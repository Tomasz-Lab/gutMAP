#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
         M E T A G E N O M I C S - P I P E L I N E
         =========================================
         """

Channel
    .fromPath(params.input)
    .splitCsv(header:true)
    .map { row -> tuple(row.sample, row.sra_id) }
    .set { ch_samples }

// ========================================================================================
//  WORKFLOW

workflow {
    // --- 1. Get Reads ---
    ENA_DOWNLOAD(ch_samples)
    ch_raw_reads = ENA_DOWNLOAD.out.reads

    // --- 2. Quality Control & Decontamination ---
    FASTP_QC(ch_raw_reads)
    REMOVE_HUMAN_READS(FASTP_QC.out.reads)
    ch_non_human_reads = REMOVE_HUMAN_READS.out.reads

    // --- 3. Metagenomic Assembly (Shared Step) ---
    MEGAHIT_ASSEMBLY(ch_non_human_reads)
    ch_assembly = MEGAHIT_ASSEMBLY.out.assembly

    // =======================================================
    // --- BRANCH 1: Taxonomic Annotation with Sylph ---
    SYLPH_TAXONOMY(ch_non_human_reads)
    // =======================================================

    // =======================================================
    // --- BRANCH 2: Gene Catalog and Abundance ---
    BAKTA_ANNOTATION(ch_assembly)

    // Create a combined channel with the required inputs for catalog creation
    ch_for_catalog_creation = ch_non_human_reads
        .join(BAKTA_ANNOTATION.out.nucleotides)
        .join(BAKTA_ANNOTATION.out.tsv)
        .map { sample_id, reads, ffn, tsv -> [ [id: sample_id], reads, ffn, tsv ] }

    CREATE_CATALOG_AND_INDEX(ch_for_catalog_creation)

    // Combine reads with the newly created index for alignment
    ch_for_alignment = ch_non_human_reads
        .combine(CREATE_CATALOG_AND_INDEX.out.index)
        .map { sample_id, reads, meta, index_files -> [ meta, reads, index_files ] }

    ALIGN_AND_QUANTIFY_READS(ch_for_alignment)

    // Combine the read counts with the original tsv for final annotation
    ch_for_annotation = ALIGN_AND_QUANTIFY_READS.out.idxstats
        .combine(BAKTA_ANNOTATION.out.tsv)
        .map { meta, idxstats, sample_id, tsv -> [ meta, idxstats, tsv ] }

    CALCULATE_TPM_AND_ANNOTATE(ch_for_annotation)
    // =======================================================

    // =======================================================
    // --- BRANCH 3: Metagenome-Assembled Genome (MAG) Creation ---
    ch_input_for_mapping = ch_non_human_reads.join(ch_assembly)
    MAP_FOR_BINNING(ch_input_for_mapping)
    METABAT2_BINNING(MAP_FOR_BINNING.out.bam)
    CHECKM_QA(METABAT2_BINNING.out.bins)
    joined_for_annotation = METABAT2_BINNING.out.bins.join(CHECKM_QA.out.summary)
    FILTER_AND_ANNOTATE(joined_for_annotation)
    // =======================================================

    // --- FINAL STEP: Collect all specified results into one directory ---
    COLLECT_RESULTS(
        FASTP_QC.out.html_report.collect(),
        FASTP_QC.out.json_report.collect(),
        SYLPH_TAXONOMY.out.collect(),
        MEGAHIT_ASSEMBLY.out.assembly.map{ it[1] }.collect(),
        BAKTA_ANNOTATION.out.annotation_dir.map{ it[1] }.collect(),
        METABAT2_BINNING.out.bins.map{ it[1] }.collect(),
        CHECKM_QA.out.checkm_dir.map{ it[1] }.collect()
    )
}

// ========================================================================================
//  PROCESSES
// ========================================================================================

process ENA_DOWNLOAD {
    tag "ENA Download: $sample_id"
    publishDir "$params.outdir/published/01_raw_reads/$sample_id", mode: 'symlink', pattern: "*.fastq.gz"
    memory '5.GB'
    cpus 2

    input:
        tuple val(sample_id), val(sra_id)

    output:
        tuple val(sample_id), path("*.fastq.gz"), emit: reads

    script:
    """
    #!/bin/bash
    set -e

    urls=\$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${sra_id}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n +2 | cut -f2 | tr ';' '\n' | sed 's|^|ftp://|')

    for url in \$urls; do
        echo "Downloading \$url ..."
        wget -c "\$url" -O "\$(basename \$url)" &
    done

    wait
    """
}

process FASTP_QC {
    tag "fastp: $sample_id"
    publishDir "$params.outdir/published/02_fastp_qc/$sample_id", mode: 'symlink'
    conda "fastp"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("*.trimmed.fastq.gz"), emit: reads
        path "*.html", emit: html_report
        path "*.json", emit: json_report

    script:
    def (r1, r2) = reads
    """
    fastp --in1 ${r1} --in2 ${r2} --out1 ${sample_id}_1.trimmed.fastq.gz --out2 ${sample_id}_2.trimmed.fastq.gz --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json --thread ${task.cpus} --detect_adapter_for_pe
    """
}

process REMOVE_HUMAN_READS {
    tag "Bowtie2 Decontam: $sample_id"
    publishDir "$params.outdir/published/03_non_human_reads/$sample_id", mode: 'symlink'
    conda "bowtie2"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.nonhuman.fastq.{1,2}.gz"), emit: reads

    script:
    def (r1, r2) = reads
    """
    bowtie2 -x ${params.bowtie2_hg38_index} -1 ${r1} -2 ${r2} --threads ${task.cpus} --very-sensitive-local --un-conc-gz ${sample_id}.nonhuman.fastq.gz -S /dev/null
    """
}

process MEGAHIT_ASSEMBLY {
    tag "MEGAHIT: $sample_id"
    publishDir "$params.outdir/published/04_shared_assembly/$sample_id", mode: 'symlink'
    conda "megahit"
    memory '32.GB'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_megahit_out"), emit: assembly

    script:
    def (r1, r2) = reads
    """
    megahit -1 ${r1} -2 ${r2} -o ${sample_id}_megahit_out -t ${task.cpus} --min-contig-len 1500
    """
}

process SYLPH_TAXONOMY {
    tag "Sylph: $sample_id"
    publishDir "$params.outdir/published/b1_taxonomy/$sample_id", mode: 'symlink'
    conda "sylph"

    input:
        tuple val(sample_id), path(reads)

    output:
        path "${sample_id}.sylph_profile.tsv"

    script:
    def (r1, r2) = reads
    """
    sylph profile ${params.sylph_db} -1 ${r1} -2 ${r2} -o ${sample_id}.sylph_profile.tsv -p ${task.cpus}
    """
}

process BAKTA_ANNOTATION {
    tag "Bakta: $sample_id"
    publishDir "$params.outdir/published/b2_annotation/$sample_id", mode: 'symlink'
    conda "bakta=1.9.0 ncbi-amrfinderplus" // bakta is updated to 1.11.x

    input:
        tuple val(sample_id), path(contigs_fa)

    output:
        tuple val(sample_id), path("${sample_id}.bakta/*.faa"), emit: proteins
        tuple val(sample_id), path("${sample_id}.bakta/*.ffn"), emit: nucleotides
        tuple val(sample_id), path("${sample_id}.bakta/*.tsv"), emit: tsv
        tuple val(sample_id), path("${sample_id}.bakta"), emit: annotation_dir

    script:
    """
    MEMFS_DB=\$MEMFS/db
    echo \$MEMFS/db
    cp ${params.bakta_db} \$MEMFS_DB -r

    bakta \
    --db \$MEMFS_DB \
    --output ${sample_id}.bakta \
    --prefix ${sample_id} \
    --threads ${task.cpus} \
    --debug \
    --verbose \
    --skip-trna \
    --skip-tmrna \
    --skip-rrna \
    --skip-ncrna \
    --skip-ncrna-region \
    --skip-crispr \
    --skip-pseudo \
    --skip-gap \
    --skip-ori \
    --skip-filter \
    --skip-plot \
    ${contigs_fa}/final.contigs.fa

    rm -rf \$MEMFS/db
    """
}

process CREATE_CATALOG_AND_INDEX {
    tag "Catalog & Index for ${meta.id}"
    conda "bioconda::cd-hit=4.8.1 bioconda::bwa-mem2=2.2.1"
    publishDir "$params.outdir/published/b2_gene_catalog/${meta.id}", mode: 'symlink'

    input:
        tuple val(meta), path(reads), path(ffn), path(tsv)

    output:
        tuple val(meta), path("${meta.id}_clustered_catalog.fna*"), emit: index

    script:
    def sampleid = meta.id
    def catalog_fna = "${sampleid}_clustered_catalog.fna"
    """
    # Step 1: Cluster at 100%
    cd-hit-est -i ${ffn} -o ${catalog_fna} -c 1.0 -T ${task.cpus} -d 0

    # Step 2: Index the gene catalogue
    bwa-mem2 index ${catalog_fna}
    """
}

process ALIGN_AND_QUANTIFY_READS {
    tag "Align & Count for ${meta.id}"
    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.19.2"
    publishDir "$params.outdir/published/b2_gene_quantification/${meta.id}", mode: 'symlink'

    input:
        tuple val(meta), path(reads), path(index_files)

    output:
        tuple val(meta), path("${meta.id}.idxstats.txt"), emit: idxstats

    script:
    def sampleid = meta.id
    def sorted_bam = "${sampleid}.bam"
    // Find the .fna file in the index file list to use as the base for alignment
    def index_base = index_files.find { it.name.endsWith('.fna') }
    """
    # Step 3: Align reads to the catalogue and create a sorted BAM file
    bwa-mem2 mem -t ${task.cpus} ${index_base} ${reads[0]} ${reads[1]} | \\
        samtools view -bS -@ ${task.cpus} - | \\
        samtools sort -@ ${task.cpus} -o ${sorted_bam}

    # Step 4: Index the sorted BAM file and get read counts
    samtools index ${sorted_bam}
    samtools idxstats ${sorted_bam} > ${sampleid}.idxstats.txt
    """
}

process CALCULATE_TPM_AND_ANNOTATE {
    tag "TPM & Annotate for ${meta.id}"

    input:
        tuple val(meta), path(idxstats), path(tsv_file)

    output:
        tuple val(meta), path("${meta.id}_gene_quantification.txt"), emit: quantification

    script:
    def sampleid = meta.id
    def abundances_file = "${sampleid}_gene_abundances.txt"
    """
    python3 - <<'EOF'
    import sys

    gene_lengths = {}
    read_counts = {}
    rpk_sum = 0.0

    with open("${idxstats}", 'r') as f:
        for line in f:
            parts = line.strip().split('\\t')
            gene_id, length, mapped_reads, _ = parts
            if gene_id == '*' or int(length) == 0:
                continue

            length = int(length)
            mapped_reads = int(mapped_reads)
            gene_lengths[gene_id] = length
            read_counts[gene_id] = mapped_reads
            rpk = mapped_reads / (length / 1000.0)
            rpk_sum += rpk

    scaling_factor = rpk_sum / 1_000_000.0 if rpk_sum > 0 else 0

    with open("${abundances_file}", 'w') as out_file:
        out_file.write("gene_id\\tlength\\tread_count\\ttpm\\n")
        for gene_id in gene_lengths:
            length = gene_lengths[gene_id]
            reads = read_counts[gene_id]
            rpk = reads / (length / 1000.0)
            tpm = rpk / scaling_factor if scaling_factor > 0 else 0
            out_file.write(f"{gene_id}\\t{length}\\t{reads}\\t{tpm}\\n")
    EOF

    # Join TPM values with gene descriptions
    awk 'BEGIN{FS=OFS="\\t"} FNR==NR{if(FNR>1)desc[\$6]=\$8; next} {if(FNR==1)print "gene_id","description","tpm"; else print \$1,desc[\$1],\$4}' ${tsv_file} ${abundances_file} > ${sampleid}_gene_quantification.txt
    """
}

// =======================================================
// --- BRANCH 3 PROCESSES ---

process MAP_FOR_BINNING {
    tag "Map for Binning: $sample_id"
    publishDir "$params.outdir/published/b3_mags/01_mapping/$sample_id", mode: 'symlink'
    conda "bowtie2 samtools"

    input:
        tuple val(sample_id), path(reads), path(megahit_dir)

    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam"), path(megahit_dir), emit: bam


    script:
    def (r1, r2) = reads
    def contigs = "${megahit_dir}/final.contigs.fa"
    """
    bowtie2-build ${contigs} ${sample_id}.assembly.idx
    bowtie2 \\
        -x ${sample_id}.assembly.idx \\
        -1 ${r1} \\
        -2 ${r2} \\
        --threads ${task.cpus} \\
        -S temp.sam

    samtools view -bS -F 4 temp.sam | samtools sort - > ${sample_id}.sorted.bam
    rm temp.sam
    """
}

process METABAT2_BINNING {
    tag "MetaBAT2 Binning: $sample_id"
    publishDir "$params.outdir/published/b3_mags/02_binning/$sample_id", mode: 'symlink'
    conda "metabat2"

    input:
        tuple val(sample_id), path(bam), path(megahit_dir)

    output:
        tuple val(sample_id), path("${sample_id}_bins"), emit: bins

    script:
    def contigs = "${megahit_dir}/final.contigs.fa"
    """
    jgi_summarize_bam_contig_depths --outputDepth ${sample_id}.depth.txt ${bam}
    mkdir -p ${sample_id}_bins
    metabat2 -i ${contigs} -a ${sample_id}.depth.txt -o ${sample_id}_bins/bin -t ${task.cpus} -m 1500 -v
    """
}

process CHECKM_QA {
    tag "CheckM: $sample_id"
    publishDir "$params.outdir/published/b3_mags/03_checkm_qa/$sample_id", mode: 'symlink'
    conda "checkm2"

    input:
        tuple val(sample_id), path(bins_dir)

    output:
        tuple val(sample_id), path("${sample_id}_checkm2_out"), emit: checkm_dir
        tuple val(sample_id), path("${sample_id}_checkm2_out/diamond_output/DIAMOND_RESULTS.tsv"), emit: summary


    script:
    """
    checkm2 predict --input ${bins_dir} --output-directory ${sample_id}_checkm2_out --threads ${task.cpus} -x fa
    """
}

process FILTER_HQ_MAGS {
    tag "Filter HQ MAGs"
    publishDir "$params.outdir/published/b3_mags/04_high_quality_mags", mode: 'symlink'

    input:
        path checkm_summaries

    output:
        path "high_quality_mags.txt"

    script:
    """
    echo -e "Bin Id\\tMarker lineage\\t# genomes\\t# markers\\t# marker sets\\t0\\t1\\t2\\t3\\t4\\t5+\\tCompleteness\\tContamination\\tStrain heterogeneity" > combined_summary.tsv
    for f in ${checkm_summaries.join(' ')}; do
        tail -n +2 "\$f" >> combined_summary.tsv
    done
    awk -F'\\t' 'BEGIN{OFS=FS} { if (\$12 > 90 && \$13 < 5) { print \$1 } }' combined_summary.tsv > high_quality_mags.txt
    """
}

process FILTER_AND_ANNOTATE {
    tag "Filter & Annotate GTDB-Tk: $sample_id"
    conda "/net/afscra/people/plgpkica/metagenome_proj/conda/gtdbtk_and_pandas"

    input:
        tuple val(sample_id), path(bins_dir), path(checkm_summary)

    output:
        tuple val(sample_id), path("${sample_id}_gtdbtk_out"), emit: gtdbtk_results

    script:
    def min_completeness = 90
    def max_contamination = 5
    """
    mkdir hq_bins -p

    # Use a here-document to create a python filtering script on-the-fly
    cat << EOF > filter_bins.py
    #!/usr/bin/env python
    import pandas as pd
    import os
    import shutil
    import sys
    import ast

    checkm_file = sys.argv[1]
    bins_dir = sys.argv[2]
    out_dir = sys.argv[3]
    completeness = float(sys.argv[4])
    contamination = float(sys.argv[5])

    # This list will hold the parsed data from each line
    parsed_data = []

    # Open and read the file line by line
    with open(checkm_file, 'r') as f:
        for line in f:
            # Skip any empty lines
            if not line.strip():
                continue

            # Split the line into two parts at the first whitespace
            # Part 1: bin_name (e.g., 'bin.15')
            # Part 2: dict_string (e.g., "{'marker lineage': ...}")
            try:
                bin_name, dict_string = line.strip().split(None, 1)

                # Safely evaluate the string to convert it into a Python dictionary
                data_dict = ast.literal_eval(dict_string)

                # Add the bin name to the dictionary
                data_dict['Bin Id'] = bin_name

                # Add the complete dictionary to our list
                parsed_data.append(data_dict)

            except (ValueError, SyntaxError) as e:
                print(f"Skipping malformed line: {line.strip()} - Error: {e}")

    # Create the DataFrame from the list of dictionaries
    df = pd.DataFrame(parsed_data)

    # Filter the DataFrame
    hq_df = df[(df['Completeness'] >= completeness) & (df['Contamination'] <= contamination)]

    hq_bin_ids = hq_df['Bin Id'].tolist()

    # Copy the corresponding FASTA files
    for bin_id in hq_bin_ids:
            source_path = os.path.join(bins_dir, f"{bin_id}.fa")
            dest_path = os.path.join(out_dir, f"{bin_id}.fa")
            if os.path.exists(source_path):
                shutil.copy(source_path, dest_path)
    EOF

    # Step 1: Execute the python script to populate the hq_bins directory
    python filter_bins.py \\
        "${checkm_summary}" \\
        "${bins_dir}" \\
        "hq_bins" \\
        "${min_completeness}" \\
        "${max_contamination}"
    # Step 2: Run GTDB-Tk, but only if high-quality bins were actually found
    if [ -n "\$(find hq_bins -name '*.fa')" ]; then
        gtdbtk classify \\
        --genome_dir hq_bins \\
        --align_dir hq_align \\
        --out_dir ${sample_id}_gtdbtk_out \\
        --extension fa \\
        --cpus ${task.cpus}
    else
        mkdir ${sample_id}_gtdbtk_out
        touch ${sample_id}_gtdbtk_out/NO_HQ_BINS_FOUND.txt
    fi
    """
}

process COLLECT_RESULTS {
    tag "Collecting all results"
    publishDir "${params.outdir}/collected", mode: 'move'

    input:
        path fastp_html_reports
        path fastp_json_reports
        path sylph_profiles
        path megahit_assemblies
        path bakta_annotation_dirs
        path metabat_bin_dirs
        path checkm_output_dirs

    output:
        path "qc_reports"
        path "taxonomy_profiles"
        path "assemblies"
        path "annotation"
        path "mags_and_quality"

    script:
    """
    mkdir -p qc_reports
    mkdir -p taxonomy_profiles
    mkdir -p assemblies
    mkdir -p annotation
    mkdir -p mags_and_quality/raw_bins
    mkdir -p mags_and_quality/checkm_output

    echo "Collecting FASTP QC reports..."
    cp ${fastp_html_reports.join(' ')} qc_reports/
    cp ${fastp_json_reports.join(' ')} qc_reports/

    echo "Collecting SYLPH taxonomy profiles..."
    cp ${sylph_profiles.join(' ')} taxonomy_profiles/

    echo "Collecting MEGAHIT assemblies..."
    cp -r ${megahit_assemblies.join(' ')} assemblies/

    echo "Collecting Bakta annotations..."
    cp -r ${bakta_annotation_dirs.join(' ')} annotation/

    echo "Collecting MetaBAT2 raw bins..."
    cp -r ${metabat_bin_dirs.join(' ')} mags_and_quality/raw_bins/

    echo "Collecting CheckM output directories..."
    cp -r ${checkm_output_dirs.join(' ')} mags_and_quality/checkm_output/

    echo "All results collected successfully!"
    """
}