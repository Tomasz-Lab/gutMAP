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
    // --- BRANCH 1: Taxonomic Annotation with Sylph, GTDB version 226 ---
    SYLPH_TAXONOMY(ch_non_human_reads)
    // =======================================================

    // =======================================================
    // --- BRANCH 2: Gene Catalog and Gene Abundance ---
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
    // --- BRANCH 3: Metagenome-Assembled Genomes (MAG) Creation ---
    ch_input_for_mapping = ch_non_human_reads.join(ch_assembly)
    MAP_FOR_BINNING(ch_input_for_mapping)
    METABAT2_BINNING(MAP_FOR_BINNING.out.bam)
    CHECKM_QA(METABAT2_BINNING.out.bins)
    joined_for_annotation = METABAT2_BINNING.out.bins.join(CHECKM_QA.out.checkm_summary)
    FILTER_AND_ANNOTATE(joined_for_annotation)
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

// =======================================================
// --- BRANCH 1 PROCESSES ---

process SYLPH_TAXONOMY {
    tag "Sylph: $sample_id"
    publishDir "$params.outdir/published/branch_1/01_taxonomy/$sample_id", mode: 'symlink'
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

// =======================================================
// --- BRANCH 2 PROCESSES ---

process BAKTA_ANNOTATION {
    tag "Bakta: $sample_id"
    publishDir "$params.outdir/published/branch_2/01_annotation/$sample_id", mode: 'symlink'
    conda "bakta ncbi-amrfinderplus"

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
    publishDir "$params.outdir/published/branch_2/02_gene_catalog/${meta.id}", mode: 'symlink'
    conda "bioconda::cd-hit=4.8.1 bioconda::bwa-mem2=2.2.1"

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
    publishDir "$params.outdir/published/branch_2/03_gene_quantification/${meta.id}", mode: 'symlink'
    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.19.2"

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
    publishDir "$params.outdir/published/branch_2/04_tpm_and_annotate/${meta.id}", mode: 'symlink'

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
    publishDir "$params.outdir/published/branch_3/01_mapping/$sample_id", mode: 'symlink'
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
    publishDir "$params.outdir/published/branch_3/02_binning/$sample_id", mode: 'symlink'
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
    publishDir "$params.outdir/published/branch_3/03_checkm_qa/$sample_id", mode: 'symlink'
    conda "checkm2"

    input:
        tuple val(sample_id), path(bins_dir)

    output:
        tuple val(sample_id), path("${sample_id}_checkm2_out"), emit: checkm_dir
        tuple val(sample_id), path("${sample_id}_checkm2_out/quality_report.tsv"), emit: checkm_summary


    script:
    """
    checkm2 predict --input ${bins_dir} --output-directory ${sample_id}_checkm2_out --threads ${task.cpus} -x fa
    """
}

process FILTER_AND_ANNOTATE {
    tag "Filter & Annotate GTDB-Tk: $sample_id"
    publishDir "$params.outdir/published/branch_3/04_filter_and_annotate/$sample_id", mode: 'symlink'
    conda "python=3.12 gtdbtk pandas"

    input:
        tuple val(sample_id), path(bins_dir), path(checkm_summary)

    output:
        tuple val(sample_id), path("${sample_id}_gtdbtk_out"), emit: gtdbtk_results

    script:
    def min_completeness = 90
    def max_contamination = 5
    """
    # Use a here-document to create a python filtering script on-the-fly
    cat << EOF > filter_bins.py
    #!/usr/bin/env python
    import pandas as pd
    import os
    import sys
    import shutil

    checkm_file = sys.argv[1]
    bins_dir = sys.argv[2]
    out_dir = sys.argv[3]
    completeness = float(sys.argv[4])
    contamination = float(sys.argv[5])

    df = pd.read_csv(checkm_file, sep='\t', header=0, index_col=None, on_bad_lines='skip').dropna(subset=['Completeness', 'Contamination'])

    hq_df = df[(df['Completeness'] >= completeness) & (df['Contamination'] <= contamination)]
    hq_bin_ids = hq_df['Name'].tolist()

    os.makedirs(out_dir, exist_ok=True)
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
        gtdbtk classify_wf \\
        --genome_dir hq_bins \\
        --out_dir ${sample_id}_gtdbtk_out \\
        --cpus ${task.cpus} \\
        --extension fa
    else
        mkdir ${sample_id}_gtdbtk_out
        touch ${sample_id}_gtdbtk_out/NO_HQ_BINS_FOUND.txt
    fi
    """
}
