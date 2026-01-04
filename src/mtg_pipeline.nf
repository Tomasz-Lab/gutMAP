#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
         M E T A G E N O M I C S - P I P E L I N E
         =========================================
         """

// ========================================================================================
//  WORKFLOW

workflow {
    // --- 1. Get Reads ---
    ch_inputs = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> row.sra_id }
        .map { sra_id ->
            def sra_file = file("${params.indir}/${sra_id}/${sra_id}.sra")
            def fq1      = file("${params.indir}/${sra_id}/${sra_id}_1.fastq.gz")
            def fq2      = file("${params.indir}/${sra_id}/${sra_id}_2.fastq.gz")
            def fq       = file("${params.indir}/${sra_id}/${sra_id}.fastq.gz")

            if (sra_file.exists()) {
                return tuple(sra_id, sra_file, 'sra')
            } else if (fq1.exists()) {
                return tuple(sra_id, [fq1, fq2], 'fastq')
            } else if (fq.exists()) {
                log.warn "Single fastq input file ${sra_id} - not implemented"
            } else {
                log.warn "Missing input files for ${sra_id}"
            }
        }
        .branch {
            to_convert: it[2] == 'sra'
            fastq_gz:  it[2] == 'fastq'
        }


    ch_sra_input = ch_inputs.to_convert.map { id, file, type -> tuple(id, file) }

    CHECK_AND_CONVERT_READS(ch_sra_input)
    ch_fastq_gz = ch_inputs.fastq_gz.map { id, files, type -> tuple(id, files) }

    ch_ready_for_qc = ch_fastq_gz.mix(CHECK_AND_CONVERT_READS.out.fastq)

    // --- 2. Quality Control & Decontamination ---
    FASTP_QC(ch_ready_for_qc)
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

    CREATE_CATALOG_AND_INDEX(ch_for_catalog_creation)

    // Combine reads with the newly created index for alignment
    ch_for_alignment = ch_non_human_reads
        .join(CREATE_CATALOG_AND_INDEX.out.index)

    ALIGN_AND_QUANTIFY_READS(ch_for_alignment)

    // Combine the read counts with the original tsv for final annotation
    ch_for_annotation = ALIGN_AND_QUANTIFY_READS.out.idxstats
        .join(BAKTA_ANNOTATION.out.tsv)

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
    tag "ENA Download: $sra_id"
    conda "kingfisher"

    input:
        val(sra_id)

    output:
        tuple val(sra_id), path("*.fastq.gz"), emit: reads

    script:
    """
    kingfisher get -r $sra_id -t ${task.cpus} -f fastq.gz -m ena-ftp prefetch aws-http
    """
}

process CHECK_AND_CONVERT_READS {
    conda "kingfisher"  // kingfisher includes fasterq-dump

    input:
        tuple val(sra_id), path(sra_file)

    output:
        tuple val(sra_id), path("*.fastq"), emit: fastq

    script:
    """
    fasterq-dump "${sra_file}" --threads ${task.cpus}
    """
}

process FASTP_QC {
    tag "fastp: $sra_id"
    publishDir "$params.outdir/published/02_fastp_qc/$sra_id", mode: 'copy', pattern: "*.json"
    conda "fastp"

    input:
        tuple val(sra_id), path(reads)

    output:
        tuple val(sra_id), path("*.trimmed.fastq"), emit: reads
        tuple val(sra_id), path("*.html"), emit: html_report
        tuple val(sra_id), path("*.json"), emit: json_report

    script:
    def (r1, r2) = reads
    """
    fastp --in1 ${r1} --in2 ${r2} --out1 ${sra_id}_1.trimmed.fastq --out2 ${sra_id}_2.trimmed.fastq \\
        --html ${sra_id}.fastp.html --json ${sra_id}.fastp.json --thread ${task.cpus} --detect_adapter_for_pe
    """
}

process REMOVE_HUMAN_READS {
    tag "Bowtie2 Decontam: $sra_id"
    conda "bowtie2"

    input:
        tuple val(sra_id), path(reads)

    output:
        tuple val(sra_id), path("${sra_id}.nonhuman.{1,2}.fastq"), emit: reads

    script:
    def (r1, r2) = reads
    """
    bowtie2 -x ${params.bowtie2_hg38_index} -1 ${r1} -2 ${r2} --threads ${task.cpus} --very-sensitive-local \\
        --un-conc ${sra_id}.nonhuman.fastq -S /dev/null
    """
}

process MEGAHIT_ASSEMBLY {
    tag "MEGAHIT: $sra_id"
    publishDir "$params.outdir/published/04_shared_assembly/$sra_id", mode: 'copy', pattern: "final.contigs.fa"
    conda "megahit"

    input:
        tuple val(sra_id), path(reads)

    output:
        tuple val(sra_id), path("${sra_id}_megahit_out"), emit: assembly
        tuple val(sra_id), path("final.contigs.fa")

    script:
    def (r1, r2) = reads
    """
    megahit -1 ${r1} -2 ${r2} -o ${sra_id}_megahit_out -t ${task.cpus} --min-contig-len 1500
    cp ${sra_id}_megahit_out/final.contigs.fa final.contigs.fa
    """
}

// =======================================================
// --- BRANCH 1 PROCESSES ---

process SYLPH_TAXONOMY {
    tag "Sylph: $sra_id"
    publishDir "$params.outdir/published/branch_1/01_taxonomy/$sra_id", mode: 'copy', pattern: "${sra_id}.sylph_profile.tsv"
    conda "sylph"

    input:
        tuple val(sra_id), path(reads)

    output:
        tuple val(sra_id), path("${sra_id}.sylph_profile.tsv")

    script:
    def (r1, r2) = reads
    """
    sylph profile ${params.sylph_db} -1 ${r1} -2 ${r2} -o ${sra_id}.sylph_profile.tsv -p ${task.cpus}
    """
}

// =======================================================
// --- BRANCH 2 PROCESSES ---

process BAKTA_ANNOTATION {
    tag "Bakta: $sra_id"
    publishDir "$params.outdir/published/branch_2/01_annotation/$sra_id", mode: 'copy', pattern: "${sra_id}.{faa,ffn,tsv,gff3,txt}"
    conda "bakta ncbi-amrfinderplus"

    input:
        tuple val(sra_id), path(contigs_fa)

    output:
        tuple val(sra_id), path("${sra_id}.bakta/*.faa"), emit: proteins
        tuple val(sra_id), path("${sra_id}.bakta/*.ffn"), emit: nucleotides
        tuple val(sra_id), path("${sra_id}.bakta/${sra_id}.tsv"), emit: tsv
        tuple val(sra_id), path("${sra_id}.bakta"), emit: annotation_dir
        tuple val(sra_id), path("*.{faa,ffn,tsv,gff3,txt}")

    script:
    """
    MEMFS_DB=\$MEMFS/db
    echo \$MEMFS/db
    cp ${params.bakta_db} \$MEMFS_DB -r

    bakta \
    --db \$MEMFS_DB \
    --output ${sra_id}.bakta \
    --prefix ${sra_id} \
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
    cp ${sra_id}.bakta/* .
    """
}

process CREATE_CATALOG_AND_INDEX {
    tag "Catalog & Index for ${sra_id}"
    publishDir "$params.outdir/published/branch_2/02_gene_catalog/${sra_id}", mode: 'copy', pattern: "${sra_id}_clustered_catalog.fna"
    conda "bioconda::cd-hit=4.8.1 bioconda::bwa-mem2=2.2.1"

    input:
        tuple val(sra_id), path(reads), path(ffn), path(tsv)

    output:
        tuple val(sra_id), path("${sra_id}_clustered_catalog.fna*"), emit: index

    script:
    def catalog_fna = "${sra_id}_clustered_catalog.fna"
    """
    # Step 1: Cluster at 100%
    cd-hit-est -i ${ffn} -o ${catalog_fna} -c 1.0 -T ${task.cpus} -d 0 -M 0

    # Step 2: Index the gene catalogue
    bwa-mem2 index ${catalog_fna}
    """
}

process ALIGN_AND_QUANTIFY_READS {
    tag "Align & Count for ${sra_id}"
    publishDir "$params.outdir/published/branch_2/03_gene_quantification/${sra_id}", mode: 'copy', pattern: "${sra_id}.idxstats.txt"
    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.19.2"

    input:
        tuple val(sra_id), path(reads), path(index_files)

    output:
        tuple val(sra_id), path("${sra_id}.idxstats.txt"), emit: idxstats

    script:
    def sorted_bam = "${sra_id}.bam"
    // Find the .fna file in the index file list to use as the base for alignment
    def index_base = index_files.find { it.name.endsWith('.fna') }
    """
    # Step 3: Align reads to the catalogue and create a sorted BAM file
    bwa-mem2 mem -t ${task.cpus} ${index_base} ${reads[0]} ${reads[1]} | \\
        samtools view -b - | samtools sort --threads ${task.cpus} -o ${sorted_bam}

    # Step 4: Index the sorted BAM file and get read counts
    samtools index --threads ${task.cpus} ${sorted_bam}
    samtools idxstats --threads ${task.cpus} ${sorted_bam} > ${sra_id}.idxstats.txt
    """
}

process CALCULATE_TPM_AND_ANNOTATE {
    tag "TPM & Annotate for ${sra_id}"
    publishDir "$params.outdir/published/branch_2/04_tpm_and_annotate/${sra_id}", mode: 'copy', pattern: "${sra_id}_gene_quantification.tsv"
    conda "pandas"

    input:
        tuple val(sra_id), path(idxstats), path(tsv_file)

    output:
        tuple val(sra_id), path("${sra_id}_gene_quantification.tsv"), emit: quantification

    script:
        def merged_output = "${sra_id}_gene_quantification.tsv"

        """
        python3 - <<'EOF'
        import pandas as pd
        from io import StringIO

        # -----------------------------
        # Read idxstats file and compute TPM
        # -----------------------------
        gene_lengths = {}
        read_counts = {}
        rpk_sum = 0.0

        with open("${idxstats}", "r") as f:
            for line in f:
                parts = line.strip().split("\\t")
                gene_id, length, mapped_reads, _ = parts

                if gene_id == "*" or int(length) == 0:
                    continue

                length = int(length)
                mapped_reads = int(mapped_reads)

                gene_lengths[gene_id] = length
                read_counts[gene_id] = mapped_reads

                rpk = mapped_reads / (length / 1000.0)
                rpk_sum += rpk

        scaling_factor = rpk_sum / 1_000_000.0 if rpk_sum > 0 else 0

        abundances_data = []
        for gene_id in gene_lengths:
            length = gene_lengths[gene_id]
            reads = read_counts[gene_id]
            rpk = reads / (length / 1000.0)
            tpm = rpk / scaling_factor if scaling_factor > 0 else 0

            abundances_data.append({
                "gene_id": gene_id,
                "length": length,
                "read_count": reads,
                "tpm": tpm
            })

        abundances_df = pd.DataFrame(abundances_data)

        # -----------------------------
        # Read annotation TSV
        # (keep only last '#' line as header)
        # -----------------------------
        with open("${tsv_file}", "r") as f:
            lines = f.readlines()

        header_idx = max(i for i, line in enumerate(lines) if line.startswith("#"))
        header_line = lines[header_idx].lstrip("#").strip()
        data_lines = lines[header_idx + 1:]

        tsv_str = header_line + "\\n" + "".join(data_lines)
        tsv_df = pd.read_csv(StringIO(tsv_str), sep="\\t", skipinitialspace=True)

        # -----------------------------
        # Merge and write output
        # -----------------------------
        merged_df = tsv_df.merge(
            abundances_df,
            left_on="Locus Tag",
            right_on="gene_id",
            how="inner"
        )

        merged_df.to_csv("${merged_output}", sep="\\t", index=False)
        EOF
        """
}

// =======================================================
// --- BRANCH 3 PROCESSES ---

process MAP_FOR_BINNING {
    tag "Map for Binning: $sra_id"
    conda "bowtie2 samtools"

    input:
        tuple val(sra_id), path(reads), path(megahit_dir)

    output:
        tuple val(sra_id), path("${sra_id}.sorted.bam"), path(megahit_dir), emit: bam


    script:
    def (r1, r2) = reads
    def contigs = "${megahit_dir}/final.contigs.fa"
    """
    bowtie2-build --threads ${task.cpus} -f ${contigs} ${sra_id}.assembly.idx
    bowtie2 -x ${sra_id}.assembly.idx -1 ${r1} -2 ${r2} --threads ${task.cpus} -S temp.sam

    samtools view -u -b -F 4 temp.sam | samtools sort --threads ${task.cpus} -m 4G -o ${sra_id}.sorted.bam
    rm temp.sam
    """
}

process METABAT2_BINNING {
    tag "MetaBAT2 Binning: $sra_id"
    conda "metabat2"

    input:
        tuple val(sra_id), path(bam), path(megahit_dir)

    output:
        tuple val(sra_id), path("${sra_id}_bins"), emit: bins

    script:
    def contigs = "${megahit_dir}/final.contigs.fa"
    """
    jgi_summarize_bam_contig_depths --outputDepth ${sra_id}.depth.txt ${bam}
    mkdir -p ${sra_id}_bins
    metabat2 -i ${contigs} -a ${sra_id}.depth.txt -o ${sra_id}_bins/bin -t ${task.cpus} -m 1500 -v
    """
}

process CHECKM_QA {
    tag "CheckM: $sra_id"
    publishDir "$params.outdir/published/branch_3/03_checkm_qa/$sra_id", mode: 'copy', pattern: "{quality_report.tsv,DIAMOND_RESULTS.tsv}"
    conda "checkm2"

    input:
        tuple val(sra_id), path(bins_dir)

    output:
        tuple val(sra_id), path("quality_report.tsv"), emit: checkm_summary
        tuple val(sra_id), path("DIAMOND_RESULTS.tsv")


    // cp is to make pattern in publishDir work correctly, does not work when publishing from any subdirectory
    script:
    """
    checkm2 predict --input ${bins_dir} --output-directory ${sra_id}_checkm2_out --threads ${task.cpus} -x fa
    cp ${sra_id}_checkm2_out/quality_report.tsv quality_report.tsv
    cp ${sra_id}_checkm2_out/diamond_output/DIAMOND_RESULTS.tsv DIAMOND_RESULTS.tsv
    """
}

process FILTER_AND_ANNOTATE {
    tag "Filter & Annotate GTDB-Tk: $sra_id"
    publishDir "$params.outdir/published/branch_3/04_filter_and_annotate/$sra_id", mode: 'copy', pattern: '{gtdbtk.bac120.summary.tsv,hq_bins}'
    conda "python=3.12 gtdbtk pandas"

    input:
        tuple val(sra_id), path(bins_dir), path(checkm_summary)

    output:
        tuple val(sra_id), path("${sra_id}_gtdbtk_out"), emit: gtdbtk_results
        tuple val(sra_id), path("hq_bins"), emit: hq_bins
        tuple val(sra_id), path("gtdbtk.bac120.summary.tsv"), emit: gtdbtk_summary

    script:
    def min_completeness = 90
    def max_contamination = 5
    """
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
        --out_dir ${sra_id}_gtdbtk_out \\
        --cpus ${task.cpus} \\
        --extension fa
        cp ${sra_id}_gtdbtk_out/classify/gtdbtk.bac120.summary.tsv gtdbtk.bac120.summary.tsv
    else
        touch NO_HQ_BINS_FOUND.tsv
    fi
    """
}
