#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def withLayout(meta, reads) {
    def names = [reads].flatten()*.name
    def single_end = !names.any { it ==~ /.*_2\.fastq(\.gz)?/ }
    if (!single_end && names.size() > 2) {
        log.warn "${meta.id}: unpaired reads are not used: ${names.join(' ')}"
    }
    tuple(meta + [single_end: single_end], reads)
}

// ========================================================================================
//  WORKFLOW

workflow {
    log.info """
         M E T A G E N O M I C S - P I P E L I N E
         =========================================
         """

    ch_inputs = Channel
        .fromPath(params.input)
        .splitCsv(header: true, limit: params.limit)
        .map { row ->
            def sra_id = row.sra_id
            def meta = [id: sra_id]
            def sra_file = file("${params.indir}/${sra_id}/${sra_id}.sra")
            def fq1 = file("${params.indir}/${sra_id}/${sra_id}_1.fastq.gz")
            def fq2 = file("${params.indir}/${sra_id}/${sra_id}_2.fastq.gz")
            def fq  = file("${params.indir}/${sra_id}/${sra_id}.fastq.gz")

            if (sra_file.exists()) {
                return tuple(meta, sra_file, 'sra')
            }
            else if (fq1.exists() && fq2.exists()) {
                return tuple(meta, [fq1, fq2], 'fastq')
            }
            else if (fq.exists()) {
                tuple(meta, [fq], 'fastq')
            }
            else {
                log.warn "Bad/missing input files for ${sra_id}"
            }
        }
        .branch {
            to_convert: it[2] == 'sra'
            fastq:  it[2] == 'fastq'
        }


    ch_sra_input = ch_inputs.to_convert.map { meta, file, type -> tuple(meta, file) }

    CHECK_AND_CONVERT_READS(ch_sra_input)
    ch_fastq = ch_inputs.fastq.map { meta, files, type -> tuple(meta, files) }

    ch_ready_for_qc = ch_fastq.mix(CHECK_AND_CONVERT_READS.out.fastq)
        .map { meta, reads -> withLayout(meta, reads) }

    // --- 2. Quality Control & Decontamination ---
    FASTP_QC(ch_ready_for_qc)
    REMOVE_HUMAN_READS(FASTP_QC.out.reads)
    ch_non_human_reads = REMOVE_HUMAN_READS.out.reads.filter { _meta, reads -> [reads].flatten().every { f -> f.size() > 0 } }

    // --- 3. Metagenomic Assembly (Shared Step) ---
    MEGAHIT_ASSEMBLY(ch_non_human_reads)
    ch_assembly = MEGAHIT_ASSEMBLY.out.contigs.filter { _meta, contigs -> contigs.size() > 0 }

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

process CHECK_AND_CONVERT_READS {
    tag "${meta.id}"
    afterScript "D=${params.outdir}/published/01_check_and_convert/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} \$D/ 2>/dev/null || true"
    conda "kingfisher"  // kingfisher includes fasterq-dump
    cpus 8
    memory "20 GB"

    input:
        tuple val(meta), path(sra_file)

    output:
        tuple val(meta), path("*.fastq"), emit: fastq

    script:
    """
    fasterq-dump "${sra_file}" --threads ${task.cpus}
    """
}

process FASTP_QC {
    tag "${meta.id}"
    publishDir "$params.outdir/published/02_fastp_qc/${meta.id}", mode: 'copy', pattern: "{*.json,*.html}"
    afterScript "D=${params.outdir}/published/02_fastp_qc/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} \$D/ 2>/dev/null || true"
    conda "fastp"

    cpus 6
    memory "30 GB"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*.trimmed.fastq"), emit: reads
        tuple val(meta), path("*.html"), emit: html_report
        tuple val(meta), path("*.json"), emit: json_report

    script:
    def in_args  = meta.single_end ? "-i ${reads}" : "--in1 *_1.fastq* --in2 *_2.fastq* --detect_adapter_for_pe"
    def out_args = meta.single_end ? "-o ${meta.id}.trimmed.fastq" : "--out1 ${meta.id}_1.trimmed.fastq --out2 ${meta.id}_2.trimmed.fastq"
    def converted = [reads].flatten().findAll { f -> f.name.endsWith('.fastq') }.join(' ')
    """
    fastp ${in_args} ${out_args} \
          --html ${meta.id}.fastp.html \
          --json ${meta.id}.fastp.json \
          --thread ${task.cpus}

    # The FASTQ files of CHECK_AND_CONVERT_READS have no other consumer. Their removal frees disk early and breaks -resume.
    # Downloaded inputs end in .fastq.gz and stay.
    for f in ${converted}; do rm -f "\$(readlink -f "\$f")"; done
    """
}

process REMOVE_HUMAN_READS {
    tag "${meta.id}"
    afterScript "D=${params.outdir}/published/03_non_human_reads/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} \$D/ 2>/dev/null || true"
    conda "bowtie2"

    cpus 16
    memory "32 GB"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.id}.nonhuman*.fastq"), emit: reads

    script:
    def reads_args = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def unaligned  = meta.single_end ? "--un" : "--un-conc"
    """
    bowtie2 -x ${params.bowtie2_hg38_index} \
            ${reads_args} \
            --threads ${task.cpus} \
            --very-sensitive-local \
            ${unaligned} ${meta.id}.nonhuman.fastq \
            -S /dev/null

    # The trimmed reads of FASTP_QC have no other consumer. Their removal frees disk early and breaks -resume.
    rm -f \$(readlink -f ${reads})
    """
}

process MEGAHIT_ASSEMBLY {
    tag "${meta.id}"
    publishDir "$params.outdir/published/04_shared_assembly/${meta.id}", mode: 'copy', pattern: "final.contigs.fa"
    afterScript "D=${params.outdir}/published/04_shared_assembly/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} log options.json \$D/ 2>/dev/null || true"
    conda "megahit"

    cpus 20
    memory "40 GB"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("final.contigs.fa"), emit: contigs

    script:
    def reads_args = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    megahit ${reads_args} -o ${meta.id}_megahit_out --tmp-dir \${MEMFS:-\$PWD} -t ${task.cpus} --min-contig-len 1500
    mv ${meta.id}_megahit_out/{final.contigs.fa,log,options.json} .
    """
}

// =======================================================
// --- BRANCH 1 PROCESSES ---

process SYLPH_TAXONOMY {
    tag "${meta.id}"
    publishDir "$params.outdir/published/branch_1/01_taxonomy/${meta.id}", mode: 'copy', pattern: "${meta.id}.sylph_profile.tsv"
    afterScript "D=${params.outdir}/published/branch_1/01_taxonomy/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} \$D/ 2>/dev/null || true"
    conda "sylph"

    cpus 2
    memory "40 GB"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.id}.sylph_profile.tsv"), emit: sylph_profile

    script:
    def reads_args = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    sylph profile ${params.sylph_db} ${reads_args} -o ${meta.id}.sylph_profile.tsv -p ${task.cpus}
    """
}

// =======================================================
// --- BRANCH 2 PROCESSES ---

process BAKTA_ANNOTATION {
    tag "${meta.id}"
    publishDir "$params.outdir/published/branch_2/01_annotation/${meta.id}", mode: 'copy', pattern: "${meta.id}.{faa,ffn,tsv,gff3,txt}"
    afterScript "D=${params.outdir}/published/branch_2/01_annotation/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} ${meta.id}.log \$D/ 2>/dev/null || true"
    conda "bioconda::bakta=1.12.2 bioconda::ncbi-amrfinderplus=4.2.7"
    cpus 5
    memory "20 GB"
    time 10.h

    input:
        tuple val(meta), path(contigs)

    output:
        tuple val(meta), path("${meta.id}.ffn"), emit: nucleotides
        tuple val(meta), path("${meta.id}.tsv"), emit: tsv
        tuple val(meta), path("${meta.id}.{faa,gff3,txt}")

    script:
    """
    bakta \\
        --db ${params.bakta_db} \\
        --prefix ${meta.id} \\
        --threads ${task.cpus} \\
        --skip-trna \\
        --skip-tmrna \\
        --skip-rrna \\
        --skip-ncrna \\
        --skip-ncrna-region \\
        --skip-crispr \\
        --skip-pseudo \\
        --skip-gap \\
        --skip-ori \\
        --skip-filter \\
        --skip-plot \\
        --meta \\
        --tmp-dir \$PWD \\
        ${contigs}
    """
}

process CREATE_CATALOG_AND_INDEX {
    tag "${meta.id}"
    publishDir "$params.outdir/published/branch_2/02_gene_catalog/${meta.id}", mode: 'copy', pattern: "${meta.id}_clustered_catalog.fna"
    afterScript "D=${params.outdir}/published/branch_2/02_gene_catalog/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} \$D/ 2>/dev/null || true"
    conda "bioconda::cd-hit=4.8.1 bioconda::bwa-mem2=2.2.1"

    cpus 4
    memory "20 GB"

    input:
        tuple val(meta), path(reads), path(ffn), path(tsv)

    output:
        tuple val(meta), path("${meta.id}_clustered_catalog.fna*"), emit: index

    script:
    def catalog_fna = "${meta.id}_clustered_catalog.fna"
    """
    # Step 1: Cluster at 100%
    cd-hit-est -i ${ffn} -o ${catalog_fna} -c 1.0 -T ${task.cpus} -d 0 -M 0

    # Step 2: Index the gene catalogue
    bwa-mem2 index ${catalog_fna}
    """
}

process ALIGN_AND_QUANTIFY_READS {
    tag "${meta.id}"
    publishDir "$params.outdir/published/branch_2/03_gene_quantification/${meta.id}", mode: 'copy', pattern: "${meta.id}.idxstats.txt"
    afterScript "D=${params.outdir}/published/branch_2/03_gene_quantification/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} \$D/ 2>/dev/null || true"
    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.19.2"

    cpus 16
    memory "40 GB"

    input:
        tuple val(meta), path(reads), path(index_files)

    output:
        tuple val(meta), path("${meta.id}.idxstats.txt"), emit: idxstats

    script:
    // Find the .fna file in the index file list to use as the base for alignment
    def index_base = index_files.find { it.name.endsWith('.fna') }
    def sorted_bam = "${meta.id}.bam"
    """
    bwa-mem2 mem -t ${task.cpus} ${index_base} ${reads} | samtools view -b - | samtools sort --threads ${task.cpus} -o ${sorted_bam}
    samtools index --threads ${task.cpus} ${sorted_bam}
    samtools idxstats --threads ${task.cpus} ${sorted_bam} > ${meta.id}.idxstats.txt
    """
}

process CALCULATE_TPM_AND_ANNOTATE {
    tag "${meta.id}"
    publishDir "$params.outdir/published/branch_2/04_tpm_and_annotate/${meta.id}", mode: 'copy', pattern: "${meta.id}_gene_quantification.tsv"
    afterScript "D=${params.outdir}/published/branch_2/04_tpm_and_annotate/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} \$D/ 2>/dev/null || true"
    conda "pandas"

    cpus 1
    memory "10 GB"

    input:
        tuple val(meta), path(idxstats), path(tsv_file)

    output:
        tuple val(meta), path("${meta.id}_gene_quantification.tsv"), emit: quantification

    script:
    """
    tpm_annotate.py ${idxstats} ${tsv_file} ${meta.id}_gene_quantification.tsv
    """
}

// =======================================================
// --- BRANCH 3 PROCESSES ---

process MAP_FOR_BINNING {
    tag "${meta.id}"
    afterScript "D=${params.outdir}/published/branch_3/01_mapping/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} \$D/ 2>/dev/null || true"
    conda "bowtie2 samtools"

    cpus 30
    memory "60 GB"

    input:
        tuple val(meta), path(reads), path(contigs)

    output:
        tuple val(meta), path("${meta.id}.sorted.bam"), path(contigs), emit: bam

    script:
    def reads_args = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    bowtie2-build --threads ${task.cpus} -f ${contigs} ${meta.id}.assembly.idx
    bowtie2 -x ${meta.id}.assembly.idx ${reads_args} --threads ${task.cpus} | samtools view -u -b -F 4 | samtools sort --threads ${task.cpus} -o ${meta.id}.sorted.bam
    """
}

process METABAT2_BINNING {
    tag "${meta.id}"
    afterScript "D=${params.outdir}/published/branch_3/02_binning/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} \$D/ 2>/dev/null || true"
    conda "metabat2"

    cpus 8
    memory "20 GB"

    input:
        tuple val(meta), path(bam), path(contigs)

    output:
        tuple val(meta), path("${meta.id}_bins"), emit: bins

    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth ${meta.id}.depth.txt ${bam}
    mkdir -p ${meta.id}_bins
    metabat2 -i ${contigs} -a ${meta.id}.depth.txt -o ${meta.id}_bins/bin -t ${task.cpus} -m 1500 -v > metabat2.log 2>&1 \
        || grep -q "There were no large target contigs" metabat2.log || { cat metabat2.log; exit 1; }
    cat metabat2.log
    """
}

process CHECKM_QA {
    tag "${meta.id}"
    publishDir "$params.outdir/published/branch_3/03_checkm_qa/${meta.id}", mode: 'copy', pattern: "{quality_report.tsv,DIAMOND_RESULTS.tsv}"
    afterScript "D=${params.outdir}/published/branch_3/03_checkm_qa/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} checkm2.log \$D/ 2>/dev/null || true"
    conda "checkm2"

    cpus 12
    memory "50 GB"

    input:
        tuple val(meta), path(bins_dir)

    output:
        tuple val(meta), path("quality_report.tsv"), emit: checkm_summary
        tuple val(meta), path("DIAMOND_RESULTS.tsv"), emit: diamond_results

    script:
    """
    N_BINS=\$(find -L ${bins_dir} -maxdepth 1 -type f -name "*.fa" | wc -l)
    echo "\$(date -Iseconds) CheckM input ${bins_dir}: \$N_BINS bin(s)"
    if [ "\$N_BINS" -eq 0 ]; then
        echo "No bins in ${bins_dir}: MetaBAT2 made none, so CheckM2 does not run."
        printf 'Name\\tCompleteness\\tContamination\\n' > quality_report.tsv
        touch DIAMOND_RESULTS.tsv
        exit 0
    fi

    checkm2 predict --input ${bins_dir} --output-directory ${meta.id}_checkm2_out --threads ${task.cpus} -x fa
    mv ${meta.id}_checkm2_out/{quality_report.tsv,checkm2.log} .
    # With more than 250 bins, CheckM2 writes DIAMOND_RESULTS_0.tsv, _1.tsv and so on.
    cat ${meta.id}_checkm2_out/diamond_output/DIAMOND_RESULTS*.tsv > DIAMOND_RESULTS.tsv
    """
}

process FILTER_AND_ANNOTATE {
    tag "${meta.id}"
    publishDir "$params.outdir/published/branch_3/04_filter_and_annotate/${meta.id}", mode: 'copy', pattern: '{*.tsv,hq_bins}' // move mode does not work and throws java.nio.file.DirectoryNotEmptyException for hq_bins
    afterScript "D=${params.outdir}/published/branch_3/04_filter_and_annotate/${meta.id}/logs; mkdir -p \$D && cp -f ${task.workDir}/.command.{sh,err,out,log} gtdbtk.log gtdbtk.json \$D/ 2>/dev/null || true"
    conda "python=3.12 gtdbtk pandas"

    cpus { task.attempt == 1 ? 20 : 35 }
    memory { task.attempt == 1 ? "80 GB" : "140 GB" }

    errorStrategy { task.attempt <=1 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
        tuple val(meta), path(bins_dir), path(checkm_summary)

    output:
        tuple val(meta), path("hq_bins"), emit: hq_bins
        tuple val(meta), path("gtdbtk.bac120.summary.tsv"), emit: gtdbtk_summary, optional: true
        tuple val(meta), path("NO_HQ_BINS_FOUND.tsv"), emit: NO_HQ_BINS_FOUND, optional: true

    script:
    def min_completeness = 90
    def max_contamination = 5
    """
    filter_bins.py ${checkm_summary} ${bins_dir} hq_bins ${min_completeness} ${max_contamination}
    
    if [ -n "\$(find hq_bins -name '*.fa')" ]; then
        gtdbtk classify_wf \\
            --genome_dir hq_bins \\
            --out_dir ${meta.id}_gtdbtk_out \\
            --cpus ${task.cpus} \\
            --extension fa
        mv ${meta.id}_gtdbtk_out/{gtdbtk.log,gtdbtk.json} .
        if [ -e ${meta.id}_gtdbtk_out/classify/gtdbtk.bac120.summary.tsv ]; then
            mv ${meta.id}_gtdbtk_out/classify/gtdbtk.bac120.summary.tsv .
        fi
    else
        touch NO_HQ_BINS_FOUND.tsv
    fi
    """
}
