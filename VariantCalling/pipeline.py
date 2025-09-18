#!/usr/bin/env python3

import os
import argparse
import subprocess
import json

def run(cmd, log_file):
    print(f"[RUNNING] {cmd}")
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, shell=True, stdout=log, stderr=log)
        if result.returncode != 0:
            print(f"[ERROR] Komanda nije uspela: {cmd}")
            exit(1)

def file_exists_and_not_empty(filepath):
    return os.path.exists(filepath) and os.path.getsize(filepath) > 0

def bam_is_valid(bam_file):
    try:
        result = subprocess.run(f"samtools quickcheck {bam_file}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.returncode == 0
    except Exception:
        return False

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq1", required=True)
    parser.add_argument("--fastq2", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. FASTQC
    fastqc_dir = os.path.join(args.outdir, "01_fastqc")
    os.makedirs(fastqc_dir, exist_ok=True)
    fastqc_log = os.path.join(fastqc_dir, "fastqc.log")
    if not file_exists_and_not_empty(fastqc_log):
        run(f"fastqc {args.fastq1} {args.fastq2} -o {fastqc_dir}", fastqc_log)
    else:
        print(f"[SKIPPING] FastQC već izvršen.")

    # 2. BWA MEM + SAMTOOLS SORT
    mapping_dir = os.path.join(args.outdir, "02_mapping")
    os.makedirs(mapping_dir, exist_ok=True)
    sam_file = os.path.join(mapping_dir, "aligned.sam")
    bam_file = os.path.join(mapping_dir, "aligned.bam")
    sorted_bam = os.path.join(mapping_dir, "sorted.bam")

    bwa_log = os.path.join(mapping_dir, "bwa.log")
    if not file_exists_and_not_empty(sam_file):
        run(f"bwa mem -t 2 -K 10000000 -R '@RG\\tID:sample1\\tSM:sample1\\tPL:ILLUMINA' {args.reference} {args.fastq1} {args.fastq2} > {sam_file}", bwa_log)
    else:
        print(f"[SKIPPING] {sam_file} već postoji.")

    view_log = os.path.join(mapping_dir, "samtools_view.log")
    if not file_exists_and_not_empty(bam_file):
        run(f"samtools view -Sb {sam_file} > {bam_file}", view_log)
    else:
        print(f"[SKIPPING] {bam_file} već postoji.")

    sort_log = os.path.join(mapping_dir, "samtools_sort.log")
    if not file_exists_and_not_empty(sorted_bam) or not bam_is_valid(sorted_bam):
        run(f"samtools sort {bam_file} -o {sorted_bam}", sort_log)
    else:
        print(f"[SKIPPING] {sorted_bam} već postoji i validan je.")

    index_log = os.path.join(mapping_dir, "samtools_index.log")
    if not file_exists_and_not_empty(sorted_bam + ".bai"):
        run(f"samtools index {sorted_bam}", index_log)
    else:
        print(f"[SKIPPING] BAM indeks već postoji.")

    # 3. MarkDuplicates (Picard)
    dedup_dir = os.path.join(args.outdir, "03_dedup")
    os.makedirs(dedup_dir, exist_ok=True)
    dedup_bam = os.path.join(dedup_dir, "dedup.bam")
    metrics_file = os.path.join(dedup_dir, "dup_metrics.txt")
    picard_log = os.path.join(dedup_dir, "markduplicates.log")

    if not file_exists_and_not_empty(dedup_bam):
        picard_cmd = (
            f"java -jar ~/tools/picard.jar MarkDuplicates "
            f"I={sorted_bam} "
            f"O={dedup_bam} "
            f"M={metrics_file} "
            f"REMOVE_DUPLICATES=false "
            f"CREATE_INDEX=true"
        )
        run(picard_cmd, picard_log)
    else:
        print(f"[SKIPPING] {dedup_bam} već postoji.")

    # 4. BQSR (Base Quality Score Recalibration)
    bqsr_dir = os.path.join(args.outdir, "04_bqsr")
    os.makedirs(bqsr_dir, exist_ok=True)
    recal_table = os.path.join(bqsr_dir, "recal_data.table")
    recal_bam = os.path.join(bqsr_dir, "recalibrated.bam")
    bqsr_base_log = os.path.join(bqsr_dir, "bqsr_base_recalibrator.log")
    bqsr_apply_log = os.path.join(bqsr_dir, "bqsr_apply.log")

    known_sites = os.path.join(os.path.expanduser("~"), "projekat2", "resources", "known_sites_fixed.vcf")

    if not file_exists_and_not_empty(recal_table):
        run(f"gatk BaseRecalibrator -R {args.reference} -I {dedup_bam} --known-sites {known_sites} -O {recal_table}", bqsr_base_log)
    else:
        print(f"[SKIPPING] {recal_table} već postoji.")

    if not file_exists_and_not_empty(recal_bam):
        run(f"gatk ApplyBQSR -R {args.reference} -I {dedup_bam} --bqsr-recal-file {recal_table} -O {recal_bam}", bqsr_apply_log)
    else:
        print(f"[SKIPPING] {recal_bam} već postoji.")

    # 5. HaplotypeCaller
    haplotype_dir = os.path.join(args.outdir, "05_haplotypecaller")
    os.makedirs(haplotype_dir, exist_ok=True)
    raw_vcf = os.path.join(haplotype_dir, "raw_variants.vcf")
    haplotype_log = os.path.join(haplotype_dir, "haplotypecaller.log")

    if not file_exists_and_not_empty(raw_vcf):
        run(f"gatk HaplotypeCaller -R {args.reference} -I {recal_bam} -O {raw_vcf} --native-pair-hmm-threads 2", haplotype_log)
    else:
        print(f"[SKIPPING] {raw_vcf} već postoji.")

    # 6. Variant Filtering
    filtering_dir = os.path.join(args.outdir, "06_filtering")
    os.makedirs(filtering_dir, exist_ok=True)
    filtered_vcf = os.path.join(filtering_dir, "filtered_variants.vcf")
    filtering_log = os.path.join(filtering_dir, "variant_filtering.log")

    if not file_exists_and_not_empty(filtered_vcf):
        run(f"gatk VariantFiltration -R {args.reference} -V {raw_vcf} --filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0\" --filter-name \"basic_snp_filter\" -O {filtered_vcf}", filtering_log)
    else:
        print(f"[SKIPPING] {filtered_vcf} već postoji.")

    # 7. JSON zapis rezultata
    results = {
        "fastqc": {
            "log": fastqc_log
        },
        "mapping": {
            "sam": sam_file,
            "bam": bam_file,
            "sorted_bam": sorted_bam,
            "log": bwa_log
        },
        "deduplication": {
            "bam": dedup_bam,
            "metrics": metrics_file,
            "log": picard_log
        },
        "bqsr": {
            "table": recal_table,
            "bam": recal_bam,
            "log": bqsr_apply_log
        },
        "haplotypecaller": {
            "vcf": raw_vcf,
            "log": haplotype_log
        },
        "filtering": {
            "filtered_vcf": filtered_vcf,
            "log": filtering_log
        }
    }

    with open(os.path.join(args.outdir, "results.json"), "w") as f:
        json.dump(results, f, indent=4)

if __name__ == "__main__":
    main()
