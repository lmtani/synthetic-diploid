#!/usr/bin/env python

import pysam
import sys
import coloredlogs, logging
coloredlogs.install(
    level=logging.INFO,
    fmt='%(asctime)s,%(msecs)03d %(hostname)s %(levelname)s %(message)s'
)

SAMPLE_1 = 0
SAMPLE_2 = 1

VCF_1 = 0
VCF_2 = 1


def load_variants(vcf_reader):
    result = {}
    logging.info(f"📦 Loading variants from {vcf_reader.filename}")
    sample_names = [s for s in vcf_reader.header.samples]
    logging.info(f"  Sample: {' '.join(sample_names)}")
    if len(sample_names) > 1:
        raise ValueError("More than one sample found in the same VCF file. Only single sample VCFs are supported")
    sample_name = sample_names[0]

    indels = 0
    snvs = 0
    heterozigous = 0
    multiallelicsite = 0
    ref_hom = 0
    total = 0
    for r in vcf_reader:
        total += 1
        if len(r.alts) > 1:
            multiallelicsite += 1
            continue

        gt = r.samples[sample_name]["GT"]
        if gt == (0,) or gt == (None, ):
            ref_hom += 1
            continue

        if len(gt) != 1:
            if gt == (0, 0) or gt == (None, None):
                ref_hom += 1
                continue
            if gt[0] != gt[1]:
                heterozigous += 1
                continue

        if len(r.ref) > 1 or len(r.alts[0]) > 1:
            indels += 1
        else:
            snvs += 1

        key = (r.chrom, r.pos, r.ref, r.alts[0])
        result.setdefault(key, []).append(r)

    if not result:
        raise ValueError("No variants found in the VCF file")

    logging.info(f"  Total variants: {total}")
    logging.info(f"    - SNVs: {snvs}")
    logging.info(f"    - Indels: {indels}")
    logging.info(f"    - Reference homozygous (dropped): {ref_hom}")

    logging.info(f"  Dropped because more than one alternative allele: {multiallelicsite}")
    logging.info(f"  Dropped because heterozigous genotypes: {heterozigous}")
    return sample_name, result


def merge_dictionaries(records_1, records_2):
    result = {}
    for key in records_1:
        if key in records_2:
            result[key] = records_1[key] + records_2[key]
        else:
            result[key] = records_1[key]

    for key in records_2:
        if key not in records_1:
            result[key] = records_2[key]

    return result


def create_vcf_header(vcf_reader, sample_name):
    new_header = pysam.VariantHeader()
    new_header.add_sample(sample_name)
    new_header.add_line('##FORMAT=<ID=OQ,Number=R,Type=Float,Description="Original genotyping quality">')
    new_header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    new_header.add_line('##FORMAT=<ID=SD,Number=R,Type=Integer,Description="Original genotyping depth">')
    new_header.add_line('##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Samples where the variant was found">')
    for contig in vcf_reader.header.contigs.values():
        new_header.add_line(f"##contig=<ID={contig.name},length={contig.length}>")
    return new_header


def extract_vcf_format_field(sample1, sample2, vars_in_position, field: str):
    val1 = vars_in_position[SAMPLE_1].samples[sample1][field]
    val2 = vars_in_position[SAMPLE_2].samples[sample2][field]
    return val1, val2


def parse_fields(vars_in_position: list, phase: dict):
    if len(vars_in_position) == 2:
        first_sample = 0  # We don't support more than one sample per VCF
        sample1 = vars_in_position[VCF_1].samples.keys()[first_sample]
        sample2 = vars_in_position[VCF_2].samples.keys()[first_sample]

        hap1, hap2 = extract_vcf_format_field(sample1, sample2, vars_in_position, "GT")
        oq1, oq2 = extract_vcf_format_field(sample1, sample2, vars_in_position, "GQ")
        dp1, dp2 = extract_vcf_format_field(sample1, sample2, vars_in_position, "DP")
        haploid_sanity_check(hap1, hap2)

        # We keep only with haploid genotypes or homozigous diploid genotypes, so we can get the first one
        first_gt_field = 0
        gt = phase["both"](hap1[first_gt_field], hap2[first_gt_field])
        oq = phase["both"](oq1, oq2)
        dp = phase["both"](dp1, dp2)
        joined_samples = ",".join(vars_in_position[0].samples.keys() + vars_in_position[1].samples.keys())
    elif len(vars_in_position) == 1:
        sample = vars_in_position[0].samples.keys()[0]
        gt = phase[sample](vars_in_position[0].samples[sample]["GT"][0])
        oq = phase[sample](vars_in_position[0].samples[sample]["GQ"])
        dp = phase[sample](vars_in_position[0].samples[sample]["DP"])
        joined_samples = ",".join(vars_in_position[-1].samples.keys())
    else:
        raise ValueError("More than 2 variants found in the same position")
    return gt, oq, dp, joined_samples


def haploid_sanity_check(hap1, hap2):
    if len(hap1) > 1:
        if hap1[0] != hap1[1]:
            logging.warning("Looks like this VCF is not represented as haploid.")
            raise ValueError("We don't support heterozigoud diploid VCFs")
    if len(hap2) > 1:
        if hap2[0] != hap2[1]:
            logging.warning("Looks like this VCF is not represented as haploid.")
            raise ValueError("We don't support heterozigoud diploid VCFs")


def create_vcf_record(new_header, original_variant, oq, dp, gt, joined_samples, sample_name):
    logging.debug("Creating VCF record")
    logging.debug("GT: %s", gt)
    logging.debug("OQ: %s", oq)
    logging.debug("DP: %s", dp)
    logging.debug("Samples: %s", joined_samples)
    r = new_header.new_record()
    r.contig = str(original_variant.chrom)
    r.pos = int(original_variant.pos)
    r.id = original_variant.id
    r.ref = original_variant.ref
    r.alts = original_variant.alts
    new_qualities = min(oq)
    r.qual = max(oq) if gt != (1, 1) else new_qualities  # if variant only in one sample, we use the original quality because zero is attributed to the other sample
    r.info["SAMPLES"] = joined_samples
    r.samples[sample_name]["OQ"] = oq
    r.samples[sample_name]['SD'] = dp
    r.samples[sample_name]['GT'] = gt
    r.samples[sample_name].phased = True
    return r


def main(path_1, path_2, sample_name, output_name):
    vcf = pysam.VariantFile(path_1)
    vcf2 = pysam.VariantFile(path_2)
    sample_name_1, haplotype_1 = load_variants(vcf)
    sample_name_2, haplotype_2 = load_variants(vcf2)
    variants = merge_dictionaries(haplotype_1, haplotype_2)
    header = create_vcf_header(vcf, sample_name)
    haplotype_phase = {
        sample_name_1: lambda n: (n, 0),
        sample_name_2: lambda n: (0, n),
        "both": lambda n, m: (n, m)
    }
    heterozigous = 0
    homozigous = 0
    # sort variants by contig and position
    variants = dict(sorted(variants.items(), key=lambda item: (item[0][0], item[0][1])))
    with pysam.VariantFile(output_name, 'w', header=header) as vcf_out:
        for var in variants:
            original_record = variants[var]
            genotype, original_qualities, read_depth, samples = parse_fields(variants[var], haplotype_phase)

            if genotype == (1, 1):
                homozigous += 1
            elif genotype == (0, 1) or genotype == (1, 0):
                heterozigous += 1
            else:
                logging.error("Unexpected genotype at position %s %s: %s",  original_record[0].contig, original_record[0].pos, genotype)
                raise ValueError("Unexpected genotype")

            # we junst want positions, ref, alt, so we can get the first original record
            record = create_vcf_record(
                header, original_record[0], original_qualities, read_depth, genotype, samples, sample_name
            )

            vcf_out.write(record)
    logging.info(f"🖫 Wrote {len(variants)} variants to {output_name}")
    logging.info(f"  Homozigous: {homozigous}")
    logging.info(f"  Heterozigous: {heterozigous}")

    return variants, homozigous, heterozigous


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python make_diploid.py <vcf_file> <vcf_file>")
        sys.exit(1)

    vcf_path_1 = sys.argv[1]
    vcf_path_2 = sys.argv[2]
    synthetic_sample_name = "synthetic-diploid"
    outfile = "diploid_genotypes.vcf"

    main(vcf_path_1, vcf_path_2, synthetic_sample_name, outfile)
