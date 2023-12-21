import pysam
import sys
import coloredlogs, logging
logging.basicConfig(level=logging.INFO)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
coloredlogs.install(level=logging.INFO)

SAMPLE_1 = 0
SAMPLE_2 = 1

VCF_1 = 0
VCF_2 = 1


def load_variants(vcf_reader):
    result = {}
    logging.info(f"Loading variants from {vcf_reader.filename}")
    sample_names = [s for s in vcf_reader.header.samples]
    logging.info(f"Samples: {' '.join(sample_names)}")
    if len(sample_names) > 1:
        raise ValueError("More than one sample found in the same VCF file. Only single sample VCFs are supported")

    diploid_warning = False
    for r in vcf_reader:
        variant_key = f"{r.chrom}:{r.pos} {r.ref} {r.alts}"
        if len(r.alts) > 1:
            logging.warning("Skipping this variant because it has more than one alternative allele: %s", variant_key)
            continue

        gt = r.samples[sample_names[0]]["GT"]
        if len(gt) != 1:
            if not diploid_warning:
                logging.warning("Looks like this VCF is not represented as haploid.")
                diploid_warning = True

            if gt[0] != gt[1]:
                logging.warning(f"Skipping this variant because it can't be represented as haploid: {variant_key}")
                continue

        key = (r.chrom, r.pos, r.ref, r.alts[0])
        result.setdefault(key, []).append(r)
    return result


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


def parse_fields(vars_in_position: list):
    if len(vars_in_position) == 2:
        first_sample = 0  # We don't support more than one sample per VCF
        sample1 = vars_in_position[VCF_1].samples.keys()[first_sample]
        sample2 = vars_in_position[VCF_2].samples.keys()[first_sample]

        hap1, hap2 = extract_vcf_format_field(sample1, sample2, vars_in_position, "GT")
        oq1, oq2 = extract_vcf_format_field(sample1, sample2, vars_in_position, "GQ")
        dp1, dp2 = extract_vcf_format_field(sample1, sample2, vars_in_position, "DP")
        haploid_sanity_check(hap1, hap2)

        first_gt_field = 0  # We keep only with haploid genotypes or homozigous diploid genotypes, so we can get the first one
        gt = (hap1[first_gt_field], hap2[first_gt_field])
        oq = (oq1, oq2)
        dp = (dp1, dp2)
        joined_samples = ",".join(vars_in_position[0].samples.keys() + vars_in_position[1].samples.keys())
    elif len(vars_in_position) == 1:
        logging.info(
            f"Only one vcf has a variant in this position: {vars_in_position[0].samples.keys()[0]} - {vars_in_position[0].pos}"
        )

        sample = vars_in_position[0].samples.keys()[0]
        gt = (0, vars_in_position[0].samples[sample]["GT"][0])
        oq = (0, vars_in_position[0].samples[sample]["GQ"])
        dp = (0, vars_in_position[0].samples[sample]["DP"])
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


def create_vcf_record(original_variant, new_qual, oq, dp, gt, joined_samples):
    logging.debug("Creating VCF record")
    logging.debug("GT: %s", gt)
    logging.debug("OQ: %s", oq)
    logging.debug("DP: %s", dp)
    logging.debug("Samples: %s", joined_samples)
    r = header.new_record()
    r.contig = str(original_variant.chrom)
    r.pos = int(original_variant.pos)
    r.id = original_variant.id
    r.ref = original_variant.ref
    r.alts = original_variant.alts
    r.qual = new_qual
    r.info["SAMPLES"] = joined_samples
    r.samples[synthetic_sample_name]["OQ"] = oq
    r.samples[synthetic_sample_name]['SD'] = dp
    r.samples[synthetic_sample_name]['GT'] = gt
    r.samples[synthetic_sample_name].phased = False
    return r


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python make_diploid.py <vcf_file> <vcf_file>")
        sys.exit(1)

    vcf_path_1 = sys.argv[1]
    vcf_path_2 = sys.argv[2]
    synthetic_sample_name = "synthetic-diploid"
    outfile = "diploid_genotypes.vcf"

    vcf = pysam.VariantFile(vcf_path_1)
    vcf2 = pysam.VariantFile(vcf_path_2)

    haplotype_1 = load_variants(vcf)
    haplotype_2 = load_variants(vcf2)
    variants = merge_dictionaries(haplotype_1, haplotype_2)
    header = create_vcf_header(vcf, synthetic_sample_name)

    with pysam.VariantFile(outfile, 'w', header=header) as vcf_out:
        for var in variants:
            original_record = variants[var]
            genotype, original_qualities, read_depth, samples = parse_fields(variants[var])

            new_qualities = min(original_qualities)
            # we junst want positions, ref, alt, so we can get the first original record
            record = create_vcf_record(
                original_record[0], new_qualities, original_qualities, read_depth, genotype, samples
            )

            vcf_out.write(record)
