# Synthetic Diploid Maker

## Description

This script takes two VCF files (one for each haplotype) and creates a new VCF file with a synthetic diploid genome.

### Limitations

* Does not process multisample VCF files
* Removes reference calls (e.g.: 0/0)
* Assume that both VCF files follow the same pattern of variant representation (e.g.: left-aligned, normalized, etc.)
* If a VCF file has diploid representation (e.g.: 0/1), the script will discard positions with heterozygous variants

## Usage

It requires python 3.7 or higher.


```bash
pip install requirements.txt
python synthetic_diploid_maker.py path-to-vcf-1 path-to-vcf-2
```

## Output

The output is a VCF file (diploid_genotypes.vcf) with the synthetic diploid genome.

This will include the following INFO field:
- SAMPLES: Samples where the variant was found

And the following FORMAT fields:
- GT: Genotype
- SD (Sample Depth): Depth of coverage for the sample
- OQ (Original Quality): Quality of the variant in the original VCF file
