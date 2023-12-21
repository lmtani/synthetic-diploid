import pytest
import pysam


@pytest.fixture()
def vcf_1():
    return "test/data/SSC39-MAT-1.deepvariant.vcf.gz"


@pytest.fixture()
def vcf_2():
    return "test/data/SSC39-MAT-2.clair3.vcf.gz"


@pytest.fixture()
def vcf_3():
    return "test/data/multisample.vcf.gz"


@pytest.fixture()
def vcf_reader_1(vcf_1):
    return pysam.VariantFile(vcf_1)


@pytest.fixture()
def vcf_reader_2(vcf_2):
    return pysam.VariantFile(vcf_2)


@pytest.fixture()
def vcf_reader_multisample(vcf_3):
    return pysam.VariantFile(vcf_3)


@pytest.fixture()
def single_variant_record(vcf_reader_1):
    return [i for i in vcf_reader_1.fetch("chr1", 113456, 113457)][0]
