import pytest
from make_diploid import load_variants


def test_load_variants_should_raise_exception_if_more_than_one_sample(vcf_reader_multisample):
    with pytest.raises(ValueError) as e:
        load_variants(vcf_reader_multisample)
    assert str(e.value) == "More than one sample found in the same VCF file. Only single sample VCFs are supported"


def test_load_variants_should_return_sample_name_and_variants(vcf_reader_2):
    sample_name, variants = load_variants(vcf_reader_2)

    assert sample_name == "SSC39-MAT-2-DNA"
    assert len(variants) == 2

    a = [i for i in vcf_reader_2.fetch("chr1", 113456, 113457)]
    b = [i for i in vcf_reader_2.fetch("chr1", 113466, 113467)]
    expected = {
        ('chr1', 113457, 'G', 'GC'): a,
        ('chr1', 113467, 'G', 'T'): b
    }
    key1 = ('chr1', 113457, 'G', 'GC')
    assert variants[key1] == expected[key1]
    key2 = ('chr1', 113467, 'G', 'T')
    assert variants[key2] == expected[key2]
