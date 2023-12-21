from make_diploid import create_vcf_record, create_vcf_header


def test_create_vcf_record_should_return_a_vcf_record(vcf_reader_1, single_variant_record):
    # arrange
    oq = (0, 1)
    dp = (0, 20)
    gt = (0, 1)
    sample_name = "diploid"
    header = create_vcf_header(vcf_reader_1, sample_name)
    joined_samples = "one_sample"

    # call
    result = create_vcf_record(header, single_variant_record, oq, dp, gt, joined_samples, sample_name)

    # assert
    assert result.pos == single_variant_record.pos
    assert result.ref == single_variant_record.ref
    assert result.alts == single_variant_record.alts
    assert result.qual == 1
    assert result.samples[sample_name]["OQ"] == oq
    assert result.samples[sample_name]["SD"] == dp
    assert result.samples[sample_name]["GT"] == gt
    assert result.info["SAMPLES"] == (joined_samples, )
