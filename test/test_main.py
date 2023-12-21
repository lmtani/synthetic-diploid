from make_diploid import main


def test_main_should_return_correct_number_of_variants_het_hom(vcf_1, vcf_2):
    # call
    total, hom, het = main(vcf_1, vcf_2, "new_sample_name", "temp-output.vcf")

    # assert
    assert len(total) == 3
    assert hom == 1
    assert het == 2

