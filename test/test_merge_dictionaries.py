from make_diploid import merge_dictionaries


def test_merge_dictionaries_should_merge_two_dictionaries():
    # arrange
    records_1 = {
        ('chr1', 113457, 'G', 'GC'): ["a"],
        ('chr1', 113467, 'G', 'T'): ["b"]
    }
    records_2 = {
        ('chr1', 113457, 'G', 'GC'): ["c"],
        ('chr1', 113467, 'G', 'T'): ["d"]
    }

    # call
    result = merge_dictionaries(records_1, records_2)

    # assert
    assert len(result) == 2
    assert result[('chr1', 113457, 'G', 'GC')] == ["a", "c"]
    assert result[('chr1', 113467, 'G', 'T')] == ["b", "d"]


def test_merge_dictionaries_should_merge_two_dictionaries_with_different_keys():
    # arrange
    records_1 = {
        ('chr1', 113457, 'G', 'GC'): ["a"],
        ('chr1', 113467, 'G', 'T'): ["b"]
    }
    records_2 = {
        ('chr1', 113458, 'G', 'GC'): ["c"],
        ('chr1', 113468, 'G', 'T'): ["d"]
    }

    # call
    result = merge_dictionaries(records_1, records_2)

    # assert
    assert len(result) == 4
    assert result[('chr1', 113457, 'G', 'GC')] == ["a"]
    assert result[('chr1', 113467, 'G', 'T')] == ["b"]
    assert result[('chr1', 113458, 'G', 'GC')] == ["c"]
    assert result[('chr1', 113468, 'G', 'T')] == ["d"]


def test_merge_dictionaries_should_merge_two_dictionaries_with_equal_and_different_keys():
    # arrange
    records_1 = {
        ('chr1', 113457, 'G', 'GC'): ["a"],
        ('chr1', 113467, 'G', 'T'): ["b"]
    }
    records_2 = {
        ('chr1', 113457, 'G', 'GC'): ["c"],
        ('chr1', 113468, 'G', 'T'): ["d"]
    }

    # call
    result = merge_dictionaries(records_1, records_2)

    # assert
    assert len(result) == 3
    assert result[('chr1', 113457, 'G', 'GC')] == ["a", "c"]
    assert result[('chr1', 113467, 'G', 'T')] == ["b"]
    assert result[('chr1', 113468, 'G', 'T')] == ["d"]
