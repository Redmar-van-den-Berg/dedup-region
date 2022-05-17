from transcript_cov import Depth, get_exon_positions


def test_create_depth():
    d = Depth('chr1', 1, 10)
    assert True


def test_get_exon_positions():
    exon = {'start': '100', 'end': '105', 'strand': '+'}
    expected = [100, 101, 102, 103, 104, 105]

    assert list(get_exon_positions(exon)) == expected
