from transcript_cov import Depth, get_exon_positions, get_transcript_positions


def test_create_depth():
    d = Depth('chr1', 1, 10)
    assert True


def test_get_exon_positions_pos_strand():
    exon = {'start': '100', 'end': '105', 'strand': '+'}
    expected = [100, 101, 102, 103, 104, 105]

    assert list(get_exon_positions(exon)) == expected


def test_get_exon_positions_neg_strand():
    exon = {'start': '100', 'end': '105', 'strand': '-'}
    expected = [105, 104, 103, 102, 101, 100]

    assert list(get_exon_positions(exon)) == expected


def test_get_transcript_positions_pos_strand():
    # A transcript is simply a list of exons
    exon1 = {'start': '100', 'end': '105', 'strand': '+'}
    exon2 = {'start': '110', 'end': '113', 'strand': '+'}
    transcript = [exon1, exon2]

    expected = [100, 101, 102, 103, 104, 105, 110, 111, 112, 113]

    assert list(get_transcript_positions(transcript)) == expected


def test_get_transcript_positions_neg_strand():
    # A transcript is simply a list of exons
    # For negative strand exons:
    # 1. The exons within the transcript are in reverse order on the chromosome
    # 2. The start positions is always before the end positions on the
    #    chromosome
    #
    # To get the positions for the transcript (in order), the positions should
    # also be *decreasing* relative to the chromsome position
    exon1 = {'start': '110', 'end': '113', 'strand': '-'}
    exon2 = {'start': '100', 'end': '105', 'strand': '-'}
    transcript = [exon1, exon2]

    expected = [113, 112, 111, 110, 105, 104, 103, 102, 101, 100]

    assert list(get_transcript_positions(transcript)) == expected
