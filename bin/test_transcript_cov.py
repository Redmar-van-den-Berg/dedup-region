import pytest

from transcript_cov import Depth
from transcript_cov import get_exon_positions, get_transcript_positions
from transcript_cov import get_transcript_depth
from transcript_cov import get_avg_exon_depth
from transcript_cov import get_avg_transcript_depth


@pytest.fixture
def depths():
    # Create list of depths
    # (pos, depth)
    depths = [
            (100, 5),
            (101, 4),
            (102, 3),
            (103, 4),
            (104, 5),
            (105, 5),
            (106, 5),
            (107, 6),
            (108, 3),
            (109, 3),
            (110, 3)
    ]
    # Just prepend chr1 here, since it is constant
    return [Depth('chr1', *x) for x in depths]


def test_create_depth():
    d = Depth('chr1', 1, 10)
    assert True


def test_get_exon_positions_pos_strand():
    exon = {'start': '100', 'end': '105', 'strand': '+', 'seqname': 'chr1'}
    expected = [100, 101, 102, 103, 104, 105]

    assert list(get_exon_positions(exon)) == expected


def test_get_exon_positions_neg_strand():
    exon = {'start': '100', 'end': '105', 'strand': '-', 'seqname': 'chr1'}
    expected = [105, 104, 103, 102, 101, 100]

    assert list(get_exon_positions(exon)) == expected


def test_get_transcript_positions_pos_strand():
    # A transcript is simply a list of exons
    exon1 = {'start': '100', 'end': '105', 'strand': '+', 'seqname': 'chr1'}
    exon2 = {'start': '110', 'end': '113', 'strand': '+', 'seqname': 'chr1'}
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
    exon1 = {'start': '110', 'end': '113', 'strand': '-', 'seqname': 'chr1'}
    exon2 = {'start': '100', 'end': '105', 'strand': '-', 'seqname': 'chr1'}
    transcript = [exon1, exon2]

    expected = [113, 112, 111, 110, 105, 104, 103, 102, 101, 100]

    assert list(get_transcript_positions(transcript)) == expected


def test_get_transcript_depth_pos_strand(depths):
    # A transcript is simply a list of exons
    exon1 = {'start': '100', 'end': '103', 'strand': '+', 'seqname': 'chr1'}
    exon2 = {'start': '107', 'end': '108', 'strand': '+', 'seqname': 'chr1'}
    transcript = [exon1, exon2]

    expected = [
            Depth('chr1', 100, 5),
            Depth('chr1', 101, 4),
            Depth('chr1', 102, 3),
            Depth('chr1', 103, 4),
            Depth('chr1', 107, 6),
            Depth('chr1', 108, 3)
    ]

    assert list(get_transcript_depth(depths, transcript)) == expected


def test_get_transcript_depth_neg_strand(depths):
    # A transcript is simply a list of exons
    exon1 = {'start': '107', 'end': '108', 'strand': '-', 'seqname': 'chr1'}
    exon2 = {'start': '100', 'end': '103', 'strand': '-', 'seqname': 'chr1'}
    transcript = [exon1, exon2]

    expected = [
            Depth('chr1', 108, 3),
            Depth('chr1', 107, 6),
            Depth('chr1', 103, 4),
            Depth('chr1', 102, 3),
            Depth('chr1', 101, 4),
            Depth('chr1', 100, 5)
    ]

    assert list(get_transcript_depth(depths, transcript)) == expected


def test_get_transcript_depth_no_depth():
    # A transcript is simply a list of exons
    exon1 = {'start': '107', 'end': '108', 'strand': '-', 'seqname': 'chr1'}
    exon2 = {'start': '100', 'end': '103', 'strand': '-', 'seqname': 'chr1'}
    transcript = [exon1, exon2]

    expected = [
            Depth('chr1', 108, 0),
            Depth('chr1', 107, 0),
            Depth('chr1', 103, 0),
            Depth('chr1', 102, 0),
            Depth('chr1', 101, 0),
            Depth('chr1', 100, 0)
    ]

    assert list(get_transcript_depth([], transcript)) == expected


def test_get_avg_exon_depth(depths):
    # A transcript is simply a list of exons
    exon1 = {'start': '107', 'end': '108', 'strand': '-', 'seqname': 'chr1'}
    exon2 = {'start': '107', 'end': '108', 'strand': '+', 'seqname': 'chr1'}

    # The average depth is the same for + and - strand exons
    assert get_avg_exon_depth(depths, exon1) == 4.5
    assert get_avg_exon_depth(depths, exon2) == 4.5


def test_get_avg_transcript_depth(depths):
    # A transcript is simply a list of exons
    exon1 = {'start': '107', 'end': '108', 'strand': '-', 'seqname': 'chr1'}
    exon2 = {'start': '100', 'end': '103', 'strand': '-', 'seqname': 'chr1'}
    transcript = [exon1, exon2]

    expected = [4.5, 4]

    assert get_avg_transcript_depth(depths, transcript) == expected


def test_avg_exon_depth_no_depth():
    """Test that we receive all zero coverage when no depth is defined

    This will happen sometimes with the output of samtools depth
    """
    exon1 = {'start': '107', 'end': '108', 'strand': '-', 'seqname': 'chr1'}
    assert get_avg_exon_depth([], exon1) == 0

def test_get_avg_transcript_depth_no_depth():
    """Test that we receive all zero coverage when no depth is defined

    This will happen sometimes with the output of samtools depth
    """
    # A transcript is simply a list of exons
    exon1 = {'start': '107', 'end': '108', 'strand': '-', 'seqname': 'chr1'}
    exon2 = {'start': '100', 'end': '103', 'strand': '-', 'seqname': 'chr1'}
    transcript = [exon1, exon2]

    expected = [0, 0]

    assert get_avg_transcript_depth(list(), transcript) == expected
