import pytest
from transcript_to_bed import transcript_to_bed, Bed


@pytest.fixture
def transcript():
    """Fake transcript object"""
    return {
            'seqname': 'chr1',
            'start': '1',
            'end': '10',
            'attribute': {'transcript_id': 'ENST00000375687'}
    }


def test_transcript_to_bed_positions(transcript):
    bed = transcript_to_bed(transcript)
    assert bed.start == 0
    assert bed.end == 9


def test_transcript_to_bed_positions(transcript):
    bed = transcript_to_bed(transcript)
    assert bed.start == 0
    assert bed.end == 9


def test_transcript_to_bed_chrom(transcript):
    bed = transcript_to_bed(transcript)
    assert bed.chrom == 'chr1'


def test_transcript_to_bed_name(transcript):
    bed = transcript_to_bed(transcript)
    assert bed.name == 'ENST00000375687'


def test_print_bed(transcript):
    bed = transcript_to_bed(transcript)
    assert str(bed) == 'chr1\t0\t9\tENST00000375687'


def test_order_bed_records():
    bed1 = Bed('chr1', 0, 10, 'ENTS01')
    bed2 = Bed('chr1', 10, 20, 'ENTS01')

    l = [bed2, bed1]

    assert sorted(l) == [bed1, bed2]

