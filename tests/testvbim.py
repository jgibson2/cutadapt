from __future__ import print_function, division, absolute_import
from nose.tools import raises, assert_raises

from cutadapt.seqio import Sequence
from cutadapt.vbim import (VBIMSeq, VBIMMatch, ANYWHERE)


def test_sequence():
    seq = VBIMSeq(
        sequence='CCACCATGGATTACAAGGATGACGACGATAAGAATTCTT',
        where=ANYWHERE,
        max_error_rate=0.1,
        read_wildcards=False,
        adapter_wildcards=False)
    read = Sequence(name="test1", sequence='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCACCATGGATTACAAGGATGACGACGATAAGAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT')
    m = seq.match_to(read)
    assert m.trimmed().sequence == 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTT'