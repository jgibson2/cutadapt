# coding: utf-8
"""
VBIM sequence finding and trimming classes

The ...Seq classes are responsible for finding adapters.
The ...Match classes trim the reads.
"""
from __future__ import print_function, division, absolute_import
import sys
import re
from collections import defaultdict
from cutadapt import align, colorspace
from cutadapt.seqio import ColorspaceSequence, FastaReader, Sequence

# Constants for the find_best_alignment function.
# The function is called with SEQ1 as the adapter, SEQ2 as the read.
# TODO get rid of those constants, use strings instead
BACK = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ1
FRONT = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.START_WITHIN_SEQ1
PREFIX = align.STOP_WITHIN_SEQ2
SUFFIX = align.START_WITHIN_SEQ2
ANYWHERE = align.SEMIGLOBAL
LINKED = 'linked'


def parse_braces(sequence):
    """
    Replace all occurrences of ``x{n}`` (where x is any character) with n
    occurrences of x. Raise ValueError if the expression cannot be parsed.

    >>> parse_braces('TGA{5}CT')
    TGAAAAACT
    """
    # Simple DFA with four states, encoded in prev
    result = ''
    prev = None
    for s in re.split('(\{|\})', sequence):
        if s == '':
            continue
        if prev is None:
            if s == '{':
                raise ValueError('"{" must be used after a character')
            if s == '}':
                raise ValueError('"}" cannot be used here')
            prev = s
            result += s
        elif prev == '{':
            prev = int(s)
            if not 0 <= prev <= 10000:
                raise ValueError('Value {} invalid'.format(prev))
        elif isinstance(prev, int):
            if s != '}':
                raise ValueError('"}" expected')
            result = result[:-1] + result[-1] * prev
            prev = None
        else:
            if s != '{':
                raise ValueError('Expected "{"')
            prev = '{'
    # Check if we are in a non-terminating state
    if isinstance(prev, int) or prev == '{':
        raise ValueError("Unterminated expression")
    return result


class VBIMParser(object):
    """
    Factory for Adapter classes that all use the same parameters (error rate,
    indels etc.). The given **kwargs will be passed to the Adapter constructors.
    """

    def __init__(self, **kwargs):
        self.constructor_args = kwargs
        self.adapter_class = VBIMSeq

    def _parse_no_file(self, spec, name=None, cmdline_type='anywhere'):
        """
        Parse an adapter specification not using ``file:`` notation and return
        an object of an appropriate Adapter class. The notation for anchored
        5' and 3' adapters is supported. If the name parameter is None, then
        an attempt is made to extract the name from the specification
        (If spec is 'name=ADAPTER', name will be 'name'.)


        """
        if name is None:
            name, spec = self._extract_name(spec)
        orig_spec = spec
        types = dict(back=BACK, front=FRONT, anywhere=ANYWHERE)
        if cmdline_type not in types:
            raise ValueError('cmdline_type cannot be {0!r}'.format(cmdline_type))
        where = types[cmdline_type]

        front_anchored, back_anchored = False, False
        if spec.startswith('^'):
            spec = spec[1:]
            front_anchored = True
        if spec.endswith('$'):
            spec = spec[:-1]
            back_anchored = True

        sequence1, middle, sequence2 = spec.partition('...')
        if where == ANYWHERE:
            if front_anchored or back_anchored:
                raise ValueError("'anywhere' (-b) VBIM sequences may not be anchored")
            if middle == '...':
                raise ValueError("'anywhere' (-b) VBIM sequences may not be linked")
            return self.adapter_class(sequence=spec, where=where, name=name, **self.constructor_args)

        assert where == FRONT or where == BACK
        if middle == '...':
            raise ValueError("VBIM sequences cannot be linked")
        if front_anchored and back_anchored:
            raise ValueError('Trying to use both "^" and "$" in VBIM seq specification {!r}'.format(orig_spec))
        if front_anchored:
            if where == BACK:
                raise ValueError("Cannot anchor the 3' VBIM seq at its 5' end")
            where = PREFIX
        elif back_anchored:
            if where == FRONT:
                raise ValueError("Cannot anchor 5' VBIM seq at 3' end")
            where = SUFFIX

        return self.adapter_class(sequence=spec, where=where, name=name, **self.constructor_args)

    def parse(self, spec, cmdline_type='anywhere'):
        """
        Parse an adapter specification and yield appropriate Adapter classes.
        This works like the _parse_no_file() function above, but also supports the
        ``file:`` notation for reading adapters from an external FASTA
        file. Since a file can contain multiple adapters, this
        function is a generator.
        """

        if spec.startswith('file:'):
            # read adapter sequences from a file
            with FastaReader(spec[5:]) as fasta:
                for record in fasta:
                    name = record.name.split(None, 1)[0]
                    yield self._parse_no_file(record.sequence, name, cmdline_type)
        else:
            name, spec = self._extract_name(spec)
            yield self._parse_no_file(spec, name, cmdline_type)

    def _extract_name(self, spec):
        """
        Parse an adapter specification given as 'name=adapt' into 'name' and 'adapt'.
        """
        fields = spec.split('=', 1)
        if len(fields) > 1:
            name, spec = fields
            name = name.strip()
        else:
            name = None
        spec = spec.strip()
        return name, spec

    def parse_multi(self, back, anywhere, front):
        """
        Parse all three types of commandline options that can be used to
        specify adapters. back, anywhere and front are lists of strings,
        corresponding to the respective commandline types (-a, -b, -g).

        Return a list of appropriate Adapter classes.
        """
        adapters = []
        for specs, cmdline_type in (back, 'back'), (anywhere, 'anywhere'), (front, 'front'):
            for spec in specs:
                adapters.extend(self.parse(spec, cmdline_type))
        return adapters


class VBIMMatch(object):
    """
    Representation of a single adapter matched to a single read.

    TODO creating instances of this class is relatively slow and responsible for quite some runtime.
    """
    __slots__ = ['astart', 'astop', 'rstart', 'rstop', 'matches', 'errors',
                 'vbim', 'read', 'length', '_trimmed_read', 'adjacent_base']

    def __init__(self, astart, astop, rstart, rstop, matches, errors, vbim, read):
        self.astart = astart
        self.astop = astop
        self.rstart = rstart
        self.rstop = rstop
        self.matches = matches
        self.errors = errors
        self.vbim = vbim
        self.read = read
        self._trim()
        # Number of aligned characters in the adapter. If there are
        # indels, this may be different from the number of characters
        # in the read.
        self.length = self.astop - self.astart
        assert self.length > 0
        assert self.errors / self.length <= self.vbim.max_error_rate
        assert self.length - self.errors > 0

    def __repr__(self):
        return 'Match(astart={0}, astop={1}, rstart={2}, rstop={3}, matches={4}, errors={5})'.format(
            self.astart, self.astop, self.rstart, self.rstop, self.matches, self.errors)

    def wildcards(self, wildcard_char='N'):
        """
        Return a string that contains, for each wildcard character,
        the character that it matches. For example, if the adapter
        ATNGNA matches ATCGTA, then the string 'CT' is returned.

        If there are indels, this is not reliable as the full alignment
        is not available.
        """
        wildcards = [self.read.sequence[self.rstart + i] for i in range(self.length)
                     if self.vbim.sequence[self.astart + i] == wildcard_char and
                     self.rstart + i < len(self.read.sequence)]
        return ''.join(wildcards)

    def rest(self):
        """
        Return parts of the sequence not containing the VBIM sequence
        """
        #return self.read.sequence[:self.rstart] + self.read.sequence[self.rstop:]
        return self.read.sequence[self.rstop:]

    def get_info_record(self):
        seq = self.read.sequence
        qualities = self.read.qualities
        info = (
            self.read.name,
            self.errors,
            self.rstart,
            self.rstop,
            seq[0:self.rstart],
            seq[self.rstart:self.rstop],
            seq[self.rstop:],
            self.vbim.name
        )
        if qualities:
            info += (
                qualities[0:self.rstart],
                qualities[self.rstart:self.rstop],
                qualities[self.rstop:]
            )
        else:
            info += ('', '', '')

        return info

    def trimmed(self):
        return self._trimmed_read

    def _trim(self):
        """Compute the trimmed read"""
        """
        #Only removes VBIM tag
        self._trimmed_read = Sequence(self.read.name, self.read.sequence[:self.rstart] + self.read.sequence[self.rstop:],
                                      qualities=self.read.qualities[:self.rstart] + self.read.qualities[self.rstop:]
                                      if self.read.qualities else None, second_header=self.read.second_header, match=self)
        """
        self._trimmed_read = Sequence(self.read.name,
                                      self.read.sequence[self.rstop:],
                                      qualities=self.read.qualities[self.rstop:] if self.read.qualities else None,
                                      second_header=self.read.second_header,
                                      match=self)
        adjacent_base = self.read.sequence[self.rstart - 1]
        if adjacent_base not in 'ACGT':
            adjacent_base = ''
        self.adjacent_base = adjacent_base

    def update_statistics(self, statistics):
        """Update AdapterStatistics in place"""
        # this isn't used at all?
        statistics.errors[self.rstop][self.errors] += 1
        #statistics.adjacent_bases[self.adjacent_base] += 1


def _generate_adapter_name(_start=[1]):
    name = str(_start[0])
    _start[0] += 1
    return name


class VBIMSeq(object):
    """
    This class can find a single adapter characterized by sequence, error rate,
    type etc. within reads.

    where --  One of the BACK, FRONT, PREFIX, SUFFIX or ANYWHERE constants.
        This influences where the adapter is allowed to appear within in the
        read and also which part of the read is removed.

    sequence -- The adapter sequence as string. Will be converted to uppercase.
        Also, Us will be converted to Ts.

    max_error_rate -- Maximum allowed error rate. The error rate is
        the number of errors in the alignment divided by the length
        of the part of the alignment that matches the adapter.

    minimum_overlap -- Minimum length of the part of the alignment
        that matches the adapter.

    read_wildcards -- Whether IUPAC wildcards in the read are allowed.

    adapter_wildcards -- Whether IUPAC wildcards in the adapter are
        allowed.

    name -- optional name of the adapter. If not provided, the name is set to a
        unique number.
    """

    def __init__(self, sequence, where, max_error_rate=0.1, min_overlap=35,
                 read_wildcards=False, adapter_wildcards=True, name=None, indels=True):
        self.debug = False
        self.name = _generate_adapter_name() if name is None else name
        self.sequence = parse_braces(sequence.upper().replace('U', 'T'))  # TODO move away
        if not self.sequence:
            raise ValueError('Sequence is empty')
        self.where = where
        self.max_error_rate = max_error_rate
        self.min_overlap = min(min_overlap, len(self.sequence))
        self.indels = indels
        iupac = frozenset('XACGTURYSWKMBDHVN')
        if adapter_wildcards and not set(self.sequence) <= iupac:
            for c in self.sequence:
                if c not in iupac:
                    raise ValueError('Character {!r} in adapter sequence {!r} is '
                                     'not a valid IUPAC code. Use only characters '
                                     'XACGTURYSWKMBDHVN.'.format(c, self.sequence))
        # Optimization: Use non-wildcard matching if only ACGT is used
        self.adapter_wildcards = adapter_wildcards and not set(self.sequence) <= set('ACGT')
        self.read_wildcards = read_wildcards

        self.aligner = align.Aligner(self.sequence, self.max_error_rate,
                                     flags=self.where, wildcard_ref=self.adapter_wildcards,
                                     wildcard_query=self.read_wildcards)
        self.aligner.min_overlap = self.min_overlap
        if not self.indels:
            # TODO
            # When indels are disallowed, an entirely different algorithm
            # should be used.
            self.aligner.indel_cost = 100000



    def __repr__(self):
        return '<Adapter(name={name!r}, sequence={sequence!r}, where={where}, ' \
               'max_error_rate={max_error_rate}, min_overlap={min_overlap}, ' \
               'read_wildcards={read_wildcards}, ' \
               'adapter_wildcards={adapter_wildcards}, ' \
               'indels={indels})>'.format(**vars(self))

    def enable_debug(self):
        """
        Print out the dynamic programming matrix after matching a read to an
        adapter.
        """
        self.debug = True
        self.aligner.enable_debug()

    def match_to(self, read, match_class=VBIMMatch):
        """
        Attempt to match this adapter to the given read.

        Return a Match instance if a match was found;
        return None if no match was found given the matching criteria (minimum
        overlap length, maximum error rate).
        """
        read_seq = read.sequence.upper()  # temporary copy
        pos = -1
        # try to find an exact match first unless wildcards are allowed
        if not self.adapter_wildcards:
            if self.where == PREFIX:
                pos = 0 if read_seq.startswith(self.sequence) else -1
            elif self.where == SUFFIX:
                pos = (len(read_seq) - len(self.sequence)) if read_seq.endswith(self.sequence) else -1
            else:
                pos = read_seq.find(self.sequence)
        if pos >= 0:
            match = match_class(
                0, len(self.sequence), pos, pos + len(self.sequence),
                len(self.sequence), 0, self, read)
        else:
            # try approximate matching
            if not self.indels and self.where in (PREFIX, SUFFIX):
                if self.where == PREFIX:
                    alignment = align.compare_prefixes(self.sequence, read_seq,
                                                       wildcard_ref=self.adapter_wildcards,
                                                       wildcard_query=self.read_wildcards)
                else:
                    alignment = align.compare_suffixes(self.sequence, read_seq,
                                                       wildcard_ref=self.adapter_wildcards,
                                                       wildcard_query=self.read_wildcards)
                astart, astop, rstart, rstop, matches, errors = alignment
                if astop - astart >= self.min_overlap and errors / (astop - astart) <= self.max_error_rate:
                    match = match_class(*(alignment + (self, read)))
                else:
                    match = None
            else:
                alignment = self.aligner.locate(read_seq)
                if self.debug:
                    print(self.aligner.dpmatrix)  # pragma: no cover
                if alignment is None:
                    match = None
                else:
                    astart, astop, rstart, rstop, matches, errors = alignment
                    match = match_class(astart, astop, rstart, rstop, matches, errors, self, read)

        if match is None:
            return None
        assert match.length > 0 and match.errors / match.length <= self.max_error_rate, match
        assert match.length >= self.min_overlap
        return match

    def __len__(self):
        return len(self.sequence)

    def random_match_probabilities(self, gc_content):
        """
        Estimate probabilities that this adapter matches a
        random sequence. Indels are not taken into account.

        Returns a list p, where p[i] is the probability that
        i bases of this adapter match a random sequence with
        GC content gc_content.
        """
        seq = self.sequence
        allowed_bases = 'CGRYSKMBDHVN' if self.adapter_wildcards else 'GC'
        p = 1
        probabilities = [p]
        for i, c in enumerate(seq):
            if c in allowed_bases:
                p *= gc_content / 2.
            else:
                p *= (1 - gc_content) / 2
            probabilities.append(p)
        return probabilities

