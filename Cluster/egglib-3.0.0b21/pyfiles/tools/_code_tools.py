"""
This module contains the single function :function:`translate` which is
available from the package :mod:`tools`.
"""

__license__ = """
    Copyright 2008-2012,2015-2016 Stephane De Mita, Mathieu Siol

    This file is part of EggLib.

    EggLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EggLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EggLib.  If not, see <http://www.gnu.org/licenses/>.
"""

import re
from .. import _eggwrapper, _interface
from . import _reading_frame, _seq_manip

########################################################################

# preload all genetic codes (at the level of the class)
_codes = {}
_stops = {}
_starts = {}
for i in xrange(_eggwrapper.GeneticCode.num_codes()):
    code = _eggwrapper.GeneticCode(i)
    key = code.get_code()
    _codes[key] = code

    # record the list of start/stop codons
    _starts[key] = []
    _stops[key] = []
    for a in 'ACGT':
        for b in 'ACGT':
            for c in 'ACGT':
                codon = _eggwrapper.GeneticCode.codon(*map(ord, (a, b, c)))
                if _codes[key].aminoacid(codon) == '*': _stops[key].append(''.join([a, b, c]))
                elif _codes[key].start(codon): _starts[key].append(''.join([a, b, c]))

########################################################################

def int2codon(codon):

    """
    Return the three-character string of a codon encoded with a single
    integer (such as those obtained from :class:`CodingDiversity`).
    Return ``"???"`` if the encoded codon is out of expected range.
    """

    if codon < 0 or codon > 63: return '???'
    return ''.join((_eggwrapper.GeneticCode.base(codon, 0),
                    _eggwrapper.GeneticCode.base(codon, 1),
                    _eggwrapper.GeneticCode.base(codon, 2)))

########################################################################

class Translator(object):

    """
    Class providing methods to translate nucleotide (DNA) sequences to
    proteins.

    :param code: genetic code identifier (see :ref:`genetic-codes`).
        Required to be an integer among the valid values. The default
        value is the standard genetic code.

    :param smart: "smart" translation.

    :param delete_truncated: if ``True``, codons that are truncated
        (either because the reading frame is 5', internally, or 3'
        partial, or if the number of bases is not a multiple of three)
        are skipped when translated. By default, truncated codons are
        translated as ``X``.

    .. versionadded:: 3.0.0
    """

    ####################################################################

    def __init__(self, code=1, smart=False, delete_truncated=False):
        if code not in _codes: raise ValueError, 'unknown genetic code: {0}'.format(code)
        self._code = _codes[code]
        self._smart = smart
        self._delete_truncated = delete_truncated

    ####################################################################

    @property
    def delete_truncated(self):
        """
        Value of the *delete_truncated* option. The value can be
        modified.
        """
        return self._delete_truncated

    @delete_truncated.setter
    def delete_truncated(self, value): self._delete_truncated = value

    ####################################################################

    def translate_codon(self, first, second, third):

        """
        Translate a single codon based on the genetic code defined at
        construction time.

        :param first: first base of the codon as a one-character string.
        :param second: second base of the codon as a one-character
            string.
        :param third: third base of the codon as a one-character string.
        :return: The one-letter amino acid code if the codon can be
            translated, '-' if the codon is '---', ``X`` otherwise
            (including all cases with invalid nucleotides).
        """

        return chr(self._code.translate(ord(first), ord(second), ord(third), self._smart))

    ####################################################################

    def _help_alt(self, base1, base2, base3):

        if base1 == 45 and base2 == 45 and base3 == 45: return 0    # gap
        if base1 == 65 and base2 == 84 and base3 == 71: return 1    # normal start
        codon = self._code.codon(base1, base2, base3)
        if codon == _eggwrapper.UNKNOWN: return 1                   # missing data
        if self._code.start(codon): return 2 # alternative start
        raise ValueError, 'all sequences are required to start by a valid initiation codon'

    def translate_align(self, align, frame=None, allow_alt=False, in_place=False):

        """
        Translate a :class:`~.Align` instance.

        :param align: a :class:`~.Align` containing DNA sequences.

        :param frame: a :class:`~.ReadingFrame` instance providing the
            exon positions in the correct frame. By default, a
            non-segmented frame covering all sequences is assumed (in
            case the provided alignment is the coding region).

        :param allow_alt: a boolean telling whether alternative start
            (initiation) codons should be considered. If ``False``,
            codons are translated as a methionine (M) if, and only if,
            there are among the alternative start codons for the considered
            genetic code and they appear at the first position for the
            considered sequence (excluding all triplets of gap symbols
            appearing at the 5' end of the sequence). With this option,
            it is required that all sequences start by a valid initiation
            codon unless the first codon is partial or contains invalid
            data (in such cases, it is ignored).

        :param in_place: place translated sequences into the original
            :class:`~.Align` instance (this discards original data). By
            default, returns a new instance.

        :return: By default, an original :class:`~.Align` instance
            containing translated (protein) sequences. If *in_place* was
            ``True``, return ``None``.
        """

        # check input
        if not isinstance(align, _interface.Align): raise TypeError, 'argument must be an `Align` instance'
        ls = align.ls
        ns = align.ns
        no = align.no

        # create frame
        if frame == None:
            if ls == 0: frame = _reading_frame.ReadingFrame([])
            else: frame = _reading_frame.ReadingFrame([(0, ls)])
        elif ls < frame.num_needed_bases:
            raise ValueError, 'reading frame is extending past the end of the alignment'
        if self._delete_truncated: naa = frame.num_full_codons
        else: naa = frame.num_codons

        # create destination
        if not in_place:
            dest = _interface.Align(nsam=ns, nout=no, nsit=naa)
            ng = align.ng
            dest.ng = ng
            for i in xrange(ns):
                dest.set_name(i, align.get_name(i))
                for j in xrange(ng):
                    dest.set_label(i, j, align.get_label(i, j))
            for i in xrange(no):
                dest.set_name_o(i, align.get_name_o(i))
                dest.set_label_o(i, align.get_label_o(i))
        else:
            dest = align

        # detect alternative start codons
        if allow_alt:
            alt_list_i = []
            for idx in xrange(ns):
                i = 0
                while i < frame.num_codons:
                    a, b, c = frame.codon_bases(i)
                    x = self._help_alt(align._obj.get_i(idx, a),
                                       align._obj.get_i(idx, b),
                                       align._obj.get_i(idx, c))
                    if x == 0: i += 1
                    elif x == 1: break
                    else:
                        alt_list_i.append((idx, i))
                        break

            alt_list_o = []
            for idx in xrange(no):
                i = 0
                while i < frame.num_codons:
                    a, b, c = frame.codon_bases(i)
                    x = self._help_alt(align._obj.get_o(idx, a),
                                       align._obj.get_o(idx, b),
                                       align._obj.get_o(idx, c))
                    if x == 0: i += 1
                    elif x == 1: break
                    else:
                        alt_list_o.append((idx, i))
                        break

        # translate
        for i, (a,b,c) in enumerate(frame.iter_codons(skip_partial=self._delete_truncated)):
            for idx in xrange(ns):
                if None in (a, b, c):
                    dest._obj.set_i(idx, i, ord('X'))
                else:
                    A, B, C = align._obj.get_i(idx, a), align._obj.get_i(idx, b), align._obj.get_i(idx, c)
                    dest._obj.set_i(idx, i, self._code.translate(A, B, C, self._smart))
            for idx in xrange(no):
                if None in (a, b, c):
                    dest._obj.set_o(idx, i, ord('X'))
                else:
                    A, B, C = align._obj.get_o(idx, a), align._obj.get_o(idx, b), align._obj.get_o(idx, c)
                    dest._obj.set_o(idx, i, self._code.translate(A, B, C, self._smart))

        # fix alternative start codons
        if allow_alt:
            for i, j in alt_list_i: dest._obj.set_i(i, j, 77)
            for i, j in alt_list_o: dest._obj.set_o(i, j, 77)

        # return
        if in_place: dest._obj.set_nsit(naa)
        else: return dest

    ####################################################################

    def translate_container(self, container, allow_alt=False, in_place=False):

        """
        Translate a :class:`~.Container` instance.

        :param align: a :class:`~.Container` containing DNA sequences.

        :param allow_alt: a boolean telling whether alternative start
            (initiation) codons should be considered. If ``False``,
            codons are translated as a methionine (M) if, and only if,
            there are among the alternative start codons for the considered
            genetic code and they appear at the first position for the
            considered sequence. With this option, it is required that
            all sequences start by a valid initiation
            codon unless the first codon is missing, partial, or contains invalid
            data (in such cases, it is ignored).

        :param in_place: place translated sequences into the original
            :class:`~.Container` instance (this discards original data).
            By default, returns a new instance.

        :return: By default, an original :class:`~.Container` instance
            containing translated (protein) sequences. If *in_place* was
            ``True``, return ``None``.
        """

        # check input
        if not isinstance(container, _interface.Container): raise TypeError, 'argument must be an `Container` instance'
        ns = container.ns
        no = container.no

        # create destination
        if not in_place:
            dest = _interface.Container()
            dest._ns = ns
            dest._obj.set_nsam_i(ns)
            dest._ns_o = no
            dest._obj.set_nsam_o(no)
            ng = container.ng
            dest.ng = ng
            for i in xrange(ns):
                dest.set_name(i, container.get_name(i))
                for j in xrange(ng):
                    dest.set_label(i, j, container.get_label(i, j))
            for i in xrange(no):
                dest.set_name_o(i, container.get_name_o(i))
                dest.set_label_o(i, container.get_label_o(i))
        else:
            dest = container

        # translate
        frame = _reading_frame.ReadingFrame([])
        for idx in xrange(ns):
            ls = container.ls(idx)
            if ls == 0: frame.process([])
            else: frame.process([(0, ls)])
            naa = frame.num_full_codons if self._delete_truncated else frame.num_codons
            if not in_place: dest._obj.set_nsit_i(idx, naa)
            if allow_alt and frame.num_codons > 0:
                a, b, c = frame.codon_bases(0)
                alt_flag = None not in (a, b, c) and self._help_alt(
                    container._obj.get_i(idx, a),
                    container._obj.get_i(idx, b),
                    container._obj.get_i(idx, c)) == 2
            for i, (a,b,c) in enumerate(frame.iter_codons(skip_partial=self._delete_truncated)):
                if None in (a, b, c):
                    dest._obj.set_i(idx, i, ord('X'))
                else:
                    A, B, C = container._obj.get_i(idx, a), container._obj.get_i(idx, b), container._obj.get_i(idx, c)
                    dest._obj.set_i(idx, i, self._code.translate(A, B, C, self._smart))
            if in_place: dest._obj.set_nsit_i(idx, naa)
            if allow_alt and alt_flag: dest._obj.set_i(idx, 0, 77)
        for idx in xrange(no):
            ls = container.lo(idx)
            if ls == 0: frame.process([])
            else: frame.process([(0, ls)])
            naa = frame.num_full_codons if self._delete_truncated else frame.num_codons
            if not in_place: dest._obj.set_nsit_o(idx, naa)
            if allow_alt and frame.num_codons > 0:
                a, b, c = frame.codon_bases(0)
                alt_flag = None not in (a, b, c) and self._help_alt(
                    container._obj.get_o(idx, a),
                    container._obj.get_o(idx, b),
                    container._obj.get_o(idx, c)) == 2
            for i, (a,b,c) in enumerate(frame.iter_codons(skip_partial=self._delete_truncated)):
                if None in (a, b, c):
                    dest._obj.set_o(idx, i, ord('X'))
                else:
                    A, B, C = container._obj.get_o(idx, a), container._obj.get_o(idx, b), container._obj.get_o(idx, c)
                    dest._obj.set_o(idx, i, self._code.translate(A, B, C, self._smart))
            if in_place: dest._obj.set_nsit_o(idx, naa)
            if allow_alt and alt_flag: dest._obj.set_o(idx, 0, 77)

        # return
        if not in_place: return dest

    ####################################################################

    def translate_sequence(self, sequence, frame=None, allow_alt=False):

        """
        Translate a sequence.

        :param sequence: a :class:`str`, :class:`~.SequenceView` or
            compatible instance containing DNA sequences. All sequence
            types yielding one-character strings are acceptable.

        :param frame: a :class:`~.ReadingFrame` instance providing the
            exon positions in the correct frame. By default, a
            non-segmented frame covering all sequences is assumed (in
            case the provided alignment is the coding region).

        :param allow_alt: a boolean telling whether alternative start
            (initiation) codons should be considered. If ``False``,
            codons are translated as a methionine (M) if, and only if,
            there are among the alternative start codons for the considered
            genetic code and they appear at the first position for the
            sequence. With this option, it is required that the sequence starts by a valid initiation
            codon unless the first codon is missing, partial, or contains invalid
            data (in such cases, it is ignored).

        :return: A new :class:`str` instance containing translated
            (protein) sequences.

        .. note::
           Character mapping is ignored for this method (sequences must
           be provided as is, with ``A``, ``C``, ``G`` and ``T`` and
           ambiguity characters if relevant; invalid characters cause
           translation to ``X`` in resulting amino acid sequences).
        """

        if not isinstance(sequence, _interface.SequenceView): sequence = map(ord, sequence)
        ls = len(sequence)
        if frame is None:
            if ls == 0: frame = _reading_frame.ReadingFrame([])
            else: frame = _reading_frame.ReadingFrame([(0, ls)])
        elif frame.num_needed_bases > ls:
            raise ValueError, 'reading frame is extending past the end of the sequence'
        naa = frame.num_full_codons if self._delete_truncated else frame.num_codons
        dest = [None] * naa
        for i, (a,b,c) in enumerate(frame.iter_codons(skip_partial=self._delete_truncated)):
            if None in (a, b, c):
                dest[i] = 'X'
            else:
                dest[i] = chr(self._code.translate(sequence[a], sequence[b], sequence[c], self._smart))
        if allow_alt and naa > 0:
            a, b, c = frame.codon_bases(0)
            if None not in (a, b, c) and self._help_alt( sequence[a], sequence[b], sequence[c]) == 2: dest[0] = 'M'
        return ''.join(dest)

########################################################################

def translate(seq, frame=None, code=1, smart=False, delete_truncated=False, allow_alt=False, in_place=False):

    """
    Translates DNA nucleotide sequences to proteins. See the class
    :class:`~.Translator` for more details. This is a convenience method
    allowing to translate nucleotide sequences in a single call. For
    repeatitive calls, direct use of :class:`~.Translator` can be more
    efficient.

    :param seq: input DNA nucleotide sequences. Accepted types are
        :class:`~.Align`, :class:`~.Container`, :class:`~.SequenceView`
        or :class:`str` (or also types compatible with :class:`str`).

    :param frame: reading frame as a :class:`~.ReadingFrame` instance
        (see :meth:`~.Translator.translate_align` for details). Not
        allowed if *seq* is a :class:`~.Container` instance.

    :param code: genetic code identifier (see :ref:`genetic-codes`).
        Required to be an integer among the valid values. The default
        value is the standard genetic code.

    :param smart: "smart" translation.

    :param delete_truncated: skip codons that are truncated (by default,
        they are retained and translated as ``X``. See
        :class:`~.Translator` for details.

    :param allow_alt: a boolean telling whether alternative start
        (initiation) codons should be considered. If ``False``,
        codons are translated as a methionine (M) if, and only if,
        there are among the alternative start codons for the considered
        genetic code and they appear at the first position for the
        sequence. With this option, it is required that the sequence starts
        by a valid initiation codon unless the first codon is missing,
        partial, or contains invalid data (in such cases, it is ignored).
        If *seq* is an :class:`.Align`, leading gaps are ignored as long as
        they are a multiple of 3 (fully missing triplets).

    :param in_place: place translated sequences in the provided
        :class:`~.Align` or :class:`~.Container` instance, overwritting
        initial data. Not allowed if *seq* is not of one of these two
        types. See :meth:`~.Translator.translate_align` for details).

    :return: Protein sequences as an :class:`~.Align` or a
        :class:`~.Container` if either of these types have been provided
        as *seq*, or a :class:`str` otherwise. If *in_place* has been
        provided, returns ``None``.

    .. versionchanged:: 3.0.0
       Added options *code*, *frame* and the options to define
       nucleotide values (default values are backward-compatible).
       Option *strip* is removed. Functionality is moved to class
       :class:`~.Translator`.
    """

    # check input type

    if isinstance(seq, _interface.Align):
        T = Translator(code=code, smart=smart, delete_truncated=delete_truncated)
        return T.translate_align(seq, frame=frame, allow_alt=allow_alt, in_place=in_place)
    elif isinstance(seq, _interface.Container):
        if frame is not None: raise ValueError, 'cannot specify a reading frame for translating a `Container`'
        T = Translator(code=code, smart=smart, delete_truncated=delete_truncated)
        return T.translate_container(seq, allow_alt=allow_alt, in_place=in_place)
    else:
        if in_place == True: raise ValueError, 'cannot translate sequence in place if it is not provided as a `Container` or an `Align`'
        T = Translator(smart=smart, delete_truncated=delete_truncated)
        return T.translate_sequence(seq, allow_alt=allow_alt, frame=frame)

########################################################################

def orf_iter(sequence, code=1, smart=False, min_length=1,
             forward_only=False, force_start=True, start_ATG=False,
             force_stop=True):

    """
    Return an iterator over non-segmented open reading frames (ORFs)
    found in a provided DNA sequence in any of the six possible frames.

    :param sequence: a :class:`str` representing a DNA sequence.

    :param code: genetic code identifier (see :ref:`genetic-codes`).
        Required to be an integer among the valid values. The default
        value is the standard genetic code.

    :param smart: "smart" translation.

    :param min_length: minimum length of returned ORFs. This value must
        be at least 1. It is understood as the length of the encoded
        peptide (therefore, selected ORFs will have a length of least
        three times this value).

    :param forward_only: consider only the three forward frames (do not
        consider the reverse strand).

    :param force_start: if ``True``, all returned ORFs are required to
        start with a start codon. Otherwise, all non-stop codons are
        included in ORFs. Note that alternative start codons (``CTG``
        and ``TTG`` for the standard genetic code) are also supported.

    :param start_ATG: consider ``ATG`` as the only start codon (ignored
        if *force_start* is ``False``.

    :param force_stop: if ``True``, all returned ORFs are required to
        end with a stop codon (this only excludes 3'-partial ORFs, that
        is ORFs that end with the end of the provided sequennce).

    :return: An iterator over all detected ORFs. Each ORF is represented
        by a ``(start, stop, length, frame)`` tuple where *start* is the
        start position of the ORF and *stop* the stop position (such as
        ``sequence[start:stop]`` returns the ORF sequence or its reverse
        complement ), *length* is the ORF length and *frame* is the
        reading frame on which it was found: +1, +2, +3 are the frames
        on the forward strand (starting respectively at the first,
        second, and third base), and -1, -2, -3 are the frames on the
        reverse strand (starting respectively at the last, last but one,
        and last but two base).

    .. versionadded:: 3.0.0
        Take over the basic functionality of :func:`.longest_orf`;
        turned into an iterator. This new version comes with a better
        implementation, a changed signature (with a new *forward_only*
        option) and returns ORF positions instead of sequences.
    """

    # get/check arguments
    if code not in _codes: raise ValueError, 'unknown genetic code: {0}'.format(code)
    if min_length < 1: raise ValueError, 'value for `min_length` is too small'
    min_length *= 3
    ls = len(sequence)
    if force_start:
        if start_ATG: list_starts = ['ATG']
        else: list_starts = _starts[code]

    # initialize dict with position of stop codons in six (or three) frames
    frames = { +1: [], +2: [], +3: []}
    if not forward_only: frames.update({ -1: [], -2: [], -3: []})

    # get the positions of all stop codons in all six (or three) frames
    for stop in _stops[code]:

        # get positions of stops in the three forward frames
        for mo in re.finditer('(?={0})'.format(stop), sequence):
            pos = mo.start()
            frames[pos%3+1].append(pos)

        if not forward_only:
            for mo in re.finditer('(?={0})'.format(_seq_manip.rc(stop)), sequence):
                pos = mo.start()
                rpos = (ls - 1) - (pos + 2)
                frames[-(rpos%3+1)].append(pos)

    # add the segment between the last stop and the end of the sequence (only forward)
    if not force_stop:
        for frame_id in [+1, +2, +3]:
            # add the position of the last full codon (in the proper frame)
            frames[frame_id].append( (ls-3) - (ls%3 - (frame_id-1)) % 3)

    # process all ORFs
    for idx, (frame_id, frame) in enumerate(frames.iteritems()):
        frame.sort()

        # trick to skip the segment between the start of sequence and the last stop (only reverse)
        if frame_id < 0 and force_stop == True and len(frame) > 0: cur = frame[0]
        else: cur = 0

        for pos in frame:

            # pos is the position of the stop codon
            if frame_id > 0: stop = pos + 3 # stop codon included
            else: stop = pos # stop codon excluded (it is stop of next ORF)

            # adjust length to set the start on the right frame (remove the remainder)
            length = (stop - cur) - (stop-cur)%3
            start = stop - length

            # start position for next orf is immediately past this stop codon
            # need to be set before "seek start" block because this can shift stop
            cur = stop

            # seek start codon if requested
            if force_start == True:
                if frame_id > 0:
                    while sequence[start:start+3] not in list_starts and length > 0:
                        start += 3
                        length -= 3
                else:
                    while _seq_manip.rc(sequence[stop-3:stop]) not in list_starts and length > 0:
                        stop -= 3
                        length -= 3

            # yield the orf if long enough
            if length >= min_length:
                yield start, stop, length, frame_id

########################################################################

def longest_orf(*args, **kwargs):

    """
    Detect the longest open reading frame 

    Arguments are identical to :func:`.orf_iter`.

    :return: A ``(start, stop, length, frame)`` tuple (see the return
        value of :func:`.orf_iter` for details), or ``None`` if no open
        reading frame fits the requirements (typically, the minimum
        length).

    An :exc:`exceptions.ValueError` is raised if two or more open
    reading frames have the largest length.

    .. versionchanged:: 2.0.1
        Added options; return the trailing stop codon when appropriate.

    .. versionchanged:: 2.1.0
        Added option *mini*. The behaviour of previous versions is
        reproduced by setting *mini* to 0.

    .. versionchanged:: 3.0.0
        Most of the functionality is moved to :func:`.orf_iter` with an
        updated interface.
    """

    max_length = 0
    candidate = None

    for start, stop, length, frame in orf_iter(*args, **kwargs):
        if length == max_length:
            candidate = None
        elif length > max_length:
            max_length = length
            candidate = start, stop, length, frame

    if max_length == 0: return None
    if candidate == None: raise ValueError, 'several equally long ORF found'
    return candidate

########################################################################

def backalign(nucl, aln, code=1, smart=False, ignore_names=False, ignore_mismatches=False, fix_stop=False):
    
    """
    Align coding sequence based on the corresponding protein alignment.

    :param nucl: a :class:`.Container` or :class:`.Align` instance
        containing coding DNA sequences that should be aligned. Codons
        containing an ambiguity IPUAC character or ``?`` (see
        :ref:`iupac-nomenclature`) are translated as ``X``. All
        alignment gaps (character ``-``) will be stripped from the
        sequences.

    :param aln: a :class::`.Align` instance containing an alignment of
        the protein sequences encoded by the coding sequences provided
        as value for the *nucl* argument. If there is an outgroup in
        *nucl*, it must also be present in *aln*. Group labels of *aln*
        are not taken into account.

    :param code: genetic code identifier (see :ref:`genetic-codes`).
        Required to be an integer among the valid values. The default
        value is the standard genetic code.

    :param smart: "smart" translation.

    :param ignore_names: if ``True``, ignore names for matching
        sequences in the protein alignment to coding sequences.
        Sequences will be matched using their rank and the names in the
        returned alignment will be taken from *nucl*.

    :param fix_stop: if ``True``, support a single trailing stop codon
        in coding sequences not represented by a ``*`` in the provided
        protein alignment (such as if the final stop codons have been
        stripped during alignment). If found, this stop codon will be
        flushed as left as possible (immediately after the last non-gap
        character) in the returned coding alignment.

    :param ignore_mismatches: if ``True``, do not generate any exception
        if a predicted protein does not match the provided protein
        sequence (if the lengths differ, an exception is always raised).

    :return: A :class:`.Align` instance containing aligned coding DNA
        sequence (including the outgroup)

    If a mismatch is detected between a protein from *aln* and the
    corresponding prediction from *nucl*, an instance of
    :exc:`.BackalignError` (a subclass of :exc:`.exceptions.ValueError`)
    is raised. The attribute :attr:`.BackalignError.alignment` can be
    used to help identify the reason of the error. Mismatches (but not
    differences of length) can be ignored with the option
    *ignore_mismatches*.
    """

    # get code
    translator = Translator(code=code, smart=smart, delete_truncated=False)

    # initialize the return object
    ns = aln.ns
    no = aln.no
    ls = aln.ls
    ng = nucl.ng
    result = _eggwrapper.DataHolder(True)
    result.set_nsam_i(ns)
    result.set_nsam_o(no)
    result.set_nsit(ls * 3 + (3 if fix_stop else 0))
    result.set_ngroups(ng)

    # set 3 trailing stops at the end of all sequences to allow 1 or more additional stop codons
    if fix_stop:
        for i in xrange(ns):
            for j in xrange(ls * 3, ls * 3 + 3): result.set_i(i, j, 45)
        for i in xrange(no):
            for j in xrange(ls * 3, ls * 3 + 3): result.set_o(i, j, 45)
        last_used = False # set to True if at least one sequence uses the last three bases

    # get mapping of indexes
    nucl = nucl._obj
    aln = aln._obj
    if nucl.get_nsam_i() != ns: raise ValueError, 'numbers of sequences don\'t match'
    if nucl.get_nsam_o() != no: raise ValueError, 'numbers of outgroup sequences don\'t match'
    names_i = [nucl.get_name_i(i) for i in xrange(ns)]
    names_o = [nucl.get_name_o(i) for i in xrange(no)]
    if not ignore_names:
        if len(set(names_i)) < ns: raise ValueError, 'duplicate name found'
        if len(set(names_o)) < no: raise ValueError, 'duplicate name found in outgroup'
        aln_names_i = [aln.get_name_i(i) for i in xrange(ns)]
        aln_names_o = [aln.get_name_o(i) for i in xrange(no)]
        rank_i = []
        rank_o = []
        for name in names_i:
            try: rank_i.append(aln_names_i.index(name))
            except ValueError: 'sequence name not found in protein alignment: {0}'.format(name)
        for name in names_o:
            try: rank_o.append(aln_names_o.index(name))
            except ValueError: 'outgroup sequence name not found in protein alignment: {0}'.format(name)
    else:
        rank_i = range(ns)
        rank_o = range(no)

    # process all ingroup samples
    for i, j in enumerate(rank_i):
        result.set_name_i(i, names_i[i])
        for k in xrange(ng):
            result.set_group_i(i, k, nucl.get_group_i(i, k))
        exitcode = _backalign_helper(nucl.get_i, aln.get_i, result.set_i, result.get_i, locals())
        if exitcode == 0: raise BackalignError(names_i[i], nucl.get_i, aln.get_i, i, j, ls, nucl.get_nsit_i(i), translator)
        if fix_stop: last_used |= (exitcode == 2)

    # same for outgroup samples
    for i, j in enumerate(rank_o):
        result.set_names_o(i, name_o[i])
        result.set_group_o(i, nucl.get_group_o(i))
        exitcode = _backalign_helper(nucl.get_o, aln.get_o, result.set_o, result.get_o, locals())
        if exitcode == 0: raise BackalignError(names_o[i], nucl.get_o, aln.get_o, i, j, ls, nucl.get_nsit_i(i), translator)
        if fix_stop: last_used |= (exitcode == 2)

    # if the last codon position has not been used (by fix_stop), remove them
    if fix_stop and not last_used:
        result.set_nsit(ls * 3) # decrease by three

    # that's all, folks!
    return _interface.Align._create_from_data_holder(result)

########                                                        ########

def _backalign_helper(fnucl, faln, fset, fget, variables):

    # import main variables from backalign()'s namespace
    i_nucl = variables['i']
    i_aln = variables['j']
    ls_nucl = variables['nucl'].get_nsit_i(i_nucl)
    translator = variables['translator']

    # safety
    if ls_nucl < 1: raise ValueError, 'empty nucleotide sequence'

    # process all aligned aa positions
    c = 0
    for idx in xrange(variables['ls']):
        aa = faln(i_aln, idx)

        # process gaps
        if aa == 45:
            fset(i_nucl, 3*idx, 45)
            fset(i_nucl, 3*idx+1, 45)
            fset(i_nucl, 3*idx+2, 45)

        # non-gap
        else:
            codon = []

            # get codon
            for i in xrange(3):
                while True:
                    if c == ls_nucl: return 0
                    ch = fnucl(i_nucl, c)
                    c += 1
                    if ch != 45: break
                codon.append(ch)

            # check mismatch
            if (not variables['ignore_mismatches'] and
                translator.translate_codon(* map(chr, codon)) != chr(aa)):
                    return 0

            # set codon
            fset(i_nucl, 3*idx, codon[0])
            fset(i_nucl, 3*idx+1, codon[1])
            fset(i_nucl, 3*idx+2, codon[2])

    # if there is additional nucleotides at the end
    if c != ls_nucl:

        # if requested, check that it is a stop codon
        if variables['fix_stop']:
            codon = [None, None, None]

            # get next three bases
            for i in xrange(3):
                while True:
                    if c == ls_nucl: return 0 # if <3 bases: error
                    ch = fnucl(i_nucl, c)
                    c += 1
                    if ch != 45: break
                codon[i] = ch

            # fix stop codon (only if all bases have been read)
            if c == ls_nucl and translator.translate_codon(* map(chr, codon)) == '*':
                i = variables['ls'] * 3 # there is necessarily room there, filled by gaps
                while i > 0 and fget(i_nucl, i-3) == 45:
                    i-=3
                fset(i_nucl, i, codon[0])
                fset(i_nucl, i+1, codon[1])
                fset(i_nucl, i+2, codon[2])
                if i == variables['ls'] * 3: return 2 # used last codon position
                else: return 1

        # if no caught by stop-fixing block: error
        return 0

    # no errors, no stop codon fixing
    return 1

########                                                        ########

class BackalignError(ValueError):

    """
    Subclass of :class:`~exceptions.ValueError` used to report errors
    occurring during the use of :func:`.backalign` because of mismatches
    between the provided alignment and predicted proteins.
    """

    ####################################################################

    def __init__(self, name, fnuc, faln, i_nuc, i_aln, ls_aln, ls_nuc, translator):

        self._name = name
        message = 'mismatch between provided and predicted proteins for {0}'.format(name)

        prov = ''.join(map(chr, [faln(i_aln, i) for i in xrange(ls_aln)])).translate(None, '-')
        pred = ''.join(map(chr, [fnuc(i_nuc, i) for i in xrange(ls_nuc)])).translate(None, '-')
        if len(pred) % 3 != 0: message = 'length of nucleotide sequence not a multiple of 3 for {0}'.format(name)
        pred = ''.join([translator.translate_codon(* codon) for codon in zip(pred[::3], pred[1::3], pred[2::3])])

        n = max([len(prov), len(pred)])

        mid = []
        for i in xrange(n):
            if i >= len(prov):
                prov += '-'
                mid.append('~')
            elif i >= len(pred):
                pred += '-'
                mid.append('~')
            elif prov[i] != pred[i]:
                mid.append('#')
            else:
                mid.append('|')
        mid = ''.join(mid)

        c = 0
        self._alignment = []
        while True:
            A = prov[c:c+60]
            A = ' '.join([A[i*10:i*10+10] for i in range(6)])
            B = mid[c:c+60]
            B = ' '.join([B[i*10:i*10+10] for i in range(6)])
            C = pred[c:c+60]
            C = ' '.join([C[i*10:i*10+10] for i in range(6)])
            self._alignment.append('[provided]  {0}'.format(A) + '[{0}]\n'.format(c+60).rjust(7))
            self._alignment.append('            {0}\n'.format(B))
            self._alignment.append('[predicted] {0}'.format(C) + '[{0}]\n'.format(c+60).rjust(7))
            c += 40
            if c < n: self._alignment.append('\n')
            else: break

        self._alignment = ''.join(self._alignment)
        super(ValueError, self).__init__(message)

    ####################################################################

    @property
    def name(self):

        """
        Name of sequence for which the error occurred.
        """

        return self._name

    ####################################################################

    @property
    def alignment(self):

        """
        String representing the alignment of the provided and predicted
        proteins, incorporating the following codes in the middle line:

        * ``|``: match.
        * ``#``: mismatch.
        * ``~``: one protein shorter.
          differ).
        """

        return self._alignment

########################################################################

def trailing_stops(align, frame=None, action=0, code=1, smart=False,
            include_outgroup=False, gap=45, replacement=(63, 63, 63)):

    """
    trailing_stops(align, frame=None, action=0, code=None, include_outgroup=False, gap='-', replacement='???')

    Detect and (optionally fix) stop codons at the end of the sequences.
    The last three non-gap data of a sequence must form a single codon
    in the specified frame, meaning that if the final codon is
    interrupted by a gap of shifted out of frame, it will not be
    detected as a stop codon. If the last codon is truncated, it will
    also be considered as non gap. If the last codon is fully falling in
    a gap, the previous one is considered (and so on).

    :param align: a :class:`.Align` containing aligned coding DNA
        sequences.

    :param frame: a :class:`~.ReadingFrame` instance providing the exon
        positions in the correct frame. By default, a non-segmented
        frame covering all sequences is assumed (in case the provided
        alignment is the coding region).

    :param action: an integer specifying what should be done if a stop
        codon is found at the end of a given sequence. Possible actions
        are listed in the following table:

        +------+---------------------------------------------------+
        | Code | Action                                            +
        +======+===================================================+
        | 0    | Nothing (just count them).                        +
        +------+---------------------------------------------------+
        | 1    | Replace them by gaps, and delete the final three  +
        |      | positions if they are made by gaps only.          +
        +------+---------------------------------------------------+
        | 2    | Replace them by the value given as *replacement*. +
        +------+---------------------------------------------------+

        Note that using ``action=1`` is not stricly equivalent to using
        ``action=2, replacement='---'`` because the former deletes the
        last three positions of the alignment if needed while the latter
        does not.

    :param code: genetic code identifier (see :ref:`genetic-codes`).
        Required to be an integer among the valid values. The default
        value is the standard genetic code.

    :param smart: "smart" translation.

    :param include_outgroup: if ``True``, process both ingroup and
        outgroup samples; if ``False``, process only ingroup.

    :param gap: the character representing gaps. It is allowed to pass
        the allele value as a single-character :class:`str`, a
        single-character :class:`unicode`, and as a :class:`int`.

    :param replacement: if *action* is set to 2, provide the three
        values that should be used to replace stop codons. This value
        must be a three-character :class:`str` or a three-item sequence
        of integers. By default, replace final stop codons by
        uncharacterized bases.

    :return: The number of sequences that had a trailing stop codons
        among the considered sequences (including outgroup if
        *include_outgroup* is ``True``).
    """

    # argument check
    if not isinstance(align, _interface.Align): raise TypeError, '`align` must be an `Align` instance'
    if code not in _codes: raise ValueError, 'unknown genetic code: {0}'.format(code)
    code = _codes[code]
    ns = align.ns
    no = align.no
    ls = align.ls
    if isinstance(gap, basestring):
        if len(gap) != 1: raise ValueError, '`gap` must be a single-character string or an int'
        gap = ord(gap)
    if len(replacement) != 3: raise ValueError, '`replacement` must contain three items'
    if isinstance(replacement, basestring):
        replacement = map(ord, replacement)

    # create frame
    if frame == None:
        if ls == 0: frame = _reading_frame.ReadingFrame([])
        else: frame = _reading_frame.ReadingFrame([(0, ls)])
    elif ls < frame.num_needed_bases:
        raise ValueError, 'reading frame is extending past the end of the alignment'
    if frame.num_codons == 0: return 0

    cnt = 0
    num_gapped = 0

    # process ingroup
    for i in xrange(ns):
        cur = frame.num_codons - 1
        while True:
            a, b, c = frame.codon_bases(cur)
            if a == None or c == None:
                aa = 0
                break
            aa = code.translate(align._obj.get_i(i, a), align._obj.get_i(i, b), align._obj.get_i(i, c), smart)
            if aa != gap: break
            cur -= 1
            if cur < 0: break
        if aa == 42:
            cnt += 1
            if action == 1:
                align._obj.set_i(i, a, 45)
                align._obj.set_i(i, b, 45)
                align._obj.set_i(i, c, 45)
                num_gapped += 1
            elif action == 2:
                align._obj.set_i(i, a, replacement[0])
                align._obj.set_i(i, b, replacement[1])
                align._obj.set_i(i, c, replacement[2])

    # process outgroup
    if include_outgroup:
        for i in xrange(no):    
            cur = frame.num_codons - 1
            while True:
                a, b, c = frame.codon_bases(cur)
                if a == None or c == None:
                    aa = 0
                    break
                aa = code.translate(align._obj.get_o(i, a), align._obj.get_o(i, b), align._obj.get_o(i, c), smart)
                if aa != gap: break
                cur -= 1
                if cur < 0: break
            if aa == 42:
                cnt += 1
                if action == 1:
                    align._obj.set_o(i, a, 45)
                    align._obj.set_o(i, b, 45)
                    align._obj.set_o(i, c, 45)
                    num_gapped += 1
                elif action == 2:
                    align._obj.set_o(i, a, replacement[0])
                    align._obj.set_o(i, b, replacement[1])
                    align._obj.set_o(i, c, replacement[2])
            elif action == 1:
                flag |= False

    # remove final 3 gaps (only if ALL samples where replaced by gaps)
    if num_gapped == ns + no: # implies action==1
        align._obj.set_nsit(ls-3)

    return cnt

########################################################################

def iter_stops(align, frame=None, code=1, smart=False, include_outgroup=False):

    """
    Return an iterator providing the coordinates of all stop codons in
    the alignment (over all sequences). Only stop codons in the
    specified frame are detected, excluding all those that are segmented
    by a gap or are in a shifted frame. Each iteration returns a
    ``(sample, position, flag)`` where *sample* is the sample index,
    *position* is the position of the first base of the stop codon, and
    *flag* is ``True`` if the sample belongs to the ingroup and
    ``False`` if the sample belongs to the outgroup.

    :param align: a :class:`.Align` containing aligned coding DNA
        sequences. Values should consist in A, C, G and T (and possibly
        IUPAC ambiguity codes if *smart* is toggled). Other values are
        treated as missing data.

    :param frame: a :class:`~.ReadingFrame` instance providing the exon
        positions in the correct frame. By default, a non-segmented
        frame covering all sequences is assumed (in case the provided
        alignment is the coding region).

    :param code: genetic code identifier (see :ref:`genetic-codes`).
        Required to be an integer among the valid values. The default
        value is the standard genetic code.

    :param smart: "smart" translation.

    :param include_outgroup: if ``True``, process both ingroup and
        outgroup samples; if ``False``, process only ingroup.

    :return: An iterator over the `(sample, position, flag)`` tuples
        corresponding to each stop codon found in the alignment (see
        above).
    """

    # argument check
    if not isinstance(align, _interface.Align): raise TypeError, '`align` must be an `Align` instance'
    ns = align.ns
    no = align.no
    if code not in _codes: raise ValueError, 'unknown genetic code: {0}'.format(code)
    code = _codes[code]

    # create frame
    if frame == None:
        if align.ls == 0: frame = _reading_frame.ReadingFrame([])
        else: frame = _reading_frame.ReadingFrame([(0, align.ls)])
    elif align.ls < frame.num_needed_bases:
        raise ValueError, 'reading frame is extending past the end of the alignment'

    # process ingroup
    for i in xrange(ns):
        for a, b, c in frame.iter_codons(skip_partial=True):
            if code.translate(align._obj.get_i(i, a), align._obj.get_i(i, b), align._obj.get_i(i, c), smart) == 42:
                yield (i, a, True)

    # process outgroup
    if include_outgroup:
        for i in xrange(no):    
            for a, b, c in frame.iter_codons(skip_partial=True):
                if code.translate(align._obj.get_o(i, a), align._obj.get_o(i, b), align._obj.get_o(i, c), smart) == 42:
                    yield (i, a, False)

########################################################################

def has_stop(align, frame=None, code=1, smart=False, include_outgroup=False):

    """
    Return ``True`` if the alignment contains at least one codon stop
    at any position in any sequence, and ``False`` otherwise. Only sto
    codons in the specified frame are detected, excluding all those tha
    are segmented by a gap or are in a shifted frame.

    :param align: a :class:`.Align` containing aligned coding DNA
        sequences.

    :param frame: a :class:`~.ReadingFrame` instance providing the exon
        positions in the correct frame. By default, a non-segmented
        frame covering all sequences is assumed (in case the provided
        alignment is the coding region).

    :param frame: a :class:`~.ReadingFrame` instance providing the exon
        positions in the correct frame. By default, a non-segmented
        frame covering all sequences is assumed (in case the provided
        alignment is the coding region).

    :param code: genetic code identifier (see :ref:`genetic-codes`).
        Required to be an integer among the valid values. The default
        value is the standard genetic code.

    :param include_outgroup: if ``True``, process both ingroup and
        outgroup samples; if ``False``, process only ingroup.

    :return: A boolean.
    """

    for stop in iter_stops(align, frame=frame, code=code, smart=smart, include_outgroup=include_outgroup): break
    else: return False
    return True
