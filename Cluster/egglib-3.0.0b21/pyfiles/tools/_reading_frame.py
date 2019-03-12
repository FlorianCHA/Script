"""
This module contains the single class :class:`ReadingFrame` which is
available from the package :mod:`tools`.
"""

__license__ = """
    Copyright 2008-2012, 2015 Stephane De Mita, Mathieu Siol

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

########################################################################

class ReadingFrame(object):
    
    """
    Handles reading frame positions. The reading frame positions can be
    loaded as constructor argument or using the method
    :meth:`~.ReadingFrame.process`. By default, builds an instance with
    no exons. If the argument is specified, it must be identical to the
    argument to the method :meth:`~.ReadingFrame.process`.

    .. versionchanged:: 3.0.0
       Previously, bases from truncated codons were discarded; they are
       not included as part of partial codons. Functionality extended.
    """

    def __init__(self, frame = None):

        if frame == None: frame = []
        self.process(frame)

    ####################################################################

    def process(self, frame):

        """
        Load a reading frame. All previously loaded data are discarded.

        :param frame: the reading frame specification must be a sequence
            of ``(start, stop[, codon_start])`` pairs or triplets where
            *start* and *stop* give the limits of an exon, such that
            ``sequence[start:stop]`` returns the exon sequence, and
            *codon_start*, if specified, can be:
            
            * 1 if the first position of the exon is the first position
              of a codon (e.g. ``ATG ATG``),

            * 2 if the first position of the segment is the second
              position of a codon (e.g. ``TG ATG``),

            * 3 if the first position of the segment is the third
              position a of codon (e.g. ``G ATG``),

            * ``None`` if the reading frame is continuing the previous
              exon.

        If *codon_start* of the first segment is ``None``, 1 will be
        assumed. If *codon_start* of any non-first segment is not
        ``None``, the reading frame is supposed to be interupted. This
        means that if any codon was not completed at the end of the
        previous exon, it will remain incomplete.
        """

        cur = -1
        self._starts = []
        self._stops = []
        self._codon_starts = []
        for exon in frame:
            if len(exon) == 2:
                start, stop = exon
                codon_start = None
            elif len(exon) == 3:
                start, stop, codon_start  = exon
                if codon_start not in set([1, 2, 3, None]): raise ValueError, 'invalid value for `codon_start`: {0}'.format(codon_start)
            if start < cur: raise ValueError, 'exon limits must in increasing order'
            if stop <= start: raise ValueError, 'exon limits must in increasing order'
            cur = stop
            self._starts.append(start)
            self._stops.append(stop)
            self._codon_starts.append(codon_start)

        self._num_exons = len(self._starts)
        if self._num_exons > 0:
            self._num_tot_bases = self._stops[-1] - self._starts[0]
            self._num_exon_bases = sum([stop-start for (start, stop) in zip(self._starts, self._stops)])
            if self._codon_starts[0] == None: self._codon_starts[0] = 1 # important
        else:
            self._num_tot_bases = 0
            self._num_exon_bases = 0

        self._bases = {} # keys are bases, values are (exon, codon) indexes pairs
        self._codons = []

        for idx, (start, stop, codon_start) in enumerate(zip(self._starts, self._stops, self._codon_starts)):

            # if codon_start is specified, complete the previous codon
            if codon_start != None:
                if len(self._codons):
                    self._codons[-1] += [None] * (3 - len(self._codons[-1]))

                # add non codon with shift if needed
                if codon_start == 1: self._codons.append([])
                elif codon_start == 2: self._codons.append([None])
                elif codon_start == 3: self._codons.append([None, None])
                else: raise ValueError, 'invalid value for codon start'

            # loop over bases, adding codons when the current is filled
            # it is guaranteed that there will always be one codon (because 1st exon is in frame 1)
            for i in xrange(start, stop):
                if len(self._codons[-1]) == 3: self._codons.append([])
                self._codons[-1].append(i)
                self._bases[i] = (idx, len(self._codons) - 1)

        # final processing & complete last codon if needed
        self._num_codons = len(self._codons)
        if self._num_codons > 0: self._codons[-1] += [None] * (3 - len(self._codons[-1]))
        self._num_full_codons = sum([1 for i in self._codons if None not in i])
        self._codons = map(tuple, self._codons)
        if self._num_exons > 0: self._needed_bases = self._stops[-1] - 1
        else: self._needed_bases = 0

    ####################################################################

    @property
    def num_needed_bases(self):

        """
        Number of bases needed for a sequence to apply this reading
        frame. In practice, the value equals to the end of the last exon
        plus one, or zero if there is no exon. If the reading frame is
        used with a shorted sequence, it can lead to errors.
        """

        return self._needed_bases

    ####################################################################

    @property
    def num_tot_bases(self):

        """
        Total number of bases (starting from the start of the first exon
        up to end of the last one).
        """

        return self._num_tot_bases

    ####################################################################

    @property
    def num_exon_bases(self):

        """
        Number of bases in exons.
        """

        return self._num_exon_bases

    ####################################################################

    @property
    def num_exons(self):

        """
        Number of exons.
        """

        return self._num_exons

    ####################################################################

    @property
    def num_codons(self):

        """
        Number of codons (including truncated codons).
        """

        return self._num_codons

    ####################################################################

    @property
    def num_full_codons(self):

        """
        Number of full (non-truncated) codons.
        """

        return self._num_full_codons

    ####################################################################

    def exon_index(self, base):
        
        """
        Find the exon in which a given base falls.

        :param base: any base index.
        :return: The index of the corresponding exon, or ``None`` if the
            base does not fall in any exon.
        """

        if base in self._bases: return self._bases[base][0]
        else: return None

    ####################################################################

    def codon_index(self, base):
        
        """
        Find the codon in which a given base falls.

        :param base: any base index.
        :return: The index of the corresponding codon, or ``None`` if
            the base does not fall in any codon.
        """

        if base in self._bases: return self._bases[base][1]
        else: return None

    ####################################################################

    def codon_position(self, base):
        
        """
        Tell if the given base is the 1st, 2nd or 3rd position of the
        codon in which it falls.

        :param base: any base index.
        :return: The index of the base in the codon (0, 1 or 3), or
            ``None`` if the base does not fall in any codon.
        """

        if base in self._bases: return self._codons[self._bases[base][1]].index(base)
        else: return None

    ####################################################################

    def codon_bases(self, codon):
        
        """
        Give the position of the three bases of a given codon. One or
        two positions (but never the middle one alone) will be ``None``
        if the codon is truncated (beginning/end of an exon without
        coverage of the previous/next one).
        
        :param codon: any codon index.
        :return: A tuple with the three base positions, potentially
            containing one or two ``None``, or, instead of the tuple,
            ``None`` if the codon index is out of range.
        """

        try: return self._codons[codon]
        except IndexError: return None

    ####################################################################

    def iter_exon_bounds(self):

        """
        This iterator returns ``(start, stop)`` tuples of the positions
        of the limits of each exon.
        """

        for start, stop in zip(self._starts, self._stops): yield start, stop

    ####################################################################

    def iter_codons(self, skip_partial=False):

        """
        This iterator returns ``(first, second, third)`` tuples of the
        positions of the three bases of each codon. If *skip_partial* is
        ``False``, partial codons (containing one or two ``None`` for
        bases that are in non-covered exons) are included.

        :param skip_partial: tells if codons containing one or two
            non-represented bases should be included.
        """

        for codon in self._codons:
            if (not skip_partial) or (None not in codon): yield codon

########################################################################
