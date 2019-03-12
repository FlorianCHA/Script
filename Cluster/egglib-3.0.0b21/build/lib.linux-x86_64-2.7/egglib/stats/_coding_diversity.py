"""
This module contains the :class:`.CodingDiversity` class.
"""

__license__ = """
    Copyright 2015-2016 Stephane De Mita, Mathieu Siol

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

from .. import _eggwrapper, _interface, tools, misc
from ..tools import _code_tools
from . import _filter, _site

########################################################################

class CodingDiversity(object):

    """
    This class processes alignments with a reading frame specification
    in order to detect synonymous and non-synonymous variable positions.
    It provides basic statistics, but it can also filter data to let the
    user compute all other statistics on synonymous-only, or
    non-synonymous-only variation (e.g. :math:`\pi` or ``D``).

    The constructor takes optional arguments. By default, build an empty
    instance. If arguments are passed, they must match the signature of
    :meth:`~.CodingDiversity.process` that will be called.

    The method :meth:`~.CodingDiversity.process` does all the work. Once
    it is called, data are available through the different instance
    attributes, and it is possible to generate alignments containing
    only codon sites with either one synonymous or one non-synonymous
    mutation. It is also possible to iterate over sites of both kind. In
    both cases, the generated data contains only codons, where each
    codon is represented by a single integer (see the methods
    :func:`.tools.int2codon`). These data can be analysed in the module
    :mod:`.stats` using the pre-defined :class:`.Filter` instance
    :attr:`.stats.filter_codon`.

    Note that, currently, the outgroup is ignored.
    """

    ####################################################################

    def __init__(self, *args, **kwargs):

        # internal processing helpers
        self._coding_site_stock = [_eggwrapper.CodingSite() for i in xrange(12)] # create a dozen sites to start with
        self._site1 = _eggwrapper.SiteHolder(0)
        self._site2 = _eggwrapper.SiteHolder(0)
        self._site3 = _eggwrapper.SiteHolder(0)
        self._sites_S = []
        self._sites_NS = []
        self._alleles_S = []
        self._alleles_NS = []
        self._static_frame = tools.ReadingFrame()
        self._filter = _filter.filter_dna._obj

        # run process
        if len(args) + len(kwargs) > 0: self.process(*args, **kwargs)

        # set internal variables to default values
        else:
            self._num_tot = 0
            self._num_eff = 0
            self._num_NS = 0.0
            self._num_S = 0.0
            self._num_pol_single = 0
            self._num_multiple_hits = 0
            self._num_multiple_alleles = 0
            self._num_pol_NS = 0
            self._num_pol_S = 0
            self._num_stop = 0
            self._groups = 0, [], []

    ####################################################################

    def _get_coding_site(self):
        if len(self._coding_site_stock) == 0: return _eggwrapper.CodingSite()
        else: return self._coding_site_stock.pop()

    ####################################################################

    def process(self, align, frame=None, struct=None, code=1,
            max_missing=0.0, consider_outgroup_missing=False,
            skipstop=True, raise_stop=False,
            multiple_alleles=False, multiple_hits=False,
            allow_alt=False):

        """
        Process an alignment. It this instance already had data in
        memory, they will all be erased.

        :param align: a :class:`.Align` instance containing the coding
            sequence to process. All sequences must be proper nucleotide
            sequences (only upper case characters). Standard ambiguity
            characters and alignment gaps are supported (but are treated
            as missing data) and unrecognized characters cause an error.
            DNA sequences (A, C, G and T as exploitable characters) are
            expected, as encoded using the ASCII table (which is the
            default for data imported from fasta files).

        :param frame: a :class:`.tools.ReadingFrame` instance specifying
            the position of coding sequence fragments (that is, exons
            without UTR) and the reading frame, or ``None``. It is legal
            to pass a reading frame specification with either terminal
            or internal gaps (only full codons will be considered). If
            ``None``, assume that coding sequences have been provided
            (the whole sequence can be translated in the first reading
            frame without interruption).

        :param struct: a :class:`.Structure` instance defining the
            structure to analyze. If ``None``, no structure is used (all
            samples are placed in a single population). If not ``None``,
            the passed :class:`Structure` instance must contain a
            population level or an individual level, or both. If both,
            they are required to be nested (individuals are mapped to
            populations). If individuals are specified, the outgroup
            samples must be mapped also, otherwise the ougroup will be
            ignored.

        :param code: genetic code identifier (see :ref:`genetic-codes`).
            Required to be an integer among the valid values. The
            default value is the standard genetic code.

        :param max_missing: maximum proportion of missing data to allow
            (including stop codons if *skipstop* if ``true``). By
            default, all codons with any missing data are excluded.

        :param consider_outgroup_missing: if ``True``, outgroup samples are
            included in the count for missing data (by default, outgroup
            samples are not considered). Only considered if an :class:`.Align`
            is passed.

        :param skipstop: if ``True``, stop codons are treated as missing
            data and skipped. If so, potential mutations to stop codons
            are not taken into account when estimating the number of
            non-synonymous sites. Warning (this may be
            counter-intuitive): it actually assumes that stop codons are
            not biologically plausible and considers them as missing
            data. On the other hand, if *skipstop* is ``False``, it
            considers stop codons as if they were valid amino acids.
            This option has no effect if *raise_stop* is ``True``.

        :param raise_stop: raise a :exc:`~exceptions.ValueError` if a
            stop codon is met. If ``True``, *skipstop* has no effect.
            This method does allow to ensure that stop codons are not
            present in an alignment, but only that stop codons are not
            present in any of the considered sites.

        :param multiple_alleles: include coding sites with more than
            two alleles (regardless of whether mutations hit the same
            position within the triplet). If there are more than two
            different codons, they must either encode for the same amino acid
            or all encode for different amino acids (otherwise the site is
            rejected). If more than one of the three codon position are
            mutated, the option *multiple_hits* is considered.

        :param multiple_hits: include coding sites for which more than
           one of the three positions has changed (regardless of the
           number of alleles). If there are more than two alleles, the
           option *multiple_alleles* is also considered.

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
        """

        # check code
        if code not in _code_tools._codes: raise ValueError, 'unknown genetic code: {0}'.format(code)
        code = _code_tools._codes[code]

        # if raise_stop is on: force skipstop to be false to make them appear
        if raise_stop and skipstop: skipstop = False

        # return all CodingSite's to stock
        self._coding_site_stock.extend(self._sites_S)
        self._coding_site_stock.extend(self._sites_NS)
        del self._sites_S[:]
        del self._sites_NS[:]
        del self._alleles_S[:]
        del self._alleles_NS[:]

        # ensure basic consistency
        if frame is None:
            frame = tools.ReadingFrame([(0, align.ls, 1)])
        elif frame.num_needed_bases > align.ls:
            raise ValueError, 'alignment is not matching the reading frame (not long enough)'

        # get structure (in case NS- or S-Align are requested)
        if struct == None:
            self._pl = 1
            self._ni = align.ns
            self._no = align.no
        else:
            struct = struct._obj
            if struct.get_ns_req() > align._obj.get_nsam_i(): raise ValueError, 'structure does not match the alignment'
            if struct.get_no_req() > align._obj.get_nsam_o(): raise ValueError, 'structure does not match the alignment'
            self._pl = struct.get_ploidy()
            self._ni = struct.num_indiv_ingroup()
            self._no = struct.num_indiv_outgroup()
        tot = self._ni + self._no

        # get the first coding site
        current = self._get_coding_site()

        # initialize variables
        self._num_tot = frame.num_full_codons
        self._num_eff = 0
        self._num_stop = 0
        self._num_NS = 0.0
        self._num_S = 0.0
        self._num_pol_single = 0
        self._num_multiple_hits = 0
        self._num_multiple_alleles = 0
        self._num_pol_NS = 0
        self._num_pol_S = 0

        if consider_outgroup_missing: max_missing = int(max_missing * self._pl * (self._ni+self._no))
        else: max_missing = int(max_missing * self._pl * self._ni)

        if allow_alt:
            first_codon = {}
            had_first_codon = {}
            for i in xrange(self._ni):
                first_codon[i] = dict.fromkeys(xrange(self._pl), False)
                had_first_codon[i] = dict.fromkeys(xrange(self._pl), False)

        # process all codon sites
        for idx, (i1, i2, i3) in enumerate(frame.iter_codons(True)):

            # generate three Site's (don't care about missing data)
            self._site1.reset(self._pl)
            self._site2.reset(self._pl)
            self._site3.reset(self._pl)
            self._site1.process_align(align._obj, i1, struct, self._filter, tot, True)
            self._site2.process_align(align._obj, i2, struct, self._filter, tot, True)
            self._site3.process_align(align._obj, i3, struct, self._filter, tot, True)

            # detect first non-missing codon if we are to allow alternative starts
            if allow_alt:
                for i in xrange(self._ni):
                    for j in xrange(self._pl):
                        if not had_first_codon[i][j] and (
                            self._site1.get_allele(self._site1.get_i(i, j)) != 45 or
                            self._site2.get_allele(self._site1.get_i(i, j)) != 45 or
                            self._site3.get_allele(self._site1.get_i(i, j)) != 45):
                                had_first_codon[i][j] = True
                                first_codon[i][j] = True

            # assemble the three positions (and check the number of missing data)
            good = current.process(self._site1, self._site2, self._site3, code, skipstop, max_missing)

            # check stop codon
            if current.nstop() > 0:
                if raise_stop == True: raise ValueError, 'stop codon found in sequences'
                self._num_stop += 1 # incremented even if skipstop is False

            # skip if too many missing data (but still count stop codon)
            if not good: continue

            # check number of alleles (exclude if no alleles)
            na = current.codons().get_nall_ing()
            if na < 1:
                continue

            if na == 2:
                self._num_pol_single += 1 # decremented below if the two alleles differ at more than one position

            if na > 2:
                self._num_multiple_alleles += 1
                good &= multiple_alleles

            # check for multiple hits
            if (self._site1.get_nall_ing() > 1) + (self._site2.get_nall_ing() > 1) + (self._site3.get_nall_ing() > 1) > 1:
                self._num_multiple_hits += 1
                good &= multiple_hits
                if na == 2: self._num_pol_single -= 1

            # exclude site if too multiple alleles/hits
            if not good: continue

            # save number of sites (because if we keep sites, we overwrite current
            na_tot = current.codons().get_nall()
            NSsites = current.NSsites()
            Ssites = current.Ssites()

            # process alleles
            if na > 1:
                aas = set()
                for i in xrange(self._ni): # cannot iterate over alleles because we need to know about first codon (for alt start)
                    for j in xrange(self._pl):
                        cdn = current.codons().get_allele(current.codons().get_i(i, j))
                        if cdn != _eggwrapper.MISSINGDATA:
                            if allow_alt and first_codon[i][j]:
                                first_codon[i][j] = False
                                if cdn != 35 and code.start(cdn):
                                    aas.add('M')
                                else:
                                    aas.add(code.aminoacid(cdn))
                            else:
                                aas.add(code.aminoacid(cdn))
                if len(aas) == 1: SYN = True
                elif len(aas) == na: SYN = False
                elif len(aas) < na: continue # skip because ambiguous syn/non-syn
                else: raise RuntimeError, 'unexpected error in CodingDiversity'

                if SYN:
                    self._num_pol_S += 1
                    self._sites_S.append(current)
                    self._alleles_S.append([current.codons().get_allele(i) for i in xrange(na_tot)])
                else:
                    self._num_pol_NS += 1
                    self._sites_NS.append(current)
                    self._alleles_NS.append([current.codons().get_allele(i) for i in xrange(na_tot)])

                # generate a new `current` site
                current = self._get_coding_site()

            # count site as exploitable
            self._num_eff += 1
            self._num_NS += NSsites
            self._num_S += Ssites

    ####################################################################

    @property
    def num_codons_tot(self):

        """
        Total number of considered codon sites that have been co. Only
        complete codons have been considered, but this value includes
        codons that have been rejected because of missing data.
        """

        return self._num_tot

    ####################################################################

    @property
    def num_codons_eff(self):

        """
        Number of codon sites that have been analysed (like
        :attr:`.num_codons_tot` but excluding sites rejected because of
        missing data).
        """

        return self._num_eff

    ####################################################################

    @property
    def num_codons_stop(self):

        """
        Number of codon sites (among those that have been analysed) with
        at least one codon stop in them.
        """

        return self._num_stop

    ####################################################################

    @property
    def num_sites_NS(self):

        """
        Estimated number of non-synonymous sites. Note that the total
        number of sites per codon is always 3.
        """

        return self._num_NS

    ####################################################################

    @property
    def num_sites_S(self):

        """
        Estimated number of synonymous sites. Note that the total number
        of sites per codon is always 3.
        """

        return self._num_S

    ####################################################################

    @property
    def num_pol_single(self):

        """
        Number of polymorphic coding sites with only one mutation. All these
        sites are always included.
        """

        return self._num_pol_single

    ####################################################################

    @property
    def num_multiple_alleles(self):

        """
        Number of polymorphic coding sites with more than two alleles. These
        sites are included only if *multiple_alleles* is ``True`` except those
        who mix synonymous and non-synonymous changes (they can be rejected
        if there are more than two alleles in total as well).
        """

        return self._num_multiple_alleles

    ####################################################################

    @property
    def num_multiple_hits(self):

        """
        Number of polymorphic codons for which more than one position is
        changed. These sites are included only if *multiple_hits* is ``True``
        and depdenting on the total number of alleles.
        """

        return self._num_multiple_hits

    ####################################################################

    @property
    def num_pol_NS(self):

        """
        Number of polymorphic codon sites with only one non-synonymous
        mutation.
        """

        return self._num_pol_NS

    ####################################################################

    @property
    def num_pol_S(self):

        """
        Number of polymorphic codon sites with only one synonymous
        mutation.
        """

        return self._num_pol_S

    ####################################################################

    def mk_align_S(self):

        """
        Create an :class:`.Align` instance with only synonymous codon
        sites. The alignment contains the same number of ingroup and
        outgroup samples as the original alignment, and a number of
        sites equal to :attr:`.num_pol_S`. Note that the returned
        alignment does not have group labels. These data can be analysed in the module
        :mod:`.stats` using the pre-defined :class:`.Filter` instance
        :attr:`.stats.filter_codon`.
        """

        align = _interface.Align(self._ni * self._pl, self._no * self._pl, self._num_pol_S)
        for sit, (site, alleles) in enumerate(zip(self._sites_S, self._alleles_S)):
            codons = site.codons()
            for idv in xrange(self._ni):
                for chm in xrange(self._pl):
                    X = codons.get_i(idv, chm)
                    if X == _eggwrapper.MISSING: X = 64
                    else: X = alleles[X]
                    align.set_i(idv * self._pl + chm, sit, X)
            for idv in xrange(self._no):
                for chm in xrange(self._pl):
                    X = codons.get_o(idv, chm)
                    if X == _eggwrapper.MISSING: X = 64
                    else: X = alleles[X]
                    align.set_o(idv * self._pl + chm, sit, X)

        return align

    ####################################################################

    def mk_align_NS(self):

        """
        Create an :class:`.Align` instance with only non-synonymous
        codon sites. Mostly similar to the method :meth:`.mk_align_S`.
        """

        align = _interface.Align(self._ni * self._pl, self._no * self._pl, self._num_pol_NS)
        for sit, (site, alleles) in enumerate(zip(self._sites_NS, self._alleles_NS)):
            codons = site.codons()
            for idv in xrange(self._ni):
                for chm in xrange(self._pl):
                    X = codons.get_i(idv, chm)
                    if X == _eggwrapper.MISSING: X = 64
                    else: X = alleles[X]
                    align.set_i(idv * self._pl + chm, sit, X)
            for idv in xrange(self._no):
                for chm in xrange(self._pl):
                    X = codons.get_o(idv, chm)
                    if X == _eggwrapper.MISSING: X = 64
                    else: X = alleles[X]
                    align.set_o(idv * self._pl + chm, sit, X)

        return align

    ####################################################################

    def iter_S(self):

        """
        Iterate over synonymous sites. Proposed as a more performant
        alternative to :meth:`.mk_align_S`. This method returns a
        generator that can be used in expressions such as::

            cs = egglib.stats.ComputeStats()
            cs.add_stat('thetaW')
            thetaW = 0.0
            cdiv = egglib.stats.CodingDiversity(align, frame)
            for site in cdiv.iter_S():
                stats = cs.process_site(site)
                thetaW += stats['thetaW']

        .. warning::

            This method returns repetively the same
            :class:`.SiteFrequency` with information updated at each
            round. The reason for this is performance. Never keep a
            reference to the iterator variable (that's ``site`` in the
            example above) outside the loop body.
        """

        site = _site.Site()
        obj = site._obj
        for coding_site in self._sites_S:
            site._obj = coding_site.codons()
            yield site
        site._obj = obj

    ####################################################################

    def iter_NS(self):

        """
        Iterate over non-synonymous sites. Mostly similar to the method
        :meth:`.iter_S` (see important warning about the fact that
        returned values are actually always the same
        :class:`.SiteFrequency` instance that is updated at each
        iteration round).
        """


        site = _site.Site()
        obj = site._obj
        for coding_site in self._sites_NS:
            site._obj = coding_site.codons()
            yield site
        site._obj = obj
