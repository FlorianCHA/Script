__license__ = """
    Copyright 2015-2018 Stephane De Mita, Mathieu Siol

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

from .. import _eggwrapper
from . import _filter, _freq, _structure, _site
from ..io import _vcf

class ComputeStats(object):

    """
    This class allows customizable and efficient analysis of diversity
    data. It is designed to minimize redundancy of underlying analyses
    so it is best to compute as many statistics as possible with a
    single instance. It also takes advantage of the object reuse policy,
    improving the efficiency of analyses when several (and especially
    many) datasets are examined in a row by the same instance of
    :class:`~.ComputeStats`.

    The constructor takes arguments that are automatically passed to the
    :meth:`configure` method.

    Statistics to compute are set using the method :meth:`~.ComputeStats.add_stats`
    which allows specifying several statistics at once and which can
    also be called several times to add more statistics to compute.
    """

    def list_stats(self):

        """
        Returns a list of tuples giving, for each available statistic,
        its code and a short description.
        """

        return [(k, v[0]) for (k, v) in self._stats.iteritems()]

    def _get_site_from_store(self):
        if len(self._site_store): return self._site_store.pop()
        else: return _eggwrapper.SiteHolder(0)

    def __init__(self, *args, **kwargs):

        self._sd = _eggwrapper.SiteDiversity()
        self._as = _eggwrapper.AlleleStatus()
        self._d1 = _eggwrapper.Diversity1()
        self._d2 = _eggwrapper.Diversity2()
        self._h = _eggwrapper.Haplotypes()
        self._struct_local = _eggwrapper.StructureHolder()
        self._struct_auxiliary = _eggwrapper.StructureHolder()
        self._rd = _eggwrapper.Rd()
        self._ld = _eggwrapper.MatrixLD()
        self._cv = _eggwrapper.ComputeV()
        self._frq = _freq.Freq()
        self._struct = _eggwrapper.StructureHolder()
        self._site = _eggwrapper.SiteHolder(0)
        self._site_cache = []
        self._site_store = []

        self._stats = {
            # code        # description                                                 # flag  # toggler
            'ns_site':    ('Number of analyzed samples per site',                       'sd',   []), 
            'Aing':       ('Number of alleles in ingroup',                              'sd',   []),
            'Aotg':       ('Number of alleles in outgroup',                             'sd',   []),
            'Atot':       ('Number of alleles in whole dataset',                        'sd',   []),
            'As':         ('Number of singleton alleles',                               'sd',   []),
            'Asd':        ('Number of singleton alleles (derived)',                     'sd',   []),
            'R':          ('Allelic richness',                                          'sd',   []),
            'thetaIAM':   ('Theta estimator based on He & IAM model',                   'sd',   []),
            'thetaSMM':   ('Theta estimator based on He & SMM model',                   'sd',   []),
            'He':         ('Expected heterozygosity',                                   'sd',   []),
            'Ho':         ('Observed heterozygosity',                                   'sd',   []),
            'Hi':         ('Inter-individual heterozygosity',                           'sd',   []),
            'Fis':        ('Inbreeding coefficient',                                    'sd',   []),
            'MAF':        ('Relative frequency of second most frequent allele',         'sd',   []),
            'MAF_pop':    ('MAF for each pop, as a list (check order of pops)',         'sd',   []),

            'WCst':       ('Weir and Cockerham for haploid data',                       'sd',   [self._sd.toggle_fstats_haplo]),
            'WCist':      ('Weir and Cockerham for diploid data',                       'sd',   [self._sd.toggle_fstats_diplo]),

            'WCisct':     ('Weir and Cockerham for hierarchical structure',             'sd',   [self._sd.toggle_fstats_hier]),

            'Dj':         ('Jost\'s D',                                                 'sd',   [self._sd.toggle_hstats]),
            'Hst':        ('Hudson\'s Hst',                                             'sd',   [self._sd.toggle_hstats]),
            'Gst':        ('Nei\'s Gst',                                                'sd',   [self._sd.toggle_hstats]),
            'Gste':       ('Hedrick\'s Gst\'',                                          'sd',   [self._sd.toggle_hstats]),

            'numSp':      ('Number of population-specific alleles',                     'as',   []),
            'numSpd':     ('Number of population-specific derived alleles',             'as',   []),
            'numShA':     ('Number of shared alleles',                                  'as',   []),
            'numShP':     ('Number of shared segregating alleles',                      'as',   []),
            'numFxA':     ('Number of fixed alleles',                                   'as',   []),
            'numFxD':     ('Number of fixed differences',                               'as',   []),
            'numSp*':     ('Sites with at least 1 pop-specific allele',                 'as',   []),
            'numSpd*':    ('Sites with at least 1 pop-specific derived allele',         'as',   []),
            'numShA*':    ('Sites with at least 1 shared allele',                       'as',   []),
            'numShP*':    ('Sites with at least 1 shared segregating allele',           'as',   []),
            'numFxA*':    ('Sites with at least 1 fixed allele',                        'as',   []),
            'numFxD*':    ('Sites with at least 1 fixed difference',                    'as',   []),

            'lseff':      ('Number of analysed sites',                                  'd1',   []),
            'nsmax':      ('Maximal number of available samples per site',              'd1',   []),
            'S':          ('Number of segregating sites',                               'd1',   []),
            'Ss':         ('Number of sites with only one singleton allele',            'd1',   []),
            'eta':        ('Minimal number of mutations',                               'd1',   []),
            'Pi':         ('Nucleotide diversity',                                      'd1',   []),

            'sites':      ('Index of polymorphic sites',                                'd1',   [self._d1.toggle_site_lists]),
            'singl':      ('Index of sites with at least one singleton allele',         'd1',   [self._d1.toggle_site_lists]),

            'lseffo':     ('Number of analysed orientable sites',                       'd1',   [self._d1.toggle_ori_site]),
            'nsmaxo':     ('Maximal number of available samples per orientable site',   'd1',   [self._d1.toggle_ori_site]),
            'sites_o':    ('Index of orientable polymorphic sites',                     'd1',   [self._d1.toggle_site_lists, self._d1.toggle_ori_site]),
            'singl_o':    ('Index of sites with at least one singleton allele',         'd1',   [self._d1.toggle_site_lists, self._d1.toggle_ori_site]),
            'So':         ('Number of segregating orientable sites',                    'd1',   [self._d1.toggle_ori_site]),
            'Sso':        ('Number of orientable sites with only one singleton allele', 'd1',   [self._d1.toggle_ori_site]),
            'nsingld':    ('Number of derived singletons',                              'd1',   [self._d1.toggle_ori_site]),
            'etao':       ('Minimal number of mutations are orientable sites',          'd1',   [self._d1.toggle_ori_site]),
            'nM':         ('Number of sites available for MFDM test',                   'd1',   [self._d1.toggle_ori_site]),
            'pM':         ('P-value of MDFM test',                                      'd1',   [self._d1.toggle_ori_site]),

            'nseffo':     ('Average number of exploitable samples at orientable sites', 'd1',   [self._d1.toggle_ori_div]),
            'thetaPi':    ('Pi using orientable sites',                                 'd1',   [self._d1.toggle_ori_div]),
            'thetaH':     ('Fay and Wu\'s estimator of theta',                          'd1',   [self._d1.toggle_ori_div]),
            'thetaL':     ('Zeng et al.\'s estimator of theta',                         'd1',   [self._d1.toggle_ori_div]),
            'Hns':        ('Fay and Wu\'s H (unstandardized)',                          'd1',   [self._d1.toggle_ori_div]),
            'Hsd':        ('Fay and Wu\'s H (standardized)',                            'd1',   [self._d1.toggle_ori_div]),
            'E':          ('Zeng et al.\'s E',                                          'd1',   [self._d1.toggle_ori_div]),
            'Dfl':        ('Fu and Li\'s D',                                            'd1',   [self._d1.toggle_ori_div]),
            'F':          ('Fu and Li\'s F',                                            'd1',   [self._d1.toggle_ori_div]),

            'nseff':      ('Average number of exploitable samples',                     'd1',   [self._d1.toggle_basic]),
            'thetaW':     ('Watterson\'s estimator of theta',                           'd1',   [self._d1.toggle_basic]),
            'Dxy':        ('Pairwise distance (if two populations)',                    'd1',   [self._d1.toggle_basic]),
            'Da':         ('Net pairwise distance (if two populations)',                'd1',   [self._d1.toggle_basic]),
            'D':          ('Tajima\'s D',                                               'd1',   [self._d1.toggle_basic]),
            'Deta':       ('Tajima\'s D using eta instead of S',                        'd1',   [self._d1.toggle_basic]),
            'D*':         ('Fu and Li\'s D*',                                           'd1',   [self._d1.toggle_basic]),
            'F*':         ('Fu and Li\'s F*',                                           'd1',   [self._d1.toggle_basic]),

            'R2':         ('Ramos-Onsins and Rozas\'s R2 (using singletons)',           'd2',   [self._d2.toggle_singletons]),
            'R3':         ('Ramos-Onsins and Rozas\'s R3 (using singletons)',           'd2',   [self._d2.toggle_singletons]),
            'R4':         ('Ramos-Onsins and Rozas\'s R4 (using singletons)',           'd2',   [self._d2.toggle_singletons]),
            'Ch':         ('Ramos-Onsins and Rozas\'s Ch (using singletons)',           'd2',   [self._d2.toggle_singletons]),

            'R2E':        ('Ramos-Onsins and Rozas\'s R2E (using external singletons)', 'd2',   [self._d2.toggle_singletons]),
            'R3E':        ('Ramos-Onsins and Rozas\'s R3E (using external singletons)', 'd2',   [self._d2.toggle_singletons]),
            'R4E':        ('Ramos-Onsins and Rozas\'s R4E (using external singletons)', 'd2',   [self._d2.toggle_singletons]),
            'ChE':        ('Ramos-Onsins and Rozas\'s ChE (using external singletons)', 'd2',   [self._d2.toggle_singletons]),

            'B':          ('Wall\'s B statistic',                                       'd2',   [self._d2.toggle_partitions]),
            'Q':          ('Wall\'s Q statistic',                                       'd2',   [self._d2.toggle_partitions]),

            'K':          ('Number of haplotypes',                                      'h',    []),
            'Ke':         ('Number of haplotypes (only ingroup)',                       'h',    []),
            'Fst':        ('Hudson\'s Fst',                                             'h',    [self._toggle_haplotype_stats]),
            'Kst':        ('Hudson\'s Kst',                                             'h',    [self._toggle_haplotype_stats]),
            'Snn':        ('Hudson\'s nearest nearest neighbour statistic',             'h',    [self._toggle_haplotype_Snn]),

            'Fs':         ('Fu\'s Fs',                                                  'Fs',   []),

            'rD':         ('R_bar{d} statistic',                                        'rD',   []),

            'Rmin':       ('Minimal number of recombination events',                    'LD',   [self._ld.toggle_Rmin]),
            'RminL':      ('Number of sites used to compute Rmin',                      'LD',   [self._ld.toggle_Rmin]),
            'Rintervals': ('List of start/end positions of recombination intervals',    'LD',   [self._ld.toggle_Rmin]),
            'nPairs':     ('Number of allele pairs used for ZnS, Z*nS, and Z*nS*',      'LD',   [self._ld.toggle_stats]),
            'nPairsAdj':  ('Allele pairs at adjacent sites (used for ZZ and Za)',       'LD',   [self._ld.toggle_stats]),
            'ZnS':        ('Kelly et al.\'s ZnS',                                       'LD',   [self._ld.toggle_stats]),
            'Z*nS':       ('Kelly et al.\'s Z*nS',                                      'LD',   [self._ld.toggle_stats]),
            'Z*nS*':      ('Kelly et al.\'s Z*nS*',                                     'LD',   [self._ld.toggle_stats]),
            'Za':         ('Rozas et al.\'s Za',                                        'LD',   [self._ld.toggle_stats]),
            'ZZ':         ('Rozas et al.\'s ZZ',                                        'LD',   [self._ld.toggle_stats]),

            'V':          ('Allele size variance',                                      'V',    [])
        }

        self.clear_stats()
        self.configure(*args, **kwargs)

    def reset(self):

        """
        Reset all currenctly computed statistics (but keep the list of
        statistics to compute).
        """

        self._sd.reset()
        self._as.reset()
        self._d1.reset_stats()
        self._d2.reset()
        self._h.setup(None)
        self._ld.reset()
        self._cv.reset()
        self._rd.reset()
        self._site_index = 0
        self._static_sites = True # set to False if sites cannot be considered to be static (required for LD stats)
        self._site_store.extend(self._site_cache)
        del self._site_cache[:]
        self._d1.set_option_ns_set(_eggwrapper.UNKNOWN)
        self._structure = None # structure passed as argument to configure()

    def configure(self,
                  struct = None,
                  multi = False,
                  multi_hits = False,
                  ignore_ns = False,
                  set_ns = None,
                  LD_min_n = 2,
                  LD_max_maj = 1.0,
                  LD_multiallelic = 0,
                  LD_min_freq = 0,
                  Rmin_oriented = False):

        """
        Configure the instance. The values provided for parameters will
        affect all subsequent analyses. Calling this method resets all
        statistics.

        :param struct: a :class:`~Structure` instance describing the
            structure of objects that will be passed. This argument can
            be bypassed by the *struct* arguments of subsequents
            methods.

        :param multi: process several alignments, set of sites, or sites
            in a row and only yield statistics when :meth:`.results` is
            called (each of the methods considered will return ``None``).

        :param multi_hits: allow multiple mutations at the same site.

        :param ignore_ns: do not require that all sites have the same
            number of samples for computing Fay and Wu's statistics.

        :param set_ns: set the total number of samples (useful to compute
            Fay and Wu's statistics, only in the case when data contain
            missing data and data are provided as sites or as several
            alignments with different numbers of samples. This option is
            ignored if *ignore_ns* is ``False``.

        :param LD_min_n: Minimal number of non-missing samples. Allows to
            specify a more stringent filter than *max_missing*.
            Only considered for calculating Rozas *et al.*'s and Kelly's
            statistics (the most stringent of the two criteria applies).

        :param LD_max_maj: Maximal relative frequency of the main allele.
            Only considered for calculating Rozas *et al.*'s and Kelly's
            statistics.

        :param LD_multiallelic: One of 0 (ignore them), 1 (use main
            allele only), and 2 (use all possible pairs of alleles).
            Defines what is done for pairs of sites of which one or
            both have more than two alleles (while computing linkage
            disequilibrium). In case of option 2, a filter can be
            applied with option *LD_min_freq*. Only considered for
            calculating Rozas *et al.*'s and Kelly's statistics.

        :param LD_min_freq: Only considered if option 2 is used for
            *LD_multiallelic*. Only consider alleles that are in
            absolute frequency equal to or larger than the given value.
            Only considered for calculating Rozas *et al.*'s and Kelly's
            statistics

        :param Rmin_oriented: Only for computing Rmin: use only
            orientable sites.
        """

        self.reset()
        self._structure = struct
        if struct is not None:
            self._h.setup(struct._obj)
        if LD_min_n < 0: raise ValueError, 'LD_min_n is too small'
        self._LD_min_n = LD_min_n
        if LD_max_maj <= 0.0 or LD_max_maj > 1.0: raise ValueError, 'LD_max_maj out of range'
        self._LD_max_maj = LD_max_maj
        if LD_multiallelic < 0: raise ValueError, 'LD_multiallelic cannot be negative'
        try: self._LD_multiallelic = [_eggwrapper.MatrixLD.ignore,
                                      _eggwrapper.MatrixLD.use_main,
                                      _eggwrapper.MatrixLD.use_all][LD_multiallelic]
        except IndexError: raise ValueError, 'LD_multiallelic out of range'
        if LD_min_freq < 0: raise ValueError, 'LD_min_freq out of range'
        self._LD_min_freq = LD_min_freq
        self._Rmin_oriented = Rmin_oriented
        self._multi_hits = bool(multi_hits)
        self._multi = bool(multi)
        self._d1.set_option_multiple(multi_hits)
        self._d2.set_option_multiple(multi_hits)
        self._ignore_ns = bool(ignore_ns)
        if set_ns is None:
            self._set_ns = None
        else:
            if set_ns < 1: raise ValueError, 'set_ns out of range'
            self._set_ns = int(set_ns)

    def add_stats(self, *stats):

        """
        add_stats(stat, ...)

        Add one or more statistics to compute. Every statistic identifier must be among
        the list of available statistics, regardless of what data is to
        be analyzed. If statistics cannot be computed, they will be
        returned as ``None``.
        """

        for stat in stats:
            if stat not in self._stats:
                raise ValueError, 'invalid statistic: {0}'.format(stat)
            self._wanted_stats.add(stat)
            self._flags[self._stats[stat][1]] = True;
            for f in self._stats[stat][2]: f()

        if self._flags['Fs']:
            self._flags['d1'] = True
            self._flags['h'] = True

        if (self._flags['as'] or self._flags['d1'] or self._flags['d2']
                              or self._flags['rD'] or self._flags['LD']
                              or self._flags['h']):
            self._flags['sd'] = True

    def all_stats(self):

        """
        Add all possible statistics. Those who cannot be computed will
        be reported as ``None``. Also reset all currently computed statistics (if any).
        """

        self.add_stats(*self._stats)

    def clear_stats(self):

        """
        Clear the list of statistics to compute.
        """

        self._wanted_stats = set()
        self._flags = dict.fromkeys([v[1] for v in self._stats.itervalues()], False)
        self._sd.toggle_off()
        self._d1.toggle_off()
        self._ld.toggle_off()
        self._toggle_off()

    def _toggle_off(self):
        self._flag_haplotype_dist = False
        self._flag_haplotype_stats = False

    def _toggle_haplotype_stats(self):
        self._flag_haplotype_dist = True
        self._flag_haplotype_stats = True

    def _toggle_haplotype_Snn(self):
        self._flag_haplotype_dist = True

    def process_freq(self, frq, alleles=None, position=None):

        """
        Analyze a site based on already computed frequencies.

        :param frq: a :class:`.Freq` instance.
        :params alleles: list of allele values (as integers). Only used
            if statistic ``V`` is required and ignored otherwise.
            Statistic ``V`` is not computed if this argument is not
            specified. The length of the list must match the total
            number of alleles.
        :param position: position of the site. Must be an integer value.
            By default, use the site loading index.
        :return: A dictionary of statistics, unless *multi* is set.
        """

        if not self._multi:
            stats = dict.fromkeys(self._wanted_stats, None)

        if position is None: position = self._site_index

        if self._flags['sd'] == True:
            res = self._sd.process(frq._obj)
            if not self._multi: self._get_sd_stats(res, stats)

        if self._flags['as'] == True and (res & 1024) != 0 and self._sd.npop_eff1() > 1:
            self._as.process(frq._obj)
            if not self._multi: self._get_as_stats(stats)

        if self._flags['d1'] == True and (res & 2) != 0:
            if self._ignore_ns == False and self._d1.lt() == 0:
                if self._set_ns is None: self._d1.set_option_ns_set(frq.nseff(cpt=frq.ingroup))
                else: self._d1.set_option_ns_set(self._set_ns)
            self._d1.load(frq._obj, self._sd, position)

        if self._flags['V'] == True and alleles is not None:
            if len(alleles) != frq._obj.num_alleles():
                raise ValueError, 'invalid number of alleles'
            self._cv.setup_alleles(len(alleles))
            for i,v in enumerate(alleles): self._cv.set_allele(i, v)
            V = self._cv.compute(frq._obj.frq_ingroup())
            if not self._multi and V != _eggwrapper.UNDEF: stats['V'] = V

        self._site_index += 1
        if not self._multi: return stats

    def process_site(self, site, position=None, struct=None):

        """
        Analyze a site.

        :param site: a :class:`.Site` instance.
        :param position: position of the site. Must be an integer value.
            By default, use the site loadingdef _results index.
        :struct: a :class:`.Structure` instance describing population
            structure. By default, doesn't apply structure. It is not
            allowed to pass a structure several times (with *multi* set
            to ``True``). In that case, pass the structure to
            :meth:`configure`.
        :return: A dictionary of statistics, unless *multi* is set.
        """

        if site._obj.get_ploidy() == 0: raise ValueError, 'cannot process site with ploidy 0'

        if position is None: position = self._site_index
        if not self._multi: stats = dict.fromkeys(self._wanted_stats, None)

        if self._multi:
            if struct is not None:
                raise ValueError, 'cannot pass a structure if `multi` is active (use `configure()` instead)'
            struct = self._structure

        self._frq.process_site(site, struct=struct)

        if self._flags['sd'] == True:
            res = self._sd.process(self._frq._obj)
            if not self._multi: self._get_sd_stats(res, stats)
        if self._flags['as'] == True and (res & 1024) != 0 and self._sd.npop_eff1() > 1:
            self._as.process(self._frq._obj)
            if not self._multi: self._get_as_stats(stats)
        if self._flags['d1'] == True and (res & 2) != 0:
            if self._ignore_ns == False and self._d1.lt() == 0:
                if self._set_ns is None: self._d1.set_option_ns_set(site._obj.get_ning() * site._obj.get_ploidy())
                else: self._d1.set_option_ns_set(self._set_ns)
            self._d1.load(self._frq._obj, self._sd, position)
        if self._flags['d2'] == True and (res & 2) != 0:
            self._d2.load(site._obj, self._sd, self._frq._obj)
        if self._flags['h'] == True and self._sd.Aglob() > 1 and (self._sd.Aing() == 2 or self._multi_hits == True):
            self._h.load(site._obj)
        if self._flags['rD'] == True:
            self._rd.load(site._obj)
        if self._flags['V'] == True:
            self._cv.setup_alleles_from_site(site._obj)
            V = self._cv.compute(self._frq._obj.frq_ingroup())
            if not self._multi and V != _eggwrapper.UNDEF: stats['V'] = V
        self._site_index += 1

        if self._multi: self._static_sites = False
        else: return stats

    def process_sites(self, sites, positions=None, struct=None):

        """
        Analyze a list of sites.

        :param sites: a list, or other sequence of :class:`.Site`
            instance. Or a sliding window as object of the :class:`.VcfWindow`
            with a number of sites greater than 0.
        :param positions: a list, or other sequence of positions for all
            sites, or ``None``. If ``None``, use the index of each site. Otherwise,
            must be a sequences of integer values (length matching the
            length of *sites*).
        :struct: a :class:`.Structure` instance describing population
            structure. By default, doesn't apply structure. It is not
            allowed to pass a structure several times (with *multi* set
            to ``True``). In that case, pass the structure to
            :meth:`configure`.
        :return: A dictionary of statistics, unless *multi* is set.
        """

        if not self._multi: self.reset()
        if positions is None:
            positions = (self._site_index + i for i in xrange(len(sites)))
        elif len(positions) != len(sites):
            raise ValueError, 'number of positions must match the number of sites'

        if self._multi:
            if struct is not None:
                raise ValueError, 'cannot pass a structure if `multi` is active (use `configure()` instead)'
            struct = self._structure

        elif self._flags['h'] == True and struct is not None:
            self._h.setup(struct._obj)

        if len(sites) > 0:
            if self._flags['d1'] == True and self._ignore_ns == False and self._d1.lt() == 0:
                if self._set_ns is None: self._d1.set_option_ns_set(sites[0]._obj.get_ning() * sites[0]._obj.get_ploidy())
                else: self._d1.set_option_ns_set(self._set_ns)

            for site, pos in zip(sites, positions):
                self._frq.process_site(site, struct=struct)
                if self._flags['sd'] == True: res = self._sd.process(self._frq._obj)  # don't require consistency
                if self._flags['as'] == True and (res & 1024) != 0 and self._sd.npop_eff1() > 1: self._as.process(self._frq._obj)  # don't require consistency
                if self._flags['d1'] == True and (res & 2) != 0: self._d1.load(self._frq._obj, self._sd, pos)  # don't require consistency
                if self._flags['d2'] == True and (res & 2) != 0: self._d2.load(site._obj, self._sd, self._frq._obj)
                if self._flags['h'] == True and self._sd.Aglob() > 1 and (self._multi_hits == True or self._sd.Aing() < 3): self._h.load(site._obj)
                if self._flags['rD'] == True and (res & 2) != 0 and self._sd.Aing() > 1: self._rd.load(site._obj)
                if self._flags['LD'] == True and (self._sd.Aing() == 2 or (self._sd.Aing() > 2 and self._LD_multiallelic != _eggwrapper.MatrixLD.ignore)): self._ld.load(site._obj, pos)
                if self._flags['V'] == True:
                    self._cv.setup_alleles_from_site(site._obj)
                    V = self._cv.compute(self._frq._obj.frq_ingroup())

        self._site_index += len(sites)
        if self._multi: self._static_sites = False # prevent attempting using sites (which may be invalidated)
        else: return self.results()

    def process_align(self, align, positions=None, struct=None,
        filtr=None, max_missing=0.0, consider_outgroup_missing=False):

        """
        Analyze an alignment.

        :param align: an :class:`.Align` instance.
        :param positions: a list, or other sequence of positions for all
            sites, or ``None``. If ``None``, use the index of each site. Otherwise,
            must be a sequences of integer values (length matching the
            number of sites).
        :param struct: a :class:`.Structure` instance describing population
            structure.  *struct* can also be an integer specifying the level
            of population labels in the alignment labels table. Only
            populations can be specified this way. If ``None`` (default),
            don't use structure. It is not allowed to pass a structure
            several times (with *multi* set to ``True``). In that case,
            pass the structure to :meth:`configure`.
        :param filtr: :class:`~.stats.Filter` instance determining
            what allelic values are acceptable and what ones should be
            considered as missing (all other values causing an exception).
            By default, nucleotide sequences are accepted (including
            IUPAC ambiguity characters as missing data), with
            case-independent matching. The user may want to use a custom
            instance or use one of the predefined instances in the
            :mod:`~.stats` module).
        :param max_missing: Maximum proportion of missing data. The
            default is to exclude all sites wdef _resultsith missing data. Missing
            data include all ambiguity characters and alignment gaps.
            Sites not passing this threshold are ignored for all
            analyses.
        :param consider_outgroup_missing: if ``True``, take the outgroup
            into account in the max_missing threshold (by default, only
            ingroup samples are considered).
        :return: A dictionary of statistics, unless *multi* is set.
        """

        # reset if needed
        if not self._multi: self.reset()

        # check max_missing
        max_missing = float(max_missing)
        if max_missing < 0.0 or max_missing >= 1.0:
            raise ValueError, 'invalid value for `max_missing` argument: {0}'.format(max_missing)
        otg_missing = bool(consider_outgroup_missing)

        # get number of sites and positions
        ls = align._obj.get_nsit()
        if positions is None:
            positions = (self._site_index + i for i in xrange(ls))
        elif len(positions) != ls:
            raise ValueError, 'number of positions must match the number of sites'

        # get filtr
        if filtr is None: filtr = _filter.filter_dna

        # process structure
        if self._multi:
            if struct is not None:
                raise ValueError, 'cannot pass a structure if `multi` is active (use `configure()` instead)'
            struct = self._structure

        if struct is None:
            pl = 1
        else:
            if isinstance(struct, int):
                if struct < 0 or struct >= align._obj.get_ngroups(): raise ValueError, 'invalid group level index'
                self._struct_local.reset()
                self._struct_local.get_structure(align._obj, _eggwrapper.UNKNOWN, struct, _eggwrapper.UNKNOWN, 1, False)
                struct = self._struct_local
                pl = 1
            else:
                struct = struct._obj
                if struct.get_ns_req() > align._obj.get_nsam_i(): raise ValueError, 'structure does not match alignment'
                if struct.get_no_req() > align._obj.get_nsam_o(): raise ValueError, 'structure does not match alignment'
                if otg_missing: max_missing_num = max_missing * struct.get_ns()
                else: max_missing_num = max_missing * (struct.get_ns() + struct.get_no())
                pl = struct.get_ploidy()
            _structure.Structure._make_auxiliary(struct, self._struct_auxiliary)

        # set the structure of Haplotypes (after processing potential index-based struct)
        if self._flags['h']  and self._multi == False and struct is not None:
            self._h.setup(self._struct_auxiliary)

        # get max_missing number
        if otg_missing: max_missing_num = max_missing * (align._obj.get_nsam_i() + align._obj.get_nsam_o())
        else: max_missing_num = max_missing * align._obj.get_nsam_i()
        max_missing_num = int(max_missing_num)

        # get numbers of samples
        if struct is None:
            ni = align._obj.get_nsam_i()
            no = align._obj.get_nsam_o()
        else:
            ni = struct.num_indiv_ingroup()
            no = struct.num_indiv_outgroup()

        # initialize ns_set
        if self._flags['d1'] == True and self._ignore_ns == False and self._d1.lt() == 0:
            if self._set_ns is None: self._d1.set_option_ns_set(align._obj.get_nsam_i())
            else: self._d1.set_option_ns_set(self._set_ns)

        # process sites
        for idx, pos in enumerate(positions):
            self._site.reset(pl)
            if not self._site.process_align(align._obj, idx, struct, filtr._obj, max_missing_num, otg_missing):
                continue

            if struct is None:
                self._frq._obj.setup_raw(1, 1, no, pl)
                self._frq._obj.setup_pop(0, 0, 0, ni)
            else:
                self._frq._obj.setup_structure(self._struct_auxiliary, pl)

            self._frq._obj.process_site(self._site)
            if self._flags['sd'] == True: res = self._sd.process(self._frq._obj)
            if self._flags['as'] == True and (res & 1024) != 0 and self._sd.npop_eff1() > 1: self._as.process(self._frq._obj)
            if self._flags['d1'] == True and (res & 2) != 0: self._d1.load(self._frq._obj, self._sd, pos)
            if self._flags['d2'] == True and (res & 2) != 0: self._d2.load(self._site, self._sd, self._frq._obj)
            if self._flags['h'] == True and self._sd.Aglob() > 1 and (self._multi_hits == True or self._sd.Aing() < 3): self._h.load(self._site)
            if self._flags['rD'] == True and (res & 2) != 0 and self._sd.Aing() > 1: self._rd.load(self._site)
            if self._flags['LD'] == True and (self._sd.Aing() == 2 or (self._sd.Aing() > 2 and self._LD_multiallelic != _eggwrapper.MatrixLD.ignore)): self._ld.load(self._site, pos)
            if self._flags['V'] == True:
                self._cv.setup_alleles_from_site(self._site)
                V = self._cv.compute(self._frq._obj.frq_ingroup())

            # if need to keep sites, put it in cache and claim new
            if self._flags['LD']:
                self._site_cache.append(self._site)
                self._site = self._get_site_from_store()

        self._site_index += ls
        if not self._multi: return self.results()

    def results(self):

        """
        Return the value of statistics for all sites since the last call
        to this method, to :meth:`reset`, or any addition of statistics,
        or the object creation,
        whichever is most recent. For statistics that can not be
        computed, ``None`` is returned.
        """

        stats = dict.fromkeys(self._wanted_stats, None)

        if self._flags['sd'] == True:
            res = self._sd.average()
            self._get_sd_stats(res, stats)

        if self._flags['as'] == True:
            self._as.total()
            if (res & 1024) != 0 and self._sd.npop_eff1() > 1:
                self._get_as_stats(stats)
                if self._as.nsites() > 0:
                    if 'numSp*' in stats: stats['numSp*'] = self._as.Sp_T1()
                    if 'numSpd*' in stats and self._as.nsites_o() > 0: stats['numSpd*'] = self._as.Spd_T1()
                    if 'numShP*' in stats: stats['numShP*'] = self._as.ShP_T1()
                    if 'numShA*' in stats: stats['numShA*'] = self._as.ShA_T1()
                    if 'numFxD*' in stats: stats['numFxD*'] = self._as.FxD_T1()
                    if 'numFxA*' in stats: stats['numFxA*'] = self._as.FxA_T1()

        if self._flags['d1'] == True:
            res = self._d1.compute()
            if 'lseff'  in stats: stats['lseff'] =  self._d1.ls()
            if (res & 1) != 0:
                if 'nsmax'   in stats: stats['nsmax'] =    self._d1.nsmax()
                if 'S'       in stats: stats['S'] =        self._d1.S()
                if 'Ss'      in stats: stats['Ss'] =       self._d1.Ss()
                if 'sites'   in stats: stats['sites'] =    [self._d1.site(i) for i in xrange(self._d1.S())]
                if 'singl'   in stats: stats['singl'] =    [self._d1.singl(i) for i in xrange(self._d1.Ss())]
                if 'eta'     in stats: stats['eta'] =      self._d1.eta()
                if 'Pi'      in stats: stats['Pi'] =       self._d1.Pi()
            if (res & 4) != 0:
                if 'lseffo'  in stats: stats['lseffo'] =   self._d1.lso()
            if (res & 8) != 0:
                if 'nsmaxo'  in stats: stats['nsmaxo'] =   self._d1.nsmaxo()
                if 'So'      in stats: stats['So'] =       self._d1.So()
                if 'Sso'     in stats: stats['Sso'] =      self._d1.Sso()
                if 'nsingld' in stats: stats['nsingld'] =  self._d1.nsingld()
                if 'sites_o' in stats: stats['sites_o'] =  [self._d1.site_o(i) for i in xrange(self._d1.So())]
                if 'singl_o' in stats: stats['singl_o'] =  [self._d1.singl_o(i) for i in xrange(self._d1.Sso())]
                if 'etao'    in stats: stats['etao'] =     self._d1.etao()
                if 'nM'      in stats: stats['nM'] =       self._d1.nM()
            if (res & 16) != 0:
                if 'pM'      in stats: stats['pM'] =       self._d1.pM()
            if (res & 32) != 0:
                if 'nseffo'  in stats: stats['nseffo'] =   self._d1.nseffo()
                if 'thetaPi' in stats: stats['thetaPi'] =  self._d1.thetaPi()
                if 'thetaH'  in stats: stats['thetaH'] =   self._d1.thetaH()
                if 'thetaL'  in stats: stats['thetaL'] =   self._d1.thetaL()
                if 'Hns'     in stats: stats['Hns'] =      self._d1.Hns()
            if (res & 1024) != 0:
                if 'Hsd'     in stats: stats['Hsd'] =      self._d1.Hsd()
            if (res & 2048) != 0:
                if 'E'       in stats: stats['E'] =        self._d1.E()
            if (res & 4096) != 0:
                if 'Dfl'     in stats: stats['Dfl'] =      self._d1.Dfl()
            if (res & 8192) != 0:
                if 'F'       in stats: stats['F'] =        self._d1.F()
            if (res & 128) != 0:
                if 'nseff'   in stats: stats['nseff'] =    self._d1.nseff()
                if 'thetaW'  in stats: stats['thetaW'] =   self._d1.thetaW()
            if (res & 16384) != 0:
                if 'Dxy'     in stats: stats['Dxy'] =      self._d1.Dxy()
                if 'Da'      in stats: stats['Da'] =       self._d1.Da()
            if (res & 256) != 0:
                if 'F*'      in stats: stats['F*'] =       self._d1.Fstar()
            if (res & 512) != 0:
                if 'D'       in stats: stats['D'] =        self._d1.D()
                if 'Deta'    in stats: stats['Deta'] =     self._d1.Deta()
                if 'D*'      in stats: stats['D*'] =       self._d1.Dstar()

        if self._flags['d2'] == True:
            res = self._d2.compute()
            if (res & 1) != 0:
                if (res & 2) != 0: pass
                elif (res & 4) != 0: raise ValueError, 'inconsistent number of samples across sites'
                elif (res & 8) != 0: raise ValueError, 'inconsistent ploidy across sites'
                elif (res & 16) != 0: raise ValueError, 'inconsistent frequencies set provided for site'
                elif (res & 32) != 0: raise ValueError, 'inconsistent site diversity provided'
                else: raise RuntimeError, 'unknown error occurred but error flag is unknown (status flag is: {0})'.format(res)
            else:
                if (res & 256) != 0:
                    if 'R2' in stats: stats['R2'] = self._d2.R2()
                    if 'R3' in stats: stats['R3'] = self._d2.R3()
                    if 'R4' in stats: stats['R4'] = self._d2.R4()
                    if 'Ch' in stats: stats['Ch'] = self._d2.Ch()
                if (res & 512) != 0:
                    if 'R2E' in stats: stats['R2E'] = self._d2.R2E()
                    if 'R3E' in stats: stats['R3E'] = self._d2.R3E()
                    if 'R4E' in stats: stats['R4E'] = self._d2.R4E()
                    if 'ChE' in stats: stats['ChE'] = self._d2.ChE()
                if (res & 1024) != 0:
                    if 'B' in stats: stats['B'] = self._d2.B()
                    if 'Q' in stats: stats['Q'] = self._d2.Q()

        if self._flags['h'] == True and self._h.n_sites() > 0:
            self._h.cp_haplotypes()
            if 'K' in stats and (self._h.ne_ing() > 0 or self._h.ne_otg() > 0): stats['K'] = self._h.ng_hapl()
            if self._h.ne_ing() > 0 and 'Ke' in stats: stats['Ke'] = self._h.ni_hapl()
            if self._flag_haplotype_dist == True:
                self._h.cp_dist()
                if self._flag_haplotype_stats == True:
                    res = self._h.cp_stats()
                    if 'Fst' in stats and res > 0: stats['Fst'] = self._h.Fst()
                    if 'Kst' in stats and res > 1: stats['Kst'] = self._h.Kst()
                if 'Snn' in stats and self._h.nstot() > 1: stats['Snn'] = self._h.Snn()

        if 'Fs' in stats and self._h.n_sites() and self._h.ni_hapl() > 0:
            stats['Fs'] = _eggwrapper.Fs(self._h.ne_ing(), self._h.ni_hapl(), self._d1.Pi())

        if self._static_sites == True and 'rD' in stats:
            stats['rD'] = self._rd.compute()
            if stats['rD'] == _eggwrapper.UNDEF: stats['rD'] = None

        if self._static_sites == True and self._flags['LD'] == True:
            flag = self._ld.process(self._LD_min_n, self._LD_max_maj,
                                    self._LD_multiallelic, self._LD_min_freq,
                                    self._Rmin_oriented)

            if 'nPairs' in stats: stats['nPairs'] = self._ld.num_allele_pairs()
            if 'nPairsAdj' in stats: stats['nPairsAdj'] = self._ld.num_allele_pairs_adj()
            if 'RminL' in stats: stats['RminL'] = self._ld.Rmin_num_sites()
            if (flag&1) != 0:
                if 'ZnS' in stats: stats['ZnS'] = self._ld.ZnS()
                if 'Z*nS' in stats: stats['Z*nS'] = self._ld.ZnS_star1()
                if 'Z*nS*' in stats: stats['Z*nS*'] = self._ld.ZnS_star2()
            if (flag&2) != 0:
                if 'Za' in stats: stats['Za'] = self._ld.Za()
                if 'ZZ' in stats: stats['ZZ'] = self._ld.ZZ()
            if (flag&4) != 0:
                if 'Rmin' in stats: stats['Rmin'] = self._ld.Rmin()
                if 'Rintervals' in stats: stats['Rintervals'] = [(self._ld.Rmin_left(i), self._ld.Rmin_right(i)) for i in xrange(self._ld.Rmin())]

        if 'V' in stats and self._cv.num_sites() > 0:
            stats['V'] = self._cv.average()

        self.reset()
        return stats

    def _get_sd_stats(self, flag, stats):
        if (flag & 1) != 0 and 'ns_site' in stats:
            if (flag & 512) != 0: stats['ns_site'] = int(self._sd.ns())
            else: stats['ns_site'] = self._sd.ns()
        if (flag & 2) != 0:
            if 'Aing' in stats:
                if (flag & 512) != 0: stats['Aing'] = int(self._sd.Aing())
                else: stats['Aing'] = self._sd.Aing()
            if 'Atot' in stats:
                if (flag & 512) != 0: stats['Atot'] = int(self._sd.Aglob())
                else: stats['Atot'] = self._sd.Aglob()
            if 'As' in stats:
                if (flag & 512) != 0: stats['As'] = int(self._sd.S())
                else: stats['As'] = self._sd.S()
            if 'R' in stats: stats['R'] = self._sd.R()
            if 'He' in stats: stats['He'] = self._sd.He()

        if (flag & 4) != 0:
            if 'thetaIAM' in stats: stats['thetaIAM'] = self._sd.thetaIAM()
            if 'thetaSMM' in stats: stats['thetaSMM'] = self._sd.thetaSMM()

        if (flag & 8) != 0:
            if 'Ho' in stats: stats['Ho'] = self._sd.Ho()
            if 'Hi' in stats: stats['Hi'] = self._sd.Hi()
            if 'Fis' in stats and self._sd.He() > 0.0:
                stats['Fis'] = 1 - self._sd.Ho() / self._sd.Hi()

        if (flag & 16) != 0:
            if 'Asd' in stats:
                if (flag & 512) != 0: stats['Asd'] = int(self._sd.Sd())
                else: stats['Asd'] = self._sd.Sd()
        if (flag & 2048) != 0:
            if 'Aotg' in stats:
                if (flag & 512) != 0: stats['Aotg'] = int(self._sd.Aout())
                else: stats['Aotg'] = self._sd.Aout()

        if (flag & 32) != 0 and 'WCst' in stats:
            n = self._sd.n()
            d = self._sd.d()
            stats['WCst'] = n/d if d > 0.0 else None

        if (flag & 64) != 0 and 'WCist' in stats:
            a = self._sd.a()
            b = self._sd.b()
            c = self._sd.c()
            stats['WCist'] = ( 1.0 - c/(b+c) if (b+c) > 0.0 else None,
                a/(a+b+c), 1.0 - c/(a+b+c)) if (a+b+c) > 0.0 else None

        if (flag & 128) != 0 and 'WCisct' in stats:
            a0 = self._sd.a0()
            b2 = self._sd.b2()
            b1 = self._sd.b1()
            c0 = self._sd.c0()
            stats['WCisct'] = ( 1.0 - c0/(b1+c0) if (b1+c0) > 0.0 else None,
                    (a0+b2) / (a0+b2+b1+c0), a0/(a0+b2+b1+c0),
                    1.0 - c0/(a0+b2+b1+c0)) if (a0+b2+b1+c0) > 0.0 else None

        if (flag & 256) != 0:
            if 'Dj' in stats: stats['Dj'] = self._sd.D()
            if 'Hst' in stats and (flag & 2048) != 0: stats['Hst'] = self._sd.Hst()
            if 'Gst' in stats and (flag & 4096) != 0: stats['Gst'] = self._sd.Gst()
            if 'Gste' in stats and (flag & 8192) != 0: stats['Gste'] = self._sd.Gste()

        if (flag & 1024) != 0:
            if 'MAF' in stats: stats['MAF'] = self._sd.MAF()
            if 'MAF_pop' in stats:
                stats['MAF_pop'] = [None if i == _eggwrapper.UNDEF else i for i in [self._sd.MAF_pop(i) for i in xrange(self._sd.k())]]

    def _get_as_stats(self, stats):
        if self._as.nsites() > 0:
            if 'numSp' in stats: stats['numSp'] = self._as.Sp()
            if 'numSpd' in stats and self._as.nsites_o() > 0: stats['numSpd'] = self._as.Spd()
            if 'numShP' in stats: stats['numShP'] = self._as.ShP()
            if 'numShA' in stats: stats['numShA'] = self._as.ShA()
            if 'numFxD' in stats: stats['numFxD'] = self._as.FxD()
            if 'numFxA' in stats: stats['numFxA'] = self._as.FxA()
