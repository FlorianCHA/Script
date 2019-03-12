"""
Contains the :function:`.haplotypes` function.
"""

__license__ = """
    Copyright 2015 Stephane De Mita, Mathieu Siol

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

from .. import _eggwrapper, _interface
from . import _site, _filter

########################################################################

def haplotypes(sites, impute_threshold=0, struct=None, filtr=None,
            max_missing=0.0, consider_outgroup_missing=False, dest=None,
            multiple=False):

    """
    Identify haplotypes from sites provided as either an :class:`.Align`
    instance or a :class:`list` of :class:`.Site` instances, and return data as a
    single :class:`.Site` instance containing one sample for each sample
    of the original data. Alleles in the returned site are representing
    all identified haplotypes (or missing data when the haplotypes could
    not be derived.

    .. note::
        There must be at least one site with at least two alleles (overall,
        including the outgroup), otherwise the produced site only contains
        missing data.

    :param sites: an :class:`.Align` instance, or a :class:`list` of
        :class:`.Site` instances.

    :param impute_threshold: by default, all samples with a least one
        occurrence of missing data will be treated as missing data. If
        this argument is more than 0, the provided value will be used as
        maximum number of missing data. All samples with as many or less
        missing data will be processed to determine which extant of
        haplotype they might belong (to which they are identical save
        for missing data). If there is only one such haplotype, the
        corresponding samples will be treated as a repetition of this
        haplotype. This option will never allow detecting new
        haplotypes. Only small values of this option make sense.

    :param struct: a :class:`.Structure` instance defining the samples
        to process. Only valid if *sites* is a :class:`.Align`. The
        population and cluster structures are not used. If the ploidy is
        larger than 1, the individuals are used, and sites are assumed
        to be phased.

    :param filtr: a :class:`.Filter` instance controlling what allelic
        values are supported. By default, assume DNA sequences. Only
        allowed if an :class:`.Align` is passed.

    :param max_missing: maximum proportion of missing data to process a
        site. Only considered if an :class:`.Align` is passed.

    :param consider_outgroup_missing: if ``True``, outgroup samples are
        included in the count for missing data (by default, outgroup
        samples are not considered). Only considered if an :class:`.Align`
        is passed.

    :param dest: a :class:`.Site` instance that will be reset and in
        which data will be placed. If specified, this function returns
        nothing. By default, the function returns a new :class:`.Site`
        instance.

    :param multiple: allow sites with more than two alleles in the
        ingroup.

    :return: A :class:`.Site` instance, (if *dest* is ``None``) or ``None``
        (otherwise).
    """

    # check arguments
    if filtr is None: filtr = _filter.filter_dna
    elif not isinstance(sites, _interface.Align): raise ValueError, 'cannot pass a value for `filtr` argument if `sites` is not an Align'
    if max_missing < 0.0 or max_missing > 1.0: raise ValueError, 'max_missing argument out of range'
    if impute_threshold < 0: raise ValueError, 'invalid value for impute_threshold argument'

    obj = _eggwrapper.Haplotypes()
    frq = _eggwrapper.FreqHolder()
    impute_idx = []

    # pass sites if they come from an alignement
    if isinstance(sites, _interface.Align):
        if struct is not None:
            struct = struct._obj
            pl = struct.get_ploidy()
            ni = struct.num_indiv_ingroup()
            no = struct.num_indiv_outgroup()
            if struct.get_ns_req() > sites._obj.get_nsam_i(): raise ValueError, 'structure object does not match alignment'
            if struct.get_no_req() > sites._obj.get_nsam_o(): raise ValueError, 'structure object does not match alignment'
        else:
            pl = 1
            ni = sites._obj.get_nsam_i()
            no = sites._obj.get_nsam_o()

        obj.setup(struct)
        frq.setup_raw(1, 1, no, pl)
        frq.setup_pop(0, 0, 0, ni)
        site = _eggwrapper.SiteHolder(pl)

        if consider_outgroup_missing: max_missing = int(max_missing * pl * (ni+no))
        else: max_missing = int(max_missing * pl * ni)

        for i in xrange(sites._obj.get_nsit()):
            site.reset(pl)
            ret = site.process_align(sites._obj, i, struct, filtr._obj, max_missing, consider_outgroup_missing)
            if not ret: continue
            frq.process_site(site)
            nall = frq.frq_ingroup().num_alleles_eff()
            if site.get_nall() > 1 and (multiple == True or nall < 3):
                obj.load(site)
                if impute_threshold > 0: impute_idx.append(i)

        if obj.n_sites() == 0:
            if dest is None:
                return _site.site_from_list([[None]*pl for i in xrange(ni)], [[None]*pl for i in xrange(no)])
            else:
                dest.process_list([[None]*pl for i in xrange(ni)], [[None]*pl for i in xrange(no)], reset=True)
                return None

    # process list of Site
    else:
        if struct is not None: raise ValueError, 'cannot pass a value for `struct` argument if `sites` is not an Align'
        if len(sites) > 0:
            pl = sites[0]._obj.get_ploidy()
            ni = sites[0]._obj.get_ning()
            no = sites[0]._obj.get_nout()
        obj.setup(None)

        frq.setup_raw(1, 1, no, pl)
        frq.setup_pop(0, 0, 0, ni)

        if consider_outgroup_missing: max_missing = int(max_missing * pl * (ni+no))
        else: max_missing = int(max_missing * pl * ni)

        for i, site in enumerate(sites):
            if pl != site._obj.get_ploidy(): raise ValueError, 'ploidy must be consistent between sites'
            if ni != site._obj.get_ning(): raise ValueError, 'number of samples must be consistent between sites'
            if no != site._obj.get_nout(): raise ValueError, 'number of outgroup samples must be consistent between sites'
            frq.process_site(site._obj)
            nall = frq.frq_ingroup().num_alleles_eff()
            if site._obj.get_nall() > 1 and (multiple == True or nall < 3):
                obj.load(site._obj)
                if impute_threshold > 0: impute_idx.append(i)

    # complete haplotype detection
    obj.cp_haplotypes()

    # impute if requested
    if impute_threshold > 0:
        obj.prepare_impute(impute_threshold)
        for i in impute_idx:
            if isinstance(sites, _interface.Align):
                site.reset(pl)
                site.process_align(sites._obj, i, struct, filtr._obj, _eggwrapper.MAX, True)
                obj.resolve(site)
            else:
                obj.resolve(sites[i]._obj)
        obj.impute()

    # generate site
    if dest is None: dest = _site.Site()
    obj.get_site(dest._obj)
    return dest
