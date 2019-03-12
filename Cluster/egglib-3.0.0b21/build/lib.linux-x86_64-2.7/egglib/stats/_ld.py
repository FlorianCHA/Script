"""
This module contains linkage disequilibirum tools.
"""

__license__ = """
    Copyright 2015-2017 Stephane De Mita, Mathieu Siol

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
from . import _filter
from . import _structure
from . import _site

########################################################################

def _get_stats(pld, multiple_policy, min_freq):

    if pld.num_alleles1() > 2 or pld.num_alleles2() > 2:
        if multiple_policy == 'forbid':
            raise ValueError, 'at least one locus has more than two alleles and this is not allowed'
        elif multiple_policy == 'main':
            freq1 = [pld.freq1(i) for i in xrange(pld.num_alleles1())]
            freq2 = [pld.freq2(i) for i in xrange(pld.num_alleles2())]
            a1 = [freq1.index(max(freq1))]
            a2 = [freq2.index(max(freq2))]
        elif multiple_policy == 'average':
            a1 = [i for i in xrange(pld.num_alleles1()) if pld.freq1(i) >= min_freq]
            a2 = [i for i in xrange(pld.num_alleles2()) if pld.freq2(i) >= min_freq]
        else:
            raise ValueError, 'invalid value for `multiple_policy` option: {0}'.format(multiple_policy)
    else:
        a1 = [0]
        a2 = [0]

    n = 0
    D = 0.0
    Dp = 0.0
    r = 0.0
    rsq = 0.0
    for i in a1:
        for j in a2:
            n += 1
            pld.compute(i, j)
            D += pld.D()
            Dp += pld.Dp()
            r += pld.r()
            rsq += pld.rsq()
    if n == 0:
        return {'D': None, 'Dp': None, 'r': None, 'rsq': None}
    if n > 1:
        D /= n
        Dp /= n
        r /= n
        rsq /= n
    return {'D': D, 'Dp': Dp, 'r': r, 'rsq': rsq, 'n': n}

########################################################################

def pairwise_LD(locus1, locus2, struct=None, multiple_policy='main', min_freq=0):

    """
    This function computes linkage disequilibium between a pair of loci.
    
    :class:`.Site` instances must be used to describe loci. Only the order of
    samples is considered and both loci must have the same number of
    samples. If sites have a ploidy larger than one, genotypes are considered.

    :param locus1: A :class:`.Site` instance.

    :param locus2: A :class:`.Site` instance.

    :param struct: A :class:`.Structure` instance determining which samples
        to use. By default, use all samples.

    :param multiple_policy: Determine what is done if either input locus
        has more than two alleles. Possible values are, ``"forbid"``
        (raise an exception if thise occurs), ``"main"`` (take
        the most frequent allele of each locus) and ``"average"`` (
        compute the unweighted average over all possible pair of
        alleles). More options might be added in future versions. This
        option is ignored if both loci have less than three alleles. If
        ``"main"`` and there are several equally most frequent alleles,
        the first-occurring one is used (arbitrarily).

    :param min_freq: Only used if at least one site has more than two
        alleles and *multiple_policy* is set to *average*. Set the
        minimum absolute frequency to consider an allele.

    :return: A dictionary of linkage disequilibrium statistics. In case
        statistics cannot be computed (either site fixed, or less than
        two samples with non-missing data at both sites), computed
        values are replaced by ``None``. ``n`` gives the number of pairs
        of alleles considered.
    """

    if locus1.ns != locus2.ns: raise ValueError, 'the number of samples must match in the two loci'
    if locus1.ploidy != locus2.ploidy: raise ValueError, 'the ploidy must match in the two loci'
    if locus1.ns < 2: raise ValueError, 'the number of samples must be at least 2'
    if struct is not None:
        if struct.req_ns > locus1.ns: raise ValueError, 'invalid structure (largest sample index is out of range)'
        if struct.ploidy != 1: raise ValueError, 'invalid structure (ploidy must be 1)'
        struct = struct._obj
    ld = _eggwrapper.PairwiseLD()
    if not ld.process(locus1._obj, locus2._obj, 0, 1.0, struct):
        return {'D': None, 'Dp': None, 'r': None, 'rsq': None}
    return _get_stats(ld, multiple_policy, min_freq)

########################################################################

def matrix_LD(align, stats, multiple_policy='main', min_freq=0,
            min_n=2, max_maj=1.0, positions=None, filtr=None):

    """
    Compute the matrix of linkage disequilibrium statistics between all
    pairs of sites of the provided alignment. The computed statistics
    are selected by an argument of this function. Return a matrix (as a
    nested :class:`list`) of the requested statistics. In all cases, all
    pairs of sites are present in the returned matrices. If statistics
    cannot be computed, they are replaced by :data:`None`.

    The available statistics are:

    * ``d`` -- Distance between sites of the pairs.
    * ``D`` -- Standard linkage disequilibrium.
    * ``Dp`` -- Lewontin's D'.
    * ``r`` -- Correlation coefficient.
    * ``rsq`` -- Equivalent to r\ :sup:`2`.

    :param align: A :class:`.Align` instance.

    :param stats: Requested statistic or statistics (see list of
        available statistics above, as a single string or as a list of
        one or more of these statistics (in any order).

    :param multiple_policy: Specify what is done for pairs of sites for
        which at least one locus has only one allele. See
        :func:`pairwise_LD` for further description.

    :param min_freq: Only used if at least one site has more than two
        alleles and depending on the value of *multiple_policy*.  See
        :func:`pairwise_LD` for further description.

    :param min_n: Minimum number of samples used (this value must
        always be larger than 1). Sites not fulfilling this criterion
        will be dropped.

    :param max_maj: Maximum relative frequency of the majority allele.
        Sites not fulfilling this criterion will be dropped.

    :param positions: A sequence of positions, whose length must match
        the number of sites of the provided alignment. Used in the
        return value to describe the used sites, and, if requested, to
        compute the distance between sites. By default, the position of
        sites in the original alignment is used.

    :param filtr: A :class:`.Filter` instance providing the list
        of valid (including missing) allelic values.

    :return: Returns a tuple with two items: first is the list of
        positions of sites used in the matrix (a subset of the sites of
        the provided alignment), with positions provided by the
        corresponding argument (by default, the index of sites); second
        is the matrix, as the nested lower half matrix. The matrix
        contains items for all i and j indexes with 0 <= j <= i < n
        where n is the number of retained sites. The content of the
        matrix is represented by a single value (if a single statistic
        has been requested) or as a list of 1 or more values (if a list
        of 1 or more, accordingly, statistics have been requested), or
        ``None`` for the diagonal or if the pairwise comparison was
        dropped for any reason.
    """

    # initialize local variables
    n = align.ls
    mLD = _eggwrapper.MatrixLD()
    retained = []

    # check arguments
    min_n = int(min_n)
    if min_n < 2: raise ValueError, 'too small value for `min_n` argument: {0}'.format(min_n)
    max_maj = float(max_maj)
    if max_maj < 0.5 or max_maj > 1.0: raise ValueError, 'invalid value for `max_maj` argument: {0}'.format(max_maj)

    if filtr == None: filtr = _filter.Filter()
    elif not isinstance(filtr, _filter.Filter): raise TypeError, 'invalid type provided for `filtr` argument'

    if positions == None:
        positions = map(float, range(n))
    elif len(positions) != n:
        raise ValueError, 'provided list of positions does not have the right number of items'

    if stats in ['D', 'Dp', 'r', 'rsq']:
        multi = False
        stats = [stats]
    else:
        multi = True

    # make a list of variable sites
    sites = []
    final_positions = []
    sd = _eggwrapper.SiteDiversity()
    frq = _eggwrapper.FreqHolder()
    frq.setup_raw(1, 1, align.no, 1)
    frq.setup_pop(0, 0, 0, align.ns)
    for i in xrange(align.ls):
        site = _site.site_from_align(align, i, filtr)
        frq.process_site(site._obj)
        if (sd.process(frq)&2) != 0:
            if sd.Aing() > 1.0:
                sites.append(site)
                final_positions.append(positions[i])
                mLD.load(site._obj, positions[i])

    matrix = [[None for j in xrange(i+1)] for i in xrange(len(sites))]

    # compute LD
    mLD.computeLD(min_n, max_maj)

    # extract requested values
    for i in xrange(mLD.num_pairs()):
        pld = mLD.pairLD(i)
        idx2 = mLD.index1(i)
        idx1 = mLD.index2(i) # reverse indexes
        computed_stats = _get_stats(pld, multiple_policy, min_freq)
        matrix[idx1][idx2] = []
        for stat in stats:
            if stat == 'd': matrix[idx1][idx2].append( final_positions[idx1] - final_positions[idx2] )
            elif stat == 'D': matrix[idx1][idx2].append( computed_stats['D'] )
            elif stat == 'Dp': matrix[idx1][idx2].append( computed_stats['Dp'] )
            elif stat == 'r': matrix[idx1][idx2].append( computed_stats['r'] )
            elif stat == 'rsq': matrix[idx1][idx2].append( computed_stats['rsq'] )
            else: raise ValueError, 'invalid statistic code: `{0}`'.format(stat)
        if not multi:
            matrix[idx1][idx2] = matrix[idx1][idx2][0]

    # return
    return final_positions, matrix
