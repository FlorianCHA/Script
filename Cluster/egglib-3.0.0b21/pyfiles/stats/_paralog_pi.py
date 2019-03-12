__license__ = """
    Copyright 2016 Stephane De Mita, Mathieu Siol

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

def paralog_pi(align, struct_p, struct_i, max_missing=0, filtr=None):

    """
    Compute diversity statistics for a gene family (Innan 2003). An
    estimate of genetic diversity is provided for every paralog and for
    every pair of paralogs, provided that enough non-missing data is
    available (at least 2 samples are required). Note that sites with
    more than two alleles are always considered.

    :param align: an :class:`.Align` containing the sequence of all
        available paralog for all samples. The outgroup is ignored.

    :param struct_p: a :class:`.Structure` providing the organisation of
        sequences in paralogs. This structure must have a ploidy of 1 (no
        individual structure). Clusters, if defined, are ignored. The
        sequences of all individuals for a given paralog should be grouped
        together in that structure. There might be a different number of
        sequences per paralog due to missing data.

    :param struct_i: a :class:`.Structure` providing the organisation of
        sequences in individuals. This structure must have a ploidy of 1
        (no individual structure). Clusters, if defined, are ignored. The
        sequences of all paralogs for a given individual should be grouped
        together in that structure. There might be a different number of
        sequences per individual due to missing data.

    :param max_missing: maximum proportion of missing data (if there are
        more missing data at a site, the site is ignored).

    :param filtr: a :class:`.Filter` instance providing the list of
        exploitable allelic values and those that should be considered
        as missing. By default, use :data:`.filter_dna`.

    :return: A new :class:`.ParalogPi` instance which provides methods
        to access the number of used sites and the diversity for each
        paralog/paralog pair.
    """

    pp = ParalogPi()
    pp.setup(struct_p, struct_i)
    pp.process_align(align, max_missing, filtr)
    return pp

########################################################################

class ParalogPi(object):

    """
    Class computing Innan's within- and between-paralog diversity
    statistics. See :func:`.paralog_pi` for more details. This class
    can be used directly (1) to analyse data with more efficiency (by
    reusing the same instance) or (2) to combine data from different
    alignments, or (3) for pass individual sites. Do first call
    :meth:`.setup`.
    """

    def __init__(self):
        self._obj = _eggwrapper.ParalogPi()
        self._req = 0
        self._iter = _site.Iterator()

    ##############

    def setup(self, struct_p, struct_i):

        """
        Specify the structure in paralog and individuals. The two arguments
        are :class:`.Structure` instances as described for :func:`.paralog_pi`.
        Only this method resets the instance.
        """

        if struct_p.ploidy != 1: raise ValueError, 'ploidy is required to be 1'
        if struct_i.ploidy != 1: raise ValueError, 'ploidy is required to be 1'
        self._req = struct_p.req_ns
        self._obj.reset(struct_p._obj, struct_i._obj)

    ##############

    def process_align(self, aln, max_missing=0, filtr=None):

        """
        Process an alignment matching the structure passed to :meth:`.setup`.
        Diversity estimates are incremented (no reset).

        :param aln: an :class:`.Align` instance.
        :param max_missing: maximum proportion of missing data.
        :param filtr: a class:`.Filter` instance.
        """

        if not isinstance(aln, _interface.Align): raise TypeError, 'expect an Align instance'
        if aln.ns < 2: return
        if aln.ns > self._req: raise ValueError, 'unsufficient number of samples in alignment'
        if filtr is None: filtr = _filter.filter_dna
        for site in self._iter.iter(aln, None, filtr, max_missing, False):
            self._obj.load(site)

    ##############

    def process_site(self, site):

        """
        Process a site matching the structure passed to :meth:`.setup`.
        Diversity estimates are incremented (no reset).

        :param site: a :class:`.Site` instance.
        """

        if site.ns < 2: return
        if site.ns > self._req: raise ValueError, 'unsufficient number of samples in site'
        if filtr is None: filtr = _filter.filter_dna
        self._obj.load(site._obj)

    ##############

    def num_sites(self, *args):

        """
        num_sites([i[, j]])

        Number of sites with any data (without arguments), with data for
        paralog *i* (if only *i* specified), or with data for the pair
        of paralogs *i* and *j* (if both specified).
        """

        if len(args) == 0:
            return self._obj.num_sites_tot()
        elif len(args) == 1:
            if args[0] < 0 or args[0] >= self._obj.num_paralogs(): raise IndexError, 'invalid paralog index'
            return self._obj.num_sites_paralog(args[0])
        elif len(args) == 2:
            if args[0] < 0 or args[0] >= self._obj.num_paralogs(): raise IndexError, 'invalid paralog index'
            if args[1] < 0 or args[1] >= self._obj.num_paralogs() or args[1] == args[0]: raise IndexError, 'invalid paralog index'
            return self._obj.num_sites_pair(args[0], args[1])
        else:
            raise ValueError, 'invalid number of arguments'

    ##############

    def Piw(self, i):

        """
        Get within-paralog diversity for paralog *i*.
        """

        if i<0 or i>=self._obj.num_paralogs(): raise IndexError, 'invalid paralog index'
        if self._obj.num_sites_paralog(i) < 1: return None
        return self._obj.Piw(i)

    ##############

    def Pib(self, i, j):

        """
        Get between-paralog diversity for paralogs *i* and *j*.
        """

        if i<0 or i>=self._obj.num_paralogs(): raise IndexError, 'invalid paralog index'
        if j<0 or j>=self._obj.num_paralogs() or j==i: raise IndexError, 'invalid paralog index'
        if self._obj.num_sites_pair(i, j) < 1: return None
        return self._obj.Pib(i, j)
