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

from .. import _eggwrapper

########################################################################

def get_structure(data, lvl_clust=None, lvl_pop=None, lvl_indiv=None,
                    pop_filter=None, ploidy=None, skip_outgroup=False):

    """
    Create a new :class:`.Structure` instance based on the group labels
    of a :class:`Align` or :class:`Container` instance.

    :param data: an :class:`.Align` or :class:`.Container` instance
        containing the grouping levels to be processed.

    :param lvl_clust: index of the grouping level containing cluster
        labels. If ``None``, all populations are placed in a single
        cluster with label 0.

    :param lvl_pop: index of the grouping level containing population
        labels. If ``None``, all individuals of a cluster are placed in
        a single cluster with the same label as their cluster.

    :param lvl_indiv: index of the grouping level containing individual
        labels. If ``None``, individuals are not specified and each
        sample is placed in a haploid individual, for both the ingroup
        and the outgroup (unless *skip_outgroup* is ``True``), meaning
        that outgroup group labels are ignored.

    :param pop_filter: process only the populations bearing the label or
        labels provided in the provided list. If ``None``, all
        populations are processed. An empty list is the same as
        ``None``. It is allowed to include repeated labels in the list,
        as well as labels that are not actually present in the data.

    :param ploidy: indicate the ploidy. Must be a positive number and,
        if specified, data must match the value. If not specified,
        ploidy will be detected automatically (it must still be
        consistent over all ingroup and outgroup individuals). Ploidy is
        ignored if *lvl_indiv* is ``None``.

    :param skip_outgroup: specify if outgroup samples should be skipped.
        No effect if there are no outgroup samples. If *lvl_indiv* is
        not ``None``, the group label of outgroup samples are considered
        to be individual labels and the outgroup individuals will be
        recorded (requiring a consistent ploidy as well). If *lvl_indiv*
        is ``None``, outgroup samples are imported in one-sample
        individuals just like ingroup samples. If the *skip_outgroup*
        flag is set, outgroup samples are not imported at all.

    :return: A new :class:`.Structure`.
    """

    obj = Structure()
    obj.get_structure(data, lvl_clust=lvl_clust, lvl_pop=lvl_pop,
            lvl_indiv=lvl_indiv, pop_filter=pop_filter, ploidy=ploidy,
            skip_outgroup=skip_outgroup)
    return obj

########################################################################

def make_structure(ingroup, outgroup):

    """
    Create a new :class:`.Structure` instance based based on the
    structure provided as dictionaries. The two arguments must match the
    format of the return value of the :meth:`.Structure.as_dict` method
    (see the documentation). Either argument can be replaced by ``None``
    which is equivalent to an empty dictionary (no samples). Not that
    all keys must be positive integers.

    :param ingroup: a three-fold nested dictionary of ingroup samples
        indexes, or ``None``.

    :param outgroup: a dictionary of outgroup samples indexes, or
        ``None``.

    :return: A new :class:`.Structure`.
    """

    obj = Structure()
    obj.make_structure(ingroup, outgroup)
    return obj

########################################################################

class Structure(object):

    """
    Describe the organisation of samples in individuals, populations,
    and clusters of populations. The structure is necessarily
    hierarchical and all levels are always defined but it is possible to
    bypass any level if information for this level is not available or
    irrelevant. The number of individuals per population can vary, but
    the number of samples per individuals (that is, the ploidy) must be
    constant.

    New objects are created using the fonctions :func:`.get_structure`
    and :func:`.make_structure`, and old objects can be recycled with
    the corresponding methods :meth:`.get_structure` and
    :meth:`make_structure`.
    """

    ####################################################################

    def __init__(self):
        self._obj = _eggwrapper.StructureHolder()

    ####################################################################

    def get_structure(self, data, lvl_clust=None, lvl_pop=None,
                lvl_indiv=None, pop_filter=None, ploidy=None, skip_outgroup=False):

        """
        Reset the instance as if it was built using :func:`.get_structure`.
        The definitions of arguments are identical.
        """

        # must be done before loading filter
        self._obj.reset()

        # convert arguments to C++ friendly + checking
        ng = data.ng

        if lvl_clust == None: lvl_clust = _eggwrapper.MISSING
        elif lvl_clust >= ng: raise ValueError, 'invalid value for `lvl_clust`: invalid index'

        if lvl_pop == None: lvl_pop = _eggwrapper.MISSING
        elif lvl_pop >= ng: raise ValueError, 'invalid value for `lvl_pop`: invalid index'

        if lvl_indiv == None: lvl_indiv = _eggwrapper.MISSING
        elif lvl_indiv >= ng: raise ValueError, 'invalid value for `lvl_indiv`: invalid index'

        if ploidy == None: ploidy = _eggwrapper.MISSING
        elif ploidy < 0: raise ValueError, 'ploidy must be strictly positive'

        if pop_filter is not None:
            self._obj.reserve_filter(len(pop_filter))
            for lbl in pop_filter: self._obj.add_pop_filter(lbl)

        # import data
        self._obj.get_structure(data._obj, lvl_clust, lvl_pop, lvl_indiv, ploidy, skip_outgroup)
        return

    ####################################################################

    def make_structure(self, ingroup, outgroup):

        """
        Reset the instance as if it was built using :func:`.make_structure`.
        The definitions of arguments are identical.
        """

        self._obj.reset()
        pop_cache = set()
        idv_cache = set()
        idx_cache = set()

        if ingroup is not None:
            for clt, clt_d in ingroup.iteritems():
                if not isinstance(clt, int) or clt<0: raise TypeError, 'cluster label must be an integer >=0'
                clt_o = self._obj.add_cluster(clt)

                for pop, pop_d in clt_d.iteritems():
                    if not isinstance(pop, int) or pop<0: raise TypeError, 'population label must be an integer >=0'
                    if pop in pop_cache: raise ValueError, 'structure not hierarchical, duplicated population label: {0}'.format(pop)
                    pop_cache.add(pop)
                    pop_o = self._obj.add_population(pop, clt_o)

                    for idv, samples in pop_d.iteritems():
                        if not isinstance(idv, int) or idv<0: raise TypeError, 'individual label must be an integer >=0'
                        if idv in idv_cache: raise ValueError, 'structure not hierarchical, duplicated individual label: {0}'.format(idv)
                        idv_cache.add(idv)
                        idv_o = self._obj.add_individual_ingroup(idv, clt_o, pop_o)

                        for idx in samples:
                            if idx in idx_cache: raise ValueError, 'ingroup index found several times: {0}'.format(idx)
                            idx_cache.add(idx)
                            self._obj.add_sample_ingroup(idx, clt_o, pop_o, idv_o)

        idx_cache = set()

        if outgroup is not None:
            for idv, samples in outgroup.iteritems():
                if not isinstance(idv, int) or idv<0: raise TypeError, 'individual label must be an integer >=0'
                if idv in idv_cache: raise ValueError, 'structure not hierarchical, duplicated individual label: {0}'.format(idv)
                idv_cache.add(idv)
                idv_o = self._obj.add_individual_outgroup(idv)

                for idx in samples:
                    if idx in idx_cache: raise ValueError, 'outgroup index found several times: {0}'.format(idx)
                    idx_cache.add(idx)
                    self._obj.add_sample_outgroup(idx, idv_o)

        self._obj.check_ploidy()

    ####################################################################

    def as_dict(self):

        """
        Return a tuple of two :class:`dict` representing, respectively,
        the ingroup and outgroup structure.

        The ingroup dictionary is a three-fold nested dictionary
        (meaning it is a dictionary of dictionaries of dictionaries)
        holding lists of sample indexes. The keys are, respectively,
        cluster, population, and individual labels. Based on how the
        instance was created, there may be just one item or even none at
        all in any dictionary. In practice, if ``d` is the ingroup
        :class:`dict` and ``clt``, ``pop`` and ``idv`` are, respectively,
        cluster, population, and individual labels, the expression
        ``d[clt][pop][idv]`` will yield a list of sample indexes.

        The outgroup dictionary is a non-nested dictionary with
        individual labels as keys and lists of sample indexes as values.
        """

        ingroup = {}
        for idx_clust in xrange(self._obj.num_clust()):
            clust = self._obj.get_cluster(idx_clust)
            ingroup[clust.get_label()] = {}
            for idx_pop in xrange(clust.num_pop()):
                pop = clust.get_population(idx_pop)
                ingroup[clust.get_label()][pop.get_label()] = {}
                for idx_indiv in xrange(pop.num_indiv()):
                    indiv = pop.get_indiv(idx_indiv)
                    ingroup[clust.get_label()][pop.get_label()][indiv.get_label()] = [
                        indiv.get_sample(idx_sam) for idx_sam in xrange(indiv.num_samples())]
        outgroup = {}
        for idx_indiv in xrange(self._obj.num_indiv_outgroup()):
            indiv = self._obj.get_indiv_outgroup(idx_indiv)
            outgroup[indiv.get_label()] = [
                indiv.get_sample(idx_sam) for idx_sam in xrange(indiv.num_samples())]
        return ingroup, outgroup

    ####################################################################

    def make_auxiliary(self):

        """
        Return a new :class:`.Structure` instance describing the organisation
        of individuals in clusters and populations (ignoring the intra-individual
        level) using the rank of individuals as indexes, the individuals
        being ranked in the order of increasing cluster and population
        indexes.
        """

        ret = Structure()
        self._make_auxiliary(self._obj, ret._obj)
        return ret

    @staticmethod
    def _make_auxiliary(source, dest):
        dest.reset()
        cur = 0
        for c in xrange(source.num_clust()):
            clu = source.get_cluster(c)
            for p in xrange(clu.num_pop()):
                pop = clu.get_population(p)
                for i in xrange(pop.num_indiv()):
                    dest.process_ingroup(cur, clu.get_label(), pop.get_label(), pop.get_indiv(i).get_label())
                    cur += 1
        cur = 0
        for i in xrange(source.num_indiv_outgroup()):
            dest.process_outgroup(cur, source.get_indiv_outgroup(i).get_label())
            cur += 1
        dest.check_ploidy(1)

    ####################################################################

    @property
    def ns(self):

        """
        Number of processed ingroup samples.
        """

        return self._obj.get_ns()

    ####################################################################

    @property
    def no(self):

        """
        Number of processed outgroup samples.
        """

        return self._obj.get_no()

    ####################################################################

    @property
    def req_ns(self):

        """
        Required number of ingroup sample index in objects using this
        structure (equal to the largest index overall plus one).
        """

        return self._obj.get_ns_req()

    ####################################################################

    @property
    def req_no(self):

        """
        Required number of ougroup sample index in objects using this
        structure (equal to the largest index overall plus one).
        """

        return self._obj.get_no_req()

    ####################################################################

    @property
    def num_clust(self):

        """
        Number of clusters.
        """

        return self._obj.num_clust()

    ####################################################################

    @property
    def num_pop(self):

        """
        Total number of populations.
        """

        return self._obj.num_pop()

    ####################################################################

    @property
    def num_indiv_ingroup(self):

        """
        Total number of ingroup individuals.
        """

        return self._obj.num_indiv_ingroup()

    ####################################################################

    @property
    def num_indiv_outgroup(self):

        """
        Number of outgroup individuals.
        """

        return self._obj.num_indiv_outgroup()

    ####################################################################

    @property
    def ploidy(self):

        """
        Ploidy.
        """

        return self._obj.get_ploidy()

    ####################################################################

    def get_samples(self):
        
        """
        Return a new :class:`list` with the index of all samples from
        the ingroup (all clusters, populations, and individuals stored
        in this object, excluding the outgroup). 
        """

        return [self._obj.get_indiv_ingroup(i). get_sample(j)
                    for i in xrange(self._obj.num_indiv_ingroup())
                            for j in xrange(self._obj.get_ploidy())]
