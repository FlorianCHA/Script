"""
Main manager of the coalesce simulator module.
"""

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

import sys
from .. import _eggwrapper, _interface, _tree
from ..stats import _filter
from . import _param_helpers

########################################################################

class Simulator(object):

    """
    Manager of the coalescent simulator. The constructor takes arguments
    controlling the demographic and mutation models used for simulation.
    Only the number of populations is required at the time of
    construction. Once it is set, it can be never modified. Other
    constructor arguments are all optional and can be set or modified
    later using the :attr:`.params` instance attribute (either using its
    :meth:`~.ParamDict.update` method or the ``[]`` operator, such as in
    ``simulator.params['theta'] = 2.85``. Keyword arguments are passed as is to
    :meth:`~.ParamDict.update`. List-based parameters can be set if all
    values are provided in a sequence. Matrix-based parameters cannot be
    set here.

    :param num_pop: number of populations.
    :param migr: migration rate.

    Other keyword arguments can be any of the parameters defined in the
    table below.

    +-----------------+------------------------------------+-----------------+
    |Parameter        | Definition                         | Default         |
    +=================+====================================+=================+
    |``num_pop``      | Number of populations              | None, required  |
    +-----------------+------------------------------------+-----------------+
    |``num_sites``    | Number of sites                    | 0 (ISM)         |
    +-----------------+------------------------------------+-----------------+
    |``num_mut``      | Fixed number of mutations          | 0               |
    +-----------------+------------------------------------+-----------------+
    |``theta``        | :math:`4N_0\mu` parameter          | 0.0             |
    +-----------------+------------------------------------+-----------------+
    |``recomb``       | :math:`4N_0c` parameter            | 0.0             |
    +-----------------+------------------------------------+-----------------+
    |``mut_model``    | Mutation model (among ``KAM``,     | ``KAM``         |
    |                 | ``IAM``, ``SMM`` and ``TPM``)      |                 |
    +-----------------+------------------------------------+-----------------+
    |``TPM_proba``    | Probability parameter of TPM       | 0.5             |
    +-----------------+------------------------------------+-----------------+
    |``TPM_param``    | Shape parameter of TPM             | 0.5             |
    +-----------------+------------------------------------+-----------------+
    |``num_alleles``  | Number of alleles for KAM          | 2               |
    +-----------------+------------------------------------+-----------------+
    |``rand_start``   | Pick start allele randomly         | False           |
    |                 | for KAM (boolean)                  |                 |
    +-----------------+------------------------------------+-----------------+
    |``num_chrom``    | Number of sampled chromosomes      | 0 for all       |
    |                 | (per population)                   |                 |
    +-----------------+------------------------------------+-----------------+
    |``num_indiv``    | Number of sampled individuals      | 0 for all       |
    |                 | (per population)                   |                 |
    +-----------------+------------------------------------+-----------------+
    |``N``            | Population size, expressed         | 1.0 for all     |
    |                 | relatively to :math:`N_0`          |                 |
    |                 | (per population)                   |                 |
    +-----------------+------------------------------------+-----------------+
    |``G``            | Exponential growth/decline rate,   | 0.0 for all     |
    |                 | negative values mean decline       |                 |
    |                 | (per population)                   |                 |
    +-----------------+------------------------------------+-----------------+
    |``s``            | Population selfing probability     | 0.0 for all     |
    |                 | (per population)                   |                 |
    +-----------------+------------------------------------+-----------------+
    |``outgroup``     | Label of samples to place in the   | None            |
    |                 | outgroup (population index for     |                 |
    |                 | normal samples, label for delayed) |                 |
    +-----------------+------------------------------------+-----------------+
    |``site_pos``     | Site position, as values           | Equally spread  |
    |                 | between 0.0 and 1.0                |                 |
    |                 | (per site)                         |                 |
    +-----------------+------------------------------------+-----------------+
    |``site_weight``  | Site mutation weight, controlling  | 1.0 for all     |
    |                 | the relative probability of sites  |                 |
    |                 | (per site)                         |                 |
    +-----------------+------------------------------------+-----------------+
    |``migr_matrix``  | Pairwise migration rate matrix     | 0.0 for all     |
    |                 | (the diagonal cannot be set)       |                 |
    +-----------------+------------------------------------+-----------------+
    |``trans_matrix`` | Matrix of transition weights       | 1.0 for all     |
    |                 | between pairs of alleles           |                 |
    +-----------------+------------------------------------+-----------------+
    |``events``       | List of events added by the user   | Empty           |
    +-----------------+------------------------------------+-----------------+
    |``seed``         | Random number generator seed       | Based on clock  |
    +-----------------+------------------------------------+-----------------+
    |``max_iter``     | Maximum number of iterations       | 100,000         |
    +-----------------+------------------------------------+-----------------+

    The following table presents the categories of events that can be
    added to the ``events`` list using :meth:`~.EventList.add`.

    +------------------+-------------------------+----------------------------------+
    |Event code        | Description             | Parameters                       |
    +==================+=========================+==================================+
    |``size``          | Change population size  | ``T`` -- date (1)                |
    |                  |                         +----------------------------------+
    |                  |                         | ``N`` -- new size                |
    |                  |                         +----------------------------------+
    |                  |                         | ``idx`` -- population index (2)  |
    +------------------+-------------------------+----------------------------------+
    |``migr``          | Change all migration    | ``T`` -- date (1)                |
    |                  | rate                    +----------------------------------+
    |                  |                         | ``M`` -- migration rate (all     |
    |                  |                         | pairwise migration rates are     |
    |                  |                         | set to ``M/(num_pop-1)``)        |
    +------------------+-------------------------+----------------------------------+
    |``pair_migr``     | Change pairwise         | ``T`` -- date (1)                |
    |                  | migration rate          +----------------------------------+
    |                  |                         | ``M`` -- pairwise migration      |
    |                  |                         | rate                             |
    |                  |                         +----------------------------------+
    |                  |                         | ``src`` -- source population     |
    |                  |                         | index                            |
    |                  |                         +----------------------------------+
    |                  |                         | ``idx`` -- destination           |
    |                  |                         | population index                 |
    +------------------+-------------------------+----------------------------------+
    |``growth``        | Change population       | ``T`` -- date (1)                |
    |                  | exponential             +----------------------------------+
    |                  | growth/decline rate     | ``G`` -- new rate                |
    |                  |                         +----------------------------------+
    |                  |                         | ``idx`` -- population index (2)  |
    +------------------+-------------------------+----------------------------------+
    |``selfing``       | Change population       | ``T`` -- date (1)                |
    |                  | self-fertilization      +----------------------------------+
    |                  | rate                    | ``s`` -- new rate                |
    |                  |                         +----------------------------------+
    |                  |                         | ``idx`` -- population index (2)  |
    +------------------+-------------------------+----------------------------------+
    |``recombination`` | Change recombination    | ``T`` -- date (1)                |
    |                  | rate                    +----------------------------------+
    |                  |                         | ``R`` -- new rate                |
    +------------------+-------------------------+----------------------------------+
    |``bottleneck``    | Apply a bottleneck      | ``T`` -- date (1)                |
    |                  |                         +----------------------------------+
    |                  |                         | ``S`` -- bottleneck strength (3) |
    |                  |                         +----------------------------------+
    |                  |                         | ``idx`` -- population index (2)  |
    +------------------+-------------------------+----------------------------------+
    |``admixture``     | Move lineages from one  | ``T`` -- date (1)                |
    |                  | population to another   +----------------------------------+
    |                  |                         | ``proba`` -- migration           |
    |                  |                         | probability (in [0,1] range      |
    |                  |                         +----------------------------------+
    |                  |                         | ``src`` -- source population     |
    |                  |                         | index                            |
    |                  |                         +----------------------------------+
    |                  |                         | ``dst`` -- destination           |
    |                  |                         | population index                 |
    +------------------+-------------------------+----------------------------------+
    |``merge``         | Merge a population to   | ``T`` -- date (1)                |
    |                  | another (take all       +----------------------------------+
    |                  | lineages from ``src``,  | ``src`` -- source population     |
    |                  | move them to ``dst``,   | index                            |
    |                  | and remove ``src``      +----------------------------------+
    |                  |                         | ``dst`` -- destination           |
    |                  |                         | population index                 |
    +------------------+-------------------------+----------------------------------+
    |``sample``        | Perform a delayed       | ``T`` -- date (1)                |
    |                  | a some point in the     +----------------------------------+
    |                  | past in one of the      | ``idx`` -- population index      |
    |                  | populations             |                                  |
    |                  |                         +----------------------------------+
    |                  |                         | ``label`` -- group label (4)     |
    |                  |                         +----------------------------------+
    |                  |                         | ``num_chrom`` -- number of       |
    |                  |                         | sampled chromosomes              |
    |                  |                         +----------------------------------+
    |                  |                         | ``num_indiv`` -- number of       |
    |                  |                         | sampled individuals              |
    +------------------+-------------------------+----------------------------------+

    1. Time is expressed in units of :math:`4N_0` generations.
    2. If ``idx`` is omitted, the event is applied to all populations 
       at once.
    3. Bottleneck strength is expressed in time units. A bottleneck is
       implemented as a period of time during which coalescences are the
       only event allowed to occur.
    4. The label of delayed sample can be set to the same value than the
       populations index (in this case, delayed samples will have the
       same label than normal samples from the same population), or to a
       different label, as the user's option. Based on the value of the
       *outgroup* option, this label controls whether delayed samples
       should go to the outgroup.

    The parameters ``num_chrom``, ``num_indiv``, ``N``, ``G``, ``s``,
    ``site_pos``, and ``site_weight`` are represented by
    :class:`.ParamList` instances that behave like lists (except that
    their length cannot be changed). In particular, :class:`.ParamList`
    support subscript indexing and it can also be initialized by passing
    a sequence. The parameters ``migr_matrix`` and ``trans_matrix`` are
    represented by :class:`.ParamMatrix` instances that support a
    double-index subscript system to read/changes values (as in
    ``params['migr_matrix'][i,j]`` to access the value at row *i* and
    column *j*. Diagonal values are read as ``None`` and cannot be
    changed. Finally, ``events`` is represented by a :class:`.EventList`
    instance that exhibit limited list functionality and provides
    methods to add, read, and modify events. See this class for more
    information about out editing the list of events.
    """

    ####################################################################

    def __init__(self, num_pop, migr=0.0, **kwargs):
        self._params = _param_helpers.ParamDict(num_pop)
        self._params.set_migr(migr)
        self._params.update(**kwargs)
        self._coalesce = _eggwrapper.Coalesce()
        self._align = _interface.Align._create_from_data_holder(self._coalesce.data())

    ####################################################################

    @property
    def params(self):

        """
        Simulation parameters of this instance. This object is a
        instance of :class:`.ParamDict`, which is a specialized version
        of :class:`dict` that does not let the users add or remove
        parameters, but lets them modify values of parameters
        (similarly, the number of items of per-population or per-site
        parameters cannot be modified).
        """

        return self._params

    ####################################################################

    def simul(self):

        """
        Perform a single simulation based on the current value of
        parameters stored in :attr:`.params`.

        :return: A new :class:`.Align` instance containing the simulated
            data unless *dest* is specified (otherwise ``None``).

        .. note::

            The method :meth:`.iter_simul` can be much more efficient.
            For performance-critical applications, its use is
            recommended.
        """

        # run the simulation
        self._coalesce.simul(self._params._params,
                             self._params._random,
                             True)

        # reset the local align instance
        self._align._ns, self._align._ns_o = self._params._params.get_nsam(self._align._ns, self._align._ns_o)
        self._align._ng = 2

        # return a deep copy of the alignment
        return _interface.Align.create(_interface.Align._create_from_data_holder(self._coalesce.data()))

    ####################################################################

    def iter_simul(self, nrepet, cs=None, cs_filter=None, cs_structure=None, dest_trees=None, **feed_params):

        """
        Performs simulations based on the current value of parameters
        stored in :attr:`.params`. Return an iterator that will loop
        over the requested number of simulations. The simulated
        alignments are always available at each iteration loop as the
        instance attribute :attr:`.align`. By default, all iterations
        return a reference to this instance, but if *cs* is specified a
        dictionary of statistics is returned at each simulation.

        :param nrepet: number of simulations. If ``None``, iterate for
            ever (you are required to use ``break`` in loop with your
            own stopping criterion is reach).
        :param cs: :class:`.ComputeStats` instance properly configured
            to compute statistics from simulated data. If *cs* is specified,
            each iteration round will yield the dictionary of statistics
            obtained from :meth:`.ComputeStats.process_align`. See also
            *cs_filter* and *cs_struct*.
        :param cs_filter: a :class:`.Filter` instance to be passed to
            :class:`.ComputeStats` (ignored if *cs* is ``None``). By default,
            use :data:`~stats.filter_default`.
        :param cs_structure: a :class:`.Structure` instance to be passed to
            :class:`.ComputeStats` (ignored if *cs* is ``None``).
        :param dest_trees: a :class:`list` in which simulated trees will
            be appended. Since each simulation can yield several trees
            (because of recombination), a new sub-list will be appended
            for each simulation, with one item for each tree. Each tree
            covered a defined region of the simulated chromosome, so
            each item of each sub-list will be ``(tree, start, stop)``
            tuple with ``tree`` as a new :class:`.Tree` instance, and
            ``start`` and ``stop`` as the start and stop positions (real
            numbers comprises between 0 and 1). If recombination
            occurred, trees will not be sorted within their sub-list. If
            recombination did not occur, each simulation will be
            represented by a list with a single tuple. Any previously
            data contained in *dest_tree* is left untouched.
        :param feed_params: other arguments provide sequences of values
            for any parameter (except ``num_pop`` and ``seed``),
            allowing to modify any set of parameters between
            simulations. All changes are permanent and affect any later
            simulation. All additional options must be in the
            ``key=value`` format, with have one of the parameters as
            ``key``, and a sequence of values as ``value``. All
            sequences must be of length at least equal to *nrepet*. If
            longer, additional values are ignored. For list-based
            parameters, each value must be a ``(index, value)`` tuple
            and, for matrix-based parameters, each value must be a
            ``(index1, index2, value)`` tuple.

        :return: An iterator, yielding a reference to :attr:`.align` if
            *cs* is ``None``, or otherwise a dictionary of computed
            statistics.
        """

        # reset the local align instance
        self._align._ns, self._align._ns_o = self._params._params.get_nsam(self._align._ns, self._align._ns_o)
        self._align._ng = 2

        # safety checking
        for key, values in feed_params.iteritems():
            if key not in self._params._keys: raise ValueError, 'invalid parameter: `{0}`'.format(key)
            if len(values) < nrepet: ValueError, 'not enough values provided for: `{0}`'.format(key)

        # take default filter
        if cs is not None and cs_filter is None:
            cs_filter = _filter.filter_default

        # main loop
        i = 0
        stop = -1 if nrepet is None else nrepet
        while i != nrepet:

            # set the provided params
            for key, values in feed_params.iteritems():
                value = values[i]
                self._params._set(key, value)

            # run the simulation
            self._coalesce.simul(self._params._params,
                                 self._params._random,
                                 False)

            # export tree if needed
            if dest_trees != None:
                dest_trees.append([
                    ( _tree.Tree._from_coalesce(self._coalesce, i, str), 
                      self._coalesce.tree_start(i),
                      self._coalesce.tree_stop(i) )
                        for i in xrange(self._coalesce.number_of_trees()) ])

            # mutate
            self._coalesce.mutate()

            # return what is requested
            if cs != None: yield cs.process_align(self._align, filtr=cs_filter, struct=cs_structure)
            else: yield self._align

            i += 1

        # reset object
        self._align.reset()

    ####################################################################

    @property
    def align(self):

        """
        A :class:`.Align` instance containing simulated data for the
        current iteration of :meth:`.iter_simul`. The data in this
        instance will be updated at each iteration round, and deleted at
        the end of the iteration. If it must be copied, a deep copy is
        required (typically with :meth:`.Align.create`).
        """

        return self._align
