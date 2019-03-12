"""
This module contains the single class :class:`Discretize` which is
available from the package :mod:`tools`.
"""

__license__ = """
    Copyright 2013-2015 Stephane De Mita, Mathieu Siol

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

import operator

########################################################################

class Discretize(object):

    """
    This class discretizes a continuous distribution based on a given
    number of samples with an arbitrary number of dimensions.

    The default constructor builds an empty instance which can
    process data using :meth:`~.Discretize.process`. It is possible to
    process data directly at object construction time by passing
    arguments to the constructor. In that case, arguments must match the
    syntax of the :meth:`~.Discretize.process` method.

    .. versionadded:: 3.0.0
    """

    ####################################################################


    def __init__(self, *args, **kwargs):

        self.reset()
        if len(args) + len(kwargs) > 0:
            self.process(*args, **kwargs)


    ####################################################################

    def reset(self):

        """
        Clear all data and reset the object to default settings.
        """

        self._ndim = None
        self._ncat = None
        self._ndata = None
        self._noob = None
        self._bounds = None
        self._check_bounds = None
        self._dist = []
        self._cache_idx = {}

    ####################################################################

    def process(self, data, ncat, bounds=None, allow_oob=False):

        """
        Process a data set and generate a discretized distribution.

        :param data: a sequence of fixed-length sequences of numeric
            values. Typically, one will pass a list of *d*-length lists
            or tuples containing floats (but integers are also
            supported), where *d* is the number of dimensions of the
            distribution. It is required that all items of *data* have
            the same length, and that this length is at least 1. The
            number of dimensions is read from the first item of *data*.

        :param ncat: numbers of categories used for binarization. The
           value can be a single integer or a sequence of integers
           whose length matches the number of dimensions of the data
           set. If an integer is passed, it must be a strictly positive
           integer and it will be used as the number of categories for
           all dimensions. Otherwise, each value describes the number of
           categories for each corresponding dimension. There is no
           default; this argument is required.

        :param bounds: bounds of the distribution given as minimum and
           maximum values for all dimensions. The position of category
           limits are determined by these bounds. The value can either
           be (1) ``None``, (2) a sequence of two integers, or (3) a
           sequence of sequences of two integers. In case (1), the
           actual minimum and maximum values of all dimensions are used
           a distribution bounds. In case (2), the values are used as
           minimum and maximum (respectively) for all dimensions. In
           case (3), the length of the sequence must be equal to the
           number of dimensions. Each item is used to define the minimum
           and maximum (respectively) for the corresponding dimension.
           In addition, it is possible to replace any item of this
           sequence ``None`` to specify that the actual minimum and
           maximum values must be used as distribution bounds for the
           corresponding dimensions. By default, the actual extreme
           values of the data are used. In that case, it is guaranteed
           that there will be no out of bounds values.

        :param allow_oob: a boolean indicating what to do if a value is
          out of bounds for any of the dimensions. This can happen only
          if the argument *bounds* has been set to a non-default value.
          If out of bounds values are allowed, they will be ignored and
          a counter will be incremented. Otherwise, a
          :class:`ValueError` will be raised on the first occurrence of
          an out of bound value.
        """

        self.reset()

        # get data set dimensions

        if len(data) == 0: raise ValueError, 'cannot discretize data: no data'
        if len(data[0]) == 0: raise ValueError, 'cannot discretize data: first item is empty'
        self._ndim = len(data[0])

        # get number of categories

        if isinstance(ncat, int):
            if ncat <= 0: raise ValueError, 'cannot discretize: null or negative number of categories: {0}'.format(ncat)
            self._ncat = [ncat] * self._ndim
        else:
            for i, v in enumerate(ncat):
                if v <= 0: raise ValueError, 'cannot discretize: null or negative number of categories for dimension {0}: {0}'.format(i+1, v)
            self._ncat = list(ncat) # deep copy

        # get bounds

        if bounds == None: # none
            self._bounds = [None] * self._ndim
            self._check_bounds = False
        else:
            try:
                MIN, MAX = map(float, bounds)
            except (ValueError, TypeError): # ndim-length of pairs of values
                if len(bounds) != self._ndim: raise ValueError, 'cannot discretize: invalid number of items in `bounds` argument: {0}'.format(len(bounds))
                self._bounds = [None] * self._ndim
                self._check_bounds = False # unless there is at least one non-None
                for i, bound in enumerate(bounds):
                    if bound != None: 
                        try: MIN, MAX = map(float, bound)
                        except ValueError: raise ValueError, 'cannot discretize: invalid item in `bounds` argument: {0}'.format(bound)
                        if MIN >= MAX: raise ValueError, 'cannot discretize: invalid bounds (MIN>=MAX)'
                        self._bounds[i] = MIN, MAX, MAX-MIN
                        self._check_bounds = True
            else: # single pair of floats (no exception above)
                if MIN >= MAX: raise ValueError, 'cannot discretize: invalid bounds (MIN>=MAX)'
                self._bounds = [(MIN, MAX, MAX-MIN) for i in xrange(self._ndim)]
                self._check_bounds = True

        # pre-screen data to determine bounds if needed

        if None in self._bounds:
            pos = filter(lambda x: self._bounds[x] == None, xrange(self._ndim))
            for p in pos: self._bounds[p] = [data[0][p]] * 2
            for item in data:
                for p in pos:
                    if item[p] < self._bounds[p][0]: self._bounds[p][0] = item[p]
                    if item[p] > self._bounds[p][1]: self._bounds[p][1] = item[p]
            for p in pos:
                MIN, MAX = self._bounds[p]
                if MIN==MAX: raise ValueError, 'cannot discretize: no variation on dimension {0}'.format(p+1)
                self._bounds[p] = MIN, MAX, MAX-MIN

        # initialize counters

        self._ndata = 0
        self._noob = 0
        
        # initialize the distribution as a single-dimension list 

        self._dist = [0] * reduce(operator.mul, self._ncat)

        # process all values

        for idx, item in enumerate(data):

            # check bounds

            if self._check_bounds:
                skip = False
                for i, v in enumerate(item):
                    if v < self._bounds[i][0] or v > self._bounds[i][1]:
                        if not allow_oob: raise ValueError, 'cannot discretize: value out of bound: {0}'.format(item)
                        skip = True
                if skip:
                    self._noob += 1
                    continue

            # compute index tuple as (X-min)/(max-min) with special case if X=MAX

            cat = tuple(ncat-1 if val==MAX else int(ncat * (val-MIN)/SPAN)
                    for (val, (MIN, MAX, SPAN), ncat) in zip(item, self._bounds, self._ncat))

            self._dist[self._index(cat)] += 1
            self._ndata += 1

    ####################################################################

    @property
    def bounds(self):

        """
        Bounds of the final distribution. By default (if no data have
        been processed), returns ``None``. If data have been processed,
        returns a list of length matching :attr:`~.ndim` containing
        ``(min, max)`` tuples giving the extreme values for each
        dimension. If :attr:`~.bounds` have been set, returns a copy of
        the passed value regardless of whether data have been processed
        or not.
        """

        if self._ndim == None: raise ValueError, 'cannot access `Discretize` data: no data'
        return [(MIN, MAX) for (MIN, MAX, SPAN) in self._bounds]

    ####################################################################

    @property
    def ndim(self):

        """
        Number of dimensions of the processed data set. The value is
        ``None`` if no data set has been processed.
        """

        return self._ndim

    ####################################################################

    @property
    def ncat(self):

        """
        Number of categories of the processed data set. The value is
        ``None`` if no data set has been processed. Otherwise the value
        is necessarily at least 1.
        """

        return self._ncat

    ####################################################################

    @property
    def ndata(self):

        """
        Number of data set items of the processed data set. The value is
        ``None`` if no data set has been processed. If out of bounds
        values were allowed, and did actually occur, this value will be
        the number of data set items that were included and used.
        """

        return self._ndata

    ####################################################################

    @property
    def nskipped(self):

        """
        Number of data set items that were skipped because they were
        out of bounds (only if non-default bounds were used). The value
        is ``None`` is no data set has been processed.
        """

        return self._noob

    ####################################################################

    def _index(self, idx):

        """
        This helpers converts a tuple of indexes (one for each
        dimension) to a single index for the one-dimensional data array.
        Should only be called if valid data have been loaded.
        """

        if idx in self._cache_idx: return self._cache_idx[idx]

        self._cache_idx[idx] = 0
        acc = 1
        for dim in xrange(self._ndim):
            self._cache_idx[idx] += idx[-1-dim] * acc
            acc *= self._ncat[-1-dim]
        return self._cache_idx[idx]

    ####################################################################

    def get(self, *idx):

        """
        Get a frequency value from discretized data. Positional
        arguments should be the categories indexes (one index per
        dimension). Passing invalid indexes or an invalid number of
        index values, or calling this method on an instance that does
        not contain data, result in a :class:`ValueError`. Indexes can
        be negative (to count from the end). Slicing is not supported.
        """

        if self._ndim is None: raise ValueError, 'cannot access `Discretize` data: no data'
        if len(idx) != self._ndim: raise ValueError, 'cannot access `Discretize` data: invalid number of indexes (expected {0}, got {1})'.format(self._ndim, len(idx))
        for i,x in enumerate(idx):
            y = x if x>=0 else self._ncat[i] + x
            if y<0 or y>=self._ncat[i]: raise ValueError, 'cannot access `Discretize` data: invalid index for dimension {0}: {1}'.format(i+1, x)
        return self._dist[self._index(idx)]

    ####################################################################

    def mget(self, *idx):

        """
        Get marginal values. This method is equivalent to
        :meth:`~.Discretize.get` except that it can integrate values
        over any dimension. To obtain marginal sums for a given set of
        dimensions (one or more) one should replace all other indexes by
        ``None``. This method requires, as :meth:`~.Discretize.get`
        does, one argument per dimension. If argument are provided they
        must be proper indexes. It is possible to replace any index by
        ``None``, including several at the same time, all of them
        (equivalent to :attr:`~.ndata`) and none of them (equivalent to
        :meth:`~.Discretize.get`).
        """

        if None not in idx: return self.get(*idx)
        dim = idx.index(None)
        ret = 0
        for i in range(self._ncat[dim]):
            idx2 = list(idx)
            idx2[dim] = i
            ret += self.mget(*idx2)
        return ret

    ####################################################################

    def mid(self, dim):

        """
        Get the list of category midpoints for a given dimension. Raises
        a :class:`ValueError` if no data are present in the instance or
        if the index is out of bounds.
        """

        if self._ndim == None: raise ValueError, 'cannot get category data: no data in instance'
        if dim < 0: raise ValueError, 'cannot get category data: negative indexes are not supported'
        if dim >= self._ndim: raise ValueError, 'cannot get category data: index out of bounds: {0}'.format(dim)

        return [self._bounds[dim][0] + self._bounds[dim][2]*(0.5+i)/self._ncat[dim] for i in range(self._ncat[dim])]

    ####################################################################

    def left(self, dim):

        """
        Get the list of category left points (minimum of each category)
        for a given dimension. Raises a :class:`ValueError` if no data
        are present in the instance or if the index is out of bounds.
        """

        if self._ndim == None: raise ValueError, 'cannot get category data: no data in instance'
        if dim < 0: raise ValueError, 'cannot get category data: negative indexes are not supported'
        if dim >= self._ndim: raise ValueError, 'cannot get category data: index out of bounds: {0}'.format(dim)

        return [self._bounds[dim][0] + self._bounds[dim][2]*(1.0*i)/self._ncat[dim] for i in range(self._ncat[dim])]

    ####################################################################

    def right(self, dim):

        """
        Get the list of category right points (maximum of each category)
        for a given dimension. Raises a :class:`ValueError` if no data
        are present in the instance or if the index is out of bounds.
        """

        if self._ndim == None: raise ValueError, 'cannot get category data: no data in instance'
        if dim < 0: raise ValueError, 'cannot get category data: negative indexes are not supported'
        if dim >= self._ndim: raise ValueError, 'cannot get category data: index out of bounds: {0}'.format(dim)

        return [self._bounds[dim][0] + self._bounds[dim][2]*(1.0+i)/self._ncat[dim] for i in range(self._ncat[dim])]
