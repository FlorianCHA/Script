"""
Class to perform machine learning fitting using neural networks.
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

from .. import _eggwrapper
import os, sys, math

########################################################################

class Nnet(object):

    """
    Training and prediction with neural networks. This classes
    implements the back propagation algorithm for training neural
    network. In the documentation, the term *parameters* is used to
    design input variables, and *statistics* for output variables.

    Data can be loaded with methods :meth:`.parse` and/or :meth:`.load`.
    """

    def __init__(self):
        self._network = _eggwrapper.Network()
        self._data = _eggwrapper.Data()
        self._predict_data = _eggwrapper.Data()
        self._random = _eggwrapper.Random()
        self.reset()

    ####################################################################

    def reset(self):

        """
        Clear all loaded data.
        """

        self._data.set_num_patterns(0)
        self._data.set_num_input(0)
        self._data.set_num_output(0)
        self._init = False

    ####################################################################

    @staticmethod
    def _process_line(line):

        """
        Read a single line. Return two lists: one with the parameter
        values and one with the statistics values.
        """

        if line == '': raise ValueError, 'empty file: {0}'.format(fname)
        if line.count('#') != 1: raise ValueError, 'cannot read data from this filee: {0}'.format(fname)
        params, stats = line.split('#')
        try:
            params = map(float, params.split())
            stats = map(float, stats.split())
        except ValueError:
            raise ValueError, 'cannot read data from this file: {0}'.format(fname)
        return params, stats

    ####################################################################

    def parse(self, fname, num_max):

        """
        Import data from a provided file. All previously loaded data are
        discarded.

        :param fname: The name of an input file, containing data to be
            used for fitting. The data are formatted as follow: one line
            per pattern; each line contains space- or tab-separated
            values; the values are all parameter values, followed by a
            hash (``#``) as separator, followed by all statistic values.
            The order of parameters and statistics is required to be
            consistent over all lines.
        :param num_max: Maximal number of patterns. If the file contains
            less patterns than this number, all patterns are imported
            (and this does not count as an error). If the file contains
            more patterns, only the first *num_max* of them are
            imported.

        :return: The number of patterns actually imported.
        """

        # init
        self.reset()
        if num_max <= 0: raise ValueError, 'the value of `num_max` must be at least 1'
        if not os.path.isfile(fname): raise ValueError, 'cannot open this file: {0}'.format(fname)
        f = open(fname)

        # read a line and set the data tables sizes
        params, stats = self._process_line(f.readline())
        np = len(params)
        ns = len(stats)
        if np == 0: raise ValueError, 'invalid file (number of parameters is 0): {0}'.format(fname)
        if ns == 0: raise ValueError, 'invalid file (number of statistics is 0): {0}'.format(fname)
        self._data.set_num_input(ns)
        self._data.set_num_output(np)
        self._data.set_num_patterns(num_max) # it can be more efficient to set this last

        # process the first line
        for i,v in enumerate(stats): self._data.set_input(0, i, v)
        for i,v in enumerate(params): self._data.set_output(0, i, v)
        cnt = 1

        # process other lines
        for line in f:
            params, stats = self._process_line(line)
            if len(params) != np: raise ValueError, 'invalid line {0} of file {1}: inconsistent number of parameters'.format(cnt+1, fname)
            if len(stats) != ns: raise ValueError, 'invalid line {0} of file {1}: inconsistent number of stats'.format(cnt+1, fname)
            for i,v in enumerate(stats): self._data.set_input(cnt, i, v)
            for i,v in enumerate(params): self._data.set_output(cnt, i, v)
            cnt += 1
            if cnt == num_max: break

        # reset dataset size if not all memory has been used
        if num_max > cnt: self._data.set_num_patterns(cnt)

        # finish
        f.close()
        self._data.normalize()
        return cnt

    ####################################################################

    def load(self, patterns):

        """
        Pass data. All previously loaded data are discarded.

        :param patterns: A list of ``(params, stats)`` tuples. All sets
            of parameters must have the same, non-null, number of items,
            and the same goes for sets of statistics.
        """

        # init
        self.reset()
        if len(patterns) == 0: raise ValueError, 'can\'t set the number of patterns to 0'
        if len(patterns[0]) != 2: raise ValueError, 'invalid list of patterns: first pattern does not have two (params, stats) items'
        np = len(patterns[0][0])
        ns = len(patterns[0][1])
        if np == 0: raise ValueError, 'invalid list of patterns (number of parameters is 0)'
        if ns == 0: raise ValueError, 'invalid list of patterns (number of statistics is 0)'

        # read a line and set the data tables sizes
        self._data.set_num_input(ns)
        self._data.set_num_output(np)
        self._data.set_num_patterns(len(patterns)) # it can be more efficient to set this last

        # process other lines
        for idx, pattern in enumerate(patterns):
            if len(pattern) != 2: raise ValueError, 'invalid pattern {0}: does not have two (params, stats) items'.format(idx+1)
            params, stats = pattern
            if len(params) != np: raise ValueError, 'invalid pattern {0}: inconsistent number of parameters'.format(idx+1, fname)
            if len(stats) != ns: raise ValueError, 'invalid pattern {0}: inconsistent number of stats'.format(idx+1, fname)
            for i,v in enumerate(stats): self._data.set_input(idx, i, v)
            for i,v in enumerate(params): self._data.set_output(idx, i, v)

        # normalize
        self._data.normalize()

    ####################################################################

    @property
    def dim(self):

        """
        A three-item tuple providing, in that order: the number of
        loaded patterns; the number of parameters; the number of
        statistics.
        """

        return self._data.get_num_input(), self._data.get_num_output()

    ####################################################################

    def init_network(self, num_train=None, num_pred=None, num_hidden=5,
              fun1='log', fun2='linear', rate1=0.1, rate2=0.1,
              bound=100.0, momentum=0.001, init1=1.0, init2=1.0):

        """
        Setup training parameters. The default argument values give the
        default parameters.

        :param num_train: Number of patterns to use for training (the
            rest will be used for prediction). If ``None`` (default),
            the value will be set to half the number of patterns.

        :param num_pred: Number of patterns to use for prediction. If
            ``None`` (default), the value will be set to the rest of
            patterns.

        :param num_hidden: Number of neurons in the hidden layer of the
            network.

        :param fun1: Activation function to use for neurons of the
            hidden layer. Must be one of ``'linear'``, ``'log'`` and
            ``'tanh'``.

        :param fun2: Activation function to use for output neurons. Must
            be one of ``'linear'``, ``'log'`` and ``'tanh'``.

        :param rate1: Rate for weights from input to hidden layer.

        :param rate2: Rate for weights form hidden layer to output.

        :param bound: Extreme value (for both signs) for weight values.

        :param momentum: Proportion of the previous weight change to
            apply to all weight changes (use 0.0 to skip momentum).

        :param init1: bound for initialization values of weights from
            input to hidden layer.

        :param init2: bound for initialization values of weights from
            hidden layer to output.
        """

        if num_train != None and num_train < 1: raise ValueError, '`num_train` cannot be <1'
        if num_pred != None and num_pred < 0: raise ValueError, '`num_pred` cannot be <0'
        if num_hidden < 1: raise ValueError, 'invalid number of hidden layer neurons'
        if fun1 == 'linear': fun1 = _eggwrapper.linear
        elif fun1 == 'log': fun1 = _eggwrapper.log
        elif fun1 == 'tanh': fun1 = _eggwrapper.tanh
        else: raise ValueError, 'invalid activation function: {0}'.format(fun1)
        if fun2 == 'linear': fun2 = _eggwrapper.linear
        elif fun2 == 'log': fun2 = _eggwrapper.log
        elif fun2 == 'tanh': fun2 = _eggwrapper.tanh
        else: raise ValueError, 'invalid activation function: {0}'.format(fun2)
        if rate1 <= 0.0: raise '`rate1` argument is too low'
        if rate2 <= 0.0: raise '`rate2` argument is too low'

        # check and convert arguments
        if self._data.get_num_patterns() == 0: raise ValueError, 'cannot train network: not enough data loaded'
        if self._data.get_num_input() == 0: raise ValueError, 'cannot train network: not enough data loaded'
        if self._data.get_num_output() == 0: raise ValueError, 'cannot train network: not enough  data loaded'
        if num_train == None: num_train = int(math.ceil(self._data.get_num_patterns() / 2.0))
        elif num_train >= self._data.get_num_patterns(): raise ValueError, 'invalid number of training patterns'
        if num_pred == None: num_pred = self._data.get_num_patterns() - num_train
        elif num_train+num_pred >= self._data.get_num_patterns(): raise ValueError, 'invalid number of prediction patterns'

        # create & initialize network
        self._network.setup(self._data, 0, num_train, self._data,
                num_train, num_train + num_pred, num_hidden, fun1,
                fun2, rate1, rate2, bound, momentum, self._random)
        self._network.init_weights(init1, init2)
        self._init = True
        self._num_hidden = num_hidden

    ####################################################################

    def train(self, n):

        """
        Train network using loaded data and parameter values specified
        using :meth:`.setup` (or default values).

        :param n: Number of training round.

        :return: A tuple containing the training/testing errors.
        """

        if not self._init: raise ValueError, 'it is required to call `init_network()` before training'
        self._network.train(n)
        return self._network.training_error(), self._network.testing_error()

    ####################################################################

    def get_weights(self):

        """
        Yields current values of weights. The return value a tuple of
        two lists (one for the input-to-hidden layer weights, and the
        other for the hidden-to-output layer weights). Each list
        contains sublists giving, for each input and hidden node
        (respectively), the weights for each of the nodes of the next
        layer. For example, for getting the weight of the connexion of
        node *i* of the input layer to node *j* of the hidden layer, use
        the subscript ``[0][i][j]`` and for the weight of the connextion
        of node *i* of the hidden layer to node *j* of the output layer,
        use the subscript ``[1][i][j]``. There is an addition index for
        the bias neuron at both levels.
        """

        if not self._init: raise ValueError, 'it is required to call `init_network()` before accessing weights'

        return ([[self._network.weight_hidden(j, i)
                    for j in range(self._num_hidden)]
                        for i in range(self._data.get_num_input()+1)],
                [[self._network.weight_output(j, i)
                    for j in range(self._data.get_num_output())]
                        for i in range(self._num_hidden+1)])

    ####################################################################

    def predict(self, data):

        """
        Predict output values from provided input values using the
        current state of the network (whatever level training has been
        previously achieved).
        """

        if not self._init: raise ValueError, 'it is required to call `init_network()` before predicting data'

        n = len(data)
        ni = self._data.get_num_input()
        no = self._data.get_num_output()

        self._predict_data.set_num_input(ni)
        self._predict_data.set_num_output(no)
        self._predict_data.set_num_patterns(n)

        for i, item in enumerate(data):
            if len(item) != ni: raise ValueError, 'number of input variables is not matching'
            for j, v in enumerate(item):
                self._predict_data.set_input(i, j, v)

        for i in range(ni):
            self._predict_data.normalize_input(i, self._data.mean_input(i), self._data.std_input(i))

        for i in range(n):
            self._network.predict(self._predict_data, i, False)
            for j in range(no):
                self._predict_data.set_output(i, j, self._network.prediction(j))

            self._predict_data.unnormalize_output(i, self._data.mean_output(i), self._data.std_output(i))

        return [[self._predict_data.get_output(i, j) for j in range(no)] for i in range(n)]
