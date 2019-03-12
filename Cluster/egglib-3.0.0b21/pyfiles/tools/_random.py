"""
This module contains the single class :class:`Random` which is
available from the package :mod:`tools`.
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

########################################################################

class Random(object):

    """
    Pseudo-random number generator.

    This class implements the Mersenne Twister algorithm for
    pseudo-random number generation. It is based on work by `Makoto
    Matsumoto and Takuji Nishimura
    <http://www.math.sci.hiroshima-u.ac.jp/~m-MAT/MT/emt.html>`_ and
    `Jasper Bedaux <http://www.bedaux.net/mtrand/>`_ for the core
    generator, and the :cpp:class:`Random` class of Egglib up to 2.2 for
    conversion to other laws than uniform.

    Note that different instances of the class have independent chains
    of pseudo-random numbers. If several instances have the same seed,
    they will generate the exact same chain of pseudo-random numbers.
    Note that this applies if the default constructor is used and that
    instances are created within the same second.

    All non-uniform distribution laws generators are based either on the
    :meth:`~.integer_32bit` or the standard (half-open, 32 bit)
    :meth:`~.uniform` methods.

    :param seed: The constructor accepts an optional integer value to
        seed the pseudo-random number generator sequence. By default,
        the current system clock value will be used, which means that
        all instances created within the same second will generate
        strictly identical sequences. Favor large, high-complexity
        seeds. When using different instances of this class in a
        program, or different instances of the same program launched
        simultaneously, ensure they are all seeded using different
        seeds. 
    """

    ####################################################################

    def __init__(self, seed = None):
        if seed == None: self._obj = _eggwrapper.Random()
        elif isinstance(seed, (int, long)): self._obj = _eggwrapper.Random(seed)
        else: raise TypeError, 'invalid type for `seed`'

    ####################################################################

    def boolean(self):

        """
        Draw a boolean with equal probabilities (*p* = 0.5).

        :return: A boolean.
        """

        return self._obj.brand()

    ####################################################################

    def bernoulli(self, p):

        """
        Draw a boolean with given probability.

        :param p: probability of returning ``True``.
        :return: A boolean.
        """

        if p < 0 or p > 1: raise ValueError, 'probability `p` must be >=0 and <=1'
        return self._obj.uniform() < p

    ####################################################################

    def binomial(self, n, p):

        """
        Draw a value from a binomial distribution.

        :param n: Number of tests (>=0).
        :param p: Test probability (>=0 and <=1).
        :return: An :class:`int` (number of successes).
        """

        if n < 0: raise ValueError, 'number of tests `n` must be >=0'
        if p < 0 or p > 1: raise ValueError, 'probability `p` must be >=0 and <=1'
        return self._obj.binomrand(n, p)

    ####################################################################

    def exponential(self, expectation):

        """
        Draw a value from an exponential distribution.

        :param expectation: Distribution's mean (equal to
            1/:math:`\lambda` , if :math:`\lambda` is the rate
            parameter). Required to be >0.
        :return: A :class:`long`.
        """

        if expectation <=0: raise ValueError, 'argument `expectation` must be strictly positive'
        return self._obj.erand(expectation)

    ####################################################################

    def geometric(self, p):

        """
        Draw a value from a geometric distribution.

        :param p: Geometric law parameter (>0 and <=1).
        :return: A positive :class:`int`.
        """

        if p<=0 or p>1: raise ValueError, 'probability `p` must be >0 and <=1'
        return self._obj.grand(p)

    ####################################################################

    def normal(self):

        """
        Draw a value from the normal distribution with expectation 0
        and variance 1. The expression ``rand.normal() * sd + m`` can be
        used to rescale the drawn value to a normal distribution with
        expectation *m* and standard deviation *sd*.

        :return: A :class:`float`.
        """

        return self._obj.nrand()

    ####################################################################

    def poisson(self, p):

        """
        Draw a value from a Poisson distribution.

        :param p: Poisson distribution parameter (usually noted
            :math:`\lambda`). Required to be >0
        :return: A positive :class:`int`.
        """

        if p <= 0: raise ValueError, 'argument `p` must be strictly positive'
        return self._obj.prand(p)

    ####################################################################

    def integer(self, n):

        """
        Draw an integer from a uniform distribution.

        :param n: Number of possible values). Note that this number is
            excluded and will never be returned. Required to be a
            stricly positive integer.
        :return: An `int` in range ``[0, n-1]``.
        """

        if n <=0: raise ValueError, 'argument `n` must be strictly positive'
        return self._obj.irand(n)

    ####################################################################

    def integer_32bit(self):

        """
        Generate a 32-bit random integer.

        :return: An :class:`long` in the interval [0, 4294967295] (that
            is in the interval [0, 2^32-1].
        """

        return self._obj.rand_int32()

    ####################################################################

    def uniform(self):

        """
        Generate a :class:`float` in the half-open interval [0,1) with
        default 32-bit precision. The value 1 is not included.
        """

        return self._obj.uniform()

    ####################################################################

    def uniform_53bit(self):

        """
        Generate a :class:`float` in the half-open interval [0,1) with
        increased 53-bit precision. The value 1 is not included.

        The increased precision increases the number of possible values
        (2^53 = 9007199254740992 instead of 2^32 = 4294967296). This
        comes with the cost of increased computing time.
        """

        return self._obj.uniform53()

    ####################################################################

    def uniform_closed(self):

        """
        Generate a :class:`float` in the closed interval [0,1] with
        default 32-bit precision. Both limits are included.
        """

        return self._obj.uniformcl()

    ####################################################################

    def uniform_open(self):

        """
        Generate a :class:`float` in the open interval (0,1) with
        default 32-bit precision. Both limits are excluded.
        """

        return self._obj.uniformop()

    ####################################################################

    def get_seed(self):

        """
        Get the seed value (value used at object creation, or value used
        to reset the instance with :meth:`.set_seed`.

        :return: A :class:`long`.
        """

        return self._obj.get_seed()

    ####################################################################

    def set_seed(self, seed):

        """
        Reset the instance by providing a seed value. See the comments
        on the seed value in the class description.
        """

        if not isinstance(seed, (int, long)): raise TypeError, 'invalid type for `seed`'
        self._obj.set_seed(seed)

########################################################################
