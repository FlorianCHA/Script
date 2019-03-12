"""
This module contains the function to estimate to misorientation
probability based on Baudry and Depaulis (2003).
"""

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
from . import _site, _freq, _filter

########################################################################

class ProbaMisoriented(object):

    """
    Estimate the misorientation probability from a user-provided set
    of sites with a fixed outgroup. Only sites that are either variable
    within the ingroup or have a fixed difference with respect to the
    outgroup are considered. Sites with more than two different alleles
    in the ingroup, or more than one allele in the outgroup, are
    ignored.

    This function is an implementation of the method mentioned in Baudry
    and Depaulis (2003), allowing to estimate the probability that a
    site oriented using the provided outgroup have be misoriented due to
    a homoplasic mutation in the branch leading to the outgroup. Note
    that this estimation neglects the probability of shared
    polymorphism.

    Reference: Baudry, E. & F. Depaulis. 2003. Effect of misoriented
    sites on neutrality tests with outgroup.
    *Genetics*\ \xa0\  **165**\ : 1619-1622.

    :param align: a :class:`.Align` containing the sites to analyse.

    If the instance is created with an alignment as constructor argument,
    then the statistics are computed. The method :meth:`.load_align`
    does the same from an existing :class:`.ProbaMisoriented` instance.
    Otherwise, individual sites can be loaded with :meth:`.load_site`,
    and then the statistics can be computed using :meth:`.compute` (the
    latter is preferable if generate :class:`.Freq` instances for an
    other use).

    .. versionadded:: 3.0.0
    """

    ####################################################################

    def __init__(self, align=None):

        self._sd = _eggwrapper.SiteDiversity()
        self._site = _site.Site()
        self._freq = _freq.Freq()
        self.reset()
        if align is not None: self.load_align(align)

    ####################################################################

    def reset(self):

        """
        Clear all loaded or computed data.
        """

        self._S = 0
        self._D = 0
        self._Ti = 0
        self._Ti_cnt = 0
        self._TiTv = None
        self._pM = None
        self._M = 0
        self._M_cnt = 0

    ####################################################################

    def load_align(self, align):

        """
        Load all sites of align that meet criteria. If there are previously
        loaded data, they are discarded. This method computes statistics.
        Data are required to be DNA sequences.

        :param align: an :class:`.Align` instance.
        """

        self.reset()
        for i in xrange(align._obj.get_nsit()):
            self._site.process_align(align, i, filtr=_filter.filter_dna, struct=None, reset=True)
            self._freq.process_site(self._site)
            self.load_site(self._freq)
        self.compute()

    ####################################################################

    def load_site(self, freq):

        """
        Load a single site. If there are previously loaded data, they
        are retained. To actually compute the misorientation
        probability, the user must call :meth:`.compute`.

        :param site: a :class:`.Freq` instance.
        """

        flag = self._sd.process(freq._obj)
        if (flag&2) == 0: return
        if self._sd.Aglob() < 2: return # skip sites that are fixed overall
        if self._sd.Aing() > 2: return # skip sites with 3+ alleles ingroup
        if self._sd.Aout() > 1: return # skip sites with 2+ alleles outgroup

        if self._sd.Aing() == 2:
            self._S += 1
            if freq._obj.frq_outgroup().nseff() > 0:
                self._M_cnt += 1
                if self._sd.Aglob() > self._sd.Aing():
                    self._M += 1
        else:
            self._D += 1

        if self._sd.Aglob() == 2:
            self._Ti_cnt += 1
            if     self._sd.global_allele(0) == 65 or self._sd.global_allele(0) == 97:
                if self._sd.global_allele(1) == 71 or self._sd.global_allele(1) == 103: self._Ti += 1
            elif   self._sd.global_allele(0) == 71 or self._sd.global_allele(0) == 103:
                if self._sd.global_allele(1) == 65 or self._sd.global_allele(1) == 97:  self._Ti += 1
            elif   self._sd.global_allele(0) == 67 or self._sd.global_allele(0) == 99:
                if self._sd.global_allele(1) == 84 or self._sd.global_allele(1) == 116: self._Ti += 1
            elif   self._sd.global_allele(0) == 84 or self._sd.global_allele(0) == 116:
                if self._sd.global_allele(1) == 67 or self._sd.global_allele(1) == 99:  self._Ti += 1

    ####################################################################

    def compute(self):

        """
        Compute :attr:`.pM` and :attr:`.TiTv` statistics. Requires that
        sites have been loaded using :meth:`.load_site`. This method
        does not reset the instance.
        """

        self._pM = None
        self._TiTv = None

        if self._Ti_cnt > 0 and self._M_cnt > 0:
            a = float(self._Ti) / self._Ti_cnt
            b = (1-a)/2
            if b > 0:
                self._TiTv = a/b
                self._pM = ((float(self._M)/self._M_cnt) * a**2 + 2*b**2) / (2*b*(2*a+b))

    ####################################################################

    @property
    def S(self):

        """
        Number of loaded polymorphic sites (within the ingroup).
        """

        return self._S

    ####################################################################

    @property
    def D(self):

        """
        Number of loaded sites with a fixed difference with respect to
        the outgroup.
        """

        return self._D

    ####################################################################

    @property
    def TiTv(self):

        """
        Ratio of transition and transversion rates ratio. ``None`` if
        the value cannot be computed (no loaded data or null
        transversion rate). Requires that :meth:`.compute` has been
        called.
        """

        return self._TiTv

    ####################################################################

    @property
    def pM(self):

        """
        Probability of misorientation. ``None`` if the value cannot be
        computed (no loaded data, no valid polymorphism, null
        transversion rate). Requires that :meth:`.compute` has been
        called.
        """

        return self._pM
