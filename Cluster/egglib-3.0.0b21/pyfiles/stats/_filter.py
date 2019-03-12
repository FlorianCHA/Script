__license__ = """
    Copyright 2015-2016 Stephane De Mita, Mathieu Siol

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

class Filter(object):

    """
    Holds lists of valid (exploitable and missing) data codes. The
    constructor allows to immediately instanciate the class with all
    desired information. It is also possible to call the modifier
    methods to add information. Note that this module provides a few
    pre-defined instances with common settings that the user may use and
    possibly extend.

    .. note::
        If no codes are specified (which is the default behaviour if no
        arguments are specified), the filter will treat all codes as
        acceptable.

    :param exploitable: Value or values used as exploitable data.
        Possible values are ``None``, a single integer, a list (or other
        iterable) of integers, or a string. An iterable of one-character
        strings is also possible. Negative integers are supported.

    :param rng: A range of integer values taken as exploitable. Must
        contain two values, providing the minimum and maximum values.

    :param missing: Value or values used as missing data. Specifications
        as for *exploitable*.

    :param lower_alias: For all *exploitable* and *missing* (but not
        *rng*) values, use the lower case character as alias. Ignored
        for all integer values.

    :param exploitable_alias: Provide the alias values for exploitable
        data. The number of values and their order must match the
        argument *exploitable*.

    :param missing_alias: Provide the alias values for missing data. The
        number of values and their order must match the argument
        *missing*.
    """

    ####################################################################

    @staticmethod
    def _process(array):
        if isinstance(array, int): return [array]
        try: return [i if isinstance(i, int) else ord(i) for i in array]
        except TypeError, e: raise ValueError, 'invalid data: invalid type for exploitable/missing value found: {0}'.format(e.message)

    ####################################################################

    def __init__(self, exploitable=None, rng=None, missing=None,
                 exploitable_alias=None, missing_alias=None,
                 lower_alias=False):

        self._obj = _eggwrapper.Filter()

        if exploitable != None:
            self.exploitable(exploitable, alias=exploitable_alias, lower_alias=lower_alias)
        elif exploitable_alias != None:
            raise ValueError, 'cannot specify `exploitable_alias` if `exploitable` is `None`'

        if missing != None:
            self.missing(missing, alias=missing_alias, lower_alias=lower_alias)
        elif missing_alias != None:
            raise ValueError, 'cannot specify `missing_alias` if `missing` is `None`'

        if rng != None:
            try:
                if len(rng) != 2: raise ValueError, 'invalid `rng` argument: must have two items'
            except TypeError: raise TypeError, 'invalid `rng` argument'
            self.rng(*rng)

    ####################################################################

    def exploitable(self, values, alias=None, lower_alias=False):  # almost identical copies . . .

        """
        Add exploitable values. See the class documentation for more
        details.
        """


        if lower_alias == True and alias != None:
            raise ValueError, 'cannot specify both `alias` and `lower_alias` at the same time'

        values = self._process(values)
        if alias == None:
            for a in values:
                if lower_alias == True and a>=65 and a<=90:
                    self._obj.add_exploitable_with_alias(a, a+32)
                else:
                    self._obj.add_exploitable(a)
        else:
            alias = self._process(alias)
            if len(values) != len(alias): raise ValueError, 'length of `values` and `alias` don\'t match'
            for a, b in zip(values, alias): self._obj.add_exploitable_with_alias(a, b)

    ####################################################################

    def missing(self, values, alias=None, lower_alias=False): # almost identical copies . . .

        """
        Add missing values. See the class documentation for more
        details.
        """

        if lower_alias == True and alias != None:
            raise ValueError, 'cannot specify both `alias` and `lower_alias` at the same time'

        values = self._process(values)
        if alias == None:
            for a in values:
                if lower_alias == True and a>=65 and a<=90:
                    self._obj.add_missing_with_alias(a, a+32)
                else:
                    self._obj.add_missing(a)
        else:
            alias = self._process(alias)
            if len(values) != len(alias): raise ValueError, 'length of `values` and `alias` don\'t match'
            for a,b in zip(values, alias): self._obj.add_missing_with_alias(a, b)

    ####################################################################

    def rng(self, mini, maxi):

        """
        Add a range of exploitable data. See the class documentation for
        more details.
        """

        if mini > maxi: raise ValueError, 'minimum value of a range of exploitable values must not be greater than maximum'
        mini, maxi = self._process([mini, maxi])
        self._obj.add_exploitable_range(mini, maxi)

########################################################################

#: All possible allele values are accepted as valid.
#: This filter is case-sensitive.
filter_default = Filter()

#: Used for DNA sequences (case-independent).
#: Exploitable: ``ACGT``. Missing: ``RYSWKMBDHVN?-``.
#: Lower-case characters are mapped to the corresponding upper-case character.
filter_dna = Filter(exploitable='ACGT', missing='RYSWKMBDHVN?-', lower_alias=True)

#: Used for RNA sequences (case-independent).
#: Exploitable: ``ACGU``. Missing: ``RYSWKMBDHVN?-``. The IUPAC codes for
#: missing data are identical in meaning to those for DNA sequences.
#: Lower-case characters are mapped to the corresponding upper-case character.
filter_rna = Filter(exploitable='ACGU', missing='RYSWKMBDHVN?-', lower_alias=True)

#: Paranoid filter for DNA sequences. Exploitable: ``ACGT``. Missing: ``N``.
#: Lower-case characters are not allowed.
filter_strict = Filter(exploitable='ACGT', missing='N', lower_alias=False)

#: Used for strictly positive SSR alleles.
#: Exploitable: range 1--999. Missing: -1.
filter_ssr = Filter(rng=(1, 999), missing=-1)

#: Used for numerical values.
#: Exploitable: range -999--999. Missing: none.
filter_num = Filter(rng=(-999, 999))

#: Used for amino acid sequences. Exploitable: ``ACDEFGHIKLMNPQRSTVWY``.
#: Missing: ``X-``. 
#: Lower-case characters are mapped to the corresponding upper-case character.
filter_amino = Filter(exploitable='ACDEFGHIKLMNPQRSTVWY', missing='X-', lower_alias=True)

#: Used for numerical codon representation. Exploitable: range 0--63. Missing: 64.
filter_codon = Filter(rng=(0, 63), missing=64)
