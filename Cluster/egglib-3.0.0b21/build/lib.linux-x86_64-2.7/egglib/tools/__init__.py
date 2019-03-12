"""
This module contains pure python classes or functions that can be of
general or specific use.
"""

__license__ = """
    Copyright 2008-2012,2014-2015 Stephane De Mita, Mathieu Siol

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

from ._discretize import Discretize
from ._random import Random
from ._concat import concat
from ._code_tools import Translator, translate, int2codon, orf_iter, \
                        longest_orf, backalign, BackalignError, \
                        trailing_stops, iter_stops, has_stop
from ._reading_frame import ReadingFrame
from ._seq_manip import rc, compare, regex, motif_iter, ungap
