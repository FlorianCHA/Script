"""
EggLib Python Package
=====================

EggLib is a Python package for population genetics and genomics
incorporating sequence alignment management tools and wrappers for a few
external applications. The EggLib C++ library is wrapped within the
package and can be used through secure wrappers.
"""


__license__ = """
    Copyright 2008-2016 Stephane De Mita, Mathieu Siol

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

version = '3.0.0b21'                         #EGGVERSION#
""" egglib version number. """

from . import _eggwrapper
from . import tools
from . import stats
from . import io
from . import fit
from . import coalesce
from . import wrappers
from ._interface import Align, Container, SampleView, SequenceView, GroupView
from ._tree import Tree, Node
