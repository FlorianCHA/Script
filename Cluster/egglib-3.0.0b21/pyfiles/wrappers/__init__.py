"""
This module contains wrappers to external applications.
"""

__license__ = """
    Copyright 2009-2011,2016 Stephane De Mita, Mathieu Siol

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

from ._phyml import phyml
from ._codeml import codeml
from ._muscle import muscle
from ._clustal import clustal
from ._utils import paths
from ._blast import blastn, blastp, tblastn, tblastx, blastx, makeblastdb, blast_config
paths.load() # must be after loading all application modules
