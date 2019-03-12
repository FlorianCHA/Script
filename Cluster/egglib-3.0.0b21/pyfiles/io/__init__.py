"""
This module contains tools for import/export of several form of data in
several file formats.
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

from _export import to_ms
from _fasta import fasta_iter, from_fasta
from _gff3 import GFF3, GFF3Feature
from _vcf import VcfParser, Variant
from _genbank import GenBank, GenBankFeature, GenBankFeatureLocation
from _legacy import from_clustal, from_staden, from_genalys, get_fgenesh
