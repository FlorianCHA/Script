"""
This module contains tools to analyze compute statistics from genetic
variation data.
"""

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

from _ld import pairwise_LD, matrix_LD
from _haplotypes import haplotypes
from _coding_diversity import CodingDiversity
from _ehh import EHH
from _baudry import ProbaMisoriented
from _structure import Structure, get_structure, make_structure
from _site import Site, site_from_align, site_from_list, site_from_vcf
from _freq import Freq, freq_from_site, freq_from_list, freq_from_vcf
from _cstats import ComputeStats
from _paralog_pi import paralog_pi, ParalogPi
from _filter import (Filter, filter_default, filter_dna, filter_rna,
                             filter_strict, filter_ssr, filter_amino,
                             filter_codon, filter_num)

