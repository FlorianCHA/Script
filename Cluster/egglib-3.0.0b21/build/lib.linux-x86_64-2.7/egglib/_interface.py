"""
This module aims to provide a thin interface to :mod:`_eggwrapper`, the
wrapped C++ library, that is with a little extra functionality added but
with secured accessors.
"""

__license__ = """
    Copyright 2008-2011,2013-2016 Stephane De Mita, Mathieu Siol

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

import re, operator
import _eggwrapper, misc, tools

########################################################################

_random_pool = misc.Pool(tools.Random, (), None)
_intersperse_pool = misc.Pool(_eggwrapper.IntersperseAlign, (), None)
_fasta_formatter_pool = misc.Pool(_eggwrapper.FastaFormatter, (), _eggwrapper.FastaFormatter.close)

CONSENSUS_VALID = set('ACGTURYSWKMBDHVN-acgturyswkmbdhvn')
CONSENSUS_MAPPING = {
     ('B', 'N'): 'N', ('K', 'd'): 'D', ('S', 'W'): 'N', ('y', 'b'): 'B',
     ('n', 'C'): 'N', ('G', 'G'): 'G', ('w', 'M'): 'H', ('N', 'h'): 'N',
     ('c', 'u'): 'Y', ('G', 'b'): 'B', ('d', 'U'): 'D', ('c', 'g'): 'S',
     ('K', 'u'): 'K', ('r', 'w'): 'D', ('a', 'H'): 'H', ('h', 'D'): 'N',
     ('m', 'y'): 'H', ('k', 'b'): 'B', ('n', 'H'): 'N', ('A', 'w'): 'W',
     ('Y', 'r'): 'N', ('W', 'v'): 'N', ('s', 'A'): 'V', ('A', 'N'): 'N',
     ('A', 'Y'): 'H', ('?', 's'): '?', ('k', 'N'): 'N', ('?', 'K'): '?',
     ('V', 'N'): 'N', ('Y', 'V'): 'N', ('T', 'g'): 'K', ('c', 'R'): 'V',
     ('T', 'r'): 'D', ('r', '?'): '?', ('a', 'N'): 'N', ('M', 'u'): 'H',
     ('u', 'w'): 'W', ('H', 'B'): 'N', ('u', 'A'): 'W', ('-', 'd'): '?',
     ('k', 'y'): 'B', ('M', '?'): '?', ('s', 'H'): 'N', ('m', 'k'): 'N',
     ('K', 'V'): 'N', ('C', 'Y'): 'Y', ('H', 'y'): 'H', ('h', 'V'): 'N',
     ('Y', 't'): 'Y', ('-', 'v'): '?', ('G', 'D'): 'D', ('r', 'T'): 'D',
     ('u', 'v'): 'N', ('T', 'N'): 'N', ('W', 'W'): 'W', ('m', 'a'): 'M',
     ('K', 'r'): 'D', ('w', 'd'): 'D', ('?', 'W'): '?', ('h', 'S'): 'N',
     ('t', 'k'): 'K', ('b', 'R'): 'N', ('S', 'S'): 'S', ('N', 'r'): 'N',
     ('v', 'Y'): 'N', ('M', 'w'): 'H', ('K', 'C'): 'B', ('g', 'R'): 'R',
     ('?', 'c'): '?', ('n', 'Y'): 'N', ('U', 'y'): 'Y', ('t', 'W'): 'W',
     ('n', 'a'): 'N', ('a', 'c'): 'M', ('?', 'V'): '?', ('a', 's'): 'V',
     ('-', 'Y'): '?', ('?', 'R'): '?', ('b', 'S'): 'B', ('r', 'M'): 'V',
     ('A', 'K'): 'D', ('Y', '-'): '?', ('D', 'w'): 'D', ('A', 'B'): 'N',
     ('t', 'b'): 'B', ('m', 'h'): 'H', ('y', 'n'): 'N', ('s', 'M'): 'V',
     ('W', 'U'): 'W', ('y', 'm'): 'H', ('m', '?'): '?', ('m', 'S'): 'V',
     ('U', 'T'): 'T', ('Y', 'g'): 'B', ('n', 't'): 'N', ('D', 'S'): 'N',
     ('C', 'k'): 'B', ('C', 'V'): 'V', ('v', 'B'): 'N', ('V', 'T'): 'N',
     ('h', '?'): '?', ('g', 'Y'): 'B', ('w', 'G'): 'D', ('R', 'a'): 'R',
     ('u', 's'): 'B', ('k', 'w'): 'D', ('Y', 'K'): 'B', ('W', 'T'): 'W',
     ('D', '?'): '?', ('t', 'U'): 'T', ('t', 'v'): 'N', ('d', 'D'): 'D',
     ('N', 'c'): 'N', ('a', 'u'): 'W', ('s', 'D'): 'N', ('b', 'U'): 'B',
     ('C', 'U'): 'Y', ('g', 't'): 'K', ('k', 'S'): 'B', ('D', 'u'): 'D',
     ('R', 'b'): 'N', ('u', 'r'): 'D', ('c', 'C'): 'C', ('W', 'k'): 'D',
     ('m', 'v'): 'V', ('S', 'a'): 'V', ('U', 'V'): 'N', ('d', 'N'): 'N',
     ('y', 's'): 'B', ('w', 'C'): 'H', ('u', 'H'): 'H', ('b', 'V'): 'N',
     ('v', 'd'): 'N', ('n', 'K'): 'N', ('s', 'u'): 'B', ('v', 'U'): 'N',
     ('M', 's'): 'V', ('R', 'c'): 'V', ('g', 'V'): 'V', ('U', 'r'): 'D',
     ('m', 'W'): 'H', ('d', 'B'): 'N', ('a', 'w'): 'W', ('b', 'W'): 'N',
     ('-', 'U'): '?', ('A', 't'): 'W', ('C', 'h'): 'H', ('s', 'y'): 'B',
     ('V', 'W'): 'N', ('T', 'C'): 'Y', ('d', 'y'): 'N', ('R', 'd'): 'D',
     ('T', 'H'): 'H', ('s', 'g'): 'S', ('V', 'a'): 'V', ('S', 'h'): 'N',
     ('a', 'n'): 'N', ('h', 'M'): 'H', ('C', 'g'): 'S', ('C', 'R'): 'V',
     ('T', 't'): 'T', ('H', 'T'): 'H', ('c', 't'): 'Y', ('M', 'm'): 'M',
     ('w', 'R'): 'D', ('W', 'h'): 'H', ('s', 'R'): 'V', ('B', 'r'): 'N',
     ('n', 'b'): 'N', ('r', '-'): '?', ('a', 'B'): 'N', ('t', 'R'): 'D',
     ('v', 'w'): 'N', ('b', 'Y'): 'B', ('m', 'T'): 'H', ('H', 'a'): 'H',
     ('-', '-'): '-', ('N', 'y'): 'N', ('T', 'A'): 'W', ('c', 'K'): 'B',
     ('T', 'V'): 'N', ('k', 'M'): 'N', ('c', 'n'): 'N', ('B', 's'): 'B',
     ('R', '-'): '?', ('Y', 'm'): 'H', ('a', 't'): 'W', ('T', 'm'): 'H',
     ('t', 's'): 'B', ('S', 'K'): 'B', ('M', 't'): 'H', ('b', 'w'): 'N',
     ('g', 'S'): 'S', ('w', 'Y'): 'H', ('R', 'g'): 'R', ('K', 'K'): 'K',
     ('w', 'U'): 'W', ('N', 'T'): 'N', ('y', 'd'): 'N', ('y', 'B'): 'B',
     ('B', 't'): 'B', ('m', 'm'): 'M', ('a', 'D'): 'D', ('G', 'y'): 'B',
     ('A', 'H'): 'H', ('T', 'k'): 'K', ('C', 'd'): 'N', ('r', 'U'): 'D',
     ('A', 'C'): 'M', ('V', 'S'): 'V', ('y', 'U'): 'Y', ('S', 'v'): 'V',
     ('k', 'D'): 'D', ('R', 'h'): 'N', ('D', 'd'): 'D', ('k', 's'): 'B',
     ('A', 'U'): 'W', ('?', 'b'): '?', ('W', 'm'): 'H', ('S', 'd'): 'N',
     ('y', 'H'): 'H', ('B', 'u'): 'B', ('m', 'R'): 'V', ('s', 'U'): 'B',
     ('V', 'B'): 'N', ('y', 'y'): 'Y', ('c', 'b'): 'B', ('U', 'g'): 'K',
     ('C', 'c'): 'C', ('H', 'S'): 'N', ('v', 'S'): 'V', ('w', 'V'): 'N',
     ('Y', 'S'): 'B', ('H', 'N'): 'N', ('s', 'N'): 'N', ('B', 'v'): 'N',
     ('u', 'U'): 'T', ('N', 'k'): 'N', ('h', '-'): '?', ('r', 'N'): 'N',
     ('b', 'H'): 'N', ('c', 'M'): 'M', ('?', 'u'): '?', ('y', 'W'): 'H',
     ('c', 'N'): 'N', ('d', 'c'): 'N', ('K', 'H'): 'N', ('D', 'b'): 'N',
     ('d', 'h'): 'N', ('N', 'W'): 'N', ('W', 'c'): 'H', ('B', 'w'): 'N',
     ('u', 'T'): 'T', ('S', 'Y'): 'B', ('h', 'b'): 'N', ('W', 'N'): 'N',
     ('g', 'H'): 'N', ('D', 'Y'): 'N', ('K', 'Y'): 'B', ('r', 'c'): 'V',
     ('S', 'G'): 'S', ('n', 'S'): 'N', ('g', 'W'): 'D', ('k', 'v'): 'N',
     ('M', 'k'): 'N', ('R', 'k'): 'D', ('b', 't'): 'B', ('g', 'N'): 'N',
     ('H', '-'): '?', ('?', 'a'): '?', ('G', 'R'): 'R', ('t', 'r'): 'D',
     ('t', 'n'): 'N', ('r', 'G'): 'R', ('h', 'c'): 'H', ('M', '-'): '?',
     ('k', 'R'): 'D', ('u', '?'): '?', ('B', '?'): '?', ('A', 'G'): 'R',
     ('Y', 'B'): 'B', ('S', 'r'): 'V', ('g', 'D'): 'D', ('d', 'a'): 'D',
     ('M', 'V'): 'V', ('V', 'y'): 'N', ('W', 'a'): 'W', ('r', 'B'): 'N',
     ('B', 'y'): 'B', ('N', 'm'): 'N', ('v', 'H'): 'N', ('?', 'm'): '?',
     ('T', 'w'): 'W', ('H', 'W'): 'H', ('y', 'r'): 'N', ('V', 'h'): 'N',
     ('g', 'U'): 'K', ('r', 'H'): 'N', ('R', 'm'): 'V', ('u', 'G'): 'K',
     ('v', 't'): 'N', ('u', 'S'): 'B', ('H', 'r'): 'N', ('t', 't'): 'T',
     ('-', 'T'): '?', ('G', 's'): 'S', ('-', 'S'): '?', ('b', 'A'): 'N',
     ('Y', 'D'): 'N', ('w', 'v'): 'N', ('U', 'A'): 'W', ('T', '?'): '?',
     ('R', 'n'): 'N', ('c', 'W'): 'H', ('a', 'K'): 'D', ('W', 'g'): 'D',
     ('K', 'W'): 'D', ('K', 'b'): 'B', ('S', 'U'): 'B', ('?', 'C'): '?',
     ('Y', 'u'): 'Y', ('W', 'B'): 'N', ('T', 'u'): 'T', ('U', 'm'): 'H',
     ('b', 'B'): 'B', ('S', 'C'): 'S', ('W', '-'): '?', ('g', 'K'): 'K',
     ('M', 'g'): 'V', ('K', 's'): 'B', ('g', 'B'): 'B', ('Y', 'Y'): 'Y',
     ('G', 'V'): 'V', ('m', 'A'): 'M', ('C', 'n'): 'N', ('v', 'r'): 'V',
     ('h', 'g'): 'N', ('n', 'N'): 'N', ('m', 'n'): 'N', ('c', 'G'): 'S',
     ('?', 't'): '?', ('V', 'k'): 'N', ('A', '-'): '?', ('S', 'n'): 'N',
     ('U', 'C'): 'Y', ('A', 'R'): 'R', ('M', 'R'): 'V', ('G', '?'): '?',
     ('V', 'u'): 'N', ('C', 'm'): 'M', ('U', 'd'): 'D', ('Y', 'w'): 'H',
     ('T', 'R'): 'D', ('D', 'C'): 'N', ('t', 'y'): 'Y', ('b', 'D'): 'N',
     ('D', 'H'): 'N', ('V', 'd'): 'N', ('r', 'r'): 'R', ('M', 'a'): 'M',
     ('u', 'C'): 'Y', ('k', 'g'): 'K', ('W', 'd'): 'D', ('H', 'v'): 'N',
     ('t', 'h'): 'H', ('K', 'T'): 'K', ('d', 't'): 'D', ('N', 's'): 'N',
     ('G', 'w'): 'D', ('v', 'm'): 'V', ('r', 'V'): 'V', ('s', 't'): 'B',
     ('h', 'v'): 'N', ('n', '?'): '?', ('H', 'm'): 'H', ('k', 'C'): 'B',
     ('Y', '?'): '?', ('s', '?'): '?', ('d', 'k'): 'D', ('R', 'r'): 'R',
     ('u', 'B'): 'B', ('c', 'S'): 'S', ('h', 't'): 'H', ('?', 'h'): '?',
     ('?', 'D'): '?', ('a', 'C'): 'M', ('w', 'H'): 'H', ('D', 'A'): 'D',
     ('K', 'A'): 'D', ('r', 'k'): 'D', ('S', '?'): '?', ('M', 'h'): 'H',
     ('M', 'c'): 'M', ('R', 's'): 'V', ('N', 't'): 'N', ('h', 'u'): 'H',
     ('w', 'u'): 'W', ('n', 'A'): 'N', ('D', '-'): '?', ('C', '?'): '?',
     ('K', 'm'): 'N', ('d', 'r'): 'D', ('h', 'k'): 'N', ('G', 'u'): 'K',
     ('A', 'D'): 'D', ('b', 'G'): 'B', ('h', 'd'): 'N', ('n', 'h'): 'N',
     ('V', 'g'): 'V', ('T', 'S'): 'B', ('m', 'd'): 'N', ('A', 'V'): 'V',
     ('R', 't'): 'D', ('M', 'N'): 'N', ('a', 'A'): 'A', ('W', 'y'): 'H',
     ('m', 's'): 'V', ('B', 'a'): 'N', ('h', 'n'): 'N', ('V', 'V'): 'V',
     ('G', 't'): 'K', ('r', 'd'): 'D', ('c', 'v'): 'V', ('C', 'w'): 'H',
     ('k', 'u'): 'K', ('H', 'D'): 'N', ('R', 'u'): 'D', ('d', '-'): '?',
     ('Y', 'h'): 'H', ('V', '?'): '?', ('U', '-'): '?', ('s', 'B'): 'B',
     ('B', 'b'): 'B', ('u', 'Y'): 'Y', ('w', 'n'): 'N', ('a', 'R'): 'R',
     ('G', 'k'): 'K', ('-', 'G'): '?', ('D', 'N'): 'N', ('m', 'g'): 'V',
     ('K', 'N'): 'N', ('c', 'A'): 'M', ('t', 'g'): 'K', ('R', 'v'): 'V',
     ('u', 'N'): 'N', ('y', '?'): '?', ('N', 'n'): 'N', ('B', 'c'): 'B',
     ('m', 'u'): 'H', ('S', 'M'): 'V', ('g', 'A'): 'R', ('n', 'U'): 'N',
     ('U', 'u'): 'T', ('t', 'C'): 'Y', ('-', 'r'): '?', ('A', 'm'): 'M',
     ('M', 'd'): 'N', ('g', 'C'): 'S', ('C', '-'): '?', ('R', 'w'): 'D',
     ('r', 'A'): 'R', ('h', 'y'): 'H', ('?', 'g'): '?', ('G', 'N'): 'N',
     ('k', 'h'): 'N', ('B', 'd'): 'N', ('m', 't'): 'H', ('y', 'h'): 'H',
     ('d', 'G'): 'D', ('?', 'N'): '?', ('n', 'V'): 'N', ('c', 'H'): 'H',
     ('-', 'A'): '?', ('h', 'w'): 'H', ('C', 't'): 'Y', ('s', 'V'): 'V',
     ('A', 'S'): 'V', ('V', 'c'): 'V', ('h', 'H'): 'H', ('U', 'K'): 'K',
     ('n', 'D'): 'N', ('u', '-'): '?', ('D', 'T'): 'D', ('n', 'r'): 'N',
     ('S', 'T'): 'B', ('R', '?'): '?', ('V', 'R'): 'V', ('G', 'h'): 'N',
     ('y', 'A'): 'H', ('D', 'K'): 'D', ('U', 'w'): 'W', ('t', 'A'): 'W',
     ('C', 's'): 'S', ('H', 'C'): 'H', ('y', 't'): 'Y', ('s', 'C'): 'S',
     ('H', 'H'): 'H', ('R', 'y'): 'N', ('u', 'K'): 'K', ('y', 'C'): 'Y',
     ('n', 'B'): 'N', ('d', 'w'): 'D', ('t', '-'): '?', ('D', 'v'): 'N',
     ('-', 'H'): '?', ('a', 'V'): 'V', ('-', 'C'): '?', ('-', 'M'): '?',
     ('B', '-'): '?', ('U', 'B'): 'B', ('k', 'K'): 'K', ('y', 'g'): 'B',
     ('D', 'm'): 'N', ('y', 'M'): 'H', ('K', '-'): '?', ('D', 'R'): 'D',
     ('a', 'G'): 'R', ('n', 'G'): 'N', ('W', 's'): 'N', ('B', 'g'): 'B',
     ('a', 'T'): 'W', ('U', 'n'): 'N', ('h', 'R'): 'N', ('T', '-'): '?',
     ('r', 'Y'): 'N', ('r', 's'): 'V', ('A', 'a'): 'A', ('t', 'd'): 'D',
     ('n', 'c'): 'N', ('g', 'G'): 'G', ('Y', 'n'): 'N', ('G', 'B'): 'B',
     ('d', 'u'): 'D', ('B', 'h'): 'N', ('C', 'b'): 'B', ('T', 'D'): 'D',
     ('r', 'W'): 'D', ('a', 'h'): 'H', ('G', 'm'): 'V', ('c', 'D'): 'N',
     ('s', 'k'): 'B', ('k', 'B'): 'B', ('b', 'K'): 'B', ('w', 'b'): 'N',
     ('A', 'W'): 'W', ('Y', 'R'): 'N', ('S', 'b'): 'B', ('-', 'K'): '?',
     ('w', 'w'): 'W', ('a', 'Y'): 'H', ('b', '-'): '?', ('g', '?'): '?',
     ('C', 'a'): 'M', ('w', 'c'): 'H', ('V', 'n'): 'N', ('y', 'V'): 'N',
     ('H', 'G'): 'N', ('W', '?'): '?', ('v', '?'): '?', ('u', 'W'): 'W',
     ('N', 'u'): 'N', ('H', 'b'): 'N', ('T', 'B'): 'B', ('-', 'D'): '?',
     ('k', 'Y'): 'B', ('G', 'c'): 'S', ('v', 'a'): 'V', ('-', '?'): '?',
     ('s', 'h'): 'N', ('?', 'U'): '?', ('K', 'v'): 'N', ('c', 'Y'): 'Y',
     ('Y', 'a'): 'H', ('Y', 'T'): 'Y', ('T', 'Y'): 'Y', ('v', 'N'): 'N',
     ('b', 'N'): 'N', ('u', 'V'): 'N', ('T', 'n'): 'N', ('N', 'v'): 'N',
     ('W', 'w'): 'W', ('B', 'k'): 'B', ('K', 'G'): 'K', ('b', 'c'): 'B',
     ('w', 'D'): 'D', ('?', 'w'): '?', ('W', 'R'): 'D', ('k', 't'): 'K',
     ('?', 'r'): '?', ('t', 'K'): 'K', ('b', 'r'): 'N', ('a', '?'): '?',
     ('V', 'M'): 'V', ('N', 'R'): 'N', ('b', 'T'): 'B', ('v', 'G'): 'V',
     ('K', 'c'): 'B', ('g', 'r'): 'R', ('m', 'D'): 'N', ('h', 'a'): 'H',
     ('d', 'M'): 'N', ('n', 'y'): 'N', ('c', '?'): '?', ('?', 'v'): '?',
     ('G', 'a'): 'R', ('v', 'c'): 'V', ('M', 'Y'): 'H', ('D', 'G'): 'D',
     ('b', 's'): 'B', ('U', 'H'): 'H', ('Y', 'c'): 'Y', ('?', '?'): '?',
     ('D', 'W'): 'D', ('U', 'S'): 'B', ('t', 'm'): 'H', ('m', 'H'): 'H',
     ('M', 'B'): 'N', ('r', 'b'): 'N', ('W', 'u'): 'W', ('B', 'm'): 'N',
     ('H', 'U'): 'H', ('?', 'M'): '?', ('U', 't'): 'T', ('Y', 'G'): 'B',
     ('d', 'S'): 'N', ('C', 'K'): 'B', ('H', 'K'): 'N', ('N', 'g'): 'N',
     ('V', 't'): 'N', ('g', 'y'): 'B', ('u', 'd'): 'D', ('c', 'r'): 'V',
     ('k', 'W'): 'D', ('y', 'K'): 'B', ('W', 't'): 'W', ('s', 'v'): 'V',
     ('B', 'n'): 'N', ('K', 'D'): 'D', ('S', 'w'): 'N', ('d', 'd'): 'D',
     ('N', 'C'): 'N', ('G', 'g'): 'G', ('w', 'm'): 'H', ('s', 'd'): 'N',
     ('w', 'A'): 'W', ('b', 'u'): 'B', ('c', 'U'): 'Y', ('?', 'y'): '?',
     ('D', 'U'): 'D', ('w', 'S'): 'N', ('K', 'U'): 'K', ('u', 'R'): 'D',
     ('c', 'c'): 'C', ('G', '-'): '?', ('v', 'R'): 'V', ('S', 'A'): 'V',
     ('U', 'v'): 'N', ('W', 'V'): 'N', ('A', 'n'): 'N', ('b', 'v'): 'N',
     ('A', 'y'): 'H', ('n', 'k'): 'N', ('?', 'S'): '?', ('k', 'n'): 'N',
     ('u', 'M'): 'H', ('g', 'v'): 'V', ('c', '-'): '?', ('Y', 'v'): 'N',
     ('T', 'G'): 'K', ('w', 't'): 'W', ('d', 'b'): 'N', ('M', 'U'): 'H',
     ('A', 'T'): 'W', ('?', 'T'): '?', ('C', 'H'): 'H', ('V', 'w'): 'N',
     ('T', 'c'): 'Y', ('u', 'a'): 'W', ('T', 'h'): 'H', ('S', 'H'): 'N',
     ('C', 'y'): 'Y', ('H', 'Y'): 'H', ('n', '-'): '?', ('G', 'd'): 'D',
     ('r', 't'): 'D', ('C', 'G'): 'S', ('H', 't'): 'H', ('-', 'k'): '?',
     ('m', '-'): '?', ('K', 'R'): 'D', ('v', 'v'): 'V', ('s', 'r'): 'V',
     ('B', 'R'): 'N', ('S', 's'): 'S', ('k', 'A'): 'D', ('a', 'b'): 'N',
     ('v', 'y'): 'N', ('M', 'W'): 'H', ('b', 'n'): 'N', ('b', 'y'): 'B',
     ('N', 'Y'): 'N', ('T', 'a'): 'W', ('U', 'Y'): 'Y', ('t', 'w'): 'W',
     ('T', 'v'): 'N', ('a', 'S'): 'V', ('-', 'y'): '?', ('b', 'k'): 'B',
     ('B', 'S'): 'B', ('r', 'm'): 'V', ('A', 'k'): 'D', ('Y', 'M'): 'H',
     ('A', 'b'): 'N', ('t', 'S'): 'B', ('D', 'g'): 'D', ('t', 'D'): 'D',
     ('g', 's'): 'S', ('K', 'k'): 'K', ('N', 'A'): 'N', ('n', 'T'): 'N',
     ('D', 's'): 'N', ('U', '?'): '?', ('B', 'T'): 'B', ('C', 'v'): 'V',
     ('v', 'b'): 'N', ('?', 'Y'): '?', ('G', 'Y'): 'B', ('w', 'g'): 'D',
     ('R', 'A'): 'R', ('b', 'M'): 'N', ('C', 'D'): 'N', ('Y', 'k'): 'B',
     ('V', 's'): 'V', ('y', 'u'): 'Y', ('S', 'V'): 'V', ('t', 'u'): 'T',
     ('u', 'm'): 'H', ('D', 'D'): 'D', ('a', 'U'): 'W', ('w', '-'): '?',
     ('S', 'D'): 'N', ('B', 'U'): 'B', ('C', 'u'): 'Y', ('g', 'T'): 'K',
     ('w', 'a'): 'W', ('V', 'b'): 'N', ('D', 'n'): 'N', ('R', 'B'): 'N',
     ('C', 'C'): 'C', ('H', 's'): 'N', ('W', 'K'): 'D', ('m', 'V'): 'V',
     ('s', 'c'): 'S', ('?', 'd'): '?', ('H', 'n'): 'N', ('s', 'n'): 'N',
     ('B', 'V'): 'N', ('N', 'K'): 'N', ('a', '-'): '?', ('y', 'N'): 'N',
     ('M', 'S'): 'V', ('R', 'C'): 'V', ('c', 'm'): 'M', ('U', 'R'): 'D',
     ('y', 'w'): 'H', ('c', 'k'): 'B', ('m', 'w'): 'H', ('K', 'h'): 'N',
     ('D', 'B'): 'N', ('a', 'W'): 'W', ('n', 'W'): 'N', ('-', 'u'): '?',
     ('m', 'B'): 'N', ('B', 'W'): 'N', ('s', 'Y'): 'B', ('h', 'B'): 'N',
     ('W', 'n'): 'N', ('d', 'Y'): 'N', ('R', 'D'): 'D', ('K', 'y'): 'B',
     ('r', 'C'): 'V', ('s', 'G'): 'S', ('V', 'A'): 'V', ('s', '-'): '?',
     ('n', 's'): 'N', ('g', 'w'): 'D', ('k', 'V'): 'N', ('?', 'G'): '?',
     ('g', 'n'): 'N', ('h', 'm'): 'H', ('N', '-'): '?', ('t', 'N'): 'N',
     ('C', 'r'): 'V', ('T', 'T'): 'T', ('h', 'C'): 'H', ('s', 'w'): 'N',
     ('c', 'T'): 'Y', ('M', 'M'): 'M', ('h', 'h'): 'H', ('w', 'r'): 'D',
     ('W', 'H'): 'H', ('S', 'R'): 'V', ('d', 'V'): 'N', ('n', 'M'): 'N',
     ('M', 'v'): 'V', ('-', 'w'): '?', ('B', 'Y'): 'B', ('H', 'A'): 'H',
     ('N', 'D'): 'N', ('y', 'T'): 'Y', ('H', 'w'): 'H', ('k', 'm'): 'N',
     ('g', 'u'): 'K', ('u', 'h'): 'H', ('N', 'B'): 'N', ('T', 'M'): 'H',
     ('S', 'k'): 'B', ('M', 'T'): 'H', ('G', 'S'): 'S', ('w', 'y'): 'H',
     ('R', 'G'): 'R', ('m', 'c'): 'M', ('b', 'a'): 'N', ('N', 'H'): 'N',
     ('y', 'D'): 'N', ('-', 'n'): '?', ('-', 'B'): '?', ('U', 'a'): 'W',
     ('c', 'w'): 'H', ('H', '?'): '?', ('a', 'k'): 'D', ('A', 'h'): 'H',
     ('g', '-'): '?', ('m', 'Y'): 'H', ('K', 'w'): 'D', ('r', 'u'): 'D',
     ('A', 'c'): 'M', ('w', 'T'): 'W', ('Y', 'U'): 'Y', ('W', 'b'): 'N',
     ('r', 'y'): 'N', ('k', 'd'): 'D', ('R', 'H'): 'N', ('b', 'b'): 'B',
     ('n', 'd'): 'N', ('A', 'u'): 'W', ('?', 'B'): '?', ('W', 'M'): 'H',
     ('g', 'k'): 'K', ('v', '-'): '?', ('m', 'r'): 'V', ('g', 'b'): 'B',
     ('-', 's'): '?', ('?', 'k'): '?', ('c', 'B'): 'B', ('U', 'G'): 'K',
     ('C', 'N'): 'N', ('s', 'K'): 'B', ('h', 'G'): 'N', ('n', 'n'): 'N',
     ('v', 's'): 'V', ('v', 'D'): 'N', ('m', 'N'): 'N', ('m', 'M'): 'M',
     ('g', 'a'): 'R', ('Y', 's'): 'B', ('S', 'N'): 'N', ('U', 'c'): 'Y',
     ('u', 'u'): 'T', ('n', 'w'): 'N', ('M', 'r'): 'V', ('a', 'm'): 'M',
     ('r', 'n'): 'N', ('C', 'M'): 'M', ('m', 'K'): 'N', ('Y', 'W'): 'H',
     ('d', 'C'): 'N', ('-', 'V'): '?', ('t', 'Y'): 'Y', ('b', 'd'): 'N',
     ('d', 'H'): 'N', ('N', 'w'): 'N', ('W', 'C'): 'H', ('r', 'R'): 'R',
     ('u', 't'): 'T', ('S', 'y'): 'B', ('k', 'G'): 'K', ('T', 'K'): 'K',
     ('D', 'y'): 'N', ('t', 'H'): 'H', ('K', 't'): 'K', ('S', 'g'): 'S',
     ('N', 'S'): 'N', ('G', 'W'): 'D', ('v', 'M'): 'V', ('M', 'K'): 'N',
     ('R', 'K'): 'D', ('n', 'm'): 'N', ('v', 'n'): 'N', ('G', 'r'): 'R',
     ('?', '-'): '?', ('r', 'g'): 'R', ('c', 's'): 'S', ('h', 'T'): 'H',
     ('N', '?'): '?', ('-', 'm'): '?', ('k', 'r'): 'D', ('b', 'C'): 'B',
     ('b', '?'): '?', ('A', 'g'): 'R', ('Y', 'b'): 'B', ('k', '?'): '?',
     ('d', 'A'): 'D', ('K', 'a'): 'D', ('r', 'K'): 'D', ('t', 'B'): 'B',
     ('V', 'Y'): 'N', ('W', 'A'): 'W', ('s', 'm'): 'V', ('t', '?'): '?',
     ('d', 'v'): 'N', ('w', 'k'): 'D', ('h', 'U'): 'H', ('T', 'W'): 'W',
     ('d', 'm'): 'N', ('t', 'V'): 'N', ('k', '-'): '?', ('V', 'H'): 'N',
     ('h', 'K'): 'N', ('G', 'U'): 'K', ('w', 's'): 'N', ('R', 'M'): 'V',
     ('u', 'g'): 'K', ('b', 'g'): 'B', ('a', 'd'): 'D', ('H', 'R'): 'N',
     ('T', 's'): 'B', ('s', 'W'): 'N', ('-', 't'): '?', ('M', 'n'): 'N',
     ('a', 'a'): 'A', ('s', 'a'): 'V', ('B', 'A'): 'N', ('Y', 'd'): 'N',
     ('V', 'v'): 'V', ('G', 'T'): 'K', ('r', 'D'): 'D', ('R', 'N'): 'N',
     ('C', 'W'): 'H', ('k', 'U'): 'K', ('h', 's'): 'N', ('W', 'G'): 'D',
     ('H', 'd'): 'N', ('a', 'M'): 'M', ('K', 'B'): 'B', ('S', 'u'): 'B',
     ('g', 'd'): 'D', ('Y', 'H'): 'H', ('T', 'U'): 'T', ('U', 'M'): 'H',
     ('s', 'b'): 'B', ('B', 'B'): 'B', ('S', 'c'): 'S', ('?', 'A'): '?',
     ('N', 'b'): 'N', ('a', 'r'): 'R', ('G', 'K'): 'K', ('M', 'G'): 'V',
     ('K', 'S'): 'B', ('K', 'n'): 'N', ('c', 'a'): 'M', ('Y', 'y'): 'Y',
     ('G', 'v'): 'V', ('t', 'G'): 'K', ('N', 'N'): 'N', ('w', '?'): '?',
     ('d', '?'): '?', ('B', 'C'): 'B', ('m', 'U'): 'H', ('S', '-'): '?',
     ('V', 'K'): 'N', ('h', 'N'): 'N', ('n', 'u'): 'N', ('A', 'r'): 'R',
     ('-', 'R'): '?', ('A', 'M'): 'M', ('V', 'U'): 'N', ('g', 'c'): 'S',
     ('U', 'D'): 'D', ('h', 'Y'): 'H', ('y', 'a'): 'H', ('D', 'c'): 'N',
     ('k', 'H'): 'N', ('B', 'D'): 'N', ('s', 'S'): 'S', ('D', 'h'): 'N',
     ('V', 'D'): 'N', ('?', 'n'): '?', ('n', 'v'): 'N', ('c', 'h'): 'H',
     ('M', 'A'): 'M', ('u', 'c'): 'Y', ('C', 'T'): 'Y', ('w', 'N'): 'N',
     ('W', 'D'): 'D', ('H', 'V'): 'N', ('U', 'k'): 'K', ('d', 'T'): 'D',
     ('V', '-'): '?', ('r', 'v'): 'V', ('s', 'T'): 'B', ('-', 'N'): '?',
     ('v', 'W'): 'N', ('H', 'M'): 'H', ('k', 'c'): 'B', ('V', 'r'): 'V',
     ('G', 'H'): 'N', ('v', 'h'): 'N', ('d', 'K'): 'D', ('R', 'R'): 'R',
     ('u', 'b'): 'B', ('C', 'S'): 'S', ('H', 'c'): 'H', ('?', 'H'): '?',
     ('H', 'h'): 'H', ('w', 'h'): 'H', ('v', 'k'): 'N', ('y', 'c'): 'Y',
     ('D', 'a'): 'D', ('D', 'V'): 'N', ('M', 'H'): 'H', ('a', 'v'): 'V',
     ('M', 'C'): 'M', ('R', 'S'): 'V', ('b', 'm'): 'N', ('U', 'b'): 'B',
     ('D', 'M'): 'N', ('w', 'W'): 'W', ('K', 'M'): 'N', ('d', 'R'): 'D',
     ('a', 'g'): 'R', ('n', 'g'): 'N', ('A', 'd'): 'D', ('B', 'G'): 'B',
     ('V', 'G'): 'V', ('A', 'v'): 'V', ('R', 'T'): 'D', ('m', 'C'): 'M',
     ('r', 'S'): 'V', ('A', 'A'): 'A', ('W', 'Y'): 'H', ('v', 'g'): 'V',
     ('g', 'g'): 'G', ('N', 'a'): 'N', ('Y', 'N'): 'N', ('h', 'W'): 'H',
     ('c', 'V'): 'V', ('B', 'H'): 'N', ('C', 'B'): 'B', ('T', 'd'): 'D',
     ('k', 'a'): 'D', ('G', 'M'): 'V', ('c', 'd'): 'N', ('w', 'K'): 'D',
     ('R', 'U'): 'D', ('y', '-'): '?', ('w', 'B'): 'N', ('y', 'R'): 'N',
     ('r', 'h'): 'N', ('S', 'B'): 'B', ('u', 'y'): 'Y', ('a', 'y'): 'H',
     ('-', 'g'): '?', ('v', 'T'): 'N', ('-', 'W'): '?', ('C', 'A'): 'M',
     ('b', 'h'): 'N', ('y', 'v'): 'N', ('R', 'V'): 'V', ('u', 'n'): 'N',
     ('H', 'g'): 'N', ('K', '?'): '?', ('S', 'm'): 'V', ('g', 'M'): 'V',
     ('N', 'U'): 'N', ('U', 'U'): 'T', ('t', 'c'): 'Y', ('g', 'h'): 'N',
     ('T', 'b'): 'B', ('M', 'D'): 'N', ('v', 'V'): 'V', ('G', 'C'): 'S',
     ('v', 'A'): 'V', ('R', 'W'): 'D', ('g', 'm'): 'V', ('r', 'a'): 'R',
     ('c', 'y'): 'Y', ('Y', 'A'): 'H', ('G', 'n'): 'N', ('T', 'y'): 'Y',
     ('y', 'S'): 'B', ('-', 'b'): '?', ('N', 'V'): 'N', ('-', 'a'): '?',
     ('-', 'c'): '?', ('B', 'K'): 'B', ('K', 'g'): 'K', ('A', 's'): 'V',
     ('V', 'C'): 'V', ('W', 'r'): 'D', ('d', 'g'): 'D', ('k', 'T'): 'K',
     ('D', 't'): 'D', ('s', 's'): 'S', ('V', 'm'): 'V', ('n', 'R'): 'N',
     ('S', 't'): 'B', ('m', 'b'): 'N', ('N', 'M'): 'N', ('h', 'A'): 'H',
     ('D', 'k'): 'D', ('U', 'W'): 'W', ('t', 'a'): 'W', ('G', 'A'): 'R',
     ('v', 'C'): 'V', ('M', 'y'): 'H', ('R', 'Y'): 'N', ('u', 'k'): 'K',
     ('U', 'h'): 'H', ('Y', 'C'): 'Y', ('d', 'W'): 'D', ('U', 's'): 'B',
     ('t', 'M'): 'H', ('-', 'h'): '?', ('d', 'n'): 'N', ('M', 'b'): 'N',
     ('t', 'T'): 'T', ('y', 'Y'): 'Y', ('v', 'u'): 'N', ('m', 'G'): 'V',
     ('N', 'd'): 'N', ('B', 'M'): 'N', ('H', 'u'): 'H', ('k', 'k'): 'K',
     ('y', 'G'): 'B', ('d', 's'): 'N', ('D', 'r'): 'D', ('H', 'k'): 'N',
     ('v', 'K'): 'N', ('N', 'G'): 'N', ('W', 'S'): 'N', ('u', 'D'): 'D',
     ('A', '?'): '?', ('U', 'N'): 'N', ('h', 'r'): 'N', ('y', 'k'): 'B'
}

########################################################################

class SampleView(object):
    
    """
    This class manages the name, sequence and groups of an item of the
    :class:`.Align` and :class:`.Container` classes.
    :class:`~.SampleView` objects allow iteration and general
    manipulation of (large) data sets without unnecessary extraction of
    full sequences. Modifications of :class:`~.SampleView` objects are
    immediately applied to the underlying data holder object.
    :class:`~.SampleView` objects are iterable and allow indexing (the
    values are: name, sequence, group, in that order).

    In principle, only :class:`.Align` and :class:`.Container` instances
    are supposed to build :class:`~.SampleView` instances.

    :param parent: a :class:`.Align` or :class:`.Container` instance.
    :param index: an index within the parent instance.
    :param outgroup: a boolean indicating whether the instance
        represents an outgroup (rather than ingroup) sample.

    .. versionadded:: 2.0.1

    .. versionchanged:: 3.0.0
       Renamed from :class:`.SequenceItem`. No more index check at
       construction. No more string formatting. Supports multiple
       groups. Sequence managed by :class:`.SequenceView` and groups by
       :class:`.GroupView`.
    """
    
    ####################################################################

    def __init__(self, parent, index, outgroup):

        self._parent = parent
        self._index = index
        self._outgroup = bool(outgroup)
        self._sequence = SequenceView(self._parent, self._index, self._outgroup)
        self._group = GroupView(self._parent, self._index, self._outgroup)

    ####################################################################

    @property
    def outgroup(self):

        """
        ``True`` if the sample is part of the outgroup. This value
        cannot be modifed.
        """

        return self._outgroup

    ####################################################################

    @property
    def ls(self):

        """
        Length of the sequence for this sample.
        """

        if self._parent._is_matrix:
            return self._parent._obj.get_nsit()
        elif self._outgroup:
            return self._parent.lo(self._index)
        else:
            return self._parent.ls(self._index)

    ####################################################################

    def __iter__(self):

        """
        Iteration over :class:`~.SampleView` instances yields
        successively the name (as a string), the sequence (as a
        :class:`~.SequenceView` instance), and the group labels (as a
        :class:`~.GroupView` instance).
        """

        for i in (self.name, self.sequence, self.group):
            yield i

    ####################################################################

    def __len__(self):

        """
        The length of :class:`~.SampleView` instances is always 3.
        """

        return 3

    ####################################################################

    def __getitem__(self, index):

        """
        The name (as a string), the sequence (as a
        :class:`~.SequenceView` instance), and the group labels (as a
        :class:`~.GroupView` instance) can be accessed by index (in that
        order).

        :param index: the only supported values are 0, 1 and 2.
        """

        if index==0: return self.name
        elif index==1: return self.sequence
        elif index==2: return self.group
        else: raise IndexError, 'SampleView instances contain only three item'

    ####################################################################

    def __setitem__(self, index, value):

        """
        Only the name and sequence can be set directly.

        :param index: the only supported values are 0 and 1.
        :param value: a value appropriate for setting the name (index 0)
            or the sequence (index 1).
        """

        if index==0: self.name = value
        elif index==1: self._parent._set_sequence(self._index, value, self._outgroup)
        else: raise IndexError, 'invalid index; only 0 and 1 are supported.'

    ####################################################################

    @property
    def index(self):

        """
        Index of the sample in the parent  (:class:`.Align` or
        :class:`.Container` instance containing this sample).
        """

        return self._index

    ####################################################################

    @property
    def parent(self):

        """
        Reference of the parent instance (:class:`.Align` or
        :class:`.Container` instance containing this sample).
        """

        return self._parent

    ####################################################################

    @property
    def name(self):

        """
        Get or set the sample name.
        """

        if self._outgroup: return self._parent._obj.get_name_o(self._index)
        else: return self._parent._obj.get_name_i(self._index)

    @name.setter
    def name(self, value):

        if self._outgroup: self._parent._obj.set_name_o(self._index, value)
        else: self._parent._obj.set_name_i(self._index, value)

    ####################################################################

    @property
    def sequence(self):

        """
        Access to data entries. This attribute is represented by a
        :class:`~.SequenceView` instance which allow modifying the
        underlying sequence container object. It is possible to set this
        attribute using either a string, a list of integers (it can be
        any iterable with a length and it may mix single-character
        strings and integers to represent individual data entries). When
        modifying an :class:`.Align` instance, all sequences must have
        the same length as the current alignment length.
        """

        return self._sequence

    @sequence.setter
    def sequence(self, value):

        self._parent._set_sequence(self._index, value, self._outgroup)

    ####################################################################

    @property
    def group(self):

        """
        Access to group labels. Returns a :class:`~.GroupView` instance
        that can be modified. The values used for assignment may be a
        :class:`~.GroupView`, a list of integers or a single integer. If
        there are less values defined than the number of group levels,
        the values after the last one are left unchanged.
        """

        return self._group

    @group.setter
    def group(self, value):

        if isinstance(value, int):
            value = [value]

        if self._outgroup:
            if len(value) != 1:
                raise ValueError, 'outgroup sample can have only one group label'
            self._parent.set_label_o(self._index, value[0])
            
        else:
            for i, v in enumerate(value):
                self._parent.set_label(self._index, i, v)

########################################################################

class SequenceView(object):

    """
    This class manages the sequence of an item of the :class:`.Align`
    and :class:`.Container` classes. Supports iteration and random
    access (including with slices) which is the preferred way if
    sequences are large because it prevents extraction of the full
    sequence.

    In principle, only :class:`~.SampleView`, :class:`.Align` and
    :class:`.Container` instances are supposed to build
    :class:`~.SequenceView` instances.

    :param parent: a :class:`.Align` or :class:`.Container` instance.
    :param index: an index within the parent instance.
    :param outgroup: a boolean indicating whether the instance
        represents an outgroup (rather than ingroup) sample.

    .. versionadded:: 3.0.0
    """

    ####################################################################

    def __init__(self, parent, index, outgroup):

        self._parent = parent
        self._index = index
        self._outgroup = outgroup

    ####################################################################

    def __iter__(self):

        """
        Iterate over data entries. The returned values are signed
        integers.
        """

        if self._outgroup:
            for i in xrange(len(self)):
                yield self._parent.get_o(self._index, i)
        else:
            for i in xrange(len(self)):
                yield self._parent.get_i(self._index, i)

    ####################################################################

    def string(self):

        """
        Generate a string from all data entries. All allele values must
        be >= 0. Even if this condition is met, the generated string
        might be non printable.
        """

        try:
            return ''.join(map(chr, self))
        except ValueError:
            raise ValueError, 'cannot convert allele value to string'

    ####################################################################
#
#    def __str__(self):
#
#        """
#        Equivalent to the :meth:`~.string` method.
#        """
#
#        return self.string()
#
    ####################################################################

    def __getitem__(self, index):

        """
        Get a given data entry or a range of data entries (if the passed
        index is a :class:`slice`). Returned values are signed integers.

        :param index: index of the value to access, or :class:`slice`.
        """

        if isinstance(index, slice):
            if self._outgroup:
                return [self._parent.get_o(self._index, i) for i in xrange(* index.indices(len(self)))]
            else:
                return [self._parent.get_i(self._index, i) for i in xrange(* index.indices(len(self)))]

        elif isinstance(index, int):
            if index < 0:
                index = len(self) + index
                if index < 0: raise IndexError, 'data index out of range'
            if self._outgroup:
                return self._parent.get_o(self._index, index)
            else:
                return self._parent.get_i(self._index, index)

        else:
            raise TypeError, 'invalid index type (expects slice or int)'

    ####################################################################

    def __setitem__(self, index, value):

        """
        Set data entries.

        :param index: index of the value to set. :class:`slice` are
            supported but, in the case of an :class:`.Align` instance,
            the slice length must be equal to the length of the *value*
            argument. For a :class:`.Container`, a segment may be
            replaced by a longer or short sequence with way.
        :param value: a signed integer or a one-character string (if
            index is a single value) or a list of integers or a string
            with length matching the number of positions in the slice.
            It is allowed to mix integers and one-character strings in a
            list or use other types provided that they are iterable and
            have a length.
        """

        if isinstance(index, int):

            # all checking is done by set/set_o

            if self._outgroup:
                self._parent.set_o(self._index, index, value)
            else:
                self._parent.set_i(self._index, index, value)

        elif isinstance(index, slice):

            start, stop, step = index.indices(len(self))
            num = stop - start

            # allow a non-None argument only if the number values matches

            if step != 1:
                indices = range(start, stop, step)
                if len(indices) != len(value): raise ValueError, 'the number of values must match when a step is used'

                if self._outgroup:
                    for i,j in zip(indices, value):
                        self._parent.set_o(self._index, i, j)
                else:
                    for i,j in zip(indices, value):
                        self._parent.set_i(self._index, i, j)

                # from now, the step may be ignored
            # conservative set (except those with step)

            elif num == len(value):
                if self._outgroup:
                    for i, j in zip(xrange(start, stop), value):
                        self._parent.set_o(self._index, i, j)
                else:
                    for i, j in zip(xrange(start, stop), value):
                        self._parent.set_i(self._index, i, j)

            # non conservative set
    
            else:
                if self._parent._is_matrix:
                    raise ValueError, 'the length of the sequence cannot be changed -- provided {0} data items (required: {1})'.format(len(value), num)

                # insert

                if len(value) > num:
                    if self._outgroup:
                        self._parent.insert_sites_o(self._index, stop, value[num:])
                        for i,v in enumerate(value[:num]): self._parent.set_o(self._index, start+i, v)
                    else:
                        self._parent.insert_sites(self._index, stop, value[num:])
                        for i,v in enumerate(value[:num]): self._parent.set_i(self._index, start+i, v)

                # delete

                elif len(value) < num:
                    if self._outgroup:
                        self._parent.del_sites_o(self._index, start+len(value), num - len(value))
                        for i,v in enumerate(value): self._parent.set_o(self._index, start+i, v)
                    else:
                        self._parent.del_sites(self._index, start+len(value), num - len(value))
                        for i,v in enumerate(value): self._parent.set_i(self._index, start+i, v)

                else:
                    raise RuntimeError, 'Wait.  What am I doing here?'

        else:
            raise TypeError, 'invalid index type (expects slice or int)'

    ####################################################################

    def __delitem__(self, sites):

        """
        Delete data entries.

        This method is only available for samples belonging to a
        :class:`.Container` instance. For :class:`.Align` instances, or
        for inserting data entries to all ingroup and outgroup samples,
        use the method :meth:`~.DataBase.del_site`.

        :param sites: index of the site or sites to remove (to remove
            several sites, use a slice).
        """

        if self._parent._is_matrix:
            raise ValueError, 'cannot delete sites from an Align\'s sequence'

        if isinstance(sites, int):
            if self._outgroup: self._parent.del_sites_o(self._index, sites, 1)
            else: self._parent.del_sites(self._index, sites, 1)
 
        elif isinstance(sites, slice):
            start, stop, step = sites.indices(len(self))
            if step == 1:
                if self._outgroup: self._parent.del_sites_o(self._index, start, stop-start)
                else: self._parent.del_sites(self._index, start, stop-start)
            else:
                if step < 0:
                    raise ValueError, 'negative step is not supported for deleting' # because we will process sites backward (to avoid site shifting and be more efficient)
                if self._outgroup:
                    for i in range(start, stop, step)[::-1]:
                        self._parent.del_sites_o(self._index, i, 1)
                else: 
                    for i in range(start, stop, step)[::-1]:
                        self._parent.del_sites(self._index, i, 1)

        else:
            raise TypeError, 'invalid index type (expects slice or int)'
 
    ####################################################################

    def __len__(self):

        """
        Get the number of sites for this sample.
        """

        if self._parent._is_matrix:
            return self._parent._obj.get_nsit()
        elif self._outgroup:
            return self._parent.lo(self._index)
        else:
            return self._parent.ls(self._index)

    ####################################################################

    def insert(self, position, values):

        """
        Insert data entries.

        This method is only available for samples belonging to a
        :class:`~.Container` instance. For :class:`~.Align` instances,
        (for which it is possible to insert data entries to all ingroup
        and outgroup samples) use the method :meth:`.insert_columns`.

        :param position: the position at which to insert sites. The new
            sites are inserted before the specified index. Use 0 to add
            sites at the beginning of the sequence, and the current
            number of sites for this sample to add sites at the end. If
            the value is larger than the current number of sites for
            this sample, or if ``None`` is provided, new sites are added
            at the end of the sequence.
        :param values: a list of integers or a string containing the
            data entries to insert in the sequence. It is allowed to mix
            integers and one-character strings in a list or use other
            types provided that they are iterable and have a length
        """

        if self._parent._is_matrix:
            raise ValueError, 'cannot insert sites in sequences from an Align'
        if self._outgroup:
            self._parent.insert_sites_o(self._index, position, values)
        else:
            self._parent.insert_sites(self._index, position, values)

    ####################################################################

    def find(self, motif, start=0, stop=None):

        """
        Locate the first instance of a motif.

        Returns the index of the first exact hit to a given substring.
        The returned value is the position of the first base of the hit.
        Only exact matches are implemented. To use regular expression
        (for example to find degenerated motifs), one should extract the
        string for the sequence and use a tool such as the regular
        expression module (:mod:`re`).

        :param motif: a list of integer, or one-character strings (or
            mixing both) or a string constituting the motif to search.
        :param start: position at which to start searching. The method
            will never return a value smaller than *start*. By default,
            search from the start of the sequence.
        :param stop: position at which to stop search (the motif
            cannot overlap this position). No returned value will be
            larger than stop-*n*. By default, or if *stop* is equal to
            or larger than the length of the sequence, search until the
            end of the sequence.
        """

        if self._outgroup: return self._parent.find_motif_o(self._index, motif, start, stop)
        else: return self._parent.find_motif(self._index, motif, start, stop)

    ####################################################################

    def to_upper(self):

        """
        Converts all allele values of this sample to upper case. More
        specifically, transforms all values lying in the range ``a-z``
        to their equivalent in the range ``A-Z``. All other allele
        values are ignored. The underlying object data is modified and
        this method returns ``None``.
        """

        if self._outgroup and self._index >= self._parent._obj.get_nsam_o():
            raise IndexError, 'sample index out of range'

        if not self._outgroup and self._index >= self._parent._obj.get_nsam_i():
            raise IndexError, 'sample index out of range'

        for i in xrange(len(self)):
            if self._outgroup:
                x = self._parent._obj.get_o(self._index, i)
                if x >= 97 and x <= 122:
                    self._parent._obj.set_o(self._index, i, x-32)
            else:
                x = self._parent._obj.get_i(self._index, i)
                if x >= 97 and x <= 122:
                    self._parent._obj.set_i(self._index, i, x-32)

    ####################################################################

    def to_lower(self):

        """
        Converts all allele values of this sample to lower case. More
        specifically, transforms all values lying in the range ``A-Z``
        to their equivalent in the range ``a-z``. All other allele
        values are ignored. The underlying object data is modified and
        this method returns ``None``.
        """

        if self._outgroup and self._index >= self._parent._obj.get_nsam_o():
            raise IndexError, 'sample index out of range'

        if not self._outgroup and self._index >= self._parent._obj.get_nsam_i():
            raise IndexError, 'sample index out of range'

        for i in xrange(len(self)):
            if self._outgroup:
                x = self._parent._obj.get_o(self._index, i)
                if x >= 65 and x <= 90:
                    self._parent._obj.set_o(self._index, i, x+32)
            else:
                x = self._parent._obj.get_i(self._index, i)
                if x >= 65 and x <= 90:
                    self._parent._obj.set_i(self._index, i, x+32)

    ####################################################################

    def strip(self, values, left=True, right=True):

        """
        Delete leading and/or trailing occurrences of any characters
        given in the *values* argument. The underlying object is
        modified and this method returns ``None``.

        :param values: *values* should be a list of integers but can
            also be a string.
        :param left: A bolean indicating whether left-side characters
            should be stripped.
        :param right: A bolean indicating whether left-side characters
            should be stripped.
        """

        if self._parent._is_matrix:
            raise TypeError, 'cannot strip sequences from an Align'

        ls = len(self) # check index of self
        fget = self._parent._obj.get_o if self._outgroup else self._parent._obj.get_i

        if isinstance(values, basestring):
            values = map(ord, values)
        else:
            values = [ord(i) if isinstance(i, basestring) else i for i in values]

        if left:
            i = 0
            while i < ls and fget(self._index, i) in values: i+= 1
            if i > 0:
                del self[:i]
                if i == ls: return
                ls = len(self)

        if right:
            i = ls - 1
            while i >= 0 and fget(self._index, i) in values: i -= 1
            if i < ls - 1: del self[i+1:]

########################################################################

class GroupView(object):

    """
    This class manages the list of group labels of an item of the
    :class:`.Align` and :class:`.Container` classes. This class can be
    represented a list of unsigned integers, is iterable and allows
    random access (read and write) but not slices.

    In principle, only :class:`~.SampleView`, :class:`.Align` and
    :class:`.Container` instances are supposed to build
    :class:`~.GroupView` instances.

    :param parent: a :class:`.Align` or :class:`.Container` instance.
    :param index: an index within the parent instance.
    :param outgroup: a boolean indicating whether the instance
        represents an outgroup (rather than ingroup) sample.

    .. versionadded:: 3.0.0
    """

    ####################################################################

    def __init__(self, parent, index, outgroup):

        self._parent = parent
        self._index = index
        self._outgroup = outgroup

    ####################################################################

    @property
    def outgroup(self):

        """
        ``True`` if the sample is part of the outgroup.
        """

        return self._outgroup

    ####################################################################

    def __len__(self):

        """
        Get the number of available groups.
        """

        if self._outgroup: return 1
        else: return self._parent.ng

    ####################################################################

    def __iter__(self):

        """
        Iterate over values.
        """

        if not self._outgroup:
            for i in xrange(self._parent.ng):
                yield self._parent.get_label(self._index, i)
        else:
            for i in [self._parent.get_label_o(self._index)]:
                yield i

    ####################################################################

    def __getitem__(self, index):

        """
        Get a value. For an outgroup sample, only the value 0 is
        allowed. Negative indexes are not supported.
        """

        if not self._outgroup:
            return self._parent.get_label(self._index, index)

        else:
            return self._parent.get_label_o(self._index)

    ####################################################################

    def __setitem__(self, index, value):

        """
        Set an item value. For an outgroup sample, only the value 0 is
        available. Negative indexes are not supported.
        """

        if not self._outgroup:
            return self._parent.set_label(self._index, index, value)
        else:
            return self._parent.set_label_o(self._index, value)

########################################################################

class DataBase(object):

    """
    Base class for :class:`.Align` and :class:`.Container`.

    This base class cannot be instanciated. Attempting to do so will
    raise an exception.
    """

    ####################################################################

    @classmethod
    def _create_from_data_holder(cls, obj):

        if isinstance(obj, _eggwrapper.DataHolder):

            new_instance = cls.__new__(cls)

            if cls == Align:
                if not obj.get_is_matrix():
                    ls = None
                    for i in xrange(obj.get_nsam_i()):
                        if ls == None: ls = obj.get_nsit_i(i)
                        elif obj.get_nsit_i(i) != ls: raise ValueError, 'cannot convert non-matrix `DataHolder` to `Align`: sequence lengths don\'t match'
                    for i in xrange(obj.get_nsam_o()):
                        if ls == None: ls = obj.get_nsit_o(i)
                        elif obj.get_nsit_o(i) != ls: raise ValueError, 'cannot convert non-matrix `DataHolder` to `Align`: sequence lengths don\'t match'
                obj.set_is_matrix(True)
                new_instance._is_matrix = True
                new_instance._window = None
            else:
                obj.set_is_matrix(False)
                new_instance._is_matrix = False

            new_instance._ns = obj.get_nsam_i()
            new_instance._ns_o = obj.get_nsam_o()
            new_instance._ng = obj.get_ngroups()
            new_instance._motif = _eggwrapper.VectorInt()
            new_instance._obj = obj

            return new_instance

        else:
            raise TypeError, 'unsupported type: {0}'.format(type(obj))

    ####################################################################

    @classmethod
    def create(cls, obj):

        """
        Create a new instance by copying date from the data container
        passed as *obj*. The object *obj* can be:

        * an :class:`~.Align`,
        * a :class:`~.Container` (but all sequences are required to have
          the same length if the target type is :class:`~.Align`),
        * any iterable type yielding a acceptable items (see the
          documentation for :meth:`~.DataBase.add_samples` for more details).

        .. versionadded:: 2.0.1
        """

        new_instance = cls.__new__(cls)
        new_instance.__init__()

        if isinstance(obj, (Align, Container)):
            new_instance.ng = obj._ng
            for sample in obj.iter_outgroup():
                new_instance.add_outgroup(*sample)
        else:
            ng = 0
            for item in obj:
                if len(item) == 3:
                    if isinstance(item[2], int): ng_i = 1
                    else: ng_i = len(item[2])
                    if ng_i > ng: ng = ng_i
            new_instance.ng = ng

        new_instance.add_samples(obj)

        return new_instance

    ####################################################################

    def __init__(self):

        raise 'cannot create a DataBase instance--use Align or Container'

    ####################################################################

    def to_fasta(self, fname=None, first=0, last=_eggwrapper.MAX,
        mapping=None, groups=False, shift_labels=False,
        include_outgroup=True, linelength=50):

        """
        Export alignment in the fasta format.

        :param fname: Name of the file to export data to. By default, the
            file is created (or overwritten if it already exists). If the
            option *append* is ``True``, data is appended at the end of the
            file (and it must exist). If *fname* is ``None`` (default), no
            file is created and the formatted data is returned as a
            :class:`str`. In the alternative case, nothing is returned.

        :param first: If only part of the sequences should be exported:
            index of the first sequence to export.

        :param last: If only part of the sequences should be exported: index
            of the last sequence to export. If the value is larger than the
            index of the last sequence, all sequences are exported until the
            last (this is the default). If *last*<*first*, no sequences are
            exported.

        :param mapping: A string providing the character mapping. Use the
            specified list of characters to map integer allelic values. If
            a non-empty string is provided, the length of the string must be
            larger than the largest possible allele values. In that case,
            the allele values will be used as indexes in order to determine
            which character from this string must be used for outputting. In
            the case that this method is used with an empty string, the
            mapping will not be used and the allele values will be casted
            directly to characters.

        :param groups: A boolean indicating whether group labels
            should be exported, or ignored.

        :param shift_labels: A boolean indicating whether group labels
            should be incremented of one unit when exporting.

        :param include_outgroup: A boolean indicating whether the outgroup
            should be exported. It will be exported at the end of the
            ingroup, and without discriminating label if *include_labels* is
            ``False``.

        :param linelength: The length of lines for internal breaking of
            sequences.
        """

        fasta_formatter = _fasta_formatter_pool.get()

        if fname != None:
            if fasta_formatter.open_file(fname) == False:
                raise ValueError, 'cannot open {0}'.format(fname)
        else:
            fasta_formatter.to_str()

        if mapping == None: fasta_formatter.set_mapping('')
        else: fasta_formatter.set_mapping(mapping)

        fasta_formatter.set_first(first)
        fasta_formatter.set_last(last)
        fasta_formatter.set_labels(groups)
        fasta_formatter.set_shift_groups(shift_labels)
        fasta_formatter.set_outgroup(include_outgroup)

        if linelength < 1: raise ValueError, 'too small value for `linelength` argument'
        fasta_formatter.set_linelength(linelength)

        fasta_formatter.write(self._obj)

        if fname == None: ret = fasta_formatter.get_str()
        _fasta_formatter_pool.put(fasta_formatter)
        if fname == None: return ret

    ####################################################################

    def clear(self):

        """
        Clear the instance and release all memory. In most cases, it is
        preferable to use the method :meth:`~.DataBase.reset`.
        """

        self._obj.clear(self._is_matrix)
        self._ns = 0
        self._ns_o = 0
        self._ng = 0
        self._motif.clear()

    ####################################################################

    def reset(self):

        """
        Reset the instance.
        """

        self._obj.reset(self._is_matrix)
        self._ns = 0
        self._ns_o = 0
        self._ng = 0

    ####################################################################

    @property
    def is_matrix(self):

        """
        ``True`` if the instance is an :class:`.Align`, and ``False`` if
        it is a :class:`.Container`.
        """

        return self._is_matrix

    ####################################################################

    def __len__(self):

        """
        The :py:func:`len` method is a synonym for the :attr:`.ns`
        attribute (note that the outgroup is not taken into account.
        """

        return self._ns

    ####################################################################

    @property
    def ns(self):

        """
        Current number of samples, not considering the outgroup (cannot
        be modified directly).
        """

        return self._ns

    ####################################################################

    @property
    def no(self):

        """
        Current number of outgroup samples (cannot be modified
        directly).
        """

        return self._ns_o

    ####################################################################

    @property
    def ng(self):

        """
        Number of group levels for ingroup (for outgroup, this number is
        always one). It possible to change the value directly.
        """

        return self._ng

    ####################################################################

    @ng.setter
    def ng(self, value):

        if value < 0:
            raise ValueError, 'cannot set number of group levels to a negative value!'

        old_ng = self._ng
        self._ng = value
        self._obj.set_ngroups(value)
        for i in self:
            for j in xrange(old_ng, self._ng):
                i.group[j] = 0

    ####################################################################

    def add_sample(self, name, data, groups=None):

        """
        Add a sample to the instance.

        :param name: name of the new sample.
        :param data: an iterable (may be a string or a list) containing
            the data to set for the new sample. For an :class:`.Align`
            instance and if the sample is not the first, the number of
            data must fit any previous ones (including outgroup samples,
            if any).
        :param groups: if not None, must be an iterable with an
            unsigned integer for each group levels of the instance (if
            None, all group labels are set to 0), or a single integer
            value if only one level needs to be specified.
        """

        self._add(name, data, groups, False)

    ###################################################################

    def add_outgroup(self, name, data, group=None):

        """
        Add an outgroup sample to the instance.

        :param name: name of the new sample.
        :param data: an iterable (may be a string or a list) containing
            the data to set for the new sample. For an :class:`.Align`
            instance and if the sample is not the first, the number of
            data must fit any previous ones (ingroup or outgroup).
        :param group: if not None, must be an unsigned integer (the
            default value is 0).
        """

        self._add(name, data, group, True)

    ####################################################################

    def _add(self, name, data, group, outgroup):

        # increase number of sites if needed

        if self._is_matrix and self._ns==0 and self._ns_o==0:
            self._obj.set_nsit(len(data))

        # increase number of samples

        if outgroup:
            index = self._ns_o
            self._ns_o += 1
            self._obj.set_nsam_o(self._ns_o)
        else:
            index = self._ns
            self._ns += 1
            self._obj.set_nsam_i(self._ns)

        # set the sample

        self._set_sample(index, name, data, group, outgroup)

    ###################################################################

    def add_samples(self, items):

        """
        Add several samples at the end of the instance.

        :param items: items must be an iterable that have a length (for
            example a list, an :class:`.Align` instance or a
            :class:`.Container` instance. Each item of *items* must be
            of length 2 or 3. The first item must be the sample name
            string, the second item is the data values (as a list of
            signed integers or a single string) and the third item, if
            provided, is the list of group labels (one unsigned integer
            for each level). See the method :meth:`~.Align.add_sample`
            for more details about each added item. If the current
            instance is an :class:`.Align`, all added items must have
            the same length which must be the same as any items
            currently present in the instance. For a
            :class:`.Container`, the items may have different lengths.
            The number of group levels is set by the sample with the
            larger number of group labels, and all samples with less
            labels are completed by e default value (0).

        .. versionadded:: 2.0.1
            Original name is :meth:`.addSequences`.

        .. versionchanged:: 3.0.0
            Renamed as ``add_samples``. All items are added in one shot
            (rather than calling the one-sample add method iteratively).
        """

        # increase the number of sites if needed

        if self._is_matrix and self._ns==0 and self._ns_o==0 and len(items) > 0 and len(items[0])>1:
            self._obj.set_nsit(len(items[0][1]))

        # incease the number of samples

        incr = len(items)
        cur = self._ns
        self._ns += incr
        self._obj.set_nsam_i(self._ns)

        # set values

        for i,v in enumerate(items):
            try:
                if len(v) == 2: self._set_sample(cur+i, v[0], v[1], None, False)
                elif len(v) == 3: self._set_sample(cur+i, v[0], v[1], v[2], False)
                else: raise ValueError, 'invalid number of items for sample'
            except:
                self._ns -= incr
                self._obj.set_nsam_i(self._ns)
                raise

    ####################################################################

    def __iter__(self):

        """
        Iterator: each iteration yields an :class:`~.SampleView`
        instance representing the sample at the corresponding index.
        """

        for i in xrange(self._ns):
            yield SampleView(self, i, False)

    ####################################################################

    def iter_outgroup(self):

        """
        Iterator over the outgroup.
        """

        for i in xrange(self._ns_o):
            yield SampleView(self, i, True)

    ####################################################################

    def iter_samples(self, ingrp, outgrp):

        """
        Iterator over either ingroup samples or outgroup samples, or
        both ingroup and outgroup (in that order) samples. Provided for
        flexibility.

        :param ingrp: include ingroup samples.
        :param outgroup: include outgroup samples.
        """

        if ingrp:
            for i in xrange(self._ns):
                yield SampleView(self, i, False)
        if outgrp:
            for i in xrange(self._ns_o):
                yield SampleView(self, i, True)

    ####################################################################

    def _sample(self, index, outgroup):

        if isinstance(index, slice):
            raise ValueError, 'slices are not supported'

        if index < 0:

            if outgroup:
                index = self._ns_o + index
            else:
                index = self._ns + index

            if index < 0: raise IndexError, 'sample index out of range'

        if outgroup:
            if index >= self._ns_o: raise IndexError, 'sample index out of range'
        else:
            if index >= self._ns: raise IndexError, 'sample index out of range'

        return index

    ####################################################################

    def _site(self, index, outgroup, sample):

        if isinstance(index, slice):
            raise ValueError, 'slices are not supported'

        if index < 0:
            if self._is_matrix: index = self._obj.get_nsit() + index
            elif outgroup: index = self._obj.get_nsit_o(sample) + index
            else: index = self._obj.get_nsit_i(sample) + index
            if index < 0: raise IndexError, 'site index out of range'

        if self._is_matrix:
            if index >= self._obj.get_nsit(): raise IndexError, 'site index out of range'

        else:
            if outgroup:
                if index >= self._obj.get_nsit_o(sample): raise IndexError, 'site index out of range'
            else:
                if index >= self._obj.get_nsit_i(sample): raise IndexError, 'site index out of range'

        return index

    ####################################################################

    def _level(self, index):

        if isinstance(index, slice):
            raise ValueError, 'slices are not supported'
        if index < 0:
            index = self._ng + index
            if index < 0: raise IndexError, 'invalid group index (not enough group levels?)'
        if index >= self._ng: raise IndexError, 'invalid group index (not enough group levels?)'

        return index

    ####################################################################

    def get_name(self, index):

        """
        Get the name of a sample.

        :param index: sample index.
        """

        return self._obj.get_name_i(self._sample(index, False))

    ####################################################################

    def get_name_o(self, index):

        """
        Get the name of an outgroup sample.

        :param index: outgroup sample index.
        """

        return self._obj.get_name_o(self._sample(index, True))

    ####################################################################

    def set_name(self, index, name):

        """
        Set the name of a sample.

        :param index: index of the sample
        :param name: new name value.
        """

        if not isinstance(name, str): raise TypeError, 'only str type supported for name'
        self._obj.set_name_i(self._sample(index, False), name)

    ####################################################################

    def set_name_o(self, index, name):

        """
        Set the name of an outgroup sample.

        :param index: index of the outgroup sample
        :param name: new name value.
        """

        if not isinstance(name, str): raise TypeError, 'only str type supported for name'
        self._obj.set_name_o(self._sample(index, True), name)

    ####################################################################

    def get_sequence(self, index):

        """
        Access to the data entries of a given ingroup sample. Returns a
        :class:`~.SequenceView` instance which allow modifying the
        underlying sequence container object.

        :param index: ingroup sample index.
        """

        return SequenceView(self, self._sample(index, False), False)

    ####################################################################

    def get_sequence_o(self, index):

        """
        Access to the data entries of a given outgroup sample. Returns a
        :class:`~.SequenceView` instance which allow modifying the
        underlying sequence container object.

        :param index: outgroup sample index.
        """

        return SequenceView(self, self._sample(index, True), True)

    ####################################################################

    def _set_sequence(self, index, value, outgroup):

        # check index

        index = self._sample(index, outgroup)

        # check input sequence

        if isinstance(value, basestring):
            value = map(ord, value)
        else:
            value = [ord(i) if isinstance(i, basestring) else i for i in value]
        n = len(value)

        # get current length of the sequence

        if self._is_matrix:
            ls = self._obj.get_nsit()
            if n != ls: raise ValueError, 'cannot change length of a sequence for an Align'
        elif outgroup: ls = self._obj.get_nsit_o(index)
        else: ls = self._obj.get_nsit_i(index)

        # change sequence length as needed

        if n != ls:
            if outgroup: self._obj.set_nsit_o(index, n)
            else: self._obj.set_nsit_i(index, n)

        # set all values

        if outgroup: fset = self._obj.set_o
        else: fset = self._obj.set_i

        for i, v in enumerate(value):
            fset(index, i, v)

    ####################################################################

    def set_sequence(self, index, value):

        """

        Replace all data entries for a given ingroup sample by new
        values.

        :param index: ingroup sample index.
        :param value: can be a :class:`~.SequenceView` instance, a
            string or a list of integers (or any iterable with a length
            which may mix single-character strings and integers). When
            modifying an :class:`.Align` instance, all sequences must
            have the same length as the current alignment length.
        """

        self._set_sequence(index, value, False)

    ####################################################################

    def set_sequence_o(self, index, value):

        """

        Replace all data entries for a given outgroup sample by new
        values.

        :param index: outgroup sample index.
        :param value: can be a :class:`~.SequenceView` instance, a
            string or a list of integers (or any iterable with a length
            which may mix single-character strings and integers). When
            modifying an :class:`.Align` instance, all sequences must
            have the same length as the current alignment length.
        """

        self._set_sequence(index, value, True)

    ####################################################################

    def __getitem__(self, index):

        """
        Get a sample. See :meth:`~.get_sample`.
        """

        return self.get_sample(index)

    ####################################################################

    def get_sample(self, index):

        """
        Get the :class:`~.SampleView` instance corresponding to the
        requested index. The returned object allows to modify the
        underlying data.

        :param index: index of the sample to access.
        """

        return SampleView(self, self._sample(index, False), False)

    ####################################################################

    def get_outgroup(self, index):

        """
        Get the :class:`~.SampleView` instance corresponding to the
        requested index for the outgroup. The returned object allows to
        modify the underlying data.

        :param index: index of the outgroup sample to access.
        """

        return SampleView(self, self._sample(index, True), True)

    ####################################################################

    def __delitem__(self, index):

        """
        Delete an item of the instance.

        :param index: index of the sample to delete.
        """

        self.del_sample(index)

    ####################################################################

    def del_sample(self, index):

        """
        Delete an item of the instance.

        :param index: index of the sample to delete.
        """

        self._obj.del_sample_i(self._sample(index, False))
        self._ns -= 1

    ####################################################################

    def del_outgroup(self, index):

        """
        Delete an item of the instance.

        :param index: index of the sample to delete.
        """

        self._obj.del_sample_o(self._sample(index, True))
        self._ns_o -= 1

    ####################################################################

    def __setitem__(self, index, sample):

        """
        Set the values of a sample. See :meth:`~.set_sample`.
        """

        self.set_sample(index, *sample)

    ####################################################################

    def set_sample(self, index, name, data, groups=None):

        """
        Set the values for the sample corresponding to the requested
        index.

        :param index: index of the sample to access (slices are not
            permitted).
        :param name: new name of the sample.
        :param data: string or list of integers given the new values to
            set. In case of an :class:`.Align`, it is required to pass a
            sequence with length matching the number of sites of the
            instance.
        :param groups: a list of integer label values, or a single
            integer value, to set as group labels. The class will ensure
            that all samples have the same number of group labels,
            padding with 0 as necessary. The default corresponds to an
            empty list.
        """

        self._set_sample(index, name, data, groups, False)

    ####################################################################

    def set_outgroup(self, index, name, data, group=None):

        """
        Set the values for a sample of the outgroup

        :param index: index of the sample to access (slices are not
            permitted).
        :param name: new name of the sample.
        :param data: string or list of integers given the new values to
            set. In case of an :class:`.Align`, it is required to pass a
            sequence with length matching the number of sites of the
            instance.
        :param groups: a positive integer to serve as group label (the
            default value is 0.
        """

        self._set_sample(index, name, data, group, True)

    ####################################################################

    def _set_sample(self, index, name, data, groups, outgroup):

        # check index
        self._sample(index, outgroup)

        # set name
        if outgroup:
            self._obj.set_name_o(index, name)
        else:
            self._obj.set_name_i(index, name)

        # set sequence
        ls = len(data)

        if self._is_matrix:
            # trust that the caller has increased ls if the object is an Align
            if ls != self._obj.get_nsit():
                raise ValueError, 'sequence length must match the alignment length'
            if outgroup:
                for i,v in enumerate(data): self.set_o(index, i, v)
            else:
                for i,v in enumerate(data): self.set_i(index, i, v)
        elif outgroup:
            self._obj.set_nsit_o(index, ls)
            for i,v in enumerate(data): self.set_o(index, i, v)
        else:
            self._obj.set_nsit_i(index, ls)
            for i,v in enumerate(data): self.set_i(index, i, v)

        # set group label
        if outgroup:
            if groups == None:
                    self._obj.set_group_o(index, 0)
            elif isinstance(groups, GroupView):
                self._obj.set_group_o(index, groups[0])
            elif not isinstance(groups, int):
                raise TypeError, 'an integer is expected as group label for an outgroup sample'
            elif groups < 0:
                raise ValueError, 'group labels must be >= 0'
            else:
                self._obj.set_group_o(index, groups)

        else:
            if groups == None:
                for i in xrange(self._ng):
                    self._obj.set_group_i(index, i, 0)
            else:
                if isinstance(groups, int): groups = [groups]  # supports single value
                ng = len(groups)

                # not allowed to add new group levels
                if ng > self._ng:
                    raise ValueError, 'not enough group levels specified'

                # set group levels
                for i, g in enumerate(groups):
                    if g < 0:
                        raise ValueError, 'group labels must be >= 0'
                    else:
                        self._obj.set_group_i(index, i, g)

                # set other (non defined) group levels with 0's
                for i in xrange(ng, self._ng):
                    self._obj.set_group_i(index, i, 0)

    ####################################################################

    def get_i(self, sample, site):

        """
        Get a data entry.

        :param sample: sample index.
        :param site: site index.
        """

        return self._obj.get_i(self._sample(sample, False), self._site(site, False, sample))

    ####################################################################

    def get_o(self, sample, site):

        """
        Get an outgroup data entry.

        :param sample: outgroup sample index.
        :param site: site index.
        """

        sample = self._sample(sample, True)
        site = self._site(site, True, sample)
        return self._obj.get_o(sample, site)

    ####################################################################

    def set_i(self, sample, site, value):

        """
        Set an ingroup data entry. The value must be a signed integer or a
        one-character string. In the latter case, :py:func:`ord` is
        called on the behalf of the user.

        :param sample: sample index.
        :param site: site index.
        :param value: allele value.
        """

        if isinstance(value, basestring): value = ord(value)
        self._obj.set_i(self._sample(sample, False), self._site(site, False, sample), value)

    ####################################################################

    def set_o(self, sample, site, value):

        """
        Set an outgroup data entry. The value must be a signed integer
        or a one-character string. In the latter case, :py:func:`ord` 
        is called on the behalf of the user.

        :param sample: sample index.
        :param site: site index.
        :param value: allele value.
        """

        if isinstance(value, basestring): value = ord(value)
        self._obj.set_o(self._sample(sample, True), self._site(site, True, sample), value)

    ####################################################################

    def get_label(self, sample, level):

        """
        Get a group label.

        :param sample: sample index.
        :param level: level index.
        """

        return self._obj.get_group_i(self._sample(sample, False), self._level(level))

    ####################################################################

    def get_label_o(self, sample):

        """
        Get an outgroup label.

        :param sample: sample index.
        """

        return self._obj.get_group_o(self._sample(sample, True))

    ####################################################################

    def set_label(self, sample, level, value):

        """
        Set a group label.

        :param sample: sample index.
        :param level: level index.
        :param value: new group value (unsigned integer).
        """

        if value < 0: raise ValueError, 'negative group labels are not supported'

        self._obj.set_group_i(self._sample(sample, False), self._level(level), value)

    ####################################################################

    def set_label_o(self, sample, value):

        """
        Get an outgroup label.

        :param sample: sample index.
        :param value: new group value (unsigned integer).
        """

        if value < 0: raise ValueError, 'negative group labels are not supported'

        self._obj.set_group_o(self._sample(sample, True), value)

    ####################################################################

    def reserve(self, nsam=0, nout=0, lnames=0, ngrp=0, nsit=0):

        """
        Pre-allocate memory. This method can be used when the size of
        arrays is known a priori, in order to speed up memory
        allocations. It is not necessary to set all values. Values less
        than 0 are ignored.

        :param nsam: number of samples in the ingroup.
        :param nout: number of samples in the outgroup.
        :param lnames: length of sample names.
        :param ngrp: number of group labels.
        :param nsit: number of sites.
        """

        self._obj.reserve(max(0, nsam), max(0, nout), max(0, lnames),
                           max(0, ngrp), max(0, nsit))

    ####################################################################

    def to_outgroup(self, index, label=0):

        """
        Transfer a sample to the outgroup.

        :param index: index of the sample to move to the outgroup.
        :param label: label to assign to the sample (usually, this label
            assigned the sample to an individual).
        """

        if index < 0:
            index = self.ns + index
            if index < 0: raise IndexError, 'sample index out of range'
        if index >= self.ns: raise IndexError, 'sample index out of range'

        self._obj.to_outgroup(index, label)
        self._ns -= 1
        self._ns_o += 1

    ####################################################################

    def find(self, name, include_outgroup=False, regex=False, multi=False, flags=None, index=False):

        """
        Find a sample by its name.

        :param name: name of sample to identify.
        :param include_outgroup: a boolean indicating whether the
            outgroup should be considered.
        :param regex: a boolean indicating whether the value passed a
            *name* is a regular expression. If so, the string is passed
            as is to the re module (using function
            :py:func:`re.search`).
            Otherwise, only exact matches will be considered.
        :param multi: a boolean indicating whether all hits should be
            returned. If so, a list of :class:`~.SampleView` instances
            is always returned (the list will be empty in case of no
            hits). Otherwise, a single :class:`~.SampleView` instance
            (or its index)  will be returned for the first hit, or
            ``None`` in case of no hits.
        :param flags: list of flags to be passed to :py:func:`re.search`
            (ignored if *regex* is ``False``). For example, when looking
            for samples containing the term "sample", for being case
            insensitive, use the following syntax:
            ``align.find("sample", regex=True, flags=[re.I])``. By
            default (``None``) no further argument is passed.
        :param index: boolean indicating whether the index of the sample
            should be returned. In that case return values for hits are
            :class:`int` (by default, :class:`~.SampleView`) instances.
            Warning: it is not allowed to set both *include_outgroup*
            and *index* to ``True`` as there would be no way to
            distinguish between ingroup and outgroup indexes.

        :return: ``None`` if no hits were found, a (potentially empty)
            list of :class:`~.SampleView` instances or class:`int` if
            *multi* is ``True``, or a single :class:`~.SampleView` or
            :class:`int` otherwise.

        """

        if index and include_outgroup: raise ValueError, 'cannot toggle `include_outgroup` and `index` together'
        if multi: ret = []
        if flags==None: flags = []
        for item in self.iter_samples(True, include_outgroup):
            if ((regex==True and re.search(name, item.name, reduce(operator.or_, flags, 0))) or
                (regex==False and name==item.name)):
                    if multi: ret.append(item.index if index else item)
                    else: return item.index if index else item

        if multi: return ret
        else: return None

    ####################################################################

    def find_motif(self, sample, motif, start=0, stop=None):

        """
        Locate the first instance of a motif for an ingroup sample

        Returns the index of the first exact hit to a given substring.
        The returned value is the position of the first base of the hit.
        Only exact matches are implemented. To use regular expression
        (for example to find degenerated motifs), one should extract the
        string for the sequence and use a tool such as the regular
        expression module (:mod:`re`).

        :param sample: sample index.
        :param motif: a list of integer, or one-character strings (or
            mixing both) or a string constituting the motif to search.
        :param start: position at which to start searching. The method
            will never return a value smaller than *start*. By default,
            search from the start of the sequence.
        :param stop: position at which to stop search (the motif
            cannot overlap this position). No returned value will be
            larger than stop-*n*. By default, or if *stop* is equal to
            or larger than the length of the sequence, search until the
            end of the sequence.
        """

        return self._find_motif(sample, motif, start, stop, False)

    ####################################################################

    def find_motif_o(self, sample, motif, start=0, stop=None):

        """
        Locate the first instance of a motif for an outgroup sample

        Returns the index of the first exact hit to a given substring.
        The returned value is the position of the first base of the hit.
        Only exact matches are implemented. To use regular expression
        (for example to find degenerated motifs), one should extract the
        string for the sequence and use a tool such as the regular
        expression module (:mod:`re`).

        :param sample: outgroup sample index.
        :param motif: a list of integer, or one-character strings (or
            mixing both) or a string constituting the motif to search.
        :param start: position at which to start searching. The method
            will never return a value smaller than *start*. By default,
            search from the start of the sequence.
        :param stop: position at which to stop search (the motif
            cannot overlap this position). No returned value will be
            larger than *stop* minus the length of the motif. By
            default, or if *stop* is equal to or larger than the length
            of the sequence, search until the end of the sequence.
        """

        return self._find_motif(sample, motif, start, stop, True)

    ####################################################################

    def _find_motif(self, sample, motif, start, stop, outgroup):

        # adjust indexes

        sample = self._sample(sample, outgroup)
        start = self._site(start, outgroup, sample)

        # adjust stop position

        if self._is_matrix: ls = self.ls
        elif outgroup: ls = self.lo(sample)
        else: ls = self.ls(sample)
        if stop == None or stop > ls: stop = ls
        else: stop = self._site(stop, outgroup, sample)

        # set motif sequence

        self._motif.set_num_values(len(motif))

        if isinstance(motif, basestring):
            for i,v in enumerate(motif):
                self._motif.set_item(i, ord(v))
        else:
            for i,v in enumerate(motif):
                if isinstance(v, basestring): self._motif.set_item(i, ord(v))
                else: self._motif.set_item(i, v)

        ret = self._obj.find(sample, outgroup, self._motif, start, stop)
        if ret == _eggwrapper.MAX: return None
        else: return ret

    ####################################################################

    def names(self):

        """
        Generate the list of ingroup sample names.
        """

        return [item.name for item in self]

    ####################################################################

    def names_outgroup(self):

        """
        Generate the list of outgroup sample names.
        """

        return [item.name for item in self.iter_outgroup()]

    ####################################################################

    def __contains__(self, name):
        
        """
        Tests whether a sequence name occurs in the instance.
        Corresponds to expressions such as ``name in fasta`` and
        returns a boolean. This operation is equivalent the following
        statement ``find(name) != None``.
        """

        return self.find(name) != None

    ####################################################################

    def name_mapping(self):

        """
        Generates a dictionary mapping names to :class:`~.SampleView`
        instances representing all samples of this instance. This method
        is most useful when several sequences have the same name. It may
        be used to detect and process duplicates. It processes the
        ingroup only.
        """

        res = {}
        for i in self:
            if i.name not in res: res[i.name] = []
            res[i.name].append(i)
        return res

    ####################################################################

    def group_mapping(self, level=0, indices=False):

        """
        Generates a dictionary mapping group labels to either
        :class:`~.SampleView` instances (by default) or their indexes
        representing all samples of this instance. It can process and
        ingroup label level or the outgroup.

        :param level: index of the level to consider. If ``None``,
            processes the outgroup.
        :param indices: if ``True``, represent samples by their index
            (within the ingroup or the outgroup depending on *level*
            being non-``None`` or ``None``) instead of
            :class:`~SampleView` instances.
        """

        res = {}
        if level == None:
            for i, item in enumerate(self.iter_outgroup()):
                label = item.group[0]
                if label not in res: res[label] = []
                res[label].append(i if indices else item)
        else:
            level = self._level(level)
            for i, item in enumerate(self):
                label = item.group[level]
                if label not in res: res[label] = []
                res[label].append(i if indices else item)
        return res

    ####################################################################

    def remove_duplicates(self):

        """
        Remove all duplicates, based on name exact matching. For all
        pairs of samples with identical name, only the one occurring
        first is conserved. The current instance is modified and this
        method returns ``None``.
        """

        names = set()
        i = 0
        while i < self._ns:
            name = self.get_sample(i).name
            if name in names:
                self.del_sample(i)
            else:
                names.add(name)
                i+=1

    ###################################################################

    _key_init = 'ABCDEDGHIJKLMNOPQRSTUVWXYZ'
    _key_code = _key_init + _key_init.lower() + '0123456789_'

    def encode(self, nbits=10, random=None, include_outgroup=False):

        """
        Renames all sequences using a random mapping of names to unique
        keys of length *nbits*.

        :param nbits: length of the keys (encoded names). This value
            must be >= 4 and <= 63.

        :param random: a :class:`.Random` instance to be used a random
            generator. By default, use a default instance.

        :param include_outgroup: a :class:`bool` indicating whether the
            outgroup samples should also be considered (in this case, a
            single mapping is returned including both ingroup and
            outgroup samples).

        :return: A dictionary mapping all the generated keys to the
            actual sequence names. The keys are case-dependent and
            guaranteed not to start with a number.
        
        The returned mapping can be used to restore the original names
        using :meth:`~.DataBase.rename`. This method is not affected by
        the presence of sequences with identical names in the original
        instance (and :meth:`~.DataBase.rename` will also work properly
        in that case).

        .. versionadded:: 2.0.1

        .. versionchanged:: 3.0.0
            Keys are forced to start with a capital letter. Now take
            a :class:`.Random` instance. Use library's own random
            number generator. Added an option for the outgroup.
        """

        if nbits<4 or nbits>len(self._key_code):
            raise ValueError, 'invalid value for `nbits`'

        if random == None: rnd = _random_pool.get()
        else: rnd = random

        mapping = {}
        for item in self.iter_samples(True, include_outgroup):
            while True:
                key = (self._key_init[rnd.integer(len(self._key_init))]
                    + ''.join([ self._key_code[rnd.integer(len(self._key_code))] for i in xrange(nbits-1) ]))
                if key not in mapping: break
            mapping[key] = item.name
            item.name = key

        if random == None: _random_pool.put(rnd)

        return mapping

    ###################################################################

    def rename(self, mapping, liberal=False, include_outgroup=False):
        
        """
        Rename sequences of the instance using the provided mapping.

        :param mapping: a :class:`dict` providing the mapping of old
            names (as keys) to new names (which may, if needed, contain
            duplicated).

        :param liberal: if this argument is `False` and a name does not
            appear in *mapping*, a :exc:`~exceptions.ValueError` is
            raised. If *liberal* is ``True``, names that don't appear in
            *mapping* are left unchanged.

        :param include_outgroup: if ``True``, consider both ingroup and
            outgroup samples together (they should be provided together
            in the same mapping).

        :return: The number of samples that have been actually renamed,
            overall.

        .. versionadded:: 2.0.1

        .. versionchanged:: 3.0.0
            Added an option for the outgroup. Added return value.
        """

        cnt = 0
        for item in self.iter_samples(True, include_outgroup):
            name = item.name
            if name in mapping:
                item.name = mapping[name]
                cnt += 1
            else:
                if not liberal:
                    raise ValueError, 'cannot rename sequence: {0}'.format(name)
        return cnt

    ###################################################################

    def subset(self, samples, outgroup=None):

        """
        Generate and return a copy of the instance with only a specified
        list of samples. It is possible to select ingroup and/or
        outgroup samples and the sample indexes are not required to be
        consecutive.

        :param samples: a list (or other iterable type with a length) of
            sample indexes giving the list of ingroup samples that must
            be exported to the return value object. If ``None``, do not
            export any ingroup samples.
        :param outgroup: a list (or other iterable type with a length)
            of sample indexes giving the list of outgroup samples that
            must be exported to the return value object. If ``None``, do
            not export any outgroup samples.

        .. versionadded:: 3.0.0
            Established as a method for :class:`.Align` and
            :class:`.Container`.
        """

        if samples == None: samples = []
        if outgroup == None: outgroup = []

        if self._is_matrix:
            ret = Align()
        else:
            ret = Container()

        ret._ns = len(samples)
        ret._obj.set_nsam_i(ret._ns)
        ret._ns_o = len(outgroup)
        ret._obj.set_nsam_o(ret._ns_o)
        if ret._ns > 0:
            ret._ng = self._ng
            ret._obj.set_ngroups(self._ng)

        if self._is_matrix:
            ls = self._obj.get_nsit()
            ret._obj.set_nsit(ls)
        
        for i, v in enumerate(samples):
            ret.set_name(i, self.get_name(v))
            for j in xrange(self._ng):
                ret._obj.set_group_i(i, j, self._obj.get_group_i(v, j))
            if not self._is_matrix:
                ls = self._obj.get_nsit_i(v)
                ret._obj.set_nsit_i(i, ls)
            for j in xrange(ls):
                ret._obj.set_i(i, j, self.get_i(v, j)) # must use Align.get_i() to have checks on v
        for i, v in enumerate(outgroup):
            ret.set_name_o(i, self.get_name_o(v))
            ret._obj.set_group_o(i, self._obj.get_group_o(v))
            if not self._is_matrix:
                ls = self._obj.get_nsit_o(v)
                ret._obj.set_nsit_o(i, ls)
            for j in xrange(ls):
                ret._obj.set_o(i, j, self.get_o(v, j))

        return ret

    ###################################################################

    def subgroup(self, groups, outgroup_all=False):

        """
        Generate and return a copy of the instance with only samples
        from the specified groups (identified by their group labels,
        including outgroup).

        :param groups: an integer or a dictionary providing the labels
            of the groups that are selected. If an integer, it is
            understood as an ingroup label corresponding to the first
            grouping level (if more than one). If a dictionary, the keys
            of this dictionary must be integer corresponding to grouping
            levels (or ``None`` for the outgroup) and the values must be
            lists containing the requested labels. It is not required to
            include all levels and the outgroup in the dictionary. It is
            allowed to specify labels that are actually not represented
            in the data.
        :param outgroup_all: if ``True``, always include all outgroup
            samples in the returned instance, regardless of whether
            their group labels are specified in the first argument.

        .. versionadded:: 3.0.0
        """

        if self._is_matrix: ret = Align()
        else: ret = Container()
        ret.ng = self._ng

        if isinstance(groups, int): groups = {0: [groups]}

        if None in groups:
            outgroup = groups[None]
            del groups[None]
        else:
            outgroup = []

        for k in groups:
            if k < 0 or k >= self._ng: raise ValueError, 'invalid group level: {0}'.format(k)

        for i in self:
            for k, v in groups.iteritems():
                if i.group[k] in v: ret.add_sample(*i)

        for i in self.iter_outgroup():
            if outgroup_all or i.group[0] in outgroup: ret.add_outgroup(*i)

        return ret

    ###################################################################

    def shuffle(self, level=0, random=None):
        
        """
        Shuffle group labels.

        Randomly reassigns group labels. Modifies the current instance
        and returns ``None``. Only the specified level is affected, and
        only the group labels are modified (the order of samples is not
        changed.

        :param level: index of the group level to shuffle. To shuffle
            the outgroup's labels, set *level* to ``None``.

        :param random: A :class:`.Random` instance to be used a random
            generator. By default, use a default instance.

        .. versionchanged:: 3.0.0
           The outgroup is necessarily processed separately from the
           rest of the instance.  Allow pass a :class:`.Random`
           instance. Use library's own random number generator.
        """

        # extract labels

        if level == None:
            labels1 = [self._obj.get_group_o(i) for i in xrange(self._ns_o)]
        else:
            level = self._level(level)
            labels1 = [self._obj.get_group_i(i, level) for i in xrange(self._ns)]

        # shuffle labels

        if random == None: rnd = _random_pool.get()
        else: rnd = random

        n = len(labels1)
        labels2 = []
        while n>1:
            x = rnd.integer(n)
            labels2.append(labels1.pop(x))
            n-=1
        labels2.append(labels1[0])
        del labels1

        if random == None: _random_pool.put(rnd)

        # reaffect labels

        if level == None:
            for i,v in enumerate(labels2):
                self._obj.set_group_o(i, v)
        else:
            for i,v in enumerate(labels2):
                self._obj.set_group_i(i, level, v)

########################################################################

class Align(DataBase):

    """
    Holds a data set with associated sample names and group information.
    The data consists of given numbers of ingroup and outgroup samples,
    each with the same number of sites. There can be any number of group
    levels (but this number must be the same for all samples), meaning
    that samples can be described by several group labels in addition to
    their name. Group labels are not group indices (they do not need to
    be consecutive). There is a separate data set for samples belonging
    to the outgroup. There can be any number of outgroup samples.
    Outgroup samples always have one level of group labels that should
    be used to specify individuals (when appropriate). All data are
    represented by signed integers

    By default, the constructor generates an empty instance (0
    samples and  0 sites).

    By default, samples have empty names and no group levels are
    defined.

    :param num_sam: number of ingroup samples.
    :param num_out: number of outgroup samples.
    :param num_sit: number of sites.
    :param init: initial values for all data entries (may be a
        signed integer or a single-character string; ignored if
        *num_sam* or *num_sit* is 0).

    .. versionadded:: 3.0.0
       Reimplementation of the Align class.
    """

    ####################################################################

    def __init__(self, nsam=0, nout=0, nsit=0, init=0):

        self._ng = 0
        self._is_matrix = True
        self._window = None
        self._motif = _eggwrapper.VectorInt()
        self._obj = _eggwrapper.DataHolder(True)

        if nsam != 0 or nout != 0 or nsit != 0:
            self._ns = nsam
            self._ns_o = nout
            self._obj.set_nsam_i(nsam)
            self._obj.set_nsam_o(nout)
            self._obj.set_nsit(nsit)
            for i in xrange(nsam):
                for j in xrange(nsit):
                    self.set_i(i, j, init)
            for i in xrange(nout):
                for j in xrange(nsit):
                    self.set_o(i, j, init)
                self._obj.set_group_o(i, 0)

        else:
            self._ns = 0
            self._ns_o = 0

    ####################################################################

    @property
    def ls(self):

        """
        Alignment length. This value cannot be set or modified directly.
        It is not possible to get the number of samples of a single
        ingroup or ingroup sample as they are all equal to this value.
        """

        return self._obj.get_nsit()

    ####################################################################

    def del_columns(self, site, num=1):

        """
        Delete full columns for all ingroup and outgroup samples.

        By default (if ``num=1``), remove a single site. If *num* is
        larger than 1, remove a range of sites.

        :param site: index of the (first) site to remove. This site must
            be a valid index.
        :param num: maximal number of sites to remove. The value cannot
            be negative.
        """

        if num < 0: raise ValueError, 'cannot delete a negative number of columns'
        site = self._site(site, False, None)
        self._obj.del_sites(site, site+num)

    ####################################################################

    def insert_columns(self, position, values):

        """
        Insert sites at a given position to an alignment.

        :param position: the position at which to insert sites. Sites
            are inserted *before* the specified position, so the user
            can use 0 to insert sites at the beginning of the sequence.
            To insert sites at the end of the sequence, pass the current
            length of the alignment, or ``None``. If *position* is
            larger than the length of the sequence or ``None``, new
            sites are inserted at the end of the alignment. The position
            might be negative to count from the end. Warning: the
            position -1 means *before* the last position.
        :param values: a list of signed integers, or a string, providing
            data to insert into the instance. The same sequence will be
            inserted for all ingroup and outgroup samples.
        """

        num = len(values)
        if isinstance(values, basestring):
            values = map(ord, values)

        ls = self._obj.get_nsit()
        if position == None or position > ls: position = ls
        if position < 0:
            position = ls + position
            if position < 0: raise IndexError, 'invalid index (negative value out of bound)'

        self._obj.insert_sites(position, num)
        for i in xrange(self._ns):
            for j, v in enumerate(values):
                self.set_i(i, position+j, v)
        for i in xrange(self._ns_o):
            for j, v in enumerate(values):
                self.set_o(i, position+j, v)

    ####################################################################

    def extract(self, *args):

        """
        Extract given positions (or columns) of the alignment and
        returns a new alignment.

        The two possible ways to call this method are:

            ``extract(start, stop)`` to extract a continuous range of
            sites,

        and 

            ``extract(indexes)`` to extract a random list of positions
            (in any order).

        :param start: first position to extract. This position must be
            a valid index for this alignment.
        :param stop: stop position for the range to extract. **This
            position is not extracted.** If this position is equal to or
            smaller than *start*, empty sequences are extracted. If this
            position is equal to or larger than the length of the
            alignment, or if it is equal to ``None``, all positions
            until the end of the alignment are extracted.
        :param indexes: a list (or other iterable type with a length) of
            alignment positions (or column indexes). This list may
            contain repetitions and does not need to be sorted. The
            positions will be extracted in the specified order.

        Keyword arguments are not supported.

        .. versionadded:: 2.0.1
    """

        # initialize the output object

        ret = Align()
        ret._ns = self._ns
        ret._ns_o = self._ns_o
        ret._ng = self._ng
        ret._obj.set_nsam_i(self._ns)
        ret._obj.set_nsam_o(self._ns_o)
        ret._obj.set_ngroups(self._ng)

        # load the names and group labels

        for i in xrange(self._ns):
            ret.set_name(i, self.get_name(i))
            for j in xrange(self._ng):
                ret._obj.set_group_i(i, j, self._obj.get_group_i(i, j))

        for i in xrange(self._ns_o):
            ret.set_name_o(i, self.get_name_o(i))
            ret._obj.set_group_o(i, self._obj.get_group_o(i))

        # load the sequences

        if len(args) == 1:
            ls = len(args[0])
            ret._obj.set_nsit(ls)
            for j,v in enumerate(args[0]):
                for i in xrange(self._ns):
                    ret._obj.set_i(i, j, self.get_i(i, v))  # use self.get_i() to check index
                for i in xrange(self._ns_o):
                    ret._obj.set_o(i, j, self.get_o(i, v))  # use self.get_o() to check index

        elif len(args) == 2:
            start = self._site(args[0], None, None)
            if args[1] == None:
                stop = self._obj.get_nsit()
            elif args[1] < 0:
                stop = self._site(args[1], None, None)
                if stop < 0: raise IndexError, 'stop position is out of range'
            else:
                stop = min(args[1], self._obj.get_nsit())
            if stop < start: stop = start
            ls = stop - start
            ret._obj.set_nsit(ls)
            for i in xrange(self._ns):
                for j in xrange(ls):
                    ret._obj.set_i(i, j, self._obj.get_i(i, start+j))
            for i in xrange(self._ns_o):
                for j in xrange(ls):
                    ret._obj.set_o(i, j, self._obj.get_o(i, start+j))

        # capture type error

        else:
            raise ValueError, 'extract() expects 1 or 2 arguments, got {0}'.format(len(args))

        return ret

    ####################################################################

    def fix_ends(self):
        
        """
        Designed for nucleotide sequence alignments. Replaces all
        leading and trailing occurrence of the alignment gap symbol
        (the numeric equivalent of ``-``) by missing data symbols
        (``?``). Internal alignment gaps (those having at least one
        character other than ``-`` and ``?`` at each side) are left
        unchanged.

        .. versionchanged:: 3.0.0
           Renamed from :meth:`.fix_gap_ends` to :meth:`.fix_ends`.
        """

        ls = self._obj.get_nsit()

        for ns, fget, fset in [ (self._ns,   self._obj.get_i, self._obj.set_i),
                                (self._ns_o, self._obj.get_o, self._obj.set_o) ]:

            for i in xrange(ns):

                # clean left

                j = 0
                k = ls - 1

                while j < ls:
                    if fget(i, j) == 45: fset(i, j, 63)
                    elif fget(i, j) != 63: break
                    j += 1
                else:
                    k = -1 # prevent right clean

                # clean right (unless reached end of previous while)

                while k >= 0:
                    if fget(i, k) == 45: fset(i, k, 63)
                    elif fget(i, k) != 63: break
                    k -= 1

    ####################################################################

    def column(self, index, ingroup=True, outgroup=True):

        """
        Extract the allele values of a site at a given position.

        :param index: the index of a site within the alignment.

        :param ingroup: A boolean indicating whether ingroup samples
            should be extracted.

        :param outgroup: A boolean indicating whether outgroup samples
            should be extracted.

        Returns one or two lists of integers providing the allele for
        all samples. If both *ingroup* and *outgroup* are ``True``, two
        lists are returned, respectively for the ingroup and for the
        outgroup. If only one of *ingroup* and *outgroup* is ``True``,
        a single list is returned with data for the selected group. If
        neither is ``True``, the returned value is ``None`` (no error).
        """

        index = self._site(index, None, None)

        if ingroup and outgroup: return [self._obj.get_i(i, index) for i in xrange(self._ns)], [self._obj.get_o(i, index) for i in xrange(self._ns_o)]
        elif ingroup: return [self._obj.get_i(i, index) for i in xrange(self._ns)]
        elif outgroup: return [self._obj.get_o(i, index) for i in xrange(self._ns_o)]
        else: return None

    ###################################################################

    def nexus(self, prot=False):

        """
        Generates a simple nexus-formatted string. If *prot* is
        ``True``, adds ``datatype=protein`` in the file, allowing it to
        be imported as proteins (but doesn't perform further checking).

        Returns a nexus-formatted string. Note: any spaces and tabs in
        sequence names are replaced by underscores. This nexus
        implementation is minimal but will normally suffice to export
        sequences to programs expecting nexus.

        Note: only the ingroup is exported. The data must be exportable
        as strings.
        """

        string = ['#NEXUS\n']
        string += ['begin data;\n']
        string += ['dimensions ntax={0} nchar={1};\n'.format(self.ns, self.ls)]
        if prot: type = 'prot'
        else: type = 'dna'
        string += ['format datatype={0} interleave=no gap=-;\n'.format(type)]
        string += ['matrix\n']
        L = 0
        for i in self:
            if (len(i.name) > L): L = len(i.name)
        for i in self:
            string += ['%s  %s\n' %(i.name.ljust(L).replace(' ', '_').replace('\t', '_'), i.sequence.string())]
        string += [';\nend;\n']
        return ''.join(string)

    ####################################################################

    def filter(self, ratio, valid='ACGTacgt', ingroup=True, outgroup=True, relative=True):

        """
        Removes the sequences with too few valid sites. This method
        modifies the current instance and returns ``None``.

        :param ratio: limit threshold, expressed as a proportion of
            either the maximum number of valid data over all processed
            samples (if the *relative* argument is ``True``) or the
            alignment length (otherwise).
        :param valid: a string or an interable of one-character strings
            or integers giving the allelic values considered to be valid
            (note that the comparisons are case-dependent).
        :param ingroup: A boolean indicating whether ingroup samples
            must be processed.
        :param outgroup: A boolean indicating whether outgroup samples
            must be processed.

        If the length of the alignment is 0, or if both *ingroup* and
        *outgroup* are ``False``, nothing is done.
        """

        if self.ls == 0 or (ingroup == False and outgroup == False): return

        # corrects the valid list

        if isinstance(valid, basestring):
            valid = map(ord, valid)
        else:
            valid = [ord(i) if isinstance(i, basestring) else i for i in valid]

        # get the threshold

        if ratio < 0.0 or ratio > 1.0:
            raise ValueError, 'invalid value for the argument ratio'

        if not relative:
            limit = ratio * self.ls
        else:
            if ingroup and outgroup:
                limit = 1.0 * ratio * max([0] + [sum([self._obj.get_i(i, j) in valid for j in xrange(self.ls)]) for i in xrange(self._ns)] +
                                                [sum([self._obj.get_o(i, j) in valid for j in xrange(self.ls)]) for i in xrange(self._ns_o)])
            elif ingroup:
                limit = 1.0 * ratio * max([0] + [sum([self._obj.get_i(i, j) in valid for j in xrange(self.ls)]) for i in xrange(self._ns)])
            elif outgroup:
                limit = 1.0 * ratio * max([0] + [sum([self._obj.get_o(i, j) in valid for j in xrange(self.ls)]) for i in xrange(self._ns_o)])
            else:
                limit = None

        # fills the kill list

        if ingroup:
            kill1 = []
            for i in xrange(self._ns):
                if sum([self._obj.get_i(i, j) in valid for j in xrange(self.ls)]) < limit:
                    kill1.append(i)
            kill1.reverse()

        if outgroup:
            kill2 = []
            for i in xrange(self._ns_o):
                if sum([self._obj.get_o(i, j) in valid for j in xrange(self.ls)]) < limit:
                    kill2.append(i)
            kill2.reverse()

        # delete the samples in kill list, last first

        if ingroup:
            for i in kill1:
                self._obj.del_sample_i(i)
                self._ns -= 1

        if outgroup:
            for i in kill2:
                self._obj.del_sample_o(i)
                self._ns_o -= 1

    ####################################################################

    def phyml(self, ingroup=True, outgroup=False, strict=True, dtype=None):

        """
        Returns a phyml-formatted string representing the content of the
        instance. The phyml format is suitable as input data for the
        PhyML and PAML programmes. Raises a
        :exc:`~exceptions.ValueError` if any name of the instance
        contains at least one character in the following list:
        ``()[]{},;`` or a space, tab, newline or linefeed.
        Group information is never exported.

        The sequences must all be convertible into characters.

        :param ingroup: A boolean indicating whether the ingroup must be
            exported.

        :param outgroup: A boolean indicating whether the outgroup must
            be exported.

        :param strict: enforce that all characters within the instances
            are valid (assuming either nucleotide or amino acid sequences,
            see *dtype*). If ``False``, characters are not checked. Also,
            enforce that their is no blank character, round bracket,
            colon, or comma in sequence names. When checking, both ingroup
            and outgroup are always checked.

        :param dtype: one of ``None`` (default), ``nt`` (nucleotides), or ``aa``
            (amino acids). Type of data assumed (only if *strict* is set to
            ``True``). By default, allow for either nucleotides or amino acids
            (but not a combination of both).

        .. versionchanged:: 3.0.0
            Added the *ingroup*, *outgroup*, *strict*, and *dtype* options.
        """

        if strict:
            if dtype != None and dtype != 'aa' and dtype != 'nt':
                raise ValueError, 'invalid data type'
            if (
                (dtype is None and self._obj.valid_phyml_nt() == False and self._obj.valid_phyml_aa() == False)
             or (dtype is 'nt' and self._obj.valid_phyml_nt() == False)
             or (dtype is 'aa' and self._obj.valid_phyml_aa() == False)):
                    raise ValueError, 'cannot perform phyml conversion: invalid character in sequences'
            if not self._obj.valid_phyml_names():
                    raise ValueError, 'cannot perform phyml conversion: invalid character in names or empty name'

        if ingroup and outgroup: ns = self._ns + self._ns_o
        elif ingroup: ns = self._ns
        elif outgroup: ns = self._ns_o
        else: ns = 0
        
        lines = ['{0} {1}'.format(ns, self.ls)]

        if ingroup: lines += ['{0}  {1}'.format(i.name, i.sequence.string()) for i in self]
        if outgroup: lines += ['{0}  {1}'.format(i.name, i.sequence.string()) for i in self.iter_outgroup()]

        return '\n'.join(lines)

    ####################################################################

    def phylip(self, format='I', ingroup=True, outgroup=False):
        
        """
        Returns a phyml-formatted string representing the content of the
        instance. The phyml format is suitable as input data for PhyML
        and PAML software. Raises a :exc:`~exceptions.ValueError` if any
        name of the instance contains at least one character of the
        following list: ``()[]{},;`` or a space, tab, newline or
        linefeed. Group labels are never exported. Sequence names cannot
        be longer than 10 characters. A :exc:`~exceptions.ValueError`
        will be raised if a longer name is met. *format* must be 'I' or
        'S' (case-independent), indicating whether the data should be
        formatted in the sequential (S) or interleaved (I) format (see
        PHYLIP's documentation for definitions). The user is responsible
        of ensuring that all names are unique. If not, the exported file
        my cause subsequent programs to fail.

        The sequences must all be convertible into characters.

        :param ingroup: A boolean indicating whether the ingroup must be
            exported.

        :param outgroup: A boolean indicating whether the outgroup must
            be exported.

        .. versionchanged:: 3.0.0
            Added the *ingroup* and *outgroup* options.
        """
        
        BLOCK = 10
        NBLOCKS = 6

        if ingroup and outgroup: ns = self._ns + self._ns_o
        elif ingroup: ns = self._ns
        elif outgroup: ns = self._ns_o
        else: ns = 0

        ls = self.ls

        format = format.upper()
        if format not in set('IS'): 
            raise ValueError, 'unknown value for option `format`: %s' %str(format)

        if format.upper()=='I':
            lines = ['  {0:d} {1:d} I'.format(ns, ls)]
        else:
            lines = ['  {0:d} {1:d}'.format(ns, ls)]

        c = 0
        ci = 0

        if ingroup:
            for i in self:
                line = []
                name = i.name
                if len(set('(){}[],; \n\r\t').intersection(name)):
                    raise ValueError, 'phylip format conversion error, invalid character in name: %s' %n
                if len(name)>10:
                    raise ValueError, 'phylip format conversion error, this name is too long: %s' %n
                ci = c
                line.append('{0}{1}'.format(name.ljust(10), ''.join(map(chr, i.sequence[:BLOCK-10]))))
                ci += BLOCK-10
                n = 0
                while n < (NBLOCKS-1) and ci<ls:
                    line.append(' {0}'.format(''.join(map(chr, i.sequence[ci:ci+BLOCK]))))
                    ci += BLOCK
                    if format=='I': n += 1
                lines.append(''.join(line))

        if outgroup:
            for i in self.iter_outgroup():
                line = []
                name = i.name
                if len(set('(){}[],; \n\r\t').intersection(name)):
                    raise ValueError, 'phylip format conversion error, invalid character in name: %s' %n
                if len(name)>10:
                    raise ValueError, 'phylip format conversion error, this name is too long: %s' %n
                ci = c
                line.append('{0}{1}'.format(name.ljust(10), ''.join(map(chr, i.sequence[:BLOCK-10]))))
                ci += BLOCK-10
                n = 0
                while n < (NBLOCKS-1) and ci<ls:
                    line.append(' {0}'.format(''.join(map(chr, i.sequence[ci:ci+BLOCK]))))
                    ci += BLOCK
                    if format=='I': n += 1
                lines.append(''.join(line))

        c = ci

        # if sequential, c should be full

        while c < ls:
            lines.append('')
            line = []
            for i in self.iter_samples(ingroup, outgroup):
                ci = c
                n = 0
                while n < NBLOCKS and ci < ls:
                    if n != 0: line.append(' ')
                    line.append('{0}'.format(''.join(map(chr, i.sequence[ci:ci+BLOCK]))))
                    ci += BLOCK
                    n += 1
                lines.append(''.join(line))
                line = []
            c = ci
            
        return '\n'.join(lines)

    ####################################################################

    def slider(self, wwidth, wstep):
        
        """
        Provides a means to perform sliding-windows analysis over the
        alignment. This method returns a generator that can be used as
        in ``for window in align.slider(wwidth, wstep)``, where each
        step *window* of the iteration will be the reference to a
        :class:`.Align` instance of length *wwidth* (or less if not
        enough sequence is available near the end of the alignment).
        Each step moves forward following the value of *wstep*.

        .. versionchanged:: 3.0.0
            The returned :class:`.Align` is actually a reference to the
            same object, which is stored within the instance and
            repeatively returned.
        """

        if self._window == None:
            self._window = Align()

        cache = 0
        for i in xrange(0, self.ls, wstep):

            # avoids redundant windows
            if min(self.ls, i+wwidth) == cache: break

            self._window.reset()
            for seq in self: self._window.add_sample(seq.name, seq.sequence[i:i+wwidth], seq.group)
            for seq in self.iter_outgroup(): self._window.add_outgroup(seq.name, seq.sequence[i:i+wwidth], seq.group)
            yield self._window

            # record the new position
            cache = min(self.ls, i+wwidth)

    ####################################################################

    def random_missing(self, rate, ch='N', valid='ACGTacgt', ingroup=True, outgroup=True, random=None):

        """
        Randomly introduces missing data in the current instance. Random
        positions of the alignment are changed to missing data. Only
        data that are currently non-missing data are considered.

        :param rate: probability that a non-mssing (as defined after the
            *valid* argument) data is turned into missing data. 
        :param ch: missing data character, to be used for all
            replacements (as a single-character string or an integer).
        :param valid: a string or list of integers (or any iterable,
            possibly missing integers and one-character strings) given
            the allele values that may be turned into missing data.
        :param ingroup: A boolean indicating whether the ingroup must be
            processed.
        :param outgroup: A boolean indicating whether the outgroup must
            be processed.
        :param random: a :class:`~.tools.Random` instance. If ``None``,
            the method will use a default instance.

        .. versionchanged:: 2.1.0
           Restricted to :class:`.Align` instances.

        .. versionchanged:: 3.0.0
           Takes the *ch*, *valid*, *ingroup*, *outgroup* and *random*
           arguments. Reimplementation dropping dependency on the
           non-default ``numpy`` package.
        """

        # make random

        if random == None: rnd = _random_pool.get()

        # check that rate is valid

        if rate < 0 or rate > 1:
            raise ValueError, 'invalid value for the argument rate'

        # corrects the valid list

        if isinstance(valid, basestring):
            valid = map(ord, valid)
        else:
            valid = [ord(i) if isinstance(i, basestring) else i for i in valid]

        # process ingroup

        if ingroup:
            for i in self:
                indexes = [j for (j,v) in enumerate(i.sequence) if v in valid]
                n = rnd.binomial(len(indexes), rate)
                for j in xrange(n):
                    x = rnd.integer(len(indexes))
                    pos = indexes.pop(x)
                    i.sequence[pos] = ch

        # process outgroup

        if outgroup:
            for i in self.iter_outgroup():
                indexes = [j for (j,v) in enumerate(i.sequence) if v in valid]
                n = rnd.binomial(len(indexes), rate)
                for j in xrange(n):
                    x = rnd.integer(len(indexes))
                    pos = indexes.pop(x)
                    i.sequence[pos] = ch

        if random == None: _random_pool.put(rnd)

    ####################################################################

    def consensus(self, ingroup=True, outgroup=False, ignore=False):
        
        """
        Generates the consensus of the object, assuming nucleotide
        sequences. The consensus is generated based on standard
        ambiguity (IUPAC) codes. The consensus is returned as a string
        object, of length matching the alignment length. The input
        alignment can contain nucleotide bases (``A``, ``C``, ``G`` and
        ``T``), and all ambiguity codes (``V``, ``H``, ``M``, ``D``,
        ``R``, ``W``, ``B``, ``S``, ``Y``, ``K`` and ``N``). ``N``
        stands for any of ``A``, ``C``, ``G`` and ``T``. ``B`` stands
        for any of ``C``, ``G`` and ``T`` and so on. Case is ignored.
        ``U`` is treated exactly as ``T``. Gaps are represented by ``-``
        and missing data by ``?``. Any other value will result in a
        ValueError. If a site is not variable, the fixed value is
        incorporated in the consensus in all cases.

        :param ingroup: A boolean indicating whether the ingroup must be
            exported.
        :param outgroup: A boolean indicating whether the outgroup must
            be exported.
        :param ignore: a boolean indicating whether missing data should
            be ignored. If missing data are not ignored (by default),
            any missing data (including gaps) cause the whole site to
            have ? as consensus. If *ignore* is ``True``, missing data
            are skipped and the consensus is based on the other
            (non-missing) data. If there is no non-missing data, the
            consensus is ``N``.

        .. versionchanged:: 3.0.0
            Options are added, and implementation is modified.
        """

        # return empty string if no data and check argument what

        if ingroup and outgroup:
            if self._ns + self._ns_o == 0: return ''
        elif ingroup:
            if self._ns == 0: return ''
        elif outgroup:
            if self._ns_o == 0: return ''
        else:
            return ''

        if self.ls == 0: return ''

        # do the thing

        return ''.join([self._consensus_site(i, ingroup, outgroup, ignore) for i in xrange(self.ls)])

    ####################################################################

    def _consensus_site(self, pos, ingroup, outgroup, ignore):

        cur = None
        if ingroup:
            for i in xrange(self._ns):
                a = chr(self._obj.get_i(i, pos))
                if ignore and a not in CONSENSUS_VALID: pass
                elif cur == None: cur = a
                elif cur == a: pass
                elif (cur, a) not in CONSENSUS_MAPPING: raise ValueError, 'unsupported character: {0}'.format(a)
                else: cur = CONSENSUS_MAPPING[(cur, a)]
        if outgroup:
            for i in xrange(self._ns_o):
                a = chr(self._obj.get_o(i, pos))
                if ignore and a not in CONSENSUS_VALID: pass
                elif cur == None: cur = a
                elif cur == a: pass
                elif (cur, a) not in CONSENSUS_MAPPING: raise ValueError, 'unsupported character: {0}'.format(a)
                else: cur = CONSENSUS_MAPPING[(cur, a)]
        if ignore and cur == None: return 'N'
        else: return cur

    ####################################################################

    def intersperse(self, length, positions=None, alleles='A', random=None):

        """
        Insert non-varying sites within the alignment. The current
        object is permanently modified.

        :param length: Desired length of the final alignment. If the
            value is smaller than the original (current) alignment
            length, nothing is done and the alignment is unchanged.

        :param positions: List of positions of sites of this alignment.
            The number of positions must be equal to the number of sites
            of the alignment (before interspersing). The argument value
            must be either a sequence of positive integers or a sequence
            of real numbers comprised between 0 and 1. In either case,
            values must be in increasing order. In the former case, the
            last (maximal) value must be smaller than the desired length
            of the final alignment. In the latter case, values are
            expressed relatively, and they will be converted to integer
            indexes by the method. In that case, if site positioning is
            non-trivial (typically, if conversion of positions to
            integer yield identical position for different conscutive
            sites), it will be resolved randomly. By default (if
            ``None``), sites are placed regularly along the final
            alignment length. If :class:`int` and :class:`float` types
            are mixed, the first occurring type will condition what will
            happen.

        :param alleles: String (or any sequence of one-character
            strings) providing the alleles to be used to fill
            non-varying positions of the resulting alignment. If there
            is more than one allele, the allele will be picked randomly
            for each site, independently for each interted site.

        :param random: A :class:`.Random` instance to be used a random
            generator. By default, use a default instance.
        """

        # escape if length = 0 (nothing will be done, and it can cause internal errors

        if length < 1: return
        intersperse = _intersperse_pool.get()
        intersperse.set_length(length)

        # get random object

        if random == None:
            rnd = _random_pool.get()
            intersperse.set_random(rnd._obj)
        else:
            intersperse.set_random(random._obj)

        # load alignment

        intersperse.load(self._obj)

        # get positions

        ls = self._obj.get_nsit()

        if positions == None:
            if ls == 1: positions = [0.5]
            else: positions = [float(i)/(ls-1) for i in xrange(ls)]

        if len(positions) != ls: raise ValueError, 'lengths of `positions` must be equal to the original alignment length'
        if isinstance(positions[0], int): need_rounding = False
        elif isinstance(positions[0], float): need_rounding = True
        else: raise TypeError, 'argument `positions` only support `int` and `float` items'

        for idx, pos in enumerate(positions)    :
            cache = -1
            if pos<cache: raise ValueError, 'item in `positions` is not larger than the previous one: {0}'.format(pos)
            cache = pos

            if need_rounding:
                if pos<0 or pos>1: raise ValueError, 'invalid value provided in `positions`: {0}'.format(pos)
                intersperse.set_position(idx, pos)

            else:
                if pos<0 or pos>=length: raise ValueError, 'invalid value provided in `positions`: {0}'.format(pos)
                intersperse.set_round_position(idx, pos)

        # get alleles

        intersperse.set_num_alleles(len(alleles))
        for idx, allele in enumerate(alleles):
            intersperse.set_allele(idx, ord(allele))

        # process

        intersperse.intersperse(need_rounding)

        # return objects to pool
        if random == None: _random_pool.put(rnd)
        _intersperse_pool.put(intersperse)

########################################################################

class Container(DataBase):

    """
    Holds a data set with associated sample names and group information.
    The data consists of given numbers of ingroup and outgroup samples,
    and each of those may have a different number of sites. There can be
    any number of group levels (but this number must be the same for all
    samples), meaning that samples can be described by several group
    labels in addition to their name. Group labels are not group indices
    (they do not need to be consecutive). There is a separate data set
    for samples belonging to the outgroup. There can be any number of
    outgroup samples. Outgroup samples always have one level of group
    labels that should be used to specify individuals (when
    appropriate). All data are represented by signed integers

    Default instance is empty (0 samples and 0 sites).

    .. versionadded:: 3.0.0
       Reimplementation of the Container class.
    """

    ####################################################################

    def __init__(self):

        self._ns = 0
        self._ns_o = 0
        self._ng = 0
        self._is_matrix = False
        self._obj = _eggwrapper.DataHolder(False)
        self._motif = _eggwrapper.VectorInt()

    ####################################################################

    def ls(self, index):

        """
        Get the number of sites of an ingroup sample.

        :param index: sample index.
        """

        return self._obj.get_nsit_i(self._sample(index, False))

    ####################################################################

    def lo(self, index):

        """
        Get the number of sites of an outgroup sample.

        :param index: sample index.
        """

        return self._obj.get_nsit_o(self._sample(index, True))

    ####################################################################

    def _del_sites(self, sample, site, num, outgroup):

        if num < 0:
            raise ValueError, 'the number of sites to remove cannot be negative'

        sample = self._sample(sample, outgroup)
        site = self._site(site, outgroup, sample)

        if outgroup:
            self._obj.del_sites_o(sample, site, site+num)
        else:
            self._obj.del_sites_i(sample, site, site+num)

    ####################################################################

    def del_sites(self, sample, site, num=1):

        """
        Delete data entries from an ingroup sample.

        By default (if ``num=1``), remove a single site. If *num* is
        larger than 1, remove a range of sites.

        :param sample: ingroup sample index.
        :param site: index of the (first) site to remove. This site must
            be a valid index.
        :param num: maximal number of sites to remove. The value cannot
            be negative.
        """

        self._del_sites(sample, site, num, False)

    ####################################################################

    def del_sites_o(self, sample, site, num=1):

        """
        Delete data entries from an outgroup sample.

        By default (if just ``num=1``), remove a single site. If *num*
        is larger than 1, remove a range of sites.

        :param sample: outgroup sample index.
        :param site: index of the (first) site to remove. This site must
            be a valid index.
        :param num: maximal number of sites to remove. The value cannot
            be negative.
        """

        self._del_sites(sample, site, num, True)

    ###################################################################

    def _insert_sites(self, sample, position, values, outgroup):

        num = len(values)
        if isinstance(values, basestring):
            values = map(ord, values)
        sample = self._sample(sample, outgroup)
        ls = self.lo(sample) if outgroup else self.ls(sample)
        if position == None or position >= ls:
            position = ls
        else:
            position = self._site(position, outgroup, sample)

        if outgroup:
            self._obj.insert_sites_o(sample, position, num)
            for i,v in enumerate(values):
                self._obj.set_o(sample, position+i, values[i])
                    
        else:
            self._obj.insert_sites_i(sample, position, num)
            for i,v in enumerate(values):
                self._obj.set_i(sample, position+i, values[i])

    ###################################################################

    def insert_sites(self, sample, position, values):

        """
        Insert sites at a given position for a given ingroup sample

        :param sample: index of the ingroup sample to which insert
            sites.
        :param position: the position at which to insert sites. Sites
            are inserted *before* the specified position, so the user
            can use 0 to insert sites at the beginning of the sequence.
            To insert sites at the end of the sequence, pass the current
            length of the sequence. If *position* is larger than the
            length of the sequence or ``None``, new sites are inserted
            at the end of the sequence. The position might be negative
            Warning: the position -1 means *before* the last position.
        :param values: a list of signed integers, or a string, providing
            data to insert into the instance.
        """

        self._insert_sites(sample, position, values, False)

    ###################################################################

    def insert_sites_o(self, sample, position, values):

        """
        Insert sites at a given position for a given outgroup samples

        :param sample: index of the outgroup sample to which insert
            sites.
        :param position: the position at which to insert sites. Sites
            are inserted *before* the specified position, so the user
            can use 0 to insert sites at the beginning of the sequence.
            To insert sites at the end of the sequence, pass the current
            length of the sequence. If *position* is larger than the
            length of the sequence or ``None``, new sites are inserted
            at the end of the sequence. The position might be negative
            Warning: the position -1 means *before* the last position.
        :param values: a list of signed integers, or a string, providing
            data to insert into the instance.
        """

        self._insert_sites(sample, position, values, True)

    ###################################################################

    def equalize(self, value = '?'):

        """
        Extend sequences such as they all have the length of the longest
        sequence (over both ingroup and outgroup).

        :param value: the value to use to extend sequences, as an
            integer or a single-character string.
        """

        ls = 0
        for i in xrange(self._ns):
            if self._obj.get_nsit_i(i) > ls:
                ls = self._obj.get_nsit_i(i)
        for i in xrange(self._ns_o):
            if self._obj.get_nsit_o(i) > ls:
                ls = self._obj.get_nsit_o(i)

        if isinstance(value, basestring):
            value = [ord(value)]
        else:
            value = [value]

        for i in xrange(self._ns):
            if self._obj.get_nsit_i(i) < ls:
                n = ls - self._obj.get_nsit_i(i)
                self.insert_sites(i, None, value * n)

        for i in xrange(self._ns_o):
            if self._obj.get_nsit_o(i) < ls:
                n = ls - self._obj.get_nsit_o(i)
                self.insert_sites_o(i, None, value * n)
