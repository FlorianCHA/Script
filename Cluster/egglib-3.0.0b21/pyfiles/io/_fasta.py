"""
This module contains the fasta parser.
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

from .. import _eggwrapper
from .. import _interface
from .. import misc

########################################################################

_parser_pool = misc.Pool(_eggwrapper.FastaParser, (), _eggwrapper.FastaParser.close, None)
_cont_pool = misc.Pool(_interface.Container, (), _interface.Container.reset, None)

########################################################################

def from_fasta(source, groups=False, string=False, cls=None):

    """
    Create a new instance of either :class:`.Align` or
    :class:`.Container` from data read from the file whose name is
    provided as argument (or, if *string* is set to ``True``, from
    the string passed as first argument).

    :param source: name of a fasta-formatted sequence file. If the
        *string* argument is ``True``, read *source* as a
        fasta-formatted string. If the returned type if
        :class:`.Align`, the sequence are required to be aligned.
    :param groups: boolean indicating whether group labels should be
        imported. If so, they are not actually required to be
        present for each (or any) sequence. If not, the labels are
        ignored and considered to be part of the sequence.
    :param string: boolean indicating whether the first argument is
        an explicit fasta-formatted string (by default, it is taken
        as the name of a fasta-formatted file).
    :param cls: type that should be generated. Possible values are:
        :class:`.Align` (then, data must be aligned),
        :class:`.Container` or ``None``. In the latter case, an
        :class:`!Align` is returned if data are found to be aligned or
        if the data set is empty, and otherwise a :class:`.Container` is
        returned.

    :return: A new :class:`.Container` or :class:`.Align` instance
        depending on the value of the *cls* option.
    """

    fasta_parser = _parser_pool.get()

    if string: fasta_parser.set_string(source)
    else: fasta_parser.open_file(source)
    obj = _eggwrapper.DataHolder(False)
    fasta_parser.read_all(groups, obj)

    if cls is _interface.Align or cls is None:
        ns = set([obj.get_nsit_i(i) for i in xrange(obj.get_nsam_i())] + [obj.get_nsit_o(i) for i in xrange(obj.get_nsam_o())])
        if len(ns) == 0:
            ns = 0
            if cls is None: cls = _interface.Align
        elif len(ns) == 1:
            ns = ns.pop()
            if cls is None: cls = _interface.Align
        else:
            if cls is _interface.Align:
                raise ValueError, 'cannot create `Align` from {0}: lengths of sequences do not match'.format('string' if string else '`'+source+'`')
            cls = _interface.Container

    elif cls is not _interface.Container:
        raise ValueError, 'invalid value provided for `cls`'

    try:
        if cls is _interface.Container: return _interface.Container._create_from_data_holder(obj)
        else: return _interface.Align._create_from_data_holder(obj)
    finally:
        _parser_pool.put(fasta_parser)

########################################################################

class fasta_iter(object):

    """
    Iterative sequence-by-sequence fasta parser. Return an object that
    can be iterated over:

    .. code-block:: python

       for item in egglib.io.fasta_iter(fname):
           do things

    :class:`.fasta_iter` objects support the ``with`` statement:

    .. code-block:: python

       with egglib.io.fasta_iter(fname) as f:
           for item in f:
               do things

    Each iteration yields a :class:`.SampleView` instance (which is
    valid only during the iteration round, see the warning below). It is
    also possible to iterate manually using :meth:`~fasta_iter.next`.
    The number of groups is defined by the current sample (if the number
    of defined groups varies among samples, it is reset at each
    iteration).

    .. warning::
        The aim of this iterator is to iterator over large fasta files
        without actually storing all data in memory at the same time.
        :class:`.SampleView` provided for each iteration are a proxy to
        a local :class:`.Container` instance that is recycled at each
        iteration step. They should be used immediately and never stored
        as this. If one want to store data accessible through any
        :class:`.SampleView`, they should copy this data to another data
        structure (typically using :meth:`.Container.add_sample`).

    :param fname: name of a fasta-formatted file.
    :param groups: if ``True``, import group labels from sequence names
        (by default, they are considered as part of the name).

    .. versionadded:: 3.0.0
    """

    ####################################################################

    def __init__(self, fname, groups=False):
        self._parser = _parser_pool.get()
        self._parser.open_file(fname)
        self._cont = _cont_pool.get()
        self._groups = groups

    ####################################################################

    def __del__(self):
        if self._parser != None: _parser_pool.put(self._parser)
        if self._cont != None: _cont_pool.put(self._cont)

    ####################################################################

    def __enter__(self):
        return self

    ####################################################################

    def __exit__(self, exc_type, exc_val, exc_tb):
        _parser_pool.put(self._parser)
        _cont_pool.put(self._cont)
        self._parser = None
        self._cont = None
        return False

    ####################################################################

    def __iter__(self):
        return self

    ####################################################################

    def next(self):

        """
        Perform an iteration round. Raise a
        :exc:`~exceptions.StopIteration` exception if the file is
        exhausted. The normal usage of this type of objects is with the
        ``for`` statement.
        """

        if not self._parser.good(): raise StopIteration
        self._cont.reset()
        self._parser.read_sequence(self._groups, self._cont._obj)
        self._cont._ns = self._cont._obj.get_nsam_i()
        self._cont._ns_o = self._cont._obj.get_nsam_o()
        self._cont._ng = self._cont._obj.get_ngroups()

        if self._cont._ns == 1: outgroup = False
        else: outgroup = True
        return _interface.SampleView(self._cont, 0, outgroup)
