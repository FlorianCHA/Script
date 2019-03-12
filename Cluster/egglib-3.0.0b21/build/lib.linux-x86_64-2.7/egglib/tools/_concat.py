"""
This module contains the single function :function:`concat` which is
available from the package :mod:`tools`.
"""

__license__ = """
    Copyright 2008-2012,2015-2016 Stephane De Mita, Mathieu Siol

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

from .. import _interface

########################################################################

def concat(*aligns, **kwargs):

    """ egglib.tools.concat(align1, align2, ..., spacer=0, ch='?', group_check=True, no_missing=False, ignore_names=False, dest=None)

    Concatenates sequence alignments provided as :class:`.Align`
    instances passed as arguments to this function. A unique
    :class:`.Align` is produced. All different sequences from all passed
    alignments are represented in the final alignment. Sequences whose
    name match are matching are concatenated. In case several sequences
    have the same name in a given segment, the first one is considered
    and others are discarded. In case a sequence is missing for a
    particular segment, a stretch of non-varying characters is inserted
    to replace the unknown sequence.

    All options (excluding the alignements to be concatenated) must be
    specified as keyword arguments, otherwise they will be treated as
    alignments, which may generate an error.

    :param align1, align2: Two or more :class:`.Align` instances (their
        order is used for concatenation. It is not allowed to specify
        them using the keyword syntax.

    :param spacer: Length  of  unsequenced stretches (represented by
        non-varying characters) between concatenated alignments. If
        *spacer* is a positive integer, the length of all stretches will
        be identical. If *spacer* is an iterable containing integers,
        each specifying the interval between two consecutive alignments
        (if *aligns* contains ``n`` alignments, *spacer* must be of
        length ``n-1``).

    :param ch: Character to used for conserved stretches and for missing
        segments.

    :param group_check: If ``True``, an exception will be raised in case
        of a mismatch between group labels of different sequence
        segments bearing the same name. Otherwise, the group of the
        first segment found will be used as group label of the final
        sequence.

    :param no_missing: If ``True``, an exception will be raised in case
        the list of samples differs between :class:`.Align` instances.
        Then, the number of samples must always be the same and all
        samples must always be present (although it is possible that
        they consist in missing data only). Ignored if *ignore_names* is
        ``True``.

    :param ignore_names: Don't consider sample names and concatenate
        sequences based on they order in the instance. Then the value of
        the option *no_missing* is ignored and the number of samples is
        required to be constant over alignments.

    :param dest: An optional :class:`.Align` instance where to place
        results. This instance is automatically reset, ignoring all
        data previously loaded. If this argument is not ``None``, the
        function returns nothing and the passed instance is modified.
        Allows to recycle the same object in intensive applications.

    :return: If *dest* is ``None``, a new :class:`.Align` instance.
        If *dest* is ``None``, this function returns ``None``.

    .. note::
        Outgroup samples are processed like ingroup samples, but they
        are processed independently. In particular, they may have the
        same name as ingroup samples but will not be skipped.

    .. versionadded:: 2.0.1
        The arguments allowing to customize function's behaviour.

    .. versionchanged:: 3.0.0
        Major interface change: the alignments are not passed as a list
        as before, but as a variable-length list of positional
        arguments, implying that options must be specified as keyword
        arguments only. Besides, added options *no_missing*,
        *ignore_names* and *dest*. Renamed option *groupCheck* as
        *group_check*. Removed option *strict* (now name comparison is
        always strict).
    """

    # import default value of options

    spacer = kwargs.get('spacer', 0)
    ch = kwargs.get('ch', '?')
    group_check = kwargs.get('group_check', True)
    no_missing = kwargs.get('no_missing', False)
    ignore_names = kwargs.get('ignore_names', False)
    dest = kwargs.get('dest', None)
    for key in kwargs:
        if key not in ['spacer', 'ch', 'group_check', 'no_missing',
                       'ignore_names', 'dest']:
            raise ValueError, 'this key is not a valid argument for concat: `{0}`'.format(key)

    # process arguments

    for aln in aligns:
        if not isinstance(aln, _interface.Align):
            raise TypeError, 'all positional arguments are required to be Align instances'

    nloc = len(aligns)
    if isinstance(spacer, int):
        spacer = [spacer] * (nloc -1) # supports nloc==0
        if spacer < 0: raise ValueError, '`spacer` argument must not be negative'
    elif min(spacer) < 0: raise ValueError, '`spacer` argument must not be negative'
    elif len(spacer) == 0 and nloc == 0:
        pass
    elif len(spacer) != nloc-1:
        raise ValueError, '`spacer` argument does not have the right number of items'
    spacer.append(0) # convenience to avoid having an "if" in the main loop
    if len(ch) != 1:
        raise ValueError, '`ch` argument must be a single character'

    # get the total length

    ls = sum([aln.ls for aln in aligns]) + sum(spacer)

    # get the list of samples

    if not ignore_names:

        # get the list of names of each alignment

        names_in = []
        names_out = []
        for aln in aligns:
            names_in.append(aln.names())
            names_out.append(aln.names_outgroup())

        # get the total list of names (as a dict without values)

        samples_in = dict.fromkeys(set().union(*names_in))
        samples_out = dict.fromkeys(set().union(*names_out))

        # numer of samples

        nsi = len(samples_in)
        nso = len(samples_out)

        # check that list is constant if requested

        if no_missing:
            for aln in aligns:
                if aln.ns != nsi or aln.no != nso:
                    raise ValueError, 'the list of samples is required to be constant over alignments (`no_missing` option is ``True``)'

        # get the index of each sample

        for name in samples_in:
            samples_in[name] = []
            for aln, names in zip(aligns, names_in):
                if name in names: samples_in[name].append(names.index(name))
                else: samples_in[name].append(None)

        for name in samples_out:
            samples_out[name] = []
            for aln, names in zip(aligns, names_out):
                if name in names: samples_out[name].append(names.index(name))
                else: samples_out[name].append(None)

    else:
        nsi = set([aln.ns for aln in aligns])
        nso = set([aln.no for aln in aligns])
        if len(nsi) != 1 or len(nso) != 1:
            raise ValueError, 'all alignments are required to have the same number of samples (`ignore_names` option is ``True``)'
        nsi = nsi.pop()
        nso = nso.pop()

        samples_in = dict([(i, [i]*nloc) for i in xrange(nsi)])
        samples_out = dict([(i, [i]*nloc) for i in xrange(nso)])

    # create or reset destination

    if dest == None:
        conc = _interface.Align(nsam=nsi, nout=nso, nsit=0)
        if not ignore_names:
            for i, v in enumerate(samples_in): conc.set_name(i, v)
            for i, v in enumerate(samples_out): conc.set_name_o(i, v)
    else:
        dest.reset()
        conc = dest
        if ignore_names:
            for v in samples_in: conc.add_sample('', [])
            for v in samples_out: conc.add_outgroup('', [])
        else:
            for v in samples_in: conc.add_sample(v, [])
            for v in samples_out: conc.add_outgroup(v, [])

    conc._obj.set_nsit(ls) # this doesn't initialize new values

    # set names

    names_in = sorted(samples_in)
    names_out = sorted(samples_out)

    if not ignore_names:
        for i, name in enumerate(names_in): conc.set_name(i, name)
        for i, name in enumerate(names_out): conc.set_name_o(i, name)

    # process ingroup groups

    group_mapping = {}
    for name in samples_in:
        groups = [aligns[i].get_sample(j).group for (i,j) in enumerate(samples_in[name]) if j != None]
        if len(groups) > 0:
            if group_check:
                for i in xrange(len(groups)-1):
                    if len(groups[0]) != len(groups[i+1]):
                        raise ValueError, 'the group labels are required to be consistent (`group_check` option if ``True``)'
                    for a, b in zip(groups[0], groups[i+1]):
                        if a != b:
                            raise ValueError, 'the group labels are required to be consistent (`group_check` option if ``True``)'
            group_mapping[name] = groups[0]
        else:
            group_mapping[name] = []

    ng = max(map(len, group_mapping.itervalues()))
    conc.ng = ng

    for idx, name in enumerate(names_in):
        grp = conc.get_sample(idx).group
        for i, g in enumerate(group_mapping[name]):
            grp[i] = g

    # process outgroup groups

    for idx, name in enumerate(names_out):
        groups = [aligns[i].get_outgroup(j).group[0] for (i,j) in enumerate(samples_out[name]) if j != None]
        if len(groups) > 0:
            if group_check:
                for i in xrange(len(groups)-1):
                    if groups[0] != groups[i+1]:
                        raise ValueError, 'the group labels are required to be consistent (`group_check` option if ``True``)'
            conc.get_outgroup(idx).group[0] = groups[0]

    # process the sequences themselves

    curr = 0

    for align_idx, (align, spc) in enumerate(zip(aligns, spacer)):

        # add sequence of this align + spacer

        ls = align.ls

        for main_idx, name in enumerate(names_in):
            sample_idx = samples_in[name][align_idx]
            if sample_idx != None:
                conc.get_sample(main_idx).sequence[curr:curr+ls] = align.get_sample(sample_idx).sequence
                conc.get_sample(main_idx).sequence[curr+ls:curr+ls+spc] = ch * spc
            else:
                conc.get_sample(main_idx).sequence[curr:curr+ls+spc] = ch * (ls+spc)

        for main_idx, name in enumerate(names_out):
            sample_idx = samples_out[name][align_idx]
            if sample_idx != None:
                conc.get_outgroup(main_idx).sequence[curr:curr+ls] = align.get_outgroup(sample_idx).sequence
                conc.get_outgroup(main_idx).sequence[curr+ls:curr+ls+spc] = ch * spc
            else:
                conc.get_outgroup(main_idx).sequence[curr:curr+ls+spc] = ch * (ls+spc)

        curr += ls + spc

    # return if needed

    if dest == None: return conc
