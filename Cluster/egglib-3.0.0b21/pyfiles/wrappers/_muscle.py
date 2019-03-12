"""
The :func:`.muscle` function.
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

import subprocess, re, os
from .. import _interface
from ..io import _fasta
import _utils

########################################################################

class _Muscle(_utils._App):

    @_utils._protect_run
    def _check_path(self, path):

        # test the "help" option to ensure that a phyml exist
        cmd = (path, '-version')
        try:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()
            if len(stderr): return stderr
        except OSError, e:
            return e.strerror

        mo = re.match('MUSCLE v(\d+)\.(\d+)\.(\d+)', stdout)
        if mo is None: return 'cannot read version number from Muscle'
        M = int(mo.group(1))
        m = int(mo.group(2))
        if M != 3 or m < 8: return 'Muscle version 3.8 or higher is required'
        return None

    ####################################################################

_app = _Muscle(key='muscle', default='muscle')
_utils.paths._add(_app)

########################################################################

_args = {
    # these ones cannot passed as is to muscle
    'diags':           (bool,  None),
    'diags1':          (bool,  None),
    'diags2':          (bool,  None),
    'anchors':         (bool,  None),
    'dimer':           (bool,  None),
    'brenner':         (bool,  None),
    'cluster':         (bool,  None),
    'teamgaps4':       (bool,  lambda x: x is True),
    'nt_profile':      (str,   ['spn']),
    'aa_profile':      (str,   ['le', 'sp', 'sv']),

    # these ones should be able to be passed as is
    'anchorspacing':   (int,   lambda x: x>0),
    'center':          (float, lambda x: x<0),
    'cluster1':        (str,   ['upgma', 'upgmb', 'neighborjoining']),
    'diagbreak':       (int,   lambda x: x>0),
    'diaglength':      (int,   lambda x: x>0),
    'diagmargin':      (int,   lambda x: x>0),
    'diagmargin':      (int,   lambda x: x>0),
    'distance1':       (str,   ['kmer6_6', 'kmer20_3', 'kmer20_4', 'kbit20_3', 'kmer4_6']),
    'distance2':       (str,   ['pctidkimura', 'pctidlog']),
    'gapopen':         (float, lambda x: x<0),
    'hydro':           (int,   lambda x: x>0),
    'hydrofactor':     (float, lambda x: x>0),
    'maxiters':        (int,   lambda x: x>0),
    'maxtrees':        (int,   lambda x: x>0),
    'minbestcolscore': (float, lambda x: x>=0),
    'minsmoothscore':  (float, lambda x: x>=0),
    'objscore':        (str,   ['sp', 'ps', 'dp', 'xp', 'spf', 'spm']),
    'refinewindow':    (int,   lambda x: x>0),
    'root1':           (str,   ['pseudo', 'midlongestspan', 'minavgleafdist']),
    'seqtype':         (str,   ['protein', 'dna', 'rna', 'auto']),
    'smoothscoreceil': (float, lambda x: x>0),
    'SUEFF':           (float, lambda x: x>=0.0 and x<=1.0),
    'weight1':         (str,   ['none', 'henikoff', 'henikoffpb', 'gsc', 'clustalw', 'threeway']),
    'weight2':         (str,   ['none', 'henikoff', 'henikoffpb', 'gsc', 'clustalw', 'threeway'])
}



########################################################################

@_utils._protect_run
def muscle(source, ref=None, include_outgroup=True, verbose=False, **kwargs):

    """ muscle(source, ref=None, include_outgroup=True, verbose=False, \
                **kwargs)

    Performs multiple alignment using `Muscle <http://www.drive5.com/muscle/>`_.

    Muscle's default options tend to produce high-quality alignments but
    may be long to run on large data set. Muscle's author recommends
    using the option ``maxiters=2`` for large data sets, and, for fast
    alignment (in particular of closely related sequences):
    ``maxiters=1 diags=True aa_profile='sv' distance1='kbit20_3'`` (for
    amino acid sequences) and ``maxiters=1 diags=True`` (for nucleotide
    sequences).

    :param source: a :class:`.Container` or :class:`.Align` containing
        sequences to align. If an :class:`!Align` is provided, sequences
        are assumed to be already aligned and alignment will be refined
        (using the ``-refine`` option of Muscle), unless an alignment is
        also provided as *ref*. In the latter case,  the two alignments
        are preversed (their columns are left unchanged), and they
        aligned with respect to each other.

    :param ref: an :class:`.Align` instance providing an alignement that
        should be aligned with respect to the alignment provided as
        *source*. If *ref* is provided, it is required both *source* and
        *ref* are :class:`!Align` instances.

    :param include_outgroup: boolean indicating whether sequences from
        the outgroup should be included in the alignment.

    :param verbose: display Muscle's console output.

    :param kwargs: other keyword arguments are passed to Muscle. The
        available options are listed below:

        .. include:: muscle_arguments.txt

        For a description of options, see the `Muscle manual <http://www.drive5.com/muscle/manual/options.html>`_.
        Most of Muscle's options are available. Note that function takes
        no flag option, and Muscle's flag options are passed as boolean
        keyword arguments (except options relative to the amino acid or
        nucleotide profile score options, that are passed as string as
        ``aa_profile`` and ``nt_profile``, respectively. The order of
        options is preserved.

    :return: An :class:`.Align` containing aligned sequences.

    .. versionchanged:: 3.0.0

        Added support for most options.
    """

    # check that program is available
    path = _app.get_path()
    if path is None:
        raise RuntimeError, 'Muscle program not available -- please configure path'
    command_line = [path, '-out', 'o']
    if not verbose: command_line.append('-quiet')

    # mapping for sample reference
    mapping = {}

    # write source
    if not isinstance(source, (_interface.Container, _interface.Align)): raise TypeError, 'invalid type for `source` argument'
    _utils._write(source, include_outgroup, 'source', mapping)

    # then save ref
    if ref != None:
        if not isinstance(ref, _interface.Align): raise TypeError, 'invalid type for `ref` argument'
        if not isinstance(source, _interface.Align): raise TypeError, 'if `ref` is provided, `source` must be an Align'
        _utils._write(ref, include_outgroup, 'ref', mapping)
        command_line.extend(['--profile', '-in1', 'ref', '-in2', 'source'])

    elif isinstance(source, _interface.Align):
        command_line.extend(['-in', 'source', '-refine'])

    else: # the default case
        command_line.extend(['-in', 'source'])

    # process kwargs
    for opt, value in kwargs.iteritems():
        if opt not in _args: raise ValueError, 'invalid option: {0}'.format(opt)
        type_, test = _args[opt]
        try: value = type_(value)
        except (TypeError, ValueError): raise ValueError, 'invalid value for argument {0}'.format(opt)
        if type_ == str:
            if value not in test: raise ValueError, 'invalid value for argument {0}'.format(opt)
        elif test is not None and not test(value): raise ValueError, 'invalid value for argument {0}'.format(opt)

        if opt == 'diags':
            if value == True: command_line.append('-diags')
        elif opt == 'diags1':
            if value == True: command_line.append('-diags1')
        elif opt == 'diags2':
            if value == True: command_line.append('-diags2')
        elif opt == 'anchors':
            if value == True: command_line.append('-anchors')
            else: command_line.append('-noanchors')
        elif opt == 'dimer':
            if value == True: command_line.append('-dimer')
        elif opt == 'brenner':
            if value == True: command_line.append('-brenner')
        elif opt == 'cluster':
            if value == True: command_line.append('-cluster')
        elif opt == 'teamgaps4':
            if value == True: command_line.append('-teamgaps4')
        elif opt == 'nt_profile':
            command_line.append('-' + value)
        elif opt == 'aa_profile':
            command_line.append('-' + value)
        else:
            command_line.extend(['-' + opt, str(value)])

    # run the program
    p = subprocess.Popen(command_line, stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=(None if verbose else subprocess.PIPE))
    stdout, stderr = p.communicate('Y\n')

    # check error
    if not os.path.isfile('o'):
        if stderr is not None:
            lines = stderr.splitlines()
            if len(lines) and 'ERROR' in lines[-1]:
                raise RuntimeError, 'error while running Muscle: {0}'.format(lines[-1].strip())
        raise RuntimeError, 'unknown error while running Muscle (try running in verbose mode)'

    # get alignment
    aln = _fasta.from_fasta('o', groups=False, cls=_interface.Align)
    if ref != None and ref.ng > source.ng: aln.ng = ref.ng
    else: aln.ng = source.ng

    # move outgroup samples to outgroup
    for code, ref in mapping.iteritems():
        if ref.outgroup:
            aln.to_outgroup(aln.find(code, index=True), ref.group[0])

    # set proper names and groups
    for sam in aln.iter_samples(True, True):
        ref = mapping[sam.name]
        sam.name = ref.name
        if not ref.outgroup:
            for i,v in enumerate(ref.group): sam.group[i] = v

    # return
    return aln
