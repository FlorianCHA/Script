"""
The :func:`.phyml` function.
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

import os, re, subprocess
from .. import _interface, _tree
from ._utils import _App, _protect_run, paths

########################################################################

_example_file = """ 4 47
one       ACGGCGCTATAAATCGAGCGTTAGCGCGCGGAGAGAGACCTCTCTAG
two       AGGGCGCTATATATCGCGCGTTTGCGCGCGGAGAGAGAGCTCTCTAG
three     AGGGCGGTTTATATCGGGCGTTTTCGCCCGGTCCGAGAGCACTCTAG
four      ACGGCGCTTTATATCCGCCGTTAACGCCCGGTCCGAGGGCACTCTAG
"""

########################################################################

class _Phyml(_App):

    @_protect_run
    def _check_path(self, path):

        # test the "help" option to ensure that a phyml exist
        cmd = (path, '--help')
        try:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()
            if len(stderr): return stderr
        except OSError, e:
            return e.strerror

        # test that phyml works
        f = open('input.phy', 'w')
        f.write(_example_file)
        f.close()

        cmd = (path, '-i', 'input.phy', '-d', 'nt', '-b', '10',
               '-m', 'HKY85', '-c', '4', '-a', 'e')

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if len(stderr): return stderr
        if p.returncode != 0: return 'PhyML version is probably wrong (version 3.1 or higher is required)'
        if not os.path.isfile('input.phy_phyml_tree')       and not os.path.isfile('input.phy_phyml_tree.txt'):       return 'output file not generated'
        if not os.path.isfile('input.phy_phyml_stats')      and not os.path.isfile('input.phy_phyml_stats.txt'):      return 'output file not generated'
        if not os.path.isfile('input.phy_phyml_boot_stats') and not os.path.isfile('input.phy_phyml_boot_stats.txt'): return 'output file not generated'
        if not os.path.isfile('input.phy_phyml_boot_trees') and not os.path.isfile('input.phy_phyml_boot_trees.txt'): return 'output file not generated'

        # check the version
        cmd = (path, '--version')
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if len(stderr): return stderr
        mo = re.search('This is PhyML version 3\.2\.(\d+)', stdout)
        if mo == None:
            mo = re.search('This is PhyML version (\d+)', stdout)
            if mo == None:
                return 'cannot check version (invalid version string)'
            version = int(mo.group(1))
            if version < 20120000: return 'PhyML version 3.1 or higher is required'

        return None

    ####################################################################

_app = _Phyml(key='phyml', default='phyml')
paths._add(_app)

########################################################################

_nt_models = ['JC69', 'K80', 'F81', 'HKY85', 'F84', 'TN93', 'GTR']
_aa_models = ['LG', 'WAG', 'JTT', 'MtREV', 'Dayhoff', 'DCMut', 'RtREV',
              'CpREV', 'VT', 'Blosum62', 'MtMam', 'MtArt', 'HIVw',
              'HIVb']

########################################################################

@_protect_run
def phyml(align, model, include_outgroup=True, export_labels=False,
           rates=1, boot=0, start_tree='nj', fixed_topology=False,
           fixed_brlens=False, search='NNI', num_rand_start=5,
           freq=None, TiTv=4.0, pinv=0.0, alpha=None, use_median=False,
           free_rates=False, seed=None, verbose=False):

    """ phyml(align, model, include_outgroup=True, export_labels=True, \
              rates=1, boot=0, start_tree='nj', fixed_topology=False, \
              fixed_brlens=False, search='NNI', num_rand_start=5, \
              freq=None, TiTv=4.0, pinv=0.0, alpha=None, \
              use_median=False,  free_rates=False, seed=None, \
              verbose=False)

    .. include:: <isoamsa.txt>

    Reconstruct maximum-likelihood phylogeny using `PhyML <http://www.atgc-montpellier.fr/phyml/>`_.

    PhyML is a program performing maximum-likelihood phylogeny
    estimation using nucleotide or amino acid sequence alignments.
    **Ref.:** Guindon S., J.-F. Dufayard, V. Lefort, M. Animisova, W.
    Hordijk, and O. Gascuel. 2010. New algorithms and methods to
    estimate maximum-likelihood phylogenies: assessing the performance
    of PhyML 3.0. *Syst. Biol.* **59**\ : 307-321.

    :param align: input sequence alignment as an :class:`.Align`
        instance.
    :param model: substitution model to use (see list below).
    :param include_ougroup: boolean indicating whether sequences from
        the outgroup should be used for the phylogeny estimation.
    :param export_labels: boolean indicating whether the group labels
        should be included in the names of sequences (such as they will
        appear in the resulting tree).
    :param rates: number of discrete categories of evolutionary rate. If
        different of 1, fits a gamma distribution of rates.
    :param boot: number of bootstrap repetitions. Values of -1, -2 and
        -4 activate one the test-based branch support evaluation methods
        that provide faster alternatives to bootstrap repetitions (-1:
        aLRT statistics, -2: Chi2-based parametric tests, -4: Shimodaira
        and Hasegawa-like statistics). A value of 0 provides no branch
        support at all.
    :param start_tree: starting topology used by the program. Possible
        values are the string ``nj`` (neighbour-joining tree), the
        string ``pars`` (maximum-parsimony tree), the string ``rand``
        (random tree, only if *search* is ``SPR``), and a
        :class:`.Tree` instance containing a user-provided topology. In
        the latter case, the names of leaves of the tree must match the
        names of the input alignment (without group labels), implying
        that names cannot be repeated.
    :param fixed_topology: boolean indicating whether the topology
        provided in the :class:`.Tree` instance passed as *start_tree*
        argument should be NOT be improved.
    :param fixed_brlens: boolean indicating whether the branch lengths
        provided in the :class:`.Tree` instance passed as *start_tree*
        argument should be NOT be improved. All branch lengths must be
        specified. Automatically sets *fixed_topology* to ``True``.
    :param search: search method.  Available modes are: ``NNI``
        (nearest- neighbour interchange, fastest), ``SPR`` (subtree
        pruning and regrafting), and ``BEST`` (best of both methods).
    :param num_rand_start: number of random start trees (only considered
        if *search* is ``SPR`` and *start_tree* is ``rand``).
    :param freq: nucleotide or amino acid frequencies. Possible values
        are the string ``o`` (observed, frequencies measured from the
        data), the string ``m`` (estimated by maximum likelihood for
        nucleotides, or retrieved from the substitution model for amino
        acids), or a four-item tuples provided the relative frequencies
        of A, C, G and T respectively (only for nucleotides). By
        default, use ``o``  for nucleotides and ``m`` for amino acids.
    :param TiTv: transition/transversion ratio. If ``None``, estimated
        by maximum likelihood. Otherwise, must be a stricly positive
        value. Ignored if data are not nucleotides or if the model does
        not support 
    :param pinv: proportion of invariable sites. If ``None`` estimated
        by maximum likelihood. Otherwise, must be in the range [0, 1].
    :param alpha: gamma shape parameter. If ``None``, estimated by
        maximum likelihood. Otherwise, must be a strictly positive
        value. Ignored if *rates* is 1 or if *free_rates* is ``True``.
    :param use_median: boolean indicating whether the median (instead of
        the mean) should be use to report values for rate from the
        discretized gamma distribution. Ignored if *rates* is 1 or if
        *free_rates* is ``True``.
    :param free_rates: boolean indicating whether a mixture model should
        be used for substitution rate categories instead of the
        discretized gamma. In this case all reates and their frequencies
        will be estimated. Requires that *rates* is larger than 1.
    :param seed: pseudo-random number generator seed. Must be a stricly
        positive integer, preferably large.
    :param verbose: boolean indicating whether standard output of PhyML
        should be displayed.

    :return: A ``(tree, stats)`` tuple where *tree* is a :class:`.Tree`
        instance and *stats* is a :class:`dict` containing the following
        statistics or estimated parameters:

            * ``lk`` -- log-likelihood of the returned tree and model.
            * ``pars`` -- parsimony score of the returned tree.
            * ``size`` --  length of the returned tree.
            * ``rates`` -- relative substitution rates, as a list
              providing values in the following order: A\ |harr|\ C,
              A\ |harr|\ G, A\ |harr|\ T, C\ |harr|\ G, C\ |harr|\ T,
              and G\ |harr|\ T (only available if model is ``GTR`` or
              custom, only for nucleotide sequences).
            * ``alpha`` -- gamma shape parameter (only if the number of
              rate categories is larger than 1 and if ``free_rates`` was
              ``False``).
            * ``cats`` -- list of ``(rate, proportion)`` tuples for each
              discrete rate category (only if ``free_rates`` was
              ``True``, implying that the number of rates was larger
              than 1).
            * ``freq`` -- list of the relative base frequencies, in the
              following order: A, C, G, and T (only for nucleotide
              sequences).
            * ``ti/tv`` -- transition/transversion ratio (available for
              the models ``K80``, ``HKY85``, ``F84``, and ``TN93``). For
              the ``TN93`` model, the resulting value is a pair of
              transition/transversion ratios, one for purines and one
              for pyrimidines (in that order).
            * ``pinv`` -- proportion of invariable sites (only if the
              corresponding option was not set to 0).

    The choice of the model defines the type of data that are expected.
    The available models are:

    * Nucleotides:

        ========= ============================= ===== ================
        Code      Full name                     Rates Base frequencies
        ========= ============================= ===== ================
        ``JC69``  Jukes and Cantor 1969         one   equal
        ``K80``   Kimura 1980                   two   equal
        ``F81``   Felsenstein 1981              one   unequal
        ``HKY85`` Hasegawa, Kishino & Yano 1985 two   unequal
        ``F84``   Felsenstein 1984              two   unequal
        ``TN93``  Tamura and Nei 1993           three unequal
        ``GTR``   general time reversible       six   unequal
        ========= ============================= ===== ================

        In addition, custom nucleotide substitution models can be
        specified. In that case, *model* must be a six-character
        strings of numeric characters specifying which of the six
        (reversable) substitution rates are allowed to vary. The
        one-rate model is specified by the string ``000000``, the
        two-rate model (separate transition and transversion rates)
        is specified by ``010010``, and the GTR model is specified
        by ``012345``. The substitution rates are specified in the
        following order: A\ |harr|\ C, A\ |harr|\ G, A\ |harr|\ T,
        C\ |harr|\ G, C\ |harr|\ T, and G\ |harr|\ T.

    * Amino acids:

        ============ ======================================================================
        Code         Authors
        ============ ======================================================================
        ``LG``       Le & Gascuel (*Mol. Biol. Evol.* 2008)
        ``WAG``      Whelan & Goldman (*Mol. Biol. Evol.* 2001)
        ``JTT``      Jones, Taylor & Thornton (*CABIOS* 1992)
        ``MtREV``    Adachi & Hasegawa (*in* *Computer Science Monographs* 1996)
        ``Dayhoff``  Dayhoff *et al.* (*in* *Atlas of Protein Sequence and Structure* 1978)
        ``DCMut``    Kosiol & Goldman (*Mol. Biol. Evol.* 2004)
        ``RtREV``    Dimmic *et al.* (*J. Mol. Evol.* 2002)
        ``CpREV``    Adachi *et al.* (*J. Mol. Evol.* 2000)
        ``VT``       Muller & Vingron (*J. Comput. Biol.* 2000)
        ``Blosum62`` Henikoff & Henikoff (*PNAS* 1992)
        ``MtMam``    Cao *et al.* (*J. Mol. Evol.* 1998)
        ``MtArt``    Abascal, Posada & Zardoya (*Mol. Biol. Evol.* 2007)
        ``HIVw``     Nickle *et al.* (*PLoS One* 2007)
        ``HIVb``     *ibid.*
        ============ ======================================================================

    .. versionchanged:: 3.0.0

        No more default value for *model* option. Added custom model for
        nucleotides. Changed SH pseudo-bootstrap option flag from -3 to
        -4. *quiet* function replaced by *verbose*. Several additional
        options are added. The syntax for input a user tree is modified.
        The second item in the returned tuple is a dictionary of
        statistics.
    """

    # check that program is available
    path = _app.get_path()
    if path is None:
        raise RuntimeError, 'PhyML program not available -- please configure path'
    command_line = [path]

    # check that alignment has enough sequences
    if not isinstance(align, _interface.Align):
        raise TypeError, '`align` should be an `Align` instance'
    n = align.ns
    if include_outgroup: n += align.no
    if n < 3:
        raise ValueError, 'alignment does not contain enough sequences to run PhyML'
    if align.ls < 1:
        raise ValueError, 'alignment does not contain enough data to run PhyML'

    # encode and save alignment
    f = open('i', 'w')
    f.write('{0} {1}\n'.format(n, align.ls))
    mapping1 = {} # without labels
    mapping2 = {} # with labels if export_labels==True (same as mapping1 otherwise)
    for i, item in enumerate(align):
        name = 'in_{0}'.format(i+1)
        if export_labels:
            mapping1[name] = item.name
            mapping2[name] = '{0} @{1}'.format(item.name, ','.join(map(str, item.group)))
        else:
            mapping1[name] = item.name
        f.write('{0}  {1}\n'.format(name, item.sequence.string()))
    if include_outgroup:
        for i, item in enumerate(align.iter_outgroup()):
            name = 'out_{0}'.format(i+1)
            if export_labels:
                mapping1[name] = item.name
                mapping2[name] = '{0} @#{1}'.format(item.name, item.group[0])
            else:
                mapping1[name] = item.name
            f.write('{0}  {1}\n'.format(name, item.sequence.string()))
    f.close()
    if not export_labels: mapping2 = mapping1

    command_line.extend(['-i', 'i'])

    # get model and determines data type
    custom_flag = False
    if model in _aa_models: dt = 'aa'
    elif model in _nt_models: dt = 'nt'
    else:
        if not isinstance(model, basestring): raise TypeError, '`model` must be a string'
        if len(model) != 6: raise ValueError, 'invalid model name or string'
        try: x = map(int, model)
        except ValueError: raise ValueError, 'invalid model name or string'
        if range(len(set(x))) != sorted(set(x)): raise ValueError, 'invalid custom model specification'
        custom_flag = True
        dt = 'nt'
    command_line.extend(['-d', dt, '-m', model])

    # check boot option
    if not isinstance(boot, int): raise TypeError, '`boot` argument must be an int'
    if  boot < -4 and boot != -3: raise ValueError, 'invalid value for `boot` argument'
    command_line.extend(['-b', str(boot)])

    # check ti/tv
    if TiTv == None: TiTv = 'e'
    elif TiTv <= 0.0: raise ValueError, 'invalid value for `TiTv` argument'
    command_line.extend(['-t', str(TiTv)])

    # manage freq option
    if freq == None:
        if dt == 'nt': freq = 'o'
        else: freq = 'm'
    elif freq == 'o': pass
    elif freq == 'm': pass
    else:
        if len(freq) != 4: raise ValueError, 'invalid value for `freq` argument'
        if abs(sum(freq) - 1.0) > 0.000000001: raise ValueError, 'invalid value for `freq` argument'
        freq = ','.join(map(str, freq))
    command_line.extend(['-f', str(freq)])

    # check pinv option
    if pinv == None: pinv = 'e'
    elif pinv < 0.0 or pinv > 1.0: raise ValueError, 'invalid value for `pinv` argument'
    command_line.extend(['-v', str(pinv)])

    # set gamma distribution parameters
    if not isinstance(rates, int): raise TypeError, '`rates` argument must be an int'
    if rates < 1: raise ValueError, 'invalid value for `rates` argument'
    if alpha == None: alpha = 'e'
    elif not isinstance(alpha, float): raise TypeError, '`alpha` argument must be a float'
    if alpha <= 0.0: raise ValueError, 'invalid value for `alpha` argument'
    command_line.extend(['-c', str(rates), '-a', str(alpha)])
    if free_rates:
        if rates == 1: raise ValueError, 'cannot use `free_rates` if `rates` < 2'
        command_line.append('--free_rates')

    # input tree
    user_tree_flag = False
    if start_tree == 'nj': pass
    elif start_tree == 'pars': command_line.append('--pars')
    elif start_tree == 'rand': command_line.append('--rand_start')
    else:
        # if input tree is Tree, make a copy and rename all leaves
        if not isinstance(start_tree, _tree.Tree): raise ValueError, 'invalid value for `start_tree` argument'

        # try to reverse mapping (will not work if sequences are ambiguous in input)
        rmap = {}
        for k, v in mapping1.iteritems():
            if v in rmap: raise ValueError, 'cannot use a user tree because of replicates'
            rmap[v] = k

        # rename all leaves
        for leaf in start_tree.iter_leaves():
            lbl = leaf.label
            if lbl not in rmap: raise ValueError, 'invalid user tree: leaves do not match alignment'
            leaf.label = rmap[lbl]
            del rmap[lbl]
        if len(rmap) > 0: raise ValueError, 'invalid user tree: not enough leaves'

        # convert tree to string (check that branch lengths are there if required
        if fixed_brlens:
            for node in start_tree.depth_iter():
                for child in node.children():
                    if node.branch_to(child) is None: raise ValueError, '`fixed_brlens` requires branch length values'
            f = open('t', 'w')
            f.write(start_tree.newick(labels=False, brlens=True) + '\n')
            f.close()
        else:
            f = open('t', 'w')
            f.write(start_tree.newick(labels=False, brlens=False) + '\n')
            f.close()

        # set argument
        command_line.extend(['-u', 't'])
        user_tree_flag = True

    # optimization option
    if fixed_topology or fixed_brlens and not user_tree_flag: raise ValueError, '`fixed_topology` option requires a user tree'
    if fixed_brlens: command_line.extend(['-o', 'r'])
    elif fixed_topology: command_line.extend(['-o', 'lr'])
    else: command_line.extend(['-o', 'tlr'])

    # search option
    if search not in ('NNI', 'SPR', 'BOTH'): raise ValueError, 'invalid value for `search` argument'
    command_line.extend(['--search', search])
    if search == 'SPR' and start_tree == 'rand':
        if not isinstance(num_rand_start, int): raise TypeError, '`num_rand_start` argument must be an int'
        if num_rand_start < 1: raise ValueError, 'invalid value for `num_rand_start` argument'
        command_line.extend(['-n_rand_starts', str(num_rand_start)])

    # random seed
    if seed is not None:
        if not isinstance(seed, int): raise TypeError, '`seed` argument must be an int'
        if seed < 1: raise ValueError, 'invalid value for `seed` argument'
        command_line.extend(['--r_seed_num', str(seed)])

    # set other options
    if use_median: command_line.append('--use_median')
    if not verbose: command_line.append('--quiet')
    command_line.append('--no_memory_check')

    # run the program
    p = subprocess.Popen(command_line, stdin=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        stdout=(None if verbose else subprocess.PIPE))

    stdout, stderr = p.communicate('Y\n')

    # check error
    stderr = stderr.strip()
    if len(stderr):
        raise RuntimeError, 'error while running PhyML: {0}'.format(stderr)

    if os.path.isfile('i_phyml_stats.txt'): statsf = 'i_phyml_stats.txt'
    elif os.path.isfile('i_phyml_stats'): statsf = 'i_phyml_stats'
    else: raise RuntimeError, 'unknown error while running PhyML (try running in verbose mode)'

    if os.path.isfile('i_phyml_tree.txt'): treef = 'i_phyml_tree.txt'
    elif os.path.isfile('i_phyml_tree'): treef = 'i_phyml_tree'
    else: raise RuntimeError, 'unknown error while running PhyML (try running in verbose mode)'

    # get resulting tree
    try:
        tree = _tree.Tree(fname=treef)
    except ValueError:
        raise RuntimeError, 'cannot import PhyML output (try running in verbose mode)'

    # decode tree leaves
    for leaf in tree.iter_leaves():
        leaf.label = mapping2[leaf.label]

    # import statistics
    stats = {}
    f = open(statsf)
    string = f.read()
    f.close()

    # likelihood
    mo = re.search('\. Log-likelihood\:[ \t]+([\-\dEe\.]+)', string)
    if mo:
        try: stats['lk'] = float(mo.group(1))
        except ValueError: raise RuntimeError, 'cannot parse PhyML output (invalid log-likelihood expression)'
    else:
        raise RuntimeError, 'no log-likelihood found in PhyML output'

    # parsimony score
    mo = re.search('\. Parsimony\:[ \t]+([\d]+)', string)
    if mo:
        try: stats['pars'] = int(mo.group(1))
        except ValueError: raise RuntimeError, 'cannot parse PhyML output (invalid parsimony score)'
    else:
        raise RuntimeError, 'no parsimony score found in PhyML output'

    # tree size
    mo = re.search('\. Tree size\:[ \t]+([\.\dEe\-]+)', string)
    if mo:
        try: stats['size'] = float(mo.group(1))
        except ValueError: raise RuntimeError, 'cannot parse PhyML output (invalid tree size)'
    else:
        raise RuntimeError, 'no tree size found in PhyML output'

    # nucleotide frequencies
    if dt == 'nt':
        stats['freqs'] = []
        for base in 'A', 'C', 'G', 'T':
            mo = re.search('f\({0}\)= ([\.\dEe\-]+)'.format(base), string)
            if mo:
                try: stats['freqs'].append(float(mo.group(1)))
                except ValueError: raise RuntimeError, 'cannot parse PhyML output (invalid nucleotide frequency)'
            else:
                raise RuntimeError, 'missing frequencies in PhyML output'

    # nucleotide substitution rate
    if model == 'GTR' or custom_flag:
        stats['rates'] = []
        for a, b in [('A', 'C'), ('A', 'G'), ('A', 'T'), ('C', 'G'), ('C', 'T'), ('G', 'T')]:
            mo = re.search('{0} <-> {1} [ \t]+([\.\dEe\-]+)'.format(a, b), string)
            if mo:
                try: stats['rates'].append(float(mo.group(1)))
                except ValueError: raise RuntimeError, 'cannot parse PhyML output (invalid nucleotide substitution rate)'
            else:
                raise RuntimeError, 'missing substitution rates in PhyML output'

    # alpha parameter or discrete rates
    if rates > 1:
        if not free_rates:
            mo = re.search('- Gamma shape parameter:[ \t]+([\d\.Ee\-]+)', string)
            if mo:
                try: stats['alpha'] = float(mo.group(1))
                except ValueError: raise RuntimeError, 'cannot parse PhyML output (invalid gamma shape parameter)'
            else:
                raise RuntimeError, 'missing gamma shape parameter in PhyML output'
        else:
            stats['cats'] = []
            for i in xrange(rates):
                mo = re.search('\- Relative rate in class {0}:[ \t]+([\.\dEe\-]+) \[prop\=([\.\dEe\-]+)\]'.format(i+1), string)
                if mo:
                    r, p = mo.groups()
                    try: stats['cats'].append((float(r), float(p)))
                    except ValueError: raise RuntimeError, 'cannot parse PhyML output (invalid discrete rates parameter)'
                else:
                    raise RuntimeError, 'missing discrete rate'

    # ti/tv ratio
    if model in ['K80', 'HKY85', 'F84']:
        mo = re.search('\. Transition/transversion ratio:[ \t]+([\.\dEe\-]+)', string)
        if mo:
            try: stats['ti/tv'] = float(mo.group(1))
            except ValueError: raise RuntimeError, 'cannot parse PhyML output (invalid Ti/Tv ratio)'
        else:
            raise RuntimeError, 'missing Ti/Tv ratio in PhyML output'
    if model == 'TN93':
        stats['ti/tv'] = [None, None]
        mo = re.search('\. Transition/transversion ratio for purines:[ \t]+([\.\dEe\-]+)', string)
        if mo:
            try: stats['ti/tv'][0] = float(mo.group(1))
            except ValueError: raise RuntimeError, 'cannot parse PhyML output (invalid Ti/Tv ratio)'
        else:
            raise RuntimeError, 'missing Ti/Tv ratio in PhyML output'

        mo = re.search('\. Transition/transversion ratio for pyrimidines:[ \t]+([\.\dEe\-]+)', string)
        if mo:
            try: stats['ti/tv'][1] = float(mo.group(1))
            except ValueError: raise RuntimeError, 'cannot parse PhyML output (invalid Ti/Tv ratio)'
        else:
            raise RuntimeError, 'missing Ti/Tv ratio in PhyML output'

    # proportion of fixed sites
    if pinv != 0.0: # including None
        mo = re.search('\. Proportion of invariant:[ \t]+([\.\dEe\-]+)', string)
        if mo:
            try: stats['pinv'] = float(mo.group(1))
            except ValueError: raise RuntimeError, 'cannot parse PhyML output (invalid fixed sites proportion)'
        else:
            raise RuntimeError, 'missing fixed sites proportion in PhyML output'

    # return
    return tree, stats
