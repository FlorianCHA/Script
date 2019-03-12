"""
Class to hold (and verify) parameters for coalescent simulations.
"""

__license__ = """
    Copyright 2015-2018 Stephane De Mita, Mathieu Siol

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

from .. import _eggwrapper, stats

########################################################################

class ParamDict(object):

    """
    Emulator of :class:`dict` for managing parameters. Do most of what a
    dictionary does except for adding and removing keys. The order of
    parameters is fixed and consistent.

    In addition to the methods documented below, :class:`.ParamDict`
    instances support the following operations (where ``params`` in one
    instance):

    +------------------------+--------------------------------------------------------+
    |Expression              | Action                                                 |
    +========================+========================================================+
    |``len(params)``         | Number of parameters                                   |
    +------------------------+--------------------------------------------------------+
    |``params[key]``         | Get the value of parameter ``key``                     |
    +------------------------+--------------------------------------------------------+
    |``params[key] = value`` | Assign ``value`` to parameter ``key`` (see note below) |
    +------------------------+--------------------------------------------------------+
    |``for key in params``   | Same as ``for key in params.iterkeys()``               |
    +------------------------+--------------------------------------------------------+
    |``reversed(params)``    | Return a reversed iterator                             |
    +------------------------+--------------------------------------------------------+
    |``key in params``       | Check if ``key`` in a parameter name                   |
    +------------------------+--------------------------------------------------------+
    |``str(params)``         | Representation of the instance as a :class:`dict`      |
    +------------------------+--------------------------------------------------------+

    Note that the ``params[key] = value`` expression is straightforward
    only for parameters that have a single value. For parameters
    represented by a :class:`.ParamList` instance, this expression is
    only supported if the right-hand operand is a sequence of matching
    length (in order to set all values at once). For parameters
    represented by a  :class:`.ParamMatrix` instance, this expression is
    only supported if the right-hand operand is another
    :class:`.ParamMatrix` instance of matching dimension. For
    ``events``, this expression is not supported at all. In all cases
    where ``params[key]`` returns :class:`.ParamList`,
    :class:`.ParamMatrix` or a :class:`.EventList` instance, the
    returned value can be modified using its own methods.
    """

    ####################################################################

    _keys = ['num_pop', 'num_sites', 'recomb', 'theta', 'num_mut',
             'mut_model', 'TPM_proba', 'TPM_param', 'num_alleles',
             'rand_start', 'num_chrom', 'num_indiv', 'N', 'G', 's',
             'site_pos', 'site_weight', 'migr_matrix', 'trans_matrix',
             'events', 'seed', 'outgroup', 'max_iter']

    # sub-lists with only parameters that can be changed between
    # simulations

    _keys_scalar = ['num_sites', 'recomb', 'theta', 'num_mut',
             'mut_model', 'TPM_proba', 'TPM_param', 'num_alleles',
             'rand_start', 'outgroup']

    _keys_list = ['num_chrom', 'num_indiv', 'N', 'G', 's', 'site_pos',
                  'site_weight']

    _keys_matrix = ['migr_matrix', 'trans_matrix']

    _random_seeds = set()

    ####################################################################

    def __init__(self, npop):

        self._random = _eggwrapper.Random()
        self._seed = self._random.get_seed() # first because __del__ requires it

        if npop < 1: raise ValueError, 'invalid value for `npop`'

        self._params = _eggwrapper.Params(npop, 0.0)
        self._migr = self._params.M()
        self._params.set_L(0)
        self._npop = npop
        while self._seed in self._random_seeds:
            self._seed = self._random.irand(10**6) + 10**4
            self._random.set_seed(self._seed)
        self._random_seeds.add(self._seed)

        self._num_chrom = ParamList(getter=self._params.get_n1,
                                    setter=self._params.set_n1,
                                    num_values=npop,
                                    check_item=lambda x: x>=0)

        self._num_indiv = ParamList(getter=self._params.get_n2,
                                    setter=self._params.set_n2,
                                    num_values=npop,
                                    check_item=lambda x: x>=0)

        self._N = ParamList(getter=self._params.get_N,
                            setter=self._params.set_N,
                            num_values=npop,
                            check_item=lambda x: x>0.0)

        self._G = ParamList(getter=self._params.get_G,
                            setter=self._params.set_G,
                            num_values=npop,
                            check_item=lambda x: True)

        self._s = ParamList(getter=self._params.get_s,
                            setter=self._params.set_s,
                            num_values=npop,
                            check_item=lambda x: x>=0.0 and x<=1.0)

        self._site_pos = ParamList(getter=self._params.get_sitePos,
                                   setter=self._params.set_sitePos,
                                   num_values=0,
                                   check_item=lambda x: x>=0.0 and x<=1.0)

        self._site_weight = ParamList(getter=self._params.get_siteW,
                                      setter=self._params.set_siteW,
                                      num_values=0,
                                      check_item=lambda x: x>0.0)

        self._migr_matrix = ParamMatrix(getter=self._migr.get_pair,
                                        setter=self._migr.set_pair,
                                        num_values=npop,
                                        check_item=lambda x:x>=0.0,
                                        check_active=lambda: True,
                                        set_active=None)

        self._trans_matrix = ParamMatrix(getter=self._params.get_transW_pair,
                                         setter=self._params.set_transW_pair,
                                         num_values=2,
                                         check_item=lambda x:x>0.0,
                                         check_active=lambda: self._params.get_transW_matrix(),
                                         set_active=self._params.set_transW_matrix)

        self._events = EventList(npop=npop,
                                 add=self._params.addChange,
                                 clear=self._params.clearChanges)

    ####################################################################

    def __del__(self):
        self._random_seeds.discard(self._seed)

    ####################################################################

    def disable_trans_matrix(self):

        """
        Disable the matrix of weights of transition between alleles (it
        is disabled by default, and automatically activated if any value
        is set). This sets all weights to 1.0.
        """

        self._params.set_transW_matrix(False)

    ####################################################################

    def get_values(self, other):

        """
        Update the dictionary with values from another
        :class:`.ParamDict` instance (with the same numbers of
        populations and sites).
        """

        if not isinstance(other, ParamDict): raise TypeError, '`other` should be a `ParamDict` instance'
        self['num_sites'] = other['num_sites']
        self['num_alleles'] = other['num_alleles']
        for k in self._keys:
            if k not in ['num_pop', 'num_sites', 'num_alleles', 'events']: self[k] = other[k]
        self._events.replace(other._events)

    ####################################################################

    def summary(self):

        """
        Return a string displaying all current values of parameters at
        the level of the C++ simulator. Use for debugging only, as the
        format is not guaranteed to be stable.
        """

        return self._params.summary()

    ####################################################################

    def _set(self, key, values):
        # alternative setter

        if key in self._keys_scalar: self[key] = values
        elif key in self._keys_list: self[key][values[0]] = values[1]
        elif key in self._keys_matrix: self[key][values[0], values[1]] = values[2]
        else: raise KeyError, 'this parameter cannot be modified between simulations'

    ####################################################################

    def set_migr(self, value):

        """
        Set all pairwise migration rates to ``value/(num_pop-1)``.
        """

        if value < 0.0: raise ValueError, 'invalid value for `migr`'
        self._migr.set_all(value)

    ####################################################################

    def mk_structure(self, skip_indiv=False):

        """
        Export a :class:`.stats.Structure` instance containing the
        structure information corresponding to currently loaded
        simulation parameters.

        :param skip_indiv: this argument determines whether the individuals
            level should be skipped (if ``True``, alleles from each given
            individual are treated as belonging to separate haploid individuals.

            .. warning::
                The individual level cannot be processed if the ploidy is
                not constant (mixing sampled chromosomes and sampled
                individuals).

        :return: A new :class:`.stats.Structure` instance.
        """

        # record outgroup label
        otg = self._params.get_outgroup()

        # process populations
        struct = {0: {}}
        outgroup = {}
        cnt_sam_ing = 0
        cnt_sam_otg = 0
        cnt_idv_ing = 0
        cnt_idv_otg = 0
        for k in xrange(self._npop):

            # get list of sample indexes (as a list of tuples/per individual)
            if k == otg: sam_start = cnt_sam_otg
            else: sam_start = cnt_sam_ing

            idv = [(sam_start+i*2, sam_start+i*2+1) for i in xrange(self._params.get_n2(k))]
            sam_start += 2*len(idv)
            idv.extend([(sam_start+i,) for i in xrange(self._params.get_n1(k))])

            if k == otg: cnt_sam_otg += 2*self._params.get_n2(k) + self._params.get_n1(k)
            else: cnt_sam_ing += 2*self._params.get_n2(k) + self._params.get_n1(k)

            # process individuals if outgroup
            if k == otg:
                if not skip_indiv: # do not merge the two "if"'s
                    outgroup.update(dict(zip(xrange(cnt_idv_otg, cnt_idv_otg+len(idv)),idv)))
                    cnt_idv_otg += len(idv)

            # process population
            else:
                if skip_indiv:
                    idv = [(j,) for i in idv for j in i]
                idv = dict(zip(xrange(cnt_idv_ing, cnt_idv_ing+len(idv)), idv))
                struct[0][k] = idv
                cnt_idv_ing += len(idv)

        # process delayed samples
        if self._params.nDSChanges() > 0:
            for event in self._events:
                if event['cat'] == 'sample':
                    k = event['label']

                    # make a copy of the right counters
                    if k == otg:
                        if skip_indiv: continue
                        cnt_sam = cnt_sam_otg
                    else:
                        cnt_sam = cnt_sam_ing

                    idv = [(cnt_sam+i*2, cnt_sam+i*2+1) for i in xrange(event['num_indiv'])]
                    idv.extend([(cnt_sam+i,) for i in xrange(event['num_chrom'])])

                    if k == otg: cnt_sam_otg += 2 * event['num_indiv'] + event['num_chrom']
                    else: cnt_sam_ing += 2 * event['num_indiv'] + event['num_chrom']

                    # process individuals if outgroup
                    if k == otg:
                        outgroup.update(dict(zip(xrange(cnt_idv_otg, cnt_idv_otg+len(idv)), idv)))
                        cnt_idv_otg += len(idv)

                    # process population otherwise
                    else:
                        if k not in struct[0]: struct[0][k] = {}
                        if skip_indiv: idv = [(j,) for i in idv for j in i]
                        idv = dict(zip(xrange(cnt_idv_ing, cnt_idv_ing+len(idv)), idv))
                        struct[0][k].update(idv)
                        cnt_idv_ing += len(idv)
 
        # check that mapping is non-empty
        if len(struct[0]) == 0: raise ValueError, 'cannot generate structure (no population left in outgroup)'

        # generate structure
        return stats.make_structure(struct, dict([(k+cnt_idv_ing, v) for k,v in outgroup.iteritems()]))

    ####################################################################
    # below: dict-like methods
    ####################################################################

    def keys(self):

        """
        Return a copy of the list of parameter names.
        """

        return list(self._keys)

    ####################################################################

    def values(self):

        """
        Return a copy of the list of parameter values.
        """

        return [self[key] for key in self._keys]

    ####################################################################

    def items(self):

        """
        Return a copy of the list of ``(key, value)`` parameter name and
        value pairs.
        """

        return [(key, self[key]) for key in self._keys]

    ####################################################################

    def has_key(self, key):

        """
        Return a boolean indicating if the passed name is one of the
        parameters.
        """

        return key in self._keys

    ####################################################################

    def get(self, key, default=None):

        """
        Return the value for *key* if *key* is one of the parameters,
        else return the value passed as *default* (by default,
        ``None``). This method therefore never raises a
        :exc:`~exceptions.KeyError`.
        """

        if key in self._keys: return self[key]
        return default

    ####################################################################

    def __iter__(self):
        for key in self._keys: yield key

    ####################################################################

    def iterkeys(self):

        """
        Return an iterator over the parameter names.
        """

        for key in self._keys: yield key

    ####################################################################

    def itervalues(self):

        """
        Return an iterator over the parameter values.
        """

        for key in self._keys: yield self[key]

    ####################################################################

    def iteritems(self):

        """
        Return an iterator over the ``(key, value)`` parameter name and
        value pairs.
        """

        for key in self._keys: yield key, self[key]

    ####################################################################

    def copy(self):

        """
        Return a shallow copy of this instance. All modifications of the
        values of either copy will modify the other one (even
        modification of integer, float and string parameters).
        """

        ret = ParamDict(self._npop)
        ret._params = self._params
        return ret

    ####################################################################

    def update(self, other=None, **kwargs):

        """
        Update the dictionary with the ``(key, value)`` parameter name
        and value pairs from the object passed as *other*, overwriting
        existing keys and raising a :exc:`~exceptions.KeyError`
        exception if an unknown parameter name is met.

        Accepts :class:`dict` instances, or else any iterable of
        ``(key, value)`` pairs. If keyword arguments are specified, the
        instance is then updated with those, as in
        ``param_dict.update(theta=5.0, recomb=2.5)``.
        """

        if isinstance(other, ParamDict):
            self.get_values(other)

        elif other is not None:
            if 'num_sites' in other: self['num_sites'] = other['num_sites']
            if 'num_alleles' in other: self['num_alleles'] = other['num_alleles']
            try:
                iter_ = other.iteritems()
            except AttributeError: iter_ = iter(other)
            for k, v in iter_: self[k] = v

        if 'num_sites' in kwargs:
            self['num_sites'] = kwargs['num_sites']
            del kwargs['num_sites']
        if 'num_alleles' in kwargs:
            self['num_alleles'] = kwargs['num_alleles']
            del kwargs['num_alleles']

        for k, v in kwargs.iteritems(): self[k] = v

    ####################################################################

    def __len__(self):
        return len(self._keys)

    ####################################################################

    def __reversed__(self):
        for key in reversed(self._keys):
            yield key

    ####################################################################

    def __contains__(self, key):
        return key in self._keys

    ####################################################################

    def __str__(self):
        return str(dict(self.iteritems()))

    ####################################################################

    def __getitem__(self, key):
        if key == 'num_pop': return self._npop
        elif key == 'num_sites': return self._params.get_L()
        elif key == 'recomb': return self._params.get_R()
        elif key == 'theta': return self._params.get_theta()
        elif key == 'num_mut': return self._params.get_fixed()
        elif key == 'mut_model':
            model = self._params.get_mutmodel()
            if model == _eggwrapper.Params.KAM: return 'KAM'
            elif model == _eggwrapper.Params.IAM: return 'IAM'
            elif model == _eggwrapper.Params.SMM: return 'SMM'
            elif model == _eggwrapper.Params.TPM: return 'TPM'
        elif key == 'TPM_proba': return self._params.get_TPMproba()
        elif key == 'TPM_param': return self._params.get_TPMparam()
        elif key == 'num_alleles': return self._params.get_K()
        elif key == 'rand_start': return self._params.get_random_start_allele()
        elif key == 'num_chrom': return self._num_chrom
        elif key == 'num_indiv': return self._num_indiv
        elif key == 'N': return self._N
        elif key == 'G': return self._G
        elif key == 's': return self._s
        elif key == 'site_pos': return self._site_pos
        elif key == 'site_weight': return self._site_weight
        elif key == 'migr_matrix': return self._migr_matrix
        elif key == 'trans_matrix': return self._trans_matrix
        elif key == 'seed': return self._seed
        elif key == 'events': return self._events
        elif key == 'max_iter': return self._params.get_max_iter()
        elif key == 'outgroup': return self._params.get_outgroup() if self._params.get_outgroup() != _eggwrapper.MISSING else None
        else: raise KeyError, 'invalid parameter name: {0}'.format(key)

    ####################################################################

    def __setitem__(self, key, value):

        if key == 'num_pop': raise ValueError, 'the number of populations is read-only (it must be set at construction time)'
        elif key == 'num_sites':
            if value < 0: raise ValueError, 'invalid value for `{0}`'.format(key)
            self._params.set_L(value)
            self._site_pos._num_values = value
            self._site_weight._num_values = value
            if value > 0: self._params.autoSitePos()
        elif key == 'recomb':
            if value < 0.0: raise ValueError, 'invalid value for `{0}`'.format(key)
            self._params.set_R(value)
        elif key == 'theta':
            if value < 0: raise ValueError, 'invalid value for `{0}`'.format(key)
            if value > 0.0 and self._params.get_fixed() > 0.0: raise ValueError, 'it is not allowed to set both theta and the number of mutations to non-zero'
            self._params.set_theta(value)
        elif key == 'num_mut':
            if value < 0: raise ValueError, 'invalid value for `{0}`'.format(key)
            if value > 0.0 and self._params.get_theta() > 0.0: raise ValueError, 'it is not allowed to set both theta and the number of mutations to non-zero'
            self._params.set_fixed(value)
        elif key == 'mut_model':
            if value == 'KAM': self._params.set_mutmodel(_eggwrapper.Params.KAM)
            elif value == 'IAM': self._params.set_mutmodel(_eggwrapper.Params.IAM)
            elif value == 'SMM': self._params.set_mutmodel(_eggwrapper.Params.SMM)
            elif value == 'TPM': self._params.set_mutmodel(_eggwrapper.Params.TPM)
            else: raise ValueError, 'invalid value for `{0}`'.format(key)
        elif key == 'TPM_proba':
            if value < 0 or value > 1: raise ValueError, 'invalid value for `{0}`'.format(key)
            self._params.set_TPMproba(value)
        elif key == 'TPM_param':
            if value < 0 or value > 1: raise ValueError, 'invalid value for `{0}`'.format(key)
            self._params.set_TPMparam(value)
        elif key == 'num_alleles':
            if value < 2: raise ValueError, 'invalid value for `{0}`'.format(key)
            self._params.set_K(value)
            self._trans_matrix._num_values = value
        elif key == 'rand_start':
            self._params.set_random_start_allele(value)
        elif key == 'max_iter':
            if value < 0: raise ValueError, 'maximum number of iterations cannot be negative'
            self._params.set_max_iter(value)
        elif key == 'seed':
            if value < 0: raise ValueError, 'random number generator seed cannot be negative'
            self._seed = value
            self._random.set_seed(value)
        elif key == 'num_chrom': self._num_chrom[:] = value
        elif key == 'num_indiv': self._num_indiv[:] = value
        elif key == 'N': self._N[:] = value
        elif key == 'G': self._G[:] = value
        elif key == 's': self._s[:] = value
        elif key == 'site_pos': self._site_pos[:] = value
        elif key == 'site_weight': self._site_weight[:] = value
        elif key == 'migr_matrix': self._migr_matrix.get_values(value)
        elif key == 'trans_matrix': self._trans_matrix.get_values(value)
        elif key == 'events': raise ValueError, 'cannot replace the list of events as a whole'
        elif key == 'outgroup':
            if value is None:
                self._params.set_outgroup(_eggwrapper.MISSING)
            elif value < 0:
                raise IndexError, 'invalid population index: {0}'.format(value)
            else:
                self._params.set_outgroup(value)
        else: raise KeyError, 'invalid parameter name: `{0}`'.format(key)

    ####################################################################

    def add_event(self, cat, T, **params):

        """
        Add an event to the list. See the documentation for the class
        :class:`.Simulator` for more details on what is expected as
        arguments.

        :param cat: category of the event.
        :param T: date of the event.
        :param params: all needed parameters for this category of
            events.
        """

        self._events.add(cat, T, **params)

########################################################################

class ParamList(object):

    """
    Emulator of :class:`list` for managing values for a parameter. Do
    most of what a list does except for adding and removing values. The
    methods :meth:`~list.reverse` and :meth:`~list.sort` are also not
    provided as they make little sense. Expressions involving the ``+``
    and ``*`` arithmetic operators are supported, but return genuine
    :class:`list` instances that are disconnected from the original
    parameter holder.

    In addition to the methods documented below, :class:`.ParamList`
    instances support the following operations (where ``items`` is a
    :class:`.ParamList` instance):

    +--------------------------+----------------------------------------------------------------+
    |Expression                | Action                                                         |
    +==========================+================================================================+
    |``len(items)``            | Number of items                                                |
    +--------------------------+----------------------------------------------------------------+
    |``items[i]``              | Get *i*\ th item                                               |
    +--------------------------+----------------------------------------------------------------+
    |``items[i:j]``            | Slice from *i* to *j*                                          |
    +--------------------------+----------------------------------------------------------------+
    |``items[i:j:k]``          | Slice from *i* to *j* with step *k*                            |
    +--------------------------+----------------------------------------------------------------+
    |``items[i] = value``      | Assign ``value`` to parameter ``key``                          |
    +--------------------------+----------------------------------------------------------------+
    |``items[i:j] = values``   | Replace a slice of values (with a sequence of the same length) |
    +--------------------------+----------------------------------------------------------------+
    |``items[i:j:k] = values`` | Replace a slice of values (with a sequence of the same length) |
    +--------------------------+----------------------------------------------------------------+
    |``for item in items``     | Iterates over items                                            |
    +--------------------------+----------------------------------------------------------------+
    |``reversed(items)``       | Return a reversed iterator                                     |
    +--------------------------+----------------------------------------------------------------+
    |``item in items``         | Check if ``item`` is present among items                       |
    +--------------------------+----------------------------------------------------------------+
    |``items + other``         | Same as ``list(items) + other`` (returns a :class:`list`)      |
    +--------------------------+----------------------------------------------------------------+
    |``items * n``             | Same as ``list(items) * n`` (returns a :class:`list`)          |
    +--------------------------+----------------------------------------------------------------+
    |``str(items)``            | Representation of the instance as a list                       |
    +--------------------------+----------------------------------------------------------------+

    The indexing operator (``[]``) supports negative indexes (to count
    from the end) and slices operators, just as for the built-in type
    :class:`list`. However, all operators (such as ``del``) or methods
    (such as :meth:`~.list.append`, :meth:`~.list.extend`,
    :meth:`~.list.remove`) that would change the length of the list are
    not available.
    """

    ####################################################################

    def __init__(self, getter, setter, num_values, check_item):
        self._getter = getter
        self._setter = setter
        self._num_values = num_values
        self._check_item = check_item

    ####################################################################

    def __getitem__(self, index):
        if isinstance(index, slice):
            return map(self.__getitem__, xrange(*index.indices(self._num_values)))
        else:
            if index < 0: index = self._num_values + index
            if index >= self._num_values: raise ValueError, 'list index out of range'
            return self._getter(index)

    ####################################################################

    def __setitem__(self, index, value):
        if isinstance(index, slice):
            rng = xrange(*index.indices(self._num_values))
            if len(rng) != len(value): raise ValueError, 'number of values is required to match number of indexes for a slice assignment'
            for i, v in zip(rng, value):
                if self._check_item(v) == False: raise ValueError, 'value out of range'
                self._setter(i, v)
        else:
            if index < 0: index = self._num_values + index
            if index >= self._num_values: raise ValueError, 'list index out of range'
            if self._check_item(value) == False: raise ValueError, 'value out of range'
            self._setter(index, value)

    ####################################################################

    def __iter__(self):
        for i in xrange(self._num_values):
            yield self._getter(i)

    ####################################################################

    def count(self, value):

        """
        Return the number of items that are equal to *value*.
        """

        cnt = 0
        for i in xrange(self._num_values):
            if self._getter(i) == value: cnt += 1
        return cnt

    ####################################################################

    def index(self, value, imin=0, imax=None):

        """
        Return the index of the first value matching *value*. The
        returned value is ``>=imin`` and ``<imax`` if they are provided.
        Raise a :exc:`~exceptions.ValueError` if the value is not found.
        """

        if imax == None: imax = self._num_values
        for i in xrange(imin, imax):
            if self._getter(i) == value: return i
        else:
            raise ValueError, '{0} is not in list'.format(value)

    ####################################################################

    def __add__(self, other):
        return list(self) + other

    def __radd__(self, other):
        return other + list(self)

    def __mul__(self, num):
        return list(self) * num

    def __rmul__(self, num):
        return num * list(self)

    ####################################################################

    def __contains__(self, value):
        for i in xrange(self._num_values):
            if self._getter(i) == value: return True
        else:
            return False

    ####################################################################

    def __len__(self):
        return self._num_values

    ####################################################################

    def __str__(self):
        return str(list(self))

    ####################################################################

    def __repr__(self):
        if self._num_values <= 10: return str(list(self))
        else: return '<list of {0} values>'.format(self._num_values)

    ####################################################################

    def __reversed__(self):
        for i in reversed(xrange(self._num_values)):
            yield self._getter(i)

########################################################################

class ParamMatrix(object):

    """
    Emulator of :class:`list` for managing values for a parameter as a
    matrix. Similar to :class:`.ParamList` except that it holds a
    matrix. It therefore implements less methods (no support for ``+``
    and ``*`` operators, no slices) and use double indexing system as in
    ``matrix[i, j] = x``, allowing both setting and getting values. The
    iterator, such as in ``for i in matrix``, flattens the matrix,
    returning all values from the first row, then those from the second
    tow, and so on (replacing the diagonal values by ``None``).
    ``len(matrix)`` returns the dimension (number of rows, which is the
    same as the number of columns).

    In addition to the methods documented below, :class:`.ParamMatrix`
    instances support the following operations (where ``matrix`` is a
    :class:`.ParamMatrix` instance):

    +--------------------------+----------------------------------------------------------------------+
    |Expression                | Action                                                               |
    +==========================+======================================================================+
    |``len(matrix)``           | Size of the matrix (as number of rows/colums)                        |
    +--------------------------+----------------------------------------------------------------------+
    |``matrix[i,j]``           | Get *j*\ th item of the *i*\ th row (``None`` for diagonal)          |
    +--------------------------+----------------------------------------------------------------------+
    |``matrix[i,j] = value``   | Assign ``value`` to *(i,j)*\ th item (no diagonal)                   |
    +--------------------------+----------------------------------------------------------------------+
    |``for item in matrix``    | Iterates over items (row-first *flat* iteration)                     |
    +--------------------------+----------------------------------------------------------------------+
    |``reversed(matrix)``      | Return a reversed iterator                                           |
    +--------------------------+----------------------------------------------------------------------+
    |``item in matrix``        | Check if ``item`` is present among items                             |
    +--------------------------+----------------------------------------------------------------------+
    |``str(matrix)``           | Representation of the instance as a nested list (including diagonal) |
    +--------------------------+----------------------------------------------------------------------+

    It is not possible to slice-set ranges of values in the matrix, but
    one can set all values at once from a compatible source with
    :meth:`~.ParamMatrix.get_values`.
    """

    ####################################################################

    def __init__(self, getter, setter, num_values, check_item, check_active, set_active):
        self._getter = getter
        self._setter = setter
        self._num_values = num_values
        self._check_item = check_item
        self._check_active = check_active
        self._set_active = set_active

    ####################################################################

    def get_values(self, other):

        """
        Get all values from *other*, which must be another
        :class:`.ParamMatrix` or a nested sequence (such as a
        :class:`list` or :class:`list`). The dimension of *other* must
        be the same as the current instance. Not that in the input
        object is not a :class:`.ParamMatrix`, any value it can have on
        the diagonal will be ignored (but they must be present to avoid
        shifting indexes).
        """

        if isinstance(other, ParamMatrix):
            if other._num_values != self._num_values: raise ValueError, '`other` must have the same dimension'

            if not other._check_active():
                if self._check_active(): self._set_active(False)
                return

            if not self._check_active():
                self._set_active(True)

            for i in xrange(self._num_values):
                for j in xrange(self._num_values):
                    if i != j: self._setter(i, j, other._getter(i, j))
                        # use of other._getter requires that number has been checked
        else:
            if not self._check_active(): self._set_active(True)
            if len(other) != self._num_values: raise ValueError, '`other` must have the same dimension'
            for i, row in enumerate(other):
                if len(row) != self._num_values: raise ValueError, '`other` must have the same dimension'
                for j, value in enumerate(row):
                    if i != j:
                        if self._check_item(value) == False: raise ValueError, 'value out of range'
                        self._setter(i, j, value)

    ####################################################################

    def __getitem__(self, (i, j)):
        if i < 0: i = self._num_values + i
        if j < 0: j = self._num_values + j
        if i >= self._num_values: raise ValueError, 'matrix index out of range'
        if j >= self._num_values: raise ValueError, 'matrix index out of range'
        if i == j: return None
        return self._getter(i, j)

    ####################################################################

    def __setitem__(self, (i, j), value):
        if i < 0: i = self._num_values + i
        if j < 0: j = self._num_values + j
        if i >= self._num_values: raise ValueError, 'matrix index out of range'
        if j >= self._num_values: raise ValueError, 'matrix index out of range'
        if i == j:
            if value is not None: raise ValueError, 'cannot set matrix diagonal to a value different than `None`'
            else: return
        if self._check_item(value) == False: raise ValueError, 'value out of range'
        if not self._check_active(): self._set_active(True)
        self._setter(i, j, value)

    ####################################################################

    def __iter__(self):
        for i in xrange(self._num_values):
            for j in xrange(self._num_values):
                if i == j: yield None
                else: yield self._getter(i, j)

    ####################################################################

    def count(self, value):

        """
        Return the number of items that are equal to *value*.
        """

        cnt = 0
        for i in xrange(self._num_values):
            for j in xrange(self._num_values):
                if i!=j and self._getter(i, j) == value: cnt += 1
        return cnt

    ####################################################################

    def index(self, value, imin=0, imax=None, jmin=0, jmax=None):

        """
        Return the ``(i, j)`` tuple of indexes of the first value
        (iterating first over rows and then over columns) matching
        *value*. The returned row index is >=imin and <imax and the
        returns column index is >=jmin and <jmax if they are provided.
        Raise a :exc:`~exceptions.ValueError` if the value is not found.
        The diagonal is not considered.
        """

        if imax == None: imax = self._num_values
        if jmax == None: jmax = self._num_values
        for i in xrange(imin, imax):
            for j in xrange(jmin, jmax):
                if i!=j and self._getter(i, j) == value: return (i, j)
        else:
            raise ValueError, '{0} is not in matrix'.format(value)

    ####################################################################

    def __contains__(self, value):
        for i in xrange(self._num_values):
            for j in xrange(self._num_values):
                if i!=j and self._getter(i, j) == value: return True
        else:
            return False

    ####################################################################

    def __len__(self):
        return self._num_values

    ####################################################################

    def __str__(self):
        return str([[self[i,j] for j in xrange(self._num_values)] for i in xrange(self._num_values)])

    ####################################################################

    def __repr__(self):
        if self._num_values <= 4: return str([[self[i,j] for j in xrange(self._num_values)] for i in xrange(self._num_values)])
        else: return '<matrix of {0}*{0} values>'.format(self._num_values)

    ####################################################################

    def __reversed__(self):
        for i in reversed(xrange(self._num_values)):
            for j in reversed(xrange(self._num_values)):
                if i == j: yield None
                else: yield self._getter(i, j)

########################################################################

class EventList(object):

    """
    Class storing the list of demographic events. Even if events appear
    as unsorted, they are internally sorted so that the user does not
    have to care about their order. This class has limited
    functionality: add events (but not removing them), accessing and
    modifying parameters.

    In addition to the methods documented below, :class:`.EventList`
    instances support the following operations (where ``events`` is an
    :class:`.EventList` instance):

    +--------------------------+----------------------------------------------------------------------+
    |Expression                | Action                                                               |
    +==========================+======================================================================+
    | ``len(events)``          | Number of events currently loaded                                    |
    +--------------------------+----------------------------------------------------------------------+
    | ``events[i]``            | Return a dictionary with parameters of the *i*\ th event (see note)  |
    +--------------------------+----------------------------------------------------------------------+
    | ``for event in events``  | Same as ``for i in range(len(events)): event = events[i]``           |
    +--------------------------+----------------------------------------------------------------------+
    | ``str(events)```         | Representation of the instance, roughly as a list                    |
    +--------------------------+----------------------------------------------------------------------+

    The string representation is such as a string at the first level,
    but each event is represented as an angle-bracketed delimited string
    containing comma-separated ``key=value`` pairs with the parameter as
    ``key`` and its value  as ``value`` (with an additional key,
    ``event_index``, providing the index of the event in the list).

    .. note::

        There is no way to modify content of the instance using the
        indexing operator ``events[i]``. One must use
        :meth:`~.EventList.update`.

    More information on the list of events and their parameters is
    available in the documentation for class :class:`.Simulator`.
    """

    ####################################################################

    _map_enum = {
        'size':          _eggwrapper.Event.change_N,
        'migr':          _eggwrapper.Event.change_M,
        'pair_migr':     _eggwrapper.Event.change_Mp,
        'growth':        _eggwrapper.Event.change_G,
        'selfing':       _eggwrapper.Event.change_s,
        'recombination': _eggwrapper.Event.change_R,
        'bottleneck':    _eggwrapper.Event.bottleneck,
        'admixture':     _eggwrapper.Event.admixture,
        'sample':        _eggwrapper.Event.delayed,
        'merge':         None
    }

    _required_parameters = {  # not: T is enforced in add's signature
        'size':          set(['N']),
        'migr':          set(['M']),
        'pair_migr':     set(['src', 'dst', 'M']),
        'growth':        set(['G']),
        'selfing':       set(['s']),
        'recombination': set(['R']),
        'bottleneck':    set(['S']),
        'admixture':     set(['src', 'dst', 'proba']),
        'sample':        set(['idx', 'label', 'num_chrom', 'num_indiv']),
        'merge':         set(['src', 'dst'])
    }

    _optional_parameters = {
        'size':       set(['idx']),
        'growth':     set(['idx']),
        'selfing':    set(['idx']),
        'bottleneck': set(['idx'])
    }

    _all_parameters = {}
    for i in _required_parameters:
        _all_parameters[i] = _required_parameters[i] | _optional_parameters.get(i, set())
        _all_parameters[i].add('T')

    ####################################################################

    def __init__(self, npop, add, clear):
        self._npop = npop
        self._add = add
        self._clear = clear
        self._events = []

    ####################################################################

    def replace(self, other):

        """
        Replace own list of events with the one in the
        :class:`.EventList` instance passed as *other*. The current list
        of events is dropped.
        """

        self.clear() # let self.clear() call self._clear() for security
        for event in other: self.add(** event)

    ####################################################################

    def clear(self):

        """
        Clear list of events.
        """

        self._clear() # do this first or the Event objects might be garbage-collected!
        del self._events[:]

    ####################################################################

    def __len__(self):
        return len(self._events)

    ####################################################################

    def __iter__(self):
        for params, changes in self._events:
            yield dict(params)

    ####################################################################

    def __str__(self):
        return ('[' + ', '.join([
                '<event_index={0};'.format(i)
                    + ';'.join(['{0}={1}'.format(k,v) for (k,v) in d.iteritems()]) + '>'
                        for i, (d, changes) in enumerate(self._events)]) + ']')

    ####################################################################

    def __repr__(self):
        return str(self)

    ####################################################################

    def __getitem__(self, i):
        if i >= len(self._events): raise ValueError, 'invalid event index'
        return dict(self._events[i][0])

    ####################################################################

    def add(self, cat, T, **params):

        """
        Add an event to the list. See the documentation for the class
        :class:`.Simulator` for more details on what is expected as
        arguments.

        :param cat: category of the event.
        :param T: date of the event.
        :param params: all needed parameters for this category of
            events.
        """

        # check that category is valid
        if cat not in self._all_parameters:
            raise ValueError, 'invalid event category: `{0}`'.format(cat)

        # check that all required parameters are present
        if not self._required_parameters[cat].issubset(params):
            raise ValueError, 'event `{0}` requires: {1}'.format(cat, ', '.join(self._required_parameters[cat] - set(params.keys())))

        # create the required backend objects
        if cat == 'merge':
            changes = [_eggwrapper.Event(self._map_enum['admixture'], T)]
            if (params['src'] < 0 or params['src'] >= self._npop or
                params['dst'] < 0 or params['dst'] >= self._npop or
                params['src'] == params['dst']):
                    raise ValueError, 'invalid population indexes provided for `merge` event'
            changes[0].set_index(params['src'])
            changes[0].set_dest(params['dst'])
            changes[0].set_param(1.0)
            for i in xrange(self._npop):
                if i != params['src']:
                    changes.append(_eggwrapper.Event(self._map_enum['pair_migr'], T))
                    changes[-1].set_param(0.0)
                    changes[-1].set_index(i)
                    changes[-1].set_dest(params['src'])
        else:
            changes = [_eggwrapper.Event(self._map_enum[cat], T)]

        # add the event to the internal list
        params['cat'] = cat
        params['T'] = T
        self._events.append((params, changes))
        params_update = dict(params)

        # use the update method to set parameters (delete the event in case of an error)
        del params_update['cat']
        del params_update['T']
        try:
            self.update(len(self._events) - 1, **params_update)
        except ValueError:
            del self._events[-1]
            raise

        # actually add changes to the lower-level Params
        for change in changes: self._add(change)

    ####################################################################

    def update(self, event_index, **params):

        """
        Modify any parameter from one of the event of the list. If an
        event's date is modified, sorting will be updated automatically.

        :param event_index: index of the event to modify (based on the
            order in which events have been specified with :meth:`.add`,
            which is the same order as events appear when representing
            the instance or iterating).
        :param params: keyword arguments specifying what parameters to
            modify. Only parameters that have to be changed should be
            specified.
        """

        # check index, get category
        if event_index < 0 or event_index >= len(self._events): raise IndexError, 'invalid event index'
        cat = self._events[event_index][0]['cat']
        change0 = self._events[event_index][1][0] # not used if complex change

        # initialize bool to perform ad-hoc tests
        test1 = False
        test2 = False

        # process complex events
        if cat == 'merge':
            if 'T' in params:
                for change in self._events[event_index][1]:
                    change.move(params['T'])
            if 'src' in params:
                test2 = True
                if params['src'] < 0 or params['src'] > self._npop: raise ValueError, 'event `{0}`, population index out of range'.format(cat)
                self._events[event_index][1][0].set_index(params['src'])
                for change in self._events[event_index][1][1:]:
                    change.set_dest(params['src'])

            if 'dst' in params:
                test2 = True
                if params['dst'] < 0 or params['dst'] > self._npop: raise ValueError, 'event `{0}`, population index out of range'.format(cat)
                self._events[event_index][1][0].set_dest(params['dst'])

        # process simple events
        else:

            # process all parameters and check+set them
            for key, value in params.iteritems():
                if key not in self._all_parameters[cat]:
                    raise ValueError, 'event `{0}` unknown parameter: {1}'.format(cat, key)

                if key == 'T':
                    if value < 0.0: raise ValueError, 'event `{0}`, date must be positive'.format(cat)
                    change0.move(value)
                elif key == 'N':
                    if value <= 0.0: raise ValueError, 'event `{0}`, size must be strictly positive'.format(cat)
                    change0.set_param(value)
                elif key == 'M':
                    if value < 0.0: raise ValueError, 'event `{0}`, size cannot be negative'.format(cat)
                    change0.set_param(value)
                elif key == 'G':
                    change0.set_param(value)
                elif key == 's':
                    if value < 0.0 or value > 1.0: raise ValueError, 'event `{0}`, selfing rate must be between 0 and 1'.format(cat)
                    change0.set_param(value)
                elif key == 'R':
                    if value < 0.0: raise ValueError, 'event `{0}`, recombination rate must be positive'.format(cat)
                    change0.set_param(value)
                elif key == 'S':
                    if value < 0.0: raise ValueError, 'event `{0}`, size must be positive'.format(cat)
                    change0.set_param(value)
                elif key == 'proba':
                    if value < 0.0 or value > 1.0: raise ValueError, 'event `{0}`, probability must be between 0 and 1'.format(cat)
                    change0.set_param(value)
                elif key == 'idx':
                    if value < 0 or value > self._npop: raise ValueError, 'event `{0}`, population index out of range'.format(cat)
                    change0.set_index(value)
                elif key == 'src':
                    test2 = True
                    if value < 0 or value > self._npop: raise ValueError, 'event `{0}`, population index out of range'.format(cat)
                    change0.set_index(value)
                elif key == 'dst':
                    test2 = True
                    if value < 0 or value > self._npop: raise ValueError, 'event `{0}`, population index out of range'.format(cat)
                    change0.set_dest(value)
                elif key == 'label':
                    if value < 0: raise ValueError, 'event `{0}`, label must be positive'.format(cat)
                    change0.set_dest(value)
                elif key == 'num_chrom':
                    test1 = True
                    if value < 0: raise ValueError, 'event `{0}`, number of samples must be positive'.format(cat)
                    change0.set_number1(value)
                elif key == 'num_indiv':
                    test1 = True
                    if value < 0: raise ValueError, 'event `{0}`, number of samples must be positive'.format(cat)
                    change0.set_number2(value)
                else:
                    raise RuntimeError, 'unexpected error'

        # perform ad-hoc tests
        if test1 and change0.get_number1() + change0.get_number2() == 0:
            raise ValueError, 'event `{0}`, there must be at least one sample'.format(cat)

        if test2 and change0.get_index() == change0.get_dest():
            raise ValueError, 'event `{0}`, cannot source/destination populations cannot be the same'.format(cat)

        # copy all parameters
        self._events[event_index][0].update(params)
