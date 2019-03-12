"""
This module contains the GFF3 parser.
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

from .. import _eggwrapper

########################################################################

class GFF3(object):

    """
    Read General Feature Format (GFF)-formatted (version 3) genome
    annotation data from a file specified by name or from a provided
    string.

    See the description of the GFF3 format `<http://www.sequenceontology.org/gff3.shtml>`_.

    This class supports segmented features but only if they are
    consecutive in the file. All features are loaded into memory and can
    be processed interactively.

    :param source: name of GFF3-formatted data file, or a GFF3-formatted
        string is *from_string* is ``True``.
    :param from_string: if ``True``, the first argument is a
        GFF3-formatted string; if ``False``, the first argument is a
        file name.
    :param liberal: if ``True``, support some violations of the GFF3
        format.
    :param only_genes: if ``True``, only index top-level features that
        have ``gene`` as type. This does not reduce memory usage, but it
        reduces access time.

    The ``liberal`` argument allows to support a few violations from the
    canonical GFF3 format. The current list of violations is:

    * CDS features may lack a phase.

    The list of supported violations may change in the future between
    minor versions of EggLib. It is recommended to use ``liberal=True``
    only for exploratory analyses. Otherwise, it is better to fix the
    original file so that it complies with the format.

    .. versionadded:: 3.0.0
    """

    _thinning_factor = 10 ** 6

    ####################################################################

    def __init__(self, source, from_string=False, liberal=False, only_genes=False):

        self._obj = _eggwrapper.GFF3()
        self._obj.liberal(liberal)

        if from_string: self._obj.parse_string(source)
        else: self._obj.parse(source)

        self._metadata = [(self._obj.metadata_key(i), self._obj.metadata_value(i))
                                    for i in xrange(self._obj.num_metadata())]

        self._idx_mapping = {}
        self._num_top_features = 0
        self._types = {}

        for i in xrange(self._obj.num_features()):
            feat = self._obj.feature(i)
            if feat.get_num_parents() > 0: continue
            ftype = feat.get_type()
            if only_genes and ftype != 'gene': continue
            self._num_top_features += 1
            seqid = feat.get_seqid()
            if seqid not in self._idx_mapping: self._idx_mapping[seqid] = {}
            start = feat.get_start(0)
            startkey = start // self._thinning_factor
            if ftype not in self._types: self._types[ftype] = 1
            else: self._types[ftype] += 1
            if startkey not in self._idx_mapping: self._idx_mapping[seqid][startkey] = []
            self._idx_mapping[seqid][startkey].append( (
                            i,
                            start,
                            feat.get_end(feat.get_num_fragments() - 1), 
                            feat.get_ID(),
                            ftype))

        self._feature_cache = {}

    ####################################################################

    @property
    def metadata(self):

        """
        Metadata of the imported file (information present in the file
        header). Metadata are available as a :class:`list` of
        ``(key, value)`` :func:`tuple` 's. It is not possible to replace
        the list by another object, but it is allowed to modify the
        returned object.
        """

        return self._metadata

    ####################################################################

    @property
    def seqid(self):

        """
        All seqid values present in the imported file, as a
        :class:`frozenset` instance.
        """

        return frozenset(self._idx_mapping)

    ####################################################################

    @property
    def types(self):

        """
        Number of top features of each of the encountered types. This
        does not include lower-level features (that is, those who have
        a parent). Value is a :class:`dict` (there is no good reaon to
        modify it).
        """

        return self._types

    ####################################################################

    @property
    def num_top_features(self):

        """
        Number of features that are directly accessible (that is, those
        who don't have a parent).
        """

        return self._num_top_features

    ####################################################################

    @property
    def num_tot_features(self):

        """
        Total number of features of the input file, included lower-level
        features and features that have not been indexed, if any.
        """

        return self._obj.num_features()

    ####################################################################

    def feature_iter(self, seqid, start, end, feat_type = None):

        """
        Return an iterator over features, returned as
        :class:`.GFF3Feature` instances. Only top features complying
        with arguments are considered.

        :param seqid: seqid identifier (only features associated with
            this seqid are yielded.
        :param start: start position on the considered seqid (only
            features whose start position is >= this value are yielded).
        :param end: end position on the considered seqid (only features
            whose end position is <= this value are yielded. One can use
            ``None`` to process features until the end of the seqid.
        :param feat_type: process only features whose type is equal to
            this value. By default, process all features within range.

        If *seqid* is not valid for this data file, raise a
        :exc:`~exceptions.ValueError`.
        """

        if seqid not in self._idx_mapping: raise ValueError, 'invalid `seqid`: "{0}"'.format(seqid)
        startkey = start // self._thinning_factor
        if end == None: endkey = max(self._idx_mapping[seqid])
        else: endkey = end // self._thinning_factor
        for key in self._idx_mapping[seqid]:
            if key >= startkey and key <= endkey:
                for item in self._idx_mapping[seqid][key]:
                    idx, item_start, item_end, item_id, item_type = item
                    if item_start >= start and (end == None or item_end <= end) and (feat_type == None or item[3] == feat_type):
                        if item_id not in self._feature_cache:
                            self._feature_cache[item_id] = GFF3Feature._create(self, self._obj.feature(idx))
                        yield self._feature_cache[item_id]

########################################################################

class GFF3Feature(object):

    """
    Provide information related to a given feature. Currently, instances
    cannot be created by the user and are read-only. Data are available
    as read-only properties; some of these are lists and can be modified
    but it does not make any sense to do so.
    """

    ####################################################################

    def __init__(self):
        raise NotImplementedError, '`GFF3Feature` instances cannot be created by the use'

    ####################################################################

    @classmethod
    def _create(cls, gff3, feat):

        obj = cls.__new__(cls)
        obj._gff3 = gff3

        obj._aliases = [feat.get_Alias(i) for i in xrange(feat.get_num_Alias())]

        obj._attributes = [(feat.get_attribute_key(i),
            [feat.get_attribute_value(i, j) for j in xrange(feat.get_num_items_attribute(i))])
                            for i in xrange(feat.get_num_attributes())]

        obj._dbxref = [feat.get_Dbxref(i) for i in xrange(feat.get_num_Dbxref())]

        obj._derives_from = feat.get_Derives_from()
        if obj._derives_from == '': obj._derives_from = None

        obj._pos = [(feat.get_start(i), feat.get_end(i)) for i in xrange(feat.get_num_fragments())]
        if len(obj._pos) == 0: raise ValueError, 'feature must have at least one fragment'

        obj._gap = feat.get_Gap()
        if obj._gap == '': obj._gap = None

        obj._ID = feat.get_ID()
        if obj._ID == '': obj._ID = None

        obj._is_circular = feat.get_Is_circular()

        obj._name = feat.get_Name()
        if obj._name == '': obj._name = None

        obj._note = [feat.get_Note(i) for i in xrange(feat.get_num_Note())]

        obj._ontology_terms = [feat.get_Ontology_term(i) for i in xrange(feat.get_num_Ontology_term())]

        obj._parents = [feat.get_Parent(i) for i in xrange(feat.get_num_Parent())]

        obj._parts = [feat.get_part(i) for i in xrange(feat.get_num_parts())]

        obj._phase = feat.get_phase()
        if obj._phase == feat.zero: obj._phase = 0
        elif obj._phase == feat.one: obj._phase = 1
        elif obj._phase == feat.two: obj._phase = 2
        elif obj._phase == feat.no_phase: obj._phase = None
        else: raise ValueError, 'invalid value for phase'

        obj._score = feat.get_score()
        if obj._score == _eggwrapper.UNDEF: obj._score = None

        obj._seqid = feat.get_seqid()
        if obj._seqid == '': obj._seqid = None

        obj._source = feat.get_source()
        if obj._source == '': obj._source = None

        obj._target = feat.get_Target()
        if obj._target == '': obj._target = None

        obj._type = feat.get_type()
        if obj._type == '': obj._type = None

        return obj

    ####################################################################

    @property
    def aliases(self):

        """
        List of Alias attributes.
        """

        return self._aliases

    ####################################################################

    @property
    def attributes(self):

        """
        List of non-predefined attributes, as ``(key, items)`` tuples,
        where ``items`` is itself a list.
        """

        return self._attributes

    ####################################################################

    @property
    def dbxref(self):

        """
        List of Dbxref attributes.
        """

        return self._dbxref

    ####################################################################

    @property
    def derives_from(self):

        """
        Value of the Derives_from attribute, or ``None`` if this
        attribute was not defined.
        """

        return self._derives_from

    ####################################################################

    @property
    def gap(self):

        """
        Value of the Gap attribute, or ``None`` if this attribute was
        not defined.
        """

        return self._gap

    ####################################################################

    @property
    def ID(self):

        """
        Value of the ID attribute, or ``None`` if this attribute was not
        defined.
        """

        return self._ID

    ####################################################################

    @property
    def is_circular(self):

        """
        Value of the Is_circular attribute, as a boolean.
        """

        return self._is_circular

    ####################################################################

    @property
    def name(self):

        """
        Value of the Name attribute, or ``None`` if this attribute was
        not defined.
        """

        return self._name

    ####################################################################

    @property
    def notes(self):

        """
        List of Note attributes.
        """

        return self._note

    ####################################################################

    @property
    def ontology_terms(self):

        """
        List of Ontology_term attributes.
        """

        return self._ontology_terms

    ####################################################################

    @property
    def phase(self):

        """
        Value of the Phase attribute. Possible values are: 0, 1, 2 and
        ``None`` (if undefined).
        """

        return self._phase

    ####################################################################

    @property
    def score(self):

        """
        Value of the Score attribute, or ``None`` if this attribute was
        not defined.
        """

        return self._score

    ####################################################################

    @property
    def seqid(self):

        """
        Value of seqid for this feature.
        """

        return self._seqid

    ####################################################################

    @property
    def source(self):

        """
        Source of this feature.
        """

        return self._source

    ####################################################################

    @property
    def target(self):

        """
        Value of The Target attribute, or ``None`` if this attribute was
        not defined. The target value is no processed and is provided as
        a single string.
        """

        return self._target

    ####################################################################

    @property
    def type(self):

        """
        Type of this feature.
        """

        return self._type

    ####################################################################

    @property
    def start(self):

        """
        Start position.
        """

        return self._pos[0][0]

    ####################################################################

    @property
    def end(self):

        """
        End position.
        """

        return self._pos[-1][1]

    ####################################################################

    @property
    def num_fragments(self):

        """
        Number of fragments.
        """

        return len(self._pos)

    ####################################################################

    @property
    def positions(self):

        """
        List of start/end positions of all fragments.
        """

        return self._pos

    ####################################################################

    @property
    def num_parents(self):

        """
        Number of parents.
        """

        return len(self._parents)

    ####################################################################

    @property
    def num_parts(self):

        """
        Number of parts (descending features).
        """

        return len(self._parts)

    ####################################################################

    def get_parent(self, idx):

        """
        Get a parent, as a :class:`.GFF3Feature` instance. The index
        should be within range.
        """

        if idx < 0 or idx >= len(self._parents): raise IndexError, 'feature parent index out of range'
        item_id = self._parents[idx].get_ID()
        if item_id not in self._gff3._feature_cache:
            self._gff3._feature_cache[item_id] = self._create(self._gff3, self._parents[idx])
        return self._gff3._feature_cache[item_id]

    ####################################################################

    def get_part(self, idx):

        """
        Get a part (descending feature), as a :class:`.GFF3Feature`
        instance. The index should be within range.
        """

        if idx < 0 or idx >= len(self._parts): raise IndexError, 'feature part index out of range'
        item_id = self._parts[idx].get_ID()
        if item_id not in self._gff3._feature_cache:
            self._gff3._feature_cache[item_id] = self._create(self._gff3, self._parts[idx])
        return self._gff3._feature_cache[item_id]
