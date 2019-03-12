__license__ = """
    Copyright 2008-2011,2016 Stephane De Mita, Mathieu Siol

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

import re
from ..tools import _seq_manip
import StringIO

########################################################################

class GenBankFeatureLocation(object):

    """
    Hold the location of a GenBank feature. Supports various forms of
    location as defined in the GenBank format specification. The
    constructor contains a parser working from a GenBank-formatted
    string. By default, features are on the forward strand and segmented
    features are ranges (not orders).

    In addition to methods documented below, the following operations
    are supported for ``loc`` if it is a
    :class:`.GenBankFeatureLocation` instance:

    +-----------------------------+-----------------------------------------------------------------+
    |Expression                   | Action                                                          |
    +=============================+=================================================================+
    |``len(loc)``                 | Number of segments                                              |
    +-----------------------------+-----------------------------------------------------------------+
    |``loc[index]``               | Return the ``(fist, last)`` tuple for the corresponding segment |
    +-----------------------------+-----------------------------------------------------------------+
    |``for (first, last) in loc`` | Iterator over segments                                          |
    +-----------------------------+-----------------------------------------------------------------+
    |``str(params)``              | Generate a GenBank representation                               |
    +-----------------------------+-----------------------------------------------------------------+

    :class:`~.GenBankFeatureLocation` supports iteration and allows to
    iterate over ``(first,last)`` segments regardless of their types
    (for a single-base segment a position *position*, the tuple
    ``(position,position)`` is returned; similar 2-item tuples are
    returned for other types of segment as well).
    """

    ####################################################################

    def __init__(self, string=None):
        self._complement = False # True or False
        self._range = True # True or False
        self._type = '' # string of S|R|B|C for Single, Range, Between and Choice
        self._pos = [] # list of start/stop positions, always in increasing number
        self._ends = [] # list of boolean tuples (True if the segment is 5'/3' partial)

        if string!=None:
            self._parse(string)

    ####################################################################

    def __iter__(self):
        for pos in self._pos: yield pos
        
    ####################################################################

    def __len__(self):
        return len(self._pos)

    ####################################################################

    def __getitem__(self, index):
        return self._pos[index]

    ####################################################################

    def __str__(self):

        if not len(self._pos):
            raise RuntimeError, 'cannot format Genbank feature\'s location: no positions were loaded'

        string = []

        for i in range(len(self._pos)):
            if self._type[i]=='S':
                string.append( str(self._pos[i][0]  +1  ) )
            elif self._type[i]=='R' or self._type[i]=='C':
                item = ''
                if self._ends[i][0]:
                    item += '<'
                item += str(self._pos[i][0]  +1  )
                item += '.'
                if self._type[i]=='R':
                    item +='.'
                if self._ends[i][1]:
                    item += '>'
                item += str(self._pos[i][1]  +1  )
                string.append( item )
            elif self._type[i]=='B':
                string.append(  '%d^%d' %(self._pos[i][0] +1 , self._pos[i][1] +1) )
            else: raise RuntimeError, 'internal error in GenBankFeatureLocation.__str__'

        string = ','.join(string)

        if len(self._pos)>1:
            if self._range:
                string = 'join(%s)' %string
            else:
                string = 'order(%s)' %string

        if self._complement:
            string = 'complement(%s)' %string

        return string

    ####################################################################

    def copy(self):

        """
        Return a deep copy of the current instance.
        """

        gbfl = GenBankFeatureLocation()
        gbfl._complement = self._complement
        gbfl._range = self._range
        gbfl._type = str(self._type)
        gbfl._pos = [(i,j) for i,j in self._pos]
        gbfl._ends = [(i,j) for i,j in self._ends]

        return gbfl

    ####################################################################

    def _parse(self, string):

        # process top-level location groups

        # complement not incompatible with the rest
        if string[:11]=='complement(':
            if string[-1]!=')':
                raise IOError, 'invalid GenBank feature location: %s' %string
            self._complement = True
            string = string[11:-1]

        # range or order (only remove it)
        if string[:5] == 'join(':
            if string[-1]!=')':
                raise IOError, 'invalid GenBank feature location: %s' %string
            string = string[5:-1]
        elif string[:6] == 'order(':
            if string[-1]!=')':
                raise IOError, 'invalid GenBank feature location: %s' %string
            self._range = True
            string = string[6:-1]

        # extracts the segment(s)

        segments = string.split(',')

        for segment in segments:

            # position formats
            single = re.match('^(\d+)$', segment)
            rangeChoice   = re.match('^(\<?)(\d+)(\<?)(\.\.?)(\>?)(\d+)(\>?)$', segment)
            between = re.match('^(\d+)\^(\d+)$', segment)

            if single:
                try:
                    pos = int(single.group(1))  -1
                except ValueError:
                    raise IOError, 'invalid position in GenBankFeature location: %s' %segment
                self.add_single_base(pos)

            elif rangeChoice:

                # gets positions
                try:
                    start = int(rangeChoice.group(2))  -1
                    stop = int(rangeChoice.group(6))   -1
                except ValueError:
                    raise IOError, 'invalid position in GenBankFeature location: %s' %segment

                # gets partial marks
                Ls = rangeChoice.group(1)+rangeChoice.group(3)
                if Ls=='':
                    Lp = False
                elif Ls=='<':
                    Lp = True
                else:
                    raise IOError, 'invalid specification in GenBankFeature location: %s' %segment
                Rs = rangeChoice.group(5)+rangeChoice.group(7)
                if Rs=='':
                    Rp = False
                elif Rs=='>':
                    Rp = True
                else:
                    raise IOError, 'invalid specification in GenBankFeature location: %s' %segment
    
                # loads appropriate feature
                if rangeChoice.group(4)=='..':
                    self.add_base_range(start, stop, Lp, Rp)
                elif rangeChoice.group(4)=='.':
                    self.add_base_choice(start, stop, Lp, Rp)
                else:
                    raise RuntimeError, 'error in GenBankFeatureLocation constructor'

            elif between:
                try:
                    start = int(between.group(1))  -1
                    stop = int(between.group(2))   -1
                except ValueError:
                    raise IOError, 'invalid position in GenBankFeature location: %s' %segment
                if stop!=start+1:
                    raise IOError, 'invalid between-base feature specification: %s' %segment
                self.add_between_base(start)
            
            else:
                raise IOError, 'invalid GenBank feature segment positions: %s' %segment
        
    ####################################################################

    def set_complement(self):

        """
        Place the feature on the complement strand.
        """

        self._complement = True

    ####################################################################

    def set_forward(self):

        """
        Place the feature on the forward (not complement) strand, which
        is the default.
        """

        self._complement = False

    ####################################################################

    def is_complement(self):

        """
        ``True`` if the feature is on the complement strand.
        """

        return self._complement

    ####################################################################
    
    def as_order(self):

        """
        Define the feature as an order instead of a range.
        """

        self._range = False

    ####################################################################

    def as_range(self):

        """
        Define the features as a range, with is the default.
        """

        self._range = True

    ####################################################################

    def is_range(self):

        """
        ``True`` if the feature is a range (the default), ``False`` if
        it is an order.
        """

        return self._range

    ####################################################################

    def shift(self, offset):

        """
        Shift all positions according to the (positive of negative)
        argument.
        """

        self._pos = map(lambda x: [x[0]+offset, x[1]+offset], self._pos)

    ####################################################################

    def add_single_base(self, position):

        """
        Add a single-base segment to the feature. If no segments were
        entered previously, set the unique segment location. *position*
        must be an integer. All entered positions must be larger than
        any positions entered previously.
        """

        if not isinstance(position, int):
            raise TypeError, 'GenBankFeatureLocation positions must be of type int'
        if len(self._pos) and position<self._pos[-1][1]:
            raise ValueError, 'GenBankFeatureLocation positions must be entered in increasing order'

        self._pos.append((position,position))
        self._type += 'S'
        self._ends.append((False,False))

    ####################################################################

    def add_between_base(self, position):

        """
        Add a segment lying between two consecutive  bases. If no
        segments were entered previously, set the unique segment
        location. *position* must be an integer. The feature will be set
        between *position* and *position* + 1. If the feature is
        intended to be placed on the complement strand between
        positions, say, 1127 and 1128, one must use
        ``add_between_base(1127)`` in combination with
        ``set_complement()``. All entered positions must be larger than
        any positions entered previously.
        """

        if not isinstance(position, int):
            raise TypeError, 'GenBankFeatureLocation positions must be of type int'
        if len(self._pos) and position<self._pos[-1][1]:
            raise ValueError, 'GenBankFeatureLocation positions must be entered in increasing order'

        self._pos.append((position,position+1))
        self._type += 'B'
        self._ends.append((False,False))

    ####################################################################

    def add_base_range(self, first, last, left_partial=False, right_partial=False):

        """
        Add a base range the feature. If no segments were previously
        enter, set the unique segment location. *first* and *last* must
        be integers. The feature will be set between *first* and *last*
        positions, including both limits. If the feature is intended to
        be placed on the complement strand between positions, say, 1127
        and 1482, one must use ``add_base_range(1127,1482)`` in
        combination with ``set_complement()``. All entered positions
        must be larger than any positions entered previously and *last*
        must be larger than *first* (but can be equal). *left_partial*
        and/or *right_partial* must be set to ``True`` if, respectively,
        the real start of the segment lies 5' of *first* and/or the real
        end of the segment lies beyond *last* (relatively to the forward
        strand and consistently with the numbering system).
        """

        if not isinstance(first, int) or not isinstance(last, int):
            raise TypeError, 'GenBankFeatureLocation positions must be of type int'
        if len(self._pos) and first<self._pos[-1][1]:
            raise ValueError, 'GenBankFeatureLocation positions must be entered in increasing order'
        if (last<first):
            raise ValueError, 'GenBankFeatureLocation positions must be entered in increasing order'

        self._pos.append((first,last))
        self._type += 'R'
        self._ends.append(( bool(left_partial), bool(right_partial) ))

    ####################################################################

    def add_base_choice(self, first, last, left_partial=False, right_partial=False):

        """
        Add a segment corresponding to a single base chosen within a
        base range. If no segments were previously enter, set the unique
        segment location. *first* and *last* must be integers. The
        feature will be set between *first* and *last* positions,
        including both limits. If the feature is intended to be placed
        on the complement strand between positions, say, 1127 and 1482,
        one must use ``add_base_choice(1127,1482)`` in combination with
        ``set_complement()``. All entered positions must be larger than
        any positions entered previously and *last* must be strictly
        larger than *first*. *left_partial* and/or *right_partial* must
        be set to ``True`` if, respectively, the real start of the
        segment lies 5' of *first* and/or the real end of the segment
        lies beyond *last* (relatively to the forward strand and
        consistently with the numbering system).
        """

        if not isinstance(first, int) or not isinstance(last, int):
            raise TypeError, 'GenBankFeatureLocation positions must be of type int'
        if len(self._pos) and first<self._pos[-1][1]:
            raise ValueError, 'GenBankFeatureLocation positions must be entered in increasing order'
        if (last<first):
            raise ValueError, 'GenBankFeatureLocation positions must be entered in increasing order'
        if (first==last):
            raise ValueError, 'invalid use of GenBankFeatureLocation.add_base_choice: start=stop (use add_single_base instead)'

        self._pos.append((first,last))
        self._type += 'C'
        self._ends.append(( bool(left_partial), bool(right_partial) ))

    ####################################################################

    def rc(self, length):

        """
        Reverse the feature positions: positions are modified to be
        counted from the end. The *length* of the complete sequence must
        be passed.
        """

        self._complement = not self._complement
        self._pos = self._pos[::-1]
        self._ends = self._ends[::-1]
        self._type = self._type[::-1]
        for i in range(len(self._pos)):
            a,b = self._pos[i]
            self._pos[i] = (length-b-1, length-a-1)
            self._ends[i] = ( self._ends[i][1], self._ends[i][0] )

########################################################################

class GenBankFeature(object):
        
    """
    Instances of this class represent features associated to a GenBank
    instance. They should not be instantiated or used separatedly of a
    :class:`.GenBank` instance. The constructor creates an empty
    instance (although a :class:`.GenBank` instance must be passed as
    *parent*) and either :meth:`~.GenBankFeature.set` or :meth:`~.GenBankFeature.parse` must be used
    subsequently.

    :param parent: a :class:`.GenBank` instance to which the feature
        should be attached.

    In addition to methods documented below, the following operations
    are supported for ``feat`` if it is a :class:`.GenBankFeature`
    instance:

    +--------------+---------------------------------------+
    |Expression    | Action                                |
    +==============+=======================================+
    |``str(feat)`` | GenBank representation of the feature |
    +--------------+---------------------------------------+
    """

    ####################################################################

    def __init__(self, parent):
        self._parent = parent
        self._type = ''
        self._location = GenBankFeatureLocation()
        self._qualifiers = []

    ####################################################################

    def get_type(self):
        
        """
        Return the type string of the instance.
        """

        return self._type

    ####################################################################

    def qualifiers(self):

        """
        Return a dictionary with all qualifier values. This method
        cannot be used to change data within the instance. Note that
        changes of the returned dictionary don't affect data contained
        in the instance.

        .. versionchanged:: 2.1.0
           Meaning changed.
        """

        return dict(self._qualifiers)

    ####################################################################

    def add_qualifier(self, key, value):

        """
        Add a qualifier to the instance's qualifiers.
        """

        self._qualifiers.append((key, value))

    ####################################################################

    def update(self, feat_type, location, ** qualifiers):

        """
        Update feature information.

        :param feat_type: a string identifying the feature type (such as
            ``"gene"``, ``"CDS"``, ``"misc_feature"``, etc.). All
            strings are acceppted.
        :param location: a :class:`.GenBankFeatureLocation` instance
            giving the feature's location.
        :param qualifiers: other qualifiers must be passed as keyword
            arguments. It is not allowed to use ``"type"`` as a
            qualifier keyword.
        """

        if 'type' in qualifiers:
            raise ValueError, 'cannot use "type" as custom qualifier in `GenBankFeature`\'s constructor'

        self._qualifiers = [(i,qualifiers[i]) for i in qualifiers]
        self._location = location
        self._type = feat_type

    ####################################################################

    def parse(self, string):

        """
        Update feature information from information read in a
        GenBank-formatted string.
        """

        try:
            self._type = string.split()[0]
            locstring = string.split()[1]
            string = string.split()[2:]
        except IndexError:
            raise IOError, 'invalid GenBank feature string'
        while len(string) and string[0][0]!='/':
            locstring += string.pop(0)
        qualifiers=string

        # now we have the feature's position in a genuine string (locstring)
        self._location = GenBankFeatureLocation(locstring)

        # and the rest of qualifiers in a list (qualifiers)
        self._qualifiers = []
        if not len(qualifiers):
            return

        # fuses qualifiers where needed
        i=1
        while i<len(qualifiers):
            if qualifiers[i][0]!='/' or '=' not in qualifiers[i]:
                qualifiers[i-1]+=' '+qualifiers[i]
                del qualifiers[i]
            else:
                i+=1
        
        # processes the qualifiers
        for qualifier in qualifiers:
            if qualifier[0]!='/':
                raise IOError, 'invalid qualifier string: %s' %qualifier
            try:
                pos = qualifier.index('=')
            except ValueError:
                raise IOError, 'invalid qualifier string: %s' %qualifier
            key = qualifier[1:pos]
            value = qualifier[pos+1:]
            if key=='translation':
                value = ''.join(value.split())

            self._qualifiers.append((key,value))

    ####################################################################

    def get_sequence(self):

        """
        Return the string corresponding to this feature. If the
        positions pass beyond the end of the parent's sequence, a
        :exc:`~exceptions.RuntimeError` (and not an
        :exc:`~exceptions.IndexError`) is raised.
        """

        seq = ''
        for i,j in self._location:
            if j>len(self._parent):
                raise RuntimeError, 'GenBank feature exceeds sequence length'
            seq += self._parent._sequence[i:j+1]
        if self._location.is_complement(): return _seq_manip.rc(seq)
        else: return seq

    ####################################################################

    def get_start(self):
        
        """
        First position of the first (or unique) segment, in such a way
        that ``start()`` is always smaller than ``stop()``.
        """

        return self._location[0][0]

    ####################################################################

    def get_stop(self):

        """
        Last position of the last (or unique) segment, in such a way
        that ``start()`` is always smaller than ``stop()``.
        """

        return self._location[-1][1]

    ####################################################################

    def copy(self, genbank):
        
        """
        Return a copy of the current instance, connected to the
        :class:`.GenBank` instance *genbank*.
        """

        feature = GenBankFeature(genbank)
        feature.update(
            self._type, self._location.copy(), **dict(self._qualifiers))

        return feature

    ####################################################################

    def shift(self, offset):

        """
        Shift all positions according to the (positive of negative)
        argument.
        """

        self._location.shift(offset)

    ####################################################################

    def __str__(self):
        string = '     %s %s\n'%(self._type.ljust(15), str(self._location))
        string = self._parent._wrap(string, 21) + '\n'
        
        for key, value in self._qualifiers:
            if '/' in value:
                if value[0]+value[-1]!='""' and value[0]+value[-1]!='\'\'':
                    value = '"%s"' %value
            string+= self._parent._wrap('                     /%s=%s' %(key,value), 21)
            string+= '\n'

        return string

    ####################################################################

    def rc(self, length=None):

        """
        Reverse-complement the feature: apply it to the complement
        strand and reverse positions counting from the end. The *length*
        argument specifies the length of the complete sequence and is
        usually not required.
        """

        if length==None:
            length = len(self._parent)
        self._location.rc(length)

########################################################################

class GenBank(object):

    """
    This class represents a GenBank-formatted DNA sequence record.

    :param fname: input file name.
    :param string: GenBank-formatted string.

    Only one of the two arguments *fname* and *string* can be
    non-``None``. If both are ``None``, the constructor generates an
    empty instance with sequence of length 0. If *fname* is
    non-``None``, a GenBank record is read from the file with this name.
    If ``string`` is non-``None``, a GenBank record is read directly
    from this string. The following variables are read from the parsed
    input if present: *accession*, *definition*, *title*, *version*,
    *GI*, *keywords*, *source*, *references* (which is a list), *locus*
    and *others*. Their default value is ``None`` except for
    *references* and *others* for which default is an empty list.
    *source* is a (*description*, *species*, *taxonomy*) tuple. Each of
    *references* is a (*header*, *raw reference*) tuple and each of
    *others* is a (*key*, *raw*) tuple.

    In addition to methods documented below, the following operations
    are supported for ``gb`` if it is a :class:`.GenBank` instance:

    +-------------------+-------------------------------------------------+
    |Expression         | Action                                          |
    +===================+=================================================+
    |``len(gb)``        | Length of the sequence attached to this record  |
    +-------------------+-------------------------------------------------+
    |``str(gb)``        | GenBank representation of the record            |
    +-------------------+-------------------------------------------------+
    |``for feat in gb`` | Iterate over :class:`.GenBankFeature` instances |
    |                   | of this record                                  |
    +-------------------+-------------------------------------------------+

    """

    ####################################################################

    def __init__(self, fname=None, string=None):
        self._sequence = ''
        self._features = []
        self.accession = None
        self.definition = None
        self.title = None
        self.locus = None
        self.version = None
        self.GI = None
        self.source = (None,None,None)
        self.references = []
        self.keywords = None
        self.others = []

        if fname!=None and string!=None:
            raise ValueError, 'GenBank constructor expects at most one argument'

        stream = None

        if fname:
            stream = open(fname)
        if string:
            stream = StringIO.StringIO(string)

        if stream:
            self._parse(stream)

    ####################################################################

    def add_feature(self, feature):
        
        """
        Add a feature to the instance. The argument *feature* must be a
        well-formed :class:`.GenBankFeature` instance.
        """

        self._features.append(feature)

    ####################################################################

    def number_of_features(self):

        """
        Give the number of features contained in the instance.
        """

        return len(self._features)

    ####################################################################

    @property
    def sequence(self):

        """
        Sequence string (can be modified). Note that changing the
        record's string might obsolete the features (meaning that the
        setting an invalid sequence might cause the features to point to
        incorrect or out-of-bounds regions of the sequence).
        """

        return self._sequence

    @sequence.setter
    def sequence(self, string): self._sequence = string

    ####################################################################

    def __iter__(self):
        for feature in self._features:
            yield feature

    ####################################################################

    def extract(self, from_pos, to_pos):
        
        """
        Return a new :class:`.GenBank` instance representing a subset of
        the current instance, from position *from_pos*  to *to_pos*. All
        (and only) features that are completely included in the
        specified range are exported.
        """

        if from_pos < 0 or to_pos>=len(self):
            raise ValueError, 'invalid positions for GenBank extraction'

        gb = GenBank()
        gb._sequence = self._sequence[from_pos:to_pos+1]
        for feature in self:
            if feature.get_start() >= from_pos and feature.get_stop() <= to_pos:
                    clone = feature.copy(gb)
                    clone.shift(-from_pos)
                    gb._features.append( clone )
        return gb

    ####################################################################
    
    def _parse(self, stream):
        
        # identifies blocks marked by a capitalized word at the very
        # beginning of the line, and send them to the dynamic parser
        
        key = None  # current block key
        buffer = [] # block being read
        
        while True:

            line = stream.readline()
            if not len(line): break

            if line.strip()=='ORIGIN':
                match = re.match('(.+)()', line.strip())
            else:
                match = re.match('^([A-Z]+) (.+)', line)
            if match:
                if key:
                    self._parse_block(key, buffer)
                key = match.group(1)
                buffer = [ match.group(2) ]
                if key=='ORIGIN':
                    break
            else:
                buffer.append(line)
                
        if key!='ORIGIN':
            raise IOError, 'GenBank records lacks sequence'
        
        self._parse_sequence(stream)

    ####################################################################
    
    def _parse_block(self, key, content):
        
        # alias "the dynamic parser". Processes the different GenBank
        # blocks appropriately
        
        def merge(items):
            return ' '.join([' '.join(i.split()) for i in items])

        if key=='LOCUS':
            if len(content)>1:
                raise IOError, 'GenBank record exhibits an invalid LOCUS line'
            self.locus = content[0].strip()

        elif key=='DEFINITION':
            self.definition = merge(content)

        elif key=='ACCESSION':
            self.accession = merge(content)

        elif key=='VERSION':
            content = merge(content)
            match = re.match('([^ ]+) +GI\:(\d+)', content)
            if not match:
                raise IOError, 'Incorrect VERSION/GI line in GenBank record'
            self.version = match.group(1)
            self.GI = match.group(2)

        elif key=='KEYWORDS':
            self.keywords = merge(content)

        elif key=='SOURCE':
            if len(content)<3:
                raise IOError, 'Incorrect SOURCE section in GenBank record'
            else:
                source= content[0].strip()
                organism= content[1].strip()
                if organism[:8] != 'ORGANISM':
                    raise IOError, 'Incorrect SOURCE section in GenBank record'
                organism= organism[10:]
                taxonomy= merge(content[2:])
            self.source = source, organism, taxonomy

        elif key=='REFERENCE':
            if len(content)<2:
                raise IOError, 'Incorrect REFERENCE section in GenBank record'
            self.references.append((content[0].strip(), ''.join(content[1:])))

        elif key=='FEATURES':
            self._parse_features(content[1:])

        else:
            self.others.append(( key, content[0]+'\n'+''.join(content[1:]).rstrip() ))

    ####################################################################

    def _parse_features(self, lines):

        buffer = ''
        for line in lines:
            # new feature
            if not re.match('^ {8}', line):
                if len(buffer):
                    feature = GenBankFeature(self)
                    feature.parse(buffer)
                    self._features.append(feature)
                    buffer=''
            buffer+=line

        if len(buffer):
            feature = GenBankFeature(self)
            feature.parse(buffer)
            self._features.append(feature)

    ####################################################################

    def _parse_sequence(self, stream):

        self._sequence = ''

        while True:
            line= stream.readline()
            if not len(line):
                raise IOError, 'GenBank sequence doesn\'t terminated appropriately'
            line=line.strip()
            match = re.match(' *[0-9]+ ([A-Za-z ]+)', line)
            if match:
                fragment = ''.join( match.group(1).split())
                self._sequence += fragment.upper()
            elif line=='//':
                break
            else:
                raise IOError, 'invalid sequence line in GenBank record (reproduced below)\n%s' %line

    ####################################################################

    def __len__(self):
        return len(self._sequence)

    ####################################################################

    def write(self, fname):

        """
        Create a file named *fname* and write the formatted record in.
        """

        f = open(fname, 'w')
        try:
            self.write_stream(f)
        finally:
            f.close()

    ####################################################################

    def __str__(self):
        stream = StringIO.StringIO()
        self.write_stream(stream)
        return stream.getvalue()

    ####################################################################

    def write_stream(self, stream):

        """
        Writes the content of the instance as a Genbank-formatted
        string within the passed file (or file-compatible) stream.
        """

        try: stream.write('LOCUS       %s\n' %self.locus)
        except AttributeError: raise TypeError, '`stream` must be a file-type object'
        if self.definition!=None:
            stream.write(self._wrap('DEFINITION  %s' %self.definition, 12) + '\n')
        if self.accession!=None:
            stream.write('ACCESSION   %s\n' %self.accession)
        if self.version!=None or self.GI!=None:
            stream.write('VERSION     %s  GI:%s\n' %(self.version, self.GI))
        if self.keywords!=None:
            stream.write(self._wrap('KEYWORDS    %s' %self.keywords, 12)+ '\n')
        if (self.source[0]!=None or self.source[1]!=None or self.source[2]!=None):
            stream.write(self._wrap('SOURCE      %s' %str(self.source[0]), 12)+ '\n')
        if self.accession!=None:
            stream.write(self._wrap('  ORGANISM  %s' %str(self.source[1]), 12)+ '\n')
            stream.write(self._wrap('            %s' %str(self.source[2]), 12)+ '\n')

        for reference in self.references:
            stream.write('REFERENCE   %s\n%s' %reference)

        for other in self.others:
            # ignores "BASE COUNT"
            if other[0]=='BASE' and other[1].split()[0]=='COUNT':
                continue
            # otherwise writes
            stream.write('%s %s\n' %other)    

        stream.write('FEATURES             Location/Qualifiers\n')
        for feature in self:
            stream.write(str( feature ))

        stream.write('ORIGIN')
        
        if len(self) > 999999999:
            raise IOError, 'cannot export a GenBank instance with a sequence longest than 999999999 base pairs'

        c=0
        while c<len(self._sequence):
            stream.write('\n' +str(c+1).rjust(9))
            for i in range(6):
                stream.write(' '+self._sequence[c:c+10].lower())
                c+=10
                if c>=len(self._sequence):
                    break
        stream.write('\n//\n')

    ####################################################################

    _LINEWIDTH = 80
    _MAXBREAK = 40

    def _wrap(self, string, spacing):
        buffer = ''.join([i for i in string])
        space=0
        comma=0
        c=0
        res = ''
        while c<len(buffer):
            if buffer[c]==' ':
                space=c
            if buffer[c]==',':
                comma=c
            c+=1
            if c==(self._LINEWIDTH-1):
                if c<len(buffer) and buffer[c]==' ':
                    space = c
                if space>=self._MAXBREAK:
                    res+=buffer[:space]
                    res+= '\n'
                    buffer = ''.join([' ']*spacing) + buffer[space+1:]
                    c=0
                    space=0
                    comma=0
                elif comma>=self._MAXBREAK:
                    res+=buffer[:comma+1]
                    res+= '\n'
                    buffer = ''.join([' ']*spacing) + buffer[comma+1:]
                    c=0
                    space=0
                    comma=0
                else:
                    res+=buffer[:c]
                    res+= '\n'
                    buffer = ''.join([' ']*spacing) + buffer[c:]
                    c=0
                    comma=0
                    space=0
        res += buffer
        return res.rstrip()

    ####################################################################

    def rc(self):

        """
        Reverse-complement the instance (in place). All features
        positions and the sequence will be reverted and applied to the
        complementary strand. The features will be sorted in increasing
        start position (after reverting). This method should be applied
        only on genuine nucleotide sequences.
        """

        self._sequence= _seq_manip.rc(self._sequence)
        for feature in self: feature.rc(len(self))
        self._features.sort(lambda x,y: cmp(x.get_start(), y.get_start()))
