"""
This module contains the VCF parser.
"""

__license__ = """
    Copyright 2015-2018 Stephane De Mita, Mathieu Siol, Thomas Coudoux

    This file is part of EggLib.

    EggLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EggLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License@@@
    along with EggLib.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
from .. import _eggwrapper
from ..stats import _site

########################################################################

class VcfParser(object):

    """
    Read Variant Call Format (VCF)-formatted data for genomic
    polymorphism information from a file specified by name or from a
    strings.

    :param fname: name of a properly formatted VCF file. The header
        section will be processed upon instance creation, and lines will
        be read later, when the user iterates over the instance (or call
        :meth:`~.VcfParser.next`).

    :param allow_X: if ``True``, the characters ``X`` and ``x`` can be
        used instead of a base in alternate alleles. This is not allowed
        in the VCF specification but some software has actually used it.
        If ``X`` is allowed and one is found, the alternate type will be
        set to an *ad hoc* type and the corresponding allele string will
        be ``X`` (regardless of the original case).

    :param allow_gap: if ``True``, the gap symbol ``-`` is accepted as a
        valid base for the specification of both reference and alternate
        alleles. This is not allowed in the VCF specification which
        follows a different convention to represent insertions and
        deletions.

    :param find_index: if ``True``, the :class:`.VcfParser` instance will 
        find a index file linked to the Vcf data (by a default name), and
        will save it, in the attribut "_index" as `:class:VcfIndex` instance.
        If no file was found or the file found doesn't match with the Vcf
        data, the attribut "_index" will be "None".

    :param threshold_PL: this parameter controls how genotype calling (GT field) is
        performed from the PL (phred-scaled genotype likelihood) field.
        By default (``None``), this conversion is never performed. It is
        only performed if GT is not available and PL is available. The
        genotype with the lowest PL value is selected. The parameter gives the
        minimum acceptable gap between the best genotype and the second one.
        If the second genotype has a too good score, based on this parameter,
        the genotype is called unknown. The parameter must be at least 1.

    :param outgroup: list of samples (identified by their index in the
        sample list of this VCF file) that are supposed to belong to the
        outgroup. By default (``None``), same the as an empty list.

    :param outgroup_AA: a boolean indicating whether the ``AA`` field of
        the VCF file should be used as outgroup when populating :class:`Site`
        object. If the field is absent, the outgroup is set to missing
        data. Otherwise, its value is imported as outgroup. The outgroup is
        assumed to be homozygote for the allele indicated by ``AA``. When
        this option is used, it always sets it before all samples specified
        by *outgroup* for a given VCF file (there may be several VCF files
        loaded in the same :class:`Site` object in some situations).

    See the `description of the VCF format. <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41>`_

    There are two ways to process VCF data: one is using a static file
    that is iteratively parsed, using the standard constructor
    ``VcfParser(fname)`` and then iterate over lines in a ``for`` loop
    (alternatively, one can use ``VcfParser.next()`` directly), and the
    other way is to use the class factory method
    ``VcfParser.from_header(string)`` and then feed manually each line
    as a string using ``VcfParser.read_line(string)``.

    :class:`.VcfParser` instances are iterable (with support for the
    ``for`` statement and the :meth:`~.VcfParser.next` method) only if
    they are created with a file to process. Otherwise they must be fed
    line-by-line with the :meth:`.read_line` method. Every loop in a
    ``for`` loop or call to :meth:`~.VcfParser.next` or
    :meth:`.read_line` yields a ``(chromosome, position, num_all)``
    tuples that allows the user to determines if the variant is of
    interest. Note that the position is considered as an index and therefore
    has been decremented compared with the value found in the file.
     If the variant is of interest, the :class:`.VcfParser` object provides methods to
    extract all data for this variant (which can be time-consuming and
    should be restricted to pre-filtered lines to improve efficiency.

    .. versionadded:: 3.0.0
    """

    ####################################################################

    def __init__(self, fname, allow_X=False, allow_gap=False, find_index=True, 
                    threshold_PL=None, threshold_GL=None, outgroup=None,
                    outgroup_AA=False):

        self._parser = _eggwrapper.VcfParser()
        self._parser.open_file(fname)
        self._parser.allow_X(allow_X)
        self._parser.allow_gap(allow_gap)
        if outgroup is not None: self._set_outgroup(outgroup)
        if outgroup_AA: self._parser.set_outgroup_AA()
        self._fname = fname
        self._bed = None
        if find_index: self._index = self._find_index(fname)
        else: self._index= None

        if threshold_PL is None and threshold_GL is None:
            self._parser.set_threshold_PL(_eggwrapper.UNKNOWN)
            self._parser.set_threshold_GL(_eggwrapper.UNKNOWN)
        elif threshold_PL is not None and threshold_GL is None :
            if threshold_PL < 1: raise ValueError, 'threshold_PL must be at least 1'
            else: self._parser.set_threshold_PL(threshold_PL) 
        elif threshold_GL is not None and threshold_PL is None :
            if threshold_GL < 1: raise ValueError, 'threshold_GL must be at least 1'
            else: self._parser.set_threshold_GL(threshold_GL)            
        else:
            raise ValueError, 'Cannot use threshold_PL and threshold_GL at the same time'
   
    ####################################################################

    @classmethod
    def from_header(cls, string, allow_X=False, allow_gap=False,
        outgroup=None, outgroup_AA=False):

        """
        Create and return a new :class:`.VcfParser` instance reading
        the header passed as the *string* argument.

        :param string: single string including system-consistent line
            endings, the first line being the file format specification
            and the last line being the header line (starting with
            ``#CHROM``). This function allows leading and trailing white
            spaces (spaces, tabs, empty lines).
        :param allow_X: see class description.
        :param allow_gap: see class description.
        :param outgroup: list of samples (identified by their index in the
            sample list of this VCF file) that are supposed to belong to the
            outgroup. By default (``None``), same the as an empty list.
        :param outgroup_AA: a boolean indicating whether the ``AA`` field of
            the VCF file should be used as outgroup when populating :class:`Site`
            object. If the field is absent,
            the outgroup is set to missing data. Otherwise, its value is
            imported as outgroup. If data should be exported as genotypes,
            the outgroup is assumed to be homozygote for the allele indicated
            by ``AA``. Otherwise, only one outgroup genotype is created. When
            this option is used, it always sets the first outgroup sample.
        """

        obj = cls.__new__(cls)
        string = string.strip()
        obj._parser = _eggwrapper.VcfParser()
        obj._parser.read_header(string)
        obj._parser.allow_X(allow_X)
        obj._parser.allow_gap(allow_gap)
        if outgroup is not None: obj._set_outgroup(outgroup)
        if outgroup_AA: obj._parser.set_outgroup_AA()
        obj._fname = None
        obj._index= None
        return obj

    ####################################################################

    def _set_outgroup(self, indexes):
        for idx in indexes:
            if not isinstance(idx, int): raise TypeError, 'invalid type provided for outgroup index'
            if idx < 0: idx = self._parser.num_samples() + idx
            if idx < 0: raise IndexError, 'invalid index provided for outgroup index'
            self._parser.set_outgroup(idx)

    ####################################################################

    file_format = property(lambda self: self._parser.file_format(), doc='File format present in the read header.')
    num_info = property(lambda self: self._parser.num_info(), doc='Number of defined INFO fields.')
    num_format = property(lambda self: self._parser.num_format(), doc='Number of defined FORMAT fields.')
    num_filter = property(lambda self: self._parser.num_filter(), doc='Number of defined FILTER fields.')
    num_alt = property(lambda self: self._parser.num_alt(), doc='Number of defined ALT fields.')
    num_meta = property(lambda self: self._parser.num_meta(), doc='Number of defined META fields.')
    num_samples = property(lambda self: self._parser.num_samples(), doc='Number of samples read from header.')

    ####################################################################

    def get_sample(self, idx):

        """
        Get the name of a sample read from the header. The passed index
        must be smaller than :attr:`~.VcfParser.num_samples`.
        """

        if idx < 0 or idx >= self._parser.num_samples(): raise IndexError, 'sample index out of range'
        return self._parser.get_sample(idx)

    ####################################################################

    def get_info(self, idx):

        """
        Get data for a given INFO field defined in the VCF header. The
        passed index must be smaller than :attr:`.num_info`. Return a
        :class:`dict` containing the following data:

        * ``id``: ID string.
        * ``type``: one of ``"Integer"``, ``"Float"``, ``"Flag"``,
          ``"Character"``, and ``"String"``.
        * ``description``: description string.
        * ``number``: expected number of items. Special values are
          ``None`` (if undefined), ``"NUM_GENOTYPES"`` (number matching
          the number of genotypes for any particular variant),
          ``"NUM_ALTERNATE"`` (number matching the number of alternate
          alleles for any particular variant), and ``"NUM_ALLELES"``
          (number matching the number of alleles--including the
          reference--for any particular variant).
        * ``extra``: all extra qualifiers, presented as a :class:`list`
          of ``(key, value)`` tuples.
        """

        if idx < 0 or idx >= self._parser.num_info(): raise IndexError, 'INFO index out of range'
        info = self._parser.get_info(idx)

        n = info.get_number()
        if n == _eggwrapper.UNKNOWN: n = None
        elif n == _eggwrapper.NUM_ALTERNATE: n = 'NUM_ALTERNATE'
        elif n == _eggwrapper.NUM_GENOTYPES: n = 'NUM_GENOTYPES'
        elif n == _eggwrapper.NUM_POSSIBLE_ALLELES: n = 'NUM_ALLELES'

        t = info.get_type()
        if t == _eggwrapper.Info.Integer: t = 'Integer'
        elif t == _eggwrapper.Info.Float: t = 'Float'
        elif t == _eggwrapper.Info.Flag: t = 'Flag'
        elif t == _eggwrapper.Info.Character: t = 'Character'
        elif t == _eggwrapper.Info.String: t = 'String'
        else: raise ValueError, 'invalid VCF metainformation type found'

        return {
            'id': info.get_ID(),
            'description': info.get_description(),
            'number': n,
            'extra': [(info.get_extra_key(j), info.get_extra_value(j))
                                for j in xrange(info.get_num_extra())],
            'type': t
        }

    ####################################################################

    def get_format(self, idx):

        """
        Get data for a given FORMAT field defined in the VCF header. The
        passed index must be smaller than :attr:`.num_format`. Return a
        :class:`dict` containing the following data:

        * ``id``: ID string.
        * ``type``: one of ``"Integer"``, ``"Float"``, ``"Character"``,
          and ``"String"``.
        * ``description``: description string.
        * ``number``: expected number of items. Special values are
          ``None`` (if undefined), ``"NUM_GENOTYPES"`` (number matching
          the number of genotypes for any particular variant),
          ``"NUM_ALTERNATE"`` (number matching the number of alternate
          alleles for any particular variant), and ``"NUM_ALLELES"``
          (number matching the number of alleles--including the
          reference--for any particular variant).
        * ``extra``: all extra qualifiers, presented as a :class:`list`
          of ``(key, value)`` tuples.
        """

        if idx < 0 or idx >= self._parser.num_format(): raise IndexError, 'FORMAT index out of range'
        format_ = self._parser.get_format(idx)

        n = format_.get_number()
        if n == _eggwrapper.UNKNOWN: n = None
        elif n == _eggwrapper.NUM_ALTERNATE: n = 'NUM_ALTERNATE'
        elif n == _eggwrapper.NUM_GENOTYPES: n = 'NUM_GENOTYPES'
        elif n == _eggwrapper.NUM_POSSIBLE_ALLELES: n = 'NUM_ALLELES'

        t = format_.get_type()
        if t == _eggwrapper.Info.Integer: t = 'Integer'
        elif t == _eggwrapper.Info.Float: t = 'Float'
        elif t == _eggwrapper.Info.Character: t = 'Character'
        elif t == _eggwrapper.Info.String: t = 'String'
        else: raise ValueError, 'invalid VCF metainformation type found'

        return {
            'id': format_.get_ID(),
            'description': format_.get_description(),
            'number': n,
            'extra': [(format_.get_extra_key(j), format_.get_extra_value(j))
                                for j in xrange(format_.get_num_extra())],
            'type': t
        }

    ####################################################################

    def get_filter(self, idx):

        """
        Get data for a given FILTER field defined in the VCF header. The
        passed index must be smaller than :attr:`.num_filter`. Return a
        :class:`dict` containing the following data:

        * ``id``: ID string.
        * ``description``: description string.
        * ``extra``: all extra qualifiers, presented as a :class:`list`
          of ``(key, value)`` tuples.
        """

        if idx < 0 or idx >= self._parser.num_filter(): raise IndexError, 'FILTER index out of range'
        filter_ = self._parser.get_filter(idx)

        return {
            'id': filter_.get_ID(),
            'description': filter_.get_description(),
            'extra': [(filter_.get_extra_key(j), filter_.get_extra_value(j))
                                for j in xrange(filter_.get_num_extra())]
        }

    ####################################################################

    def get_alt(self, idx):

        """
        Get data for a given ALT field defined in the VCF header. The
        passed index must be smaller than :attr:`.num_alt`. Return a
        :class:`dict` containing the following data:

        * ``id``: ID string.
        * ``description``: description string.
        * ``extra``: all extra qualifiers, presented as a :class:`list`
          of ``(key, value)`` tuples.
        """

        if idx < 0 or idx >= self._parser.num_alt(): raise IndexError, 'ALT index out of range'
        alt = self._parser.get_alt(idx)

        return {
            'id': alt.get_ID(),
            'description': alt.get_description(),
            'extra': [(alt.get_extra_key(j), alt.get_extra_value(j))
                                for j in xrange(alt.get_num_extra())]
        }

    ####################################################################

    def get_meta(self, idx):

        """
        Get data for a given META field defined in the VCF header. The
        passed index must be smaller than :attr:`.num_meta`. Return a
        :class:`tuple` containing the key and the value of the META
        field.
        """

        if idx < 0 or idx >= self._parser.num_meta(): raise IndexError, 'META index out of range'
        meta = self._parser.get_meta(idx)

        return ( meta.get_key(), meta.get_value() )

    ####################################################################

    @property
    def good(self):

        """
        Tell if the file is good for reading (available valid stream,
        and not end of file).
        """

        return self._parser.good()

    ####################################################################

    def __iter__(self):
        return self

    def next(self):

        """
        Read one variant. Raise a :exc:`~exceptions.StopIteration`
        exception if no data is available.

        :return: The same as an iteration loop (see class description).
        """

        if not self._parser.good(): raise StopIteration

        self._parser.read()
        return self._parser.chromosome(), self._parser.position(), self._parser.num_alternate() + 1

    ####################################################################

    def read_line(self, string):

        """
        Read one variant from a user-provided single line. The string
        should contain a single line of VCF-formatted data (no header).
        All field specifications and sample information should be
        consistent with the information contained in the header that
        has been provided at creation-time to this instance (whichever
        it was read from a file or also provided as a string).

        :return: The same as an iteration loop (see class description).
        """

        self._parser.read_line(string)
        return self._parser.chromosome(), self._parser.position(), self._parser.num_alternate() + 1

    ####################################################################

    def last_variant(self):

        """
        Return a :class:`Variant` instance containing all data available
        for the last variant processed by this instance. It is required
        that a variant has been effectively processed.
        """

        return Variant._make(self)

    ####################################################################

    def get_genotypes(*parsers, **args):

        """
        get_genotypes([parser1[, parser2]], ..., dest=None)

        Process genotype data loaded into one or more
        :class:`.VcfParser` instances and return them as a single
        :class:`.Site` instance.

        There are two ways to use this method:

        1. As an instance method, to process a single parser, as in:
           ``parser.get_genotypes()`` (if ``parser`` is a
           :class:`.VcfParser` instance).

        2. As a class method, to process several parsers, as in:
           ``VcfParser.get_genotypes(parser1, parser2, parser3)`` (where
           ``parser1``, ``parser2`` and ``parser3`` are three
           :class:`.VcfParser` instances).

        This method requires that all parsers have been used to process
        valid data.

        :param parser: a :class:`.VcfParser` instance (only required if
            used as a class method). This argument can be repeated
            several times, but can not be passed as a keyword argument.

        :param dest: if specified, it must be a :class:`.Site`
            instance that will be recycled and used to place results.

        :return: A :class:`.Site` instance by default, or
            ``None`` if *dest* was specified.
        """

        # check parsers
        if len(parsers) == 0: raise ValueError, 'at least one `VcfParser` instance is required'
        for i, p in enumerate(parsers):
            if not isinstance(p, VcfParser): raise TypeError, 'a `VcfParser` instance is expected'
            if not p._parser.has_data(): raise ValueError, 'data must have been read from VCF parser'
            if not p._parser.has_GT(): raise ValueError, '`VcfParser` instance must have `GT` data'
        pl = parsers[0]._parser.ploidy()
        if pl < 1: raise ValueError, 'cannot import VCF parser: GT ploidy is 0'

        # create a new site instance or recycle if needed
        site = None
        for arg in args:
            if arg == 'dest':
                siteobj = args['dest']._obj
                siteobj.reset(pl)
            else:
                raise ValueError, 'invalid argument received: ' + arg

        if site is None:
            site = _site.Site()
            siteobj = site._obj
            siteobj.reset(pl)

        # process all parsers
        for parser in parsers:
            if pl != parser._parser.ploidy(): raise ValueError, 'inconsistent ploidy across VCF parsers'
            siteobj.process_vcf(
                        parser._parser,
                        0,                            # start
                        parser._parser.num_samples(), # stop
                        _eggwrapper.MAX)              # max_missing

        # create/recycle provided instance
        if 'dest' not in args:
            return site
  
    ####################################################################

    def _find_index(self, fname):

        # this helper looks for the index file corresponding to the provided VCF file name
        # returns a VcfIndex instance or None if not index file is found

        ifname=_eggwrapper._index_from_fname(fname)
        if(os.path.exists(ifname)):
            if(self._parser.file_end()==_eggwrapper.get_index_eof(ifname)): 
                index=_eggwrapper.VcfIndex()
                index.set_index_list(ifname)
                return index
            else:
                return None
        else:
            return None

    ####################################################################

    def read_index(self, fname):

        """
        An index file allows to increase the execution's speed of the progression of the :class:VcfParser` instance
        on the variants of a Vcf file. In Fact this binary file contains all start line position of each variant,
        linked to the chromosome and chromosomal position. That's allows to move the parser at a specific position
        according to a chromosome or/and chromosomal position given.

        This method allows to find, read and load as an VcfIndex an Index file linked to the current
        :class:VcfParser` instance.  The search of the index file, is done with a default name generated from the name
        of the Vcf file loaded in the current VcfParser.

        :param fname: name of a properly formatted Index file. The Index file must be linked to the vcf
            data of the :class:VcfParser` instance. In fact in the header on an index file, there are printed the "EOF"
            index of the vcf file behind it. The "EOF"  indexes, of the index file and of the Vcf file
            passed in the VcfParser, must match.
            EOF*: End Of File.
        """

        if(os.path.exists(fname)):
            if(self._parser.file_end()==_eggwrapper.get_index_eof(fname)): 
                self._index=_eggwrapper.VcfIndex()
                self._index.set_index_list(fname)
                self._fname=fname
            else:
                raise ValueError, 'your index file doesnt match with your object VcfParser'
        else:
            raise IOError, 'The file doesnt exist. Create a index_file from your vcf_file with the method "make index"'

    ####################################################################

    @property
    def num_index(self):

        """
        Gets the number of index loaded in the :class:VcfParser` instance.
        """

        if(self._index != None):
            return self._index.n_index()
        else:
            return None

    ####################################################################

    @property 
    def has_index(self):

        """
        Checks if the :class:VcfParser` instance has an index file linked loaded.
        """

        if(self._index != None):
               return True
        else:
            return False

    ####################################################################

    def get_position_index(self, chromosome, position):

        """
        Gets the index of a variant according to a chromosome and a chromosomal position.
        :param chromosome: This function expects the name of a desired chromosome (string)
        :param position: This function expects an chromosomal position, linked to the chromosome passed as argument
           according the Index file. (int)
        :beware: An index file must be loaded in the object "VcfIndex", with the function "set_index_file"
           before use the method "get_contigu_index	
        """

        if(self._index != None):
            return self._index.get_position_index(chromosome, position)
        else:
            raise ValueError, 'There are no index in the object VcfParser. Use the method read_index() to load an index file'

    ####################################################################

    def get_contig_index(self, chromosome):

        """
        Gets the index of the first variant of a chromosome inf the index file.
        :param chromosome: This function expects the name of a desired chromosome (string)
        Beware: an index file must be loaded in :class:`.VcfIndex` instance, with the function "set_index_file"
           before use the method "get_contigu_index	
        """

        if(self._index != None):
            return self._index.get_contigu_index(chromosome)
        else:
            raise ValueError, 'There are no index in the object VcfParser. Use the method read_index() to load an index file'


    ####################################################################

    def get_last_contig_index(self, chromosome):

        """
        Get the start stream position of the first variant linked to the last chromosome.
        :beware : A index file must be loaded in the object "VcfIndex", with the function "set_index_file"
           before use the method "get_contigu_index
        """

        if(self._index != None):
            return self._index.get_last_position_contig_index(chromosome)
        else:
            raise ValueError, 'There are no index in the object VcfParser. Use the method read_index() to load an index file'

    ####################################################################

    def get_variant_index(self, line):

        """
        Gets the index of a variant according to an indice (variant's position in the vcf file).
        :param line: This function expects a line number smaller of equal of the number
           of indice in the file index loaded. (int)
        :beware : A index file must be loaded in :class:`.VcfParser` instance, with the function
           "read_index"	before use the method 'get_indice_index'. 	
        """

        if(self._index == None):
            raise ValueError, 'There is no index in the object VcfParser. Use the method read_index() to load an index file'
        if (line > self._index.n_index()):
            raise ValueError, 'The indice passed as argument is greater than the number of indices in the index file loaded in the VcfParser' 
        return self._index.get_indice_index(line)

    ####################################################################"

    def make_index(self, fname=None, load=False):

        """
        This method allows to create a Index File.
        :param output: Name of the index file created.If the output argument is None,
           the created Index file will be named with a default name. The extension of
           an Index file is "vcfi".
        :param load: if ``True``, the VcfParser will loaded the index file created in the :class:`.VcfIndex` instance,
            else no data will be loaded after its creation. 
            call :meth:`.VcfWindow.next`
        """

        if(self._fname==None):
            raise NotImplementedError, 'No vcf file has been parsed in the "VCFParser" instance'

        if(fname==None):
            fname=_eggwrapper._index_from_fname(self._fname)

        _eggwrapper.index_file_from_vcf(self._parser ,fname)

        if(load==True):
            self.read_index(fname)
        return 0

   ####################################################################

    def goto(self, chromosome, position = None): #position=0):

        """
        This method allows to move  a :class:`~.VcfParser` instance  at a specific position
        in the vcf file according to a chromosome and a chromosomal position.
        :param chromosome: a string with a name of the chromosome desired
        :param position: a int with a chromosomal position linked to the chromosome.
        If only a chromosome is passed as argument or the position argument is None, the
        :class:`~.VcfParser` instance will moved at the position of the first variant of this 
        chromosome in the VcfFile.

        Load an index file in a :class:`~.VcfParser` instance allows to increase the execution's speed
        of the method :meth:`.VcfWindow.goto`
        """

        r_position = None

        if(self._index != None):
            r_index=0
            c_line=0
            if(position != None): #if(position!=0):
                r_index=self._index.get_position_index(chromosome, position)
                c_line=self._index.get_index_indice(r_index)-1
            else:
                r_index=self._index.get_contigu_index(chromosome)
                c_line=self._index.get_index_indice(r_index)-1

            self._parser.set_index(r_index)
            self._parser.set_currline(c_line)

            self._parser.read()
            r_position = self._parser.position()
            self._parser.unread()
            return r_position

        else:
            self._parser.rewind()
            while self._parser.good():
                self._parser.read_chromosome()
		if(position == None and self._parser.chromosome() == chromosome):
           		r_position = self._parser.position()
			self._parser.unread()
            		return r_position
		elif( position != None and self._parser.chromosome() == chromosome and self._parser.position() == position):
           		r_position = self._parser.position()
			self._parser.unread()
            		return r_position
			
	    if(r_position == None):
		raise ValueError, "The chromosome and chromosomal position desired , don't match with VcfParser"

    #################################################################### 

    def get_index(self):

        """
        This method allows to get the index of the current position of a :class:`~.VcfParser` instance.
        """

        return self._parser.get_index()
        
    #################################################################### 

    def unread(self):

        """
        This methods allows to unread the last variant read by the VcfParser with
        the method 'next'
        """

        self._parser.unread()


    ####################################################################   

    def rewind(self):

        """
        This method allows to move a :class:`~.VcfParser` instance to the first variant of the Vcffile loaded. 
        """

        self._parser.rewind()

    ####################################################################

    def slider(self, size, step, as_bp=True,
            start=0, stop=None, max_missing=0):

        """
        This method allows to create a sliding window from the current position.

        :param size: size of the sliding window (by default, in base pairs).
        :param step: increment of the sliding window (by default, in base pairs).
        :param as_bp: toggle size/step from base pairs to number of sites.
        :param start: start position of the sliding window.
        :param stop: stop position of the sliding window.
        :param max_missing: maximum number of missing alleles.

        :return: A :class:`~._vcf.VcfWindow` instance.
        """

        return VcfWindow._make(self, size, step, as_bp, start, stop, max_missing)

    ####################################################################   

    def get_bed_file(self,fname):

        """
        This method allows to load a ".bed" file in an object "VcfParser". 
        In a "VcfParser" a bed file is used to get somes variants in a ".vcf"
        file according to chromosomals coordinates saved in a ".bed" file loaded.
        """

        self._bed = _eggwrapper.BedParser()
        self._bed.get_bed_file(fname)


    ####################################################################

    def bed_slider(self, miss_pv, fill=False, start_pv=0, end_pv=None):
        
        """
        :param start_pv: if a subset of samples must be considered,
            index of the first sample to consider (by default,
            all samples are considered).
        :param end_pv: if a subset of samples must be considered,
            index of the last sample to consider (by default,
            all samples are considered).
        :param miss_pv: maximum number of missing alleles.
            If this proportion is processing is stopped and
            get_missing() returns max_missing + 1. Only missing
            data in this data set are considered.
            The starting point of the sliding window depends of the
            stream of the Vcf Parser passed in parameter.
        :param ploidy: the ploidy must be a strictly positive number.
        :param fill: if ``True``, the sliding window will read the first window,
            else no data will be loaded in the :class:`.VcfWindow` instance before 
            call :meth:`.VcfWindow.next`

        To use the method bed_slider, you must load a bed file with the method
        'get_bed_file' before.
        """

        if(self._bed == None): raise ValueError, 'There are no bed file in the object VcfParser. Use the method get_bed_file() to load an bed file'
        return BedWindow._make(self, miss_pv,  start_pv=0, end_pv=None,  fill=False)

########################################################################

class Variant(object):

    """
    Represent a single variant (one line from a VCF-formatted data
    file). The user cannot create instances of this class himself
    (instances are generated by :class:`VcfParser`) and instances are
    not modifiable in principle (however, some attributes provide
    mutable objects, as mentioned).

    .. note::

        The ``AA`` (ancestral allele), ``AN`` (allele number),
        ``AC`` (allele count), and ``AF`` (allele frequency) INFO fields
        as well as the ``GT`` (deduced genotype) FORMAT are
        automatically extracted if they are present in the the file and
        if their definition matches the format specification (meaning
        that they were not re-defined with different number/type) in
        the header. If present, they are available through the dedicated
        attributes :attr:`.AN`, :attr:`.AA`, :attr:`.AC`, :attr:`.AF`,
        :attr:`.GT`, :attr:`.ploidy` and :attr:`.GT_phased`. However,
        they are still available in the respective :attr:`.info` and
        :attr:`.samples` (sub)-dictionaries.
    """

    ####################################################################

    # the following are provided for comparison to content of variant.alternate_types
    alt_type_default = _eggwrapper.Default #: Explicit alternate allele (the string represents the nucleotide sequence of the allele).
    alt_type_referred = _eggwrapper.Referred #: Alternate allele referring to a pre-defined allele (the string provides the ID of the allele).
    alt_type_breakend = _eggwrapper.Breakend #: Alternate allele symbolizing a breakend (see VCF description for more details).

    ####################################################################

    def __init__(self):
        raise NotImplementedError, 'cannot create `Variant` instances directly'

    ####################################################################

    @staticmethod
    def _filter_missing(v, type_):
        if ((type_==int and v==_eggwrapper.MISSINGDATA) or
            (type_==float and v==_eggwrapper.UNDEF) or
            (type_==str and v[0]==_eggwrapper.MAXCHAR)): return None
        else: return v

    ####################################################################

    @classmethod
    def _make(cls, parser):

        if parser._parser.len_reference() == 0: raise ValueError, 'cannot generate `Variant` instance: no VCF line has been parsed (or the length of the reference allele is null)'

        obj = cls.__new__(cls)

        # get chromosome
        obj._chrom = parser._parser.chromosome()
        if obj._chrom == '': obj._chrom = None

        # get position
        obj._pos = parser._parser.position()
        if obj._pos == _eggwrapper.MISSING: obj._pos = None
        elif obj._pos == _eggwrapper.UNKNOWN: obj._pos = -1

        # get IDs
        obj._id = tuple(parser._parser.ID(i) for i in xrange(parser._parser.num_ID()))

        # get reference+alternate alleles
        obj._alleles = [parser._parser.reference()]
        if obj._alleles[0] == '': obj._alleles[0] = None
        obj._alternate_types = []
        obj._num_alternate = parser._parser.num_alternate()
        obj._num_alleles = obj._num_alternate + 1
        for i in xrange(obj._num_alternate):
            obj._alleles.append(parser._parser.alternate(i))
            obj._alternate_types.append(parser._parser.alternate_type(i))
        obj._alleles = tuple(obj._alleles)
        obj._alternate_types = tuple(obj._alternate_types)

        # get quality
        obj._quality = parser._parser.quality()
        if obj._quality == _eggwrapper.UNDEF: obj._quality = None

        # get test info
        n = parser._parser.num_failed_tests()
        if n == _eggwrapper.UNKNOWN: obj._failed_tests = None
        else: obj._failed_tests = tuple(parser._parser.failed_test(i) for i in xrange(n))

        # get info for the whole variant (site)
        obj._info = {}
        for i in xrange(parser._parser.num_FlagInfo()):
            item = parser._parser.FlagInfo(i)
            obj._info[item.get_ID()] = ()
        for i in xrange(parser._parser.num_IntegerInfo()):
            item = parser._parser.IntegerInfo(i)
            if item.get_expected_number() == 1:
                obj._info[item.get_ID()] = cls._filter_missing(item.item(0), int)
            else:
                obj._info[item.get_ID()] = tuple(cls._filter_missing(item.item(j), int) for j in xrange(item.num_items()))

        for i in xrange(parser._parser.num_FloatInfo()):
            item = parser._parser.FloatInfo(i)
            if item.get_expected_number() == 1:
                obj._info[item.get_ID()] = cls._filter_missing(item.item(0), float)
            else:
                obj._info[item.get_ID()] = tuple(cls._filter_missing(item.item(j), float) for j in xrange(item.num_items()))

        for i in xrange(parser._parser.num_CharacterInfo()):
            item = parser._parser.CharacterInfo(i)
            if item.get_expected_number() == 1:
                obj._info[item.get_ID()] = cls._filter_missing(item.item(0), str)
            else:
                obj._info[item.get_ID()] = tuple(cls._filter_missing(item.item(j), str) for j in xrange(item.num_items()))
        for i in xrange(parser._parser.num_StringInfo()):
            item = parser._parser.StringInfo(i)
            if item.get_expected_number() == 1:
                obj._info[item.get_ID()] = cls._filter_missing(item.item(0), str)
            else:
                obj._info[item.get_ID()] = tuple(cls._filter_missing(item.item(j), str) for j in xrange(item.num_items()))

        # get predefined AN/AA/AC/AF info fields if they are matching definition
        if parser._parser.has_AN(): obj._AN = parser._parser.AN()
        else: obj._AN = None
        if parser._parser.has_AA(): obj._AA = parser._parser.AA_string()
        else: obj._AA = None
        if parser._parser.has_AC(): obj._AC = tuple(parser._parser.AC(i) for i in xrange(parser._parser.num_AC()))
        else: obj._AC = None
        if parser._parser.has_AF(): obj._AF = tuple(parser._parser.AF(i) for i in xrange(parser._parser.num_AF()))
        else: obj._AF = None

        # pre-process the format column to allow match the ID's to the proper accessor method of SampleInfo
        fields = []
        for i in xrange(parser._parser.num_fields()):
            format_ = parser._parser.field(i)
            type_ = format_.get_type()
            if type_ == _eggwrapper.Info.String:
                f1 = _eggwrapper.SampleInfo.num_StringItems
                f2 = _eggwrapper.SampleInfo.StringItem
                t = str
            elif type_ == _eggwrapper.Info.Float:
                f1 = _eggwrapper.SampleInfo.num_FloatItems
                f2 = _eggwrapper.SampleInfo.FloatItem
                t = float
            elif type_ == _eggwrapper.Info.Integer:
                f1 = _eggwrapper.SampleInfo.num_IntegerItems
                f2 = _eggwrapper.SampleInfo.IntegerItem
                t = int
            elif type_ == _eggwrapper.Info.Character:
                f1 = _eggwrapper.SampleInfo.num_CharacterItems
                f2 = _eggwrapper.SampleInfo.CharacterItem
                t = str
            fields.append((
                format_.get_ID(), # the ID
                parser._parser.field_rank(i), # the index within range
                f1, # SampleInfo's method to get number of items
                f2, # SampleInfo's method to get a given item
                t   # type (required to convert missing data to None)
            ))
        obj._format_fields = frozenset(i[0] for i in fields)

        # get all SampleInfo data
        obj._num_samples = parser._parser.num_samples()
        obj._samples = []
        for sam in xrange(obj._num_samples):
            sample_info = parser._parser.sample_info(sam)
            sample_data = {}
            for ID, idx, f1, f2, t in fields:
                sample_data[ID] = (
                    None if f1(sample_info, idx)==0 else
                    tuple(cls._filter_missing(f2(sample_info, idx, i), t) for i in xrange(f1(sample_info, idx))))
                    # the above line get all items for a given FORMAT (using method pointers)
            obj._samples.append(sample_data)

        # get ploidy/num genotpes
        obj._ploidy = parser._parser.ploidy()
        obj._num_genotypes = parser._parser.num_genotypes()

        # get genotypes
        if parser._parser.has_GT():
            obj._gt = tuple(
                tuple((None if parser._parser.GT(i,j) == _eggwrapper.UNKNOWN
                            else obj._alleles[parser._parser.GT(i,j)])
                                    for j in xrange(obj._ploidy))
                                    for i in xrange(obj._num_samples))
            obj._gt_phased = tuple(parser._parser.GT_phased(i) for i in xrange(obj._num_samples))

            obj._gt_field = []
            gt_f = None
            for i in xrange(obj._num_samples):
                aa_f = []
                for j in xrange(obj._ploidy):
                    if parser._parser.GT(i,j) == _eggwrapper.UNKNOWN:
                        aa_f = None
                    else:
                        aa_f.append(parser._parser.GT(i,j))
                if(aa_f != None):
                    aa_fs = map(str, aa_f)  
                    if parser._parser.GT_phased(i) == False:
                        gt_f = '/'.join(aa_fs)
                    else : gt_f = '|'.join(aa_fs)
                else:
                    gt_f = '.'
                obj._gt_field.append(gt_f)
            obj._gt_field= tuple(obj._gt_field)

        else:
            obj._gt = None
            obj._gt_phased = None

        # get PL
        if parser._parser.has_PL():
            obj._pl = tuple(
                    tuple(parser._parser.PL(i, j) for j in xrange(obj._num_genotypes))
                            for i in xrange(obj._num_samples))
        else:
            obj._pl = None

        # get GL
        if parser._parser.has_GL():
            obj._gl = tuple(
                    tuple(parser._parser.GL(i, j) for j in xrange(obj._num_genotypes))
                            for i in xrange(obj._num_samples))
        else:
            obj._gl = None


        return obj

    ####################################################################
    # block of accessors below: they don't DO anything (it's only doc who makes is cumbersome)

    chromosome = property(lambda self: self._chrom, doc='Chromosome name (``None`` if missing).')
    position = property(lambda self: self._pos, doc='Position (as an index; first value is 0) (``None`` if missing).')
    ID = property(lambda self: self._id, doc='Tuple containing all IDs (even if just one or none).')
    num_alleles = property(lambda self: len(self._alleles), doc='Number of alleles (including the reference in all cases).')
    num_alternate = property(lambda self: self._num_alternate, doc='Number of alternate. Equal to :attr:`~.Variant.num_alleles` minus 1.')
    alleles = property(lambda self: self._alleles, doc='Variant alleles (the first is the reference and is not guaranteed to be present in samples), as a tuple.')
    alternate_types = property(lambda self: self._alternate_types, doc="""Alternate alleles types, as a tuple.
One value is provided for each alternate allele. The provided values are
integers whose values should always be compared to class attributes
:attr:`.alt_type_default`, :attr:`.alt_type_referred` and
:attr:`.alt_type_breakend`, as in (for the type of the first alternate
allele)::

    type_ = variant.alternate_types[0]
    if type_ == variant.alt_type_default:
        allele = variant.allele(0)rewind""")
    quality = property(lambda self: self._quality, doc='Variant quality (``None`` if missing).')
    failed_tests = property(lambda self: self._failed_tests, doc='Named of filters at which this variant failed, as a tuple (``None`` if no filters applied).')
    AA = property(lambda self: self._AA, doc='Value of the AA info field (``None`` if missing).')
    AN = property(lambda self: self._AN, doc='Value of the AN info field (``None`` if missing).')
    AC = property(lambda self: self._AC, doc='Value of the AC info field, as a tuple (``None`` if missing).')
    AF = property(lambda self: self._AF, doc='Value of the AF info field, as a tuple (``None`` if missing).')
    info = property(lambda self: self._info, doc="""Dictionary of INFO
fields for this variant. Keys are ID of INFO fields available for this
variant, and values are always a tuple of items. For flag INFO types,
the value is always an empty tuple.

.. note::

    This :class:`dict` is mutable, which enables the user to modify the
    data contained in the instance. Note that this will modify the data
    contained in this :class:`.Variant` instance, although not in the
    related :class:`.VcfParser` instance.""")
    format_fields = property(lambda self: self._format_fields, doc='Available FORMAT fields ID\'s available for each sample, as a :class:`frozenset` (empty if no sample data is available).')
    num_samples = property(lambda self: self._num_samples, doc='Number of samples (equivalent to ``len(Variant.samples)``).')
    samples = property(lambda self: self._samples, doc="""List of
information available for each sample (empty list if no samples are
defined). The list contains one :class:`dict` for each sample: keys of
these dictionary are FORMAT fields ID (the keys are always the same as
the content of :attr:`.format_fields`), and their values are tuples in
all cases.

.. note::

    This :class:`list` and the :class:`dict` instances it contains are
    all mutable, which enables the user to modify the data contained in
    the instance. Note that this will modify the data contained in this
    :class:`.Variant` instance, although not in the related
    :class:`.VcfParser` instance.""")
    ploidy = property(lambda self: self._ploidy, doc='Ploidy among genotypes (always 2 if GT is not available).')
    num_genotypes = property(lambda self: self._num_genotypes, doc='Number of genotypes.')
    GT_phased = property(lambda self: self._gt_phased, doc='Boolean indicating whether the genotype for each sample is phased (``None`` if GT is not available).')
    GT = property(lambda self: self._gt, doc="""Genotypes from GT
fields (only if this format field is available), provided as a tuple of
sub-tuples. The number of sub-tuples is equal to the number of samples
(:attr:`~.Variant.num_samples`). The number of items within each
sub-tuples is equal to the ploidy (:attr:`.ploidy`). These items are
allele expression (as found in :attr:`.alleles`), or ``None`` (for
missing values). This attribute is ``None`` if GT is not available.""")

    GT_vcf = property(lambda self: self._gt_field, doc = """GT field as written in a vcf file""")
    PL = property(lambda self: self._pl, doc = """PL values for all samples""")
    GL = property(lambda self: self._gl, doc = """GL values for all samples""")

########################################################################

class VcfWindow(object):
        
    """
    This class manages sliding Windows on :class:`.VcfParser` 
    instances. This class cannot be instanciated directly: instances are
    returned by the method :meth:`.VcfParser.slider`.
    
    :class:`~._vcf.VcfWindow` instances support the ``for`` statement).
    The iteration is concern the sites of the current sliding window. Each 
    itteration step returns a :class:`~.stats._site.Site` instance.
    :class:`~._vcf.VcfWindow` instances also support access using the ``[]`` operator
    as well as ``len()``.
    """

    ####################################################################

    def __init__(self):
        raise NotImplementedError, 'cannot create `VCFWindow` instances directly'

    ####################################################################

    @classmethod
    def _make(cls, parser, *args, **kwargs):
        obj = cls.__new__(cls)
        obj._wdw = _eggwrapper.VcfWindow()
        obj._parser = parser._parser
        obj.configure(*args, **kwargs)
        return obj

    ####################################################################

    def configure(self, size, step, as_bp, start, stop, max_missing):

        """
        Configure a sliding window (see :meth:`VcfParser.slider` for arguments)
        """

        if not self._parser.good(): raise ValueError, 'parser reached end of file'
        if self._parser.get_index() == self._parser.file_end(): raise ValueError, 'your parser reached the end of file, you cannot reconfigure your sliding window from this position. Use the method :meth:`.VcfParser.rewind` to be at the beginning of file or :meth:`.VcfParser.goto` to be at a higher position on the vcf file.'
        if max_missing > self._parser.num_samples(): raise ValueError, 'invalid value for `max_missing`'
        if stop is None: stop = _eggwrapper.UNKNOWN
        if start < 0: raise ValueError, 'invalid start value'
        if stop < start: raise ValueError, 'invalid stop value'

        ns = self._parser.num_samples()
        self._wdw.setup(self._parser, size, step, as_bp, start, stop, max_missing)

    ####################################################################

    @property
    def good(self):

        """
        Checks if :class:`~._vcf.VcfWindow` instance can continue on the 
        :class:`.VcfParser` instance
        """

        return self._wdw.good()

    ####################################################################

    num_sites = property(lambda self: self._wdw.num_sites(), doc='get the number of sites in the current :class:`~._vcf.VcfWindow` instance')
    chromosome = property(lambda self: self._wdw.chromosome(), doc='get the chromosome in the current :class:`~._vcf.VcfWindow` instance')
    num_samples = property(lambda self: self._num_samples(), doc='number of samples used by this window')

    @property
    def first_site(self):
        """ Position of the first actual site in the current window (``None`` if no sites) """
        if self._wdw.num_sites() == 0: return None
        else: return self._wdw.first_site_pos()

    @property
    def last_site(self):
        """ Position of the last actual site in the current window (``None`` if no sites) """
        if self._wdw.num_sites() == 0: return None
        else: return self._wdw.last_site_pos()

    @property
    def start_bound(self):
        """ First position of the current window (first site if the unit is not bp) """
        if self._wdw.first_pos() == _eggwrapper.UNKNOWN: return None
        else: return self._wdw.first_pos()

    @property
    def end_bound(self):
        """ Last position of the current window (last site if the unit is not bp) """
        if self._wdw.last_pos() == _eggwrapper.UNKNOWN: return None
        else: return self._wdw.last_pos()

    ####################################################################

    def next(self):

        """
        Allows to move forward the sliding window, from the initial position of the
        vcfparser to a chromosomial position passed as parameter 'stop' in the configuration
        method, until the last variant of the :class:`.VcfParser` instance or of 
        the chromosome read by the sliding window.
        This progression depends of the size window and the increment.	
        To start the sliding window at a specific position, use the method :meth:`.VcfParser.goto`
        method of the :class:`.VcfParser` before call the the method :meth:`.VcfParser.slider` 
        of the :class:`.VcfParser`
        """

        if self._wdw.good(): self._wdw.next_window()
        else: raise ValueError, 'sliding window has completed (check `good` attribute to avoid this error)'

    ####################################################################

    def __getitem__(self, idx):
        if idx >= self._wdw.num_sites() or idx < 0:
            raise ValueError, 'invalid site index: {0} is out of bounds [0, {1}['.format(idx, self._wdw.num_sites())
        else:
            return _site.Site._from_site_holder(self._wdw.get_site(idx).site())

    ####################################################################

    def __iter__(self):
        site = self._wdw.first_site()
        while site is not None:
            yield _site.Site._from_site_holder(site.site())
            site = site.next()

    def __len__(self):
        return self._wdw.num_sites()

########################################################################

class BedWindow(object):
   
    """
    This class allows to create a bed sliding Window on a :class:`.VcfParser` 
    instance. The creation of the bed sliding window is done according to a ".bed" file.
    In fact, this class allows to get a set of variants from a VcfParser instance and
    depending on the chromosomes coordinates of the ".bed" file linked.
 
    But this class cannot be called directly. To create a bed sliding window you must call
    the method :meth:`.VcfParser.slider` of the :class:`VcfParser` after a ".bed" file 
    has been loaded in the instance :class:`VcfParser`.
    
    :class:`~._vcf.BedWindow` instances are iterable by the special 
    :meth:`~._vcf.BedWindow.__iter__` (with support 'for' statement).
    The iteration is done on the sites of the current sliding window. Each 
    itteration return a :class:`.stats._site.Site` instance.
    :class:`~._vcf.BedWindow` instances are slices objects (usable by the next
    command: self[key]). The method :meth:`~._vcf.BedWindow.__getitem__` allows 
    the access by index. This method gets the site at given index. Negatives
    values are not allowed, else returns a 'ValueError'.
     
    """  

    ####################################################################

    def __init__(self):
        raise NotImplementedError, 'cannot create `BedWindow` instances directly'

    ####################################################################

    @classmethod
    def _make(cls, parser, miss_pv, start_pv=0, end_pv=None, fill=False):
        """
        Configure a bed sliding window.
        """

        if (parser._bed == None or parser._bed.n_bed_data() == 0): raise ValueError, 'You cannot use the class bEDWindow without load a bed file in the VcfParser'
        if(parser._fname == None): raise ValueError, 'No variant loaded, you cannot use cannot create a sliding window on the VCFPaser Instance'
        if not parser._parser.good(): raise ValueError,  'your parser reached the end of file, you cannot create a sliding window from this position. Use the method :meth:`.VcfParser.rewind` to be at the beginning of file or :meth:`.VcfParser.goto` to be at a higher position on the vcf file.'	
        if( miss_pv > parser.num_samples ): raise ValueError, "The value of the parameter 'miss_pv' cannot be greater than the number of samples"

        n_sample=parser._parser.num_samples()
        if start_pv < 0 or end_pv > n_sample: raise IndexError, 'invalid start index :arg:start_pv'
        if end_pv == None: end_pv = n_sample

        obj = cls.__new__(cls)
        obj._wbed = _eggwrapper.BedWindow()
        obj._parser = parser
        obj._wbed.configure(parser._parser, parser._bed, start_pv, end_pv, miss_pv)

        if(fill == True): obj.next()
        return obj

    def next(self):

        """
        This function allows the bed sliding window to position at a next chromosome coordinate according 
        to a ".bed" file linked.
        """

        if not self._wbed.good(): raise ValueError, 'you raised the last value of the vector of bed data loaded in the vcf parser' 
        self._wbed.next()


    ####################################################################	
    def at(self, indice):
        """
	This function allows the bed sliding window to position at a chromosome coordinate according 
	to a ".bed" and an indice. In fact chromosomes coordinates of a ".bed" file are saved in a vector
	and each chromosomes coordinates can be accessible by an indice smallest than the numbers total of 
	chromosomes coordinates in the vector. 
        """
        if indice > self._parser._bed.n_bed_data() or indice < 0: 
           raise ValueError, 'The indice passed as argument is greater than the number of bed data'
        self._wbed.at(indice)

    ####################################################################
    def __getitem__(self, indice):
        """
        Get the site at given index (used in expression: self[indice]). 

        """ 

        if indice > self.num_sites or indice < 0:
           raise ValueError, 'The indice passed as argument is greater than the number of site in the window'
        else:
           site=_site.Site._from_site_holder(self._wbed.get_site(indice))
           return site

    ####################################################################
    def __iter__(self):
        """
        Implementation of the iterable interface of the :class:`~._vcf.BedWindow`
	The iteration is done on the sites. 
        So the :class:`~._vcf.BedWindow` instances are iterable and can
	be use with support 'for' 
        """
        for i in xrange (self._wbed.get_n_site()):
            site=_site.Site._from_site_holder(self._wbed.get_site(i))
            yield site


    ####################################################################
    def __len__(self):
        """
        get the number of sites in the current :class:`~._vcf.BedWindow`
        instance
        this method is usable in this way: len(self)
        """
        return self._wbed.get_n_site()


    ####################################################################
    @property
    def good(self):
        """
        Checks if :class:`~._vcf.BedWindow` instance can continue to grow on the 
        :class:`.VcfParser` instance and according the number of data in the 
	BedParser reference.
        
        """
	if(self._wbed.good()):
	   return True
	else:
	   return False 
    ####################################################################

    start_variant= property(lambda self: self._wbed.get_start_w(), doc='start position of the current :class:`~._vcf.BedWindow` instance.')
    end_variant= property(lambda self: self._wbed.get_end_w(),doc='end position of the current :class:`~._vcf.BedWindow` instance.')
    num_sites = property(lambda self: self._wbed.get_n_site(), doc='get the number of sites in the current :class:`~._vcf.BedWindow` instance')
    chromosome = property(lambda self: self._wbed.get_chrom_w(), doc='get the chromosome in the current :class:`~._vcf.BedWindow` instance')

#########################################################################
