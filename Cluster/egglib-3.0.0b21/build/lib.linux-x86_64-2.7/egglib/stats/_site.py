__license__ = """
    Copyright 2016-2018 Stephane De Mita, Mathieu Siol

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

def site_from_align(align, index, filtr, struct=None):

    """
    Import allelic and genotypic data from a position of the provided
    :class:`Align` instance. The *struct* argument allows to process
    only a subset of the samples, and also controls the genotypic
    struct and the ploidy. This means that the same :class:`.Structure`
    must be used again to further process the resulting :class:`.Site`
    instance.

    :param align: a :class:`.Align` instance.
    :param index: the index of a valid (not out of range) position of
        the alignment.
    :param filtr: a :class:`.Filter` instance providing the list of
        alleles to be considered as valid or missing data.
    :param struct: a :class:`.Structure` instance that will be used
        to group samples in individuals. Note that if the :class:`.Structure`
        describes only a subset of the samples of the alignment, only
        those samples will be included in the resulting :class:`.Site`
        instance.  By default, all sampled are processed in haploid
        individuals.

    :return: A new :class:`.Site` instance. The numbers of ingroup and
        outgroup items of this instance are defined by the :class:`.Structure`
        instance passed as *struct* and can be smaller than the number
        of samples of the original alignment.
    """

    obj = Site()
    obj.process_align(align, index, filtr, struct)
    return obj

def site_from_list(ingroup, outgroup, flat=False):

    """
    Import allelic and genotpyic data from provided lists. Input data
    are equivalent to the return value of :meth:`.Site.as_list`.
    *ingroup* and *outgroup* provide data for the ingroup and outgroup
    respectively. Either can be replaced by ``None`` which is equivalent
    to an empty list. They are supposed to be lists of tuples if *flat*
    is ``False``, and lists of integers otherwise, but lists and tuples
    can be replaced by other sequence types.

    :param ingroup: if *flat* is ``False`` (default): a list of
        tuples (one tuple per individual) which contain a constant
        number of allelic values. All allelic values are taken as
        valid alleles except ``None`` which represent missing data.
        Alleles can be any integers. The number of items of all
        individuals of both ingroup and outgroup is required to be
        constant. If *flat* is ``True``: a lists of integers (``None``
        for missing data) giving allelic values of all samples.

    :param outgroup: same format as *ingroup*.

    :param flat: determines if genotypic (default) or allelic data
        are provided.
    """

    obj = Site()
    obj.process_list(ingroup, outgroup, flat)
    return obj

def site_from_vcf(vcf, start=0, stop=None, flat=False):

    """
    Import allelic and genotypic data from a VCF parser. The VCF parser
    must have processed a variant and the variant is required to have
    genotypic data available as the GT format field. An exception is
    raised otherwise.

    .. warning::
        VCF genotypes are exported as allele indexes (0 for the
        reference allele). This function treats them as allele values,
        meaning that they might be shifted (0 is the first allele found,
        which is not necessarily the reference allele). Be aware of this
        fact when appending VCF data to data from other sources (including
        other VCF files) in the same site data, or when processing
        individual alleles.

    :param vcf: a :class:`.VcfParser` instance containing data. There
        must at least one sample and the GT format field must be
        available. It is not required to extract variant data manually.

    :param start: index of the first sample to process. Index is
        required to be within available bounds (it must be at least 0
        and smaller than the number of samples in the VCF data). Note
        that in a VCF file a sample corresponds to a genotype.

    :param stop: sample index at which processing must be stopped (this
        sample is not processed). Index is required to be within
        available bounds (if must be at least equal to *start* and not
        larger than the number of samples in the VCF data). Note that in
        a VCF file, a sample corresponds to a genotype.

    :param flat: ignore individuals and load data as haploid genotypes.
        By default, genotypes are imported based on the ploidy defined
        in the data.
    """

    obj = Site()
    obj.process_vcf(vcf, start, stop, flat)
    return obj

class Site(object):

    """
    Store allelic and genotypic data at a single site. Allelic values
    are represented by integers. In case of complex, such as structural,
    variation, one can use the index of variants as allelic values, and
    some objects can lack allelic value information (a number of alleles
    of 0 will be reported).
    
    Loading data is incremental: processed new data will be added
    to previously loaded data (requiring that ploidy matches). To avoid
    this, use :meth:`~.Site.reset`.
    """

    @classmethod
    def _from_site_holder(cls, obj):
        ret = object.__new__(cls)
        ret._obj = obj
        return ret

    def __init__(self):
        self._obj = _eggwrapper.SiteHolder(0)

    def reset(self):

        """
        Clear all data from the instance.
        """

        self._obj.reset(0)

    @property
    def ploidy(self):

        """
        Current value of the ploidy (to change it, one needs to reset
        the instance).
        """

        return self._obj.get_ploidy()

    @property
    def ns(self):

        """
        Total number of samples (alleles) in the ingroup.
        """

        return self._obj.get_ploidy() * self._obj.get_ning()

    @property
    def no(self):

        """
        Total number of samples (alleles) in the outgroup.
        """

        return self._obj.get_ploidy() * self._obj.get_nout()

    @property
    def num_ingroup(self):

        """
        Total number of individuals in the ingroup.
        """

        return self._obj.get_ning()

    @property
    def num_outgroup(self):

        """
        Total number of individuals in the outroup.
        """

        return self._obj.get_nout()

    @property
    def num_missing(self):

        """
        Number of missing data (expressed in number of samples).
        """

        return self._obj.get_tot_missing()

    @property
    def num_missing_ingroup(self):

        """
        Number of missing data (expressed in number of samples) in the
        ingroup only.
        """

        return self._obj.get_missing_ing()

    @property
    def num_missing_outgroup(self):

        """
        Number of missing data (expressed in number of samples) in the
        outgroup only.
        """

        return self._obj.get_missing_otg()

    @property
    def num_alleles(self):

        """
        Number of alleles.
        """

        return self._obj.get_nall()

    def alleles(self):

        """
        Generate the list of alleles
        """

        return [self._obj.get_allele(i) for i in xrange(self._obj.get_nall())]

    def as_list(self, flat=False, skip_outgroup=False):

        """
        Generate one or two lists containing data from the instance.

        :param flat: flatten return lists and skip the individual
            level (by default, return lists are a list of tuples which
            represent individuals).

        :param skip_outgroup: return only data for the ingroup (by
            default, return two lists for respectively ingroup and
            outgroup data).

        :return: lists of allelic index (or ``None`` for missing data).
            the return value is either a single or two lists (based on
            the value of *skip_outgroup*), and the list or lists
            contain tuples representing individuals unless *flat* is
            ``True``.
        """

        if skip_outgroup: return self._as_list_helper(flat, self._obj.get_i, self._obj.get_ning())
        else: return (
            self._as_list_helper(flat, self._obj.get_i, self._obj.get_ning()),
            self._as_list_helper(flat, self._obj.get_o, self._obj.get_nout()))

    def _as_list_helper(self, flat, getter, num):
        if flat: return [getter(i, j) if getter(i, j) != _eggwrapper.MISSING else None for i in xrange(num) for j in xrange(self._obj.get_ploidy())]
        else: return [tuple(getter(i, j) if getter(i, j) != _eggwrapper.MISSING else None for j in xrange(self._obj.get_ploidy())) for i in xrange(num)]

    def ingroup(self, indiv, chrom=None):

        """ingroup(indiv[, chrom])

        Get a given ingroup genotype or allele, as allele indexes.

        :param indiv: individual index.
        :param chrom: chromosome index.

        If *chrom* is omitted, return the genotype as a tuple. Otherwise,
        return the specific allele. Setting *chrom* to ``None`` is the
        same as omitting it. *chrom* can be safely omitted for haploid
        data.
        """

        if indiv < 0 or indiv >= self._obj.get_ning():
            raise IndexError, 'invalid individual index'
        if chrom is None:
            return tuple(self._obj.get_i(indiv, j) if self._obj.get_i(indiv, j) != _eggwrapper.MISSING
                                      else None
                                        for j in xrange(self._obj.get_ploidy()))
        else:
            if chrom < 0 or chrom >= self._obj.get_ploidy():
                raise IndexError, 'allele index is inconsistent with ploidy'
            return self._obj.get_i(indiv, chrom) if self._obj.get_i(indiv, chrom) != _eggwrapper.MISSING else None

    def outgroup(self, indiv, chrom=None):

        """outroup(indiv[, chrom])

        Get a given outgroup genotype or allele, as allele indexes.

        :param indiv: individual index.
        :param chrom: chromosome index.

        If *chrom* is omitted, return the genotype as a tuple. Otherwise,
        return the specific allele. Setting *chrom* to ``None`` is the
        same as omitting it. *chrom* can be safely omitted for haploid
        data.
        """

        if indiv < 0 or indiv >= self._obj.get_nout():
            raise IndexError, 'invalid individual index'
        if chrom is None:
            return tuple(self._obj.get_o(indiv, j) if self._obj.get_o(indiv, j) != _eggwrapper.MISSING
                                      else None
                                        for j in xrange(self._obj.get_ploidy()))
        else:
            if chrom < 0 or chrom >= self._obj.get_ploidy():
                raise IndexError, 'allele index is inconsistent with ploidy'
            return self._obj.get_o(indiv, chrom) if self._obj.get_o(indiv, chrom) != _eggwrapper.MISSING else None

    def process_align(self, align, index, filtr, struct=None, reset=True):

        """
        Import data from the provided :class:`.Align` to the data
        currently held by the instance. Arguments are identical to the
        function :func:`.site_from_align`, expect *reset*.

        :param reset: if ``True``, reset the instance as if newly
            created. If ``False``, append the data to current data, if any.

        If *reset* is ``False`` and this instance currently holds data,
        the ploidy defined by the *struct* argument is required to match
        the current value. If *struct* is ``None``, the implied ploidy
        is 1 and is still required to match.
        """

        if struct == None:
            ploidy = 1
        else:
            struct = struct._obj
            if struct.get_ns_req() > align.ns: raise ValueError, 'structure does not match alignment'
            if struct.get_no_req() > align.no: raise ValueError, 'structure does not match alignment'
            ploidy = struct.get_ploidy()

        if reset or self._obj.get_ploidy() == 0: self._obj.reset(ploidy)
        elif ploidy != self._obj.get_ploidy():
            raise ValueError, 'ploidy of loaded data does not match previously current value'

        if index >= align.ls: raise IndexError, 'invalid site index'

        self._obj.process_align(align._obj, index, struct, filtr._obj, _eggwrapper.MAX, False)

    def process_list(self, ingroup, outgroup, flat=False, reset=True):

        """
        Import data from the provided lists to the data currently held
        by the instance. Arguments are identical to the function
        :func:`.site_from_list`, expect *reset*.

        :param reset: if ``True``, reset the instance as if newly
            created. If ``False``, append the data to current data, if any.

        If *reset* is ``False`` and this instance currently holds data,
        the ploidy defined by the
        input data is required to match the current value. If *flat* is
        ``True``, the implied ploidy is 1 and is still required to match.
        """

        if ingroup is None: ingroup = []
        if outgroup is None: outgroup = []
        if len(ingroup) == 0 and len(outgroup) == 0: 
            self._obj.reset()
            return # still was reset
        
        if flat == True:
            ploidy = 1
        else:
            ploidy = set(map(len, ingroup)) | set(map(len, outgroup))
            if len(ploidy) != 1: raise ValueError, 'ploidy of loaded data is not consistent'
            ploidy = ploidy.pop()

        if reset or self._obj.get_ploidy() == 0: self._obj.reset(ploidy)
        elif ploidy != self._obj.get_ploidy():
            raise ValueError, 'ploidy of loaded data does not match previously current value'
        
        if flat == True: _iter_ = lambda x: ([i] for i in x)
        else: _iter_ = iter

        cur = self._obj.get_ning()
        self._obj.add_ing(len(ingroup))
        for i, v in enumerate(_iter_(ingroup)):
            for j, k in enumerate(v):
                if k is None: k = _eggwrapper.MISSINGDATA
                if isinstance(k, basestring): k = ord(k)
                self._obj.load_ing(cur+i, j, k)

        cur = self._obj.get_nout()
        self._obj.add_otg(len(outgroup))
        for i, v in enumerate(_iter_(outgroup)):
            for j, k in enumerate(v):
                if k is None: k = _eggwrapper.MISSINGDATA
                if isinstance(k, basestring): k = ord(k)
                self._obj.load_otg(cur+i, j, k)

    def process_vcf(self, vcf, start=0, stop=None, flat=False, reset=True):

        """
        Import data from the provided VCF parser to the data currently
        held by the instance. Arguments are identical to the function
        :func:`.site_from_vcf`, expect *reset*.

        :param reset: if ``True``, reset the instance as if newly
            created. If ``False``, append the data to current data, if any.

        If *reset* is ``False`` and this instance currently holds data,
        the ploidy defined by the
        input data is required to match the current value. If *flat* is
        ``True``, the implied ploidy is 1 and is still required to match.

        .. warning::
            VCF genotypes are exported as allele indexes (0 for the
            reference allele). This function treats them as allele values,
            meaning that they might be shifted (0 is the first allele found,
            which is not necessarily the reference allele). Be aware of this
            fact when appending VCF data to data from other sources (including
            other VCF files) in the same site data, or when processing
            individual alleles.
        """

        if vcf._parser.has_data() == False: raise ValueError, 'data must have been read from VCF parser'
        if vcf._parser.has_GT() == False: raise ValueError, 'VCF parser does not have GT data'

        if flat == True: ploidy = 1
        else:
            ploidy = vcf._parser.ploidy()

        if reset or self._obj.get_ploidy() == 0: self._obj.reset(ploidy)
        elif ploidy != self._obj.get_ploidy():
            raise ValueError, 'ploidy of loaded data does not match previously current value'

        n = vcf._parser.num_samples()
        if start < 0 or start >= n: raise IndexError, 'invalid start index'
        if stop == None: stop = n
        elif stop < start or stop > n: raise IndexError, 'invalid stop index'

        self._obj.process_vcf(vcf._parser, flat, start, stop, _eggwrapper.MAX)

class Iterator(object):

    """
    Generate a iterator over sites of an alignment, the sites with too
    many missing data being skipped.
    """

    def __init__(self):
        self._site = _eggwrapper.SiteHolder(1)

    def iter(self, aln, struct, filtr, max_missing, consider_outgroup_missing):

        """
        Return the iterator over sites of an alignment. Warning, the
        returned site is static (it is the same object all along and it
        is reset at the end of the iteration).

        :param aln: an :class:`.Align` instance.
        :param struct: a :class:`.Structure` instance (or ``None`` for
            unstructured, haploid data).
        :param filtr: a :class:`.Filter` instance.
        :param max_missing: maximum proporiton of missing data
        :param consider_outgroup_missing: whether include the outgroup to
            evaluate the maximum number of missing data.

        :return: An iterator yielding a :class:`.Site` instance which is
            updated for each site of the alignment (skipping those with
            too many missing data).
        """


        if struct is None:
            pl = 1
        else:
            struct = struct._obj
            pl = struct.get_ploidy()

        if max_missing < 0 or max_missing > 1: raise ValueError, 'max_missing value is out of range'
        max_missing = int(round(max_missing * aln._obj.get_nsam_i()))

        for i in xrange(aln._obj.get_nsit()):
            self._site.reset(pl)
            if not self._site.process_align(aln._obj, i, struct,
                    filtr._obj, max_missing, consider_outgroup_missing):
                continue
            yield self._site

        self._site.reset(1)
