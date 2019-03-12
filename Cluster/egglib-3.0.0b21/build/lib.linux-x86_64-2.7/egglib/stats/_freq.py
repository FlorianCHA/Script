__license__ = """
    Copyright 2016 Stephane De Mita, Mathieu Siol

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
from . import _structure

########################################################################

def freq_from_site(site, struct=None):

    """
    Create a new :class:`.Freq` instance based on data of the provided
    site.

    :param site: a :class:`.Site` instance.
    :param struct: this argument can be:

        * A :class:`.Structure` instance
          with ploidy equal to 1 (since the individual level must already
          be implemented in the provided :class:`.Site`), allowing to select a
          subset of samples and/or define the structure. This is the only way
          to specify a hierachical (with clusters) structure.
        * A :class:`list` (or compatible) of at least one integer providing
          the sample size (as numbers of individuals) of all populations,
          assuming the individuals are organized in the corresponding order
          (all individuals of a given population grouped together).
        * ``None`` (no structure, all individuals placed in a single population).

    :return: A new :class:`.Freq` instance.
    """

    obj = Freq()
    obj.process_site(site, struct)
    return obj

########################################################################

def freq_from_list(ingroup, outgroup, geno_list=None):

    """
    Create a new :class:`.Freq` instance based on already computed
    frequency data.

    :param ingroup: a nested list of genotype or allele frequencies (based
        on the value of *geno_list*). The list must have three levels:
        (i) clusters, (ii) populations, and (iii) alleles/genotypes. If
        clusters/populations are not defined, a single-item list should
        be provided for this level. Empty lists are also allowed. The
        frequencies must be null or positive integers. The number of
        frequencies per population is required to be constant for all
        populations (corresponding to the number of alleles or genotypes).
        If *geno_list* is defined, data must be the frequencies of the
        provided genotypes, in the same order. Otherwise, data must be
        allelic frequencies, in the order of increasing allele index (in
        the latter case, data will be loaded as haploid).
    :param outgroup: a list of allele/genotype frequencies for the
        outgroup. The number of alleles or genotypes is also required to
        match. If ``None``, no outgroup samples (equivalent to a list of
        zeros).
    :param geno_list: list of genotypes. Genotypes must be provided as
        tuples or lists. Their length is equal to the ploidy and is required to
        be at least two and constant for all genotypes. Order
        of alleles within genotypes is significant. If ``None``, data are
        loaded as haploid alleles.

    :return: A new :class:`.Freq` instance.

    Note that it is required that there is at least one cluster and one
    population.
    """

    obj = Freq()
    obj.process_list(ingroup, outgroup, geno_list)
    return obj

########################################################################

def freq_from_vcf(vcf):

    """
    Import allelic frequencies from a VCF parser. The VCF parser
    must have processed a variant and the variant is required to have
    frequency data available as the AC format field along with the AN
    field. An exception is raised otherwise.

    This function only imports haploid allele frequencies in the ingroup
    (without structure). The first allele is the reference, by
    construction, then all the alternate alleles in the order in which
    they are provided in the VCF file.

    :param vcf: a :class:`.VcfParser` instance containing data. There
        must at least one sample and the AN/AC format fields must be
        available. It is not required to extract variant data manually.
    """

    obj = Freq()
    obj.process_vcf(vcf)
    return obj

########################################################################

class Freq(object):

    """
    Hold allelic and genotypic frequencies for a single site. :class:`!.Freq`
    instances can be created using the three functions :func:`.freq_from_site`,
    :func:`.freq_from_list`, and :func:`.freq_from_vcf`, or using the
    default constructor. After it is created by any way, instances can be
    re-used (which is faster), using their methods :meth:`~.Freq.process_site`,
    :meth:`~.Freq.process_list`, and :meth:`~.Freq.process_vcf`.
    """

    ingroup = 0
    outgroup = 1
    cluster = 2
    population = 3

    ####################################################################

    def __init__(self):
        self._obj = _eggwrapper.FreqHolder()

    ####################################################################

    def process_site(self, site, struct=None):

        """
        Reset the instance as if it had been created using :func:`.freq_from_site`.
        Arguments are identical to this function.
        """

        if struct is None:
            self._obj.setup_raw(1, 1, site._obj.get_nout(), site._obj.get_ploidy())
            self._obj.setup_pop(0, 0, 0, site._obj.get_ning())
        elif isinstance(struct, _structure.Structure):
            struct = struct._obj
            if struct.get_ploidy() != 1: raise ValueError, 'invalid structure (ploidy is required to be 1)'
            if struct.get_ns_req() > site._obj.get_ning(): raise ValueError, 'invalid structure (ingroup sample index out of range)'
            if struct.get_no_req() > site._obj.get_nout(): raise ValueError, 'invalid structure (outgroup sample index out of range)'
            self._obj.setup_structure(struct, site._obj.get_ploidy())
        else:
            if len(struct) < 1: raise ValueError, 'there must be at least one population size'
            if sum(struct) != site._obj.get_ning(): raise ValueError, 'invalid structure (sample size is required to match)'
            self._obj.setup_raw(1, len(struct), site._obj.get_nout(), site._obj.get_ploidy())
            for i, ns in enumerate(struct): self._obj.setup_pop(i, 0, i, ns)
            flag = True
        self._obj.process_site(site._obj)

    ####################################################################

    def process_list(self, ingroup, outgroup, geno_list=None):

        """
        Reset the instance as if it had been created using :func:`.freq_from_list`.
        Arguments are identical to this function.
        """

        # get structure properties
        nc = len(ingroup)
        ns = [sum(pop) for clu in ingroup for pop in clu]
        np = len(ns)
        ni = sum(ns)

        # get number of alleles/genotypes
        na = set(map(len, [j for i in ingroup for j in i]))
        if outgroup is not None: na.add(len(outgroup))
        if len(na) != 1: raise ValueError, 'number of frequencies must be the same for all populations and outgroup'
        na = na.pop()
        if na < 1: raise ValueError, 'there must be at least one allele'
        if outgroup is None: outgroup = [0] * na
        no = sum(outgroup)
        if geno_list is None:
            pl = 1
            ng = 0
        else:
            ng = na
            if ng != len(geno_list): raise ValueError, 'invalid number of genotypes'
            pl = set(map(len, geno_list))
            if len(pl) != 1: raise ValueError, 'ploidy is not constant in genotypes'
            pl = pl.pop()
            if pl < 2: raise ValueError, 'ploidy must be at least 2'
            geno_set = set()
            for g in geno_list:
                if g in geno_set: raise ValueError, 'genotype {0} is repeated'.format(g)
                geno_set.add(g)
            na = max([j for i in geno_list for j in i]) + 1

        # setup the instance
        self._obj.setup_raw(nc, np, no, pl)
        idx = 0
        for i, clu in enumerate(ingroup):
            for j, pop in enumerate(clu):
                self._obj.setup_pop(idx, i, j, sum(pop))
                idx += 1
        self._obj.set_nall(na, ng)
        if pl > 1:
            for i, g in enumerate(geno_list):
                for j, v in enumerate(g):
                    self._obj.set_genotype_item(i, j, v)

        # load data
        pop_idx = 0
        ing_frq = self._obj.frq_ingroup()
        for clu_idx, clu_data in enumerate(ingroup):
            clu_frq = self._obj.frq_cluster(clu_idx)
            for pop_data in clu_data:
                pop_frq = self._obj.frq_population(pop_idx)
                for i, n in enumerate(pop_data):
                    if pl > 1:
                        ing_frq.incr_genotype(i, n)
                        clu_frq.incr_genotype(i, n)
                        pop_frq.incr_genotype(i, n)
                        for j in geno_list[i]:
                            ing_frq.incr_allele(j, n)
                            clu_frq.incr_allele(j, n)
                            pop_frq.incr_allele(j, n)
                    else:
                        ing_frq.incr_allele(i, n)
                        clu_frq.incr_allele(i, n)
                        pop_frq.incr_allele(i, n)
                pop_idx += 1

        otg_frq = self._obj.frq_outgroup()
        for i, n in enumerate(outgroup):
            if pl > 1:
                otg_frq.incr_genotype(i, n)
                for j in geno_list[i]:
                    otg_frq.incr_allele(j, n)
            else:
                otg_frq.incr_allele(i, n)

        # load heterozygote genotypes
        if pl > 1:
            for i, g in enumerate(geno_list):
                g = set(g)
                if len(g) > 1:
                    for j in g:
                        self._obj.frq_ingroup().tell_het(i, j)
                        self._obj.frq_outgroup().tell_het(i, j)
                        for k in xrange(self._obj.num_clusters()):
                            self._obj.frq_cluster(k).tell_het(i, j)
                        for k in xrange(self._obj.num_populations()):
                            self._obj.frq_population(k).tell_het(i, j)

    ####################################################################

    def process_vcf(self, vcf):

        """
        Reset the instance as if it had been created using :func:`.freq_from_vcf`.
        Argument is identical to this function.
        """

        if vcf._parser.has_data() == False: raise ValueError, 'data must have been read from VCF parser'
        if vcf._parser.has_AC() == False or vcf._parser.has_AN() == False: raise ValueError, 'VCF data must have AC and AN info fields'

        self._obj.process_vcf(vcf._parser)

    ####################################################################

    @property
    def ploidy(self):

        """
        Ploidy.
        """

        return self._obj.ploidy()

    ####################################################################

    @property
    def num_alleles(self):

        """
        Number of alleles in the whole site.
        """

        return self._obj.num_alleles()

    ####################################################################

    @property
    def num_genotypes(self):

        """
        Number of genotypes in the whole site.
        """

        return self._obj.num_genotypes()

    ####################################################################

    @property
    def num_clusters(self):

        """
        Number of clusters.
        """

        return self._obj.num_clusters()

    ####################################################################

    @property
    def num_populations(self):

        """
        Number of populations.
        """

        return self._obj.num_populations()

    ####################################################################

    @property
    def ploidy(self):

        """
        Ploidy
        """

        return self._obj.ploidy()

    ####################################################################

    def genotype(self, idx):

        """
        Get a genotype, as a tuple of allele indexes.
        """

        if idx<0 or idx>=self._obj.num_genotypes():
            raise IndexError, 'invalid genotype index'
        return tuple(self._obj.genotype_item(idx, i) for i in xrange(self._obj.ploidy()))

    ####################################################################

    def _getter(self, cpt, idx):
        if cpt == self.ingroup: return self._obj.frq_ingroup()
        elif cpt == self.outgroup: return self._obj.frq_outgroup()
        elif cpt == self.cluster:
            if idx is None:
                raise ValueError, 'cluster index is required'
            if idx<0 or idx>= self._obj.num_clusters():
                raise IndexError, 'invalid cluster index'
            return self._obj.frq_cluster(idx)
        elif cpt == self.population:
            if idx is None:
                raise ValueError, 'population index is required'
            if idx<0 or idx>= self._obj.num_populations():
                raise IndexError, 'invalid population index'
            return self._obj.frq_population(idx)
        else: raise ValueError, 'invalid compartment identifier'

    ####################################################################

    def freq_allele(self, allele, cpt=ingroup, idx=None):

        """
        Get the frequency of an allele.

        :param allele: allele index.
        :param cpt: compartment identifier.
        :param idx: compartment index (required for clusters and
            populations, ignored otherwise).
        """

        if allele<0 or allele>=self._obj.num_alleles(): raise IndexError, 'invalid allele index'
        return self._getter(cpt, idx).frq_all(allele)

    ####################################################################

    def freq_genotype(self, genotype, cpt=ingroup, idx=None):

        """
        Get the frequency of an genotype.

        :param genotype: genotype index.
        :param cpt: compartment identifier.
        :param idx: compartment index (required for clusters and
            populations, ignored otherwise).
        """

        if genotype<0 or genotype>=self._obj.num_genotypes(): raise IndexError, 'invalid genotype index'
        return self._getter(cpt, idx).frq_gen(genotype)

    ####################################################################

    def nieff(self, cpt=ingroup, idx=None):

        """
        Get the number of individuals within a given compartment. In
        the haploid case, this method is identical to :meth:`.nseff`.

        :param cpt: compartment identifier.
        :param idx: compartment index (required for clusters and
            populations, ignored otherwise).
        """

        return self._getter(cpt, idx).nieff()

    ####################################################################

    def nseff(self, cpt=ingroup, idx=None):

        """
        Get the number of samples within a given compartment.

        :param cpt: compartment identifier.
        :param idx: compartment index (required for clusters and
            populations, ignored otherwise).
        """

        return self._getter(cpt, idx).nseff()
