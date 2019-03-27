/*
    Copyright 2016 Stéphane De Mita, Mathieu Siol

    This file is part of the EggLib library.

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
*/

#ifndef EGGLIB_FREQ_HOLDER_HPP
#define EGGLIB_FREQ_HOLDER_HPP

namespace egglib {

    class SiteHolder;
    class StructureHolder;
    class VcfParser;

    /// Class holding frequencies in a given compartment
    class FreqSet {

        private:

            FreqSet(const FreqSet& src) {}
            FreqSet& operator=(const FreqSet& src) { return * this; }

            unsigned int _nall;
            unsigned int _nall_c;
            unsigned int _nall_eff;
            unsigned int _gen_eff;
            unsigned int _ngen;
            unsigned int _ngen_c;
            unsigned int _ngen_eff;
            unsigned int _nsam;
            unsigned int _nind;
            unsigned int _nhet;
            unsigned int * _frq_all; // _nall_c
            unsigned int * _frq_het; // _nall_c
            unsigned int * _frq_gen; // _ngen_c
            bool * _gen_het;         // _ngen_c

            void init();
            void free();

        public:

            FreqSet(); ///< Constructor (all empty)
            ~FreqSet(); ///< Destructor
            void setup(); ///< Setup
            void reset(unsigned int); ///< Set number of alleles (set nsam/ngen to 0)
            void add_genotypes(unsigned int num); ///< Add genotypes
            void incr_allele(unsigned int all_idx, unsigned int num); ///< Increment frequency of a given allele
            void incr_genotype(unsigned int gen_idx, unsigned int num); ///< Increment frequency of a given genotype

           /** \brief Tell the class that genotype i is heterozygote for allele a
            * call it several times!
            * don't change frequencies after that!
            */
            void tell_het(unsigned int i, unsigned int a);

            unsigned int num_alleles() const; ///< Number of alleles (equal to user-provided value)
            unsigned int num_alleles_eff() const; ///< Number of alleles with non-null frequency
            unsigned int num_genotypes() const; ///< Number of genotypes (user-provided)
            unsigned int num_genotypes_eff() const; ///< Number of genotypes with non-null frequency
            unsigned int nseff() const; ///< Total frequency (number of samples)
            unsigned int nieff() const; ///< Total frequency (number of individuals) (0 if haploid)
            unsigned int frq_all(unsigned int) const; ///< Get an allele frequency
            unsigned int frq_gen(unsigned int) const; ///< Get an genotype frequency
            unsigned int frq_het(unsigned int) const; ///< Frequency of heterozygotes have >= 1 copies of allele
            unsigned int tot_het() const; ///< Total frequency of heterozygotes
    };

   /** \brief Class holding frequencies for all compartments for a site
    *
    * Possible uses of this class:
    *
    * * Process a site with structure stored in a StructureHolder:
    *
    *     * setup_structure(structure, ploidy, flag) and keep structure available
    *     * process_site()
    *
    * * Process a site without structure or with a manual structure:
    *
    *     * setup_raw(nc, np, no, ploidy, flag)
    *     * setup_pop(i, cluster, relative, ns) for all populations
    *     * process_site() assuming all individuals are consecutive
    *
    * * Enter frequencies manually:
    *
    *     * setup_raw(nc, np, no, ploidy, flag)
    *     * setup_pop(i, cluster, relative, ns) for all populations
    *     * set_nall(na, ng)
    *     * FreqSet.incr_allele()
    *     * [FreqSet.ncr_genotype()]
    *     * [FreqSet.tell_het()]
    *
    * * Process data from a VCF parser:
    *
    *     * process_vcf()
    *
    */
    class FreqHolder {

        private:
            FreqSet _frq_ing;
            FreqSet _frq_otg;
            FreqSet ** _frq_clu;    // _nclu_c
            FreqSet ** _frq_pop;    // _npop_c
            const StructureHolder * _structure;

            unsigned int _npop;
            unsigned int _npop_c;
            unsigned int _nclu;
            unsigned int _nclu_c;

            unsigned int * _clu_idx;      // _npop_c
            unsigned int * _rel_pop_idx;  // _npop_c
            unsigned int * _pop_ns;       // _npop_c
            unsigned int _outgroup_ns;

            unsigned int _nall;
            unsigned int _nall_c;
            int * _alleles;         // _nall_c
            unsigned int _pl;

            unsigned int _ngen;
            unsigned int _ngen_c;
            unsigned int * _gen_c2;

            unsigned int ** _genotypes;
            bool * _matched;
            unsigned int _matched_c;
            bool * _gen_het;

            FreqHolder(const FreqHolder& src) {}
            FreqHolder& operator=(const FreqHolder& src) {return *this;}

            unsigned int _find_genotype(const unsigned int * alleles);
            void _add_genotypes(unsigned int);
            void _set_frq(unsigned int nc, unsigned int np);
            void _set_nall(unsigned int);

        public:

            FreqHolder(); ///< Constructor
            ~FreqHolder(); ///< Destructor

            void setup_structure(const StructureHolder * structure, unsigned int ploidy); ///< Set up based on provided structure (no individual level)
            void setup_raw(unsigned int nc, unsigned int np, unsigned int no, unsigned int ploidy); ///< Setup manual structure
            void setup_pop(unsigned int i, unsigned int clu_idx, unsigned int rel_idx, unsigned int ns); ///< Follows setup_raw() (for all pops)
            void set_nall(unsigned int na, unsigned int ng); ///< Before loading frequencies manually
            void process_site(const SiteHolder& site); ///< Compute frequencies (structure cannot have individual level)
            void process_vcf(const VcfParser& vcf); ///< Get frequencies from VCF

            const FreqSet& frq_ingroup() const; ///< Get frequencies in whole ingroup
            const FreqSet& frq_outgroup() const; ///< Get frequencies in outgroup
            unsigned int num_clusters() const; ///< Get number of clusters
            unsigned int num_populations() const; ///< Get number of populations
            const FreqSet& frq_cluster(unsigned int) const; ///< Get frequencies in a cluster
            const FreqSet& frq_population(unsigned int) const; ///< Get frequencies in a population
            unsigned int num_alleles() const; ///< Number of alleles
            int allele(unsigned int) const; ///< Get an allele value
            unsigned int ploidy() const; ///< Ploidy
            unsigned int num_genotypes() const; ///< Number of genotypes with non-null frequency (0 if haploid)
            bool genotype_het(unsigned int) const; ///< True if genotype is heterozygote
            const unsigned int * genotype(unsigned int) const; ///< Get a genotype (as array of allele indexes) (none if haploid)
            unsigned int genotype_item(unsigned int, unsigned int) const; ///< Get part of a genotype
            void set_genotype_item(unsigned int i, unsigned int j, unsigned int a); ///< Set part of a genotype
            unsigned int cluster_index(unsigned int) const; ///< Get the index of the cluster of a given population
    };

/* on the way to removal but kept here because not committed
    /// Hold a list of genotype indexes
    class GenotypesTable {
        private:
            GenotypesTable(const GenotypesTable& src) {}
            GenotypesTable& operator=(const GenotypesTable& src) { return *this; }
            unsigned int _num_genotypes;
            unsigned int _sz_genotypes;
            unsigned int _num_pop;
            unsigned int _notg;
            unsigned int _sz_pop;
            unsigned int _sz_otg;
            unsigned int * _nind; // _sz_pop
            unsigned int * _sz_ind; // _sz_pop
            unsigned int ** _genotypes_i; // _sz_pop * _sz_ind[i]
            unsigned int * _genotypes_o; // _sz_otg
            bool * _heter; // _sz_genotypes

            void _init();
            void _free();

        public:
            GenotypesTable(); /// create
            ~GenotypesTable(); /// destroy
            void reset(); /// set number of genotypes to 0
            void set_num_pop(unsigned int num_pop); /// set number of pop (no reset)
            void set_num_idv(unsigned int pop, unsigned int num); /// set number of indiv / pop (no default)
            void set_num_otg(unsigned int num); /// set number of outgroup individuals (no default)
            void add_genotype(); /// add a genotype
            unsigned int num_genotypes(); /// number of genotypes
            void set_heterozygote(unsigned int genotype, bool heterozygote); /// tell if genotype is heterozygote (no default)
            void set_genotype_i(unsigned int pop, unsigned int idv, unsigned int genotype); /// set genotype index for individual
            void set_genotype_o(unsigned int idv, unsigned int genotype); /// set genotype index for individual (outgroup
            unsigned int num_pop() const; /// number of populations
            unsigned int num_idv(unsigned int) const; /// number of individuals in a population
            unsigned int num_otg() const; /// number of individuals in outgroup
            unsigned int genotype_i(unsigned int, unsigned int) const; /// genotype index of an ingroup individual
            unsigned int genotype_o(unsigned int) const; /// genotype index of an outgroup individual
            bool genotype_heter(unsigned int) const; /// tell if genotype is heterozygote
            bool ing_heter(unsigned int, unsigned int); /// tell if individual is heterozygote
            bool otg_heter(unsigned int); /// tell if outgroup individual is heterozygote
    };
*/
}
#endif
