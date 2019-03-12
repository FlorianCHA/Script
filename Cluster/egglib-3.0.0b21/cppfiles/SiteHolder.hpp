/*
    Copyright 2012-2018 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_SITE_HOLDER_HPP
#define EGGLIB_SITE_HOLDER_HPP

namespace egglib {

    class DataHolder;
    class VcfParser;
    class StructureHolder;
    class Filter;

   /** \brief Holds data for a site for diversity analysis
    *
    * Usage of this class: first set the ploidy. Then, either load an
    * alignment or data from a VCF, or individuals manually. Before
    * loading individuals manually, it is required to pre-set the number
    * such as the indexes will exist. If you don't set all samples
    * manually with load_ing() or load_otg(), you must force set the
    * alleles to the default value. Note: the instance is not reset
    * unless you ask it. Data will add up.
    *
    * \ingroup diversity
    *
    */
    class SiteHolder {

        private:

            SiteHolder(SiteHolder& src) {}
            SiteHolder& operator=(SiteHolder& src) { return * this; }

        protected:

            unsigned int _nall;
            unsigned int _nall_ing;
            unsigned int _nall_c;
            int * _alleles;
            unsigned int _ni;
            unsigned int _no;
            unsigned int _ni_c;
            unsigned int _no_c;
            unsigned int * _ing; // size: _ns_c * ploidy // get_pi depends on this structure
            unsigned int * _otg; // size: _ns_c * ploidy // get_po depends on this structure
            unsigned int _ploidy;
            unsigned int _missing;
            unsigned int _missing_ing;
            unsigned int _missing_otg;
            unsigned int _tot_missing;
            unsigned int _consider_outgroup_missing;

            void _init(unsigned int);
            unsigned int _process_allele(int, bool ing);
            void _load_helper_i(unsigned int sam_idx, unsigned int idv, unsigned int chrom);
            void _load_helper_o(unsigned int sam_idx, unsigned int idv, unsigned int chrom);
            const DataHolder * _data_pointer;
            const Filter * _filter_pointer;
            unsigned int _site_index;
            bool _flag;

        public:

            SiteHolder(); ///< Constructor (ploidy = 1)
            SiteHolder(unsigned int ploidy); ///< Constructor
            virtual ~SiteHolder(); ///< Destructor

            void reset(unsigned int ploidy); ///< Reset all to defaults
            void add_ing(unsigned int num); ///< Add ingroup individuals
            void add_otg(unsigned int num); ///< Add outgroup individuals
            void load_ing(unsigned int idv, unsigned int chrom, int allele); ///< Analyze an ingroup allele (MISSINGDATA for missing data)
            void load_otg(unsigned int idv, unsigned int chrom, int allele); ///< Analyze an outgroup allele (MISSINGDATA for missing data)
            unsigned int get_ning() const; ///< Get number of ingroup individuals
            unsigned int get_nout() const; ///< Get number of outgroup individuals
            unsigned int get_ploidy() const;  ///< Get ploidy
            unsigned int get_nall() const; ///< Number of alleles
            unsigned int get_nall_ing() const; ///< Number of alleles in ingroup only
            int get_allele(unsigned int) const; ///< Get an allele (MISSINGDATA for MISSING)
            void set_nall(unsigned int all, unsigned int ing); ///< Set number of alleles (up to you that all is consistent)
            void set_allele(unsigned int, int) const; ///< Set an allele value
            unsigned int get_i(unsigned int idv, unsigned int chrom) const; ///< Get allele index for ingroup (MISSING for missing data)
            unsigned int get_o(unsigned int idv, unsigned int chrom) const; ///< Get allele index for outgroup (MISSING for missing data)
            unsigned int get_straight_i(unsigned int sam) const; ///< Get allele index for ingroup (MISSING for missing data)
            unsigned int get_straight_o(unsigned int sam) const; ///< Get allele index for outgroup (MISSING for missing data)
            const unsigned int * get_pi(unsigned int idv) const; ///< get_i as pointer
            const unsigned int * get_po(unsigned int idv) const; ///< get_o as pointer
            void set_i(unsigned int idv, unsigned int chrom, unsigned int all); ///< Set allele index for ingroup (MISSING for missing data)
            void set_o(unsigned int idv, unsigned int chrom, unsigned int all); ///< Set allele index for outgroup (MISSING for missing data)
            unsigned int get_missing() const; ///< Number of missing alleles found in the last processed data
            unsigned int get_tot_missing() const; ///< Total number of missing alleles
            unsigned int get_missing_ing() const; ///< Total number of missing alleles in ingroup
            unsigned int get_missing_otg() const; ///< Total number of missing alleles in outgroup

           /** \brief Process an alignment.
            *
            * Does not reset instance! Ploidy must be defined before!
            *
            * \param data an alignment.
            * \param idx index of the site to process.
            * \param struc the structure to use (NULL can be passed to
            *   process all samples as haploid individuals).
            * \param filtr the allele filter.
            * \param max_missing maximum number of missing alleles. If there
            *   are more missing data, stop processing and return false.
            *   Then, the instance should absolutely not be used further.
            * \param consider_outgroup_missing if true, consider missing
            *   data for the max_missing argument, otherwise, only count
            *   missing data of the ingroup.
            *
            * \return true if the number of missing data was not exceeded.
            *
            */
            bool process_align(const DataHolder& data, unsigned int idx,
                    const StructureHolder * struc,
                    const Filter& filtr, unsigned int max_missing,
                    bool consider_outgroup_missing);

           /** \brief Import allelic data and compute frequencies from VCF data
            *
            * Beware: this method does not reset the instance.
            *
            * \param data a VcfParser reference containing data and
            *        having the GT format field filled.
            * \param start index of the first sample to consider.
            * \param stop index of the last sample to consider.
            * \param max_missing maximum number of missing alleles.
            *        If this proportion is processing is stopped and
            *        get_missing() returns max_missing + 1. Only missing
            *        data in this data set, and in the ingroup, are considered.
            *
            * \return A boolean specifying whether processing was
            *         completed.
            *
            */
            bool process_vcf(VcfParser& data,
                             unsigned int start, unsigned int stop,
                             unsigned int max_missing);
    };


   /** \brief SiteHolder subclass to transform regular data to genotypic */
    class SiteGeno : public SiteHolder {
        private:
            unsigned int _sz_geno;
            unsigned int _ploidy; // real ploidy: the ploidy of the object is always 1
            unsigned int _sz_ploidy1;
            unsigned int * _sz_ploidy; // _sz_geno
            bool * _homoz; // _sz_geno
            unsigned int ** _geno; // _sz_geno * _sz_ploidy[i]
            bool * _geno_matched; // _sz_ploidy1

            unsigned int _analyse_genotype(const unsigned int *);

            SiteGeno(const SiteGeno& src) : SiteHolder(1) {}
            SiteGeno& operator=(const SiteGeno& src) { return * this; }
        
        public:
            SiteGeno();
            ~SiteGeno();
            void process(const SiteHolder&); ///< reset and get data from a site
            bool homoz(unsigned int genotype) const; ///< tell if a genotype is homozygote
    };
}

#endif
