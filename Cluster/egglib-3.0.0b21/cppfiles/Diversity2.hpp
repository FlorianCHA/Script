/*
    Copyright 2008-2009,2013,2015-2016 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_DIVERSITY2_HPP
#define EGGLIB_DIVERSITY2_HPP

namespace egglib {

    class SiteHolder;
    class SiteDiversity;
    class FreqHolder;

   /** \brief Compute population summary statistics from an array of sites
    *
    * \ingroup diversity
    *
    * Diversity instances cannot be copied. This class is designed to
    * allow reuse of objects without unnecessary memory reallocation.
    *
    * This class computes statistics that require access to the full
    * Site instance (and therefore to the individual allele of each
    * individual). Sites with missing data are ignored when computing
    * Wall's statistics.
    *
    * Meaning of flag:
    * * flag&1      an error occurred (+ one of 2, 4, 8, 16, 32)
    * * flag&2      error: less than 2 samples (including missing)
    * * flag&4      error: inconsistent number of samples
    * * flag&8      error: inconsistent ploidy
    * * flag&16     error: inconsistent frequency holder (sample size not checked)
    * * flag&32     error: provided SiteDiversity does not have proper data
    * * flag&64     at least 1 polymorphic site with at least 2 non-missing samples
    * * flag&128    at least 1 polymorphic, orientable site with at least 2 non-missing samples
    * * flag&256    computed R2, R3, R4, and Ch
    * * flag&512    computed R2E, R3E, R4E, and ChE
    * * flag&1024   computed B and Q (at least 2 sites with no missing data)
    *
    * Header: <egglib-cpp/Diversity.hpp>
    *
    */
    class Diversity2 {

        private:
            Diversity2(const Diversity2& src) {}
            Diversity2& operator=(const Diversity2& src) {return *this;}

            void init();
            void free();
            void _compute_singletons();
            void _compute_partitions();
            void _alloc_samples();

            unsigned int _num_sites;
            unsigned int _num_clear;
            unsigned int _num_siteso;
            unsigned int _ploidy;
            unsigned int _num_samples;
            unsigned int _num_sequences;
            unsigned int _res_sequences;
            unsigned int * _singletons;      // size: res_samples
            unsigned int * _extsingletons;   // size: res_samples
            double _k;
            double _ko;
            double _R2;
            double _R3;
            double _R4;
            double _Ch;
            double _R2E;
            double _R3E;
            double _R4E;
            double _ChE;
            unsigned int _Bp;
            double _B;
            double _Q;
            unsigned int _num_dihap;
            unsigned int _res_dihap;
            unsigned int ** _dihap;          // size: res_dihap * 2
            unsigned int * _site_cache;      // size: res_samples
            unsigned int _num_part;
            unsigned int _res_part;
            unsigned int * _res_part2;       // size: _res_part
            unsigned int ** _part;           // size: _res_part * _res_part2[i]
            unsigned int _flag;
            bool _flag_singletons;
            bool _flag_partitions;
            bool _flag_multiple;

        public:

            Diversity2(); ///< Constructor
            ~Diversity2(); ///< Destructor
            void reset(); ///< Restore all variables to the default state (except toggled flags)
            void toggle_off(); ///< Cancel flags
            void toggle_singletons(); ///< Activate computation of Rx/Ch RxE/ChE stats
            void toggle_partitions(); ///< Activate computation of B and Q stats (must be set before load()
            void load(const SiteHolder&,
                const SiteDiversity&, const FreqHolder&); ///< Load site with its matching structure (requires basic stats)
            void _check_partition(); // helper to load()
            void _add_dihap(); // idem
            void set_option_multiple(bool); ///< Toggle option for multiple alleles

            unsigned int num_sites() const; ///< Number of loaded sites (only polymorphic)
            unsigned int num_orientable() const; ///< Number of orientable sites
            unsigned int num_clear() const;  ///< Number of sites with 0 missing data (Wall stats)
            unsigned int num_sequences() const;  ///< Number of sequences

            unsigned int compute(); ///< Compute singletons and/or partitions stats, return flag

            double k() const; ///< Average number of differences
            double ko() const; ///< Average number of differences at orientable sites
            double R2() const; ///< Ramos-Onsins and Rozas's statistic
            double R3() const; ///< Ramos-Onsins and Rozas's statistic
            double R4() const; ///< Ramos-Onsins and Rozas's statistic
            double Ch() const; ///< Ramos-Onsins and Rozas's statistic
            double R2E() const; ///< Ramos-Onsins and Rozas's statistic
            double R3E() const; ///< Ramos-Onsins and Rozas's statistic
            double R4E() const; ///< Ramos-Onsins and Rozas's statistic
            double ChE() const; ///< Ramos-Onsins and Rozas's statistic
            double B() const; ///< Wall's statistic
            double Q() const; ///< Wall's statistic
    };
}

#endif
