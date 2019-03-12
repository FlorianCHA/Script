/*
    Copyright 2013,2015-2016 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_ALLELESTATUS_HPP
#define EGGLIB_ALLELESTATUS_HPP

namespace egglib {

    class FreqHolder;

   /** \brief Classify alleles and site for frequencies with several populations
    *
    * \ingroup diversity
    *
    * Statistics:
    *   Sp  -- population-specific alleles
    *   Spd  -- population-specific derived alleles
    *   ShP -- number of alleles segregating in at least one pair of populations
    *   ShA -- number of alleles in non-null frequencies in at least one pair of populations
    *   FxA -- number of alleles fixed in at least one population
    *   FxD -- number of fixed differences (two different alleles fixed in a pair of populations)
    *
    * The user must ensure that all passed sites are polymorphic. The user
    * should also probably exclude populations with low sample sizes if they
    * are interested in the number of fixed alleles (populations with no
    * samples are automatically skipped).
    *
    * The statistics are computed for each site. Sums for multi-sites
    * are available as Sp_T and Sp_T1 (and similarly for other statistics).
    * T1 is such as each site is counted only once for any statistic.
    *
    * Header: <egglib-cpp/AlleleStatus.hpp>
    *
    */
    class AlleleStatus {

        public:

            AlleleStatus(); ///< Constructor
            ~AlleleStatus(); ///< Destructor
            void reset(); ///< Reset sums (but keep toggle flag)
            void process(const FreqHolder& freqs); ///< Analyze a site
            void total(); ///< Copy all sums to director accessors

            unsigned int Sp() const; ///< Pop-specific alleles
            unsigned int Sp_T1() const; ///< Pop-specific alleles

            unsigned int Spd() const; ///< Pop-specific derived alleles
            unsigned int Spd_T1() const; ///< Pop-specific derived alleles

            unsigned int ShP() const; ///< Shared polymorphisms
            unsigned int ShP_T1() const; ///< Shared polymorphisms

            unsigned int ShA() const; ///< Shared alleles
            unsigned int ShA_T1() const; ///< Shared alleles

            unsigned int FxD() const; ///< Fixed differences
            unsigned int FxD_T1() const; ///< Fixed differences

            unsigned int FxA() const; ///< Fixed alleles
            unsigned int FxA_T1() const; ///< Fixed alleles

            unsigned int nsites() const; ///< Number of sites with valid data
            unsigned int nsites_o() const; ///< Number of orientable sites with valid data

        private:

            AlleleStatus(const AlleleStatus& src) {}
            AlleleStatus& operator=(const AlleleStatus& src) {return *this;}

            unsigned int _Sp, _Sp_T, _Sp_T1,
                         _Spd, _Spd_T, _Spd_T1,
                         _ShP, _ShP_T, _ShP_T1,
                         _ShA, _ShA_T, _ShA_T1,
                         _FxA, _FxA_T, _FxA_T1,
                         _FxD, _FxD_T, _FxD_T1;
            unsigned int _npop, _nall;
            unsigned int _nsites, _nsites_o;

            void _check_Sp(const FreqHolder&);
            void _check_Spd(const FreqHolder&);
            void _check_ShP(const FreqHolder&);
            void _check_ShA(const FreqHolder&);
            void _check_FxD(const FreqHolder&);
            void _check_FxA(const FreqHolder&);
    };
}

#endif
