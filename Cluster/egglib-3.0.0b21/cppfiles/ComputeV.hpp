/*
    Copyright 2016 Stéphane De Mita

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

#ifndef EGGLIB_COMPUTEV_HPP
#define EGGLIB_COMPUTEV_HPP

namespace egglib {

    class FreqSet;
    class SiteHolder;

   /** \brief Compute allele size variance
    *
    * \ingroup diversity
    *
    */
    class ComputeV {

        private:
            unsigned int _num_all;
            unsigned int _res_all;
            int * _alleles;
            unsigned int _num_sites;
            double _acc_V;

        public:
            ComputeV(); ///< Constructor
            ~ComputeV(); ///< Destructor
            void reset(); ///< Reset
            void setup_alleles_from_site(const SiteHolder&); ///< Get alleles directly from site
            void setup_alleles(unsigned int); ///< Specify number of alleles
            void set_allele(unsigned int, int); ///< Set an allele value
            double compute(const FreqSet&); ///< Compute V (UNDEF if not computable)
            double average() const; ///< Get average V (UNDEF if no computed values)
            unsigned int num_sites() const; ///< Number of sites with computed V
    };
}

#endif
