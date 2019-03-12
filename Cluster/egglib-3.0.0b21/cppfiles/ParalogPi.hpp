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

#ifndef EGGLIB_PARALOG_PI_HPP
#define EGGLIB_PARALOG_PI_HPP

namespace egglib {

    class SiteHolder;
    class StructureHolder;

   /** \brief Specifically designed to compute Innan 2003 statistics for a gene family
    * 
    * \ingroup diversity
    *
    * This class computes the within- and between-paralog Pi of Innan
    * (2003). The within-paralog Pi is the same as the standard Pi,
    * except that it is not unbiased. The between-paralog Pi is the same
    * as Dxy, taking the paralogs as populations, except that one pair
    * of genes (paralogs from the same sample) is not considered.
    *
    * Setup provides the two structure objects describing respectively
    * the structure in paralogs and the structure in samples (two
    * different structure objects are required because they are necessarily
    * non-nested). It is required that the structures of interest are
    * loaded as population levels. Cluster levels are ignored. It is required
    * that ploidy is 1 (if you have genotypic data, skip the individual
    * level). Required tests: ploidy of both structure and of all sites
    * is equal to 1, the maximum index of the paralog structure is
    * represented in all sites (other disagrements are treated as missing
    * data).
    *
    */
    class ParalogPi {

        public:

            ParalogPi(); ///< Constructor (default: 0 pop)
            ~ParalogPi(); ///< Destructor

            void reset(const StructureHolder& str_prl, const StructureHolder& str_idv); ///< Reset and setup structure
            void load(const SiteHolder&); ///< Load a site
            unsigned int num_sites_tot() const; ///< Total number of analyzed sites
            unsigned int num_sites_paralog(unsigned int) const; ///< Number of analyzed sites for a paralog
            unsigned int num_sites_pair(unsigned int, unsigned int) const; ///< Number of analyzed sites for a pair of paralogs
            unsigned int num_samples() const; ///< Number of samples (or: size of each pop)
            unsigned int num_paralogs() const; ///< Number of copies (or: number of pops)
            double Piw(unsigned int) const; ///< Within-paralog Pi
            double Pib(unsigned int, unsigned int) const; ///< Between-paralog Pi

        private:

            ParalogPi(const ParalogPi& src) {}
            ParalogPi& operator=(const ParalogPi& src) { return *this; }

            unsigned int _ls;
            unsigned int _ns;
            unsigned int _np;

            unsigned int _np_c1;
            unsigned int * _np_c2;      // size: np_c1
            double * _piw;              // size: np_c1
            unsigned int * _lsi;        // idem
            double ** _pib;             // size: np_c1 * np_c2[i]-i-1
            unsigned int ** _lsij;      // idem

            unsigned int * _ns_c;       // size: np_c1
            unsigned int ** _samples;   // size: np_c1 * _ns_c[i]
        };
}

#endif
