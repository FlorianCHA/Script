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

#ifndef EGGLIB_DIVERSITY_FS_HPP
#define EGGLIB_DIVERSITY_FS_HPP

namespace egglib {

   /** \brief Compute Fu's Fs
    *
    * This function computes Fu's Fs statistic using haplotype
    * statistics (that should have been computed using the Haplotypes
    * class) and, as a theta estimator, pi provided by the Diversity1
    * class. The values must have been computed using the same data
    * set.
    *
    * Warning: this function is not available for values of n (number
    * of samples) larger than MAX_N_STIRLING. k must be >= 1 and <= n.
    *
    * \param n number of exploited samples for determining the number
    * of haplotypes.
    *
    * \param K number of haplotypes obtained with the same data.
    *
    * \param pi average number of pairwise differences (as theta
    * estimator, per gene).
    *
    * The behaviour of the function is not defined if K < 0. The
    * function returns UNDEF if the value cannot be computed, which
    * can happen:
    *
    *   \li if n is larger than MAX_N_STIRLING;
    *   \li if the sum of probabilities of k values >= K is too close
    * of 0 or 1 (based on the computer's precision);
    *   \li if pi is 0 (no polymorphism);
    *   \li if K > n (which is an error).
    *
    * Header: <egglib-cpp/Fs.hpp>
    *
    */
    double Fs(unsigned int n, unsigned int K, double pi);
}

#endif
