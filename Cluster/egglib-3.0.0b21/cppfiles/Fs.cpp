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

#include <cmath>
#include "egglib.hpp"
#include "stirling.hpp"
#include "Fs.hpp"

namespace egglib {

    double Fs(unsigned int n, unsigned int K, double pi) {

        // escape if the number of samples is too large to be supported or if no polymorphism

        if (n > MAX_N_STIRLING) return UNDEF;
        if (pi == 0.0) return UNDEF;

        // compute Ewens's recursion

        double logPi = log(pi);
        double Ew = logPi;
        for (unsigned int i=1; i<n; i++) Ew += log(pi + i);

        // compute S'

        double Sp = 0.;

        for (unsigned int k=K; k<=n; k++) Sp += exp(stirling_table(n, k) + k * logPi - Ew);

        // reject extreme values of S'

        if (Sp <= 0.0) return UNDEF;
        if (1.0 - Sp <= 0.0) return UNDEF;

        // compute Fs

        return log(Sp) - log(1.0-Sp);
    }
}
