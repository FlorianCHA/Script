/*
    This class is an implementation of the Mersenne Twister
    pseudo-random number generator that has period 2^(19937)-1. The
    algorithm and original C source code have been designed by M.
    Matsumoto and T. Nishimura: Mersenne Twister: A 623-Dimensionally
    Equidistributed Uniform Pseudo-Random Number Generator, ACM
    Transactions on Modeling and Computer Simulation, Vol. 8, No. 1.
    <http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html>.

    Copyright 1997-2002 Makoto Matsumoto and Takuji Nishimura, all
    rights reserved

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:

    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in
       the documentation and other materials provided with the
       distribution.

    3. The names of its contributors may not be used to endorse or
       promote products derived from this software without specific
       prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
    COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
    INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
    BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
    ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.

    Any feedback is very welcome.
    http://www.math.keio.ac.jp/matumoto/emt.html
    email: matumoto@math.keio.ac.jp

    The generators returning floating point numbers are based on a
    version by Isaku Wada.

    Copyright 2002 Isaku Wada

    The Mersenne Twister random number generator has been ported to C++
    by Jesper Bedaux <http://www.bedaux.net/mtrand/>.

    Copyright 2003 Jasper Bedaux

    Feedback about the C++ port should be sent to Jasper Bedaux,
    see http://www.bedaux.net/mtrand/ for e-mail address and info.

    The interface was revised and non-uniform random number generators
    have been included from the Random class of EggLib v2.1.3

    Copyright 2008, 2009, 2012 Stéphane De Mita, Mathieu Siol

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

    The binomrand function contains code adapted from NumPy version
    v1.8.0.

    Copyright 2005 Robert Kern

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation files
    (the "Software"), to deal in the Software without restriction,
    including without limitation the rights to use, copy, modify, merge,
    publish, distribute, sublicense, and/or sell copies of the Software,
    and to permit persons to whom the Software is furnished to do so,
    subject to the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.
*/

#include "egglib.hpp"
#include "Random.hpp"
#include <cstdlib>
#include <ctime>
#include <cmath>

namespace egglib {

    // Default constructor
    Random::Random() {
        state = (unsigned long*) malloc(n * sizeof(unsigned long));
        if (!state) throw EGGMEM;
        for (unsigned int i=0; i<n; i++) state[i] = 0;
        pos = 0;
        set_seed(time(NULL));
        b_ncached = false;
        v_ncached = 0.;
        _binom_cache = false;
    }


    // Seeded constructor
    Random::Random(unsigned long s) {
        state = (unsigned long*) malloc(n * sizeof(unsigned long));
        if (!state) throw EGGMEM;
        for (unsigned int i=0; i<n; i++) state[i] = 0;
        pos = 0;
        set_seed(s);
        b_ncached = false;
        v_ncached = 0.;
        _binom_cache = false;
    }    


    // Seeding method
    void Random::set_seed(unsigned long s) {
        state[0] = s & 0xFFFFFFFFUL;
        for (unsigned int i = 1; i < n; ++i) {
            state[i] = 1812433253UL * (state[i - 1] ^ (state[i - 1] >> 30)) + i;
            state[i] &= 0xFFFFFFFFUL;
        }
        pos = n;
        b_ncached = false;
        v_ncached = 0.;
        _binom_cache = false;
        _seed = s;
    }

    // Seed accessor
    unsigned long Random::get_seed() const {
        return _seed;
    }

    // Core generator
    unsigned long Random::rand_int32() {
        if (pos == n) gen_state();
        unsigned long x = state[pos++];
        x ^= (x >> 11);
        x ^= (x << 7) & 0x9D2C5680UL;
        x ^= (x << 15) & 0xEFC60000UL;
        return x ^ (x >> 18);
    }


    // Empty destructor
    Random::~Random() {
        if (state) free(state);
    }


    // Boolean
    bool Random::brand() {
        return rand_int32() < 2147483648;
            // true if rand int < 2^32 / 2
    }
    

    // [0, 1) uniform
    double Random::uniform() {
        return static_cast<double>(rand_int32()) * (1. / 4294967296.);
            // rand int  / 2^32
    }


    // [0, 1] uniform
    double Random::uniformcl() {
        return static_cast<double>(rand_int32()) * (1. / 4294967295.);
            // rand int / (2^32 - 1)
    }


    // (0, 1) uniform
    double Random::uniformop() {
        return (static_cast<double>(rand_int32()) + .5) * (1. / 4294967296.);
            // rand int half-shifted right / 2^32
    }


    // [0, 1) uniform, 53 bits
    double Random::uniform53() {
        return (static_cast<double>(rand_int32() >> 5) * 67108864. + 
          static_cast<double>(rand_int32() >> 6)) * (1. / 9007199254740992.);
    }


    // Helper
    unsigned long Random::twiddle(unsigned long u, unsigned long v) {
        return (((u & 0x80000000UL) | (v & 0x7FFFFFFFUL)) >> 1) ^ ((v & 1UL) ? 0x9908B0DFUL : 0x0UL);
    }


    // Generate new state
    void Random::gen_state() {
        for (unsigned int i = 0; i < (n - m); ++i) {
            state[i] = state[i + m] ^ twiddle(state[i], state[i + 1]);
        }
        for (unsigned int i = n - m; i < (n - 1); ++i) {
            state[i] = state[i + m - n] ^ twiddle(state[i], state[i + 1]);
        }
        state[n - 1] = state[m - 1] ^ twiddle(state[n - 1], state[0]);
        pos = 0;
    }


    // Exponential real 
    double Random::erand(double esperance) {
        double tp = 0.0;
        while (tp==0.0) tp = uniform();
        return ( -(esperance)*log(tp));
    }


    // Uniform integer in specified range
    unsigned int Random::irand(unsigned int ncards) {
        return (unsigned int) (uniform()*ncards);
    }


    // Poisson real
    unsigned int Random::prand(double mean) {
        unsigned int i=0;
        double cumul;
        cumul= (-1/mean)*log(uniform());
        while (cumul<1) {
            cumul += (-1/mean)*log(uniform());
            i++;
        }
        return i;
    }


    // Geometric integer
    unsigned int Random::grand(double param) {
        if (param==1.) return 1;
        double X = 1.-uniform();
        return (unsigned int) ceil(log(X)/log(1.-param));
    }


    // Normal real
    double Random::nrand() {

        // return cached value, if so

        if (b_ncached) {
            b_ncached = false;
            return v_ncached;
        }
        
        // polar form of the Box-Muller transformation
        // implementation taken as is from http://www.taygeta.com/random/gaussian.html Nov 10th 2010
        
        float x1, x2, w, y1, y2;
 
        do {
            x1 = 2.0 * uniform() - 1.0;
            x2 = 2.0 * uniform() - 1.0;
            w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = sqrt( (-2.0 * log( w ) ) / w );
         y1 = x1 * w;
         y2 = x2 * w;

        // cache one value and return the other

        b_ncached = true;
        v_ncached = y2;

        return y1;
    }

    // Binomial (from numpy 1.8.0 numpy/random/mtrand/distributions.c)

    unsigned long Random::binomrand(long n, double p) {
        if (n < 0) throw EggArgumentValueError("`binomrand' cannot take a negative value for `n'");
        if (p <= 0.5) {
            if (p * n <= 30.0) return _binomrand_inversion(n, p);
            else return _binomrand_btpe(n, p);
        }
        else {
            double q = 1.0 - p;
            if (q * n <= 30.0) return n - _binomrand_inversion(n, q);
            else return n - _binomrand_btpe(n, q);
        }
    }

    // helper

    unsigned long Random::_binomrand_btpe(long n, double p) {
        double r, q, fm, p1, xm, xl, xr, c, laml, lamr, p2, p3, p4;
        double a, u, v, s, F, rho, t, A, nrq, x1, x2, f1, f2, z, z2, w, w2, x;
        long m, y, k, i;

        if (_binom_cache == false || _binom_n != n || _binom_p != p) {
            _binom_n = n;
            _binom_p = p;
            _binom_cache = true;
            _binom_r = r = p < 0.5 ? p : 1.0 - p;
            _binom_q = q = 1.0 - r;
            _binom_fm = fm = n * r + r;
            _binom_m = m = (long) floor(fm);
            _binom_p1 = p1 = floor(2.195 * sqrt(n*r*q) - 4.6*q) + 0.5;
            _binom_xm = xm = m + 0.5;
            _binom_xl = xl = xm - p1;
            _binom_xr = xr = xm + p1;
            _binom_c = c = 0.134 + 20.5/(15.3 + m);
            a = (fm - xl) / (fm - xl*r);
            _binom_laml = laml = a * (1.0 + a/2.0);
            a = (xr - fm)/(xr * q);
            _binom_lamr = lamr = a * (1.0 + a/2.0);
            _binom_p2 = p2 = p1 * (1.0 + 2.0*c);
            _binom_p3 = p3 = p2 + c/laml;
            _binom_p4 = p4 = p3 + c/lamr;
        }
        else {
            r = _binom_r;
            q = _binom_q;
            fm = _binom_fm;
            m = _binom_m;
            p1 = _binom_p1;
            xm = _binom_xm;
            xl = _binom_xl;
            xr = _binom_xr;
            c = _binom_c;
            laml = _binom_laml;
            lamr = _binom_lamr;
            p2 = _binom_p2;
            p3 = _binom_p3;
            p4 = _binom_p4;
        }

        /* the while loop below replaces goto-based code  ... */

        while (true) {

            // Step10

            nrq = n*r*q;
            u = uniform()*p4;
            v = uniform();
            if (u <= p1) {
                y = (long) floor(xm - p1*v + u);
                break; // ... goto Step60
            }
            // ... goto Step20

            // Step20

            if (u > p2) { // ... goto Step30

                // Step30
                if (u > p3) { // ... goto Step40

                    // Step40
                    y = (long) floor(xr - log(v)/lamr);
                    if (y > n) continue; // ... goto Step40
                    v = v*(u-p3)*lamr;
                    // .. goto Step50
                }

                // (still Step30)
                else {
                    y = (long) floor(xl + log(v)/laml);
                    if (y < 0) continue; // ... goto Step10
                    v = v*(u-p2)*laml;
                    // .. goto Step50
                }
            }

            // (still Step20)
            else {
                x = xl + (u - p1)/c;
                v = v*c + 1.0 - fabs(m - x + 0.5)/p1;
                if (v > 1.0) continue; // .. goto Step10
                y = (long) floor(x);
                // ... goto Step50
            }

            // Step50
            k = fabs(y - m);
            if (!((k > 20) && (k < ((nrq)/2.0 - 1)))) {

                // ... not gotoing Step52
                s = r/q;
                a = s*(n+1);
                F = 1.0;
                if (m < y) {
                    for (i=m; i<=y; i++) F *= (a/i - s);
                }
                else if (m > y) {
                    for (i=y; i<=m; i++) F /= (a/i - s);
                }
                else {
                    if (v > F) continue; // ... goto Step10
                    else break; // ... goto Step60
                }
            }

            // Step52

            rho = (k/(nrq))*((k*(k/3.0 + 0.625) + 0.16666666666666666)/nrq + 0.5);
            t = -k*k/(2*nrq);
            A = log(v);
            if (A < (t - rho)) break; // ... goto Step60
            if (A > (t + rho)) continue; // ... goto Step10

            x1 = y+1;
            f1 = m+1;
            z = n+1-m;
            w = n-y+1;
            x2 = x1*x1;
            f2 = f1*f1;
            z2 = z*z;
            w2 = w*w;
            if (A > (xm*log(f1/x1)
                    + (n-m+0.5)*log(z/w)
                    + (y-m)*log(w*r/(x1*q))
                    + (13680.-(462.-(132.-(99.-140./f2)/f2)/f2)/f2)/f1/166320.
                    + (13680.-(462.-(132.-(99.-140./z2)/z2)/z2)/z2)/z/166320.
                    + (13680.-(462.-(132.-(99.-140./x2)/x2)/x2)/x2)/x1/166320.
                    + (13680.-(462.-(132.-(99.-140./w2)/w2)/w2)/w2)/w/166320.))
            {
                continue; // .. goto Step10
            }
        }

        // Step60
        if (p > 0.5) y = n - y;
        return y;
    }

    // helper

    unsigned long Random::_binomrand_inversion(long n, double p) {
        double q, qn, np, px, U;
        long X, bound;

        if (_binom_cache == false || _binom_n != n || _binom_p != p) {
            _binom_n = n;
            _binom_p = p;
            _binom_cache = true;
            _binom_q = q = 1.0 - p;
            _binom_r = qn = exp(n * log(q));
            _binom_c = np = n*p;
            bound = np + 10.0 * sqrt(np*q + 1);
            if (n < bound) bound = n;
            _binom_m = bound;
        }
        else {
            q = _binom_q;
            qn = _binom_r;
            np = _binom_c;
            bound = _binom_m;
        }
        X = 0;
        px = qn;
        U = uniform();
        while (U > px) {
            X++;
            if (X > bound) {
                X = 0;
                px = qn;
                U = uniform();
            }
            else {
                U -= px;
                px  = ((n-X+1) * p * px)/(X*q);
            }
        }
        return X;
    }
}
