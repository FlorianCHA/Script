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

#ifndef EGGLIB_RANDOM_HPP
#define EGGLIB_RANDOM_HPP

namespace egglib {

   /** \brief Pseudo-random number generator
    *
    * \ingroup core
    *
    * This class implements the Mersenne Twister algorithm for
    * pseudo-random number generation. It is based on work by Makoto
    * Matsumoto and Takuji Nishimura (see <http://www.math.sci.hiroshima-u.ac.jp/~m-MAT/MT/emt.html>)
    * and Jasper Bedaux (see <http://www.bedaux.net/mtrand/>) for the
    * core generator, and the Random class of Egglib up to 2.2 for
    * conversion to other laws than uniform.
    *
    * Note that different instances of the class have independent chain
    * of pseudo-random numbers. If several instances have the same seed,
    * they will generate the exact same chain of pseudo-random numbers.
    * Note that this applies if the default constructor is used and that
    * instances are created within the same second.
    *
    * All non-uniform distribution laws generators are based either on
    * the rand_int32() or the standard (half-open, 32 bit) uniform()
    * methods.
    *
    * Header: <egglib-cpp/Random.hpp>
    *
    */
    class Random {

        public:

           /** \brief Constructor with default seed
            *
            * Uses the current system clock second as seed.
            *
            */
            Random();

           /** \brief Constructor with custom seed
            *
            * Favor large, high-complexity seeds. When using different
            * instances of Random in a program, or different processes
            * using Random, ensure they are all seeded using different
            * seeds.
            *
            */
            Random(unsigned long s);

           /** \brief Re-seed an instance
            *
            * Favor large, high-complexity seeds. When using different
            * instances of Random in a program, or different processes
            * using Random, ensure they are all seeded using different
            * seeds.
            *
            */
            void set_seed(unsigned long s);

            /** \brief Get seed value
             *
             * Return the value of the seed that was used to initiate
             * the instance. If the generator was re-seeded, return the
             * seed value passed at that point.
             *
             */
             unsigned long get_seed() const;

           /** \brief Generate a 32-bit random integer
            *
            * Returns an integer in the range [0, 4294967295] (that is
            * in the range [0, 2^32-1].
            *
            */
            unsigned long rand_int32();

           /** \brief Destructor
            *
            */
            virtual ~Random();

           /** \brief Generate a real in the half-open interval [0,1)
            *
            * 0 is included but not 1.
            *
            */
            double uniform();

           /** \brief Generate a real in the closed interval [0,1]
            *
            * Both 0 and 1 are included.
            *
            */
            double uniformcl();

           /** \brief Generate a real in the open interval (0,1)
            *
            * Neither 0 nor 1 is included.
            *
            */
            double uniformop();

           /** \brief Generate a 53-bit real
            *
            * The value has increased precision: even uniform integer
            * pseudo-random numbers can take a finite number of values
            * (2^32 of them, that is). This method increases the
            * complexity of return values, with a cost as increased
            * computing time.
            *
            */
            double uniform53();

           /** \brief Boolean integer
            *
            * Return true with probability 0.5.
            *
            */
            bool brand();

            /** \brief Draws a uniform integer
             *
             * The argument is the number of values that can be
             * generated. Returns an integer in the range [0, ncards-1].
             * Therefore, ncards is not included in the range.
             *
             */
             unsigned int irand(unsigned int ncards);

            /** \brief Draws a number from an exponential distribution
             *
             * Beware, the argument is the distribution's mean (and is
             * also 1/lambda where lambda is the rate parameter).
             *
             */
             double erand(double expectation);

            /** \brief Draws an integer from a Poisson distribution
             *
             * The argument is the Poisson distribution parameters.
             *
             */
             unsigned int prand(double p);

            /** \brief Draws a number from a normal distribution
             *
             * Return a normal variation with expectation 0 and standard
             * deviation 1. The algorithm used is the polar form of the
             * Box-Muller algorithm. A draw is performed every two calls
             * unless the instance is re-seeded.
             *
             */
             double nrand();

            /** \brief Draws a number from a geometric law
             *
             * The argument is the geometric law parameter.
             *
             */
             unsigned int grand(double);

            /** \brief Draws a number from a binomial law
             *
             * \param n number of tests (must be >=0).
             * \param p test probability.
             *
             */
            unsigned long binomrand(long n, double p);

        private:

            // Seed used at creation or re-seeding
            unsigned long _seed;

            // Mersenne Twister parameter
            static const unsigned int n = 624;

            // Mersenne Twister parameter
            static const unsigned int m = 397;

            // Stored bits
            unsigned long* state;

            // Current bit position
            unsigned int pos;

            // Internal helper
            unsigned long twiddle(unsigned long u, unsigned long v);

            // Generate a new state
            void gen_state();

            // Normal value cached boolean
            bool b_ncached;

            // Cached normal value
            double v_ncached;

            // binomrand helper
            unsigned long _binomrand_inversion(long n, double p);

            // binomrand helper
            unsigned long _binomrand_btpe(long n, double p);

            // binomial data cached boolean
            bool _binom_cache;

            // cached binomial data
            double _binom_p;
            long _binom_n;
            double _binom_r;
            double _binom_q;
            double _binom_fm;
            unsigned long _binom_m;
            double _binom_p1;
            double _binom_xm;
            double _binom_xl;
            double _binom_xr;
            double _binom_c;
            double _binom_laml;
            double _binom_lamr;
            double _binom_p2;
            double _binom_p3;
            double _binom_p4;

            /// \brief Copy constructor is disabled and not available
            Random(const Random& src) {}

            /// \brief Copy assignment operator disabled and not available
            Random& operator=(const Random& src) { return *this; }
    };
}

#endif
