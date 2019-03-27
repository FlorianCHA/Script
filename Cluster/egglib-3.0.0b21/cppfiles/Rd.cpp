/*
    Copyright 2008-2018 Stéphane De Mita, Mathieu Siol

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

#include "egglib.hpp"
#include <cstdlib>
#include <cmath>
#include "SiteHolder.hpp"
#include "Rd.hpp"

namespace egglib {

    Rd::Rd() {
        _res_loci = 0;
        _res_ploidy = 0;
        _res_pairs = 0;
        _var = NULL;
        _diff = NULL;
        _diff_n = NULL;
        _flags = NULL;
        reset();
    }

    void Rd::reset() {
        _ploidy = 0;
        _Ve = 0.0;
        _num_loci = 0;
        _num_indiv = 0;
        _npt = 0;
        _no_data = true;
        _invalid = false;
    }

    unsigned int Rd::num_loci() const {
        return _num_loci;
    }

    unsigned int Rd::ploidy() const {
        return _ploidy;
    }

    unsigned int Rd::num_indiv() const {
        return _num_indiv;
    }

    void Rd::_configure(unsigned int num_indiv, unsigned int ploidy) {
        reset();
        _num_indiv = num_indiv;
        _ploidy = ploidy;
        _npt = _num_indiv * (_num_indiv - 1) / 2;
        if (_npt > _res_pairs) {
            _diff = (unsigned int *) realloc(_diff, _npt * sizeof(unsigned int));
            if (!_diff) throw EGGMEM;
            _diff_n = (unsigned int *) realloc(_diff_n, _npt * sizeof(unsigned int));
            if (!_diff_n) throw EGGMEM;
            _res_pairs = _npt;
        }
        if (_ploidy > _res_ploidy) {
            _flags = (bool *) realloc(_flags, _ploidy * sizeof(bool));
            if (!_flags) throw EGGMEM;
            _res_ploidy = _ploidy;
        }
        for (unsigned int i=0; i<_npt; i++) {
            _diff[i] = 0;
            _diff_n[i] = 0;
        }
        _no_data = false;
    }

    Rd::~Rd() {
        if (_var) free(_var);
        if (_diff) free(_diff);
        if (_diff_n) free(_diff_n);
        if (_flags) free(_flags);
    }

    void Rd::load(const SiteHolder& site) {

        if (_no_data) _configure(site.get_ning(), site.get_ploidy());
        else if (site.get_ning() != _num_indiv || site.get_ploidy() != _ploidy) {
            _num_loci = 0;
            _invalid = true;
        }
        if (_invalid) return;

        // initialise
        unsigned int sum_d = 0;  // sum of pairwise diff
        unsigned int sum_d2 = 0; // sum of squared pairwise diff
        unsigned int c = 0;      // current pair index
        unsigned int npi = 0;    // number of pairs with data

        // compute number of differences for all pairs
        for (unsigned int i=0; i<_num_indiv-1; i++) {
            for (unsigned int j=i+1; j<_num_indiv; j++) {
                unsigned int d = _cmp_diff(site, i, j);
                if (d != UNKNOWN) {
                    _diff[c] += d;
                    _diff_n[c]++;
                    sum_d += d;
                    sum_d2 += d * d;
                    npi++;
                }
                c++;
            }
        }

        // only process sites if >0 pair of samples
        if (npi > 0) {

            // allocate new site
            _num_loci++;
            if (_num_loci > _res_loci) {
                _var = (double *) realloc(_var, _num_loci * sizeof(double));
                if (!_var) throw EGGMEM;
                _res_loci = _num_loci;
            }

            _var[_num_loci-1] = (sum_d2 - sum_d * sum_d / static_cast<double>(npi)) / npi;
            _Ve += _var[_num_loci-1];
        }
    }

    unsigned int Rd::_cmp_diff(const SiteHolder& site, unsigned int a, unsigned int b) {
        for (unsigned int i=0; i<_ploidy; i++) _flags[i] = true;
        unsigned int d = 0;
        for (unsigned int i=0; i<_ploidy; i++) {
            unsigned int j;
            for (j=0; j<_ploidy; j++) {
                if (_flags[j]) {
                    if (site.get_i(a, i) == MISSING) return UNKNOWN;
                    if (site.get_i(b, j) == MISSING) return UNKNOWN;
                    if (site.get_i(a, i) == site.get_i(b, j)) {
                        _flags[j] = false;
                        break;
                    }
                }
            }
            if (j == _ploidy) d++;
        }
        return d;
    }

    double Rd::compute() {
        double acc = 0.0;
        for (unsigned int i=0; i<_num_loci; i++) {
            for (unsigned int j=i+1; j<_num_loci; j++) {
                acc += sqrt(_var[i] * _var[j]);
            }
        }

        unsigned int npe = 0;
        double sum_D = 0.0;
        double sum_D2 = 0.0;
        for (unsigned int i=0; i<_npt; i++) {
            if (_diff_n[i] > 0) {
                double d = static_cast<double>(_diff[i]) / _diff_n[i] * _num_loci;
                sum_D += d;
                sum_D2 += d*d;
                npe++;
            }
        }

        if (npe > 0) {
            double Vo = (sum_D2 - sum_D * sum_D / npe) / npe;
            if (acc > 0.0) {
                return (Vo - _Ve) / ( 2 * acc);
            }
        }
        return UNDEF;
    }
}
