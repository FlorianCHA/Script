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

#include <cstdlib>
#include "egglib.hpp"
#include "ComputeV.hpp"
#include "FreqHolder.hpp"
#include "SiteHolder.hpp"

namespace egglib {

    ComputeV::ComputeV() {
        _res_all = 0;
        _num_all = 0;
        _alleles = NULL;
        reset();
    }

    void ComputeV::reset() {
        _num_sites = 0;
        _acc_V = 0.0;
    }

    ComputeV::~ComputeV() {
        if (_alleles) free(_alleles);
    }

    void ComputeV::setup_alleles_from_site(const SiteHolder& site) {
        setup_alleles(site.get_nall());
        for (unsigned int i=0; i<_num_all; i++) {
            _alleles[i] = site.get_allele(i);
        }
    }

    void ComputeV::setup_alleles(unsigned int n) {
        _num_all = n;
        if (n > _res_all) {
            _alleles = (int *) realloc(_alleles, n * sizeof(int));
            if (!_alleles) throw EGGMEM;
            _res_all = n;
        }
    }

    void ComputeV::set_allele(unsigned int i, int a) {
        _alleles[i] = a;
    }

    double ComputeV::compute(const FreqSet& frq) {
        if (frq.nseff() < 2) return UNDEF;
        double _sum = 0.0;
        double _sum2 = 0.0;
        double a;
        for (unsigned int i=0; i<_num_all; i++) {
            a = static_cast<double>(_alleles[i]);
            _sum += frq.frq_all(i) * a;
            _sum2 += frq.frq_all(i) * a * a;
        }
        _sum /= frq.nseff();
        double V = _sum2/frq.nseff() - _sum * _sum;
        V *= frq.nseff() / (frq.nseff()-1.0); 
        _num_sites++;
        _acc_V += V;
        return V;
    }

    unsigned int ComputeV::num_sites() const {
        return _num_sites;
    }

    double ComputeV::average() const {
        if (_num_sites == 0) return UNDEF;
        return _acc_V / _num_sites;
    }
}
