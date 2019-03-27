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

#include <cstdlib>
#include <cmath>
#include "egglib.hpp"
#include "Diversity2.hpp"
#include "SiteHolder.hpp"
#include "SiteDiversity.hpp"
#include "FreqHolder.hpp"

namespace egglib {

    Diversity2::Diversity2() {
        init();
    }

    //////

    Diversity2::~Diversity2() {
        free();
    }

    //////

    void Diversity2::reset() {
        _num_sites = 0;
        _num_clear = 0;
        _num_siteso = 0;
        _k = 0.;
        _ko = 0.;
        _num_samples = 0;
        _num_sequences = 0;
        _R2 = 0.;
        _R3 = 0.;
        _R4 = 0.;
        _Ch = 0.;
        _R2E = 0.;
        _R3E = 0.;
        _R4E = 0.;
        _ChE = 0.;
        _Bp = 0;
        _B = 0.0;
        _Q = 0.0;
        _flag = 0;
        _num_part = 0;
    }

    //////

    void Diversity2::toggle_off() {
        _flag_partitions = false;
        _flag_singletons = false;
    }

    //////

    void Diversity2::toggle_singletons() {
        _flag_singletons = true;
    }

    //////

    void Diversity2::toggle_partitions() {
        _flag_partitions = true;
    }

    //////

    void Diversity2::set_option_multiple(bool b) {
        _flag_multiple = b;
    }

    //////

    void Diversity2::init() {
        _res_sequences = 0;
        _singletons = NULL;
        _extsingletons = NULL;
        reset();
        toggle_off();
        _res_dihap = 0;
        _dihap = NULL;
        _res_part = 0;
        _res_part2 = NULL;
        _part = NULL;
        _site_cache = NULL;
        _flag_multiple = false;
    }

    //////

    void Diversity2::free() {
        if (_singletons) ::free(_singletons);
        if (_extsingletons) ::free(_extsingletons);
        if (_res_dihap) {
            for (unsigned int i=0; i<_res_dihap; i++) {
                if (_dihap[i]) ::free(_dihap[i]);
            }
            ::free (_dihap);
        }
        if (_site_cache) ::free(_site_cache);
        if (_res_part2) ::free(_res_part2);
        if (_part) {
            for (unsigned int i=0; i<_res_part; i++) {
                if (_part[i]) ::free(_part[i]);
            }
            ::free(_part);
        }
    }

    //////

    void Diversity2::_alloc_samples() {
        if (_num_sequences > _res_sequences) {
            _singletons = (unsigned int *) realloc(_singletons, _num_sequences * sizeof(unsigned int));
            if (!_singletons) throw EGGMEM;
            _extsingletons = (unsigned int *) realloc(_extsingletons, _num_sequences * sizeof(unsigned int));
            if (!_extsingletons) throw EGGMEM;
            _site_cache = (unsigned int *) realloc(_site_cache, _num_sequences * sizeof(unsigned int));
            if (!_site_cache) throw EGGMEM;
            _res_sequences = _num_sequences;
        }

        for (unsigned int i=0; i<_num_sequences; i++) {
            _singletons[i] = 0;
            _extsingletons[i] = 0;
        }
    }

    //////

    void Diversity2::load(const SiteHolder& site, const SiteDiversity& div, const FreqHolder& frq) {

        // skip if previous error
        if ((_flag&1) != 0) return;

        // get or check number of samples
        if (_num_samples == 0) {
            _ploidy = site.get_ploidy();
            _num_samples = site.get_ning();
            _num_sequences = _num_samples * _ploidy;
            if (_num_sequences < 2) { _flag = 3; return; }
            _alloc_samples();
        }
        else {
            if (_num_samples != site.get_ning()) { _flag = 5; return; }
            if (_ploidy != site.get_ploidy()) { _flag = 9; return; }
        }

        // check freq holder
        const FreqSet& ing = frq.frq_ingroup();
        const FreqSet& otg = frq.frq_outgroup();
        if (ing.nseff() < 2 || ing.num_alleles_eff() < 2) return;
        if (ing.num_alleles_eff() > 2 && _flag_multiple == false) return;
        if (ing.num_alleles() != site.get_nall()) { _flag = 17; return; }

        // require diversity data
        if ((div.flag()&2) == 0) {
            if ((div.flag()&1) == 0) _flag = 33;
            return;
        }
        if (div.ns() != frq.frq_ingroup().nseff()) {
            _flag = 33;
            return;
        }

        _num_sites++;
        _flag |= 64;
        _k += div.pairdiff();

        if (div.orientable()) {
            _num_siteso++;
            _flag |= 128;
            _ko += div.pairdiff();
        }

        // detect singletons
        if (div.S()) {
            for (unsigned int i=0; i<_num_samples; i++) {
                for (unsigned int j=0; j<_ploidy; j++) {
                    if (site.get_i(i, j) != MISSING && ing.frq_all(site.get_i(i, j)) == 1) {
                        _singletons[i*_ploidy+j] ++;
                        if (div.orientable() && otg.frq_all(site.get_i(i, j)) == 0) _extsingletons[i*_ploidy+j] ++;
                    }
                }
            }
        }

        // check partitions only if no missing data
        if (_flag_partitions && div.ns() == _num_samples * _ploidy) {
            _num_clear++;
            _flag |= 128;

            // check if site is congruent with previous
            if (_num_clear > 1) {
                unsigned int k;
                _num_dihap = 0;
                for (unsigned int i=0; i<_num_samples; i++) {
                    for (unsigned int j=0; j<_ploidy; j++) {
                        bool found = false;
                        for (k=0; k<_num_dihap; k++) {
                            if (_dihap[k][0] == _site_cache[i*_ploidy+j] && _dihap[k][1] == site.get_i(i, j)) {
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            _add_dihap();
                            _dihap[_num_dihap-1][0] = _site_cache[i*_ploidy+j];
                            _dihap[_num_dihap-1][1] = site.get_i(i, j);
                        }
                    }
                }
            }

            // increment B' if pair is congruent
            if (_num_clear > 1 && _num_dihap == 2) {
                _Bp++;
                _check_partition();
            }

            else {
                // copy cache (skipped if site is identical)
                for (unsigned int i=0; i<_num_samples; i++) {
                    for (unsigned int j=0; j<_ploidy; j++) {
                        _site_cache[i*_ploidy+j] = site.get_i(i, j);
                    }
                }
            }
        }
    }

    //////

    void Diversity2::_check_partition() {
        unsigned int n = 0;
        for (; n<_num_part; n++) {
            _num_dihap = 0;
            for (unsigned int i=0; i<_num_sequences; i++) {
                bool found = false;
                for (unsigned int k=0; k<_num_dihap; k++) {
                    if (_dihap[k][0] == _site_cache[i] && _dihap[k][1] == _part[n][i]) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    _add_dihap();
                    _dihap[_num_dihap-1][0] = _site_cache[i];
                    _dihap[_num_dihap-1][1] = _part[n][i];
                }
            }
            if (_num_dihap == 2) break;
        }

        if (n == _num_part) {
            _num_part++;
            if (_num_part > _res_part) {
                _res_part2 = (unsigned int *) realloc(_res_part2, _num_part * sizeof(unsigned int));
                if (!_res_part2) throw EGGMEM;
                _part = (unsigned int **) realloc(_part, _num_part * sizeof(unsigned int *));
                if (!_part) throw EGGMEM;
                _res_part2[n] = 0;
                _part[n] = NULL;
            }
            if (_num_sequences > _res_part2[n]) {
                _part[n] = (unsigned int *) realloc(_part[n], _num_sequences * sizeof(unsigned int));
                if (!_part[n]) throw EGGMEM;
                _res_part2[n] = _num_sequences;
            }
            for (unsigned int i=0; i<_num_sequences; i++) {
                _part[n][i] = _site_cache[i];
            }
        }
    }

    //////

    void Diversity2::_add_dihap() {
        _num_dihap++;
        if (_num_dihap > _res_dihap) {
            _dihap = (unsigned int **) realloc(_dihap, _num_dihap * sizeof(unsigned int *));
            if (!_dihap) throw EGGMEM;
            _dihap[_num_dihap-1] = (unsigned int *) malloc(2 * sizeof(unsigned int));
            if (!_dihap[_num_dihap-1]) throw EGGMEM;
            _res_dihap = _num_dihap;
        }
    }

    //////

    unsigned int Diversity2::num_sequences() const {
        return _num_sequences;
    }

    //////

    unsigned int Diversity2::compute() {
        if (_flag_singletons) _compute_singletons();
        if (_flag_partitions) _compute_partitions();
        return _flag;
    }

    //////
   
    void Diversity2::_compute_singletons() {
        if (_num_sites > 0) {
            _flag |= 256;
            _R2 = 0.0;
            _R3 = 0.0;
            _R4 = 0.0;
            unsigned int U = 0;
            double kp = _k * _num_sequences / (_num_sequences-1);
            for (unsigned int i=0; i<_num_sequences; i++) {
                _R2 += pow(_singletons[i] - kp / 2.0, 2);
                _R3 += pow(_singletons[i] - kp / 2.0, 3);
                _R4 += pow(_singletons[i] - kp / 2.0, 4);
                U += _singletons[i];
            }

            _R2 /= _num_sequences;
            _R3 /= _num_sequences;
            _R4 /= _num_sequences;
            _R2 = pow(_R2, 1/2.0);
            if (_R3 < 0.0) _R3 = - pow(- _R3, 1/3.0);
            else _R3 = pow(_R3, 1/3.0);
            _R4 = pow(_R4, 1/4.0);
                // todo check ERANGE after calls to pow()
            _R2 /= _num_sites;
            _R3 /= _num_sites;
            _R4 /= _num_sites;
            double m = _k * _num_sequences / (_num_sequences - 1);
            _Ch = pow(U - m, 2) * _num_sites / (m * (_num_sites - m));

            if (_num_siteso > 0) {
                _flag |= 512;
                _R2E = 0.0;
                _R3E = 0.0;
                _R4E = 0.0;
                unsigned int V = 0;
                for (unsigned int i=0; i<_num_sequences; i++) {
                    _R2E += pow(_extsingletons[i] - _ko / 2.0, 2);
                    _R3E += pow(_extsingletons[i] - _ko / 2.0, 3);
                    _R4E += pow(_extsingletons[i] - _ko / 2.0, 4);
                    V += _extsingletons[i];
                }
                _R2E /= _num_sequences;
                _R3E /= _num_sequences;
                _R4E /= _num_sequences;

                _R2E = pow(_R2E, 1/2.0);
                if (_R3E < 0.0) _R3E = - pow(- _R3E, 1/3.0);
                else _R3E = pow(_R3E, 1/3.0);
                _R4E = pow(_R4E, 1/4.0);

                _R2E /= _num_siteso;
                _R3E /= _num_siteso;
                _R4E /= _num_siteso;

                double mo = _ko * _num_sequences / (_num_sequences - 1);
                _ChE = pow(V - mo, 2) * _num_siteso / (mo * (_num_siteso - mo));
            }
        }
    }

    //////

    void Diversity2::_compute_partitions() {
        if (_num_clear > 2) {
            _flag |= 1024;
            _B = (double) _Bp / (_num_clear - 1);
            _Q = static_cast<double>(_Bp + _num_part) / _num_clear;
        }
    }

    //////

    unsigned int Diversity2::num_sites() const { return _num_sites; }
    unsigned int Diversity2::num_clear() const { return _num_clear; }
    unsigned int Diversity2::num_orientable() const { return _num_siteso; }
    double Diversity2::k() const { return _k; }
    double Diversity2::ko() const { return _ko; }
    double Diversity2::R2() const { return _R2; }
    double Diversity2::R3() const { return _R3; }
    double Diversity2::R4() const { return _R4; }
    double Diversity2::Ch() const { return _Ch; }
    double Diversity2::R2E() const { return _R2E; }
    double Diversity2::R3E() const { return _R3E; }
    double Diversity2::R4E() const { return _R4E; }
    double Diversity2::ChE() const { return _ChE; }
    double Diversity2::B() const { return _B; }
    double Diversity2::Q() const { return _Q; }
}
