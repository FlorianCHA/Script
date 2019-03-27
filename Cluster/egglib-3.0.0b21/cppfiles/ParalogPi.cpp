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
#include "egglib.hpp"
#include "ParalogPi.hpp"
#include "SiteHolder.hpp"
#include "Structure.hpp"

namespace egglib {

    ParalogPi::ParalogPi() {
        _np = 0;
        _ns = 0;
        _ls = 0;
        _np_c1 = 0;
        _np_c2 = NULL;
        _ns_c = NULL;
        _piw = NULL;
        _pib = NULL;
        _lsi = NULL;
        _lsij = NULL;
        _samples = NULL;
    }

    ///////

    ParalogPi::~ParalogPi() {
        if (_piw) free(_piw);
        if (_lsi) free(_lsi);
        for (unsigned int i=0; i<_np_c1; i++) {
            if (_pib[i]) free(_pib[i]);
            if (_lsij[i]) free(_lsij[i]);
            if (_samples[i]) free(_samples[i]);
        }
        if (_np_c2) free(_np_c2);
        if (_ns_c) free(_ns_c);
        if (_pib) free(_pib);
        if (_lsij) free(_lsij);
        if (_samples) free(_samples);
    }

    ///////

    void ParalogPi::reset(const StructureHolder& str_prl, const StructureHolder& str_idv) {
        _np = str_prl.num_pop();
        _ns = str_idv.num_pop();
        _ls = 0;

        // allocate the per-paralog tables
        if (_np > _np_c1) {
            _piw = (double *) realloc(_piw, _np * sizeof(double));                        if (!_piw) throw EGGMEM;
            _pib = (double **) realloc(_pib, _np * sizeof(double *));                     if (!_pib) throw EGGMEM;
            _lsi = (unsigned int *) realloc(_lsi, _np * sizeof(double));                  if (!_lsi) throw EGGMEM;
            _lsij = (unsigned int **) realloc(_lsij, _np * sizeof(unsigned int *));       if (!_lsij) throw EGGMEM;
            _samples = (unsigned int **) realloc(_samples, _np * sizeof(unsigned int *)); if (!_samples) throw EGGMEM;
            _np_c2 = (unsigned int *) realloc(_np_c2, _np * sizeof(unsigned int));        if (!_np_c2) throw EGGMEM;
            _ns_c = (unsigned int *) realloc(_ns_c, _np * sizeof(unsigned int));          if (!_ns_c) throw EGGMEM;
            for (unsigned int i=_np_c1; i<_np; i++) {
                _np_c2[i] = 0;
                _ns_c[i] = 0;
                _pib[i] = NULL;
                _lsij[i] = NULL;
                _samples[i] = NULL;
            }
            _np_c1 = _np;
        }

        // allocate the second level of np*np tables
        for (unsigned int i=0; i<_np; i++) {
            _piw[i] = 0.0;
            _lsi[i] = 0;
            if ((_np-i-1) > _np_c2[i]) {
                _pib[i] = (double *) realloc(_pib[i], (_np-i-1) * sizeof(double));               if (!_pib[i]) throw EGGMEM;
                _lsij[i] = (unsigned int *) realloc(_lsij[i], (_np-i-1) * sizeof(unsigned int)); if (!_lsij[i]) throw EGGMEM;
                _np_c2[i] = (_np-i-1);
            }
            for (unsigned int j=0; j<(_np-i-1); j++) {
                _pib[i][j] = 0.0;
                _lsij[i][j] = 0;
            }
        }

        // allocate the second level of samples table
        for (unsigned int i=0; i<_np; i++) {
            if (_ns > _ns_c[i]) {
                _samples[i] = (unsigned int *) realloc(_samples[i], _ns * sizeof(unsigned int));
                if (!_samples[i]) throw EGGMEM;
                _ns_c[i] = _ns;
            }
        }

        // set up the sample index table
        for (unsigned int i=0; i<_np; i++) {
            for (unsigned int j=0; j<_ns; j++) {
                _samples[i][j] = MISSING; // initialize because they are allowed to be missing
            }
        }
        for (unsigned int paralog=0; paralog<_np; paralog++) {
            const StructurePopulation& pop = str_prl.get_population(paralog);
            for (unsigned int i=0; i<pop.num_indiv(); i++) {
                unsigned int idx = pop.get_indiv(i).get_sample(0); // assume ploidy=1
                unsigned int indiv = str_idv.get_pop_index(idx);
                if (indiv != MISSING) _samples[paralog][indiv] = idx;
            }
        }
    }

    ///////

    void ParalogPi::load(const SiteHolder& site) {

        bool any_data = false;

        // compute Piw for each paralog and increment lsi's
        unsigned int i1, i2, a1, a2, d, n;
        for (unsigned int p=0; p<_np; p++) {
            d = 0;
            n = 0;
            for (unsigned int i=0; i<_ns-1; i++) {
                i1 = _samples[p][i];        if (i1 == MISSING) continue;
                a1 = site.get_i(i1, 0);     if (a1 == MISSING) continue; // assume ploidy=1
                for (unsigned int j=i+1; j<_ns; j++) {
                    i2 = _samples[p][j];    if (i2 == MISSING) continue;
                    a2 = site.get_i(i2, 0); if (a2 == MISSING) continue; // assume ploidy=1
                    n++;
                    if (a1 != a2) d++;
                }
            }
            if (n > 0) {
                any_data |= 1;
                _lsi[p]++;
                _piw[p] += 1.0 * d / n;
            }
        }

        // doing the same for each pair of paralogs
        for (unsigned int p=0; p<_np; p++) {
            for (unsigned int q=p+1; q<_np; q++) {
                d = 0;
                n = 0;
                for (unsigned int i=0; i<_ns; i++) {
                    i1 = _samples[p][i];        if (i1 == MISSING) continue;
                    a1 = site.get_i(i1, 0);     if (a1 == MISSING) continue; // assume ploidy=1
                    for (unsigned int j=0; j<_ns; j++) {
                        if (i==j) continue;
                        i2 = _samples[q][j];        if (i2 == MISSING) continue;
                        a2 = site.get_i(i2, 0);     if (a2 == MISSING) continue; // assume ploidy=1
                        n++;
                        if (a1 != a2) d++;
                    }
                }
                if (n > 0) {
                    any_data |= 1;
                    _lsij[p][q-p-1]++;
                    _pib[p][q-p-1] += 1.0 * d / n;  // that would be n+1 to follow eq. (5), but only if no missing data
                }
            }
        }

        if (any_data) _ls++;
    }

    ///////

    unsigned int ParalogPi::num_sites_tot() const {
        return _ls;
    }

    ///////

    unsigned int ParalogPi::num_sites_paralog(unsigned int index) const {
        return _lsi[index];
    }

    ///////

    unsigned int ParalogPi::num_sites_pair(unsigned int index1, unsigned int index2) const {
        if (index1 < index2) return _lsij[index1][index2-index1-1];
        else return _lsij[index2][index1-index2-1];
    }

    ///////

    unsigned int ParalogPi::num_samples() const {
        return _ns;
    }

    ///////

    unsigned int ParalogPi::num_paralogs() const {
        return _np;
    }

    ///////

    double ParalogPi::Piw(unsigned int index) const {
        return _piw[index];
    }

    ///////

    double ParalogPi::Pib(unsigned int index1, unsigned int index2) const {
        if (index1 < index2) return _pib[index1][index2-index1-1];
        else return _pib[index2][index1-index2-1];
    }
}
