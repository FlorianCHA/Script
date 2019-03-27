/*
    Copyright 2012-2018 Stéphane De Mita, Mathieu Siol

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
#include "SiteHolder.hpp"
#include "Structure.hpp"
#include "Haplotypes.hpp"

namespace egglib {

    Haplotypes::Haplotypes() {
        _init();
    }

    Haplotypes::~Haplotypes() {
        _free();
    }

    void Haplotypes::_free() {
        if (_c_hapl2) free(_c_hapl2);
        if (_c_hapl4) free(_c_hapl4);
        if (_c_sites) free(_c_sites);
        if (_freq_i) free(_freq_i);
        if (_freq_o) free(_freq_o);
        if (_map) free(_map);
        if (_map_cache) free(_map_cache);
        if (_n_missing) free(_n_missing);
        if (_mis_idx) free(_mis_idx);
        if (_Ki) free(_Ki);
        if (_ni) free(_ni);
        if (_Kd) free(_Kd);
        if (_pop_i) free(_pop_i);
        if (_hapl) {
            for (unsigned int i=0; i<_c_hapl; i++) if (_hapl[i]) free(_hapl[i]);
            free(_hapl);
        }
        if (_potential) {
            for (unsigned int i=0; i<_n_mis; i++) if (_potential[i]) free(_potential[i]);
            free(_potential);
        }
        if (_dist) {
            for (unsigned int i=0; i<_c_hapl3; i++) if (_dist[i]) free(_dist[i]);
            free(_dist);
        }
    }

    void Haplotypes::_init() {
        _c_sam = 0;
        _c_hapl = 0;
        _c_hapl2 = NULL;
        _c_hapl3 = 0;
        _c_hapl4 = NULL;
        _c_sites = NULL;
        _map = NULL;
        _map_cache = NULL;
        _hapl = NULL;
        _freq_i = NULL;
        _freq_o = NULL;
        _c_mis = 0;
        _mis_idx = NULL;
        _potential = NULL;
        _n_missing = NULL;
        _dist = NULL;
        _c_pop = 0;
        _Ki = NULL;
        _ni = NULL;
        _Kd = NULL;
        _pop_i = NULL;
        _pop_o = NULL;
        setup(NULL);
    }

    void Haplotypes::setup(const StructureHolder * struc) {
        _invalid = false;
        _pl = 0;
        _ni_tot = 0;
        _no_tot = 0;
        _nsi = 0;
        _nso = 0;
        _n_mis = 0;
        _nt_hapl = 0;
        _ni_hapl = 0;
        _ng_hapl = 0;
        _n_sites = 0;
        _Fst = UNDEF;
        _Kst = UNDEF;
        _n_pop = 0;
        _ne_pop = 0;
        _nstot = 0;
        _struct = struc;
    }

    void Haplotypes::_set_sample_size(unsigned int pl, unsigned int ni, unsigned int no) {
        _ni_tot = ni;
        _no_tot = no;
        _pl = pl;
        _nsi = ni*pl;
        _nso = no*pl;
        _n_sam = _nsi + _nso;
        if (_n_sam > _c_sam) {
            _map = (unsigned int *) realloc(_map, _n_sam * sizeof(unsigned int));
            if (!_map) throw EGGMEM;
            _map_cache = (unsigned int *) realloc(_map_cache, _n_sam * sizeof(unsigned int));
            if (!_map_cache) throw EGGMEM;
            _n_missing = (unsigned int *) realloc(_n_missing, _n_sam * sizeof(unsigned int));
            if (!_n_missing) throw EGGMEM;
            _pop_i = (unsigned int *) realloc(_pop_i, _n_sam * sizeof(unsigned int));
            if (!_n_sam) throw EGGMEM;
            _c_sam = _n_sam;
        }
        for (unsigned int i=0; i<_n_sam; i++) {
            _map[i] = 0;
            _n_missing[i] = 0;
        }
        _pop_o = _pop_i + _nsi;

        if (_struct == NULL) {
            for (unsigned int i=0; i<_nsi; i++) _pop_i[i] = 0;
            for (unsigned int i=0; i<_nso; i++) _pop_o[i] = 0;
        }

        else {
            // initialize all individuals to "ignored"
            for (unsigned int i=0; i<_nsi; i++) _pop_i[i] = UNKNOWN;
            for (unsigned int i=0; i<_nso; i++) _pop_o[i] = UNKNOWN;

            _n_pop = _struct->num_pop();
            if (_n_pop > _c_pop) {
                _Ki = (unsigned int *) realloc(_Ki, _n_pop * sizeof(unsigned int));
                if (!_Ki) throw EGGMEM;
                _ni = (unsigned int *) realloc(_ni, (_n_pop+1) * sizeof(unsigned int));
                if (!_ni) throw EGGMEM;
                _Kd = (unsigned int *) realloc(_Kd, _n_pop * _n_pop * sizeof(unsigned int));
                if (!_Kd) throw EGGMEM;
                _c_pop = _n_pop;
            }

            for (unsigned int i=0; i<_n_pop; i++) {
                const StructurePopulation& pop = _struct->get_population(i);
                for (unsigned int j=0; j<pop.num_indiv(); j++) {
                    unsigned int sample = pop.get_indiv(j).get_sample(0);
                    for (unsigned int k=0; k<_pl; k++) {
                        _pop_i[sample*_pl+k] = i;
                    }
                }
            }

            for (unsigned int i=0; i<_struct->num_indiv_outgroup(); i++) {
                unsigned int sample = _struct->get_indiv_outgroup(i).get_sample(0);
                for (unsigned int j=0; j<_pl; j++) {
                    _pop_o[sample*_pl+j] = 9999; // any non-UNKNOWN value
                }
            }
        }
    }

    void Haplotypes::load(const SiteHolder& site) {

        if (_invalid) return;

        if (_n_sites == 0) {
            _set_sample_size(site.get_ploidy(), site.get_ning(), site.get_nout());
        }
        else {
            if (_ni_tot != site.get_ning() ||
                _no_tot != site.get_nout() ||
                _pl != site.get_ploidy()) {
                    setup(NULL);
                    _invalid = true;
                    return;
                }
        }

        _n_sites++;

        // record the haplotype of all samples at previous site
        for (unsigned int i=0; i<_n_sam; i++) {
            _map_cache[i] = _map[i];
        }

        // add a column with N's to all haplotypes
        for (unsigned int i=0; i<_nt_hapl; i++) {
            if (_n_sites > _c_sites[i]) {
                _hapl[i] = (unsigned int *) realloc(_hapl[i], _n_sites * sizeof(unsigned int));
                if (!_hapl[i]) throw EGGMEM;
            }
            _hapl[i][_n_sites-1] = MISSING;
        }

        // process all samples
        for (unsigned int i=0; i<_nsi; i++) {
            if (_pop_i[i] != UNKNOWN) {
                _process(site.get_i(i/_pl, i%_pl), i);
            }
        }
        for (unsigned int i=0; i<_nso; i++) {
            if (_pop_o[i] != UNKNOWN) {
                _process(site.get_o(i/_pl, i%_pl), _nsi+i);
            }
        }
    }

    void Haplotypes::_add_hapl() {
        _nt_hapl++;
        if (_nt_hapl > _c_hapl) {
            _c_sites = (unsigned int *) realloc(_c_sites, _nt_hapl * sizeof(unsigned int));
            if (!_c_sites) throw EGGMEM;
            _hapl = (unsigned int **) realloc(_hapl, _nt_hapl * sizeof(unsigned int *));
            if (!_hapl) throw EGGMEM;
            _freq_i = (unsigned int *) realloc(_freq_i, _nt_hapl * sizeof(unsigned int));
            if (!_freq_i) throw EGGMEM;
            _freq_o = (unsigned int *) realloc(_freq_o, _nt_hapl * sizeof(unsigned int));
            if (!_freq_o) throw EGGMEM;
            _c_sites[_nt_hapl-1] = 0;
            _hapl[_nt_hapl-1] = NULL;
            _c_hapl = _nt_hapl;
        }
        if (_n_sites > _c_sites[_nt_hapl-1]) {
            _hapl[_nt_hapl-1] = (unsigned int *) realloc(_hapl[_nt_hapl-1], _n_sites * sizeof(unsigned int));
            _c_sites[_nt_hapl-1] = _n_sites;
        }
    }

    void Haplotypes::_process(unsigned int a, unsigned int idx) {
        // if missing data, they need to be counted
        if (a == MISSING) {
            if (_map[idx] != MISSING) _map[idx] = MISSING;
            _n_missing[idx]++;
            return;
        }

        // only consider sample if it is not marked as missing
        if (_map[idx] == MISSING) return;

        // for first haplotype (necessarily first site)
        if (_nt_hapl == 0) {
            _add_hapl();
            _hapl[0][0] = a;
            return;
            // note: all samples are pre-set to be haplotype 0
        }

        // if first sample of its haplotype, set its base as base of the corresponding haplotype
        if (_hapl[_map[idx]][_n_sites-1] == MISSING) {
            _hapl[_map[idx]][_n_sites-1] = a;
            return;
        }

        // find if the sample = one of the samples that used to be in the same haplotype (at previous site)
        for (unsigned int i=0; i<idx; i++) {
            if (_map[i] != MISSING && _map_cache[i] == _map_cache[idx] && a == _hapl[_map[i]][_n_sites-1]) {
                if (_map[idx] != _map[i]) _map[idx] = _map[i];
                return;
            }
        }

        // if this sample differ from all previous samples of the same haplotype, break haplotype
        _add_hapl();
        for (unsigned int i=0; i<_n_sites-1; i++) _hapl[_nt_hapl-1][i] = _hapl[_map[idx]][i];
        _hapl[_nt_hapl-1][_n_sites-1] = a;
        _map[idx] = _nt_hapl-1;
    }

    void Haplotypes::cp_haplotypes() {
        // initialize the frequencies
        _ne_ing = 0;
        _ne_otg = 0;
        for (unsigned int i=0; i<_nt_hapl; i++) {
            _freq_i[i] = 0;
            _freq_o[i] = 0;
        }
        _ne_ing = 0;
        _ne_otg = 0;

        // compute frequencies
        for (unsigned int i=0; i<_nsi; i++) {
            if (_pop_i[i] == UNKNOWN) _map[i] = MISSING;
            else {
                if (_map[i] != MISSING) {
                    _ne_ing++;
                    _freq_i[_map[i]]++;
                    _nstot++;
                }
            }
        }
        for (unsigned int i=0; i<_nso; i++) {
            if (_pop_o[i] == UNKNOWN) _map[_nsi+i] = MISSING;
            else {
                if (_map[_nsi+i] != MISSING) {
                    _ne_otg++;
                    _freq_o[_map[i+_nsi]]++;
                }
            }
        }

        // compute effective number of haplotypes
        _ng_hapl = 0;
        _ni_hapl = 0;
        for (unsigned int i=0; i<_nt_hapl; i++) {
            if (_freq_i[i] > 0) {
                _ni_hapl++;
                _ng_hapl++;
            }
            else if (_freq_o[i] > 0) _ng_hapl++;
        }
    }

    void Haplotypes::prepare_impute(unsigned int max_missing) {
        _n_mis = 0;
        for (unsigned int i=0; i<_n_sam; i++) {
            if (_pop_i[i] != UNKNOWN) {
                if (_n_missing[i] > 0 && _n_missing[i] <= max_missing) {
                    _n_mis++;
                    if (_n_mis > _c_mis) {
                        _c_hapl2 = (unsigned int *) realloc(_c_hapl2, _n_mis * sizeof(unsigned int));
                        if (!_c_hapl2) throw EGGMEM;
                        _potential = (unsigned int **) realloc(_potential, _n_mis * sizeof(unsigned int *));
                        if (!_potential) throw EGGMEM;
                        _mis_idx = (unsigned int *) realloc(_mis_idx, _n_mis * sizeof(unsigned int));
                        if (!_mis_idx) throw EGGMEM;
                        _c_hapl2[_n_mis-1] = 0;
                        _potential[_n_mis-1] = NULL;
                        _c_mis = _n_mis;
                    }
                    if (_nt_hapl > _c_hapl2[_n_mis-1]) {
                        _potential[_n_mis-1] = (unsigned int *) realloc(_potential[_n_mis-1], (_nt_hapl+1) * sizeof(unsigned int));
                        if (!_potential[_n_mis-1]) throw EGGMEM;
                        _c_hapl2[_n_mis-1] = _nt_hapl;
                    }
                    for (unsigned int k=0; k<_nt_hapl; k++) {
                        if (_freq_i[k] + _freq_o[k] > 0) _potential[_n_mis-1][k] = 1;
                        else _potential[_n_mis-1][k] = 0;
                    }
                    _mis_idx[_n_mis-1] = i;
                    _potential[_n_mis-1][_nt_hapl] = _ng_hapl;
                }
            }
        }
        _site_index = 0;
    }

    void Haplotypes::resolve(const SiteHolder& site) {
        unsigned int indiv, chrom, allele;
        unsigned int (SiteHolder::*getter)(unsigned int, unsigned int) const;
        for (unsigned int i=0; i<_n_mis; i++) {
            if (_mis_idx[i] >= _nsi) {
                indiv = (_mis_idx[i]-_nsi) / _pl;
                chrom = (_mis_idx[i]-_nsi) % _pl;
                getter = &SiteHolder::get_o;
            }
            else {
                indiv = _mis_idx[i] / _pl;
                chrom = _mis_idx[i] % _pl;
                getter = &SiteHolder::get_i;
            }
            allele = (site.*getter)(indiv, chrom);
            if (allele != MISSING) {
                for (unsigned int h=0; h<_nt_hapl; h++) {
                    if (_hapl[h][_site_index] != allele) { // this will catch the case when the haplotype base is missing ("ghost haplotypes")
                        if (_potential[i][h] == 1) {
                            _potential[i][_nt_hapl]--;
                            _potential[i][h] = 0;
                        }
                    }
                }
            }
        }
        _site_index++;
    }

    void Haplotypes::impute() {
        unsigned int c;
        for (unsigned int i=0; i<_n_mis; i++) {
            c = 0;
            for (unsigned int j=0; j<_nt_hapl && c<_potential[i][_nt_hapl]; j++) {
                if (_potential[i][j] == 1) _potential[i][c++] = j;
            }
            if (_potential[i][_nt_hapl] == 1) {
                if (_mis_idx[i] < _nsi) {
                    _ne_ing++;
                    _freq_i[_potential[i][0]]++;
                    _map[_mis_idx[i]] = _potential[i][0];
                }
                else {
                    _ne_otg++;
                    _freq_o[_potential[i][0]]++;
                    _map[_mis_idx[i]] = _potential[i][0];
                }
            }
        }
    }

    void Haplotypes::get_site(SiteHolder& site) {
        site.reset(_pl);
        site.add_ing(_nsi/_pl);
        site.add_otg(_nso/_pl);
        site.set_nall(_ng_hapl, _ni_hapl);
        unsigned int i=0;
        for (unsigned int a=0; a<_nt_hapl; a++) {
            if (_freq_i[a] + _freq_o[a] > 0) site.set_allele(i++, a);
        }
        if (i!=_ng_hapl) throw EggRuntimeError("unexpected case in Haplotypes::get_site()");

        for (unsigned int i=0; i<_nsi; i++) site.set_i(i/_pl, i%_pl, _map[i]);
        for (unsigned int i=0; i<_nso; i++) site.set_o(i/_pl, i%_pl, _map[_nsi+i]);
    }

    void Haplotypes::cp_dist() {
        // alloc
        if (_nt_hapl > _c_hapl3) {
            _c_hapl4 = (unsigned int *) realloc(_c_hapl4, _nt_hapl * sizeof(unsigned int));
            if (!_c_hapl4) throw EGGMEM;
            _dist = (unsigned int **) realloc(_dist, _nt_hapl * sizeof(unsigned int *));
            if (!_dist) throw EGGMEM;
            for (unsigned int i=_c_hapl3; i<_nt_hapl; i++) {
                _dist[i] = NULL;
                _c_hapl4[i] = 0;
            }
            _c_hapl3 = _nt_hapl;
        }
        for (unsigned int i=0; i<_nt_hapl; i++) {
            if (i > _c_hapl4[i]) {
                _dist[i] = (unsigned int *) realloc(_dist[i], i * sizeof(unsigned int));
                if (!_dist[i]) throw EGGMEM;
                _c_hapl4[i] = _nt_hapl;
            }
        }

        // compute (truncated haplotypes are set to 0)
        for (unsigned int i=0; i<_nt_hapl; i++) {
            for (unsigned int j=0; j<i; j++) {
                _dist[i][j] = 0;
                if ((_freq_i[i]>0 || _freq_o[i]>0) && (_freq_i[j]>0 || _freq_o[j]>0)) {
                    for (unsigned int k=0; k<_n_sites; k++) {
                        if (_hapl[i][k] != _hapl[j][k]) _dist[i][j]++;
                    }
                }
            }
        }
    }

    unsigned int Haplotypes::cp_stats() {

        if (_n_pop == 0) return 0;

        // local variables
        double Kt = 0.0;
        for (unsigned int i=0; i<_n_pop; i++) {
            _Ki[i] = 0;
            _ni[i] = 0;
        }
        for (unsigned int i=0; i<_n_pop*_n_pop; i++) {
            _Kd[i] = 0;
        }
        _ni[_n_pop] = 0;
        unsigned int d;

        // process all pairs of samples
        for (unsigned int i=0; i<_nsi; i++) {
            unsigned int a = _map[i];
            if (a == MISSING) continue;
            _ni[_pop_i[i]]++;
            _ni[_n_pop]++;
            for (unsigned int j=i+1; j<_nsi; j++) {
                unsigned int b = _map[j];
                if (b == MISSING) continue;
                if (a == b) d = 0;
                else if (a > b) d = _dist[a][b];
                else d = _dist[b][a];
                Kt += d;
                if (_pop_i[i] == _pop_i[j]) {
                    _Ki[_pop_i[i]] += d;
                }
                else {
                    _Kd[_n_pop * _pop_i[i] + _pop_i[j]] += d;
                }
            }
        }

        // compute statistics
        unsigned int flag = 0;
        double Ks = 0.0;
        double Hw = 0.0;
        unsigned int nw = 0;
        _ne_pop = 0;
        unsigned int n_tot = 0;
        for (unsigned int i=0; i<_n_pop; i++) {
            if (_ni[i] >= 2) {
                Ks += _ni[i] * (static_cast<double>(_Ki[i])) * 2 / (_ni[i] * (_ni[i] - 1));
                _ne_pop++;
                n_tot += _ni[i];
            }
            Hw += (static_cast<double>(_Ki[i])) * 2 / (_ni[i] * (_ni[i] - 1));
            nw++;
        }

        double Hb = 0.0;
        unsigned int nb = 0;
        for (unsigned int i=0; i<_n_pop; i++) {
            for (unsigned int j=i+1; j<_n_pop; j++) {
                Hb += static_cast<double>((_Kd[i*_n_pop + j] + _Kd[j*_n_pop + i])) / (_ni[i] * _ni[j]);
                nb++;
            }
        }

        if (nw > 0 && Hb > 0 && nb > 0) {
            flag += 1;
            _Fst = 1 - (Hw/nw) / (Hb/nb);
        }

        if (Kt > 0.0 && n_tot>0) {
            flag += 1;
            // note: if Kt is 0 (no variation), then Hb is necessarily 0 (no between-pop variation)
            // then if flag is 1 it is only because Fst cannot be computed
            Ks /= n_tot;
            Kt /= 0.5 * n_tot * (n_tot-1);
            _Kst = 1 - Ks/Kt;
        }

        return flag;
    }

    double Haplotypes::Snn() const {

        unsigned int curr;
        unsigned int nni;
        unsigned int nnt;
        unsigned int d;
        double Snn = 0.0;

        for (unsigned int i=0; i<_nsi; i++) {
            if (_map[i] == MISSING) continue;

            curr = UNKNOWN;
            nni = 0;
            nnt = 0;

            for (unsigned int j=0; j<_nsi; j++) {
                if (_map[j] == MISSING) continue;
                if (i==j) continue;

                if (_map[i] == _map[j]) d = 0;
                else if (_map[i] > _map[j]) d = _dist[_map[i]][_map[j]];
                else d = _dist[_map[j]][_map[i]];
                if (d > curr) continue;

                if (d == curr) {
                    nnt += 1;
                    if (_pop_i[i] == _pop_i[j]) nni += 1;
                }

                else {
                    nnt = 1;
                    if (_pop_i[i] == _pop_i[j]) nni = 1;
                    else nni = 0;
                    curr = d;
                }
            }
            Snn += (double)nni/nnt;
        }
        return Snn / _nstot;
    }

    unsigned int Haplotypes::n_ing() const { return _nsi; }
    unsigned int Haplotypes::n_otg() const { return _nso; }
    unsigned int Haplotypes::n_sam() const { return _n_sam; }
    unsigned int Haplotypes::nt_hapl() const { return _nt_hapl; }
    unsigned int Haplotypes::hapl(unsigned int i, unsigned int j) const { return _hapl[i][j]; }
    unsigned int Haplotypes::map_sample(unsigned int i) const { return _map[i]; }
    unsigned int Haplotypes::n_sites() const { return _n_sites; }
    unsigned int Haplotypes::ne_ing() const { return _ne_ing; }
    unsigned int Haplotypes::ne_otg() const { return _ne_otg; }
    unsigned int Haplotypes::ni_hapl() const { return _ni_hapl; }
    unsigned int Haplotypes::ng_hapl() const { return _ng_hapl; }
    unsigned int Haplotypes::freq_i(unsigned int i) const { return _freq_i[i]; }
    unsigned int Haplotypes::freq_o(unsigned int i) const { return _freq_o[i]; }
    unsigned int Haplotypes::n_missing(unsigned int i) const { return _n_missing[i]; }
    unsigned int Haplotypes::n_mis() const { return _n_mis; }
    unsigned int Haplotypes::mis_idx(unsigned int i) const { return _mis_idx[i]; }
    unsigned int Haplotypes::n_potential(unsigned int i) const { return _potential[i][_nt_hapl]; }
    unsigned int Haplotypes::potential(unsigned int i, unsigned int j) const { return _potential[i][j]; }
    unsigned int Haplotypes::dist(unsigned int i, unsigned int j) const { return _dist[i][j]; }
    double Haplotypes::Fst() const { return _Fst; }
    double Haplotypes::Kst() const { return _Kst; }
    unsigned int Haplotypes::n_pop() const { return _n_pop; }
    unsigned int Haplotypes::ne_pop() const { return _ne_pop; }
    unsigned int Haplotypes::nstot() const { return _nstot; }
}
