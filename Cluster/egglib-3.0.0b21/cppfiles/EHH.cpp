/*
    Copyright 2013-2017 Stéphane De Mita, Mathieu Siol

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
#include "EHH.hpp"
#include <cmath>
#include <cstdlib>
#include "SiteHolder.hpp"

namespace egglib {

    EHH::EHH() {

        _par_EHH_thr = 0.0;
        _par_EHHc_thr = 0.0;
        _par_EHHS_thr = 0.0;
        _par_EHHG_thr = 0.0;
        _par_min_sam = 0;
        _opt_genotypes = false;
        _opt_crop_EHHS = false;

        _sz_K_core = 0;
        _sz_K_cur = 0;
        _sz_sam = 0;
        _sz_branches = 0;

        _K_core = 0;
        _K_cur = 0;

        _pl = 0;
        _ni = 0;
        _nsam = 0;
        _num_tot = 0;
        _ncur_tot = 0;
        _ncur_cur = NULL;
        _ncur_core = NULL;
        _num_core = NULL;
        _hap_core = NULL;
        _hap_cur = NULL;
        _hap_allele = NULL;
        _hap_origin = NULL;
        _homoz_core = NULL;
        _homoz_cur = NULL;
        _branches = NULL;

        _EHHS = 0.0;
        _EHHG = 0.0;
        _dEHHS = 0.0;
        _dEHHG = 0.0;
        _iES = 0.0;
        _iEG = 0.0;
        _flag_dEHHS = false;
        _flag_dEHHG = false;

        _EHH = NULL;
        _EHHc = NULL;
        _IHH = NULL;
        _IHHc = NULL;
        _dEHH = NULL;
        _dEHHc = NULL;
        _num_dEHH = 0;
        _num_dEHHc = 0;
        _dist = 0.0;
    }

    /* ************************************************************** */

    EHH::~EHH() {
        if (_ncur_cur) free(_ncur_cur);
        if (_ncur_core) free(_ncur_core);
        if (_num_core) free(_num_core);
        if (_hap_core) free(_hap_core);
        if (_hap_cur) free(_hap_cur);
        if (_hap_allele) free(_hap_allele);
        if (_hap_origin) free(_hap_origin);
        if (_homoz_cur) free(_homoz_cur);
        if (_homoz_core) free(_homoz_core);
        if (_EHH) free(_EHH);
        if (_EHHc) free(_EHHc);
        if (_IHH) free(_IHH);
        if (_IHHc) free(_IHHc);
        if (_dEHH) free(_dEHH);
        if (_dEHHc) free(_dEHHc);
        if (_branches) {
            for (unsigned int i=0; i<_sz_branches; i++) if (_branches[i]) free(_branches[i]);
            free(_branches);
        }
    }

    /* ************************************************************** */

    void EHH::set_core(const SiteHolder * input_site, bool genotypes,
                        double EHH_thr, double EHHc_thr,
                        double EHHS_thr, double EHHG_thr,
                        unsigned int min_freq, unsigned int min_sam,
                        bool crop) {

        // set arguments
        _par_EHHS_thr = EHHS_thr;
        _par_EHHG_thr = EHHG_thr;
        _par_EHH_thr = EHH_thr;
        _par_EHHc_thr = EHHc_thr;
        _par_min_sam = min_sam;
        _opt_genotypes = genotypes;
        _opt_crop_EHHS = crop;

        // get genotypized site (or reference to original)
        const SiteHolder * working_site;
        if (genotypes) {
            _site_geno.process(*input_site);
            working_site = & _site_geno;
        }
        else {
            working_site = input_site;
        }

        // get number of total samples (but only ingroup)
        _pl = working_site->get_ploidy();
        _ni = working_site->get_ning();
        _nsam = _pl * _ni;
        if (_nsam > _sz_sam) {
            _hap_core = (unsigned int *) realloc(_hap_core, _nsam * sizeof(unsigned int)); if (!_hap_core) throw EGGMEM;
            _hap_cur  = (unsigned int *) realloc(_hap_cur,  _nsam * sizeof(unsigned int)); if (!_hap_cur)  throw EGGMEM;
            _sz_sam = _nsam;
        }

        // process site
        _frq.setup_raw(1, 1, 0, _pl);
        _frq.setup_pop(0, 0, 0, _ni);
        _frq.process_site(*working_site);

        // analyze core alleles
        _K_core = 0;
        _K_cur = 0;
        _ncur_tot = 0;
        unsigned int allele_idx, haplo;
        for (unsigned int i=0; i<_nsam; i++) {
            allele_idx = working_site->get_straight_i(i);
            if (allele_idx == MISSING || _frq.frq_ingroup().frq_all(allele_idx) < min_freq) {
                _hap_core[i] = MISSING;
                _hap_cur[i] = MISSING;
            }
            else {

                // get index of core allele
                for (haplo=0; haplo<_K_core; haplo++) {
                    if (allele_idx == _hap_allele[haplo]) break; // existing haplotype
                }

                // new core allele
                if (haplo == _K_core) {
                    _K_core++;

                    if (_K_core > _sz_K_core) {
                        _ncur_core = (unsigned int *) realloc(_ncur_core, _K_core * sizeof(unsigned int));    if (!_ncur_core) throw EGGMEM;
                        _num_core = (unsigned int *) realloc(_num_core, _K_core * sizeof(unsigned int));      if (!_num_core) throw EGGMEM;
                        _EHH = (double *) realloc(_EHH, _K_core * sizeof(double));                            if (!_EHH) throw EGGMEM;
                        _EHHc = (double *) realloc(_EHHc, _K_core * sizeof(double));                          if (!_EHHc) throw EGGMEM;
                        _IHH = (double *) realloc(_IHH, _K_core * sizeof(double));                            if (!_IHH) throw EGGMEM;
                        _IHHc = (double *) realloc(_IHHc, _K_core * sizeof(double));                          if (!_IHHc) throw EGGMEM;
                        _dEHH = (double *) realloc(_dEHH, _K_core * sizeof(double));                          if (!_dEHH) throw EGGMEM;
                        _dEHHc = (double *) realloc(_dEHHc, _K_core * sizeof(double));                        if (!_dEHHc) throw EGGMEM;
                        _sz_K_core = _K_core;
                    }

                    // K_cur is initialized to be equal to K_core
                    _K_cur++;
                    if (_K_cur > _sz_K_cur) {
                        _ncur_cur = (unsigned int *) realloc(_ncur_cur, _K_cur * sizeof(unsigned int));       if (!_ncur_cur) throw EGGMEM;
                        _hap_allele = (unsigned int *) realloc(_hap_allele, _K_cur * sizeof(unsigned int));   if (!_hap_allele) throw EGGMEM;
                        _hap_origin = (unsigned int *) realloc(_hap_origin, _K_cur * sizeof(unsigned int));   if (!_hap_origin) throw EGGMEM;
                        _homoz_core = (bool *) realloc(_homoz_core, _K_cur * sizeof(bool));                   if (!_homoz_core) throw EGGMEM;
                        _homoz_cur = (bool *) realloc(_homoz_cur, _K_cur * sizeof(bool));                     if (!_homoz_cur) throw EGGMEM;
                        _sz_K_cur = _K_cur;
                    }

                    _ncur_core[haplo] = 0;
                    _ncur_cur[haplo] = 0;
                    _hap_allele[haplo] = allele_idx;
                    _hap_origin[haplo] = haplo;
                }

                // increment frequency
                _ncur_tot++;
                _ncur_core[haplo]++;
                _ncur_cur[haplo]++;
                _hap_core[i] = haplo;
                _hap_cur[i] = haplo;
            }
        }

        // initialise statistics
        for (unsigned int haplo=0; haplo<_K_core; haplo++) {
            _num_tot = _ncur_tot;
            _num_core[haplo] = _ncur_core[haplo];
            _EHH[haplo] = UNDEF;
            _EHHc[haplo] = UNDEF;
            _dEHH[haplo] = UNDEF;
            _dEHHc[haplo] = UNDEF;
            if (genotypes) {
                _homoz_core[haplo] = _site_geno.homoz(_hap_allele[haplo]);
                _homoz_cur[haplo] = _homoz_core[haplo];
            }
            _IHH[haplo] = (_ncur_core[haplo] == 0) ? UNDEF : 0.0;
            _IHHc[haplo] = (_ncur_core[haplo] == _ncur_tot) ? UNDEF : 0.0;
        }

        _dist = 0.0;
        _iES = 0.0;
        _iEG = 0.0;
        _EHHS = 1.0;
        _dEHHS = UNDEF;
        _dEHHG = UNDEF;
        _flag_dEHHS = false;
        _flag_dEHHG = false;
        _num_dEHH = 0;
        _num_dEHHc = 0;

        // process the core site as if it was distant to compute statistics
        load_distant(input_site, 0.0);
    }

    /* ************************************************************** */

    void EHH::load_distant(const SiteHolder * input_site, double distance) {

        if (distance < _dist) throw EggArgumentValueError("distant sites must be loaded with increasing distance");

        // get genotypized site (or reference to original)
        const SiteHolder * working_site;
        if (_opt_genotypes) {
            _site_geno.process(*input_site);
            working_site = & _site_geno;
        }
        else {
            working_site = input_site;
        }

        // check consistency of passed instance
        if (_pl != working_site->get_ploidy()) throw EggArgumentValueError("distant site is not consistent with core site");
        if (_ni != working_site->get_ning()) throw EggArgumentValueError("distant site is not consistent with core site");

        // initialize current haplotype frequency
        for (unsigned int i=0; i<_K_cur; i++) _ncur_cur[i] = 0;

        // process site
        _frq.process_site(*working_site);

        // analyse haplotypic structure
        unsigned int nbr = 0;
        for (unsigned int i=0; i<_nsam; i++) {

            // don't consider individuals that are already dropped
            if (_hap_cur[i] == MISSING) continue;
            unsigned int allele_idx = working_site->get_straight_i(i);

            // if missing allele, decrement it
            if (allele_idx == MISSING) {
                _hap_cur[i] = MISSING;
                _ncur_core[_hap_core[i]]--;
                _ncur_tot--;
                continue;
            }

            // if not missing allele
            unsigned int h = _hap_cur[i];

            // first allele at this site for the previous haplotype
            if (_ncur_cur[h] == 0) {
                _ncur_cur[h] = 1;
                _hap_allele[h] = allele_idx;
                if (_opt_genotypes && _homoz_cur[h] == true && _site_geno.homoz(allele_idx) == false) _homoz_cur[h] = false;
            }

            else {
                // not first, but same allele
                if (allele_idx == _hap_allele[h]) {
                    _ncur_cur[h]++;
                }

                // not same allele than reference (ie first) individual for this haplotype
                else {

                    // maybe one of already known haplotypes (screen breaks)
                    unsigned int j;
                    for (j=0; j<nbr; j++) {
                        if (_branches[j][0] == h && _branches[j][1] == allele_idx) {
                            h = _branches[j][2];
                            _ncur_cur[h]++;
                            _hap_cur[i] = h;
                            break;
                        }
                    }

                    // new haplotype
                    if (j == nbr) {
                        _K_cur++;
                        if (_K_cur > _sz_K_cur) {
                            _ncur_cur = (unsigned int *) realloc(_ncur_cur, _K_cur * sizeof(unsigned int));        if (!_ncur_cur) throw EGGMEM;
                            _hap_allele = (unsigned int *) realloc(_hap_allele, _K_cur * sizeof(unsigned int));    if (!_hap_allele) throw EGGMEM;
                            _hap_origin = (unsigned int *) realloc(_hap_origin, _K_cur * sizeof(unsigned int));    if (!_hap_origin) throw EGGMEM;
                            _homoz_core = (bool *) realloc(_homoz_core, _K_cur * sizeof(bool));                    if (!_homoz_core) throw EGGMEM;
                            _homoz_cur = (bool *) realloc(_homoz_cur, _K_cur * sizeof(bool));                      if (!_homoz_cur) throw EGGMEM;
                            _sz_K_cur = _K_cur;
                        }
                        _hap_allele[_K_cur-1] = allele_idx;
                        _hap_origin[_K_cur-1] = _hap_origin[h];
                        _ncur_cur[_K_cur-1] = 1;
                        _hap_cur[i] = _K_cur-1;
                        if (_opt_genotypes) {
                            _homoz_core[_K_cur-1] = _homoz_core[_hap_origin[h]];
                            if (_homoz_cur[h] && _site_geno.homoz(allele_idx)) _homoz_cur[_K_cur-1] = true;
                            else _homoz_cur[_K_cur-1] = false;
                        }

                        // add breakpoint
                        nbr++;
                        if (nbr > _sz_branches) {
                            _branches = (unsigned int **) realloc(_branches, nbr * sizeof(unsigned int *));
                            if (!_branches) throw EGGMEM;
                            _branches[j] = (unsigned int *) malloc(3 * sizeof(unsigned int));
                            if (!_branches[j]) throw EGGMEM;
                            _sz_branches = nbr;
                        }
                        _branches[j][0] = h;
                        _branches[j][1] = allele_idx;
                        _branches[j][2] = _K_cur-1;
                    }
                }
            }
        }

        // compute EHH/EHHc and IHH/IHHc
        for (unsigned int i=0; i<_K_core; i++) {
            double ehh = 0.0;

            // only consider combination if samples available (otherwise, let the default value)
            if (_ncur_core[i] >= _par_min_sam) {

                // compute EHH value
                for (unsigned int k=0; k<_K_cur; k++) {
                    if (_hap_origin[k] == i && _ncur_cur[k] > 0) {
                        ehh += _ncur_cur[k] * (_ncur_cur[k] - 1);
                    }
                }
                ehh /= (_ncur_core[i] * (_ncur_core[i] - 1));

                if (_EHH[i] != UNDEF && ehh > _EHH[i]) throw EggRuntimeError("increasing EHH");

                // increment IHH if threshold has not been reached yet
                if (_dEHH[i] == UNDEF) {
                    if (ehh < _par_EHH_thr) {
                        _IHH[i] += (distance - _dist)*(_EHH[i]-_par_EHH_thr)/(_EHH[i]-ehh) * (_EHH[i] - _par_EHH_thr) / 2.0;
                        _dEHH[i] = distance;
                        _num_dEHH++;
                    }
                    else {
                        _IHH[i] +=  (distance - _dist) * (ehh + _EHH[i] - 2 * _par_EHH_thr) / 2.0;
                    }
                }
                _EHH[i] = ehh;
            }
            else {
                _EHH[i] = UNDEF;
            }

            // compute EHHc
            unsigned int n = _ncur_tot - _ncur_core[i];
            if (n >= _par_min_sam) {

                // compute EHHc using all other samples
                double ehhc = 0.0;
                for (unsigned int k=0; k<_K_cur; k++) {
                    if (_hap_origin[k] != i && _ncur_cur[k] > 0) ehhc += _ncur_cur[k] * (_ncur_cur[k] - 1);
                }
                ehhc /= n * (n - 1);
                if (_EHHc[i] != UNDEF && ehhc > _EHHc[i]) throw EggRuntimeError("increasing EHHc");

                // increment IHHc and monitor dEHHc
                if (_dEHHc[i] == UNDEF) {

                    // monitor dEHHc and integrate if passing threshold
                    if (ehhc < _par_EHHc_thr) {
                        _IHHc[i] += (distance - _dist)*(_EHHc[i]-_par_EHHc_thr)/(_EHHc[i]-ehhc) * (_EHHc[i] - _par_EHHc_thr) / 2.0;
                        _dEHHc[i] = distance;
                        _num_dEHHc++;
                    }

                    // if threshold not reached, integrate
                    else {
                        _IHHc[i] +=  (distance - _dist) * (ehhc + _EHHc[i] - 2 * _par_EHHc_thr) / 2.0;
                    }
                }
                _EHHc[i] = ehhc;
            }
            else {
                _EHHc[i] = UNDEF;
            }
        }

        // compute EHHS and iES
        if (_ncur_tot >= _par_min_sam) {
            double ehhs = 0.0;
            for (unsigned int k=0; k<_K_cur; k++) {
                ehhs += _ncur_cur[k] * _ncur_cur[k];
            }
            ehhs = 1 - ehhs / (_ncur_tot * _ncur_tot);
            ehhs = 1 - ehhs * _ncur_tot / (_ncur_tot - 1);

            double den = 0.0;
            for (unsigned int k=0; k<_K_core; k++) {
                den += _ncur_core[k] * _ncur_core[k];
            }
            den = 1 - _ncur_tot / (_ncur_tot - 1.0) * (1 - den / (_ncur_tot * _ncur_tot));

            ehhs /= den;

            if (_opt_crop_EHHS && ehhs < _par_EHHS_thr) ehhs = 0.0;
            if (_EHHS != UNDEF && ehhs > _EHHS) throw EggRuntimeError("increasing EHHS");

            // increment iES and monitor dEHHS
            if (_dEHHS == UNDEF) {
                if (ehhs < _par_EHHS_thr) {
                    _iES += (distance - _dist)*(_EHHS-_par_EHHS_thr)/(_EHHS-ehhs) * (_EHHS - _par_EHHS_thr) / 2.0;
                    _dEHHS = distance;
                    _flag_dEHHS = true;
                }
                else {
                    _iES += (distance - _dist) * (ehhs + _EHHS - 2 * _par_EHHS_thr) / 2.0;
                }
            }
            _EHHS = ehhs;
        }
        else {
            _EHHS = UNDEF;
        }

        // compute EHHS for genotypic data
        if (_opt_genotypes) {
            if (_ncur_tot >= _par_min_sam) {
                double ehhg = 0.0;
                for (unsigned int i=0; i<_K_cur; i++) {
                    if (_homoz_cur[i]) ehhg += _ncur_cur[i];
                }
                double den = 0.0;
                unsigned int curcoretot = 0;
                for (unsigned int i=0; i<_K_core; i++) {
                    if (_homoz_core[i]) den += _ncur_core[i];
                    curcoretot += _ncur_core[i];
                }

                ehhg = ehhg / _ncur_tot / (den / curcoretot);
                if (_opt_crop_EHHS && _EHHG < _par_EHHG_thr) ehhg = 0.0;
    
                // increment iES and monitor dEHHS
                if (_dEHHG == UNDEF) {
                    if (ehhg < _par_EHHG_thr) {
                        _iEG += (distance - _dist)*(_EHHG-_par_EHHG_thr)/(_EHHG-ehhg) * (_EHHG - _par_EHHG_thr) / 2.0;
                        _dEHHG = distance;
                        _flag_dEHHG = true;
                    }
                    else {
                        _iEG += (distance - _dist) * (ehhg + _EHHG - 2 * _par_EHHG_thr) / 2.0;
                    }
                }
                _EHHG = ehhg;

            }
        }
        _dist = distance;
    }

    /* ************************************************************** */

    unsigned int EHH::K_core() const {
        return _K_core;
    }

    /* ************************************************************** */

    unsigned int EHH::K_cur() const {
        return _K_cur;
    }

    /* ************************************************************** */

    unsigned int EHH::num_avail_tot() const {
        return _ncur_tot;
    }

    /* ************************************************************** */

    unsigned int EHH::num_avail_core(unsigned int i) const {
        return _ncur_core[i];
    }

    /* ************************************************************** */

    unsigned int EHH::num_avail_cur(unsigned int i) const {
        return _ncur_cur[i];
    }

    /* ************************************************************** */

    double EHH::EHHS() const {
        return _EHHS;
    }

    /* ************************************************************** */

    double EHH::EHHG() const {
        return _EHHG;
    }

    /* ************************************************************** */

    double EHH::EHHi(unsigned int haplotype) const {
        return _EHH[haplotype];
    }

    /* ************************************************************** */

    double EHH::EHHc(unsigned int haplotype) const {
        return _EHHc[haplotype];
    }

    /* ************************************************************** */

    double EHH::rEHH(unsigned int haplotype) const {
        if (_EHH[haplotype] == UNDEF || _EHHc[haplotype] == UNDEF || _EHHc[haplotype] == 0.0) return UNDEF;
        else return _EHH[haplotype] / _EHHc[haplotype];
    }

    /* ************************************************************** */

    double EHH::IHH(unsigned int haplotype) const {
        return _IHH[haplotype];
    }

    /* ************************************************************** */

    double EHH::IHHc(unsigned int haplotype) const {
        return _IHHc[haplotype];
    }

    /* ************************************************************** */

    double EHH::iHS(unsigned int haplotype) const {
        if (_IHHc[haplotype] == 0 ||
            _IHH[haplotype] == 0 ||
            _IHHc[haplotype] == UNDEF ||
            _IHH[haplotype] == UNDEF) return UNDEF;
        else return log(_IHHc[haplotype] / _IHH[haplotype]);
    }

    /* ************************************************************** */

    double EHH::iES() const {
        return _iES;
    }

    /* ************************************************************** */

    double EHH::iEG() const {
        return _iEG;
    }

    /* ************************************************************** */

    unsigned int EHH::num_EHH_done() const {
        return _num_dEHH;
    }

    /* ************************************************************** */

    unsigned int EHH::num_EHHc_done() const {
        return _num_dEHHc;
    }

    /* ************************************************************** */

    bool EHH::flag_EHHS_done() const {
        return _flag_dEHHS;
    }

    /* ************************************************************** */

    bool EHH::flag_EHHG_done() const {
        return _flag_dEHHG;
    }

    /* ************************************************************** */

    double EHH::dEHHS() const {
        return _dEHHS;
    }

    /* ************************************************************** */

    double EHH::dEHHG() const {
        return _dEHHG;
    }

    /* ************************************************************** */

    double EHH::dEHH(unsigned int haplotype) const {
        return _dEHH[haplotype];
    }

    /* ************************************************************** */

    double EHH::dEHHc(unsigned int haplotype) const {
        return _dEHHc[haplotype];
    }

    /* ************************************************************** */

    double EHH::dEHH_max() const {
        if (_K_core == 0) return UNDEF;
        double max = -1.0;
        for (unsigned int i=0; i<_K_core; i++) {
            if (_dEHH[i] == UNDEF) return UNDEF;
            if (_dEHH[i] > max) max = _dEHH[i];
        }
        return max;
    }

    /* ************************************************************** */

    double EHH::dEHH_mean() const {
        if (_K_core == 0 || _ncur_tot == 0) return UNDEF;
        double acc = 0.0;
        for (unsigned int i=0; i<_K_core; i++) {
            if (_dEHH[i] == UNDEF) return UNDEF;
            acc += _dEHH[i] * _num_core[i];
        }
        return acc / _num_tot;
    }
}
