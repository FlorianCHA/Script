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

#include <cmath>
#include <cstdlib>
#include "egglib.hpp"
#include "Diversity1.hpp"
#include "FreqHolder.hpp"
#include "SiteDiversity.hpp"

namespace egglib {

    Diversity1::Diversity1() {
        _init();
        toggle_off();
    }

    //////

    Diversity1::~Diversity1() {
        _free();
    }

    //////

    void Diversity1::reset_stats() {
        _lt = 0;
        _ls = 0;
        _lso = 0;
        _nseff = 0.;
        _nseffo = 0.;
        _nsmax = 0;
        _nsmaxo = 0;
        _S = 0;
        _Ss = 0;
        _So = 0;
        _Sso = 0;
        _eta = 0;
        _etao = 0;
        _k = 0.0;
        _ko = 0.0;
        _D = 0.0;
        _Dxy = 0.0;
        _Da = 0.0;
        _Deta = 0.;
        _Dstar = 0.0;
        _Fstar = 0.0;
        _thetaW = 0.0;
        _Pi = 0.0;
        _PiForF = 0.0;
        _Pi0 = 0.0;
        _Pi1 = 0.0;
        _nsingl = 0;
        _nsingld = 0;
        _thetaPi = 0.0;
        _thetaH = 0.0;
        _thetaL = 0.0;
        _Hns = 0.0;
        _Hsd = 0.0;
        _E = 0.0;
        _Dfl = 0.0;
        _F = 0.0;
        _pM = 999.999;
        _nM = 0;
        _flag = 0;
    }

    //////

    void Diversity1::_init() {
        _esse_c = 0;
        _esse = NULL;
        _sites = NULL;
        _sites_o = NULL;
        _singl = NULL;
        _singl_o = NULL;
        _c_sites = 0;
        _c_sites_o = 0;
        _c_singl = 0;
        _c_singl_o = 0;
        _option_multiple = false;
        _option_ns_set = UNKNOWN;
        reset_stats();
    }

    //////

    void Diversity1::_free() {
        if (_esse) free(_esse);
        if (_sites) free(_sites);
        if (_sites_o) free(_sites_o);
        if (_singl) free(_singl);
        if (_singl_o) free(_singl_o);
    }

    //////

    void Diversity1::toggle_off() {
        _flag_ori_site = false;
        _flag_ori_div = false;
        _flag_basic = false;
        _flag_site_lists = false;
    }

    //////

    void Diversity1::toggle_ori_site() {
        _flag_ori_site = true;
    }

    //////

    void Diversity1::toggle_ori_div() {
        _flag_ori_site = true;
        _flag_ori_div = true;
    }

    //////

    void Diversity1::toggle_basic() {
        _flag_basic = true;
    }

    //////

    void Diversity1::toggle_site_lists() {
        _flag_site_lists = true;
    }

    //////

    void Diversity1::set_option_multiple(bool b) {
        _option_multiple = b;
    }

    //////

    void Diversity1::set_option_ns_set(unsigned int n) {
        _option_ns_set = n;
        if (n + 1 > _esse_c) {
            _esse = (unsigned int *) realloc(_esse, (n + 1) * sizeof(unsigned int));
            if (!_esse) throw EGGMEM;
            _esse_c = n + 1;
        }
        for (unsigned int i=0; i<n+1; i++) _esse[i] = 0;
    }

    //////

    void Diversity1::load(const FreqHolder& freqs, const SiteDiversity& div, unsigned int position) {

        // general counter
        _lt++;

        // get number of samples
        unsigned int ns = freqs.frq_ingroup().nseff();
        if (ns < 2) return;

        // exclude sites with >2 alleles if requested
        unsigned int nall = div.Aing();
        if (nall > 2 && _option_multiple == false)  return;

        // this is analyzable site
        _flag |= 1;
        _nseff += ns;
        if (ns > _nsmax) _nsmax = ns;
        _ls++;

        // compute distance
        if (freqs.num_populations() == 2
            && freqs.frq_population(0).nseff() > 0
            && freqs.frq_population(1).nseff() > 0)
        {
            _flag |= 16384;
            _Dxy += div.pairdiff_inter(0, 1);
            _Pi0 += div.He_pop(0);
            _Pi1 += div.He_pop(1);
        }

        // polymorphic site
        bool pol = nall > 1;
        if (pol) {
            _eta += nall-1;

            _S++;
            _Pi += div.He();
            _k += div.pairdiff();
            _nsingl += div.S();

            if (_flag_site_lists) {
                if (_S > _c_sites) {
                    _sites = (unsigned int *) realloc(_sites, _S * sizeof(unsigned int));
                    if (!_sites) throw EGGMEM;
                    _c_sites = _S;
                }
                _sites[_S-1] = position;
            }

            if (div.S() == 1 && nall == 2) {
                _Ss++;
                if (_flag_site_lists) {
                    if (_Ss > _c_singl) {
                        _singl = (unsigned int *) realloc(_singl, _Ss * sizeof(unsigned int));
                        if (!_singl) throw EGGMEM;
                        _c_singl = _Ss;
                    }
                    _singl[_Ss-1] = position;
                }
            }
        }

        // process orientable site
        if (_flag_ori_site) {
            _flag |= 4;
            
            if (div.orientable()) {
                _flag |= 8;
                _lso++;
                _nseffo += ns;
                if (ns > _option_ns_set) throw EggArgumentValueError("value used to set the number of samples is invalid");

                // if number of samples may vary
                if (ns > _nsmaxo) {
                    if (_option_ns_set != UNKNOWN) {
                        if (ns + 1 > _esse_c) {
                            _esse = (unsigned int *) realloc(_esse, (ns + 1) * sizeof(unsigned int));
                            if (!_esse) throw EGGMEM;
                            _esse_c = ns + 1;
                        }
                        for (unsigned int i=_nsmaxo; i<ns+1; i++) _esse[i] = 0;
                    }
                    _nsmaxo = ns;
                }

                // process site if polymorphic
                if (pol==true) {
                    _So++;
                    _etao += nall - 1;
                    _ko += div.pairdiff();
                    _nsingld += div.Sd();
                    _PiForF += div.He();

                    if (_flag_site_lists) {
                        if (_So > _c_sites_o) {
                            _sites_o = (unsigned int *) realloc(_sites_o, _So * sizeof(unsigned int));
                            if (!_sites_o) throw EGGMEM;
                            _c_sites_o = _So;
                        }
                        _sites_o[_So-1] = position;
                    }

                    if (div.S() == 1 && nall == 2) {
                        _Sso++;
                        if (_flag_site_lists) {
                            if (_Sso > _c_singl_o) {
                                _singl_o = (unsigned int *) realloc(_singl_o, _Sso * sizeof(unsigned int));
                                if (!_singl_o) throw EGGMEM;
                                _c_singl_o = _Sso;
                            }
                            _singl_o[_Sso-1] = position;
                        }
                    }

                    if (_option_ns_set != UNKNOWN) {
                        for (unsigned int i=0; i<div.num_derived(); i++) {
                            _esse[static_cast<unsigned int>(0.5 + _option_ns_set * div.derived(i) / ns)]++;
                        }
                    }
                    else {
                        for (unsigned int i=0; i<div.num_derived(); i++) {
                            _esse[static_cast<unsigned int>(div.derived(i))]++;
                        }
                    }
                }
            }

            // MFDM
            unsigned int max_der = 0;
            for (unsigned int i=0; i<div.num_derived(); i++) {
                if (div.derived(i) > max_der) max_der = div.derived(i);
            }
            double half = static_cast<double>(ns) / 2;
            if (max_der == half) {
                if (1.0 < _pM) _pM = 1.0;
                _nM++;
                _flag |= 16;
            }
            else {
                if (max_der > half) {
                    double pval = 2.0 * (ns - max_der) / (ns - 1);
                    if (pval < _pM) _pM = pval;
                    _nM++;
                    _flag |= 16;
                }
            }
        }
    }

    //////

    unsigned int Diversity1::compute() {
        if (_flag_basic) _basic();
        if (_flag_ori_div) _oriented();
        return _flag;
    }

    //////

    void Diversity1::_basic() {

        // requires at least one site with valid data
        if ((_flag&1) == 0) return;
        _flag |= 128;

        // compute Dxy/Da
        if ((_flag&16384) != 0) {
            _Dxy /= _ls;
            _Da = _Dxy - (_Pi0/_ls + _Pi1/_ls) / 2;
        }

        // compute the effective number of samples
        _nseff /= _ls;

        // round to an int for use in formulae
        unsigned int ns = (unsigned int)(_nseff + 0.5);

        // computes theta Watterson
        double a1 = 0.0;
        double a2 = 0.0;
        for (unsigned int i=1; i<ns; i++) {
            a1 += 1.0/i;
            a2 += 1.0/(i*i);
        }

        double b1 = (ns+1) / (3.0 * (ns-1));
        double b2 = 2.0 * (ns * ns + ns + 3) / (9.0 * ns * (ns-1));
        double c1 = b1 - 1.0 / a1;
        double c2 = b2 - (ns+2) / (a1 * ns) + a2 / (a1*a1);
        double e1 = c1 / a1;
        double e2 = c2 / (a1*a1 + a2);

        _thetaW = _S / a1;

        // require at least one polymorphic site
        if (_S == 0 || ns < 3) return;
        _flag |= 256;

        // compute D of Tajima
        double V = e1 * _S + e2 * _S * (_S-1);
        double Veta = e1 * _eta + e2 * _eta * (_eta-1);

        // compute D* of Fu and Li
            // note: an = a1, an+1 = a1+1/n, bn = a2

        double cn = (ns <= 2) ? 1 : (2.0 * (ns * a1 - 2 * (ns - 1)) / ((ns - 1) * (ns - 2)));
        double dn = cn + (ns-2.0) / ((ns-1) * (ns-1)) + (2.0 / (ns-1)) * (1.5 - (2 * (a1 + 1.0/ns)-3) / (ns-2.0) - 1.0/ns);
        double vD = ((ns * ns / ((ns-1.0) * (ns-1.0))) * a2 + a1 * a1 * dn - 2 * ns * a1 * (a1+1) / ((ns-1.0) * (ns-1.0))) / (a1 * a1 + a2);
        double uD = (ns / (ns-1.0)) * (a1 - ns / (ns-1.0)) - vD;
        if (ns > 3) {
            _flag |= 512;
            _D = (_Pi - _thetaW) / sqrt(V);
            _Deta = (_Pi - _eta / a1) / sqrt(Veta);
            _Dstar = ((ns / (ns-1.0)) * _eta - a1 * _nsingl) / sqrt(uD * _eta + vD * _eta * _eta);
        }

        // compute F* of Fu and Li
        double vF = ((2.0*ns*ns*ns + 110*ns*ns - 255*ns + 153)
                      / (9.0*ns*ns*(ns-1))
                            + (2.0*(ns-1)*a1)/(ns*ns) - 8.0*a2/ns) / (a1*a1+a2);
        double uF = ((4.0*ns*ns + 19*ns + 3 - 12*(ns+1)*(a1+1.0/ns))/(3*ns*(ns-1)))/a1 - vF;
            /* vF and uF from:
                    SIMONSEN, K. L., G. A. CHURCHILL and C. F. AQUADRO. (1995).
                    Properties of statistical tests of neutrality for DNA
                    polymorphism data. Genetics 141: 413-429.
            */
        _Fstar = (_Pi - (ns - 1.0) * _nsingl / ns) / sqrt(uF * _eta + vF * _eta * _eta);
    }

    //////

    void Diversity1::_oriented() {

        // requires at least one site with valid data and 3 samples
        if ((_flag&8) == 0) return;
        _flag |= 32;

        // compute the effective number of samples
        _nseffo /= _lso;

        // compute old school H
        unsigned int nso = (unsigned int)(_nseffo + 0.5);
        double a1 = 0.0;
        double bn = 0.0;
        for (unsigned int i=1; i<nso; i++) {
            a1 += 1.0 / i;
            bn += 1.0 / (i*i);
        }
        double bnp1 = 0.0;
        bnp1 = bn + 1.0 / (nso * nso);

        _thetaH = 0.0;
        _thetaPi = 0.0;
        _thetaL = 0.0;
        for (unsigned int i=1; i<_nsmaxo; i++) {
            _thetaPi += i * (_nsmaxo - i) * _esse[i];
            _thetaH += i * i * _esse[i];
            _thetaL += i * _esse[i];
        }
        _thetaPi *= 2.0 / (_nsmaxo * (_nsmaxo-1));  // different of Pi (maybe because random effect of non-orientable site)
        _thetaH *= 2.0 / (_nsmaxo * (_nsmaxo-1));
        _thetaL /= _nsmaxo - 1;

        _Hns = _thetaPi - _thetaH;

        // other require So>0
        if (_So == 0 or _nseffo < 3.0) return;

        // compute standardized H
        double theta = _etao / a1;
        double theta2 =  _etao*(_etao - 1.0) / (a1 * a1 + bn);
        double VarZ = theta*(_nseffo-2)/(6*(_nseffo-1)) +
                theta2*(18*_nseffo*_nseffo*(3*_nseffo+2)*bnp1 - 
                (88*_nseffo*_nseffo*_nseffo+9*_nseffo*_nseffo-13*_nseffo+6))/
                    (9*_nseffo*(_nseffo-1)*(_nseffo-1));

        if (VarZ > 0.0) {
            _flag |= 1024;
            _Hsd = (_thetaPi - _thetaL) / sqrt(VarZ);
        }

        // compute E
        double VarE = theta*(_nseffo/(2*(_nseffo-1))-1/a1)
             +  theta2*(
                   bn/(a1*a1) + 2*bn*(_nseffo/(_nseffo-1))*(_nseffo/(_nseffo-1))
                 - 2*(_nseffo*bn-_nseffo+1)/((_nseffo-1)*a1)
                 - (3*_nseffo+1)/(_nseffo-1));

        if (VarE > 0.0) {
            _flag |= 2048;
            _E = (_thetaL - theta) / sqrt(VarE);
        }

        // compute D of Fu and Li
        double cn = (_nseffo <= 2) ? 1 :
            2.0 * (_nseffo * a1 - 2 * (_nseffo - 1)) / ((_nseffo - 1) * (_nseffo - 2));

        double vD = 1 + (a1 * a1 / (bn + a1 * a1)) * (cn - (_nseffo + 1) / (_nseffo - 1));
        double uD = a1 - 1 - vD;

        double VarDfl = uD * _etao + vD * _etao * _etao;
        if (VarDfl > 0.0) {
            _flag |= 4096;
            _Dfl = (_etao - a1 * _nsingld) / sqrt(VarDfl);
        }

        // compute F of Fu and Li
        double vF = (cn + 2 * (_nseffo * _nseffo + _nseffo + 3) / (9 * _nseffo * (_nseffo-1))
                        - 2 / (_nseffo-1))
                                / (a1 * a1 + bn);
        double uF = (1 + (_nseffo+1) / (3 * (_nseffo-1))
                    - (4.0 * (_nseffo+1) / ((_nseffo-1) * (_nseffo-1)))
                            * (a1+1.0/nso - 2 * _nseffo / (_nseffo+1))) / a1 - vF;
        double VarF = uF * _etao + vF * _etao * _etao;
        if (VarF > 0.0) {
            _flag |= 8192;
            _F = (_PiForF - _nsingld) / sqrt(VarF);
        }
    }

    //////

    unsigned int Diversity1::lt() const                    { return _lt;         }
    unsigned int Diversity1::ls() const                    { return _ls;         }
    unsigned int Diversity1::lso() const                   { return _lso;        }
    unsigned int Diversity1::S() const                     { return _S;          }
    unsigned int Diversity1::Ss() const                    { return _Ss;         }
    unsigned int Diversity1::So() const                    { return _So;         }
    unsigned int Diversity1::Sso() const                   { return _Sso;        }
    unsigned int Diversity1::nsingld() const               { return _nsingld;    }
    unsigned int Diversity1::eta() const                   { return _eta;        }
    unsigned int Diversity1::etao() const                  { return _etao;       }
    double Diversity1::D() const                           { return _D;          }
    double Diversity1::Deta() const                        { return _Deta;       }
    double Diversity1::Dstar() const                       { return _Dstar;      }
    double Diversity1::Fstar() const                       { return _Fstar;      }
    double Diversity1::Pi() const                          { return _Pi;         }
    double Diversity1::thetaW() const                      { return _thetaW;     }
    double Diversity1::Dxy() const                         { return _Dxy;        }
    double Diversity1::Da() const                          { return _Da;         }
    double Diversity1::nseff() const                       { return _nseff;      }
    double Diversity1::nseffo() const                      { return _nseffo;     }
    unsigned int Diversity1::nsmax() const                 { return _nsmax;      }
    unsigned int Diversity1::nsmaxo() const                { return _nsmaxo;     }
    double Diversity1::thetaH() const                      { return _thetaH;     }
    double Diversity1::thetaPi() const                     { return _thetaPi;    }
    double Diversity1::thetaL() const                      { return _thetaL;     }
    double Diversity1::Hns() const                         { return _Hns;        }
    double Diversity1::Hsd() const                         { return _Hsd;        }
    double Diversity1::E() const                           { return _E;          }
    double Diversity1::Dfl() const                         { return _Dfl;        }
    double Diversity1::F() const                           { return _F;          }
    double Diversity1::pM() const                          { return _pM;         }
    unsigned int Diversity1::nM() const                    { return _nM;         }
    unsigned int Diversity1::site(unsigned int i) const    { return _sites[i];   }
    unsigned int Diversity1::site_o(unsigned int i) const  { return _sites_o[i]; }
    unsigned int Diversity1::singl(unsigned int i) const   { return _singl[i];   }
    unsigned int Diversity1::singl_o(unsigned int i) const { return _singl_o[i]; }
}
