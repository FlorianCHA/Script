/*
    Copyright 2013,2015-2016 Stéphane De Mita, Mathieu Siol

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
#include "AlleleStatus.hpp"
#include "FreqHolder.hpp"
#include <cstdlib>

namespace egglib {

    AlleleStatus::AlleleStatus() {
        reset();
    }

    //////

    AlleleStatus::~AlleleStatus() {
    }

    //////

    void AlleleStatus::reset() {
        _Sp = 0;
        _Sp_T = 0;
        _Sp_T1 = 0;
        _Spd = 0;
        _Spd_T = 0;
        _Spd_T1 = 0;
        _ShP = 0;
        _ShP_T = 0;
        _ShP_T1 = 0;
        _ShA = 0;
        _ShA_T = 0;
        _ShA_T1 = 0;
        _FxD = 0;
        _FxD_T = 0;
        _FxD_T1 = 0;
        _FxA = 0;
        _FxA_T = 0;
        _FxA_T1 = 0;
        _nsites = 0;
        _nsites_o = 0;
    }

    //////

    void AlleleStatus::process(const FreqHolder& freqs) {

        _npop = freqs.num_populations();
        _nall = freqs.num_alleles();

        // initializes  per-site counters
        _Sp = 0;
        _Spd = 0;
        _ShP = 0;
        _ShA = 0;
        _FxA = 0;
        _FxD = 0;

        _nsites++;

        // determine if site can be oriented (must be exactly 1 outgroup allele)
        bool ori = false;
        for (unsigned int all=0; all<_nall; all++) {
            if (freqs.frq_outgroup().frq_all(all) > 0) {
                if (ori == true) {
                    ori = false;
                    break;
                }
                ori = true;
                _nsites_o++;
            }
        }

        _check_Sp(freqs);
        if (ori) _check_Spd(freqs);
        _check_ShP(freqs);
        _check_ShA(freqs);
        _check_FxA(freqs);
        _check_FxD(freqs);

        _Sp_T += _Sp;               if (_Sp > 0) _Sp_T1++;
        _Spd_T += _Spd;             if (_Spd > 0) _Spd_T1++;
        _ShP_T += _ShP;             if (_ShP > 0) _ShP_T1++;
        _ShA_T += _ShA;             if (_ShA > 0) _ShA_T1++;
        _FxA_T += _FxA;             if (_FxA > 0) _FxA_T1++;
        _FxD_T += _FxD;             if (_FxD > 0) _FxD_T1++;
    }

    //////

    void AlleleStatus::total() {
        _Sp = _Sp_T;
        _Spd = _Spd_T;
        _ShP = _ShP_T;
        _ShA = _ShA_T;
        _FxA = _FxA_T;
        _FxD = _FxD_T;
    }

    //////

    unsigned int AlleleStatus::nsites() const {
        return _nsites;
    }

    //////

    unsigned int AlleleStatus::nsites_o() const {
        return _nsites_o;
    }

    //////

    void AlleleStatus::_check_Sp(const FreqHolder& freqs) {
        for (unsigned int all=0; all<_nall; all++) {
            for (unsigned int pop=0; pop<_npop; pop++) {

                if (freqs.frq_population(pop).frq_all(all) > 0 &&
                    freqs.frq_population(pop).frq_all(all) == freqs.frq_ingroup().frq_all(all))
                {
                    _Sp++;
                    break; // cannot be in another population
                }
            }
        }
    }

    //////

    void AlleleStatus::_check_Spd(const FreqHolder& freqs) {
        for (unsigned int all=0; all<_nall; all++) {
            for (unsigned int pop=0; pop<_npop; pop++) {

                if (freqs.frq_population(pop).frq_all(all) > 0 &&
                    freqs.frq_population(pop).frq_all(all) == freqs.frq_ingroup().frq_all(all)
                 && freqs.frq_outgroup().frq_all(all) == 0) 
                {
                    _Spd++;
                    break; // cannot be in another population
                }
            }
        }
    }

    //////

    void AlleleStatus::_check_FxA(const FreqHolder& freqs) {
        for (unsigned int all=0; all<_nall; all++) {
            unsigned int c = 0;
            for (unsigned int pop=0; pop<_npop; pop++) {
                if (freqs.frq_population(pop).nseff() > 0) {
                    if (freqs.frq_population(pop).frq_all(all) == 0) c |= 1;
                    if (freqs.frq_population(pop).frq_all(all) == freqs.frq_population(pop).nseff()) c |= 2;
                }
            }
            if (c == 3) {
                _FxA++;
            }
        }
    }

    //////

    void AlleleStatus::_check_ShA(const FreqHolder& freqs) {
        bool flag;
        for (unsigned int all=0; all<_nall; all++) {
            flag = false;
            for (unsigned int pop1=0; pop1<_npop; pop1++) {
                for (unsigned int pop2=pop1+1; pop2<_npop; pop2++) {

                    if (freqs.frq_population(pop1).frq_all(all) > 0 &&
                        freqs.frq_population(pop2).frq_all(all) > 0)
                    {
                        _ShA++;
                        flag = true;
                        break;
                    }
                }
                if (flag) break; // don't consider other pairs
            }
        }
    }

   //////

    void AlleleStatus::_check_ShP(const FreqHolder& freqs) {
        bool flag;
        for (unsigned int pop1=0; pop1<_npop; pop1++) {
            if (freqs.frq_population(pop1).nseff() == 0) continue;
            for (unsigned int pop2=pop1+1; pop2<_npop; pop2++) {
                if (freqs.frq_population(pop2).nseff() == 0) continue;
                for (unsigned int all1=0; all1<_nall; all1++) {
                    flag = false;
                    for (unsigned int all2=all1+1; all2<_nall; all2++) {
                        if (freqs.frq_population(pop1).frq_all(all1) > 0 &&
                            freqs.frq_population(pop1).frq_all(all1) < freqs.frq_population(pop1).nseff() &&
                            freqs.frq_population(pop2).frq_all(all1) > 0 &&
                            freqs.frq_population(pop2).frq_all(all1) < freqs.frq_population(pop2).nseff() &&
                            freqs.frq_population(pop1).frq_all(all2) > 0 &&
                            freqs.frq_population(pop1).frq_all(all2) < freqs.frq_population(pop1).nseff() &&
                            freqs.frq_population(pop2).frq_all(all2) > 0 &&
                            freqs.frq_population(pop2).frq_all(all2) < freqs.frq_population(pop2).nseff())
                        {
                            _ShP++;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) break;
                }
            }
        }
    }

   //////

    void AlleleStatus::_check_FxD(const FreqHolder& freqs) {
        for (unsigned int pop1=0; pop1<_npop; pop1++) {
            if (freqs.frq_population(pop1).nseff() == 0) continue;
            for (unsigned int pop2=pop1+1; pop2<_npop; pop2++) {
                if (freqs.frq_population(pop2).nseff() == 0) continue;
                for (unsigned int all1=0; all1<_nall; all1++) {
                    for (unsigned int all2=all1+1; all2<_nall; all2++) {

                        if ((freqs.frq_population(pop1).frq_all(all1) == freqs.frq_population(pop1).nseff() &&
                             freqs.frq_population(pop2).frq_all(all2) == freqs.frq_population(pop2).nseff()) ||
                            (freqs.frq_population(pop1).frq_all(all2) == freqs.frq_population(pop1).nseff() &&
                             freqs.frq_population(pop2).frq_all(all1) == freqs.frq_population(pop2).nseff()))
                        {
                            _FxD++;
                        }
                    }
                }
            }
        }
    }

    //////

    unsigned int AlleleStatus::Sp()     const { return _Sp; }
    unsigned int AlleleStatus::Sp_T1()  const { return _Sp_T1; }
    unsigned int AlleleStatus::Spd()    const { return _Spd; }
    unsigned int AlleleStatus::Spd_T1() const { return _Spd_T1; }
    unsigned int AlleleStatus::ShP()    const { return _ShP; }
    unsigned int AlleleStatus::ShP_T1() const { return _ShP_T1; }
    unsigned int AlleleStatus::ShA()    const { return _ShA; }
    unsigned int AlleleStatus::ShA_T1() const { return _ShA_T1; }
    unsigned int AlleleStatus::FxD()      const { return _FxD; }
    unsigned int AlleleStatus::FxD_T1()   const { return _FxD_T1; }
    unsigned int AlleleStatus::FxA()      const { return _FxA; }
    unsigned int AlleleStatus::FxA_T1()   const { return _FxA_T1; }
}
