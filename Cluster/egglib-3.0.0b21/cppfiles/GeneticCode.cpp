/*
    Copyright 2013,2016 Stéphane De Mita, Mathieu Siol

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
#include "GeneticCode.hpp"
#include "SiteHolder.hpp"

namespace egglib {

    unsigned int GeneticCode::codon(char first, char second, char third) {

        unsigned int a, b, c;

        switch (first) {
            case 'T': a = 0; break;
            case 'C': a = 1; break;
            case 'A': a = 2; break;
            case 'G': a = 3; break;
            case 't': a = 0; break;
            case 'c': a = 1; break;
            case 'a': a = 2; break;
            case 'g': a = 3; break;
            default: return UNKNOWN;
        }

        switch (second) {
            case 'T': b = 0; break;
            case 'C': b = 1; break;
            case 'A': b = 2; break;
            case 'G': b = 3; break;
            case 't': b = 0; break;
            case 'c': b = 1; break;
            case 'a': b = 2; break;
            case 'g': b = 3; break;
            default: return UNKNOWN;
        }

        switch (third) {
            case 'T': c = 0; break;
            case 'C': c = 1; break;
            case 'A': c = 2; break;
            case 'G': c = 3; break;
            case 't': c = 0; break;
            case 'c': c = 1; break;
            case 'a': c = 2; break;
            case 'g': c = 3; break;
            default: return UNKNOWN;
        }

        return int2codon(a, b, c);
    }

    //////

    const char * GeneticCode::name() const {
        return _names[_index];
    }

    //////

    GeneticCode::GeneticCode() {
        _shift = 0;
        _index = 0;
    }

    //////


    GeneticCode::GeneticCode(unsigned int index) {
        set_code(index);
    }

    //////

    unsigned int GeneticCode::get_code() const {
        return _code;
    }

    //////

    void GeneticCode::set_code(unsigned int index) {
        _index = index;
        _code = _codes[index];
        _shift = index * 64;
    }

    //////

    char GeneticCode::aminoacid(unsigned int codon) const {
        if (codon < 64) return _aa[codon + _shift];
        else return 'X';
    }

    //////

    bool GeneticCode::start(unsigned int codon) const {
        if (codon < 64) return _start[codon + _shift] == 'M';
        else return false;
    }

    //////

    double GeneticCode::NSsites(unsigned int codon, bool ignorestop) const {
        return ignorestop? _NS2[codon + _shift] : _NS1[codon + _shift];
    }

    //////

    double GeneticCode::Ssites(unsigned int codon, bool ignorestop) const {
        return ignorestop? _S2[codon + _shift] : _S1[codon + _shift];
    }

    //////

    double GeneticCode::NSsites(const SiteHolder& site1, const SiteHolder& site2, const SiteHolder& site3, unsigned int& num_samples, bool ignorestop, int A, int C, int G, int T) const {
        num_samples = 0;
        unsigned int ni = site1.get_ning();
        unsigned int pl = site1.get_ploidy();
        double NS = 0.0;
        unsigned int a, b, c;
        int all;

        for (unsigned int i=0; i<ni; i++) {
            for (unsigned int j=0; j<pl; j++) {
                all = site1.get_i(i, j);
                if (all == T) a = 0;
                else if (all == C) a = 1;
                else if (all == A) a = 2;
                else if (all == G) a = 3;
                else continue;
                all = site2.get_i(i, j);
                if (all == T) b = 0;
                else if (all == C) b = 1;
                else if (all == A) b = 2;
                else if (all == G) b = 3;
                else continue;
                all = site3.get_i(i, j);
                if (all == T) c = 0;
                else if (all == C) c = 1;
                else if (all == A) c = 2;
                else if (all == G) c = 3;
                else continue;
                unsigned int codon =  int2codon(a, b, c) + _shift;
                if (ignorestop && _aa[codon] == '*') continue;
                num_samples++;
                NS += ignorestop? _NS2[codon] : _NS1[codon];
            }
        }

        return NS;
    }

    //////

    double GeneticCode::Ssites(const SiteHolder& site1, const SiteHolder& site2, const SiteHolder& site3, unsigned int& num_samples, bool ignorestop, int A, int C, int G, int T) const {
        num_samples = 0;
        unsigned int ni = site1.get_ning();
        unsigned int pl = site1.get_ploidy();
        double S = 0.0;
        unsigned int a, b, c;
        int all;

        for (unsigned int i=0; i<ni; i++) {
            for (unsigned int j=0; j<pl; j++) {
                all = site1.get_i(i, j);
                if (all == T) a = 0;
                else if (all == C) a = 1;
                else if (all == A) a = 2;
                else if (all == G) a = 3;
                else continue;
                all = site2.get_i(i, j);
                if (all == T) b = 0;
                else if (all == C) b = 1;
                else if (all == A) b = 2;
                else if (all == G) b = 3;
                else continue;
                all = site3.get_i(i, j);
                if (all == T) c = 0;
                else if (all == C) c = 1;
                else if (all == A) c = 2;
                else if (all == G) c = 3;
                else continue;
                unsigned int codon =  int2codon(a, b, c) + _shift;
                if (ignorestop && _aa[codon] == '*') continue;
                num_samples++;
                S += ignorestop? _S2[codon] : _S1[codon];
            }
        }

        return S;
    }

    //////

    double GeneticCode::NSsites(const SiteHolder& codons, unsigned int& num_samples, bool ignorestop, int A, int C, int G, int T) const {
        num_samples = 0;
        unsigned int ni = codons.get_ning();
        unsigned int pl = codons.get_ploidy();
        double NS = 0.0;
        int allele;

        for (unsigned int i=0; i<ni; i++) {
            for (unsigned int j=0; j<pl; j++) {
                allele = codons.get_i(i, j);
                if (allele < 0 || allele > 63) continue;
                if (ignorestop && _aa[allele + _shift] == '*') continue;
                num_samples++;
                NS += ignorestop? _NS2[allele + _shift] : _NS1[allele + _shift];
            }
        }

        return NS;
    }

    //////

    double GeneticCode::Ssites(const SiteHolder& codons, unsigned int& num_samples, bool ignorestop, int A, int C, int G, int T) const {
        num_samples = 0;
        unsigned int ni = codons.get_ning();
        unsigned int pl = codons.get_ploidy();
        double S = 0.0;
        int allele;

        for (unsigned int i=0; i<ni; i++) {
            for (unsigned int j=0; j<pl; j++) {
                allele = codons.get_i(i, j);
                if (allele < 0 || allele > 63) continue;
                if (ignorestop && _aa[allele + _shift] == '*') continue;
                num_samples++;
                S += ignorestop? _S2[allele + _shift] : _S1[allele + _shift];
            }
        }

        return S;
    }

    //////

    void GeneticCode::_dev_implied(char ch, unsigned int& num, int dest[]) {
        switch (ch) {
            case 'T': case 't':
                num = 1; dest[0] = 0; break;
            case 'C': case 'c':
                num = 1; dest[0] = 1; break;
            case 'A': case 'a':
                num = 1; dest[0] = 2; break;
            case 'G': case 'g':
                num = 1; dest[0] = 3; break;
            case 'Y': case 'y':
                num = 2; dest[0] = 0; dest[1] = 1; break;
            case 'W': case 'w':
                num = 2; dest[0] = 0; dest[1] = 2; break;
            case 'K': case 'k':
                num = 2; dest[0] = 0; dest[1] = 3; break;
            case 'M': case 'm':
                num = 2; dest[0] = 1; dest[1] = 2; break;
            case 'S': case 's':
                num = 2; dest[0] = 1; dest[1] = 3; break;
            case 'R': case 'r':
                num = 2; dest[0] = 2; dest[1] = 3; break;
            case 'H': case 'h':
                num = 3; dest[0] = 0; dest[1] = 1; dest[2] = 2; break;
            case 'B': case 'b':
                num = 3; dest[0] = 0; dest[1] = 1; dest[2] = 3; break;
            case 'D': case 'd':
                num = 3; dest[0] = 0; dest[1] = 2; dest[2] = 3; break;
            case 'V': case 'v':
                num = 3; dest[0] = 1; dest[1] = 2; dest[2] = 3; break;
            case 'N': case 'n':
                num = 4; dest[0] = 0; dest[1] = 1; dest[2] = 2; dest[3] = 3; break;
            default:
                num = 0;
        }
    }

    //////

    int GeneticCode::translate(int first, int second, int third, bool smart) {

        if (first == 45 && second == 45 && third == 45) return 45;

        if (smart) {
            _dev_implied(first, _num_implied_1st, _implied_1st);
            _dev_implied(second, _num_implied_2nd, _implied_2nd);
            _dev_implied(third, _num_implied_3rd, _implied_3rd);

            int cache = 'X';
            for (unsigned int i=0; i<_num_implied_1st; i++) {
                for (unsigned int j=0; j<_num_implied_2nd; j++) {
                    for (unsigned int k=0; k<_num_implied_3rd; k++) {
                        if (cache == 'X') cache = _aa[_shift + int2codon(_implied_1st[i], _implied_2nd[j], _implied_3rd[k])];
                        else if (cache != _aa[_shift + int2codon(_implied_1st[i], _implied_2nd[j], _implied_3rd[k])]) return 'X';
                    }
                }
            }
            return cache;
        }

        else {
            unsigned int a, b, c;
            switch (first) {
                case 'T': case 't': a = 0; break;
                case 'C': case 'c': a = 1; break;
                case 'A': case 'a': a = 2; break;
                case 'G': case 'g': a = 3; break;
                default: return 'X';
            }
            switch (second) {
                case 'T': case 't': b = 0; break;
                case 'C': case 'c': b = 1; break;
                case 'A': case 'a': b = 2; break;
                case 'G': case 'g': b = 3; break;
                default: return 'X';
            }
            switch (third) {
                case 'T': case 't': c = 0; break;
                case 'C': case 'c': c = 1; break;
                case 'A': case 'a': c = 2; break;
                case 'G': case 'g': c = 3; break;
                default: return 'X';
            }
            return _aa[_shift + int2codon(a, b, c)];
        }
    }

    //////

#include "GeneticCode.epp"

}
