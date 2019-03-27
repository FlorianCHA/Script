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

#include <cstdlib>
#include "egglib.hpp"
#include "CodingSite.hpp"
#include "GeneticCode.hpp"
#include "SiteHolder.hpp"

namespace egglib {

    CodingSite::CodingSite() {
        _init();
    }

    //////

    CodingSite::~CodingSite() {
        _free();
    }

    //////

    void CodingSite::_init() {
        _ni = 0;
        _no = 0;
        _pl = 0;
        _nseff = 0;
        _nseffo = 0;
        _NSsites = 0.0;
        _Ssites = 0.0;
        _nstop = 0;
        flags = NULL;
        _nall_c = 0;
        _nall = 0;
        _aminoacids = (SiteHolder *) new(std::nothrow) SiteHolder(0);
        if (!_aminoacids) throw EGGMEM;
        _codons = (SiteHolder *) new(std::nothrow) SiteHolder(0);
        if (!_codons) throw EGGMEM;
    }

    //////

    void CodingSite::_free() {
        if (_nall_c > 0) for (unsigned int i=0; i<_nall_c-1; i++) if (flags[i]) free(flags[i]);
        if (flags) free(flags);
        if (_codons) delete _codons;
        if (_aminoacids) delete _aminoacids;
    }

    //////

    unsigned int CodingSite::_check_allele(int allele) {
        switch (allele) {
            case 'T':
            case 't': return 0;
            case 'C':
            case 'c': return 1;
            case 'A':
            case 'a': return 2;
            case 'G':
            case 'g': return 3;
        }
        return MISSING;
    }

    //////

    bool CodingSite::_process_helper(const GeneticCode& code, bool skipstop, int a1, int a2, int a3) {
        unsigned int all1 = _check_allele(a1);
        unsigned int all2 = _check_allele(a2);
        unsigned int all3 = _check_allele(a3);
        if (all1 == MISSING || all2 == MISSING || all3 == MISSING) return false;
        _codon1 = code.int2codon(all1, all2, all3);
        _amino = code.aminoacid(_codon1);
        if (_amino == '*') {
            _nstop++;
            if (skipstop) return false;
        }
        _nseff++;
        _NSsites += code.NSsites(_codon1, skipstop);
        _Ssites += code.Ssites(_codon1, skipstop);
        return true;
    }

    void CodingSite::_process_helper2(const GeneticCode& code, unsigned int all) {
        unsigned int codon2;
        if (all == _nall) {
            _nall++;
            if (_nall > _nall_c) {
                flags = (int **) realloc(flags, (_nall-1) * sizeof(int *));
                if (!flags) throw EGGMEM;
                for (unsigned int i=_nall_c>0?_nall_c-1:0; i<_nall-1; i++) flags[i] = NULL;
                for (unsigned int i=0; i<_nall-1; i++) {
                    flags[i] = (int *) realloc(flags[i], (_nall-i-1) * sizeof(int));
                    if (!flags[i]) throw EGGMEM;
                }
                _nall_c = _nall;
            }

            for (unsigned int k=0; k<_nall-1; k++) {
                codon2 = _codons->get_allele(k);
                flags[k][_nall-k-2] = (code.aminoacid(_codon1) != code.aminoacid(codon2)) ? 1 : 0;
                unsigned int tmp = 0;
                if (GeneticCode::diff1(_codon1, codon2)) {
                    flags[k][_nall-k-2] |= 2;
                    tmp += 1;
                }
                if (GeneticCode::diff2(_codon1, codon2)) {
                    flags[k][_nall-k-2] |= 4;
                    tmp += 1;
                }
                if (GeneticCode::diff3(_codon1, codon2)) {
                    flags[k][_nall-k-2] |= 8;
                    tmp += 1;
                }
                flags[k][_nall-k-2] |= (tmp << 4);
            }
        }
    }

    bool CodingSite::process(const SiteHolder& site1, const SiteHolder& site2,
                             const SiteHolder& site3, const GeneticCode& code,
                             bool skipstop, unsigned int max_missing) {

        _pl = site1.get_ploidy();
        _ni = site1.get_ning();
        _no = site1.get_nout();
        _nseff = 0;
        _nseffo = 0;
        unsigned int missing = 0;
        _NSsites = 0;
        _Ssites = 0;
        _nstop = 0;
        _nall = 0;
        _codons->reset(_pl);
        _codons->add_ing(_ni);
        _codons->add_otg(_no);
        _aminoacids->reset(_pl);
        _aminoacids->add_ing(_ni);
        _aminoacids->add_otg(_no);

        // process ingroup
        for (unsigned int i=0; i<_ni; i++) {
            for (unsigned int j=0; j<_pl; j++) {
                if (_process_helper(
                        code, skipstop,
                        site1.get_allele(site1.get_i(i, j)),
                        site2.get_allele(site2.get_i(i, j)),
                        site3.get_allele(site3.get_i(i, j))) == false) {
                    missing++;
                    if (missing > max_missing) return false;
                    _codons->load_ing(i, j, MISSINGDATA);
                    _aminoacids->load_ing(i, j, MISSINGDATA);
                    continue;
                }
                _codons->load_ing(i, j, static_cast<int>(_codon1));
                _aminoacids->load_ing(i, j, static_cast<int>(_amino));
                _process_helper2(code, _codons->get_i(i, j));
            }
        }

        // process outgroup
        for (unsigned int i=0; i<_no; i++) {
            for (unsigned int j=0; j<_pl; j++) {
                if (_process_helper(
                        code, skipstop,
                        site1.get_allele(site1.get_o(i, j)),
                        site2.get_allele(site2.get_o(i, j)),
                        site3.get_allele(site3.get_o(i, j))) == false) {
                    missing++;
                    if (missing > max_missing) return false;
                    _codons->load_otg(i, j, MISSINGDATA);
                    _aminoacids->load_otg(i, j, MISSINGDATA);
                    continue;
                }
                _codons->load_otg(i, j, static_cast<int>(_codon1));
                _aminoacids->load_otg(i, j, static_cast<int>(_amino));
                _process_helper2(code, _codons->get_o(i, j));
            }
        }

        if (_nseff > 0) {
            _NSsites /= _nseff;
            _Ssites /= _nseff;
        }
        return true;
    }

    //////

    const SiteHolder& CodingSite::codons() const {
        return * _codons;
    }

    //////

    const SiteHolder& CodingSite::aminoacids() const {
        return * _aminoacids;
    }

    //////

    unsigned int CodingSite::ni() const { return _ni; }
    unsigned int CodingSite::no() const { return _no; }
    unsigned int CodingSite::pl() const { return _pl; }

    //////

    unsigned int CodingSite::nseff() const {
        return _nseff;
    }

    //////

    double CodingSite::NSsites() const {
        return _NSsites;
    }

    //////

    double CodingSite::Ssites() const {
        return _Ssites;
    }

    //////

    unsigned int CodingSite::nstop() const {
        return _nstop;
    }

    //////

    bool CodingSite::mutated(unsigned int i, unsigned int j, unsigned int pos) const {
        if (i < j) return (flags[i][j-i-1]) & (pos==0?2:pos==1?4:8);
        else return (flags[j][i-j-1]) & (pos==0?2:pos==1?4:8);
    }

    //////

    unsigned int CodingSite::ndiff(unsigned int i, unsigned int j) const {
        if (i < j) return (flags[i][j-i-1]) >> 4;
        else return (flags[j][i-j-1]) >> 4;
    }

    //////

    bool CodingSite::NS(unsigned int i, unsigned int j) const {
        if (i < j) return (flags[i][j-i-1]) & 1;
        else return (flags[j][i-j-1]) & 1;
    }
}
