/*
    Copyright 2012-2013 Stéphane De Mita, Mathieu Siol
    
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
#include "Filter.hpp"
#include <cstdlib>

namespace egglib {

    Filter::Filter() {
        init();
    }

    //////

    Filter::Filter(const Filter& src) {
        init();
        copy(src);
    }

    //////

    Filter& Filter::operator=(const Filter& src) {
        copy(src);
        return *this;
   }

    //////

    Filter::~Filter() {
        free();
    }

    //////

    void Filter::reserve(unsigned int num_expl, unsigned int num_expl_rg, unsigned int num_missing, unsigned int num_synonyms, unsigned int num_missing_synonyms) {

        if (num_expl > _ne_c) {
            _exploitable = (int *) realloc(_exploitable, num_expl * sizeof(int));
            if (!_exploitable) throw EGGMEM;
            _ne_c = num_expl;
        }

        if (num_expl_rg > _ne_c_rg) {
            _exploitable_rg = (int *) realloc(_exploitable_rg, 2 * num_expl_rg * sizeof(int));
            if (!_exploitable_rg) throw EGGMEM;
            _ne_c_rg = num_expl_rg;
        }

        if (num_missing > _nm_c) {
            _missing = (int *) realloc(_missing, num_missing * sizeof(int));
            if (!_missing) throw EGGMEM;
            _nm_c = num_missing;
        }
        
        if (num_synonyms > _nsyn_c) {
            _syn = (int *) realloc(_syn, 2 * num_synonyms * sizeof(int));
            if (!_syn) throw EGGMEM;
            _nsyn_c = num_synonyms;
        }

        if (num_missing_synonyms > _n_missing_syn_c) {
            _missing_syn = (int *) realloc(_missing_syn, 2 * num_missing_synonyms * sizeof(int));
            if (!_missing_syn) throw EGGMEM;
            _n_missing_syn_c = num_missing_synonyms;
        }
    }

    //////

    void Filter::clear() {
        free();
        init();
    }

    //////

    void Filter::add_exploitable(int code) {
        _ne++;
        if (_ne > _ne_c) {
            _exploitable = (int *) realloc(_exploitable, _ne * sizeof(int));
            if (!_exploitable) throw EGGMEM;
            _ne_c = _ne;
        }
        _exploitable[_ne-1] = code;
    }

    //////

    void Filter::add_exploitable_with_alias(int code, int alias) {
        add_exploitable(code);
        _nsyn++;
        if (_nsyn > _nsyn_c) {
            _syn = (int *) realloc(_syn, 2 * _nsyn * sizeof(int));
            if (!_syn) throw EGGMEM;
            _nsyn_c = _nsyn;
        }
        _syn[2 * (_nsyn-1)] = alias;
        _syn[2 * (_nsyn-1) + 1] = code;
    }
    
    //////

    void Filter::add_exploitable_range(int first, int last) {
        _ne_rg++;
        if (_ne_rg > _ne_c_rg) {
            _exploitable_rg = (int *) realloc(_exploitable_rg, 2 * _ne_rg * sizeof(int));
            if (!_exploitable_rg) throw EGGMEM;
            _ne_c_rg = _ne_rg;
        }
        _exploitable_rg[2*(_ne_rg-1)] = first;
        _exploitable_rg[2*(_ne_rg-1) + 1] = last;
    }

    //////

    void Filter::add_missing(int code) {
        _nm++;
        if (_nm > _nm_c) {
            _missing = (int*) realloc(_missing, _nm*sizeof(int));
            if (!_missing) throw EGGMEM;
            _nm_c = _nm;
        }
        _missing[_nm-1] = code;
    }

    //////

    void Filter::add_missing_with_alias(int code, int alias) {
        add_missing(code);
        _n_missing_syn++;
        if (_n_missing_syn > _n_missing_syn_c) {
            _missing_syn = (int *) realloc(_missing_syn, 2 * _n_missing_syn * sizeof(int));
            if (!_missing_syn) throw EGGMEM;
            _n_missing_syn_c = _n_missing_syn;
        }
        _missing_syn[2 * (_n_missing_syn-1)] = alias;
        _missing_syn[2 * (_n_missing_syn-1) + 1] = code;
    }
    
    //////

    bool Filter::is_exploitable(int& code) const {

        // liberal case
        if (_ne+_ne_rg==0) return true;
        
        // check if match in lonely codes
        for (unsigned int i=0; i<_ne; i++) {
            if (code == _exploitable[i]) return true;
        }

        // check if match in code ranges
        for (unsigned int i=0; i<_ne_rg; i++) {
            if (code >= _exploitable_rg[2*i] && code <= _exploitable_rg[2*i+1]) return true;
        }
        
        // check if match in synonyms
        for (unsigned int i=0; i<_nsyn; i++) {
            if (code == _syn[2*i]) {
                code = _syn[2*i + 1];
                return true;
            }
        }

        // no match
        return false;
    }

    //////

    bool Filter::is_missing(int& code) const {
        
        // check if match
        for (unsigned int i=0; i<_nm; i++) {
            if (code == _missing[i]) return true;
        }
        
        // check if match in synonyms
        for (unsigned int i=0; i<_n_missing_syn; i++) {
            if (code == _missing_syn[2*i]) {
                code = _missing_syn[2*i + 1];
                return true;
            }
        }
        
        // no match
        return false;
    }

    //////

    int Filter::check(int code, bool& flag) const {
        if (_ne+_ne_rg==0) return code;
        for (unsigned int i=0; i<_ne; i++) if (code == _exploitable[i]) return code;
        for (unsigned int i=0; i<_ne_rg; i++) if (code >= _exploitable_rg[2*i] && code <= _exploitable_rg[2*i+1]) return code;
        for (unsigned int i=0; i<_nsyn; i++) if (code == _syn[2*i]) return _syn[2*i + 1];
        for (unsigned int i=0; i<_nm; i++) if (code == _missing[i]) return MISSINGDATA;
        for (unsigned int i=0; i<_n_missing_syn; i++) if (code == _missing_syn[2*i]) return MISSINGDATA;
        flag = true;
        return MISSINGDATA;
    }

    //////

    void Filter::init() {
        _exploitable = NULL;
        _ne = 0;
        _ne_c = 0;
        _exploitable_rg = NULL;
        _ne_rg = 0;
        _ne_c_rg = 0;
        _missing = NULL;
        _nm = 0;
        _nm_c = 0; 
        _syn = NULL;
        _nsyn = 0;
        _nsyn_c = 0;
        _missing_syn = NULL;
        _n_missing_syn = 0;
        _n_missing_syn_c = 0;
    }

    //////

    void Filter::copy(const Filter& src) {

        reserve(src._ne, src._ne_rg, src._nm, src._nsyn, src._n_missing_syn);
        _ne = src._ne;
        _ne_rg = src._ne_rg;
        _nm = src._nm;
        _nsyn = src._nsyn;
        _n_missing_syn = src._n_missing_syn;

        for (unsigned int i=0; i<src._ne; i++) {
            _exploitable[i] = src._exploitable[i];
        }
        
        for (unsigned int i=0; i<src._ne_rg; i++) {
            _exploitable_rg[2*i] = src._exploitable_rg[2*i];
            _exploitable_rg[2*i+1] = src._exploitable_rg[2*i+1];
        }

        for (unsigned int i=0; i<src._nm; i++) {
            _missing[i] = src._missing[i];
        }
        
        for (unsigned int i=0; i<src._nsyn; i++) {
            _syn[2*i] = src._syn[2*i];
            _syn[2*i+1] = src._syn[2*i+1];
        }

        for (unsigned int i=0; i<src._n_missing_syn; i++) {
            _missing_syn[2*i] = src._missing_syn[2*i];
            _missing_syn[2*i+1] = src._missing_syn[2*i+1];
        }
    }
            
    //////

    void Filter::free() {
        if (_exploitable) ::free(_exploitable);
        if (_exploitable_rg) ::free(_exploitable_rg);
        if (_missing) ::free(_missing);
        if (_syn) ::free(_syn);
        if (_missing_syn) ::free(_missing_syn);
    }
}
