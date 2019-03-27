/*
    Copyright 2012-2016 Stéphane De Mita, Mathieu Siol

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
#include "DataHolder.hpp"
#include "Random.hpp"
#include "Coalesce.hpp"
#include <cstdlib>
#include <cstring>
#include <cmath>

namespace egglib {

    //~~~~~~~~~//

    void VectorInt::copy(const VectorInt& src) {
        set_num_values(src._num);
        for (unsigned int i=0; _num; i++) _values[i] = src._values[i];
    }

    //~~~~~~~~~//

    VectorInt::VectorInt() {
        _num = 0;
        _res = 0;
        _values = NULL;
    }

    //~~~~~~~~~//

    VectorInt::~VectorInt() {
        if (_values) free(_values);
    }

    //~~~~~~~~~//

    VectorInt::VectorInt(const VectorInt& src) {
        _res = 0;
        _values = NULL;
        copy(src);
    }

    //~~~~~~~~~//

    VectorInt& VectorInt::operator=(const VectorInt& src) {
        copy(src);
        return * this;
    }

    //~~~~~~~~~//

    void VectorInt::set_num_values(unsigned int n) {
        _num = n;
        if (_num > _res) {
            _values = (int *) realloc(_values, _num * sizeof(int));
            if (!_values) throw EGGMEM;
        }
    }

    //~~~~~~~~~//

    unsigned int VectorInt::get_num_values() const {
        return _num;
    }

    //~~~~~~~~~//

    void VectorInt::set_item(unsigned int i, int value) {
        _values[i] = value;
    }

    //~~~~~~~~~//

    int VectorInt::get_item(unsigned int i) const {
        return _values[i];
    }

    //~~~~~~~~~//

    void VectorInt::clear() {
        if (_values) free(_values);
        _res = 0;
        _num = 0;
        _values = NULL;
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    void DataHolder::_init() {
        _nsi = 0;
        _nso = 0;
        _ls = 0;
        _ng = 0;
        _nsi_r = 0;
        _nso_r = 0;
        _lsi = NULL;
        _lso = NULL;
        _lsi_r = NULL;
        _lso_r = NULL;
        _lni = NULL;
        _lni_r = NULL;
        _lno = NULL;
        _lno_r = NULL;
        _ng_r = NULL;
        _datai = NULL;
        _datao = NULL;
        _namesi = NULL;
        _nameso = NULL;
        _groups_i = NULL;
        _groups_o = NULL;
    }

    //~~~~~~~~~//

    void DataHolder::_free() {

        // clean 1-dimensional tables

        if (_lni) free(_lni);
        if (_lni_r) free(_lni_r);
        if (_lno) free(_lno);
        if (_lno_r) free(_lno_r);
        if (_lsi) free(_lsi);
        if (_lso) free(_lso);
        if (_lsi_r) free(_lsi_r);
        if (_lso_r) free(_lso_r);
        if (_ng_r) free(_ng_r);

        // clean ingroup table and names

        for (unsigned int i=0; i<_nsi_r; i++) {
            if (_datai[i]) free(_datai[i]);
            if (_namesi[i]) free(_namesi[i]);
        }
        if (_datai) free(_datai);
        if (_namesi) free(_namesi);

        // clean group labels

        for (unsigned int i=0; i<_nsi_r; i++) {
            free(_groups_i[i]);
        }
        if (_groups_i) free(_groups_i);
        if (_groups_o) free(_groups_o);

        // clean outgroup table and names

        for (unsigned int i=0; i<_nso_r; i++) {
            if (_datao[i]) free(_datao[i]);
            if (_nameso[i]) free(_nameso[i]);
        }
        if (_datao) free(_datao);
        if (_nameso) free(_nameso);
    }

    //~~~~~~~~~//

    void DataHolder::_alloc_nsi(unsigned int nsi) {

        if (nsi > _nsi_r) {

            _lsi = (unsigned int * ) realloc(_lsi, nsi * sizeof(unsigned int));
            if (!_lsi) throw EGGMEM;
            _lni = (unsigned int * ) realloc(_lni, nsi * sizeof(unsigned int));
            if (!_lni) throw EGGMEM;
            _datai = (int ** ) realloc(_datai, nsi * sizeof(int * ));
            if (!_datai) throw EGGMEM;
            _namesi = (char ** ) realloc(_namesi, nsi * sizeof(char * ));
            if (!_namesi) throw EGGMEM;
            _groups_i = (unsigned int ** ) realloc(_groups_i, nsi * sizeof(unsigned int * ));
            if (!_groups_i) throw EGGMEM;

            _lsi_r = (unsigned int * ) realloc(_lsi_r, nsi * sizeof(unsigned int));
            if (!_lsi_r) throw EGGMEM;
            _lni_r = (unsigned int * ) realloc(_lni_r, nsi * sizeof(unsigned int));
            if (!_lni_r) throw EGGMEM;
            _ng_r = (unsigned int * ) realloc(_ng_r, nsi * sizeof(unsigned int));
            if (!_ng_r) throw EGGMEM;

            for (unsigned int i=_nsi_r; i<nsi; i++) {
                _datai[i] = NULL;
                _lsi_r[i] = 0;
                _namesi[i] = NULL;
                _lni_r[i] = 0;
                _groups_i[i] = NULL;
                _ng_r[i] = 0;
            }

            _nsi_r = nsi;
        }
    }

    //~~~~~~~~~//

    void DataHolder::_alloc_nso(unsigned int nso) {

        if (nso > _nso_r) {

            _lso = (unsigned int * ) realloc(_lso, nso * sizeof(unsigned int));
            if (!_lso) throw EGGMEM;
            _lso_r = (unsigned int * ) realloc(_lso_r, nso * sizeof(unsigned int));
            if (!_lso_r) throw EGGMEM;
            _lno = (unsigned int * ) realloc(_lno, nso * sizeof(unsigned int));
            if (!_lno) throw EGGMEM;
            _lno_r = (unsigned int * ) realloc(_lno_r, nso * sizeof(unsigned int));
            if (!_lno_r) throw EGGMEM;
            _datao = (int ** ) realloc(_datao, nso * sizeof(int * ));
            if (!_datao) throw EGGMEM;
            _nameso = (char ** ) realloc(_nameso, nso * sizeof(char * ));
            if (!_nameso) throw EGGMEM;
            _groups_o = (unsigned int * ) realloc(_groups_o, nso * sizeof(unsigned int));
            if (!_groups_o) throw EGGMEM;
            for (unsigned int i=_nso_r; i<nso; i++) {
                _datao[i] = NULL;
                _lso_r[i] = 0;
                _nameso[i] = NULL;
                _lno_r[i] = 0;
            }

            _nso_r = nso;
        }
    }

    //~~~~~~~~~//

    void DataHolder::_alloc_ng(unsigned int ng) {

        for (unsigned int i=0; i<_nsi; i++) {
            if (ng > _ng_r[i]) {
                _groups_i[i] = (unsigned int * ) realloc(_groups_i[i], ng * sizeof(unsigned int));
                if (!_groups_i[i]) throw EGGMEM;
                _ng_r[i] = ng;
            }
        }
    }

    //~~~~~~~~~//

    void DataHolder::_alloc_ls(unsigned int ls) {
        for (unsigned int i=0; i<_nsi; i++) _alloc_lsi(i, ls);
        for (unsigned int i=0; i<_nso; i++) _alloc_lso(i, ls);
    }

    //~~~~~~~~~//

    void DataHolder::_alloc_lsi(unsigned int i, unsigned int ls) {
        if (ls > _lsi_r[i]) {
            _datai[i] = (int * ) realloc(_datai[i], ls * sizeof(int));
            if (!_datai[i]) throw EGGMEM;
            _lsi_r[i] = ls;
        }
    }

    //~~~~~~~~~//

    void DataHolder::_alloc_lso(unsigned int i, unsigned int ls) {
        if (ls > _lso_r[i]) {
            _datao[i] = (int * ) realloc(_datao[i], ls * sizeof(int));
            if (!_datao[i]) throw EGGMEM;
            _lso_r[i] = ls;
        }
    }

    //~~~~~~~~~//

    void DataHolder::_alloc_lni(unsigned int i, unsigned int ln) {
        if (ln > _lni_r[i]) {
            _namesi[i] = (char * ) realloc(_namesi[i], ln * sizeof(char));
            if (!_namesi[i]) throw EGGMEM;
            _lni_r[i] = ln;
        }
    }

    //~~~~~~~~~//

    void DataHolder::_alloc_lno(unsigned int i, unsigned int ln) {
        if (ln > _lno_r[i]) {
            _nameso[i] = (char * ) realloc(_nameso[i], ln * sizeof(char));
            if (!_nameso[i]) throw EGGMEM;
            _lno_r[i] = ln;
        }
    }

    //~~~~~~~~~//

    void DataHolder::_copy(const DataHolder& src) {

        _is_matrix = src._is_matrix;

        // allocate arrays

        _nsi = src._nsi;
        _nso = src._nso;
        _ng = src._ng;
        _alloc_nsi(_nsi);
        _alloc_nso(_nso);
        _alloc_ng(_ng);

        // get names and groups

        for (unsigned int i=0; i<_nsi; i++) {
            set_name_i(i, src._namesi[i]);
            for (unsigned int j=0; j<_ng; j++) _groups_i[i][j] = src._groups_i[i][j];
        }
        for (unsigned int i=0; i<_nso; i++) {
            set_name_o(i, src._nameso[i]);
            _groups_o[i] = src._groups_o[i];
        }

        // get data

        if (_is_matrix) {
            _ls = src._ls;
            _alloc_ls(_ls);
            for (unsigned int j=0; j<_ls; j++) {
                for (unsigned int i=0; i<_nsi; i++) _datai[i][j] = src._datai[i][j];
                for (unsigned int i=0; i<_nso; i++) _datao[i][j] = src._datao[i][j];
            }
        }
        else {
            for (unsigned int i=0; i<_nsi; i++) {
                _lsi[i] = src._lsi[i];
                _alloc_lsi(i, src._lsi[i]);
                for (unsigned int j=0; j<_lsi[i]; j++) _datai[i][j] = src._datai[i][j];
            }
            for (unsigned int i=0; i<_nso; i++) {
                _lso[i] = src._lso[i];
                _alloc_lso(i, src._lso[i]);
                for (unsigned int j=0; j<_lso[i]; j++) _datao[i][j] = src._datao[i][j];
            }
        }
    }

    //~~~~~~~~~//

    DataHolder::DataHolder(bool is_matrix) {
        _init();
        _is_matrix = is_matrix;
    }

    //~~~~~~~~~//

    DataHolder::DataHolder(const DataHolder& src) {
        _init();
        _copy(src);
    }

    //~~~~~~~~~//

    DataHolder& DataHolder::operator=(const DataHolder& src) {
        _copy(src);
        return *this;
    }

    //~~~~~~~~~//

    DataHolder::~DataHolder() {
        _free();
    }

    //~~~~~~~~~//

    void DataHolder::set_is_matrix(bool flag) {

        // matrix -> non-matrix

        if (_is_matrix == true && flag == false) {
            for (unsigned int i=0; i<_nsi; i++) _lsi[i] = _ls;
            for (unsigned int i=0; i<_nso; i++) _lso[i] = _ls;
            _is_matrix = false;
        }

        // non-matrix -> matrix

        if (_is_matrix == false && flag == true) {
            if (_nsi > 0) _ls = _lsi[0];
            else {
                if (_nso > 0) _ls = _lso[0];
                else _ls = 0;
            }
            _is_matrix = true;
        }

        // nothing to do otherwise
    }

    //~~~~~~~~~//

    bool DataHolder::get_is_matrix() const {
        return _is_matrix;
    }

    //~~~~~~~~~//

    void DataHolder::reserve(unsigned int nsi, unsigned int nso, unsigned int ln, unsigned int ng, unsigned int ls) {
        _alloc_nsi(nsi);
        _alloc_nso(nso);
        _alloc_ng(ng);
        _alloc_ls(ls);
        for (unsigned int i=0; i<nsi; i++) _alloc_lni(i, ln);
        for (unsigned int i=0; i<nso; i++) _alloc_lno(i, ln);
    }

    //~~~~~~~~~//

    unsigned int DataHolder::get_nsam_i() const {
        return _nsi;
    }

    //~~~~~~~~~//

    void DataHolder::set_nsam_i(unsigned int nsam) {

        _alloc_nsi(nsam);

        for (unsigned int i=_nsi; i<nsam; i++) {

            set_name_i(i, "");

            if (_is_matrix) {
                _alloc_lsi(i, _ls);
                _lsi[i] = _ls;
            }
            else {
                _lsi[i] = 0;
            }

            if (_ng > _ng_r[i]) {
                _groups_i[i] = (unsigned int * ) realloc(_groups_i[i], _ng * sizeof(int));
                if (!_groups_i[i]) throw EGGMEM;
                _ng_r[i] = _ng;
            }
        }

        _nsi = nsam;
    }

    //~~~~~~~~~//

    unsigned int DataHolder::get_nsam_o() const {
        return _nso;
    }

    //~~~~~~~~~//

    void DataHolder::set_nsam_o(unsigned int nsam) {

        _alloc_nso(nsam);

        for (unsigned int i=_nso; i<nsam; i++) {

            set_name_o(i, "");

            if (_is_matrix) {
                _alloc_lso(i, _ls);
                _lso[i] = _ls;
            }

            else _lso[i] = 0;
        }

        _nso = nsam;
    }

    //~~~~~~~~~//

    unsigned int DataHolder::get_nsit() const {
        return _ls;
    }

    //~~~~~~~~~//

    unsigned int DataHolder::get_nsit_i(unsigned int sam) const {
        return _lsi[sam];
    }

    //~~~~~~~~~//

    unsigned int DataHolder::get_nsit_o(unsigned int sam) const {
        return _lso[sam];
    }

    //~~~~~~~~~//

    void DataHolder::set_nsit(unsigned int val) {
        _ls = val;
        _alloc_ls(_ls);
    }

    //~~~~~~~~~//

    void DataHolder::set_nsit_i(unsigned int sam, unsigned int val) {
        _lsi[sam] = val;
        _alloc_lsi(sam, val);
    }

    //~~~~~~~~~//

    void DataHolder::set_nsit_o(unsigned int sam, unsigned int val) {
        _lso[sam] = val;
        _alloc_lso(sam, val);
    }

    //~~~~~~~~~//

    void DataHolder::insert_sites(unsigned int pos, unsigned int num) {

        if (num == 0) return;          // nothing to do (prevents side effect with j iterator)

        // increase the data table (add num to all arrays)

        if (_is_matrix) {
            _ls += num;
            _alloc_ls(_ls);
        }
        else {
            for (unsigned int i=0; i<_nsi; i++) {
                _lsi[i] += num;
                _alloc_lsi(i, _lsi[i]);
            }
            for (unsigned int i=0; i<_nso; i++) {
                _lso[i] += num;
                _alloc_lso(i, _lso[i]);
            }
        }

        // shift sites (only if pos is not MAX)

        if (pos != MAX) {

            unsigned int ls = _ls; // changed if not matrix

            for (unsigned int i=0; i<_nsi; i++) {
                if (!_is_matrix) ls = _lsi[i];
                for (unsigned int j=ls-1; j>=pos+num; j--) {
                    _datai[i][j] =  _datai[i][j-num];
              }
            }

            for (unsigned int i=0; i<_nso; i++) {
                if (!_is_matrix) ls = _lso[i];
                for (unsigned int j=ls-1; j>=pos+num; j--) _datao[i][j] = _datao[i][j-num];
            }
        }
    }

    //~~~~~~~~~//

    void DataHolder::insert_sites_i(unsigned int sam, unsigned int pos, unsigned int num) {
        if (num == 0) return;
        _lsi[sam] += num;
        _alloc_lsi(sam, _lsi[sam]);
        for (unsigned int i=_lsi[sam]; i>pos+num; i--) _datai[sam][i-1] = _datai[sam][i-1-num];
    }

    //~~~~~~~~~//

    void DataHolder::insert_sites_o(unsigned int sam, unsigned int pos, unsigned int num) {
        if (num == 0) return;
        _lso[sam] += num;
        _alloc_lso(sam, _lso[sam]);
        for (unsigned int i=_lso[sam]; i>pos+num; i--) _datao[sam][i-1] = _datao[sam][i-1-num];
    }

    //~~~~~~~~~//

    unsigned int DataHolder::get_ngroups() const {
        return _ng;
    }

    //~~~~~~~~~//

    void DataHolder::set_ngroups(unsigned int ngrp) {
        _ng = ngrp;
        _alloc_ng(_ng);
    }

    //~~~~~~~~~//

    int DataHolder::get_i(unsigned int sam, unsigned int sit) const {
        return _datai[sam][sit];
    }

    //~~~~~~~~~//

    void DataHolder::set_i(unsigned int sam, unsigned int sit, int value) {
        _datai[sam][sit] = value;
    }

    //~~~~~~~~~//

    int DataHolder::get_o(unsigned int sam, unsigned int sit) const {
        return _datao[sam][sit];
    }

    //~~~~~~~~~//

    void DataHolder::set_o(unsigned int sam, unsigned int sit, int value) {
        _datao[sam][sit] = value;
    }

    //~~~~~~~~~//

    unsigned int DataHolder::get_group_i(unsigned int sam, unsigned int lvl) const {
        return _groups_i[sam][lvl];
    }

    //~~~~~~~~~//

    unsigned int DataHolder::get_group_o(unsigned int sam) const {
        return _groups_o[sam];
    }

    //~~~~~~~~~//

    void DataHolder::set_group_i(unsigned int sam, unsigned int lvl, unsigned int label) {
        _groups_i[sam][lvl] = label;
    }

    //~~~~~~~~~//

    void DataHolder::set_group_o(unsigned int sam, unsigned int label) {
        _groups_o[sam] = label;
    }

    //~~~~~~~~~//

    const char * DataHolder::get_name_i(unsigned int sam) const {
        return _namesi[sam];
    }

    //~~~~~~~~~//

    void DataHolder::set_name_i(unsigned int sam, const char * name) {
        _lni[sam] = strlen(name) + 1; // be careful, at some point we check empty names by _lni[i]==1 (valid_phyml_names)
        _alloc_lni(sam, _lni[sam]);
        strcpy(_namesi[sam], name);
    }

    //~~~~~~~~~//

    void DataHolder::name_appendch_i(unsigned int sam, char ch) {
        _lni[sam]++;
        _alloc_lni(sam, _lni[sam]);
        _namesi[sam][_lni[sam]-2] = ch;
        _namesi[sam][_lni[sam]-1] = '\0';
    }

    //~~~~~~~~~//

    void DataHolder::name_append_i(unsigned int sam, const char * str) {
        unsigned int cur = _lni[sam];
        _lni[sam] += strlen(str);          // might overflow because strlen returns a size_t
        _alloc_lni(sam, _lni[sam]);
        strcpy(_namesi[sam]+cur-1, str);
        _namesi[sam][_lni[sam]-1] = '\0';
    }

    //~~~~~~~~~//

    const char * DataHolder::get_name_o(unsigned int sam) const {
        return _nameso[sam];
    }

    //~~~~~~~~~//

    void DataHolder::set_name_o(unsigned int sam, const char * str) {
        _lno[sam] = strlen(str) + 1;
        _alloc_lno(sam, _lno[sam]);
        strcpy(_nameso[sam], str);
    }

    //~~~~~~~~~//

    void DataHolder::name_appendch_o(unsigned int sam, char ch) {
        _lno[sam]++;
        _alloc_lno(sam, _lno[sam]);
        _nameso[sam][_lno[sam]-2] = ch;
        _nameso[sam][_lno[sam]-1] = '\n';
    }

    //~~~~~~~~~//

    void DataHolder::name_append_o(unsigned int sam, const char * str) {
        unsigned int cur = _lno[sam];
        _lno[sam] += strlen(str);
        _alloc_lno(sam, _lno[sam]);
        strcpy(_nameso[sam]+cur-1, str);
        _nameso[sam][_lno[sam]-1] = '\n';
    }

    //~~~~~~~~~//

    void DataHolder::to_outgroup(unsigned int sam, unsigned int label) {

        // save data that must be transfered

        char * name = _namesi[sam];
        int * data = _datai[sam];
        unsigned int ln = _lni[sam];
        unsigned int ls_r = _lsi_r[sam];
        unsigned int ln_r = _lni_r[sam];

        unsigned int ls = _lsi[sam]; // only used if !_is_matrix

        // free groups (not transferred)

        if (_ng_r[sam]) ::free(_groups_i[sam]);

        // copy all pointers at previous index

        for (unsigned int i = sam; i < _nsi-1; i++) {
            _namesi[i] = _namesi[i+1];
            _datai[i] = _datai[i+1];
            _groups_i[i] = _groups_i[i+1];
            _lni[i] = _lni[i+1];
            _lni_r[i] = _lni_r[i+1];
            _lsi_r[i] = _lsi_r[i+1];
            _ng_r[i] = _ng_r[i+1];
            if (!_is_matrix) _lsi[i] = _lsi[i+1];
        }

        // set new size

        _nsi--;
        _nsi_r--;

        // allocate new outgroup sample if needed

        if (_nso == _nso_r) _alloc_nso(_nso+1);

        // replace the last item (which has been duplicated) by the one that will be overwritten by transferred data

        _namesi[_nsi] = _nameso[_nso];
        _datai[_nsi] = _datao[_nso];
        _groups_i[_nsi] = NULL;
        _lni_r[_nsi] = _lno_r[_nso];
        _lsi_r[_nsi] = _lso_r[_nso];
        _ng_r[_nsi] = 0;

        // add an outgroup item

        _nso++;

        _nameso[_nso-1] = name;
        _datao[_nso-1] = data;
        _lno[_nso-1] = ln;
        _lno_r[_nso-1] = ln_r; 
        _lso_r[_nso-1] = ls_r;

        if (!_is_matrix) _lso[_nso-1] = ls;

        // set the label

        _groups_o[_nso-1] = label;
    }

    //~~~~~~~~~//

    void DataHolder::del_sample_i(unsigned int sam) {

        if (sam == _nsi-1) {
            _nsi--;
        }

        else {
            char * name = _namesi[sam];
            int * data = _datai[sam];
            unsigned int * group = _groups_i[sam];
            unsigned int ls = _lsi[sam]; // only used if not matrix
            unsigned int ln_r = _lni_r[sam];
            unsigned int ls_r = _lsi_r[sam];
            unsigned int ng_r = _ng_r[sam];

            for (unsigned int i = sam; i < _nsi-1; i++) {
                _namesi[i] = _namesi[i+1];
                _datai[i] = _datai[i+1];
                _groups_i[i] = _groups_i[i+1];
                _lni[i] = _lni[i+1];
                _lni_r[i] = _lni_r[i+1];
                _lsi_r[i] = _lsi_r[i+1];
                _ng_r[i] = _ng_r[i+1];
            }

            _nsi--;

            _namesi[_nsi] = name;
            _datai[_nsi] = data;
            _groups_i[_nsi] = group;
            _lni_r[_nsi] = ln_r;
            _lsi_r[_nsi] = ls_r;
            _ng_r[_nsi] = ng_r;
            if (!_is_matrix) _lsi[_nsi] = ls;
        }

        if (_nsi == 0 && _nso == 0) _ls = 0;
    }

    //~~~~~~~~~//

    void DataHolder::del_sample_o(unsigned int sam) {

        if (sam == _nso-1) {
            _nso--;
        }

        else {
            char * name = _nameso[sam];
            int * data = _datao[sam];
            unsigned int group = _groups_o[sam];
            unsigned int ls = _lso[sam]; // only used if not matrix
            unsigned int ln_r = _lno_r[sam];
            unsigned int ls_r = _lso_r[sam];
            
            for (unsigned int i = sam; i < _nso-1; i++) {
                _nameso[i] = _nameso[i+1];
                _datao[i] = _datao[i+1];
                _groups_o[i] = _groups_o[i+1];
                _lno[i] = _lno[i+1];
                _lno_r[i] = _lno_r[i+1];
                _lso_r[i] = _lso_r[i+1];
            }

            _nso--;

            _nameso[_nso] = name;
            _datao[_nso] = data;
            _groups_o[_nso] = group;
            _lno_r[_nso] = ln_r;
            _lso_r[_nso] = ls_r;
            if (!_is_matrix) _lso[_nso] = ls;
        }

        if (_nsi == 0 && _nso == 0) _ls = 0;
    }

    //~~~~~~~~~//

    void DataHolder::del_sites(unsigned int start, unsigned int stop) {

        if (start>=stop) return;

        unsigned int ls;

        if (_is_matrix) {
            ls = _ls;
            if (start >= ls) return;
            if (stop > ls) stop = ls;
            _ls -= (stop-start);
        }

        for (unsigned int i=0; i<_nsi; i++) {
            if (!_is_matrix) del_sites_i(i, start, stop);
            else for (unsigned int j=0; stop+j<_ls; j++) _datai[i][start+j] = _datai[i][stop+j];
        }

        for (unsigned int i=0; i<_nso; i++) {
            if (!_is_matrix) del_sites_o(i, start, stop);
            else for (unsigned int j=0; stop+j<_ls; j++) _datao[i][start+j] = _datao[i][stop+j];
        }
    }

    //~~~~~~~~~//

    void DataHolder::del_sites_i(unsigned int sam, unsigned int start, unsigned int stop) {

        unsigned int ls = _lsi[sam];
        if (start >= ls) return;
        if (stop > ls) stop = ls;
        _lsi[sam] -= (stop-start);
            
        for (unsigned int i=0; stop+i<ls; i++) _datai[sam][start+i] = _datai[sam][stop+i];
    }

    //~~~~~~~~~//

    void DataHolder::del_sites_o(unsigned int sam, unsigned int start, unsigned int stop) {

        unsigned int ls = _is_matrix? _ls : _lso[sam];
        if (start >= ls) return;
        if (stop > ls) stop = ls;
        _lso[sam] -= (stop-start);
            
        for (unsigned int i=0; stop+i<ls; i++) _datao[sam][start+i] = _datao[sam][stop+i];
    }

    //~~~~~~~~~//

    void DataHolder::reset(bool is_matrix) {
        _nsi = 0;
        _nso = 0;
        if (is_matrix) _ls = 0;
        _ng = 0;
        _is_matrix = is_matrix;
    }

    //~~~~~~~~~//

    void DataHolder::clear(bool is_matrix) {
        _free();
        _init();
        _is_matrix = is_matrix;
    }

    //~~~~~~~~~//

    unsigned int DataHolder::find(unsigned int sam, bool of_outgroup, VectorInt& motif, unsigned int start, unsigned int stop) const {

        if (_is_matrix) {
            if (stop > _ls) stop = _ls;
        }
        else {
            if (of_outgroup) {
                if (stop > _lso[sam]) stop = _lso[sam];
            }
            else {
                if (stop > _lsi[sam]) stop = _lsi[sam];
            }
        }

        unsigned int n = motif.get_num_values();

        if (start >= stop || n==0) return MAX;

        unsigned int i=start;
        unsigned int j;

        while (true) {
            j = 0;
            while (true) {
                if (i+j == stop) return MAX;
                if (of_outgroup) {
                    if (_datao[sam][i+j] != motif.get_item(j)) break;
                }
                else {
                    if (_datai[sam][i+j] != motif.get_item(j)) break;
                }
                j++;
                if (j==n) return i;
            }
            i++;
        }
        
        throw EggRuntimeError("this point should not be reached (DataHolder::find)");
        return MAX;
    }

    //~~~~~~~~~//

    bool DataHolder::is_equal() const {
        unsigned int n;
        if (_nsi > 0) n = _lsi[0];
        else {
            if (_nsi > 0) n = _lso[0];
            else return true;
        }
        for (unsigned int i=0; i<_nsi; i++) if (_lsi[i] != n) return false;
        for (unsigned int i=0; i<_nso; i++) if (_lso[i] != n) return false;
        return true;
    }

    //~~~~~~~~~//

    bool DataHolder::valid_phyml_nt() const {
        for (unsigned int j=0; j<_ls; j++) {
            for (unsigned int i=0; i<_nsi; i++) {
                switch (_datai[i][j]) {
                    case 'A': case 'C': case 'G': case 'T':
                    case 'U': case 'M': case 'R': case 'W':
                    case 'S': case 'Y': case 'K': case 'B':
                    case 'D': case 'H': case 'V':
                    case '-': case 'N': case 'X': case '?': break;
                    default: return false;
                }
            }
            for (unsigned int i=0; i<_nso; i++) {
                switch (_datao[i][j]) {
                    case 'A': case 'C': case 'G': case 'T':
                    case 'U': case 'M': case 'R': case 'W':
                    case 'S': case 'Y': case 'K': case 'B':
                    case 'D': case 'H': case 'V':
                    case '-': case 'N': case 'X': case '?': break;
                    default: return false;
                }
            }
        }
        return true;
    }

    //~~~~~~~~~//

    bool DataHolder::valid_phyml_aa() const {
        for (unsigned int j=0; j<_ls; j++) {
            for (unsigned int i=0; i<_nsi; i++) {
                switch (_datai[i][j]) {
                    case 'A': case 'R': case 'N': case 'B':
                    case 'D': case 'C': case 'Q': case 'Z':
                    case 'E': case 'G': case 'H': case 'I':
                    case 'L': case 'K': case 'M': case 'F':
                    case 'P': case 'S': case 'T': case 'W':
                    case 'Y': case 'V': case '-': case 'X': case '?': break;
                    default: return false;
                }
            }
            for (unsigned int i=0; i<_nso; i++) {
                switch (_datai[i][j]) {
                    case 'A': case 'R': case 'N': case 'B':
                    case 'D': case 'C': case 'Q': case 'Z':
                    case 'E': case 'G': case 'H': case 'I':
                    case 'L': case 'K': case 'M': case 'F':
                    case 'P': case 'S': case 'T': case 'W':
                    case 'Y': case 'V': case '-': case 'X': case '?': break;
                    default: return false;
                }
            }
        }
        return true;
    }

    //~~~~~~~~~//

    bool DataHolder::valid_phyml_names() const {
        for (unsigned int i=0; i<_nsi; i++) {
            if (_lni[i] == 1) return false;
            for (unsigned int j=0; j<_lni[i]; j++) {
                switch (_namesi[i][j]) {
                    case ' ': case '\t': case '\r': case '\n':
                    case ',': case ':': case '(': case ')':
                        return false;
                }
            }
        }
        for (unsigned int i=0; i<_nso; i++) {
            if (_lno[i] == 1) return false;
            for (unsigned int j=0; j<_lno[i]; j++) {
                switch (_nameso[i][j]) {
                    case ' ': case '\t': case '\r': case '\n':
                    case ',': case ':': case '(': case ')':
                        return false;
                }
            }
        }
        return true;
    }

    //~~~~~~~~~//

    IntersperseAlign::IntersperseAlign() {
        _alleles = (int *) malloc(1 * sizeof(int));
        if (!_alleles) throw EGGMEM;
        _alleles[0] = static_cast<int>('A');
        _num_alleles = 1;
        _res_alleles = 1;

        _positions = NULL;
        _round_positions = NULL;
        _offset = NULL;
        _res_positions = 0;

        _data = NULL;
        _nsites = 0;
        _length = 0;
        _random = NULL;
    }

    //~~~~~~~~~//

    IntersperseAlign::~IntersperseAlign() {
        if (_alleles) free(_alleles); // in principle, test is not needed
        if (_positions) free(_positions);
        if (_round_positions) free(_round_positions);
        if (_offset) free(_offset);
    }

    //~~~~~~~~~//

    void IntersperseAlign::load(DataHolder& data) {
        _data = & data;
        _nsites = data.get_nsit();
        if (_nsites > _res_positions) {
            _positions = (double *) realloc(_positions, _nsites * sizeof(double));
            if (!_positions) throw EGGMEM;
            _round_positions = (unsigned int *) realloc(_round_positions, _nsites * sizeof(unsigned int));
            if (!_round_positions) throw EGGMEM;
            _offset = (unsigned int *) realloc(_offset, (_nsites + 1) * sizeof(unsigned int));
            if (!_offset) throw EGGMEM;
            _res_positions = _nsites;
        }
    }

    //~~~~~~~~~//

    void IntersperseAlign::set_length(unsigned int length) {
        _length = length;
    }

    //~~~~~~~~~//

    void IntersperseAlign::set_position(unsigned int index, double position) {
        _positions[index] = position;
    }

    //~~~~~~~~~//

    void IntersperseAlign::set_round_position(unsigned int index, unsigned int position) {
        _round_positions[index] = position;
    }

    //~~~~~~~~~//

    void IntersperseAlign::get_positions(const Coalesce& coalesce) {
        for (unsigned int i=0; i<_nsites; i++) {
            _positions[i] = coalesce.site_position(i); // positions are guaranteed to be increasing
        }
    }

    //~~~~~~~~~//

    void IntersperseAlign::set_num_alleles(unsigned int num) {
        if (num > _res_alleles) {
            _alleles = (int *) realloc(_alleles, num * sizeof(int));
            if (!_alleles) EGGMEM;
            _res_alleles = num;
        }
        _num_alleles = num;
    }

    //~~~~~~~~~//

    void IntersperseAlign::set_allele(unsigned int index, int allele) {
        _alleles[index] = allele;
    }

    //~~~~~~~~~//

    void IntersperseAlign::set_random(Random * random) {
        _random = random;
    }

    //~~~~~~~~~//

    void IntersperseAlign::intersperse(bool round_positions) {

        // don't do anything if no sites need to be inserted

        if (_nsites >= _length) return;

        // convert positions to integers

        if (round_positions) {
            for (unsigned int i=0; i<_nsites; i++) {
                _round_positions[i] = floor(0.5 + _positions[i] * (_length-1));
            }
        }

        // determines offsets

        unsigned int tot_insert = 0;

        // first site (if any)

        if (_nsites > 0) {
            _offset[0] = _round_positions[0];
            tot_insert += _offset[0];
        }

        // other sites

        for (unsigned int i=1; i<_nsites; i++) {
            _offset[i] =  (_round_positions[i] == _round_positions[i-1]) ? 0 : _round_positions[i] - _round_positions[i-1] - 1;
            tot_insert += _offset[i];
        }

        // after last site

        if (_nsites > 0) _offset[_nsites] = _length - _round_positions[_nsites-1] - 1;
        else _offset[_nsites] = _length;
        tot_insert += _offset[_nsites];

        // correct error randomly (not enough or too many inserted sites)

        int error = tot_insert + _nsites - _length;

        double X;

        while (error != 0) {

            // pick a random interval

            X = _random->uniform() * tot_insert;

            unsigned int i;
            unsigned int acc = 0;
            for (i=0; i<_nsites+1; i++) {
                acc += _offset[i];
                if (X < acc) break;
            }

            // remove/add a site there

            if (error > 0) {
                _offset[i]--;
                tot_insert--;
                error--;
            }

            else {
               _offset[i]++;
               tot_insert++; 
               error++;
           }
        }

        // insert sites in the alignment

        unsigned int pos = 0;
        int allele;
        for (unsigned int site=0; site<_nsites+1; site++) {
            _data->insert_sites(pos, _offset[site]);
            for (unsigned int y=pos; y<pos+_offset[site]; y++) {
                if (_num_alleles == 1) allele = _alleles[0];
                else allele = _alleles[_random->irand(_num_alleles)];
                for (unsigned int x=0; x<_data->get_nsam_i(); x++)  _data->set_i(x, y, allele);
                for (unsigned int x=0; x<_data->get_nsam_o(); x++)  _data->set_o(x, y, allele);
            }
            pos += _offset[site] + 1;
        }
    }
}
