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

#include <cstring>
#include <cmath>
#include <string>
#include <iostream>
#include <cctype>
#include <new>
#include "egglib.hpp"
#include "VCF.hpp"
#include "Index.hpp"
using namespace std;

namespace egglib {

    // seems they must be defined here because however NUM_ALTERNATE is... equal to MAX (meaning MAX is 0)
    const double UNDEF = - std::numeric_limits<double>::max();
    const int MISSINGDATA = std::numeric_limits<int>::max();
    const char MAXCHAR = std::numeric_limits<char>::max();
    const unsigned int UNKNOWN = std::numeric_limits<unsigned int>::max();
    const unsigned int MAX = UNKNOWN;
    const unsigned int OUTGROUP = UNKNOWN;
    const unsigned int MISSING = UNKNOWN;
    const unsigned long int BEFORE = std::numeric_limits<unsigned long int>::max() - 2; // don't use max() in case unsigned long int == unsigned int

    namespace vcf {

        const unsigned int NUM_ALTERNATE = MAX - 1;
        const unsigned int NUM_GENOTYPES = MAX - 2;
        const unsigned int NUM_POSSIBLE_ALLELES = MAX - 3;

        //////

        const char * _index_from_fname (const char * filename){
            string file(filename);
            string n_file;
            char * r_file;
            string path;
            string f_file;
            size_t none=-1;
            size_t slash = file.find_last_of("/\\");
            if(slash != none){ 
                path = file.substr(0, slash);
                file = file.substr(slash+1);
            }

            if(file.substr(file.find_last_of(".") + 1) == "vcf"){
                std::size_t ext = file.find_last_of(".");
                file.replace(ext, sizeof(".vcfi")-1 ,".vcfi");
            } else {
                file=file+".vcfi";
            }

            if (slash != none) {
                f_file=path+"/"+file;
            } else {
                f_file=file;
            }
            r_file=(char *) malloc(f_file.length()+1 * sizeof(char));
            strcpy(r_file, f_file.c_str());
            return (const char *)r_file;
        }

        //////

        Filter::Filter() {
            init();
        }

        //////

        Filter::Filter(const char * id, const char * descr) {
            init();
            set_ID(id);
            set_description(descr);
        }

        //////

        void Filter::update(const char * id, const char * descr) {
            set_ID(id);
            set_description(descr);
            _extra_n = 0;
        }

        //////

        Filter::Filter(const Filter& src) {
            copy(src);
        }

        //////

        Filter& Filter::operator=(const Filter& src) {
            clear();
            copy(src);
            return * this;
        }

        //////

        void Filter::copy(const Filter& src) {
            init();
            set_ID(src._ID);
            set_description(src._descr);
            for (unsigned int i=0; i<src._extra_n; i++) set_extra(src._extra_key[i], src._extra_val[i]);
        }

        //////

        Filter::~Filter() {
            free();
        }

        //////

        const char * Filter::get_ID() const {
            return _ID;
        }

        //////

        void Filter::set_ID(const char * id) {
            unsigned int n = strlen(id) + 1;
            if (n > _ID_r) {
                _ID = (char *) realloc(_ID, n * sizeof(char));
                if (!_ID) throw EGGMEM;
                _ID_r = n;
            }
            strcpy(_ID, id);
        }

        //////

        const char * Filter::get_description() const {
            return _descr;
        }

        //////

        void Filter::set_description(const char * descr) {
            unsigned int n = strlen(descr) + 1;
            if (n > _descr_r) {
                _descr = (char *) realloc(_descr, n * sizeof(char));
                if (!_descr) throw EGGMEM;
                _descr_r = n;
            }
            strcpy(_descr, descr);
        }

        //////

        void Filter::clear() {
            free();
            init();
        }

        //////

        void Filter::init() {
            _ID_r = 0;
            _descr_r = 0;
            _ID = NULL;
            _descr = NULL;
            _extra_n = 0;
            _extra_r = 0;
            _extra_key_r = 0;
            _extra_val_r = 0;
            _extra_key = NULL;
            _extra_val = NULL;
        }

        //////

        void Filter::free() {
            if (_ID) ::free(_ID);
            if (_descr) ::free(_descr);
            for (unsigned int i=0; i<_extra_r; i++) {
                if (_extra_key[i]) ::free(_extra_key[i]);
                if (_extra_val[i]) ::free(_extra_val[i]);
            }
            if (_extra_key) ::free(_extra_key);
            if (_extra_val) ::free(_extra_val);
            if (_extra_key_r) ::free(_extra_key_r);
            if (_extra_val_r) ::free(_extra_val_r);
        }

        //////

        void Filter::set_extra(const char * key, const char * value) {
            _extra_n++;
            if (_extra_n > _extra_r) {
                _extra_key_r = (unsigned int *) realloc(_extra_key_r, _extra_n * sizeof(unsigned int));
                if (!_extra_key_r) throw EGGMEM;

                _extra_val_r = (unsigned int *) realloc(_extra_val_r, _extra_n * sizeof(unsigned int));
                if (!_extra_val_r) throw EGGMEM;

                _extra_key = (char **) realloc(_extra_key, _extra_n * sizeof(char *));
                if (!_extra_key) throw EGGMEM;

                _extra_val = (char **) realloc(_extra_val, _extra_n * sizeof(char *));
                if (!_extra_val) throw EGGMEM;

                _extra_key_r[_extra_n-1] = 0;
                _extra_val_r[_extra_n-1] = 0;
                _extra_key[_extra_n-1] = NULL;
                _extra_val[_extra_n-1] = NULL;

                _extra_r = _extra_n;
            }

            unsigned int n = strlen(key) + 1;
            if (n > _extra_key_r[_extra_n-1]) {
                _extra_key[_extra_n-1] = (char *) realloc(_extra_key[_extra_n-1], n * sizeof(char));
                if (!_extra_key[_extra_n-1]) throw EGGMEM;
                _extra_key_r[_extra_n-1] = n;
            }
            strcpy(_extra_key[_extra_n-1], key);

            n = strlen(value) + 1;
            if (n > _extra_val_r[_extra_n-1]) {
                _extra_val[_extra_n-1] = (char *) realloc(_extra_val[_extra_n-1], n * sizeof(char));
                if (!_extra_val[_extra_n-1]) throw EGGMEM;
                _extra_val_r[_extra_n-1] = n;
            }
            strcpy(_extra_val[_extra_n-1], value);
        }

        //////

        const char * Filter::get_extra_key(unsigned int idx) const {
            return _extra_key[idx];
        }

        //////

        const char * Filter::get_extra_value(unsigned int idx) const {
            return _extra_val[idx];
        }

        //////

        unsigned int Filter::get_num_extra() const {
            return _extra_n;
        }

        ////////////////////////////////////////////////////////////////////

        Info::Info(const char * id, unsigned int num, Info::Type t, const char * descr) {
            init();
            Filter::update(id, descr);
            set_type(t);
            set_number(num);
        }

        //////

        void Info::update(const char * id, unsigned int num, Info::Type t, const char * descr) {
            Filter::update(id, descr);
            set_type(t);
            set_number(num);
        }

        //////

        Info::Info(const Info& src) {
            init();
            Filter::copy(src);
            set_type(src._type);
            set_number(src._number);
        }

        //////

        Info& Info::operator=(const Info& src) {
            clear();
            Filter::copy(src);
            set_type(src._type);
            set_number(src._number);
            return *this;
        }

        //////

        unsigned int Info::get_number() const {
            return _number;
        }

        //////

        void Info::set_number(unsigned int num) {
            _number = num;
        }

        //////

        Info::Type Info::get_type() const {
            return _type;
        }

        //////

        void Info::set_type(Info::Type t) {
            _type = t;
        }

        ////////////////////////////////////////////////////////////////////

        Format::Format(const char * id, unsigned int num, Info::Type t, const char * descr) {
            init();
            Filter::update(id, descr);
            set_number(num);
            set_type(t);
        }

        //////

        void Format::update(const char * id, unsigned int num, Info::Type t, const char * descr) {
            Filter::update(id, descr);
            set_number(num);
            set_type(t);
        }

        //////

        Format::Format(const Format& src) {
            init();
            copy(src);
            set_type(src._type);
            set_number(src._number);
        }

        //////

        Format::Format(const Info& src) {
            init();
            Filter::update(src.get_ID(), src.get_description());
            set_type(src.get_type());
            set_number(src.get_number());
        }

        //////

        Format& Format::operator=(const Format& src) {
            clear();
            copy(src);
            set_type(src._type);
            set_number(src._number);
            return *this;
        }

        //////

        Format& Format::operator=(const Info& src) {
            clear();
            set_ID(src.get_ID());
            Filter::update(src.get_ID(), src.get_description());
            set_type(src.get_type());
            set_number(src.get_number());
            return *this;
        }

        //////

        void Format::set_type(Info::Type t) {
            if (t == Info::Flag) throw EggArgumentValueError("flag is not allowed as FORMAT type");
            _type = t;
        }

        //////

        Info::Type Format::get_type() const {
            return _type;
        }

        ////////////////////////////////////////////////////////////////////

        Meta::Meta() {
            init();
        }

        //////

        Meta::Meta(const char * k, const char * v) {
            init();
            set_key(k);
            set_value(v);
        }

        //////

        void Meta::update(const char * k, const char * v) {
            set_key(k);
            set_value(v);
        }

        //////

        Meta::Meta(Meta& src) {
            init();
            set_key(src._key);
            set_value(src._val);
        }

        //////

        Meta& Meta::operator=(Meta& src) {
            clear();
            set_key(src._key);
            set_value(src._val);
            return *this;
        }

        //////

        Meta::~Meta() {
            free();
        }

        //////

        void Meta::set_key(const char * k) {
            unsigned int n = strlen(k) + 1;
            if (n > _key_r) {
                _key = (char *) realloc(_key, n * sizeof(char));
                if (!_key) throw EGGMEM;
                _key_r = n;
            }
            strcpy(_key, k);
        }
            
        //////

        void Meta::set_value(const char * v) {
            unsigned int n = strlen(v) + 1;
            if (n > _val_r) {
                _val = (char *) realloc(_val, n * sizeof(char));
                if (!_val) throw EGGMEM;
                _val_r = n;
            }
            strcpy(_val, v);
        }

        //////

        const char * Meta::get_key() const {
            return _key;
        }

        //////

        const char * Meta::get_value() const {
            return _val;
        }

        //////

        void Meta::clear() {
            free();
            init();
        }

        //////

        void Meta::init() {
            _key_r = 0;
            _key = NULL;
            _val_r = 0;
            _val = NULL;
        }

        //////

        void Meta::free() {
            if (_key) ::free(_key);
            if (_val) ::free(_val);
        }

        ////////////////////////////////////////////////////////////////////

        Alt& Alt::operator=(Alt& src) {
            Filter::operator=(src);
            return *this;
        }

        //////

        Alt& Alt::operator=(Filter& src) {
            Filter::operator=(src);
            return *this;
        }

        ////////////////////////////////////////////////////////////////////

        FlagInfo::FlagInfo() {
            _res_ID = 0;
            _ID = NULL;
        }

        //////

        FlagInfo::FlagInfo(const FlagInfo& src) {
            _res_ID = 0;
            _ID = NULL;
            copy(src);
        }

        //////

        FlagInfo& FlagInfo::operator=(const FlagInfo& src) {
            copy(src);
            return *this;
        }

        //////

        FlagInfo::~FlagInfo() {
            if (_ID) free(_ID);
        }

        //////

        void FlagInfo::set_ID(const char * id) {
            unsigned int n = strlen(id) + 1;
            if (n > _res_ID) {
                _ID = (char *) realloc(_ID, n * sizeof(char));
                if (!_ID) throw EGGMEM;
                _res_ID = n;
            }
            strcpy(_ID, id);
        }

        //////

        const char * FlagInfo::get_ID() const {
            return _ID;
        }

        //////

        void FlagInfo::copy(const FlagInfo& src) {
            unsigned int n = strlen(src._ID) + 1;
            if (n > _res_ID) {
                _ID = (char *) realloc(_ID, n * sizeof(char));
                if (!_ID) throw EGGMEM;
                _res_ID = n;
            }
            strcpy(_ID, src._ID);
           
        }
        
        ////////////////////////////////////////////////////////////////////

        StringInfo::StringInfo() {
            _res_ID = 0;
            _ID = NULL;
            _num_items = 0;
            _res_items = 0;
            _res_len_items = NULL;
            _items = NULL;
        }

        //////

        StringInfo::StringInfo(const StringInfo& src) {
            _res_ID = 0;
            _ID = NULL;
            _num_items = 0;
            _res_items = 0;
            _res_len_items = NULL;
            _items = NULL;
            copy(src);
        }

        //////

        StringInfo& StringInfo::operator=(const StringInfo& src) {
            copy(src);
            return * this;
        }
        
        //////
        
        StringInfo::~StringInfo() {
            for (unsigned int i=0; i<_res_items; i++) if (_items[i]) free(_items[i]);
            if (_res_len_items) free(_res_len_items);
        }
        
       //////
        
        void StringInfo::add() {
            _num_items++;
            if (_num_items > _res_items) {
                _items = (char **) realloc(_items, _num_items * sizeof(char*));
                if (!_items) throw EGGMEM;
                _items[_num_items-1] = (char *) malloc(1 * sizeof(char));
                if (!_items[_num_items-1]) throw EGGMEM;
                _items[_num_items-1][0] = '\0';
                _res_len_items = (unsigned int *) realloc(_res_len_items, _num_items * sizeof(unsigned int));
                if (!_res_len_items) throw EGGMEM;
                _res_items = _num_items;
            }
            _res_len_items[_num_items-1] = 1;
        }

        //////
        
        void StringInfo::copy(const StringInfo& src) {
            set_ID(src._ID);
            _num_items = src._num_items;
            if (_num_items > _res_items) {
                _items = (char **) realloc(_items, _num_items * sizeof(char *));
                if (!_items) throw EGGMEM;
                _res_len_items = (unsigned int *) realloc(_res_len_items, _num_items * sizeof(unsigned int));
                if (!_res_len_items) throw EGGMEM;
                for (unsigned int i=_res_items; i<_num_items; i++) {
                    _items[i] = NULL;
                    _res_len_items[i] = 0;
                }
                _res_items = _num_items;
            }
            for (unsigned int i=0; i<_num_items; i++) {
                unsigned int n = strlen(src._items[i]) + 1;
                if (n > _res_len_items[i]) {
                    _items[i] = (char *) realloc(_items[i], n * sizeof(char));
                    if (!_items[i]) throw EGGMEM;
                    _res_len_items[i] = n;
                }
                strcpy(_items[i], src._items[i]);
            }
        }
        
        //////

        void StringInfo::change(unsigned int item, unsigned int position, char value) {
            _items[item][position] = value;
        }

        ////////////////////////////////////////////////////////////////////

        SampleInfo::SampleInfo() {
            init();
        }

        //////

        SampleInfo::SampleInfo(const SampleInfo& src) {
            init();
            copy(src);
        }

        //////

        SampleInfo& SampleInfo::operator=(const SampleInfo& src) {
            copy(src);
            return *this;
        }

        //////

        SampleInfo::~SampleInfo() {
            free();
        }
            
        //////

        void SampleInfo::reset() {
            _num_IntegerEntries = 0;
            _num_FloatEntries = 0;
            _num_CharacterEntries = 0;
            _num_StringEntries = 0;
        }

        //////
            
        void SampleInfo::clear() {
            free();
            init();
        }
            
        //////

        unsigned int SampleInfo::num_IntegerEntries() const {
            return _num_IntegerEntries;
        }

        //////
            
        unsigned int SampleInfo::num_IntegerItems(unsigned int i) const {
            return _num_IntegerItems[i];
        }

        //////
            
        int SampleInfo::IntegerItem(unsigned int i, unsigned int j) const {
            return _IntegerItems[i][j];
        }

        //////
            
        unsigned int SampleInfo::num_FloatEntries() const {
            return _num_FloatEntries;
        }

        //////
            
        unsigned int SampleInfo::num_FloatItems(unsigned int i) const {
            return _num_FloatItems[i];
        }

        //////
            
        double SampleInfo::FloatItem(unsigned int i, unsigned int j) const {
            return _FloatItems[i][j];
        }

        //////
            
        unsigned int SampleInfo::num_CharacterEntries() const {
            return _num_CharacterEntries;
        }

        //////
            
        unsigned int SampleInfo::num_CharacterItems(unsigned int i) const {
            return _num_CharacterItems[i];
        }

        //////
            
        char SampleInfo::CharacterItem(unsigned int i, unsigned int j) const {
            return _CharacterItems[i][j];
        }

        //////
            
        unsigned int SampleInfo::num_StringEntries() const {
            return _num_StringEntries;
        }

        //////
            
        unsigned int SampleInfo::num_StringItems(unsigned int i) const {
            return _num_StringItems[i];
        }

        //////
            
        const char * SampleInfo::StringItem(unsigned int i, unsigned int j) const {
            return _StringItems[i][j];
        }

        //////

        void SampleInfo::init() {
            _num_IntegerEntries = 0;
            _res_IntegerEntries = 0;
            _num_IntegerItems = NULL;
            _res_IntegerItems = NULL;
            _IntegerItems = NULL;
            
            _num_FloatEntries = 0;
            _res_FloatEntries = 0;
            _num_FloatItems = NULL;
            _res_FloatItems = NULL;
            _FloatItems = NULL;
            
            _num_CharacterEntries = 0;
            _res_CharacterEntries = 0;
            _num_CharacterItems = NULL;
             _res_CharacterItems = NULL;
            _CharacterItems = NULL;
            
            _num_StringEntries = 0;
            _res_StringEntries = 0;
            _num_StringItems = NULL;
             _res_StringItems = NULL;
            _res_len_StringItems = NULL;
            _StringItems = NULL;
        }
        
        //////
        
        void SampleInfo::copy(const SampleInfo& src) {
        
            // copy Integer data
        
            _num_IntegerEntries = src._num_IntegerEntries;
            if (_num_IntegerEntries > _res_IntegerEntries) {
                _num_IntegerItems = (unsigned int *) realloc(_num_IntegerItems, _num_IntegerEntries * sizeof(unsigned int));
                if (!_num_IntegerItems) throw EGGMEM;
                _res_IntegerItems = (unsigned int *) realloc(_res_IntegerItems, _num_IntegerEntries * sizeof(unsigned int));
                if (!_res_IntegerItems) throw EGGMEM;
                _IntegerItems = (int **) realloc(_IntegerItems, _num_IntegerEntries * sizeof(int *));
                if (!_IntegerItems) throw EGGMEM;
                for (unsigned int i=_res_IntegerEntries; i<_num_IntegerEntries; i++) {
                    _num_IntegerItems[i] = 0;
                    _res_IntegerItems[i] = 0;
                    _IntegerItems[i] = NULL;
                }
                _res_IntegerEntries = _num_IntegerEntries;
            }
            for (unsigned int i=0; i<_num_IntegerEntries; i++) {
                _num_IntegerItems[i] = src._num_IntegerItems[i];
                if (_num_IntegerItems[i] > _res_IntegerItems[i]) {
                    _IntegerItems[i] = (int *) realloc(_IntegerItems[i], _num_IntegerItems[i] * sizeof(int));
                    if (!_IntegerItems[i]) throw EGGMEM;
                    _res_IntegerItems[i] = _num_IntegerItems[i];
                }
                for (unsigned int j=0; j<_num_IntegerItems[i]; j++) _IntegerItems[i][j] = src._IntegerItems[i][j];
            }

            // copy Float data
        
            _num_FloatEntries = src._num_FloatEntries;
            if (_num_FloatEntries > _res_FloatEntries) {
                _num_FloatItems = (unsigned int *) realloc(_num_FloatItems, _num_FloatEntries * sizeof(unsigned int));
                if (!_num_FloatItems) throw EGGMEM;
                _res_FloatItems = (unsigned int *) realloc(_res_FloatItems, _num_FloatEntries * sizeof(unsigned int));
                if (!_res_FloatItems) throw EGGMEM;
                _FloatItems = (double **) realloc(_FloatItems, _num_FloatEntries * sizeof(double *));
                if (!_FloatItems) throw EGGMEM;
                for (unsigned int i=_res_FloatEntries; i<_num_FloatEntries; i++) {
                    _num_FloatItems[i] = 0;
                    _res_FloatItems[i] = 0;
                    _FloatItems[i] = NULL;
                }
                _res_FloatEntries = _num_FloatEntries;
            }
            for (unsigned int i=0; i<_num_FloatEntries; i++) {
                _num_FloatItems[i] = src._num_FloatItems[i];
                if (_num_FloatItems[i] > _res_FloatItems[i]) {
                    _FloatItems[i] = (double *) realloc(_FloatItems[i], _num_FloatItems[i] * sizeof(double));
                    if (!_FloatItems[i]) throw EGGMEM;
                    _res_FloatItems[i] = _num_FloatItems[i];
                }
                for (unsigned int j=0; j<_num_FloatItems[i]; j++) _FloatItems[i][j] = src._FloatItems[i][j];
            }

            // copy Char data
        
            _num_CharacterEntries = src._num_CharacterEntries;
            if (_num_CharacterEntries > _res_CharacterEntries) {
                _num_CharacterItems = (unsigned int *) realloc(_num_CharacterItems, _num_CharacterEntries * sizeof(unsigned int));
                if (!_num_CharacterItems) throw EGGMEM;
                _res_CharacterItems = (unsigned int *) realloc(_res_CharacterItems, _num_CharacterEntries * sizeof(unsigned int));
                if (!_res_CharacterItems) throw EGGMEM;
                _CharacterItems = (char **) realloc(_CharacterItems, _num_CharacterEntries * sizeof(char *));
                if (!_CharacterItems) throw EGGMEM;
                for (unsigned int i=_res_CharacterEntries; i<_num_CharacterEntries; i++) {
                    _num_CharacterItems[i] = 0;
                    _res_CharacterItems[i] = 0;
                    _CharacterItems[i] = NULL;
                }
                _res_CharacterEntries = _num_CharacterEntries;
            }
            for (unsigned int i=0; i<_num_CharacterEntries; i++) {
                _num_CharacterItems[i] = src._num_CharacterItems[i];
                if (_num_CharacterItems[i] > _res_CharacterItems[i]) {
                    _CharacterItems[i] = (char *) realloc(_CharacterItems[i], _num_CharacterItems[i] * sizeof(char));
                    if (!_CharacterItems[i]) throw EGGMEM;
                    _res_CharacterItems[i] = _num_CharacterItems[i];
                }
                for (unsigned int j=0; j<_num_CharacterItems[i]; j++) _CharacterItems[i][j] = src._CharacterItems[i][j];
            }

            // copy String data
        
            _num_StringEntries = src._num_StringEntries;
            if (_num_StringEntries > _res_StringEntries) {
                _num_StringItems = (unsigned int *) realloc(_num_StringItems, _num_StringEntries * sizeof(unsigned int));
                if (!_num_StringItems) throw EGGMEM;
                _res_StringItems = (unsigned int *) realloc(_res_StringItems, _num_StringEntries * sizeof(unsigned int));
                if (!_res_StringItems) throw EGGMEM;
                _res_len_StringItems = (unsigned int **) realloc(_res_len_StringItems, _num_StringEntries * sizeof(unsigned int *));
                if (!_res_StringItems) throw EGGMEM;
                _StringItems = (char ***) realloc(_StringItems, _num_StringEntries * sizeof(char **));
                if (!_StringItems) throw EGGMEM;
                for (unsigned int i=_res_StringEntries; i<_num_StringEntries; i++) {
                    _num_StringItems[i] = 0;
                    _res_StringItems[i] = 0;
                    _res_len_StringItems[i] = NULL;
                    _StringItems[i] = NULL;
                }
                _res_StringEntries = _num_StringEntries;
            }
            for (unsigned int i=0; i<_num_StringEntries; i++) {
                _num_StringItems[i] = src._num_StringItems[i];
                if (_num_StringItems[i] > _res_StringItems[i]) {
                    _StringItems[i] = (char **) realloc(_StringItems[i], _num_StringItems[i] * sizeof(char *));
                    if (!_StringItems[i]) throw EGGMEM;
                    _res_len_StringItems[i] = (unsigned int *) realloc(_res_len_StringItems[i], _num_StringItems[i] * sizeof(unsigned int *));
                    if (!_res_len_StringItems[i]) throw EGGMEM;
                    for (unsigned int j=_res_StringItems[i]; j<_num_StringItems[i]; j++) {
                        _res_len_StringItems[i][j] = 1;
                        _StringItems[i][j] = (char *) malloc(1 * sizeof(char));
                        if (!_StringItems[i][j]) throw EGGMEM;
                        _StringItems[i][j][0] = '\0';
                    }
                    _res_StringItems[i] = _num_StringItems[i];
                }
                for (unsigned int j=0; j<_num_StringItems[i]; j++) {
                    unsigned int n = strlen(src._StringItems[i][j]) + 1;
                    if (n > _res_len_StringItems[i][j]) {
                        _StringItems[i][j] = (char *) realloc(_StringItems[i][j], n * sizeof(char));
                        if (!_StringItems[i][j]) throw EGGMEM;
                        _res_len_StringItems[i][j] = n;
                    }
                    strcpy(_StringItems[i][j], src._StringItems[i][j]);
                }
            }
        }
        
        //////
        
        void SampleInfo::free() {

            for (unsigned int i=0; i<_res_IntegerEntries; i++) {
                if (_IntegerItems[i]) ::free(_IntegerItems[i]);
            }
            if (_IntegerItems) ::free(_IntegerItems);
            if (_num_IntegerItems) ::free(_num_IntegerItems);
            if (_res_IntegerItems) ::free(_res_IntegerItems);

            for (unsigned int i=0; i<_res_FloatEntries; i++) {
                if (_FloatItems[i]) ::free(_FloatItems[i]);
            }
            if (_FloatItems) ::free(_FloatItems);
            if (_num_FloatItems) ::free(_num_FloatItems);
            if (_res_FloatItems) ::free(_res_FloatItems);

            for (unsigned int i=0; i<_res_CharacterEntries; i++) {
                if (_CharacterItems[i]) ::free(_CharacterItems[i]);
            }
            if (_CharacterItems) ::free(_CharacterItems);
            if (_num_CharacterItems) ::free(_num_CharacterItems);
            if (_res_CharacterItems) ::free(_res_CharacterItems);

            for (unsigned int i=0; i<_res_StringEntries; i++) {
                if (_res_len_StringItems[i]) ::free(_res_len_StringItems[i]);
                for (unsigned int j=0; j<_res_StringItems[i]; j++) {
                    if (_StringItems[i][j]) ::free(_StringItems[i][j]);
                }
                if (_StringItems[i]) ::free(_StringItems[i]);
            }
            if (_num_StringItems) ::free(_num_StringItems);
            if (_res_StringItems) ::free(_res_StringItems);
            if (_res_len_StringItems) ::free(_res_len_StringItems);
            if (_StringItems) ::free(_StringItems);
        }

        //////
            
        void SampleInfo::addIntegerEntry() {
            _num_IntegerEntries++;
            if (_num_IntegerEntries > _res_IntegerEntries) {
                _num_IntegerItems = (unsigned int *) realloc(_num_IntegerItems, _num_IntegerEntries * sizeof(unsigned int));
                if (!_num_IntegerItems) throw EGGMEM;
                _res_IntegerItems = (unsigned int *) realloc(_res_IntegerItems, _num_IntegerEntries * sizeof(unsigned int));
                if (!_res_IntegerItems) throw EGGMEM;
                _res_IntegerItems[_num_IntegerEntries-1] = 0;
                _IntegerItems = (int **) realloc(_IntegerItems, _num_IntegerEntries * sizeof(int*));
                if (!_IntegerItems) throw EGGMEM;
                _IntegerItems[_num_IntegerEntries-1] = NULL;
                _res_IntegerEntries = _num_IntegerEntries;
            }
            _num_IntegerItems[_num_IntegerEntries-1] = 0;
        }

        //////
            
        void SampleInfo::addIntegerItem() {
            _num_IntegerItems[_num_IntegerEntries-1]++;
            if(_num_IntegerItems[_num_IntegerEntries-1] > _res_IntegerItems[_num_IntegerEntries-1]) {
                _IntegerItems[_num_IntegerEntries-1] = (int *) realloc(_IntegerItems[_num_IntegerEntries-1], _num_IntegerItems[_num_IntegerEntries-1] * sizeof(int));
                if (!_IntegerItems[_num_IntegerEntries-1]) throw EGGMEM;
                _res_IntegerItems[_num_IntegerEntries-1] = _num_IntegerItems[_num_IntegerEntries-1];
            }
        }

        //////
            
        void SampleInfo::addFloatEntry() {
            _num_FloatEntries++;
            if (_num_FloatEntries > _res_FloatEntries) {
                _num_FloatItems = (unsigned int *) realloc(_num_FloatItems, _num_FloatEntries * sizeof(unsigned int));
                if (!_num_FloatItems) throw EGGMEM;
                _res_FloatItems = (unsigned int *) realloc(_res_FloatItems, _num_FloatEntries * sizeof(unsigned int));
                if (!_res_FloatItems) throw EGGMEM;
                _res_FloatItems[_num_FloatEntries-1] = 0;
                _FloatItems = (double **) realloc(_FloatItems, _num_FloatEntries * sizeof(double*));
                if (!_FloatItems) throw EGGMEM;
                _FloatItems[_num_FloatEntries-1] = NULL;
                _res_FloatEntries = _num_FloatEntries;
            }
            _num_FloatItems[_num_FloatEntries-1] = 0;
        }

        //////
            
        void SampleInfo::addFloatItem() {
            _num_FloatItems[_num_FloatEntries-1]++;
            if(_num_FloatItems[_num_FloatEntries-1] > _res_FloatItems[_num_FloatEntries-1]) {
                _FloatItems[_num_FloatEntries-1] = (double *) realloc(_FloatItems[_num_FloatEntries-1], _num_FloatItems[_num_FloatEntries-1] * sizeof(double));
                if (!_FloatItems[_num_FloatEntries-1]) throw EGGMEM;
                _res_FloatItems[_num_FloatEntries-1] = _num_FloatItems[_num_FloatEntries-1];
            }
        }

        //////

        void SampleInfo::addCharacterEntry() {
            _num_CharacterEntries++;
            if (_num_CharacterEntries > _res_CharacterEntries) {
                _num_CharacterItems = (unsigned int *) realloc(_num_CharacterItems, _num_CharacterEntries * sizeof(unsigned int));
                if (!_num_CharacterItems) throw EGGMEM;
                _res_CharacterItems = (unsigned int *) realloc(_res_CharacterItems, _num_CharacterEntries * sizeof(unsigned int));
                if (!_res_CharacterItems) throw EGGMEM;
                _res_CharacterItems[_num_CharacterEntries-1] = 0;
                _CharacterItems = (char **) realloc(_CharacterItems, _num_CharacterEntries * sizeof(char*));
                if (!_CharacterItems) throw EGGMEM;
                _CharacterItems[_num_CharacterEntries-1] = NULL;
                _res_CharacterEntries = _num_CharacterEntries;
            }
            _num_CharacterItems[_num_CharacterEntries-1] = 0;
        }

        //////
            
        void SampleInfo::addCharacterItem() {
            _num_CharacterItems[_num_CharacterEntries-1]++;
            if(_num_CharacterItems[_num_CharacterEntries-1] > _res_CharacterItems[_num_CharacterEntries-1]) {
                _CharacterItems[_num_CharacterEntries-1] = (char *) realloc(_CharacterItems[_num_CharacterEntries-1], _num_CharacterItems[_num_CharacterEntries-1] * sizeof(char));
                if (!_CharacterItems[_num_CharacterEntries-1]) throw EGGMEM;
                _res_CharacterItems[_num_CharacterEntries-1] = _num_CharacterItems[_num_CharacterEntries-1];
            }
        }

        //////
            
        void SampleInfo::addStringEntry() {
            _num_StringEntries++;
            if (_num_StringEntries > _res_StringEntries) {
                _num_StringItems = (unsigned int *) realloc(_num_StringItems, _num_StringEntries * sizeof(unsigned int));
                if (!_num_StringItems) throw EGGMEM;
                _res_StringItems = (unsigned int *) realloc(_res_StringItems, _num_StringEntries * sizeof(unsigned int));
                if (!_res_StringItems) throw EGGMEM;
                _res_StringItems[_num_StringEntries-1] = 0;
                _res_len_StringItems = (unsigned int **) realloc(_res_len_StringItems, _num_StringEntries * sizeof(unsigned int*));
                if (!_res_len_StringItems) throw EGGMEM;
                _res_len_StringItems[_num_StringEntries-1] = NULL;
                _StringItems = (char ***) realloc(_StringItems, _num_StringEntries * sizeof(char **));
                if (!_StringItems) throw EGGMEM;
                _StringItems[_num_StringEntries-1] = NULL;            
                _res_StringEntries = _num_StringEntries;
            }
            _num_StringItems[_num_StringEntries-1] = 0;
        }


        //////
            
        void SampleInfo::addStringItem() {
            _num_StringItems[_num_StringEntries-1]++;
            if (_num_StringItems[_num_StringEntries-1] > _res_StringItems[_num_StringEntries-1]) {
                _res_len_StringItems[_num_StringEntries-1] = (unsigned int *) realloc(_res_len_StringItems[_num_StringEntries-1], _num_StringItems[_num_StringEntries-1] * sizeof(unsigned int));
                if (!_res_len_StringItems[_num_StringEntries-1]) throw EGGMEM;
                _res_len_StringItems[_num_StringEntries-1][_num_StringItems[_num_StringEntries-1]-1] = 1;  // last item of last entry
                _StringItems[_num_StringEntries-1] = (char **) realloc(_StringItems[_num_StringEntries-1], _num_StringItems[_num_StringEntries-1] * sizeof(char *));
                if (!_StringItems[_num_StringEntries-1]) throw EGGMEM;
                _StringItems[_num_StringEntries-1][_num_StringItems[_num_StringEntries-1]-1] = (char *) malloc(1 * sizeof(char));  // last item of last entry
                if (!_StringItems[_num_StringEntries-1][_num_StringItems[_num_StringEntries-1]-1]) throw EGGMEM;
                _StringItems[_num_StringEntries-1][_num_StringItems[_num_StringEntries-1]-1][0] = '\0';
                _res_StringItems[_num_StringEntries-1] = _num_StringItems[_num_StringEntries-1];
            }
        }
    }

    ////////////////////////////////////////////////////////////////////

    VcfParser::VcfParser() {
        init();
    }

    VcfParser::~VcfParser() {
        free();
    }

    void VcfParser::open_file(const char* fname) {

        // reset memory

        reset();
        
        // save file name

        unsigned int lfname = strlen(fname) + 1;

        if (lfname > _res_fname) {
            _fname = (char *) realloc(_fname, lfname * sizeof(char));
            if (!_fname) throw EGGMEM;
            _res_fname = lfname;
        }

        strcpy(_fname, fname);
        
        // open file
        
        _stream = &_localstream;
        _localstream.open(fname);
        if (!_localstream.is_open()) {
            throw EggOpenFileError(fname);
        }

        // parse header
        header();
        _idx_frt_sample=get_index();
        _idx_end_file = file_end();
    }

    void VcfParser::set_stream(std::istream& stream) {

        // reset memory

        reset();

        // save file name

        if (_res_fname < 9) {
            _fname = (char *) realloc(_fname, 9 * sizeof(char));
            if (!_fname) throw EGGMEM;
            _res_fname = 9;
        }

        strcpy(_fname, "<stream>");

        // get stream
        
        _stream = &stream;

        // parse header

        header();
        _idx_frt_sample=get_index();
        _idx_end_file = file_end();
    }

    void VcfParser::clear_stream(){
        _stream->clear();
    }

    long long int VcfParser::get_index() const{
        if(_stream->eof()) _stream->clear();
        return _stream->tellg();
    }

    void VcfParser::set_currline(unsigned int line){
        _currline=_start_currline+line;
    }

    void VcfParser::set_index(long long int index){
        if(_stream->eof()) _stream->clear();
        _stream->seekg(index,_stream->beg);
    }

    void VcfParser::set_previous_next(){
        if(_stream->good()) _previous_next=_stream->tellg();
        else _previous_next=_idx_frt_sample;
    }

    long long int VcfParser::file_end(){
        long long int cur_stream=_stream->tellg();
        _stream->seekg (0, _stream->end);
        long long int end_stream=_stream->tellg();
        _stream->seekg(cur_stream,_stream->beg);
        return end_stream;
    }

    void VcfParser::unread(){
        if((_previous_next!=0) && (_position != UNKNOWN))
            set_index(_previous_next);
        else
            throw EggException("unexpected use of method: you can't call the unread() method without have called the read() method before");
    }

    void VcfParser::read_index(const char *fname){
        unsigned int lfname=strlen(fname)+1;
        ifstream f_stream(fname);
        string line;
        if(f_stream.is_open()){
            std::getline(f_stream,line);
            if (!strcmp(line.c_str(), "#Egglib_index_file\n")) throw EggFormatError(fname, 1, "index file", "the first line of the index file is not correct");
                if(file_end() == egglib::get_index_eof(fname)){
                    _index=(char *) realloc (_index, lfname * sizeof(char*));
            if (!_index) throw EGGMEM;
                strcpy(_index,fname);
                _has_index=true;
            }else{
                throw EggFormatError(fname, 1, "index file", "end file false");
            }
        }else{
            throw EggRuntimeError("cannot open the index file");
        }
        f_stream.close();
    }

    bool VcfParser::has_index() {
        return _has_index;
    }

    const char * VcfParser::index_fname() {
        if(_index){ //if(_has_index(){
            return _index;
        }else{
            throw EggRuntimeError("This VcfParser have no index file linked. Use the method 'read_index' to get or load an index file corresponding to the VcfParser");
        }    
    }

    void VcfParser::rewind(){
        set_index(_idx_frt_sample);
    }

    long long int VcfParser::get_idx_frt_sample(){
        return _idx_frt_sample;
    }

    void VcfParser::read_header(const char * string) {

        // reset memory
        reset();

        // save file name
        if (_res_fname < 9) {
            _fname = (char *) realloc(_fname, 9 * sizeof(char));
            if (!_fname) throw EGGMEM;
            _res_fname = 9;
        }
        strcpy(_fname, "<string>");

        // activate string stream
        _stringstream.clear();
        _stringstream.str(string);
        _stream = &_stringstream;

        // parse header
        header();

        // invalidate stream
        _stream = NULL;
    }
    
    void VcfParser::predefine() {
        add_info()->update("AA",         1,                  vcf::Info::String,  "ancestral allele");
        add_info()->update("AC",         vcf::NUM_ALTERNATE, vcf::Info::Integer, "allele count in genotypes, for each ALT allele, in the same order as listed");
        add_info()->update("AF",         vcf::NUM_ALTERNATE, vcf::Info::Float,   "allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes");
        add_info()->update("AN",         1,                  vcf::Info::Integer, "total number of alleles in called genotypes");
        add_info()->update("BQ",         1,                  vcf::Info::Integer, "RMS base quality at this position");
        add_info()->update("CIGAR",      vcf::NUM_ALTERNATE, vcf::Info::String,  "cigar string describing how to align an alternate allele to the reference allele");
        add_info()->update("DB",         0,                  vcf::Info::Flag,    "dbSNP membership");
        add_info()->update("DP",         1,                  vcf::Info::Integer, "combined depth across samples, e.g. DP=154");
        add_info()->update("END",        1,                  vcf::Info::Integer, "end position of the variant described in this record (for use with symbolic alleles)");
        add_info()->update("H2",         0,                  vcf::Info::Flag,    "membership in hapmap2");
        add_info()->update("H3",         0,                  vcf::Info::Flag,    "membership in hapmap3");
        add_info()->update("MQ",         1,                  vcf::Info::Integer, "RMS mapping quality, e.g. MQ=52");
        add_info()->update("MQ0",        1,                  vcf::Info::Integer, "number of MAPQ == 0 reads covering this record");
        add_info()->update("NS",         1,                  vcf::Info::Integer, "number of samples with data");
        add_info()->update("SB",         1,                  vcf::Info::Float,   "strand bias at this position");
        add_info()->update("SOMATIC",    0,                  vcf::Info::Flag,    "indicates that the record is a somatic mutation, for cancer genomics");
        add_info()->update("VALIDATED",  0,                  vcf::Info::Flag,    "validated by follow-up experiment");
        add_info()->update("1000G",      0,                  vcf::Info::Flag,    "membership in 1000 Genomes");
        add_info()->update("IMPRECISE",  0,                  vcf::Info::Flag,    "imprecise structural variation");
        add_info()->update("NOVEL",      0,                  vcf::Info::Flag,    "Indicates a novel structural variation");
// NO REDEFINE! add_info()->update("END",1,                  vcf::Info::Integer, "end position of the variant described in this record");
        add_info()->update("SVTYPE",     1,                  vcf::Info::String,  "type of structural variant");
        add_info()->update("SVLEN",      UNKNOWN,            vcf::Info::Integer, "difference in length between REF and ALT alleles");
        add_info()->update("CIPOS",      2,                  vcf::Info::Integer, "confidence interval around POS for imprecise variants");
        add_info()->update("CIEND",      2,                  vcf::Info::Integer, "confidence interval around END for imprecise variants");
        add_info()->update("HOMLEN",     UNKNOWN,            vcf::Info::Integer, "length of base pair identical micro-homology at event breakpoints");
        add_info()->update("HOMSEQ",     UNKNOWN,            vcf::Info::String,  "sequence of base pair identical micro-homology at event breakpoints");
        add_info()->update("BKPTID",     UNKNOWN,            vcf::Info::String,  "ID of the assembled alternate allele in the assembly file");
        add_info()->update("MEINFO",     4,                  vcf::Info::String,  "mobile element info of the form NAME,START,END,POLARITY");
        add_info()->update("METRANS",    4,                  vcf::Info::String,  "mobile element transduction info of the form CHR,START,END,POLARITY");
        add_info()->update("DGVID",      1,                  vcf::Info::String,  "ID of this element in Database of Genomic Variation");
        add_info()->update("DBVARID",    1,                  vcf::Info::String,  "ID of this element in DBVAR");
        add_info()->update("DBRIPID",    1,                  vcf::Info::String,  "ID of this element in DBRIP");
        add_info()->update("MATEID",     UNKNOWN,            vcf::Info::String,  "ID of mate breakends");
        add_info()->update("PARID",      1,                  vcf::Info::String,  "ID of partner breakend");
        add_info()->update("EVENT",      1,                  vcf::Info::String,  "ID of event associated to breakend");
        add_info()->update("CILEN",      2,                  vcf::Info::Integer, "confidence interval around the length of the inserted material between breakends");
// NO REDEFINE! add_info()->update("DP", 1,                  vcf::Info::Integer, "read depth of segment containing breakend");
        add_info()->update("DPADJ",      UNKNOWN,            vcf::Info::Integer, "read depth of adjacency");
        add_info()->update("CN",         1,                  vcf::Info::Integer, "copy number of segment containing breakend");
        add_info()->update("CNADJ",      UNKNOWN,            vcf::Info::Integer, "copy number of adjacency");
        add_info()->update("CICN",       2,                  vcf::Info::Integer, "confidence interval around copy number for the segment");
        add_info()->update("CICNADJ",    UNKNOWN,            vcf::Info::Integer, "confidence interval around copy number for the adjacency");

        add_format()->update("GT",   1,                  vcf::Info::String,  "genotype, encoded as allele values separated by either of \"/\" or \"|\". The allele values are 0 for the reference allele (what is in the REF field), 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on. For diploid calls examples could be 0/1, 1|0, or 1/2, etc. For haploid calls, e.g. on Y, male non-pseudoautosomal X, or mitochondrion, only one allele value should be given; a triploid call might look like 0/0/1. If a call cannot be made for a sample at a given locus, \".\" should be specified for each missing allele in the GT field (for example \"./.\" for a diploid genotype and \".\" for haploid genotype). The meanings of the separators are as follows (see the PS field below for more details on incorporating phasing information into the genotypes): \"/\" : genotype unphased, \"|\" : genotype phased");
        add_format()->update("DP",   1,                  vcf::Info::Integer, "read depth at this position for this sample");
        add_format()->update("FT",   1,                  vcf::Info::String,  "sample genotype filter indicating if this genotype was \"called\" (similar in concept to the FILTER field). Again, use PASS to indicate that all filters have been passed, a semi-colon separated list of codes for filters that fail, or \".\" to indicate that filters have not been applied. These values should be described in the meta-information in the same way as FILTERs (no white-space or semi-colons permitted)");
        add_format()->update("GL",   UNKNOWN,            vcf::Info::Float,   "genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields. In presence of the GT field the same ploidy is expected and the canonical order is used; without GT field, diploidy is assumed. If A is the allele in REF and B,C,... are the alleles as ordered in ALT, the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j.  In other words, for biallelic sites the ordering is: AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc.  For example: GT:GL 0/1:-323.03,-99.29,-802.53");
        add_format()->update("GLE",  1,                  vcf::Info::String,  "genotype likelihoods of heterogeneous ploidy, used in presence of uncertain copy number. For example: GLE=0:-75.22,1:-223.42,0/0:-323.03,1/0:-99.29,1/1:-802.53");
        add_format()->update("PL",   UNKNOWN,            vcf::Info::Integer, "the phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely as the GL field)");
        add_format()->update("GP",   UNKNOWN,            vcf::Info::Float,   "the phred-scaled genotype posterior probabilities (and otherwise defined precisely as the GL field); intended to store imputed genotype probabilities");
        add_format()->update("GQ",   1,                  vcf::Info::Integer, "conditional genotype quality, encoded as a phred quality -10log_10p(genotype call is wrong, conditioned on the site's being variant)");
        add_format()->update("HQ",   2,                  vcf::Info::Integer, "haplotype qualities, two comma separated phred qualities");
        add_format()->update("PS",   1,                  vcf::Info::Integer, "phase set.  A phase set is defined as a set of phased genotypes to which this genotype belongs.  Phased genotypes for an individual that are on the same chromosome and have the same PS value are in the same phased set.  A phase set specifies multi-marker haplotypes for the phased genotypes in the set.  All phased genotypes that do not contain a PS subfield are assumed to belong to the same phased set.  If the genotype in the GT field is unphased, the corresponding PS field is ignored.  The recommended convention is to use the position of the first variant in the set as the PS identifier (although this is not required). (Non-negative 32-bit Integer)");
        add_format()->update("PQ",   1,                  vcf::Info::Integer, "phasing quality, the phred-scaled probability that alleles are ordered incorrectly in a heterozygote (against all other members in the phase set).  We note that we have not yet included the specific measure for precisely defining \"phasing quality\"; our intention for now is simply to reserve the PQ tag for future use as a measure of phasing quality.");
        add_format()->update("EC",   vcf::NUM_ALTERNATE, vcf::Info::Integer, "comma separated list of expected alternate allele counts for each alternate allele in the same order as listed in the ALT field (typically used in association analyses)");
        add_format()->update("MQ",   1,                  vcf::Info::Integer, "RMS mapping quality, similar to the version in the INFO field.");
        add_format()->update("CN",   1,                  vcf::Info::Integer, "copy number genotype for imprecise events");
        add_format()->update("CNQ",  1,                  vcf::Info::Float,   "copy number genotype quality for imprecise events");
        add_format()->update("CNL",  UNKNOWN,            vcf::Info::Float,   "copy number genotype likelihood for imprecise events");
        add_format()->update("NQ",   1,                  vcf::Info::Integer, "phred style probability score that the variant is novel with respect to the genome's ancestor");
        add_format()->update("HAP",  1,                  vcf::Info::Integer, "unique haplotype identifier");
        add_format()->update("AHAP", 1,                  vcf::Info::Integer, "unique identifier of ancestral haplotype");
        
        add_alt()->update("DEL",        "deletion relative to the reference");
        add_alt()->update("INS",        "insertion of novel sequence relative to the reference");
        add_alt()->update("DUP",        "region of elevated copy number relative to the reference");
        add_alt()->update("INV",        "inversion of reference sequence");
        add_alt()->update("CNV",        "copy number variable region (may be both deletion and duplication)");
        add_alt()->update("DUP:TANDEM", "tandem duplication");
        add_alt()->update("DEL:ME",     "deletion of mobile element relative to the reference");
        add_alt()->update("INS:ME",     "insertion of a mobile element relative to the reference");
    }

    void VcfParser::reset() {
        if (_localstream.is_open()) _localstream.close();
        _stream = NULL;
        _idx_end_file = 0;
        _index[0]= '\0';
        _has_index=false;
        _idx_frt_sample=0;
        _fname[0] = '\0';
        _currline = 0;
        _start_currline = 0;
        _ff[0] = '\0';
        _num_info = 0;
        _num_filter = 0;
        _num_format = 0;
        _num_meta = 0;
        _num_samples = 0;
        _num_alt = 0;
        curr_ch = 0;
        prev_ch = 0;
        predefine();
        _allow_X = false;
        _allow_gap = false;
        _threshold_PL = UNKNOWN;
        _threshold_GL = UNKNOWN;
        reset_variant();
    }
    
    void VcfParser::reset_variant() {
        _chrom[0] = '\0';
        _position = UNKNOWN;
        _num_ID = 0;
        _len_reference = 0;
        _num_alternate = 0;
        _num_genotypes = UNKNOWN;
        _ploidy = UNKNOWN;
        _quality = UNDEF;
        _num_failed_test = UNKNOWN;
        _num_FlagInfo = 0;
        _num_CharacterInfo = 0;
        _num_IntegerInfo = 0;
        _num_FloatInfo = 0;
        _num_StringInfo = 0;
        _num_formatEntries = 0;
        _has_AN = false;
        _has_AC = false;
        _has_AC_ss = false;
        _has_AF = false;
        _has_AA = false;
        _has_GT = false;
        _has_PL = false;
        _has_GL = false;
        _GT_all_phased = true;
        _has_data = false;
    }

    void VcfParser::clear() {
        free();
        init();
    }
    
    void VcfParser::header() {

        // read fileformat

        next();    if (curr_ch != '#') throw EggFormatError(_fname, _currline+1, "VCF", "first character of file is not \"#\" as expected", curr_ch);
        next();    if (curr_ch != '#') throw EggFormatError(_fname, _currline+1, "VCF", "second character of file is not \"#\" as expected", curr_ch);
        get_string(_buffer, _res_buffer, &VcfParser::check_letter, &VcfParser::stop_equal, false);
        if (strcmp(_buffer, "fileformat")) throw EggFormatError(_fname, _currline+1, "VCF", "first meta-information of the file is not \"fileformat\"");

        get_string(_ff, _res_ff, &VcfParser::check_sign_and_space, &VcfParser::stop_line, false);
        _currline++;
            
        // read all meta-information until the header line (detected as `#C`)

        while (true) {

            next();    if (curr_ch != '#') throw EggFormatError(_fname, _currline+1, "VCF", "first character of line is not \"#\" as expected", curr_ch);
            next();
            if (curr_ch != '#') {
                if (curr_ch == 'C') break;
                throw EggFormatError(_fname, _currline+1, "VCF", "second character of line is not \"#\" as expected", curr_ch);
            }

            //get_string(_buffer, _res_buffer, &VcfParser::check_letter, &VcfParser::stop_equal, false);
            get_string(_buffer, _res_buffer, &VcfParser::check_alphanumericunderscore, &VcfParser::stop_equal, false);

            // check is not new fileformat

            if (!strcmp(_buffer, "fileformat")) throw EggFormatError(_fname, _currline+1, "VCF", "fileformat meta-information cannot be specified more than once");

            // INFO

            if (!strcmp(_buffer, "INFO")) {
                get_4_fields(false);
                continue;
            }
        
            // FORMAT

            if (!strcmp(_buffer, "FORMAT")) {
                get_4_fields(true);
                continue;
            }

            // FILTER

            if (!strcmp(_buffer, "FILTER")) {
                get_2_fields(false);
                continue;
            }

            // ALT

            if (!strcmp(_buffer, "ALT")) {
                get_2_fields(true);
                continue;
            }

            // otherwise it is a Meta
            
            vcf::Meta& meta = *add_meta();
            meta.set_key(_buffer);
            get_string(_buffer, _res_buffer, &VcfParser::check_sign_and_space, &VcfParser::stop_line, false);
            meta.set_value(_buffer);
            _currline++;
        }

        // check the the current line is indeed the header
        
        get_string(_buffer, _res_buffer, &VcfParser::check_sign_and_space, &VcfParser::stop_tab, false);
        if (strcmp(_buffer, "HROM")) throw EggFormatError(_fname, _currline+1, "VCF", "invalid header line");

        get_string(_buffer, _res_buffer, &VcfParser::check_sign_and_space, &VcfParser::stop_tab, false);
        if (strcmp(_buffer, "POS")) throw EggFormatError(_fname, _currline+1, "VCF", "invalid header line");

        get_string(_buffer, _res_buffer, &VcfParser::check_sign_and_space, &VcfParser::stop_tab, false);
        if (strcmp(_buffer, "ID")) throw EggFormatError(_fname, _currline+1, "VCF", "invalid header line");

        get_string(_buffer, _res_buffer, &VcfParser::check_sign_and_space, &VcfParser::stop_tab, false);
        if (strcmp(_buffer, "REF")) throw EggFormatError(_fname, _currline+1, "VCF", "invalid header line");

        get_string(_buffer, _res_buffer, &VcfParser::check_sign_and_space, &VcfParser::stop_tab, false);
        if (strcmp(_buffer, "ALT")) throw EggFormatError(_fname, _currline+1, "VCF", "invalid header line");

        get_string(_buffer, _res_buffer, &VcfParser::check_sign_and_space, &VcfParser::stop_tab, false);
        if (strcmp(_buffer, "QUAL")) throw EggFormatError(_fname, _currline+1, "VCF", "invalid header line");

        get_string(_buffer, _res_buffer, &VcfParser::check_sign_and_space, &VcfParser::stop_tab, false);
        if (strcmp(_buffer, "FILTER")) throw EggFormatError(_fname, _currline+1, "VCF", "invalid header line");

        get_string(_buffer, _res_buffer, &VcfParser::check_sign_and_space, &VcfParser::stop_linetabEOF, false);
        if (strcmp(_buffer, "INFO")) throw EggFormatError(_fname, _currline+1, "VCF", "invalid header line");

        _currline++;

       // read the putative format

        if (_stream->eof()) return;
        if (curr_ch == '\n') return;
        if (curr_ch != '\t') throw EggFormatError(_fname, _currline, "VCF", "invalid header line");
        
        get_string(_buffer, _res_buffer, &VcfParser::check_sign_and_space, &VcfParser::stop_tab, false);
        if (strcmp(_buffer, "FORMAT")) throw EggFormatError(_fname, _currline, "VCF", "invalid header line");
        
       // read the samples (by default: _num_samples = 0)
        while (true) {
            
            unsigned int n = get_string(_buffer, _res_buffer, &VcfParser::check_sign_and_space, &VcfParser::stop_linetabEOF, false);
            n++; // for termination character
            _num_samples++;
            
            if (_num_samples > _res_samples) {
                _res_len_samples = (unsigned int *) realloc(_res_len_samples, _num_samples * sizeof(unsigned int));
                if (!_res_len_samples) throw EGGMEM;
                _res_len_samples[_num_samples-1] = 0;
                _samples = (char **) realloc(_samples, _num_samples * sizeof(char *));
                if (!_samples) throw EGGMEM;
                _samples[_num_samples-1] = NULL;
                _sampleInfo = (vcf::SampleInfo **) realloc(_sampleInfo, _num_samples * sizeof(vcf::SampleInfo *));
                if (!_sampleInfo) throw EGGMEM;
                _sampleInfo[_num_samples-1] = new(std::nothrow) vcf::SampleInfo;
                if (!_sampleInfo[_num_samples-1]) throw EGGMEM;
                _res_samples = _num_samples;
                _GT = (unsigned int **) realloc(_GT, _num_samples * sizeof(unsigned int *));
                if (!_GT) throw EGGMEM;
                _GT[_num_samples-1] = NULL;
                _PL = (unsigned int **) realloc(_PL, _num_samples * sizeof(unsigned int *));
                if (!_PL) throw EGGMEM;
                _PL[_num_samples-1] = NULL;
                _res_GT = (unsigned int *) realloc(_res_GT, _num_samples * sizeof(unsigned int));
                if (!_res_GT) throw EGGMEM;
                _res_GT[_num_samples-1] = 0;
                _res_PL = (unsigned int *) realloc(_res_PL, _num_samples * sizeof(unsigned int));
                if (!_res_PL) throw EGGMEM;
                _res_PL[_num_samples-1] = 0;
                _GL = (double **) realloc(_GL, _num_samples * sizeof(double *));
                if (!_GL) throw EGGMEM;
                _GL[_num_samples-1] = NULL;
                _res_GL = (unsigned int *) realloc(_res_GL, _num_samples * sizeof(unsigned int));
                if (!_res_GL) throw EGGMEM;
                _res_GL[_num_samples-1] = 0;
                _GT_phased = (bool *) realloc(_GT_phased, _num_samples * sizeof(bool));
                if (!_GT_phased) throw EGGMEM;
                _outgroup_map = (bool *) realloc(_outgroup_map, _num_samples * sizeof(bool));
                if (!_outgroup_map) throw EGGMEM;
                _outgroup_map[_num_samples - 1] = false;
            }
            
            if (n > _res_len_samples[_num_samples-1]) {
                _samples[_num_samples-1] = (char *) realloc(_samples[_num_samples-1], n * sizeof(char));
                if (!_samples[_num_samples-1]) throw EGGMEM;
                _res_len_samples[_num_samples - 1] = n;
            }
            
            strcpy(_samples[_num_samples-1], _buffer);

            if (_stream->eof()) break;
            if (curr_ch == '\n') break;
            if (curr_ch != '\t') throw EggFormatError(_fname, _currline, "VCF", "invalid header line");
        }

        _start_currline = _currline-1;
    }

    bool VcfParser::is_outgroup(unsigned int x) const {
        return _outgroup_map[x];
    }

    void VcfParser::set_outgroup(unsigned int x) {
        _outgroup_map[x] = true;
    }

    void VcfParser::allow_X(bool flag) {
        _allow_X = flag;
    }

    void VcfParser::allow_gap(bool flag) {
        _allow_gap = flag;
    }

    void VcfParser::read_line(const char * string) {
        _stringstream.clear();
        _stringstream.str(string);
        _stream = &_stringstream;
        read();
        _stream = NULL;
    }

    void VcfParser::read() {

        _previous_next=get_index();
        reset_variant();

        // checking
        if (_stream == NULL || (!_stream->good())) throw EggFormatError(_fname, _currline+1, "VCF", "error in file or invalid stream");

        // get chromosome
        get_string(_chrom, _res_chrom, &VcfParser::check_sign, &VcfParser::stop_tab, false);

        // get position
        unsigned int n = get_string(_buffer, _res_buffer, &VcfParser::check_integer, &VcfParser::stop_tab, true);
        if (_buffer[0] == MAXCHAR) throw EggFormatError(_fname, _currline+1, "VCF", "position cannot be missing");

        _position = atol(_buffer);
        if (_position == 0) _position = BEFORE;
        else _position--;

        // get ID(s)
        while (true) {
            _num_ID++;
            if (_num_ID > _res_ID) {

                _ID = (char **) realloc(_ID, _num_ID * sizeof(char *));
                if (!_ID) throw EGGMEM;
                
                _ID[_num_ID-1] = NULL;
                _res_len_ID = (unsigned int *) realloc(_res_len_ID, _num_ID * sizeof(unsigned int));
                _res_len_ID[_num_ID-1] = 0;
                
                _res_ID = _num_ID;
            }
            get_string(_ID[_num_ID-1], _res_len_ID[_num_ID-1], &VcfParser::check_sign_and_space, &VcfParser::stop_tabsemicolon, false);

            if (curr_ch == '\t') break;            
        }
        
        // get reference allele
        _len_reference = get_string(_reference, _res_reference, &VcfParser::check_bases, &VcfParser::stop_tab, false);
        
        // get alternate allele(s) (defines _num_alternate)
        while (true) {

            // read next allele
            unsigned int n = get_string(_buffer, _res_buffer, &VcfParser::check_sign_and_space, &VcfParser::stop_tabcomma, false);

            // is it a missing value?
            if (!strcmp(_buffer, ".")) {
                if (_num_alternate != 0) throw EggFormatError(_fname, _currline+1, "VCF", "one out of several alternate alleles declared as missing: this is illegal");
                if (curr_ch != '\t') throw EggFormatError(_fname, _currline+1, "VCF", "one out of several alternate alleles declared as missing: this is illegal");
                break;
            }

            // alloc new allele in all other cases
            _num_alternate++;
            if (_num_alternate > _res_alternate) {
                _type_alternate = (vcf::AltType *) realloc(_type_alternate, _num_alternate * sizeof(vcf::AltType));
                if (!_type_alternate) throw EGGMEM;
                
                _res_len_alternate = (unsigned int *) realloc(_res_len_alternate, _num_alternate * sizeof(unsigned int));
                if (!_res_len_alternate) throw EGGMEM;
                _res_len_alternate[_num_alternate-1] = 0;
                
                _alternate = (char **) realloc(_alternate, _num_alternate * sizeof(char *));
                if (!_alternate) throw EGGMEM;
                _alternate[_num_alternate-1] = NULL;
                
                _res_alternate = _num_alternate;
            }

            if (_buffer[0] == '<' && _buffer[n-1] == '>') {
                _type_alternate[_num_alternate-1] = vcf::Referred;
                if (n-1 > _res_len_alternate[_num_alternate-1]) {
                    _alternate[_num_alternate-1] = (char *) realloc(_alternate[_num_alternate-1], (n-1) * sizeof(char));
                    if (!_alternate[_num_alternate-1]) throw EGGMEM;
                    _res_len_alternate[_num_alternate-1] = n-1;
                }
                strncpy(_alternate[_num_alternate-1], _buffer+1, n-2);
                _alternate[_num_alternate-1][n-2] = '\0';

                if (curr_ch == '\t') break;
                continue;
            }

            // copy full buffer in all other cases
            if ((n+1) > _res_len_alternate[_num_alternate-1]) {
                _alternate[_num_alternate-1] = (char *) realloc(_alternate[_num_alternate-1], (n+1) * sizeof(char));
                if (!_alternate[_num_alternate-1]) throw EGGMEM;
                _res_len_alternate[_num_alternate-1] = n+1;
            }
            strcpy(_alternate[_num_alternate-1], _buffer);

            // if this a breakend?   recognized by two `]` or two `[`, one `:`
            unsigned int colon = 0;
            unsigned int open = 0;
            unsigned int close = 0;
            
            for (unsigned int i=0; i<n; i++) {
                if (_buffer[i] == '[') open++;
                if (_buffer[i] == ']') close++;
                if (_buffer[i] == ':') colon++;
            }
            
            if (colon==1 && ( (open==2 && close==0) || (open==0 && close==2) )) {
                _type_alternate[_num_alternate-1] = vcf::Breakend;
                if (curr_ch == '\t') break;
                continue;
            }

            // if this a simple X? (only if allow_X has been set)
            if (_allow_X == true && n == 1 && (_buffer[0] == 'X' or _buffer[0] == 'x')) {
                _type_alternate[_num_alternate-1] = vcf::X;
                if (_buffer[0] == 'x') _alternate[_num_alternate-1][0] = 'X'; // ported to upper case
                if (curr_ch == '\t') break;
                continue;
            }

            // finally, it should be an explicit allele
            for (unsigned int i=0; i<n; i++) {
                switch (_buffer[i]) {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                    case 'N':
                        continue;
                    case '*':
                        _alternate[_num_alternate-1][i] = '-'; // overlapping gaps are recoded as alignment gap
                        continue;
                    case 'a':
                    case 'c':
                    case 'g':
                    case 't':
                    case 'n':
                        _alternate[_num_alternate-1][i] = toupper(_buffer[i]); // all is ported to upper case
                        continue;
                    case '-':
                        if (_allow_gap) continue;
                    default:
                        throw EggFormatError(_fname, _currline+1, "VCF", "the following alternate allele is invalid: ", '\0', _alternate[_num_alternate-1]);
                }
            }
            
            _type_alternate[_num_alternate-1] = vcf::Default;
            if (curr_ch == '\t') break;
        }

        // get quality
        n = get_string(_buffer, _res_buffer, &VcfParser::check_float, &VcfParser::stop_tab, true);

        if (_buffer[0] == MAXCHAR) {
            // leave quality to default value (UNDEF)
        }

        else {
            unsigned int exp = 0;
            unsigned int dot = 0;
            for (unsigned int i=0; i<n; i++) {
                if (_buffer[i] == 'e' || _buffer[i] == 'E') exp++;
                if (_buffer[i] == '.') dot++;
                if (_buffer[i] == '-' && i != 0 && _buffer[i-1] != 'E' && _buffer[i-1] != 'e') throw EggFormatError(_fname, _currline+1, "VCF", "invalid numeric value: ", '\0', _buffer);
            }
            if (exp > 1 || dot > 1) throw EggFormatError(_fname, _currline+1, "VCF", "invalid numeric value: ", '\0', _buffer);
            
            _quality = atof(_buffer);
        }

        // get filter results
        _num_failed_test = 0;
        while (true) {

            // read next test
            unsigned int n = get_string(_buffer, _res_buffer, &VcfParser::check_sign, &VcfParser::stop_tabsemicolon, true);

            // check special case of no tests or no failed tests
            if (_buffer[0] == MAXCHAR) {
                if (_num_failed_test != 0) throw EggFormatError(_fname, _currline+1, "VCF", "one out of several filters is the missing value (\".\"): this is illegal");
                if (curr_ch != '\t') throw EggFormatError(_fname, _currline+1, "VCF", "one out of several filters is the missing value (\".\"): this is illegal", curr_ch);
                _num_failed_test = UNKNOWN;
                break;
            }

            if (!strcmp(_buffer, "PASS")) {
                if (_num_failed_test != 0) throw EggFormatError(_fname, _currline+1, "VCF", "one out of several filters is PASS: this is illegal");
                if (curr_ch != '\t') throw EggFormatError(_fname, _currline+1, "VCF", "one out of several filters is PASS: this is illegal", curr_ch);
                break;
            }

            // alloc new filter in all other cases
            _num_failed_test++;
            if (_num_failed_test > _res_failed_test) {
                _res_len_failed_test = (unsigned int *) realloc(_res_len_failed_test, _num_failed_test * sizeof(unsigned int));
                if (!_res_len_failed_test) throw EGGMEM;
                _res_len_failed_test[_num_failed_test-1] = 0;
                
                _failed_test = (char **) realloc(_failed_test, _num_failed_test * sizeof(char *));
                if (!_failed_test) throw EGGMEM;
                _failed_test[_num_failed_test-1] = NULL;
                
                _res_failed_test = _num_failed_test;
            }

            if ((n+1) > _res_len_failed_test[_num_failed_test-1]) {
                _failed_test[_num_failed_test-1] = (char *) realloc(_failed_test[_num_failed_test-1], (n+1) * sizeof(char));
                if (!_failed_test[_num_failed_test-1]) throw EGGMEM;
                _res_len_failed_test[_num_failed_test-1] = n+1;
            }
            strcpy(_failed_test[_num_failed_test-1], _buffer);

            if (curr_ch == '\t') break;
        }

        // get info fields
        while (true) {

            // get key
            get_string(_buffer, _res_buffer, &VcfParser::check_sign, &VcfParser::stop_equalsemicolontablineEOF, true);

            // detect missing block
            if (_buffer[0] == MAXCHAR) {
                if (curr_ch == '\t' || curr_ch == '\n') break;
                else throw EggFormatError(_fname, _currline+1, "VCF", "INFO ID cannot be a single period", '\0');
            }

            // locate it amongst info specifications
            vcf::Info * info = find_info(_buffer);

            // default values (if unspecified)
            vcf::Info::Type type = vcf::Info::String;
            unsigned int number = UNKNOWN;

            if (info) {
                type = info->get_type();
                number = info->get_number();
            }

            vcf::FlagInfo * next = 0;
            switch(type) {
                case vcf::Info::Flag:
                    next = add_FlagInfo(number);
                    break;
                case vcf::Info::Integer:
                    next = add_IntegerInfo(number);
                    break;
                case vcf::Info::Float:
                    next = add_FloatInfo(number);
                    break;
                case vcf::Info::Character:
                    next = add_CharacterInfo(number);
                    break;
                case vcf::Info::String:
                    next = add_StringInfo(number);
                    break;
            }

            next->set_ID(_buffer);

            // if type if Flag, there should be no values, otherwise, must be value(s)
            if (type==vcf::Info::Flag) {
                if (curr_ch == '\t') break;
                if (curr_ch == '\n') break;
                if (_stream->eof()) break;
                if (curr_ch != ';') throw EggFormatError(_fname, _currline+1, "VCF", "no values allowed for Flag-type INFO - ID: ", '\0', _buffer);
                continue;
            }

            if (curr_ch != '=') throw EggFormatError(_fname, _currline+1, "VCF", "values are expected for non-Flag-type INFO ", '\0', _buffer);

            // get value(s)
            unsigned int n = 0;
            while (true) {

                // receptacles to simplify handling of pointers
                vcf::StringInfo * s;
                vcf::TypeInfo<int> * i;
                vcf::TypeInfo<char> * c;
                vcf::TypeInfo<double> * f;

                switch (type) {
                    case vcf::Info::String:
                        s = (vcf::StringInfo *) next;
                        s->add();
                        get_string(s->_items[n], s->_res_len_items[n], &VcfParser::check_sign, &VcfParser::stop_tabsemicoloncommalineEOF, true);
                        break;
                        
                    case vcf::Info::Integer:
                        i = (vcf::TypeInfo<int> *) next;
                        i->add();
                        get_string(_buffer, _res_buffer, &VcfParser::check_integer, &VcfParser::stop_tabsemicoloncommalineEOF, true);
                        if (_buffer[0] == MAXCHAR) i->_items[n] = MISSINGDATA;
                        else i->_items[n] = atoi(_buffer);
                        break;
                        
                    case vcf::Info::Float:
                        f = (vcf::TypeInfo<double> *) next;
                        f->add();
                        get_string(_buffer, _res_buffer, &VcfParser::check_float, &VcfParser::stop_tabsemicoloncommalineEOF, true);
                        if (_buffer[0] == MAXCHAR) f->_items[n] = UNDEF;
                        else f->_items[n] = atof(_buffer);
                        break;
                    
                    case vcf::Info::Character:
                        c = (vcf::TypeInfo<char> *) next;
                        c->add();
                        get_string(_buffer, _res_buffer, &VcfParser::check_integer, &VcfParser::stop_tabsemicoloncommalineEOF, true);
                        if (strlen(_buffer) != 1) throw EggFormatError(_fname, _currline+1, "VCF", "invalid value for Character-type INFO: ", '\0', _buffer);
                        c->_items[n] = _buffer[0];
                        break;

                    case vcf::Info::Flag:
                        throw EggRuntimeError("unexpected error in egglib::VcfParser::read() [1]");
                }

                n++;
                
                if (curr_ch == ';') break;
                if (curr_ch == '\t') break;
                if (curr_ch == '\n') break;
                if (_stream->eof()) break;
            }

            // check that the number of values is correct (don't modify number since it will check as special value)
            unsigned int actual_number = number;
            if (actual_number == vcf::NUM_ALTERNATE) actual_number = _num_alternate;
            if (actual_number == vcf::NUM_POSSIBLE_ALLELES) actual_number = _num_alternate + 1;
            if (actual_number == vcf::NUM_GENOTYPES) actual_number = (_num_alternate + 1) * (_num_alternate + 2) / 2; // at this point number of genotypes assuming diploidy
            if (actual_number != UNKNOWN && n != actual_number) throw EggFormatError(_fname, _currline+1, "VCF", "incorrect number of arguments for INFO entry: ", '\0', next->get_ID());
            
            // intercept AN, AC, AF and AA
            if (!strcmp(next->get_ID(), "AN") && type==vcf::Info::Integer && number==1) {
                _has_AN = true;
                _AN = ((vcf::TypeInfo<int> *) next)->item(0);
                
                if (_has_AC_ss) {
                    _has_AC = true;
                }
            }

            if (!strcmp(next->get_ID(), "AC") && type==vcf::Info::Integer && number==vcf::NUM_ALTERNATE) {

                _has_AC_ss = true;

                _num_AC = _num_alternate;
                if (_num_AC > _res_AC) {
                    _AC = (unsigned int *) realloc(_AC, _num_AC * sizeof(unsigned int));
                    if (!_AC) throw EGGMEM;
                    _res_AC = _num_AC;
                }

                for (unsigned int i=0; i<_num_alternate; i++) {
                    _AC[i] = ((vcf::TypeInfo<int> *) next)->item(i);
                }

                if (_has_AN) {
                    _has_AC = true;
                }

            }

            if (!strcmp(next->get_ID(), "AF") && type==vcf::Info::Float && number==vcf::NUM_ALTERNATE) {
                _has_AF = true;

                _num_AF = _num_alternate;
                if (_num_AF > _res_AF) {
                    _AF = (double *) realloc(_AF, _num_AF * sizeof(double));
                    if (!_AF) throw EGGMEM;
                    _res_AF = _num_AF;
                }
                
                for (unsigned int i=0; i<_num_alternate; i++) {
                    _AF[i] = ((vcf::TypeInfo<double> *) next)->item(i);
                }
            }

            if (!strcmp(next->get_ID(), "AA") && type==vcf::Info::String && number==1) {

                _has_AA = true;
                _AA_string = ((vcf::StringInfo *) next)->item(0);

                // if AA is missing
                if (strlen(_AA_string)==1 && (_AA_string[0]=='.' || _AA_string[0]=='N' || _AA_string[0]=='n')) {
                    _AA_index = UNKNOWN;
                    _AA_string = _AA_missing;
                }

                // if AA non-missing, get its index
                else {
                    _AA_index = UNKNOWN;

                    // uppercase the AA string
                    unsigned int n = strlen(_AA_string);
                    if (_AA_string[0] != '<' && (!(_AA_string[0] == '[' || _AA_string[0] == ']' || _AA_string[n-1] == '[' || _AA_string[n-1] == ']'))) {
                        unsigned int d = 'a' - 'A';
                        for (unsigned int i=0; i<n; i++) {
                            if (_AA_string[i] >= 'a' && _AA_string[i] <='z') ((vcf::StringInfo *) next)->change(0, i, _AA_string[i] - d);
                        }
                    }

                    // find its index in list of alleles
                    _AA_index = UNKNOWN;
                    if (!strcmp(_AA_string, _reference)) _AA_index = 0;
                    else {
                        for (unsigned int i=0; i<_num_alternate; i++) {
                            if (!strcmp(_AA_string, _alternate[i])) {
                                _AA_index = i + 1;
                                break;
                            }
                        }
                    }
                    if (_AA_index == UNKNOWN) _AA_index = _num_alternate + 1;
                }
            }

            // finish
            if (curr_ch == '\t') break;
            if (curr_ch == '\n') break;
            if (_stream->eof()) break;
        }

        // check that FORMAT is present if >0 samples, and is not otherwise
        if (_num_samples == 0) {
            if (curr_ch != '\n' && _stream->eof()==false) throw EggFormatError(_fname, _currline+1, "VCF", "expect end of line (and no FORMAT string)", curr_ch);
            _currline++;
            _stream->peek();
            return;
            
        }

        if (curr_ch == '\n' || _stream->eof())  throw EggFormatError(_fname, _currline+1, "VCF", "expect FORMAT string (found end of line)");

        // read FORMAT specification
        unsigned int rankI = 0;
        unsigned int rankF = 0;
        unsigned int rankC = 0;
        unsigned int rankS = 0;        

        // default values (overrided if GT or PL is found)
        _ploidy = 2;
        _num_genotypes = (_num_alternate + 1) * (_num_alternate + 2) / 2;
        _has_GT = false;
        unsigned int PL_index = UNKNOWN;
        unsigned int GL_index = UNKNOWN;

        while (true) {
            // preallocate new format specification
            _num_formatEntries++;

            if (_num_formatEntries > _res_formatEntries) {
                
                _formatEntries = (vcf::Format **) realloc(_formatEntries, _num_formatEntries * sizeof(vcf::Format *));
                if (!_formatEntries) throw EGGMEM;
                // don't set it (has to be done)
                
                _formatRank = (unsigned int *) realloc(_formatRank, _num_formatEntries * sizeof(unsigned int));
                if (!_formatRank) throw EGGMEM;
                // don't set it (has to be done)
                
                _res_formatEntries = _num_formatEntries;
            }

            // get the format ID
            get_string(_buffer, _res_buffer, &VcfParser::check_sign, &VcfParser::stop_tabcolon, false);

            // identifies the Format concerned
            _formatEntries[_num_formatEntries-1] = find_format(_buffer);
            if (_formatEntries[_num_formatEntries-1] == NULL) EggFormatError(_fname, _currline+1, "VCF", "undefined FORMAT type: ", '\0', _buffer);            

            // record where GT can be found (only CANONICAL GT!)
            if (!strcmp(_formatEntries[_num_formatEntries-1]->get_ID(), "GT")
                && _formatEntries[_num_formatEntries-1]->get_number() == 1
                && _formatEntries[_num_formatEntries-1]->get_type() == vcf::Info::String) {
                    if (_num_formatEntries != 1) throw EggFormatError(_fname, _currline+1, "VCF", "GT must be used as first FORMAT entry");
                    _has_GT = true;
            }

            // record where PL can be found
            if (!strcmp(_formatEntries[_num_formatEntries-1]->get_ID(), "PL")
                && _formatEntries[_num_formatEntries-1]->get_type() == vcf::Info::Integer ) {
                    _has_PL = true;
                    PL_index = _num_formatEntries - 1;
            }

            // record where GL can be found
            if (!strcmp(_formatEntries[_num_formatEntries-1]->get_ID(), "GL")
                && _formatEntries[_num_formatEntries-1]->get_type() == vcf::Info::Float ) {
                    _has_GL = true;
                    GL_index = _num_formatEntries - 1;
            }


            // get format rank (within its type)
            switch (_formatEntries[_num_formatEntries-1]->get_type()) {
                
                case vcf::Info::Integer:
                    _formatRank[_num_formatEntries-1] = rankI++;
                    break;

                case vcf::Info::Float:
                    _formatRank[_num_formatEntries-1] = rankF++;
                    break;

                case vcf::Info::Character:
                    _formatRank[_num_formatEntries-1] = rankC++;
                    break;

                case vcf::Info::String:
                    _formatRank[_num_formatEntries-1] = rankS++;
                    break;
                
                case vcf::Info::Flag:
                    throw EggRuntimeError("unexpected error in egglib::VcfParser::read() [2]");
            }
            
            if (curr_ch=='\t') break;
        }

        // flag to set _has_GT after all samples are processed
        bool will_have_GT = false;

        // read all samples
        for (unsigned int i=0; i<_num_samples; i++) {
            // reset the sample info receptacle
             _sampleInfo[i]->reset();
             
            // process all entries
            unsigned int j;
            for (j=0; j<_num_formatEntries; j++) {

                // add an entry of the appropriate type
                if (_formatEntries[j] == NULL || _formatEntries[j]->get_type() == vcf::Info::String) {
                    _sampleInfo[i]->addStringEntry();
                }
                else {
                    switch (_formatEntries[j]->get_type()) {
                        case vcf::Info::Character:
                            _sampleInfo[i]->addCharacterEntry();
                            break;

                        case vcf::Info::Integer:
                            _sampleInfo[i]->addIntegerEntry();
                            break;

                        case vcf::Info::Float:
                            _sampleInfo[i]->addFloatEntry();
                            break;

                        case vcf::Info::String:
                        case vcf::Info::Flag:
                            throw EggRuntimeError("unexpected error in egglib::VcfParser::read() [3]");
                    }
                }

                // read data item(s)
                unsigned int n = 0;
                unsigned int num_missing = 0;

                while (true) { // read sample-specific values (matching FORMAT entries)

                    // detect an item of the appropriate type (and manage missing data)
                    switch (_formatEntries[j]->get_type()) {
                        case vcf::Info::String:
                            _sampleInfo[i]->addStringItem();
                            get_string(_sampleInfo[i]->_StringItems[_sampleInfo[i]->_num_StringEntries-1][n], _sampleInfo[i]->_res_len_StringItems[_sampleInfo[i]->_num_StringEntries-1][n], &VcfParser::check_sign, &VcfParser::stop_tabcoloncommalineEOF, true);
                            if (_sampleInfo[i]->_StringItems[_sampleInfo[i]->_num_StringEntries-1][n][0] == MAXCHAR) num_missing++;
                            break;
                        
                        case vcf::Info::Integer:
                            _sampleInfo[i]->addIntegerItem();
                            get_string(_buffer, _res_buffer, &VcfParser::check_integer, &VcfParser::stop_tabcoloncommalineEOF, true);

                            if (_buffer[0] == MAXCHAR) {
                                num_missing++;
                                _sampleInfo[i]->_IntegerItems[_sampleInfo[i]->_num_IntegerEntries-1][n] = MISSINGDATA;
                            }
                            else {
                                _sampleInfo[i]->_IntegerItems[_sampleInfo[i]->_num_IntegerEntries-1][n] = atoi(_buffer);
                            }
                            break;

                        case vcf::Info::Float:
                            _sampleInfo[i]->addFloatItem();
                            get_string(_buffer, _res_buffer, &VcfParser::check_float, &VcfParser::stop_tabcoloncommalineEOF, true);
                            if (_buffer[0] == MAXCHAR) {
                                num_missing++;
                                _sampleInfo[i]->_FloatItems[_sampleInfo[i]->_num_FloatEntries-1][n] = UNDEF;
                            }
                            else {
                                _sampleInfo[i]->_FloatItems[_sampleInfo[i]->_num_FloatEntries-1][n] = atof(_buffer);
                            }
                            break;

                        case vcf::Info::Character:
                            _sampleInfo[i]->addCharacterItem();
                            if (get_string(_buffer, _res_buffer, &VcfParser::check_float, &VcfParser::stop_tabcoloncommalineEOF, true) != 1) throw EggFormatError(_fname, _currline+1, "VCF", "expect a single character instead of ", '\0', _buffer);
                            if (_buffer[0] == MAXCHAR) num_missing++;
                            _sampleInfo[i]->_CharacterItems[_sampleInfo[i]->_num_CharacterEntries-1][n] = _buffer[0];
                            break;

                        case vcf::Info::Flag:
                            throw EggRuntimeError("unexpected error in egglib::VcfParser::read() [4]");
                    }
 
                    n++;
                    if (curr_ch == ':' || curr_ch == '\t' || curr_ch == '\n' || _stream->eof()) break;
                }

                // if only one missing data, set the number of items to 0
                if (num_missing == 1 && n == 1) {
                    switch (_formatEntries[j]->get_type()) {
                        case vcf::Info::String:
                            _sampleInfo[i]->_num_StringItems[_sampleInfo[i]->_num_StringEntries-1]--;
                            break;
                        case vcf::Info::Integer:
                            _sampleInfo[i]->_num_IntegerItems[_sampleInfo[i]->_num_IntegerEntries-1]--;
                            break;
                        case vcf::Info::Float:
                            _sampleInfo[i]->_num_FloatItems[_sampleInfo[i]->_num_FloatEntries-1]--;
                            break;
                        case vcf::Info::Character:
                            _sampleInfo[i]->_num_CharacterItems[_sampleInfo[i]->_num_CharacterEntries-1]--;
                            break;
                        case vcf::Info::Flag:
                            throw EggRuntimeError("unexpected error in egglib::VcfParser::read() [4]");
                    }
                }

                else {
                    // ensure the number of items is correct
                    unsigned int number = _formatEntries[j]->get_number();
                    if (number == vcf::NUM_ALTERNATE) number = _num_alternate;
                    if (number == vcf::NUM_POSSIBLE_ALLELES) number = _num_alternate + 1;
                    if (number == vcf::NUM_GENOTYPES) number = _num_genotypes;
                    if (number != UNKNOWN && n != number) throw EggFormatError(_fname, _currline+1, "VCF", "incorrect number of arguments for sample: ", '\0', _samples[i]);
                }

                // process GT string
                if (_has_GT && j == 0) { // if canonical GT was found, it must be at index 0

                    unsigned int n = strlen(_sampleInfo[i]->StringItem(0, 0));

                    // copy GT in buffer
                    if (n > _res_buffer) { // technically it is probably impossible
                        _buffer = (char *) realloc(_buffer, n * sizeof(char));
                        if (!_buffer) throw EGGMEM;
                        _res_buffer = n;
                    }

                    strcpy(_buffer, _sampleInfo[i]->StringItem(0, 0));

                    // detect ploidy
                    unsigned int nb = 0;
                    unsigned int ns = 0;

                    if (_buffer[0] != MAXCHAR) {
                        for (unsigned int k=0; k<n; k++) {
                            if (_buffer[k] == '|') { nb++; continue; }
                            if (_buffer[k] == '/') { ns++; continue; }
                            if (_buffer[k] == '.') { continue; }
                            if (_buffer[k] < '0' || _buffer[k] > '9') throw EggFormatError(_fname, _currline+1, "VCF", "invalid GT field: ", '\0', _buffer);
                        }
                    }

                    if (nb && ns) throw EggFormatError(_fname, _currline+1, "VCF", "invalid GT field (mixing phased/unphased codes): ", '\0', _buffer);

                    if (nb) {
                        _GT_phased[i] = true;
                    }
                    else {
                        _GT_phased[i] = false;
                        _GT_all_phased = false;
                    }

                    _ploidy = nb + ns + 1;
                    if (_ploidy == 0) throw EggFormatError(_fname, _currline+1, "VCF", "invalid genotype (GT) specification: ploidy is 0"); //don't believe that this can actually occur
                    for (unsigned int k=0; k<_num_samples; k++) {
                        if (_ploidy > _res_GT[k]) {
                            _GT[k] = (unsigned int *) realloc(_GT[k], _ploidy * sizeof(unsigned int));
                            if (!_GT[k]) throw EGGMEM;
                            _res_GT[k] = _ploidy;
                        }
                    }

                    // get all values
                    unsigned int start = 0;
                    unsigned int stop;
                    for (unsigned int k=0; k<_ploidy; k++) {
                        if (_buffer[0] == MAXCHAR) {
                            _GT[i][k] = UNKNOWN;
                            continue;
                        }
                        stop = start + 1;
                        while (true) {
                            if (_buffer[stop] == '/' || _buffer[stop] == '|') { _buffer[stop] = '\0'; break; }
                            if (_buffer[stop] == '\0') break;   
                            stop++;
                            if (stop == n) throw EggFormatError(_fname, _currline+1, "VCF", "invalid GT field: ", '\0', _sampleInfo[i]->StringItem(0, 0));
                        }
                        
                        if (!strcmp(_buffer + start, ".")) {
                            _GT[i][k] = UNKNOWN;
                        }
                        else {
                            _GT[i][k] = atoi(_buffer + start);
                            if (_GT[i][k] > _num_alternate /* counting ref */) throw EggFormatError(_fname, _currline+1, "VCF", "invalid GT field: ", '\0', _sampleInfo[i]->StringItem(0, 0));
                        }
                        start = stop + 1;
                    }

                    // compute number of genotypes
                    double x = 1.0;
                    unsigned int na = _num_alternate+1;
                    for (unsigned int i=1; i<na; i++) x *= (_ploidy + na - i) / static_cast<double>(i);
                    _num_genotypes = static_cast<unsigned int>(x);
                }

                // process PL string
                if (_has_PL && j == PL_index) {

                    if (_num_genotypes > _res_PL[i]) {
                        _PL[i] = (unsigned int *) realloc(_PL[i], _num_genotypes * sizeof(unsigned int));
                        if (!_PL[i]) throw EGGMEM;   
                        _res_PL[i] = _num_genotypes;
                    }


                    // check consistent ploidy
                    if (_sampleInfo[i]->num_IntegerItems(field_rank("PL")) != _num_genotypes) {
                        throw EggFormatError(_fname, _currline+1, "VCF", "invalid PL field, inconsistent number of genotypes in PL: ", '\0', _buffer);
                    }

                    // import PL values
                    for (unsigned int j=0; j<_num_genotypes; j++){
                        _PL[i][j] = _sampleInfo[i]->IntegerItem(field_rank("PL"), j);
                    }
                }
  

                // process gL string
                if (_has_GL && j == GL_index) {
                    if (_num_genotypes > _res_GL[i]) {
                        _GL[i] = (double *) realloc(_GL[i], _num_genotypes * sizeof(double));
                        if (!_GL[i]) throw EGGMEM;   
                        _res_GL[i] = _num_genotypes;

                    }

                    // check consistent ploidy
                    if (_sampleInfo[i]->num_FloatItems(field_rank("GL")) != _num_genotypes) {
                        throw EggFormatError(_fname, _currline+1, "VCF", "invalid GL field, inconsistent number of genotypes in GL: ", '\0', _buffer);
                    }

                    // import GL values
                    for (unsigned int j=0; j<_num_genotypes; j++){
                        _GL[i][j] = _sampleInfo[i]->FloatItem(field_rank("GL"), j);
                    }
                }


                // if GT is missing, attempt to convert from PL
                if (!_has_GT && _has_PL && _threshold_PL != UNKNOWN && _threshold_GL == UNKNOWN  && !strcmp(_formatEntries[j]->get_ID(), "PL")) {
                    if (_sampleInfo[i]->num_IntegerItems(field_rank("PL")) != _num_genotypes) throw EggFormatError(_fname, _currline+1, "VCF", "invalid PL field: invalid number of genotypes -- ", '\0', _buffer);

                    // set arbitrarily flags
                    _GT_phased[i] = false;
                    _GT_all_phased = false;
                    will_have_GT = true;

                    // allocate the GT table
                    if (_ploidy > _res_GT[i]) {
                        _GT[i] = (unsigned int *) realloc(_GT[i], _ploidy * sizeof(unsigned int));
                        if (!_GT[i]) throw EGGMEM;
                        _res_GT[i] = _ploidy;
                    }

                    // find best (smallest) PL value (ex aequo are unimportant)
                    unsigned int smallest=MAX;
                    unsigned int smallest_idx=UNKNOWN;
                    unsigned int second_smallest=MAX;
                    for (unsigned int k=0; k<_num_genotypes; k++) {
                        if (_PL[i][k] < smallest) {
                            smallest = _PL[i][k];
                            smallest_idx = k;
                        }
                    }

                    // find the second best PL = best while ignoring the best index
                    // if there is an ex aequo, it will be the same value
                    for (unsigned int k=0; k<_num_genotypes; k++) {
                        if (k!=smallest_idx && _PL[i][k] < second_smallest) {
                            second_smallest = _PL[i][k];
                        }
                    }

                    if (second_smallest - smallest >= _threshold_PL) {
                        // find the alleles matching the best genotype
                        _find_alleles(_num_alternate+1, _GT[i], smallest_idx);
                    }
                    else {
                        for (unsigned int k=0; k<_ploidy; k++) {
                            _GT[i][k] = UNKNOWN;
                        }
                    }
                }


                // if GT is missing, attempt to convert from GL
                if (!_has_GT && _has_GL && _threshold_PL == UNKNOWN && _threshold_GL != UNKNOWN  && !strcmp(_formatEntries[j]->get_ID(), "GL")) {
                   if (_sampleInfo[i]->num_FloatItems(field_rank("GL")) != _num_genotypes) throw EggFormatError(_fname, _currline+1, "VCF", "invalid GL field: invalid number of genotypes -- ", '\0', _buffer); 

                   // set arbitrarily flags
                    _GT_phased[i] = false;
                    _GT_all_phased = false;
                    will_have_GT = true;

                    // allocate the GT table

                    if (_ploidy > _res_GT[i]) {
                        _GT[i] = (unsigned int *) realloc(_GT[i], _ploidy * sizeof(unsigned int));
                        if (!_GT[i]) throw EGGMEM;
                        _res_GT[i] = _ploidy;
                    }

                    // find best (smallest) GL value (ex aequo are unimportant)
                    unsigned int smallest=MAX;
                    unsigned int smallest_idx=UNKNOWN;
                    unsigned int second_smallest=MAX;
                    //unsigned int tmp_PL;
                    double tmp_PL;
                    double tmp_GL;
                    for (unsigned int k=0; k<_num_genotypes; k++) {
                        //tmp_PL= (unsigned int)round(-10*_GL[i][k]); // convert GL(log10(L)) to PL(-10log10(L)) 
                        tmp_PL = -10*_GL[i][k];
                        if (tmp_PL < smallest) {
                            smallest = tmp_PL;
                            tmp_GL = _GL[i][k];
                            smallest_idx = k;
                        }
                    }
                    

                    // find the second best GL = best while ignoring the best index
                    // if there is an ex aequo, it will be the same value
                    for (unsigned int k=0; k<_num_genotypes; k++) {
                        //tmp_PL=  (unsigned int)round(-10*_GL[i][k]);// convert GL(log10(L)) to PL(-10log10(L)) 
                        tmp_PL = -10*_GL[i][k];
                        if (k!=smallest_idx && tmp_PL < second_smallest) {
                            second_smallest = tmp_PL;
                        }
                    }

                    if (second_smallest - smallest >= _threshold_GL) {
                        // find the alleles matching the best genotype
                        _find_alleles(_num_alternate+1, _GT[i], smallest_idx);
                    }
                    else {
                        for (unsigned int k=0; k<_ploidy; k++) {
                            _GT[i][k] = UNKNOWN;
                        }
                    }
                        
                }


                // finish
                if (curr_ch == '\t' || curr_ch == '\n' || _stream->eof()) break;  // as a result of this break
                                                                // the stop condition is never used
                                                                // (unless there are two many fields)
            }

            // if entries are missing, add them anyway except if GT
            for (j++; j<_num_formatEntries; j++) {  // the j++ start condition is there because we skip the last j++ of the previous loop

                if (!strcmp(_formatEntries[j]->get_ID(), "GT")) throw EggFormatError(_fname, _currline+1, "VCF", "it is not allowed to skip GT FORMAT field - for sample: ", '\0', _samples[i]);

                if (_formatEntries[j] == NULL || _formatEntries[j]->get_type() == vcf::Info::String) {
                    _sampleInfo[i]->addStringEntry();
                }
                else switch (_formatEntries[j]->get_type()) {

                    case vcf::Info::Character:
                        _sampleInfo[i]->addCharacterEntry();
                        break;
                    
                    case vcf::Info::Integer:
                        _sampleInfo[i]->addIntegerEntry();
                        break;
                        
                    case vcf::Info::Float:
                        _sampleInfo[i]->addFloatEntry();
                        break;
                        
                    case vcf::Info::String:
                    case vcf::Info::Flag:
                        throw EggRuntimeError("unexpected error in egglib::VcfParser::read() [5]");
                }
            }
            if (curr_ch != '\n' && curr_ch != '\t' && _stream->eof()==false) throw EggFormatError(_fname, _currline+1, "VCF", "0:expect a tab or an end of line there", curr_ch);
        }

        if (curr_ch != '\n' && _stream->eof()==false) throw EggFormatError(_fname, _currline+1, "VCF", "1:expect an end of line there", curr_ch);

        // conclude
        _has_GT |= will_have_GT;
        _currline++;
        _has_data = true;
    }

    void VcfParser::read_chromosome() {
        _previous_next=get_index();
        reset_variant();

        // checking
        if (_stream == NULL || (!_stream->good())) throw EggFormatError(_fname, _currline+1, "VCF", "error in file or invalid stream");

        // get chromosome
        get_string(_chrom, _res_chrom, &VcfParser::check_sign, &VcfParser::stop_tab, false);

        // get position
        get_string(_buffer, _res_buffer, &VcfParser::check_integer, &VcfParser::stop_tab, true);
        if (_buffer[0] == MAXCHAR) throw EggFormatError(_fname, _currline+1, "VCF", "position cannot be missing");

        _position = atol(_buffer);
        if (_position == 0) _position = BEFORE;
        else _position--;

        _stream->ignore(numeric_limits<streamsize>::max(), '\n'); // go to end of line
        _stream->peek();
    }

    bool VcfParser::has_data() const {
        return _has_data;
    }
     
    bool VcfParser::good() const {
        return (_stream != NULL && _stream->good() && _stream->tellg() != _idx_end_file);
    }

    const char * VcfParser::file_format() const {
        return _ff;
    }

    void VcfParser::init() {
        curr_ch = 0;
        prev_ch = 0;
        _stream = NULL;
        _idx_end_file = 0L;
        _previous_next=0;
        _idx_frt_sample=0;
        _res_fname = 1;
        _fname = (char *) malloc(1 * sizeof(char));
        if (!_fname) throw EGGMEM;
        _fname[0] = '\0';
        _currline = 0;
        _index= (char *) malloc(1 * sizeof(char));;
        _has_index=false;
        _res_buffer = 0;
        _buffer = NULL;
        _buffer_float = NULL;
        _res_buffer2 = 0;
        _buffer2 = NULL;

        _res_ff = 1;
        _ff = (char *) malloc(1 * sizeof(char));
        if (!_ff) throw EGGMEM;
        _ff[0] = '\0';

        _num_samples = 0;
        _res_samples = 0;
        _res_len_samples = NULL;
        _samples = NULL;
        _sampleInfo = NULL;
        _outgroup_map = NULL;
        _outgroup_AA = false;
        
        _num_filter = 0;
        _num_info = 0;
        _num_format = 0;
        _num_meta = 0;
        _num_alt = 0;
        _res_filter = 0;
        _res_info = 0;
        _res_format = 0;
        _res_meta = 0;
        _res_alt = 0;
        _filter = NULL;
        _info = NULL;
        _format = NULL;
        _meta = NULL;
        _alt = NULL;
        
        _res_chrom = 1;
        _chrom = (char *) malloc(1 * sizeof(char));
        if (!_chrom) throw EGGMEM;
        _chrom[0] = '\0';
        _position = UNKNOWN;
        _num_ID = 0;
        _res_ID = 0;
        _res_len_ID = NULL;
        _ID = NULL;
        _len_reference = 0;
        _res_reference = 1;
        _reference = (char *) malloc(1 * sizeof(char));
        if (!_reference) throw EGGMEM;
        _reference[0] = '\0';
        _num_alternate = 0;
        _res_alternate = 0;
        _type_alternate = NULL;
        _res_len_alternate = NULL;
        _alternate = NULL;
        _num_failed_test = UNKNOWN;
        _res_failed_test = 0;
        _res_len_failed_test = NULL;
        _failed_test = NULL;
        
        _num_FlagInfo = 0;
        _res_FlagInfo = 0;
        _FlagInfo = NULL;
        _num_CharacterInfo = 0;
        _res_CharacterInfo = 0;
        _num_IntegerInfo = 0;
        _res_IntegerInfo = 0;
        _IntegerInfo = NULL;
        _CharacterInfo = NULL;
        _num_FloatInfo = 0;
        _res_FloatInfo = 0;
        _FloatInfo = NULL;
        _num_StringInfo = 0;
        _res_StringInfo = 0;
        _StringInfo = NULL;

        _num_formatEntries = 0;
        _res_formatEntries = 0;
        _formatEntries = NULL;
        _formatRank = NULL;

        _has_data = false;
        _has_AN = false;
        _AN = 0;
        _has_AC = false;
        _has_AC_ss = false;
        _num_AC = 0;
        _res_AC = 0;
         _AC = NULL;
        _has_AF = false;
        _num_AF = 0;
        _res_AF = 0;
        _AF = NULL;
        _has_AA = false;
        _AA_index = UNKNOWN;
        _AA_string = NULL;
        _AA_missing = (char *) malloc(2 * sizeof(char));
        if (!_AA_missing) throw EGGMEM;
        strcpy(_AA_missing, "?");

        _has_GT = false;
        _res_GT = NULL;
        _res_PL = NULL;
        _res_GL = NULL;
        _GT = NULL;
        _GT_phased = NULL;
        _GT_all_phased = false;
        _has_PL = false;
        _PL = NULL;
        _has_GL = false;
        _GL = NULL;
        _num_genotypes = UNKNOWN;
        _ploidy = UNKNOWN;

        _genotype_idx_helper_size = 2;
        _genotype_idx_helper = (unsigned int **) malloc(_genotype_idx_helper_size * sizeof(unsigned int *));
        if (!_genotype_idx_helper) throw EGGMEM;
        _genotype_idx_helper[0] = NULL;
        for (unsigned int i=1; i<_genotype_idx_helper_size; i++) {
            unsigned int ng = i*(i+1)/2;
            _genotype_idx_helper[i] = (unsigned int *) malloc(2 * ng * sizeof(unsigned int));
            if (!_genotype_idx_helper[i]) throw EGGMEM;
            for (unsigned int a=0; a<i; a++) {
                for (unsigned int b=a; b<i; b++) {
                    unsigned int idx = b*(b+1)/2 + a;
                    _genotype_idx_helper[i][2*idx] = a;
                    _genotype_idx_helper[i][2*idx+1] = b;
                }
            }
        }
    }

    void VcfParser::free() {

         if (_fname) ::free(_fname);
         if (_index) ::free(_index);
         if (_buffer) ::free(_buffer);
         if(_buffer_float) ::free(_buffer_float);
         if (_ff) ::free(_ff);

         for (unsigned int i=0; i<_res_filter; i++) {
             if (_filter[i]) delete (_filter[i]);
         }
         if (_filter) ::free(_filter);

         for (unsigned int i=0; i<_res_info; i++) {
             if (_info[i]) delete (_info[i]);
         }
         if (_info) ::free(_info);

         for (unsigned int i=0; i<_res_format; i++) {
             if (_format[i]) delete (_format[i]);
         }
         if (_format) ::free(_format);

         for (unsigned int i=0; i<_res_meta; i++) {
             if (_meta[i]) delete (_meta[i]);
         }
         if (_meta) ::free(_meta);
         
         for (unsigned int i=0; i<_res_alt; i++) {
             if (_alt[i]) delete (_alt[i]);
         }
         if (_alt) ::free(_alt);

         for (unsigned int i=0; i<_res_samples; i++) {
             if (_samples[i]) ::free(_samples[i]);
             if (_sampleInfo[i]) delete _sampleInfo[i];
         }
         if (_res_len_samples) ::free(_res_len_samples);
         if (_samples) ::free(_samples);
         if (_sampleInfo) ::free(_sampleInfo);
         if (_outgroup_map) ::free(_outgroup_map);
         
         if (_chrom) ::free(_chrom);
         for (unsigned int i=0; i<_res_ID; i++) {
             if (_ID[i]) ::free(_ID[i]);
         }
         if (_res_len_ID) ::free(_res_len_ID);
         if (_ID) ::free(_ID);

         if (_reference) ::free(_reference);

         for (unsigned int i=0; i<_res_alternate; i++) {
             if (_alternate[i]) ::free(_alternate[i]);
         }
         if (_alternate) ::free(_alternate);
         if (_type_alternate) ::free(_type_alternate);
         if (_res_len_alternate) ::free(_res_len_alternate);

         for (unsigned int i=0; i<_res_failed_test; i++) {  // DO NOT USE _num_failed_test, it can be egglib::UNKNOWN!!
             if (_failed_test[i]) ::free(_failed_test[i]);
         }
         if (_failed_test) ::free(_failed_test);
         if (_res_len_failed_test) ::free(_res_len_failed_test);

        for (unsigned int i=0; i<_res_FlagInfo; i++) {
            if (_FlagInfo[i]) delete _FlagInfo[i];
        }
        if (_FlagInfo) ::free(_FlagInfo);

        for (unsigned int i=0; i<_res_CharacterInfo; i++) {
            if (_CharacterInfo[i]) delete _CharacterInfo[i];
        }
        if (_CharacterInfo) ::free(_CharacterInfo);

        for (unsigned int i=0; i<_res_IntegerInfo; i++) {
            if (_IntegerInfo[i]) delete _IntegerInfo[i];
        }
        if (_IntegerInfo) ::free(_IntegerInfo);

        for (unsigned int i=0; i<_res_FloatInfo; i++) {
            if (_FloatInfo[i]) delete _FloatInfo[i];
        }
        if (_FloatInfo) ::free(_FloatInfo);

        for (unsigned int i=0; i<_res_StringInfo; i++) {
            if (_StringInfo[i]) delete _StringInfo[i];
        }
        if (_StringInfo) ::free(_StringInfo);

        if (_formatEntries) ::free(_formatEntries);
        if (_formatRank) ::free(_formatRank);
        
        if (_AC) ::free(_AC);
        if (_AF) ::free(_AF);
        
        if (_GT_phased) ::free(_GT_phased);
        for (unsigned int i=0; i<_res_samples; i++) {
            if (_GT[i]) ::free(_GT[i]);
        }
        if (_GT) ::free(_GT);

        for (unsigned int i=0; i<_res_samples; i++) {
           if (_PL[i]) ::free(_PL[i]);
        }
        if (_PL) ::free(_PL);

        for (unsigned int i=0; i<_res_samples; i++) {
           if (_GL[i]) ::free(_GL[i]);
        }
        if (_GL) ::free(_GL);

        if (_res_GT) ::free(_res_GT);
        if (_res_PL) ::free(_res_PL);
        if (_res_GL) ::free(_res_GL);
        if (_AA_missing) ::free(_AA_missing);

        for (unsigned int i=0; i<_genotype_idx_helper_size; i++) {
            if (_genotype_idx_helper[i]) ::free(_genotype_idx_helper[i]);
        }
        if (_genotype_idx_helper) ::free(_genotype_idx_helper);
    }

    void VcfParser::get_4_fields(bool format) {
        next();
        if (curr_ch != '<') throw EggFormatError(_fname, _currline+1, "VCF", "first character of meta-information specification is not \"<\" as expected", curr_ch);

        unsigned int flag = 0;  // should be 1 + 2 + 4 + 8 = 15
        unsigned int nterms = 0; // should 4
        vcf::Info info;

        while (true) {

            if (curr_ch == '>') {
                next();
                if (curr_ch == '\r') next();
                if (curr_ch != '\n') throw EggFormatError(_fname, _currline+1, "VCF", "expect end of line after INFO or FORMAT specification", curr_ch);
                _currline++;
                break;
            }
                    
            get_string(_buffer, _res_buffer, &VcfParser::check_letter, &VcfParser::stop_equal, false);
            
            if (!strcmp(_buffer, "ID")) {
                get_string(_buffer, _res_buffer, &VcfParser::check_alphanumericunderscore, &VcfParser::stop_field, false);
                if (_buffer[0] >= '0' && _buffer[0] <= '9') throw EggFormatError(_fname, _currline+1, "VCF", "invalid INFO or FORMAT ID (cannot start by a number)");
                info.set_ID(_buffer);
                flag |= 1;
                nterms++;
                continue;
            }

            if (!strcmp(_buffer, "Number")) {
                switch (_stream->peek()) {
                    case 'A':
                        info.set_number(vcf::NUM_ALTERNATE);
                        next();
                        next();
                        if (curr_ch != ',' && curr_ch != '>') throw EggFormatError(_fname, _currline+1, "VCF", "expect a comma after `Number` specification", curr_ch);
                        break;
                    case 'R':
                        info.set_number(vcf::NUM_POSSIBLE_ALLELES);
                        next();
                        next();
                        if (curr_ch != ',' && curr_ch != '>') throw EggFormatError(_fname, _currline+1, "VCF", "expect a comma after `Number` specification", curr_ch);
                        break;
                    case 'G':
                        info.set_number(vcf::NUM_GENOTYPES);
                        next();
                        next();
                        if (curr_ch != ',' && curr_ch != '>') throw EggFormatError(_fname, _currline+1, "VCF", "expect a comma after `Number` specification", curr_ch);
                        break;
                    case '.':
                        info.set_number(UNKNOWN);
                        next();
                        next();
                        if (curr_ch != ',' && curr_ch != '>') throw EggFormatError(_fname, _currline+1, "VCF", "expect a comma after `Number` specification", curr_ch);
                        break;
                    default:
                        get_string(_buffer, _res_buffer, &VcfParser::check_integer, &VcfParser::stop_field, false);
                        info.set_number(atoi(_buffer));
                }
                flag |= 2;
                nterms++;
                continue;
            }

            if (!strcmp(_buffer, "Type")) {
                get_string(_buffer, _res_buffer, &VcfParser::check_letter, &VcfParser::stop_field, false);

                if (!strcmp(_buffer, "Integer")) {
                    info.set_type(vcf::Info::Integer);
                    flag |= 4;
                    nterms++;
                    continue;
                }

                if (!strcmp(_buffer, "Float")) {
                    info.set_type(vcf::Info::Float);
                    flag |= 4;
                    nterms++;
                    continue;
                }

                if (!strcmp(_buffer, "String")) {
                    info.set_type(vcf::Info::String);
                    flag |= 4;
                    nterms++;
                    continue;
                }

                if (!strcmp(_buffer, "Character")) {
                    info.set_type(vcf::Info::Character);
                    flag |= 4;
                    nterms++;
                    continue;
                }

                if (!format && !strcmp(_buffer, "Flag")) {
                    info.set_type(vcf::Info::Flag);
                    flag |= 4;
                    nterms++;
                    continue;
                }

                throw EggFormatError(_fname, _currline+1, "VCF", "invalid Type for INFO for FORMAT specification");
            }

            if (!strcmp(_buffer, "Description")) {
                next();
                if (curr_ch != '"') throw  EggFormatError(_fname, _currline+1, "VCF", "invalid description string for INFO or FORMAT (expect a double quote `\"`)", curr_ch);
                get_string(_buffer, _res_buffer, &VcfParser::check_string, &VcfParser::stop_quote, false);
                info.set_description(_buffer);
                next();
                if (curr_ch != '>' && curr_ch != ',') throw  EggFormatError(_fname, _currline+1, "VCF", "malformed INFO or FORMAT specification: expect a field separator \",\" or \">\"", curr_ch);
                flag |= 8;
                nterms++;
                continue;
            }

            next();
            if (curr_ch != '"') throw  EggFormatError(_fname, _currline+1, "VCF", "invalid extra field string for INFO or FORMAT (expect a double quote `\"`)", curr_ch);
            get_string(_buffer2, _res_buffer2, &VcfParser::check_string, &VcfParser::stop_quote, false);
            info.set_extra(_buffer, _buffer2);
            next();
            if (curr_ch != '>' && curr_ch != ',') throw  EggFormatError(_fname, _currline+1, "VCF", "malformed INFO or FORMAT specification: expect a field separator \",\" or \">\"", curr_ch);
            continue;
        }
        
        // check consistency
        
        if (flag != 15) throw  EggFormatError(_fname, _currline+1, "VCF", "malformed INFO or FORMAT specification");
        if (nterms != 4) throw  EggFormatError(_fname, _currline+1, "VCF", "malformed INFO or FORMAT specification");
        if (info.get_type() == vcf::Info::Flag && info.get_number() != 0) throw  EggFormatError(_fname, _currline+1, "VCF", "malformed INFO specification: if Type is Flag, Number should be 0");
        if (info.get_type() != vcf::Info::Flag && info.get_number() == 0) throw  EggFormatError(_fname, _currline+1, "VCF", "malformed INFO specification: if Type is not Flag, Number should not be 0");
        
        // copy to store (ensure no duplicate)

        vcf::Info * dest;
        if (!format) {
            dest = find_info(info.get_ID());
            if (dest == NULL) dest = add_info();
        }
        else {
            dest = find_format(info.get_ID());
            if (dest == NULL) dest = (vcf::Info *) add_format();
        }
        dest->update(info.get_ID(), info.get_number(), info.get_type(), info.get_description());
    }

    void VcfParser::get_2_fields(bool alt) {
        next();
        if (curr_ch != '<') throw EggFormatError(_fname, _currline+1, "VCF", "first character of meta-information specification is not \"<\" as expected", curr_ch);

        unsigned int flag = 0;  // should be 1 + 2 = 3
        unsigned int nterms = 0; // should be 2
        vcf::Filter filter;
        
        while (true) {

            if (curr_ch == '>') {
                next();
                if (curr_ch == '\r') next();
                if (curr_ch != '\n') throw EggFormatError(_fname, _currline+1, "VCF", "expect end of line after FILTER or ALT specification", curr_ch);
                _currline++;
                break;
            }

            get_string(_buffer, _res_buffer, &VcfParser::check_letter, &VcfParser::stop_equal, false);
            
            if (!strcmp(_buffer, "ID")) {
                get_string(_buffer, _res_buffer, &VcfParser::check_alphanumericunderscore, &VcfParser::stop_field, false);
                if (_buffer[0] >= '0' && _buffer[0] <= '9') throw EggFormatError(_fname, _currline+1, "VCF", "invalid FILTER or ALT ID (cannot start by a number)");
                filter.set_ID(_buffer);
                flag |= 1;
                nterms++;
                continue;
            }

            if (!strcmp(_buffer, "Description")) {
                next();
                if (curr_ch != '"') throw  EggFormatError(_fname, _currline+1, "VCF", "invalid description string for FILTER or ALT (expect a double quote `\"`)", curr_ch);
                get_string(_buffer, _res_buffer, &VcfParser::check_string, &VcfParser::stop_quote, false);
                filter.set_description(_buffer);
                next();
                if (curr_ch != '>' && curr_ch != ',') throw  EggFormatError(_fname, _currline+1, "VCF", "malformed FILTER or ALT specification: expect a field separator \",\" or \">\"", curr_ch);
                flag |= 2;
                nterms++;
                continue;
            }

            next();
            if (curr_ch != '"') throw  EggFormatError(_fname, _currline+1, "VCF", "invalid extra field string for INFO or FORMAT (expect a double quote `\"`)", curr_ch);
            get_string(_buffer2, _res_buffer2, &VcfParser::check_string, &VcfParser::stop_quote, false);
            filter.set_extra(_buffer, _buffer2);
            next();
            if (curr_ch != '>' && curr_ch != ',') throw  EggFormatError(_fname, _currline+1, "VCF", "malformed INFO or FORMAT specification: expect a field separator \",\" or \">\"", curr_ch);
            continue;
        }

        // check consistency

        if (flag != 3) throw  EggFormatError(_fname, _currline, "VCF", "malformed FILTER or ALT specification");
        if (nterms != 2) throw  EggFormatError(_fname, _currline, "VCF", "malformed FILTER or ALT specification");

        // copy to storage (ensure no duplicate)

        vcf::Filter * dest;
        if (!alt) {
            dest = find_filter(filter.get_ID());
            if (dest == NULL) dest = add_filter();
        }
        else {
            dest = (vcf::Filter *) find_alt(filter.get_ID());
            if (dest == NULL) dest = add_alt();
        }
        dest->update(filter.get_ID(), filter.get_description());
    }

    unsigned int VcfParser::get_string(char *& where, unsigned int& _res_, bool (VcfParser::* check)(unsigned int), bool (VcfParser::* stop)(), bool catch_missing) {

        unsigned int n = 0;

        bool dot_escaped = false; // allow a single dot whatever check method

        while (true) {
            next();

            if ((this->*stop)()) break;
            if (!(this->*check)(n)) {
                if (curr_ch == '.' && dot_escaped == false) {
                    dot_escaped = true;
                }else{
                    throw EggFormatError(_fname, _currline+1, "VCF", "invalid character found", curr_ch);
                }
            }

            n++;
            if ((n+1) > _res_) {
                where = (char *) realloc(where, (n+1) * sizeof(char));
                if (!where) throw EGGMEM;
                _res_ = (n+1);
            }
            where[n-1] = curr_ch;
        }

        where[n] = '\0';
        if (n < 1) throw EggFormatError(_fname, _currline+1, "VCF", "empty field or specification here");

        if (catch_missing && !strcmp(where, ".")) where[0] = MAXCHAR;

        return n;
    }

    bool VcfParser::stop_equal() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, _currline+1, "VCF", "file truncated [code: 1]");
        if (curr_ch == '=') return true;
        if (curr_ch == '\r') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return");
        if (curr_ch == '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected end of line");
        return false;
    }

    bool VcfParser::stop_equalsemicolontablineEOF() {
        if (_stream->gcount()==0 && _stream->eof()) {
            return true;
        }
        if (curr_ch == '=') return true;
        if (curr_ch == ';') return true;
        if (curr_ch == '\t') return true;
        if (curr_ch == '\n') return true;
        if (curr_ch == '\r') {
            next();
            if (curr_ch != '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return (not followed by a new line)");
            return true;
        }
        return false;
    }

    bool VcfParser::stop_colon() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, _currline+1, "VCF", "file truncated [code: 2]");
        if (curr_ch == ':') return true;
        if (curr_ch == '\r') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return");
        if (curr_ch == '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected end of line");
        return false;
    }

    bool VcfParser::stop_tab() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, _currline+1, "VCF", "file truncated [code: 3]");
        if (curr_ch == '\t') return true;
        if (curr_ch == '\r') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return");
        if (curr_ch == '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected end of line");
        return false;
    }

    bool VcfParser::stop_tabsemicolon() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, _currline+1, "VCF", "file truncated [code: 4]");
        if (curr_ch == '\t') return true;
        if (curr_ch == ';') return true;
        if (curr_ch == '\r') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return");
        if (curr_ch == '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected end of line");
        return false;
    }

    bool VcfParser::stop_tabcolon() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, _currline+1, "VCF", "file truncated [code: 5]");
        if (curr_ch == '\t') return true;
        if (curr_ch == ':') return true;
        if (curr_ch == '\r') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return");
        if (curr_ch == '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected end of line");
        return false;
    }

    bool VcfParser::stop_tabcomma() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, _currline+1, "VCF", "file truncated [code: 6]");
        if (curr_ch == '\t') return true;
        if (curr_ch == ',') return true;
        if (curr_ch == '\r') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return");
        if (curr_ch == '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected end of line");
        return false;
    }

    bool VcfParser::stop_tabsemicoloncomma() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, _currline+1, "VCF", "file truncated [code: 7]");
        if (curr_ch == '\t') return true;
        if (curr_ch == ';') return true;
        if (curr_ch == ',') return true;
        if (curr_ch == '\r') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return");
        if (curr_ch == '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected end of line");
        return false;
    }

    bool VcfParser::stop_tabsemicoloncommalineEOF() {
        if (_stream->gcount()==0 && _stream->eof()) {
            return true;
        }
        if (curr_ch == '\t') return true;
        if (curr_ch == ';') return true;
        if (curr_ch == ',') return true;
        if (curr_ch == '\n') return true;
        if (curr_ch == '\r') {
            next();
            if (curr_ch != '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return (not followed by a new line)");
            return true;
        }
        return false;
    }

    bool VcfParser::stop_tabcoloncommalineEOF() {
        if (_stream->gcount()==0 && _stream->eof()) {
            return true;
        }
        if (curr_ch == '\t') return true;
        if (curr_ch == ':') return true;
        if (curr_ch == ',') return true;
        if (curr_ch == '\n') return true;
        if (curr_ch == '\r') {
            next();
            if (curr_ch != '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return (not followed by a new line)");
            return true;
        }
        return false;
    }

    bool VcfParser::stop_quote() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, _currline+1, "VCF", "file truncated [code: 8]");
        if (curr_ch == '\\' && _stream->peek() == '"') {
            _stream->get(curr_ch);
            return false;
        }
        if (curr_ch == '\\' && _stream->peek() == '\\') {
            _stream->get(curr_ch);
            return false;
        }
        if (curr_ch == '"') return true;
        if (curr_ch == '\r') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return");
        if (curr_ch == '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected end of line");
        return false;
    }

    bool VcfParser::stop_field() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, _currline+1, "VCF", "file truncated [code: 9]");
        if (curr_ch == ',') return true;
        if (curr_ch == '>') return true;
        if (curr_ch == '\r') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return");
        if (curr_ch == '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected end of line");
        return false;
    }

    bool VcfParser::stop_line() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, _currline+1, "VCF", "file truncated [code: 10s]");
        if (curr_ch == '\n') return true;
        if (curr_ch == '\r') {
            next();
            if (curr_ch != '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return (not followed by a new line)");
            return true;
        }
        return false;
    }

    bool VcfParser::stop_linetabEOF() {
        if (_stream->gcount()==0 && _stream->eof()) {
            return true;
        }
        if (curr_ch == '\t') return true;
        if (curr_ch == '\n') return true;
        if (curr_ch == '\r') {
            next();
            if (curr_ch != '\n') throw EggFormatError(_fname, _currline+1, "VCF", "unexpected carriage return (not followed by a new line)");
            return true;
        }
        return false;
    }

    bool VcfParser::check_letter(unsigned int idx) {
        if (curr_ch >= 'A' && curr_ch <= 'Z') return true;
        if (curr_ch >= 'a' && curr_ch <= 'z') return true;
        if (curr_ch == '_') return true;
        return false;
    }

    bool VcfParser::check_sign(unsigned int idx) {
        if (curr_ch < '!' || curr_ch > '~') return false;
        return true;
    }

    bool VcfParser::check_sign_and_space(unsigned int idx) {
        if (curr_ch < ' ' || curr_ch > '~') return false;
        return true;
    }

    bool VcfParser::check_string(unsigned int idx) {
        if (curr_ch == ' ' || curr_ch == '\t') return true;
        if (curr_ch < '!' || curr_ch > '~') return false;
        return true;
    }
    
    bool VcfParser::check_alphanumericunderscore(unsigned int idx) {
        if (curr_ch >= 'A' && curr_ch <= 'Z') return true;
        if (curr_ch >= 'a' && curr_ch <= 'z') return true;
        if (curr_ch >= '0' && curr_ch <= '9') return true;
        if (curr_ch == '*') return true;
        if (curr_ch == '_') return true;
        if (curr_ch == '-') return true;
        if (curr_ch == '+') return true;
        if (curr_ch == '.') return true;
        return false;
    }
    
    bool VcfParser::check_integer(unsigned int idx) {
        if (curr_ch >= '0' && curr_ch <= '9') return true;
        if (idx == 0 && (curr_ch == '-' || curr_ch == '+')) return true;
        return false;
    }
    
    bool VcfParser::check_float(unsigned int idx) {
        if (curr_ch >= '0' && curr_ch <= '9') return true;
        if (curr_ch == '.') return true;
        if (curr_ch == '-') return true; // - and + can appear after a "e" (even if not idx==0)
        if (curr_ch == '+') return true;
        if (curr_ch == 'E') return true;
        if (curr_ch == 'e') return true;
        return false;
    }
    
    bool VcfParser::check_bases(unsigned int idx) {

        switch (curr_ch) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'N':
                return true;
            case 'a':
            case 'c':
            case 'g':
            case 't':
            case 'n':
                curr_ch = toupper(curr_ch);
                return true;
            case '-':
                return _allow_gap;
            default:
                return false;
        }
    }

    void VcfParser::next() {
        prev_ch = curr_ch;
        _stream->get(curr_ch);
    }

    vcf::Filter * VcfParser::add_filter() {
        _num_filter++;
        if (_num_filter > _res_filter) {
            _filter = (vcf::Filter **) realloc(_filter, _num_filter * sizeof(vcf::Filter *));
            if (!_filter) throw EGGMEM;
            _filter[_num_filter-1] = new(std::nothrow) vcf::Filter;
            if (!_filter[_num_filter-1]) throw EGGMEM;
            _res_filter = _num_filter;
        }
        return _filter[_num_filter-1];
    }

    vcf::Info * VcfParser::add_info() {
        _num_info++;
        if (_num_info > _res_info) {
            _info = (vcf::Info **) realloc(_info, _num_info * sizeof(vcf::Info *));
            if (!_info) throw EGGMEM;
            _info[_num_info-1] = new(std::nothrow) vcf::Info;
            if (!_info[_num_info-1]) throw EGGMEM;
            _res_info = _num_info;
        }
        return _info[_num_info-1];
    }

    vcf::Format * VcfParser::add_format() {
        _num_format++;
        if (_num_format > _res_format) {
            _format = (vcf::Format **) realloc(_format, _num_format * sizeof(vcf::Format *));
            if (!_format) throw EGGMEM;
            _format[_num_format-1] = new(std::nothrow) vcf::Format;
            if (!_format[_num_format-1]) throw EGGMEM;
            _res_format = _num_format;
        }
        return _format[_num_format-1];
    }

    vcf::Meta * VcfParser::add_meta() {
        _num_meta++;
        if (_num_meta > _res_meta) {
            _meta = (vcf::Meta **) realloc(_meta, _num_meta * sizeof(vcf::Meta *));
            if (!_meta) throw EGGMEM;
            _meta[_num_meta-1] = new(std::nothrow) vcf::Meta;
            if (!_meta[_num_meta-1]) throw EGGMEM;
            _res_meta = _num_meta;
        }
        return _meta[_num_meta-1];
    }

    vcf::Alt * VcfParser::add_alt() {
        _num_alt++;
        if (_num_alt > _res_alt) {
            _alt = (vcf::Alt **) realloc(_alt, _num_alt * sizeof(vcf::Alt *));
            if (!_alt) throw EGGMEM;
            _alt[_num_alt-1] = new(std::nothrow) vcf::Alt;
            if (!_alt[_num_alt-1]) throw EGGMEM;
            _res_alt = _num_alt;
        }
        return _alt[_num_alt-1];
    }

    void VcfParser::add_filter(const char * id, const char * descr) {
        vcf::Filter * dest = find_filter(id);
        if (dest == NULL) dest = add_filter();
        dest->update(id, descr);
    }
        
    void VcfParser::add_alt(const char * id, const char * descr) {
        vcf::Alt * dest = find_alt(id);
        if (dest == NULL) dest = add_alt();
        dest->update(id, descr);
    }

    void VcfParser::add_info(const char * id, unsigned int num, vcf::Info::Type type, const char * descr) {
        vcf::Info * dest = find_info(id);
        if (dest == NULL) dest = add_info();
        dest->update(id, num, type, descr);
    }

    void VcfParser::add_format(const char * id, unsigned int num, vcf::Info::Type type, const char * descr) {
        vcf::Format * dest = find_format(id);
        if (dest == NULL) dest = add_format();
        dest->update(id, num, type, descr);
    }
    
    void VcfParser::add_meta(const char * key, const char * val) {
        vcf::Meta * dest = find_meta(key);
        if (dest == NULL) dest = add_meta();
        dest->update(key, val);
    }

    vcf::FlagInfo * VcfParser::add_FlagInfo(unsigned int expected_number) {
        _num_FlagInfo++;
        if (_num_FlagInfo > _res_FlagInfo) {
            _FlagInfo = (vcf::FlagInfo **) realloc(_FlagInfo, _num_FlagInfo * sizeof(vcf::FlagInfo *));
            if (!_FlagInfo) throw EGGMEM;
            _FlagInfo[_num_FlagInfo-1] = new(std::nothrow) vcf::FlagInfo;
            if (!_FlagInfo[_num_FlagInfo-1]) throw EGGMEM;
            _res_FlagInfo = _num_FlagInfo;
        }
        return _FlagInfo[_num_FlagInfo-1];
    }

    vcf::TypeInfo<char> * VcfParser::add_CharacterInfo(unsigned int expected_number) {
        _num_CharacterInfo++;
        if (_num_CharacterInfo > _res_CharacterInfo) {
            _CharacterInfo = (vcf::TypeInfo<char> **) realloc(_CharacterInfo, _num_CharacterInfo * sizeof(vcf::TypeInfo<char> *));
            if (!_CharacterInfo) throw EGGMEM;
            _CharacterInfo[_num_CharacterInfo-1] = new(std::nothrow) vcf::TypeInfo<char>;
            if (!_CharacterInfo[_num_CharacterInfo-1]) throw EGGMEM;
            _res_CharacterInfo = _num_CharacterInfo;
        }
        else {
            _CharacterInfo[_num_CharacterInfo-1]->reset();
        }
        _CharacterInfo[_num_CharacterInfo-1]->set_expected_number(expected_number);
        return _CharacterInfo[_num_CharacterInfo-1];
    }

    vcf::TypeInfo<int> * VcfParser::add_IntegerInfo(unsigned int expected_number) {
        _num_IntegerInfo++;
        if (_num_IntegerInfo > _res_IntegerInfo) {
            _IntegerInfo = (vcf::TypeInfo<int> **) realloc(_IntegerInfo, _num_IntegerInfo * sizeof(vcf::TypeInfo<int> *));
            if (!_IntegerInfo) throw EGGMEM;
            _IntegerInfo[_num_IntegerInfo-1] = new(std::nothrow) vcf::TypeInfo<int>;
            if (!_IntegerInfo[_num_IntegerInfo-1]) throw EGGMEM;
            _res_IntegerInfo = _num_IntegerInfo;
        }
        else {
            _IntegerInfo[_num_IntegerInfo-1]->reset();
        }
        _IntegerInfo[_num_IntegerInfo-1]->set_expected_number(expected_number);
        return _IntegerInfo[_num_IntegerInfo-1];
    }

    vcf::TypeInfo<double> * VcfParser::add_FloatInfo(unsigned int expected_number) {
        _num_FloatInfo++;
        if (_num_FloatInfo > _res_FloatInfo) {
            _FloatInfo = (vcf::TypeInfo<double> **) realloc(_FloatInfo, _num_FloatInfo * sizeof(vcf::TypeInfo<double> *));
            if (!_FloatInfo) throw EGGMEM;
            _FloatInfo[_num_FloatInfo-1] = new(std::nothrow) vcf::TypeInfo<double>;
            if (!_FloatInfo[_num_FloatInfo-1]) throw EGGMEM;
            _res_FloatInfo = _num_FloatInfo;
        }
        else {
            _FloatInfo[_num_FloatInfo-1]->reset();
        }
        _FloatInfo[_num_FloatInfo-1]->set_expected_number(expected_number);
        return _FloatInfo[_num_FloatInfo-1];
    }

    vcf::StringInfo * VcfParser::add_StringInfo(unsigned int expected_number) {
        _num_StringInfo++;
        if (_num_StringInfo > _res_StringInfo) {
            _StringInfo = (vcf::StringInfo **) realloc(_StringInfo, _num_StringInfo * sizeof(vcf::StringInfo *));
            if (!_StringInfo) throw EGGMEM;
            _StringInfo[_num_StringInfo-1] = new(std::nothrow) vcf::StringInfo;
            if (!_StringInfo[_num_StringInfo-1]) throw EGGMEM;
            _res_StringInfo = _num_StringInfo;
        }
        else {
            _StringInfo[_num_StringInfo-1]->reset();
        }
        _StringInfo[_num_StringInfo-1]->set_expected_number(expected_number);
        return _StringInfo[_num_StringInfo-1];
    }

    unsigned int VcfParser::num_filter() const {
        return _num_filter;
    }
    
    const vcf::Filter * VcfParser::get_filter(unsigned int i) const {
        return _filter[i];
    }

    unsigned int VcfParser::num_info() const {
        return _num_info;
    }
    
    const vcf::Info * VcfParser::get_info(unsigned int i) const {
        return _info[i];
    }

    unsigned int VcfParser::num_format() const {
        return _num_format;
    }
    
    const vcf::Format * VcfParser::get_format(unsigned int i) const {
        return _format[i];
    }

    unsigned int VcfParser::num_meta() const {
        return _num_meta;
    }
    
    const vcf::Meta * VcfParser::get_meta(unsigned int i) const {
        return _meta[i];
    }

    unsigned int VcfParser::num_alt() const {
        return _num_alt;
    }
    
    const vcf::Alt * VcfParser::get_alt(unsigned int i) const {
        return _alt[i];
    }

    unsigned int VcfParser::num_samples() const {
        return _num_samples;
    }
    
    const char * VcfParser::get_sample(unsigned int i) const {
        return _samples[i];
    }

    const char * VcfParser::chromosome() const {
        return _chrom;
    }

    unsigned long VcfParser::position() const {
        return _position;
    }

    unsigned int VcfParser::num_ID() const {
        return _num_ID;
    }

    const char * VcfParser::ID(unsigned int i) const {
        return _ID[i];
    }

    unsigned int VcfParser::len_reference() const {
        return _len_reference;
    }

    const char * VcfParser::reference() const {
        return _reference;
    }
    
    vcf::Filter * VcfParser::find_filter(const char * id) {
        for (unsigned int i=0; i<_num_filter; i++) {
            if (!strcmp(id, _filter[i]->get_ID())) return _filter[i];
        }
        return NULL;
    }

    vcf::Format * VcfParser::find_format(const char * id) {
        for (unsigned int i=0; i<_num_format; i++) {
            if (!strcmp(id, _format[i]->get_ID())) return _format[i];
        }
        return NULL;
    }

    vcf::Info * VcfParser::find_info(const char * id) {
        for (unsigned int i=0; i<_num_info; i++) {
            if (!strcmp(id, _info[i]->get_ID())) return _info[i];
        }
        return NULL;
    }

    vcf::Meta * VcfParser::find_meta(const char * id) {
        for (unsigned int i=0; i<_num_meta; i++) {
            if (!strcmp(id, _meta[i]->get_key())) return _meta[i];
        }
        return NULL;
    }
        
    vcf::Alt * VcfParser::find_alt(const char * id) {
        for (unsigned int i=0; i<_num_alt; i++) {
            if (!strcmp(id, _alt[i]->get_ID())) return _alt[i];
        }
        return NULL;
    }

    unsigned int VcfParser::num_alternate() const {
        return _num_alternate;
    }
   
    vcf::AltType VcfParser::alternate_type(unsigned int i) const {
        return _type_alternate[i];
    }

    const char * VcfParser::alternate(unsigned int i) const {
        return _alternate[i];
    }
    
    double VcfParser::quality() const {
        return _quality;
    }

    unsigned int VcfParser::num_failed_tests() const {
        return _num_failed_test;
    }

    const char * VcfParser::failed_test(unsigned int i) const {
        return _failed_test[i];
    }

    unsigned int VcfParser::num_FlagInfo() const {
        return _num_FlagInfo;
    }

    const vcf::FlagInfo VcfParser::FlagInfo(unsigned int i) const {
        return *_FlagInfo[i];
    }

    unsigned int VcfParser::num_IntegerInfo() const {
        return _num_IntegerInfo;
    }

    const vcf::TypeInfo<int>& VcfParser::IntegerInfo(unsigned int i) const {
        return *_IntegerInfo[i];
    }

    unsigned int VcfParser::num_FloatInfo() const {
        return _num_FloatInfo;
    }

    const vcf::TypeInfo<double>& VcfParser::FloatInfo(unsigned int i) const {
        return *_FloatInfo[i];
    }
    
    unsigned int VcfParser::num_CharacterInfo() const {
        return _num_CharacterInfo;
    }

    const vcf::TypeInfo<char>& VcfParser::CharacterInfo(unsigned int i) const {
        return *_CharacterInfo[i];
    }

    unsigned int VcfParser::num_StringInfo() const {
        return _num_StringInfo;
    }

    const vcf::StringInfo& VcfParser::StringInfo(unsigned int i) const {
        return *_StringInfo[i];
    }

    unsigned int VcfParser::num_fields() const {
        return _num_formatEntries;
    }
    
    const vcf::Format& VcfParser::field(unsigned int i) const {
        return *_formatEntries[i];
    }
        
    unsigned int VcfParser::field_index(const char * ID) const {
        for (unsigned int i=0; i<_num_formatEntries; i++) {
            if (!strcmp(ID, _formatEntries[i]->get_ID())) return i;
        }
        return UNKNOWN;
    }

    unsigned int VcfParser::field_rank(unsigned int i) const {
        return _formatRank[i];
    }

    unsigned int VcfParser::field_rank(const char * ID) const {
        return _formatRank[field_index(ID)];
    }

    const vcf::SampleInfo& VcfParser::sample_info(unsigned int i) const {
        return *_sampleInfo[i];
    }

    bool VcfParser::has_AN() const {
        return _has_AN;
    }

    unsigned int VcfParser::AN() const {
        return _AN;
    }

    bool VcfParser::has_AA() const {
        return _has_AA;
    }

    const char * VcfParser::AA_string() const {
        return _AA_string;
    }

    unsigned int VcfParser::AA_index() const {
        return _AA_index;
    }

    bool VcfParser::has_AC() const {
        return _has_AC;
    }

    unsigned int VcfParser::num_AC() const {
        return _num_AC;
    }

    unsigned int VcfParser::AC(unsigned int i) const {
        return _AC[i];
    }

    bool VcfParser::has_AF() const {
        return _has_AF;
    }

    unsigned int VcfParser::num_AF() const {
        return _num_AF;
    }

    double VcfParser::AF(unsigned int i) const {
        return _AF[i];
    }

    bool VcfParser::has_GT() const {
        return _has_GT;
    }

    unsigned int VcfParser::ploidy() const {
        return _ploidy;
    }

    unsigned int VcfParser::num_genotypes() const {
        return _num_genotypes;
    }

    bool VcfParser::GT_phased(unsigned int i) const {
        return _GT_phased[i];
    }

    bool VcfParser::GT_phased() const {
        return _GT_all_phased;
    }

    unsigned int VcfParser::GT(unsigned int sample, unsigned int allele) const {
        return _GT[sample][allele];
    }

    bool VcfParser::has_PL() const {
        return _has_PL;
    }

    bool VcfParser::has_GL() const {
        return _has_GL;
    }

    unsigned int VcfParser::num_missing_GT(){
        unsigned int result = 0;
            if(!_has_GT) throw EggRuntimeError("The current variant has not genotypes data");
            for (unsigned int i=0; i < _num_samples; ++i){
            for(unsigned int j=0; j < _ploidy; ++j){
                    if(_GT[i][j] == UNKNOWN) ++result;
                }
        }
        return result;
    }

    void VcfParser::set_threshold_PL(unsigned int v) {
        _threshold_PL = v;
    }

    unsigned int VcfParser::get_threshold_PL() const {
        return _threshold_PL;
    }

    void VcfParser::set_threshold_GL(unsigned int v) {
        _threshold_GL = v;
    }

    unsigned int VcfParser::get_threshold_GL() const {
        return _threshold_GL;
    }

    void VcfParser::_find_alleles(unsigned int na, unsigned int * dest, unsigned int idx) {
        if (na >= _genotype_idx_helper_size) {
            _genotype_idx_helper = (unsigned int **) realloc(_genotype_idx_helper, (na+1) * sizeof(unsigned int *));
            if (!_genotype_idx_helper) throw EGGMEM;
            for (unsigned int i=_genotype_idx_helper_size; i<na+1; i++) {
                unsigned int ng = i*(i+1)/2;
                _genotype_idx_helper[i] = (unsigned int *) malloc(2 * ng * sizeof(unsigned int));
                if (!_genotype_idx_helper[i]) throw EGGMEM;
                for (unsigned int a=0; a<i; a++) {
                    for (unsigned int b=a; b<i; b++) {
                        unsigned int idx = b*(b+1)/2 + a;
                        _genotype_idx_helper[i][2*idx] = a;
                        _genotype_idx_helper[i][2*idx+1] = b;
                    }
                }
            }
            _genotype_idx_helper_size = na + 1;
        }

        dest[0] = _genotype_idx_helper[na][2*idx];
        dest[1] = _genotype_idx_helper[na][2*idx+1];
    }

    unsigned int VcfParser::PL(unsigned int sample, unsigned int genotype) const {
        return _PL[sample][genotype];
    }

    double VcfParser::GL(unsigned int sample, unsigned int genotype) const {
        return _GL[sample][genotype];
    }

    bool VcfParser::outgroup_AA() const {
        return _outgroup_AA;
    }
    void VcfParser::set_outgroup_AA() {
        _outgroup_AA = true;
    }

    long long int file_end(const char * fname){
        long long int end_stream;
        egglib::VcfParser VCF;
        VCF.open_file(fname);
        end_stream=VCF.file_end();
        return end_stream;
    }
}
