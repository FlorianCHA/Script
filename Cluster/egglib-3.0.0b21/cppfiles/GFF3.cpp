/*
    Copyright 2012,2014 Stéphane De Mita, Mathieu Siol

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
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <new>
#include "egglib.hpp"
#include "Fasta.hpp"
#include "GFF3.hpp"

namespace egglib {

    Feature::Feature() {
        init();
    }

    //////

    Feature::Feature(const Feature& src) {
        init();
        copy(src);
    }

    //////

    Feature& Feature::operator=(const Feature& src) {
        reset();
        copy(src);
        return *this;
    }

    //////

    Feature::~Feature() {
        free();
    }

    //////

    void Feature::reset() {
        _seqid[0] = '\0';
        _source[0] = '\0';
        _type[0] = '\0';
        _num_frag = 0;
         _strand = no_strand;
        _phase = no_phase;
        _ID[0] = '\0';
        _Name[0] = '\0';
        _num_Alias = 0;
        _num_Parent = 0;
        _Target[0] = '\0';
        _Gap[0] = '\0';
        _Derives_from[0] = '\0';
        _num_Note = 0;
        _num_Dbxref = 0;
        _num_Ontology_term = 0;
        _Is_circular = 0;
        _num_attributes = 0;
        _num_parents = 0;
        _num_parts = 0;
    }

    //////

    void Feature::clear() {
        free();
        init();
    }

    //////

    const char * Feature::get_seqid() const {
        return _seqid;
    }

    //////

    void Feature::set_seqid(const char * str) {
        unsigned int n = strlen(str) + 1;
        if (n > _res_seqid) {
            _seqid = (char *) realloc(_seqid, n * sizeof(char));
            if (!_seqid) throw EGGMEM;
            _res_seqid = n;
        }
        strcpy(_seqid, str);
    }

    //////

    const char * Feature::get_source() const {
        return _source;
    }

    //////

    void Feature::set_source(const char * str) {
        unsigned int n = strlen(str) + 1;
        if (n > _res_source) {
            _source = (char *) realloc(_source, n * sizeof(char));
            if (!_source) throw EGGMEM;
            _res_source = n;
        }
        strcpy(_source, str);
    }

    //////

    const char * Feature::get_type() const {
        return _type;
    }

    //////

    void Feature::set_type(const char * str) {
        unsigned int n = strlen(str) + 1;
        if (n > _res_type) {
            _type = (char *) realloc(_type, n * sizeof(char));
            if (!_type) throw EGGMEM;
            _res_type = n;
        }
        strcpy(_type, str);
    }

    //////

    unsigned int Feature::get_num_fragments() const {
        return _num_frag;
    }

    //////
    
    void Feature::set_num_fragments(unsigned int num) {
        if (num > _res_frag) {
            _start = (unsigned int *) realloc(_start, num * sizeof(unsigned int));
            if (!_start) throw EGGMEM;

            _end = (unsigned int *) realloc(_end, num * sizeof(unsigned int));
            if (!_end) throw EGGMEM;
            
            _res_frag = num;
        }
        
        for (unsigned int i=_num_frag; i<num; i++) {
            _start[i] = 0;
            _end[i] = 0;
        }
        
        _num_frag = num;
    }

    //////

    unsigned int Feature::get_start(unsigned int i) const {
        return _start[i];
    }

    //////

    void Feature::set_start(unsigned int i, unsigned int val) {
        _start[i] = val;
    }

    //////

    unsigned int Feature::get_end(unsigned int i) const {
        return _end[i];
    }

    //////

    void Feature::set_end(unsigned int i, unsigned int val) {
        _end[i] = val;
    }

    //////

    double Feature::get_score() const {
        return _score;
    }

    //////

    void Feature::set_score(double d) {
        _score = d;
    }

    //////

    Feature::STRAND Feature::get_strand() const {
        return _strand;
    }

    //////

    void Feature::set_strand(Feature::STRAND s) {
        _strand = s;
    }

    //////

    Feature::PHASE Feature::get_phase() const {
        return _phase;
    }

    //////

    void Feature::set_phase(Feature::PHASE p) {
        _phase = p;
    }

    //////

    const char * Feature::get_ID() const {
        return _ID;
    }

    //////

    void Feature::set_ID(const char * str) {
        unsigned int n = strlen(str) + 1;
        if (n > _res_ID) {
            _ID = (char *) realloc(_ID, n * sizeof(char));
            if (!_ID) throw EGGMEM;
            _res_ID = n;
        }
        strcpy(_ID, str);
    }
   
    //////

    const char * Feature::get_Name() const {
        return _Name;
    }

    //////

    void Feature::set_Name(const char * str) {
        unsigned int n = strlen(str) + 1;
        if (n > _res_Name) {
            _Name = (char *) realloc(_Name, n * sizeof(char));
            if (!_Name) throw EGGMEM;
            _res_Name = n;
        }
        strcpy(_Name, str);
    }

    //////

    unsigned int Feature::get_num_Alias() const {
        return _num_Alias;
    }

    //////

    void Feature::set_num_Alias(unsigned int num) {

        // realloc if needed
        
        if (num > _res_Alias) {
            _res_len_Alias = (unsigned int *) realloc(_res_len_Alias, num * sizeof(unsigned int));
            if (!_res_len_Alias) throw EGGMEM;
            _Alias = (char **) realloc(_Alias, num * sizeof(char *));
            if (!_Alias) throw EGGMEM;
            
            // initialize new values
            
            for (unsigned int i=_res_Alias; i<num; i++) {
                _res_len_Alias[i] = 1;
                _Alias[i] = (char *) malloc(1 * sizeof(char));
                if (!_Alias[i]) throw EGGMEM;
            }
            _res_Alias = num;
        }

        // (re-)initialize (pseudo-)new values
        
        for (unsigned int i=_num_Alias; i<num; i++) {
            _Alias[i][0] ='\0';
        }
        
        _num_Alias = num;
    }
    
    //////

    const char * Feature::get_Alias(unsigned int i) const {
        return _Alias[i];
    }
    
    //////

    void Feature::set_Alias(unsigned int i, const char * s) {
        unsigned int n = strlen(s) + 1;
        
        if (n > _res_len_Alias[i]) {
            _Alias[i] = (char *) realloc(_Alias[i], n * sizeof(char));
            if (!_Alias[i]) throw EGGMEM;
            _res_len_Alias[i] = n;
        }
        
        strcpy(_Alias[i], s);
    }

    //////

    unsigned int Feature::get_num_Parent() const {
        return _num_Parent;
    }
    
    //////

    void Feature::set_num_Parent(unsigned int num) {

        // realloc if needed
        
        if (num > _res_Parent) {
            _res_len_Parent = (unsigned int *) realloc(_res_len_Parent, num * sizeof(unsigned int));
            if (!_res_len_Parent) throw EGGMEM;
            _Parent = (char **) realloc(_Parent, num * sizeof(char *));
            if (!_Parent) throw EGGMEM;
            
            // initialize new values
            
            for (unsigned int i=_res_Parent; i<num; i++) {
                _res_len_Parent[i] = 1;
                _Parent[i] = (char *) malloc(1 * sizeof(char));
                if (!_Parent[i]) throw EGGMEM;
            }
            _res_Parent = num;
        }
        
        // (re-)initialize (pseudo-)new values
        
        for (unsigned int i=_num_Parent; i<num; i++) {
            _Parent[i][0] ='\0';
        }
        
        _num_Parent = num;
    }
    
    //////

    const char * Feature::get_Parent(unsigned int i) const {
        return _Parent[i];
    }
    
    //////

    void Feature::set_Parent(unsigned int i, const char * s) {
        unsigned int n = strlen(s) + 1;
        
        if (n > _res_len_Parent[i]) {
            _Parent[i] = (char *) realloc(_Parent[i], n * sizeof(char));
            if (!_Parent[i]) throw EGGMEM;
            _res_len_Parent[i] = n;
        }
        
        strcpy(_Parent[i], s);
    }

    //////

    const char * Feature::get_Target() const {
        return _Target;
    }

    //////

    void Feature::set_Target(const char * str) {
        unsigned int n = strlen(str) + 1;
        if (n > _res_Target) {
            _Target = (char *) realloc(_Target, n * sizeof(char));
            if (!_Target) throw EGGMEM;
            _res_Target = n;
        }
        strcpy(_Target, str);
    }

    //////

    const char * Feature::get_Gap() const {
        return _Gap;
    }

    //////

    void Feature::set_Gap(const char * str) {
        unsigned int n = strlen(str) + 1;
        if (n > _res_Gap) {
            _Gap = (char *) realloc(_Gap, n * sizeof(char));
            if (!_Gap) throw EGGMEM;
            _res_Gap = n;
        }
        strcpy(_Gap, str);
    }

    //////

    const char * Feature::get_Derives_from() const {
        return _Derives_from;
    }

    //////

    void Feature::set_Derives_from(const char * str) {
        unsigned int n = strlen(str) + 1;
        if (n > _res_Derives_from) {
            _Derives_from = (char *) realloc(_Derives_from, n * sizeof(char));
            if (!_Derives_from) throw EGGMEM;
            _res_Derives_from = n;
        }
        strcpy(_Derives_from, str);
    }

    //////

    unsigned int Feature::get_num_Note() const {
        return _num_Note;
    }

    //////

    void Feature::set_num_Note(unsigned int num) {
        
        // realloc if needed
        
        if (num > _res_Note) {
            _res_len_Note = (unsigned int *) realloc(_res_len_Note, num * sizeof(unsigned int));
            if (!_res_len_Note) throw EGGMEM;
            _Note = (char **) realloc(_Note, num * sizeof(char *));
            if (!_Note) throw EGGMEM;
            
            // initialize new values
            
            for (unsigned int i=_res_Note; i<num; i++) {
                _res_len_Note[i] = 1;
                _Note[i] = (char *) malloc(1 * sizeof(char));
                if (!_Note[i]) throw EGGMEM;
            }
            _res_Note = num;
        }
        
        // (re-)initialize (pseudo-)new values
        
        for (unsigned int i=_num_Note; i<num; i++) {
            _Note[i][0] ='\0';
        }
        
        _num_Note = num;
    }

    //////

    const char * Feature::get_Note(unsigned int i) const {
        return _Note[i];
    }

    //////

    void Feature::set_Note(unsigned int i, const char * s) {
        unsigned int n = strlen(s) + 1;
        
        if (n > _res_len_Note[i]) {
            _Note[i] = (char *) realloc(_Note[i], n * sizeof(char));
            if (!_Note[i]) throw EGGMEM;
            _res_len_Note[i] = n;
        }
        
        strcpy(_Note[i], s);
    }

    //////

    unsigned int Feature::get_num_Dbxref() const {
        return _num_Dbxref;
    }

    //////

    void Feature::set_num_Dbxref(unsigned int num) {

        // realloc if needed
        
        if (num > _res_Dbxref) {
            _res_len_Dbxref = (unsigned int *) realloc(_res_len_Dbxref, num * sizeof(unsigned int));
            if (!_res_len_Dbxref) throw EGGMEM;
            _Dbxref = (char **) realloc(_Dbxref, num * sizeof(char *));
            if (!_Dbxref) throw EGGMEM;
            
            // initialize new values
            
            for (unsigned int i=_res_Dbxref; i<num; i++) {
                _res_len_Dbxref[i] = 1;
                _Dbxref[i] = (char *) malloc(1 * sizeof(char));
                if (!_Dbxref[i]) throw EGGMEM;
            }
            _res_Dbxref = num;
        }
        
        // (re-)initialize (pseudo-)new values
        
        for (unsigned int i=_num_Dbxref; i<num; i++) {
            _Dbxref[i][0] ='\0';
        }
        
        _num_Dbxref = num;
    }

    //////

    const char * Feature::get_Dbxref(unsigned int i) const {
        return _Dbxref[i];
    }

    //////

    void Feature::set_Dbxref(unsigned int i, const char * s) {
        unsigned int n = strlen(s) + 1;
        
        if (n > _res_len_Dbxref[i]) {
            _Dbxref[i] = (char *) realloc(_Dbxref[i], n * sizeof(char));
            if (!_Dbxref[i]) throw EGGMEM;
            _res_len_Dbxref[i] = n;
        }
        
        strcpy(_Dbxref[i], s);
    }

    //////

    unsigned int Feature::get_num_Ontology_term() const {
        return _num_Ontology_term;
    }

    //////

    void Feature::set_num_Ontology_term(unsigned int num) {

        // realloc if needed
        
        if (num > _res_Ontology_term) {
            _res_len_Ontology_term = (unsigned int *) realloc(_res_len_Ontology_term, num * sizeof(unsigned int));
            if (!_res_len_Ontology_term) throw EGGMEM;
            _Ontology_term = (char **) realloc(_Ontology_term, num * sizeof(char *));
            if (!_Ontology_term) throw EGGMEM;
            
            // initialize new values
            
            for (unsigned int i=_res_Ontology_term; i<num; i++) {
                _res_len_Ontology_term[i] = 1;
                _Ontology_term[i] = (char *) malloc(1 * sizeof(char));
                if (!_Ontology_term[i]) throw EGGMEM;
            }
            _res_Ontology_term = num;
        }
        
        // (re-)initialize (pseudo-)new values
        
        for (unsigned int i=_num_Ontology_term; i<num; i++) {
            _Ontology_term[i][0] ='\0';
        }
        
        _num_Ontology_term = num;
    }

    //////

    const char * Feature::get_Ontology_term(unsigned int i) const {
        return _Ontology_term[i];
    }

    //////

    void Feature::set_Ontology_term(unsigned int i, const char * s) {
        unsigned int n = strlen(s) + 1;
        
        if (n > _res_len_Ontology_term[i]) {
            _Ontology_term[i] = (char *) realloc(_Ontology_term[i], n * sizeof(char));
            if (!_Ontology_term[i]) throw EGGMEM;
            _res_len_Ontology_term[i] = n;
        }
        
        strcpy(_Ontology_term[i], s);
    }

    //////

    bool Feature::get_Is_circular() const {
        return _Is_circular;
    }

    //////

    void Feature::set_Is_circular(bool b) {
        _Is_circular = b;
    }
    
    //////

    unsigned int Feature::get_num_attributes() const {
        return _num_attributes;
    }

    //////

    void Feature::set_num_attributes(unsigned int num) {
        
        // realloc if needed
        
        if (num > _res_attributes) {
           
            _num_attributes_items = (unsigned int *) realloc(_num_attributes_items, num * sizeof(unsigned int));
            if (!_num_attributes_items) throw EGGMEM;

            _res_attributes_items = (unsigned int *) realloc(_res_attributes_items, num * sizeof(unsigned int));
            if (!_res_attributes_items) throw EGGMEM;

            _res_len_attributes_key = (unsigned int *) realloc(_res_len_attributes_key, num * sizeof(unsigned int));
            if (!_res_len_attributes_key) throw EGGMEM;

            _res_len_attributes_val = (unsigned int **) realloc(_res_len_attributes_val, num * sizeof(unsigned int *));
            if (!_res_len_attributes_val) throw EGGMEM;

            _attributes_key = (char **) realloc(_attributes_key, num * sizeof(char *));
            if (!_attributes_key) throw EGGMEM;

            _attributes_val = (char ***) realloc(_attributes_val, num * sizeof(char **));
            if (!_attributes_val) throw EGGMEM;
            
            // initialize new values
            
            for (unsigned int i=_res_attributes; i<num; i++) {
                _res_len_attributes_key[i] = 1;
                _attributes_key[i] = (char *) malloc(1 * sizeof(char));
                if (!_attributes_key[i]) throw EGGMEM;

                _num_attributes_items[i] = 0;
                _res_attributes_items[i] = 0;
                _res_len_attributes_val[i] = NULL;
                _attributes_val[i] = NULL;
                
            }
            _res_attributes = num;
        }

        // (re-)initialize (pseudo-)new values

        for (unsigned int i=_num_attributes; i<num; i++) {
                _attributes_key[i][0] ='\0';
                _num_attributes_items[i] = 0;
        }
        
        _num_attributes = num;
    }

    //////

    unsigned int Feature::get_num_items_attribute(unsigned int i) const {
        return _num_attributes_items[i];
    }

    //////

    void Feature::set_num_items_attribute(unsigned int i, unsigned int num) {
        
        // realloc if needed

        if (num > _res_attributes_items[i]) {
            _res_len_attributes_val[i] = (unsigned int *) realloc(_res_len_attributes_val[i], num * sizeof(unsigned int));
            if (!_res_len_attributes_val[i]) throw EGGMEM;
            _attributes_val[i] = (char **) realloc(_attributes_val[i], num * sizeof(char *));
            if (!_attributes_val[i]) throw EGGMEM;
            
            // initialize new values

            for (unsigned int j=_res_attributes_items[i]; j<num; j++) {
                _res_len_attributes_val[i][j] = 1;
                _attributes_val[i][j] = (char *) malloc(1 * sizeof(char *));
                if (!_attributes_val[i][j]) throw EGGMEM;
            }
            _res_attributes_items[i] = num;
        }
        
        // (re-)initialize (pseudo-)new values
        
        for (unsigned int j=_num_attributes_items[i]; j<num; j++) {
            _attributes_val[i][j][0] = '\0';

        }
        
        _num_attributes_items[i] = num;
        
    }

    //////

    const char * Feature::get_attribute_key(unsigned int i) const {
        return _attributes_key[i];
    }

    //////

    void Feature::set_attribute_key(unsigned int i, const char * str) {
        unsigned int n = strlen(str) + 1;
        
        if (n > _res_len_attributes_key[i]) {
            _attributes_key[i] = (char *) realloc(_attributes_key[i], n * sizeof(char));
            if (!_attributes_key[i]) throw EGGMEM;
            _res_len_attributes_key[i] = n;
        }
        
        strcpy(_attributes_key[i], str);
    }

    //////

    const char * Feature::get_attribute_value(unsigned int attr, unsigned int item) const {
        return _attributes_val[attr][item];
    }

    //////

    void Feature::set_attribute_value(unsigned int attr, unsigned int item, const char * str) {
        unsigned int n = strlen(str) + 1;
        
        if (n > _res_len_attributes_val[attr][item]) {
            _attributes_val[attr][item] = (char *) realloc(_attributes_val[attr][item], n * sizeof(char));
            if (!_attributes_val[attr][item]) throw EGGMEM;
            _res_len_attributes_val[attr][item] = n;
        }
        
        strcpy(_attributes_val[attr][item], str);
    }

    //////

    unsigned int Feature::get_num_parents() const {
        return _num_parents;
    }

    //////

    void Feature::set_num_parents(unsigned int num) {
        if (num > _res_parents) {
            _parents = (Feature **) realloc(_parents, num * sizeof(Feature *));
            if (!_parents) throw EGGMEM;
            _res_parents = num;
        }
        for (unsigned int i=_num_parents; i<num; i++) _parents[i] = NULL;
        _num_parents = num;
    }

    //////

    Feature * Feature::get_parent(unsigned int i) const {
        return _parents[i];
    }

    //////

    void Feature::set_parent(unsigned int i, Feature * feat) {
        _parents[i] = feat;
    }

    //////

    unsigned int Feature::get_num_parts() const {
        return _num_parts;
    }

    //////

    void Feature::set_num_parts(unsigned int num) {
        if (num > _res_parts) {
            _parts = (Feature **) realloc(_parts, num * sizeof(Feature *));
            if (!_parts) throw EGGMEM;
            _res_parts = num;
        }
        for (unsigned int i=_num_parts; i<num; i++) _parts[i] = NULL;
        _num_parts = num;
    }

    //////

    Feature * Feature::get_part(unsigned int i) const {
        return _parts[i];
    }

    //////

    void Feature::set_part(unsigned int i, Feature * feat) {
        _parts[i] = feat;
    }

    //////

    void Feature::init() {
        _res_seqid = 1;
        _seqid = (char *) malloc(sizeof(char));
        if (!_seqid) throw EGGMEM;
        _seqid[0] = '\0';
        _res_source = 1;
        _source = (char *) malloc(sizeof(char));
        if (!_source) throw EGGMEM;
        _source[0] = '\0';
        _res_type = 1;
        _type = (char *) malloc(sizeof(char));
        if (!_type) throw EGGMEM;
        _type[0] = '\0';
        _num_frag = 0;
        _res_frag = 0;
        _start = NULL;
        _end  = NULL;
        _score = UNDEF;
        _strand = no_strand;
        _phase = no_phase;
        _res_ID = 1;
        _ID = (char *) malloc(sizeof(char));
        if (!_ID) throw EGGMEM;
        _ID[0] = '\0';
        _res_Name = 1;
        _Name = (char *) malloc(sizeof(char));
        if (!_Name) throw EGGMEM;
        _Name[0] = '\0';
        _num_Alias = 0;
        _res_Alias = 0 ;
        _res_len_Alias = NULL;
        _Alias = NULL;
        _num_Parent = 0;
        _res_Parent = 0;
        _res_len_Parent = NULL;
        _Parent = NULL;
        _res_Target = 1;
        _Target = (char *) malloc(sizeof(char));
        if (!_Target) throw EGGMEM;
        _Target[0] = '\0';
        _res_Gap = 1;
        _Gap = (char *) malloc(sizeof(char));
        if (!_Gap) throw EGGMEM;
        _Gap[0] = '\0';
        _res_Derives_from = 1;
        _Derives_from = (char *) malloc(sizeof(char));
        if (!_Derives_from) throw EGGMEM;
        _Derives_from[0] = '\0';
        _num_Note = 0;
        _res_Note = 0;
        _res_len_Note = NULL;
        _Note = NULL;
        _num_Dbxref = 0;
        _res_Dbxref = 0;
        _res_len_Dbxref = NULL;
        _Dbxref = NULL;
        _num_Ontology_term = 0;
        _res_Ontology_term = 0;
        _res_len_Ontology_term = NULL;
        _Ontology_term = NULL;
        _Is_circular = false;
        _num_attributes = 0;
        _res_attributes = 0;
        _num_attributes_items = NULL;
        _res_attributes_items = NULL;
        _res_len_attributes_key = NULL;
        _res_len_attributes_val = NULL;
        _attributes_key = NULL;
        _attributes_val = NULL;
        _num_parents = 0;
        _res_parents = 0;
        _parents = NULL;
        _num_parts = 0;
        _res_parts = 0;
        _parts = NULL;
    }

    //////

    void Feature::free() {
        if (_seqid) ::free(_seqid);
        if (_source) ::free(_source);
        if (_type) ::free(_type);
        if (_start) ::free(_start);
        if (_end) ::free(_end);
        if (_ID) ::free(_ID);
        if (_Name) ::free(_Name);
        if (_res_len_Alias) ::free(_res_len_Alias);
        for (unsigned int i=0; i<_res_Alias; i++) if (_Alias[i]) ::free(_Alias[i]);
        if (_Alias) ::free(_Alias);
        if (_res_len_Parent) ::free(_res_len_Parent);
        for (unsigned int i=0; i<_res_Parent; i++) if (_Parent[i]) ::free(_Parent[i]);
        if (_Parent) ::free(_Parent);
        if (_Target) ::free(_Target);
        if (_Gap) ::free(_Gap);
        if (_Derives_from) ::free(_Derives_from);
        if (_res_len_Note) ::free(_res_len_Note);
        for (unsigned int i=0; i<_res_Note; i++) if (_Note[i]) ::free(_Note[i]);
        if (_Note) ::free(_Note);
        if (_res_len_Dbxref) ::free(_res_len_Dbxref);
        for (unsigned int i=0; i<_res_Dbxref; i++) if (_Dbxref[i]) ::free(_Dbxref[i]);
        if (_Dbxref) ::free(_Dbxref);
        if (_res_len_Ontology_term) ::free(_res_len_Ontology_term);
        for (unsigned int i=0; i<_res_Ontology_term; i++) if (_Ontology_term[i]) ::free(_Ontology_term[i]);
        if (_Ontology_term) ::free(_Ontology_term);
        for (unsigned int i=0; i<_res_attributes; i++) {
            for (unsigned int j=0; j<_res_attributes_items[i]; j++) {
                if (_attributes_val[i][j]) ::free(_attributes_val[i][j]);
            }
            if (_attributes_key[i]) ::free(_attributes_key[i]);
            if (_res_len_attributes_val[i]) ::free(_res_len_attributes_val[i]);
            if (_attributes_val[i]) ::free(_attributes_val[i]);
        }
        if (_res_len_attributes_key) ::free(_res_len_attributes_key);
        if (_res_len_attributes_val) ::free(_res_len_attributes_val);
        if (_attributes_key) ::free(_attributes_key);
        if (_attributes_val) ::free(_attributes_val);
        if (_num_attributes_items) ::free(_num_attributes_items);
        if (_res_attributes_items) ::free(_res_attributes_items);
        if (_parents) ::free(_parents);
        if (_parts) ::free(_parts);
    }

    //////

    void Feature::copy(const Feature& src) {
        
        set_seqid(src._seqid);
        set_source(src._source);
        set_type(src._type);

        set_num_fragments(src._num_frag);
        for (unsigned int i=0; i<_num_frag; i++) {
            _start[i] = src._start[i];
            _end[i] = src._end[i];
        }

        _score = src._score;
        _strand = src._strand;
        _phase = src._phase;
        set_ID(src._ID);
        set_Name(src._Name);

        set_num_Alias(src._num_Alias);
        for (unsigned int i=0; i<_num_Alias; i++) set_Alias(i, src._Alias[i]);

        set_num_Parent(src._num_Parent);
        for (unsigned int i=0; i<_num_Parent; i++) set_Parent(i, src._Parent[i]);
        
        set_Target(src._Target);
        set_Gap(src._Gap);
        set_Derives_from(src._Derives_from);

        set_num_Note(src._num_Note);
        for (unsigned int i=0; i<_num_Note; i++) set_Note(i, src._Note[i]);
        
        set_num_Dbxref(src._num_Dbxref);
        for (unsigned int i=0; i<_num_Dbxref; i++) set_Dbxref(i, src._Dbxref[i]);

        set_num_Ontology_term(src._num_Ontology_term);
        for (unsigned int i=0; i<_num_Ontology_term; i++) set_Ontology_term(i, src._Ontology_term[i]);
        
        _Is_circular = src._Is_circular;
        
        set_num_attributes(src._num_attributes);
        for (unsigned int i=0; i<_num_attributes; i++) {
            set_attribute_key(i, src._attributes_key[i]);
            set_num_items_attribute(i, src._num_attributes_items[i]);
            for (unsigned int j=0; j<_num_attributes_items[i]; j++) {
                set_attribute_value(i, j, src._attributes_val[i][j]);
            }
        }

        set_num_parents(src._num_parents);
        for (unsigned int i=0; i<_num_parents; i++) _parents[i] = src._parents[i];
        
        set_num_parts(src._num_parts);
        for (unsigned int i=0; i<_num_parts; i++) _parts[i] = src._parts[i];
    }

////////////////////////////////////////////////////////////////////////

    GFF3::GFF3() {
        init();
    }

    //////

    GFF3::~GFF3() {
        free();
    }

    //////

    void GFF3::liberal(bool flag) {
        _liberal = flag;
    }

    //////

    void GFF3::parse(const char* fname) {
        unsigned int n = strlen(fname) + 1;
        if (n > _res_fname) {
            _fname = (char *) realloc(_fname, n * sizeof(char));
            if (!_fname) throw EGGMEM;
            _res_fname = n;
        }
        strcpy(_fname, fname);
        std::ifstream fstream(fname);
        if (!fstream.is_open()) throw EggOpenFileError(fname);
        _parse(fstream);
        _fname[0] = '\0';
    }

    //////

    void GFF3::parse_string(std::string& string) {
        unsigned int n = 9;
        if (n > _res_fname) {
            _fname = (char *) realloc(_fname, n * sizeof(char));
            if (!_fname) throw EGGMEM;
            _res_fname = n;
        }
        strcpy(_fname, "<string>");
        std::istringstream sstream(string);
        _parse(sstream);
        _fname[0] = '\0';
    }

    //////

    void GFF3::parse(std::istream& stream) {
        unsigned int n = 9;
        if (n > _res_fname) {
            _fname = (char *) realloc(_fname, n * sizeof(char));
            if (!_fname) throw EGGMEM;
            _res_fname = n;
        }
        strcpy(_fname, "<stream>");
        _parse(stream);
        _fname[0] = '\0';
    }

    //////
    
    void GFF3::_parse(std::istream& stream) {

        reset();
        _stream = &stream;

        // check format line
        
        get_string(buffer, res_buffer, &GFF3::check_string, &GFF3::stop_tabspace);
        if (strcmp(buffer, "##gff-version")) throw EggFormatError(_fname, currline+1, "GFF3", "first directive should be \"gff-version\", but found: ", '\0', buffer);

        get_string(buffer, res_buffer, &GFF3::check_string, &GFF3::stop_lineEOF, true);
        if (strcmp(buffer, "3")) throw EggFormatError(_fname, currline+1, "GFF3", "this GFF version is not supported: ", '\0', buffer);

        currline++;
        
        // read the rest of the file...
        
        while (true) {
            
            if (!_stream->good()) throw EggFormatError(_fname, currline+1, "GFF3", "cannot read file");
            
            switch (_stream->peek()) {

                case '>':
                    get_fasta();
                    break;

                case '#':
                    _stream->get(curr_ch);
                    if (_stream->peek() == '#') get_directive();
                    else {
                        if (!_stream->good()) throw EggFormatError(_fname, currline+1, "GFF3", "truncated data or problem with file");
                        skip_line();
                    }
                    break;

                default:
                    get_annotation();
            }

            currline++;
            _stream->peek();
            if (_stream->eof()) break;
        }
    }

    //////

    void GFF3::get_fasta() {
        FastaParser fp;
        try {
            fp.set_stream(*_stream);
            fp.read_all(false, _sequences);
        }
        catch (EggFormatError e) {
            throw EggFormatError(_fname, e.line()+currline, "GFF3 (FASTA section)", e.m(), e.character(), e.info());
        }
        if (_sequences.get_nsam_i() == 0) throw EggFormatError(_fname, currline, "GFF3 (FASTA section)", "no sequences found");
    }

    //////

    void GFF3::get_directive() {
       
        _num_metadata++;
        
        if (_num_metadata > _res_metadata) {
            
            _res_len_metadata_key = (unsigned int *) realloc(_res_len_metadata_key, _num_metadata * sizeof(unsigned int));
            if (!_res_len_metadata_key) throw EGGMEM;
            _res_len_metadata_key[_num_metadata-1] = 0;

            _res_len_metadata_val = (unsigned int *) realloc(_res_len_metadata_val, _num_metadata * sizeof(unsigned int));
            if (!_res_len_metadata_val) throw EGGMEM;
            _res_len_metadata_val[_num_metadata-1] = 0;
            
            _metadata_key = (char **) realloc(_metadata_key, _num_metadata * sizeof(char *));
            if (!_metadata_key) throw EGGMEM;
            _metadata_key[_num_metadata-1] = NULL;

            _metadata_val = (char **) realloc(_metadata_val, _num_metadata * sizeof(char *));
            if (!_metadata_val) throw EGGMEM;
            _metadata_val[_num_metadata-1] = NULL;

            _res_metadata = _num_metadata;
        }

        _stream->get(curr_ch);
        if (curr_ch != '#') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected error; a \"#\" is expected here", curr_ch);

        get_string(_metadata_key[_num_metadata-1], _res_len_metadata_key[_num_metadata-1], &GFF3::check_string, &GFF3::stop_tabspacelineEOF);
    
        if (!strcmp(_metadata_key[_num_metadata-1], "#")) {
            if (curr_ch != '\n' && _stream->eof() != true) throw EggFormatError(_fname, currline+1, "GFF3", "don't expect data after \"###\" directive on the same line");
            _num_metadata--;
            mark = _num_features > 0 ? _num_features - 1 : 0;
            return;
        }

        if (!strcmp(_metadata_key[_num_metadata-1], "FASTA")) {
            if (curr_ch != '\n') throw EggFormatError(_fname, currline+1, "GFF3", "don't expect data after \"##FASTA\" directive on the same line");
            if (_stream->eof() == true) throw EggFormatError(_fname, currline+1, "GFF3", "no data after \"##FASTA\" directive (expects sequences)");
            currline++;
            get_fasta();
            _num_metadata--;
            return;
        }

        get_string(_metadata_val[_num_metadata-1], _res_len_metadata_val[_num_metadata-1], &GFF3::check_string, &GFF3::stop_lineEOF, true);
    }

    //////

    void GFF3::skip_line() {
        while (true) {
            if (_stream->eof()) break;
            if (!_stream->good()) throw EggFormatError(_fname, currline+1, "GFF3", "cannot read file");
            if (curr_ch == '\r') {
                _stream->get(curr_ch);
                if (curr_ch != '\n') throw EggFormatError(_fname, currline+1, "GFF3", "expect \"\n\" after \"\r\"", curr_ch);
            }
            if (curr_ch == '\n') break;
            _stream->get(curr_ch);
        }
    }

    //////

    void GFF3::get_annotation() {
        
        unsigned int n;

        // add a feature

        _num_features++;
        
        if (_num_features > _res_features) {
            
            _features = (Feature **) realloc(_features, _num_features * sizeof(Feature *));
            if (!_features) throw EGGMEM;
            
            _features[_num_features-1] = new(std::nothrow) Feature;
            if (!_features[_num_features-1]) throw EGGMEM;
            
            _res_features = _num_features;
        }
        else {
            _features[_num_features-1]->reset();            
        }
        
        // set the number of fragments to 1 by default
        
        Feature * f = _features[_num_features-1];  // only for readability
        f->set_num_fragments(1);

        // get seqid

        get_string(f->_seqid,  f->_res_seqid,  &GFF3::check_stringESC, &GFF3::stop_tab);

        // get source

        get_string(f->_source, f->_res_source, &GFF3::check_stringESC, &GFF3::stop_tab);

        // get type

        get_string(f->_type,   f->_res_type,   &GFF3::check_stringESC, &GFF3::stop_tab);

        // get start

        get_string(buffer, res_buffer, &GFF3::check_integer, &GFF3::stop_tab);
        f->set_start(0, atoi(buffer) - 1);

        // get end

        get_string(buffer, res_buffer, &GFF3::check_integer, &GFF3::stop_tab);
        f->set_end(0, atoi(buffer) - 1);
        if (f->get_start(0) > f->get_end(0)) throw EggFormatError(_fname, currline+1, "GFF3", "start must be >= end");

        // get score

        if (_stream->peek() == '.') {  // looks like a missing data
            _stream->get(curr_ch);

            // recover numbers like .001

            char c = _stream->peek();
            if (c >= '0' && c <= '9') {
                _stream->putback(c);
                get_string(buffer, res_buffer, &GFF3::check_float, &GFF3::stop_tab);
                f->set_score(atof(buffer));
            }

            // really a missing data

            else {
                _stream->get(curr_ch);
                if (curr_ch != '\t') throw EggFormatError(_fname, currline+1, "GFF3", "expect a tabulation after missing (\".\") score", curr_ch);
                // leave f->score() to default
            }
        }
        else {
            get_string(buffer, res_buffer, &GFF3::check_float, &GFF3::stop_tab);
            if (buffer[0] != 'N' && buffer[0] != 'n') f->set_score(atof(buffer));
                    // if buffer[0] is N/n, must be NaN, then leave to Gff3Undefined
        }

        // get strand
        
        n = get_string(buffer, res_buffer, &GFF3::check_string, &GFF3::stop_tab);
        if (n != 1) throw EggFormatError(_fname, currline+1, "GFF3", "invalid strand specification: ", '\0', buffer);

        switch (buffer[0]) {
            case '+':  f->set_strand(Feature::plus);     break;
            case '-':  f->set_strand(Feature::minus);    break;
            case '.':  /* default (no_strand) */     break;
            default:   throw EggFormatError(_fname, currline+1, "GFF3", "invalid strand specification: ", buffer[0]);
        }

        // get phase

        n = get_string(buffer, res_buffer, &GFF3::check_string, &GFF3::stop_tab);
        if (n != 1) throw EggFormatError(_fname, currline+1, "GFF3", "invalid pahse specification: ", '\0', buffer);
        switch (buffer[0]) {
            case '0':  f->set_phase(Feature::zero);   break;
            case '1':  f->set_phase(Feature::one);    break;
            case '2':  f->set_phase(Feature::two);    break;
            case '.':  /* default (no_phase) */   break;
            default:   throw EggFormatError(_fname, currline+1, "GFF3", "invalid phase specification: ", buffer[0]);
        }

        if (!strcmp(f->get_type(), "CDS")) {
            if (_liberal == false && f->get_phase() == Feature::no_phase) throw EggFormatError(_fname, currline+1, "GFF3", "CDS feature must have a phase");
        }
        else {
            if (f->get_phase() != Feature::no_phase) throw EggFormatError(_fname, currline+1, "GFF3", "this type of feature cannot have a phase: ", '\0', f->get_type());
        }

        // get attributes

        int flag = 0;
            // flag to check that attributes are not checked more than once
            // to set bit i: flag |= 1, 2, ..., 1024;
            // to get bit i: flag & 1, 2, ..., 1024 ?;

        while (true) {

            // support trailing ; character (for GFF3 generated by lazy programs)

            if (curr_ch==';') {
                if (_stream->peek() == '\n') {
                    _stream->get(curr_ch);
                    break;
                }
                if (_stream->peek() == '\r') {
                    _stream->get(curr_ch);
                    _stream->get(curr_ch);
                    if (curr_ch != '\n') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected carriage return (not followed by a new line)");
                    break;
                }
                if (_stream->eof()) break;
            }

            // support empty attribute (directly; ) -- also for lazy programs

            if (_stream->peek() == ';') {
                _stream->get(curr_ch);
                continue;
            }

            // get attribute key

            get_string(buffer, res_buffer, &GFF3::check_stringESC, &GFF3::stop_equalsemicolonlineEOF);

            if (curr_ch!='=') {
                if (strcmp(buffer, "Is_circular")) throw EggFormatError(_fname, currline+1, "GFF3", "expect a \"=\" after ", '\0', buffer);
                if (flag & 1024) throw EggFormatError(_fname, currline+1, "GFF3", "attribute \"Is_circular\" defined more than once");
                flag |= 1024;
                f->set_Is_circular(true);
                if (curr_ch == ';') continue;
                break;
            }
            if (!strcmp(buffer, "Is_circular")) throw EggFormatError(_fname, currline+1, "GFF3", "attribute Is_circular doesn't expect arguments");

            //  read value

            if (buffer[0] >= 'A' && buffer[0] <= 'Z') predefined_attribute(f, flag);
            else custom_attribute(f);
    
            if (curr_ch == '\n' || _stream->eof()) break;
        }

        // if ID is defined, check if it is a new fragment of the previous feature
        
        if (flag & 1 && _num_features > 1 && !strcmp(f->_ID, _features[_num_features-2]->_ID)) {

            // check that type is matching (other fields/attributes are no checked)

            if (strcmp(f->get_type(), _features[_num_features-2]->get_type()))   throw EggFormatError(_fname, currline+1, "GFF3", "a feature ID is used several times but type values are inconsistent; ID: ", '\0', f->_ID);

            // check positions
                    
            if (f->get_start(0) <= _features[_num_features-2]->get_end(_features[_num_features-2]->get_num_fragments()-1)) throw EggFormatError(_fname, currline+1, "GFF3", "invalid start position with respect to the end position of the previous fragment; ID: ", '\0', f->_ID);
                    
            // add the fragment to the base feature
                
            _features[_num_features-2]->set_num_fragments(_features[_num_features-2]->get_num_fragments()+1);
            _features[_num_features-2]->set_start(_features[_num_features-2]->get_num_fragments()-1, f->get_start(0));
            _features[_num_features-2]->set_end(_features[_num_features-2]->get_num_fragments()-1, f->get_end(0));

            // forget the new feature and return

            _num_features--;
            return;
        }

        // store all genes in a special array

        if (!strcmp(f->_type, "gene")) {
            _num_genes++;
            if (_num_genes > _res_genes) {
                _genes = (Feature**) realloc(_genes, _num_genes * sizeof(Feature *));
                if (!_genes) throw EGGMEM;
                _res_genes = _num_genes;
            }
            _genes[_num_genes-1] = f;
        }
    }

    //////

    void GFF3::predefined_attribute(Feature * f, int& flag) {
        
        if (!strcmp(buffer, "ID")) {
            if (flag & 1) throw EggFormatError(_fname, currline+1, "GFF3", "attribute \"ID\" defined more than once");
            flag |= 1;
            get_string(f->_ID, f->_res_ID, &GFF3::check_stringESC, &GFF3::stop_semicolonlineEOF);
            return;
        }

        if (!strcmp(buffer, "Name")) {
            if (flag & 2) throw EggFormatError(_fname, currline+1, "GFF3", "attribute \"Name\" defined more than once");
            flag |= 2;
            get_string(f->_Name, f->_res_Name, &GFF3::check_stringESC, &GFF3::stop_semicolonlineEOF);
            return;
        }

        if (!strcmp(buffer, "Target")) {
            if (flag & 4) throw EggFormatError(_fname, currline+1, "GFF3", "attribute \"Target\" defined more than once");
            flag |= 4;
            get_string(f->_Target, f->_res_Target, &GFF3::check_stringESC, &GFF3::stop_semicolonlineEOF);
            return;
        }
        
        if (!strcmp(buffer, "Gap")) {
            if (flag & 8) throw EggFormatError(_fname, currline+1, "GFF3", "attribute \"Gap\" defined more than once");
            flag |= 8;
            get_string(f->_Gap, f->_res_Gap, &GFF3::check_stringESC, &GFF3::stop_semicolonlineEOF);
            return;
        }

        if (!strcmp(buffer, "Derives_from")) {
            if (flag & 16) throw EggFormatError(_fname, currline+1, "GFF3", "attribute \"Derives_from\" defined more than once");
            flag |= 16;
            get_string(f->_Derives_from, f->_res_Derives_from, &GFF3::check_stringESC, &GFF3::stop_semicolonlineEOF);
            return;
        }

        if (!strcmp(buffer, "Alias")) {
            if (flag & 32) throw EggFormatError(_fname, currline+1, "GFF3", "attribute \"Alias\" defined more than once");
            flag |= 32;
            get_items(f, &Feature::set_num_Alias, f->_res_len_Alias, f->_Alias);
            return;
        }

        if (!strcmp(buffer, "Parent")) {
            if (flag & 64) throw EggFormatError(_fname, currline+1, "GFF3", "attribute \"Parent\" defined more than once");
            flag |= 64;
            get_items(f, &Feature::set_num_Parent, f->_res_len_Parent, f->_Parent);

            // locates all parents

            if (_num_features < 2) throw EggFormatError(_fname, currline+1, "GFF3", "this feature cannot have parents because it is the first one");

            f->set_num_parents(f->_num_Parent);

            Feature * parent;
            for (unsigned int i=0; i<f->_num_Parent; i++) {
                parent = NULL;
                unsigned int j = _num_features - 2;
                while (true) {
                    if (!strcmp(_features[j]->_ID, f->_Parent[i])) {
                        parent = _features[j];
                        break;
                    }
                    if (j == mark) break;
                    j--;
                }

                if (parent == NULL) throw EggFormatError(_fname, currline+1, "GFF3", "this parent cannot be identified: ", '\0', f->_Parent[i]);

                f->set_parent(i, parent);
                unsigned int nparts = parent->get_num_parts();
                parent->set_num_parts(nparts + 1);
                parent->set_part(nparts, f);
            }
            return;
        }

        if (!strcmp(buffer, "Note")) {
            if (flag & 128) throw EggFormatError(_fname, currline+1, "GFF3", "attribute \"Note\" defined more than once");
            flag |= 128;
            get_items(f, &Feature::set_num_Note, f->_res_len_Note, f->_Note);
            return;
        }

        if (!strcmp(buffer, "Dbxref")) {
            if (flag & 256) throw EggFormatError(_fname, currline+1, "GFF3", "attribute \"Dbxref\" defined more than once");
            flag |= 256;
            get_items(f, &Feature::set_num_Dbxref, f->_res_len_Dbxref, f->_Dbxref);
            return;
        }

        if (!strcmp(buffer, "Ontology_term")) {
            if (flag & 512) throw EggFormatError(_fname, currline+1, "GFF3", "attribute \"Ontology_term\" defined more than once");
            flag |= 512;
            get_items(f, &Feature::set_num_Ontology_term, f->_res_len_Ontology_term, f->_Ontology_term);
            return;
        }
        
        // if we got there, it means the key is not known
        throw EggFormatError(_fname, currline+1, "GFF3", "capitalized attribute names are reserved - unknown attribute: ", '\0', buffer);
    }

    //////

    void GFF3::custom_attribute(Feature * f) {

        unsigned int attr = f->get_num_attributes();
        f->set_num_attributes(attr + 1);
        f->set_attribute_key(attr, buffer);

        unsigned int n = 0;
        
        while (true) {
            f->set_num_items_attribute(attr, n + 1);
            get_string(f->_attributes_val[attr][n], f->_res_len_attributes_val[attr][n], &GFF3::check_stringESC, &GFF3::stop_semicoloncommalineEOF);
            n++;
            if (curr_ch == ',') continue;
            break; // if ; \t \n or EOF
        }
    }

    //////

    void GFF3::get_items(Feature * f, void (Feature::* num_items)(unsigned int), unsigned int *& res_len, char **& _items) {

        unsigned int n = 0;
        
        while (true) {
            (f->*(num_items))(n+1);
            get_string(_items[n], res_len[n], &GFF3::check_stringESC, &GFF3::stop_semicoloncommalineEOF);
            n++;
            if (curr_ch == ',') continue;
            break; // if ; \t \n or EOF
        }
    }

    //////

    unsigned int GFF3::get_string(char *& where, unsigned int& _res_, bool (GFF3::* check)(), bool (GFF3::* stop)(), bool skip_initial_spaces) {

        curr_pos = 0;

        while (true) {
            _stream->get(curr_ch);
            if (curr_ch == ' ') {
                if (skip_initial_spaces == true) continue;
            }
            else skip_initial_spaces = false; // skip only leading spaces (not interanl
            if ((this->*stop)()) break;
            if (!(this->*check)()) throw EggFormatError(_fname, currline+1, "GFF3", "invalid character found", curr_ch);

            curr_pos++;
            if ((curr_pos+1) > _res_) {
                where = (char *) realloc(where, (curr_pos+1) * sizeof(char));
                if (!where) throw EGGMEM;
                _res_ = (curr_pos+1);
            }
            where[curr_pos-1] = curr_ch;
        }

        where[curr_pos] = '\0';
        if (curr_pos < 1) throw EggFormatError(_fname, currline+1, "GFF3", "empty field or specification here");
        return curr_pos;
    }

    //////

    bool GFF3::stop_equalsemicolonlineEOF() {
        if (_stream->gcount()==0 && _stream->eof()) return true;
        if (curr_ch == '\t') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected tabulation");
        if (curr_ch == '=') return true;
        if (curr_ch == ';') return true;
        if (curr_ch == '\n') return true;
        if (curr_ch == '\r') {
            _stream->get(curr_ch);
            if (curr_ch != '\n') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected carriage return (not followed by a new line)");
            return true;
        }
        return false;
    }

    //////

    bool GFF3::stop_semicolonlineEOF() {
        if (_stream->gcount()==0 && _stream->eof()) return true;
        if (curr_ch == '\t') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected tabulation");
        if (curr_ch == ';') return true;
        if (curr_ch == '\n') return true;
        if (curr_ch == '\r') {
            _stream->get(curr_ch);
            if (curr_ch != '\n') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected carriage return (not followed by a new line)");
            return true;
        }
        return false;
    }

    //////

    bool GFF3::stop_semicoloncommalineEOF() {
        if (_stream->gcount()==0 && _stream->eof()) return true;
        if (curr_ch == '\t') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected tabulation");
        if (curr_ch == ';') return true;
        if (curr_ch == ',') return true;
        if (curr_ch == '\t') return true;
        if (curr_ch == '\n') return true;
        if (curr_ch == '\r') {
            _stream->get(curr_ch);
            if (curr_ch != '\n') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected carriage return (not followed by a new line)");
            return true;
        }
        return false;
    }

    //////

    bool GFF3::stop_lineEOF() {
        if (_stream->gcount()==0 && _stream->eof()) return true;
        if (curr_ch == '\n') return true;
        if (curr_ch == '\r') {
            _stream->get(curr_ch);
            if (curr_ch != '\n') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected carriage return (not followed by a new line)");
            return true;
        }
        return false;
    }

    //////

    bool GFF3::stop_tabspacelineEOF() {
        if (_stream->gcount()==0 && _stream->eof()) return true;
        if (curr_ch == '\t') return true;
        if (curr_ch == ' ') return true;
        if (curr_ch == '\n') return true;
        if (curr_ch == '\r') {
            _stream->get(curr_ch);
            if (curr_ch != '\n') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected carriage return (not followed by a new line)");
            return true;
        }
        return false;
    }

    //////

    bool GFF3::stop_tabspace() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, currline+1, "GFF3", "file truncated");
        if (curr_ch == '\t') return true;
        if (curr_ch == ' ') return true;
        if (curr_ch == '\n') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected end of line");
        if (curr_ch == '\r') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected carriage return");
        return false;
    }

    //////

    bool GFF3::stop_tab() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, currline+1, "GFF3", "file truncated");
        if (curr_ch == '\t') return true;
        if (curr_ch == '\n') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected end of line");
        if (curr_ch == '\r') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected carriage return");
        return false;
    }

    //////

    bool GFF3::stop_equal() {
        if (_stream->gcount()==0 && _stream->eof()) throw EggFormatError(_fname, currline+1, "GFF3", "file truncated");
        if (curr_ch == '=') return true;
        if (curr_ch == '\n') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected end of line");
        if (curr_ch == '\r') throw EggFormatError(_fname, currline+1, "GFF3", "unexpected carriage return");
        return false;
    }

    //////

    bool GFF3::check_string() {
        if (curr_ch < ' ' || curr_ch > '~') return false;
        return true;
    }

    //////

    bool GFF3::check_stringESC() {
        
        if (curr_ch == '%') {
            _stream->get(curr_ch);
            if (curr_ch < '0' || (curr_ch > '9' && curr_ch < 'A') || curr_ch > 'F') {
                _stream->putback(curr_ch);
            }
            else {
                buffer_ESC[0] = curr_ch;
                _stream->get(curr_ch);
                if (curr_ch < '0' || (curr_ch > '9' && curr_ch < 'A') || curr_ch > 'F') {
                    _stream->putback(curr_ch);
                }
                else {
                    buffer_ESC[1] = curr_ch;
                    curr_ch = (char) strtol(buffer_ESC, NULL, 16);
                    return true;
                }
            }
        }
        
        if (curr_ch < ' ' || curr_ch > '~') return false;
        return true;
    }

    //////

    bool GFF3::check_integer() {
        if (curr_ch >= '0' && curr_ch <= '9') return true;
        return false;
    }

    //////

    bool GFF3::check_float() {
        if (curr_ch >= '0' && curr_ch <= '9') return true;
        if (curr_ch == '-') return true;
        if (curr_ch == '.') return true;
        if (curr_ch == 'E') return true;
        if (curr_ch == 'e') return true;
        if (curr_ch == 'N' || curr_ch == 'n') {
            if (curr_pos == 0) return true;
            if (curr_pos == 2 && (buffer[1] == 'A' || buffer[1] == 'a')) return true;
        }
        if ((curr_ch == 'A' || curr_ch == 'a') && curr_pos == 1 && (buffer[0] == 'N' || buffer[0] == 'n')) return true;
        return false;
    }

    //////

    void GFF3::reset() {
        currline = 0;
        _num_metadata = 0;
        curr_ch = '\0';
        _num_features = 0;
        _num_genes = 0;
        mark = 0;
        _sequences.reset(false);
    }

    //////

    void GFF3::clear() {
        free();
        init();
    }

    //////

    void GFF3::init() {
        _fname = (char *) malloc(sizeof(char));
        if (!_fname) throw EGGMEM;
        _fname[0] = '\0';
        _res_fname = 1;
        _stream = NULL;
        buffer = NULL;
        res_buffer = 0;
        buffer_ESC = (char *) malloc(3 * sizeof(char));
        if (!buffer_ESC) throw EGGMEM;
        buffer_ESC[2] = '\0';
        _num_metadata = 0;
        _res_metadata = 0;
        _res_len_metadata_key = NULL;
        _res_len_metadata_val = NULL;
        _metadata_key = NULL;
        _metadata_val = NULL;
        _num_features = 0;
        _res_features = 0;
        _features = NULL;
        _num_genes = 0;
        _res_genes = 0;
        _genes = NULL;
        mark = 0;
        _liberal = false;
    }

    //////

    void GFF3::free() {
        if (_fname) ::free(_fname);
        if (buffer) ::free(buffer);
        if (buffer_ESC) :: free(buffer_ESC);
        for (unsigned int i=0; i<_res_metadata; i++) {
            if (_metadata_key[i]) ::free(_metadata_key[i]);
            if (_metadata_val[i]) ::free(_metadata_val[i]);
        }
        if (_metadata_key) ::free(_metadata_key);
        if (_metadata_val) ::free(_metadata_val);
        if (_res_len_metadata_key) ::free(_res_len_metadata_key);
        if (_res_len_metadata_val) ::free(_res_len_metadata_val);
        for (unsigned int i=0; i<_res_features; i++) {
            if (_features[i]) delete _features[i];
        }
        if (_features) ::free(_features);
        if (_genes) ::free(_genes);
    }

    //////

    unsigned int GFF3::num_metadata() const {
        return _num_metadata;
    }
        
    //////

    const char * GFF3::metadata_key(unsigned int i) const {
        return _metadata_key[i];
    }

    //////

    const char * GFF3::metadata_value(unsigned int i) const {
        return _metadata_val[i];
    }

    //////

    unsigned int GFF3::num_genes() const {
        return _num_genes;
    }

    //////

    Feature& GFF3::gene(unsigned int i) {
        return *_genes[i];
    }

    //////
        
    unsigned int GFF3::num_features() const {
        return _num_features;
    }
        
    //////

    Feature& GFF3::feature(unsigned int i) {
        return *_features[i];
    }

    //////

    const DataHolder& GFF3::sequences() const {
        return _sequences;
    }
}
