/*
    Copyright 2008-2009,2012-2014 Stéphane De Mita, Mathieu Siol

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
#include "Fasta.hpp"
#include <cstdlib>
#include <sstream>
#include <cstring>

namespace egglib {

    //////
    // constructors helpers and memory management

    FastaParser::FastaParser() {
        init();
    }

    //////

    FastaParser::~FastaParser() {
        close();
        free();
    }

    //////

    void FastaParser::clear() {
        free();
        init();
    }

    //////

    void FastaParser::init() {
        _lname = 1;
        _lname_r = 1;
        _name = (char *) malloc(1 * sizeof(char));
        if (!_name) throw EGGMEM;
        _name[0] = '\0';
        _lseq = 0;
        _lseq_r = 0;
        _seq = NULL;
        _ngrp = 0;
        _ngrp_r = 0;
        _grp = NULL;
        _outgroup = false;
        _grp_o = 0;
        _good = false;
        _stream = NULL;
        _lfname_r = 1;
        _fname = (char *) malloc(1 * sizeof(char));
        if (!_fname) throw EGGMEM;
        _fname[0] = '\0';
        _currline = 0;
    }

    //////

    void FastaParser::free() {
        if (_name) ::free(_name);
        if (_seq) ::free(_seq);
        if (_grp) ::free(_grp);
        if (_fname) ::free(_fname);
    }

    //////

    void FastaParser::reset_sequence() {
        _lname = 1;
        _name[0] = '\0';
        _lseq = 0;
        _ngrp = 0;
        _outgroup = false;
        _grp_o = 0;
    }

    //////

    void FastaParser::reserve(unsigned int ln, unsigned int ls, unsigned int ng, unsigned int lf) {

        if ((ln+1) > _lname_r) {
            _name = (char *) realloc(_name, (ln+1) * sizeof(char));
            if (!_name) throw EGGMEM;
            _lname_r = ln+1;
        }

        if (ls > _lseq_r) {
            _seq = (char *) realloc(_seq, ls * sizeof(char));
            if (!_seq) throw EGGMEM;
            _lseq_r = ls;
        }

        if (ng > _ngrp_r) {
            _grp = (unsigned int *) realloc(_grp, ng * sizeof(unsigned int));
            if (!_grp) throw EGGMEM;
            _ngrp_r = ng;
        }

        if ((lf+1) > _lfname_r) {
            _fname = (char *) realloc(_fname, (lf+1) * sizeof(char));
            if (!_fname) throw EGGMEM;
            _lfname_r = lf+1;
        }
    }

    //////
    // stream management

    void FastaParser::open_file(const char * fname) {

        // reset memory

        close();
        reset_sequence();

        // save file name

        unsigned int lfname = strlen(fname);
        if ((lfname+1) > _lfname_r) {
            _fname = (char *) realloc(_fname, (lfname+1) * sizeof(char));
            if (!_fname) throw EGGMEM;
            _lfname_r = lfname+1;
        }

        strcpy(_fname, fname);

        // open stream

        _stream = &_fstream;
        _fstream.open(fname);
        if (!_fstream.is_open()) {
            throw EggOpenFileError(fname);
        }

        // check that first character is '>'

        check();
    }

    //////

    void FastaParser::set_stream(std::istream& stream) {

        // reset memory

        close();
        reset_sequence();

        // save file name

        if (9 > _lfname_r) {
            _fname = (char *) realloc(_fname, 9 * sizeof(char));
            if (!_fname) throw EGGMEM;
            _lfname_r = 9;
        }

        strcpy(_fname, "<stream>");

        // check stream

        _stream = &stream;
        if (! _stream->good()) {
            throw EggArgumentValueError("FastaParser: invalid stream (not good for reading)");
        }

        // check that first character is '>'

        check();
    }

    //////

    void FastaParser::set_string(const char * str) {

        // reset memory

        close();
        reset_sequence();

        // save file name

        if (9 > _lfname_r) {
            _fname = (char *) realloc(_fname, 9 * sizeof(char));
            if (!_fname) throw EGGMEM;
            _lfname_r = 9;
        }

        strcpy(_fname, "<string>");

        // check stream

        _sstream.clear();
        _sstream.str(str);
        _stream = & _sstream;
        if (! _stream->good()) {
            throw EggArgumentValueError("FastaParser: invalid stream (cannot read string for some reasong)");
        }

        // check that first character is '>'

        check();
    }

    //////

    void FastaParser::check() {

        char c;
        _stream->get(c);

        if (_stream->eof()) {
            _good = false;
            return;
        }

        if (_stream->fail()) {
            throw EggFormatError(_fname, _currline+1, "fasta", "cannot read data from file");
        }

        if (c!='>') {
            throw EggFormatError(_fname, _currline+1, "fasta", "a '>' character is expected here");
        }

        _good = true;
    }

    //////

    void FastaParser::close() {
        reset_sequence();
        _stream = NULL;
        if (_fstream.is_open()) _fstream.close();
        _fstream.clear();
        _good = false;
        _fname[0] = '\0';
        _currline = 0;
    }

    //////

    bool FastaParser::good() const {
        return _good;
    }

    //////
    // reading data

    void FastaParser::read_sequence(bool groups, DataHolder * dest) {

        // add a sequence

        unsigned int index = 0;

        if (dest) {
            index = dest->get_nsam_i();
            dest->set_nsam_i(index + 1);
        }
        else reset_sequence();

        // gets the name until a group mark or the end of the line is met

        char ch;
        bool readinglabel = false;

        while (true) {
            _stream->get(ch);
            if (_stream->eof()) throw EggFormatError(_fname, _currline+1, "fasta", "unexpected end of file - file might be truncated");
            if (_stream->fail()) throw EggFormatError(_fname, _currline+1, "fasta", "invalid header (stream error)");

            if (ch=='\n') {
                _currline++;
                break;
            }
            if (ch=='\r') continue;

            if (ch==' ' && groups && _stream->peek() == '@') {
                continue;
            }

            if (ch=='@' && groups) {
                readinglabel = true;
                break;
            }

            if (ch != '\r') {
                if (dest) dest->name_appendch_i(index, ch);
                else name_append(ch);
            }
        }

        // process group label

        unsigned int nlevels = 0;
        bool outgroup = false;

        if (readinglabel) {

            unsigned int label;

            // read all items

            while (true) {
                ch = _stream->peek();

                // the only non-number character allowed is # (outgroup)

                if (ch < '0' || ch > '9') {
                    if (ch == '#') {
                        outgroup = true;
                        if (nlevels > 0) {
                            throw EggFormatError(_fname, _currline+1, "fasta", "invalid group label (outgroup is incompatible with groups)");
                        }
                        _stream->get(ch);  // the # is read

                        // if there is a number next, get it as label

                        ch = _stream->peek();
                        if (ch >= '0' && ch <= '9') (* _stream) >> _grp_o;
                        else _grp_o = 0;
                    }
                    else {
                        throw EggFormatError(_fname, _currline+1, "fasta", "invalid group label (unexpected character found instead of integer)");
                    }
                }

                // if not outgroup, reads label

                else {

                    if (outgroup) {
                        throw EggFormatError(_fname, _currline+1, "fasta", "invalid group label (outgroup is incompatible with groups)");
                    }

                    (* _stream) >> label;
                    if (_stream->fail()) throw EggFormatError(_fname, _currline+1, "fasta", "invalid group label (stream error)");
                    nlevels++;

                    if (dest) {

                        // if the DataHolder has less levels, set all
                        // previous samples to zero at the nevel

                        if (nlevels > dest->get_ngroups()) {
                            dest->set_ngroups(nlevels);
                            for (unsigned int i=0; i<index; i++) dest->set_group_i(i, nlevels-1, 0);
                        }

                        // set the group label

                        dest->set_group_i(index, nlevels-1, label);
                    }

                    else group_append(label);
                }

                // finally: ignores all \r characters

                ch = '\r';
                while (ch == '\r') _stream->get(ch);

                if (ch == '\n') {
                    _currline++;
                    break;
                }
                if (ch == ',') continue;

                throw EggFormatError(_fname, _currline+1, "fasta", "invalid group label (end of line expected here)");
            }
        }

        // gets the sequence

        unsigned int ls = 0;

        while (true) {

            _stream->get(ch);

            if (_stream->eof()) { _good = false; break; }
            if (_stream->fail()) throw EggFormatError(_fname, _currline+1, "fasta", "sequence (stream error)");
            if (ch == '>') break;

            if (ch != '\n' && ch != '\r') {
                ls++;
                if (dest) {
                    dest->set_nsit_i(index, ls);
                    dest->set_i(index, ls-1, (int)ch);
                }
                else seq_append(ch);
            }
            if (ch == '\n') _currline++;
        }

        // move to outgroup if needed

        if (outgroup) {
            if (dest) {
                dest->to_outgroup(index, _grp_o);
                _grp_o = 0;
            }
            else {
                _outgroup = true;
            }
        }

        // set non-defined (or non-parsed) group levels to 0

        else if (dest && nlevels < dest->get_ngroups()) {
            unsigned int ngrp = dest->get_ngroups();
            for (unsigned int i=nlevels; i<ngrp; i++) {
                dest->set_group_i(index, i, 0);
            }
        }
    }

    //////

    void FastaParser::read_all(bool groups, DataHolder& dest) {
        while (_good) {
            read_sequence(groups, &dest);
        }
    }

    //////

    void FastaParser::name_append(char c) {
        _lname++;
        if (_lname > _lname_r) {
            _name = (char *) realloc(_name, _lname * sizeof(char));
            if (!_name) throw EGGMEM;
            _lname_r = _lname;
        }
        _name[_lname-2] = c;
        _name[_lname-1] = '\0';
    }

    //////

    void FastaParser::seq_append(char c) {
        _lseq++;
        if (_lseq > _lseq_r) {
            _seq = (char *) realloc(_seq, _lseq * sizeof(char));
            if (!_seq) throw EGGMEM;
            _lseq_r = _lseq;
        }
        _seq[_lseq-1] = c;
    }

    //////

    void FastaParser::group_append(unsigned int i) {
        _ngrp++;
        if (_ngrp > _ngrp_r) {
            _grp = (unsigned int *) realloc(_grp, _ngrp * sizeof(unsigned int));
            if (!_grp) throw EGGMEM;
            _ngrp_r = _ngrp;
        }
        _grp[_ngrp-1] = i;
    }

    //////
    // Accessors

    const char * FastaParser::name() const {
        return _name;
    }

    //////

    unsigned int FastaParser::ls() const {
        return _lseq;
    }

    //////

    char FastaParser::ch(unsigned int index) const {
        return _seq[index];
    }

    //////

    unsigned int FastaParser::ngroups() const {
        return _ngrp;
    }

    //////

    unsigned int FastaParser::group(unsigned int index) const {
        return _grp[index];
    }

    //////

    bool FastaParser::outgroup() const {
        return _outgroup;
    }

    //////

    unsigned int FastaParser::group_o() const {
        return _grp_o;
    }

    //////
    // Independent parsing functions

    void read_fasta_file(const char * fname, bool groups, DataHolder& dest) {
        FastaParser fp;
        fp.open_file(fname);
        fp.read_all(groups, dest);
    }

    //////

    void read_fasta_string(const std::string str, bool groups, DataHolder& dest) {
        std::istringstream stream(str);
        FastaParser fp;
        fp.set_stream(stream);
        fp.read_all(groups, dest);
    }

    /******************************************************************/
    // Formatting class

    FastaFormatter::FastaFormatter() : BaseFormatter() {

        // load default config parameters
        defaults();

        _mapping = NULL;
        _res_mapping = 0;
    }

    /**********/

    void FastaFormatter::defaults() {
        set_first();
        set_last();
        set_outgroup();
        set_labels();
        set_linelength();
        set_shift_groups();
        set_mapping();
    }

    /**********/

    FastaFormatter::~FastaFormatter() {
        if (_mapping) free(_mapping);
    }

    /**********/

    void FastaFormatter::set_first(unsigned int first) {
        _first = first;
    }

    /**********/

    void FastaFormatter::set_last(unsigned int last) {
        _last = last;
    }

    /**********/

    void FastaFormatter::set_outgroup(bool outgroup) {
        _outgroup = outgroup;
    }

    /**********/

    void FastaFormatter::set_labels(bool labels) {
        _labels = labels;
    }

    /**********/

    void FastaFormatter::set_linelength(unsigned int linelength) {
        _linelength = linelength;
    }

    /**********/

    void FastaFormatter::set_shift_groups(bool shift_groups) {
        _shift_groups = shift_groups ? 1 : 0;
    }

    /**********/

    void FastaFormatter::set_mapping(const char * mapping) {
        unsigned int len = strlen(mapping);
        if (len > 0) {
            _bool_mapping = true;
            if (len > _res_mapping) {
                _mapping = (char * ) realloc(_mapping, len * sizeof(char));
                if (!_mapping) throw EGGMEM;
                _res_mapping = len;
            }
            _len_mapping = len;
            strncpy(_mapping, mapping, len);
        }
        else {
            _bool_mapping = false;
        }
    }

    /**********/

    std::string FastaFormatter::write_string(const DataHolder& src) {
        _cache_stream = _stream;
        _stream = & _sstream;
        _sstream.str("");
        write(src);
        _stream = _cache_stream;
        return _sstream.str();
    }

    /**********/

    char FastaFormatter::_convert(int allele) {
        if (_bool_mapping) {
            if (allele < 0 || allele >= static_cast<int>(_len_mapping)) throw EggInvalidCharacterError(allele);
            return _mapping[allele];
        }
        else {
            return static_cast<char>(allele);
        }
    }

    /**********/

    void FastaFormatter::write(const DataHolder& src) {

        if (!_stream->good()) {
            throw EggRuntimeError("cannot export data data: invalid stream (not good for writing)");
        }

        // variables

        unsigned int nsam = src.get_nsam_i();
        unsigned int ngroups = src.get_ngroups();
        unsigned int c;
        unsigned int ls = 0;
        bool is_matrix = src.get_is_matrix();
        if (is_matrix) ls = src.get_nsit();

        // export specified sequences of the ingroup

        for (unsigned int i=_first; i<nsam && i<=_last; i++) {

            // export header

            (*_stream) << ">" << src.get_name_i(i);

            if (_labels) {
                if (ngroups > 0) {
                    (*_stream) << " @";
                    (*_stream) << src.get_group_i(i, 0) + _shift_groups;
                    for (unsigned int j=1; j<ngroups; j++) {
                        (*_stream) << "," << src.get_group_i(i, j) + _shift_groups;
                    }
                }
            }

            (*_stream) << std::endl;

            // export sequence

            c = 0;
            if (!is_matrix) ls = src.get_nsit_i(i);

            for (unsigned int j=0; j<ls; j++) {
                (*_stream) << _convert(src.get_i(i, j));
                if (_linelength!=0) {
                    c++;
                    if (c==_linelength) {
                        c = 0;
                        (*_stream) << std::endl;
                    }
                }
            }
            if (c!=0 || _linelength==0) (*_stream) << std::endl;
        }

        // export outgroup if requested

        if (_outgroup) {

            nsam = src.get_nsam_o();

            for (unsigned int i=0; i<nsam; i++) {

                // export header

                (*_stream) << ">" << src.get_name_o(i);
                if (_labels) (*_stream) << " @#" << src.get_group_o(i) + _shift_groups;
                (*_stream) << std::endl;

                // export sequence

                c = 0;
                if (!is_matrix) ls = src.get_nsit_o(i);
                for (unsigned int j=0; j<ls; j++) {
                    (*_stream) << _convert(src.get_o(i, j));
                    if (_linelength!=0) {
                        c++;
                        if (c==_linelength) {
                            c = 0;
                            (*_stream) << std::endl;
                        }
                    }
                }
                if (c!=0 || _linelength==0) (*_stream) << std::endl;
            }
        }

        if (!_stream->good()) {
            throw EggRuntimeError("error while writing fasta data");
        }
    }
}
