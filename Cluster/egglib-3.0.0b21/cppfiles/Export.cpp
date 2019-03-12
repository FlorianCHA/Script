/*
    Copyright 2014 Stéphane De Mita, Mathieu Siol

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
#include "Export.hpp"
#include <cstdlib>
#include <ios>
#include <iomanip>

namespace egglib {

    BaseFormatter::BaseFormatter() {
        _is_file = false;
        _stream = & std::cout;
    }

    /**********/

    BaseFormatter::~BaseFormatter() {
        if (_fstream.is_open()) _fstream.close();
    }

    /**********/

    void BaseFormatter::close() {
        if (_fstream.is_open()) _fstream.close();
        _stream = & std::cout;
    }

    /**********/

    void BaseFormatter::write(const char * bit, bool eol) {
        *(_stream) << bit;
        if (eol) *(_stream) << std::endl;
    }

    /**********/

    void BaseFormatter::flush() {
        _stream->flush();
    }

    /**********/

    bool BaseFormatter::open_file(const char * fname) {
        if (_fstream.is_open()) _fstream.close();
        _fstream.open(fname);
        if (_fstream.fail()) {
            _fstream.clear();
            return false;
        }
        _is_file = true;
        _stream = & _fstream;
        return true;
    }

    /**********/

    void BaseFormatter::to_str() {
        if (_fstream.is_open()) _fstream.close();
        _sstream.str("");
        _is_file = false;
        _stream = & _sstream;
    }

    /**********/

    std::string BaseFormatter::get_str() {
        return _sstream.str();
    }

    /**********/

    void BaseFormatter::to_cout() {
        _is_file = false;
        _sstream.str("");
        _stream = & std::cout;
    }

    /**********/

    Export::Export() : BaseFormatter() {
        _ms_res_positions = 0;
        _ms_positions = NULL;
    }

    /**********/

    Export::~Export() {
        if (_ms_positions) free(_ms_positions);
    }

    /**********/

    void Export::newick(const Tree& tree, bool blen, bool eol) {
        const Node * node = tree.root();
        if (node->is_terminal()) {
            (*_stream) << node->label();
        }
        else {
            (*_stream) << '(';
            _newick(tree, tree.node(node->son1()), blen);
            (*_stream) << ',';
            _newick(tree, tree.node(node->son2()), blen);
            (*_stream) <<  ')';
        }
        (*_stream) << ';';
        if (eol) (*_stream) << std::endl;
    }

    /**********/

    void Export::_newick(const Tree& tree, const Node * node, bool blen) {
        if (node->is_terminal()) {
            (*_stream) << node->label();
        }
        else {
            (*_stream) << '(';
            _newick(tree, tree.node(node->son1()), blen);
            (*_stream) << ',';
            _newick(tree, tree.node(node->son2()), blen);
            (*_stream) <<  ')';
        }
        if (blen) (*_stream) << ':' << node->get_L();
    }

    /**********/

    void Export::ms_num_positions(unsigned int n) {
        if (n > _ms_res_positions) {
            _ms_positions = (double *) realloc(_ms_positions, n * sizeof(double));
            if (!_ms_positions) throw EGGMEM;
        }
    }

    /**********/

    void Export::ms_position(unsigned int site, double position) {
        _ms_positions[site] = position;
    }

    /**********/

    void Export::ms_auto_positions(unsigned int n) {
        ms_num_positions(n);
        if (n == 1) {
            _ms_positions[0] = 0.5;
        }
        else {
            for (unsigned int i=0; i<n; i++) {
                _ms_positions[i] = static_cast<double>(i) / (n-1);
            }
        }
    }

    /**********/

    void Export::ms(const DataHolder& data, bool spacer, bool outgroup) {
        unsigned int nsi = data.get_nsam_i();
        unsigned int nso = data.get_nsam_o();
        unsigned int ls = data.get_nsit();
        (*_stream) << "//" << std::endl;
        (*_stream) << "segsites: " << ls << std::endl;
        if (ls > 0) {
            (*_stream) << "positions:";
            (*_stream) << std::fixed;
            (*_stream) << std::setprecision(5);
            for (unsigned int i=0; i<ls; i++) (*_stream) << ' ' << std::setprecision(5) << _ms_positions[i];
            (*_stream) << std::endl;
            _stream->unsetf(std::ios_base::floatfield);
            for (unsigned int i=0; i<nsi; i++) {
                (*_stream) << data.get_i(i, 0);
                for (unsigned int j=1; j<ls; j++) {
                    if (spacer) (*_stream) << ' ';
                    (*_stream) << data.get_i(i, j);
                }
                (*_stream) << std::endl;
            }
            if (outgroup) {
                for (unsigned int i=0; i<nso; i++) {
                    (*_stream) << data.get_o(i, 0);
                    for (unsigned int j=1; j<ls; j++) {
                        if (spacer) (*_stream) << ' ';
                        (*_stream) << data.get_o(i, j);
                    }
                    (*_stream) << std::endl;
                }
            }
        }
        else (*_stream) << std::endl;
    }
}
