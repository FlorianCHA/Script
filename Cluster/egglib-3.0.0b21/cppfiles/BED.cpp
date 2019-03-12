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

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <sstream>
#include <string.h>
#include <limits>
#include "egglib.hpp"
#include "BED.hpp"
using namespace std;

//const char MAXCHAR = std::numeric_limits<char>::max();

namespace egglib {

        BedParser::BedParser(){

            //prev_ch=0;
            _chrom = (char *) malloc(1 * sizeof(char));
            if (!_chrom) throw EGGMEM;
            _chrom[0] = '\0';
            _buffer = (char *) malloc(1 * sizeof(char));
            if (!_buffer) throw EGGMEM;
            _buffer[0] = '\0';
            _fname = (char *) malloc(1 * sizeof(char));
            if (!_fname) throw EGGMEM;
            _fname[0] = '\0';

            _name = (char *) malloc(1 * sizeof(char));
            if (!_name) throw EGGMEM;
            _name[0] = MAXCHAR;

            _score = (char *) malloc(1 * sizeof(char));
            if (!_score) throw EGGMEM;
            _score[0] = MAXCHAR;
        
            _strand = MAXCHAR;


            curr_ch=0;
            _start_position=0;
            _end_position=0;
            _res_chrom = 0;
            _res_buffer = 0;  
            _res_name = 0;
            _res_score = 0;
            _res_strand = 0;

            _currline = 0;    
        }
        
        BedParser::~BedParser(){
            curr_ch=0;
            //prev_ch=0;
            _start_position=0;
            _end_position=0;
            _res_chrom = 0;
            _res_buffer = 0;
            _res_name = 0;
            _res_score = 0;
            _res_strand = 0;
            _currline = 0;


            if(_chrom) free(_chrom);
            if(_buffer) free(_buffer);
            if(_name) free (_name);
            if(_score) free (_score);
            if(_fname) free (_fname);
            if(_Bed_list.size()!=0) _Bed_list.clear();
            if(_stream.is_open()) _stream.close();

        }
    
        void BedParser::open_file(const char * filename){
            _stream.open(filename);
            if(!_stream.is_open()) throw EggOpenFileError(filename);


            unsigned int lfname = strlen(filename) + 1;
            _fname = (char *) realloc(_fname, lfname * sizeof(char));
            if (!_fname) throw EGGMEM;
            //_res_fname = lfname;

            strcpy(_fname, filename);
            pass_header();


        }


        void BedParser::next() {
            _stream.get(curr_ch);
        }

        bool BedParser::check_sign(unsigned int idx) {
            if (curr_ch < '!' || curr_ch > '~') return false;
            return true;
        }


        bool BedParser::stop_spacetabreturn() {
            if (_stream.gcount()==0 && _stream.eof()) { 
                return true;
            }
            if (curr_ch == '\n') return true;
            if (std::isspace(curr_ch)){
                while(std::isspace(curr_ch)){
                    if(curr_ch == '\n') return true;
                    next();
                }
                _stream.unget();
	        return true;
            }
            if (curr_ch == '\t') return true;
            if (curr_ch == '\r') {
                next();
                if (curr_ch != '\n') throw  EggFormatError(_fname, _currline+1, "BED", "unexpected carriage return (not followed by a new line)");
                return true;
            }
            return false;
        }
            

        bool BedParser::check_integer(unsigned int idx) {
            if (curr_ch >= '0' && curr_ch <= '9') return true;
            if (idx == 0 && (curr_ch == '-' || curr_ch == '+')) return true;
            return false;
        }

        bool BedParser::check_strand(unsigned int idx) {
            if (curr_ch == '-') return true;
            if (curr_ch == '+') return true;
            if (curr_ch == '.') return true;
            return false;
        }


        unsigned int BedParser::get_string(char *& where, unsigned int& _res_, bool (BedParser::* check)(unsigned int), bool (BedParser::* stop)(), bool catch_missing) {

            unsigned int n = 0;
            bool dot_escaped = false;
            while (true) {
                next();
                if ((this->*stop)()) break;

                if (!(this->*check)(n)) {
                    if (curr_ch == '.' && dot_escaped == false) {
                        dot_escaped = true;
                    }else{
                        throw EggFormatError(_fname, _currline+1, "BED", "invalid character found", curr_ch);
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
            if (n < 1) throw EggFormatError(_fname, _currline+1, "BED", "empty field or specification here");

            if (catch_missing && !strcmp(where, ".") ) where[0] = MAXCHAR; //

            return n;
        }


        unsigned int BedParser::get_char(char & where, unsigned int& _res_, bool (BedParser::* check)(unsigned int), bool (BedParser::* stop)(), bool catch_missing) {

            //where = MAXCHAR;
            unsigned int n = 0;
            bool dot_escaped = false;
            while (true) {
                next();
                if ((this->*stop)()) break;
                if (!(this->*check)(n)) {
                    if (curr_ch == '.' && dot_escaped == false) {
                        dot_escaped = true;
                    }else{
                        throw EggFormatError(_fname, _currline+1, "BED", "invalid character found", curr_ch);
                    }
                }
                n++;
                where = curr_ch;
            }
		
            if (n < 1) throw EggFormatError(_fname, _currline+1, "BED", "empty field or specification here");

            if (catch_missing && where != '.') where = MAXCHAR; //

            return n;
        }



        void BedParser::reach_eol(){
            while(_stream.get(curr_ch)){
                if (curr_ch == '\n') break;
            }
        }


        void BedParser::read() {

            if (!_stream.is_open() || (!_stream.good()))  throw EggFormatError(_fname, _currline+1, "BED", "error in file or invalid stream");

            unsigned int n =0;
            n = get_string(_chrom, _res_chrom, &BedParser::check_sign, &BedParser::stop_spacetabreturn, true);
            if (_chrom[0] == MAXCHAR) throw EggFormatError(_fname, _currline+1, "BED", "chromosome cannot be missing");

            n = get_string(_buffer, _res_buffer, &BedParser::check_integer, &BedParser::stop_spacetabreturn, true);
            if (_buffer[0] == MAXCHAR) throw EggFormatError(_fname, _currline+1, "BED", "start position cannot be missing");
            _start_position = atol(_buffer);

            n = get_string(_buffer, _res_buffer, &BedParser::check_integer, &BedParser::stop_spacetabreturn, true);
            if (_buffer[0] == MAXCHAR) throw EggFormatError(_fname, _currline+1, "BED", "end position cannot be missing");
            _end_position = atol(_buffer);

            if(curr_ch == '\n') {
                _stream.peek();
                _currline++;
                return;
            }

            n = get_string(_name , _res_name, &BedParser::check_sign, &BedParser::stop_spacetabreturn, true);
            if(curr_ch == '\n') {
                _stream.peek();
                _currline++;
                return;
            }

            n = get_string(_score, _res_score, &BedParser::check_sign, &BedParser::stop_spacetabreturn, true);
            if(curr_ch == '\n') {
                _stream.peek();
                _currline++;
                return;
            }

            n = get_char(_strand, _res_strand, &BedParser::check_strand, &BedParser::stop_spacetabreturn, false);
            if(curr_ch == '\n') {
                _stream.peek();
                _currline++;
                return;
            }

            if(curr_ch != '\n') reach_eol();

            _stream.peek();
            _currline++;
        }

        bool BedParser::good() const {
            if (!_stream.is_open()) return false; // _stream cannot be NULL it is not a pointer
            return _stream.good();
        }


        void BedParser::get_bed_file(const char * filename){
            open_file(filename);
            _bed_t * bed_data=NULL;
            while(good()){
                bed_data = (_bed_t*) realloc(bed_data ,(1*sizeof(_bed_t)));
                if (!bed_data) throw EGGMEM;        
                bed_data->chromosome = (char *) malloc((sizeof(_chrom)-1) * sizeof(char));
                if (!bed_data->chromosome) throw EGGMEM;
                bed_data->name = (char *) malloc((sizeof(_name)-1) * sizeof(char));
                if (!bed_data->name) throw EGGMEM;
                bed_data->score = (char *) malloc((sizeof(_score)-1) * sizeof(char));
                if (!bed_data->score) throw EGGMEM;

                read();

                strncpy(bed_data->chromosome, _chrom, (sizeof(_chrom)-1));
                bed_data->chromosome[sizeof(_chrom)-1] = '\0';
                bed_data->start_pos=_start_position;
                bed_data->end_pos=_end_position;

                if(_name[0] != MAXCHAR){
                    strncpy(bed_data->name, _name, (sizeof(_name)-1));
                    bed_data->name[sizeof(_name)-1] = '\0';
                    memset(_name, 0, strlen(_name));
                    _name[0] = MAXCHAR;
		}else{
                    bed_data->name[0] = MAXCHAR;
		}

                if(_score[0] != MAXCHAR){
                    strncpy(bed_data->score, _score, (sizeof(_score)-1));
                    bed_data->score[sizeof(_score)-1] = '\0';
                    memset(_score, 0, strlen(_score));
                    _score[0] = MAXCHAR;
                }else{
                    bed_data->score[0] = MAXCHAR;
		}

                if(_strand != MAXCHAR){
                    bed_data->strand = _strand; 
                    _strand = MAXCHAR;
                }else{
                    bed_data->strand = MAXCHAR; 
		}

                _Bed_list.push_back((*bed_data));
            }


            free(bed_data);

        }

    
        bool BedParser::header(){
            if (_stream.is_open()==false || (!_stream.good()))  throw EggFormatError(_fname, _currline+1, "BED", "error in file or invalid stream");

            int end_h = _stream.tellg();
            get_string(_buffer, _res_buffer, &BedParser::check_sign, &BedParser::stop_spacetabreturn, true);
            if (strcmp(_buffer, "browser") == 0 || strcmp(_buffer, "track") == 0 || _buffer[0] == '#' ) {
                reach_eol();
                _stream.peek();
                end_h = _stream.tellg();
                return true;
            }else{
                _stream.seekg(end_h, _stream.beg);
                return false;
            }    
        }


        void BedParser::pass_header(){
            bool head = true;
            while(true){
                head = header();
                if(!head) break;
                if(!good()) break;
            }
        }


        vector< _bed_t > BedParser::get_bed_list(){
            if(_Bed_list.size() > 0){
                return _Bed_list;
            }else{
                throw EggRuntimeError("you must load an BED file in the object");
            }
        
        }


        _bed_t BedParser::get_bed_data(int i){
            if(_Bed_list.size() == 0) throw EggRuntimeError("you must load an BED file in the object");
            return _Bed_list.at(i);
        }


        unsigned int BedParser::n_bed_data(){
            return _Bed_list.size();
        }


}
