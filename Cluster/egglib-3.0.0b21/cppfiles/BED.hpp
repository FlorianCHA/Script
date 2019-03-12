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

#ifndef EGGLIB_BED_HPP
#define EGGLIB_BED_HPP
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <limits>
#include "egglib.hpp"
using namespace std;



namespace egglib {

    typedef struct _bed_struct{
        char  * chromosome;
        unsigned long start_pos;
        unsigned long end_pos;
        char * name;
        char * score;
        char strand;
    } _bed_t;


    class BedParser {
        private:
            ifstream _stream;
            char curr_ch;
            char * _fname;
            unsigned int _currline;

            char * _buffer;
            unsigned int _res_chrom;
            unsigned int _res_name;
            unsigned int _res_score;
            unsigned int _res_strand;
            unsigned int _res_buffer;
        
            /* \brief This method allows to check if a line in a bed file belongs to the header
            */
            bool header();

            /* \brief This method allows the parser to pass at another data's line of the ".bed" file.
            *          
            */
            void next();
    
            bool check_sign(unsigned int idx); // all printing characters (`!` to `~`)
            bool stop_spacetabreturn(); // stop at `\t`, '\n', and '\s'
            bool check_integer(unsigned int idx); // `0` to `9`
            bool check_strand(unsigned int idx);
            unsigned int get_string(char *& where, unsigned int& _res_, bool (BedParser::* check)(unsigned int), bool (BedParser::* stop)(), bool catch_missing);
            unsigned int get_char(char & where, unsigned int& _res_, bool (BedParser::* check)(unsigned int), bool (BedParser::* stop)(), bool catch_missing);

            /* \brief this method allows the parser on the ".bed" file to reached the end of the line.
            *  In fact this parser gets only the first three fields (chromosome, startr position and end position),
                        *  so if the data line contains other fields this method allow to skip them, after had got the first
                        *  three field and pass at the next line.       
            */    
            void reach_eol();

        public:
            vector< _bed_t > _Bed_list;
            char * _chrom;
            unsigned long _start_position;
            unsigned long _end_position;
            char * _name;
            char * _score;
            char _strand;
    
            BedParser();
        
            ~BedParser();

            /* \brief this method allows to open a '.bed' file in reading mode  
            *  \param filename path of a '.bed' file   
            */        
            void open_file(const char * filename);
        
            /* \brief this method allows to read data line in a '.bed' file
            *     
            */
            void read();

            /* \brief this method checks if the BedParser reached the end of the '.bed' file
            */
            bool good() const; 
    
            /* \brief
            *    This function allows to read and saves alls chromosomial coordinates of a ".bed" 
            *    file in a "BedParser" vector.
            */
            void get_bed_file(const char * filename);
            
            /* \brief
            *    This function allows the stream on the ".bed" file to pass the header of an bed file 
            *    and reached lines with chromosomal coordinates.
            */    
            void pass_header();
    
            /* \brief
            *     This function allows to get the vector containing alls chromosomal coordinantes got 
            *     from a ".bed" file
            */
            vector< _bed_t > get_bed_list();

            /* \brief this method allows to get a bed_data object containing chromosomial coordinates saved in the   
            *   vector the  current BedParser instance
            */
            _bed_t get_bed_data(int i);

            /* brief this method gets the number of bed_data object saved in the current BedParser instance
            */
            unsigned int n_bed_data();
    };

}

#endif
