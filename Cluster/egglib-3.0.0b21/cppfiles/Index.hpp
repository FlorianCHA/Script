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

#ifndef EGGLIB_INDEX_HPP
#define EGGLIB_INDEX_HPP
using namespace std;
#include <vector>
#include <list>    

namespace egglib {

    typedef struct _index_struct{
        long long int index;
        char chromosome[50];
        unsigned long position;
        unsigned int indice;
    } _index_t;


    class VcfIndex {
        private:
            char * _file_loaded;

            vector< _index_t > _Idx_list;
            
            long long int _index_eof;

            /* \brief gets the vector of index struct loaded in the vcfIndex object
            */
            vector< _index_t > get_index_list();

        public:
            VcfIndex(); ///< Constructor

            ~VcfIndex();

            /* \brief gets the name of the index file loaded in the vcfIndex object
            */
            char * get_file_name();

            /* \brief gets and stores all information of an index file in a list of vector
            */
            void set_index_list(const char * filename);

            /* \brief gets and stores the indexes of all variants of a vcf file in an index file ".vcfi".
            *
            *   This function expects a VCF file in input and the name of the index file desired in output
            *   static void index_file(const char * input, const char * output);
            */

            /* \brief gets the index of the first variant for an desired chromosome 
            *    
            *   This function expects the name of the desired chromosome (const char)
            *   beware: A index file must be loaded in the object "VcfIndex", with the function "set_index_file"
            *          before use the method "get_contigu_index    
            */
            long long int  get_contigu_index(const char * chromosome);

            /* \brief gets the index of a variant for a desired chromosome position
            *  
            *   This function expects a chromosomal position (int)
            *   beware: A index file must be loaded in the object "VcfIndex", with the function "set_index_file"
            *          before use the method "get_position_index
            */
            long long int get_position_index(const char * chromosome, unsigned long position);
            
            /* \brief gets the index of a variant for a desired position variant in the vcf file
            *     (line number without considering the header) 
            *
            *   This function expects a line number  (int).
            *   beware: A index file must be loaded in the object "VcfIndex", with the function "set_index_file"
            *           before use the method "get_indice_index    
            */
            long long int get_indice_index(unsigned int indice);

            unsigned int get_index_indice(long long int index);

            /* brief gets the index of the EOF of the Vcf file, whriting in the Index file. 
            */
            long long int get_index_eof();
        
            /* brief gets the number of index load in the VcfIndex 
            */
            unsigned int n_index();
        
            /* brief gets the index of the last variant of the crhomomome given
            */
            long long int get_last_position_contig_index(const char * chromosome);
    };


    /* \brief gets the EOF of an VCF file.
    */
    long long int get_index_eof(const char * fname);

    /* \brief creates an index file from an VcfParser
    */
    void index_file_from_vcf(egglib::VcfParser& VCF, const char * output);

}

#endif
