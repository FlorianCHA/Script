#include <stdlib.h>
#include <utility> 
#include <iostream>
#include <iterator>
#include <string>
#include <stdio.h>
#include <cstring>
#include <vector>
#include "egglib.hpp"
#include "VCF.hpp"
#include "Index.hpp"
#include <list>  
using namespace std;

/*
  long long int n;
  string a = "65862598745";
  const char * a_cstchar = a.c_str ();
  n=strtoll(a_cstchar, NULL, 10);


*/
namespace egglib {

    VcfIndex::VcfIndex(){
        _file_loaded=NULL;
        _index_eof=0;
    }

    VcfIndex::~VcfIndex(){
        if (_file_loaded) free(_file_loaded);
        if(_Idx_list.size()!=0) _Idx_list.clear();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    vector< _index_t > VcfIndex::get_index_list(){ return _Idx_list; }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void VcfIndex::set_index_list(const char * filename){
        if(_Idx_list.size()!=0) _Idx_list.clear();
        fstream f_stream(filename, ios::in | ios::binary);
        if (!f_stream.is_open()) {
            throw EggOpenFileError(filename);
        }

        unsigned int n_element;
        f_stream.read(reinterpret_cast<char*>(&_index_eof),sizeof(long long int));
        f_stream.read(reinterpret_cast<char*>(&n_element),sizeof(unsigned int));
        _Idx_list.resize(n_element);
        f_stream.read(reinterpret_cast<char*>(&_Idx_list[0]), n_element*sizeof(_index_t));
        
        int s_file=strlen(filename)+2;
        _file_loaded=(char*) realloc(_file_loaded, s_file * sizeof(char));
        strcpy(_file_loaded, filename);    

        f_stream.close();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    long long int  VcfIndex::get_contigu_index(const char * chromosome){    
        if(_Idx_list.size()==0) throw EggRuntimeError("you must load an index file in the object \"VcfIndex\" before use the method \"get_contigu_index\"");
        long long int found=UNKNOWN;

        for (vector< _index_t>::iterator it = _Idx_list.begin(); it != _Idx_list.end(); ++it){
            if(strcmp(chromosome, (*it).chromosome)==0){
                found=(*it).index;
                break;
            }
        }
        if(found == UNKNOWN){
            throw EggInvalidChromosomeIdxError(chromosome, _file_loaded);
        }
        return found;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    long long int  VcfIndex::get_position_index(const char * chromosome, unsigned long position){
        if(_Idx_list.size()==0) throw EggRuntimeError("you must load an index file in the object \"VcfIndex\" before use the method \"get_contigu_index\"");
        long long int found=UNKNOWN;
        if(position != 0){
            for (vector< _index_t>::iterator it = _Idx_list.begin(); it != _Idx_list.end(); ++it){
                if((*it).position == position && strcmp(chromosome,  (*it).chromosome)==0){
                    found=(*it).index;
                    break;
                }
            }
        }else{
            found=get_contigu_index(chromosome);        
        }

        if(found == UNKNOWN){
            throw EggInvalidPositionIdxError(chromosome, position, _file_loaded);
        }

        return found;

    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    long long int  VcfIndex::get_indice_index(unsigned int indice){
    
        if(_Idx_list.size()==0) throw EggRuntimeError("you must load an index file in the object \"VcfIndex\" before use the method \"get_contigu_index\"");
        
        long long int found=UNKNOWN;

        for (vector< _index_t>::iterator it = _Idx_list.begin(); it != _Idx_list.end(); ++it){
            if((*it).indice == indice){
                found=(*it).index;
                break;
            }
        }

        if(found == UNKNOWN){
            throw EggInvalidLineIdxError(indice, _file_loaded);
        }
        return found;

    }


    unsigned int VcfIndex::get_index_indice(long long int index){
        if(_Idx_list.size()==0) throw EggRuntimeError("you must load an index file in the object \"VcfIndex\" before use the method \"get_contigu_index\"");
        unsigned int found=UNKNOWN;

        for (vector< _index_t>::iterator it = _Idx_list.begin(); it != _Idx_list.end(); ++it){
            if((*it).index == index){
                found=(*it).indice;
                break;
            }
        }
        if(found == UNKNOWN){
            throw EggInvalidLineIdxError(index, _file_loaded);
        }
        return found;

    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    char * VcfIndex::get_file_name(){
        if(_file_loaded){
            return _file_loaded;
        }else{
            throw EggRuntimeError("you must load an index file in the object \"VcfIndex\" before use the method \"get_file_name\"");
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    long long int VcfIndex::get_index_eof() {
        if(_index_eof){        
            return _index_eof;
        }else{
            throw EggRuntimeError("you must load an index file in the object \"VcfIndex\" before use the method \"get_index_eof\"");
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    long long int get_index_eof(const char * fname){
        long long int idx_eof;
        fstream f_stream(fname, ios::in | ios::binary);
        if (!f_stream.is_open()) {
            throw EggOpenFileError(fname);
        }
        f_stream.read(reinterpret_cast<char*>(&idx_eof),sizeof(long long int));
        f_stream.close();
        return idx_eof;
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
    long long int VcfIndex::get_last_position_contig_index(const char * chromosome){

        if(_Idx_list.size()==0) throw EggRuntimeError("you must load an index file in the object \"VcfIndex\" before use the method \"get_last_position_contig_index\"");

        long long int found = UNKNOWN;
        unsigned long tmp_max = 0;

        for (vector< _index_t>::reverse_iterator it = _Idx_list.rbegin(); it != _Idx_list.rend(); ++it){
            if(strcmp(chromosome,  (*it).chromosome)==0){
                    found=(*it).index;
                    break;
            }
        }

        if(found == UNKNOWN){
            throw EggInvalidChromosomeIdxError(chromosome, _file_loaded);
        }

        return found;

    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    unsigned int VcfIndex::n_index(){
        return _Idx_list.size();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void index_file_from_vcf(egglib::VcfParser& VCF, const char * output){
 
            unsigned int i=1;
            long long int idx_vrt;
            long long int end_file=VCF.file_end();
            long long int curr_stream=VCF.get_index();
            VCF.set_index(VCF.get_idx_frt_sample());
            vector< _index_t > index_list;

            _index_t * index_data=NULL;

            while(1){
                idx_vrt=VCF.get_index();
                //VCF.read();
                VCF.read_chromosome();
                index_data = (_index_t*) realloc(index_data ,(1*sizeof(_index_t)));
                index_data->index=idx_vrt;
                strcpy(index_data->chromosome, VCF.chromosome());
                index_data->position=VCF.position();//+1;  according with the vcf file
                index_data->indice=i;
                index_list.push_back((*index_data));
                ++i;

                if(!VCF.good()) break;            
            }
            unsigned int n_element=index_list.size();

            fstream idx_stream (output, ios::out | ios::binary);
            if (!idx_stream.is_open()) {
                throw EggOpenFileError(output);
            }
            
            idx_stream.write(reinterpret_cast<const char*>(&end_file),sizeof(end_file));
            idx_stream.write(reinterpret_cast<const char*>(&n_element),sizeof(n_element));
            idx_stream.write(reinterpret_cast<const char*>(&index_list[0]), index_list.size() * sizeof(index_list[0]));
            idx_stream.close();
            free(index_data);
            VCF.set_index(curr_stream);
    }




}
