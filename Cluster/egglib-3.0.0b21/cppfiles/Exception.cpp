/*
    Copyright 2009,2011,2012,2016 Stéphane De Mita, Mathieu Siol

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
#include <sstream>
#include <iostream>

namespace egglib {

    EggException::EggException() {
        message = "";
    }

    EggException::EggException(const char * m) {
        message = m;
    }
            
    const char * EggException::what() const  throw() {
        return message.c_str();
    }

    EggMemoryError::EggMemoryError(unsigned int line, const char * file) {
        std::ostringstream stream;
        stream << "error while allocating memory";
        stream << " [for debugging: line " << line << " of file " << file << "]";
        message = stream.str();
    }

    EggArgumentValueError::EggArgumentValueError(const char * m) {
        message = "invalid value for function argument: ";
        message += m;
    }

    EggRuntimeError::EggRuntimeError(const char * m) {
        message = m;
    }

    EggFormatError::EggFormatError(const char * fileName, unsigned int line, const char * expectedFormat, const char * m, char c, const char * paste_end) {
        std::ostringstream stream;
        stream << "cannot parse \"";
        stream << expectedFormat;
        stream << "\" data from \"";
        stream << fileName;
        stream << "\" (line ";
        stream << line;
        stream << "): ";
        stream << m;
        stream << paste_end;
        if (c != 0) {
            stream << " [char: ";
            stream << c;
            stream << "]";
        }
        message = stream.str();
        this->c = c;
        this->paste_end = paste_end;
        _line = line;
        _message = m;
    }

    const char * EggFormatError::m() {
        return _message.c_str();
    }
    
    unsigned int EggFormatError::line() {
        return _line;
    }

    char EggFormatError::character() {
        return c;
    }

    const char * EggFormatError::info() {
        return paste_end.c_str();
    }
    
    EggOpenFileError::EggOpenFileError(const char * fileName) {
        message = "error while opening this file: ";
        message+= fileName;
    }

    EggUnalignedError::EggUnalignedError() {
        message = "sequence doesn't match the alignment length";
    }

    EggInvalidAlleleError::EggInvalidAlleleError(int c, unsigned int seqIndex, unsigned int posIndex) {
        std::ostringstream stream;
        stream << "invalid allele value found: " << c;
        if (c>=32 && c<=126) stream << " (character: " << static_cast<char>(c) << ")";
        stream << " at position " << posIndex+1 << " for sample number " << seqIndex+1;
        message = stream.str();
    }

    EggInvalidCharacterError::EggInvalidCharacterError(int value) {
        std::ostringstream stream;
        stream << "cannot export this allelic value: " << value;
        message = stream.str();
    }

    EggPloidyError::EggPloidyError() {
        std::ostringstream stream;
        stream << "invalid ploidy";
        message = stream.str();
    }

    EggNonHierarchicalStructure::EggNonHierarchicalStructure(bool indiv_flag, unsigned int label) {
        std::ostringstream stream;
        stream << "structure is not hierarchical: ";
        if (indiv_flag) stream << "individual";  else stream << "population";
        stream << " " << label << " found in different ";
        if (indiv_flag) stream << "populations";  else stream << "clusters";
        message = stream.str();
    }

    EggNonHierarchicalStructure::EggNonHierarchicalStructure(unsigned int label) {
        std::ostringstream stream;
        stream << "structure is not hierarchical: individual " << label << " found in both ingroup and outgroup";
        message = stream.str();
    }

    EggInvalidChromosomeIdxError::EggInvalidChromosomeIdxError(const char * chromosome, const char * file){
	std::ostringstream stream;
        stream << "The desired chromosome:"<< chromosome <<" doesn't match with chromosomes in the index file:"<< file;
        message = stream.str();
    }

    EggInvalidPositionIdxError::EggInvalidPositionIdxError(const char * chromosome, unsigned int position, const char * file){
	std::ostringstream stream;
        stream << "The desired chromosomal position:"<< position <<" doesn't match with the chromosome:" << chromosome << " in the index file: "<< file;
        message = stream.str();
    }

    EggInvalidLineIdxError::EggInvalidLineIdxError(unsigned int line, const char * file){
	std::ostringstream stream;
        stream << "The desired line:"<< line <<" doesn't match with the lines numbers of the index file:"<< file;
        message = stream.str();
    }

}
