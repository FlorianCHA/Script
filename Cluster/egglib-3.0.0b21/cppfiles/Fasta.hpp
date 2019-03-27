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

#ifndef EGGLIB_FASTA_HPP
#define EGGLIB_FASTA_HPP

#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <limits>
#include "egglib.hpp"
#include "DataHolder.hpp"
#include "Export.hpp"

namespace egglib {

   /** \brief Sequence-by-sequence Fasta parser
    *
    * \ingroup parsers
    *
    * Read fasta-formatted sequence data from a file specified by name
    * or from an open stream. See the description of the format below.
    *
    *    - Each sequence is preceded by a header limited to a single
    *      line and starting by a ">" character.
    *
    *    - The header length is not limited and all characters are
    *      allowed but white spaces and special characters are
    *      discouraged. The header is terminated by a newline character.
    *
    *    - Group labels are specified a special markup system placed at
    *      the end of the header line. The labels are specified by an
    *      at sign ("@" followed by any integer value ("\@0", "\@1",
    *      "\@2" and so on). It is allowed to define several group
    *      labels for any sequence. In that case, integer values must be
    *      enter consecutively after the at sign, separated by commas,
    *      as in "\@1,3,2" for a sequence belonging to groups 1, 3 and 2
    *      in three different grouping levels. Multiple grouping levels
    *      can be used to specify hierarchical structure, but not only
    *      (independent grouping structure can be freely specified).
    *      The markup "\@#" (at sign or hash sign) specifies an outgroup
    *      sequence. The hash sign may be followed by a single integer
    *      to specify a unique group label. Multiple grouping levels
    *      are not allowed for the outgroup. The group labels of the
    *      ingroup and the outgroup are independent, so the same labels
    *      may be used. The at sign can be preceded by a unique space.
    *      In that case, the parser automatically discards one space
    *      before the at sign (both ">name@1" and ">name @1" are read as
    *      "name") but if there are more than one space, additional
    *      spaces are considered to be part of the name. By default, no
    *      grouping structure is assumed and all sequences are assumed
    *      to be part of the ingroup.
    *
    *    - Group indices are ignored unless specifically specified in a
    *      parser's options.
    *
    *    - The sequence itself continues on following lines until the
    *      next ">" character or the end of the file.
    *
    *    - White spaces, tab and carriage returns are allowed at any
    *      position. They are ignored unless for terminating the header
    *      line. There is no limitation in length and different
    *      sequences can have different lengths.
    *
    *    - Characters case is preserved and imported. Note that, when
    *      *groups* is true and that sequences are placed in a
    *      DataHolder instance, their position in the original fasta
    *      file is lost. Exporting to fasta will automatically place
    *      them at the end of the file.
    *
    * Header: <egglib-cpp/Fasta.hpp>
    *
    */
    class FastaParser {

      public:

       /** \brief Constructor
        *
        * The constructor does not generate an object ready for use.
        * Call to open or set methods is needed before starting to
        * parse data.
        *
        */
        FastaParser();

       /** \brief Destructor
        *
        */
        ~FastaParser();

       /** \brief Reserve memory to speed up data loading
        *
        * This method does not change the size of the data set
        * contained in the instance, but reserves memory in order to
        * speed up subsequent loading of data. The passed values are not
        * required to be accurate. In case the instance has allocated
        * more memory than what is requested, nothing is done (this
        * applies to all parameters independently). It is always valid
        * to use 0 for any values (in that case, no memory is pre
        * allocated for the corresponding array, and memory will be
        * allocated when needed). Notethat one character is always
        * pre-allocated for all names.
        *
        * \param ln expected length of name.
        * \param ls expected length of sequence.
        * \param ng expected number of groups.
        * \param lf expected length of file name.
        *
        */
        void reserve(unsigned int ln, unsigned int ls, unsigned int ng, unsigned int lf);

       /** \brief Open a file for reading
        *
        * This method attempts to open the specified file and to read a
        * single character. If the file cannot be open, an
        * EggOpenFileError exception is thrown; if the read character is
        * not '>', an EggFormatError exception is thrown; if the file is
        * empty, no exception is thrown.
        *
        * In case the instance was already processing a stream, it will
        * be dismissed. The stream created by this method will be closed
        * if another stream is created or set (call to open_file() or set_stream()
        * methods), if the close() method is called or upon object
        * destruction.
        *
        * \param fname name of the fasta-formatted file to open.
        *
        */
        void open_file(const char * fname);

       /** \brief Pass an open stream for reading
        *
        * This method sets the passed stream (which is supposed to have
        * been opened for reading) and attempts to read a single
        * character. If the stream is not open or if data cannot be read
        * from it, an EggArgumentValueError (and not EggOpenFileError)
        * exception is thrown; if the read character is not '>', an
        * EggFormatError exception is thrown; if no data is found, no
        * exception is thrown.
        *
        * In case the instance was already processing a stream, it will
        * be dismissed. The stream passed by this method not be closed
        * by the class even when calling close().
        *
        * \param stream open stream to read fasta-formatted sequences from.
        */
        void set_stream(std::istream& stream);

       /** \brief Pass a string for reading
        *
        * This method opens a reading stream initialized on the passed
        * string and attempts to read a single character. If data cannot
        * be read, an EggArgumentValueError (and not EggOpenFileError)
        * exception is thrown; if the read character is not '>', an
        * EggFormatError exception is thrown; if no data is found, no
        * exception is thrown.
        *
        * In case the instance was already processing a stream, it will
        * be dismissed.
        *
        * \param str a string to be read.
        */
        void set_string(const char * str);

       /** \brief Read a single sequence
        *
        * If the argument *dest* is NULL (default):
        *
        * Read a sequence from the stream and load it in the object
        * memory. Read data can be accessed using name(), ch() and
        * group() methods (plus outgroup() and group_o() for an outgroup
        * sequence). Note that memory allocated for storing data is
        * retained after subsequent calls to read() (but not data
        * themselves). This means that subsequent sequences will be read
        * faster. It also means that, after reading a long sequence,
        * memory will be used until destruction of the object or call to
        * clear() method. Note that read data will be lost as soon as
        * the current stream is dismissed (using the close() method), or
        * a new stream is opened or set, or clear() is caller, or read()
        * is called again with a NULL *dest* argument, but not if read()
        * is called with a non-NULL *dest* argument.
        *
        * If the argument *dest* is **not** NULL:
        *
        * Read a sequence from the stream and load it into the passed
        * DataHolder instance. This will result in the addition of one
        * sequence to the DataHolder. If the argument *groups* is true,
        * the destination DataHolder might be modified. New group levels
        * will be added as needed to accomodate the label(s) found in
        * the sequence. In case the destination instance already
        * contained samples, they will be assumed to belong to group 0
        * for all levels where they were not specified. Warning:
        * *dest* must absolutely be a non-matrix.
        *
        * In either case:
        *
        * If no data can be read (no open stream, stream closed or
        * reached end of file), an EggRuntimeError exception will be
        * thrown.
        *
        * \param groups if false, any group labels found in sequence
        *  headers will be ignored.
        *
        * \param dest if not NULL, destination where to place read data
        * (otherwise, data are stored within the current instance).
        *
        * \return The number of read characters in sequence.
        *
        */
        void read_sequence(bool groups, DataHolder * dest = NULL);

       /** \brief Read a multiple sequences into a DataHolder
        *
        * This method calls read() repetitively passing the DataHolder
        * reference which is filled incrementally, until the end of the
        * fasta stream is reached. If the DataHolder instance already
        * contains sequences, new sequences are appended at the end.
        * Warning: *dest* must absolutely be a non-matrix.
        *
        * \param groups if false, any group labels found in sequence
        * headers will be ignored.
        *
        * \param dest destination where to place read data.
        *
        */
        void read_all(bool groups, DataHolder& dest);

       /** \brief Close the opened file
        *
        * This method closes the file that was opened using the
        * open_file() method. If the file was open using the open_file() method of
        * the same instance, it is actually closed. If the file was
        * passed as a stream using set_stream(), it is forgotten but not
        * closed. If no stream is present, this method does nothing.
        *
        */
        void close();

       /** \brief Check if the instance is good for reading
        *
        * Return true if an open stream is available and if the last
        * reading operation (or by default opening) found that the next
        * character is a '>'.
        *
        */
        bool good() const;

       /** \brief Get the last read name
        *
        * Return a c-string containing the name of the last sequence
        * read by the read() method. By default, an empty string is
        * returned.
        *
        */
        const char * name() const;

       /** \brief Get the length of the last read sequence
        *
        * Return the length of the last sequence read by the read()
        * method. By default, the value is 0.
        *
        */
        unsigned int ls() const;

       /** \brief Get a character of the last read sequence
        *
        * Get the value of a specified index of the last sequence read
        * by the read() method. The index must be valid.
        *
        */
        char ch(unsigned int index) const;

       /** \brief Get the number of group labels specified for the last read sequence
        *
        * Return the number group labels specified for the last sequence
        * read by the read() method. By default, the value is 0.
        *
        */
        unsigned int ngroups() const;

       /** \brief Get a group label of the last read sequence
        *
        * Get one of the group label specified for the last sequence
        * read by the read() method. The index must be valid.
        *
        */
        unsigned int group(unsigned int index) const;

       /** \brief Check if the last read sequence is part of the outgroup
        *
        * Return true if the last sequence read by the read() method is
        * labelled as outgroup. If this method returns true, it is
        * necessary that ngroups() returns 0.
        *
        */
        bool outgroup() const;

       /** \brief Get the outgroup's label
        *
        * Undefined if outgroup() returns false. The default value (in
        * case the label was "@#") is 0.
        *
        */
        unsigned int group_o() const;

       /** \brief Actually clears the memory of the instance
        *
        * Actually frees the memory of the instance. This is useful if
        * a large sequence have been read, in order to really free
        * memory.
        *
        */
        void clear();

      private:

       /** \brief Copy constructor is disabled */
        FastaParser(const FastaParser& src) {}

       /** \brief Copy assignment operator is disabled */
        FastaParser& operator=(const FastaParser& src) { return * this; }

        void init();
        void free();
        void reset_sequence();
        void name_append(char c);
        void seq_append(char c);
        void group_append(unsigned int i);
        void check();

        unsigned int _lname;
        unsigned int _lname_r;
        char * _name;
        unsigned int _lseq;
        unsigned int _lseq_r;
        char * _seq;
        unsigned int _ngrp;
        unsigned int _ngrp_r;
        unsigned int * _grp;
        unsigned int _grp_o;
        bool _outgroup;

        bool _good;
        std::ifstream _fstream;
        std::istringstream _sstream;
        std::istream* _stream;
        unsigned int _lfname_r;
        char * _fname;
        unsigned int _currline;
    };

   /** \brief Multi-sequence fasta parser (from file)
    *
    * \ingroup parsers
    *
    * Read fasta-formatted sequence data from a file specified by name.
    * For format specification, see the documentation of the class
    * FastaParser, which is used behind the scenes.
    *
    * Note that, for optimal performance, the read_multi() method of
    * FastaParser requires only one FastaParser instance (the best is
    * to re-use a single DataHolder instance to take advantage of memory
    * caching).
    *
    * \param fname name of the fasta-formatted file to read.
    * \param groups boolean specifying whether group labels should be
    * imported or ignored. If true, group labels are stripped from names
    * and missing labels are replaced by 0.
    * \param dest reference to the instance where to place sequences. If
    * the object already contains sequences, new sequences will be
    * appended to it. In any case, the destination object must always be
    * a non-matrix.
    *
    * Header: <egglib-cpp/Fasta.hpp>
    *
    */
    void read_fasta_file(const char * fname, bool groups, DataHolder& dest);

   /** \brief Multi-sequence fasta parser (from string)
    *
    * \ingroup parsers
    *
    * Read fasta-formatted sequence data from a raw string. For format
    * specification, see the documentation of the class FastaParser,
    * which is used behind the scenes.
    *
    * \param str string containing fasta-formatted sequences.
    * \param groups boolean specifying whether group labels should be
    * imported or ignored. If true, group labels are stripped from names
    * and missing labels are replaced by 0.
    * \param dest reference to the instance where to place sequences. If
    * the object already contains sequences, new sequences will be
    * appended to it. In any case, the destination object must always be
    * a non-matrix.
    *
    * Header: <egglib-cpp/Fasta.hpp>
    *
    */
    void read_fasta_string(const std::string str, bool groups, DataHolder& dest);

   /** \brief Fasta formatter
    *
    * \ingroup parsers
    *
    * Write genetic data to a file, a string or standard output using
    * the fasta format (formally described in FastaParser). It is
    * required that all exported allele values are exportable as
    * characters (in particular, negative values are never allowed). See
    * the methods documentation for more details, in particular
    * set_mapping() to understand how to map alleles to user-specified
    * characters (such as mapping 0, 1, 2, 3 to *A*, *C*, *G*, *T*, for
    * example).
    *
    * Header: <egglib-cpp/Fasta.hpp>
    *
    */
    class FastaFormatter : public BaseFormatter {

        public:

           /** \brief Constructor
            *
            * Parametrization of the instance can be performed using the
            * setter methods (their names start with *set_*). The
            * default values of the argument of these method represent
            * the default value of the options By default, output is
            * sent to the standard output. A (new) output file can be
            * created at any time using open_file().
            *
            * Concerning the arguments *first* and *last*, please note
            * that sequences will be imported if *last* < *first*, or if
            * *first* is larger than the last index. If *last* is larger
            * than the last index, all last sequences are exported and
            * no error is caused. The default values of *first* and
            * *last* ensure that all sequences are exported.
            *
            */
            FastaFormatter();

           /** \brief Destructor
            *
            * Destroys the object. If an output file is currently open,
            * it is closed at this point.
            *
            */
            virtual ~FastaFormatter();

           /** \brief Sets all parameters to defaults
            *
            */
            void defaults();

           /** \brief Sets index of the first sample to export
            *
            */
            void set_first(unsigned int first = 0);

           /** \brief Sets index of the last sample to export
            *
            */
            void set_last(unsigned int last = MAX);

           /** \brief Specifies whether the outgroup should be exported
            *
            */
            void set_outgroup(bool outgroup = true);

           /** \brief Specifies whether the group labels should be exported
            *
            */
            void set_labels(bool labels = true);

           /** \brief Sets the length of sequence line.
            *
            * If zero, the whole sequence is written on a single line.
            *
            */
            void set_linelength(unsigned int linelength = 50);

           /** \brief Specifies whether group labels should be shifted.
            *
            * If `true`, all group labels are shift of one unit
            * (positively) when exporting.
            *
            */
            void set_shift_groups(bool shift_groups = false);

           /** \brief Specifies whether a character mapping should be used.
            *
            * Use the specified list of characters to map integer
            * allelic values. If a non-empty string is provided, the
            * length of the string must be larger than the largest
            * possible allele values. In that case, the allele values
            * will be used as indexes in order to determine which char
            * from this string must be used for outputting. In the case
            * that this method is used with an empty string, the mapping
            * will not be used and the allele values will be casted
            * directly to characters.
            *
            */
            void set_mapping(const char * mapping = "");

           /** \brief Write fasta-formatted data
            *
            * The parameters specified by the last call to config() (or
            * the defaults) apply). If an output file has be opened with
            * open_file(), data are written into this file. Otherwise, data
            * are written to the standard output.
            *
            */
            void write(const DataHolder& src);

           /** \brief Write fasta-formatted data to string
            *
            * As to_stream() but generates and and returns a string. If
            * there is an open file, it is not touched.
            *
            */
            std::string write_string(const DataHolder& src);

        private:

           /** \brief It is not allowed to copy FastaFormatter instances
            *
            */
            FastaFormatter(const FastaFormatter& src) {}

           /** \brief It is not allowed to copy FastaFormatter instances
            *
            */
            FastaFormatter& operator=(const FastaFormatter& src) { return * this;}

            char _convert(int allele);

            bool _labels;
            bool _outgroup;
            unsigned int _first;
            unsigned int _last;
            unsigned int _linelength;
            unsigned int _shift_groups;
            bool _bool_mapping;
            unsigned int _len_mapping;
            unsigned int _res_mapping;
            char * _mapping;
    };
}

#endif
