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

#ifndef EGGLIB_VCF_HPP
#define EGGLIB_VCF_HPP
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "egglib.hpp"

namespace egglib {

    class VcfParser;

   /** \brief Classes for parsing VCF files
    *
    * This namespace gathers the helpers of VcfParser (which is itself
    * outside of the namespace to be used conveniently).
    *
    */
    namespace vcf {

       /** \brief Enum for alternate allele specification types
        *
        * The enum is used to specify what the alternate allele value
        * is representing.
        *
        */
        enum AltType {
            Default,    ///< Explicit alternate allele
            X,          ///< The letter "X" instead of a base (requires allow_X())
            Referred,   ///< ID (without angle brackets) of a Alt specification
            Breakend    ///< A breakend specification, provided as is
        };

        ////////////////////////////////////////////////////////////////////

        const char *  _index_from_fname (const char * filename);

        ////////////////////////////////////////////////////////////////////

       /** \brief Large positive value used for a special case
        *
        * This value should means that the variable should be
        * re-evaluated to the number of alternate alleles.
        *
        */
        extern const unsigned int NUM_ALTERNATE;

       /** \brief Large positive value used for a special case
        *
        * This value should means that the variable should be
        * re-evaluated to the number of genotypes.
        *
        */
        extern const unsigned int NUM_GENOTYPES;

       /** \brief Large positive value used for a special case
        *
        * This value should means that the variable should be
        * re-evaluated to the number of alternate alleles + 1.
        *
        */
        extern const unsigned int NUM_POSSIBLE_ALLELES;

        ////////////////////////////////////////////////////////////////////

       /** \brief Class to handle VCF FILTER specifications
        *
        * \ingroup parsers
        *
        * By default, string accessors return null pointers.
        *
        * Header: <egglib-cpp/VCF.hpp>
        *
        */
        class Filter {

          public:

           /** \brief Default constructor */
            Filter();

           /** \brief Initialization constructor */
            Filter(const char * id, const char * descr);

           /** \brief Copy constructor */
            Filter(const Filter& src);

           /** \brief Copy assignment operator */
            Filter& operator=(const Filter& src);

           /** \brief Destructor */
            virtual ~Filter();

           /** \brief Setter (reset extra fields) */
            void update(const char * id, const char * descr);

           /** \brief Get ID string */
            const char * get_ID() const;

           /** \brief Set ID string */
            void set_ID(const char * id);

           /** \brief Get description string */
            const char * get_description() const;

           /** \brief Set description string */
            void set_description(const char * descr);

           /** \brief Actually free memory */
            void clear();

           /** \brief Set extra field */
            void set_extra(const char * key, const char * value);

           /** \brief Get extra field key */
            const char * get_extra_key(unsigned int idx) const;

           /** \brief Get extra field value */
            const char * get_extra_value(unsigned int idx) const;

           /** \brief Get number of extra fields */
            unsigned int get_num_extra() const;

          protected:

            void init();
            void free();
            void copy(const Filter& src);

            unsigned int _ID_r;
            char * _ID;

            unsigned int _descr_r;
            char * _descr;

            unsigned int _extra_n;
            unsigned int _extra_r;
            unsigned int * _extra_key_r;    // size: _extra_r
            unsigned int * _extra_val_r;    // size: _extra_r
            char ** _extra_key;             // size: _extra_r x _extra_key_r[i]
            char ** _extra_val;             // size: _extra_r x _extra_val_r[i]
        };

        ////////////////////////////////////////////////////////////////////

       /** \brief Class to handle VCF INFO specifications
        *
        * \ingroup parsers
        *
        * By default, string accessors return null pointers and other
        * accessors return undefined values. The number of values might
        * be egglib::UNKNOWN (unspecified number of values, represented
        * by the character "." in files) or vcf::NUM_ALTERNATE (match
        * the number of ALT variants, represented by the character "A"
        * in files) or vcf::NUM_GENOTYPES (match the number of possible
        * genotypes, represented by the character "G" in files). Note
        * that UNKNOWN, vcf::NUM_ALTERNATE and vcf::NUM_GENOTYPES are
        * all large positive values. One more special value:
        * vcf::NUM_POSSIBLE_ALLELES (like vcf::NUM_ALTERNATE but
        * including the reference).
        *
        * Header: <egglib-cpp/VCF.hpp>
        *
        */
        class Info : public Filter {

          public:

           /** \brief Meta-information types
            *
            * This enum is used to specify FORMAT and INFO types (Flag is not
            * accepted only for FORMAT).
            *
            */
            enum Type { Integer, Float, Flag, Character, String };

           /** \brief Default constructor */
            Info() : Filter() {}

           /** \brief Initialization constructor */
            Info(const char * id, unsigned int num, Info::Type t, const char * descr);

           /** \brief Copy constructor */
            Info(const Info& src);

           /** \brief Copy assignment operator */
            Info& operator=(const Info& src);

           /** \brief Setter */
            void update(const char * id, unsigned int num, Info::Type t, const char * descr);

           /** \brief Get number of values */
            unsigned int get_number() const;

           /** \brief Set number of values */
            void set_number(unsigned int num);

           /** \brief Get type */
            virtual Info::Type get_type() const;

           /** \brief Set type */
            virtual  void set_type(Info::Type t);

          protected:

            unsigned int _number;
            Info::Type _type;
        };

        ////////////////////////////////////////////////////////////////////

       /** \brief Class to handle VCF FORMAT specifications
        *
        * Format is identical in structure to Info except that the
        * type Flag is not allowed (an EggArgumentValueError is caused when
        * attempting to set Flag as type).
        *
        * \ingroup parsers
        *
        * Header: <egglib-cpp/VCF.hpp>
        *
        */
        class Format : public Info {

          public:

           /** \brief Default constructor */
            Format() : Info() {}

           /** \brief Initialization constructor */
            Format(const char * id, unsigned int num, Info::Type t, const char * descr);

           /** \brief Copy constructor */
            Format(const Format& src);

           /** \brief Copy constructor */
            Format(const Info& src);

           /** \brief Copy assignment operator */
            Format& operator=(const Format& src);

           /** \brief Copy assignment operator */
            Format& operator=(const Info& src);

           /** \brief Setter */
            void update(const char * id, unsigned int num, Info::Type t, const char * descr);

           /** \brief Get type */
            virtual Info::Type get_type() const;

           /** \brief Set type */
            virtual void set_type(Info::Type t);
        };

        ////////////////////////////////////////////////////////////////////

       /** \brief Class to handle VCF meta-information entries
        *
        * \ingroup parsers
        *
        * By default, string accessors return null pointers.
        *
        * Header: <egglib-cpp/VCF.hpp>
        *
        */
        class Meta {

          public:

           /** \brief Default constructor */
            Meta();

           /** \brief Initialization constructor */
            Meta(const char * k, const char * v);

           /** \brief Copy constructor */
            Meta(Meta& src);

           /** \brief Copy assignemnt operator */
            Meta& operator=(Meta& src);

           /** \brief Destructor */
            ~Meta();

           /** \brief Setter */
            void update(const char * k, const char * v);

           /** \brief Set key */
            void set_key(const char * k);

           /** \brief Set value */
            void set_value(const char * v);

           /** \brief Get key */
            const char * get_key() const;

           /** \brief Get value */
            const char * get_value() const;

           /** \brief Actually free memory */
            void clear();

          private:

            void init();
            void free();

            unsigned int _key_r;
            unsigned int _val_r;
            char * _key;
            char * _val;
        };

        ////////////////////////////////////////////////////////////////////

       /** \brief Class to handle VCF alternate allele definitions
        *
        * \ingroup parsers
        *
        * By default, string accessors return null pointers.
        *
        * Header: <egglib-cpp/VCF.hpp>
        *
        */
        class Alt : public Filter {

          public:

           /** \brief Default constructor */
            Alt() : Filter() {}

           /** \brief Initialization constructor */
            Alt(const char * id, const char * descr) : Filter(id, descr) {}

           /** \brief Copy constructor */
            Alt(Alt& src) : Filter(src) {}

           /** \brief Copy constructor */
            Alt(Filter& src) : Filter(src) {}

           /** \brief Copy assignemnt operator */
            Alt& operator=(Alt& src);

           /** \brief Copy assignemnt operator */
            Alt& operator=(Filter& src);

           /** \brief Destructor */
            ~Alt() {}
        };

        ////////////////////////////////////////////////////////////////////

       /** \brief Class representing a Flag-type INFO field
        *
        * \ingroup parsers
        *
        * Header: <egglib-cpp/VCF.hpp>
        *
        */
        class FlagInfo {

          public:

            friend class egglib::VcfParser;

           /** \brief Constructor */
            FlagInfo();

           /** \brief Copy constructor */
            FlagInfo(const FlagInfo& src);

           /** \brief Copy assignment operator */
            FlagInfo& operator=(const FlagInfo& src);

           /** \brief Destructor */
            virtual ~FlagInfo();

           /** \brief Set ID string */
            void set_ID(const char * id);

           /** \brief Get ID string */
            const char * get_ID() const;

          protected:

            virtual void copy(const FlagInfo& src);

            unsigned int _res_ID;
            char * _ID;
        };

        ////////////////////////////////////////////////////////////////////

       /** \brief Template for Character (char), Integer (int) and Float (double)-type INFO fields
        *
        * \ingroup parsers
        *
        * Header: <egglib.cpp/VCF.hpp>
        *
        */
        template <class T> class TypeInfo : public FlagInfo {

          public:

            friend class egglib::VcfParser;

           /** \brief Constructor */
            TypeInfo() {
                _res_ID = 0;
                _ID = NULL;
                _num_items = 0;
                _res_items = 0;
                _items = NULL;
                _expected_number = UNKNOWN;
            }

           /** \brief Copy constructor */
            TypeInfo(const TypeInfo<T>& src) {
                _res_ID = 0;
                _ID = NULL;
                _num_items = 0;
                _res_items = 0;
                _items = NULL;
                copy(src);
            }

           /** \brief Copy assignment operator */
            TypeInfo<T>& operator=(const TypeInfo<T>& src) {
                if (_ID) ::free(_ID);
                if (_items) ::free(_items);
                _res_ID = 0;
                _ID = NULL;
                _num_items = 0;
                _res_items = 0;
                _items = NULL;
                copy(src);
            }

           /** \brief Destructor */
            virtual ~TypeInfo() {
                if (_items) ::free(_items);
            }

           /** \brief Get number of items available in the instance */
            unsigned int num_items() const {
                return _num_items;
            }

           /** \brief Get an item (missing data are encoded by type-specific special values) */
            const T& item(unsigned int i) const {
                return _items[i];
            }

           /** \brief Reset instance */
            void reset() {
                _num_items = 0;
            }

           /** \brief Get expected number of items */
            unsigned int get_expected_number() const {
                return _expected_number;
            }

           /** \brief Set expected number of items */
            void set_expected_number(unsigned int n) {
                _expected_number = n;
            }

          protected:

            virtual void copy(const TypeInfo& src) {
                set_ID(src._ID);
                _num_items = src._num_items;
                if (_num_items > _res_items) {
                    _items = (T *) realloc(_items, _num_items * sizeof(T));
                    if (!_items) throw EGGMEM;
                    _res_items = _num_items;
                }
                for (unsigned int i=0; i<_num_items; i++) {
                    _items[i] = src._items[i];
                }
                _expected_number = src._expected_number;
            }

           /** \brief Add a new item and return its reference */
            void add() {
                _num_items++;
                if (_num_items > _res_items) {
                    _items = (T *) realloc(_items, _num_items * sizeof(T));
                    if (!_items) throw EGGMEM;
                    _res_items = _num_items;
                }
            }

            unsigned int _expected_number;
            unsigned int _num_items;
            unsigned int _res_items;
            T * _items;
        };

        ////////////////////////////////////////////////////////////////////

        /** \brief Class from String-type INFO fields
         *
         * \ingroup parsers
         *
         * Header: <egglib-cpp/VCF.hpp>
         *
         */
        class StringInfo : public TypeInfo<char *> {

          public:

            friend class egglib::VcfParser;

           /** \brief Constructor */
            StringInfo();

           /** \brief Copy constructor */
            StringInfo(const StringInfo& src);

           /** \brief Copy assignment operator */
            StringInfo& operator=(const StringInfo& src);

           /** \brief Destructor */
            virtual ~StringInfo();

           /** \brief Set a character (must fit in length) */
            void change(unsigned int item, unsigned int position, char value);

          protected:

           /** \brief Add a new item and return its reference */
            void add();

            virtual void copy(const StringInfo& src);
            unsigned int * _res_len_items;
        };

        ////////////////////////////////////////////////////////////////////

       /** \brief Class storing information fields for a sample
        *
        * \ingroup parsers
        *
        * Header: <egglib-cpp/VCF.hpp>
        *
        */
        class SampleInfo {

          public:

            friend class egglib::VcfParser;

           /** \brief Constructor */
            SampleInfo();

           /** \brief Copy constructor */
            SampleInfo(const SampleInfo& src);

           /** \brief Copy assignment operator */
            SampleInfo& operator=(const SampleInfo& src);

           /** \brief Destructor */
            ~SampleInfo();

           /** \brief Restore the object to the initial state
            *
            * This method does not free allocated memory, which is reserved
            * for latter use.
            *
            */
            void reset();

           /** \brief Actually free all memory allocated by this instance */
            void clear();

           /** \brief Number of Integer-type entries */
            unsigned int num_IntegerEntries() const;

           /** \brief Number of items for an Integer-type entry */
            unsigned int num_IntegerItems(unsigned int i) const;

           /** \brief Get an item for an Integer-type entry */
            int IntegerItem(unsigned int i, unsigned int j) const;

           /** \brief Number of Float-type entries */
            unsigned int num_FloatEntries() const;

           /** \brief Number of items for an Float-type entry */
            unsigned int num_FloatItems(unsigned int i) const;

           /** \brief Get an item for an Float-type entry */
            double FloatItem(unsigned int i, unsigned int j) const;

           /** \brief Number of Character-type entries */
            unsigned int num_CharacterEntries() const;

           /** \brief Number of items for an Character-type entry */
            unsigned int num_CharacterItems(unsigned int i) const;

           /** \brief Get an item for an Character-type entry */
            char CharacterItem(unsigned int i, unsigned int j) const;

           /** \brief Number of String-type entries */
            unsigned int num_StringEntries() const;

           /** \brief Number of items for an String-type entry */
            unsigned int num_StringItems(unsigned int i) const;

           /** \brief Get an item for an String-type entry */
            const char * StringItem(unsigned int i, unsigned int j) const;

          private:

            void init();
            void copy(const SampleInfo& src);
            void free();

            void addIntegerEntry();
            void addIntegerItem();
            void addFloatEntry();
            void addFloatItem();
            void addCharacterEntry();
            void addCharacterItem();
            void addStringEntry();
            void addStringItem();

            unsigned int _num_IntegerEntries;
            unsigned int _res_IntegerEntries;
            unsigned int * _num_IntegerItems;
            unsigned int * _res_IntegerItems;
            int ** _IntegerItems;

            unsigned int _num_FloatEntries;
            unsigned int _res_FloatEntries;
            unsigned int * _num_FloatItems;
            unsigned int * _res_FloatItems;
            double ** _FloatItems;

            unsigned int _num_CharacterEntries;
            unsigned int _res_CharacterEntries;
            unsigned int * _num_CharacterItems;
            unsigned int * _res_CharacterItems;
            char ** _CharacterItems;

            unsigned int _num_StringEntries;
            unsigned int _res_StringEntries;
            unsigned int * _num_StringItems;
            unsigned int * _res_StringItems;
            unsigned int ** _res_len_StringItems;
            char *** _StringItems;
        };
    }
    ////////////////////////////////////////////////////////////////////

   /** \brief Line-by-line variant call format (VCF) parser
    *
    * \ingroup parsers
    *
    * Read VCF-formatted variation data from a file specified by name
    * or from an open stream.
    *
    * This parser supports the 4.1 specification of the variant call
    * format as described at this address:
    *
    * http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
    *
    * Upon opening a file (or passing an open stream), VcfParser objects
    * automatically read all meta-information until the header
    * (included) and stop before reading first data item. To open a file
    * use the open_file() method, and to pass an open stream to VCF data use
    * set_stream() (take care in the latter case that the file must be open and
    * that the first line to be read must be the "fileformat"
    * specification line). Both methods import information from the
    * header. It is possible to pass a string containing the header only
    * with read_header() but then it will be necessary to passed all
    * lines separately with read_line().
    *
    * After reading the header, several items of information are made
    * available through methods: file format (string identifying the
    * format used to encode data), which is required, and optional
    * meta-information fields which can be multiple and are identified
    * by their ID: INFO (specifying information fields relative to each
    * variable position), FORMAT (specifying information fields relative
    * to each sample for a given variable position), FILTER (identifies
    * criteria used to filter variable positons) and ALT (identifies
    * pre-defined alternate alleles). The fileformat string is available
    * through the file_format() method. For INFO specifications,
    * num_info() gives the number of INFO specifications and info()
    * gives access to a given index (the equivalent methods exist for
    * FORMAT, FILTER and ALT). The accessors return dedicated classes.
    * If it also possible to use find_info() and equivalent methods who
    * look to a specification by its ID. Before parsing the header, the
    * object loads a number of pre-defined INFO, FORMAT and ALT
    * specifications as defined in the VCF 4.1 format definition. If the
    * file contains a specification matching an existing one (either
    * from pre-defined specifications or from the file itself), it
    * overwrites it. In addition, it is possible to use the methods
    * add_info() and similar to add user-specified definitions. Beware,
    * then, that specification loaded from the file might overwrite
    * them. All specifications are reset upon opening or setting a new
    * file. Meta-information lines that do not use the "fileformat",
    * "INFO", "FORMAT", "FILTER" and "ALT" keys fall into the default
    * "meta" category and can be accessed using the num_meta(), meta()
    * and find_meta() methods, and set/modified using add_meta() method.
    * Finally, the header line must define the number of samples (if
    * any). The number of samples and their names are accessible using
    * the methods num_samples() and sample().
    *
    * The method allow_X() switches support (be default no support) for
    * using ``X`` or ``x`` instead of a base in alternate alleles. This
    * is not allowed in VCF specification format 4.1 but some software
    * does actually use it.  If X is allowed and one is found, the
    * alternate type will be set to vcf::X (not vcf::Default) and the
    * corresponding allele string will be ``X`` (regardless of the
    * original case).
    *
    * The method allow_gap() switches support (be default no support)
    * for gap (``-``) as a valid base (allow both reference and
    * alternate alleles to either contain or be a gap symbol. This is
    * not allowed in VCF specification format 4.1 and its use is
    * discouraged.
    * 
    * Each call to read() or read_line() processes a single variant
    * position. Each further call invalidates data stored from the last
    * read operation. At any moment, the method bool() tells whether the
    * underlying file stream is good for reading, but does not guarantee
    * that the next read operation will succeed. If no data is left to
    * read, the method read() returns false. Many formatting errors will
    * be intercepted and result in a EggFormatError exception specifying
    * the line number and as much information as possible.After reading
    * a line, a number of information items are available. Note that
    * allocated memory is never freed unless explicitly requested,
    * therefore speeding up the processing of large files (this also
    * applies when several files are processed in a row).
    *
    * List of information available after reading a VCF line:
    *
    *   - Chromosome (or other molecule) name: chromosome().
    *
    *   - Position of the variant on the chromosome: position(). Note
    * that the first position is 0, and that the first telomere is
    * represented by the constant egglib::BEFORE.
    *
    *   - The list of IDs defined for this variant can be accessed using
    * num_ID() and ID(). The number of IDs can be 0, 1 or more.
    *
    *   - Reference allele: reference(). The reference allele is
    * represented by one or more bases A, C, G, T or N. The number of
    * bases is directly accessible throught len_reference().
    *
    *   - Alternate alleles: num_alternate(), alternate_type() and
    * alternate(). There must be at least one alternate allele. The
    * alternate alleles can be represented by different types (see the
    * documentation).
    *
    *   - Variant quality score: quality().
    *
    *   - The list of failed tests can be analyzed using
    * num_failed_tests() and failed_test(). If all tests passed, the
    * number of failed tests is 0. If no tests were performed, the
    * number of failed tests is set to egglib::UNKNOWN (which is a very
    * large positive value).
    *
    *   - An arbitrary number (including none) of INFO fields can be
    * available. These INFO are separated between Flag, Integer, Float,
    * Character and String. The number of INFO items falling into each
    * category and each item can be accessed using num_FlagInfo() and
    * FlagInfo() methods, respectively (and equivalent for other types).
    * The types used are FlagInfo, TypeInfo<int> (for Integer),
    * TypeInfo<double> (for Float), TypeInfo<char> (for Character)
    * and StringInfo for String. If a non-specified INFO field is
    * used, its type is set to String.
    *
    *   - If the INFO fields AN, AC, AF and AA are defined and match
    * if their definition is conform to the standard definitions, their
    * value is directly accessible using dedicated members. The booleans
    * has_AN(), has_AC(), has_AF() and has_AA() allow to test if data
    * are available. The counters num_AC() and num_AF() return the
    * number of entries (which must be equal to the number of alternate
    * alleles). The value are accessible through AN(), AC(), AF() and
    * AA().
    *
    *   - If, and only if, more than 0 samples are defined,
    * sample-specific description fields are available. The fields are
    * described by IDs that normally are defined in the header or in
    * pre-defined types (as FORMAT specification). Undefined types are
    * not allowed. The methods num_field() and field() allow to explore
    * the IDs of FORMAT fields used. As for INFO fields, they are sorted
    * by type (except that there is not Flag type for FORMAT
    * specifications). To get the index of a FORMAT field amongst fields
    * of its types, use the overloaded methods field_rank(). The method
    * sample_info() provides an object of the class SampleInfo that
    * contains all FORMAT fields for a given sample (identified by its
    * index). The type, and the type-specific index returned by
    * field_rank(), are required to get the value corresponding to a
    * given FORMAT specification.
    *
    * Header: <egglib-cpp/VCF.hpp>
    *
    */
    class VcfParser {

      public:

       /** \brief Constructor
        *
        * The constructor does not generate an object ready for use.
        * Call to open or set methods is needed before starting to
        * parse data. The constructor automatically imports a set of
        * pre-defined INFO and FORMAT specification as specified by the
        * format standard.
        *
        */
        VcfParser();

       /** \brief Destructor
        *
        */
        virtual ~VcfParser();

       /** \brief Open a file for reading
        *
        * This method attempts to open the specified file and to read a
        * the VCF header. If the file cannot be open, an
        * EggOpenFileError exception is thrown; if the header is invalid
        * an EggFormatError exception is thrown.
        *
        * In case the instance was already processing a stream, it will
        * be dismissed. The stream created by this method will be closed
        * if another stream is created or set (call to open_file() or set_stream()
        * methods), if the close() method is called or upon object
        * destruction.
        *
        * \param fname name of the VCF-formatted file to open.
        *
        */
        void open_file(const char* fname);

       /** \brief Pass an open stream for reading
        *
        * This method set the passed stream (which is supposed to have
        * been opened for reading). If the stream is not good for
        * reading, an EggArgumentValueError (and not EggOpenFileError)
        * exception is thrown.
        *
        * In case the instance was already processing a stream, it will
        * be dismissed. The stream passed by this method not be closed
        * by the class even when calling close().
        *
        * \param stream open stream to read VCF data from.
        */
        void set_stream(std::istream& stream);

       /** \brief Pass header string for reading
        *
        * In case the instance was already processing a stream, it will
        * be dismissed. This function opens no stream, and it will not
        * be able read any further line using using read_line().
        *
        * \param string string containing the header.
        */
        void read_header(const char * string);

       /** \brief Switch support for ``X`` as an alternate allele
        *
        * This call will affect all subsequent reading operations but
        * the default value will be restored at the next call to set_stream(),
        * open_file(), read_header(), or reset(). The default is ``false``.
        *
        */
        void allow_X(bool flag);

       /** \brief Switch support for ``-`` as a valid base in reference/alternate alleles
        *
        * This call will affect all subsequent reading operations but
        * the default value will be restored at the next call to set_stream(),
        * open_file(), read_header(), or reset(). The default is ``false``.
        *
        */
        void allow_gap(bool flag);

       /** \brief Read a single variant
        *
        */
        void read();

       /** \brief Unread a single variant
        *
        *    Don't call this if no single variant has been read.
        */
        void unread();

       /** \brief set the attribut previous_next 
        *
        * Don't call this if no stream set
        */
        void set_previous_next();

       /** \brief set position in stream 
        *
        * Don't call this if no stream set
        */
        void set_index(long long int index);

       /** \brief set _currline variable 
        *
        */
        void set_currline(unsigned int line);

       /** \brief The current value of the flags is overwritten
        */
        void clear_stream();

       /** \brief get current position in stream
        *
        * Don't call this if no stream set
        */
        long long int get_index() const;

       /** \brief Read a single variant from a provided string
        *
        */
        void read_line(const char *);

       /** \brief Check if object is ready to parse
        *
        * \return true if the object has a valid stream and the stream
        * is ready to parse data and not end of file.
        *
        */
        bool good() const;

        bool has_data() const; ///< True if any data has been read

       /** \brief Reset this objet
        *
        * This method closes the file that was opened using the
        * open_file() method. If the file was open using the open_file() method of
        * the same instance, it is actually closed. If the file was
        * passed as a stream using set_stream(), it is forgotten but not
        * closed. If no stream is present, this method does nothing.
        *
        */
        void reset();

       /** \brief Forget information from last read variant
        *
        * It is not necessary to call this method before calling read(),
        * even if variants have been read previously.
        *
        */
        void reset_variant();

       /** \brief Actually clears the memory of the instance
        *
        * Actually frees the memory of the instance. This method must
        * not be used while reading a file.
        *
        */
        void clear();

       /** \brief File format string of the current file
        *
        * By default, if no VCF file has been set, returns an empty but
        * valid string.
        *
        */
        const char * file_format() const;

       /** \brief Get the number of filter entries of the instance
        *
        */
        unsigned int num_filter() const;

       /** \brief Get a specific filter entry
        *
        */
        const vcf::Filter * get_filter(unsigned int i) const;

       /** \brief Get the number of info entries of the instance
        *
        */
        unsigned int num_info() const;

       /** \brief Get a specific info entry
        *
        */
        const vcf::Info * get_info(unsigned int i) const;

       /** \brief Get the number of format entries of the instance
        *
        */
        unsigned int num_format() const;

       /** \brief Get a specific format entry
        *
        */
        const vcf::Format * get_format(unsigned int i) const;

       /** \brief Get the number of meta-information entries of the instance
        *
        * This excludes any FILTER, INFO, FORMAT, ALT and the fileformat.
        *
        */
        unsigned int num_meta() const;

       /** \brief Get a specific meta-information entry
        *
        */
        const vcf::Meta * get_meta(unsigned int i) const;

       /** \brief Get the number of alternative allele entries of the instance
        *
        */
        unsigned int num_alt() const;

       /** \brief Get a specific alternative allele entry
        *
        */
        const vcf::Alt * get_alt(unsigned int i) const;

       /** \brief Add a filter entry
        *
        * If a filter with this ID already exists, it will be
        * overwritten.
        *
        */
        void add_filter(const char * id, const char * descr);

       /** \brief Add an alternative allele entry
        *
        * If an alternative allele with this ID already exists, it will
        * be overwritten.
        *
        */
        void add_alt(const char * id, const char * descr);

       /** \brief Add an info entry
        *
        * If an info allele with this ID already exists, it will be
        * overwritten.
        *
        */
        void add_info(const char * id, unsigned int num, vcf::Info::Type type, const char * descr);

       /** \brief Add a format entry
        *
        * Flag is not permitted as a type.
        *
        * If a format with this ID already exists, it will be
        * overwritten.
        *
        */
        void add_format(const char * id, unsigned int num, vcf::Info::Type type, const char * descr);

       /** \brief Add a meta-information entry
        *
        * If a meta-information with this key already exists, it will be
        * overwritten.
        *
        */
        void add_meta(const char * key, const char * val);

       /** \brief Get the number of samples read from the header */
        unsigned int num_samples() const;

       /** \brief Get a sample name */
        const char * get_sample(unsigned int i) const;

       /** \brief Find a filter specification
        *
        * Return the address of the vcf::Filter with the specified ID. If
        * none is found, return NULL.
        *
        */
        vcf::Filter * find_filter(const char * id);

       /** \brief Find a format specification
        *
        * Return the address of the vcf::Format with the specified ID. If
        * none is found, return NULL.
        *
        */
        vcf::Format * find_format(const char * id);

       /** \brief Find an info specification
        *
        * Return the address of the vcf::Info with the specified ID. If
        * none is found, return NULL.
        *
        */
        vcf::Info * find_info(const char * id);

       /** \brief Find a meta-information specification
        *
        * Return the address of the vcf::Meta with the specified key. If
        * none is found, return NULL.
        *
        */
        vcf::Meta * find_meta(const char * key);

       /** \brief Find an alternate allele specification
        *
        * Return the address of the vcf::Alt with the specified ID. If
        * none is found, return NULL.
        *
        */
        vcf::Alt * find_alt(const char * id);

       /** \brief Get the last read record chromosome
        *
        */
        const char * chromosome() const;

       /** \brief Get the last read variant position
        *
        * If the position was the one before first, the constant
        * BEFORE is returned.
        *
        */
        unsigned long position() const;

       /** \brief Get the number of IDs of the last read variant
        *
        */
        unsigned int num_ID() const;

       /** \brief Get an ID from the last read record
        *
        */
        const char * ID(unsigned int i) const;

       /** \brief Get the length of the reference allele from the last read variant
        *
        * The default is 0.
        *
        */
        unsigned int len_reference() const;

       /** \brief Get the reference allele from the last read variant
        *
        * The default is an empty string.
        *
        */
        const char * reference() const;

       /** \brief Get the number of alternate alleles for the last read variant
        *
        * The value is 0 is no variants were provided.
        *
        */
        unsigned int num_alternate() const;

       /** \brief Get the type of one of the alternate alleles for the last read variant
        *
        */
        vcf::AltType alternate_type(unsigned int i) const;

       /** \brief Get one of the alternate alleles for the last read variant
        *
        * If alternate_type(i) is Default, the returned string is the
        * allele itself. If alternate_type(i) is Referred, it is its
        * ID, without angle brackets `< >` (it must be present in the
        * meta-information). If alternate_type(i) is Breakend, then it
        * is a breakend specification, reproduced as is.
        *
        */
        const char * alternate(unsigned int i) const;

       /** \brief Get the quality of the last read variant
        *
        * Returns the phred-scaled quality of the variant (or no
        * variation), or egglib::UNDEF in case of missing data. Note
        * that UNDEF is a large negative value.
        *
        */
        double quality() const;

       /** \brief Get the number of failed tests for the last read variants
        *
        * The value is 0 is all tests passed, and egglib::UNKNOWN if no
        * tests were performed (missing value in file).
        *
        */
        unsigned int num_failed_tests() const;

       /** \brief Get the ID of one of the failed tests for the last read variant
        *
        * If alternate_type(i) is Default, the returned string is the
        * allele itself. If alternate_type(i) is Referred, it is its
        * ID, without angle brackets `< >` (it must be present in the
        * meta-information). If alternate_type(i) is Breakend, then it
        * is a breakend specification, reproduced as is.
        *
        * Warning, do not iterate over the value returned by
        * num_failed_tests() without checking that it is not equal to
        * egglib::UNKNOWN.
        *
        */
        const char * failed_test(unsigned int i) const;

       /** \brief Get the number of Flag-type INFO entries for the last read variant
        *
        */
        unsigned int num_FlagInfo() const;

       /** \brief Get a Flag-type INFO entry for the last read variant
        *
        */
        const vcf::FlagInfo FlagInfo(unsigned int) const;

       /** \brief Get the number of Integer-type INFO entries for the last read variant
        *
        */
        unsigned int num_IntegerInfo() const;

       /** \brief Get a Integer-type INFO entry for the last read variant
        *
        */
        const vcf::TypeInfo<int>& IntegerInfo(unsigned int) const;

       /** \brief Get the number of Float-type INFO entries for the last read variant
        *
        */
        unsigned int num_FloatInfo() const;

       /** \brief Get a Float-type INFO entry for the last read variant 
        *
        */
        const vcf::TypeInfo<double>& FloatInfo(unsigned int) const;

       /** \brief Get the number of Character-type INFO entries for the last read variant
        *
        */
        unsigned int num_CharacterInfo() const;

       /** \brief Get a Character-type INFO entry for the last read variant
        *
        */
        const vcf::TypeInfo<char>& CharacterInfo(unsigned int) const;

       /** \brief Get the number of String-type INFO entries for the last read variant
        *
        */
        unsigned int num_StringInfo() const;

       /** \brief Get a String-type INFO entry for the last read variant
        *
        */
        const vcf::StringInfo& StringInfo(unsigned int) const;

       /** \brief Get the number of sample FORMAT fields for the last read variant
        *
        */
        unsigned int num_fields() const;

       /** \brief Get a sample FORMAT field specification for the last read variant
        *
        * Return the corresponding FORMAT specification as a reference
        * to a vcf::Format instance.
        *
        */
        const vcf::Format& field(unsigned int i) const;

       /** \brief Get the index of a sample FORMAT field specification by its ID for the last read variant
        *
        * This is equivalent to looping over the field(unsigned int)
        * method until the type() method of the returned value matches
        * the specified ID. Returns egglib::UNKNOWN if the ID is not
        * found.
        *
        */
        unsigned int field_index(const char * ID) const;

       /** \brief Get the rank of a sample FORMAT field for the last read variant
        *
        * The index returned by this method gives the rank of the
        * corresponding FORMAT field among fields of the same type. See
        * sample_info(unsigned int) to understand why it is useful. This
        * methods is faster than using field_rank(const char *).
        *
        */
        unsigned int field_rank(unsigned int i) const;

       /** \brief Get the rank of a sample FORMAT field for the last read variant using its ID
        *
        * The index returned by this method gives the rank of the
        * corresponding FORMAT field among fields of the same type. See
        * sample_info(unsigned int) to understand why it is useful. The
        * behaviour is not defined if the ID is not found. This method
        * performs a search operation, and using field_rank(unsigned
        * int) is faster.
        *
        */
        unsigned int field_rank(const char * ID) const;

       /** \brief Get all FORMAT fields for a sample for the last read variant
        *
        * The returned instance contains all FORMAT fields for the
        * specified sample. The class vcf::SampleInfo provides methods to
        * access data which are based on the type and on the field rank
        * among of the fields of the same type.
        *
        * Assuming you know that the type of the data field you wish to
        * extract is, say, a String and its ID, you should first get
        * its type-based rank using field_rank(const char *) using its
        * ID and then call vcf::SampleInfo::StringItem(unsigned int,
        * unsigned int) to get the value of an item. Take care that
        * missing values result in fields with 0 items.
        *
        */
        const vcf::SampleInfo& sample_info(unsigned int i) const;

        /** \brief Check if the INFO field AN (number of called alleles) is available
         *
         * The method returns false if the AN field is not defined for
         * the last variant, or if its definition does not match the
         * standard.
         *
         */
         bool has_AN() const;

        /** \brief Get the AN value if defined
         *
         * If AN is not available, returns an undefined value.
         *
         */
         unsigned int AN() const;

        /** \brief Check if the INFO field AA (ancestral allele) is available
         *
         * The method returns false if the AA field is not defined for
         * the last variant or if its definition does not match the
         * standard.
         *
         */
         bool has_AA() const;

        /** \brief Get the index of the allele given as AA if defined
         *
         * If AA is not available, returns an undefined value. If AA is
         * available but the ancestral allele is not determined
         * (provided as a missing value), returns UNKNOWN. If AA is
         * available but not one of the alleles for this site, returns
         * the next valid index (num_alternate + 1).
         *
         */
         unsigned int AA_index() const;

        /** \brief Get the INFO field AA if available
         *
         * If AA is not available, returns an undefined value. If AA is
         * available but the ancestral allele is not determined
         * (provided as a missing value), return "?".
         *
         */
         const char * AA_string() const;

        /** \brief Check if the INFO field AC (allele absolute frequencies) is available
         *
         * The method returns false if the AC field is not defined for
         * the last variant, or if its definition does not match the
         * expectation, or if AC is not available.
         *
         */
         bool has_AC() const;

        /** \brief Get the number of AN values, if defined
         *
         * The number of AC values is equal to the number of alternate
         * alleles. If AC is not available, returns an undefined value.
         *
         */
         unsigned int num_AC() const;

        /** \brief Get an AC value if defined
         *
         * If AC is not available, or if the index is over the value
         * returned by num_AC(), this method might cause a crash.
         *
         */
         unsigned int AC(unsigned int i) const;

        /** \brief Check if the INFO field AF (allele relative frequencies) is available
         *
         * The method returns false if the AF field is not defined for
         * the last variant, or if its definition does not match the
         * expectation.
         *
         */
         bool has_AF() const;

        /** \brief Get the number of AF values, if defined
         *
         * The number of AF values is equal to the number of alternate
         * alleles. If AF is not available, returns an undefined value.
         *
         */
         unsigned int num_AF() const;

        /** \brief Get an AF value if defined
         *
         * If AF is not available, or if the index is over the value
         * returned by num_AF(), this method might cause a crash.
         *
         */
         double AF(unsigned int i) const;

        /** \brief Check if the FORMAT field GT (genotype of each sample) is available
         *
         * The method returns false if the GT field is not defined for
         * the last variant, or if its definition does not match the
         * expectation.
         *
         */
         bool has_GT() const;

        /** \brief Returns the ploidy of the last variant
         *
         * The ploidy must be a strictly positive number. If has_GT()
         * returns false, the returned value is 2.
         *
         */
         unsigned int ploidy() const;

        /// \brief Number of genotypes
         unsigned int num_genotypes() const;

        /** \brief Check if the genotype of a given sample is phased
         *
         * If has_GT() returns false, or if ploidy() returns 1, the
         * returned value is undefined.
         *
         */
         bool GT_phased(unsigned int i) const;

        /** \brief Check if the genotype of all samples are phased
         *
         * If has_GT() returns false, or if ploidy() returns 1, the
         * returned value is undefined. If the number of samples is 0,
         * returns true.
         *
         */
         bool GT_phased() const;

        /** \brief Get a genotype value
         *
         * If has_GT() returns false, the behaviour of this method is
         * undefined.
         *
         * \param sample sample index (must be < num_samples()).
         * \param allele allele, or chromosome, index (must be <
         * ploidy()).
         * \return The index of the allele carried by this sample at
         * this ploidy index, or egglib::UNKNOWN if value is missing. It
         * is NOT guaranteed that if one allele is missing, all are. The
         * method returns 0 if the carried allele is the reference.
         * Otherwise, it returns 1 + the index of the allele within the
         * list given by alternate alleles.
         *
         */
         unsigned int GT(unsigned int sample, unsigned int allele) const;

        /// \brief Get a PL value
         unsigned int PL(unsigned int sample, unsigned int genotype) const;
        
        /// \brief Get a GL value
         double GL(unsigned int sample, unsigned int genotype) const ;

        /** \brief Check if the Variant read has PL data.
         *
         */
         bool has_PL() const;

        /** \brief Check if the Variant read has GL data.
         *
         */
         bool has_GL() const;

        /** \brief Check the index of the EOF of the VcfFile loaded
         *   \return the stream position of the end of file.
         */
         long long int file_end();

        /** \brief this method allows to find, read and load as an VcfIndex an Index
         *       file linked to the current Vcf Parser.  The search of the index file, 
         *       is done with a default name generated from the name of the Vcf file 
         *       loaded in the current VcfParser.
         *   \param fname path of an index file
         *  
         */
         void read_index(const char *fname);

        /** \brief checks if the VcfParser has an index file loaded in a object VcfIndex.
         */
         bool has_index();

        /** \brief gets a default name for an index file according the name of a Vcf File 
         * loaded in the VcfParser
         */
        const char * index_fname();

        /** \brief gets the index of the first variant in the VcfParser
         *
         */
        long long int get_idx_frt_sample();

        /** \brief this method allows to move the VcfParser to the first variant of the Vcffile loaded. 
         *
         */
        void rewind();

        /// \brief get only chromosome name
        void read_chromosome();

        unsigned int num_missing_GT();

        /// \brief Set the threshold for extracting GT from PL (UNKNOWN to prevent, otherwise 1 or more, never 0)
        void set_threshold_PL(unsigned int);

        /// \brief Get the threshold used for extracting GT from PL
        unsigned int get_threshold_PL() const;


        /// \brief Set the threshold for extracting GT from PL (UNKNOWN to prevent, otherwise 1 or more, never 0)
        void set_threshold_GL(unsigned int);

        /// \brief Get the threshold used for extracting GT from PL
        unsigned int get_threshold_GL() const;

        bool is_outgroup(unsigned int) const; ///< Tells if a sample is set to be outgroup
        void set_outgroup(unsigned int); ///< Set a sample to be outgroup
        bool outgroup_AA() const; ///< Tells if AA data should be picked to generate sites
        void set_outgroup_AA(); ///< Tells to add AA data as first sample of outgroup when sites are generated

      private:

       /** \brief Copy constructor is disabled */
        VcfParser(const VcfParser& src) {}

       /** \brief Copy assignment operator is disabled */
        VcfParser& operator=(const VcfParser& src) { return *this; }

        void init();
        void free();
        void header();
        void _check_eof();

        unsigned int get_string(char *& where, unsigned int& _res_, bool (VcfParser::* check)(unsigned int), bool (VcfParser::* stop)(), bool catch_missing);

        // EOF supported wherever \n is allowed!
        bool stop_equal();                      // stop at `=`, error if `\r` or `\n`
        bool stop_equalsemicolontablineEOF();   // stop at `=`, `;`, `\t`, `\n` or EOF, expect `n` if `\r`
        bool stop_colon();                      // stop at ':', error if `\r` or `\n`
        bool stop_tab();                        // stop at `\t`, error if `\r` or `\n`
        bool stop_tabcolon();                   // stop at `\t` or ':', error if `\r` or `\n`
        bool stop_tabsemicolon();               // stop at `\t` or ';', error if `\r` or `\n`
        bool stop_tabcomma();                   // stop at `\t` or ',', error if `\r` or `\n`
        bool stop_tabsemicoloncomma();          // stop at `\t` or ';', error if `\r` or `\n`
        bool stop_tabsemicoloncommalineEOF();   // stop at `\t`, `,`, `;`, `\n` or EOF, expect `\n` if `\r`
        bool stop_tabcoloncommalineEOF();       // stop at `\t`, `,`, `:`, `\n` or EOF, expect `\n` if `\r`
        bool stop_line();                       // stop at `\n`, expect `\n` if `\r`
        bool stop_linetabEOF();                 // stop at `\n`, EOF or '\t', expect `\n` if `\r`
        bool stop_field();                      // stop at `,` or `>`, error if `\r or `\n`
        bool stop_quote();                      // stop at `"` if previous is not `\`, error if `\r` or `\n`

        bool check_letter(unsigned int i);                    // `A` to `Z` or `a` to `z`
        bool check_sign(unsigned int i);                      // all printing characters (`!` to `~`)
        bool check_sign_and_space(unsigned int i);            // all printing characters (`!` to `~`) plus space
        bool check_string(unsigned int i);                    // all printing characters (`!` to `~`) plus space and tab
        bool check_alphanumericunderscore(unsigned int i);    // `0` to `9` or `A` to `Z` or `a` to `z` + the signs in `_+-.*`
        bool check_integer(unsigned int i);                   // `0` to `9`
        bool check_float(unsigned int i);                     // `0` to `9` or `-` or `E` or `e` or `.`
        bool check_bases(unsigned int i);                     // `A`, `C`, `G`, `T` or `N` (and force upper)

        vcf::Filter * add_filter();
        vcf::Info   * add_info();
        vcf::Format * add_format();
        vcf::Meta   * add_meta();
        vcf::Alt    * add_alt();

        vcf::FlagInfo         * add_FlagInfo(unsigned int expected_number);
        vcf::TypeInfo<char>   * add_CharacterInfo(unsigned int expected_number);
        vcf::TypeInfo<int>    * add_IntegerInfo(unsigned int expected_number);
        vcf::TypeInfo<double> * add_FloatInfo(unsigned int expected_number);
        vcf::StringInfo       * add_StringInfo(unsigned int expected_number);

        void predefine();

        void next();
        void get_4_fields(bool format);
        void get_2_fields(bool alt);

        void _find_alleles(unsigned int na, unsigned int *, unsigned int idx);
        unsigned int ** _genotype_idx_helper;
        unsigned int _genotype_idx_helper_size;

        std::ifstream _localstream;
        std::istringstream _stringstream; // used only within functions (not persistent)
        std::istream* _stream;
        unsigned long _previous_next;
        unsigned int _currline;
        unsigned int _start_currline;
        long long int  _idx_frt_sample;
        long long int _idx_end_file;
        unsigned int _res_fname;
        char * _fname;
        unsigned int _threshold_PL;
        unsigned int _threshold_GL;

        double * _buffer_float;
        unsigned int _res_buffer;
        char * _buffer;
        unsigned int _res_buffer2;
        char * _buffer2;

        char curr_ch;
        char prev_ch;

        unsigned int _res_ff;
        char * _ff;

        unsigned int _num_filter;
        unsigned int _num_info;
        unsigned int _num_format;
        unsigned int _num_meta;
        unsigned int _num_alt;
        unsigned int _res_filter;
        unsigned int _res_info;
        unsigned int _res_format;
        unsigned int _res_meta;
        unsigned int _res_alt;
        vcf::Filter ** _filter;
        vcf::Info   ** _info;
        vcf::Format ** _format;
        vcf::Meta   ** _meta;
        vcf::Alt    ** _alt;

        unsigned int _num_samples;
        unsigned int _res_samples;
        unsigned int * _res_len_samples;
        char ** _samples;
        vcf::SampleInfo ** _sampleInfo;
        bool * _outgroup_map;
        bool _outgroup_AA;

        unsigned int _res_chrom;
        char * _chrom;
        unsigned long _position;
        unsigned int _num_ID;
        unsigned int _res_ID;
        unsigned int * _res_len_ID;
        char ** _ID;
        unsigned int _len_reference;
        unsigned int _res_reference;
        char * _reference;
        unsigned int _num_alternate;
        unsigned int _res_alternate;
        vcf::AltType * _type_alternate;
        unsigned int * _res_len_alternate;
        char ** _alternate;
        unsigned int _ploidy;
        unsigned int _num_genotypes;
        double _quality;
        unsigned int _num_failed_test;
        unsigned int _res_failed_test;
        unsigned int * _res_len_failed_test;
        char ** _failed_test;

        unsigned int _num_FlagInfo;
        unsigned int _res_FlagInfo;
        vcf::FlagInfo ** _FlagInfo;

        unsigned int _num_CharacterInfo;
        unsigned int _res_CharacterInfo;
        vcf::TypeInfo<char> ** _CharacterInfo;

        unsigned int _num_IntegerInfo;
        unsigned int _res_IntegerInfo;
        vcf::TypeInfo<int> ** _IntegerInfo;

        unsigned int _num_FloatInfo;
        unsigned int _res_FloatInfo;
        vcf::TypeInfo<double> ** _FloatInfo;

        unsigned int _num_StringInfo;
        unsigned int _res_StringInfo;
        vcf::StringInfo ** _StringInfo;

        unsigned int _num_formatEntries;
        unsigned int _res_formatEntries;
        vcf::Format ** _formatEntries; // addresses not freeable
        unsigned int * _formatRank;

        bool _has_data;

        bool _has_AN;
        unsigned int _AN;

        bool _has_AC;
        bool _has_AC_ss;
        unsigned int _num_AC;
        unsigned int _res_AC;
        unsigned int * _AC;

        bool _has_AF;
        unsigned int _num_AF;
        unsigned int _res_AF;
        double * _AF;

        bool _has_AA;
        unsigned int _AA_index;
        const char * _AA_string;
        char * _AA_missing;

        bool _has_GT;
        unsigned int * _res_GT; // size: _num_samples
        unsigned int * _res_PL; // size: _num_samples
        unsigned int * _res_GL; // size: _num_samples

        unsigned int ** _GT; // size : _num_samples * _res_GT[i], the values are UNKNOWN, _num_alternate (if _reference), or an alternate allele index (not as in VCF)
        bool * _GT_phased;  // size: _num_samples
        bool _GT_all_phased;

        bool _has_PL;
        unsigned int ** _PL;  // size: _num_samples * res_PL[i]


        bool _has_GL;
        double ** _GL;  // size: _num_samples * res_GL[i]

        bool _allow_X;
        bool _allow_gap;
        bool _has_index;
        char * _index;
    };

   /* \brief gets the EOF position of an VCF file.
    *  \param fname path of a VCF file.
    */
    long long int file_end( const char * fname);
}

#endif
