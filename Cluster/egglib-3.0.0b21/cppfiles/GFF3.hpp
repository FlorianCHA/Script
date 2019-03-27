/*
    Copyright 2012,2014-2015 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_GFF3_HPP
#define EGGLIB_GFF3_HPP

#include <istream>
#include "DataHolder.hpp"

namespace egglib {

    class GFF3;

   /** \brief Annotation feature
    *
    * \ingroup parsers
    *
    * Objects of this class are used to describe annotation features in
    * the %GFF3 format read by the GFF3 class.
    *
    * Header: <egglib-cpp/GFF3.hpp>
    *
    */
    class Feature {

      public:

        friend class GFF3;

       /** \brief Enum for reading frame specification */
        enum PHASE {
            zero, ///< Codon starts at first base
            one, ///< Codon starts at second base
            two, ///< Codon starts at third base
            no_phase ///< Not defined (irrelevant)
        };

       /** \brief Enum for strand specification */
        enum STRAND {
            plus, ///< Forward strand
            minus, ///< Reverse strand
            no_strand ///< Not defined (irrelevant)
        };

       /** \brief Constructor */
        Feature();

       /** \brief Copy constructor
        *
        * Note that "parent" and "parts" members (which are links to
        * other Feature objects) are shallow-copied. As a result, the
        * copy object will point to the same parents and/or parts as the
        * original.
        *
        */
        Feature(const Feature& src);

       /** \brief Copy assignment operator
        *
        * Note that "parent" and "parts" members (which are links to
        * other Feature objects) are shallow-copied. As a result, the
        * copy object will point to the same parents and/or parts as the
        * original.
        *
        */
        Feature& operator=(const Feature& src);

       /** \brief Destructor */
        virtual ~Feature();

       /** \brief Reset instance (but retain allocated memory) */
        void reset();

       /** \brief Actually release memory of the instance*/
        void clear();

       /** \brief Get "seqid" field (empty string by default) */
        const char * get_seqid() const;

       /** \brief Set "seqid" field */
        void set_seqid(const char * str);

       /** \brief Get "source" field (empty string by default)*/
        const char * get_source() const;

       /** \brief Set "source" field */
        void set_source(const char * str);

       /** \brief Get "type" field (empty string by default) */
        const char * get_type() const;

       /** \brief Set "type" field */
        void set_type(const char * str);

       /** \brief Get number of fragments (number of "start" and "end" fields; 0 by default) */
        unsigned int get_num_fragments() const;

       /** \brief Set number of fragments (number of "start" and "end" fields)
        *
        * If the value is larger, new values are initialized to 0 (both
        * "start" and "end"). If the new value is smaller, last values
        * are lost.
        *
        */
        void set_num_fragments(unsigned int num);

       /** \brief Get "start" for a given fragment */
        unsigned int get_start(unsigned int i) const;

       /** \brief Set "start" for a given fragment */
        void set_start(unsigned int i, unsigned int val);

       /** \brief Get "end"  for a given fragment */
        unsigned int get_end(unsigned int i) const;

       /** \brief Set "end" for a given fragment */
        void set_end(unsigned int i, unsigned int val);

       /** \brief Get "score" field (UNDEF by default)
        *
        * egglib::UNDEF stands for undetermined. It equals to a very
        * large negative value. If the file had a "NaN" (or any case
        * combination), UNDEF is also returned.
        *
        */
        double get_score() const;

       /** \brief Set "score" field */
        void set_score(double d);

       /** \brief Get "strand" field (no_strand by default) */
        STRAND get_strand() const;

       /** \brief Set "strand" field */
        void set_strand(STRAND s);

       /** \brief Get "phase" field (no_phase by default) */
        PHASE get_phase() const;

       /** \brief Set "phase" field */
        void set_phase(PHASE p);

       /** \brief Get "ID" attribute (empty string if missing) */
        const char * get_ID() const;

       /** \brief Set "ID" attribute (use an empty string to skip) */
        void set_ID(const char * str);

       /** \brief Get "Name" attribute (empty string if missing) */
        const char * get_Name() const;

       /** \brief Set "Name" attribute (use an empty string to skip) */
        void set_Name(const char * str);

       /** \brief Get the number of "Alias" attributes */
        unsigned int get_num_Alias() const;

       /** \brief Set the number of "Alias" attributes
        *
        * If the new value is larger, new attributes are set to empty
        * strings. If the new value is smaller, last attributes are
        * lost.
        *
        */
        void set_num_Alias(unsigned int num);

       /** \brief Get the "Alias" attribute at the specified index */
        const char * get_Alias(unsigned int i) const;

       /** \brief Set the "Alias" attribute at the specified index */
        void set_Alias(unsigned int i, const char * s);

       /** \brief Get the number of "Parent" attributes */
        unsigned int get_num_Parent() const;

       /** \brief Set the number of "Parent" attributes
        *
        * If the new value is larger, new attributes are set to empty
        * strings. If the new value is smaller, last attributes are
        * lost.
        *
        */
        void set_num_Parent(unsigned int num);

       /** \brief Get the "Parent" attribute at the specified index */
        const char * get_Parent(unsigned int i) const;

       /** \brief Set the "Parent" attribute at the specified index */
        void set_Parent(unsigned int i, const char * s);

       /** \brief Get "Target" attribute (empty string if missing) */
        const char * get_Target() const;

       /** \brief Set "Target" attribute (use an empty string to skip) */
        void set_Target(const char * str);

       /** \brief Get "Gap" attribute (empty string if missing) */
        const char * get_Gap() const;

       /** \brief Set "Gap" attribute (use an empty string to skip) */
        void set_Gap(const char * str);

       /** \brief Get "Derives_from" attribute (empty string if missing) */
        const char * get_Derives_from() const;

       /** \brief Set "Derives_from" attribute (use an empty string to skip) */
        void set_Derives_from(const char * str);

       /** \brief Get the number of "Note" attributes */
        unsigned int get_num_Note() const;

       /** \brief Set the number of "Note" attributes
        *
        * If the new value is larger, new attributes are set to empty
        * strings. If the new value is smaller, last attributes are
        * lost.
        *
        */
        void set_num_Note(unsigned int num);

       /** \brief Get the "Note" attribute at the specified index */
        const char * get_Note(unsigned int i) const;

       /** \brief Set the "Note" attribute at the specified index */
        void set_Note(unsigned int i, const char * s);

       /** \brief Get the number of "Dbxref" attributes */
        unsigned int get_num_Dbxref() const;

       /** \brief Set the number of "Dbxref" attributes
        *
        * If the new value is larger, new attributes are set to empty
        * strings. If the new value is smaller, last attributes are
        * lost.
        *
        */
        void set_num_Dbxref(unsigned int num);

       /** \brief Get the "Dbxref" attribute at the specified index */
        const char * get_Dbxref(unsigned int i) const;

       /** \brief Set the "Dbxref" attribute at the specified index */
        void set_Dbxref(unsigned int i, const char * s);

       /** \brief Get the number of "Ontology_term" attributes */
        unsigned int get_num_Ontology_term() const;

       /** \brief Set the number of "Ontology_term" attributes
        *
        * If the new value is larger, new attributes are set to empty
        * strings. If the new value is smaller, last attributes are
        * lost.
        *
        */
        void set_num_Ontology_term(unsigned int num);

       /** \brief Get the "Ontology_term" attribute at the specified index */
        const char * get_Ontology_term(unsigned int i) const;

       /** \brief Set the "Ontology_term" attribute at the specified index */
        void set_Ontology_term(unsigned int i, const char * s);

       /** \brief Get "Is_circular" attribute (true if present) */
        bool get_Is_circular() const;

       /** \brief Set "Is_circular" attribute (false to skip) */
        void set_Is_circular(bool b);

       /** \brief Number of custom attributes
        *
        * Does not take into account pre-defined attributes which are
        * identified by a first capital letter, and which are accessible
        * through specific methods.
        *
        */
        unsigned int get_num_attributes() const;

       /** \brief Set the number of attributes
        *
        * Does not take into account pre-defined attributes (which are
        * identified by a first capital letter). If the new value is
        * larger, new attributes are set to an empty string for key, and
        * a number of items of 0. If the new value is smaller, last
        * attributes are lost.
        *
        */
        void set_num_attributes(unsigned int num);

       /** \brief Number of items of a custom attribute
        *
        */
        unsigned int get_num_items_attribute(unsigned int i) const;

       /** \brief Set the number of items of a custom attribute
        *
        * If the new value is larger, new items are set to empty
        * strings. If the new value is smaller, last items are lost.
        *
        */
        void set_num_items_attribute(unsigned int i, unsigned int num);

       /** \brief Get a custom attribute key
        *
        */
        const char * get_attribute_key(unsigned int i) const;

       /** \brief Set a custom attribute key
        *
        * Legal keys don't start with a capital letter (otherwise you
        * must use one of the pre-defined attributes).
        *
        */
        void set_attribute_key(unsigned int i, const char * str);

       /** \brief Get the value of a custom attribute item
        *
        */
        const char * get_attribute_value(unsigned int attr, unsigned int item) const;

       /** \brief Set a custom attribute value
        *
        */
        void set_attribute_value(unsigned int attr, unsigned int item, const char * str);

       /** \brief Get the number of parents
        *
        */
        unsigned int get_num_parents() const;

       /** \brief Set the number of parents
        *
        * If the new value is larger, new items are set to NULL. If the
        * new value is smaller, last items are lost.
        *
        */
        void set_num_parents(unsigned int num);

       /** \brief Get a parent
        *
        */
        Feature * get_parent(unsigned int i) const;

       /** \brief Set a parent
        *
        */
        void set_parent(unsigned int i, Feature * feat);

       /** \brief Get the number of parts
        *
        */
        unsigned int get_num_parts() const;

       /** \brief Set the number of parts
        *
        * If the new value is larger, new items are set to NULL. If the
        * new value is smaller, last items are lost.
        *
        */
        void set_num_parts(unsigned int num);

       /** \brief Get a part
        *
        */
        Feature * get_part(unsigned int i) const;

       /** \brief Set a part
        *
        */
        void set_part(unsigned int i, Feature * feat);

      private:

        // utilities

        void init();
        void free();
        void copy(const Feature& src);

        // standard fields

        unsigned int _res_seqid;
        char * _seqid;
        unsigned int _res_source;
        char * _source;
        unsigned int _res_type;
        char * _type;
        unsigned int _num_frag;
        unsigned int _res_frag;
        unsigned int * _start;
        unsigned int * _end;
        double _score;
        STRAND _strand;
        PHASE _phase;

        // reserved attributes

        unsigned int _res_ID;
        char * _ID;

        unsigned int _res_Name;
        char * _Name;

        unsigned int _num_Alias;
        unsigned int _res_Alias;
        unsigned int * _res_len_Alias;
        char ** _Alias;

        unsigned int _num_Parent;
        unsigned int _res_Parent;
        unsigned int * _res_len_Parent;
        char ** _Parent;

        unsigned int _res_Target;
        char * _Target;

        unsigned int _res_Gap;
        char * _Gap;

        unsigned int _res_Derives_from;
        char * _Derives_from;

        unsigned int _num_Note;
        unsigned int _res_Note;
        unsigned int * _res_len_Note;
        char ** _Note;

        unsigned int _num_Dbxref;
        unsigned int _res_Dbxref;
        unsigned int * _res_len_Dbxref;
        char ** _Dbxref;

        unsigned int _num_Ontology_term;
        unsigned int _res_Ontology_term;
        unsigned int * _res_len_Ontology_term;
        char ** _Ontology_term;

        bool _Is_circular;

        // custom attributes

        unsigned int _num_attributes;
        unsigned int _res_attributes;
        unsigned int * _num_attributes_items;
        unsigned int * _res_attributes_items;
        unsigned int * _res_len_attributes_key;
        unsigned int ** _res_len_attributes_val;
        char ** _attributes_key;
        char *** _attributes_val;

        // links to parents / parts instances

        unsigned int _num_parents;
        unsigned int _res_parents;
        Feature ** _parents;

        unsigned int _num_parts;
        unsigned int _res_parts;
        Feature ** _parts;
    };


////////////////////////////////////////////////////////////////////////

   /** \brief %GFF3 parser
    *
    * \ingroup parsers
    *
    * Read %GFF3-formatted genome annotation data from a file specified
    * by name or from an open stream.
    *
    * The description of the %GFF3 format:
    *
    * http://www.sequenceontology.org/gff3.shtml
    *
    * This class supports segmented features but only if they are
    * consecutive in the file. All features are loaded into memory and
    * can be processed interatively. Two accesors are provided: one,
    * feature() and num_features(), allows to process all imported
    * features in the order in which they were loaded. Each provides
    * access to its own parents and parts; and the second, gene() and
    * num_genes(), allows to process the subset of the latter that are
    * of type gene.
    *
    * Header: <egglib-cpp/GFF3.hpp>
    *
    */
    class GFF3 {

      public:

       /** \brief Build an empty object */
        GFF3();

       /** \brief Destructor */
        virtual ~GFF3();

       /** \brief Set the liberal flag
        * 
        * The liberal feature allows a few violations of the format as
        * specified in the standard (more violations might be added in
        * the future:
        *     - CDS features may lack a phase (then the no_phase value
        *       is used.
        *
        * By default, the parsers are strict. The new value affects all
        * consecutive parsing operations.
        *
        */
        void liberal(bool flag);

       /** \brief Parse a GFF3-formatted file */
        void parse(const char* fname);

       /** \brief Parse an open GFF3-formatted stream */
        void parse(std::istream& stream);

       /** \brief Parse an open GFF3-formatted string */
        void parse_string(std::string& string);

       /** \brief Clear data stored in the object (but don't free memory) */
        void reset();

       /** \brief Like reset() but actually free memory */
        void clear();

       /** \brief Get the number of defined meta-data */
        unsigned int num_metadata() const;

       /** \brief Get a meta-data key */
        const char * metadata_key(unsigned int i) const;

       /** \brief Get a meta-data value */
        const char * metadata_value(unsigned int i) const;

       /** \brief Number of gene features */
        unsigned int num_genes() const;

       /** \brief Get a gene feature */
        Feature& gene(unsigned int i);

       /** \brief Total number of features */
        unsigned int num_features() const;

       /** \brief Get a feature */
        Feature& feature(unsigned int i);

       /** \brief Get sequences (in case they were present in file) as a non-matrix object */
        const DataHolder& sequences() const;

      private:

       /** \brief Copy constructor is disabled */
        GFF3(const GFF3& src) {}

       /** \brief Copy assignment operator is disabled */
        GFF3& operator=(const GFF3& src) { return *this; }

        // utilities

        void init();
        void free();
        void copy(const GFF3& src);

        // parsing method

        void _parse(std::istream& stream);
        void get_fasta();
        void get_directive();
        void skip_line();
        void get_annotation();
        void predefined_attribute(Feature * f, int& flag);
        void custom_attribute(Feature * f);
        void get_items(Feature * f, void (Feature::* num_items)(unsigned int), unsigned int *& res_len, char **& _items);

        // parsing helpers

        unsigned int get_string(char *& where, unsigned int& _res_, bool (GFF3::* check)(), bool (GFF3::* stop)(), bool skip_initial_spaces=false);

        bool stop_lineEOF();                  // stop at `\n` or EOF, expect `\n` if `\r`
        bool stop_equalsemicolonlineEOF();    // stop at `=`, `;`, `\n` or EOF, expect `\n` if `\r`, error if `\t` (for attribute key or flag)
        bool stop_semicolonlineEOF();         // stop at `;`, `\n` or EOF, expect `\n` if `\r`, error if `\t` (for unique item)
        bool stop_semicoloncommalineEOF();    // stop at `;`, `,`, `\n` or EOF, expect `\n` if `\r`, error if `\t` (for item or items)
        bool stop_tabspacelineEOF();          // stop at `\t` or ` ` or '\n' or EOF, expect `\n` if `\r`
        bool stop_tabspace();                 // stop at `\t` or ` `, error if `\n`, `\r` or EOF
        bool stop_tab();                      // stop at `\t`, error if `\n`, `\r` or EOF
        bool stop_equal();                    // stop at `=`, error if `\n`, `\r` or EOF

        bool check_string();                  // all printing characters (`!` to `~`)
        bool check_stringESC();               // all printing characters (`!` to `~`), supported HEX escape
        bool check_integer();                 // `0` to `9`
        bool check_float();                   // `0` to `9` or `-` or `E` or `e` or `.`, also tolerate `N` or `n` (if first or third assuming that second is `a` or `A`) and `A` or `a` if second and first was `A` and `a`

        // parsing members

        std::istream * _stream;
        char * _fname;   // "" while not parsing
        unsigned int _res_fname;
        unsigned int currline;
        char * buffer_ESC;
        char * buffer;
        unsigned int res_buffer;
        char curr_ch;
        unsigned int curr_pos;
        unsigned int mark; // first feature allowed to be completed
        bool _liberal;

        // data

        unsigned int _num_metadata;
        unsigned int _res_metadata;
        unsigned int * _res_len_metadata_key;
        unsigned int * _res_len_metadata_val;
        char ** _metadata_key;
        char ** _metadata_val;

        unsigned int _num_features;
        unsigned int _res_features;
        Feature ** _features;

        unsigned int _num_genes;
        unsigned int _res_genes;
        Feature ** _genes;      // do not delete objects

        DataHolder _sequences;
    };
}

#endif
