/*
    Copyright 2012-2013,2016 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_FILTER_HPP
#define EGGLIB_FILTER_HPP

namespace egglib {

   /** \brief Holds lists of valid (exploitable and missing) data codes
    *
    * \ingroup diversity
    *
    * This class holds two list of integer values: one corresponding to
    * data codes that should be treated as valid and exploitable, and
    * one corresponding to data codes that should be treated as valid
    * but missing. If the list of exploitable data is empty, all data
    * are considered to be exploitable (even those who are missing).
    *
    * Header: <egglib-cpp/Filter.hpp>
    *
    */
    class Filter {

        public:

            /** \brief Constructor */
             Filter();

            /** \brief Copy constructor */
             Filter(const Filter& src);

            /** \brief Copy assignment operator */
             Filter& operator=(const Filter& src);

            /** \brief Destructor */
             virtual ~Filter();

            /** \brief Reserve memory
             *
             * The method pre-allocates data arrays in order to speed up
             * subsequent loading operations (up to the numbers passed).
             * The instance is not formally changed by this method, and
             * it is absolutely not required to call this method prior
             * setting valid or missing data codes.
             *
             * \param num_expl expected number of exploitable data
             *        codes.
             * \param num_expl_ranges expected number of ranges of
             *        exploitable data codes.
             * \param num_missing expected number of missing data codes.
             * \param num_synonyms expected number of synonyms (it is
             * not required that all exploitable codes have a synonym).
             * \param num_missing_synonyms expected number of synonym
             * for missing data (it is not required that all missing
             * codes have a synonym).
             */
             void reserve(unsigned int num_expl, unsigned int
                          num_expl_ranges, unsigned int num_missing,
                          unsigned int num_synonyms, unsigned int
                          num_missing_synonyms);

            /** \brief Reset instance
             *
             * The instance is delivered as newly created, and memory
             * is actually released.
             *
             */
             void clear();

            /** \brief Add a exploitable data code
             *
             * By default (that is, if no exploitable data codes have
             * been entered),  all data values are considered to be
             * exploitable, even if they appear in missing. If one or
             * more data have been set as exploitable, only those will
             * be considered to be exploitable.
             *
             */
             void add_exploitable(int code);

            /** \brief Add a exploitable data code with synonym character
             *
             * Similar to add_exploitable(int), except that a second
             * parameter is passed to specify a synonym for the main
             * code. If a data item matching `alias` is passed to
             * is_exploitable(), it will be replaced by `code` and the
             * function will return `true`. For a large range of
             * continguous codes, do not use this method within a loop
             * and use the appropriate and more efficient method
             * add_exploitable_range(int, int).
             *
             */
             void add_exploitable_with_alias(int code, int alias);

            /** \brief Add a range of exploitable data codes
             *
             * Like add_exploitable(int), but add all values included
             * between first and last (both included). It is not
             * required at all that all exploitable codes have a
             * synonym.
             *
             */
             void add_exploitable_range(int first, int last);

            /** \brief Add a missing data code
             *
             */
             void add_missing(int code);

            /** \brief Add a missing data code with synonym character
             *
             * Similar to add_missing(int), except that a second
             * parameter is passed to specify a synonym for the main
             * code. If a data item matching `alias` is passed to
             * is_missing(int), it will be replaced by `code` and the
             * function will return `true`.
             *
             */
             void add_missing_with_alias(int code, int alias);

            /** \brief Check if a data value is exploitable
             *
             * If the list of exploitable data is empty, then this
             * method always returns `true`. Otherwise, it checks if the
             * value matches one of the value in the exploitable list.
             * If the value does not match any of the exploitable list,
             * but matches one of the synonyms entered using
             * add_exploitable(int, int), then the method returns `true`
             * and modifies the passed value (so that the synonym is
             * modified to the reference exploitable code). Otherwise,
             * it never modified the value
             *
             */
             bool is_exploitable(int& code) const;

            /** \brief Check if a data value is missing
             *
             * Returns `true` if the code matches one of the missing
             * data codes. Otherwise, returns `true` and modifies the
             * value if the code matches one of the synonyms. If also
             * not, return `false`.
             *
             */
             bool is_missing(int& code) const;

            /** \brief Check data value.
             *
             * Return: (i) the allelic value if it is exploitable as is;
             * (ii) the main value if the provided value is a synonym
             * for an exploitable value; (iii) MISSINGDATA if it is one
             * of the missing data codes or synonyms; (iv) MISSINGDATA
             * if the value is invalid (and in addition the flag argument
             * is set to true).
             *
             */
             int check(int code, bool& flag) const;


        private:

           /** \brief Initializes members (bare instance) */
            void init();

           /** \brief Copy data (over not-so-necessarily empty instance) */
            void copy(const Filter& src);

           /** \brief Free memory (does not re-initialize) */
            void free();

           /** \brief Exploitable data array */
            int * _exploitable;

           /** \brief Number of exploitable data */
            unsigned int _ne;

           /** \brief Number of exploitable data (cache) */
            unsigned int _ne_c;

           /** \brief Exploitable data array ranges */
            int * _exploitable_rg;

           /** \brief Number of exploitable data ranges */
            unsigned int _ne_rg;

           /** \brief Number of exploitable data ranges (cache) */
            unsigned int _ne_c_rg;

           /** \brief Synonym mapping (actual size: 2*_nsyn_c) */
            int * _syn;

           /** \brief Number of synonyms */
            unsigned int _nsyn;

           /** \brief Number of synonyms (cache) */
            unsigned int _nsyn_c;

           /** \brief Missing synonyms mapping (actual size: 2*_n_missing_syn_c) */
            int * _missing_syn;

           /** \brief Number of missing synonyms */
            unsigned int _n_missing_syn;

           /** \brief Number of missing synonyms (cache) */
            unsigned int _n_missing_syn_c;

           /** \brief Missing data array */
            int* _missing;

           /** \brief Number of missing data */
            unsigned int _nm;

           /** \brief Number of missing data (cache) */
            unsigned int _nm_c;
    };
}

#endif
