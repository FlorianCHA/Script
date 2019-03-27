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

#ifndef EGGLIB_DATAHOLDER_HPP
#define EGGLIB_DATAHOLDER_HPP

#include "egglib.hpp"

namespace egglib {

    class Coalesce;
    class Random;

   /** \brief Minimal reimplementation of a vector<int>
    *
    * \ingroup core
    * 
    * Header: <egglib-cpp/DataHolder.hpp>
    *
    */
    class VectorInt {

        private:

            int * _values;
            unsigned int _num;
            unsigned int _res;

            void copy(const VectorInt& src);

        public:

            /// Constructor (default: 0 values)
            VectorInt();

            /// Destructor
            virtual ~VectorInt();

            /// Copy constructor
            VectorInt(const VectorInt& src);

            /// Copy assignment operator
            VectorInt& operator=(const VectorInt& src);

            /// Set the number of vqlues (values are not initialized)
            void set_num_values(unsigned int n);

            /// Get the number of values
            unsigned int get_num_values() const;

            /// Set a value
            void set_item(unsigned int i, int value);

            /// Get a value
            int get_item(unsigned int i) const;

            /// Release memory
            void clear();
    };

   /** \brief Integer data set
    *
    * \ingroup core
    *
    * Holds a data set with associated sample names and group
    * information. The data consists of given numbers of ingroup and
    * outgroup samples, which can all have a different number of sites,
    * unless the object is configured to be a matrix. In that cases, it
    * is assumed that all samples have the same number of sites as the
    * first loaded sample. There can be any number of group levels (but
    * this number must be the same for all samples), meaning that
    * samples can be described by several group labels in addition to
    * their name. Group labels are not group indices (they do not need
    * to be consecutive). There is a separate data set for sample
    * belonging to the outgroup. There can be any number of outgroup
    * samples. The outgroup has always one level of group labels, but
    * the labels are not initialized. All data are represented by signed
    * integers. Note that none of the accessors performs out-of-bound
    * checking. The user is responsible to provide valid indices. This
    * class follows a memory caching system: allocated memory is never
    * freed with the aim of efficiently reusing the same object.
    *
    * Header: <egglib-cpp/DataHolder.hpp>
    *
    */
    class DataHolder {

        protected:

            void _init();
            void _alloc_nsi(unsigned int nsi);
            void _alloc_nso(unsigned int nso);
            void _alloc_ng(unsigned int ng);
            void _alloc_ls(unsigned int ls);
            void _alloc_lsi(unsigned int i, unsigned int ls);
            void _alloc_lso(unsigned int i, unsigned int ls);
            void _alloc_lni(unsigned int i, unsigned int ln);
            void _alloc_lno(unsigned int i, unsigned int ln);
            void _free();
            void _copy(const DataHolder& src);
            void _nsam_i_helper(unsigned int nsam);
            void _nsam_o_helper(unsigned int nsam);

            bool _is_matrix;
            unsigned int _nsi;
            unsigned int _nsi_r;
            unsigned int _nso;
            unsigned int _nso_r;
            unsigned int _ls;
            unsigned int * _lsi;            // size: _nsi_r
            unsigned int * _lsi_r;          // size: _nsi_r
            unsigned int * _lso;            // size: _nso_r
            unsigned int * _lso_r;          // size: _nso_r
            unsigned int * _lni;            // size: _nsi_r
            unsigned int * _lni_r;          // size: _nsi_r
            unsigned int * _lno;            // size: _nso_r
            unsigned int * _lno_r;          // size: _nso_r
            unsigned int _ng;
            unsigned int * _ng_r;           // size: _nsi_r

            int ** _datai;                  // size: _nsi_r x _lsi_r[i]
            int ** _datao;                  // size: _nso_r x _lso_r[i]
            char ** _namesi;                // size: _nsi_r x _lni_r[i]
            char ** _nameso;                // size: _nso_r x _lso_r[i]
            unsigned int ** _groups_i;      // size: _nsi_r x _ng_r[i]
            unsigned int * _groups_o;       // size: _nso_r

        public:

           /** \brief Default constructor
            *
            * Create an empty matrix. The object will be usable only
            * when resizers will have been called.
            *
            * \param is_matrix determines if the object is configured to
            *   a matrix (that is, where all ingroup and outgroup
            *   samples have the same number of sites). If so, the user
            *   is responsible to ensure that all loaded samples are
            *   consistent. This value can be changed using the method
            *   is_matrix(bool).
            *
            */
            DataHolder(bool is_matrix = false);

           /** \brief Copy constructor
            *
            * The reserved memory of the source is not copied.
            *
            */
            DataHolder(const DataHolder& src);

           /** \brief Assignment operator
            *
            * The reserved memory of the source is not copied. The
            * reserved memory of the current object is retained.
            *
            */
            DataHolder& operator=(const DataHolder& src);

           /** \brief Destructor
            *
            */
            virtual ~DataHolder();

            /// True if all names are got for PhyMl
            bool valid_phyml_names() const;

            /// True if all data are nucleotides expected by PhyML (for alignment)
            bool valid_phyml_nt() const;

            /// True if all data are amino acids expected by PhyML (for alignment)
            bool valid_phyml_aa() const;

           /** \brief Configure the object (not) to be a matrix
            *
            * If a non-matrix object is converted to a matrix, the user
            * is responsible of ensuring that all samples (including
            * outgroup) have the same number of samples. The method will
            * assume that all samples have the same number of sites as
            * the first sample among the ingroup and outgroup. There is
            * no requirement for converting a matrix to a non-matrix.
            *
            */
            void set_is_matrix(bool flag);

           /** \brief Check if the object is configured to be a matrix
            *
            */
            bool get_is_matrix() const;

           /** \brief Reserve memory to speed up data loading
            *
            * This method does not change the size of the data set
            * contained in the instance, but reserves memory in order to
            * speed up incremental loading of data. The passed values
            * are not required to be accurate. In case the instance has
            * allocated more memory than what is requested, nothing is
            * done (this applies to both dimensions independently). It
            * is always valid to use 0 for any values (in that case,
            * nothing is done). Note that one character is always
            * pre-allocated for all names.
            *
            * \param nsi expected number of ingroup samples.
            * \param nso expected number of outgroup samples.
            * \param ln expected length of names.
            * \param ng expected number of groups.
            * \param ls expected number of sites (the same for all
            *   ingroup and outgroup samples, whichever the object is a
            *   matrix or not).
            *
            */
            void reserve(unsigned int nsi, unsigned int nso, unsigned int ln,
                         unsigned int ng, unsigned int ls);

           /** \brief Get the number of ingroup samples
            *
            */
            unsigned int get_nsam_i() const;

           /** \brief Set the number of ingroup samples
            *
            * Perform memory allocation as needed but does not
            * initialize new values (except names). If the object is a
            * matrix, new samples are set to have the current number of
            * sites. Otherwise, new samples have no sites. Set the
            * number of samples to a smaller value equals to remove the
            * last samples.
            *
            */
            void set_nsam_i(unsigned int nsam);

           /** \brief Get the number of outgroup samples
            *
            */
            unsigned int get_nsam_o() const;

           /** \brief Set the number of outgroup samples
            *
            * As nsam(unsigned int) but for the outgroup data table.
            *
            */
            void set_nsam_o(unsigned int nsam);

           /** \brief Get the number of sites
            *
            * This method **may not** be called on a non-matrix object.
            *
            */
            unsigned int get_nsit() const;

           /** \brief Get the number of sites for an ingroup sample
            *
            * This method **may not** be called on a matrix object.
            *
            */
            unsigned int get_nsit_i(unsigned int sam) const;

           /** \brief Get the number of sites for an outgroup sample
            *
            * This method **may not** be called on a matrix object.
            *
            */
            unsigned int get_nsit_o(unsigned int sam) const;

           /** \brief Set the number of sites
            *
            * Perform memory allocation as needed but does not
            * initialize new values. It is possible to use this method
            * for both matrix and non-matrix objects. In both cases,
            * the effective result is that all ingroup and outgroup
            * samples are resized to the specified number of sites.
            *
            */
            void set_nsit(unsigned int val);

           /** \brief Set the number of sites for an ingroup sample
            *
            * Similar to nsit(unsigned int) but for only one sample
            * of the ingroup. Available only for non-matrix objects.
            *
            */
            void set_nsit_i(unsigned int sam, unsigned int val);

           /** \brief Set the number of sites for an outgroup sample
            *
            * Similar to nsit(unsigned int) but for only one sample of
            * the outgroup. Available only for non-matrix objects.
            *
            */
            void set_nsit_o(unsigned int sam, unsigned int val);

           /** \brief Insert sites at a given position
            *
            * Increase the number of sites for all samples. This method
            * may be used for matrix or non-matrix objects. Note that
            * the inserted sites are not initialized.
            *
            * \param pos the position at which to insert sites. The new
            *   sites are inserted before the specified index. Use 0 to
            *   add sites at the beginning of the sequence, and the
            *   current number of sites to add sites at the end. If the
            *   value is larger than the current number of sites, sites
            *   are added at the end of the sequence. Therefore it is
            *   possible to use egglib::MAX as the position to specify
            *   that new sites must be inserted at the end.
            * \param num number of sites at add.
            *
            */
            void insert_sites(unsigned int pos, unsigned int num);

           /** \brief Insert sites at a given position for an ingroup sample
            *
            * As insert_sites(unsigned int, unsigned int, int) but for
            * only one sample of the ingroup. Available only for
            * non-matrix objects.
            *
            */
            void insert_sites_i(unsigned int sam, unsigned int pos, unsigned int num);

           /** \brief Insert sites at a given position for an outgroup sample
            *
            * As insert_sites(unsigned int, unsigned int, int) but for
            * only one sample of the outgroup. Available only for
            * non-matrix objects.
            *
            */
            void insert_sites_o(unsigned int sam, unsigned int pos, unsigned int num);

           /** \brief Get the number of group levels
            *
            */
            unsigned int get_ngroups() const;

           /** \brief Set the number of group levels
            *
            * Perform memory allocation as needed but does not
            * initialize new values.
            *
            */
            void set_ngroups(unsigned int ngrp);

           /** \brief Get an ingroup data entry
            *
            * \param sam sample index.
            * \param sit site index.
            *
            * The indices must be valid, otherwise a segmentation fault
            * or aberrant behaviour will occur.
            *
            */
            int get_i(unsigned int sam, unsigned int sit) const;

           /** \brief Set an ingroup data entry
            *
            * \param sam sample index.
            * \param sit site index.
            * \param value allele value.
            *
            * The indices must be valid, otherwise a segmentation fault
            * or aberrant behaviour will occur.
            *
            */
            void set_i(unsigned int sam, unsigned int sit, int value);

           /** \brief Get an outgroup data entry
            *
            * \param sam sample index.
            * \param sit site index.
            *
            * The indices must be valid, otherwise a segmentation fault
            * or aberrant behaviour will occur.
            *
            */
            int get_o(unsigned int sam, unsigned int sit) const;

           /** \brief Set an outgroup data entry
            *
            * \param sam sample index.
            * \param sit site index.
            * \param value allele value.
            *
            * The indices must be valid, otherwise a segmentation fault
            * or aberrant behaviour will occur.
            *
            */
            void set_o(unsigned int sam, unsigned int sit, int value);

           /** \brief Get a group label
            *
            * \param sam sample index.
            * \param lvl group level index.
            *
            * The indices must be valid, otherwise a segmentation fault
            * or aberrant behaviour will occur.
            *
            */
            unsigned int get_group_i(unsigned int sam, unsigned int lvl) const;

           /** \brief Get the group label for an outgroup sample
            *
            * \param sam sample index.
            *
            * The indexe must be valid, otherwise a segmentation fault
            * or aberrant behaviour will occur. There is necessarily one
            * group level for outgroups, and the default value is 0.
            *
            */
            unsigned int get_group_o(unsigned int sam) const;

           /** \brief Set a group label
            *
            * \param sam sample index.
            * \param lvl group level index.
            * \param label group label.
            *
            * The indices must be valid, otherwise a segmentation fault
            * or aberrant behaviour will occur.
            *
            */
            void set_group_i(unsigned int sam, unsigned int lvl, unsigned int label);

           /** \brief Set the group label for an outgroup sample
            *
            * \param sam sample index.
            * \param label group label.
            *
            * The indexe must be valid, otherwise a segmentation fault
            * or aberrant behaviour will occur. There is necessarily one
            * group level for outgroups.
            *
            */
            void set_group_o(unsigned int sam, unsigned int label);

           /** \brief Get an ingroup name
            * 
            */
            const char * get_name_i(unsigned int sam) const;

           /** \brief Set an ingroup name
            * 
            */
            void set_name_i(unsigned int sam, const char * name);

           /** \brief Add a character to the specified ingroup name
            * 
            */
            void name_appendch_i(unsigned int sam, char ch);

           /** \brief Add character at the end of specified ingroup name
            * 
            */
            void name_append_i(unsigned int sam, const char * ch);

           /** \brief Get an outgroup name
            * 
            */
            const char * get_name_o(unsigned int sam) const;

           /** \brief Set an outgroup name
            * 
            */
            void set_name_o(unsigned int sam, const char * name);

           /** \brief Add a character to the specified outgroup name
            * 
            */
            void name_appendch_o(unsigned int sam, char ch);

           /** \brief Add characters at the end of specified ingroup name
            * 
            */
            void name_append_o(unsigned int sam, const char * ch);

           /** \brief Move a sample to the outgroup
            *
            * The specified sample is moved to the outgroup and its
            * group labels are discarded. Obviously, this decreases the
            * ingroup size by 1, and increases the outgroup size
            * accordingly.
            *
            * \param sam ingroup sample index.
            * \param label group label to assign to the sample after it
            *   it is moved to the outgroup (use 0 if not relevant).
            *
            */
            void to_outgroup(unsigned int sam, unsigned int label);

           /** \brief Delete a sample
            *
            * Delete the specified sample and decrease index of all
            * subsequent samples by one. If there is only one sequence
            * in the instance, set the number of sites (the maximal
            * number of sites for non-matrix objects) to 0.
            *
            */
            void del_sample_i(unsigned int sam);

           /** \brief Delete an outgroup sample
            *
            * Delete the specified sample and decrease index of all
            * subsequent samples by one. If there is only one sequence
            * in the instance, set the number of sites (the maximal
            * number of sites for non-matrix objects) to 0.
            *
            */
            void del_sample_o(unsigned int sam);

           /** \brief Delete a given range of sites
            *
            * The sites are removed for all ingroup and outgroup
            * samples. The index must be valid. This method may be used
            * for both matrix and  non-matrix objects.
            *
            * \param start start position of the range to remove.
            * \param stop stop position of the range to remove (this
            *   site IS NOT removed).
            *
            * If the stop argument is larger than the number of sites
            * (the number of sites for this sample in the case of a
            * non-matrix object), then sites are removed until the end
            * of the sequence. If the start argument is larger or equal
            * to the number of sites (the number of sites for this
            * sample in the case of a non-matrix object), then nothing
            * is done.
            *
            */
            void del_sites(unsigned int start, unsigned int stop);

           /** \brief Delete a given range of sites for an ingroup sample
            *
            * As del_sites(unsigned int, unsigned int) but for a single
            * sample. This method **may not** be called on a matrix
            * object.
            *
            */
            void del_sites_i(unsigned int sam, unsigned int start, unsigned int stop);

           /** \brief Delete a given range of sites for an outgroup sample
            *
            * As del_sites(unsigned int, unsigned int) but for a single
            * sample. This method **may not** be called on a matrix
            * object.
            *
            */
            void del_sites_o(unsigned int sam, unsigned int start, unsigned int stop);

           /** \brief Restore object to initial state
            *
            * This method is designed to allow reusing the object and
            * reusing previously allocated memory. All data contained in
            * the instance is considered to be lost, but allocated
            * memory is actually retained to speed up later resizing
            * operations.
            *
            */
            void reset(bool is_matrix);

           /** \brief Clear object
            *
            * Actually clears all memory stored by the object (including
            * cache). All memory vector data are effectively lost and
            * memory is released.
            *
            */
            void clear(bool is_matrix);

           /** \brief Find the start position of the first match of a motif
            *
            * \param sam sample index.
            * \param of_outgroup specifies whether the sample is in the
            *   outgroup.
            * \param motif the list of integers representing the motif
            *   to be found.
            * \param start at which to start search. No returned value
            *   will be smaller than this value.
            * \param stop position at which to stop search (the motif
            *   cannot overlap this position). No returned value will be
            *   larger than *stop* - *n*.
            * \return The index of the start positon of the first exact
            *   match for the passed set of values, or egglib::MAX if
            *   no match was found (within the specified region).
            *
            */
            unsigned int find(unsigned int sam, bool of_outgroup, VectorInt& motif, unsigned int start=0, unsigned int stop=MAX) const;

           /** \brief Test if all sequences have the same length
            *
            * True if no sequences at all. Only valid for containers.
            *
            */
            bool is_equal() const;
    };

   /** \brief Insert non-varying sites within alignments
    *
    * \ingroup core
    *
    * This class allows to add non-varying sites within an alignment at
    * given positions. The procedure below must be strictly followed:
    *
    * - Create an instance. The constructors takes no arguments.
    *
    * - Load a DataHolder instance using load(). It is required that it
    *   is an alignment (a matrix). The instance will create an array of
    *   positions internally.
    *
    * - Specify the desired length of the final alignment using
    *   set_length().
    *
    * - Specify the position of all sites of the original alignment.
    *   This can be achieved by three ways:
    *    + Specify manually all positions as real numbers using
    *      set_position().
    *    + Specify manually all positions as extant indexes (as integer
    *      values) using set_round_position() for all positions. If this
    *      approach is used, it is necessary to set the round option of
    *      intersperse to false.
    *    + Pass the reference to the Coalesce instance that has
    *      generated the alignment (assuming it is a simulation) and let
    *      the instance find itself the site positions, with the method
    *      get_positions().
    *
    * - Specify the list of alleles values used for non-varying
    *   positions using set_num_alleles() and then set_allele() as many
    *   times as needed. If there is more than one allele, non-varying
    *   alleles will be picked randomly. This method can be skipped (by
    *   default, the value corresponding to `A` will be used).
    *
    * - Provide the address of a random number generator using
    *   set_random() (it is always needed).
    *
    * - Call intersperse(). This will change the DataHolder originally
    *   loaded.
    *
    * Header: <egglib-cpp/DataHolder.hpp>
    *
    */
    class IntersperseAlign {

        public:

           /** \brief Constructor.
            *
            */
            IntersperseAlign();

           /** \brief Destructor.
            *
            */
            ~IntersperseAlign();

           /** \brief Loads a data set.
            *
            * The loaded data set may contain any number of sites (even
            * zero). 
            *
            * This method does not reset the random number generator,
            * the final alignment length, the position of sites (unless
            * the new alignment has a different number of sites compared
            * with the previous one) or the number and value of
            * non-varying alleles.
            *
            */
            void load(DataHolder& data);

           /** \brief Specifies the desired length of the final alignment.
            *
            * If this method is skipped, the default value is 0 (the
            * final alignment is identical to the original one), or the
            * previously specified value (if set_length() has been
            * called previously).
            *
            */
            void set_length(unsigned int length);

           /** \brief Sets the position of one of the sites of the original alignment.
            *
            * A DataHolder reference must have been loaded using load().
            * This method allows to specify the position of each of the
            * sites of the passed DataHolder instance. Note that the
            * position of all sites must be specified, that positions
            * must always be increasing (consecutive positions might be
            * equal), and all positions must be at least 0 and at most
            * 1.
            *
            */
            void set_position(unsigned int index, double position);

           /** \brief Sets the position of one of the sites of the original alignment.
            *
            * A DataHolder reference must have been loaded using load().
            * This method allows to specify the position of each of the
            * sites of the passed DataHolder instance. Note that the
            * position of all sites must be specified, that positions
            * must always be increasing (consecutive positions might be
            * equal), and all positions must be at least 0 and at most
            * ls-1 where ls is the length of the final alignment.
            *
            * If you use this method, you must use it for all sites and
            * then set the argument of intersperse() to false.
            *
            */
            void set_round_position(unsigned int index, unsigned int position);

           /** \brief Gets automatically the positions of sites of the original alignment.
            *
            * A DataHolder reference must have been loaded using load(),
            * and this DataHolder object must be the last one simulated
            * using the Coalesce object whose reference is passed. This
            * method will load the positions of all sites as provided by
            * Coalesce.
            * 
            */
            void get_positions(const Coalesce& coalesce);

           /** \brief Specifies the number of possible alleles at inserted positions.
            *
            * This method may be called at any time. The value must be
            * at least one. If more than one, alleles at inserted
            * positions will be picked randomly (they will be fixed
            * among samples). All alleles must be specified using
            * set_allele(). The default value is one and the default
            * first allele is `A`. It is possible not to specify the
            * first allele even if the number of alleles is increased
            * (the `A` value will be retained).
            *
            */
            void set_num_alleles(unsigned int num);

           /** \brief Sets an allele for inserted positions.
            *
            * The number of alleles must have been fixed using
            * set_num_alleles(). The default value for the first allele
            * is `A`.
            *
            */
            void set_allele(unsigned int index, int allele);

           /** \brief Provides a random number generator.
            *
            * It is always required to provide a random number
            * generator, even if set_num_alleles() is one.
            *
            */
            void set_random(Random * random);

           /** \brief Insert non-varying sites.
            *
            * This method modifies the DataHolder instance that has been
            * loaded. It is required to have loaded one, and to have
            * specified the position of each of its sites. It is also
            * logical (but not formally required) to have specified the
            * desired length of the alignment. It is possible to specify
            * more than one alleles for inserted positions. It is
            * required to have passed a random number generator.
            *
            * After call to this method, the loaded DataHolder instance
            * will have a length equal to the value specified using
            * set_length(), unless the original DataHolder was longer
            * (in such case, it is not changed at all).
            *
            * \param round_positions a boolean indicating if site
            *   positions must be rounded. Set it to false if already
            *   rounded positions have been provided.
            *
            */
            void intersperse(bool round_positions = true);

        protected:

            DataHolder * _data;
            unsigned int _nsites;
            unsigned int _length;
            unsigned int _res_positions;
            double * _positions;
            unsigned int * _round_positions;
            unsigned int * _offset;
            unsigned int _final_offset;
            unsigned int _num_alleles;
            unsigned int _res_alleles;
            int * _alleles;
            Random * _random;

        private:

            IntersperseAlign(const IntersperseAlign& src) {}
            IntersperseAlign& operator=(const IntersperseAlign& src) { return * this; }
    };
}

#endif
