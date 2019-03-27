/*
    Copyright 2015-2017 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_STRUCTURE_HPP
#define EGGLIB_STRUCTURE_HPP

#include "egglib.hpp"

namespace egglib {

    class DataHolder;
    class StructureCluster;
    class StructurePopulation;
    class StructureIndiv;

   /** \brief Manage hierarchical group structure.
    *
    * \ingroup core
    *
    */
    class StructureHolder {

        private:

            unsigned int _ns;
            unsigned int _no;
            unsigned int _required_ns;
            unsigned int _required_no;
            unsigned int _ploidy;
            unsigned int _num_clust;
            unsigned int _num_pop;
            unsigned int _num_indiv_i;
            unsigned int _num_indiv_o;
            StructureCluster ** _clusters;
            StructurePopulation ** _pops;
            StructureIndiv ** _indivs_i;
            StructureIndiv ** _indivs_o;

            unsigned int _num_clust_c;
            unsigned int _num_pop_c;
            unsigned int _num_indiv_i_c;
            unsigned int _num_indiv_o_c;

            unsigned int _num_filter;
            unsigned int _num_filter_c;
            unsigned int * _filter;

            unsigned int _last_indiv;

            void init();
            void free();

            StructureHolder(const StructureHolder& src) {}
            StructureHolder& operator=(const StructureHolder& src) { return *this; }

        public:

            StructureHolder();
            ~StructureHolder();

           /** \brief Reset to defaults.
            *
            */
            void reset();

           /** \brief Reset pop filter only.
            *
            */
            void reset_filter();

           /** \brief Pre-alloc filter table.
            *
            */
            void reserve_filter(unsigned int howmany);

           /** \brief Add a population label to filter.
            *
            * If at least one is passed, process only those passed. By
            * default, include all populations. Use reset_filter() to
            * reset to default.
            *
            */
            void add_pop_filter(unsigned int lbl);

           /** \brief Ensure ploidy is consistent and optionally equal to passed value.
            *
            * Automatically called by get_structure(). Need to be called
            * if process_ingroup() and/or process_outgroup() is used.
            *
            * Value must be >0.
            *
            */
            void check_ploidy(unsigned int value=UNKNOWN);

           /** \brief Get ploidy.
            *
            * Default is UNKNOWN.
            *
            */
            unsigned int get_ploidy() const;

           /** \brief Process labels from a DataHolder.
            *
            * Use UNKNOWN for any level to skip (but skipping
            * individuals is not the same as skipping clusters/pops).
            * You must set pop filter separately.
            *
            */
            void get_structure(DataHolder& data, unsigned int lvl_clust,
                    unsigned int lvl_pop, unsigned int lvl_indiv,
                    unsigned int ploidy, bool skip_outgroup);

           /** \brief Add a cluster with no samples in it
            *
            */
            StructureCluster * add_cluster(unsigned int label);

           /** \brief Add a population with no samples in it
            *
            */
            StructurePopulation * add_population(unsigned int label,
                StructureCluster * cluster);

           /** \brief Add an ingroup individual with no samples in it
            *
            */
            StructureIndiv * add_individual_ingroup(unsigned int label,
                StructureCluster * cluster, StructurePopulation * population);

           /** \brief Add an outgroup individual with no samples in it
            *
            */
            StructureIndiv * add_individual_outgroup(unsigned int label);

           /** \brief Add one ingroup sample
            *
            */
            void add_sample_ingroup(unsigned int sam_idx,
                                StructureCluster * cluster,
                                StructurePopulation * population,
                                StructureIndiv * indiv);

           /** \brief Add one outgroup sample
            *
            */
            void add_sample_outgroup(unsigned int sam_idx,
                                StructureIndiv * indiv);

           /** \brief Process one ingroup sample.
            *
            */
            void process_ingroup(unsigned int idx, unsigned int lbl_clust,
                unsigned int lbl_pop, unsigned int lbl_indiv);

           /** \brief Process one outgroup sample.
            *
            */
            void process_outgroup(unsigned int idx, unsigned int lbl_indiv);

           /** \brief Number of clusters. */
            unsigned int num_clust() const;

           /** \brief Number of populations (total). */
            unsigned int num_pop() const;

           /** \brief Number of ingroup individuals (total). */
            unsigned int num_indiv_ingroup() const;

           /** \brief Number of outgroup individuals. */
            unsigned int num_indiv_outgroup() const;

           /** \brief Get a cluster. */
            const StructureCluster& get_cluster(unsigned int idx) const;

           /** \brief Get a population. */
            const StructurePopulation& get_population(unsigned int idx) const;

           /** \brief Get an ingroup individual. */
            const StructureIndiv& get_indiv_ingroup(unsigned int idx) const;

           /** \brief Get an outgroup individual. */
            const StructureIndiv& get_indiv_outgroup(unsigned int idx) const;

           /** \brief Get number of ingroup samples. */
            unsigned int get_ns() const;

           /** \brief Get number of outgroup samples. */
            unsigned int get_no() const;

           /** \brief Get required number of ingroup samples. */
            unsigned int get_ns_req() const;

           /** \brief Get required number of outgroup samples. */
            unsigned int get_no_req() const;

           /** \brief Index of the population containing this sample (default: MISSING). */
            unsigned int get_pop_index(unsigned int) const;

            /// \brief Copy data frome source object
            void copy(const StructureHolder& source);
    };

   /** \brief Manage a cluster.
    *
    */
    class StructureCluster {

        private:

            StructureHolder * _parent;
            unsigned int _label;
            unsigned int _ns;
            unsigned int _num_pop;
            unsigned int _num_pop_c;
            unsigned int _num_indiv;
            unsigned int _num_indiv_c;
            StructurePopulation ** _pops;
            StructureIndiv ** _indivs;

            void init();
            void free();

            StructureCluster() {}
            StructureCluster(const StructureCluster& src) {}
            StructureCluster& operator=(const StructureCluster& src) { return *this; }

        public:

            StructureCluster(StructureHolder * holder, unsigned int label);
            ~StructureCluster();

           /** \brief Restore defaults. */
            void reset(StructureHolder * holder, unsigned int label);

           /** \brief Add and create a population. */
            StructurePopulation * add_pop(unsigned int label);

           /** \brief Add and create an individual. */
            StructureIndiv * add_indiv(StructurePopulation * pop, unsigned int label);

           /** \brief Add a sample. */
            void add_sample();

           /** \brief Number of populations. */
            unsigned int num_pop() const;

           /** \brief Number of individuals (total for this cluster). */
            unsigned int num_indiv() const;

           /** \brief Get a population. */
            const StructurePopulation& get_population(unsigned int idx) const;

           /** \brief Get an individual. */
            const StructureIndiv& get_indiv(unsigned int idx) const;

           /** \brief Get label. */
            unsigned int get_label() const;
    };

   /** \brief Manage a population.
    *
    */
    class StructurePopulation {

        private:

            StructureHolder * _holder;
            StructureCluster * _cluster;
            unsigned int _label;
            unsigned int _ns;
            unsigned int _num_indiv;
            unsigned int _num_indiv_c;
            StructureIndiv ** _indivs;

            void init();
            void free();

            StructurePopulation() {}
            StructurePopulation(const StructurePopulation& src) {}
            StructurePopulation& operator=(const StructurePopulation& src) { return *this; }

        public:

            StructurePopulation(StructureHolder * holder,
                                StructureCluster * cluster,
                                unsigned int label);
            ~StructurePopulation();

           /** \brief Restore defaults. */
            void reset(StructureHolder * holder,
                       StructureCluster * cluster,
                       unsigned int label);

           /** \brief Add and create an individual. */
            StructureIndiv * add_indiv(unsigned int label);

           /** \brief Add a sample. */
            void add_sample();

           /** \brief Number of individuals. */
            unsigned int num_indiv() const;

           /** \brief Get an individual. */
            const StructureIndiv& get_indiv(unsigned int idx) const;

           /** \brief Get label. */
            unsigned int get_label() const;

           /** \brief Get containing cluster. */
            StructureCluster * get_cluster();
    };

   /** \brief Manage an individual.
    *
    */
    class StructureIndiv {

        private:

            StructureHolder * _holder;
            StructureCluster * _cluster;
            StructurePopulation * _population;
            unsigned int _label;
            unsigned int _num_sam;
            unsigned int _num_sam_c;
            unsigned int * _samples;

            void init();
            void free();

            StructureIndiv() {}
            StructureIndiv(const StructureIndiv& src) {}
            StructureIndiv& operator=(const StructureIndiv& src) { return *this; }

        public:

            StructureIndiv(StructureHolder * holder,
                           StructureCluster * cluster,
                           StructurePopulation * population,
                           unsigned int label);
            ~StructureIndiv();

           /** \brief Restore defaults. */
            void reset(StructureHolder * holder,
                       StructureCluster * cluster,
                       StructurePopulation * population,
                       unsigned int label);

           /** \brief Number of samples. */
            unsigned int num_samples() const;

           /** \brief Add a sample. */
            void add_sample(unsigned int index);

           /** \brief Get a sample. */
            unsigned int get_sample(unsigned int idx) const;

           /** \brief Get label. */
            unsigned int get_label() const;

           /** \brief Get containing cluster (NULL if outgroup). */
            StructureCluster * get_cluster();

           /** \brief Get containing population (NULL if outgroup). */
            StructurePopulation * get_population();
    };
}

#endif
