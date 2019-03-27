/*
    Copyright 2012-2014 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_TREE_HPP
#define EGGLIB_TREE_HPP

#include "DataHolder.hpp"

namespace egglib {

    class Params;
    class Random;

   /** \brief Implements a node of bifurcating tree
    *
    * \ingroup coalesce
    *
    * Companion of the class Tree. A given Node instance is connected to
    * either no or two other instances ("sons"). Sons addresses are lost
    * by calling the destructor or the assignment operator. The sons are
    * never destroyed.
    *
    * The member default values are:
    *  \li t: 0
    *  \li L: 0
    *  \li nmut: 0
    *  \li label: egglib::UNKNOWN
    *  \li son1: egglib::UNKNOWN
    *  \li son2: egglib::UNKNOWN
    *
    * Header: <egglib-cpp/Tree.hpp>
    *
    */
    class Node {

        public:

           /** \brief Default constructor
            *
            * Initialize members to default values.
            *
            */
            Node();

           /** \brief Copy constructor
            *
            * Copy members values and sons addresses.
            *
            */
            Node(const Node& src);

           /** \brief Assignment operator
            *
            * Copy members values and sons addresses.
            *
            */
            Node& operator=(const Node& src);

           /** \brief Destructor */
            virtual ~Node();

           /** \brief Specify an internal node
            *
            * Other members than those specified by argument are not
            * changed (except boolean is_terminal).
            *
            * \param t node time.
            * \param son1 first descending node.
            * \param son2 second descending node.
            *
            */
            void set_internal(double t, unsigned int son1, unsigned int son2);

           /** \brief Specify a terminal node
            *
            * Other members than those specified by argument are not
            * changed (except boolean is_terminal).
            *
            * \param t node time (>=0).
            * \param label sample label.
            *
            */
            void set_terminal(double t, unsigned int label);

            /** \brief Internal or terminal node */
            bool is_terminal() const;

            /** \brief Gets label */
            unsigned int label() const;

           /** \brief Get first son
            *
            * If the node is terminal, UNKNOWN is returned.
            *
            */
            unsigned int son1() const;

           /** \brief Get second son
            *
            * If the node is terminal, UNKNOWN is returned.
            *
            */
            unsigned int son2() const;

           /** \brief Get node date */
            double date() const;

           /** \brief Get branch length */
            double get_L() const;

           /** \brief Set branch length
            *
            * The branch is going from the ancestror node to this node.
            * As a result, the root has a branch length of 0.
            *
            */
            void set_L(double value);

           /** \brief Number of mutations borne by this node */
            unsigned int nmut() const;

           /** \brief Get the site index of a mutation
            *
            * Warning: no out-of-bound check is performed!
            *
            */
            unsigned int mutationSite(unsigned int mut) const;

           /** \brief Get the site position of a mutation
            *
            * Warning: no out-of-bound check is performed!
            *
            * TODO check if this method is needed anyway
            *
            */
            double mutationPos(unsigned int mut) const;

           /** \brief Add a mutation
            *
            * \param site site index
            * \param pos site position (between 0 and 1)
            *
            */
            void addMutation(unsigned int site, double pos);

           /** \brief Restore initial state
            *
            * Deliver object as new, but retain reserved memory.
            *
            */
            void reset();

           /** \brief Remove already stored mutations
            *
            */
            void clear_mutations();

        protected:

            void copy(const Node& src); // allocate to src._nmut
            void init();  // like reset, but initalize reserved memory
            void alloc(unsigned int nmut); // allocate to nmut and set nmut to _nmut

            bool _is_terminal;
            double _date;
            double _length;
            unsigned int _son1;
            unsigned int _son2;
            unsigned int _nmut;
            unsigned int _reserved;
            unsigned int * _mutationSite;
            double * _mutationPos;
            unsigned int _label;
    };

   /** \brief This class handles a bifurcating tree
    *
    * \ingroup coalesce
    *
    * Each tree corresponds to a given segment of chromosome. The class
    * provides a method for generating allelic data for a given site,
    * and exposes to root node, allowing recursive tree exploration.
    *
    * Header: <egglib-cpp/Tree.hpp>
    *
    */
    class Tree {

        public:

           /** \brief Constructor
            *
            * The constructor builds a preliminary tree consisting of
            * only leaves, all available for coalescence. Leaves are
            * labelled from 0 to (numberOfLeaves - 1).
            *
            * \param numberOfLeaves number of leaves (unconnected).
            * \param start start of the region covered by the tree.
            * \param stop end of the region covered by the tree
            *
            */
            Tree(unsigned int numberOfLeaves, double start, double stop);

           /** \brief Copy constructor
            *
            * Deep copy of the tree structure.
            */
            Tree(const Tree& tree);                                                                                // copy constructor

           /** \brief Copy assignment operator
            *
            * Deep copy of the tree structure.
            *
            */
            Tree& operator=(Tree tree);                                                                            // assignement

           /** \brief Destructor */
            virtual ~Tree();

           /** \brief Generate a new tree by recombination
            *
            * It is assumed that the recombination point lies between
            * the start and stop positions of this tree. This tree is
            * copied to the new tree whose address is passed. The
            * current tree's interval is shrunk to the what is at the
            * left of the recombination point (from start to point), and
            * the complement to the new tee.
            *
            */
            void recomb(double point, Tree * new_tree);

           /** \brief %Coalesce two lineages
            *
            * \param nodeIndex1 first node to coalesce.
            * \param nodeIndex2 second node to coalesce.
            * \param time absolute time point of the event.
            * \return The index of the parent node.
            *
            * The two indices must be different and be valid indices for
            * this Tree object.
            *
            */
            unsigned int coal(unsigned int nodeIndex1, unsigned int nodeIndex2, double time);                                        // coalesce 2 nodes, return ancestor index

           /** \brief Generate genetic data for a given site
            *
            * \param site index of the site (that is, column of the
            * DataHolder) for which mutations should be generated.
            * \param data DataHolder object in which mutations should be
            * placed.
            * \param params parameters holder (where parameters should be
            * read).
            * \param random a random number generator.
            * \return The number of mutations (>=0) found at this site.
            *
            * The DataHolder instance should have the correct dimensions but
            * its values do not need to be initialized.
            *
            */
            unsigned int mutate(unsigned int site, DataHolder& data, const Params *params, Random *random);                // mutate a site (calls _mutate)

           /** \brief Start position of the segment covered by the tree */
            double start() const;

           /** \brief End position of the segment covered by the tree */
            double stop() const;

           /** \brief Length of the segment covered by the tree
            *
            * The tree length is not taken into account; this is exactly
            * stop() - start().
            *
            */
            double cov() const;

           /** \brief Total tree length */
            double L() const;                                                                                      // total tree length

           /** \brief Get the root of the tree
            *
            * It is mandatory that the tree is completed (all coalescence
            * has been performed, with only one lineage left) before using
            * this method. In this case, recursing over node descendants
            * allows to explore the tree.
            *
            */
            Node * root() const;

           /** \brief Get the number of nodes */
            unsigned int nnodes() const;

           /** \brief Get a node address by its index
            *
            * The index must be within range.
            *
            */
            Node * const node(unsigned int i) const;

           /** \brief Add a node and return its index */
            unsigned int addNode(double t, unsigned int label);

           /** \brief Reset to initial state
            *
            * Restore the object as newly created, but retain allocated
            * memory arrays. Leaves are automatically labelled from 0 up
            * to (numberOfLeaves - 1).
            *
            * \param numberOfLeaves number of leaves (unconnected).
            * \param start start of the region covered by the tree.
            * \param stop end of the region covered by the tree
            *
            */
            void reset(unsigned int numberOfLeaves, double start, double stop);

           /** \brief Remove already stored mutations
            *
            */
            void clear_mutations();

        protected:

            void init();  // members to zero
            void free();  // free memory (that's all)
            void set(unsigned int numberOfLeaves, double start, double stop);  // (re)initialize instance (don't reinitialize nodes as terminal)
            void copy(const Tree& tree);  // import data from tree
            void realloc(unsigned int size);   // resize Node array, create objects and initialize them
            void alloc_from(unsigned int size, Node ** source);  // resize Node array, create objects and copy them from source

            int next_allele(int allele);   // draw the next allele (based on params instance)
            void r_mutate(unsigned int site, Node * node, int allele);  // mutate helper

            double _start;   // segment start
            double _stop;    // segment end
            double _cov;     // stop - start
            double _length;  // total branch length
            unsigned int n; // number of nodes (<= r)
            unsigned int r; // reserved number of nodes
            Node ** nodes;   // arrays of nodes (size = r)
            const Params * params;   // used only while recursing in mutate
            Random * random;   // used only while recursing in mutate
            unsigned int nmut;  // number of mutations occurring while a mutate() call
            DataHolder * data;   // used only while recursing in mutate

        private:

           /** \brief No default constructor available */
            Tree() {}
    };
}

#endif
