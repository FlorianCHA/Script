/*
    Copyright 2012 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_LINEAGE_HPP
#define EGGLIB_LINEAGE_HPP

namespace egglib {

  /** \brief Implements a single ARG lineage
    *
    * \ingroup coalesce
    *
    * This class implements a (recombinable and coalescable) lineage at
    * the ancestral recombination graph (ARG) level. It essentially
    * points to a node in each of the trees used to model the ARG.
    *
    * Header: <egglib-cpp/Lineage.hpp>
    *
    */
    class Lineage {

        public:

           /** \brief Constructor (no default constructor)
            *
            * Initializes the size of the node mapping array (node index
            * in each of the trees of the ARG). The array cells are not
            * initalized.
            *
            * \param ntrees number of Trees of the ARG (if the ARG has
            * more trees, excess ones are considered not to be covered
            * by this lineage.
            *
            */
            Lineage(unsigned int ntrees);

           /** \brief Destructor
            *
            */
            virtual ~Lineage();

           /** \brief Add a cell to the node mapping array
            *
            * \param node node index (with respect to the tree).
            * \param cov coverage of the tree (increment the total
            * coverage of the lineage).
            *
            */
            void addTree(unsigned int node, double cov);

           /** \brief Access node mapping cell
            *
            * This methods returns UNKNOWN if the considered tree is not
            * covered. The index must be in range.
            *
            */
            unsigned int get_node(unsigned int treeIndex) const;


           /** \brief Set node mapping cell
            *
            * \param index tree index.
            * \param node node index (with respect to the tree).
            * \param cov coverage of the tree (increment the total
            * coverage of the lineage).
            *
            * Indices must be valid, except that node must be UNKNOWN
            * for trees that are not covered.
            *
            */
            void set_node(unsigned int index, unsigned int node, double cov);

           /** \brief Get total coverage
            *
            * This is the sum of tree coverage (stop-start) for all
            * trees where this lineage is represented.
            *
            */
            double cov() const;

           /** \brief Reset the number of tree
            *
            * It allocates the node mapping array (only allocate if no
            * reserved memory available). This method does not
            * initialize the content of the table. It sets the coverage
            * to 0.
            *
            */
            void reset(unsigned int ntrees);

        protected:

            void init();  // like reset, but initalize reserved memory
            void alloc(unsigned int ntrees);

            unsigned int _ntrees;
            unsigned int _reserved;
            unsigned int *nodeMapping; //makes the connection between the nodes and the trees
            double _cov;

        private:

           /** \brief Default constructor not available*/
            Lineage(); //Default constructor (just to forbid default initialisation)

           /** \brief Copy constructor not available */
            Lineage(const Lineage& src);

           /** \brief Assignment operator not available */
            Lineage& operator=(const Lineage& src);
    };
}

#endif
