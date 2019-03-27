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

#ifndef EGGLIB_COALESCE_HPP
#define EGGLIB_COALESCE_HPP

#include "Params.hpp"
#include "DataHolder.hpp"
#include "Lineage.hpp"
#include "Tree.hpp"
#include "Random.hpp"

namespace egglib {

   /** \brief Coalescence simulation manager
    *
    * \ingroup coalesce
    *
    * The coalesce object must take a Params address. The Params object
    * will not modified. Coalesce objects cannot be created without
    * passing a valid Params address (NULL is out of question) and
    * cannot be copied either. Not that all action must be effected
    * using the simul method, including setting parameters. The other
    * methods such as coalescence, delayedSample and so on are meant to
    * be called by change types and should not be used directly. This
    * does not apply, of course, to DataHolder and tree getters.
    *
    * Header: <egglib-cpp/Coalesce.hpp>
    *
    */
    class Coalesce {

        public:

           /** \brief Constructor */
            Coalesce();

           /** \brief Destructor */
            ~Coalesce();

           /** \brief Perform a simulation
            *
            * \param params the address of a Params instance to
            * provide parameters. The Params object will be restored
            * after completion of the simulation.
            * \param random address of a random number generator.
            * \param mutate if this boolean is false, mutation phase
            * will be skipped and the DataHolder member will keep any
            * previous value (or will be empty otherwise). In all cases
            * the trees will be generated.
            *
            */
            void simul(Params *params, Random *random, bool mutate=true);

           /** \brief Perform mutation only
            *
            * The simul() method must have been called previous (with
            * the *mutate* boolean true or false). The DataHolder will
            * be regenerated.
            *
            */
            void mutate();

           /** \brief Perform a single coalescence event
            *
            * \param pop population index.
            * \param i index of the first lineage to coalesce.
            * \param j index of the second lineage to coalesce.
            *
            * It is not required that i<j. However it is required that
            * i!=j and that all indices are within their respective
            * ranges.
            *
            */
            void coalescence(unsigned int pop, unsigned int i, unsigned int j);

           /** \brief Perform a single migration event
            *
            * \param source source population.
            * \param i index of the lineage to migrate.
            * \param dest destination population.
            *
            * It is required that all indices are within their repestive
            * ranges.
            *
            */
            void migrate(unsigned int source, unsigned int i, unsigned int dest);

           /** \brief Admixture
            *
            * Instant migration of lineages from a given population to
            * another.
            *
            * \param source index of the population providing migrants.
            * \param dest index of the population accepting migrants.
            * \param proba migration probability.
            *
            * It is legal to provide a probability equal to 0 or 1, and
            * to use a source population which does not actually have
            * migrants. It is not leage to use equal source and dest
            * indexes.
            *
            */
            void admixt(unsigned int source, unsigned int dest, double proba);

           /** \brief Bottleneck
            *
            * Perform coalescence in the specified population during
            * a given amount of time.
            *
            */
            void bottleneck(unsigned int pop, double duration);

           /** \brief Number of trees
            *
            * This method provides access to the number of tree, or
            * number of recombination segments (equal to the number
            * of recombination events plus one). If no recombination
            * event has occurred, the value is 1. If no simulation has
            * been performed, the value is 0.
            *
            */
            unsigned int number_of_trees() const;

           /** \brief Add a delayed sample to the simulation
            *
            * This method ought to be called by the class DelayedSample.
            *
            * \param date date of the event (should be larger than the
            * current date).
            * \param pop index of the population (should be smaller
            * than the current number of populations).
            * \param n1 number of single samples.
            * \param n2 number of double samples.
            *
            */
            void delayedSample(double date, unsigned int pop, unsigned int n1, unsigned int n2);

           /** \brief Get a tree
            *
            * This method provides access to a given tree. The passed
            * index must be smaller than the value returned by the
            * method number_of_trees(). Obviously, a simulation must
            * have been performed. The returned pointer refers to an
            * instance which is stored within the Coalesce object and
            * managed by it. The object can not be modified by any
            * means. Out of range values might result in program crash
            * or even worse.  Beware that the address to the internally
            * stored object is returned, and the address must be
            * considered invalid after the next call to the simul()
            * method.
            *
            */
            Tree const * const tree(unsigned int i) const;

           /** \brief Get the start position of a tree
            *
            */
            double tree_start(unsigned int i) const;

           /** \brief Get the stop position of a tree
            *
            */
            double tree_stop(unsigned int i) const;

           /** \brief Get simulated genotypes
            *
            * Get the DataHolder (as a matrix object) corresponding to
            * the last simulated data set. Only meaningful if the
            * simul() method has been called with mutate=true.
            * Otherwise, the returned instance will correspond to the
            * previous simulations, or will be empty. In all cases,
            * there will be two levels of structure, the first
            * corresponding to populations (from 0 to k-1) where k is
            * the number of populations, and the second corresponding to
            * individual indices (ranging from 0 to the total number of
            * sampled individuals over all populations, minus one). An
            * individual is represented by either one or two consecutive
            * samples. The allele values are depending on the mutation
            * model. For FSS, they are restricted to the range [0, K-1]
            * where K is the number of allowed alleles, while for other
            * models unrestricted positive and negative values are
            * allowed. In the infinitely many sites model (that is, if
            * the sites parameter is 0), only the variable sites are
            * exported, and all exported sites are variables. If the
            * number of sites is finite, all sites are exported and
            * there can be invariable sites. Beware that the address to
            * the internally stored object is returned, and the address
            * must be considered invalid after the next call to the
            * simul method.
            *
            */
            DataHolder const * const data() const;

           /** \brief Get the position of a mutation
            *
            * This method must not be called if the number of sites was
            * fixed to a non-zero value in the Params instances used to
            * configure simulations. This method returns the position of
            * a given mutation when the mutation was randomly drawn. It
            * is guaranteed that positions are in increasing order. It
            * is required that the index is smaller than the number of
            * mutations drawn from the last simulation.
            *
            */
            double site_position(unsigned int mut) const;

        protected:

            void alloc_pop();   // update the Lineage* and counter arrays
            void alloc_pop(unsigned int pop, unsigned int n);
            void add_one_lineage(unsigned int pop); // increment pops array but don't create lineage
            unsigned int alloc_stack(unsigned int incr);   // create new Lineage objects as needed / new objects have ntrees trees
            void alloc_one_tree();   // the tree must be reset in all cases
            ////-void alloc_mut_mapping(unsigned int nmut); // allocates nmut positions, and tree addresses

            void diploid();
            void label();
            void tcoal();
            double tcoal(unsigned int pop);
            void tmigr();
            void trec();
            void tevent();
            void do_coal();
            void do_migr();
            void do_rec();

            char nextW;
            double nextT;
            unsigned int nextP;
            double nextM;

            DataHolder _data;
            Params * params;
            Random * random;

            unsigned int npop;
            unsigned int npop_r;
            double * crec;
            unsigned int * popsize;
            unsigned int * popsize_r;
            Lineage *** pops;

            unsigned int stack;
            unsigned int stack_r;
            Lineage ** lineages;

            double time;
            unsigned int ns; // lineages at start
            unsigned int remaining;  // remaining lineages

            unsigned int ntrees;
            unsigned int ntrees_r;

            unsigned int site_mutation_c;
            unsigned int * site_mutation; // UNKNOWN if no mutation, index otherwise
            unsigned int mutation_site_c;
            unsigned int * mutation_site;
            unsigned int site_tree_c;
            Tree ** site_tree;
            unsigned int site_pos_c;
            double * site_pos;

            Tree ** trees;

        private:

            /** \brief Copy constructor not available */
            Coalesce(const Coalesce& coalesce) {}

            /** \brief Copy assignment operator not available */
            Coalesce& operator=(const Coalesce& coalesce) {return *this;}
    };
}

#endif
