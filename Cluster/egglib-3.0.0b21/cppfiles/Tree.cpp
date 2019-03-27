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

#include "egglib.hpp"
#include "Tree.hpp"
#include "Random.hpp"
#include <cstdlib>
#include <new>
#include "Params.hpp"

namespace egglib {

    /* *********************** */

    Node::Node() {
        init();
    }

    /* *********************** */

    Node::~Node() {
        if (_reserved) {
            free(_mutationSite);
            free(_mutationPos);
        }
    }

    /* *********************** */

    Node::Node(const Node& src) {
        init();
        copy(src);
    }

    /* *********************** */

    Node& Node::operator=(const Node& src) {
        reset();
        copy(src);
        return *this;
    }

    /* *********************** */

    void Node::set_internal(double t, unsigned int son1, unsigned int son2) {
        this->_date = t;
        this->_son1 = son1;
        this->_son2 = son2;
        _is_terminal = false;
    }

    /* *********************** */

    void Node::set_terminal(double t, unsigned int label) {
        this->_date = t;
        this->_label = label;
        _is_terminal = true;
    }

    /* *********************** */

    bool Node::is_terminal() const {
        return _is_terminal;
    }

    /* *********************** */

    unsigned int Node::label() const {
        return _label;
    }

    /* *********************** */

    unsigned int Node::son1() const {
        return _son1;
    }

    /* *********************** */

    unsigned int Node::son2() const {
        return _son2;
    }

    /* *********************** */

    double Node::date() const {
        return _date;
    }

    /* *********************** */

    double Node::get_L() const {
        return _length;
    }

    /* *********************** */

    void Node::set_L(double value) {
        _length = value;
    }

    /* *********************** */

    unsigned int Node::nmut() const {
        return _nmut;
    }

    /* *********************** */

    unsigned int Node::mutationSite(unsigned int mut) const {
        return _mutationSite[mut];
    }

    /* *********************** */

    double Node::mutationPos(unsigned int mut) const {
        return _mutationPos[mut];
    }

    /* *********************** */

    void Node::addMutation(unsigned int site, double pos) {
        alloc(_nmut+1); // changes _nmut
        _mutationSite[_nmut-1] = site;
        _mutationPos[_nmut-1] = pos;
    }

    /* *********************** */

    void Node::reset() { // note there is a small redundancy with this method being followed by set_XXX normally
        _date = 0.;
        _length = 0.;
        _son1 = UNKNOWN;
        _son2 = UNKNOWN;
        _nmut = 0;
        _label = UNKNOWN;
        _is_terminal = true;
        // _reserved don't change
        // memory allocated under _mutationSite / _mutationPos is retained
    }

    /* *********************** */

    void Node::init() {
        reset();  // set other members
        _reserved = 0;
        _mutationSite = NULL;
        _mutationPos = NULL;
    }

    /* *********************** */

    void Node::copy(const Node& src) {
        _date = src._date;
        _length = src._length;
        _son1 = src._son1;
        _son2 = src._son2;
        _label = src._label;
        _is_terminal = src._is_terminal;
        alloc(src._nmut);  // changes _nmut
        for (unsigned int i=0; i<_nmut; i++) {
            _mutationSite[i] = src._mutationSite[i];
            _mutationPos[i] = src._mutationPos[i];
        }
    }

    /* *********************** */

    void Node::alloc(unsigned int nmut) {
        _nmut = nmut;
        if (_nmut > _reserved) {
            _mutationSite = (unsigned int *) realloc(_mutationSite, _nmut * sizeof(unsigned int));
            _mutationPos = (double *) realloc(_mutationPos, _nmut * sizeof(double));
            if (!_mutationSite || !_mutationPos) throw EGGMEM;
            _reserved = _nmut;
        }
    }

    /* *********************** */

    void Node::clear_mutations() {
        _nmut = 0;
    }

    /* Tree class */

    void Tree::init() {
        n = 0;
        r = 0;
        nodes = NULL;
        _length = 0.0;
    }

    /* *********************** */

    void Tree::free() {
        if (r > 0) {
            for (unsigned int i=0; i<r; i++) delete nodes[i];
            ::free(nodes);
        }
    }

    /* *********************** */

    void Tree::copy(const Tree& tree) {
        alloc_from(tree.n, tree.nodes);
        _length = tree._length;
        _start = tree._start;
        _stop = tree._stop;
        _cov = tree._cov;
    }

    /* *********************** */

    void Tree::realloc(unsigned int size) {

        unsigned int old = r;

        if (size > r) {

            // extend the table as needed

            nodes = (Node **) ::realloc(nodes, size * sizeof(Node *));
            if (!nodes) throw EGGMEM;

            // create the new objects

            for (unsigned int i=r; i<size; i++) {
                nodes[i] = new (std::nothrow) Node;
                if (!nodes[i]) throw EGGMEM;
            }
            r = size;
        }

        // reset recycled instances

        for (unsigned int i=n; i<old && i<size; i++) {
            nodes[i]->reset();
        }

        // record new value

        n = size;
    }

    /* *********************** */

    void Tree::alloc_from(unsigned int size, Node **source) {

        unsigned int old = r;

        if (size > r) {

            // extend the table as needed

            nodes = (Node **) ::realloc(nodes, size * sizeof(Node *));
            if (!nodes) throw EGGMEM;

            r = size;
        }

        // copy over existing nodes

        for (unsigned int i=0; i<old; i++) *(nodes[i]) = *(source[i]);

        // copy-create others

        for (unsigned int i=old; i<size; i++) {
            nodes[i] = new (std::nothrow) Node(*(source[i]));
            if (!nodes[i]) throw EGGMEM;
        }

        n = size;
    }

    /* *********************** */

    void Tree::set(unsigned int numberOfLeaves, double start, double stop) {
        _start = start;
        _stop = stop;
        _cov = stop - start;
        realloc(numberOfLeaves);
    }

    /* *********************** */

    Tree::Tree(unsigned int numberOfLeaves, double start, double stop) {
        init();
        set(numberOfLeaves, start, stop);
    }

    /* *********************** */

    Tree::Tree(const Tree& tree) {
        init();
        copy(tree);
    }

    /* *********************** */

    Tree& Tree::operator=(Tree tree) {
        copy(tree);
        return *this;
    }

    /* *********************** */

    Tree::~Tree() {
        free();
    }

    /* *********************** */

    void Tree::recomb(double point, Tree *new_tree) {

        // copy own nodes to the new_tree

        new_tree->alloc_from(n, nodes);

        // set up other members

        new_tree->_length = _length;
        new_tree->_start = point;
        new_tree->_stop = _stop;
        new_tree->_cov = _stop - point;

        // shrunk itself

        _stop = point;
        _cov = point - _start;
    }

    /* *********************** */

    unsigned int Tree::coal(unsigned int nodeIndex1, unsigned int nodeIndex2, double date) {

        // add one lineage (its index is n-1)

        realloc(n + 1);

        // set relationship and branch lengths

        nodes[n-1]->set_internal(date, nodeIndex1, nodeIndex2);
        nodes[nodeIndex1]->set_L( date - nodes[nodeIndex1]->date() );
        nodes[nodeIndex2]->set_L( date - nodes[nodeIndex2]->date() );

        // increment tree length

        _length += nodes[nodeIndex1]->get_L();
        _length += nodes[nodeIndex2]->get_L();

        return n - 1;
    }

    /* *********************** */

    unsigned int Tree::mutate(unsigned int site, DataHolder& data, const Params *params, Random *random) {

        nmut = 0;

        this->data = &data;
        this->params = params;
        this->random = random;

        // determine start allele

        unsigned int allele;

        if (params->get_random_start_allele()) {
            allele = random->irand( params->get_K() );
        }
        else {
            allele = 0;
        }

        // recurse (call directly both root's descendant, as the root node cannot have mutations (length 0)

        r_mutate(site, nodes[nodes[n-1]->son1()], allele);
        r_mutate(site, nodes[nodes[n-1]->son2()], allele);

        return nmut;
    }

    /* *********************** */

    double Tree::start() const {
        return _start;
    }

    /* *********************** */

    double Tree::stop() const {
        return _stop;
    }

    /* *********************** */

    double Tree::cov() const {
        return _cov;
    }

    /* *********************** */

    double Tree::L() const {
        return _length;
    }

    /* *********************** */

    Node* Tree::root() const {
        return nodes[n-1];
    }

    /* *********************** */

    unsigned int Tree::nnodes() const {
        return n;
    }

    /* *********************** */

    Node * const Tree::node(unsigned int i) const {
        return nodes[i];
    }

    /* *********************** */

    unsigned int Tree::addNode(double t, unsigned int label) {

        // claim new Node instance

        realloc(n + 1);

        // set it

        nodes[n-1]->set_terminal(t, label);

        // return index

        return n-1;
    }

    /* *********************** */

    void Tree::reset(unsigned int numberOfLeaves, double start, double stop) {
        n = 0;
        _length = 0.0;
        set(numberOfLeaves, start, stop);
    }

    /* *********************** */

    void Tree::clear_mutations() {
        for (unsigned int i=0; i<n; i++) nodes[i]->clear_mutations();
    }

    /* *********************** */

    void Tree::r_mutate(unsigned int site, Node *node, int allele) {

        // process all node's mutations

        bool mutated = false;

        for (unsigned int imut = 0; imut < node->nmut(); imut++) {

            // process the mutation only if it hits this site
            if (node->mutationSite(imut) == site) {

                // if IAM, one or more mutations are the same (only process the first)
                if (params->get_mutmodel() != Params::IAM || mutated == false) allele = next_allele(allele);

                nmut++;
                mutated = true;
            }
        }

        // if the node is terminal set, set DataHolder entries

        if (node->is_terminal()) {
            data->set_i(node->label(), site, allele);
        }

        // otherwise process descending nodes

        else {
            r_mutate(site, nodes[node->son1()], allele);
            r_mutate(site, nodes[node->son2()], allele);
        }
    }

    /* *********************** */

    int Tree::next_allele(int allele) {

        double X;
        unsigned int sign, step;

        switch (params->get_mutmodel()) {

            case (Params::KAM) :

                /* draw a random in the range defined by transition weights
                   note : in KAM the allele index is always positive */

                X = random->uniform() * params->get_transW_row((unsigned int) allele);

                // identify the picked allele

                for (unsigned int i=0; i<params->get_K(); i++) {

                    // skip the identity cell
                    if (i == (unsigned int) allele) continue;

                    X -= params->get_transW_pair((unsigned int) allele, i);
                    if (X<0) return i;
                }

                break; // this will ensure that the Exception is raised in case of problem


            case (Params::IAM) :

                return nmut;


            case (Params::SMM) :

                if (random->uniform()<0.5) return allele - 1;
                else return allele + 1;

                return allele + (random->uniform()<0.5) ? -1 : 1;


            case (Params::TPM) :

                // draw sign

                sign = (random->uniform()<0.5) ? -1 : 1;

                // draw step size

                if (random->uniform() < params->get_TPMproba()) {
                    step = random->grand(params->get_TPMparam());
                }
                else step = 1;

                return allele + sign * step;

        }

        throw EggRuntimeError("an unexpected error happened in Tree:nextallele - please report this bug");

        return -999999;
    }
}
