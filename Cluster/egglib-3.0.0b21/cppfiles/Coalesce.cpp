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

#include "Coalesce.hpp"
#include "egglib.hpp"
#include <new>
#include <cstdlib>
#include <cmath>

namespace egglib {

    // ************************************* Constructor / destructor //

    Coalesce::Coalesce() {

        time = 0.;

        npop = 0;
        npop_r = 0;
        crec = NULL;
        popsize = NULL;
        popsize_r = NULL;
        pops = NULL;

        stack = 0;
        stack_r = 0;
        lineages = NULL;

        ntrees = 0;
        ntrees_r = 0;
        trees = NULL;

        site_mutation_c = 0;
        site_mutation = NULL;
        mutation_site_c = 0;
        mutation_site = NULL;
        site_tree_c = 0;
        site_tree = NULL;
        site_pos_c = 0;
        site_pos = NULL;

        _data.set_is_matrix(true);
    }

    // *********** //

    Coalesce::~Coalesce() {

        if (npop_r) {
            free(crec);
            for (unsigned int i=0; i<npop_r; i++) {
                if (pops[i]) free(pops[i]);
            }
            free(pops);
            free(popsize);
            free(popsize_r);
        }

        if (stack_r) {
            for (unsigned int i=0; i<stack_r; i++) {
                delete lineages[i];
            }
            free(lineages);
        }

        if (ntrees_r) {
            for (unsigned int i=0; i<ntrees_r; i++) {
                delete trees[i];
            }
            free(trees);
        }

        if (site_mutation) free(site_mutation);
        if (mutation_site) free(mutation_site);
        if (site_tree) free(site_tree);
        if (site_pos) free(site_pos);
    }

    // *********************************************** Main function //

    void Coalesce::simul(Params *params, Random *random, bool mutate) {

        this->params = params;
        this->random = random;

        // setup

        time = 0.;
        ntrees = 0;
        ns = 0;
        npop = 0;
        stack = 0;

        // allocate memory

        alloc_one_tree();  // ntrees++
        alloc_pop(); // also setup npop, ns, crec and popsize
        remaining = ns;

        // raise an error if ns too small
        if (ns < 2) {
            bool flag = false;
            if (params->nDSChanges() > 0) {
                const Event * e = params->firstChange();
                while (e) {
                    if (e->event_type() == Event::delayed && (e->get_number1() > 0 || e->get_number2() > 0)) {
                        flag = true;
                        break;
                    }
                }
            }
            if (!flag) throw EggArgumentValueError("at least 2 samples overall are required");
        }

        trees[0]->reset(ns, 0, 1);

        // setup the lineages (cover the unique tree)

        unsigned int c = 0;

        for (unsigned int i=0; i<npop; i++) {

            for (unsigned int j=0; j<popsize[i]; j++) {
                pops[i][j]->set_node(0, c, 1.);
                trees[0]->node(c)->set_terminal(0.0, c);
                c++;
            }
        }

        // diploid lineages

        diploid();

        // main loop

        unsigned int iter_counter=0;
        while (remaining > 1 || params->nDSChanges() > 0) {

            // get event time

            nextW = '0';
            nextT = UNDEF;

            tcoal();
            tmigr();
            trec();
            tevent();

            if (nextW == '0') throw EggRuntimeError("infinite coalescent time (unconnected populations or excessive ancestral population size)");

            // apply event

            time += nextT;

            switch (nextW) {

                case 'c':
                    do_coal();
                    break;

                case 'm':
                    do_migr();
                    break;

                case 'r':
                    do_rec();
                    break;

                case 'e':
                    params->nextChangeDo(this);
                    break;

                default:
                    throw EggRuntimeError("infinite coalescent event - please report this bug");
            }

            // prevent infinite loop
            if (++iter_counter > params->get_max_iter()) {
                throw EggRuntimeError("failed to complete coalescent tree: two lineages might be trapped to unconnected populations (if you are sure your model is correct, increase the parameter `max_iter`)");
            }
        }

        // mutate if needed

        if (mutate) Coalesce::mutate();

        // restore parameters holder

        params->restore();
    }

    // ******************************************** Allocator helpers //

    void Coalesce::alloc_pop() {

        // get new number of pop
        unsigned int new_npop = params->k();

        // if number of pop increases (compared to cache), extend array
        if (new_npop > npop_r) {

            // actually extend arrays
            pops = (Lineage ***) realloc(pops, new_npop * sizeof(Lineage **));
            if (pops == NULL) throw EGGMEM;

            popsize = (unsigned int *) realloc(popsize, new_npop * sizeof(unsigned int));
            if (popsize == NULL) throw EGGMEM;

            crec = (double *) realloc(crec, new_npop * sizeof(double));
            if (popsize == NULL) throw EGGMEM;

            popsize_r = (unsigned int *) realloc(popsize_r, new_npop * sizeof(unsigned int));
            if (popsize_r == NULL) throw EGGMEM;

            // initialize new data
            for (unsigned int i=npop_r; i<new_npop; i++) {
                popsize_r[i] = 0;
                pops[i] = NULL;
            }

            // update reserve counter
            npop_r = new_npop;
        }

        // initialize all population to (officially) 0
        for (unsigned int i=0; i<new_npop; i++) {
            popsize[i] = 0;
        }

        // allocate sub-arrays
        unsigned int n;
        for (unsigned int i=npop; i<new_npop; i++) {

            // get new pop size
            n = 2 * params->get_n2(i) + params->get_n1(i);
            ns += n;

            // set it
            alloc_pop(i, n);
            crec[i] = n;
        }

        // update operational counter
        npop = new_npop;
    }

    // *********** //

    void Coalesce::alloc_pop(unsigned int pop, unsigned int n) {

        // if pop size increased (compared to cache), extend array
        if (n > popsize_r[pop]) {
            pops[pop] = (Lineage **) realloc(pops[pop], n * sizeof(Lineage *));
            if (pops[pop] == NULL) throw EGGMEM;

            // update reserve counter
            popsize_r[pop] = n;
        }

        // claim Lineage objects
        unsigned int x = alloc_stack(n - popsize[pop]);

        // populate new entries with addresses
        for (unsigned int i=0; i<n-popsize[pop]; i++) {
            pops[pop][popsize[pop] + i] = lineages[x+i];
        }

        // set the newly claimed objects to the right number of trees
            // (redundant for newly created Lineage but we cannot tell)
            // (and it is not actually costy)
        for (unsigned int i=popsize[pop]; i<n; i++) {
            pops[pop][i]->reset(ntrees);
        }

        // update operational counter
        popsize[pop] = n;
    }

    // *********** //

    void Coalesce::add_one_lineage(unsigned int pop) {

        // need real allocation

        if (popsize[pop] + 1 > popsize_r[pop]) {

            pops[pop] = (Lineage **) realloc(pops[pop], (popsize[pop] + 1) * sizeof(Lineage *));
            if (!pops[pop]) throw EGGMEM;

            popsize_r[pop] = popsize[pop] + 1;
        }

        // update operational counter

        popsize[pop]++;
    }

    // *********** //

    unsigned int Coalesce::alloc_stack(unsigned int incr) {

        // if allocation needed for all or part of new objects

        if (stack + incr > stack_r) {

            // allocate the table

            lineages = (Lineage **) realloc(lineages, (stack + incr) * sizeof(Lineage *));
            if (lineages == NULL) throw EGGMEM;

            // create the new objects as needed

            for (unsigned int i=stack_r; i<stack + incr; i++) {
                lineages[i] = new (std::nothrow) Lineage(ntrees);
                if (lineages[i] == NULL) throw EGGMEM;
            }

            // update reserved counter

            stack_r = stack + incr;
        }

        // update operation counter
        stack += incr;

        // return index of the first new object
        return stack - incr;
    }

    // *********** //

    void Coalesce::alloc_one_tree() {

        // realloc the Tree* and create the object if needed

        if ((ntrees+1) > ntrees_r) {
            trees =  (Tree **) realloc(trees, (ntrees+1) * sizeof(Tree*));
            if (trees == NULL) throw EGGMEM;

            trees[ntrees] = new (std::nothrow) Tree(0,0,0);
            if (trees[ntrees] == NULL) throw EGGMEM;

            ntrees_r = ntrees + 1;
        }

        // otherwise just recycle an old tree

        else {
            trees[ntrees]->reset(0,0,0);
        }

        ntrees++;
    }

    // ********************************************** heavy modifiers //

    void Coalesce::diploid() {
        unsigned int c;
        double CRIT;
        for (unsigned int i=0; i<npop; i++) {
            c = 0;
            CRIT = params->get_s(i)/(2.-params->get_s(i));
            for (unsigned int j=0; j<params->get_n2(i); j++) {
                if (random->uniform()<CRIT) {
                    coalescence(i, c, c+1);
                }
                else {
                    c += 2;
                }
            }
        }
    }

    // *********** //

    void Coalesce::admixt(unsigned int source, unsigned int dest, double proba) {
        unsigned int i=0;
        while (i<popsize[source]) {
            if (random->uniform() < proba) {
                migrate(source, i, dest);
            }
            else i++;
        }
    }

    // *********** //

    void Coalesce::bottleneck(unsigned int pop, double duration) {

        // define local variables

        double t = 0.0;
        double incr = 0.0;
        unsigned int i1, i2;

        // iterate not more than duration

        while (true) {

            // pick a time

            incr = tcoal(pop);
            if (incr == UNDEF) break;  // no possibility of coalescence
            if (t + incr > duration) break; // done

            // pick two lineages

            i2 = i1 = random->irand(popsize[pop]);
            while (i1 == i2) {
                i2 = random->irand(popsize[pop]);
            }

            // perform coalesce

            coalescence(pop, i1, i2);

            // increment only local time

            t += incr;
        }
    }

    // *********************************************** mutation function //

    void Coalesce::mutate() {

        double X;         // for random values
        double Wt = 0.0;  // sum of weights * tree length at each sites (only if nsites>0)

        // reset trees (in case mutations are already stored)

        for (unsigned int i=0; i<ntrees; i++) {
            trees[i]->clear_mutations();
        }

        // total branch length

        double ARG_length = 0;
        for (unsigned int i=0; i<ntrees; i++) {
            ARG_length += trees[i]->L() * trees[i]->cov();
        }

        // defines the number of mutations

        unsigned int nmut;
        if (params->get_theta()==0.) {
            nmut = params->get_fixed();
        }
        else {
            nmut = random->prand(params->get_theta() * ARG_length);
        }

        /*
         * Affects mutation to positions in the case of fixed number of
         * sites.
         *
         * Warning: Params::L() and Tree::L() are different parameters.
         *
         */

        unsigned int nsites = params->get_L();

        if (nmut > 0 && nsites > 0) {

            // allocate the array to store Tree addresses of each site

            if (nsites > site_tree_c) {
                site_tree = (Tree **) realloc(site_tree, nsites * sizeof(Tree *));
                if (!site_tree) throw EGGMEM;
                site_tree_c = nsites;
            }

            unsigned int c;   // for current site

            /*
             * Record all tree's total length.
             * Also record correspond tree's address (will then be
             * mapped to a tree mapping).
             *
             */

            for (unsigned int tr=0; tr < ntrees; tr++) {

                // process all sites falling in this segment
                c = 0;

                while (c < nsites && params->get_sitePos(c) < trees[tr]->start()) c++;

                if (c == nsites) break;
                    // move on until site is after start of this tree's segment

                // while current site is within this tree's segment
                while (c < nsites && params->get_sitePos(c) <= trees[tr]->stop()) {

                    // record this site's tree
                    Wt += params->get_siteW(c) * trees[tr]->L();
                    site_tree[c] = trees[tr];

                    // next site
                    c++;
                }
            }

            /*
             * Draw a site for each mutation. Records the tree, index
             * and position corresponding to each mutation.
             *
             */

            // allocate the arrays for mapping mutations to sites and all the way around

            if (nsites > site_mutation_c) {
                site_mutation = (unsigned int *) realloc(site_mutation, nsites * sizeof(unsigned int));
                if (!site_mutation) throw EGGMEM;
                site_mutation_c = nsites;
            }
            for (unsigned int i=0; i<nsites; i++) site_mutation[i] = UNKNOWN; // not UNKNOWN only if this site is mutated

            if (nmut > mutation_site_c) {
                mutation_site = (unsigned int *) realloc(mutation_site, nmut * sizeof(unsigned int));
                if (!mutation_site) throw EGGMEM;
                mutation_site_c = nmut;
            } // no need to initialize (all mutations have a site)

            for (unsigned int imut=0; imut<nmut; imut++) {

                X = random->uniform() * Wt;

                for (unsigned int sit=0; sit<nsites; sit++) {
                    X -= params->get_siteW(sit) * site_tree[sit]->L();
                    if (X < 0) {                     // this assumes that uniform() is [0; 1)
                        site_mutation[sit] = imut;
                        mutation_site[imut] = sit;
                        break;
                    }
                }
            }
            // from now on site_mapping_tree is not needed
        }

        /*
         * Affects mutation to positions in the case of infinite number
         * of sites. First, the position of mutations in the [0,1]
         * interval are drawn. It is way in such a way that they are in
         * order (the intervals between them are drawn). Then, map a
         * tree to each mutation.
         *
         */

        if (nmut > 0 && nsites == 0) {

            // allocate the array to store Tree addresses of each mutation (therefore each site)

            if (nmut > site_tree_c) {
                site_tree = (Tree **) realloc(site_tree, nmut * sizeof(Tree *));
                if (!site_tree) throw EGGMEM;
                site_tree_c = nmut;
            }

            // and the array for site positions

            if (nmut > site_pos_c) {
                site_pos = (double *) realloc(site_pos, nmut * sizeof(double));
                if (!site_pos) throw EGGMEM;
                site_pos_c = nmut;
            }

            // draw S+1 intervals and store cumul

            double c = 0.;
            for (unsigned int i=0; i<nmut; i++) {
                c += random->uniform();
                site_pos[i] = c;
            }

            // draw a last interval (between last mutation and end of chromosome)

            c += random->uniform();

            // rescale to [0;1]

            for (unsigned int i=0; i<nmut; i++) site_pos[i] /= c;

            /*
             * Map a tree to each mutation. The principle
             * is the same as for mapping trees to sites.
             *
             */

            unsigned int mut = 0;
            unsigned int tr = 0;

            while (mut < nmut) {

                // move forward until current mut is within current tree's interval

                while (site_pos[mut] < trees[tr]->start() || site_pos[mut] >= trees[tr]->stop()) {
                    tr++;
                    if (tr == ntrees) tr = 0;
                }

                // store addresses while current mut is within current tree's interval

                while (mut < nmut && (tr == ntrees || site_pos[mut] < trees[tr]->stop())) {
                                // note: the first term forces the last mutation(s) to be within the last tree
                                // this is expect to occur exceptionnally (pos = 1) or through errors
                    site_tree[mut] = trees[tr];
                    mut++;
                }
            }
        }

        /*
         * Mutations are now placed on the chromosome. We have the
         * following arrays available whatever nsites:
         *      site_tree (Tree*)
         *      site_mutation (unsigned int) ; only if nsites>0
         *      site_pos (double) ; only if nsite=0
         *
         */

        bool site_flag; // true if all sites are mutated by definition
        if (nsites == 0) {
            nsites = nmut;
            site_flag = true;
        }
        else {
            site_flag = false;
        }

        /*
         * Now, place the mutations on tree nodes (each node has an array)
         *
         */

        Tree *tree;
        unsigned int nnodes;

        for (unsigned int imut=0; imut<nmut; imut++) {

            // the tree

            tree = site_flag ? site_tree[imut] : site_tree[mutation_site[imut]];
            nnodes = tree->nnodes();

            // draw a random in tree length coverage range

            X = random->uniform() * tree->L();

            // find the node corresponding to this X

            for (unsigned int n=0; n<nnodes; n++) {

                X -= tree->node(n)->get_L();

                if (X<0) {    // the < (instead <=) assumes 1 cannot be drawn

                    // add the mutation

                    tree->node(n)->addMutation(site_flag ? imut : mutation_site[imut],
                                               site_flag ? site_pos[imut] : params->get_sitePos(mutation_site[imut]));
                    break;
                }
            }
        }

        /*
         * Apply mutations. This is mostly done by the dedicated method
         * of Tree. The loop is over sites ( = nmum if sites is 0). The
         * DataHolder is not initialized because the method necessarily
         * sets all entries of a given column
         *
         */

        _data.set_nsam_i(ns);
        _data.set_nsam_o(0);
        _data.set_ngroups(2);
        _data.set_nsit(nsites);

        double check_site_pos = 0.0; // used to check that site positions are in order

        for (unsigned sit = 0; sit < nsites; sit++) {

            // check that site position are in increasing order
            if (!site_flag) {
                if (params->get_sitePos(sit) < check_site_pos) throw EggArgumentValueError("site positions are not sorted");
                check_site_pos = params->get_sitePos(sit);
            }

            /*
             * Perform mutations for all mutated sites
             *                              = all sites if nsites == 0
             *                              = flagged sites otherwise
             *
             */

            if (site_flag || (nmut>0 && site_mutation[sit] != UNKNOWN)) {
                site_tree[sit]->mutate(sit, _data, params, random);
            }

            /*
             * The mutation counter is incremented by the number of
             * mutations found at this site.
             *
             */

            /*
             * If no mutations, it is necessary to set alleles forcily
             * 0 by default, but a random value if KAM is used with the
             * random start allele option.
             *
             */

            else {
                int all = (params->get_random_start_allele()) ? random->irand(params->get_K()) : 0;
                for (unsigned int i=0; i<ns; i++) _data.set_i(i, sit, all);
             }
        }

        label();
    }

    // *********** //

    void Coalesce::label() {

        // data has already 2 group levels (constructor)

        unsigned int c = 0;	   // counter for samples
        unsigned int I = 0;    // counter for individuals

        // process all pops

        for (unsigned int pop = 0; pop < params->k(); pop++) {

            // label doubles

            for (unsigned int i = 0; i < params->get_n2(pop); i++) {
                _data.set_group_i(c, 0, pop);
                _data.set_group_i(c, 1, I);
                c++;
                _data.set_group_i(c, 0, pop);
                _data.set_group_i(c, 1, I);
                c++;
                I++;
            }

            // label singles

            for (unsigned int i = 0; i < params->get_n1(pop); i++) {
                _data.set_group_i(c, 0, pop);
                _data.set_group_i(c, 1, I);
                c++;
                I++;
            }
        }

        // label delayed samples

        Event * cur = params->firstChange();
        while (cur != NULL) {
            if (cur->event_type() == Event::delayed) {

                // label doubles

                for (unsigned int i = 0; i < cur->get_number2(); i++) {
                    _data.set_group_i(c, 0, cur->get_dest());
                    _data.set_group_i(c, 1, I);
                    c++;
                    _data.set_group_i(c, 0, cur->get_dest());
                    _data.set_group_i(c, 1, I);
                    c++;
                    I++;
                }

                // label singles

                for (unsigned int i = 0; i < cur->get_number1(); i++) {
                    _data.set_group_i(c, 0, cur->get_dest());
                    _data.set_group_i(c, 1, I);
                    c++;
                    I++;
                }
            }
            cur = cur->next();
        }

        /** todo: find a way to skip the above if the params hasn't changed */

        // send the designed samples to the outgroup

        if (params->get_outgroup() != MISSING) {
            unsigned int lbl = params->get_outgroup();
            unsigned int num = _data.get_nsam_i();
            for (unsigned int idx=0; idx<num; ) {
                if (lbl == _data.get_group_i(idx, 0)) {
                    _data.to_outgroup(idx, _data.get_group_i(idx, 1));
                    num--;
                }
                else idx++;
            }
        }
    }

    // ********************************************** get events time //

    void Coalesce::tcoal() {

        double t;
        for (unsigned int i=0; i<npop; i++) {

            t = tcoal(i);
            if (t == UNDEF) continue;
            if (nextT == UNDEF || t < nextT) {
                nextT = t;
                nextW = 'c';
                nextP = i  ;
            }
        }
    }

    // *********** //

    double Coalesce::tcoal(unsigned int pop) {

        unsigned int n = popsize[pop];
        if (n<2) return UNDEF;
        double expect = params->get_N(pop)*(2.-params->get_s(pop)) / (2.*n*(n-1));
        double t;
        if (params->get_G(pop)==0.) t = random->erand(expect);
        else {
            double arg  = (1. + params->get_G(pop) * exp(-params->get_G(pop) *
                      (time - params->lastChange(pop))) * random->erand(expect));
            if (arg > 0) t = log(arg)/params->get_G(pop);
            else return UNDEF;
        }

        return t;
    }

    // *********** //

    void Coalesce::trec() {

        if (params->get_R() == 0.) return;

        double t;
        for (unsigned int i=0; i<npop; i++) {

            // because of rounding error somewhere,
            // I sometimes have crec (slightly) < 0 for empty deme
            // screwing up my nice algorithm
            if (popsize[i]==0) continue;

            double R = params->get_R() * (1 - params->get_s(i)/(2.-params->get_s(i)));

            if (R==0.) continue;
            R *= 1. * crec[i];
            t = random->erand(1./R);

            if (nextT == UNDEF && params->nextChangeDate() == UNDEF) {
                throw EggRuntimeError("infinite coalescent time: simulation will never complete (unconnected populations or excessive ancestral population size: the only possible event is recombination and no demographic change planned)");
            }

            if (nextT == UNDEF || t < nextT) {
                nextT = t;
                nextW = 'r';
                nextP = i;
            }
        }
    }

    // *********** //

    void Coalesce::tmigr() {

        if (npop<2) return;

        nextM = 0.;
        for (unsigned i=0; i<npop; i++) {
            nextM += popsize[i] * params->M().get_row(i);
        }

        if (nextM < 0.00000000000000001) return;

        double t = random->erand(1./nextM);

        if (nextT == UNDEF || t < nextT) {
            nextT = t;
            nextW = 'm';
        }
    }

    // ********** //

    void Coalesce::tevent() {

        if (params->nextChangeDate() == UNDEF) return;

        double te = params->nextChangeDate();

        te -= time;
        if (te < -0.000000000000001) throw EggRuntimeError("negative time to next change: invalid change date or invalid order of changes");
        if (te < 0) te = 0;
        if (nextT == UNDEF || te < nextT) {
            nextW = 'e';
            nextT = te;
        }
    }

    // ************************************ applies coalescent events //

    void Coalesce::do_coal() {

        // draw a lineage in the coalescing pop
        unsigned int a = random->irand(popsize[nextP]);

        // draw a second (different) lineage
        unsigned int b;
        do b = random->irand(popsize[nextP]);
        while (a==b);

        // coalesce them
        coalescence(nextP, a, b);
    }

    // ********** //

    void Coalesce::do_migr() {

        double x;

        // draw one pair of populations

        x = random->uniform() * nextM;

        for (unsigned i=0; i<npop; i++) {
            for (unsigned int j=0; j<npop; j++) {

                if (i==j) continue;
                x -= popsize[i] *  params->M().get_pair(i, j);

                // picked a pair

                if (x < 0) {

                    // pick a lineage and migrate it

                    unsigned int lin = random->irand(popsize[i]);
                    migrate(i, lin, j);
                    return;
                }
            }
        }
    }

    // ********** //

    void Coalesce::do_rec() {

        // pick a recombination point (lineage + point in chromosome)

        unsigned int lin1 = UNKNOWN;
        unsigned int treeIndex = UNKNOWN;
        double point = random->uniform()*crec[nextP];

        for (unsigned int i=0; i<popsize[nextP]; i++) {
            for (unsigned int j=0; j<ntrees; j++) {
                if (pops[nextP][i]->get_node(j) != UNKNOWN) {
                    point -= trees[j]->cov();
                    if (point < 0) {
                        lin1 = i;
                        treeIndex = j;
                        break;
                    }
                }
            }
            if (point < 0) break;
        }

        if (lin1 == UNKNOWN) throw EggRuntimeError("bug in Coalesce::do_rec - please report this error");

        point += trees[treeIndex]->stop();  // converts the (negative) remainder of X to a position

        // claim a new tree

        alloc_one_tree();   // the new tree is at index ntrees-1

        // make the recombination

        trees[treeIndex]->recomb(point, trees[ntrees-1]);

        // add the new tree to all lineages

        for (unsigned int i=0; i<npop; i++) {
            for (unsigned int j=0; j<popsize[i]; j++) {
                pops[i][j]->addTree(
                    pops[i][j]->get_node(treeIndex), 0.); // the node index is the same as for the old tree
                                                          // note that there is no coverage increment
            }
        }

        // create a new lineage and add it to the population

        unsigned int lin2 = alloc_stack(1);  // index of new lineage in lineages array
        lineages[lin2]->reset(ntrees);
        add_one_lineage(nextP);

        pops[nextP][ popsize[nextP] - 1 ] = lineages[lin2];
        lin2 = popsize[nextP] - 1;  // index in pops[] (as for lin1)

        // map trees to these two lineages

        for (unsigned int i=0; i<ntrees; i++) {

            // what is not covered in the original lineage will not be covered by new lineage

            if (pops[nextP][lin1]->get_node(i) == UNKNOWN) {
                pops[nextP][lin2]->set_node(i, UNKNOWN, 0.);
                continue;
            }

            // what is at the left of point is for lin1
            // set it as non-covered for lin2 (remains in lin1)

            if (trees[i]->stop() <= point) {
                pops[nextP][lin2]->set_node(i, UNKNOWN, 0.);
            }

            // what is at the right of point is for lin2
            // set it as covered for lin2 (takes the node of lin1)
            // set it as non-covered for lin1 (need to decrement the coverage)

            else {
                pops[nextP][lin2]->set_node(i,
                        pops[nextP][lin1]->get_node(i), trees[i]->cov());
                pops[nextP][lin1]->set_node(i, UNKNOWN, - trees[i]->cov());
            }
        }

        // one more lineage, damn it!

        remaining++;

        // (note that crec does not change at all)
    }

    // ************************* coalescent structure main operations //

    void Coalesce::coalescence(unsigned int pop, unsigned int i, unsigned int j) {

        // claim a new lineage (and get its address)

        unsigned int pos = alloc_stack(1);
        Lineage* parent = lineages[pos]; // merging these two lines led to a segmentation fault
        parent->reset(ntrees);

        // get reference of two coalesced lineages

        Lineage* son1 = pops[pop][i];
        Lineage* son2 = pops[pop][j];

        // coalesce (if needed) tree branches for all trees

        for (unsigned int k=0; k<ntrees; k++) {

            // get nodes of the lineages for that tree

            unsigned int n1 = son1->get_node(k);
            unsigned int n2 = son2->get_node(k);

            // only if both nodes are present, coalesce them

            if (n1 != UNKNOWN && n2 != UNKNOWN) {

                // coalesce nodes of this tree / set result for parent

                parent->set_node( k, trees[k]->coal(n1, n2, time), trees[k]->cov() );

                // reduces recombination opportunities

                if (crec[pop] > 0) {
                    crec[pop] -= trees[k]->cov();
                }
            }

            // otherwise the parent takes the node of whichever lineage is covered

            else {

                if (n1 != UNKNOWN) parent->set_node(k, n1, trees[k]->cov());
                else {
                    if (n2 != UNKNOWN) parent->set_node(k, n2, trees[k]->cov());
                    else parent->set_node(k, UNKNOWN, 0.);
                }
            }
        }

        // remove the two lineages from pop (but they are not moved away from lineages)

        unsigned int k, r;
        for (k=0, r=0; /* ! */ ; k++, r++) {
            while (r==i || r==j) r++;  // second curser goes to first non-removable lineage
            if (r == popsize[pop]) break; // if none, quit
            if (r>k) pops[pop][k] = pops[pop][r];  // copies only greater addresses
        }

        // reduce virtually the array size and add the parent

        popsize[pop]--;
        pops[pop][popsize[pop]-1] = parent;

        // decrement general counter

        remaining--;

    }

    // ********** //

    void Coalesce::migrate(unsigned int source, unsigned int i, unsigned int dest) {

        // update recombination coefficients

        double cov = pops[source][i]->cov();

        crec[source] -= cov;
        crec[dest] += cov;

        // add the lineage to its new population

        add_one_lineage(dest); // memory allocation
        pops[dest][popsize[dest]-1] = pops[source][i];

        // remove the  lineage from its old population

        if (i < (popsize[source] - 1)) {   // no need to do anything if lineage is last
            for (unsigned int j=i+1; j<popsize[source]; j++) {  // iterator is the next lineage
                pops[source][j-1] = pops[source][j];
            }
        }

        // update (virtually) the pop array size

        popsize[source]--;
    }

    // ********** //

    void Coalesce::delayedSample(double date, unsigned int pop, unsigned int n1, unsigned int n2) {

        // increase counters

        unsigned int n = 2 * n2 + n1;
        remaining += n;

        crec[pop] += n;

        // record the first index of new lineages

        unsigned int index = popsize[pop];

        // increase population size (also creates lineages)

        alloc_pop(pop, popsize[pop]+n);

        // add each new lineage to all trees (creates new nodes as needed)

        for (unsigned int i=0; i<n; i++) {
            for (unsigned int j=0; j<ntrees; j++) {
                pops[pop][index+i]->set_node(j, trees[j]->addNode(date, ns), trees[j]->cov());
            }
            ns++;
        }

        // manage diploid samples

        double CRIT = params->get_s(pop)/(2.-params->get_s(pop));
        unsigned int c = 0;
        for (unsigned int i=0; i<n2; i++) {
            if (random->uniform()<CRIT) coalescence(pop, c, c+1);
            else c += 2;
        }
    }

    // **************************************************** accessors //

    unsigned int Coalesce::number_of_trees() const {
        return ntrees;
    }

    // ********** //

    Tree const * const Coalesce::tree(unsigned int i) const {
        return trees[i];
    }

    // ********** //

    double Coalesce::tree_start(unsigned int i) const {
        return trees[i]->start();
    }

    // ********** //

    double Coalesce::tree_stop(unsigned int i) const {
        return trees[i]->stop();
    }

    // ********** //

    DataHolder const * const Coalesce::data() const {
        return &_data;
    }

    // ********** //

    double Coalesce::site_position(unsigned int mut) const {
        return site_pos[mut];
    }
}
