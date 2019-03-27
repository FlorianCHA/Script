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

#include "egglib.hpp"
#include "Structure.hpp"
#include "DataHolder.hpp"
#include <cstdlib>
#include <new>

namespace egglib {

    void StructureHolder::init() {
        _clusters = NULL;
        _pops = NULL;
        _indivs_i = NULL;
        _indivs_o = NULL;
        _num_clust_c = 0;
        _num_pop_c = 0;
        _num_indiv_i_c = 0;
        _num_indiv_o_c = 0;
        _num_filter_c = 0;
        _filter = NULL;
        reset();
    }

    //~~~~~~~~~//

    void StructureHolder::reset() {
        _ns = 0;
        _no = 0;
        _required_ns = 0;
        _required_no = 0;
        _ploidy = UNKNOWN;
        _num_clust = 0;
        _num_pop = 0;
        _num_indiv_i = 0;
        _num_indiv_o = 0;
        _num_filter = 0;
        _last_indiv = UNKNOWN;
    }

    //~~~~~~~~~//

    void StructureHolder::free() {
        if (_pops) ::free(_pops);
        if (_indivs_i) ::free(_indivs_i);
        if (_indivs_o) {
            for (unsigned int i=0; i<_num_indiv_o_c; i++) if (_indivs_o[i]) delete(_indivs_o[i]);
            ::free(_indivs_o);
        }
        if (_clusters) {
            for (unsigned int i=0; i<_num_clust_c; i++) {
                if (_clusters[i]) delete _clusters[i];
            }
            ::free(_clusters);
        }
        if (_filter) ::free(_filter);
    }

    //~~~~~~~~~//

    StructureHolder::StructureHolder() {
        init();
    }

    //~~~~~~~~~//

    StructureHolder::~StructureHolder() {
        free();
    }

    //~~~~~~~~~//

    void StructureHolder::reset_filter() {
        _num_filter = 0;
    }

    //~~~~~~~~~//

    void StructureHolder::copy(const StructureHolder& src) {
        reset();
        _ploidy = src._ploidy;
        const StructureCluster * src_clu;
        const StructurePopulation * src_pop;
        const StructureIndiv * src_idv;
        StructureCluster * dst_clu;
        StructurePopulation * dst_pop;
        StructureIndiv * dst_idv;
        for (unsigned int i=0; i<src.num_clust(); i++) {
            src_clu = & src.get_cluster(i);
            dst_clu = add_cluster(src_clu->get_label());
            for (unsigned int j=0; j<src_clu->num_pop(); j++) {
                src_pop = & (src_clu->get_population(j));
                dst_pop = add_population(src_pop->get_label(), dst_clu);
                for (unsigned int k=0; k<src_pop->num_indiv(); k++) {
                    src_idv = & (src_pop->get_indiv(k));
                    dst_idv = add_individual_ingroup(src_idv->get_label(), dst_clu, dst_pop);
                    for (unsigned int m=0; m<_ploidy; m++) {
                        add_sample_ingroup(src_idv->get_sample(m), dst_clu, dst_pop, dst_idv);
                    }
                }
            }
        }
        for (unsigned int i=0; i<src.num_indiv_outgroup(); i++) {
            src_idv = & src.get_indiv_outgroup(i);
            dst_idv = add_individual_outgroup(src_idv->get_label());
            for (unsigned int j=0; j<_ploidy; j++) {
                add_sample_outgroup(src_idv->get_sample(j), dst_idv);
            }
        }
    }

    //~~~~~~~~~//

    void StructureHolder::add_pop_filter(unsigned int lbl) {
        reserve_filter(_num_filter+1);
        _num_filter++;
        _filter[_num_filter-1] = lbl;
    }

    //~~~~~~~~~//

    void StructureHolder::reserve_filter(unsigned int howmany) {
        if (howmany > _num_filter_c) {
            _filter = (unsigned int *) realloc(_filter, howmany * sizeof(unsigned int));
            if (!_filter) throw EGGMEM;
            _num_filter_c = howmany;
        }
        }

    //~~~~~~~~~//

    unsigned int StructureHolder::get_ploidy() const {
        return _ploidy;
    }

    //~~~~~~~~~//

    void StructureHolder::get_structure(DataHolder& data, unsigned int lvl_clust,
        unsigned int lvl_pop, unsigned int lvl_indiv, unsigned int ploidy, bool skip_outgroup) {

        unsigned int lbl_clust;
        unsigned int lbl_pop;
        unsigned int lbl_indiv;

        // process all ingroup samples
        for (unsigned int idx_sam=0; idx_sam<data.get_nsam_i(); idx_sam++) {

            // get labels
            if (lvl_clust == UNKNOWN) lbl_clust = 0;
            else lbl_clust = data.get_group_i(idx_sam, lvl_clust);

            if (lvl_pop == UNKNOWN) lbl_pop = lbl_clust;
            else lbl_pop = data.get_group_i(idx_sam, lvl_pop);

            if (lvl_indiv == UNKNOWN) lbl_indiv = idx_sam;
            else lbl_indiv = data.get_group_i(idx_sam, lvl_indiv);

            // check if filter applies
            if (_num_filter > 0) {
                unsigned int i;
                for (i=0; i<_num_filter; i++) {
                    if (lbl_pop == _filter[i]) break;
                }
                if (i == _num_filter) continue; // skip it
            }

            // process sample
            process_ingroup(idx_sam, lbl_clust, lbl_pop, lbl_indiv);
        }

        // process outgroup samples
        if (skip_outgroup == false) {
            for (unsigned int idx_sam=0; idx_sam<data.get_nsam_o(); idx_sam++) {
                if (lvl_indiv == UNKNOWN) lbl_indiv = data.get_nsam_i() + idx_sam;
                else lbl_indiv = data.get_group_o(idx_sam);
                process_outgroup(idx_sam, lbl_indiv);
            }
        }

        // check and set ploidy
        if (lvl_indiv == UNKNOWN) _ploidy = 1;
        else check_ploidy(ploidy);
    }

    //~~~~~~~~~//

    void StructureHolder::process_ingroup(unsigned int sam_idx, unsigned int lbl_clust, unsigned int lbl_pop, unsigned int lbl_indiv) {

        unsigned int idx_clust;
        unsigned int idx_pop;
        unsigned int idx_indiv;

        // find cluster #############
        for (idx_clust=0; idx_clust<_num_clust; idx_clust++) {
            if (lbl_clust == _clusters[idx_clust]->get_label()) break;
        }

        // if cluster not found, create it
        if (idx_clust == _num_clust) add_cluster(lbl_clust);

        // find pop #############
        for (idx_pop=0; idx_pop<_num_pop; idx_pop++) {
            if (lbl_pop == _pops[idx_pop]->get_label()) break;
        }

        // if pop found, check it is in right cluster (enforce hierarchical)
        if (idx_pop < _num_pop) {
            if (_pops[idx_pop]->get_cluster() != _clusters[idx_clust]) {
                throw EggNonHierarchicalStructure(0, lbl_pop);
            }
        }

        // if pop not found, create it
        else {
            add_population(lbl_pop, _clusters[idx_clust]);
        }

        // find indiv #############
        // check last individual first (usually samples from the same individuals are together)
        if (_last_indiv < _num_indiv_i && _indivs_i[_last_indiv]->get_label() == lbl_indiv) {
            idx_indiv = _last_indiv;
            _last_indiv = UNKNOWN; // forget it (in most cases we'll have a diploid)
        }

        else {
            // check all ingroup
            for (idx_indiv=0; idx_indiv<_num_indiv_i; idx_indiv++) {
                if (_indivs_i[idx_indiv]->get_label() == lbl_indiv) break;
            }

            // if not found, check that NOT in outgroup
            if (idx_indiv == _num_indiv_i) {
                for (unsigned int i=0; i<_num_indiv_o; i++) {
                    if (_indivs_o[i]->get_label() == lbl_indiv) {
                        throw EggNonHierarchicalStructure(lbl_indiv);
                    }
                }

                // create individual
                add_individual_ingroup(lbl_indiv, _clusters[idx_clust], _pops[idx_pop]);
            }
        }

        // add sample
        add_sample_ingroup(sam_idx, _clusters[idx_clust], _pops[idx_pop], _indivs_i[idx_indiv]);
    }

    //~~~~~~~~~//

    void StructureHolder::process_outgroup(unsigned int sam_idx, unsigned int lbl_indiv) {

        // find indiv
        unsigned int idx_indiv;
        for (idx_indiv=0; idx_indiv<_num_indiv_o; idx_indiv++) {
            if (_indivs_o[idx_indiv]->get_label() == lbl_indiv) {
                break;
            }
        }

        if (idx_indiv == _num_indiv_o) {

            // if not found, check that not in ingroup
            for (unsigned int i=0; i<_num_indiv_i; i++) {
                if (_indivs_i[i]->get_label() == lbl_indiv) {
                    throw EggNonHierarchicalStructure(lbl_indiv);
                }
            }

            // create individual
            add_individual_outgroup(lbl_indiv);
        }

        // add sample
        add_sample_outgroup(sam_idx, _indivs_o[idx_indiv]);
    }

    //~~~~~~~~~//

    void StructureHolder::add_sample_ingroup(unsigned int sam_idx,
                                                StructureCluster * cluster,
                                                StructurePopulation * population,
                                                StructureIndiv * indiv) {
        _ns++;
        if (sam_idx+1 > _required_ns) _required_ns = sam_idx + 1;
        cluster->add_sample();
        population->add_sample();
        indiv->add_sample(sam_idx);
    }

    //~~~~~~~~~//

    void StructureHolder::add_sample_outgroup(unsigned int sam_idx,
                                                StructureIndiv * indiv) {
        _no++;
        if (sam_idx+1 > _required_no) _required_no = sam_idx + 1;
        indiv->add_sample(sam_idx);
    }

    //~~~~~~~~~//

    StructureCluster * StructureHolder::add_cluster(unsigned int label) {
        _num_clust++;
        if (_num_clust > _num_clust_c) {
            _clusters = (StructureCluster **) realloc(_clusters, _num_clust * sizeof(StructureCluster *));
            if (!_clusters) throw EGGMEM;
            _clusters[_num_clust-1] = new(std::nothrow) StructureCluster(this, label);
            if (!_clusters[_num_clust-1]) throw EGGMEM;
            _num_clust_c = _num_clust;
        }
        else {
            _clusters[_num_clust-1]->reset(this, label);
        }
        return _clusters[_num_clust-1];
    }

    //~~~~~~~~~//

    StructurePopulation * StructureHolder::add_population(unsigned int label,
                                                    StructureCluster * cluster) {
        _num_pop++;
        if (_num_pop > _num_pop_c) {
            _pops = (StructurePopulation **) realloc(_pops, _num_pop * sizeof(StructurePopulation *));
            if (!_pops) throw EGGMEM;
            _num_pop_c = _num_pop;
        }
        _pops[_num_pop-1] = cluster->add_pop(label);
        return _pops[_num_pop-1];
    }

    //~~~~~~~~~//

    StructureIndiv * StructureHolder::add_individual_ingroup(unsigned int label,
                    StructureCluster * cluster, StructurePopulation * population) {
        _num_indiv_i++;
        if (_num_indiv_i > _num_indiv_i_c) {
            _indivs_i = (StructureIndiv **) realloc(_indivs_i, _num_indiv_i * sizeof(StructureIndiv *));
            if (!_indivs_i) throw EGGMEM;
            _num_indiv_i_c = _num_indiv_i;
        }
        _indivs_i[_num_indiv_i-1] = cluster->add_indiv(population, label);
        return _indivs_i[_num_indiv_i-1];
    }

    //~~~~~~~~~//

    StructureIndiv * StructureHolder::add_individual_outgroup(unsigned int label) {
        _num_indiv_o++;
        if (_num_indiv_o > _num_indiv_o_c) {
            _indivs_o = (StructureIndiv **) realloc(_indivs_o, _num_indiv_o * sizeof(StructureIndiv *));
            if (!_indivs_o) throw EGGMEM;
            _indivs_o[_num_indiv_o-1] = new(std::nothrow) StructureIndiv(this, NULL, NULL, label);
            if (!_indivs_o[_num_indiv_o-1]) throw EGGMEM;
            _num_indiv_o_c = _num_indiv_o;
        }
        else {
            _indivs_o[_num_indiv_o-1]->reset(this, NULL, NULL, label);
        }
        return _indivs_o[_num_indiv_o-1];
    }

    //~~~~~~~~~//

    void StructureHolder::check_ploidy(unsigned int value) {
        _ploidy = value;
        if (_ploidy == UNKNOWN) {
            if (_num_indiv_i > 0) _ploidy = _indivs_i[0]->num_samples();
            else if (_num_indiv_o > 0) _ploidy = _indivs_o[0]->num_samples();
        }
        for (unsigned int i=0; i<_num_indiv_i; i++) {
            if (_indivs_i[i]->num_samples() != _ploidy) throw EggPloidyError();
        }
        for (unsigned int i=0; i<_num_indiv_o; i++) {
            if (_indivs_o[i]->num_samples() != _ploidy) throw EggPloidyError();
        }
    }

    //~~~~~~~~~//

    unsigned int StructureHolder::num_clust() const {
        return _num_clust;
    }

    //~~~~~~~~~//

    unsigned int StructureHolder::num_pop() const {
        return _num_pop;
    }

    //~~~~~~~~~//

    unsigned int StructureHolder::num_indiv_ingroup() const {
        return _num_indiv_i;
    }

    //~~~~~~~~~//

    unsigned int StructureHolder::num_indiv_outgroup() const {
        return _num_indiv_o;
    }

    //~~~~~~~~~//

    const StructureCluster& StructureHolder::get_cluster(unsigned int idx) const {
        return * _clusters[idx];
    }

    //~~~~~~~~~//

    const StructurePopulation& StructureHolder::get_population(unsigned int idx) const {
        return * _pops[idx];
    }

    //~~~~~~~~~//

    const StructureIndiv& StructureHolder::get_indiv_ingroup(unsigned int idx) const {
        return * _indivs_i[idx];
    }

    //~~~~~~~~~//

    const StructureIndiv& StructureHolder::get_indiv_outgroup(unsigned int idx) const {
        return * _indivs_o[idx];
    }

    //~~~~~~~~~//

    unsigned int StructureHolder::get_ns() const {
        return _ns;
    }

    //~~~~~~~~~//

    unsigned int StructureHolder::get_no() const {
        return _no;
    }

    //~~~~~~~~~//

    unsigned int StructureHolder::get_ns_req() const {
        return _required_ns;
    }

    //~~~~~~~~~//

    unsigned int StructureHolder::get_no_req() const {
        return _required_no;
    }

    //~~~~~~~~~//

    unsigned int StructureHolder::get_pop_index(unsigned int idx) const {
        for (unsigned int i=0; i<_num_pop; i++) {
            for (unsigned int j=0; j<_pops[i]->num_indiv(); j++) {
                for (unsigned int k=0; k<_ploidy; k++) {
                    if (idx == _pops[i]->get_indiv(j).get_sample(k)) return i;
                }
            }
        }
        return MISSING;
    }

    //~~~~~~~~~//     //~~~~~~~~~//     //~~~~~~~~~//     //~~~~~~~~~//

    void StructureCluster::init() {
        _num_pop_c = 0;
        _num_indiv_c = 0;
        _pops = NULL;
        _indivs = NULL;
    }

    //~~~~~~~~~//

    void StructureCluster::reset(StructureHolder * parent, unsigned int label) {
        _parent = parent;
        _label = label;
        _ns = 0;
        _num_pop = 0;
        _num_indiv = 0;
    }

    //~~~~~~~~~//

    void StructureCluster::free() {
        if (_pops) {
            for (unsigned int i=0; i<_num_pop_c; i++) {
                if (_pops[i]) delete _pops[i];
            }
            ::free(_pops);
        }
        if (_indivs) ::free(_indivs);
    }

    //~~~~~~~~~//

    StructureCluster::StructureCluster(StructureHolder * _parent, unsigned int label) {
        init();
        reset(_parent, label);
    }

    //~~~~~~~~~//

    StructureCluster::~StructureCluster() {
        free();
    }

    //~~~~~~~~~//

    StructurePopulation * StructureCluster::add_pop(unsigned int label) {
        _num_pop++;
        if (_num_pop > _num_pop_c) {
            _pops = (StructurePopulation **) realloc(_pops, _num_pop * sizeof(StructurePopulation *));
            if (!_pops) throw EGGMEM;
            _pops[_num_pop-1] = new(std::nothrow) StructurePopulation(_parent, this, label);
            if (!_pops[_num_pop-1]) throw EGGMEM;
            _num_pop_c = _num_pop;
        }
        else {
            _pops[_num_pop-1]->reset(_parent, this, label);
        }
        return _pops[_num_pop-1];
    }

    //~~~~~~~~~//

    StructureIndiv * StructureCluster::add_indiv(StructurePopulation * pop, unsigned int label) {
        _num_indiv++;
        if (_num_indiv > _num_indiv_c) {
            _indivs = (StructureIndiv **) realloc(_indivs, _num_indiv * sizeof(StructureIndiv *));
            if (!_indivs) throw EGGMEM;
            _num_indiv_c = _num_indiv;
        }
        _indivs[_num_indiv-1] = pop->add_indiv(label);
        return _indivs[_num_indiv-1];
    }

    //~~~~~~~~~//

    void StructureCluster::add_sample() {
        _ns++;
    }

    //~~~~~~~~~//

    unsigned int StructureCluster::num_pop() const {
        return _num_pop;
    }

    //~~~~~~~~~//

    unsigned int StructureCluster::num_indiv() const {
        return _num_indiv;
    }

    //~~~~~~~~~//

    const StructurePopulation& StructureCluster::get_population(unsigned int idx) const {
        return * _pops[idx];
    }

    //~~~~~~~~~//

    const StructureIndiv& StructureCluster::get_indiv(unsigned int idx) const {
        return * _indivs[idx];
    }

    //~~~~~~~~~//

    unsigned int StructureCluster::get_label() const {
        return _label;
    }

    //~~~~~~~~~//    //~~~~~~~~~//    //~~~~~~~~~//    //~~~~~~~~~//

    void StructurePopulation::init() {
        _num_indiv_c = 0;
        _indivs = NULL;
    }

    //~~~~~~~~~//

    void StructurePopulation::reset(StructureHolder * holder,
                                    StructureCluster * cluster,
                                    unsigned int label) {
        _holder = holder;
        _cluster = cluster;
        _label = label;
        _ns = 0;
        _num_indiv = 0;
    }

    //~~~~~~~~~//

    void StructurePopulation::free() {
        if (_indivs) {
            for (unsigned int i=0; i<_num_indiv_c; i++) {
                if (_indivs[i]) delete _indivs[i];
            }
            ::free(_indivs);
        }
    }

    //~~~~~~~~~//

    StructurePopulation::StructurePopulation(StructureHolder * holder,
                                             StructureCluster * cluster,
                                             unsigned int label) {
        init();
        reset(holder, cluster, label);
    }

    //~~~~~~~~~//

    StructurePopulation::~StructurePopulation() {
        free();
    }

    //~~~~~~~~~//

    StructureIndiv * StructurePopulation::add_indiv(unsigned int label) {
        _num_indiv++;
        if (_num_indiv > _num_indiv_c) {
            _indivs = (StructureIndiv **) realloc(_indivs, _num_indiv * sizeof(StructureIndiv *));
            if (!_indivs) throw EGGMEM;
            _indivs[_num_indiv-1] = new(std::nothrow) StructureIndiv(_holder, _cluster, this, label);
            if (!_indivs[_num_indiv-1]) throw EGGMEM;
            _num_indiv_c = _num_indiv;
        }
        else {
            _indivs[_num_indiv-1]->reset(_holder, _cluster, this, label);
        }
        return _indivs[_num_indiv-1];
    }

    //~~~~~~~~~//

    void StructurePopulation::add_sample() {
        _ns++;
    }

    //~~~~~~~~~//

    unsigned int StructurePopulation::num_indiv() const {
        return _num_indiv;
    }

    //~~~~~~~~~//

    const StructureIndiv& StructurePopulation::get_indiv(unsigned int idx) const {
        return * _indivs[idx];
    }

    //~~~~~~~~~//

    StructureCluster * StructurePopulation::get_cluster() {
        return _cluster;
    }

    //~~~~~~~~~//

    unsigned int StructurePopulation::get_label() const {
        return _label;
    }

    //~~~~~~~~~//    //~~~~~~~~~//    //~~~~~~~~~//    //~~~~~~~~~//

    void StructureIndiv::init() {
        _num_sam_c = 0;
        _samples = NULL;
    }

    //~~~~~~~~~//

    void StructureIndiv::reset(StructureHolder * holder,
                           StructureCluster * cluster,
                           StructurePopulation * population,
                           unsigned int label) {
        _holder = holder;
        _cluster = cluster;
        _population = population;
        _label = label;
        _num_sam = 0;
    }

    //~~~~~~~~~//

    void StructureIndiv::free() {
        if (_samples) ::free(_samples);
    }

    //~~~~~~~~~//

    StructureIndiv::StructureIndiv(StructureHolder * holder,
                           StructureCluster * cluster,
                           StructurePopulation * population,
                           unsigned int label) {
        init();
        reset(holder, cluster, population, label);
    }

    //~~~~~~~~~//

    StructureIndiv::~StructureIndiv() {
        free();
    }

    //~~~~~~~~~//

    void StructureIndiv::add_sample(unsigned int index) {
        _num_sam++;
        if (_num_sam > _num_sam_c) {
            _samples = (unsigned int *) realloc(_samples, _num_sam * sizeof(unsigned int));
            if (!_samples) throw EGGMEM;
            _num_sam_c = _num_sam;
        }
        _samples[_num_sam-1] = index;
    }

    //~~~~~~~~~//

    unsigned int StructureIndiv::num_samples() const {
        return _num_sam;
    }

    //~~~~~~~~~//

    unsigned int StructureIndiv::get_sample(unsigned int idx) const {
        return _samples[idx];
    }

    //~~~~~~~~~//

    StructureCluster * StructureIndiv::get_cluster() {
        return _cluster;
    }

    //~~~~~~~~~//

    StructurePopulation * StructureIndiv::get_population() {
        return _population;
    }

    //~~~~~~~~~//

    unsigned int StructureIndiv::get_label() const {
        return _label;
    }
}
