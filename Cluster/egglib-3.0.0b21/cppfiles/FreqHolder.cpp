/*
    Copyright 2016 Stéphane De Mita, Mathieu Siol

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

#include <cstdlib>
#include <cstring>
#include "egglib.hpp"
#include "FreqHolder.hpp"
#include "SiteHolder.hpp"
#include "Structure.hpp"
#include "VCF.hpp"

namespace egglib {

    FreqSet::FreqSet() {
        init();
    }

    ////////////

    FreqSet::~FreqSet() {
        free();
    }

    ////////////

    void FreqSet::init() {
        _nall_c = 0;
        _ngen_c = 0;
        _frq_all = NULL;
        _frq_het = NULL;
        _frq_gen = NULL;
        _gen_het = NULL;
        setup();
        reset(0);
    }

    ////////////

    void FreqSet::free() {
        if (_frq_all) ::free(_frq_all);
        if (_frq_het) ::free(_frq_het);
        if (_frq_gen) ::free(_frq_gen);
        if (_gen_het) ::free(_gen_het);
    }

    ////////////

    void FreqSet::setup() { 
        _nall = 0;
        _nall_eff = 0;
        _ngen = 0;
        _ngen_eff = 0;
        _nsam = 0;
        _nind = 0;
        _nhet = 0;
    }

    ////////////

    void FreqSet::reset(unsigned int na) {
        _nall = na;
        if (_nall > _nall_c) {
            _frq_all = (unsigned int *) realloc(_frq_all, _nall * sizeof(unsigned int));
            if (!_frq_all) throw EGGMEM;
            _frq_het = (unsigned int *) realloc(_frq_het, _nall * sizeof(unsigned int));
            if (!_frq_het) throw EGGMEM;
            _nall_c = _nall;
        }
        for (unsigned int i=0; i<_nall; i++) {
            _frq_all[i] = 0;
            _frq_het[i] = 0;
        }
        _ngen = 0;
        _ngen = 0;
        _nsam = 0;
        _nhet = 0;
        _nall_eff = 0;
        _ngen_eff = 0;
    }

    ////////////

    void FreqSet::add_genotypes(unsigned int num) {
        _ngen += num;
        if (_ngen > _ngen_c) {
            _frq_gen = (unsigned int *) realloc(_frq_gen, _ngen * sizeof(unsigned int));
            if (!_frq_gen) throw EGGMEM;
            _gen_het = (bool *) realloc(_gen_het, _ngen * sizeof(bool));
            if (!_gen_het) throw EGGMEM;
            _ngen = _ngen;
        }
        for (unsigned int i=0; i<num; i++) {
            _frq_gen[_ngen - 1 - i] = 0;
            _gen_het[_ngen - 1 - i] = false;
        }
    }

    ////////////

    void FreqSet::incr_allele(unsigned int all_idx, unsigned int num) {
        _nsam += num;
        if (_frq_all[all_idx] == 0 && num > 0) _nall_eff++;
        _frq_all[all_idx] += num;
    }

    ////////////

    void FreqSet::incr_genotype(unsigned int gen_idx, unsigned int num) {
        _nind += num;
        if (_frq_gen[gen_idx] == 0 && num > 0) _ngen_eff++;
        _frq_gen[gen_idx] += num;
    }

    ////////////

    void FreqSet::tell_het(unsigned int i, unsigned int a) {
        _frq_het[a] += _frq_gen[i];
        if (!_gen_het[i]) {
            _gen_het[i] = true;
            _nhet += _frq_gen[i];
        }
    }

    ////////////

    unsigned int FreqSet::num_alleles() const {
        return _nall;
    }

    ////////////

    unsigned int FreqSet::num_alleles_eff() const {
        return _nall_eff;
    }

    ////////////

    unsigned int FreqSet::num_genotypes() const {
        return _ngen;
    }

    ////////////

    unsigned int FreqSet::num_genotypes_eff() const {
        return _ngen_eff;
    }

    ////////////

    unsigned int FreqSet::nseff() const {
        return _nsam;
    }

    ////////////

    unsigned int FreqSet::nieff() const {
        return _nind;
    }

    ////////////

    unsigned int FreqSet::frq_all(unsigned int i) const {
        return _frq_all[i];
    }

    ////////////

    unsigned int FreqSet::frq_het(unsigned int all) const {
        return _frq_het[all];
    }

    ////////////

    unsigned int FreqSet::tot_het() const {
        return _nhet;
    }

    ////////////

    unsigned int FreqSet::frq_gen(unsigned int i) const {
        return _frq_gen[i];
    }

    ////////////////////////////////////////////////////////////////////

    FreqHolder::FreqHolder() {
        _npop_c = 0;
        _nclu_c = 0;
        _frq_clu = NULL;
        _frq_pop = NULL;
        _clu_idx = NULL;
        _rel_pop_idx = NULL;
        _pop_ns = NULL;
        _nall_c = 0;
        _ngen_c = 0;
        _gen_c2 = NULL;
        _genotypes = NULL;
        _matched_c = 0;
        _matched = NULL;
        _alleles = NULL;
        _gen_het = NULL;
        _structure = NULL;

        _npop = 0;
        _nclu = 0;
        _outgroup_ns = 0;
        _nall = 0;
        _pl = 0;
        _ngen = 0;
    }

    ////////////

    FreqHolder::~FreqHolder() {
        if (_frq_clu) {
            for (unsigned int i=0; i<_nclu_c; i++) {
                if (_frq_clu[i]) delete _frq_clu[i];
            }
            free(_frq_clu);
        }
        if (_frq_pop) {
            for (unsigned int i=0; i<_npop_c; i++) {
                if (_frq_pop[i]) delete _frq_pop[i];
            }
            free(_frq_pop);
        }
        if (_genotypes) {
            for (unsigned int i=0; i<_ngen_c; i++) {
                if (_genotypes[i]) free(_genotypes[i]);
            }
            free(_genotypes);
        }
        if (_gen_c2) free(_gen_c2);
        if (_matched) free(_matched);
        if (_alleles) free(_alleles);
        if (_gen_het) free(_gen_het);
        if (_clu_idx) free(_clu_idx);
        if (_rel_pop_idx) free(_rel_pop_idx);
        if (_pop_ns) free(_pop_ns);
    }

    ////////////

    void FreqHolder::_set_frq(unsigned int nc, unsigned int np) {
        _nclu = nc;
        if (_nclu > _nclu_c) {
            _frq_clu = (FreqSet **) realloc(_frq_clu, _nclu * sizeof(FreqSet *));
            if (!_frq_clu) throw EGGMEM;
            for (unsigned int i=_nclu_c; i<_nclu; i++) {
                _frq_clu[i] = new(std::nothrow) FreqSet;
                if (!_frq_clu[i]) throw EGGMEM;
            }
            _nclu_c = _nclu;
        }

        _npop = np;
        if (_npop > _npop_c) {
            _frq_pop = (FreqSet **) realloc(_frq_pop, _npop * sizeof(FreqSet *));
            if (!_frq_pop) throw EGGMEM;
            for (unsigned int i=_npop_c; i<_npop; i++) {
                _frq_pop[i] = new(std::nothrow) FreqSet;
                if (!_frq_pop[i]) throw EGGMEM;
            }
            _clu_idx = (unsigned int *) realloc(_clu_idx, _npop * sizeof(unsigned int));
            if (!_clu_idx) throw EGGMEM;
            _rel_pop_idx = (unsigned int *) realloc(_rel_pop_idx, _npop * sizeof(unsigned int));
            if (!_rel_pop_idx) throw EGGMEM;
            _pop_ns = (unsigned int *) realloc(_pop_ns, _npop * sizeof(unsigned int));
            if (!_pop_ns) throw EGGMEM;
            _npop_c = _npop;
        }
    }

    ////////////

    void FreqHolder::_set_nall(unsigned int na) {
        _nall = na;
        if (_nall > _nall_c) {
            _alleles = (int *) realloc(_alleles, _nall * sizeof(int));
            if (!_alleles) throw EGGMEM;
            _nall_c = _nall;
        }

        if (_nall > _matched_c) {
            _matched = (bool *) realloc(_matched, _nall * sizeof(bool));
            if (!_matched) throw EGGMEM;
            _matched_c = na;
        }
        _frq_ing.reset(na);
        _frq_otg.reset(na);
        for (unsigned int i=0; i<_nclu; i++) _frq_clu[i]->reset(na);
        for (unsigned int i=0; i<_npop; i++) _frq_pop[i]->reset(na);
    }

    ////////////

    void FreqHolder::_add_genotypes(unsigned int num) {
        _ngen += num;
        if (_ngen > _ngen_c) {
            _genotypes = (unsigned int **) realloc(_genotypes, _ngen * sizeof(unsigned int *));
            if (!_genotypes) throw EGGMEM;
            _gen_c2 = (unsigned int *) realloc(_gen_c2, _ngen * sizeof(unsigned int));
            if (!_gen_c2) throw EGGMEM;
            _gen_het = (bool *) realloc(_gen_het, _ngen * sizeof(bool));
            if (!_gen_het) throw EGGMEM;
            for (unsigned int i=_ngen_c; i<_ngen; i++) {
                _gen_c2[i] = 0;
                _genotypes[i] = NULL;
                _gen_het[i] = false;
            }
            _ngen_c = _ngen;
        }
        for (unsigned int i=0; i<num; i++) {
            if (_pl > _gen_c2[_ngen-1-i]) {
                _genotypes[_ngen-1-i] = (unsigned int *) realloc(_genotypes[_ngen-1-i], _pl * sizeof(unsigned int));
                if (!_genotypes[_ngen-1-i]) throw EGGMEM;
                _gen_c2[_ngen-1-i] = _pl;
            }
        }
        _frq_ing.add_genotypes(num);
        _frq_otg.add_genotypes(num);
        for (unsigned int i=0; i<_nclu; i++) _frq_clu[i]->add_genotypes(num);
        for (unsigned int i=0; i<_npop; i++) _frq_pop[i]->add_genotypes(num);
    }

    ////////////

    unsigned int FreqHolder::_find_genotype(const unsigned int * alleles) {
        for (unsigned int i=0; i<_pl; i++) {
            if (alleles[i] == MISSING) return MISSING;
        }
        unsigned int idx, all_idx, all_idx2;
        for (idx=0; idx<_ngen; idx++) {
            for (all_idx=0; all_idx<_pl; all_idx++) _matched[all_idx] = false;
            for (all_idx=0; all_idx<_pl; all_idx++) {
                for (all_idx2=0; all_idx2<_pl; all_idx2++) {
                    if (_matched[all_idx2]==false && _genotypes[idx][all_idx2] == alleles[all_idx]) {
                        _matched[all_idx2] = true;
                        break; // found this allele
                    }
                }
                if (all_idx2 == _pl) break; // did not find a match for alleles[all_idx]
            }
            if (all_idx == _pl) break; // all alleles are matching
        }
        if (idx == _ngen) { // no genotype is matching
            _add_genotypes(1);
            for (all_idx=0; all_idx<_pl; all_idx++) {
                _genotypes[idx][all_idx] = alleles[all_idx];
            }
        }
        return idx;
    }

    ////////////

    void FreqHolder::setup_structure(const StructureHolder * structure, unsigned int ploidy) {
        _pl = ploidy;
        _ngen = 0;
        _nall = 0;
        _structure = structure;

        // setup global FreqSet objects
        _frq_ing.setup();
        _frq_otg.setup();
        _outgroup_ns = structure->num_indiv_outgroup();

        // create needed FreqSet objects
        _set_frq(structure->num_clust(), structure->num_pop());

        // setup FreqSet objects
        unsigned int cur = 0;
        for (unsigned int i=0; i<_nclu; i++) {
            _frq_clu[i]->setup();
            for (unsigned int j=0; j<structure->get_cluster(i).num_pop(); j++) {
                _frq_pop[cur]->setup();
                _clu_idx[cur] = i;
                _rel_pop_idx[cur] = j;
                _pop_ns[cur] = structure->get_cluster(i).get_population(j).num_indiv();
                cur++;
            }
        }
    }

    ////////////

    void FreqHolder::setup_raw(unsigned int nc, unsigned int np, unsigned int no, unsigned int ploidy) {
        _pl = ploidy;
        _ngen = 0;
        _nall = 0;
        _structure = NULL;
        _outgroup_ns = no;

        _set_frq(nc, np);
        _frq_ing.setup();
        _frq_otg.setup();
        for (unsigned int i=0; i<nc; i++) _frq_clu[i]->setup();
    }

    ////////////

    void FreqHolder::setup_pop(unsigned int idx, unsigned int cluster, unsigned int rel_idx, unsigned int ns) {
        _frq_pop[idx]->setup();
        _clu_idx[idx] = cluster;
        _rel_pop_idx[idx] = rel_idx;
        _pop_ns[idx] = ns;
    }

    ////////////

    void FreqHolder::set_nall(unsigned int na, unsigned int ng) {
        _set_nall(na);
        _ngen = 0;
        _add_genotypes(ng);
    }

    ////////////

    void FreqHolder::process_site(const SiteHolder& site) {

        // set variables
        _ngen = 0;
        _set_nall(site.get_nall());
        for (unsigned int i=0; i<_nall; i++) _alleles[i] = site.get_allele(i);

        // set ingroup freq
        unsigned int gen_idx;
        unsigned int indiv_idx = 0;
        for (unsigned int p=0; p<_npop; p++) {
            for (unsigned int i=0; i<_pop_ns[p]; i++) {
                if (_structure != NULL) indiv_idx = _structure->get_cluster(_clu_idx[p]).get_population(_rel_pop_idx[p]).get_indiv(i).get_sample(0);
                if (_pl > 1) {
                    gen_idx = _find_genotype(site.get_pi(indiv_idx));
                    if (gen_idx != MISSING) {
                        _frq_ing.incr_genotype(gen_idx, 1);
                        _frq_clu[_clu_idx[p]]->incr_genotype(gen_idx, 1);
                        _frq_pop[p]->incr_genotype(gen_idx, 1);
                    }
/*                    if (geno_table != NULL) {
                        if (gen_idx == _ngen-1) {
                            geno_table->add_genotype();
                            unsigned int j;
                            for (j=1; j<_pl; j++) if (site.get_i(indiv_idx, j) != site.get_i(indiv_idx, 0)) {
                                geno_table->set_heterozygote(gen_idx, true);
                                break;
                            }
                            if (j == _pl) geno_table->set_heterozygote(gen_idx, false);
                        }
                        geno_table->set_genotype_i(p, i, gen_idx);
                    }*/
                }
                for (unsigned int j=0; j<_pl; j++) {
                    unsigned int all = site.get_i(indiv_idx, j);
                    if (all != MISSING) {
                        _frq_ing.incr_allele(all, 1);
                        _frq_clu[_clu_idx[p]]->incr_allele(all, 1);
                        _frq_pop[p]->incr_allele(all, 1);
                    }
                }
                if (_structure == NULL) indiv_idx++;
            }
        }

        // set outgroup freq
        if (_structure == NULL) indiv_idx = 0;
        for (unsigned int i=0; i<_outgroup_ns; i++) {
            if (_structure != NULL) indiv_idx = _structure->get_indiv_outgroup(i).get_sample(0);
            if (_pl > 1) {
                gen_idx = _find_genotype(site.get_po(indiv_idx));
                if (gen_idx != MISSING) {
                    _frq_otg.incr_genotype(gen_idx, 1);
                }
/*                if (geno_table != NULL) {
                    if (gen_idx == _ngen-1) {
                        geno_table->add_genotype();
                        unsigned int j;
                        for (j=1; j<_pl; j++) if (site.get_o(indiv_idx, j) != site.get_o(indiv_idx, 0)) {
                            geno_table->set_heterozygote(gen_idx, true);
                            break;
                        }
                        if (j == _pl) geno_table->set_heterozygote(gen_idx, false);
                    }
                    geno_table->set_genotype_o(i, gen_idx);
                }*/
            }
            for (unsigned int j=0; j<_pl; j++) {
                unsigned int all = site.get_o(indiv_idx, j);
                if (all != MISSING) {
                    _frq_otg.incr_allele(all, 1);
                }
            }
            if (_structure == NULL) indiv_idx++;
        }

        // process heterozygote genotypes
        if (_pl > 1) {
            for (unsigned int i=0; i<_ngen; i++) {
                unsigned int j;
                for (j=1; j<_pl; j++) {
                    if (_genotypes[i][j] != _genotypes[i][0]) break;
                }
                // if heterozygote
                if (j < _pl) {
                    _gen_het[i] = true;
                    for (j=0; j<_nall; j++) _matched[j] = false;
                    for (j=0; j<_pl; j++) _matched[_genotypes[i][j]] = true;
                    for (j=0; j<_nall; j++) {
                        if (_matched[j]) {
                            _frq_ing.tell_het(i, j);
                            _frq_otg.tell_het(i, j);
                            for (unsigned int k=0; k<_nclu; k++) _frq_clu[k]->tell_het(i, j);
                            for (unsigned int k=0; k<_npop; k++) _frq_pop[k]->tell_het(i, j);
                        }
                    }
                }
            }
        }
    }

    ////////////

    void FreqHolder::process_vcf(const VcfParser& vcf) {
        unsigned int AN = vcf.AN();
        if (AN == UNKNOWN) throw EggArgumentValueError("cannot import VCF data: AN is missing");
        unsigned int acc = 0;
        setup_raw(1, 1, 0, 1);
        setup_pop(0, 0, 0, vcf.num_AC());
        _ngen = 0;
        _set_nall(vcf.num_AC()+1);

        if (strlen(vcf.reference()) == 1) _alleles[0] = static_cast<int>(vcf.reference()[0]);
        else _alleles[0] = 0;
        for (unsigned int i=1; i<_nall; i++) {
            if (vcf.alternate_type(i-1) == vcf::Default && strlen(vcf.alternate(i-1)) == 1) _alleles[i] = static_cast<int>(vcf.alternate(i-1)[0]);
            else _alleles[i] = i;
        }

        unsigned int AC;
        for (unsigned int i=1; i<_nall; i++) {
            AC = vcf.AC(i-1);
            if (AC == UNKNOWN) AC = 0;
            _frq_ing.incr_allele(i, AC);
            _frq_clu[0]->incr_allele(i, AC);
            _frq_pop[0]->incr_allele(i, AC);
            acc += AC;
        }
        if (acc > AN) throw EggRuntimeError("invalid VCF data: sum of AC fields is > AN");
        acc = AN - acc;
        _frq_ing.incr_allele(0, acc);
        _frq_clu[0]->incr_allele(0, acc);
        _frq_pop[0]->incr_allele(0, acc);
    }

    ////////////

    const FreqSet& FreqHolder::frq_ingroup() const {
        return _frq_ing;
    }

    ////////////

    const FreqSet& FreqHolder::frq_outgroup() const {
        return _frq_otg;
    }

    ////////////

    const FreqSet& FreqHolder::frq_cluster(unsigned int i) const {
        return * _frq_clu[i];
    }

    ////////////

    const FreqSet& FreqHolder::frq_population(unsigned int i) const {
        return * _frq_pop[i];
    }

    ////////////

    unsigned int FreqHolder::ploidy() const {
        return _pl;
    }

    ////////////

    unsigned int FreqHolder::num_alleles() const {
        return _nall;
    }

    ////////////

    unsigned int FreqHolder::num_genotypes() const {
        return _ngen;
    }

    ////////////

    const unsigned int * FreqHolder::genotype(unsigned int i) const {
        return _genotypes[i];
    }

    ////////////

    bool FreqHolder::genotype_het(unsigned int i) const {
        return _gen_het[i];
    }

    ////////////

    unsigned int FreqHolder::genotype_item(unsigned int i, unsigned int j) const {
        return _genotypes[i][j];
    }

    ////////////

    void FreqHolder::set_genotype_item(unsigned int i, unsigned int j, unsigned int a) {
        _genotypes[i][j] = a;
    }

    ////////////

    unsigned int FreqHolder::num_clusters() const {
        return _nclu;
    }

    ////////////

    unsigned int FreqHolder::num_populations() const {
        return _npop;
    }

    ////////////

    int FreqHolder::allele(unsigned int i) const {
        return _alleles[i];
    }

    ////////////

    unsigned int FreqHolder::cluster_index(unsigned int i) const {
        return _clu_idx[i];
    }

/* on the way to removal but kept here because not committed

    ////////////

    void GenotypesTable::_init() {
        _num_genotypes = 0;
        _sz_genotypes = 0;
        _num_pop = 0;
        _sz_pop = 0;
        _nind = NULL;
        _notg = 0;
        _sz_ind = NULL;
        _sz_otg = 0;
        _genotypes_i = NULL;
        _genotypes_o = NULL;
        _heter = NULL;
    }

    ////////////

    void GenotypesTable::_free() {
        if (_genotypes_i) {
            for (unsigned int i=0; i<_sz_pop; i++) if (_genotypes_i[i]) free(_genotypes_i[i]);
            free(_genotypes_i);
        }
        if (_genotypes_o) free(_genotypes_o);
        if (_nind) free(_nind);
        if (_sz_ind) free(_sz_ind);
        if (_heter) free(_heter);
    }

    ////////////

    GenotypesTable::GenotypesTable() {
        _init();
    }

    ////////////

    GenotypesTable::~GenotypesTable() {
        _free();
    }

    ////////////

    void GenotypesTable::reset() {
        _num_genotypes = 0;
    }

    ////////////

    void GenotypesTable::set_num_pop(unsigned int num_pop) {
        _num_pop = num_pop;
        if (_num_pop > _sz_pop) {
            _nind = (unsigned int *) realloc(_nind, _num_pop * sizeof(unsigned int));
            if (!_nind) throw EGGMEM;
            _sz_ind = (unsigned int *) realloc(_sz_ind, _num_pop * sizeof(unsigned int));
            if (!_sz_ind) throw EGGMEM;
            _genotypes_i = (unsigned int **) realloc(_genotypes_i, _num_pop * sizeof(unsigned int *));
            if (!_genotypes_i) throw EGGMEM;
            for (unsigned int i=_sz_pop; i<_num_pop; i++) {
                _sz_ind[i] = 0;
                _genotypes_i[i] = NULL;
            }
            _sz_pop = _num_pop;
        }
    }

    ////////////

    void GenotypesTable::set_num_idv(unsigned int pop, unsigned int num) {
        _nind[pop] = num;
        if (num > _sz_ind[pop]) {
            _genotypes_i[pop] = (unsigned int *) realloc(_genotypes_i[pop], num * sizeof(unsigned int));
            if (!_genotypes_i[pop]) throw EGGMEM;
            _sz_ind[pop] = num;
        }
    }

    ////////////

    void GenotypesTable::set_num_otg(unsigned int num) {
        _notg = num;
        if (num > _sz_otg) {
            _genotypes_o = (unsigned int *) realloc(_genotypes_o, num * sizeof(unsigned int));
            if (_genotypes_o) throw EGGMEM;
            _sz_otg = num;
        }
    }

    ////////////

    void GenotypesTable::add_genotype() {
        _num_genotypes++;
        if (_num_genotypes > _sz_genotypes) {
            _heter = (bool *) realloc(_heter, _num_genotypes * sizeof(bool));
            if (!_heter) throw EGGMEM;
        }
    }

    ////////////

    unsigned int GenotypesTable::num_genotypes() {
        return _num_genotypes;
    }

    ////////////

    void GenotypesTable::set_heterozygote(unsigned int genotype, bool heterozygote) {
        _heter[genotype] = heterozygote;
    }

    ////////////

    void GenotypesTable::set_genotype_i(unsigned int pop, unsigned int idv, unsigned int genotype) {
        _genotypes_i[pop][idv] = genotype;
    }

    ////////////

    void GenotypesTable::set_genotype_o(unsigned int idv, unsigned int genotype) {
        _genotypes_o[idv] = genotype;
    }

    ////////////

    unsigned int GenotypesTable::num_pop() const {
        return _num_pop;
    }

    ////////////

    unsigned int GenotypesTable::num_idv(unsigned int pop) const {
        return _nind[pop];
    }

    ////////////

    unsigned int GenotypesTable::num_otg() const {
        return _notg;
    }

    ////////////

    unsigned int GenotypesTable::genotype_i(unsigned int i, unsigned int j) const {
        return _genotypes_i[i][j];
    }

    ////////////

    unsigned int GenotypesTable::genotype_o(unsigned int i) const {
        return _genotypes_o[i];
    }

    ////////////

    bool GenotypesTable::genotype_heter(unsigned int i) const {
        return _heter[i];
    }

    ////////////

    bool GenotypesTable::ing_heter(unsigned int i, unsigned int j) {
        return _heter[_genotypes_i[i][j]];
    }

    ////////////

    bool GenotypesTable::otg_heter(unsigned int i) {
        return _heter[_genotypes_o[i]];
    }
*/

}
