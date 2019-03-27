/*
    Copyright 2012-2018 Stéphane De Mita, Mathieu Siol

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
#include <cstdlib>
#include <cstring>
#include "DataHolder.hpp"
#include "Structure.hpp"
#include "SiteHolder.hpp"
#include "Filter.hpp"
#include "VCF.hpp"

namespace egglib {

    SiteHolder::SiteHolder() {
        _init(1);
    }

    SiteHolder::SiteHolder(unsigned int ploidy) {
        _init(ploidy);
    }

    void SiteHolder::_init(unsigned int ploidy) {
        _nall_c = 0;
        _alleles = NULL;
        _ni_c = 0;
        _no_c = 0;
        _ing = NULL;
        _otg = NULL;
        _flag = false;
        reset(ploidy);
    }

    SiteHolder::~SiteHolder() {
        if (_alleles) ::free(_alleles);
        if (_ing) ::free(_ing);
        if (_otg) ::free(_otg);
    }

    void SiteHolder::reset(unsigned int ploidy) {
        _flag = false;
        _nall = 0;
        _nall_ing = 0;
        _ni = 0;
        _no = 0;
        _ploidy = ploidy;
        _missing = 0;
        _tot_missing = 0;
        _missing_ing = 0;
        _missing_otg = 0;
    }

    void SiteHolder::add_ing(unsigned int num) {
        _ni += num;
        if (_ni * _ploidy > _ni_c) {
            _ing = (unsigned int *) realloc(_ing, (_ni * _ploidy) * sizeof(unsigned int));
            if (!_ing) throw EGGMEM;
            _ni_c = _ni * _ploidy;
        }

    }

    void SiteHolder::add_otg(unsigned int num) {
        _no += num;
        if (_no * _ploidy > _no_c) {
            _otg = (unsigned int *) realloc(_otg, (_no * _ploidy) * sizeof(unsigned int));
            if (!_otg) throw EGGMEM;
            _no_c = _no * _ploidy;
        }

    }

    unsigned int SiteHolder::get_ning() const {
        return _ni;
    }

    unsigned int SiteHolder::get_nout() const {
        return _no;
    }

    unsigned int SiteHolder::get_ploidy() const {
        return _ploidy;
    }

    unsigned int SiteHolder::get_nall() const {
        return _nall;
    }

    unsigned int SiteHolder::get_nall_ing() const {
        return _nall_ing;
    }

    int SiteHolder::get_allele(unsigned int i) const {
        if (i == MISSING) return MISSINGDATA;
        return _alleles[i];
    }

    void SiteHolder::set_nall(unsigned int n, unsigned int ning) {
        _nall = n;
        _nall_ing = ning;
        if (_nall > _nall_c) {
            _alleles = (int *) realloc(_alleles, _nall * sizeof(int));
            if (!_alleles) throw EGGMEM;
            _nall_c = _nall;
        }
    }

    void SiteHolder::set_allele(unsigned int i, int a) const {
        _alleles[i] = a;
    }

    unsigned int SiteHolder::get_i(unsigned int idv, unsigned int chrom) const {
        return _ing[idv*_ploidy + chrom];
    }

    unsigned int SiteHolder::get_o(unsigned int idv, unsigned int chrom) const {
        return _otg[idv*_ploidy + chrom];
    }

    unsigned int SiteHolder::get_straight_i(unsigned int sam) const {
        return _ing[sam];
    }

    unsigned int SiteHolder::get_straight_o(unsigned int sam) const {
        return _otg[sam];
    }

    void SiteHolder::set_i(unsigned int idv, unsigned int chrom, unsigned int all) {
        if (all == MISSING) {
            _tot_missing++;
            _missing_ing++;
        }
        _ing[idv*_ploidy + chrom] = all;
    }

    void SiteHolder::set_o(unsigned int idv, unsigned int chrom, unsigned int all) {
        if (all == MISSING) {
            _tot_missing++;
            _missing_otg++;
        }
        _otg[idv*_ploidy + chrom] = all;
    }

    void SiteHolder::load_ing(unsigned int idv, unsigned int chrom, int allele) {
        set_i(idv, chrom, _process_allele(allele, true));
    }

    void SiteHolder::load_otg(unsigned int idv, unsigned int chrom, int allele) {
        set_o(idv, chrom, _process_allele(allele, false));
    }

    const unsigned int * SiteHolder::get_pi(unsigned int idv) const {
        return _ing + idv*_ploidy;
    }

    const unsigned int * SiteHolder::get_po(unsigned int idv) const {
        return _otg + idv*_ploidy;
    }

    unsigned int SiteHolder::get_missing() const {
        return _missing;
    }

    unsigned int SiteHolder::get_missing_ing() const {
        return _missing_ing;
    }

    unsigned int SiteHolder::get_missing_otg() const {
        return _missing_otg;
    }

    unsigned int SiteHolder::get_tot_missing() const {
        return _tot_missing;
    }

    bool SiteHolder::process_align(const DataHolder& data, unsigned int idx,
            const StructureHolder * struc, const Filter& filtr,
            unsigned int max_missing, bool consider_outgroup_missing) {

        if (data.get_is_matrix() == false) throw EggArgumentValueError("argument must be an alignment");
        const StructureCluster * clu;
        const StructurePopulation * pop;
        const StructureIndiv * idv;
        _data_pointer = & data;
        _filter_pointer = & filtr;
        _site_index = idx;
        _missing = 0;
        _consider_outgroup_missing = consider_outgroup_missing;

        unsigned int cur_i = _ni;
        unsigned int cur_o = _no;

        if (struc) {
            unsigned int cnt = 0;
            add_ing(struc->num_indiv_ingroup());
            add_otg(struc->num_indiv_outgroup());
            for (unsigned int i=0; i<struc->num_clust(); i++) {
                clu = & struc->get_cluster(i);
                for (unsigned int j=0; j<clu->num_pop(); j++) {
                    pop = & clu->get_population(j);
                    for (unsigned int k=0; k<pop->num_indiv(); k++) {
                        idv = & pop->get_indiv(k);
                        for (unsigned int m=0; m<_ploidy; m++) _load_helper_i(idv->get_sample(m), cur_i+cnt, m);
                        cnt++;
                        if (_missing > max_missing) return false;
                    }
                }
            }
            for (unsigned int i=0; i<struc->num_indiv_outgroup(); i++) {
                idv = & struc->get_indiv_outgroup(i);
                for (unsigned int j=0; j<_ploidy; j++) _load_helper_o(idv->get_sample(j), cur_o+i, j);
                if (consider_outgroup_missing && _missing > max_missing) return false;
            }
        }
        else {
            add_ing(data.get_nsam_i());
            add_otg(data.get_nsam_o());
            for (unsigned int i=0; i<data.get_nsam_i(); i++) {
                _load_helper_i(i, cur_i+i, 0);
                if (_missing > max_missing) return false;
            }
            for (unsigned int i=0; i<data.get_nsam_o(); i++) {
                _load_helper_o(i, cur_o+i, 0);
                if (consider_outgroup_missing && _missing > max_missing) return false;
            }
        }

        _tot_missing += _missing;
        return true;
    }

    void SiteHolder::_load_helper_i(unsigned int sam_idx, unsigned int idv, unsigned int chrom) {
        _ing[idv*_ploidy+chrom] = _process_allele(
            _filter_pointer->check(_data_pointer->get_i(sam_idx, _site_index), _flag), true);
        if (_ing[idv*_ploidy+chrom] == MISSING) {
            _missing++;
            _missing_ing++;
        }
        if (_flag) throw EggInvalidAlleleError(_data_pointer->get_i(sam_idx, _site_index), sam_idx, _site_index);
    }

    void SiteHolder::_load_helper_o(unsigned int sam_idx, unsigned int idv, unsigned int chrom) {
        _otg[idv*_ploidy+chrom] = _process_allele(
            _filter_pointer->check(_data_pointer->get_o(sam_idx, _site_index), _flag), false);
        if (_otg[idv*_ploidy+chrom] == MISSING && _consider_outgroup_missing) {
            _missing++;
            _missing_otg++;
        }
        if (_flag) throw EggInvalidAlleleError(_data_pointer->get_o(sam_idx, _site_index), sam_idx, _site_index);
    }

    unsigned int SiteHolder::_process_allele(int x, bool ing) {
        if (x == MISSINGDATA) {
            return MISSING;
        }

        for (unsigned int i=0; i<_nall; i++) {
            if (x == _alleles[i]) return i;
        }
        _nall++;

        if (ing) _nall_ing++;
        if (_nall > _nall_c) {
            _alleles = (int *) realloc(_alleles, _nall * sizeof(int));
            if (!_alleles) throw EGGMEM;
            _nall_c = _nall;
        }
        _alleles[_nall-1] = x;
        return _nall - 1;
    }

    bool SiteHolder::process_vcf(VcfParser& data,
                     unsigned int start, unsigned int stop,
                     unsigned int max_missing) {

        unsigned int cur = _ni;
        unsigned int cur_o = _no;

        _missing = 0;
        unsigned int allele;
        unsigned int input_ploidy = data.ploidy();
        int recoded;

        if (input_ploidy != _ploidy) throw EggArgumentValueError("ploidy mismatch in SiteHolder::process_vcf()");
        unsigned int num_outgroup = 0;
        for (unsigned int i=start; i<stop; i++) if (data.is_outgroup(i)) num_outgroup++;
        add_ing(stop-start-num_outgroup); // if AA is included in outgroup, do not count it here!
        if (data.outgroup_AA()) num_outgroup++;
        add_otg(num_outgroup);

        // get AA as outgroup
        unsigned int otg_idx = 0;
        if (data.outgroup_AA()) {
            recoded = MISSINGDATA;
            if (data.has_AA() and strcmp(data.AA_string(), "?") != 0) {
                if (strcmp(data.AA_string(), data.reference()) == 0) recoded = 0;
                for (unsigned int i=0; i<data.num_alternate(); i++) {
                    if (strcmp(data.AA_string(), data.alternate(i)) == 0) recoded = 1 + i;
                }
            }
            // set the alleles to all alleles of the genotypes (homozygote assumed)
            for (unsigned int i=0; i<input_ploidy; i++) _otg[cur_o+otg_idx++] = _process_allele(recoded, false);
        }

        unsigned int ing_idx = 0;
        for (unsigned int i=start; i<stop; i++) {
            for (unsigned int j=0; j<input_ploidy; j++) {
                allele = data.GT(i, j);
                if (allele == UNKNOWN) {
                    _missing++;
                    if (!data.is_outgroup(i)) _missing_ing++;
                    if (_missing_ing > max_missing) return false;
                    recoded = MISSINGDATA;
                }
                else {
                    recoded = static_cast<int>(allele);
                }

                if (data.is_outgroup(i)) {
                    _otg[cur_o+otg_idx++] = _process_allele(recoded, false);
                }
                else {
                    _ing[cur+ing_idx++] = _process_allele(recoded, true);
                }
            }
        }

        _tot_missing += _missing;

        return true;
    }

    SiteGeno::SiteGeno() : SiteHolder(1) {
        _sz_geno = 0;
        _ploidy = 0;
        _sz_ploidy1 = 0;
        _sz_ploidy = NULL;
        _homoz = NULL;
        _geno = NULL;
        _geno_matched = NULL;
    }

    SiteGeno::~SiteGeno() {
        if (_sz_ploidy) free(_sz_ploidy);
        if (_homoz) free(_homoz);
        if (_geno) {
            for (unsigned int i=0; i<_sz_geno; i++) if (_geno[i]) free(_geno[i]);
            free(_geno);
        }
        if (_geno_matched) free(_geno_matched);
    }

    bool SiteGeno::homoz(unsigned int genotype) const {
        return _homoz[genotype];
    }

    void SiteGeno::process(const SiteHolder& src) {
        reset(1);
        add_ing(src.get_ning());
        add_otg(src.get_nout());

        _nall = 0;
        _ploidy = src.get_ploidy();

        if (_ploidy > _sz_ploidy1) {
            _geno_matched = (bool *) realloc(_geno_matched, _ploidy * sizeof(bool));
            if (!_geno_matched) throw EGGMEM;
            _sz_ploidy1 = _ploidy;
        }

        for (unsigned int i=0; i<_ni; i++) _ing[i] = _analyse_genotype(src.get_pi(i));
        _nall_ing = _nall;
        for (unsigned int i=0; i<_no; i++) _otg[i] = _analyse_genotype(src.get_po(i));
    }

    unsigned int SiteGeno::_analyse_genotype(const unsigned int * p) {

        // check if any missing data
        unsigned int chr;
        for (chr=0; chr<_ploidy; chr++) {
            if (p[chr] == MISSING) return MISSING;
        }

        unsigned int idx, all_idx, all_idx2;
        for (idx=0; idx<_nall; idx++) {
            for (all_idx=0; all_idx<_ploidy; all_idx++) _geno_matched[all_idx] = false;
            for (all_idx=0; all_idx<_ploidy; all_idx++) {
                for (all_idx2=0; all_idx2<_ploidy; all_idx2++) {
                    if (_geno_matched[all_idx2]==false && _geno[idx][all_idx2] == p[all_idx]) {
                        _geno_matched[all_idx2] = true;
                        break; // found this allele
                    }
                }
                if (all_idx2 == _ploidy) break; // did not find a match for p[all_idx]
            }
            if (all_idx == _ploidy) {
                return idx; // all alleles are matching
            }
        }

        // add memory space for new genotype
        _nall++;
        if (_nall > _nall_c) {
            _alleles = (int *) realloc(_alleles, _nall * sizeof(int));
            if (!_alleles) throw EGGMEM;
            _nall_c = _nall;
        }
        if (_nall > _sz_geno) {
            _homoz = (bool *) realloc(_homoz, _nall * sizeof(bool));
            if (!_homoz) throw EGGMEM;
            _sz_ploidy = (unsigned int *) realloc(_sz_ploidy, _nall * sizeof(unsigned int));
            if (!_sz_ploidy) throw EGGMEM;
            _geno = (unsigned int **) realloc(_geno, _nall * sizeof(unsigned int *));
            if (!_geno) throw EGGMEM;
            _sz_ploidy[_sz_geno] = 0;
            _geno[_sz_geno] = NULL;
            _sz_geno = _nall;
        }

        // this is the second dimension
        if (_ploidy > _sz_ploidy[_nall-1]) {
            _geno[_nall-1] = (unsigned int *) realloc(_geno[_nall-1], _ploidy * sizeof(unsigned int));
            if (!_geno[_nall-1]) throw EGGMEM;
            _sz_ploidy[_nall-1] = _ploidy;
        }

        // copy genotype
        for (chr=0; chr<_ploidy; chr++) _geno[_nall-1][chr] = p[chr];
        _alleles[_nall-1] = static_cast<int>(_nall - 1); // hugh

        // check homozygosity and return
        for (chr=1; chr<_ploidy; chr++) {
            if (p[chr] != p[0]) {
                _homoz[_nall-1] = false;
                return idx;
            }
        }

        // default
        _homoz[_nall-1] = true;
        return idx;
    }
}
