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

#ifndef EGGLIB_HAPLOTYPES_HPP
#define EGGLIB_HAPLOTYPES_HPP

namespace egglib {

    class SiteHolder;
    class StructureHolder;

   /** \brief Identifies haplotypes from a set of sites
    *
    * \ingroup diversity
    *
    * How to use this class:
    *
    * * 1) Setup (or reset) and an optional structure.
    * * 2) Load all sites with load(site). Haplotypes are computed and all
    *      samples with at least one missing data are marked as missing.
    *      While loading sites, you can monitor the values of:
    *      * n_ing()
    *      * n_otg()
    *      * n_sam()
    *      * nt_hapl()
    *      * hapl()
    *      * map()
    *      * n_sites()
    * * 3) Call cp_haplotypes() to finalize haplotype processing. After
    *      that you may use:
    *      * ne_ing()
    *      * ne_otg()
    *      * ni_hapl()
    *      * ng_hapl()
    *      * freq_i()
    *      * freq_o()
    *      * nstot()
    * * 4) If you wish (but you don't have to), you can try and guess the
    *      haplotype of samples with missing data. For this you need first
    *      to call prepare_impute(). After that, you may use:
    *      * n_mis()
    *      * mis_idx()
    * * 5) If you impute, load again all sites with solve(). This may
    *      change the value of:
    *      * map()
    * * 6) If you impute, you are required to call impute() (and
    *      otherwise you can't). If you do, the following values will be
    *      updated [note: maybe impute and cp_haplotypes do exactly the same]
    *      * ni_hapl()
    *      * ng_hapl()
    *      * freq_i()
    *      * freq_o()
    *      and these values will be set:
    *      * num_potential()
    *      * potential()
    * * 7) After whether or not you performed 4-6, you can now call this
    *      to make a site:
    *      * get_site()
    * * 8) Whether or not you performed 4-6 and/or 7, you can now call
    *      cp_dist() that will let you access to the distance matrix:
    *      * dist()
    * * 9) To compute stats, call cp_stats() (you must have set a structure
    *       and called 8). Then you may use:
    *       * Fst()
    *       * Kst()
    *       The function cp_stats() returns a flag: 0 (no stats computed),
    *       1 (Fst computed), or 2 (both Fst and Kst computed).
    * * 10) You can also call this method (requires a structure and 8 but not 9, and
    *       valid if nstot is >1):
    *       * Snn()
    *
    * Header: <egglib-cpp/Haplotypes.hpp>
    *
    */
    class Haplotypes {

        private:
            Haplotypes(const Haplotypes& src) {}
            Haplotypes& operator=(const Haplotypes& src) {return *this;}
            void _init();
            void _free();
            void _add_hapl();
            void _process(unsigned int allele, unsigned int index);
            void _set_sample_size(unsigned int, unsigned int, unsigned int);
            bool _invalid;
            const StructureHolder * _struct;
            unsigned int _pl;
            unsigned int _ni_tot;
            unsigned int _no_tot;
            unsigned int _nsi;
            unsigned int _nso;
            unsigned int _ne_ing;
            unsigned int _ne_otg;
            unsigned int _n_sam;
            unsigned int _c_sam;
            unsigned int _nt_hapl;
            unsigned int   _c_hapl;
            unsigned int * _c_hapl2;     // max size: _n_mis
            unsigned int   _c_hapl3;
            unsigned int * _c_hapl4;     // max size: _c_hapl3
            unsigned int _ng_hapl;
            unsigned int _ni_hapl;
            unsigned int _n_sites;
            unsigned int * _c_sites;     // max size: _c_hapl
            unsigned int * _freq_i;      // max size: _c_hapl
            unsigned int * _freq_o;      // max size: _c_hapl
            unsigned int * _map;         // max size: _c_sam
            unsigned int * _map_cache;   // max size: _c_sam
            unsigned int ** _hapl;       // max size: _c_hapl * _c_sites[i]
            unsigned int _n_mis;
            unsigned int _c_mis;
            unsigned int * _mis_idx;     // max_size: _c_mis
            unsigned int ** _potential;  // max size: _c_mis * _c_hapl2[i] (truly allocated=_n_hapl+1)
            unsigned int * _n_missing;   // max size: _c_sam
            unsigned int ** _dist;       // max size: _c_hapl3 * _c_hapl4[i] (->i)
            unsigned int _n_pop;
            unsigned int _ne_pop;
            unsigned int _c_pop;
            unsigned int * _Ki;          // max size: _c_pop;
            unsigned int * _ni;          // max size: _c_pop (+1 for sum);
            unsigned int * _Kd;          // max size: _c_pop * (_c_pop - 1);
            unsigned int * _pop_i;       // max size: _c_sam
            unsigned int * _pop_o;       // pointer to _pop[i][_n_ing]
            unsigned int _nstot;
            double _Fst;
            double _Kst;
            unsigned int _site_index;

        public:
            Haplotypes(); ///< Constructor
            ~Haplotypes(); ///< Destructor
            void setup(const StructureHolder * struc = NULL); ///< Setup/reset instance (with optional structure object)
            void load(const SiteHolder&); ///< Process a site (implies call to setup())
            void prepare_impute(unsigned int); ///< Initialize required tables for imputing
            void resolve(const SiteHolder&); ///< Second pass: try to resolve missing data
            void cp_haplotypes(); ///< Finalize haplotype estimation.
            void impute(); ///< Try to guess haplotype of samples with missing data.
            void cp_dist(); ///< Compute distance matrix
            unsigned int cp_stats(); ///< Compute differentiation stats
            unsigned int n_ing() const; ///< Total number of ingroup samples
            unsigned int n_otg() const; ///< Total number of outgroup samples
            unsigned int n_sam() const; ///< Total number of samples
            unsigned int nstot() const; ///< Number of samples (non-missing and not ignored by structure)
            unsigned int nt_hapl() const; ///< Total number of haplotypes (including truncated ones)
            unsigned int hapl(unsigned int, unsigned int) const; ///< Get site j of haplotype i
            unsigned int map_sample(unsigned int) const; ///< Get haplotype index of a sample
            unsigned int n_sites() const; ///< Number of sites
            unsigned int ne_ing() const; ///< Non-missing number of ingroup samples
            unsigned int ne_otg() const; ///< Non-missing number of outgroup samples
            unsigned int n_missing(unsigned int) const; ///< Number of missing data per sample
            unsigned int n_mis() const; ///< Number of samples with missing data
            unsigned int mis_idx(unsigned int) const; ///< Index of one of the samples with missing data
            unsigned int ni_hapl() const; ///< Number of haplotypes at non-zero ingroup frequency
            unsigned int ng_hapl() const; ///< Number of haplotypes at non-zero frequency overall
            unsigned int freq_i(unsigned int) const; ///< Frequency of haplotype in intgroup
            unsigned int freq_o(unsigned int) const; ///< Frequency of haplotype in outgroup
            unsigned int n_potential(unsigned int) const; ///< Number of compatible haplotypes for a sample with missing data
            unsigned int potential(unsigned int, unsigned int) const; ///< Index of a potential haplotype
            unsigned int dist(unsigned int, unsigned int) const; ///< Distance matrix entry (0<=j<i<_nt_hapl)
            unsigned int n_pop() const; ///< Number of populations
            unsigned int ne_pop() const; ///< Number of populations with 
            double Fst() const; ///< Fst value
            double Kst() const; ///< Gst value
            double Snn() const; ///< Snn value computed on the fly
            void get_site(SiteHolder&); ///< Get haplotypic data as a site (recycle passed object)
    };
}

#endif
