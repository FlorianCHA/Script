/*
    Copyright 2012-2017 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_SITEDIVERSITY_HPP
#define EGGLIB_SITEDIVERSITY_HPP

namespace egglib {

    class FreqHolder;

   /** \brief %Diversity analyses at the level of a site
    *
    * \ingroup diversity
    *
    * Computes standard diversity indexes for a unique site or marker.
    *
    * process() and average() return a composite flag.
    *
    * Statistics:
    *   * If fstats_diplo is called
    *           * npop_eff2 (pops with >= 1 indiv)
    *   * If fstats_haplo is called
    *           * npop_eff3 (pops with >= 1 sample)
    *   * If fstats_hier is called
    *           * nclu
    *           * nclu_eff (>= 1 pops each with >= 1 indiv)
    *           * npop_eff2 (same as for fstats_diplo)
    *   * flag&1 (always on for process()):
    *           * ns
    *   * flag&2:
    *           * npop
    *           * Aglob
    *           * Aing
    *           * Stot
    *           * pairdiff
    *           * He
    *           * R
    *           * npop_eff1 (pops with >= 2 samples)
    *           * He[pop] (for pops with >= 2 samples)
    *           * pairdiff_pop[pop1][pop2] (for pops with >= 2 samples and pop2 != pop1)
    *   * flag&4:
    *           * thetaIAM
    *           * thetaSMM
    *   * flag&8:
    *           * Ho
    *           * Hi
    *   * flag&16:
    *           * Sder
    *           * der
    *   * flag&2048:
    *           * Aout
    *   * flag&32:
    *           * n
    *           * d
    *   * flag&64:
    *           * a
    *           * b
    *           * c
    *   * flag&128:
    *           * a0
    *           * b1
    *           * b2
    *           * c0
    *   * flag&256:
    *           * JostD
    *   * flag&2048:
    *           * Hst
    *   * flag&4096:
    *           * Gst
    *   * flag&8192:
    *           * Gste
    *   * flag&512: ns, Aglob, Aing, Aout, Stot, Sder, and der are actually integers
    *   * flag&1024: site is polymorphic / there is at least one polymorphic site
    *           * MAF
    *           * MAF_pop
    *
    * Fit = 1 - c/(a+b+c)
    * Fst = a/(a+b+c)
    * Fis = 1 - c/(b+c)
    *
    * Fst = n/d
    *
    * Fit = 1 - c0/(a0+b2+b1+c0)
    * Fst = (a0+b2)/(a0+b2+b1+c0)
    * Fct = a0/(a0+b2+b1+c0)
    * Fis = 1 - c0/(b1+c0)
    *
    * Hst = 1 - Hs / He
    * Gst = 1 - Hs / Httilde
    * Gste = 1 - Hse / Hte
    *
    * Requires: stats()
    *
    */
    class SiteDiversity {

        public:

            SiteDiversity();                  ///< Constructor
            virtual ~SiteDiversity();         ///< Destructor
            void toggle_off();                ///< Set all flags to off
            void toggle_fstats_diplo();       ///< Toggle F-statistics
            void toggle_fstats_haplo();       ///< Toggle F-statistics
            void toggle_fstats_hier();        ///< Toggle F-statistics
            void toggle_hstats();             ///< Toggle H-statistics
            unsigned int flag() const;        ///< Get flag value
            unsigned int process(const FreqHolder&);  ///< Compute toggled statistics
            unsigned int average();           ///< Compute the average of all stats (except those per pop)
            double ns() const;                ///< Number of analyzed samples (stats)
            unsigned int k() const;           ///< Number of populations (stats)
            unsigned int npop_eff1() const;   ///< Number of populations with >= 2 samples (stats)
            unsigned int npop_eff2() const;   ///< Number of populations with >= 1 indiv (fstats_diplo + fstats_hier)
            unsigned int npop_eff3() const;   ///< Number of populations with >= 1 sample (fstats_haplo)
            unsigned int nclu_eff() const;    ///< Number of clusters with >= 1 pop with >= 1 indiv (fstats_hier)
            double Aglob() const;             ///< Number of alleles (including outgroup-specific alleles) (stats)
            double Aing() const;              ///< Number of alleles excluding outgroup-specific alleles (stats)
            double Aout() const;              ///< Number of different alleles in the outgroup (stats)
            int global_allele(unsigned int) const;    ///< Get one of the global alleles
            double S() const;                 ///< Number of alleles at frequency one (singletons) (stats)
            double Sd() const;                ///< Number of derived singletons (stats) (requires outgroup)
            double R() const;                 ///< Allelic richness (stats)
            double pairdiff() const;          ///< Average number of pairwise differences (stats)
            double He() const;                ///< Unbiased heterozygosity (averaged if relevant) stats()
            double Ho() const;                ///< Frequency of heterozygotes (stats)
            double Hi() const;                ///< Avg number of differents between individuals (stats)
            double He_pop(unsigned int pop) const;  ///< Unbiased heterozygosity for a population stats()
            double pairdiff_inter(unsigned int pop1, unsigned int pop2) const;   ///< Average number of differences between a pair of population (stats)
            bool orientable() const;          ///< True if the site is orientable (stats)
            unsigned int num_derived() const; /// < Number of derived alleles (stats+outgroup)
            double derived(unsigned int) const; ///< Derived allele frequency (stats+outgroup)
            double a() const;              ///< Computed by fstats_diplo()
            double b() const;              ///< Computed by fstats_diplo()
            double c() const;              ///< Computed by fstats_diplo()
            double a0() const;             ///< Computed by fstats_hier()
            double b1() const;             ///< Computed by fstats_hier()
            double b2() const;             ///< Computed by fstats_hier()
            double c0() const;             ///< Computed by fstats_hier()
            double n() const;              ///< Computed by fstats_haplo()
            double d() const;              ///< Computed by fstats_haplo()
            double Hst() const;            ///< Computed by hstats()
            double Gst() const;            ///< Computed by hstats()
            double Gste() const;           ///< Computed by hstats()
            double D() const;              ///< Computed by hstats()
            double thetaIAM() const;       ///< Requires stats()
            double thetaSMM() const;       ///< Requires stats()
            void reset();                  ///< Reset stats sums to 0 (keep toggled flags)
            unsigned int nsites1() const;  ///< For average ns
            unsigned int nsites2() const;  ///< For average Aglob, Aing, Atot, Stot, pairdiff, He, thetaIAM, and thetaSMM.
            unsigned int nsites3() const;  ///< For average thetaIAM and thetaSMM
            unsigned int nsites4() const;  ///< For average Ho and Hi
            unsigned int nsites5() const;  ///< For average derived and Sd
            unsigned int nsites6() const;  ///< For average n and d
            unsigned int nsites7() const;  ///< For average a, b, and c
            unsigned int nsites8() const;  ///< For average c0, b1, b2, a0
            unsigned int nsites9() const;  ///< For average D
            unsigned int nsites10() const; ///< For average Hst
            unsigned int nsites11() const; ///< For average Gst
            unsigned int nsites12() const; ///< For average Gste
            double MAF() const;         ///< Frequency of second most frequent allele
            double MAF_pop(unsigned int) const;    ///< MAF for a population

        private:

            void _stats(const FreqHolder& freqs);
            void _fstats_diplo(const FreqHolder& freq);
            void _fstats_haplo(const FreqHolder& freqs);
            void _fstats_hier(const FreqHolder& freqs);
            void _hstats(const FreqHolder& freqs);

            SiteDiversity(const SiteDiversity& src) {}
            SiteDiversity& operator=(const SiteDiversity& src) {return *this;}

            void init();
            void free();
            void alloc();
            unsigned int _npop;
            unsigned int _nclu;
            unsigned int _npop_c;
            unsigned int _npop_eff1;
            unsigned int _npop_eff2;
            unsigned int _npop_eff3;
            unsigned int _nclu_eff;
            double _ns;
            double _A_tot;
            double _A_ing;
            double _A_glo;
            double _A_out;
            int * _glob_alleles;
            unsigned int _glob_c;
            double _S_tot;
            double _S_der;
            double _R;
            double _pairdiff;
            double ** _pairdiff_inter;
            double _He;
            double _thetaIAM;
            double _thetaSMM;
            double _Ho;
            double _Hi;
            double * _He_pop;
            double * _der;
            unsigned int _n_der;
            unsigned int _c_der;
            double _a;
            double _b;
            double _c;
            double _a0;
            double _b1;
            double _b2;
            double _c0;
            double _n;
            double _d;
            unsigned int * _ni;
            bool * _clu_flags;
            double * _pi;
            unsigned int _r_c;
            bool _orientable;
            double _Hst;
            double _Gst;
            double _Gste;
            double _JostD;
            bool _flag_fdip;
            bool _flag_fhap;
            bool _flag_fhie;
            bool _flag_h;
            unsigned int _flag;
            unsigned int _nsites1;
            unsigned int _nsites2;
            unsigned int _nsites3;
            unsigned int _nsites4;
            unsigned int _nsites5;
            unsigned int _nsites6;
            unsigned int _nsites7;
            unsigned int _nsites8;
            unsigned int _nsites9;
            unsigned int _nsites10;
            unsigned int _nsites11;
            unsigned int _nsites12;
            double _ns_acc;
            double _Aglob_acc;
            double _Aing_acc;
            double _Aout_acc;
            double _S_acc;
            double _R_acc;
            double _pairdiff_acc;
            double _He_acc;
            double _Ho_acc;
            double _Hi_acc;
            double _Sd_acc;
            double _a_acc;
            double _b_acc;
            double _c_acc;
            double _a0_acc;
            double _b1_acc;
            double _b2_acc;
            double _c0_acc;
            double _n_acc;
            double _d_acc;
            double _Hst_acc;
            double _Gst_acc;
            double _Gste_acc;
            double _D_acc;
            double _thetaIAM_acc;
            double _thetaSMM_acc;
            double _MAF;
            double * _MAF_pop;
    };
}

#endif
