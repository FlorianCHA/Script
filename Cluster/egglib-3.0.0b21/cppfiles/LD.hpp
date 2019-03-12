/*
    Copyright 2009-2017 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_LINKAGE_DISEQUILIBRIUM_HPP
#define EGGLIB_LINKAGE_DISEQUILIBRIUM_HPP

#include "SiteHolder.hpp"
#include "Structure.hpp"

namespace egglib {

   /** \brief Analyzes linkage disequilibrium for a pair of polymorphic sites
    *
    * \ingroup diversity
    *
    * This class considers a single pair of polymorphic sites at a
    * time. The first method, process(), detects alleles at both sites
    * under consideration and determines whether the pairwise
    * comparison is fit for analysis (based on the presence of
    * polymorphism, and allele frequencies). Statistics are computed
    * by compute() for a given pair of alleles. Letting the user
    * filter out sites for which that are more than two alleles and,
    * if necessary, process multiple pairs of alleles.
    *
    * One should first process a pair of sites with process(). If the
    * return value is false, one should not process data further.
    * Otherwise, one can access data with num_alleles1(), num_alleles2(),
    * index1(), index2(), freq1(), freq2(), freq(), and nsam(), and can
    * also compute LD with compute() (for a given pair of alleles) and
    * then access to LD estimates.
    *
    */
    class PairwiseLD {

        private:
            PairwiseLD(const PairwiseLD& src) {}
            PairwiseLD& operator=(const PairwiseLD& src) {return * this;}

            double _D;
            double _Dp;
            double _r;
            double _rsq;
            SiteGeno _site1_geno;
            SiteGeno _site2_geno;
            const SiteHolder * __site1__;
            const SiteHolder * __site2__;

            void _alloc(unsigned int a1i, unsigned a2i);
            void _helper(unsigned int);

            unsigned int _a1;        // number of alleles site 1 (initial)
            unsigned int _a2;        // number of alleles site 2 (initial)
            unsigned int _a1e;       // number of alleles site 1 (effective)
            unsigned int _a2e;       // number of alleles site 2 (effective)
            unsigned int _neff;      // effective number of samples
            unsigned int * _map1;    // (effective) allele indexes for site 1
            unsigned int * _map2;    // (effective) allele indexes for site 2

            unsigned int _a1c;       // allocated number of alleles site 1
            unsigned int _a2c;       // allocated number of alleles site 2
            unsigned int * _pc;      // allocated size per row in **p
            unsigned int * _p1;      // frequency of alleles site 1
            unsigned int * _p2;      // frequency of alleles site 2
            unsigned int ** _p;      // frequency of haplotypes

        public:

            PairwiseLD(); ///< Default constructor
            ~PairwiseLD(); ///< Destructor

           /** \brief Analyze a pair of sites
            *
            * The method takes two sites as argument. The two sites
            * must be taken from the same data set. In particular, the
            * sample sizes must be identical. %Structure and outgroup
            * are ignored. The indexes of samples must be matching
            * over the two sites. Samples which are missing in either
            * of the samples are skipped. If the remaining samples are
            * less than the argument min_n, the whole computation is
            * dropped. Genotypes are ignored (only alleles are
            * considered).
            *
            * \param site1 first site.
            *
            * \param site2 second site.
            *
            * \param min_n minimum number of samples used (this value
            * must always be larger than 1).
            *
            * \param max_maj maximum relative frequency of the
            * majority allele (if any allele at either site has a
            * frequency larger than this value, the pairwise
            * comparison is dropped).
            *
            * \param stru a Structure object (used to process a subset of
            * samples). By default, use all samples.
            *
            * \return true if computations have been performed, false
            * if the sites fall in one the following cases: not enough
            * samples (based on the min_n argument); either site is
            * fixed; the allele frequencies are too unbalanced with at
            * least one allele at a frequency larger than max_maj.
            *
            * \note Due to missing data, a site that is initially
            * polymorphic might appear to be fixed when considering
            * only samples that are not missing for the other site,
            * causing this method to drop the pairwise comparison.
            * Conversely, a site that has more than two alleles might
            * have only two when considering only samples that are not
            * missing for the other site. For this reason, it is not
            * trivial to filter out sites before calling this method,
            * and sites might not be consistently included or
            * rejected.
            *
            * If this method returns true, statistics might be
            * computed for a given pair of alleles using the compute()
            * method. The number of alleles available for analysis is
            * available at either site using num_alleles1() and
            * num_alleles2(). When returning false, this method stops
            * as early as possible, and the state of the object might
            * be inconsistent. In this case, no accessor must be used
            * and compute() must not be called.
            *
            */
            bool process(const SiteHolder& site1, const SiteHolder& site2, unsigned min_n = 2,
                        double max_maj = 1.0, const StructureHolder * stru=NULL);

           /** \brief Get the actual number of alleles at the first site
            *
            * The method process() must have been executed and must
            * have returned true.
            *
            * Gives the number of different alleles at the first site,
            * considering only samples for which both sites have
            * exploitable data.
            *
            */
            unsigned int num_alleles1() const;

           /** \brief Get the actual number of alleles at the second site
            *
            * The method process() must have been executed and must
            * have returned true.
            *
            * Gives the number of different afirst lleles at the
            * second site, considering only samples for which both
            * sites have exploitable data.
            *
            */
            unsigned int num_alleles2() const;

           /** \brief Get the index of an allele for the first site
            *
            * The method process() must have been executed and must
            * have returned true.
            *
            * For a given allele, get its index within the original
            * SiteHolder instance. The indexes can be shifted by process()
            * due to missing data.
            *
            */
            unsigned int index1(unsigned int allele) const;

           /** \brief Get the index of an allele for the second site
            *
            * The method process() must have been executed and must
            * have returned true.
            *
            * For a given allele, get its index within the original
            * SiteHolder instance. The indexes can be shifted by process()
            * due to missing data.
            *
            */
            unsigned int index2(unsigned int allele) const;

           /** \brief Get the frequency of an allele for the first site
            *
            * The method process() must have been executed and must
            * have returned true.
            *
            * The index must be smaller than the value returned by
            * num_alleles1().
            *
            */
            unsigned int freq1(unsigned int allele) const;

           /** \brief Get the frequency of an allele for the second site
            *
            * The method process() must have been executed and must
            * have returned true.
            *
            * The index must be smaller than the value returned by
            * num_alleles2().
            *
            */
            unsigned int freq2(unsigned int allele) const;

           /** \brief Get the frequency of a genotype
            *
            * The method process() must have been executed and must
            * have returned true.
            *
            * The indexes must be smaller than the value returned by
            * num_alleles1() and num_alleles2() respectively.
            *
            */
            unsigned int freq(unsigned int allele1, unsigned int allele2) const;

           /** \brief Get the number of analyzed samples
            *
            * The method process() must have been executed and must have
            * returned true.
            *
            * The returned value might be smaller than the initial
            * number of samples due to missing data.
            *
            */
            unsigned int nsam() const;

           /** \brief Compute D, D', r and r^2 statistics for a given pair of alleles
            *
            * The method process() must have been executed and must
            * have returned true.
            *
            * Statistics are computed only for a given pair of
            * alleles. If there are only two alleles, all allele pairs
            * result in consistent results. Otherwise, some
            * multi-allele summarizing methodology has to be applied.
            *
            * allele1 and allele2 are the allele indexes at the first
            * and second site, respectively.
            *
            */
            void compute(unsigned int allele1, unsigned int allele2);

           /** \brief Get the D statistic
            *
            * This value is reset to 0 upon call to process().
            *
            * Requires compute().
            *
            */
            double D() const;

           /** \brief Get the D' statistic
            *
            * This value is reset to 0 upon call to process().
            *
            * Requires compute();
            *
            */
            double Dp() const;

           /** \brief Get the r statistic
            *
            * This value is reset to 0 upon call to process().
            *
            * Requires compute().
            *
            */
            double r() const;

           /** \brief Get the r^2 statistic
            *
            * Same as r()*r(). This value is reset to 0 upon call to
            * process().
            *
            *
            * Requires compute().
            *
            */
            double rsq() const;

           /** \brief Reset all values to default
            *
            * Call to this method is usually not necessary since process()
            * automatically resets the instance.
            *
            */
            void reset();
    };

    //////////////////////////////////////////////////////////////////

   /** \brief Analyzes linkage disequilibrium between pairs of sites
    *
    * \ingroup diversity
    *
    * This class processes a set of SiteHolder instances and computes
    * linkage disequilibrium for all pairs of sites. A PairwiseLD
    * instance is provided for all comparison, skipping all pairs for
    * which LD cannot be computed (there are several criteria). The
    * approach consists in first calling load() by providing a set
    * of SiteHolder instances. The method computeLD() computes
    * the LD for each pair and  computeStats() computes the statistics
    * of Kelly (1997) and Rozas et al. (2001). These statistics are
    * based on the average of pairwise linkage disequilibrium
    * statistics. In addition, computeRmin() computes Rm of Hudson and
    * Kaplan (1985) and does not generate nor use PairwiseLD instances,
    * and it can be used independently.
    *
    */
    class MatrixLD {

        private:
            MatrixLD(const MatrixLD& src) {}
            MatrixLD& operator=(const MatrixLD& src) {return *this;}

            unsigned int _nsites;
            unsigned int _nsites_c;
            const SiteHolder ** _sites;
            double * _positions;
            unsigned int * _index1;
            unsigned int * _index2;
            PairwiseLD ** _linkage;
            unsigned int * _dist;
            bool * _adjacent;
            double _ZnS;
            double _ZnS_star1;
            double _ZnS_star2;
            double _Za;
            double _ZZ;
            unsigned int _num_allele_pairs;
            unsigned int _num_allele_pairs_adj;
            unsigned int _ntot;
            unsigned int _npairs;
            unsigned int _np_c;
            unsigned int _nalls;
            unsigned int _Rmin;
            unsigned int _Rmin_res_sites;
            unsigned int _Rmin_res_intervals;
            unsigned int _Rmin_num_sites;
            bool * _Rmin_sites;
            unsigned int * _Rmin_left;
            unsigned int * _Rmin_right;
            bool * _Rmin_bool;
            unsigned int _nall_c;
            bool * _all_flags;
            bool _toggle_stats;
            bool _toggle_Rmin;
            unsigned int _ns;
            unsigned int _pl;
            bool _mismatch;

        public:

           /** \brief Flags for processing multiallelic sites
            *
            * This enum is used to specify what should be done with
            * pairs of sites for which at least one site has more than
            * two alleles.
            *
            */
            enum MultiAllelic {
                ignore,     ///< Only process pairs of sites with exactly two alleles
                use_main,   ///< Use the allele with highest frequency
                use_all     ///< Use all alleles
            };

            MatrixLD(); ///< Constructor
            ~MatrixLD(); ///< Destructor
            void reset(); ///< Reset to defaults

            void toggle_off(); ///< Toggle all off
            void toggle_Rmin(); ///< Toggle Rmin
            void toggle_stats(); ///< Toggle summary statistics

           /** \brief Load a site
            *
            * \param site One of the sites.
            * 
            * \param position All sites must have a valid position.
            * Positions are required to be increasing. For computing Rmin,
            * positions are ignored (they only are fed back if interval
            * limits are required).
            *
            */
            void load(const SiteHolder& site, double position);

           /** \brief Compute LD between all pairs of sites
            *
            * Use sites loaded using load() and process all possible
            * pairs. Each pairwise comparison is retained only if all
            * filters are passed (see arguments of this method). After
            * call of this method, the number of pairs can be accessed
            * using num_tot() (it is equal to n(n-1)/2 where n is the
            * number of loaded sites); the number of analyzed pairs
            * can be accessed using num_pairs(); the total number of
            * analyzed allele pairs can be accessed using
            * num_alleles(). Then the method compute() can be called
            * to compute Kelly's statistics.
            *
            * \param min_n minimum number of samples used (this value
            * must always be larger than 1).
            *
            * \param max_maj maximum relative frequency of the
            * majority allele (if any allele at either site has a
            * frequency larger than this value, the pairwise
            * comparison is dropped).
            *
            * \note Due to missing data, it is not trivial to predict
            * whether a pairwise comparison will be dropped. See the
            * documentation of PairwiseLD::process().
            *
            */
            void computeLD(unsigned min_n = 2, double max_maj = 1.0);

           /** \brief Get the total number of pairs of sites
            *
            * Requires that site pairs have been processed using
            * computeLD().
            *
            */
            unsigned int num_tot() const;

           /** \brief Get the number of analyzed pairs of sites
            *
            * Requires computeLD(). The returned value excludes all
            * pairwise comparisons with no polymorphism failing any
            * other criterion (see the computeLD() method).
            *
            */
            unsigned int num_pairs() const;

           /** \brief Get the total number of allele pairs
            *
            * Requires computeLD().  Returns the sum of allele pairs
            * over all analyzed sites (see num_pairs()). The minimum
            * value is twice num_pairs() (since, by definition, there
            * must be at least two alleles at each retained site).
            *
            */
            unsigned int num_alleles() const;

           /** \brief Get linkage disequilibrium for a given pair of sites
            *
            * Requires that enough pairs have been loaded using load()
            * and that requested index must be smaller than
            * num_pairs(). Use methods index1() and index2() to obtain
            * the corresponding site indexes.
            *
            */
            const PairwiseLD& pairLD(unsigned int index) const;

           /** \brief Index of first site for a given pair.
            *
            * See pairLD().
            *
            */
            unsigned int index1(unsigned int allele) const;

           /** \brief Index of second site for a given pair.
            *
            * See pairLD(). Note that index2 is always > index1.
            *
            */
            unsigned int index2(unsigned int allele) const;

           /** \brief Get the distance between sites for a given pair
            *
            * Requires that enough pairs have been loaded using load()
            * and that requested index must be smaller than
            * num_pairs(). The distance is returned as an absolute
            * value.  It is not possible to determine to which sites
            * the pair index corresponds. If you need it, you might
            * want to use PairwiseLD directly.
            *
            */
            unsigned int distance(unsigned int index) const;

           /** \brief Call {computeLD() and computeStats()} and/or computeRmin() based on toggled flags
            *
            * Arguments are like for the three methods. Enter anything if they are not used.
            *
            * Return value is a flag with the following bits:
            * * 0: ZnS, Z*nS, and Z*nS* are computed.
            * * 1: Za and ZZ are computed.
            * * 2: Rmin was computed.
            *
            */
            unsigned int process(unsigned min_n, double max_maj, MultiAllelic multiallelic, unsigned int min_freq,
                bool oriented);

           /** \brief Computes Kelly's and Rozas et al.'s statistics
            *
            * Computes ZnS, Z*nS and Z*nS* (Kelly 1997), and Za and ZZ
            * (Rozas et al. 2001) on the basis of analyzed site pairs
            * (requires computeLD()). The number of alleles pairs used
            * for computing ZnS, Z*nS and Z*nS* is available as
            * num_allele_pairs(), and the number of allele pairs used
            * for computing Za and ZZ is available as
            * num_allele_pairs_adj(). The statistics must not be used if
            * the corresponding number of allele pairs is 0.
            *
            * If multiallelic equals to MatrixLD::ignore, only pairs
            * of sites for which both sites have exactly two alleles
            * are processed. In this case, the first allele of each
            * site is considered. If multiallelic is
            * MatrixLD::use_main, the alleles with highest frequency
            * are considered (even if one or both sites have only two
            * alleles). In case of equality, the first allele is
            * considered. If multiallelic is MatrixLD::use_all, then
            * all alleles of all sites are used, and the final
            * statistics are averaged over num_alleles() (rather than
            * num_pairs).
            *
            * \param multiallelic modifies the behaviour of the method
            * (see above).
            *
            * \param min_freq this flag has an effect only if used in
            * conjunction with MatrixLD::use_all (it is ignored
            * otherwise); if larger than 0, rather than using all
            * alleles, use only those that have a frequency equal to
            * or larger than the given value.
            *
            */
            void computeStats(MultiAllelic multiallelic = ignore, unsigned int min_freq = 0);

           /** \brief Number of allele pairs used to compute Kelly's statistics
            * 
            * Return the number of allele pairs used by computeStats()
            * to compute Kelly's ZnS, Z*nS and Z*nS* statistics. If
            * multiallelic equals ignore, this value equals the number
            * of pairs of sites with exactly two alleles each (at most,
            * num_pairs()); if multiallelic was use_main, this value
            * equals num_pairs(); if multiallelic was use_all, this
            * value equals num_alleles(). If the returned value is 0 (no
            * loaded pairs of sites, or no pairs of diallelic sites, if
            * multiallelic was set to MatrixLD::ignore), Kelly's
            * statistics have been reset to 0 but should then be
            * considered as not computable. If num_allele_pairs() is
            * null, none of the Kelly's and Rozas et al.'s statistics
            * can be computed.
            *
            */
            unsigned int num_allele_pairs() const;

           /** \brief Number of allele pairs used to compute Rozas et al.'s statistics
            * 
            * Return the number of allele pairs used by computeStats()
            * to compute Rozas et al.'s Za and ZZ statistics. See the
            * documentation of num_allele_pairs() for reference. The
            * meaning of this value is similar, except that it applies
            * only to adjacent polymorphic sites (the value hence can
            * only be smaller, or equal is limit cases). If this vallue
            * is 0, Rozas et al.'s statistics have been reset to 0 but
            * should then be considered as not computable.
            *
            */
            unsigned int num_allele_pairs_adj() const;

           /** \brief Get Kelly's ZnS statistic
            *
            * Requires computeStats(). See documation of this method to
            * known when this value is defined.
            *
            */
            double ZnS() const;

           /** \brief Get Kelly's Z*nS statistic
            *
            * Requires computeStats(). See documation of this method to
            * known when this value is defined.
            *
            */
            double ZnS_star1() const;

           /** \brief Get Kelly's Z*nS* statistic
            *
            * Requires computeStats(). See documation of this method to
            * known when this value is defined.
            *
            */
            double ZnS_star2() const;

           /** \brief Get Rozas et al.s Za statistic
            *
            * Requires computeStats(). See documation of this method to
            * known when this value is defined.
            *
            */
            double Za() const;

           /** \brief Get Rozas et al.s ZZ statistic
            *
            * Requires computeStats(). See documation of this method to
            * known when this value is defined.
            *
            */
            double ZZ() const;

           /** \brief Computes Hudson and Kaplan's Rmin
            *
            * To be used, this method requires that sites have been
            * loaded in increasing position order. Only sites with
            * exactly two alleles and no missing data at all are used.
            * In addition, if the oriented argument is set to true, only
            * orientable sites are considered. Sites with more than two
            * alleles and sites with any missing data, and sites not
            * orientable if oriented is set to true, are ignored.
            *
            * When the method has finished, a few methods provide access
            * to the results. Rmin_num_sites() give the number of sites
            * considered for the analysis. If the value is less than
            * two, the statistic itself bears no signification. Rmin()
            * gives the minimum number of recombination events. Finally,
            * all of the non-reductible intervals containing a
            * recombination even can be accessed using the two methods
            * Rmin_left(unsigned int) and Rmin_right(unsigned int). The
            * number of intervals is always Rmin().
            *
            * \param oriented if true, consider only orientable sites
            * and apply the three- (instead four-) gametes rule. If
            * false, ignore all outgroup data and include orientable and
            * non-orientable sites.
            *
            */
            void computeRmin(bool oriented=false);

           /** \brief Minimal number of recombination events
            * 
            * Requires computeRmin().
            * 
            */
            unsigned int Rmin() const;

           /** \brief Number of sites used for computing Rmin
            *
            * 
            * Requires computeRmin().
            *
            * Fixed to 0 if there are less than two sites in total (no
            * computation is performed by computeRmin() in that case).
            * If Rmin_num_sites() is less than 2, Rmin() is not defined
            * and fixed to 0.
            *
            */
            unsigned int Rmin_num_sites() const;

           /** \brief Left bound of a recombination interval
            * 
            * Requires computeRmin() and that Rmin_num_sites() is at
            * least 2.
            *
            * \param i interval index (there are Rmin() intervals).
            *
            * \return Position of the site at the 5' end of the given
            * interval (as provided to load()).
            *
            */
            unsigned int Rmin_left(unsigned int i) const;

           /** \brief Right bound of a recombination interval
            * 
            * Requires computeRmin() and that Rmin_num_sites() is at
            * least 2.
            *
            * \param i interval index (there are Rmin() intervals).
            *
            * \return Position of the site at the 5' end of the given
            * interval (as provided to load()).
            *
            */
            unsigned int Rmin_right(unsigned int i) const;
    };
}

#endif
