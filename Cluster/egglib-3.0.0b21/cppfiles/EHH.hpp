/*
    Copyright 2013-2017 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_EHH_HPP

#include "egglib.hpp"
#include "FreqHolder.hpp"
#include "SiteHolder.hpp"

namespace egglib {

    class SiteHolder;

   /** \brief Compute Extended Haplotype Homozygosity statistics
    *
    * \ingroup diversity
    *
    * Compute statistics described in Sabeti et al. (Nature 2002),
    * Voight et al. (PLoS Biology 2006), Ramirez-Soriano et al.
    * (Genetics 2008) and Tang et al. (PLoS Biology 2007).
    *
    * The user must first load the core haplotype or site using the
    * set_core() method which also allows to specify option values, and
    * then all needed distant sites using load_distant(). Distant sites
    * must be loaded for one side only and with always increasing
    * distance relatively to the core. To load sites of the other side,
    * the user needs to call set_core() again with the same core site in
    * order to reset statistics. Statistics are automatically computed
    * and updated at each loaded distant site. It is required to load at
    * least one valid core site before using accessors.
    *
    * Header: <egglib-cpp/EHH.hpp>
    *
    */
    class EHH {

        private:

            // disabled methods
            EHH(const EHH& ehh) {}
            EHH& operator=(const EHH& ehh) {return * this;}

            // parameters
            double _par_EHH_thr;            // limit for IHH and dEHH
            double _par_EHHc_thr;           // limit for IHHc and dEHHc
            double _par_EHHS_thr;           // limit for IES
            double _par_EHHG_thr;           // limit for IES (genotypes)
            unsigned int _par_min_sam;      // minimum number of non-missing samples
            bool _opt_genotypes;            // whether to use genotypes
            bool _opt_crop_EHHS;            // whether to set EHHS values below threshold to 0 (both standard and genotypes)

            // variables
            SiteGeno _site_geno;            // object managing genotypes
            FreqHolder _frq;                // object used to analyse each site
            unsigned int _pl;               // ploidy (read from core site)
            unsigned int _ni;               // number of individuals (read from core site)
            unsigned int _nsam;             // total number of samples
            unsigned int _ncur_tot;         // current number of samples (excluding missing)
            unsigned int _num_tot;         // current number of samples (excluding missing)
            unsigned int _K_core;           // number of core haplotypes
            unsigned int _K_cur;            // current number of haplotypes
            unsigned int * _hap_core;       // core haplotype index for each sample // size: _sz_nsam
            unsigned int * _hap_cur;        // current haplotype index for each sample // size: _sz_nsam
            bool * _homoz_core;             // true if homozygote at core // size: _sz_cur
            bool * _homoz_cur;              // true if homozygote over [core-cur] // size: _sz_cur
            unsigned int * _num_core;       // number of non-missing samples per core haplotype, at core site // size: _sz_K_core
            unsigned int * _ncur_core;      // current number of samples per core haplotype // size: _sz_K_core
            unsigned int * _ncur_cur;       // current number of samples per cur haplotype // size: _sz_K_cur
            unsigned int * _hap_allele;     // allele index of each current haplotype (MISSING for missing) // size: _sz_K_cur
            unsigned int * _hap_origin;     // core haplotype index of each current haplotype // size: _K_cur
            unsigned int ** _branches;      // mapping of old-to-new haplotypes at breakpoints [old hap / allele_idx / new_hap] // size: _sz_branches * 3

            // cache
            unsigned int _sz_sam;
            unsigned int _sz_K_core;
            unsigned int _sz_K_cur;
            unsigned int _sz_branches;

            // statistics
            double _EHHS;                   // site EHHS
            double _EHHG;                   // site EHHS (genotypes)
            double _dEHHS;                  // EHHS decay
            double _dEHHG;                  // EHHS (genotypes) decay
            double _iES;                    // integrated EHHS
            double _iEG;                    // integrated EHHS (genotypes)
            bool _flag_dEHHS;               // true if EHHS reached decay
            bool _flag_dEHHG;               // true if EHHS (genotypes) reached decay
            double * _EHH;                  // EHH for each core haplotype
            double * _EHHc;                 // complementary EHH
            double * _IHH;                  // integrated EHH
            double * _IHHc;                 // integrated EHHC
            double * _dEHH;                 // EHH decay
            double * _dEHHc;                // EHHC decay
            unsigned int _num_dEHH;         // number of EHH reaching decay
            unsigned int _num_dEHHc;        // number of EHHc reaching decay
            double _dist;                   // previous distance

        public:

           /** \brief Constructor */
            EHH();

           /** \brief Destructor */
            virtual ~EHH();

           /** \brief Load the core site or region
            *
            * This method automatically resets the instance (clear all
            * previously computed data and reallocate arrays to proper
            * sizes). The Site instance passed as core is only used by
            * this method. All counters will be incremented, until the
            * next call to set_core(), or eventual destruction of
            * object. All thresholds are understood as either EHH or
            * EHHS values and therefore must lie between 0.0 and 1.0.
            *
            * \param site core site or region. If a region, haplotypes
            *        within the core region must have been identified
            *        previously and should be loaded as a Site instance.
            *        The site may contain missing data. The samples
            *        containing missing data at the core site will be
            *        ignored for all subsequently loaded distant site.
            * \param genotypes if true, consider that data are entered
            *        as unphased genotypes (the Site instance must have
            *        consistent data).
            * \param EHH_thr threshold EHH value.
            * \param EHHc_thr threshold EHHc value.
            * \param EHHS_thr threshold EHHS value.
            * \param EHHG_thr threshold EHHS (genotypes) value.
            * \param min_freq minimal absolute frequency for haplotypes
            *        (haplotypes with lower frequencies are ignored).
            *        Required to be strictly larger than zero.
            * \param min_sam minimal number of samples to continue computing
            *        (applied both within core haplotypes and for the total).
            * \param crop if True, set values of EHHS that are below the
            *        threshold to 0 to emulate the behaviour of the R
            *        package rehh (also affects iES).
            *
            */
            void set_core(const SiteHolder * site, bool genotypes,
                          double EHH_thr, double EHHc_thr,
                          double EHHS_thr, double EHHG_thr, unsigned int min_freq, unsigned int min_sam, bool crop);

           /** \brief Load a distant site
            *
            * For each core haplotype, compute or update all statistics.
            *
            * \param site the distant site to be loaded. The method will only
            *        throw an exception if the number of samples differ.
            * \param distance between the core and the distant site. The
            *        nature of the distance metrics is up to the user
            *        but must be consistent over sites. Distances must
            *        be >=0 this must correspond to the positions of
            *        distant sites relatively to the core region or
            *        site. Site must be loaded such as distances must
            *        always be increasing.
            *
            */
            void load_distant(const SiteHolder * site, double distance);

           /** \brief Number of used haplotypes of the core */
            unsigned int K_core() const;

           /** \brief Current number of haplotypes */
            unsigned int K_cur() const;

           /** \brief Current number of non-missing samples */
            unsigned int num_avail_tot() const;

           /** \brief Current number of non-missing samples for a core haplotype */
            unsigned int num_avail_core(unsigned int) const;

           /** \brief Current number of non-missing samples for a current haplotype */
            unsigned int num_avail_cur(unsigned int) const;

           /** \brief Get an EHHS value */
            double EHHS() const;

           /** \brief Get an EHHG value */
            double EHHG() const;

           /** \brief Get an EHH value */
            double EHHi(unsigned int haplotype) const;

           /** \brief Get an EHHc value */
            double EHHc(unsigned int haplotype) const;

           /** \brief Get an rEHH value */
            double rEHH(unsigned int haplotype) const;

           /** \brief Get an iES value */
            double iES() const;

           /** \brief Get an iEG value */
            double iEG() const;

           /** \brief Get an IHH value */
            double IHH(unsigned int haplotype) const;

           /** \brief Get an IHHc value */
            double IHHc(unsigned int haplotype) const;

           /** \brief Get an iHS value */
            double iHS(unsigned int haplotype) const;

           /** \brief Number of haplotypes for which computation of dEHH and IHH has been completed */
            unsigned int num_EHH_done() const;

           /** \brief Number of haplotypes for which computation of dEHHc and IHHc has been completed */
            unsigned int num_EHHc_done() const;

           /** \brief Tell if decay has been reached for EHHS */
            bool flag_EHHS_done() const;

           /** \brief Tell if decay has been reached for EHHG */
            bool flag_EHHG_done() const;

           /** \brief Get an EHH decay value */
            double dEHH(unsigned int haplotype) const;

           /** \brief Get an EHHc decay value */
            double dEHHc(unsigned int haplotype) const;

           /** \brief Get the maximum EHH decay value (computed on the fly) */
            double dEHH_max() const;

           /** \brief Get the average EHH decay value (on the fly) */
            double dEHH_mean() const;

           /** \brief Get an EHHS decay value */
            double dEHHS() const;

           /** \brief Get an EHHS (genotypes) decay value  */
            double dEHHG() const;
    };
}

#endif
