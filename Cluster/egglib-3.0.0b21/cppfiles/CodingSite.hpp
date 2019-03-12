/*
    Copyright 2013,2016 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_CODINGSITE_HPP
#define EGGLIB_CODINGSITE_HPP

#include "egglib.hpp"

namespace egglib {

    class GeneticCode;
    class SiteHolder;

   /** \brief Holds data for a coding site (as a triplet of sites)
    *
    * This class assists the detection of polymorphism at codon sites,
    * although diversity analyses themselves have to be performed using
    * SiteDiversity itself. The class perform analyses through the
    * process() method, which automatically resets all previously
    * stored data.
    *
    * Header: <egglib-cpp/CodongSite.hpp>
    *
    */
    class CodingSite {

        private:
            CodingSite(const CodingSite& src) {}
            CodingSite& operator=(const CodingSite& src) {return *this;}
            void _init();
            void _free();
            unsigned int _check_allele(int);
            unsigned int _ni;
            unsigned int _no;
            unsigned int _pl;
            unsigned int _nseff;
            unsigned int _nseffo;
            double _NSsites;
            double _Ssites;
            unsigned int _nstop;
            int ** flags; // for each pairwise value:
                                // bit 0: non synonymous
                                // bit 1: first position changed
                                // bit 2: second position changed
                                // bit 3: third position changed
                                // bits 4 and 5: number of changes
            unsigned int _nall;
            unsigned int _nall_c;
            SiteHolder * _aminoacids;
            SiteHolder * _codons;
            char _amino;
            unsigned int _codon1;

            bool _process_helper(const GeneticCode&, bool, int, int, int);
            void _process_helper2(const GeneticCode&, unsigned int);

        public:
            CodingSite(); ///< Constructor
            ~CodingSite(); ///< Destructor

           /** \brief Analyzes a codon site
            *
            * The three codon positions must be loaded as Site instances
            * containing nucleotides encoded as integer values. All
            * values except values equal to A, C, G and T
            * (case-independent) are treated as missing data. Obviously,
            * the three sites must have the same number of samples and
            * also the same number of populations (and matching
            * affectation of samples to populations). Upon processing,
            * the class generates and holds a Site instance (available
            * as codons()) containing data from the three sites merged,
            * and another (available as aminoacids()) with the same data
            * translated).
            *
            *   \param site1 first nucleotide position of the codon.
            *   \param site2 second nucleotide position of the codon.
            *   \param site3 third nucleotide position of the codon.
            *   \param code GeneticCode instance representing the code
            *   to be used for treating this codon.
            *   \param skipstop if true, stop codons are treated as
            *   missing data and skipped. If set to true, potential
            *   mutations to stop codons are not taken into account when
            *   estimating the number of non-synonymous sites. Warning
            *   (this may be counter-intuitive). It actually assumes
            *   that stop codons are not biologically plausible and
            *   considers them as missing data. On the other hand, if
            *   skipstop is false, it takes stop codons as if they were
            *   valid amino acids.
            *   \param max_missing maximum number of missing data to
            *   allow (including stop codons if skipstop if true).
            *   \return A boolean indicating whether analysis was
            *   completed. If False, data contained in the object should
            *   not be used since stored objects will not have been
            *   filled. Even ni() and nieff() will be invalid.
            *
            * \note It is not allowed to use egglib::UNKNOWN for any of
            * A, C, G and T argument.
            *
            */
            bool process(const SiteHolder& site1, const SiteHolder& site2,
                const SiteHolder& site3, const GeneticCode& code,
                bool skipstop=true, unsigned int max_missing=MAX);

           /** \brief Get access to merged codon data
            *
            * Requires that a codon site has been analyzed using
            * process().
            *
            * This instance contains codon alleles. See GeneticCode for
            * information about encoding of codons.
            *
            */
            const SiteHolder& codons() const;

           /** \brief Get access to amino acid data
            *
            * Requires that a codon site has been analyzed using
            * process().
            *
            * This instance contains amino acids. It contains integer
            * data representing amino acids ('*' for stop codons). The
            * GeneticCode instance passed to process() determines the
            * translation.
            *
            */
            const SiteHolder& aminoacids() const;

            unsigned int ni() const; ///< Number of ingroup indiv
            unsigned int no() const; ///< Number of outgroup indiv
            unsigned int pl() const; ///< Ploidy

           /** \brief Number of analyzed samples
            *
            * Requires that a codon site has been analyzed using
            * process().
            *
            * Value bound by 0 and ns(). Depends on the number of
            * missing data (and stop codons, if skipstop was set to
            * true).
            *
            */
            unsigned int nseff() const;

            unsigned int nseffo() const; ///< Analyzed samples for outgroup

           /** \brief Estimated number of nonsynonymous sites
            *
            * Requires that a codon site has been analyzed using
            * process().
            *
            */
            double NSsites() const;

           /** \brief Estimated number of synonymous sites
            *
            * Requires that a codon site has been analyzed using
            * process().
            *
            */
            double Ssites() const;

           /** \brief Number of stop codons met during processing of codon site
            *
            * Requires that a codon site has been analyzed using
            * process().
            *
            * Not affected by the skipstop option.
            *
            */
            unsigned int nstop() const;

           /** \brief Number of nucleotides differences between two given codon alleles
            *
            * Give, for two given codon alleles (see alleles()), the
            * number of nucleotide differences between them.  Possible
            * values are 1, 2 or 3.
            *
            * Requires that a codon site has been analyzed using
            * process().
            *
            * The indexes must both be < nall() but they may be passed
            * in any order.
            *
            */
            unsigned int ndiff(unsigned int codon_allele1, unsigned int codon_allele2) const;

           /** \brief Check if two given codon alleles differ at a given position
            *
            * Tells, for two given codon alleles (see alleles()), if
            * they differ (at least) at the specified position. Possible
            * values are 0, 1 or 2.
            *
            * Requires that a codon site has been analyzed using
            * process().
            *
            * The indexes must both be < nall() but they may be passed
            * in any order.
            *
            */
            bool mutated(unsigned int codon_allele1, unsigned int codon_allele2, unsigned int pos) const;

           /** \brief True if the two given codon alleles encode different aminoacids
            *
            * Requires that a codon site has been analyzed using
            * process().
            *
            * The indexes must both be < nall() but they may be passed
            * in any order.
            *
            */
            bool NS(unsigned int codon_allele1, unsigned int codon_allele2) const;
    };
}

#endif
