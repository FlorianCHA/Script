/*
    Copyright 2008-2018 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_DIVERSITY1_HPP
#define EGGLIB_DIVERSITY1_HPP

namespace egglib {

    class FreqHolder;
    class SiteDiversity;

   /** \brief Compute population summary statistics from allele frequencies at several sites
    *
    * \ingroup diversity
    *
    * Diversity1 instances cannot be copied. This class is designed to
    * allow reuse of objects without unnecessary memory reallocation.
    *
    * This class computes statistics that does not require access to a
    * full Site instance and for which only frequencies are needed. The
    * frequency for all sites that must be analyzed should be loaded.
    *
    * Statistics:
    *
    * code     | requirement                 | flag  | toggle flag
    * =========|=============================|=======|============
    * lt       | -                           | -     | -
    * ls       | -                           | -     | -
    * nsmax    | ls>0                        | 1     | -
    * S        | ls>0                        | 1     | -
    * Ss       | ls>0                        | 1     | -
    * eta      | ls>0                        | 1     | -
    * Pi       | ls>0                        | 1     | -
    * lso      | -                           | 4     | ori_site
    * nsmaxo   | lso>0                       | 8     | ori_site
    * So       | lso>0                       | 8     | ori_site
    * Sso      | lso>0                       | 8     | ori_site
    * etao     | lso>0                       | 8     | ori_site
    * lM       | lso>0                       | 8     | ori_site
    * pM       | lM>0                        | 16    | ori_site
    * nseffo   | lso>0                       | 32    | ori_div
    * thetaH   | lso>0                       | 32    | ori_div
    * thetaL   | lso>0                       | 32    | ori_div
    * Hns      | lso>0                       | 32    | ori_div
    * Hsd      | So>0 & nseffo>=3 & varZ>0   | 1024  | ori_div
    * E        | So>0 & nseffo>=3 & varE>0   | 2048  | ori_div
    * Dfl      | So>0 & nseffo>=3 & varDfl>0 | 4096  | ori_div
    * F        | So>0 & nseffo>=3 & varF>0   | 8192  | ori_div
    * nseff    | ls>0                        | 128   | basic
    * thetaW   | ls>0                        | 128   | basic
    * Dxy      | ls>0 npop=2                 | 16384 | basic
    * Da       | ls>0 npop=2                 | 16384 | basic
    * Fstar    | S>0 & ns>2                  | 256   | basic
    * D        | S>0 & ns>3                  | 512   | basic
    * Deta     | S>0 & ns>3                  | 512   | basic
    * Dstar    | S>0 & ns>3                  | 512   | basic
    * sites    | i<S                         | -     | site_lists
    * sites_o  | i<So                        | -     | site_lists
    * singl    | i<Ss                        | -     | site_lists
    * singl_o  | i<Sso                       | -     | site_lists
    *
    * Header: <egglib-cpp/Diversity.hpp>
    *
    */
    class Diversity1 {

        private:

            void _init();
            void _free();

            unsigned int _lt;
            unsigned int _ls;
            unsigned int _lso;

            double _nseff;
            double _nseffo;
            unsigned int _nsmax;
            unsigned int _nsmaxo;

            unsigned int _S;
            unsigned int _Ss;
            unsigned int _So;
            unsigned int _Sso;
            unsigned int _eta;
            unsigned int _etao;
            double _Pi;
            double _PiForF;
            double _Pi0;
            double _Pi1;
            double _thetaW;
            double _Dxy;
            double _Da;
            double _k;
            double _ko;
            unsigned int _nsingl;
            unsigned int _nsingld;

            double _D;
            double _Deta;
            double _Dstar;
            double _Fstar;

            unsigned int _option_ns_set;
            unsigned int * _esse;
            unsigned int _esse_c;

            double _thetaPi;
            double _thetaH;
            double _thetaL;
            double _Hns;
            double _Hsd;
            double _E;
            double _Dfl;
            double _F;
            double _pM;
            unsigned int _nM;

            void _basic();
            void _oriented();

            unsigned int _flag;
            bool _flag_ori_site;
            bool _flag_ori_div;
            bool _flag_basic;
            bool _flag_site_lists;
            bool _option_multiple;

            unsigned int * _sites;
            unsigned int * _sites_o;
            unsigned int * _singl;
            unsigned int * _singl_o;
            unsigned int _c_sites;
            unsigned int _c_sites_o;
            unsigned int _c_singl;
            unsigned int _c_singl_o;

            Diversity1(const Diversity1& src) {}
            Diversity1& operator=(const Diversity1& src) { return *this; }

        public:

            Diversity1(); ///< Constructor
            ~Diversity1(); ///< Destructor

            void load(const FreqHolder& freqs, const SiteDiversity& div, unsigned int position); ///< Analyze a site
            unsigned int compute(); ///< Compute statistics, return flag but does not reset
            void reset_stats(); ///< Reset counters to 0
            void toggle_off(); ///< Cancel all flags
            void toggle_ori_site(); ///< Activate per-site oriented
            void toggle_ori_div(); ///< Activate per-gene oriented
            void toggle_basic(); ///< Activate basic per-gene
            void toggle_site_lists(); ///< Activate lists of site positions
            void set_option_multiple(bool); ///< Set multiple option (default: False)
            void set_option_ns_set(unsigned int); ///< Set maximum number of samples, for H and co. (default: UNKNOWN)

            unsigned int lt() const; ///< Number of loaded sites (total)
            unsigned int ls() const; ///< Number of loaded sites (with >=2 valid data)
            unsigned int lso() const;  ///< Number of loaded orientable sites (with valid data)
            unsigned int S() const; ///< Number of polymorphic sites
            unsigned int Ss() const; ///< Number of polymorphic sites with =1 singleton
            unsigned int So() const; ///< Number of polymorphic orientable sites
            unsigned int Sso() const; ///< Number of polymorphic orientable sites with =1 singleton
            unsigned int nsingld() const; ///< Number derived singletons
            unsigned int eta() const; ///< eta
            unsigned int etao() const; ///< eta for orientable sites
            double D() const; ///< Tajima's D
            double Deta() const;  ///< Tajima's D using eta instead of S
            double Dstar() const; ///< Fu and Li's D*
            double Fstar() const; ///< Fu and Li's F*
            double thetaW() const; ///< Theta estimator based on S
            double Pi() const;  ///< Sum of He
            double Dxy() const; ///< Pairwise distance for 1st pair
            double Da() const; ///< Net pairwise distance for 1st pair
            double nseff() const; ///< Average number of used samples
            unsigned int nsmax() const; ///< Largest number of used samples
            double nseffo() const; ///< Average number of used samples for orientable sites
            unsigned int nsmaxo() const; ///< Largest number of used samples for orientable sites
            double thetaPi() const; ///< thetaPi estimator (using orientable sites)
            double thetaH() const; ///< ThetaH estimator
            double thetaL() const; ///< ThetaL estimator
            double Hns() const; ///< Unstandardized Fay and Wu's H
            double Hsd() const; ///< Fay and Wu's H standardized by Zeng et al.
            double E() const; ///< Zeng et al.'s E
            double Dfl() const; ///< Fu and Li's D
            double F() const; ///< Fu and Li's F
            double pM() const; ///< Li's MFDM test p value (large positive value by default)
            unsigned int nM() const; ///< Sites available for MFDM test
            unsigned int site(unsigned int) const; ///< Get position of polymorphic site
            unsigned int site_o(unsigned int) const; ///< Get position of polymorphic orientable site
            unsigned int singl(unsigned int) const; ///< Get position of site with a singleton
            unsigned int singl_o(unsigned int) const; ///< Get position of site with an orientable singleton
    };
}

#endif
