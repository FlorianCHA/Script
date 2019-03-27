/*
    Copyright 2016-2018 Stéphane De Mita, Mathieu Siol, Thomas Coudoux

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

#ifndef EGGLIB_WINDOW_HPP
#define EGGLIB_WINDOW_HPP

using namespace std;
#include "egglib.hpp"
#include "SiteHolder.hpp"
#include "VCF.hpp"
#include "BED.hpp"


// remove this somehow
#include <vector>


namespace egglib {

    class WPool;

    ///< Class for sites within a double-linked list
    class WSite {
        private:
            WSite * _next;
            WSite * _prev;
            unsigned int _pos;
            SiteHolder _site;
            WPool * _pool;

            WSite() {}
            WSite(const WSite& src) {}
            WSite& operator=(const WSite* src) { return *this; }

        public:
            WSite(WPool *); ///< Constructor
            ~WSite(); ///< Destructor
            SiteHolder& site(); ///< Get the included site object
            void set_pos(unsigned int); ///< Set site position
            unsigned int get_pos() const; ///< Get site position
            WSite * next(); ///< Get next site
            WSite * prev(); ///< Get previous site
            WSite * push_back(WSite *); ///< Add site to the end, return new end
            WSite * pop_front(); ///< Disconnect first, return it to pool, return its follower
            WSite * pop_back(); ///< Disconnect last, return it to pool, return its predecessor
            void reset(unsigned int pl); ///< Reset values
            void init(); ///< Set pointers to NULL
    };

    ///< WSite pool
    class WPool {
        private:
            WSite ** _cache;
            WSite ** _pool;
            unsigned int _c_cache;
            unsigned int _c_pool;
            unsigned int _n_pool;
            WPool(const WPool& src) {}
            WPool& operator=(const WPool& src) { return *this; }

        public:
            WPool();
            ~WPool();
            WSite * get();
            void put(WSite *);
    };

    ///< Class for sliding window
    class VcfWindow {
        private:
            WPool _pool;
            WSite * _first_site;            // first site of current window
            WSite * _last_site;             // last site of current window
            unsigned int _num;              // number of sites

            VcfParser * _vcf;
            unsigned int _wsize;       // window size
            unsigned int _wstep;       // window step
            bool _unit_bp;             // true if wsize/wstep expressed in genomic bp
            unsigned int _max_missing; // SiteHolder argument

            unsigned int _start_pos;   // where to start the sliding window
            unsigned int _stop_pos;    // where to stop the sliding window (not included)
            unsigned int _cur_pos;     // current position (if unit_bp)

            unsigned int _first_pos;   // start of current window
            unsigned int _last_pos;    // end of current window (included)
            bool _good;                // there will be a new window

            char * _chrom;
            unsigned int _c_chrom;

            bool _read(unsigned int);  // read one site
                                       // return false if no site available,
                                       //   specified position or next chromosome is reached
                                       // if return false, nothing is read and set _good to false
            void _add();

            VcfWindow(const VcfWindow& src) {}
            VcfWindow& operator=(const VcfWindow& src) { return *this; }

        public:
            VcfWindow(); ///< Constructor
            ~VcfWindow(); ///< Destructor
            void setup(VcfParser&, unsigned int wsize, unsigned int wstep, bool unit_bp,
                    unsigned int start_pos, unsigned int stop_pos,
                    unsigned int max_missing); ///< Setup a new sliding window
            void next_window(); ///< Load the next window
            const char * chromosome() const; ///< Get chromosome
            unsigned int num_samples() const; ///< Number of samples
            unsigned int num_sites() const; ///< Number of actual sites (from vcf)
            unsigned int first_pos() const; ///< Window first position (UNKNOWN if unit is not bp and no site at all)
            unsigned int last_pos() const; ///< Window last position (included) (UNKNOWN if unit is not bp and no site at all)
            const WSite * first_site() const; ///< Get first site (NULL if no site at all)
            const WSite * last_site() const; ///< Get last site (NULL if no site at all)
            unsigned int first_site_pos() const; ///< Position of the first site
            unsigned int last_site_pos() const; ///< Position of the last site
            bool good() const; ///< False if sliding has completed
            const WSite * get_site(unsigned int) const; ///< Get a random site (slower -- there must be enough sites)
};





    class BedWindow {

        private:
            VcfParser * _VCF;
            BedParser * _Bed;
            vector<egglib::SiteHolder*> Site_list;
            char * _chrom_window;
            unsigned int start_w;
            unsigned int end_w;
            unsigned int _ploidy_pv;
            unsigned int _end_bound;
            unsigned int _indice;

            /* \brief edits the chromosome of the sliding window
            *   \param chrom: chromosome.
            */
            void set_chromosome_window(const char * chrom);

            /* \brief edits the start chromosomal position of the sliding window
            *   \param position:start chromosomal position.     
            */
            void set_start_window(unsigned int position);

            /* \brief edits the end chromosomal position of the sliding window
            *   \param position:end chromosomal position.     
            */
            void set_end_window(unsigned int position);


            /* \brief gets the start chromosomal position of the current sliding window.
            *
            */
            void get_start_variant_vcf(unsigned int indice);

            /* \brief checks if the sliding windows has reached the chromosome end position or the last variant of 
                        *         the chromosome read by the sliding window.
            */
            bool before_bound() const;


            /* \brief this method, allows to read and store variants read by the sliding_window.
            *          The method vcf_process of the class SiteHolder allows to extract alleclic data 
            *          and  compute frequencies from each variant read by the sliding window in an 
                        *          object of the class SiteHolder. Under this form, each variants read are stored in a list. 
                        *    Beware: This list is emptied after each calls of the method "get_sliding_window" and fills
            *            by new variants.
            */
            void get_sliding_window();


        public:

            unsigned int b_size;
            unsigned int _start_pv;
            unsigned int _end_pv;
            unsigned int _miss_pv;


            BedWindow(); //<-- constructor

            ~BedWindow(); //<--

            /* \brief get the list of SiteHolder object extract by the sliding window.
            *
            */
            std::vector<egglib::SiteHolder*> get_site_list();

            /* \brief set all variables used to configure the bed sliding window 
            *
            *   \param VCF a VcfParser reference containing data and having 
            *          the GT format field filled 
            *   \param Bed a BedParser reference containing coordinates chromosomal
            *          according to the VcfParser loaded.
            *   \param start if a subset of samples must be considered,
            *          index of the first sample to consider (by default,
            *          all samples are considered).
            *   \param stop if a subset of samples must be considered,
            *          index of the last sample to consider (by default,
            *          all samples are considered).
            *   \param max_missing maximum number of missing alleles.
            *          If this proportion is processing is stopped and
            *          get_missing() returns max_missing + 1. Only missing
            *          data in this data set are considered.
            */
            void configure(egglib::VcfParser * VCF, BedParser * Bed, unsigned int start_pv, unsigned int end_pv, unsigned int miss_pv);


            /* \brief This method allows emptying the list of site, filled by the sliding window. it allows to 
            *         dellocate the memory used to create sites.
            */
            void clear_window();

            /* \brief This function allows the bed sliding window to position at the next chromosome coordinate 
            *         of the ".bed" file linked.
            *
            */
            void next();

            /* \brief resets all variables used to configure the bed sliding window
            *
            */
            void init();

            /* \brief checks if the sliding window can continue to grow on the file vcf according the number of 
            *   of data in the bed object loaded.          
            */
            bool good();

            /* \brief gets the sliding window chromosomal starting position
            *
            */
            unsigned int get_start_w();

            /* \brief gets the sliding window chromosomal ending position
            *
            */
            unsigned int get_end_w();

            /* \brief gets the sliding window chromosome
            *
            */
            char * get_chrom_w();

            /* \brief gets the sliding window chromosome
            *   \param i: indice of the site searched. 
            *    \beware : The indice must be smaller than the total number of sites
            *              contained in the current sliding window. 
            */
            egglib::SiteHolder& get_site(unsigned int i);

            /* \brief gets the number of site saved in the current sliding window
            */
            int unsigned get_n_site();


            void at(unsigned int indice);

    };

}

#endif

