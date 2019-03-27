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

#ifndef EGGLIB_GENETICCODE_HPP
#define EGGLIB_GENETICCODE_HPP

namespace egglib {

    class SiteHolder;

   /** \brief Hold genetic code tables
    *
    * \ingroup core
    *
    * Handle genetic codes. All genetic codes defined by the National
    * Centor for Biotechnology Information are supported.
    *
    * Header: <egglib-cpp/GeneticCode.hpp>
    *
    */
    class GeneticCode {

        public:

           /** \brief Constructor
            *
            * Build an instance of the genetic code.
            *
            */
            GeneticCode(unsigned int index);

           /** \brief Default constructor
            *
            * Like GeneticCode(int), except that the code used is 1
            * (standard).
            *
            */
            GeneticCode();

           /** \brief Get the current genetic code
            *
            */
            unsigned int get_code() const;

           /** \brief Set the genetic code
            *
            */
            void set_code(unsigned int index);

           /** \brief Returns the translation of a codon
            *
            * The codon should be passed as an integer code (see
            * codon()). This methods returns the single amino acid code
            * for valid codons (represented by integers in the range
            * 0-63), and 'X' for any other integer.
            *
            */
            char aminoacid(unsigned int codon) const;

           /** \brief Tells if a codon is an initiation codon
            *
            * The codon should be passed as an integer code (see
            * codon()). This methods returns `true` if the codon in an
            * initiation codon (including any of the alternative
            * initiation codons known for the genetic code of the set
            * for the current object), and `false` otherwise (including
            * for invalid codon codes).
            */
            bool start(unsigned int codon) const;

           /** \brief Give the number of non-synonymous sites of a codon
            *
            * The number is in the range 0-3 (3 is all changes at all of
            * the three positions would lead to a non-synonymous
            * change).
            *
            * \param codon codon integer code (see codon()).
            *
            * \param ignorestop if true, potential changes to stop codons
            * are excluded and all stop codons return 0; if false, changes
            * to stop codons are considered to be non-synonymous.
            *
            */
            double NSsites(unsigned int codon, bool ignorestop = true) const;

           /** \brief Give the number of synonymous sites of a codon
            *
            * The number is in the range 0-3 (3 is all changes at all of
            * the three positions would lead to a non-synonymous
            * change).
            *
            * \param codon codon integer code (see codon()).
            *
            * \param ignorestop if true, potential changes to stop codons
            * are excluded and all stop codons return 0; if false, changes
            * to stop codons are considered to be non-synonymous.
            *
            */
            double Ssites(unsigned int codon, bool ignorestop = true) const;

           /** \brief Give the number of non-synonymous sites of a codon site
            *
            * The number is in the range 0-3 (3 is all changes at all of
            * the three positions would lead to a non-synonymous
            * change). This is the same as NSistes(unsigned int, bool),
            * but average over all samples based on provided Site
            * instances.
            *
            * \param site1 first position of the codon site.
            * \param site2 second position of the codon site.
            * \param site3 third position of the codon site.
            * \param num_samples variable used to provide the number of
            * samples analyzed by the method (that is, number of samples
            * minus number of samples containing at least one missing
            * data). The original value of the variable is ignored and
            * is modified by the instance. If 0, the return value should
            * be ignored.
            * \param ignorestop if true, potential changes to stop codons
            * are excluded and all stop codons are treated as missing
            * data; if false, changes to stop codons are considered to
            * be non-synonymous.
            * \param A integer value representing bases A.
            * \param C integer value representing bases C.
            * \param G integer value representing bases G.
            * \param T integer value representing bases T.
            *
            * All three codon positions must have the same number of
            * samples such as the ith nucleotides at the three sites
            * give the codon for the ith sample. (Same ploidy as well.)
            *
            * All nucleotides represented by an integer allele values
            * not matching either of the four values passed as the A, C,
            * G and T arguments are considered as missing data.
            *
            * \warning This method will not accept coding sequences
            * mixing upper and lower case characters. It is however
            * possible to configure how the four nucleotides are
            * represented.
            *
            */
            double NSsites(const SiteHolder& site1, const SiteHolder& site2, const SiteHolder& site3, unsigned int& num_samples, bool ignorestop = true, int A='A', int C='C', int G='G', int T='T') const;

            /// Like NSsites, for synonmous sites
            double Ssites(const SiteHolder& site1, const SiteHolder& site2, const SiteHolder& site3, unsigned int& num_samples, bool ignorestop = true, int A='A', int C='C', int G='G', int T='T') const;

           /** \brief Give the number of non-synonymous sites of a codon site
            *
            * See NSsites(const SiteHolder&, const SiteHolder&, const SiteHolder&, unsigned int&, bool, int, int C, int, int).
            * Do the same, except that this version takes a single SiteHolder
            * reference instance of three, and the SiteHolder reference passed
            * to this function contains integer alleles representing
            * codons, that is in the range 0-63. Other values are
            * considered to be missing data.
            *
            */
            double NSsites(const SiteHolder& codons, unsigned int& num_samples, bool ignorestop = true, int A='A', int C='C', int G='G', int T='T') const;

            /// Like NSsites, but for synonymous sites
            double Ssites(const SiteHolder& codons, unsigned int& num_samples, bool ignorestop = true, int A='A', int C='C', int G='G', int T='T') const;

           /** \brief Return the integer code for a codon
            *
            * The first, second and third bases of the codon must be
            * passed as character arguments. The case of characters is
            * ignored. Returns and integer in the range [0, 63] for the
            * 64 codons (see the table below).
            *
            * \warning The base 'U', although biologically relevant, is
            * treated as an invalid base.
            *
            * \note As a static method, this method can be called as
            * `GeneticCode::codon(base1, base2, base3)` directly (it
            * does not require instanciation of an object) and it is not
            * dependent on any genetic code specification.
            *
            * The codons are identified by single integers as given by
            * the table below:
            *
            * Integer code | First base | Second base | Third base
            * :----------: | :--------: | :---------: | :--------:
            * 0            | T          | T           | T
            * 1            | T          | T           | C
            * 2            | T          | T           | A
            * 3            | T          | T           | G
            * 4            | T          | C           | T
            * 5            | T          | C           | C
            * 6            | T          | C           | A
            * 7            | T          | C           | G
            * 8            | T          | A           | T
            * 9            | T          | A           | C
            * 10           | T          | A           | A
            * 11           | T          | A           | G
            * 12           | T          | G           | T
            * 13           | T          | G           | C
            * 14           | T          | G           | A
            * 15           | T          | G           | G
            * 16           | C          | T           | T
            * 17           | C          | T           | C
            * 18           | C          | T           | A
            * 19           | C          | T           | G
            * 20           | C          | C           | T
            * 21           | C          | C           | C
            * 22           | C          | C           | A
            * 23           | C          | C           | G
            * 24           | C          | A           | T
            * 25           | C          | A           | C
            * 26           | C          | A           | A
            * 27           | C          | A           | G
            * 28           | C          | G           | T
            * 29           | C          | G           | C
            * 30           | C          | G           | A
            * 31           | C          | G           | G
            * 32           | A          | T           | T
            * 33           | A          | T           | C
            * 34           | A          | T           | A
            * 35           | A          | T           | G
            * 36           | A          | C           | T
            * 37           | A          | C           | C
            * 38           | A          | C           | A
            * 39           | A          | C           | G
            * 40           | A          | A           | T
            * 41           | A          | A           | C
            * 42           | A          | A           | A
            * 43           | A          | A           | G
            * 44           | A          | G           | T
            * 45           | A          | G           | C
            * 46           | A          | G           | A
            * 47           | A          | G           | G
            * 48           | G          | T           | T
            * 49           | G          | T           | C
            * 50           | G          | T           | A
            * 51           | G          | T           | G
            * 52           | G          | C           | T
            * 53           | G          | C           | C
            * 54           | G          | C           | A
            * 55           | G          | C           | G
            * 56           | G          | A           | T
            * 57           | G          | A           | C
            * 58           | G          | A           | A
            * 59           | G          | A           | G
            * 60           | G          | G           | T
            * 61           | G          | G           | C
            * 62           | G          | G           | A
            * 63           | G          | G           | G
            *
            * All other triplets: egglib.UNKNOWN.
            *
            */
            static unsigned int codon(char first, char second, char third);

           /** \brief Return the integer code for a codon
            *
            * The first, second and third bases of the codon must be
            * passed as integer code, according to the following
            * mapping: 0 for T, 1 for C, 2 for G and 3 for T. This code
            * must absolutely be following and no other value may be
            * passed. Returns and integer in the range [0, 63] for the
            * 64 codons (see the documentation of the method
            * codon(char, char, char)).
            *
            * \note As a static method, this method can be called as
            * `GeneticCode::int2codon(base1, base2, base3)` directly (it
            * does not require instanciation of an object) and it is not
            * dependent on any genetic code specification.
            *
            */
            static inline unsigned int int2codon(unsigned int base1, unsigned int base2, unsigned int base3) {
                return (base1<<4|base2<<2|base3);
            }

           /** \brief Returns one of the base of a codon
            *
            * This method can be called on the class directly (as in
            * `GeneticCode::base(0, 0)` and it is not dependent on the
            * specification of a genetic code).
            *
            * \param codon the codon should be passed as an integer code
            * (see codon()). Only values in the range 0-63 are
            * supported.
            * \param index index of the base to extract (only 0, 1 and 2
            * are accepted; other values will result in aberrant
            * outcome).
            * \return Returns the character at the specified position of
            * the codon (as an upper-case character).
            *
            * \warning The methods returns '?' if it cannot perform
            * base extraction, but it is not guaranteed that all invalid
            * arguments will be detected properly.
            *
            */
            static inline char base(unsigned int codon, unsigned int index) {

               /*
                * index is       0, 1 or 2
                * 2-index is     2, 1 or 0 respectively
                * 2*(2-index) is 4, 2 or 0 respectively
                * 3 is (in binary) 000011
                * 3 << (2-index) is 48, 12 or 3, that is, in binary,
                * 110000, 001100 or 000011 and is used as mask
                *
                * codon & (3 << (2-index)) extracts bits for an index
                * and >> (2*(2-index)) shifts the bits to unit level
                *
                * in the end we have 00, 01, 10 or 11 (in binary), that
                * is 0, 1, 2 or 4 respectively
                *
                */

                switch ((codon & (3 << (2*(2-index)))) >> (2*(2-index))) {
                    case 0: return 'T';
                    case 1: return 'C';
                    case 2: return 'A';
                    case 3: return 'G';
                    default: return '?';
                }
            }

           /** \brief Returns the number of nucleotide differences between two codons
            *
            * This method can be called on the class directly (as in
            * `GeneticCode::base(0, 0)` and it is not dependent on the
            * specification of a genetic code).
            *
            * This method is only valid if both arguments are less than
            * 64. Returns only 0, 1, 2 or 3.
            *
            */
            static inline unsigned int ndiff(unsigned int codon1, unsigned int codon2) {
                return diff1(codon1, codon2) + diff2(codon1, codon2) + diff3(codon1, codon2);
            }

           /** \brief Check if the first position of two codons is identical
            *
            */
            static inline bool diff1(unsigned int codon1, unsigned int codon2) {
                return ((codon1 & 48) >> 4) != ((codon2 & 48) >> 4);
            }

           /** \brief Check if the second position of two codons is identical
            *
            */
            static inline bool diff2(unsigned int codon1, unsigned int codon2) {
                return ((codon1 & 12) >> 2) != ((codon2 & 12) >> 2);
            }

           /** \brief Check if the third position of two codons is identical
            *
            */
            static inline bool diff3(unsigned int codon1, unsigned int codon2) {
                return (codon1 & 3) != (codon2 & 3);
            }

           /** \brief Translate a codon directly
            *
            * \param first first codon position.
            * \param second second codon position.
            * \param third third codon position.
            * \param smart smart translation.
            * \return ASCII-coded aminoacid.
            *
            * Codon positions should be ASCII-coded. Return 'X' if
            * missing data or invalid nucleotides. For fourfold
            * degenerate positions. Codons including non-ambiguity
            * characters always return 'X' (even '?' at a fourfold
            * degenerate position), except if the codon is '---' (in
            * that case, '-' is returned).
            *
            */
            int translate(int first, int second, int third, bool smart);

            /// Get the name of the genetic code
            const char * name() const;

            /// Get the number of available codes
            static unsigned int num_codes();

        private:

           /** \brief Copy constructor not available */
            GeneticCode(GeneticCode& src) {}

           /** \brief Copy assignment operator not available */
            GeneticCode& operator=(GeneticCode& src) { return *this; }

            void _dev_implied(char ch, unsigned int& num, int dest[]);
            int _implied_1st[4];
            int _implied_2nd[4];
            int _implied_3rd[4];
            unsigned int _num_implied_1st;
            unsigned int _num_implied_2nd;
            unsigned int _num_implied_3rd;

            static const char _aa[];
            static const char _start[];
            static const double _NS1[]; // stop codons considered
            static const double _NS2[]; // ignorestop
            static const double _S1[];
            static const double _S2[];
            static const char * _names[];
            static const unsigned int _codes[];
            unsigned int _code;
            unsigned int _index;
            unsigned int _shift;
    };
}

#endif
