/*
    Copyright 2014 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_EXPORT_HPP
#define EGGLIB_EXPORT_HPP

#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <limits>
#include "egglib.hpp"
#include "DataHolder.hpp"
#include "Tree.hpp"

namespace egglib {

   /** \brief Abstract Base class for exporters
    *
    */
    class BaseFormatter {

        public:

           /** \brief Constructor
            *
            * Sets to console output.
            *
            */
            BaseFormatter();

           /** \brief Destructor
            *
            */
            virtual ~BaseFormatter();

           /** \brief Open a file for writing
            *
            * All subsequent formatting operations will be done using
            * this file as output, until close() is called or open_file()
            * again (to create a new file). By default (if open_file() is
            * never called), output goes to the standard output.
            *
            * If an output file was already open, it is closed prior
            * opening the new one (it is not necessary to call close()).
            *
            * \return A boolean indicating whether the file has been
            * opened sucessfully.
            * 
            */
            bool open_file(const char * fname);

           /** \brief Write to a string buffer
            *
            * All subsequent formatting operations will be done into
            * an internally-stored string buffer which can be accessed
            * using string(). This closes the file if one is open. This
            * clears the string buffer if one has already been opened.
            *
            */
            void to_str();

           /** \brief Access the internal string buffer
            *
            * Returns the current string stored in the output string
            * buffer. This does not reset the output buffer (one must
            * call to_str() again to initialize a new output string.
            *
            */
            std::string get_str();

           /** \brief Restore the default (export to standard output)
            *
            * This closes the file if one is open.
            *
            */
            void to_cout();

           /** \brief Write a line as is
            *
            * A newline is automatically if the argument ``eol`` is
            * ``true``. It is legal to pass an empty string.
            *
            *
            */
            void write(const char * bit, bool eol);

           /** \brief Flush the current output
            *
            */
            void flush();

           /** \brief Close the output file
            *
            * If there is no open file, nothing is done. Otherwise the
            * current output file is closed and any subsequent output
            * will be directed to the standard output. It is not
            * required to call this method between to consecutive calls
            * to open_file() (in order to change the output file). The file
            * is also closed properly when this object is destroyed.
            *
            */
            void close();

        protected:

            const char * _fname;
            std::ostream * _stream;
            std::ofstream _fstream;
            std::ostringstream _sstream;
            std::ostream * _cache_stream;
            bool _is_file;

    };

   /** \brief Holder class for ms-type and newick tree formatting methods.
    *
    * \ingroup parsers
    *
    * This class cannot be built. Only method can be called.
    *
    * Header: <egglib-cpp/Export.hpp>
    *
    */
    class Export : public BaseFormatter {

        public:

           /** \brief Constructor. */
            Export();

           /** \brief Destructor. */
            ~Export();

           /** \brief Write a newick-formatted tree.
            *
            * Write the data as a single line, complete with its newline
            * character. It is required that leaves are properly and
            * hierarchically connected up to the root; non-network
            * structure).
            *
            * \param tree a completed genealogical tree with labelled
            * leaves.
            *
            * \param blen whether to export the value of branch lengths.
            *
            * \param eol whether to print a newline character after the
            * tree.
            *
            */
            void newick(const Tree& tree, bool blen, bool eol);

           /** \brief Specifies the number of positions for ms exporting
            *
            * This method should be called before calling ms() in order
            * to specify the number of sites. This value must match the
            * number of sites in the DataHolder that will be passed to
            * ms(), and all positions must be specified with
            * set_positions() after call to this method.
            *
            */
            void ms_num_positions(unsigned int n);

           /** \brief Specifies the position of a site for ms exporting
            *
            * This method should be called before calling ms() in order
            * to specify the position of each site. The number of sites
            * must have been specified previously using
            * ms_num_positions(). The position must be <0 and >1.
            *
            */
            void ms_position(unsigned int site, double position);

           /** \brief Assign default positions to all sites
            *
            * This method should be called before calling ms() in order
            * to specify the number of sites, if the position are not
            * defined or irrelevant. The argument must be the number of
            * sites and must match the number of sites in the DataHolder
            * that will be passed to ms(). The positions will be evenly
            * spread between 0 and 1.
            *
            */
            void ms_auto_positions(unsigned int n);

           /** \brief %Export data as ms format
            *
            * The number of sites must have been previously (and to the
            * correct value) using either ms_num_positions() (with all
            * positions specified with ms_position()) or
            * ms_auto_positions(). The format is as follow: one line
            * with two slashes; one line with the number of sites; one
            * line with the positions, or an empty line if the number of
            * sites is zero; then the matrix of genotypes (one line per
            * sample), only if the number of sites is larger than zero.
            *
            * \param data the data set to export.
            *
            * \param spacer if ``True``, insert a space between each
            * genotype value.
            *
            * \param outgroup if ``True``, outgroup data is exported
            * after the main group.
            *
            */
            void ms(const DataHolder& data, bool spacer, bool outgroup=true);

        private:

           /** \brief This call cannot be copied. */
            Export(const Export& src) {}

           /** \brief This call cannot be copied. */
            Export& operator=(const Export& src) { return * this;}

            void _newick(const Tree& tree, const Node * node, bool blen);

            unsigned int _ms_res_positions;
            double * _ms_positions;
    };
}

#endif
