/*
    Copyright 2011,2013,2015 Stéphane De Mita, Mathieu Siol

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

#include "ABC.hpp"
#include "egglib.hpp"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>

#ifndef __woABC
#include "gsl/gsl_multifit.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_errno.h"
#endif

namespace egglib {

    void egghandler(const char * reason, const char * file, int line, int gsl_errno) {
        std::ostringstream stream;
        #ifndef __woABC
        stream << "[GSL error] " << gsl_strerror(gsl_errno) << " - " << reason;
        #endif
        throw EggRuntimeError(stream.str().c_str());
    }

    ////

    ABC::ABC() {
        #ifdef __woABC
        throw EggRuntimeError("EggLib has been compiled without support for the ABC class");
        #endif
        init();
    }

    ////

    ABC::~ABC() {
        free();
    }

    ////

    void ABC::init() {
        _sd = NULL;
        _obs = NULL;
        _sum = NULL;
        _sum2 = NULL;
        _nstats_c = 0;
        _fnames = NULL;
        _nparams = NULL;
        _lnames = NULL;
        _nsam_part = NULL;
        _nfiles_c = 0;
        #ifndef __woABC
        gsl_set_error_handler(&egghandler);
        #endif
        reset();
    }

    ////

    void ABC::free() {
        if (_sd) ::free(_sd);
        if (_obs) ::free(_obs);
        if (_sum) ::free(_sum);
        if (_sum2) ::free(_sum2);
        for (unsigned int i=0; i<_nfiles_c; i++) if (_fnames[i]) ::free(_fnames[i]);
        if (_fnames) ::free(_fnames);
        if (_lnames) ::free(_lnames);
        if (_nparams) ::free(_nparams);
        if (_nsam_part) ::free(_nsam_part);
    }

    ////

    void ABC::reset() {
        _nstats = 0;
        _nsam = 0;
        _threshold = -1;
        _nfiles = 0;
    }
     
    ////

    void ABC::number_of_statistics(unsigned int ns) {

        reset();
        _nstats = ns;

        if (_nstats > _nstats_c) {
            _sd = (double *) realloc(_sd, ns * sizeof(double));
            if (!_sd) throw EGGMEM;
            _obs = (double *) realloc(_obs, ns * sizeof(double));
            if (!_obs) throw EGGMEM;
            _sum = (double *) realloc(_sum, ns * sizeof(double));
            if (!_sum) throw EGGMEM;
            _sum2 = (double *) realloc(_sum2, ns * sizeof(double));
            if (!_sum2) throw EGGMEM;

            _nstats_c = _nstats;
        }

        for (unsigned int i=0; i<ns; i++) {
            _sd[i] = 0.;
            // sum and sum2 are initialized later
        }
    }

    ////

    void ABC::add_fname(const char * fname, unsigned int np) {

        _nfiles++;
        if (_nfiles > _nfiles_c) {
            _fnames = (char **) realloc(_fnames, _nfiles * sizeof(char *));
            if (!_fnames) throw EGGMEM;
            _lnames = (unsigned int *) realloc(_lnames, _nfiles * sizeof(char *));
            if  (!_lnames) throw EGGMEM;
            _nparams = (unsigned int *) realloc(_nparams, _nfiles * sizeof(unsigned int));
            if (!_nparams) throw EGGMEM;
            _nsam_part = (unsigned int *) realloc(_nsam_part, _nfiles * sizeof(unsigned int));
            if (!_nsam_part) throw EGGMEM;
            _fnames[_nfiles-1] = NULL;
            _lnames[_nfiles-1] = 0;
            _nsam_part[_nfiles-1] = 0;
            _nfiles_c++;
        }

        _nparams[_nfiles-1] = np;
        unsigned int len = strlen(fname);
        if (len > _lnames[_nfiles-1]) {
            _fnames[_nfiles-1] = (char *) realloc(_fnames[_nfiles-1], (len+1) * sizeof(char));
            if (!_fnames[_nfiles-1]) throw EGGMEM;
            _lnames[_nfiles-1] = len;
        }
        strcpy(_fnames[_nfiles-1], fname);
    }

    ////

    double ABC::sd(unsigned int index) const {
        return _sd[index];
    }

    ////

    void ABC::obs(unsigned int index, double value) {
        _obs[index] = value;
    }

    ////

    void ABC::get_threshold(double tolerance) {
        
        // common variables

        std::string line;
        std::string token;
        double value;
        unsigned int linenum;

        for (unsigned int i=0; i<_nstats; i++) {
            _sum[i] = 0;
            _sum2[i] = 0;
        }

        // gets data for standard deviations
        
        for (unsigned int i=0; i<_nfiles; i++) {

            std::ifstream fstream(_fnames[i]);
            linenum = 0;
            if (!fstream.good()) {
                throw EggOpenFileError(_fnames[i]);
            }

            while (!fstream.eof()) {
                if (!fstream.good()) {
                    throw EggFormatError(_fnames[i], linenum, "ABC sample data", "expected a new line");
                }

                getline(fstream, line);
                linenum++;
                std::istringstream sstream(line);

                if (line.size() == 0) break;

                for (unsigned int j=0; j<_nparams[i]; j++) {
                    if (!sstream.good()) {
                        throw EggFormatError(_fnames[i], linenum, "ABC sample data", "invalid line");
                    }
                    sstream >> token;
                }

                if (!sstream.good()) {
                    throw EggFormatError(_fnames[i], linenum, "ABC sample data", "invalid line");
                }

                sstream >> token;
                if (token != "#") {
                    throw EggFormatError(_fnames[i], linenum, "ABC sample data", "invalid line (invalid number of parameters?)");
                }

                for (unsigned int j=0; j<_nstats; j++) {
                    if (!sstream.good()) {
                        throw EggFormatError(_fnames[i], linenum, "ABC sample data", "invalid line");
                    }
                    sstream >> token;
                    value = atof(token.c_str());
                    _sum[j] += value;
                    _sum2[j] += (value*value);
                }

                if (sstream.good()) throw EggFormatError(_fnames[i], linenum, "ABC sample data", "line longer than expected");

                _nsam++;
                _nsam_part[i]++;
            }

            fstream.close();
        }

        if (_nsam==0) throw EggArgumentValueError("ABC: no sample data found");

        // computes standard deviations

        double m, m2;

        for (unsigned int i=0; i<_nstats; i++) {
            m = _sum[i] / _nsam;
            m2 = _sum2[i] / _nsam;
            _sd[i] = sqrt(m2 - m*m);
        }

        // computes Euclidean distances

        _euclid.resize(_nsam);
        unsigned int cnt = 0;

        for (unsigned int i=0; i<_nfiles; i++) {

            std::ifstream fstream(_fnames[i]);
            linenum = 0;
            if (!fstream.good()) throw EggOpenFileError(_fnames[i]);

            while (!fstream.eof()) {

                if (!fstream.good()) throw EggFormatError(_fnames[i], linenum, "ABC sample data", "expected a new line");

                getline(fstream, line);
                std::istringstream sstream(line);
                linenum++;

                if (line.size() == 0) break;

                if (cnt>=_nsam) throw EggArgumentValueError("ABC: trying to compute more euclidean distances than the number of sample data used for computed statistics standard deviatons");

                for (unsigned int j=0; j<_nparams[i]; j++) {
                    if (!sstream.good()) {
                        throw EggFormatError(_fnames[i], linenum, "ABC sample data", "invalid line");
                    }
                    sstream >> token;
                }

                if (!sstream.good()) {
                    throw EggFormatError(_fnames[i], linenum, "ABC sample data", "invalid line");
                }
                
                sstream >> token;
                if (token != "#") {
                    throw EggFormatError(_fnames[i], linenum, "ABC sample data", "invalid line (invalid number of parameters?)");
                }

                _euclid[cnt] = 0.;

                for (unsigned int j=0; j<_nstats; j++) {
                    if (!sstream.good()) throw EggFormatError(_fnames[i], linenum, "ABC sample data", "invalid line");
                    sstream >> token;
                    value = atof(token.c_str());
                    if (_sd[j]>0) {
                        value /= _sd[j];
                        _euclid[cnt] += ( (value-_obs[j]/_sd[j]) * 
                                          (value-_obs[j]/_sd[j]) );
                    }
                }

                if (sstream.good()) throw EggFormatError(_fnames[i], linenum, "ABC sample data", "line longer than expected");

                _euclid[cnt] = sqrt(_euclid[cnt]);
                cnt++;
            }
            fstream.close();
        }

        if (cnt!=_nsam) throw EggRuntimeError("ABC: inconsistent number of data between two steps to data access");

        unsigned int position = (unsigned int)ceil(_nsam * tolerance);

        if (position<2) {
            throw EggArgumentValueError("tolerance is not large enough to catch at least two points");
        }

        std::sort(_euclid.begin(), _euclid.end());
        _threshold = _euclid[position];
    }

    ////

    unsigned int ABC::rejection(const char* outfname, bool exportlabels, bool strip) {

        unsigned int n = 0;
        unsigned int accept = 0;
        double euclid;
        double value;
        double w;
        std::string line;
        std::string token;
        unsigned int linenum;

        std::ofstream outfstream(outfname);
        if (!outfstream.good()) throw EggOpenFileError(outfname);

        for (unsigned int i=0; i<_nfiles; i++) {

            std::ifstream infstream(_fnames[i]);
            if (!infstream.good()) throw EggOpenFileError(_fnames[i]);
            linenum = 0;

            while (!infstream.eof()) {

                if (!infstream.good()) throw EggFormatError(_fnames[i], linenum, "ABC sample data", "expected a new line");

                getline(infstream, line);
                std::istringstream sstream(line);
                linenum++;

                if (line.size() == 0) break;

                for (unsigned int j=0; j<_nparams[i]; j++) {
                    if (!sstream.good()) {
                        throw EggFormatError(_fnames[i], linenum, "ABC sample data", "invalid line");
                    }
                    sstream >> token;
                }

                if (!sstream.good()) {
                    throw EggFormatError(_fnames[i], linenum, "ABC sample data", "invalid line");
                }

                sstream >> token;
                if (token != "#") {
                    throw EggFormatError(_fnames[i], linenum, "ABC sample data", "invalid line (invalid number of parameters?)");
                }

                euclid = 0.;

                for (unsigned int j=0; j<_nstats; j++) {
                    if (!sstream.good()) throw EggFormatError(_fnames[i], linenum, "ABC sample data", "invalid line");
                    sstream >> token;
                    value = atof(token.c_str());
                    if (_sd[j]>0) {
                        value /= _sd[j];
                        euclid += ( (value-_obs[j]/_sd[j]) * 
                                    (value-_obs[j]/_sd[j]) );
                    }
                }

                euclid = sqrt(euclid);
                n++;
                if (euclid < _threshold) {
                    accept++;
                    w = 1-((euclid*euclid)/(_threshold*_threshold));
                    if (exportlabels) outfstream << "[" << i+1 << "] ";
                    if (strip) outfstream << line.substr(0, line.find('#')-1) << std::endl;
                    else outfstream << line << " # " << w << std::endl;
                }
            }

            infstream.close();
        }

        if (n != _nsam) {
            throw EggRuntimeError("number of samples do not match between ABC::get_sd() and ABC::threshold() calls");
        }

        outfstream.close();
        return accept;        
    }

    //////

    unsigned int ABC::regression(const char* infname, const char* outfname,
                                 TransformMode mode, const char* header) {

        // first pass to get the number of points and checks file consistency
        
        std::ifstream fstream(infname);
        if (!fstream.good()) throw EggOpenFileError(infname);

        std::string line;
        std::string token;
        unsigned int linenum = 0;

        unsigned int nparams = 0;
        unsigned int npoints = 0;
        
        while (!fstream.eof()) {

            if (!fstream.good()) throw EggFormatError(infname, linenum, "ABC sample data", "expected a new line");

            getline(fstream, line);
            linenum++;
            if (line.size()==0) break;
            if (line[0]=='[') throw EggFormatError(infname, linenum, "ABC sample data", "model labels are not supported for regression");

            std::istringstream sstream(line);
            unsigned int i=0;

            while (1) {
                if (!sstream.good()) throw EggFormatError(infname, linenum, "ABC sample data", "invalid line structure  (3)");
                sstream >> token;
                if (token == "#") break;
                i++;
            }

            if (i==0) throw EggFormatError(infname, linenum, "ABC sample data", "no params?");
            if (nparams!=0 && i!=nparams) throw EggFormatError(infname, linenum, "ABC sample data", "inconsistent number of parameters");
            if (nparams==0) nparams = i;
            npoints++;
        }

        if (npoints==0) throw EggRuntimeError("cannot perform ABC: no simulations");

        fstream.close();

        // allocates memory

        #ifndef __woABC
        gsl_multifit_linear_workspace* workspace =
                    gsl_multifit_linear_alloc(npoints, _nstats+1);  // +1 for intercept

        gsl_matrix* mstats = gsl_matrix_alloc(npoints, _nstats+1); // +1 for intercept
        gsl_vector* vweights = gsl_vector_alloc(npoints);
        gsl_matrix* mparams = gsl_matrix_alloc(npoints, nparams);
        gsl_vector* vparams = gsl_vector_alloc(npoints);
        gsl_vector* vcoefs = gsl_vector_alloc(_nstats+1); // +1 for intercept
        gsl_matrix* vcov = gsl_matrix_alloc(_nstats+1, _nstats+1);  // +1 for intercept

        if (!mstats || !workspace || !vweights ||
            !vparams || !mparams || !vcoefs || !vcov) {
            if (mstats) gsl_matrix_free(mstats);
            if (vweights) gsl_vector_free(vweights);
            if (mparams) gsl_matrix_free(mparams);
            if (vparams) gsl_vector_free(vparams);
            if (vcoefs) gsl_vector_free(vcoefs);
            if (vcov) gsl_matrix_free(vcov);
            if (workspace) gsl_multifit_linear_free(workspace);
            throw EGGMEM;
        }
        #endif

        // loads data

        fstream.open(infname);
        if (!fstream.good()) throw EggOpenFileError(infname);
        linenum = 0;

        unsigned int point=0;

        while (!fstream.eof()) {
            if (!fstream.good()) throw EggFormatError(infname, linenum, "ABC sample data", "expected a new line");

            getline(fstream, line);
            if (line.size() == 0) break;

            std::istringstream sstream(line);

            // gets params

            for (unsigned int param=0; param<nparams; param++) {
                if (!sstream.good()) throw EggFormatError(infname, linenum, "ABC sample data", "invalid line structure  (6)");
                sstream >> token;
                #ifndef __woABC
                gsl_matrix_set(mparams, point, param, atof(token.c_str()));
                #endif
            }

            sstream >> token;
            if (token != "#") throw EggFormatError(infname, linenum, "ABC sample data", "invalid line structure  (7)");

            // gets stats

            #ifndef __woABC
            gsl_matrix_set(mstats, point, 0, 1.); // for intercept
            #endif

            for (unsigned int stat=0; stat<_nstats; stat++) {
                if (!sstream.good()) throw EggFormatError(infname, linenum, "ABC sample data", "invalid line structure  (8)");
                sstream >> token;
                if (_sd[stat]==0) throw EggRuntimeError("at least one statistic not variable in the local region (or ABC class not used properly)");
                #ifndef __woABC
                gsl_matrix_set(mstats, point, stat+1, atof(token.c_str())/_sd[stat]);
                #endif
            }
            sstream >> token;
            if (token != "#") throw EggFormatError(infname, linenum, "ABC sample data", "invalid line structure  (9)");

            // gets weight

            if (!sstream.good()) throw EggFormatError(infname, linenum, "ABC sample data", "invalid line structure  (10)");
            sstream >> token;
            #ifndef __woABC
            gsl_vector_set(vweights, point, atof(token.c_str()));
            #endif

            point++;
        }

        if (point!=npoints) throw EggRuntimeError("inconsistent ABC file parsing...");

        fstream.close();

        // fits all parameters separately

        #ifndef __woABC
        double chisq;
        double mean;
        #endif

        for (unsigned int param=0; param<nparams; param++) {

            // loads the explained variable

            for (unsigned int point=0; point<npoints; point++) {
                #ifndef __woABC
                gsl_vector_set(vparams, point,
                                gsl_matrix_get(mparams, point, param));
                #endif
            }

            // performs transformation

            #ifndef __woABC
            double min=0.0;
            double max=999999999.0;
            double x=0.;
            #endif

            if (mode==TAN) {
                #ifndef __woABC
                min = gsl_vector_get(vparams, 0);
                max = gsl_vector_get(vparams, 0);
                for (unsigned int point=1; point<npoints; point++) {
                    x = gsl_vector_get(vparams, point);
                    if (x<min) min=x;
                    if (x>max) max=x;
                }
                min -= 0.00000001;
                max += 0.00000001;
                for (unsigned int point=0; point<npoints; point++) {
                    x = gsl_vector_get(vparams, point);
                    x = -log(1./(tan( ((x-min)/(max-min))*(M_PI/2.) )));
                    gsl_vector_set(vparams, point, x);
                }
                #endif
            }

            if (mode==LOG) {
                #ifndef __woABC
                for (unsigned int point=0; point<npoints; point++) {
                    x = gsl_vector_get(vparams, point);
                    if (x<=0.) throw EggRuntimeError("ABC transformation LOG not available for parameter values <= 0");
                    gsl_vector_set(vparams, point, log(x));
                }
                #endif
            }

            // performs regression

            #ifndef __woABC
            gsl_multifit_wlinear(mstats, vweights, vparams,
                                    vcoefs, vcov, &chisq, workspace);
            #endif

            // computes predicted mean

            #ifndef __woABC
            mean = 1. * gsl_vector_get(vcoefs, 0);   // intercept
            for (unsigned int stat=0; stat<_nstats; stat++) {
                mean += _obs[stat]/_sd[stat] * gsl_vector_get(vcoefs, stat+1);
            }
            #endif

            // computes residuals and predicted params (overwrites params table)

            #ifndef __woABC
            for (unsigned int point=0; point<npoints; point++) {
                double observed = gsl_vector_get(vparams, point);
                double predicted = 0.;
                for (unsigned int stat=0; stat<_nstats+1; stat++) {
                    predicted += gsl_matrix_get(mstats, point, stat)
                                                * gsl_vector_get(vcoefs, stat);
                }

                x = mean + (observed - predicted);

                // untransforms immediately

                if (mode==LOG) x = exp(x);
                if (mode==TAN) x = min + (2./M_PI)*(max-min)*atan(exp(x));

                // and sets final value

                gsl_matrix_set(mparams, point, param, x);
            }
            #endif
        }

        // exports data

        std::ofstream out(outfname);
        if (!out.good()) throw EggOpenFileError(outfname);

        if (strlen(header)) out << header << std::endl;

        #ifndef __woABC
        for (unsigned int point=0; point<npoints; point++) {
            out << gsl_matrix_get(mparams, point, 0);
            for (unsigned int param=1; param<nparams; param++) {
                out << " " << gsl_matrix_get(mparams, point, param);
            }
            out << std::endl;
        }
        #endif

        out.close();

        // frees memory

        #ifndef __woABC
        gsl_multifit_linear_free(workspace);
        gsl_matrix_free(mstats);
        gsl_vector_free(vweights);
        gsl_matrix_free(mparams);
        gsl_vector_free(vparams);
        gsl_vector_free(vcoefs);
        gsl_matrix_free(vcov);
        #endif

        return npoints;        
    }

    ////

    unsigned int ABC::number_of_samples() const {
        return _nsam;
    }

    ////

    unsigned int ABC::number_of_samples_part(unsigned int i) const {
        return _nsam_part[i];
    }

    ////

    double ABC::threshold() const {
        return _threshold;
    }
}
