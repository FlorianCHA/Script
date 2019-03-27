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

#include "egglib.hpp"
#include "Params.hpp"
#include "Coalesce.hpp"
#include "Params.hpp"
#include <cstdlib>
#include <sstream>
#include <cmath>

namespace egglib {

    Migration::Migration(unsigned int n, double M) {
        npop = n;
        npop_reserved = n;
        matrix = (double**) malloc(n * sizeof(double*));
        if (!matrix) throw EGGMEM;
        for (unsigned int i=0; i<n; i++) {
            matrix[i] = (double*) malloc(n * 2 * sizeof(double));
            if (!matrix[i]) throw EGGMEM;
        }
        set_all(M);
    }

    Migration::~Migration() {
        if (matrix) {
            for (unsigned int i=0; i<npop_reserved; i++) if (matrix[i]) free(matrix[i]);
            free(matrix);
        }
    }

    unsigned int Migration::n() const {
        return npop;
    }

    void Migration::set_all(double M) {
        for (unsigned int i=0; i<npop; i++) {
            matrix[i][2*i] = M;
            matrix[i][2*i+1] = M; // cached value
            for (unsigned int j=i+1; j<npop; j++) {
                matrix[i][2*j] = M/(npop-1);
                matrix[j][2*i] = M/(npop-1);
                matrix[i][2*j+1] = M/(npop-1);  // cache
                matrix[j][2*i+1] = M/(npop-1);  // cache
            }
        }
    }

    void Migration::set_all_R(double M) {
        for (unsigned int i=0; i<npop; i++) {
            matrix[i][2*i] = M; // set only current value
            for (unsigned int j=i+1; j<npop; j++) { // set only current values
                matrix[i][2*j] = M/(npop-1);
                matrix[j][2*i] = M/(npop-1);
            }
        }
    }

    void Migration::set_row(unsigned int i, double M) {
        matrix[i][2*i] = M;
        matrix[i][2*i+1] = M;
        for (unsigned int j=0; j<npop; j++) {
            if (i==j) continue;
            matrix[i][2*j] = M/(npop-1);
            matrix[i][2*j+1] = M/(npop-1);
        }
    }

    void Migration::set_row_R(unsigned int i, double M) {
        matrix[i][2*i] = M;
        for (unsigned int j=0; j<npop; j++) {
            if (i==j) continue;
            matrix[i][2*j] = M/(npop-1);
        }
    }

    void Migration::set_pair(unsigned int i, unsigned int j, double m) {
        matrix[i][2*i] -= matrix[i][2*j];
        matrix[i][2*i] += m;
        matrix[i][2*i+1] = matrix[i][2*i];  // cache
        matrix[i][2*j] = m;
        matrix[i][2*j+1] = m;  // cache
    }

    void Migration::set_pair_R(unsigned int i, unsigned int j, double m) {
        matrix[i][2*i] -= matrix[i][2*j];
        matrix[i][2*i] += m;
        matrix[i][2*j] = m;
    }

    double Migration::get_row(unsigned int i) const {
        return matrix[i][2*i];
    }

    double Migration::get_pair(unsigned int i, unsigned int j) const {
        return matrix[i][2*j];
    }

    void Migration::restore() {
        for (unsigned int i=0; i<npop; i++) {
            matrix[i][2*i] = matrix[i][2*i+1];
            for (unsigned int j=i+1; j<npop; j++) {
                matrix[i][2*j] = matrix[i][2*j+1];
                matrix[j][2*i] = matrix[j][2*i+1];
            }
        }
    }

    const double Event::_small = 0.00000000001;

    Event::Event(Type type, double date) {
        _prev = NULL;
        _next = NULL;
        _type = type;
        _date = date;
        _index = MAX;
        _dest = 0;
        _param = 0.0;
        _number1 = 0;
        _number2 = 0;
    }

    void Event::copy(const Event& src) {
        _prev = NULL;
        _next = NULL;
        _type = src._type;
        _date = src._date;
        _index = src._index;
        _dest = src._dest;
        _param = src._param;
        _number1 = src._number1;
        _number2 = src._number2;
    }

    Event::Type Event::event_type() const {
        return _type;
    }

    double Event::date() const {
        return _date;
    }

    double Event::get_param() const {
        return _param;
    }

    unsigned int Event::get_index() const {
        return _index;
    }

    unsigned int Event::get_dest() const {
        return _dest;
    }

    unsigned int Event::get_number1() const {
        return _number1;
    }

    unsigned int Event::get_number2() const {
        return _number2;
    }

    void Event::insert(Event * event) {
        if (_next == NULL) {
            _next = event;
            event->_prev = this;
            event->_next = NULL;
            return;
        }

        if (_next->_date - event->_date > _small) {
            event->_next = _next;
            _next->_prev = event;
            _next = event;
            event->_prev = this;
            return;
        }

        _next->insert(event);
    }

    void Event::disconnect() {
        if (_next != NULL) {
            _next->_prev = NULL;
            _next->disconnect();
            _next = NULL;
        }
    }

    Event * Event::next() {
        return _next;
    }

    Event * Event::prev() {
        return _prev;
    }

    void Event::move(double date) {

        // store date

        _date = date;

        // if this is earlier than previous, move itself up

        if (_prev != NULL && _prev->_date - date > _small) {

            // connect previous to next
            _prev->_next = _next;
            if (_next != NULL) _next->_prev = _prev;

            // reconnect earlier (will overwrite links)
            _prev->_insert_up(this);
        }

        // if this is later than next, move itself down (similar)

        if (_next != NULL && date - _next->_date > _small) {
            _prev->_next = _next;
            if (_next != NULL) _next->_prev = _prev;
            _next->insert(this);
        }
    }

    void Event::set_index(unsigned int i) {
        _index = i;
    }

    void Event::set_param(double p) {
        _param = p;
    }

    void Event::set_dest(unsigned int d) {
        _dest = d;
    }

    void Event::set_number1(unsigned int n) {
        _number1 = n;
    }

    void Event::set_number2(unsigned int n) {
        _number2 = n;
    }

    void Event::_insert_up(Event * event) {  // event must be before this (and there must be an event before)
        if (_prev->_date - event->_date > _small) {
            _prev->_insert_up(event);
        }
        else {
            _prev->_next = event;
            event->_prev = _prev;
            _prev = event;
            event->_next = this;
        }
    }

    void Event::perform(Params * params, Coalesce * coal) {

        switch (_type) {

            case delayed:
                coal->delayedSample(_date, _index, _number1, _number2);
                break;

            case change_N:
                if (_index == MAX) for (unsigned int i=0; i<params->k(); i++) params->N_R(i, _param, _date);
                else params->N_R(_index, _param, _date);
                break;

            case change_G:
                if (_index == MAX) for (unsigned int i=0; i<params->k(); i++) params->G_R(i, _param, _date);
                else params->G_R(_index, _param, _date);
                break;

            case change_R:
                params->R_R(_param);
                break;

            case change_s:
                if (_index == MAX) for (unsigned int i=0; i<params->k(); i++) params->s_R(i, _param, _date);
                else params->s_R(_index, _param, _date);
                break;

            case change_M:
                params->M().set_all_R(_param);
                break;

            case change_Mp:
                params->M().set_pair_R(_index, _dest, _param);
                break;

            case admixture:
                coal->admixt(_index, _dest, _param);
                break;

            case bottleneck:
                if (_index == MAX) for (unsigned int i=0; i<params->k(); i++) coal->bottleneck(i, _param);
                else coal->bottleneck(_index, _param);
                break;

            case none:
                throw EggRuntimeError("cannot use \"none\" as an actual historical parameter change in a simulation");
        }
    }

    Params::Params() {
        _init(1, 0.0);
    }

    Params::Params(unsigned int npop, double migr) {
        _init(npop, migr);
    }

    Params::~Params() {
        if (_base_change) delete _base_change;
        if (_migr) delete _migr;
        if (_n1) free(_n1);
        if (_n2) free(_n2);
        if (_popsize) free(_popsize);
        if (_popsize_cache) free(_popsize_cache);
        if (_growthrate) free(_growthrate);
        if (_growthrate_cache) free(_growthrate_cache);
        if (_selfing) free(_selfing);
        if (_selfing_cache) free(_selfing_cache);
        if (_lastChange) free(_lastChange);
        for (unsigned int i=0; i<_nalleles_r; i++) {
            if (_transW[i]) free(_transW[i]);
        }
        if (_transW) free(_transW);
        if (_sitePos) free(_sitePos);
        if (_siteW) free(_siteW);
    }

    std::string Params::summary() const {

        std::ostringstream stream;

        stream << "Number of populations: " << _npop << std::endl;
        for (unsigned int i=0; i<_npop; i++) {
            stream << "    Population " << i+1 << ":" << std::endl;
            stream << "       Single samples: " << _n1[i] << std::endl;
            stream << "       Double samples: " << _n2[i] << std::endl;
            stream << "       Relative size: " << _popsize[i] << std::endl;
            stream << "       Growth rate: " << _growthrate[i] << std::endl;
            stream << "       Selfing rate: " << _selfing[i] << std::endl;
        }
        stream << "Recombination rate: " << _recomb << std::endl;

        stream << "Migration matrix:" << std::endl;

        double m;
        std::ostringstream ss;
        ss.precision(3);
        std::string s;
        unsigned int ln;

        for (unsigned int i=0; i<_npop; i++) {
            stream << "  ";
            for (unsigned int j=0; j<_npop; j++) {
                if (i==j) m = _migr->get_row(i);
                else m = _migr->get_pair(i, j);
                ss.str("");
                ss.clear();
                ss << m;
                s = ss.str();
                ln = (9 - s.size()) / 2;
                for (unsigned int k=0; k<ln; k++) stream << " ";
                stream << s;
                for (unsigned int k = ln+s.size(); k<9; k++) stream << " ";
            }
            stream << std::endl;
        }

        stream << "Mutation rate: " << _theta << std::endl;
        stream << "Fixed number of alleles: " << _fixed << std::endl;
        stream << "Mutation model: ";
        switch (_mutmodel) {
            case KAM:
                stream << "KAM";
                break;
            case IAM:
                stream << "IAM";
                break;
            case SMM:
                stream << "SMM";
                break;
            case TPM:
                stream << "TPM";
                break;
        }
        stream << std::endl;

        stream << "Number of alleles: " << _nalleles << std::endl;
        stream << "Random start allele: " << (_random_start_allele) << std::endl;

        if (_mutmodel == TPM) {
            stream << "TPM probability parameter: " << _TPMproba << std::endl;
            stream << "TPM shape parameter: " << _TPMparam << std::endl;
        }

        stream << "Custom transition matrix: " << (_transW_matrix) << std::endl;

        double r;

        for (unsigned int i=0; i<_nalleles; i++) {
            stream << "  ";
            for (unsigned int j=0; j<_nalleles; j++) {
                if (i==j) r = get_transW_row(i);
                else r = get_transW_pair(i, j);
                ss.str("");
                ss.clear();
                ss << r;
                s = ss.str();
                ln = (9 - s.size()) / 2;
                for (unsigned int k=0; k<ln; k++) stream << " ";
                stream << s;
                for (unsigned int k = ln+s.size(); k<9; k++) stream << " ";
            }
            stream << std::endl;
        }

        stream << "Number of mutable sites: " << _nsites << std::endl;
        for (unsigned int i=0; i<_nsites; i++) {
            stream << "    Site " << i+1 << ":" << std::endl;
            stream << "        Position: " << _sitePos[i] << std::endl;
            stream << "        Site weight: " << _siteW[i] << std::endl;
        }

        stream << "Number of changes: " << _num_changes << std::endl;

        Event * cur = _base_change->next();
        unsigned int c = 0;

        while (cur != NULL) {
            switch (cur->event_type()) {

                case Event::delayed:
                    stream << "    Change " << c+1 << ": Delayed sample" << std::endl;
                    stream << "        Date: "           << cur->date() << std::endl;
                    stream << "        Population: "     << cur->get_index() << std::endl;
                    stream << "        Label: "          << cur->get_dest() << std::endl;
                    stream << "        Single samples: " << cur->get_number1() << std::endl;
                    stream << "        Double samples: " << cur->get_number2() << std::endl;
                    break;

                case Event::change_N:
                    stream << "    Change " << c+1 << ": Population size change" << std::endl;
                    stream << "        Date: "       << cur->date() << std::endl;
                    stream << "        Population: " << cur->get_index() << std::endl;
                    stream << "        Size: "       << cur->get_param() << std::endl;
                    break;

                case Event::change_G:
                    stream << "    Change " << c+1 << ": Growth rate change" << std::endl;
                    stream << "        Date: "       << cur->date() << std::endl;
                    stream << "        Population: " << cur->get_index() << std::endl;
                    stream << "        Rate: "       << cur->get_param() << std::endl;
                    break;

                case Event::change_M:
                    stream << "    Change " << c+1 << ": All migration rates change" << std::endl;
                    stream << "        Date: "       << cur->date() << std::endl;
                    stream << "        Rate: "       << cur->get_param() << std::endl;
                    break;

                case Event::bottleneck:
                    stream << "    Change " << c+1 << ": Bottleneck" << std::endl;
                    stream << "        Date: "       << cur->date() << std::endl;
                    stream << "        Population: " << cur->get_index() << std::endl;
                    stream << "        Strength: "   << cur->get_param() << std::endl;
                    break;

                case Event::change_s:
                    stream << "    Change " << c+1 << ": All selfing rates change" << std::endl;
                    stream << "        Date: "   << cur->date() << std::endl;
                    stream << "        Rate: "   << cur->get_param() << std::endl;
                    break;

                case Event::admixture:
                    stream << "    Change " << c+1 << ": Admixture" << std::endl;
                    stream << "        Date: "             << cur->date() << std::endl;
                    stream << "        Population: "       << cur->get_index() << std::endl;
                    stream << "        Other population: " << cur->get_dest() << std::endl;
                    stream << "        Probability: "      << cur->get_param() << std::endl;
                    break;

                case Event::change_Mp:
                    stream << "    Change " << c+1 << ": Pairwise migration rate change" << std::endl;
                    stream << "        Date: "         << cur->date() << std::endl;
                    stream << "        Source: "       << cur->get_index() << std::endl;
                    stream << "        Destination: "  << cur->get_dest() << std::endl;
                    stream << "        Rate: "         << cur->get_param() << std::endl;
                    break;

                case Event::change_R:
                    stream << "    Change " << c+1 << ": Recombination rate change" << std::endl;
                    stream << "        Date: "         << cur->date() << std::endl;
                    stream << "        Rate: "         << cur->get_param() << std::endl;
                    break;

                case Event::none:
                    throw EggRuntimeError("cannot use \"none\" as an actual historical parameter change in a simulation");
            }
            cur = cur->next();
            c++;
        }

        return stream.str();
    }

    void Params::set_outgroup(unsigned int lbl) {
        _outgroup = lbl;
    }

    unsigned int Params::get_outgroup() const {
        return _outgroup;
    }

    unsigned int Params::k() const {
        return _npop;
    }

    unsigned int Params::get_K() const {
        return _nalleles;
    }

    void Params::set_K(unsigned int value) {
        if (_transW_matrix) alloc_nalleles(value);
        else _nalleles = value;
    }

    unsigned long Params::get_max_iter() const {
        return _max_iter;
    }

    void Params::set_max_iter(unsigned long x) {
        _max_iter = x;
    }

    bool Params::get_random_start_allele() const {
        return _random_start_allele;
    }

    void Params::set_random_start_allele(bool value) {
        _random_start_allele = value;
    }

    void Params::set_R(double value) {
        _recomb = value;
        _recomb_cache = value;
    }

    void Params::R_R(double value) {
        _recomb = value;
    }

    double Params::get_R() const {
        return _recomb;
    }

    void Params::set_theta(double value) {
        _theta = value;
    }

    double Params::get_theta() const {
        return _theta;
    }

    void Params::set_fixed(unsigned int value) {
        _fixed = value;
    }

    unsigned int Params::get_fixed() const {
        return _fixed;
    }

    void Params::set_mutmodel(MutationModel value) {
        _mutmodel = value;
    }

    Params::MutationModel Params::get_mutmodel() const {
        return _mutmodel;
    }

    void Params::set_TPMproba(double value) {
        _TPMproba = value;
    }

    double Params::get_TPMproba() const {
        return _TPMproba;
    }

    void Params::set_TPMparam(double value) {
        _TPMparam = value;
    }

    double Params::get_TPMparam() const {
        return _TPMparam;
    }

    void Params::set_n1(unsigned int pop, unsigned int value) {
        _n1[pop] = value;
    }

    unsigned int Params::get_n1(unsigned int pop) const {
        return _n1[pop];
    }

    void Params::set_n2(unsigned int pop, unsigned int value) {
        _n2[pop] = value;
    }

    unsigned int Params::get_n2(unsigned int pop) const {
        return _n2[pop];
    }

    void Params::set_N(unsigned int pop, double value) {
        _popsize[pop] = value;
        _popsize_cache[pop] = value;
    }

    void Params::N_R(unsigned int pop, double value, double t) {
        _popsize[pop] = value;
        _lastChange[pop] = t;
    }

    double Params::get_N(unsigned int pop) const {
        return _popsize[pop];
    }

    void Params::set_G(unsigned int pop, double value) {
        _growthrate[pop] = value;
        _growthrate_cache[pop] = value;
    }

    void Params::G_R(unsigned int pop, double value, double t) {
        _popsize[pop] = _popsize[pop] * exp( - _growthrate[pop] *
                                            (t - _lastChange[pop]));
        _growthrate[pop] = value;
        _lastChange[pop] = t;
    }

    double Params::get_G(unsigned int pop) const {
        return _growthrate[pop];
    }

    void Params::set_s(unsigned int pop, double value) {
        _selfing[pop] = value;
        _selfing_cache[pop] = value;        
    }

    void Params::s_R(unsigned int pop, double value, double t) {
        _selfing[pop] = value;
    }

    double Params::get_s(unsigned int pop) const {
        return _selfing[pop];
    }

    unsigned int Params::get_L() const {
        return _nsites;
    }

    void Params::set_L(unsigned int value) {
        alloc_nsites(value);
    }

    void Params::autoSitePos() {

        // if no sites, does nothing

        if (_nsites == 0) return;

        // if one site, at the middle

        if (_nsites == 1) {
            _sitePos[0] = 0.5;
            return;
        }

        // if >=2 sites - define the gap between two sites

        double d = 1. / (_nsites-1);

        // set positions

        _sitePos[0] = 0.0;
        _sitePos[_nsites-1] = 1.0;

        for (unsigned int i=1; i<_nsites-1; i++) {
            _sitePos[i] = i * d;
        }
    }

    void Params::set_sitePos(unsigned int site, double value) {
        _sitePos[site] = value;
    }

    double Params::get_sitePos(unsigned int site) const {
        return _sitePos[site];
    }

    void Params::set_siteW(unsigned int site, double value) {
        _totalSiteW -= _siteW[site];
        _siteW[site] = value;
        _totalSiteW += _siteW[site];
    }

    double Params::get_siteW(unsigned int site) const {
        return _siteW[site];
    }

    double Params::totalSiteW() const {
        return _totalSiteW;
    }

    Migration& Params::M() {
        return * _migr;
    }

    void Params::set_transW_matrix(bool flag) {
        _transW_matrix = flag;
        if (flag) alloc_nalleles(_nalleles); // alloc the matrix if needed
    }

    bool Params::get_transW_matrix() const {
        return _transW_matrix;
    }

    void Params::set_transW_pair(unsigned int i, unsigned int j, double value) {
        _transW[i][i] -= _transW[i][j];
        _transW[i][i] += value;
        _transW[i][j] = value;
    }

    double Params::get_transW_pair(unsigned int i, unsigned int j) const {
        if (_transW_matrix) return _transW[i][j];
        return 1.;
    }

    double Params::get_transW_row(unsigned int i) const {
        if (_transW_matrix) return _transW[i][i];
        return _nalleles - 1;
    }

    void Params::addChange(Event * e) {
        _num_changes++;
        _base_change->insert(e);
        if (e->event_type() == Event::delayed) {
            _num_DSchanges++;
            _num_DSchanges_cache++;
        }
    }

    void Params::clearChanges() {
        _num_changes = 0;
        _num_DSchanges = 0;
        _num_DSchanges_cache = 0;
        _base_change->disconnect();
    }

    Event * Params::firstChange() {
        return _base_change->next();
    }

    double Params::nextChangeDate() const {
        return _cur_change->next() != NULL ? _cur_change->next()->date() : UNDEF;
    }

    void Params::nextChangeDo(Coalesce * coal) {
        _cur_change->next()->perform(this, coal);
        if (_cur_change->next()->event_type() == Event::delayed) _num_DSchanges--;
        _cur_change = _cur_change->next();
    }

    unsigned int Params::numChanges() const {
        return _num_changes;
    }

    unsigned int Params::nDSChanges() const {
        return _num_DSchanges;
    }

    double Params::lastChange(unsigned int pop) const {
        return _lastChange[pop];
    }

    void Params::restore() {
        if (_num_changes > 0) {
            _migr->restore();
            _recomb = _recomb_cache;
            for (unsigned int i=0; i<_npop; i++) {
                _popsize[i] = _popsize_cache[i];
                _growthrate[i] = _growthrate_cache[i];
                _selfing[i] = _selfing_cache[i];
                _lastChange[i] = 0.;
            }
            _cur_change = _base_change;
            _num_DSchanges = _num_DSchanges_cache;
        }
    }

    void Params::_init(unsigned int npop, double migr) {
        _recomb = 0.0;
        _recomb_cache = 0.0;

        _migr = NULL;

        // initialize, alloc and set /population arrays

        _npop_r = 0;
        _n1 = NULL;
        _n2 = NULL;
        _popsize = NULL;
        _popsize_cache = NULL;
        _growthrate = NULL;
        _growthrate_cache = NULL;
        _selfing = NULL;
        _selfing_cache = NULL;
        _lastChange = NULL;

        alloc_npop(npop);

        for (unsigned int i=0; i<npop; i++) {
            _n1[i] = 0;
            _n2[i] = 0;
            _popsize[i] = 1.0;
            _popsize_cache[i] = 1.0;
            _growthrate[i] = 0.0;
            _growthrate_cache[i] = 0.0;
            _selfing[i] = 0.0;
            _selfing_cache[i] = 0.0;
            _lastChange[i] = 0.0;
        }

        // set variables

        _max_iter = 100000;
        _theta = 0.0;
        _mutmodel = KAM;
        _nalleles = 2;
        _nalleles_r = 0;
        _random_start_allele = false;
        _fixed = 0;
        _TPMproba = 0.5;
        _TPMparam = 0.5;
        _outgroup = MISSING;

        // initialize /allele array

        _transW_matrix = false;
        _transW = NULL;

        // initialize /site arrays & variables

        _nsites = 0;
        _nsites_r = 0;
        _sitePos = NULL;
        _siteW = NULL;
        _totalSiteW = 0.0;

        // initialize /change arrays & variables

        _base_change = new (std::nothrow) Event(Event::none, UNDEF);
        if (_base_change == NULL) throw EGGMEM;
        _num_changes = 0;
        _cur_change = _base_change;
        _num_DSchanges = 0;
        _num_DSchanges_cache = 0;

        // initialize migration matrix

        _migr = new (std::nothrow) Migration(npop, migr);
        if (_migr == NULL) throw EGGMEM;
    }

    void Params::alloc_npop(unsigned int size) {

        if (size > _npop_r) {

            _n1 = (unsigned int*) realloc(_n1, size * sizeof(unsigned int));
            if (!_n1) throw EGGMEM;

            _n2 = (unsigned int*) realloc(_n2, size * sizeof(unsigned int));
            if (!_n2) throw EGGMEM;

            _popsize = (double*) realloc(_popsize, size * sizeof(double));
            if (!_popsize) throw EGGMEM;

            _popsize_cache = (double*) realloc(_popsize_cache, size * sizeof(double));
            if (!_popsize_cache) throw EGGMEM;

            _growthrate = (double*) realloc(_growthrate, size * sizeof(double));
            if (!_growthrate) throw EGGMEM;

            _growthrate_cache = (double*) realloc(_growthrate_cache, size * sizeof(double));
            if (!_growthrate_cache) throw EGGMEM;

            _selfing = (double*) realloc(_selfing, size * sizeof(double));
            if (!_n1) throw EGGMEM;

            _selfing_cache = (double*) realloc(_selfing_cache, size * sizeof(double));
            if (!_selfing_cache) throw EGGMEM;

            _lastChange = (double*) realloc(_lastChange, size * sizeof(double));
            if (!_lastChange) throw EGGMEM;

            _npop_r = size;
        }
        _npop = size;
    }

    void Params::alloc_nalleles(unsigned int size) {

        // allocate the table if needed

        if (size > _nalleles_r) {

            // first dimension

            _transW = (double **) realloc(_transW, size * sizeof(double *));
            if (!_transW) throw EGGMEM;

            // second dimension for reused cells

            for (unsigned int i=0; i<_nalleles_r; i++) {
                _transW[i] = (double *) realloc(_transW[i], size * sizeof(double));
                if (!_transW[i]) throw EGGMEM;
            }

            // second dimension for new cells

            for (unsigned int i=_nalleles_r; i<size; i++) {
                _transW[i] = (double *) malloc(size * sizeof(double));
                if (!_transW[i]) throw EGGMEM;
            }

            // set reserved size

            _nalleles_r = size;
        }

        // record operation size

        _nalleles = size;

        // initialize values whichever reused or not

        for (unsigned int i=0; i<_nalleles; i++) {
            for (unsigned int j=0; j<_nalleles; j++) {
                if (i==j) _transW[i][i] = _nalleles - 1.0;
                else _transW[i][j] = 1.0;
            }
        }
    }

    void Params::alloc_nsites(unsigned int size) {

        if (size > _nsites_r) {

            _sitePos = (double*) realloc(_sitePos, size * sizeof(double));
            if (!_sitePos) throw EGGMEM;

            _siteW = (double*) realloc(_siteW, size * sizeof(double));
            if (!_siteW) throw EGGMEM;

            _nsites_r = size;
        }

        for (unsigned int i=_nsites; i<size; i++) _siteW[i] = 1.;
        _totalSiteW = size;
        _nsites = size;
    }

    void Params::validate() const {

        // check main parameters

        if (_theta < 0.0) throw EggArgumentValueError("invalid theta (must be >= 0)");
        if (_theta > 0.0 && _fixed > 0) throw EggArgumentValueError("it is not allowed to set both theta and fixed to non-null values");
        if (_recomb < 0.0) throw EggArgumentValueError("invalid recombination rate (must be >= 0)"); 

        // check migration matrix consistency and values

        if (_npop != _migr->n()) throw EggArgumentValueError("invalid migration matrix (number of populations)"); 
        for (unsigned int i=0; i<_npop; i++) {
            for (unsigned int j=0; j<_npop; j++) {
                if (i==j) continue;
                if (_migr->get_pair(i, j) < 0.0)  throw EggArgumentValueError("invalid migration rate (must be >= 0)");
            }
        }

        // check population parameters

        for (unsigned int i=0; i<_npop; i++) {
            if (_popsize[i] <= 0.0) throw EggArgumentValueError("invalid population size (must be > 0)");
            if (_selfing[i] < 0.0) throw EggArgumentValueError("invalid selfing rate (must be >= 0)");
            if (_selfing[i] > 1.0) throw EggArgumentValueError("invalid selfing rate (must be <= 1)");
        }

        // check TPM parameters

        if (_TPMproba < 0.0) throw EggArgumentValueError("invalid TPM probability (must be >= 0");
        if (_TPMproba > 1.0) throw EggArgumentValueError("invalid TPM probability (must be <= 1");
        if (_TPMproba < 0.0) throw EggArgumentValueError("invalid TPM shape parameter (must be >= 0");
        if (_TPMproba > 1.0) throw EggArgumentValueError("invalid TPM shape parameter (must be <= 1");

        // check transition weights

        if (_transW_matrix) {
            for (unsigned int i=0; i<_nalleles; i++) {
                for (unsigned int j=0; j<_nalleles; j++) {
                    if (i==j) continue;
                    if (_transW[i][j] < 0.0) throw EggArgumentValueError("invalid transition weight (must be >= 0)");
                }
            }
        }

        // check site weights

        double tot = 0.;
        for (unsigned int i=0; i<_nsites; i++) {
            if (_siteW[i] < 0.0) throw EggArgumentValueError("invalid site weight (must be >= 0)");
            tot += _siteW[i];
        }
        if (tot != _totalSiteW) throw EggArgumentValueError("invalid site weights (sum does not match)");

        // check site positions

        if (_nsites > 0) {
            if (_sitePos[0] < 0.0) throw EggArgumentValueError("invalid site position (must be >= 0)");
            if (_sitePos[_nsites - 1] > 1.0) throw EggArgumentValueError("invalid site position (must be <= 1)");
        }
        for (unsigned int i=1; i<_nsites; i++) {
            if (_sitePos[i] <= _sitePos[i-1]) throw EggArgumentValueError("invalid site position (must be > previous)");
        }

        // check changes

        double t = 0.;
        Event * c = _base_change->next();

        while (c != NULL) {

            if (c->date() < t) throw EggArgumentValueError("invalid change: invalid date (must be >= 0 and >= previous)");

            if ((c->event_type() == Event::change_N || c->event_type() == Event::change_Mp ||
                c->event_type() == Event::change_G || c->event_type() == Event::change_s ||
                c->event_type() == Event::delayed ||  c->event_type() == Event::admixture ||
                c->event_type() == Event::bottleneck) && (c->get_index() >= _npop)) throw EggArgumentValueError("invalid change: invalid population index (too large)");

            if (c->event_type() == Event::bottleneck && c->get_param() < 0.0) throw EggArgumentValueError("invalid change: invalid bottleneck strength (must be >= 0)");

            if (c->event_type() == Event::delayed && c->get_number1() + c->get_number2() == 0) throw EggArgumentValueError("invalid change: delayed sample with no samples");

            if (c->event_type() == Event::change_Mp) {
                if (c->get_param() < 0) throw EggArgumentValueError("invalid change: invalid pairwise migration rate (must be >= 0)");
                if (c->get_index() == c->get_dest()) throw EggArgumentValueError("invalid change: invalid population indices (must be different)");
                if (c->get_dest() >= _npop) throw EggArgumentValueError("invalid change: invalid population index (too large)");
            }

            if (c->event_type() == Event::change_M && c->get_param() < 0) throw EggArgumentValueError("invalid change: invalid migration rate (must be >= 0)");

            if (c->event_type() == Event::change_s) {
                if (c->get_param() < 0.0) throw EggArgumentValueError("invalid change: invalid selfing rate (must be >= 0");
                if (c->get_param() > 1.0) throw EggArgumentValueError("invalid change: invalid selfing rate (must be <= 1");
            }

            if (c->event_type() == Event::change_N && c->get_param() < 0.0) throw EggArgumentValueError("invalid change: invalid population size (must be >=0)");

            if (c->event_type() == Event::admixture) {
                if (c->get_param() < 0.0) throw EggArgumentValueError("invalid change: invalid migration probability (must be >= 0");
                if (c->get_param() > 1.0) throw EggArgumentValueError("invalid change: invalid migration probability (must be <= 1");
                if (c->get_dest() >= _npop) throw EggArgumentValueError("invalid change: invalid population index (too large)");
            }

            if (c->event_type() == Event::change_R && c->get_param() < 0.0) throw EggArgumentValueError("invalid change: invalid recombination rate (must be >= 0)");

            c = c->next();
        }
    }

    void Params::get_nsam(unsigned int& ns, unsigned int& no) {
        ns = 0;
        no = 0;
        for (unsigned int i=0; i<_npop; i++) {
            if (i == _outgroup) no += _n1[i] + _n2[i] * 2;
            else ns += _n1[i] + _n2[i] * 2;
        }
        Event * cur = _base_change->next();
        while (cur != NULL) {
            if (cur->event_type() == Event::delayed) {
                if (cur->get_dest() == _outgroup) no += cur->get_number1() + cur->get_number2() * 2;
                else ns += cur->get_number1() + cur->get_number2() * 2;
            }
            cur = cur->next();
        }
    }
}
