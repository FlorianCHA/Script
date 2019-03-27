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

#include <cstdlib>
#include <cstring>
#include <new>
#include "Window.hpp"

namespace egglib {

    WSite::WSite(WPool * p) {
        _pool = p;
        reset(1);
    }

    void WSite::reset(unsigned int pl) {
        _pos = UNKNOWN;
        _site.reset(pl);
    }

    void WSite::init() {
        _prev = NULL;
        _next = NULL;
    }

    WSite::~WSite() {
    }

    SiteHolder& WSite::site() {
        return _site;
    }

    void WSite::set_pos(unsigned int p) {
        _pos = p;
    }

    unsigned int WSite::get_pos() const {
        return _pos;
    }

    WSite * WSite::next() {
        return _next;
    }

    WSite * WSite::prev() {
        return _prev;
    }

    WSite * WSite::push_back(WSite * ws) {
        _next = ws;
        ws->_prev = this;
        ws->_next = NULL;
        return ws;
    }

    WSite * WSite::pop_front() {
        if (_next != NULL) _next->_prev = NULL;
        _pool->put(this);
        return _next;
    }

    WSite * WSite::pop_back() {
        if (_prev != NULL) _prev->_next = NULL;
        _pool->put(this);
        return _prev;
    }

    WPool::WPool() {
        _cache = NULL;
        _pool = NULL;
        _c_cache = 0;
        _c_pool = 0;
        _n_pool = 0;
    }

    WPool::~WPool() {
        for (unsigned int i=0; i<_c_cache; i++) {
            if (_cache[i]) delete _cache[i];
        }
        if (_cache) free(_cache);
        if (_pool) free(_pool);
    }

    WSite * WPool::get() {
        if (_n_pool == 0) {
            _c_cache++;
            _cache = (WSite **) realloc(_cache, _c_cache * sizeof(WSite *));
            if (!_cache) throw EGGMEM;
            _cache[_c_cache-1] = new(std::nothrow) WSite(this);
            if (!_cache[_c_cache-1]) throw EGGMEM;
            return _cache[_c_cache-1];
        }
        else {
            return _pool[--_n_pool];
        }
    }

    void WPool::put(WSite * p) {
        _n_pool++;
        if (_n_pool > _c_pool) {
            _pool = (WSite **) realloc(_pool, _n_pool * sizeof(WSite *));
            if (!_pool) throw EGGMEM;
            _c_pool = _n_pool;
        }
        _pool[_n_pool-1] = p;
    }

    VcfWindow::VcfWindow() {
        _num = 0;
        _chrom = (char *) malloc(1 * sizeof(char));
        if (!_chrom) throw EGGMEM;
        _c_chrom = 1;
        _cur_pos = UNKNOWN;
        _start_pos = UNKNOWN;
        _stop_pos = UNKNOWN;
        _first_site = NULL;
        _last_site = NULL;
        _first_pos = UNKNOWN;
        _last_pos = UNKNOWN;
        _good = false;
    }

    VcfWindow::~VcfWindow() {
        if (_chrom) free(_chrom);
    }

    void VcfWindow::setup(VcfParser& vcf, unsigned int wsize, unsigned int wstep, bool unit_bp,
                unsigned int start_pos, unsigned int stop_pos,
                unsigned int max_missing) {

        while (_first_site != NULL) _first_site = _first_site->pop_front();
        _last_site = NULL;

        _vcf = &vcf;
        _wsize = wsize;
        _wstep = wstep;
        _unit_bp = unit_bp;
        _max_missing = max_missing;
        _num = 0;
        _start_pos = start_pos;
        _stop_pos = stop_pos;
        if (unit_bp) _cur_pos = _start_pos;
        _first_site = NULL;
        _last_site = NULL;
        _first_pos = UNKNOWN;
        _last_pos = UNKNOWN;
        _good = vcf.good();

        // get chromosome name
        if (_good) {
            vcf.read();
            const char * chr = vcf.chromosome();
            if (strlen(chr) > _c_chrom) {
                _chrom = (char *) realloc(_chrom, (strlen(chr)+1) * sizeof(char));
                if (!_chrom) throw EGGMEM;
            }
            strcpy(_chrom, chr);
            vcf.unread();
        }
        else {
            _chrom[0] = '\0';
        }

        // advance until start position is met
        while (_read(_start_pos));
    }

    void VcfWindow::next_window() {

        if (_unit_bp) {
            _first_pos = _cur_pos;
            _last_pos = _cur_pos + _wsize - 1;
            if (_last_pos >= _stop_pos) {
                _last_pos = _stop_pos - 1;
                _good = false;
            }
            _cur_pos += _wstep;
            if (_cur_pos >= _stop_pos) {
                _good = false;
            }
        }

        // pop sites from start of window
        if (_unit_bp) {
            while (_first_site != NULL && _first_site->get_pos() < _first_pos) {
                _first_site = _first_site->pop_front();
                if (!_first_site) _last_site = NULL;
                _num--;
            }
            if (!_first_site) _last_site = NULL;
        }
        else {
            for (unsigned int i=0; i<_wstep && i<_wsize; i++) {
                if (_first_site == NULL) break;
                _first_site = _first_site->pop_front();
                if (!_first_site) _last_site = NULL;
                _num--;
            }
        }

        // add sites at the end of window
        if (_unit_bp) {
            while (_read(_last_pos+1)) _add();
        }
        else while (_num < _wsize && _read(_stop_pos)) {
            _add();
        }

        // set start/stop positions of current window
        if (!_unit_bp) {
            if (_num > 0) {
                _first_pos = _first_site->get_pos();
                _last_pos = _last_site->get_pos();
            }
            else {
                _first_pos = UNKNOWN;
                _last_pos = UNKNOWN;
            }

            // check that there is at least one site left
            if (_read(_stop_pos)) _vcf->unread();
            else _good = false;
        }

        // go to start of next window is step>size
        if (_unit_bp) {
            while (_read(_cur_pos));
        }
        else {
            for (unsigned int i=_wsize; i<_wstep; i++) {
                if (!_read(_stop_pos)) break;
            }
        }
    }

    bool VcfWindow::_read(unsigned int limit) {
        // no data at all
        if (!_vcf->good()) {
            _good = false;
            return false;
        }

        // read one
        _vcf->read();

        // is next chromosome
        if (strcmp(_vcf->chromosome(), _chrom)) {
            _vcf->unread();
            _good = false;
            return false;
        }

        // reached limit
        if (_vcf->position() >= limit) {
            _vcf->unread();
            return false;
        }
        return true;
    }

    void VcfWindow::_add() {
        if (!_last_site) {
            _last_site = _pool.get();
            _first_site = _last_site;
            _first_site->init();
        }
        else _last_site = _last_site->push_back(_pool.get());
        _last_site->reset(_vcf->ploidy());
        _last_site->set_pos(_vcf->position());

        if (!_last_site->site().process_vcf(*_vcf, 0, _vcf->num_samples(), _max_missing)) {
            _last_site = _last_site->pop_back();
            if (!_last_site) _first_site = NULL;
        }
        else {
            _num++;
        }
    }

    const char * VcfWindow::chromosome() const { return _chrom; }
    unsigned int VcfWindow::first_pos() const { return _first_pos; }
    unsigned int VcfWindow::last_pos() const { return _last_pos; }
    unsigned int VcfWindow::num_sites() const { return _num; }
    const WSite * VcfWindow::first_site() const { return _first_site; }
    const WSite * VcfWindow::last_site() const { return _last_site; }
    bool VcfWindow::good() const { return _good; }
    unsigned int VcfWindow::num_samples() const { return _vcf->num_samples(); }

    const WSite * VcfWindow::get_site(unsigned int idx) const {
        WSite * site = _first_site;
        for (unsigned int i=0; i<idx; i++) site = site->next();
        return site;
    }

    unsigned int VcfWindow::first_site_pos() const {
        return _first_site->get_pos();
    }

    unsigned int VcfWindow::last_site_pos() const {
        return _last_site->get_pos();
    }

//////////

    BedWindow::BedWindow(){
        init();
    }

    BedWindow::~BedWindow(){
        _VCF = NULL;
        _Bed = NULL;
        if(Site_list.size()!=0) clear_window();
        free(_chrom_window); 
        _indice = 0;
        b_size = 0;
        _start_pv = 0;
        _end_pv = 0;
        _miss_pv =0;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void BedWindow::set_start_window(unsigned int position) { start_w = position; }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void BedWindow::set_end_window(unsigned int position) { end_w = position; }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void BedWindow::set_chromosome_window(const char * chrom) {
        int s_chrom=0;
        s_chrom=strlen(chrom)+2;
        _chrom_window= (char*) realloc (_chrom_window, s_chrom *sizeof(char));
        strcpy(_chrom_window, chrom);    
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    unsigned int BedWindow::get_end_w() { return end_w; }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    unsigned int BedWindow::get_start_w() { return start_w; }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    char *  BedWindow::get_chrom_w() { return _chrom_window; }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<egglib::SiteHolder*> BedWindow::get_site_list() { return Site_list;}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void BedWindow::init(){
        _VCF = NULL;
        _Bed = NULL;
        _chrom_window = (char *) malloc(1 * sizeof(char));
        if (!_chrom_window) throw EGGMEM;
        _chrom_window[0] = '\0';
        _indice = 0;
        start_w=0;
        end_w=0;
    }        
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool BedWindow::good(){
        return (_indice < _Bed->get_bed_list().size());    
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
    void BedWindow::configure(egglib::VcfParser * VCF, BedParser * Bed, unsigned int start_pv, unsigned int end_pv, unsigned int miss_pv){
        _VCF = VCF;
        _Bed = Bed;
        b_size = _Bed->get_bed_list().size();
        _start_pv = start_pv;
        _end_pv = end_pv;
        _miss_pv = miss_pv;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void BedWindow::get_start_variant_vcf(unsigned int indice){
        _VCF->rewind();
        while(_VCF->good()){ //true
            _VCF->read();
            if((strcmp(_VCF->chromosome(),(_Bed->get_bed_list()[indice].chromosome))==0) && (_VCF->position() >= _Bed->get_bed_list()[indice].start_pos) && (_VCF->position() < _Bed->get_bed_list()[indice].end_pos)){
                if(_VCF->has_data() && _VCF->has_GT() && (_VCF->num_missing_GT() <= _miss_pv)){               
                start_w = _VCF->position();
                    set_chromosome_window(_VCF->chromosome());
                    _VCF->unread();    
                    break;
                }
            }    
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void BedWindow::clear_window(){
        for (unsigned int i = 0; i < Site_list.size(); i++ ){           
            delete Site_list[i];    
        }
        Site_list.clear(); 
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool BedWindow::before_bound() const{
        return ((strcmp(_chrom_window,_VCF->chromosome())==0) && ((_VCF->position()) <= _end_bound));
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void BedWindow::get_sliding_window(){    

        if(Site_list.size()!=0) clear_window();
        int unsigned i=1;
        bool b_site = 0;
        while(_VCF->good()){
            _VCF->read();
            if(!before_bound()){     
                _VCF->unread();
                break;
            }

            _ploidy_pv = _VCF->ploidy();

            egglib::SiteHolder* site = NULL;
            site = new(std::nothrow) egglib::SiteHolder(_ploidy_pv);
            if (!site) throw EGGMEM;
            b_site=(*site).process_vcf(*_VCF, _start_pv, _end_pv, _miss_pv);

            if (b_site) {
                Site_list.push_back(site);
            end_w = _VCF->position();
                ++i;
            }
            else {

                delete site;
            }
            
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////            

    void BedWindow::next(){
        if(!good()) {
            throw EggArgumentValueError("The Bed Window reached the last data of the bed file.");
        }

        get_start_variant_vcf(_indice);
        _end_bound = _Bed->get_bed_list()[_indice].end_pos;
        get_sliding_window();
        _indice ++;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    egglib::SiteHolder& BedWindow::get_site(unsigned int i){
        if(i > Site_list.size()) throw EggArgumentValueError("you can't pass in parameter a indice greater than the total number of sites in the list of sites");
        return * Site_list[i];    
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int unsigned BedWindow::get_n_site(){ return Site_list.size(); }


    void BedWindow::at(unsigned int indice){
        get_start_variant_vcf(indice);
        _end_bound = _Bed->get_bed_list()[indice].end_pos;
        get_sliding_window();
    }
}
