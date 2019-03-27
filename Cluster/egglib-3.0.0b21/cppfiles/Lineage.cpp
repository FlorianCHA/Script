/*
    Copyright 2012 Stéphane De Mita, Mathieu Siol
    
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
#include "egglib.hpp"
#include "Lineage.hpp"

namespace egglib {
    
    // constructor to given size - no initialization
    Lineage::Lineage(unsigned int ntrees) {
        _ntrees = 0;
        _reserved = 0;
        nodeMapping = NULL;
        _cov = 0.;
        reset(ntrees);
    }
    
    
    // destructor
    Lineage::~Lineage() {
        if (_reserved) {
            free(nodeMapping);
        }
    }
    
    
    // add a tree to the array making the connection between the nodes and the trees
    void Lineage::addTree(unsigned int node, double cov) {
        alloc(_ntrees+1);
        nodeMapping[_ntrees-1] = node;
        _cov += cov;
    }
    
    
    // get a node index
    unsigned int Lineage::get_node(unsigned int index) const {
        return nodeMapping[index];
    }
    
    
    // set a Node
    void Lineage::set_node(unsigned int index, unsigned int node, double cov) {
        nodeMapping[index] = node;
        _cov += cov;
    }


   
    // reset
    void Lineage::reset(unsigned int ntrees) {
        alloc(ntrees);
        _cov = 0.;
    }

    // allocation or reallocation (if needed)
    void Lineage::alloc(unsigned int ntrees) {
        if (ntrees > _reserved) {
            nodeMapping = (unsigned int *) realloc(nodeMapping, ntrees * sizeof(unsigned int));
            if (!nodeMapping) throw EGGMEM;
            _reserved = ntrees;
        }
        _ntrees = ntrees;
    }
    
    // get coverage
    double Lineage::cov() const {
        return _cov;
    }
}
