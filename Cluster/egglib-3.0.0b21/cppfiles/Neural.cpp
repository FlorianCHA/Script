/*
    Copyright 2013 Stéphane De Mita, Mathieu Siol

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
#include "Neural.hpp"
#include <cmath>
#include <cstdlib>

namespace egglib {
    namespace nnet {

        Neuron::Neuron() {
            _fun = linear;
            _num_input = 0;
            _cache_input = 0;
            _input = NULL;
            _weight = NULL;
            _change_cache = NULL;
            _delta = NULL;
            _output = 0.0;
            _deriv = 0.0;
        }

        //////

        Neuron::~Neuron() {
            if (_input) free(_input);
            if (_weight) free(_weight);
            if (_change_cache) free(_change_cache);
            if (_delta) free(_delta);
        }

        //////

        void Neuron::config(ActivationFunction fun, unsigned int num_input) {
            _fun = fun;
            if (num_input > _cache_input) {
                _input = (double *) realloc(_input, num_input * sizeof(double));
                if (!_input) throw EGGMEM;
                _weight = (double *) realloc(_weight, num_input * sizeof(double));
                if (!_weight) throw EGGMEM;
                _change_cache = (double *) realloc(_change_cache, num_input * sizeof(double));
                if (!_change_cache) throw EGGMEM;
                _delta = (double *) realloc(_delta, num_input * sizeof(double));
                if (!_delta) throw EGGMEM;
                _cache_input = num_input;
            }
            _num_input = num_input;

            for (unsigned int i=0; i<_num_input; i++) _change_cache[i] = 0.0;
        }

        //////

        void Neuron::set_input(unsigned int index, double value) {
            _input[index] = value;
        }

        //////

        void Neuron::set_weight(unsigned int index, double value) {
            _weight[index] = value;
        }

        //////

        double Neuron::get_weight(unsigned int index) const {
            return _weight[index];
        }

        //////

        void Neuron::activate() {
            _output = 0.0;
            for (unsigned int i=0; i<_num_input; i++) {
                _output += _input[i] * _weight[i];
            }
            if (_fun == log) {
                _output = 1. / (1 + exp(-_output));
                _deriv = _output * (1 - _output);
            }
            else {
                if (_fun == tanh) {
                    _output = ::tanh(_output);
                    _deriv = 1.0 - _output * _output;
                }
                else {
                    _deriv = 1.0;
                }
            }
        }

        //////

        void Neuron::propagate(unsigned int index, double val) {
            _delta[index] = _deriv * val;
        }

        //////

        double Neuron::delta(unsigned int index) const {
            return _delta[index];
        }

        //////

        void Neuron::update(double rate, double bound, double momentum) {

            double change;

            for (unsigned int i=0; i<_num_input; i++) {
                change = momentum * _change_cache[i] + rate * _delta[i] * _input[i];
                _change_cache[i] = change;
                _weight[i] += change;
                if (_weight[i] < -bound) _weight[i] = -bound;
                else if (_weight[i] > bound) _weight[i] = bound;
            }
        }

        //////

        double Neuron::get_output() const {
            return _output;
        }

        //************************************************************//

        void Data::copy(const Data& src) {
            set_num_input(src._num_input);
            set_num_output(src._num_output);
            set_num_patterns(src._num_patterns);

            for (unsigned int i=0; i<_num_patterns; i++) {
                for (unsigned int j=0; j<_num_input; j++) {
                    _input[i][j] = src._input[i][j];
                }
                for (unsigned int j=0; j<_num_output; j++) {
                    _output[i][j] = src._output[i][j];
                }
            }

            if (src._normalized) {

                _normalized = true;

                if (_num_input > _cache_input_meanstd) {
                    _mean_input = (double *) realloc(_mean_input, _num_input * sizeof(double));
                    if (!_mean_input) throw EGGMEM;
                    _std_input = (double *) realloc(_std_input, _num_input * sizeof(double));
                    if (!_std_input) throw EGGMEM;
                    _cache_input_meanstd = _num_input;
                }

                if (_num_output > _cache_output_meanstd) {
                    _mean_output = (double *) realloc(_mean_output, _num_output * sizeof(double));
                    if (!_mean_output) throw EGGMEM;
                    _std_output = (double *) realloc(_std_output, _num_output * sizeof(double));
                    if (!_std_output) throw EGGMEM;
                    _cache_output_meanstd = _num_output;
                }

                for (unsigned int i=0; i<_num_input; i++) {
                    _mean_input[i] = src._mean_input[i];
                    _std_input[i] = src._std_input[i];
                }

                for (unsigned int i=0; i<_num_output; i++) {
                    _mean_output[i] = src._mean_output[i];
                    _std_output[i] = src._std_output[i];
                }
            }
        }

        //////

        void Data::init() {
            _num_patterns = 0;
            _num_input = 0;
            _num_output = 0;
            _cache_patterns = 0;
            _cache_input = 0;
            _cache_output = 0;
            _input = NULL;
            _output = NULL;
            _normalized = false;
            _mean_input = NULL;
            _std_input = NULL;
            _mean_output = NULL;
            _std_output = NULL;
            _cache_input_meanstd = 0;
            _cache_output_meanstd = 0;
        }

        //////

        void Data::free() {
            for (unsigned int i=0; i<_cache_patterns; i++) {
                if (_input[i]) ::free(_input[i]);
                if (_output[i]) ::free(_output[i]);
            }
            if (_input) ::free(_input);
            if (_output) ::free(_output);
            if (_cache_input) ::free(_cache_input);
            if (_cache_output) ::free(_cache_output);
            if (_mean_input) ::free(_mean_input);
            if (_std_input) ::free(_std_input);
            if (_mean_output) ::free(_mean_output);
            if (_std_output) ::free(_std_output);
        }
            
        //////

        Data::Data() {
            init();
        }

        //////

        Data::Data(const Data& src) {
            init();
            copy(src);
        }

        //////

        Data& Data::operator=(const Data& src) {
            copy(src);
            return * this;
        }

        //////

        Data::~Data() {
            free();
        }

        //////

        void Data::set_num_patterns(unsigned int num) {

            _normalized = false;

            // get value

            _num_patterns = num;

            // allocate new lines

            if (_num_patterns > _cache_patterns) {

                _input = (double **) realloc(_input, _num_patterns * sizeof(double *));
                if (!_input) throw EGGMEM;
                    
                _output = (double **) realloc(_output, _num_patterns * sizeof(double *));
                if (!_output) throw EGGMEM;

                _cache_input = (unsigned int *) realloc(_cache_input, _num_patterns * sizeof(unsigned int));
                if (!_cache_input) throw EGGMEM;
                    
                _cache_output = (unsigned int *) realloc(_cache_output, _num_patterns * sizeof(unsigned int));
                if (!_cache_output) throw EGGMEM;

                // initialize lines

                for (unsigned int i=_cache_patterns; i<_num_patterns; i++) {
                    _input[i] = NULL;
                    _cache_input[i] = 0;
                    _output[i] = NULL;
                    _cache_output[i] = 0;
                }

                _cache_patterns = _num_patterns;
            }

            // ensure all used lines have enough length

            for (unsigned int i=0; i<_num_patterns; i++) {

                if (_num_input > _cache_input[i]) {
                    _input[i] = (double *) realloc(_input[i], _num_input * sizeof(double));
                    if (!_input[i]) throw EGGMEM;
                    _cache_input[i] = _num_input;
                }

                if (_num_output > _cache_output[i]) {
                    _output[i] = (double *) realloc(_output[i], _num_output * sizeof(double));
                    if (!_output[i]) throw EGGMEM;
                    _cache_output[i] = _num_output;
                }
            }
        }

        //////

        unsigned int Data::get_num_patterns() const {
            return _num_patterns;
        }

        //////

        void Data::set_num_input(unsigned int num) {

            _normalized = false;

            for (unsigned int i=0; i<_num_patterns; i++) {
                if (num > _cache_input[i]) {
                    _input[i] = (double *) realloc(_input[i], num * sizeof(double));
                    if (!_input[i]) throw EGGMEM;
                    _cache_input[i] = num;
                }
            }
            _num_input = num;
        }

        //////

        unsigned int Data::get_num_input() const {
            return _num_input;
        }

        //////

        void Data::set_num_output(unsigned int num) {

            _normalized = false;

            for (unsigned int i=0; i<_num_patterns; i++) {
                if (num > _cache_output[i]) {
                    _output[i] = (double *) realloc(_output[i], num * sizeof(double));
                    if (!_output[i]) throw EGGMEM;
                    _cache_output[i] = num;
                }
            }
            _num_output = num;
        }

        //////

        unsigned int Data::get_num_output() const {
            return _num_output;
        }

        //////

        void Data::set_input(unsigned int pattern, unsigned int variable, double value) {
            _input[pattern][variable] = value;
        }

        //////

        void Data::set_output(unsigned int pattern, unsigned int variable, double value) {
            _output[pattern][variable] = value;
        }

        //////

        double Data::get_input(unsigned int pattern, unsigned int variable) const {
            return _input[pattern][variable];
        }

        //////

        double Data::get_output(unsigned int pattern, unsigned int variable) const {
            return _output[pattern][variable];
        }

        //////

        void Data::normalize() {

            // alloc for mean/std table if needed

            _normalized = true;

            if (_num_input > _cache_input_meanstd) {
                _mean_input = (double *) realloc(_mean_input, _num_input * sizeof(double));
                if (!_mean_input) throw EGGMEM;
                _std_input = (double *) realloc(_std_input, _num_input * sizeof(double));
                if (!_std_input) throw EGGMEM;
                _cache_input_meanstd = _num_input;
            }

            if (_num_output > _cache_output_meanstd) {
                _mean_output = (double *) realloc(_mean_output, _num_output * sizeof(double));
                if (!_mean_output) throw EGGMEM;
                _std_output = (double *) realloc(_std_output, _num_output * sizeof(double));
                if (!_std_output) throw EGGMEM;
                _cache_output_meanstd = _num_output;
            }

            // compute average and std deviation for each variable

            double sum;
            double sum2;

            for (unsigned int i=0; i<_num_input; i++) {
                sum = 0.0;
                sum2 = 0.0;
                for (unsigned int j=0; j<_num_patterns; j++) {
                    sum += _input[j][i];
                    sum2 += _input[j][i] * _input[j][i];
                }
                _mean_input[i] = sum / _num_patterns;
                _std_input[i] = sqrt(sum2 / _num_patterns - _mean_input[i] * _mean_input[i]);
            }

            for (unsigned int i=0; i<_num_output; i++) {
                sum = 0.0;
                sum2 = 0.0;
                for (unsigned int j=0; j<_num_patterns; j++) {
                    sum += _output[j][i];
                    sum2 += _output[j][i] * _output[j][i];
                }
                _mean_output[i] = sum / _num_patterns;
                _std_output[i] = sqrt(sum2 / _num_patterns - _mean_output[i] * _mean_output[i]);
            }

            // normalize data

            for (unsigned int i=0; i<_num_patterns; i++) {
                for (unsigned int j=0; j<_num_input; j++) {
                    if (_std_input[j] > 0.0) _input[i][j] = (_input[i][j] - _mean_input[j]) / _std_input[j];
                }
                for (unsigned int j=0; j<_num_output; j++) {
                    if (_std_output[j] > 0.0) _output[i][j] = (_output[i][j] - _mean_output[j]) / _std_output[j];
                }
            }
        }

        //////

        void Data::normalize_input(unsigned int index, double mean, double std) {
            if (std > 0.0) {
                for (unsigned int i=0; i<_num_patterns; i++) {
                    _input[i][index] = (_input[i][index] - mean) / std;
                }
            }
        }

        void Data::normalize_output(unsigned int index, double mean, double std) {
            if (std > 0.0) {
                for (unsigned int i=0; i<_num_patterns; i++) {
                    _output[i][index] = (_output[i][index] - mean) / std;
                }
            }
        }

        void Data::unnormalize_output(unsigned int index, double mean, double std) {
            for (unsigned int i=0; i<_num_patterns; i++) {
                _output[i][index] = _output[i][index] * std + mean;
            }
        }

        //////

        double Data::mean_input(unsigned int index) const {
            return _mean_input[index];
        }

        //////

        double Data::std_input(unsigned int index) const {
            return _std_input[index];
        }

        //////

        double Data::mean_output(unsigned int index) const {
            return _mean_output[index];
        }

        //////

        double Data::std_output(unsigned int index) const {
            return _std_output[index];
        }

        //************************************************************//

        Network::Network() {
            _neurons = NULL;
            _cache_neurons = 0;
            _num_input = 0;
            _num_hidden = 0;
            _num_output = 0;
            _num_patterns = 0;

            _training_data = NULL;
            _testing_data = NULL;
            _training_start = 0;
            _training_stop = MAX;
            _testing_start = 0;
            _testing_stop = 0;
            _fun1 = linear;
            _fun2 = linear;
            _random = NULL;
            _rate1 = 0.0;
            _rate2 = 0.0;
            _bound = 0.0;

            _predictions = NULL;
            _error = NULL;
            _cache_output = 0;

            _training_error = 0;
            _testing_error = 0;
            _pattern_error = 0;

            _num_iter = 0;
        }

        //////

        Network::~Network() {
            for (unsigned int i=0; i<_cache_neurons; i++) {
                if (_neurons[i]) delete _neurons[i];
            }
            if (_neurons) free(_neurons);
            if (_predictions) free(_predictions);
            if (_error) free(_error);
        }

        //////

        void Network::setup(const Data& training_data,
                    unsigned int training_start, unsigned int training_stop,
                    const Data& testing_data,
                    unsigned int testing_start, unsigned int testing_stop,
                    unsigned int num_hidden,
                    ActivationFunction fun1, ActivationFunction fun2,
                    double rate1, double rate2, double bound,
                    double momentum, Random * random) {

            // set training & testing data

            _training_data = & training_data;
            _training_start = training_start;
            _training_stop = training_stop < training_data.get_num_patterns() ? training_stop : training_data.get_num_patterns();

            _testing_data = & testing_data;
            _testing_start = testing_start;

            _testing_stop = testing_stop < testing_data.get_num_patterns() ? testing_stop : testing_data.get_num_patterns();

            // set network dimensions

            _num_input = training_data.get_num_input();
            _num_hidden = num_hidden;
            _num_output = training_data.get_num_output();
            _num_patterns = _training_stop - _training_start;
            unsigned int num_neurons = _num_hidden + _num_output;

            // claim neurons objects

            if (num_neurons > _cache_neurons) {
                _neurons = (Neuron **) realloc(_neurons, num_neurons * sizeof(Neuron *));
                if (!_neurons) throw EGGMEM;

                for (unsigned int i=_cache_neurons; i<num_neurons; i++) {
                    _neurons[i] = new (std::nothrow) Neuron;
                    if (!_neurons[i]) throw EGGMEM;
                }                
                _cache_neurons = num_neurons;
            }

            // set network connexions

            for (unsigned int i=0; i<_num_hidden; i++) _neurons[i]->config(fun1, _num_input + 1);               // + 1 for bias
            for (unsigned int i=0; i<_num_output; i++) _neurons[_num_hidden+i]->config(fun2, _num_hidden + 1);  // + 1 for bias

            // allocate output arrays

            if (_num_output > _cache_output) {
                _predictions = (double *) realloc(_predictions, _num_output * sizeof(double));
                if (!_predictions) throw EGGMEM;

                _error = (double *) realloc(_error, _num_output * sizeof(double));
                if (!_error) throw EGGMEM;

                _cache_output = _num_output;
            }

            // set parameters

            _num_iter = 0;
            _fun1 = fun1;
            _fun2 = fun2;
            _rate1 = rate1;
            _rate2 = rate2;
            _bound = bound;
            _momentum = momentum;
            _random = random;
        }

        //////

        void Network::init_weights(double range_input, double range_output) {

            for (unsigned int i=0; i<_num_hidden; i++) {
                for (unsigned int j=0; j<_num_input + 1; j++) { // + 1 for bias
                    _neurons[i]->set_weight(j, (_random->uniformcl() - 0.5) * 2 * range_input);
                }
            }

            for (unsigned int i=0; i<_num_output; i++) {
                for (unsigned int j=0; j<_num_hidden + 1; j++) { // + 1 for bias
                    _neurons[_num_hidden+i]->set_weight(j, (_random->uniformcl() - 0.5) * 2 * range_output);
                }
            }
        }

        //////

        void Network::train(unsigned int num_iter) {
            for (unsigned int i=0; i<num_iter; i++) {
                for (unsigned int j=0; j<_num_patterns; j++) {
                    unsigned int rnd = _training_start + _random->irand(_num_patterns);
                    predict(* _training_data, rnd, true);
                    backpropagate();
                }
            }
            compute_error();
            _num_iter += num_iter;
        }

        //////

        void Network::compute_error() {
            _training_error = 0.0;
            for (unsigned int i=_training_start; i<_training_stop; i++) {
                predict(* _training_data, i, true);
if (_pattern_error > 1000000000) exit(0);
                _training_error += _pattern_error;
            }
            _training_error = sqrt(_training_error);

            _testing_error = 0.0;
            for (unsigned int i=_testing_start; i<_testing_stop; i++) {
                predict(* _testing_data, i, true);
                _testing_error += _pattern_error;
            }
            _testing_error = sqrt(_testing_error);
        }

        //////

        void Network::predict(const Data& data, unsigned int pattern, bool compute_error) {

            // pass hidden layer

            for (unsigned int i=0; i<_num_hidden; i++) {
                for (unsigned int j=0; j<_num_input; j++) {
                    _neurons[i]->set_input(j, data.get_input(pattern, j));
                }
                _neurons[i]->set_input(_num_input, 1.0); // + 1 for bias
                _neurons[i]->activate();
            }

            // pass output layer and compute error

            if (compute_error) {
                _pattern_error = 0.0;
            }

            for (unsigned int i=0; i<_num_output; i++) {
                for (unsigned int j=0; j<_num_hidden; j++) {
                    _neurons[_num_hidden+i]->set_input(j, _neurons[j]->get_output());
                }
                _neurons[_num_hidden+i]->set_input(_num_hidden, 1.0); // + 1 for bias
                _neurons[_num_hidden+i]->activate();
                _predictions[i] = _neurons[_num_hidden+i]->get_output();
                if (compute_error) {
                    _error[i] = data.get_output(pattern, i) - _predictions[i];
                    _pattern_error += (_predictions[i] - data.get_output(pattern, i)) * (_predictions[i] - data.get_output(pattern, i));
                }
            }
        }

        //////

        void Network::backpropagate() {

            // propagate errors to the output layer

            for (unsigned int i=0; i<_num_output; i++) {
                for (unsigned int j=0; j<_num_hidden + 1; j++) { // + 1 for bias
                    _neurons[_num_hidden+i]->propagate(j, _error[i]);
                }
            }

            // propagate errors to the hidden layer

            double err;

            for (unsigned int i=0; i<_num_hidden; i++) {

                // get the total error received from the output layer

                err = 0.0;
                for (unsigned int j=0; j<_num_output; j++) {
                    err += _neurons[_num_hidden+j]->delta(i) * _neurons[_num_hidden+j]->get_weight(i);
                }

                // send this error to all weights connecting the input layer

                for (unsigned int j=0; j<_num_input + 1; j++) { // + 1 for bias
                    _neurons[i]->propagate(j, err);
                }
            }

            // apply all changes

            for (unsigned int i=0; i<_num_hidden; i++) _neurons[i]->update(_rate1, _bound, _momentum);
            for (unsigned int i=0; i<_num_output; i++) _neurons[_num_hidden+i]->update(_rate2, _bound, _momentum);
        }

        //////

        double Network::prediction(unsigned int output_var) const {
            return _predictions[output_var];
        }

        //////

        double Network::training_error() const {
            return _training_error;
        }

        //////

        double Network::testing_error() const {
            return _testing_error;
        }

        //////

        double Network::weight_hidden(unsigned int i, unsigned int j) const {
            return _neurons[i]->get_weight(j);
        }

        //////

        double Network::weight_output(unsigned int i, unsigned int j) const {
            return _neurons[_num_hidden+i]->get_weight(j);
        }

        //////

        void Network::weight_hidden(unsigned int i, unsigned int j, double value) {
            _neurons[i]->set_weight(j, value);
        }

        //////

        void Network::weight_output(unsigned int i, unsigned int j, double value) {
            _neurons[_num_hidden+i]->set_weight(j, value);
        }

        //////

        unsigned int Network::num_iter() const {
            return _num_iter;
        }
    }
}
