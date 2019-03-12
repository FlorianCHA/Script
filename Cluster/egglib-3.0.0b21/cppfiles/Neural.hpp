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

#ifndef EGGLIB_NNET_HPP
#define EGGLIB_NNET_HPP

#include "egglib.hpp"
#include "Random.hpp"

namespace egglib {

   /** \brief Classes for neural networks
    *
    * The classes in this namespace are designed to train and use neural
    * networks. The main class is Network. Refer to this class's
    * documentation for more details.
    *
    */
    namespace nnet {

       /** \brief Activation functions */
        enum ActivationFunction {
            linear,     ///< linear function: f(x) = x
            log,        ///< logistic function: f(x) = 1 / [1+e^(-x)]
            tanh        ///< hyperbolic tangent function: f(x) = [1-e^(-2x)] / [1+e^(-2x)]
        };

        //************************************************************//
        
       /** \brief %Neuron, that is a node of a neural network
        *
        * \ingroup core
        *
        * Header: <egglib-cpp/Neural.hpp>
        *
        */
        class Neuron {

            private:

                ActivationFunction _fun;            // default: linear
                unsigned int _num_input;
                double * _input;
                double * _weight;
                double * _change_cache;
                double * _delta;
                unsigned int _cache_input;
                double _output;
                double _deriv;

               /** \brief It makes no sense to copy neurons */
                Neuron(const Neuron& n) {}

               /** \brief It makes no sense to copy neurons */
                Neuron& operator=(const Neuron& n) { return *this; }

            public:

               /** \brief Create a neuron
                *
                * The neuron is created naked, empty and unusable. The
                * user must use config() before doing anything with it.
                * To reuse a neuron after it has been used (e.g. for
                * training again the same network), it is required to
                * call config() again.
                * 
                */
                Neuron();

               /** \brief Delete a neuron
                *
                */
                virtual ~Neuron();

               /** \brief Set up the neuron
                *
                * The object is reset, but previously data may not have
                * been reinitialized.
                *
                * \param fun function used for activation of this
                *       neuron.
                * \param num_input number of incoming connexions (either
                *      neurons of the previous layer or input variables)
                *      of this neuron.
                *
                */
                void config(ActivationFunction fun, unsigned int num_input);

               /** \brief Set an input value
                *
                * Load a value for a given incoming connexion. Until
                * changed, the value will be remembered. Changing the
                * value does not update the output value of this neuron.
                *
                */
                void set_input(unsigned int index, double value);

               /** \brief Set a weight
                *
                * Load the value of the weight to be applied to a given
                * incoming connexion. Until changed, the value will be
                * remembered. Changing the value does not update the
                * output value of this neuron.
                *
                */
                void set_weight(unsigned int index, double value);

               /** \brief Get a weight
                *
                * Get the current valueo of the weight applied to a
                * given incoming connexion.
                *
                */
                double get_weight(unsigned int index) const;

               /** \brief Activate the neuron
                *
                * Process all incoming connexions, apply weights and the
                * activation function to generate the output.
                *
                */
                void activate();

               /** \brief Collect the output
                *
                * Get the output value. Calling this method does not
                * update the output if any input data has changed (make
                * sure to call the activate() method for this).
                *
                */
                double get_output() const;

               /** \brief Propagate error for a given weight
                *
                * This method can only be used with neurons that have
                * loaded input and have been activated. After the error
                * has been propagated, the delta value can be obtained
                * using delta(), typically for propagating errors to the
                * previous layer. The delta value is the value passed as
                * val to this method multiplied by the value derivative
                * of the activation function at the current input value
                * of this neuron. Note that the neuron's weight are not
                * modified by this method, but only when update() is
                * called.
                *
                * \param index index of one of the weights of this
                *   neuron.
                * \param val propagated error for this weight (for an
                *   output neuron: this neuron's error times the output
                *   value of the corresponding neuron of the hidden
                *   layer; for a hidden neuron: the sum of delta values
                *   for all output neurons, weighted by the
                *   corresponding weights.
                *
                */
                void propagate(unsigned int index, double val);

               /** \brief Get the delta value for a given weight
                *
                * The propagate() method must have been called, which
                * itself requires that the neuron had been previously
                * loaded with all needed data and activated. This method
                * returns the change of the specified weight value based
                * on the propagated error. The delta value is computed
                * as f'(I) * E[index] where f' is the derivative
                * of the activation function, I is the input value for
                * this neuron and E[index] is the error associated to
                * the weight in question.
                *
                */
                double delta(unsigned int index) const;

               /** \brief Update all weights
                *
                * The propagate() method must have been called for all
                * weights. This method applies all delta values.
                *
                * \param rate learning rate.
                * \param bound absolute limit value: weights are bound
                *   to the range [-bound, +bound].
                * \param momentum proportion of the previous weight
                *   change to apply to the new change.
                *
                */
                void update(double rate, double bound, double momentum);
        };

        //************************************************************//

       /** \brief Holds input and output data for training neural networks
        *
        * When creating a data set, the user must first set the number
        * of patterns with num_patterns(), of input variables with
        * num_input() and of output variables with num_output(). Only
        * then it is possible to load data using input() and output().
        * Make sure to load every declared slots as data are not
        * initialized.
        *
        * \note In case a same object must be reused will smaller values
        * of num_input and/or num_output more larger num_patterns, it is
        * more efficient to call the methods num_input() and
        * num_output() before num_patterns().
        * 
        * \ingroup core
        *
        * Header: <egglib-cpp/Neural.hpp>
        *
        */
        class Data {

            private:

                unsigned int _num_patterns;
                unsigned int _num_input;
                unsigned int _num_output;
                double ** _input;         // num_patterns * num_input
                double ** _output;        // num_patterns * num_output
                unsigned int _cache_patterns;
                unsigned int * _cache_input;    // num_patterns
                unsigned int * _cache_output;   // num_patterns
                bool _normalized;
                double * _mean_input;     // num_input
                double * _std_input;      // num_input
                double * _mean_output;    // num_output
                double * _std_output;     // num_output
                unsigned int _cache_input_meanstd;
                unsigned int _cache_output_meanstd;

                void copy(const Data& source);
                void init();
                void free();

            public:

               /** \brief Constructor */
                Data();
            
               /** \brief Copy constructor */
                Data(const Data& src);

               /** \brief Copy assignment operator */
                Data& operator=(const Data& src);

               /** \brief Destructor */
                virtual ~Data();

               /** \brief Set number of patterns
                *
                * Invalidate all means and standard deviations if the
                * data was previously normalized.
                *
                */
                void set_num_patterns(unsigned int num);

               /** \brief Get number of patterns */
                unsigned int get_num_patterns() const;

               /** \brief Set number of input variables
                *
                * Invalidate all means and standard deviations if the
                * data was previously normalized.
                *
                */
                void set_num_input(unsigned int num);

               /** \brief Get number of input variables */
                unsigned int get_num_input() const;

               /** \brief Set number of output variables
                *
                * Invalidate all means and standard deviations if the
                * data was previously normalized.
                *
                */
                void set_num_output(unsigned int num);

               /** \brief Get number of output variables */
                unsigned int get_num_output() const;

               /** \brief Load an input data item */
                void set_input(unsigned int pattern, unsigned int variable, double value);

               /** \brief Load an output data item */
                void set_output(unsigned int pattern, unsigned int variable, double value);

               /** \brief Get an input data item */
                double get_input(unsigned int pattern, unsigned int variable) const;

               /** \brief Get an output data item */
                double get_output(unsigned int pattern, unsigned int variable) const;

               /** \brief Normalize data for all input and output variables
                *
                * All data are modified permanently and the average and
                * standard deviation for each input and output variable
                * are saved and remain accessible using the methods
                * mean_input(), mean_output(), std_input() and
                * std_output(), until normalize() is called again or the
                * number of input or output variables or the number of
                * patterns is modified.
                *
                */
                void normalize();

               /** \brief Normalize data for an input variable
                *
                * All data are modified permanently. The passed mean and
                * standard deviation are not saved.
                *
                */
                void normalize_input(unsigned int index, double mean, double std);

               /** \brief Normalize data for an output variable
                *
                * All data are modified permanently. The passed mean and
                * standard deviation are not saved.
                *
                */
                void normalize_output(unsigned int index, double mean, double std);

               /** \brief Unnormalize data for an output variable
                *
                * All data are modified permanently. The passed mean and
                * standard deviation are not saved.
                *
                */
                void unnormalize_output(unsigned int index, double mean, double std);

               /** \brief Get the mean for an input variable
                *
                * If data have been normalized using the method
                * normalize(), this method returns the mean of a given
                * input variable in the original data.
                *
                */
                double mean_input(unsigned int index) const;

               /** \brief Get the standard deviation for an input variable
                *
                * If data have been normalized using the method
                * normalize(), this method returns the standard
                * deviation of a given input variable in the original
                * data.
                *
                */
                double std_input(unsigned int index) const;

               /** \brief Get the mean for an output variable
                *
                * If data have been normalized using the method
                * normalize(), this method returns the mean of a given
                * output variable in the original data.
                *
                */
                double mean_output(unsigned int index) const;

               /** \brief Get the standard deviation for an output variable
                *
                * If data have been normalized using the method
                * normalize(), this method returns the standard
                * deviation of a given output variable in the original
                * data.
                *
                */
                double std_output(unsigned int index) const;
        };

        //************************************************************//

       /** \brief Training and prediction with neural networks
        *
        * This classes implements the back propagation algorithm for
        * training neural network.
        *
        * The constructor of Neural generates a bare, unusable instance.
        * It is required set up the network with the setup() method and
        * initialize weights (normally using the init_weights() method,
        * or by setting manually all weights) before calling train().
        * More information is available in the documentation of the
        * setup() method.
        *
        * \ingroup core
        *
        * Header: <egglib-cpp/Neural.hpp>
        *
        */
        class Network {

            private:

                Neuron ** _neurons;
                unsigned int _cache_neurons;
                unsigned int _num_input;
                unsigned int _num_hidden;
                unsigned int _num_output;
                unsigned int _num_patterns;

                const Data * _training_data;
                const Data * _testing_data;
                unsigned int _training_start;
                unsigned int _training_stop;
                unsigned int _testing_start;
                unsigned int _testing_stop;
                ActivationFunction _fun1;
                ActivationFunction _fun2;
                Random * _random;
                double _rate1;
                double _rate2;
                double _bound;
                double _momentum;

                double * _predictions;
                double * _error;
                unsigned int _cache_output;

                double _training_error;
                double _testing_error;
                double _pattern_error;

                unsigned int _num_weights;
                unsigned int _cache_weights;
                double * _change_cache;

                unsigned int _num_iter;

                void backpropagate();
                void compute_error();

               /** \brief It is not possible to copy neural networks */
                Network(const Network& n) {}

               /** \brief It is not possible to copy neural networks */
                Network& operator=(const Network& n) { return *this; }

            public:

               /** \brief Constructor */
                Network();

               /** \brief Destructor */
                virtual ~Network();

               /** \brief Setup the network
                *
                * Upon call to this method, the network is set up based
                * on the number of input and output variables of the
                * provided data set, and the specified number of neurons
                * in the hidden layer. Note that all neurons (from both
                * the hidden and output layers) are automatically
                * connected to an additional neuron generating a
                * constant input of 1.0 (the bias).
                *
                * Loading a training data set (using the train_data) is
                * required. The first and last-plus-one indexes must be
                * passed to indicate which patterns must be processed
                * for training. Note that the pattern corresponding to
                * the train_stop argument is not included. It is
                * possible to use a train_stop value larger than the
                * number of patterns in train_data. To use all patterns,
                * use train_start=0 and train_stop=egglib::MAX.
                *
                * Loading a test data set is not mandatory but
                * advisable. The test data set is not used for training
                * but allows to evaluate the predictive ability of the
                * network. To skip this option, pass any data set as
                * test_data (for example, the same as train_data) and
                * set test_start >= test_stop.
                *
                * \param training_data a training data set with at least
                *   one input variable, at least one output variable, at
                *   least one pattern and all data loaded. The object
                *   passed must not modified until training is finished.
                * \param training_start index of the first pattern to
                *   process for training.
                * \param training_stop index of the pattern immediately
                *   after the last pattern to process for training.
                * \param testing_data a data set to use for testing the
                *   predictive ability of the network (not used for
                *   fitting). The data set must have the same number of
                *   input and output variables as train_data (preferably
                *   it will the same object, this non-overlapping ranges
                *   of patterns to process).
                * \param testing_start index of the first pattern to
                *   process for testing.
                * \param testing_stop index of the pattern immediately
                *   after the last pattern to process for test.
                * \param num_hidden number of neurons in the hidden
                *   layer of the network.
                * \param fun1 activation function to use for neurons of
                *   the hidden layer.
                * \param fun2 activation function to use for output
                *   neurons.
                * \param rate1 rate for weights from input to hidden
                *   layer.
                * \param rate2 rate for weights form hidden layer to
                *   output.
                * \param bound extreme value (for both signs) for weight
                *   values.
                * \param momentum proportion of the previous weight
                *   change to apply to all weight changes (use 0.0 to
                *   skip momentum).
                * \param random Random object to use for generating
                *   random numbers (only used for initial weights).
                *
                */
                void setup(const Data& training_data,
                    unsigned int training_start, unsigned int training_stop,
                    const Data& testing_data,
                    unsigned int testing_start, unsigned int testing_stop,
                    unsigned int num_hidden,
                    ActivationFunction fun1, ActivationFunction fun2,
                    double rate1, double rate2, double bound,
                    double momentum, Random * random);

               /** \brief Initialize weights
                *
                * It is required to call this method before calling
                * training, otherwise the weights have undefined values.
                * All weights are initialized to random values from a
                * uniform distribution in the range [-X, X] with X is
                * the range_input argument for all weights connecting
                * input variables to neurons of the hidden layer, and
                * the range_output arguments for all weights connecting
                * neurons of the hidden layer to the output neurons.
                *
                * \param range_input bound for first-level weights.
                * \param range_output bound for second-level weights.
                *
                */
                void init_weights(double range_input = 0.01, double range_output = 0.01);

               /** \brief Train the neural network
                *
                * Train the network for a fixed number of iterations.
                * This method must be called iteratively until a given
                * stop criterion is fulfilled. After each call to
                * train(), it is possible to access the the current
                * values of weights, the error for training data and, if
                * testing data have been loaded, the error for testing
                * data.
                *
                */
                void train(unsigned int num_iter);

               /** \brief Use the neural network to predict output
                *
                * The network must have been trained using a Data
                * instance with the same number of input and output
                * variables. This method generates an array of predicted
                * values based on the current values of the weight
                * (after training) that can be accessed using the
                * prediction() method. If the compute_error argument is
                * true, the error is accessed using error().
                *
                * \param data data set with the correct number of input
                *   (always) and output variables (unless compute_errors
                *   is false).
                * \param pattern pattern to process.
                * \param compute_error if false, don't compute errors,
                *   and output variables are not considered.
                *
                */ 
                void predict(const Data& data, unsigned int pattern, bool compute_error);

               /** \brief Get a predicted output
                *
                * The predict() method must have been called. Get the
                * predicted value for one of the output variables.
                * Warning: training (using the method train()) modifies
                * the predicted values and invalidates the results of
                * the method predict().
                *
                * \param output_var index of the output variable (must
                *   be smalled than the number of output variables
                *   defined by both the Data instances used for training
                *   and for prediction).
                * 
                */
                double prediction(unsigned int output_var) const;

               /** \brief Get total error for training data
                *
                * The error is computed sqrt(sum((pred[i]-obs[i])^2))
                * where sqrt is the square root function, sum is the
                * sum over all output nodes, pred[i] is the predicted
                * value provided by output neuron i and obs[i] is the
                * observed value for output variable i. Computed using
                * the training data only. Requires train() but also
                * modified if predict() is called with
                * compute_error=true.
                *
                */
                double training_error() const;

               /** \brief Get total error for testing data
                *
                * As training_error() but for the testing data set. Only
                * valid if testing data has been passed to setup().
                *
                */
                double testing_error() const;

               /** \brief Value of a weight to a hidden neuron
                *
                * Access to the current value of the weight of the
                * connexion of a hidden neuron to an input variable or
                * to the bias neuron. The bias weight's index is equal
                * to the number of input variables.
                *
                * \param i index of the hidden neuron.
                * \param j index of the input variable.
                *
                */
                double weight_hidden(unsigned int i, unsigned int j) const;

               /** \brief Value of a weight to an output neuron
                *
                * Access to the current value of the weight of the
                * connexion of an output neuron to hidden neuron or to
                * the bias neuron. The bias weight's index is equal to
                * the number of hidden neurons.
                *
                * \param i index of the output neuron.
                * \param j index of the hidden neuron.
                *
                */
                double weight_output(unsigned int i, unsigned int j) const;

               /** \brief Set the value of a weight to a hidden neuron
                *
                * Set the value of the weight of the connexion of a
                * hidden neuron to an input variable or to the bias
                * neuron. The bias weight's index is equal to the number
                * of input variables.
                *
                * \param i index of the hidden neuron.
                * \param j index of the input variable.
                * \param value weight value.
                *
                */
                void weight_hidden(unsigned int i, unsigned int j, double value);

               /** \brief Set the value of a weight to an output neuron
                *
                * Set the value of the weight of the connexion of an
                * output variable to a hidden neuron or to the bias
                * neuron. The bias weight's index is equal to the number
                * of hidden neurons.
                *
                * \param i index of the output neuron.
                * \param j index of the hidden neuron.
                * \param value weight value.
                *
                */
                void weight_output(unsigned int i, unsigned int j, double value);

               /** \brief Get the number of training iterations
                *
                * The internal counter is reset by setup().
                * 
                */
                unsigned int num_iter() const;
        };
    }
}

#endif
