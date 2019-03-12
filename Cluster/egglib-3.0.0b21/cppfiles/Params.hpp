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

#ifndef EGGLIB_PARAMS_HPP
#define EGGLIB_PARAMS_HPP

#include "egglib.hpp"
#include <string>

namespace egglib {

    class Coalesce;
    class Params;

   /** \brief Arbitrary migration matrix type
    *
    * \ingroup coalesce
    *
    * This class represents an arbitrary migration model where any
    * pairwise can be specified independently.
    *
    * Header: <egglib-cpp/Params.hpp>
    *
    */
    class Migration {

        public:

           /** \brief Constructor
            *
            * \param n number of populations.
            * \param M migration rate.
            *
            * The migration rate argument is per population. By default,
            * all pairwise migration rates are M/(n-1).
            *
            */
            Migration(unsigned int n, double M);

           /** \brief Destructor
            *
            */
            ~Migration();

           /** \brief Number of populations
            *
            */
            unsigned int n() const;

           /** \brief Set all migration rates
            *
            * The migration rate argument is per population. Change all
            * pairwise migration rates to the value M/(n-1). This method
            * is definitive - the previous value is lost.
            *
            */
            void set_all(double M);

           /** \brief Set all migration rates (recursive)
            *
            * The migration rate argument is per population. Change all
            * pairwise migration rates to the value M/(n-1). This method
            * is reversible by the restore method.
            *
            */
            void set_all_R(double M);

           /** \brief Set migration rates of a given source population
            *
            * Change all pairwise migration rates originating from the
            * specified population. The migration rate argument is per
            * population. Change all pairwise migration rates to the
            * value M/(n-1). population should be smaller than n. If
            * this condition is not met, the behaviour is not defined.
            * This method is definitive - the previous value is lost.
            *
            */
            void set_row(unsigned int i, double M);

           /** \brief Set migration rates from a pop (recursive)
            *
            * Change all pairwise migration rates originating from the
            * specified population. The migration rate argument is per
            * population. Change all pairwise migration rates to the
            * value M/(n-1). population should be smaller than n. If
            * this condition is not met, the behaviour is not defined.
            * This method is reversible by the restore method. In
            * M/(n-1).
            *
            */
            void set_row_R(unsigned int i, double M);

           /** \brief Set a given pairwise migration rate
            *
            * The migration rate argument is the pairwise value. The
            * diagonal sum is automatically updated. i and j should both
            * be smaller than n, and should not be equal. If any of
            * these conditions is not met, the behaviour is not defined
            * and is likely to be unpleasant. The user is prompted to be
            * particularly vigilant regarding the fact that i and j
            * must be different. This method is definitive - the
            * previous value is lost.
            *
            */
            void set_pair(unsigned int i, unsigned int j, double m);

           /** \brief Set a pairwise migration rate (reversible)
            *
            * The migration rate argument is the pairwise value. The
            * diagonal sum is automatically updated. i and j should both
            * be smaller than n, and should not be equal. If any of
            * these conditions is not met, the behaviour is not defined
            * and is likely to be unpleasant. The user is prompted to be
            * particularly vigilant regarding the fact that i and j
            * must be different. This method is definitive - the
            * previous value is lost.  This method is reversible by the
            * restore method.
            *
            */
            void set_pair_R(unsigned int i, unsigned int j, double m);

           /** \brief Get a population migration rate
            *
            * Return the sum of pairwise migration rates from this
            * population. If i>=n1*n2, the behaviour is undefined.
            *
            */
            double get_row(unsigned int i) const;

           /** \brief Get a pairwise migration rate
            *
            * Return the pairwise migration from population i to
            * population j. Both arguments should be smaller than n
            * and i should be different of j. It is not guaranteed
            * that all errors will properly generate an exception.
            *
            */
            double get_pair(unsigned int i, unsigned int j) const;

           /** \brief Restore to initial state
            *
            *
            */
            void restore();

        protected:

            unsigned int npop;
            unsigned int npop_reserved;
            double **matrix; // cache is included (all values doubled)

        private:

           /** \brief Copy constructor not available */
            Migration(const Migration& src) {}

           /** \brief Copy assignment operator not available*/
            Migration& operator=(const Migration& src) { return *this; }
    };

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

   /** \brief Generic class for all historical events
    *
    * \ingroup coalesce
    *
    * All types of events are implemented by this class. To create an
    * event of a given type, pass the date and the appropriate flag to
    * the constructor and then set the correct parameters. The table
    * below gives the meaning of parameters for all types of events
    * (blank when the parameter is undefined).
    *
    * \verbatim
    * | type       | date | param             | index       | dest     | number1         | number2         |
    * | :--------- | :--: | :---------------: | :---------: | :------: | :-------------: | :-------------: |
    * | none       | +    |                   |             |          |                 |                 |
    * | change_N   | +    | new size          | population* |          |                 |                 |
    * | change_M   | +    | new popwise rate  |             |          |                 |                 |
    * | change_Mp  | +    | new pairwise rate | population  | dest pop |                 |                 |
    * | change_G   | +    | new rate          | population* |          |                 |                 |
    * | change_s   | +    | new rate          | population* |          |                 |                 |
    * | change_R   | +    | new rate          |             |          |                 |                 |
    * | delayed    | +    |                   | population  | label    | haploid samples | diploid samples |
    * | admixture  | +    | migr proba        | source pop  | dest pop |                 |                 |
    * | bottleneck | +    | strength          | population* |          |                 |                 |
    * \endverbatim
    *
    * (*) for those events (and them only), it is possible to use
    * ``egglib::MAX`` as population index to specify "all populations".
    *
    */
    class Event {

        public:

            /// Type of the event
            enum Type {
                none,        ///< no event (use only as list head)
                change_N,    ///< change population size
                change_M,    ///< change all migration rates
                change_Mp,   ///< change a pairwise migration rate
                change_G,    ///< change exponential growth rate/decline
                change_s,    ///< change autofertilization rate
                change_R,    ///< change recombination rate
                bottleneck,  ///< cause immediate coalescences in a population
                admixture,   ///< move lineages from a population to another
                delayed      ///< delayed sample
            };

        private:

            /// Copy not available
            Event(const Event& src) {}

            /// Copy not available
            Event& operator=(const Event& src) {return * this;}

            void _insert_up(Event * event);

            static const double _small;
            Event * _prev;
            Event * _next;
            Type _type;
            double _date;
            unsigned int _index;
            unsigned int _dest;
            double _param;
            unsigned int _number1;
            unsigned int _number2;


        public:

           /** \brief Constructor
            *
            * Depending on the object type, it is required to specify
            * at least one other parameter in addition to the date.
            *
            */
            Event(Type type, double date);

            /// \brief Copy all members (except links) from a source
            void copy(const Event& src);

            /// \brief Get type
            Type event_type() const;

            /// \brief Get date
            double date() const;

            /// \brief Get parameter
            double get_param() const;

            /// \brief Get locus or main population index
            unsigned int get_index() const;

            /// \brief Get destination population index or label for delayed events
            unsigned int get_dest() const;

            /// \brief Get first number
            unsigned int get_number1() const;

            /// \brief Get second number
            unsigned int get_number2() const;

           /** \brief Apply the event
            *
            * The parameter set is modified and (for some events only)
            * the coalesce instance as well. It is illegal to call this
            * on an event of none type.
            *
            */
            void perform(Params * param, Coalesce * coal);

           /** \brief Insert an event in the chain
            *
            * The current previous and next links value of the passed
            * instance (if any) will be ignored and overwritten. This
            * method should be called only on the head item of a chain.
            * It is illegal to load a none type.
            *
            */
            void insert(Event * event);

           /** \brief Disconnect all events in chain (this and down)
            *
            */
            void disconnect();

           /** \brief Next event (NULL if end)
            *
            */
            Event * next();

           /** \brief previous event (NULL if end)
            *
            */
            Event * prev();

           /** \brief Change the date and update position
            *
            * The event is moved either up or down the chain of events
            * until its date is at the right location.
            *
            * It is possible to call this method on an unconnected
            * instance.
            *
            */
            void move(double date);

           /** \brief Set the event index
            *
            * Sets the population (for event types that support it) that
            * should be affected by the event. All events except
            * recombination rate changes and global migration rate
            * changes  must have a value for this parameter (except the
            * type "none" that must not be used as a genuine event). The
            * value egglib::MAX may be used (but only for some kinds of
            * events; see class description) to specify "all
            * populations".
            *
            */
            void set_index(unsigned int i);

           /** \brief Set the event parameter
            *
            * Most types have a value for this, but the meaning varies:
            * it is generally a new rate (migration, recombination,
            * selfing) or a new size (population size) and also the
            * strength for bottleneck and instant emigration for
            * admixture.
            *
            */
            void set_param(double p);

           /** \brief Set the event destination population
            *
            * For admixture and pairwise migration rate change events:
            * sets the index of the destination population. For delayed
            * sample events: set the value used as population label for
            * this sample.
            *
            */
            void set_dest(unsigned int d);

           /** \brief Set the number of haploid samples
            *
            * Number of samples with one sampled allele (for delayed
            * sample events).
            *
            */
            void set_number1(unsigned int n);

           /** \brief Set the number of diploid samples
            *
            * Number of samples with both sampled alleles (for delayed
            * sample events).
            *
            */
            void set_number2(unsigned int n);
    };

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

   /** \brief Coalescence simulation parameters holder
    *
    * \ingroup coalesce
    *
    * Holds parameters for a coalescence simulation. Default values are:
    *  \li k: 1
    *  \li L: 0
    *  \li R: 0.
    *  \li theta: 0.
    *  \li fixed: 0
    *  \li mutmodel: KAM
    *  \li TPMproba: 0.5
    *  \li TPMparam: 0.5
    *  \li K: 2
    *  \li random_start_allele: false
    *  \li n1: 0 for any population
    *  \li n2: 0 for any population
    *  \li N: 1. for any population
    *  \li G: 0. for any population
    *  \li s: 0. for any population
    *  \li sitePos: not initialized (beware!)
    *  \li siteW: 1. for any site
    *  \li M: all rates at 0.0
    *  \li transW_matrix: false
    *  \li transW: 1. for any combination
    *  \li changes: none
    *  \li max_iter: 100000
    *
    * See the corresponding setter methods for more details on the
    * signification of parameters. Note that every population has an
    * independent size, growth/decline rate and selfing rate (in
    * addition to numbers of samples, obviously).
    *
    * Header: <egglib-cpp/Params.hpp>
    *
    */
    class Params {

        public:

           /** \brief Specification for mutation models
            *
            * The possible values are:
            *  \li KAM: fixed finite number of alleles.
            *  \li IAM: infinite number of alleles.
            *  \li SMM: stepwise mutation model.
            *  \li TPM: two phase mutation model.
            *
            * All models are constrained by the mutation parameter
            * (theta or, if theta is 0, fixed). TPM is the only model
            * that requires additional parameters (TPMproba and TPMparam
            * see below).
            *
            * Nucleotide sequences can be simulated using the KAM with
            * K=2 or 4 (the latter should be used with models with a
            * finite number of sites, to control homoplasy). When using
            * an infinite site model, each site cannot experience more
            * than one mutation. Therefore, it makes little sense to use
            * any other mutation model than KAM with K=2.
            *
            * In the KAM, each mutation cause a transition to one of the
            * other alleles. Transition biases can be controlled using
            * the tranW matrix (if tranW_matrix is true). The allele
            * values are integers in the range 0, 1, ..., K-1.
            *
            * In the IAM, each mutation generate a new allele. Alleles
            * values are integers in the range 0, 1, ..., S-1 where S is
            * the number of mutations having occurred for a particular
            * site. The allele value indicates the rank of occurrence
            * of the mutation, and should not be taken as having any
            * biological signification.
            *
            * In the SMM, each mutation change the allele value by a
            * step of -1 or +1 (with equal probabilities). Alleles
            * values can therefore have any negative or positive values.
            *
            * In the TPM, each mutation change the allele value by a
            * step of -n or +n (with equal probabilities). n is either
            * 1 or drawn from a geometric distribution. TPMproba gives
            * the probability that n follows a geometric distribution
            * (that is, of not following a SMM) and TPMparam gives the
            * shape parameter of the geometric distribution. The GSM
            * (generalized stepwise model) can be obtained using the TPM
            * by setting TPMproba to 1.0.
            *
            */
            enum MutationModel { KAM, IAM, SMM, TPM };

           /** \brief Default constructor
            *
            * Parameters are initialized to default values. The number
            * of populations is 1.
            *
            */
            Params();

           /** \brief Standard constructor
            *
            * Except the number of populations, the parameters are
            * initialized to default values. The number of populations
            * must be at least 1. All pairwise migration rates can be
            * specified later (the optional argument is the
            * per-population migration rate to be applied to all
            * populations.
            *
            */
            Params(unsigned int npop, double migr=0.0);

           /** \brief Destructor
            *
            */
            ~Params();

           /** \brief String summary of the parameter values
            *
            * The string is generated each time this method is called.
            *
            */
            std::string summary() const;

           /** \brief Get the number of populations
            *
            */
            unsigned int k() const;

           /** \brief Get the number of alleles
            *
            */
            unsigned int get_K() const;

           /** \brief Set the number of alleles
            *
            * The value must be at least 2. If the transition matrix
            * flag is on, calling this method will reset all weights to
            * default value (all weights to 1, sum to number of alleles
            * minus one). Also if the new value is smaller than or equal
            * to the previous one.
            *
            */
            void set_K(unsigned int value);

            /// Maximum number of iterations
            unsigned long get_max_iter() const;

            /// Maximum number of iterations
            void set_max_iter(unsigned long);

           /** \brief Get random start allele boolean
            *
            */
            bool get_random_start_allele() const;

           /** \brief Set random start allele boolean
            *
            * If true, the coalescent simulator will use as start allele
            * a random value in the range [0, K-1]. If false, it will
            * use 0.
            *
            */
            void set_random_start_allele(bool value);

           /** \brief Set the recombination rate
            *
            * The value must be at least 0.
            *
            */
            void set_R(double value);

           /** \brief Set the recombination rate (reversible)
            *
            */
            void R_R(double value);

           /** \brief Get the recombination rate
            *
            */
            double get_R() const;

           /** \brief Set the mutation rate
            *
            * If theta is 0, the value of fixed is considered. The value
            * cannot be less than 0.
            *
            */
            void set_theta(double value);

           /** \brief Get the mutation rate
            *
            * If theta is 0, the value of fixed is considered.
            *
            */
            double get_theta() const;

           /** \brief Set the fixed number of mutations
            *
            * This parameter is ignored if theta is different of 0.
            * The value can be 0.
            *
            */
            void set_fixed(unsigned int value);

           /** \brief Get the fixed number of mutations
            *
            * This parameter is ignored if theta is different of 0.
            *
            */
            unsigned int get_fixed() const;

           /** \brief Set the mutation model
            *
            * Detailed information is provided within the documentation
            * of the enum type MutationModel.
            *
            */
            void set_mutmodel(MutationModel value);

           /** \brief Get the mutation model value
            *
            * Detailed information is provided within the documentation
            * of the enum type MutationModel.
            *
            */
            MutationModel get_mutmodel() const;

           /** \brief Set the TPM probability parameter
            *
            * This parameter is ignored if the mutation model is any
            * other than TPM. The probability gives the probability that
            * mutation steps are drawn from a geometric distribution
            * (otherwise, they are set to 1). The value must be at least
            * 0. and at most 1.
            *
            */
            void set_TPMproba(double value);

           /** \brief Get the TPM probability parameter
            *
            * This parameter is ignored if the mutation model is any
            * other than TPM. The probability gives the probability that
            * mutation steps are drawn from a geometric distribution
            * (otherwise, they are set to 1).
            *
            */
            double get_TPMproba() const;

           /** \brief Set the shape parameter of the TPM distribution
            *
            * This parameter is ignored if the mutation model is any
            * other than TPM. This parameter sets the shape of the
            * geometric distribution from which mutation steps are drawn
            * with probability given by TPMproba. The passed value
            * should be >= 0 and <= 1.
            *
            */
            void set_TPMparam(double value);

           /** \brief Get the shape parameter of the TPM distribution
            *
            * This parameter is ignored if the mutation model is any
            * other than TPM. This parameter sets the shape of the
            * geometric distribution from which mutation steps are drawn
            * with probability given by TPMproba.
            *
            */
            double get_TPMparam() const;

           /** \brief Set the number of single samples in a population
            *
            * \param pop population index.
            * \param value number of single samples.
            *
            * The population index must be smaller than the last number
            * of populations specified. Single samples are individuals
            * in which only one allele is sampled. If the selfing rate
            * is 0, there is no difference between 2n single samples and
            * n double samples (or any combination summing up to a total
            * of 2n alleles).
            *
            */
            void set_n1(unsigned int pop, unsigned int value);

           /** \brief Get the number of single samples in a population
            *
            */
            unsigned int get_n1(unsigned int pop) const;

           /** \brief Set the number of double samples in a population
            *
            * \param pop population index.
            * \param value number of double samples.
            *
            * The population index must be smaller than the last number
            * of populations specified. Double samples are individuals
            * in which both alleles are sampled. They therefore actually
            * represent two samples (which are correlated whenever the
            * selfing rate is larger than 0). If the selfing rate is 0,
            * there is no difference between 2n single samples and n
            * double samples (or any combination summing up to a total
            * of 2n alleles).
            *
            */
            void set_n2(unsigned int pop, unsigned int value);

           /** \brief Get the number of double samples in a population
            *
            */
            unsigned int get_n2(unsigned int pop) const;

           /** \brief Set a population relative size
            *
            * \param pop population index.
            * \param value population size.
            *
            * The population index must be smaller than the last number
            * of population specified. The population size must be
            * larger than 0.
            *
            */
            void set_N(unsigned int pop, double value);

           /** \brief Set a population relative size (reversible)
            *
            * \param pop population index.
            * \param value population size.
            * \param t time of the change.
            *
            */
            void N_R(unsigned int pop, double value, double t);

           /** \brief Get a population relative size
            *
            */
            double get_N(unsigned int pop) const;

           /** \brief Set a population relative growth/decline rate
            *
            * \param pop population index.
            * \param value growth/decline rate.
            *
            * The population index must be smaller than the last number
            * of population specified. If the rate is > 0, then the
            * population is experiencing an exponential growth (was
            * smaller in the past). If the rate is < 0, then the
            * population is experiencing an exponential decline (was
            * bigger in the past). Caution is required when using
            * exponential population decline. Going back in time,
            * population size can grow to excessive size and cause
            * infinite coalescence time (this will cause a runtime
            * error). To prevent this, always ensure that population
            * decline stops at some point in the past and that
            * coalescence of remaining lineages can occur.
            *
            */
            void set_G(unsigned int pop, double value);

           /** \brief Set a population relative growth/decline rate(reversible)
            *
            * \param pop population index.
            * \param value growth/decline rate.
            * \param t time of the change.
            *
            */
            void G_R(unsigned int pop, double value, double t);

           /** \brief Get a population growth/decline rate
            *
            */
            double get_G(unsigned int pop) const;

           /** \brief Set a population selfing rate
            *
            * \param pop population index.
            * \param value selfing rate.
            *
            * The population index must be smaller than the last number
            * of population specified. The selfing rate must be >= 0 and
            * <= 1.
            *
            */
            void set_s(unsigned int pop, double value);

           /** \brief Set a population selfing rate (reversible)
            *
            * \param pop population index.
            * \param value selfing rate.
            * \param t time of the change.
            *
            */
            void s_R(unsigned int pop, double value, double t);

           /** \brief Get a population selfing rate
            *
            */
            double get_s(unsigned int pop) const;

           /** \brief Set label identifying outgroup samples.
            *
            * MISSING for none.
            *
            */
            void set_outgroup(unsigned int lbl);

           /** \brief Get label identifying outgroup samples.
            *
            * MISSING for none.
            *
            */
            unsigned int get_outgroup() const;

           /** \brief Get number of mutable sites
            *
            * A value of 0 denotes an infinite number of sites.
            *
            */
            unsigned int get_L() const;

           /** \brief Set number of mutable sites
            *
            * A value of 0 denotes an infinite number of sites. The site
            * weights are initialized to 1, but the site positions are
            * not initialized.
            *
            * \param value new number.
            *
            */
            void set_L(unsigned int value);

           /** \brief Automatically set site positions
            *
            * Set all site positions such as they use all possible span
            * over the simulated chromosome. If there only one site, it
            * is set at position 0.5. If there are two sites, there are
            * set a positions 0.0 and 1.0. If there are three sites, at
            * 0.0, 0.5, 1.0, and so on.
            *
            */
            void autoSitePos();

           /** \brief Set a site position
            *
            * \param site site index.
            * \param value site position.
            *
            * The site index must be smaller than L. The site position
            * must be >= 0 and <= 1. It is assumed that sites are ranked
            * by increased site position, but this is not enforced.
            *
            */
            void set_sitePos(unsigned int site, double value);

           /** \brief Get a site position
            *
            */
            double get_sitePos(unsigned int site) const;

           /** \brief Set a site weight
            *
            * Site weights give the relative mutation probability of
            * sites. The higher the site weight, the more likely that a
            * mutation will hit this site. Note that the sum of weights
            * is not required to sum up to one. Weights don't affect the
            * total mutation rate. Once a mutation occurs, they
            * determine which site is affected.
            *
            * \param site site index.
            * \param value site weight.
            *
            * The site index must be smaller than L. The site weight can
            * be any strictly positive value.
            *
            */
            void set_siteW(unsigned int site, double value);

           /** \brief Get a site weight
            *
            */
            double get_siteW(unsigned int site) const;

           /** \brief Get the sum of site weights
            *
            */
            double totalSiteW() const;

           /** \brief Get the migration matrix
            *
            */
            Migration & M();

           /** \brief %Allele transition matrix flag
            *
            * If true, transition weights between alleles are taken
            * from a matrix (values must be set using method transW).
            * Setting the flag to true will reset all weights to default
            * values: all weights to 1, sum to number of alleles
            * minus one. Also if the new value is smaller than or equal
            * to the previous one (even if the flag was already true).
            *
            */
            void set_transW_matrix(bool flag);

           /** \brief %Allele transition matrix flag
            *
            * If true, transition weights between alleles are taken
            * from a matrix (values must be set using method transW).
            *
            */
            bool get_transW_matrix() const;

           /** \brief Set an allele transition matrix entry
            *
            * Set transition weights from allele i to allele j to value.
            * Transition weights are not required to sum up to 1. If the
            * sum is different of 1, it will not affect the overall
            * mutation rate. If a mutation occur and if the current
            * allele is the one given by i, the value will determine the
            * relative probability of mutating to the allele given by j.
            * The value must be strictly larger than 0. Both i and j
            * must be smaller than the number of alleles (K), and they
            * must be different to each other. It is not allowed to call
            * this method if transW_matrix() is false (and it will
            * likely result in a crash).
            *
            */
            void set_transW_pair(unsigned int i, unsigned int j, double value);

           /** \brief Get an allele transition matrix entry
            *
            * This method should not be used to access to diagonal
            * items. If transW_matrix() is false, this method always
            * return 1.
            *
            */
            double get_transW_pair(unsigned int i, unsigned int j) const;

           /** \brief Get transition weight sum for an allele
            *
            * If transW_matrix() is false, this method always return
            * the number of alleles minus one (the weight is always 1).
            *
            */
            double get_transW_row(unsigned int i) const;

           /** \brief Add a demographic change
            *
            * The argument must be an instance of Event (but not of the
            * none type). Each of the types represents a given type of
            * demographic change (meaning is broad). All changes have a
            * date parameter, specifying at what date the change should
            * occur (in backward, coalescent time). The unit is 4N
            * generations, where N is the number of diploi individuals
            * of a population. Multiple changes can be specified (just
            * call this method repetitively). Changes can be entered in
            * any order (they are automatically sorted). However, if
            * two changes appear to occur to the same date, the first
            * entered will also be the first to be applied. Note that
            * many parameters imply constraints on values: all dates
            * must be >= 0, population indices must refer to a valid
            * population, pairwise migration rate population indices
            * must be different, selfing rate (and other probabilities
            * must be >= 0 and <= 1, migration rates must be > 0. Beware
            * that final coalescence of all sample must be allowed
            * (beware of null migration rates and exponential decline
            * rates). As a rule, the change, Params and Coal classes DO
            * NOT perform sanity checks, so most invalid values will
            * result in program crash, incorrect results or exceptions.
            * The user should carefully check passed values or ensure
            * passed values are valid for all change type used. Note that
            * Params does not take ownership of the Event pointer and
            * will not delete the object at destruction time.
            *
            */
            void addChange(Event * e);

           /** \brief Clear the event list
            *
            */
            void clearChanges();

           /** \brief Get the number of changes
            *
            */
            unsigned int numChanges() const;

           /** \brief Get the first change
            *
            * Returns the first event of the queue (irrespective of the
            * internal event counter). The value is NULL if there is no
            * event loaded.
            *
            */
            Event * firstChange();

           /** \brief Get the date of the next change
            *
            * Returns the date of the next change to be applied. The
            * number will be >= 0 if there is at least one change, and
            * egglib::UNDEF otherwise.
            *
            */
            double nextChangeDate() const;

           /** \brief Apply next change
            *
            * The argument must be the address of the Coalesce instance
            * that should be modified (note that not all change types
            * will require affecting the Coalesce instance directly).
            *
            */
            void nextChangeDo(Coalesce * coal);

           /** \brief Number of delayed sample changes left to apply
            *
            */
            unsigned int nDSChanges() const;

            /** \brief Date of the last size change for a given pop
             *
             */
             double lastChange(unsigned int pop) const;

            /** \brief Restore object to the initial state
             *
             * All variable parameters (that is, those who have a
             * reversible method, and the number of populations) will
             * be reset to the last value set using a normal setter (or
             * to the default value).
             *
             */
            void restore();

           /** \brief Check consistency of entered parameters
            *
            * This methods ensures that the number of populations of the
            * migration matrix is consistent with the one of the current
            * Params instance, and that all changes have valid parameter
            * values (increase date, valid indices, etc.). The method
            * throws an EggArgumentValueError for the first encounterer
            * problem.
            *
            */
            void validate() const;

            /// \brief Gives the number of ingroup/outgroup samples (including delayed samples)
            void get_nsam(unsigned int& ns, unsigned int& no);

        protected:

            void _init(unsigned int npop, double migr);  // initialize members
            void alloc_npop(unsigned int size);
            void alloc_nalleles(unsigned int size);
            void alloc_nsites(unsigned int size);

            double _recomb;
            double _recomb_cache;
            Migration * _migr;

            unsigned long _max_iter;
            unsigned int _npop;
            unsigned int _npop_r;
            unsigned int * _n1;
            unsigned int * _n2;
            unsigned int _outgroup;
            double * _popsize;
            double * _popsize_cache;
            double * _growthrate;
            double * _growthrate_cache;
            double * _selfing;
            double * _selfing_cache;
            double * _lastChange;

            double _theta;
            MutationModel _mutmodel;
            unsigned int _nalleles;
            unsigned int _nalleles_r;
            bool _random_start_allele;
            unsigned int _fixed;
            double _TPMproba;
            double _TPMparam;

            bool _transW_matrix;
            double ** _transW;  // size _nalleles * _nalleles

            unsigned int _nsites;
            unsigned int _nsites_r;
            double * _sitePos;
            double * _siteW;
            double _totalSiteW;

            Event * _base_change;
            Event * _cur_change;
            unsigned int _num_changes;
            unsigned int _num_DSchanges;
            unsigned int _num_DSchanges_cache;

        private:

           /** \brief Copy constructor not available */
            Params(const Params& src) {}

           /** \brief Copy assignment operator not available */
            Params& operator=(const Params& src) { return *this; }
    };
}

#endif
