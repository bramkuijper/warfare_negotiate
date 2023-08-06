#ifndef IBM_WARFARE_H
#define IBM_WARFARE_H

#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "patch.h"
#include "parameters.h"

// struct with default parameters for the class constructor

class IBM_Warfare
{
    private:
        // a data file containing the results of the simulation
        std::ofstream data_file;

        // keep track of the time step of the simulation
        long unsigned time_step = 0;

        // random device which is used to generate
        // proper random seeds
        std::random_device rd;

        // store the random seed
        // we need to store this so that we can output the
        // random seed, so that we could 'replay' the exact
        // same sequence of random numbers for debugging purposes etc
        unsigned int seed;

        // variable for stats purposes: 
        // total fecundity in population
        double global_fecundity_total;

        // random number generator
        std::mt19937 rng_r;

        // uniform distribution to compare against probabilities
        std::uniform_real_distribution<double> uniform;

        // uniform distribution to get random patch
        std::uniform_int_distribution<int> patch_sampler;

        // parameter object
        // containing all the parameters for this run
        Parameters par;

        // save number of fights
        int n_escalated_fights = 0;

    public:
        // metapopulation of patches
        std::vector<Patch> metapop;

        // the class constructor
        IBM_Warfare(const Parameters &params);

        void reproduce();
        double Cy(double const yphen);
        double Cx(double const xphen);
        void attack();
        void attack_sealed_bid();

        void escalated_contest(int const attacking_patch_idx
                ,int const defending_patch_idx);

        void produce_offspring();
        void replace_breeders_by_juvs();

        void write_parameters();
        void write_data();
        void write_data_headers();

}; // end class IBM_Warfare

#endif
