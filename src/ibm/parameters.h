#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <string>
struct Parameters {

    double d = {0.3}; // dispersal rates of species 1 and 2
    int max_time_steps{50000}; // number of generations the simulation is supposed to last
    unsigned int npatches{250}; // number of patches

    // number of individuals per patch
    // of each species
    int npp = 5;

    // initial value of belligerence and bravery
    double init_belligerence = 0.5;
    double init_bravery = 0.5;

    double omega = 0.5;

    // the name of the file to which data is written
    std::string base_name{"sim_ibm_warfare"};

    // interval in generations at which data is written to the data
    int data_interval = 10;

    // standard deviation in mutational effect size
    double sdmu = 0.02;

    // mutational probability baseline belligerence
    double mu_a_bel = 0.01;
    // mutational probability slope belligerence
    double mu_b_bel = 0.01;
    // mutational probability baseline bravery
    double mu_a_brav = 0.01;
    // mutational probability slope bravery
    double mu_b_brav = 0.01;
};

#endif
