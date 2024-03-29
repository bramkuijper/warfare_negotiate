#ifndef _PATCH_H_
#define _PATCH_H_

#include <random>
#include "individual.h"

class Patch
{
    public:
        // two dimensional list of individuals in each patch
        // which reflect the breeders of species 1 and species 2
        std::vector < Individual > breeders;
        std::vector < Individual > juveniles;

        std::vector < int > patches_attack_to;
        std::vector < int > patches_attack_from;

        std::vector <double> local_fecundity;
        double local_fecundity_total{0.0};

        double a_brav_phen_group{0.0};
        double a_bel_phen_group{0.0};

        int conquered_by{0};

        // the default constructor 
        Patch(int const n, Parameters const &params);

        // copy constructor
        Patch(Patch const &other);

        Parameters par; 
        
        void operator=(Patch const &other);

        // typedef
        typedef std::vector<Individual>::iterator ind_iter;
}; // end patch class declaration
#endif
