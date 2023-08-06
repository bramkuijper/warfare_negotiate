#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#include <random>
#include "parameters.h"

class Individual
{
    public:
        // belligerence reaction norm
        double a_bel{0.0};
        double b_bel{0.0};

        double bel_phen{0.0};

        // bravery reaction norm
        double a_brav{0.0};
        double b_brav{0.0};

        double brav_phen{0.0};

        Individual();

        // init constructor
        Individual(double belligerence, double bravery);

        // copy ctor 
        Individual(Individual const &other);

        // birth constructor
        Individual(Individual const &mother
                ,std::mt19937 &rng
                ,Parameters const &params);

        void operator=(Individual const &other);
};

#endif
