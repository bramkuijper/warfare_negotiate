#include <random>
#include <algorithm>
#include <set>
#include "individual.h"
#include "parameters.h"

Individual::Individual(double belligerence, double bravery) :
    a_bel{belligerence},
    b_bel{0.0},
    bel_phen{belligerence},
    a_brav{bravery},
    b_brav{0.0},
    brav_phen{bravery}
{}

// copy constructor
Individual::Individual(Individual const &other) :
    a_bel{other.a_bel}
    ,b_bel{other.b_bel}
    ,bel_phen{other.bel_phen}
    ,a_brav{other.a_brav}
    ,b_brav{other.b_brav}
    ,brav_phen{other.brav_phen}
{}

// birth constructor
Individual::Individual(Individual const &parent
                ,std::mt19937 &rng
                ,Parameters const &params) :
    a_bel{parent.a_bel}
    ,b_bel{parent.b_bel}
    ,bel_phen{parent.bel_phen}
    ,a_brav{parent.a_brav}
    ,b_brav{parent.b_brav}
    ,brav_phen{parent.brav_phen}
{
    // set up random number distributions to realize
    // mutation rates...
    std::uniform_real_distribution<> uniform{0.0,1.0};

    // ... and mutational effect sizes
    std::normal_distribution<> normal{0,params.sdmu};

    if (uniform(rng) < params.mu_a_bel)
    {
        a_bel = a_bel + normal(rng);
    }

    if (uniform(rng) < params.mu_b_bel)
    {
        b_bel = b_bel + normal(rng);
    }

    // for now phenotype similar to a_bel:
    bel_phen = a_bel;
    
    if (uniform(rng) < params.mu_a_brav)
    {
        a_brav = a_brav + normal(rng);
    }

    if (uniform(rng) < params.mu_b_brav)
    {
        b_brav = b_brav + normal(rng);
    }
    
    // for now phenotype similar to a_brav:
    brav_phen = a_brav;
} // end parent constructor

void Individual::operator=(Individual const &other)
{
    a_bel = other.a_bel;
    b_bel = other.b_bel;
    bel_phen = other.bel_phen;
    a_brav = other.a_brav;
    b_brav = other.b_brav;
    brav_phen = other.brav_phen;
}
