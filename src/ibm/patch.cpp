#include "patch.h"
#include "individual.h"
#include <random>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

Patch::Patch(
        int const nbreeders
        ,Parameters const &params
        ) :
    par{params}
    ,breeders(nbreeders, Individual(
            params.init_belligerence
            ,params.init_bravery))
    ,juveniles{0}
    ,local_fecundity{0}
    ,local_fecundity_total{0}
    ,patches_attack_to{0}
    ,patches_attack_from{0}
{
} // end constructor

// copy constructor
Patch::Patch(Patch const &other) :
    par{other.par}
    ,breeders{other.breeders}
    ,juveniles{other.juveniles}
    ,local_fecundity{other.local_fecundity}
    ,local_fecundity_total{other.local_fecundity_total}
    ,patches_attack_to{other.patches_attack_to}
    ,patches_attack_from{other.patches_attack_from}
{
}

// assignment operator
void Patch::operator=(Patch const &other)
{
    par = other.par;
    breeders = other.breeders;
    juveniles = other.juveniles;
    local_fecundity = other.local_fecundity;
    local_fecundity_total = other.local_fecundity_total;
    patches_attack_to = other.patches_attack_to;
    patches_attack_from = other.patches_attack_from;
}
