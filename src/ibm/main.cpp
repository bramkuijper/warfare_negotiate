#include <cstdlib>
#include <string>
#include "ibm_mutualism.h"

int main(int argc, char **argv)
{
    Parameters parms{};

    parms.base_name = argv[1];

    parms.npp[0] = atof(argv[2]);
    parms.npp[1] = atof(argv[2]);
    
    parms.mu_fec_h = atof(argv[3]);
    parms.mu_surv_h = atof(argv[4]);

    parms.initial_d[0] = atof(argv[5]);
    parms.initial_d[1] = atof(argv[5]);

    // parms.between_species = argv[2] == std::string{"true"};
    // parms.death_birth = argv[3] == std::string{"true"};
    // parms.partner_mechanism = atof(argv[4]);  

    // parms.initial_fec_h[0] = atof(argv[4]);
    // parms.initial_fec_h[1] = atof(argv[5]);

    // parms.initial_surv_h[0] = atof(argv[6]);
    // parms.initial_surv_h[1] = atof(argv[7]);

    IBM_Mutualism sim{parms};

    return 0;
}
