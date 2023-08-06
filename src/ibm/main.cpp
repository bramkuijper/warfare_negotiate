#include <cstdlib>
#include <string>
#include "ibm_warfare.h"

int main(int argc, char **argv)
{
    Parameters parms{};

    parms.base_name = argv[1];

    parms.d = std::stod(argv[2]);
    parms.h = std::stod(argv[3]);
    parms.npp = std::stoi(argv[4]);
    parms.npatches = std::stoi(argv[5]);
    parms.omega = std::stod(argv[6]);
    parms.mu_b_bel = std::stod(argv[7]);
    parms.mu_b_brav = std::stod(argv[8]);

    IBM_Warfare sim{parms};

    return 0;
}
