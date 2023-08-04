#include <stdexcept>
#include <random>
#include <cassert>
#include <vector>
#include <array>
#include "ibm_warfare.h"
#include "parameters.h"

// constructor expects a parameter object
IBM_Warfare::IBM_Warfare(Parameters const &params) : // constructors first initialize data members
    rd{} // initialize random device, see *.h file
    ,seed{rd()} // initialize seed
    ,rng_r{seed} // initialize the random number generator
    ,uniform{0.0,1.0} // initialize the uniform distribution
    ,patch_sampler{0,(int)params.npatches - 1} // initialize uniform distribution to sample patch indices from, as we are counting from 0 this can be any number from 0 to npatches - 1
    ,data_file{params.base_name.c_str()} // initialize the data file by giving it a name
    ,par{params} // initialize the parameter data member with the constructor argument
    ,metapop{par.npatches, Patch(par.npp, par)} // initialize a meta population each with n1 individuals of species 1 and n2 individuals of species 2
{
    // write the output
    write_data_headers();

    for (time_step = 0; 
            time_step <= par.max_time_steps; 
            ++time_step)
    {
        produce_offspring();
        replace_breeders_by_juvs();
        attack();

        if (time_step % par.data_interval == 0 ||
                time_step == par.max_time_steps)
        {
            write_data();
        }
    }

    write_parameters();
} // end IBM_Warfare::IBM_Warfare()

double IBM_Warfare::Cx(double const xphen)
{
    return(std::pow(xphen,4.0));
}

double IBM_Warfare::Cy(double const yphen)
{
    return(std::pow(yphen,4.0));
}

// see which patches attack other patches
double IBM_Warfare::attack_sealed_bid()
{
    for (int patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        metapop[patch_idx].patches_attack_to.clear();
        metapop[patch_idx].patches_attack_from.clear();
    }


    // Determine which patches are going to be attacked by each group.
    // The per-patch attack probability is then given by a/(npatch - 1)
    // hence the total number of patches being attacked is 
    // Binom(a/(npatch-1),npatch-1)
    //
    int nd_minus_1 = metapop.size() - 1.0;

    int npatches_to_attack, attacked_patch_idx;

    // first determine number of attacks of each patch
    for (int patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        // sample number of patches subject to attach
        std::binomial_distribution sample_patch_to_attack(
                nd_minus_1
                ,metapop[patch_idx].a_bel_phen_group / nd_minus_1);

        npatches_to_attack = sample_patch_to_attack(rng_r);

        for (int to_attack_patch_counter = 0; 
                to_attack_patch_counter < npatches_to_attack; 
                ++to_attack_patch_counter)
        {
            // sample index of patch
            do {
                attacked_patch_idx = patch_sampler(rng_r);
            } 
            while(attacked_patch_idx == patch_idx);

            // ok found victim patch, noted
            metapop[patch_idx].patches_attack_to.push_back(attacked_patch_idx);
            metapop[attacked_patch_idx].patches_attack_from.push_back(patch_idx);
        }
    } // end int patch_idx

    // now that all attacks have been noted let's pick one per patch 
    // for escalation
    for (int patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        n_patches_attacked_from = 
            metapop[patch_idx].patches_attacked_from.size();

        std::uniform_discrete_distribution patch_contest(0, 
                metapop[patch_idx].patches_attacked_from.size() - 1);


    }
} // double IBM_Warfare::attack_sealed_bid


// produce offspring
void IBM_Warfare::produce_offspring()
{
    for (int patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        metapop[patch_idx].juveniles.clear();
    }

    // some auxiliary variables
    double bel_phen, brav_phen, fecundity;

    std::vector <double> global_fecundity{0};
    std::vector <int> global_patch_id{0};
    std::vector <int> global_individual_id{0};
    double global_fecundity_total = 0.0;

    // now produce offspring dependent on an individual's level of x, y
    for (int patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        metapop[patch_idx].local_fecundity.clear();
        metapop[patch_idx].local_fecundity_total = 0.0;

        for (int breeder_idx = 0; 
                breeder_idx < metapop[patch_idx].breeders.size(); 
                ++breeder_idx)
        {
            bel_phen = metapop[patch_idx].breeders[breeder_idx].bel_phen;
            brav_phen = metapop[patch_idx].breeders[breeder_idx].brav_phen;

            fecundity = (1.0 - Cx(bel_phen)) * (1.0 - Cy(brav_phen));

            // add fecundity to global stack of fecundities
            // this is for later sampling of immigrant juveniles
            // (i.e., more fecund patches more likely to be chosen)
            global_fecundity.push_back(fecundity);
            global_patch_idx.push_back(patch_idx);
            global_individual_idx.push_back(breeder_idx);

            // add to global total of fecundities
            // this is for later sampling for immigrant juveniles
            global_fecundity_total += fecundity;

            // add fecundity to local stack of fecundities
            // this is for later sampling of native juveniles
            // (i.e., more fecund breeders more likely to be chosen)
            metapop[patch_idx].local_fecundity.push_back(fecundity);

            // add fecundity to local total of fecundities
            metapop[patch_idx].local_fecundity_total += fecundity;
        }
    } // all fecundities calculated now sample

    // make global fecundity sampler 
    // sampling from this will produce an index that can be used to produce
    // patch and breeder indices
    std::discrete_distribution <int> global_fecundity_sampler(
            global_fecundity.begin(), global_fecundity.end());

    double local_feq_conqueree;
    double local_feq_conqueror;

    std::array <double> fecundity_origin_list(3,0.0);

    int conquered_by_patch_idx;
        

    // generate a sampling distribution of any remote fecundities
    std::discrete_distribution <int> remote_fecundity_sampler(
            global_fecundity.begin()
            ,global_fecundity.end());

    // now produce offspring dependent on an individual's level of x, y
    for (int patch_idx = 0; 
            patch_idx < metapop.size(); 
            ++patch_idx)
    {
        conquered_by_patch_idx = metapop[patch_idx].conquered_by;

        // make a sampler for either from 
        // 1. local
        // 2. local, conqueror (if conquered)
        // 3. remote
        //
        fecundity_origin_list[0] = 
            (1.0 - par.d) * par.h * metapop[patch_idx].local_fecundity;

        fecundity_origin_list[1] = conquered_by_patch_idx >= 0 ?
            (1.0 - par.d) * (1.0 - par.h) * 
                metapop[conquered_by_patch_idx].local_fecundity_total
                :
                0.0;

        fecundity_origin_list[2] = par.d * global_fecundity_total;

        // generate a sampling distribution for sampling the origin of any 
        // breeder
        std::discrete_distribution <int> origin_sampler(fecundity_origin_list.begin()
                ,fecundity_origin_list.end());

        // aux variables to sample new breeder
        int origin_id, parent_idx, patch_parent_idx;
            
        std::discrete_distribution <int> local_fecundity_sampler(
                    metapop[patch_idx].local_fecundity.begin()
                    ,metapop[patch_idx].local_fecundity.end());

        std::discrete_distribution <int> conqueror_fecundity_sampler();

        // TODO
        if (conquered_by_patch_idx >= 0)
        {
            conqueror_fecundity_sampler
        }

        // clear existing juvs
        metapop[patch_idx].juveniles.clear();

        // sample new individuals to take up breeder positions
        for (int breeder_idx = 0; 
                breeder_idx < params.npp; 
                ++breeder_idx)
        {
            // sample origin id (0,1,2)
            origin_id = origin_sampler(rng_r);

            switch(origin_id)
            {
                case 0: // sample from local patch
                    {
                        patch_parent_idx = patch_idx;
                        parent_idx = local_fecundity_sampler(rng_r);
                        break;
                    }
                case 1: // sample from conqueror's patch
                    {
                        patch_parent_idx = conquered_by_patch_idx;
                        parent_idx = conqueror_fecundity_sampler(rng_r);
                        break;
                    }
                case 2: // sample from global patch
                    {
                        // store index from the sampler, this can be used to 
                        // lookup patch idx and breeder idx globally
                        int global_fecundity_idx =  global_fecundity_sampler(rng_r);
                        patch_parent_idx = global_patch_idx[global_fecundity_idx];
                        parent_idx = global_individual_idx[global_fecundity_idx];
                        break;
                    }
                case 3:
                    // TODO throw a range error
                    std::range_error("Wrong range when sampling origin id.");
                    break;
            }

            assert(patch_parent_idx >= 0);
            assert(patch_parent_idx < metapop.size());

            assert(parent_idx >= 0);
            assert(parent_idx < metapop[patch_parent_idx].breeders.size());

            // make new individual from parent
            Individual Kid(
                metapop[patch_parent_idx].breeders[parent_idx]
                ,rng_r
                ,par);

            metapop[patch_idx].juveniles.push_back(Kid);
        } // end for breeder idx
    } // end patch_idx
    
} // end IBM_Warfare::reproduce()

// TODO: if one has a larger patch make sure replacement also accommodates for this
void IBM_Warfare::replace_breeders_by_juvs()
{
    // now produce offspring dependent on an individual's level of x, y
    for (int patch_idx = 0; 
            patch_idx < metapop.size(); 
            ++patch_idx)
    {
        assert(metapop[patch_idx].juveniles.size() == par.npp);

        for (int breeder_idx = 0; 
                breeder_idx < params.npp; 
                ++breeder_idx)
        {
            metapop[patch_idx].breeders[breeder_idx] = 
                metapop[patch_idx].juveniles[breeder_idx];
        }
    } // end patch idx
}


//  write parameters to the output file
void IBM_Warfare::write_parameters()
{
    data_file << std::endl
        << std::endl;

    // write parameters to the file for each species
    for (int species_idx = 0; species_idx < 2; ++species_idx)
    {
        data_file << "d" << (species_idx + 1) << ";"
                    << par.initial_d[species_idx] << std::endl

                << "npp" << (species_idx + 1) << ";"
                    << par.npp[species_idx] << std::endl

                << "baseline_survival" << (species_idx + 1) << ";"
                    << par.baseline_survival[species_idx] << std::endl

                << "baseline_fecundity" << (species_idx + 1) << ";"
                    << par.baseline_fecundity[species_idx] << std::endl

                << "strength_survival" << (species_idx + 1) << ";"
                    << par.strength_survival[species_idx] << std::endl

                << "fecundity_help" << (species_idx +1 ) << ";"
                    << par.initial_fec_h[species_idx] << std::endl

                << "survival_help" << (species_idx + 1) << ";"
                    << par.initial_surv_h[species_idx] << std::endl

                << "survival_cost_of_surv_help" << (species_idx + 1) << ";"
                    << par.survival_cost_of_surv_help[species_idx] << std::endl

                << "survival_cost_of_fec_help" << (species_idx + 1) << ";"
                    << par.survival_cost_of_fec_help[species_idx] << std::endl

                << "fecundity_cost_of_surv_help" << (species_idx + 1) << ";"
                    << par.fecundity_cost_of_surv_help[species_idx] << std::endl

                << "fecundity_cost_of_fec_help" << (species_idx + 1) << ";"
                    << par.fecundity_cost_of_fec_help[species_idx] << std::endl

                << "juvenile_survival_weight" << (species_idx + 1) << ";"
                    << juvenile_survival_weight[species_idx] << std::endl;
    }

    data_file << "mu_fec_h;" << par.mu_fec_h << std::endl
                << "mu_surv_h;" << par.mu_surv_h << std::endl
                << "mu_d;" << par.mu_disp << std::endl
                << "sdmu;" << par.sdmu << std::endl
                << "between_species;" << par.between_species << std::endl
                << "death_birth;" << par.death_birth << std::endl
                << "partner_mechanism;" << par.partner_mechanism << std::endl
                << "fidelity_prob;" << par.fidelity_prob << std::endl
                << "negotiate_once;" << par.negotiate_once << std::endl
                << "sd_pcerr;" << par.sd_pcerr << std::endl
                << "seed;" << seed << std::endl;

} // end write_parameters()

void IBM_Warfare::write_data_headers()
{
    data_file << "time_step;";

    for (int species_idx = 1; species_idx <= 2; ++species_idx)
    {
        data_file
            << "mean_disp" << species_idx << ";"
            << "mean_fec_h" << species_idx << ";"
            << "mean_surv_h" << species_idx << ";"
            << "mean_given_fec_h" << species_idx << ";"
            << "mean_given_surv_h" << species_idx << ";"

            << "var_disp" << species_idx << ";"
            << "var_fec_h" << species_idx << ";"
            << "var_surv_h" << species_idx << ";"
            << "var_given_fec_h" << species_idx << ";"
            << "var_given_surv_h" << species_idx << ";"

            << "mean_surv_prob" << species_idx << ";"
            << "mean_surv_help_per_ind" << species_idx << ";"
            << "patch_occupancy" << species_idx << ";"
            << "nsurvivors" << species_idx << ";"
            << "mean_offspring" << species_idx << ";"
            
            << "mean_adult_survival_weight" << species_idx << ";"
            << "mean_juvenile_survival_weight" << species_idx << ";";
    }

    data_file << std::endl;
} // end IBM_Warfare::write_data_headers()

// write the data to an output file
void IBM_Warfare::write_data()
{
    // store means and sums of squares
    // (latter for variance calculations)
    double mean_disp[2]         = {0.0,0.0};
    double mean_fec_h[2]        = {0.0,0.0};
    double mean_surv_h[2]       = {0.0,0.0};
    double mean_given_fec_h[2]  = {0.0,0.0};
    double mean_given_surv_h[2] = {0.0,0.0};

    double ss_disp[2]           = {0.0,0.0};
    double ss_fec_h[2]          = {0.0,0.0};
    double ss_surv_h[2]         = {0.0,0.0};
    double ss_given_fec_h[2]    = {0.0,0.0};
    double ss_given_surv_h[2]   = {0.0,0.0};

    // aux variables to store trait values
    double disp,f,s,gf,gs;

    // array for counts of population sizes
    // (although for now this should be simply metapop.size() * npp[species_idx])
    int n_events[2] = {0,0};

    // go through all patches
    for (int patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        // go through the 2 species
        for (int species_idx = 0; species_idx < 2; ++species_idx)
        {
            for (ind_iter individual_iter =
                    metapop[patch_idx].breeders[species_idx].begin();
                    individual_iter != metapop[patch_idx].breeders[species_idx].end();
                    ++individual_iter)
            {
                // obtain allelic values
                disp    = (individual_iter->d[0]        + individual_iter->d[1]) / 2;
                f       = individual_iter->fec_h[0]     + individual_iter->fec_h[1];
                s       = individual_iter->surv_h[0]    + individual_iter->surv_h[1];
                gf      = individual_iter->given_fec_h;
                gs      = individual_iter->given_surv_h;

                mean_disp[species_idx]          += disp;
		        ss_disp[species_idx]            += disp * disp;

                mean_fec_h[species_idx]         += f;
                ss_fec_h[species_idx]           += f * f;

                mean_surv_h[species_idx]        += s;
                ss_surv_h[species_idx]          += s * s;

                mean_given_fec_h[species_idx]   += gf;
                ss_given_fec_h[species_idx]     += gf *gf;

                mean_given_surv_h[species_idx]  += gs;
                ss_given_surv_h[species_idx]    += gs * gs;

                // update population stats
                // (although for now, it should be the same as
                ++n_events[species_idx];
            }
        } // end for species_idx
    } // end for patch_idx


    // output the time step for starters
    data_file << time_step << ";";

    // calculate means variances of traits for each species
    for (int species_idx = 0; species_idx < 2; ++species_idx)
    {
        mean_disp[species_idx]          /= n_events[species_idx];
        mean_fec_h[species_idx]         /= n_events[species_idx];
        mean_surv_h[species_idx]        /= n_events[species_idx];
        mean_given_fec_h[species_idx]   /= n_events[species_idx];
        mean_given_surv_h[species_idx]  /= n_events[species_idx];

        // var = E[x^2] - E[x]^2
        double var_disp = ss_disp[species_idx] / n_events[species_idx] -
	        mean_disp[species_idx] * mean_disp[species_idx];

        double var_fec_h = ss_fec_h[species_idx] / n_events[species_idx] -
            mean_fec_h[species_idx] * mean_fec_h[species_idx];

        double var_surv_h = ss_surv_h[species_idx] / n_events[species_idx] -
            mean_surv_h[species_idx] * mean_surv_h[species_idx];

        double var_given_fec_h = ss_given_fec_h[species_idx] / n_events[species_idx] -
            mean_given_fec_h[species_idx] * mean_given_fec_h[species_idx];

        double var_given_surv_h = ss_given_surv_h[species_idx] / n_events[species_idx] -
            mean_given_surv_h[species_idx] * mean_given_surv_h[species_idx];

        data_file
                    << mean_disp[species_idx] << ";"
                    << mean_fec_h[species_idx] << ";"
                    << mean_surv_h[species_idx] << ";"
                    << mean_given_fec_h[species_idx] << ";"
                    << mean_given_surv_h[species_idx] << ";"

                    << var_disp << ";"
                    << var_fec_h << ";"
                    << var_surv_h << ";"
                    << var_given_fec_h << ";"
                    << var_given_surv_h << ";"

                    << mean_surv_prob[species_idx] << ";"
                    // TODO: missing mean_fec_help_per_individual?
                    << mean_surv_help_per_individual[species_idx] << ";"
                    << patch_occupancy[species_idx] << ";"
                    << nsurvivors[species_idx] << ";"
                    << mean_offspring[species_idx] << ";"

                    << mean_adult_survival_weight[species_idx] << ";"
                    << mean_juvenile_survival_weight[species_idx] << ";";
    } // end for species_idx

    data_file << std::endl;

} // end IBM_Warfare::write_data()
