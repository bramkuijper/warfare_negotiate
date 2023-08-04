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
    // clear list of attacked patches
    for (int patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        metapop[patch_idx].patches_attack_to.clear();
        metapop[patch_idx].patches_attack_from.clear();
        metapop[patch_idx].conquered_by = -1;
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
        std::uniform_discrete_distribution patch_contest_sampler(0, 
                metapop[patch_idx].patches_attacked_from.size() - 1);

        int contest_idx = patch_contest_sampler(rng_r);

        int contest_patch_idx = metapop[patch_idx].patches_attacked_from[contest_idx];

        // range checking
        assert(contest_patch_idx >= 0);
        assert(contest_patch_idx < params.npatches);

        escalated_contest(contest_patch_idx, patch_idx);
    }
} // double IBM_Warfare::attack_sealed_bid

void IBM_Warfare::escalated_contest(
        int const attacking_patch_idx
        ,int const defending_patch_idx)
{
    double g_attacker = exp(
            metapop[attacking_patch_idx].breeders.size() *
            metapop[attacking_patch_idx].a_brav_phen_group);

    double g_defender = exp(
            metapop[defending_patch_idx].breeders.size() *
            metapop[defending_patch_idx].a_brav_phen_group);

    double prob_winning = par.omega * g_attacker /
        (par.omega * g_attacker + (1.0 - par.omega) * g_defender);

    if (uniform(rng_r) < prob_winning)
    {
        metapop[defending_patch_idx].conquered_by = attacking_patch_idx;
    }
} // end escalated_contest






// produce offspring
void IBM_Warfare::produce_offspring()
{
    // reset juveniles
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
        // collect patch id of patch that has conquered focal
        // if no conquest took place, this is negative
        conquered_by_patch_idx = metapop[patch_idx].conquered_by;

        // make a sampler for where to pick juveniles from
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
                0.0; // if not conquered by others this is 0

        fecundity_origin_list[2] = par.d * global_fecundity_total;

        // generate a sampling distribution for sampling the origin of any 
        //  breeder delivering a juvenile
        std::discrete_distribution <int> origin_event_sampler(
                fecundity_origin_list.begin()
                ,fecundity_origin_list.end());

        // aux variables to sample new breeder
        int origin_id, parent_idx, patch_parent_idx;
           
        // sampler for the local fecundity distribution
        std::discrete_distribution <int> local_fecundity_sampler(
                    metapop[patch_idx].local_fecundity.begin()
                    ,metapop[patch_idx].local_fecundity.end());

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
                        assert(conquered_by_patch_idx >= 0);
                        assert(conquered_by_patch_idx < params.npatches);

                        // make fecundity distribution of globals
                        // there is no way to initialize this distribution
                        // dynamically, so that we have to do this here... meh
                        std::discrete_distribution <int> conqueror_fecundity_sampler(
                                metapop[conquered_by_patch_idx].local_fecundity.begin()
                                metapop[conquered_by_patch_idx].local_fecundity.end()
                                );

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
                    {
                    // TODO throw a range error
                    std::range_error("Wrong range when sampling origin id.");
                    break;
                    }
            } // end switch

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

    data_file << "d;" << par.d << std::endl <<
        "omega;" << par.omega << std::endl <<
        "npp;" << par.npp << std::endl <<
        "sdmu;" << par.sdmu << std::endl <<
        "mu_a_bel;" << par.mu_a_bel << std::endl <<
        "mu_b_bel;" << par.mu_b_bel << std::endl <<
        "mu_a_brav;" << par.mu_a_brav << std::endl <<
        "mu_b_brav;" << par.mu_b_brav << std::endl <<
        "init_belligerence;" << par.init_belligerence << std::endl <<
        "init_bravery;" << par.init_bravery << std::endl; 
} // end write_parameters()

void IBM_Warfare::write_data_headers()
{
    data_file << "time_step;"
        << "mean_a_bel;" 
        << "mean_b_bel;"
        << "mean_belligerence;"
        << "mean_a_brav;"
        << "mean_b_brav;"
        << "mean_bravery;"
        << "var_a_bel;" 
        << "var_b_bel;"
        << "var_belligerence;"
        << "var_a_brav;"
        << "var_b_brav;"
        << "var_bravery;"
        << std::endl;
} // end IBM_Warfare::write_data_headers()

// write the data to an output file
void IBM_Warfare::write_data()
{
    double mean_a_bel{0.0};
    double ss_a_bel{0.0};
    double mean_b_bel{0.0};
    double ss_b_bel{0.0};
    double mean_bel{0.0};
    double ss_bel{0.0};

    double mean_a_brav{0.0};
    double ss_a_brav{0.0};
    double mean_b_brav{0.0};
    double ss_b_brav{0.0};
    double mean_brav{0.0};
    double ss_brav{0.0};

    int a_bel, b_bel, bel, a_brav, b_brav, brav;

    int nbreeder = 0;

    for (int patch_idx = 0; patch_idx < params.npatches; ++patch_idx)
    {
        nbreeder += metapop[patch_idx].breeders.size();
        for (int breeder_idx = 0; 
                breeder_idx < metapop[patch_idx].breeders.size(); ++breeder_idx)
        {
            a_bel = metapop[patch_idx].breeders[breeder_idx].a_bel;
            b_bel = metapop[patch_idx].breeders[breeder_idx].b_bel;
            bel = metapop[patch_idx].breeders[breeder_idx].bel;
            a_brav = metapop[patch_idx].breeders[breeder_idx].a_brav;
            b_brav = metapop[patch_idx].breeders[breeder_idx].b_brav;
            brav = metapop[patch_idx].breeders[breeder_idx].brav;

            mean_a_bel += a_bel;
            ss_a_bel += a_bel * a_bel;
            
            mean_b_bel += b_bel;
            ss_b_bel += b_bel * b_bel;
            
            mean_bel += bel;
            ss_bel += bel * bel;
            
            mean_a_brav += a_brav;
            ss_a_brav += a_brav * a_brav;
            
            mean_b_brav += b_brav;
            ss_b_brav += b_brav * b_brav;
            
            mean_brav += brav;
            ss_brav += brav * brav;
        }
    }

    mean_a_bel /= nbreeder;
    double var_a_bel = ss_a_bel / nbreeder - mean_a_bel * mean_a_bel;
    
    mean_b_bel /= nbreeder;
    double var_b_bel = ss_b_bel / nbreeder - mean_b_bel * mean_b_bel;
    
    mean_bel /= nbreeder;
    double var_bel = ss_bel / nbreeder - mean_bel * mean_bel;
    
    mean_a_brav /= nbreeder;
    double var_a_brav = ss_a_brav / nbreeder - mean_a_brav * mean_a_brav;
    
    mean_b_brav /= nbreeder;
    double var_b_brav = ss_b_brav / nbreeder - mean_b_brav * mean_b_brav;
    
    mean_brav /= nbreeder;
    double var_brav = ss_brav / nbreeder - mean_brav * mean_brav;

    data_file << time_step << ";"
        << mean_a_bel << ";"
        << mean_b_bel << ";"
        << mean_bel << ";"
        << mean_a_brav << ";"
        << mean_b_brav << ";"
        << mean_brav << ";"
        << var_a_bel << ";"
        << var_b_bel << ";"
        << var_bel << ";"
        << var_a_brav << ";"
        << var_b_brav << ";"
        << var_brav << ";"
        << std::endl;
} // end IBM_Warfare::write_data()
