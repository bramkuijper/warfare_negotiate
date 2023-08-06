#!/usr/bin/env python3

import datetime
import numpy as np

d = [0.1,0.5,0.9]

npp = [5]
npatches = [1000]

h = [0.05,0.1,0.25,0.5,0.75,0.9, 0.95]

omega = [ 0.5, 0.9] 

mu_b_bel = [0.0]
mu_b_brav = [0.0]

date = datetime.datetime.now()
base_name = "sim_stress_" +\
        f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"

exe_name = "WarfareIBM"  

ctr = 0
nrep=5
for rep_i in range(0,nrep):
    for d_i in d:
        for h_i in h:
            for npp_i in npp:
                for npatches_i in npatches:
                    for omega_i in omega:
                        for mu_b_bel_i in mu_b_bel:
                            for mu_b_brav_i in mu_b_brav:
                                base_name_i = base_name + str(ctr)
                                str_out = f"./{exe_name} {base_name_i} {d_i} {h_i} {npp_i} {npatches_i} {omega_i} {mu_b_bel_i} {mu_b_brav_i}"
                                ctr += 1
                                print(str_out)
