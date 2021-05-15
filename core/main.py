import numpy as np      #TODO import only numpy functions used
from pandas import DataFrame

from conf.settings import *
from core.mcmc import *
from core.ow import *
from core.rpm import *
from core.stoch_sim import *
from core.plotting import *


top = SIMUL_PARAS['top']
base = SIMUL_PARAS['base']
sr = SIMUL_PARAS['sr']
nfacies = SIMUL_PARAS['nfacies']

P = np.matrix([[SIMUL_PARAS['p11'], SIMUL_PARAS['p12'], SIMUL_PARAS['p13']],
[SIMUL_PARAS['p21'], SIMUL_PARAS['p22'], SIMUL_PARAS['p23']],
[SIMUL_PARAS['p31'], SIMUL_PARAS['p32'], SIMUL_PARAS['p33']]])

code = SIMUL_PARAS['code']
color = SIMUL_PARAS['color']

mux = FACIES_PARAS['mux']
muy = FACIES_PARAS['muy']
covf1 = FACIES_PARAS['covf1']
covf2 = FACIES_PARAS['covf2']
covf3 = FACIES_PARAS['covf3']
cov = [covf1, covf2, covf3]

lh1 = STOCH_PARAS['lh1']
lh2 = STOCH_PARAS['lh2']
lh3 = STOCH_PARAS['lh3']

ow = FINAL_PARAS['ow']
so = FINAL_PARAS['so']
rho_oil = FINAL_PARAS['rho_oil']
k_oil = FINAL_PARAS['k_oil']
rho_water = FINAL_PARAS['rho_water']
k_water = FINAL_PARAS['k_water']
rho_min = FINAL_PARAS['rho_min']
rpm = FINAL_PARAS['rpm']
rpm_p = FINAL_PARAS['rpm_p']
rpm_k = FINAL_PARAS['rpm_k']
rpm_g = FINAL_PARAS['rpm_g']
phi_c = FINAL_PARAS['phi_c']
coord_n = FINAL_PARAS['coord_n']
c_error = FINAL_PARAS['c_error']
r_error = FINAL_PARAS['r_error']

#TODO configure the script for the user feel free to type n facies


def run():

    depth = np.arange(top, base, sr)

    facies, code_facies = mcmc(nfacies, depth, P, code, sf=1)
    transition_matrix_show(P)

    f1 = (facies == 1)
    f2 = (facies == 2)
    f3 = (facies == 3)

    rv = []
    rv_pdf = []

    for i in range(3):

        rv.append(gaussian_gen(mux[i], muy[i], cov[i])[0])
        rv_pdf.append(gaussian_gen(mux[i], muy[i], cov[i])[1])
    
    f1_p = sampling_from_gaussian(rv[0], depth, f1)
    f2_p = sampling_from_gaussian(rv[1], depth, f2)
    f3_p = sampling_from_gaussian(rv[2], depth, f3)


    m1 = np.empty(len(depth))
    m2 = np.empty(len(depth))

    m1[f1], m2[f1] = f1_p[0], f1_p[1]
    m1[f2], m2[f2] = f2_p[0], f2_p[1]
    m1[f3], m2[f3] = f3_p[0], f3_p[1]

    model_f1 = syn_variogram(depth, f1_p[2][0, 0], lh1)
    model_f2 = syn_variogram(depth, f2_p[2][0, 0], lh2)
    model_f3 = syn_variogram(depth, f3_p[2][0, 0], lh3)

    sim_x_f1, sim_y_f1 = simulations(model_f1, f1_p[2], m1, m2)
    sim_x_f2, sim_y_f2 = simulations(model_f2, f2_p[2], m1, m2)
    sim_x_f3, sim_y_f3 = simulations(model_f3, f3_p[2], m1, m2)

    sim_x = np.empty(len(depth))
    sim_y = np.empty(len(depth))

    sim_x[f1], sim_y[f1] = sim_x_f1[f1], sim_y_f1[f1]
    sim_x[f2], sim_y[f2] = sim_x_f2[f2], sim_y_f2[f2]
    sim_x[f3], sim_y[f3] = sim_x_f3[f3], sim_y_f3[f3]

    phi = sim_x
    vclay = sim_y

    sw = sw_gen(depth, ow, so)
    #TODO identify if one of lithologies is Shale and define SW = 1.0 in these intervals

    rhob = rhob_gen(depth, phi, ow, rho_min, rho_oil, rho_water)

    Ksoft, Gsoft = soft_sand(rpm_k, rpm_g, rho_min, phi, phi_c, coord_n, rpm_p)
    Ksoft[f1], Gsoft[f1] = soft_sand(K=21.0 * 10**9, G=7.0 * 10**9, rho=2.58, phi=phi[f1], phic=phi_c, n=11.0, P=rpm_p)   # Using other rpm parameters for some lithologies.

    Gsat = Gsoft
    Ksat = np.empty(len(Gsat))
    Ksat[(depth <= ow)] = gassmann(Ksoft[(depth <= ow)], rpm_k, k_oil, phi[(depth <= ow)])
    Ksat[(depth > ow)] = gassmann(Ksoft[(depth > ow)], rpm_k, k_water, phi[(depth > ow)])

    Msat = Ksat + (Gsat * 4./3.)
    Vp = (Msat/10**9/rhob)**0.5
    Vs = (Gsat/10**9/rhob)**0.5

    #TODO introduce spatial correlated error

    logplots(depth, vclay, phi, sw, Vp, Vs, rhob, f1, f2, f3)
    crossplot_vpvs(Vp, Vs, z=vclay, zlabel='VCLAY [dec]')

    DataFrame(data=np.array([depth, facies, code_facies, vclay, phi, sw, Vp, Vs, rhob]).T,
              columns=['DEPTH', 'CODE', 'FACIES', 'VCLAY', 'PHI', 'SW', 'VP', 'VS', 'RHOB']).to_csv('results/pseudowell.las',
                                                                                                     sep = ',', index=False)
