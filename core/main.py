import numpy as np

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

def run():

    depth = np.arange(top, base, sr)

    facies, code_facies = mcmc(nfacies, depth, P, code, sf=1)

    shale = (facies == 1)
    silty = (facies == 2)
    sand = (facies == 3)

    rv = []
    rv_pdf = []

    for i in range(3):

        rv.append(gaussian_gen(mux[i], muy[i], cov[i])[0])
        rv_pdf.append(gaussian_gen(mux[i], muy[i], cov[i])[1])
    
    shale_p = sampling_from_gaussian(rv[0], depth, shale)
    silty_p = sampling_from_gaussian(rv[1], depth, silty)
    sand_p = sampling_from_gaussian(rv[2], depth, sand)


    m1 = np.empty(len(depth))
    m2 = np.empty(len(depth))

    m1[shale], m2[shale] = shale_p[0], shale_p[1]
    m1[silty], m2[silty] = silty_p[0], silty_p[1]
    m1[sand], m2[sand] = sand_p[0], sand_p[1]

    model_shale = syn_variogram(depth, shale_p[2][0, 0], lh1)
    model_silty = syn_variogram(depth, silty_p[2][0, 0], lh2)
    model_sand = syn_variogram(depth, sand_p[2][0, 0], lh3)

    sim_x_shale, sim_y_shale = simulations(model_shale, shale_p[2], m1, m2)
    sim_x_silty, sim_y_silty = simulations(model_silty, silty_p[2], m1, m2)
    sim_x_sand, sim_y_sand = simulations(model_silty, sand_p[2], m1, m2)

    sim_x = np.empty(len(depth))
    sim_y = np.empty(len(depth))

    sim_x[shale], sim_y[shale] = sim_x_shale[shale], sim_y_shale[shale]
    sim_x[silty], sim_y[silty] = sim_x_silty[silty], sim_y_silty[silty]
    sim_x[sand], sim_y[sand] = sim_x_sand[sand], sim_y_sand[sand]

    phi = sim_x
    vclay = sim_y

    sw = sw_gen(depth, ow, so)

    rhob = rhob_gen(depth, phi, ow, rho_min, rho_oil, rho_water)

    Ksoft, Gsoft = soft_sand(rpm_k, rpm_g, rho_min, phi, phi_c, coord_n, rpm_p)

    Gsat = Gsoft
    Ksat = gassmann(Ksoft, rpm_k, k_oil, phi)

    Msat = Ksat + (Gsat * 4./3.)
    Vp = (Msat/10**9/rhob)**0.5
    Vs = (Gsat/10**9/rhob)**0.5

    plot = logplots(depth, vclay, phi, sw, Vp, Vs, rhob, shale, silty, sand)