import numpy as np      #TODO import only numpy functions used
import pandas as pd
pd.options.mode.chained_assignment = None

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
[SIMUL_PARAS['p31'], SIMUL_PARAS['p32'], SIMUL_PARAS['p33']]]).T

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
lh = [lh1, lh2, lh3]

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
boolean = FINAL_PARAS['use_dif_params']


def run():

    depth = np.arange(top, base, sr)
    df = pd.DataFrame(data=depth.T, columns=['Depth'])

    for i in range(nfacies):
        if (len(P) == nfacies) & (len(code) == nfacies):
            df['facies'], df['code_facies'] = mcmc(nfacies, depth, P, code, sf=0)
        else:
            raise ValueError('The number of facies is different from the codes of facies or matrix P dimensions. Check your parameters')


    transition_matrix_show(P)
    
    for i in range(nfacies):
        df[f'fbool{i}'] = (df['facies'] == i)

    rv = []
    rv_pdf = []

    for i in range(nfacies):
        rv.append(gaussian_gen(mux[i], muy[i], cov[i])[0])
        rv_pdf.append(gaussian_gen(mux[i], muy[i], cov[i])[1])

    m1 = np.empty(len(depth))
    m2 = np.empty(len(depth))
    covariance = [0, 0, 0]


    for i in range(nfacies):
        m1[df[f'fbool{i}']], m2[df[f'fbool{i}']], covariance[i] = sampling_from_gaussian(rv[i], depth, df[f'fbool{i}'])

    df['m1'] = m1
    df['m2'] = m2


    for i in range(nfacies):
        df[f'model_f{i}'] = syn_variogram(depth, covariance[i][0, 0], lh[i])


    for i in range(nfacies):
        df[f'sim_x_f{i}'], df[f'sim_y_f{i}'], c_ex = simulations(df[f'model_f{i}'], covariance[i], df['m1'], df['m2']) 


    sim_x = np.empty(len(depth))
    sim_y = np.empty(len(depth))


    for i in range(nfacies):
        sim_x[df[f'fbool{i}']], sim_y[df[f'fbool{i}']] = df[f'sim_x_f{i}'][df[f'fbool{i}']], df[f'sim_y_f{i}'][df[f'fbool{i}']]

    df['phi'] = sim_x
    df['vclay'] = sim_y


    df['sw'] = sw_gen(depth, ow, so)
    i = 0
    for words in code:
        if words.lower() == 'shale':
            index = i
            df['sw'][df[f'fbool{index}']] = 1.0
            break
        else:
            i = i+1

    df['rhob'] = rhob_gen(depth, df['phi'], ow, rho_min, rho_oil, rho_water)
    
    df['Ksoft'], df['Gsoft'] = soft_sand(rpm_k, rpm_g, rho_min, df['phi'], phi_c, coord_n, rpm_p)
    print(df.head(50))
    
    if boolean == True:
        df['Ksoft'][df[f'fbool{index}']], df['Gsoft'][df[f'fbool{index}']] = soft_sand(K=21.0 * 10**9, G=7.0 * 10**9, rho=2.58,
                                                                               phi=df['phi'][df[f'fbool{index}']], phic=phi_c,
                                                                               n=11.0, P=rpm_p)    
    
    rpm_plot(np.linspace(0, 1, 100),
             soft_sand(rpm_k, rpm_g, rho_min, np.linspace(0, 1, 100), phi_c, coord_n, rpm_p)[0],
             np.linspace(0, 1, 100),
             soft_sand(K=21.0 * 10**9, G=7.0 * 10**9, rho=2.58, phi=np.linspace(0, 1, 100), phic=phi_c, n=11.0, P=rpm_p)[0],
             color[0],
             color[2],
             code[0],
             code[2])

    df['Gsat'] = df['Gsoft']
    df['Ksat'] = np.empty(len(df['Gsat']))
    df['Ksat'][depth <= ow] = gassmann(df['Ksoft'][depth <= ow], rpm_k, k_oil, df['phi'][depth <= ow])
    df['Ksat'][depth > ow] = gassmann(df['Ksoft'][depth > ow], rpm_k, k_water, df['phi'][depth > ow])

    df['Msat'] = df['Ksat'] + (df['Gsat'] * 4./3.)
    df['Vp'] = (df['Msat']/10**9/df['rhob'])**0.5
    df['Vs'] = (df['Gsat']/10**9/df['rhob'])**0.5
    
    error_var = syn_variogram(depth, 1.0, 1.0)

    df['Vp'] = simulations_error(error_var, df['Vp'], 0.4)
    df['Vs'] = simulations_error(error_var, df['Vs'], 0.1)
    df['rhob'] = simulations_error(error_var, df['rhob'], 0.1)

    #print(np.shape(df.filter(regex='^fbool', axis=1).values))
    X, Y = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
    bivariate_plot(X, Y, rv_pdf, color)
    logplots(depth, df['vclay'], df['phi'], df['sw'], df['Vp'], df['Vs'], df['rhob'], code, df.filter(regex='^fbool', axis=1).values, color)
    crossplot_vpvs(df['Vp'], df['Vs'], z=df['vclay'], zlabel='VCLAY [dec]')

    #pd.DataFrame(data=np.array([depth, df['facies'], df['code_facies'], df['vclay'], df['phi'], df['sw'], df['Vp'], df['Vs'], df['rhob']]).T,
    #          columns=['DEPTH', 'CODE', 'FACIES', 'VCLAY', 'PHI', 'SW', 'VP', 'VS', 'RHOB']).to_csv('results/pseudowell.las',
    #                                                                                                 sep = ',', index=False)
    