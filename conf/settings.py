import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


OUTPUT = {
    'synthetic': os.path.join(BASE_DIR, 'synthetic.las'),      # path of synthetic
    'figures': os.path.join(BASE_DIR, 'figures')        # path of data figures
}

SIMUL_PARAS = {
    'top': 2000,                # top of interval
    'base': 2200,               # base of interval
    'sr': 0.5,                   # sample rate
    'nfacies': 3,               # number of facies (only works for 3)
    'p11': 0.9,                # facies 1 thickness
    'p12': 0.05,                # probability of transition from facies 1 to 2   
    'p13': 0.05,
    'p22': 0.93,
    'p21': 0.0,
    'p23': 0.07,
    'p33': 0.95,
    'p31': 0.05,
    'p32': 0.0,
    'code': ['Shale', 'Silty', 'Sand'],
    'color': ['green', 'brown', 'yellow']
}

FACIES_PARAS = {
    'mux': [0.05, 0.18, 0.25],                      # porosity average [Mean_f1, Mean_f2, Mean_f3]                                                  
    'muy': [0.5, 0.2, 0.1],                         # clay volume average
    'covf1': [[0.002, -0.002], [-0.002, 0.02]],     # covariance matrix for facies 1
    'covf2': [[0.0013, -0.001], [-0.001, 0.0032]],
    'covf3': [[0.0019, -0.0035], [-0.0035, 0.01]]  
}

STOCH_PARAS = {
    'lh1': 7.5,                  # correlation range for facies 1
    'lh2': 7.5,
    'lh3': 7.5
}

FINAL_PARAS = {
    'ow': 2100.0,                  # depth of oil-water contact
    'so': 0.9,                   # Oil saturation
    'rho_oil': 0.63,
    'k_oil': 0.4 * 10**9,
    'rho_water': 1.01,
    'k_water': 2.18 * 10**9,
    'rho_min': 2.65,
    'rpm': 'soft_sand',
    'rpm_p': 25 * 10**6,
    'rpm_k': 36 * 10**9,
    'rpm_g': 45 * 10**9, 
    'phi_c': 0.4,
    'coord_n': 4.0,
    'c_error': 0.5,
    'r_error': 1.0
}