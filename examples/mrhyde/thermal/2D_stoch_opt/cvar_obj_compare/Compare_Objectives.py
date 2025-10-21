import numpy as np
import pdb

alpha = 0.5
beta = 0.7

lofi_objs = np.loadtxt('sample_obj_vals_lofi.dat')
N = lofi_objs.shape[0]
T = int(np.round(N*(1-beta)))
lofi_objs.sort()
lofi_obj = alpha * lofi_objs[-T:].mean() + (1-alpha) * lofi_objs.mean()

update_objs = np.loadtxt('sample_obj_vals_mean_update.dat')
update_objs.sort()
update_obj = alpha * update_objs[-T:].mean() + (1-alpha) * update_objs.mean()

hifi_objs = np.loadtxt('sample_obj_vals_hifi.dat')
hifi_objs.sort()
hifi_obj = alpha * hifi_objs[-T:].mean() + (1-alpha) * hifi_objs.mean()

filename = 'cvar_objectives.txt'
with open(filename, 'w') as file:
    file.write(f'CVaR objective at low-fidelity: {lofi_obj}\n')
    file.write(f'CVaR objective at mean update: {update_obj}\n')
    file.write(f'CVaR objective at high-fidelity: {hifi_obj}\n')

