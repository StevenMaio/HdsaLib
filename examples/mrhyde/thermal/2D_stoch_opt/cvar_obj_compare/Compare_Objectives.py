import numpy as np
import pdb

lofi_objs = np.loadtxt('sample_obj_vals_lofi.dat')
lofi_objs.sort()
lofi_cvar = lofi_objs[-4:].mean()

update_objs = np.loadtxt('sample_obj_vals_mean_update.dat')
update_objs.sort()
update_cvar = update_objs[-4:].mean()

hifi_objs = np.loadtxt('sample_obj_vals_hifi.dat')
hifi_objs.sort()
hifi_cvar = hifi_objs[-4:].mean()

filename = 'cvar_objectives.txt'
with open(filename, 'w') as file:
    file.write(f'CVaR objective at low-fidelity: {lofi_cvar}\n')
    file.write(f'CVaR objective at mean update: {update_cvar}\n')
    file.write(f'CVaR objective at high-fidelity: {hifi_cvar}\n')

