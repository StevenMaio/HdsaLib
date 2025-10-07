import numpy as np
import pdb

lofi_objs = np.loadtxt('sample_obj_vals_lofi.dat')
lofi_objs.sort()
lofi_obj = 0.5 * lofi_objs[-4:].mean() + 0.5 * lofi_objs.mean()

update_objs = np.loadtxt('sample_obj_vals_mean_update.dat')
update_objs.sort()
update_obj = 0.5 * update_objs[-4:].mean() + 0.5 * update_objs.mean()

hifi_objs = np.loadtxt('sample_obj_vals_hifi.dat')
hifi_objs.sort()
hifi_obj = 0.5 * hifi_objs[-4:].mean() + 0.5 * hifi_objs.mean()

filename = 'cvar_objectives.txt'
with open(filename, 'w') as file:
    file.write(f'CVaR objective at low-fidelity: {lofi_obj}\n')
    file.write(f'CVaR objective at mean update: {update_obj}\n')
    file.write(f'CVaR objective at high-fidelity: {hifi_obj}\n')

