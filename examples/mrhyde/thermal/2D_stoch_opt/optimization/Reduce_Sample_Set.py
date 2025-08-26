import numpy as np

print("Reducing the sample set to active samples")

tau = 1e-4  
weights = np.loadtxt('sample_weights.dat')
sample_set = np.loadtxt('sample_set.dat')

mask = weights >= tau
reduced_weights = weights[mask]
reduced_sample_set = sample_set[mask]

reduced_weights = reduced_weights/reduced_weights.sum()

np.savetxt('sample_weights.dat', reduced_weights)
np.savetxt('sample_set.dat', reduced_sample_set)
