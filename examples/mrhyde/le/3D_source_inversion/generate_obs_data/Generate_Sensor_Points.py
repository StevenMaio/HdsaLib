import numpy as np

# Number of sensor points
num_sensors = 100

xy = np.random.rand(num_sensors, 2)
xy[:,1] = .4 * xy[:,1]
z = np.ones((num_sensors, 1))
points0 = np.hstack((xy, 0*z))
points1 = np.hstack((xy, .1*z))
points = np.vstack((points0,points1))

# Write to file
np.savetxt("sensor_points.dat", points, fmt="%.8f")

print(f"Wrote {num_sensors} sensor points to sensor_points.dat")