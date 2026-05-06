import numpy as np

# Number of sensor points
num_sensors = 10

# Generate random points in (0,1)^2
xy = np.random.rand(num_sensors, 2)

# Append z = 0 for every point
z = np.zeros((num_sensors, 1))
points = np.hstack((xy, z))

# Write to file
np.savetxt("sensor_points.dat", points, fmt="%.8f")

print(f"Wrote {num_sensors} sensor points to sensor_points.dat")