import os
import numpy as np
import pyvista as pv

# ----------------------------------------------------------------------
# File names
# ----------------------------------------------------------------------
exodus_file = "output.exo"
points_file = "sensor_points.dat"
output_prefix = "sensor_data"
noise_level = 0.0  # percent

# Name of the vector field in the Exodus file
vector_field_name = "d"

# Output names for the three components
component_labels = ["dx", "dy", "dz"]

# Optional: set for reproducible noise
# np.random.seed(12345)

# ----------------------------------------------------------------------
# Read sensor points
# Expecting one point per line: x y z
# ----------------------------------------------------------------------
if not os.path.exists(points_file):
    raise FileNotFoundError(f"Points file not found: {points_file}")

points = np.loadtxt(points_file)

# Handle single-point file
if points.ndim == 1:
    if points.size != 3:
        raise ValueError(f"{points_file} must contain rows of x y z")
    points = points.reshape((1, 3))

if points.shape[1] != 3:
    raise ValueError(f"{points_file} must contain exactly 3 columns")

# ----------------------------------------------------------------------
# Read Exodus mesh
# ----------------------------------------------------------------------
if not os.path.exists(exodus_file):
    raise FileNotFoundError(f"Exodus file not found: {exodus_file}")

mesh = pv.read(exodus_file)

# Exodus may come in as a multiblock; combine if needed
if isinstance(mesh, pv.MultiBlock):
    mesh = mesh.combine()

# ----------------------------------------------------------------------
# Check vector field
# ----------------------------------------------------------------------
point_fields = list(mesh.point_data.keys())
cell_fields = list(mesh.cell_data.keys())

if vector_field_name not in point_fields:
    raise ValueError(
        f"Vector field '{vector_field_name}' not found in point data.\n"
        f"Available point fields: {point_fields}\n"
        f"Available cell fields: {cell_fields}"
    )

# ----------------------------------------------------------------------
# Probe/sample vector field at sensor points
# ----------------------------------------------------------------------
sensor_points = pv.PolyData(points)
sampled = sensor_points.sample(mesh)

if vector_field_name not in sampled.point_data:
    raise RuntimeError(f"Sampled output does not contain field '{vector_field_name}'")

values = np.asarray(sampled.point_data[vector_field_name], dtype=float)

if values.ndim != 2 or values.shape[1] < 3:
    raise ValueError(
        f"Field '{vector_field_name}' must have at least 3 components, "
        f"but got shape {values.shape}"
    )

# Keep only the first 3 components as x,y,z
values = values[:, :3]

# ----------------------------------------------------------------------
# Split into dx, dy, dz and write separate files with noise
# ----------------------------------------------------------------------
for i, field_name in enumerate(component_labels):
    component_values = values[:, i]
    output_file = f"{output_prefix}_{field_name}.dat"

    # Add mean-zero Gaussian noise
    avg_magnitude = np.mean(np.abs(component_values))
    sigma = (noise_level / 100.0) * avg_magnitude
    noisy_values = component_values + np.random.normal(
        loc=0.0,
        scale=sigma,
        size=component_values.shape
    )

    np.savetxt(output_file, noisy_values, header=field_name, comments="")

    print(
        f"Successfully wrote '{output_file}' "
        f"(noise sigma = {sigma:.6e}, noise_level = {noise_level}%)"
    )


# import os
# import numpy as np
# import pyvista as pv

# # ----------------------------------------------------------------------
# # File names
# # ----------------------------------------------------------------------
# exodus_file = "output.exo"
# points_file = "sensor_points.dat"
# output_prefix = "sensor_data"
# noise_level = 0.0  # percent
# field_names = ["dx", "dy", "dz"]

# # Optional: set for reproducible noise
# # np.random.seed(12345)

# # ----------------------------------------------------------------------
# # Read sensor points
# # Expecting one point per line: x y z
# # ----------------------------------------------------------------------
# if not os.path.exists(points_file):
#     raise FileNotFoundError(f"Points file not found: {points_file}")

# points = np.loadtxt(points_file)

# # Handle single-point file
# if points.ndim == 1:
#     if points.size != 3:
#         raise ValueError(f"{points_file} must contain rows of x y z")
#     points = points.reshape((1, 3))

# if points.shape[1] != 3:
#     raise ValueError(f"{points_file} must contain exactly 3 columns")

# # ----------------------------------------------------------------------
# # Read Exodus mesh
# # ----------------------------------------------------------------------
# if not os.path.exists(exodus_file):
#     raise FileNotFoundError(f"Exodus file not found: {exodus_file}")

# mesh = pv.read(exodus_file)

# # Exodus may come in as a multiblock; combine if needed
# if isinstance(mesh, pv.MultiBlock):
#     mesh = mesh.combine()

# # ----------------------------------------------------------------------
# # Check fields
# # ----------------------------------------------------------------------
# point_fields = list(mesh.point_data.keys())
# cell_fields = list(mesh.cell_data.keys())

# missing = [name for name in field_names if name not in point_fields]
# if missing:
#     raise ValueError(
#         f"These fields were not found in point data: {missing}\n"
#         f"Available point fields: {point_fields}\n"
#         f"Available cell fields: {cell_fields}"
#     )

# # ----------------------------------------------------------------------
# # Probe/sample fields at sensor points
# # ----------------------------------------------------------------------
# sensor_points = pv.PolyData(points)
# sampled = sensor_points.sample(mesh)

# for field_name in field_names:
#     if field_name not in sampled.point_data:
#         raise RuntimeError(f"Sampled output does not contain field '{field_name}'")

#     values = np.asarray(sampled.point_data[field_name], dtype=float)
#     output_file = f"{output_prefix}_{field_name}.dat"

#     # ------------------------------------------------------------------
#     # Add mean-zero Gaussian noise
#     # std dev = (noise_level/100) * average magnitude of the sampled data
#     # ------------------------------------------------------------------
#     avg_magnitude = np.mean(np.abs(values))
#     sigma = (noise_level / 100.0) * avg_magnitude
#     noisy_values = values + np.random.normal(loc=0.0, scale=sigma, size=values.shape)

#     # ------------------------------------------------------------------
#     # Write output
#     # Scalar field: one value per line
#     # Multi-component field: one row per point, multiple columns
#     # ------------------------------------------------------------------
#     if noisy_values.ndim == 1:
#         np.savetxt(output_file, noisy_values, header=field_name, comments="")
#     else:
#         header = " ".join(f"{field_name}_{i}" for i in range(noisy_values.shape[1]))
#         np.savetxt(output_file, noisy_values, header=header, comments="")

#     print(
#         f"Successfully wrote '{output_file}' "
#         f"(noise sigma = {sigma:.6e}, noise_level = {noise_level}%)"
#     )