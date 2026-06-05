# OedLib

A C++ library implementing optimal experimental design (OED) for Bayesian inverse problem.

## File Structure

### Installation

`OedLib` depends on the C++ header library
[`Eigen`](https://libeigen.gitlab.io/).
Follow the instructions to install `Eigen` above and then create a symbolic link to the `Eigen` directory within the `Eigen` project called `Eigen` inside the include directory.
This can be accomplished with the command `ln -s /path/to/eigen/Eigen /path/to/OedLib/include/Eigen`.
