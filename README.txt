The repository is organized as follows.

Under the "Old" folder is an old version of the codes which we are keeping in the repository for now to aid the migration of some old code into the current version.

We currently do not have a CMake facility to build hdsalib. Rather, we are utilizing the CMake capabilities of the application driving the analysis. Under "src" we have the hdsalib source code, all of which are ".hpp" files.

Under "examples" we have subfolders corresponding to the applications driving the analysis. For instance, under "hdsalib/examples/rol" we have examples that execute from within ROL. This is done by creating a symbolic link to place folders under "hdsalib/examples/rol" within the directory "Trilinos/packages/rol/examples/". From there, we use the "CMakeList" within "Trilinos/packages/rol/examples/" to compile the example. All hdsa source code is included via the "#include "../../../src/source_file.hpp" within the example file.

Under "tests" we have some basic regression tests. As with the ROL examples above, we create a symbolic link to place the "hdsalib/tests" folder under "Trilinos/packages/rol/examples/" and utilize the ROL's CMakeList to include the tests as examples.

As we develop examples within MrHyDE, we will keep the interface source code in "hdsalib/interfaces/mrhyde" and the examples in "hdsalib/examples/mrhyde". We will follow a similar convention with Albany.
