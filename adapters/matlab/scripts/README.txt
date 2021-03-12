I need to compile Trilinos using the configure-trilinos included in this directory in order for the matlab adapater codes to link with Trilinos correctly.

I need to execute the set_library_paths script (using source) in order to put the correct matlab dynamic libraries in scope for the codes to execute.

compile_script contains the command to compile the matlab adapter examples. For each example I need to update the EXAMPLE_NAME in two places in compile_script
