# HdsaLib (Hyper-differential Sensitivity Analysis Library)

This C++ library implements various Hyper-differential Sensitivity Analysis algorithms, with a focus on model discrepancy analysis.

## File Structure

- [**src/**](./src/): Source code defining class interfaces and codes for analysis.
- [**examples/**](./examples/): Applications / examples with real physics.
- [**tests/**](./tests/): Unit tests and toy examples for verification.
- [**scripts/**](./scripts/): Example scripts for configuration

### Installation

HdsaLib depends on Trilinos; see "scripts/trilinos-configure-macos-clang" for an example of a Trilinos build.
Trilinos must be installed first and provided to the HdsaLib configuration script; see 
"scripts/hdsalib-configure-macos-clang" for an example of this. HdsaLib can be built with its own examples
and tests. Also, HdsaLib has an interface to MrHyDE. To use this, first build HdsaLib and then 
provide it to the MrHyDE build configuration; see "scripts/mrhyde-configure-macos-clang" for an example.

