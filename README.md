# 2PJIT

Two-phase 3D Jet Instability Tool in Cylindrical Coordinates

This tool consists of 3 matlab scripts. 'twoPJIT.m' is the script is run by the user. The properties of the jet are entered in this script along with other inputs which are prompted while running the script. This script calls 'Cylindrical_3D_solution.m' to solve the discretized equations for each input frequency. The 'twoPJIT.m' script also calls 'perturbation.m' to solve for perturbation flow fields for the frequency corresponding to the highest growth rate.

For further description, refer to Documentation.pdf.
