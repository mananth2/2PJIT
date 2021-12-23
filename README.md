# 2PJIT

Two-phase 3D Jet Instability Tool in Cylindrical Coordinates

This tool consists of 3 main matlab scripts. 'twoPJIT.m' is the script is run by the user. The properties of the jet are entered in this script along with other inputs which are prompted while running the script. This script calls 'Cylindrical_3D_solution.m' to solve the discretized equations for each input frequency. The 'twoPJIT.m' script also calls 'perturbation.m' to solve for perturbation flow fields for the frequency corresponding to the highest growth rate.

For further description, refer to Documentation.pdf.

'Validation_Lin_Gordillo.m' is a separate script that has to be run separately. This scripts outputs the dispersion plots which can be validated with the plots in the following papers; Lin and Chen (Journal of Fluid Mechanics, 1998) and Gordillo and P-Saborid (Journal of Fluid Mechanics, 2005).
