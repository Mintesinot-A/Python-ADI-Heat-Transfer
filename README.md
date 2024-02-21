ADI Scheme Implementation in Python
This repository contains a Python implementation of the Alternating-Direction Implicit (ADI) scheme for solving two-dimensional parabolic partial differential equations (PDEs). 
The specific implementation is based on the example presented in Steven C. Chapra's book "Numerical Methods for Engineers: Finite Difference: Parabolic Equations".
Purpose:
This code provides a straightforward way to apply the ADI scheme to a specific heat transfer problem with given boundary conditions. 
By adjusting the input parameters and initial conditions, the code can be adapted to solve other similar parabolic PDEs.
How to use:
Install required libraries: No external libraries are necessary for this code to run.
Edit parameters: Modify the following variables in the adi_scheme.py file to match your specific problem:
initial_temp: Initial temperature distribution inside the domain (excluding boundaries).
r_temp, l_temp, top_temp, bottom_temp: Boundary temperatures at right, left, top, and bottom, respectively.
time_step: Time step for the numerical solution.
step_size: Spatial step size for the grid discretization.
length: Length of the domain in the spatial direction.
time: Total simulation time.
diffusivity: Thermal diffusivity of the material.

Run the code: Execute the adi_scheme.py file. The output will be a list of temperature distributions at each time step, stored in the temp_over_time variable.
Further notes:
The code currently solves a specific case with fixed boundary conditions. It can be extended to handle more general scenarios by modifying the boundary condition 
calculations within the code.
This is a basic implementation and may not be optimized for speed or memory efficiency. For more complex problems, consider exploring optimized libraries or specialized software packages.
This is just a basic example, and you can customize it further to fit your specific needs and provide more detailed instructions or usage examples. 
