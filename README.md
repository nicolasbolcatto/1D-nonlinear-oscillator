# 1D-nonlinear-oscillator
MATLAB Code from my Master's Thesis Report

This is the code I wrote and used to create a 1D-nonlinear oscillator model, designed to simulate the behaviour of a soil sample at the Resonant Column Test. This model was part of my Master's Thesis Report to fulfill the requirements for my Master's Degree in Geotechnical Engineering, back in 2019. You can check the full report if you like, it's included in the repository.

The model is based on input variables such as the initial shear modulus of the soil, sample height and diameter, soil density, load timing function, load amplitude and more. The model takes this inputs and iterates by modifying soil shear modulus and hysteretical damping ratio with a Newton-Raphson scheme. Time integration is performed via Newmark finite difference method with mean acceleration.

Code was implemented in MATLAB, as it was the propietary language used at the University. I plan to refactor it to Python or maybe JavaScript in the future, and perhaps work out a UI to make it more appealing. Feel free to contact me if you want to take part in the quest!

