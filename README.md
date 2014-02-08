AirfoilOpt
==========

A grid-less flow solver routine in MATLAB based on the Hess-Smith Panel Method which I coupled with an evolutionary optimization algorithm to optimize airfoil shapes for desired criteria.

Details
=======

Airfoil shapes are parametrized as per NACA 4-digit scheme. The solver takes in the 4 digit input, along with flow details such as Reynold's Number and angle of attack to give the Pressure/Velocity Distribution around the airfoil. This can be used to calculate some important specifications such as Lift, Drag and other derived quantities.

The code is written in such a way that alternative forms of parametrization can be used without changing the rest of the code. Same with the flow solvers. The code was written to be used as a cost function to an evolutionary optimization algorithm which aimed to generate ideal airfoil shapes for desired criteria as specified by the user. The user can also modify the cost function to suit the needs.
