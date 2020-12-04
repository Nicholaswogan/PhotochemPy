# PhotochemPy
PhotochemPy will be a photochemical model of terrestrial atmospheres. But its not done yet.

## History
PhotochemPy is an updated version of the `Atmos` photochemical model, originally developed by Jim Kasting and Kevin Zahnle and further developed by many of their students and colleges. The code was originally written in Fortran 77. This version of the code is written in Fortran 90. It has been restructured so that all varaibles are dynamically allocated. Also, it can be easily compliled as a Python package using `f2py` if you like to work in Python. Although, if you prefer Fortran and don't like Python, you can just use the code in fortran too. I'll include examples for both.

The code is structured so that the photochemical model can easily be evolved with different ODE integration methods. My plan is evolve the photochemical model with CVODE, so that the atmosphere is actually resolved in time.

I will also implement a simple backward Euler integration scheme. This will be used to integrate the system of ODEs to steady-state. Also, I will try to see if optimization methods can be used to find steady-state atmospheres.
