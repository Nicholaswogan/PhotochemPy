Welcome to PhotochemPy's documentation
======================================

PhotochemPy is a 1-D photochemical model of rocky planet's atmospheres. Given inputs, like the stellar UV flux, the atmospheric temperature structure, etc., this code will find the steady-state chemical composition of an atmosphere.

Basic Usage
^^^^^^^^^^^

.. code-block:: Python

   from PhotochemPy import PhotochemPy # import

   # Load files defining your problem
   pc = PhotochemPy('species.dat', \
                    'reactions.rx', \
                    'planet.dat', \
                    'input_photchem.dat', \
                    'atmosphere.txt', \
                    'star_flux.txt')
   pc.integrate() # Finds photochemical equilibrium

   # Output atmospheric mixing ratio
   out = pc.out_dict()

   # print ground level methane mixing ratio
   print('CH4 concentration =','%.2e'%out['CH4'][0])

A more complete example is available in the Tutorial linked in the User Guide.

.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   install
   Tutorial
   PhotochemPy
   todo
