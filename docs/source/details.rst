How photochemical models work
=============================

Here I roughly describe how photochemical models work. I follow Appendix B in `Catling and Kasting 2017 <https://www.cambridge.org/us/academic/subjects/physics/planetary-systems-and-astrobiology/atmospheric-evolution-inhabited-and-lifeless-worlds?format=HB&isbn=9780521844123>`_.

In a nut shell
^^^^^^^^^^^^^^

PhotochemPy is a one-dimensional (1-D) steady-state photochemical model. It is 1-D because it assumes the atmosphere has no horizontal heterogeneity, and that gas composition only changes with altitude. The model is steady-state because it seeks a state of photochemical equilibrium, i.e. where production of each gas at each altitude is equal to loss mechanisms. PhotochemPy can not track the time-evolution of an atmosphere.

PhotochemPy is governed by a large system of ordinary differential equations (ODEs) which define the production rate and loss rate of each species at each altitude due to chemistry to transport. Photochemical equilibrium is found by starting with some given atmospheric composition, then evolving the atmosphere forward in time according to the ODEs. After evolving the atmosphere for several billion years, we assume a steady-state is achieved.

Photochemical Model Equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The photochemical model is governed by a system of partial differential equations (PDEs) of the form:

.. math:: \frac{dn_i}{dt}= P_i - l_i n_i - \frac{d\Phi}{dz}
   :label: cont

Where,

.. math:: \Phi_i = Kn\frac{df_i}{dz}-D_in\left(\frac{1}{n_i} \frac{dn_i}{dz} + \frac{1}{H_i} + \frac{1+\alpha_{Ti}}{T}\frac{dT}{dz}\right)
   :label: flux

Here,

- :math:`t` = time (seconds)
- :math:`z` = altitude (cm)
- :math:`n_i` = number density (:math:`\text{molcules/cm}^3`) of species :math:`i`
- :math:`l_i` = chemical loss frequency (1/s)
- :math:`P_i` = chemical production rate (:math:`\mathrm{molcules/cm^3/s}`)
- :math:`\Phi_i` = flux of species :math:`i` (:math:`\mathrm{molcules/cm^2/s}`)
- :math:`f_i` = :math:`n_i/n` = mixing ratio
- :math:`D_i` = diffusion coefficient between species :math:`i` and the background atmosphere
- :math:`K` = eddy diffusion coefficient (:math:`\mathrm{cm^2/s}`)
- :math:`H_i` = :math:`kT/m_ig` = scale height of species :math:`i`
- :math:`\alpha_{Ti}` = thermal diffusion coefficient of species :math:`i` with respect to the background atmosphere
- :math:`k` = Boltzmann's constant (cgs units)
- :math:`m_i` = molecular mass of species :math:`i` (g/molecule)
- :math:`m` = mean molecular mass (g/molecule)
- :math:`g` = gravitational acceleration (:math:`\mathrm{cm/s^2}`)

The molecullar diffusion term in Equation :eq:`flux` is only valid for a minor constituent diffusing through a background atmosphere (e.g. :math:`\mathrm{H_2}` diffusion through  :math:`\mathrm{N_2}` in the modern Earth atmosphere).

For convenience, we recast Equation :eq:`cont` in terms of mixing ratios:

.. math::

   \frac{df_i}{dz} = \frac{1}{n}\frac{dn_i}{dz}-\frac{n_i}{n^2}\frac{dn}{dz}

Differentiation the ideal gas law (:math:`p=nkT`) with respect to :math:`z` gives:

.. math::

   \frac{dp}{dz} = k \left(n  \frac{dT}{dz}+T \frac{dn}{dz} \right)

Dividing each side by :math:`nkT`:

.. math::

   \frac{dp}{dz} = k \left(n  \frac{dT}{dz}+T \frac{dn}{dz} \right)
