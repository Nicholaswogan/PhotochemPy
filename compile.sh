



f2py -c Photochem.f90 -m Photochem only: allocate_memory right_hand_side \
jacobian read_species read_reactions read_atmosphere photgrid rates gridw readflux \
initphoto initmie
