

# improvements

- Put global data in different modules organized by how much they might change. Then explicitly bring into scope all global data in every subroutine (completed)
- Increase number of arguments in each subroutine, to remove as much global data as possible. Ideally would remove `photochem_wrk` and `photochem_vars`. This should only be done for "lower" subroutines. Subroutines like `right_hand_side` will need a bunch of global data.


We need a subroutine that checks for nonlinearities in SL species!



notes:
- initphoto - depends on usol via the column depth for O2 and NO xsection
- rainout - depends on usol, T, den, fco2
- lightning - DO THIS NEXT


NEXT, lets 
- add to read_reactions so that SL can not photolyze.
- Write subroutine which determines the dimension of the problem

Following catagories

1. Stuff that never changes after run once
  - read_species
  -

2. Depends on T, Den, edd, other things in `photochem_data`

3. Depends on things that change during integration (`photochem_wrk`)
