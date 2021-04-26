

# improvements

- Put global data in different modules organized by how much they might change. Then explicitly bring into scope all global data in every subroutine (completed)
- Increase number of arguments in each subroutine, to remove as much global data as possible. Ideally would remove `photochem_wrk` and `photochem_vars`. This should only be done for "lower" subroutines. Subroutines like `right_hand_side` will need a bunch of global data.
