

# improvements


- Increase number of arguments in each subroutine, to remove as much global data as possible.
- I think I will treat T as something that can not change after being read in. This is necesary to allow whole parallel runs. Otherwise, parallel reading of data files would need to happen, which sounds hard.
- 
