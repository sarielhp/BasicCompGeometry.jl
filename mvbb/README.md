# MVBB 

This package was previously called gdiam.

## Whats in here?

This project contains source code in C++ for

- Approximating the minimum volume enclosing bounding box (MVBB) of a
  point set in 3d.

- Approximating the diameter.

- Ground breaking revolutionary code for computing convex-hull in 2d
  (joking, joking [but it does contain 2d convex-hull code])

The code seems to be fast and reasonably stable numerically. The code
is somewhat a weird mix of all style C++ and new style. At some point
it would be good to rewrite the code to move completely to new style
c++ (i.e., STL). It seems to compile/run perfectly under
g++/linux. Should work also on other operating systems without any
change.

I wrote the code 25 years ago, and it can definitely use more
documentation. If I have the energy I would write such doucmentation,
or hmm, get some AI tool to do so ;).

---------------------------------------------------------

## Refs 
Source code implementing the algorithms described in two papers

- A Practical Approach for Computing the Diameter of a Point-Set
  Sariel Har-Peled
  https://arxiv.org/abs/2505.11317

- Gill Barequet, Sariel Har-Peled: 
  Efficiently Approximating the Minimum-Volume Bounding Box of a
  Point Set in Three Dimensions. J. Algorithms 38(1): 91-109 (2001)

## Copyright

Copyright 2001 Sariel Har-Peled. This program is free software; you
can redistribute it and/or modify it under the terms of either:

 * the GNU General Public License as published by the Free Software
   Foundation; either version 2, or (at your option) any later version.

or

 * the GNU Lesser General Public License as published by the Free
   Software Foundation; either version 2.1, or (at your option)
   any later version.

or

 * MIT License https://opensource.org/licenses/MIT
 
If you need the source code under a different open-source license, let
me know - I don't really care.

Program is provided without any guarantees. Use at your own risk.

## A competing/followup program

Another implementation of this and similar algorithm is available here:
(https://github.com/gabyx/ApproxMVBB)

## Thanks 

    Many people over the years modified the code and added stuff, fixed stuff, etc. 
    See history below for details.

## Random comment

I am working on a Julia implementation of some parts of this code. I
have the diameter code working, in case anybody is interested.



## History

- 12/12/25
  + Rewrote some code, added some comments, and posted on github.
   
  + Added namespace MVBB for the package. Much cleaner than
     without. (Its a breaking change.)
     
  + Rewrote the 2d convex-hull code. Renamed some packages, added some
    minor documentation. 

- 3/6/18 
  + Tobias Stohr reported & fixed a bug wis missing constructor for ???

- 8/18/16 
  + Apply various code tweaks from BRL-CAD (C. Yapp)

- 12/22/14 
  + Apply GSoC patch switching the convex hull algorithm to
    monotone chain (P. Amidon)
 
- 8/7/13 
  + Add DLL import/export logic for Windows (C. Yapp)

- 8/4/13 
  + Add get_vertex method for low-level data translation. (C. Yapp)

- 8/3/13 -
  + Updated licensing - can now use either GPLv2 or LGPLv2.1.

- 3/28/01 
  + Original code updated to be more robust. It should now
    handle really abnoxious inputs well (i.e., points with equal
    coordinates, etc.


