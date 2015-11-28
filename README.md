# nbody
Tiny experimental high precision n-body simulation for small n.

## Features

* Uses RK4 integration with variable time step, determined in the first acceleration computation step.

* Optionally uses the KahanAdder class from OpenSAGE, available at http://opensage.googlecode.com/svn/core_code/trunk/src/c++/include/numerics/kahan.h.

  Kahan summation ought to greatly reduce the amount of roundoff error introduced by the many small time steps taken in close encounters,   besides the overall long term numerical drift due to the chaotic system.

* Orbit computation is done in a separate thread to allow smooth display while computing. Orbit data is protected with a mutex.

* Some very primitive camera tracking code.

## Planned additions

* Free camera, a la WASD + mouse.

* Arbitrary precision arithmetic, so paths can be "verified" by increasing the number of bits used and checking how long until the paths diverge.

* Export of path data in CSV and OBJ formats (either as path or via mesh generation).

* Lua scripting for setting up the initial conditions.
