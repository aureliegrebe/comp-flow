# comp-flow

## About

This library contains functions for basic compressible flow relations. 
The naming of the API is heavily inspired by the python library [`compflow`](https://github.com/jb753/compflow).

The included functions have no input checking or error handling whatsoever.
Invalid (non-physical) inputs such as mach < 1 for a shock relation or gamma < 1
may produce non-sensical outputs.

## Features

 - Isentropic relations calculated from Mach numbers.
 - Inverse isentropic relations calculating Mach numbers.
 - Normal shock relations.
 - Weak oblique shock relations.

## To Do

 - Fanno Flow
 - Rayleigh Flow
 - Oblique Shock Relations

## Contributing

Contributions are welcome, particularily optimizing the functions that make
use of Newton's method. If a particulary relation that you use is not included,
feel free to add it.
