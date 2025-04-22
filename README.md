# Introduction

**TODO**

# Quickstart

From a Sage >= 10.2 installation, install the package with :
```Bash
pip install git+https://github.com/ea-icj/moment_cone
```
eventually prefixed by `sage -m`, then test it using for example the `moment_cone` script:
```Bash
$ moment_cone kron 4 4 4 1 --formats none
Configuration:
	representation: Kronecker
	N: [4, 4, 4, 1]
	[...]

MomentConeStep: ...                                                              {MomentConeStep}
    [...]
    BirationalityStep: ...                                                       {BirationalityStep}
        ineq_candidates: Dataset(#pending=0, #validated=47)                      {BirationalityStep}                                                                                                               
        Done (Wall: 2561.108ms, CPU: 2561.020ms (100%))                          {BirationalityStep}
    Done (Wall: 4186.895ms, CPU: 4186.689ms (100%))                              {MomentConeStep}
```
or using the Python/Sage interpreter:
```Python
>>> from moment_cone import *
>>> V = KroneckerRepresentation((4, 4, 4, 1))
>>> ineqs = moment_cone(V)
MomentConeStep: ...                                                              {MomentConeStep}
    [...]
    Done (Wall: 4224.499ms, CPU: 4224.381ms (100%))                              {MomentConeStep}
>>> ineqs
Dataset(#pending=0, #validated=47)
```
# Testing

First install the development packages by specifying the optional `[dev]` dependency:
```Bash
pip install git+https://github.com/ea-icj/moment_cone[dev]
```

Then :
- for the unit tests:
    ```bash
    make unittest
    ```
- testing the examples from the documentation:
    ```bash
    make doctest
    ```
- static type checker:
    ```bash
    make mypy
    ```
or all these tests by simply removing the target:
```bash
make
```
