# Introduction

**TODO**

# Installation

The **first installation method** is starting from a [SageMath](https://www.sagemath.org/) >= 10.2 [installation](https://doc.sagemath.org/html/en/installation/index.html). It can be done using [`conda`](https://www.anaconda.com/docs/getting-started/miniconda/main):
```bash
conda create -n moment_cone -c conda-forge "sage>=10.2" "python>=3.11"
```
or using the given environment file:
```bash
conda env create -f sage_environment.yml
```

Then the `moment_cone` package can be installed using `pip`:
```bash
pip install git+https://github.com/ea-icj/moment_cone
```
eventually prefixed by `sage -m`.

The **second installation method** relies on a modularized fork of SageMath, named [passagemath](https://github.com/passagemath/passagemath), that only requires a Python>=3.11 installation (eg `conda create -n moment_cone "python=3.11"`) and the following `pip` command:
```bash
pip install "moment_cone[passagemath] @ git+https://github.com/ea-icj/moment_cone"
```

**TODO:** PyPI

# Quickstart

Once installed, the package can be tested using the `moment_cone` **command-line** script:
```bash
$ moment_cone kron 4 4 4 --formats py
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
with computed inequalities in generated `ineq_Kronecker_4_4_4.py` Python file (Normaliz format also available).


It can also be done through the Python/Sage **interpreter**:
```python
>>> from moment_cone import *
>>> V = KroneckerRepresentation((4, 4, 4, 1))
>>> ineqs = moment_cone(V)
MomentConeStep: ...                                                              {MomentConeStep}
    [...]
    Done (Wall: 4224.499ms, CPU: 4224.381ms (100%))                              {MomentConeStep}
>>> ineqs.validated()
[Inequality(tau  = 0 0 0 0 | 0 0 0 0 | 1 1 1 0 | -1,
            w    = 0 1 2 3 | 0 1 2 3 | 0 1 2 3 | 0,
            wtau = 0 0 0 0 | 0 0 0 0 | 1 1 1 0 | -1),
 Inequality(tau  = 1 0 0 0 | 1 0 0 0 | 1 1 1 0 | -2,
            w    = 0 1 2 3 | 0 1 2 3 | 1 2 3 0 | 0,
            wtau = 1 0 0 0 | 1 0 0 0 | 0 1 1 1 | -2),
 ...]
```

Computation can also be done in **parallel** on a multi-core computer:
```bash
$ moment_cone kron 5 5 5 --formats py --parallel
Configuration:
	representation: Kronecker
	N: [5, 5, 5, 1]
	[...]

MomentConeStep: ...                                                              {MomentConeStep}
    [...]
    BirationalityStep: ...                                                       {BirationalityStep}
        ineq_candidates: ListDataset(#pending=0, #validated=462)                 {BirationalityStep}
        Done (Wall: 24563.861ms, CPU: 274287.168ms (1117%))                      {BirationalityStep}
    [...]
    Done (Wall: 42339.791ms, CPU: 382788.186ms (904%))                           {MomentConeStep}                             {MomentConeStep}
```
or from a Python interpreter:
```python
>>> from moment_cone import *
>>> Parallel.configure()
>>> V = KroneckerRepresentation((4, 4, 4, 1))
>>> ineqs = moment_cone(V)
MomentConeStep: ...                                                              {MomentConeStep}
    [...]
Done (Wall: 43905.266ms, CPU: 388104.087ms (884%))                           {MomentConeStep}
>>> ineqs.validated()
[Inequality(tau  = 1 0 0 0 0 | 1 0 0 0 0 | 1 1 1 1 0 | -2,
            w    = 0 1 2 3 4 | 0 1 2 3 4 | 1 2 3 4 0 | 0,
            wtau = 1 0 0 0 0 | 1 0 0 0 0 | 0 1 1 1 1 | -2),
 Inequality(tau  = 1 0 0 0 0 | 1 0 0 0 0 | 1 1 1 1 0 | -2,
            w    = 0 1 2 3 4 | 2 0 1 3 4 | 0 1 3 4 2 | 0,
            wtau = 1 0 0 0 0 | 0 0 1 0 0 | 1 1 0 1 1 | -2),
 ...]
```

# Testing

The unit tests can be launched using the `dev` optional dependency:
```bash
pip install "moment_cone[dev] @ git+https://github.com/ea-icj/moment_cone"
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
