# timevar-poly

A lightweight Python engine for **time-variant polynomials**.

## Features
- Coefficients can be constants or functions of time (`TimeFunction`).
- Automatic constant reduction (optimizes `0`, `1`, and constant operations).
- Symbolic-numeric differentiation (spatial `diff_x()` and temporal `diff_t()`).
- No external dependencies.

## Quick Start
```python
from timevar_poly import TimeVariantPolynomial, TimeFunction
import math

# P(x, t) = (sin t)x + 5
P = TimeVariantPolynomial({1: TimeFunction(math.sin), 0: 5})
print(P(x=2, t=math.pi/2)) # Output: 7.0
