# 📈 timevar-poly

**timevar-poly** is a lightweight Python engine designed to manipulate polynomials whose coefficients vary over time. It bridges the gap between symbolic structure and numerical evaluation, allowing for complex algebraic operations and derivatives.

## ✨ Key Features

* **Dynamic Coefficients:** Use any Python function (sin, cos, exp, etc.) or constants as polynomial coefficients.
* **Hybrid Computation:** Combines the structural benefits of polynomials with the flexibility of time-dependent functions.
* **Partial Evaluation (Currying):** Intelligently "freeze" either the spatial variable ($x$) or the temporal variable ($t$) to obtain specialized objects.
* **Advanced Derivatives:** Full support for spatial derivatives ($\frac{\partial F}{\partial x}$), temporal derivatives ($\frac{\partial F}{\partial t}$), and mixed partial derivatives.

## 🚀 Installation

```bash
pip install timevar-poly
```

## 🛠 Quick Start

### 1. Define a Polynomial
Let's define $F(x, t) = \cos(t)x + \sin(t)$:

```python
import math
from timevar_poly import TimeFunction, TimeVariantPolynomial

F = TimeVariantPolynomial({
    1: TimeFunction(math.cos),
    0: TimeFunction(math.sin)
}, label="F")

print(F) # Output: (cos(t))x + (sin(t))
```

### 2. The Power of Currying
You can evaluate one dimension while keeping the other as a variable:

```python
# Fix time: Returns a standard Polynomial in x
P_at_t1 = F(t=1)
print(P_at_t1) # Outputs numerical coefficients evaluated at t=1

# Fix space: Returns a pure TimeFunction of t
f_at_x2 = F(x=2)
print(f_at_x2) # Output: F(2, t)
print(f_at_x2(0)) # Calculates the value at t=0
```

### 3. Calculus
```python
# Derivative with respect to x
df_dx = F.diff_x()

# Derivative with respect to time (t)
df_dt = F.diff_t()

# Mixed second-order derivative (Schwarz's Theorem)
d2f_dtdx = F.diff_x().diff_t()
d2f_dxdt = F.diff_t().diff_x()

assert math.isclose(d2f_dtdx(2, 1), d2f_dxdt(2, 1), rel_tol=1e-9) # Should be True
```


## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👤 Author
**Ifèdé Assogba** - *Creator of the timevar-poly engine*
