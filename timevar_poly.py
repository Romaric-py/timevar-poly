"""
timevar-poly: A lightweight engine for time-variant polynomials.
Supports algebraic operations and temporal/spatial derivatives.
"""

__version__ = "0.1.1"
__author__ = "Ifèdé Assogba"


import inspect
import operator


class TimeFunction:
    def __init__(self, func):
        sig = inspect.signature(func)
        if len(sig.parameters) != 1:
            raise ValueError("TimeFunction must be a one-parameter function")
        self.func = func
        self._made_from_constant = False
        self._constant_value = None

    def __repr__(self):
        return f"TimeFunction({self.func.__name__})"

    def __str__(self):
        return self.func.__name__

    def __call__(self, x):
        return self.func(x)

    @classmethod
    def from_constant(cls, value):
        def const_func(x):
            return value

        const_func.__name__ = str(value)
        time_function = cls(const_func)
        time_function._made_from_constant = True
        time_function._constant_value = value
        return time_function

    def _ensure_time_func(self, obj):
        if isinstance(obj, (int, float)):
            return self.from_constant(obj)
        if isinstance(obj, TimeFunction):
            return obj
        raise TypeError(f"Cannot operate with {type(obj)}")

    def _combine(self, other, op, op_symbol):
        other = self._ensure_time_func(other)

        if self._made_from_constant and other._made_from_constant:
            return self.from_constant(op(self._constant_value, other._constant_value))

        def combined_func(x):
            return op(self.func(x), other.func(x))

        combined_func.__name__ = (
            f"({self.func.__name__}{op_symbol}{other.func.__name__})"
        )
        return TimeFunction(combined_func)

    def _is_constant(self, obj, value):
        return (
            isinstance(obj, (int, float))
            and obj == value
            or (
                isinstance(obj, TimeFunction)
                and obj._made_from_constant
                and obj._constant_value == value
            )
        )

    def __add__(self, other):
        if self._is_constant(other, 0):
            return self
        return self._combine(other, operator.add, "+")

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self._combine(other, operator.sub, "-")

    def __rsub__(self, other):
        return self.from_constant(other).__sub__(self)

    def __mul__(self, other):
        if self._is_constant(other, 0):
            return self.from_constant(0)
        if self._is_constant(other, 1):
            return self
        return self._combine(other, operator.mul, "*")

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if self._is_constant(other, 0):
            raise ZeroDivisionError("Cannot divide by zero")
        if self._is_constant(other, 1):
            return self
        return self._combine(other, operator.truediv, "/")

    def __rtruediv__(self, other):
        return self.from_constant(other)._combine(self, operator.truediv, "/")

    def __neg__(self):
        if self._is_constant(self, 0):
            return self

        def neg_func(x):
            return -self.func(x)

        neg_func.__name__ = f"(-{self.func.__name__})"
        return TimeFunction(neg_func)

    def __pow__(self, power):
        if power == 0:
            return self.from_constant(1)

        if power == 1:
            return self

        if self._is_constant(self, 0) and power > 0:
            return self.from_constant(0)
        if self._is_constant(self, 1):
            return self.from_constant(1)

        def pow_func(x):
            return self.func(x) ** power

        pow_func.__name__ = f"({self.func.__name__}**{power})"
        return TimeFunction(pow_func)

    def __rpow__(self, base):
        def rpow_func(x):
            return base ** self.func(x)

        rpow_func.__name__ = f"({base}**{self.func.__name__})"
        return TimeFunction(rpow_func)

    def derivative(self, h=1e-5):
        if self._made_from_constant:
            return self.from_constant(0)

        def deriv_func(x):
            return (self.func(x + h) - self.func(x - h)) / (2 * h)

        deriv_func.__name__ = f"derivative({self.func.__name__})"
        return TimeFunction(deriv_func)


class TimeVariantPolynomial:
    def __init__(self, data):
        self.data = {deg: v for deg, v in data.items() if self._is_not_zero(v)}

    def _is_not_zero(self, val):
        if isinstance(val, (int, float)):
            return val != 0
        return True

    def __repr__(self):
        if not self.data:
            return "0"

        terms = []
        for deg in sorted(self.data.keys(), reverse=True):
            coeff = self.data[deg]
            c_str = f"({coeff}(t))" if isinstance(coeff, TimeFunction) else str(coeff)

            if deg == 0:
                terms.append(c_str)
            elif deg == 1:
                terms.append(f"{c_str}x")
            else:
                terms.append(f"{c_str}x^{deg}")
        return " + ".join(terms).replace("+ -", "- ")

    def __call__(self, x, t=None):

        res = 0
        for deg, coeff in self.data.items():
            c_val = coeff(t) if isinstance(coeff, TimeFunction) else coeff
            res += c_val * (x**deg)
        return res

    def __add__(self, other):
        if isinstance(other, (int, float)):
            other = TimeVariantPolynomial({0: other})

        all_degs = set(self.data.keys()) | set(other.data.keys())
        new_data = {}
        for d in all_degs:
            v1 = self.data.get(d, 0)
            v2 = other.data.get(d, 0)
            new_data[d] = v1 + v2
        return TimeVariantPolynomial(new_data)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return TimeVariantPolynomial({d: v * other for d, v in self.data.items()})

        new_data = {}
        for deg1, c1 in self.data.items():
            for deg2, c2 in other.data.items():
                d_res = deg1 + deg2

                new_data[d_res] = new_data.get(d_res, 0) + (c1 * c2)
        return TimeVariantPolynomial(new_data)

    def diff_x(self):
        new_data = {d - 1: v * d for d, v in self.data.items() if d > 0}
        return TimeVariantPolynomial(new_data)

    def diff_t(self):
        new_data = {}
        for d, v in self.data.items():
            v_func = v if isinstance(v, TimeFunction) else TimeFunction.from_constant(v)
            deriv = v_func.derivative()

            if not (deriv._made_from_constant and deriv._constant_value == 0):
                new_data[d] = deriv

        return TimeVariantPolynomial(new_data)


if __name__ == "__main__":

    def identity(x):
        return x

    def sqrt(x):
        return x**2

    def cube(x):
        return x**3

    F = TimeVariantPolynomial({2: TimeFunction(sqrt), 1: 3, 0: TimeFunction(identity)})

    print("F(x, t) =", F)
    print("F(2, 0) =", F(2, 0))  # sqrt(0) * 2^2 + 3*2 + 0 = 0*4 + 6 + 0 = 6
    print("F(2, 1) =", F(2, 1))  # sqrt(1) * 2^2 + 3*2 + 1 = 1*4 + 6 + 1 = 11
    print("dF/dx =", F.diff_x())  # 2*sqrt(t)*x + 3
    print("dF/dt =", F.diff_t())  # (1/(2*sqrt(t)))*x^2 + 1

    G = TimeVariantPolynomial({1: TimeFunction(cube), 0: 5})

    print("G(x, t) =", G)  # cube(t)*x + 5
    print("F + G =", F + G)  # (sqrt(t)*x^2 + (3 + cube(t))*x + (identity(t) + 5))
    print(
        "F * G =", F * G
    )  # (sqrt(t)*cube(t)*x^3 + (3*cube(t) + 5*sqrt(t))*x^2 + (15 + cube(t)*identity(t))*x + 5*identity(t))

    a = TimeFunction.from_constant(2)
    b = TimeFunction.from_constant(3)
    print("a =", a)  # 2
    print("b =", b)  # 3
    print("a + b =", a + b)  # 5
    print("da/dt =", a.derivative())  # 0
    print("a * b =", a * b)  # 6
