from latex2sympy2 import latex2sympy
import sympy
from sympy import latex
def show(expr,title):
    e = sympy.simplify(expr)
    print(title,'\t',e,'\t',latex(e))

INPUT_STRING="""
x_1 - x_2 + 2x_1^2 +2x_1x_2 +x_2^2
"""

d = sympy.Symbol('d')
sym = latex2sympy(INPUT_STRING)
deriv = sympy.Matrix([sympy.diff(sym,'x_1'),sympy.diff(sym,'x_2')])
phi = sym.subs('x_1',d).subs('x_2',-d )

print (sym)
print(deriv)
print(sympy.latex(deriv.subs('x_1',0).subs('x_2',0)))
print("phi", phi,latex(phi))

nabla2 = deriv.subs({'x_1':-1,'x_2':1})
print(nabla2)

phi2 = sym.subs({'x_1':-1-d,'x_2':1-d})
phi2 = sympy.simplify(phi2)
print("phi2",'\t',phi2,'\t',latex(phi2))

x_d = sympy.Matrix([-1-d,1-d])
x_2 = x_d.subs('d',-0.1)
show(x_2,"x_fin")