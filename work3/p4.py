from latex2sympy2 import latex2sympy
import sympy
from sympy import latex,Matrix,Rational
from dataclasses import dataclass,field
def show(expr,title):
    e = sympy.simplify(expr)
    print(title,'\t',e,'\t',latex(e))

x,y,d = sympy.symbols("x y d")
expr = 2*x**2 + 2*y**2 -2*x*y - 4*x -6*y

nx = sympy.diff(expr,x)
ny = sympy.diff(expr,y)

show(nx,"nx")
show(ny,"ny")

grad = Matrix([nx,ny])

g0 = grad.subs({x:0,y:0})
show(g0,"g0")    


r1 = x+y-2
r2 = x+5*y-5
r3 = -x
r4 = -y
R = [r1,r2,r3,r4]

@dataclass
class CheckRes:
    limits : list = field(default_factory=list)
    fails: list = field(default_factory=list)

def check_restraint(x_v,y_v):
    print(f"  on ({x_v}, {y_v})")
    res = CheckRes()
    for (i,r) in enumerate(R):
        v= r.subs({x:x_v,y:y_v})
        if v > 0:
            print(f"  fail: r{i+1}")
            res.fails.append(i)
        elif v == 0:
            print(f"  limit: r{i+1}")
            res.limits.append(i)
        else:
            print(f"  pass: r{i+1}")
    return res

lims = check_restraint(0,0)
print(lims)

phi1 = expr.subs({x:2*d,y:3*d})
show(phi1,"phi1")
R1 = [r.subs({x:2*d,y:3*d})for r in R]
print(latex(Matrix(R1)))
lims = check_restraint(10/17,15/17)
g1 = grad.subs({x:Rational(10,17),y:Rational(15,17)})
show(g1,"g1")

sub_2 = {x:Rational(10,17) + 5*d,y:Rational(15,17)-d}
phi2 = expr.subs(sub_2)
show(phi2,"phi2")

x3 = Rational(10,17)+ 5*d
y3 = Rational(15,17)-d
x3_1 = x3.subs(d,Rational(57,527))
y3_1 = y3.subs(d,Rational(57,527))
check_restraint(x3.subs(d,Rational(57,527)),y3.subs(d,Rational(57,527)))
show(x3_1,"x_3")
show(y3_1,"y_3")
show(grad.subs({x:x3_1,y:y3_1}),"grad")

