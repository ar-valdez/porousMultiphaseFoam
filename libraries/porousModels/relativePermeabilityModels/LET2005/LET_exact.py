 #!usr/bin/python
from sympy import *
x = Symbol("x")
y = Symbol("y")

Sw = Symbol("Sw")
Sw_max = Symbol("Sw_max")
Sw_min = Symbol("Sw_min")
Kw = Symbol("Kw")
Ko = Symbol("Ko")

Lw = Symbol("Lw")
Ew = Symbol("Ew")
Tw = Symbol("Tw")
Lo = Symbol("Lo")
Eo = Symbol("Eo")
To = Symbol("To")

Swe = (Sw - Sw_min) / (Sw_max - Sw_min)

#Swe = Symbol("Swe")
Krw = Kw * (Swe**Lw) / (Swe**Lw + Ew*(1-Swe)**Tw)
Kro = Ko * ((1-Swe)**Lo) / ((1-Swe)**Lo + Eo*(Swe)**To)

print('Water perm',Krw)
print('')
print('Oil perm',Kro)
print('')

print('Kw der',simplify(diff(Krw,Sw)))
print('')
print('Ko der',simplify(diff(Kro,Sw)))
print('')

