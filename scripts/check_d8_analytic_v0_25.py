#!/usr/bin/env python3
import math
import sympy as sp

y,c=sp.symbols('y c', real=True)
Iq=sp.integrate((2-y)**2*(1-y)**2,(y,0,1))
Iqb=sp.integrate((2-y)**2,(y,0,1))
assert sp.simplify(Iq-sp.Rational(31,30))==0
assert sp.simplify(Iqb-sp.Rational(7,3))==0

# Fixed-helicity q-channel angular polynomial after undoing the
# unpolarized initial-quark average in the DIS cross section.
Fq=(3+c)*(1+c)/4
d1=(1+c)/2
d2=(1+c)*(2*c-1)/2
assert sp.simplify(Fq-(sp.Rational(7,4)*d1+sp.Rational(1,4)*d2))==0

# Strict |Re a_J| <= 1/2 limits at shat=Lambda^2.
# q channel: a1=7 C/(96 pi), so |C|<=48 pi/7.
# anti-q channel: M=C(3+c)=3 C P0 + C P1,
#                 a0=3 C/(16 pi), so |C|<=8 pi/3.
Cmax_q=48*math.pi/7
Cmax_qbar=8*math.pi/3
print(f'Iq={Iq} Iqbar={Iqb} ratio={float(Iqb/Iq):.9f}')
print(f'Cmax_q_J1_strict={Cmax_q:.9f}')
print(f'Cmax_qbar_J0_strict={Cmax_qbar:.9f}')
print('binding=anti-q J=0')
print('positivity_single_coefficient_branch=C<=0')
