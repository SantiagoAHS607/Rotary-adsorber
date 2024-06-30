function q_H2O=isotherm_H2O_silicasol_honeycomb(C,T)
R=8.3145; %[J/(mol K)]
p=C*R*T;
psat=10^6*exp(9.3876-3826.36/(T-45.47));
x=p/psat;

cm=1.1829;
cG=57.0632;
Kads=0.8818;

q_H2O=cm*cG*Kads*x/(1-Kads*x)/(1+(cG-1)*Kads*x);