function q=isotherm_CO2_silicasol_honeycomb(C,T) % C in [mol/m3], T in [K]
T0=273.15;
R=8.3145;
p=C*R*T; %[Pa]

b01=118.0789;
b02=0.0084;
ns01_1=0.2398;
ns01_2=1.5546;
ns02_1=5.3635;
ns02_2=25.1742;
xeta1_1=4.1941;
xeta1_2=3.5142;
xeta2_1=2.92;
xeta2_2=5.2444;
t01=0.3826;
t02=0.099;
alpha1=9.4569;
alpha2=4.7045;
delta_HToth1=124940;
delta_HToth2=120000;

if T<348
    ns1=ns01_1*exp(xeta1_1*(1-T0/T));
else
    ns1=ns01_2*exp(xeta1_2*(1-T/T0));
end
b1=b01*exp(delta_HToth1/R/T0*(T0/T-1));
t1=t01*exp(alpha1*(T/T0-1));

if T<328
    ns2=ns02_1*exp(xeta2_1*(1-T0/T));
else
    ns2=ns02_2*exp(xeta2_2*(1-T/T0));
end
b2=b02*exp(delta_HToth2/R/T0*(T0/T-1));
t2=t02*exp(alpha2*(T/T0-1));

q=ns1*b1*p/(1+(b1*p)^t1)^(1/t1)+ns2*b2*p/(2+(b2*p)^t2)^(2/t2); %[mmol/g]
