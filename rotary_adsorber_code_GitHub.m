tic
Tatm=273.15+35; %[K]
Catm_CO2=0.0192; %[mol/m3]
Catm_H2O=0.6742; %ambient H2O concentration [mol/m3]
Tr_in=273.15+75; %inlet air temperature in regeneration section [K]
Cr_in_CO2=Catm_CO2; %inlet CO2 molar concentration in regeneration section [mol/m3]
Cr_in_H2O=Catm_H2O; %inlet H2O molar concentration in regeneration section [mol/m3]
Tp_in=Tatm; %inlet air temperature in process section [K]
Cp_in_CO2=Catm_CO2; %inlet air H2O concentration in process section [K]
Cp_in_H2O=Catm_H2O;
vr=1; %inlet air velocity in regeneration section [m/s]
vp=3; %inlet air velocity in adsorption section [m/s]

L=0.3; %wheel thickness [m]
Dout=0.4; %wheel outer diameter [m]
Din=0.02; %wheel inner diameter [m]
yita=0.8; %effective area coefficient
yita_e=0.9; %electric heater efficiency [-]
m=20;
dz=L/m;
ncycle=10000;
nr=ncycle/4;
np=ncycle-nr;
S_wheel=yita*pi*((Dout/2)^2-(Din/2)^2); %wheel cross-section area
Sr=S_wheel*nr/ncycle;
Sp=S_wheel*np/ncycle;
Fr=vr*Sr*3600; %regeneration air flowrate [m3/h]
Fp=vp*Sp*3600; %process air flowrate [m3/h]

r_speed=7; %wheel rotation speed
tcycle=1/r_speed*3600; %cycle duration [s]
dt=tcycle/ncycle;

T_rot=150; %motor torque [N m]
i_rot=20; %transmission ratio [-]
n_rot=(60/tcycle)*i_rot; %motor angular velocity [r/min]
P_motor=T_rot*2*pi*n_rot/60; %motor power [W]

lad=0.001896;%adsorbent linear density %[kg/m]
cp_m=0;
po=0.0325; %porosity [-]
a=0.001; %channel half height [m]
b=0.0015; %channel half width [m]

A=2*a*b; %channel cross-section area [m2]
R_channel=sqrt(A/pi); %channel equivalent diameter [m]
P=2*b+2*sqrt(b^2+a^2*pi^2)*(3+4*b^2/a^2/pi^2)/(4+4*b^2/a^2/pi^2); %channel perimeter [m]
f=0.8;
Nu=1.1791*(1+2.7701*(a/b)-3.1901*(a/b)^2+1.9975*(a/b)^3-0.4966*(a/b)^4);

vis_a=18.37e-6; %[Pa s]
den_a=1.168; %[kg/m3]
th_a=0.0321; %[W/(m K)]
cp_a=1004; %[J/(kg K)]
cp_c=846; %[J/(kg K)]
cp_v=1864; %[J/(kg K)]
th_ad=1.25; %[W/(m K)]
cp_ad=960; %[J/(kg K)]
Mc=0.044; %[kg/mol]
Mh=0.018; %[kg/mol]

DM_CO2=1.6e-5;
DK_CO2=1.2622e-6;
DS_CO2=4.6804e-11;
DP_CO2=1/(1/DM_CO2+1/DK_CO2);
DM_H2O=1.27e-4;
DK_H2O=1.9734e-6;
DS_H2O=3.1072e-11;
DP_H2O=1/(1/DM_H2O+1/DK_H2O);
rkd=1.4;
h=rkd*Nu*th_a*P/4/A;
Sh=Nu;
hm_CO2=Sh*DM_CO2*P/4/A;
hm_H2O=Sh*DM_H2O*P/4/A;
kA_CO2=0.0011; %[s-1]
kD_CO2=0.0123; %[s-1]
kA_H2O=0.005; %[s-1]
kD_H2O=0.01; %[s-1]

H_CO2=zeros(20, 1); %CO2 adosrption heats [J/mol]
H_H2O=zeros(20, 1); %H2O adosrption heats [J/mol]
qe_CO2=zeros(20, 1); %equilibrium CO2 adsorption capacity [mol/kg]
qe_H2O=zeros(20, 1); %equilibrium H2O adsorption capacity [mol/kg]
a1=zeros(160,160); %coefficient matrix of equations [-]
a2=zeros(160,1); %constant at the right side of equations [-]
a3=zeros(160,1); %output at previous time [-]
x=zeros(160,1); %for data exchanging [-]

w1_CO2=hm_CO2*P/f/A;
w1_H2O=hm_H2O*P/f/A;
w2=lad/po/(1-f)/A;
w3_CO2=hm_CO2*P/po/(1-f)/A;
w3_H2O=hm_H2O*P/po/(1-f)/A;
w7=(1-f)*(1-po)*A*th_ad;

w4=zeros(20, 1);
w5=zeros(20, 1);
w6=zeros(20, 1);

fid=fopen('adsorber_wheel_outlet_profile.csv','a+');
fid_channel=fopen('adsorber_channel_profile.csv','a+');

fprintf(fid,'Desorption gas inlet temperature[oC],'); fprintf(fid,'%f',Tr_in-273.15); fprintf(fid,'\r\n');
fprintf(fid,'Desorption gas inlet CO2 concentration[ppm],'); fprintf(fid,'%f',Cr_in_CO2*22.4*1000); fprintf(fid,'\r\n');
fprintf(fid,'Desorption gas inlet H2O concentration[o/o],'); fprintf(fid,'%f',Cr_in_H2O*22.4/10); fprintf(fid,'\r\n');
fprintf(fid,'Desorption gas velocity[m/s],'); fprintf(fid,'%f',vr); fprintf(fid,'\r\n');
fprintf(fid,'Desorption gas flowrate[m3/h],'); fprintf(fid,'%f',Fr); fprintf(fid,'\r\n');
fprintf(fid,'Process gas inlet temperature[oC],'); fprintf(fid,'%f',Tp_in-273.15); fprintf(fid,'\r\n');
fprintf(fid,'Process gas inlet CO2 concentration[ppm],'); fprintf(fid,'%f',Cp_in_CO2*22.4*1000); fprintf(fid,'\r\n');
fprintf(fid,'Process gas inlet H2O concentration[o/o],'); fprintf(fid,'%f',Cp_in_H2O*22.4/10); fprintf(fid,'\r\n');
fprintf(fid,'Process gas velocity[m/s],'); fprintf(fid,'%f',vp); fprintf(fid,'\r\n');
fprintf(fid,'Process gas flowrate[m3/h],'); fprintf(fid,'%f',Fp); fprintf(fid,'\r\n');
fprintf(fid,'CO2 Adsorption kinetic coefficient[1/s],'); fprintf(fid,'%f',kA_CO2); fprintf(fid,'\r\n');
fprintf(fid,'CO2 Desorption kinetic coefficient[1/s],'); fprintf(fid,'%f',kD_CO2); fprintf(fid,'\r\n');
fprintf(fid,'H2O Adsorption kinetic coefficient[1/s],'); fprintf(fid,'%f',kA_H2O); fprintf(fid,'\r\n');
fprintf(fid,'H2O Desorption kinetic coefficient[1/s],'); fprintf(fid,'%f',kD_H2O); fprintf(fid,'\r\n');
fprintf(fid,'Nu[-],'); fprintf(fid,'%f',Nu); fprintf(fid,'\r\n');
fprintf(fid,'Convection heat transfer coefficient[W/m2/oC],'); fprintf(fid,'%f',h); fprintf(fid,'\r\n');
fprintf(fid,'CO2 mass transfer coefficient[-],'); fprintf(fid,'%f',hm_CO2); fprintf(fid,'\r\n');
fprintf(fid,'H2O mass transfer coefficient[-],'); fprintf(fid,'%f',hm_H2O); fprintf(fid,'\r\n');
fprintf(fid,'Adsorbent linear density[kg/m],'); fprintf(fid,'%f',lad); fprintf(fid,'\r\n');
fprintf(fid,'Wheel thickness[m],'); fprintf(fid,'%f',L); fprintf(fid,'\r\n');
fprintf(fid,'Wheel outer diameter[m],'); fprintf(fid,'%f',Dout); fprintf(fid,'\r\n');
fprintf(fid,'Wheel inner diameter[m],'); fprintf(fid,'%f',Din); fprintf(fid,'\r\n');
fprintf(fid,'Channel half width[m],'); fprintf(fid,'%f',a); fprintf(fid,'\r\n');
fprintf(fid,'Channel half height[m],'); fprintf(fid,'%f',b); fprintf(fid,'\r\n');
fprintf(fid,'Rotation speed[r/h],'); fprintf(fid,'%f',3600/tcycle); fprintf(fid,'\r\n');

%initial conditions
q0_CO2=isotherm_CO2_silicasol_honeycomb(Catm_CO2, Tatm);
q0_H2O=isotherm_H2O_silicasol_honeycomb(Catm_H2O, Tatm);

for i=1:160
    if(i<=20)
        a3(i)=Catm_CO2;
    elseif(i>20&&i<=40)
        a3(i)=Catm_H2O;
    elseif(i>40&&i<=60)
        a3(i)=Catm_CO2;
    elseif(i>60&&i<=80)
        a3(i)=Catm_H2O;
    elseif(i>80&&i<=100)
        a3(i)=Tatm;
    elseif(i>100&&i<=120)
        a3(i)=Tatm;
    elseif(i>120&&i<=140)
        a3(i)=q0_CO2;
    else
        a3(i)=q0_H2O;
    end
end

for s=1:6 %s is the number of revolution of the wheel
    Cr_ave_CO2=0;
    Cr_ave_H2O=0;
    Tr_ave=0;
    qr_ave_CO2=0;
    qr_ave_H2O=0;
    Cp_ave_CO2=0;
    Cp_ave_H2O=0;
    Tp_ave=0;
    qp_ave_CO2=0;
    qp_ave_H2O=0;
    r_outlet_CO2_molar_amount=0;
    r_des_CO2_molar_amount=0;
    p_inlet_CO2_molar_amount=0;
    PR=0; 
    EC_thermal=0;
    EC_thermal_ideal=0;
    EC_blower=0; 

    %for regeneration section, control equations
    for k=1:nr
        for i=1:20
            qe_CO2(i)=isotherm_CO2_silicasol_honeycomb(a3(i+40),a3(i+100));
            qe_H2O(i)=isotherm_H2O_silicasol_honeycomb(a3(i+60),a3(i+100));
            H_CO2(i)=-295573*qe_CO2(i)^5+292760*qe_CO2(i)^4-111580*qe_CO2(i)^3+20260*qe_CO2(i)^2-1666.5*qe_CO2(i)-17.9296;
            H_H2O(i)=-43300;
            w4(i)=h*P/((cp_c*Mc*a3(i)+cp_v*Mh*a3(i+20)+cp_a*den_a)*f*A);
            w5(i)=(cp_c*hm_CO2*P*Mc*(a3(i)-a3(i+40))+cp_v*hm_H2O*P*Mh*(a3(i+20)-a3(i+60)))/((cp_c*Mc*a3(i)+cp_v*Mh*a3(i+20)+cp_a*den_a)*f*A);
            w6(i)=cp_ad*lad+cp_c*lad*Mc*a3(i+120)+cp_v*lad*Mh*a3(i+140)+cp_c*po*(1-f)*Mc*a3(i+40)*A+cp_v*po*(1-f)*Mh*a3(i+60)*A+cp_a*po*(1-f)*den_a*A;
        end

        %input parameters on the left side of equations
        for i=1:160
            if(i<=20)
                if(i==1)
                    a1(i,i)=1/dt+vr/dz+w1_CO2;
                    a1(i,i+40)=-w1_CO2;
                else
                    a1(i,i-1)=-vr/dz;
                    a1(i,i)=1/dt+vr/dz+w1_CO2;
                    a1(i,i+40)=-w1_CO2;
                end
            elseif(i>20&&i<=40)
                if(i==21)
                    a1(i,i)=1/dt+vr/dz+w1_H2O;
                    a1(i,i+40)=-w1_H2O;
                else
                    a1(i,i-1)=-vr/dz;
                    a1(i,i)=1/dt+vr/dz+w1_H2O;
                    a1(i,i+40)=-w1_H2O;
                end
            elseif(i>40&&i<=60)
                if(i==41)
                    a1(i,i-40)=-w3_CO2;
                    a1(i,i)=1/dt+DP_CO2/dz^2+w3_CO2;
                    a1(i,i+1)=-DP_CO2/dz^2;
                    a1(i,i+80)=w2/dt+w2*DS_CO2/dz^2;
                    a1(i,i+80+1)=-w2*DS_CO2/dz^2;
                elseif(i==60)
                    a1(i,i-40)=-w3_CO2;
                    a1(i,i-1)=-DP_CO2/dz^2;
                    a1(i,i)=1/dt+DP_CO2/dz^2+w3_CO2;
                    a1(i,i+80-1)=-w2*DS_CO2/dz^2;
                    a1(i,i+80)=w2/dt+w2*DS_CO2/dz^2;
                else
                    a1(i,i-40)=-w3_CO2;
                    a1(i,i-1)=-DP_CO2/dz^2;
                    a1(i,i)=1/dt+2*DP_CO2/dz^2+w3_CO2;
                    a1(i,i+1)=-DP_CO2/dz^2;
                    a1(i,i+80-1)=-w2*DS_CO2/dz^2;
                    a1(i,i+80)=w2/dt+2*w2*DS_CO2/dz^2;
                    a1(i,i+80+1)=-w2*DS_CO2/dz^2;
                end
            elseif(i>60&&i<=80)
                if(i==61)
                    a1(i,i-40)=-w3_H2O;
                    a1(i,i)=1/dt+DP_H2O/dz^2+w3_H2O;
                    a1(i,i+1)=-DP_H2O/dz^2;
                    a1(i,i+80)=w2/dt+w2*DS_H2O/dz^2;
                    a1(i,i+80+1)=-w2*DS_H2O/dz^2;
                elseif(i==80)
                    a1(i,i-40)=-w3_H2O;
                    a1(i,i-1)=-DP_H2O/dz^2;
                    a1(i,i)=1/dt+DP_H2O/dz^2+w3_H2O;
                    a1(i,i+80-1)=-w2*DS_H2O/dz^2;
                    a1(i,i+80)=w2/dt+w2*DS_H2O/dz^2;
                else
                    a1(i,i-40)=-w3_H2O;
                    a1(i,i-1)=-DP_H2O/dz^2;
                    a1(i,i)=1/dt+2*DP_H2O/dz^2+w3_H2O;
                    a1(i,i+1)=-DP_H2O/dz^2;
                    a1(i,i+80-1)=-w2*DS_H2O/dz^2;
                    a1(i,i+80)=w2/dt+2*w2*DS_H2O/dz^2;
                    a1(i,i+80+1)=-w2*DS_H2O/dz^2;
                end
            elseif(i>80&&i<=100)
                if(i==81)
                    a1(i,i)=1/dt+vr/dz+w4(i-80)+w5(i-80);
                    a1(i,i+20)=-w4(i-80)-w5(i-80);
                else
                    a1(i,i-1)=-vr/dz;
                    a1(i,i)=1/dt+vr/dz+w4(i-80)+w5(i-80);
                    a1(i,i+20)=-w4(i-80)-w5(i-80);
                end
            elseif(i>100&&i<=120)
                if(i==101)
                    a1(i,i-20)=-h*P-cp_c*hm_CO2*P*Mc*(a3(i-100)-a3(i-60))-cp_v*hm_H2O*P*Mh*(a3(i-80)-a3(i-40));
                    a1(i,i)=w6(i-100)/dt+w7/dz^2+h*P+cp_c*hm_CO2*P*Mc*(a3(i-100)-a3(i-60))+cp_v*hm_H2O*P*Mh*(a3(i-80)-a3(i-40));
                    a1(i,i+1)=-w7/dz^2;
                    a1(i,i+20)=lad*H_CO2(i-100)/dt;
                    a1(i,i+40)=lad*H_H2O(i-100)/dt;
                elseif(i==120)
                    a1(i,i-20)=-h*P-cp_c*hm_CO2*P*Mc*(a3(i-100)-a3(i-60))-cp_v*hm_H2O*P*Mh*(a3(i-80)-a3(i-40));
                    a1(i,i-1)=-w7/dz^2;
                    a1(i,i)=w6(i-100)/dt+w7/dz^2+h*P+cp_c*hm_CO2*P*Mc*(a3(i-100)-a3(i-60))+cp_v*hm_H2O*P*Mh*(a3(i-80)-a3(i-40));
                    a1(i,i+20)=lad*H_CO2(i-100)/dt;
                    a1(i,i+40)=lad*H_H2O(i-100)/dt;
                else
                    a1(i,i-20)=-h*P-cp_c*hm_CO2*P*Mc*(a3(i-100)-a3(i-60))-cp_v*hm_H2O*P*Mh*(a3(i-80)-a3(i-40));
                    a1(i,i-1)=-w7/dz^2;
                    a1(i,i)=w6(i-100)/dt+2*w7/dz^2+h*P+cp_c*hm_CO2*P*Mc*(a3(i-100)-a3(i-60))+cp_v*hm_H2O*P*Mh*(a3(i-80)-a3(i-40));
                    a1(i,i+1)=-w7/dz^2;
                    a1(i,i+20)=lad*H_CO2(i-100)/dt;
                    a1(i,i+40)=lad*H_H2O(i-100)/dt;
                end
            elseif(i>120&&i<=140)
                a1(i,i)=1/dt;
            elseif(i>140&&i<=160)
                a1(i,i)=1/dt;
            end
        end

        %input parameters on the right side of equations
        for i=1:160
            if(i<=20)
                if(i==1)
                    a2(i)=1/dt*a3(i)+vr/dz*Cr_in_CO2;  
                else
                    a2(i)=1/dt*a3(i);
                end
            elseif(i>20&&i<=40)
                if(i==21)
                    a2(i)=1/dt*a3(i)+vr/dz*Cr_in_H2O;
                else
                    a2(i)=1/dt*a3(i);
                end
            elseif(i>40&&i<=60)
                a2(i)=1/dt*a3(i)+w2/dt*a3(i+80);
            elseif(i>60&&i<=80)
                a2(i)=1/dt*a3(i)+w2/dt*a3(i+80);
            elseif(i>80&&i<=100)
                if(i==81)
                    a2(i)=1/dt*a3(i)+vr/dz*Tr_in;
                else
                    a2(i)=1/dt*a3(i);
                end
            elseif(i>100&&i<=120)
                a2(i)=w6(i-100)/dt*a3(i)+lad*H_CO2(i-100)/dt*a3(i+20)+lad*H_H2O(i-100)/dt*a3(i+40);
            elseif(i>120&&i<=140)
                if(a3(i)<qe_CO2(i-120))
                    a2(i)=kA_CO2*(qe_CO2(i-120)-a3(i))+1/dt*a3(i);
                else
                    a2(i)=-kD_CO2*(a3(i)-qe_CO2(i-120))+1/dt*a3(i);
                end
            elseif(i>140&&i<=160)
                if(a3(i)<qe_H2O(i-140))
                    a2(i)=kA_H2O*(qe_H2O(i-140)-a3(i))+1/dt*a3(i);
                else
                    a2(i)=-kD_H2O*(a3(i)-qe_H2O(i-140))+1/dt*a3(i);
                end
            end
        end

        a3=a1\a2;

        for i=1:20
            if a3(i)<0
                a3(i)=0;%error("CO2 concentration on the gas side of the regeneration section<0");
            end
            if a3(i+20)<0
                a3(i+20)=0;%error("H2O concentration on the gas side of the regeneration section<0");
            end
            if a3(i+40)<0
                a3(i+40)=0;%error("CO2 concentration on the solid side of the regeneration section<0");
            end
            if a3(i+60)<0
                a3(i+60)=0;%error("H2O concentration on the solid side of the regeneration section<0");
            end
            if a3(i+80)<0
                a3(i+80)=0;%error("temperature on the gas side of the regeneration section<0");
            end
            if a3(i+100)<0
                a3(i+100)=0;%error("temperature on the solid side of the regeneration section<0");
            end
        end

        %calculation of evaluation indicators
        qr_channel_ave_CO2=0;
        qr_channel_ave_H2O=0;
        for i=1:m
            qr_channel_ave_CO2=qr_channel_ave_CO2+a3(i+120); 
            qr_channel_ave_H2O=qr_channel_ave_H2O+a3(i+140);
        end
        qr_channel_ave_CO2=qr_channel_ave_CO2/m; 
        qr_channel_ave_H2O=qr_channel_ave_H2O/m; 

        Cr_ave_CO2=Cr_ave_CO2+a3(20); 
        Cr_ave_H2O=Cr_ave_H2O+a3(40); 
        Tr_ave=Tr_ave+a3(100); 
        qr_ave_CO2=qr_ave_CO2+qr_channel_ave_CO2; 
        qr_ave_H2O=qr_ave_H2O+qr_channel_ave_H2O; 
        r_outlet_CO2_molar_amount=r_outlet_CO2_molar_amount+a3(20)*vr*Sr; 
        r_des_CO2_molar_amount=r_des_CO2_molar_amount+(a3(20)-Cr_in_CO2)*vr*Sr; 
        PR=PR+((a3(20)-Cr_in_CO2)*vr*f*A)/(lad*L); 
        EC_thermal=EC_thermal+den_a*vr*Sr*cp_a*(Tr_in-Tatm); 
        EC_thermal_ideal=EC_thermal_ideal+den_a*vr*Sr*cp_a*(Tr_in-a3(100)); 
        EC_blower=EC_blower+8*L*vis_a/R_channel^2*vr*Sr; 

        if k==1
            fprintf(fid,'No.');fprintf(fid,'%f,',s);
            fprintf(fid,'rotation degree,');
            fprintf(fid,'C_CO2_a(ppm),');
            fprintf(fid,'C_H2O_a(percent),');
            fprintf(fid,'C_CO2_s(ppm),');
            fprintf(fid,'C_H2O_s(percent),');
            fprintf(fid,'T_a(oC),');
            fprintf(fid,'T_s(oC),');
            fprintf(fid,'q_ave_CO2(mmol/g),');
            fprintf(fid,'q_ave_H2O(mmol/g)');
            fprintf(fid,'\r\n');
            fprintf(fid_channel,'No.');fprintf(fid_channel,'%f',s);
            fprintf(fid_channel,'\r\n');
        end

        if(mod(k,100)==0) %output wheel outlet parameter profiles to Excel
            if k==100
                fprintf(fid,'Enter regeneration section,');
                wr=[3.6*k/100 a3(20)*22.4*1000 a3(40)*22.4/10 a3(60)*22.4*1000 a3(80)*22.4/10 a3(100)-273.15 a3(120)-273.15 qr_channel_ave_CO2 qr_channel_ave_H2O];
                fprintf(fid,'%f, %f, %f, %f, %f, %f, %f, %f, %f',wr);
                fprintf(fid,'\r\n');
            else
                wr=[3.6*k/100 a3(20)*22.4*1000 a3(40)*22.4/10 a3(60)*22.4*1000 a3(80)*22.4/10 a3(100)-273.15 a3(120)-273.15 qr_channel_ave_CO2 qr_channel_ave_H2O];
                fprintf(fid,',%f, %f, %f, %f, %f, %f, %f, %f, %f',wr);
                fprintf(fid,'\r\n');
            end
        end

        if (mod(k,100)==0) %output channel parameter profiles to Excel
            if k==100
                fprintf(fid_channel, 'Enter regeneration section');fprintf(fid_channel, '\r\n');
                fprintf(fid_channel,'Gas side CO2 concentration,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'Gas side H2O concentration,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'Solid side CO2 concentration,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'Solid side H2O concentration,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'Gas side temperature,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'Solid side temperature,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'CO2 uptake,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'H2O uptake');
                fprintf(fid_channel,'\r\n');
            end
            for i=1:160
                if (i<=20)||(i>40&&i<=60)
                    fprintf(fid_channel, '%f,', a3(i)*22.4*1000);
                elseif (i>20&&i<=40)||(i>60&&i<=80)
                    fprintf(fid_channel, '%f,', a3(i)*22.4/10);
                elseif i>80&&i<=120
                    fprintf(fid_channel, '%f,', a3(i)-273.15);
                else
                    fprintf(fid_channel, '%f,', a3(i));
                end
                if (mod(i,20)==0)
                    fprintf(fid_channel,',');
                end
            end
            fprintf(fid_channel,'\r\n');
        end

        for i=1:160
            a2(i)=0;
            for j=1:160
                a1(i,j)=0;
            end
        end

    end %end current time step

    for i=1:160 %reverse the data, as the regeneration and adsorption gases are arranged in a countercurrent manner
        x(i)=a3(i);
    end
    for i=1:20
        a3(i)=x(21-i);
        a3(i+20)=x(41-i);
        a3(i+40)=x(61-i);
        a3(i+60)=x(81-i);
        a3(i+80)=x(101-i);
        a3(i+100)=x(121-i);
        a3(i+120)=x(141-i);
        a3(i+140)=x(161-i);
    end

    %for adsorption section，control equations
    for k=1:np
        for i=1:20
            qe_CO2(i)=isotherm_CO2_silicasol_honeycomb(a3(i+40),a3(i+100));
            qe_H2O(i)=isotherm_H2O_silicasol_honeycomb(a3(i+60),a3(i+100));
            H_CO2(i)=-295570*qe_CO2(i)^5+292760*qe_CO2(i)^4-111580*qe_CO2(i)^3+20260*qe_CO2(i)^2-1666.5*qe_CO2(i)-17.9296;
            H_H2O(i)=-43300;
            w4(i)=h*P/((cp_c*Mc*a3(i)+cp_v*Mh*a3(i+20)+cp_a*den_a)*f*A);
            w5(i)=(cp_c*hm_CO2*P*Mc*(a3(i)-a3(i+40))+cp_v*hm_H2O*P*Mh*(a3(i+20)-a3(i+60)))/((cp_c*Mc*a3(i)+cp_v*Mh*a3(i+20)+cp_a*den_a)*f*A);
            w6(i)=cp_ad*lad+cp_c*lad*Mc*a3(i+120)+cp_v*lad*Mh*a3(i+140)+cp_c*po*(1-f)*Mc*a3(i+40)*A+cp_v*po*(1-f)*Mh*a3(i+60)*A+cp_a*po*(1-f)*den_a*A;
        end

        %input parameters on the left side of equations
        for i=1:160
            if(i<=20)
                if(i==1)
                    a1(i,i)=1/dt+vp/dz+w1_CO2;
                    a1(i,i+40)=-w1_CO2;
                else
                    a1(i,i-1)=-vp/dz;
                    a1(i,i)=1/dt+vp/dz+w1_CO2;
                    a1(i,i+40)=-w1_CO2;
                end
            elseif(i>20&&i<=40)
                if(i==21)
                    a1(i,i)=1/dt+vp/dz+w1_H2O;
                    a1(i,i+40)=-w1_H2O;
                else
                    a1(i,i-1)=-vp/dz;
                    a1(i,i)=1/dt+vp/dz+w1_H2O;
                    a1(i,i+40)=-w1_H2O;
                end
            elseif(i>40&&i<=60)
                if(i==41)
                    a1(i,i-40)=-w3_CO2;
                    a1(i,i)=1/dt+DP_CO2/dz^2+w3_CO2;
                    a1(i,i+1)=-DP_CO2/dz^2;
                    a1(i,i+80)=w2/dt+w2*DS_CO2/dz^2;
                    a1(i,i+80+1)=-w2*DS_CO2/dz^2;
                elseif(i==60)
                    a1(i,i-40)=-w3_CO2;
                    a1(i,i-1)=-DP_CO2/dz^2;
                    a1(i,i)=1/dt+DP_CO2/dz^2+w3_CO2;
                    a1(i,i+80-1)=-w2*DS_CO2/dz^2;
                    a1(i,i+80)=w2/dt+w2*DS_CO2/dz^2;
                else
                    a1(i,i-40)=-w3_CO2;
                    a1(i,i-1)=-DP_CO2/dz^2;
                    a1(i,i)=1/dt+2*DP_CO2/dz^2+w3_CO2;
                    a1(i,i+1)=-DP_CO2/dz^2;
                    a1(i,i+80-1)=-w2*DS_CO2/dz^2;
                    a1(i,i+80)=w2/dt+2*w2*DS_CO2/dz^2;
                    a1(i,i+80+1)=-w2*DS_CO2/dz^2;
                end
            elseif(i>60&&i<=80)
                if(i==61)
                    a1(i,i-40)=-w3_H2O;
                    a1(i,i)=1/dt+DP_H2O/dz^2+w3_H2O;
                    a1(i,i+1)=-DP_H2O/dz^2;
                    a1(i,i+80)=w2/dt+w2*DS_H2O/dz^2;
                    a1(i,i+80+1)=-w2*DS_H2O/dz^2;
                elseif(i==80)
                    a1(i,i-40)=-w3_H2O;
                    a1(i,i-1)=-DP_H2O/dz^2;
                    a1(i,i)=1/dt+DP_H2O/dz^2+w3_H2O;
                    a1(i,i+80-1)=-w2*DS_H2O/dz^2;
                    a1(i,i+80)=w2/dt+w2*DS_H2O/dz^2;
                else
                    a1(i,i-40)=-w3_H2O;
                    a1(i,i-1)=-DP_H2O/dz^2;
                    a1(i,i)=1/dt+2*DP_H2O/dz^2+w3_H2O;
                    a1(i,i+1)=-DP_H2O/dz^2;
                    a1(i,i+80-1)=-w2*DS_H2O/dz^2;
                    a1(i,i+80)=w2/dt+2*w2*DS_H2O/dz^2;
                    a1(i,i+80+1)=-w2*DS_H2O/dz^2;
                end
            elseif(i>80&&i<=100)
                if(i==81)
                    a1(i,i)=1/dt+vp/dz+w4(i-80)+w5(i-80);
                    a1(i,i+20)=-w4(i-80)-w5(i-80);
                else
                    a1(i,i-1)=-vp/dz;
                    a1(i,i)=1/dt+vp/dz+w4(i-80)+w5(i-80);
                    a1(i,i+20)=-w4(i-80)-w5(i-80);
                end
            elseif(i>100&&i<=120)
                if(i==101)
                    a1(i,i-20)=-h*P-cp_c*hm_CO2*P*Mc*(a3(i-100)-a3(i-60))-cp_v*hm_H2O*P*Mh*(a3(i-80)-a3(i-40));
                    a1(i,i)=w6(i-100)/dt+w7/dz^2+h*P+cp_c*hm_CO2*P*Mc*(a3(i-100)-a3(i-60))+cp_v*hm_H2O*P*Mh*(a3(i-80)-a3(i-40));
                    a1(i,i+1)=-w7/dz^2;
                    a1(i,i+20)=lad*H_CO2(i-100)/dt;
                    a1(i,i+40)=lad*H_H2O(i-100)/dt;
                elseif(i==120)
                    a1(i,i-20)=-h*P-cp_c*hm_CO2*P*Mc*(a3(i-100)-a3(i-60))-cp_v*hm_H2O*P*Mh*(a3(i-80)-a3(i-40));
                    a1(i,i-1)=-w7/dz^2;
                    a1(i,i)=w6(i-100)/dt+w7/dz^2+h*P+cp_c*hm_CO2*P*Mc*(a3(i-100)-a3(i-60))+cp_v*hm_H2O*P*Mh*(a3(i-80)-a3(i-40));
                    a1(i,i+20)=lad*H_CO2(i-100)/dt;
                    a1(i,i+40)=lad*H_H2O(i-100)/dt;
                else
                    a1(i,i-20)=-h*P-cp_c*hm_CO2*P*Mc*(a3(i-100)-a3(i-60))-cp_v*hm_H2O*P*Mh*(a3(i-80)-a3(i-40));
                    a1(i,i-1)=-w7/dz^2;
                    a1(i,i)=w6(i-100)/dt+2*w7/dz^2+h*P+cp_c*hm_CO2*P*Mc*(a3(i-100)-a3(i-60))+cp_v*hm_H2O*P*Mh*(a3(i-80)-a3(i-40));
                    a1(i,i+1)=-w7/dz^2;
                    a1(i,i+20)=lad*H_CO2(i-100)/dt;
                    a1(i,i+40)=lad*H_H2O(i-100)/dt;
                end
            elseif(i>120&&i<=140)
                a1(i,i)=1/dt;
            elseif(i>140&&i<=160)
                a1(i,i)=1/dt;
            end
        end

        %input parameters on the right side of equations
        for i=1:160
            if(i<=20)
                if(i==1)
                    a2(i)=1/dt*a3(i)+vp/dz*Cp_in_CO2;
                else
                    a2(i)=1/dt*a3(i);
                end
            elseif(i>20&&i<=40)
                if(i==21)
                    a2(i)=1/dt*a3(i)+vp/dz*Cp_in_H2O;
                else
                    a2(i)=1/dt*a3(i);
                end
            elseif(i>40&&i<=60)
                a2(i)=1/dt*a3(i)+w2/dt*a3(i+80);
            elseif(i>60&&i<=80)
                a2(i)=1/dt*a3(i)+w2/dt*a3(i+80);
            elseif(i>80&&i<=100)
                if(i==81)
                    a2(i)=1/dt*a3(i)+vp/dz*Tp_in;
                else
                    a2(i)=1/dt*a3(i);
                end
            elseif(i>100&&i<=120)
                a2(i)=w6(i-100)/dt*a3(i)+lad*H_CO2(i-100)/dt*a3(i+20)+lad*H_H2O(i-100)/dt*a3(i+40);
            elseif(i>120&&i<=140)
                if(a3(i)<qe_CO2(i-120))
                    a2(i)=kA_CO2*(qe_CO2(i-120)-a3(i))+1/dt*a3(i);
                else
                    a2(i)=-kD_CO2*(a3(i)-qe_CO2(i-120))+1/dt*a3(i);
                end
            elseif(i>140&&i<=160)
                if(a3(i)<qe_H2O(i-140))
                    a2(i)=kA_H2O*(qe_H2O(i-140)-a3(i))+1/dt*a3(i);
                else
                    a2(i)=-kD_H2O*(a3(i)-qe_H2O(i-140))+1/dt*a3(i);
                end
            end
        end

        a3=a1\a2;

        for i=1:20
            if a3(i)<0
                a3(i)=0; %error("CO2 concentration on the gas side of the adsorption section<0");
            end
            if a3(i+20)<0
                a3(i+20)=0; % error("H2O concentration on the gas side of the adsorption section<0");
            end
            if a3(i+40)<0
                a3(i+40)=0; % error("CO2 concentration on the solid side of the adsorption section<0");
            end
            if a3(i+60)<0
                a3(i+60)=0; % error("H2O concentration on the solid side of the adsorption section<0");
            end
            if a3(i+80)<0
                a3(i+80)=0;%error("temperature on the gas side of the adsorption section<0");
            end
            if a3(i+100)<0
                a3(i+100)=0;%error("temperature on the solid side of the adsorption section<0");
            end
        end

        qp_channel_ave_CO2=0;
        qp_channel_ave_H2O=0;
        for i=1:m
            qp_channel_ave_CO2=qp_channel_ave_CO2+a3(i+120); 
            qp_channel_ave_H2O=qp_channel_ave_H2O+a3(i+140); 
        end
        qp_channel_ave_CO2=qp_channel_ave_CO2/m; 
        qp_channel_ave_H2O=qp_channel_ave_H2O/m; 
        
        Cp_ave_CO2=Cp_ave_CO2+a3(20); 
        Cp_ave_H2O=Cp_ave_H2O+a3(40); 
        Tp_ave=Tp_ave+a3(100); 
        qp_ave_CO2=qp_ave_CO2+qp_channel_ave_CO2; 
        qp_ave_H2O=qp_ave_H2O+qp_channel_ave_H2O; 
        p_inlet_CO2_molar_amount=p_inlet_CO2_molar_amount+Cp_in_CO2*vp*Sp; 
        EC_blower=EC_blower+8*L*vis_a/R_channel^2*vp*Sp; 

        if(mod(k,100)==0) %output wheel outlet parameter profiles to Excel
            if k==100
                wr=[3.6*(k+nr)/100 a3(20)*22.4*1000 a3(40)*22.4/10 a3(60)*22.4*1000 a3(80)*22.4/10 a3(100)-273.15 a3(120)-273.15 qp_channel_ave_CO2 qp_channel_ave_H2O];
                fprintf(fid,'Enter adsorption section,');
                fprintf(fid,'%f, %f, %f, %f, %f, %f, %f, %f, %f',wr);
                fprintf(fid,'\r\n');
            else
                wr=[3.6*(k+nr)/100 a3(20)*22.4*1000 a3(40)*22.4/10 a3(60)*22.4*1000 a3(80)*22.4/10 a3(100)-273.15 a3(120)-273.15 qp_channel_ave_CO2 qp_channel_ave_H2O];
                fprintf(fid,',%f, %f, %f, %f, %f, %f, %f, %f, %f',wr);
                fprintf(fid,'\r\n');
            end
        end

        if (mod(k,100)==0) %output channel parameter profiles to Excel
            if k==100
                fprintf(fid_channel, 'Enter adsorption section');fprintf(fid_channel, '\r\n');
                fprintf(fid_channel,'Gas side CO2 concentration,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'Gas side H2O concentration,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'Solid side CO2 concentration,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'Solid side H2O concentration,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'Gas side temperature,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'Solid side temperature,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'CO2 uptake,,,,,,,,,,,,,,,,,,,,,');
                fprintf(fid_channel,'H2O uptake');
                fprintf(fid_channel,'\r\n');
            end
            for i=1:160
                if (i<=20)||(i>40&&i<=60)
                    fprintf(fid_channel, '%f,', a3(i)*22.4*1000);
                elseif (i>20&&i<=40)||(i>60&&i<=80)
                    fprintf(fid_channel, '%f,', a3(i)*22.4/10);
                elseif i>80&&i<=120
                    fprintf(fid_channel, '%f,', a3(i)-273.15);
                else
                    fprintf(fid_channel, '%f,', a3(i));
                end
                if (mod(i,20)==0)
                    fprintf(fid_channel,',');
                end
            end
            fprintf(fid_channel,'\r\n');
        end

        for i=1:160
            a2(i)=0;
            for j=1:160
                a1(i,j)=0;
            end
        end
    end %end current time step

    for i=1:160 %reverse the data
        x(i)=a3(i);
    end
    for i=1:20
        a3(i)=x(21-i);
        a3(i+20)=x(41-i);
        a3(i+40)=x(61-i);
        a3(i+60)=x(81-i);
        a3(i+80)=x(101-i);
        a3(i+100)=x(121-i);
        a3(i+120)=x(141-i);
        a3(i+140)=x(161-i);
    end

    %calculation of evaluation indicators
    Cr_ave_CO2=Cr_ave_CO2/nr; %the average CO2 concentration at the outlet of the regeneration section
    Cr_ave_H2O=Cr_ave_H2O/nr; %the average H2O concentration at the outlet of the regeneration section
    Tr_ave=Tr_ave/nr; %the average temperature at the outlet of the regeneration section
    qr_ave_CO2=qr_ave_CO2/nr; %the average CO2 uptake of the regeneration section
    qr_ave_H2O=qr_ave_H2O/nr; %the average H2O uptake of the regeneration section
    Cp_ave_CO2=Cp_ave_CO2/np; %the average CO2 concentration at the outlet of the adsorption section
    Cp_ave_H2O=Cp_ave_H2O/np; %the average H2O concentration at the outlet of the adsorption section
    Tp_ave=Tp_ave/np; %the average temperature at the outlet of the adsorption section
    qp_ave_CO2=qp_ave_CO2/np; %the average CO2 uptake of the adsorption section
    qp_ave_H2O=qp_ave_H2O/np; %the average H2O uptake of the adsorption section
    CCR=1-Cp_ave_CO2/Cp_in_CO2; %carbon capture rate
    CP=Cr_ave_CO2; %product gas CO2 concentration
    PR=PR/(nr+np); %CO2 productivity
    EC_thermal=EC_thermal/r_des_CO2_molar_amount/1000/yita_e; %thermal energy consumption
    EC_thermal_ideal=EC_thermal_ideal/r_des_CO2_molar_amount/1000/yita_e; %ideal thermal energy consumption
    EC_blower=EC_blower/1000/r_outlet_CO2_molar_amount; %blower energy consumption
    EC_motor=P_motor*tcycle/1000/r_outlet_CO2_molar_amount; %motor energy consumption

    %output evaluation indicators
    fprintf(fid,'CCR(percent), CP(ppm), PR(mol/kg/day), EC_thermal(kJ/mol), EC_thermal_ideal(kJ/mol), EC_motor(kJ/mol), EC_blower(kJ/mol), qr_ave_CO2(mmol/g), qr_ave_H2O(mmol/g), qp_ave_CO2(mmol/g), qp_ave_H2O(mmol/g), r_des_CO2_molar_amount(mol)');fprintf(fid,'\r\n');
    fprintf(fid,'%f',CCR*100);fprintf(fid,',');
    fprintf(fid,'%f',CP*22.4*1000);fprintf(fid,',');
    fprintf(fid,'%f',PR*3600*24);fprintf(fid,',');
    fprintf(fid,'%f',EC_thermal);fprintf(fid,',');
    fprintf(fid,'%f',EC_thermal_ideal);fprintf(fid,',');
    fprintf(fid,'%f',EC_motor);fprintf(fid,',');
    fprintf(fid,'%f',EC_blower);fprintf(fid,',');
    fprintf(fid,'%f',qr_ave_CO2);fprintf(fid,',');
    fprintf(fid,'%f',qr_ave_H2O);fprintf(fid,',');
    fprintf(fid,'%f',qp_ave_CO2);fprintf(fid,',');
    fprintf(fid,'%f',qp_ave_H2O);fprintf(fid,',');
    fprintf(fid,'%f',r_des_CO2_molar_amount);fprintf(fid,'\r\n');
    fprintf(fid,'CO2_r_ave_out(ppm), H2O_r_ave_out(percent), T_r_ave_out(oC), CO2_p_ave_out(ppm), H2O_p_ave_out(percent), T_p_ave_out(oC)');fprintf(fid,'\r\n');
    fprintf(fid,'%f,',Cr_ave_CO2*22.4*1000); fprintf(fid,'%f,',Cr_ave_H2O*22.4/10); fprintf(fid,'%f,',Tr_ave-273.15);
    fprintf(fid,'%f,',Cp_ave_CO2*22.4*1000); fprintf(fid,'%f,',Cp_ave_H2O*22.4/10); fprintf(fid,'%f,',Tp_ave-273.15); fprintf(fid,'\r\n');
end %end current cycle
fclose(fid);
fclose(fid_channel);
toc