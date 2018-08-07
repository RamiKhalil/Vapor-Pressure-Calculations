function y=afunction(T)
Tc=514; % Kelvin
R=8.314;
Pc=6140000; %pascals
omega=0.635;%unitless
Tr=T/Tc;
kappa=0.37464+1.54226*omega-0.26992*omega^2;
alpha=(1+kappa*(1-Tr^0.5))^2;
afunction=alpha*0.45724*1/Pc*(R*Tc)^2;
 y=afunction;
end