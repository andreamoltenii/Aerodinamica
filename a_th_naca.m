% calcola angolo di theodorsen analiticamente per profili NACA

NACA = '2414';

m = str2num(NACA(1))/100;
p = str2num(NACA(2))/10;

phi = acos(1-2*p);

alpha_th = 2*m/pi/(1-p)^2*((p-.5)*(pi-phi)-.5*sin(phi));
if p ~= 0
    alpha_th = alpha_th + 2*m/pi/p^2*(p*phi-(phi-sin(phi))/2);
end
alpha_th = 180*alpha_th/pi

