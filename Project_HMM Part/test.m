clear all;

STATES_NUM = 49;
Quantify = 100;

E_G = zeros(STATES_NUM,Quantify);
for i=1:STATES_NUM
    E_G(i,i*2-1) = 0.1;
    E_G(i,i*2) = 0.2;
    E_G(i,i*2+1) = 0.4;
    E_G(i,i*2+2) = 0.2;
    E_G(i,i*2+3) = 0.1;
end
E_G = E_G(:,1:100);