% Program uporabljen pri seminarju 4, >>ZvezSist - 03<< , Višja dinamika  2015/16
clc; clear all; close all;

%% w(x,t) = W(x) * T(t)   Zapišimo funkcijo kraja W(x)
syms x beta L
C=sym('C',[4,1]);
W=C(1)*cos(beta*x)+C(2)*sin(beta*x)+C(3)*cosh(beta*x)+C(4)*sinh(beta*x);

% iz prvih dveh RP na desni strani dobimo C(1)=C(3) in C(2)=C(4)
W=subs(W,C(3),C(1));
W=subs(W,C(4),C(2));

%za robna pogoja na levi strani pri masni tocki rabimo drugi in tretji odvod
ddW=diff(W,x,x);
dddW=diff(ddW,x);
%za x vstavimo L
ddW_L=simplify(subs(ddW,x,L));
dddW_L=simplify(subs(dddW,x,L));
W_L=simplify(subs(W,x,L));

%zapisemo robne pogoje, koeficienti so še v enacbi
syms EI mk c
RP=[ddW_L;  EI * dddW_L+mk*W_L *(c^2) * beta^4 ];

for i=1:2
    [clen,koeficient]=coeffs(RP(i),C);   %razclenili bomo koeficiente
    A(i,:)=clen;
end
detA=simplify(det(A));

%% Zapisemo podatke
mk=4.56;
L=0.88;
a=0.041;
b=0.033;
EI= 1e10 * pi/4 * (a/2)^3 *b/2;
rho=900;
Area=pi* a/2 * b/2;
c=sqrt(EI/(rho*Area));
% podatke vnesemo v funkcijo
funkcija=eval(detA);

%Numericno resevanje
options = optimset('TolFun',1E-12,'Display','off');
eval(['funkcija= @(beta)' char(funkcija)]);

for i=0:20
    [x,fval]=fsolve(funkcija,i,options);
    bete(1+i,1)=x;
end
f = bete.^2*c/2/pi;

%% Spreminjamo maso macole, vpliv na prvo lastno frekvenco
m_m=(0:0.1:10);
for j=1:length(m_m)
    mk=m_m(j);
    funkcija=detA;
    funkcija=eval(detA);
    eval(['funkcija= @(beta)' char(funkcija)]);
    [x,fval]=fsolve(funkcija,6,options);
    bete_m(j,1)=x;
end
f_m = bete_m.^2*c/2/pi;
mk=4.56;
%% Izris in izpis
figure(1)
plot(m_m,f_m,'LineWidth',2)
hold on
plot(4.56,f(6),'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','r')
grid on
xlabel('Masa jeklenega dela macole [kg]'); ylabel('Prva lastna frekvenca  f _{01} [Hz]')
clc
% Najdemo prve 4 netrivialne frekvence
fprintf('w0_1 = %1.2f Hz \n',f(6))
fprintf('w0_2 = %1.2f Hz \n',f(9))
fprintf('w0_3 = %1.2f Hz \n',f(12)) 
fprintf('w0_4 = %1.2f Hz \n',f(16)) 

% %% Spremenimo se dolzino
% L_L=(0:0.1:10);
% for j=1:length(m_m)
%     L=L_L(j);
%     funkcija=detA;
%     funkcija=eval(detA);
%     eval(['funkcija= @(beta)' char(funkcija)]);
%     [x,fval]=fsolve(funkcija,6,options);
%     bete_m(j,1)=x;
% end
% f_L = bete_m.^2*c/2/pi;
% L=0.88;
% %% Izris za dolzino
% figure(2)
% plot(L_L,f_L,'LineWidth',2)
% hold on
% plot(0.88,f(6),'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','r')
% grid on
% xlabel('Dolžina macole [m]'); ylabel('Prva lastna frekvenca  \omega_{01} [Hz]')