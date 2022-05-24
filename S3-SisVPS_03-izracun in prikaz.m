% Program uporabljen pri predmetu Višja dinamika 2015/16, Seminar 3, SisVPS-03
clc;clear all; close all;

%% Podatki:
N=6;                                        %stevilo razdelkov
m=223;                                   %skupna masa nosilca [kg]
L=7.2;                                     %dolžina nosilca [L]
EI=2.1e11 * 3.95e-5;            %modul * vstrajnostni moment  [Pa * m^4 = Nm^2]
A=0.0049;                              %prerez [m^2]
F0=311.8;                              %amplituda sile [N]
g=9.81;                                   %težnostni pospešek [m/s^2]

n=N-1;                                     %št. masnih tock
mi=m/n;                                   %masa tocke
ku=1.61*EI/(2*L/N)^3 * (L/N)^2;         %nadomestna togost

%% 2. del  -> Masna in togostna matrika
M=zeros(n);                             %najprej masna matrika sistema
for i=1:n
    M(i,i)=mi;                              %masna matrika je diagonalna
end

K=zeros(n);                              %togostna matrika
K(1,1:2)=[2, -1];                       %primer na zgornjem levem delu matrike
K(n,n-1:n)=[-1,2];                      %primer na spodnjem desnem delu matrike
for i=2:n-1                                 
    K(i,(i-1):(i+1))=[-1,2,-1];       %zapisemo vse primere, ki so vmes
end
K=ku*K;                                     %matriko pomnozimo se z nadomestno togostjo

%% 3. del  -> Dolocitev lastnih frekvenc, oblik nosilca
Amat=M\K;                                %uporabili bomo metodo eigenvalues, ki zahteva tak vhod
[oblike,frek]=eig(Amat);           %oblike predstavlja matriko lastnih vrednosti, frek pa matriko lastnih kotnih frekvenc
nicle=zeros(1,N-1);                      %zraven lastnih vrednosti dodamo se nicle, ki predstavljajo vpetje
oblike1=[nicle;oblike;nicle];
barve={'k','b','r','g','m'};
figure(1)
x=0:N;
for i=1:N-1
    plot(x,oblike1(:,i)/oblike1(2,i),barve{i},'LineWidth',3)      %tu normiramo na prvo masno tocko
    hold on; grid on;
    xlabel('Tocke na nosilcu (N+1)');
    ylabel('Amplituda normiranega odmika posamezne tocke')
end
legend('X_{1}','X_{2}','X_{3}','X_{4}','X_{5}')
frek=sqrt(frek)                           % Lastna kotna frekvenca [rad/s]
f=frek/(2*pi);                              % Lastne frekvence [Hz]

%% 4. del -> Prehod v modalne koordinate
Fi=zeros(N-1,N-1);
for i=1:N-1
    Fi(:,i)=oblike(:,i)/oblike(1,i);                %normiramo lastni vektor
end
Mmod=transpose(Fi) * M * Fi;                                  %Naša modalna masna matrika
Kmod=transpose(Fi) * K * Fi;                                   %Naša modalna togostna matrika

%% 5. del  -> Dolocitev frekvence vzbujanja
w_vzb=1.2*frek(3,3);
f_vzb = w_vzb/(2*pi);
T=1/f_vzb;

%% 6. del -> Doloèitev odziva sistema
     % 6.1)  Popis sile
%izracun koeficientov ki jih bomo uporabili pri vsoti (da ne traja predolgo)
% syms n t F0 T w_vzb
% a0=(2/T)*int(2*F0/T*t,t,0,T/2)+(2/T)*int(-F0,t,T/2,T);
% aj=simplify((2/T)*int(2*F0/T*t*cos(n*w_vzb*t),t,0,T/2)+(2/T)*int(-F0*cos(n*w_vzb*t),t,T/2,T));
% bj=simplify((2/T)*int(2*F0/T*t*sin(n*w_vzb*t),t,0,T/2)+(2/T)*int(-F0*sin(n*w_vzb*t),t,T/2,T));     
Cleni = 500;                              %popis sile s 500 cleni              
a0=-F0/2;
aj=zeros(1,Cleni);                     %priprava vektorjev za vsoto
bj=zeros(1,Cleni);
for n=1:Cleni
    aj(n)=-(2*F0*(4*sin((T*n*w_vzb)/4)^2 + T*n*w_vzb*sin(T*n*w_vzb) - 2*T*n*w_vzb*sin((T*n*w_vzb)/2)))/(T^2*n^2*w_vzb^2);
    bj(n)=(2*F0*(2*sin((T*n*w_vzb)/2) - T*n*w_vzb - 2*T*n*w_vzb*cos((T*n*w_vzb)/2) + 2*T*n*w_vzb*cos((T*n*w_vzb)/2)^2))/(T^2*n^2*w_vzb^2);
end
cas=0:T/5000:T;                     % x-os
F=a0/2;                                      %F je vzbujevalna sila enega motorja
for i=1:Cleni
    F=F+aj(i)*cos(i*w_vzb*cas)+bj(i)*sin(i*w_vzb*cas);
end
Ffilt=sgolayfilt(F,0,71);                       %uporabimo filter, da odstranimo nelinearnosti pri popisu sile
figure(2)
plot(cas,F,'LineWidth',2)
hold on
plot(cas,Ffilt,'r','LineWidth',2)
legend('Fouriejeva vrsta','Filter')
xlabel('Cas[s]');ylabel('F [N]')
grid on

sila=[1;0;1;0;1];
g=transpose(Fi) *sila;                                       %Vektor sile v modalnih koordinatah

    % 6.2) Izracun odziva
for i=1:(N-1)
    c(i,1)=(g(i,1)/Mmod(i,i))^-1;            %normirana sila v posamezni tocki
    w_i(i,1)=sqrt(Kmod(i,i)/Mmod(i,i));
end


delta=0;                                                %Nimamo dusenja. Razmernik dusenja=0
p=1;
%etac=zeros(N-1,5001);
for t=0:T/5000:T                                    %razdelimo na isto delov kokr smo silo
    for j=1:Cleni                                       %stevilo clenov fouerjeve vrste
        for i=[1 3 5]
            fi(i,j) = atan2(2*delta*j*w_vzb/w_i(i,1), 1-(j*w_vzb/w_i(i,1))^2);               %fazni zamik ni vedno niè
            beta(i,j) = 1/(1-(j*w_vzb / w_i(i,1))^2);
            eta1(i,j) = aj(j) / (c(i,1) * (w_i(i,1))^2) * beta(i,j) * cos(j*w_vzb*t-fi(i,j));
            eta2(i,j) = bj(j) / (c(i,1) * (w_i(i,1))^2) * beta(i,j) * sin(j*w_vzb*t-fi(i,j));
        end
    end
    for i=[1 3 5]
        etac(i,p)=a0/(2*(w_i(i,1))^2 * c(i,1)) + sum(eta1(i,:)) + sum(eta2(i,:));
    end
    p=p+1;
end
y= Fi * etac;                                                       %realne koordinate

%% 7. del -> Prikaz odziva sistema
tocka={'Odziv prve masne tocke','Odziv druge masne tocke','Odziv tretje masne tocke','Odziv cetrte masne tocke','Odziv pete masne tocke'};
for i=1:N-1
    figure(i+2)
    plot(cas,y(i,:)*1000,'LineWidth',2)                                       %izris odziva posamezne toèke
    grid on;
    xlabel('t [s]'); ylabel('Pomik [mm]'); title(tocka{i})
    axis([0 T -0.3 0.15]);
end
%% ANIMACIJA
steviloperiod = 1;                                                  %stevilo period
Naslov=sprintf('Animacija odziva sistema v ustaljenem stanju, stevilo period T= %d ',steviloperiod);
for i=1:steviloperiod                      
figure(8)
    for t=1:30:length(y)
        set(gcf,'color','w');
        pomik=[0 y(1,t)*1000 y(2,t)*1000 y(3,t)*1000 y(4,t)*1000 y(5,t)*1000 0];
        pozicija=[0 L/6 2*L/6 3*L/6 4*L/6 5*L/6 L];
        plot(pozicija,pomik,'LineWidth',2,'Marker','s')  
        %axis([0 L -0.3 0.15]);
        axis([0 L -1 1]);
        grid on;
        xlabel('x [m]'); ylabel('Pomik [mm]'); title(Naslov )
        drawnow
    end
end
%% 8. del -> Lastne frekvence pri veèji diskretizaciji
novefrek=zeros(1,5);
inkrement=1;
segmenti=[6 24 60 120 600];
for i=1:length(segmenti)
    razdelki=segmenti(i);
    tocke=razdelki-1;
    kun=1.61*EI/(2*L/razdelki)^3 * (L/razdelki)^2;          %nova togost vsakic
    Mnova=m/tocke *eye(tocke,tocke);                            %nova masna matrika (diagonalna)
    Knova=zeros(tocke,tocke);                                              %togostna matrika
    Knova(1,1:2)=[2, -1];                       
    Knova(tocke,tocke-1:tocke)=[-1,2];                      
    for i=2:tocke-1                                 
        Knova(i,(i-1):(i+1))=[-1,2,-1];                                     %zapisemo vse primere, ki so vmes
    end
    Knova=kun*Knova;   
    Amatn=Mnova\Knova;                                                  %uporabili bomo metodo eigenvalues, ki zahteva tak vhod
    [obliken,frekn]=eig(Amatn);
    novefrek(inkrement)=sqrt(frekn(3,3));
    inkrement=inkrement+1;
end
figure(9)
plot(segmenti,novefrek,'LineWidth',2,'Marker','s')
grid on
xlabel('Stevilo segmentov [/]'); ylabel('\omega_{03} [rad/s]'); title('Graf konvergence')