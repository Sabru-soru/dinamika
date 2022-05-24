clc; clear all; close all
% Podatki 
k = 890; L = 0.5; m = 9.3; g = 9.81;
F = linspace(0,1260,101);
syms fi
for i=1:101
Fi = eval(solve(F(i)*L*cos(fi) - m*g*L*cos(fi) - 2*k*L^2*(1 - cos(fi))*sin(fi) == 0, fi,'Real',true)); %samo realne rešitve
Fi1 = (Fi(2))*180/pi     %v stopinjah
kot(i)=Fi1;
end
fprintf('Kot fi znaša %1.2f ° v statièni legi, èe zunanja sila ni prisotna. \n',kot(1))
plot(F,kot,'s','MarkerSize',2,'MarkerFaceColor','b')
grid on; xlabel('Sila F [N]'); ylabel('Kot fi [°]'); title('Graf odvisnosti kota fi od sile F, fi(F) ')
%Izraèun mejne sile
fi=0;
clearvars F; syms F
Fm=eval(solve(F*L*cos(fi) - m*g*L*cos(fi) - 2*k*L^2*(1 - cos(fi))*sin(fi) == 0, F));
fprintf('Mejna sila, kjer bo kot fi enak niè znaša F=%1.3f N. \n',Fm)