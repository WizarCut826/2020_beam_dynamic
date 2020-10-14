% Simple shaft bending FRF without supports
clear all
%% Frequency Range
W = 1:0.1:100;W = W*2*pi;
% W = 1:200;W = W*2*pi;
%% Initializing shaft material property,    
n = 3;
E = 2.1e11;eta = 0.0;
E = E.*(1+1i*eta);E = E*ones(n,1);
V = 0.3;
G = 0.5*E./(1.+V);    % shear module

R = [0.01 0.02 0.01];D = 2*R;
r = [0 0 0]; d = 2*r;

l = [1 1 1];
rho = 7850;
Rho = rho*ones(1,n);


I = pi*(D.^4-d.^4)/64;
A = pi*(D.^2-d.^2)/4;

kappa = 6*(1.+V).*(1+(d./D).^2).^2 ./...
    ((7+6.*V).*(1+(d./D).^2).^2+(20+12.*V).*(d./D).^2) ;  % kappa: Timoshenko shear coefficient.
Prop = [I;A;l;kappa;Rho];                                 % Prop:Physical property matrix
%% Support,BC，mass and Forces
sup = zeros(n,2);
sup(1:2,2)=1e7*ones(1,2);

sup = [zeros(1,2);sup];
supr = zeros(n+1,1);
supa = zeros(n+1,1);

M = zeros(n+1,1); J = zeros(n+1,1); Jt = zeros(n+1,1);
%M:集中质量 J:集中转动惯量 Jr：扭转集中转动惯量
M(1)= 40; M(4) = 50;
Jt(1) = 4; Jt(4) = 10;
%% Node Set
ShaftNodeSet=1:n+1;

NON = length(ShaftNodeSet);
Hshaft{length(W)}=zeros(NON);
Hrshaft{length(W)}=zeros(NON);
Hashaft{length(W)}=zeros(NON);

%% Frequency Response Function
h = waitbar(0,'please wait');
for ww = 1:length(W) %
    %% Transverse
 	Hshaft{ww}=nam(ShaftNodeSet,E,G,W(ww),Prop,sup,M,J);
    Hrshaft{ww}=namr(ShaftNodeSet,G,W(ww),Prop,supr,Jt);
    Hashaft{ww}=nama(ShaftNodeSet,E,W(ww),Prop,supr,M);
    %% waitbar
    str = ['...',num2str(ww/length(W)*100),'%'];
    waitbar(ww/length(W),h,str);
end
delete(h);
%% Plot
for i = 1:length(W)
    disp1(i) = Hashaft{i}(1,1);
    disp2(i) = Hashaft{i}(1,2);
    disp3(i) = Hashaft{i}(1,3);
    disp4(i) = Hashaft{i}(1,4);
end
figure;
plot(W/2/pi,(disp1),'linewidth',1);hold on;
plot(W/2/pi,(disp2),'linewidth',2);hold on;
plot(W/2/pi,(disp3),'linewidth',1);hold on;
plot(W/2/pi,(disp4),'linewidth',1);hold on;
% set(gca,'yscale','log')
% plot(W/2/pi,abs(disp1));
% semilogy(W/2/pi,abs(disp))