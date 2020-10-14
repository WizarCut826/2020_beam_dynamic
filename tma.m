function [Acl,Acr]=tma(E,w,Prop,supa,M)
n=size(Prop,2);
% Prop = [I;A;l;kappa;Rho]
A = Prop(2,:);
l = Prop(3,:);
Rho = Prop(5,:);    


lambda = zeros(n,1);
for i=1:n
    lambda(i) = w * sqrt( Rho(i) / E(i) );
end
Ktr = supa(end);Mr = M(end); 
Kbtr = -(Ktr - Mr*w^2); 
Acl{n}=zeros(2);Acr{n}=zeros(2);
for i = 1:n    
    Acl{i} = [1                 0;
            -supa(i)+w^2*M(i) E(i)*A(i)*lambda(i)];
    Acr{i} = [cos(lambda(i) * l(i))                         sin(lambda(i) * l(i));
              -E(i)*A(i)*lambda(i)*sin(lambda(i) * l(i))    E(i)*A(i)*lambda(i)*cos(lambda(i) * l(i))];
end
if i == n
    Acr{i}(2,:)=[-E(i)*A(i)*lambda(i)*sin(lambda(i) * l(i)) - Kbtr*cos(lambda(i)*l(i)) ...
                E(i)*A(i)*lambda(i)*cos(lambda(i) * l(i)) - Kbtr*sin(lambda(i)*l(i))];
end