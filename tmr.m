function [Acl,Acr]=tmr(G,w,Prop,supr,Jr)
n=size(Prop,2);
% Prop = [I;A;l;kappa;Rho]
I = Prop(1,:)*2;% rotational inertia *2
A = Prop(2,:);
l = Prop(3,:);
Rho = Prop(5,:);    


lambda = zeros(n,1);
for i=1:n
    lambda(i) = w * sqrt( Rho(i) / G(i) );
end
Ktr = supr(end);J = Jr(end); 
Kbtr = -(Ktr - J*w^2); 
Acl{n}=zeros(2);Acr{n}=zeros(2);
for i = 1:n    
    Acl{i} = [1                 0;
            -supr(i)+w^2*Jr(i) G(i)*I(i)*lambda(i)];
    Acr{i} = [cos(lambda(i) * l(i))                         sin(lambda(i) * l(i));
              -G(i)*I(i)*lambda(i)*sin(lambda(i) * l(i))    G(i)*I(i)*lambda(i)*cos(lambda(i) * l(i))];
end
if i == n
    Acr{i}(2,:)=[-G(i)*I(i)*lambda(i)*sin(lambda(i) * l(i)) - Kbtr*cos(lambda(i)*l(i)) ...
                G(i)*I(i)*lambda(i)*cos(lambda(i) * l(i)) - Kbtr*sin(lambda(i)*l(i))];
end


