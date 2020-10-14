%function tm: derive lateral transfer matrix of beam section
% Acl,Acr: left & right point transfer matrix of a beam section

function [Acl,Acr] = tm(E,G,w,Prop,sup,M,J)
n=size(Prop,2); 
I = Prop(1,:);      A = Prop(2,:);
l = Prop(3,:);      kappa = Prop(4,:);
rho = Prop(5,:);


alpha = zeros(n,1);beta = zeros(n,1);gamma = zeros(n,1);
lambda = zeros(n,1);lambda_b = zeros(n,1);
q = zeros(n,1);q_b = zeros(n,1);
a1 = zeros(n,1);a2 = zeros(n,1);a3 = zeros(n,1);a4 = zeros(n,1);
for i=1:n
    alpha(i) =  rho(i) * w^2 / E(i);
    beta(i)  =  rho(i) * w^2 / kappa(i) / G(i);
    gamma(i) =  rho(i) * A(i) * w^2 / E(i) / I(i);
    lambda(i) = (sqrt(((alpha(i)-beta(i))/2)^2+gamma(i))-(alpha(i)+beta(i))/2)^0.5;
    lambda_b(i) = (sqrt(((alpha(i)-beta(i))/2)^2+gamma(i))+(alpha(i)+beta(i))/2)^0.5;
    q(i) = rho(i) * w^2 / (kappa(i) * G(i) * lambda(i)) + lambda(i);
    q_b(i) = rho(i) * w^2 / (kappa(i) * G(i) * lambda_b(i)) - lambda_b(i);
    a1(i) = E(i)*I(i)*lambda(i)*q(i);                a2(i) = E(i)*I(i)*lambda_b(i)*q_b(i);
    a3(i) = kappa(i)*G(i)*A(i)*(q(i)-lambda(i));     a4(i) = kappa(i)*G(i)*A(i)*(q_b(i)+lambda_b(i));
end

% [Ktr,Krr,Jr,Mr] = bc(2,:);
Ktr = sup(end,2);Krr = sup(end,1);
Jr = J(end); Mr = M(end);

Kbrr = (Krr - Jr*w^2); Kbtr = (Ktr - Mr*w^2); 
for i=1:n
    Acl{i}(:,:)= [ 1        0       1        0;
                   0        q(i)    0       -q_b(i);
                   a1(i)    0       a2(i)    0;
                   0        a3(i)   0       -a4(i)];
    Acl{i}(3,2) = -(-w^2*J(i)+sup(i,1))*q(i);
    Acl{i}(3,4) = (-w^2*J(i)+sup(i,1))*q_b(i);
    Acl{i}(4,1) = -w^2*M(i) + sup(i,2);
    Acl{i}(4,3) = -w^2*M(i) + sup(i,2);
    
    Acr{i}(:,:)=[cosh(lambda(i)*l(i))           sinh(lambda(i)*l(i))         cos(lambda_b(i)*l(i))        sin(lambda_b(i)*l(i));
                q(i)*sinh(lambda(i)*l(i))       q(i)*cosh(lambda(i)*l(i))    q_b(i)*sin(lambda_b(i)*l(i)) -q_b(i)*cos(lambda_b(i)*l(i));
                a1(i)*cosh(lambda(i)*l(i))      a1(i)*sinh(lambda(i)*l(i))   a2(i)*cos(lambda_b(i)*l(i))  a2(i)*sin(lambda_b(i)*l(i));
                a3(i)*sinh(lambda(i)*l(i))      a3(i)*cosh(lambda(i)*l(i))   a4(i)*sin(lambda_b(i)*l(i))  -a4(i)*cos(lambda_b(i)*l(i))];
    if i == n
        Acr{i}(3:4,:)=...
            [a1(i)*cosh(lambda(i)*l(i))+Kbrr*q(i)*sinh(lambda(i)*l(i))   a1(i)*sinh(lambda(i)*l(i))+Kbrr*q(i)*cosh(lambda(i)*l(i)) ...
            a2(i)*cos(lambda_b(i)*l(i))+Kbrr*q(i)*sin(lambda_b(i)*l(i))  a2(i)*sin(lambda_b(i)*l(i))+Kbrr*q(i)*cos(lambda_b(i)*l(i));
            a3(i)*sinh(lambda(i)*l(i))-Kbtr*cosh(lambda(i)*l(i))         a3(i)*cosh(lambda(i)*l(i))-Kbtr*sinh(lambda(i)*l(i)) ... 
            a4(i)*sin(lambda_b(i)*l(i))-Kbtr*cos(lambda_b(i)*l(i))      -a4(i)*cos(lambda_b(i)*l(i))-Kbtr*sin(lambda_b(i)*l(i))];
    end    
end



