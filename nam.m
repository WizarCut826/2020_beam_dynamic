% nam: return FRF matrix at frequency W
% input: ShaftNodeSet: calculate node index
%       E: Young's module  G:shear module
%       W: frequency       Prop: Property matrix of beam 
%       sup:support(stiffness)
%       M: lumped mass     J; lumped innertia
function  FRF_matrix=nam(ShaftNodeSet,E,G,W,Prop,sup,M,J)
n = size(Prop,2); %% number of beam sections
% n+1: number of beam points
NON = length(ShaftNodeSet);
% max of NON:n+1.
FRF_matrix = zeros(NON);

%% NAM matrix: T, weighting matrix: Add
[Acl,Acr] = tm(E,G,W,Prop,sup,M,J);
T = zeros(4*n);

Al = Acl{1}(1:2,:);Ar = Acr{n}(1:2,:);
Pl = -Acl{1}(3:4,:);Pr = Acr{n}(3:4,:);
for j = 1:n-1
    T(4*j-1:4*j+2,4*j-3:4*j) = Acr{j};
    T(4*j-1:4*j+2,4*j+1:4*j+4) = -Acl{j+1};
end
T(1:2,1:4) = Pl;
T(end-1:end,end-3:end) = Pr;
for j = 1:4*n
    temp_A(j) = sqrt(1/ sum(T(j,:).^2));
end
Add = diag(temp_A);

for Nodes = 1:NON
   %% Force input
    F_input = zeros(4*n,1);
    if ShaftNodeSet(Nodes)==1
        F_input(2) = 1;
    elseif ShaftNodeSet(Nodes)==n+1
        F_input(end) = 1;
    else 
        F_input(4*(ShaftNodeSet(Nodes)-1)+2)=1;
    end
    %% NAM Matrix
    C = (Add*T)\Add*F_input;
    C_d = zeros(4,n);
    for j = 1:n
        C_d(:,j) = C(4*j-3:4*j);
    end
    
    for j = 1:Nodes
        if ShaftNodeSet(j) == n+1
            FRF_matrix(Nodes,j)=Ar(1,:)*C_d(:,end);
        else
            FRF_matrix(Nodes,j)=Acl{ShaftNodeSet(j)}(1,:)*C_d(:,ShaftNodeSet(j));
        end
    end
end

FRF_matrix = FRF_matrix+tril(FRF_matrix,-1)';
end