function  FRFa_matrix=nama(ShaftNodeSet,E,W,Prop,supa,M)
n = size(Prop,2); %% number of beam sections
% n-1: number of inner points
% n+1: number of beam points
NON = length(ShaftNodeSet);
% max of NON:n+1
FRFa_matrix = zeros(NON);

%% NAM matrix: T, weighting matrix: Add
[aAcl,aAcr] = tma(E,W,Prop,supa,M);
aPl = -aAcl{1}(2,:);aPr = aAcr{n}(2,:);
aAl = aAcl{1}(1,:);aAr = aAcr{n}(1,:);
Ta = zeros(2*n);
for j = 1:n-1
    Ta(2*j:2*j+1,2*j-1:2*j) = aAcr{j};
    Ta(2*j:2*j+1,2*j+1:2*j+2) = -aAcl{j+1};
end
Ta(1,1:2) = aPl;
Ta(2*n,2*n-1:2*n) = aPr;
for j = 1:2*n
    atemp_A(j) = sqrt(1/sum(Ta(j,:).^2));
end
aAdd = diag(atemp_A);

for Nodes = 1:NON
   %% Force input
    Fa_input = zeros(2*n,1);
    if ShaftNodeSet(Nodes)==1
        Fa_input(1) = 1;
    elseif ShaftNodeSet(Nodes)==n+1
        Fa_input(end) = 1;
    else 
        Fa_input(2*(ShaftNodeSet(Nodes)-1)+1)=1;
    end
    %% NAM Matrix
    aC = (aAdd*Ta)\aAdd*Fa_input;
    aC_d = zeros(2,n);
    for j = 1:n
        aC_d(:,j) = aC(2*j-1:2*j);
    end
    
    for j = 1:Nodes
        if ShaftNodeSet(j) == n+1
            FRFa_matrix(Nodes,j)=aAr(1,:)*aC_d(:,end);
        else
            FRFa_matrix(Nodes,j)=aAcl{ShaftNodeSet(j)}(1,:)*aC_d(:,ShaftNodeSet(j));
        end
    end
end

FRFa_matrix = FRFa_matrix+tril(FRFa_matrix,-1)';
end