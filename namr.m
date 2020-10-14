function  FRFr_matrix=namr(ShaftNodeSet,G,W,Prop,supr,Jt)
n = size(Prop,2); %% number of beam sections
% n-1: number of inner points
% n+1: number of beam points
NON = length(ShaftNodeSet);
% max of NON:n+1
FRFr_matrix = zeros(NON);

%% NAM matrix: T, weighting matrix: Add
[tAcl,tAcr] = tmr(G,W,Prop,supr,Jt);
tPl = -tAcl{1}(2,:);tPr = tAcr{n}(2,:);
tAl = tAcl{1}(1,:);tAr = tAcr{n}(1,:);
Tr = zeros(2*n);
for j = 1:n-1
    Tr(2*j:2*j+1,2*j-1:2*j) = tAcr{j};
    Tr(2*j:2*j+1,2*j+1:2*j+2) = -tAcl{j+1};
end
Tr(1,1:2) = tPl;
Tr(2*n,2*n-1:2*n) = tPr;
for j = 1:2*n
    ttemp_A(j) = sqrt(1/sum(Tr(j,:).^2));
end
tAdd = diag(ttemp_A);

for Nodes = 1:NON
   %% Force input
    Fr_input = zeros(2*n,1);
    if ShaftNodeSet(Nodes)==1
        Fr_input(1) = 1;
    elseif ShaftNodeSet(Nodes)==n+1
        Fr_input(end) = 1;
    else 
        Fr_input(2*(ShaftNodeSet(Nodes)-1)+1)=1;
    end
    %% NAM Matrix
    tC = (tAdd*Tr)\tAdd*Fr_input;
    tC_d = zeros(2,n);
    for j = 1:n
        tC_d(:,j) = tC(2*j-1:2*j);
    end
    
    for j = 1:Nodes
        if ShaftNodeSet(j) == n+1
            FRFr_matrix(Nodes,j)=tAr(1,:)*tC_d(:,end);
        else
            FRFr_matrix(Nodes,j)=tAcl{ShaftNodeSet(j)}(1,:)*tC_d(:,ShaftNodeSet(j));
        end
    end
end

FRFr_matrix = FRFr_matrix+tril(FRFr_matrix,-1)';
end