function out=fitness(Gjb,Gbj,Gib,Gbi,Gji,Gij,Gii,rhou, rhod,pj,pi,pbs, gammaireq, gammajrequ,gammajreqd)
%%%%%% Calculate constriatints
%function out=fitness(Gjb,Gbj,Gib,Gbi,Gji,Gij,Gii,rhou, rhod,pj,pi,pbs, gammaireq, gammajrequ,gammajreqd,pimax,pjmax,pbsmax)


count=0;
test=true;

M=length(Gjb);%number of cellular users  and iterated using j
N=length(Gbi);%number of DUEs and iterated using i
gammau=zeros(1,M);%gamma uplink for cellular
gammad=zeros(1,M);%gamma downlink for cellular
gammai=zeros(1,N);%gamma for DUE
N0=0.01;%variance of white Gaussian Noise
for i=(1:M)
    deno=0;
    deno1=0;
    for j=(1:N)
        deno=deno+rhou(j,i)*pi(j)*Gib(j);
        deno1=deno1+rhod(j,i)*pi(j)*Gij(j,i);
    end
    gammau(i)=pj(i)*Gjb(i)/(deno+N0);
    gammad(i)=pj(i)*pbs(i)*Gbj(i)/(deno1+N0);
end

for i=(1:N)
    deno=0;
    for j=(1:M)
        deno=deno+rhou(i,j)*pj(j)*Gji(i,j)+rhod(i,j)*pbs(j)*Gbi(i);
    end
    gammai(i)=pi(i)*Gii(i)/(deno+N0);
end

%%%Condition 1
for j=(1:M)
    if(gammau(j)<gammajrequ || gammad(j)<gammajreqd)
        test=false;
        count=count+1;
        ;
        break;
        
    end
end
%%%Condition 2
for i=(1:N)
    if(gammai(i)<gammaireq) %% i is for devices j for 
        test=false;
        count=count+1;
        ;
        break;
    end
end
%%%Condition 3
% for i=(1:N)
%     if(pi(i)<0 || pi(i)>pimax)
%         test=false;
%         count=count+1;
%         ;
%         break;
%     end
% end
%%%Condition 4
% for j=(1:M)
%     if(pj(j)<0 || pj(j)>pjmax)
%         test=false;
%         count=count+1;
%         ;
%         break;
%     end
% end
%%%Condition 5
% for j=(1:M)
%     if(pbs(j)<0 || pbs(j)>pbsmax)
%         test=false;
%         count=count+1;
%         ;
%         break;
%     end
% end
%%%Condition 6
for j=(1:M)
    sumrhod=0;
    sumrhou=0;
    for i=(1:N)
        sumrhod=sumrhod+rhod(i,j);
        sumrhou=sumrhou+rhou(i,j);
    end

    if( (sumrhod>1 || sumrhod<0 ) || (sumrhou>1 || sumrhou<0 ) )
        test=false;
        count=count+1;
        ;
        break;
    end
end
%%%Condition 7
for i=(1:N)
    sumrhod=0;
    sumrhou=0;
    for j=(1:M)
        sumrhod=sumrhod+rhod(i,j);
        sumrhou=sumrhou+rhou(i,j);
    end

    if( sumrhod*sumrhou ~=0 )
        test=false;
        count=count+1;
        ;
        break;
    end
end
%%%Condition 8


%%%%% Calculation of R


Ru=log2(1+gammau);%uplink data rate cellular
Rd=log2(1+gammad);%downlink data rate cellular
Ri=log2(1+gammai);%data rate DUE
Fit=-1*(sum(Ru)+sum(Rd)+sum(Ri));%data rate
if(count == 0)
    
    X=['In run Fit: ',num2str(Fit),'count: ',num2str(count)];
    %disp(X)
end
out = [Fit count];
end