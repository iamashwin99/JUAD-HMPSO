function out=fitness(Gjb,Gbj,Gib,Gbi,Gji,Gij,Gii,rhou, rhod,pj,pi,pbs, gammaireq, gammajrequ,gammajreqd)
%Hi shashwat
    M=length(Gjb);%number of cellular users
    N=length(Gbi);%number of DUEs
    gammau=zeros(M);%gamma uplink for cellular
    gammad=zeros(N);%gamma downlink for cellular
    gammai=zeros(N);%gamma for DUE
    N0=1;%variance of white Gaussian Noise
    for i=(1:M)
        deno=0;
        deno1=0;
        for j=(1:N)
            deno=deno+rhou(j,i)*pi(j)*Gib(j);
            deno1=deno1+rhod(j,i)*pi(j)*Gij(j,i);
        end
        gammau(i)=pj(i)*Gjb(i)/(deno+N0);
        gammad(i)=pj(i)*pbs(i)*Gbj(j)/(deno1+N0);
    end
    for i=(1:N)
        deno=0;
        for j=(1:M)
            deno=deno+rhou(i,j)*pj(j)*Gji(i,j)+rhod*pbs(j)*Gbi(i);
        end
        gammai(i)=pi(i)*Gii(i)/(deno+N0);
    end
    Ru=log2(1+gammau);%uplink data rate cellular
    Rd=log2(1+gammad);%downlink data rate cellular
    Ri=log2(1+gammai);%data rate DUE
    out=sum(Ru)+sum(Rd)+sum(Ri);%data rate
    count=0;
    test=true;
    
    for j=(1:M)
        if(gammau(j)<gammajrequ(j) || gammad(j)<gammajreqd(j))
            test=false;
            count=count+1;
            break;
        end
    end
    for i=(1:N)
        if(gamma(i)<gammaireq(i)) %% i is for devices j for 
            test=false;
            count=count+1;
            break;
        end
    end
    
    
    
    
end