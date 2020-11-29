function fit = inter2(solvec,M,N)
%{
...
Solution is of len 4m+n m-cue n-due
1st set of m - givesrhou s(j) = i +> rhou (i,j) = 1 (ith device and jth cel)
              round(i) ;add constrint i<m+1
2nd set of m - gives rhod s(j) = i +> rhod (i,j) = 1 (ith device and jth cel)
              round(i) ;add constrint i<m+1

3rd set of m - Pj  s(j) = a +> Pj(j) = a  ///for cue

4th set of m - Pbj  s(j) = a +> Pbj(j) = a ///for cue

last set of n Pi s(j) = a +> Pi(j) = a ///Due


Gjb, channel gain between the BS and CUE j
Gbj, between DUE i and the BS
Gib, between DUE i and the BS
Gbi, between the BS and DUE i
Gji, between CUE j and DUE i
Gij, between DUE i and CUE j
Gii  between the transmitter and receiver of DUE i



...
...
%}
rhou=zeros(N,M);
rhod=zeros(N,M);
pj=zeros(1,M);
pbs=zeros(1,M);
pi=zeros(1,N);
for iter=(1:M)
    if solvec(iter)<0.5
        ;
    elseif solvec(iter)>N+0.5
           ;
    else
        K=round(solvec(iter));
        rhou(K,iter)=1;
    end
    if solvec(iter+M)<0.5
        ;
    elseif solvec(iter+M)>N+0.5
            ;
    else
        K=round(solvec(iter+M));
        rhod(K,iter)=1;
    end
    if solvec(2*M+iter)<0
        pj(iter)=0;
    else
        pj(iter)=solvec(2*M+iter);
    end
    if solvec(3*M+iter)<0
        pbs(iter)=0;
    else
        pbs(iter)=solvec(3*M+iter);
    end
end
for iter=(1:N)
    if solvec(4*M+iter)<0
        pi(iter)=0;
    else
        pi(iter)=solvec(4*M+iter);
    end
end
Gjb=ones(1,M);
Gbj=Gjb;
Gib=ones(1,N);
Gbi=ones(1,N);
Gji=ones(N,M);
Gij=ones(N,M);
Gii=ones(1,N);
gammaireq=0.01;
gammajrequ=0.01;
gammajreqd=0.01;



fit=fitness(Gjb,Gbj,Gib,Gbi,Gji,Gij,Gii,rhou, rhod,pj,pi,pbs, gammaireq, gammajrequ,gammajreqd);

end

