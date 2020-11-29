function fit = inter(solvecArray,M,N)
[r,c]=size(solvecArray);
for iter=(1:r)
    runiter = inter2(solvecArray(iter,:),M,N);
    f1(iter)=runiter(1);
    f2(iter)=runiter(2);
end
fit = [f1' f2'];
end

