function [result]=simple(K,s,alpha,beta,gamma)
% s is the true class label.

[m,n]=size(K);
Z=eye(n);
c=length(unique(s));
A=inv(K+2*gamma*eye(n));
%options = optimset( 'Algorithm','interior-point-convex','Display','off');
for i=1:200
    Zold=Z;

    D = diag(sum(Z));
    L = D-Z;

    [F, temp, ev]=eig1(L, c, 0);


    for ij=1:n
        for ji=1:n
            all(ji)=(norm(F(ij,:)-F(ji,:)))^2;
        end

       Z(:,ij)=A*(alpha*K(ij,:)'-beta/2*all');
       
% we use the free package to solve quadratic equation: http://sigpromu.org/quadprog/index.html
    %    [Z(:,ij),err,lm] = qpas(H,ff,[],[],ones(1,n),1,zeros(n,1),ones(n,1));
% Z(:,ij)=quadprog(H,(beta/2*all'-(alpha-2)*K(:,ij))',[],[],ones(1,n),1,zeros(n,1),ones(n,1),Z(:,ij),options);
    end
    Z(find(Z<0))=0;
    Z= (Z+Z')/2;
    if i>5 &((norm(Z-Zold)/norm(Zold))<1e-5)
        break
    end

end


actual_ids= kmeans(F, c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
[result] = ClusteringMeasure( actual_ids,s);
