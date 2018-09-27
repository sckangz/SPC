function [result]=multiple3(K,s,alpha,beta,gamma,mu)
% s is the true class label.

[m,n,r]=size(K);
Z=eye(n);
c=length(unique(s));
g = ones(1,r)/r;
D = zeros(n);

 for j = 1:r
        D = D + g(j)*K(:,:,j);
    end
%options = optimset( 'Algorithm','interior-point-convex','Display','off');
for i=1:200
    Zold=Z;

    dd = diag(sum(Z));
    L = dd-Z;

    [F, temp, ev]=eig1(L, c, 0);


    for ij=1:n
        for ji=1:n
            all(ji)=(norm(F(ij,:)-F(ji,:)))^2;
        end

       Z(:,ij)=(2*gamma*eye(n)+D)\(alpha*D(ij,:)'-beta/2*all');
       
% we use the free package to solve quadratic equation: http://sigpromu.org/quadprog/index.html
    %    [Z(:,ij),err,lm] = qpas(H,ff,[],[],ones(1,n),1,zeros(n,1),ones(n,1));
% Z(:,ij)=quadprog(H,(beta/2*all'-(alpha-2)*K(:,ij))',[],[],ones(1,n),1,zeros(n,1),ones(n,1),Z(:,ij),options);
    end
    Z(find(Z<0))=0;
    Z= (Z+Z')/2;
    D=zeros(n);
    for j = 1:r
       D= D + g(j)*K(:,:,j);
    end 
   
    h=zeros(12,1);
    for j=1:12
        h(j)=.5*trace(K(:,:,j)-2*alpha*K(:,:,j)+Z'*K(:,:,j)*Z);  
    end
    for j=1:12
        g(j)=(h(j)*sum(1./h))^(-2);
    end
    
    
    
    if i>5 &((norm(Z-Zold)/norm(Zold))<1e-5)
        break
    end

end


actual_ids= kmeans(F, c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
[result] = ClusteringMeasure( actual_ids,s);
