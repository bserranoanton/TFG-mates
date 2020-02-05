function [phi] = Q9analytic(k1,k2onu,k3onu2,N)
    
    phi = zeros(1,51);
    phi(1) = 0.03;
    phi(2) = phi(1)*((k1*N)/(k1+k2onu*(N-1)+k3onu2*(N-1)*(N-2)));
   

    %recursion
    for n = 2:49
        coeffn = k1*N + 2*k2onu*n*(N-n) + k3onu2*n*(N-n)*(N-2);
        coeffnminus1 = k1*(N-n+1)+k2onu*(n-1)*(N-n+1)+k3onu2*(n-1)*(n-2)*(N-n+1); 
        coeffnplus1 = k1*(n+1)+k2onu*(n+1)*(N-n-1)+k3onu2*(n+1)*(N-n-1)*(N-n-2);
        phi(n+1)=  (phi(n)*coeffn - phi(n-1)*coeffnminus1)/coeffnplus1;
    end
    %normalise
  % phi=phi/sum(phi);

end
