function[MM]=gammaMatrix(alpha,K,N)
  MM=zeros(N);
  for i=1:K
      for j=1:i
          MM(i,j)=((-1)^(i-j))*Cgamma(alpha,i-j);
      end
  end
  for i=K+1:N
      for j=i-K+1:i
          MM(i,j)=((-1)^(i-j))*Cgamma(alpha,i-j);
      end
  end 
end