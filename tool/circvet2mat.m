
function N = circvet2mat(x)
n=length(x);
N=zeros(n);
for j=1:n
N(j,:)=circshift(x,[0,j-1]);
end