clc;
clear all;
disp('Coefficient matrix of set of equations:')
A = [5 1.2 0.3 -0.6;
     1.2 6 -0.4 0.9; 
     0.3 -0.4 8 1.7; 
   -0.6 0.9 1.7 10]
B = [0.0630, -0.6358, 0.5937, -0.1907]';
N = length( A );
L = zeros( N, N );
for i=1:N
    L(i,i)=sqrt(A(i,i)-L(i,:)*L(i,:)');

    for j=(i + 1):N
        L(j,i)=(A(j,i)-L(i,:)*L(j,:)')/L(i,i);
    end
end
disp('Lower and Upper triangular matrices after decomposition:')
L
U = L'
B
disp('Forward substitution on L matrix')
disp("Intermittant solution obtaioned is 'Y'")
%y = L \ b
Y = zeros(N,1);
Y(1) = B(1)/L(1,1);
for k = 2:N
    Y(k) = (B(k) -L(k,1:k-1)*Y(1:k-1))/L(k,k);
end
Y
disp('Backward substitution on U matrix')
%x = U \ y
a = [U , Y]
[N,M] = size(a);
x=zeros(N,1);
for i = N:-1:1
    sum=0;
    for k=2:N
        sum=sum+a(i,k)*x(k);
    end
    x(i)=(a(i,M)-sum)/a(i,i);
end
disp('Solution set by cholesky method is:')
x