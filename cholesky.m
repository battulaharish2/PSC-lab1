% % Choleski's Factorization
clc;
clear all;
A = [5 1.2 0.3 -0.6; 
     1.2 6 -0.4 0.9; 
     0.3 -0.4 8 1.7; 
   -0.6 0.9 1.7 10];
N = length(A);
L = zeros(N, N);
U = zeros(N, N);
L(1,1) = sqrt(A(1,1));
U(1,1) = L(1,1);
for a= 2:N
    L(a,1) = A(a,1)/L(1,1);
    U(1,a) = A(1,a)/L(1,1);
end
for i = 2:N
    for j = 1:N
        if i==j
            L(j,i)=sqrt(A(j,i)-L(j,1:i-1)*U(1:i-1,i));
            U(j,i) = L(j,i);
        else
            L(j,i) = (A(j,i)-L(j,1:i-1)*U(1:i-1,i))/L(i,i);
        end
    end
    for k=i+1:N
        U(i,k) = (A(i,k)-L(i,1:i-1)*U(1:i-1,k))/L(i,i);
    end
end
L
U
B = [0.063; -0.6358; 0.5937; -0.1907];
y = L \ B;
x = U \ y;
disp('Solution set by cholesky method is:')
Answer = x'