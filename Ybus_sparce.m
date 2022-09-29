clc;
clear all;

% sl no fbus tbus r(pu) x(pu) b(pu) turn ratio
Ldata = [
    1 1 4 0      0.0576 0      1
    2 2 7 0      0.0625 0      1
    3 3 9 0      0.0586 0      1
    4 4 5 0.010  0.085  0.088  1
    5 4 6 0.017  0.092  0.079  1
    6 5 7 0.03   0.161  0.153  1
    7 6 9 0.039  0.170  0.179  1
    8 7 8 0.0085 0.072  0.0745 1
    9 8 9 0.0119 0.1008 0.1045 1
    ];

fbus=Ldata(:,2);
tbus=Ldata(:,3);
r=Ldata(:,4);
x=Ldata(:,5);
b=Ldata(:,6);
b;
t=Ldata(:,7);
z=r+i*x;
y=1./z;
y;
b=i*b;

nbus=max(max(fbus),max(tbus));
Y=zeros(nbus,nbus);
nline = length(Ldata);
%%Form NLCOUNT
NLCOUNT = zeros(1,nbus);
for k = 1:nline;
    P = fbus(k); 
    Q = tbus(k);
    NLCOUNT(P)= NLCOUNT(P)+1;
    NLCOUNT(Q)= NLCOUNT(Q)+1;
end
%NLCOUNT=[2 4 3 3 3 2 3];
%%Form ITAGF and ITAGT
ITAGF(1)= 1; 
ITAGT(1) = NLCOUNT(1);
for i = 2:nbus;
    ITAGF(i)= ITAGT(i-1)+1;
    ITAGT(i)= ITAGT(i-1)+NLCOUNT(i);
end
%ITAGF = [1 3 7 10 13 16 18];
%ITAGT = [2 6 9 12 15 17 20];
%%Form YPP
YPP = zeros(1,nbus);
for k = 1 : nline;
    P = fbus(k);
    Q = tbus(k);
    YPP(P) = YPP(P)+t(k)^2*y(k)+i*b(k)/2;
    YPP(Q) = YPP(Q)+y(k)+i*b(k)/2;
end
%%Form YPQ, ADJQ and ADJL
YPQ=zeros(1,2*nline);
ADJQ=zeros(1,2*nline);
ADJL=zeros(1,2*nline);
IBUS=zeros(1,nbus);
for k = 1 : nbus;
P = fbus(k); 
Q = tbus(k);
LPQ = ITAGF(P)+IBUS(P);
LQP = ITAGF(Q)+IBUS(Q);
YPQ(LPQ) = -t(k)*y(k);
YPQ(LQP) = -t(k)*y(k);
ADJQ(LPQ) = Q;
ADJQ(LQP) = P;
ADJL(LPQ) = k;
ADJL(LQP) = k;
IBUS(P) = IBUS(P)+1;
IBUS(Q) = IBUS(Q)+1;
end
YPP
YPQ
NLCOUNT
ITAGF
ITAGT