% sl no fbus tbus r(pu) x(pu) b(pu) tap
Ldata=[
1 1 2 0.042 0.168 0.082 1
2 1 5 0.031 0.126 0.062 1
3 2 3 0.031 0.126 0.062 1
4 3 4 0.024 0.136 0.164 1
5 3 5 0.053 0.21 0.102 1
6 4 5 0.063 0.252 0.122 1
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
%%Form ITAGF and ITAGT
ITAGF(1)= 1; 
ITAGT(1) = NLCOUNT(1);
for i = 2:nbus;
    ITAGF(i)= ITAGT(i-1)+1;
    ITAGT(i)= ITAGT(i-1)+NLCOUNT(i);
end
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
