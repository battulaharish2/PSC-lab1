 clear all;
 clc;
 display('----------Ybus formation by direct inspection--------');

% display('Ybus for the sample problem considred in theoretical calculation is')
% Ldata = [
%     1 1 2 0.06 0.18 0.05
%     2 1 3 0.02 0.06 0.06
%     3 2 3 0.04 0.12 0.05
%     ];
display('Ybus for the IEEE 9-bus system is')
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
fb=Ldata(:,2);
tb=Ldata(:,3);
r=Ldata(:,4);
x=Ldata(:,5);
b=Ldata(:,6);
%tap = Ldata(:,7);
z=r+i*x;
y=1./z;
b=i*b;
nbus=max(max(fb),max(tb));
nline = size(Ldata,1);
Y=zeros(nbus,nbus);
% off diagonal element

for k=1:nline
    Y(fb(k),tb(k))=Y(fb(k),tb(k))-y(k);
    Y(tb(k),fb(k))=Y(fb(k),tb(k));
end

% diagonal element

for m=1:nbus
    for n=1:nline
        if fb(n)==m 
            Y(m,m)=Y(m,m)+y(n)+b(n);
        elseif tb(n)==m
            Y(m,m)=Y(m,m)+y(n)+b(n);
        end
    end
end
Y