clc
clear all
% sl sb rb R x b tap
linedata=[1 1 2 0.042 0.168 0.082 1;
2 1 5 0.031 0.126 0.062 1;
3 2 3 0.031 0.126 0.062 1;
4 3 4 0.024 0.136 0.164 1;
5 3 5 0.053 0.210 0.102 1;
6 4 5 0.063 0.252 0.122 1];
% bus type V d Pg Qg Pl Ql Qmin() Qmax
busdata = [ 1 1 1.00 0 0.0 0 0 0 0 0;
2 2 1.00 0 50 0 0 0 0 500;
3 2 1.00 0 100 0 0 0 -500 500;
4 3 1.00 0 0.0 0 115 60 -500 500;
5 3 1.00 0 0.0 0 85 40 -500 500;
];
nline = linedata(:,1);
sb = linedata(:,2); %Sending End Bus column
rb = linedata(:,3); %Receiving End Bus column
R = linedata(:,4); %Resistance(pu) in lines
X = linedata(:,5); %Reactance(pu) in lines
hc = 1i*linedata(:,6); %Capacitive Susceptance(pu)
t = linedata(:,7); %nominal tap ratio
Z = R + 1i*X;y=(1./Z);
sb1=max(sb);sb2=max(rb);
bus = max(sb1,sb2); %number of buses
line = length(nline); %number of lines
ybus=zeros(bus,bus);
for k=1:line %loop for tap ratio
if t(k)~=1
t1 = (1-1/t(k));
t2 = -t1/t(k);
hc(k) = t2*y(k);
y(k) = y(k)/t(k);
end
end
for k=1:line %loop for Ybus
p=sb(k);q=rb(k);
ybus(p,p)=ybus(p,p)+y(k)+hc(k)/2;
ybus(q,q)=ybus(q,q)+y(k)+hc(k)/2;
ybus(p,q)=ybus(p,q)-y(k);
ybus(q,p)=ybus(p,q);
end
busnumber=busdata(:,1); %Bus number
type = busdata(:,2);%bus type
Vmag = busdata(:,3);%Voltage magnitude
V = busdata(:,3).*exp(1i*busdata(:,4));
ANG = busdata(:,4);%bus angle
PG = busdata(:,5);%generated real power
QG = busdata(:,6);%generated reactive power
PL = busdata(:,7);%load real power
QL = busdata(:,8);%load reactive power
BaseMVA =100;
Load = PL + 1i*QL;
Qmin = busdata(:,9)/100; % Minimum Reactive Power Limit
Qmax = busdata(:,10)/100; % Maximum Reactive Power Limit
G = real(ybus); % Conductance
B = imag(ybus); % Susceptance
Psp=(PG-PL)/100;
Qsp=(QG-QL)/100;
Ssp= Psp +1i*Qsp;
PV= find(type==2);
PQ= find(type==3);
NS= [PV;PQ];
M=length(PV); %no. of PV bus
N=length(PQ);%no. of PQ bus
O=length(NS);
dP=zeros(O,1);
dQ=zeros(N,1);
MV = [dP; dQ];
iter=1;
tol=1;
st=clock; % stop the iteration time clock
while (tol>1e-8)
S = V.*conj(ybus*V);
S1 = S-Ssp;
Qc = imag(S);
dP=real(S1(NS));
s=0; %set pv violation
for i = 1:bus
if type(i) == 2
if (Qc(i) > Qmax(i)) || (Qc(i) < Qmin(i))
if Qc(i) < Qmin(i) % Whether violated the lower limit.
Qsp(i) = Qmin(i);
else % No, violated the upper limit.
Qsp(i) = Qmax(i);
end
type(i) = 3; % If Violated, change PV bus to PQ bus.
end
s=s+1;
end
end
PV= find(type==2);
PQ= find(type==3);
NS= [PV;PQ];
M=length(PV); %no. of PV bus
N=length(PQ);%no. of PQ bus
O=length(NS);
dQ=imag(S1(PQ));
B1= -imag(ybus(NS,NS));
B11= -imag(ybus(PQ,PQ));
invB1 = inv(B1);
invB11 = inv(B11);
dth=-invB1*dP; %correction vector : theta
dV=-invB11*dQ; %correction vector : Voltage
MV=[dP; dQ];
ANG(NS) = ANG(NS)+ dth;
Vmag(PQ) = Vmag(PQ)+ dV;
V=Vmag.*exp(1i*ANG);
iter = iter+1;
tol =max(abs(MV));
end
PVnew= find(type==2);
PQnew= find(type==3);
dell = ANG*180/pi;
for m=1:bus
for n=1:bus
I(m,n)=-(V(m)-V(n))*ybus(m,n);
end
I(m,m)=ybus(m,:)*V;
end
for m=1:bus
for n=1:bus
Pow_Flow(m,n)=V(m)*I(m,n)';
end
Pow_Flow(m,m)=V(m)*I(m,m)';
end
s1=find(type==1);
Slack_bus_power=Pow_Flow(s1,s1)*BaseMVA;
for m=1:bus
for n=1:bus
Loss(m,n)=Pow_Flow(m,n)+Pow_Flow(n,m);
end
Loss(m,m)=0;
end
ste=clock; % stop the iteration time clock
lf = [type Vmag dell real(S).*100 Qc.*100];
Losses = real(Loss);
Total_loop_loss=sum((sum(real(Loss))))/2*BaseMVA;
disp('--------------------------------------------------------------------------')
disp(' FDLF Load-Flow Study');
disp('--------------------------------------------------------------------------')
disp(date);
fprintf('Number of iterations : %d \n', iter);
fprintf('Error tolerance considered : %d \n', tol);
fprintf('Solution time : %g sec.\n\n',etime(ste,st));
fprintf(' Ybus Matrix:\n');
disp('==================');
disp(ybus);
disp(" B'");
disp('========');
disp(B1)
disp(' B"');
disp('========');
disp(B11);
fprintf(' BUS Data:\n');
disp('============');
disp(' bus type v delta Pg Qg Pl Q1 Qmin Qmax');
disp(' ----- ----- --- ----- ---- ---- ----- --- ------ ------');
disp(busdata);
fprintf('Load Flow:\n');
disp('=============');
disp(' Bus Type Vpu angle P(MW) Q(MVAr)');
disp(' -------- -------- -------- -------- -------');
disp(lf);
fprintf('Power loss:%d MW\n \n',Total_loop_loss);