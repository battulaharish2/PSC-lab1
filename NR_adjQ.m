clc;
clear;
%% data
load('busd14');
load('ldata14');
%% ybus formation
l = ldata14(:,1);
nline=length(l);
nbus=length(busd14(:,1));
fb = ldata14(:,2);
tb = ldata14(:,3);
R = ldata14(:,4);
X = ldata14(:,5);
hc = ldata14(:,6);
t = ldata14(:,7);
Z = complex(R,X);
y=(1./Z);
sb1=max(fb);
sb2=max(tb);
bus = max(sb1,sb2);
ybus=zeros(bus,bus);
for k=1:nline
    p=fb(k);
    q=tb(k);
    ybus(p,p)=ybus(p,p)+y(k)+1i*hc(k);
    ybus(q,q)=ybus(q,q)+y(k)+1i*hc(k);
    ybus(p,q)=ybus(p,q)-y(k);
    ybus(q,p)=ybus(p,q);
end
Y=ybus;
%% reading bus data
BMVA=100;
bus=busd14(:,1);
type=busd14(:,2);
v=busd14(:,3);
dell=busd14(:,4);
del=dell*pi/180;
Pg=busd14(:,5)/BMVA;
Qg=busd14(:,6)/BMVA;
Pl=busd14(:,7)/BMVA;
Ql=busd14(:,8)/BMVA;
Qmin=busd14(:,9)/BMVA;
Qmax=busd14(:,10)/BMVA;
P=Pg-Pl;
Q=Qg-Ql;
Psp=P;
Qsp=Q;
G=real(Y);
B=imag(Y);
PQ=find(type==3);
npq=length(PQ);
Tol=1e-3;
err = 1;
Iter=0;
vv=0;
%% NR load flow
st=clock; % start the iteration time clock
while err>Tol
    P=zeros(nbus,1);
    Q=zeros(nbus,1);
    %caluclate p and q
    for i=1:nbus
        for k=1:nbus
            P(i)=P(i)+v(i)*v(k)*(G(i,k)*cos(del(i)-del(k))+B(i,k)*sin(del(i)-del(k)));
            Q(i)=Q(i)+v(i)*v(k)*(G(i,k)*sin(del(i)-del(k))-B(i,k)*cos(del(i)-del(k)));
        end
    end
    % checking q-limit violation
    for i=2:nbus
        if type(i) == 2
            if Q(i) < Qmin(i)         % Whether violated the lower limit.
                Q(i) = Qmin(i);
                type(i) = 3;
                vv = i;
            end
            if Q(i) > Qmax(i)         % No, violated the upper limit
                Q(i) = Qmax(i);
                type(i) = 3;          % If Violated, change PV bus to PQ bus.
                vv = i;
            end
        end
        type;
        PQ=find(type==3);
        npq=length(PQ);
    end
    
    %calculate MISMATCH: change from specified value
    dPa=Psp-P;
    dQa=Qsp-Q;
    k=1;
    dQ=zeros(npq,1);
    for i=1:nbus
        if type(i)==3
            dQ(k,1)=dQa(i);
            k=k+1;
        end
    end
    dP=dPa(2:nbus);
    M=[dP;dQ];
    %jacobian
    %J1- dervative of real power injection with angles
    J1= zeros(nbus-1,nbus-1);
    for i=1:nbus-1
        m=i+1;
        for k=1:nbus-1
            n=k+1;
            if n==m
                for n=1:nbus
                    J1(i,k)=J1(i,k)+v(m)*v(n)*(-G(m,n)*sin(del(m)-del(n))+B(m,n)*cos(del(m)-del(n)));
                end
                J1(i,k)=J1(i,k)-v(m)^2*B(m,m);
            else
                J1(i,k)=v(m)*v(n)*(G(m,n)*sin(del(m)-del(n))-B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    % J2 derivative of real power injection with v
    J2= zeros(nbus-1,npq);
    for i=(1:nbus-1)
        m=i+1;
        for k=1:npq
            n=PQ(k);
            if n==m
                for n=1:nbus
                    J2(i,k)=J2(i,k)+v(n)*(G(m,n)*cos(del(m)-del(n))+B(m,n)*sin(del(m)-del(n)));
                end
                J2(i,k)=J2(i,k)+v(m)^2*G(m,m);
            else
                J2(i,k)=v(m)*(G(m,n)*cos(del(m)-del(n))+B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    %J3 derviative of reactive of reactive power injection with
    %angles
    J3= zeros(npq,nbus-1);
    for i=1:npq
        m=PQ(i);
        for k=1:(nbus-1)
            n=k+1;
            if n==m
                for n=1:nbus
                    J3(i,k)=J3(i,k)+v(m)*v(n)*(G(m,n)*cos(del(m)-del(n))+B(m,n)*sin(del(m)-del(n)));
                end
                J3(i,k)=J3(i,k)-v(m)^2*G(m,m);
            else
                J3(i,k)=v(m)*v(n)*(-G(m,n)*cos(del(m)-del(n))-B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    % J4 derivative of reative power injection with v
    J4= zeros(npq,npq);
    for i=1:npq
        m=PQ(i);
        for k=1:npq
            n=PQ(k);
            if n==m
                for n=1:nbus
                    J4(i,k)=J4(i,k)+v(m)*v(n)*(G(m,n)*sin(del(m)-del(n))-B(m,n)*cos(del(m)-del(n)));
                end
                J4(i,k)=J4(i,k)-v(m)^2*B(m,m);
            else
                J4(i,k)=v(m)*v(n)*(G(m,n)*sin(del(m)-del(n))-B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    J=[J1 J2; J3 J4];
    x=J\M;
    dth=x(1:nbus-1);
    dv=x(nbus:end);
    %updating state vectors
    del(2:nbus)=dth+del(2:nbus);
    k=1;
    for i=2:nbus
        if type(i)==3
            v(i)=dv(k)+v(i);
            k=k+1;
        end
    end
    for i=1:nbus
        E(i)=v(i)*(cos(del(i))+sin(del(i))*1j);
    end
    Iter=Iter+1;
    err=max(abs(M));
end
ste=clock; % stop the iteration time clock
lf = [type abs(E)' rad2deg(angle(E))' P.*100 Q.*100];
%% printing results
disp('--------------------------------------------------------------------------')
disp('                      Newton Rapson Load-Flow Study');
disp('--------------------------------------------------------------------------')
disp(date);
fprintf('Number of iterations       : %d \n', Iter);
fprintf('Error tolerance considered : %d \n', Tol);
fprintf('Solution time              : %g sec.\n',etime(ste,st));
fprintf('Maximum mismatch           : %d\n\n',max(abs(M)));
% fprintf('Jacobian Matrix:\n');
% disp('==================');
% disp(J);
% fprintf('BUS Data:\n');
% disp('============');
% disp('   bus   type   v   delta    Pg   Qg     Pl   Q1  Qmin   Qmax');
% disp('  ----- -----  ---  -----   ---- ----  ----- --- ------ ------')
% disp(busd);
% fprintf('Load Flow:\n');
% disp('=============');
% fprintf('Q-limit voilated at bus: %d\n\n', vv);
disp('   Bus Type    Vpu      angle      P(MW)   Q(MVAr)')
disp('   --------  --------  --------  --------  -------');
disp(lf)
power_loss = sum(lf(:,4));
fprintf('Power loss:%d MW\n \n',power_loss);