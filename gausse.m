A = [5 1.2 0.3 -0.6;1.2 6 -0.4 0.9;0.3 -0.4 8 1.7;-0.6 0.9 1.7 10];
B = [0.063; -0.6358; 0.5937; -0.1907];
%disp('Agumented matrix is:');
a = [A,B];
[m,n]=size(a);
%disp('Agumented matrix during iterations is:');
for j=1:m-1
    for z=2:m
        if a(j,j)==0
            t=a(j,:);
            a(j,:)=a(z,:);
            a(z,:)=t;
        end
    end
    for i=j+1:m
        a(i,:)=a(i,:)-a(j,:)*(a(i,j)/a(j,j));
        
    end
end
x=zeros(1,m);
%disp('By back substituion we get the solution as');
for i=m:-1:1
    sum=0;
    for k=2:m
        sum=sum+a(i,k)*x(k);
    end
    x(i)=(a(i,n)-sum)/a(i,i);
    
end
disp('Solution of Gauss elimination method:');
x'

