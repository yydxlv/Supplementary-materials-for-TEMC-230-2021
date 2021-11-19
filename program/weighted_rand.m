
function a=weighted_rand(N,Scope,num)

a=zeros(N,1);
n=0;
f1=@(t) sin(t);

tt=linspace(0,Scope,1001);
ff=f1(deg2rad(tt));
s=trapz(tt,ff);  
ff=ff/s;   
 
 
while n<N
    t=rand(1)*Scope;
    f=f1(deg2rad(t))/s;
    r=rand(1);  
    if r<=f   
        n=n+1;
        a(n)=t;
    end
end
 
num=100;        
[x,c]=hist(a,num);    
dc=Scope/num;      
x=x/N/dc;        
