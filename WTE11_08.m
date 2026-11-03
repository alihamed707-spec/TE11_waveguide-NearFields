% This program computes Ex and Ey for TEmn mode waveguide antenna, m=1, n=1
% at height z = 10*lambda in the plane y = 0 (0<x<10*lambda)
% the radius of the antenna is a = 5**lambda
% ph is the angle in azimuth direction



clear all
clc

eps0=8.85418781761e-12;
mu0=pi*4e-7;
c0=1/sqrt(eps0*mu0);
f=c0;
lambda=c0/f;
beta=2*pi/lambda;

m=1;n=1;Xmn=1.841;

M=40;N=10;

a=5*lambda;
z=10*lambda;
y=0;

ph=pi/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



xo=eps:0.1:10;


for t=1:length(xo)
    
    x=xo(t);
    
    rho=sqrt(x^2+y^2);
    r=sqrt(z^2+rho^2);
    
    Ra=sqrt(r^2+a^2);
    
    %Y00=-(1.841/2)*(besselh(0.5,2,beta*r)*sqrt(pi/(2*beta*r))-besselh(0.5,2,beta*Ra)*sqrt(pi/(2*beta*Ra)))/(beta*a);
    
    sumqc=0;
    sumqs=0;
    
    for q=0:N
        
        Fq=((-1)^q)*((Xmn/(beta*a))^(2*q+m))/(gamma(q+1)*gamma(q+m));
        
        
        sumpc=0;
        sumps=0;
        
        for p=0:M
            
            
            Cpq=((Xmn/a)^2)*(rho^2)/((q+m)*(q+m+1)*(p+m)*(p+m+1));
            
            Fp=((beta*rho)^(2*p+m-1))/(gamma(p+1)*gamma(p+m));
            
            L=p+q+m;
            
            FLa=2^(L);
            
            Da0=(gamma(L)/FLa);
            Ga0=((beta*r)^(q-p));
            Ha0=besselh(p-q+0.5,2,beta*r)*sqrt(pi/(2*beta*r));
            
            suma=Da0*Ga0*Ha0;
            
            for va=0:L-1
                
                Da=gamma(L)/gamma(L-va);
                Ea=((beta*a/2)^(2*q+m-va-1))/(2^(va+1));
                Ga= ((a/(2*Ra))^(2*p+m-va-1));
                Ha=besselh(2*p+m-va-0.5,2,beta*Ra)*sqrt(pi/(2*beta*Ra));
                
                
                suma=suma-Da*Ea*Ga*Ha;
                
            end
            
            FLb=2^(L+2); 
            
            Db0=gamma(L+2)/FLb;
            Gb0=(beta*r)^(q-p);
            Hb0=besselh(p-q+0.5,2,beta*r)*sqrt(pi/(2*beta*r));
            
            sumb=Db0*Gb0*Hb0;
            
            for vb=0:L+1
                
                Db=gamma(L+2)/gamma(L-vb+2);
                Eb=((beta*a/2)^(2*q+m-vb+1))/(2^(vb+1));
                Gb=((a/(2*Ra))^(2*p+m-vb+1));
                Hb=besselh(2*p+m-vb+1.5,2,beta*Ra)*sqrt(pi/(2*beta*Ra));
                
                sumb=sumb-Db*Eb*Gb*Hb;
                
            end
            
            Fcpa=cos((m-1)*ph)*suma;
            Fcpb=cos((m+1)*ph)*Cpq*sumb;
            
            Fspa=sin((m-1)*ph)*suma;
            Fspb=sin((m+1)*ph)*Cpq*sumb;
            
            
            sumpc=sumpc+Fp*(Fcpa-Fcpb);
            sumps=sumps+Fp*(Fspa+Fspb);
            
            
            
        end
        
        sumqc=sumqc+Fq*sumpc;
        sumqs=sumqs+Fq*sumps;
        
    end
    
    Ex(t)=j*beta^2*z*sumqs;
    Ey(t)=j*beta^2*z*sumqc;
    
    
end

%save e:\JUNE2023\datasam\WTEmn\TE11020005.mat Ey xo a z ph Xmn m n q p  %From WT1108
%save e:\JUNE2023\datasam\WTEmn\TE11029005.mat Ey xo a z ph Xmn m n q p  %From WT1108

%save e:\JUNE2023\datasam\WTEmn\TE21050010.mat Ey xo a z ph Xmn m n q p  %From WT1108
%save e:\JUNE2023\datasam\WTEmn\TE21059010.mat Ex xo a z ph Xmn m n q p  %From WT1108

%save e:\JUNE2023\datasam\WTEmn\TE21050005.mat Ey xo a z ph Xmn m n q p  %From WT1108
%save e:\JUNE2023\datasam\WTEmn\TE21059005.mat Ex xo a z ph Xmn m n q p  %From WT1108

plot(xo,abs(Ex)/max(abs(Ex)),'r.-')
        
