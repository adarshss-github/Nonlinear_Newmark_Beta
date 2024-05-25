function [ ud fs ta] = NBNonLin(m,zhi,Tn,uga,dt,fy,tol)

n = length(uga) ;
ta = 0:dt:(n-1)*dt ;

wn = (2*pi)/Tn ;
c = 2*m*wn*zhi;
k = m*(wn)^2 ;
alpha = 0.5 ;
beta = 1/6 ;
ud = zeros(n,1) ;
uv = zeros(n,1) ;
ua = zeros(n,1) ;
fs = zeros(n,1);
p = zeros(n,1) ;
p = -m.*9.81.*uga ;
ko = k ;

ud(1,1) = 0 ; % Assumed
uv(1,1) = 0 ; % Assumed
fs(1,1) = statedeterEPP(k,fy,0,0,0) ; %Assumed


ua(1,1) = ( p(1,1)-c*uv(1,1) - fs(1,1) )/m ;
a = (1/(beta*dt))*m + (alpha/beta)*c ;
b = (1/(2*beta))*m + dt*((alpha/(2*beta)) -1 )*c ;
kdash = a/dt ;
punb = 0 ;

for i = 1:1:(n-1)
    
    dp = p(i+1) - p(i) ;
    punb = dp + a*uv(i,1) + b*ua(i,1) + punb ;
    kt = tanstiffEPP(k,fy,fs(i,1)) ;
    khat = kt + kdash ;
    udtrial = 0 ;
    fdtrial = 0 ;
    dud = 0 ;
    
    while abs(punb)> tol
        
        dud = dud + (punb/khat) ;
        udtrial = dud + ud(i,1) ;
        duv = (alpha/(beta*dt))*dud - (alpha/beta)*uv(i,1) + dt*(1- (alpha/(2*beta)) )*ua(i,1) ;
        dua = (1/(beta*dt^2))*dud - (1/(beta*dt))*uv(i,1) - (1/(2*beta))*ua(i,1) ;
        fdtrial =  statedeterEPP(k,fy,fs(i,1),ud(i,1),udtrial) ;
        punb = p(i+1) - ( m*(ua(i,1)+ dua) + c*(uv(i,1) + duv) + fdtrial ) ;
        
    end
    
    ud(i+1,1) = ud(i,1) + dud ;
    uv(i+1,1) = uv(i,1) + duv ;
    ua(i+1,1) = ua(i,1) + dua ;
    fs(i+1) = statedeterEPP(k,fy,fs(i,1),ud(i,1),ud(i+1,1)) ;
    
end














end

