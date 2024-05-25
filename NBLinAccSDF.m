function [ ud , t] = NBLinAccSDF(uga,dt,wn,zhi)

% SI  units followed

m = 1 ; % 1 kg assumed, because this code is primarily intended for SDoFs subjected to ground motion
k = wn^2 ;
c = 2*m*wn*zhi;
n = length(uga) ;
t = 0:1:(n-1) ;
t = t' ;
t = t.*dt ; %Time axis
p = zeros(n,1) ;
p = -m.*9.81.*uga ;

ud = zeros(n,1) ;
uv = zeros(n,1) ;
ua = zeros(n,1) ;

alpha = 0.5 ;
beta = 1/6 ;

%Initial caculations (Change if necessary)

ud(1,1) = 0 ; % Assumed
uv(1,1) = 0 ; % Assumed

ua(1,1) = (p(1,1) - c*uv(1,1) - k*ud(1,1))/m ;
khat = k + (alpha/(beta*dt))*c + (1/(beta*dt^2))*m ;
a = (1/(beta*dt))*m + (alpha/beta)*c ;
b = (1/(2*beta))*m + dt*((alpha/(2*beta)) -1 )*c ;

for i = 1:1:(n-1)
    
    dphat = (p(i+1)-p(i,1)) + a*uv(i,1) + b*ua(i,1) ;
    dud = dphat/khat ;
    duv = (alpha/(beta*dt))*dud - (alpha/beta)*uv(i,1) + dt*(1- (alpha/(2*beta)) )*ua(i,1) ;
    dua = (1/(beta*dt^2))*dud - (1/(beta*dt))*uv(i,1) - (1/(2*beta))*ua(i,1) ;
    ud(i+1) = ud(i,1) + dud ;
    uv(i+1) = uv(i,1) + duv ;
    ua(i+1) = ua(i,1) + dua ;
    
    
    
end








end

