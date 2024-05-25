clear all
close all
clc

U = importdata('ElCentro.txt') ;
t = U(:,1) ;
dtor = 0.02 ;
ugor = U(:,2) ;
oldsize = length(ugor) ;
dtnew = 0.005 ;
div = dtor/dtnew ;
newsize = oldsize + (oldsize-1)*(div-1) ;
tnew = [0:dtnew:(t(oldsize))]' ;
ugnew = zeros(newsize,1) ;

i = div + 1 ;

for j = 1:1:(oldsize-1)
    
    slope =  ( ugor(j+1,1) - ugor(j,1) )/dtor ;
    ugnew(i-4,1) = ugor(j,1) ;
    ugnew(i-3,1) = ugnew(i-4,1) + slope*dtnew ;
    ugnew(i-2,1) = ugnew(i-4,1) + 2*slope*dtnew ;
    ugnew(i-1,1) = ugnew(i-4,1) + 3*slope*dtnew ;
    ugnew(i,1) = ugor(j+1,1) ;
    i = i + div ;
    
end

plot(tnew,ugnew)
Unew = [tnew ugnew] ;
save('UpElcentro.txt','Unew','-ascii')
