function [ ff ] = statedeterEPP(k,fy,fc,uc,uf)

du = uf - uc ;
ftrial = fc + du*k ;

if abs(ftrial) > fy
    
    ff = sign(ftrial)*fy ;
    
    
else
    
    ff = ftrial ;
    
    
end


end

