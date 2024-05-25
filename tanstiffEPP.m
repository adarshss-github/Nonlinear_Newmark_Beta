function [ kt ] = tanstiffEPP(k,fy,fi)

if abs(fi)>fy
    
    kt = 0 ;
    
else
    
    kt = k ;
    
end


end

