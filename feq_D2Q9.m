function [feq]= feq_D2Q9(k,rho,u,v,w,e,cs2)
    t1=u.*u+v.*v;
    t2=u.*e(k,1)+v.*e(k,2);
    feq=rho.*w(k).*(1.0+t2/cs2+t2.*t2/(2*cs2*cs2)-t1/(2*cs2));
end