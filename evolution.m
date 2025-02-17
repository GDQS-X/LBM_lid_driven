for k = 1:1:9    
    f(k,:,:)= circshift(f(k,:,:),[0,round(e(k,1)/cc),round(e(k,2)/cc)]);
end