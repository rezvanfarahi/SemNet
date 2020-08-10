for f=1:4;
    for s=1:17
    for t=1:240
        mat(s,:,:,f,t)=squeeze(plv_matc(s,:,:,f,t))+squeeze(plv_matc(s,:,:,f,t))';
    end
    end
end

 for s=1:17
     smat=squeeze(mat(s,:,:,:,:));
     smat(smat>prctile(smat(:),70))=1;
     smat(smat<=prctile(smat(:),70))=0;
     mat2(s,:,:,:,:)=smat;
 end
     
  

figure; imagesc(squeeze(mean(mean(mat2(:,:,:,1,:),1),2))); colorbar