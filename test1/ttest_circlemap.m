load('/imaging/rf02/TypLexMEG/connmat_abstract.mat')
load('/imaging/rf02/TypLexMEG/connmat_concrete.mat')
datheta=reshape(Dmat_abs(:,:,1,:),size(Dmat_abs,1),size(Dmat_abs,2),size(Dmat_abs,4));
dctheta=reshape(Dmat_cncrt(:,:,1,:),size(Dmat_cncrt,1),size(Dmat_abs,2),size(Dmat_abs,4));
for ii=1:size(datheta,1)
for jj=1:size(datheta,2)
if jj>ii
datheta(ii,jj,:)=datheta(jj,ii,:);
dctheta(ii,jj,:)=dctheta(jj,ii,:);
end
end
end

% h=[];p=[];
% for ii=1:size(datheta,1); for jj=1:size(datheta,2)
% [h(ii,jj),p(ii,jj)]=ttest(datheta(ii,jj,:),dctheta(ii,jj,:));
% end
% end

dabst=datheta(datheta~=0);
dabst=reshape(sum(reshape(dabst,71,72,17),1),72,17);
dcrtt=dctheta(dctheta~=0);
dcrtt=reshape(sum(reshape(dcrtt,71,72,17),1),72,17);

h=[];p=[];
for ii=1:size(dabst,1)
[h(ii),p(ii)]=ttest(dabst(ii,:),dcrtt(ii,:));
end



