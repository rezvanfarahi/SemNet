load('/imaging/rf02/TypLexMEG/Xa_pca.mat')
load('/imaging/rf02/TypLexMEG/Xc_pca.mat')
expVc=[];
expVa=[];
Xs_pca=Xc_pca-Xa_pca;
for K=2
K    
[IDXc,Cc,sumdc,Dc]=kmeans(Xs_pca',K);
[IDXa,Ca,sumda,Da]=kmeans(Xa_pca',K);
[xxc,yyc]=min(Dc,[],2);
unexpVc=sum(xxc);

[xxa,yya]=min(Da,[],2);
unexpVa=sum(xxa);

pdc=pdist(Xc_pca');
totVc=sum(pdc.^2)/17;

pda=pdist(Xa_pca');
totVa=sum(pda.^2)/17;

expVc(K)=100*(totVc-unexpVc)/totVc;
expVa(K)=100*(totVa-unexpVa)/totVa;
end