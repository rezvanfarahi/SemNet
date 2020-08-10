clear
clc
cd('/imaging/rf02/TypLexMEG')
load('./lefthem_winner_vertices.mat')
unver_left=unver((unver>=0 & unver<10242));
unver=double(unver_left);
[x,y,z]=textread('fsaverage_rh_MNIxyz.txt','%f %f %f', 'headerlines',1);
MNIrh=[x(1:10242),y(1:10242),z(1:10242)];
[x,y,z]=textread('fsaverage_lh_MNIxyz.txt','%f %f %f', 'headerlines',1);
MNIlh=[-x(1:10242),y(1:10242),z(1:10242)];
for ii=double(unver_left)
    ii
    this_vertex=repmat(MNIlh(ii+1,:),size(MNIrh,1),1);
    this_dist=sqrt(sum((MNIrh-this_vertex).^2,2));
    [min_dist,w_min_dist]=min(this_dist);
    while sum(unver==w_min_dist+10241)==1
        this_dist(w_min_dist)=20000;
        [min_dist,w_min_dist]=min(this_dist);
    end
    unver=[unver,w_min_dist+10241];
end
% unver=[-100,unver];
save('bothhem_winner_vertices.mat','unver')
    