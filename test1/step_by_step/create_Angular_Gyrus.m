clear
clc
cd('/imaging/rf02/TypLexMEG/fsaverage')
load('./MNI_mat_FG2.mat')
% unver_left=unver((unver>=0 & unver<10242));
MNI_AG=double(MNI_mat_FG2);
cd('/imaging/rf02/TypLexMEG')
[x,y,z]=textread('fsaverage_rh_MNIxyz.txt','%f %f %f', 'headerlines',1);
MNIrh=[x,y,z];%[x(1:10242),y(1:10242),z(1:10242)];
[x,y,z]=textread('fsaverage_lh_MNIxyz.txt','%f %f %f', 'headerlines',1);
MNIlh=[x,y,z];%[-x(1:10242),y(1:10242),z(1:10242)];
MNI_both=[MNIlh;MNIrh];
left_verts=[];
right_verts=[];
min_dist=[];
for ii=1:size(MNI_AG,1)%double(unver_left)
    ii
    this_vertex=repmat(MNI_AG(ii,:),size(MNI_both,1),1);
    this_dist=sqrt(sum((MNI_both-this_vertex).^2,2));
    [min_dist(ii),w_min_dist]=min(this_dist);
    %     while sum(unver==w_min_dist+10241)==1
    %         this_dist(w_min_dist)=20000;
    %         [min_dist,w_min_dist]=min(this_dist);
    %     end
    if (MNI_AG(ii,1)<0 && w_min_dist<=size(MNIlh,1) && min_dist(ii)<=13)
        left_verts=[left_verts,w_min_dist];
    elseif (MNI_AG(ii,1)>0 && w_min_dist>size(MNIlh,1) && min_dist(ii)<=13)
        right_verts=[right_verts,w_min_dist-size(MNIlh,1)];
    end
end
left_verts_u=unique(left_verts);
right_verts_u=unique(right_verts);
% unver=[-100,unver];
% save('bothhem_winner_vertices.mat','unver')
out_path='/imaging/rf02/TypLexMEG/fsaverage/FG2_left.mat';
save(out_path,'left_verts_u')
out_path='/imaging/rf02/TypLexMEG/fsaverage/FG2_right.mat';
save(out_path,'right_verts_u')