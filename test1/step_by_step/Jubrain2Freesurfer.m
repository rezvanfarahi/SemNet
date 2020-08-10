clear
addpath('C:\Users\rf02\Documents\Rezvan\PhDproject\semloc\April16 report\Toolbox_22c\NIfTI_20140122')
nii = load_nii('C:\Users\rf02\Documents\Rezvan\PhDproject\semloc\April16 report\Toolbox_22c\Anatomy\PMaps\IPL_PGa.nii');
ni=nii.img;
ni(ni<0.5)=0;
nii.img=ni;
ff1=find(ni>0);
[aa,bb,cc]=ind2sub(size(ni),ff1);
coor_mat=[aa,bb,cc];
orig_mat=repmat(nii.hdr.hist.originator(1:3),size(coor_mat,1),1);
MNI_mat_PGa=coor_mat-orig_mat;
out_path='\\cbsu\data\Imaging\rf02\TypLexMEG\fsaverage\MNI_mat_PGa.mat';
save(out_path, 'MNI_mat_PGa')

nii = load_nii('C:\Users\rf02\Documents\Rezvan\PhDproject\semloc\April16 report\Toolbox_22c\Anatomy\PMaps\IPL_PGp.nii');
ni=nii.img;
ni(ni<0.5)=0;
nii.img=ni;
ff2=find(ni>0);
[aa,bb,cc]=ind2sub(size(ni),ff2);
coor_mat=[aa,bb,cc];
orig_mat=repmat(nii.hdr.hist.originator(1:3),size(coor_mat,1),1);
MNI_mat_PGp=coor_mat-orig_mat;
out_path='\\cbsu\data\Imaging\rf02\TypLexMEG\fsaverage\MNI_mat_PGp.mat';
save(out_path, 'MNI_mat_PGp')


nii = load_nii('C:\Users\rf02\Documents\Rezvan\PhDproject\semloc\April16 report\Toolbox_22c\Anatomy\PMaps\IPL_PF.nii');
ni=nii.img;
ni(ni<0.5)=0;
nii.img=ni;
ff3=find(ni>0);
[aa,bb,cc]=ind2sub(size(ni),ff3);
coor_mat=[aa,bb,cc];
orig_mat=repmat(nii.hdr.hist.originator(1:3),size(coor_mat,1),1);
MNI_mat_PF=coor_mat-orig_mat;
out_path='\\cbsu\data\Imaging\rf02\TypLexMEG\fsaverage\MNI_mat_PF.mat';
save(out_path, 'MNI_mat_PF')

ff=[ff1;ff2;ff3];
ffu=unique(ff);
[aa,bb,cc]=ind2sub(size(ni),ffu);
coor_mat=[aa,bb,cc];
orig_mat=repmat(nii.hdr.hist.originator(1:3),size(coor_mat,1),1);
MNI_mat_AGSMG=coor_mat-orig_mat;
out_path='\\cbsu\data\Imaging\rf02\TypLexMEG\fsaverage\MNI_mat_AGSMG.mat';
save(out_path, 'MNI_mat_AGSMG')

nii = load_nii('C:\Users\rf02\Documents\Rezvan\PhDproject\semloc\April16 report\Toolbox_22c\Anatomy\PMaps\Visual_FG2.nii');
ni=nii.img;
ni(ni<0.5)=0;
nii.img=ni;
ff4=find(ni>0);
ff4=unique(ff4);
[aa,bb,cc]=ind2sub(size(ni),ff4);
coor_mat=[aa,bb,cc];
orig_mat=repmat(nii.hdr.hist.originator(1:3),size(coor_mat,1),1);
MNI_mat_FG2=coor_mat-orig_mat;
out_path='\\cbsu\data\Imaging\rf02\TypLexMEG\fsaverage\MNI_mat_FG2.mat';
save(out_path, 'MNI_mat_FG2')