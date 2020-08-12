
exp1  = load('/imaging/rf02/Semnet/semnet4semloc/dcm/SemLoc_SD_80wins.mat');
exp2a = load('/imaging/rf02/Semnet/semnet4semloc/dcm/SemNet_SD_80wins.mat');
exp2b = load('/imaging/rf02/Semnet/semnet4semloc/dcm/SemNet_LD_80wins.mat');

Nr = size(exp1.Datafm_wins,1)
Nt = size(exp1.Datafm_wins,2)
xaxis=[50:5:50+79*5];
Ns(1) = size(exp1.Datafm_wins,4);
Ns(2) = size(exp2a.Datafm_wins,4)

X = [ones(Ns(1),1) zeros(Ns(1),2)];
X = [X; zeros(2*Ns(2),1) kron(eye(2),ones(Ns(2),1))]
X = [X [zeros(Ns(1),Ns(2)); eye(Ns(2)); eye(Ns(2))]];
figure,imagesc(X)

Fcom = []; Pcom = []; Fdif = []; Pdif = [];
for r = 1:Nr
    for t = 1:Nt
        y =     squeeze(exp1.Datafm_wins(r,t,:,:))'  * [1 -1]'; % contrast of conditions
        y = [y; squeeze(exp2a.Datafm_wins(r,t,:,:))' * [1 -1]']; % contrast of conditions
        y = [y; squeeze(exp2b.Datafm_wins(r,t,:,:))' * [1 -1]']; % contrast of conditions
        %[Tval,Fcom(r,t),Pcom(r,t),df,R2,cR2,B,res,aR2,iR2,Bcov] = glm(y,X,[1 1 1 zeros(1,Ns(2))]');
        [Tval,Fcom(r,t),Pcom(r,t),df,R2,cR2,B,res,aR2,iR2,Bcov] = glm(y,X,[1 1 1 ones(1,Ns(2))*2/Ns(2)]'); % http://www.sbirc.ed.ac.uk/cyril/download/Contrast_Weighting_Glascher_Gitelman_2008.pdf
        [Tval,Fdif(r,t),Pdif(r,t),df,R2,cR2,B,res,aR2,iR2,Bcov] = glm(y,X,[detrend(eye(3),0) zeros(3,Ns(2))]');
    end
end

alpha = 0.05 % uncorrected
%alpha = alpha/(Nr*Nt) % extreme Bonferonni, ignoring that time windows correlated
%alpha = alpha/Nr % compromise!

figure,imagesc(Pcom),colorbar,xlabel('Twin'),ylabel('ROI'),set(gca,'XTick',[50:5:50+79*5],'YTick',[1:5]),colormap('gray'),caxis([0 alpha]),title('Common condition effect across exps')
figure,imagesc(Pdif),colorbar,xlabel('Twin'),ylabel('ROI'),set(gca,'XTick',[50:5:50+79*5],'YTick',[1:5]),colormap('gray'),caxis([0 alpha]),title('Interaction between exp and condition')
