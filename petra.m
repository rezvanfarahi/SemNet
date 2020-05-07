clear
clc
close all
fname='/home/rf02/Documents/petra/RAexercise.csv';
M=csvread(fname,1,0);
Mt=M(:,4:end);
Mt(Mt==0)=NaN;
M(:,4:end)=Mt;
M(159,3)=NaN;
Mr=zeros(size(M,1),5,8);
for ii=1:8
    Mr(:,1:4,ii)=M(:,1:4);
    Mr(:,5,ii)=M(:,4+ii);
end
p=[];
stats={};
pall=zeros(6,8);
pall_corr=zeros(6,8);
results1=[];
for ii=1:8
    g1=squeeze(Mr(:,2,ii));
    g2=squeeze(Mr(:,3,ii));
    g3=squeeze(Mr(:,4,ii));
    [p(:,ii),tb,stats{ii}]=anovan(squeeze(Mr(:,5,ii)),{g1 g2 g3},'display','off');%,'model','full'
    figure,
    [results,means] = multcompare(stats{ii},'CType','lsd');
    results1=[results1,results(:,[4,3,5])];
    pall(:,ii)=results(:,6);
    
end;
close all
for jj=1:6
    pall_corr(jj,:)=mafdr(pall(jj,:));
end

% pall_corr(pall_corr>=0.05)=1;
% pall_corr(pall_corr<1)=2;
figure, imagesc(pall_corr), colormap(hot)
figure, 
for ii=1:8
subplot(2,4,ii),boxplot(squeeze(Mr(:,5,ii)),squeeze(Mr(:,2,ii)))
end



C1=zeros(8,8,4);
kv=[0,1,2,3];
for kk=1:4
    Mr11=M(M(:,2)==kv(kk),5:12);
for ii=1:8
    for jj=1:8
        if ii>jj
        v1=Mr11(:,ii); nv1=isnan(v1);v1(nv1)=[];
        v2=Mr11(:,jj); v2(nv1)=[];nv2=isnan(v2);v2(nv2)=[];v1(nv2)=[];
        [bootstat,bootsam] = bootstrp(1000,@corr,v1,v2);
        bss=sign(sort(bootstat));
        if sum(bss==1)>=950 || sum(bss==-1)>=950
                    c=corrcoef([v1,v2]);
        C1(ii,jj,kk)=c(1,2);
        else
            C1(ii,jj,kk)=0;
        end
%         c=corrcoef([v1,v2]);
%         C1(ii,jj,kk)=c(1,2);
    end
end
end
end
figure,
kk=0;
for ii=1:4
    for jj=1:4
        kk=kk+1;
        if jj>ii
        subplot(4,4,kk), imagesc(C1(:,:,ii)+C1(:,:,ii)'-C1(:,:,jj)-C1(:,:,jj)',[-0.5,0.5]), colormap(hot)
        end
    end
end
figure,imagesc(C1(:,:,1)+C1(:,:,1)',[-0.5,0.5]), colormap(hot)
figure,imagesc(C1(:,:,2)+C1(:,:,2)',[-0.5,0.5]), colormap(hot)
figure,imagesc(C1(:,:,3)+C1(:,:,3)',[-0.5,0.5]), colormap(hot)
figure,imagesc(C1(:,:,4)+C1(:,:,4)',[-0.5,0.5]), colormap(hot)


% C1=corr(Mr1(:,2:end));
Mr1=M(M(:,2)==0,[2,5:12]);
Mr2=M(M(:,2)==1,[2,5:12]);
Mr3=M(M(:,2)==2,[2,5:12]);
Mr4=M(M(:,2)==3,[2,5:12]);
Mr1(sum(isnan(Mr1),2)>0,:)=[];
Mr2(sum(isnan(Mr2),2)>0,:)=[];
Mr3(sum(isnan(Mr3),2)>0,:)=[];
Mr4(sum(isnan(Mr4),2)>0,:)=[];


Indices1 = crossvalind('Kfold', size(Mr1,1), 5);
Indices2 = crossvalind('Kfold', size(Mr2,1), 5);
Indices3 = crossvalind('Kfold', size(Mr3,1), 5);
Indices4 = crossvalind('Kfold', size(Mr4,1), 5);
pp=pall_corr;
pp(pall_corr>=0.05)=1;
pp(pall_corr<0.05)=0;
pp=1-pp;
pp=[zeros(6,1),pp];
acc1=[];
for ii=1:5
    cms=pp(1,:)>0;
 svmStruct = svmtrain([Mr1(Indices1~=ii,cms);Mr2(Indices2~=ii,cms)],[Mr1(Indices1~=ii,1);Mr2(Indices2~=ii,1)],'ShowPlot',false);%;Mr3(Indices3~=ii,1);Mr4(Indices4~=ii,1)%;Mr3(Indices3~=ii,2:end);Mr4(Indices4~=ii,2:end)]
Group = svmclassify(svmStruct,[Mr1(Indices1==ii,cms);Mr2(Indices2==ii,cms)]);
refg=[zeros(sum(Indices1==ii),1);ones(sum(Indices2==ii),1)];
acc1(ii)=sum(abs(Group-refg)==0)/length(Group);
end

%%%
acc2=[];
for ii=1:5
    cms=pp(2,:)>0;
 svmStruct = svmtrain([Mr1(Indices1~=ii,cms);Mr3(Indices3~=ii,cms)],[Mr1(Indices1~=ii,1);Mr3(Indices3~=ii,1)],'ShowPlot',false);%;Mr3(Indices3~=ii,1);Mr4(Indices4~=ii,1)%;Mr3(Indices3~=ii,2:end);Mr4(Indices4~=ii,2:end)]
Group = svmclassify(svmStruct,[Mr1(Indices1==ii,cms);Mr3(Indices3==ii,cms)]);
refg=[zeros(sum(Indices1==ii),1);2*ones(sum(Indices3==ii),1)];
acc2(ii)=sum(abs(Group-refg)==0)/length(Group);
end


%%%
acc3=[];
for ii=1:5
    cms=pp(3,:)>0;
 svmStruct = svmtrain([Mr1(Indices1~=ii,cms);Mr4(Indices4~=ii,cms)],[Mr1(Indices1~=ii,1);Mr4(Indices4~=ii,1)],'ShowPlot',false);%;Mr3(Indices3~=ii,1);Mr4(Indices4~=ii,1)%;Mr3(Indices3~=ii,2:end);Mr4(Indices4~=ii,2:end)]
Group = svmclassify(svmStruct,[Mr1(Indices1==ii,cms);Mr4(Indices4==ii,cms)]);
refg=[zeros(sum(Indices1==ii),1);3*ones(sum(Indices4==ii),1)];
acc3(ii)=sum(abs(Group-refg)==0)/length(Group);
end

%%%
acc4=[];
for ii=1:5
    cms=pp(4,:)>0;
 svmStruct = svmtrain([Mr2(Indices2~=ii,cms);Mr3(Indices3~=ii,cms)],[Mr2(Indices2~=ii,1);Mr3(Indices3~=ii,1)],'ShowPlot',false);%;Mr3(Indices3~=ii,1);Mr4(Indices3~=ii,1)%;Mr3(Indices3~=ii,2:end);Mr4(Indices3~=ii,2:end)]
Group = svmclassify(svmStruct,[Mr2(Indices2==ii,cms);Mr3(Indices3==ii,cms)]);
refg=[ones(sum(Indices2==ii),1);2*ones(sum(Indices3==ii),1)];
acc4(ii)=sum(abs(Group-refg)==0)/length(Group);
end

%%%
acc5=[];
for ii=1:5
    cms=pp(5,:)>0;
 svmStruct = svmtrain([Mr2(Indices2~=ii,cms);Mr4(Indices4~=ii,cms)],[Mr2(Indices2~=ii,1);Mr4(Indices4~=ii,1)],'ShowPlot',false);%;Mr3(Indices3~=ii,1);Mr4(Indices4~=ii,1)%;Mr3(Indices3~=ii,2:end);Mr4(Indices4~=ii,2:end)]
Group = svmclassify(svmStruct,[Mr2(Indices2==ii,cms);Mr4(Indices4==ii,cms)]);
refg=[ones(sum(Indices2==ii),1);3*ones(sum(Indices4==ii),1)];
acc5(ii)=sum(abs(Group-refg)==0)/length(Group);
end

%%%
acc6=[];
for ii=1:5
    cms=pp(6,:)>0;
 svmStruct = svmtrain([Mr3(Indices3~=ii,cms);Mr4(Indices4~=ii,cms)],[Mr3(Indices3~=ii,1);Mr4(Indices4~=ii,1)],'ShowPlot',false);%;Mr3(Indices3~=ii,1);Mr4(Indices4~=ii,1)%;Mr3(Indices3~=ii,2:end);Mr4(Indices4~=ii,2:end)]
Group = svmclassify(svmStruct,[Mr3(Indices3==ii,cms);Mr4(Indices4==ii,cms)]);
refg=[2*ones(sum(Indices3==ii),1);3*ones(sum(Indices4==ii),1)];
acc6(ii)=sum(abs(Group-refg)==0)/length(Group);
end

