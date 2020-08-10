clear
close all
t=-pi:0.01:pi;
S2=sin(t+pi/4);
S1=0.5*sin(t)+sin(2*t+pi/3);
C=corrcoef(S1,S2); C=C(1,2);
ii=0;jj=0;
for A=0:100
    ii=ii+1; jj=0;
    for B=0:100
        jj=jj+1;
        S1est=A*S1+B*S2; %estimated S1
        S2est=A*S2+B*S1;
        rat(ii,jj)=A/B;
        Rp=corrcoef(S1est,S2);R_uncorrected(ii,jj)=Rp(1,2); 
        Spnew=S1est-S1est.*S2;
        Rp=corrcoef(Spnew,S2);R_corrected(ii,jj)=Rp(1,2); 
    end
end
X=0:100;Y=0:100;
figure,pcolor(X,Y,R_uncorrected)
colorbar()
xlabel('B')
ylabel('A')
title(['correlation before leakage correction, Real Correlation: ', num2str(C)])

figure,pcolor(X,Y,R_corrected)
colorbar()
xlabel('B')
ylabel('A')
title(['correlation after leakage correction, Real Correlation: ', num2str(C)])

