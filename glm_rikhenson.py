#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 23:52:52 2020

@author: rezvanh
"""

import numpy as np
from np.linalg import pinv as pinv
from np.linalg.matrix_rank import matrix_rank as matrix_rank
def glm_rik(y,X,c,pflag=None):
    if pflag == None:
        pflag = 0
    
    B=pinv(X)@y
    
    Y = X@B
    
    r = y - Y
    
    cE = np.eye(len(r))
    Bcov = (pinv(X)@cE)@pinv(X).T
    
    if c.shape[1] > c.shape[0]:
        print('Transposing c!')
        c = c.T;
        
    l = c.shape[0]
    if l < X.shape[1]:
        c = np.vstack((c,np.zeros((X.shape[1]-l,c.shape[1]))))
        print('Padding c with zeros!')
    
    df = len(y) - matrix_rank(X)
    
    if c.shape[1]==1:
        s = (r.T@r) / df
        t = (c.T@B) / np.sqrt(s@c.T@pinv(X.T@X)@c)
        
    



if size(c,2)==1
  s = r'*r / df;
  t = c'*B / sqrt(s*c'*pinv(X'*X)*c);
  p = t2p(t,df);
  F = t.^2; 
  cR2 = F2R(F,1,df);
  R2 = 1 - (r'*r) / (y'*y); 
  aR2 = 1 - (r'*r/df(end)) / (y'*y/(length(y)-1));
  if pflag~=-1
     fprintf('T(%d)=%f, p(one-tailed)=%f  (R2=%3.2f; overall R2=%3.2f (adjusted R2=%3.2f))\n',df,t,p,cR2,R2,aR2);
  end
else
  c_0 = eye(size(X,2)) - c*pinv(c);
  X_0 = X*c_0;
  R   = eye(size(X,1)) - X*pinv(X);
  R_0 = eye(size(X,1)) - X_0*pinv(X_0);
  M = R_0 - R;
  df = [rank(X)-rank(X_0) size(X,1)-rank(X)];
  F  = ((B'*X'*M*X*B)/df(1)) / ((y'*R*y)/df(2)); 
  p  = F2p(F,df(1),df(2));
  t = [];
  cR2 = F2R(F,df(1),df(2));
  R2 = 1 - (r'*r) / (y'*y); 
  aR2 = 1 - (r'*r/df(end)) / (y'*y/(length(y)-1));
  if pflag~=-1
      fprintf('F(%d,%d)=%f, p(two-tailed)=%f  (R2=%3.2f; overall R2=%3.2f (adjusted R2=%3.2f))\n',df(1),df(2),F,p,cR2,R2,aR2); 
  end
end

%Below doesn't work if constant not included
%R2 = 1 - (r'*r) / (y'*y);   %R2 = (Y'*Y)/(y'*y) equivalent
%disp(sprintf('Overall R2 = %f',R2))
%aR2 = 1 - (r'*r/df(end)) / (y'*y/(length(y)-1));
%disp(sprintf('Overall Adjusted R2 (assumes a constant in model) = %f',R2))

% Hacky way to calculate unique contribution of contrast 
c_0 = eye(size(X,2)) - c*pinv(c);
X_0 = X*c_0;
y_0 = X_0*(pinv(X_0)*y);
r_0 = y - y_0;
iR2 = R2 - (1 - (r_0'*r_0) / (y'*y));


if pflag==1
    figure(pflag), clf
    Ne = size(c,2);
    for e = 1:Ne
        subplot(1,Ne,e), hold on
        Xc = X*c(:,e);
        Yc = Xc*(c(:,e)'*B);
        yc = Yc+r;
        plot(Xc,yc,'r.')
        plot(Xc,Yc,'b-')
    end
end


return

%%%%%%%%%%%%%%%%%%
function p=t2p(t,df);

% one-tailed p-value

try
    p=tcdf(t,df);
catch
    try 
        p=spm_Tcdf(t,df);
    catch
        error('Need Matlab Stats toolbox (tcdf) or SPM on path')
    end
end

f=find(p>0.5);
p(f)=1-p(f);

return

%%%%%%%%%%%%%%%%%%%%

function p=F2p(F,df1,df2);
% one-tailed p-value

try
    p=1-fcdf(F,df1,df2);
catch
    try 
        p=1-spm_Fcdf(F,df1,df2);
    catch
        error('Need Matlab Stats toolbox (fcdf) or SPM on path')
    end
end

return

%%%%%%%%%%%%%%%%%%%%
function R2=F2R(F,df1,df2);

R2=df1*F/(df2+df1*F);

return
