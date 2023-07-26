% Generate the coefficients for stationary echoes filtering
% The filter implementation is based in the work of:
% S. M. Torres and D. S. Zrnic “Ground clutter cancelling with a regression
% filter” J. Atmos. Ocean. Technol, vol. 16, no. 10, pp. 1364-1372, Oct. 
% 1999, DOI: 10.1175/1520-0426(1999)016<1364:GCCWAR>2.0.CO;2
% 
%
% Input:
% Mf: length of the regressin filter (number of temporal samples)
% p: regression polynomial order 
% m e n: staggered PRT ratio T1/T2=m/n

% Output
% F: regression filter matrix 
% B: Basis matrix
% tm: time vector used for generating the basis matrix
function [F,B,tm] = PolyRegFilterLeg(Mf,p,m,n)

F=zeros(Mf,Mf+1);% initialize F based in Mf
B=zeros(p,Mf+1); % initialize B based in p and Mf
tm=zeros(1,Mf+1);% create a time vector
for i=1:Mf+1 % sweep from 1 to Mf+1 to create the time vector tm based on PRT ratio
    if i==1
        tm(i)=0;
    else
        if mod(i,2)==0 % even means T1 (to match with the data line drop pattern) 
            tm(i)=tm(i-1)+m;
        else         % odd meand T2
            tm(i)=tm(i-1)+n;
        end
    end
end
tm=2*tm/tm(end)-1; %normalizes tm to be between -1 to +1 due to lengedre function input criterion
for i=0:p % create the Basis matrix using lengendre polynomial function 
    temp=legendre(i,tm);
    B(i+1,:)=temp(1,:)/sqrt(temp(1,:)*temp(1,:)');% normalize B
end
F=eye(Mf+1)-B'*B;% compute equation (6) of reference
end



