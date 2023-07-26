% Script is adapted on the work
% Altube P, Bech J, Argemí O, Rigo T, Pineda N, Collis S, Helmus J, (2017) 
% "Correction of dual-PRF Doppler velocity outliers in the presence of aliasing," 
% J. Atmos. Oceanic Technol., vol. 34, no. 7, pp. 1529-1543, doi: 10.1175/JTECH-D-16-0065.1
% It was adapted to Doppler method for industrial pipe flow

function [v1dc,v2dc,Betah,Betal,outlier1,outlier2] = OutlierCorrection(M,v1d,v2d,mratio,nratio,va1,va2)

if M==0 % no correction output=input!
    v1dc=v1d;
    v2dc=v2d;
    return %exit
end
%% Error Detection

vae=nratio*va2;

% Step 1
alpha1=v1d*pi/vae;
alpha2=v2d*pi/vae;

% Step 2
alpha1l=alpha1*mratio;
alpha2l=alpha2*nratio;

% Step 3
passo=M;central_gate=floor(passo/2)+1;
meio_passo=floor(passo/2);
Betah=zeros(size(v1d));
Betal=zeros(size(v2d));

alpha1ltemp= alpha1l(2:size(v1d,1)-1,:); %crop walls
alpha2ltemp= alpha2l(2:size(v1d,1)-1,:);

alpha1ltemp=[alpha1ltemp(1:central_gate-1,:); alpha1ltemp ;alpha1ltemp(end:-1:end-central_gate-1,:)]; %append
alpha2ltemp=[alpha2ltemp(1:central_gate-1,:); alpha2ltemp ;alpha2ltemp(end:-1:end-central_gate-1,:)];
indices=1:passo;

clear Betah;
clear Betal;
for j=1:size(v1d,2) % sweep emission by emission (time direction)
indices_sc=[indices(1:meio_passo) indices(central_gate+1:end)];m=2;
    for i=central_gate:size(v1d,1)-1+central_gate-1-1 
          Betah(m,j)=atan2(sum(sin(alpha1ltemp(indices_sc,j)))/passo,sum(cos(alpha1ltemp(indices_sc,j)))/passo );
          Betal(m,j)=atan2(sum(sin(alpha2ltemp(indices_sc,j)))/passo,sum(cos(alpha2ltemp(indices_sc,j)))/passo );
          indices_sc=indices_sc+1;
          m=m+1;
    end
end
Betal(m,:)=0;
Betah(m,:)=0;
% imagesc(Betah),c1=colorbar;

% Step 4
alpharef=Betal-Betah;
alpharef(alpharef<0)=alpharef(alpharef<0)+2*pi;

% Step 5
vout1=((alpha1-alpharef))*vae/pi;
vout2=((alpha2-alpharef))*vae/pi;

outlier1=abs(vout1)>(va1);
outlier2=abs(vout2)>(va2);
% for debug purposes
% figure; % imagesc(outlier1),c1=colorbar;title('outlier1')
% figure; imagesc(outlier2),c1=colorbar;title('outlier2')
% figure;imagesc(v1d),c1=colorbar;title('v1d')
% figure;imagesc(v2d),c1=colorbar;title('v2d')
% figure;imagesc(vout1),c1=colorbar;title('vout1')
% figure;imagesc(vout2),c1=colorbar;title('vout2')
% figure;imagesc(Betal),c1=colorbar;title('Betal')
% figure;imagesc(Betah),c1=colorbar;title('Betah')
%% Error correction 

v1dtemp= v1d(2:size(v1d,1)-1,:); %crop walls
v2dtemp= v2d(2:size(v2d,1)-1,:);
outlier1temp=outlier1(2:size(v2d,1)-1,:);
outlier2temp=outlier2(2:size(v2d,1)-1,:);

v1dtemp=[v1dtemp(1:central_gate-1,:); v1dtemp ;v1dtemp(end:-1:end-central_gate-1,:)]; %append
v2dtemp=[v2dtemp(1:central_gate-1,:); v2dtemp ;v2dtemp(end:-1:end-central_gate-1,:)];

outlier1temp=[zeros(length(1:central_gate-1),size(outlier1temp,2)); outlier1temp ;zeros(size(outlier1temp(end:-1:end-central_gate-1,:),1),size(outlier1temp,2))]; %append
outlier2temp=[zeros(length(1:central_gate-1),size(outlier1temp,2)); outlier2temp ;zeros(size(outlier1temp(end:-1:end-central_gate-1,:),1),size(outlier1temp,2))];
outlier1tempmask=outlier1temp<1;
outlier2tempmask=outlier2temp<1;

indices=1:passo;
clear localmeanh;
clear localmeanl;
for j=1:size(v1d,2)
    indices_sc=[indices(1:meio_passo) indices(central_gate+1:end)];m=2;
    for i=central_gate:size(v1d,1)-1+central_gate-1-1
        localmeanh(m,j)=sum(v1dtemp(indices_sc,j).*outlier1tempmask(indices_sc,j))/sum(outlier1tempmask(indices_sc,j));
        localmeanl(m,j)=sum(v2dtemp(indices_sc,j).*outlier2tempmask(indices_sc,j))/sum(outlier2tempmask(indices_sc,j));
        indices_sc=indices_sc+1;
        m=m+1;
    end
end
localmeanh(m,:)=0;
localmeanl(m,:)=0;

i=1;
clear vcorr1;
clear vcorr2;
clear vdif1;
clear vdif2;

for mc=-mratio:1:mratio
    vcorr1(:,:,i)=(outlier1.*v1d)-2*mc*va1;
    vcorr2(:,:,i)=(outlier2.*v2d)-2*mc*va2;
    vdif1(:,:,i)=outlier1.*(localmeanh)-outlier1.*(vcorr1(:,:,i));
    vdif2(:,:,i)=outlier2.*localmeanl-outlier2.*vcorr2(:,:,i);
    i=i+1;
end
mc=-mratio:1:mratio;
v1dc=v1d;
v2dc=v2d;
for i=2:size(vdif1,1)-1%crop walls
    for j=1:size(vdif1,2)
        itemp=(i-1)+central_gate-1;
        if (outlier1temp(itemp,j)==1) && (sum(outlier1temp(itemp-meio_passo:itemp+meio_passo,j))<(passo-1))
            [~,ind]=min(abs(vdif1(i,j,:)));
            v1dc(i,j)=vcorr1(i,j,ind);%+2*mc(ind)*va1;
        end
        if (outlier2temp(itemp,j)==1) && (sum(outlier2temp(itemp-meio_passo:itemp+meio_passo,j))<(passo-1))
             [~,ind]=min(abs(vdif2(i,j,:)));
            v2dc(i,j)=vcorr2(i,j,ind);
        end
    end
end
% figure;imagesc(v1dc),c1=colorbar;title('V1dc')
% figure;imagesc(v2dc),c1=colorbar;title('V2dc')
end