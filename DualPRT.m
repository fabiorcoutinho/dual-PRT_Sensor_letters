function [v1d,v2d,v1,v2,va1,va2,W] = DualPRT(iq,c,fc,T1,T2,Ns,Nc,m,n,F)
% Modification of the autocorrelation method or phase shift estimator
% to estimate speed using the staggered trigger technique or
% dual PRT or staggered PRT. At the end, implement the defined rules
% by Torres et al. in the article "Design, Implementation, and Demonstration
%of a Staggered PRT Algorithm for the WSR-88D", Journal of Atmospheric and Oceanic
% Technology, 2004, volume 21, 1389-1399, to perform the dealiasing of the
% estimated velocity.
% 

% input:
% iq: ultrasound data matriz (IQ demodulated) lines= ultrasound emissions
% c soundspeed
% fc : transducer central frequency
% T1:  pulse% repetition time 1 ou PRT1) - T1  < e T2
% T2:  pulse repetition time 2 ou PRT2)
% Ns: number of spatial samples -range gate, 
% Nc: number os emissions
% m e n: staggered PRT ratio m/n

% Outuput
% v1d: dealiased spatiotemporal velocity map regarding T1 
% v2d: dealiased spatiotemporal velocity map regarding T2 
% v1: aliased spatiotemporal velocity map regarding T1 
% v2: aliased spatiotemporal velocity map regarding T2 
% va1: Nyquist's maximum velocity for PRT1 or T1
% va2: Nyquist's maximum velocity for PRT2 or T2
% W: spatiotemporal velocity energy


%adapts the num of spatial samples to be windowed  with Ns
% by adding zero values ​​at the end of the array (zero pads)
pad_n=Ns-rem(size(iq,1)-Ns,Ns);
pad=zeros(pad_n,size(iq,2));
iq=[iq ; pad]; % acrescenta zeros na profundidade no fim do vetor

% calculates the number of spatial channels already with the pad of zeros at the end 
nchannels=(size(iq,1)-Ns)/(Ns);
% calculates the number of time channels
nchannelst=floor((size(iq,2)-Nc)/(Nc))+0;

%initializes the algorithm output result matrices with zeros
v1d=zeros(nchannels,nchannelst); % spatiotemporal velocity map (dealiased)  T1
v2d=zeros(nchannels,nchannelst); % spatiotemporal velocity map (dealiased) radial T2
v12=zeros(nchannels,nchannelst); % spatiotemporal velocity map - v1-v2
v1=zeros(nchannels,nchannelst); % spatiotemporal velocity map aliased T1
v2=zeros(nchannels,nchannelst); % spatiotemporal velocity map aliased T2
W=zeros(nchannels,nchannelst); % spatiotemporal energy

for i=1:nchannelst % scans all time channels
    p_i=(i-1)*Nc+1;
    p_f=p_i+Nc;
    for j=1:nchannels % scans all spatial channels
        ps_i=(j-1)*Ns+1;
        ps_f=ps_i+Ns-1;
        data=iq(ps_i:ps_f,p_i:p_f);
        vdata=sum(data,1); % integrates all values ​​of a sample range
        if F ~=0 % F not zero indicates to filter data for eliminate stationary echoes
            vdatar=F*real(vdata)';
            vdatar=vdatar';%*F;%
            vdatai=F*imag(vdata)';
            vdatai=vdatai';%*F;%
            vdata=complex(vdatar,vdatai);
        end
        Nc1=length(vdata(2:2:end));
        Nc2=length(vdata(3:2:end));
        %  Calculate the autocorrelation of R(T1) and R(T2)
        auto1T1= (vdata(2:2:end)* vdata(1:2:end-1)') / (Nc1-1); % R(T1)
        auto1T2= (vdata(3:2:end)* vdata(2:2:end-1)') / (Nc2-1);% R(T2)
        %  Calculate the autocorrelation of RT1(0) and RT2(0)
        auto0T1=(sum(abs(vdata(2:2:end)).^2)+sum(abs(vdata(1:2:end)).^2)) /   (2*(Nc1-1));  % R(0)=W
        auto0T2=(sum(abs(vdata(3:2:end)).^2)+sum(abs(vdata(2:2:end-1)).^2)) / (2*(Nc2-1));  % R(0)=W
        W(j,i)=(auto0T1+auto0T2)/2;
                
        %  Calculate the phase angle from the autocorrelation of lag 1
        phi_est1 =  atan2(imag(auto1T1),real(auto1T1));
        phi_est2 =  atan2(imag(auto1T2),real(auto1T2));

        % calculate the axial velocity
        v1(j,i) = -c*(1/(T1))/(4*pi*(fc)) *  phi_est1;
        v2(j,i) = -c*(1/(T2))/(4*pi*(fc)) *  phi_est2;
        v12(j,i)= v1(j,i)-v2(j,i); % velocity difference
    end

end

%Dealiasing
% calculates maximum conventional velocity for each period (T1, T2)
va1=c/(4*T1*fc);
va2=c/(4*T2*fc);
va=m*va1*1;
dC=2*va/(m*n);

% constant values for Dealising rules
C0=0;
C1=C0+2*va2;
C2=C1-2*va1;
C3=C2+2*va2;
C4=C3-2*va1;
C5=C4+2*va2;
C6=C5-2*va1;


% Aapply dealising rules according to:
% M. Torres, Y. Dubel, D. S. Zrnic "Design, implementation, and demonstration
% of a staggered PRT algorithm for the WSR-88D," J. Atmos. Ocean. Technol, 
% vol. 21, no. 9, pp. 1389-1399, Sep. 2004, 
% DOI: 10.1175/1520-0426(2004)021<1389:DIADOA>2.0.CO;2

% m=0
rule=(v12>(C0-dC/2))&(v12<(C0+dC/2));
P=0;Q=0;
v1d(rule)=v1(rule)+2*P*va1;
v2d(rule)=v2(rule)+2*Q*va2;


if m>=1
    rule1p=(v12>(C1-dC/2))&(v12<(C1+dC/2));
    P=0;Q=1;
    v1d(rule1p)=v1(rule1p)+2*P*va1;
    v2d(rule1p)=v2(rule1p)+2*Q*va2;
    %%%
    rule1n=(v12<-(C1-dC/2))&(v12>-(C1+dC/2));
    P=0;Q=-1;
    v1d(rule1n)=v1(rule1n)+2*P*va1;
    v2d(rule1n)=v2(rule1n)+2*Q*va2;
end

if m>=2
    rule2p=(v12>(C2-dC/2))&(v12<(C2+dC/2));
    P=1;Q=1;
    v1d(rule2p)=v1(rule2p)+2*P*va1;
    v2d(rule2p)=v2(rule2p)+2*Q*va2;
    %%%%
    rule2n=(v12<-(C2-dC/2))&(v12>-(C2+dC/2));
    P=-1;Q=-1;
    v1d(rule2n)=v1(rule2n)+2*P*va1;
    v2d(rule2n)=v2(rule2n)+2*Q*va2;
end

if m>=3
    rule=(v12>(C3-dC/2))&(v12<(C3+dC/2));
    P=1;Q=2;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
    %%%
    rule=(v12<-(C3-dC/2))&(v12>-(C3+dC/2));
    P=-1;Q=-2;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
end

if m>=4
    rule=(v12>(C4-dC/2))&(v12<(C4+dC/2));
    P=2;Q=2;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
    %%%
    rule=(v12<-(C4-dC/2))&(v12>-(C4+dC/2));
    P=-2;Q=-2;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
end

if m>=5
    rule=(v12>(C5-dC/2))&(v12<(C5+dC/2));
    P=2;Q=3;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
    % % % % %
    rule=(v12<-(C5-dC/2))&(v12>-(C5+dC/2));
    P=-2;Q=-3;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
end

if m>=6
    rule=(v12>(C6-dC/2))&(v12<(C6+dC/2));
    P=3;Q=3;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
    %%%
    rule=(v12<-(C6-dC/2))&(v12>-(C6+dC/2));
    P=-3;Q=-3;
    v1d(rule)=v1(rule)+2*P*va1;
    v2d(rule)=v2(rule)+2*Q*va2;
end

end