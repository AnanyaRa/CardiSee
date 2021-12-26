clear
%close all
% Objective of hrv1b: access Physionet data OR read local URochester THEW
% files and calculate various metrics
%% Links of various Physionet Databases
% https://physionet.org/physiobank/database/mitdb/
% https://physionet.org/physiobank/database/cudb/
% https://physionet.org/physiobank/database/nsr2db/
% https://www.physionet.org/physiobank/database/chf2db/
% https://www.physionet.org/physiobank/database/chfdb/
%%--
%%
load('/Users/ananyar/Documents/MATLAB/Add-Ons/Apps/HRVTool/code/data/HLT/10106_ann.mat');
Fs=200;

%Ann = rdann('nsrdb/19830','atr');
%Ann = rdann('nsr2db/nsr008','ecg');
%Fs=128;
%Ann = rdann('chf2db/chf229','ecg');
%Fs = 128;
%Ann = rdann('chfdb/chf002','ecg');
%Fs = 250;
%Ann = rdann('cudb/cu03','atr');
%Fs = 250;
%Ann = rdann('mitdb/234','atr');
%Fs = 360;
%Ann = rdann('svdb/828','atr');
%Fs = 128;
%load('/Users/ananyar/Documents/MATLAB/Add-Ons/Apps/HRVTool/code/data/100m.mat')
%Fs=360;
%Ann = rdann('afdb/00735','qrs');
%Fs = 250;
%Ann = rdann('ltafdb/110','qrs');
%Fs = 128;
%--
%2015 Physionet Challenge Data
%load('/Users/ananyar/Downloads/training/a103l.mat')
%[tm,signal,Fs,siginfo]= rdmat('/Users/ananyar/Downloads/training/a103l');
%sqrs('/Users/ananyar/Downloads/training/a103l',[],[],2,1000)
%Fs=250;

Ann=Ann/Fs;
Ann = Ann(0.1*3600:end);

RR = [NaN; diff(Ann)];
%RR = [0; diff(Ann)];

RR = HRV.RRfilter(RR,90);
%RR = RR(RR>=0.6 & RR <=1.2);

%       % Corresponding relative RR intervals:

rr = HRV.rrx(RR);
% length(rr(~isnan(rr)))    % number of intervals without NaN
sdnn_RR = HRV.SDNN(RR,0,1)*1000;
pnn50_RR = HRV.pNN50(RR,0,0)*100;
%sdsd_RR = HRV.SDSD(RR,0,1);
%sd1psd2=0.7071*(0.7071*sdnn_RR + sqrt(2*sdnn_RR^2-0.5*sdsd_RR^2));

% Out_boundary_calculation
% rr>0.05 returns logical (TRUE/FALSE) array of elements that meet condition
% rr(rr>0.05) returns elements ONLY with TRUE value
%val_out5p=(sum(rr(rr>0.05))+abs(sum(rr(rr<-0.05))))/(length(rr(rr>0.05))+length(rr(rr<-0.05)))*100

val_out5p=sum(abs(rr(abs(rr)>0.05)))/length(rr(abs(rr)>0.05))*100;
val_out10p=sum(abs(rr(abs(rr)>0.10)))/length(rr(abs(rr)>0.10))*100;
val_out20p=sum(abs(rr(abs(rr)>0.20)))/length(rr(abs(rr)>0.20))*100;

%pct_5p=length(rr(abs(rr)>0.05))/size(rr,1)*100;
%pct_10p=length(rr(abs(rr)>0.10))/size(rr,1)*100;
%pct_20p=length(rr(abs(rr)>0.20))/size(rr,1)*100;

% consecutive +20/-20 swings
rra20p = rr > 0.20;
rrd20p = -1*(rr < -0.20);
rrs1_20p=[rra20p+rrd20p; NaN];
rrs2_20p=[NaN; rra20p+rrd20p];
rrs12_20p=rrs1_20p.*rrs2_20p; % -1 determines consecutive +20/-20 swings
z3e20p = (sum(rra20p == 1) + sum(rrd20p == -1))/size(rr,1)*100;
cm_20p=sum(rrs12_20p == -1)/size(rr,1)*100;

% consecutive +10%/-10% swings
rra10p = rr > 0.10;
rrd10p = -1*(rr < -0.10);
rrs10p=rra10p+rrd10p;
rrs1_10p=[rra10p+rrd10p; NaN];
rrs2_10p=[NaN; rra10p+rrd10p];
rrs12_10p=rrs1_10p.*rrs2_10p; %-1 determines consecutive +10/-10 swings
%cm_10p = (sum(rra10p == 1) + sum(rrd10p == -1))/size(rr,1)*100;
cm_10p=sum(rrs12_10p == -1)/size(rr,1)*100;

% consecutive +5%/-5% swings
rra5p = rr > 0.05;
rrd5p = -1*(rr < -0.05);
rrs1_5p=[rra5p+rrd5p; NaN];
rrs2_5p=[NaN; rra5p+rrd5p];
rrs12_5p=rrs1_5p.*rrs2_5p; %-1 determines consecutive +5/-5 swings
%cm_5p = (sum(rra5p == 1) + sum(rrd5p == -1))/size(rr,1)*100;
cm_5p=sum(rrs12_5p == -1)/size(rr,1)*100;

% consecutive +2%/-2% swings
rra2p = rr > 0.02;
rrd2p = -1*(rr < -0.02);
rrs1_2p=[rra2p+rrd2p; NaN];
rrs2_2p=[NaN; rra2p+rrd2p];
rrs12_2p=rrs1_2p.*rrs2_2p; %-1 determines consecutive +5/-5 swings
%cm_2p = (sum(rra2p == 1) + sum(rrd2p == -1))/size(rr,1)*100;
cm_2p=sum(rrs12_2p == -1)/size(rr,1)*100;

[pLF,pHF,LFHFratio,VLF,LF,HF] = HRV.fft_val_fun(RR(~isnan(RR)),Fs);

%- Bin information
edges=[-1 -0.2:0.05:0.2 1];
figure(1)
hist_out=histogram(rr,edges);
hist_bin_values=hist_out.Values;
hist_bin_values=hist_bin_values/sum(hist_bin_values)*100;
%figure(1); 
%histogram(rr);
set(gca, 'YScale', 'log');
grid on
%--
% Separate into 3 zones (+1,0,-1)
% ACCOUNTING for quantization noise
%--
qnt=1*0.5/(Fs*nanmedian(RR));

rrc=rr(~isnan(rr));
rrp1=rrc > qnt;
rrn1=-1*(rrc < -qnt);
rr1_3z=rrp1+rrn1;
total=size(rr1_3z(:,1),1);      % can use length(rr1_3z) as a simpler alternative
% OPTIONAL substitution/smoothing to prevent fragmentation due to quantization noise
for i=2:total-1             
   if ((rr1_3z(i-1,1)==rr1_3z(i+1,1)) & abs(rrc(i,1))<qnt)
        rr1_3z(i,1)=rr1_3z(i-1,1);
    end
end
% of changes introduced by smoothing
count_smoothed_data=numel(find((rrp1+rrn1)~=rr1_3z)); % count_smoothed_data=sum((rrp1+rrn1)~=rr1_3z)
pct_smoothed_data=count_smoothed_data/total*100;
% END - Separate into 3 zones (+1,0,-1)
%--
% COUNTING of inflexion points
count=0;
for i=2:total
    count=count+sum(rr1_3z((i),1)~=rr1_3z((i-1),1));
end
pct_inflex_point=count/length(rr(~isnan(rr)))*100; % count of inflexion points (any transition between 1,0,-1 zones)
%--
%e-figure(2)
%e-subplot(3,1,1)
%e-bar(rr(~isnan(rr)))
%e-subplot(3,1,2)
%e-bar(rr1_3z)
%--
% Counting of consecutive intervals
bp1 = diff([0 rrp1' 0] ==1);
cp1 = find(bp1==-1) - find(bp1==1);
tabulate(cp1) % tabulation of "+1" consecutive intervals and their occurence (count)
bn1 = diff([0 rrn1' 0] ==-1);
cn1 = find(bn1==-1) - find(bn1==1);
tabulate(cn1) % tabulation of "-1" consecutive intervals and their occurence (count)
b0 = diff([0 not(rr1_3z)' 0] ==1);
c0 = find(b0==-1) - find(b0==1);
pct_zero_time=sum(c0)/total*100;
tabulate(c0) % tabulation of "0" consecutive intervals and their occurence (count)
fprintf('SDNN=%.4f, pNN50=%.4f\n',sdnn_RR,pnn50_RR);
fprintf('Val>5%%=%.4f, Val>10%%=%.4f, Val>20%%=%.4f\n',val_out5p,val_out10p,val_out20p);
%fprintf('rr>5%%=%.4f, rr>10%%=%.4f, rr>20%%=%.4f\n',pct_5p,pct_10p,pct_20p);
fprintf('cm>2%%=%.4f, cm>5%%=%.4f, cm>10%%=%.4f, cm>20%%=%.4f\n',cm_2p, cm_5p,cm_10p,cm_20p);
fprintf('Pct Smoothed Date=%.4f, Pct Inflexion Points=%.4f\n',pct_smoothed_data,pct_inflex_point);
fprintf('Histogram Bin Value = %.2f\n',hist_bin_values);
fprintf('# of Valid rr intervals=%d\n',length(rr(~isnan(rr))));
fprintf('LFHFr=%.4f, LF=%.4f, HF=%.4f\n',LF/HF,LF,HF);
figure(2);
scatter([RR;NaN],[NaN;RR])
grid on;
figure(4)
plot([rr;NaN],[NaN;rr]);
grid on;
fprintf('z3e20p=%.4f\n',z3e20p);
