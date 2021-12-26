%Interesting files: A00002;A00015;A00018;A00019;A00023;A00036;
clear

% -- Section to readin 2017 Physionet Data stored locally
[tm,ecg,fs,siginfo]=rdmat('A00036');
input_ecg = ecg;

% -- Artifically add noise (sine for baseline wander; rand white noise
ecg = (0.01*sin(2*pi*60*tm') + 0.25*rand(size(tm'))) + ecg;
input_ecg = ecg;

% -- Section to load OHSU test data
%load('700000.mat');
%lead7 = ECG12Lead(:,7)/1000;
%input_ecg=lead7;

% -- Section to create fft/visual inspection of freq. components
tm=1/fs;
taxis=tm.*[0:1:length(input_ecg)-1]';
figure(1)
subplot(2,1,1)
plot(taxis,input_ecg)
m=length(input_ecg);    % s6.m total number of sample points in time domain
n=1024;
%n=pow2(nextpow2(m));
y=fft(input_ecg,n);
f=(0:n-1)*fs/n;         % frequency vector
power=abs(y).^2/n;      % power spectrum
subplot(2,1,2)
hold on
plot(f(1:floor(n/2)),10*log10(power(1:floor(n/2))))
xlim([0.001 100])
set(gca, 'XScale', 'log')

% -- Filtering of input signal
[b,a] = butter(12,(30/(fs/2)),'low'); % LPF; 12th order, 30Hz as corner-frequency
ecg = filtfilt(b,a,input_ecg);
[b,a] = butter(1,[0.5/(300/2) 30/(300/2)]); %BandpassF; 1st order; 0.4~30Hz pass-band
% b=[1 -1]; a=[1 -0.99];                      % dc notch filter/recursive filter
ecg = filtfilt(b,a,input_ecg);
raw_ecg = ecg;          % raw_ecg is copy of filtered input signal

% -- FFT of filtered signal
y1=fft(raw_ecg,n);
power1=abs(y1).^2/n;
plot(f(1:floor(n/2)),10*log10(power1(1:floor(n/2))))
hold off

% -- Shreyasi Data function
% ecg_noisecancellation function does the following
% 1st: calls ecg_points.m: rloc results are used subsequently (not p,q,s,t)
% ecg_points using Joachin Behar's qrs_detect2.m are used for rloc calc.
% 2nd: spectogram method to identify locations of noise 
% 3rd: identifies RR segments which are noisy and removes them
% 4th: returns clean RR segments without noise
[ ecg,tm1 ] = ecg_noisecancellation( ecg, fs );
% cleaned [ecg] is now used for p,q,r,s,t location calc.
[ P_index, Q_index, R_index, S_index, T_index] = ecg_points( ecg, fs );

figure(2);
subplot(2,1,1)          
plot(taxis,input_ecg)     % plot of ecg prior to S.Datta spectogram cleaning
subplot(2,1,2)
plot(tm1,ecg)           % plot of ecg post S.Datta spectogram cleaning

% Following code is flat (non-function call) copy of above to understand
% how spectogram method works
    [ ploc, Q_index, rloc, S_index, T_index] = ecg_points( raw_ecg, fs );
    [s,f,t] = spectrogram(raw_ecg,fs/10,0.9*fs/10,10000,fs); % t is mid-point of each time segment
    s1 = [];
    med_rr = median(diff(rloc));    % med_rr not used in this fn
    mov_mean_len = 100;             % mov_mean_len * diff(t) = 1second; this is moving average window
    s1=movmean(sum(abs(s(101:end,:))),mov_mean_len);
    thres = 1.25*median(s1);
    g1= find(s1>thres);
    figure(4)
    plot(t,s1);
    line([0 max(t)],[thres thres]);     % yline(thres) is available is 2019a release
    
    p=find(diff(g1)==1);            % 1 x m; JayR: index of segments of consecutive g1>threshold (EXCEPT LAST ENTRY)
    q=[p;p+1];                      % 2 x m
    n1 = g1(q);                     % 2 x m,this plottgives all the pairs of consecutive numbers
    [rows columns] = size(n1);
    if(rows~=2)
        n1 = n1';
    end
    fact = floor(length(raw_ecg)/length(s1)); % JayR: length(raw_ecg) :-)
    n1 = n1 * fact;                 % JayR: align s1 segments to ecg time samples

    n2 = [];                        % Create continuous segment #s for noisy segments
    if(~isempty(n1))
        for i=1:length(n1(1,:))
            n2 = [n2 n1(1,i):1:n1(2,i)];
        end
    end
    n2 = unique(n2);                % noisy time-stamps in ecg sample numbers :-)
    
    noise_rloc = [];
    if(length(rloc)>=2)
        for i=1:length(rloc)        % for every rloc
            %i
            %pause
            if( ~isempty (intersect(n2, rloc(i))))  %common elements in n2 & rloc
                if(i==1)
                    noise_rloc = [noise_rloc; [rloc(1) rloc(2)]]
                elseif i==length(rloc)
                    noise_rloc = [noise_rloc; [rloc(end-1) rloc(end)]]
                else
                    noise_rloc = [noise_rloc; [rloc(i-1) rloc(i)]; [rloc(i) rloc(i+1)]];
                end
            end    
        end
    end
    
    n3=[];
    if (~isempty(noise_rloc))
        length(noise_rloc)
        for i=1:length(noise_rloc(:,1))
            %i
            n3 = [n3 noise_rloc(i,1):1:noise_rloc(i,2)];
            %pause
        end
    end
    
    more_noise_rloc = [];
    if(~isempty(n2))
        n3=setdiff(n2,noise_rloc);  % JayR: extract noisy time-stamps that are not in noise_rloc
        for i=1:length(n3)
           r=find(rloc<n3(i));      % JayR: r is list of r that is less than noisy-time stamp
           if(isempty(r))
               continue;
           end
           start = rloc(r); start=start(end); % JayR: start is last r-location before noise-time stamp
           r= find(rloc>n3(i));     % JayR: r is list of r that is greater than noisy-time stamp
           if(isempty(r))
               continue;
           end
           stop  = rloc(r); stop=stop(1); % JayR: stop of is first location after noisy-time stamp
           if(isempty(more_noise_rloc))
               more_noise_rloc = [start stop]; % JayR: first entry of more_noise_rloc
           elseif(~find(more_noise_rloc==start)) % JayR: subsequent entries of more_noise_loc
                more_noise_rloc = [more_noise_rloc; [start stop]];
           end;
        end
    end
    noise_rloc = [noise_rloc; more_noise_rloc]; 
    noise_rloc1 = unique(noise_rloc);
    % https://www.mathworks.com/help/signal/examples/signal-smoothing.html
    