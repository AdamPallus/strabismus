function b=spikesort(~,~)

h=evalin('base','handles');
b=evalin('base','b');
data=evalin('base','data');

threshold=str2double(h.thresh.String);

flip=0; %in case of fibers with positive spikes- use negative threshold
if threshold<0
    flip=1;
end

if flip
%     data={b.Unit.values*-1};
%     assignin('base','data',data)
    data{1}=data{1}*-1;
    threshold=threshold*-1;
end

if abs(threshold) < 1 || abs(threshold) >10
    display('Invalid Threshold. Using default of 3.3');
    threshold=3.3;
end
    
addpath(genpath('C:\Users\setup\Documents\GitHub\UltraMegaSort')) 

Fs=50000; %sampling rate of waveform recorder

extrafilter=1;

if extrafilter
    Wp = [ 800  8000] * 2 / Fs;
    Ws = [ 600 10000] * 2 / Fs;
    [N,Wn] = buttord( Wp, Ws, 3, 20);
    [B,A] = butter(N,Wn);
    for j = 1:length(data)
        data2{j} = filtfilt( B, A, data{j} );
    end
    data=data2;
end

% if flip
%     data{1}=data{1}*-1;
% end

spikes = ss_default_params(Fs,'thresh',threshold);
spikes = ss_detect(data,spikes);
spikes = ss_align(spikes);
spikes = ss_kmeans(spikes);
spikes = ss_energy(spikes);
spikes = ss_aggregate(spikes);

% main tool
xxx=figure;
splitmerge_tool(spikes)
close(xxx)
