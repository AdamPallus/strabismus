function resaveascsv2018()
% [filename, filepath]=uigetfile('.mat','Select Data File',...
%     'C:\Users\setup\Desktop\NRTP Vergence\OMN');
[filenames, filepath]=uigetfile('*.mat','Select Data Files',...
    'C:\Users\setup\Desktop\NRTP Vergence',...
    'multiselect','on');

if ~iscell(filenames)
    filenames={filenames};
end

h=waitbar(0,'Converting to .csv');

for i = 1:length(filenames)
    waitbar(i/length(filenames),h,['Converting ' filenames{i},'...'])
    b=load([filepath filenames{i}]);
    
    if ~isfield(b.s,'SpikeTimes')
        display(['No Spikes in' filenames{i}])
        continue %skip to next
    end
    
    monkey=strtok(filenames{i},'_');
    %Exotropes: Patos, Pilchuck
    if strcmp(monkey,'Patos') || strcmp(monkey,'Pilchuck')
        exa=1;
    else
        exa=0;
    end
    
        %create spike density function
    [~, rasters]=makesdf(b,20);
    %Determine which channel is the right eye:
    try
        verg=b.H_Eye2.values-b.H_Eye.values;
    catch
        verg=mean(b.H_Eye2.values)-mean(b.H_Eye.values);
    end
    
    if mean(verg) < 0  && ~exa %Probably Eye2 is RIGHT
        lep=b.H_Eye.values;%horizontal right eye position
        rep=b.H_Eye2.values;%horizontal right eye position
        lepV=b.V_Eye.values;%vertical right eye position
        repV=b.V_Eye2.values;%vertical right eye position
    else
        rep=b.H_Eye.values;%horizontal right eye position
        lep=b.H_Eye2.values;%horizontal right eye position
        repV=b.V_Eye.values;%vertical right eye position
        lepV=b.V_Eye2.values;%vertical right eye position
    end
    
    %apply Hamming Filter (from Cullen)
    hamming=fir1(51,0.25);
    rep=conv(rep,hamming,'same');
    repV=conv(repV,hamming,'same');
    lep=conv(lep,hamming,'same');
    lepV=conv(lepV,hamming,'same');
    
    thp=b.H_Targ.values; %horizontal target position
    tvp=b.V_Targ.values; %vertical target position
    if isfield(b,'H2_Targ')
        thp2=b.H2_Targ.values; %horizontal target position
        tvp2=b.V2_Targ.values; %vertical target position
    else
        thp2=thp;
        tvp2=tvp;
    end
    
    if abs(length(tvp)-length(tvp2))>1000
        thp2=thp;
        tvp2=tvp;
    end
    rev=parabolicdiff(smooth(rep,15),5);%horizontal right eye velocity
    revV=parabolicdiff(smooth(repV,15),5);%vertical right eye velocity
    lev=parabolicdiff(smooth(lep,15),5);%horizontal right eye velocity
    levV=parabolicdiff(smooth(lepV,15),5);%vertical right eye velocity
    
    variablenames={'rasters','rep','rev','repV','revV'...
        'lep','lev','lepV','levV','thp','tvp','thp2','tvp2'};
    for j = 1:length(variablenames)
        xx(j)=eval(['length(', variablenames{j},')']);
    end
    
    trimlength=min(xx);
    if max(xx) > trimlength
        sprintf('Trimming: %d',max(xx)-min(xx))
        rep=rep(1:trimlength);
        repV=repV(1:trimlength);
        thp=thp(1:trimlength);
        tvp=tvp(1:trimlength);
        thp2=thp2(1:trimlength);
        tvp2=tvp2(1:trimlength);
        lep=lep(1:trimlength);
        lepV=lepV(1:trimlength);
        rev=rev(1:trimlength);
        revV=revV(1:trimlength);
        lev=lev(1:trimlength);
        levV=levV(1:trimlength);
        rasters=rasters(1:trimlength);
    end
    t=table(rasters',rep,rev,repV,revV,...
        lep,lev,lepV,levV,thp,tvp,thp2,tvp2,...
        'variablenames',{'rasters','rep','rev','repV','revV'...
        'lep','lev','lepV','levV','thp','tvp','thp2','tvp2'});
    

    if exa
        savename=[filenames{i}(1:end-11), 'EXO','.csv'];
    else
        savename=[filenames{i}(1:end-11) '.csv'];
    end
    display([filepath savename])
    writetable(t,[filepath savename])
    
end
close(h)

function [sdf, rasters]=makesdf(b, stdsize)
datalength=length(b.H_Eye.values);
rasters=zeros([1,datalength]);
rasters(floor(b.s.SpikeTimes))=1;
rasters=rasters(1:datalength);
gaus=fspecial('gaussian',[1 stdsize*10],stdsize)*1000;
sdf=conv(rasters,gaus,'same')';
