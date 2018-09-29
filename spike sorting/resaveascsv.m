function resaveascsv()

[filenames, filepath]=uigetfile('*.mat','multiselect','on');
if ~iscell(filenames)
    filenames={filenames};
end

h=waitbar(0,'Converting to .csv');

for i = 1:length(filenames)
    waitbar(i/length(filenames),h,['Converting ' filenames{i},'...'])
    b=load([filepath filenames{i}]);
    if ~isfield(b,'spiketimes')
        display(['No Spikes in' filenames{i}])
        continue
    end
    %create spike density function
    [sdf, rasters]=makesdf(b,20);
    rep=b.H_Eye.values;%horizontal right eye position
    lep=b.H_Eye2.values;%horizontal right eye position
    repV=b.V_Eye.values;%vertical right eye position
    lepV=b.V_Eye2.values;%vertical right eye position
    
    thp=b.H_Targ.values; %horizontal target position
    tvp=b.V_Targ.values; %vertical target position
    
    if length(rep)>length(lep)
        rep=rep(1:length(lep));
        repV=repV(1:length(lep));
    elseif length(lep)<length(rep)
        lep=lep(1:length(rep));
        lepV=lepV(1:length(rep));
    end
    if length(sdf)>length(lep)
        sdf=sdf(1:length(lep));
    end
    
    rev=parabolicdiff(smooth(rep,15),5);%horizontal right eye velocity
    revV=parabolicdiff(smooth(repV,15),5);%vertical right eye velocity
    lev=parabolicdiff(smooth(lep,15),5);%horizontal right eye velocity
    levV=parabolicdiff(smooth(lepV,15),5);%vertical right eye velocity
    
    %sdf=[sdf zeros[1 length(rep)-length(sdf)]] %pad sdf
    
    % t=table(sdf,rep,rev,repV,revV,...
    %     lep,lev,lepV,levV,...
    %     'variablenames',{'sdf','rep','rev','repV','revV'...
    %     'lep','lev','lepV','levV'});
    
    % t=table(sdf,rep,rev,repV,revV,...
    %     lep,lev,lepV,levV,thp,tvp,...
    %     'variablenames',{'sdf','rep','rev','repV','revV'...
    %     'lep','lev','lepV','levV','thp','tvp'});
    
    t=table(rasters',rep,rev,repV,revV,...
        lep,lev,lepV,levV,thp,tvp,...
        'variablenames',{'rasters','rep','rev','repV','revV'...
        'lep','lev','lepV','levV','thp','tvp'});
    
    savename=[filenames{i}(1:end-11) '.csv'];
    display([filepath savename])
    writetable(t,[filepath savename])
    
end
close(h)

function [sdf, rasters]=makesdf(b, stdsize)
datalength=length(b.H_Eye.values);
rasters=zeros([1,datalength]);
rasters(floor(b.spiketimes/50))=1;
rasters=rasters(1:datalength);
gaus=fspecial('gaussian',[1 stdsize*10],stdsize)*1000;
sdf=conv(rasters,gaus,'same')';
