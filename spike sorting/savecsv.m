function savecsv(~,~)
b=evalin('base','b');

hasspikes=isfield(b,'Unit');

try
b.spiketimes= evalin('base','spiketimes');
end

if hasspikes
    savespikes=questdlg('Save Spiketimes?','Save?');
    savelocation=questdlg('Save .csv to same folder?','Save?');
else
    savespikes=0;
    savelocation=questdlg(sprintf('No unit channel\nSave .csv to same folder?'),'Save?');
end

monkey=strtok(b.filename,'_');
%Exotropes: Patos, Pilchuck 
exa= strcmp(questdlg(['Is ',monkey,' an Exotrope?'],'EXO?','No','Yes','No'),'Yes');

if strcmp(savespikes,'Yes')
    b.spikes=evalin('base','spikes');
    save([b.filepath, b.filename(1:end-4), '-sorted.mat'],'-struct','b')
end

if hasspikes
    %create spike density function
    [sdf, rasters]=makesdf(b,20);
end
%Determine which channel is the right eye:
try
    verg=b.H_Eye2.values-b.H_Eye.values;
catch
    verg=mean(b.H_Eye2.values)-mean(b.H_Eye.values);
end
% if mean(verg) > 0  && exa %Probably Eye2 is RIGHT
%     rep=b.H_Eye.values;%horizontal right eye position
%     lep=b.H_Eye2.values;%horizontal right eye position
%     repV=b.V_Eye.values;%vertical right eye position
%     lepV=b.V_Eye2.values;%vertical right eye position
% else
%     lep=b.H_Eye.values;%horizontal right eye position
%     rep=b.H_Eye2.values;%horizontal right eye position
%     lepV=b.V_Eye.values;%vertical right eye position
%     repV=b.V_Eye2.values;%vertical right eye position
% end  
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
trimlength=min(length(rep),length(lep));

if length(rep) ~= length(lep)
    
    rep=rep(1:trimlength);
    repV=repV(1:trimlength);
    thp=thp(1:trimlength);
    tvp=tvp(1:trimlength);    
    thp2=thp2(1:trimlength);
    tvp2=tvp2(1:trimlength);
    lep=lep(1:trimlength);
    lepV=lepV(1:trimlength);
end

% if length(sdf)> trimlength
%     sdf=sdf(1:trimlength);
% end

if length(tvp)>trimlength || length(thp2)> trimlength
    thp=thp(1:trimlength);
    tvp=tvp(1:trimlength);
    thp2=thp2(1:trimlength);
    tvp2=tvp2(1:trimlength);
end

if hasspikes
    if length(rasters) > trimlength
        rasters=rasters(1:trimlength);
    end
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

if hasspikes
    t=table(rasters,rep,rev,repV,revV,...
        lep,lev,lepV,levV,thp,tvp,thp2,tvp2,...
        'variablenames',{'rasters','rep','rev','repV','revV'...
        'lep','lev','lepV','levV','thp','tvp','thp2','tvp2'});
else
    t=table(rep,rev,repV,revV,...
        lep,lev,lepV,levV,thp,tvp,thp2,tvp2,...
        'variablenames',{'rep','rev','repV','revV'...
        'lep','lev','lepV','levV','thp','tvp','thp2','tvp2'});
end

if exa
    defaultname=[b.filepath, b.filename(1:end-4),'-EXO','.csv'];
else
    defaultname=[b.filepath, b.filename(1:end-4), '.csv'];
end

ww= waitbar(0,'SAVING...');
if strcmp(savelocation,'Yes')
    writetable(t,defaultname)
else
    [filename, filepath]=uiputfile('*.csv','Save Table','~/data');
    display([filepath filename])
    % assignin('base','t',t)
    writetable(t,[filepath filename])
end
close(ww)

function [sdf, rasters]=makesdf(b, stdsize)
rasters=zeros([1,length(b.H_Eye.values)]);
rasters(floor(b.spiketimes/50)+1)=1;
gaus=fspecial('gaussian',[1 stdsize*10],stdsize)*1000;
sdf=conv(rasters,gaus,'same')';
rasters=rasters';
