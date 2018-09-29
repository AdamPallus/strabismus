
function savecsvbehavior(b)

% % [filename, filepath]=uigetfile('.mat','Select Data File','C:\Users\setup\Desktop\Nucleus Prepositus Hypoglossi');
% [filename, filepath]=uigetfile('.mat','Select Data File','C:\Users\setup\Desktop\NRTP Vergence\Raw Data');
% % [filename, filepath]=uigetfile('.mat','Select Data File','C:\Users\setup\Desktop\SOA from Mark');
% % [filename, filepath]=uigetfile('.mat','Select Data File','C:\Users\setup\Desktop\INC Recording');
% 
% 
% if filename == 0
%     return
% end
% 
% b=load([filepath filename]);
b.filepath=filepath;
b.filename=filename;


savelocation=questdlg('Save .csv to same folder?','Save?');
exa= strcmp(questdlg('Are these data from an Exotrope?','EXO?','No','Yes','No'),'Yes');


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
thp2=b.H2_Targ.values; %horizontal target position
tvp2=b.V2_Targ.values; %vertical target position

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


if length(tvp)>trimlength
    thp=thp(1:trimlength);
    tvp=tvp(1:trimlength);
    thp2=thp2(1:trimlength);
    tvp2=tvp2(1:trimlength);
end

    
rev=parabolicdiff(smooth(rep,15),5);%horizontal right eye velocity
revV=parabolicdiff(smooth(repV,15),5);%vertical right eye velocity
lev=parabolicdiff(smooth(lep,15),5);%horizontal right eye velocity
levV=parabolicdiff(smooth(lepV,15),5);%vertical right eye velocity


t=table(rep,rev,repV,revV,...
    lep,lev,lepV,levV,thp,tvp,thp2,tvp2,...
    'variablenames',{'rep','rev','repV','revV'...
    'lep','lev','lepV','levV','thp','tvp','thp2','tvp2'});

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


