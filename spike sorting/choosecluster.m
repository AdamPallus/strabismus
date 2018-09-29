function choosecluster(~,~)
h=evalin('base','handles');
b=evalin('base','b');

if ~isfield(b,'spikes')
    
    try
        spikes=evalin('base','spikes');
    catch
        errordlg('Save Spikes from Sort GUI','Saves Spikes First')
        return
    end
else
    spikes=b.spikes;
end

% clusters=unique(spikes.assigns);
clusters=spikes.labels(spikes.labels(:,2)~=4,1);
clear xx
for i =1:length(clusters)
    xx{i}=num2str(clusters(i));
end

[choice cancel]=listdlg('ListString',xx,...
    'PromptString','Choose a Cluster','selectionmode','single');
if cancel ~=0
    spiketimes=spikes.spiketimes(spikes.assigns==str2double(xx{choice}));
    spiketimes=spiketimes*spikes.params.Fs;
else
    display('Cancelled')
end

doplot= questdlg('Analyze Choice','Plot?');

if strcmp(doplot,'Yes')
    b=evalin('base','b');
    a=figure;
    plot(b.Unit.values)
    hold on
    plot(spiketimes,-0.07,'^r')
    title(['Cluster: ',xx{choice},' #spikes: ',num2str(length(spiketimes))])
    a.Position=[35 672 1852 420];
    waitfor(a)
    keep= questdlg('Keep Choice','Keep?');
else
    keep='Yes';
end


if strcmp(keep,'Yes')
    assignin('base','spiketimes',spiketimes)
    h.savecsv.Enable='on';
end



