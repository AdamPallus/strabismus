function [sdf, rasters]=makesdf(b, stdsize)
datalength=length(b.H_Eye.values);
rasters=zeros([1,datalength]);
rasters(floor(b.spiketimes/50))=1;
rasters=rasters(1:datalength);
gaus=fspecial('gaussian',[1 stdsize*10],stdsize)*1000;
sdf=conv(rasters,gaus,'same')';