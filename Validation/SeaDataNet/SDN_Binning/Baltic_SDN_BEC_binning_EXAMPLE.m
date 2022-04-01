

folder_in = (['/Volumes/Rogue/Data/SSS/Baltic/BEC/Validation/indata/SDN/BINNED/FLAG#49/monthly/']);

load([folder_in 'Baltic_SDN_BINNED_25KM_20110116_20110124.mat'])

% remove zeros as NaN (to improve pcolor)
count_bin(count_bin==0)=NaN;

figure;
fillmap_baltic; hold on
pcolorm(lat_bin,lon_bin,count_bin);
fillmap_baltic; hold on
colorbar

plotm(lat_raw,lon_raw,'ko','MarkerFaceColor','w','markersize',2)