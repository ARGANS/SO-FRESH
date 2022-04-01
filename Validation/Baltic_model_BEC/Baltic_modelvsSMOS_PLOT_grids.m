


figure;

plot(lon(:),ylim,'-','color',[1 1 1].*0.5,'linewidth',2); hold on
plot(xlim,lat(:),'-','color',[1 1 1].*0.5,'linewidth',2); hold on

h1 = plot(lon,lat,'ok','MarkerFaceColor','k','MarkerSize',2); hold on
h2 = plot(lon_irange(:),lat_irange(:),'ro','MarkerSize',5);


xlim([5 xmax]); ylim([ymin 68]);

xlabel('longitude (\circ)','fontsize',24)
ylabel('latitude (\circ)','fontsize',24)

fillmap


lg = legend ([h1(1);h2(1)],...
    {'Model gridpoints: 5220','Sampled grid points: 212'},...
    'fontsize',18,'Location','SouthEast');