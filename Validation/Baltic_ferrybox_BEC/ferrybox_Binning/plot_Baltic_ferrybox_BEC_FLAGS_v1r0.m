% snippet [ferrybox]: plot SALT and TEMP time series with and without flag

PRODS = {'RAW',['FLAG ~= ' num2str(flag_good)]};


figure; 
subplot(2,1,1); 
h1 = plot(time_number,SALT,'k.'); hold on; 
h2 = plot(time_number(ind),SALT(ind),'ro');

lg = legend([h1,h2],PRODS,...
    'fontsize',18,'location','SouthWest');


title({...
    [REF_str ' flag good data (FLAG ~= ' num2str(flag_good) ')'];...
    [time_str1 '-' time_str2]});


grid on
grid minor
box on
datetick('x', 'mm/dd/yy')
%xticklabels(datestr(time_number,'dd-mm-yyyy'))
xtickangle(45)


subplot(2,1,2); 
h1 = plot(time_number,TEMP,'k.'); 
hold on; 
h2 = plot(time_number(ind),TEMP(ind),'ro');

lg = legend([h1,h2],PRODS,...
    'fontsize',18,'location','SouthWest');


grid on
grid minor
box on
datetick('x', 'mm/dd/yy')
%xticklabels(datestr(time_number,'dd-mm-yyyy'))
xtickangle(45)
