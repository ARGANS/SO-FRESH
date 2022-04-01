% Syntax: COMPUTE_STATS_Baltic_argo_BEC_colocations
%
% Description
%
% Master script compute stats of dSSS as computed using function
% Baltic_argo_BEC_colocations_stats.m
%
%
% Options for idata_type
%              [4] Nominal,
%              [5] Nodal Sampling,
%              [6] Global BEC (v001);
%              [7] CCI+SSS (v01.07);
%              [8] Model (NEMO) data in Baltic
%
% iregion: There are four Baltic+ study regions (see Baltic+ DUM, p. 28):
% [1] Arkona Basin           [ArB] (55°31'30.22"N,  16°16'15.06"E)
% [2] Bothian Sea            [BOS] (61°55'50.19"N,  19°14'46.24"E)
% [3] Gulf of Finland        [GOF] (59°35'55.07"N,  23°07'27.96"E)
% [4] Northern Baltic Proper [NBP] (57°36'10.05"N,  19°48'25.78"E)
% [5] ALL-Baltic             [ALL]
%
% iships: There are 7-ship route data set (see Baltic+ DUM, p. 15):
% Ship routes are available within the time range shown
% [1] BalticQueen   [2015:2018]
% [2] FinnMaid      [2011:2018] *
% [3] Romantica     [2013:2015] *
% [4] SiljaSernade  [2010, 2014:2018] *
% [5] Transpaper    [2011:2018] *
% [6] Victoria      [2015:2016] *
%           * Ship routes within the time range Baltic+ (2011-2013)
%       
%
%
% =========================================================================
% Author: rcatany
%
% =========================================================================

% set path to data
path_root = ('/Volumes/Rogue/Data/');



%% Script inputs
ibasin = 9;     % Basin number: 9 (Baltic)

[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

iDATA_TYPE = [4:8];
iREGION = [5];
flag_on = 1; % [flag_on= 1] apply flags as (qc = 49 "good value data"); or not [flag_on = 0]
bin_data = 1; % [bin_data = 1] bin REF, or not [bin_data = 0]
REF_str = 'ferrybox'; %input the string of the reference dataset (eg. SDN, ARGO, TARA, etc.)

ROUTES = 5;%[1,2,3,4,5];

r = 0; % (km) radius from ref-2-sat
r_str = num2str(r);

% ship routes
SHIPS = {'BalticQueen','FinnMaid','Romantika',...
    'SiljaSernade','Transpaper','Victoria'};


%% set flag_str to notice filename
if flag_on == 1
    flag_str = 'FLAG#ON';
elseif flag_on == 0
    flag_str = 'FLAG#00';
end

%% set bin_str to notice in filename
if bin_data == 0
    bin_str = 'NOBIN';
elseif bin_data == 1
    bin_str = 'BINNED';
end



for SH = 1:length(ROUTES)
    iships = ROUTES(SH);
    
    if iships == 1 % BalticQueen
        ships_str = SHIPS{1};
        iyear = 2015:2018;
        
    elseif iships == 2 % FinnMaid
        ships_str = SHIPS{2};
        %iyear = 2011:2018;
        iyear = 2011:2013;
        
    elseif iships == 3 % Romantica
        ships_str = SHIPS{3};
        % iyear = 2013:2015;
        iyear = 2013;
        
    elseif iships == 4 % SiljaSernade
        ships_str = SHIPS{4};
        % iyear = 2010:2018;
        iyear = 2011:2013;
        
    elseif iships == 5 % Transpaper
        ships_str = SHIPS{5};
        %iyear = 2011:2018;
        iyear = 2011:2013;
        
    elseif iships == 6 % Victoria
        ships_str = SHIPS{6};
        iyear = 2015:2016;
    end
    
    % get year strings to call filenames
    year1_str = num2str(iyear(1));
    year2_str = num2str(iyear(end));
    


for dd = 1:length(iDATA_TYPE)
    idata_type = iDATA_TYPE(dd);
    
    for rr = 1:length(iREGION)
        iregion = iREGION(rr);
        
        % string of the satellite version
        if idata_type == 1
            version_str = 'plusv3r0';
            
        elseif idata_type == 3
            version_str = 'v2r0';
            
        elseif idata_type == 4
            version_str = 'v1r0_NOM';
            
        elseif idata_type == 5
            version_str = 'v1r0_NS';
            
        elseif idata_type == 6
            version_str = 'v001';
            
        elseif idata_type == 7
            version_str = 'CCI+SSSv1.7';
            
        elseif idata_type == 8
            version_str = 'REANALYSIS_PHY_003_011';
            
        end
        
        
        % Colocatiion fn_out for each Regional study in the Baltic
        if strcmpi(iregion,'ARB') || iregion == 1     % Arkona Basin
            region = 'ARB';
            
        elseif strcmpi(iregion,'BOS')|| iregion == 2  % Bothian Sea
            region = 'BOS';
            
        elseif strcmpi(iregion,'GOF')|| iregion == 3 % Gulf of Finland
            region = 'GOF';
            
        elseif strcmpi(iregion,'NBP')|| iregion == 4 % North Baltic Proper
            region = 'NBP';
            
        elseif strcmpi(iregion,'ALL') || iregion == 5 % all-Baltic region stats
            region = 'ALL';
            
        end
        
        
        %% Statistical computation of dSSS
        
        folder_data =([path_root...
            'SSS/' basin_str '/BEC/Validation/indata/' REF_str '/Colocations/'...
            bin_str '/' flag_str '/cat_files/'...
            version_str '/']);
        
        folder_out = [path_root...
            'SSS/' basin_str...
            '/BEC/Validation/indata/' REF_str '/Colocations/'...
            bin_str '/netCDF/' flag_str '/cat_files/' version_str '/'];
        foldercheck_raf(folder_out);
        
        fn = [basin_str '_' ships_str '_' REF_str...
            '_R' r_str   '_10m_COLOCATIONS_'...
            year1_str '_' year2_str '_' version_str]; 
        
        fn_in = [folder_data fn];
        fn_exist = exist([fn_in '.mat'],'file');
        
        if fn_exist == 2
            
            [vars_out] = Baltic_ferrybox_BEC_colocations_stats(fn_in);
            
            
            %% Convert matlab file with stats to netCDF
            fn = [basin_str '_' ships_str '_ferrybox_R25_10m_COLOCATIONS_'...
                year1_str '_' year2_str '_' version_str '_STATS'];
            
            fn_in = [folder_data fn '.mat'];
            fn_out = ([folder_out fn '.nc']);
            
            %write_netcdf_ferrybox_Colocations_STATS(fn_in,fn_out)
            
        end
        
    end
end

end
