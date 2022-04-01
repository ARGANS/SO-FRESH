% Syntax: COMPUTE_STATS_Baltic_argo_BEC_colocations
%
% Description
%
% Master script compute stats of dSSS as computed using function
% Baltic_argo_BEC_colocations_stats.m
% 
% Input
%
% idata_type
%           [4] Nominal, [5] Nodal Sampling, 
%           [6] Global BEC (v001); 
%           [7] CCI+SSS (v01.07); 
%           [8] Model (NEMO) data in Baltic
%
% % iregion: There are four Baltic+ study regions (see Baltic+ DUM, p. 28):
%                       [1] Arkona Basin           [ArB] (55°31'30.22"N,  16°16'15.06"E)
%                       [2] Bothian Sea            [BOS] (61°55'50.19"N,  19°14'46.24"E)
%                       [3] Gulf of Finland        [GOF] (59°35'55.07"N,  23°07'27.96"E)
%                       [4] Northern Baltic Proper [NBP] (57°36'10.05"N,  19°48'25.78"E)
%                       [5] ALL-Baltic             [ALL] 
% 
% * for more informaiton about Baltic regions type "help Baltic_studyregion.m"
%
% =========================================================================

%% Script inputs
iyear = 2013%2011:2013;
iregion = 5;    
ibasin = 9;  % Basin number: 9 ([9] Baltic, [7] Arctic)
bin_data = 1; % [bin_data = 1] bin REF, or not [bin_data = 0] - go to section 1.4

r = 0; % Radius distance platform (in Km)
r_str = num2str(r);


%% Bin stiring to store colocatios
if bin_data == 1
    bin_str = 'BINNED';
else
    bin_str = 'NOBIN';
end

% Define first and last year of timeseries (stings)
year1_str = num2str(iyear(1));
year2_str = num2str(iyear(end));

% Colocatiion fn_out for each Regional study in the Baltic
if strcmpi(iregion,'ARB') || iregion == 1% Arkona Basin
    region = 'ARB';
    
elseif strcmpi(iregion,'BOS') || iregion == 2 % Bothian Sea
    region = 'BOS';
    
elseif strcmpi(iregion,'GOF') || iregion == 3 % Gulf of Finland
    region = 'GOF';
    
elseif strcmpi(iregion,'NBP') || iregion == 4 % Gulf of Finland
    region = 'NBP';
    
elseif strcmpi(iregion,'ALL') || iregion == 5 % all-Baltic region stats
    region = 'ALL';
    
end

[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

iDATA_TYPE = [4,5,7,8];
iREGION = [5];

for dd = 1:length(iDATA_TYPE)
    idata_type = iDATA_TYPE(dd);
    
    for rr = 1:length(iREGION)
        iregion = iREGION(rr);        
        
        path_root = ('/Volumes/Rogue/Data/');
        
        % define version string
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
        
        
        % define region
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
        
        % Define folder input/output
        folder_data =([path_root...
            'SSS/' basin_str '/BEC/Validation/indata/Argo/Colocations/' bin_str '/cat_files/'...
            version_str '/']);
        
        folder_out = [path_root...
            'SSS/' basin_str...
            '/BEC/Validation/indata/argo/Colocations/netCDF/' bin_str '/cat_files/' version_str '/'];
        foldercheck_raf(folder_out);
        
        % Define filename
        fn = [basin_str '_' region '_argo_R' r_str '_10m_COLOCATIONS_' year1_str '_' year2_str '_' version_str];
        
        fn_in = [folder_data fn];
        
        fn_exist = exist([fn_in '.mat'],'file');
        
        if fn_exist == 2
            
            [vars_out] = Baltic_argo_BEC_colocations_stats_v2r0(fn_in);
            
            
            %% Convert matlab file with stats to netCDF
            
            fn = [basin_str '_' region '_argo_R' r_str '_10m_COLOCATIONS_' year1_str '_' year2_str '_' version_str '_STATS'];
            
            fn_in = [folder_data fn '.mat'];
            
            fn_out = ([folder_out fn '.nc']);
            
            %write_netcdf_ARGO_Colocations_STATS(fn_in,fn_out)
            
        end
        
    end
end

