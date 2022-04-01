% Syntax: COMPUTE_Baltic_argo_BEC_colocations
%
% Description
%
% Master script compute colocations of Argo (ref) to satellite grid 
% as computed using function Baltic_argo_BEC_colocations.m
clc; clear

%% Script inputs
iyear = 2013%2011:2013;
imonth = 1:12;
bin_data = 1; % [bin_data = 1] bin REF, or not [bin_data = 0] - go to section 1.4

%% Bin stiring to store colocatios
if bin_data == 1
    bin_str = 'BINNED';
else
    bin_str = 'NOBIN';
end



%iregion = 1;    % [1] Atlantic; [2] Pacific; [3] ALL-Arctic
for idata_type = [4:8]
% Options for idata_type
%              [4] Nominal, 
%              [5] Nodal Sampling, 
%              [6] Global BEC (v001); 
%              [7] CCI+SSS (v01.07); 
%              [8] Model (NEMO) data in Baltic
ibasin = 9;  % Basin number: 9 ([9] Baltic, [7] Arctic)



%% Get SSS-satellite data within a given distance from each Argo float
r = 0; % Radius distance platform (in Km)
lat_factor = 1.25; % Areas at high-lats are about 25% bigger than at equator

r_str = num2str(r);

r2 = r*lat_factor; % multiply lat_factor to get a more accurate sampling radius

depth_ref = 10; % Reference depth (m), by gral. assumtion 10 m


[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

iREGION = [5];

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

% Colocation fn_out for each Regional study in the Arctic North NA and PA
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


path_root = ('/Volumes/Rogue/Data/');

folder_data = [path_root...
    'SSS/' basin_str...
    '/BEC/Validation/indata/argo/Colocations/' bin_str '/monthly/' version_str '/'];

folder_out = [path_root...
    'SSS/' basin_str...
    '/BEC/Validation/indata/argo/Colocations/' bin_str '/cat_files/' version_str '/'];
foldercheck_raf(folder_out);


for yy = 1:length(iyear)
    for mm = 1:length(imonth)
        
        iYEAR = iyear(yy);
        iMONTH = imonth(mm);
        
        %         [vars_out] =...
        %             Baltic_argo_BEC_colocations_v2r0(iYEAR,iMONTH,iregion,idata_type);
        [vars_out] =...
            Baltic_argo_BEC_colocations_v3r0(iYEAR,iMONTH,iregion,idata_type);

    
    end
end


year_start_str = num2str(iyear(1));
year_end_str = num2str(iyear(end));


% Output variables to save
argo_vars_out = {...
    'platform_irange','ID_irange',...
    'PRES_irange','SALT_irange','TEMP_irange',...
    'lon_irange','lat_irange','JULD_irange','time_irange'...
    'TEMP_intp','SALT_intp','PRES_intp'};

sat_vars_out = {...
    'sss_irange',...
    'sss_error_irange',...
    'lon_sss_irange','lat_sss_irange'};

colocated_vars_out = {...
    'sss_irange_mn',...
    'sss_irange_md',...
    'sss_error_irange_mn'};

% variables to concatenate
cat_vars = [argo_vars_out sat_vars_out colocated_vars_out];

save_vars_out = [cat_vars {'version_str'}];

ncount = 0; % set counter to zero
for yy = 1: length(iyear)
    for mm = 1:length(imonth)
        
        ncount = ncount+1;
        
        iYEAR = iyear(yy);
        iMONTH = imonth(mm);
        
        fn_out = [folder_out basin_str '_' region '_argo_R' r_str '_'...
            sprintf('%02.0f',depth_ref) 'm_COLOCATIONS_'...
            year_start_str '_' year_end_str '_'...
            version_str '.mat'];
        
        fn_out_exist = exist(fn_out,'file');
        
        
        if ncount == 1
            
            fn1 = [folder_data basin_str '_' region '_argo_R' r_str '_'...
                sprintf('%02.0f',depth_ref) 'm_COLOCATIONS_'...
                sprintf('%02.0f',iYEAR) sprintf('%02.0f',iMONTH) '_'...
                version_str '.mat'];
            
            if exist(fn1,'file') == 2
                A = load(fn1);
                
            else
                error([fn1 'does not exist'])
            end
            
            for nn = 1:length(cat_vars)
                this_param = cat_vars{nn};
                eval([this_param ' = A.' this_param ';']);
                
            end
            
        elseif ncount > 1
            fn2 = [folder_data basin_str '_' region '_argo_R' r_str '_'...
                sprintf('%02.0f',depth_ref) 'm_COLOCATIONS_'...
                sprintf('%02.0f',iYEAR) sprintf('%02.0f',iMONTH) '_'...
                version_str '.mat'];
            
            A = load(fn_out);
            
            if exist(fn2,'file') == 2
                B = load(fn2);
            else
                error([fn2 'does not exist'])
            end
            
            
            for nn = 1:length(cat_vars)
                this_param = cat_vars{nn};
                eval([this_param ' = cat(2,A.' this_param  ',B.' this_param ');']);
                
            end
        end
        
        save(fn_out,save_vars_out{:});
        
    end
end


end



end % loop idata_type




%% Compute STATS of the colocations
% run script for same period as the colocations
% COMPUTE_STATS_Baltic_argo_BEC_colocations.m
