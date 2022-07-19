% ----------------------------------------------------------------------------
%         Read electrode coordinates
%         Convert Allen CCF indices to Franklin-Paxinos labels
%         Export table with anatomical label of each electrode
%
%   Based on data from Chon et al. Enhanced and unified anatomical labeling 
%   for a common mouse brain atlas (2020).
%
% ----------------------------------------------------------------------------


%% Start with coordinates within the Allen CCF mouse brain atlas
% in the form [AP1, DV1, ML1
%              AP2, DV2, ML2]

output_folder = 'C:\Users\Kacper\Desktop\LDTg\histology\KK013_CCF\processed';
electrode_points_file = 'C:\Users\Kacper\Desktop\LDTg\histology\KK013_CCF\processed\electrode_points.mat';
load(electrode_points_file);
probes = [1 2 3 4]; % probes (==shanks) tracked on earlier steps
          
% directory of reference files
annotation_volume_location = 'C:\Users\Kacper\Desktop\AllenAtlas\annotation_volume_10um_by_index.npy'; % from the allen inst (see readme)
structure_tree_location = 'C:\Users\Kacper\Desktop\AllenAtlas\structure_tree_safe_2017.csv'; % located in github repo
CCF_to_FP_location =  'C:\Users\Kacper\Desktop\AllenAtlas\CCF_to_FP.csv'; % located in github repo
FP_table_location = 'C:\Users\Kacper\Desktop\AllenAtlas\FP_table_Chon_2020.csv'; % located in github repo
chon_images_loc = 'C:\Users\Kacper\Desktop\AllenAtlas\Suppl_File1_Labels'; % from chon et al (supplementary data 4, https://www.nature.com/articles/s41467-019-13057-w)

% generate values for pixel-to-coordinate transformation
bregma = allenCCFbregma(); % estimated bregma position in reference data space
atlas_resolution = 0.010; % pixels to mm

% should the brain image be dark or light
black_brain = true;
brain_points_color = [.5 .5 1];


%% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end
if ~exist('CCFtoFPtable','var') || ~exist('FPtable','var')
    disp('loading CCF-FP lookup tables...')
    CCFtoFPtable = loadCCFtoFP(CCF_to_FP_location);
    FPtable = loadFPtable(FP_table_location);
end

%% plot empty wire frame brain
fwireframe = plotBrainGrid([], [], [], black_brain); hold on; 
fwireframe.InvertHardcopy = 'off';
figure(fwireframe); hold on

%% Loop through all the shanks and save all results as a table
full_table = [];

for shank = probes
    
brain_points = electrodePoints{shank}.brain_coord;
from_tip = electrodePoints{shank}.from_tip;
    
% to read the manually selected points, not electrodes from each shank:
%brain_points = pointList.pointList{shank,1}(:,[3, 2, 1]); 
        
%% initialize array of region annotations
annotation_CCF = cell(size(brain_points,1),3);    
annotation_FP = cell(size(brain_points,1),3); 

%% process data
% loop through every point to get ROI locations and region annotations
for point = 1:size(brain_points,1)

    % find the annotation, name, and acronym of the current point from
    % Allen CCF data
    ann = av(brain_points(point,1),brain_points(point,2),brain_points(point,3));
    name = st.safe_name{ann};
    acr = st.acronym{ann};

    annotation_CCF{point,1} = ann;
    annotation_CCF{point,2} = name;
    annotation_CCF{point,3} = acr;

    % find the annotation, name, and acronym of the current ROI pixel
    % using Chon et al data synthesizing CCF and Franklin-Paxinos
    [ann_FP, name_FP, acr_FP] = CCF_to_FP(brain_points(point,1), brain_points(point,2), brain_points(point,3), ...
                                          CCFtoFPtable, FPtable, chon_images_loc);

    annotation_FP{point,1} = ann_FP;
    annotation_FP{point,2} = name_FP;
    annotation_FP{point,3} = acr_FP;
    
end

% get coordinates relative to bregm\
ap = -(brain_points(:,1)-bregma(1))*atlas_resolution;
dv = (brain_points(:,2)-bregma(2))*atlas_resolution;
ml = (brain_points(:,3)-bregma(3))*atlas_resolution;

% generate table
shank_tab = repelem(shank, length(from_tip))';
data_table = table(shank_tab, from_tip, annotation_CCF(:,2),annotation_CCF(:,3), annotation_FP(:,2),annotation_FP(:,3),...
                        ap,dv,ml, annotation_CCF(:,1), annotation_FP(:,1),...
         'VariableNames', {'Probe', 'Pos_from_tip', 'CCF_name', 'CCF_abbrv', 'FP_name', 'FP_abbrv', 'AP_location', 'DV_location', 'ML_location', 'CCF_index', 'FP_index'});

full_table = [full_table; data_table];

%% plot results
hp = plot3(brain_points(:,1), brain_points(:,3), brain_points(:,2), '.','linewidth',2, 'color',brain_points_color,'markers',10);   

end

%% Display and save the results
disp(full_table);
writetable(full_table, fullfile(output_folder, 'brain_structures.csv'));
disp('Results saved as .csv table');