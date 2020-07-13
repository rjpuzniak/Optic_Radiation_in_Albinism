%% Run AFQ on groups
%% Initialization of paths, pipeline flags and metadata
root_folder='/home/rjp/1_OVGU/1_Connectivity_in_albinism/';
data_folder = strcat(root_folder,'Data');

analysis_folder=strcat(root_folder,'1_OR_albinism/');
tracks_main_folder=strcat(analysis_folder,'Tracking_results/');

afq_folder=strcat(analysis_folder,'AFQ/');
afq_base_folder=strcat(afq_folder,'AFQ_base/');

% Analysis variants
variants_lgn={'default','thalamic'};
variants_tracking={'Mrtrix', 'Contrack'};

% Flags
afq_creation_required=0;
preparation_required=1;

% Subjects
subjects = {'CON_1', 'CON_2','CON_3','CON_4','CON_5','CON_6','CON_7','CON_8','ALB_1','ALB_2','ALB_3','ALB_4','ALB_5','ALB_6','ALB_7','ALB_8','ALB_9'};
sub_group = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1];

sub_dirs = strcat(data_folder,'/sub-',subjects,'/dtiinit');
sub_dirs128 = strcat(data_folder,'/sub-',subjects,'/dtiinit/dti128trilin');

%% Creation of the AFQ structure with new default fibers
for i=1:length(variants_lgn)
    variant_lgn=variants_lgn{i};
    
    for j=1:length(variants_tracking)
        variant_tracking=variants_tracking{j};
        
        if afq_creation_required
            % Add here a code running default AFQ structure creation    
        else
            load(sprintf('%s/AFQ_default_afq_structure.mat', afq_base_folder));
            patient_data= load(sprintf('%s/AFQ_default_patients_data.mat', afq_base_folder));
            control_data = load(sprintf('%s/AFQ_default_control_data.mat', afq_base_folder));
            norms = load(sprintf('%s/AFQ_default_norms.mat', afq_base_folder));
            abn = load(sprintf('%s/AFQ_default_abn.mat', afq_base_folder));
            abnTracts = load(sprintf('%s/AFQ_default_abnTracts', afq_base_folder));       
        end
        
        if preparation_required
            %% Preparation of LGN ROIs (the V1 are shared)
            if(strcmp(variant_lgn,'default'))
                new_lgn_rois = {'LGN_L','LGN_R'};
                old_lgn_rois = {'_LGN_lh_30','_LGN_rh_30'};
            elseif(strcmp(variant_lgn,'thalamic'))
                new_lgn_rois = {'LGN_L','LGN_R'};
                old_lgn_rois = {'ROIlgn_L','ROIlgn_R'};
            end
            
            for k = 1:length(subjects)
                subj = subjects{k};
                % Converting LGN ROIs .nii.gz to .mat and pasting in the
                % dtiinit128trilin/ROIs folder
                for l = 1:length(new_lgn_rois)
                    
                    if(strcmp(variant_lgn,'default'))
                        tmp = niftiRead(sprintf('%s/ROIs_%s/%s/%s%s.nii.gz',analysis_folder,variant_lgn,subj,subj,old_lgn_rois{l}));  
                    elseif(strcmp(variant_lgn,'thalamic'))
                        tmp_path = dir(fullfile('/home/rjp/1_OVGU/1_Connectivity_in_albinism/1_OR_albinism', strcat('ROIs_',variant_lgn), strcat('sub-',subj),'**/rois',strcat(old_lgn_rois{l},'.nii.gz')));
                        tmp = niftiRead(fullfile(tmp_path.folder, tmp_path.name));         
                    end
                        
                    roisToMake=unique(tmp.data);
                    roisToMake=roisToMake(roisToMake~=0);
                    thisRoi = find(tmp.data==roisToMake);
                    [x1,y1,z1] = ind2sub(size(tmp.data), thisRoi);    
                    roiMask = dtiNewRoi(new_lgn_rois{l}, rand(1, 3));
                    roiMask.coords = mrAnatXformCoords(tmp.qto_xyz, [x1,y1,z1]);
                    clear x1 y1 z1;
                    roiMask = dtiRoiClean(roiMask, [], [0 1 0]);
                    dtiWriteRoi(roiMask,sprintf('%s/ROIs/%s.mat',sub_dirs128{k},new_lgn_rois{l}))
                    clear tmp roiMask  
                end
            end          
                
            %% Preparation of tracks
            if(strcmp(variant_tracking,'Mrtrix'))
                tracks_folder=strcat(tracks_main_folder,variant_tracking,'_',variant_lgn,'_clean');   
                fibers_input_names={'Left_LGN_fov','Right_LGN_fov','Left_LGN_parafov','Right_LGN_parafov','Left_LGN_medperi','Right_LGN_medperi','Left_LGN_peri','Right_LGN_peri'};
            elseif(strcmp(variant_tracking,'Contrack'))
                tracks_folder=strcat(tracks_main_folder,variant_tracking,'_',variant_lgn); 
                fibers_input_names = {'Left_Optic_Radiation_Foveal','Right_Optic_Radiation_Foveal','Left_Optic_Radiation_Parafoveal', 'Right_Optic_Radiation_Parafoveal','Left_Optic_Radiation_Medperipheral','Right_Optic_Radiation_Medperipheral','Left_Optic_Radiation_Peripheral','Right_Optic_Radiation_Peripheral'};               
            end
            fibers_output_names = {'Left_Optic_Radiation_Foveal','Right_Optic_Radiation_Foveal','Left_Optic_Radiation_Parafoveal', 'Right_Optic_Radiation_Parafoveal','Left_Optic_Radiation_Medperipheral','Right_Optic_Radiation_Medperipheral','Left_Optic_Radiation_Peripheral','Right_Optic_Radiation_Peripheral'};               
                        
            for k=1:length(subjects)
                subj = subjects{k};
                for l=1:length(fibers_input_names)
                    
                    if(strcmp(variant_tracking,'Contrack'))
                        copyfile(strcat(tracks_folder,'/',subj,'/',fibers_input_names{l},'.mat'),sprintf('%s/sub-%s/dtiinit/dti128trilin/fibers/%s.mat',data_folder,subj, fibers_output_names{l}))
                    elseif(strcmp(variant_tracking,'Mrtrix'))
                        new = dtiImportFibersMrtrix(sprintf('%s/%s/%s_%s.tck',tracks_folder,subj, subj, fibers_input_names{l}));
                        new.fibers = cellfun(@double,new.fibers,'UniformOutput',false); 
                        new.name = fibers_output_names{l};
                        fgWrite(new, sprintf('%s/sub-%s/dtiinit/dti128trilin/fibers/%s',data_folder,subj, fibers_output_names{l}),'mat');         
                    end
          
                end
            end 
            
        end
        
        %% Add new group to the AFQ structure
        % Set parameters for fibers exclusion
        %afq.params.maxDist =5
        %afq.params.maxLen = 4

        % ROIs and fibers already stored in Mrdif folder due to contrack
        % requirements
        lgn_rois={'LGN_L','LGN_R'};
        v1_roi_names={'V1_FOV_L','V1_FOV_R','V1_PARAFOV_L','V1_PARAFOV_R','V1_MEDPERI_L','V1_MEDPERI_R','V1_PERI_L','V1_PERI_R'};
        fibers_output_names = {'Left_Optic_Radiation_Foveal','Right_Optic_Radiation_Foveal','Left_Optic_Radiation_Parafoveal','Right_Optic_Radiation_Parafoveal','Left_Optic_Radiation_Medperipheral','Right_Optic_Radiation_Medperipheral','Left_Optic_Radiation_Peripheral','Right_Optic_Radiation_Peripheral'};               
            
        for m = 1:length(fibers_output_names)
            afq = AFQ_AddNewFiberGroup(afq, fibers_output_names{m}, lgn_rois{rem(m+1,2)+1}, v1_roi_names{m}, [], [], [], strcat(fibers_output_names{m},'.mat'), true);
        end
        
        %% Save tract cores
        cores_output_folder=strcat(tracks_main_folder,'A_core_',variant_tracking,'_',variant_lgn);
        
        % Load dummy mrtrix fiber
        dummy_mrtrix = dtiImportFibersMrtrix('/home/rjp/1_OVGU/1_Connectivity_in_albinism/1_OR_albinism/Tracking_results/Mrtrix_default/CON_1/CON_1_Left_LGN_fov.tck');

        for r = 1:length(subjects)
            subject = subjects{r};
            mkdir(strcat(cores_output_folder,'/',subject));

            for p = 1:length(fibers_output_names)
                fiber = fibers_output_names{p};                
                or = fgRead(sprintf('%s/fibers/%s_clean_D5_L4.mat',sub_dirs128{r}, fiber));
                or.params = dummy_mrtrix.params;
                dtiExportFibersMrtrix(or, strcat(cores_output_folder,'/',subject,'/',subject,'_',fiber,'.tck'));
            end
        end
      
        %% Recompute norms with 2.5 and 97.5% thresholds
        afq.params.cutoff = [2.5 97.5];
        [norms, patient_data, control_data, afq] = AFQ_ComputeNorms(afq);

        %% Identify Patients With Abnormal Diffusion Measurements
        property = 'FA'; % can be also MD, RD, AD
        afq.params.cutoff = [2.5 97.5];
        [abn, abnTracts] = AFQ_ComparePatientsToNorms(patient_data, norms, afq.params.cutoff, property, 'profile'); %switch between profile and mean

        %% Save the results of the AFQ with custom tracts       
        % Define and create output folder
        output_folder=strcat(afq_folder,variant_tracking,'_',variant_lgn);
        mkdir(output_folder);
        
        save(sprintf('%s/afq_%s_%s.mat',output_folder, variant_tracking, variant_lgn),'afq')
        save(sprintf('%s/patient_data_%s_%s.mat',output_folder, variant_tracking, variant_lgn),'patient_data')
        save(sprintf('%s/control_data_%s_%s.mat',output_folder, variant_tracking, variant_lgn),'control_data')
        save(sprintf('%s/norms_%s_%s.mat',output_folder, variant_tracking, variant_lgn),'norms')
        save(sprintf('%s/abn_%s_%s.mat',output_folder, variant_tracking, variant_lgn),'abn')
        save(sprintf('%s/abnTracts_%s_%s.mat',output_folder, variant_tracking, variant_lgn),'abnTracts')

    end
end      
%%
%{


% Flags
Default_AFQ_Done = 1

%Segmentation_type='Segmentation_half' % Segmentation_khazar Segmentation_yoshimine

Segmentation_types={'Segmentation_yoshimine'}

for i = 1:length(Segmentation_types)
    Segmentation_type=Segmentation_types{i}

if ~Default_AFQ_Done
    [AFQbase AFQdata AFQfunc AFQutil AFQdoc AFQgui] = AFQ_directories;

    % Declare relevant folders
    data_folder='/home/rjp/1_OVGU/1_Connectivity_in_albinism/Data';

    % Declare subjects and groups they belong to
    


    % Declare path to subjects' data and images
    sub_dirs = strcat(data_folder,'/sub-',subjects,'/dtiinit/dti128trilin');
    images = cellfun(@(x)[x '/bin/b0.nii.gz'],sub_dirs,'uni', false);

    % Prepare FG files and declare their paths
    fibers_converted = 1

    if ~fibers_converted
        for i=1:length(subjects)
            subj = subjects{i};
            disp(subj)

            wb = dtiImportFibersMrtrix(sprintf('%s/sub-%s/%s_tractogram_300k.tck',data_folder,subj,subj));
            wb.fibers = cellfun(@double,wb.fibers,'UniformOutput',false);
            fgWrite(wb, sprintf('%s/sub-%s/%s_tractogram_300k',data_folder,subj, subj),'mat'); % changed from 300k to 5k

        end
    end
    fibers = strcat(data_folder,'/sub-',subjects,'/',subjects,'_tractogram_300k.mat');


    % Create AFQ
    afq = AFQ_Create('normalization','ants');

    % Check presence of ANTS
    AFQ_get(afq, 'use ANTS');

    % Set AFQ values
    afq = AFQ_set(afq, 'sub_dirs', sub_dirs);
    afq = AFQ_set(afq, 'sub_group',sub_group);
    afq = AFQ_set(afq, 'images', images);

    % Set AFQ values in the loop through all of the subjects
    for j = 1:length(subjects)

        % Make sure there's ROI folder in dtiinit/dti128trilin so that ANTS can
        % output warped ROIs
        mkdir(sprintf('%s/ROIs',sub_dirs{j}))

        subj = subjects{j}

        afq = AFQ_set(afq, 'wholebrain fg path', 'subnum',j,fibers{j});

        afq = AFQ_set(afq, 'ants warp', sprintf('%s/bin/b0Warp.nii.gz',sub_dirs{j}) ,j);

        afq = AFQ_set(afq, 'ants inverse warp', sprintf('%s/bin/b0InverseWarp.nii.gz',sub_dirs{j}) ,j);  

        afq.overwrite.fibers.wholebrain(j)=0;
        afq.overwrite.fibers.segmented(j)=1;
        afq.overwrite.fibers.clean(j) = 1;
        afq.overwrite.vals(j)=1;
    end

    % Recomputation of ROIs will be obtained via modification of
    % AFQ_SegmentFiberGroups inside AFQ_Run - e.g. [2 0] for 2mm distance and
    % no recomputation

    % Run AFQ
    
    % Alternative run for Khazar's segmentation (AFQ_run limited to alb
    % patients with indices 1,4,6,7,8,9)
    %afq = AFQ_set(afq, 'runsubjects',[1 2 3 4 5 6 7 8 9 12 14 15 16 17])
    
    [afq patient_data control_data norms abn abnTracts] = AFQ_run(sub_dirs, sub_group, afq)
        
else
    

end


%}