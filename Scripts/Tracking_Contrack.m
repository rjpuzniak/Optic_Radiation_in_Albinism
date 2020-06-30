%% Initialization
subjects = {'CON_1','CON_2','CON_3','CON_4','CON_5','CON_6','CON_7','CON_8','ALB_1','ALB_2','ALB_3','ALB_4','ALB_5','ALB_6','ALB_7','ALB_8','ALB_9'};

project_folder='/home/rjp/1_OVGU/1_Connectivity_in_albinism/';
analysis_folder = strcat(project_folder,'1_OR_albinism');
data_folder = strcat(project_folder,'Data');

sub_dirs = strcat(data_folder,'/sub-',subjects,'/dtiinit');
sub_dirs128 = strcat(data_folder,'/sub-',subjects,'/dtiinit/dti128trilin');

%subjects_fibers_folder = strcat(data_folder,'/sub-',subjects,'/dtiinit/dti128trilin/fibers');
%% Import dummy mrtrix tracks to copy header
dummy_mrtrix = dtiImportFibersMrtrix('/home/rjp/1_OVGU/1_Connectivity_in_albinism/7_OR_segmentation/Segmentation_yoshimine/OR_tracts/CON_1/CON_1_Left_LGN_fov.tck');

%% Flags

rois_sets = {'default', 'thalamic'};

%% Contrack - preparing shared ROIs for tractography

% Copy V1 ROIs to dtiinit128trilin/ROIs_folder - shared by both methods
new_cort_rois = {'V1_INF_L','V1_INF_R','V1_SUP_L','V1_SUP_R','V1_FOV_L','V1_FOV_R','V1_PARAFOV_L','V1_PARAFOV_R','V1_MEDPERI_L','V1_MEDPERI_R','V1_PERI_L','V1_PERI_R'};
old_cort_rois = {'lh.V1_inf_vol', 'rh.V1_inf_vol','lh.V1_sup_vol', 'rh.V1_sup_vol','lh.V1_fov_vol', 'rh.V1_fov_vol','lh.V1_parafov_vol', 'rh.V1_parafov_vol','lh.V1_medperi_vol', 'rh.V1_medperi_vol','lh.V1_peri_vol', 'rh.V1_peri_vol'};

for j = 1:length(subjects)
   subj = subjects{j};

   for l = 1:length(old_cort_rois)

        tmp = niftiRead(sprintf('/home/rjp/1_OVGU/1_Connectivity_in_albinism/1_OR_albinism/ROIs_default/%s/%s_%s.nii.gz',subj,subj,old_cort_rois{l})); 

        roisToMake=unique(tmp.data);
        roisToMake=roisToMake(roisToMake~=0);
        thisRoi = find(tmp.data==roisToMake);
        [x1,y1,z1] = ind2sub(size(tmp.data), thisRoi);    
        roiMask = dtiNewRoi(new_cort_rois{l}, rand(1, 3));
        roiMask.coords = mrAnatXformCoords(tmp.qto_xyz, [x1,y1,z1]);
        clear x1 y1 z1;
        roiMask = dtiRoiClean(roiMask, [], [0 1 0]);

        dtiWriteRoi(roiMask,sprintf('%s/ROIs/%s.mat',sub_dirs128{j},new_cort_rois{l}));
        clear tmp roiMask
   end
end

% Convert inclusion and exclusion ROIs - shared by both methods
excl_roi_names={'MNI_template_OR_exclusion','MNI_template_OR_inclusion'};

for i = 1:length(subjects)
    subj = subjects{i};

    for j = 1:length(excl_roi_names)

        tmp = niftiRead(sprintf(sprintf('/home/rjp/1_OVGU/1_Connectivity_in_albinism/1_OR_albinism/ROIs_exclusion/%s/%s_%s.nii.gz',subj, subj, excl_roi_names{j}))); 

        roisToMake=unique(tmp.data);
        roisToMake=roisToMake(roisToMake~=0);
        thisRoi = find(tmp.data==roisToMake);
        [x1,y1,z1] = ind2sub(size(tmp.data), thisRoi);    
        roiMask = dtiNewRoi('EXCLUSION', rand(1, 3));
        roiMask.coords = mrAnatXformCoords(tmp.qto_xyz, [x1,y1,z1]);
        clear x1 y1 z1;
        roiMask = dtiRoiClean(roiMask, [], [0 1 0]);

        dtiWriteRoi(roiMask,sprintf(sprintf('/home/rjp/1_OVGU/1_Connectivity_in_albinism/1_OR_albinism/ROIs_exclusion/%s/%s_%s.mat',subj, subj, excl_roi_names{j})));

    end
end

%% Contrack - preparing ROIs specific to given method
% Seed LGN ROIs and all following steps vary between methods
for i=1:length(rois_sets)
    
    roi_set = rois_sets{i};

    % Preparation of ROIs
    switch roi_set
        
        case 'default'
            
            new_lgn_rois = {'LGN_L','LGN_R'};
            old_lgn_rois = {'_LGN_lh_30','_LGN_rh_30'};
           
            for j = 1:length(subjects)

                subj = subjects{j};

                % Converting LGN ROIs .nii.gz to .mat and pasting in the
                % dtiinit128trilin/ROIs folder
                for k = 1:length(new_lgn_rois)

                    tmp = niftiRead(sprintf('/home/rjp/1_OVGU/1_Connectivity_in_albinism/1_OR_albinism/ROIs_%s/%s/%s%s.nii.gz',roi_set,subj,subj,old_lgn_rois{k}));                    
                    roisToMake=unique(tmp.data);
                    roisToMake=roisToMake(roisToMake~=0);
                    thisRoi = find(tmp.data==roisToMake);
                    [x1,y1,z1] = ind2sub(size(tmp.data), thisRoi);    
                    roiMask = dtiNewRoi(new_lgn_rois{k}, rand(1, 3));
                    roiMask.coords = mrAnatXformCoords(tmp.qto_xyz, [x1,y1,z1]);
                    clear x1 y1 z1;
                    roiMask = dtiRoiClean(roiMask, [], [0 1 0]);

                    dtiWriteRoi(roiMask,sprintf('%s/ROIs/%s.mat',sub_dirs128{i},new_lgn_rois{k}))
                    clear tmp roiMask

                end
            end
                                 
        case 'thalamic'
            
            new_lgn_rois = {'LGN_L','LGN_R'};
            old_lgn_rois = {'ROIlgn_L','ROIlgn_R'};
            
            for j = 1:length(subjects)

                subj = subjects{j};
                
                for k = 1:length(new_lgn_rois)
                   
                    tmp_path = dir(fullfile('/home/rjp/1_OVGU/1_Connectivity_in_albinism/1_OR_albinism', strcat('ROIs_',roi_set), strcat('sub-',subj),'**/rois',strcat(old_lgn_rois{k},'.nii.gz')));
                    tmp = niftiRead(fullfile(tmp_path.folder, tmp_path.name));               
                    roisToMake=unique(tmp.data);
                    roisToMake=roisToMake(roisToMake~=0);
                    thisRoi = find(tmp.data==roisToMake);
                    [x1,y1,z1] = ind2sub(size(tmp.data), thisRoi);    
                    roiMask = dtiNewRoi(new_lgn_rois{k}, rand(1, 3));
                    roiMask.coords = mrAnatXformCoords(tmp.qto_xyz, [x1,y1,z1]);
                    clear x1 y1 z1;
                    roiMask = dtiRoiClean(roiMask, [], [0 1 0]);

                    dtiWriteRoi(roiMask,sprintf('%s/ROIs/%s.mat',sub_dirs128{i},new_lgn_rois{k}))
                    clear tmp roiMask
                    
                end
            end           
    end
    
    
    %% Contrack - tractography

    % Initalize Contrack batch parameters
    ctrParams = ctrInitBatchParams;

    ctrParams.projectName  = strcat('OR_',roi_set)
    ctrParams.baseDir      = '/home/rjp/1_OVGU/1_Connectivity_in_albinism/Data';
    ctrParams.subs         = strcat('sub-',subjects);
    ctrParams.dtDir        = 'dtiinit/dti128trilin';
    ctrParams.roiDir       = 'dtiinit/dti128trilin/ROIs';
    ctrParams.roi1         = {'LGN_L','LGN_R','LGN_L','LGN_R','LGN_L','LGN_R','LGN_L','LGN_R'};
    ctrParams.roi2         = {'V1_FOV_L','V1_FOV_R', 'V1_PARAFOV_L', 'V1_PARAFOV_R', 'V1_MED-PERI_L', 'V1_MED-PERI_R','V1_PERI_L','V1_PERI_R'};
    ctrParams.nSamples     = 12500;
    ctrParams.maxNodes     = 250;
    ctrParams.minNodes     = 10;
    ctrParams.stepSize     = 0.5;
    ctrParams.multithread  = 1;

    % Create command and run tracking for Contrack
    [cmd infoFile] = ctrInitBatchTrack(ctrParams);
    system(cmd)

    % Score Contrack tracks
    batchFileName = ctrInitBatchScore(infoFile, 2500, 1)
    system(batchFileName)
    
    %% Cleaning tractography results
    subjects_contrack_folder=strcat(data_folder,'/sub-',subjects,'/dtiinit/dti128trilin/fibers/conTrack/', strcat('OR_',roi_set));
    v1_rois = {'LGN_L_V1_FOV_L','LGN_R_V1_FOV_R','LGN_L_V1_PARAFOV_L','LGN_R_V1_PARAFOV_R','LGN_L_V1_MED-PERI_L','LGN_R_V1_MED-PERI_R','LGN_L_V1_PERI_L','LGN_R_V1_PERI_R'};
    
    fiber_labels={'Left_Optic_Radiation_Foveal','Right_Optic_Radiation_Foveal','Left_Optic_Radiation_Parafoveal','Right_Optic_Radiation_Parafoveal','Left_Optic_Radiation_Medperipheral','Right_Optic_Radiation_Medperipheral','Left_Optic_Radiation_Peripheral','Right_Optic_Radiation_Peripheral'};
    
    % Loop performing segmentation and cleaning of fibers
    for j = 1:length(subjects)

        subject = subjects{j};

        for k = 1:length(v1_rois)
            
            or = mtrImportFibers(sprintf('%s/scoredFG_%s_%s_top2500.pdb', subjects_contrack_folder{j},strcat('OR_',roi_set), v1_rois{k}))

            % remove fibers crossing to the other hemisphere
            or_cleaned = dtiCleanFibers(or, [0 nan nan]); 

            % Load exclusion ROI
            excl_roi = dtiReadRoi(sprintf(sprintf('/home/rjp/1_OVGU/1_Connectivity_in_albinism/1_OR_albinism/ROIs_exclusion/%s/%s_MNI_template_OR_exclusion.mat',subject, subject)));
            incl_roi = dtiReadRoi(sprintf(sprintf('/home/rjp/1_OVGU/1_Connectivity_in_albinism/1_OR_albinism/ROIs_exclusion/%s/%s_MNI_template_OR_inclusion.mat',subject, subject)));

            [or_cleaned_tmp,~,~] = dtiIntersectFibersWithRoi([],{'not'},0.5, excl_roi, or_cleaned);
            [or_cleaned_selected] = dtiIntersectFibersWithRoi([],{'and'},0.5, incl_roi, or_cleaned_tmp);

            [or_cleaned_selected, cont_fg, keep] = dtiIntersectFibersWithRoi([],{'and'},1, roi ,or_cleaned_selected);
            or_cleaned_selected = AFQ_removeFiberOutliers(or_cleaned_selected, 1, 1, 100);
            or_cleaned_selected.name = fiber_labels{k};
            
            % Create output folder depending on approach
            output_folder = strcat(analysis_folder, '/Results/Contrack_',roi_set,'/',subject)
            mkdir(output_folder)

            fgWrite(or_cleaned_selected, sprintf('%s/%s.mat',output_folder, fiber_labels{k}));

            % saving fibers in .tck format
            %system(strcat('mkdir -p',32,project_folder,'1_OR_albinism/OR_tracts_contrack/',subject));
            or_cleaned_selected.params = dummy_mrtrix.params;

            dtiExportFibersMrtrix(or_cleaned_selected, strcat(output_folder,'1_OR_albinism/OR_tracts_contrack/',subject,'/',subject,'_',fiber_labels{k},'.tck'));
        end
    end    
end