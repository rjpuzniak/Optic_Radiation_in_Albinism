% Analysis of the AFQ results

% Initialize folder names and other data
Segmentation_type='Segmentation_yoshimine'

subjects = {'CON_1','CON_2','CON_3','CON_4','CON_5','CON_6','CON_7','CON_8','ALB_1','ALB_2','ALB_3','ALB_4','ALB_5','ALB_6','ALB_7','ALB_8','ALB_9'}

calculate_correlations = 0
group_plot=0
render_single_fibers = 0
render_overlaid_fibers=1
tract_core_voxels= 0

analysis_folder=sprintf('/home/rjp/1_OVGU/1_Connectivity_in_albinism/7_OR_segmentation/%s/Results_contrack',Segmentation_type);
mkdir(sprintf('%s/Figures',analysis_folder))

% Parameters
afq.params.cutoff = [2.5 97.5];    
properties={'FA','MD','RD','AD'};

% Load the data
load(sprintf('%s/afq.mat',analysis_folder));    
load(sprintf('%s/patient_data.mat',analysis_folder));
load(sprintf('%s/control_data.mat',analysis_folder));
load(sprintf('%s/norms.mat',analysis_folder));
load(sprintf('%s/abn.mat',analysis_folder));
load(sprintf('%s/abnTracts.mat',analysis_folder)); 
    
% Recompute various supporting parameters
[norms, patient_data, control_data, afq] = AFQ_ComputeNorms(afq); % this uses afq.params.cutoff only


%% Perform statistical testing
for h=1:length(properties)
        
    property = properties{h};
    mkdir(sprintf('%s/Figures/%s/',analysis_folder,property))

    % Identify Patients With Abnormal Diffusion Measurements
    [abn, abnTracts] = AFQ_ComparePatientsToNorms(patient_data, norms, afq.params.cutoff, property, 'profile'); %switch between profile and mean

    % Run a t-test to compare given property (FA, MD, RD, AD) along each fiber tract for patients vs. controls
    for jj = 1:length(afq.fgnames)
        [stath(jj,:),statp(jj,:),~,Tstats(jj)] = ttest2(afq.patient_data(jj).(property),afq.control_data(jj).(property));
    end
    
    
    % Crop the results by 10 first and 10 last nodes
    stath_cropped = stath(:,11:90);
    statp_cropped = statp(:,11:90);
    
    %% Correlate results with other measures
    if calculate_correlations
    
        % Initialize extracted data, 1st column is controls, 2nd is patients 
         outcome = array2table(nan(length(afq.sub_dirs),10));
         outcome.Properties.VariableNames=afq.fgnames(21:30);
         outcome.Properties.RowNames = subjects;

        % For each patient and each track find mean property value in the
        % biggest abnormal cluster
        for a = 21:length(afq.fgnames)

                tmp = diff([0,stath(a,:)]);

                idstart=find(tmp==1);
                idend=find(tmp==-1);

                max_length = max(idend - idstart);

                max_id = find((idend - idstart) == max_length);

                cluster_indices = idstart(max_id):(idend(max_id)-1);

                outcome{1:8,a-20} = mean(afq.control_data(a).(property)(:,[cluster_indices]),2);
                outcome{9:17,a-20} = mean(afq.patient_data(a).(property)(:,[cluster_indices]),2);

        end

        % Read the decussation index data
        OC_idx = load('/home/rjp/1_OVGU/1_Connectivity_in_albinism/Scripts/Data_from_OC_study/CSD_seed_weights.mat');
        OC_idx = OC_idx.ids_table;
        
        % Translate names from BL to our standard
        subjects_alt={'la21','lw37','nb30','ow93','ta14','uf97','xn78','xs62','ce04','fe21','ib57','kw99','nh50','rx88','sj22','tq63','uh47'}; % [controls albinism]
        group =[1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0];
        M = containers.Map(subjects, subjects_alt);        
        outcome.Properties.RowNames = values(M, (outcome.Properties.RowNames).');
        
        % For each tract plot one vs another
        figure('units','normalized','outerposition',[0 0 1 1])
        for b=1:length(outcome.Properties.VariableNames)
            for c=1:length(OC_idx.Properties.VariableNames)
                
                subplot(length(OC_idx.Properties.VariableNames),length(outcome.Properties.VariableNames),b+(c-1)*length(outcome.Properties.VariableNames))
            
                plot_name = strcat(afq.fgnames{20+b},OC_idx.Properties.VariableNames{c});

                I_d = OC_idx{subjects_alt,c};
                I_F_A = outcome{subjects_alt,b};
                gscatter( I_d, I_F_A, group);
                hold on
                
                if(b==1 && c==1)
                    legend('Albisnim', 'Controls')
                else
                    legend('hide')
                end
                
                if(b==1)
                    ylabel(OC_idx.Properties.VariableNames{c})
                end
                
                if(c==1)
                    title(erase(afq.fgnames{20+b},'Optic_Radiation_'), 'Interpreter','None')
                end
                
                % It's linear regression time!
                
                % For controls
                I_d_corr=I_d(logical(group));
                I_d_corr=[ones(length(I_d_corr),1) I_d_corr];
                
                I_F_A_corr=I_F_A(logical(group));
                
                regres = I_d_corr \ I_F_A_corr;
                yCalc2 = I_d_corr * regres;
                
                Rsq = 1 - sum((I_F_A_corr - yCalc2).^2)/sum((I_F_A_corr - mean(I_F_A_corr)).^2);
                
                
                if(Rsq>0.70)              
                    plot(I_d(logical(group)),yCalc2,'b--');
                else
                    plot(I_d(logical(group)),yCalc2,'c--'); 
                end
                clear I_d_corr I_F_A_corr regres yCalc2 Rsq
                
                 % For albinism
                I_d_corr=I_d(~logical(group));
                I_d_corr=[ones(length(I_d_corr),1) I_d_corr];
                
                I_F_A_corr=I_F_A(~logical(group));
                
                regres = I_d_corr \ I_F_A_corr;
                yCalc2 = I_d_corr * regres;
                
                Rsq = 1 - sum((I_F_A_corr - yCalc2).^2)/sum((I_F_A_corr - mean(I_F_A_corr)).^2)
                
                yCalc2 = I_d_corr * regres;
                if(Rsq>0.70)              
                    plot(I_d(~logical(group)),yCalc2,'b--');
                else
                    plot(I_d(~logical(group)),yCalc2,'r--'); 
                end
                clear I_d_corr I_F_A_corr regres yCalc2 Rsq               
            end
        end 
        print(gcf, '-dpng',sprintf('%s/Figures/%s_weights_correlations',analysis_folder,property))
        close()
    end
         
    %% Plot group tract profiles with information about statistically significant segments
    if group_plot
        
        gnames = {'control', 'patient', 'p<0.05'};
                
        data{1} = AFQ_get(afq,'control_data'); %control
        data{2} = AFQ_get(afq,'patient_data'); %patient
        
        tracts = 1:length(data{2});
        c = lines(length(data)); % define colors
        
        % generate random numbers for the figure windowss
        fignums = ceil(rand(1,max(tracts)).*10000);
    
        % Loop over all the tracts
        for jj = 1 : length(tracts)
            % For each tract loop over the number of groups and plot each
            % on the same plot
            for ii = 1 : length(data)
                figure(fignums(jj)); hold on;
                
                % collect the value of interest
                switch(property)
                    case 'FA'
                        vals = data{ii}(tracts(jj)).FA;
                        label = 'Fractional Anisotropy';
                    case 'RD'
                        vals = data{ii}(tracts(jj)).RD;
                        label = 'Radial Diffusivity';
                    case 'AD'
                        vals = data{ii}(tracts(jj)).AD;
                        label = 'Axial Diffusivity';
                    case 'MD'
                        vals = data{ii}(tracts(jj)).MD;
                        label = 'Mead Diffusivity';
                end
                
                % crop the vals
                vals = vals(:,11:90);
                
                % number of subjects with measurements for this tract
                n  = sum(~isnan(vals(:,1)));
                % group mean diffusion profile
                m  = nanmean(vals);
                % standard deviation at each node
                sd = nanstd(vals);
                % standard error of the mean at each node
                se = sd./sqrt(n);
                % plot the mean
                h(ii) = plot(m,'-','Color',c(ii,:),'linewidth',3);
                % plot the confidence interval
                plot(m+se,'--','Color',c(ii,:));
                plot(m-se,'--','Color',c(ii,:));
                % label the axes etc.
                xlabel('Location','fontName','Times','fontSize',12);
                ylabel(label,'fontName','Times','fontSize',12);
                title(afq.fgnames{tracts(jj)},'fontName','Times','fontSize',12, 'Interpreter', 'None');
                set(gca,'fontName','Times','fontSize',12);

            end            
            % Highlight segments where mean values are of statistical significant difference                        
            y = ylim;
            stat_data=stath_cropped(jj,:);
            histogram('BinEdges',1:81,'BinCounts',stat_data,'DisplayStyle','bar','EdgeColor','none') 
            ylim(y)
            
            % add a legend to the plot
            legend(h,gnames);
            
            % save plot
            print(gcf, '-dpng',sprintf('%s/Figures/%s/%s',analysis_folder,property,afq.fgnames{jj}))
            close()                        
        end    
    end

    %% Render fibers and tract cores
    
        if render_single_fibers
        
        subjects = {'CON_1', 'CON_2','CON_3','CON_4','CON_5','CON_6','CON_7','CON_8','ALB_1','ALB_2','ALB_3','ALB_4','ALB_5','ALB_6','ALB_7','ALB_8','ALB_9'};
        analysis_data_folder='/home/rjp/1_OVGU/1_Connectivity_in_albinism/7_OR_segmentation';

        for i= 1:length(subjects)
            
            % Create folder
            mkdir(sprintf('%s/%s/Results_contrack/Figures/Subjects/%s',analysis_data_folder,Segmentation_type,subjects{i}));
            
            % Load tract profiles for given subject
            TractProfile = afq.TractProfiles(i,:);  
            
            % Load fg
            fg = dtiReadFibers(fullfile(afq.sub_dirs{i},'fibers','MoriGroups_clean_D5_L4.mat'));
            for k=21:length(patient_data)
                tmp=dtiReadFibers(fullfile(afq.sub_dirs{i},'fibers',strcat(afq.fgnames{k},'.mat')));
                tmp = rmfield(tmp,'pathwayInfo');
                fg(k) = tmp;
                clear tmp
            end
            
            % Add the pvalues and T statistics from the group comparison to the tract
            % profile. This same code could be used to add correlation coeficients or
            % other statistics
            for jj = 1:length(patient_data)
                TractProfile(jj) = AFQ_TractProfileSet(TractProfile(jj),'vals','pval',statp(jj,:)<0.05); % p-values are equal 0, fov is 21-22, medperi is 34-24 and peri is 25-26                
                %TractProfile(jj) = AFQ_TractProfileSet(TractProfile(jj),'vals','Tstat',Tstats(jj).tstat);
            
            cmap = 'jet';
            crange = [0.05 0.95];
            numfibers = 400;
            
            if isempty(TractProfile(jj).nfibers)
                continue
            end
        
        AFQ_RenderFibers(fg(jj),'color',[.8 .8 1],'tractprofile',TractProfile(jj),...
            'val','pval','numfibers',numfibers,'cmap',cmap,...
            'radius',[0.5 4]);
      
        %b0 = readFileNifti(dt.files.b0);
        b0 = readFileNifti(afq.files.images.path{i});
        AFQ_AddImageTo3dPlot(b0,[35,0,0]);
        
        print(gcf, '-dpng',sprintf('%s/%s/Figures/Subjects/%s/%s',analysis_data_folder,Segmentation_type,subjects{i},afq.fgnames{jj}))
        close()
        
            end
                  
        end
    end
    
    if render_overlaid_fibers
        
        subjects = {'CON_1', 'CON_2','CON_3','CON_4','CON_5','CON_6','CON_7','CON_8','ALB_1','ALB_2','ALB_3','ALB_4','ALB_5','ALB_6','ALB_7','ALB_8','ALB_9'};
        analysis_data_folder='/home/rjp/1_OVGU/1_Connectivity_in_albinism/7_OR_segmentation';

        for i= 1:length(subjects)
            
            % Create folder
            mkdir(sprintf('%s/%s/Results_contrack/Figures/Subjects/%s',analysis_data_folder,Segmentation_type,subjects{i}));
            
            % Load tract profiles for given subject
            TractProfile = afq.TractProfiles(i,:); 
            
            % Load fg
            fg = dtiReadFibers(fullfile(afq.sub_dirs{i},'fibers','MoriGroups_clean_D5_L4.mat'));
            for k=21:length(patient_data)
                tmp=dtiReadFibers(fullfile(afq.sub_dirs{i},'fibers',strcat(afq.fgnames{k},'.mat')));
                tmp = rmfield(tmp,'pathwayInfo');
                fg(k) = tmp;
                clear tmp
            end
            
            % Add the pvalues and T statistics from the group comparison to the tract
            % profile. This same code could be used to add correlation coeficients or
            % other statistics
            for jj = 1:length(patient_data)
                TractProfile(jj) = AFQ_TractProfileSet(TractProfile(jj),'vals','pval',(statp(jj,:)>0.05)*jj); % p-values are equal 0, fov is 21-22, medperi is 34-24 and peri is 25-26                
                %TractProfile(jj) = AFQ_TractProfileSet(TractProfile(jj),'vals','Tstat',Tstats(jj).tstat);
            end
            
            % Render the tract profiles of p-values 
             crange = [0 32];

            % Create new color map          
            newmap = jet(32);              
            newmap(1:20,:)=repmat([255/255 20/255 147/255],20,1);
            newmap(26:27,:)=repmat([0 0 1],2,1);
            newmap(28:29,:) = repmat([1 1 0],2,1);
            newmap(30:31,:) = repmat([1 0 0],2,1);
            newmap(32:end,:) = repmat([0 1 0],1,1);

            % Set the number of fibers to render. More fibers takes longer
            numfibers = 100;
 
            if isempty(TractProfile(jj).nfibers)
                continue
            end
            
            % Render eccenctricity of the left hemisphere           
            fg_to_render= 'Left_Optic_Radiation_Foveal';
            match = reshape(strcmp({fg.name}, fg_to_render ), size(fg));
            idx_to_render = find(match);         
            [lightH. fiberMesh, h]  = AFQ_RenderFibers(fg(idx_to_render),'color',[0 0 1],'tractprofile',TractProfile(idx_to_render),'val','pval','numfibers',numfibers,'cmap',newmap,'crange',crange,'radius',[0.25 2],'jittershading', 0.2, 'alpha', 0.1);
            
            fg_to_render= 'Left_Optic_Radiation_Parafoveal';
            match = reshape(strcmp({fg.name}, fg_to_render ), size(fg));
            idx_to_render = find(match);         
            [lightH. fiberMesh, h]  = AFQ_RenderFibers(fg(idx_to_render),'color',[1 1 0],'tractprofile',TractProfile(idx_to_render),'val','pval','numfibers',numfibers,'cmap',newmap,'crange',crange,'radius',[0.25 2],'jittershading', 0.2, 'alpha', 0.1, 'newfig', false);
            
            fg_to_render= 'Left_Optic_Radiation_Medperipheral';
            match = reshape(strcmp({fg.name}, fg_to_render ), size(fg));
            idx_to_render = find(match);     
            try
                AFQ_RenderFibers(fg(idx_to_render),'color',[1 0 0],'tractprofile',TractProfile(idx_to_render),'val','pval','numfibers',numfibers,'cmap',newmap,'crange',crange,'radius',[0.25 2],'jittershading', 0.1, 'alpha', 0.1, 'newfig', false);
            catch
                disp('Problem!')
            end
            
            fg_to_render= 'Left_Optic_Radiation_Peripheral';
            match = reshape(strcmp({fg.name}, fg_to_render ), size(fg));
            idx_to_render = find(match);     
            try
                AFQ_RenderFibers(fg(idx_to_render),'color',[0 1 0],'tractprofile',TractProfile(idx_to_render),'val','pval','numfibers',numfibers,'cmap',newmap,'crange',crange,'radius',[0.25 2],'jittershading', 0.1, 'alpha', 0.1, 'newfig', false);
            catch
                disp('Problem!')
            end
            
            b0 = readFileNifti(afq.files.images.path{i});
            AFQ_AddImageTo3dPlot(b0,[35,0,0]);
            
            print(gcf, '-dpng',sprintf('%s/%s/Figures/Subjects/%s/Left_eccentricity',analysis_data_folder,Segmentation_type,subjects{i}))
            close()
            
            % Render eccenctricity of the right hemisphere                    
            fg_to_render= 'Right_Optic_Radiation_Foveal';
            match = reshape(strcmp({fg.name}, fg_to_render ), size(fg));
            idx_to_render = find(match);         
            [lightH. fiberMesh, h]  = AFQ_RenderFibers(fg(idx_to_render),'color',[0 0 1],'tractprofile',TractProfile(idx_to_render),'val','pval','numfibers',numfibers,'cmap',newmap,'crange',crange,'radius',[0.25 2],'jittershading', 0.2, 'alpha', 0.1);
            
            fg_to_render= 'Right_Optic_Radiation_Parafoveal';
            match = reshape(strcmp({fg.name}, fg_to_render ), size(fg));
            idx_to_render = find(match);         
            [lightH. fiberMesh, h]  = AFQ_RenderFibers(fg(idx_to_render),'color',[1 1 0],'tractprofile',TractProfile(idx_to_render),'val','pval','numfibers',numfibers,'cmap',newmap,'crange',crange,'radius',[0.25 2],'jittershading', 0.2, 'alpha', 0.1, 'newfig', false);
            
            fg_to_render= 'Right_Optic_Radiation_Medperipheral';
            match = reshape(strcmp({fg.name}, fg_to_render ), size(fg));
            idx_to_render = find(match);     
            try
                AFQ_RenderFibers(fg(idx_to_render),'color',[1 0 0],'tractprofile',TractProfile(idx_to_render),'val','pval','numfibers',numfibers,'cmap',newmap,'crange',crange,'radius',[0.25 2],'jittershading', 0.1, 'alpha', 0.1, 'newfig', false);
            catch
                disp('Problem!')
            end
            
            fg_to_render= 'Right_Optic_Radiation_Peripheral';
            match = reshape(strcmp({fg.name}, fg_to_render ), size(fg));
            idx_to_render = find(match);     
            try
                AFQ_RenderFibers(fg(idx_to_render),'color',[0 1 0],'tractprofile',TractProfile(idx_to_render),'val','pval','numfibers',numfibers,'cmap',newmap,'crange',crange,'radius',[0.25 2],'jittershading', 0.1, 'alpha', 0.1, 'newfig', false);
            catch
                disp('Problem!')
            end
            
            b0 = readFileNifti(afq.files.images.path{i});
            AFQ_AddImageTo3dPlot(b0,[35,0,0]);
            
            print(gcf, '-dpng',sprintf('%s/%s/Figures/Subjects/%s/Right_eccentricity',analysis_data_folder,Segmentation_type,subjects{i}))
            close()
      
        end
    end
    
    if tract_core_voxels
    
        
    
        dtiCreateRoiFromFibers
    end
    
end    
    %% TO DO
    
    % Find vulerable regions
    
    % Plot corelations between VA, misrouting, ID
    
    % All plots only for those 8 or 10 OR tracts
% Add FD, FDC and log(FC) values from Mrtrix
%https://community.mrtrix.org/t/calculating-average-fba-metrics-of-specific-tracts/1805

%https://community.mrtrix.org/t/how-to-restrict-fixelcfestats-to-subset-of-tracts/845/2
