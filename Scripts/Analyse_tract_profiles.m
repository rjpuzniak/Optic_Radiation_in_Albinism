%% Initialization of paths, pipeline flags and metadata
% Subjects
controls = {'CON_1','CON_2','CON_3','CON_4','CON_5','CON_6','CON_7','CON_8'};
albinism = {'ALB_1','ALB_2','ALB_3','ALB_4','ALB_5','ALB_6','ALB_7','ALB_8','ALB_9'};
subjects = {controls, albinism};
sub_group = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1];

root_folder='/home/rjp/1_OVGU/1_Connectivity_in_albinism/';

% Data files structured as Vistasoft and AFQ requires
data_folder = strcat(root_folder,'Data');
sub_dirs = strcat(data_folder,'/sub-',subjects,'/dtiinit');
sub_dirs128 = strcat(data_folder,'/sub-',subjects,'/dtiinit/dti128trilin');

% Folders specific to analysis
analysis_folder=strcat(root_folder,'1_OR_albinism/');
tracks_main_folder=strcat(analysis_folder,'Tracking_results/');
afq_folder=strcat(analysis_folder,'AFQ/');
afq_base_folder=strcat(afq_folder,'AFQ_base/');

% Analysis variants
variants_lgn={'default','thalamic'};
variants_tracking={'Mrtrix', 'Contrack'};

% Flags for types of analyses to be done
overview=1;
clusters=1;
averaging=1;

visualize_mrtrix_tracks=0;
caculate_voxels_overlap=0;
calculate_correlations = 0;
group_plot=0;
render_single_fibers = 0;
render_overlaid_fibers=0;
tract_core_voxels= 0;

calculate_correlations = 0;
group_plot=0;
render_single_fibers = 0;
render_overlaid_fibers=0;
tract_core_voxels= 0;


% Analysis-specific parameters
percentiles = [2.5 97.5];    
properties={'FA','MD','RD','AD'};

for a=1:length(variants_lgn)
    variant_lgn=variants_lgn{a};
    for b=1:length(variants_tracking)
        variant_tracking=variants_tracking{b};
        
        % Define directory with results
        results_folder=strcat(afq_folder,variant_tracking,'_',variant_lgn);
        mkdir(results_folder);
        
        % Load and recompute the AFQ data accordingly to the new params
        load(sprintf('%s/afq_%s_%s.mat',results_folder, variant_tracking, variant_lgn));    
        load(sprintf('%s/patient_data_%s_%s.mat',results_folder, variant_tracking, variant_lgn));
        load(sprintf('%s/control_data_%s_%s.mat',results_folder, variant_tracking, variant_lgn));
        load(sprintf('%s/norms_%s_%s.mat',results_folder, variant_tracking, variant_lgn));
        load(sprintf('%s/abn_%s_%s.mat',results_folder, variant_tracking, variant_lgn));
        load(sprintf('%s/abnTracts_%s_%s.mat',results_folder, variant_tracking, variant_lgn)); 
        
        afq.params.cutoff = percentiles;

        % Recompute various supporting parameters
        [norms, patient_data, control_data, afq] = AFQ_ComputeNorms(afq); % this uses afq.params.cutoff only
         
        %% Custom_afq with averaged profiles
        custom_afq = afq;

        new_fg = {'Averaged_Optic_Radiation_Foveal', 'Averaged_Optic_Radiation_Parafoveal','Averaged_Optic_Radiation_Medperipheral','Averaged_Optic_Radiation_Peripheral'};
        custom_afq.fgnames = horzcat(afq.fgnames, new_fg);

        in_idx = [21:2:28];
        out_idx = [29:32];

        for b = 1:length(in_idx)
            for c =1:length(properties)   
                prop = properties{c};
                custom_afq.patient_data(out_idx(b)).(prop) = nanmean(cat(3, afq.patient_data(in_idx(b)).(prop), afq.patient_data(in_idx(b)+1).(prop)),3);
                custom_afq.control_data(out_idx(b)).(prop) = nanmean(cat(3, afq.control_data(in_idx(b)).(prop), afq.control_data(in_idx(b)+1).(prop)),3);
            end    
        end

        tract_idx = {[21:28];[29:32]};
              
        %% Statistical tests along tract profiles (property-specific)
        for h=1:length(properties)
            property = properties{h};
        
            figures_folder=sprintf('%s/Figures/%s',results_folder,property);
            mkdir(figures_folder);

            % Identify Patients With Abnormal Diffusion Measurements    
            %[abn, abnTracts] = AFQ_ComparePatientsToNorms(patient_data, norms, afq.params.cutoff, property, 'profile'); %switch between profile and mean

            % Run a t-test to compare given property (FA, MD, RD, AD) along each fiber tract for patients vs. controls
            for jj = 1:length(custom_afq.fgnames)
                [stath(jj,:),statp(jj,:),~,Tstats(jj)] = ttest2(custom_afq.patient_data(jj).(property),custom_afq.control_data(jj).(property), 'tail','both', 'alpha', 0.01);
            end
     
            % Crop the results by 10 first and 10 last nodes
            stath_cropped = stath(:,11:90);
            statp_cropped = statp(:,11:90);         
            
            %% Plot overview
            if overview
    
                output_folder=sprintf('%s/1_Overview',figures_folder);
                mkdir(output_folder);

                for i = 1:size(subjects,2) % 1-controls, 2-patients

                    subj_group = subjects{1,i};

                    if i ==1
                        group_name = 'control';
                    elseif i == 2
                        group_name = 'patient';
                    end

                    for j =1:size(tract_idx,1)

                        tract_group = tract_idx{j};

                            for k = 1:length(tract_group)

                                fg_idx = tract_group(k);
                                fg_name = custom_afq.fgnames{fg_idx};

                                figure('Position', get(0, 'Screensize'));
                                hold on;

                                % Extract mean for control group
                                mean_val = nanmean(custom_afq.control_data(fg_idx).(property)(:,:));
                                std_val = nanstd(custom_afq.control_data(fg_idx).(property)(:,:),1);

                                % Plot mean for control
                                plot(mean_val,'-','Color','k','linewidth',4);
                                plot(mean_val + std_val, '--', 'Color', 'k', 'linewidth', 2);
                                plot(mean_val - std_val, '--', 'Color', 'k', 'linewidth', 2);

                                color = jet(length(subj_group));

                                for l=1:length(subj_group)
                                    % Plot tract profiles for group of interest
                                    plot(custom_afq.(strcat(group_name,'_data'))(fg_idx).(property)(l,:), '-', 'Color',color(l,:), 'linewidth',1.5);
                                end   

                                title(sprintf('%s in %s group',fg_name, group_name), 'Interpreter', 'None');

                                if property=='FA'
                                    ylim([0 1]);
                                end

                                ylabel(property);
                                xlabel('% of length');
                                legend(horzcat({'mean in control group', 'mean+std', 'mean-std'}, subjects{i}), 'Interpreter', 'None', 'FontSize', 14);

                                saveas(gcf, sprintf('%s/%s_%s_%s.jpg', output_folder, property, group_name, fg_name));
                                close                                       
                        end
                    end
                end
            end
    
        %% Plot clusters 
        if clusters

            output_folder=sprintf('%s/2_Clusters',figures_folder);
            mkdir(output_folder);

            % Loop through all fiber groups
            % Run a t-test to compare given property (FA, MD, RD, AD) along each fiber tract for patients vs. controls
            for jj = 1:length(custom_afq.fgnames)
                fg_name = custom_afq.fgnames{jj};

                % Plot mean values for both groups
                figure('Position', get(0, 'Screensize'))
                hold on;

                c=lines(2);

                ctr_mean = nanmean(custom_afq.control_data(jj).(property)(:,:));
                ctr_std =   nanstd(custom_afq.control_data(jj).(property)(:,:),1);

                pat_mean = nanmean(custom_afq.patient_data(jj).(property)(:,:),1);
                pat_std =   nanstd(custom_afq.patient_data(jj).(property)(:,:),1);

                plot(ctr_mean,'-','Color',c(1,:),'linewidth',2)
                plot(ctr_mean + ctr_std, '--', 'Color', c(1,:), 'linewidth', 1)
                plot(ctr_mean - ctr_std, '--', 'Color', c(1,:), 'linewidth', 1)

                plot(pat_mean,'-','Color',c(2,:),'linewidth',2)
                plot(pat_mean + pat_std, '--', 'Color', c(2,:)', 'linewidth', 1)
                plot(pat_mean - pat_std, '--', 'Color', c(2,:), 'linewidth', 1)

                stat_data=stath(jj,:);
                try
                    histogram('BinEdges',1:101,'BinCounts',stat_data,'DisplayStyle','bar','EdgeColor','none') 
                end

                if strcmp(property,'FA')
                    ylim([0 1]);
                end

                title(sprintf('Clusters for group-averaged %s',fg_name), 'Interpreter', 'None')
                ylabel(property)
                xlabel('% of length')
                legend(horzcat({'mean in control group', 'mean+std', 'mean-std', 'mean in albinism', 'mean+std', 'mean-std', 'Pval < 0.01'}), 'Interpreter', 'None', 'FontSize', 14)

                fg_name = strrep(fg_name, ' ','_');

                saveas(gcf, sprintf('%s/%s_%s.jpg', output_folder, property, fg_name));
                close                                       

            end
        end
        
        %% Plot averaging process
        if averaging 

           output_folder=sprintf('%s/3_Averaging',figures_folder);
           mkdir(output_folder);

            group = [21, 22, 29; 
                     23, 24, 30;
                     25, 26, 31;
                     27, 28, 32];

            for i =1:size(group,1)
                descp=custom_afq.fgnames{28+i}

                figure('Position', get(0, 'Screensize'))

                for j = 1:size(group,2)

                    subplot(size(group,2), 1, j);
                    ylim([0 1]);

                    ctr_mean = nanmean(custom_afq.control_data(group(i,j)).(property)(:,:));
                    ctr_std =   nanstd(custom_afq.control_data(group(i,j)).(property)(:,:));

                    pat_mean = nanmean(custom_afq.patient_data(group(i,j)).(property)(:,:));
                    pat_std =   nanstd(custom_afq.patient_data(group(i,j)).(property)(:,:));

                    x = 1:100;

                    x1 = [x, fliplr(x)];
                    x2 = [ctr_mean + ctr_std, fliplr(ctr_mean - ctr_std)];
                    x3 = [pat_mean + pat_std, fliplr(pat_mean - pat_std)];

                    hold on
                    ylim([0 1]);

                    plot(ctr_mean, 'ob');
                    plot(pat_mean, 'or');

                    h1 = fill(x1, x2, '-.b');
                    h1.FaceAlpha = 0.25;

                    h2 = fill(x1, x3, '-.r');
                    h2.FaceAlpha = 0.25;

                    % Add histogram
                    try
                        histogram('BinEdges',1:101,'BinCounts',stath(group(i,j),:),'DisplayStyle','bar','EdgeColor','none', 'FaceColor', 'm') 
                    end

                    ylim([0 1])

                    title(sprintf('Group-average of %s',custom_afq.fgnames{group(i,j)}), 'Interpreter', 'None')
                    ylabel(property)
                    xlabel('% of length')

                    if j==1
                        legend(horzcat({'Mean - control', 'Mean - albinism', '+- std - control', '+- std - albinism', 'p < 0.01', 'Pval < 0.01'}), 'Interpreter', 'None', 'FontSize', 14);
                    end
                end

                % savefig
                saveas(gcf, sprintf('%s/%s_%s.jpg', output_folder, property, descp));
                close 
                
            end
        end
        
        %% Render fibers
        
        %% Correlate
        
        
        end
    end
end

%{

    %{
    if render
               
        subjects = {'CON_1', 'CON_2','CON_3','CON_4','CON_5','CON_6','CON_7','CON_8','ALB_1','ALB_2','ALB_3','ALB_4','ALB_5','ALB_6','ALB_7','ALB_8','ALB_9'};
        mkdir(strcat('/home/rjp/1_OVGU/1_Connectivity_in_albinism/7_OR_segmentation/Segmentation_yoshimine/Results_contrack/Figures/4_Renders/',property))      
                
       for jj = 1:length(custom_afq.fgnames)
            [stath(jj,:),statp(jj,:),~,~] = ttest2(custom_afq.patient_data(jj).(property),custom_afq.control_data(jj).(property), 'tail','both', 'alpha', 0.01);
       end 
         
        for i= 1:length(subjects)
                        
            % Load default tract profiles for given subject
            TractProfile = custom_afq.TractProfiles(i,:); 
            
            % Load fg
            fg = {}; %= dtiReadFibers(fullfile(custom_afq.sub_dirs{i},'fibers','MoriGroups_clean_D5_L4.mat'));
            for k=21:length(custom_afq.patient_data)-6
                fg{k}=dtiReadFibers(fullfile(custom_afq.sub_dirs{i},'fibers',strcat(custom_afq.fgnames{k},'.mat')));
            end
            
            % Add the pvalues and T statistics from the group comparison to the tract
            % profile. This same code could be used to add correlation coeficients or
            % other statistics
            for jj = 1:length(patient_data)
                TractProfile(jj) = AFQ_TractProfileSet(TractProfile(jj),'vals','pval',(statp(jj,:)>0.05)*jj); % p-values are equal 0, fov is 21-22, medperi is 34-24 and peri is 25-26                
                %TractProfile(jj) = AFQ_TractProfileSet(TractProfile(jj),'vals','Tstat',Tstats(jj).tstat);
            end
            
            % Render the tract profiles of p-values     
            crange = [0.05 0.95];
            
            % Create new color map          
            newmap = jet(30);              
            newmap(1:20,:)=repmat([1 1 0],20,1);
            newmap(21:23,:)=repmat([0 0 1],3,1);
            newmap(24:26,:) = repmat([1 0 0],3,1);
            newmap(26:30,:) = repmat([0 1 0],5,1);

            % Set the number of fibers to render. More fibers takes longer
            numfibers = 100;
 
            %if isempty(TractProfile(jj).nfibers)
            %    continue
            %end
            
            % Render eccenctricity of the left hemisphere           
            fg_to_render= 'Left_Optic_Radiation_Foveal';
            match = reshape(strcmp({fg.name}, fg_to_render ), size(fg));
            idx_to_render = find(match);         
            [lightH. fiberMesh, h]  = AFQ_RenderFibers(fg{21},'color',[0 0 1],'tractprofile',TractProfile(21),'val','pval','numfibers',numfibers,'cmap',newmap,'crange',crange,'radius',[0.25 2],'jittershading', 0.2, 'alpha', 0.1);
            
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
            
            b0 = readFileNifti(dt.files.b0);
            AFQ_AddImageTo3dPlot(b0,[35,0,0]);
            
            print(gcf, '-dpng',sprintf('%s/%s/Figures/Subjects/%s/Left_eccentricity',analysis_data_folder,Segmentation_type,subjects{i}))
            close()
            
            % Render eccenctricity of the right hemisphere           
            fg_to_render= 'Right_Optic_Radiation_Foveal';
            match = reshape(strcmp({fg.name}, fg_to_render ), size(fg));
            idx_to_render = find(match);         
            [lightH. fiberMesh, h]  = AFQ_RenderFibers(fg(idx_to_render),'color',[0 0 1],'tractprofile',TractProfile(idx_to_render),'val','pval','numfibers',numfibers,'cmap',newmap,'crange',crange,'radius',[0.25 2],'jittershading', 0.2, 'alpha', 0.1);
            
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
            AFQ_RenderFibers(fg(idx_to_render),'color',[0 1 0],'tractprofile',TractProfile(idx_to_render),'val','pval','numfibers',numfibers,'cmap',newmap,'crange',crange,'radius',[0.25 2],'jittershading', 0.1, 'alpha', 0.1, 'newfig', false);
  
            b0 = readFileNifti(dt.files.b0);
            AFQ_AddImageTo3dPlot(b0,[35,0,0]);
            
            print(gcf, '-dpng',sprintf('%s/%s/Figures/Subjects/%s/Right_eccentricity',analysis_data_folder,Segmentation_type,subjects{i}))
            close()
      
        end
    end
   %}    
        
        
end
        
  
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
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


}%