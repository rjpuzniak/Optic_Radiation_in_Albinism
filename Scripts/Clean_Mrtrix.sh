# Folders
data_folder=/home/rjp/1_OVGU/1_Connectivity_in_albinism/Data
analysis_folder=/home/rjp/1_OVGU/1_Connectivity_in_albinism/1_OR_albinism
scripts_folder=$analysis_folder/Scripts

participants='CON_1 CON_2 CON_3 CON_4 CON_5 CON_6 CON_7 CON_8 ALB_1 ALB_2 ALB_3 ALB_4 ALB_5 ALB_6 ALB_7 ALB_8 ALB_9'

variants='default thalamic'

for variant in $variants; do
	for subj in $participants; do

		# Create tmp folder
		mkdir -p $analysis_folder/Tracking_results/Tmp_Mrtrix_$variant/$subj

		# Apply exclusion and inclusion ROIs
		cd $analysis_folder/Tracking_results/Mrtrix_$variant/$subj

		for file in *; do
			tckedit $file -exclude $analysis_folder/ROIs_exclusion/$subj/$subj\_MNI_template_OR_exclusion.nii.gz -include $analysis_folder/ROIs_exclusion/$subj/$subj\_MNI_template_OR_inclusion.nii.gz $analysis_folder/Tracking_results/Tmp_Mrtrix_$variant/$subj/${file::-4}.tck -force
		done

		# Calculate mean and std of tracks length
		mkdir -p $analysis_folder/Tracking_results/Mrtrix_$variant\_clean/$subj
		cd $analysis_folder/Tracking_results/Tmp_Mrtrix_$variant\_clean/$subj 
		for file in *; do 	
			mean=$(tckstats $file -output mean -quiet)
			std=$(tckstats $file -output std -quiet)

			meanmin=$(echo "($mean-$std)"| bc -l)
			meanplus=$(echo "($mean+$std)"| bc -l)

			tckedit $file -minlength $meanmin -maxlength $meanplus $analysis_folder/Tracking_results/Mrtrix_$variant\_clean/$subj/${file::-4}.tck -force

		done
	done

	# Remove tmp folder
	rm -rf $analysis_folder/Tracking_results/Tmp_Mrtrix_$variant
done

cd $scripts_folder

: <<'END'
END
