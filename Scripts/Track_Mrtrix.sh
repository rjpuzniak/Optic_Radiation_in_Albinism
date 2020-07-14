### Initialization
#participants='CON_1 CON_2 CON_3 CON_4 CON_5 CON_6 CON_7 CON_8 ALB_1 ALB_2 ALB_3 ALB_4 ALB_5 ALB_6 ALB_7 ALB_8 ALB_9'
participants='CON_1'

main_folder=/home/rjp/1_OVGU/1_Connectivity_in_albinism
data_folder=$main_folder/Data
scripts=/home/rjp/1_OVGU/1_Connectivity_in_albinism/1_OR_albinism/Scripts

# Tracking parameters
lmax="6 8 10"
FAs="0.03 0.06"
Curv="30 45 60"

# 'Default' approach
echo 'Default approach'

ROIs_folder=$main_folder/1_OR_albinism/ROIs_default
tracts_folder=$main_folder/1_OR_albinism/Tracking_results/Mrtrix_default
mkdir $tracts_folder

for sub in $participants; do

	echo $sub

	output=$tracts_folder/$sub
	mkdir -p $output

	for i in $lmax; do
		for j in $FAs; do
			for k in $Curv; do

					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Left_LGN_fov_$i\_$j\_$k.tck -seed_image $ROIs_folder/$sub/$sub\_LGN_lh_30.nii.gz -include $ROIs_folder/$sub/$sub\_lh.V1_fov_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;
					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Left_LGN_parafov_$i\_$j\_$k.tck -seed_image $ROIs_folder/$sub/$sub\_LGN_lh_30.nii.gz -include $ROIs_folder/$sub/$sub\_lh.V1_parafov_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;
					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Left_LGN_medperi_$i\_$j\_$k.tck -seed_image $ROIs_folder/$sub/$sub\_LGN_lh_30.nii.gz -include $ROIs_folder/$sub/$sub\_lh.V1_medperi_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;
					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Left_LGN_peri_$i\_$j\_$k.tck -seed_image $ROIs_folder/$sub/$sub\_LGN_lh_30.nii.gz -include $ROIs_folder/$sub/$sub\_lh.V1_peri_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;

					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Right_LGN_fov_$i\_$j\_$k.tck -seed_image $ROIs_folder/$sub/$sub\_LGN_rh_30.nii.gz -include $ROIs_folder/$sub/$sub\_rh.V1_fov_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;
					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Right_LGN_parafov_$i\_$j\_$k.tck -seed_image $ROIs_folder/$sub/$sub\_LGN_rh_30.nii.gz -include $ROIs_folder/$sub/$sub\_rh.V1_parafov_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;
					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Right_LGN_medperi_$i\_$j\_$k.tck -seed_image $ROIs_folder/$sub/$sub\_LGN_rh_30.nii.gz -include $ROIs_folder/$sub/$sub\_rh.V1_medperi_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;
					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Right_LGN_peri_$i\_$j\_$k.tck -seed_image $ROIs_folder/$sub/$sub\_LGN_rh_30.nii.gz -include $ROIs_folder/$sub/$sub\_rh.V1_peri_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;

				done
			done
		done

	cd $output

	tckedit $sub\_Left_LGN_fov_* $sub\_Left_LGN_fov.tck
	tckedit $sub\_Left_LGN_parafov_* $sub\_Left_LGN_parafov.tck
	tckedit $sub\_Left_LGN_medperi_* $sub\_Left_LGN_medperi.tck
	tckedit $sub\_Left_LGN_peri_* $sub\_Left_LGN_peri.tck

	tckedit $sub\_Right_LGN_fov_* $sub\_Right_LGN_fov.tck
	tckedit $sub\_Right_LGN_parafov_* $sub\_Right_LGN_parafov.tck
	tckedit $sub\_Right_LGN_medperi_* $sub\_Right_LGN_medperi.tck
	tckedit $sub\_Right_LGN_peri_* $sub\_Right_LGN_peri.tck

	rm

	rm $sub\_Left_LGN_fov_*
	rm $sub\_Left_LGN_parafov_*
	rm $sub\_Left_LGN_medperi_* 
	rm $sub\_Left_LGN_peri_* 

	rm $sub\_Right_LGN_fov_* 
	rm $sub\_Right_LGN_parafov_* 
	rm $sub\_Right_LGN_medperi_* 
	rm $sub\_Right_LGN_peri_* 


	cd $scripts

done

# 'Thalamic' approach
echo 'Thalamic approach'

ROIs_folder=$main_folder/1_OR_albinism/ROIs_default
tracts_folder=$main_folder/1_OR_albinism/Tracking_results/Mrtrix_thalamic
mkdir $tracts_folder

for sub in $participants; do

	echo $sub

	output=$tracts_folder/$sub
	mkdir -p $output

	L_LGN=$main_folder/1_OR_albinism/ROIs_thalamic/sub-$sub/*/rois/ROIlgn_L.nii.gz
	R_LGN=$main_folder/1_OR_albinism/ROIs_thalamic/sub-$sub/*/rois/ROIlgn_R.nii.gz

	for i in $lmax; do
		for j in $FAs; do
			for k in $Curv; do

					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Left_LGN_fov_$i\_$j\_$k.tck -seed_image $L_LGN -include $ROIs_folder/$sub/$sub\_lh.V1_fov_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;
					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Left_LGN_parafov_$i\_$j\_$k.tck -seed_image $L_LGN -include $ROIs_folder/$sub/$sub\_lh.V1_parafov_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;
					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Left_LGN_medperi_$i\_$j\_$k.tck -seed_image $L_LGN -include $ROIs_folder/$sub/$sub\_lh.V1_medperi_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;
					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Left_LGN_peri_$i\_$j\_$k.tck -seed_image $L_LGN -include $ROIs_folder/$sub/$sub\_lh.V1_peri_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;

					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Right_LGN_fov_$i\_$j\_$k.tck -seed_image $R_LGN -include $ROIs_folder/$sub/$sub\_rh.V1_fov_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;
					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Right_LGN_parafov_$i\_$j\_$k.tck -seed_image $R_LGN -include $ROIs_folder/$sub/$sub\_rh.V1_parafov_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;
					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Right_LGN_medperi_$i\_$j\_$k.tck -seed_image $R_LGN -include $ROIs_folder/$sub/$sub\_rh.V1_medperi_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;
					tckgen $data_folder/sub-$sub/$sub\_FOD_lmax_$i.mif $output/$sub\_Right_LGN_peri_$i\_$j\_$k.tck -seed_image $R_LGN -include $ROIs_folder/$sub/$sub\_rh.V1_peri_vol.nii.gz -act $data_folder/sub-$sub/$sub\_5tt.nii.gz -select 50 -max_attempts_per_seed 10000 -stop;

				done
			done
		done

	cd $output

	tckedit $sub\_Left_LGN_fov_* $sub\_Left_LGN_fov.tck
	tckedit $sub\_Left_LGN_parafov_* $sub\_Left_LGN_parafov.tck
	tckedit $sub\_Left_LGN_medperi_* $sub\_Left_LGN_medperi.tck
	tckedit $sub\_Left_LGN_peri_* $sub\_Left_LGN_peri.tck

	tckedit $sub\_Right_LGN_fov_* $sub\_Right_LGN_fov.tck
	tckedit $sub\_Right_LGN_parafov_* $sub\_Right_LGN_parafov.tck
	tckedit $sub\_Right_LGN_medperi_* $sub\_Right_LGN_medperi.tck
	tckedit $sub\_Right_LGN_peri_* $sub\_Right_LGN_peri.tck

	rm

	rm $sub\_Left_LGN_fov_*
	rm $sub\_Left_LGN_parafov_*
	rm $sub\_Left_LGN_medperi_* 
	rm $sub\_Left_LGN_peri_* 

	rm $sub\_Right_LGN_fov_* 
	rm $sub\_Right_LGN_parafov_* 
	rm $sub\_Right_LGN_medperi_* 
	rm $sub\_Right_LGN_peri_* 


	cd $scripts

done
