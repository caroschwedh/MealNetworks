********************************************************************************
* Project: 		Meal patterns and BMI/weight - PEAS.
* Program: 		meal_PEAS_preg
* Author: 		Carolina Schwedhelm
* Created: 		09/18/2019
* Last Changes: 10/11/2019, by Carolina Schwedhelm
* Origin: 		self made
* Category: 	creation of dataset
* Program before: 	Autoexec
* Short descr.: Creating datasets for main analysis (PREGNANCY):
				PREG (24 datasets) --> 	normal BMI vs. high BMI
		 								adequate GWG vs excessive GWG
										low HEI vs. high HEI
*******************************************************************************;

*infitems_meals_wide --> meal data

*********************************************************************************
********************************************************************************;


*********************************************************************************
* ADD VARIABLES (BMI, GWG, HEI, EPW)
********************************************************************************;

data base; 																					* 458 women;
	set 'N:\DIPHRHBB\PEAS\Data - Updated\Analytic\Main\PEAS-pregnancy\base.sas7bdat';		* 0 BMI missing;
	keep PATID BMIgroup gwg_iom EPW_ever HEI2015_TOT_PREG HEI2015_TOT_POST;					* gwg 91 missing;
run;																						* EPW 184 missing;
*check missings;																			* HEI_preg 93 missing;
proc means data=base n nmiss;																* HEI_post 191 missing;
	var BMIgroup gwg_iom EPW_ever HEI2015_TOT_PREG HEI2015_TOT_POST;
run;
*make HEI tertiles;
proc rank data= base groups=3 out=base;
	var HEI2015_TOT_PREG HEI2015_TOT_POST;
	ranks HEI_PREG_TERT HEI_POST_TERT;
run;
*rename tertiles;
data base; set base;
	if HEI_PREG_TERT = 2 then HEI_PREG_TERT = 3;
	if HEI_PREG_TERT = 1 then HEI_PREG_TERT = 2;
	if HEI_PREG_TERT = 0 then HEI_PREG_TERT = 1;

	if HEI_POST_TERT = 2 then HEI_POST_TERT = 3;
	if HEI_POST_TERT = 1 then HEI_POST_TERT = 2;
	if HEI_POST_TERT = 0 then HEI_POST_TERT = 1;
run;	
/*
*check tertiles;
proc sort data=base; by HEI_PREG_TERT;
proc means data=base;
var HEI2015_TOT_PREG;
by HEI_PREG_TERT;
run;
proc sort data=base; by HEI_POST_TERT;
proc means data=base;
var HEI2015_TOT_POST;
by HEI_POST_TERT;
run;
*/

proc sort data=base; by patid;
proc sort data=infitems_meals_wide; by patid;
data PEAS_meal_networks;
	merge infitems_meals_wide (in=a) base;
	by patid;
	if a;
run;

*export for co-authors;
libname outloc 'V:\Staff Subdirectories\Carolina Schwedhelm\Datasets for Manuscript _1_PEAS networks/';
proc datasets library=work;
	copy in=work out=outloc;
	select PEAS_meal_networks;
run;
quit;
proc export data=PEAS_meal_networks dbms=csv
	outfile="V:\Staff Subdirectories\Carolina Schwedhelm\Datasets for Manuscript _1_PEAS networks/PEAS_meal_networks.csv"
	replace;
run;


*********************************************************************************
* create PREG + BMI datasets;
********************************************************************************;

*breakfast;
data PREG_BMI1_occ1; set PEAS_meal_networks;
	where pregpost = 1 and occ_name = 1 and bmigroup = 1;
run;
data PREG_BMI2_3_occ1; set PEAS_meal_networks;
	where pregpost = 1 and occ_name = 1 and bmigroup in (2,3);
run;


*lunch;
data PREG_BMI1_occ3; set PEAS_meal_networks;
	where pregpost = 1 and occ_name = 3 and bmigroup = 1;
run;
data PREG_BMI2_3_occ3; set PEAS_meal_networks;
	where pregpost = 1 and occ_name = 3 and bmigroup in (2,3);
run;


*dinner;
data PREG_BMI1_occ4; set PEAS_meal_networks;
	where pregpost = 1 and occ_name in (4,5) and bmigroup = 1;
run;
data PREG_BMI2_3_occ4; set PEAS_meal_networks;
	where pregpost = 1 and occ_name in (4,5) and bmigroup in (2,3);
run;

/*
*Nr participants;
proc means data=PREG_BMI2_3_occ4 noprint;
	class patid;
	ways 1;
	var tot_foodamt;
	output out=check mean=;
run;
*/

*********************************************************************************
* create PREG + GWG datasets;
********************************************************************************;

*breakfast;
data PREG_GWG2_occ1; set PEAS_meal_networks;
	where pregpost = 1 and occ_name = 1 and gwg_iom = 2;
run;
data PREG_GWG3_occ1; set PEAS_meal_networks;
	where pregpost = 1 and occ_name = 1 and gwg_iom = 3;
run;

*lunch;
data PREG_GWG2_occ3; set PEAS_meal_networks;
	where pregpost = 1 and occ_name = 3 and gwg_iom = 2;
run;
data PREG_GWG3_occ3; set PEAS_meal_networks;
	where pregpost = 1 and occ_name = 3 and gwg_iom = 3;
run;

*dinner;
data PREG_GWG2_occ4; set PEAS_meal_networks;
	where pregpost = 1 and occ_name in (4,5) and gwg_iom = 2;
run;
data PREG_GWG3_occ4; set PEAS_meal_networks;
	where pregpost = 1 and occ_name in (4,5) and gwg_iom = 3;
run;
/*
*Nr participants;
proc means data=PREG_GWG3_occ4 noprint;
	class patid;
	ways 1;
	var tot_foodamt;
	output out=check mean=;
run;
*/

*********************************************************************************
* create PREG + HEI datasets;
********************************************************************************;
/*
* check distribution of HEI;
proc univariate data=base;
	var HEI2015_TOT_PREG hei2015_tot_post;
	histogram HEI2015_TOT_PREG hei2015_tot_post;			* 57.7 = median PREG;
run;														* 58.2 = median POST;
*/
*breakfast;
data PREG_HEI_T1_occ1; set PEAS_meal_networks;
	where pregpost = 1 and occ_name = 1 and HEI_PREG_TERT = 1;
run;
data PREG_HEI_T3_occ1; set PEAS_meal_networks;
	where pregpost = 1 and occ_name = 1 and HEI_PREG_TERT  = 3;
run;

*lunch;
data PREG_HEI_T1_occ3; set PEAS_meal_networks;
	where pregpost = 1 and occ_name = 3 and HEI_PREG_TERT  = 1;
run;
data PREG_HEI_T3_occ3; set PEAS_meal_networks;
	where pregpost = 1 and occ_name = 3 and HEI_PREG_TERT = 3;
run;

*dinner;
data PREG_HEI_T1_occ4; set PEAS_meal_networks;
	where pregpost = 1 and occ_name in (4,5) and HEI_PREG_TERT  = 1;
run;
data PREG_HEI_T3_occ4; set PEAS_meal_networks;
	where pregpost = 1 and occ_name in (4,5) and HEI_PREG_TERT = 3;
run;

/*
*Nr participants;
proc means data=PREG_HEI_T3_occ1 noprint;
	class patid;
	ways 1;
	var tot_foodamt;
	output out=check mean=;
run;
*/


********************************************************************************;
*save data set file to a local folder on server for importing to sas;
libname out 'C:\Users\schwedhelmramc2\Documents\PEAS data_M1/PREG/';
********************************************************************************;

*preg BMI - BREAKFAST;
data out.PREG_BMI1_occ1;
	set PREG_BMI1_occ1;
	drop 	patid pregpost visit recallno occ_no occ_name _FREQ_ tot_foodamt 
			BMIgroup gwg_iom EPW_ever HEI2015_TOT_PREG HEI2015_tot_post HEI_PREG_TERT HEI_POST_TERT
			food_17 food_40;
	rename 	food_1 = Milk_drinks					food_2 = Milk_desserts					food_3 = Cheese
			food_4 = Poultry						food_5 = Fish_shellfish					food_6 = Meat
			food_7 = Cured_meat						food_8 = Eggs							food_9 = Legumes_nuts_seeds
			food_10 = Breads_refined				food_11 = Breads_whole					food_12 = Cakes_cookies
			food_13 = Pancakes_other				food_14 = Savory_pies					food_15 = Sandwiches
			food_16 = Patties_loaves														food_18 = Tortilla_based			 	
			food_19 = Mayonnaise_salads				food_20 = Pastas_rice					food_21 = Crackers_salty_snacks
			food_22 = Breakfast_cereals_high_sugar	food_23 = Breakfast_cereals_low_sugar	food_24 = Fruit_juice
			food_25 = Whole_fruits					food_26 = Potatoes						food_27 = Potatoes_fried
			food_28 = DG_vegetables					food_29 = RO_vegetables					food_30 = Other_vegetables
			food_31 = Solid_fats					food_32 = Oils							food_33 = Sauces
			food_34 = Soups							food_35 = Sugars_sweets					food_36 = Coffee_tea
			food_37 = Sugar_sweetened_drinks		food_38 = Alcoholic_drinks				food_39 = Water;
run;

data PREG_BMI2_3_occ1;
	set PREG_BMI2_3_occ1;
	rename 	food_1 = Milk_drinks					food_2 = Milk_desserts					food_3 = Cheese
			food_4 = Poultry						food_5 = Fish_shellfish					food_6 = Meat
			food_7 = Cured_meat						food_8 = Eggs							food_9 = Legumes_nuts_seeds
			food_10 = Breads_refined				food_11 = Breads_whole					food_12 = Cakes_cookies
			food_13 = Pancakes_other				food_14 = Savory_pies					food_15 = Sandwiches
			food_16 = Patties_loaves				food_17 = Pasta_based					food_18 = Tortilla_based			 	
			food_19 = Mayonnaise_salads				food_20 = Pastas_rice					food_21 = Crackers_salty_snacks
			food_22 = Breakfast_cereals_high_sugar	food_23 = Breakfast_cereals_low_sugar	food_24 = Fruit_juice
			food_25 = Whole_fruits					food_26 = Potatoes						food_27 = Potatoes_fried
			food_28 = DG_vegetables					food_29 = RO_vegetables					food_30 = Other_vegetables
			food_31 = Solid_fats					food_32 = Oils							food_33 = Sauces
			food_34 = Soups							food_35 = Sugars_sweets					food_36 = Coffee_tea
			food_37 = Sugar_sweetened_drinks												food_39 = Water;
run;
data out.PREG_BMI2_3_occ1;
	set PREG_BMI2_3_occ1;
	drop 	patid pregpost visit recallno occ_no occ_name _FREQ_ tot_foodamt 
			BMIgroup gwg_iom EPW_ever HEI2015_TOT_PREG HEI2015_tot_post HEI_PREG_TERT HEI_POST_TERT
			food_38 food_40;
run;
* dataset with frequency eaten;
proc sort data=PREG_BMI2_3_occ1; by patid visit recallno occ_no;
proc transpose data= PREG_BMI2_3_occ1 out = long1;
	by patid visit recallno occ_no;
	var 	Milk_drinks						Milk_desserts					Cheese
			Poultry							Fish_shellfish					Meat
			Cured_meat						Eggs							Legumes_nuts_seeds
			Breads_refined					Breads_whole					Cakes_cookies
			Pancakes_other					Savory_pies						Sandwiches
			Patties_loaves					Pasta_based						Tortilla_based			 	
			Mayonnaise_salads				Pastas_rice						Crackers_salty_snacks
			Breakfast_cereals_high_sugar	Breakfast_cereals_low_sugar		Fruit_juice
			Whole_fruits					Potatoes						Potatoes_fried
			DG_vegetables					RO_vegetables					Other_vegetables
			Solid_fats						Oils							Sauces
			Soups							Sugars_sweets					Coffee_tea
			Sugar_sweetened_drinks											Water;
run;
data long1; set long1;
	if COL1 = 0 then freq = 0;
	else if COL1 ne 0 then freq = 1;
run;
proc means data=long1 noprint;
	class _NAME_ _LABEL_;
	ways 2;
	var freq;
	output out=PREG_BMI2_3_occ1_FREQ sum=;
run;
data PREG_BMI2_3_occ1_FREQ;
	set PREG_BMI2_3_occ1_FREQ;
	drop _TYPE_ _FREQ_;
run;
%ds2csv (data=PREG_BMI2_3_occ1_FREQ, runmode=b, csvfile=C:\Users\schwedhelmramc2\Documents\PEAS data_M1/PREG/BMI2_3_FREQ.csv);



*preg BMI - LUNCH;
libname out 'C:\Users\schwedhelmramc2\Documents\PEAS data_M1/PREG/';
data out.PREG_BMI1_occ3;
	set PREG_BMI1_occ3;
	drop 	patid pregpost visit recallno occ_no occ_name _FREQ_ tot_foodamt 
			BMIgroup gwg_iom EPW_ever HEI2015_TOT_PREG HEI2015_TOT_POST
			food_40;
	rename 	food_1 = Milk_drinks					food_2 = Milk_desserts					food_3 = Cheese
			food_4 = Poultry						food_5 = Fish_shellfish					food_6 = Meat
			food_7 = Cured_meat						food_8 = Eggs							food_9 = Legumes_nuts_seeds
			food_10 = Breads_refined				food_11 = Breads_whole					food_12 = Cakes_cookies
			food_13 = Pancakes_other				food_14 = Savory_pies					food_15 = Sandwiches
			food_16 = Patties_loaves				food_17 = Pasta_based					food_18 = Tortilla_based			 	
			food_19 = Mayonnaise_salads				food_20 = Pastas_rice					food_21 = Crackers_salty_snacks
			food_22 = Breakfast_cereals_high_sugar	food_23 = Breakfast_cereals_low_sugar	food_24 = Fruit_juice
			food_25 = Whole_fruits					food_26 = Potatoes						food_27 = Potatoes_fried
			food_28 = DG_vegetables					food_29 = RO_vegetables					food_30 = Other_vegetables
			food_31 = Solid_fats					food_32 = Oils							food_33 = Sauces
			food_34 = Soups							food_35 = Sugars_sweets					food_36 = Coffee_tea
			food_37 = Sugar_sweetened_drinks		food_38 = Alcoholic_drinks				food_39 = Water;
run;

data out.PREG_BMI2_3_occ3;
	set PREG_BMI2_3_occ3;
	drop 	patid pregpost visit recallno occ_no occ_name _FREQ_ tot_foodamt 
			BMIgroup gwg_iom EPW_ever HEI2015_TOT_PREG HEI2015_TOT_POST
			food_38 food_40;
	rename 	food_1 = Milk_drinks					food_2 = Milk_desserts					food_3 = Cheese
			food_4 = Poultry						food_5 = Fish_shellfish					food_6 = Meat
			food_7 = Cured_meat						food_8 = Eggs							food_9 = Legumes_nuts_seeds
			food_10 = Breads_refined				food_11 = Breads_whole					food_12 = Cakes_cookies
			food_13 = Pancakes_other				food_14 = Savory_pies					food_15 = Sandwiches
			food_16 = Patties_loaves				food_17 = Pasta_based					food_18 = Tortilla_based			 	
			food_19 = Mayonnaise_salads				food_20 = Pastas_rice					food_21 = Crackers_salty_snacks
			food_22 = Breakfast_cereals_high_sugar	food_23 = Breakfast_cereals_low_sugar	food_24 = Fruit_juice
			food_25 = Whole_fruits					food_26 = Potatoes						food_27 = Potatoes_fried
			food_28 = DG_vegetables					food_29 = RO_vegetables					food_30 = Other_vegetables
			food_31 = Solid_fats					food_32 = Oils							food_33 = Sauces
			food_34 = Soups							food_35 = Sugars_sweets					food_36 = Coffee_tea
			food_37 = Sugar_sweetened_drinks												food_39 = Water;
run;


*preg BMI - DINNER;
libname out 'C:\Users\schwedhelmramc2\Documents\PEAS data_M1/PREG/';
data out.PREG_BMI1_occ4;
	set PREG_BMI1_occ4;
	drop 	patid pregpost visit recallno occ_no occ_name _FREQ_ tot_foodamt 
			BMIgroup gwg_iom EPW_ever HEI2015_TOT_PREG HEI2015_TOT_POST
			food_40;
	rename 	food_1 = Milk_drinks					food_2 = Milk_desserts					food_3 = Cheese
			food_4 = Poultry						food_5 = Fish_shellfish					food_6 = Meat
			food_7 = Cured_meat						food_8 = Eggs							food_9 = Legumes_nuts_seeds
			food_10 = Breads_refined				food_11 = Breads_whole					food_12 = Cakes_cookies
			food_13 = Pancakes_other				food_14 = Savory_pies					food_15 = Sandwiches
			food_16 = Patties_loaves				food_17 = Pasta_based					food_18 = Tortilla_based			 	
			food_19 = Mayonnaise_salads				food_20 = Pastas_rice					food_21 = Crackers_salty_snacks
			food_22 = Breakfast_cereals_high_sugar	food_23 = Breakfast_cereals_low_sugar	food_24 = Fruit_juice
			food_25 = Whole_fruits					food_26 = Potatoes						food_27 = Potatoes_fried
			food_28 = DG_vegetables					food_29 = RO_vegetables					food_30 = Other_vegetables
			food_31 = Solid_fats					food_32 = Oils							food_33 = Sauces
			food_34 = Soups							food_35 = Sugars_sweets					food_36 = Coffee_tea
			food_37 = Sugar_sweetened_drinks		food_38 = Alcoholic_drinks				food_39 = Water;
run;
data out.PREG_BMI2_3_occ4;
	set PREG_BMI2_3_occ4;
	drop 	patid pregpost visit recallno occ_no occ_name _FREQ_ tot_foodamt 
			BMIgroup gwg_iom EPW_ever HEI2015_TOT_PREG HEI2015_TOT_POST
			food_38 food_40;
	rename 	food_1 = Milk_drinks					food_2 = Milk_desserts					food_3 = Cheese
			food_4 = Poultry						food_5 = Fish_shellfish					food_6 = Meat
			food_7 = Cured_meat						food_8 = Eggs							food_9 = Legumes_nuts_seeds
			food_10 = Breads_refined				food_11 = Breads_whole					food_12 = Cakes_cookies
			food_13 = Pancakes_other				food_14 = Savory_pies					food_15 = Sandwiches
			food_16 = Patties_loaves				food_17 = Pasta_based					food_18 = Tortilla_based			 	
			food_19 = Mayonnaise_salads				food_20 = Pastas_rice					food_21 = Crackers_salty_snacks
			food_22 = Breakfast_cereals_high_sugar	food_23 = Breakfast_cereals_low_sugar	food_24 = Fruit_juice
			food_25 = Whole_fruits					food_26 = Potatoes						food_27 = Potatoes_fried
			food_28 = DG_vegetables					food_29 = RO_vegetables					food_30 = Other_vegetables
			food_31 = Solid_fats					food_32 = Oils							food_33 = Sauces
			food_34 = Soups							food_35 = Sugars_sweets					food_36 = Coffee_tea
			food_37 = Sugar_sweetened_drinks												food_39 = Water;
run;


/*
	rename 	food_1 = Milk_drinks					food_2 = Milk_desserts					food_3 = Cheese
			food_4 = Poultry						food_5 = Fish_shellfish					food_6 = Meat
			food_7 = Cured_meat						food_8 = Eggs							food_9 = Legumes_nuts_seeds
			food_10 = Breads_refined				food_11 = Breads_whole					food_12 = Cakes_cookies
			food_13 = Pancakes_other				food_14 = Savory_pies					food_15 = Sandwiches
			food_16 = Patties_loaves				food_17 = Pasta_based					food_18 = Tortilla_based			 	
			food_19 = Mayonnaise_salads				food_20 = Pastas_rice					food_21 = Crackers_salty_snacks
			food_22 = Breakfast_cereals_high_sugar	food_23 = Breakfast_cereals_low_sugar	food_24 = Fruit_juice
			food_25 = Whole_fruits					food_26 = Potatoes						food_27 = Potatoes_fried
			food_28 = DG_vegetables					food_29 = RO_vegetables					food_30 = Other_vegetables
			food_31 = Solid_fats					food_32 = Oils							food_33 = Sauces
			food_34 = Soups							food_35 = Sugars_sweets					food_36 = Coffee_tea
			food_37 = Sugar_sweetened_drinks		food_38 = Alcoholic_drinks				food_39 = Water
			food_40 = Nutritional_drinks;
*/