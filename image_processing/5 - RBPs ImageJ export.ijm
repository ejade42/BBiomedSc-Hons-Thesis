//path=getDirectory("Choose a data folder"); print(path); exit;
input_path = "/Volumes/ScotterLab/Evelyn Jade/Experiments/240903 - Exp 13B - RBPs/"
input_folder = "images_for_analysis/"

output_path = "/Volumes/ScotterLab/Evelyn Jade/Experiments/241007 - Exp 16 - Microscopy images for thesis/5 - RBPs/"
output_folder = "imagej_export/"
output_format = "TIFF"



flag = newArray("40X/r02c05f05p01-ch2sk1fk1fl1.tiff", "40X/r02c08f01p01-ch2sk1fk1fl1.tiff",
				"40X/r04c04f01p01-ch2sk1fk1fl1.tiff", "40X/r04c08f04p01-ch2sk1fk1fl1.tiff",
				"40X/r05c04f01p01-ch2sk1fk1fl1.tiff", "40X/r05c09f01p01-ch2sk1fk1fl1.tiff",
				"40X/r06c04f01p01-ch2sk1fk1fl1.tiff", "40X/r06c08f02p01-ch2sk1fk1fl1.tiff",
				"40X/r07c04f01p01-ch2sk1fk1fl1.tiff", "40X/r07c08f05p01-ch2sk1fk1fl1.tiff");
		
flag_names = newArray("TDP43-Mock-F", "TDP43-GGC98-F",
					  "Fus-Mock-F", "Fus-GGC98-F",
					  "Matrin-Mock-F", "Matrin-GGC98-F",
					  "hnRNP-Mock-F", "hnRNP-GGC98-F",
					  "ATXN2-Mock-F", "ATXN2-GGC98-F");
					  

RBPs = newArray("40X/r02c05f05p01-ch3sk1fk1fl1.tiff", "40X/r02c08f01p01-ch3sk1fk1fl1.tiff",
				"40X/r04c04f01p01-ch3sk1fk1fl1.tiff", "40X/r04c08f04p01-ch3sk1fk1fl1.tiff",
				"40X/r05c04f01p01-ch3sk1fk1fl1.tiff", "40X/r05c09f01p01-ch3sk1fk1fl1.tiff",
				"40X/r06c04f01p01-ch3sk1fk1fl1.tiff", "40X/r06c08f02p01-ch3sk1fk1fl1.tiff",
				"40X/r07c04f01p01-ch3sk1fk1fl1.tiff", "40X/r07c08f05p01-ch3sk1fk1fl1.tiff");
		
RBPs_names = newArray("TDP43-Mock-R", "TDP43-GGC98-R",
					  "Fus-Mock-R", "Fus-GGC98-R",
					  "Matrin-Mock-R", "Matrin-GGC98-R",
					  "hnRNP-Mock-R", "hnRNP-GGC98-R",
					  "ATXN2-Mock-R", "ATXN2-GGC98-R");



p62 = newArray("40X/r02c05f05p01-ch4sk1fk1fl1.tiff", "40X/r02c08f01p01-ch4sk1fk1fl1.tiff",
			   "40X/r04c04f01p01-ch4sk1fk1fl1.tiff", "40X/r04c08f04p01-ch4sk1fk1fl1.tiff",
			   "40X/r05c04f01p01-ch4sk1fk1fl1.tiff", "40X/r05c09f01p01-ch4sk1fk1fl1.tiff",
			   "40X/r06c04f01p01-ch4sk1fk1fl1.tiff", "40X/r06c08f02p01-ch4sk1fk1fl1.tiff",
			   "40X/r07c04f01p01-ch4sk1fk1fl1.tiff", "40X/r07c08f05p01-ch4sk1fk1fl1.tiff");
		
p62_names = newArray("TDP43-Mock-P", "TDP43-GGC98-P",
					 "Fus-Mock-P", "Fus-GGC98-P",
					 "Matrin-Mock-P", "Matrin-GGC98-P",
					 "hnRNP-Mock-P", "hnRNP-GGC98-P",
					 "ATXN2-Mock-P", "ATXN2-GGC98-P");
					 
					 
hoechst = newArray("40X/r02c05f05p01-ch1sk1fk1fl1.tiff", "40X/r02c08f01p01-ch1sk1fk1fl1.tiff",
				   "40X/r04c04f01p01-ch1sk1fk1fl1.tiff", "40X/r04c08f04p01-ch1sk1fk1fl1.tiff",
				   "40X/r05c04f01p01-ch1sk1fk1fl1.tiff", "40X/r05c09f01p01-ch1sk1fk1fl1.tiff",
				   "40X/r06c04f01p01-ch1sk1fk1fl1.tiff", "40X/r06c08f02p01-ch1sk1fk1fl1.tiff",
				   "40X/r07c04f01p01-ch1sk1fk1fl1.tiff", "40X/r07c08f05p01-ch1sk1fk1fl1.tiff");
		
hoechst_names = newArray("TDP43-Mock-H", "TDP43-GGC98-H",
					     "Fus-Mock-H", "Fus-GGC98-H",
					     "Matrin-Mock-H", "Matrin-GGC98-H",
					     "hnRNP-Mock-H", "hnRNP-GGC98-H",
					     "ATXN2-Mock-H", "ATXN2-GGC98-H");

// Function for cropping
function crop_and_save(slice_number, left_x, top_y, names_list, do_exit) {
	setSlice(slice_number);
	run("Duplicate...", "title=temp");
	makeRectangle(left_x, top_y, 280, 210);
	if (do_exit == true) {exit};
	run("Crop");
	saveAs(output_format, output_path+output_folder+"cropped_"+names_list[slice_number-1]);
	close();
}


// Make stack of FLAG
for (i=0; i<flag.length; i++) {
	open(input_path + input_folder + flag[i]);
	rename(flag_names[i]);
}

run("Images to Stack", "name=[flag] use");
setMinAndMax(100, 5000);
run("Apply LUT", "stack");

// Export cropped regions
crop_and_save(2, 618, 442, flag_names, false);
crop_and_save(4, 650, 350, flag_names, false);
crop_and_save(6, 250, 590, flag_names, false);
crop_and_save(8, 750, 365, flag_names, false);
crop_and_save(10, 700, 575, flag_names, false);

// Export uncropped regions
for (i=1; i<=flag.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+flag_names[i-1]);
	close();
}




// Make stack of RBPs - different levels
for (i=0; i<2; i++) {
	open(input_path + input_folder + RBPs[i]);
	rename(RBPs_names[i]);
}
run("Images to Stack", "name=[RBPs_1] use");
setMinAndMax(50, 200);
run("Apply LUT", "stack");

for (i=2; i<4; i++) {
	open(input_path + input_folder + RBPs[i]);
	rename(RBPs_names[i]);
}
run("Images to Stack", "name=[RBPs_2] use");
setMinAndMax(50, 750);
run("Apply LUT", "stack");

for (i=4; i<6; i++) {
	open(input_path + input_folder + RBPs[i]);
	rename(RBPs_names[i]);
}
run("Images to Stack", "name=[RBPs_3] use");
setMinAndMax(50, 1250);
run("Apply LUT", "stack");

for (i=6; i<10; i++) {
	open(input_path + input_folder + RBPs[i]);
	rename(RBPs_names[i]);
}
run("Images to Stack", "name=[RBPs_4] use");
setMinAndMax(50, 750);
run("Apply LUT", "stack");

run("Concatenate...", "title=RBPs image1=RBPs_1 image2=RBPs_2 image3=RBPs_3 image4=RBPs_4");

// Export cropped regions
crop_and_save(2, 618, 442, RBPs_names, false);
crop_and_save(4, 650, 350, RBPs_names, false);
crop_and_save(6, 250, 590, RBPs_names, false);
crop_and_save(8, 750, 365, RBPs_names, false);
crop_and_save(10, 700, 575, RBPs_names, false);

// Export uncropped regions
for (i=1; i<=RBPs.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+RBPs_names[i-1]);
	close();
}




// Make stack of p62
for (i=0; i<p62.length; i++) {
	open(input_path + input_folder + p62[i]);
	rename(p62_names[i]);
}
run("Images to Stack", "name=[p62] use");
setMinAndMax(50, 500);
run("Apply LUT", "stack");

// Export cropped regions
crop_and_save(2, 618, 442, p62_names, false);
crop_and_save(4, 650, 350, p62_names, false);
crop_and_save(6, 250, 590, p62_names, false);
crop_and_save(8, 750, 365, p62_names, false);
crop_and_save(10, 700, 575, p62_names, false);

// Export uncropped regions
for (i=1; i<=p62.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+p62_names[i-1]);
	close();
}



// Make stack of Hoechst
for (i=0; i<hoechst.length; i++) {
	open(input_path + input_folder + hoechst[i]);
	rename(hoechst_names[i]);
}

run("Images to Stack", "name=[hoechst] use");
setMinAndMax(300, 10000);
run("Apply LUT", "stack");

// Export cropped regions
crop_and_save(2, 618, 442, hoechst_names, false);
crop_and_save(4, 650, 350, hoechst_names, false);
crop_and_save(6, 250, 590, hoechst_names, false);
crop_and_save(8, 750, 365, hoechst_names, false);
crop_and_save(10, 700, 575, hoechst_names, false);

// Export uncropped regions
for (i=1; i<=hoechst.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+hoechst_names[i-1]);
	close();
}


run("Merge Channels...", "c1=RBPs c2=flag c3=hoechst c5=p62 create");
exit

close("*");
