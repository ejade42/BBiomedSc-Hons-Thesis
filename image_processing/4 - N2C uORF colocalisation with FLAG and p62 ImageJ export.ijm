//path=getDirectory("Choose a data folder"); print(path); exit;
input_path = "/Volumes/ScotterLab/Evelyn Jade/Experiments/240903 - Exp 13 - Further Notch testing/"
input_folder = "images_for_analysis_plate_a/"

output_path = "/Volumes/ScotterLab/Evelyn Jade/Experiments/241007 - Exp 16 - Microscopy images for thesis/4 - NOTCH2NLC colocalisation with FLAG and p62/"
output_folder = "imagej_export/"
output_format = "TIFF"



green = newArray("40X/r02c04f02p01-ch2sk1fk1fl1.tiff", "40X/r02c08f02p01-ch2sk1fk1fl1.tiff",
			     "40X/r04c04f02p01-ch2sk1fk1fl1.tiff", "40X/r04c08f05p01-ch2sk1fk1fl1.tiff",
			     "40X/r05c04f04p01-ch2sk1fk1fl1.tiff", "40X/r05c08f02p01-ch2sk1fk1fl1.tiff",
			     "40X/r06c04f03p01-ch2sk1fk1fl1.tiff", "40X/r06c09f05p01-ch2sk1fk1fl1.tiff");

green_names = newArray("FLAG-Mock-G", "FLAG-GGC98-G",
					   "4C4-Mock-G", "4C4-GGC98-G",
					   "4D12-Mock-G", "4D12-GGC98-G",
					   "Rabbit-Mock-G", "Rabbit-GGC98-G");


red = newArray("40X/r02c04f02p01-ch3sk1fk1fl1.tiff", "40X/r02c08f02p01-ch3sk1fk1fl1.tiff",
			   "40X/r04c04f02p01-ch3sk1fk1fl1.tiff", "40X/r04c08f05p01-ch3sk1fk1fl1.tiff",
			   "40X/r05c04f04p01-ch3sk1fk1fl1.tiff", "40X/r05c08f02p01-ch3sk1fk1fl1.tiff",
			   "40X/r06c04f03p01-ch3sk1fk1fl1.tiff", "40X/r06c09f05p01-ch3sk1fk1fl1.tiff");

red_names = newArray("FLAG-Mock-R", "FLAG-GGC98-R",
					 "4C4-Mock-R", "4C4-GGC98-R",
					 "4D12-Mock-R", "4D12-GGC98-R",
					 "Rabbit-Mock-R", "Rabbit-GGC98-R");


p62 = newArray("40X/r02c04f02p01-ch4sk1fk1fl1.tiff", "40X/r02c08f02p01-ch4sk1fk1fl1.tiff",
			   "40X/r04c04f02p01-ch4sk1fk1fl1.tiff", "40X/r04c08f05p01-ch4sk1fk1fl1.tiff",
			   "40X/r05c04f04p01-ch4sk1fk1fl1.tiff", "40X/r05c08f02p01-ch4sk1fk1fl1.tiff",
			   "40X/r06c04f03p01-ch4sk1fk1fl1.tiff", "40X/r06c09f05p01-ch4sk1fk1fl1.tiff");

p62_names = newArray("FLAG-Mock-P", "FLAG-GGC98-P",
					 "4C4-Mock-P", "4C4-GGC98-P",
					 "4D12-Mock-P", "4D12-GGC98-P",
					 "Rabbit-Mock-P", "Rabbit-GGC98-P");
					 
					 
hoechst = newArray("40X/r02c04f02p01-ch1sk1fk1fl1.tiff", "40X/r02c08f02p01-ch1sk1fk1fl1.tiff",
			       "40X/r04c04f02p01-ch1sk1fk1fl1.tiff", "40X/r04c08f05p01-ch1sk1fk1fl1.tiff",
				   "40X/r05c04f04p01-ch1sk1fk1fl1.tiff", "40X/r05c08f02p01-ch1sk1fk1fl1.tiff",
				   "40X/r06c04f03p01-ch1sk1fk1fl1.tiff", "40X/r06c09f05p01-ch1sk1fk1fl1.tiff");

hoechst_names = newArray("FLAG-Mock-H", "FLAG-GGC98-H",
					     "4C4-Mock-H", "4C4-GGC98-H",
					     "4D12-Mock-H", "4D12-GGC98-H",
					     "Rabbit-Mock-H", "Rabbit-GGC98-H");

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


// Make stack of green
for (i=0; i<2; i++) {
	open(input_path + input_folder + green[i]);
	rename(green_names[i]);
}
run("Images to Stack", "name=[green_flag] use");
setMinAndMax(100, 10000);
run("Apply LUT", "stack");

for (i=2; i<8; i++) {
	open(input_path + input_folder + green[i]);
	rename(green_names[i]);
}
run("Images to Stack", "name=[green_notch] use");
setMinAndMax(250, 1250);
run("Apply LUT", "stack");

run("Concatenate...", "title=green image1=green_flag image2=green_notch");

// Export cropped regions
crop_and_save(2, 478, 210, green_names, false);
crop_and_save(4, 1000, 330, green_names, false);
crop_and_save(6, 125, 35, green_names, false);
crop_and_save(8, 550, 500, green_names, false);

// Export uncropped regions
for (i=1; i<=green.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+green_names[i-1]);
	close();
}




// Make stack of red
for (i=0; i<red.length; i++) {
	open(input_path + input_folder + red[i]);
	rename(red_names[i]);
}

run("Images to Stack", "name=[red] use");
setMinAndMax(250, 750);
run("Apply LUT", "stack");

// Export cropped regions
crop_and_save(2, 478, 210, red_names, false);
crop_and_save(4, 1000, 330, red_names, false);
crop_and_save(6, 125, 35, red_names, false);
crop_and_save(8, 550, 500, red_names, false);

// Export uncropped regions
for (i=1; i<=red.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+red_names[i-1]);
	close();
}


// Make stack of p62
for (i=0; i<p62.length; i++) {
	open(input_path + input_folder + p62[i]);
	rename(p62_names[i]);
}
run("Images to Stack", "name=[p62] use");
setMinAndMax(100, 10000);
run("Apply LUT", "stack");

// Export cropped regions
crop_and_save(2, 478, 210, p62_names, false);
crop_and_save(4, 1000, 330, p62_names, false);
crop_and_save(6, 125, 35, p62_names, false);
crop_and_save(8, 550, 500, p62_names, false);

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
setMinAndMax(500, 10000);
run("Apply LUT", "stack");

// Export cropped regions
crop_and_save(2, 478, 210, hoechst_names, false);
crop_and_save(4, 1000, 330, hoechst_names, false);
crop_and_save(6, 125, 35, hoechst_names, false);
crop_and_save(8, 550, 500, hoechst_names, false);

// Export uncropped regions
for (i=1; i<=hoechst.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+hoechst_names[i-1]);
	close();
}


run("Merge Channels...", "c1=red c2=green c3=p62 c4=hoechst create");
exit

close("*");
