//path=getDirectory("Choose a data folder"); print(path); exit;
input_path = "/Volumes/ScotterLab/Evelyn Jade/Experiments/240910 - Exp 14 - p62 + ubiquitin testing/"
input_folder = "images_for_analysis/"

output_path = "/Volumes/ScotterLab/Evelyn Jade/Experiments/241007 - Exp 16 - Microscopy images for thesis/6 - Ubiquitin and p62 further testing/"
output_folder = "imagej_export/"
output_format = "TIFF"

// Consider FLAG green, NOTCH red, Ub/p62 yellow?

hoechst = newArray("40X/r02c02f01p01-ch1sk1fk1fl1.tiff", "40X/r02c06f05p01-ch1sk1fk1fl1.tiff",
			       "40X/r03c02f01p01-ch1sk1fk1fl1.tiff", "40X/r03c06f04p01-ch1sk1fk1fl1.tiff",  //consider field 4 or 5
				   "40X/r04c02f01p01-ch1sk1fk1fl1.tiff", "40X/r04c06f05p01-ch1sk1fk1fl1.tiff",
				   "40X/r05c02f01p01-ch1sk1fk1fl1.tiff", "40X/r05c06f04p01-ch1sk1fk1fl1.tiff",
				   "40X/r06c02f03p01-ch1sk1fk1fl1.tiff", "40X/r06c06f04p01-ch1sk1fk1fl1.tiff");

hoechst_names = newArray("K48-Mock-H", "K48-GGC98-H",
					     "Ub-Mock-H", "Ub-GGC98-H",
					     "K63-Mock-H", "K63-GGC98-H",
					     "mP62-Mock-H", "mP62-GGC98-H",
					     "gpP62-Mock-H", "gpP62-GGC98-H");


green = newArray("40X/r02c02f01p01-ch2sk1fk1fl1.tiff", "40X/r02c06f05p01-ch2sk1fk1fl1.tiff",
			     "40X/r03c02f01p01-ch2sk1fk1fl1.tiff", "40X/r03c06f04p01-ch2sk1fk1fl1.tiff",  //consider field 4 or 5
				 "40X/r04c02f01p01-ch2sk1fk1fl1.tiff", "40X/r04c06f05p01-ch2sk1fk1fl1.tiff",
				 "40X/r05c02f01p01-ch2sk1fk1fl1.tiff", "40X/r05c06f04p01-ch2sk1fk1fl1.tiff",
				 "40X/r06c02f03p01-ch2sk1fk1fl1.tiff", "40X/r06c06f04p01-ch2sk1fk1fl1.tiff");

green_names = newArray("K48-Mock-G", "K48-GGC98-G",
					   "Ub-Mock-G", "Ub-GGC98-G",
					   "K63-Mock-G", "K63-GGC98-G",
					   "mP62-Mock-G", "mP62-GGC98-G",
					   "gpP62-Mock-G", "gpP62-GGC98-G");
					   
					   
red = newArray("40X/r02c02f01p01-ch3sk1fk1fl1.tiff", "40X/r02c06f05p01-ch3sk1fk1fl1.tiff",
			   "40X/r03c02f01p01-ch3sk1fk1fl1.tiff", "40X/r03c06f04p01-ch3sk1fk1fl1.tiff",  //consider field 4 or 5
			   "40X/r04c02f01p01-ch3sk1fk1fl1.tiff", "40X/r04c06f05p01-ch3sk1fk1fl1.tiff",
			   "40X/r05c02f01p01-ch3sk1fk1fl1.tiff", "40X/r05c06f04p01-ch3sk1fk1fl1.tiff",
			   "40X/r06c02f03p01-ch3sk1fk1fl1.tiff", "40X/r06c06f04p01-ch3sk1fk1fl1.tiff");

red_names = newArray("K48-Mock-R", "K48-GGC98-R",
					 "Ub-Mock-R", "Ub-GGC98-R",
					 "K63-Mock-R", "K63-GGC98-R",
					 "mP62-Mock-R", "mP62-GGC98-R",
					 "gpP62-Mock-R", "gpP62-GGC98-R");
					 
					 
guineapig = newArray("40X/r02c02f01p01-ch4sk1fk1fl1.tiff", "40X/r02c06f05p01-ch4sk1fk1fl1.tiff",
			  		 "40X/r03c02f01p01-ch4sk1fk1fl1.tiff", "40X/r03c06f04p01-ch4sk1fk1fl1.tiff",  //consider field 4 or 5
			   		 "40X/r04c02f01p01-ch4sk1fk1fl1.tiff", "40X/r04c06f05p01-ch4sk1fk1fl1.tiff",
			   		 "40X/r05c02f01p01-ch4sk1fk1fl1.tiff", "40X/r05c06f04p01-ch4sk1fk1fl1.tiff",
			   		 "40X/r06c02f03p01-ch4sk1fk1fl1.tiff", "40X/r06c06f04p01-ch4sk1fk1fl1.tiff");

guineapig_names = newArray("K48-Mock-647", "K48-GGC98-647",
					 	   "Ub-Mock-647", "Ub-GGC98-647",
					 	   "K63-Mock-647", "K63-GGC98-647",
					 	   "mP62-Mock-647", "mP62-GGC98-647",
					 	   "gpP62-Mock-647", "gpP62-GGC98-647");

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
for (i=0; i<green.length; i++) {
	open(input_path + input_folder + green[i]);
	rename(green_names[i]);
}
run("Images to Stack", "name=[green] use");
setMinAndMax(250, 4000);
run("Apply LUT", "stack");

// Export cropped regions
crop_and_save(2, 300, 575, green_names, false);
crop_and_save(4, 600, 225, green_names, false);
crop_and_save(6, 250, 700, green_names, false);
crop_and_save(8, 125, 807, green_names, false);
crop_and_save(10, 1050, 482, green_names, false);

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
setMinAndMax(250, 2500);
run("Apply LUT", "stack");


// Export cropped regions
crop_and_save(2, 300, 575, red_names, false);
crop_and_save(4, 600, 225, red_names, false);
crop_and_save(6, 250, 700, red_names, false);
crop_and_save(8, 125, 807, red_names, false);
crop_and_save(10, 1050, 482, red_names, false);

// Export uncropped regions
for (i=1; i<=red.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+red_names[i-1]);
	close();
}




// Make stack of alexa 647
for (i=0; i<guineapig.length; i++) {
	open(input_path + input_folder + guineapig[i]);
	rename(guineapig_names[i]);
}
run("Images to Stack", "name=[647] use");
setMinAndMax(250, 2500);
run("Apply LUT", "stack");

// Export cropped regions
crop_and_save(2, 300, 575, guineapig_names, false);
crop_and_save(4, 600, 225, guineapig_names, false);
crop_and_save(6, 250, 700, guineapig_names, false);
crop_and_save(8, 125, 807, guineapig_names, false);
crop_and_save(10, 1050, 482, guineapig_names, false);

// Export uncropped regions
for (i=1; i<=guineapig.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+guineapig_names[i-1]);
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
crop_and_save(2, 300, 575, hoechst_names, false);
crop_and_save(4, 600, 225, hoechst_names, false);
crop_and_save(6, 250, 700, hoechst_names, false);
crop_and_save(8, 125, 807, hoechst_names, false);
crop_and_save(10, 1050, 482, hoechst_names, false);

// Export uncropped regions
for (i=1; i<=hoechst.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+hoechst_names[i-1]);
	close();
}

run("Merge Channels...", "c1=red c2=green c3=hoechst c5=647 create");
exit

close("*");
