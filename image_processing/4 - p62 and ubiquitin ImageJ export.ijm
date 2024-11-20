//path=getDirectory("Choose a data folder"); print(path); exit;
input_path_13 = "/Volumes/ScotterLab/Evelyn Jade/Experiments/240903 - Exp 13 - Further Notch testing/"
input_folder_13 = "images_for_analysis_plate_a/"

input_path_14 = "/Volumes/ScotterLab/Evelyn Jade/Experiments/240910 - Exp 14 - p62 + ubiquitin testing/"
input_folder_14 = "images_for_analysis/"


output_path = "/Volumes/ScotterLab/Evelyn Jade/Experiments/241007 - Exp 16 - Microscopy images for thesis/4 - Aggregation associated proteins/"
output_folder = "imagej_export/"
output_format = "TIFF"



alexa_488 = newArray("40X/r06c04f03p01-ch2sk1fk1fl1.tiff", "40X/r06c09f05p01-ch2sk1fk1fl1.tiff",
					 "40X/r05c02f01p01-ch2sk1fk1fl1.tiff", "40X/r05c06f04p01-ch2sk1fk1fl1.tiff",
					 "40X/r06c02f03p01-ch2sk1fk1fl1.tiff", "40X/r06c06f04p01-ch2sk1fk1fl1.tiff",
					 "40X/r02c02f01p01-ch2sk1fk1fl1.tiff", "40X/r02c06f05p01-ch2sk1fk1fl1.tiff",
					 "40X/r04c02f01p01-ch2sk1fk1fl1.tiff", "40X/r04c06f05p01-ch2sk1fk1fl1.tiff",
					 "40X/r03c02f01p01-ch2sk1fk1fl1.tiff", "40X/r03c06f05p01-ch2sk1fk1fl1.tiff");

alexa_488_names = newArray("gpP62-Mock-488", "gpP62-GGC98-488",
						   "mP62-Mock-488", "mP62-GGC98-488",
						   "FLAG-p62-Mock-488", "FLAG-p62-GGC98-488",
						   "K48-Mock-488", "K48-GGC98-488",
						   "K63-Mock-488", "K63-GGC98-488",
						   "Pan-Mock-488", "Pan-GGC98-488");


alexa_594 = newArray("40X/r06c04f03p01-ch3sk1fk1fl1.tiff", "40X/r06c09f05p01-ch3sk1fk1fl1.tiff",
					 "40X/r05c02f01p01-ch3sk1fk1fl1.tiff", "40X/r05c06f04p01-ch3sk1fk1fl1.tiff",
					 "40X/r06c02f03p01-ch3sk1fk1fl1.tiff", "40X/r06c06f04p01-ch3sk1fk1fl1.tiff",
					 "40X/r02c02f01p01-ch3sk1fk1fl1.tiff", "40X/r02c06f05p01-ch3sk1fk1fl1.tiff",
					 "40X/r04c02f01p01-ch3sk1fk1fl1.tiff", "40X/r04c06f05p01-ch3sk1fk1fl1.tiff",
					 "40X/r03c02f01p01-ch3sk1fk1fl1.tiff", "40X/r03c06f05p01-ch3sk1fk1fl1.tiff");

alexa_594_names = newArray("gpP62-Mock-594", "gpP62-GGC98-594",
						   "mP62-Mock-594", "mP62-GGC98-594",
						   "FLAG-p62-Mock-594", "FLAG-p62-GGC98-594",
						   "K48-Mock-594", "K48-GGC98-594",
						   "K63-Mock-594", "K63-GGC98-594",
						   "Pan-Mock-594", "Pan-GGC98-594");


alexa_647 = newArray("40X/r06c04f03p01-ch4sk1fk1fl1.tiff", "40X/r06c09f05p01-ch4sk1fk1fl1.tiff",
					 "40X/r05c02f01p01-ch4sk1fk1fl1.tiff", "40X/r05c06f04p01-ch4sk1fk1fl1.tiff",
					 "40X/r06c02f03p01-ch4sk1fk1fl1.tiff", "40X/r06c06f04p01-ch4sk1fk1fl1.tiff",
					 "40X/r02c02f01p01-ch4sk1fk1fl1.tiff", "40X/r02c06f05p01-ch4sk1fk1fl1.tiff",
					 "40X/r04c02f01p01-ch4sk1fk1fl1.tiff", "40X/r04c06f05p01-ch4sk1fk1fl1.tiff",
					 "40X/r03c02f01p01-ch4sk1fk1fl1.tiff", "40X/r03c06f05p01-ch4sk1fk1fl1.tiff");

alexa_647_names = newArray("gpP62-Mock-647", "gpP62-GGC98-647",
						   "mP62-Mock-647", "mP62-GGC98-647",
						   "FLAG-p62-Mock-647", "FLAG-p62-GGC98-647",
						   "K48-Mock-647", "K48-GGC98-647",
						   "K63-Mock-647", "K63-GGC98-647",
						   "Pan-Mock-647", "Pan-GGC98-647");
					 
					 
hoechst = newArray("40X/r06c04f03p01-ch1sk1fk1fl1.tiff", "40X/r06c09f05p01-ch1sk1fk1fl1.tiff",
				   "40X/r05c02f01p01-ch1sk1fk1fl1.tiff", "40X/r05c06f04p01-ch1sk1fk1fl1.tiff",
				   "40X/r06c02f03p01-ch1sk1fk1fl1.tiff", "40X/r06c06f04p01-ch1sk1fk1fl1.tiff",
				   "40X/r02c02f01p01-ch1sk1fk1fl1.tiff", "40X/r02c06f05p01-ch1sk1fk1fl1.tiff",
				   "40X/r04c02f01p01-ch1sk1fk1fl1.tiff", "40X/r04c06f05p01-ch1sk1fk1fl1.tiff",
				   "40X/r03c02f01p01-ch1sk1fk1fl1.tiff", "40X/r03c06f05p01-ch1sk1fk1fl1.tiff");

hoechst_names = newArray("gpP62-Mock-H", "gpP62-GGC98-H",
						 "mP62-Mock-H", "mP62-GGC98-H",
						 "FLAG-p62-Mock-H", "FLAG-p62-GGC98-H",
						 "K48-Mock-H", "K48-GGC98-H",
						 "K63-Mock-H", "K63-GGC98-H",
						 "Pan-Mock-H", "Pan-GGC98-H");

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


// Make stack of Alexa 488
for (i=0; i<2; i++) {
	open(input_path_13 + input_folder_13 + alexa_488[i]);
	rename(alexa_488_names[i]);
}
run("Images to Stack", "name=[alexa_488_1] use");
setMinAndMax(500, 2000);
run("Apply LUT", "stack");

for (i=2; i<4; i++) {
	open(input_path_14 + input_folder_14 + alexa_488[i]);
	rename(alexa_488_names[i]);
}
run("Images to Stack", "name=[alexa_488_2] use");
setMinAndMax(100, 800);
run("Apply LUT", "stack");

for (i=4; i<8; i++) {
	open(input_path_14 + input_folder_14 + alexa_488[i]);
	rename(alexa_488_names[i]);
}
run("Images to Stack", "name=[alexa_488_3] use");
setMinAndMax(300, 3000);
run("Apply LUT", "stack");

for (i=8; i<10; i++) {
	open(input_path_14 + input_folder_14 + alexa_488[i]);
	rename(alexa_488_names[i]);
}
run("Images to Stack", "name=[alexa_488_4] use");
setMinAndMax(100, 2000);
run("Apply LUT", "stack");

for (i=10; i<12; i++) {
	open(input_path_14 + input_folder_14 + alexa_488[i]);
	rename(alexa_488_names[i]);
}
run("Images to Stack", "name=[alexa_488_5] use");
setMinAndMax(100, 2000);
run("Apply LUT", "stack");

run("Concatenate...", "title=alexa_488 image1=alexa_488_1 image2=alexa_488_2 image3=alexa_488_3 image4=alexa_488_4 image5=alexa_488_5");


// Export cropped regions
crop_and_save(2, 495, 65, alexa_488_names, false);
crop_and_save(4, 100, 690, alexa_488_names, false);
crop_and_save(6, 1050, 480, alexa_488_names, false);
crop_and_save(8, 250, 150, alexa_488_names, false);
crop_and_save(10, 65, 500, alexa_488_names, false);
crop_and_save(12, 250, 265, alexa_488_names, false);

// Export uncropped regions
for (i=1; i<=alexa_488.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+alexa_488_names[i-1]);
	close();
}




// Make stack of Alexa 594
for (i=0; i<2; i++) {
	open(input_path_13 + input_folder_13 + alexa_594[i]);
	rename(alexa_594_names[i]);
}
run("Images to Stack", "name=[alexa_594_1] use");

for (i=2; i<6; i++) {
	open(input_path_14 + input_folder_14 + alexa_594[i]);
	rename(alexa_594_names[i]);
}
run("Images to Stack", "name=[alexa_594_2] use");

for (i=6; i<8; i++) {
	open(input_path_14 + input_folder_14 + alexa_594[i]);
	rename(alexa_594_names[i]);
}
run("Images to Stack", "name=[alexa_594_3] use");

for (i=8; i<12; i++) {
	open(input_path_14 + input_folder_14 + alexa_594[i]);
	rename(alexa_594_names[i]);
}
run("Images to Stack", "name=[alexa_594_4] use");



run("Concatenate...", "title=alexa_594 image1=alexa_594_1 image2=alexa_594_2 image3=alexa_594_3 image4=alexa_594_4");
setMinAndMax(350, 2500);
run("Apply LUT", "stack");

// Export cropped regions
crop_and_save(2, 495, 65, alexa_594_names, false);
crop_and_save(4, 100, 690, alexa_594_names, false);
crop_and_save(6, 1050, 480, alexa_594_names, false);
crop_and_save(8, 250, 150, alexa_594_names, false);
crop_and_save(10, 65, 500, alexa_594_names, false);
crop_and_save(12, 250, 265, alexa_594_names, false);


// Export uncropped regions
for (i=1; i<=alexa_594.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+alexa_594_names[i-1]);
	close();
}



// Make stack of Alexa 647
for (i=0; i<2; i++) {
	open(input_path_13 + input_folder_13 + alexa_647[i]);
	rename(alexa_647_names[i]);
}
run("Images to Stack", "name=[alexa_647_1] use");

for (i=2; i<12; i++) {
	open(input_path_14 + input_folder_14 + alexa_647[i]);
	rename(alexa_647_names[i]);
}
run("Images to Stack", "name=[alexa_647_2] use");

run("Concatenate...", "title=alexa_647 image1=alexa_647_1 image2=alexa_647_2");
setMinAndMax(400, 4000);
run("Apply LUT", "stack");


// Export cropped regions
crop_and_save(2, 495, 65, alexa_647_names, false);
crop_and_save(4, 100, 690, alexa_647_names, false);
crop_and_save(6, 1050, 480, alexa_647_names, false);
crop_and_save(8, 250, 150, alexa_647_names, false);
crop_and_save(10, 65, 500, alexa_647_names, false);
crop_and_save(12, 250, 265, alexa_647_names, false);

// Export uncropped regions
for (i=1; i<=alexa_647.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+alexa_647_names[i-1]);
	close();
}



// Make stack of Hoechst
for (i=0; i<2; i++) {
	open(input_path_13 + input_folder_13 + hoechst[i]);
	rename(hoechst_names[i]);
}
for (i=2; i<hoechst.length; i++) {
	open(input_path_14 + input_folder_14 + hoechst[i]);
	rename(hoechst_names[i]);
}

run("Images to Stack", "name=[hoechst] use");
setMinAndMax(500, 10000);
run("Apply LUT", "stack");

// Export cropped regions
crop_and_save(2, 495, 65, hoechst_names, false);
crop_and_save(4, 100, 690, hoechst_names, false);
crop_and_save(6, 1050, 480, hoechst_names, false);
crop_and_save(8, 250, 150, hoechst_names, false);
crop_and_save(10, 65, 500, hoechst_names, false);
crop_and_save(12, 250, 265, hoechst_names, false);

// Export uncropped regions
for (i=1; i<=hoechst.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+hoechst_names[i-1]);
	close();
}


run("Merge Channels...", "c1=alexa_594 c2=alexa_488 c3=hoechst c5=alexa_647 create");
exit

close("*");
