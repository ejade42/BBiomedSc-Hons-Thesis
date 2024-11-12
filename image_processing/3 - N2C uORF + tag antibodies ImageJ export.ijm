//path=getDirectory("Choose a data folder"); print(path); exit;
input_path = "/Volumes/ScotterLab/Evelyn Jade/Experiments/240902 - Exp 12 - NOTCH2NLC antibody testing take 2/"
input_folder = "images_for_analysis/"

output_path = "/Volumes/ScotterLab/Evelyn Jade/Experiments/241007 - Exp 16 - Microscopy images for thesis/3 - NOTCH2NLC + tag antibodies with NOTCH plasmid/"
output_folder = "imagej_export/"
output_format = "TIFF"



alexa_488 = newArray("40X/r02c02f01p01-ch1sk1fk1fl1.tiff", "40X/r02c04f01p01-ch1sk1fk1fl1.tiff",
					 "40X/r03c02f02p01-ch1sk1fk1fl1.tiff", "40X/r03c04f01p01-ch1sk1fk1fl1.tiff",
					 "40X/r04c02f03p01-ch1sk1fk1fl1.tiff", "40X/r04c04f03p01-ch1sk1fk1fl1.tiff",
					 "40X/r05c02f01p01-ch1sk1fk1fl1.tiff", "40X/r05c04f01p01-ch1sk1fk1fl1.tiff",
					 "40X/r06c02f01p01-ch1sk1fk1fl1.tiff", "40X/r06c04f01p01-ch1sk1fk1fl1.tiff",
					 "40X/r07c02f01p01-ch1sk1fk1fl1.tiff", "40X/r07c04f01p01-ch1sk1fk1fl1.tiff");

alexa_488_names = newArray("4C4-Mock-G", "4C4-GGC98-G",
					 	   "4D12-Mock-G", "4D12-GGC98-G",
					 	   "Rabbit-Mock-G", "Rabbit-GGC98-G",
					 	   "FLAG-Mock-G", "FLAG-GGC98-G",
					 	   "HA-Mock-G", "HA-GGC98-G",
					 	   "Myc-Mock-G", "Myc-GGC98-G");


hoescht = newArray("40X/r02c02f01p01-ch2sk1fk1fl1.tiff", "40X/r02c04f01p01-ch2sk1fk1fl1.tiff",
				   "40X/r03c02f02p01-ch2sk1fk1fl1.tiff", "40X/r03c04f01p01-ch2sk1fk1fl1.tiff",
				   "40X/r04c02f03p01-ch2sk1fk1fl1.tiff", "40X/r04c04f03p01-ch2sk1fk1fl1.tiff",
				   "40X/r05c02f01p01-ch2sk1fk1fl1.tiff", "40X/r05c04f01p01-ch2sk1fk1fl1.tiff",
				   "40X/r06c02f01p01-ch2sk1fk1fl1.tiff", "40X/r06c04f01p01-ch2sk1fk1fl1.tiff",
				   "40X/r07c02f01p01-ch2sk1fk1fl1.tiff", "40X/r07c04f01p01-ch2sk1fk1fl1.tiff");

hoescht_names = newArray("4C4-Mock-H", "4C4-GGC98-H",
					 	 "4D12-Mock-H", "4D12-GGC98-H",
					 	 "Rabbit-Mock-H", "Rabbit-GGC98-H",
					 	 "FLAG-Mock-H", "FLAG-GGC98-H",
					 	 "HA-Mock-H", "HA-GGC98-H",
					 	 "Myc-Mock-H", "Myc-GGC98-H");


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
for (i=0; i<6; i++) {
	open(input_path + input_folder + alexa_488[i]);
	rename(alexa_488_names[i]);
}
run("Images to Stack", "name=[alexa_488_notch] use");
setMinAndMax(50, 400);
run("Apply LUT", "stack");

for (i=6; i<12; i++) {
	open(input_path + input_folder + alexa_488[i]);
	rename(alexa_488_names[i]);
}
run("Images to Stack", "name=[alexa_488_tags] use");
setMinAndMax(50, 600);
run("Apply LUT", "stack");

run("Concatenate...", "title=alexa_488 image1=alexa_488_notch image2=alexa_488_tags");

// Export cropped regions
crop_and_save(2, 575, 675, alexa_488_names, false);
crop_and_save(4, 622, 550, alexa_488_names, false);
crop_and_save(6, 450, 750, alexa_488_names, false);
crop_and_save(8, 325, 525, alexa_488_names, false);
crop_and_save(10, 975, 350, alexa_488_names, false);
crop_and_save(12, 540, 407, alexa_488_names, false);

// Export uncropped regions
for (i=1; i<=alexa_488.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+alexa_488_names[i-1]);
	close();
}



// Make stack of Hoescht
for (i=0; i<hoescht.length; i++) {
	open(input_path + input_folder + hoescht[i]);
	rename(hoescht_names[i]);
}

run("Images to Stack", "name=[hoescht] use");
setMinAndMax(1000, 10000);
run("Apply LUT", "stack");

// Export cropped regions
crop_and_save(2, 575, 675, hoescht_names, false);
crop_and_save(4, 622, 550, hoescht_names, false);
crop_and_save(6, 450, 750, hoescht_names, false);
crop_and_save(8, 325, 525, hoescht_names, false);
crop_and_save(10, 975, 350, hoescht_names, false);
crop_and_save(12, 540, 407, hoescht_names, false);

// Export uncropped regions
for (i=1; i<=hoescht.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+hoescht_names[i-1]);
	close();
}

close("*");
