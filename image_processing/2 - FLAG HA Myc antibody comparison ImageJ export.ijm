//path=getDirectory("Choose a data folder"); print(path); exit;
input_path_9 = "/Volumes/ScotterLab/Evelyn Jade/Experiments/240729 - Exp 9 - FLAG:HA:Myc transfection attempt 2/"
input_path_10 = "/Volumes/ScotterLab/Evelyn Jade/Experiments/240819 - Exp 10 - Myc re-testing/"
input_folder = "images_for_analysis/"

output_path = "/Volumes/ScotterLab/Evelyn Jade/Experiments/241007 - Exp 16 - Microscopy images for thesis/2 - FLAG HA Myc antibody comparison/"
output_folder = "imagej_export/"
output_format = "TIFF"


// Hoescht = Ch 2
// Alexa 488 = Ch 1

// FLAG = Exp 9 R3 C4 F3
// HA = Exp 9 R4,5 C6-7       Goat: 46-4, Mouse = 56-2
// Myc = Exp 10 R4,5 C5 F3


hoescht = newArray("003002-3-001001002.tif", "003004-3-001001002.tif", 
				   "004002-4-001001002.tif", "004006-4-001001002.tif", 
				   "005002-1-001001002.tif", "005007-4-001001002.tif", 
				   "004003-3-001001002.tif", "004005-3-001001002.tif", 
				   "005003-3-001001002.tif", "005005-3-001001002.tif");

hoescht_names = newArray("mFLAG-Mock-H", "mFLAG-Positive-H",
						 "gHA-Mock-H", "gHA-Positive-H",
						 "mHA-Mock-H", "mHA-Positive-H",
						 "9E10-Mock-H", "9E10-Positive-H",
						 "9E11-Mock-H", "9E11-Positive-H");


alexa_488 = newArray("003002-3-001001001.tif", "003004-3-001001001.tif", 
					 "004002-4-001001001.tif", "004006-4-001001001.tif", 
					 "005002-1-001001001.tif", "005007-4-001001001.tif", 
					 "004003-3-001001001.tif", "004005-3-001001001.tif", 
					 "005003-3-001001001.tif", "005005-3-001001001.tif");

alexa_488_names = newArray("mFLAG-Mock-G", "mFLAG-Positive-G",
						   "gHA-Mock-G", "gHA-Positive-G",
						   "mHA-Mock-G", "mHA-Positive-G",
						   "9E10-Mock-G", "9E10-Positive-G",
						   "9E11-Mock-G", "9E11-Positive-G");


// Make stack of Hoescht
for (i=0; i<6; i++) {
	open(input_path_9 + input_folder + hoescht[i]);
	rename(hoescht_names[i]);
}
run("Images to Stack", "name=[hoescht_1] use");
setMinAndMax(2000, 10000);
run("Apply LUT", "stack");

for (i=6; i<hoescht.length; i++) {
	open(input_path_10 + input_folder + hoescht[i]);
	rename(hoescht_names[i]);
}
run("Images to Stack", "name=[hoescht_2] use");
setMinAndMax(1000, 15000);
run("Apply LUT", "stack");

run("Concatenate...", "title=hoescht image1=hoescht_1 image2=hoescht_2");

for (i=1; i<=hoescht.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+hoescht_names[i-1]);
	close();
}



// Make stack of Alexa 488
for (i=0; i<2; i++) {
	open(input_path_9 + input_folder + alexa_488[i]);
	rename(alexa_488_names[i]);
}
run("Images to Stack", "name=[alexa_488_1] use");
setMinAndMax(1000, 17500);
run("Apply LUT", "stack");

for (i=2; i<6; i++) {
	open(input_path_9 + input_folder + alexa_488[i]);
	rename(alexa_488_names[i]);
}
run("Images to Stack", "name=[alexa_488_2] use");
setMinAndMax(150, 3000);
run("Apply LUT", "stack");

for (i=6; i<alexa_488.length; i++) {
	open(input_path_10 + input_folder + alexa_488[i]);
	rename(alexa_488_names[i]);
}
run("Images to Stack", "name=[alexa_488_3] use");
setMinAndMax(250, 2500);
run("Apply LUT", "stack");

run("Concatenate...", "title=alexa_488 image1=alexa_488_1 image2=alexa_488_2 image3=alexa_488_3");

for (i=1; i<=alexa_488.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+alexa_488_names[i-1]);
	close();
}

close("*");
