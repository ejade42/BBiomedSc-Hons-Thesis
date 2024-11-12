//path=getDirectory("Choose a data folder"); print(path); exit;
input_path = "/Volumes/ScotterLab/Evelyn Jade/Experiments/240722 - Exp 8 - EGFP-TDP-43-WT test transfection/"
input_folder = "images_for_analysis_v2/"

output_path = "/Volumes/ScotterLab/Evelyn Jade/Experiments/241007 - Exp 16 - Microscopy images for thesis/1 - Transfection optimisation/"
output_folder = "imagej_export/"
output_format = "TIFF"



gfp = newArray("plate_1/r03c03f03p01-ch2sk1fk1fl1.tiff",
			   "plate_1/r03c04f03p01-ch2sk1fk1fl1.tiff",
			   "plate_2/r03c05f03p01-ch2sk1fk1fl1.tiff",
			   "plate_2/r03c06f03p01-ch2sk1fk1fl1.tiff",
			   "plate_2/r03c03f03p01-ch2sk1fk1fl1.tiff",
			   "plate_2/r03c04f03p01-ch2sk1fk1fl1.tiff");

gfp_names = newArray("Lipo24-Mock-G", "Lipo24-GFP-G",
					 "Lipo48-Mock-G", "Lipo48-GFP-G",
					 "Fugene-Mock-G", "Fugene-GFP-G");


hoescht = newArray("plate_1/r03c03f03p01-ch1sk1fk1fl1.tiff",
				   "plate_1/r03c04f03p01-ch1sk1fk1fl1.tiff",
				   "plate_2/r03c05f03p01-ch1sk1fk1fl1.tiff",
				   "plate_2/r03c06f03p01-ch1sk1fk1fl1.tiff",
				   "plate_2/r03c03f03p01-ch1sk1fk1fl1.tiff",
				   "plate_2/r03c04f03p01-ch1sk1fk1fl1.tiff");

hoescht_names = newArray("Lipo24-Mock-H", "Lipo24-GFP-H",
						 "Lipo48-Mock-H", "Lipo48-GFP-H",
						 "Fugene-Mock-H", "Fugene-GFP-H");


// Make stack of GFP
for (i=0; i<gfp.length; i++) {
	open(input_path + input_folder + gfp[i]);
	rename(gfp_names[i]);
}

run("Images to Stack", "name=[gfp] use");
setMinAndMax(200, 6000);
run("Apply LUT", "stack");
run("Scale Bar...", "width=1000 thickness=7 hide overlay");

for (i=1; i<=gfp.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+gfp_names[i-1]);
	close();
}


// Make stack of Hoescht
for (i=0; i<hoescht.length; i++) {
	open(input_path + input_folder + hoescht[i]);
	rename(hoescht_names[i]);
}

run("Images to Stack", "name=[hoescht] use");
setMinAndMax(1250, 12500);
run("Apply LUT", "stack");
run("Scale Bar...", "width=1000 thickness=7 hide overlay");

for (i=1; i<=hoescht.length; i++) {
	setSlice(i);
	run("Duplicate...", "title=temp");
	saveAs(output_format, output_path+output_folder+hoescht_names[i-1]);
	close();
}

close("*");
