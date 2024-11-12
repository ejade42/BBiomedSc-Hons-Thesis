//path=getDirectory("Choose a data folder"); print(path); exit;
input_path = "U:/Evelyn Jade/Experiments/241008 - Exp 17 - First NIID IHC/"

output_path = "U:/Evelyn Jade/Experiments/241007 - Exp 16 - Microscopy images for thesis/7 - NIID IHC/"
output_folder = "imagej_export/"
output_format = "TIFF"


// Slide 2 region 1
// makeRectangle(7660, 6950, 420, 315)

// Slide 2 region 2
// makeRectangle(6335, 890, 420, 315)

// Levels: Hoechst 0/250, 488 120/750, 647 120/2000.


// Tick crop on import, set width to 15000 otherwise won't load.
open(input_path + "241017_NOTCH_20x_Slide2.nd2");

setSlice(1);
setMinAndMax(0, 750);
run("Apply LUT");   // click "no" on apply to all
run("Duplicate...", "title=temp");
makeRectangle(7660, 6950, 420, 315);
run("Crop");
saveAs(output_format, output_path+output_folder+"slide_2_region_1_hoechst");
close();

run("Duplicate...", "title=temp");
makeRectangle(6335, 890, 420, 315);
run("Crop");
saveAs(output_format, output_path+output_folder+"slide_2_region_2_hoechst");
close();


setSlice(2);
setMinAndMax(120, 750);
run("Apply LUT");
run("Duplicate...", "title=temp");
makeRectangle(7660, 6950, 420, 315);
run("Crop");
saveAs(output_format, output_path+output_folder+"slide_2_region_1_notch");
close();

run("Duplicate...", "title=temp");
makeRectangle(6335, 890, 420, 315);
run("Crop");
saveAs(output_format, output_path+output_folder+"slide_2_region_2_notch");
close();



setSlice(3);
setMinAndMax(120, 2000);
run("Apply LUT");
run("Duplicate...", "title=temp");
makeRectangle(7660, 6950, 420, 315);
run("Crop");
saveAs(output_format, output_path+output_folder+"slide_2_region_1_p62");
close();

run("Duplicate...", "title=temp");
makeRectangle(6335, 890, 420, 315);
run("Crop");
saveAs(output_format, output_path+output_folder+"slide_2_region_2_p62");
close();
close();
