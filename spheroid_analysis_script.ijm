// Spheroid Area Measurement of a Single Plate

// Author: Aneesh Dalvi
// Glasgow Lab, University of California, San Diego
// Date: October 1, 2024

// Spheroid area measurement of U87, U251, LN229, and LN319 cell lines. Returns annotated images with 
// spheroid perimeter outlined in blue.

// Directory Structure:
// Main/
// ├── Plate folder/
// │   ├── adjusted_images/
// │   ├── raw_images/

// Output Directory Structure:
// Main/
// ├── Plate folder/
// │   ├── adjusted_images/
// │   ├── raw_images/
// │   ├── ROI_images/


#@ File(label="Input directory", style="directory") path
path = path + "/"
path = correct_path(path);

main();

function main() {
	// Collect list of all images in plate 
	image_folders = getFileList(path);
		
	// Reset ROI images folder 
	roi_images_path = path + "ROI_images/";
	roi_files = getFileList(roi_images_path);
	if (roi_files.length > 0) {
		for (j = 0; j < roi_files.length; j++) {
			File.delete(roi_images_path + roi_files[j]);	
		}
	}
	// Initialize settings
	settings();
	
	// Iterate through images 
	images = getFileList(path + "raw_images/");
	for (k = 0; k < images.length; k++) {
		image_path = path + "raw_images/" + images[k];
		if (endsWith(image_path, "TIF")) {
			measure_area(image_path, roi_images_path, images[k]);
		}
	}
	
	print("Spheroid analysis complete.");
}

function settings() {
	setBatchMode(true);
	setOption("BlackBackground", true);
	run("Set Measurements...", "area display redirect=None decimal=0");
}

// Macro to measure spheroid area
function measure_area(image_path, roi_images_path, sample_name) {
	open(image_path);	
	
	Image.removeScale();
	roiManager("Reset");
	
	image = getImageID();
	Overlay.remove;

	// Filtering
	run("Duplicate...", " ");
	run("Gamma...", "value=0.75");
	//run("Subtract Background...", "rolling=25 light sliding");
	//run("Multiply...", "value=1.4");
	run("Median...", "radius=5");
	run("Enhance Contrast...", "saturated=0.35");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Fill Holes (Binary/Gray)");
	run("EDM Binary Operations", "iterations=6 operation=erode");
	run("EDM Binary Operations", "iterations=6 operation=dilate");
	run("Morphological Filters", "operation=Closing element=Disk radius=12");
	run("Fill Holes (Binary/Gray)");
	run("Adjustable Watershed", "tolerance=2");
	image2 = getImageID();
	Image.removeScale();
	
	run("Analyze Particles...", "size=2000-2000000 circularity=0.20-1.00 show=Masks include");
	image3 = getImageID();
	Image.removeScale();
	run("Erode");
	run("Dilate");
	
	imageCalculator("OR create", image2, image3);
	Image.removeScale();
	
	// Connected component analysis
	run("Analyze Particles...", "size=10000-2000000 circularity=0.20-1.00 include add composite");
	close("Mask");
	
	// Measure on the original image
	selectImage(image);
	Image.removeScale();
	roiManager("Set Color", "cyan");
	roiManager("Set Line Width", 6);
	roiManager("Show all");
	run("Labels...", "color=yellow font=6 show bold");
	run("Flatten");
	roiManager("Measure");
	
	// Save Image
	saveAs("png", roi_images_path + sample_name);
}


// Fix output from path parameter prompt, replacing "\" with "/"
function correct_path(incorrect_path) {
	return replace(incorrect_path, "\\", "/");
}