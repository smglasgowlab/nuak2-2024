// Script for IF Analysis

// Author: Aneesh Dalvi
// Glasgow Lab, University of California, San Diego
// Date: July 2, 2024

// Performs batch image analysis of cell counting of cells fluoresced with DAPI, TXRed, and GFP
// Returns ImageJ overlay of cell counts, composite_images folder containing overlapped TXRed and GFP channels
// image_annotations folder containing ROIs of individual channels

// Directory Structure:
// Main/
// ├── Image folders/
// │   ├── stack1/
// │   │   ├── image.ets
// ├── corresponding_file.vsi

// Output Directory Structure:
// Main/
// ├── Image folders/
// │   ├── stack1/
// │   │   ├── image.ets
// ├── composite_images/
// ├── image_annotations/
// ├── corresponding_file.vsi


#@ File(label="Input directory", style="directory") path
path = path + "/"
path = correct_path(path);


sample_files = getFileList(path);
sample_files = Array.sort(sample_files);

setBatchMode(true);

for (i = 0; i < sample_files.length/2; i++) {
	
	if (sample_files[i] != "annotated_images") {
	
		// Access nested stack file
		stack_folder = getFileList(path + sample_files[i]);
		stack_file = getFileList(path + sample_files[i] + stack_folder[0]);
		stack_file_path = path + sample_files[i] + stack_folder[0] + stack_file[0];
		
		// Bioformat open settings
		run("Bio-Formats Importer", 
		"open=[" + stack_file_path + "] autoscale color_mode=Colorized rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT series_1");
		
		
		// Perform single-positive analysis for each frame
		for (j = 0; j < 3; j++) {
			single_count_IF(j, sample_files[i]);
		}
		
		// Perform double-positive analysis for each frame
		double_count_IF(sample_files[i]);
			
		move_annotations();
		move_composites();
	}
}


// Image analysis macro. Using IsoData dark threshold, watershed, then analyzing particles with size=20-25000 and circularity 0.40-1.00.
function single_count_IF(count, sample_name) { 
	if (count == 0) {
		// Image analysis for DAPI
		selectImage(stack_file_path + " - DAPI, GFP, TXRed - C=" + count);
		run("Duplicate...", " ");
		selectImage(stack_file_path + " - DAPI, GFP, TXRed - C=" + count + "-1");
		run("8-bit");
		run("Auto Local Threshold", "method=Phansalkar radius=15 parameter_1=0 parameter_2=0 white");
		run("Fill Holes");
		run("Analyze Particles...", "size=150-10000 pixel circularity=0.20-1.00 show=[Overlay Masks] summarize");		
	} else {
		// Image analysis
		selectImage(stack_file_path + " - DAPI, GFP, TXRed - C=" + count);
		run("Duplicate...", " ");
		selectImage(stack_file_path + " - DAPI, GFP, TXRed - C=" + count + "-1");
		run("16-bit");
		setAutoThreshold("IsoData dark");
		run("Threshold...");
		run("Convert to Mask");
		run("Watershed");
		run("Analyze Particles...", "size=150-25000 pixel circularity=0.20-1.00 show=[Overlay Masks] exclude summarize");
	}
	
	
	// Save annotated image
	selectImage(stack_file_path + " - DAPI, GFP, TXRed - C=" + count + "-1");
	sub_sample_name = substring(sample_name, 0, sample_name.length - 2);
	
	if (count == 0) {
		saveAs("PNG", path + "annotated_" + sub_sample_name + "_DAPI");
	} else if (count == 1) {
		saveAs("PNG", path + "annotated_" + sub_sample_name + "_GFP");
	} else if (count == 2) {
		saveAs("PNG", path + "annotated_" + sub_sample_name + "_TXRed");
	}
}


function double_count_IF(sample_name) {
// Create composite image
	run("Merge Channels...", "c1=[" + stack_file_path + " - DAPI, GFP, TXRed - C=2] c2=[" + stack_file_path + " - DAPI, GFP, TXRed - C=1] create keep");
	
	// Image analysis
	selectImage(stack_file_path + " - DAPI, GFP, TXRed - C=1");
	run("16-bit");
	setAutoThreshold("IsoData dark");
	run("Convert to Mask");
	
	selectImage(stack_file_path + " - DAPI, GFP, TXRed - C=2");
	run("16-bit");
	setAutoThreshold("IsoData dark");
	run("Convert to Mask");
	
	imageCalculator("AND create", stack_file_path + " - DAPI, GFP, TXRed - C=1", stack_file_path + " - DAPI, GFP, TXRed - C=2");
	selectImage("Result of " + stack_file_path + " - DAPI, GFP, TXRed - C=1");
	run("Watershed");

	run("Analyze Particles...", "size=150-25000 pixel circularity=0.20-1.00 show=[Overlay Masks] exclude summarize");
	
	
	// Save composite image
	sub_sample_name = substring(sample_name, 0, sample_name.length - 2);
	selectImage("Composite");
	saveAs("PNG", path + "composite" + sub_sample_name);
	
	// Save annotated image
	selectImage("Result of " + stack_file_path + " - DAPI, GFP, TXRed - C=1");
	saveAs("PNG", path + "annotated_DP" + sub_sample_name);
}


// Moves all annotated images to a separate "image_annotations" folder
function move_annotations() {
	File.makeDirectory(path + "image_annotations");
	files = getFileList(path);
	for (i = 0; i < files.length; i++) {
		if (substring(files[i], 0, 3) == "ann") {
			path1 = path + files[i];
			path2 = path + "image_annotations/" + files[i];
			File.rename(path1, path2);
		}
	}
}


// Moves all annotated images to a separate "image_annotations" folder
function move_composites() {
	File.makeDirectory(path + "composite_images");
	files = getFileList(path);
	for (i = 0; i < files.length; i++) {
		if (substring(files[i], 0, 3) == "com") {
			path1 = path + files[i];
			path2 = path + "composite_images/" + files[i];
			File.rename(path1, path2);
		}
	}
}


// Fix output from path parameter prompt, replacing "\" with "/"
function correct_path(incorrect_path) {
	return replace(incorrect_path, "\\", "/");
}