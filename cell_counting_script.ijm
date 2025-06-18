// Script for IF Analysis

// Author: Aneesh Dalvi
// Glasgow Lab, University of California, San Diego
// Date: July 2, 2024
// Contact: adalvi@ucsd.edu

// Folders containing the .ets file should have the same name as the corresponding .vsi files with a leading and trailing _
// Example: _U251_WT_1_ is a directory containing the nested .ets file 
// 			U251_WT_1.vsi is the corresponding vsi file
// No other files should be present in the main directory

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
setBatchMode(true);

main();

function main() {
	sample_files = getFileList(path);
	sample_files = Array.sort(sample_files);
	
	// Prompt map channel numbers 
	showMessage("Three Channel IHC Cell Counter", 
    "For the following prompts, enter -1 if not applicable.\n" +
    "Channel numbers for images can be viewed by opening the images individually.\n" +
    "The number in the part of the image title 'C=_' corresponds to the channel number.");
	num_ch = getNumber("How many channels are there? (Up to 3)", 3);
	DAPI_ch = getNumber("Enter the DAPI channel:", -1);
	GFP_ch = getNumber("Enter the GFP channel:", -1);
	TXRed_ch = getNumber("Enter the TXRed channel:", -1);
	if (num_ch > 1) {
		count_db = getBoolean("Would you like to perform double-positive cell counting?");
		if (count_db) {
			dp_ch_1 = getNumber("Select the first channel for double-positive counting (0, 1, or 2):", -1);
			dp_ch_2 = getNumber("Select the second channel for double-positive counting (0, 1, or 2):", -1);
		}
	}

	for (i = 0; i < sample_files.length/2; i++) {
		// Access nested stack file
		stack_folder = getFileList(path + sample_files[i]);
		stack_file = getFileList(path + sample_files[i] + stack_folder[0]);
		stack_file_path = path + sample_files[i] + "stack1/" + stack_file[0];
			
		// Bioformat open settings
		run("Bio-Formats Importer", 
		"open=[" + stack_file_path + "] autoscale color_mode=Colorized rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT series_1");
		image_names = getList("image.titles");
			
		// Perform single-positive analysis for each frame
		channel_names = newArray(num_ch);  // Initialize an empty array
		for (j = 0; j < num_ch; j++) {
		    if (j == DAPI_ch) {
		        DAPI = single_count_IF(j, sample_files[i]);
		        channel_names[DAPI_ch] = DAPI;
		    } else if (j == GFP_ch) {
		        GFP = single_count_IF(j, sample_files[i]);
		        channel_names[GFP_ch] = GFP;
		    } else if (j == TXRed_ch) {
		        TXRed = single_count_IF(j, sample_files[i]);	
		        channel_names[TXRed_ch] = TXRed;
		    }
		}
		
		roiManager("Reset");
		
		// Perform double-positive analysis for each frame
		if (num_ch > 1 && count_db == true) {
			dp_name_1 = channel_names[dp_ch_1];
			dp_name_2 = channel_names[dp_ch_2];
			double_count_IF(sample_files[i], dp_name_1, dp_name_2);
		}
				
		// Move annotations and close all open windows
		move_annotations();
		close("*");
	}
	
	print("Cell counting complete.");
}
	
// Image analysis macro. Using IsoData dark threshold, watershed, then analyzing particles with size=150-25000 and circularity 0.20-1.00.
function single_count_IF(count, sample_name) { 
	roiManager("Reset");
	image_name = image_names[count];
	selectImage(image_name);
	image1 = getTitle();
	Overlay.remove;
		
	run("Duplicate...", " ");
	image2 = getTitle();
	run("Multiply...", "value=1.7");
	
	name = getTitle();
	run("Morphological Filters", "operation=[Internal Gradient] element=Disk radius=2");

	run("8-bit");
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("EDM Binary Operations", "iterations=5 operation=close");
	run("Fill Holes (Binary/Gray)");
	run("Adjustable Watershed", "tolerance=0.5");
	image3 = getTitle();
	
	// Save the processed images for GFP, TXRed, rename
	// Call doublecount and use imagecalculator to save the overlaps between images
	if (count == 1) {
		run("Duplicate...", " ");
		name2 = getTitle();
		name2 = name2 + "-1";
		rename(name2);
	} else if (count == 0) {
		run("Duplicate...", " ");
		name2 = getTitle();
		name2 = name2 + "-0";
		rename(name2);
	} else {
		name2 = "None";
	}
	
	selectImage(image3);
	rename(name);
	run("Analyze Particles...", "size=50-Infinity pixel circularity=0.20-1.00 display exclude include summarize add");
	selectImage(image1);
	roiManager("Set Color", "cyan");
	roiManager("Set Line Width", 2);
	roiManager("Show all");
	run("Labels...", "color=yellow font=6 show bold");
	run("Flatten");
	
	// Save annotated image
	sub_sample_name = substring(sample_name, 0, sample_name.length - 2);
	channels = newArray("_DAPI", "_TXRed", "_TXRed");
	saveAs("PNG", path + "annotated_" + sub_sample_name + channels[count]);
	
	return name2;
}


function double_count_IF(sample_name, dp_name_1, dp_name_2) {
	// Create composite image
	run("Merge Channels...", "c1=[" + image_names[dp_ch_1] + "] c2=[" + image_names[dp_ch_2] + "] create keep");
	image4 = getTitle();
	
	// Image Analysis
	roiManager("Reset");
	Overlay.remove;
	imageCalculator("AND create", dp_name_1, dp_name_2);
	run("Analyze Particles...", "size=50-Infinity pixel circularity=0.20-1.00 display exclude include summarize add");
	
	selectImage(image4);
	roiManager("Set Color", "yellow");
	roiManager("Set Line Width", 2);
	roiManager("Show all");
	run("Labels...", "color=yellow font=6 show bold");
	run("Flatten");
	
	// Save annotated image
	sub_sample_name = substring(sample_name, 0, sample_name.length - 2);
	saveAs("PNG", path + "annotated_DP" + sub_sample_name);
}


// Moves all annotated images to a separate "image_annotations" folder
function move_annotations() {
	if (File.exists(path + "image annotations")) {
		File.delete(path + "image_annotations");
	}
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


// Fix output from path parameter prompt, replacing "\" with "/"
function correct_path(incorrect_path) {
	return replace(incorrect_path, "\\", "/");
}