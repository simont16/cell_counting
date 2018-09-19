from ij import IJ
from glob import glob
import os, sys, re

def day_finder(folder_loc, experiment, virus, days):
    day_find1 = glob (folder_loc + '/' + experiment + '/' + virus +'/')
    day_find2 = glob(day_find1[0]+'/*/')    
    for ii in range(len(day_find2)):
            first_split = re.split(r'/', day_find2[ii])
            second_split = re.split(r'_', first_split[9])
           #print(second_split)
            third_split = second_split[1].replace(month, '')
            print(third_split)
            days.append(third_split)

def image_finder(file_paths,image_name):
    for iii in range(len(file_paths)):
        first_file_split = re.split(r'\\', file_paths[iii])
        second_file_split = first_file_split[1]
        image_name.append(second_file_split)
    return image_name

folder_loc = 'C:/Users/SWC/Dropbox (UCL - SWC)/Data Files/Simon/Data'
experiment = 'HSV Rabies/HSV_HEK_Transfection_18_Aug2018'
virus = 'N2c'
file_type = '_I_'
month = '2018'
experiment_type = 'HSVRepeat_EVOS10x'
days = ['2908', '3108' , '0209']

extension = ['*RFP.tif', '*GFP.tif']
threshold = ['Li', 'Triangle']
plate = ['P1', 'P2']

#day_finder(folder_loc, experiment, virus, days)

for g in range(len(plate)):
	for b in range(len(threshold)):
		for f in range(len(extension)):
			for p in range(len(days)):
			    ## select folder of interest
			    image_folder = folder_loc  +'/'+ experiment +'/'+ virus +'/'+ 'ST_' + month + days[p] + '_' + plate[g] + file_type + virus + '_' + experiment_type + "/"
			    ## folder to which results will be saved
			    print(image_folder)
			    result_fol_step_1 = folder_loc +'/' + experiment +'/' + 'results' + '/'
			    result_fol_step_2 = result_fol_step_1  + month + days[p] + file_type + virus + '_' + experiment_type +'/' 
			    result_fol_step_3 = result_fol_step_2 +  threshold[b] + "/"
			    #print('working folder: {0} \nsaving to: {1}'.\
			          #format(image_folder,result_fol_step_3 ))
			     #check if results folder exisits, else create it
			    if not os.path.isdir(result_fol_step_1):
			        #print('results folder not found \ncreating...')
			        os.mkdir(result_fol_step_1)
			    if not os.path.isdir(result_fol_step_2):
			        #print('results folder not found \ncreating...')
			        os.mkdir(result_fol_step_2)
			    if not os.path.isdir(result_fol_step_3):
			        #print('results folder not found \ncreating...')
			        os.mkdir(result_fol_step_3)
			    
			    result_folder = result_fol_step_3
			    
			     ## file extension
			    # extension = '*.jpeg'
			    file_paths = glob(image_folder +'/'+ extension[f])
			    #print(file_paths)
			    file_num = len(file_paths)
			    #print file_num
			    image_name = []
			    image_finder(file_paths, image_name)
			
			    for q in range(len(image_name)):
			    	#print file_paths[q]
			    	imp = IJ.openImage(file_paths[q])
			    	imp.show()
			    	IJ.run("Set Measurements...", "mean redirect=None decimal=3");
			    	
			    	IJ.run(imp, "Measure", "");
			    	
			    	IJ.run("32-bit");
			    	IJ.setAutoThreshold(imp, threshold[b]);
			    	IJ.run("Convert to Mask");
			    	IJ.run("Invert");
			    	IJ.run("Erode");
			    	IJ.run("Erode");
			    	IJ.run("Watershed");
			    	IJ.run("Analyze Particles...", "size=40 -Infinity show=Ellipses summarize");
			    	IJ.selectWindow( 'ST_' + month + days[p]+ '_' + plate[g] + file_type + virus + '_' + experiment_type + "\\" + image_name[q]);
			    	IJ.run("Close");
			    	IJ.selectWindow("Drawing of " +  'ST_' + month + days[p] + '_' + plate[g] + file_type + virus + '_' + experiment_type + "\\" + image_name[q]);
			    	IJ.run("Invert");
			
			    	imp2 = IJ.openImage(file_paths[q])
			    	imp2.show()
			    	IJ.run(imp2, "Enhance Contrast...", "saturated=0")
			    	IJ.selectWindow( 'ST_' + month + days[p] + '_' + plate[g] + file_type + virus + '_' + experiment_type + "\\" + image_name[q])
			    	IJ.run("Add Image...", "image=[Drawing of " +  'ST_' + month + days[p] + '_' + plate[g] + file_type + virus + '_' + experiment_type + "\\" + image_name[q] + "] x=0 y=0 opacity=90 zero");
			    	IJ.selectWindow( 'ST_' + month + days[p] + '_' + plate[g] + file_type + virus + '_' + experiment_type + "\\" + image_name[q])	
			    	imp3 = imp2.flatten();
			    	IJ.saveAs(imp3, "Tiff", result_folder + days[p] + image_name[q])
			    	
			    	IJ.run("Close All", "");
			    		    		    	
IJ.selectWindow("Summary");	
IJ.saveAs("Results", result_folder + "cell_count" + '_' + extension[f].replace("*", "") + ".csv");
#IJ.deleteRows(0,527);
IJ.run("Close");
IJ.selectWindow("Results");	
IJ.saveAs("Results", result_folder + "mean_grey" + '_' + extension[f].replace("*", "") + ".csv");
#IJ.deleteRows(0,527);
IJ.run("Close");
		   
			
	    