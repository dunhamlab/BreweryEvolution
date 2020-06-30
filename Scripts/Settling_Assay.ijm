Directory = "/Users/chris/Google Drive/Dunham Lab/Postdoc_Brewing/SettlingAssay/"

list = getFileList(Directory + "20180621_Images/");
setBatchMode(false);
//print(list);
for(i = 0; i < list.length - 1; i++)
	OpenAndProcess(list[i]);

function OpenAndProcess(FileName) {
	//print(Directory + "20180621_Images/" + list[i]);
	open(Directory + "20180621_Images/" + FileName);
	title_JPG=getTitle();
	run("Channels Tool...");
	run("Make Composite");
	Stack.setDisplayMode("grayscale");
	setTool("rectangle");
	waitForUser("Pause","Select rectangle for cropping"); 
	run("Crop");
	
	saveAs("Tiff", Directory + "20180621_Images/" + FileName + ".tiff");
	title_Tiff=getTitle();
	for(j=0; j < 3; j++)
		PlotProfile(FileName);
	selectWindow(title_Tiff);
	close();
}

function PlotProfile(FileName) {
	print(j);
	setTool("line");
	waitForUser("Pause","Select line for plot profile"); 
	run("Plot Profile");
	Plot.getValues(x, y);
	openFile = File.open(Directory + "MinimumValues_20180621/" + FileName + "" + j + ".txt");
	for (k=0; k<x.length; k++)
		print(openFile, x[k] + "\t" + y[k]);
			
	File.close(openFile);
	selectWindow("Plot of " + title_JPG);
	close();
}
