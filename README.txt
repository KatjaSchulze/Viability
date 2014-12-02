This plugin determines the vitality of Synechocystis cells via red chlorophyll fluorescence (indicates live cells) and an unspecific green fluorescence (indicates dead cells)

For the use ImageJ is needed and can be dowloaded from http://rsbweb.nih.gov/ij/. Copy the jar-file in the ImageJ plugin folder and restart ImageJ. 

Usage:
1 ) Start the plugin under Plugin>LivingDead
2 ) In a first dialog the parameters for the concentration calculation can be adjusted (height and width of the image, depth of the used counting chamber). 
    Additionally it can be specified if images should be shown during the anlysis and if the result should be stored.
3 ) In a second dialog you have to specify the folder which contains the images that should be analyzed (images must be tif, jpeg or bmp).
4 ) The automated analysis will start for all images in the specified folder 
5 ) Results are displayed as a normalized histogram for the mean green intensity. The log shows the total cell count of all, red(live) and green(dead) cells, the percentage of red(live) and green(dead) cells in the sample and the calculated total concentration in cells/ml.  
