import java.io.BufferedReader;
import java.io.FileReader;

import com.google.gson.Gson

import qupath.lib.objects.PathAnnotationObject;
import qupath.lib.roi.ROIs;
import qupath.lib.roi.EllipseROI;
import qupath.lib.measurements.Measurement;





//#############################################
//#  Following Parameters need to be changed  #
//#############################################
//1a. Create a folder where your scalefactors_json.json and tissue_positions_list.csv are stored. DO NOT ADD "/" at the end
//1b. For windows users, change all the '/' symbol to '\'
//2. Meassure the actual fiducial diameter on the image. You can use Line Functions in QuPath 

visium_spot_drawer('/path/to/sample/folder/',100)


//##############################################################################
def visium_spot_drawer(String file_directory, double fiducial_diameter_qupath){
//Call two essential files
def json_file = new File(file_directory + '/scalefactors_json.json').text
def csvfile = file_directory + '/tissue_positions_list.csv'

//Align the pixel length
def gson = new Gson()
def json = gson.fromJson(json_file, Map.class)
def fiducial_diameter_fullres = json.fiducial_diameter_fullres.toDouble()
def spot_diameter_fullres = json.spot_diameter_fullres.toDouble();
double spot_diameter_qupath = spot_diameter_fullres* fiducial_diameter_qupath/  fiducial_diameter_fullres;




// Create BufferedReader
def imageData = getCurrentImageData();
def csvReader = new BufferedReader(new FileReader(csvfile));
def plane = ImagePlane.getDefaultPlane();

// Create an empty list to store the Data
def spots = [] 

// Loop through all the rows of the CSV file.
while ((row = csvReader.readLine()) != null) {

    def rowContent = row.split(",")
    
    String barcode = Arrays.toString(rowContent[0]);
    barcode = barcode.replaceAll("[\\[\\]\"]", "");
    int  in_tissue = rowContent[1] as int;
    int  array_row = rowContent[2] as int;
    int  array_col = rowContent[3] as int;
    double cy = rowContent[4] as double;
    double cx = rowContent[5] as double;
    
    int first_time = 1;
    // Create annotation
    if (in_tissue) {
        double px = cx - spot_diameter_qupath/2; 
        double py = cy - spot_diameter_qupath/2; 
        def roi = new EllipseROI(px, py, spot_diameter_qupath, spot_diameter_qupath, plane)
        
        def annotation = new PathAnnotationObject(roi, PathClassFactory.getPathClass("Spot"));
        annotation.setName(barcode);
        
        if (first_time) {
            annotation.getMeasurementList().addMeasurement("array_row", array_row);
            annotation.getMeasurementList().addMeasurement("array_col", array_col);
            annotation.getMeasurementList().addMeasurement("cx", cx);
            annotation.getMeasurementList().addMeasurement("cy", cy);
            first_time = 0;
        } else {
            annotation.getMeasurementList().putMeasurement("array_row", array_row);
            annotation.getMeasurementList().putMeasurement("array_col", array_col);
            annotation.getMeasurementList().putMeasurement("cx", cx);
            annotation.getMeasurementList().putMeasurement("cy", cy);
        }  
         
         spots << annotation
        
    }
    addObjects(spots)
}
}
