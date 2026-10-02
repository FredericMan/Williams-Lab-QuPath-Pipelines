/**
 * ============================================================
 * TVA/Rabies two-pass cell detection - v9
 * ============================================================
 *
 * QuPath 0.7.0
 *
 * MULTI-ROI VERSION
 *
 * Select ONE annotation/ROI before running the script.
 *
 *
 * PASS 1
 * -------
 * FITC/GFP nuclear detection
 *       ↓
 * cell expansion
 *       ↓
 * measure TRITC/mCherry
 *       ↓
 * GFP+ / mCherry+  = Starter
 * GFP+ / mCherry-  = GFP-only
 *
 *
 * PASS 2
 * -------
 * TRITC/mCherry detection
 *       ↓
 * compare mCherry detection centroid with saved GFP Starter
 * detections
 *       ↓
 * matched GFP Starter = Starter
 * no GFP Starter match = RV+
 *
 *
 * FINAL
 * -----
 * Cell (Starter)
 * Cell (GFP-only)
 * Cell (RV+)
 *
 *
 * IMPORTANT
 * ----------
 * Only detections belonging to the CURRENT selected annotation
 * are modified or removed.
 *
 * Detections belonging to previously analyzed ROIs are preserved.
 *
 * Pass-1 GFP information is saved before Pass 2.
 *
 * After Pass 2, ONLY GFP-only cells are recreated.
 *
 * GFP Starter objects are NOT recreated because the corresponding
 * Pass-2 mCherry detections are the final Starter objects.
 *
 * ============================================================
 */


// ============================================================
// IMPORTS
// ============================================================

import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import static qupath.lib.gui.scripting.QPEx.*


// ============================================================
// USER PARAMETERS
// ============================================================


// ------------------------------------------------------------
// CHANNEL NAMES
// ------------------------------------------------------------

def GFP_CHANNEL = "FITC"
def MCHERRY_CHANNEL = "TRITC"


// ------------------------------------------------------------
// CURRENT ROI / ANNOTATION
// ------------------------------------------------------------

def analysisParent = getSelectedObject()

if (analysisParent == null || !analysisParent.isAnnotation()) {

    throw new Exception(
        "Please select exactly ONE annotation/ROI before running the script."
    )
}

def ROI_NAME = analysisParent.getName()

if (ROI_NAME == null || ROI_NAME.trim().isEmpty()) {
    ROI_NAME = "Unnamed_ROI"
}


// ============================================================
// PASS 1 - GFP DETECTION PARAMETERS
// ============================================================

def GFP_PIXEL_SIZE = 1.0

def GFP_BACKGROUND_RADIUS = 10.0
def GFP_MEDIAN_RADIUS = 2.0
def GFP_SIGMA = 5.0

def GFP_MIN_AREA = 6.0
def GFP_MAX_AREA = 600.0

def GFP_THRESHOLD = 200.0

def GFP_CELL_EXPANSION = 3.0


// ============================================================
// PASS 1 - mCHERRY STARTER THRESHOLD
// ============================================================

def MCHERRY_STARTER_THRESHOLD = 2000.0


// ============================================================
// PASS 2 - mCHERRY DETECTION PARAMETERS
// ============================================================

def MCHERRY_PIXEL_SIZE = 1.0

def MCHERRY_BACKGROUND_RADIUS = 10.0
def MCHERRY_MEDIAN_RADIUS = 2.0
def MCHERRY_SIGMA = 5.0

def MCHERRY_MIN_AREA = 10.0
def MCHERRY_MAX_AREA = 600.0

def MCHERRY_THRESHOLD = 600.0

def MCHERRY_CELL_EXPANSION = 3.0


// ============================================================
// STARTER SPATIAL MATCHING
// ============================================================

def STARTER_MATCH_DISTANCE_MICRONS = 5.0


// ============================================================
// EXPORT
// ============================================================

def EXPORT_RESULTS = true

def EXPORT_FILENAME = "TVA_Rabies_Cell_Analysis_v9"


// ============================================================
// INITIAL CHECK
// ============================================================

def imageData = getCurrentImageData()

if (imageData == null) {
    throw new Exception("No image is currently open.")
}


println ""
println "================================================"
println " TVA/Rabies two-pass cell analysis - v9"
println " QuPath 0.7.0"
println "================================================"
println ""

println "Analysis ROI:    " + ROI_NAME
println "GFP channel:     " + GFP_CHANNEL
println "mCherry channel: " + MCHERRY_CHANNEL
println ""


// ============================================================
// FUNCTION:
// GET ONLY DETECTIONS DIRECTLY UNDER CURRENT ROI
// ============================================================
//
// This is deliberately NOT getDetectionObjects().
//
// getDetectionObjects() returns detections globally and therefore
// is unsuitable when several analyzed ROIs coexist in the image.
//

def getCurrentROIDetections = {

    return analysisParent.getChildObjects().findAll {
        it.isDetection()
    }
}


// ============================================================
// REMOVE EXISTING DETECTIONS FROM CURRENT ROI ONLY
// ============================================================

def existingDetections = getCurrentROIDetections()

if (!existingDetections.isEmpty()) {

    println(
        "Removing existing detections from current ROI: " +
        existingDetections.size()
    )

    analysisParent.removeChildObjects(
        existingDetections
    )

    fireHierarchyUpdate()

} else {

    println "No existing detections found in current ROI."
}

println ""

println "Detections in other ROIs are preserved."
println ""


// ============================================================
// PASS 1
// GFP DETECTION
// ============================================================

println "----------------------------------------------"
println "PASS 1: GFP detection"
println "ROI: " + ROI_NAME
println "----------------------------------------------"


runPlugin(
    'qupath.imagej.detect.cells.PositiveCellDetection',

    """{
        "detectionImage":"${GFP_CHANNEL}",
        "requestedPixelSizeMicrons":${GFP_PIXEL_SIZE},
        "backgroundRadiusMicrons":${GFP_BACKGROUND_RADIUS},
        "backgroundByReconstruction":true,
        "medianRadiusMicrons":${GFP_MEDIAN_RADIUS},
        "sigmaMicrons":${GFP_SIGMA},
        "minAreaMicrons":${GFP_MIN_AREA},
        "maxAreaMicrons":${GFP_MAX_AREA},
        "threshold":${GFP_THRESHOLD},
        "watershedPostProcess":true,
        "cellExpansionMicrons":${GFP_CELL_EXPANSION},
        "includeNuclei":true,
        "smoothBoundaries":true,
        "makeMeasurements":true,
        "thresholdCompartment":"Cell: ${MCHERRY_CHANNEL} mean",
        "thresholdPositive1":${MCHERRY_STARTER_THRESHOLD},
        "thresholdPositive2":${MCHERRY_STARTER_THRESHOLD},
        "thresholdPositive3":${MCHERRY_STARTER_THRESHOLD},
        "singleThreshold":true
    }"""
)


println "GFP detection completed."
println ""


// ============================================================
// GET PASS-1 DETECTIONS FROM CURRENT ROI ONLY
// ============================================================

def gfpDetections = getCurrentROIDetections()

println(
    "Number of GFP-derived detections in " +
    ROI_NAME + ": " +
    gfpDetections.size()
)

println ""


// ============================================================
// DATA STRUCTURES FOR PASS 1
// ============================================================

def gfpStarterData = []

def gfpOnlyData = []


// ============================================================
// DIAGNOSTIC OUTPUT + SAVE PASS-1 DATA
// ============================================================

println "GFP detection measurements:"
println "----------------------------------------------"


gfpDetections.eachWithIndex { detection, index ->

    def ml = detection.getMeasurementList()


    // --------------------------------------------------------
    // Measurements
    // --------------------------------------------------------

    def gfpNucleus =
        ml.get("Nucleus: ${GFP_CHANNEL} mean")

    def gfpCell =
        ml.get("Cell: ${GFP_CHANNEL} mean")

    def mcherryNucleus =
        ml.get("Nucleus: ${MCHERRY_CHANNEL} mean")

    def mcherryCell =
        ml.get("Cell: ${MCHERRY_CHANNEL} mean")


    // --------------------------------------------------------
    // Cell ROI
    // --------------------------------------------------------

    def cellROI =
        detection.getROI()


    def centroidX =
        cellROI.getCentroidX()

    def centroidY =
        cellROI.getCentroidY()


    // --------------------------------------------------------
    // Nucleus ROI
    // --------------------------------------------------------

    def nucleusROI = null

    try {

        if (detection.metaClass.respondsTo(
                detection,
                "getNucleusROI")) {

            nucleusROI =
                detection.getNucleusROI()
        }

    } catch (Exception e) {

        nucleusROI = null
    }


    // --------------------------------------------------------
    // Diagnostic output
    // --------------------------------------------------------

    println(
        "${index + 1}: " +
        "GFP nuclear mean = ${gfpNucleus}, " +
        "mCherry cell mean = ${mcherryCell}"
    )


    // --------------------------------------------------------
    // Biological classification
    // --------------------------------------------------------

    def classification

    if (
        mcherryCell != null &&
        !Double.isNaN(mcherryCell) &&
        mcherryCell >= MCHERRY_STARTER_THRESHOLD
    ) {

        classification = "Starter"

    } else {

        classification = "GFP-only"
    }


    // --------------------------------------------------------
    // Save all information required to reconstruct GFP-only
    // --------------------------------------------------------

    def record = [

        cellROI: cellROI,

        nucleusROI: nucleusROI,

        centroidX: centroidX,
        centroidY: centroidY,

        gfpNucleus: gfpNucleus,
        gfpCell: gfpCell,

        mcherryNucleus: mcherryNucleus,
        mcherryCell: mcherryCell,

        cellArea: cellROI.getArea(),

        classification: classification,

        roiName: ROI_NAME
    ]


    if (classification == "Starter") {

        gfpStarterData << record

    } else {

        gfpOnlyData << record
    }


    // --------------------------------------------------------
    // Temporarily classify the Pass-1 object
    // --------------------------------------------------------

    detection.setClassification(
        classification
    )
}


println ""

println(
    "Starter cells from GFP pass: " +
    gfpStarterData.size()
)

println(
    "GFP-only cells:              " +
    gfpOnlyData.size()
)

println ""


// ============================================================
// SAVE PASS-1 DATA CONFIRMATION
// ============================================================

println "----------------------------------------------"
println "Saving Pass-1 data"
println "----------------------------------------------"

println(
    "GFP-only objects saved: " +
    gfpOnlyData.size()
)

println(
    "GFP Starter objects saved for spatial matching: " +
    gfpStarterData.size()
)

println ""


// ============================================================
// REMOVE PASS-1 DETECTIONS FROM CURRENT ROI ONLY
// ============================================================
//
// IMPORTANT:
//
// We do NOT use removeDetections().
//
// This removes ONLY the children belonging to analysisParent.
//
// All previously analyzed ROIs remain untouched.
//

def pass1CurrentROIDetections =
    getCurrentROIDetections()

if (!pass1CurrentROIDetections.isEmpty()) {

    analysisParent.removeChildObjects(
        pass1CurrentROIDetections
    )

    fireHierarchyUpdate()
}


println "Pass-1 detections removed from current ROI."
println "Pass-1 data retained in memory."

println "Detections in other ROIs remain untouched."
println ""


// ============================================================
// PASS 2
// mCHERRY DETECTION
// ============================================================

println "----------------------------------------------"
println "PASS 2: mCherry detection"
println "ROI: " + ROI_NAME
println "----------------------------------------------"


runPlugin(
    'qupath.imagej.detect.cells.PositiveCellDetection',

    """{
        "detectionImage":"${MCHERRY_CHANNEL}",
        "requestedPixelSizeMicrons":${MCHERRY_PIXEL_SIZE},
        "backgroundRadiusMicrons":${MCHERRY_BACKGROUND_RADIUS},
        "backgroundByReconstruction":true,
        "medianRadiusMicrons":${MCHERRY_MEDIAN_RADIUS},
        "sigmaMicrons":${MCHERRY_SIGMA},
        "minAreaMicrons":${MCHERRY_MIN_AREA},
        "maxAreaMicrons":${MCHERRY_MAX_AREA},
        "threshold":${MCHERRY_THRESHOLD},
        "watershedPostProcess":true,
        "cellExpansionMicrons":${MCHERRY_CELL_EXPANSION},
        "includeNuclei":true,
        "smoothBoundaries":true,
        "makeMeasurements":true,
        "thresholdCompartment":"Cell: ${GFP_CHANNEL} mean",
        "thresholdPositive1":1500.0,
        "thresholdPositive2":1500.0,
        "thresholdPositive3":1500.0,
        "singleThreshold":true
    }"""
)


println "mCherry detection completed."
println ""


// ============================================================
// GET PASS-2 DETECTIONS FROM CURRENT ROI ONLY
// ============================================================

def mcherryDetections =
    getCurrentROIDetections()

println(
    "Number of mCherry-derived detections in " +
    ROI_NAME + ": " +
    mcherryDetections.size()
)

println ""


// ============================================================
// PASS 2 CLASSIFICATION
// ============================================================

def rvCount = 0
def starterCount = 0


// ------------------------------------------------------------
// Pixel calibration
// ------------------------------------------------------------

def pixelCalibration =
    imageData.getServer().getPixelCalibration()


def pixelSizeMicrons =
    pixelCalibration.getAveragedPixelSizeMicrons()


if (
    pixelSizeMicrons == null ||
    Double.isNaN(pixelSizeMicrons) ||
    pixelSizeMicrons <= 0
) {

    throw new Exception(
        "Could not determine image pixel size."
    )
}


def starterMatchDistancePixels =
    STARTER_MATCH_DISTANCE_MICRONS /
    pixelSizeMicrons


println(
    "Starter matching distance: " +
    STARTER_MATCH_DISTANCE_MICRONS +
    " µm"
)

println(
    "Approximate pixel distance: " +
    starterMatchDistancePixels
)

println ""


// ============================================================
// MATCH EACH mCHERRY DETECTION TO SAVED GFP STARTERS
// ============================================================

mcherryDetections.eachWithIndex {
    detection, index ->


    def roi =
        detection.getROI()


    def mcherryX =
        roi.getCentroidX()


    def mcherryY =
        roi.getCentroidY()


    def nearestDistancePixels =
        Double.POSITIVE_INFINITY


    def matchedStarter =
        null


    // --------------------------------------------------------
    // Compare against every GFP Starter
    // --------------------------------------------------------

    gfpStarterData.each { starter ->


        def dx =
            mcherryX -
            starter.centroidX


        def dy =
            mcherryY -
            starter.centroidY


        def distancePixels =
            Math.sqrt(
                dx * dx +
                dy * dy
            )


        if (
            distancePixels <
            nearestDistancePixels
        ) {

            nearestDistancePixels =
                distancePixels

            matchedStarter =
                starter
        }
    }


    // --------------------------------------------------------
    // Classification
    // --------------------------------------------------------

    if (
        matchedStarter != null &&
        nearestDistancePixels <=
            starterMatchDistancePixels
    ) {


        detection.setClassification(
            "Starter"
        )


        starterCount++


        println(
            "mCherry detection ${index + 1}: " +
            "STARTER " +
            "(distance = " +
            (
                nearestDistancePixels *
                pixelSizeMicrons
            ) +
            " µm)"
        )


    } else {


        detection.setClassification(
            "RV+"
        )


        rvCount++


        println(
            "mCherry detection ${index + 1}: RV+"
        )
    }


    // --------------------------------------------------------
    // Add ONLY numeric analysis metadata.
    //
    // DO NOT store analysisParent.getID() here.
    // QuPath MeasurementList accepts numeric measurements,
    // whereas getID() returns a UUID.
    // --------------------------------------------------------

    def ml =
        detection.getMeasurementList()


    ml.put(
        "TVA_Rabies_Analysis",
        1.0
    )
}


println ""

println(
    "Starter cells identified in Pass 2: " +
    starterCount
)

println(
    "RV+ cells identified in Pass 2: " +
    rvCount
)

println ""


// ============================================================
// RESTORE GFP-ONLY CELLS
// ============================================================
//
// GFP-only cells are recreated as TRUE CELL OBJECTS.
//
// The QuPath API uses:
//
// createCellObject(cellROI, nucleusROI, PathClass)
//
// Therefore they should appear in the hierarchy as:
//
// Cell (GFP-only)
//
// rather than:
//
// Detection (GFP-only)
//

println "----------------------------------------------"
println "Restoring GFP-only cells"
println "----------------------------------------------"


def restoredGFPOnlyObjects = []


def GFPOnlyPathClass =
    PathClass.fromString("GFP-only")


gfpOnlyData.each { record ->


    def restoredObject


    // --------------------------------------------------------
    // Create a genuine Cell object
    // --------------------------------------------------------

    restoredObject =
        PathObjects.createCellObject(
            record.cellROI,
            record.nucleusROI,
            GFPOnlyPathClass
        )


    // --------------------------------------------------------
    // Restore measurements
    // --------------------------------------------------------

    def ml =
        restoredObject.getMeasurementList()


    if (
        record.gfpNucleus != null &&
        !Double.isNaN(record.gfpNucleus)
    ) {

        ml.put(
            "Nucleus: ${GFP_CHANNEL} mean",
            record.gfpNucleus
        )
    }


    if (
        record.gfpCell != null &&
        !Double.isNaN(record.gfpCell)
    ) {

        ml.put(
            "Cell: ${GFP_CHANNEL} mean",
            record.gfpCell
        )
    }


    if (
        record.mcherryNucleus != null &&
        !Double.isNaN(record.mcherryNucleus)
    ) {

        ml.put(
            "Nucleus: ${MCHERRY_CHANNEL} mean",
            record.mcherryNucleus
        )
    }


    if (
        record.mcherryCell != null &&
        !Double.isNaN(record.mcherryCell)
    ) {

        ml.put(
            "Cell: ${MCHERRY_CHANNEL} mean",
            record.mcherryCell
        )
    }


    // --------------------------------------------------------
    // Numeric metadata
    // --------------------------------------------------------

    ml.put(
        "TVA_Rabies_Pass1_GFPonly",
        1.0
    )


    restoredGFPOnlyObjects <<
        restoredObject
}


// ============================================================
// ADD GFP-ONLY CELLS UNDER ORIGINAL ROI
// ============================================================

if (!restoredGFPOnlyObjects.isEmpty()) {

    analysisParent.addChildObjects(
        restoredGFPOnlyObjects
    )

    fireHierarchyUpdate()
}


println(
    "GFP-only cells restored: " +
    restoredGFPOnlyObjects.size()
)

println(
    "Restored under annotation: " +
    ROI_NAME
)

println ""


// ============================================================
// GET FINAL DETECTIONS FROM CURRENT ROI ONLY
// ============================================================

def finalDetections =
    getCurrentROIDetections()


// ============================================================
// FINAL COUNTS
// ============================================================

def finalStarter =
    finalDetections.count {

        def pc =
            it.getPathClass()

        pc != null &&
        pc.getName() == "Starter"
    }


def finalGFPOnly =
    finalDetections.count {

        def pc =
            it.getPathClass()

        pc != null &&
        pc.getName() == "GFP-only"
    }


def finalRV =
    finalDetections.count {

        def pc =
            it.getPathClass()

        pc != null &&
        pc.getName() == "RV+"
    }


// ============================================================
// FINAL RESULTS
// ============================================================

println ""
println "================================================"
println " FINAL TVA/RABIES RESULTS"
println "================================================"
println ""

println "ROI:                 " + ROI_NAME
println ""

println "Starter:             " + finalStarter
println "GFP-only:            " + finalGFPOnly
println "RV+:                 " + finalRV

println ""

println "Total detections:    " +
        finalDetections.size()

println "================================================"
println ""


// ============================================================
// EXPORT
// ============================================================
//
// A separate CSV is created for each ROI.
//
// Example:
//
// TVA_Rabies_Cell_Analysis_v9_SUB.csv
// TVA_Rabies_Cell_Analysis_v9_CA1.csv
// TVA_Rabies_Cell_Analysis_v9_CA3.csv
//
// This prevents later ROI analyses from overwriting earlier
// results.
//

if (EXPORT_RESULTS) {


    println "Exporting measurements..."
    println ""


    def safeROIName =
        ROI_NAME.replaceAll(
            "[^A-Za-z0-9_-]",
            "_"
        )


    def outputFile =
        buildFilePath(
            PROJECT_BASE_DIR,
            EXPORT_FILENAME +
            "_" +
            safeROIName +
            ".csv"
        )


    def writer =
        new File(outputFile).newPrintWriter()


    // --------------------------------------------------------
    // Header
    // --------------------------------------------------------

    writer.println(
        "ROI," +
        "Classification," +
        "Centroid_X_px," +
        "Centroid_Y_px," +
        "Cell_Area_px2," +
        "GFP_Nucleus_Mean," +
        "GFP_Cell_Mean," +
        "mCherry_Nucleus_Mean," +
        "mCherry_Cell_Mean"
    )


    // --------------------------------------------------------
    // Rows
    // --------------------------------------------------------

    finalDetections.each { detection ->


        def ml =
            detection.getMeasurementList()


        def roi =
            detection.getROI()


        def x =
            roi.getCentroidX()


        def y =
            roi.getCentroidY()


        def area =
            roi.getArea()


        def gfpNucleus =
            ml.get(
                "Nucleus: ${GFP_CHANNEL} mean"
            )


        def gfpCell =
            ml.get(
                "Cell: ${GFP_CHANNEL} mean"
            )


        def mcherryNucleus =
            ml.get(
                "Nucleus: ${MCHERRY_CHANNEL} mean"
            )


        def mcherryCell =
            ml.get(
                "Cell: ${MCHERRY_CHANNEL} mean"
            )


        def pathClass =
            detection.getPathClass()


        def classification =
            pathClass == null ?
            "" :
            pathClass.getName()


        writer.println(
            '"' + ROI_NAME + '",' +
            '"' + classification + '",' +
            x + ',' +
            y + ',' +
            area + ',' +
            gfpNucleus + ',' +
            gfpCell + ',' +
            mcherryNucleus + ',' +
            mcherryCell
        )
    }


    writer.close()


    println "Results exported to:"
    println outputFile
    println ""
}


// ============================================================
// COMPLETE
// ============================================================

println ""
println "================================================"
println "Analysis complete."
println "ROI: " + ROI_NAME
println "================================================"