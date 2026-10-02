/**
 * ============================================================
 * TVA/Rabies two-pass cell detection - v7
 * ============================================================
 *
 * QuPath 0.7.0
 *
 * MULTI-ROI VERSION
 *
 * Select ONE annotation/ROI before running the script.
 *
 * The script analyzes only that ROI and preserves detections
 * belonging to all other ROIs.
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
 * are temporarily removed.
 *
 * Detections in other annotations/ROIs are not touched.
 *
 * GFP-only cells are saved before Pass 2 and reconstructed as
 * CELL objects under the original annotation.
 *
 * ============================================================
 */


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

println ""
println "Analysis annotation: " + ROI_NAME
println ""


// ------------------------------------------------------------
// PASS 1 - GFP DETECTION
// ------------------------------------------------------------

def GFP_PIXEL_SIZE = 1.0

def GFP_BACKGROUND_RADIUS = 10.0
def GFP_MEDIAN_RADIUS = 2.0
def GFP_SIGMA = 5.0

def GFP_MIN_AREA = 6.0
def GFP_MAX_AREA = 600.0

def GFP_THRESHOLD = 200.0

def GFP_CELL_EXPANSION = 3.0


// ------------------------------------------------------------
// PASS 1 - mCherry threshold for Starter
// ------------------------------------------------------------

def MCHERRY_STARTER_THRESHOLD = 2000.0


// ------------------------------------------------------------
// PASS 2 - mCherry DETECTION
// ------------------------------------------------------------

def MCHERRY_PIXEL_SIZE = 1.0

def MCHERRY_BACKGROUND_RADIUS = 10.0
def MCHERRY_MEDIAN_RADIUS = 2.0
def MCHERRY_SIGMA = 5.0

def MCHERRY_MIN_AREA = 10.0
def MCHERRY_MAX_AREA = 600.0

def MCHERRY_THRESHOLD = 600.0

def MCHERRY_CELL_EXPANSION = 3.0


// ------------------------------------------------------------
// SPATIAL MATCHING
// ------------------------------------------------------------

def STARTER_MATCH_DISTANCE_MICRONS = 5.0


// ------------------------------------------------------------
// EXPORT
// ------------------------------------------------------------

def EXPORT_RESULTS = true

def EXPORT_FILENAME = "TVA_Rabies_Cell_Analysis_v7"


// ============================================================
// INITIAL CHECK
// ============================================================

def imageData = getCurrentImageData()

if (imageData == null) {
    throw new Exception("No image is currently open.")
}


println ""
println "================================================"
println " TVA/Rabies two-pass cell analysis - v7"
println " QuPath 0.7.0"
println "================================================"
println ""

println "Analysis ROI:    " + ROI_NAME
println "GFP channel:     " + GFP_CHANNEL
println "mCherry channel: " + MCHERRY_CHANNEL
println ""


// ============================================================
// REMOVE ONLY EXISTING DETECTIONS IN CURRENT ROI
// ============================================================
//
// IMPORTANT:
// Do NOT use removeDetections() here because that removes
// detections outside the current ROI as well.
//
// We remove only detection children of analysisParent.
//

def existingDetections =
        analysisParent.getChildObjects().findAll {
            it.isDetection()
        }

if (!existingDetections.isEmpty()) {

    println "Removing existing detections from current ROI: " +
            existingDetections.size()

    analysisParent.removeChildObjects(existingDetections)

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
// GET ONLY PASS-1 DETECTIONS FROM CURRENT ROI
// ============================================================

def gfpDetections =
        analysisParent.getChildObjects().findAll {
            it.isDetection()
        }

println "Number of GFP-derived detections in " +
        ROI_NAME + ": " +
        gfpDetections.size()

println ""


// ============================================================
// SAVE PASS-1 DATA
// ============================================================

def pass1Data = []

def gfpStarterData = []

def gfpOnlyData = []


println "GFP detection measurements:"
println "----------------------------------------------"


gfpDetections.eachWithIndex { detection, index ->

    def ml = detection.getMeasurementList()


    def gfpNucleus =
            ml.get("Nucleus: ${GFP_CHANNEL} mean")

    def gfpCell =
            ml.get("Cell: ${GFP_CHANNEL} mean")

    def mcherryNucleus =
            ml.get("Nucleus: ${MCHERRY_CHANNEL} mean")

    def mcherryCell =
            ml.get("Cell: ${MCHERRY_CHANNEL} mean")


    def roi = detection.getROI()

    def nucleusROI = null

    // Try to preserve the nucleus ROI.
    // If unavailable, the script continues with null.
    try {
        nucleusROI = detection.getNucleusROI()
    } catch (Exception e) {
        nucleusROI = null
    }


    def centroidX = roi.getCentroidX()
    def centroidY = roi.getCentroidY()


    println(
        "${index + 1}: " +
        "GFP nuclear mean = ${gfpNucleus}, " +
        "mCherry cell mean = ${mcherryCell}"
    )


    // --------------------------------------------------------
    // Biological classification
    // --------------------------------------------------------

    def classification

    if (mcherryCell != null &&
            !Double.isNaN(mcherryCell) &&
            mcherryCell >= MCHERRY_STARTER_THRESHOLD) {

        classification = "Starter"

    } else {

        classification = "GFP-only"
    }


    // --------------------------------------------------------
    // Save complete Pass-1 information
    // --------------------------------------------------------

    def record = [

        roi: roi,

        nucleusROI: nucleusROI,

        classification: classification,

        centroidX: centroidX,
        centroidY: centroidY,

        gfpNucleus: gfpNucleus,
        gfpCell: gfpCell,

        mcherryNucleus: mcherryNucleus,
        mcherryCell: mcherryCell,

        cellArea: roi.getArea(),

        roiName: ROI_NAME
    ]


    pass1Data << record


    if (classification == "Starter") {

        gfpStarterData << record

    } else {

        gfpOnlyData << record
    }


    detection.setClassification(classification)
}


println ""

println "Starter cells from GFP pass: " +
        gfpStarterData.size()

println "GFP-only cells:              " +
        gfpOnlyData.size()

println ""


// ============================================================
// SAVE GFP-ONLY DATA BEFORE PASS 2
// ============================================================

println "----------------------------------------------"
println "Saving Pass-1 GFP-only objects"
println "----------------------------------------------"

println "GFP-only objects saved: " +
        gfpOnlyData.size()

println "GFP Starter objects saved for spatial matching: " +
        gfpStarterData.size()

println ""


// ============================================================
// REMOVE ONLY PASS-1 DETECTIONS FROM CURRENT ROI
// ============================================================
//
// This is the critical multi-ROI change.
//
// Do NOT use:
//     removeDetections()
//
// because that can affect detections belonging to other ROIs.
//

def pass1CurrentROIdetections =
        analysisParent.getChildObjects().findAll {
            it.isDetection()
        }

if (!pass1CurrentROIdetections.isEmpty()) {

    analysisParent.removeChildObjects(
        pass1CurrentROIdetections
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
// GET ONLY PASS-2 DETECTIONS FROM CURRENT ROI
// ============================================================

def mcherryDetections =
        analysisParent.getChildObjects().findAll {
            it.isDetection()
        }

println "Number of mCherry-derived detections in " +
        ROI_NAME + ": " +
        mcherryDetections.size()

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

def starterMatchDistancePixels =
        STARTER_MATCH_DISTANCE_MICRONS /
        pixelSizeMicrons


println "Starter matching distance: " +
        STARTER_MATCH_DISTANCE_MICRONS +
        " µm"

println "Approximate pixel distance: " +
        starterMatchDistancePixels

println ""


// ============================================================
// MATCH EACH mCHERRY DETECTION TO SAVED GFP STARTERS
// ============================================================

mcherryDetections.eachWithIndex { detection, index ->

    def roi = detection.getROI()

    def mcherryX = roi.getCentroidX()
    def mcherryY = roi.getCentroidY()


    def nearestDistancePixels =
            Double.POSITIVE_INFINITY

    def matchedStarter = null


    // --------------------------------------------------------
    // Compare with saved GFP Starters
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


        if (distancePixels <
                nearestDistancePixels) {

            nearestDistancePixels =
                    distancePixels

            matchedStarter =
                    starter
        }
    }


    // --------------------------------------------------------
    // Classification
    // --------------------------------------------------------

    if (matchedStarter != null &&
            nearestDistancePixels <=
                    starterMatchDistancePixels) {

        detection.setClassification("Starter")

        starterCount++


        println(
            "mCherry detection ${index + 1}: " +
            "STARTER " +
            "(distance = " +
            (nearestDistancePixels *
             pixelSizeMicrons) +
            " µm)"
        )

    } else {

        detection.setClassification("RV+")

        rvCount++


        println(
            "mCherry detection ${index + 1}: RV+"
        )
    }


    // --------------------------------------------------------
    // Add ROI metadata to Pass-2 object
    // --------------------------------------------------------

    def ml = detection.getMeasurementList()

    ml.put(
        "TVA_Rabies_ROI_ID",
        analysisParent.getID()
    )

    ml.put(
        "TVA_Rabies_Analysis",
        1.0
    )
}


println ""

println "Starter cells identified in Pass 2: " +
        starterCount

println "RV+ cells identified in Pass 2: " +
        rvCount

println ""


// ============================================================
// RESTORE GFP-ONLY OBJECTS
// ============================================================
//
// IMPORTANT:
// Create CELL objects rather than generic DETECTION objects.
//
// This should make the hierarchy display:
//
//     Cell (GFP-only)
//
// rather than:
//
//     Detection (GFP-only)
//

println "----------------------------------------------"
println "Restoring GFP-only cells"
println "----------------------------------------------"


def restoredGFPOnlyObjects = []


gfpOnlyData.each { record ->


    def restoredObject


    // --------------------------------------------------------
    // Try to recreate a true Cell object using both the cell
    // ROI and the original nucleus ROI.
    // --------------------------------------------------------

    if (record.nucleusROI != null) {

        restoredObject =
                PathObjects.createCellObject(
                    record.roi,
                    record.nucleusROI,
                    PathClass.fromString("GFP-only")
                )

    } else {

        // Fallback if no nucleus ROI is available.
        restoredObject =
                PathObjects.createCellObject(
                    record.roi,
                    null,
                    PathClass.fromString("GFP-only")
                )
    }


    // --------------------------------------------------------
    // Restore measurements
    // --------------------------------------------------------

    def ml =
            restoredObject.getMeasurementList()


    if (record.gfpNucleus != null &&
            !Double.isNaN(record.gfpNucleus)) {

        ml.put(
            "Nucleus: ${GFP_CHANNEL} mean",
            record.gfpNucleus
        )
    }


    if (record.gfpCell != null &&
            !Double.isNaN(record.gfpCell)) {

        ml.put(
            "Cell: ${GFP_CHANNEL} mean",
            record.gfpCell
        )
    }


    if (record.mcherryNucleus != null &&
            !Double.isNaN(record.mcherryNucleus)) {

        ml.put(
            "Nucleus: ${MCHERRY_CHANNEL} mean",
            record.mcherryNucleus
        )
    }


    if (record.mcherryCell != null &&
            !Double.isNaN(record.mcherryCell)) {

        ml.put(
            "Cell: ${MCHERRY_CHANNEL} mean",
            record.mcherryCell
        )
    }


    // --------------------------------------------------------
    // ROI metadata
    // --------------------------------------------------------

    ml.put(
        "TVA_Rabies_ROI_ID",
        analysisParent.getID()
    )

    ml.put(
        "TVA_Rabies_Pass1_GFPonly",
        1.0
    )


    restoredGFPOnlyObjects << restoredObject
}


// ============================================================
// RESTORE GFP-ONLY CELLS UNDER CURRENT ROI
// ============================================================

analysisParent.addChildObjects(
    restoredGFPOnlyObjects
)

fireHierarchyUpdate()


println "GFP-only cells restored: " +
        restoredGFPOnlyObjects.size()

println "Restored under annotation: " +
        ROI_NAME

println ""


// ============================================================
// FINAL DETECTIONS IN CURRENT ROI
// ============================================================

def finalDetections =
        analysisParent.getChildObjects().findAll {
            it.isDetection()
        }


def finalStarter =
        finalDetections.count {
            it.getClassification() == "Starter"
        }


def finalGFPOnly =
        finalDetections.count {
            it.getClassification() == "GFP-only"
        }


def finalRV =
        finalDetections.count {
            it.getClassification() == "RV+"
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

println "Starter:             " +
        finalStarter

println "GFP-only:            " +
        finalGFPOnly

println "RV+:                 " +
        finalRV

println ""

println "Total detections in ROI: " +
        finalDetections.size()

println "================================================"
println ""


// ============================================================
// EXPORT
// ============================================================
//
// Rather than relying exclusively on saveDetectionMeasurements(),
// create a CSV that explicitly contains the ROI name.
//
// This makes the file immediately useful for combining results
// from SUB, CA1, CA3, DG, etc.
//

if (EXPORT_RESULTS) {

    println "Exporting measurements..."
    println ""


    def outputFile =
            buildFilePath(
                PROJECT_BASE_DIR,
                EXPORT_FILENAME + "_" +
                ROI_NAME.replaceAll(
                    "[^A-Za-z0-9_-]",
                    "_"
                ) +
                ".csv"
            )


    def writer =
            new File(outputFile).newPrintWriter()


    // --------------------------------------------------------
    // CSV header
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
    // CSV rows
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


        def classification =
                detection.getClassification()


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

println "Analysis complete."
println ""

println "IMPORTANT:"
println "Previous ROI detections were preserved."
println "Current ROI: " + ROI_NAME
println ""