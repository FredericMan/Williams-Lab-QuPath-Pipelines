/**
 * ============================================================
 * TVA/Rabies two-pass cell detection - v6
 * ============================================================
 *
 * QuPath 0.7.0
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
 * overlaps GFP Starter = Starter
 * no GFP Starter overlap = RV+
 *
 *
 * FINAL
 * -----
 * Starter
 * GFP-only
 * RV+
 *
 *
 * IMPORTANT
 * ----------
 * Pass-1 detections are saved in memory before Pass 2.
 *
 * All Pass-1 detections are then removed so that QuPath can
 * perform the second PositiveCellDetection normally.
 *
 * After Pass 2, the saved GFP-only detections are recreated
 * and added back to the hierarchy.
 *
 * The original GFP Starter detections are NOT recreated,
 * because the corresponding mCherry detections are retained
 * as the final Starter objects.
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
// PASS 1 - GFP DETECTION
// ------------------------------------------------------------

def GFP_PIXEL_SIZE = 1.0

def GFP_BACKGROUND_RADIUS = 10.0
def GFP_MEDIAN_RADIUS = 2.0
def GFP_SIGMA = 5.0

def GFP_MIN_AREA = 10.0
def GFP_MAX_AREA = 600.0

def GFP_THRESHOLD = 200.0

def GFP_CELL_EXPANSION = 3.0


// ------------------------------------------------------------
// PASS 1 - mCherry threshold for Starter
// ------------------------------------------------------------
//
// A GFP-derived cell is considered a Starter if its
// cellular TRITC mean is >= this value.
//

def MCHERRY_STARTER_THRESHOLD = 1500.0


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
//
// Maximum distance between the centroid of a mCherry-derived
// detection and the centroid of a saved GFP Starter detection.
//
// This is expressed in microns.
//
// Because GFP is nuclear and mCherry is cellular, we do not
// require the two centroids to be identical.
//

def STARTER_MATCH_DISTANCE_MICRONS = 5.0


// ------------------------------------------------------------
// EXPORT
// ------------------------------------------------------------

def EXPORT_RESULTS = true

def EXPORT_FILENAME = "TVA_Rabies_Cell_Analysis_v5"


// ============================================================
// INITIAL CHECK
// ============================================================

def imageData = getCurrentImageData()

if (imageData == null) {
    throw new Exception("No image is currently open.")
}


println ""
println "================================================"
println " TVA/Rabies two-pass cell analysis - v5"
println " QuPath 0.7.0"
println "================================================"
println ""

println "GFP channel:     " + GFP_CHANNEL
println "mCherry channel: " + MCHERRY_CHANNEL
println ""


// ============================================================
// REMOVE ANY EXISTING DETECTIONS
// ============================================================

removeDetections()

println "Existing detections cleared."
println ""


// ============================================================
// PASS 1
// GFP DETECTION
// ============================================================

println "----------------------------------------------"
println "PASS 1: GFP detection"
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
// GET PASS-1 DETECTIONS
// ============================================================

def gfpDetections = getDetectionObjects()

println "Number of GFP-derived detections: " +
        gfpDetections.size()

println ""


// ============================================================
// SAVE PASS-1 DATA
// ============================================================
//
// We create independent records containing the information
// we need after Pass 2.
//
// This is important because the actual QuPath detection
// objects will be removed before Pass 2.
//

def pass1Data = []

def gfpStarterData = []

def gfpOnlyData = []


// ------------------------------------------------------------
// Diagnostic output
// ------------------------------------------------------------

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

    def centroidX = roi.getCentroidX()
    def centroidY = roi.getCentroidY()


    println(
        "${index + 1}: " +
        "GFP nuclear mean = ${gfpNucleus}, " +
        "mCherry cell mean = ${mcherryCell}"
    )


    // --------------------------------------------------------
    // Determine biological classification
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

        classification: classification,

        centroidX: centroidX,
        centroidY: centroidY,

        gfpNucleus: gfpNucleus,
        gfpCell: gfpCell,

        mcherryNucleus: mcherryNucleus,
        mcherryCell: mcherryCell,

        cellArea: roi.getArea()
    ]


    pass1Data << record


    if (classification == "Starter") {

        gfpStarterData << record

    } else {

        gfpOnlyData << record
    }


    // --------------------------------------------------------
    // Assign classification to original object
    // --------------------------------------------------------

    detection.setClassification(classification)
}


println ""


println "Starter cells from GFP pass: " +
        gfpStarterData.size()

println "GFP-only cells:              " +
        gfpOnlyData.size()

println ""


// ============================================================
// IMPORTANT: SAVE GFP-ONLY DATA BEFORE PASS 2
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
// REMOVE PASS-1 DETECTIONS
// ============================================================
//
// We now let QuPath run the mCherry detection normally.
//
// The GFP information is safely stored in pass1Data,
// gfpStarterData and gfpOnlyData.
//

removeDetections()

println "Pass-1 detections removed from hierarchy."
println "Pass-1 data retained in memory."
println ""


// ============================================================
// PASS 2
// mCHERRY DETECTION
// ============================================================

println "----------------------------------------------"
println "PASS 2: mCherry detection"
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
// GET PASS-2 DETECTIONS
// ============================================================

def mcherryDetections = getDetectionObjects()

println "Number of mCherry-derived detections: " +
        mcherryDetections.size()

println ""


// ============================================================
// PASS 2 CLASSIFICATION
// ============================================================
//
// For each mCherry detection:
//
// 1. Determine whether it is spatially associated with a
//    previously detected GFP Starter.
//
// 2. If yes:
//       Starter
//
// 3. If no:
//       RV+
//
// IMPORTANT:
// GFP-only cells are NOT used for the Starter matching step.
// They were GFP+ / mCherry- in Pass 1.
//

def rvCount = 0
def starterCount = 0


// ------------------------------------------------------------
// Pixel calibration
// ------------------------------------------------------------
//
// Convert the user-defined matching distance in microns to
// pixels.
//

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


    def nearestDistancePixels = Double.POSITIVE_INFINITY

    def matchedStarter = null


    // --------------------------------------------------------
    // Compare with every saved GFP Starter
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
            (nearestDistancePixels * pixelSizeMicrons) +
            " µm)"
        )

    } else {

        detection.setClassification("RV+")

        rvCount++


        println(
            "mCherry detection ${index + 1}: " +
            "RV+"
        )
    }
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
// The mCherry detection pass replaced the previous detection
// objects. We now recreate ONLY the GFP-only cells saved from
// Pass 1.
//
// The original Starter objects are deliberately NOT restored.
// The corresponding mCherry detections are now the final
// Starter objects.
//

println "----------------------------------------------"
println "Restoring GFP-only objects"
println "----------------------------------------------"


def restoredGFPOnlyObjects = []


gfpOnlyData.each { record ->

    def restoredObject =
            PathObjects.createDetectionObject(
                record.roi,
                PathClass.fromString("GFP-only")
            )


    // --------------------------------------------------------
    // Restore important measurements
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


    ml.put(
        "TVA_Rabies_Pass1_GFPonly",
        1.0
    )


    restoredGFPOnlyObjects << restoredObject
}


// ------------------------------------------------------------
// Add restored objects to hierarchy
// ------------------------------------------------------------

addObjects(restoredGFPOnlyObjects)

fireHierarchyUpdate()


println "GFP-only objects restored: " +
        restoredGFPOnlyObjects.size()

println ""


// ============================================================
// FINAL COUNTS
// ============================================================

def finalDetections =
        getDetectionObjects()


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

println "Starter:             " +
        finalStarter

println "GFP-only:            " +
        finalGFPOnly

println "RV+:                 " +
        finalRV

println ""

println "Total detections:    " +
        finalDetections.size()

println "================================================"
println ""


// ============================================================
// EXPORT
// ============================================================

if (EXPORT_RESULTS) {

    println "Exporting measurements..."
    println ""


    def outputFile =
            buildFilePath(
                PROJECT_BASE_DIR,
                EXPORT_FILENAME + ".csv"
            )


    saveDetectionMeasurements(
        outputFile
    )


    println "Results exported to:"
    println outputFile
    println ""
}


// ============================================================
// COMPLETE
// ============================================================

println "Analysis complete."