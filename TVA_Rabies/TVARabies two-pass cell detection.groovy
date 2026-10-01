// Script name: TVA_Rabies_pipeline_v1
// Author: Fred Manseau
// Laboratory: Sylvain Williams Lab
// Institution: Douglas Research Centre / McGill University
// Created: 2026-08-20
// Last modified:
// Description: QuPath pipeline for analysing TVA-Rabies histological tracing experiments

/**
 * TVA/Rabies two-pass cell detection
 *
 * QuPath 0.7.0
 *
 * PASS 1:
 *   GFP nuclear detection
 *   -> cell expansion
 *   -> measure cellular mCherry
 *   -> GFP+ / mCherry+ = Starter
 *   -> GFP+ / mCherry- = GFP-only
 *
 * PASS 2:
 *   mCherry-positive cell detection
 *   -> cell expansion
 *   -> measure GFP
 *   -> mCherry+ / GFP- = RV+
 *   -> mCherry+ / GFP+ = Starter candidate
 *
 * The script keeps the two detection populations separate.
 *
 * IMPORTANT:
 *   Adjust the channel names and thresholds below to match
 *   the image server exactly.
 */


// ============================================================
// USER PARAMETERS
// ============================================================

// ---------- CHANNEL NAMES ----------

def GFP_CHANNEL = "GFP"
def MCHERRY_CHANNEL = "mCherry"


// ---------- PASS 1: GFP DETECTION ----------

def GFP_PIXEL_SIZE = 1.0
def GFP_BACKGROUND_RADIUS = 10.0
def GFP_MEDIAN_RADIUS = 2.0
def GFP_SIGMA = 5.0
def GFP_MIN_AREA = 10.0
def GFP_MAX_AREA = 600.0
def GFP_THRESHOLD = 600.0
def GFP_CELL_EXPANSION = 3.0


// mCherry threshold for GFP-derived cells
def MCHERRY_STARTER_THRESHOLD = 1500.0


// ---------- PASS 2: mCherry DETECTION ----------

def MCHERRY_PIXEL_SIZE = 1.0
def MCHERRY_BACKGROUND_RADIUS = 10.0
def MCHERRY_MEDIAN_RADIUS = 2.0
def MCHERRY_SIGMA = 5.0
def MCHERRY_MIN_AREA = 10.0
def MCHERRY_MAX_AREA = 600.0
def MCHERRY_THRESHOLD = 600.0
def MCHERRY_CELL_EXPANSION = 3.0


// GFP threshold used to distinguish
// mCherry+ / GFP+ from mCherry+ / GFP-
def GFP_POSITIVE_THRESHOLD = 1500.0


// ---------- EXPORT ----------

def EXPORT_RESULTS = true

// Change this to wherever you want the results saved.
// The script will create a filename automatically.
def EXPORT_FILENAME = "TVA_Rabies_Cell_Analysis"


// ============================================================
// CHECK IMAGE
// ============================================================

def imageData = getCurrentImageData()

if (imageData == null) {
    throw new Exception("No image is currently open.")
}

println ""
println "=============================================="
println " TVA/Rabies two-pass cell analysis"
println " QuPath 0.7.0"
println "=============================================="
println ""

println "GFP channel:     " + GFP_CHANNEL
println "mCherry channel: " + MCHERRY_CHANNEL
println ""


// ============================================================
// REMOVE EXISTING DETECTIONS
// ============================================================
//
// This prevents previous detection runs from being mixed with
// the current analysis.
//

clearDetections()

println "Existing detections cleared."
println ""


// ============================================================
// PASS 1
// GFP NUCLEAR DETECTION
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
// SAVE PASS 1 DETECTIONS
// ============================================================
//
// We need to retain these objects before starting the second
// detection pass.
//

def gfpDetections = getDetectionObjects()

println "Number of GFP-derived detections: " + gfpDetections.size()
println ""


// ============================================================
// CLASSIFY PASS 1
// ============================================================
//
// Positive cells are GFP+ / mCherry+ = Starter
//
// Negative cells are GFP+ / mCherry- = GFP-only
//
// We don't rely exclusively on QuPath's Positive/Negative
// classification; we explicitly assign our own biological
// classifications.
//

def starterCountPass1 = 0
def gfpOnlyCount = 0

gfpDetections.each { detection ->

    def ml = detection.getMeasurementList()

    def mcherry = ml.get("Cell: ${MCHERRY_CHANNEL} mean")

    if (mcherry != null &&
            !Double.isNaN(mcherry) &&
            mcherry >= MCHERRY_STARTER_THRESHOLD) {

        detection.setClassification(
                PathClass.fromString("Starter")
        )

        starterCountPass1++

    } else {

        detection.setClassification(
                PathClass.fromString("GFP-only")
        )

        gfpOnlyCount++
    }
}


println "Starter cells from GFP pass: " + starterCountPass1
println "GFP-only cells:              " + gfpOnlyCount
println ""


// ============================================================
// STORE PASS 1
// ============================================================
//
// We cannot simply leave these detections in place while
// performing a second independent detection.
//
// Store them temporarily.
//

def pass1Objects = new ArrayList(gfpDetections)


// ============================================================
// CLEAR DETECTIONS FOR PASS 2
// ============================================================

clearDetections()

println "GFP detections temporarily removed."
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
            "thresholdPositive1":${GFP_POSITIVE_THRESHOLD},
            "thresholdPositive2":${GFP_POSITIVE_THRESHOLD},
            "thresholdPositive3":${GFP_POSITIVE_THRESHOLD},
            "singleThreshold":true
        }"""
)

println "mCherry detection completed."
println ""


def mcherryDetections = getDetectionObjects()

println "Number of mCherry-derived detections: " +
        mcherryDetections.size()

println ""


// ============================================================
// CLASSIFY PASS 2
// ============================================================
//
// mCherry+ / GFP- = RV+
//
// mCherry+ / GFP+ = Starter candidate
//
// IMPORTANT:
// These Starter candidates are NOT added to the Starter
// population from Pass 1. This prevents double counting.
//

def rvCount = 0
def starterCandidateCount = 0

mcherryDetections.each { detection ->

    def ml = detection.getMeasurementList()

    def gfp = ml.get("Cell: ${GFP_CHANNEL} mean")

    if (gfp != null &&
            !Double.isNaN(gfp) &&
            gfp >= GFP_POSITIVE_THRESHOLD) {

        detection.setClassification(
                PathClass.fromString("Starter_candidate")
        )

        starterCandidateCount++

    } else {

        detection.setClassification(
                PathClass.fromString("RV+")
        )

        rvCount++
    }
}


println "RV+ cells:                  " + rvCount
println "Starter candidates:         " + starterCandidateCount
println ""


// ============================================================
// ADD A POPULATION MEASUREMENT
// ============================================================
//
// We add a numerical code that can be useful when exporting
// the measurements:
//
// 1 = Starter
// 2 = GFP-only
// 3 = RV+
// 4 = Starter candidate
//

pass1Objects.each { detection ->

    def ml = detection.getMeasurementList()

    def classification = detection.getPathClass()

    if (classification != null &&
            classification.toString() == "Starter") {

        ml.put("TVA_Rabies_Population", 1.0)

    } else {

        ml.put("TVA_Rabies_Population", 2.0)
    }
}


mcherryDetections.each { detection ->

    def ml = detection.getMeasurementList()

    def classification = detection.getPathClass()

    if (classification != null &&
            classification.toString() == "RV+") {

        ml.put("TVA_Rabies_Population", 3.0)

    } else {

        ml.put("TVA_Rabies_Population", 4.0)
    }
}


// ============================================================
// RESTORE BOTH POPULATIONS
// ============================================================
//
// This is the difficult part of a two-pass workflow.
//
// We need both sets of detections available in the image
// hierarchy after the second detection pass.
//
// NOTE:
// QuPath detection objects belong to the hierarchy, so simply
// storing them in a Groovy list does not necessarily restore
// them to the hierarchy. We therefore add them back to the
// current parent object.
//

def parent = getCurrentParent()

if (parent == null) {
    throw new Exception(
            "Please select the annotation/parent object containing the analysis region."
    )
}


// Remove the current mCherry detections from the hierarchy
clearDetections()


// Add both detection populations
addObjects(pass1Objects)
addObjects(mcherryDetections)

fireHierarchyUpdate()


// ============================================================
// FINAL COUNTS
// ============================================================

def allDetections = getDetectionObjects()

def finalStarter = allDetections.count {
    it.getPathClass()?.toString() == "Starter"
}

def finalGFPOnly = allDetections.count {
    it.getPathClass()?.toString() == "GFP-only"
}

def finalRV = allDetections.count {
    it.getPathClass()?.toString() == "RV+"
}

def finalStarterCandidates = allDetections.count {
    it.getPathClass()?.toString() == "Starter_candidate"
}


println ""
println "=============================================="
println " FINAL TVA/RABIES RESULTS"
println "=============================================="
println ""
println "Starter:             " + finalStarter
println "GFP-only:            " + finalGFPOnly
println "RV+:                 " + finalRV
println "Starter candidate:   " + finalStarterCandidates
println ""
println "Total detections:    " + allDetections.size()
println "=============================================="
println ""


// ============================================================
// EXPORT
// ============================================================

if (EXPORT_RESULTS) {

    println "Exporting measurements..."

    def outputFile = buildFilePath(
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


print "Analysis complete."