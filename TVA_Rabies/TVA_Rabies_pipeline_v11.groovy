// Script name: TVA_Rabies_pipeline_v11
// Author: Fred Manseau
// Laboratory: Sylvain Williams Lab
// Institution: Douglas Research Centre / McGill University
// Created: 2026-08-19
// Last modified: 2026-09-30
// Description: QuPath pipeline for analysing TVA-Rabies histological tracing experiments
// Development history and version tracking:
// See the associated Git repository:


/**
 * ============================================================
 * TVA/Rabies two-pass cell analysis - v11
 * ============================================================
 *
 * QuPath 0.7.0
 *
 * MULTI-ROI / MULTI-SLICE VERSION
 *
 * Select ONE annotation/ROI before running the script.
 *
 * v11 CHANGELOG (relative to v10_chan_name)
 * -------------------------------------------------------
 *   1. ANIMAL_ID added for unambiguous per-animal traceability,
 *      ported from the companion Synaptag pipeline
 *      (mRuby2_SypEGFP_two-pass_v14_Atlas_Hemisphere) with the same
 *      column name, position, and derivation. Inferred from the
 *      folder containing the QuPath project's base directory (NOT
 *      from PROJECT_BASE_DIR as a string, which can hold an
 *      unresolved token - see the ANIMAL ID section below). Printed
 *      as a standing diagnostic every run, and added as a new column
 *      to all three CSV exports (object, ROI summary, and the
 *      TVA_Rabies-specific Excluded_Artifacts file) immediately
 *      after Image_Name, plus the run-parameters log.
 *   2. Sectioning-angle tilt correction reset to identity: the CCFv3
 *      -> stereotaxic rotation (was 5 deg) and DV squeeze (was
 *      0.9434) were identified as a mistaken, redundant correction -
 *      they were originally entered to compensate for an X/Y
 *      cutting-angle tilt that ABBA's own Atlas Slicing -> X and Y
 *      Rotation (deg) step already corrects, in 3D, before
 *      coordinates are ever exported. Both are now 0/1.0 by default,
 *      renamed/grouped into a dedicated SECTIONING-ANGLE CORRECTION
 *      section (aliased back to the original constant names so nothing
 *      downstream needed to change), clearly separated from the
 *      still-unverified X/Y/Z origin offsets, which are unchanged.
 *      TRANSFORM_VERSION updated accordingly.
 *
 * UNDOCUMENTED FIXES CARRIED OVER FROM v10_chan_name (retroactively
 * logged here, since that file never received its own changelog entry)
 * -------------------------------------------------------
 *   - Channel name resolution: GFP_CHANNEL/MCHERRY_CHANNEL are now
 *     resolved against the actual channel names of whichever image
 *     is open (exact match, then case-insensitive prefix match)
 *     before Pass 1 begins, with a loud fail-fast exception if
 *     neither matches. Fixes a real bug where QuPath's automatic
 *     "(Cn)" channel-name disambiguation suffix (applied per-image,
 *     inconsistently) caused Pass 2 to silently detect 0 cells for
 *     at least one image, with the failure only visible in QuPath's
 *     own Log panel (runPlugin() catches this internally rather than
 *     throwing a catchable Groovy exception).
 *   - The current ROI is now explicitly re-selected
 *     (setSelectedObject(analysisParent)) immediately before EACH of
 *     the two detection passes, with a console line confirming the
 *     selection, since runPlugin() acts on whatever is currently
 *     selected in the viewer rather than taking an explicit target.
 *
 * v10 CHANGELOG (relative to v9 + interactive QC v2/v3)
 * -------------------------------------------------------
 *   1. Starter reconciliation: GFP+/mCherry+ candidates from Pass 1
 *      that have no matching independent Pass-2 detection are no
 *      longer silently dropped - they are restored as Starter cells
 *      (Starter_Source_Code = 2), and a mismatch triggers an
 *      explicit console warning.
 *   2. Pass-2 matching is now a one-to-one greedy assignment
 *      (closest pairs claimed first), preventing one GFP Starter
 *      from being claimed by two Pass-2 detections.
 *   3. The CCFv3->stereotaxic transform is isolated into a single
 *      function with a documented (and currently UNVERIFIED)
 *      provenance, an automatic per-detection plausibility check
 *      (Stereotaxic_Plausible), and a TRANSFORM_VERSION tag written
 *      to every export. A standalone companion script,
 *      CCFv3_Stereotaxic_Transform_Calibration_v1.groovy, is
 *      provided to refit these constants from landmark data.
 *   4. A complete run-parameters log (every tunable constant, plus
 *      per-ROI signal/background QC statistics and warnings) is now
 *      written alongside the object/summary CSVs. Manual QC now
 *      supports marking a detection "Artifact" (recording its prior
 *      classification) instead of deleting it outright; artifacts
 *      are excluded from the main results but documented in a
 *      separate *_Excluded_Artifacts.csv.
 *
 * WORKFLOW
 * --------
 * PASS 1:
 *   FITC/GFP nuclear detection
 *      -> cell expansion
 *      -> measure TRITC/mCherry
 *      -> Starter or GFP-only
 *
 * PASS 2:
 *   TRITC/mCherry detection
 *      -> spatially match to saved GFP Starter cells (one-to-one)
 *      -> Starter or RV+
 *
 * RESTORATION:
 *   Restore GFP-only cells, AND any unmatched GFP Starter
 *   candidates (see v10 changelog item 1), as true QuPath Cell
 *   objects under the original selected annotation.
 *
 * COORDINATES:
 *   Calculate Allen CCFv3 X/Y/Z
 *   Calculate stereotaxic AP/ML/DV using the CCFv3 -> stereotaxic
 *   transformation (see ccfv3ToStereotaxic() and its provenance
 *   warning - constants are currently UNVERIFIED, inherited from
 *   earlier scripts).
 *
 * HEMISPHERE:
 *   Determined from stereotaxic ML sign.
 *   The positive-ML convention is configurable below.
 *
 * QC VERSION NOTE:
 *   Pre-QC automatic counts use the names automaticStarter,
 *   automaticGFPOnly and automaticRV. Final post-QC counts use
 *   finalStarter, finalGFPOnly and finalRV. This avoids duplicate
 *   top-level/local variable declarations in Groovy. Detections
 *   marked "Artifact" during QC are excluded from all of the above.
 *
 * IMPORTANT:
 *   Only detections directly belonging to the CURRENT ROI
 *   are removed or modified. Detections in other ROIs are
 *   preserved.
 *
 * ============================================================
 */

// ============================================================
// IMPORTS
// ============================================================

import net.imglib2.RealPoint
import qupath.ext.biop.abba.AtlasTools
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass

import javax.swing.JButton
import javax.swing.JDialog
import javax.swing.JLabel
import javax.swing.JOptionPane
import javax.swing.JPanel
import javax.swing.WindowConstants
import java.awt.BorderLayout
import java.awt.Dialog
import java.awt.FlowLayout
import java.awt.GridLayout

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
// ML -> HEMISPHERE CONVENTION
// ------------------------------------------------------------
//
// PositiveML_R:
//     positive ML = R
//     negative ML = L
//
// PositiveML_L:
//     positive ML = L
//     negative ML = R
//
// This must remain explicit because ABBA/Allen orientation can
// differ between projects/images.
//

def ML_POSITIVE_HEMISPHERE = "PositiveML_R"

if (!(ML_POSITIVE_HEMISPHERE in ["PositiveML_R", "PositiveML_L"])) {
    throw new Exception(
        "ML_POSITIVE_HEMISPHERE must be exactly " +
        "\"PositiveML_R\" or \"PositiveML_L\". Current value: " +
        ML_POSITIVE_HEMISPHERE
    )
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


// ------------------------------------------------------------
// GFP-derived Starter criterion
// ------------------------------------------------------------
//
// A GFP-derived cell is a Starter if its expanded CELL ROI has
// TRITC/mCherry mean >= this value.
//
// Otherwise it is GFP-only.
//

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
//
// A Pass-2 mCherry detection is classified as Starter when its
// centroid is within this distance of a saved Pass-1 GFP Starter.
//

def STARTER_MATCH_DISTANCE_MICRONS = 5.0


// ============================================================
// ATLAS / STEREOTAXIC TRANSFORM PARAMETERS
// ============================================================
//
// PROVENANCE WARNING (as of v10, revised v11):
//   The constants below (axis offsets, rotation angle, DV squeeze
//   factor) were inherited from earlier scripts:
//     import_coord_Atlas_Stereo_QP7_correct.groovy
//     mRuby2_SypEGFP_two-pass_v4_Atlas_Hemisphere.groovy
//   Their original derivation was NOT documented anywhere in this
//   codebase or its history.
//
//   v11 UPDATE: the rotation and DV squeeze specifically (NOT the
//   X/Y/Z offsets) have since been identified as a mistaken,
//   redundant correction and reset to identity (0 degrees, 1.0x).
//   They were originally entered to compensate for a slight X/Y
//   sectioning angle in some brains - but ABBA's own
//   Atlas Slicing -> X and Y Rotation (deg) step already corrects
//   for that same tilt, in 3D, on the atlas volume itself, BEFORE
//   any CCFv3 coordinate is ever exported. This script's rotation
//   was a separate, simpler 2D correction applied AFTER the fact to
//   coordinates that had already been corrected once by ABBA -
//   effectively double-compensating for the same tilt. See the
//   dedicated SECTIONING-ANGLE CORRECTION section below for where a
//   future tilted sample's correction belongs instead, if needed.
//
//   The axis offsets (OFFSET_X/Y/Z) are UNCHANGED by this and remain
//   of unverified origin - see the ORIGIN OFFSET section below.
//
//   RECOMMENDED FOLLOW-UP (not yet done, needs landmark data):
//     Refit the offsets from >= 3-5 independently-known stereotaxic
//     landmarks (e.g. bregma itself, a midline commissure landmark,
//     one or two hippocampal/cortical landmarks with published
//     Paxinos coordinates) against their corresponding CCFv3
//     coordinates, via least squares, rather than relying on a
//     single example point. Once refit, bump TRANSFORM_VERSION
//     below and update this comment block with the new derivation
//     and its residual error. The standalone companion script
//     CCFv3_Stereotaxic_Transform_Calibration_v1.groovy does this.
//
//   Until that refit happens, PLAUSIBLE_*_RANGE_MM below provide
//   an automatic, per-run sanity check: any computed coordinate
//   outside the expected adult mouse brain envelope is flagged
//   (see "STEREOTAXIC PLAUSIBILITY CHECK" further below). This
//   would catch gross errors (sign flips, axis swaps, unit
//   mistakes) automatically, even without a full recalibration.
//

def TRANSFORM_VERSION =
    "offsets_inherited_unverified__rotation_squeeze_reset_2026-09-30"


// ------------------------------------------------------------
// ORIGIN OFFSET (CCFv3 -> stereotaxic origin alignment)
// ------------------------------------------------------------
//
// Pure translation - aligns the CCFv3 origin with stereotaxic
// bregma. Unrelated to sectioning angle; of unverified origin (see
// provenance warning above). Do not confuse with the sectioning-
// angle correction below - these solve different problems.
//

def CCFV3_TO_STEREO_OFFSET_X = 5.40
def CCFV3_TO_STEREO_OFFSET_Y = 0.44
def CCFV3_TO_STEREO_OFFSET_Z = 5.70


// ------------------------------------------------------------
// SECTIONING-ANGLE CORRECTION (X/Y cutting-angle tilt)
// ------------------------------------------------------------
//
// DEFAULT = IDENTITY (0 degrees, 1.0x) - assumes a negligible X/Y
// cutting angle, matching sections whose tilt was already corrected
// in ABBA (Atlas Slicing -> X and Y Rotation (deg)) at registration
// time, before coordinates were ever exported. For such sections,
// applying anything here would double-compensate for the same tilt.
//
// Only touch these two values for a specific sample with a genuine,
// uncorrected X/Y sectioning angle - and even then, note that this
// is NOT a direct substitute for ABBA's own 3D correction: ABBA
// applies two separate 3D tilt angles (X rotation and Y rotation) to
// the atlas volume itself during registration, whereas this is a
// single 2D in-plane rotation (plus a one-axis DV scale) applied
// afterward, downstream of ABBA, to already-extracted coordinates.
// The two are only approximately equivalent for small angles. The
// correct place to compensate for a real sectioning tilt is
// normally in ABBA itself, not here - these two constants exist so
// that a downstream approximation CAN be entered and clearly
// tracked if it is ever genuinely needed, not as the intended first
// line of correction.
//
// ROTATION_DEGREES: rotates CCF X/Y (-> ML/DV) only; never touches
// Z (-> AP). DV_SQUEEZE: scales the DV axis only, after rotation.
//

def SECTIONING_ANGLE_XY_ROTATION_DEGREES = 0.0
def SECTIONING_ANGLE_DV_SQUEEZE = 1.0

// Kept as aliases so the rest of the script (ccfv3ToStereotaxic(),
// the run-parameters log, etc.) does not need to change: everything
// downstream reads CCFV3_TO_STEREO_ROTATION_DEGREES / _DV_SQUEEZE,
// which now simply mirror the two named constants above.
def CCFV3_TO_STEREO_ROTATION_DEGREES = SECTIONING_ANGLE_XY_ROTATION_DEGREES
def CCFV3_TO_STEREO_DV_SQUEEZE = SECTIONING_ANGLE_DV_SQUEEZE

// Generous adult mouse brain envelope (bregma-relative, mm),
// intended only to catch gross transform errors - NOT a precise
// anatomical boundary, and NOT a substitute for real validation.
def PLAUSIBLE_AP_RANGE_MM = [-8.0, 5.0]
def PLAUSIBLE_ML_RANGE_MM = [-6.5, 6.5]
def PLAUSIBLE_DV_RANGE_MM = [-8.5, 1.0]


// ============================================================
// PER-ROI SIGNAL/BACKGROUND QC PARAMETERS
// ============================================================
//
// These do NOT change any classification or detection - they only
// control automatic warnings, computed per ROI, meant to catch an
// ROI whose staining/signal looks unusual before you trust its
// results. All are heuristic defaults, adjustable per experiment;
// none of them represent a validated biological or statistical
// expectation.
//

// Flag an ROI if it produces fewer detections than this in either
// pass - likely a failed detection, wrong channel, empty tissue, or
// a misplaced ROI.
def MIN_EXPECTED_DETECTIONS_PER_ROI = 5

// A detection's relevant intensity measurement is considered
// "borderline" if it falls within +/- this fraction of the
// threshold that was used to detect/classify it (e.g. 0.20 means
// within 20% of the threshold value in either direction).
def BORDERLINE_THRESHOLD_MARGIN_FRACTION = 0.20

// Flag an ROI if more than this fraction of its detections are
// "borderline" (see above) for a given threshold - suggests that
// threshold is marginal for this image's actual signal level, and
// that small threshold changes would meaningfully change counts.
def BORDERLINE_FRACTION_WARNING_LIMIT = 0.15

// Flag an ROI if the fraction of GFP+ cells classified Starter
// (vs GFP-only) falls outside this range - an all-or-nothing split
// often indicates a threshold miscalibration for that image rather
// than genuine biology.
def EXTREME_STARTER_FRACTION_LOW = 0.02
def EXTREME_STARTER_FRACTION_HIGH = 0.98


// ============================================================
// EXPORT
// ============================================================

def SCRIPT_VERSION = "TVA_Rabies_pipeline_v11"

def EXPORT_RESULTS = true

def EXPORT_OBJECT_FILENAME =
    "TVA_Rabies_Cell_Analysis_v11"

def EXPORT_SUMMARY_FILENAME =
    "TVA_Rabies_ROI_Summary_v11"

def EXPORT_PARAMETERS_FILENAME =
    "TVA_Rabies_Run_Parameters_v11"


// ============================================================
// CURRENT IMAGE / ROI
// ============================================================

def imageData = getCurrentImageData()

if (imageData == null) {
    throw new Exception("No image is currently open.")
}


def analysisParent = getSelectedObject()

if (analysisParent == null || !analysisParent.isAnnotation()) {
    throw new Exception(
        "Please select ONE annotation/ROI before running the script."
    )
}


def ROI_NAME = analysisParent.getName()

if (ROI_NAME == null || ROI_NAME.trim().isEmpty()) {
    ROI_NAME = "Unnamed_ROI"
}


// ------------------------------------------------------------
// Image / slice identification
// ------------------------------------------------------------

def imageName =
    imageData.getServer().getMetadata().getName()

if (imageName == null || imageName.trim().isEmpty()) {
    imageName = "Unnamed_Image"
}

def SLICE_IDENTIFIER = imageName

def safeImageName =
    imageName.replaceAll(
        "[^A-Za-z0-9._-]",
        "_"
    )

def safeROIName =
    ROI_NAME.replaceAll(
        "[^A-Za-z0-9._-]",
        "_"
    )


// ============================================================
// ANIMAL ID (inferred from folder structure)
// ============================================================
//
// There is currently no explicit per-animal identifier anywhere in
// the project or image metadata. As a best current approximation,
// the animal ID is taken from the name of the folder that CONTAINS
// the QuPath project's own base folder - i.e. one level above the
// .qpproj file itself. This assumes each animal has its own QuPath
// project folder, nested inside a parent folder named for that
// animal (e.g. D:\Synaptag\GABA\PV\PV_2F\qupath -> Animal_ID =
// "PV_2F"). If that folder-naming convention isn't followed for a
// given project, this will silently produce the wrong value, so
// it's worth a quick visual check the first time this runs on a
// new project.
//
// NOTE: this does NOT match against any specific folder name - it
// simply takes whatever the parent folder is called, whatever that
// is, and it does not depend on where QuPath itself is installed.
//
// IMPORTANT: this deliberately does NOT read PROJECT_BASE_DIR as a
// plain string. PROJECT_BASE_DIR can hold a literal, unresolved
// token ("{%PROJECT}") rather than an actual path - it appears to
// be a convenience value meant to be resolved BY buildFilePath()
// itself (which is why file export via
// buildFilePath(PROJECT_BASE_DIR, ...) still works fine elsewhere in
// this script), not a string safe to parse directly. Instead,
// getProject().getBaseDirectory() is used, which is QuPath's actual
// Project API and always returns a real, resolved File. This exact
// approach - and this exact bug - was already found and fixed in
// the companion Synaptag pipeline
// (mRuby2_SypEGFP_two-pass_v12_Atlas_Hemisphere onward).
//

def qpProject = getProject()

def qupathBaseDirFile = null

if (qpProject != null) {
    qupathBaseDirFile = qpProject.getBaseDirectory()
}

def animalFolderFile =
    qupathBaseDirFile == null ?
    null :
    qupathBaseDirFile.getParentFile()

def ANIMAL_ID = "Unknown_Animal"

if (animalFolderFile != null) {
    def animalFolderName = animalFolderFile.getName()

    if (animalFolderName != null && !animalFolderName.trim().isEmpty()) {
        ANIMAL_ID = animalFolderName
    }
}

// Always printed as a standing visible confirmation of what
// Animal_ID resolved to on every run, not only when extraction fails
// (matching the companion Synaptag pipeline's v14 behavior).
println "----------------------------------------------"
println "Animal_ID diagnostic"
println "PROJECT_BASE_DIR (raw, unused for this):  " + PROJECT_BASE_DIR
println(
    "QuPath base dir (from Project API):       " +
    (qupathBaseDirFile == null ? "<no project open>" : qupathBaseDirFile.getPath())
)
println(
    "Animal folder (resolved):                 " +
    (animalFolderFile == null ? "<none - no parent found>" : animalFolderFile.getPath())
)
println "Animal_ID:                                 " + ANIMAL_ID
println "----------------------------------------------"

if (ANIMAL_ID == "Unknown_Animal") {
    println(
        "WARNING: Could not infer Animal_ID from the folder " +
        "structure above the QuPath project. See the diagnostic " +
        "lines above for the exact paths that were checked. " +
        "Animal_ID will be recorded as \"Unknown_Animal\"."
    )
}


// ============================================================
// INITIAL INFORMATION
// ============================================================

println ""
println "================================================"
println " TVA/Rabies two-pass cell analysis - v11"
println " QuPath 0.7.0"
println "================================================"
println ""

println "Analysis ROI:    " + ROI_NAME
println "Image / slice:   " + imageName
println "Animal_ID:       " + ANIMAL_ID
println "GFP channel:     " + GFP_CHANNEL
println "mCherry channel: " + MCHERRY_CHANNEL
println "ML convention:   " + ML_POSITIVE_HEMISPHERE
println ""


// ============================================================
// FUNCTION:
// GET ONLY DETECTIONS DIRECTLY UNDER CURRENT ROI
// ============================================================
//
// Deliberately do NOT use global getDetectionObjects() here.
// This allows several analyzed ROIs to coexist in the same image.
//

def getCurrentROIDetections = {
    return analysisParent.getChildObjects().findAll {
        it.isDetection()
    }
}


// ============================================================
// REMOVE EXISTING DETECTIONS FROM CURRENT ROI ONLY
// ============================================================
//
// QuPath 0.7.0: use removeDetections()/removeChildObjects,
// not deprecated clearDetections().
//
// We remove the current ROI's children explicitly so other ROIs
// remain untouched.
//

def existingDetections = getCurrentROIDetections()

if (!existingDetections.isEmpty()) {

    println(
        "Removing existing detections from current ROI: " +
        existingDetections.size()
    )

    // Explicitly remove only the detections belonging to this ROI.
    analysisParent.removeChildObjects(existingDetections)

    fireHierarchyUpdate()

} else {

    println "No existing detections found in current ROI."
}

println ""
println "Detections in other ROIs are preserved."
println ""


// ============================================================
// SIGNAL/BACKGROUND QC RESULTS ACCUMULATOR
// ============================================================
//
// Populated by the Pass-1 and Pass-2 QC blocks below, then written
// out in full alongside the other parameters in the run-parameters
// log file (see EXPORT section near the end of the script).
//

def signalQCResults = [:]


// ============================================================
// CHANNEL NAME RESOLUTION (fail-fast + per-image tolerant match)
// ============================================================
//
// Different images in this dataset can carry slightly different
// channel name metadata even for "the same" staining panel - e.g.
// QuPath/Bio-Formats appends an automatic "(Cn)" disambiguation
// suffix to a channel's name for some images but not others,
// depending on what triggers a naming collision in that specific
// file. A hardcoded exact channel name (e.g. "TRITC") that matches
// most images can silently fail to match one image whose real
// channel name is "TRITC (C3)".
//
// This matters more than it might seem: runPlugin() catches that
// kind of failure INSIDE QuPath's own plugin-running machinery, not
// as a Groovy exception this script can catch - it logs an ERROR
// line and simply does nothing. The result is a silent zero
// detections for that image, with no crash and nothing for this
// script's own error handling to react to.
//
// This resolves each configured channel name against the ACTUAL
// channel names of the CURRENTLY OPEN image before Pass 1 even
// starts: exact match first, falling back to a case-insensitive
// "starts with" match (so "TRITC" resolves to "TRITC (C3)"). If
// neither matches, it fails loudly and immediately, with the full
// list of available channels, rather than letting detection
// silently do nothing partway through the script.
//
// GFP_CHANNEL and MCHERRY_CHANNEL are REASSIGNED to the resolved
// names here. Every later ${GFP_CHANNEL} / ${MCHERRY_CHANNEL}
// interpolation (detection calls, measurement key lookups) is a
// fresh expression evaluated at that point in the script, so all of
// them automatically pick up the corrected value from this point on -
// no other part of the script needs to change.
//

def resolveChannelName = { String configuredName, List actualChannelNames ->

    if (configuredName in actualChannelNames) {
        return configuredName
    }

    def startsWithMatch =
        actualChannelNames.find {
            it.toLowerCase().startsWith(configuredName.toLowerCase())
        }

    if (startsWithMatch != null) {
        println(
            "NOTE: configured channel '" + configuredName +
            "' matched to actual channel '" + startsWithMatch +
            "' for this image (name differs per-image; handled " +
            "automatically)."
        )
        return startsWithMatch
    }

    throw new Exception(
        "Configured channel '" + configuredName + "' was not found " +
        "in this image, and no channel starting with that name " +
        "exists either. Available channels: " + actualChannelNames +
        ". Check GFP_CHANNEL / MCHERRY_CHANNEL at the top of the " +
        "script against this image's actual channel names."
    )
}

def actualChannelNames =
    imageData.getServer().getMetadata().getChannels().collect {
        it.getName()
    }

println "----------------------------------------------"
println "Channel name resolution"
println "----------------------------------------------"
println "Actual channel names in this image: " + actualChannelNames

GFP_CHANNEL = resolveChannelName(GFP_CHANNEL, actualChannelNames)
MCHERRY_CHANNEL = resolveChannelName(MCHERRY_CHANNEL, actualChannelNames)

println "Resolved GFP_CHANNEL:     " + GFP_CHANNEL
println "Resolved MCHERRY_CHANNEL: " + MCHERRY_CHANNEL
println ""


// ============================================================
// PASS 1
// GFP DETECTION
// ============================================================

println "----------------------------------------------"
println "PASS 1: GFP detection"
println "ROI: " + ROI_NAME
println "----------------------------------------------"


// ------------------------------------------------------------
// Explicitly (re-)select the ROI before running detection.
// ------------------------------------------------------------
//
// runPlugin() operates on whatever is CURRENTLY SELECTED in the
// viewer at the moment it runs - it does not take analysisParent
// as an explicit target. analysisParent was only captured once, at
// script start; if anything changes the active selection before
// this point, Pass 1 would silently run against the wrong object
// (or nothing) with no error. This costs nothing and removes that
// entire class of failure.
//

setSelectedObject(analysisParent)

println(
    "Selected object immediately before Pass 1 detection: " +
    (getSelectedObject() == analysisParent ?
        "confirmed = analysisParent (" + ROI_NAME + ")" :
        "MISMATCH - selection is NOT analysisParent! " +
        "Detection below will not run on the intended ROI.")
)
println ""


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

def gfpDetections =
    getCurrentROIDetections()

println(
    "Number of GFP-derived detections in " +
    ROI_NAME + ": " +
    gfpDetections.size()
)

println ""


// ============================================================
// SAVE PASS-1 DATA
// ============================================================
//
// Two groups are retained in memory:
//
//   gfpStarterData
//       Used only for spatial matching during Pass 2.
//
//   gfpOnlyData
//       Reconstructed as true Cell objects after Pass 2.
//
// The original cell and nucleus ROIs are retained.
//

def gfpStarterData = []
def gfpOnlyData = []


def starterCountPass1 = 0
def gfpOnlyCountPass1 = 0


println "GFP detection measurements:"
println "----------------------------------------------"


gfpDetections.eachWithIndex { detection, index ->

    def ml =
        detection.getMeasurementList()


    // --------------------------------------------------------
    // GFP measurements
    // --------------------------------------------------------

    def gfpNucleus =
        ml.get("Nucleus: ${GFP_CHANNEL} mean")

    def gfpCell =
        ml.get("Cell: ${GFP_CHANNEL} mean")


    // --------------------------------------------------------
    // mCherry measurements
    // --------------------------------------------------------

    def mcherryNucleus =
        ml.get("Nucleus: ${MCHERRY_CHANNEL} mean")

    def mcherryCell =
        ml.get("Cell: ${MCHERRY_CHANNEL} mean")


    // --------------------------------------------------------
    // Cell ROI / centroid
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
                "getNucleusROI"
        )) {

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
    // Classification
    // --------------------------------------------------------

    def classification

    if (
        mcherryCell != null &&
        !Double.isNaN(mcherryCell) &&
        mcherryCell >= MCHERRY_STARTER_THRESHOLD
    ) {

        classification = "Starter"
        starterCountPass1++

    } else {

        classification = "GFP-only"
        gfpOnlyCountPass1++
    }


    // --------------------------------------------------------
    // Save reconstruction / matching data
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

        roiName: ROI_NAME,
        imageName: imageName
    ]


    if (classification == "Starter") {

        gfpStarterData << record

    } else {

        gfpOnlyData << record
    }


    // --------------------------------------------------------
    // Temporary Pass-1 classification
    // --------------------------------------------------------

    detection.setClassification(
        classification
    )
}


println ""

println(
    "Starter cells from GFP pass: " +
    starterCountPass1
)

println(
    "GFP-only cells:              " +
    gfpOnlyCountPass1
)

println ""


// ============================================================
// SAVE PASS-1 DATA
// ============================================================

println "----------------------------------------------"
println "Saving Pass-1 GFP-only objects"
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
// PASS-1 SIGNAL/BACKGROUND QC
// ============================================================
//
// See "PER-ROI SIGNAL/BACKGROUND QC PARAMETERS" near the top of
// the script. This does not change any classification - it only
// computes descriptive statistics and prints warnings for values
// this ROI's GFP/mCherry signal takes that are worth a manual look
// before trusting the results.
//

// --------------------------------------------------------
// Small descriptive-statistics helper (mean/median/SD/min/max).
// Returns null if given no valid (non-null, non-NaN) values.
// --------------------------------------------------------

def computeDistributionStats = { List rawValues ->

    def values =
        rawValues.findAll {
            it != null && !Double.isNaN(it)
        }

    if (values.isEmpty()) {
        return null
    }

    def n = values.size()
    def mean = values.sum() / n

    def sorted = values.sort(false)
    def median =
        (n % 2 == 1) ?
            sorted[(int) (n / 2)] :
            (sorted[(int) (n / 2) - 1] + sorted[(int) (n / 2)]) / 2.0

    def variance =
        values.collect { (it - mean) ** 2 }.sum() / n

    def sd = Math.sqrt(variance)

    return [
        n: n,
        mean: mean,
        median: median,
        sd: sd,
        min: sorted.first(),
        max: sorted.last()
    ]
}

def printStats = { String label, Map stats ->

    if (stats == null) {
        println "  " + label + ": no valid measurements"
        return
    }

    println(
        String.format(
            "  %-28s n=%d  mean=%.1f  median=%.1f  sd=%.1f  min=%.1f  max=%.1f",
            label, stats.n, stats.mean, stats.median,
            stats.sd, stats.min, stats.max
        )
    )
}

println "----------------------------------------------"
println "Pass-1 signal/background QC (GFP-derived cells)"
println "----------------------------------------------"

def pass1Records = gfpStarterData + gfpOnlyData

def pass1GFPStats =
    computeDistributionStats(pass1Records.collect { it.gfpCell })

def pass1MCherryStats =
    computeDistributionStats(pass1Records.collect { it.mcherryCell })

printStats("Cell: ${GFP_CHANNEL} mean (Pass 1)", pass1GFPStats)
printStats("Cell: ${MCHERRY_CHANNEL} mean (Pass 1)", pass1MCherryStats)
println ""

def pass1Warnings = []

if (pass1Records.size() < MIN_EXPECTED_DETECTIONS_PER_ROI) {
    pass1Warnings << (
        "Only " + pass1Records.size() + " GFP+ detection(s) found " +
        "in this ROI (expected at least " +
        MIN_EXPECTED_DETECTIONS_PER_ROI +
        "). Check ROI placement, GFP channel selection, and " +
        "GFP_THRESHOLD for this image."
    )
}

def gfpBorderlineMargin =
    GFP_THRESHOLD * BORDERLINE_THRESHOLD_MARGIN_FRACTION

def gfpBorderlineCount =
    pass1Records.count {
        it.gfpCell != null && !Double.isNaN(it.gfpCell) &&
        Math.abs(it.gfpCell - GFP_THRESHOLD) <= gfpBorderlineMargin
    }

def gfpBorderlineFraction =
    pass1Records.isEmpty() ? 0.0 :
        gfpBorderlineCount / (double) pass1Records.size()

if (gfpBorderlineFraction > BORDERLINE_FRACTION_WARNING_LIMIT) {
    pass1Warnings << (
        String.format(
            "%.0f%% of GFP+ detections (%d / %d) have a " +
            "Cell: %s mean within %.0f%% of GFP_THRESHOLD (%.0f). " +
            "Detection in this ROI may be sensitive to small " +
            "changes in GFP_THRESHOLD.",
            gfpBorderlineFraction * 100, gfpBorderlineCount,
            pass1Records.size(), GFP_CHANNEL,
            BORDERLINE_THRESHOLD_MARGIN_FRACTION * 100, GFP_THRESHOLD
        )
    )
}

def mcherryStarterBorderlineMargin =
    MCHERRY_STARTER_THRESHOLD * BORDERLINE_THRESHOLD_MARGIN_FRACTION

def mcherryStarterBorderlineCount =
    pass1Records.count {
        it.mcherryCell != null && !Double.isNaN(it.mcherryCell) &&
        Math.abs(it.mcherryCell - MCHERRY_STARTER_THRESHOLD) <=
            mcherryStarterBorderlineMargin
    }

def mcherryStarterBorderlineFraction =
    pass1Records.isEmpty() ? 0.0 :
        mcherryStarterBorderlineCount / (double) pass1Records.size()

if (mcherryStarterBorderlineFraction > BORDERLINE_FRACTION_WARNING_LIMIT) {
    pass1Warnings << (
        String.format(
            "%.0f%% of GFP+ detections (%d / %d) have a " +
            "Cell: %s mean within %.0f%% of MCHERRY_STARTER_THRESHOLD " +
            "(%.0f). The Starter/GFP-only split for this ROI may be " +
            "sensitive to small changes in that threshold.",
            mcherryStarterBorderlineFraction * 100,
            mcherryStarterBorderlineCount, pass1Records.size(),
            MCHERRY_CHANNEL, BORDERLINE_THRESHOLD_MARGIN_FRACTION * 100,
            MCHERRY_STARTER_THRESHOLD
        )
    )
}

def starterFractionPass1 =
    pass1Records.isEmpty() ? 0.0 :
        starterCountPass1 / (double) pass1Records.size()

if (!pass1Records.isEmpty() &&
    (starterFractionPass1 < EXTREME_STARTER_FRACTION_LOW ||
     starterFractionPass1 > EXTREME_STARTER_FRACTION_HIGH)) {

    pass1Warnings << (
        String.format(
            "Starter fraction among GFP+ cells is %.1f%% (%d / %d) " +
            "- an all-or-nothing split like this often indicates " +
            "MCHERRY_STARTER_THRESHOLD is miscalibrated for this " +
            "image rather than genuine biology.",
            starterFractionPass1 * 100, starterCountPass1,
            pass1Records.size()
        )
    )
}

if (pass1Warnings.isEmpty()) {
    println "No Pass-1 signal/background flags raised."
} else {
    pass1Warnings.each { println "WARNING: " + it }
}

println ""

signalQCResults["Pass1_N_GFP_Positive_Cells"] = pass1Records.size()
signalQCResults["Pass1_GFP_Mean_Of_Cell_Means"] =
    pass1GFPStats?.mean
signalQCResults["Pass1_GFP_Median_Of_Cell_Means"] =
    pass1GFPStats?.median
signalQCResults["Pass1_GFP_SD_Of_Cell_Means"] =
    pass1GFPStats?.sd
signalQCResults["Pass1_mCherry_Mean_Of_Cell_Means"] =
    pass1MCherryStats?.mean
signalQCResults["Pass1_mCherry_Median_Of_Cell_Means"] =
    pass1MCherryStats?.median
signalQCResults["Pass1_mCherry_SD_Of_Cell_Means"] =
    pass1MCherryStats?.sd
signalQCResults["Pass1_Starter_Fraction"] = starterFractionPass1
signalQCResults["Pass1_GFP_Borderline_Fraction"] = gfpBorderlineFraction
signalQCResults["Pass1_mCherryStarter_Borderline_Fraction"] =
    mcherryStarterBorderlineFraction
signalQCResults["Pass1_Warning_Count"] = pass1Warnings.size()
signalQCResults["Pass1_Warnings"] =
    pass1Warnings.isEmpty() ? "none" : pass1Warnings.join(" | ")


// ============================================================
// REMOVE PASS-1 DETECTIONS FROM CURRENT ROI ONLY
// ============================================================
//
// Pass-1 objects are removed from the hierarchy temporarily,
// but their data remain in gfpStarterData / gfpOnlyData.
//
// Other ROIs are untouched.
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


// ------------------------------------------------------------
// Explicitly (re-)select the ROI before running detection.
// ------------------------------------------------------------
//
// See identical comment above Pass 1's runPlugin call. This is
// the more important of the two re-selection points: everything
// between script start and here (Pass 1 detection, classification,
// restoring GFP-only cells, removeChildObjects(), fireHierarchyUpdate())
// is an opportunity for the viewer's active selection to have
// silently changed. If it has, Pass 2 detection would silently
// run against the wrong object (or nothing), reporting 0
// detections with no error - indistinguishable from "there is
// nothing here to detect" unless explicitly checked.
//

setSelectedObject(analysisParent)

println(
    "Selected object immediately before Pass 2 detection: " +
    (getSelectedObject() == analysisParent ?
        "confirmed = analysisParent (" + ROI_NAME + ")" :
        "MISMATCH - selection is NOT analysisParent! " +
        "Detection below will not run on the intended ROI.")
)
println ""


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
        "makeMeasurements":true
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
// PASS-2 SIGNAL/BACKGROUND QC
// ============================================================
//
// Same rationale as the Pass-1 QC block above, applied to the
// independent Pass-2 mCherry detections and MCHERRY_THRESHOLD
// (the Pass-2 detection/segmentation threshold, distinct from
// MCHERRY_STARTER_THRESHOLD used in Pass 1).
//

println "----------------------------------------------"
println "Pass-2 signal/background QC (mCherry-derived cells)"
println "----------------------------------------------"

def pass2MCherryValues =
    mcherryDetections.collect {
        it.getMeasurementList().get("Cell: ${MCHERRY_CHANNEL} mean")
    }

def pass2MCherryStats =
    computeDistributionStats(pass2MCherryValues)

printStats("Cell: ${MCHERRY_CHANNEL} mean (Pass 2)", pass2MCherryStats)
println ""

def pass2Warnings = []

if (mcherryDetections.size() < MIN_EXPECTED_DETECTIONS_PER_ROI) {
    pass2Warnings << (
        "Only " + mcherryDetections.size() + " mCherry+ " +
        "detection(s) found in this ROI (expected at least " +
        MIN_EXPECTED_DETECTIONS_PER_ROI +
        "). Check ROI placement, mCherry channel selection, and " +
        "MCHERRY_THRESHOLD for this image."
    )
}

def mcherryDetectionBorderlineMargin =
    MCHERRY_THRESHOLD * BORDERLINE_THRESHOLD_MARGIN_FRACTION

def mcherryDetectionBorderlineCount =
    pass2MCherryValues.count {
        it != null && !Double.isNaN(it) &&
        Math.abs(it - MCHERRY_THRESHOLD) <= mcherryDetectionBorderlineMargin
    }

def mcherryDetectionBorderlineFraction =
    pass2MCherryValues.isEmpty() ? 0.0 :
        mcherryDetectionBorderlineCount / (double) pass2MCherryValues.size()

if (mcherryDetectionBorderlineFraction > BORDERLINE_FRACTION_WARNING_LIMIT) {
    pass2Warnings << (
        String.format(
            "%.0f%% of Pass-2 mCherry detections (%d / %d) have a " +
            "Cell: %s mean within %.0f%% of MCHERRY_THRESHOLD (%.0f). " +
            "Detection in this ROI may be sensitive to small changes " +
            "in MCHERRY_THRESHOLD.",
            mcherryDetectionBorderlineFraction * 100,
            mcherryDetectionBorderlineCount, pass2MCherryValues.size(),
            MCHERRY_CHANNEL, BORDERLINE_THRESHOLD_MARGIN_FRACTION * 100,
            MCHERRY_THRESHOLD
        )
    )
}

if (pass2Warnings.isEmpty()) {
    println "No Pass-2 signal/background flags raised."
} else {
    pass2Warnings.each { println "WARNING: " + it }
}

println ""

signalQCResults["Pass2_N_mCherry_Positive_Cells"] =
    mcherryDetections.size()
signalQCResults["Pass2_mCherry_Mean_Of_Cell_Means"] =
    pass2MCherryStats?.mean
signalQCResults["Pass2_mCherry_Median_Of_Cell_Means"] =
    pass2MCherryStats?.median
signalQCResults["Pass2_mCherry_SD_Of_Cell_Means"] =
    pass2MCherryStats?.sd
signalQCResults["Pass2_mCherry_Borderline_Fraction"] =
    mcherryDetectionBorderlineFraction
signalQCResults["Pass2_Warning_Count"] = pass2Warnings.size()
signalQCResults["Pass2_Warnings"] =
    pass2Warnings.isEmpty() ? "none" : pass2Warnings.join(" | ")


// ============================================================
// PIXEL CALIBRATION FOR SPATIAL MATCHING
// ============================================================

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
// PASS 2 CLASSIFICATION
// ============================================================
//
// One-to-one nearest-neighbor matching between Pass-2 mCherry
// detections and Pass-1 GFP Starter candidates.
//
// Rationale (v10 change):
//   The previous version matched each Pass-2 detection
//   independently to its single nearest GFP Starter, with no
//   check on whether that Starter had already been claimed by
//   another Pass-2 detection. If a bright/over-segmented soma
//   produced two nearby Pass-2 objects, both could be matched to
//   the same GFP Starter and counted twice.
//
//   Here, ALL valid (detection, starter) pairs within
//   STARTER_MATCH_DISTANCE_MICRONS are collected, sorted by
//   distance, and assigned greedily so that each Pass-2 detection
//   claims at most one GFP Starter, and each GFP Starter is
//   claimed by at most one Pass-2 detection.
//
// Pass-2 detections that match no GFP Starter -> RV+.
//
// GFP Starter candidates that no Pass-2 detection ever claims are
// NOT discarded. Because Pass 1 and Pass 2 use independently
// parameterized segmentations of essentially the same underlying
// mCherry signal, a true double-positive cell can clear the Pass-1
// "Cell: mCherry mean" threshold while failing to form its own
// passing object in the independent Pass-2 segmentation (e.g. a
// diffuse signal that does not meet the Pass-2 area/threshold
// requirements for a discrete detection). Silently dropping such
// cells would understate the Starter count that this whole
// experiment is built around. Unmatched candidates are therefore
// restored further below as "Starter" cells using their own Pass-1
// ROI, with provenance recorded via Starter_Source_Code so they
// remain fully traceable and can be flagged for extra scrutiny
// during manual QC.
//

def starterCount = 0
def rvCount = 0

// Numeric provenance codes (see also restoration step below):
//   0 = not a Starter (RV+ or GFP-only)
//   1 = Starter confirmed by independent Pass-2 mCherry detection
//   2 = Starter restored from Pass-1 GFP-based data only
//       (no matching independent Pass-2 detection found)
def STARTER_SOURCE_NONE = 0.0
def STARTER_SOURCE_PASS2_MATCHED = 1.0
def STARTER_SOURCE_PASS1_UNMATCHED = 2.0

// Give every Pass-1 GFP Starter candidate a stable index so we can
// track which ones get claimed during matching.
gfpStarterData.eachWithIndex { starter, i -> starter.starterIndex = i }

// --------------------------------------------------------
// Collect all candidate (mCherry detection, GFP starter) pairs
// within the allowed matching distance.
// --------------------------------------------------------

def candidatePairs = []

mcherryDetections.eachWithIndex { detection, mIndex ->

    def roi = detection.getROI()
    def mcherryX = roi.getCentroidX()
    def mcherryY = roi.getCentroidY()

    gfpStarterData.each { starter ->

        def dx = mcherryX - starter.centroidX
        def dy = mcherryY - starter.centroidY
        def distancePixels = Math.sqrt(dx * dx + dy * dy)

        if (distancePixels <= starterMatchDistancePixels) {
            candidatePairs << [
                mIndex: mIndex,
                starterIndex: starter.starterIndex,
                distancePixels: distancePixels
            ]
        }
    }
}

// --------------------------------------------------------
// Greedy one-to-one assignment, closest pairs claimed first.
// --------------------------------------------------------

candidatePairs.sort { it.distancePixels }

def assignedStarterForMcherry = [:]   // mIndex -> starterIndex
def claimedStarterIndices = new HashSet()

candidatePairs.each { pair ->

    if (!assignedStarterForMcherry.containsKey(pair.mIndex) &&
        !claimedStarterIndices.contains(pair.starterIndex)) {

        assignedStarterForMcherry[pair.mIndex] = pair.starterIndex
        claimedStarterIndices << pair.starterIndex
    }
}

// --------------------------------------------------------
// Apply classification to each Pass-2 detection.
// --------------------------------------------------------

mcherryDetections.eachWithIndex { detection, index ->

    def ml = detection.getMeasurementList()

    if (assignedStarterForMcherry.containsKey(index)) {

        def starterIndex = assignedStarterForMcherry[index]

        def matchedPair = candidatePairs.find {
            it.mIndex == index && it.starterIndex == starterIndex
        }

        def matchedDistanceMicrons =
            matchedPair.distancePixels * pixelSizeMicrons

        detection.setClassification("Starter")
        starterCount++

        ml.put("Starter_Source_Code", STARTER_SOURCE_PASS2_MATCHED)
        ml.put("Starter_Match_Distance_um", matchedDistanceMicrons)

        println(
            "mCherry detection ${index + 1}: Starter " +
            "(distance = ${matchedDistanceMicrons} µm, " +
            "matched to GFP Starter #${starterIndex + 1})"
        )

    } else {

        detection.setClassification("RV+")
        rvCount++

        ml.put("Starter_Source_Code", STARTER_SOURCE_NONE)

        println("mCherry detection ${index + 1}: RV+")
    }

    // ----------------------------------------------------
    // Numeric analysis flag
    // ----------------------------------------------------
    //
    // DO NOT put analysisParent.getID() into MeasurementList:
    // QuPath's FloatList accepts numeric values, not UUIDs.
    //

    ml.put("TVA_Rabies_Analysis", 1.0)
}


println ""

println(
    "Pass-2 mCherry detections matched to a GFP Starter: " +
    starterCount
)

println(
    "Pass-2 mCherry detections classified RV+:            " +
    rvCount
)

println ""


// ============================================================
// STARTER RECONCILIATION: PASS 1 vs PASS 2
// ============================================================
//
// Identify GFP Starter candidates (from Pass 1) that were never
// claimed by any Pass-2 mCherry detection. These are restored
// later in the script alongside the GFP-only cells.
//

def unmatchedGFPStarters =
    gfpStarterData.findAll {
        !claimedStarterIndices.contains(it.starterIndex)
    }

println "----------------------------------------------"
println "Starter reconciliation (Pass 1 vs Pass 2)"
println "----------------------------------------------"

println(
    "GFP Starter candidates identified in Pass 1:  " +
    gfpStarterData.size()
)

println(
    "Matched to an independent Pass-2 detection:   " +
    claimedStarterIndices.size()
)

println(
    "Unmatched (no corresponding Pass-2 object):   " +
    unmatchedGFPStarters.size()
)

if (!unmatchedGFPStarters.isEmpty()) {

    println ""
    println(
        "WARNING: " + unmatchedGFPStarters.size() +
        " GFP+/mCherry+ candidate(s) from Pass 1 had no matching " +
        "independent detection in Pass 2 within " +
        STARTER_MATCH_DISTANCE_MICRONS + " um."
    )
    println(
        "These are being RESTORED below as Starter cells " +
        "(Starter_Source_Code = " +
        STARTER_SOURCE_PASS1_UNMATCHED +
        ") using their Pass-1 ROI, instead of being silently " +
        "dropped from the final counts."
    )
    println "Please inspect these specifically during manual QC."
}

println ""


// ============================================================
// RESTORE GFP-ONLY CELLS
// ============================================================
//
// IMPORTANT:
// GFP-only cells are recreated AFTER Pass 2.
//
// They are created with createCellObject(), so they are genuine
// QuPath Cell objects rather than generic Detection objects.
//
// They are added directly to analysisParent, ensuring they remain
// under the original ROI in the hierarchy.
//

println "----------------------------------------------"
println "Restoring GFP-only cells"
println "----------------------------------------------"


def restoredGFPOnlyObjects = []

def GFPOnlyPathClass =
    PathClass.fromString("GFP-only")


gfpOnlyData.each { record ->

    def restoredObject =
        PathObjects.createCellObject(
            record.cellROI,
            record.nucleusROI,
            GFPOnlyPathClass
        )


    def ml =
        restoredObject.getMeasurementList()


    // --------------------------------------------------------
    // Restore GFP measurements
    // --------------------------------------------------------

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


    // --------------------------------------------------------
    // Restore mCherry measurements
    // --------------------------------------------------------

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

    ml.put(
        "Starter_Source_Code",
        STARTER_SOURCE_NONE
    )

    ml.put(
        "TVA_Rabies_Analysis",
        1.0
    )


    restoredGFPOnlyObjects <<
        restoredObject
}


// ------------------------------------------------------------
// Add restored cells UNDER THE ORIGINAL ROI
// ------------------------------------------------------------

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


// ============================================================
// RESTORE UNMATCHED PASS-1 GFP STARTER CANDIDATES
// ============================================================
//
// See "STARTER RECONCILIATION" section above for the rationale.
//
// These are candidates that satisfied the Pass-1 GFP+/mCherry+
// criterion but were not claimed by any independent Pass-2
// mCherry detection. Rather than being dropped, they are restored
// as genuine "Starter" Cell objects using their ORIGINAL Pass-1
// cell/nucleus ROI (i.e. the GFP-detection-derived shape, not an
// mCherry-detection-derived shape, since no Pass-2 object exists
// for them).
//
// Starter_Source_Code = 2.0 on every one of these objects makes
// them fully distinguishable downstream (CSV export, QC review)
// from Starter cells confirmed independently in Pass 2
// (Starter_Source_Code = 1.0).
//

println ""
println "----------------------------------------------"
println "Restoring unmatched Pass-1 GFP Starter candidates"
println "----------------------------------------------"


def restoredUnmatchedStarterObjects = []

def StarterPathClass =
    PathClass.fromString("Starter")


unmatchedGFPStarters.each { record ->

    def restoredObject =
        PathObjects.createCellObject(
            record.cellROI,
            record.nucleusROI,
            StarterPathClass
        )


    def ml =
        restoredObject.getMeasurementList()


    // --------------------------------------------------------
    // Restore GFP measurements
    // --------------------------------------------------------

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


    // --------------------------------------------------------
    // Restore mCherry measurements
    // --------------------------------------------------------

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
    // Provenance / numeric metadata
    // --------------------------------------------------------

    ml.put(
        "Starter_Source_Code",
        STARTER_SOURCE_PASS1_UNMATCHED
    )

    ml.put(
        "TVA_Rabies_Analysis",
        1.0
    )


    restoredUnmatchedStarterObjects <<
        restoredObject
}


// ------------------------------------------------------------
// Add restored cells UNDER THE ORIGINAL ROI
// ------------------------------------------------------------

if (!restoredUnmatchedStarterObjects.isEmpty()) {

    analysisParent.addChildObjects(
        restoredUnmatchedStarterObjects
    )

    fireHierarchyUpdate()
}


println(
    "Unmatched GFP Starter candidates restored as Starter cells: " +
    restoredUnmatchedStarterObjects.size()
)

if (!restoredUnmatchedStarterObjects.isEmpty()) {

    println(
        "NOTE: these carry Starter_Source_Code = " +
        STARTER_SOURCE_PASS1_UNMATCHED +
        " and should be given extra attention during manual QC, " +
        "since their shape is Pass-1 (GFP-derived) rather than " +
        "confirmed by an independent Pass-2 mCherry detection."
    )
}

println ""

println ""


// ============================================================
// FINAL DETECTIONS
// ============================================================
//
// This retrieves only the direct detections belonging to the
// current ROI. It includes:
//
//   Cell (Starter)   -> Pass 2
//   Cell (RV+)       -> Pass 2
//   Cell (GFP-only)  -> restored Pass 1
//

def finalDetections =
    getCurrentROIDetections()


// ============================================================
// FINAL CLASS COUNTS
// ============================================================

def automaticStarter =
    finalDetections.count {

        def pc =
            it.getPathClass()

        pc != null &&
        pc.getName() == "Starter"
    }


def automaticGFPOnly =
    finalDetections.count {

        def pc =
            it.getPathClass()

        pc != null &&
        pc.getName() == "GFP-only"
    }


def automaticRV =
    finalDetections.count {

        def pc =
            it.getPathClass()

        pc != null &&
        pc.getName() == "RV+"
    }


// ============================================================
// ATLAS + STEREOTAXIC COORDINATES
// ============================================================
//
// This reproduces the coordinate functionality from:
//
//   import_coord_AtlasXYZ_QP7.groovy
//   import_coord_Atlas_Stereo_QP7_correct.groovy
//
// and the supplied:
//   mRuby2_SypEGFP_two-pass_v4_Atlas_Hemisphere.groovy
//
// Coordinate convention:
//
//   Atlas X -> stereotaxic ML
//   Atlas Y -> stereotaxic DV
//   Atlas Z -> stereotaxic AP
//
// The established transformation is:
//
//   x_stereo = CCF_X - 5.40
//   y_stereo = CCF_Y - 0.44
//   z_stereo = CCF_Z - 5.70
//
//   rotate X/Y by +5 degrees
//   squeeze DV by 0.9434
//
// Final signed coordinates:
//
//   ML = rotated X
//   DV = -rotated Y
//   AP = -Z
//
// This preserves the AP/ML/DV orientation correction previously
// established.
//

println "----------------------------------------------"
println "Calculating atlas / stereotaxic coordinates"
println "----------------------------------------------"


def pixelToAtlasTransform =
    AtlasTools
        .getAtlasToPixelTransform(imageData)
        .inverse()


// ------------------------------------------------------------
// CCFv3 -> stereotaxic transform function
// ------------------------------------------------------------
//
// Isolated into a single function so that, once the constants
// above are validated/refit (see provenance warning where they
// are defined), this is the ONLY place that needs to change.
//
// Returns a Map with keys: ap, ml, dv (all in mm).
//

def ccfv3ToStereotaxic = { double x_ccfv3, double y_ccfv3, double z_ccfv3 ->

    // ----------------------------------------------------
    // CCFv3 -> stereotaxic offset
    // ----------------------------------------------------

    def x_stereo = x_ccfv3 - CCFV3_TO_STEREO_OFFSET_X
    def y_stereo = y_ccfv3 - CCFV3_TO_STEREO_OFFSET_Y
    def z_stereo = z_ccfv3 - CCFV3_TO_STEREO_OFFSET_Z

    // ----------------------------------------------------
    // Rotate X/Y
    // ----------------------------------------------------

    def angleCorrection =
        CCFV3_TO_STEREO_ROTATION_DEGREES / 180.0 * Math.PI

    def rot_x_stereo =
        x_stereo * Math.cos(angleCorrection) -
        y_stereo * Math.sin(angleCorrection)

    def rot_y_stereo =
        x_stereo * Math.sin(angleCorrection) +
        y_stereo * Math.cos(angleCorrection)

    // ----------------------------------------------------
    // Compress DV axis
    // ----------------------------------------------------

    def corrected_y_stereo =
        rot_y_stereo * CCFV3_TO_STEREO_DV_SQUEEZE

    // ----------------------------------------------------
    // Established axis/sign convention
    // ----------------------------------------------------
    //
    // CCF X -> ML, CCF Y -> DV, CCF Z -> AP.
    //
    // The signs are chosen so that the Subiculum example
    // (CCF approximately X=7.4, Y=1.6, Z=8.9) produces
    // approximately AP=-3.2, ML=+2.0, DV=-1.2, matching the
    // expected AP/ML/DV orientation. (These example ML/DV figures
    // assume the default SECTIONING_ANGLE_XY_ROTATION_DEGREES=0 /
    // SECTIONING_ANGLE_DV_SQUEEZE=1.0 - see that section above if
    // you have changed them for a specific tilted sample.) See the
    // provenance warning above the constants: the X/Y/Z offsets
    // remain of unverified origin regardless of this.
    //

    return [
        ap: -z_stereo,
        ml: rot_x_stereo,
        dv: -corrected_y_stereo
    ]
}


// ------------------------------------------------------------
// Stereotaxic plausibility check
// ------------------------------------------------------------
//
// Cheap, automatic per-detection sanity check that does NOT
// require landmark data: flags any AP/ML/DV outside a generous
// adult mouse brain envelope. Meant to catch gross transform bugs
// (sign flips, axis swaps, unit errors) immediately, on every run,
// while the underlying transform constants remain unverified.
//

def isWithinPlausibleRange = { double value, List range ->
    !Double.isNaN(value) && value >= range[0] && value <= range[1]
}

def outOfRangeCoordinateCount = 0


// ------------------------------------------------------------
// Hemisphere function
// ------------------------------------------------------------

def determineHemisphere = { double mlValue ->

    if (Double.isNaN(mlValue)) {
        return "Undefined"
    }

    if (mlValue > 0) {

        return ML_POSITIVE_HEMISPHERE ==
               "PositiveML_R" ?
               "R" :
               "L"
    }

    if (mlValue < 0) {

        return ML_POSITIVE_HEMISPHERE ==
               "PositiveML_R" ?
               "L" :
               "R"
    }

    return "Undefined"
}


// ------------------------------------------------------------
// Coordinate results stored for CSV export
// ------------------------------------------------------------

def coordinateResults = [:]


finalDetections.each { detection ->

    def atlasCoordinates =
        new RealPoint(3)


    // --------------------------------------------------------
    // Detection centroid in pixel coordinates
    // --------------------------------------------------------

    atlasCoordinates.setPosition([

        detection.getROI().getCentroidX(),

        detection.getROI().getCentroidY(),

        0

    ] as double[])


    // --------------------------------------------------------
    // Pixel -> Allen CCFv3
    // --------------------------------------------------------

    pixelToAtlasTransform.apply(
        atlasCoordinates,
        atlasCoordinates
    )


    def x_ccfv3 =
        atlasCoordinates.getDoublePosition(0)

    def y_ccfv3 =
        atlasCoordinates.getDoublePosition(1)

    def z_ccfv3 =
        atlasCoordinates.getDoublePosition(2)


    def ml =
        detection.getMeasurementList()


    // --------------------------------------------------------
    // Atlas measurements
    // --------------------------------------------------------

    ml.put(
        "Atlas_X",
        x_ccfv3
    )

    ml.put(
        "Atlas_Y",
        y_ccfv3
    )

    ml.put(
        "Atlas_Z",
        z_ccfv3
    )


    ml.put(
        "Atlas_Allen_CCFv3_X (mm)",
        x_ccfv3
    )

    ml.put(
        "Atlas_Allen_CCFv3_Y (mm)",
        y_ccfv3
    )

    ml.put(
        "Atlas_Allen_CCFv3_Z (mm)",
        z_ccfv3
    )


    // --------------------------------------------------------
    // CCFv3 -> stereotaxic
    // --------------------------------------------------------
    //
    // See "ATLAS / STEREOTAXIC TRANSFORM PARAMETERS" near the top
    // of the script for the provenance warning covering the
    // constants used inside ccfv3ToStereotaxic().
    //

    def stereo =
        ccfv3ToStereotaxic(
            x_ccfv3,
            y_ccfv3,
            z_ccfv3
        )

    def ap_coordinates_mm = stereo.ap
    def ml_coordinates_mm = stereo.ml
    def dv_coordinates_mm = stereo.dv


    // --------------------------------------------------------
    // Store stereotaxic measurements
    // --------------------------------------------------------

    ml.put(
        "AP (mm)",
        ap_coordinates_mm
    )

    ml.put(
        "ML (mm)",
        ml_coordinates_mm
    )

    ml.put(
        "DV (mm)",
        dv_coordinates_mm
    )


    // --------------------------------------------------------
    // Stereotaxic plausibility check
    // --------------------------------------------------------
    //
    // Does not confirm the transform is CORRECT - only flags
    // results that are anatomically IMPOSSIBLE, which would
    // indicate a sign/axis/unit error. See closure definition
    // above for details.
    //

    def coordinatesPlausible =
        isWithinPlausibleRange(ap_coordinates_mm, PLAUSIBLE_AP_RANGE_MM) &&
        isWithinPlausibleRange(ml_coordinates_mm, PLAUSIBLE_ML_RANGE_MM) &&
        isWithinPlausibleRange(dv_coordinates_mm, PLAUSIBLE_DV_RANGE_MM)

    ml.put(
        "Stereotaxic_Plausible",
        coordinatesPlausible ? 1.0 : 0.0
    )

    if (!coordinatesPlausible) {

        outOfRangeCoordinateCount++

        println(
            "WARNING: implausible stereotaxic coordinate " +
            "(AP=" + ap_coordinates_mm +
            ", ML=" + ml_coordinates_mm +
            ", DV=" + dv_coordinates_mm +
            " mm) for a detection at pixel (" +
            detection.getROI().getCentroidX() + ", " +
            detection.getROI().getCentroidY() +
            "). This may indicate a registration problem for " +
            "this section, or an error in the CCFv3->stereotaxic " +
            "transform. Flagged via Stereotaxic_Plausible = 0."
        )
    }


    // --------------------------------------------------------
    // Hemisphere
    // --------------------------------------------------------

    def hemisphere =
        determineHemisphere(
            ml_coordinates_mm
        )


    // Numeric helper:
    //
    // R = +1
    // L = -1
    // Undefined = 0
    //

    ml.put(
        "Hemisphere_Code",

        hemisphere == "R" ? 1.0 :

        hemisphere == "L" ? -1.0 :

        0.0
    )


    coordinateResults[detection] = [

        atlasX: x_ccfv3,
        atlasY: y_ccfv3,
        atlasZ: z_ccfv3,

        ap: ap_coordinates_mm,
        ml: ml_coordinates_mm,
        dv: dv_coordinates_mm,

        hemisphere: hemisphere,
        stereotaxicPlausible: coordinatesPlausible
    ]
}


fireHierarchyUpdate()


println(
    "Coordinates added to " +
    finalDetections.size() +
    " final detections."
)

println(
    "Transform version: " +
    TRANSFORM_VERSION
)

println ""

println "----------------------------------------------"
println "Stereotaxic plausibility check"
println "----------------------------------------------"

println(
    "Detections with implausible AP/ML/DV: " +
    outOfRangeCoordinateCount +
    " / " +
    finalDetections.size()
)

if (outOfRangeCoordinateCount > 0) {
    println(
        "See WARNING lines above for the affected pixel " +
        "coordinates. Check per-section ABBA registration quality " +
        "for this image before trusting AP/ML/DV values from it, " +
        "and treat the transform constants themselves as still " +
        "unverified (see provenance warning near the top of the " +
        "script)."
    )
}

println ""


// NOTE: v3 preserves the coordinates/hemisphere calculated before QC.
//       The post-QC export filters the existing coordinateResults table
//       instead of recalculating coordinates from the surviving objects.
//
// ============================================================
// INTERACTIVE QUALITY CONTROL
// ============================================================
//
// IMPORTANT:
//   At this point all automatic detection/classification and
//   coordinate calculations have been completed.
//
//   The script now pauses the EXPORT workflow and gives the user
//   an opportunity to inspect the detections in QuPath manually.
//
//   The dialog is modeless so the QuPath viewer remains usable.
//   Delete any erroneous/unwanted detections directly in QuPath.
//   The CSV export is performed only after Continue is pressed,
//   and the export routine re-reads the current ROI detections.
//
//   Therefore, deleted detections are NOT included in the CSV.
//
// ============================================================

println "================================================"
println " INTERACTIVE QUALITY CONTROL"
println "================================================"
println ""
println "Review the detections in QuPath before CSV export."
println "Delete any erroneous or unwanted detections."
println ""


def getCurrentQCDCounts = {
    def detections = getCurrentROIDetections()

    def starter = detections.count {
        def pc = it.getPathClass()
        pc != null && pc.getName() == "Starter"
    }

    def gfpOnly = detections.count {
        def pc = it.getPathClass()
        pc != null && pc.getName() == "GFP-only"
    }

    def rv = detections.count {
        def pc = it.getPathClass()
        pc != null && pc.getName() == "RV+"
    }

    def artifact = detections.count {
        def pc = it.getPathClass()
        pc != null && pc.getName() == "Artifact"
    }

    return [
        starter: starter,
        gfpOnly: gfpOnly,
        rv: rv,
        artifact: artifact,
        total: detections.size()
    ]
}


// ------------------------------------------------------------
// "Mark as Artifact" support
// ------------------------------------------------------------
//
// Rather than deleting an erroneous/unwanted detection outright
// (which leaves no record of it ever having existed), the user can
// select it in the QuPath viewer and click "Mark Selected as
// Artifact" in the QC dialog below. This reclassifies it to
// "Artifact" instead of removing it: the object, its coordinates,
// and its measurements are preserved, but it is excluded from the
// Starter/GFP-only/RV+ counts and from the main results CSV. It is
// instead written to a separate *_Excluded_Artifacts CSV, together
// with the classification it held before being marked, so QC
// decisions remain fully documented rather than silently
// disappearing.
//
// Numeric code stored under "QC_Original_Classification_Code":
//   0 = unknown/none, 1 = Starter, 2 = GFP-only, 3 = RV+
//

// NOTE: setClassification() below is called with the String form
// "Artifact", consistent with how Starter/GFP-only/RV+ are already
// set elsewhere in this script (no separate PathClass object needed).

def classificationNameToCode = { String name ->
    if (name == "Starter") return 1.0
    if (name == "GFP-only") return 2.0
    if (name == "RV+") return 3.0
    return 0.0
}

def classificationCodeToName = { Number code ->
    if (code == null || Double.isNaN(code.doubleValue())) return "Unknown"
    if (code == 1.0) return "Starter"
    if (code == 2.0) return "GFP-only"
    if (code == 3.0) return "RV+"
    return "Unknown"
}

def markSelectedAsArtifact = {

    def currentROIDetectionSet = getCurrentROIDetections() as Set

    def selected =
        getSelectedObjects().findAll {
            it.isDetection() && currentROIDetectionSet.contains(it)
        }

    if (selected.isEmpty()) {
        println(
            "No detections belonging to this ROI are currently " +
            "selected in the viewer - select one or more detections " +
            "first, then click 'Mark Selected as Artifact'."
        )
        return
    }

    def markedCount = 0

    selected.each { detection ->

        def pc = detection.getPathClass()

        if (pc != null && pc.getName() == "Artifact") {
            // Already an artifact - leave its original-classification
            // record untouched.
            return
        }

        def currentName = pc == null ? null : pc.getName()
        def ml = detection.getMeasurementList()

        ml.put(
            "QC_Original_Classification_Code",
            classificationNameToCode(currentName)
        )

        detection.setClassification("Artifact")
        markedCount++
    }

    fireHierarchyUpdate()

    println(
        "Marked " + markedCount + " selected detection(s) as " +
        "Artifact (" + (selected.size() - markedCount) +
        " already were)."
    )
}

def unmarkSelectedArtifact = {

    def currentROIDetectionSet = getCurrentROIDetections() as Set

    def selected =
        getSelectedObjects().findAll {
            it.isDetection() &&
            currentROIDetectionSet.contains(it) &&
            it.getPathClass() != null &&
            it.getPathClass().getName() == "Artifact"
        }

    if (selected.isEmpty()) {
        println(
            "No detections currently classified 'Artifact' are " +
            "selected - select one or more first, then click " +
            "'Unmark Artifact'."
        )
        return
    }

    selected.each { detection ->

        def ml = detection.getMeasurementList()
        def originalCode = ml.get("QC_Original_Classification_Code")
        def originalName = classificationCodeToName(originalCode)

        if (originalName == "Unknown") {
            println(
                "Could not restore original classification for a " +
                "detection (none was recorded) - leaving it as " +
                "Artifact. Reclassify it manually if needed."
            )
        } else {
            detection.setClassification(originalName)
        }
    }

    fireHierarchyUpdate()

    println("Restored original classification for " + selected.size() +
            " detection(s).")
}


// ------------------------------------------------------------
// CSV export is placed in a closure so it can be called AFTER
// the user finishes the manual QC step.
// ------------------------------------------------------------
def exportFinalResults = {

    try {

        println ""
        println "----------------------------------------------"
        println "FINALIZING RESULTS AFTER QUALITY CONTROL"
        println "----------------------------------------------"

        // --------------------------------------------------------
        // Re-read detections AFTER manual QC
        // --------------------------------------------------------
        def rawDetectionsAfterQC = getCurrentROIDetections()

        println "Detections remaining after QC: " +
                rawDetectionsAfterQC.size()
        println ""

        // --------------------------------------------------------
        // Split out Artifact-classified detections
        // --------------------------------------------------------
        //
        // Detections marked "Artifact" during the QC step (see
        // markSelectedAsArtifact() above) are NOT deleted, so they
        // are still present in rawDetectionsAfterQC. They must be
        // excluded from finalDetectionsAfterQC so that every count,
        // density, and CSV row below this point - none of which is
        // otherwise modified - automatically excludes them, exactly
        // as if they had been deleted, while still being documented
        // separately (see EXCLUDED ARTIFACTS EXPORT further below).
        //

        def finalArtifacts =
            rawDetectionsAfterQC.findAll {
                def pc = it.getPathClass()
                pc != null && pc.getName() == "Artifact"
            }

        def finalDetectionsAfterQC =
            rawDetectionsAfterQC.findAll {
                def pc = it.getPathClass()
                pc == null || pc.getName() != "Artifact"
            }

        println "Of which marked as Artifact (excluded, logged separately): " +
                finalArtifacts.size()
        println ""

        // --------------------------------------------------------
        // IMPORTANT: Do NOT recalculate coordinates after QC.
        //
        // Coordinates and hemisphere were calculated once, before
        // the manual QC step, and stored in coordinateResults and
        // on each detection's MeasurementList.  QC only removes
        // detections.  The surviving detections therefore retain
        // exactly the coordinate/hemisphere values generated by
        // the original working coordinate calculation.
        // --------------------------------------------------------

        // Keep only coordinate records belonging to detections that
        // are still present after manual QC.
        def qcDetectionSet = finalDetectionsAfterQC as Set
        def coordinateResultsAfterQC = [:]
        qcDetectionSet.each { detection ->
            if (coordinateResults.containsKey(detection)) {
                coordinateResultsAfterQC[detection] = coordinateResults[detection]
            }
        }

        // Replace the working export table with the QC-filtered table.
        coordinateResults.clear()
        coordinateResults.putAll(coordinateResultsAfterQC)

        fireHierarchyUpdate()

        // --------------------------------------------------------
        // Recalculate classification counts AFTER QC
        // --------------------------------------------------------
        def finalStarter = finalDetectionsAfterQC.count {
            def pc = it.getPathClass()
            pc != null && pc.getName() == "Starter"
        }

        def finalGFPOnly = finalDetectionsAfterQC.count {
            def pc = it.getPathClass()
            pc != null && pc.getName() == "GFP-only"
        }

        def finalRV = finalDetectionsAfterQC.count {
            def pc = it.getPathClass()
            pc != null && pc.getName() == "RV+"
        }

        // --------------------------------------------------------
        // Stereotaxic plausibility count AFTER QC
        // --------------------------------------------------------
        //
        // Recomputed from the QC-filtered coordinateResults (not
        // the pre-QC outOfRangeCoordinateCount) so this reflects
        // only the detections actually present in the final export.
        //

        def outOfRangeCoordinateCountAfterQC =
            coordinateResults.values().count {
                it.stereotaxicPlausible == false
            }

        // --------------------------------------------------------
        // Determine ROI hemisphere AFTER QC
        // --------------------------------------------------------
        def detectedHemispheresAfterQC =
            coordinateResults.values()
                .collect { it.hemisphere }
                .findAll { it == "R" || it == "L" }
                .unique()

        def ROI_HEMISPHERE_AFTER_QC =
            detectedHemispheresAfterQC.size() == 1 ?
                detectedHemispheresAfterQC[0] :
            detectedHemispheresAfterQC.size() > 1 ?
                "Mixed" :
                "Undefined"

        if (ROI_HEMISPHERE_AFTER_QC == "Mixed") {
            println "WARNING: Final QC-approved detections include BOTH R and L coordinates."
            println "Check the selected ROI and ABBA orientation."
            println ""
        }

        // --------------------------------------------------------
        // ROI-level measurements AFTER QC
        // --------------------------------------------------------
        def roiArea = analysisParent.getROI().getArea()
        def roiAreaMm2 = roiArea / 1000000.0

        def totalStarterArea =
            finalDetectionsAfterQC
                .findAll {
                    def pc = it.getPathClass()
                    pc != null && pc.getName() == "Starter"
                }
                .sum { it.getROI().getArea() } ?: 0.0

        def totalGFPOnlyArea =
            finalDetectionsAfterQC
                .findAll {
                    def pc = it.getPathClass()
                    pc != null && pc.getName() == "GFP-only"
                }
                .sum { it.getROI().getArea() } ?: 0.0

        def totalRVArea =
            finalDetectionsAfterQC
                .findAll {
                    def pc = it.getPathClass()
                    pc != null && pc.getName() == "RV+"
                }
                .sum { it.getROI().getArea() } ?: 0.0

        def totalCellArea =
            totalStarterArea +
            totalGFPOnlyArea +
            totalRVArea

        def totalCellDensity =
            roiAreaMm2 > 0 ?
            finalDetectionsAfterQC.size() / roiAreaMm2 :
            Double.NaN

        def starterDensity =
            roiAreaMm2 > 0 ?
            finalStarter / roiAreaMm2 :
            Double.NaN

        def gfpOnlyDensity =
            roiAreaMm2 > 0 ?
            finalGFPOnly / roiAreaMm2 :
            Double.NaN

        def rvDensity =
            roiAreaMm2 > 0 ?
            finalRV / roiAreaMm2 :
            Double.NaN

        // --------------------------------------------------------
        // FINAL RESULTS CONSOLE SUMMARY
        // --------------------------------------------------------
        println "================================================"
        println " FINAL TVA/RABIES RESULTS AFTER QC"
        println "================================================"
        println ""
        println "ROI:                 " + ROI_NAME
        println "Image / slice:       " + imageName
        println "Hemisphere:          " + ROI_HEMISPHERE_AFTER_QC
        println "ML convention:       " + ML_POSITIVE_HEMISPHERE
        println ""
        println "Starter:             " + finalStarter
        println "GFP-only:            " + finalGFPOnly
        println "RV+:                 " + finalRV
        println ""
        println "Total detections:    " + finalDetectionsAfterQC.size()
        println "================================================"
        println ""

        // --------------------------------------------------------
        // CSV EXPORT
        // --------------------------------------------------------
        if (EXPORT_RESULTS) {

            println "Exporting QC-approved measurements..."
            println "PROJECT_BASE_DIR: " + PROJECT_BASE_DIR
            println ""

            // ----------------------------------------------------
            // RUN PARAMETERS LOG
            // ----------------------------------------------------
            //
            // A complete, human-readable record of every tunable
            // parameter and transform constant used to produce this
            // run's output, plus the per-ROI signal/background QC
            // results computed earlier. Written alongside the
            // per-object and per-ROI summary CSVs so any output file
            // can be traced back to exactly what produced it, without
            // relying on remembering which script version was used.
            //

            def runParameters = [:]

            runParameters["Script_Version"] = SCRIPT_VERSION
            runParameters["Run_Timestamp"] =
                new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ss")
                    .format(new Date())
            runParameters["Image_Name"] = imageName
            runParameters["Animal_ID"] = ANIMAL_ID
            runParameters["ROI_Name"] = ROI_NAME
            runParameters["Hemisphere_After_QC"] = ROI_HEMISPHERE_AFTER_QC
            runParameters["ML_Positive_Hemisphere_Convention"] =
                ML_POSITIVE_HEMISPHERE

            runParameters["GFP_Channel"] = GFP_CHANNEL
            runParameters["mCherry_Channel"] = MCHERRY_CHANNEL

            runParameters["GFP_Pixel_Size_um"] = GFP_PIXEL_SIZE
            runParameters["GFP_Background_Radius_um"] = GFP_BACKGROUND_RADIUS
            runParameters["GFP_Median_Radius_um"] = GFP_MEDIAN_RADIUS
            runParameters["GFP_Sigma_um"] = GFP_SIGMA
            runParameters["GFP_Min_Area_um2"] = GFP_MIN_AREA
            runParameters["GFP_Max_Area_um2"] = GFP_MAX_AREA
            runParameters["GFP_Threshold"] = GFP_THRESHOLD
            runParameters["GFP_Cell_Expansion_um"] = GFP_CELL_EXPANSION
            runParameters["mCherry_Starter_Threshold_Pass1"] =
                MCHERRY_STARTER_THRESHOLD

            runParameters["mCherry_Pixel_Size_um"] = MCHERRY_PIXEL_SIZE
            runParameters["mCherry_Background_Radius_um"] =
                MCHERRY_BACKGROUND_RADIUS
            runParameters["mCherry_Median_Radius_um"] = MCHERRY_MEDIAN_RADIUS
            runParameters["mCherry_Sigma_um"] = MCHERRY_SIGMA
            runParameters["mCherry_Min_Area_um2"] = MCHERRY_MIN_AREA
            runParameters["mCherry_Max_Area_um2"] = MCHERRY_MAX_AREA
            runParameters["mCherry_Threshold_Pass2"] = MCHERRY_THRESHOLD
            runParameters["mCherry_Cell_Expansion_um"] = MCHERRY_CELL_EXPANSION

            runParameters["Starter_Match_Distance_um"] =
                STARTER_MATCH_DISTANCE_MICRONS

            runParameters["Transform_Version"] = TRANSFORM_VERSION
            runParameters["CCFv3_To_Stereo_Offset_X"] =
                CCFV3_TO_STEREO_OFFSET_X
            runParameters["CCFv3_To_Stereo_Offset_Y"] =
                CCFV3_TO_STEREO_OFFSET_Y
            runParameters["CCFv3_To_Stereo_Offset_Z"] =
                CCFV3_TO_STEREO_OFFSET_Z
            runParameters["CCFv3_To_Stereo_Rotation_Degrees"] =
                CCFV3_TO_STEREO_ROTATION_DEGREES
            runParameters["CCFv3_To_Stereo_DV_Squeeze"] =
                CCFV3_TO_STEREO_DV_SQUEEZE
            runParameters["Plausible_AP_Range_mm"] =
                PLAUSIBLE_AP_RANGE_MM.join(" to ")
            runParameters["Plausible_ML_Range_mm"] =
                PLAUSIBLE_ML_RANGE_MM.join(" to ")
            runParameters["Plausible_DV_Range_mm"] =
                PLAUSIBLE_DV_RANGE_MM.join(" to ")
            runParameters["Coordinates_Out_Of_Plausible_Range_Count"] =
                outOfRangeCoordinateCountAfterQC

            runParameters["Min_Expected_Detections_Per_ROI"] =
                MIN_EXPECTED_DETECTIONS_PER_ROI
            runParameters["Borderline_Threshold_Margin_Fraction"] =
                BORDERLINE_THRESHOLD_MARGIN_FRACTION
            runParameters["Borderline_Fraction_Warning_Limit"] =
                BORDERLINE_FRACTION_WARNING_LIMIT
            runParameters["Extreme_Starter_Fraction_Low"] =
                EXTREME_STARTER_FRACTION_LOW
            runParameters["Extreme_Starter_Fraction_High"] =
                EXTREME_STARTER_FRACTION_HIGH

            runParameters.putAll(signalQCResults)

            runParameters["Automatic_Starter_Count_Pre_QC"] = automaticStarter
            runParameters["Automatic_GFPOnly_Count_Pre_QC"] = automaticGFPOnly
            runParameters["Automatic_RV_Count_Pre_QC"] = automaticRV
            runParameters["Final_Starter_Count_Post_QC"] = finalStarter
            runParameters["Final_GFPOnly_Count_Post_QC"] = finalGFPOnly
            runParameters["Final_RV_Count_Post_QC"] = finalRV

            def parametersOutputFile =
                buildFilePath(
                    PROJECT_BASE_DIR,
                    EXPORT_PARAMETERS_FILENAME +
                    "_" +
                    safeImageName +
                    "_" +
                    safeROIName +
                    "_" +
                    ROI_HEMISPHERE_AFTER_QC +
                    ".csv"
                )

            def parametersWriter =
                new File(parametersOutputFile).newPrintWriter()

            parametersWriter.println("Parameter,Value")

            runParameters.each { key, value ->
                parametersWriter.println(
                    '"' + key + '","' +
                    (value == null ? "NA" : value.toString()) +
                    '"'
                )
            }

            parametersWriter.close()

            println "Run parameters log written to:"
            println "  " + parametersOutputFile
            println ""

            def objectOutputFile =
                buildFilePath(
                    PROJECT_BASE_DIR,
                    EXPORT_OBJECT_FILENAME +
                    "_" +
                    safeImageName +
                    "_" +
                    safeROIName +
                    "_" +
                    ROI_HEMISPHERE_AFTER_QC +
                    ".csv"
                )

            def objectWriter =
                new File(objectOutputFile).newPrintWriter()

            objectWriter.println(
                "ROI," +
                "Image_Name," +
                "Animal_ID," +
                "Slice," +
                "Object_Type," +
                "Hemisphere," +
                "Classification," +
                "Starter_Source_Code," +
                "Starter_Source_Label," +
                "Starter_Match_Distance_um," +
                "Centroid_X_px," +
                "Centroid_Y_px," +
                "Area_px2," +
                "GFP_Nucleus_Mean," +
                "GFP_Cell_Mean," +
                "mCherry_Nucleus_Mean," +
                "mCherry_Cell_Mean," +
                "Atlas_X," +
                "Atlas_Y," +
                "Atlas_Z," +
                "Atlas_Allen_CCFv3_X (mm)," +
                "Atlas_Allen_CCFv3_Y (mm)," +
                "Atlas_Allen_CCFv3_Z (mm)," +
                "AP (mm)," +
                "ML (mm)," +
                "DV (mm)," +
                "Stereotaxic_Plausible," +
                "Transform_Version"
            )

            finalDetectionsAfterQC.each { detection ->

                def ml = detection.getMeasurementList()
                def roi = detection.getROI()

                def x = roi.getCentroidX()
                def y = roi.getCentroidY()
                def area = roi.getArea()

                def pathClass = detection.getPathClass()
                def classification =
                    pathClass == null ?
                    "" :
                    pathClass.getName()

                def objectType = "Cell"

                def result = coordinateResults[detection]

                def objectHemisphere =
                    result == null ?
                    "Undefined" :
                    result.hemisphere

                def gfpNucleus =
                    ml.get("Nucleus: ${GFP_CHANNEL} mean")

                def gfpCell =
                    ml.get("Cell: ${GFP_CHANNEL} mean")

                def mcherryNucleus =
                    ml.get("Nucleus: ${MCHERRY_CHANNEL} mean")

                def mcherryCell =
                    ml.get("Cell: ${MCHERRY_CHANNEL} mean")

                // --------------------------------------------------------
                // Starter provenance (see PASS 2 CLASSIFICATION /
                // STARTER RECONCILIATION sections):
                //   0 = not a Starter (RV+ or GFP-only)
                //   1 = Starter confirmed by an independent Pass-2
                //       mCherry detection
                //   2 = Starter restored from Pass-1 GFP-based data only
                //       (no matching independent Pass-2 detection found;
                //       flag for extra scrutiny during manual QC)
                // --------------------------------------------------------

                def starterSourceCode =
                    ml.get("Starter_Source_Code")

                def starterSourceLabel

                if (starterSourceCode == null || Double.isNaN(starterSourceCode)) {
                    starterSourceLabel = "NA"
                } else if (starterSourceCode == 1.0) {
                    starterSourceLabel = "Pass2_confirmed"
                } else if (starterSourceCode == 2.0) {
                    starterSourceLabel = "Pass1_only_unmatched"
                } else {
                    starterSourceLabel = "Not_a_starter"
                }

                def starterMatchDistanceUm =
                    ml.get("Starter_Match_Distance_um")

                objectWriter.println(
                    '"' + ROI_NAME + '",' +
                    '"' + imageName + '",' +
                    '"' + ANIMAL_ID + '",' +
                    '"' + SLICE_IDENTIFIER + '",' +
                    '"' + objectType + '",' +
                    '"' + objectHemisphere + '",' +
                    '"' + classification + '",' +
                    starterSourceCode + ',' +
                    '"' + starterSourceLabel + '",' +
                    starterMatchDistanceUm + ',' +
                    x + ',' +
                    y + ',' +
                    area + ',' +
                    gfpNucleus + ',' +
                    gfpCell + ',' +
                    mcherryNucleus + ',' +
                    mcherryCell + ',' +
                    ml.get("Atlas_X") + ',' +
                    ml.get("Atlas_Y") + ',' +
                    ml.get("Atlas_Z") + ',' +
                    ml.get("Atlas_Allen_CCFv3_X (mm)") + ',' +
                    ml.get("Atlas_Allen_CCFv3_Y (mm)") + ',' +
                    ml.get("Atlas_Allen_CCFv3_Z (mm)") + ',' +
                    ml.get("AP (mm)") + ',' +
                    ml.get("ML (mm)") + ',' +
                    ml.get("DV (mm)") + ',' +
                    ml.get("Stereotaxic_Plausible") + ',' +
                    '"' + TRANSFORM_VERSION + '"'
                )
            }

            objectWriter.close()

            // ----------------------------------------------------
            // ROI summary CSV
            // ----------------------------------------------------
            def summaryOutputFile =
                buildFilePath(
                    PROJECT_BASE_DIR,
                    EXPORT_SUMMARY_FILENAME +
                    "_" +
                    safeImageName +
                    "_" +
                    safeROIName +
                    "_" +
                    ROI_HEMISPHERE_AFTER_QC +
                    ".csv"
                )

            def summaryWriter =
                new File(summaryOutputFile).newPrintWriter()

            summaryWriter.println(
                "ROI," +
                "Image_Name," +
                "Animal_ID," +
                "Slice," +
                "Hemisphere," +
                "ML_Convention," +
                "ROI_Area_px2," +
                "ROI_Area_mm2," +
                "Starter_Count," +
                "GFP_only_Count," +
                "RVplus_Count," +
                "Total_Cell_Count," +
                "Starter_Density_per_mm2," +
                "GFP_only_Density_per_mm2," +
                "RVplus_Density_per_mm2," +
                "Total_Cell_Density_per_mm2," +
                "Starter_Total_Area_px2," +
                "GFP_only_Total_Area_px2," +
                "RVplus_Total_Area_px2," +
                "Total_Cell_Area_px2," +
                "Transform_Version," +
                "Coordinates_Out_Of_Plausible_Range_Count"
            )

            summaryWriter.println(
                '"' + ROI_NAME + '",' +
                '"' + imageName + '",' +
                '"' + ANIMAL_ID + '",' +
                '"' + SLICE_IDENTIFIER + '",' +
                '"' + ROI_HEMISPHERE_AFTER_QC + '",' +
                '"' + ML_POSITIVE_HEMISPHERE + '",' +
                roiArea + ',' +
                roiAreaMm2 + ',' +
                finalStarter + ',' +
                finalGFPOnly + ',' +
                finalRV + ',' +
                finalDetectionsAfterQC.size() + ',' +
                starterDensity + ',' +
                gfpOnlyDensity + ',' +
                rvDensity + ',' +
                totalCellDensity + ',' +
                totalStarterArea + ',' +
                totalGFPOnlyArea + ',' +
                totalRVArea + ',' +
                totalCellArea + ',' +
                '"' + TRANSFORM_VERSION + '",' +
                outOfRangeCoordinateCountAfterQC
            )

            summaryWriter.close()

            println "Object results exported to:"
            println objectOutputFile
            println ""
            println "ROI summary exported to:"
            println summaryOutputFile
            println ""

            // ----------------------------------------------------
            // EXCLUDED ARTIFACTS EXPORT
            // ----------------------------------------------------
            //
            // Detections marked "Artifact" during manual QC (see
            // markSelectedAsArtifact() above) are excluded from the
            // main object/summary CSVs above, but are NOT deleted -
            // they remain in the QuPath hierarchy and are documented
            // here, together with the classification they held
            // before being marked, so QC exclusions have a permanent
            // record instead of simply disappearing.
            //

            def artifactsOutputFile =
                buildFilePath(
                    PROJECT_BASE_DIR,
                    EXPORT_OBJECT_FILENAME +
                    "_" +
                    safeImageName +
                    "_" +
                    safeROIName +
                    "_" +
                    ROI_HEMISPHERE_AFTER_QC +
                    "_Excluded_Artifacts.csv"
                )

            def artifactsWriter =
                new File(artifactsOutputFile).newPrintWriter()

            artifactsWriter.println(
                "ROI," +
                "Image_Name," +
                "Animal_ID," +
                "Slice," +
                "Original_Classification," +
                "Hemisphere," +
                "Centroid_X_px," +
                "Centroid_Y_px," +
                "Area_px2," +
                "GFP_Nucleus_Mean," +
                "GFP_Cell_Mean," +
                "mCherry_Nucleus_Mean," +
                "mCherry_Cell_Mean," +
                "AP (mm)," +
                "ML (mm)," +
                "DV (mm)"
            )

            finalArtifacts.each { detection ->

                def ml = detection.getMeasurementList()
                def roi = detection.getROI()

                def originalCode = ml.get("QC_Original_Classification_Code")
                def originalName = classificationCodeToName(originalCode)

                def apMm = ml.get("AP (mm)")
                def mlMm = ml.get("ML (mm)")
                def dvMm = ml.get("DV (mm)")

                def hemisphere = determineHemisphere(mlMm)

                artifactsWriter.println(
                    '"' + ROI_NAME + '",' +
                    '"' + imageName + '",' +
                    '"' + ANIMAL_ID + '",' +
                    '"' + SLICE_IDENTIFIER + '",' +
                    '"' + originalName + '",' +
                    '"' + hemisphere + '",' +
                    roi.getCentroidX() + ',' +
                    roi.getCentroidY() + ',' +
                    roi.getArea() + ',' +
                    ml.get("Nucleus: ${GFP_CHANNEL} mean") + ',' +
                    ml.get("Cell: ${GFP_CHANNEL} mean") + ',' +
                    ml.get("Nucleus: ${MCHERRY_CHANNEL} mean") + ',' +
                    ml.get("Cell: ${MCHERRY_CHANNEL} mean") + ',' +
                    apMm + ',' +
                    mlMm + ',' +
                    dvMm
                )
            }

            artifactsWriter.close()

            println(
                "Excluded artifacts (" + finalArtifacts.size() +
                ") logged to:"
            )
            println artifactsOutputFile
            println ""
        }

        println "================================================"
        println "Analysis complete AFTER manual QC."
        println "ROI: " + ROI_NAME
        println "Image / slice: " + imageName
        println "Hemisphere: " + ROI_HEMISPHERE_AFTER_QC
        println "Starter: " + finalStarter
        println "GFP-only: " + finalGFPOnly
        println "RV+: " + finalRV
        println "Total: " + finalDetectionsAfterQC.size()
        println "================================================"

    } catch (Exception e) {
        println ""
        println "ERROR during finalization / CSV export:"
        e.printStackTrace()

        // A silently-failed export is easy to miss if only printed to
        // the Log panel - show it in a dialog too, so it cannot be
        // mistaken for "nothing happened".
        JOptionPane.showMessageDialog(
            null,
            "CSV export failed:\n\n" +
            e.class.simpleName + ": " + e.getMessage() +
            "\n\nFull details were printed to the QuPath Log panel " +
            "(View > Show log).",
            "TVA/Rabies export error",
            JOptionPane.ERROR_MESSAGE
        )
    }
}


// ============================================================
// BUILD THE QC DIALOG
// ============================================================

def qcDialog = new JDialog()
qcDialog.setTitle("TVA/Rabies v11 - Quality Control")
qcDialog.setModal(false)
qcDialog.setModalityType(Dialog.ModalityType.MODELESS)
qcDialog.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE)


def instructionLabel = new JLabel(
    "<html>Review detections in QuPath. Select unwanted detection(s) and click " +
    "<b>Mark Selected as Artifact</b> (preferred - keeps a record), or delete " +
    "them directly if you prefer. Then click Continue.</html>"
)


def countsLabel = new JLabel()


def updateQCCounts = {
    def counts = getCurrentQCDCounts()

    countsLabel.setText(
        "<html>" +
        "<b>Current ROI: " + ROI_NAME + "</b><br>" +
        "Starter: " + counts.starter +
        " &nbsp;&nbsp; GFP-only: " + counts.gfpOnly +
        " &nbsp;&nbsp; RV+: " + counts.rv +
        " &nbsp;&nbsp; Artifact: " + counts.artifact +
        " &nbsp;&nbsp; <b>Total: " + counts.total + "</b>" +
        "</html>"
    )
}

updateQCCounts()


def refreshButton = new JButton("Refresh counts")
refreshButton.addActionListener {
    updateQCCounts()
}


def markArtifactButton = new JButton("Mark Selected as Artifact")
markArtifactButton.addActionListener {
    markSelectedAsArtifact()
    updateQCCounts()
}


def unmarkArtifactButton = new JButton("Unmark Artifact")
unmarkArtifactButton.addActionListener {
    unmarkSelectedArtifact()
    updateQCCounts()
}


def continueButton = new JButton("Continue to CSV Export")
continueButton.addActionListener {
    continueButton.setEnabled(false)
    refreshButton.setEnabled(false)
    markArtifactButton.setEnabled(false)
    unmarkArtifactButton.setEnabled(false)
    qcDialog.dispose()
    exportFinalResults()
}


def cancelButton = new JButton("Cancel - no CSV export")
cancelButton.addActionListener {
    qcDialog.dispose()
    println ""
    println "CSV export cancelled by user."
    println "Detections remain in the current ROI."
}


def buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER))
buttonPanel.add(refreshButton)
buttonPanel.add(markArtifactButton)
buttonPanel.add(unmarkArtifactButton)
buttonPanel.add(continueButton)
buttonPanel.add(cancelButton)


def mainPanel = new JPanel(new BorderLayout(10, 10))
mainPanel.setBorder(javax.swing.BorderFactory.createEmptyBorder(12, 12, 12, 12))
mainPanel.add(instructionLabel, BorderLayout.NORTH)
mainPanel.add(countsLabel, BorderLayout.CENTER)
mainPanel.add(buttonPanel, BorderLayout.SOUTH)

qcDialog.setContentPane(mainPanel)
qcDialog.setSize(820, 200)
qcDialog.setLocationRelativeTo(null)
qcDialog.setAlwaysOnTop(true)
qcDialog.setVisible(true)

println "QC dialog opened."
println "Select unwanted detections and click 'Mark Selected as Artifact' (preferred), or delete them directly."
println "Click Continue when done. No CSV will be written until Continue is pressed."

// The script intentionally ends here.  The modeless QC dialog remains
// available in QuPath, and Continue triggers exportFinalResults(),
// which re-reads the detections AFTER manual QC.
