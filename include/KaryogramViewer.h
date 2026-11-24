#ifndef KARYOGRAM_VIEWER_H
#define KARYOGRAM_VIEWER_H

#include <string>
#include <vector>
#include <filesystem>
#include <utility>

namespace karyoke
{

struct KaryogramPopupSection                    // One part of a chromosome's popup, e.g. "Homolog 1".
{
    std::string heading;                        // Shown above the counts; may be empty when the chromosome has only one homolog.
    std::vector<std::pair<std::string, int>> counts;   // Mutation type and count, only the types that are present.
};

struct KaryogramHotspot                         // A clickable area of the karyogram (a chromosome label) that opens a popup.
{
    std::string title;                          // Popup title, e.g. "Chromosome 1".
    double x, y, width, height;                 // Area in image pixels, from Karyogram::getLabelRegions().
    std::vector<KaryogramPopupSection> sections;
};

struct KaryogramDetailSegment                   // One stretch of a drawn chromosome piece, at its true size, for the zoomed homolog view.
{
    double startMbp, endMbp;                    // Position along the piece as drawn, measured from its top end.
    double r, g, b;                             // Fill colour (0-1), the colour of the chromosome the material came from.
    bool gap = false;                           // White deletion gap.
    bool reversed = false;                      // Inverted material.
    bool secondHomolog = false;                 // Material from a second homolog: its stretch of outline is dashed.
    double centromereStartMbp = -1.0;           // Centromere inside this stretch, along the piece; -1 when there is none.
    double centromereEndMbp = -1.0;
    double originStartMbp = -1.0;               // Span of these bases on the chromosome they came from (for a white gap: the deleted span); -1 when not known.
    double originEndMbp = -1.0;
    std::string originName;                     // Homolog the material came from, e.g. "Chromosome 5, homolog 2"; empty for white gaps or when not known.
};

struct KaryogramDetailMarker                    // A feature too small to see at most zoom levels, e.g. an indel, marked with a pointer.
{
    double positionMbp;                         // Along the piece as drawn.
    std::string text;                           // Shown next to the pointer, e.g. "Indel: 47 bp deleted".
    double deletedMbp = 0.0;                    // Bases deleted right below the position, drawn as a white gap.
    bool inserted = false;                      // Bases inserted at the position, drawn as a line across the piece.
};

struct KaryogramDetailDamage                    // One SDD damage site with a DSB or SSBs, for the zoomed homolog view (hovering shows its breaks).
{
    double positionMbp;                         // Along the chromatid, from the p-arm end.
    bool doubleStrandBreak = false;             // Drawn as a black line.
    int singleStrandBreaks = 0;                 // Drawn as a white line that gets longer with the count.
};

struct KaryogramDetailPiece                     // One shape drawn in a homolog's slot (SDR: the chromosome, an excised piece, an ecDNA or ring circle; SDD: a chromatid).
{
    double x, y, width, height;                 // Area on the karyogram in image pixels, clickable to open the zoomed view.
    std::string caption;                        // Shown above the piece in the zoomed view, e.g. "Ring"; may be empty.
    double lengthMbp = 0.0;                     // True length of the piece as drawn, white gaps included.
    std::vector<KaryogramDetailSegment> segments;   // Top to bottom.
    std::vector<std::pair<double, double>> chevrons;    // Inverted stretches (start, end in Mbp along the piece): one chevron each.
    std::vector<std::pair<double, double>> movedArrows; // Pieces moved within their own chromosome: from where they were cut out to their middle (Mbp).
    std::vector<int> dashedBoundaries;          // Boundary k lies between segments k-1 and k: drawn dashed (a moved piece's new place).
    std::vector<int> hiddenBoundaries;          // Boundaries drawn without a line (where a moved piece was cut out).
    std::vector<KaryogramDetailMarker> markers;
    std::vector<KaryogramDetailDamage> damages; // SDD only: damage sites on this chromatid.
    bool roundTopCap = true;                    // Telomere end (rounded) or cut end (flat).
    bool roundBottomCap = true;
    bool outline = true;                        // Whole chromosomes have a black outline, excised pieces do not.
    bool circular = false;                      // An ecDNA or ring, drawn as a circle in the zoomed view.
};

struct KaryogramDetailSlot                      // Everything drawn for one homolog, opened together when any of its pieces is clicked.
{
    std::string title;                          // E.g. "Chromosome 5, homolog 1".
    std::vector<KaryogramDetailPiece> pieces;
};

struct KaryogramImage                           // One karyogram to show in the web viewer.
{
    std::string label;                          // Shown in the cell list, e.g. "Cell 0".
    std::string filename;                       // Name used for downloads from the viewer, e.g. "cell0_karyogram.png" (the extension is swapped to .svg for SVG downloads).
    std::vector<unsigned char> svgBytes;        // Vector image (KaryogramFormat::SVG): shown in the viewer, stays sharp at any zoom level, offered as SVG download.
    std::vector<unsigned char> pngBytes;        // Raster image (KaryogramFormat::PNG): offered as PNG download, and shown only when svgBytes is empty.
                                                // Either one may be left empty; the viewer then only offers the format it has.
    std::vector<KaryogramHotspot> hotspots;     // Clickable chromosome labels with mutation or damage summaries (may be empty).
    std::string popupNote;                      // Small print at the bottom of each popup (may be empty).
    std::string popupEmptyText = "No mutations detected";   // Shown for a homolog with nothing to list, e.g. "No damage detected" for SDD.
    std::vector<KaryogramDetailSlot> detailSlots;   // Clickable homologs with their pieces at true scale, from Karyogram::getDetailSlots() (may be empty).
};


// Writes a single self-contained HTML file with every karyogram embedded in it (no server or extra files needed).
bool writeKaryogramViewer(
    const std::vector<KaryogramImage>& images,  // Karyograms to show, in order.
    const std::string& title,                   // Page title, e.g. the input file name.
    const std::filesystem::path& htmlPath       // Where to write the HTML file.
);

// Writes several files (e.g. every karyogram of an input file) into one .zip archive, stored without compression since PNG and
// SVG output is already compact enough. Entries are {file name inside the zip, file contents}. Returns false if it could not write.
bool writeZipArchive(
    const std::vector<std::pair<std::string, std::vector<unsigned char>>>& entries,
    const std::filesystem::path& zipPath
);

// Returns false when no browser can be shown on this machine, e.g. over SSH, on an HPC cluster or in a container,
// or in WSL with Windows interop switched off. Used to decide where to write the viewer and whether to open it.
bool canOpenBrowser();

// Opens a file in the system's default web browser (xdg-open on Linux, open on macOS, start on Windows, the Windows
// browser from WSL). Returns false if it could not, in which case the caller should tell the user where the file is.
bool openInBrowser(const std::filesystem::path& htmlPath);

} // namespace karyoke

#endif
