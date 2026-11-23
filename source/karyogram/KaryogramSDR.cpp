#include <cairo/cairo.h>
#include <cmath>
#include <limits>
#include <algorithm>
#include <map>
#include <string>
#include <set>
#include <iostream>

#include "Karyogram.h"
#include "SDRutilities.h"





namespace karyoke
{

// Gap between a chromosome and a column of circular pieces (rings, ecDNA) beside it. Wider than the 6 px stackGap between linear pieces,
// because a round edge only comes closest to the chromosome at one point and reads as touching it at 6 px.
static constexpr double circleColumnGap = 12.0;

// A linear record that mixes home-strand fragments with foreign (exchanged) ones can still hold lost material: two ADJACENT fragments of the SAME strand (home or foreign) that do not touch
// leave a gap (an independent deletion elsewhere on the strand, e.g. next to a deletion-translocation or deletion-insertion). Returns the record's drawn length in Mbp with
// those gaps included, and sets hasGap when at least one exists. Mirrors the gap rule of buildMixedStrandSegmentsWithGaps().
static double mixedRecordLengthWithHomeGaps(const SDRdataRecord& record, int homeOldStrandID, bool& hasGap)
{
    const double tolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)
    double totalLengthMbp = 0.0;
    hasGap = false;

    for (std::size_t i = 0; i < record.fragments.size(); ++i)
    {
        const SDRfragment& fragment = record.fragments[i];
        totalLengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);

        if (i == 0 || fragment.oldStrandID != record.fragments[i - 1].oldStrandID)
        {
            continue;
        }

        const SDRfragment& previous = record.fragments[i - 1];

        // Foreign material is only trusted to be in ascending order, a shuffled (chromothripsis-style) foreign pair says nothing about lost material.
        if (fragment.oldStrandID != homeOldStrandID && (previous.oldStartPosition > previous.oldEndPosition || fragment.oldStartPosition > fragment.oldEndPosition || previous.oldEndPosition > fragment.oldStartPosition))
        {
            continue;
        }
        const double prevLo = std::min(previous.oldStartPosition, previous.oldEndPosition);
        const double prevHi = std::max(previous.oldStartPosition, previous.oldEndPosition);
        const double curLo = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
        const double curHi = std::max(fragment.oldStartPosition, fragment.oldEndPosition);

        if (prevHi <= curLo && !approxEqual(prevHi, curLo, tolerance))
        {
            totalLengthMbp += curLo - prevHi;
            hasGap = true;
        }
        else if (curHi <= prevLo && !approxEqual(curHi, prevLo, tolerance))
        {
            totalLengthMbp += prevLo - curHi;
            hasGap = true;
        }
    }

    return totalLengthMbp;
}

static bool hasPiecesMovedWithinStrand(const std::vector<SDRfragment>& fragments);

// Did pieces of the home strand of a mixed record move within the chromosome (see hasPiecesMovedWithinStrand())? Foreign pieces are left out.
// Only a record that carries the home strand's centromere (the body of the chromosome) qualifies: its home pieces are laid out against the whole
// chromosome, with the stretches it lacks drawn as deleted. An acentric piece of the chromosome (e.g. the partner of a dicentric) holds only part of
// it, and the rest usually sits in another record rather than being deleted, so it keeps the plain mixed layout.
static bool homePiecesMoved(const SDRdataRecord& record, int homeOldStrandID)
{
    std::vector<SDRfragment> home;
    bool homeCentromere = false;
    for (const SDRfragment& fragment : record.fragments)
    {
        if (fragment.oldStrandID == homeOldStrandID)
        {
            home.push_back(fragment);
            homeCentromere = homeCentromere || fragment.hasCentromere;
        }
    }
    return homeCentromere && home.size() != record.fragments.size() && home.size() >= 2 && hasPiecesMovedWithinStrand(home);
}

static bool isSingleMixedRecordWithHomeGap(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID)
{
    if (records.size() != 1 || !records[0]->linear)
    {
        return false;
    }

    bool hasForeign = false;
    for (const SDRfragment& fragment : records[0]->fragments)
    {
        hasForeign = hasForeign || fragment.oldStrandID != homeOldStrandID;
    }

    bool hasGap = false;
    mixedRecordLengthWithHomeGaps(*records[0], homeOldStrandID, hasGap);

    return hasForeign && hasGap;
}



bool Karyogram::generateSDRkaryogram(					// Writes the SDR karyogram to a PNG (default) or SVG file.
    const SDRmasterHeader& masterHeader,
    const SDRsubHeader& subHeader,
    bool humanGenome,
    const std::string& outputFilename,
    KaryogramFormat format)
{
    return renderSDRkaryogram(masterHeader, subHeader, humanGenome, outputFilename, nullptr, format);
}


bool Karyogram::generateSDRkaryogram(					// Returns the SDR karyogram as PNG or SVG bytes in memory.
    const SDRmasterHeader& masterHeader,
    const SDRsubHeader& subHeader,
    bool humanGenome,
    std::vector<unsigned char>& imageBytes,
    KaryogramFormat format)
{
    return renderSDRkaryogram(masterHeader, subHeader, humanGenome, std::string(), &imageBytes, format);
}


bool Karyogram::renderSDRkaryogram(					// Function to draw a karyogram of the rearranged chromosome segments on the karyogram
    const SDRmasterHeader& masterHeader,				// Access SDR master header for intact chromosome sizes
    const SDRsubHeader& rawSubHeader, 					// Access SDR subheader for cell-by-cell damage/misrepair data
    bool humanGenome,							// Check if human genome specified
    const std::string& outputFilename,					// Output file name of the karyogram, concatenated with the start of the SDR filename (when imageBytes is null).
    std::vector<unsigned char>* imageBytes,				// If non-null, receives the image in memory instead.
    KaryogramFormat format)						// PNG or SVG.
{

    labelRegions.clear();						// Label positions are recorded fresh for every karyogram
    detailSlots.clear();						// So are the homolog slots of the zoomed view
    detailStrandNames.clear();
    currentDetailSlot = -1;
    detailIndels.clear();
    detailMasterHeader = &masterHeader;
    detailHumanGenome = humanGenome;

    // Check for absent data in SDR headers/subHeaders
    if (rawSubHeader.dataRecords.empty())
    {
        std::cerr << "ERROR: No SDR data records to draw for cell " << rawSubHeader.cellID << ".\n";
        return false;
    }

    // Chromosome sizes must be passed in the masterHeader, otherwise Karyogram cannot be drawn and mutations cannot be checked.
    const std::vector<double>& sizes = masterHeader.intactChromosomeSizes;
    if (sizes.empty())
    {
        std::cerr << "ERROR: SDR master header has no chromosome sizes.\n";
        return false;
    }

    // Indels (< 100 bp) are only counted, never drawn: the chromosomes are drawn from the cell with its indels taken out (see removeIndels), so an
    // inserted base or a 5 bp gap neither adds a sliver of another chromosome's colour nor splits a piece. The summary box still counts them.
    // The indels taken out are kept, to be marked in the zoomed homolog view.
    const SDRsubHeader subHeader = removeIndels(rawSubHeader, static_cast<int>(sizes[0]), &detailIndels);


    // All data records grouped by their determined home strand (the centromere-bearing fragment when unambiguous, 
    // overridden for dicentric/acentric pairs, falling back to the first listed fragment otherwise)
    std::map<int, std::vector<const SDRdataRecord*>> recordsByOriginalStrand;
    // Use to merge aberrant fragments onto a single chromosome in karyogram
    std::vector<SDRdataRecord> mergedRecordStorage;
    mergedRecordStorage.reserve(2 * subHeader.dataRecords.size()); 			// Safe upper bound (clustered bars plus the acentric slots after X and Y).

    // Store maps of chromosome IDs for dicentric and acentric chromosomes to be drawn in their corresponding remaining home slots.
    std::map<int, int> dicentricAcentricOverrides = findDicentricAcentricHomeOverrides(subHeader, static_cast<int>(sizes[0]), masterHeader);

    // Two-chromosome chromothripsis: its multi-fragment records interleave both strands, so they are drawn as concatenated bars, one per strand slot,
    // whether or not the excised single-fragment pieces exist in the SDR file. The dicentric/acentric pairing above needs those excised pieces to tile the
    // strands, so when they are absent it finds nothing, and every record would fall back to the first fragment's strand and pile into ONE slot,
    // leaving the other strand drawn untouched. Give each multi-fragment record of the cluster its own strand slot here instead.
    // Complex genomic rearrangements are found on the same kind of cluster, so they are drawn the same way.
    std::vector<SDRchromothripsisEvent> shatteredClusters = detectChromothripsis(subHeader, static_cast<int>(sizes[0]));
    const std::vector<SDRcomplexRearrangementEvent> complexRearrangements = detectComplexGenomicRearrangements(subHeader, static_cast<int>(sizes[0]));
    shatteredClusters.insert(shatteredClusters.end(), complexRearrangements.begin(), complexRearrangements.end());

    for (const SDRchromothripsisEvent& event : shatteredClusters)
    {
        if (event.involvedStrandIDs.size() != 2)
        {
            continue;
        }

        std::set<int> takenSlots;
        for (const SDRdataRecord& record : subHeader.dataRecords)
        {
            if (record.fragments.size() < 2 || std::find(event.contributingNewStrandIDs.begin(), event.contributingNewStrandIDs.end(), record.newStrandID) == event.contributingNewStrandIDs.end())
            {
                continue;
            }

            const auto existing = dicentricAcentricOverrides.find(record.newStrandID);
            int slot = (existing != dicentricAcentricOverrides.end()) ? existing->second : determineHomeStrandID(record);

            if (takenSlots.count(slot))
            {
                for (int strandID : event.involvedStrandIDs)
                {
                    if (!takenSlots.count(strandID))
                    {
                        slot = strandID;
                        break;
                    }
                }
            }

            takenSlots.insert(slot);
            dicentricAcentricOverrides[record.newStrandID] = slot;
        }
    }


    // Chromoplexy-specific shape: an original strand fully consumed by fusing onto two OTHER strands' centromeres (one piece each) gets hidden 
    // from its own slot entirely, and its two partners' own leftover acentric pieces are drawn together in one new slot after X and Y.
    std::vector<SDRchromoplexyAcentricPlacement> chromoplexyPlacements = findChromoplexyAcentricRecombinations(subHeader, static_cast<int>(sizes[0]), masterHeader);

    // Any other linear record without a centromere that joins pieces of two or more chromosomes (an acentric fragment left by a multi-chromosome
    // rearrangement) is drawn in the slot of one of those chromosomes that would otherwise be empty, its first one in the record. When every one of its chromosomes already has something
    // to draw, it gets its own slot after X and Y, like the chromoplexy leftovers. Acentric partners of a dicentric keep their slot.
    {
        std::set<int> placed;
        for (const SDRchromoplexyAcentricPlacement& placement : chromoplexyPlacements)
        {
            placed.insert(placement.leftoverNewStrandIDs.begin(), placement.leftoverNewStrandIDs.end());
        }

        std::vector<const SDRdataRecord*> candidates;
        std::map<int, int> recordsPerHome;						// Slot strand -> records drawn there (candidates left out)

        for (const SDRdataRecord& record : subHeader.dataRecords)
        {
            std::set<int> strands;
            bool centric = false;
            for (const SDRfragment& fragment : record.fragments)
            {
                strands.insert(fragment.oldStrandID);
                centric = centric || fragment.hasCentromere;
            }

            if (placed.count(record.newStrandID) || record.fragments.empty())
            {
                continue;
            }

            if (record.linear && !centric && strands.size() >= 2 && !dicentricAcentricOverrides.count(record.newStrandID))
            {
                candidates.push_back(&record);
                continue;
            }

            const auto overrideIt = dicentricAcentricOverrides.find(record.newStrandID);
            recordsPerHome[(overrideIt != dicentricAcentricOverrides.end()) ? overrideIt->second : determineHomeStrandID(record)]++;
        }

        for (const SDRdataRecord* record : candidates)
        {
            int emptySlot = -1;
            for (const SDRfragment& fragment : record->fragments)
            {
                if (recordsPerHome[fragment.oldStrandID] == 0)
                {
                    emptySlot = fragment.oldStrandID;
                    break;
                }
            }

            if (emptySlot != -1)
            {
                dicentricAcentricOverrides[record->newStrandID] = emptySlot;
                recordsPerHome[emptySlot]++;
            }
            else
            {
                chromoplexyPlacements.push_back({-1, {record->newStrandID}});
            }
        }
    }

    std::set<int> chromoplexyConsumedStrands;
    std::set<int> chromoplexyRecombinedNewStrandIDs; 					// newStrandIDs now drawn ONLY in the dedicated chromoplexy slot, not in their own chromosome's slot.
    for (const SDRchromoplexyAcentricPlacement& placement : chromoplexyPlacements)
    {
        chromoplexyConsumedStrands.insert(placement.consumedStrandID);
        for (int leftoverNewStrandID : placement.leftoverNewStrandIDs)
        {
            chromoplexyRecombinedNewStrandIDs.insert(leftoverNewStrandID);
        }
    }


    // ecDNA assembled from excised fragments of TWO OR MORE different original strands is acentric and doesn't "belong" to any one chromosome. 
    // Like the chromoplexy leftover pieces above, it gets its own dedicated slot instead of being grouped under whichever
    // strand determineHomeStrandID() happens to fall back to. A same-origin ecDNA (everyfragment from the same strand, however many gaps) is 
    // left alone here: isECDNAshape() already draws those correctly next to their home chromosome.
    const std::vector<SDRecDNAevent> ecDNAevents = detectECDNA(subHeader, static_cast<int>(sizes[0]), masterHeader, false);	// false: drawing wants every circular piece, including those scored under a chromothripsis

    std::set<int> multiOriginECDNAnewStrandIDs; 					// excisedStrandID of every multi-origin ecDNA event.
    for (const SDRecDNAevent& event : ecDNAevents)
    {
    	const std::set<int> distinctOrigins(event.oldStrandIDs.begin(), event.oldStrandIDs.end());

        if (distinctOrigins.size() >= 2)
        {
            multiOriginECDNAnewStrandIDs.insert(event.excisedStrandID);
        }
    }


    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (record.fragments.empty())
        {
            continue;
        }
	// This record's material is drawn only in the dedicated chromoplexy acentric slot below, or the dedicated multi-foreign fragment ecDNA 
	// acentric slot, leave it out of its own chromosome's normal slot grouping so it isn't drawn twice.
	else if (chromoplexyRecombinedNewStrandIDs.count(record.newStrandID) || multiOriginECDNAnewStrandIDs.count(record.newStrandID))
	{
	    continue;
	}
        else
        {
	    // A dicentric/acentric pair needs an explicit override, since the acentric side's correct home strand depends on the
            // PAIRED dicentric record's second fragment (information determineHomeStrandID() can't see on its own).
            const auto overrideIt = dicentricAcentricOverrides.find(record.newStrandID);
            const int homeStrandID = (overrideIt != dicentricAcentricOverrides.end()) ? overrideIt->second : determineHomeStrandID(record);

	    // Grouping drawings to home strand which is whichever strand has the centromere, not the first fragment in the list, unless no centromere in the record, 
	    // then falls back to first fragment in the list.
            recordsByOriginalStrand[homeStrandID].push_back(&record);
        }
    }



    // Determining the chromosome layout, based on how chromosome sizes were passed
    const int chromosomeCount = static_cast<int>(sizes[0]);				// First entry in SDR intact chromosome sizes header is the number of chromosomes listed

    const ChromosomeLayout chromosomeLayout = determineChromosomeLayout(sizes);		// The order of how the chromosome sizes are passed alters the karyogram drawing code.
    const bool hasHomologs = chromosomeLayout != ChromosomeLayout::NON_HOMOLOGOUS;	// Default chromosome layout is non-homologous: [count, 1,2,3,...,N]
    const int homologousPairs = hasHomologs ? (chromosomeCount - 2) / 2 : 0;		// Number of homologous pairs excludes the X and Y sex chromosomes and then divided by two for pairs.
    const int drawableGroups = hasHomologs ? homologousPairs : chromosomeCount - 2;	// Determines the number of chromosome groups to draw in the karyogram, the X and Y chromosomes are excluded and will be drawn at the end.


    // Storage for synthesized baseline records, for original strand that have NO data records at all in the file.
    std::vector<SDRdataRecord> intactRecordStorage;
    intactRecordStorage.reserve(static_cast<std::size_t>(chromosomeCount));


    // Scale: the largest intact chromosome size in the master header is drawn maxRenderHeight tall, in every cell of the file.
    // Why: the scale used to come from the longest record in each cell. A cell whose records were all shorter than the largest
    // chromosome (e.g. a single deletion, with the intact strands not listed as records) got a larger scale, so its intact chromosomes
    // grew into the next row and the legend; a cell with a long derivative got a smaller one. Chromosome sizes now look the same in
    // every karyogram, and a derivative longer than the largest chromosome makes its own row taller instead (see rowTop below).
    double maxLengthMbp = 0.0;

    for (std::size_t i = 1; i < sizes.size(); ++i)
    {
        maxLengthMbp = std::max(maxLengthMbp, sizes[i]);
    }

    // Total length of a record's fragments, in Mbp.
    const auto recordLengthMbp = [](const SDRdataRecord& record)
    {
        double lengthMbp = 0.0;

        for (const SDRfragment& fragment : record.fragments)
        {
            lengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
        }

        return lengthMbp;
    };

    if (maxLengthMbp <= 0.0)							// No usable sizes in the header: fall back to the longest record.
    {
        for (const SDRdataRecord& record : subHeader.dataRecords)
        {
            maxLengthMbp = std::max(maxLengthMbp, recordLengthMbp(record));
        }
    }

    if (maxLengthMbp <= 0.0)
    {
        std::cerr << "ERROR: SDR data records have no positive-length fragments.\n";
        return false;
    }


    // --------------------------------------------------
    // Karyogram dimensions
    // --------------------------------------------------

    // Fixed layout constants for the whole karyogram image. imgWidth/columns/rowHeight/startY control
    // the grid every slot is positioned within (changing them requires re-checking the row/column math
    // throughout this function); maxRenderHeight/chromosomeWidth/homologGap are purely visual sizing and
    // safe to adjust independently.
    const int imgWidth = 1000;
    const int columns = 4;
    const double colWidth = static_cast<double>(imgWidth) / columns;
    const double maxRenderHeight = 180.0;					// Height of the largest intact chromosome, safe to tune.
    detailMbpPerPixel = maxLengthMbp / maxRenderHeight;
    const double rowHeight = 250.0;						// Default row height; a row holding a longer derivative grows (see rowTop below).
    const double rowSpacing = rowHeight - maxRenderHeight;			// Room under the tallest bar of a row for its label and the gap to the next row.
    const double startY = 288.0;						// Below the summary box (5 rows of counts, bottom at y = 218).
    const int rows = (drawableGroups + columns - 1) / columns;

    // Check if X and Y can fit in last row, if not, put in their own row
    const int sexChromosomeRow = (drawableGroups > 0) ? (drawableGroups - 1) / columns : 0;
    const int lastRowUsedCols = drawableGroups - sexChromosomeRow * columns;
    const bool sexChromosomesFitInLastRow = (columns - lastRowUsedCols) >= 2;

    // Chromoplexy acentric-recombination slots and multi-origin ecDNA slots share ONE row sequence after X and Y. ecDNA slots continue filling any 
    // columns the chromoplexy slots left open in their last row before wrapping to a new row.
    const int chromoplexySlotCount = static_cast<int>(chromoplexyPlacements.size());
    const int multiOriginECDNAslotCount = static_cast<int>(multiOriginECDNAnewStrandIDs.size());
    // Combine acentric rows for both ecDNA and chromoplexy to pack into the same rows both mutations if there is space, otherwise 
    // move on to new row.
    const int combinedAcentricSlotCount = chromoplexySlotCount + multiOriginECDNAslotCount;
    const int combinedAcentricRows = (combinedAcentricSlotCount + columns - 1) / columns;
    const int autosomeAndSexRows = sexChromosomesFitInLastRow ? rows : rows + 1;
    const int acentricsStartRow = autosomeAndSexRows;
    const int sexRow = sexChromosomesFitInLastRow ? sexChromosomeRow : rows;
    const int totalRows = autosomeAndSexRows + combinedAcentricRows;

    // Row heights. Every row is rowHeight tall, except a row holding a derivative longer than the largest chromosome, which grows just enough
    // to fit it, so it doesn't run into the row below or the legend. Only that row grows; the others keep the default height.
    std::vector<double> rowLongestMbp(static_cast<std::size_t>(totalRows), 0.0);	// Longest strand drawn in each row

    const auto noteSlotStrand = [&](int row, int oldStrandID)
    {
        const auto it = recordsByOriginalStrand.find(oldStrandID);

        if (row < 0 || row >= totalRows || it == recordsByOriginalStrand.end() || chromoplexyConsumedStrands.count(oldStrandID))
        {
            return;
        }

        for (const SDRdataRecord* record : it->second)
        {
            rowLongestMbp[row] = std::max(rowLongestMbp[row], recordLengthMbp(*record));
        }
    };

    for (int i = 0; i < drawableGroups; ++i)					// Same slot-to-strand mapping as the drawing loop below
    {
        if (chromosomeLayout == ChromosomeLayout::SPLIT_HOMOLOGS)
        {
            noteSlotStrand(i / columns, i + 1);
            noteSlotStrand(i / columns, i + 1 + homologousPairs);
        }
        else if (chromosomeLayout == ChromosomeLayout::ADJACENT_HOMOLOGS)
        {
            noteSlotStrand(i / columns, 2 * i + 1);
            noteSlotStrand(i / columns, 2 * i + 2);
        }
        else
        {
            noteSlotStrand(i / columns, i + 1);
        }
    }

    noteSlotStrand(sexRow, chromosomeCount - 1);					// Y
    noteSlotStrand(sexRow, chromosomeCount);					// X

    for (std::size_t p = 0; p < chromoplexyPlacements.size(); ++p)		// Chromoplexy leftovers are drawn joined into one bar
    {
        double combinedMbp = 0.0;

        for (const SDRdataRecord& record : subHeader.dataRecords)
        {
            const std::vector<int>& leftovers = chromoplexyPlacements[p].leftoverNewStrandIDs;

            if (std::find(leftovers.begin(), leftovers.end(), record.newStrandID) != leftovers.end())
            {
                combinedMbp += recordLengthMbp(record);
            }
        }

        const int row = acentricsStartRow + static_cast<int>(p) / columns;
        rowLongestMbp[row] = std::max(rowLongestMbp[row], combinedMbp);
    }

    std::vector<double> rowTop(static_cast<std::size_t>(totalRows) + 1, startY);	// rowTop[r] = top of row r; rowTop[totalRows] = bottom of the last row

    for (int r = 0; r < totalRows; ++r)
    {
        rowTop[r + 1] = rowTop[r] + std::max(rowHeight, maxRenderHeight * (rowLongestMbp[r] / maxLengthMbp) + rowSpacing);
    }

    const double chromosomeWidth = 20.0;					// Visual width of a single bar/column, safe to tune.
    const double homologGap = 18.0;						// Visual gap between left/right homolog slots, safe to tune.

    // Draw karyogram layout. Everything except the legend is drawn on a recording surface first: the image height is only known once
    // the lowest label has been drawn, since the legend starts a fixed gap below it (see writeKaryogramImage).
    cairo_surface_t* content = cairo_recording_surface_create(CAIRO_CONTENT_COLOR_ALPHA, nullptr);
    cairo_t* cr = cairo_create(content);

    cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);

    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12.0);


    // --------------------------------------------------
    // Summary Box
    // --------------------------------------------------
    drawSDRsummary(cr, masterHeader, rawSubHeader);				// The raw cell, so the box can count its indels.


    // Resolves the records to draw for a given original strand's slot: its own rearranged records if any exist, or a synthesized intact baseline if it was
    // never touched, or nothing if consumed by chromoplexy rearrangement. This check has to run before the map lookup below,
    // a fully-consumed strand still has its own trivial baseline entry in recordsByOriginalStrand (nothing else references that strand 
    // to trigger filterBaselineIfMutated), so checking 'found in the map' first would let that baseline slip through as a false intact chromosome.
    auto resolveSlotRecords = [&](int oldStrandID) -> std::vector<const SDRdataRecord*>
    {
	if (chromoplexyConsumedStrands.count(oldStrandID))
    	{
            return {};
    	}

        const auto it = recordsByOriginalStrand.find(oldStrandID);

        if (it != recordsByOriginalStrand.end())
        {
            return clusterRecordsForDrawing(filterBaselineIfMutated(it->second, chromosomeCount), oldStrandID, mergedRecordStorage, masterHeader);
        }

        // No record of its own, but its centromere sits in another chromosome's record: the whole chromosome was used up by rearrangements (its pieces
        // all sit in other chromosomes' records), so its slot stays empty instead of showing an intact chromosome. Bases
        // merely copied into an indel insertion do not count (the cell without its indels is checked).
        for (const SDRdataRecord& record : subHeader.dataRecords)
        {
            for (const SDRfragment& fragment : record.fragments)
            {
                if (fragment.oldStrandID == oldStrandID && fragment.hasCentromere)	// Its centromere went elsewhere: not still intact
                {
                    return {};
                }
            }
        }

        return synthesizeIntactRecordIfMissing(oldStrandID, subHeader.cellID, masterHeader, intactRecordStorage);
    };



    // Name of every original strand, as its slot is titled, so the zoomed view can say which homolog a translocated segment came from
    // (slots are drawn one after another, so the names are all set before any slot is drawn).
    for (int i = 0; i < drawableGroups; ++i)
    {
        const std::string chromosome = "Chromosome " + std::to_string(i + 1);

        if (chromosomeLayout == ChromosomeLayout::SPLIT_HOMOLOGS)
        {
            detailStrandNames[i + 1] = chromosome + ", homolog 1";
            detailStrandNames[i + 1 + homologousPairs] = chromosome + ", homolog 2";
        }
        else if (chromosomeLayout == ChromosomeLayout::ADJACENT_HOMOLOGS)
        {
            detailStrandNames[2 * i + 1] = chromosome + ", homolog 1";
            detailStrandNames[2 * i + 2] = chromosome + ", homolog 2";
        }
        else
        {
            detailStrandNames[i + 1] = chromosome;
        }
    }
    detailStrandNames[chromosomeCount - 1] = "Chromosome Y";			// Y and X are always the last two strands, in that order
    detailStrandNames[chromosomeCount] = "Chromosome X";

    // ----------------------------------------------------------------------------------------------------
    // Draw each numbered chromosome slot depending on chromosome sizes layout, X and Y will be drawn later
    // ----------------------------------------------------------------------------------------------------
    for (int i = 0; i < drawableGroups; ++i)
    {
        int firstOldStrandID = 0;
        int secondOldStrandID = -1; // -1 = no homolog

        if (chromosomeLayout == ChromosomeLayout::SPLIT_HOMOLOGS)		// [count, 1,2,3...,1,2,3...]
        {
            firstOldStrandID = i + 1;
            secondOldStrandID = i + 1 + homologousPairs;
        }
        else if (chromosomeLayout == ChromosomeLayout::ADJACENT_HOMOLOGS)	// [count,1,1,2,2,3,3...]
        {
            firstOldStrandID = 2 * i + 1;
            secondOldStrandID = 2 * i + 2;
        }
        else									// [count,1,2,3,...]
        {
            firstOldStrandID = i + 1;
        }


        const int col = i % columns;
        const int row = i / columns;
        const double groupCenterX = (col * colWidth) + (colWidth / 2.0);
        const double posY = rowTop[row];

        double labelHeight = maxRenderHeight;                   		// Fallback, if slot has nothing to draw.

	// Homologous chromosome sizes layout specified in SDR header
        if (hasHomologs)
        {

	    // Cluster records for drawing only on mutated records, such that mutations are grouped according to their original strand location
	    // resolveSlotRecords() also leaves a slot empty (rather than synthesizing a false intact baseline) when a
	    // chromoplexy fusion elsewhere has fully consumed that strand.
	    const std::vector<const SDRdataRecord*> leftRecords = resolveSlotRecords(firstOldStrandID);	// First homolog
            const std::vector<const SDRdataRecord*> rightRecords = resolveSlotRecords(secondOldStrandID);	// Second homolog


	    // Adjust the gap of the left and right homologs and their corresponding deletions/ecDNA depending on the actual width drawStackedMutations() will use.
            double leftColumn2Width = chromosomeWidth;
            double rightColumn2Width = chromosomeWidth;
	    // Adjust slot width to account for drawing of deleted segments and ecDNA beside the original chromosomes
            const double leftStackWidth = computeSlotWidth(leftRecords, firstOldStrandID, chromosomeWidth, maxLengthMbp, maxRenderHeight, masterHeader, leftColumn2Width);
            const double rightStackWidth = computeSlotWidth(rightRecords, secondOldStrandID, chromosomeWidth, maxLengthMbp, maxRenderHeight, masterHeader, rightColumn2Width);
            double leftSlotCenterX = groupCenterX - leftStackWidth / 2.0 - homologGap / 2.0;
            double rightSlotCenterX = groupCenterX + rightStackWidth / 2.0 + homologGap / 2.0;

            // Moved-segment arrows are drawn on the LEFT of a chromosome, so on the second homolog they reach into the gap towards the first homolog and
            // its circles. The second homolog is drawn once off-screen to measure how far its drawing reaches left of its slot, and the
            // gap is widened by that much, so the arrows keep the usual homologGap from the first homolog's slot, like the first homolog's circles do.
            if (!rightRecords.empty())
            {
                cairo_surface_t* probeSurface = cairo_recording_surface_create(CAIRO_CONTENT_COLOR_ALPHA, nullptr);
                cairo_t* probe = cairo_create(probeSurface);
                drawStackedMutations(probe, rightRecords, rightSlotCenterX, posY, chromosomeWidth, maxLengthMbp, maxRenderHeight, humanGenome, masterHeader, secondOldStrandID);
                cairo_destroy(probe);

                double inkX = 0.0, inkY = 0.0, inkWidth = 0.0, inkHeight = 0.0;
                cairo_recording_surface_ink_extents(probeSurface, &inkX, &inkY, &inkWidth, &inkHeight);
                cairo_surface_destroy(probeSurface);

                const double overhang = (rightSlotCenterX - rightStackWidth / 2.0) - inkX;	// How far the drawing reaches left of the slot

                if (inkWidth > 0.0 && overhang > 2.0)					// Outlines alone stick out about a pixel; only arrows reach further
                {
                    leftSlotCenterX -= overhang / 2.0;
                    rightSlotCenterX += overhang / 2.0;
                }
            }

	    // After checking for the presence of long deletions and ecDNA drawn beside the original chromosomes, draw them for the left and right homologs.
	    // And draw other mutations stacked on top of the original intact chromosome
            if (!leftRecords.empty())
            {
                beginDetailSlot("Chromosome " + std::to_string(i + 1) + ", homolog 1", firstOldStrandID);
                drawStackedMutations(cr, leftRecords, leftSlotCenterX, posY, chromosomeWidth, maxLengthMbp, maxRenderHeight, humanGenome, masterHeader, firstOldStrandID);
                endDetailSlot();
            }

            if (!rightRecords.empty())
            {
                beginDetailSlot("Chromosome " + std::to_string(i + 1) + ", homolog 2", secondOldStrandID);
                drawStackedMutations(cr, rightRecords, rightSlotCenterX, posY, chromosomeWidth, maxLengthMbp, maxRenderHeight, humanGenome, masterHeader, secondOldStrandID);
                endDetailSlot();
            }

	    // Add a label
            if (!leftRecords.empty() || !rightRecords.empty())
            {
                labelHeight = std::max(
                    computeMaxBarHeight(leftRecords, maxLengthMbp, maxRenderHeight, firstOldStrandID, masterHeader),
                    computeMaxBarHeight(rightRecords, maxLengthMbp, maxRenderHeight, secondOldStrandID, masterHeader)
                );
            }

        }

	// Non-homologous chromosome sizes layout specified by the user in the SDR header (simpler, do not need to adjust gap between homologs for possible long deletions)
        else
        {
	    // resolveSlotRecords() clusters mutated records, falls back to a synthesized intact baseline for untouched strands, and leaves the slot 
	    // empty for a strand fully consumed by a chromoplexy fusion elsewhere.
            const std::vector<const SDRdataRecord*> filtered = resolveSlotRecords(firstOldStrandID);

            if (!filtered.empty())						// If data present, draw the mutation and label it
            {
                beginDetailSlot("Chromosome " + std::to_string(i + 1), firstOldStrandID);
                drawStackedMutations(cr, filtered, groupCenterX, posY, chromosomeWidth, maxLengthMbp, maxRenderHeight, humanGenome, masterHeader, firstOldStrandID);
                endDetailSlot();
                labelHeight = computeMaxBarHeight(filtered, maxLengthMbp, maxRenderHeight, firstOldStrandID, masterHeader);
            }
        }


        // Group label (chromosome number).
        cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 16.0);
        std::string label = "Chr " + std::to_string(i + 1);
        std::vector<int> slotStrandIDs = {firstOldStrandID};		// Strands summarized when the label is clicked in the web viewer
        if (hasHomologs)
        {
            slotStrandIDs.push_back(secondOldStrandID);
        }
        drawChromosomeLabel(cr, label, groupCenterX - 20.0, posY + labelHeight + 25.0, slotStrandIDs);

    }




    // -------------------------------------------------------------------------------------------------------------
    // Drawing X and Y groups, placed at the end of last row if there's room, otherwise they're given their own row.
    // -------------------------------------------------------------------------------------------------------------
    if (chromosomeCount >= 2)
    {
	// Determine positions of X and Y chromosomes depeding on the number of chromosomes specified.
        const double sexChromPosY = rowTop[sexRow];
        int yChromCol;
        int xChromCol;

        if (sexChromosomesFitInLastRow)					// Put sex chromosomes in last two positions of last row if they fit
        {
            yChromCol = lastRowUsedCols;
            xChromCol = lastRowUsedCols + 1;
        }
        else								// Put sex chromosomes in first two positions of a new last row if they don't fit
        {
            yChromCol = 0;
            xChromCol = 1;
        }

	// Y and X are always the last two entries in the chromosome list, in that order, this ordering convention is required, not incidental 
	// (see humanCentromeres in KaryogramLayout.cpp)
        int yOldStrandID = chromosomeCount - 1;
        int xOldStrandID = chromosomeCount;

        double yCenterPosX = (yChromCol * colWidth) + (colWidth / 2.0);
        double xCenterPosX = (xChromCol * colWidth) + (colWidth / 2.0);


	// resolveSlotRecords() clusters mutated records, falls back to a synthesized intact baseline for an untouched sex chromosome, and leaves the 
	// slot empty if a chromoplexy fusion elsewhere fully consumed it.
        const std::vector<const SDRdataRecord*> yRecords = resolveSlotRecords(yOldStrandID);

        double yLabelHeight = maxRenderHeight;
        if (!yRecords.empty())
        {
            beginDetailSlot("Chromosome Y", yOldStrandID);
            drawStackedMutations(cr, yRecords, yCenterPosX, sexChromPosY, chromosomeWidth, maxLengthMbp, maxRenderHeight, humanGenome, masterHeader, yOldStrandID);
            endDetailSlot();
            yLabelHeight = computeMaxBarHeight(yRecords, maxLengthMbp, maxRenderHeight, yOldStrandID, masterHeader);
        }

        cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 16.0);
        drawChromosomeLabel(cr, "Y", yCenterPosX - 10.0, sexChromPosY + yLabelHeight + 25.0, {yOldStrandID});

        const std::vector<const SDRdataRecord*> xRecords = resolveSlotRecords(xOldStrandID);

        double xLabelHeight = maxRenderHeight;

        if (!xRecords.empty())
        {
            beginDetailSlot("Chromosome X", xOldStrandID);
            drawStackedMutations(cr, xRecords, xCenterPosX, sexChromPosY, chromosomeWidth, maxLengthMbp, maxRenderHeight, humanGenome, masterHeader, xOldStrandID);
            endDetailSlot();
            xLabelHeight = computeMaxBarHeight(xRecords, maxLengthMbp, maxRenderHeight, xOldStrandID, masterHeader);
        }

        cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 16.0);
        drawChromosomeLabel(cr, "X", xCenterPosX - 10.0, sexChromPosY + xLabelHeight + 25.0, {xOldStrandID});
    }


    // ------------------------------------------------------------------------------------------------------------------------------------
    // Draw chromoplexy acentric-recombination slots, one per placement, packed into rows after X and Y exactly like the autosome slots are.
    // ------------------------------------------------------------------------------------------------------------------------------------
    for (std::size_t p = 0; p < chromoplexyPlacements.size(); ++p)
    {
        const SDRchromoplexyAcentricPlacement& placement = chromoplexyPlacements[p];

        // Gather the leftover record(s) this placement points to, in the order it lists them.
        std::vector<const SDRdataRecord*> leftoverRecords;
        for (int leftoverNewStrandID : placement.leftoverNewStrandIDs)
        {
            for (const SDRdataRecord& record : subHeader.dataRecords)
            {
                if (record.newStrandID == leftoverNewStrandID)
                {
                    leftoverRecords.push_back(&record);
                    break;
                }
            }
        }

        if (leftoverRecords.empty())
        {
	    // Shouldn't happen, findChromoplexyAcentricRecombinations() only returns newStrandIDs it found in this same subHeader.
            continue;
        }

        // Concatenate every leftover record's fragments into one synthetic record. When the SDR file already combined both partners' leftover pieces 
	// into a single record, this is just that record's own fragments (a concatenation); when they're still two separate records,
        // this joins them. Either way, these are the loose, non-centromeric ends of the same 3+-way rearrangement re-joining each other.
        SDRdataRecord combinedRecord{};
        combinedRecord.cellID = subHeader.cellID;
        combinedRecord.newStrandID = leftoverRecords.front()->newStrandID; 		// Representative ID, used only for the drawn label.
        combinedRecord.linear = true;

        for (const SDRdataRecord* leftoverRecord : leftoverRecords)
        {
            combinedRecord.linear = combinedRecord.linear && leftoverRecord->linear;
            combinedRecord.fragments.insert(combinedRecord.fragments.end(), leftoverRecord->fragments.begin(), leftoverRecord->fragments.end());
        }

        mergedRecordStorage.push_back(combinedRecord);
        const std::vector<const SDRdataRecord*> combinedRecords = { &mergedRecordStorage.back() };

        const int col = static_cast<int>(p) % columns;
        const int row = static_cast<int>(p) / columns;
        const double slotCenterX = (col * colWidth) + (colWidth / 2.0);
        const double slotPosY = rowTop[acentricsStartRow + row];

        // homeOldStrandID is only used by drawStackedMutations()/computeMaxBarHeight() to test for named single-strand shapes (deletion, ecDNA, etc.). 
	// This record mixes two different strands' fragments, so none of those shapes can match regardless of which ID is passed.
        const int chromoplexyDetailSlot = beginDetailSlot("", combinedRecord.fragments.front().oldStrandID);	// Titled below, once the label is known
        drawStackedMutations(cr, combinedRecords, slotCenterX, slotPosY, chromosomeWidth, maxLengthMbp, maxRenderHeight, humanGenome, masterHeader, combinedRecord.fragments.front().oldStrandID);
        endDetailSlot();
        const double slotLabelHeight = computeMaxBarHeight(combinedRecords, maxLengthMbp, maxRenderHeight, combinedRecord.fragments.front().oldStrandID, masterHeader);

        //  Label acentric strand with the fragments that appear on the drawn strand, trailing "*" marking it as acentric.
	std::string chromoplexyLabel;
	for (std::size_t f = 0; f < combinedRecord.fragments.size(); ++f)
	{
    	    if (f > 0)
    	    {
        	chromoplexyLabel += "-";
    	    }
    	    chromoplexyLabel += "chr" + std::to_string(combinedRecord.fragments[f].oldStrandID);
	}
	chromoplexyLabel += "*";

	if (chromoplexyDetailSlot >= 0 && chromoplexyDetailSlot < static_cast<int>(detailSlots.size()))
	{
	    detailSlots[chromoplexyDetailSlot].title = "Acentric piece " + chromoplexyLabel;
	}

	cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_set_font_size(cr, 13.0);

	cairo_text_extents_t labelExtents;
	cairo_text_extents(cr, chromoplexyLabel.c_str(), &labelExtents);
	cairo_move_to(cr, slotCenterX - labelExtents.width / 2.0, slotPosY + slotLabelHeight + 25.0);
	cairo_show_text(cr, chromoplexyLabel.c_str());
    }



    // ------------------------------------------------------------------------------------------------------------------------------------------
    // Draw multi-origin ecDNA slots, one per event, packed into rows after the chromoplexy slots exactly like the autosome/chromoplexy slots are.
    // ------------------------------------------------------------------------------------------------------------------------------------------
    std::size_t ecDNAslotIndex = 0;

    for (const SDRecDNAevent& event : ecDNAevents)
    {
        const std::set<int> distinctOrigins(event.oldStrandIDs.begin(), event.oldStrandIDs.end());

        if (distinctOrigins.size() < 2)
        {
            continue; 									// Same-origin ecDNA is already drawn next to its home chromosome above.
        }

        const SDRdataRecord* circularRecord = nullptr;

        for (const SDRdataRecord& record : subHeader.dataRecords)
        {
            if (record.newStrandID == event.excisedStrandID)
            {
                circularRecord = &record;
                break;
            }
        }

        if (circularRecord == nullptr)
        {
            continue; 									// Shouldn't happen, excisedStrandID came from detectECDNA() on this same subHeader.
        }

        const int slotPosition = chromoplexySlotCount + static_cast<int>(ecDNAslotIndex);
        const int col = slotPosition % columns;
        const int row = slotPosition / columns;
        const double slotCenterX = (col * colWidth) + (colWidth / 2.0);
        const double slotPosY = rowTop[acentricsStartRow + row];
        ++ecDNAslotIndex;

        double totalLengthMbp = 0.0;

        for (const SDRfragment& fragment : circularRecord->fragments)
        {
            totalLengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
        }

        const double diameter = std::min(computeSDRbarHeight(totalLengthMbp, maxLengthMbp, maxRenderHeight), chromosomeWidth * 3.0);
        const double centerY = slotPosY + diameter / 2.0;


	// Sum each contributing strand's total length across however many segments it supplied (a strand with two gaps recombining into this ecDNA gets ONE wedge sized by their
        // combined length, not two separate wedges), then build one wedge per strand, in the order each strand first appears on the circular record, same order the label below
        // uses, so the wedge position and the label token line up.
        std::map<int, double> lengthByStrand;

        for (const SDRfragment& fragment : circularRecord->fragments)
        {
            lengthByStrand[fragment.oldStrandID] += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
        }

        std::vector<std::pair<RGB, double>> wedges;
        std::set<int> seenStrands;

        for (const SDRfragment& fragment : circularRecord->fragments)
        {
            if (seenStrands.count(fragment.oldStrandID))
            {
                continue;
            }

            seenStrands.insert(fragment.oldStrandID);
            wedges.push_back({ getColorForOriginalStrand(fragment.oldStrandID, masterHeader), lengthByStrand[fragment.oldStrandID] });
        }

        drawPieChartFragment(cr, slotCenterX, centerY, diameter, wedges);
        const int ecDNAdetailSlot = beginDetailSlot("", -1);			// Titled below, once the label is known
        recordDetailCircle(*circularRecord, slotCenterX, centerY, diameter);
        endDetailSlot();


        // Label reads e.g. "chrA-chrB*", one "chr<oldStrandID>" token per DISTINCT contributing strand, in the order each strand first appears on the 
	// circular record, trailing "*" marking it as acentric. This names which chromosomes are involved, not how many pieces each gave.
        std::string ecDNAlabel;
        std::set<int> labeledStrands;

        for (int strandID : event.oldStrandIDs)
        {
            if (labeledStrands.count(strandID))
            {
                continue;
            }

            labeledStrands.insert(strandID);

            if (!ecDNAlabel.empty())
            {
                ecDNAlabel += "-";
            }

            ecDNAlabel += "chr" + std::to_string(strandID);
        }

        ecDNAlabel += "*";

        if (ecDNAdetailSlot >= 0 && ecDNAdetailSlot < static_cast<int>(detailSlots.size()))
        {
            detailSlots[ecDNAdetailSlot].title = "ecDNA " + ecDNAlabel;
        }

        cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 13.0);

        cairo_text_extents_t ecDNAlabelExtents;
        cairo_text_extents(cr, ecDNAlabel.c_str(), &ecDNAlabelExtents);
        cairo_move_to(cr, slotCenterX - ecDNAlabelExtents.width / 2.0, slotPosY + diameter + 25.0);
        cairo_show_text(cr, ecDNAlabel.c_str());
    }


    // ---------------------------------------------
    // Draw legend at bottom, a fixed gap below the lowest label, and write the image
    // ---------------------------------------------
    cairo_destroy(cr);

    return writeKaryogramImage(content, imgWidth, &Karyogram::drawSDRlegend, outputFilename, imageBytes, format);
}















// Returns the segments with their fractions remapped so every non-empty segment spans at least minPixels of a bar height pixels tall. The
// extra room is taken proportionally from the larger segments, and centromere fractions move with their segment, so the bar keeps its height.
std::vector<PaintedSegment> Karyogram::applyMinimumSegmentHeight(const std::vector<PaintedSegment>& segments, double height, double minPixels)
{
    if (height <= 0.0)
    {
        return segments;
    }

    const double minFraction = minPixels / height;

    // Boundaries of every segment, sorted, define intervals; each interval inside a segment is a piece to stretch or shrink.
    std::vector<double> bounds;
    for (const PaintedSegment& segment : segments)
    {
        bounds.push_back(std::min(segment.startFraction, segment.endFraction));
        bounds.push_back(std::max(segment.startFraction, segment.endFraction));
    }
    std::sort(bounds.begin(), bounds.end());
    bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());

    if (bounds.size() < 2)
    {
        return segments;
    }

    double tinyTotal = 0.0;
    double restTotal = 0.0;
    std::vector<bool> isTiny(bounds.size() - 1, false);

    for (std::size_t i = 0; i + 1 < bounds.size(); ++i)
    {
        const double length = bounds[i + 1] - bounds[i];

        // Only a stretch some segment paints is enlarged. Pieces listed on neighbouring bases (e.g. 120.000000 then 120.000001) leave a
        // sliver no segment covers; enlarging it to minPixels opened an unpainted band and shifted everything below it.
        bool covered = false;
        for (const PaintedSegment& segment : segments)
        {
            covered = covered || (std::min(segment.startFraction, segment.endFraction) <= bounds[i] && std::max(segment.startFraction, segment.endFraction) >= bounds[i + 1]);
        }

        if (covered && length > 0.0 && length < minFraction)
        {
            isTiny[i] = true;
            tinyTotal += minFraction;
        }
        else
        {
            restTotal += length;
        }
    }

    const double span = bounds.back() - bounds.front();

    if (tinyTotal == 0.0 || restTotal <= 0.0 || tinyTotal >= span)		// Nothing to enlarge, or not enough room to do it.
    {
        return segments;
    }

    const double restScale = (span - tinyTotal) / restTotal;
    std::vector<double> newBounds(bounds.size());
    newBounds[0] = bounds[0];

    for (std::size_t i = 0; i + 1 < bounds.size(); ++i)
    {
        newBounds[i + 1] = newBounds[i] + (isTiny[i] ? minFraction : (bounds[i + 1] - bounds[i]) * restScale);
    }

    // Piecewise-linear map from an old fraction to its new position.
    const auto remap = [&](double fraction)
    {
        if (fraction <= bounds.front())
        {
            return newBounds.front() + (fraction - bounds.front());
        }

        for (std::size_t i = 0; i + 1 < bounds.size(); ++i)
        {
            if (fraction <= bounds[i + 1])
            {
                const double length = bounds[i + 1] - bounds[i];
                const double t = (length > 0.0) ? (fraction - bounds[i]) / length : 0.0;
                return newBounds[i] + t * (newBounds[i + 1] - newBounds[i]);
            }
        }

        return newBounds.back() + (fraction - bounds.back());
    };

    std::vector<PaintedSegment> result = segments;

    for (PaintedSegment& segment : result)
    {
        segment.startFraction = remap(segment.startFraction);
        segment.endFraction = remap(segment.endFraction);

        if (segment.hasCentromere)
        {
            segment.centromereStartFraction = remap(segment.centromereStartFraction);
            segment.centromereEndFraction = remap(segment.centromereEndFraction);
        }
    }

    return result;
}



void Karyogram::drawPaintedFragment(			// Function to draw the individual fragments making up the chromosomes on the karyogram and the excised linear fragments (circular fragments handled later)
    cairo_t* cr,
    double x, 						// X position of the painted segment
    double y, 						// Y position of the painted segment
    double height, 					// Height of the chromosome
    double width, 					// Width of the chromosome
    const std::vector<PaintedSegment>& inputSegments,	// Painted segment vector to add the aberrant strands on the original chromosomes to depict structural variations
    bool roundTopCap,					// Whether the top edge gets a rounded telomere cap or a flat cut edge
    bool roundBottomCap,				// Whether the bottom edge gets a rounded telomere cap or a flat cut edge
    bool drawOutline)					// Whether the shape should have a black outline, default to true
{
    if (inputSegments.empty())
    {
        return;
    }

    // Tiny segments (e.g. a 0.1 Mbp translocated piece) would round away to nothing, so every segment gets a minimum height. 3 px rather than 1,
    // because a 1 px band at a rounded end hides under the 1.5 px outline.
    const double minSegmentPixels = 3.0;
    const std::vector<PaintedSegment> segments = applyMinimumSegmentHeight(inputSegments, height, minSegmentPixels);

    const double capRadius = width / 2.0;



    // --------------------------------------------------
    // Determine centromere locations
    // --------------------------------------------------
    struct CentromereSpan									// Local struct only available for the SDR chromosome drawings
    {
        double centromereTopY;									// Centromere top height in pixels
        double centromereBottomY;								// Centromere bottom height in pixels
    };

    std::vector<CentromereSpan> centromereSpans;
    for (const PaintedSegment& segment : segments)
    {
        if (segment.hasCentromere)								// Check if fragments/strands have a centromere
        {
            CentromereSpan span{};								// Initialize the CentromereSpan vector span
            span.centromereTopY = y + height * segment.centromereStartFraction;			// Top of centromere in pixels
            span.centromereBottomY = y + height * segment.centromereEndFraction;		// Bottom of centromere in pixels
            centromereSpans.push_back(span);							// Store this centromere's geometry
        }
    }

    // Sort centromere spans in order of height (top before bottom)
    std::sort(centromereSpans.begin(), centromereSpans.end(), [](const CentromereSpan& a, const CentromereSpan& b)
    {
        return a.centromereTopY < b.centromereTopY;
    });

    // Draw constriction of the chromosome near the centromere center.
    const double constrictionHeight = 7.0;
    const double constrictionAmount = 2.0;


    // Trace the outer capsule outline with a constriction notch at each centromere position along the multiple fragments making up a mutated chromosome
    const auto traceCapsulePath = [&]()
    {
        cairo_new_path(cr);

        if (roundTopCap)
        {
            cairo_arc(cr, x + capRadius, y + capRadius, capRadius, M_PI, 2 * M_PI);
        }
        else
        {
            cairo_move_to(cr, x, y);
            cairo_line_to(cr, x + width, y);
        }

        for (const CentromereSpan& span : centromereSpans)
        {
            const double centromereY = (span.centromereTopY + span.centromereBottomY) / 2.0;
            const double constrictionTop = centromereY - constrictionHeight / 2.0;
            const double constrictionBottom = centromereY + constrictionHeight / 2.0;

            cairo_line_to(cr, x + width, constrictionTop);
            cairo_line_to(cr, x + width - constrictionAmount, centromereY);
            cairo_line_to(cr, x + width, constrictionBottom);
        }

        const double bottomRightY = roundBottomCap ? (y + height - capRadius) : (y + height);
        cairo_line_to(cr, x + width, bottomRightY);

        if (roundBottomCap)
        {
            cairo_arc(cr, x + capRadius, y + height - capRadius, capRadius, 0, M_PI);
        }
        else
        {
            cairo_line_to(cr, x, y + height);
        }

        for (auto it = centromereSpans.rbegin(); it != centromereSpans.rend(); ++it)
        {
            const double centromereY = (it->centromereTopY + it->centromereBottomY) / 2.0;
            const double constrictionTop = centromereY - constrictionHeight / 2.0;
            const double constrictionBottom = centromereY + constrictionHeight / 2.0;

            cairo_line_to(cr, x, constrictionBottom);
            cairo_line_to(cr, x + constrictionAmount, centromereY);
            cairo_line_to(cr, x, constrictionTop);
        }

        cairo_close_path(cr);
    };


    traceCapsulePath(); 									// build it once, for the clip


    // --------------------------------------------------------
    // Fill each painted segment, clipped to the capsule shape
    // --------------------------------------------------------

    // Draw white segments on chromosomes where deletions occurred instead of truncating the chromosome, for clarity.
    const auto isWhiteColor = [](const RGB& color)
    {
        const double tolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)
        return std::fabs(color.r - 1.0) < tolerance && std::fabs(color.g - 1.0) < tolerance && std::fabs(color.b - 1.0) < tolerance;
    };

    const double minWhiteGapHeight = 1.0;	// Keep a deletion's white gap from shrinking below a visible size, however small the actual deleted region is.

    cairo_save(cr);
    cairo_clip_preserve(cr);

    for (const PaintedSegment& segment : segments)
    {
        const double segmentTop = y + height * segment.startFraction;
        const double segmentBottom = y + height * segment.endFraction;

        cairo_set_source_rgb(cr, segment.color.r, segment.color.g, segment.color.b);
        cairo_rectangle(cr, x, segmentTop, width, segmentBottom - segmentTop);
        cairo_fill(cr);
    }

    // A white (deletion-gap) segment can shrink to a fraction of a pixel for a tiny deletion and effectively disappear. Re-painting it, re-centered on
    // its original midpoint and padded out to the minimum height, in a SEPARATE pass after every segment has already been filled once above ensures the
    // padding isn't immediately painted back over by whichever neighboring segment gets drawn next in the main loop. Non-white segments are untouched,
    // and nothing here changes the fraction boundaries other segments or the centromere markers use, only how tall this one gap's fill rectangle is.
    for (const PaintedSegment& segment : segments)
    {
        if (!isWhiteColor(segment.color))
        {
            continue;
        }

        const double segmentTop = y + height * segment.startFraction;
        const double segmentBottom = y + height * segment.endFraction;

        if ((segmentBottom - segmentTop) >= minWhiteGapHeight)
        {
            continue;
        }

        const double gapCenter = (segmentTop + segmentBottom) / 2.0;
        const double paddedTop = gapCenter - minWhiteGapHeight / 2.0;
        const double paddedBottom = gapCenter + minWhiteGapHeight / 2.0;

        cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
        cairo_rectangle(cr, x, paddedTop, width, paddedBottom - paddedTop);
        cairo_fill(cr);
    }

    std::vector<std::pair<double, double>> chevronsPx;					// Top and bottom of each chevron, for the zoomed view

    // Inversion chevrons, drawn in their own pass after every fill (and white-gap padding) so no neighbouring segment can paint over them. Each chevron is black and spans the
    // full height of its inverted segment and the full width of the chromosome, so a short inversion gives a flatter and flatter arrow, ending as a single horizontal line.
    // Still inside the capsule clip, so the arm ends are trimmed to the chromosome outline.
    // A long deletion inside an inverted segment (reversed piece, white gap, reversed piece) gets ONE chevron over the whole segment, drawn on top of the white gap.
    for (std::size_t i = 0; i < segments.size(); ++i)
    {
        const PaintedSegment& segment = segments[i];

        if (!segment.isReversed || isWhiteColor(segment.color))
        {
            continue;
        }

        double endFraction = segment.endFraction;

        while (i + 2 < segments.size() && segments[i + 1].insideInversion && segments[i + 2].isReversed && !isWhiteColor(segments[i + 2].color))
        {
            i += 2;
            endFraction = segments[i].endFraction;
        }

        drawInversionChevron(cr, x + width / 2.0, y + height * segment.startFraction, y + height * endFraction, width);
        chevronsPx.emplace_back(y + height * segment.startFraction, y + height * endFraction);
    }

    cairo_restore(cr);


    // ------------------------------------------------
    // Pieces moved within their own chromosome (intra-chromosomal translocations)
    // ------------------------------------------------
    // Pieces of one original strand are grouped into blocks: an inverted stretch (a reversed piece and the pieces that keep reading downward after it, across
    // white gaps too) is one block, so plain and nested inversions never count as moves. The blocks of each strand that stay in their original order are the heaviest increasing run of their
    // lowest positions (the most material left in place); every other block was moved. Each moved block gets a red hooked arrow on the left of the
    // chromosome (deletion pieces and ecDNA are drawn on the right) from the spot it was cut out of, to its own middle, and its
    // borders at its new place are drawn dashed (the place it was cut out of gets no line at all). Only for whole chromosomes (drawOutline), not for excised pieces.
    struct PieceBlock
    {
        int strandID;
        double lowMbp;
        double highMbp;
        std::size_t first;								// First and last segment index of the block
        std::size_t last;
        bool swappedWithTrailingPiece = false;					// An inverted block whose trailing forward piece(s) were split off (see below)
        double keyMbp = 0.0;								// Position the block is ordered by (its lowest, or the place it took, see below)
    };

    std::vector<bool> dashedBoundary(segments.size() + 1, false);			// Boundary k lies between segments k-1 and k
    std::vector<bool> hiddenBoundary(segments.size() + 1, false);			// Where a moved piece was cut out: no line at all
    std::vector<PieceBlock> movedBlocks;
    std::vector<std::pair<double, double>> movedArrows;					// fromY, toY in pixels

    if (drawOutline)
    {
        std::vector<PieceBlock> blocks;

        // Is any other piece of this strand lying between two positions? Then two pieces next to each other on the drawing were not neighbours originally.
        const auto pieceLiesBetween = [&](int strandID, double lowMbp, double highMbp, std::size_t skipA, std::size_t skipB)
        {
            for (std::size_t k = 0; k < segments.size(); ++k)
            {
                if (k == skipA || k == skipB || segments[k].oldStrandID != strandID || isWhiteColor(segments[k].color))
                {
                    continue;
                }

                const double otherLow = std::min(segments[k].oldStartPosition, segments[k].oldEndPosition);
                const double otherHigh = std::max(segments[k].oldStartPosition, segments[k].oldEndPosition);

                if (otherHigh > lowMbp + svToleranceMbp && otherLow < highMbp - svToleranceMbp)
                {
                    return true;
                }
            }
            return false;
        };

        for (std::size_t i = 0; i < segments.size(); ++i)
        {
            const PaintedSegment& segment = segments[i];

            if (segment.oldStrandID < 0 || isWhiteColor(segment.color))
            {
                continue;
            }

            const double lowMbp = std::min(segment.oldStartPosition, segment.oldEndPosition);
            const double highMbp = std::max(segment.oldStartPosition, segment.oldEndPosition);

            // An inverted block starts with a reversed piece and keeps reading downward: each next piece lies wholly below the previous one with nothing of the
            // strand in between (a piece flipped back inside the inversion, a nested inversion, stays part of it).
            if (!blocks.empty())
            {
                PieceBlock& previous = blocks.back();
                const PaintedSegment& previousSegment = segments[previous.last];
                const double previousLow = std::min(previousSegment.oldStartPosition, previousSegment.oldEndPosition);
                bool onlyWhiteBetween = true;

                for (std::size_t j = previous.last + 1; j < i; ++j)
                {
                    onlyWhiteBetween = onlyWhiteBetween && isWhiteColor(segments[j].color);
                }

                if (onlyWhiteBetween && previous.strandID == segment.oldStrandID && segments[previous.first].isReversed
                    && highMbp <= previousLow + svToleranceMbp && !pieceLiesBetween(segment.oldStrandID, highMbp, previousLow, i, previous.last))
                {
                    previous.lowMbp = std::min(previous.lowMbp, lowMbp);
                    previous.highMbp = std::max(previous.highMbp, highMbp);
                    previous.last = i;
                    continue;
                }
            }

            blocks.push_back({segment.oldStrandID, lowMbp, highMbp, i, i});
        }

        // An inverted block also ends with a reversed piece: forward pieces trailing it were not flipped back inside an inversion but
        // swapped with it (the forward piece lies below the reversed one, so it moved). They become
        // blocks of their own.
        for (std::size_t b = 0; b < blocks.size(); ++b)
        {
            PieceBlock& block = blocks[b];

            if (!segments[block.first].isReversed || segments[block.last].isReversed)
            {
                continue;
            }

            std::size_t newLast = block.last;
            while (newLast > block.first && (!segments[newLast].isReversed || isWhiteColor(segments[newLast].color)))
            {
                --newLast;
            }

            std::vector<PieceBlock> trailing;
            for (std::size_t k = newLast + 1; k <= block.last; ++k)
            {
                if (!isWhiteColor(segments[k].color))
                {
                    const double lowMbp = std::min(segments[k].oldStartPosition, segments[k].oldEndPosition);
                    const double highMbp = std::max(segments[k].oldStartPosition, segments[k].oldEndPosition);
                    trailing.push_back({segments[k].oldStrandID, lowMbp, highMbp, k, k});
                }
            }

            block.last = newLast;
            block.swappedWithTrailingPiece = !trailing.empty();
            block.lowMbp = std::numeric_limits<double>::max();
            block.highMbp = std::numeric_limits<double>::lowest();
            for (std::size_t k = block.first; k <= block.last; ++k)
            {
                if (!isWhiteColor(segments[k].color))
                {
                    block.lowMbp = std::min({block.lowMbp, segments[k].oldStartPosition, segments[k].oldEndPosition});
                    block.highMbp = std::max({block.highMbp, segments[k].oldStartPosition, segments[k].oldEndPosition});
                }
            }

            blocks.insert(blocks.begin() + static_cast<long>(b) + 1, trailing.begin(), trailing.end());
            b += trailing.size();
        }

        for (PieceBlock& block : blocks)
        {
            block.keyMbp = block.lowMbp;

            // An inverted block whose own original place is a clean cut between two pieces drawn next to each other was cut out from there: it is the
            // moved one (see movedFragmentFlags()).
            if (block.swappedWithTrailingPiece)
            {
                std::size_t previousPiece = segments.size();
                for (std::size_t k = 0; k < segments.size(); ++k)
                {
                    if (isWhiteColor(segments[k].color))
                    {
                        continue;
                    }

                    if (previousPiece < segments.size() && !(previousPiece >= block.first && previousPiece <= block.last))
                    {
                        const PaintedSegment& above = segments[previousPiece];
                        const PaintedSegment& below = segments[k];
                        const double aboveLow = std::min(above.oldStartPosition, above.oldEndPosition);
                        const double aboveHigh = std::max(above.oldStartPosition, above.oldEndPosition);
                        const double belowLow = std::min(below.oldStartPosition, below.oldEndPosition);
                        const double belowHigh = std::max(below.oldStartPosition, below.oldEndPosition);

                        if (above.oldStrandID == block.strandID && below.oldStrandID == block.strandID
                            && ((approxEqual(aboveHigh, block.lowMbp, svToleranceMbp) && approxEqual(belowLow, block.highMbp, svToleranceMbp))
                                || (approxEqual(aboveLow, block.highMbp, svToleranceMbp) && approxEqual(belowHigh, block.lowMbp, svToleranceMbp))))
                        {
                            block.swappedWithTrailingPiece = false;
                            break;
                        }
                    }

                    previousPiece = k;
                }
            }
        }

        // An inverted block drawn right before a forward piece of its strand that originally lay below it (not joined above because other pieces of the
        // strand lie in between) also swapped places with it (see movedFragmentFlags()): it keeps the place it took, which
        // orders it, and the forward piece gets the arrow.
        for (std::size_t b = 0; b + 1 < blocks.size(); ++b)
        {
            PieceBlock& block = blocks[b];
            const PieceBlock& next = blocks[b + 1];
            const PaintedSegment& nextSegment = segments[next.first];

            if (block.swappedWithTrailingPiece || !segments[block.first].isReversed || !segments[block.last].isReversed || next.first != next.last
                || next.first != block.last + 1 || nextSegment.isReversed || next.strandID != block.strandID || next.highMbp > block.lowMbp + svToleranceMbp)
            {
                continue;
            }

            // Not when the inversion's own place is a clean cut between two pieces drawn next to each other (it was moved out from between them).
            bool leftCleanCut = false;
            std::size_t previousPiece = segments.size();
            for (std::size_t k = 0; k < segments.size() && !leftCleanCut; ++k)
            {
                if (isWhiteColor(segments[k].color))
                {
                    continue;
                }

                if (previousPiece < segments.size())
                {
                    const PaintedSegment& above = segments[previousPiece];
                    const PaintedSegment& below = segments[k];
                    leftCleanCut = above.oldStrandID == block.strandID && below.oldStrandID == block.strandID && !above.isReversed && !below.isReversed
                                   && approxEqual(std::max(above.oldStartPosition, above.oldEndPosition), block.lowMbp, svToleranceMbp)
                                   && approxEqual(std::min(below.oldStartPosition, below.oldEndPosition), block.highMbp, svToleranceMbp);
                }

                previousPiece = k;
            }

            // And only when the strand reads on upward from the forward piece up to the inversion's upper end (see movedFragmentFlags()).
            bool readsOnToInversion = false;
            std::size_t current = next.first;
            for (std::size_t k = next.first + 1; k < segments.size(); ++k)
            {
                if (isWhiteColor(segments[k].color))
                {
                    continue;
                }

                const PaintedSegment& above = segments[current];
                const PaintedSegment& following = segments[k];
                const double followingLow = std::min(following.oldStartPosition, following.oldEndPosition);

                if (following.oldStrandID != above.oldStrandID || following.isReversed || followingLow < std::max(above.oldStartPosition, above.oldEndPosition) - svToleranceMbp)
                {
                    break;
                }

                if (approxEqual(followingLow, block.highMbp, svToleranceMbp))
                {
                    readsOnToInversion = true;
                    break;
                }

                if (std::max(following.oldStartPosition, following.oldEndPosition) > block.lowMbp + svToleranceMbp)
                {
                    break;
                }

                current = k;
            }

            if (leftCleanCut || !readsOnToInversion)
            {
                continue;
            }

            block.swappedWithTrailingPiece = true;
            block.keyMbp = next.lowMbp;
        }

        std::map<int, std::vector<std::size_t>> blocksByStrand;
        for (std::size_t b = 0; b < blocks.size(); ++b)
        {
            blocksByStrand[blocks[b].strandID].push_back(b);
        }

        const auto blockMiddleY = [&](const PieceBlock& block)
        {
            return y + height * (segments[block.first].startFraction + segments[block.last].endFraction) / 2.0;
        };

        for (const auto& [strandID, indices] : blocksByStrand)
        {
            const std::size_t n = indices.size();

            if (n < 2)
            {
                continue;
            }

            // Heaviest increasing run of lowest positions, weighted by block length (plus one base, so single-base pieces count too)
            const double oneBaseMbp = 0.000001;
            std::vector<double> best(n, 0.0);
            std::vector<int> previousInRun(n, -1);

            // An inverted block swapped with the forward piece trailing it keeps its place, and the forward piece is the moved one: the inversion already shows
            // with its chevron, so the arrow marks the forward piece (originally before the inverted one).
            const double invertedBlockBonusMbp = 1000.0;
            const auto blockWeight = [&](const PieceBlock& block)
            {
                return block.highMbp - block.lowMbp + oneBaseMbp + (block.swappedWithTrailingPiece ? invertedBlockBonusMbp : 0.0);
            };

            for (std::size_t a = 0; a < n; ++a)
            {
                const PieceBlock& blockA = blocks[indices[a]];
                best[a] = blockWeight(blockA);

                for (std::size_t c = 0; c < a; ++c)
                {
                    const PieceBlock& blockC = blocks[indices[c]];
                    const double candidate = best[c] + blockWeight(blockA);

                    if (blockC.keyMbp < blockA.keyMbp && candidate > best[a])
                    {
                        best[a] = candidate;
                        previousInRun[a] = static_cast<int>(c);
                    }
                }
            }

            std::vector<bool> inPlace(n, false);
            int runEnd = static_cast<int>(std::max_element(best.begin(), best.end()) - best.begin());
            while (runEnd != -1)
            {
                inPlace[runEnd] = true;
                runEnd = previousInRun[runEnd];
            }

            for (std::size_t a = 0; a < n; ++a)
            {
                if (inPlace[a])
                {
                    continue;
                }

                const PieceBlock& moved = blocks[indices[a]];

                // The block it used to sit next to: the one starting closest to where it ended (it sat above it), else the one ending closest to where it started.
                int originBlock = -1;
                double closest = 0.0;

                for (std::size_t c = 0; c < n; ++c)
                {
                    if (c == a)
                    {
                        continue;
                    }

                    const PieceBlock& other = blocks[indices[c]];
                    const double distanceBelow = std::fabs(other.lowMbp - moved.highMbp);
                    const double distanceAbove = std::fabs(moved.lowMbp - other.highMbp);

                    if (originBlock == -1 || distanceBelow < closest)
                    {
                        originBlock = static_cast<int>(c);
                        closest = distanceBelow;
                    }

                    if (distanceAbove < closest)
                    {
                        originBlock = static_cast<int>(c);
                        closest = distanceAbove;
                    }
                }

                const PieceBlock& origin = blocks[indices[originBlock]];

                // The arrow starts exactly where the piece was cut out: the end of the former neighbour that touched it (top of a segment for its listed
                // start, bottom for its listed end).
                double fromY = blockMiddleY(origin);
                double closestEnd = -1.0;

                for (std::size_t k = origin.first; k <= origin.last; ++k)
                {
                    const PaintedSegment& piece = segments[k];

                    if (isWhiteColor(piece.color))
                    {
                        continue;
                    }

                    const double startDistance = std::min(std::fabs(piece.oldStartPosition - moved.lowMbp), std::fabs(piece.oldStartPosition - moved.highMbp));
                    const double endDistance = std::min(std::fabs(piece.oldEndPosition - moved.lowMbp), std::fabs(piece.oldEndPosition - moved.highMbp));

                    if (closestEnd < 0.0 || startDistance < closestEnd)
                    {
                        closestEnd = startDistance;
                        fromY = y + height * piece.startFraction;
                    }

                    if (endDistance < closestEnd)
                    {
                        closestEnd = endDistance;
                        fromY = y + height * piece.endFraction;
                    }
                }

                movedArrows.emplace_back(fromY, blockMiddleY(moved));
                movedBlocks.push_back(moved);

                dashedBoundary[moved.first] = true;
                dashedBoundary[moved.last + 1] = true;
            }
        }
    }


    // The spot a moved piece was cut out of: the two pieces that used to flank it now meet there (one ends where the moved piece started, the other starts
    // where it ended). That border is left without a line, so only the new place is marked.
    for (std::size_t k = 1; k < segments.size(); ++k)
    {
        const PaintedSegment& above = segments[k - 1];
        const PaintedSegment& below = segments[k];

        if (dashedBoundary[k] || above.oldStrandID < 0 || above.oldStrandID != below.oldStrandID || isWhiteColor(above.color) || isWhiteColor(below.color))
        {
            continue;
        }

        const double gapLow = std::min(above.oldEndPosition, below.oldStartPosition);
        const double gapHigh = std::max(above.oldEndPosition, below.oldStartPosition);

        for (const PieceBlock& moved : movedBlocks)
        {
            if (moved.strandID == above.oldStrandID && approxEqual(moved.lowMbp, gapLow, svToleranceMbp) && approxEqual(moved.highMbp, gapHigh, svToleranceMbp))
            {
                hiddenBoundary[k] = true;
            }
        }
    }


    // ------------------------------------------------
    // Segment boundaries
    // ------------------------------------------------

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_set_line_width(cr, 1.0);

    // Skip boundary lines touching a white (deletion gap) segment, a thin gap with black lines on both edges reads as solid black instead of white. Inverted segments are
    // marked by their chevron alone, no outlining lines.
    // The borders of a moved piece at its new place are dashed, also next to an inverted piece or a white deletion gap.
    const double junctionDash[2] = { 3.0, 2.0 };

    // Boundary lines are clipped to the chromosome's own shape, so a junction inside a rounded end (a short piece next to the telomere) or at the
    // centromere pinch never sticks out past the outline.
    cairo_save(cr);
    traceCapsulePath();
    cairo_clip(cr);

    for (std::size_t i = 0; i + 1 < segments.size(); i++)
    {
        const bool dashed = dashedBoundary[i + 1];

        // Next to a white gap only a moved piece's dashed border is drawn (between two white segments never)
        if ((isWhiteColor(segments[i].color) || isWhiteColor(segments[i + 1].color)) && (!dashed || (isWhiteColor(segments[i].color) && isWhiteColor(segments[i + 1].color))))
        {
            continue;
        }

        if (!dashed && (segments[i].isReversed || segments[i + 1].isReversed || hiddenBoundary[i + 1]))
        {
            continue;
        }

        const double boundaryY = y + height * segments[i].endFraction;

        if (dashed)
        {
            cairo_set_dash(cr, junctionDash, 2, 0.0);
        }

        cairo_move_to(cr, x, boundaryY);
        cairo_line_to(cr, x + width, boundaryY);
        cairo_stroke(cr);

        if (dashed)
        {
            cairo_set_dash(cr, nullptr, 0, 0.0);
        }
    }

    cairo_restore(cr);


    // --------------------------------------
    // Chromosome outline
    // --------------------------------------
    if (drawOutline)
    {
        // Outline style follows the material: stretches from the second homolog are dashed, everything else solid, so material exchanged between homologs keeps the
        // outline of the chromosome it came from. White gap segments carry no origin, so they take the style of the nearest preceding real segment (else the next one).
        std::vector<bool> dashed(segments.size(), false);
        for (std::size_t i = 0; i < segments.size(); i++)
        {
            if (!isWhiteColor(segments[i].color))
            {
                dashed[i] = segments[i].isSecondHomolog;
                continue;
            }

            bool found = false;
            for (std::size_t j = i; j-- > 0 && !found;)
            {
                if (!isWhiteColor(segments[j].color)) { dashed[i] = segments[j].isSecondHomolog; found = true; }
            }
            for (std::size_t j = i + 1; j < segments.size() && !found; j++)
            {
                if (!isWhiteColor(segments[j].color)) { dashed[i] = segments[j].isSecondHomolog; found = true; }
            }
        }

        // Merge consecutive same-style segments into runs, then stroke the whole outline once per run, clipped to that run's vertical span (first/last run reach past the caps).
        const double dashPattern[2] = { 4.0, 3.0 };
        const double big = 1.0e6;
        std::size_t runStart = 0;

        while (runStart < segments.size())
        {
            std::size_t runEnd = runStart;
            while (runEnd + 1 < segments.size() && dashed[runEnd + 1] == dashed[runStart])
            {
                runEnd++;
            }

            const double clipTop = (runStart == 0) ? -big : y + height * segments[runStart].startFraction;
            const double clipBottom = (runEnd + 1 == segments.size()) ? big : y + height * segments[runEnd].endFraction;

            cairo_save(cr);
            cairo_rectangle(cr, x - 10.0, clipTop, width + 20.0, clipBottom - clipTop);
            cairo_clip(cr);

            traceCapsulePath(); // Rebuild. The previous stroke / clip consumed the path

            cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
            cairo_set_line_width(cr, 1.5);

            if (dashed[runStart])
            {
                cairo_set_dash(cr, dashPattern, 2, 0.0);
            }

            cairo_stroke(cr);
            cairo_restore(cr);

            runStart = runEnd + 1;
        }
    }


    // --------------------------------------
    // Draw centromere(s) per strand
    // --------------------------------------
    const double ellipseWidth = width + 2.0;

    for (const CentromereSpan& span : centromereSpans)
    {
        const double centromereEllipseCenterY = (span.centromereTopY + span.centromereBottomY) / 2.0;
        const double ellipseHeight = std::max(5.0, span.centromereBottomY - span.centromereTopY);

        cairo_save(cr);
        cairo_translate(cr, x + width / 2.0, centromereEllipseCenterY);
        cairo_scale(cr, ellipseWidth / 2.0, ellipseHeight / 2.0);
        cairo_arc(cr, 0.0, 0.0, 1.0, 0.0, 2.0 * M_PI);
        cairo_restore(cr);

        cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
        cairo_fill(cr);
    }


    // --------------------------------------
    // Arrows of pieces moved within the chromosome, drawn last so nothing covers them
    // --------------------------------------
    for (std::size_t i = 0; i < movedArrows.size(); ++i)
    {
        const double fromY = movedArrows[i].first;
        const double toY = movedArrows[i].second;
        const double bulge = std::max(11.0, 0.6 * std::fabs(toY - fromY)) + 6.0 * static_cast<double>(i);	// Further arrows reach a little further out

        drawMovedSegmentArrow(cr, x, fromY, toY, bulge, 1.3, 5.5, 3.6);
    }

    if (currentDetailSlot >= 0)
    {
        recordDetailBar(x, y, height, width, inputSegments, segments, chevronsPx, movedArrows, dashedBoundary, hiddenBoundary, roundTopCap, roundBottomCap, drawOutline);
    }
}















// Helper to draw rearranged fragments on their original strands, stacked on top or beside in the same chromosome slot.
void Karyogram::drawStackedMutations(
    cairo_t* cr,
    const std::vector<const SDRdataRecord*>& records, 					// Access SDR data records to determine fragment origins and sizes
    double slotCenterX, 								// Determine the slot center to group the intact chromosomes and their associated mutations
    double posY, 									// Generic Y position variable to determine Y height of the chromosome group
    double chromosomeWidth, 								// Chromosome width shout match fragment width
    double maxLengthMbp, 								// Max length to normalize the heights of all the chromosomes drawn
    double maxRenderHeight, 								// Put a cap on the height of the chromosome heights
    bool humanGenome, 									// Check if user specified a human genome.
    const SDRmasterHeader& masterHeader,						// Check SDR master header for intact chromosome sizes layout
    int homeOldStrandID)								// The original chromosome this slot belongs to
{
    if (records.empty())
    {
        return;
    }


    // A long deletion's pair of records gets special treatment: the remaining piece is drawn at the ORIGINAL chromosome's full
    // length with a white gap marking the deletion, and the excised piece gets flat (non-rounded) caps on whichever end isn't a
    // real chromosome edge.
    // Unified gap shape: a home-strand-only slot with one remaining record carrying gap(s), plus any OPTIONAL excised pieces (linear) or ecDNA (circular).
    // Covers long deletions with or without an excised record, deletion-inversions, same-origin ecDNA and lone gaps (deletion-insertion donors).
    const SDRdataRecord* gapRemainingRecord = nullptr;
    std::vector<const SDRdataRecord*> gapPieceRecords;
    const bool isGapShape = analyzeHomeStrandGapShape(records, homeOldStrandID, masterHeader, gapRemainingRecord, gapPieceRecords);
    const bool isChromosomeRing = isChromosomeRingShape(records, homeOldStrandID, masterHeader);

    const std::size_t sizeIndex = static_cast<std::size_t>(homeOldStrandID);
    const double originalSizeMbp = (sizeIndex < masterHeader.intactChromosomeSizes.size()) ? masterHeader.intactChromosomeSizes[sizeIndex] : 0.0;

    // Determine width of a given homologous pair drawing slot
    const double stackGap = 6.0;
    double column2Width = chromosomeWidth;					// Default for shapes with only one column; overwritten below for deletion-like/ecDNA shapes with a second column.
    const double totalWidth = computeSlotWidth(records, homeOldStrandID, chromosomeWidth, maxLengthMbp, maxRenderHeight, masterHeader, column2Width);


    if (isGapShape && originalSizeMbp > 0.0)
    {
        // The remaining record is drawn at the ORIGINAL chromosome's full length, a white gap wherever material is missing (internal or terminal), reversed
        // fragments carry the inversion chevron. Every gap is drawn whether or not an excised partner exists in the file.
        const double remainingBarHeight = computeSDRbarHeight(originalSizeMbp, maxLengthMbp, maxRenderHeight);
        const double remainingBarX = slotCenterX - totalWidth / 2.0;
        const std::vector<PaintedSegment> segments = buildDeletionRemainingSegments(*gapRemainingRecord, homeOldStrandID, humanGenome, masterHeader);

        drawPaintedFragment(cr, remainingBarX, posY, remainingBarHeight, chromosomeWidth, segments);

        // Optional second column: ONLY the excised pieces actually present in the file. Linear pieces get flat caps unless they sit at a real chromosome end,
        // circular pieces are ecDNA circles. No pieces -> no second column and nothing drawn beside the chromosome.
        const bool hasCircularPiece = std::any_of(gapPieceRecords.begin(), gapPieceRecords.end(), [](const SDRdataRecord* piece) { return !piece->linear; });
        const double excisedColumnCenterX = remainingBarX + chromosomeWidth + (hasCircularPiece ? circleColumnGap : stackGap) + column2Width / 2.0;
        const double verticalGap = 6.0;
        const double delTolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)
        double excisedY = posY;

        for (const SDRdataRecord* piece : gapPieceRecords)
        {
            if (!piece->linear)
            {
                double pieceLengthMbp = 0.0;

                for (const SDRfragment& fragment : piece->fragments)
                {
                    pieceLengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
                }

                // Diameter scaled to the fragment's length (not circumference), capped at 3x chromosome width, floored so a tiny circle stays visible.
                const double minPieceSize = 5.0;			// Matches the floor computeSlotWidth() applies when it sizes this slot's column.
                const double diameter = std::max(minPieceSize, std::min((pieceLengthMbp / originalSizeMbp) * remainingBarHeight, chromosomeWidth * 3.0));
                const double centerY = excisedY + diameter / 2.0;

                // A circular piece that kept the centromere is a ring chromosome, not ecDNA.
                drawSingleOriginCircularPiece(cr, *piece, excisedColumnCenterX, centerY, diameter, homeOldStrandID, humanGenome, masterHeader);

                excisedY += diameter + verticalGap;
            }
            else
            {
                const SDRfragment& fragment = piece->fragments[0];
                const double fragLower = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
                const double fragHigher = std::max(fragment.oldStartPosition, fragment.oldEndPosition);

                const bool roundTopCap = approxEqual(fragLower, 0.0, delTolerance);
                const bool roundBottomCap = approxEqual(fragHigher, originalSizeMbp, delTolerance);

                // Same fraction of the remaining bar height as its length is of the original chromosome, so the piece always matches its white gap exactly.
                const double minPieceSize = 5.0;			// Keep a tiny excised deletion bar from shrinking below a visible size.
                const double excisedBarHeight = std::max(minPieceSize, ((fragHigher - fragLower) / originalSizeMbp) * remainingBarHeight);
                const double excisedBarX = excisedColumnCenterX - chromosomeWidth / 2.0;
                const std::vector<PaintedSegment> pieceSegments = buildPaintedSegments(*piece, humanGenome, masterHeader);

                drawPaintedFragment(cr, excisedBarX, excisedY, excisedBarHeight, chromosomeWidth, pieceSegments, roundTopCap, roundBottomCap, false);

                excisedY += excisedBarHeight + verticalGap;
            }
        }

        return;
    }


    if (isChromosomeRing && originalSizeMbp > 0.0)
    {
        const SDRdataRecord* ringRecord = nullptr;
        std::vector<const SDRdataRecord*> excisedRecords;

        for (const SDRdataRecord* record : records)
        {
            if (!record->linear)
            {
                ringRecord = record;
            }
            else
            {
                excisedRecords.push_back(record);
            }
        }

        if (ringRecord == nullptr)
        {
            return;
        }

        // Column 1: the two excised ends, stacked separately at the
        // original chromosome's slot position.
        const double barX = slotCenterX - totalWidth / 2.0;
        const double verticalGap = 6.0;
        double excisedY = posY;
        const double delTolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)

        for (const SDRdataRecord* record : excisedRecords)
        {
            const SDRfragment& fragment = record->fragments[0];
            const double fragLower = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
            const double fragHigher = std::max(fragment.oldStartPosition, fragment.oldEndPosition);
            const double pieceLengthMbp = fragHigher - fragLower;

            const bool roundTopCap = approxEqual(fragLower, 0.0, delTolerance);
            const bool roundBottomCap = approxEqual(fragHigher, originalSizeMbp, delTolerance);

            const double pieceBarHeight = computeSDRbarHeight(pieceLengthMbp, maxLengthMbp, maxRenderHeight);
            const std::vector<PaintedSegment> segments = buildPaintedSegments(*record, humanGenome, masterHeader);

            drawPaintedFragment(cr, barX, excisedY, pieceBarHeight, chromosomeWidth, segments, roundTopCap, roundBottomCap, false);

            excisedY += pieceBarHeight + verticalGap;
        }

        // Column 2: the ring itself, sized to its own length.
        double ringLengthMbp = 0.0;

        for (const SDRfragment& fragment : ringRecord->fragments)
        {
            ringLengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
        }

	const double remainingBarHeight = computeSDRbarHeight(originalSizeMbp, maxLengthMbp, maxRenderHeight);
        const double equivalentLinearHeight = (ringLengthMbp / originalSizeMbp) * remainingBarHeight;
        const double minPieceSize = 5.0;				// Matches the floor computeSlotWidth() already applied when it sized this slot's column.
        const double ringDiameter = std::max(minPieceSize, std::min(equivalentLinearHeight, chromosomeWidth * 3.0));

        const double ringColumnCenterX = barX + chromosomeWidth + circleColumnGap + column2Width / 2.0;
        const double ringCenterY = posY + ringDiameter / 2.0;

        drawSingleOriginCircularPiece(cr, *ringRecord, ringColumnCenterX, ringCenterY, ringDiameter, homeOldStrandID, humanGenome, masterHeader);

        return;
    }


    if (isDeletionTranslocationDonorShape(records, homeOldStrandID))
    {
        // The recombined molecule carries no gap, the deleted region simply isn't part of it, which is already correct. It draws
        // exactly like an ordinary translocation derivative: plain fragments, scaled to its own actual length, no white-gap or
        // full-original-length treatment. The excised piece gets the same flat-cap treatment as any other deletion's excised piece.
        const SDRdataRecord* recombinedRecord = nullptr;
        const SDRdataRecord* excisedRecord = nullptr;

	// Tell the two records apart by fragment origin: the one carrying at least one foreign
	// (non-home-strand) fragment is the recombined derivative; the one made up entirely of home-strand
	// fragments is the excised piece.
        for (const SDRdataRecord* record : records)
        {
            bool hasForeign = false;

	    // The deletion translocation donor must contain at least one foreign fragment
            for (const SDRfragment& fragment : record->fragments)
            {
                if (fragment.oldStrandID != homeOldStrandID)			// The origin of the fragment must be different than the location it was drawn in.
                {
                    hasForeign = true;
                    break;
                }
            }

            if (hasForeign)
            {
                recombinedRecord = record;
            }
            else
            {
                excisedRecord = record;
            }
        }

	// Deletion-translocation donor requires a recombined record containing original and foreign fragments, and an excised record.
	if (recombinedRecord == nullptr || excisedRecord == nullptr)
        {
            return;
        }

	// Determine the horizontal positions to draw the left and right homologs of the chromosome
        const double columnX1 = slotCenterX - totalWidth / 2.0;
        const double columnX2 = columnX1 + chromosomeWidth + stackGap;

	// Calculate change in chromosome new length after the translocation occurs
        double recombinedLengthMbp = 0.0;
	for (const SDRfragment& fragment : recombinedRecord->fragments)
        {
            recombinedLengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
        }

	// Convert bar height of a chromosome in Mbp to a height in pixels on the image
	const double recombinedBarHeight = computeSDRbarHeight(recombinedLengthMbp, maxLengthMbp, maxRenderHeight);
	const std::vector<PaintedSegment> recombinedSegments = buildPaintedSegments(*recombinedRecord, humanGenome, masterHeader);

	// Draw all the original and foreign fragments at a given chromosome slot
	drawPaintedFragment(cr, columnX1, posY, recombinedBarHeight, chromosomeWidth, recombinedSegments);

        const SDRfragment& excisedFragment = excisedRecord->fragments[0];
        const double fragLower = std::min(excisedFragment.oldStartPosition, excisedFragment.oldEndPosition);
        const double fragHigher = std::max(excisedFragment.oldStartPosition, excisedFragment.oldEndPosition);
        const double excisedLengthMbp = fragHigher - fragLower;

        const double delTolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)
	// Booleans to decide if the deleted segment was at a chromosome end, and if it should draw a rounded cap to indicate the deletion occurred at an end.
        const bool roundTopCap = approxEqual(fragLower, 0.0, delTolerance);
        const bool roundBottomCap = (originalSizeMbp > 0.0) ? approxEqual(fragHigher, originalSizeMbp, delTolerance) : false;

	// Calculate the height of the deletion in pixels and build the deleted fragment beside the chromosome it originated from
        const double excisedBarHeight = computeSDRbarHeight(excisedLengthMbp, maxLengthMbp, maxRenderHeight);
        const std::vector<PaintedSegment> excisedSegments = buildPaintedSegments(*excisedRecord, humanGenome, masterHeader);

	// Draw the deleted segment at the set location beside its original chromosome
        drawPaintedFragment(cr, columnX2, posY, excisedBarHeight, chromosomeWidth, excisedSegments, roundTopCap, roundBottomCap, false);

        return;

    }


    // Multiple records, none matching a named shape - shattered/collapsed-concatenation case (chromothripsis-style). One base
    // record (most fragments) drawn plainly, every other piece sized exactly to its own length (no floor) and stacked vertically,
    // wrapping into additional columns rather than growing past maxRenderHeight into the next karyogram row.
    if (records.size() > 1 || isSingleMixedRecordWithHomeGap(records, homeOldStrandID))
    {
        const SDRdataRecord* baseRecord = nullptr;
        std::vector<const SDRdataRecord*> otherRecords;

        for (const SDRdataRecord* record : records)
        {
            if (baseRecord == nullptr || record->fragments.size() > baseRecord->fragments.size())
            {
                if (baseRecord != nullptr)
                {
                    otherRecords.push_back(baseRecord);
                }

                baseRecord = record;
            }
            else
            {
                otherRecords.push_back(record);
            }
        }

        const double baseX = slotCenterX - totalWidth / 2.0;

        if (baseRecord != nullptr)
        {
	    bool baseIsSingleStrand = true;

            for (const SDRfragment& fragment : baseRecord->fragments)
            {
                if (fragment.oldStrandID != homeOldStrandID)
                {
                    baseIsSingleStrand = false;
                    break;
                }
            }

            std::vector<PaintedSegment> baseSegments;
            double baseLengthMbp = 0.0;

	    if (baseIsSingleStrand && baseRecord->fragments.size() >= 2)
            {
                // Pure single-strand case - reuse the same function plain multi-deletion already uses, position-sorted
                // and robust to file-listing order.
                baseSegments = buildDeletionRemainingSegments(*baseRecord, homeOldStrandID, humanGenome, masterHeader);

                // buildDeletionRemainingSegments scales by the chromosome's full original size, not just the fragment+gap sum.
                // Match that here so the bar height stays consistent.
                if (sizeIndex < masterHeader.intactChromosomeSizes.size())
                {
                    baseLengthMbp = masterHeader.intactChromosomeSizes[sizeIndex];
                }
            }

	    else
            {
                // Mixed-strand case (material exchanged between two chromosomes), fragment positions live in different
                // coordinate systems, so buildDeletionRemainingSegments' sorting assumption doesn't hold. Walks fragments in
                // file order instead, inserting white gaps only between immediately-adjacent HOME-strand fragments.
                // Home pieces moved within the chromosome: laid out like a single-strand record with
                // moved pieces, so deleted stretches sit where they belong and moved pieces get their arrows.
                baseSegments = buildMixedStrandSegmentsWithGaps(*baseRecord, homeOldStrandID, humanGenome, masterHeader, baseLengthMbp, homePiecesMoved(*baseRecord, homeOldStrandID));
            }

            const double baseHeight = computeSDRbarHeight(baseLengthMbp, maxLengthMbp, maxRenderHeight);

            drawPaintedFragment(cr, baseX, posY, baseHeight, chromosomeWidth, baseSegments);
        }

        const double verticalGap = 6.0;
        const double delTolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)
        std::vector<double> pieceHeights;
        computeExcisedColumnLayout(otherRecords, maxLengthMbp, maxRenderHeight, maxRenderHeight, verticalGap, pieceHeights);

        double columnX = baseX + chromosomeWidth + stackGap;
        double columnY = posY;


	for (std::size_t i = 0; i < otherRecords.size(); ++i)
        {
            const SDRdataRecord* record = otherRecords[i];
            const double pieceHeight = pieceHeights[i];

            if (i > 0)
            {
                const double previousHeight = pieceHeights[i - 1];

                if ((columnY - posY) + previousHeight + verticalGap + pieceHeight > maxRenderHeight)
                {
                    columnX += chromosomeWidth + stackGap;
                    columnY = posY;
                }
                else
                {
                    columnY += previousHeight + verticalGap;
                }
            }

            // A piece that is PURELY home-strand material represents a simple deleted/excised segment (nothing exchanged with
            // another chromosome), gets the same "no outline" treatment as a plain deletion's excised piece. A piece
            // carrying foreign material represents an exchange, not a deletion, so it keeps its outline.
            bool isPureHomeStrand = true;

            for (const SDRfragment& fragment : record->fragments)
            {
                if (fragment.oldStrandID != homeOldStrandID)
                {
                    isPureHomeStrand = false;
                    break;
                }
            }

            const bool drawOutline = !isPureHomeStrand;

            if (!record->linear)
            {
                // Circular (ecDNA-like) piece, draw as a circle, same diameter convention as the plain ecDNA branch (scaled
                // directly to fragment length, capped at 3x chromosomeWidth).
                const double diameter = std::min(pieceHeight, chromosomeWidth * 3.0);
                const double centerY = columnY + diameter / 2.0;

		// Correct because genuinely mixed-origin circular records are filtered out upstream via multiOriginECDNAnewStrandIDs before reaching this fallback.
                // A record that slipped through with foreign fragments would be
		// miscolored here, since nothing local to this line checks the record's actual fragment origins.
                // A circular piece that kept the centromere is a ring chromosome, not ecDNA.
                drawSingleOriginCircularPiece(cr, *record, columnX + chromosomeWidth / 2.0, centerY, diameter, homeOldStrandID, humanGenome, masterHeader);
            }
            else
            {
                bool roundTopCap = false;
                bool roundBottomCap = false;

                if (record->fragments.size() == 1 && originalSizeMbp > 0.0)
                {
                    const SDRfragment& fragment = record->fragments[0];
                    const double fragLower = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
                    const double fragHigher = std::max(fragment.oldStartPosition, fragment.oldEndPosition);

                    roundTopCap = approxEqual(fragLower, 0.0, delTolerance);
                    roundBottomCap = approxEqual(fragHigher, originalSizeMbp, delTolerance);
                }

                const std::vector<PaintedSegment> segments = buildPaintedSegments(*record, humanGenome, masterHeader);

                drawPaintedFragment(cr, columnX, columnY, pieceHeight, chromosomeWidth, segments, roundTopCap, roundBottomCap, drawOutline);
            }
        }

        return;
    }


    // ---------------------------------------------------------------------------
    // Everything else keeps the existing one-column-per-record horizontal layout.
    // ---------------------------------------------------------------------------

    // Determine where to draw the abberrant strands and how large to draw them
    double barX = slotCenterX - totalWidth / 2.0;                          			// Aberrant strand X position

    for (const SDRdataRecord* record : records)
    {
        double lengthMbp = 0.0;                                            			// Convert Mbp lengths to pixel sizes
        for (const SDRfragment& fragment : record->fragments)
        {
            lengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);	// Determine new rearranged strand length
        }

        double barHeight = computeSDRbarHeight(lengthMbp, maxLengthMbp, maxRenderHeight); 	// Determine the height of the aberrant strand
        const std::vector<PaintedSegment> segments = buildPaintedSegments(*record, humanGenome, masterHeader); // Store all the segments associated with a given data record to be built and drawn

        drawPaintedFragment(cr, barX, posY, barHeight, chromosomeWidth, segments);        	// Draw the chromosome with the rearranged fragments
        barX += chromosomeWidth + stackGap;                                			// Adjust the position of the aberrant strand fragment x position
    }

}












// Function to build the painted segment for one new strand after it has been rearranged in the SDR file
std::vector<PaintedSegment> Karyogram::buildPaintedSegments(
    const SDRdataRecord& record, 							// Access the data record for the fragment start and end locations
    bool humanGenome, 									// Check if human genome is passed to adjust centromeres on fragments
    const SDRmasterHeader& masterHeader)						// Access SDR master header for chromosome sizes
{

    std::vector<PaintedSegment> segments;						// Initializing segments vector to hold PaintedSegment objects

    // Determining new strand length by summing the lengths of the fragments that make it up
    double totalLengthMbp = 0.0;
    for (const SDRfragment& fragment : record.fragments)
    {
        totalLengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
    }
    if (totalLengthMbp <= 0)
    {
        return segments;
    }

    // Temporary variable to hold locations of the segment breaks
    double runningPositionMbp = 0.0;

    for (const SDRfragment& fragment : record.fragments)
    {
        const double fragmentLengthMbp = std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);	// Data field 3 fragment length is the difference in fragment start and end locations in Mbp

        PaintedSegment segment{};							// Instantiate the PaintedSegment object segment{}
        segment.startFraction = runningPositionMbp / totalLengthMbp;			// Find start location of fragment as a fraction of the total chromosome length
        segment.endFraction = (runningPositionMbp + fragmentLengthMbp) / totalLengthMbp;// Find end location of fragment as a fraction of the total chromosome length
        segment.color = getColorForOriginalStrand(fragment.oldStrandID, masterHeader);	// Obtain color of the original chromosome from which this fragment originated
        segment.isSecondHomolog = isSecondHomologStrand(fragment.oldStrandID, masterHeader);
        segment.isReversed = fragment.oldStartPosition > fragment.oldEndPosition;	// If the start location is > the end location of the fragment, an inversion occurred.
        segment.oldStrandID = fragment.oldStrandID;
        segment.oldStartPosition = fragment.oldStartPosition;
        segment.oldEndPosition = fragment.oldEndPosition;
        segment.hasCentromere = false;							// By default, segments/fragments do not contain centromeres, unless hasCentromere = true.

        if (fragment.hasCentromere && fragmentLengthMbp > 0.0)				// If fragment does have a centromere
        {
            double centromereStartBP = 0.0;						// Determine centromere range
            double centromereEndBP = 0.0;

            if (getCentromereForOriginalStrand(fragment.oldStrandID, humanGenome, masterHeader, centromereStartBP, centromereEndBP)) // If original strand has centromere
            {
                const double centromereStartMbp = centromereStartBP / 1000000.0;
                const double centromereEndMbp = centromereEndBP / 1000000.0;

		// Ensure the centromere range stays within the span of a given fragment
                const double fragMinMbp = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
                const double fragMaxMbp = std::max(fragment.oldStartPosition, fragment.oldEndPosition);
                // Clamp the centromere span to the portion actually contained within this fragment.
                const double overlapStartMbp = std::max(centromereStartMbp, fragMinMbp);
                const double overlapEndMbp = std::min(centromereEndMbp, fragMaxMbp);

                double localFractionStart;
                double localFractionEnd;

                if (overlapEndMbp > overlapStartMbp)
                {
                    if (!segment.isReversed)						// No balanced inversion
                    {
			// Centromere range on the local fragment as a fraction of the fragment length
                        localFractionStart = (overlapStartMbp - fragMinMbp) / fragmentLengthMbp;
                        localFractionEnd = (overlapEndMbp - fragMinMbp) / fragmentLengthMbp;
                    }
                    else
                    {
                        // Reversed fragment - genomic direction is flipped relative to the visual (top-to-bottom) direction of the drawn segment.
                        localFractionStart = 1.0 - (overlapEndMbp - fragMinMbp) / fragmentLengthMbp;
                        localFractionEnd = 1.0 - (overlapStartMbp - fragMinMbp) / fragmentLengthMbp;
                    }
                }
                else
                {
                    // hasCentromere is set but our lookup doesn't actually overlap this fragment's range (e.g. table/data mismatch). Fall back to a
                    // point at the fragment's midpoint rather than drawing nothing.
                    localFractionStart = 0.5;
                    localFractionEnd = 0.5;
                }

		// Adjust the centromere location relative to the fragment length and associate it to the entire mutated chromosome length
		// To be converted to pixels and drawn on the karyogram
                segment.hasCentromere = true;
                segment.centromereStartFraction = segment.startFraction + localFractionStart * (segment.endFraction - segment.startFraction);
                segment.centromereEndFraction = segment.startFraction + localFractionEnd * (segment.endFraction - segment.startFraction);
            }
        }

        segments.push_back(segment);					// Store this segment in the segments vector to be drawn in the karyogram
        runningPositionMbp += fragmentLengthMbp;			// Move on to the next fragment in a given data record
    }

    return segments;							// Return geometrical information of the mutated segments/fragments to be drawn on the karyogram.
}













// Are pieces of one strand listed out of their original order, beyond inversions? Pieces are grouped into blocks first: an inverted stretch (a reversed piece
// and the pieces that keep reading downward right after it, nothing else of the strand lying in between) is one block, so plain and nested inversions don't count.
// True when a block starts below one listed before it. Same block rule as the moved-piece arrows in drawPaintedFragment().
static bool hasPiecesMovedWithinStrand(const std::vector<SDRfragment>& fragments)
{
    const auto lowOf = [](const SDRfragment& fragment) { return std::min(fragment.oldStartPosition, fragment.oldEndPosition); };
    const auto highOf = [](const SDRfragment& fragment) { return std::max(fragment.oldStartPosition, fragment.oldEndPosition); };

    std::map<int, double> highestBlockLow;						// Per strand, the highest block start seen so far
    std::size_t i = 0;

    while (i < fragments.size())
    {
        const SDRfragment& first = fragments[i];
        double blockLow = lowOf(first);
        std::size_t last = i;

        if (first.oldStartPosition > first.oldEndPosition)
        {
            while (last + 1 < fragments.size())
            {
                const SDRfragment& previous = fragments[last];
                const SDRfragment& next = fragments[last + 1];

                if (next.oldStrandID != first.oldStrandID || highOf(next) > lowOf(previous) + svToleranceMbp)
                {
                    break;
                }

                // An inverted block ends with a reversed piece: a forward piece is part of it only when a reversed piece of the block follows it
                // (a piece flipped back inside the inversion), otherwise it was swapped with the inversion (see drawPaintedFragment()).
                if (next.oldStartPosition <= next.oldEndPosition)
                {
                    bool reversedFollows = false;
                    for (std::size_t k = last + 2; k < fragments.size() && fragments[k].oldStrandID == first.oldStrandID; ++k)
                    {
                        if (fragments[k].oldStartPosition > fragments[k].oldEndPosition)
                        {
                            reversedFollows = highOf(fragments[k]) <= lowOf(next) + svToleranceMbp;
                            break;
                        }
                    }

                    if (!reversedFollows)
                    {
                        break;
                    }
                }

                bool somethingBetween = false;
                for (std::size_t k = 0; k < fragments.size() && !somethingBetween; ++k)
                {
                    somethingBetween = k != last && k != last + 1 && fragments[k].oldStrandID == first.oldStrandID
                                       && highOf(fragments[k]) > highOf(next) + svToleranceMbp && lowOf(fragments[k]) < lowOf(previous) - svToleranceMbp;
                }

                if (somethingBetween)
                {
                    break;
                }

                ++last;
                blockLow = std::min(blockLow, lowOf(next));
            }
        }

        const auto seen = highestBlockLow.find(first.oldStrandID);
        if (seen != highestBlockLow.end() && blockLow < seen->second)
        {
            return true;
        }

        highestBlockLow[first.oldStrandID] = (seen == highestBlockLow.end()) ? blockLow : std::max(seen->second, blockLow);
        i = last + 1;
    }

    return false;
}



// Which pieces of a single-chromosome record moved, by the same rule drawPaintedFragment() uses for the arrows: pieces are grouped into
// blocks (an inverted block starts and ends with a reversed piece and reads downward with nothing of the strand in between; forward pieces trailing it were
// swapped with it), and the heaviest increasing run of blocks by length stays in place (an inverted block swapped with a trailing piece is preferred to stay).
static std::vector<bool> movedFragmentFlags(const std::vector<SDRfragment>& fragments)
{
    const auto lowOf = [](const SDRfragment& fragment) { return std::min(fragment.oldStartPosition, fragment.oldEndPosition); };
    const auto highOf = [](const SDRfragment& fragment) { return std::max(fragment.oldStartPosition, fragment.oldEndPosition); };
    const auto reversed = [](const SDRfragment& fragment) { return fragment.oldStartPosition > fragment.oldEndPosition; };

    struct Block { std::size_t first; std::size_t last; double low; double high; bool swapped; double key = 0.0; };
    std::vector<Block> blocks;

    for (std::size_t i = 0; i < fragments.size(); ++i)
    {
        if (!blocks.empty())
        {
            Block& previous = blocks.back();
            const SDRfragment& previousFragment = fragments[previous.last];
            bool somethingBetween = false;

            for (std::size_t k = 0; k < fragments.size() && !somethingBetween; ++k)
            {
                somethingBetween = k != i && k != previous.last && fragments[k].oldStrandID == fragments[i].oldStrandID
                                   && highOf(fragments[k]) > highOf(fragments[i]) + svToleranceMbp && lowOf(fragments[k]) < lowOf(previousFragment) - svToleranceMbp;
            }

            if (fragments[previous.first].oldStrandID == fragments[i].oldStrandID && reversed(fragments[previous.first])
                && highOf(fragments[i]) <= lowOf(previousFragment) + svToleranceMbp && !somethingBetween)
            {
                previous.last = i;
                previous.low = std::min(previous.low, lowOf(fragments[i]));
                previous.high = std::max(previous.high, highOf(fragments[i]));
                continue;
            }
        }

        blocks.push_back({i, i, lowOf(fragments[i]), highOf(fragments[i]), false});
    }

    for (std::size_t b = 0; b < blocks.size(); ++b)					// Split off forward pieces trailing an inverted block
    {
        Block& block = blocks[b];

        if (!reversed(fragments[block.first]) || reversed(fragments[block.last]))
        {
            continue;
        }

        std::size_t newLast = block.last;
        while (newLast > block.first && !reversed(fragments[newLast]))
        {
            --newLast;
        }

        std::vector<Block> trailing;
        for (std::size_t k = newLast + 1; k <= block.last; ++k)
        {
            trailing.push_back({k, k, lowOf(fragments[k]), highOf(fragments[k]), false});
        }

        block.last = newLast;
        block.swapped = true;
        block.low = lowOf(fragments[block.first]);
        block.high = highOf(fragments[block.first]);
        for (std::size_t k = block.first; k <= block.last; ++k)
        {
            block.low = std::min(block.low, lowOf(fragments[k]));
            block.high = std::max(block.high, highOf(fragments[k]));
        }

        blocks.insert(blocks.begin() + static_cast<long>(b) + 1, trailing.begin(), trailing.end());
        b += trailing.size();
    }

    for (Block& block : blocks)
    {
        block.key = block.low;

        // An inverted block whose own original place is a clean cut between two pieces listed next to each other (one ends where it started, the next
        // starts where it ended) was cut out from there: it is the moved one, not the forward piece it was swapped with.
        if (block.swapped)
        {
            for (std::size_t k = 0; k + 1 < fragments.size(); ++k)
            {
                const SDRfragment& above = fragments[k];
                const SDRfragment& below = fragments[k + 1];
                const int strandID = fragments[block.first].oldStrandID;

                if (k >= block.first && k <= block.last)
                {
                    continue;
                }

                if (above.oldStrandID == strandID && below.oldStrandID == strandID
                    && ((approxEqual(highOf(above), block.low, svToleranceMbp) && approxEqual(lowOf(below), block.high, svToleranceMbp))
                        || (approxEqual(lowOf(above), block.high, svToleranceMbp) && approxEqual(highOf(below), block.low, svToleranceMbp))))
                {
                    block.swapped = false;
                    break;
                }
            }
        }
    }

    // An inverted block listed right before a forward piece of its strand that originally lay below it, with other pieces of the strand in between (so the
    // two were not joined above), also swapped places with that piece: the inversion took the forward piece's place and keeps it, the forward piece
    // is the moved one. The inversion is ordered by the place it took, so the pieces after it stay in place.
    for (std::size_t b = 0; b + 1 < blocks.size(); ++b)
    {
        Block& block = blocks[b];
        const Block& next = blocks[b + 1];
        const SDRfragment& nextFragment = fragments[next.first];

        if (block.swapped || !reversed(fragments[block.first]) || !reversed(fragments[block.last]) || next.first != next.last || reversed(nextFragment)
            || nextFragment.oldStrandID != fragments[block.first].oldStrandID || highOf(nextFragment) > block.low + svToleranceMbp)
        {
            continue;
        }

        // Not when the inversion's own place is a clean cut: two pieces listed next to each other that flanked it originally (it was moved out from
        // between them, so the inversion is the moved piece).
        bool leftCleanCut = false;
        for (std::size_t k = 0; k + 1 < fragments.size() && !leftCleanCut; ++k)
        {
            leftCleanCut = fragments[k].oldStrandID == nextFragment.oldStrandID && fragments[k + 1].oldStrandID == nextFragment.oldStrandID
                           && !reversed(fragments[k]) && !reversed(fragments[k + 1])
                           && approxEqual(highOf(fragments[k]), block.low, svToleranceMbp) && approxEqual(lowOf(fragments[k + 1]), block.high, svToleranceMbp);
        }

        // And only when the strand reads on upward from the forward piece, up to the inversion's own upper end (deleted stretches may lie in between): the
        // inversion and that run simply traded places.
        bool readsOnToInversion = false;
        for (std::size_t k = next.first; k + 1 < fragments.size(); ++k)
        {
            const SDRfragment& current = fragments[k];
            const SDRfragment& following = fragments[k + 1];

            if (following.oldStrandID != current.oldStrandID || reversed(following) || lowOf(following) < highOf(current) - svToleranceMbp)
            {
                break;
            }

            if (approxEqual(lowOf(following), block.high, svToleranceMbp))
            {
                readsOnToInversion = true;
                break;
            }

            if (highOf(following) > block.low + svToleranceMbp)
            {
                break;
            }
        }

        if (leftCleanCut || !readsOnToInversion)
        {
            continue;
        }

        block.swapped = true;
        block.key = next.low;
    }

    const std::size_t n = blocks.size();
    const double oneBaseMbp = 0.000001;
    const auto weight = [&](const Block& block) { return block.high - block.low + oneBaseMbp + (block.swapped ? 1000.0 : 0.0); };
    std::vector<double> best(n, 0.0);
    std::vector<int> previousInRun(n, -1);

    for (std::size_t a = 0; a < n; ++a)
    {
        best[a] = weight(blocks[a]);

        for (std::size_t c = 0; c < a; ++c)
        {
            if (fragments[blocks[c].first].oldStrandID == fragments[blocks[a].first].oldStrandID && blocks[c].key < blocks[a].key && best[c] + weight(blocks[a]) > best[a])
            {
                best[a] = best[c] + weight(blocks[a]);
                previousInRun[a] = static_cast<int>(c);
            }
        }
    }

    std::vector<bool> moved(fragments.size(), true);
    int runEnd = n ? static_cast<int>(std::max_element(best.begin(), best.end()) - best.begin()) : -1;

    while (runEnd != -1)
    {
        for (std::size_t k = blocks[runEnd].first; k <= blocks[runEnd].last; ++k)
        {
            moved[k] = false;
        }
        runEnd = previousInRun[runEnd];
    }

    return moved;
}


// Builds painted segments for the "remaining" (2-or-more-fragment) half of a long deletion, scaled to the ORIGINAL chromosome's full length
// rather than this record's own reduced total - so the deleted region shows up as a distinct white gap rather than disappearing.
std::vector<PaintedSegment> Karyogram::buildDeletionRemainingSegments(
    const SDRdataRecord& record,
    int homeOldStrandID,
    bool humanGenome,
    const SDRmasterHeader& masterHeader)
{
    std::vector<PaintedSegment> segments;


    // This function must only ever be called on a record already confirmed by isDeletionShape/isECDNAshape/isDeletionInversionShape
    // to have EVERY fragment from the home strand. It assumes all fragment positions share one coordinate system and sorts them
    // directly. A record mixing in foreign translocated fragments must never be routed here.

    if (record.fragments.empty())
    {
        return segments;
    }

    const std::size_t sizeIndex = static_cast<std::size_t>(homeOldStrandID);

    if (sizeIndex >= masterHeader.intactChromosomeSizes.size())
    {
        return segments;
    }

    const double originalSizeMbp = masterHeader.intactChromosomeSizes[sizeIndex];

    if (originalSizeMbp <= 0.0)
    {
        return segments;
    }

    // Pieces moved within the chromosome would be put back in their original order by the position sort below, hiding the move (only a border line
    // was left). Such records keep their file order instead, with each deleted stretch still drawn as a white gap next to the piece it bordered.
    if (hasPiecesMovedWithinStrand(record.fragments))
    {
        double totalLengthMbp = 0.0;
        return buildMixedStrandSegmentsWithGaps(record, homeOldStrandID, humanGenome, masterHeader, totalLengthMbp, true);
    }


    std::vector<SDRfragment> fragments = record.fragments;


    std::sort(fragments.begin(), fragments.end(), [](const SDRfragment& a, const SDRfragment& b)
    {
        return std::min(a.oldStartPosition, a.oldEndPosition) < std::min(b.oldStartPosition, b.oldEndPosition);
    });



    const RGB homeColor = getColorForOriginalStrand(homeOldStrandID, masterHeader);


    auto addFragmentSegment = [&](const SDRfragment& fragment)
    {
        const double lowerPos = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
        const double higherPos = std::max(fragment.oldStartPosition, fragment.oldEndPosition);

        PaintedSegment segment{};
        segment.startFraction = lowerPos / originalSizeMbp;
        segment.endFraction = higherPos / originalSizeMbp;
        segment.color = homeColor;
        segment.isSecondHomolog = isSecondHomologStrand(homeOldStrandID, masterHeader);
        segment.isReversed = fragment.oldStartPosition > fragment.oldEndPosition;
        segment.hasCentromere = false;

        if (fragment.hasCentromere)
        {
            double centromereStartBP = 0.0;
            double centromereEndBP = 0.0;

            if (getCentromereForOriginalStrand(homeOldStrandID, humanGenome, masterHeader, centromereStartBP, centromereEndBP))
            {
                const double centStartMbp = centromereStartBP / 1000000.0;
                const double centEndMbp = centromereEndBP / 1000000.0;

                const double lowerOverlap = std::max(lowerPos, centStartMbp);
                const double higherOverlap = std::min(higherPos, centEndMbp);

                if (higherOverlap > lowerOverlap)
                {
                    segment.hasCentromere = true;
                    segment.centromereStartFraction = lowerOverlap / originalSizeMbp;
                    segment.centromereEndFraction = higherOverlap / originalSizeMbp;
                }
            }
        }

        segments.push_back(segment);
    };


    const double delTolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)


    // Leading white gap if the first fragment doesn't start at position 0, a terminal deletion removed the chromosome's start.
    const double firstFragmentStart = std::min(fragments.front().oldStartPosition, fragments.front().oldEndPosition);

    if (firstFragmentStart > 0.0 && !approxEqual(firstFragmentStart, 0.0, delTolerance))
    {
        PaintedSegment gapSegment{};
        gapSegment.startFraction = 0.0;
        gapSegment.endFraction = firstFragmentStart / originalSizeMbp;
        gapSegment.color = RGB{1.0, 1.0, 1.0};
        gapSegment.isReversed = false;
        gapSegment.hasCentromere = false;
        segments.push_back(gapSegment);
    }


    // Only a balanced inversion (not reaching either chromosome end) is drawn as one chevron across a deletion inside it; a terminal one keeps a chevron per piece.
    const std::vector<SDRfragment> joinedInversions = joinInvertedPiecesAcrossGaps(record.fragments, delTolerance);
    auto insideTerminalInversion = [&](double gapStart, double gapEnd)
    {
        for (const SDRfragment& fragment : joinedInversions)
        {
            const double low = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
            const double high = std::max(fragment.oldStartPosition, fragment.oldEndPosition);

            if (fragment.oldStrandID == homeOldStrandID && low <= gapStart && gapEnd <= high)
            {
                return approxEqual(low, 0.0, delTolerance) || approxEqual(high, originalSizeMbp, delTolerance);
            }
        }
        return false;
    };


    // Draw central white gaps for central deleted segments
    for (std::size_t i = 0; i < fragments.size(); i++)
    {
        addFragmentSegment(fragments[i]);

        if (i + 1 < fragments.size())
        {
            const double gapStart = std::max(fragments[i].oldStartPosition, fragments[i].oldEndPosition);
            const double gapEnd = std::min(fragments[i + 1].oldStartPosition, fragments[i + 1].oldEndPosition);

            // Same tolerance as the leading/trailing gaps and the detectors (svToleranceMbp): pieces listed on neighbouring bases (e.g. 120.000000 then
            // 120.000001) touch, so no white sliver is drawn between them.
            if (gapEnd > gapStart && !approxEqual(gapStart, gapEnd, delTolerance))
            {
                PaintedSegment gapSegment{};
                gapSegment.startFraction = gapStart / originalSizeMbp;
                gapSegment.endFraction = gapEnd / originalSizeMbp;
                gapSegment.color = RGB{1.0, 1.0, 1.0};
                gapSegment.isReversed = false;
                gapSegment.hasCentromere = false;
                gapSegment.insideInversion = recordHasGapInsideInvertedSegment(record, homeOldStrandID, gapStart, gapEnd, delTolerance)
                                             && !insideTerminalInversion(gapStart, gapEnd);
                segments.push_back(gapSegment);
            }
        }
    }


    // Trailing white gap if the last fragment doesn't end at the chromosome's true full length, a terminal deletion removed the end.
    const double lastFragmentEnd = std::max(fragments.back().oldStartPosition, fragments.back().oldEndPosition);

    if (originalSizeMbp > lastFragmentEnd && !approxEqual(lastFragmentEnd, originalSizeMbp, delTolerance))
    {
        PaintedSegment gapSegment{};
        gapSegment.startFraction = lastFragmentEnd / originalSizeMbp;
        gapSegment.endFraction = 1.0;
        gapSegment.color = RGB{1.0, 1.0, 1.0};
        gapSegment.isReversed = false;
        gapSegment.hasCentromere = false;
        segments.push_back(gapSegment);
    }

    return segments;
}





















// Builds painted segments for a base record that may mix home-strand fragments with foreign (exchanged) ones. Walks fragments in FILE
// ORDER (not sorted by position, since foreign fragments live in a different chromosome's coordinate system entirely). A white gap is
// inserted ONLY between two immediately-adjacent fragments of the SAME strand (home or foreign) that don't touch, real lost material. No gap is ever inserted
// around a foreign fragment, since that material is still physically present on the molecule, just relocated from elsewhere. Also
// allows the pure single-strand case (every fragment home) with the same logic, since every adjacent pair is then a home-home check.
std::vector<PaintedSegment> Karyogram::buildMixedStrandSegmentsWithGaps(
    const SDRdataRecord& record,
    int homeOldStrandID,
    bool humanGenome,
    const SDRmasterHeader& masterHeader,
    double& outTotalLengthMbp,								// Total length of chromosome plus the new fragment (can be larger than the original chromosome length)
    bool movedPiecesLayout)								// Single-strand record with pieces moved within the chromosome (see buildDeletionRemainingSegments)
{
    std::vector<PaintedSegment> segments;
    outTotalLengthMbp = 0.0;

    if (record.fragments.empty())
    {
        return segments;
    }

    const double tolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)

    struct PlannedPiece
    {
        bool isGap;
        const SDRfragment* fragment; 							// valid if !isGap
        double lengthMbp;
        bool insideInversion = false;							// Gaps only, see PaintedSegment::insideInversion
    };

    std::vector<PlannedPiece> plan;
    double totalLengthMbp = 0.0;

    if (movedPiecesLayout)
    {
        // Pieces in file order. Every stretch of the chromosome found in no piece (deleted) is a white gap: the chromosome-end ones at the top/bottom,
        // an inner one next to a piece it borders (see below), so material that merely moved
        // elsewhere on the chromosome is never drawn as a gap. Scaled to the original chromosome length like the sorted deletion layout.
        const std::size_t sizeIndex = static_cast<std::size_t>(homeOldStrandID);
        const double originalSizeMbp = (sizeIndex < masterHeader.intactChromosomeSizes.size()) ? masterHeader.intactChromosomeSizes[sizeIndex] : 0.0;

        std::vector<std::pair<double, double>> covered;
        for (const SDRfragment& fragment : record.fragments)
        {
            if (fragment.oldStrandID != homeOldStrandID)				// A mixed record's foreign pieces (exchanged material) are not part of the home chromosome
            {
                continue;
            }
            covered.emplace_back(std::min(fragment.oldStartPosition, fragment.oldEndPosition), std::max(fragment.oldStartPosition, fragment.oldEndPosition));
        }
        std::sort(covered.begin(), covered.end());

        std::vector<std::pair<double, double>> deleted;					// Stretches no piece covers, from 0 to the original length
        double reached = 0.0;
        for (const auto& [low, high] : covered)
        {
            if (low > reached && !approxEqual(low, reached, tolerance))
            {
                deleted.emplace_back(reached, low);
            }
            reached = std::max(reached, high);
        }
        if (originalSizeMbp > reached && !approxEqual(originalSizeMbp, reached, tolerance))
        {
            deleted.emplace_back(reached, originalSizeMbp);
        }

        std::vector<bool> moved = movedFragmentFlags(record.fragments);
        for (std::size_t i = 0; i < record.fragments.size(); ++i)
        {
            moved[i] = moved[i] && record.fragments[i].oldStrandID == homeOldStrandID;
        }
        std::vector<std::vector<std::pair<double, double>>> gapsBefore(record.fragments.size());
        std::vector<std::vector<std::pair<double, double>>> gapsAfter(record.fragments.size());
        std::vector<std::pair<double, double>> leadingGaps;
        std::vector<std::pair<double, double>> trailingGaps;

        for (const auto& gap : deleted)
        {
            // An end of the chromosome exchanged for another chromosome's material (a foreign piece at that end of the record) was not deleted: no gap.
            if (approxEqual(gap.first, 0.0, tolerance))
            {
                if (record.fragments.front().oldStrandID == homeOldStrandID)
                {
                    leadingGaps.push_back(gap);
                }
                continue;
            }
            if (approxEqual(gap.second, originalSizeMbp, tolerance))
            {
                if (record.fragments.back().oldStrandID == homeOldStrandID)
                {
                    trailingGaps.push_back(gap);
                }
                continue;
            }

            // Deleted stretches stay next to a piece that did not move: only a moved piece leaves its place. Among pieces in place, the
            // gap goes next to the piece whose upper end borders it (after a forward piece, before a reversed one, e.g. above an inversion), else next to the
            // piece whose lower end borders its upper end (before a forward piece, after a reversed one), but never at the spot a moved piece was cut out
            // of, which stays plain (there the white gap follows the moved piece instead). Moved pieces are used only when nothing else fits.
            bool placed = false;

            const auto lowOf = [](const SDRfragment& fragment) { return std::min(fragment.oldStartPosition, fragment.oldEndPosition); };
            const auto highOf = [](const SDRfragment& fragment) { return std::max(fragment.oldStartPosition, fragment.oldEndPosition); };

            // Boundary k (between pieces k-1 and k) is a cut-out spot when those two pieces are the original neighbours of a moved piece.
            const auto isCutSite = [&](std::size_t k)
            {
                if (k == 0 || k >= record.fragments.size())
                {
                    return false;
                }

                const SDRfragment& above = record.fragments[k - 1];
                const SDRfragment& below = record.fragments[k];

                for (std::size_t m = 0; m < record.fragments.size(); ++m)
                {
                    if (!moved[m] || m == k - 1 || m == k)
                    {
                        continue;
                    }

                    const SDRfragment& piece = record.fragments[m];
                    const bool aboveLowSide = approxEqual(highOf(above), lowOf(piece), tolerance) && approxEqual(lowOf(below), highOf(piece), tolerance);
                    const bool aboveHighSide = approxEqual(lowOf(above), highOf(piece), tolerance) && approxEqual(highOf(below), lowOf(piece), tolerance);

                    if (aboveLowSide || aboveHighSide)
                    {
                        return true;
                    }
                }

                return false;
            };

            for (int pass = 0; pass < 2 && !placed; ++pass)				// pass 0: pieces in place only, pass 1: any piece
            {
                for (std::size_t i = 0; i < record.fragments.size() && !placed; ++i)
                {
                    const SDRfragment& fragment = record.fragments[i];
                    const bool forward = fragment.oldStartPosition <= fragment.oldEndPosition;

                    if ((pass == 1 || (!moved[i] && !isCutSite(forward ? i + 1 : i))) && approxEqual(highOf(fragment), gap.first, tolerance))
                    {
                        (forward ? gapsAfter[i] : gapsBefore[i]).push_back(gap);
                        placed = true;
                    }
                }

                for (std::size_t i = 0; i < record.fragments.size() && !placed; ++i)
                {
                    const SDRfragment& fragment = record.fragments[i];
                    const bool forward = fragment.oldStartPosition <= fragment.oldEndPosition;

                    if ((pass == 1 || (!moved[i] && !isCutSite(forward ? i : i + 1))) && approxEqual(lowOf(fragment), gap.second, tolerance))
                    {
                        (forward ? gapsBefore[i] : gapsAfter[i]).push_back(gap);
                        placed = true;
                    }
                }
            }
        }

        const auto addGap = [&](const std::pair<double, double>& gap)
        {
            PlannedPiece gapPiece{true, nullptr, gap.second - gap.first};
            gapPiece.insideInversion = recordHasGapInsideInvertedSegment(record, homeOldStrandID, gap.first, gap.second, tolerance);
            plan.push_back(gapPiece);
            totalLengthMbp += gapPiece.lengthMbp;
        };

        for (const auto& gap : leadingGaps) addGap(gap);

        for (std::size_t i = 0; i < record.fragments.size(); ++i)
        {
            for (const auto& gap : gapsBefore[i]) addGap(gap);

            const double fragLengthMbp = std::fabs(record.fragments[i].oldEndPosition - record.fragments[i].oldStartPosition);
            plan.push_back({false, &record.fragments[i], fragLengthMbp});
            totalLengthMbp += fragLengthMbp;

            for (const auto& gap : gapsAfter[i]) addGap(gap);
        }

        for (const auto& gap : trailingGaps) addGap(gap);
    }

    for (std::size_t i = 0; i < record.fragments.size() && !movedPiecesLayout; ++i)
    {
        const SDRfragment& fragment = record.fragments[i];
        const double fragLengthMbp = std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);

        if (i > 0)
        {
            const SDRfragment& previous = record.fragments[i - 1];

            if (fragment.oldStrandID == previous.oldStrandID && (fragment.oldStrandID == homeOldStrandID || (previous.oldStartPosition <= previous.oldEndPosition && fragment.oldStartPosition <= fragment.oldEndPosition && previous.oldEndPosition <= fragment.oldStartPosition)))
            {
                const double prevLo = std::min(previous.oldStartPosition, previous.oldEndPosition);
                const double prevHi = std::max(previous.oldStartPosition, previous.oldEndPosition);
                const double curLo = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
                const double curHi = std::max(fragment.oldStartPosition, fragment.oldEndPosition);

                double gapStart = 0.0;
                double gapEnd = 0.0;
                bool hasGap = false;

                if (prevHi <= curLo && !approxEqual(prevHi, curLo, tolerance))
                {
                    gapStart = prevHi;
                    gapEnd = curLo;
                    hasGap = true;
                }
                else if (curHi <= prevLo && !approxEqual(curHi, prevLo, tolerance))
                {
                    gapStart = curHi;
                    gapEnd = prevLo;
                    hasGap = true;
                }

                if (hasGap)
                {
                    const double gapLengthMbp = gapEnd - gapStart;
                    plan.push_back({true, nullptr, gapLengthMbp});
                    totalLengthMbp += gapLengthMbp;
                }
            }
        }

        plan.push_back({false, &fragment, fragLengthMbp});
        totalLengthMbp += fragLengthMbp;
    }

    if (totalLengthMbp <= 0.0)
    {
        return segments;
    }

    outTotalLengthMbp = totalLengthMbp;

    double runningPositionMbp = 0.0;

    for (const PlannedPiece& piece : plan)
    {
        if (piece.isGap)
        {
            PaintedSegment gapSegment{};
            gapSegment.startFraction = runningPositionMbp / totalLengthMbp;
            gapSegment.endFraction = (runningPositionMbp + piece.lengthMbp) / totalLengthMbp;
            gapSegment.color = RGB{1.0, 1.0, 1.0};
            gapSegment.isReversed = false;
            gapSegment.hasCentromere = false;
            gapSegment.insideInversion = piece.insideInversion;
            segments.push_back(gapSegment);

            runningPositionMbp += piece.lengthMbp;
            continue;
        }

        const SDRfragment& fragment = *piece.fragment;
        const double fragmentLengthMbp = piece.lengthMbp;

        PaintedSegment segment{};
        segment.startFraction = runningPositionMbp / totalLengthMbp;
        segment.endFraction = (runningPositionMbp + fragmentLengthMbp) / totalLengthMbp;
        segment.color = getColorForOriginalStrand(fragment.oldStrandID, masterHeader);
        segment.isSecondHomolog = isSecondHomologStrand(fragment.oldStrandID, masterHeader);
        segment.isReversed = fragment.oldStartPosition > fragment.oldEndPosition;
        segment.oldStrandID = fragment.oldStrandID;
        segment.oldStartPosition = fragment.oldStartPosition;
        segment.oldEndPosition = fragment.oldEndPosition;
        segment.hasCentromere = false;

        if (fragment.hasCentromere && fragmentLengthMbp > 0.0)
        {
            double centromereStartBP = 0.0;
            double centromereEndBP = 0.0;

            if (getCentromereForOriginalStrand(fragment.oldStrandID, humanGenome, masterHeader, centromereStartBP, centromereEndBP))
            {
                const double centromereStartMbp = centromereStartBP / 1000000.0;
                const double centromereEndMbp = centromereEndBP / 1000000.0;

                const double fragMinMbp = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
                const double fragMaxMbp = std::max(fragment.oldStartPosition, fragment.oldEndPosition);
                const double overlapStartMbp = std::max(centromereStartMbp, fragMinMbp);
                const double overlapEndMbp = std::min(centromereEndMbp, fragMaxMbp);

		// Clamp the centromere span to the portion actually contained within this fragment.
                if (overlapEndMbp > overlapStartMbp)
                {
                    double localFractionStart;
                    double localFractionEnd;

                    if (!segment.isReversed)
                    {
			// Centromere range on the local fragment as a fraction of the fragment length
                        localFractionStart = (overlapStartMbp - fragMinMbp) / fragmentLengthMbp;
                        localFractionEnd = (overlapEndMbp - fragMinMbp) / fragmentLengthMbp;
                    }
                    else
                    {
			// Reversed fragment, genomic direction is flipped relative to the visual (top-to-bottom) direction of the drawn segment.
                        localFractionStart = 1.0 - (overlapEndMbp - fragMinMbp) / fragmentLengthMbp;
                        localFractionEnd = 1.0 - (overlapStartMbp - fragMinMbp) / fragmentLengthMbp;
                    }

                    segment.hasCentromere = true;
                    segment.centromereStartFraction = segment.startFraction + localFractionStart * (segment.endFraction - segment.startFraction);
                    segment.centromereEndFraction = segment.startFraction + localFractionEnd * (segment.endFraction - segment.startFraction);
                }
            }
        }

        segments.push_back(segment);
        runningPositionMbp += fragmentLengthMbp;
    }

    return segments;
}













void Karyogram::drawMovedSegmentArrow(					// Red hooked arrow of a piece moved within its own chromosome, drawn left of the chromosome edge
    cairo_t* cr,
    double edgeX,								// Left edge of the chromosome in pixels
    double fromY,								// Middle of the piece the moved piece used to sit next to
    double toY,									// Middle of the moved piece, where the head points
    double bulge,								// How far left of the edge the curve reaches
    double lineWidth,
    double headLength,
    double headHalfWidth)
{
    // One rounded C curve that leaves and enters the chromosome horizontally, so the head (drawn as part of the line, the line ends inside it) sits on the curve.
    const double tipX = edgeX - 1.0;
    const double headBaseX = tipX - headLength;

    cairo_save(cr);
    cairo_set_source_rgb(cr, 0.75, 0.22, 0.17);

    cairo_new_path(cr);
    cairo_arc(cr, tipX, fromY, lineWidth * 1.4, 0.0, 2.0 * M_PI);			// Small dot where the arrow starts
    cairo_fill(cr);

    cairo_set_line_width(cr, lineWidth);
    cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
    cairo_move_to(cr, tipX, fromY);
    cairo_curve_to(cr, edgeX - bulge, fromY, edgeX - bulge, toY, headBaseX + 1.0, toY);
    cairo_stroke(cr);

    cairo_move_to(cr, tipX, toY);
    cairo_line_to(cr, headBaseX, toY - headHalfWidth);
    cairo_line_to(cr, headBaseX + 1.5, toY);
    cairo_line_to(cr, headBaseX, toY + headHalfWidth);
    cairo_close_path(cr);
    cairo_fill(cr);

    cairo_restore(cr);
}



void Karyogram::drawInversionChevron(					// Function to draw a marker to denote a balanced inversion occurred, spanning the inverted fragment
    cairo_t* cr,
    double centerX, 							// X position of the chevron's apex (the chromosome's horizontal center) in pixels
    double topY, 							// Top of the inverted fragment in pixels, where the apex sits
    double bottomY, 							// Bottom of the inverted fragment in pixels, where the two arms end
    double width)							// Width of the chromosome in pixels, the arms reach its left and right edges
{
    // Black chevron "^": apex at the top-centre of the segment, arms to the bottom-left and bottom-right corners. Its height is the segment's own height, so the
    // smaller the inversion the flatter the arrow gets; a segment of ~0 height collapses to one horizontal line the width of the chromosome.
    const double lineWidth = 1.5;
    const double halfWidth = width / 2.0;

    cairo_save(cr);

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_set_line_width(cr, lineWidth);
    cairo_set_line_cap(cr, CAIRO_LINE_CAP_BUTT);
    cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);			// Round join keeps the apex tidy as the arrow flattens, no long miter spike

    cairo_new_path(cr);
    cairo_move_to(cr, centerX - halfWidth, bottomY);
    cairo_line_to(cr, centerX, topY);
    cairo_line_to(cr, centerX + halfWidth, bottomY);
    cairo_stroke(cr);

    cairo_restore(cr);
}












// Draws an ecDNA fragment as a filled circle rather than a bar visually distinct from every other (linear) piece on the karyogram.
void Karyogram::drawCircularFragment(
    cairo_t* cr, 
    double centerX, 
    double centerY, 
    double diameter, 
    RGB color)
{
    const double radius = diameter / 2.0;

    cairo_new_path(cr);
    cairo_arc(cr, centerX, centerY, radius, 0.0, 2.0 * M_PI);

    cairo_set_source_rgb(cr, color.r, color.g, color.b);
    cairo_fill_preserve(cr);

    drawCircularFragmentInnerRing(cr, centerX, centerY, diameter);

}










// Draws one wedge per {color, fraction} entry in order, each spanning fraction * 360 degrees of the circle, so the assembled circle looks like a pie chart of each contributing strand's own
// color, proportional to how much of the fragment's total length that strand contributed. The wedge outlines against each other are always drawn (so the split between colors is visible even
// when they're similar).
void Karyogram::drawPieChartFragment(
    cairo_t* cr,
    double centerX,
    double centerY,
    double diameter,
    const std::vector<std::pair<RGB, double>>& wedges)
{
    const double radius = diameter / 2.0;

    // A single-wedge (single-origin) circle needs no pie slicing at all, draw it as a plain filled circle, same as drawCircularFragment(), to avoid a stray radius line across it.
    if (wedges.size() <= 1)
    {
        const RGB color = wedges.empty() ? RGB{0.5, 0.5, 0.5} : wedges.front().first;
        drawCircularFragment(cr, centerX, centerY, diameter, color);
        return;
    }

    double totalFraction = 0.0;

    for (const auto& wedge : wedges)
    {
        totalFraction += wedge.second;
    }

    if (totalFraction <= 0.0)
    {
        return;
    }

    double startAngle = -M_PI / 2.0; 							// Start at the top of the circle (12 o'clock), like a conventional pie chart.

    for (const auto& [color, fraction] : wedges)
    {
        const double sweepAngle = (fraction / totalFraction) * 2.0 * M_PI;
        const double endAngle = startAngle + sweepAngle;

        cairo_new_path(cr);
        cairo_move_to(cr, centerX, centerY);
        cairo_arc(cr, centerX, centerY, radius, startAngle, endAngle);
        cairo_close_path(cr);

        cairo_set_source_rgb(cr, color.r, color.g, color.b);
        cairo_fill_preserve(cr);

        // Thin line between wedges so adjacent similar colors stay visually distinguishable.
        cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
        cairo_set_line_width(cr, 0.75);
        cairo_stroke(cr);

        startAngle = endAngle;
    }

    drawCircularFragmentInnerRing(cr, centerX, centerY, diameter);

}











// Draws the small white, black-outlined inner circle that creates the "ring" appearance for ecDNA fragments - half the diameter of the
// outer circle, centered on top of it.
void Karyogram::drawCircularFragmentInnerRing(cairo_t* cr, double centerX, double centerY, double outerDiameter)
{
    const double innerRadius = outerDiameter / 4.0; 				// Half the outer DIAMETER = half the outer RADIUS.

    cairo_new_path(cr);
    cairo_arc(cr, centerX, centerY, innerRadius, 0.0, 2.0 * M_PI);

    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_fill(cr);
}












// Draws a circular piece made of one chromosome's material. If any fragment kept a centromere (the file's centromere flag, the same
// test detectChromosomeRings() uses), it is a ring chromosome and gets the ring shape with its centromere mark, dashed for a second
// homolog. Otherwise it is acentric ecDNA and gets the plain ecDNA circle. Both shapes have the same diameter, so slot widths are unchanged.
void Karyogram::drawSingleOriginCircularPiece(
    cairo_t* cr,
    const SDRdataRecord& record,
    double centerX,
    double centerY,
    double diameter,
    int homeOldStrandID,							// The slot's chromosome, which sets the color and the dashed outline
    bool humanGenome,								// Human centromere table or generic centromere at 50% of the chromosome length
    const SDRmasterHeader& masterHeader)
{
    const RGB color = getColorForOriginalStrand(homeOldStrandID, masterHeader);

    if (currentDetailSlot >= 0)
    {
        recordDetailCircle(record, centerX, centerY, diameter);
    }

    double ringLengthMbp = 0.0;
    const SDRfragment* centromereFragment = nullptr;
    double centromereFragmentOffsetMbp = 0.0;					// Where the centromere-bearing fragment starts along the ring

    for (const SDRfragment& fragment : record.fragments)
    {
        if (fragment.hasCentromere && centromereFragment == nullptr)
        {
            centromereFragment = &fragment;
            centromereFragmentOffsetMbp = ringLengthMbp;
        }

        ringLengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
    }

    if (centromereFragment == nullptr)
    {
        drawCircularFragment(cr, centerX, centerY, diameter, color);
        return;
    }

    // The fractions only size the centromere mark, which drawRingFragment() always places at 12 o'clock.
    double centromereStartFraction = 0.0;
    double centromereEndFraction = 0.0;
    double centromereStartBP = 0.0;
    double centromereEndBP = 0.0;

    if (ringLengthMbp > 0.0 && getCentromereForOriginalStrand(centromereFragment->oldStrandID, humanGenome, masterHeader, centromereStartBP, centromereEndBP))
    {
        const double fragLow = std::min(centromereFragment->oldStartPosition, centromereFragment->oldEndPosition);

        centromereStartFraction = (centromereFragmentOffsetMbp + centromereStartBP / 1000000.0 - fragLow) / ringLengthMbp;
        centromereEndFraction = (centromereFragmentOffsetMbp + centromereEndBP / 1000000.0 - fragLow) / ringLengthMbp;
    }

    // A ring formed from a second homolog gets dashed outlines, like that homolog's linear chromosome.
    drawRingFragment(cr, centerX, centerY, diameter, color, centromereStartFraction, centromereEndFraction,
                     isSecondHomologStrand(homeOldStrandID, masterHeader));
}



// Draws a chromosome ring: a colored circle (same fill and inner white
// ring as drawCircularFragment(), but WITH a visible black outline, so
// the centromere constriction notch actually reads against it) with a
// pinched constriction at 12 o'clock marking the centromere, and a
// grey centromere ellipse centered there. The SAME shape
// drawPaintedFragment() draws for linear chromosomes, just rotated 90
// degrees: its long axis runs radially (up/down at 12 o'clock) instead
// of tangentially (along the bar).
void Karyogram::drawRingFragment(
    cairo_t* cr,
    double centerX,
    double centerY,
    double diameter,
    RGB color,
    double centromereStartFraction,
    double centromereEndFraction,
    bool dashedOutline)						// Rings from a second homolog get dashed inner and outer outlines, matching the dashed outline of second homologs.
{
    const double radius = diameter / 2.0;
    const double dashPattern[2] = { 4.0, 3.0 };				// Same dash as drawPaintedFragment()'s second-homolog outline.
    const double innerRadius = diameter / 4.0;
    const double topAngle = -M_PI / 2.0; 					// 12 o'clock, same convention as drawPieChartFragment().

    const double constrictionArcLength = 7.0; 					// Matches drawPaintedFragment()'s constrictionHeight.
    const double constrictionAmount = 2.0; 					// Matches drawPaintedFragment()'s inward pinch depth.
    const double notchHalfAngle = (constrictionArcLength / 2.0) / radius; 	// Shared by both circles, so the two notches line up at the same angular position.

    const auto traceNotchedCircle = [&](double circleRadius, double pinchRadius)
    {
        cairo_new_path(cr);
        cairo_arc(cr, centerX, centerY, circleRadius, topAngle + notchHalfAngle, topAngle - notchHalfAngle + 2.0 * M_PI);

        const double pinchX = centerX + pinchRadius * std::cos(topAngle);
        const double pinchY = centerY + pinchRadius * std::sin(topAngle);

        cairo_line_to(cr, pinchX, pinchY); 					// In to the pinch point...
        cairo_close_path(cr); 							// ...and back out to close the notch.
    };

    // Outer colored band, pinched INWARD (toward center) at 12 o'clock.
    traceNotchedCircle(radius, radius - constrictionAmount);

    cairo_set_source_rgb(cr, color.r, color.g, color.b);
    cairo_fill_preserve(cr);

    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
    cairo_set_line_width(cr, 1.5);
    if (dashedOutline)
    {
        cairo_set_dash(cr, dashPattern, 2, 0.0);
    }
    cairo_stroke(cr);

    const double innerPinchRadius = std::min(radius, innerRadius + constrictionAmount); 	// Pinches OUTWARD into the colored band, meeting the outer ring's own inward pinch partway, clamped so it can't cross past the outer boundary.

    traceNotchedCircle(innerRadius, innerPinchRadius);

    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_fill_preserve(cr);

    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
    cairo_set_line_width(cr, 1.5);
    cairo_stroke(cr);							// Still dashed here when dashedOutline is set.
    cairo_set_dash(cr, nullptr, 0, 0.0);				// Back to solid for everything drawn afterwards.

    // Grey centromere mark, spanning the full width of the colored band
    // (outer edge to inner edge), centered at the band's radial
    // midpoint - still at 12 o'clock, still rotated 90 degrees from the
    // linear version (long axis radial instead of tangential).
    const double bandThickness = radius - innerRadius;
    const double bandMidRadius = (radius + innerRadius) / 2.0;

    const double centromereFraction = std::max(0.0, centromereEndFraction - centromereStartFraction);
    const double radialHeight = bandThickness;

    // Safeguard: keep the mark narrower than it is tall, so it always reads as a radial (vertical at 12 o'clock) centromere. On a small ring
    // the centromere is a large share of the ring (e.g. a 14.6 Mbp ring holding the 2.7 Mbp generic centromere of -ko), and
    // sizing the width by that share alone made the mark wider than the band and turned it horizontal.
    const double maxTangentialWidth = 0.6 * radialHeight;
    const double tangentialWidth = std::min(std::max(5.0, centromereFraction * M_PI * diameter), maxTangentialWidth);

    const double markCenterX = centerX + bandMidRadius * std::cos(topAngle);
    const double markCenterY = centerY + bandMidRadius * std::sin(topAngle);

    cairo_save(cr);
    cairo_translate(cr, markCenterX, markCenterY);
    cairo_scale(cr, tangentialWidth / 2.0, radialHeight / 2.0);
    cairo_arc(cr, 0.0, 0.0, 1.0, 0.0, 2.0 * M_PI);
    cairo_restore(cr);

    cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
    cairo_fill(cr);
}














// Computes per-excised-piece heights (no floor, unlike computeSDRbarHeight) and how many vertical-stack columns are needed so no column exceeds
// maxColumnHeight. Shared between computeSlotWidth() (spacing) and drawStackedMutations() (actual drawing) so both always agree on layout.
int Karyogram::computeExcisedColumnLayout(
    const std::vector<const SDRdataRecord*>& excisedRecords,
    double maxLengthMbp,
    double maxRenderHeight,
    double maxColumnHeight,
    double verticalGap,
    std::vector<double>& outHeights)
{
    outHeights.clear();
    outHeights.reserve(excisedRecords.size());

    for (const SDRdataRecord* record : excisedRecords)
    {
        double lengthMbp = 0.0;

        for (const SDRfragment& fragment : record->fragments)
        {
            lengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
        }

	const double height = (maxLengthMbp > 0.0) ? (maxRenderHeight * (lengthMbp / maxLengthMbp)) : 0.0;
        const double minPieceSize = 5.0;				// Keep tiny excised/ecDNA pieces (bars and circles alike) from shrinking below a visible size.
        outHeights.push_back(std::max(minPieceSize, height));

    }

    if (outHeights.empty())
    {
        return 0;
    }

    // Greedy first-fit packing: keep adding pieces to the current column as long as they fit within
    // maxColumnHeight; once the next piece would overflow it, start a new column instead. Pieces are
    // packed in the order given (not re-sorted), so column balance depends on caller ordering.
    int columns = 1;
    double columnUsedHeight = 0.0;
    bool columnHasContent = false;

    for (double pieceHeight : outHeights)
    {
        const double additionalHeight = columnHasContent ? (verticalGap + pieceHeight) : pieceHeight;

	// Increase the number of columns for excised/acentric fragments if they exist, and increase it by the required additional height 
        if (columnHasContent && (columnUsedHeight + additionalHeight) > maxColumnHeight)
        {
            ++columns;
            columnUsedHeight = 0.0;
            columnHasContent = false;
        }

        columnUsedHeight += columnHasContent ? (verticalGap + pieceHeight) : pieceHeight;
        columnHasContent = true;
    }

    return columns;
}


















// The same height-scaling formula drawStackedMutations uses per bar. Pulled out so label positioning can reuse it without duplicating
// the clamp logic.
double Karyogram::computeSDRbarHeight(
    double lengthMbp, 								// Length of the bar fragment
    double maxLengthMbp, 							// Max length of all chromosomes for normalizing chromosome heights
    double maxRenderHeight)							// Max render height of the entire karyogram image
{
    double barHeight = maxRenderHeight * (lengthMbp / maxLengthMbp);		// Calculate segment height and on karyogram in pixels

    if (barHeight < 30.0)							// Prevent segment from being too small
    {
        barHeight = 30.0;
    }

    return barHeight;
}









// Tallest bar among a slot's (possibly stacked) records, used to position that slot's label just below whatever actually got drawn.
double Karyogram::computeMaxBarHeight(
    const std::vector<const SDRdataRecord*>& records, 				// Access the fragment sizes
    double maxLengthMbp, 							// Max length in Mbp of the chromosomes for normalization
    double maxRenderHeight,							// Max height of the entire karyogram for scaling
    int homeOldStrandID,							// Adjust label height to the home slot it was designated for
    const SDRmasterHeader& masterHeader)					// Access intact chromosome sizes
{
    if (records.empty())
    {
	return 0.0;
    }


    // These are the same special shapes that can cause the remaining chromosome to be drawn at its ORIGINAL height
    // rather than its reduced fragment length.
    const SDRdataRecord* gapRemainingRecord = nullptr;
    std::vector<const SDRdataRecord*> gapPieceRecords;
    const bool isGapShape = analyzeHomeStrandGapShape(records, homeOldStrandID, masterHeader, gapRemainingRecord, gapPieceRecords);


    // For deletion-like chromosome drawings, the remaining chromosome is drawn using the full original chromosome
    // length, with the deleted material represented as a white gap.
    if (isGapShape)
    {
        const std::size_t sizeIndex = static_cast<std::size_t>(homeOldStrandID);

        if (sizeIndex < masterHeader.intactChromosomeSizes.size())
        {
            const double originalSizeMbp = masterHeader.intactChromosomeSizes[sizeIndex];

            if (originalSizeMbp > 0.0)
            {
                return computeSDRbarHeight(originalSizeMbp, maxLengthMbp, maxRenderHeight);
            }
        }
    }


    // Normal case: determine the tallest actual record.
    double tallest = 0.0;

    for (const SDRdataRecord* record : records)
    {
        double lengthMbp = 0.0;

        for (const SDRfragment& fragment : record->fragments)
        {
            lengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
        }

        if (record->linear)							// A mixed record's lost home-strand gaps are drawn as white space, so they count towards its height.
        {
            bool hasHomeGap = false;
            const double lengthWithGaps = mixedRecordLengthWithHomeGaps(*record, homeOldStrandID, hasHomeGap);
            lengthMbp = hasHomeGap ? lengthWithGaps : lengthMbp;

            // Safeguard: a record made only of home-strand material is never drawn taller than the intact chromosome. When its fragments are shuffled
            // or inverted (e.g. 1/132.6/131.7, 1/124.2/104.9, 1/124.2/131.7, 1/104.9/0), neighbouring fragments don't touch and the gap rule above
            // counts the same stretch twice, which pushed the label far below the chromosome (and into the legend on the last row).
            const bool homeStrandOnly = std::all_of(record->fragments.begin(), record->fragments.end(),
                [homeOldStrandID](const SDRfragment& fragment) { return fragment.oldStrandID == homeOldStrandID; });
            const std::size_t sizeIndex = static_cast<std::size_t>(homeOldStrandID);

            if (homeStrandOnly && sizeIndex < masterHeader.intactChromosomeSizes.size() && masterHeader.intactChromosomeSizes[sizeIndex] > 0.0)
            {
                lengthMbp = std::min(lengthMbp, masterHeader.intactChromosomeSizes[sizeIndex]);
            }
        }

        const double height =computeSDRbarHeight(lengthMbp, maxLengthMbp, maxRenderHeight);

        if (height > tallest)
        {
            tallest = height;
        }
    }

    return tallest;
}











// Function to determine how wide a slot's drawing will actually be, used both by drawStackedMutations() itself (to position things)
// and by renderSDRkaryogram() (to space homologs apart correctly). Returns the total width; outColumn2Width is set to the second
// column's width for deletion-like/ecDNA shapes (chromosomeWidth otherwise, unused by the caller in that case).
double Karyogram::computeSlotWidth(
    const std::vector<const SDRdataRecord*>& records,
    int homeOldStrandID,
    double chromosomeWidth,
    double maxLengthMbp,
    double maxRenderHeight,
    const SDRmasterHeader& masterHeader,
    double& outColumn2Width)
{
    outColumn2Width = chromosomeWidth;

    if (records.empty())
    {
        return 0.0;
    }

    const SDRdataRecord* gapRemainingRecord = nullptr;
    std::vector<const SDRdataRecord*> gapPieceRecords;
    const bool isGapShape = analyzeHomeStrandGapShape(records, homeOldStrandID, masterHeader, gapRemainingRecord, gapPieceRecords);

    const bool isChromosomeRing = isChromosomeRingShape(records, homeOldStrandID, masterHeader);


    const std::size_t sizeIndex = static_cast<std::size_t>(homeOldStrandID);
    const double originalSizeMbp = (sizeIndex < masterHeader.intactChromosomeSizes.size()) ? masterHeader.intactChromosomeSizes[sizeIndex] : 0.0;

    const double stackGap = 6.0;


    if (isGapShape && originalSizeMbp > 0.0)
    {
        // One column for the remaining bar. A second column exists ONLY if excised pieces are present in the file, sized to the widest ecDNA circle (if any).
        if (gapPieceRecords.empty())
        {
            return chromosomeWidth;
        }

        const double remainingBarHeight = computeSDRbarHeight(originalSizeMbp, maxLengthMbp, maxRenderHeight);
        double maxDiameter = 0.0;

        for (const SDRdataRecord* piece : gapPieceRecords)
        {
            if (piece->linear)
            {
                continue;
            }

            double pieceLengthMbp = 0.0;

            for (const SDRfragment& fragment : piece->fragments)
            {
                pieceLengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
            }

            const double minPieceSize = 5.0;			// Same floor the drawing applies to each circle
            maxDiameter = std::max(maxDiameter, std::max(minPieceSize, std::min((pieceLengthMbp / originalSizeMbp) * remainingBarHeight, chromosomeWidth * 3.0)));
        }

        outColumn2Width = std::max(chromosomeWidth, maxDiameter);

        return chromosomeWidth + (maxDiameter > 0.0 ? circleColumnGap : stackGap) + outColumn2Width;	// Must match drawStackedMutations()'s gap choice.
    }


    if (isChromosomeRing && originalSizeMbp > 0.0)
    {
        const SDRdataRecord* ringRecord = nullptr;

        for (const SDRdataRecord* record : records)
        {
            if (!record->linear)
            {
                ringRecord = record;
                break;
            }
        }

        double ringLengthMbp = 0.0;

        if (ringRecord != nullptr)
        {
            for (const SDRfragment& fragment : ringRecord->fragments)
            {
                ringLengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
            }
        }

        const double remainingBarHeight = computeSDRbarHeight(originalSizeMbp, maxLengthMbp, maxRenderHeight);
        const double equivalentLinearHeight = (ringLengthMbp / originalSizeMbp) * remainingBarHeight;
        const double minPieceSize = 5.0;				// Keep a tiny ring from shrinking below a visible size.
        const double ringDiameter = std::max(minPieceSize, std::min(equivalentLinearHeight, chromosomeWidth * 3.0));

        outColumn2Width = std::max(chromosomeWidth, ringDiameter);

        return chromosomeWidth + circleColumnGap + outColumn2Width;
    }


    if (isDeletionTranslocationDonorShape(records, homeOldStrandID))
    {
        // Two plain columns, both chromosomeWidth-wide, the recombined piece isn't scaled to full-original-length, so no ecDNA-style
        // diameter accounting is needed here.
        return 2 * chromosomeWidth + stackGap;
    }


    // Multiple records, none matching a named shape, this is the shattered/collapsed-concatenation case (chromothripsis-style):
    // one base record (most fragments) plus N unrelated leftover pieces, none of which decompose into a recognized shape.
    if (records.size() > 1)
    {
        const SDRdataRecord* baseRecord = nullptr;
        std::vector<const SDRdataRecord*> otherRecords;

        for (const SDRdataRecord* record : records)
        {
            if (baseRecord == nullptr || record->fragments.size() > baseRecord->fragments.size())
            {
                if (baseRecord != nullptr)
                {
                    otherRecords.push_back(baseRecord);
                }

                baseRecord = record;
            }
            else
            {
                otherRecords.push_back(record);
            }
        }

        const double verticalGap = 6.0;
        std::vector<double> pieceHeights;
        const int numColumns = computeExcisedColumnLayout(otherRecords, maxLengthMbp, maxRenderHeight, maxRenderHeight, verticalGap, pieceHeights);

        outColumn2Width = chromosomeWidth; 					// Each stacked column is one chromosomeWidth wide.

        return chromosomeWidth + stackGap + static_cast<double>(numColumns) * chromosomeWidth + static_cast<double>(std::max(0, numColumns - 1)) * stackGap;
    }

    const std::size_t count = records.size();
    double totalWidth = count * chromosomeWidth;

    if (count > 0)
    {
        totalWidth += (count - 1) * stackGap;
    }

    return totalWidth;

}












// Synthesizes a plain, unmutated baseline record for an original strand that has NO data records at all in the SDR file, e.g. a
// file that only lists mutated strands, omitting trivial "nothing happened here" entries for every untouched chromosome. Without
// this, such a strand's slot would be left completely blank on the karyogram instead of showing a correctly-sized intact chromosome.
std::vector<const SDRdataRecord*> Karyogram::synthesizeIntactRecordIfMissing(
    int oldStrandID,
    int cellID,
    const SDRmasterHeader& masterHeader,
    std::vector<SDRdataRecord>& intactRecordStorage)
{
    const std::size_t sizeIndex = static_cast<std::size_t>(oldStrandID);
    const double sizeMbp = (sizeIndex < masterHeader.intactChromosomeSizes.size()) ? masterHeader.intactChromosomeSizes[sizeIndex] : 0.0;

    if (sizeMbp <= 0.0)
    {
        return {};
    }

    SDRdataRecord record{};
    record.cellID = cellID;
    record.newStrandID = oldStrandID;				 // Matches the file's own baseline convention (newStrandID == oldStrandID when unmutated)
    record.linear = true;

    SDRfragment fragment{};
    fragment.oldStrandID = oldStrandID;
    fragment.oldStartPosition = 0.0;
    fragment.oldEndPosition = sizeMbp;
    fragment.hasCentromere = true; 				// Matches real baseline records, which always flag the whole chromosome as centromere-bearing

    record.fragments.push_back(fragment);

    intactRecordStorage.push_back(record);

    return { &intactRecordStorage.back() };
}











// If a slot contains any genuine mutation-derived record (newStrandID > numOriginalStrands), drop any leftover baseline/unmutated record
// for that same original strand - only the resulting mutated strand(s) should be drawn, not the untouched original alongside them.
std::vector<const SDRdataRecord*> Karyogram::filterBaselineIfMutated(
    const std::vector<const SDRdataRecord*>& records,					// Access SDR data records to assess which records contain mutated strands and which are intact
    int numOriginalStrands)								// Use to assess when the mutated SDR data records begin
{

    bool hasMutation = false;

    for (const SDRdataRecord* record : records)
    {
        if (record->newStrandID > numOriginalStrands)					// New strand IDs > 46 for human chromosomes are aberrant/mutated strands
        {
            hasMutation = true;
            break;
        }
    }

    if (!hasMutation)									// If no mutated strands found, return the SDR records of the intact strands
    {
        return records;
    }

    std::vector<const SDRdataRecord*> filtered;						// Vector 'filtered' to store unique data records as being either completely intact, or partially mutated, not both intact and mutated.

    for (const SDRdataRecord* record : records)
    {
        if (record->newStrandID <= numOriginalStrands)					// If new strand ID > 46 for human chromosomes, this is a mutated/aberrant fragment strand
        {
            continue; 									// drop the baseline (intact) data entry now that a real mutation exists for this strand
        }

        filtered.push_back(record);
    }

    return filtered;
}
















// The first record's full fragment list is the base structure. For every subsequent record, its own home-strand fragment(s) 
// tell us what portion of the original chromosome that record actually retains.
// The base structure gets truncated down to fit within that retained range (keeping only the overlapping part of each base fragment, dropping any base 
// fragment with no overlap at all). The subsequent record's FOREIGN fragments are then tacked on at the end.
SDRdataRecord Karyogram::concatenateRecordsForDrawing(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID)
{
    SDRdataRecord merged{};
    merged.cellID = records.front()->cellID;
    merged.newStrandID = records.front()->newStrandID; 					// Representative ID, used only for the drawn label.
    merged.linear = records.front()->linear;

    std::vector<SDRfragment> baseFragments(records.front()->fragments.begin(), records.front()->fragments.end());

    for (std::size_t i = 1; i < records.size(); ++i)
    {
        for (const SDRfragment& fragment : records[i]->fragments)
        {
            if (fragment.oldStrandID != homeOldStrandID)
            {
                continue;
            }

            const double retainLo = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
            const double retainHi = std::max(fragment.oldStartPosition, fragment.oldEndPosition);

            std::vector<SDRfragment> truncated;

            for (SDRfragment base : baseFragments)
            {
                const double baseLo = std::min(base.oldStartPosition, base.oldEndPosition);
                const double baseHi = std::max(base.oldStartPosition, base.oldEndPosition);
                const bool baseReversed = base.oldStartPosition > base.oldEndPosition;

                const double overlapLo = std::max(baseLo, retainLo);
                const double overlapHi = std::min(baseHi, retainHi);

                if (overlapHi <= overlapLo)
                {
                    continue; 							// No overlap with the retained range, fully dropped.
                }

                // Truncate to just the overlapping portion, preserving orientation.
                if (baseReversed)
                {
                    base.oldStartPosition = overlapHi;
                    base.oldEndPosition = overlapLo;
                }
                else
                {
                    base.oldStartPosition = overlapLo;
                    base.oldEndPosition = overlapHi;
                }

                truncated.push_back(base);
            }

            baseFragments = truncated;
        }
    }

    for (const SDRfragment& fragment : baseFragments)
    {
        merged.fragments.push_back(fragment);
    }

    for (std::size_t i = 1; i < records.size(); ++i)
    {
        for (const SDRfragment& fragment : records[i]->fragments)
        {
            if (fragment.oldStrandID != homeOldStrandID)
            {
                merged.fragments.push_back(fragment);
            }
        }
    }

    return merged;
}















// Decides whether a slot's multiple records should be merged into one synthetic drawing (via concatenateRecordsForDrawing) or left as separate pieces. 
// Named shapes (deletion, ecDNA, deletion-inversion, deletion-translocation-donor) are returned untouched, since drawStackedMutations() has its own 
// dedicated layout for each. Everything else attempts the merge, but falls back to the original records if the result comes out degenerate.
std::vector<const SDRdataRecord*> Karyogram::clusterRecordsForDrawing(
    const std::vector<const SDRdataRecord*>& records,
    int homeOldStrandID,
    std::vector<SDRdataRecord>& mergedStorage,
    const SDRmasterHeader& masterHeader)
{
    // Nothing to merge. Built element by element rather than copied with "return records;": GCC 12-14 at -O2 wrongly
    // warns (-Warray-bounds) that copying a vector known to hold at most one pointer reads past its end.
    if (records.empty())
    {
        return {};
    }

    if (records.size() == 1)
    {
        return {records.front()};
    }

    const SDRdataRecord* gapRemainingRecord = nullptr;
    std::vector<const SDRdataRecord*> gapPieceRecords;

    if (analyzeHomeStrandGapShape(records, homeOldStrandID, masterHeader, gapRemainingRecord, gapPieceRecords)
        || isDeletionTranslocationDonorShape(records, homeOldStrandID)
	|| isChromosomeRingShape(records, homeOldStrandID, masterHeader))
    {
        return records;
    }

    SDRdataRecord merged = concatenateRecordsForDrawing(records, homeOldStrandID);

    double mergedLengthMbp = 0.0;

    for (const SDRfragment& fragment : merged.fragments)
    {
        mergedLengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
    }

    double firstRecordLengthMbp = 0.0;

    for (const SDRfragment& fragment : records.front()->fragments)
    {
        firstRecordLengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
    }

    // If concatenation collapsed to (near) nothing relative to the base record it started from, the merge assumptions this function was
    // built for don't hold for this group - fall back to drawing every record separately instead of returning a degenerate result.
    if (firstRecordLengthMbp > 0.0 && (mergedLengthMbp / firstRecordLengthMbp) < 0.1)
    {
        return records;
    }

    mergedStorage.push_back(merged);

    return { &mergedStorage.back() };
}



// --------------------------------------------------------------------------------------------------------------------------------------------
// Zoomed homolog view: every piece drawn in a homolog's slot, recorded at its true size for the web viewer (see getDetailSlots()).
// --------------------------------------------------------------------------------------------------------------------------------------------

namespace
{
    // A recorded segment plus where its bases came from, used to place the cell's indels on it.
    struct DetailSegmentOrigin
    {
        int strandID;					// Original strand, -1 for white gaps or when not known
        double lowMbp;					// Original positions covered by the segment
        double highMbp;
    };

    // Adds a marker to the piece for every indel whose site lies on one of its segments, at the site's position along the piece.
    // Indels are taken out before drawing (see removeIndels()), so without a marker they could not be found in the zoomed view.
    void addIndelMarkers(KaryogramDetailPiece& piece, const std::vector<DetailSegmentOrigin>& origins, const std::vector<SDRindelEvent>& indels)
    {
        const double siteTolerance = 0.0000005;		// Half a base: the site is the end of the piece left of it

        for (const SDRindelEvent& indel : indels)
        {
            for (std::size_t k = 0; k < piece.segments.size() && k < origins.size(); ++k)
            {
                const KaryogramDetailSegment& segment = piece.segments[k];
                const DetailSegmentOrigin& origin = origins[k];

                if (segment.gap || origin.strandID != indel.oldStrandID
                    || indel.position < origin.lowMbp - siteTolerance || indel.position > origin.highMbp + siteTolerance)
                {
                    continue;
                }

                const double offsetMbp = segment.reversed ? origin.highMbp - indel.position : indel.position - origin.lowMbp;
                const double positionMbp = std::min(segment.endMbp, std::max(segment.startMbp, segment.startMbp + offsetMbp));

                std::string text = "Indel: ";
                if (indel.deletedBases > 0)
                {
                    text += std::to_string(indel.deletedBases) + " bp deleted";
                }
                if (indel.insertedBases > 0)
                {
                    text += std::string(indel.deletedBases > 0 ? ", " : "") + std::to_string(indel.insertedBases) + " bp inserted";
                }

                KaryogramDetailMarker marker;
                marker.positionMbp = positionMbp;
                marker.text = text;
                marker.deletedMbp = indel.deletedBases / 1000000.0;	// The joined piece still spans the deleted bases, right below the site
                marker.inserted = indel.insertedBases > 0;
                piece.markers.push_back(marker);
                break;							// One marker per indel, even where two segments touch at its site
            }
        }
    }

    bool isWhiteSegment(const RGB& color)
    {
        return std::fabs(color.r - 1.0) < 1e-6 && std::fabs(color.g - 1.0) < 1e-6 && std::fabs(color.b - 1.0) < 1e-6;
    }
}


// Opens a slot for the zoomed view: pieces drawn from now on until endDetailSlot() belong to it. Alters detailSlots and the current slot state.
int Karyogram::beginDetailSlot(const std::string& title, int homeOldStrandID)
{
    KaryogramDetailSlot slot;
    slot.title = title;
    detailSlots.push_back(slot);

    currentDetailSlot = static_cast<int>(detailSlots.size()) - 1;
    detailHomeStrandID = homeOldStrandID;
    return currentDetailSlot;
}


// Closes the current slot; a slot that recorded no piece is removed again so it never shows up as an empty zoomed view.
void Karyogram::endDetailSlot()
{
    if (currentDetailSlot >= 0 && currentDetailSlot == static_cast<int>(detailSlots.size()) - 1 && detailSlots.back().pieces.empty())
    {
        detailSlots.pop_back();
    }

    currentDetailSlot = -1;
    detailHomeStrandID = -1;
}


// Records one bar drawn by drawPaintedFragment() in the current slot, at its true size. The drawing stretches tiny segments to a minimum height,
// so the pixel positions of chevrons and arrows are mapped back through the drawn and true segment boundaries to Mbp along the piece.
// Alters detailSlots (adds a piece to the current slot).
void Karyogram::recordDetailBar(double x, double y, double height, double width,
                                const std::vector<PaintedSegment>& trueSegments, const std::vector<PaintedSegment>& drawnSegments,
                                const std::vector<std::pair<double, double>>& chevronsPx, const std::vector<std::pair<double, double>>& arrowsPx,
                                const std::vector<bool>& dashedBoundary, const std::vector<bool>& hiddenBoundary,
                                bool roundTopCap, bool roundBottomCap, bool drawOutline)
{
    if (currentDetailSlot < 0 || currentDetailSlot >= static_cast<int>(detailSlots.size()) || trueSegments.empty() || height <= 0.0)
    {
        return;
    }

    // Mbp per unit of fraction along the bar: from the segments that carry their original positions, otherwise (a deletion drawn on its
    // own chromosome, where fractions are positions over the chromosome size) the home chromosome's size, otherwise the karyogram scale.
    double knownMbp = 0.0;
    double knownFraction = 0.0;

    for (const PaintedSegment& segment : trueSegments)
    {
        if (segment.oldStrandID >= 0 && !isWhiteSegment(segment.color) && segment.endFraction > segment.startFraction)
        {
            knownMbp += std::fabs(segment.oldEndPosition - segment.oldStartPosition);
            knownFraction += segment.endFraction - segment.startFraction;
        }
    }

    double mbpPerFraction = height * detailMbpPerPixel;
    bool fractionsArePositions = false;					// Fractions are home-chromosome positions over its size

    if (knownFraction > 1e-12 && knownMbp > 0.0)
    {
        mbpPerFraction = knownMbp / knownFraction;
    }
    else if (detailMasterHeader != nullptr && detailHomeStrandID > 0
             && static_cast<std::size_t>(detailHomeStrandID) < detailMasterHeader->intactChromosomeSizes.size()
             && detailMasterHeader->intactChromosomeSizes[detailHomeStrandID] > 0.0)
    {
        mbpPerFraction = detailMasterHeader->intactChromosomeSizes[detailHomeStrandID];
        fractionsArePositions = true;
    }

    // Drawn fraction -> true fraction, through the boundaries both versions share (the minimum-height stretch is piecewise linear).
    std::vector<std::pair<double, double>> boundaryPairs;
    for (std::size_t k = 0; k < trueSegments.size() && k < drawnSegments.size(); ++k)
    {
        boundaryPairs.emplace_back(drawnSegments[k].startFraction, trueSegments[k].startFraction);
        boundaryPairs.emplace_back(drawnSegments[k].endFraction, trueSegments[k].endFraction);
    }
    std::sort(boundaryPairs.begin(), boundaryPairs.end());

    const auto pixelToMbp = [&](double pixelY)
    {
        const double drawnFraction = (pixelY - y) / height;
        double trueFraction = drawnFraction;

        if (!boundaryPairs.empty())
        {
            if (drawnFraction <= boundaryPairs.front().first)
            {
                trueFraction = boundaryPairs.front().second;
            }
            else if (drawnFraction >= boundaryPairs.back().first)
            {
                trueFraction = boundaryPairs.back().second;
            }
            else
            {
                for (std::size_t i = 0; i + 1 < boundaryPairs.size(); ++i)
                {
                    const auto& [drawnLow, trueLow] = boundaryPairs[i];
                    const auto& [drawnHigh, trueHigh] = boundaryPairs[i + 1];

                    if (drawnFraction <= drawnHigh)
                    {
                        const double t = (drawnHigh > drawnLow) ? (drawnFraction - drawnLow) / (drawnHigh - drawnLow) : 0.0;
                        trueFraction = trueLow + t * (trueHigh - trueLow);
                        break;
                    }
                }
            }
        }

        return trueFraction * mbpPerFraction;
    };

    KaryogramDetailPiece piece;
    piece.x = x;
    piece.y = y;
    piece.width = width;
    piece.height = height;
    piece.caption = drawOutline ? "" : "Excised piece";
    piece.lengthMbp = mbpPerFraction;
    piece.roundTopCap = roundTopCap;
    piece.roundBottomCap = roundBottomCap;
    piece.outline = drawOutline;

    std::vector<DetailSegmentOrigin> origins;

    for (const PaintedSegment& segment : trueSegments)
    {
        KaryogramDetailSegment detail;
        detail.startMbp = segment.startFraction * mbpPerFraction;
        detail.endMbp = segment.endFraction * mbpPerFraction;
        detail.r = segment.color.r;
        detail.g = segment.color.g;
        detail.b = segment.color.b;
        detail.gap = isWhiteSegment(segment.color);
        detail.reversed = segment.isReversed;
        detail.secondHomolog = segment.isSecondHomolog;

        if (segment.hasCentromere)
        {
            detail.centromereStartMbp = segment.centromereStartFraction * mbpPerFraction;
            detail.centromereEndMbp = segment.centromereEndFraction * mbpPerFraction;
        }

        piece.segments.push_back(detail);

        DetailSegmentOrigin origin{ -1, 0.0, 0.0 };
        if (!detail.gap && segment.oldStrandID >= 0)
        {
            origin = { segment.oldStrandID, std::min(segment.oldStartPosition, segment.oldEndPosition), std::max(segment.oldStartPosition, segment.oldEndPosition) };
        }
        else if (!detail.gap && fractionsArePositions)
        {
            origin = { detailHomeStrandID, detail.startMbp, detail.endMbp };
        }
        origins.push_back(origin);
    }

    // Original span of every segment, shown when hovering over it in the zoomed view. A white gap's deleted span lies between the ends of the
    // pieces of the same chromosome that border it; on a deletion drawn on its own chromosome the gap's position is already that span.
    for (std::size_t k = 0; k < piece.segments.size(); ++k)
    {
        KaryogramDetailSegment& detail = piece.segments[k];

        if (!detail.gap)
        {
            if (origins[k].strandID >= 0)
            {
                detail.originStartMbp = origins[k].lowMbp;
                detail.originEndMbp = origins[k].highMbp;
                detail.foreign = detailHomeStrandID > 0 && origins[k].strandID != detailHomeStrandID;

                const auto name = detailStrandNames.find(origins[k].strandID);
                if (name != detailStrandNames.end())
                {
                    detail.originName = name->second;
                }
            }
            continue;
        }

        if (fractionsArePositions)
        {
            detail.originStartMbp = detail.startMbp;
            detail.originEndMbp = detail.endMbp;
            continue;
        }

        if (k == 0 || k + 1 >= piece.segments.size() || origins[k - 1].strandID < 0 || origins[k - 1].strandID != origins[k + 1].strandID)
        {
            continue;
        }

        const double aboveEdge = piece.segments[k - 1].reversed ? origins[k - 1].lowMbp : origins[k - 1].highMbp;
        const double belowEdge = piece.segments[k + 1].reversed ? origins[k + 1].highMbp : origins[k + 1].lowMbp;
        detail.originStartMbp = std::min(aboveEdge, belowEdge);
        detail.originEndMbp = std::max(aboveEdge, belowEdge);
    }

    for (const auto& [topY, bottomY] : chevronsPx)
    {
        piece.chevrons.emplace_back(pixelToMbp(topY), pixelToMbp(bottomY));
    }

    for (const auto& [fromY, toY] : arrowsPx)
    {
        piece.movedArrows.emplace_back(pixelToMbp(fromY), pixelToMbp(toY));
    }

    for (std::size_t k = 0; k < dashedBoundary.size(); ++k)
    {
        if (dashedBoundary[k])
        {
            piece.dashedBoundaries.push_back(static_cast<int>(k));
        }
    }

    for (std::size_t k = 0; k < hiddenBoundary.size(); ++k)
    {
        if (hiddenBoundary[k])
        {
            piece.hiddenBoundaries.push_back(static_cast<int>(k));
        }
    }

    addIndelMarkers(piece, origins, detailIndels);
    detailSlots[currentDetailSlot].pieces.push_back(std::move(piece));
}


// Records a circular piece (ecDNA or ring) in the current slot, clickable over the drawn circle: its segments are laid end to end from its
// first fragment, which the zoomed view wraps clockwise around a circle.
// Alters detailSlots (adds a piece to the current slot).
void Karyogram::recordDetailCircle(const SDRdataRecord& record, double centerX, double centerY, double diameter)
{
    if (currentDetailSlot < 0 || currentDetailSlot >= static_cast<int>(detailSlots.size()) || detailMasterHeader == nullptr)
    {
        return;
    }

    const std::vector<PaintedSegment> segments = buildPaintedSegments(record, detailHumanGenome, *detailMasterHeader);
    double lengthMbp = 0.0;
    bool hasCentromere = false;

    for (const SDRfragment& fragment : record.fragments)
    {
        lengthMbp += std::fabs(fragment.oldEndPosition - fragment.oldStartPosition);
        hasCentromere = hasCentromere || fragment.hasCentromere;
    }

    if (segments.empty() || lengthMbp <= 0.0)
    {
        return;
    }

    KaryogramDetailPiece piece;
    piece.x = centerX - diameter / 2.0;
    piece.y = centerY - diameter / 2.0;
    piece.width = diameter;
    piece.height = diameter;
    piece.caption = hasCentromere ? "Ring" : "ecDNA";
    piece.lengthMbp = lengthMbp;
    piece.circular = true;
    piece.roundTopCap = false;							// A circle has no telomeres
    piece.roundBottomCap = false;

    std::vector<DetailSegmentOrigin> origins;

    for (const PaintedSegment& segment : segments)
    {
        KaryogramDetailSegment detail;
        detail.startMbp = segment.startFraction * lengthMbp;
        detail.endMbp = segment.endFraction * lengthMbp;
        detail.r = segment.color.r;
        detail.g = segment.color.g;
        detail.b = segment.color.b;
        detail.reversed = segment.isReversed;
        detail.secondHomolog = segment.isSecondHomolog;

        if (segment.hasCentromere)
        {
            detail.centromereStartMbp = segment.centromereStartFraction * lengthMbp;
            detail.centromereEndMbp = segment.centromereEndFraction * lengthMbp;
        }

        if (detail.reversed)
        {
            piece.chevrons.emplace_back(detail.startMbp, detail.endMbp);
        }

        detail.originStartMbp = std::min(segment.oldStartPosition, segment.oldEndPosition);
        detail.originEndMbp = std::max(segment.oldStartPosition, segment.oldEndPosition);
        piece.segments.push_back(detail);
        origins.push_back({ segment.oldStrandID, detail.originStartMbp, detail.originEndMbp });
    }

    addIndelMarkers(piece, origins, detailIndels);
    detailSlots[currentDetailSlot].pieces.push_back(std::move(piece));
}

} // namespace karyoke
