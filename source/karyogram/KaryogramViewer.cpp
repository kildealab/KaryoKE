#include "KaryogramViewer.h"
#include "generalUtilities.h"
#include "KaryoKEversion.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>




namespace karyoke
{

namespace
{
    // Encode PNG bytes as base64 so they can be embedded in the HTML as data: URLs.
    std::string base64Encode(const std::vector<unsigned char>& bytes)
    {
        static const char table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

        std::string encoded;
        encoded.reserve(((bytes.size() + 2) / 3) * 4);

        size_t i = 0;
        for (; i + 2 < bytes.size(); i += 3)                                    // Every 3 input bytes become 4 output characters.
        {
            const unsigned int chunk = (bytes[i] << 16) | (bytes[i + 1] << 8) | bytes[i + 2];
            encoded += table[(chunk >> 18) & 0x3F];
            encoded += table[(chunk >> 12) & 0x3F];
            encoded += table[(chunk >> 6) & 0x3F];
            encoded += table[chunk & 0x3F];
        }

        const size_t remaining = bytes.size() - i;                             // Pad the last 1 or 2 bytes with '='.
        if (remaining == 1)
        {
            const unsigned int chunk = bytes[i] << 16;
            encoded += table[(chunk >> 18) & 0x3F];
            encoded += table[(chunk >> 12) & 0x3F];
            encoded += "==";
        }
        else if (remaining == 2)
        {
            const unsigned int chunk = (bytes[i] << 16) | (bytes[i + 1] << 8);
            encoded += table[(chunk >> 18) & 0x3F];
            encoded += table[(chunk >> 12) & 0x3F];
            encoded += table[(chunk >> 6) & 0x3F];
            encoded += '=';
        }

        return encoded;
    }

    // Escape text placed inside the HTML (file names can contain characters such as & or <).
    std::string htmlEscape(const std::string& text)
    {
        std::string escaped;
        for (char c : text)
        {
            switch (c)
            {
                case '&':  escaped += "&amp;";  break;
                case '<':  escaped += "&lt;";   break;
                case '>':  escaped += "&gt;";   break;
                case '"':  escaped += "&quot;"; break;
                case '\'': escaped += "&#39;";  break;
                default:   escaped += c;
            }
        }
        return escaped;
    }

    // Write text as a JavaScript string literal (with quotes) for use inside the viewer's <script>.
    // Why: file names and titles come from the input file name, and could contain quotes, backslashes or "</script>",
    // which would otherwise break the script. HTML escaping does not apply inside <script>, so JavaScript escaping is used.
    std::string jsString(const std::string& text)
    {
        std::string quoted = "\"";
        for (char c : text)
        {
            switch (c)
            {
                case '\\': quoted += "\\\\"; break;
                case '"':  quoted += "\\\""; break;
                case '<':  quoted += "\\u003c"; break;
                case '\n': quoted += "\\n"; break;
                case '\r': quoted += "\\r"; break;
                default:   quoted += c;
            }
        }
        return quoted + "\"";
    }

    // Write a number for the viewer's <script> with enough digits to keep single bases apart on the longest chromosomes (e.g. 248.956422 Mbp).
    std::string jsNumber(double value)
    {
        char buffer[32];
        std::snprintf(buffer, sizeof(buffer), "%.10g", value);
        return buffer;
    }

    // Write a 0-1 RGB colour as a CSS hex colour, e.g. "#3f8f6b".
    std::string cssColor(double r, double g, double b)
    {
        const auto channel = [](double value) { return static_cast<int>(std::lround(std::min(1.0, std::max(0.0, value)) * 255.0)); };
        char buffer[8];
        std::snprintf(buffer, sizeof(buffer), "#%02x%02x%02x", channel(r), channel(g), channel(b));
        return buffer;
    }

    // Write the clickable homologs of one karyogram as a JavaScript array for the zoomed homolog view: per slot its title, whether its pieces are
    // sister chromatids (SDD replicated homolog) and its pieces, each piece
    // with its clickable area on the karyogram and its segments, chevrons, arrows, boundary styles and indel markers in Mbp along the piece.
    std::string detailSlotsJs(const std::vector<KaryogramDetailSlot>& slots)
    {
        std::string js = "[";

        for (const KaryogramDetailSlot& slot : slots)
        {
            js += "\n      { title: " + jsString(slot.title) + (slot.sisterChromatids ? ", sisters: true" : "") + ", pieces: [";

            for (const KaryogramDetailPiece& piece : slot.pieces)
            {
                js += "\n        { x: " + jsNumber(piece.x) + ", y: " + jsNumber(piece.y) + ", w: " + jsNumber(piece.width) + ", h: " + jsNumber(piece.height)
                    + ", caption: " + jsString(piece.caption) + ", len: " + jsNumber(piece.lengthMbp)
                    + ", roundTop: " + (piece.roundTopCap ? "true" : "false") + ", roundBottom: " + (piece.roundBottomCap ? "true" : "false")
                    + ", outline: " + (piece.outline ? "true" : "false") + ", circular: " + (piece.circular ? "true" : "false") + ",\n          segs: [";

                for (const KaryogramDetailSegment& segment : piece.segments)		// [start, end, colour, gap, reversed, second homolog, centromere start, centromere end, original start, original end, origin name]
                {
                    js += "[" + jsNumber(segment.startMbp) + "," + jsNumber(segment.endMbp) + ",\"" + cssColor(segment.r, segment.g, segment.b) + "\","
                        + (segment.gap ? "1" : "0") + "," + (segment.reversed ? "1" : "0") + "," + (segment.secondHomolog ? "1" : "0") + ","
                        + jsNumber(segment.centromereStartMbp) + "," + jsNumber(segment.centromereEndMbp) + ","
                        + jsNumber(segment.originStartMbp) + "," + jsNumber(segment.originEndMbp) + "," + jsString(segment.originName) + "],";
                }

                js += "],\n          chevrons: [";
                for (const auto& [startMbp, endMbp] : piece.chevrons)
                {
                    js += "[" + jsNumber(startMbp) + "," + jsNumber(endMbp) + "],";
                }

                js += "], arrows: [";
                for (const auto& [fromMbp, toMbp] : piece.movedArrows)
                {
                    js += "[" + jsNumber(fromMbp) + "," + jsNumber(toMbp) + "],";
                }

                js += "], dashed: [";
                for (int boundary : piece.dashedBoundaries)
                {
                    js += std::to_string(boundary) + ",";
                }

                js += "], hidden: [";
                for (int boundary : piece.hiddenBoundaries)
                {
                    js += std::to_string(boundary) + ",";
                }

                js += "], markers: [";
                for (const KaryogramDetailMarker& marker : piece.markers)
                {
                    js += "[" + jsNumber(marker.positionMbp) + "," + jsString(marker.text) + "," + jsNumber(marker.deletedMbp) + "," + (marker.inserted ? "1" : "0") + "],";
                }

                js += "], damages: [";
                for (const KaryogramDetailDamage& damage : piece.damages)		// [position, DSB, number of SSBs]
                {
                    js += "[" + jsNumber(damage.positionMbp) + "," + (damage.doubleStrandBreak ? "1" : "0") + "," + std::to_string(damage.singleStrandBreaks) + "],";
                }

                js += "] },";
            }

            js += " ] },";
        }

        return js + " ]";
    }

#if !defined(_WIN32) && !defined(__APPLE__)
    // Detect Windows Subsystem for Linux (WSL).
    // Why: a program built inside WSL is a Linux program, so it would normally use xdg-open. But WSL usually has no
    // Linux browser installed, so xdg-open fails silently and nothing appears. In WSL we instead hand the file to the
    // Windows browser. WSL sets WSL_DISTRO_NAME, and its kernel version string mentions "microsoft" as a backup check.
    bool isRunningInWSL()
    {
        if (std::getenv("WSL_DISTRO_NAME") != nullptr)
        {
            return true;
        }

        std::ifstream version("/proc/version");
        std::string line;
        std::getline(version, line);
        return toLower(line).find("microsoft") != std::string::npos;
    }

    // Check that WSL can launch Windows programs.
    // Why: "Windows interop" can be switched off in WSL (common on locked-down work laptops). Then wslpath or
    // explorer.exe are missing, and we must fall back to telling the user where the file is.
    bool wslCanLaunchWindowsPrograms()
    {
        return std::system("command -v wslpath > /dev/null 2>&1 && command -v explorer.exe > /dev/null 2>&1") == 0;
    }
#endif

#if defined(__APPLE__)
    // Check whether this session is connected over SSH (only needed on macOS, Linux checks for a display instead).
    // Why: on a remote server, HPC cluster or container reached over SSH, the browser would have to open on the
    // remote machine, which has no screen. SSH sets these variables in every remote session.
    bool isRemoteSSHsession()
    {
        return std::getenv("SSH_CONNECTION") != nullptr || std::getenv("SSH_TTY") != nullptr;
    }
#endif

#if !defined(_WIN32)
    // Quote a path for a POSIX shell (Linux, macOS, WSL).
    // Why: the path goes into a std::system() command line. Input file names can contain spaces, $, quotes or other
    // characters the shell would otherwise interpret, which could break the command or run something unintended.
    // Inside single quotes the shell treats everything literally; an embedded ' is written as '\'' (close, escaped
    // quote, reopen).
    std::string shellQuote(const std::string& path)
    {
        std::string quoted = "'";
        for (char c : path)
        {
            quoted += (c == '\'') ? std::string("'\\''") : std::string(1, c);
        }
        quoted += "'";
        return quoted;
    }
#endif
}


bool writeKaryogramViewer(
    const std::vector<KaryogramImage>& images,
    const std::string& title,
    const std::filesystem::path& htmlPath)
{
    if (images.empty())
    {
        std::cerr << "ERROR: No karyograms to show in the web viewer.\n";
        return false;
    }

    std::ofstream html(htmlPath, std::ios::binary);

    if (!html)
    {
        std::cerr << "ERROR: Failed to create web viewer file " << htmlPath << "\n";
        return false;
    }

    const std::string safeTitle = htmlEscape(title);

    // ------------------------------------------
    // Page header and styling
    // ------------------------------------------
    // The karyogram sits inside a scrollable "viewport". Zooming changes the image's real width and height (not a CSS
    // scale transform), so the browser redraws the SVG at the new size and lines and text stay sharp at every zoom level.
    html << "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
         << "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
         << "<title>" << safeTitle << " karyograms</title>\n"
         << R"HTML(<style>
html, body { height: 100%; margin: 0; }
body { display: flex; flex-direction: column; font-family: sans-serif; background: #e9e9e9; color: #222; }
#titlebar { display: grid; grid-template-columns: 1fr auto 1fr; align-items: center; gap: 12px; padding: 6px 16px 0; background: #fff; }
#version { font-size: 13px; color: #666; white-space: nowrap; }
#titlebar h1 { font-size: 16px; margin: 0; text-align: center; overflow-wrap: anywhere; }
header { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; padding: 6px 16px 8px; background: #fff; border-bottom: 1px solid #ccc; }
.group { display: flex; gap: 4px; align-items: center; }
.spacer { flex: 1; }
button, select { font-size: 14px; padding: 4px 10px; cursor: pointer; }
#counter, #zoomLevel { color: #555; font-size: 14px; min-width: 52px; text-align: center; }
#viewport { flex: 1; overflow: auto; position: relative; cursor: grab; user-select: none; touch-action: none; }
#viewport.dragging { cursor: grabbing; }
#stage { position: relative; }
#frame { position: absolute; background: #fff; box-shadow: 0 1px 6px rgba(0,0,0,0.25); }
#karyogram { display: block; width: 100%; height: 100%; }
.hotspot { position: absolute; cursor: pointer; border-radius: 3px; }
.hotspot:hover { background: rgba(30,110,220,0.15); outline: 1px solid rgba(30,110,220,0.6); }
#popup { position: fixed; z-index: 10; min-width: 200px; max-width: 320px; background: #fff; border: 1px solid #bbb; border-radius: 6px; box-shadow: 0 4px 16px rgba(0,0,0,0.25); padding: 10px 14px; font-size: 14px; }
#popup[hidden] { display: none; }
#popup h2 { font-size: 15px; margin: 0 24px 6px 0; }
#popup h3 { font-size: 13px; margin: 10px 0 4px; color: #444; }
#popup table { border-collapse: collapse; width: 100%; }
#popup td { padding: 2px 0; }
#popup td.count { text-align: right; padding-left: 16px; font-variant-numeric: tabular-nums; }
#popup .none { color: #666; font-style: italic; }
#popup .note { margin-top: 8px; font-size: 11px; color: #777; }
#popupClose { position: absolute; top: 4px; right: 6px; border: none; background: none; font-size: 18px; padding: 0 4px; }
#hint { padding: 4px 16px; font-size: 12px; color: #666; background: #fff; border-top: 1px solid #ccc; }
.piecespot { position: absolute; cursor: zoom-in; border-radius: 4px; }
.piecespot:hover { outline: 2px solid rgba(30,110,220,0.7); outline-offset: 2px; }
#detail { position: fixed; inset: 0; z-index: 20; background: rgba(0,0,0,0.35); display: flex; justify-content: center; padding: 20px; }
#detail[hidden] { display: none; }
#detailPanel { flex: 0 1 auto; width: 100%; min-width: 0; display: flex; flex-direction: column; background: #fff; border-radius: 8px; box-shadow: 0 6px 24px rgba(0,0,0,0.3); overflow: hidden; }
#detailBar { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; padding: 8px 12px; border-bottom: 1px solid #ccc; }
#detailTitle { font-size: 15px; font-weight: bold; }
#detailCell, #detailScale { color: #555; font-size: 13px; }
#detailScale { min-width: 120px; text-align: right; font-variant-numeric: tabular-nums; }
#detailView { flex: 1; min-height: 0; overflow-x: auto; overflow-y: hidden; cursor: grab; user-select: none; touch-action: none; }
#detailView.dragging { cursor: grabbing; }
#detailSvg { display: block; }
#detailSvg text { font-family: sans-serif; font-size: 12px; fill: #333; }
#detailSvg .tick { font-size: 11px; fill: #555; font-variant-numeric: tabular-nums; }
#detailSvg .caption { font-size: 12px; fill: #555; }
#detailSvg .marker { cursor: pointer; }
#detailSvg .marker text { fill: #1e6edc; }
#detailSvg .inverted { font-size: 11px; fill: #333; }
#detailTip { position: fixed; z-index: 30; pointer-events: none; background: #fff; border: 1px solid #bbb; border-radius: 6px; box-shadow: 0 4px 16px rgba(0,0,0,0.25); padding: 6px 10px; font-size: 13px; line-height: 1.45; max-width: 280px; }
#detailTip[hidden] { display: none; }
#detailTip .site { font-weight: bold; font-variant-numeric: tabular-nums; }
#detailTip .site + div, #detailTip .more { color: #333; }
#detailTip .tipsite + .tipsite { margin-top: 6px; }
#detailTip .head { color: #555; font-style: italic; margin-bottom: 4px; }
#detailHint { padding: 4px 12px; font-size: 12px; color: #666; border-top: 1px solid #ccc; }
</style>
</head>
<body>
)HTML";

    // ------------------------------------------
    // Toolbar: switch between karyograms, zoom controls
    // ------------------------------------------
    // Title bar: KaryoKE version on the left, input file name in the middle.
    html << "<div id=\"titlebar\"><span id=\"version\">" << htmlEscape(KARYOKE_VERSION_TEXT) << "</span><h1>" << safeTitle << "</h1><span></span></div>\n"
         << "<header>\n"
         << "<div class=\"group\">\n"
         << "<button id=\"prev\" title=\"Previous karyogram\">&larr; Previous</button>\n"
         << "<select id=\"cell\" title=\"Choose a karyogram\">\n";

    for (size_t i = 0; i < images.size(); i++)
    {
        html << "<option value=\"" << i << "\">" << htmlEscape(images[i].label) << "</option>\n";
    }

    html << R"HTML(</select>
<button id="next" title="Next karyogram">Next &rarr;</button>
<span id="counter"></span>
</div>
<div class="spacer"></div>
<div class="group">
<button id="zoomOut" title="Zoom out (-)">&minus;</button>
<span id="zoomLevel"></span>
<button id="zoomIn" title="Zoom in (+)">+</button>
<button id="fit" title="Fit the whole karyogram in the window (0)">Fit</button>
<button id="actual" title="Show at original size (1)">100%</button>
</div>
<div class="group">
<select id="format" title="File format for downloads">
<option value="png">PNG</option>
<option value="svg">SVG</option>
</select>
<button id="download" title="Download the karyogram shown">Download</button>
<button id="downloadAll" title="Download every karyogram in one .zip file">Download all (.zip)</button>
</div>
</header>
<div id="viewport"><div id="stage"><div id="frame"><img id="karyogram" alt="Karyogram" draggable="false"><div id="hotspots"></div></div></div></div>
<div id="popup" hidden role="dialog"><button id="popupClose" title="Close (Esc)">&times;</button><div id="popupBody"></div></div>
<div id="detail" hidden role="dialog" aria-label="Zoomed homolog"><div id="detailPanel">
<div id="detailBar"><span id="detailTitle"></span><span id="detailCell"></span><div class="spacer"></div>
<div class="group"><button id="detailOut" title="Zoom out (-)">&minus;</button><button id="detailIn" title="Zoom in (+)">+</button><button id="detailWhole" title="Show the whole homolog (0)">Whole homolog</button></div>
<span id="detailScale"></span><button id="detailClose" title="Back to the karyogram (Esc)">Close &times;</button></div>
<div id="detailView"><svg id="detailSvg" xmlns="http://www.w3.org/2000/svg"></svg></div>
<div id="detailHint"></div>
</div><div id="detailTip" hidden></div></div>
<div id="hint">Click a chromosome for a zoomed view &middot; click a chromosome label for its summary &middot; scroll or pinch to zoom at the cursor &middot; drag to pan &middot; double-click to zoom in &middot; keys: + &minus; 0 (fit) 1 (100%)</div>
)HTML";

    // ------------------------------------------
    // Embedded karyogram images
    // ------------------------------------------
    // Each image is embedded as a base64 data: URL so the HTML file works on its own, with no other files or server.
    // The SVG is shown and both formats are kept for the download buttons.
    html << "<script>\nconst TITLE = " << jsString(title) << ";\nconst karyograms = [\n";

    for (const KaryogramImage& image : images)
    {
        html << "  { filename: " << jsString(image.filename)
             << ",\n    svg: \"" << (image.svgBytes.empty() ? "" : "data:image/svg+xml;base64," + base64Encode(image.svgBytes))
             << "\",\n    png: \"" << (image.pngBytes.empty() ? "" : "data:image/png;base64," + base64Encode(image.pngBytes))
             << "\",\n    hotspots: [";

        // Clickable chromosome labels: position in image pixels plus the counts shown in the popup.
        for (const KaryogramHotspot& hotspot : image.hotspots)
        {
            html << "\n      { title: " << jsString(hotspot.title)
                 << ", x: " << hotspot.x << ", y: " << hotspot.y
                 << ", w: " << hotspot.width << ", h: " << hotspot.height << ", sections: [";
            for (const KaryogramPopupSection& section : hotspot.sections)
            {
                html << "{ heading: " << jsString(section.heading) << ", counts: [";
                for (const auto& [type, count] : section.counts)
                {
                    html << "[" << jsString(type) << ", " << count << "], ";
                }
                html << "] }, ";
            }
            html << "] },";
        }

        html << " ],\n    popupNote: " << jsString(image.popupNote)
             << ",\n    popupEmptyText: " << jsString(image.popupEmptyText)
             << ",\n    detailSlots: " << detailSlotsJs(image.detailSlots) << " },\n";
    }

    // ------------------------------------------
    // Switching, zoom and pan logic
    // ------------------------------------------
    html << R"HTML(];

const viewport = document.getElementById('viewport');
const stage = document.getElementById('stage');
const frame = document.getElementById('frame');
const img = document.getElementById('karyogram');
const hotspotLayer = document.getElementById('hotspots');
const popup = document.getElementById('popup');
const popupBody = document.getElementById('popupBody');
const select = document.getElementById('cell');
const counter = document.getElementById('counter');
const zoomLevel = document.getElementById('zoomLevel');

const MIN_SCALE = 0.1, MAX_SCALE = 20, STEP = 1.25, MARGIN = 16;
let current = 0;
let scale = 1;
let fitMode = true;            // While true the image refits when the window is resized.

// Scale at which the whole karyogram fits inside the window.
function fitScale() {
  const w = img.naturalWidth || 1, h = img.naturalHeight || 1;
  return Math.min((viewport.clientWidth - 2 * MARGIN) / w, (viewport.clientHeight - 2 * MARGIN) / h);
}

// Offset of the image inside the stage: centred when smaller than the window, otherwise a small margin.
function offsets(s) {
  const w = img.naturalWidth * s, h = img.naturalHeight * s;
  return {
    x: Math.max(MARGIN, (viewport.clientWidth - w) / 2),
    y: Math.max(MARGIN, (viewport.clientHeight - h) / 2),
  };
}

// Apply a new zoom level while keeping the image point under (screenX, screenY) fixed on screen.
function setScale(newScale, screenX, screenY) {
  newScale = Math.min(MAX_SCALE, Math.max(MIN_SCALE, newScale));
  if (screenX === undefined) { screenX = viewport.clientWidth / 2; screenY = viewport.clientHeight / 2; }

  const before = offsets(scale);
  const imageX = (viewport.scrollLeft + screenX - before.x) / scale;    // Point under the cursor, in image pixels.
  const imageY = (viewport.scrollTop + screenY - before.y) / scale;

  scale = newScale;
  const after = offsets(scale);
  const w = img.naturalWidth * scale, h = img.naturalHeight * scale;
  frame.style.width = w + 'px';            // Real layout size, so the SVG is redrawn sharply at this size.
  frame.style.height = h + 'px';           // The clickable labels inside the frame are placed in %, so they follow.
  frame.style.left = after.x + 'px';
  frame.style.top = after.y + 'px';
  stage.style.width = Math.max(viewport.clientWidth, w + 2 * after.x) + 'px';
  stage.style.height = Math.max(viewport.clientHeight, h + 2 * after.y) + 'px';

  viewport.scrollLeft = imageX * scale + after.x - screenX;
  viewport.scrollTop = imageY * scale + after.y - screenY;
  zoomLevel.textContent = Math.round(scale * 100) + '%';
}

function fit() { fitMode = true; setScale(fitScale()); }
function zoomBy(factor, x, y) { fitMode = false; setScale(scale * factor, x, y); }

function show(index) {
  current = (index + karyograms.length) % karyograms.length;
  const first = !img.src;
  closePopup();
  img.onload = () => {
    placeHotspots();
    if (first || fitMode) fit(); else setScale(scale);           // Keep the zoom when switching cells, to compare them.
  };
  img.src = karyograms[current].svg || karyograms[current].png;      // Prefer the sharp vector image.
  select.value = current;
  counter.textContent = (current + 1) + ' of ' + karyograms.length;
}

// Mouse wheel and trackpad pinch (sent as wheel events with ctrlKey) zoom around the cursor.
viewport.addEventListener('wheel', (e) => {
  e.preventDefault();
  const r = viewport.getBoundingClientRect();
  const delta = e.deltaMode === 1 ? e.deltaY * 33 : e.deltaY;           // Line-based scrolling (Firefox) to pixels.
  zoomBy(Math.exp(-delta * (e.ctrlKey ? 0.01 : 0.002)), e.clientX - r.left, e.clientY - r.top);
}, { passive: false });

// Drag with the mouse (or one finger) to pan.
let drag = null;
viewport.addEventListener('pointerdown', (e) => {
  if (e.button !== 0) return;
  drag = { x: e.clientX, y: e.clientY, left: viewport.scrollLeft, top: viewport.scrollTop, moved: false };
  viewport.setPointerCapture(e.pointerId);
  viewport.classList.add('dragging');
});
viewport.addEventListener('pointermove', (e) => {
  if (!drag) return;
  if (Math.abs(e.clientX - drag.x) + Math.abs(e.clientY - drag.y) > 4) drag.moved = true;
  viewport.scrollLeft = drag.left - (e.clientX - drag.x);
  viewport.scrollTop = drag.top - (e.clientY - drag.y);
});
const endDrag = () => { drag = null; viewport.classList.remove('dragging'); };
viewport.addEventListener('pointerup', (e) => {
  // A press and release without dragging is a click. The pointer is captured by the viewport for panning, so the
  // click target is looked up from the position: a chromosome label opens its popup, anywhere else closes it.
  if (drag && !drag.moved) {
    const spot = hotspotAt(e.clientX, e.clientY);
    const slot = spot ? -1 : pieceSlotAt(e.clientX, e.clientY);
    if (spot) openPopup(spot, e.clientX, e.clientY);
    else if (slot >= 0) openDetail(slot);
    else closePopup();
  }
  endDrag();
});
viewport.addEventListener('pointercancel', endDrag);

viewport.addEventListener('dblclick', (e) => {
  if (hotspotAt(e.clientX, e.clientY) || pieceSlotAt(e.clientX, e.clientY) >= 0) return;   // Double-clicking a label or chromosome should not also zoom.
  const r = viewport.getBoundingClientRect();
  zoomBy(2, e.clientX - r.left, e.clientY - r.top);
});

document.getElementById('zoomIn').onclick = () => zoomBy(STEP);
document.getElementById('zoomOut').onclick = () => zoomBy(1 / STEP);
document.getElementById('fit').onclick = fit;
document.getElementById('actual').onclick = () => { fitMode = false; setScale(1); };
document.getElementById('prev').onclick = () => show(current - 1);
document.getElementById('next').onclick = () => show(current + 1);
select.onchange = () => show(Number(select.value));

document.addEventListener('keydown', (e) => {
  if (!detail.hidden) { detailKey(e); return; }                      // Keys act on the zoomed homolog while it is open.
  if (e.key === 'Escape') { closePopup(); return; }
  if (e.target === select || e.ctrlKey || e.metaKey || e.altKey) return;
  if (e.key === '+' || e.key === '=') zoomBy(STEP);
  else if (e.key === '-' || e.key === '_') zoomBy(1 / STEP);
  else if (e.key === '0') fit();
  else if (e.key === '1') { fitMode = false; setScale(1); }
});

// ------------------------------------------------------------
// Clickable chromosome labels and the summary popup (mutations for SDR, damages for SDD)
// ------------------------------------------------------------
// One invisible box per chromosome label, placed in % of the image size so it stays on the label at every zoom level.
function placeHotspots() {
  hotspotLayer.replaceChildren();
  const w = img.naturalWidth || 1, h = img.naturalHeight || 1;
  karyograms[current].hotspots.forEach((spot, i) => {
    const box = document.createElement('div');
    box.className = 'hotspot';
    box.dataset.index = i;
    box.title = spot.title + ': click for summary';
    box.style.left = (100 * spot.x / w) + '%';
    box.style.top = (100 * spot.y / h) + '%';
    box.style.width = (100 * spot.w / w) + '%';
    box.style.height = (100 * spot.h / h) + '%';
    hotspotLayer.appendChild(box);
  });
  // One box per drawn piece of a homolog (chromosome, excised piece, ecDNA or ring), opening that homolog's zoomed view. The two sister
  // chromatids of a replicated SDD homolog get one box around both, so the pair is outlined as one chromosome.
  (karyograms[current].detailSlots || []).forEach((slot, i) => {
    const areas = slot.sisters ? [slot.pieces.reduce((a, p) => {
      const x = Math.min(a.x, p.x), y = Math.min(a.y, p.y);
      return { x, y, w: Math.max(a.x + a.w, p.x + p.w) - x, h: Math.max(a.y + a.h, p.y + p.h) - y };
    })] : slot.pieces;
    for (const piece of areas) {
      const box = document.createElement('div');
      box.className = 'piecespot';
      box.dataset.slot = i;
      box.title = slot.title + ': click for a zoomed view';
      box.style.left = (100 * piece.x / w) + '%';
      box.style.top = (100 * piece.y / h) + '%';
      box.style.width = (100 * piece.w / w) + '%';
      box.style.height = (100 * piece.h / h) + '%';
      hotspotLayer.appendChild(box);
    }
  });
}

// Index of the homolog slot whose piece is under the screen point, or -1.
function pieceSlotAt(x, y) {
  const box = document.elementFromPoint(x, y)?.closest('.piecespot');
  return box ? Number(box.dataset.slot) : -1;
}

function hotspotAt(x, y) {
  const box = document.elementFromPoint(x, y)?.closest('.hotspot');
  return box ? karyograms[current].hotspots[Number(box.dataset.index)] : null;
}

function addElement(parent, tag, text, className) {
  const element = document.createElement(tag);
  if (text !== undefined) element.textContent = text;            // textContent, never innerHTML: names are shown as plain text.
  if (className) element.className = className;
  parent.appendChild(element);
  return element;
}

function openPopup(spot, x, y) {
  popupBody.replaceChildren();
  addElement(popupBody, 'h2', spot.title);
  for (const section of spot.sections) {
    if (section.heading) addElement(popupBody, 'h3', section.heading);
    if (section.counts.length === 0) { addElement(popupBody, 'div', karyograms[current].popupEmptyText, 'none'); continue; }
    const table = addElement(popupBody, 'table');
    for (const [type, count] of section.counts) {
      const row = addElement(table, 'tr');
      addElement(row, 'td', type);
      addElement(row, 'td', count, 'count');
    }
  }
  if (karyograms[current].popupNote) addElement(popupBody, 'div', karyograms[current].popupNote, 'note');

  popup.hidden = false;
  // Open next to the click, but keep the whole popup inside the window.
  const r = popup.getBoundingClientRect();
  popup.style.left = Math.max(8, Math.min(x + 12, window.innerWidth - r.width - 8)) + 'px';
  popup.style.top = Math.max(8, Math.min(y + 12, window.innerHeight - r.height - 8)) + 'px';
}

function closePopup() { popup.hidden = true; }
document.getElementById('popupClose').onclick = closePopup;

// ------------------------------------------------------------
// Zoomed homolog view: one homolog redrawn at its true size, with a base-pair ruler, zoomable down to single bases
// ------------------------------------------------------------
const detail = document.getElementById('detail');
const detailView = document.getElementById('detailView');
const detailSvg = document.getElementById('detailSvg');
const detailScale = document.getElementById('detailScale');
const SVG_NS = 'http://www.w3.org/2000/svg';
const D_MIN_MBP_PER_PX = 0.000001 / 12;     // Deepest zoom: 12 px per base.
const D_TOP = 44, D_BOTTOM = 24;            // Room above the pieces (captions) and below them at "whole homolog".
const D_WIDTH = 56;                         // Drawn width of a piece; the karyogram's details (pinch, dashes) are scaled from its 20 px.
const D_K = D_WIDTH / 20;
let dSlot = null, dTop = 0, dMbpPerPx = 1, dDrag = null, dClipId = 0, dColumns = [];
const detailTip = document.getElementById('detailTip');

// Opens the zoomed view of one homolog slot of the karyogram shown, fitted to the window. Alters the view state and shows the overlay.
function openDetail(slotIndex) {
  dSlot = karyograms[current].detailSlots[slotIndex];
  if (!dSlot) return;
  closePopup();
  document.getElementById('detailTitle').textContent = dSlot.title;
  document.getElementById('detailCell').textContent = select.options[current].text;
  const sdd = dSlot.pieces.some((p) => p.damages.length);
  document.getElementById('detailHint').textContent = 'Scroll or pinch to zoom at the cursor · drag to move along the homolog · '
    + (sdd ? 'point at (or tap) a damage site for its breaks · the ruler gives the position from the p-arm end'
           : 'point at any stretch of a chromosome for the chromosome it came from and its span · click a blue pointer to zoom to that indel · ecDNA and rings are drawn as circles sized by their length, not to the ruler scale · the ruler gives the position along each piece, from its top end');
  detail.hidden = false;
  detailWhole();
}

// Closes the zoomed view and returns to the karyogram.
function closeDetail() { detail.hidden = true; dSlot = null; detailSvg.replaceChildren(); hideDamageTip(); }

// Length of the longest linear piece of the open slot, in Mbp: the ruler runs from 0 to here. Circles (ecDNA, rings) are not on the ruler.
function detailLength() {
  const linear = dSlot.pieces.filter((p) => !p.circular);
  return Math.max(1e-6, ...(linear.length ? linear : dSlot.pieces).map((p) => p.len));
}

// Screen diameter of a circular piece: a 40 px minimum, so its colours always show, plus its length in pixels of the ruler at
// "whole homolog" (at most 200 px in all). Fixed while zooming, since the circle is not drawn to the ruler's scale.
function circleDiameter(piece) {
  const whole = Math.max(40, detailView.clientHeight - D_TOP - D_BOTTOM);
  return Math.min(200, 40 + whole * piece.len / detailLength());
}

// Original spans [from, to] of the material of a piece, in Mbp: stretches of one colour that continue each other (within a base) are joined.
function pieceSpans(piece) {
  const spans = [];
  for (const s of piece.segs) {
    if (s[3] || s[8] < 0) continue;
    const last = spans[spans.length - 1];
    if (last && last.color === s[2] && Math.abs(s[8] - last.to) <= 1.5e-6) last.to = s[9];
    else if (last && last.color === s[2] && Math.abs(s[9] - last.from) <= 1.5e-6) last.from = s[8];
    else spans.push({ from: s[8], to: s[9], color: s[2] });
  }
  return spans;
}

// Label lines of a circular piece: "ecDNA: X Mbp – Y Mbp" (or "Ring: ..."), the original span of its material; further spans of a circle
// made of several pieces follow on lines of their own (at most four, then how many more).
function circleCaption(piece) {
  const mbp = (value) => String(Number(value.toFixed(6))) + ' Mbp';
  const spans = pieceSpans(piece).map((s) => `${mbp(s.from)} – ${mbp(s.to)}`);
  const kind = /^Ring/.test(piece.caption) ? 'Ring' : 'ecDNA';
  if (!spans.length) return [kind];
  const lines = [`${kind}: ${spans[0]}`, ...spans.slice(1, 4)];
  if (spans.length > 4) lines.push(`and ${spans.length - 4} more`);
  return lines;
}

// Draws a circular piece (ecDNA or ring) as the karyogram does, as a ring at a fixed screen place below the captions: a coloured band with a
// white centre half its diameter, one stretch of band per segment in its colour (clockwise from the top), white for deleted bases, a thin line
// where material of two chromosomes meets, and inner and outer borders that are dashed along material from a second homolog. A ring chromosome
// is turned so its centromere sits at 12 o'clock, where both circles pinch towards each other (the karyogram's centromere constriction, scaled
// like the linear pieces' pinch) and a grey ellipse crosses the band. Alters the SVG only.
function drawDetailCircle(svg, piece, column) {
  const d = column.diameter, r = d / 2, inner = d / 4, cx = column.x + r, cy = column.y + r, len = Math.max(1e-12, piece.len);
  const cen = piece.segs.find((s) => s[6] >= 0);
  const turn = cen ? (cen[6] + cen[7]) / 2 : 0;                      // Material at the top of the ring
  const angle = (mbp) => 2 * Math.PI * (mbp - turn) / len;
  const point = (mbp, radius) => [cx + radius * Math.sin(angle(mbp)), cy - radius * Math.cos(angle(mbp))];
  const full = (a, b) => b - a >= len - 1e-9;
  // One stretch of the band, from startMbp to endMbp clockwise, between the inner and outer circle.
  const band = (startMbp, endMbp, fill) => {
    if (full(startMbp, endMbp)) {
      svgElement(body, 'path', { d: `M${cx - r},${cy} a${r},${r} 0 1 0 ${2 * r},0 a${r},${r} 0 1 0 ${-2 * r},0 Z M${cx - inner},${cy} a${inner},${inner} 0 1 1 ${2 * inner},0 a${inner},${inner} 0 1 1 ${-2 * inner},0 Z`,
        fill, 'fill-rule': 'evenodd' });
      return;
    }
    const large = (endMbp - startMbp) / len > 0.5 ? 1 : 0;
    const [ox0, oy0] = point(startMbp, r), [ox1, oy1] = point(endMbp, r), [ix0, iy0] = point(startMbp, inner), [ix1, iy1] = point(endMbp, inner);
    svgElement(body, 'path', { d: `M${ox0},${oy0} A${r},${r} 0 ${large} 1 ${ox1},${oy1} L${ix1},${iy1} A${inner},${inner} 0 ${large} 0 ${ix0},${iy0} Z`,
      fill, stroke: fill, 'stroke-width': 0.6 });
  };
  // One arc of a border circle, solid or dashed.
  const arc = (startMbp, endMbp, radius, dashed) => {
    const style = { fill: 'none', stroke: '#1a1a1a', 'stroke-width': 1.5, 'stroke-dasharray': dashed ? '6 4' : 'none' };
    if (full(startMbp, endMbp)) { svgElement(borders, 'circle', { cx, cy, r: radius, ...style }); return; }
    const [x0, y0] = point(startMbp, radius), [x1, y1] = point(endMbp, radius), large = (endMbp - startMbp) / len > 0.5 ? 1 : 0;
    svgElement(borders, 'path', { d: `M${x0},${y0} A${radius},${radius} 0 ${large} 1 ${x1},${y1}`, ...style });
  };

  // Centromere constriction of a ring: the outer circle pinches inward and the inner one outward at 12 o'clock. The band is clipped to the
  // pinched shape, and the borders leave out the top so the pinch can be drawn in their place.
  // Pinch depth and half its angle: at least the linear pieces' pinch, grown with the band and kept wider than the centromere ellipse so it shows.
  const cenWidth = cen ? Math.min(Math.max(5, 2 * Math.PI * r * (cen[7] - cen[6]) / len), 0.6 * (r - inner)) : 0;
  const depth = Math.max(2 * D_K, 0.12 * (r - inner)), notch = Math.max(3.5 * D_K, cenWidth / 2 + 2 * depth) / r;
  const at = (theta, radius) => [cx + radius * Math.sin(theta), cy - radius * Math.cos(theta)];
  const pinched = (radius, apexRadius) => {
    const [x0, y0] = at(notch, radius), [x1, y1] = at(-notch, radius);
    return `M${x0},${y0} A${radius},${radius} 0 1 1 ${x1},${y1} L${cx},${cy - apexRadius} Z`;
  };
  let body = svg, borders = svg;
  if (cen) {
    const bandClip = 'dclip' + (++dClipId), borderClip = 'dclip' + (++dClipId);
    svgElement(svgElement(svg, 'clipPath', { id: bandClip }), 'path', { d: pinched(r, r - depth) + ' ' + pinched(inner, inner + depth), 'clip-rule': 'evenodd' });
    const [bx0, by0] = at(notch, r + 4), [bx1, by1] = at(-notch, r + 4);
    svgElement(svgElement(svg, 'clipPath', { id: borderClip }), 'path', { d: `M${cx},${cy} L${bx0},${by0} A${r + 4},${r + 4} 0 1 1 ${bx1},${by1} Z` });
    body = svgElement(svg, 'g', { 'clip-path': `url(#${bandClip})` });
    borders = svgElement(svg, 'g', { 'clip-path': `url(#${borderClip})` });
  }

  const minMbp = len / 360;                                          // A deleted stretch is at least one degree wide, so it stays visible
  for (const s of piece.segs) if (!s[3]) band(s[0], s[1], s[2]);
  const gaps = piece.segs.filter((s) => s[3]).map((s) => [s[0], s[1]]).concat(piece.markers.filter((m) => m[2] > 0).map((m) => [m[0], m[0] + m[2]]));
  for (const [a, b] of gaps) { const c = (a + b) / 2, half = Math.max(b - a, minMbp) / 2; band(c - half, c + half, '#ffffff'); }
  svgElement(body, 'circle', { cx, cy, r: inner, fill: '#ffffff' });

  // Thin line across the band where the material of two chromosomes meets, as between the wedges of the karyogram's multi-chromosome circles.
  const solid = piece.segs.filter((s) => !s[3]);
  for (let k = 0; k < solid.length && solid.length > 1; k++) {
    const a = solid[k], b = solid[(k + 1) % solid.length];
    if (a[2] === b[2]) continue;
    const [ox, oy] = point(a[1] % len, r), [ix, iy] = point(a[1] % len, inner);
    svgElement(body, 'line', { x1: ox, y1: oy, x2: ix, y2: iy, stroke: '#1a1a1a', 'stroke-width': 0.75 });
  }

  // Inner and outer borders, dashed along material from a second homolog (a white gap takes the style of the material before it).
  let dashed = false;
  for (const s of piece.segs) {
    if (!s[3]) dashed = !!s[5];
    arc(s[0], s[1], r, dashed);
    arc(s[0], s[1], inner, dashed);
  }

  // Ring chromosome: the grey centromere ellipse across the band (clipped to the pinched band), then the pinch lines of both circles in the
  // style of the centromere's material.
  if (cen) {
    svgElement(body, 'ellipse', { cx, cy: cy - (r + inner) / 2, rx: cenWidth / 2, ry: (r - inner) / 2, fill: '#999999' });
    for (const [radius, apexRadius] of [[r, r - depth], [inner, inner + depth]]) {
      const [x0, y0] = at(notch, radius), [x1, y1] = at(-notch, radius);
      svgElement(svg, 'path', { d: `M${x1},${y1} L${cx},${cy - apexRadius} L${x0},${y0}`, fill: 'none', stroke: '#1a1a1a', 'stroke-width': 1.5,
        'stroke-dasharray': cen[5] ? '6 4' : 'none' });
    }
  }
  column.caption.forEach((line, i) => svgElement(svg, 'text', { x: column.x, y: column.y + d + 16 + 15 * i, class: i ? 'tick' : 'caption' }, line));
}

// Keeps the zoom between "whole homolog" and 12 px per base, and the view on the homolog. Alters dTop and dMbpPerPx.
function detailClamp() {
  const h = detailView.clientHeight, len = detailLength();
  const maxMbpPerPx = len / Math.max(40, h - D_TOP - D_BOTTOM);
  dMbpPerPx = Math.min(maxMbpPerPx, Math.max(D_MIN_MBP_PER_PX, dMbpPerPx));
  dTop = Math.min(len + D_BOTTOM * dMbpPerPx - h * dMbpPerPx, Math.max(-D_TOP * dMbpPerPx, dTop));
}

// Shows the whole homolog. Alters the view state and redraws.
function detailWhole() { dMbpPerPx = Infinity; dTop = -Infinity; detailRender(); }

// Zooms by factor around a screen y position (relative to the view), keeping the base under it in place. Alters the view state and redraws.
function detailZoomAt(factor, screenY) {
  if (screenY === undefined) screenY = detailView.clientHeight / 2;
  const at = dTop + screenY * dMbpPerPx;
  dMbpPerPx /= factor; detailClamp();
  dTop = at - screenY * dMbpPerPx;
  detailRender();
}

// Centres the view on a position and zooms in so about 400 bases fill the height. Alters the view state and redraws.
function detailJumpTo(mbp) {
  const h = detailView.clientHeight;
  dMbpPerPx = Math.max(D_MIN_MBP_PER_PX, 0.0004 / h);
  dTop = mbp - (h / 2) * dMbpPerPx;
  detailRender();
}

// Creates one SVG element with the given attributes (and text) inside parent.
function svgElement(parent, tag, attrs, text) {
  const node = document.createElementNS(SVG_NS, tag);
  for (const key in attrs) node.setAttribute(key, attrs[key]);
  if (text !== undefined) node.textContent = text;
  parent.appendChild(node);
  return node;
}

// Tick spacing in Mbp: 1, 2 or 5 times a power of ten bases, about one tick per 70 px.
function tickStep() {
  const bases = dMbpPerPx * 70 * 1e6;
  const power = Math.pow(10, Math.floor(Math.log10(bases)));
  const m = bases / power;
  return Math.max(1, (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * power) / 1e6;
}

// Formats a ruler position with the unit that suits the tick spacing: Mbp, kbp or single bases.
function formatPosition(mbp, step) {
  if (step >= 0.1) return mbp.toFixed(step >= 1 ? 0 : 1) + ' Mbp';
  if (step >= 0.0001) return (mbp * 1000).toLocaleString('en-US', { maximumFractionDigits: step >= 0.001 ? 0 : 1 }) + ' kbp';
  return Math.round(mbp * 1e6).toLocaleString('en-US') + ' bp';
}

// Formats a stretch of the piece, e.g. "50.0 – 70.0 Mbp", or in single bases when it is shorter than 1 Mbp.
function formatRange(startMbp, endMbp) {
  if (endMbp - startMbp >= 1) return `${startMbp.toFixed(1)} – ${endMbp.toFixed(1)} Mbp`;
  const bp = (mbp) => Math.round(mbp * 1e6).toLocaleString('en-US');
  return `${bp(startMbp)} – ${bp(endMbp)} bp`;
}

// Formats a length for the scale readout and labels, e.g. "1.4 Mbp", "12 kbp", "47 bp".
function formatLength(mbp) {
  const bp = mbp * 1e6;
  if (bp >= 1e6) return (bp / 1e6).toPrecision(3).replace(/\.?0+$/, '') + ' Mbp';
  if (bp >= 1e3) return (bp / 1e3).toPrecision(3).replace(/\.?0+$/, '') + ' kbp';
  return (bp >= 10 ? Math.round(bp) : bp.toPrecision(2).replace(/\.?0+$/, '')) + ' bp';
}

// Outline of a piece between screen y top and bottom: rounded or flat ends (an end off screen is flat) and the centromere pinch on both sides.
function piecePath(x, top, bottom, roundTop, roundBottom, pinches) {
  const w = D_WIDTH, r = w / 2, depth = 2 * D_K, half = 3.5 * D_K;
  const capTop = roundTop && bottom - top > w, capBottom = roundBottom && bottom - top > w;
  let d = capTop ? `M${x},${top + r} A${r},${r} 0 0 1 ${x + w},${top + r}` : `M${x},${top} L${x + w},${top}`;
  for (const cy of pinches) d += ` L${x + w},${cy - half} L${x + w - depth},${cy} L${x + w},${cy + half}`;
  d += capBottom ? ` L${x + w},${bottom - r} A${r},${r} 0 0 1 ${x},${bottom - r}` : ` L${x + w},${bottom} L${x},${bottom}`;
  for (const cy of [...pinches].reverse()) d += ` L${x},${cy + half} L${x + depth},${cy} L${x},${cy - half}`;
  return d + ' Z';
}

// Outline of one sister chromatid of a replicated SDD homolog, the zoomed form of the karyogram's X shape: at each centromere the outer side
// is pinched in deeply over a tall stretch with a smooth cosine profile, while the inner side (facing the sister) bulges by half the gap
// between them so the two sisters touch there. side is 'left' or 'right'; the ends are rounded or flat as in piecePath.
const D_SISTER_HALF = 12 * D_K, D_SISTER_OUTER = 7 * D_K;   // Half height and depth of a sister chromatid's outer pinch, in px
// How far a sister chromatid's sides are pinched at screen height py (1 at a centromere centre, 0 outside the pinch), with a cosine profile.
function sisterBump(py, pinches) {
  return Math.max(0, ...pinches.map((cy) => {
    const t = Math.abs(py - cy) / D_SISTER_HALF;
    return t >= 1 ? 0 : 0.5 * (1 + Math.cos(Math.PI * t));
  }));
}

function sisterPath(x, top, bottom, roundTop, roundBottom, pinches, side, gap) {
  const w = D_WIDTH, r = w / 2, half = D_SISTER_HALF, outer = D_SISTER_OUTER, inner = -gap / 2;
  const capTop = roundTop && bottom - top > w, capBottom = roundBottom && bottom - top > w;
  const from = capTop ? top + r : top, to = capBottom ? bottom - r : bottom;
  const bump = (py) => sisterBump(py, pinches);
  // Rows along the straight part: the two ends, then every 2 px inside each pinch.
  const rows = [from, to];
  for (const cy of pinches) for (let py = cy - half; py <= cy + half; py += 2) if (py > from && py < to) rows.push(py);
  rows.sort((a, b) => a - b);
  const rightCut = side === 'left' ? inner : outer, leftCut = side === 'left' ? outer : inner;
  let d = capTop ? `M${x},${from} A${r},${r} 0 0 1 ${x + w},${from}` : `M${x},${from} L${x + w},${from}`;
  for (const py of rows) d += ` L${x + w - rightCut * bump(py)},${py}`;
  d += capBottom ? ` A${r},${r} 0 0 1 ${x},${to}` : ` L${x},${to}`;
  for (const py of [...rows].reverse()) d += ` L${x + leftCut * bump(py)},${py}`;
  return d + ' Z';
}

// Centromeres in the zoomed view: a grey ellipse of half width rx centred on cx over each centromere, at least 5 px tall; a centromere far
// taller than the window is drawn as a band. cens holds screen tops and bottoms. Alters the SVG only.
function drawDetailCentromeres(svg, cens, cx, rx, h) {
  const far = D_WIDTH * 2, clampY = (v) => Math.max(-far, Math.min(h + far, v));
  for (const c of cens) {
    if (c.bottom < -far || c.top > h + far) continue;
    const height = Math.max(5, c.bottom - c.top), cy = (c.top + c.bottom) / 2;
    if (height > 4 * h) svgElement(svg, 'rect', { x: cx - rx, y: clampY(c.top), width: 2 * rx, height: clampY(c.bottom) - clampY(c.top), fill: '#999999' });
    else svgElement(svg, 'ellipse', { cx, cy, rx, ry: height / 2, fill: '#999999' });
  }
}

// SDD damage sites of a piece at column x, as on the karyogram: a black line past both edges for a DSB, then a white line for the SSBs that
// gets longer with their number. On a sister chromatid ({ side, gap }) a white line inside the outer pinch is moved towards the inner side,
// and shortened only if it is still wider than the chromatid there, so it never sticks out of the outline. Alters the SVG only.
function drawDetailDamages(svg, piece, x, h, sister) {
  const pinches = sister ? piece.segs.filter((s) => s[6] >= 0).map((s) => ((s[6] + s[7]) / 2 - dTop) / dMbpPerPx) : [];
  const y = (mbp) => (mbp - dTop) / dMbpPerPx;
  for (const [mbp, dsb] of piece.damages) {
    const dy = y(mbp);
    if (dsb && dy > -2 && dy < h + 2) svgElement(svg, 'line', { x1: x - 4 * D_K, y1: dy, x2: x + D_WIDTH + 4 * D_K, y2: dy, stroke: '#000', 'stroke-width': 1.5 });
  }
  for (const [mbp, , ssbs] of piece.damages) {
    const dy = y(mbp);
    if (ssbs <= 0 || dy < -2 || dy > h + 2) continue;
    let half = Math.min(D_WIDTH, (4 + 4 * (ssbs - 1)) * D_K) / 2, cx = x + D_WIDTH / 2;
    if (sister) {
      const cut = D_SISTER_OUTER * sisterBump(dy, pinches), margin = 1.5;
      const left = x + (sister.side === 'left' ? cut : 0) + margin, right = x + D_WIDTH - (sister.side === 'left' ? 0 : cut) - margin;
      half = Math.min(half, (right - left) / 2);
      cx = Math.max(left + half, Math.min(right - half, cx));
    }
    svgElement(svg, 'line', { x1: cx - half, y1: dy, x2: cx + half, y2: dy, stroke: '#fff', 'stroke-width': 2 });
  }
}

// Draws one piece of the slot at column x: fills, white gaps, boundaries, outline, centromeres, one chevron per inverted stretch,
// moved-piece arrows and indel pointers, all placed with the current zoom. A sister chromatid ({ side, gap }) gets the sister outline and
// leaves its centromere and damage sites to detailRender, which draws them once both sisters are drawn. Alters the SVG only.
function drawDetailPiece(svg, piece, x, h, sister) {
  const y = (mbp) => (mbp - dTop) / dMbpPerPx;
  const far = D_WIDTH * 2;
  const top = Math.max(-far, y(0)), bottom = Math.min(h + far, y(piece.len));
  if (bottom < -far || top > h + far) return;
  const roundTop = piece.roundTop && y(0) > -far, roundBottom = piece.roundBottom && y(piece.len) < h + far;
  const segs = piece.segs.map((s) => ({ start: s[0], end: s[1], color: s[2], gap: !!s[3], rev: !!s[4], second: !!s[5], cenStart: s[6], cenEnd: s[7] }));
  const cens = segs.filter((s) => s.cenStart >= 0).map((s) => ({ top: y(s.cenStart), bottom: y(s.cenEnd) }));
  const pinches = cens.map((c) => (c.top + c.bottom) / 2).filter((cy) => cy > top && cy < bottom);
  const outlinePath = sister ? sisterPath(x, top, bottom, roundTop, roundBottom, pinches, sister.side, sister.gap)
    : piecePath(x, top, bottom, roundTop, roundBottom, pinches);
  const bx = x - 6, bw = D_WIDTH + 12;          // Fills reach a little past the piece; the outline clip trims them (a sister's inner bulge included)
  const clipId = 'dclip' + (++dClipId);
  svgElement(svgElement(svg, 'clipPath', { id: clipId }), 'path', { d: outlinePath });
  const body = svgElement(svg, 'g', { 'clip-path': `url(#${clipId})` });
  const visible = (a, b) => b >= -far && a <= h + far;
  const clampY = (v) => Math.max(-far, Math.min(h + far, v));

  // Fills, then white gaps padded to at least 1 px so even a few deleted bases stay visible.
  for (const s of segs) {
    if (s.gap || !visible(y(s.start), y(s.end))) continue;
    svgElement(body, 'rect', { x: bx, y: clampY(y(s.start)), width: bw, height: Math.max(0.3, clampY(y(s.end)) - clampY(y(s.start))), fill: s.color });
  }
  // Indels are taken out of the karyogram drawing, so the piece runs through them: their deleted bases are drawn here as white gaps too,
  // and inserted bases as a blue line across the piece at the site.
  const gaps = segs.filter((s) => s.gap).map((s) => [s.start, s.end]).concat(piece.markers.filter((m) => m[2] > 0).map((m) => [m[0], m[0] + m[2]]));
  for (const [startMbp, endMbp] of gaps) {
    if (!visible(y(startMbp), y(endMbp))) continue;
    let a = clampY(y(startMbp)), b = clampY(y(endMbp));
    if (b - a < 1) { const c = (a + b) / 2; a = c - 0.5; b = c + 0.5; }
    svgElement(body, 'rect', { x: bx, y: a, width: bw, height: b - a, fill: '#ffffff' });
  }
  for (const [mbp, , , inserted] of piece.markers) {
    if (inserted && y(mbp) > -2 && y(mbp) < h + 2) svgElement(body, 'line', { x1: x, y1: y(mbp), x2: x + D_WIDTH, y2: y(mbp), stroke: '#1e6edc', 'stroke-width': 2 });
  }

  // One chevron per inverted stretch, stretched over the whole stretch at every zoom: apex at its top centre, arms to its bottom corners.
  // Only the part of each arm that is on screen is drawn, so very deep zooms never build huge coordinates.
  const labels = [];
  for (const [startMbp, endMbp] of piece.chevrons) {
    const yTop = y(startMbp), yBottom = y(endMbp);
    if (!visible(yTop, yBottom)) continue;
    const from = Math.max(yTop, -far), to = Math.min(yBottom, h + far), span = Math.max(1e-9, yBottom - yTop);
    const armX = (screenY) => (D_WIDTH / 2) * (screenY - yTop) / span;
    const cx = x + D_WIDTH / 2;
    svgElement(body, 'path', { d: `M${cx - armX(from)},${from} L${cx - armX(to)},${to} M${cx + armX(from)},${from} L${cx + armX(to)},${to}`,
      fill: 'none', stroke: '#000', 'stroke-width': 2, 'stroke-linejoin': 'round' });
    if (yTop < 26 && yBottom > 50) labels.push([Math.max(26, yTop) + 16, `Inverted ${formatRange(startMbp, endMbp)}`]);
  }

  // Boundaries between segments, by the karyogram's rules: none next to a white gap or an inverted piece, none where a moved piece was cut
  // out, dashed at a moved piece's new place.
  const dashedSet = new Set(piece.dashed), hiddenSet = new Set(piece.hidden);
  for (let i = 0; i + 1 < segs.length; i++) {
    const dashed = dashedSet.has(i + 1), a = segs[i], b = segs[i + 1];
    if ((a.gap || b.gap) && (!dashed || (a.gap && b.gap))) continue;
    if (!dashed && (a.rev || b.rev || hiddenSet.has(i + 1))) continue;
    const by = y(a.end);
    if (by < -2 || by > h + 2) continue;
    svgElement(body, 'line', { x1: x, y1: by, x2: x + D_WIDTH, y2: by, stroke: '#000', 'stroke-width': 1, 'stroke-dasharray': dashed ? '6 4' : 'none' });
  }

  // Outline: dashed along material from a second homolog (a white gap takes the style of the nearest piece above, else below), solid elsewhere.
  if (piece.outline) {
    const dashedRun = segs.map((s, i) => {
      if (!s.gap) return s.second;
      for (let j = i - 1; j >= 0; j--) if (!segs[j].gap) return segs[j].second;
      for (let j = i + 1; j < segs.length; j++) if (!segs[j].gap) return segs[j].second;
      return false;
    });
    for (let start = 0; start < segs.length;) {
      let end = start;
      while (end + 1 < segs.length && dashedRun[end + 1] === dashedRun[start]) end++;
      const runTop = start === 0 ? -1e4 : clampY(y(segs[start].start)), runBottom = end === segs.length - 1 ? h + 1e4 : clampY(y(segs[end].end));
      const runClip = 'dclip' + (++dClipId);
      svgElement(svgElement(svg, 'clipPath', { id: runClip }), 'rect', { x: x - 10, y: runTop, width: D_WIDTH + 20, height: Math.max(0, runBottom - runTop) });
      svgElement(svg, 'path', { d: outlinePath, fill: 'none', stroke: '#1a1a1a', 'stroke-width': 1.5, 'clip-path': `url(#${runClip})`,
        'stroke-dasharray': dashedRun[start] ? '8 6' : 'none' });
      start = end + 1;
    }
  }

  if (!sister) {
    drawDetailCentromeres(svg, cens, x + D_WIDTH / 2, D_WIDTH / 2 + 1, h);
    drawDetailDamages(svg, piece, x, h);
  }

  // Red hooked arrows of pieces moved within the chromosome, from where they were cut out to their middle, on the left of the piece.
  piece.arrows.forEach(([fromMbp, toMbp], i) => {
    const fromY = Math.max(-h, Math.min(2 * h, y(fromMbp))), toY = Math.max(-h, Math.min(2 * h, y(toMbp)));
    const tipX = x - 1, headBase = tipX - 9, bulge = Math.min(48, Math.max(16, 0.6 * Math.abs(toY - fromY))) + 6 * i;
    const g = svgElement(svg, 'g', { fill: '#bf382b', stroke: '#bf382b' });
    svgElement(g, 'circle', { cx: tipX, cy: fromY, r: 2.2, stroke: 'none' });
    svgElement(g, 'path', { d: `M${tipX},${fromY} C${x - bulge},${fromY} ${x - bulge},${toY} ${headBase + 1},${toY}`, fill: 'none', 'stroke-width': 1.8, 'stroke-linecap': 'round' });
    svgElement(g, 'path', { d: `M${tipX},${toY} L${headBase},${toY - 5} L${headBase + 2},${toY} L${headBase},${toY + 5} Z`, stroke: 'none' });
  });

  // Labels for inverted stretches whose chevron tip is above the window, so the inversion is still named.
  for (const [labelY, text] of labels) svgElement(svg, 'text', { x: x + D_WIDTH + 8, y: labelY, class: 'inverted' }, text);

  // Indels: never drawn on the karyogram, marked here with a pointer and their size; a click on one zooms to it.
  for (const [mbp, text] of piece.markers) {
    const my = y(mbp);
    if (my < -8 || my > h + 8) continue;
    const g = svgElement(svg, 'g', { class: 'marker', 'data-mbp': mbp });
    svgElement(g, 'path', { d: `M${x + D_WIDTH + 3},${my} l10,-6 v12 z`, fill: '#1e6edc' });
    svgElement(g, 'text', { x: x + D_WIDTH + 16, y: my + 4 }, text);
  }
}

const D_LEGEND_WIDTH = 240;                 // Width of the symbol legend left of the ruler.

// Draws the karyogram's symbol legend as a column in the zoomed view, fixed while zooming: for SDR chromosome, centromere, reversed segment,
// moved segment and ecDNA (as on the karyogram); for SDD chromosome, double-strand break, single-strand break and centromere. Symbols use the open homolog's colour
// and are drawn as in the zoomed view. Alters the SVG only.
function drawDetailLegend(svg, x, y) {
  const sdd = dSlot.pieces.some((p) => p.damages.length);
  const firstSeg = dSlot.pieces.flatMap((p) => p.segs).find((s) => !s[3]);
  const color = firstSeg ? firstSeg[2] : '#e69545';
  const row = 46, symbolX = x + 14, textX = x + 52;
  const items = sdd ? ['Chromosome', 'Double-Strand Break', 'Single-Strand Break', 'Centromere']
                    : ['Chromosome', 'Centromere', 'Reversed Segment', 'Moved Segment', 'ecDNA'];
  const g = svgElement(svg, 'g', {});
  svgElement(g, 'rect', { x, y, width: D_LEGEND_WIDTH, height: 34 + row * items.length, fill: '#ffffff', stroke: '#999', rx: 4 });
  svgElement(g, 'text', { x: x + 12, y: y + 22, 'font-weight': 'bold' }, 'Symbols:');
  items.forEach((name, i) => {
    const cy = y + 34 + row * i + row / 2, w = 16, cx = symbolX + w / 2;
    const bar = (top, height) => svgElement(g, 'rect', { x: symbolX, y: top, width: w, height, fill: color, stroke: '#1a1a1a', 'stroke-width': 1 });
    if (name === 'Chromosome') {
      svgElement(g, 'path', { d: piecePathSized(symbolX, cy - 18, cy + 18, w), fill: color, stroke: '#1a1a1a', 'stroke-width': 1.2 });
    } else if (name === 'Centromere') {
      svgElement(g, 'ellipse', { cx, cy, rx: w / 2 + 2, ry: 5, fill: '#999999' });
    } else if (name === 'Reversed Segment') {
      bar(cy - 12, 24);
      svgElement(g, 'path', { d: `M${symbolX},${cy + 12} L${cx},${cy - 12} L${symbolX + w},${cy + 12}`, fill: 'none', stroke: '#000', 'stroke-width': 1.5 });
    } else if (name === 'Moved Segment') {
      bar(cy - 15, 30);
      const tipX = symbolX - 1, fromY = cy + 9, toY = cy - 9, headBase = tipX - 6;
      svgElement(g, 'circle', { cx: tipX, cy: fromY, r: 1.8, fill: '#bf382b' });
      svgElement(g, 'path', { d: `M${tipX},${fromY} C${symbolX - 14},${fromY} ${symbolX - 14},${toY} ${headBase + 1},${toY}`, fill: 'none', stroke: '#bf382b', 'stroke-width': 1.5 });
      svgElement(g, 'path', { d: `M${tipX},${toY} L${headBase},${toY - 4} L${headBase + 1.5},${toY} L${headBase},${toY + 4} Z`, fill: '#bf382b' });
    } else if (name === 'ecDNA') {
      svgElement(g, 'circle', { cx, cy, r: 13, fill: color, stroke: '#1a1a1a', 'stroke-width': 1.5 });
      svgElement(g, 'circle', { cx, cy, r: 6.5, fill: '#ffffff', stroke: '#1a1a1a', 'stroke-width': 1.5 });
    } else if (name === 'Double-Strand Break') {
      bar(cy - 15, 30);
      svgElement(g, 'line', { x1: symbolX - 6, y1: cy, x2: symbolX + w + 6, y2: cy, stroke: '#000', 'stroke-width': 1.5 });
    } else if (name === 'Single-Strand Break') {
      bar(cy - 15, 30);
      svgElement(g, 'line', { x1: cx - 5, y1: cy, x2: cx + 5, y2: cy, stroke: '#fff', 'stroke-width': 2 });
    }
    svgElement(g, 'text', { x: textX, y: cy + 4 }, name);
  });
  return y + 34 + row * items.length;
}

// Counts of the open homolog, as its chromosome label's popup gives them: the popup section of that homolog (or the only section of a
// chromosome without homologs). Returns null for slots no label describes, e.g. a piece made of several chromosomes.
function slotCounts(slot) {
  const match = /^(.*), homolog (\d+)$/.exec(slot.title);
  const spot = karyograms[current].hotspots.find((h) => h.title === (match ? match[1] : slot.title));
  if (!spot) return null;
  if (!match) return spot.sections.length === 1 ? spot.sections[0].counts : null;
  const section = spot.sections.find((sec) => sec.heading === 'Homolog ' + match[2]);
  return section ? section.counts : null;
}

// Draws the open homolog's mutation (SDR) or damage (SDD) counts in a box below the symbol legend, one row per type with its count on the
// right, or the "none detected" text; SDR adds the note that events involving several chromosomes count on each. Alters the SVG only.
function drawDetailSummary(svg, x, y) {
  const counts = slotCounts(dSlot);
  if (!counts) return;
  const sdd = dSlot.pieces.some((p) => p.damages.length);
  const note = karyograms[current].popupNote ? ['Events involving several', 'chromosomes count on each.'] : [];
  const rows = counts.length ? counts.length : 1, lastRowY = y + 46 + 20 * (rows - 1), noteY = lastRowY + 26;
  const height = (note.length ? noteY + 14 * (note.length - 1) : lastRowY) + 14 - y;
  const g = svgElement(svg, 'g', {});
  svgElement(g, 'rect', { x, y, width: D_LEGEND_WIDTH, height, fill: '#ffffff', stroke: '#999', rx: 4 });
  svgElement(g, 'text', { x: x + 12, y: y + 22, 'font-weight': 'bold' }, sdd ? 'Damage:' : 'Mutations:');
  if (!counts.length) svgElement(g, 'text', { x: x + 12, y: y + 46, 'font-style': 'italic', fill: '#666' }, karyograms[current].popupEmptyText);
  counts.forEach(([type, count], i) => {
    svgElement(g, 'text', { x: x + 12, y: y + 46 + 20 * i }, type);
    svgElement(g, 'text', { x: x + D_LEGEND_WIDTH - 12, y: y + 46 + 20 * i, 'text-anchor': 'end', class: 'tick' }, String(count));
  });
  note.forEach((line, i) => svgElement(g, 'text', { x: x + 12, y: noteY + 14 * i, class: 'tick' }, line));
}

// Outline of a small rounded chromosome w px wide between y top and bottom, with a centromere pinch halfway, for the legend.
function piecePathSized(x, top, bottom, w) {
  const r = w / 2, cy = (top + bottom) / 2, depth = 2, half = 3;
  return `M${x},${top + r} A${r},${r} 0 0 1 ${x + w},${top + r} L${x + w},${cy - half} L${x + w - depth},${cy} L${x + w},${cy + half}`
    + ` L${x + w},${bottom - r} A${r},${r} 0 0 1 ${x},${bottom - r} L${x},${cy + half} L${x + depth},${cy} L${x},${cy - half} Z`;
}

// Redraws the open slot: ruler on the left, then each piece in its own column, captions on top. Alters the SVG and the scale readout.
function detailRender() {
  if (!dSlot) return;
  detailClamp();
  const h = detailView.clientHeight;
  const y = (mbp) => (mbp - dTop) / dMbpPerPx;
  dClipId = 0;
  detailSvg.replaceChildren();

  // Columns: room on the left of a piece for its arrows, on the right for indel and inversion labels.
  const legendX = 16;                                               // Symbol legend on the far left, then the ruler
  let rulerX = legendX + D_LEGEND_WIDTH + 100;
  let x = rulerX + 30;
  // Sister chromatids sit side by side with the karyogram's gap between them (scaled like the width) and no captions, as one X-shaped chromosome.
  const sisterGap = 1.5 * D_K;
  const columns = dSlot.pieces.filter((piece) => !piece.circular).map((piece, k) => {
    if (dSlot.sisters) {
      const column = { piece, x, sister: { side: k === 0 ? 'left' : 'right', gap: sisterGap } };
      x = column.x + D_WIDTH + (k === 0 ? sisterGap : 24);
      return column;
    }
    const left = piece.arrows.length ? 40 + 6 * piece.arrows.length : 0;
    const captionWidth = 7 * (piece.caption ? piece.caption.length : 0);   // Rough width of the caption
    const labelRoom = piece.markers.length ? 200 : piece.chevrons.length ? 150 : 24;   // Indel and "Inverted ..." labels go right of the piece
    const right = Math.max(labelRoom, captionWidth - D_WIDTH + 16);
    const column = { piece, x: x + left };
    x += left + D_WIDTH + right;
    return column;
  });
  // Circles (ecDNA, rings) are stacked top to bottom after the linear pieces, each with its label below it, fixed while zooming;
  // a circle that no longer fits the height starts a new stack to the right.
  let stackY = 36, stackWidth = 0;
  for (const piece of dSlot.pieces.filter((p) => p.circular)) {
    const diameter = circleDiameter(piece), caption = circleCaption(piece);
    const height = diameter + 8 + 15 * caption.length + 18;
    if (stackWidth && stackY + height > h) { x += stackWidth + 24; stackY = 36; stackWidth = 0; }
    columns.push({ piece, x, y: stackY, diameter, caption });
    stackY += height;
    stackWidth = Math.max(stackWidth, diameter, 7 * Math.max(...caption.map((line) => line.length)));
  }
  if (stackWidth) x += stackWidth + 24;
  // The panel is only as wide as its pieces need (at least wide enough for its toolbar, at most the window), so one homolog is not
  // left beside empty space. Wider content scrolls sideways.
  const panel = document.getElementById('detailPanel');
  panel.style.width = Math.min(window.innerWidth - 40, Math.max(760, x + 24)) + 'px';
  const width = Math.max(detailView.clientWidth, x);
  const shift = Math.max(0, (width - x) / 2);                       // Centre the ruler and pieces when the panel is wider than they are
  for (const column of columns) column.x += shift;
  detailSvg.setAttribute('width', width);
  detailSvg.setAttribute('height', h);
  detailSvg.setAttribute('viewBox', `0 0 ${width} ${h}`);

  for (const column of columns) {
    if (column.piece.circular) drawDetailCircle(detailSvg, column.piece, column);
    else drawDetailPiece(detailSvg, column.piece, column.x, h, column.sister);
  }
  // Sisters share one centromere ellipse across the waist where they touch (both outer pinches left out), then get their damage sites on top.
  const sisters = columns.filter((c) => c.sister);
  if (sisters.length === 2) {
    const [left] = sisters, segs = left.piece.segs.filter((s) => s[6] >= 0);
    drawDetailCentromeres(detailSvg, segs.map((s) => ({ top: y(s[6]), bottom: y(s[7]) })), left.x + D_WIDTH + sisterGap / 2, D_WIDTH - 7 * D_K + sisterGap / 2 + 1, h);
    for (const c of sisters) drawDetailDamages(detailSvg, c.piece, c.x, h, c.sister);
  }
  dColumns = columns;

  // Ruler: position along the homolog as drawn, continuous from its top end, with round ticks that follow the zoom.
  // A slot of circles only (ecDNA joining several chromosomes) has nothing on the ruler, so it gets none.
  const len = detailLength(), step = tickStep();
  rulerX += shift;
  if (dSlot.pieces.every((p) => p.circular)) rulerX = -1000;
  svgElement(detailSvg, 'line', { x1: rulerX, y1: Math.max(0, y(0)), x2: rulerX, y2: Math.min(h, y(len)), stroke: '#888' });
  const last = Math.min(len, dTop + h * dMbpPerPx);
  for (let i = Math.max(0, Math.ceil(dTop / step)); i * step <= last + 1e-12; i++) {
    const t = i * step, ty = y(t);
    svgElement(detailSvg, 'line', { x1: rulerX - 6, y1: ty, x2: rulerX, y2: ty, stroke: '#555' });
    svgElement(detailSvg, 'text', { x: rulerX - 9, y: ty + 4, 'text-anchor': 'end', class: 'tick' }, formatPosition(t, step));
  }

  // Captions (e.g. "Excised piece") stay at the top of the window; the chromosome itself has none, the ruler gives its length.
  svgElement(detailSvg, 'rect', { x: 0, y: 0, width, height: 26, fill: '#ffffff' });
  for (const column of columns) {
    if (column.piece.circular) continue;                            // Circles carry their label below them
    if (column.piece.caption && !column.sister) svgElement(detailSvg, 'text', { x: column.x, y: 17, class: 'caption' }, column.piece.caption);
  }
  drawDetailSummary(detailSvg, legendX + shift, drawDetailLegend(detailSvg, legendX + shift, 36) + 14);

  const pxPerBase = 1 / (dMbpPerPx * 1e6);
  detailScale.textContent = pxPerBase >= 1 ? `1 base = ${pxPerBase.toFixed(1)} px` : `1 px = ${formatLength(dMbpPerPx)}`;
}

// SDD damage sites within a few pixels of a screen point in the zoomed view, nearest first.
function damageSitesAt(clientX, clientY) {
  const r = detailSvg.getBoundingClientRect(), px = clientX - r.left, py = clientY - r.top;
  const linear = dColumns.filter((c) => !c.piece.circular);
  const column = linear.find((c) => px >= c.x && px <= c.x + D_WIDTH) || linear.find((c) => px >= c.x - 4 * D_K - 4 && px <= c.x + D_WIDTH + 4 * D_K + 4);
  if (!column) return [];
  const sites = column.piece.damages
    .map(([mbp, dsb, ssbs]) => ({ mbp, dsb, ssbs, distance: Math.abs((mbp - dTop) / dMbpPerPx - py) }))
    .filter((site) => site.distance <= 4)
    .sort((a, b) => a.distance - b.distance || a.mbp - b.mbp);
  // The sites the tooltip lists, to outline while they are pointed at; the box reaches past the piece as far as a DSB line does.
  if (sites.length) {
    const listed = sites.slice(0, 6).map((site) => site.mbp);
    sites.box = { x: column.x, start: Math.min(...listed), end: Math.max(...listed), padX: 4 * D_K + 4 };
  }
  return sites;
}

// Shows what the damage sites under the pointer are: per site its position, "Double Strand Break" when there is one and the number of SSBs.
// Sites closer together than the screen can separate are listed together. Alters the tooltip only.
function showDamageTip(clientX, clientY) {
  const sites = damageSitesAt(clientX, clientY);
  if (!sites.length) { hideDamageTip(); return; }
  detailTip.replaceChildren();
  if (sites.length > 1) addElement(detailTip, 'div', `${sites.length} damage sites here, zoom in to separate them`, 'head');
  for (const site of sites.slice(0, 6)) {
    const group = addElement(detailTip, 'div', undefined, 'tipsite');
    addElement(group, 'div', Math.round(site.mbp * 1e6).toLocaleString('en-US') + ' bp', 'site');
    if (site.dsb) addElement(group, 'div', 'Double Strand Break');
    if (site.ssbs > 0) addElement(group, 'div', 'Single Strand Breaks = ' + site.ssbs);
  }
  if (sites.length > 6) addElement(detailTip, 'div', `and ${sites.length - 6} more`, 'more');
  detailTip.hidden = false;
  const box = detailTip.getBoundingClientRect();
  detailTip.style.left = Math.max(8, Math.min(clientX + 14, window.innerWidth - box.width - 8)) + 'px';
  detailTip.style.top = Math.max(8, Math.min(clientY + 14, window.innerHeight - box.height - 8)) + 'px';
}

// Hides the tooltip and the outline of the stretch it described.
function hideDamageTip() { detailTip.hidden = true; showStretchBox(null); }

// Outlines the stretch or damage sites the tooltip describes with a light blue box a little wider than the piece (padX px each side) and at
// least 6 px tall, so a white gap of a few bases or a single damage site is still boxed; removes the outline when box is null. Alters the SVG only.
function showStretchBox(box) {
  const old = document.getElementById('detailStretchBox');
  if (old) old.remove();
  if (!box) return;
  const h = detailView.clientHeight, y = (mbp) => (mbp - dTop) / dMbpPerPx;
  let top = Math.max(-4, y(box.start)), bottom = Math.min(h + 4, y(box.end));
  if (bottom - top < 6) { const c = (top + bottom) / 2; top = c - 3; bottom = c + 3; }
  const padX = box.padX || 4;
  svgElement(detailSvg, 'rect', { id: 'detailStretchBox', x: box.x - padX, y: top - 2, width: D_WIDTH + 2 * padX, height: bottom - top + 4, rx: 3,
    fill: 'none', stroke: '#4fb3ff', 'stroke-width': 2.5, 'pointer-events': 'none' });
}

// Formats an original span in Mbp with up to six decimals (single bases), e.g. "40.000001 – 60 Mbp".
function formatSpan(startMbp, endMbp) {
  const mbp = (value) => String(Number(value.toFixed(6)));
  return `${mbp(startMbp)} – ${mbp(endMbp)} Mbp`;
}

// What the pointer is over in an SDR chromosome or excised piece: the chromosome its bases came from and their original span, for every
// stretch (a white gap gives the chromosome and span of the deleted bases). No mutation names: those are the detectors' reading, the tooltip
// gives only what the SDR record says. Circles (ecDNA, rings) get none, their label already gives their span. Returns the tooltip lines, the
// first one bold, or null.
function segmentTipAt(clientX, clientY) {
  const r = detailSvg.getBoundingClientRect(), px = clientX - r.left, py = clientY - r.top;
  const column = dColumns.find((c) => !c.piece.circular && px >= c.x && px <= c.x + D_WIDTH);
  if (!column || column.piece.damages.length) return null;
  const piece = column.piece, segs = piece.segs, at = dTop + py * dMbpPerPx, slack = 1.5 * dMbpPerPx;
  if (at < 0 || at > piece.len) return null;
  // A white gap first, with a little slack so even a gap of a few bases can be aimed at; otherwise the stretch under the pointer.
  let index = segs.findIndex((s) => s[3] && at >= s[0] - slack && at <= s[1] + slack);
  if (index < 0) index = segs.findIndex((s) => at >= s[0] && at <= s[1]);
  if (index < 0 || segs[index][8] < 0) return null;
  let name = segs[index][10];
  for (let k = index - 1; !name && k >= 0; k--) if (!segs[k][3]) name = segs[k][10];          // A gap belongs to the chromosome around it
  for (let k = index + 1; !name && k < segs.length; k++) if (!segs[k][3]) name = segs[k][10];
  const lines = [name || 'Unknown chromosome', formatSpan(segs[index][8], segs[index][9])];
  lines.box = { x: column.x, start: segs[index][0], end: segs[index][1] };   // The stretch to outline while it is pointed at
  return lines;
}

// Shows the tooltip for whatever is under the pointer: SDD damage sites first, otherwise the span of a deleted segment or excised piece (circles carry their span in their label).
function showDetailTip(clientX, clientY) {
  const sites = damageSitesAt(clientX, clientY);
  if (sites.length) { showDamageTip(clientX, clientY); showStretchBox(sites.box); return; }
  const lines = segmentTipAt(clientX, clientY);
  if (!lines) { hideDamageTip(); return; }
  showStretchBox(lines.box);
  detailTip.replaceChildren();
  const group = addElement(detailTip, 'div', undefined, 'tipsite');
  lines.forEach((line, i) => addElement(group, 'div', line, i === 0 ? 'site' : undefined));
  detailTip.hidden = false;
  const box = detailTip.getBoundingClientRect();
  detailTip.style.left = Math.max(8, Math.min(clientX + 14, window.innerWidth - box.width - 8)) + 'px';
  detailTip.style.top = Math.max(8, Math.min(clientY + 14, window.innerHeight - box.height - 8)) + 'px';
}

// Keyboard in the zoomed view: Esc closes, + and - zoom, 0 shows the whole homolog.
function detailKey(e) {
  if (e.key === 'Escape') closeDetail();
  else if (e.ctrlKey || e.metaKey || e.altKey) return;
  else if (e.key === '+' || e.key === '=') detailZoomAt(2);
  else if (e.key === '-' || e.key === '_') detailZoomAt(0.5);
  else if (e.key === '0') detailWhole();
}

detailView.addEventListener('wheel', (e) => {
  e.preventDefault();
  hideDamageTip();
  const delta = e.deltaMode === 1 ? e.deltaY * 33 : e.deltaY;
  detailZoomAt(Math.exp(-delta * (e.ctrlKey ? 0.01 : 0.0025)), e.clientY - detailView.getBoundingClientRect().top);
}, { passive: false });

// Drag (mouse or one finger) moves along the homolog and sideways when the pieces are wider than the window; two fingers pinch to zoom.
const dPointers = new Map();
let dPinch = null;
detailView.addEventListener('pointerdown', (e) => {
  if (e.button !== 0) return;
  detailView.setPointerCapture(e.pointerId);
  dPointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
  dDrag = { x: e.clientX, y: e.clientY, moved: false };
  detailView.classList.add('dragging');
});
detailView.addEventListener('pointermove', (e) => {
  const previous = dPointers.get(e.pointerId);
  if (!previous) { if (e.pointerType === 'mouse') showDetailTip(e.clientX, e.clientY); return; }   // Hovering (not dragging) shows a tooltip.
  hideDamageTip();
  dPointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
  if (dDrag && Math.abs(e.clientX - dDrag.x) + Math.abs(e.clientY - dDrag.y) > 4) dDrag.moved = true;
  if (dPointers.size === 2) {
    const [a, b] = [...dPointers.values()], distance = Math.hypot(a.x - b.x, a.y - b.y);
    if (dPinch) detailZoomAt(distance / dPinch, (a.y + b.y) / 2 - detailView.getBoundingClientRect().top);
    dPinch = distance || 1;
    return;
  }
  detailView.scrollLeft -= e.clientX - previous.x;
  dTop -= (e.clientY - previous.y) * dMbpPerPx;
  detailRender();
});
const endDetailDrag = (e) => {
  if (dDrag && !dDrag.moved && dPointers.size === 1) {               // A click on an indel pointer zooms to it.
    const marker = document.elementFromPoint(e.clientX, e.clientY)?.closest('.marker');
    if (marker) detailJumpTo(Number(marker.dataset.mbp));
    else showDetailTip(e.clientX, e.clientY);                        // A tap shows the tooltip, for touch screens.
  }
  dPointers.delete(e.pointerId);
  dPinch = null;
  if (!dPointers.size) { dDrag = null; detailView.classList.remove('dragging'); }
};
detailView.addEventListener('pointerup', endDetailDrag);
detailView.addEventListener('pointercancel', endDetailDrag);
detailView.addEventListener('pointerleave', (e) => { if (e.pointerType === 'mouse') hideDamageTip(); });
detailView.addEventListener('dblclick', (e) => detailZoomAt(2, e.clientY - detailView.getBoundingClientRect().top));
detail.addEventListener('pointerdown', (e) => { if (e.target === detail) closeDetail(); });   // Clicking outside the panel closes it.
document.getElementById('detailIn').onclick = () => detailZoomAt(2);
document.getElementById('detailOut').onclick = () => detailZoomAt(0.5);
document.getElementById('detailWhole').onclick = detailWhole;
document.getElementById('detailClose').onclick = closeDetail;
window.addEventListener('resize', () => { if (dSlot) detailRender(); });

// ------------------------------------------------------------
// Downloads
// ------------------------------------------------------------
const formatSelect = document.getElementById('format');
for (const option of formatSelect.options) {                       // Only offer formats that were embedded.
  option.disabled = !karyograms.some((k) => k[option.value]);
}
if (formatSelect.selectedOptions[0].disabled) formatSelect.value = 'svg';

// "cell0_karyogram_human.png" -> "cell0_karyogram_human.svg" for SVG downloads.
function downloadName(k, format) {
  return k.filename.replace(/\.(png|svg)$/i, '') + '.' + format;
}

// Decode an embedded base64 data: URL back into the original file bytes.
function dataUrlBytes(dataUrl) {
  const binary = atob(dataUrl.slice(dataUrl.indexOf(',') + 1));
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  return bytes;
}

// Hand a file to the browser's normal download (it asks where to save, or uses the Downloads folder).
function saveFile(bytes, name, type) {
  const url = URL.createObjectURL(new Blob([bytes], { type }));
  const link = document.createElement('a');
  link.href = url;
  link.download = name;
  document.body.appendChild(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

const MIME = { png: 'image/png', svg: 'image/svg+xml' };

document.getElementById('download').onclick = () => {
  const k = karyograms[current], format = formatSelect.value;
  if (!k[format]) { alert('This karyogram has no ' + format.toUpperCase() + ' version.'); return; }
  saveFile(dataUrlBytes(k[format]), downloadName(k, format), MIME[format]);
};

// Build a .zip file in the browser (stored, not compressed: PNGs are already compressed and this keeps the code short).
// Why a zip: browsers block or prompt for many separate downloads at once, a single file avoids that.
const CRC_TABLE = (() => {
  const table = new Uint32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1;
    table[n] = c >>> 0;
  }
  return table;
})();

function crc32(bytes) {
  let c = 0xFFFFFFFF;
  for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xFF] ^ (c >>> 8);
  return (c ^ 0xFFFFFFFF) >>> 0;
}

function makeZip(files) {
  const encoder = new TextEncoder();
  const chunks = [], central = [];
  let offset = 0;
  const now = new Date();
  const dosTime = (now.getHours() << 11) | (now.getMinutes() << 5) | (now.getSeconds() >> 1);
  const dosDate = ((now.getFullYear() - 1980) << 9) | ((now.getMonth() + 1) << 5) | now.getDate();

  for (const file of files) {
    const name = encoder.encode(file.name), crc = crc32(file.bytes), size = file.bytes.length;

    const local = new DataView(new ArrayBuffer(30));
    local.setUint32(0, 0x04034B50, true);   // Local file header signature
    local.setUint16(4, 20, true);           // Version needed to extract
    local.setUint16(6, 0x0800, true);       // Flags: file name is UTF-8
    local.setUint16(8, 0, true);            // Method: stored
    local.setUint16(10, dosTime, true);
    local.setUint16(12, dosDate, true);
    local.setUint32(14, crc, true);
    local.setUint32(18, size, true);        // Compressed size
    local.setUint32(22, size, true);        // Uncompressed size
    local.setUint16(26, name.length, true);
    local.setUint16(28, 0, true);           // Extra field length
    chunks.push(new Uint8Array(local.buffer), name, file.bytes);

    const entry = new DataView(new ArrayBuffer(46));
    entry.setUint32(0, 0x02014B50, true);   // Central directory signature
    entry.setUint16(4, 20, true);           // Version made by
    entry.setUint16(6, 20, true);           // Version needed
    entry.setUint16(8, 0x0800, true);
    entry.setUint16(10, 0, true);
    entry.setUint16(12, dosTime, true);
    entry.setUint16(14, dosDate, true);
    entry.setUint32(16, crc, true);
    entry.setUint32(20, size, true);
    entry.setUint32(24, size, true);
    entry.setUint16(28, name.length, true);
    entry.setUint32(42, offset, true);      // Where this file's local header starts
    central.push(new Uint8Array(entry.buffer), name);

    offset += 30 + name.length + size;
  }

  const centralSize = central.reduce((total, part) => total + part.length, 0);
  const end = new DataView(new ArrayBuffer(22));
  end.setUint32(0, 0x06054B50, true);       // End of central directory signature
  end.setUint16(8, files.length, true);
  end.setUint16(10, files.length, true);
  end.setUint32(12, centralSize, true);
  end.setUint32(16, offset, true);
  return new Blob([...chunks, ...central, new Uint8Array(end.buffer)], { type: 'application/zip' });
}

document.getElementById('downloadAll').onclick = () => {
  const format = formatSelect.value;
  const files = karyograms.filter((k) => k[format])
    .map((k) => ({ name: downloadName(k, format), bytes: dataUrlBytes(k[format]) }));
  const zipName = TITLE.replace(/\.[^.]*$/, '') + '_karyograms_' + format + '.zip';
  saveFile(makeZip(files), zipName, 'application/zip');
};

window.addEventListener('resize', () => { if (fitMode) fit(); else setScale(scale); });

show(0);
</script>
</body>
</html>
)HTML";

    html.close();

    if (!html)
    {
        std::cerr << "ERROR: Failed to write web viewer file " << htmlPath << "\n";
        return false;
    }

    return true;
}


bool writeZipArchive(
    const std::vector<std::pair<std::string, std::vector<unsigned char>>>& entries,
    const std::filesystem::path& zipPath)
{
    // CRC-32 checksum the zip format stores for every file, so unzip tools can check the contents.
    const auto crc32 = [](const std::vector<unsigned char>& bytes) -> std::uint32_t
    {
        std::uint32_t crc = 0xFFFFFFFFu;

        for (unsigned char byte : bytes)
        {
            crc ^= byte;

            for (int bit = 0; bit < 8; ++bit)
            {
                crc = (crc >> 1) ^ (0xEDB88320u & (0u - (crc & 1u)));
            }
        }

        return crc ^ 0xFFFFFFFFu;
    };

    // Zip numbers are little-endian, written byte by byte so the file is the same on every machine.
    const auto put16 = [](std::vector<unsigned char>& out, std::uint32_t value)
    {
        out.push_back(static_cast<unsigned char>(value & 0xFF));
        out.push_back(static_cast<unsigned char>((value >> 8) & 0xFF));
    };
    const auto put32 = [&put16](std::vector<unsigned char>& out, std::uint32_t value)
    {
        put16(out, value & 0xFFFF);
        put16(out, value >> 16);
    };

    // File date and time in the zip's MS-DOS format, so the extracted files show when they were made.
    const std::time_t now = std::time(nullptr);
    const std::tm* local = std::localtime(&now);
    const std::uint32_t dosTime = local ? static_cast<std::uint32_t>((local->tm_hour << 11) | (local->tm_min << 5) | (local->tm_sec / 2)) : 0;
    const std::uint32_t dosDate = local ? static_cast<std::uint32_t>(((local->tm_year - 80) << 9) | ((local->tm_mon + 1) << 5) | local->tm_mday) : 0x21;

    std::vector<unsigned char> archive;
    std::vector<unsigned char> centralDirectory;

    for (const auto& [name, data] : entries)
    {
        // Safeguard: this simple writer has no ZIP64 support, so each file and the whole archive must stay under 4 GB.
        if (data.size() > 0xFFFFFFFFu || archive.size() > 0xFFFFFFFFu || name.size() > 0xFFFFu)
        {
            std::cerr << "ERROR: " << name << " is too large to add to the zip file " << zipPath << "\n";
            return false;
        }

        const std::uint32_t crc = crc32(data);
        const std::uint32_t size = static_cast<std::uint32_t>(data.size());
        const std::uint32_t offset = static_cast<std::uint32_t>(archive.size());

        put32(archive, 0x04034B50u);                                            // Local file header
        put16(archive, 20);                                                     // Version needed to extract (2.0)
        put16(archive, 0);                                                      // Flags
        put16(archive, 0);                                                      // Stored, no compression
        put16(archive, dosTime);
        put16(archive, dosDate);
        put32(archive, crc);
        put32(archive, size);                                                   // Compressed size (same, stored)
        put32(archive, size);
        put16(archive, static_cast<std::uint32_t>(name.size()));
        put16(archive, 0);                                                      // No extra field
        archive.insert(archive.end(), name.begin(), name.end());
        archive.insert(archive.end(), data.begin(), data.end());

        put32(centralDirectory, 0x02014B50u);                                   // Central directory entry for the same file
        put16(centralDirectory, 20);                                            // Made by version 2.0
        put16(centralDirectory, 20);
        put16(centralDirectory, 0);
        put16(centralDirectory, 0);
        put16(centralDirectory, dosTime);
        put16(centralDirectory, dosDate);
        put32(centralDirectory, crc);
        put32(centralDirectory, size);
        put32(centralDirectory, size);
        put16(centralDirectory, static_cast<std::uint32_t>(name.size()));
        put16(centralDirectory, 0);                                             // No extra field
        put16(centralDirectory, 0);                                             // No comment
        put16(centralDirectory, 0);                                             // Disk number
        put16(centralDirectory, 0);                                             // Internal attributes
        put32(centralDirectory, 0);                                             // External attributes
        put32(centralDirectory, offset);
        centralDirectory.insert(centralDirectory.end(), name.begin(), name.end());
    }

    if (entries.size() > 0xFFFFu)
    {
        std::cerr << "ERROR: Too many files for the zip file " << zipPath << "\n";
        return false;
    }

    const std::uint32_t directoryOffset = static_cast<std::uint32_t>(archive.size());
    archive.insert(archive.end(), centralDirectory.begin(), centralDirectory.end());

    put32(archive, 0x06054B50u);                                                // End of central directory record
    put16(archive, 0);
    put16(archive, 0);
    put16(archive, static_cast<std::uint32_t>(entries.size()));
    put16(archive, static_cast<std::uint32_t>(entries.size()));
    put32(archive, static_cast<std::uint32_t>(centralDirectory.size()));
    put32(archive, directoryOffset);
    put16(archive, 0);                                                          // No archive comment

    std::ofstream zipFile(zipPath, std::ios::binary);
    zipFile.write(reinterpret_cast<const char*>(archive.data()), static_cast<std::streamsize>(archive.size()));

    if (!zipFile)
    {
        std::cerr << "ERROR: Failed to write zip file " << zipPath << "\n";
        return false;
    }

    return true;
}


bool canOpenBrowser()
{
#if defined(_WIN32)
    return true;                                                                // Native Windows always has a desktop and a default browser.
#elif defined(__APPLE__)
    return !isRemoteSSHsession();                                               // A Mac reached over SSH would open the browser on its own (unseen) screen.
#else
    if (isRunningInWSL())
    {
        return wslCanLaunchWindowsPrograms();                                   // WSL opens the Windows browser, which needs Windows interop.
    }

    // Plain Linux: a browser can only appear if there is a graphical display. Desktops set DISPLAY (X11) or
    // WAYLAND_DISPLAY (Wayland); SSH with X forwarding (ssh -X) also sets DISPLAY, so that case still works.
    // Why: without a display, xdg-open either fails or, worse, starts a text-mode browser (such as lynx or w3m)
    // inside the terminal, which takes over the terminal and blocks the program until the user quits it.
    return std::getenv("DISPLAY") != nullptr || std::getenv("WAYLAND_DISPLAY") != nullptr;
#endif
}


bool openInBrowser(const std::filesystem::path& htmlPath)
{
    if (!canOpenBrowser())                                                      // Safeguard for headless machines, see canOpenBrowser().
    {
        return false;
    }

    const std::string path = std::filesystem::absolute(htmlPath).string();

#if defined(_WIN32)
    // Windows file names cannot contain double quotes, so a quote here means something is wrong with the path.
    // Refusing it keeps the quoting of the "start" command below intact.
    if (path.find('"') != std::string::npos)
    {
        return false;
    }
    const std::string command = "start \"\" \"" + path + "\"";                    // The empty "" is the window title, start would otherwise treat the quoted path as the title.
#elif defined(__APPLE__)
    const std::string command = "open " + shellQuote(path);
#else
    if (isRunningInWSL())
    {
        // wslpath -w turns /tmp/x.html into \\wsl.localhost\<distro>\tmp\x.html, a path Windows programs can read,
        // and explorer.exe opens it with the user's default Windows browser.
        const std::string wslCommand =
            "explorer.exe \"$(wslpath -w " + shellQuote(path) + ")\" > /dev/null 2>&1";

        // explorer.exe returns 1 even when it succeeds, so its exit code cannot tell us anything. It is stored in a
        // variable only because some compilers warn when the result of std::system() is ignored.
        const int explorerStatus = std::system(wslCommand.c_str());
        (void)explorerStatus;
        return true;
    }

    const std::string command = "xdg-open " + shellQuote(path) + " > /dev/null 2>&1";   // Output is hidden because xdg-open can print browser noise into the CLI output.
#endif

    return std::system(command.c_str()) == 0;
}

} // namespace karyoke
