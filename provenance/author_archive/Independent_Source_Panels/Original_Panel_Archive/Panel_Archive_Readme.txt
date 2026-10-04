Original Figure Panel Archive

165 original labeled panels are preserved: 68 panels from main Figures 1-6 and 97 panels from Supplementary Figures S1-S10. Shared titles, legends, and a source footnote are retained as separate supporting elements. No analysis value, point, line, axis, color scale, symbol, annotation, or result text has been redrawn.

The independent/ directory uses original figure and panel identifiers. The reassembled/ directory follows the revised narrative numbering:
Original 1 -> revised 1
Original 2 -> revised 2
Original 3 -> revised 3
Original 4 -> revised 5
Original 5 -> revised 6
Original 6 -> revised 7
Original S1-S10 -> revised S1-S10
The new binding figure (revised Figure 4) and new virtual-perturbation figure (Figure S11) are supplied separately.

Image preservation
TIFF exports use lossless LZW compression. Original main-figure images have 300 dpi metadata; supplementary images have 1500 dpi metadata. These values are preserved. Images have not been enlarged to imply greater information content. Main-figure reassembly trims pure-white border space and relocates complete original panels without resampling.

Original Figure 3 has interleaved E/F axis labels and G/H panel letters. Original Figure 4 has an E size legend extending into the F left-margin region. These six panels use explicitly recorded neighbouring-text exclusions so their original graphical content can be separated. Excluded pixels belong to another delivered panel. Pixel-by-pixel reconstruction of the original composite is checked from all independent panels and shared elements; the reconstruction report records this check for every original figure.

Manifest and reproducibility
panel_crop_manifest.csv records the source SHA-256, exact pixel coordinates, dimensions, original dpi, output pixel SHA-256, original rectangle pixel SHA-256, any neighbouring-text exclusion rectangles, panel description, and revised figure destination. A pixel SHA-256 is calculated from decoded RGB bytes, independent of TIFF encoding.

scripts/reassemble_delivered_panels.py reconstructs the six retained main figures directly from this archive without needing the original composite images. scripts/extract_original_panels.py reproduces source extraction and main-figure reassembly. scripts/verify_panel_reconstruction.py places the panels back at their original source coordinates and checks exact equality with each source figure. The scripts expect the original extracted submission archive under work/source/转化医学投稿20260922 and panel artifacts under work/panels, relative to the workspace root. The source archive itself is unchanged.

contact_sheets/ and previews/ are review-only reduced images. Use the independent or reassembled TIFF files for figure composition and submission.
