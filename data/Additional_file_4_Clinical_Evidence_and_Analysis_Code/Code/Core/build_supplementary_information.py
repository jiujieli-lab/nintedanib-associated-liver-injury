#!/usr/bin/env python3
"""Build the journal-facing Supplementary Information PDF and Markdown.

The PDF preserves vector figure content by placing each source PDF onto a
letter-sized page.  A separate legend page follows every figure so that
full-resolution plots are never reduced to make room for long captions.
"""

from __future__ import annotations

import re
import os
import tempfile
from functools import partial
from io import BytesIO
from pathlib import Path
from xml.sax.saxutils import escape

from pypdf import PdfReader, PdfWriter, Transformation
from pypdf._page import PageObject
from pypdf.generic import DecodedStreamObject, NameObject
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A3, landscape, letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import PageBreak, Paragraph, SimpleDocTemplate, Spacer


ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "Supplementary_Figures"
OUT = ROOT
OUT.mkdir(parents=True, exist_ok=True)

TITLE = "Role-separated public-data triangulation prioritizes target–phenotype hypotheses for nintedanib-associated liver injury"
BLUE = HexColor("#2E74B5")
INK = HexColor("#0B2545")
GRAY = HexColor("#5B6573")


def figure_specs() -> list[tuple[str, Path, Path]]:
    specs: list[tuple[str, Path, Path]] = [
        ("S1", FIG / "Supplementary_Figure_S1_FAERS.pdf", FIG / "Supplementary_Figure_S1_FAERS_caption.md"),
        ("S2", FIG / "Supplementary_Figure_S2_Protein_Distributions.pdf", FIG / "Supplementary_Figure_S2_Protein_Distributions_caption.md"),
        ("S3", FIG / "Supplementary_Figure_S3_Proteomics_Robustness.pdf", FIG / "Supplementary_Figure_S3_Proteomics_Robustness_caption.md"),
        ("S4", FIG / "Supplementary_Figure_S4_Nintedanib_Perturbation_Robustness.pdf", FIG / "Supplementary_Figure_S4_Nintedanib_Perturbation_Robustness_caption.md"),
        ("S5", FIG / "Supplementary_Figure_S5_Whole_Blood_Sensitivity.pdf", FIG / "Supplementary_Figure_S5_Whole_Blood_Sensitivity_caption.md"),
        ("S6", FIG / "Supplementary_Figure_S6_GSE151374.pdf", FIG / "Supplementary_Figure_S6_GSE151374_caption.md"),
        ("S7", FIG / "Supplementary_Figure_S7_PXD052594.pdf", FIG / "Supplementary_Figure_S7_PXD052594_caption.md"),
    ]
    for label, pdf, caption in specs:
        if not pdf.exists() or not caption.exists():
            raise FileNotFoundError(f"Missing {label}: {pdf} or {caption}")
    return specs


def clean_caption(label: str, raw: str) -> str:
    text = raw.strip()
    text = re.sub(r"^#\s*", "", text)
    text = re.sub(
        r"^Supplementary Figure S?\d+\s*[.|]\s*",
        f"Supplementary Figure {label}. ",
        text,
        count=1,
        flags=re.I,
    )
    if not text.lower().startswith(f"supplementary figure {label.lower()}"):
        text = f"Supplementary Figure {label}. {text}"
    return text


def caption_for_pdf(label: str, raw: str) -> str:
    """Convert normalized Markdown caption text to ReportLab-safe markup."""
    text = escape(clean_caption(label, raw), entities={"\"": "&quot;"})
    text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<i>\1</i>", text)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    return text


def register_font() -> str:
    """Register an embedded four-face sans family for every generated text run."""
    runtime_root = os.environ.get("CODEX_PRIMARY_RUNTIME_ROOT")
    runtime_fonts = (
        Path(runtime_root)
        / "dependencies/native/libreoffice-headless/libreoffice/share/fonts/truetype"
        if runtime_root
        else None
    )
    families = [
        {
            "regular": Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
            "bold": Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
            "italic": Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf"),
            "bold_italic": Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf"),
        },
        {
            "regular": Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf"),
            "bold": Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"),
            "italic": Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Italic.ttf"),
            "bold_italic": Path("/usr/share/fonts/truetype/liberation2/LiberationSans-BoldItalic.ttf"),
        },
    ]
    if runtime_fonts is not None:
        families.insert(
            0,
            {
                "regular": runtime_fonts / "DejaVuSans.ttf",
                "bold": runtime_fonts / "DejaVuSans-Bold.ttf",
                "italic": runtime_fonts / "DejaVuSans-Oblique.ttf",
                "bold_italic": runtime_fonts / "DejaVuSans-BoldOblique.ttf",
            },
        )
    for paths in families:
        if all(path.exists() for path in paths.values()):
            pdfmetrics.registerFont(TTFont("BodySans", str(paths["regular"])))
            pdfmetrics.registerFont(TTFont("BodySans-Bold", str(paths["bold"])))
            pdfmetrics.registerFont(TTFont("BodySans-Italic", str(paths["italic"])))
            pdfmetrics.registerFont(TTFont("BodySans-BoldItalic", str(paths["bold_italic"])))
            pdfmetrics.registerFontFamily(
                "BodySans",
                normal="BodySans",
                bold="BodySans-Bold",
                italic="BodySans-Italic",
                boldItalic="BodySans-BoldItalic",
            )
            return "BodySans"
    raise FileNotFoundError("No complete embeddable sans-serif font family is available")


def text_pdf(specs: list[tuple[str, Path, Path]], font_name: str, path: Path) -> None:
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "SuppTitle", parent=styles["Title"], fontName=font_name,
        fontSize=19, leading=24, textColor=INK, alignment=TA_CENTER,
        spaceAfter=16,
    )
    sub_style = ParagraphStyle(
        "SuppSub", parent=styles["Normal"], fontName=font_name,
        fontSize=10, leading=14, textColor=GRAY, alignment=TA_CENTER,
    )
    h1 = ParagraphStyle(
        "SuppH1", parent=styles["Heading1"], fontName=font_name,
        fontSize=15, leading=19, textColor=BLUE, alignment=TA_LEFT,
        spaceAfter=10,
    )
    body = ParagraphStyle(
        "SuppBody", parent=styles["BodyText"], fontName=font_name,
        fontSize=9.5, leading=13, textColor=INK, alignment=TA_LEFT,
        spaceAfter=7,
    )
    small = ParagraphStyle(
        "SuppSmall", parent=body, fontSize=8.5, leading=12, textColor=GRAY,
    )
    doc = SimpleDocTemplate(
        str(path), pagesize=letter,
        rightMargin=0.8 * inch, leftMargin=0.8 * inch,
        topMargin=0.72 * inch, bottomMargin=0.72 * inch,
        title=f"Supplementary Information — {TITLE}",
        author="Study authors (administrative fields pending author verification)",
    )
    story = [
        Spacer(1, 0.75 * inch),
        Paragraph("Supplementary Information", title_style),
        Paragraph(TITLE, h1),
        Spacer(1, 0.18 * inch),
        Paragraph("Cell Reports Medicine — Research Article", sub_style),
        Spacer(1, 0.4 * inch),
        Paragraph("[AUTHOR LIST]", sub_style),
        Paragraph("[AFFILIATIONS]", sub_style),
        PageBreak(),
        Paragraph("Contents and interpretation boundary", h1),
    ]
    for label, _, _ in specs:
        story.append(Paragraph(f"Supplementary Figure {label} and legend", body))
    story.extend([
        Spacer(1, 0.15 * inch),
        Paragraph(
            "The figures report public-data reanalysis. FAERS disproportionality does not estimate incidence or causality; "
            "human serum proteomics is multi-drug DILI rather than nintedanib-specific; virtual knockout is unsigned model-derived "
            "network displacement; pooled tissue slices and pooled single-cell libraries retain their source-level experimental units; "
            "and public lung datasets define model-dependent pulmonary-response constraints rather than establish hepatic toxicity.",
            small,
        ),
    ])
    for label, _, caption_path in specs:
        story.extend([
            PageBreak(),
            Paragraph(f"Supplementary Figure {label} legend", h1),
            Paragraph(caption_for_pdf(label, caption_path.read_text(encoding="utf-8")), body),
        ])
    doc.build(story, canvasmaker=partial(Canvas, initialFontName=font_name))


def overlay_page(label: str, page_number: int, font_name: str, page_size: tuple[float, float]) -> PageObject:
    buf = BytesIO()
    c = Canvas(buf, pagesize=page_size, initialFontName=font_name)
    c.setFont(font_name, 8)
    c.setFillColor(GRAY)
    c.drawString(36, page_size[1] - 26, f"Supplementary Figure {label}")
    c.drawRightString(page_size[0] - 36, 22, f"Page {page_number}")
    c.save()
    buf.seek(0)
    return PdfReader(buf).pages[0]


def footer_page(page_number: int, font_name: str) -> PageObject:
    buf = BytesIO()
    c = Canvas(buf, pagesize=letter, initialFontName=font_name)
    c.setFont(font_name, 8)
    c.setFillColor(GRAY)
    c.drawRightString(letter[0] - 36, 22, f"Page {page_number}")
    c.save()
    buf.seek(0)
    return PdfReader(buf).pages[0]


def placed_figure_page(label: str, source_pdf: Path, page_number: int, font_name: str) -> PageObject:
    src_reader = PdfReader(str(source_pdf))
    if len(src_reader.pages) != 1:
        raise ValueError(f"Expected one-page source figure: {source_pdf}")
    src = src_reader.pages[0]
    sw, sh = float(src.mediabox.width), float(src.mediabox.height)
    # Portrait figures S1–S5 fit Letter at near-native scale. The two dense
    # ten-panel extensions use A3 landscape so their native 7-pt labels remain
    # print-legible instead of being reduced below 4 pt on Letter.
    page_size = landscape(A3) if label in {"S6", "S7"} else letter
    page = PageObject.create_blank_page(width=page_size[0], height=page_size[1])
    available_w = page_size[0] - 56
    available_h = page_size[1] - 76
    scale = min(available_w / sw, available_h / sh)
    tx = (page_size[0] - sw * scale) / 2
    ty = 34 + (available_h - sh * scale) / 2
    page.merge_transformed_page(src, Transformation().scale(scale).translate(tx, ty))
    page.merge_page(overlay_page(label, page_number, font_name, page_size))
    return page


def remove_unused_reportlab_base_font(page: PageObject) -> None:
    """Remove ReportLab's empty Helvetica preamble and its unused resource."""
    contents = page.get_contents()
    if contents is None:
        return
    data = contents.get_data()
    cleaned = re.sub(
        rb"\s*BT\s*/F1\s+12\s+Tf\s+14(?:\.4)?\s+TL\s+ET",
        b"",
        data,
        count=1,
    )
    if cleaned != data:
        stream = DecodedStreamObject()
        stream.set_data(cleaned)
        page.replace_contents(stream)
    if re.search(rb"/F1(?![+A-Za-z0-9_-])", cleaned) is None:
        resources = page.get("/Resources")
        fonts = resources.get("/Font") if resources is not None else None
        if fonts is not None:
            fonts = fonts.get_object()
            if NameObject("/F1") in fonts:
                del fonts[NameObject("/F1")]


def build_pdf(specs: list[tuple[str, Path, Path]], text_path: Path, out_path: Path, font_name: str) -> None:
    text_reader = PdfReader(str(text_path))
    expected_text_pages = 2 + len(specs)
    if len(text_reader.pages) != expected_text_pages:
        raise ValueError(
            f"Expected exactly {expected_text_pages} text pages (cover, contents, and one legend per figure); "
            f"found {len(text_reader.pages)}. Adjust caption typography before assembly."
        )
    # Preserve cover and contents first, then interleave a figure page before each one-page legend.
    writer = PdfWriter()
    writer.add_page(text_reader.pages[0])
    writer.pages[-1].merge_page(footer_page(1, font_name))
    remove_unused_reportlab_base_font(writer.pages[-1])
    writer.add_page(text_reader.pages[1])
    writer.pages[-1].merge_page(footer_page(2, font_name))
    remove_unused_reportlab_base_font(writer.pages[-1])
    text_index = 2
    page_number = 3
    for label, pdf, _ in specs:
        writer.add_page(placed_figure_page(label, pdf, page_number, font_name))
        page_number += 1
        if text_index >= len(text_reader.pages):
            raise ValueError(f"Missing legend page for {label}")
        writer.add_page(text_reader.pages[text_index])
        writer.pages[-1].merge_page(footer_page(page_number, font_name))
        remove_unused_reportlab_base_font(writer.pages[-1])
        text_index += 1
        page_number += 1
    assert text_index == len(text_reader.pages)
    writer.add_metadata({
        "/Title": f"Supplementary Information - {TITLE}",
        "/Subject": "Supplementary public-data analyses and figure legends",
        "/Author": "Study authors (administrative fields pending author verification)",
    })
    with out_path.open("wb") as handle:
        writer.write(handle)


def build_markdown(specs: list[tuple[str, Path, Path]], out_path: Path) -> None:
    lines = [
        "# Supplementary Information",
        "",
        f"## {TITLE}",
        "",
        "[AUTHOR LIST]",
        "",
        "[AFFILIATIONS]",
        "",
        "The supplementary figures preserve the experimental-unit and inference boundaries stated in the main manuscript.",
        "",
    ]
    for label, _, caption in specs:
        lines.extend([f"## Supplementary Figure {label}", "", clean_caption(label, caption.read_text(encoding="utf-8")), ""])
    out_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main() -> None:
    specs = figure_specs()
    font_name = register_font()
    with tempfile.TemporaryDirectory(prefix="nintedanib_supplementary_") as temp_dir:
        temp = Path(temp_dir) / "Supplementary_Information_text_pages.pdf"
        text_pdf(specs, font_name, temp)
        build_pdf(specs, temp, OUT / "Supplementary_Information.pdf", font_name)
    build_markdown(specs, OUT / "Supplementary_Information.md")
    print(f"Built Supplementary_Information.pdf with {len(specs)} figures")


if __name__ == "__main__":
    main()
