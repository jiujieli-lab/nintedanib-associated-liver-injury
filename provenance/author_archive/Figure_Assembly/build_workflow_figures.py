from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
import csv, json, io, os, hashlib
import numpy as np
from scipy.ndimage import label

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parents[2]
ASSETS = ROOT / "assets"
OUT = ROOT / "final"
OUT.mkdir(exist_ok=True)
SHEET = ASSETS / "biomedical_icon_sheet.png"
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
BOLD = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
F = lambda size, bold=False: ImageFont.truetype(BOLD if bold else FONT, size)
INK = "#243641"
TEAL = "#2C7789"
MUTED = "#52646D"
LIGHT = "#F6F9FA"
EDGE = "#DCE5E9"
ARROW = "#7C9AA7"
HEIGHT = 380

spec = {
1: [
("Clinical injury", ["FAERS signal", "13 serum proteins"]),
("Drug + liver cells", ["Pharmacology + HepG2", "5-donor liver atlas"]),
("Target priorities", ["Calibrated integration", "PDGFRA · FLT4"]),
("Independent evaluation", ["Liver exposure challenge", "Pulmonary response context"]),
],
5: [
("Liver-cell context", ["5 human donors", "Compartment + detection"]),
("Virtual knockout", ["210 archived runs", "Raw network distances"]),
("Numerical calibration", ["Zero threshold ≤10⁻¹²", "132 numerical-zero runs"]),
("Clinical projection", ["29 nonzero target runs", "Matched-null comparison"]),
],
6: [
("Liver + clinical anchors", ["Donor-balanced networks", "13 protein weights"]),
("Calibrated proximity", ["1,999 matched null sets", "Expression + detection"]),
("Evidence integration", ["Pharmacology · network", "HepG2 · STRING"]),
("Target priorities", ["17 targets · 6 algorithms", "PDGFRA · FLT4"]),
],
7: [
("Prespecified candidates", ["17 targets", "6 algorithms"]),
("Independent mouse liver", ["Healthy 3 · disease 9", "Nintedanib 10"]),
("Animal-level analysis", ["Genome-wide + target family", "92,378 allocations"]),
("Independent comparison", ["12 measured targets", "Priority vs. hepatic effect"]),
],
}

def atomic_bytes(path, data):
    tmp = path.with_name(path.name + ".partial")
    with tmp.open("wb") as fh:
        fh.write(data)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)

def image_bytes(im, fmt):
    buf = io.BytesIO()
    kw = {"dpi": (300,300)}
    if fmt == "TIFF": kw["compression"] = "tiff_lzw"
    im.save(buf, format=fmt, **kw)
    return buf.getvalue()

# The illustration grid has visually verified, unequal gutters.
row_bounds = [(0,340),(345,660),(665,950),(950,1254)]
col_bounds = [
[(0,315),(325,635),(640,939),(940,1254)],
[(0,335),(340,630),(635,937),(940,1254)],
[(0,315),(325,630),(635,935),(940,1254)],
[(0,295),(295,637),(637,935),(937,1254)],
]
sheet = Image.open(SHEET).convert("RGBA")
icons = {}
# Two bottom-row illustrations have overlapping rectangular bounds but disjoint
# alpha components. Isolate their full connected cutouts without clipping tails.
sheet_array=np.asarray(sheet)
component_map,component_count=label(sheet_array[:,:,3]>1)
component_sizes=np.bincount(component_map.ravel())
component_candidates=[]
for component_id in np.flatnonzero(component_sizes>30000):
    if component_id==0: continue
    ys,xs=np.where(component_map==component_id)
    if ys.min()>940:
        component_candidates.append((component_id,int(xs.min()),int(xs.max()),int(ys.min()),int(ys.max())))
special_ids={}
for cid,x0,x1,y0,y1 in component_candidates:
    if 290 < x0 < 320: special_ids[(3,1)]=cid
    if 605 < x0 < 630: special_ids[(3,2)]=cid
for row in range(4):
    for col in range(4):
        y0,y1 = row_bounds[row]
        x0,x1 = col_bounds[row][col]
        im = sheet.crop((x0,y0,x1,y1))
        if (row,col) in special_ids:
            cid=special_ids[(row,col)]
            arr=sheet_array.copy()
            arr[:,:,3]=np.where(component_map==cid,arr[:,:,3],0)
            im=Image.fromarray(arr,"RGBA")
        # Trim only transparent padding. Never repaint or regenerate biology.
        box = im.getchannel("A").point(lambda a: 255 if a>12 else 0).getbbox()
        if box: im=im.crop(box)
        path=ASSETS/f"icon_{row+1}_{col+1}.png"
        atomic_bytes(path,image_bytes(im,"PNG"))
        icons[row,col]=im

def center_text(draw, x, y, text, font, fill=INK):
    width=draw.textlength(text,font=font)
    draw.text((x-width/2,y),text,font=font,fill=fill)
    return width

def draw_band(fig, width, row):
    band = Image.new("RGB", (width,HEIGHT), "white")
    d=ImageDraw.Draw(band)
    d.text((0,0),"A",font=F(54,True),fill="#111111")
    left=73
    gap=29
    usable=width-left-12
    cw=(usable-3*gap)/4
    audits=[]
    for col,(title,labels) in enumerate(spec[fig]):
        x0=round(left+col*(cw+gap));x1=round(x0+cw)
        cx=(x0+x1)/2
        # A quiet open card groups each biological stage.
        d.rounded_rectangle((x0,8,x1,365),radius=17,fill=LIGHT,outline=EDGE,width=2)
        tw=center_text(d,cx,24,title,F(31,True),TEAL)
        assert tw < cw-16,(fig,title,tw,cw)
        icon=icons[row,col].copy()
        icon.thumbnail((round(cw-42),193),Image.Resampling.LANCZOS)
        ix=round(cx-icon.width/2);iy=round(169-icon.height/2)
        band.paste(icon,(ix,iy),icon)
        labelwidths=[]
        for j,label in enumerate(labels):
            emph = label=="PDGFRA · FLT4"
            lw=center_text(d,cx,281+j*36,label,F(28,emph),TEAL if emph else INK)
            assert lw < cw-14,(fig,label,lw,cw)
            labelwidths.append(lw)
        if col < 3:
            ax=x1+gap/2
            d.line((ax-9,171,ax+8,171),fill=ARROW,width=4)
            d.polygon([(ax+8,171),(ax-1,164),(ax-1,178)],fill=ARROW)
        audits.append({"title":title,"labels":labels,"card_box":[x0,8,x1,365],"icon_box":[ix,iy,ix+icon.width,iy+icon.height],"title_width":tw,"label_widths":labelwidths})
    return band,audits

manifest=list(csv.DictReader((BASE/"work/audit/panel_crop_manifest.csv").open()))
source_cache={}
qa=[]
for row,fig in enumerate([1,5,6,7]):
    lp=BASE/f"work/panels/reassembled/Figure_{fig}_layout.json"
    layout=json.loads(lp.read_text())
    srcpath=BASE/f"work/panels/reassembled/Figure_{fig}_Reassembled_Original_Panels.tiff"
    src=Image.open(srcpath);src.load();src=src.convert("RGB")
    ap=layout["placements"][0]
    band,designaudit=draw_band(fig,ap["width"],row)
    delta=HEIGHT-ap["height"]
    dst=Image.new("RGB",(src.width,src.height+delta),"white")
    # Carry the entire lower canvas forward unchanged, including panel spacing.
    cut=ap["y"]+ap["height"]
    dst.paste(src.crop((0,0,src.width,ap["y"])),(0,0))
    dst.paste(band,(ap["x"],ap["y"]))
    dst.paste(src.crop((0,cut,src.width,src.height)),(0,cut+delta))
    lower_equal=np.array_equal(np.asarray(src)[cut:,:,:],np.asarray(dst)[cut+delta:,:,:])
    assert lower_equal
    panelchecks=[]
    for p in layout["placements"][1:]:
        r=next(r for r in manifest if r["planned_figure"]==f"Figure_{fig}" and r["planned_panel"]==p["panel"])
        original=BASE/r["source_file"]
        if str(original) not in source_cache:
            oi=Image.open(original);oi.load();source_cache[str(original)]=oi.convert("RGB")
        box=tuple(int(r[k]) for k in ["crop_left","crop_top","crop_right","crop_bottom"])
        oi=source_cache[str(original)].crop(box)
        # Prior extraction excludes only documented neighboring-panel text.
        exclusions=json.loads(r["neighbour_text_exclusions"])
        for ex in exclusions:
            eb=(ex[0]-box[0],ex[1]-box[1],ex[2]-box[0],ex[3]-box[1])
            oi.paste("white",eb)
        oi=oi.crop(tuple(json.loads(r["trim_bbox"])))
        outbox=(p["x"],p["y"]+delta,p["x"]+p["width"],p["y"]+delta+p["height"])
        actual=dst.crop(outbox)
        eq=oi.size==actual.size and np.array_equal(np.asarray(oi),np.asarray(actual))
        assert eq,(fig,p["panel"],oi.size,actual.size)
        panelchecks.append({"panel":p["panel"],"original_sha256":hashlib.sha256(oi.tobytes()).hexdigest(),"output_sha256":hashlib.sha256(actual.tobytes()).hexdigest(),"pixel_identical":eq,"output_box":outbox,"prior_neighbor_text_exclusions":exclusions})
    stem=OUT/f"Current_Figure_{fig}"
    atomic_bytes(stem.with_suffix(".tiff"),image_bytes(dst,"TIFF"))
    atomic_bytes(stem.with_suffix(".png"),image_bytes(dst,"PNG"))
    # PDF keeps the native 300dpi image geometry; no artificial resolution claim.
    pdfbuf=io.BytesIO();pdf=canvas.Canvas(pdfbuf,pagesize=(dst.width*72/300,dst.height*72/300),pageCompression=1)
    pdf.setTitle(f"Current Figure {fig}")
    pdf.drawImage(ImageReader(dst),0,0,width=dst.width*72/300,height=dst.height*72/300)
    pdf.showPage();pdf.save();atomic_bytes(stem.with_suffix(".pdf"),pdfbuf.getvalue())
    atomic_bytes(OUT/f"Current_Figure_{fig}_Workflow_A.png",image_bytes(band,"PNG"))
    qa.append({"figure":fig,"source_size":src.size,"output_size":dst.size,"source_dpi":300,"output_dpi":300,"workflow_height":HEIGHT,"vertical_translation_px":delta,"entire_lower_canvas_identical":lower_equal,"quantitative_panels":panelchecks,"layout_audit":designaudit})
(ROOT/"composite_validation.json").write_text(json.dumps(qa,indent=2))
# High-resolution contact sheet of all four revised workflow bands.
contact=Image.new("RGB",(2054,4*(HEIGHT+38)),"white")
d=ImageDraw.Draw(contact)
for row,fig in enumerate([1,5,6,7]):
    im=Image.open(OUT/f"Current_Figure_{fig}_Workflow_A.png")
    y=row*(HEIGHT+38)
    d.text((14,y+2),f"Current Figure {fig}",font=F(22,True),fill=INK)
    contact.paste(im,(14,y+35))
atomic_bytes(OUT/"Workflow_Panels_Overview.png",image_bytes(contact,"PNG"))
print(json.dumps([{k:v for k,v in x.items() if k not in ["quantitative_panels","layout_audit"]} for x in qa],indent=2))
