from pathlib import Path
from PIL import Image
import shutil, json, hashlib, io, struct, gc

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent/'repaired_quantitative'
SRC=ROOT/'work/source/转化医学投稿20260922'
OUT.mkdir(exist_ok=True)
Image.MAX_IMAGE_PIXELS=None
sources=[(f'Figure_{n}',SRC/f'Figures/Figure {n}.tif') for n in (2,3)]
sources += [(f'Figure_S{n}',SRC/f'supplementary Figures_S1-S10/Figure_S{n}_1500dpi.tiff') for n in range(1,11)]
sources += [('Figure_S11',ROOT/'work/new_virtual/Figure_S11_Virtual_Perturbation_1200dpi.tiff')]
rows=[]
for name,p in sources:
    dst=OUT/f'{name}.tiff'
    with Image.open(p) as im:
        size=im.size; dpi=tuple(float(x) for x in im.info.get('dpi',(300,300)))
        mode=im.mode
        if mode=='RGB':
            shutil.copyfile(p,dst)
            mode_out='RGB'
            # Source byte identity preserves every pixel, color profile, and resolution tag.
            digest=hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
            assert digest==hashlib.file_digest(dst.open('rb'),'sha256').hexdigest()
        else:
            if mode=='RGBA':
                assert im.getchannel('A').getextrema()==(255,255)
            rgb=im.convert('RGB'); buffer=io.BytesIO()
            rgb.save(buffer,format='TIFF',compression='tiff_lzw',dpi=dpi)
            dst.with_suffix('.tmp').write_bytes(buffer.getvalue())
            dst.with_suffix('.tmp').replace(dst)
            buffer.close(); mode_out='RGB'; del rgb
        # Preview only: deliver native TIFF separately; previews are not source replacements.
        im.thumbnail((2100,2700),Image.Resampling.LANCZOS)
        im.convert('RGB').save(OUT/f'{name}.png',dpi=(300,300))
    header=dst.read_bytes()[:8]
    offset=struct.unpack('<I' if header[:2]==b'II' else '>I',header[4:8])[0]
    assert offset>0
    row={'figure':name,'source':str(p.relative_to(ROOT)),'output':str(dst.relative_to(ROOT)),
         'native_pixels':size,'native_dpi':dpi,'source_mode':mode,'output_mode':mode_out,
         'first_ifd_offset':offset,'bytes':dst.stat().st_size,'native_pixels_preserved':True,
         'source_byte_identical':mode=='RGB'}
    rows.append(row);print(name,size,dpi,flush=True);gc.collect()
(OUT/'Original_Figure_Export_Manifest.json').write_text(json.dumps(rows,indent=2))
