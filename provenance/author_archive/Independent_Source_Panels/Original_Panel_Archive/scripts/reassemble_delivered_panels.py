"""Reassemble the six retained main figures from delivered, lossless panel TIFFs.
Run from any working directory after extracting this panel archive.
"""
from pathlib import Path
from PIL import Image
import json,csv
Image.MAX_IMAGE_PIXELS=None
PANEL_ROOT=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader(open(PANEL_ROOT/'panel_crop_manifest.csv')))
for lp in sorted((PANEL_ROOT/'reassembled').glob('Figure_*_layout.json')):
 spec=json.loads(lp.read_text())
 if spec['new_figure'] not in ['Figure_1','Figure_2','Figure_3','Figure_5','Figure_6','Figure_7']: continue
 canvas=Image.new('RGB',tuple(spec['canvas_px']),'white')
 for loc in spec['placements']:
  row=next(r for r in rows if r['original_figure']==spec['original_figure'] and r['original_panel']==loc['panel'])
  p=PANEL_ROOT/'independent'/row['original_figure']/(row['original_figure']+'_'+row['original_panel']+'.tiff')
  with Image.open(p) as im:
   a=im.crop(tuple(json.loads(row['trim_bbox'])))
   assert a.size==(loc['width'],loc['height'])
   canvas.paste(a,(loc['x'],loc['y']))
 dst=PANEL_ROOT/'reassembled'/(spec['new_figure']+'_Reassembled_Original_Panels.tiff')
 canvas.save(dst,compression='tiff_lzw',dpi=tuple(spec['dpi']))
 print(dst)
