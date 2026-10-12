from pathlib import Path
import argparse,json
from PIL import Image

def main():
    p=argparse.ArgumentParser();p.add_argument('image',type=Path);a=p.parse_args()
    image=Image.open(a.image);image.load()
    if image.mode!='RGBA':raise ValueError('Expected actual RGBA transparent render')
    alpha=image.getchannel('A');bounds=alpha.getbbox()
    if bounds is None or alpha.getextrema()[0]!=0:raise ValueError('Empty render or no transparent background')
    if bounds[0]<=1 or bounds[1]<=1 or bounds[2]>=image.width-1 or bounds[3]>=image.height-1:raise ValueError('Visible image may be clipped at border')
    white=Image.new('RGBA',image.size,'white');white.alpha_composite(image)
    target=a.image.with_name(a.image.stem.replace('_transparent','')+'_white.png');white.convert('RGB').save(target)
    report={'passed':True,'transparent_image':str(a.image.resolve()),'white_image':str(target.resolve()),'size':list(image.size),'visible_bounds':list(bounds),'transparent_background':True}
    target.with_suffix('.validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
