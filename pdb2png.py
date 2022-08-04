#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 09:50:13 2022

@author: mheinzinger
"""
from pathlib import Path
import sys

if "pymol" not in  "\t".join(sys.path):
    # corresponding path to your pymol installation
    sys.path.append("/home/mheinzinger/work_space/pymol/lib/python3.9/site-packages")
    import __main__
    __main__.pymol_argv = ['pymol','-qc']
    import pymol
    pymol.finish_launching()
else:
    import pymol

def write_png(pdb_file, png_file):
    pdb_name =pdb_file.name.split('.')[0]
    pymol.cmd.load(pdb_file, pdb_name)
    pymol.cmd.disable("all")
    pymol.cmd.enable(pdb_name)
    print(pymol.cmd.get_names())
    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 0)
    pymol.cmd.color('red', 'ss h')
    pymol.cmd.color('yellow', 'ss s')
    pymol.cmd.png(str(png_file))
    #pymol.cmd.quit()

def main():
    # directory
    pdb_dir = Path("./mysite/rank_1_10char")
    png_dir = Path("./mysite/rank_1_10char_pngs")
    for pdb_p in Path(pdb_dir).glob('**/*.pdb'):
        png_p = png_dir / pdb_p.name.replace(".pdb",".png")
        write_png(pdb_p,png_p)
    
    

if __name__ == "__main__":
    main()