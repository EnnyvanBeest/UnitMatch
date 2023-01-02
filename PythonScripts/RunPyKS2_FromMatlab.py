from pathlib import Path
from pykilosort import run, add_default_handler, np1_probe, np2_probe, np2_4shank_probe, neuropixel_probe_from_metafile
#import glob
#import shutil
#import os

def RunPyKS(ThisFile):		
	print('Starting PyKS2 now')
	print(ThisFile)
	ProbeType = neuropixel_probe_from_metafile(ThisFile)
	# run(data_path, probe=ProbeType, dir_path=OutputDir) #this runs from server, not ideal
	run(ThisFile, probe=ProbeType, template_snapshots = [0.2, 0.5, 0.8]) #This runs from local
	print('DONE')
	success=1
	return success


if __name__ == '__main__':
   # ThisFile = 'D:/tmpdata/2021-02-24_EB001_g0_t0.imec0.ap.cbin' # #
    success = RunPyKS(ThisFile)