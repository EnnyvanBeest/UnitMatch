from pathlib import Path
from pykilosort import run, add_default_handler, np1_probe, np2_probe, np2_4shank_probe, neuropixel_probe_from_metafile
import glob
import os

def RunPyKS(ThisFile):
	print(ThisFile)
	data_paths = list()
	tmpdatapath = f'{ThisFile}/*.*bin'
	for name in glob.glob(tmpdatapath):
		print(Path(name))
		data_paths.append(Path(name))
		
	ProbeType = neuropixel_probe_from_metafile(data_paths[0])
	print('Starting PyKS2 now')	
	# run(data_path, probe=ProbeType, dir_path=OutputDir) #this runs from server, not ideal
	run(data_paths, probe=ProbeType) #This runs from local
	print('DONE')
	success=1
	return success


if __name__ == '__main__':
    #ThisFile = r'H:/tmpdata' #'//znas/Subjects/' # #
    success = RunPyKS(ThisFile)
