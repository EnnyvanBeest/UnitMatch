from pathlib import Path
from pykilosort.ibl import run_spike_sorting_ibl, ibl_pykilosort_params
import glob
import os

def RunPyKS(ThisFile):
	print(ThisFile)
	data_paths = list()
	tmpdatapath = f'{ThisFile}/*.*bin'
	for name in glob.glob(tmpdatapath):
		print(Path(name))
		data_paths.append(Path(name))
		

	# Path management
	scratch_dir = Path(os.path.dirname(ThisFile))
	ks_output_dir = scratch_dir.joinpath('output')	
	
	#shutil.rmtree(scratch_dir, ignore_errors=True)
	scratch_dir.mkdir(parents=True, exist_ok=True)
	ks_output_dir.mkdir(parents=True, exist_ok=True)
	
	# Load parameters
	params = ibl_pykilosort_params(data_paths)
	print('Starting PyKS2 now')	
	# Run PyKS2
	run(data_paths, dir_path=scratch_dir, output_dir=ks_output_dir, **params)

	# run(data_path, probe=ProbeType, dir_path=OutputDir) #this runs from server, not ideal
	run_spike_sorting_ibl(data_paths, delete=DELETE, log_level='INFO', params=params) #This runs from local
	print('DONE')
	success=1
	return success


if __name__ == '__main__':
    #ThisFile = r'H:/tmpdata' #'//znas/Subjects/' # #
    success = RunPyKS(ThisFile)
