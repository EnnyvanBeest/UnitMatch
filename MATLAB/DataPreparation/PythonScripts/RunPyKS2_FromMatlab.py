from pathlib import Path
from pykilosort.ibl import run, ibl_pykilosort_params
import os
#import shutil

def RunPyKS(bin_file):		
	print('Starting PyKS2 now')
	bin_file = Path(bin_file)	
	print(bin_file)

	# Path management
	scratch_dir = Path(os.path.dirname(bin_file))
	ks_output_dir = scratch_dir.joinpath('output')	
	
	#shutil.rmtree(scratch_dir, ignore_errors=True)
	scratch_dir.mkdir(parents=True, exist_ok=True)
	ks_output_dir.mkdir(parents=True, exist_ok=True)

	# load parameters
	params = ibl_pykilosort_params(bin_file)
	# Run PyKS2
	run(bin_file, dir_path=scratch_dir, output_dir=ks_output_dir, **params)
	print('DONE')
	success=1
	return success


if __name__ == '__main__':
   # bin_file = 'D:/tmpdata/2021-02-24_EB001_g0_t0.imec0.ap.cbin' # #
    success = RunPyKS(bin_file)
