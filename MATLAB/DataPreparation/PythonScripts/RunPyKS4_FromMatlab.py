from pathlib import Path
#from tqdm import tqdm
import os 
#from kilosort.utils import download_probes
from kilosort import run_kilosort

def RunKS4(bin_file,probe_file):		
	print('Starting KS4 now')
	bin_file = Path(bin_file)	
	print(bin_file)
   
    # Path management
	scratch_dir = Path(os.path.dirname(bin_file))
	
	#shutil.rmtree(scratch_dir, ignore_errors=True)
	scratch_dir.mkdir(parents=True, exist_ok=True)
	
	# Download channelmaps
	#download_probes
    # 'dminx':400
    # ,'nearest_templates':50
	settings = {'data_dir':bin_file.parent, 'n_chan_bin':385, 'probe_path':probe_file,'x_centers':10}
	ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = run_kilosort(settings=settings, filename = bin_file)
	
	print('DONE')
	success=1
	return success


if __name__ == '__main__':
	#bin_file = 'D:/tmpdata/2024-01-24_EB036_OrangeRigs_g0_t0.imec0.ap.bin'
	#probe_file = 'D:/tmpdata/2024-01-24_EB036_OrangeRigs_g0_t0.imec0_kilosortChanMap.mat'
	success = RunKS4(bin_file,probe_file)
