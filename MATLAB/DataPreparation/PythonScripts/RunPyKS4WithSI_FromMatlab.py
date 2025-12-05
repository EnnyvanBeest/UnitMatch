from pathlib import Path
import spikeinterface as si
import spikeinterface.extractors as se
import spikeinterface.sorters as ss
import spikeinterface.preprocessing as spre


def run_ks4_si(bin_file, probe_file=None):
    """
    Run SpikeInterface preprocessing + Kilosort4 on a SpikeGLX ap.bin.

    Parameters
    ----------
    bin_file : str or Path
        Path to the .ap.bin file (SpikeGLX format). The matching .meta is required.
    probe_file : str or Path, optional

    Returns
    -------
    Success : Boolean
    """
    bin_file = Path(bin_file)
    meta_file = bin_file.with_suffix(".meta")
    if not meta_file.exists():
        raise FileNotFoundError(f"Meta file not found: {meta_file}")

    output_dir = bin_file.parent / "kilosort4"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load SpikeGLX recording from the folder containing the bin/meta pair
    rec0 = se.read_spikeglx(bin_file.parent, load_sync_channel=False)  # auto-detects the matching .bin

    # Automatically have the probe from the meta; allow overriding with probe_file if provided
    rec0.get_probe().to_dataframe()
    if probe_file:
        try:
            rec0 = rec0.set_probe(probe_file)
            print(f"Loaded probe from {probe_file}")
        except Exception as exc:  # fallback to existing probe metadata
            print(f"Warning: failed to load probe file {probe_file}: {exc}")

    # Preprocessing chain
    recording = spre.phase_shift(rec0)  # correct for time delay between recording channels
    recording = spre.highpass_filter(recording)  # highpass
    bad_channel_ids, _ = spre.detect_bad_channels(recording)
    if len(bad_channel_ids) > 0:
        print("bad_channel_ids", bad_channel_ids)
        recording = recording.remove_channels(bad_channel_ids)
    else:
        print("bad_channel_ids []")
    recording = spre.common_reference(recording, operator="median", reference="global")

    # Run Kilosort 4
    print("Starting KS4 now")
    sorter_params = si.get_default_sorter_params("kilosort4")
     
    sorting = ss.run_sorter(
        "kilosort4",
        recording=recording,
        output_folder=output_dir,
        remove_existing_folder=True,
        docker_image=True,
        verbose=True,
        **sorter_params,
    )

    print("Kilosort4 output:", output_dir)
    success = sorting is not None
    return success


if __name__ == "__main__":

    # Example:
    success = run_ks4_si(r"D:/tmpdata/2025-09-22_EB053_3_g0_g0/2025-09-22_EB053_3_g0_g0_imec1/2025-09-22_EB053_3_g0_g0_t0.imec1.ap.bin",r"D:/tmpdata/2025-09-22_EB053_3_g0_g0_t0.imec1_kilosortChanMap.mat")
    
    #success = run_ks4_si(bin_file,probe_file)

