from pathlib import Path
import re
import spikeinterface as si
import spikeinterface.extractors as se
import spikeinterface.sorters as ss
import spikeinterface.preprocessing as spre


def run_ks4_si(bin_file):
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

    # Derive stream name (e.g., imec1.ap) from the bin filename
    imec_match = re.search(r"imec(\d+)", bin_file.name)
    stream_base = f"imec{imec_match.group(1)}" if imec_match else "imec0"
    if imec_match is None:
        print(f"Warning: could not parse imec index from {bin_file.name}; defaulting to {stream_base}")
    stream_name = f"{stream_base}.ap"

    # Load SpikeGLX recording. `read_spikeglx` accepts either the .ap.bin
    # full path or the folder containing the .bin/.meta pair. Try the
    # .bin first (most explicit), then fall back to the folder.
    print(bin_file)
    recording = se.read_spikeglx(bin_file.parent, stream_name=stream_name, load_sync_channel=False)
  
    # Probe handling: prefer an explicit `probe_file` if provided; otherwise
    # rely on probe metadata found in the .meta file. Be defensive: `get_probe`
    # may raise if no probe is present.
    probe = recording.get_probe()

    # Preprocessing chain
    recording = spre.phase_shift(recording)  # correct for time delay between recording channels
    bad_channel_ids, channel_labels = spre.detect_bad_channels(recording)
    recording = recording.remove_channels(bad_channel_ids)
    recording = spre.highpass_filter(recording)  # correct for time delay between recording channels

    print('bad_channel_ids', bad_channel_ids)

    # Run Kilosort 4
    print("Starting KS4 now")
    sorting_ks4 = ss.run_sorter("kilosort4", recording, output_dir, remove_existing_folder=True)
    
    print("Kilosort4 output:", output_dir)
    success = sorting_ks4 is not None
    return success


if __name__ == "__main__":

    # Example:
    # success = run_ks4_si(r"D:/tmpdata/2025-09-22_EB053_3_g0_g0/2025-09-22_EB053_3_g0_g0_imec1/2025-09-22_EB053_3_g0_g0_t0.imec1.ap.bin")
    
    success = run_ks4_si(bin_file)
