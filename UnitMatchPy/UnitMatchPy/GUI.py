from tkinter import *
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk 
import numpy as np
from matplotlib import rcParams
import os
import pickle


def precalculate_all_acgs(clus_info, param, save_path=None, bin_size=0.001, max_lag=0.05):
    """
    Pre-calculate and save autocorrelograms for all units to speed up GUI loading
    
    Parameters
    ----------
    clus_info : dict
        Cluster information dictionary containing unit metadata
    param : dict
        Parameters dictionary containing KS_dirs and other settings
    save_path : str, optional
        Path to save the ACG data. If None, saves in current directory as 'acg_cache.pkl'
    bin_size : float
        Bin size in seconds (default 1ms)
    max_lag : float
        Maximum lag in seconds (default 50ms)
        
    Returns
    -------
    acg_cache : dict
        Dictionary mapping unit_id to (autocorr, bin_centers) tuples
    """
    print("Pre-calculating autocorrelograms for all units...")
    
    if save_path is None:
        save_path = 'acg_cache.pkl'
    
    acg_cache = {}
    n_units = len(clus_info.get('session_id', []))
    
    for unit_id in range(n_units):
        try:
            # Get spike times for this unit
            spike_times = get_spike_times_for_unit_precalc(unit_id, clus_info, param)
            
            if len(spike_times) > 1:
                # Compute ACG
                autocorr, bin_centers = compute_acg_precalc(spike_times, bin_size, max_lag)
                acg_cache[unit_id] = (autocorr, bin_centers)
                
                if (unit_id + 1) % 50 == 0:  # Progress update every 50 units
                    print(f"Processed {unit_id + 1}/{n_units} units")
            else:
                acg_cache[unit_id] = (np.array([]), np.array([]))
                
        except Exception as e:
            print(f"Error computing ACG for unit {unit_id}: {e}")
            acg_cache[unit_id] = (np.array([]), np.array([]))
    
    # Save to file
    try:
        with open(save_path, 'wb') as f:
            pickle.dump(acg_cache, f)
        print(f"ACG cache saved to {save_path}")
    except Exception as e:
        print(f"Error saving ACG cache: {e}")
    
    return acg_cache


def get_spike_times_for_unit_precalc(unit_id, clus_info, param):
    """
    Get spike times for a specific unit (used in pre-calculation)
    
    Parameters
    ----------
    unit_id : int
        UnitMatch unit ID
    clus_info : dict
        Cluster information dictionary
    param : dict
        Parameters dictionary
        
    Returns
    -------
    spike_times : array
        Spike times in seconds for the unit
    """
    try:
        # Get session information for this unit
        session_id = clus_info['session_id'][unit_id]
        original_id = clus_info['original_ids'][unit_id]
        
        # Get the Kilosort directory for this session
        if 'KS_dirs' in param:
            ks_dir = param['KS_dirs'][session_id]
        else:
            return np.array([])
        
        # Load Kilosort spike times and cluster assignments
        spike_times_path = os.path.join(ks_dir, 'spike_times.npy')
        spike_clusters_path = os.path.join(ks_dir, 'spike_clusters.npy')
        
        if os.path.exists(spike_times_path) and os.path.exists(spike_clusters_path):
            # Load spike times (in samples) and cluster assignments
            spike_times_samples = np.load(spike_times_path).flatten()
            spike_clusters = np.load(spike_clusters_path).flatten()
            
            # Get sampling rate - try to load from params.py or use default
            try:
                params_path = os.path.join(ks_dir, 'params.py')
                sample_rate = 30000  # Default sample rate
                if os.path.exists(params_path):
                    # Parse params.py for sample_rate
                    with open(params_path, 'r') as f:
                        params_content = f.read()
                    # Extract sample_rate
                    for line in params_content.split('\n'):
                        if 'sample_rate' in line and '=' in line:
                            sample_rate = float(line.split('=')[1].strip())
                            break
            except:
                sample_rate = 30000  # Fallback
            
            # Get spike times for this specific unit
            unit_mask = spike_clusters == original_id
            unit_spike_times = spike_times_samples[unit_mask] / sample_rate  # Convert to seconds
            
            return unit_spike_times
        else:
            return np.array([])
            
    except Exception as e:
        return np.array([])


def compute_acg_precalc(spike_times, bin_size=0.001, max_lag=0.05):
    """
    Compute autocorrelogram for a unit (used in pre-calculation)
    
    Parameters
    ----------
    spike_times : array
        Spike times in seconds
    bin_size : float
        Bin size in seconds (default 1ms)
    max_lag : float
        Maximum lag in seconds (default 50ms)
        
    Returns
    -------
    autocorr : array
        Autocorrelogram counts
    bin_centers : array
        Bin centers in seconds
    """
    if len(spike_times) < 2:
        return np.array([]), np.array([])
    
    # Create bins centered around 0
    bins = np.arange(-max_lag, max_lag + bin_size, bin_size)
    bin_centers = bins[:-1] + bin_size/2
    
    # Calculate autocorrelogram
    autocorr = np.zeros(len(bin_centers))
    
    # For each spike, find all other spikes within max_lag
    for i, spike_time in enumerate(spike_times):
        # Find spikes within the lag window
        time_diffs = spike_times - spike_time
        valid_diffs = time_diffs[(np.abs(time_diffs) <= max_lag) & (time_diffs != 0)]
        
        # Bin the time differences
        hist, _ = np.histogram(valid_diffs, bins)
        autocorr += hist
    
    # Convert to firing rate (Hz)
    total_time = spike_times[-1] - spike_times[0] if len(spike_times) > 1 else 1
    n_spikes = len(spike_times)
    autocorr = autocorr / (bin_size * n_spikes * total_time)
    
    return autocorr, bin_centers


def load_acg_cache(cache_path='acg_cache.pkl'):
    """
    Load pre-calculated ACG cache from file
    
    Parameters
    ----------
    cache_path : str
        Path to the ACG cache file
        
    Returns
    -------
    acg_cache : dict or None
        Dictionary mapping unit_id to (autocorr, bin_centers) tuples, or None if loading failed
    """
    try:
        with open(cache_path, 'rb') as f:
            acg_cache = pickle.load(f)
        print(f"ACG cache loaded from {cache_path}")
        return acg_cache
    except Exception as e:
        print(f"Could not load ACG cache from {cache_path}: {e}")
        return None


def run_GUI():
    """
    This function runs the GUI, allowing the user to look at the result and manually curate the matches.

    Returns
    -------
    List
        lists for the manually curated matches and non-matches
    """
    global CV_tkinter
    global root
    global entry_a
    global entry_b
    global session_entry_a
    global session_entry_b
    global match_idx
    global frame_table
    global score_table
    global avg_waveform_plot
    global trajectory_plot
    global bayes_label
    global original_id_label
    global raw_waveform_plot
    global hist_plot
    global is_match
    global not_match
    global option_a
    global option_b
    global entry_frame
    global toggle_raw_val
    global toggle_UM_score_val
    global unit_legend_plot
    global hist_legend_plot
    global acg_plot
    global toggle_acg_val
    global acg_cache

    # Try to load pre-calculated ACG cache
    acg_cache = load_acg_cache()

    rcParams.update({'figure.autolayout': True})
    rcParams.update({'font.size': 14})
    rcParams.update({'font.family': 'DejaVu Sans'})
    color = 'white'
    rcParams['text.color'] = color
    rcParams['axes.labelcolor'] = color
    rcParams['xtick.color'] = color
    rcParams['ytick.color'] = color

    
    np.set_printoptions(suppress=True)
    is_match = []
    not_match = []
    root = Tk()
    # Get screen width and height
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    # Set window size to fit the screen with some padding
    window_width = screen_width - 100
    window_height = screen_height - 100
    root.geometry(f"{window_width}x{window_height}+50+50")
    
    # Configure column weights to prevent plot squashing
    for i in range(7):  # Configure columns 0-6
        root.columnconfigure(i, weight=1)
    # Configure rows to expand properly
    for i in range(10):  # Configure rows 0-9
        root.rowconfigure(i, weight=1)

    # downloaded theme from https://sourceforge.net/projects/tcl-awthemes/
    theme_path_rel = os.path.join('TkinterTheme', 'awthemes-10.4.0')
    theme_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), theme_path_rel)

    root.tk.call('lappend', 'auto_path', theme_path)
    root.tk.call('package', 'require', 'awdark')
    s = ttk.Style(root)
    s.theme_use('awdark')
    
    # Configure fonts for all TTK widgets to use DejaVu Sans
    s.configure('.', font=('DejaVu Sans', 12))
    s.configure('TLabel', font=('DejaVu Sans', 12))
    s.configure('TButton', font=('DejaVu Sans', 12))
    s.configure('TEntry', font=('DejaVu Sans', 10))
    s.configure('TCombobox', font=('DejaVu Sans', 10))
    s.configure('TCheckbutton', font=('DejaVu Sans', 12))
    s.configure('TRadiobutton', font=('DejaVu Sans', 12))
    s.configure('TLabelFrame.Label', font=('DejaVu Sans', 12, 'bold'))
    
    root.title('UMPy - Manual Curation')
    #root.geometry('800x800')

    # Construct the file path to the icon
    icon_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'GUI_icon.png')

    # Debugging print statements
    print(f"Icon path: {icon_path}")
    print(f"File exists: {os.path.exists(icon_path)}")

    # Load the icon
    try:
        icon = PhotoImage(file=icon_path)
        # Set the icon photo
        root.iconphoto(False, icon)
    except Exception as e:
        print(f"Error loading icon: {e}")


    background = ttk.Frame(root)
    background.place(x=0, y=0, relwidth=1.0, relheight=1.0)

    frame_table = ttk.LabelFrame(root)
    score_table = ttk.LabelFrame(root)
    avg_waveform_plot = Canvas(root) 
    trajectory_plot = Canvas(root)
    bayes_label = ttk.Label(root)
    original_id_label = ttk.Label(root)
    raw_waveform_plot = Canvas(root)
    hist_plot = Canvas(root)

    #Unit entry
    ######################################################################################
    #Will have Unit A - color green and unit B- color Blue
    entry_frame = ttk.LabelFrame(root, text = 'Select Units')
    label_a = ttk.Label(entry_frame, text = 'Unit A')
    label_a.configure(foreground='green')
    label_b = ttk.Label(entry_frame, text = 'Unit B')
    label_b.configure(foreground='blue')

    #select the session
    sessions_list = np.arange(1,param['n_sessions']+1).tolist()
    session_entry_a = ttk.Combobox(entry_frame, value = sessions_list, width = 2)
    session_entry_b = ttk.Combobox(entry_frame, value = sessions_list, width = 2)
    session_entry_a.set(1) #Start wiht session 1,2 if more than 1 session given
    if len(sessions_list) == 1:
        session_entry_b.set(1)
    else:
        session_entry_b.set(2)
    label_session_a = ttk.Label(entry_frame, text = 'Session No.')
    label_session_b = ttk.Label(entry_frame, text = 'Session No.')


    #select CV
    CV_options = [('Avg', 0),  ('(1,2)', 1), ('(2,1)',2)]
    CV_tkinter = IntVar()
    CV_tkinter.set(0)
    label_cv = ttk.Label(entry_frame, text = 'Select the cv option')
    for i, option in enumerate(CV_options):
        RadioCV = ttk.Radiobutton(entry_frame, text = option[0], value =option[1], variable = CV_tkinter, command = update_unit_cv).grid(row = i+1, column = 0)

    #selecting the unit
    session_a = int(session_entry_a.get())
    session_b = int(session_entry_b.get()) 
    CV = get_cv_option()
    CV_option = CV_tkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV values
        tmp_idx_a = np.argwhere((session_switch[session_a -1] <= matches_avg[:,0]) * (matches_avg[:,0] < session_switch[session_a]) == True ).squeeze()
        tmp_idx_b = np.argwhere((session_switch[session_b -1] <= matches_avg[:,1]) * (matches_avg[:,1] < session_switch[session_b]) == True).squeeze()
        in_both = np.isin(tmp_idx_a, tmp_idx_b)
        option_a = matches_avg[tmp_idx_a[in_both],:].tolist()

        entry_a = ttk.Combobox(entry_frame, values = option_a, width = 10)
        entry_a.set(option_a[0][0])

        option_b = np.flip(np.argsort(output_avg[int(entry_a.get()),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()

    else :    
        tmp_idx_a = np.argwhere((session_switch[session_a -1] <= matches_GUI[CV_option][:,0]) * ( matches_GUI[CV_option][:,0] < session_switch[session_a]) == True ).squeeze()
        tmp_idx_b = np.argwhere((session_switch[session_b -1] <= matches_GUI[CV_option][:,1]) * ( matches_GUI[CV_option][:,1] < session_switch[session_b]) == True).squeeze()
        in_both = np.isin(tmp_idx_a, tmp_idx_b)

        option_a =  matches_GUI[CV_option][tmp_idx_a[in_both],:].tolist()
        if CV_option == 0:
            option_a = sorted(option_a)

        entry_a = ttk.Combobox(entry_frame, values = option_a, width = 10)
        entry_a.set(option_a[0][0])


        option_b = np.flip(np.argsort(output_GUI[CV_option][int(entry_a.get()),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()

    ##################
    entry_b = ttk.Combobox(entry_frame, values = option_b, width = 10 )
    entry_b.set(option_b[0])
    entry_a.bind('<<ComboboxSelected>>', update_units)
    entry_b.bind('<<ComboboxSelected>>', update)

    session_entry_a.bind('<<ComboboxSelected>>', update_unit_entryA)
    session_entry_b.bind('<<ComboboxSelected>>', update_unit_entryB)

    #adding a button which swaps unit A and B
    swap_button = ttk.Button(entry_frame, text = 'Swap Units', command = swap_units)

    #Calculate the score histograms
    #for each CV pair
    #make global variable so the functions can access the histograms.
    global hist_names_avg
    global hist_names_12
    global hist_names_21
    global hist_avg
    global hist_12
    global hist_21
    global hist_matches_avg
    global hist_matches_12
    global hist_matches_21

    hist_names_avg, hist_avg, hist_matches_avg =  get_score_histograms(scores_to_include_avg, (output_avg > match_threshold))
    hist_names_12, hist_12, hist_matches_12 =  get_score_histograms(scores_to_include_GUI[0], (output_GUI[0] > match_threshold))
    hist_names_21, hist_21, hist_matches_21 =  get_score_histograms(scores_to_include_GUI[1], (output_GUI[1] > match_threshold))


    #place the widgets on the EntryFrame
    label_cv.grid(row = 0, column = 0)
    label_a.grid(row = 0, column = 1)
    label_b.grid(row = 0, column = 3)
    label_session_a.grid(row = 1, column = 1)
    session_entry_a.grid(row = 1, column = 2, padx = 15)
    label_session_b.grid(row = 1, column = 3)
    session_entry_b.grid(row = 1, column = 4, padx = 15)
    entry_a.grid(row = 2, column = 1, columnspan = 2, stick = 'WE', padx = 5)
    entry_b.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)
    swap_button.grid(row = 3, column = 2, columnspan = 2, sticky = 'WE')
    ######################################################################################


    #MatchButtons
    ######################################################################################
    match_button = ttk.Button(root, text = 'Set as Match', command = set_match)
    non_match_button = ttk.Button(root, text='Set as Non Match', command = set_not_match)

    #Toggle Plots
    ######################################################################################
    toggle_raw_val = BooleanVar()
    toggle_UM_score_val = BooleanVar()
    toggle_acg_val = BooleanVar()
    toggle_raw_val.set(False)
    toggle_UM_score_val.set(False)
    toggle_acg_val.set(True)  # Default to hiding ACGs
    toggle_raw_plot = ttk.Checkbutton(root, text ='Hide Raw Data', variable = toggle_raw_val)
    toggle_UM_score_plot = ttk.Checkbutton(root, text ='Hide UM Score Histograms', variable = toggle_UM_score_val)
    toggle_acg_plot = ttk.Checkbutton(root, text ='Hide Autocorrelograms (ACGs)', variable = toggle_acg_val)



    # Set up Key-Board shortcuts
    root.bind_all('u', update)
    root.bind_all('<Return>', update)
    root.bind_all('<Right>', next_pair)
    root.bind_all('<Left>', previous_pair)
    root.bind_all('<Up>', up_options_b_list)
    root.bind_all('<Down>', down_options_b_list)
    root.bind_all('q', set_match)
    root.bind_all('m', set_match)
    root.bind_all('e', set_not_match)
    root.bind_all('n', set_not_match)

    #Grid the units
    entry_frame.grid(row = 2, column = 0, pady=5, padx = 5)
    match_button.grid(row = 1, column = 0, sticky = 'W',  padx = 10, pady = 5)
    non_match_button.grid(row = 1, column = 1, sticky = 'W', padx = 10, pady = 5)
    toggle_UM_score_plot.grid(row = 1, column = 2, sticky = 'W',  padx = 10, pady = 5)
    toggle_raw_plot.grid(row = 1, column = 3, sticky = 'W',  padx = 10, pady = 5)
    toggle_acg_plot.grid(row = 1, column = 4, sticky = 'W',  padx = 10, pady = 5)
    
    # Create visual legend plots outside plots
    global unit_legend_plot
    global hist_legend_plot
    create_unit_legend()
    create_hist_legend()
    
    unit_legend_plot.grid(row = 0, column = 1, columnspan = 2, sticky = 'W', padx = 10, pady = 5)
    hist_legend_plot.grid(row = 0, column = 3, columnspan = 4, sticky = 'W', padx = 10, pady = 5)

    # Configure grid weights to auto-adjust subpanels
    root.grid_rowconfigure(0, weight=1)
    root.grid_columnconfigure(0, weight=1)
    entry_frame.grid_rowconfigure(0, weight=1)
    entry_frame.grid_columnconfigure(0, weight=1)

    update(None)
    match_idx = 0

    root.mainloop()

    #set default plot color back to black
    color = 'black'
    rcParams['text.color'] = color
    rcParams['axes.labelcolor'] = color
    rcParams['xtick.color'] = color
    rcParams['ytick.color'] = color

    return is_match, not_match, matches_GUI

def process_info_for_GUI(output, match_threshold_in, scores_to_include, total_score, amplitude_in, spatial_decay_in,
                         avg_centroid_in, avg_waveform_in, avg_waveform_per_tp_in, wave_idx_in, max_site_in, max_site_mean_in, 
                         waveform_in, within_session_in, channel_pos_in, clus_info_in, param_in):
    """
    This function:
    1 - passes data to the GUI
    2 - processes the data so it is in a better form for the GUI
    """
    global matches_avg
    global matches_GUI
    global output_avg
    global output_GUI
    global scores_to_include_avg
    global scores_to_include_GUI
    global amplitude
    global amplitude_avg
    global spatial_decay
    global spatial_decay_avg
    global avg_centroid
    global avg_centroid_avg
    global avg_waveform
    global avg_waveform_avg
    global avg_waveform_per_tp
    global avg_waveform_per_tp_avg
    global wave_idx
    global max_site
    global max_site_mean
    global waveform
    global clus_info
    global session_switch
    global param
    global match_threshold
    global within_session
    global channel_pos

    amplitude = amplitude_in
    spatial_decay = spatial_decay_in
    avg_centroid = avg_centroid_in
    avg_waveform = avg_waveform_in
    avg_waveform_per_tp = avg_waveform_per_tp_in
    wave_idx = wave_idx_in
    max_site = max_site_in
    max_site_mean = max_site_mean_in
    waveform = waveform_in
    clus_info = clus_info_in
    session_switch = clus_info['session_switch']
    param = param_in
    match_threshold = match_threshold_in
    within_session = within_session_in
    channel_pos = channel_pos_in

    output_threshold = np.zeros_like(output)
    output_threshold[output > match_threshold] = 1

    matches = np.argwhere(output_threshold == 1) # need all matches including same session for GUI


    #for the code it is helpful for matches to contain both (session 1, session 2) and (session 2, session 1)
    # so when changing Unit A/B session it can find matches for all permutations
    matches_12_part_1 = np.argwhere(np.tril(output) > match_threshold) 
    matches_12_part_2 = np.argwhere(np.tril(output).T > match_threshold)
    matches_12 = np.unique(np.concatenate((matches_12_part_1,matches_12_part_2)), axis = 0)


    matches_21_part_1 = np.argwhere(np.triu(output) > match_threshold) 
    matches_21_part_2 = np.argwhere(np.triu(output).T > match_threshold)
    matches_21 = np.unique(np.concatenate((matches_21_part_1,matches_21_part_2)), axis = 0)

    matches_GUI = [matches_12, matches_21]

    output_GUI1_part_1 = np.tril(output)
    output_GUI1_part_2 = np.tril(output).T
    np.fill_diagonal(output_GUI1_part_2, 0)
    output_GUI1 = output_GUI1_part_1 + output_GUI1_part_2

    output_GUI2_part_1 = np.triu(output)
    output_GUI2_part_2 = np.triu(output).T
    np.fill_diagonal(output_GUI2_part_2, 0)
    output_GUI2 = output_GUI2_part_1 + output_GUI2_part_2

    output_GUI = [output_GUI1, output_GUI2]

    scores_to_include_GUI = []
    scores_to_include_12 = {}
    scores_to_include_21 = {}
    for key, value in scores_to_include.items():
        tmp1 = np.tril(value)
        tmp2 = np.tril(value).T
        np.fill_diagonal(tmp2, 0)

        tmp3 = np.triu(value)
        tmp4 = np.triu(value).T
        np.fill_diagonal(tmp4, 0)

        scores_to_include_12[key] = tmp1 + tmp2 
        scores_to_include_21[key] = tmp3 + tmp4  

    scores_to_include_GUI = [scores_to_include_12, scores_to_include_21]


    #getting avg CV data
    # for the Scores where can do (X + X.T)/2 and take upper triangular part
    total_score_avg = np.triu( (total_score + total_score.T) / 2)

    scores_to_include_avg = {}
    for key, value in scores_to_include.items():
        scores_to_include_avg[key] = (value + value.T) / 2

    output_avg = (output_GUI[0] + output_GUI[1]) / 2

    matches_avg_part_1 = np.argwhere( output_avg > match_threshold )
    matches_avg_part_2 = np.argwhere( output_avg.T > match_threshold )
    matches_avg = np.unique(np.concatenate((matches_avg_part_1,matches_avg_part_2)), axis = 0)

    # or an simply average over both CV
    amplitude_avg = np.mean(amplitude, axis = -1)
    spatial_decay_avg = np.mean(spatial_decay, axis = -1)
    avg_centroid_avg = np.mean(avg_centroid, axis = -1)
    avg_waveform_avg = np.mean(avg_waveform, axis = -1)
    avg_waveform_per_tp_avg = np.mean(avg_waveform_per_tp, axis = -1)

def create_unit_legend():
    """Create visual legend for unit colors"""
    global unit_legend_plot
    
    fig = Figure(figsize=(6, 0.5), dpi=100)
    fig.patch.set_facecolor('#33393b')
    
    ax = fig.add_subplot(111)
    ax.set_facecolor('#33393b')
    
    # Draw legend lines and text
    ax.plot([0, 0.15], [0.5, 0.5], 'g-', lw=3)
    ax.text(0.18, 0.5, 'Unit A', color='white', fontsize=12, va='center')
    
    ax.plot([0.6, 0.75], [0.5, 0.5], 'b-', lw=3)
    ax.text(0.78, 0.5, 'Unit B', color='white', fontsize=12, va='center')
    
    ax.set_xlim(0, 1.2)
    ax.set_ylim(0, 1)
    ax.axis('off')
    
    unit_legend_plot = FigureCanvasTkAgg(fig, master=root)
    unit_legend_plot.draw()
    unit_legend_plot = unit_legend_plot.get_tk_widget()

def create_hist_legend():
    """Create visual legend for histogram colors"""
    global hist_legend_plot
    
    fig = Figure(figsize=(8, 0.5), dpi=100)
    fig.patch.set_facecolor('#33393b')
    
    ax = fig.add_subplot(111)
    ax.set_facecolor('#33393b')
    
    # Draw legend elements
    ax.plot([0, 0.1], [0.5, 0.5], 'orange', lw=3)
    ax.text(0.12, 0.5, 'All scores', color='white', fontsize=10, va='center')
    
    ax.plot([0.35, 0.45], [0.5, 0.5], 'magenta', lw=3)
    ax.text(0.47, 0.5, 'Expected matches', color='white', fontsize=10, va='center')
    
    ax.plot([0.75, 0.85], [0.5, 0.5], 'white', lw=2, linestyle='--')
    ax.text(0.87, 0.5, 'Current pair', color='white', fontsize=10, va='center')
    
    ax.set_xlim(0, 1.2)
    ax.set_ylim(0, 1)
    ax.axis('off')
    
    hist_legend_plot = FigureCanvasTkAgg(fig, master=root)
    hist_legend_plot.draw()
    hist_legend_plot = hist_legend_plot.get_tk_widget()

def compute_acg(spike_times, bin_size=0.001, max_lag=0.05):
    """
    Compute autocorrelogram for a unit
    
    Parameters
    ----------
    spike_times : array
        Spike times in seconds
    bin_size : float
        Bin size in seconds (default 1ms)
    max_lag : float
        Maximum lag in seconds (default 50ms)
        
    Returns
    -------
    autocorr : array
        Autocorrelogram counts
    bin_centers : array
        Bin centers in seconds
    """
    if len(spike_times) < 2:
        return np.array([]), np.array([])
    
    # Create bins centered around 0
    bins = np.arange(-max_lag, max_lag + bin_size, bin_size)
    bin_centers = bins[:-1] + bin_size/2
    
    # Calculate autocorrelogram
    autocorr = np.zeros(len(bin_centers))
    
    # For efficiency, subsample spikes if there are too many
    if len(spike_times) > 10000:
        indices = np.random.choice(len(spike_times), 10000, replace=False)
        spike_subset = spike_times[indices]
    else:
        spike_subset = spike_times
    
    # Calculate cross-correlation with itself  
    for i, spike_time in enumerate(spike_subset[::10]):  # Subsample further for speed
        # Find spikes within max_lag of this spike
        time_diffs = spike_times - spike_time
        valid_diffs = time_diffs[(np.abs(time_diffs) <= max_lag) & (time_diffs != 0)]
        
        if len(valid_diffs) > 0:
            hist, _ = np.histogram(valid_diffs, bins=bins)
            autocorr += hist
    
    # Convert to firing rate (spikes/sec)
    if len(spike_subset) > 0:
        recording_duration = np.max(spike_times) - np.min(spike_times) if len(spike_times) > 0 else 1
        autocorr = autocorr / (len(spike_subset) * bin_size) if recording_duration > 0 else autocorr
    
    return autocorr, bin_centers

def plot_acgs(unit_a, unit_b):
    """Plot autocorrelograms for both units overlaid in single plot"""
    global acg_plot
    global clus_info
    global acg_cache
    
    # Destroy existing plot
    if 'acg_plot' in globals() and acg_plot.winfo_exists():
        acg_plot.destroy()
    
    
    # Create figure for single overlaid ACG plot
    fig = Figure(figsize=(3, 3), dpi=100)
    fig.patch.set_facecolor('#33393b')
    
    # Create single subplot for overlaid ACGs
    ax = fig.add_subplot(1, 1, 1)
    ax.set_facecolor('#2d2d2d')
    
    max_rate = 0  # Track max rate for y-axis scaling
    
    # Try to get ACGs from cache first, otherwise compute them
    try:
        # Get ACG for Unit A
        if acg_cache is not None and unit_a in acg_cache:
            # Use cached ACG
            autocorr_a, bin_centers_a = acg_cache[unit_a]
        else:
            # Compute ACG on the fly
            spike_times_a = get_spike_times_for_unit(unit_a)
            if len(spike_times_a) > 1:
                autocorr_a, bin_centers_a = compute_acg(spike_times_a)
            else:
                autocorr_a, bin_centers_a = np.array([]), np.array([])
        
        # Plot ACG for Unit A
        if len(autocorr_a) > 0:
            positive_mask = bin_centers_a >= 0
            positive_centers = bin_centers_a[positive_mask] * 1000  # Convert to ms
            positive_autocorr = autocorr_a[positive_mask]
            
            ax.plot(positive_centers, positive_autocorr, color='green', 
                   linewidth=2, alpha=0.8, label=f'Unit A ({unit_a})')
            max_rate = max(max_rate, np.max(positive_autocorr) if len(positive_autocorr) > 0 else 0)
        
        # Get ACG for Unit B
        if acg_cache is not None and unit_b in acg_cache:
            # Use cached ACG
            autocorr_b, bin_centers_b = acg_cache[unit_b]
        else:
            # Compute ACG on the fly
            spike_times_b = get_spike_times_for_unit(unit_b)
            if len(spike_times_b) > 1:
                autocorr_b, bin_centers_b = compute_acg(spike_times_b)
            else:
                autocorr_b, bin_centers_b = np.array([]), np.array([])
        
        # Plot ACG for Unit B
        if len(autocorr_b) > 0:
            positive_mask = bin_centers_b >= 0
            positive_centers = bin_centers_b[positive_mask] * 1000  # Convert to ms
            positive_autocorr = autocorr_b[positive_mask]
            
            ax.plot(positive_centers, positive_autocorr, color='blue', 
                   linewidth=2, alpha=0.8, label=f'Unit B ({unit_b})')
            max_rate = max(max_rate, np.max(positive_autocorr) if len(positive_autocorr) > 0 else 0)
        
        # Set labels and title
        ax.set_xlabel('Time lag (ms)', fontsize=10, color='white')
        ax.set_ylabel('Rate (Hz)', fontsize=10, color='white')
        ax.set_title('Autocorrelograms', fontsize=12, color='white')
        ax.set_xlim(0, 50)  # 50ms
        if max_rate > 0:
            ax.set_ylim(0, max_rate * 1.1)
        
        # Add legend
        legend = ax.legend(fontsize=8, loc='upper right')
        legend.get_frame().set_facecolor('#33393b')
        legend.get_frame().set_alpha(0.8)
        for text in legend.get_texts():
            text.set_color('white')
        
        # Style the plot
        ax.tick_params(colors='white', labelsize=8)
        ax.spines['bottom'].set_color('white')
        ax.spines['left'].set_color('white')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
    except Exception as e:
        # If spike times not available, show placeholder
        ax = fig.add_subplot(1, 1, 1)
        ax.set_facecolor('#2d2d2d')
        ax.text(0.5, 0.5, 'ACG data\nnot available', ha='center', va='center', 
                color='white', fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
    
    # Create canvas and position to the right of waveform plot
    acg_canvas = FigureCanvasTkAgg(fig, master=root)
    acg_canvas.draw()
    acg_plot = acg_canvas.get_tk_widget()
    acg_plot.grid(row=3, column=3, rowspan=2, padx=5, pady=5, sticky='nsew')

def get_spike_times_for_unit(unit_id):
    """
    Get spike times for a specific unit from Kilosort data
    
    Parameters
    ----------
    unit_id : int
        UnitMatch unit ID
        
    Returns
    -------
    spike_times : array
        Spike times in seconds for the unit
    """
    try:
        # Get session information for this unit
        global clus_info, param
        session_id = clus_info['session_id'][unit_id]
        original_id = clus_info['original_ids'][unit_id]
        
        # Get the Kilosort directory for this session
        if 'KS_dirs' in param:
            ks_dir = param['KS_dirs'][session_id]
        else:
            # Try to infer from existing paths
            return np.array([])
        
        # Load Kilosort spike times and cluster assignments
        spike_times_path = os.path.join(ks_dir, 'spike_times.npy')
        spike_clusters_path = os.path.join(ks_dir, 'spike_clusters.npy')
        
        if os.path.exists(spike_times_path) and os.path.exists(spike_clusters_path):
            # Load spike times (in samples) and cluster assignments
            spike_times_samples = np.load(spike_times_path).flatten()
            spike_clusters = np.load(spike_clusters_path).flatten()
            
            # Get sampling rate - try to load from params.py or use default
            try:
                params_path = os.path.join(ks_dir, 'params.py')
                if os.path.exists(params_path):
                    # Parse params.py for sample_rate
                    with open(params_path, 'r') as f:
                        params_content = f.read()
                    # Extract sample_rate
                    for line in params_content.split('\n'):
                        if 'sample_rate' in line and '=' in line:
                            sample_rate = float(line.split('=')[1].strip())
                            break
                    else:
                        sample_rate = 30000.0  # Default Neuropixels sample rate
                else:
                    sample_rate = 30000.0  # Default
            except:
                sample_rate = 30000.0  # Default
            
            # Filter spikes for this unit
            unit_mask = spike_clusters == int(original_id)
            unit_spike_times = spike_times_samples[unit_mask]
            
            # Convert to seconds
            unit_spike_times_sec = unit_spike_times / sample_rate
            
            return unit_spike_times_sec
        else:
            print(f"Could not find spike_times.npy or spike_clusters.npy in {ks_dir}")
            return np.array([])
            
    except Exception as e:
        print(f"Error loading spike times for unit {unit_id}: {e}")
        return np.array([])

def update(event):
    """
    Updates the GUI.
    """
        
    unit_a = int(entry_a.get())
    unit_b = int(entry_b.get())

    CV = get_cv_option() # CV = 'Avg', if AVG is selcted else it equal [0,1] or [1,0]
    CV_option = CV_tkinter.get() - 1

    table = get_table_data(unit_a, unit_b, CV)
    MakeTable(table)

    scores_table = get_unit_score_table(unit_a, unit_b, CV_option)
    make_unit_score_table(scores_table)

    plot_avg_waveforms(unit_a, unit_b, CV)
    plot_trajectories(unit_a, unit_b, CV)

    if toggle_raw_val.get() is False:
        plot_raw_waveforms(unit_a, unit_b, CV)
    else:
        if raw_waveform_plot.winfo_exists() == 1:
            raw_waveform_plot.destroy()

    add_probability_label(unit_a, unit_b, CV_option)
    add_original_ID(unit_a, unit_b)

    if toggle_UM_score_val.get() is False:
        #plot histograms based of off the CV
        if CV_option == -1:
            plot_histograms(hist_names_avg, hist_avg, hist_matches_avg, scores_to_include_avg, unit_a, unit_b)
        if CV_option == 0:
            plot_histograms(hist_names_12, hist_12, hist_matches_12, scores_to_include_GUI[CV_option], unit_a, unit_b)
        if CV_option == 1:
            plot_histograms(hist_names_21, hist_21, hist_matches_21, scores_to_include_GUI[CV_option], unit_a, unit_b)
    else:
        if hist_plot.winfo_exists() == 1:
            hist_plot.destroy()
    
    # Plot ACGs for both units if not hidden
    if toggle_acg_val.get() is False:
        plot_acgs(unit_a, unit_b)
    else:
        if 'acg_plot' in globals() and acg_plot.winfo_exists() == 1:
            acg_plot.destroy()
    

def up_options_b_list(event):
    """  
    moves up the option B list, chose the unit with th next highest probabilty of being a match with unit A.
    """
    global entry_frame
    global option_b
    global entry_b

    tmp_entry_b = int(entry_b.get())

    tmp_list = np.asarray(option_b)
    current_idx = int(np.argwhere(tmp_list == tmp_entry_b))
    if current_idx == 0:
        return
    else:
        current_idx -=1
        if entry_b.winfo_exists() == 1:
                entry_b.destroy()

    new_entry_b = tmp_list[current_idx]

    entry_b = ttk.Combobox(entry_frame, values = option_b, width = 10 )
    entry_b.set(new_entry_b )
    entry_b.bind('<<ComboboxSelected>>', update)  

    entry_b.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)

    update(event)


def down_options_b_list(event):
    """  
    moves down the option B list, chose the unit with th next highest probabilty of being a match with unit A.
    """
    global entry_frame 
    global option_b
    global entry_b

    tmp_entry_b = int(entry_b.get())

    tmp_list = np.asarray(option_b)
    current_idx = int(np.argwhere(tmp_list == tmp_entry_b))
    if current_idx == (len(tmp_list) -1):
        return
    else:
        current_idx +=1
        if entry_b.winfo_exists() == 1:
            entry_b.destroy()

    new_entry_b = tmp_list[current_idx]

    entry_b = ttk.Combobox(entry_frame, values = option_b, width = 10 )
    entry_b.set(new_entry_b)
    entry_b.bind('<<ComboboxSelected>>', update)  

    entry_b.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)

    update(event)

#sort CV function as to be calleed as part of update
def get_cv_option():
    """
    Will read in the values of the radio button and assign an appropriate value to CV.
    In general this is a list where itis [Unit A cv, UnitB cv], however it could be the string 'Avg' 
    """
    global CV_tkinter
    ChosenOption = CV_tkinter.get()
    if ChosenOption == 0:
        CV = 'Avg'
    elif ChosenOption == 1:
        CV = [0, 1]
    elif ChosenOption == 2:
        CV = [1, 0] 
    return CV

#These are the function used to selct units, including how the CV selctition radio buttons, session selction, unitselection andmoving left and rigthfor next units
def update_unit_cv():
    """
    When updating the CV we need to do /not do the following:
    - keep the same selected units,
    - update the options in boxes for unitA and unitB as matches and likely matches units can change
    - make it so when scrolling the list #MATCHIDX AUTOMATICALLY UPDATES TO THE CORrECT POINT IN THE NEW CV
    - update the screen to show the new CV
    """
    global entry_a
    global entry_b
    global option_a
    global session_entry_a
    global session_entry_b
    global match_idx
    
    #selecting the unit
    session_a = int(session_entry_a.get())
    session_b = int(session_entry_b.get()) 

    #Keep track of the unit it was before as we dont want to change theunit viewed when changing the CV
    entry_a_tmp = int(entry_a.get())
    entry_b_tmp = int(entry_b.get())

    if entry_a.winfo_exists() == 1:
        entry_a.destroy()
    if entry_b.winfo_exists() == 1:
        entry_b.destroy()

    CV = get_cv_option()
    CV_option = CV_tkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV values
        tmp_idx_a = np.argwhere((session_switch[session_a -1] <= matches_avg[:,0]) * (matches_avg[:,0] < session_switch[session_a]) == True ).squeeze()
        tmp_idx_b = np.argwhere((session_switch[session_b -1] <= matches_avg[:,1]) * (matches_avg[:,1] < session_switch[session_b]) == True).squeeze()
        in_both = np.isin(tmp_idx_a, tmp_idx_b)
        option_a = matches_avg[tmp_idx_a[in_both],:].tolist()

        entry_a = ttk.Combobox(entry_frame, values = option_a, width = 10)
        entry_a.set(entry_a_tmp)

        option_b = np.flip(np.argsort(output_avg[int(entry_a.get()),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()

    else :

        tmp_idx_a = np.argwhere((session_switch[session_a -1] <= matches_GUI[CV_option][:,0]) * ( matches_GUI[CV_option][:,0] < session_switch[session_a]) == True).squeeze()
        tmp_idx_b = np.argwhere((session_switch[session_b -1] <= matches_GUI[CV_option][:,1]) * ( matches_GUI[CV_option][:,1] < session_switch[session_b]) == True).squeeze()
        in_both = np.isin(tmp_idx_a, tmp_idx_b)

        #So is orderby unit A
        option_a =  matches_GUI[CV_option][tmp_idx_a[in_both],:].tolist()
        if CV_option == 0:
            option_a = sorted(option_a)

        entry_a = ttk.Combobox(entry_frame, values = option_a, width = 10)
        entry_a.set(entry_a_tmp)

        ##NOT SORTED FOR CV
        option_b = np.flip(np.argsort(output_GUI[CV_option][int(entry_a.get()),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()

    ##################
    entry_b = ttk.Combobox(entry_frame, values = option_b, width = 10 )
    entry_b.set(entry_b_tmp)
    
    entry_a.bind('<<ComboboxSelected>>', update_units)
    entry_b.bind('<<ComboboxSelected>>', update)
    entry_a.grid(row = 2, column = 1, columnspan = 2, stick = 'WE', padx = 5)
    entry_b.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)

    tmp_list = [int(entry_a_tmp), int(entry_b_tmp)]
    if tmp_list in option_a:
        match_idx = option_a.index(tmp_list)

    update(None)


def update_unit_entryA(event):
    """
    This will be called when the user changes the sesion A, and the following will happen:
    - Keep the same B unit
    - Set A to the first match of the new session 
    - Set the Match IDX at 0, for the new list of match pairs
    - update the options in list A
    - update the options for list B
    - update what is in the screen
    """
    global entry_a
    global entry_b
    global option_a
    global match_idx
    global session_entry_b
    
    session_a = int(session_entry_a.get())
    session_b = int(session_entry_b.get()) 
    EntryBtmp = int(entry_b.get())

    if entry_a.winfo_exists() == 1:
        entry_a.destroy()
    if entry_b.winfo_exists() == 1:
        entry_b.destroy()

    CV = get_cv_option()
    cv_option = CV_tkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV values
        tmp_idx_a = np.argwhere((session_switch[session_a -1] <= matches_avg[:,0]) * (matches_avg[:,0] < session_switch[session_a]) == True ).squeeze()
        tmp_idx_b = np.argwhere((session_switch[session_b -1] <= matches_avg[:,1]) * (matches_avg[:,1] < session_switch[session_b]) == True).squeeze()
        in_both = np.isin(tmp_idx_a, tmp_idx_b)
        option_a = matches_avg[tmp_idx_a[in_both],:].tolist()

        entry_a = ttk.Combobox(entry_frame, values = option_a, width = 10)
        entry_a.set(option_a[0][0])

        option_b = np.flip(np.argsort(output_avg[int(entry_a.get()),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()

    else :

        tmp_idx_a = np.argwhere((session_switch[session_a -1] <= matches_GUI[cv_option][:,0]) * ( matches_GUI[cv_option][:,0] < session_switch[session_a]) == True ).squeeze()
        tmp_idx_b = np.argwhere((session_switch[session_b -1] <= matches_GUI[cv_option][:,1]) * ( matches_GUI[cv_option][:,1] < session_switch[session_b]) == True).squeeze()
        in_both = np.isin(tmp_idx_a, tmp_idx_b)
        option_a =  matches_GUI[cv_option][tmp_idx_a[in_both],:].tolist()
        if cv_option == 0:
            option_a = sorted(option_a)

        entry_a = ttk.Combobox(entry_frame, values = option_a, width = 10)
        entry_a.set(option_a[0][0])

        option_b = np.flip(np.argsort(output_GUI[cv_option][int(entry_a.get()),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()
    ##################
    entry_b = ttk.Combobox(entry_frame, values = option_b, width = 10 )
    entry_b.set(EntryBtmp)
    entry_a.bind('<<ComboboxSelected>>', update_units)
    entry_b.bind('<<ComboboxSelected>>', update)

    entry_a.grid(row = 2, column = 1, columnspan = 2, stick = 'WE', padx = 5)
    entry_b.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)
    match_idx = 0

    update(event)

def update_unit_entryB(event):
    """
    This will be called when the user changes the sesion B, and the following will happen:
    - Keep the same A unit
    - update the options in list A
    - update the options for list B
    - Set B to the highest match of unit A for the new session B 
    - Set the Match IDX at 0, for the new list of match pairs
    - update what is in the screen
    """
    global entry_a
    global entry_b
    global option_a
    global match_idx
    global session_entry_a
    
    session_a = int(session_entry_a.get())
    session_b = int(session_entry_b.get()) 
    entry_a_tmp = int(entry_a.get())

    if entry_a.winfo_exists() == 1:
        entry_a.destroy()
    if entry_b.winfo_exists() == 1:
        entry_b.destroy()

    CV = get_cv_option()
    CV_option = CV_tkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV values
        tmp_idx_a = np.argwhere((session_switch[session_a -1] <= matches_avg[:,0]) * (matches_avg[:,0] < session_switch[session_a]) == True ).squeeze()
        tmp_idx_b = np.argwhere((session_switch[session_b -1] <= matches_avg[:,1]) * (matches_avg[:,1] < session_switch[session_b]) == True).squeeze()
        in_both = np.isin(tmp_idx_a, tmp_idx_b)
        option_a = matches_avg[tmp_idx_a[in_both],:].tolist()

        entry_a = ttk.Combobox(entry_frame, values = option_a, width = 10)
        entry_a.set(entry_a_tmp)

        option_b = np.flip(np.argsort(output_avg[int(entry_a.get()),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()

    else :

        tmp_idx_a = np.argwhere((session_switch[session_a -1] <= matches_GUI[CV_option][:,0]) * ( matches_GUI[CV_option][:,0] < session_switch[session_a]) == True ).squeeze()
        tmp_idx_b = np.argwhere((session_switch[session_b -1] <= matches_GUI[CV_option][:,1]) * ( matches_GUI[CV_option][:,1] < session_switch[session_b]) == True).squeeze()
        in_both = np.isin(tmp_idx_a, tmp_idx_b)
        option_a =  matches_GUI[CV_option][tmp_idx_a[in_both],:].tolist()
        if CV_option == 0:
            option_a = sorted(option_a)

        entry_a = ttk.Combobox(entry_frame, values = option_a, width = 10)
        entry_a.set(entry_a_tmp)

        option_b = np.flip(np.argsort(output_GUI[CV_option][int(entry_a_tmp),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()
    
    ##################
    entry_b = ttk.Combobox(entry_frame, values = option_b, width = 10 )
    entry_b.set(option_b[0])
    entry_a.bind('<<ComboboxSelected>>', update_units)
    entry_b.bind('<<ComboboxSelected>>', update)  

    entry_a.grid(row = 2, column = 1, columnspan = 2, stick = 'WE', padx = 5)
    entry_b.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)
    match_idx = 0    

    update(event)

def update_units(event):
    """ 
    This function is called when a option in the Unit A dropdown is selected (is a pair of units)
    - split the list of units into two pairs for unit labels for UnitA and UnitB
    - Update the unit B dropdown box to reflect matches for the new Unit A
    - Find the new Match Idx for the new pair which is selected
    """
    global match_idx
    global entry_a
    global entry_b
    global option_a
    global session_entry_b
    

    session_b = int(session_entry_b.get())
    selected = entry_a.get()

    if entry_b.winfo_exists() == 1:
        entry_b.destroy()
    # selected is s string which is xxx yyy , where xxx and yyy are the two units seperates by a space
    tmpA = selected.split()[0]
    tmpB = selected.split()[1]
    entry_a.set(tmpA)
    

    #need to make it so the list for B updates
    CV = get_cv_option()
    CV_option = CV_tkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV value
        option_b = np.flip(np.argsort(output_avg[int(tmpA),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()
    else :
        option_b = np.flip(np.argsort(output_GUI[CV_option][int(tmpA),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()
    
    ##################
    entry_b = ttk.Combobox(entry_frame, values = option_b, width = 10 )
    entry_b.set(tmpB)
    entry_b.bind('<<ComboboxSelected>>', update)  

    entry_b.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)


    #need to get list index of the selected postion
    tmp_list = [int(tmpA), int(tmpB)]
    match_idx = option_a.index(tmp_list)
    update(event)

def next_pair(event):
    """ 
    This function is called when moving to the next unit in the mathc list:
    - get the new unit idx for unit A and unit B
    - update the Units
    - update the dropdown box for unit b to reflect the new Unit A
    """
    global match_idx
    global entry_a
    global entry_b
    global option_a
    global session_entry_b

    session_b = int(session_entry_b.get())
    match_idx +=1
    tmp_a, tmp_b = option_a[match_idx]
    entry_a.set(tmp_a)

    if entry_b.winfo_exists() == 1:
        entry_b.destroy()

    #Update B, and change it's dropbox
    CV = get_cv_option()
    cv_option = CV_tkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV value
        option_b = np.flip(np.argsort(output_avg[int(tmp_a),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()
    else :
        option_b = np.flip(np.argsort(output_GUI[cv_option][int(tmp_a),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()
    
    ##################
    entry_b = ttk.Combobox(entry_frame, values = option_b, width = 10 )
    entry_b.set(tmp_b)
    entry_b.bind('<<ComboboxSelected>>', update)  

    entry_b.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)

    update(event)

def previous_pair(event):
    """ 
    This function is called when moving to the previous unit in the mathc list:
    - get the new unit idx for unit A and unit B
    - update the Units
    - update the dropdown box for unit b to reflect the new Unit A
    """
    global match_idx
    global entry_a
    global entry_b
    global option_a
    global session_entry_b
    
    session_b = int(session_entry_b.get())
    match_idx -=1
    tmp_a, tmp_b = option_a[match_idx]
    entry_a.set(tmp_a)


    if entry_b.winfo_exists() == 1:
        entry_b.destroy()

    #Update B, and change it's dropbox
    CV = get_cv_option()
    CV_option = CV_tkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV value
        option_b = np.flip(np.argsort(output_avg[int(tmp_a),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()
    else :
        option_b = np.flip(np.argsort(output_GUI[CV_option][int(tmp_a),session_switch[session_b-1]:session_switch[session_b]]) + session_switch[session_b-1])
        option_b = option_b.tolist()
    
    ##################
    entry_b = ttk.Combobox(entry_frame, values = option_b, width = 10 )
    entry_b.set(tmp_b)
    entry_b.bind('<<ComboboxSelected>>', update)  

    entry_b.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)

    update(event)

def swap_units():
    global entry_frame
    global entry_a
    global entry_b
    global session_entry_b
    global session_entry_b
    global match_idx
    global option_a

    #get all initial info
    entry_a_tmp = int(entry_a.get())
    entry_b_tmp = int(entry_b.get())
    session_a_tmp = int(session_entry_a.get())
    session_b_tmp = int(session_entry_b.get())

    #delte old box a and b
    if entry_a.winfo_exists() == 1:
        entry_a.destroy()
    if entry_b.winfo_exists() == 1:
        entry_b.destroy()

    #swap allpairs
    session_entry_a.set(session_b_tmp)
    session_entry_b.set(session_a_tmp)
    
    #Update A and B dropsown boxs
    CV = get_cv_option()
    cv_option = CV_tkinter.get() - 1 #want 0 for cv (1,2) , want 1 for cv (2,1)
    if CV == 'Avg':
        #Use average CV values
        tmp_idx_a = np.argwhere((session_switch[session_b_tmp -1] <= matches_avg[:,0]) * (matches_avg[:,0] < session_switch[session_b_tmp]) == True ).squeeze()
        tmp_idx_b = np.argwhere((session_switch[session_a_tmp -1] <= matches_avg[:,1]) * (matches_avg[:,1] < session_switch[session_a_tmp]) == True).squeeze()
        in_both = np.isin(tmp_idx_a, tmp_idx_b)
        option_a = matches_avg[tmp_idx_a[in_both],:].tolist()

        entry_a = ttk.Combobox(entry_frame, values = option_a, width = 10)
        entry_a.set(entry_b_tmp)

        option_b = np.flip(np.argsort(output_avg[int(entry_a.get()),session_switch[session_a_tmp-1]:session_switch[session_a_tmp]]) + session_switch[session_a_tmp-1])
        option_b = option_b.tolist()

    else :

        tmp_idx_a = np.argwhere((session_switch[session_b_tmp -1] <= matches_GUI[cv_option][:,0]) * ( matches_GUI[cv_option][:,0] < session_switch[session_b_tmp]) == True).squeeze()
        tmp_idx_b = np.argwhere((session_switch[session_a_tmp -1] <= matches_GUI[cv_option][:,1]) * ( matches_GUI[cv_option][:,1] < session_switch[session_a_tmp]) == True).squeeze()
        in_both = np.isin(tmp_idx_a, tmp_idx_b)

        #So is orderby unit A
        option_a =  matches_GUI[cv_option][tmp_idx_a[in_both],:].tolist()
        if cv_option == 0:
            option_a = sorted(option_a)

        entry_a = ttk.Combobox(entry_frame, values = option_a, width = 10)
        entry_a.set(entry_b_tmp)

        ##NOT SORTED FOR CV
        option_b = np.flip(np.argsort(output_GUI[cv_option][int(entry_a.get()),session_switch[session_a_tmp-1]:session_switch[session_a_tmp]]) + session_switch[session_a_tmp-1])
        option_b = option_b.tolist()

    ##################
    entry_b = ttk.Combobox(entry_frame, values = option_b, width = 10 )
    entry_b.set(entry_a_tmp)
    
    entry_a.bind('<<ComboboxSelected>>', update_units)
    entry_b.bind('<<ComboboxSelected>>', update)
    entry_a.grid(row = 2, column = 1, columnspan = 2, stick = 'WE', padx = 5)
    entry_b.grid(row = 2, column = 3, columnspan = 2, sticky = 'WE', padx = 5)

    tmp_list = [int(entry_b_tmp), int(entry_a_tmp)]

    match_idx = option_a.index(tmp_list)
    update(None)

def get_score_histograms(scores_to_include, output_threshold):
    """  
    Scores2Include is the dictionary of all scores.
    ProbThreshold is a nUnits*nUnits array where each index is 0 (Not Match) 1 (Match)
    """

    #are lsit of length 6 each list item is 2 np arrays(bins, values)
    hist_names = []
    hist = []
    hist_matches = []
    for key, values in scores_to_include.items():
        hist_names.append(key)
        hist.append(np.histogram(values, bins = 100, density = True))
        hist_matches.append(np.histogram(values[output_threshold.astype(bool)], bins = 100, density = True))

    return hist_names, hist, hist_matches    

def add_original_ID(UnitA, UnitB):
    global original_id_label
    global root
    global clus_info

    if original_id_label.winfo_exists():
        original_id_label.destroy()

    try:
        original_id_a = int(clus_info["original_ids"][UnitA].squeeze())
        original_id_b = int(clus_info["original_ids"][UnitB].squeeze())
        original_id_label = ttk.Label(root, text=f'The Original Unit IDs are:\nUnit A: {original_id_a}   Unit B: {original_id_b}', borderwidth=2, relief='groove')
    except IndexError as e:
        print(f"Error: {e}")
        original_id_label = ttk.Label(root, text='Error: Unit ID out of bounds', borderwidth=2, relief='groove')

    original_id_label.grid(row=2, column=2, ipadx=5, ipady=5)

def add_probability_label(UnitA, UnitB, CVoption):
    global bayes_label

    if bayes_label.winfo_exists():
        bayes_label.destroy()

    if CVoption == -1:
        bayes_label = ttk.Label(root, text = f'The UM probabilty for this match is:\n {np.round(output_avg[UnitA, UnitB],5)}', borderwidth = 2 , relief= 'groove' )
        bayes_label.grid(row = 2, column = 1, ipadx = 5, ipady = 5)
    else:
        bayes_label = ttk.Label(root, text = f'The UM probabilty for this match is:\n {np.round(output_GUI[CVoption][UnitA, UnitB],5)}', borderwidth = 2 , relief= 'groove')
        bayes_label.grid(row = 2, column = 1, ipadx = 5, ipady = 5)


def set_match(event = None):
    global is_match
    unit_a = int(entry_a.get())
    unit_b = int(entry_b.get())

    is_match.append( [unit_a, unit_b] )
    is_match.append( [unit_b, unit_a] )

    

def set_not_match(event = None):
    global not_match
    unit_a = int(entry_a.get())
    unit_b = int(entry_b.get())

    not_match.append( [unit_a, unit_b] )
    not_match.append( [unit_b, unit_a] )


def MakeTable(table):
    global frame_table
    
    if frame_table.winfo_exists() == 1:
        frame_table.destroy()

    total_rows = len(table)
    total_columns = len(table[0])

    colors = ['black', 'green', 'blue']
    frame_table = ttk.LabelFrame(root, text = 'UnitData')
    for i in range(total_rows):
        for j in range(total_columns):
                
            e = ttk.Entry(frame_table, width=20)
            e.insert(END, table[i][j])
            e.configure(state='readonly')             
            e.grid(row=i, column=j)
    frame_table.grid(row = 4, column = 0, padx = 10, pady = 10)

#get table data #ADD STABILTY - prob of unit with itself accros cv
def get_table_data(UnitA, UnitB, CV):
    template = [['Unit','A','B'],
        ['Avg Centroid','tmp','tmp'],
        ['Amplitude','tmp','tmp'],
        ['Spatial Decay','tmp','tmp'],
        ['Units Matches','tmp','tmp'],
        ['Stability','tmp','tmp']]

    total_rows = len(template)
    total_columns = len(template[0])
    unit_idx_tmp =[0 , UnitA, UnitB]
    table = template

    if CV == 'Avg':
        for i in range(total_rows):
            for j in range(total_columns):
                if j == 0 :
                    continue        
                if i == 0:
                    table[i][j] = str(unit_idx_tmp[j])
                if i == 1:
                    table[i][j] = str( np.round(avg_centroid_avg[:,unit_idx_tmp[j]], 2))
                if  i == 2:
                    table[i][j] = str( np.round(amplitude_avg[unit_idx_tmp[j]], 2))
                if i == 3:
                    table[i][j] = str( np.round(spatial_decay_avg[unit_idx_tmp[j]], 3))
                if i == 4:
                    table[i][j] = str(np.argwhere( (within_session*output_avg)[unit_idx_tmp[j],:] > match_threshold)).replace('[', '').replace(']', '')
                if i ==5:
                    table[i][j] = str( np.round(output_GUI[0][unit_idx_tmp[j], unit_idx_tmp[j]], 3))
           
    else:
        for i in range(total_rows):
            for j in range(total_columns):
                if j == 0 :
                    continue        
                if i == 0:
                    table[i][j] = str(unit_idx_tmp[j])
                if i == 1:
                    table[i][j] = str( np.round(avg_centroid[:,unit_idx_tmp[j],CV[j-1]], 2))
                if  i == 2:
                    table[i][j] = str( np.round(amplitude[unit_idx_tmp[j],CV[j-1]], 2))
                if i == 3:
                    table[i][j] = str( np.round(spatial_decay[unit_idx_tmp[j],CV[j-1]], 3))
                if i == 4:
                    table[i][j] = str(np.argwhere( (within_session*output_GUI[CV[0]])[unit_idx_tmp[j],:] > match_threshold )).replace('[', '').replace(']', '')
                if i ==5:
                    table[i][j] = str( np.round(output_GUI[0][unit_idx_tmp[j], unit_idx_tmp[j]], 3))
                    
    return table

def get_unit_score_table(UnitA, UnitB, CVoption):
  
    table = [['tmp'] * 2 for i in range((len(scores_to_include_avg) +1))]

    table[0] = ['Score', f'{UnitA} and {UnitB}']

    if CVoption == -1:
        for i in range(len(scores_to_include_avg)):
            for j in range(2):
                if j ==0:
                    table[i+1][j] = list(scores_to_include_avg.keys())[i]
                else:
                    table[i+1][j] = str( np.round(scores_to_include_avg[list(scores_to_include_avg.keys())[i]][UnitA, UnitB], 3))
   
    else:
        for i in range(len(scores_to_include_avg)):
            for j in range(2):
                if j ==0:
                    table[i+1][j] = list(scores_to_include_GUI[CVoption].keys())[i]
                else:
                    table[i+1][j] = str( np.round(scores_to_include_GUI[CVoption][list(scores_to_include_GUI[CVoption].keys())[i]][UnitA, UnitB], 3))

    return table

def make_unit_score_table(table):
    global score_table
    
    if score_table.winfo_exists() == 1:
        score_table.destroy()

    total_rows = len(table)
    total_columns = len(table[0])

    colors = ['black', 'Purple']
    score_table = ttk.LabelFrame(root, text = 'UM Scores')
    for i in range(total_rows):
        for j in range(total_columns):
                
            e = ttk.Entry(score_table, width = 30)
            e.insert(END, table[i][j])
            e.configure(state='readonly')             
            e.grid(row=i, column=j)

    score_table.grid(row = 5, column = 0, columnspan=2, padx = 10, pady = 10)


def plot_avg_waveforms(UnitA, UnitB, CV):
    global avg_waveform_plot
    if avg_waveform_plot.winfo_exists() == 1:
        avg_waveform_plot.destroy()

    fig = Figure(figsize = (3, 3), 
                 dpi = 100) 
    fig.patch.set_facecolor('#33393b')
    
    plt1 = fig.add_subplot(111)
    #plt1.spines[["left", "bottom"]].set_position(("data", 0))
    plt1.spines[["bottom"]].set_position(("data", 0))
    plt1.spines[["top", "right"]].set_visible(False)
    plt1.patch.set_facecolor('#2d2d2d')
    plt1.xaxis.set_label_coords(0.9,0)


    if CV =='Avg':
        plt1.plot(avg_waveform_avg[:,UnitA], 'g', label=f'Unit A ({UnitA})')
        plt1.plot(avg_waveform_avg[:,UnitB], 'b', label=f'Unit B ({UnitB})')
        plt1.set_xlabel('Time (ms)')
        plt1.set_ylabel('Amplitude (µV)')
        # plt1.set_xlim(left = 0)
        # plt1.set_xticks([])

    else:
        plt1.plot(avg_waveform[:,UnitA,CV[0]], 'g', label=f'Unit A ({UnitA})')
        plt1.plot(avg_waveform[:,UnitB,CV[1]], 'b', label=f'Unit B ({UnitB})')
        plt1.set_xlabel('Time (ms)')
        plt1.set_ylabel('Amplitude (µV)')
        # plt1.set_xlim(left = 0)

    avg_waveform_plot = FigureCanvasTkAgg(fig, master = root)
    avg_waveform_plot.draw()
    avg_waveform_plot = avg_waveform_plot.get_tk_widget()
    avg_waveform_plot.grid(row = 3, column = 0)

def plot_trajectories(UnitA, UnitB, CV):
    global trajectory_plot
    if trajectory_plot.winfo_exists() == 1:
        trajectory_plot.destroy()

    fig = Figure(figsize = (3,3), 
                 dpi = 100) 
    fig.patch.set_facecolor('#33393b')

    
    plt2 = fig.add_subplot(111)
    plt2.patch.set_facecolor('#2d2d2d')
    plt2.set_aspect(0.5)
    plt2.spines[['right', 'top']].set_visible(False)
    

    if CV =='Avg':
        
        # AM not doing a time averaged WaveIDX (where you fins goodtimepoints), will just uses CV 0 for both
        plt2.plot(avg_waveform_per_tp_avg[1,UnitA,wave_idx[UnitA,:,0].astype(bool)], avg_waveform_per_tp_avg[2,UnitA,wave_idx[UnitA,:,0].astype(bool)], 'g', label=f'Unit A ({UnitA})')
        plt2.scatter(avg_centroid_avg[1,UnitA], avg_centroid_avg[2,UnitA], c = 'g')

        plt2.plot(avg_waveform_per_tp_avg[1,UnitB,wave_idx[UnitB,:,0].astype(bool)], avg_waveform_per_tp_avg[2,UnitB,wave_idx[UnitB,:,0].astype(bool)],'b', label=f'Unit B ({UnitB})')
        plt2.scatter(avg_centroid_avg[1,UnitB], avg_centroid_avg[2,UnitB], c = 'b')

        plt2.set_xlabel(r'X position ($\mu$m)')
        plt2.set_ylabel(r'Y position ($\mu$m)')

    else:
        plt2.plot(avg_waveform_per_tp[1,UnitA,wave_idx[UnitA,:,CV[0]].astype(bool), CV[0]], avg_waveform_per_tp[2,UnitA,wave_idx[UnitA,:,CV[0]].astype(bool), CV[0]], 'g', label=f'Unit A ({UnitA})')
        plt2.scatter(avg_centroid[1,UnitA,CV[0]], avg_centroid[2,UnitA,CV[0]], c = 'g')

        plt2.plot(avg_waveform_per_tp[1,UnitB,wave_idx[UnitB,:,CV[1]].astype(bool), CV[1]], avg_waveform_per_tp[2,UnitB,wave_idx[UnitB,:,CV[1]].astype(bool), CV[1]],'b', label=f'Unit B ({UnitB})')
        plt2.scatter(avg_centroid[1,UnitB,CV[1]], avg_centroid[2,UnitB,CV[1]], c = 'b')

        plt2.set_xlabel(r'X position ($\mu$m)')
        plt2.set_ylabel(r'Y position ($\mu$m)')

    trajectory_plot=FigureCanvasTkAgg(fig, master = root)
    trajectory_plot.draw()
    trajectory_plot = trajectory_plot.get_tk_widget()
#    TrajectoryPlot.configure(bg = '#33393b')
    trajectory_plot.grid(row = 3, column = 1, columnspan = 2)

def order_good_sites(good_sites, channel_pos, n_sessions):
    # make it so it goes from biggest to smallest
    reordered_idx = np.argsort(-channel_pos[n_sessions][good_sites,2].squeeze())
    reordered_good_sites = good_sites[reordered_idx]

    #re-arange x-axis so it goes (smaller x, bigger x)
    for i in range(9):
        a,b = channel_pos[n_sessions][reordered_good_sites[[2*i, 2*i+1]],1]

        if a > b:
            #swap order
            reordered_good_sites[[2*i + 1, 2*i]] = reordered_good_sites[[2*i, 2*i+1]]
    return reordered_good_sites

def nearest_channels(max_site, max_site_mean, channel_pos, clus_info, unit, CV):

    n_sessions = clus_info['session_id'][unit]
    if CV == 'Avg':
        maxsite = max_site_mean[unit].squeeze()
        __, x, y = channel_pos[n_sessions][maxsite,:]

        good_x_sites = np.argwhere( np.logical_and((x-50 < channel_pos[n_sessions][:,1]) == True, (channel_pos[n_sessions][:,1] < x+50) == True))
        y_values = channel_pos[n_sessions][good_x_sites,2]

        y_dist_to_max_site = np.abs(y_values - channel_pos[n_sessions][maxsite,2])
        good_sites = np.argsort(y_dist_to_max_site,axis = 0 )[:18]
        good_sites = good_x_sites[good_sites]
        reordered_good_sites = order_good_sites(good_sites, channel_pos, n_sessions)

    else:
        maxsite = max_site[unit,CV]
        __, x, y = channel_pos[n_sessions][maxsite,:]

        good_x_sites = np.argwhere( np.logical_and((x-50 < channel_pos[n_sessions][:,1]) == True, (channel_pos[n_sessions][:,1] < x+50) == True))
        y_values = channel_pos[n_sessions][good_x_sites,2]

        y_dist_to_max_site = np.abs(y_values - channel_pos[n_sessions][maxsite,2])
        good_sites = np.argsort(y_dist_to_max_site,axis = 0 )[:18]
        good_sites = good_x_sites[good_sites]
        reordered_good_sites = order_good_sites(good_sites, channel_pos, n_sessions)
        
    return reordered_good_sites

def plot_raw_waveforms(unit_a, unit_b, CV):

    session_no_a = clus_info['session_id'][unit_a]
    global raw_waveform_plot
    if raw_waveform_plot.winfo_exists() == 1:
        raw_waveform_plot.destroy()

    fig = Figure(figsize=(4,6), dpi = 100)
    fig.set_tight_layout(False)
    fig.patch.set_facecolor('#33393b')

    main_ax = fig.add_axes([0.2,0.2,0.8,0.8])
    main_ax.set_facecolor('#2d2d2d')
    main_ax_offset = 0.2
    main_ax_scale = 0.8

    if CV =='Avg':
        good_channels = nearest_channels(max_site, max_site_mean, channel_pos, clus_info, unit_a, CV)

        min_x, min_y = channel_pos[session_no_a][good_channels[-2],[1,2]].squeeze()
        max_x, maxy = channel_pos[session_no_a][good_channels[1],[1,2]].squeeze()
        delta_x = (max_x - min_x) / 2
        delta_y = (maxy - min_y) / 18

        #may want to change so it find this for both units and selects the most extreme arguments
        #however i dont think tis will be necessary
        sub_min_y = np.nanmin(waveform[unit_a,:,good_channels].mean(axis=-1))
        sub_max_y = np.nanmax(waveform[unit_a,:,good_channels].mean(axis=-1))
        # shift each waveform so 0 is at the channel site, 1/9 is width of a y waveform plot
        waveform_y_offset = (np.abs(sub_max_y) / (np.abs(sub_min_y) + np.abs(sub_max_y)) ) * 1/9
    
    else:
        good_channels = nearest_channels(max_site, max_site_mean, channel_pos, clus_info, unit_a, CV[0])

        min_x, min_y = channel_pos[session_no_a][good_channels[-2],[1,2]].squeeze()
        max_x, sub_max_y = channel_pos[session_no_a][good_channels[1],[1,2]].squeeze()
        delta_x = (max_x - min_x) / 2
        delta_y = (maxy - min_y) / 18
        #may want to change so it find this for both units and selects the most extreme arguments
        #however i dont think this will be necessary
        sub_min_y = np.nanmin(waveform[unit_a,:,good_channels, CV[0]])
        sub_max_y = np.nanmax(waveform[unit_a,:,good_channels, CV[0]])
        # shift each waveform so 0 is at the channel site, 1/9 is width of a y waveform plot
        waveform_y_offset = (np.abs(sub_max_y) / (np.abs(sub_min_y) + np.abs(sub_max_y)) ) * 1/9


    #make the main scatter positiose site as scatter with opacity 
    main_ax.scatter(channel_pos[session_no_a][good_channels,1], channel_pos[session_no_a][good_channels,2], c = 'grey', alpha = 0.3)
    main_ax.set_xlim(min_x - delta_x, max_x + delta_x)
    main_ax.set_ylim(min_y - delta_y, maxy + delta_y)

    for i in range(9):
        for j in range(2):
            #may need to change this positioning if units sizes are irregular
            if j == 0:
                #The peak in the waveform is not half way, so maths says the x axis should be starting at
                #0.1 and 0.6 so the middle is at 0.25/0.76 however chosen these values so it loks better by eye

                #
                ax =  fig.add_axes([main_ax_offset + main_ax_scale*0.25, main_ax_offset + main_ax_scale*(i/9 - 1/18 + waveform_y_offset), main_ax_scale*0.25, main_ax_scale*1/9])
            if j == 1:
                ax = fig.add_axes([main_ax_offset + main_ax_scale*0.75, main_ax_offset + main_ax_scale*(i/9 - 1/18 + waveform_y_offset), main_ax_scale*0.25, main_ax_scale*1/9])

            if CV =='Avg':
                ax.plot(waveform[unit_a,:,good_channels[i*2 + j]].mean(axis=-1).squeeze(), color = 'g')             
                ax.plot(waveform[unit_b,:,good_channels[i*2 + j]].mean(axis = -1).squeeze(), color = 'b', lw=0.8)
            else:
                ax.plot(waveform[unit_a,:,good_channels[i*2 + j], CV[0]].squeeze(), color = 'g')             
                ax.plot(waveform[unit_b,:,good_channels[i*2 + j],CV[1]].squeeze(), color = 'b', lw=0.8)                
            ax.set_ylim(sub_min_y,sub_max_y)
            ax.set_axis_off()


    main_ax.spines.right.set_visible(False)
    main_ax.spines.top.set_visible(False)
    main_ax.set_xticks([min_x, max_x])
    main_ax.set_xlabel('X position ($\mu$m)', size = 14)
    main_ax.set_ylabel('Y position ($\mu$m)', size = 14)
    


    raw_waveform_plot = FigureCanvasTkAgg(fig, master = root)
    raw_waveform_plot.draw()
    raw_waveform_plot = raw_waveform_plot.get_tk_widget()
    #RawWaveformPlot.configure(bg = '#33393b')

    raw_waveform_plot.grid(row = 3, column = 4, columnspan = 2, rowspan = 4, padx = 15, pady = 25, ipadx = 15)

def plot_histograms(hist_names, hist, hist_matched, scores_to_include, unit_a, unit_b):

    global hist_plot
    if hist_plot.winfo_exists() == 1:
        hist_plot.destroy()

    fig = Figure(figsize= (4,6), dpi = 100, layout = 'constrained')
    fig.patch.set_facecolor('#33393b')
    axs = fig.subplots(3,2, sharex='col')
    axs = axs.flat


    # Create title mapping
    title_mapping = {
        'amp_score': 'Amplitude score',
        'spatial_decay_score': 'Spatial decay score', 
        'centroid_overlord_score': "C'oid overlord score",
        'centroid_dist': "C'oid distance score",
        'waveform_score': 'Waveform score',
        'trajectory_score': 'Trajectory score'
    }
    
    #loop over indexes.. 
    for i in range(len(hist)):
        axs[i].step(hist[i][1][:-1], hist[i][0], color = 'orange', label='All scores' if i == 0 else "")
        axs[i].step(hist_matched[i][1][:-1], hist_matched[i][0], color = 'magenta', label='Expected matches' if i == 0 else "")
        axs[i].set_ylim(bottom=0)
        
        # Use improved title from mapping
        plot_title = title_mapping.get(hist_names[i], hist_names[i])
        axs[i].set_title(plot_title, fontsize = 12)
        
        # Add ylabel
        axs[i].set_ylabel('% units', fontsize=10)
        
        axs[i].axvline(scores_to_include[hist_names[i]][unit_a,unit_b], ls = '--', color = 'white', label='Current match pair' if i == 0 else "")
        axs[i].set_facecolor('#2d2d2d')
    

    hist_plot = FigureCanvasTkAgg(fig, master = root)
    hist_plot.draw()
    hist_plot = hist_plot.get_tk_widget()

    hist_plot.grid(row = 3, column = 6, columnspan = 2, rowspan = 4, padx = 5, pady = 20)

