# for array manipulation
import numpy as np

# for reading fits files
from astropy.io import fits
from astropy.table import Table

# for plotting
import matplotlib.pyplot as plt
from matplotlib import gridspec

# for filtering time-tag events by time
from costools import splittag
# for processing COS data
import calcos

# for system files
import glob
import os
from pathlib import Path

# for downloading the data
from astroquery.mast import Observations

# For supressing an unhelpful warning:
import warnings
from astropy.units import UnitsWarning

test_dir = Path('./test-NEW-A/')  # Using the pathlib style of system path
plots_dir = test_dir / 'plots-test'
reg_int_dir = test_dir / 'regular_intervals'

# Make the directories in case they don't exist
test_dir.mkdir(exist_ok=True), plots_dir.mkdir(exist_ok=True)
reg_int_dir.mkdir(exist_ok=True) # Output directory for intervals specified by start/end/length

x1d=sorted(glob.glob('*_x1d.fits'))
corrtag_a=sorted(glob.glob('*corrtag_a.fits'))
corrtag_b=sorted(glob.glob('*corrtag_b.fits'))

# Build the figure structure
#fig = plt.figure(figsize=(12, 10))
# Using gridspec to let us control panel sizes and locations
#gs = fig.add_gridspec(nrows=4, ncols=3)
# Get the time data for each exposure
my_size = 480  # Number of seconds in a bin

for i,chosen_corrtag in enumerate(corrtag_a):
# Split the file every 100 seconds
#chosen_corrtag=corrtag_a[0]
#reg_int_dir = test_dir / 'regular_intervals/'chosen_corrtag[0:9]

    reg_int_str=os.path.join(test_dir,'regular_intervals', chosen_corrtag[0:9])
    reg_int_dir=Path(reg_int_str)
    reg_int_dir.mkdir(exist_ok=True) # Output directory for intervals specified by start/end/length
    hdu=fits.open(chosen_corrtag)
    #for head in list(hdu[0].header.keys()):
    #    if hdu[0].header[head]=='COMPLETE':
    #        hdu[0].header[head]=='PERFORM'
    time=int(hdu[1].data['time'][-1])
    div=int(time/my_size)
    splittag.splittag(
        infiles=chosen_corrtag,
        outroot='./test-NEW-A/regular_intervals/'+chosen_corrtag[0:9]+'/'+chosen_corrtag[0:9],
        starttime=0.,
        increment=time/div,
        endtime=time
    )


    # Gather the split sub-exposure files:
    spec_outlist = sorted(glob.glob('./test-NEW-A/regular_intervals/'+chosen_corrtag[0:9]+'/*.fits'))
    epoch_markers=[]
    # Make histogram lightcurves as in previous plots:
    binsize = 10 # binsize of lightcurve histogram in seconds
    for splitfile in spec_outlist:  # For each of our newly split up files:
        epoch_number = os.path.basename(splitfile).split('_')[1]
        # Read in the file as a table of events:
        events_table = Table.read(splitfile, 1)
        event_times = events_table['TIME']
        epoch_markers.append((min(event_times.tolist()),max(event_times.tolist())))
        hist = plt.hist(
            event_times,
            bins=int((max(event_times)-min(event_times))/binsize),
            alpha=1, label=f"sub-exposure "+str(max(event_times.tolist()))
        )

    # Add time ranges "epochs" or "windows of time" as transparent spans:
    for epoch_time, color in zip(epoch_markers, ['C0', 'C1']):
        plt.axvspan(
            epoch_time[0], epoch_time[1],
            color=color, alpha=0.1, label=None
        )




    # This cell is for ensuring you have a valid "lref" directory of reference files.
    if not os.environ.get('lref'):
        if os.path.exists('/grp/hst/cdbs/lref/'):  # If on STScI Network/VPN
            os.environ['lref'] = '/grp/hst/cdbs/lref/'
            print("It looks like you are probably on the STScI network; setting lref to '/grp/hst/cdbs/lref/'")
        else:  # If not on STScI Network/VPN
            # PLEASE EDIT THIS PATH
            os.environ['lref'] = '/Users/hrrsergio/crds_cache/references/hst/lref'
            if not os.path.exists(os.environ['lref']):  # Check if that path exists
                print("It doesn't look like that's a valid path. Deleting it.")
                del os.environ['lref']  # delete this nonexistant path
    else:
        found_lref = Path(os.environ.get('lref'))
        print(f"You already have an lref path in your environment variables -\
     It's {found_lref}\n")

    assert os.path.exists(os.environ['lref']), f"The path to your lref directory is invalid ({os.environ['lref']})"

    for split_corrtag in spec_outlist:  # When we run CalCOS on corrtags, we must go 1-by-1
        # Define epoch as which chunk of the initial exposure. epoch 2 contains the transit.
        epoch_number = os.path.basename(split_corrtag).split('_')[1]
        print(f"Extracting file {split_corrtag} using CalCOS")
        # Make a sub-directory of output/calcos/ named for each epoch:
        cal_output_dir = f'./test-NEW-A/calcos/'+chosen_corrtag[0:9]+f'/epoch{epoch_number}/'
        os.makedirs(cal_output_dir, exist_ok=True)
        # Extract the spectrum from each of the sub-exposures:
        hdu_split=fits.open(split_corrtag)

        hdu_split[0].header['DOPPCORR']='PERFORM'
        calcos.calcos(split_corrtag, outdir=cal_output_dir, verbosity=0)
    # Print a message at the end to let us know it's finished:
    print("Done running the pipeline.")

    # Find all the `x1d` files:
    processed_files = sorted(glob.glob('./test-NEW-A/calcos/'+chosen_corrtag[0:9]+'/epoch*/*x1d.fits'))

    # Set up figure
