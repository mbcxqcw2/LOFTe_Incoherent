# NOTES

A stripped-down version of the LOFT-e incoherent beamformer for use with the LOFT-e docker: https://github.com/mbcxqcw2/LOFTe_Docker

- Code to initialise a new filterbank file uses the read_header() function developed by Scott Ransom for PRESTO: https://github.com/scottransom/presto/blob/4959ffeac8d369d8883987287179297f93aea773/python/presto/filterbank.py
- The original incoherent beamformer and associated scripts and pipelines is stored at: /share/nas1/LOFT-e/software/incoherent_beam_pipeline
- This git repository is stored at: /share/nas1/LOFT-e/software/incoherent_beam_pipeline/git_uploads/LOFTe_Incoherent/

More detailed notes to come.

---

# DEPENDENCIES

- numpy
- matplotlib
- astropy https://github.com/astropy/astropy (and dependencies)
- sigpyproc https://github.com/ewanbarr/sigpyproc (and dependencies)
- presto https://github.com/scottransom/presto (don't need to install, just need the python scripts in your python path to use the sigproc.read_hdr_val() function)

---

# USAGE

1) To start, git clone this directory and add it to your python path. 

2) To combine two or more filterbank files, run:

```
from Incoherent_9 import CombineFils
```

