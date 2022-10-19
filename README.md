# liquid-liquid-phase-seperation-simulations
## Overview
Numerical simulations the phase seperation phenomena between liquids; this repository provides numerical solutions to a variety of models tuned to exhibit this phenomenon. 

## Installation and usage (Python)
### Mandatory dependencies
A [requirements.txt](requirements.txt) file is provided, containing all the dependencies necessary for the scripts to run. Installing of these dependencies can be performed painless using pip and the command below:

```bash
pip install -r requirements.txt
```

If pip fails, below is a listing of the dependencies and links to their respective repositories. Please ensure, to avoid issues, that the installed version of each dependency matches the major version and has greater or equal minor version to the listing below:

- [matplotlib v3.5.x](https://github.com/matplotlib/matplotlib)
- [numpy v1.22.x](https://github.com/numpy/numpy)
- [scipy v1.8.x](https://github.com/scipy/scipy)

### Recommended extra ([FFmpeg](https://github.com/FFmpeg/FFmpeg))

[FFmpeg](https://github.com/FFmpeg/FFmpeg) is [matplotlib's](https://github.com/matplotlib/matplotlib) preffered tool for creating video capable files, as such its usage is highly recommended. The simplest method to install FFmpeg is through your platform's local package manager (If one is available). 

#### OS X
```bash
brew install ffmpeg
```

#### Debian/Ubuntu
```bash
sudo apt install ffmpeg
```

### Manual installation 
If a package manager is not available, then binaries for [FFmpeg](https://github.com/FFmpeg/FFmpeg) can be acquired from the [official download links](https://ffmpeg.org/download.html). For [matplotlib](https://github.com/matplotlib/matplotlib) to be able to use [FFmpeg](https://github.com/FFmpeg/FFmpeg), the `$PATH` environment variable must contain the directory to the downloaded [FFmpeg](https://github.com/FFmpeg/FFmpeg) binary.

To verify that [FFmpeg](https://github.com/FFmpeg/FFmpeg) has been installed correctly, run the command:
```bash
ffmpeg
```
If an error is reported, the installion was performed incorrectly and try again. Otherwise, if you have followed all the instructions up this step, congratulations, everything is set up and ready to be run.
