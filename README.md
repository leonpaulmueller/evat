# Electric Vehicle Auralization Toolbox (EVAT)

This repository contains Matlab code and data for the auralization of electric vehicle passages accompanying the paper:
> Leon Müller and Wolfgang Kropp, _Auralization of Electric Vehicles for the Perceptual Evaluation of Acoustic Vehicle Alerting Systems_, submitted for review and publication to Acta Acustica on 02. February 2024

## Contents
- **code**
  - `generatePassBy.m` - Main script, generates binaural EV pass-by auralization
  - `bem2SH.m` - Converts BEM pressure to a set of spherical harmonic coefficients
  - `getAvasCoeffsFromRecordings.m` - Analyzes AVAS recordings to obtain synthesis coefficients
  - functions
    - `applyAirAttenuation.m` - Applies air attenuation according to ISO 9613-1:1993 
    - `generateAvasSignal.m` - Synthesizes AVAS or Tire/Road noise source signal for a given velocity vector and coefficient set
    - `generateOutSignal.m` - Generates output signal at receiver position by applying moving Green's functions to source signal
    - `getAvasCoeffs.m` - Get AVAS synthesis coefficients from reference recording
    - `getIRs.m` - Extrapolate SH directivity to obtain Green's functions for all source positions
    - `getMinPhaseFilt.m` - Get minimum-phase representation of input IR
    - `makeBinaural.m` - Apply HRTF to input signal depending on source position to obtain binaural output.
    - `sh2p.m` - Get pressure from SH coefficients at evaluation coordinates
    - `ss2ds.m` - Converts single-sided to double-sided spectrum
    - `truncateSH.m` - Truncates SH coefficients to a lower order
- **data**
  - directivity
    - bem - BEM results for the three different vehicles (pressure in frequency domain on 131st order Lebedev sphere)
    - sh - SH Coefficients for AVAS (SH order N=64) and tire/road noise radiation (N=16)
  - synthesis
    - coeffs - AVAS and tire/road noise synthesis coefficients, obtained by analyzing reference recordings
    - settings - Pre-defined analysis settings for the three evaluated vehicles
  - measurements
    - download from [here](https://doi.org/10.5281/zenodo.10610491)
   
For consistency with the paper, we use the following aliases for the three evaluated vehicles:
- Vehicle A: Tesla Model Y 2021
- Vehicle B: Volkswagen ID.3 Pro Performance 2021
- Vehicle C: Nissan Leaf 2018


## External Dependencies
The toolbox requires the external function `getSH.m` (https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/getSH.m) to generate SH basis functions.
Additionally, the script `generatePassBy.m` uses the HRTF set `HMSII.sofa` (https://sofacoustics.org/data/database/thk/HMSII.sofa), which you can alternatively replace with your own HRTF.

## Matlab Requirements
The code was tested for Matlab versions R2021b and newer. For full functionality, the following toolboxes are required:
- Audio Toolbox
- Signal Processing Toolbox
- DSP System Toolbox
- Statistics and Machine Learning Toolbox
- Image Processing Toolbox

## Measurement Data
Due to the GitHub file size limitations, all measurement data was uploaded to Zenodo (https://doi.org/10.5281/zenodo.10610491). 
Please download the data manually and add the directories `ambience`, `avas`, `passby`, and `tires` to the `data/measurements` folder.
