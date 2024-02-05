# Electric Vehicle Auralization Toolbox (EVAT)

This repository contains Matlab code and data for the auralization of electric vehicle passages, accompanying the paper _Auralization of Electric Vehicles for the Perceptual Evaluation of Acoustic Vehicle Alerting Systems_ (submitted to Acta Acustica).

## Contents
- **code**
  - `generatePassBy.m`
  - `bem2SH.m`
  - `getAvasCoeffsFromRecordings.m`
  - functions
    - `applyAirAttenuation.m`
    - `generateAvasSignal.m`
    - `generateOutSignal.m`
    - `getAvasCoeffs.m`
    - `getIRs.m`
    - `getMinPhaseFilt.m`
    - `makeBinaural.m`
    - `sh2p.m`
    - `ss2ds.m`
    - `truncateSH.m`
- **data**
  - directivity
    - bem
    - sh
  - synthesis
    - coeffs
    - settings
  - measurements
    - download from [here](https://doi.org/10.5281/zenodo.10610491)     


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
Please download the data manually and add the directories 'ambience', 'avas', 'passby', and 'tires' to the 'data/measurements' folder.
