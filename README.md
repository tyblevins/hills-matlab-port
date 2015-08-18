# Matlab port for Hills Algorithm
## Crowdsourcing Seizure Detection Algorithms Using Kaggle and ieeg.org
### University of Pennsylvania

This repository contains programs to interface the Michael Hills algorithm with IEEG.org. Annotated datasets are downloaded from IEEG.org using the portal toolbox in matlab, the user chooses how many seizures to train on, the algorithm is trained, and then the user either chooses a custom testing segment or all of the remaining seizures are used as testing segments. Detections are uploaded to IEEG.org.

## Methods
### Hills_Detector_py2
- Wraps the python code in Matlab to preserve exact algorithm.

**Note**:
File paths in Hills_Detector_py2.m and the .sh files will need to be altered to fit your path

### Hills_Detector2
- Pure matlab port.
- Not an exact match to the python algorithm.

## Dependencies
- matlab 2014b
- pw file for ieeg toolbox (need an ieeg account)
- python 2.7
- scikit_learn-0.14.1
- numpy-1.8.1
- pandas-0.14.0
- scipy
- hickle

