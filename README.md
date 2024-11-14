# SPADE-PTM: SEC Protein Assembly Dynamics Evaluator - Post Translational Modifications

## Description
SPADE is a python tool built for analyzing SEC-MS (Size Exclusion Chromatography paired with Mass Spectrometry) data for evaluating of protein assembly state dynamics between multiple samples. The SPADE-PTM tool was built to analyze SEC-MX data - incorporating PTM (Post-Translational Modification) enriched SEC-MS data. SPADE-PTM has 4 basic functions: (a) TMT data preprocessing, (b) peak picking, (c) peak alignment between samples and PTMs, and (d) additional calculations (e.g., fold-changes, peak region evaluations,...).

- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Running Sample](#running-sample)
- [Sample Results](#sample-results)

## System Requirements

### OS
This package has been tested on the following systems:
+ macOS: Sonoma (14.4.1)
+ Windows 11

### CPU
This package has been tested with the following processors:
+ Apple (Intel, M2, M3 Max)
+ Windows (Intel)

### GPU
This package has been tested with the following memory capacities:
+ 16-36 GB

### Software
The package was tested with these versions
+ `python==3.12.3`
+ `pip==24.0`
+ `click==8.1.7`
+ `scipy==1.13.1`
+ `scikit-learn==1.5.0`
+ `matplotlib==3.8.4`
+ `numpy==1.26.4`
+ `pandas==2.2.2`

## Installation
Note: the following steps (installation and running) are meant to be run in terminal or command line. It is recommended to open a directory with at least 1GB of space. 

Installation Time Estimation: 10-30 minutes (depending on memory)
### Downloading the Repository

Clone the repository 

The repository can also be downloaded from the [GitHub](https://github.com/mjlab-Columbia/SPADE-PTM). Intalling git on Windows, refer [here](https://git-scm.com/downloads/win).

```
git clone https://github.com/mjlab-Columbia/SPADE-PTM.git
```

### Download Dependencies
Go to repository folder
```
cd SPADE-PTM
```
#### Conda Method
1. Download dependencies
```
conda env create -f environment.yaml
```

2. Activate the environment
```
conda activate spade_ptm
```
#### Pip Method
Download dependencies directly from requirements file 
```
pip install -r requirements.txt
```

### Download Sample Data
Download and unzip the [Data Inputs](https://sec-mx-example-data.s3.amazonaws.com/data_inputs.zip)

Note: if `wget` is not available, install using `brew install wget` (refer to [homebrew](https://brew.sh/)). Windows users may have difficulty using 'wget' and are recommended to use the [Data Inputs](https://sec-mx-example-data.s3.amazonaws.com/data_inputs.zip) download link.

```
wget https://sec-mx-example-data.s3.amazonaws.com/data_inputs.zip
unzip data_inputs.zip
```

## Running Sample
#### Prerequisites
+ Start in the repository folder `cd peak_alignment`
+ Check that conda environment is activated or dependencies are installed
+ Have sample data folder `/sample_data` unzipped and in the repository folder

#### Additional Information
+ Most inputs for the code are defaulted and set to run on the `sample_data` folder.
+ Space: expect results to take a few hundred MB in space
+ Time: running all steps should take 10-15 minutes

### Step 1. Create output/results directory
```
mkdir -p results
```

### Step 2. Preprocessing
The preprocessing script takes the TMT intensity tables from Spectomine and performed the following operations: (a) reshaping the data, (b) normalizing and averaging between TMT mixes, (c) median normalizing between mixes, (d) averages replicates ("-ar"), (e) exporting to SQL database.
```
python dataPreprocessing.py -o results/test_database.db -ar
```

### Step 3. Peak Picking
The peakpicking module performs the following operations on the SQL database output from the Preprocessing step: (a) smoothing using a [linear digital filter](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html) and exporting the smoothed data, (b) peak picking from each elution profile using `scipy.signal.find_peaks` function, (c) plotting peak mutliplicity of each sample, (d) exporting result to SQL database.
```
python peakPicking.py -s results/test_database.db
```

### Step 4. Peak Alignment 
In the first iteration of peak alignment, the pipeline first assigns PTM peptides to the global data and then finds the optimal alignment between conditions. For PTM alignment, each PTM peak is assigned to a global peak - allowing for multiple PTM peptides to match to a single global peak. Then a global Condition alignment is performed to find the optimal alignment given the Threshold. 

By default, the threshold is set to 3. 
```
python PeakMatching.py -s results/test_database.db -c 3
```

### Step 5. Post Alignment Editing 
After aligning elution peaks between samples, the following operations are performed: (a) calculation of the mean elution fraction of an aligned peak, (b) calcuation of fold changes of peaks between samples, (c) molecular weight estimations and peak region categorization.
```
python postMatchingEditing.py -s results/test_database.db
```
## Sample Results

Download and unzip the [Data Outputs](https://sec-mx-example-data.s3.amazonaws.com/data_outputs.zip) to see expected results.

```
wget https://sec-mx-example-data.s3.amazonaws.com/data_outputs.zip
unzip sample_outputs.zip
```
