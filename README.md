# STORM_Analysis_Parameter_Extraction

This MATLAB code analyzes point cloud images of H2B distribution in the nucleus to predict the spatial distribution of chromatin compaction and extact heterochromatin domains sizes. 

------------------------------------------------------------------------------
## System Requirement
# Operating system:
This package is supported for macOS and windows. The package has been tested on the following systems:
- macOS Ventura 13.4.1 Processor: M1 Max chip; RAM: 64 GB
- Windows 10 Processor: Intel Exon; RAM: 64 GB


# MATLAB dependencies
This package requires a MATLAB version of 2021a or newer. This package has been tested on the following MATLAB version:
- MATLAB_R2022a


## Installation guide
1. Install MATLAB (https://www.mathworks.com/products/matlab.html). The installation usually takes approximately 1 hour.
2. After installing MATLAB, users should make sure the “Machine Learning and Deep Learning” package is installed by click “APP” in the menu. This package is required for DBSCAN algorithm. We assume all basic package suggested by MATLAB have already been installed. Users may still get notified in the command window of MATLAB of all necessary apps while execute the code. If a package is suggested to be installed, follow the instruction in command window to install and then click “RUN” again.

## Demo
1. Open “Nuclear_STORM_Analysis_MATLAB_v3.m” and click “Run”. The demo code will take around 4 minutes to finish depends on the processor and data size.
2. New windows will then open to visualize the Voroni polygon of the nucleus.
3. All the measurements of the nucleus will be shown in the command window of MATLAB. All the post-analysis data can be found in the cell structure “data” in the workspace window.

## Sample Output
1. A sample "command window output" is shown in the folder "Output_ImgLib".
2. A sample png image of the Voronoi Density is shown in the folder "Output_ImgLib".

Note for “Nuclear_STORM_Analysis_MATLAB.m”: This file is configured to analyze the heterochromatin domains in the nucleus. The necessary custom written functions are store in the folder “FuncLib”. The H2B localization file named under “GSK_DEMO_1.txt” is stored in the folder “Input_LocsLib”.

## Instruction for use
To analyze more data, H2B localization (saved as a txt file) should be stored in the “Input_LocsLib” folder. Our package will automatically read all the txt file in this folder.
