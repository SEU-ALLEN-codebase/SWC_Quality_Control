# SWC_Quality_Control
A tool for batch quality control of swc files

## Project Structure
- **analyze.h/.cpp:** Contains functions that analyze the structure of swc to detect errors when the reconstruction is complete.
-  **colldetection.h/.cpp:** Contains functions that detect errors in swc.
-  **detecttask.h/.cpp:** Contains the detection thread.
-  **sort_swc.h/.cpp:** Contains functions that sort swc.
-  **others:** Contains functions that support the effective execution of quality control.

## Run
SwcQualityControl --input_swc_dir_path
