# README: User guide for SSR algorithm based on DaMo

## Author

Sijie Li, Mengqi Lu, Jing Yuan

## Introduction

Single-Spectrum Reconstruction (SSR) is a rapid super-resolution reconstruction algorithm designed for Digital Array Modulation Microscopy (DaMo), allowing reconstruction from illumination modulated images using only a single high-order spectrum. For the specific principle of SSR, refer to the article “Three-dimensional super-resolution imaging with suppressed background via digital array modulation microscopy”.

## System requirements

Our software is developed on Windows 10 / Windows 11 with MATLAB (MathWorks®). MATLAB of R2021b or a later version is preferred. As a lightweight technology, SSR necessitates no specific computational or storage resources from the underlying computer hardware.

The software has been tested on MATLAB R2024b (Windows 11, version 24H2), MATLAB R2024a (Windows 10, version 22H2), and MATLAB R2021b (Windows 11, version 23H2).

## Installation guide

1. Download the repository.
2. Extract folder "DaMo_MATLAB" and execute the main function “SSRecon.m” in it.

Total installation time is under 1 second, requiring only the extraction of the 64 MB software_DaMo folder with no dependency installations, compilation, or environment configuration (MATLAB R2021b or later must be pre-installed), which is tested on an Intel Core i7-9700 desktop (3.00 GHz base frequency) under Windows 11 Version 24H2.

## Demo

Four raw images of sample tubulin are located at: “DaMo_MATLAB\Demo\Raw_data”. Files named “tubulin\_\*.\*” use the index “\*” to denote the raw modulated images under different illumination modulation intensities.

Parameters in the code have been pre-configured for the sample data. Select the desired SSR processing mode (HC/HS/HF) by setting “SSR_mode”, then execute the code.

The SR results of three SSR modes, “tubulin_HC_SR.\*”, “tubulin_HS_SR.\*”, and “tubulin_HF_SR.\*”, will be saved at: “DaMo_MATLAB\Demo\SSR_result”.

Processing of a 2048×2048-pixel tubulin image in HS mode required 1.62 seconds on a test system comprising an Intel Core i7-9700 desktop (3.00 GHz base frequency) running MATLAB R2024b under Windows 11 Version 24H2. The demo execution time includes overheads for image loading, parameter initialization, SR reconstruction, and output generation.

## User data processing instructions

1. **Input**
   Place the raw images in the "input" folder. Ensure the raw data files are named using the following convention: "filename\_x":
   - filename: must exactly match the value of the filename variable defined in the code
   - x: must be an integer from 1 to 4, where each number corresponds to a different modulation intensity
2. **Parameters instruction**
   Configure these parameters before execution:
   - SSR_mode: select processing mode (HC, HS, or HF)
   - Blockm/Blockn: raw image dimensions in pixels
   - PixelSize: physical pixel size in nanometers
   - NA: objective numerical aperture
   - lambda: emission wavelength in micrometers
3. **Run “SSRecon.m”**
4. **Output**
   The output results will be automatically saved under the path denoted by the “pathout” variable.

## Contact

If you have any questions, please contact yuanj@hust.edu.cn for help.

## License

This software is licensed under the Apache License 2.0. See LICENSE for details.
