# PETPipelineMATLAB

## Description

The `PETPipeline` pipeline is a state-of-the-art PET preprocessing pipeline (wrapper) for MATLAB. This pipeline is designed to execute multiple preprocessing steps on BIDS structured datasets that have at least one dynamic PET scan and one anatomical MRI scan. The key steps integrated within the pipeline include:

- **Motion correction**
- **Co-registration**
- **Segmentation**
- **Partial volume correction**
- **Kinetic modeling**

## Usage

To utilize the `PETPipeline` function in your MATLAB environment:

```matlab
PETPrep(data_dir, config)
```

### Parameters:

- **data_dir (str)**: This is the path directing to the BIDS data directory.
- **config (str)**: This is the path directing to the PETPrep configuration file. This config should be placed in a 'code' directory in the main BIDS directory. Example: 'config1.json'.

**Note:** The function doesn't return any values.

## Sequential Steps Enforced by the Pipeline:

1. Transition to the BIDS data directory.
2. Integrate BIDS data into MATLAB.
3. Incorporate the configuration file.
4. Initiate derivatives directories.
5. Activate `ReconAll` for FreeSurfer reconstruction.
6. Perform GTM segmentation with `GTMSeg`.
7. Convert FreeSurfer output to BIDS structure with `ConvertFS2BIDS`.
8. Perform motion correction (various methods available).
9. Co-register with `CoReg` and plot results with `PlotCoReg`.
10. Perform GTM PVC with `GTMPVC`.
11. Convert PETsurfer output to BIDS structure with `PETsurfer2TAC`.
12. Perform kinetic modeling (various models available).

## Author

Martin Norgaard, Stanford University, 2022.

## License

[License details, if any]

## Acknowledgements

[Optional section for any acknowledgements or citations]
