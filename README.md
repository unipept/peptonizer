# Peptonizer2000 - Integration into Unipept

## Installation
In order to run the Peptonizer2000 on your own system, you should install Conda, Mamba and all of its dependencies.
Follow the installation instructions step-by-step for an explanation of what you should do.

* Make sure that Conda and Mamba are installed. If these are not yet present on your system, you can follow the instructions on their [README](https://github.com/conda-forge/miniforge).
* Go to the "workflow" directory by executing `cd workflow` from the terminal.
* Run `conda env create -f env.yml` (make sure to run this command from the workflow directory) in order to install all dependencies and create a new conda environment (which is named "peptonizer" by default).
* Run `mamba install -c conda-forge -c bioconda -n peptonizer snakemake` to install snakemake which is required to run the whole workflow from start-to-finish.
* Run `conda activate peptonizer` to switch the current Conda environment to the peptonizer environment you created earlier.
* Start the peptonizer with the command `snakemake --use-conda --cores 1`. If you have sufficient CPU and memory power available to your system, you can increase the amount of cores in order to speed up the workflow.

## Chain of operations
In order to compute the probability that a specific taxon is present in a metaproteomics sample, the Peptonizer2000 goes through a whole chain of different operations.
This section highlights each of the different steps that are performed and how these are interconnected to each other.

## Generating Peptonizer Python package
Generating and packaging the Peptonizer as a Python package can be done by executing this command in the terminal (from within the `peptonizer/peptonizer` directory): `python setup.py bdist_wheel`.

### Input
The following data is required by the Peptonizer as input and are then used during the analysis pipeline to transform the input into the required output probabilities:
* A file containing "rescored" PSMS, use MS2Rescore to transform your input to something that can be fed into this pipeline.
* 

