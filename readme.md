# pyoedometer

**pyoedometer** is a Python package for interpreting and plotting oedometer (consolidation) tests.  

## Requirements

- **Python**: 3.12 or newer  
- **Core dependencies** (see `pyproject.toml`):  
  - `numpy`  
  - `matplotlib`  
  - `pandas`  
  - `openpyxl`  
  - `tables`  
  - `PyYAML`  
  - `PyQt5`  
  - `jinja2`
  - `scipy`  

## Installation

You can install **pyoedometer** in two main ways:

### 1. From source (editable / development)

1. Clone this repository:  
   ```bash
   git clone https://github.com/tingeman/pyoedometer.git
   cd pyoedometer
   ```
2. Install package:  
   ```bash
   pip install .
   ```

To install the package locally but allow you to modify the code in place, use>

2. Install in “editable” mode:  
   ```bash
   pip install -e .
   ```


### 2. Using conda & environment file

A ready-to-go conda environment is provided in `environment.yml` (named `pyoedometer_3p12`):

```bash
conda env create -f environment.yml
conda activate pyoedometer_3p12
```

This will install Python 3.12 plus all declared dependencies via conda and pip.

---


## Usage

To interpret oedometer (consolidation) test data using the sqrt-time method (as per the ISO standard), follow these steps:

1. **Setup folder structure and organize datafiles**
   The code expects the following folder structure:

   ```text
   experiment_sample01/
   ├── sample_data/
   |   ├── hobo_data
   |   |   └── hobo_data.hdf
   │   ├── cr6_lvdt_table.dat
   │   ├── cr6_pt100_table.dat
   |   ├── hobo_data.hdf
   │   └── config_sample.yml
   └── run_sample.py
   ```

   The location of data files can be adjusted in the `config_sample.yml` file.

   In some cases, hobo data consists of multiple files. The code should read and process all hobo csv files in the dedicated directory, then save the combined time series to a `hobo_data.hdf` file (which is subsequently used in processing, if it exists).

2. **Prepare your configuration:**  
   Edit the `sample_config.yml` file to specify your test data files, sample properties, and analysis parameters. This YAML file controls all input to the module.

3. **Run the analysis:**  
   Use the command line to execute the `run_sample.py` script. For example:
   ```bash
   python run_sample.py
   ```
   The script will read the `./sample_data/sample_config.yml` file and process the data accordingly. 

   It will do basic processing for all steps (including establishing eps0 and epsf for each step), then do detailed processing of steps specified in the `config_sample.yml` file.
   
   Detailed processing may include:
   - plotting overview of strain and temperature data for entire expreiment
   - plotting overview of individual load or temperature steps
   - interpreting time curves of individual load steps
   - plotting time curves annotated with interpretation
   - saving individual step data to excel files,
   - saving interpretation to excel file
   - plotting the consolidation curve

4. **Review the output:**  
   Results (such as interpreted parameters and plots) will be saved to the locations defined in your config file.

For more details on the configuration options, see the comments in `sample_config.yml` and the code.

### General workflow

In general, the processing workflow consists of iteratively updating the `sample_config.yml` file and running the `run_sample.py` script to produce output. Processing steps are changed to work on certain elements, until an acceptable result is obtained.

Often it is most convenient to work with a single time step at a time, adjusting settings so that only that particular step is processed and plotted.

1. **Define the experiment setup**
   Modify the `sample_config.yml` configuration section, to provide appropriate paths to files, and define lvdt and temperature measurements.
   ***NB: Do not reuse calibration factors from previous experiments, unless you know what sensors they refer to or otherwise know what you are doing.***

2. **Define the experiment history**
   Setup `sample_config.yml` history section with approximate date and times for each load step start and end.

3. **Adjust for exact step start and end times**
   Use the step overview plots (use zoom functionality) to find exact start and end times and update the `sample_config.yml` file accordingly. 
   
   The step start time should be set so that the first point encountered at or after the start time is lower than `eps0`

4. **Define time range for interpretation of `epsf`**
   Inspect the step overview plots and specify the time interval for calculating `epsf`. Typically we use the final 10 minutes of the time series for each step (`dt=10`). Setting is specified for each step independently.
   
   This setting should be defined for all time curves (steps) before moving on to detailed interpretation of load steps.

5. **Interpret time curves**
   Conduct interpretation of individual time curves. This is ofte most conveniently handled one time step at a time. 

   `timec` `t1` and `t2` are used to define two points to use for fitting a line to the primary consolidation. The code will take the nearest points in the time series to `t1` and `t2` and construct a straight line (in sqrt-time) between these two points.
   `timec` `t3` and `t4` are used to define two points to use for defining the slope of the secondary consolidation (creep).

   The step start time must be correctly set for the interpretation step to work properly. If the start time is off, the sqrt-time fitting will be obscured.

   `temp`  `t0` and `t1` are used to specify the start and end time for calculating the average temperature of the time step.

   ***For load steps***, the average temperature of the entire time step should be reoprted.

   ***For temperature steps*** (no change in load), the final sample temperature should be reported, meaning that `t0` and `t1` should be set to reflect the final temperature of the sample.

   All times are given in minutes since start of the current time step.

6. **Plot consolidation curve**
   The script will generate the file `./sample_data/interpretation.xlsx` containing all interpretation parameters. 
   
   When all steps have been successfully interpreted/processed,
   rename the final `interpretation.xlsx` file to e.g. `interpretation_final.xlsx`.

   Update the `sample_config.yml` `paths` `interpretation_file` setting to point to this file.

   Open the excel file, navigate to the `plotting` sheet, and modify annotation settings in this sheet to adjust how the consolidation curve is plotted (e.g. which steps to annotate, and how to off set labels on the plot).

   When working with the final consolidation plot, processing of individual steps can be turned off, as all necessary data is stored in and read from the `interpretation_final.xlsx` file.


---

## License

Released under the GNU GPL-3.0-or-later.  

## Author

Thomas Ingeman-Nielsen (<thin@dtu.dk>)  
