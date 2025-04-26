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

Optional (or for development) dependencies are listed in `requirements.txt` and `environment.yml`.  

## Installation

You can install **pyoedometer** in two main ways:

### 2. From source (editable / development)

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

#### Alternative: requirements.txt

If you prefer a plain pip-based setup, you can install dependencies directly:

```bash
pip install -r requirements.txt
```

Then install the package:

```bash
pip install .
```

## Quick Verification

After installation, open a Python REPL and try:

```python
>>> import pyoedometer
>>> print(pyoedometer.__version__)
0.1.0
```

If no errors occur, you’re all set!  

---

## License

Released under the GNU GPL-3.0-or-later.  

## Author

Thomas Ingeman-Nielsen (<thin@dtu.dk>)  
