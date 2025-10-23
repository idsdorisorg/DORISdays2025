# [DORIS DAYS 2025] Jupyter Notebooks for Hands On Sessions

[![DORIS DAYS 2025](assets/logo_DD_2025.png)](https://ids-doris.org/ids/meetings/ids-meetings.html%20)

This repository is meant to hold source code for the "hands-on" sessions that will take 
place during the [DORIS DAYS 2025](https://ids-doris.org/ids/meetings/ids-meetings.html%20) in NTUA, Athens.

The repository contains:
 - a *core* python module named `dsoclasses`, on top of which
 - a list of [Jupyter Notebooks](https://jupyter.org/) are available, showcasing trivial DORIS data analysis tasks.

To use the Notebooks, you need to install the core python module.

## Installation of `dsoclasses` module

This guide walks you through installing the `dsoclasses` Python module from source.

0. Prerequisites

Make sure your system has the following installed:

- Python 3.8 or higher
- `git`, `pip`, and optionally `venv`

 1. Clone the repository

```bash
git clone https://github.com/idsdorisorg/DORISdays2025.git
cd DORISdays2025
```

2. Create a Virtual Environment (Optional)

**Note that depending on your OS and setup, you may need to replace `python` with `python3`.**

```bash
python -m venv .venv
source .venv/bin/activate
```

For **Windows** users the above should be replaced with:
```bash
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

3. Upgrade build tools

```bash
python -m pip install --upgrade pip
```

4. Build and Install the package

```bash
pip install -e .
```

### Updating

To fetch latest changes and/or additions to the online repository, you will need to 
run the following command from the root of the project (where the `pyproject.toml` is 
located): `git pull origin`. No other step should be needed.


## Jupyter Notebooks

The notebooks are placed under the `JupyterLab` folder. Hence, assuming jupyterlab 
is available on your system (if not, `pip install jupyterlab` would do it) the following 
command should launch a local web server and open JupyterLab in your browser 
`jupyter lab --notebook-dir=JupiterLab/` (from the top-level directory).

## A note on data

This repository comes with a few data files that are needed to run the examples presented in 
the notebooks. Alternate or updated data should be seeked at the dedicated web repositories.