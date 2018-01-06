Blood MRI Analysis Tool 
-----------------------

For measuring T1 and T2 in sets of images.

Installation
------------

The following are required:
    - python (version 2 or 3)
    - qtpy
    - pyqt, qt (version 4 or 5)
    - qtawesome
    - matplotlib
    - numpy
    - scipy
    - scikit-image
    - numexpr
    - pydicom
    - sympy
    - six

The easiest way to install these is using conda, as part of 
[miniconda](https://conda.io/miniconda.html) or
[anaconda](https://www.anaconda.com/download/).

Once the requirements are met, you can install with pip
```bash
pip install git+https://github.com/shportnoy/blood_roi_tool
```

or by cloning the repository, which will also download some example data
```bash
git clone https://github.com/shportnoy/blood_roi_tool
cd blood_roi_tool
python setup.py install
```

Thereafter, you can launch the tool with
```bash
blood_tools
```