Blood MRI Analysis Tool 
-----------------------

For measuring T1 and T2 in sets of images.

Installation
------------

The following are required:
1. python (version 2 or 3), and
2. the following packages:
* qtpy
* pyqt, qt (version 4 or 5)
* qtawesome
* matplotlib
* numpy
* scipy
* scikit-image
* numexpr
* pydicom
* sympy
* six

The easiest way to install python is by downloading [anaconda]:
(https://www.anaconda.com/download/)

Once anaconda is installed, run the following command to obtain the necessary packages (listed above):
```bash
conda install qtpy qtawesome matplotlib numpy scipy scikit-image numexpr \
    pydicom sympy six -c defaults -c conda-forge
```

Once all the requirements are met, you can install the tool (and download some example data) by cloning this repository:
```bash
git clone https://github.com/shportnoy/blood_roi_tool
cd blood_roi_tool
python setup.py install
```

Thereafter, you can launch the tool with the command:
```bash
blood_tools
```