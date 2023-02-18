# Installation

First, you need to install PyMOL v2.x on your computer. If you do not have it installed, please refer to PyMOL [website](https://pymol.org/2/).

_PyMOL2 parKVFinder Tools_ is available to use with PyMOL v2.x and has been been developed using Qt interface and Python3.

After installing PyMOL v2.x, plug-in requires the installation of [toml](https://pypi.org/project/toml/) and [pyqt5](https://pypi.org/project/PyQt5/) modules and PyMOL's native Python does not have it installed. So you need to install it:

```bash
# PyMOL v2 (python 3)
$ pip3 install toml pyqt5
# In an conda environment via Anaconda Cloud
$ conda install -c conda-forge toml
conda install -c anaconda pyqt 
```

_Note:_ [pip](https://pypi.org/project/pip/) or
[Anaconda](https://www.anaconda.com/) package management system
installation is required.

Then, to install the _PyMOL2 parKVFinder Tools_ on PyMOL v2.x, download the latest `PyMOL2-parKVFinder-Tools.zip` from [here](https://github.com/LBC-LNBio/parKVFinder/releases/latest/download/PyMOL2-parKVFinder-Tools.zip) and follow these steps:

1. Open PyMOL.
2. Go to **Plugin** menu, click on **Plugin Manager**.
3. The **Plugin Manager** window will open, go to the **Install New Plugin** tab.
4. Under **Install from local file** panel, click on **Choose file...**.
5. The **Install Plugin** window will open, select the _PyMOL2-parKVFinder-Tools_.zip` that you downloaded earliar.
6. The **Select plugin directory** window will open, select
   **/home/\<user\>/.pymol/startup** and click **OK**.
7. The **Confirm** window will open, click **OK**.
8. The **Sucess** window will appear, confirming that the plug-in has
   been installed.
9. Restart PyMOL.
10. _PyMOL2 parKVFinder Tools_ is ready to use.

Or, instead of selecting `PyMOL2-parKVFinder-Tools.zip` (Step 5), you can select `__init__.py` file on `parKVFinder/tools/PyMOL2-parKVFinder-Tools` directory tree.
