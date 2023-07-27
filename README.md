# ReadDepthVisualizer
The following package consists of a GUI implemented in Tkinter package.
The GUI enables the user to integrate read-depth signal ("RD signal") from many different samples and to visualize the similarities and 
differences in RD signal.

The package is depdendent on CNVPytor which can be found [here](https://github.com/abyzovlab/CNVpytor.git) but does not require it.
The final output of CNVPytor are HDF (Hierarchical Data Format) files which store RD signal in many different versions and according to many different binsizes.
The **ReadDepthVisualizer** rips out normalized RD signal (normalized so that the RD signal of a non-deleted and non-duplicated region is 2) and stores it into 
its own file. 

# Requirements
The package is lightweight, and requires only few extra packages. These are:
* requests >= 2.27
* h5py >=3.4.0
* numpy >= 1.22.3
* pandas >= 1.4.1
* matplotlib >= 3.4.2


The recommended minimum version of Python is 3.8.12.
If you're using CNVPytor, you likely already have everything you need.

# Installation
As of now, the quickest installation is via cloning from github
> git clone https://github.com/ivanp1994/ReadDepthVisualizer.git
> 
> pip install .

Installation is optional, since you can simply import `start_function` from the package.
Example:
```
from ReadDepthViz import start_function
start_function()
```
# Usage
Once the package is installed, simply import `start_function` from it and the GUI loop will start.

# GUI help
...gonna need to add additional stuff here
