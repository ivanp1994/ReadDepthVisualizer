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
> git cd ReadDepthVisualizer
> pip install .

Installation is optional, since you can simply import `start_function` from the package.
Example:
```
from ReadDepthViz import start_function
start_function()
```
# Usage
Once the package is installed, simply import `start_function` from it and the GUI loop will start.

# Prior requirements
The package is written to take the outputs of CNVPytor files and visualize the information within them. Briefly, CNVPytor is a read-depth based SV calling software
that works by segmenting the reference genome into contiguous segments called *bins* (or *windows*). Every read of a sample is put into its corresponding bin (it's *binned*)
and total number of reads paired with the bin is then counted. This count represents *RD signal* of a particular bin. RD signals are then normalized so that the 
average RD signal is 2 which corresponds to the diploid copy number. (Addendum: In CNVPytor's CNV signal, this RD signal is normalized that the average RD signal is 1).
The outputs of CNVPytor files are HDF (*Hierarchical Data Format*) files, and every HDF file requires a BAM (*Binary Alignment Map*) file. A BAM file is created when a
sample in FASTA format (sample) is aligned against reference FASTA (reference). It's recommended to have information about reference FASTA file in the form
of gene/exon annotations and assembly reports.

## Test data
The folder *data* has a dummy test data. The HDF files are found in *pytor* folder and belong to various samples of Mexican Tetra. The information about reference genome 
HDF files were created against is found in *assembly.csv*, *genes.csv*, and *exons.csv*. The information about samples is given in the *metadata.csv*.
Briefly, there are 3 fishes from 5 different regions (total of 15 samples) of Central America. Those 15 fishes can be further subdivided into cave and surface ecotypes, and new and old lineages.


# GUI help
The GUI interface consists of one menu titled **Data**, a bar to draw signal, and a console display. The **Data menu** consists of three segments.

* First segment envelops parameters used to read in RD signals (*Set data parameters*) plot RD signals (*Set grid parameters*)
* Second segment envelops options to add data (*Add HDF files*), add metadata (*Load Metadata*), and loading CNV data (*Load CNV data*). The latter is not yet implemented meaningfully.  
* Third segment envelops options for adding information about reference genome - the information about its assembly (*Load assembly*) and about the genomic position of genes (*Load genes*) and exons (*Load exons*)

The latter two are self-explanatory.

## Setting parameters
Clicking on *Set data parameters* or *Set grid parameters* opens up a new popup window.

### Setting parameters for data 
Popup window consists of 4 fields. First field (*Compress*) relates to how is the RD signal stored - half stores the RD signal in half (numpy type float16), etc.
Last field (*GCF*) enables user to enter GCF ID of a reference genome. Then the user can load all information about it from the NCBI data base.
The second and third field relate to which RD signal is used for plotting. Briefly, RD signal depends on selected bin size. For example, a RD signal of bin size 100 will be markedly different from the RD signal
of bin size 1000. Therefore one HDF file can have multiple RD signals. Authors of CNVPytor recomment the binsize for which the ratio of mean RD signal to its standard deviation is between 4 and 5 (corresponding to a relative standard deviation of 0.25 and 0.2). This boundary can be changed, but our dummy data has only one RD signal. 

### Setting parameters for plotting
Metadata (that is, data about our data) must be loaded for this to work. 
Popup window consists of checkable boxes for *Genes* and *Exons*, two dropdown menus for *Hue* and *Grid*, and two buttons to add and remove *Grids*.

Ticking either *Genes* or *Exons* or both adds another grid on top of all recommended where genes and exons are drawn. Genes are drawn as thin red lines, while exons are drawn as thick red lines.
Dropdown menu for *Hue* enables different colors of RD signals based on selected value, and Grid adds another plot where only RD signals of a selected value are drawn.

![Example of dropdown menus](https://github.com/ivanp1994/ReadDepthVisualizer/assets/84333373/c8cfcf26-9b3a-4979-923e-1825da4cb816)

## Drawing 
Once everything is set, RD signal can be drawn. To draw RD signal, its locus must be set - the locus consists of chromosome, start position, and end position.
An example image from our dummy data is given below.

![image](https://github.com/ivanp1994/ReadDepthVisualizer/assets/84333373/e691afc7-f646-45d3-951f-b4f0ced6ff42)

In the above example, we've colored RD signal based on lineage (new and old), drawn all surface RD signals on the top ax and cave RD signals on the bottom.
The gene we get is 103043895 which codes for a protein called *retinal G protein coupled receptor b*. Part of the gene is fully deleted in cave fishes.

### Options for drawing

Additionally, RD plots are semi-interactive. 
* Mouse movement displays real time coordinates of genomic position and diploid copy number
* Right clicking on grid axes opens up a popup menu where a sample's visiblity and be toggled on and off
* Right clicking on track (gene and exon) axes opens up a popup menu where user can toggle annotations, strandness, exons, genes, and legend on/off
* Double left clicking on grid axes prints out all genes found on this positions
* Left click, drag, and release, zooms in this area
* Double right clicking or wheel click saves the image


