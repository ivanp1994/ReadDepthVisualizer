 # -*- coding: utf-8 -*-
"""
Main wrapper
"""

import os
import logging
import tkinter as tk
from tkinter import filedialog
from tkinter.scrolledtext import ScrolledText
from functools import partial

from threading import Thread
import queue

import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

#import reading as rd
from .reading import get_samples,read_all_rdchroms,write_read_depth_file,write_cnv_call_file
from .gcf_fetcher import read_assembly_report, read_feature_table, read_exons
from .plotting import ReadDepthDrawer as rdd
from .plotting import chr_len_form
from .widgets import LociEntry, OptionMenu
from .windows import DataParameterWindow, GridParameterWindow, PopupMenu



DEFAULT_GRID = {"genes":None,
                "exon":None,
                "hue":None,
                "grids":None}


def select_file():
    """function for opening file"""
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(filetypes=(
        ("CSV files", "*.CSV"),("CSV files", "*.csv"),
        ("TSV files", "*.tsv"), ("TSV files", "*.TSV"),
        ))

    sep = None
    if file_path:
        if file_path[-3:].lower()=="csv":
            sep=","
        elif file_path[-3:].lower()=="tsv":
            sep="\t"
        else:
            file_path=None
    return file_path,sep

def chunkify_cnvs(row,size):
    """
    Breaks up bigger CNVs into smaller segments
    For example, a CNV:
        >>> [0 1000]
    will be broken up into
        >>> [0,500], [500,1000]
    """
    chrom = row.chrom
    start = row.start
    end = row.end
    
    ran = [x for x in range(start,end,size)]
    starts = ran[:-1]
    ends = ran[1:]
    chrom = ["NC_064425.1"]*len(starts)
    return list(zip(chrom,starts,ends,))

def convert_cnvs_to_list(cnv_df,lower_limit,upper_limit):
    "convers CNV dataframe to a list of [chrom,start,end]"
    #filter out small CNVS
    cnv_df = cnv_df[cnv_df["size"]>lower_limit] #87747 -> 87403
    #filter out big CNVS
    cnv_l = cnv_df[cnv_df["size"]<upper_limit][["chrom","start","end"]].values.tolist()
    #break up big CNVs
    cnv_l2 = cnv_df[cnv_df["size"]>upper_limit].apply(chunkify_cnvs,axis=1,size=upper_limit).tolist()
    return cnv_l + cnv_l2


class QueueHandler(logging.Handler):
    """Class to send logging records to a queue

    It can be used from different threads
    The ConsoleUi class polls this queue to display records in a ScrolledText widget
    """
    # Example from Moshe Kaplan: https://gist.github.com/moshekaplan/c425f861de7bbf28ef06
    # (https://stackoverflow.com/questions/13318742/python-logging-to-tkinter-text-widget) is not thread safe!
    # See https://stackoverflow.com/questions/43909849/tkinter-python-crashes-on-new-thread-trying-to-log-on-main-thread

    def __init__(self, log_queue):
        super().__init__()
        self.log_queue = log_queue

    def emit(self, record):
        self.log_queue.put(record)

class ConsoleUi:
    """Poll messages from a logging queue and display them in a scrolled text widget"""

    def __init__(self, frame,logger):
        self.frame = frame
        # Create a ScrolledText wdiget
        self.scrolled_text = ScrolledText(frame, state='disabled', height=12)
        self.scrolled_text.grid(row=0, column=0, sticky=("nsew"))
        self.scrolled_text.configure(font='TkFixedFont')
        self.scrolled_text.tag_config('INFO', foreground='black')
        self.scrolled_text.tag_config('DEBUG', foreground='gray')
        self.scrolled_text.tag_config('WARNING', foreground='orange')
        self.scrolled_text.tag_config('ERROR', foreground='red')
        self.scrolled_text.tag_config('CRITICAL', foreground='red', underline=1)
        # Create a logging handler using a queue
        self.log_queue = queue.Queue()
        self.queue_handler = QueueHandler(self.log_queue)
        formatter = logging.Formatter('%(message)s')
        self.queue_handler.setFormatter(formatter)
        logger.addHandler(self.queue_handler)
        # Start polling messages from the queue
        self.frame.after(100, self.poll_log_queue)

    def display(self, record):
        msg = self.queue_handler.format(record)
        self.scrolled_text.configure(state='normal')
        self.scrolled_text.insert(tk.END, msg + '\n', record.levelname)
        self.scrolled_text.configure(state='disabled')
        # Autoscroll to the bottom
        self.scrolled_text.yview(tk.END)

    def poll_log_queue(self):
        # Check every 100ms if there is a new message in the queue to display
        while True:
            try:
                record = self.log_queue.get(block=False)
            except queue.Empty:
                break
            else:
                self.display(record)
        self.frame.after(100, self.poll_log_queue)

# MAIN INTERFACE
class MainInterface():
    """Wrap around tkinter Tk, holding datasets, etc"""
    def quit(self):
        self.master.quit()
        self.master.destroy()
        #sys.exit(0)
        os._exit(0)

    def __init__(self):
        self.master = tk.Tk()
        self.command_menu = None

        self.data_paths = None
        self.file_path = "data.rd"
        self.cnv = None #not necessary as of now
        self.cnvindex = None #for iterating over cnv list
        self.status = None #necessary, but make it not
        self.gcf = None #optional parameter
        self.assembly = None # DataFrame
        self.exons = None # DataFrame
        self.genes = None # DataFrame
        self.plotter = None #plotting object
        self.canvas = None #canvas to draw in
        self.message = None #message displaying mouse coords
        self.popup_axes = None #popup menu to choose samples
        self.popup_track = None #popup menu to choose track
        self.canvas_frame = tk.Frame(master=self.master,width=800,height=800)
        self.coord_display = None #displays coordinates
        self.logger = logging.getLogger()

        self.parameters = {"data":None,
                           "cnv":None,
                           "rd":None,
                           "region":None
                           }

        self.commands = {"chrom":None,
                         "start":None,
                         "end":None}


        #self.buttons = dict()

        self.master.title("Read-Depth Visualizer")
        self.master.protocol("WM_DELETE_WINDOW", self.quit)


        #CONFIGURE MENU
        #and add commands
        menubar = tk.Menu(self.master)
        self.master.config(menu=menubar)
        self.data_menu = tk.Menu(menubar,tearoff=False)
        menubar.add_cascade(label="Data",menu=self.data_menu)
        self.data_menu.add_command(label="Set data parameters",command=self.set_data_params)
        self.data_menu.add_command(label="Set grid parameters",command=self.set_grid_params)

        self.data_menu.add_separator()
        for label,command in zip(["Add HDF files","Load CNV data","Load Metadata"],
                                 [self.add_data,self.add_cnv,self.add_metadata]):
            self.data_menu.add_command(label=label,command=command)
        self.data_menu.add_separator()

        commands = [(x,partial(self.gene_fromfile,what=x),partial(self.gene_fromgcf,what=x))
                    for x in ["assembly","genes","exons"]]
        for _g, _fromfile,_fromgcf in commands:
            assembly_menu = tk.Menu(self.master,tearoff=False)
            assembly_menu.add_command(label="Load from file",command=_fromfile)
            assembly_menu.add_command(label="Load from GCF",command=_fromgcf)
            self.data_menu.add_cascade(label=f"Load {_g}",menu=assembly_menu)

        #CONFIGURE LOGFRAME
        #scrollbar configuration
        log_frame = tk.Frame(master=self.master,height=30)
        #console = ConsoleUi(log_frame,self.logger)
        ConsoleUi(log_frame,self.logger)
        self.logger.setLevel(logging.INFO)

        #CONFIGURE COMMAND FRAME
        command_frame = tk.Frame(master=self.master)
        _initoptions = read_all_rdchroms(self.file_path)
        if _initoptions is None:
            self.logger.info("No chromosomes to draw - load data")
            _initoptions = ["None"]

        chrom_option = OptionMenu(command_frame,"Chrom",_initoptions)
        start_option = LociEntry(command_frame,"Start")
        end_option = LociEntry(command_frame,"End")

        chrom_option.place(0,0,sticky=None)
        start_option.place(0,2,None)
        end_option.place(0,4,None)
        start_option.entry.bind("<Return>",self.plot_region)
        end_option.entry.bind("<Return>",self.plot_region)
        self.commands = {"chrom":chrom_option,
                         "start":start_option,
                         "end":end_option}
        tk.Button(master=command_frame, command=self.plot_region, text = "Draw").grid(row=0,column=6,sticky="nse")
        tk.Button(master=command_frame,command=self.get_cnvcoords, text = "CNV").grid(row=0,column=7,sticky="nse")

        #PLACE FRAMES
        command_frame.grid(row=1,column=0,sticky="w")
        log_frame.grid(row=2, column=0,sticky="nw")
        self.canvas_frame.grid(row=2,column=1,rowspan=2)
        #self.master.grid_rowconfigure(index, kw)

# ADDING GENOME INFO
    def gene_fromfile(self,what=""):
        """read the selected from file"""
        file_path,sep = select_file()
        if file_path is None:
            return None
        data = pd.read_csv(file_path,sep=sep)
        if what=="assembly":
            self.assembly = data
        if what=="genes":
            self.genes = data
        if what=="exons":
            self.exons = data
        self.logger.warning(f"{what.upper()} successfuly set")

    def gene_fromgcf(self,what=""):
        """read the selected from gcf"""
        #check for valid gcf
        if self.gcf is None:
            self.logger.warning("Set a valid GCF accession ID to use this option")
            return None
        self.logger.warning(f"Searching NCBI for {what}...")

        if what=="assembly":
            self.assembly = read_assembly_report(self.gcf)
        if what=="genes":
            self.genes = read_feature_table(self.gcf)
        if what=="exons":
            self.exons = read_exons(self.gcf)
        self.logger.warning(f"{what.upper()} successfuly set")

    def determine_sample_col(self):
        """determines sample column from status and samples"""
        s = get_samples(self.file_path)
        status = self.status
        if s is not None and status is not None:
            di = {k:len(s.intersection(set(status[k]))) for k in status.columns}
            col = max(di,key=di.get)
            self.status = status.set_index(col)

# ADDING DATA
    def _add_data(self,**params):
        "target for a different thread"
        write_read_depth_file(**params)
        self.determine_sample_col()
        self.logger.info("Loaded Data")
        self.commands["chrom"].repopulate_options(read_all_rdchroms(self.file_path))
        self.logger.info("Added new options for chromosomes")
        #get samples
        samples = get_samples(self.file_path)
        self.logger.error(f"Added a total of {len(samples)} : {samples}")
        self.initialize_plotter()

    def add_data(self):
        "Adds data"
        root = tk.Tk()
        root.withdraw()
        self.data_paths = list(filedialog.askopenfilenames())

        params = self.parameters["data"]
        if params is None:
            params = dict()
        params.update({"pytor_files":self.data_paths,"output_file":self.file_path,"assembly_report":self.assembly})

        if self.data_paths:
            #t = Thread(target=rd.write_read_depth_file,kwargs=params)
            #rd.write_read_depth_file(**params) #add parameters
            #t.start()
            Thread(target=self._add_data,kwargs=params).start()

    def add_metadata(self):
        "Adds metadata - data with information about samples"
        file_path,sep = select_file()
        if file_path:
            self.status = pd.read_csv(file_path,sep=sep)
            self.determine_sample_col()
        self.logger.error("Added metadata")
    
    def _add_cnv(self):
        "Adds CNVs but don't know what to do with this just yet"
        params = self.parameters["data"]
        if params is None:
            params = dict()
        params.update({"pytor_files":self.data_paths,"output":None,"ret":True})

        cnvdf = write_cnv_call_file(**params)
        self.cnv = convert_cnvs_to_list(cnvdf,1000,100000)
        self.cnvindex = 0
        self.logger.info("Loaded CNVs")
    
    def add_cnv(self):
        """wrapper around _cnv for thread safe logging"""
        if self.cnv is None and self.data_paths is not None:
            Thread(target=self._add_cnv).start()
        if self.cnv is None:
            self.logger.error("Cannot add CNV data without HDF files - Add HDF files to use this functionality")
            return
            #log cannot get cnv because no data blah
        
    def get_cnvcoords(self):
        "gets CNV coordinates of the next CNV and draws it"
        if self.cnv is None:
            self.logger.error("No CNV data - load CNV data to use this functionality")
            return
        chrom,start,end = self.cnv[self.cnvindex]
        #print(chrom,start,end)
        self.cnvindex = self.cnvindex + 1
        if self.cnvindex == len(self.cnv):
            self.cnvindex = 0
        #set commands and plot region
        chrom_com = self.commands["chrom"]
        start_com = self.commands["start"]
        end_com = self.commands["end"]
        
        chrom_com.variable.set(chrom)
        for value,com in zip([start,end],[start_com,end_com]):
            com.entry.delete(0,tk.END)
            com.entry.insert(0,f"{value}")
        
        self.plot_region()
    
# SETTING PARAMETERS
    def set_data_params(self):
        "enters parameters"
        w = DataParameterWindow(self)
        w.master.protocol("WM_DELETE_WINDOW",lambda w=w:self.param_data_exit(w))

    def param_data_exit(self,w):
        "exists parameter window and loads data in"
        self.gcf = w.data_variables.pop("gcf",None).get()
        self.parameters["data"] = {k:v.get() for k,v in w.data_variables.items()}
        w.master.destroy()

    def set_grid_params(self):
        "set grid parameters"
        if self.status is None:
            self.logger.warning("Load metadata to set grid parameters")
            return None
        w = GridParameterWindow(self.status)
        w.master.protocol("WM_DELETE_WINDOW",lambda w=w:self.param_grid_exit(w))

    def param_grid_exit(self,w):
        "exits parameter window and loads grid data in"
        gridp = w.extract_parameters()
        w.master.destroy()

        genes = gridp["genes"]
        exon = gridp["exon"]
        #making this explict and as clear as possible
        if genes:
            genes = self.genes
        else:
            genes = None

        if exon:
            exon = self.exons
        else:
            exon = None

        gridp["genes"] = genes
        gridp["exon"] = exon
        self.parameters["rd"] = gridp
        self.logger.info("Loaded parameters")
        self.initialize_plotter()

# PLOTTER
    def initialize_plotter(self):
        """initializes plotter"""
        samples = get_samples(self.file_path)
        if samples is None:
            self.logger.warning("Cannot initialize plotter without data")
            return None
        rdparams = self.parameters["rd"]

        if rdparams is None:
            self.logger.warning("No meta-data, drawing with default settings")
            rdparams = DEFAULT_GRID
            rdparams["grids"] = [list(samples)]
            rdparams["hue"] = {"all":list(samples)}

        self.plotter = rdd(rdparams,samples,self.file_path)
        self.logger.info("Initialized plotter")
        fig = self.plotter.figure
        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
        if self.coord_display is not None:
            self.coord_display.destroy()

        self.canvas = FigureCanvasTkAgg(fig, master=self.canvas_frame)
        #the bottom numbers are magical to expand figure to canvas
        #TODO - test these numbers more closely
        #eh - close enough
        self.canvas.figure.subplots_adjust(left=0.042, right=0.99, top=0.988, bottom=0.071)
        self.canvas.get_tk_widget().grid(row=0, column=0)
        #connect to function
        self.canvas.mpl_connect("button_press_event",self.canvas_key_press)
        self.canvas.mpl_connect("motion_notify_event",self.mouse_motion)
        self.message = tk.StringVar()
        self.coord_display = tk.Label(master=self.canvas_frame,textvariable=self.message)
        self.coord_display.grid(row=1,column=0)

    def plot_region(self,*args):
        "plots a region"""
        
        #initialize plotter if it's not initialized
        #clear ax if plotter is initialized
        if self.plotter is not None:
            self.plotter.clear_ax()
        else:
            self.initialize_plotter()
            self.canvas.get_tk_widget().grid(row=0, column=0)

        chrom = self.commands["chrom"].variable.get()
        start_entry = self.commands["start"]
        end_entry = self.commands["end"]

        start, start_display = start_entry.get_loci()
        end, end_display = end_entry.get_loci()

        if chrom=="None":
            self.logger.warning("No chromosome selected, load data via 'Add HDF files'")
            return None

        if not start or not end:
            self.logger.warning("Set values to draw")
            return None

        if start>end:
            end,start = start,end
            start_entry.change_display(end_display)
            end_entry.change_display(start_display)
            self.logger.warning("Switching 'Start' and 'End'")

        self.plotter.draw(chrom,start,end)
        self.canvas.draw()
        #add all lines from plotter in a mapper and initialize popupmenu
        #lines = {k:v[0] for k,v in self.plotter.lines.items()} #<-this will need to be refactored
        lines = self.plotter.lines
        self.popup_axes = PopupMenu(self.master,lines,self.canvas)
        #adding popup menu for artists on track
        ###throws error if I have no track_lines so this fixes it
        if self.plotter.track_lines is not None:
            tracklines ={k:v for k,v in self.plotter.track_lines.items() if v is not None}
            self.popup_track = PopupMenu(self.master,tracklines,self.canvas)

    def _canvas_key_release(self,event,start):
        "zooms in the data"
        start = int(start)
        end = int(event.xdata)
        if start==end:
            self.canvas.mpl_disconnect(self.discid)
            return
        #print(start,end)
        start_com = self.commands["start"]
        end_com = self.commands["end"]
        for value,com in zip([start,end],[start_com,end_com]):
            com.entry.delete(0,tk.END)
            com.entry.insert(0,f"{value}")
        self.canvas.mpl_disconnect(self.discid)
        self.plot_region()

    def savefigure(self):
        "saves figure"
        figure = self.plotter.figure
        f = filedialog.asksaveasfilename(defaultextension=".png")
        if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
            return
        figure.savefig(f,bbox_inches="tight")


    def canvas_key_press(self,event,):
        """binds two method - find_genes, and pop_up based on the button clicked"""
        #make sure that the click is in axes
        if event.inaxes is None:
            return
        #if left double click - show gene
        if event.button==1 and event.dblclick:
            self.find_genes(event.xdata)

        #if left click with release
        if event.button==1 and not event.dblclick:
            ax = event.inaxes

            #add a vertical line there
            start = event.xdata
            canvas_key_release = partial(self._canvas_key_release,start=start)
            #connect to release
            self.discid = self.canvas.mpl_connect("button_release_event",canvas_key_release)
        #if right click and not double
        if event.button==3 and not event.dblclick:
            ax = event.inaxes
            gridaxes = [v for k,v in self.plotter.axes_dict.items() if k!="Track"]
            #if its in grid ax -> popup popup axes
            if ax in gridaxes:
                self.popup_axes.menu.tk_popup(event.guiEvent.x_root,event.guiEvent.y_root)
            else:
                self.popup_track.menu.tk_popup(event.guiEvent.x_root,event.guiEvent.y_root)

        #if right click and double - save
        if event.button==3 and event.dblclick or event.button==2:
            self.savefigure()

    def find_genes(self,xcoord):
        """returns gene found at a given xcoordinate"""
        chrom = self.commands["chrom"].variable.get()
        if self.genes is None:
            return
        df = self.genes.copy()
        if "genomic_accession" in df.columns:
            df = df[df.genomic_accession==chrom]
        if "chrom" in df.columns:
            df = df[df.chrom==chrom]
        df = df[(df.start<xcoord)&(df.end>xcoord)]
        df = df[df.feature=="gene"]
        for _,row in df.iterrows():
            self.logger.info(f"Gene {row.GeneID} of {row['class']} class named '{row['name']}'")

    def mouse_motion(self,event):
        """updates toolbar to show coordinate of a mouse"""
        ax = event.inaxes
        gridaxes = [v for k,v in self.plotter.axes_dict.items() if k!="Track"]
        if not ax:
            #displays nothing
            self.message.set("")
            return 
        x = chr_len_form(event.xdata,None)
        if ax not in gridaxes:
            #destroys the coordinate display
            self.message.set(f"Genome position: {x}")
            return
        
        string = f"Genome position: {x} \t\t Diploid copy number {event.ydata:.1f}"
        self.message.set(string)

def preload(minf):
    "preloads the function"
    minf.status = pd.read_csv("data/metadata.csv")[["Sample Name","dwelling","lineage","region"]]
    minf.genes = pd.read_csv("data/feature.csv")
    minf.assembly = pd.read_csv("data/assembly.csv")
    minf.exons = pd.read_csv("data/exons.csv")
    minf.determine_sample_col()

    minf.parameters["rd"] = {'genes': minf.genes,
                              'exon': minf.exons,
                              'hue': {'Rio Choy': ['Choy01', 'Choy05', 'Choy06', 'Choy09', 'Choy10', 'Choy11', 'Choy12', 'Choy13', 'Choy14'], 'Pachon': ['Pach3', 'Pach7', 'Pach8', 'Pach9', 'Pach11', 'Pach12', 'Pach14', 'Pach15', 'Pach17'], 'Molino': ['Molino2a', 'Molino7a', 'Molino9b', 'Molino10b', 'Molino11a', 'Molino12a', 'Molino13b', 'Molino14a', 'Molino15b'], 'Rascon': ['Rascon02', 'Rascon04', 'Rascon13', 'Rascon15', 'Rascon8', 'Rascon6'], 'Tinaja': ['Tinaja6', 'Tinaja12', 'TinajaB', 'Tinaja2', 'TinajaC', 'Tinaja3', 'TinajaD', 'Tinaja5', 'TinajaE']},
                              'grids': [['Choy01', 'Choy05', 'Choy06', 'Choy09', 'Choy10', 'Choy11', 'Choy12', 'Choy13', 'Choy14', 'Rascon02', 'Rascon04', 'Rascon13', 'Rascon15', 'Rascon8', 'Rascon6'],
                                        ['Pach3', 'Pach7', 'Pach8', 'Pach9', 'Pach11', 'Pach12', 'Pach14', 'Pach15', 'Pach17', 'Molino2a', 'Molino7a', 'Molino9b', 'Molino10b', 'Molino11a', 'Molino12a', 'Molino13b', 'Molino14a', 'Molino15b', 'Tinaja6', 'Tinaja12', 'TinajaB', 'Tinaja2', 'TinajaC', 'Tinaja3', 'TinajaD', 'Tinaja5', 'TinajaE']
                                        ]
                              }

def start_function():
    "starts the function"
    MainInterface().master.mainloop()

