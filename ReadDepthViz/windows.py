# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 12:29:04 2023
Supplemental toplevel windows
@author: ivanp
"""
import tkinter as tk
from .widgets import GCFcheck, OptionMenu, Slider

class DataParameterWindow():
    """
    Window for selecting parameters of data
    """
    def __init__(self,master_object):

        self.master = tk.Toplevel(master=master_object.master)
        self.logger = master_object.logger
        self.data_variables = dict()
        #
        self.master.title = "Parameter Window"
        dataframe = tk.Frame(master=self.master,)
        dataframe.grid(row=0,column=0)
        compress = OptionMenu(dataframe,"Compress",["half","normal","double"],)
        lower = Slider(dataframe,"Lower bound",4,2,8)
        upper = Slider(dataframe,"Upper bound",5,2,8)
        gcf = GCFcheck(dataframe,self.logger)
        #
        self.data_variables["compress"] = compress.variable
        self.data_variables["lower_ratio"] = lower.variable
        self.data_variables["upper_ratio"] = upper.variable
        self.data_variables["gcf"] = gcf.variable
        #
        compress.place(0,0)
        lower.place(1,0)
        upper.place(2,0)
        gcf.place(3,0)

class GridParameterWindow():
    """
    Window for setting parameters of the drawn plot
    """
    def __init__(self,status):

        self.options = self.populate_options(status)
        self.gridvalues = self.populate_values(status,self.options)
        self.hueoptions = self.populate_huevalues(status,self.options)

        self.master = tk.Toplevel()
        self.option_menus = list()
        self.show_exon = tk.IntVar()
        self.show_gene = tk.IntVar()
        self.hue_var = None

        self.extend_options()
        #frame = tk.Frame(master=self.master)
        #frame.grid(row=0,column=0)


        self.show_gene.set(1)
        tk.Checkbutton(master=self.master,text="Gene",variable=self.show_gene,onvalue=1, offvalue=0,
                       justify="left").grid(
            row=0,column=0,columnspan=2,sticky="nsew")
        tk.Checkbutton(master=self.master,text="Exon",variable=self.show_exon,onvalue=1, offvalue=0).grid(
            row=1,column=0,columnspan=2,sticky="nsew")

        hue_menu = OptionMenu(self.master,"Hue",options=["None"]+status.columns.tolist())
        hue_menu.place(r=2,c=0)
        self.hue_var = hue_menu.variable

        tk.Button(master=self.master,command=self.add_gridopt,text="Add grid").grid(
            row=3,column=0,columnspan=1,sticky="nsew")
        tk.Button(master=self.master,command=self.delete_gridopt,text="Remove grid").grid(
            row=3,column=1,columnspan=1,sticky="nsew")

        self.add_gridopt()

    def extract_parameters(self):
        "extracts parameters"
        param_dic = {"genes":bool(self.show_gene.get()),
        "exon":bool(self.show_exon.get()),
        "hue":self.hue_var.get()
        }
        param_dic["hue"] = self.hueoptions.get(param_dic["hue"],"None")

        if param_dic["hue"]=="None":
            param_dic["hue"]={"all":self.gridvalues["--All--"]}


        param_dic["grids"] = [self.gridvalues[m.variable.get()] for m in self.option_menus]
        return param_dic

    def extend_options(self):
        "Add one option that allows for the entire dataset"
        allsamps = list(self.gridvalues.values())
        allsamps = [item for sublist in allsamps for item in sublist]
        allsamps = sorted(list(set(allsamps)))
        self.options.insert(0,"--All--")
        self.gridvalues["--All--"] = allsamps

    def add_gridopt(self):
        """adds gridopt menu"""
        i = len(self.option_menus)+1
        option_menu = OptionMenu(self.master,f"Grid {i}",self.options)
        self.option_menus.append(option_menu)
        option_menu.place(r=i+3,c=0)

    def delete_gridopt(self):
        """destroys a gridopt menu"""
        if len(self.option_menus)==1:
            return None
        menu = self.option_menus.pop()
        menu.label.destroy()
        menu.optionmenu.destroy()

    @staticmethod
    def populate_options(status):
        "repopulates options"
        options = list()
        for _c in status.columns:
            options = options + [f"{_c}:{x}" for x in status[_c].unique()]
        return options

    @staticmethod
    def populate_values(status,options):
        "populates values"
        return  {x:list(status[status[x[:x.find(":")]]==x[x.find(":")+1:]].index.unique()) for x in options}

    @staticmethod
    def populate_huevalues(status,options):
        "populates values for hue"
        dic = GridParameterWindow.populate_values(status,options)
        hueoptions = {}
        for hue in status.columns:
            hueoptions[hue] = {k.replace(f"{hue}:",""):v for k,v in dic.items() if k.startswith(hue)}
        return hueoptions

class PopupMenu():
    """
    Popup menu for axes
    It contains a mapping of the type
    >>> {label[str]:artist[matplotlib.artist]}
    It's solely used to toggle the visibility of the artist
    """
    def __init__(self,master,option_dic,canvas):
        self.menu = tk.Menu(master=master,tearoff=False)
        self.label_to_line = option_dic
        self.label_to_var = dict()
        self.canvas = canvas

        for label,lines in self.label_to_line.items():
            #in instances where lines is a collection of artists
            #then first element of artists will dictate
            #the state for the rest
            if isinstance(lines,list):
                lines = lines[0]

            var = tk.BooleanVar(value=lines._visible,name=label)
            var.set(lines._visible)
            var.trace_add(mode="write",callback=self.change)
            self.menu.add_checkbutton(label=label, onvalue=1, offvalue=0, variable=var)
            self.label_to_var[label] = var
            #print(label,lines._visible,var.get())

    def change(self,var,index,mode):
        "changes the variable on/off and redraws the axis"
        #change to visible
        line = self.label_to_line[var]
        visible = self.label_to_var[var].get()
        if isinstance(line,list):
            for l in line:
                l._visible = visible
        else:
            line._visible = visible
        #self.label_to_line[var]._visible = self.label_to_var[var].get()
        self.canvas.draw()
