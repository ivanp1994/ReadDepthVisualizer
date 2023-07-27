# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 08:57:43 2023
Supplemental widgets
@author: ivanp
"""
import re
import tkinter as tk

def convert_fancy(value):
    """
    Converts the 'fancy' chromosome display to numerical coordinate
    e.g.
    >>> 20M30K163
    >>> 20 * 1e6 + 30 * 1e3 + 163

    """
    is_invalid = False
    true_value = 0
    for denom,exp in zip(["m","k"],[1e6,1e3]):
        pattern = "(?P<m>[0-9]+)" + f"[{denom}|{denom.upper()}]"
        intvalue = re.findall(pattern,value)
        if len(intvalue)>1:
            is_invalid = True
        if len(intvalue)==0:
            continue
        intvalue = intvalue[0]
        is_invalid = not intvalue.isdigit()
        if is_invalid:
            break
        intvalue = int(intvalue)
        value = value.replace(f"{intvalue}{denom}","").replace(f"{intvalue}{denom.upper()}","")
        true_value = true_value + intvalue * exp
    if not value:
        value = "0"
    if is_invalid or not value.isdigit():
        return None

    return  int(true_value + int(value))

def convert_tointeger(value:str):
    """
    Converts text based value
    to an integer.

    Takes care of:
        > commas and spaces
        > fancy formatting (e.g. M and K)
    """
    converted = value.replace(" ","").replace(",","")
    if converted.isdigit():
        return int(converted)
    converted = convert_fancy(converted)
    return converted

def convert_todisplay(converted:int):
    """
    Converts a number to better display:
        e.g. 150000 to 150 000
    (just adds spaces every 1K)
    """
    #convert to string and reverse order
    if converted is None:
        return ""
    display = str(converted)[::-1]
    #segment into list of 3 items, then reverse the list
    display = [display[i:i+3][::-1] for i in range(0, len(display), 3)][::-1]
    #join with spaces
    return " ".join(display)

class TextEntry():
    """
    Text entry
    """
    def __init__(self,master,label):
        self.variable = tk.StringVar()
        self.master = master
        self.label = tk.Label(master=master,text=label,anchor="w",relief="raised",bg="gainsboro")
        self.entry = tk.Entry(master=master)

    def place(self,r,c,sticky="nsew"):
        "places both label and entry"
        self.label.grid(row=r, column=c, sticky=sticky)
        self.entry.grid(row=r, column=c+1, sticky=sticky)

class LociEntry(TextEntry):
    """
    Specific entry for locus
    """

    def get_loci(self):
        "gets the loci in numerical value and in displayed value"
        value = self.entry.get() #get value
        converted = convert_tointeger(value) #convert to int
        display = convert_todisplay(converted) #convert to display
        self.change_display(display)
        return converted, display

    def change_display(self,display:str):
        """change display"""
        self.entry.delete(0,tk.END)
        self.entry.insert(0,display)

class GCFcheck(TextEntry):
    """
    Subclass of GCF used solely to search NCBI for a given GCF ID
    """
    def __init__(self, master, logger):
        super().__init__(master=master,label="GCF")
        self.entry.bind("<Return>",self.search_gcf)
        self.logger = logger

    def search_gcf(self,*args):
        "searches NCBI for a given GCF"
        gcf = self.entry.get()
        isvalid = gcf.check_for_gcf(gcf)
        if isvalid is not None:
            self.variable.set(gcf)
            self.logger.info(f"GCF corresponds to {isvalid}")
        else:
            self.logger.warning(f"No valid assembly for {gcf}")
            self.entry.delete(0, 'end')
            self.variable.set("")

class OptionMenu():
    """
    self.variable ->variable of menu
    self.label    ->tk.Label
    self.optionmenu ->tk.Option menu
    """
    def place(self, r, c,sticky="nsew"):
        "Places label and option menu"
        self.label.grid(row=r, column=c, sticky=sticky)
        self.optionmenu.grid(row=r, column=c+1, sticky=sticky)


    def __init__(self, master, label, options, typevar="string"):
        """
        Parameters
        ----------
        master : master for the button, a frame.
        label : string.
        options : list for choices.
        typevar : None, "bool", "int" if None then defaults to string

        Returns
        -------
        """
        label_ent = tk.Label(master=master, text=label,
                             anchor="w", relief="raised", bg="gainsboro")
        self.label = label_ent
        self.variable = tk.StringVar()
        if typevar == "bool":
            self.variable = tk.BooleanVar()
        if typevar == "int":
            self.variable = tk.IntVar()
        if typevar == "float":
            self.variable = tk.DoubleVar()

        self.variable.set(options[0])
        self.optionmenu = tk.OptionMenu(master, self.variable, *options)

    def repopulate_options(self,options):
        """
        repopulates option menu
        """

        om = self.optionmenu["menu"]
        #print(type(om))
        #print(dir(om))
        om.delete(0, "end")
        if options is None:
            return 
        if isinstance(options,str):
            options = [options]
        for string in options:
            om.add_command(label=string,command=lambda value=string: self.variable.set(value))

class Slider():
    """
    Just what it says
    """

    def place(self, r, c):
        """
        Places slider on specified row and column
        """
        self.label.grid(row=r, column=c, sticky="nsew")
        self.scale.grid(row=r, column=c+1, sticky="nsew")

    def __init__(self, master, label,default=1,lw=1,up=10):

        variable = tk.DoubleVar()
        self.variable = variable
        self.master = master
        label_ent = tk.Label(master=master, text=label,
                             anchor="w", relief="raised", bg="gainsboro")
        self.label = label_ent
        scale = tk.Scale(master=master, from_=lw, to=up,
                         variable=self.variable, orient="horizontal")
        self.variable.set(default)
        self.scale = scale
