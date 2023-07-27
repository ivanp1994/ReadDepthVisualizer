# -*- coding: utf-8 -*-
"""
Created on Fri May 19 12:22:19 2023
ReadDepth drawer - to be used in interactive tkinter session.
For details, see docs on ReadDepthDrawer
@author: Ivan
"""

import logging
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from .genotyping import get_rd_region

matplotlib.style.use("ggplot")
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

def hide_ticks(ax,axis="xy"):
    """hides ticks"""
    if "x" in axis:
        ax.xaxis.set_tick_params(length=0)
        plt.setp(ax.get_xticklabels(),visible=False)
    if "y" in axis:
        ax.yaxis.set_tick_params(length=0)
        plt.setp(ax.get_yticklabels(),visible=False)

#CONSTRUCTING DRAW_DATAFRAME
#HUE MAPPING
def get_hue_mapping(hues,loaded_samples):
    """
    Hues -> dictionary of lists containing samples
    Status -> index must be samples
    Returns
        Legend dictionary (label:color)
        Sample dictionary (sample:color)
    """
    color_legend = dict(zip(hues.keys(),COLORS))
    allsamples = list()
    allcolors = list()
    for label,samps in hues.items():
        subsamples = [x for x in samps if x in loaded_samples]
        allsamples = allsamples + subsamples
        allcolors = allcolors + [color_legend[label]]*len(subsamples)
    color_mapping = dict(zip(allsamples,allcolors))
    return color_legend, color_mapping
#EXAMPLE  #color_legend,color_mapping = get_hue_mappings(rdparams["hue"],samples)

def construct_drawdataframe(rdparams,loaded_samples):
    """
    draw dataframe:
        index -> samples
        grids -> str of Grid IDs
        color -> color
    legend:
        mapping -> color
    """
    legend, color_mapping = get_hue_mapping(rdparams["hue"],loaded_samples)
    if len(legend)==1:
        legend = None

    grid_mapping = {f"Grid {i+1}":grid for i,grid in enumerate(rdparams["grids"])}

    drawframe = list()
    for grid, samples in grid_mapping.items():
        samples = [x for x in samples if x in loaded_samples]
        df = pd.DataFrame(index=samples)
        df["grid"] = grid
        df["color"] =  df.index.to_series().map(color_mapping)
        drawframe.append(df)
    drawframe = pd.concat(drawframe,axis=0)
    return drawframe, legend
#EXAMPLE    #df, leg = construct_drawdataframe(rdparams,samples)

#DATAFRAME OPERATIONS FOR CONSTRUCTING DRAW_DATAFRAME2
def group_overlapping(data,chromosomes=None,cols=None):
    """
    Merges overlapping intervals to get merged data

    Adapted from https://stackoverflow.com/questions/57882621

    Parameters
    ----------
    data : pd.DataFrame.
    chromosomes : iterable, optional
        List of chromosomes on which to merge.
        The default is None meaning it merges on all unique
        chromosomes of found in "chrom" column
    cols : iterable, optional
        The locus parts of dataframe,
        where chromosomes are, and starts and ends are
        The default is None meaning ["chrom","start","end"]

    Returns
    -------
    merged dataframe

    """
    if cols is None:
        cols = ["chrom","start","end"]
    _chrom_col,_start_col,_end_col = cols
    if chromosomes is None:
        chromosomes = data[_chrom_col].unique().tolist()

    #sort the data - without this, nothing works
    data = data.sort_values([_chrom_col,_start_col])

    result_df = list()
    for chrom in chromosomes:
        input_df = data[data[_chrom_col]==chrom].copy()
        input_df[_start_col]=input_df[_start_col].astype(float).astype(int)
        input_df[_end_col]=input_df[_end_col].astype(float).astype(int)
        input_df["group"]=(input_df[_start_col]>input_df[_end_col].shift().cummax()).cumsum()
        for _,df in input_df.groupby("group"):
            #df["ypos"] = [i for i in range(len(df))]
            df["ypos"] = list(range(len(df)))
            result_df.append(df)

    return pd.concat(result_df,axis=0,ignore_index=True)

def find_overlaping(df,chrom,start,end):
    "finds overlap between two - used to draw exons and introns"
    #if there is no df
    if len(df)==0:
        return pd.DataFrame()

    if "genomic_accession" in df.columns:
        df = df[df.genomic_accession==chrom]
    if "chrom" in df.columns:
        df = df[df.chrom==chrom]
    return  df[(df.start<=end)&(df.end>=start)]

def get_overlapping_features(chrom,start,end,rdparams):
    """
    Gets intersecting genes and exons from a given position.
    Genes that overlap are jittered so they do not overlap each other
    on the plot. This returns dataframe with 'ypos' column that
    determines where the y position will be
    """
    #load genes and exons
    genes = rdparams.get("genes",pd.DataFrame())
    exons = rdparams.get("exon",pd.DataFrame())

    #process genes for easier acces
    if len(genes)>0:
        genes = genes[genes["feature"]=="gene"]
        if "genomic_accession" in genes.columns:
            genes = genes.rename({"genomic_accession":"chrom"},axis=1)

    #find overlaps
    fgenes = find_overlaping(genes,chrom,start,end)
    fexons = find_overlaping(exons,chrom,start,end)

    #if genes overlap - change Y position via grouping
    # group overlapping -> modifies dataframe to include "ypos" column
    if len(fgenes)>0:
        fgenes = group_overlapping(fgenes)
        # exons cannot be without genes
        if len(fexons)>0:
            fexons["ypos"] = fexons["GeneID"].map(dict(zip(fgenes["GeneID"],fgenes["ypos"])))
        final = pd.concat([fgenes,fexons],axis=0)
        #increment ypos
        final["ypos"] = final["ypos"] + 1
    #if there are no genes -> there are no exons
    #nothing to draw
    else:
        final = None
    return final

# DRAWING ONTO THE TRACK AX
def draw_quiver(row,ax):
    """
    draws arrow denoting the direction of feature
    """
    start = row.start
    end = row.end
    ypos = row.ypos
    xdir = 1
    if row.strand=="-":
        xdir = -1
    mean = (start+end)/2
    #xpos = [int(mean*x) for x in [1-_delta,1,1+_delta]]
    xpos = [int(mean)]
    xdir = [xdir] * len(xpos)
    ydir = [0] * len(xpos)
    ypos = [ypos] * len(xpos)
    return ax.quiver(xpos,ypos,xdir,ydir,zorder=1,
              headaxislength=3,headlength=3,
              alpha=0.8,minlength=0,headwidth=3)

def annotate_gene(row,ax,gene_name="GeneID"):
    """places gene name above the gene"""
    start = row.start
    end = row.end
    ypos = row.ypos
    gene_name = row[gene_name]

    mean = (start+end)/2
    return ax.annotate(gene_name,xy=(mean,ypos+0.3),weight="bold",
                annotation_clip=False,ha="center")

def log_genes(feat,start,end):
    """shoves information about genes into logging module"""
    if feat is None:
        logging.info(f"No genes within range {start}-{end}")
        return

    gns = feat.loc[feat["feature"]=="gene"]
    if len(gns)==1:
        logging.info(f"One gene is found in range {start}-{end}")
    else:
        logging.info(f"{len(gns)} genes are found in range {start}-{end}:")
    #for _,r in gns.iterrows():
    #    logging.info(f"Gene {r.GeneID} of {r['class']} class named '{r['name']}'")

def draw_track(trackax,chrom,start,end,rdparams,lw=2,gene_name="GeneID",):
    """Draws genes intersecting a range on a given trackax"""
    feat = get_overlapping_features(chrom,start,end,rdparams)
    log_genes(feat,start,end) #logs genes
    if feat is None:
        #kill the track ax
        trackax.axis("off")
        return None,None,None,None
    #make sure that track axis is on
    trackax.axis("on")
    feat["lw"] = feat["feature"].map({"gene":lw,"exon":3*lw})
    trackax.set_ylim(0,max(feat["ypos"])+1)

    genes = feat[feat["feature"]=="gene"]
    exons = feat[feat["feature"]=="exon"]
    gene_lines, strand_lines, exon_lines, annotations = None, None, None, None
    if len(genes)>0:
        #store gene lines
        gene_lines = trackax.hlines(genes["ypos"],genes["start"],genes["end"],lw=genes["lw"])
        #draw the arrow as a strand
        strand_lines = genes.apply(draw_quiver,axis=1,ax=trackax,).tolist()
        #annotate - no need to continue
        annotations = genes.apply(annotate_gene,axis=1,ax=trackax,gene_name=gene_name).tolist()
    if len(exons)>0:
        #store exon lines
        exon_lines = trackax.hlines(exons["ypos"],exons["start"],exons["end"],lw=exons["lw"])
    return gene_lines,strand_lines,exon_lines, annotations

# CREATING FIGURES
def chr_len_form(x,pos):
    "formats Chromosome view to be nicer"
    #pylint - pos is necessary for matplotlib's interface
    if x >= 1e6:
        return "%1.2f Mb" % (x*1e-6)
    if x >= 1e3:
        return "%1.2f Kb" % (x*1e-3)
    return str(x)

class ReadDepthDrawer():
    """
    Plotter for drawing read depth,
    optimized to work with interactive matplotlib backends

    It can also work within Python's script, but then some modifications are necessary
    Three main inputs are:
        rdfile -> HDF object that contains read-depth info
            (see "reading.py")
        samples -> list of samples (unique) that will be drawn
        rdparameters -> dictionary containing parameters

    "rdparameters" consists of the following values:
        "genes" : pd.DataFrame containing information about genes
        "exons" : pd.DataFrame containing information about exons
        NOTE: both of those can be None if track ax is not desired
        (see "gcf_fetcher.py" for details)
            - both of those files must contain columns :
                "feature","chrom","start","end", "GeneID"
                - valid replacement for "chrom" is "genomic_accession"
                - "GeneID" determines how gene will be annotated
                    (what will be displayed on the track axis)
                - different identifier can be selected, in that
                case just pass 'gene_name="{column}"' to .draw method of the object
        "grids" : list of samples that will be drawn, it's a list of lists, for example:
            >>> [["sampleA","sampleB"],["sampleC","sampleD"],["sampleA","sampleD"]]
            will draw three grids, with 'sampleA' and 'sampleD' being drawn twice, but
            on the different grids
        "hue" : None or dictionary, if a dictionary is passed, then it must have form
            {"label1":[samples...],"label2",[samples...]}, for example:
            >>> {"tumor":["sampleA","sampleC"],"healthy":["sampleB","sampleD"]}

    After object is initialized, call ".draw" method. Draw method takes the following arguments:
        "chrom" : str -> chromosome part of loci
        "start" int -> starting position
        "end": int -> ending position
    Additional arguments are:
        "lw": int -> thickness of lines (default is 2)
        "gene_name":str -> column in "genes" that will be displayed in the track ax


    Example usage:
        >>> rdfile = "data.rd"
        >>> samples = {'Rascon04', 'Molino10b', 'Rascon6', 'Choy06', 'Rascon02', 'Pach8', 'Molino7a', 'Pach3', 'Molino9b', 'Tinaja5', 'Choy05', 'Tinaja3', 'Choy01', 'Tinaja2', 'Pach7'}
        >>> rdparams = {'genes': pd.read_csv("data/feature.csv"),
        >>>                           'exon': pd.read_csv("data/exons.csv"),
        >>>                           'hue': {'Rio Choy': ['Choy01', 'Choy05', 'Choy06', 'Choy09', 'Choy10', 'Choy11', 'Choy12', 'Choy13', 'Choy14'], 'Pachon': ['Pach3', 'Pach7', 'Pach8', 'Pach9', 'Pach11', 'Pach12', 'Pach14', 'Pach15', 'Pach17'], 'Molino': ['Molino2a', 'Molino7a', 'Molino9b', 'Molino10b', 'Molino11a', 'Molino12a', 'Molino13b', 'Molino14a', 'Molino15b'], 'Rascon': ['Rascon02', 'Rascon04', 'Rascon13', 'Rascon15', 'Rascon8', 'Rascon6'], 'Tinaja': ['Tinaja6', 'Tinaja12', 'TinajaB', 'Tinaja2', 'TinajaC', 'Tinaja3', 'TinajaD', 'Tinaja5', 'TinajaE']},
        >>>                           'grids': [['Choy01', 'Choy05', 'Choy06', 'Choy09', 'Choy10', 'Choy11', 'Choy12', 'Choy13', 'Choy14', 'Rascon02', 'Rascon04', 'Rascon13', 'Rascon15', 'Rascon8', 'Rascon6'],
        >>>                                     ['Pach3', 'Pach7', 'Pach8', 'Pach9', 'Pach11', 'Pach12', 'Pach14', 'Pach15', 'Pach17', 'Molino2a', 'Molino7a', 'Molino9b', 'Molino10b', 'Molino11a', 'Molino12a', 'Molino13b', 'Molino14a', 'Molino15b', 'Tinaja6', 'Tinaja12', 'TinajaB', 'Tinaja2', 'TinajaC', 'Tinaja3', 'TinajaD', 'Tinaja5', 'TinajaE']]
        >>>                           }
        >>>
        >>> chrom = "NC_064408.1"
        >>> start = 41000
        >>> end = 58000
        >>> figure = ReadDepthDrawer(rdparams,samples,rdfile).draw(chrom,start,end)


    """
    def __init__(self,rdparameters,samples,rdfile):
        self.has_trackax = False
        self.rdfile = rdfile
        self.draw_dataframe, self.legend = construct_drawdataframe(rdparameters,samples)
        self.figure = None
        self.genes = rdparameters["genes"]
        self.exons = rdparameters["exon"]
        self.has_trackax = self.check_trackax(rdparameters)
        self.lines = None #this is a mapping of SAMPLE -> MATPLOTLIB LINE
        self.track_lines = None #mapping artists on TRACKAX
        self.scale = rdparameters.get("scale","log")

        if self.genes is None or self.genes is False:
            self.genes = pd.DataFrame()
        if self.exons is None or self.exons is False:
            self.exons = pd.DataFrame()

        self.init_figure(rdparameters)

    @staticmethod
    def check_trackax(rdparameters):
        """checks if there is a trackax -> trackax for genes"""
        genes = rdparameters.get("genes",None)
        exons = rdparameters.get("exons",None)

        has_trackax = True
        if genes is None and exons is None:
            has_trackax = False
        return has_trackax

    def construct_grid(self,rdparameters):
        "constructs grid -> no matplotlib initialization yet"
        grids = list()
        ratios = list()
        if self.has_trackax:
            grids.append("Track")
            ratios.append(2)
            self.has_trackax = True
        grids = grids + [f"Grid {i+1}" for i,_ in enumerate(rdparameters["grids"])]
        ratios = ratios + [5]*len(rdparameters["grids"])

        return grids, ratios

    def init_figure(self,rdparameters):
        "initializes figure"
        grids, ratios = self.construct_grid(rdparameters)
        fig, axes = plt.subplots(nrows=len(grids), ncols=1, sharex='col',
                               gridspec_kw={"height_ratios": ratios,
                                            "hspace":0.05},
                               figsize=(5*len(grids), 5))
        if len(grids)==1:
            axes = [axes]
        self.axes_dict=dict(zip(grids,axes))
        self.draw_dataframe["grid"] = self.draw_dataframe["grid"].map(self.axes_dict)

        #make Grid Axes share Y
        gridaxes = [v for x,v in self.axes_dict.items() if x!="Track"]
        gridaxes[0].get_shared_y_axes().join(gridaxes[0], *gridaxes[1:])

        self.figure = fig

    @staticmethod
    #wrapped around in a list to avoid bumping off pandas df
    def _get_rd_region(sample,hdf_file,chromosome,start,end):
        """gets rd region"""
        return [get_rd_region(sample,hdf_file,chromosome,start,end,True)]

    @staticmethod
    def draw_step(row,lw=4):
        """drawstep"""
        ax = row["grid"]
        color = row["color"]
        # values are wrapped around list to avoid bumping pandas df
        values = row["values"][0]
        return ax.step(values.index, values, color=color,linewidth=lw) #lw = 4

    @staticmethod
    def pair_sample_lines(drawdf):
        """
        pairs every sample to its line
        since same samples can be drawn on multiple axes,
        this function is necessary
        """
        lines = drawdf.reset_index().groupby("index")["lines"].agg(list)
        lines = lines.apply(lambda x: [item for sub in x for item in sub])
        return lines.to_dict()

    def draw(self,chrom,start,end,**kwargs):
        """draws - the main method"""
        #clear ax
        self.clear_ax()
        #draw on track ax
        #if self.axes_dict.get("Track",None) is not None:
        if self.has_trackax:
            glines,slines,xlines,annots = draw_track(self.axes_dict.get("Track"),
                                              chrom,start,end,
                                              {"genes":self.genes,"exon":self.exons},**kwargs)

            hide_ticks(self.axes_dict.get("Track"),"y")
            self.axes_dict.get("Track").set_ylabel("Track",rotation=90,fontsize=15)
            self.track_lines = {"genes":glines,"strand":slines,"exons":xlines,"annotations":annots}

        #draw on grid
        draw_df =  self.draw_dataframe.copy()
        draw_df["values"] = draw_df.index.to_series().apply(self._get_rd_region,
                                                            hdf_file=self.rdfile,
                                                            chromosome=chrom,
                                                            start=start, end=end)

        draw_df["lines"] = draw_df.apply(self.draw_step, axis=1,lw=2)
        self.lines = self.pair_sample_lines(draw_df)

        #beautify axes and return legend
        last_grid = list(self.axes_dict)[-1]
        legend = None
        for axname,ax in self.axes_dict.items():
            if axname!=last_grid:
                ...#hide x labels
                #ax.xaxis.set_visible(False) -> this turns off the grid
                ax.tick_params(labelbottom=False) #-> this works
            else:
                legend = self.add_legend(ax)
            self.beautify_ax(ax)
        #set limits to start and end
        ax.set_xlim(start,end)
        #store legend in tracklines
        if legend is not None and self.track_lines is not None:
            self.track_lines["legend"] = legend
        return self.figure

    @staticmethod
    def beautify_ax(ax):
        """beautifes axis"""
        ax.xaxis.set_tick_params(labelsize=15)
        ax.yaxis.set_tick_params(labelsize=15)
        #grids
        ax.set_axisbelow(True)
        ax.xaxis.grid(which="major",visible=True,zorder=0)
        ax.yaxis.grid(which="major",visible=True,zorder=0)
        #remove spines
        ax.spines[["top","right",]].set_visible(False)
        #set nicer formatter
        ax.xaxis.set_major_formatter(chr_len_form)
        #remove ticks - rely only on grid
        ax.xaxis.set_tick_params(length=0)
        ax.yaxis.set_tick_params(length=0)

    def clear_ax(self):
        "clears axes"
        for ax in self.axes_dict.values():
            ax.clear()

    def add_legend(self,ax):
        "adds legend"
        if self.legend is None:
            return None
        handles = [matplotlib.lines.Line2D([],[], color = v,label=k) for k,v in self.legend.items()]
        return ax.legend(handles=handles,loc="upper right")

    def get_grid_axes(self):
        "returns grid axes"
        return [v for k,v in self.axes_dict.items() if k!="Track"]
#%% MAIN
"""
rdfile = "data.rd"
samples = {'Rascon04', 'Molino10b', 'Rascon6', 'Choy06', 'Rascon02', 'Pach8', 'Molino7a', 'Pach3', 'Molino9b', 'Tinaja5', 'Choy05', 'Tinaja3', 'Choy01', 'Tinaja2', 'Pach7'}
rdparams = {'genes': pd.read_csv("data/feature.csv"),
                          'exon': pd.read_csv("data/exons.csv"),
                          'hue': {'Rio Choy': ['Choy01', 'Choy05', 'Choy06', 'Choy09', 'Choy10', 'Choy11', 'Choy12', 'Choy13', 'Choy14'], 'Pachon': ['Pach3', 'Pach7', 'Pach8', 'Pach9', 'Pach11', 'Pach12', 'Pach14', 'Pach15', 'Pach17'], 'Molino': ['Molino2a', 'Molino7a', 'Molino9b', 'Molino10b', 'Molino11a', 'Molino12a', 'Molino13b', 'Molino14a', 'Molino15b'], 'Rascon': ['Rascon02', 'Rascon04', 'Rascon13', 'Rascon15', 'Rascon8', 'Rascon6'], 'Tinaja': ['Tinaja6', 'Tinaja12', 'TinajaB', 'Tinaja2', 'TinajaC', 'Tinaja3', 'TinajaD', 'Tinaja5', 'TinajaE']},
                          'grids': [['Choy01', 'Choy05', 'Choy06', 'Choy09', 'Choy10', 'Choy11', 'Choy12', 'Choy13', 'Choy14', 'Rascon02', 'Rascon04', 'Rascon13', 'Rascon15', 'Rascon8', 'Rascon6'],
                                    ['Pach3', 'Pach7', 'Pach8', 'Pach9', 'Pach11', 'Pach12', 'Pach14', 'Pach15', 'Pach17', 'Molino2a', 'Molino7a', 'Molino9b', 'Molino10b', 'Molino11a', 'Molino12a', 'Molino13b', 'Molino14a', 'Molino15b', 'Tinaja6', 'Tinaja12', 'TinajaB', 'Tinaja2', 'TinajaC', 'Tinaja3', 'TinajaD', 'Tinaja5', 'TinajaE']]
                          }

chrom = "NC_064408.1"
start = 41000
end = 58000
o = ReadDepthDrawer(rdparams,samples,rdfile)
o.draw(chrom,start,end)

#ax = o.axes_dict["Grid 1"]
"""