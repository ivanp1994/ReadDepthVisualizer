# -*- coding: utf-8 -*-
"""
Created on Tue May 16 09:01:04 2023

@author: ivanp
"""
import re
import gzip
from functools import partial
from io import StringIO
import requests
import pandas as pd


#GCF = "GCF_023375975.1"
GENOME_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
LINK_PATTERN  = "<a href=\"(?P<link>[^\"]*)\">"

def get_genome_links(text,link_pattern):
    """finds all links in a given pattern"""
    return re.findall(link_pattern,text)

def joint_processing(gcf):
    """
    finding HTML page for a given GCF ID
    the function is then passed to
    "check_for_gcf" which checks if the GCF ID is valid
    or "ncbi_info" which finds all information associated with the ID
    """
    #first enter the parent directory folder
    #via gcf information
    _gcf = gcf[:gcf.find(".")]
    _gcf = _gcf.replace("_","")
    query="/".join([_gcf[i:i+3] for i in range(0, len(_gcf), 3)])
    parent_dir = GENOME_URL+query
    #find the link to the parent directory containing all information
    text = requests.get(parent_dir).text.strip("\n")
    return text,parent_dir

def check_for_gcf(gcf):
    """
    checks if the gcf is valid
    """
    text,_ = joint_processing(gcf)
    if "Error 404" in text:
        return None
    links = get_genome_links(text,LINK_PATTERN)
    #the link that starts with gcf is the ftp repo
    child_dir = [x for x in links if x.startswith(gcf)][0]
    return child_dir.replace(gcf,"")[1:-1]

def ncbi_info(gcf):
    """
    fetches information about the given NCF - to be processed 
    in downstream functions
    """
    text,parent_dir = joint_processing(gcf)
    links = get_genome_links(text,LINK_PATTERN)
    #the link that starts with gcf is the ftp repo
    child_dir = [x for x in links if x.startswith(gcf)][0]
    child_url = parent_dir + "/" + child_dir

    #enter the repository
    text = requests.get(child_url).text.strip("\n")
    links = get_genome_links(text,LINK_PATTERN)

    #the following structures are defined, starting with the name of directory they are in
    #first strip the last slash
    qname = child_dir[:-1]
    return child_url,links,qname

def _get_url(gcf,what):
    """
    common function for getting a particular file
    """
    child_url,links,qname = ncbi_info(gcf)

    ass_url = child_url +  [x for x in links if x.startswith(qname + what)][0]
    return ass_url

get_assembly_report = partial(_get_url,what="_assembly_report")
get_feature_table = partial(_get_url,what="_feature_table")
get_gff_file = partial(_get_url,what="_genomic.gff")

def load_data(url):
    """
    returns data found on a page
    as a decoded string - if data is gunzipped,
    decompresses it
    """
    data = requests.get(url).content
    if url.endswith(".gz"):
        data = gzip.decompress(data)
    return data.decode()

def read_assembly_report(gcf):
    """
    Returns AssemblyReport of a given GCF ID
    """
    url = get_assembly_report(gcf)
    ass = load_data(url).split("\n")
    #ass = requests.get(url).text.split("\n")
    i=0
    while ass[i].startswith("#"):
        i=i+1
    cols = ass[i-1].replace("#","").split("\t")
    #data = "\n".join(ass[i:])
    data = pd.read_csv(StringIO("\n".join(ass[i:])),
                       sep="\t",header=None)
    data.columns = cols
    return data

def read_feature_table(gcf):
    """reads feature table from a given GCF ID"""
    url = get_feature_table(gcf)

    data = load_data(url).split("\n")
    cols = data[0].replace("# ","").split("\t")
    data = pd.read_csv(StringIO("\n".join(data[1:])),
                            sep="\t",header=None)
    data.columns = cols
    return data

def read_exons(gcf):
    """reads only exons"""
    url = get_gff_file(gcf)
    data = load_data(url).split("\n")
    i = 0
    while data[i].startswith("#"):i=i+1
    data = pd.read_csv(StringIO("\n".join(data[i:])),
                   sep="\t",header=None)

    data.columns =  ["chrom","source","feature","start","end","score","strand","frame","attributes"]

    exons = data[data["feature"]=="exon"].copy()
    exons["GeneID"] = exons["attributes"].str.extract("Dbxref=GeneID:([0-9]*)")
    return exons.drop(["attributes"],axis=1)

#ass = get_assembly_report(GCF)
#ft = get_feature_table(GCF)
#gff = get_gff_file(GCF)
#ft = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/375/975/GCF_023375975.1_AstMex3_surface/GCF_023375975.1_AstMex3_surface_feature_table.txt.gz"
#ass = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/375/975/GCF_023375975.1_AstMex3_surface/GCF_023375975.1_AstMex3_surface_assembly_report.txt"
#gff = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/375/975/GCF_023375975.1_AstMex3_surface/GCF_023375975.1_AstMex3_surface_genomic.gff.gz"

#assembly = read_assembly_report(GCF).to_csv("data/assembly.csv")
#feature = read_feature_table(GCF).to_csv("data/feature.csv")
#exons = read_exons(GCF).to_csv("data/exons.csv")


