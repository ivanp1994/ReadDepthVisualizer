# -*- coding: utf-8 -*-
"""
This script collates read depth information from a list of
CNVPytor .hdf5 files and merges the relevant information into one .hdf5 file.
===============================================================================
To collate all files from CNVPytor's HDF files, use the function
    "write_read_depth_file"

===============================================================================
    PARAMETERS
===============================================================================
pytor_files:
    Iterable of files (paths to the files) for the CNVPytor's HDF files

output_file:
    Path of the output file

assembly_report (OPTIONAL):
    Path to the assembly report of a reference genome.
    If it's provided, two things will happen:
        - only assembled chromosomes will be recorded
        - every chromosome will get simple name in its attribute (.attrs[chrom]),
            e.g. instead of "NC_064431.1", its name will be "chrom 24"

binsize_dic (OPTIONAL):
    Python dictionary. The keys of the dictionary are the files outlined in
    pytor_files or the names of those files (without the directory and extension),
    and its values are 100 divisible integers representing histogram partition
    bin size.

    If its not passed, binsize will be detected automatically
    (see "lower_ratio" and "upper_ratio")



flag_dic (OPTIONAL):
    Python dictionary. The keys of the dictionary are the files outlined in
    pytor_files or the names of those files (without the directory and extension),
    and its values represent "flags" in CNVPytor's syntax.
    E.g. "GC" means corrected for GC.

    If its not passed, flag will be detected automatically
    as the largest suffix of the signal.

lower_ratio, upper_ratio (OPTIONAL):
    For a particular bin size, after partition, every
    bin has its own READ DEPTH. If a chromosome of length 1 000 is
    divided into bins of 100 size, there will be 10 bins. The
    relevant statistics are the MEAN and STANDARD DEVIATION of the
    READ DEPTH array. Per CNVnator's best practices, optimal binsize
    is when the ratio of MEAN to STANDARD DEVIATION is between 4 and 5.

===============================================================================
    OUTPUT
===============================================================================

The output is a structured hdf5 file.

The resulting HDF5 file's main keys are samples.

Every sample has three attributes (accessed via .attrs["x"]):
    - "flag" -> details for CNVpytor (example _GC implies that rd is GC adjusted,etc)
    - "binsize" -> the binsize used in RD plots
    - "mean" -> the mean read-depth

Within every sample, there are datasets that correspond to chromosomes.
Every sample has one attribute (accessed via .attrs["x"]):
    - "name" -> How assembly report codes this molecule - for example chrom 8

The values in chromosomes represent diploid read depth.
Meaning that the raw values in histogram partition are divided by
the mean read depth and multiplied by 2. For example:
    [200,100,400,200,...] with mean 200 becomes
    [2  ,1  ,4  ,2  ,....]

To retrieve values per chromosome, simply pass the following
    "np.array(hdf_file_buffer[sample][chromosome])"
    where:
        - hdf_file_buffer is the opened HDF5PY buffer
        - sample is a sample
        - chromosome is the chromosome

===============================================================================
     EXAMPLE OUTPUT
===============================================================================
The following code snippet will take
a chromosome and sample, and reconstruct the entire dataframe

sample = "Choy01"
chrom = "NC_064408.1"
chr_len = 134019835 #or dont pass it
chrom_df = pd.DataFrame()

with h5py.File(OUTPUT_RD_FILE) as _f:

    _bin_size = _f[sample].attrs["binsize"]
    _rd_mean = _f[sample].attrs["mean"]
    chrom_name = _f[sample][chrom].attrs["name"]

    chrom_df["NORM_RD"] = _f[sample][chrom]

chrom_df["CHROM_START"]= chrom_df.index * _bin_size
chrom_df["CHROM_END"] = (chrom_df.index + 1) * _bin_size
if chr_len:
    chrom_df.loc[len(chrom_df)-1,"CHROM_END"] = chr_len
chrom_df["RAW_RD"] = chrom_df["NORM_RD"] * _rd_mean

@author: ivanp
"""

import re
import os
from pathlib import Path
import time
import logging
import h5py as hpy
import numpy as np
import pandas as pd

COMPRESS_DICT = {"half":np.float16,
                 "norm":np.float32,
                 "double":np.float64}

def estimate_samplebinsize_logfile(sample,blogpath="",blogdir="",blogext=".binsize.log"):
    """
    Gets binsize of sample based on the LOG file created
    during CNVPytor CNV calling-

    Its a given that the LOG file corresponds to the sample name:
        e.g. sample = "polnareff"
            logfile = "polnareff.binsize.log"
    If no, you can specify the path to the logfile with "blogpath" parameter.


    Parameters
    ----------
    sample : str, sample.
    blogpath : str,
        path to the LOG file of the sample
    blogdir : str, optional
        Directory where the LOG file is found
    blogext : str, optional
        Extension of the file. The default is ".binsize.log".

    Raises
    ------
    ValueError
        When there is no logfile.

    Returns
    -------
    binsize : int
        Binsize of the sample.

    """
    if not blogpath:
        if blogdir:
            blogpath = f"{blogdir}/{sample}{blogext}"
        else:
            blogpath = f"{sample}{blogext}"

    if not os.path.isfile(blogpath):
        raise ValueError(f"No file found at {blogpath}")

    with open(blogpath, "r") as _f:
        blogstring = _f.read()
    # find binsize
    match = re.search(re.compile("binsize=(?P<b>[0-9]+)"), blogstring)
    if not match:
        raise ValueError("No Binsize information")
    binsize = int(match["b"])
    # find convergence info
    match = re.search(re.compile("converged=(?P<c>True|False)"), blogstring)
    if not match:
        raise ValueError("No Convergence information")
    return binsize

def get_sample_binsize_dict(files,blogpath="",blogdir="cnv_calls",blogext=".binsize.log"):
    """
    Gets binsize of sample based on the LOG file created
    during CNVPytor CNV calling-

    Its a given that the LOG file corresponds to the sample name:
        e.g. sample = "polnareff"
            logfile = "polnareff.binsize.log"
    If no, you can specify the path to the logfile with "blogpath" parameter.


    Parameters
    ----------
    files : iterable, list of files to pair with binsize.
    blogpath : str,
        path to the LOG file of the sample
    blogdir : str, optional
        Directory where the LOG file is found
    blogext : str, optional
        Extension of the file. The default is ".binsize.log".

    """
    return {Path(x).stem:estimate_samplebinsize_logfile(x,blogpath,blogdir,blogext) for x in files}

def convert_signal(signal):
    """
    Converts signal for read depth level
    to read depth statistics for autosomal chromosome

    For example:
        rd_level_600_GC -> rd_stat_600_auto_GC
    """
    binsize,flag = signal.replace("rd_level_","").split("_")
    return f"rd_stat_{binsize}_auto_{flag}"

def get_rd_stats(file,binsize=None,flag=None,lower=4,upper=5):
    """
    Estimates statistic in .hdf file created by CNVPytor
    The statistics returned are:
        binsize : the 100 divisible integer window size
        flag : the appending to the signal
        rd_mean : the mean read depth

    If no binsize is passed (default),
    then the best binsize is estimated using lower and upper parameters

    If no flag is passed (default),
    then the best flag is estimated as the longest flag
    """
    with hpy.File(Path(file),"r") as _f:
        #estimate best binsize
        #according to CNVNATOR parameters
        if binsize is None:
            logging.info("No binsize, estimating from the ratio")
            rd_keys = list()
            _rd_keys = [x for x in _f.keys() if x.startswith("rd_level_")]
            for _rdkey in _rd_keys:
                stat = _f[convert_signal(_rdkey)]
                binsize, flag = _rdkey.replace("rd_level_","").split("_")
                logging.info("For binsize %s, RD mean is %.2f, and RD stdev is %.2f, ratio is %.2f",binsize,stat[4],stat[5],stat[4]/stat[5])
                if lower<=(round(stat[4]/stat[5],2))<=upper:
                    rd_keys.append(_rdkey)

        else:
            logging.info("Passed binsize %f ",binsize)
            rd_keys = [x for x in _f.keys() if x.startswith(f"rd_level_{binsize}")]

        #find the best flag from the above mentioned keys
        #the best flag is the longest...duh....
        if flag is None:
            best_key = max(rd_keys, key=len)
        else:
            best_key = [x for x in rd_keys if x.endswith(flag)]
            if len(best_key)==0:
                logging.error("No signal for %s and %s for file %s",flag,binsize,file)
                raise ValueError(f"No key for specified flag: '{flag}'")
            best_key = best_key[0]

        #find the mean
        mean_rd, _ = np.array(_f[best_key])

    #get the flag and bin of the best key
    binsize,flag = best_key.replace("rd_level_","").split("_")
    return binsize,flag,mean_rd

def get_rd_chroms(file):
    """
    gets chromosomes
    """
    chromosomes = list()
    if isinstance(file,str):
        with hpy.File(file,"r") as _f:
            chromosomes = [x.decode("ascii") for x in _f["rd_chromosomes"]]
    if isinstance(file,list):
        for f in file:
            with hpy.File(f,"r") as _f:
                chromosomes += [x.decode("ascii") for x in _f["rd_chromosomes"]]
    if len(chromosomes)==0:
        logging.error("No RD chromosomes in the given file(s)")
        raise ValueError("There are no RD chromosomes in the given file(s) - pass a path or a list")
    return list(set(chromosomes))

def constrain_chromosomes_by_assembly_report(chromosomes,assembly_report):
    """
    From a given assembly report and a list of chromosomes,
    only selects those chromosomes that are actual chromosomes

    Returns a dictionary of the form:
        chromosome_name : clearer name

    """
    if assembly_report is None:
        logging.info("Not constraining by assembly report")
        return dict(zip(chromosomes,chromosomes))
    assrep = assembly_report.copy()
    try:
        assrep = assrep.loc[assrep['Assigned-Molecule-Location/Type'] == "Chromosome"]
    except KeyError as e:
        logging.error("There was a Key Error - %s",e)
        logging.error("The above key must be in the assembly report")
        raise e


    assrep["chr"] = "Chr " + assrep["Assigned-Molecule"]
    det_col = None
    for col in ["GenBank-Accn", "RefSeq-Accn"]:
        determine = set(assrep[col].tolist()) - set(chromosomes)
        if len(determine) == 0:
            det_col = col

    assrep["CHROMOSOME"] = "chrom " + assrep["Assigned-Molecule"].astype(str)


    return dict(zip(assrep[det_col], assrep["CHROMOSOME"]))

def get_array(file, signal_key):
    """
    gets array from given hdf file
    """
    with hpy.File(file,"r") as _f:
        array = _f.get(signal_key, None)
        if array:
            array = np.array(array)
    return array

def write_read_depth_file(pytor_files,output_file,assembly_report=None,compress="half",
                          binsize_dic=None,flag_dic=None,
                          lower_ratio=4,upper_ratio=5):
    """
    Writes read depth data from a list of CNVPytor-made .hdf files
    to one file found at "output_file"

    Parameters
    ----------
    pytor_files : iterable of strings
        Contains path to the CNVPytor-made .hdf files.
    output_file : str
        Where to write the file.
    assembly_report : str, optional
        Path to assembly report. Assembly report is used to
        constrain chromosomes so that only RD data of
        assembled chromosomes is stored.
    compress : str, optional
        Store in what kind of float datatype.
        "norm" for np.float32, "half" for np.float16, "double" for np.float64
        The default is "half".
    binsize_dic : dict {str:int}
       if passed, retrieves binsize from the dictionary
    flag_dic : dict {str:str}
        if passed, retrieves flagsize from the dictionary
    lower_ratio, upper_ratio:
        Float limits of the ratio of mean Read Depth and its standard deviation.
        Used to find the optimal binsize.

    Returns
    -------
    None.

    """
    #display locals
    #for _k, _v in locals().items():
    #    logging.info("VARIABLE=%s, VALUE=%s",_k,_v)

    rd_chromosomes = get_rd_chroms(pytor_files) # get all chromosomes with read depth
    #constrain them by assembly report
    chrom_map = constrain_chromosomes_by_assembly_report(rd_chromosomes,assembly_report)
    #choose data type to store - default to half
    data_type = COMPRESS_DICT.get(compress,np.float16)
    #write them
    _s = time.perf_counter()
    with hpy.File(output_file, "w") as output:

        for file in pytor_files:
            sample = Path(file).stem
            #check for samples
            if binsize_dic:
                binsize = binsize_dic.get(sample,None)
                if binsize is None:
                    binsize = binsize_dic.get(file,None)
            else:
                binsize = None


            if flag_dic:
                flag = flag_dic.get(sample,None)
                if flag is None:
                    flag = binsize_dic.get(file,None)
            else:
                flag = None

            binsize, flag, mean = get_rd_stats(file,binsize,flag,lower=lower_ratio,upper=upper_ratio)

            # create a group
            sample_group = output.create_group(sample)
            # store information in attributes

            sample_group.attrs["binsize"] = binsize
            sample_group.attrs["flag"] = flag
            sample_group.attrs["mean"] = mean
            logging.error("For sample %s, binsize is %s, flag is %s, and global RD mean is %.2f",sample,binsize,flag,mean)
            for chrom, chrom_name in chrom_map.items():
                logging.info("Writing chrom %s with name %s",chrom,chrom_name)
                # get signal and get an array
                signal_key = f"his_rd_p_{chrom}_{binsize}_{flag}"

                array = get_array(file, signal_key)
                # just in case theres no signal
                if array is None:
                    logging.warning("No data using signal %s",signal_key)
                    continue
                # this "normalizes" the array around diploid copy number
                # meaning that CN of 2 is the null
                centered_array = 2 * array / mean
                # store it as half
                centered_array = centered_array.astype(data_type)
                chrom_dataset = sample_group.create_dataset(
                    chrom, data=centered_array)
                chrom_dataset.attrs["name"] = chrom_name
            logging.error("Written data for sample %s found at %s",sample,file)
            #time.sleep(5)
    _e = time.perf_counter()
    logging.info("R/W done in %.4f seconds",(_e-_s))

def get_cnv_calls(file,chrom,binsize,flag,sample=None):
    """
    Gets the CNV dataframe from Pytor produced file
    """
    if sample is None:
        sample = Path(file).stem
    signal = f"calls_{chrom}_{binsize}_{flag}"
    array = get_array(file,signal)
    if array is None:
        logging.warning("No CNVdata using signal %s in file at %s",signal,file)
        return None
    data = pd.DataFrame(array)
    #I have no idea what the first column is
    data.pop(0)
    #remap the data
    data.columns = ["type","start","end","size","cnv","p_val","p_val_2","p_val_3","p_val_4","Q0","pN","dG"]
    data["type"] = data["type"].map({1:"duplication",-1:"deletion"})
    #broad cast start and end to integer - for nicer processing in BEDtools etc
    data["start"] = data["start"].astype(int)
    data["end"] = data["end"].astype(int)
    #insert relevant info
    data.insert(0,"chrom",chrom)
    data.insert(0,"sample",sample)

    #reshuffle it
    data = data[["chrom","start","end","size","cnv","sample","type","p_val","p_val_2","p_val_3","p_val_4","Q0","pN","dG"]]
    return data

def write_cnv_call_file(pytor_files,output,fmt="csv",assembly_report=None,
                        binsize_dic=None,flag_dic=None,ret=True,
                        lower_ratio=4,upper_ratio=5,):
    """
    Extracts CNV calls from given pytor files and writes it to a file

    pytor_files : iterable of strings
        Contains path to the CNVPytor-made .hdf files.
    output : str
        Where to write the file.
    ret : Bool,
        Return a dataframe.
    assembly_report : str, optional
        Path to assembly report. Assembly report is used to
        constrain chromosomes so that only RD data of
        assembled chromosomes is stored.
    binsize_dic : dict {str:int}
       if passed, retrieves binsize from the dictionary
    flag_dic : dict {str:str}
        if passed, retrieves flagsize from the dictionary
    lower_ratio, upper_ratio:
        Float limits of the ratio of mean Read Depth and its standard deviation.
        Used to find the optimal binsize.
    """

    #check output to determine separator and header
    if output is None:
        logging.info("Will not save file")
        fmt = None
    elif output[-4:] in [".csv",".CSV"]:
        fmt = ".csv"
        output = output[:-4]+fmt
        _sep = ","
        #write header
        with open(output,"w") as _f:
            _f.write(f"{_sep}".join(["chrom","start","end","size","cnv","sample","type","p_val","p_val_2","p_val_3","p_val_4","Q0","pN","dG"])+"\n")

    elif output[-4:] in [".bed",".BED"]:
        fmt = ".bed"
        output = output[:-4]+fmt
        #_header = False

    else:
        if fmt in ["csv","CSV",".csv",".CSV"]:

            fmt = ".csv"
            output = output+fmt
            _sep = ","
            #write header
            with open(output,"w") as _f:
                _f.write(f"{_sep}".join(["chrom","start","end","size","cnv","sample","type","p_val","p_val_2","p_val_3","p_val_4","Q0","pN","dG"])+"\n")
        if fmt in ["bed","BED",".bed",".BED"]:
            fmt = ".bed"
            output = output+fmt
            _sep = "\t"

    logging.info("Saving CNV calls at %s", output)

    #display locals
    #for _k, _v in locals().items():
    #    logging.info("VARIABLE=%s, VALUE=%s",_k,_v)
    #get chromosomes
    rd_chromosomes = get_rd_chroms(pytor_files) # get all chromosomes with read depth
    #constrain them by assembly report
    chrom_map = constrain_chromosomes_by_assembly_report(rd_chromosomes,assembly_report)
    _s = time.perf_counter()
    if ret:
        ret=list()

    for file in pytor_files:
        sample = Path(file).stem
        #check for samples
        if binsize_dic:
            binsize = binsize_dic.get(sample,None)
            if binsize is None:
                binsize = binsize_dic.get(file,None)
        else:
            binsize = None

        if flag_dic:
            flag = flag_dic.get(sample,None)
            if flag is None:
                flag = binsize_dic.get(file,None)
        else:
            flag = None

        binsize, flag, _mean = get_rd_stats(file,binsize,flag,lower=lower_ratio,upper=upper_ratio)

        for chrom, chrom_name in chrom_map.items():
            #get data
            data = get_cnv_calls(file,chrom,binsize,flag,sample)
            if data is None:
                continue
            #save data
            if output:
                data.to_csv(output,mode="a",sep=_sep,index=False,header=False)
                logging.info("Appended CNV data for %s on chromosome %s to file %s",sample,chrom,output)
            if ret is not None:
                ret.append(data)
                logging.info("Appended CNV data for %s on chromosome %s",sample,chrom)
            #
    _e = time.perf_counter()
    logging.info("R/W done in %.4f seconds",(_e-_s))
    if ret:
        return pd.concat(ret,axis=0,ignore_index=True)
    
def get_samples(rd_file):
    """Gets sample names found in RD Files"""
    try:   
        with hpy.File(rd_file,"r") as _f:
            samples = set(_f.keys())
        return samples
    except OSError:
        return None
    
def read_all_rdchroms(rd_file):
    """Gets all chromosomes in a given file"""
    try:
        with hpy.File(rd_file, "r") as f:
            chroms = [list(f[x].keys()) for x in f.keys()]
        return sorted(list(set([item for sublist in chroms for item in sublist])))
    except OSError:
        return None