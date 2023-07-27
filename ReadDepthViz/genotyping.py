# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 17:27:22 2022

This script "genotypes" (provides average copy number) for
a list of given regions.

===============================================================================
To genotype a list of given regions, use function "genotype".

===============================================================================
    PARAMETERS
===============================================================================
hdf_file:
    Pointer to a HDF file. The HDF file must be structured in the following manner.
    The topmost keys must correspond to samples. The topmost keys must
    have "binsize" attribute.
    Every topmost key must have keys corresponding to chromosomes.
    Every inner keys must have an array corresponding to
    a diploid centered read depth array.

    This hdf_file is the output of a previous module that
    concatenates outputs of CNVPytor into one smaller file.

cnv_df:
    pandas DataFrame. Must have columns ["chrom","start","end"].
    This is our list of regions that will be genotyped.

output_file:
    Where data will be saved.

max_cores (OPTIONAL):
    Number of cores for parallel processing. The default is 1.

samples (OPTIONAL):
    All keys of the sample. The default is None (iterates through all samples)

===============================================================================
    OUTPUT
===============================================================================
The output is a pandas DataFrame that is saved to a .csv file:

    - The indices of this dataframe are regions joined by a comma.
    - The columns of this dataframe are samples found in hdf_file.
    - The values of this dataframe are the average diploid copy number.


===============================================================================
    CURRENT PERFORMANCE
===============================================================================
Genotype of 42 x 15k size done:
    563.8423 seconds (single core)
    181.4731 seconds (12 cores)

For one 15k size the result is roughly 18 seconds.

===============================================================================
    INFORMATION ABOUT HOW CNVPYTOR GENOTYPES A VALUE,
    HOW IT CALCULATES MEAN, AND HOW I'VE IMPROVED ITS PERFORMANCE
===============================================================================

CNVPytor genotypes a region by taking its histogram signal, which is an array,
and takes a "modified" average. For example, suppose you have
a chromosome of size 1 500 divided into bins of 400.
1 500 // 400 = 3 meaning our array will have 3 + 1 elements:
    0 0-400 RD_BIN_0
    1 400-800 RD_BIN_1
    2 800-1200 RD_BIN_2
    3 1200-1600 RD_BIN_3

To get average read depth, you would sum all values in the third column,
then divide by the total length (assuming uniform RD distribution as the null).

The first problem lies in edge cases, for example - region of 200 - 1000.

Example 1 - a region of 200 to 1000:

    Start is 200 and the end is 1000, whilst binsize is 400.
        200 // 400 = 0 so the starting index will be 0
        1000 // 400 = 2 so the ending index will be 2
    We should then slice the array [0:2+1] (python is exclusive to end),
    but then we get the same thing as if we used 0 as a start, and 800 as
    an end.
        0 // 400 = 0 so the starting index will be 0
        1200 // 400 = 2 so the starting index will be 2+1
    The CNVPytor solves this by truncating the starting and ending indexes.
    For example, 200 is halfway into the first bin - 200/400=0.5
    So if the first bin has RD of 10, the result will be 5 (10 * 0.5).
    1000 is halfway into the third bin, so similar thing is done.

    Fractional formula:
        frac_left = 1 - (start - idx_left * binsize)/binsize
        frac_right = (end - idx_right * binsize)/binsize

    Proof:
        start = 200
        idx_left = 0

        end = 1000
        idx_right = 2
        binsize = 400

        frac_left = 1 - (200 - 0 * 400)/400
        frac_left = 1 - 0.5 = 0.5

        frac_right = (1000 - 2 * 400)/400
        frac_right = (1000 - 800)/400
        frac_right = (200)/400 = 0.5

    We solve this by the following code:

        array = entire_chromosome_array[idx_left:idx_right+1]
        array[0] = array[0] * frac_left
        array[-1] = array[-1] * frac_right

    Now we can call .sum() on the array

The second problem lies with what happens after the end of the chromosome.

Example 2 - a region of 900 to 2000
    Start is 900 and the end is 2000, whilst binsize is 400
        900 // 400 = 2 so the starting index will be 0
        2000 // 400 = 5 so the ending index will be 5

        frac_left = 1 - (900 - 2 * 400)/400 = 1 - 0.25 = 0.75
        frac_right = (2000 - 5 * 400)/400 = 0

    The first problem is that the ending index is wrong.
    But numpy's array slicing will simply ignore excess and return the entire
    array back.

    Array of size 4 sliced from 2 to 6 (5+1) will return
    array sliced from 2 to its end (3+1). The problem then is that
    frac_right defaults to 0 and will simply ignore the last index.

    To solve this, we implement the following:
        if id_right >= len_array:
            end = (len_array)*bin_size-1e-9
    Proof:
        if 5 > 4:
            end = (4) * 400 - 1e-9
            end = 1600 - 1e-9
            end = 1599.9...

    Now with our new end we have:
        900 // 400 = 2 so the starting index will be 0
        1599.9... // 400 = 3 so the ending index will be 3

        frac_left = 1 - (900 - 2 * 400)/400 = 1 - 0.25 = 0.75
        frac_right = (1599.9... - 3 * 400)/400 = 399.9.../400 = 0.9...

    Without 1e-9 part, the frac_right would simply be 0,
    as (idx * binsize)/binsize - (idx * binsize) // binsize is always 0
    if idx is a multiple of binsize.


The next problem refers to average genotype value of a region, or the average copy number (CN) of the region.
    Average CN is calculated as a mean of RD values per bins spanning the region.
    In the case above (binsize=400):
        0 0-400 RD_BIN_0
        1 400-800 RD_BIN_1
        2 800-1200 RD_BIN_2
        3 1200-1600 RD_BIN_3

    The average genotype value of the region spanned by 0 to 1600
    is (RD_BIN_0 + RD_BIN_1 + RD_BIN_2 + RD_BIN_3)/(4).

    But in the case of the overhang or underhang, for example
    for region spanned by 200 to 1400, we still have the same sum.
    But it feels wrong to divide it by 4, since we're only taking the half
    of the first and last bin. So we divide it by 3.
    The full formula for the size:
        size = len(region) - 2 + frac_left + frac_right
    For example of 200 to 1400, we have an array of 4 size.
    But since we have an overhang in the first we subtract 2 from the size,
    and then we add frac_left and frac_right, which are both equal to 0.5.
    The result is:
        size = 4 - 2 + 0.5 + 0.5 = 3

    NOTE: This is functionally identical to division of size by binsize.
    For example:
        (1400 - 200) / 400 = 1200 / 400 = 3


@author: ivanp
"""
import os
import time
import logging
from multiprocessing import Pool
from functools import partial
import h5py
import pandas as pd
import numpy as np


def infer_samples(hdf_file):
    """
    Returns all "keys" (top folders)
    of a HDF file found at "hdf_file"
    """
    with h5py.File(hdf_file) as _f:
        samples = list(_f.keys())
    return samples

def _index_on_locus(df, inplace=False):
    """
    Creates a index based on locus.
    Df in general must have columns
    "chrom","start","end" and a string type locus is created in
    the form of {chrom},{start},{end}.

    If inplace is passed, data will be reindexed inplace.
    """
    if inplace:
        cnvdf = df
    else:
        cnvdf = df.copy()

    cnvdf["region"] = cnvdf["chrom"].astype(
        str)+","+cnvdf["start"].astype(str)+","+cnvdf["end"].astype(str)
    cnvdf.set_index("region", inplace=True)
    if inplace:
        return None
    return cnvdf

def index_on_locus(df, cols = ["chrom","start","end"],inplace=False,drop=False):
    """
    Creates a index based on locus.
    Df in general must have columns
    "chrom","start","end" and a string type locus is created in
    the form of {chrom},{start},{end}.

    If inplace is passed, data will be reindexed inplace.
    """
    if inplace:
        cnvdf = df
    else:
        cnvdf = df.copy()
    #cols = []
    #if not cols:
    #    cols = ["chrom","start","end"]

    cnvdf["region"] = cnvdf[cols[0]].astype(
        str)+","+cnvdf[cols[1]].astype(str)+","+cnvdf[cols[2]].astype(str)
    cnvdf.set_index("region", inplace=True)
    if inplace:
        return None
    if drop:
        cnvdf.drop(cols,axis=1,inplace=True)
    return cnvdf

def get_rd_region(sample, hdf_file, chromosome, start, end, index=False):
    """
    Get the read depth region

    CAVEAT: When 'end' parameter exceeds
    the length of array multiplied by binsize, the array will be still
    returned correctly.
    Example
    -------
    We have a chromosome of size 1 500 divided into bins of 400

    0 : 000 - 400 : ...
    1 : 400 - 800 : ...
    2 : 800 - 1200 : ...
    3 : 1200 - 1600 : ...

    The chromosome is stored in 4 x 1 array
    If we try to access a region with end point of 1 700,
    we'll get 1 700 // 400 = 4.

    However, passing array[:3+1] and array[:4+1]
    will return the same lengths of array.
    The problem will be downstream - when genotyping.

    Since its highly unlikely that the passed end falls neatly
    into histogram divisions, read depth values are multiplied what
    percentage the region is in the histogram.

    In the above example, size 1 500 is found in bin number 3
    [1 200 - 1 600]. But its only 3/4 in the bin, so if its read-depth
    were 4 - it would be weighted to be 3.


    The full calculation for it is:
        frac_right = (end - (end//binsize)*binsize) / binsize

    In the above example of 1500,
        frac_right  = (1500 - (1500//400)*400)/400
                    = (1500 - 3 *400) / 400
                    = 300 / 400 = 3/4
        since we're going to 1500//400 = 3
        we get the value of bin[3] * 3/4

    If we had a chromosome of 1 600, it's would be:
        frac_right  = (1 600 - (1 600//400)*400)/400
                    = (1 600 - 400 * 400 )/400
                    = 0
        This would work since we're going in the ID 4 (1600//400=4).
        So the fraction of the next ID would be 0.

    But what happens when we overshoot, and use 1 700 as end point?
    Since 1 700 // 400 = 4, we'll get an entire array
    But, frac_right will be a problem:
        =(1 700 - (1700//400)*400)/400
        =(1 700 - 1600)/400
        =1/4
    The problem is that since there's next data, our
    i_right would be that 1 700//400=4, so we still have an overhang
    To solve that problem, we default on the size:
        if a given end point is greater or equal to than length of
        the array multiplied by the binsize, then the end point will simply
        be the length of the array multiplied by the binsize.
    In the case above 1 700 > 4 * 400,
    so our end size is 1600. But the problem will still persist!
    In our case, since our end_size is 1600, our fract right will go to 0.

        = (1600 - (1600//400)*400)/400
        = (1600 - 4 * 400)/400
        = 0 / 400
    To rectify this, we simply make the following calculation

    new_end = len_array * binsize - 1e-9

    The problem is that longdiv is a periodic operation, which works
    for every instance, except the end point. To solve the end point,
    we provide a hack! in the form of removing an infinitesimal number.
    So our new end will be 1599.9999. In the formula above, this will
    cause our frac_right to be 0.99999, and since we are dealing with
    half types, this will round to one.



    """
    with h5py.File(hdf_file, "r") as _f:

        bin_size = int(_f[sample].attrs["binsize"])

        id_left = int(start // bin_size)
        id_right = int(end // bin_size)

        len_array = len(_f[sample][chromosome])
        #region = _f[sample][chromosome][id_left:id_right+2]
        region = _f[sample][chromosome][id_left:id_right+1]
        #correct would be id_right+1
        # compensate when the end of array is achieved
        if id_right >= len_array:
            # remove a small piece of end to work with longdiv
            # arithmetic
            end = (len_array)*bin_size-1e-9

    if not index:
        return region, bin_size, end
    index = np.arange(start=id_left*bin_size,
                      stop=(id_right+1)*bin_size, step=bin_size)
    region = pd.Series(region, index=index)
    return region

def genotype_region(sample, hdf_file, chromosome, start, end):
    """
    Get the average genotype value of a given region
    """

    # get region
    region, bin_size, end = get_rd_region(
        sample, hdf_file, chromosome, start, end)

    # get the fractions
    frac_left = 1 - (start - int(start // bin_size) * bin_size)/bin_size
    frac_right = (end - int(end // bin_size) * bin_size)/bin_size

    # modify by the fractions
    region[0] = region[0] * frac_left
    region[-1] = region[-1] * frac_right

    # get the size - compensated by fractions
    size_of_region = len(region) - 2 + frac_left + frac_right

    # genotype
    return np.half(region.sum()/size_of_region)

def _apply_genotyping(row, sample, hdf_file):
    """
    Function for apply genotype optimized
    for application on axis=1 dataframe.
    The dataframe must have columns "chrom","start","end".

    The "sample" is the foremost "key" (folder) of
    the HDF file found at "hdf_file"
    """
    chromosome = row["chrom"]
    start = row["start"]
    end = row["end"]

    return genotype_region(sample, hdf_file, chromosome, start, end)

def genotype_sample(cnv_df, hdf_file, sample, index_locus=False):
    """
    Genotype a "sample", a key of HDF file
    found at "hdf_file" (the topmost level),
    using values found in pandas DataFrame "cnv_df"

    Parameters
    ----------
    cnv_df : pandas DataFrame
        pandas DataFrame. Must have columns ["chrom","start","end"]
    hdf_file : str
        Pointer to a HDF file.
    sample : str
        Key of a HDF file (top level folder).
    index_locus : Boolean, optional
        Whether to reindex using locus (see "index_on_locus" function).

    Returns
    -------
    pandas Series
        The name of the series is "sample" parameter, and its values
        are average read depths of a given region defined in "cnv_df".
        The average read depths are normalized around 2.

    """

    if index_locus:
        cnv_df = index_on_locus(cnv_df, False)

    _s = time.perf_counter()
    cnv_df[sample] = cnv_df.apply(_apply_genotyping, axis=1,
                                  sample=sample, hdf_file=hdf_file)
    _e = time.perf_counter()
    logging.info("Genotype of %s size for sample %s done in %.4f seconds",len(cnv_df),sample,_e - _s)
    return cnv_df[sample]

def _genotype_singlecore(func, sample_list):
    """
    Iterates one by one
    """
    result_df = list()
    for sample in sample_list:
        result_df.append(func(sample))
    dataframe = pd.concat(result_df, axis=1)
    return dataframe

def _genotype_multicore(func, sample_list, max_cores):
    """
    Uses parallel processing to genotype
    """
    with Pool(max_cores) as pool:
        result = pool.imap_unordered(func, sample_list)
        result_df = [r for r in result]
    dataframe = pd.concat(result_df, axis=1)
    return dataframe


def genotype(hdf_file, cnv_data, output_file, max_cores=1, samples=None):
    """
    Finds average diploid copy number of a region
    defined in pandas dataframe "cnv_data" for all samples
    stored in hdf_file.



    Parameters
    ----------
    hdf_file : str
        Pointer to a HDF file.
    cnv_df : pandas DataFrame
        pandas DataFrame. Must have columns ["chrom","start","end"]
    output_file : str
        Where data will be saved.
    max_cores : int, optional
        Number of cores for parallel processing. The default is 1.
    samples : list, optional
        All keys of the sample. The default is None (iterates through all samples)

    Returns
    -------
    pandas DataFrame
        The indices of this dataframe are regions joined by a comma.
        The columns of this dataframe are samples found in hdf_file.
        The values of this dataframe are the average diploid copy number.

    """
    # set maximum cores
    max_cores = min(max_cores, os.cpu_count())
    logging.info("STARTing genotyping using %s cores",max_cores)

    #display locals
    for _k, _v in locals().items():
        if isinstance(_v,pd.DataFrame):
            _v = f"pandas DataFrame of {_v.shape}"
        logging.info("VARIABLE=%s, VALUE=%s",_k,_v)

    # infer samples if not provided
    if not samples:
        samples = infer_samples(hdf_file)

    # reindex on loci
    cnv_data = index_on_locus(cnv_data)

    # freeze the function for faster processing
    genotyping_func = partial(genotype_sample, cnv_data, hdf_file)

    # start the function
    _s = time.perf_counter()
    if max_cores == 1:
        dataframe = _genotype_singlecore(genotyping_func, samples)
    else:
        dataframe = _genotype_multicore(genotyping_func, samples, max_cores)
    _e = time.perf_counter()
    logging.info("Genotype of %s size for %s samples DONE in %.4f seconds",len(cnv_data),len(samples),_e - _s)
    if output_file:
        dataframe.to_csv(output_file)
    return dataframe
