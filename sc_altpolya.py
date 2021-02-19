#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:39:47 2021

@author: beth
"""

#Load required packages
import os, subprocess, h5py, re
import pandas as pd
import numpy as np
from sklearn.mixture import BayesianGaussianMixture
from collections import defaultdict
import matplotlib.pyplot as plt
#also requires samtools

#Add directories to path so that samtools is in path
os.environ['PATH'] = '/opt/anaconda3/bin:/opt/anaconda3/condabin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin'


def get_top_genes_and_barcodes_list(file_name, n_genes=100):
    '''
    Return the top expressed genes and filtered barcodes list from 10X .h5 file

    Parameters
    ----------
    file_name : .h5 filename
        a 10x filtered .h5 file
        
    n_genes : integer
        number of top genes by counts to return. The default is 100.

    Returns
    -------
    top_genes : data frame
        the top n expressed genes and thier counts 
    barcodes : list
        the filtered barcodes list from the .h5 input file
    '''

    dset = h5py.File(file_name, 'r')
    dset = dset['matrix']
    indices = dset['indices']
    data = dset['data']
    features = dset['features']
    ids = features['id']
    names = features['name']
    barcodes = pd.Series(dset['barcodes'][()]).astype('str').str.replace('b','').str.replace("'",'')
    
    counts = pd.DataFrame({'count':data,'gene':indices})
    counts = counts.groupby('gene').sum().reset_index()
    
    ensembl = pd.DataFrame({'ids':ids,'names':names})
    ensembl = ensembl.apply(lambda x: x.astype('str').str.replace('b','')).reset_index()
    
    temp = pd.merge(left=counts,right=ensembl,how="left",left_on='gene',right_on='index')
    top_genes = temp.sort_values('count',ascending=False).reset_index()[['count','ids','names']].head(n_genes)
    
    return top_genes, barcodes


def gtf_top_genes(top_genes, gtf_file):
    '''
    Gets chrom, strand and exon information for each gene in top_genes 

    Parameters
    ----------
    top_genes : dataframe 
        top genes dataframe returned from get_top_genes_and_barcodes_list
        or any dataframe with genes of interest in column 'ids'
    gtf_file : .gtf filename
        ensembl gtf file matching genome build used in the 10X bam file

    Returns
    -------
    df : dataframe
        the top genes dataframe with column added
            exons: set of all genomic exons positions for the gene
            chrom: chromosome of gene
            strand: genomic strand of gene

    '''
    exon_list=[]
    chrom_list=[]
    strand_list=[]
    for i in top_genes['ids']:
        #print(i)
        #ensembl = !grep $i $gtf_file | cut -f1,3,4,5,7
        i=i.strip("'")
        grep_gene = subprocess.Popen(('grep', i, gtf_file), stdout=subprocess.PIPE)
        ensembl = subprocess.check_output(('cut','-f1,3,4,5,7'), stdin=grep_gene.stdout)
        ensembl = ensembl.decode("utf-8") 
        ensembl = ensembl.split('\n')
        
        try:
            exons = list(map(lambda x: ('\texon') in x, ensembl))
            ensembl = np.array(ensembl)
            ensembl = ensembl[exons]
            chrom = re.findall('(chr[^\t]*)',ensembl[0])[0]
            temp = list(map(lambda x: set(x), map(lambda x: range(int(re.findall('\t([0-9]*)',x)[1]),\
                                                      int(re.findall('\t([0-9]*)',x)[2])), ensembl)))
            exons = list(set([i for i in temp for i in i]))
            exon_list.append(exons)
            chrom_list.append(chrom)
            strand = ensembl[0][-1]
            strand_list.append(strand)
        
        #in case gene id is not in gtf
        except: 
            exon_list.append(0)
            chrom_list.append(0)
            strand_list.append(0)
            
    df  = pd.DataFrame({'gene_id':top_genes['ids'], 'name':top_genes['names'], 'chrom':chrom_list, 'exons':exon_list, 'strand':strand_list})
    df = df[df['chrom'] != 0] 
    return df


def bam_top_genes(df, filtered_barcodes, bam_file): 
    '''
    Considering gene strand, finds all reads that start between the first exon position and the last. 
    Finds read stop positions and only retains furthest 3' read-stop for each UMI group.
    Allows for read_stops to occur up to 300kbp after last exon position.

    Parameters
    ----------
    df : dataframe
        dataframe from gtf_top_genes with exon set, chrom and strand columns
    filtered_barcodes : list
        barcodes list from get_top_genes_and_barcodes_list
    bam_file : .bam filename
        10x bam file (from Cellranger pipeline)

    Returns
    -------
    df : dataframe
        input dataframe updated with columns
            read_stops: read stop positions
            barcodes: corresponding barcodes 

    '''
    
    df = df.reset_index()
    
    read_stops_list_of_list=[]
    barcodes_list_of_list = []
    
    for i in range(len(df)):
        
        if df.loc[i,'strand'] == "+":
            #print(i)
            #Get the minumum and maximum exon position
            chrom = df.loc[i,'chrom']
            start = min(df.loc[i, 'exons'])
            stop = max(df.loc[i, 'exons'])
            
            #find reads starting between these positions
            view = subprocess.run(['samtools', 'view', bam_file, chrom+":"+str(start)+"-"+str(stop)],stdout=subprocess.PIPE)
            view = view.stdout.decode("utf-8") 
            view = view.split('\n')
            view = view[:-1]

            #print('parsing view')
            #parse the samtools view for corrected UMI, barcode, strand and cigar info
            umi = list(map(lambda x: re.findall('UB:Z:([A-Z]+)', x),view))
            barcodes = list(map(lambda x: re.findall('CB:Z:([A-Z]+-[0-9]+)', x),view))
            read_strand = list(map(lambda x: x.split('\t')[1],view))
            read_starts = list(map(lambda x: x.split('\t')[3],view))
            cigars = list(map(lambda x: x.split('\t')[5],view))
            #find the read stop position using cigar string
            cigars = list(map(lambda x: re.findall('[0-9]+', x), cigars))
            cigars = list(map(lambda x: sum([int(s) for s in x]), cigars))
            read_stops = np.array(list(map(int,read_starts))) + np.array(cigars)
            
            #print('making test df')
            #make a temporary dataframe with parsed information for each read
            test = pd.DataFrame({'umi':umi,'read_starts':read_starts,'read_strand':read_strand,'barcodes':barcodes,'read_stops':read_stops})
            #only retain positive strand reads here
            test = test[test['read_strand']=='0']
            #only retain read_stops within 300k of last exon position
            test = test[(test['read_stops'] < stop + 300000)]
            #remove reads missing corrected umi or barcode
            test = test[list(map(lambda x: len(x) != 0, test['umi']))]
            test = test[list(map(lambda x: len(x) != 0, test['barcodes']))]
            #reformat/unlist the umi and barcode
            test['umi'] = list(map(lambda x: x[0], test['umi']))
            test['barcodes'] = list(map(lambda x: x[0], test['barcodes']))
            #only retain barcodes from filtered barcodes list
            test = test[test['barcodes'].isin(filtered_barcodes)]
            #sort by umi and read stops to retain only furthest 3' read stop
            test = test.sort_values(by=['umi','read_stops'],ascending=[True,False])
            test = test.drop_duplicates('umi', keep='first')
            
            #print('appending_list')
            read_stops_list_of_list.append(list(test['read_stops']))
            barcodes_list_of_list.append(list(test['barcodes']))
            
        if df.loc[i,'strand'] == "-":
            print(i)
            #Get the minumum and maximum (actual start) exon position 
            #subtract 300k from the minumum to include 300k downstream of gene
            chrom = df.loc[i,'chrom']
            actual_start = min(df.loc[i, 'exons']) 
            start = min(df.loc[i, 'exons']) - 300000
            stop = max(df.loc[i, 'exons'])
            
            #find reads starting between these positions
            view = subprocess.run(['samtools', 'view', bam_file, chrom+":"+str(start)+"-"+str(stop)],stdout=subprocess.PIPE)
            view = view.stdout.decode("utf-8") 
            view = view.split('\n')
            view = view[:-1]
            
            #print('parsing view')
            #parse the samtools view for corrected UMI, barcode, strand and cigar info
            umi = list(map(lambda x: re.findall('UB:Z:([A-Z]+)', x),view))
            barcodes = list(map(lambda x: re.findall('CB:Z:([A-Z]+-[0-9]+)', x),view))
            read_strand = list(map(lambda x: x.split('\t')[1],view))
            read_starts = list(map(lambda x: x.split('\t')[3],view))
            cigars = list(map(lambda x: x.split('\t')[5],view))
            #find the read stop position using cigar string
            #note for negative strand read stop is really the read start
            #and read start is really the read stop
            cigars = list(map(lambda x: re.findall('[0-9]+', x), cigars))
            cigars = list(map(lambda x: sum([int(s) for s in x]), cigars))
            read_stops = np.array(list(map(int,read_starts))) + np.array(cigars)

            #print('making test df')
            #make a temporary dataframe with parsed information for each read
            test = pd.DataFrame({'umi':umi,'read_starts':read_starts,'read_strand':read_strand,'barcodes':barcodes,'read_stops':read_stops})
            #only retain negative strand reads here
            test = test[test['read_strand']=='16']
            #only retain read_stops (technically read starts) within gene exons
            test = test[(test['read_stops'] > actual_start) & (test['read_stops'] < stop)]
            #remove reads missing corrected umi or barcode
            test = test[list(map(lambda x: len(x) != 0, test['umi']))]
            test = test[list(map(lambda x: len(x) != 0, test['barcodes']))]
            #reformat/unlist the umi and barcode
            test['umi'] = list(map(lambda x: x[0], test['umi']))
            test['barcodes'] = list(map(lambda x: x[0], test['barcodes']))
            #only retain barcodes from filtered barcodes list
            test = test[test['barcodes'].isin(filtered_barcodes)]
            #sort by umi and read stops to retain only furthest 3' read start (technically read stop)
            test = test.sort_values(by=['umi','read_starts'], ascending=[True,True])
            test = test.drop_duplicates('umi', keep='first')
            #print('appending_list')

            read_stops_list_of_list.append(list(test['read_starts']))
            barcodes_list_of_list.append(list(test['barcodes']))
        
    df['barcodes'] = barcodes_list_of_list
    df['read_stops'] = read_stops_list_of_list
    return df


def count_and_summarize_peaks(df2, num_peaks_exon = 3, num_peaks_other = 3):
    '''
    Uses a gaussian mixture model to find peaks up to the number specified of read_stops for 
        exons: read_stops occuring inside exons. The introns are removed and exons juxtaposed, so that
        the peaks are detected along the coding sequence of the gene.
        other: read_stops occuring outside exons. This includes introns, and genomic sequence outside the gene, up to 300kbp downstream.
    
    Parameters
    ----------
    df2 : dataframe
        dataframe after bam_top_genes() function with exons, strand, readstops, barcodes
    num_peaks_exon : int, optional
        maximum number of peaks GMM model looks for in exons. The default is 3.
    num_peaks_other : int, optional
        maximum number of peaks GMM model looks for in others. The default is 3.

    Returns
    -------
    final_counts : dataframe
        a counts table of each alt transcript for each cell barcode
    final_summary : dataframe
          a summary table with count, mean, and std information for each alt transcript found
          the mean is relative to the first exon position found from the gtf file earlier

    '''
    
    final_counts = pd.DataFrame()
    final_summary = pd.DataFrame()
    df = df2.copy()
    df = df.reset_index()
    
    for i in range(len(df)):
        #print(i)
        if df.loc[i,'strand'] == "+":
            start = min(df.loc[i,'exons'])
            exons = list(map(lambda x: x - start, df.loc[i,'exons']))
            read_stops = list(map(lambda x: int(x) - start, df.loc[i,'read_stops']))
        if df.loc[i,'strand'] == "-":
            start = max(df.loc[i,'exons'])
            exons = list(map(lambda x: -1*(x - start), df.loc[i,'exons']))
            exons = exons[::-1]
            read_stops = list(map(lambda x: -1*(int(x) - start), df.loc[i,'read_stops']))
       
        read_stops_dict = defaultdict(int)
        barcodes_dict = defaultdict(list)
        read_barcode_zip = zip(df.loc[i,'barcodes'],read_stops)
        for b in read_barcode_zip:
            read_stops_dict[b[1]] += 1
            barcodes_dict[b[1]].append(b[0])
                   
        
        #exons model
    
        read_stops_exons=defaultdict(int)
        barcodes_exons=defaultdict(list)
        for key in exons:
            if read_stops_dict.get(key):
                read_stops_exons[key] = read_stops_dict.get(key)
                barcodes_exons[key] = barcodes_dict.get(key)
                
        for_bgm = [x for y in [[exons.index(i)]*n for i, n in read_stops_exons.items()] for x in y]
        for_bgm = np.array(for_bgm).reshape([-1,1])
        result = BayesianGaussianMixture(num_peaks_exon,max_iter=2000).fit(for_bgm).predict(for_bgm)
        
        barcodes = [np.array(i).flatten() for i in barcodes_exons.values()]
        barcodes = np.concatenate(barcodes)
                
        exons_results = pd.DataFrame({'result':result,'pos':for_bgm.flatten(),'barcode':barcodes})
        exons_results['result'] = exons_results['result'].map(lambda x: df.loc[i,'name'] + '_exons_' + str(x))
        exons_results['pos'] = exons_results['pos'].map(lambda x: exons[x])
        
        #other model
        
        try:
            other_read_stops=defaultdict(int)
            other_barcodes=defaultdict(list)
            for key in read_stops_dict.keys():
                if key not in exons:
                    other_read_stops[key] = read_stops_dict.get(key)
                    other_barcodes[key] = barcodes_dict.get(key)

            for_bgm = [x for y in [[i]*n for i, n in other_read_stops.items()] for x in y]
            for_bgm = np.array(for_bgm).reshape([-1,1])
            result = BayesianGaussianMixture(num_peaks_other,max_iter=2000).fit(for_bgm).predict(for_bgm)

            barcodes = [np.array(i).flatten() for i in other_barcodes.values()]
            barcodes = np.concatenate(barcodes)

            other_results =  pd.DataFrame({'result':result,'pos':for_bgm.flatten(),'barcode':barcodes})
            other_results['result'] = other_results['result'].map(lambda x: df.loc[i,'name'] + '_others_' + str(x))
        
            gene_final = pd.concat([exons_results,other_results],axis=0,ignore_index=False)
            
        except ValueError:
            gene_final = exons_results
        
        gene_summary = gene_final.groupby(['result']).agg({'pos':['mean','std','count']})
        
        gene_counts = gene_final.groupby(['barcode','result']).count().rename(columns={'pos':'count'})
        
        final_summary = pd.concat([final_summary, gene_summary],axis=0,ignore_index=False)
        
        final_counts = pd.concat([final_counts,gene_counts],axis=0,ignore_index=False)

    
    return final_counts, final_summary


def select_alt_transcripts(final_counts, final_summary, min_per_gene = 2, max_within_gene_correlation = 1, count_greater_than_std = True):
    '''
    Filter the counts and summary dataframes. Can remove peaks/alt transcripts if less than 2 were found for a gene (a single alt transcript will be highly correlated with parent gene counts), 
    remove alt transcript if counts was less than standard deviation (wide sparse peaks are more likely to be noise) and 
    remove alt transcript if minimmum intra gene correlation was above a certain level (not informative).
    
    Parameters
    ----------
    final_counts : dataframe
        counts dataframe from count_and_summarize_peaks()
    final_summary : dataframe
        summary dataframe from count_and_summarize_peaks()
    min_per_gene : integer, optional
        Remove genes if less than this number of alt transcripts were found. The default is 2.
    max_within_gene_correlation : float between -1 and 1, optional
        Remove genes if minimum intra-gene correlation is higher than this. The default is 1.
    count_greater_than_std : boolean, optional
        Should alt transcripts with count greater than standard deviation be removed? The default is True.

    Returns
    -------
    final_counts : dataframe
        Filtered counts table
    final_summary : dataframe
        Filtered summary table

    '''
    
    #remove by counts greater than std               
    if count_greater_than_std == True:
        final_summary = final_summary[final_summary['pos']['count']>final_summary['pos']['std']]
    
    #remove by number of peaks per gene
    genes = final_summary.index.str.extract("\'(.*)\'")
    keep = genes[0].value_counts() > (min_per_gene - 1)
    keep = keep[keep == True].index
    
    final_summary = final_summary.droplevel(level=0,axis=1)
    final_summary = final_summary.reset_index()
    final_summary['name'] = genes[0]
    final_summary = final_summary[final_summary['name'].isin(keep)]
    
    final_counts = final_counts.reset_index()
    final_counts = final_counts[final_counts['result'].isin(final_summary['result'])]
    
    #remove by intra gene correlation
    min_corr = min_correlations_by_gene(final_counts, final_summary)
    keep = min_corr[min_corr['min_within_gene_correlation'] < max_within_gene_correlation]
    keep = keep['name']
    final_summary = final_summary[final_summary['name'].isin(keep)]
    final_counts = final_counts[final_counts['result'].isin(final_summary['result'])]
    
    return final_counts, final_summary
        
        
def peak_correlations_for_gene(base_gene_id, counts):
    '''
    Return correlation matrix for all alternate transcripts for a given gene.

    Parameters
    ----------
    base_gene_id : string
        The gene name/symbol for the gene of interest
    counts : dataframe
        counts matrix from seperate_into_peaks() or filtered counts from select_alt_transcripts()

    Returns
    -------
    correlation matrix

    '''
    counts = counts.loc[counts['result'].str.contains("'"+base_gene_id+"'"),].pivot(index = 'barcode', columns = "result")
    return counts.corr() 


def min_correlations_by_gene(counts, summaries):
    '''
    Find the minimum intra-gene correlation of alternate transcripts
    If all the alternate transcripts are highly correlated, there is less information 

    Parameters
    ----------
    counts : dataframe
        counts dataframe from seperate_into_peaks() either filtered or not filtered
    summaries : dataframe
        summaries dataframe from seperate_into_peaks() either filtered or not filtered


    Returns
    -------
    df : dataframe
        A dataframe with a column 'name' for gene name and 'min_within_gene_correlation' for the minimum intragene correlation 
        of alt transcripts for that gene

    '''
    min_correlations = []
    genes = []
    for gene in summaries['result'].str.extract("(\'.*\')")[0].unique():
        c = counts.loc[counts['result'].str.contains(gene),].pivot(index = 'barcode', columns = "result")
        c = c.corr().min()
        min_correlations.append(c)
        min_corr = list(map(min, min_correlations))
        genes.append(gene.strip("'"))
        df = pd.DataFrame({'name':genes, 'min_within_gene_correlation':min_corr})
    return df

def graph_exons(gene_name, df, num_peaks = 3, num_bins = 100):
    '''
    Produce 3 histogram graphs for read_stops inside exons before and after assignment by GMM algorithm. Note that the GMM model is non-deterministic 
        so sometimes the results will difer. Running this function repeatedly with a different num_peaks may give an idea of the best num_peaks parameter for \ the count_and_summarize_peaks() function. 

    Parameters
    ----------
    gene_name : string
        gene id/symbol for gene of interest
    df : dataframe
        the dataframe returned after bam_top_genes()
    num_peaks : integer, optional
        maximum number of peaks for the GMM model. The default is 3.
    num_bins : integer, optional
        number of bins for histograms. The default is 100.

    Returns
    -------
    3 graphs
        1- read_stops inside exons coding sequence histogram
        2- read_stops inside exons genomic location histogram
        3- num_peaks histograms for each of the peaks found 

    '''
    #df=df2.copy()
    df=df[df['name']=="'"+gene_name+"'"]
    df=df.reset_index()
    
    if df.loc[0,'strand'] == "+":
        start = min(df.loc[0,'exons'])
        exons = list(map(lambda x: x - start, df.loc[0,'exons']))
        read_stops = list(map(lambda x: int(x) - start, df.loc[0,'read_stops']))
    if df.loc[0,'strand'] == "-":
        start = max(df.loc[0,'exons'])
        exons = list(map(lambda x: -1*(x - start), df.loc[0,'exons']))
        exons = exons[::-1]
        read_stops = list(map(lambda x: -1*(int(x) - start), df.loc[0,'read_stops']))
    
    read_stops_dict = defaultdict(int)
    
    for stop in read_stops:
        read_stops_dict[stop] += 1
        
    read_stops_exons=defaultdict(int)
    for key in exons:
        if read_stops_dict.get(key):
            read_stops_exons[key] = read_stops_dict.get(key)
            
    for_bgm = [x for y in [[exons.index(i)]*n for i, n in read_stops_exons.items()] for x in y]
    
    plt.figure(num=1)
    plt.hist(for_bgm,bins=num_bins)
    plt.title('Exon Peaks: CDS Position' + ' ' + gene_name)
    plt.xlabel("Distance from first exon start", fontsize=12)
    plt.ylabel("Count",fontsize=12)
    
    
    for_bgm = np.array(for_bgm).reshape([-1,1])
    result = BayesianGaussianMixture(3,max_iter=2000).fit(for_bgm).predict(for_bgm)
    
    exons_results = pd.DataFrame({'result':result,'pos':for_bgm.flatten()})
    exons_results['pos'] = exons_results['pos'].map(lambda x: exons[x])
    
    plt.figure(num=2)
    plt.hist(exons_results['pos'],bins=num_bins)
    plt.title('Exon Peaks: Genomic Position' + ' ' + gene_name)
    plt.xlabel("Distance from first exon start", fontsize=12)
    plt.ylabel("Count",fontsize=12)


    fig3, ax = plt.subplots(nrows=1, ncols=num_peaks, sharex=True, sharey=True)
    for i in range(0,num_peaks):
        group = exons_results[exons_results['result'] == i]
        #num_bins=int(np.array(group['pos']).std())
        ax[i].hist(group['pos'],bins=num_bins)
    
    fig3.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", bottom=False, left=False)
    
    plt.title('Exon Peaks Found' + ' ' + gene_name)
    plt.xlabel("Distance from first exon start", fontsize=12)
    plt.ylabel("Count",fontsize=12)
    

def graph_others(gene_name,df,num_peaks = 3, num_bins = 100):
    '''
    Produce 3 histogram graphs for read_stops outside exons before and after assignment by GMM algorithm. Note that the GMM model is non-deterministic 
    so sometimes the results will difer. Running this function repeatedly with a different num_peaks may give an idea of the best num_peaks parameter for 
    the count_and_summarize_peaks() function. 

    Parameters
    ----------
    gene_name : string
        gene id/symbol for gene of interest
    df : dataframe
        the dataframe returned after bam_top_genes()
    num_peaks : integer, optional
        maximum number of peaks for the GMM model. The default is 3.
    num_bins : integer, optional
        number of bins for histograms. The default is 100.

    Returns
    -------
    3 graphs
        1- read_stops outside exons genomic location histogram
        2- read_stops outside exons genomic location histogram, x-axis limited to gene 
        3- num_peaks histograms for each of the peaks found 

    '''
    #df=df2.copy()
    df=df[df['name']=="'"+gene_name+"'"]
    df=df.reset_index()
    
    if df.loc[0,'strand'] == "+":
        start = min(df.loc[0,'exons'])
        exons = list(map(lambda x: x - start, df.loc[0,'exons']))
        read_stops = list(map(lambda x: int(x) - start, df.loc[0,'read_stops']))
    if df.loc[0,'strand'] == "-":
        start = max(df.loc[0,'exons'])
        exons = list(map(lambda x: -1*(x - start), df.loc[0,'exons']))
        exons = exons[::-1]
        read_stops = list(map(lambda x: -1*(int(x) - start), df.loc[0,'read_stops']))
    
    read_stops_dict = defaultdict(int)
    
    for stop in read_stops:
        read_stops_dict[stop] += 1
        
    other_read_stops=defaultdict(int)
        
    for key in read_stops_dict.keys():
        if key not in exons:
            other_read_stops[key] = read_stops_dict.get(key)


    for_bgm = [x for y in [[i]*n for i, n in other_read_stops.items()] for x in y]
     
    
    for_bgm = np.array(for_bgm).reshape([-1,1])
    result = BayesianGaussianMixture(3,max_iter=2000).fit(for_bgm).predict(for_bgm)

    other_results =  pd.DataFrame({'result':result,'pos':for_bgm.flatten()})
    
    plt.figure(num=1)
    plt.hist(other_results['pos'],bins=num_bins)
    plt.title('Other Peaks: Genomic Position' + ' ' + gene_name)
    plt.xlabel("Distance from first exon start", fontsize=12)
    plt.ylabel("Count",fontsize=12)
    
    plt.figure(num=2)
    plt.hist(other_results['pos'],bins=num_bins)
    plt.title('Other Peaks: Inside Gene, Genomic Position' + ' ' + gene_name)
    plt.xlabel("Distance from first exon start", fontsize=12)
    plt.ylabel("Count",fontsize=12)
    plt.xlim(0,len(exons))
    
    
    
    fig3, ax = plt.subplots(nrows=1, ncols=num_peaks, sharex=True, sharey=True)
    for i in range(0,num_peaks):
        group = other_results[other_results['result'] == i]
        #num_bins=int(np.array(group['pos']).std())
        ax[i].hist(group['pos'],bins=num_bins)
    
    fig3.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", bottom=False, left=False)
    
    plt.title('Other Peaks Found' + ' ' + gene_name)
    plt.xlabel("Distance from first exon start", fontsize=12)
    plt.ylabel("Count",fontsize=12)
        
        
def new_h5(counts,summaries,dset_file,new_dset_file):
    '''
    Make a new .h5 file containing all the information from the original .h5 file
    but with the new alternate trascript counts added.
    
    Parameters
    ----------
    counts : dataframe
        counts dataframe from seperate_into_peaks() either filtered or not filtered
    summaries : dataframe
        summaries dataframe from seperate_into_peaks() either filtered or not filtered
    dset_file : .h5 filename
        a 10x filtered .h5 file
    new_dset_file : filename
        new .h5 file to be created in current directory
    
    '''
    #load the dset_file
    dset = h5py.File(dset_file, 'r')
    dset = dset['matrix']
    indices = dset['indices']
    data = dset['data']
    
    #create the dset2 file and matrix group with barcodes and features group
    subprocess.Popen(('touch', new_dset_file))
    dset2 = h5py.File(new_dset_file, 'w')
    dset2.create_group('matrix')
    dset2['matrix'].create_dataset('barcodes', data = np.array(dset['barcodes']))
    dset2['matrix'].create_group('features')
    
    #add the new features to the old features ids
    summaries = summaries.reset_index()
    summaries['result'] = summaries['result'].str.replace("'","")
    old_genes = np.array(dset['features']['id'],dtype="|S30")
    new_genes = np.array(summaries['result'],dtype='|S30')
    dset2['matrix']['features'].create_dataset('id', data = np.append(old_genes,new_genes), dtype='|S30')
    
    #add the new features to the old features names
    old_genes_names = np.array(dset['features']['name'])
    dset2['matrix']['features'].create_dataset('name', data = np.append(old_genes_names, new_genes), dtype="|S30")
    
    #add extended genome array
    genome = dset['features']['genome'][0]
    genome = np.repeat(genome,dset2['matrix']['features']['id'].shape[0])
    dset2['matrix']['features'].create_dataset('genome', data = genome, dtype = "|S6")
    
    #add extended feature type array 
    feature_type = dset['features']['feature_type'][0]
    feature_type = np.repeat(feature_type,dset2['matrix']['features']['id'].shape[0])
    dset2['matrix']['features'].create_dataset('feature_type', data=feature_type, dtype = "|S15")
    
    #add tag
    tag = dset['features']['_all_tag_keys'][0]
    tag = np.array(tag,dtype='|S6')
    dset2['matrix']['features'].create_dataset('_all_tag_keys', data = tag, dtype = "|S6")
    
    #get a list of all features
    gene_list = list(pd.Series(dset2['matrix']['features']['id']).astype('str').str.replace('b','').str.replace("'",''))
    
    #get a list of all barcodes
    barcodes = list(pd.Series(dset['barcodes']).astype('str').str.replace('b','').str.replace("'",''))
    
    #update the indptr, data and indices arrays from the original dset by inserting new data
    indptr = np.array(dset['indptr'])
    data = np.array(dset['data'])
    indices = np.array(dset['indices'])
    counts = counts.reset_index()
    counts['result'] = counts['result'].str.replace("'",'')
    
    for i,b in enumerate(sorted(list(counts['barcode'].unique()))):
        if i % 1000 == 0:
            print(i)
        if b in barcodes:
            cell_counts = counts[counts['barcode']==b]
            cell_counts['gene_indice']=cell_counts['result'].map(lambda x: gene_list.index(x))
            cell_counts = cell_counts.sort_values(by='gene_indice')
            i = barcodes.index(b)
            pos = indptr[i]
            data = np.insert(data,pos,cell_counts['count'])              
            indices = np.insert(indices,pos,cell_counts['gene_indice'])
            indptr = np.append(indptr[:i+1],indptr[i+1:]+len(cell_counts))
        else: 
            continue
    
    #put the updated arrays in the new dset
    dset2['matrix'].create_dataset('indptr', data = indptr)
    dset2['matrix'].create_dataset('data', data = data)
    dset2['matrix'].create_dataset('indices', data = indices)
    
    #put the new matrix shape in the new dset
    new_shape =np.array([len(dset2['matrix']['features']['id']), len(dset2['matrix']['barcodes'])], dtype='int32')
    dset2['matrix'].create_dataset('shape', data = new_shape, dtype='int32')
           
               
        


if __name__ == "__main__":
    
    #Load required packages
    import os, subprocess, h5py, re
    import pandas as pd
    import numpy as np
    from sklearn.mixture import BayesianGaussianMixture
    from collections import defaultdict
    import matplotlib.pyplot as plt
    #also samtools
    
    #Add a lot of directories to path so that samtools is in path
    os.environ['PATH'] = '/opt/anaconda3/bin:/opt/anaconda3/condabin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin'

   




