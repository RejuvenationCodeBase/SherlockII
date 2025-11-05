# es_toolbox.py

# Version 18: Last Edit May 20, 2014.
# Contact:  chriskfuller@gmail.com  or chris@genome.ucsf.edu
import pdb
# Toolbox is the core collection of functions used in the Empirical Sherlock (ES) algorithm for eQTL + GWAS analysis
#import MySQLdb
import sys, os, math, numpy, string, pickle, time, tabix
from subprocess import PIPE, Popen, STDOUT
from scipy import stats
from scipy import signal
from scipy.misc import *
import scipy.stats 
from es_toolbox import *
from es_object_defs import *
from es_graphics import *
#import rpy2.robjects as robjects
#from rpy2.robjects.packages import importr
from scipy import *
from pylab import *
#from new_lib import *
import random
import bisect
random.seed()
from numpy import asarray, array
from numpy.fft import rfft, irfft
import qvalue_Storey_Tibshirani_2003

#from statsmodels.stats.multitest import fdrcorrection

####################################################################################################################
###   Functions to Process the eQTL Data, align eQTL with GWAS, and tag the matched pairs.
###   Essentially, this is the pre-analysis required before we can begin applying our statistical test to the data.
###   This approach contains a number of heuristics.  A better approach would be to use full genotypes for the
###   GWAS and eQTL data to determine linkage and identify the "true" independent blocks in the data.
###   Given that genotypes are only rarely available, we rely on a database of linkage generated from 1000 genomes.
####################################################################################################################

def get_gene_list_txt(define_gene_list, eQTL, data_folder_location):
    """Make any necessary definitions, txt version
    Define a function for obtaining the full gene list using either
    all genes in a given eQTL or imported from a text file
    """
    gene_names = []
    if define_gene_list == 'ALL':
        geneFile = data_folder_location + 'eQTL/' + eQTL + '.gene'
        # Use every gene name present in the txt file of
        # a given eQTL data set.
        with open(geneFile) as geneInput:
            # each line is one gene
            gene = geneInput.readline().strip()
            while gene:
                gene_names.append(gene)
                gene = geneInput.readline().strip()
    else:
        gene_names = define_gene_list.split(',')

    return gene_names

def update_metabolic_pickle(eQTL_datasets, define_gene_list, \
        database_cursor, data_folder_location, include_HLA_region):
    """
    update metabolic pickle, from local file (no database required)
    This funciton is called, generally, only when introduce new metabolic
    data sets
    """
    # Will take considerable time (hour or more) per eQTL dataset
    for eQTL in eQTL_datasets:
        print ('Starting analysis of eQTL: ' + eQTL)
        #### Load all gene names found in the eQTL file
        gene_names = get_gene_list_txt(define_gene_list, eQTL,
                data_folder_location)

        ### Setup the Gene List ###
        gene_list = []
        print ("Found " + str(len(gene_names)) + " genes in " + eQTL + " eQTL:")
        hugoDB = getHUGODB(database_cursor)
        geneSymbolToID = getGeneSymbolToID(database_cursor)
        geneIDPos = getGeneIDPos(database_cursor)
        allSNPs, geneSNPs = getAllSNPs_local(eQTL, data_folder_location)
        print (geneSNPs)
        #snpAnnotation was used as a dictionary to store
        #annotation for all snps, one snp one record
        snpAnnotation = getSNPAnnotation(database_cursor, allSNPs)
        #snpInLD = getSNPsInLD(database_cursor, snpAnnotation)
        snp_pos, snp_ld, snp_linked = get_eSNP_in_LD02(eQTL,
                data_folder_location)
        #get the preloaded in_ld information including position range LD>=0.2
        ld_range = get_eSNP_LD_range(data_folder_location)
        for i in range(len(gene_names)):
            tmp_Gene = Gene()
            tmp_Gene = annotate_Gene2(tmp_Gene, gene_names[i], hugoDB,
                    geneSymbolToID, geneIDPos)
            tmp_Gene.eQTL = eQTL
            gene_list.append(tmp_Gene)
            #print gene_list[i].name

        # Add SNPs from eQTL to each gene, properly sorted and listed by block
        chrom_sort_dict = {'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22, 'chrX':23, 'chrY':24}
        cntr = 0
        #database_cursor.execute("use eQTL")
        annotated_SNP = {}
        snpAnnotated = 0
        savedQuery = 0
        for gene in gene_list:
            cntr = cntr + 1
            print ("Get eQTL SNPs for gene " + str(int(cntr)) + " of " + str(len(gene_list)) + " : " + str(gene.name))
            SNPs = []
            eQTL_pvals_dict = {}
            #database_cursor.execute("use eQTL")
            #SQL_string = "SELECT * from " + eQTL + " where gene='" + gene.name + "'"
            #database_cursor.execute(SQL_string)
            #qr = database_cursor.fetchall()
            for snp1 in geneSNPs[gene.name]:
                tmp = None
                if not (snp1[0] in annotated_SNP):
                    tmp = SNP()
                    #tmp = annotate_SNP(tmp, qr[i][1], database_cursor)
                    #use the pre loaded snpAnnotation to annotate the snps
                    #add chromsome, MAF, positions etc
                    tmp = annotate_SNP2(tmp, snp1[0], snpAnnotation)
                    annotated_SNP[snp1[0]]=tmp
                    snpAnnotated += 1
                else:
                    savedQuery += 1
                    tmp = annotated_SNP[snp1[0]]

                if tmp.chrom != 'missing' and tmp.chrom in chrom_sort_dict:
                    SNPs.append(tmp)
                    eQTL_pvals_dict[snp1[0]] = snp1[1]

            if len(SNPs) != 0:  # Not all genes will be found in each eQTL data set

                # Sort SNPs by chromosome, then position
                SNPs = sort_SNPs(SNPs)

                # Remove repeat SNPs (strangely, some data sets have these!!!!). Keep only the first one.
                updated_list = []
                for i in range(len(SNPs)):
                    if i == 0:
                        updated_list.append(SNPs[i])
                    else:
                        if SNPs[i].rsID == SNPs[i-1].rsID:
                            pass
                        else:
                            updated_list.append(SNPs[i])
                SNPs = updated_list

                # Update the number of individual SNPs found for the gene
                gene.num_SNPs = len(SNPs)

                # Load SNPs into the SNP pairs
                SNP_pairs = []
                count = 0
                for item in SNPs:
                    if include_HLA_region != 'True' and item.chrom == 'chr6' and int(item.location) > 29000000 and int(item.location) < 33000000:
                        # Discard SNPs in the HLA if requested
                        pass
                    else:
                        count = count + 1
                        tmp = SNP_pair()
                        tmp.eQTL_SNP = item.rsID
                        tmp.eQTL_SNP_pval = eQTL_pvals_dict[item.rsID]
                        tmp.block_number = count
                        tmp.chrom = item.chrom
                        tmp.eQTL_SNP_location = item.location
                        tmp.eQTL_SNP_MAF = item.MAF
                        SNP_pairs.append(tmp)

                # Add these blocks to our gene object before returning
                gene.SNP_pairs = SNP_pairs

        print ("remove gene with no SNP")

        # Remove genes that don't have any SNPs (for whatever reason)
        updated_gene_list = []
        for i in range(len(gene_list)):
            if gene_list[i].SNP_pairs == 'Not Loaded':
                print ('Removing Gene with no SNP Pairs: ' + gene_list[i].name)
            else:
                updated_gene_list.append(gene_list[i])

        del gene_list
        gene_list = updated_gene_list
        del updated_gene_list

        print ("Identify SNPs in the LD")
        # Identify SNPs in the LD database
        gene_list = is_in_LD_database2(gene_list, snp_pos, ld_range)

        # Set tag SNP label to "No".  Set this to yes only for confirmed tag SNPs
        for gene in gene_list:
            for item in gene.SNP_pairs:
                item.is_tag_SNP = 'No'

        # Add gene proximity information for each SNP (either cis or trans)
        gene_list = determine_SNP_gene_proximity(gene_list)

        #below is duplicated
        # Remove genes that don't have any SNPs (for whatever reason)
        #updated_gene_list = []
        #for i in range(len(gene_list)):
        #    if gene_list[i].SNP_pairs == 'Not Loaded':
        #        print 'Removing Gene with no SNP Pairs: ' + gene_list[i].name
        #    else:
        #        updated_gene_list.append(gene_list[i])

        #del gene_list
        #gene_list = updated_gene_list
        #del updated_gene_list

        # Save to a pickle
        output_file_name = eQTL + '.pickle'
        FileWriter = open(data_folder_location + 'eQTL/' + output_file_name, 'w')
        pickle.dump(gene_list, FileWriter)
        FileWriter.close()



def update_eQTL_pickle_local(eQTL_datasets, define_gene_list, database_cursor,
        data_folder_location, include_HLA_region):
    """
    update eQTL pickle, from local file (no database required)
    This funciton is called, generally, only when introduce new eQTL data sets
    """
    # Will take considerable time (hour or more) per eQTL dataset
    for eQTL in eQTL_datasets:
        print ('Starting analysis of eQTL: ' + eQTL)

        #### Load all gene names found in the eQTL file
        gene_names = get_gene_list_txt(define_gene_list, eQTL,
                data_folder_location)

        ### Setup the Gene List ###
        gene_list = []
        print ("Found " + str(len(gene_names)) + " genes in " + eQTL + " eQTL:")
        hugoDB = getHUGODB(database_cursor)
        geneSymbolToID = getGeneSymbolToID(database_cursor)
        geneIDPos = getGeneIDPos(database_cursor)
        allSNPs, geneSNPs = getAllSNPs_local(eQTL, data_folder_location)
        #snpAnnotation was used as a dictionary to store
        #annotation for all snps, one snp one record
        snpAnnotation = getSNPAnnotation(database_cursor, allSNPs)
        #snpInLD = getSNPsInLD(database_cursor, snpAnnotation)
        snp_pos, snp_ld, snp_linked = get_eSNP_in_LD02(eQTL,
                data_folder_location)

        #get the preloaded in_ld information including position range LD>=0.2
        ld_range = get_eSNP_LD_range(data_folder_location)

        for i in range(len(gene_names)):
            tmp_Gene = Gene()
            tmp_Gene = annotate_Gene2(tmp_Gene, gene_names[i], hugoDB,
                    geneSymbolToID, geneIDPos)
            tmp_Gene.eQTL = eQTL
            gene_list.append(tmp_Gene)
            #print gene_list[i].name

        # Add SNPs from eQTL to each gene, properly sorted and listed by block
        chrom_sort_dict = {'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22, 'chrX':23, 'chrY':24}
        cntr = 0
        #database_cursor.execute("use eQTL")
        annotated_SNP = {}
        snpAnnotated = 0
        savedQuery = 0
        for gene in gene_list:
            cntr = cntr + 1
            print ("Get eQTL SNPs for gene " + str(int(cntr)) + " of " + str(len(gene_list)) + " : " + str(gene.name))
            SNPs = []
            eQTL_pvals_dict = {}
            #database_cursor.execute("use eQTL")
            #SQL_string = "SELECT * from " + eQTL + " where gene='" + gene.name + "'"
            #database_cursor.execute(SQL_string)
            #qr = database_cursor.fetchall()
            for snp1 in geneSNPs[gene.name]:
                tmp = None
                if not (snp1[0] in annotated_SNP):
                    tmp = SNP()
                    #tmp = annotate_SNP(tmp, qr[i][1], database_cursor)
                    #use the pre loaded snpAnnotation to annotate the snps
                    #add chromsome, MAF, positions etc
                    tmp = annotate_SNP2(tmp, snp1[0], snpAnnotation)
                    annotated_SNP[snp1[0]]=tmp
                    snpAnnotated += 1
                else:
                    savedQuery += 1
                    tmp = annotated_SNP[snp1[0]]

                if tmp.chrom != 'missing' and tmp.chrom in chrom_sort_dict:
                    SNPs.append(tmp)
                    eQTL_pvals_dict[snp1[0]] = snp1[1]

            if len(SNPs) != 0:  # Not all genes will be found in each eQTL data set

                # Sort SNPs by chromosome, then position
                SNPs = sort_SNPs(SNPs)

                # Remove repeat SNPs (strangely, some data sets have these!!!!). Keep only the first one.
                updated_list = []
                for i in range(len(SNPs)):
                    if i == 0:
                        updated_list.append(SNPs[i])
                    else:
                        if SNPs[i].rsID == SNPs[i-1].rsID:
                            pass
                        else:
                            updated_list.append(SNPs[i])
                SNPs = updated_list

                # Update the number of individual SNPs found for the gene
                gene.num_SNPs = len(SNPs)

                # Load SNPs into the SNP pairs
                SNP_pairs = []
                count = 0
                for item in SNPs:
                    if include_HLA_region != 'True' and item.chrom == 'chr6' and int(item.location) > 29000000 and int(item.location) < 33000000:
                        # Discard SNPs in the HLA if requested
                        pass
                    else:
                        count = count + 1
                        tmp = SNP_pair()
                        tmp.eQTL_SNP = item.rsID
                        tmp.eQTL_SNP_pval = eQTL_pvals_dict[item.rsID]
                        tmp.block_number = count
                        tmp.chrom = item.chrom
                        tmp.eQTL_SNP_location = item.location
                        tmp.eQTL_SNP_MAF = item.MAF
                        SNP_pairs.append(tmp)

                # Add these blocks to our gene object before returning
                gene.SNP_pairs = SNP_pairs

        print ("remove gene with no SNP")

        # Remove genes that don't have any SNPs (for whatever reason)
        updated_gene_list = []
        for i in range(len(gene_list)):
            if gene_list[i].SNP_pairs == 'Not Loaded':
                print ('Removing Gene with no SNP Pairs: ' + gene_list[i].name)
            else:
                updated_gene_list.append(gene_list[i])

        del gene_list
        gene_list = updated_gene_list
        del updated_gene_list

        print ("Identify SNPs in the LD")
        # Identify SNPs in the LD database
        gene_list = is_in_LD_database2(gene_list, snp_pos, ld_range)

        # Set tag SNP label to "No".  Set this to yes only for confirmed tag SNPs
        for gene in gene_list:
            for item in gene.SNP_pairs:
                item.is_tag_SNP = 'No'

        # Add gene proximity information for each SNP (either cis or trans)
        gene_list = determine_SNP_gene_proximity(gene_list)

        #below is duplicated
        # Remove genes that don't have any SNPs (for whatever reason)
        #updated_gene_list = []
        #for i in range(len(gene_list)):
        #    if gene_list[i].SNP_pairs == 'Not Loaded':
        #        print 'Removing Gene with no SNP Pairs: ' + gene_list[i].name
        #    else:
        #        updated_gene_list.append(gene_list[i])

        #del gene_list
        #gene_list = updated_gene_list
        #del updated_gene_list

        # Save to a pickle
        output_file_name = eQTL + '.pickle'
        FileWriter = open(data_folder_location + 'eQTL/' + output_file_name, 'w')
        pickle.dump(gene_list, FileWriter)
        FileWriter.close()



# This funciton is called, generally, only when we introduce new eQTL data sets into the database
def update_eQTL_pickle(eQTL_datasets, define_gene_list, database_cursor, data_folder_location, include_HLA_region):

    # Run loop for every eQTL dataset.  This will take considerable time (hour or more) per eQTL dataset
    for eQTL in eQTL_datasets:
        print ('Starting analysis of eQTL: ' + eQTL)

        #### Load all gene names found in the eQTL file
        if define_gene_list == 'ALL':
            gene_names = []
            database_cursor.execute("use eQTL")
            SQL_string = "SELECT DISTINCT gene from " + eQTL
            database_cursor.execute(SQL_string)
            qr = database_cursor.fetchall()
            if qr == ():
                print ("Nothing Found searching for eQTL gene names")
            else:
                for i in range(len(qr)):
                    gene_names.append(qr[i][0])
        else:
            gene_names = define_gene_list


        ### Setup the Gene List ###
        gene_list = []
        print ("Found " + str(len(gene_names)) + " genes in " + eQTL + " eQTL:")
        hugoDB = getHUGODB(database_cursor)
        geneSymbolToID = getGeneSymbolToID(database_cursor)
        geneIDPos = getGeneIDPos(database_cursor)
        allSNPs, geneSNPs = getAllSNPs(database_cursor, eQTL)
        #snpAnnotation was used as a dictionary to store
        #annotation for all snps, one snp one record
        snpAnnotation = getSNPAnnotation(database_cursor, allSNPs)
        snpInLD = getSNPsInLD(database_cursor, snpAnnotation)

        for i in range(len(gene_names)):
            tmp_Gene = Gene()
            tmp_Gene = annotate_Gene2(tmp_Gene, gene_names[i], hugoDB, geneSymbolToID, geneIDPos)
            tmp_Gene.eQTL = eQTL
            gene_list.append(tmp_Gene)
            #print gene_list[i].name

        # Add SNPs from eQTL to each gene, properly sorted and listed by block
        chrom_sort_dict = {'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22, 'chrX':23, 'chrY':24}
        cntr = 0
        database_cursor.execute("use eQTL")
        annotated_SNP = {}
        snpAnnotated = 0
        savedQuery = 0
        for gene in gene_list:
            cntr = cntr + 1
            print ("Get eQTL SNPs for gene " + str(int(cntr)) + " of " + str(len(gene_list)) + " : " + str(gene.name))
            SNPs = []
            eQTL_pvals_dict = {}
            #database_cursor.execute("use eQTL")
            #SQL_string = "SELECT * from " + eQTL + " where gene='" + gene.name + "'"
            #database_cursor.execute(SQL_string)
            #qr = database_cursor.fetchall()
            for snp1 in geneSNPs[gene.name]:
                tmp = None
                if not (snp1[0] in annotated_SNP):
                    tmp = SNP()
                    #tmp = annotate_SNP(tmp, qr[i][1], database_cursor)
                    #use the pre loaded snpAnnotation to annotate the snps
                    #add chromsome, MAF, positions etc
                    tmp = annotate_SNP2(tmp, snp1[0], snpAnnotation)
                    annotated_SNP[snp1[0]]=tmp
                    snpAnnotated += 1
                else:
                    savedQuery += 1
                    tmp = annotated_SNP[snp1[0]]

                if tmp.chrom != 'missing' and tmp.chrom in chrom_sort_dict:
                    SNPs.append(tmp)
                    eQTL_pvals_dict[snp1[0]] = snp1[1]

            if len(SNPs) != 0:  # Not all genes will be found in each eQTL data set

                # Sort SNPs by chromosome, then position
                SNPs = sort_SNPs(SNPs)

                # Remove repeat SNPs (strangely, some data sets have these!!!!). Keep only the first one.
                updated_list = []
                for i in range(len(SNPs)):
                    if i == 0:
                        updated_list.append(SNPs[i])
                    else:
                        if SNPs[i].rsID == SNPs[i-1].rsID:
                            pass
                        else:
                            updated_list.append(SNPs[i])
                SNPs = updated_list

                # Update the number of individual SNPs found for the gene
                gene.num_SNPs = len(SNPs)

                # Load SNPs into the SNP pairs
                SNP_pairs = []
                count = 0
                for item in SNPs:
                    if include_HLA_region != 'True' and item.chrom == 'chr6' and int(item.location) > 29000000 and int(item.location) < 33000000:
                        # Discard SNPs in the HLA if requested
                        pass
                    else:
                        count = count + 1
                        tmp = SNP_pair()
                        tmp.eQTL_SNP = item.rsID
                        tmp.eQTL_SNP_pval = eQTL_pvals_dict[item.rsID]
                        tmp.block_number = count
                        tmp.chrom = item.chrom
                        tmp.eQTL_SNP_location = item.location
                        tmp.eQTL_SNP_MAF = item.MAF
                        SNP_pairs.append(tmp)

                # Add these blocks to our gene object before returning
                gene.SNP_pairs = SNP_pairs

        print ("remove gene with no SNP")

        # Remove genes that don't have any SNPs (for whatever reason)
        updated_gene_list = []
        for i in range(len(gene_list)):
            if gene_list[i].SNP_pairs == 'Not Loaded':
                print ('Removing Gene with no SNP Pairs: ' + gene_list[i].name)
            else:
                updated_gene_list.append(gene_list[i])

        del gene_list
        gene_list = updated_gene_list
        del updated_gene_list

        print ("Identify SNPs in the LD")
        # Identify SNPs in the LD database
        gene_list = is_in_LD_database2(gene_list, snpInLD)

        # Set tag SNP label to "No".  Set this to yes only for confirmed tag SNPs
        for gene in gene_list:
            for item in gene.SNP_pairs:
                item.is_tag_SNP = 'No'

        # Add gene proximity information for each SNP (either cis or trans)
        gene_list = determine_SNP_gene_proximity(gene_list)

        #below is duplicated
        # Remove genes that don't have any SNPs (for whatever reason)
        #updated_gene_list = []
        #for i in range(len(gene_list)):
        #    if gene_list[i].SNP_pairs == 'Not Loaded':
        #        print 'Removing Gene with no SNP Pairs: ' + gene_list[i].name
        #    else:
        #        updated_gene_list.append(gene_list[i])

        #del gene_list
        #gene_list = updated_gene_list
        #del updated_gene_list

        # Save to a pickle
        output_file_name = eQTL + '.pickle'
        FileWriter = open(data_folder_location + 'eQTL/' + output_file_name, 'w')
        pickle.dump(gene_list, FileWriter)
        FileWriter.close()

def get_eSNP_LD_range(data_folder_location):
    """Return a dictionary of snps, with location, whether it is in LD database
    and its LD range"""
    ld_range = {}
    snp_LDDB_range = data_folder_location + 'eQTL/allSNPs.LDDB'
    with open(snp_LDDB_range) as LD:
        for line in LD:
            items = line.strip().split('\t')
            ld_range[items[0]]=items
            chr_pos = items[1]+":"+items[2]
            ld_range[chr_pos]=items
    return ld_range

def is_snps_in_LD_range(ld_range, snp1, snp2, snp1_pos, snp2_pos,\
        default_range=1000000):
    """Test whether two snps are in the LD range with each other
    This is called when snp1 and snp2 are on the same chromosome
    So this function does not check the chromosome
    """
    if snp1 in ld_range and ld_range[snp1][5] != 'None':
        left1 = int(ld_range[snp1][5])
        right1 = int(ld_range[snp1][6])
        if snp2_pos >=left1 and snp2_pos <= right1:
            return True
        else:
            return False
    elif snp2 in ld_range and ld_range[snp2][5] != 'None':
        left2 = int(ld_range[snp2][5])
        right2 = int(ld_range[snp2][6])
        if snp1_pos >= left2 and snp1_pos <= right2:
            return True
        else:
            return False
    else:
        if abs(snp1_pos-snp2_pos) <= default_range:
            return True
        else:
            return False

def is_snp1_in_LD_snp2(ld_range, snp1, snp2, snp1_pos, snp2_pos):
    """Test whether two snp1 is in the LD range of snp2
    This is called when snp1 and snp2 are on the same chromosome
    So this function does not check the chromosome
    This function assume snp2 is in the ld_range database therefore
    will for sure have the range.
    This is for comparing a snp not in LD database (snp1) with
    another snp in LD database (snp2)
    """
    if snp2 in ld_range and ld_range[snp2][5] != 'None':
        left2 = int(ld_range[snp2][5])
        right2 = int(ld_range[snp2][6])
        if snp1_pos >= left2 and snp1_pos <= right2:
            return True
        else:
            return False
    else:
        #should not happen
        return False


def get_eSNP_in_LD02(eQTL, data_folder_location,LD_window_cutoff=0.2):
    """Read the preprocessing eQTL snp LD02 file (LD >0.2)
    data_folder_location + 'eQTL/' + eQTL + '.LD02'
    """
    eQTL_LD_file = data_folder_location + 'eQTL/' + eQTL + '.LD02'
    print ("processing",eQTL_LD_file)
    eSNP_pos = {} #{snp:[chrom, pos]}
    eSNP_lds = {} #{snpA:{snpB:ld}}
    eSNP_linked = {} #{snpA:[snpB1,snpB2,...]}, snpB sorted by position
    i_line = 0
    #pdb.set_trace()
    with open(eQTL_LD_file) as ld:
        for line in ld:
        #line = ld.readline()
        #while line:
            items = line.strip().split("\t")
            #data format [rsid, chr, pos1, rsid1, pos2, rsid2, ld]
            #in most cases rsid == rsid1 (rsid is the snp appears in eQTL)
            eSNP_pos.setdefault(items[0],[items[1],items[2]])
            eSNP_pos.setdefault(items[5],[items[1],items[4]])
            eSNP_lds.setdefault(items[0],{})
            if float(items[6])<LD_window_cutoff: 
                continue # raise the LD window cutoff above 0.2 
            eSNP_lds[items[0]][items[5]]=float(items[6])
            eSNP_linked.setdefault(items[0],[])
            eSNP_linked[items[0]].append(items[5])
            if i_line % 10000 == 0:
                print ('processing line# ',i_line, 'in ',eQTL_LD_file)
            i_line += 1
            #line = ld.readline()
    return eSNP_pos, eSNP_lds, eSNP_linked

def get_eSNP_in_LD85(eQTL, data_folder_location):
    """Read the preprocessing eQTL snp LD85 file (LD >0.85)
    data_folder_location + 'eQTL/' + eQTL + '.LD85'
    To prepare for alignment with preload LD data

    Note: some SNPs do share the same location
    use the eSNP_linked to preserve the order of snps
    """
    eQTL_LD_file = data_folder_location + 'eQTL/' + eQTL + '.LD85'
    eSNP_pos = {} #{snp:[chrom, pos]}
    eSNP_linked = {} #{snpA:[[snpB1,ld2],[snpB2,ld2],...]} sorted by position
    with open(eQTL_LD_file) as ld:
        line = ld.readline()
        while line:
            items = line.strip().split("\t")
            #data format [rsid, chr, pos1, rsid1, pos2, rsid2, ld]
            #in most cases rsid == rsid1 (rsid is the snp appears in eQTL)
            eSNP_pos.setdefault(items[0],[items[1],items[2]])
            eSNP_pos.setdefault(items[5],[items[1],items[4]])
            ld_value = float(items[6])
            eSNP_linked.setdefault(items[0],[])
            eSNP_linked[items[0]].append([items[5], ld_value])
            line = ld.readline()
    return eSNP_pos, eSNP_linked


def read_eQTL_snps(eQTL, data_folder_location):
    """Read all eQTL_snps from local file
    """
    eQTL_file = data_folder_location + 'eQTL/' + eQTL + '.info'
    eQTL_snps = {}
    with open(eQTL_file) as info:
        line = info.readline()
        while line:
            items = line.strip().split('\t')
            eQTL_snps.setdefault(items[1],1)
            line = info.readline()
    return eQTL_snps

def getAllSNPs(database_cursor, eQTL):
    database_cursor.execute("use eQTL")
    SQL_string = "SELECT gene, snp, pval from " + eQTL
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    allSNPs = {}
    geneSNPs = {}
    index = 0
    for r in qr:
        index += 1
        geneSNPs.setdefault(r[0],[])
        geneSNPs[r[0]].append([r[1],r[2]])
        if not (r[1] in allSNPs):
            allSNPs[r[1]]=1
        if index % 100000 == 0:
            print (index)
    return allSNPs, geneSNPs

def getAllSNPs_local(eQTL, data_folder_location):
    """Read all SNPs from local .info file
    return allSNPs, geneSNPs={gene:[[snp1,pval1]...]}
    """
    eQTL_file = data_folder_location + 'eQTL/' + eQTL + '.info'
    allSNPs = {}
    geneSNPs = {}
    with open(eQTL_file) as info:
        line = info.readline()
        while line:
            items = line.split('\t')
            geneSNPs.setdefault(items[0],[])
            geneSNPs[items[0]].append([items[1],items[2]])
            allSNPs.setdefault(items[1],items[3:])
            line = info.readline()
    return allSNPs, geneSNPs

def get_eSNP_in_LD85_filter(eQTL, data_folder_location, gwas_snp, eQTL_snps,align_LD_threshold = 0.85):
    """Read the preprocessing eQTL snp LD85 file (LD >0.85)
    data_folder_location + 'eQTL/' + eQTL + '.LD85'
    To prepare for alignment with preload LD data
    only keep those that are in gwas_snp and eQTL_snps
    Note: some SNPs do share the same location
    Use eSNP_linked to preserve the order of snps
    """
    #eQTL_LD_file = data_folder_location + 'eQTL/' + eQTL + '.LD85'
    eQTL_LD_file = data_folder_location + 'eQTL/' + eQTL + '.LD50'
    eSNP_pos = {} #{snp:[chrom, pos]}
    eSNP_linked = {} #{snpA:[[snpB1,ld2],[snpB2,ld2],...]} sorted by position
    with open(eQTL_LD_file) as ld:
        #line = ld.readline()
        #while line:
        for line in ld:
            items = line.strip().split("\t")
            if items[5] in gwas_snp or items[5] in eQTL_snps: #filter snpB
                #data format [rsid, chr, pos1, rsid1, pos2, rsid2, ld]
                #in most cases rsid == rsid1 (rsid is the snp appears in eQTL)
                eSNP_pos.setdefault(items[0],[items[1],items[2]])
                eSNP_pos.setdefault(items[5],[items[1],items[4]])
                ld_value = float(items[6])
                if ld_value < align_LD_threshold: continue
                eSNP_linked.setdefault(items[0],[])
                eSNP_linked[items[0]].append([items[5], ld_value])
            #line = ld.readline()
    return eSNP_pos, eSNP_linked


# The following functions upto match_eQTL_SNPs_with_GWAS
# are called from match_eQTL_SNPs_with_GWAS
# The original setup was nested function which make the function huge
# moved them out for good
def get_all_GWAS_SNP(GWAS,database_cursor):
    """Get all GWAS snp for a GWAS(gwas_id) from database
    """
    database_cursor.execute("use GWAS")
    SQL_string = "SELECT snp, pval, chrom, location, maf from " + GWAS +";"
    #print SQL_string
    database_cursor.execute(SQL_string)
    qrs = database_cursor.fetchall()
    gwas_snp={}
    for qr in qrs:
        snp = qr[0]
        pval = qr[1]
        chrom = qr[2]
        location = qr[3]
        maf = qr[4]
        output = [snp, pval, chrom, location, maf]
        gwas_snp[snp]=output
    return qrs, gwas_snp

def get_all_GWAS_SNP_local(GWAS, gwas_file):
    """Get all gwas snp from local text file, the result should be the
    same as get_all_GWAS_SNP() which use database query
    """
    gwas_snp={}
    with open(gwas_file) as gwas_input:
        for line in gwas_input:
            tmp=line.split('\t')
            snp = tmp[0]
            pval = float(tmp[1])
            chrom = tmp[2]
            location = int(tmp[3])
		#print tmp[4]
            maf = 0
            if len(tmp)>=5 and (not (tmp[4]=='none')):
                maf = float(tmp[4])
            output = [snp, pval, chrom, location, maf]
            gwas_snp[snp]=output
    return gwas_snp

def get_all_GWAS_SNP_pos_local(GWAS, gwas_file):
    """Get all gwas snp from local text file
    use chr:pos as key instead of the rsID
    """
    gwas_snp_pos={}
    with open(gwas_file) as gwas_input:
        for line in gwas_input:
            tmp=line.split('\t')
            snp = tmp[0]
            pval = float(tmp[1])
            chrom = tmp[2]
            location = long(tmp[3])
		#print tmp[4]
            maf = 0
            if len(tmp)>=5 and (not (tmp[4]=='none')):
                maf = float(tmp[4])
            output = [snp, pval, chrom, location, maf]
            chr_pos = chrom+":"+str(location)
            #assume the location are unique for different snps
            #this is just temperory solution (need to be changed)
            gwas_snp_pos[chr_pos]=output
    return gwas_snp_pos

def get_matching_GWAS_SNP(GWAS, rsID, database_cursor):
    """Identify the matching gwas snp,
    input: GWAS, name of GWAS table(gwas_id)
    input: rsID, find the data using rsID as key the database table
    input: database_cursor
    """
    database_cursor.execute("use GWAS")
    SQL_string = "SELECT snp, pval, chrom, location, maf from " + GWAS + " WHERE snp='" + rsID + "';"
    #print SQL_string
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    if qr == ():
        output = 'none'
    else:
        snp = qr[0][0]
        pval = qr[0][1]
        chrom = qr[0][2]
        location = qr[0][3]
        maf = qr[0][4]
        output = [snp, pval, chrom, location, maf]
    return output


def get_matching_GWAS_SNP_dict(gwas_snp, rsID):
    """Identify the matching gwas snp,
    input: gwas_snp is a dictionary of {rsid: [snp,pval,chrom,location,maf]}
    input: rsID, find the data using rsID as key in gwas_snp
    """
    output = 'none'
    if rsID in gwas_snp:
        snp = gwas_snp[rsID][0]
        pval = gwas_snp[rsID][1]
        chrom = gwas_snp[rsID][2]
        location = gwas_snp[rsID][3]
        maf = gwas_snp[rsID][4]
        #use the element to create a new array
        #instead of the gwas_snp[rsID] itself, which causing interference
        #array in order of snp,pval,chrom,location,maf
        output = [snp, pval, chrom, location, maf]

    return output


def get_high_LD_GWAS_SNPs_tabix(GWAS, rsID, chrom, pos, LD_threshold,
        data_folder_location, gwas_snp, eqtlsnp_LD):
    """For a rsID (eQTL snp), find the corresponding GWAS SNPs which are in
    high LD with it.
    LD_threshold is the threshold of LD (usually 0.85 here)
    The query using tabix should be the same as using the following sql
    from hg19.LD_chr[1,2,3...]

    SQL_string1 = "SELECT SNPA, rsq, SNPA_loc from LD_" + SNP_pair.chrom +
                  "  WHERE SNPB='" + SNP_pair.eQTL_SNP + "'" + " AND rsq > " +
                  str(LD_threshold)
    SQL_string2 = "SELECT SNPB, rsq, SNPB_loc from LD_" + SNP_pair.chrom +
                  "  WHERE SNPA='" + SNP_pair.eQTL_SNP + "'" + " AND rsq > " +
                  str(LD_threshold)
    """
    # First, create a list of ALL SNPs in LD with the eQTL SNP at better than
    # the LD_rsq threshold
    SNPs_in_LD = []
    # check if rsID has been query before (if so will be stored in eqtlsnp_LD
    if rsID in eqtlsnp_LD:
        #print "already known the SNPs_in_LD for " + rsID
        SNPs_in_LD = eqtlsnp_LD[rsID]
    else:
        # using the tabix command available from the system
        # this is deprecated, use pytabix instead
        datafile=data_folder_location+'LD/'+chrom+'.LD85.tab.gz'
        position_str=chrom+":"+str(pos)+"-"+str(pos)
        proc = Popen(['tabix', datafile, position_str], stdout=PIPE)
        result = proc.stdout.read()
        lines = result.split('\n')
        qr = []
        for line in lines:
            items = line.split('\t')
            if len(items)>=6:
                # [snp, rsq, pos]
                qr.append([items[4], items[5], items[3]])

        #print "tabix1 find:" + str(len(qr)) +' lines'

        for result in qr:
            # [snp, rsq, distance_to_rsID]
            SNPs_in_LD.append([result[0], float(result[1]),
                abs(int(result[2]) - int(pos))])

        # Sort to make for an efficient search of the GWAS database
        #(i.e. use the first one that we find)
        # Sort by proximity (want closest first)
        SNPs_in_LD.sort(key=lambda x: x[2], reverse=False)
        # an alternative way will be sorted by LD (highest first
        # SNPs_in_LD.sort(key=lambda x: x[1], reverse=True)
        #print rsID+ '\tSNPs_in_LD:\t' + str(SNPs_in_LD)
        # save the result SNPs_in_LD to eqtlsnp_LD
        eqtlsnp_LD[rsID]=SNPs_in_LD

    print ('tabix\t' + rsID+ '\tSNPs_in_LD:\t' + str(SNPs_in_LD))

    # Next, search the GWAS results until we either find a SNP, or don't
    output = 'none'
    for item in SNPs_in_LD:
        matching_gwas_snp = get_matching_GWAS_SNP_dict(gwas_snp, item[0])
        if matching_gwas_snp != 'none':
            output = matching_gwas_snp
            #print 'tabix output before append:'+str(output)
            output.append(item[1])
            output.append(item[2])
            #print 'tabix output after append:'+str(output)
            #print str(output)
            break

    return output

def get_high_LD_GWAS_SNPs_pytabix(GWAS, rsID, chrom, pos, LD_threshold,
        data_folder_location, gwas_snp, eqtlsnp_LD2, tabixQuery):
    """Get GWAS snps in high LD using pytabix,
    preloaded tabixQuery (one tabix query for a merged LD85.tab.gz)
    """
    # First, create a list of ALL SNPs in LD with the eQTL SNP at better
    # than the LD_rsq threshold
    SNPs_in_LD = []
    if rsID in eqtlsnp_LD2:
        #print "already known the SNPs_in_LD for " + rsID
        SNPs_in_LD = eqtlsnp_LD2[rsID]
    else:
        position_str=chrom+":"+str(pos)+"-"+str(pos)
        records = []
        try:
            records = tabixQuery.querys(position_str) #query(chrom,pos,pos)
        except:
            print ("error query with tabix: ", position_str)

        qr = []
        for record in records:
            #print record
            if len(record)>=6:
                qr.append([record[4], record[5], record[3]]) # [snpID,LD,pos]

        if qr != ():
            for result in qr:
                SNPs_in_LD.append([result[0], float(result[1]),
                    abs(int(result[2]) - int(pos))]) # [snpID, LD, distance]

        # Sort to make for an efficient search of the GWAS database
        # (i.e. use the first one that we find)
        # Sort by distance (proximity, want closest first)
        SNPs_in_LD.sort(key=lambda x: x[2], reverse=False)
        # alternative way would be to sort by LD, (want highest first)
        # SNPs_in_LD.sort(key=lambda x: x[1], reverse=True) # sorted by LD
        #print rsID+ '\tSNPs_in_LD:\t' + str(SNPs_in_LD)
        eqtlsnp_LD2[rsID]=SNPs_in_LD

    #print 'tabix2\t' + rsID+ '\tSNPs_in_LD:\t' + str(SNPs_in_LD)
    # Next, search the GWAS results until we either find a SNP, or don't
    output = 'none'
    if len(SNPs_in_LD) > 0:
        for item in SNPs_in_LD:
            matching_gwas_snp = get_matching_GWAS_SNP_dict(gwas_snp, item[0])
            if matching_gwas_snp != 'none':
                output = matching_gwas_snp
                output.append(item[1])
                output.append(item[2])
                # output would be [[snp, pval, chrom, location, maf, LD, dist]
                #print str(output)
                break

    return output

def get_high_LD_GWAS_SNPs_preload(GWAS, rsID, chrom, pos, align_LD_threshold,
        gwas_snp, eqtlsnp_LD2, eSNP_POS85, eSNP_linked85):
    """Get GWAS snps in high LD using preload eSNP_LD85
    """
    # First, create a list of ALL SNPs in LD with the eQTL SNP at better
    # than the LD_rsq threshold
    SNPs_in_LD = []
    # chromosomal pos for rsID
    snpA_pos = int(pos)
    if rsID in eqtlsnp_LD2:
        #print "already known the SNPs_in_LD for " + rsID
        SNPs_in_LD = eqtlsnp_LD2[rsID]
    else:
        if rsID in eSNP_linked85:
            for snpB, ld in eSNP_linked85[rsID]:
                snpB_pos = int(eSNP_POS85[snpB][1]) #eSNP_POS85[snpB]=[chr,pos]
                # [snpID, LD, dist]
                SNPs_in_LD.append([snpB, ld, abs(snpB_pos - snpA_pos)])
        #else:
        #    print "error query with preload: ", rsID, chrom, pos
        # Sort to make for an efficient search of the GWAS database
        # (i.e. use the first one that we find)
        # Sort by distance (proximity, want closest first)
        SNPs_in_LD.sort(key=lambda x: x[2], reverse=False)
        # alternative way would be to sort by LD, (want highest first)
        # SNPs_in_LD.sort(key=lambda x: x[1], reverse=True) # sorted by LD
        #print rsID+ '\tSNPs_in_LD:\t' + str(SNPs_in_LD)
        eqtlsnp_LD2[rsID]=SNPs_in_LD

    #print 'tabix2\t' + rsID+ '\tSNPs_in_LD:\t' + str(SNPs_in_LD)
    # Next, search the GWAS results until we either find a SNP, or don't
    matching_gwas_snp = 'none'
    for item in SNPs_in_LD:
        if item[0] in gwas_snp:
            #array in order of snp,pval,chrom,location,maf
            matching_gwas_snp = gwas_snp[item[0]][0:5]
            #matching_gwas_snp = get_matching_GWAS_SNP_dict(gwas_snp, item[0])
            matching_gwas_snp.append(item[1])
            matching_gwas_snp.append(item[2])
            # [snp, pval, chrom, location, maf, LD, dist]
            break
    #for a small number of gwas snps the we need to match the snp_location
    if matching_gwas_snp == 'none':
        # match the gwas positions instead
        pass

    return matching_gwas_snp


def get_high_LD_GWAS_SNPs(GWAS, rsID, database_cursor, LD_threshold):
    """find high LD gwas snp -- Database version
    """
    database_cursor.execute("use hg19")
    # First, create a list of ALL SNPs in LD with the eQTL SNP at better than the LD_rsq threshold
    SQL_string = "SELECT SNPA, rsq, SNPA_loc, SNPB_loc from LD_" + SNP_pair.chrom + "  WHERE SNPB='" + SNP_pair.eQTL_SNP + "'" + " AND rsq > " + str(LD_threshold)
    database_cursor.execute(SQL_string)
    qr1 = database_cursor.fetchall()
    #print qr1
    SQL_string = "SELECT SNPB, rsq, SNPB_loc, SNPA_loc from LD_" + SNP_pair.chrom + "  WHERE SNPA='" + SNP_pair.eQTL_SNP + "'" + " AND rsq > " + str(LD_threshold)
    database_cursor.execute(SQL_string)
    qr2 = database_cursor.fetchall()
    print ('Searching: '+rsID+"\t"+str(SNP_pair.eQTL_SNP_location))
    #print qr2
    qr = qr1 + qr2
    SNPs_in_LD = []
    if qr != ():
        for result in qr:
            SNPs_in_LD.append([result[0], float(result[1]), abs(int(result[2]) - int(SNP_pair.eQTL_SNP_location))])
    if len(SNPs_in_LD)>0:
        print ('Searching for:'+ rsID+"\t"+str(SNP_pair.eQTL_SNP_location)+"\t"+str(qr[0][3]))

    # Sort to make for an efficient search of the GWAS database (i.e. use the first one that we find)
    SNPs_in_LD.sort(key=lambda x: x[2], reverse=False)     # Sort first by proximity (want closest first)
    print ('mysql\t' + rsID + '\t' + 'SNPs_in_LD:\t' + str(SNPs_in_LD))
    #SNPs_in_LD.sort(key=lambda x: x[1], reverse=True)     # Sort next by LD (want highest first)


    # Next, search the GWAS results until we either find a SNP, or don't
    output = 'none'
    if len(SNPs_in_LD) > 0:
        for item in SNPs_in_LD:
            matching_gwas_snp = get_matching_GWAS_SNP(GWAS, item[0], database_cursor)
            if matching_gwas_snp != 'none':
                output = matching_gwas_snp
                print ('mysql output before append:'+str(output))
                output.append(item[1])
                output.append(item[2])
                print ('mysql output after append:'+str(output))
                break

    return output

def get_proximal_GWAS_SNPs_tabix(GWAS, gwas_tabix, chrom, location, eQTL_SNP_MAF, proximal, MAF_match):

    # Ensure that eQTL_SNP_MAF is a float
    try:
        maf_value = float(eQTL_SNP_MAF)
    except:
        eQTL_SNP_MAF = float(0.0001)

    #print location
    #print proximal
    upper_bound = location + proximal
    lower_bound = location - proximal
    # SQL_string = "SELECT snp, pval, chrom, location, maf from " + GWAS + " WHERE chrom='" + chrom + "' AND location > " + str(lower_bound) + " AND location < " + str(upper_bound) + ";"
    position_str=chrom+":"+str(lower_bound)+"-"+str(upper_bound)
    #print position_str
    records = []
    try:
        records = gwas_tabix.querys(position_str) #query(chrom,pos,pos)
    except:
        print ("tabix query error for proximal GWAS snp:", position_str)
    #print "success query tabix"
    qr = []
    for record in records:
        #print record
        if len(record)>=5:
            qr.append([record[2], record[3], record[0], record[1], record[4]])

        #print "tabix2 find:" + str(len(qr)) +' lines'

    Proximal_SNPs = []
    if qr != ():
        #print qr
        for result in qr:
            tmp_maf = result[4]
            try:
                foo = float(tmp_maf)
                if foo == float(0):
                    foo = float(0.0001)
            except:
                foo = float(0.0001)

            Proximal_SNPs.append([result[0], result[1], result[2], result[3], foo, abs(int(result[3]) - int(location))])

    # Sort to make
    Proximal_SNPs.sort(key=lambda x: abs(float(eQTL_SNP_MAF) - float(x[4]))/(x[4]), reverse=False)     # Sort by percentage difference in MAF
    #print Proximal_SNPs

    # Next, search the GWAS results until we either find a SNP, or don't
    output = 'none'
    if len(Proximal_SNPs)>0 and float(eQTL_SNP_MAF)>0 :
        if abs(float(eQTL_SNP_MAF) - float(Proximal_SNPs[0][4]))/float(eQTL_SNP_MAF) < 0.2:
            output = Proximal_SNPs[0]

    return output


def get_proximal_GWAS_SNPs(GWAS, database_cursor, chrom, location, eQTL_SNP_MAF, proximal, MAF_match):

    # Ensure that eQTL_SNP_MAF is a float
    try:
        boo = float(eQTL_SNP_MAF)
    except:
        eQTL_SNP_MAF = float(0.0001)

    database_cursor.execute("use GWAS")
    #print location
    #print proximal
    upper_bound = location + proximal
    lower_bound = location - proximal
    SQL_string = "SELECT snp, pval, chrom, location, maf from " + GWAS + " WHERE chrom='" + chrom + "' AND location > " + str(lower_bound) + " AND location < " + str(upper_bound) + ";"
    #print SQL_string
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    Proximal_SNPs = []
    if qr != ():
        #print qr
        for result in qr:
            tmp_maf = result[4]
            try:
                foo = float(tmp_maf)
                if foo == float(0):
                    foo = float(0.0001)
            except:
                foo = float(0.0001)

            Proximal_SNPs.append([result[0], result[1], result[2], result[3], foo, abs(int(result[3]) - int(location))])

    # Sort to make
    Proximal_SNPs.sort(key=lambda x: abs(float(eQTL_SNP_MAF) - float(x[4]))/(x[4]), reverse=False)     # Sort by percentage difference in MAF
    #print Proximal_SNPs

    # Next, search the GWAS results until we either find a SNP, or don't
    output = 'none'
    if len(Proximal_SNPs) > 0:
        if abs(float(eQTL_SNP_MAF) - float(Proximal_SNPs[0][4]))/float(eQTL_SNP_MAF) < 0.2:
            output = Proximal_SNPs[0]

    return output



# Use the full list of eQTL tag SNPs to match GWAS SNPs. Unlike previous versions, we don't discard any SNPs -- simply tag them.
def match_eQTL_SNPs_with_GWAS(eQTL, gene_list, database_cursor, align_LD_threshold,
        GWAS, data_folder_location):
    """Alignment of eQTL vs GWAS
    """
    #prepare tabix
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    #url=data_folder_location+'LD/LD85.tab.gz'
    #url = '/run/shm/LD/LD85.tab.gz'
    #tabixQuery=tabix.open(url)
    #pdb.set_trace()
    eqtlsnp_LD={}
    eqtlsnp_LD2={}
    gwas_file=data_folder_location+"GWAS/"+GWAS+".txt"
    gwas_tabix_file=data_folder_location+"GWAS/"+GWAS+".tab.gz"
    gwas_snp=get_all_GWAS_SNP_local(GWAS,gwas_file)
    #create another dict with chr:pos to search for matching gwas_snps
    gwas_snp_pos = {}
    for snp_id in gwas_snp:
        chr_pos = gwas_snp[snp_id][2]+":"+str(gwas_snp[snp_id][3])
        gwas_snp_pos[chr_pos] = gwas_snp[snp_id]

    gwas_tabix=tabix.open(gwas_tabix_file)

    eQTL_snps = read_eQTL_snps(eQTL, data_folder_location)
    #preload LD85 data and filter
    time_start = time.time()
    #eSNP_POS85, eSNP_linked85 = get_eSNP_in_LD85_filter(eQTL,
    #        data_folder_location, gwas_snp, eQTL_snps)
    eSNP_POS85, eSNP_linked85 = get_eSNP_in_LD85_filter(eQTL,
            data_folder_location, gwas_snp, eQTL_snps,align_LD_threshold)

    print ('Done preload LD85. Took:',str(time.time()-time_start))
    #eSNP_POS85, eSNP_linked85 = get_eSNP_in_LD85(eQTL,
    #        data_folder_location)

    #gwas_snp_list,gwas_snp=get_all_GWAS_SNP(GWAS, database_cursor)

    t0 = time.time()
    t1 = t0
    for gene_index in range(len(gene_list)):
        time_tmp = time.time()
        time_for_gene = time_tmp-t1
        time_accum = time_tmp-t0
        t1 = time_tmp
        print ('Matching gene: ' + gene_list[gene_index].name + 'time total:'+str(time_accum) + ' time previous:'+str(time_for_gene))
        SNP_pairs = gene_list[gene_index].SNP_pairs
        SNP_pairs_out = []
        for SNP_pair in SNP_pairs:
            # First find exact matches
            #output = get_matching_GWAS_SNP(GWAS, SNP_pair.eQTL_SNP,
            #         database_cursor)
            #print str(SNP_pair.eQTL_SNP)
            output = get_matching_GWAS_SNP_dict(gwas_snp, SNP_pair.eQTL_SNP)

            if output != 'none':
                #print "matching step 1 success"
                SNP_pair.GWAS_SNP = output[0]
                SNP_pair.GWAS_SNP_pval = output[1]
                SNP_pair.GWAS_SNP_chrom = output[2]
                SNP_pair.GWAS_SNP_location = output[3]
                SNP_pair.GWAS_SNP_MAF = output[4]
                SNP_pair.LD_rsq = 1.000
                SNP_pair.separation = 0
                SNP_pairs_out.append(SNP_pair)

            # Then find matching SNPs using high LD
            else:
                #t01=time.time()
                #output = get_high_LD_GWAS_SNPs(GWAS, SNP_pair.eQTL_SNP, database_cursor, LD_threshold)
                #print "mysql search took "+str(time.time()-t01)
                #t11=time.time()
                #output1=get_high_LD_GWAS_SNPs_tabix(GWAS, SNP_pair.eQTL_SNP, SNP_pair.chrom, SNP_pair.eQTL_SNP_location, database_cursor, LD_threshold, data_folder_location, gwas_snp, eqtlsnp_LD)
                #print "tabix search 1 took:" + str(time.time()-t11)
                #t21=time.time()
                #output1=get_high_LD_GWAS_SNPs_pytabix(GWAS, SNP_pair.eQTL_SNP,
                #        SNP_pair.chrom, SNP_pair.eQTL_SNP_location,
                #        LD_threshold, data_folder_location, gwas_snp,
                #        eqtlsnp_LD2, tabixQuery)
                output=get_high_LD_GWAS_SNPs_preload(GWAS, SNP_pair.eQTL_SNP,
                        SNP_pair.chrom, SNP_pair.eQTL_SNP_location,
                        align_LD_threshold, gwas_snp,
                        eqtlsnp_LD2, eSNP_POS85, eSNP_linked85)

                #print "tabix search 2 took:" + str(time.time()-t21)


                #output1 = get_high_LD_GWAS_SNPs_tabix(GWAS, SNP_pair.eQTL_SNP, SNP_pair.chrom, SNP_pair.eQTL_SNP_location, database_cursor, LD_threshold, data_folder_location, gwas_snp, eqtlsnp_LD)
                #diff=0
                #if output != 'none':
                #    if output1 !='none':
                #        for outputIndex in range(len(output)):
                #            if not (output[outputIndex]==output1[outputIndex]):
                #                diff=1
                #                break
                #    else :
                #        diff=1
                #else:
                #    if output1 != 'none' :
                #        diff=1
                #print "output\t"+str(output)
                #print "output1\t"+str(output1)
                #if diff==1:
                #    print "ERROR: output not matched: " + SNP_pair.eQTL_SNP
                #else:
                #    print "Correct: output matched"

                #using tabix result
                if output != 'none':
                    #print "matching step 2 succeed"
                    SNP_pair.GWAS_SNP = output[0]
                    SNP_pair.GWAS_SNP_pval = output[1]
                    SNP_pair.GWAS_SNP_chrom = output[2]
                    SNP_pair.GWAS_SNP_location = output[3]
                    SNP_pair.GWAS_SNP_MAF = output[4]
                    SNP_pair.LD_rsq = output[5]
                    SNP_pair.separation = output[6]
                    SNP_pairs_out.append(SNP_pair)

            # If the above (preferred) matching appemps fail, find matches using very close proximity and similar MAF for those SNPs not in our LD database
            if output == 'none':
                output = get_proximal_GWAS_SNPs_tabix(GWAS, gwas_tabix, SNP_pair.chrom, SNP_pair.eQTL_SNP_location, SNP_pair.eQTL_SNP_MAF, proximal=250, MAF_match=0.02)
                if output != 'none':
                    #print "matching step 3 succeed"
                    SNP_pair.GWAS_SNP = output[0]
                    SNP_pair.GWAS_SNP_pval = output[1]
                    SNP_pair.GWAS_SNP_chrom = output[2]
                    SNP_pair.GWAS_SNP_location = output[3]
                    SNP_pair.GWAS_SNP_MAF = output[4]
                    SNP_pair.LD_rsq = 'NA'
                    SNP_pair.separation = output[5]
                    SNP_pairs_out.append(SNP_pair)

                # Finally, update information on the other SNPs to indicate NF: None Found
                else:
                    #print "matching fail"
                    SNP_pair.GWAS_SNP = 'NF'
                    SNP_pair.GWAS_SNP_pval = 'NA'
                    SNP_pair.GWAS_SNP_chrom = 'NA'
                    SNP_pair.GWAS_SNP_location = 'NA'
                    SNP_pair.GWAS_SNP_MAF = 'NA'
                    SNP_pair.LD_rsq = 'NA'
                    SNP_pair.separation = 'NA'
                    SNP_pairs_out.append(SNP_pair)


        # Replace the SNP_pairs with updated SNP_pairs for each gene
        gene_list[gene_index].SNP_pairs = SNP_pairs_out

    print ('End match_eQTL_SNPs_with_GWAS')

    return gene_list, eqtlsnp_LD2

# Determine which SNPs from matched SNP pairs should be used as tag SNPs.  These are the only SNPs that factor into our analysis
def assign_eQTL_tag_SNPs(gene_list, database_cursor, LD_threshold):
    # Extract eQTL SNPs for each gene and create blocks
    print ('assign_eQTL_tag_SNPs called')
    database_cursor.execute("use hg19")

    #chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

    for gene in gene_list:

        # Identify gene
        #print gene.name

        # Build blocks on a chromosome-by-chromosome basis (obviously!!)
        for chrom in chroms:

            # First make blocks using SNPs that are in the LD_DB and are on the same chromosome.
            #print gene.SNP_pairs
            SNPs_tmp = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.in_LD_DB == 'Yes' and foo.GWAS_SNP != 'NF']
            SNPs = apply_eQTL_threshold(SNPs_tmp)

             # Store linked SNPs for each chrom as a list of lists representing the block structure
            chrom_blocks = []
            for item in [str(foo.eQTL_SNP) for foo in SNPs]:
                #print 'Working SNP: ' + str(count)
                #print item
                SQL_string = "SELECT SNPB, rsq FROM LD_" + chrom + " WHERE SNPA=" + "'" + item + "'" + "AND rsq >= " + str(LD_threshold) + " AND SNPB IN (" + ', '.join(["'" + foo + "'" for foo in [str(moo.eQTL_SNP) for moo in SNPs]]) + ")"
                database_cursor.execute(SQL_string)
                qr1 = database_cursor.fetchall()
                SQL_string = "SELECT SNPA, rsq FROM LD_" + chrom + " WHERE SNPB=" + "'" + item + "'" + " AND rsq >= " + str(LD_threshold) + " AND SNPA IN (" + ', '.join(["'" + foo + "'" for foo in [str(moo.eQTL_SNP) for moo in SNPs]]) + ")"
                database_cursor.execute(SQL_string)
                qr2 = database_cursor.fetchall()

                # Identify the SNPs in high LD that have the strongest eQTL
                qr = qr1 + qr2

                #print "LD for "+item
                #print qr

                linked_SNPs = [[item, 1.00]]
                if qr != ():
                    for result in qr:
                        linked_SNPs.append([result[0], result[1]])

                # Order the SNPs in each block to ensure equivalent resutls when the eQTL p-values happen to be identical
                linked_SNPs.sort(key=lambda x: x[0])

                # Add anything we've found to the block list for this chrom
                if linked_SNPs != []:
                    chrom_blocks.append(linked_SNPs)


            # Process the blocks for each chrom
            for linked_SNPs in chrom_blocks:

                # Get the eQTL SNP p-value and MAF for these results
                comparison_list = []
                for linked_SNP in linked_SNPs:
                    #print 'linked_SNP: '+str(linked_SNP)
                    #for snp1 in gene.SNP_pairs:
                    #    print snp1.eQTL_SNP+'\t'+snp1.eQTL_SNP_pval
                    foo = next(float(foo.eQTL_SNP_pval) for foo in gene.SNP_pairs if foo.eQTL_SNP == linked_SNP[0])
                    #print 'foo:'+ str(foo)
                    comparison_list.append([linked_SNP[0], linked_SNP[1], foo])   # These are rs#, rsq, and the eQTL p-value

                    # Want to select the strongest (smallest) eQTL p-value for all the SNPs that are linked.
                    comparison_list.sort(key=lambda x: x[2], reverse=False)
                    #print comparison_list
                    tag_SNP = comparison_list[0][0]

                # Change the tag SNP flag
                for pair in gene.SNP_pairs:
                    if pair.eQTL_SNP == tag_SNP:
                        pair.is_tag_SNP = 'Yes'

                # To debug this, print everything
                #print
                #print comparison_list
                #print 'tag_SNP: ' + tag_SNP


        # Now verify that the tagged genes are not in LD with each other
        for chrom in chroms:

            # The exclude list ensures that we don't remove tag SNPs due to long-distance chaining of SNP blocks
            exclude_list = []

            SNPs = [foo.eQTL_SNP for foo in gene.SNP_pairs if foo.chrom == chrom and foo.is_tag_SNP == 'Yes']

            # Store linked SNPs for each chrom as a list of lists representing the block structure
            for item in SNPs:
                if item not in exclude_list:
                    #print 'Working SNP: ' + str(count)
                    SQL_string = "SELECT SNPB, rsq FROM LD_" + chrom + " WHERE SNPA=" + "'" + item + "'" + " AND rsq >= 0.2 AND SNPB IN (" + ', '.join(["'" + foo + "'" for foo in SNPs]) + ")"
                    database_cursor.execute(SQL_string)
                    qr1 = database_cursor.fetchall()
                    SQL_string = "SELECT SNPA, rsq FROM LD_" + chrom + " WHERE SNPB=" + "'" + item + "'" + " AND rsq >= 0.2 AND SNPA IN (" + ', '.join(["'" + foo + "'" for foo in SNPs]) + ")"
                    database_cursor.execute(SQL_string)
                    qr2 = database_cursor.fetchall()

                    # Identify the SNPs in high LD that have the strongest eQTL
                    qr = qr1 + qr2
                    #print qr

                    linked_SNPs = []
                    if qr != ():
                        # Some of the SNPs are still in LD, we will need to remove the extras.
                        for result in qr:
                            linked_SNPs.append(result[0])
                        linked_SNPs.sort()
                        for boo in linked_SNPs:
                            exclude_list.append(boo)

                    # If we find other tag SNPs in LD, we need to check their eQTL p-values and remove the large ones
                    if len(linked_SNPs) > 0:
                        # Get the eQTL SNP p-value and MAF for these results
                        comparison_list = []
                        for linked_SNP in linked_SNPs:
                            foo = next(float(foo.eQTL_SNP_pval) for foo in gene.SNP_pairs if foo.eQTL_SNP == linked_SNP)
                            comparison_list.append([linked_SNP, foo])   # These are rs# and the eQTL p-value

                            # Want to select the strongest (smallest) eQTL p-value for all the SNPs that are linked.
                            comparison_list.sort(key=lambda x: x[1], reverse=False)
                            #print comparison_list
                            tag_SNP = comparison_list[0][0]

                        # Change the tag SNP flag to No for SNPs with the larger p-values
                        for pair in gene.SNP_pairs:
                            if pair.eQTL_SNP != tag_SNP:
                                if pair.eQTL_SNP in linked_SNPs:
                                    pair.is_tag_SNP = 'No'

        # Now deal with the regions, mostly HLA, that do not have good coverage in our LD database.  Use proximity and MAF to identify tag SNPs
        for chrom in chroms:

            # The exclude list ensures that we don't remove tag SNPs due to long-distance chaining of SNP blocks
            exclude_list = []

            # Get the list of SNPs not in the LD database that have corresponding GWAS
            SNPs_tmp = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.in_LD_DB == 'No' and foo.GWAS_SNP_pval != 'NA']
            SNPs = apply_eQTL_threshold(SNPs_tmp)

            for item in [str(foo.eQTL_SNP) for foo in SNPs]:
                if item not in exclude_list:
                    item_location = next(foo.eQTL_SNP_location for foo in gene.SNP_pairs if foo.eQTL_SNP == item)
                    upper_bound = int(item_location) + 5000
                    lower_bound = int(item_location) - 5000
                    proximal_SNPs_tmp = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.in_LD_DB == 'No' and foo.GWAS_SNP_pval != 'NA' and foo.eQTL_SNP_location > lower_bound and foo.eQTL_SNP_location < upper_bound]
                    proximal_SNPs = [str(foo.eQTL_SNP) for foo in apply_eQTL_threshold(proximal_SNPs_tmp)]

                    proximal_SNPs.sort()
                    for boo in proximal_SNPs:
                        exclude_list.append(boo)

                    # Get the eQTL SNP p-value and MAF for these results (will always included at least the SNP under test (aka item))
                    comparison_list = []
                    for proximal_SNP in proximal_SNPs:
                        foo = next(float(foo.eQTL_SNP_pval) for foo in gene.SNP_pairs if foo.eQTL_SNP == proximal_SNP)
                        comparison_list.append([proximal_SNP, foo])   # These are rs# and the eQTL p-value

                        # Want to select the strongest (smallest) eQTL p-value for all the SNPs that are linked.
                        comparison_list.sort(key=lambda x: x[1], reverse=False)
                        #print comparison_list
                        tag_SNP = comparison_list[0][0]

                    #print 'Tag SNP is: ' + str(tag_SNP)
                    # Change the tag SNP flag to No for SNPs with the larger p-values
                    for pair in gene.SNP_pairs:
                        if pair.eQTL_SNP == tag_SNP:
                                pair.is_tag_SNP = 'Yes'


 # Now verify that SNPs aren't in close proximity
        min_separation = 100000
        for chrom in chroms:

            untag_snp = 'start'
            while untag_snp != 'none':
                SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.is_tag_SNP == 'Yes']
                #print SNPs
                untag_snp = 'none'
                for i in range(1,len(SNPs)):
                    if int(abs(SNPs[i-1].eQTL_SNP_location - SNPs[i].eQTL_SNP_location)) < int(min_separation):
                        if float(SNPs[i-1].eQTL_SNP_pval) <= float(SNPs[i].eQTL_SNP_pval):
                            untag_snp = SNPs[i].eQTL_SNP
                        else:
                            untag_snp = SNPs[i-1].eQTL_SNP
                        # Now update
                        for pair in gene.SNP_pairs:
                            if pair.eQTL_SNP == untag_snp:
                                #print 'Untaging SNP: ' + str(untag_snp)
                                pair.is_tag_SNP = 'No'

                        # break out of the for loop and repeate with updated gene
                        break


 # Now verify that adjacent SNPs have rather different MAFs
        min_separation = 150000
        min_MAF_diff = 0.1
        for chrom in chroms:

            untag_snp = 'start'
            while untag_snp != 'none':
                SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.is_tag_SNP == 'Yes']
                #print SNPs
                untag_snp = 'none'
                for i in range(1,len(SNPs)):
                    try:
                        if int(abs(SNPs[i-1].eQTL_SNP_location - SNPs[i].eQTL_SNP_location)) < int(min_separation) and (abs(float(SNPs[i-1].eQTL_SNP_MAF) - float(SNPs[i].eQTL_SNP_MAF)) < float(min_MAF_diff)):
                            #print "Untag for gene: " + str(gene.name)
                            if float(SNPs[i-1].eQTL_SNP_pval) <= float(SNPs[i].eQTL_SNP_pval):
                                untag_snp = SNPs[i].eQTL_SNP
                            else:
                                untag_snp = SNPs[i-1].eQTL_SNP
                            # Now update
                            for pair in gene.SNP_pairs:
                                if pair.eQTL_SNP == untag_snp:
                                    #print 'Untaging SNP: ' + str(untag_snp)
                                    pair.is_tag_SNP = 'No'

                            # break out of the for loop and repeate with updated gene
                            break
                    except:
                        # In some cases, there is no MAF in the database, so we just skip this step
                        pass


 # Do not tag chromosome X
        SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == 'chrX' and foo.is_tag_SNP == 'Yes']
        #print SNPs
        for i in range(len(SNPs)):
            for pair in gene.SNP_pairs:
                if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
                    pair.is_tag_SNP = 'No'

# Do not tag the HLA region
        SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == 'chr6' and foo.is_tag_SNP == 'Yes' and (int(foo.eQTL_SNP_location) > int(27000000) and int(foo.eQTL_SNP_location) < int(33000000))]
        #print SNPs
        for i in range(len(SNPs)):
            for pair in gene.SNP_pairs:
                if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
                    pair.is_tag_SNP = 'No'

 ## Set a maximum eQTL p-value threshold for trans SNPs
 #       SNPs = [foo for foo in gene.SNP_pairs if (float(foo.eQTL_SNP_pval) >= float(3e-6)) and foo.is_tag_SNP == 'Yes' and foo.is_cis_SNP == 'trans']
 #       #print SNPs
 #       for i in range(len(SNPs)):
 #           for pair in gene.SNP_pairs:
 #               if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
 #                   pair.is_tag_SNP = 'No'

 # Set a maximum eQTL p-value threshold for cis SNPs
        SNPs = [foo for foo in gene.SNP_pairs if (float(foo.eQTL_SNP_pval) >= float(1e-3)) and foo.is_tag_SNP == 'Yes' and foo.is_cis_SNP == 'cis']
        #print SNPs
        for i in range(len(SNPs)):
            for pair in gene.SNP_pairs:
                if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
                    pair.is_tag_SNP = 'No'


    return gene_list


# Determine which SNPs from matched SNP pairs should be used as tag SNPs.  These are the only SNPs that factor into our analysis
def assign_eQTL_tag_SNPs_tabix(gene_list, eQTLSNP_LD, LD_threshold,
        data_folder_location):
    # Extract eQTL SNPs for each gene and create blocks
    print ('assign_eQTL_tag_SNPs called')

    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
            'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
            'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
            'chr22', 'chrX', 'chrY']
    eQTL_dict={}
    #prepare tabix for LD 0.2, used merged LD02.tab.gz
    #url=data_folder_location+'LD/LD02.tab.gz'
    #testing use /run/shm/LD/LD02.tab.gz'
    url = '/run/shm/LD/LD02.tab.gz'
    tabix_LD_02 = tabix.open(url)

    for gene in gene_list:

        # Identify gene
        print (gene.name)

        # Build blocks on a chromosome-by-chromosome basis (obviously!!)
        for chrom in chroms:

            # First make blocks using SNPs that are in the LD_DB and are on the same chromosome.
            #print gene.SNP_pairs
            SNPs_tmp = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.in_LD_DB == 'Yes' and foo.GWAS_SNP != 'NF']
            SNPs = apply_eQTL_threshold(SNPs_tmp)

            # Store linked SNPs for each chrom as a list of lists
            # representing the block structure
            chrom_blocks = []

            for item in [str(foo.eQTL_SNP) for foo in SNPs]:
                # take advantage of the fact that eQTLSNP_LD has been
                # constructed from match_eQTL_SNPs_with_GWAS
                qr=[]
                if item in eQTLSNP_LD:
                    if not (item in eQTL_dict):
                        eQTL_dict[item]={}
                        for qritem in eQTLSNP_LD[item]:
                            eQTL_dict[item][qritem[0]]=qritem[1]
                    for eqtlsnp in [str(foo.eQTL_SNP) for foo in SNPs]:
                        if eqtlsnp in eQTL_dict[item]:
                            #print "found "+qritem[0]
                            qr.append([eqtlsnp,eQTL_dict[item][eqtlsnp]])
                        #else:
                        #    print qritem[0]+" is not in SNPDict"
                # Identify the SNPs in high LD that have the strongest eQTL
                #print "LD for "+item
                #print qr
                linked_SNPs = [[item, 1.00]]
                if qr != []:
                    for result in qr:
                        linked_SNPs.append([result[0], result[1]])

                # Order the SNPs in each block to ensure equivalent
                # resutls when the eQTL p-values happen to be identical
                linked_SNPs.sort(key=lambda x: x[0])

                # Add anything we've found to the block list for this chrom
                chrom_blocks.append(linked_SNPs)


            # Process the blocks for each chrom
            for linked_SNPs in chrom_blocks:

                # Get the eQTL SNP p-value and MAF for these results
                comparison_list = []
                for linked_SNP in linked_SNPs:
                    #print "linked_SNP:"+str(linked_SNP)
                    #for snp in gene.SNP_pairs:
                    #    print snp.eQTL_SNP+"\t"+snp.eQTL_SNP_pval
                    foo = next(float(foo.eQTL_SNP_pval) for foo in gene.SNP_pairs if foo.eQTL_SNP == linked_SNP[0])
                    #$print "foo:"+str(foo)
                    comparison_list.append([linked_SNP[0], linked_SNP[1], foo])   # These are rs#, rsq, and the eQTL p-value

                # Want to select the strongest (smallest) eQTL p-value for all the SNPs that are linked.
                comparison_list.sort(key=lambda x: x[2], reverse=False)
                #print comparison_list
                tag_SNP = comparison_list[0][0]

                # Change the tag SNP flag
                for pair in gene.SNP_pairs:
                    if pair.eQTL_SNP == tag_SNP:
                        pair.is_tag_SNP = 'Yes'

                # To debug this, print everything
                #print
                #print comparison_list
                #print 'tag_SNP: ' + tag_SNP


        # Now verify that the tagged genes are not in LD with each other
        for chrom in chroms:

            # The exclude list ensures that we don't remove tag SNPs due to long-distance chaining of SNP blocks
            exclude_list = {}

            #SNPs = [foo.eQTL_SNP for foo in gene.SNP_pairs if foo.chrom == chrom and foo.is_tag_SNP == 'Yes']
            # instead of creating a list, use a dictionary
            # keep snps in the original order
            sortedSNPs = []
            # SNPValues = {'snpid':[location, pval]}
            SNPValues = {}
            for item in gene.SNP_pairs:
                if item.is_tag_SNP == 'Yes' and item.chrom == chrom:
                    SNPValues[item.eQTL_SNP] = [item.eQTL_SNP_location,
                            item.eQTL_SNP_pval]
                    sortedSNPs.append(item.eQTL_SNP)

            # Store linked SNPs for each chrom as a list of lists representing the block structure
            for SNP in sortedSNPs:
                #print "checking",SNP
                #print "Excluded List",exclude_list
                if SNP not in exclude_list:
                    # SNPValues[SNP] is [pos, pval]
                    snp_location = str(SNPValues[SNP][0])
                    position_str = chrom + ":" + snp_location + \
                                           "-" + snp_location
                    records = []
                    try:
                        #query(chrom,pos,pos)
                        records = tabix_LD_02.querys(position_str)
                    except:
                        print ("error query with tabix: ", position_str)

                    qr = []
                    for record in records:
                        # [chrom,pos1,snp1,pos2,snp2,ld]
                        if len(record)>=6 and record[4] in SNPValues:
                            qr.append([record[4], record[5]])

                    # Identify the SNPs in high LD that have the strongest eQTL
                    linked_SNPs = [SNP]
                    if qr:
                        # Some of the SNPs are still in LD,
                        # we will need to remove the extras.
                        for result in qr:
                            linked_SNPs.append(result[0])
                        #linked_SNPs.sort()
                        for snp_B in linked_SNPs:
                            exclude_list.setdefault(snp_B,1)

                    # If we find other tag SNPs in LD, we need to check their eQTL p-values and remove the large ones
                    if len(linked_SNPs) > 0:
                        #print 'linked snp', linked_SNPs
                        # Get the eQTL SNP p-value and MAF for these results
                        comparison_list = [] #[SNP, float(SNPValues[SNP][1])]]
                        for linked_SNP in linked_SNPs:
                            #snp_pval = next(float(foo.eQTL_SNP_pval) for foo in gene.SNP_pairs if foo.eQTL_SNP == linked_SNP)
                            snp_pval = float(SNPValues[linked_SNP][1])
                            comparison_list.append([linked_SNP, snp_pval])   # These are rs# and the eQTL p-value

                        # Want to select the strongest (smallest) eQTL p-value for all the SNPs that are linked.
                        comparison_list.sort(key=lambda x: x[1], reverse=False)
                        #print 'comparison_list', comparison_list
                        tag_SNP = comparison_list[0][0]
                        #print 'set tag_SNP', tag_SNP

                        # Change the tag SNP flag to No for SNPs with the larger p-values
                        for pair in gene.SNP_pairs:
                            if pair.eQTL_SNP != tag_SNP:
                                if pair.eQTL_SNP in linked_SNPs:
                                    pair.is_tag_SNP = 'No'

        # Now deal with the regions, mostly HLA, that do not have good coverage in our LD database.  Use proximity and MAF to identify tag SNPs
        for chrom in chroms:

            # The exclude list ensures that we don't remove tag SNPs due to long-distance chaining of SNP blocks
            exclude_list = []

            # Get the list of SNPs not in the LD database that have corresponding GWAS
            SNPs_tmp = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.in_LD_DB == 'No' and foo.GWAS_SNP_pval != 'NA']
            SNPs = apply_eQTL_threshold(SNPs_tmp)

            for item in [str(foo.eQTL_SNP) for foo in SNPs]:
                if item not in exclude_list:
                    item_location = next(foo.eQTL_SNP_location for foo in gene.SNP_pairs if foo.eQTL_SNP == item)
                    upper_bound = int(item_location) + 5000
                    lower_bound = int(item_location) - 5000
                    proximal_SNPs_tmp = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.in_LD_DB == 'No' and foo.GWAS_SNP_pval != 'NA' and foo.eQTL_SNP_location > lower_bound and foo.eQTL_SNP_location < upper_bound]
                    proximal_SNPs = [str(foo.eQTL_SNP) for foo in apply_eQTL_threshold(proximal_SNPs_tmp)]

                    proximal_SNPs.sort()
                    for boo in proximal_SNPs:
                        exclude_list.append(boo)

                    # Get the eQTL SNP p-value and MAF for these results (will always included at least the SNP under test (aka item))
                    comparison_list = []
                    for proximal_SNP in proximal_SNPs:
                        foo = next(float(foo.eQTL_SNP_pval) for foo in gene.SNP_pairs if foo.eQTL_SNP == proximal_SNP)
                        comparison_list.append([proximal_SNP, foo])   # These are rs# and the eQTL p-value

                        # Want to select the strongest (smallest) eQTL p-value for all the SNPs that are linked.
                        comparison_list.sort(key=lambda x: x[1], reverse=False)
                        #print comparison_list
                        tag_SNP = comparison_list[0][0]

                    #print 'Tag SNP is: ' + str(tag_SNP)
                    # Change the tag SNP flag to No for SNPs with the larger p-values
                    for pair in gene.SNP_pairs:
                        if pair.eQTL_SNP == tag_SNP:
                                pair.is_tag_SNP = 'Yes'


 # Now verify that SNPs aren't in close proximity
        min_separation = 100000
        for chrom in chroms:

            untag_snp = 'start'
            while untag_snp != 'none':
                SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.is_tag_SNP == 'Yes']
                #print SNPs
                untag_snp = 'none'
                for i in range(1,len(SNPs)):
                    if int(abs(SNPs[i-1].eQTL_SNP_location - SNPs[i].eQTL_SNP_location)) < int(min_separation):
                        if float(SNPs[i-1].eQTL_SNP_pval) <= float(SNPs[i].eQTL_SNP_pval):
                            untag_snp = SNPs[i].eQTL_SNP
                        else:
                            untag_snp = SNPs[i-1].eQTL_SNP
                        # Now update
                        for pair in gene.SNP_pairs:
                            if pair.eQTL_SNP == untag_snp:
                                #print 'Untaging SNP: ' + str(untag_snp)
                                pair.is_tag_SNP = 'No'

                        # break out of the for loop and repeate with updated gene
                        break


 # Now verify that adjacent SNPs have rather different MAFs
        min_separation = 150000
        min_MAF_diff = 0.1
        for chrom in chroms:

            untag_snp = 'start'
            while untag_snp != 'none':
                SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.is_tag_SNP == 'Yes']
                #print SNPs
                untag_snp = 'none'
                for i in range(1,len(SNPs)):
                    try:
                        if int(abs(SNPs[i-1].eQTL_SNP_location - SNPs[i].eQTL_SNP_location)) < int(min_separation) and (abs(float(SNPs[i-1].eQTL_SNP_MAF) - float(SNPs[i].eQTL_SNP_MAF)) < float(min_MAF_diff)):
                            #print "Untag for gene: " + str(gene.name)
                            if float(SNPs[i-1].eQTL_SNP_pval) <= float(SNPs[i].eQTL_SNP_pval):
                                untag_snp = SNPs[i].eQTL_SNP
                            else:
                                untag_snp = SNPs[i-1].eQTL_SNP
                            # Now update
                            for pair in gene.SNP_pairs:
                                if pair.eQTL_SNP == untag_snp:
                                    #print 'Untaging SNP: ' + str(untag_snp)
                                    pair.is_tag_SNP = 'No'

                            # break out of the for loop and repeate with updated gene
                            break
                    except:
                        # In some cases, there is no MAF in the database, so we just skip this step
                        pass


 # Do not tag chromosome X
        SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == 'chrX' and foo.is_tag_SNP == 'Yes']
        #print SNPs
        for i in range(len(SNPs)):
            for pair in gene.SNP_pairs:
                if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
                    pair.is_tag_SNP = 'No'

# Do not tag the HLA region
        SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == 'chr6' and foo.is_tag_SNP == 'Yes' and (int(foo.eQTL_SNP_location) > int(27000000) and int(foo.eQTL_SNP_location) < int(33000000))]
        #print SNPs
        for i in range(len(SNPs)):
            for pair in gene.SNP_pairs:
                if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
                    pair.is_tag_SNP = 'No'

 ## Set a maximum eQTL p-value threshold for trans SNPs
 #       SNPs = [foo for foo in gene.SNP_pairs if (float(foo.eQTL_SNP_pval) >= float(3e-6)) and foo.is_tag_SNP == 'Yes' and foo.is_cis_SNP == 'trans']
 #       #print SNPs
 #       for i in range(len(SNPs)):
 #           for pair in gene.SNP_pairs:
 #               if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
 #                   pair.is_tag_SNP = 'No'

 # Set a maximum eQTL p-value threshold for cis SNPs
        SNPs = [foo for foo in gene.SNP_pairs if (float(foo.eQTL_SNP_pval) >= float(1e-3)) and foo.is_tag_SNP == 'Yes' and foo.is_cis_SNP == 'cis']
        #print SNPs
        for i in range(len(SNPs)):
            for pair in gene.SNP_pairs:
                if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
                    pair.is_tag_SNP = 'No'


    return gene_list

# Determine which SNPs from matched SNP pairs should be used as tag SNPs.  These are the only SNPs that factor into our analysis
def assign_eQTL_tag_SNPs_preload_debug(gene_list, eQTLSNP_LD, LD_threshold,
        data_folder_location, eSNP_pos, eSNP_ld, eSNP_linked):
    """This function is for debug purpose to make sure using preload data
    produce identical data compared to tabix
    """
    #pdb.set_trace()
    # Extract eQTL SNPs for each gene and create blocks
    print ('assign_eQTL_tag_SNPs called')

    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
            'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
            'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
            'chr22', 'chrX', 'chrY']
    eQTL_dict={}
    #prepare tabix for LD 0.2, used merged LD02.tab.gz
    #url=data_folder_location+'LD/LD02.tab.gz'
    #testing use /run/shm/LD/LD02.tab.gz'
    url = '/run/shm/LD/LD02.tab.gz'
    tabix_LD_02 = tabix.open(url)

    for gene in gene_list:

        # Identify gene
        # print gene.name
        time_gene = time.time()
        snp_index = {}
        eQTL_pvals = {}
        for pair_index, snp_pair in enumerate(gene.SNP_pairs):
            snp_index[snp_pair.eQTL_SNP] = pair_index
            eQTL_pvals[snp_pair.eQTL_SNP] = float(snp_pair.eQTL_SNP_pval)

        # Build blocks on a chromosome-by-chromosome basis (obviously!!)
        for chrom in chroms:

            # First make blocks using SNPs that are in the LD_DB and are on the same chromosome.
            #print gene.SNP_pairs
            SNPs_tmp = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.in_LD_DB == 'Yes' and foo.GWAS_SNP != 'NF']
            SNPs = apply_eQTL_threshold(SNPs_tmp)

            # Store linked SNPs for each chrom as a list of lists
            # representing the block structure
            chrom_blocks = []

            for item in [str(foo.eQTL_SNP) for foo in SNPs]:
                # take advantage of the fact that eQTLSNP_LD has been
                # constructed from match_eQTL_SNPs_with_GWAS
                qr=[]
                if item in eQTLSNP_LD:
                    if not (item in eQTL_dict):
                        eQTL_dict[item]={}
                        for qritem in eQTLSNP_LD[item]:
                            eQTL_dict[item][qritem[0]]=qritem[1]
                    for eqtlsnp in [str(foo.eQTL_SNP) for foo in SNPs]:
                        if eqtlsnp in eQTL_dict[item]:
                            #print "found "+qritem[0]
                            qr.append([eqtlsnp,eQTL_dict[item][eqtlsnp]])
                        #else:
                        #    print qritem[0]+" is not in SNPDict"
                # Identify the SNPs in high LD that have the strongest eQTL
                #print "LD for "+item
                #print qr
                linked_SNPs = [[item, 1.00]]
                if qr != []:
                    for result in qr:
                        linked_SNPs.append([result[0], result[1]])

                # Order the SNPs in each block to ensure equivalent
                # resutls when the eQTL p-values happen to be identical
                linked_SNPs.sort(key=lambda x: x[0])

                # Add anything we've found to the block list for this chrom
                chrom_blocks.append(linked_SNPs)


            # Process the blocks for each chrom
            for linked_SNPs in chrom_blocks:

                # Get the eQTL SNP p-value and MAF for these results
                comparison_list = []
                for linked_SNP in linked_SNPs:
                    #print "linked_SNP:"+str(linked_SNP)
                    #for snp in gene.SNP_pairs:
                    #    print snp.eQTL_SNP+"\t"+snp.eQTL_SNP_pval
                    foo = next(float(foo.eQTL_SNP_pval) for foo in gene.SNP_pairs if foo.eQTL_SNP == linked_SNP[0])
                    #$print "foo:"+str(foo)
                    # These are rs#, rsq, and the eQTL p-value
                    comparison_list.append(
                            [linked_SNP[0], linked_SNP[1], foo])
                    # These are rs#, rsq, and the eQTL p-value

                # Want to select the strongest (smallest) eQTL p-value for
                # all the SNPs that are linked.
                comparison_list.sort(key=lambda x: x[2], reverse=False)
                #print comparison_list
                tag_SNP = comparison_list[0][0]

                # Change the tag SNP flag
                for pair in gene.SNP_pairs:
                    if pair.eQTL_SNP == tag_SNP:
                        pair.is_tag_SNP = 'Yes'

                # To debug this, print everything
                #print
                #print comparison_list
                #print 'tag_SNP: ' + tag_SNP


        # Now verify that the tagged genes are not in LD with each other
        for chrom in chroms:

            # The exclude list ensures that we don't remove tag SNPs due to long-distance chaining of SNP blocks
            exclude_list = {}
            exclude_list2 = {}

            #SNPs = [foo.eQTL_SNP for foo in gene.SNP_pairs if foo.chrom == chrom and foo.is_tag_SNP == 'Yes']
            # instead of creating a list, use a dictionary
            # keep snps in the original order
            sortedSNPs = []
            # SNPValues = {'snpid':[location, pval]}
            SNPValues = {}
            for item in gene.SNP_pairs:
                if item.is_tag_SNP == 'Yes' and item.chrom == chrom:
                    SNPValues[item.eQTL_SNP] = [item.eQTL_SNP_location,
                            item.eQTL_SNP_pval]
                    sortedSNPs.append(item.eQTL_SNP)

            # Store linked SNPs for each chrom as a list of lists representing the block structure
            for SNP in sortedSNPs:
                #print "checking",SNP
                #print "Excluded List",exclude_list
                if SNP not in exclude_list2:
                    # SNPValues[SNP] is [pos, pval]
                    snp_location = str(SNPValues[SNP][0])
                    position_str = chrom + ":" + snp_location + \
                                           "-" + snp_location
                    #records = []
                    #try:
                    #    #query(chrom,pos,pos)
                    #    records = tabix_LD_02.querys(position_str)
                    #except:
                    #    print "error query with tabix: ", position_str

                    #qr = []
                    #for record in records:
                    #    # [chrom,pos1,snp1,pos2,snp2,ld]
                    #    if len(record)>=6 and record[4] in SNPValues:
                    #        qr.append([record[4], record[5]])
                    qr2 = []
                    if SNP in eSNP_linked:
                        for snp2 in eSNP_linked[SNP]:
                            if snp2 in SNPValues:
                                qr2.append([snp2, eSNP_ld[SNP][snp2]])
                                #linked_SNPs.append(snp2)
                                #exclude_list.setdefault(snp2,1)


                    # Identify the SNPs in high LD that have the strongest eQTL
                    #linked_SNPs = [SNP]
                    #if qr:
                    #    # Some of the SNPs are still in LD,
                    #    # we will need to remove the extras.
                    #    for result in qr:
                    #        linked_SNPs.append(result[0])
                    #    #linked_SNPs.sort()
                    #    for snp_B in linked_SNPs:
                    #        exclude_list.setdefault(snp_B,1)
                    linked_SNPs2 = [SNP]
                    if qr2:
                        for result in qr2:
                            linked_SNPs2.append(result[0])
                        for snp_B in linked_SNPs2:
                            exclude_list2.setdefault(snp_B,1)

                    #linked_SNPs and linked_SNP2 should be exactly the same
                    #if not (len(linked_SNPs) == len(linked_SNPs2)):
                    #    print 'Diff',gene.name, SNP, 'linked SNP not equal', position_str
                    #    print str(len(linked_SNPs)),str(len(linked_SNPs2))
                    #    print 'tabix:',linked_SNPs
                    #    print 'preload:',linked_SNPs2
                    #else:
                    #    for si, lsnp in enumerate(linked_SNPs):
                    #        lsnp2 = linked_SNPs2[si]
                    #        if not lsnp==lsnp2:
                    #            print 'DIFF', position_str
                    #            print gene.name, SNP, si,'in tabix:',lsnp
                    #            print gene.name, SNP, si,'in preload:',lsnp2


                    # If we find other tag SNPs in LD, we need to check their eQTL p-values and remove the large ones
                    if len(linked_SNPs2) > 0:
                        #print 'linked snp', linked_SNPs
                        # Get the eQTL SNP p-value and MAF for these results
                        comparison_list = [] #[SNP, float(SNPValues[SNP][1])]]
                        for linked_SNP in linked_SNPs2:
                            #snp_pval = next(float(foo.eQTL_SNP_pval) for foo in gene.SNP_pairs if foo.eQTL_SNP == linked_SNP)
                            snp_pval = float(SNPValues[linked_SNP][1])
                            comparison_list.append([linked_SNP, snp_pval])   # These are rs# and the eQTL p-value

                        # Want to select the strongest (smallest) eQTL p-value for all the SNPs that are linked.
                        comparison_list.sort(key=lambda x: x[1], reverse=False)
                        #print 'comparison_list', comparison_list
                        tag_SNP = comparison_list[0][0]
                        #print 'set tag_SNP', tag_SNP

                        # Change the tag SNP flag to No for SNPs with the larger p-values
                        #for pair in gene.SNP_pairs:
                        #    if pair.eQTL_SNP != tag_SNP:
                        #        if pair.eQTL_SNP in linked_SNPs2:
                        #            pair.is_tag_SNP = 'No'

                        # save the tag_SNP_state for the tag_SNP
                        # because the tag_SNP may already be in 'No' state
                        # in this case we are not going to tag it back to 'Yes'
                        tag_index = snp_index[tag_SNP]
                        tag_SNP_state = gene.SNP_pairs[tag_index].is_tag_SNP
                        # untag all in linkage first
                        for snp in linked_SNPs2:
                            pair_index = snp_index[snp]
                            gene.SNP_pairs[pair_index].is_tag_SNP = 'No'
                        # revert the tag_SNP back to its old state
                        gene.SNP_pairs[tag_index].is_tag_SNP = tag_SNP_state



        # Now deal with the regions, mostly HLA, that do not have good coverage in our LD database.  Use proximity and MAF to identify tag SNPs
        for chrom in chroms:

            # The exclude list ensures that we don't remove tag SNPs due to long-distance chaining of SNP blocks
            exclude_list = []

            # Get the list of SNPs not in the LD database that have corresponding GWAS
            SNPs_tmp = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.in_LD_DB == 'No' and foo.GWAS_SNP_pval != 'NA']
            SNPs = apply_eQTL_threshold(SNPs_tmp)

            for item in [str(foo.eQTL_SNP) for foo in SNPs]:
                if item not in exclude_list:
                    item_location = next(foo.eQTL_SNP_location for foo in gene.SNP_pairs if foo.eQTL_SNP == item)
                    upper_bound = int(item_location) + 5000
                    lower_bound = int(item_location) - 5000
                    proximal_SNPs_tmp = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.in_LD_DB == 'No' and foo.GWAS_SNP_pval != 'NA' and foo.eQTL_SNP_location > lower_bound and foo.eQTL_SNP_location < upper_bound]
                    proximal_SNPs = [str(foo.eQTL_SNP) for foo in apply_eQTL_threshold(proximal_SNPs_tmp)]

                    proximal_SNPs.sort()
                    for boo in proximal_SNPs:
                        exclude_list.append(boo)

                    # Get the eQTL SNP p-value and MAF for these results (will always included at least the SNP under test (aka item))
                    comparison_list = []
                    for proximal_SNP in proximal_SNPs:
                        foo = next(float(foo.eQTL_SNP_pval) for foo in gene.SNP_pairs if foo.eQTL_SNP == proximal_SNP)
                        comparison_list.append([proximal_SNP, foo])   # These are rs# and the eQTL p-value

                        # Want to select the strongest (smallest) eQTL p-value for all the SNPs that are linked.
                        comparison_list.sort(key=lambda x: x[1], reverse=False)
                        #print comparison_list
                        tag_SNP = comparison_list[0][0]

                    #print 'Tag SNP is: ' + str(tag_SNP)
                    # Change the tag SNP flag to No for SNPs with the larger p-values
                    for pair in gene.SNP_pairs:
                        if pair.eQTL_SNP == tag_SNP:
                                pair.is_tag_SNP = 'Yes'


 # Now verify that SNPs aren't in close proximity
        min_separation = 100000
        for chrom in chroms:

            untag_snp = 'start'
            while untag_snp != 'none':
                SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.is_tag_SNP == 'Yes']
                #print SNPs
                untag_snp = 'none'
                for i in range(1,len(SNPs)):
                    if int(abs(SNPs[i-1].eQTL_SNP_location - SNPs[i].eQTL_SNP_location)) < int(min_separation):
                        if float(SNPs[i-1].eQTL_SNP_pval) <= float(SNPs[i].eQTL_SNP_pval):
                            untag_snp = SNPs[i].eQTL_SNP
                        else:
                            untag_snp = SNPs[i-1].eQTL_SNP
                        # Now update
                        for pair in gene.SNP_pairs:
                            if pair.eQTL_SNP == untag_snp:
                                #print 'Untaging SNP: ' + str(untag_snp)
                                pair.is_tag_SNP = 'No'

                        # break out of the for loop and repeate with updated gene
                        break


 # Now verify that adjacent SNPs have rather different MAFs
        min_separation = 150000
        min_MAF_diff = 0.1
        for chrom in chroms:

            untag_snp = 'start'
            while untag_snp != 'none':
                SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.is_tag_SNP == 'Yes']
                #print SNPs
                untag_snp = 'none'
                for i in range(1,len(SNPs)):
                    try:
                        if int(abs(SNPs[i-1].eQTL_SNP_location - SNPs[i].eQTL_SNP_location)) < int(min_separation) and (abs(float(SNPs[i-1].eQTL_SNP_MAF) - float(SNPs[i].eQTL_SNP_MAF)) < float(min_MAF_diff)):
                            #print "Untag for gene: " + str(gene.name)
                            if float(SNPs[i-1].eQTL_SNP_pval) <= float(SNPs[i].eQTL_SNP_pval):
                                untag_snp = SNPs[i].eQTL_SNP
                            else:
                                untag_snp = SNPs[i-1].eQTL_SNP
                            # Now update
                            for pair in gene.SNP_pairs:
                                if pair.eQTL_SNP == untag_snp:
                                    #print 'Untaging SNP: ' + str(untag_snp)
                                    pair.is_tag_SNP = 'No'

                            # break out of the for loop and repeate with updated gene
                            break
                    except:
                        # In some cases, there is no MAF in the database, so we just skip this step
                        pass


 # Do not tag chromosome X
        SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == 'chrX' and foo.is_tag_SNP == 'Yes']
        #print SNPs
        for i in range(len(SNPs)):
            for pair in gene.SNP_pairs:
                if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
                    pair.is_tag_SNP = 'No'

# Do not tag the HLA region
        SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == 'chr6' and foo.is_tag_SNP == 'Yes' and (int(foo.eQTL_SNP_location) > int(27000000) and int(foo.eQTL_SNP_location) < int(33000000))]
        #print SNPs
        for i in range(len(SNPs)):
            for pair in gene.SNP_pairs:
                if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
                    pair.is_tag_SNP = 'No'

 ## Set a maximum eQTL p-value threshold for trans SNPs
 #       SNPs = [foo for foo in gene.SNP_pairs if (float(foo.eQTL_SNP_pval) >= float(3e-6)) and foo.is_tag_SNP == 'Yes' and foo.is_cis_SNP == 'trans']
 #       #print SNPs
 #       for i in range(len(SNPs)):
 #           for pair in gene.SNP_pairs:
 #               if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
 #                   pair.is_tag_SNP = 'No'

 # Set a maximum eQTL p-value threshold for cis SNPs
        SNPs = [foo for foo in gene.SNP_pairs if (float(foo.eQTL_SNP_pval) >= float(1e-3)) and foo.is_tag_SNP == 'Yes' and foo.is_cis_SNP == 'cis']
        #print SNPs
        for i in range(len(SNPs)):
            for pair in gene.SNP_pairs:
                if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
                    pair.is_tag_SNP = 'No'

        # time used by each gene
        print ('Tagged',gene.name,'used',str(time.time()-time_gene),'secs')

    return gene_list


# Determine which SNPs from matched SNP pairs should be used as tag SNPs.  These are the only SNPs that factor into our analysis
def assign_eQTL_tag_SNPs_preload(gene_list, eQTLSNP_LD, LD_threshold,
        eSNP_pos, eSNP_ld, eSNP_linked, ld_range, eQTL_threshold = 1e-5):
    """Using a preload eSNPs dictionary with LD>0.2 (included only eSNPs)

    Modified 08/18/2015: some snpID were not in LD database, or were
    represented in LD database using another rsID(sharing the same location)
    So we modified the preprocessing procedure to include the rsID in eQTL,

    Tabix approach using the chromosome position to query the SNP LD, thus
    will be able to pull down the linkage information for that chromosome pos,
    even when the snpID was missing in the LD tabix file.

    For preload to work the same way, the preload eSNP_linked table use rsID
    as keys, so we need to store the rsID used in the eQTL also as the first
    column when generating these files (see the data_preparation/eQTL/)

    Modified: 8/12/17: some rsIDs were not in the LD database or were
    represented by other rsIDs, so when filtering, need to take into acount
    of that. Fixed this problem by fix the data, in the LD data, find the snp
    that are alias to the SNPs in eQTL. And add corresponding lines to the LD
    data

    """
    # Extract eQTL SNPs for each gene and create blocks
    print ('assign_eQTL_tag_SNPs called')

    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
            'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
            'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
            'chr22', 'chrX', 'chrY']
    eQTL_dict={}

    for gene in gene_list:
        # record start time
        time_gene = time.time()
        # prebuild a dictionary from SNP to index in SNP_pairs
        snp_index = {}
        eQTL_pvals = {}
        snp_positions = {}
        for pair_index, snp_pair in enumerate(gene.SNP_pairs):
            snp_index[snp_pair.eQTL_SNP] = pair_index
            eQTL_pvals[snp_pair.eQTL_SNP] = float(snp_pair.eQTL_SNP_pval)
            snp_positions[snp_pair.eQTL_SNP] = snp_pair.eQTL_SNP_location

        # Build blocks on a chromosome-by-chromosome basis (obviously!!)
        for chrom in chroms:

            # First make blocks using SNPs that are in the LD_DB and
            # are on the same chromosome.
            SNPs_tmp = [foo for foo in gene.SNP_pairs if foo.chrom == chrom
                    and foo.in_LD_DB == 'Yes' and foo.GWAS_SNP != 'NF']
            SNPs = apply_eQTL_threshold(SNPs_tmp,eQTL_threshold)

            # Store linked SNPs for each chrom as a list of lists
            # representing the block structure
            chrom_blocks = []

            for item in [str(foo.eQTL_SNP) for foo in SNPs]:
                # take advantage of the fact that eQTLSNP_LD has been
                # constructed from match_eQTL_SNPs_with_GWAS
                qr=[]
                if item in eQTLSNP_LD:
                    if item not in eQTL_dict:
                        eQTL_dict[item]={}
                        for qritem in eQTLSNP_LD[item]:
                            eQTL_dict[item][qritem[0]]=qritem[1]
                    for eqtlsnp in [str(foo.eQTL_SNP) for foo in SNPs]:
                        if eqtlsnp in eQTL_dict[item]:
                            #print "found "+qritem[0]
                            qr.append([eqtlsnp,eQTL_dict[item][eqtlsnp]])
                        #else:
                        #    print qritem[0]+" is not in SNPDict"
                # Identify the SNPs in high LD that have the strongest eQTL
                #print "LD for "+item
                #print qr
                linked_SNPs = [[item, 1.00]]
                if qr != []:
                    for result in qr:
                        linked_SNPs.append([result[0], result[1]])

                # Order the SNPs in each block to ensure equivalent
                # resutls when the eQTL p-values happen to be identical
                linked_SNPs.sort(key=lambda x: x[0])

                # Add anything we've found to the block list for this chrom
                chrom_blocks.append(linked_SNPs)


            # Process the blocks for each chrom
            for linked_SNPs in chrom_blocks:

                # Get the eQTL SNP p-value and MAF for these results
                comparison_list = []
                for linked_SNP in linked_SNPs:
                    eQTL_p = eQTL_pvals[linked_SNP[0]]
                    comparison_list.append([linked_SNP[0], linked_SNP[1],
                        eQTL_p])   # These are rs#, rsq, and the eQTL p-value

                # Want to select the strongest (smallest) eQTL p-value for
                # all the SNPs that are linked.
                comparison_list.sort(key=lambda x: x[2], reverse=False)
                tag_SNP = comparison_list[0][0]

                # Change the tag SNP flag (use snp_index for tag_index)
                tag_index = snp_index[tag_SNP]
                gene.SNP_pairs[tag_index].is_tag_SNP = 'Yes'


        # Now verify that the tagged SNP are not in LD with each other
        for chrom in chroms:

            # excluded_SNPs ensures that we don't remove tag SNPs
            # due to long-distance chaining of SNP blocks,
            excluded_SNPs = {}

            # keep snps in the original order (sorted by position)
            sortedSNPs = []
            # SNPValues = {'snpid':[location, pval]}
            SNPValues = {}
            for item in gene.SNP_pairs:
                # identify the snps that are tagged and in given chromosome
                if item.is_tag_SNP == 'Yes' and item.chrom == chrom:
                    SNPValues[item.eQTL_SNP] = [item.eQTL_SNP_location,
                            item.eQTL_SNP_pval]
                    sortedSNPs.append(item.eQTL_SNP)

            for SNP in sortedSNPs:
                #print "checking",SNP
                #print "Excluded List",exclude_list
                if SNP not in excluded_SNPs:
                    # to preserve first SNP, initialize with linked_SNPs=[]
                    # not preserve first SNP, initialize with linked_SNPs=[SNP]
                    linked_SNPs = [SNP]
                    if SNP in eSNP_linked:
                        # use the eSNP_linked to go through the snp2 (sorted)
                        # (the order of snp2 in eSNP_lds[SNP] will be shuffled
                        # and sometimes give inconsistent results especially
                        # when two snp share the same p-values)
                        for snp2 in eSNP_linked[SNP]:
                            if snp2 in SNPValues:
                                #qr.append([snp2, eSNP_ld[SNP][snp2]])
                                linked_SNPs.append(snp2)
                                excluded_SNPs.setdefault(snp2,1)
                    #else: #means there is no snpB in eQTL that is in LD02
                    #    print "error no LD found:", SNP, position_str

                    # Identify the SNPs in high LD that have the strongest eQTL
                    # Some of the SNPs are still in LD,
                    # we will need to remove the extras.
                    # If we find other tag SNPs in LD, we need to check their
                    # eQTL p-values and remove the large ones
                    if len(linked_SNPs) > 0:
                        #print 'linked snp', linked_SNPs
                        # Get the eQTL SNP p-value and MAF for these results
                        comparison_list = [] #[SNP, float(SNPValues[SNP][1])]]
                        for linked_SNP in linked_SNPs:
                            snp_pval = eQTL_pvals[linked_SNP]
                            # comparizon_list = [[rsid,eqtl_pval]...]
                            comparison_list.append([linked_SNP, snp_pval])

                        # Want to select the strongest (smallest) eQTL p-value
                        # for all the SNPs that are linked.
                        comparison_list.sort(key=lambda x: x[1], reverse=False)
                        tag_SNP = comparison_list[0][0]

                        # Change the tag SNP flag to No for SNPs with
                        # larger p-values. Do not touch the state (Yes or No)
                        # of tag_SNP.
                        # First save the tag_SNP_state for the tag_SNP
                        # the tag_SNP may already be in 'No' state
                        # and we are not going to tag it back to 'Yes'
                        tag_index = snp_index[tag_SNP]
                        tag_SNP_state = gene.SNP_pairs[tag_index].is_tag_SNP
                        # untag all in linkage first
                        for snp in linked_SNPs:
                            pair_index = snp_index[snp]
                            gene.SNP_pairs[pair_index].is_tag_SNP = 'No'
                        # revert the tag_SNP back to its old state
                        gene.SNP_pairs[tag_index].is_tag_SNP = tag_SNP_state


        # Now deal with the regions, mostly HLA, that do not have good
        # coverage in our LD database. Use proximity and MAF to identify
        # tag SNPs
        for chrom in chroms:
            # The exclude dict ensures that we don't remove tag SNPs due to
            # long-distance chaining of SNP blocks
            exclude_list = {}

            # Get the list of SNPs not in the LD database that have
            # corresponding GWAS
            SNPs_tmp = [foo for foo in gene.SNP_pairs if foo.chrom == chrom \
                    and foo.in_LD_DB == 'No' and foo.GWAS_SNP_pval != 'NA']
            SNPs = apply_eQTL_threshold(SNPs_tmp)

            for item in [str(foo.eQTL_SNP) for foo in SNPs]:
                if item not in exclude_list:
                    #item_location = next(foo.eQTL_SNP_location for foo in gene.SNP_pairs if foo.eQTL_SNP == item)
                    snp_pos = snp_positions[item]
                    upper_bound = int(snp_pos) + 5000
                    lower_bound = int(snp_pos) - 5000
                    proximal_SNPs_tmp = [foo for foo in gene.SNP_pairs if foo.chrom == chrom and foo.in_LD_DB == 'No' and foo.GWAS_SNP_pval != 'NA' and foo.eQTL_SNP_location > lower_bound and foo.eQTL_SNP_location < upper_bound]
                    proximal_SNPs = [str(foo.eQTL_SNP) for foo in apply_eQTL_threshold(proximal_SNPs_tmp)]

                    proximal_SNPs.sort()
                    for proximal_snp in proximal_SNPs:
                        exclude_list.setdefault(proximal_snp,0)

                    # Get the eQTL SNP p-value and MAF for these results
                    # (will always included at least the SNP under
                    # test (aka item))
                    comparison_list = []
                    for proximal_SNP in proximal_SNPs:
                        snp_pval = eQTL_pvals[proximal_SNP]
                        #foo = next(float(foo.eQTL_SNP_pval) for foo in gene.SNP_pairs if foo.eQTL_SNP == proximal_SNP)
                        comparison_list.append([proximal_SNP, snp_pval])
                        # These are rs# and the eQTL p-value

                    # Want to select the strongest (smallest) eQTL p-value
                    # for all the SNPs that are linked.
                    comparison_list.sort(key=lambda x: x[1], reverse=False)
                    #print comparison_list
                    tag_SNP = comparison_list[0][0]

                    # use the snp_index to set the state of the tagged SNP
                    tag_index = snp_index[tag_SNP]
                    gene.SNP_pairs[tag_index].is_tag_SNP = 'Yes'


        # Now verify that SNPs aren't in close proximity
        min_separation = 100000
        for chrom in chroms:

            untag_snp = 'start'
            while untag_snp != 'none':
                SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == chrom \
                        and foo.is_tag_SNP == 'Yes']
                #print SNPs
                untag_snp = 'none'
                for i in range(1,len(SNPs)):
                    if int(abs(SNPs[i-1].eQTL_SNP_location -
                            SNPs[i].eQTL_SNP_location)) < int(min_separation):
                        if float(SNPs[i-1].eQTL_SNP_pval) <= \
                                float(SNPs[i].eQTL_SNP_pval):
                            untag_snp = SNPs[i].eQTL_SNP
                        else:
                            untag_snp = SNPs[i-1].eQTL_SNP
                        # Now update
                        if untag_snp in snp_index:
                            untag_index = snp_index[untag_snp]
                            gene.SNP_pairs[untag_index].is_tag_SNP = 'No'

                        # break out of the for loop and repeate
                        # with updated gene
                        break


        # Now verify that adjacent SNPs have rather different MAFs
        min_separation = 150000
        min_MAF_diff = 0.1
        for chrom in chroms:
            untag_snp = 'start'
            while untag_snp != 'none':
                SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == chrom \
                        and foo.is_tag_SNP == 'Yes']
                #print SNPs
                untag_snp = 'none'
                for i in range(1,len(SNPs)):
                    try:
                        if int(abs(SNPs[i-1].eQTL_SNP_location - SNPs[i].eQTL_SNP_location)) < int(min_separation) and (abs(float(SNPs[i-1].eQTL_SNP_MAF) - float(SNPs[i].eQTL_SNP_MAF)) < float(min_MAF_diff)):
                            #print "Untag for gene: " + str(gene.name)
                            if float(SNPs[i-1].eQTL_SNP_pval) <= float(SNPs[i].eQTL_SNP_pval):
                                untag_snp = SNPs[i].eQTL_SNP
                            else:
                                untag_snp = SNPs[i-1].eQTL_SNP
                            # Now update
                            #for pair in gene.SNP_pairs:
                            #    if pair.eQTL_SNP == untag_snp:
                                    #print 'Untaging SNP: ' + str(untag_snp)
                            #        pair.is_tag_SNP = 'No'

                            if untag_snp in snp_index:
                                untag_index = snp_index[untag_snp]
                                gene.SNP_pairs[untag_index].is_tag_SNP = 'No'


                            # break out of the for loop and repeate with updated gene
                            break
                    except:
                        # In some cases, there is no MAF in the database, so we just skip this step
                        pass

        #now verify that the ones not in LD database is far away enough from
        #the other tag snps, 9/3/17
        for chrom in chroms:
            untag_snp = 'start'
            while untag_snp != 'none':
                #list of snps with in_LD_DB=='No'
                SNPs_no_LD = [foo for foo in gene.SNP_pairs \
                        if foo.chrom == chrom \
                        and foo.is_tag_SNP == 'Yes' and foo.in_LD_DB == 'No']

                SNPs_LD = [foo for foo in gene.SNP_pairs \
                        if foo.chrom == chrom \
                        and foo.is_tag_SNP == 'Yes' and foo.in_LD_DB == 'Yes']
                #
                #print SNPs
                untag_snp = 'none'
                for i in range(len(SNPs_no_LD)-1):
                    snp_a = SNPs_no_LD[i].eQTL_SNP
                    snp_a_pos = SNPs_no_LD[i].eQTL_SNP_location

                    for j in range(i+1,len(SNPs_no_LD)):
                        if SNPs_no_LD[i].is_tag_SNP == 'Yes' and \
                                SNPs_no_LD[j].is_tag_SNP == 'Yes':
                            snp_b = SNPs_no_LD[j].eQTL_SNP
                            snp_b_pos = SNPs_no_LD[j].eQTL_SNP_location
                            linked = is_snps_in_LD_range(ld_range, snp_a,\
                                    snp_b, snp_a_pos, snp_b_pos)
                            if linked:
                                if float(SNPs_no_LD[i].eQTL_SNP_pval) <= float(SNPs_no_LD[j].eQTL_SNP_pval):
                                    untag_snp = SNPs_no_LD[j].eQTL_SNP
                                else:
                                    untag_snp = SNPs_no_LD[i].eQTL_SNP
                            # Now update
                            #for pair in gene.SNP_pairs:
                            #    if pair.eQTL_SNP == untag_snp:
                                    #print 'Untaging SNP: ' + str(untag_snp)
                            #        pair.is_tag_SNP = 'No'

                                if untag_snp in snp_index:
                                    untag_index = snp_index[untag_snp]
                                    gene.SNP_pairs[untag_index].is_tag_SNP = 'No'
                                    print ("Untag snp not in LD: " +\
                                        str(gene.name))
                                    print (gene.SNP_pairs[untag_index].eQTL_SNP)

                #then go through the SNPs_LD
                untag_snp = 'none'
                for i in range(len(SNPs_no_LD)):
                    snp_a = SNPs_no_LD[i].eQTL_SNP
                    snp_a_pos = SNPs_no_LD[i].eQTL_SNP_location
                    for j in range(len(SNPs_LD)):
                        if SNPs_no_LD[i].is_tag_SNP == 'Yes' and \
                                SNPs_LD[j].is_tag_SNP == 'Yes':
                            snp_b = SNPs_LD[j].eQTL_SNP
                            snp_b_pos = SNPs_LD[j].eQTL_SNP_location
                            #use a slightly different function
                            #to use only the range in snp_b to test
                            #snp_a is snp1, snp_b is snp2(with LD information)
                            linked = is_snp1_in_LD_snp2(ld_range, snp_a,\
                                    snp_b, snp_a_pos, snp_b_pos)
                            if linked:
                                print (snp_a, snp_a_pos, SNPs_no_LD[i].eQTL_SNP_pval)
                                print (snp_b, snp_b_pos, SNPs_LD[j].eQTL_SNP_pval)

                                if float(SNPs_no_LD[i].eQTL_SNP_pval) <= \
                                        float(SNPs_LD[j].eQTL_SNP_pval):
                                    untag_snp = SNPs_LD[j].eQTL_SNP
                                else:
                                    untag_snp = SNPs_no_LD[i].eQTL_SNP

                                if untag_snp in snp_index:
                                    untag_index = snp_index[untag_snp]
                                    gene.SNP_pairs[untag_index].is_tag_SNP = 'No'
                                    print ("Untag snp not in LD compare with snp in LD: " + str(gene.name))
                                    print (gene.SNP_pairs[untag_index].eQTL_SNP)


        # Do not tag chromosome X
        SNPs_on_X = [foo for foo in gene.SNP_pairs if foo.chrom == 'chrX' \
                and foo.is_tag_SNP == 'Yes']
        for pair_on_X in SNPs_on_X:
            if pair_on_X.eQTL_SNP in snp_index:
                untag_index = snp_index[pair_on_X.eQTL_SNP]
                gene.SNP_pairs[untag_index].is_tag_SNP = 'No'

        # Do not tag the HLA region
        SNPs = [foo for foo in gene.SNP_pairs if foo.chrom == 'chr6' and
                foo.is_tag_SNP == 'Yes' and (int(foo.eQTL_SNP_location) >
                    int(27000000) and int(foo.eQTL_SNP_location) <
                    int(33000000))]
        for pair_on_HLA in SNPs:
            if pair_on_HLA.eQTL_SNP in snp_index:
                untag_index = snp_index[pair_on_HLA.eQTL_SNP]
                gene.SNP_pairs[untag_index].is_tag_SNP = 'No'



 ## Set a maximum eQTL p-value threshold for trans SNPs
 #       SNPs = [foo for foo in gene.SNP_pairs if (float(foo.eQTL_SNP_pval) >= float(3e-6)) and foo.is_tag_SNP == 'Yes' and foo.is_cis_SNP == 'trans']
 #       #print SNPs
 #       for i in range(len(SNPs)):
 #           for pair in gene.SNP_pairs:
 #               if pair.eQTL_SNP == SNPs[i].eQTL_SNP:
 #                   pair.is_tag_SNP = 'No'

 # Set a maximum eQTL p-value threshold for cis SNPs
        SNPs = [foo for foo in gene.SNP_pairs if (float(foo.eQTL_SNP_pval) >= float(1e-3)) and foo.is_tag_SNP == 'Yes' and foo.is_cis_SNP == 'cis']
        for bad_cis_pair in SNPs:
            if bad_cis_pair.eQTL_SNP in snp_index:
                untag_index = snp_index[bad_cis_pair.eQTL_SNP]
                gene.SNP_pairs[untag_index].is_tag_SNP = 'No'


        # time used by each gene
        print ('Tagged',gene.name,'used',str(time.time()-time_gene),'secs')

    return gene_list


####################################################################################################################
###   Functions to construct GWAS PDFs and Score PDFs
###   These functions define the core of our statistical test
####################################################################################################################


### Determine the number of times that a given SNP appears in the full list of genes.  We then use this to condition the PDF for each SNP added to the score.
def assign_pleiotropic_gene_count(gene_list):
    """This method is SNP based and assigns the gene.PGC field a number
    based simply on the number of times a SNP appears in the aligned data.
    Create a dictionary of the number of times that a GWAS SNP appears as
    a tag SNP.
    """
    aligned_SNP_dict = {}
    for gene in gene_list:
        #not counting hsa00xxx for pseudo genes(kegg pathways)
        if not gene.name.startswith('hsa'):
            for pair in gene.SNP_pairs:
                if pair.is_tag_SNP == 'Yes':
                    aligned_SNP_dict.setdefault(pair.GWAS_SNP,0)
                    aligned_SNP_dict[pair.GWAS_SNP] += 1

    # Now, use the dictionary to assign PGC to each tag SNP found
    for gene in gene_list:
        for pair in gene.SNP_pairs:
            if pair.is_tag_SNP == 'Yes':
                if pair.GWAS_SNP in aligned_SNP_dict:
                    pair.PGC = aligned_SNP_dict[pair.GWAS_SNP]
                else: #when only pseudo genes are considered
                    pair.PGC = 1

    return gene_list

def pleio_bin(PGC):
    if PGC <=5:
        return PGC-1
    elif PGC <=10:
        return 5
    elif PGC <=30:
        return 6
    else:
        return 7

def get_GWAS_pdf(gene_list, GWAS, database_cursor, bin_width, include_HLA_region):
    """Obtain GWAS p-value pdf (probability distribution function)
    """
    # Here the GWAS PDF is computed contingent on SNP PGC and loaded
    # into a matrix.
    # Thus, as SNPs are added to a score, we can convolve a PDF specific
    # to the pleiotropy of the added SNP.

    # Create a list of all unique SNPs
    aligned_SNP_dict = {}
    for gene in gene_list:
        for pair in gene.SNP_pairs:
            if pair.is_tag_SNP == 'Yes':
                if int(pair.PGC) > int(800):
                    PGC = 800
                else:
                    PGC = pair.PGC
                aligned_SNP_dict[pair.GWAS_SNP] = [float(PGC), float(pair.GWAS_SNP_pval)]

    #compute the mean and var of the pval
    all_p_values = []
    for SNP in aligned_SNP_dict:
        all_p_values.append(float(aligned_SNP_dict[SNP][1]))


    #pdb.set_trace()

    all_p_values = -numpy.log10(all_p_values)
    all_p_values = all_p_values[~numpy.isinf(all_p_values)]    # exclude log_pval==inf
    mean_p = numpy.mean(all_p_values)
    var_p = numpy.var(all_p_values)

    pdf_matrix_raw = []# raw means the raw SNP pvalues, no -log10(pvals)
    for i in range(8): # 8 bins corresponding to PGC = 1,2,3,4,5,[6-10],[11-30],[30-inf]
        pdf_matrix_raw.append([])

    for SNP in aligned_SNP_dict:
        PGC = int(aligned_SNP_dict[SNP][0])
        PGC_bin = pleio_bin(PGC)
        pdf_matrix_raw[PGC_bin].append(float(aligned_SNP_dict[SNP][1]))

    # Use a gwas pdf based on the actual SNPs that align with our eQTL data
    pdf_matrix = []
    i_pdf = 0
    for pvals in pdf_matrix_raw:

        # Convert to -log10 p-values
        log_pvals = -numpy.log10(pvals)
        num_snps = len(log_pvals)

        # Threshold the pvals at genome-wide significance level
        thresh_pval = 9     # the threshold p-value in -log10 space
        log_pvals = numpy.clip(log_pvals, 0, thresh_pval)

        # Obtain histogram of pvals
        bins = numpy.arange(0, thresh_pval + bin_width, bin_width)
        print ('i_pdf=' + str(i_pdf) + 'before numpy.histogram, len(bins)=' + str(len(bins)))
        pval_hist = numpy.histogram(log_pvals, bins)
        counts = pval_hist[0]
        pval_pdf = counts/float(num_snps)
        
        bins = pval_hist[1][0:-1]   # Remove the last element since plt.bar with align='edge' wants to see only the left edge of each bar
        print ('i_pdf=' + str(i_pdf) + 'after numpy.histogram, len(bins)=' + str(len(bins)))

        # Ensure we have maximum precision
        pval_pdf = numpy.array(pval_pdf, dtype=numpy.float128)

        # Add to the matrix
        pdf_matrix.append(pval_pdf)
        i_pdf += 1
       
    # Change name to be consistent with other returns
    pval_pdf = pdf_matrix

    print ('mean_p=' + str(mean_p))
    print ('var_p='  + str(var_p))
    return pdf_matrix, thresh_pval, bins, counts, mean_p, var_p

def get_conv_length(pval_pdf, pval_pdf_headzeros, SNP_pleios):
    """Return the size of the convolved pdf"""
    pdfs = []
    ###   Step1. track the length of each SNP pdf
    for curr_pleios in SNP_pleios:
        curr_pdf = pval_pdf[curr_pleios]
        pdfs.append(curr_pdf)

    shift_count = 0
    sum_length = 0
    for n in range(0,len(pdfs)):
        sum_length += len(pdfs[n])
        shift_count = shift_count + pval_pdf_headzeros[SNP_pleios[n]]

    pdf_length = shift_count + sum_length - len(pdfs) + 1
    num_conv = len(pdfs)
    return pdf_length, num_conv

### Create a proabbility distribution function (PDF) of possible scores that we can compare a *specific* score against to determine a gene's p-value.
### This function conditions the individual SNP PDFs on their pleiotropy as it convolves the PDFs
def create_score_pdf_pleiotropic(pval_pdf, pval_pdf_headzeros, bin_width, thresh_pval, SNP_pleios, mean_p, var_p, gene_score):
    #pickle.dump( pval_pdf, open( "pval_pdf.pickle", "wb" ) )
    # Threshold to use central limit approach for large numbers of loci.
    #threshold_effective_loc = 5000000000 #set to super big to disable using centrol limit to compute pval
    threshold_effective_loc = 100
    #100 #if more than this, use the central limit
    time_start = time.time()
    conv_pdf_length, num_conv = get_conv_length(pval_pdf, pval_pdf_headzeros,\
            SNP_pleios)
    #central limit pvalue calculated no matter how many convolutions we make
    newMean = num_conv * mean_p
    newVar = num_conv * var_p
    newSigma = math.sqrt(newVar)
    central_pval = 1 - stats.norm.cdf(gene_score,newMean,newSigma)
    print ("num_conv:" , num_conv)
    print ("central_pval:", central_pval)
    print ("SNP_pleios_print_out\t" + str(len(SNP_pleios)))
    if len(SNP_pleios) < threshold_effective_loc: #5000000000:
        # Need to convolve N pdfs together to get the proper distribution.  Here, we use a different PDF each time based on pleiotropy
        score_pdf = pval_pdf[SNP_pleios[0]]     # Index of the pval_pdf matrix is the pleiotropy of the SNP + 1 (i.e. correct for index starting at zero)
        #print 'PDF'
        #print str(0)
        #print score_pdfi
        pleios_count={}
        pleios_count[SNP_pleios[0]]=1
        shift_count=pval_pdf_headzeros[SNP_pleios[0]]
        #of_gene_pleios.write(str(len(SNP_pleios)) + '\n')
        if len(SNP_pleios) > 1:
            score_pdf, shift_count, pleios_count = \
                fftconvolve_ifft_after_all_convolution(SNP_pleios,pval_pdf, \
                    shift_count,pval_pdf_headzeros,pleios_count)
        # Determine maximum score based on threshold pvalue and the number of
        # convlutions
        # Note that plt.bar with align='edge' wants to see the left edge of each bar
        #shift the score_pdf to the right by adding zeros that we removed from the head of the pdf
        if shift_count > 0 :
            score_pdf = numpy.concatenate([numpy.zeros(shift_count),score_pdf])

        # Determine maximum score based on threshold pvalue and the number of convlutions
        # Note that plt.bar with align='edge' wants to see the left edge of each bar
        max_score = thresh_pval * len(SNP_pleios)
        bins = numpy.arange(0, max_score-(bin_width*(len(SNP_pleios) - 1)), bin_width)
        bins = numpy.around(bins, decimals=int(numpy.log10(1/bin_width)))

        print ('sum score_pdf convolution: ' + str(sum(score_pdf)))
        #negative pdf should be zero
        #this is necessary because of the convolution by fft produce negative
        #values for some edge instances
        #(when there are zero probability in the middle
        score_pdf[score_pdf<0]=0
        print ('sum score_pdf after correcting negatives: '+str(sum(score_pdf)))
        print ("convoluted with the following pleiotropy")
        print ("different pleiotropy: " + str(len(pleios_count.keys())))
        print ("max_score:"+ str(max_score))
        print ("score_pdf_length "+ str(len(score_pdf)))
        print ("conv_pdf_length " + str(conv_pdf_length))
        print ("SNP_pleoids"+str(len(SNP_pleios)))
        for pleios in pleios_count:
            print (str(pleios)+":"+ str(pleios_count[pleios]))
        print ("with convolution:"+ str(time.time()-time_start) + "secs")

    else:
        #do nothing but set the
        score_pdf = None
        bins = None
    return score_pdf, bins, num_conv, central_pval
#
#    elif len(SNP_pleios) > 5000000000:
#        #central limit
#        #using the mean_p and std_p
#        max_score = thresh_pval * len(SNP_pleios)
#        bins = numpy.arange(0, max_score-(bin_width*(len(SNP_pleios) - 1)), bin_width)
#        bins = numpy.around(bins, decimals=int(numpy.log10(1/bin_width)))
#
#        # Get statistics from the base distribution:
#        N = len(SNP_pleios)
#        mean = mean_p
#        sigma = std_p
#        score_pdf = 1/(sigma * numpy.sqrt(2 * numpy.pi)) * numpy.exp( - (numpy.array(bins, dtype=numpy.float128) - mean)**2 / (2 * var))
#
#        print 'sum score_pdf central limit: ' + str(sum(score_pdf))
#        print "with centrol limit:"+ str(time.time()-time_start) + "secs"
#
#    else:
#        max_score = thresh_pval * len(SNP_pleios)
#        bins = numpy.arange(0, max_score-(bin_width*(len(SNP_pleios) - 1)), bin_width)
#        bins = numpy.around(bins, decimals=int(numpy.log10(1/bin_width)))
#
#        # Get statistics from the base distribution:
#        N = len(SNP_pleios)
#        mean = N*numpy.mean(pval_pdf[0])
#        var = N*numpy.var(pval_pdf[0])
#        sigma = numpy.sqrt(var)
#        # integrate over the bin_width
#        score_pdf = bin_width*1/(sigma * numpy.sqrt(2 * numpy.pi)) * numpy.exp( - (numpy.array(bins, dtype=numpy.float128) - mean)**2 / (2 * var))
#        print max_score, bin_width, len(SNP_pleios), mean, var, sigma
#        print 'sum score_pdf central limit: ' + str(sum(score_pdf))
#        print "with centrol limit:"+ str(time.time()-time_start) + "secs"
#


### Compute scores and pvalues for each gene in the gene list ###
def compute_scores_and_pvals(gene_list, database_cursor, GWAS, bin_width, include_HLA_region):
    def thresh_test(pval):
        # This function is used to test a gene p-value against a simple threshold
        # A false result means the p-value is BETTER than the threshold
        threshold = float(1e-3)
        if pval == 'None' or pval == 'Not Loaded':
            result = True
        elif float(pval) > threshold:
            result = True
        else:
            result = False

        return result


    def compute_pval(score_pdf, gene_score, score_bins, bin_width, gene, test):
        time_start = time.time()
        # Compute score index that works in 99% of cases
        score_index = bisect.bisect(score_bins, gene_score) - 3
        if score_index < 0:
            score_index = 0

        # Now compute the p-value
        test_pval = sum(score_pdf[score_index:len(score_pdf)])
        if test_pval == 0:
            while test_pval == 0:                           # This is required due to an accumulation of rounding error.  This requires careful consideration when updating the code!!!!
                score_index = score_index - 1               # That is, few genes (1 in 1000 or fewer) may require an offset of 3, 4, or 5.
                test_pval = sum(score_pdf[score_index:len(score_pdf)])

        # Print the scoring graph if necessary (this is useful for debugging)

        print ("compute pvalue took:"+str(time.time()-time_start))
        print ("test_pval:",test_pval, "score_pdf length:", len(score_pdf))

        return test_pval


    # Identify the threshold p-values used in the sequence of tests applied.
    # Why the heck are these hard-coded?  They are intimately linked to the eQTL data sets (and how they were created) and should not be changed.
    cis_test = [1e-3]
    trans_tests = [1e-5, 5e-6, 5e-7]

    threshold_effective_loc = 50000000 #if more than this, use the central limit

    print ('begin get_GWAS_pdf')
    time_get_GWAS_pdf = time.time()
    #pdb.set_trace()
    gwas_pdf, thresh_pval, bins, counts, mean_p, var_p = get_GWAS_pdf(gene_list, GWAS, database_cursor, bin_width, include_HLA_region)
    print ('end get_GWAS_pdf Took '+str(time.time()-time_get_GWAS_pdf))
    print ("mean_p:", mean_p)
    print ("var_p:", var_p)

    #pdb.set_trace()
    gwas_pdf2=[]
    gwas_pdf_head_zeros=[]
    for pdf in gwas_pdf:
        pdf2 = numpy.trim_zeros(pdf,'bf')
        pdf3 = numpy.trim_zeros(pdf,'f')
        gwas_pdf2.append(pdf2)
        gwas_pdf_head_zeros.append(pdf.size-pdf3.size)
    gwas_pdf=gwas_pdf2

    print ("SNP_pleios_print_out\tlen(gene_list)=" + str(len(gene_list))) 
    for gene in gene_list:

    
        best_test_pval = float(1)
        best_test_snps = []
        gene.num_blocks_used = int(0)

        SNP_set_dict = {} # used to store the SNP name sets used in each of the 4 possible tests
        ################## First perform cis test ##########################
        # Get the tag SNPs that meet our threshold criteria
        tag_SNPs = [pair for pair in gene.SNP_pairs if pair.GWAS_SNP_pval != 'NA' and pair.is_tag_SNP == 'Yes' and pair.is_cis_SNP == 'cis']
        SNP_set_dict['cis'] = set([pair.eQTL_SNP for pair in tag_SNPs])
        print ("computing for ", gene.name)
        if len(tag_SNPs) > 0:
            print ("cis snp scoring")
            # Sort the SNPs
            for pair in tag_SNPs:
                pair.GWAS_SNP_pval = float(pair.GWAS_SNP_pval)
            tag_SNPs.sort(key=lambda x: x.GWAS_SNP_pval, reverse=False)

            # Convert
            SNP_pair_probs = [-numpy.log10(float(pair.GWAS_SNP_pval)) for pair in tag_SNPs]
            SNP_pleios = [int(pair.PGC) if int(pair.PGC) < int(800) else int(800) for pair in tag_SNPs]
            #pdb.set_trace()
            print ('exact SNP_pleios=' + str(SNP_pleios))
            SNP_pleios = [pleio_bin(int(pair.PGC)) for pair in tag_SNPs]
            print ('bins of SNP_pleios=' + str(SNP_pleios))
            # How many tag SNPs are we using in the test?
            num_cis_tag_SNPs = len(SNP_pair_probs)

            # Clip the SNP_pair_probs to just better than genome wide significance
            SNP_pair_probs = numpy.clip(SNP_pair_probs, 0, 9)
            # In our method, the score is simply the sum of the -10log(GWAS p-values) for each vaid tag SNP
            gene_score = sum(SNP_pair_probs)
            # Get the score pdf for a gene with the exact number of valid tag SNP blocks
            #if gene == 'ABAT':
            #    pdb.set_trace()
            #    a = 1
            score_pdf, score_bins, num_conv, central_pval = create_score_pdf_pleiotropic(gwas_pdf, gwas_pdf_head_zeros, bin_width, thresh_pval, SNP_pleios, mean_p, var_p, gene_score)
            # set to the central limit pvalue fist
            test_pval = central_pval
            # Now that we know the gene_score and score_pdf, compute the p-value
            # if score_pdf is not None, (None if using central limit theorem)
            if score_pdf is not None :
                test_pval = compute_pval(score_pdf, gene_score, score_bins, bin_width, gene.name, 'cis')

            gene.cis_test = float(test_pval)
            best_test_pval = float(test_pval)
            gene.p_a = float(test_pval)
            gene.test_used = 'cis'
            gene.snps_used = [foo.eQTL_SNP for foo in tag_SNPs]
            gene.num_blocks_used = int(num_cis_tag_SNPs)

        else:
            print ("SNP_pleios_print_out\t-1") 
            gene.cis_test = 'None'


        ################### Next perform trans tests ########################

        for t in [1,2,0]:       # We use this order to ensure that the default test is performed last

            # Get the tag SNPs that meet our threshold criteria
            tag_SNPs = [pair for pair in gene.SNP_pairs if pair.GWAS_SNP_pval != 'NA' and pair.is_tag_SNP == 'Yes' and float(pair.eQTL_SNP_pval) <= float(trans_tests[t])]
                
            if len(tag_SNPs) > 0:
                print ("trans scoring--test type:"+str(t))
                for pair in tag_SNPs:
                    pair.GWAS_SNP_pval = float(pair.GWAS_SNP_pval)
                tag_SNPs.sort(key=lambda x: x.GWAS_SNP_pval, reverse=False)
                SNP_pair_probs = [-numpy.log10(float(pair.GWAS_SNP_pval)) for pair in tag_SNPs]
                SNP_pleios = [int(pair.PGC) if int(pair.PGC) < int(800) else int(800) for pair in tag_SNPs]
                print ('exact SNP_pleios=' + str(SNP_pleios))
                #pdb.set_trace()
                SNP_pleios = [pleio_bin(int(pair.PGC)) for pair in tag_SNPs]
                print ('bins of SNP_pleios=' + str(SNP_pleios))
                
                num_trans_tag_SNPs = len(SNP_pair_probs)

                SNP_pair_probs = numpy.clip(SNP_pair_probs, 0, 9)
                gene_score = sum(SNP_pair_probs)            # In our method, the score is simply the sum of the -10log(GWAS p-values) for each vaid tag SNP


                if t == 0:
                    SNP_set_dict['trans_test_1'] = set([pair.eQTL_SNP for pair in tag_SNPs])
                    score_pdf, score_bins, num_conv, central_pval = create_score_pdf_pleiotropic(gwas_pdf, gwas_pdf_head_zeros, bin_width, thresh_pval, SNP_pleios, mean_p, var_p, gene_score)
                    test_pval = central_pval
                    if score_pdf is not None:
                        test_pval = compute_pval(score_pdf, gene_score, score_bins, bin_width, gene.name, str(t))
                    gene.trans_test_1 = float(test_pval)
                elif t == 1:
                    SNP_set_dict['trans_test_2'] = set([pair.eQTL_SNP for pair in tag_SNPs])
                    score_pdf, score_bins , num_conv, central_pval= create_score_pdf_pleiotropic(gwas_pdf, gwas_pdf_head_zeros,  bin_width, thresh_pval, SNP_pleios, mean_p, var_p, gene_score)
                    test_pval = central_pval
                    if score_pdf is not None:
                        test_pval = compute_pval(score_pdf, gene_score, score_bins, bin_width, gene.name, str(t))
                    gene.trans_test_2 = float(test_pval)

                elif t == 2:
                    SNP_set_dict['trans_test_3'] = set([pair.eQTL_SNP for pair in tag_SNPs])
                    score_pdf, score_bins, num_conv, central_pval = create_score_pdf_pleiotropic(gwas_pdf, gwas_pdf_head_zeros,  bin_width, thresh_pval, SNP_pleios, mean_p, var_p, gene_score)
                    test_pval = central_pval
                    #set to the centrol limit first
                    #update it when score_pdf being calculated
                    if score_pdf is not None:
                        test_pval = compute_pval(score_pdf, gene_score, score_bins, bin_width, gene.name, str(t))
                    gene.trans_test_3 = float(test_pval)

                if not math.isnan(test_pval) and float(test_pval) < float(best_test_pval):
                    best_test_pval = float(test_pval)
                    gene.p_a = float(test_pval)
                    gene.test_used = 'trans'
                    gene.snps_used = [foo.eQTL_SNP for foo in tag_SNPs]
                    gene.num_blocks_used = int(num_trans_tag_SNPs)


            else:
                if t == 0:
                    gene.trans_test_1 = 'None'
                    SNP_set_dict['trans_test_1'] = set()
                elif t == 1:
                    gene.trans_test_2 = 'None'
                    SNP_set_dict['trans_test_2'] = set()
                elif t == 2:
                    gene.trans_test_3 = 'None'
                    SNP_set_dict['trans_test_3'] = set()
                print ("SNP_pleios_print_out\t-1") 

        adjustment_method = 'Number of tests'
        #adjustment_method = 'None'

        if adjustment_method == 'Number of tests':
            # Adjust the p-value based on the number of tests
            # Find the number of unique tests
            unique_tests = []
            if gene.cis_test != 'None' and float(gene.cis_test) != float(1):
                unique_tests.append(float(gene.cis_test))
            if gene.trans_test_1 != 'None' and float(gene.trans_test_1) != float(1):
                unique_tests.append(float(gene.trans_test_1))
            if gene.trans_test_2 != 'None' and float(gene.trans_test_2) != float(1):
                unique_tests.append(float(gene.trans_test_2))
            if gene.trans_test_3 != 'None' and float(gene.trans_test_3) != float(1):
                unique_tests.append(float(gene.trans_test_3))
            num_unique_tests = len(set(unique_tests))

            # Order the tests
            unique_tests.sort()

            if thresh_test(gene.cis_test) and thresh_test(gene.trans_test_2) and thresh_test(gene.trans_test_3):
                if gene.trans_test_1 != 'None':
                    best_test_pval = float(test_pval)
                    gene.p_a = float(test_pval)
                    gene.test_used = 'trans'
                    gene.snps_used = [foo.eQTL_SNP for foo in tag_SNPs]
                    gene.num_blocks_used = int(num_trans_tag_SNPs)

            # Determine if there is a significant difference between the tests before testing
            # Adjust the number of unique tests for p-values that are very close.  Find the REAL number of unique tests
            for i in range(num_unique_tests):
                try:
                    if (2*unique_tests[i]) > unique_tests[i+1]:
                        del unique_tests[i+1]
                except:
                    pass

            # Correct based on the number of tests
            print ("Number of unique tests:",unique_tests)
            num_unique_tests = len(unique_tests)
            
            print ("BEGIN SNP_set_dict")

            cis = SNP_set_dict['cis']
            trans_test_1 = SNP_set_dict['trans_test_1']
            trans_test_2 = SNP_set_dict['trans_test_2']
            trans_test_3 = SNP_set_dict['trans_test_3']

            num_unique_sets = 1 #start 1 for cis test
            if len(cis) == 0  or cis == trans_test_1 or cis == trans_test_2 or cis == trans_test_3:
                num_unique_sets = num_unique_sets - 1   # when cis is the same as one of three trans tests 
            # trans tests have diff SNP sets if and only if they have diff length 
            # because the eSNP_pval thresholds are monotonically increasing
            if len(trans_test_1) > 0 :
                num_unique_sets = num_unique_sets + 1
            if len(trans_test_2) > 0 and trans_test_2 != trans_test_1:
                num_unique_sets = num_unique_sets + 1
            if len(trans_test_3) > 0  and trans_test_3 != trans_test_1 and trans_test_3 != trans_test_2:
                num_unique_sets = num_unique_sets + 1
            for test in SNP_set_dict:
                print (str(test) + ':' + str(SNP_set_dict[test]))
             
            print ("END SNP_set_dict")
            print ('gene = ' + gene.name)
            print ('before correct, num_unique_tests='+str(num_unique_tests))
            #FIXME the key for this branch
            # use unique SNP sets instead of unique pvals to define unique tests
            num_unique_tests = num_unique_sets
            print ('num_unique_tests=num_unique_sets='+str(num_unique_tests))
 
            print ("p-value", gene.p_a)
            #print thresh_test(gene.p_a)
            if thresh_test(gene.p_a) == False and (int(num_unique_tests) > int(1)):
                print ('Changed: ' + str(gene.p_a))
                gene.p_a = float(gene.p_a)*int(num_unique_tests)
                print ('To: ' + str(gene.p_a))


        if adjustment_method == 'None':
            pass


    # Lastly, purge the genes that had no overlap (no blocks met all the alignment criteria)
    output_gene_list = []
    for gene in gene_list:
        if gene.num_blocks_used != int(0):
            output_gene_list.append(gene)

    return output_gene_list






#################################################################################################################################
###  Various Helper Functions that are called by the main functions.
#################################################################################################################################

### Get GWAS SNPs from the database.  Used, for example, in buliding the GWAS PDF.
def get_GWAS_SNPs(GWAS, database_cursor, include_HLA_region='True'):
    database_cursor.execute("use GWAS")
    SQL_string = "select pval, snp from " + GWAS
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    if qr == ():
        print ("Error: no data in GWAS query")
    else:
        pvals = [qr[i][0] for i in range(len(qr))]
        rs_nums = [qr[i][1] for i in range(len(qr))]

    # Annotate the SNPs
    GWAS_SNPs = []
    for pval, rs_num in zip(pvals, rs_nums):
        tmp = SNP()
        tmp = annotate_SNP(tmp, rs_num, database_cursor)
        tmp.pval = pval
        GWAS_SNPs.append(tmp)

    # If required, set SNPs near HLA region to NA
    GWAS_SNPs_no_HLA = []
    if include_HLA_region == 'False':
        for GWAS_SNP in GWAS_SNPs:
            if GWAS_SNP.chrom == 'chr6' and int(GWAS_SNP.location) > 27000000 and int(GWAS_SNP.location) < 33000000:
                foo = 'HLA_SNP'
            else:
                GWAS_SNPs_no_HLA.append(GWAS_SNP)
        GWAS_SNPs = GWAS_SNPs_no_HLA

    # If required, set SNPs on X chromosome to NA
    GWAS_SNPs_no_X = []
    include_X_chrom = 'False'
    if include_X_chrom == 'False':
        for GWAS_SNP in GWAS_SNPs:
            if GWAS_SNP.chrom == 'chrX':
                foo = 'X chrom'
            else:
                GWAS_SNPs_no_X.append(GWAS_SNP)
        GWAS_SNPs = GWAS_SNPs_no_X

    return GWAS_SNPs


# Define a function to stitch eQTL blocks together into "pseudo genes" that we can test for statistical significance
def construct_pseudo_genes(gene_list, pseudo_gene_subgenes, pseudo_gene_names, pseudo_gene_definitions):
    tmp_pseudo_gene_list = []
    geneToGeneList = {}
    for item in gene_list:
        geneToGeneList[item.name]=item

    print (len(pseudo_gene_subgenes))
    # First, load every eQTL block we find in the genes
    for i in range(len(pseudo_gene_subgenes)):
        tmp = Gene()
        tmp.name = pseudo_gene_names[i]
        print (i, tmp.name)
        tmp.chrom = 'pseudo'
        tmp.txStart = 'pseudo'
        tmp.txEnd = 'pseduo'
        tmp.SNP_pairs = []
        if pseudo_gene_definitions != []:
            tmp.definition = pseudo_gene_definitions[i]
        for subgene in pseudo_gene_subgenes[i]:
             if subgene in geneToGeneList:
        #for item in gene_list:
        #    if item.name in pseudo_gene_subgenes[i]:
                #print 'Found ' + item.name
                # Found a matching gene in this eQTL.  Not all eQTL will have every gene.
                for pair in geneToGeneList[subgene].SNP_pairs:
                    tmp.SNP_pairs.append(pair)
        tmp_pseudo_gene_list.append(tmp)


    print ('remove pseudo genes without snps')
    # Next, toss out pseudo genes that don't have any SNPs (i.e. none of their genes were actually in the eQTL)
    updated_tmp_pseudo_gene_list = []
    for i in range(len(tmp_pseudo_gene_list)):
        if tmp_pseudo_gene_list[i].SNP_pairs != []:
            updated_tmp_pseudo_gene_list.append(tmp_pseudo_gene_list[i])
    pseudo_gene_list = updated_tmp_pseudo_gene_list
    del tmp_pseudo_gene_list

    # Remove any duplicate SNPs
    for item in pseudo_gene_list:
        updated_SNP_pairs = []
        tmp_dict = {}
        for pair in item.SNP_pairs:
            if pair.eQTL_SNP not in tmp_dict:
                updated_SNP_pairs.append(pair)
                tmp_dict[pair.eQTL_SNP] = 1
        item.SNP_pairs = updated_SNP_pairs

    # Order the SNPs by chromosome
    for item in pseudo_gene_list:
        sort_SNPs(item.SNP_pairs)

    # Attach genes to gene list and return
    for item in pseudo_gene_list:
        gene_list.append(item)
    return gene_list


# Define a function to stitch eQTL blocks together into "pseudo genes" that we can test for statistical significance
def construct_pseudo_genes_new(gene_list, pseudo_gene_subgenes, pseudo_gene_names, pseudo_gene_definitions):
    tmp_pseudo_gene_list = []
    geneToGeneList = {}
    for item in gene_list:
        geneToGeneList[item.name]=item

    print (len(pseudo_gene_subgenes))
    # First, load every eQTL block we find in the genes
    for i in range(len(pseudo_gene_subgenes)):
        tmp = Gene()
        tmp.name = pseudo_gene_names[i]
        print (i, tmp.name)
        tmp.chrom = 'pseudo'
        tmp.txStart = 'pseudo'
        tmp.txEnd = 'pseduo'
        tmp.SNP_pairs = []
        if pseudo_gene_definitions != []:
            tmp.definition = pseudo_gene_definitions[i]
        for subgene in pseudo_gene_subgenes[i]:
             if subgene in geneToGeneList:
        #for item in gene_list:
        #    if item.name in pseudo_gene_subgenes[i]:
                #print 'Found ' + item.name
                # Found a matching gene in this eQTL.  Not all eQTL will have every gene.
                for pair in geneToGeneList[subgene].SNP_pairs:
                    tmp.SNP_pairs.append(pair)
        tmp_pseudo_gene_list.append(tmp)


    print ('remove pseudo genes without snps')
    # Next, toss out pseudo genes that don't have any SNPs (i.e. none of their genes were actually in the eQTL)
    updated_tmp_pseudo_gene_list = []
    for i in range(len(tmp_pseudo_gene_list)):
        if tmp_pseudo_gene_list[i].SNP_pairs != []:
            updated_tmp_pseudo_gene_list.append(tmp_pseudo_gene_list[i])
    pseudo_gene_list = updated_tmp_pseudo_gene_list
    del tmp_pseudo_gene_list

    # Remove any duplicate SNPs
    for item in pseudo_gene_list:
        tmp_dict = {}
        for pair in item.SNP_pairs:
            # if not already included in dict, add
            if pair.eQTL_SNP not in tmp_dict :
                tmp_dict[pair.eQTL_SNP] = pair
            # if already included, update value in the case of a more significant p-value 
            else:
                existing_pair = tmp_dict[pair.eQTL_SNP]
                if float(pair.eQTL_SNP_pval) < float(existing_pair.eQTL_SNP_pval):   
                    tmp_dict[pair.eQTL_SNP] = pair
                    print ('update ' + pair.eQTL_SNP + ' pval from ' + str(existing_pair.eQTL_SNP_pval)  + ' to ' + str(pair.eQTL_SNP_pval))
        item.SNP_pairs = tmp_dict.values()

    # Order the SNPs by chromosome
    for item in pseudo_gene_list:
        sort_SNPs(item.SNP_pairs)

    new_gene_list = []
    # Attach genes to gene list and return
    for item in pseudo_gene_list:
        new_gene_list.append(item)
    return new_gene_list



def compute_qvals_python(gene_list):
    p_values = numpy.array([float(gene.pval) for gene in gene_list])
    p_values = numpy.clip(p_values,0,1)
    qvalues = qvalue_Storey_Tibshirani_2003.estimate(p_values,verbose = True)
    for gene, qvalue in zip(gene_list, qvalues):
        gene.qval = '{:02.2e}'.format(float(qvalue))
    return gene_list
    if False: # without calculating Pi0
        p_values = [float(gene.pval) for gene in gene_list]
        # the package statsmodels.stats.multitest.fdrcorrection returns two elements
        # first being a boolean list, True if a hypothesis is rejected, False if not
        # second being the q-value, or pvalues adjusted for multiple hypothesis testing to limit FDR
        qvalues = fdrcorrection(p_values)[1].tolist() 
        pdb.set_trace()
        for gene, qvalue in zip(gene_list, qvalues):
            gene.qval = '{:02.2e}'.format(float(qvalue))
        return gene_list



### Compute q-values from the full list of p-values
def compute_qvals(gene_list):

    if len(gene_list) > 200: #1000
        try:
            qvalue = importr('qvalue')
        except:
            qvalue = importr('qvalue',lib_loc='/genomesvr1/home/jgu1/R/x86_64-pc-linux-gnu-library/3.6')
        tmp = [float(gene.pval) for gene in gene_list]
        zero_threshold = float(1e-100)                       # In R, very small p-values are getting set to zero, causing problems
        for i in range(len(tmp)):
            if tmp[i] < zero_threshold:
                tmp[i] = zero_threshold
        one_threshold = float(1)                            # Rounding error may cause p-values very slightly larger than one
        for i in range(len(tmp)):
            if tmp[i] > one_threshold:
                tmp[i] = one_threshold
        #print tmp
        #print 'Length of p-values: ' + str(len(tmp))
        p_values = robjects.FloatVector(tmp)
        print ('These are p-values:')
        print (p_values)
        #pdb.set_trace()
        # before calculating qvalue
        q_value_obj = qvalue.qvalue(p_values)
        print ("This is the q-value object")
        print (q_value_obj)
        qvalues = q_value_obj.rx2('qvalues')
        pvalues = q_value_obj.rx2('pvalues')
        lambda_ = q_value_obj.rx2('lambda')

        # For debugging
        #print 'Lambda is: ' + str(lambda_)

        pi0 = q_value_obj.rx2('pi0')
        for gene, qvalue in zip(gene_list, qvalues):
            gene.qval = '{:02.2e}'.format(float(qvalue))
            if gene.qval == '9.28e-014.09e-01':
                pdb.set_trace()
                a = 1
    else:
        for gene in gene_list:
            gene.qval = 'NA'

    return gene_list


### Compute p-values for gene enrichment based on independent gene lists (e.g. MalaCards)
def top_gene_enrichment(gene_list, independent_gene_list, number_of_top_genes_to_test):

    # Only perform this test if we have a real gene list
    if len(independent_gene_list) == 1:
        pval = 'No List Provided'
    elif len(gene_list) < 1000:
        pval = 'NA'
    else:
        # Get list of all genes in the results
        master_name_list = {}
        for gene in gene_list:
            master_name_list[gene.name]=1

        # Use only those gene names that are actually found in the results
        filtered_independent_gene_list = []
        for gene in independent_gene_list:
            if gene in master_name_list:
                filtered_independent_gene_list.append(gene)
        tmp_dict = dict(zip(filtered_independent_gene_list, filtered_independent_gene_list))
        unique_filtered_independent_gene_list = tmp_dict.keys()

        # Number of independent genes in top genes
        top_count = 0
        print ("overlap:")
        for i in range(number_of_top_genes_to_test):
            if gene_list[i].name in tmp_dict:
                top_count += 1
                print (gene_list[i].name)

        # perform hypergeometic test using R
        #Example:  pval = 1 - phyper(30,700,10000-700,100)
        if False:
            r_stats = importr('stats')
            print ("overlap x:", top_count)
            print ("m:", len(unique_filtered_independent_gene_list))
            print ("n:", len(gene_list) - len(unique_filtered_independent_gene_list))
            print ("k:", number_of_top_genes_to_test)
            tmp = r_stats.phyper(top_count-1, len(unique_filtered_independent_gene_list), len(gene_list) - len(unique_filtered_independent_gene_list), number_of_top_genes_to_test)
            #pval = '{:02.2e}'.format(1 - float(robjects.default_ri2py(tmp)[0]))
            pval = '{:02.2e}'.format(1 - float(robjects.conversion.ri2py(tmp)[0]))
        # scipy.stats.hypergeom.cdf(k,M,n,N)
        # M: total of both type
        # n: total of type 1
        # N: N<=n, how many in a draw
        # k: how many got drew 
        pval = 1 - scipy.stats.hypergeom.cdf(top_count-1,len(gene_list),len(unique_filtered_independent_gene_list),number_of_top_genes_to_test)
        pval = '{:02.2e}'.format(pval)

    #print pval
    return pval


def apply_eQTL_threshold(SNP_list,eQTL_threshold = 1e-5):
    """ This function prunes a SNP list based on the eQTL significance
    and proximity
    """
    output_list = []

    #eQTL_trans_threshold = float(1e-5)
    eQTL_trans_threshold = float(eQTL_threshold)
    eQTL_cis_threshold = float(1e-3)

    for SNP in SNP_list:
        if float(SNP.eQTL_SNP_pval) <= eQTL_trans_threshold and \
                SNP.is_cis_SNP == 'trans':
            output_list.append(SNP)
        elif float(SNP.eQTL_SNP_pval) <= eQTL_cis_threshold and \
                SNP.is_cis_SNP == 'cis':
            output_list.append(SNP)

    return output_list


def is_in_LD_database2(gene_list, snpInLD, ld_range=None):
    """snp_ld_range should store all the pre-computed snp ld information
    including the position ranges
    rs10000553  chr4    169903609   1   0   169770199   169941852
    rsid<tab>chr<tab>location<tab>in_ld_db<tab>extend<tab>left<tab>right
    """
    # Identify those SNPs in the LD database
    for gene in gene_list:
        for pair in gene.SNP_pairs:
            #first set NO
            pair.in_LD_DB = 'No'
            if pair.eQTL_SNP in snpInLD:
                pair.in_LD_DB = 'Yes'
            else:
                if ld_range:
                    if pair.eQTL_SNP in ld_range:
                        if ld_range[pair.eQTL_SNP][3]=='1':
                            pair.in_LD_DB = 'Yes'
                    else:
                        chr_pos = pair.chrom + ':'+ str(pair.eQTL_SNP_location)
                        if chr_pos in ld_range:
                            if ld_range[chr_pos][3]=='1':
                                pair.in_LD_DB = 'Yes'
    return gene_list

def is_in_LD_database(gene_list, database_cursor):
    database_cursor.execute("use hg19")
    # Identify those SNPs in the LD database
    for gene in gene_list:
        for pair in gene.SNP_pairs:
            #print pair
            SQL_string = "SELECT COUNT(*) FROM LD_" + pair.chrom + " WHERE SNPA=" + "'" + pair.eQTL_SNP + "'" + " OR SNPB=" + "'" + pair.eQTL_SNP + "'" + ";"
            database_cursor.execute(SQL_string)
            qr = database_cursor.fetchall()
            #print SQL_string
            #print qr
            # Test the number of entries
            if qr != ():
                counts = int(qr[0][0])
                #print counts
                if counts == 0:
                    pair.in_LD_DB = 'No'
                else:
                    pair.in_LD_DB = 'Yes'

    return gene_list


def sort_SNPs(input):
    chrom_sort_dict = {'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22, 'chrX':23, 'chrY':24}
    # Sort eQTL SNPs by chromosome, then position
    # Input can be EITHER a SNP or SNP_pair object

    if isinstance(input[0], SNP_pair):
        #print input[0]
        input.sort(key=lambda x: x.eQTL_SNP_location, reverse=False)
        input.sort(key=lambda x: chrom_sort_dict[x.chrom], reverse=False)

    elif isinstance(input[0], SNP):
        input.sort(key=lambda x: x.location, reverse=False)
        input.sort(key=lambda x: chrom_sort_dict[x.chrom], reverse=False)

    else:
        print ("Error in sort_SNPs")

    return input


# For every eQTL SNP in a gene's SNP pairs, determine if it is cis or trans (near of far from the gene)
def determine_SNP_gene_proximity(gene_list):
    for gene in gene_list:
        for pair in gene.SNP_pairs:
            if pair.chrom == gene.chrom:
                if (pair.eQTL_SNP_location > gene.txStart - 1000000) and (pair.eQTL_SNP_location < gene.txEnd + 1000000):
                    pair.is_cis_SNP = 'cis'
                else:
                    pair.is_cis_SNP = 'trans'
            else:
                pair.is_cis_SNP = 'trans'

    return gene_list



# Use this to find all SNPs in LD at some level with every SNP in a SNP_list
def get_SNPs_in_LD(one_SNP, database_cursor, LD_threshold):
    database_cursor.execute("use hg19")

    # Define the output list
    SNPs_in_LD = []
    rsq_of_SNPs = []

    # Search with one orientation
    SQL_string = "SELECT SNPA, rsq from LD_" + one_SNP.chrom + "  WHERE SNPB='" + one_SNP.rsID + "'" + " AND rsq > " + str(LD_threshold)
    database_cursor.execute(SQL_string)
    qr1 = database_cursor.fetchall()

    # Search with the other orientation
    SQL_string = "SELECT SNPB, rsq from LD_" + one_SNP.chrom + "  WHERE SNPA='" + one_SNP.rsID + "'" + " AND rsq > " + str(LD_threshold)
    database_cursor.execute(SQL_string)
    qr2 = database_cursor.fetchall()

    # Combine
    qr = qr1 + qr2

    # Parse
    if qr != ():
        SNPs_in_LD.append(one_SNP.rsID)
        rsq_of_SNPs.append('1.000000')
        for result in qr:
            SNPs_in_LD.append(result[0])
            rsq_of_SNPs.append(result[1])

    else:
        SNPs_in_LD = one_SNP.rsID
        rsq_of_SNPs = '1.000000'

    return SNPs_in_LD, rsq_of_SNPs



def read_GWAS_catalog(gwas_catalog_file):
    """
    Use this to read a local copy (file) of the GWAS catalog and construct
    a dictionary for all genes listed
    """
    gwas_catalog_dictionary = {}
    FileReader = open(gwas_catalog_file,errors='ignore')
    moo = FileReader.readline()	# skip the header
    for line in FileReader:
        line = line.split('\t')
        if len(line) > 13:
            items = line[13].split(',')
            items = [item.strip(' ') for item in items]
            for item in items:
                if item in gwas_catalog_dictionary:
                    num_entries = len(gwas_catalog_dictionary[item])
                    publications = [gwas_catalog_dictionary[item][i][2] \
                            for i in range(num_entries)]
                    if line[5] not in publications:
                        gwas_catalog_dictionary[item].append(
                                [line[7], line[6], line[5]])
                else:
                    gwas_catalog_dictionary[item] = [
                            [line[7], line[6], line[5]]]
                    # Order of catalog entries is 1. Disease Name
                    # 2. Publication Title  3. URL to publication
    return gwas_catalog_dictionary


def load_GWAS_catalog_entries(gene_list, catalog_file_location):
    catalog_dict = read_GWAS_catalog(catalog_file_location)
    for gene in gene_list:
        #print gene
        if gene.name in catalog_dict:
            gene.GWAS_catalog = catalog_dict[gene.name]

    return gene_list

def submit_job_SGE(command_line_string, eQTL, GWAS, email_address, log_folder_location):

    # First, create a temporary BASH script for just this job.  Add some random characters to the file name, just to ensure no collisions with other instances of this program
    job_submission_script_location = '/tmp/ES_job_submission_script_' + ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(10)) + '.qsub'

    # Store the qsub script as a string, then write it to the above temporary file
    fs = ''

    fs = fs + "#!/bin/bash" + '\n'
    fs = fs + "### Request Bourne shell as shell for job" + '\n'
    fs = fs + "#$ -S /bin/bash" + "\n"

    fs = fs + "#$ -e " + log_folder_location + '\n'
    fs = fs + "#$ -o " + log_folder_location + '\n'

    fs = fs + "### Configure email" + '\n'
    fs = fs + "#$ -M " + email_address + '\n'
    #fs = fs + "#$ -m bea" + '\n'
    #not sending mail: #$ -m n
    fs += "#$ -m n" + '\n'


    fs = fs + "### Select the host and/or que" + "\n"
    fs = fs + "#$ -l hostname=genomesvr2.ucsf.edu" + "\n"

    fs = fs + "### Name the job:" + "\n"
    fs = fs + "#$ -N " + "ES_" + eQTL + "_" + GWAS + "\n"

    fs = fs + "### Use current working directory" + "\n"
    fs = fs + "#$ -cwd" + "\n"
    fs = fs + "#$ -q highmem" + "\n"


    fs = fs + command_line_string + "\n"

    FileWriter = open(job_submission_script_location, 'w')
    FileWriter.write(fs)
    FileWriter.close()

    # Must make the file executable
    cmd = 'chmod +x ' + job_submission_script_location
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)

    # Submit the Job to sge
    command_line_string = "qsub -l mem_free=6G " + job_submission_script_location
    p = Popen(command_line_string, shell=True, stdout=PIPE, stderr=STDOUT)
    foo = p.communicate()[0]        # Provides error info. if necessary

    print ('Submitted Job: eQTL ' + eQTL + ' and GWAS ' + GWAS)

def submit_job_Mac(command_line_string, eQTL, GWAS, email_address, log_folder_location):

    # Submit the Job direction to Mac command line
    p = Popen(command_line_string + ' | tee tmp.log', shell=True, stdout=PIPE, stderr=STDOUT)
    foo = p.communicate()        # Provides error info. if necessary
    print (foo[0])



def get_eQTL_meta_txt(eQTL, data_folder_location):
    """Get eQTL meta from text file
       Output format:
       [eqtl_id, eqtl_tissue, eqtl_cohort, eqtl_study_size,
        eqtl_publication_link,eqtl_publication_title]
    """
    meta_file = data_folder_location + '/eQTL/'+eQTL+'.meta'
    meta_info = []
    if os.path.isfile(meta_file):
        with open(meta_file) as meta:
            #[eqtl_id, eqtl_tissue, eqtl_cohort, eqtl_study_size,
            #    eqtl_publication_link,eqtl_publication_title]
            meta_info = meta.readlines()
    else:
        print ('No meta data found for '+eQTL)

    if len(meta_info)<6:
        print (eQTL+ ' meta data has no enought lines '+str(len(meta_info)))
    return meta_info

def get_eQTL_meta_db(eQTL, database_cursor):
    """
    Get eQTL meta info from eQTL database
    """
    database_cursor.execute("use eQTL")
    SQL_string = "Select * from eqtl_meta_data where eqtl_id='" + eQTL + "';"
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    if qr == ():
        print ("Error: no data in eQTL meta data query")
    else:
        eqtl_id = qr[0][0]
        eqtl_tissue = qr[0][1]
        eqtl_cohort = qr[0][2]
        eqtl_study_size = qr[0][3]
        eqtl_publication_link = qr[0][4]
        eqtl_publication_title = qr[0][5]
    return [eqtl_id, eqtl_tissue, eqtl_cohort, eqtl_study_size,
            eqtl_publication_link,eqtl_publication_title]


def get_GWAS_meta_txt(GWAS, data_folder_location):
    """Get GWAS meta from text file
       Output format:
        [gwas_id, gwas_disease, gwas_search_terms, gwas_publication_link,
         gwas_publication_title,ingenuity_gene_list, malacard_gene_list]
    """
    meta_file = data_folder_location + '/GWAS/'+GWAS+'.meta'
    meta_info = []
    if os.path.isfile(meta_file):
        with open(meta_file) as meta:
            #[gwas_id, gwas_disease, gwas_search_terms, gwas_publication_link,
            # gwas_publication_title,ingenuity_gene_list, malacard_gene_list]
            meta_info = [item.strip('"') for item in meta.readlines()]
    else:
        print ('No meta data found for '+GWAS)

    if len(meta_info)<7:
        print (GWAS+ ' meta data has no enought lines: '+str(len(meta_info)))
    return meta_info

def get_GWAS_meta_db(GWAS, database_cursor):
    """
    get GWAS meta data from database
    """
    database_cursor.execute("use GWAS")
    SQL_string = "Select * from gwas_meta_data where gwas_id='" + GWAS + "';"
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    if qr == ():
        print ("Error: no data in GWAS meta data query")
        gwas_id = 'none'
        gwas_disease = 'none'
        gwas_search_terms = 'none'
        gwas_publication_link = 'none'
        gwas_publication_title = 'none'
        ingenuity_gene_list = 'none'
        malacard_gene_list = 'none'
    else:
        #records in gwas_meta_data can not be null in order for the below code
        #to function properly
        gwas_id = qr[0][0].strip('"')
        gwas_disease = qr[0][1].strip('"')
        gwas_search_terms = qr[0][2].strip('"')
        gwas_publication_link = qr[0][3].strip('"')
        gwas_publication_title = qr[0][4].strip('"')
        ingenuity_gene_list = qr[0][5].strip('"')
        malacard_gene_list = qr[0][6].strip('"')

    return [gwas_id, gwas_disease, gwas_search_terms, gwas_publication_link,
            gwas_publication_title,ingenuity_gene_list, malacard_gene_list]

    #################################################################################################################################
###   Functions used to convert the raw output pickle files into html documents with summary data.
###   In addition, it generates a series of folders containing text files for each gene in each GWAS eQTL combination detailing
###   All inputs to the calculations.
################################################################################################################################

### Load every pickle file, find the gene ranking, and report results
def place_results_online(GWAS, eQTL, database_cursor, data_folder, results_folder_location, bin_width, FDR_threshold, gwas_catalog_file):

    #hold(True)

    plot_colors = ['aqua', 'black', 'blue', 'red', 'navy', 'green', 'purple', 'yellow', 'fuchsia', 'maroon',  'olive',   'lime', 'silver', 'teal', 'gray','white']

    html_files_folder_location = results_folder_location.rstrip('/') + '/html/'

    # Get data from the GWAS catalog
    catalog_dict = read_GWAS_catalog(gwas_catalog_file)
    catalog_keys = catalog_dict.keys()

    #print 'Working ' + GWAS + ' and ' + eQTL

    # Create new subfolder for each GWAS + eQTL combination
    specific_folder_location = html_files_folder_location + GWAS + '_' + eQTL

    # Prepare folder structure to hold all results
    if not os.path.exists(specific_folder_location):
        os.makedirs(specific_folder_location)
    #for letter in string.uppercase:
    #    if not os.path.exists(specific_folder_location + '/' + letter):
    #        os.makedirs(specific_folder_location + '/' + letter)
    #if not os.path.exists(specific_folder_location + '/' + 'other'):
    #    os.makedirs(specific_folder_location + '/' + 'other')          # Need folder for oddball genes that might start not start with uppercase letter (non HUGO)

    eqtl_meta = get_eQTL_meta_txt(eQTL, data_folder)
    [eqtl_id, eqtl_tissue, eqtl_cohort, eqtl_study_size,
            eqtl_publication_link,eqtl_publication_title] = ['none']*6
    if len(eqtl_meta) >=6:
        [eqtl_id, eqtl_tissue, eqtl_cohort, eqtl_study_size,
            eqtl_publication_link,eqtl_publication_title] = eqtl_meta

    # Get info from GWAS metadata
    [gwas_id, gwas_disease, gwas_search_terms, gwas_publication_link,
            gwas_publication_title,ingenuity_gene_list, malacard_gene_list] = \
                    ['none']*7
    gwas_meta = get_GWAS_meta_txt(GWAS, data_folder)
    [gwas_id, gwas_disease, gwas_search_terms, gwas_publication_link,
            gwas_publication_title,ingenuity_gene_list, malacard_gene_list] = \
                    gwas_meta

    # Place genes into a nice list
    ingenuity_gene_list = ingenuity_gene_list.split(',')
    ingenuity_gene_list = [foo.strip() for foo in ingenuity_gene_list]
    malacard_gene_list = malacard_gene_list.split(',')
    malacard_gene_list = [foo.strip() for foo in malacard_gene_list]

    file_name = results_folder_location.rstrip('/') + '/' + GWAS + '_' + eQTL + '.pickle'
    gene_list = pickle.load(open(file_name, 'rb'))

    # Perform hypergeometic test of results using MalaCards gene list
    enrichment_pval_top_30 = top_gene_enrichment(gene_list, malacard_gene_list, 30)
    enrichment_pval_top_100 = top_gene_enrichment(gene_list, malacard_gene_list, 100)

    # Create graphs for this particular GWAS and eQTL combination
    graph_1_file_name = specific_folder_location + '/' +'ES_Results_QQ_plot.png'
    color_genes = [gene.name for gene in gene_list if (gene.qval != 'NA' and float(gene.qval) < float(FDR_threshold))]
    QQ_plot_gene_pvalues(GWAS, eQTL, database_cursor, output_file_name=graph_1_file_name, gene_list=gene_list, color_genes=color_genes)

    graph_2_file_name = specific_folder_location + '/' +'GWAS_QQ_plot.png'
    GWAS_QQ_plot_local(GWAS, eQTL, data_folder , output_file_name=graph_2_file_name)

    graph_3_file_name = specific_folder_location + '/' +'location_gene_supporting_snps.png'
    plot_gene_supporting_snps(gene_list, graph_3_file_name)

    graph_4_file_name = specific_folder_location + '/' +'PGC_pvals_plot.png'
    PGC_pvals_plot(gene_list, graph_4_file_name)

    #graph_3_file_name = specific_folder_location + '/' + 'GWAS_PDF.png'
    #create_plot_GWAS_pdf(GWAS, database_cursor, bin_width, output_file_name=graph_3_file_name)

    ### Construct a summary HTML file of the results
    file_name = specific_folder_location + '/' + GWAS + '_' + eQTL + "_summary.html"
    file_name3 = specific_folder_location + '/' + GWAS + '_' + eQTL + "_gene.txt"
    FileWriter = open(file_name, 'w')
    FileWriter3 = open(file_name3, 'w')
    file_text = []

    ### Determine the number of significant genes
    number_significant_genes = 0
    for gene in gene_list:
        if gene.qval != 'NA' and float(gene.qval) < float(FDR_threshold):
            number_significant_genes += 1

    ### Determine number of genes better than a certain (0.01) p-value
    number_genes_below_01 = 0
    for gene in gene_list:
        if gene.pval != 'NA' and float(gene.pval) < float(0.01):
            number_genes_below_01 += 1

    ### Determine number of valid genes
    number_valid_genes = len(gene_list)

    ### Determine how many genes fall outside a 95 % confidence interval

    # First compute the confidence intervals
    N = number_valid_genes
    c975 = []
    c025 = []
    for i in range(1, N+1):
        #c975.append(-numpy.log10(qbeta(0.975,i,N-i+1)[0]))
        #c025.append(-numpy.log10(qbeta(0.025,i,N-i+1)[0]))
        c975.append(-numpy.log10(scipy.stats.beta.ppf(0.975,i,N-i+1)))
        c025.append(-numpy.log10(scipy.stats.beta.ppf(0.025,i,N-i+1)))


    # Keep a running count of the ordered stats that fall outside our confidence intervals
    number_genes_outside_95_confidence = 0
    for gene, high_band, low_band in zip(gene_list, c025, c975):
        if (float(-numpy.log10(float(gene.pval))) > float(high_band)) or (float(-numpy.log10(float(gene.pval))) < float(low_band)):
            number_genes_outside_95_confidence += 1
            print (-numpy.log10(float(gene.pval)), high_band, low_band)

    ### Embed comments in file that can be used to extract results overview for the main page
    file_text.append("<!--number_significant_genes " + str(number_significant_genes) + " -->")
    file_text.append("<!--malacards_30_pval " + enrichment_pval_top_30 + " -->")
    file_text.append("<!--number_valid_genes " + str(number_valid_genes) + " -->")    # Just added 1/5/14
    file_text.append("<!--number_genes_below_01 " + str(number_genes_below_01) + " -->")    # Just added 1/5/14
    file_text.append("<!--number_genes_outside_95_confidence " + str(number_genes_outside_95_confidence) + " -->")    # Just added

    # Construct the header and start of table
    file_text.append("<h3><strong>Results summary for " + gwas_disease + " GWAS using " + eqtl_tissue + " eQTL</strong></h3>")
    file_text.append("<p></p>")
    file_text.append("<p>")
    file_text.append("<b>GWAS Data Set</b><br>")
    file_text.append("&nbsp;&nbsp;&nbsp;&nbsp;ID: " + GWAS + "<br>")
    file_text.append("&nbsp;&nbsp;&nbsp;&nbsp;Disease: " + gwas_disease + "<br>")
    file_text.append("&nbsp;&nbsp;&nbsp;&nbsp;Publication: <a href=\"" + gwas_publication_link + "\" target=\"_blank\">" + gwas_publication_title + "</a><br>" )
    file_text.append("</p>")

    file_text.append("<p>")
    file_text.append("<b>eQTL Data Set</b><br>")
    file_text.append("&nbsp;&nbsp;&nbsp;&nbsp;ID: " + eQTL + "<br>")
    file_text.append("&nbsp;&nbsp;&nbsp;&nbsp;Tissue: " + eqtl_tissue + "<br>")
    file_text.append("&nbsp;&nbsp;&nbsp;&nbsp;Cohort: " + eqtl_study_size + " " + eqtl_cohort + "<br>")
    file_text.append("&nbsp;&nbsp;&nbsp;&nbsp;Publication: <a href=\"" + eqtl_publication_link + "\" target=\"_blank\">" + eqtl_publication_title + "</a><br>" )
    file_text.append("</p>")

    file_text.append("<p>")
    file_text.append("<b>Gene Set Enrichment</b> (if available)<br>")
    file_text.append("&nbsp;&nbsp;&nbsp;&nbsp;Hypergeometric test of <a href=\"http://www.malacards.org/\" target=\"_blank\">MalaCards</a> genes in top 30 genes : " + enrichment_pval_top_30 + "<br>")
    file_text.append("&nbsp;&nbsp;&nbsp;&nbsp;Hypergeometric test of <a href=\"http://www.malacards.org/\" target=\"_blank\">MalaCards</a> genes in top 100 genes: " + enrichment_pval_top_100 + "<br>")
    file_text.append("</p>")

    file_text.append("<center><div style=\"width:1200px; height:450px; padding:5px;\">")
    file_text.append(   "<image src=\"" + 'GWAS_QQ_plot.png' + "\" style=\"width:600px; height:450px; float:left;\" />")
    file_text.append(   "<image src=\"" + 'ES_Results_QQ_plot.png' + "\" style=\"width:600px; height:450px; float:right\">")
    file_text.append("</center>")
    file_text.append("<center><div style=\"width:1200px; height:400px; padding:5px;\">")
    file_text.append(   "<image src=\"" + 'location_gene_supporting_snps.png' + "\" style=\"width:1200px; height:400px; float:center\">")
    file_text.append("</center>")
    file_text.append("<center><div style=\"width:1200px; height:400px; padding:5px;\">")
    file_text.append(   "<image src=\"" + 'PGC_pvals_plot.png' + "\" style=\"width:1200px; height:400px; float:center\">")
    file_text.append("</center>")
    file_text.append("<p></p>")
    file_text.append("<table border=\"1\" cellpadding=\"2\">")
    file_text.append("<tbody>")
    file_text.append("<tr><td><strong>Rank</strong><td><strong>Gene Symbol</strong></td><td><strong>Significant</strong></td><td><strong>Type</strong></td><td><strong>Location</strong></td><td><strong>Gene Definition</strong></td><td><strong>GWAS Catalog Gene</strong></td><td><strong>Ingenuity Gene</strong></td><td><strong>MalaCard Gene</strong></td><td><strong>p-value</strong></td><td><strong>q-value</strong></td></tr>")

    # Write to File
    FileWriter.write("\n".join(file_text))
    file_text = []

    for gene in gene_list:

        # Establish gene significance flags
        if gene.qval == 'NA':
            significant = str('No')
            font_color = str('000000')
        else:
            if float(gene.qval) < float(FDR_threshold):
                significant = str('Yes')
                font_color = str('FF0000')
            else:
                significant = str('No')
                font_color = str('000000')

        # Construct gene
        location = str(gene.chrom) + ":" + str(gene.txStart) + "-" + str(gene.txEnd)

        # Use the first letter of the HUGO gene symbol to define the folder that we store these in (don't want 20k files in one folder)
        if gene.name[0] in string.ascii_uppercase:
            gene_html_file = specific_folder_location + "/" + gene.name[0] + "/" + gene.name.replace("/", "_") + '.txt'
            gene_html_file_hyperlink = "./" + gene.name[0] + "/" + gene.name.replace("/", "_") + '.txt'
        else:
            gene_html_file = specific_folder_location + "/" + 'other' + "/" + gene.name.replace("/", "_") + '.txt'
            gene_html_file_hyperlink = "./" + 'other' + "/" + gene.name.replace("/", "_") + '.txt'

        gene_hyperlink = "http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + gene.name
        ucsc_hyperlink = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=" + location

        if gene.name in ingenuity_gene_list:
            ingenuity_gene = 'Yes'
        else:
            ingenuity_gene = 'No'

        if gene.name in catalog_keys:
            gwas_catalog_gene = 'Yes'
            gwas_catalog_diseases = [foo[0] for foo in catalog_dict[gene.name]]
            gwas_catalog_diseases = '&#10;'.join(gwas_catalog_diseases)
        else:
            gwas_catalog_gene = 'No'
            gwas_catalog_diseases = ''

        if gene.name in malacard_gene_list:
            malacard_gene = 'Yes'
        else:
            malacard_gene = 'No'

        line = "<tr><td>" + str(gene.rank) + "</td>  <td style=\"text-align: left;\"><a href=\"" + gene_hyperlink  + "\" target=\"_blank\">" + gene.name + "</a></td>  <td><FONT COLOR=\"" + font_color  + "\">" + significant + "</FONT></td> <td>" + gene.test_used + "</td> <td> <a href=\"" + ucsc_hyperlink  + "\" target=\"_blank\">" + location + "</a></td>  <td>" + gene.definition + "</td>  <td title=\"" + gwas_catalog_diseases + "\">" + gwas_catalog_gene + "</td>  <td>" + ingenuity_gene + "</td>  <td>" + malacard_gene + "</td>  <td style=\"text-align: left;\"><a href=\"" + gene_html_file_hyperlink  + "\" target=\"_blank\">" + '{:02.2e}'.format(float(gene.pval))  + "</a></td>  <td>" + gene.qval + "</td></tr>"

        #  Whatever line was created, now append it to the master file
        FileWriter.write(line + '\n')

        # Construct the individual gene results file (details all the SNPs used for each result)
        #FileWriter2 = open(gene_html_file, 'w')
        #FileWriter2.write(str(gene))
        #FileWriter2.close()
        FileWriter3.write('#Begin Gene\n')
        FileWriter3.write(str(gene))
        FileWriter3.write('\n#End Gene\n')

    # End the table
    file_text.append("</tbody>")
    file_text.append("</table>")

    # Write to File
    FileWriter.write("\n".join(file_text))
    file_text = []
    FileWriter.close()
    FileWriter3.close()

##################################################################################################################
###   Function used to calculate the pdf against which a gene with specific p-value is compared
###   Suppose this gene has n matching SNPs between eQTL and GWAS, the result pdf is the convolution of each SNP's pdf
###   Step1. track the length of each SNP pdf
###   Step2. calculate an uniform fft length 'l_uniform' that will be used in fft for each SNP pdf
###   Step3. compute the fft for each SNP pdf.
###   Step4. multiply these fft together. Because 'l_uniform' is used, all fft have the same length
###   Step5. do inverse fft
##################################################################################################################
def fftconvolve_ifft_after_all_convolution(SNP_pleios,pval_pdf,shift_count,pval_pdf_headzeros,pleios_count):
    pdfs = []
    ###   Step1. track the length of each SNP pdf
    for curr_pleios in SNP_pleios:
        curr_pdf = pval_pdf[curr_pleios]
        pdfs.append(curr_pdf)
    ###   Step2. calculate an uniform fft length 'l_uniform' that will be used in fft for each SNP pdf
    fshape,fslice = calc_fshape_fslice_arr(pdfs)
    fft = rfft(pdfs[0],fshape)
    for n in range(1,len(pdfs)):
    ###   Step3. compute the fft for each SNP pdf.
    ###   Step4. multiply these fft together. Because 'l_uniform' is used, all fft have the same length
        fft = fft * rfft(pdfs[n], fshape)
        shift_count = shift_count + pval_pdf_headzeros[SNP_pleios[n]]
        if SNP_pleios[n] in pleios_count:
            pleios_count[SNP_pleios[n]]=pleios_count[SNP_pleios[n]]+1
        else:
            pleios_count[SNP_pleios[n]]=1

    ###   Step5. do inverse fft
    convolve_result = irfft(fft,fshape)[fslice]
    return convolve_result,shift_count,pleios_count

##################################################################################################################
### Function used to calculate an uniform fft  length 'l_uniform' that each pdf will use in rfft
### Step1. Track the length of each SNP and add them together, get result 'sum'
### Step2. Calculate the next smallest regular number that is larger than 'sum'. regular number means be a multiple of 2,3 or 5. Regular Number is good for doing fft
### Step. Return the regular number 'fshape' and a slice 'fslice' which is just (0,'sum')
##################################################################################################################
def calc_fshape_fslice_arr(twoDArr):

    ### Step1. Track the length of each SNP and add them together, get result 'sum'
    shape = 1     # in the case of three arrays, it should be (s1-1) + (s2-1) + (s3-1) + 1
    for arr in twoDArr:
        shape = shape + len(arr) -1
    ### Step2. Compute the next smallest regular number that is larger than 'sum'. regular number means be a multiple of 2,3 or 5. Regular Number is good for doing fft
    fshape = _next_regular(shape)#
    fslice = [slice(0, shape)]    #
    ### Step. Return the regular number 'fshape' and a slice 'fslice' which is just (0,'sum')
    return (fshape,fslice)

##################################################################################################################
### Function used to calculate the next smallest regular number that is larger than 'sum'. regular number means be a multiple of 2,3 or 5. Regular Number is good for doing fft
##################################################################################################################
def _next_regular(target):
    """
    Find the next regular number greater than or equal to target.
    Regular numbers are composites of the prime factors 2, 3, and 5.
    Also known as 5-smooth numbers or Hamming numbers, these are the optimal
    size for inputs to FFTPACK.

    Target must be a positive integer.
    """
    if target <= 6:
        return target

    # Quickly check if it's already a power of 2
    if not (target & (target-1)):
        return target

    match = float('inf') # Anything found will be smaller
    p5 = 1
    while p5 < target:
        p35 = p5
        while p35 < target:
            # Ceiling integer division, avoiding conversion to float
            # (quotient = ceil(target / p35))
            quotient = -(-target // p35)

            # Quickly find next power of 2 >= quotient
            try:
                p2 = 2**((quotient - 1).bit_length())
            except AttributeError:
                # Fallback for Python <2.7
                p2 = 2**(len(bin(quotient - 1)) - 2)

            N = p2 * p35
            if N == target:
                return N
            elif N < match:
                match = N
            p35 *= 3
            if p35 == target:
                return p35
        if p35 < match:
            match = p35
        p5 *= 5
        if p5 == target:
            return p5
    if p5 < match:
        match = p5
    return match

def print_gene_pval_txt(data_folder_location, output_folder_location, GWAS, eQTL):
    output_file_name_root = os.path.join(output_folder_location , GWAS + '_' + eQTL)
    gene_list = pickle.load(open(output_file_name_root + '.pickle','rb'))
    HGNC_gene_set = set()
    for line in open('./lib/all_HGNC.txt','r'):
        HGNC_gene_set.add(line.strip())
    of = open(output_file_name_root + '.gene_pval.tsv','w+')
    for gene in gene_list:
        if gene.name in HGNC_gene_set:
            of.write('{}\t{}\n'.format(gene.name,gene.pval))

    # Indicate that the analysis of one eQTL + GWAS is complete
