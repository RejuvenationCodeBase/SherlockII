#! /usr/bin/env python

# Version 17: Last Edit January 23, 2014.
# Contact:  chriskfuller@gmail.com  or chris@genome.ucsf.edu

# Version v10: Last Edit M, 2013.
# This is the script run from the command line that actually launches the analysis of one eQTL and one GWAS, with requisite parameters

# Import items
import sys, os, math, random, numpy, time, configparser, pdb
sys.path.append('./lib')
import pickle,pymysql
from subprocess import PIPE, Popen, STDOUT
from scipy import stats
from scipy.misc import *
from es_toolbox import *
from es_object_defs import *
from es_graphics import *
random.seed()   # Uses system time
# Organize incomming stuff
#print arg_list
config = configparser.ConfigParser()
GWAS = sys.argv[1]
eQTL = sys.argv[2]
config.read(sys.argv[3])
data_folder_location = config.get('IO','data_folder_location') + '/'
output_folder_location = config.get('IO','output_folder_location') + '/'
if not os.path.exists(output_folder_location):
    os.makedirs(output_folder_location)
try:
    align_LD_threshold = float(config.get('parameter','align_LD_threshold'))
    LD_window_cutoff = float(config.get('parameter','LD_window_cutoff'))
    eQTL_threshold = float(config.get('parameter','eQTL_threshold'))
    bin_width = float(config.get('parameter','bin_width'))
    FDR_threshold = float(config.get('parameter','FDR_threshold'))
except ValueError:
    print('please set numeric parameters correctly')
    exit(1)

gwas_catalog_file = os.path.join(data_folder_location,'gwas_catalog_v1.0-associations_e98_r2019-12-16.tsv')

print ('GWAS = ' + str(GWAS))
##################### Begin the run for ONE eQTL and ONE GWAS: ################################
time_start = time.time()

### Obtain list of gene names for this run
gene_names = get_gene_list_txt('ALL', eQTL, data_folder_location)

### Load pre-processed eQTL blocks and keep only those in gene_names list ###
pickle_file = data_folder_location + 'eQTL/' + eQTL + '.pickle'
print ('Loading ' + pickle_file)
gene_list = pickle.load(open(pickle_file, 'rb'))
updated_gene_list = []
for i in range(len(gene_list)):
    if gene_list[i].name in gene_names:
        updated_gene_list.append(gene_list[i])
del gene_list
gene_list = updated_gene_list
del updated_gene_list

# Create Pseudo Genes and add them to the gene list
# Treat a group of gene in Kegg as a pseudo gene

### Get entries from the GWAS catalog for each gene ###
gene_list = load_GWAS_catalog_entries(gene_list,
        catalog_file_location=gwas_catalog_file)

print ('Running eQTL ' + eQTL + ' and GWAS ' + GWAS)

# Create output file name root.  This is used for three output file formats: a results gene list pickle, text representations of top hits, and QQ plot image
output_file_name_root = os.path.join(output_folder_location , GWAS + '_' + eQTL)

# Find the GWAS SNPs in highest LD with the eQTL tag SNPs.  For cases where the LD is identical (e.g. several SNPs with rsq = 1), use the closest SNPs
t_load_eQTL_LD02 = time.time()
eSNP_pos, eSNP_lds, eSNP_linked = get_eSNP_in_LD02(eQTL,
        data_folder_location,LD_window_cutoff)
print ('Completed loading LD for ' + eQTL +' Took ' + \
        str(time.time()-t_load_eQTL_LD02))

t_eQTL_GWAS_alignment = time.time()
gene_list, eQTLSNP_LD = match_eQTL_SNPs_with_GWAS(eQTL, gene_list,
        None, align_LD_threshold, GWAS, data_folder_location)
print ('Completed ' + eQTL + ' and ' + GWAS + ' alignment. Took ' + \
        str(time.time()-t_eQTL_GWAS_alignment))

print ('Length of gene_list is:' + str(len(gene_list)))

# Determine which SNPs are the tag SNPs
t_assign_eQTL_tag=time.time()
ld_range = get_eSNP_LD_range(data_folder_location)
gene_list = assign_eQTL_tag_SNPs_preload(gene_list, eQTLSNP_LD,
        align_LD_threshold, eSNP_pos, eSNP_lds, eSNP_linked, ld_range,eQTL_threshold)
#gene_list = assign_eQTL_tag_SNPs(gene_list, database_cursor, LD_threshold)
#database_cursor = pymysql.connect(host="localhost", # your host, usually localhost
#            user="root",
#            passwd="genome",
#            db="hg19",
#            unix_socket='/var/run/mysqld/mysqld.sock').cursor() 
#pdb.set_trace()
#gene_list = assign_eQTL_tag_SNPs(gene_list, database_cursor, 0.3) 
print ('Completed tag SNP identification(tabix). Took '+ str(time.time()-
        t_assign_eQTL_tag))

print ('Length of gene_list is:' + str(len(gene_list)))

time_assign_pleiotropic_gene_count = time.time()
print ('begin assign_pleiotropic_gene_count')
# Determine the pleiotripic number for each aligned SNP
gene_list = assign_pleiotropic_gene_count(gene_list)
print ('finished assign_pleiotropic_gene_count Took '+str(time.time()-time_assign_pleiotropic_gene_count))
# Compute scores and p-values for each gene in the list
t_pvalue=time.time()

gene_list = compute_scores_and_pvals(gene_list, None, GWAS, bin_width, False)
print ('Completed ' + eQTL + ' and ' + GWAS + ' scoring. Took '+str(time.time()-t_pvalue))

time_qval = time.time()

print ('Length of gene_list is:' + str(len(gene_list)))

# This is a temporary line to correct for not using pleiotropy correction
for gene in gene_list:
    gene.p_a_given_b = gene.p_a
    gene.pval = gene.p_a

print ('Length of gene_list is:' + str(len(gene_list)))

# Now rank the genes according to p-value
gene_list = sorted(gene_list, key=lambda x: x.pval, reverse=False)
for i in range(len(gene_list)):
    gene_list[i].rank = i + 1

# Now assign q-values for each gene in the list (check that we have enough genes)
gene_list = compute_qvals_python(gene_list)
print ('Completed qval, Took:' + str(time.time()-time_qval) )

# Save the full gene list to pickle file
pickle.dump(gene_list, open(output_file_name_root + '.pickle', 'wb'))
#  NOTEL: to reload, use pickle.load(open('gene_list_pickle', 'r'))
print ('Computation done. Used:' + str(time.time()-time_start) + 'seconds')

# Now use the pickle file to generate html output that's human readable
place_results_online(GWAS, eQTL, None, data_folder_location, output_folder_location, bin_width, FDR_threshold, gwas_catalog_file)



# Finally generate a minimum text output containing only gene_name and corresponding p-value 
print_gene_pval_txt(data_folder_location, output_folder_location, GWAS, eQTL)

# Indicate that the analysis of one eQTL + GWAS is complete
print ('Completed ' + eQTL + ' and ' + GWAS + ' analysis')


