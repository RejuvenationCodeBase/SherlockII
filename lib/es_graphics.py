#! /usr/bin/python

# This is the library for graphical tools used in Empirical Sherlock (ES)

# Import graphic-related items.  All other imports should happen in the master script
import matplotlib, numpy, pdb
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from matplotlib import pyplot
from matplotlib import patches
from matplotlib.pyplot import Rectangle
from es_object_defs import *
import scipy.stats
#import rpy2.robjects as robjects


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



# Create plot of top gene supporting SNP locations
def plot_gene_supporting_snps(gene_list, output_file_name):

    if len(gene_list) < 500:
        sub_gene_list = gene_list[0:len(gene_list)]
    else:
        sub_gene_list = gene_list[0:500]

    supporting_snps = []
    for gene in sub_gene_list:
        for pair in gene.SNP_pairs:
            if (pair.is_tag_SNP == 'Yes') and (float(pair.GWAS_SNP_pval) < float(1e-1)):
                supporting_snps.append(pair)

    # Sort the SNPs
    #print 'Gene Supporting SNPs for the top 500 Genes:'
    supporting_snps = sort_SNPs(supporting_snps)
    #for item in supporting_snps:
    #    print item

    bin_size = 100000

    # Scale x location of our plots based on chromosome location
    chrom_lengths = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566]
    chrom_lengths = [int(foo/float(bin_size)) for foo in chrom_lengths]
    cumulative_distance = [0]       # Chromosome 1 has no cumulative distance. Chr2 has cumulative distance of chr1.  Chr3 has chr1 + chr2, and so on.
    label_location = [int((chrom_lengths[0])/float(2))]
    for chrom in range(1,23):
        cumulative_distance.append(int(sum(chrom_lengths[0:chrom])))
        label_location.append(int(sum(chrom_lengths[0:chrom]) + chrom_lengths[chrom]/float(2)))   # Subtract half so we will know the location of the chromosome number label

    # Construct histogram. Remember to scale everything by bin_size
    array_length = int(round(sum(chrom_lengths[0:22]))) + 1
    counts = [0 for i in range(array_length)]
    for item in supporting_snps:
        if len(item.chrom) <= 5 and item.chrom != 'chrX' and item.chrom != 'chrY':
            x_coordinate = int(round(cumulative_distance[int(item.chrom.strip('chr')) - 1] + int(round(item.eQTL_SNP_location/float(bin_size)))))
            counts[x_coordinate] = counts[x_coordinate] + 1

    indices = [i for i in range(array_length)]

    plot_height = int(round(max(counts) + float(0.05)* max(counts)))

    fig, ax = pyplot.subplots()
    ax.plot(indices, counts)
    fig.set_size_inches(12,4)

    # Create a plot with bars representing the chromosome boundaries
    for i in range(1, len(cumulative_distance)):
        if i % 2 == 0:
            color = 'blue'
        else:
            color = 'green'
        width = cumulative_distance[i] - cumulative_distance[i-1]
        ax.add_patch(Rectangle((cumulative_distance[i-1],0), width, plot_height, color=color, alpha=0.2))

    # Create actual histogram of top SNP Locations
    ax.plot(indices, counts, 'k')
    pyplot.title('SNP Locations for Top ' + str(len(sub_gene_list)) + ' Genes', fontsize=12)
    pyplot.xlabel('Genomic Location (chromosomes)', fontsize=12)
    pyplot.ylabel('Number of SNPs', fontsize=12)
    pyplot.xlim([0,array_length])

    ax.set_xticks(label_location)
    ax.set_xticklabels([foo for foo in range(1,len(label_location)+1)], fontsize=8)

    # Save it
    pyplot.savefig(output_file_name)

    # Clear it
    pyplot.clf()


# Create plot of top gene supporting SNP locations
def create_gwas_manhattan_plot(gwas, database_cursor, output_file_name):

    # This will generate a manhattan plot of the average p-values in each geomic block for input GWAS
    # Get stats and full pvals from the database for specified GWAS
    database_cursor.execute("use GWAS")
    SQL_string = "select pval, chrom, location from " + gwas
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    if qr == ():
        print ("Error: no data in GWAS query")
    else:
        pvals = [qr[i][0] for i in range(len(qr))]
        chroms = [qr[i][1] for i in range(len(qr))]
        locations = [int(qr[i][2]) for i in range(len(qr))]

    print ('Number pvals: ' + str(len(pvals)))
    print ('Number chroms:' + str(len(chroms)))
    print ('Number locations:' + str(len(locations)))

    # Convert to minus log10 p-values
    pvals = -numpy.log10(pvals)

    # Scale x location of our plots based on chromosome location
    chrom_lengths = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566]
    chrom_lengths = [int(foo/float(1000000)) for foo in chrom_lengths]
    cumulative_distance = [0]       # Chromosome 1 has no cumulative distance. Chr2 has cumulative distance of chr1.  Chr3 has chr1 + chr2, and so on.
    label_location = [int((chrom_lengths[0])/float(2))]
    for chrom in range(1,23):
        cumulative_distance.append(int(sum(chrom_lengths[0:chrom])))
        label_location.append(int(sum(chrom_lengths[0:chrom]) + chrom_lengths[chrom]/float(2)))   # Subtract half so we will know the location of the chromosome number label

    # Construct histogram. Remember to scale everything by 1000000
    array_length = int(round(sum(chrom_lengths[0:22]))) + 1
    pval_array = [[] for i in range(array_length)]
    for pval, chrom, location in zip(pvals, chroms, locations):
        if len(chrom) <= 5 and chrom != 'chrX' and chrom != 'chrY':
            x_coordinate = int(round(cumulative_distance[int(chrom.strip('chr')) - 1] + int(round(location/float(1000000)))))
            pval_array[x_coordinate].append(pval)

    indices = [i for i in range(array_length)]

    # Covert the list of p-values in each bin to an AVERAGE p-value per bin
    avg_pval_array = []
    for item in pval_array:
        if len(item) == 0:
            avg_pval_array.append(0)
        else:
            avg_pval_array.append(float(sum(item))/len(item))

    plot_height = int(round(max(avg_pval_array) + float(0.5) * max(avg_pval_array)))
    #plot_height = int(round(max(avg_pval)))

    fig, ax = pyplot.subplots()
    ax.plot(indices, avg_pval_array)
    fig.set_size_inches(12,4)

    # Create a plot with bars representing the chromosome boundaries
    for i in range(1, len(cumulative_distance)):
        if i % 2 == 0:
            color = 'blue'
        else:
            color = 'green'
        width = cumulative_distance[i] - cumulative_distance[i-1]
        ax.add_patch(Rectangle((cumulative_distance[i-1],0), width, plot_height, color=color, alpha=0.2))

    # Create actual histogram of top SNP Locations
    ax.plot(indices, avg_pval_array, 'k')
    pyplot.title('Mean Bin Manhattan Plot for ' + str(gwas) + ' GWAS', fontsize=12)
    pyplot.xlabel('Genomic Location (chromosomes)', fontsize=12)
    pyplot.ylabel('Mean -log10 p-value', fontsize=12)
    pyplot.xlim([0,array_length])

    ax.set_xticks(label_location)
    ax.set_xticklabels([foo for foo in range(1,23)], fontsize=8)

    # Save it
    pyplot.savefig(output_file_name)

    # Clear it
    pyplot.clf()





# Create GWAS QQ Plot
def GWAS_QQ_plot(GWAS, eQTL, database_cursor, output_file_name):

    ##### First generate the QQ plot for input GWAS
    # Get stats and full pvals from the database for specified GWAS
    database_cursor.execute("use GWAS")
    SQL_string = "select pval from " + GWAS
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    if qr == ():
        print ("Error: no data in GWAS query")
    else:
        pvals = [qr[i][0] for i in range(len(qr))]

    # Sort based on the p-values
    foo = -numpy.log10(pvals)
    obs_pvals = numpy.flipud(numpy.sort(foo))
    threshold = 12      # Display -log(p-val) greater than threshold as threshold
    obs_pvals = numpy.clip(obs_pvals, 0, threshold)

    # Create the Null Distribution
    N = len(obs_pvals)
    null_pvals = -numpy.log10(numpy.r_[1.:(N+1)]/N)

    # Create actual GWAS QQ Plot
    pyplot.plot([0., max(null_pvals)],[0., max(null_pvals)], 'k')
    pyplot.plot(null_pvals,obs_pvals, '.b')
    pyplot.title('QQ Plot for Input GWAS', fontsize=12)
    pyplot.xlabel('Null Distribution (-log p-values)', fontsize=12)
    pyplot.ylabel('SNP Significance (-log p-values)', fontsize=12)
    pyplot.ylim(0, threshold + 1)

    # Save it
    pyplot.savefig(output_file_name)

    # Clear it
    pyplot.clf()

# Create GWAS QQ Plot
def GWAS_QQ_plot_local(GWAS, eQTL, data_folder, output_file_name):
    gwas_file=data_folder+"GWAS/"+GWAS+".txt"
    pvals=[]
    with open(gwas_file) as gwas_input:
        for line in gwas_input:
            tmp=line.split('\t')
            pval = float(tmp[1])
            pvals.append(pval)

    ##### First generate the QQ plot for input GWAS
    # Get stats and full pvals from the database for specified GWAS
    #database_cursor.execute("use GWAS")
    #SQL_string = "select pval from " + GWAS
    #database_cursor.execute(SQL_string)
    #qr = database_cursor.fetchall()
    #if qr == ():
    #    print "Error: no data in GWAS query"
    #else:
    #    pvals = [qr[i][0] for i in range(len(qr))]

    # Sort based on the p-values
    foo = -numpy.log10(pvals)
    obs_pvals = numpy.flipud(numpy.sort(foo))
    threshold = 12      # Display -log(p-val) greater than threshold as threshold
    obs_pvals = numpy.clip(obs_pvals, 0, threshold)

    # Create the Null Distribution
    N = len(obs_pvals)
    null_pvals = -numpy.log10(numpy.r_[1.:(N+1)]/N)

    # Create actual GWAS QQ Plot
    pyplot.plot([0., max(null_pvals)],[0., max(null_pvals)], 'k')
    pyplot.plot(null_pvals,obs_pvals, '.b')
    pyplot.title('QQ Plot for Input GWAS', fontsize=12)
    pyplot.xlabel('Null Distribution (-log p-values)', fontsize=12)
    pyplot.ylabel('SNP Significance (-log p-values)', fontsize=12)
    pyplot.ylim(0, threshold + 1)

    # Save it
    pyplot.savefig(output_file_name)

    # Clear it
    pyplot.clf()



# Create Empirical Sherlock Results QQ Plot
def QQ_plot_gene_pvalues(GWAS, eQTL, database_cursor, output_file_name, gene_list, color_genes):

    ##### Now generate the QQ plot for the gene list
    # Sort based on the p-values
    gene_list = sorted(gene_list, key=lambda x: x.pval, reverse=False)

    # Get the ordered list of -log10 p-values
    obs_pvals = []
    for item in gene_list:
        obs_pvals.append(-numpy.log10(item.pval))
    threshold = 6      # Display -log(p-val) greater than threshold as threshold
    obs_pvals = numpy.clip(obs_pvals, 0, threshold)

    # Create the Null Distribution
    N = len(obs_pvals)
    null_pvals = -numpy.log10(numpy.r_[1.:(N+1)]/N)

    # Create the confidence intervals
    c975 = []
    c025 = []
    for i in range(1, N+1):
        #c975.append(-numpy.log10(qbeta(0.975,i,N-i+1)[0]))
        #c025.append(-numpy.log10(qbeta(0.025,i,N-i+1)[0]))
        c975.append(-numpy.log10(scipy.stats.beta.ppf(0.975,i,N-i+1)))
        c025.append(-numpy.log10(scipy.stats.beta.ppf(0.025,i,N-i+1)))


    # Create appropriate list for specific genes we want to color
    gene_number = 0
    color_genes_x_val = []
    color_genes_y_val = []
    for item in gene_list:
        gene_number = gene_number + 1
        if item.name in color_genes:
            color_genes_x_val.append(null_pvals[gene_number - 1])
            color_genes_y_val.append(obs_pvals[gene_number -1])

    # Plot
    pyplot.plot([0., max(null_pvals)],[0., max(null_pvals)], 'k')
    pyplot.plot(null_pvals,obs_pvals, '.b')
    pyplot.plot(color_genes_x_val, color_genes_y_val, 'ro')
    pyplot.plot(null_pvals,c975, 'r')
    pyplot.plot(null_pvals,c025, 'r')
    pyplot.title('QQ Plot for Output Gene List', fontsize=12)
    pyplot.xlabel('Null Distribution (-log p-values)', fontsize=12)
    pyplot.ylabel('Gene Significance (-log p-values)', fontsize=12)
    pyplot.ylim(0, threshold + 1)

    # Show it
    pyplot.savefig(output_file_name)

    # Clear it
    pyplot.clf()

# Create Empirical Sherlock Results Uniform Plot
def uniform_plot_gene_pvalues(GWAS, eQTL, database_cursor, output_file_name, gene_list, color_genes):

    ##### Now generate a plot of the pvalues and compare to uniform distribution
    # Sort based on the p-values
    gene_list = sorted(gene_list, key=lambda x: x.pval, reverse=False)

    # Get the ordered list of -log10 p-values
    obs_pvals = []
    for item in gene_list:
        obs_pvals.append(item.pval)

    # Create the Null Distribution
    N = len(obs_pvals)
    null_pvals = numpy.r_[1.:(N+1)]/N

    # Create appropriate list for specific genes we want to color
    gene_number = 0
    color_genes_x_val = []
    color_genes_y_val = []
    for item in gene_list:
        gene_number = gene_number + 1
        if item.name in color_genes:
            color_genes_x_val.append(null_pvals[gene_number - 1])
            color_genes_y_val.append(obs_pvals[gene_number -1])

    # Plot
    pyplot.plot([0., max(null_pvals)],[0., max(null_pvals)], 'k')
    pyplot.plot(null_pvals,obs_pvals, '.b')
    pyplot.plot(color_genes_x_val, color_genes_y_val, 'ro')
    pyplot.title('QQ Plot for Output Gene List', fontsize=12)
    pyplot.xlabel('Null Distribution (-log p-values)', fontsize=12)
    pyplot.ylabel('Gene Significance (-log p-values)', fontsize=12)
    pyplot.ylim(0, 1)

    # Show it
    pyplot.savefig(output_file_name)

    # Clear it
    pyplot.clf()


# Create double QQ Plot, with color option
def QQ_plot_gene_scores(GWAS, database_cursor, eQTL, gene_list, color_genes, output_file_name):


    ##### First generate the QQ plot for input GWAS
    # Get stats and full pvals from the database for specified GWAS
    database_cursor.execute("use GWAS")
    SQL_string = "select pval from " + GWAS
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    if qr == ():
        print ("Error: no data in GWAS query")
    else:
        pvals = [qr[i][0] for i in range(len(qr))]

    # Sort based on the p-values
    foo = -numpy.log10(pvals)
    obs_pvals = numpy.flipud(numpy.sort(foo))
    threshold = 12      # Display -log(p-val) greater than threshold as threshold
    obs_pvals = numpy.clip(obs_pvals, 0, threshold)

    # Create the Null Distribution
    N = len(obs_pvals)
    null_pvals = -numpy.log10(numpy.r_[1.:(N+1)]/N)

    # Plot
    pyplot.subplot(211)
    pyplot.plot([0., max(null_pvals)],[0., max(null_pvals)], 'k')
    pyplot.plot(null_pvals,obs_pvals, '.b')
    pyplot.title('QQ Plot for GWAS: ' + GWAS, fontsize=12)
    #pyplot.xlabel('Null Distribution (-log p-values)', fontsize=12)
    pyplot.ylabel('Results (-log p-values)', fontsize=12)
    pyplot.ylim(0, threshold + 1)

    ##### Now generate the QQ plot for the gene list
    # Sort based on the p-values
    gene_list = sorted(gene_list, key=lambda x: x.pval, reverse=False)

    # Get the ordered list of -log10 p-values
    obs_pvals = []
    for item in gene_list:
        obs_pvals.append(-numpy.log10(item.pval))
    threshold = 6      # Display -log(p-val) greater than threshold as threshold
    obs_pvals = numpy.clip(obs_pvals, 0, threshold)

    # Create the Null Distribution
    N = len(obs_pvals)
    null_pvals = -numpy.log10(numpy.r_[1.:(N+1)]/N)

    # Create appropriate list for specific genes we want to color
    gene_number = 0
    color_genes_x_val = []
    color_genes_y_val = []
    for item in gene_list:
        gene_number = gene_number + 1
        if item.name in color_genes:
            color_genes_x_val.append(null_pvals[gene_number - 1])
            color_genes_y_val.append(obs_pvals[gene_number -1])

    # Plot
    pyplot.subplot(212)
    pyplot.plot([0., max(null_pvals)],[0., max(null_pvals)], 'k')
    pyplot.plot(null_pvals,obs_pvals, '.b')
    pyplot.plot(color_genes_x_val, color_genes_y_val, 'ro')
    pyplot.title('QQ Plot Gene Results for: ' + eQTL, fontsize=12)
    pyplot.xlabel('Null Distribution (-log p-values)', fontsize=12)
    pyplot.ylabel('Results (-log p-values)', fontsize=12)
    pyplot.ylim(0, threshold + 1)

    #pyplot.show()
    pyplot.savefig(output_file_name)

    # Clear it
    pyplot.clf()


def create_plot_pleio_matrix(pleio_matrix, output_file_name):

    x_data = []
    y_data = []
    for item in pleio_matrix:
        # pleio p-value is first, gene p-value is next
        y_data.append(item[0])
        x_data.append(item[1])

    # Print
    print (x_data)
    print (y_data)

    # First plot the counts
    pyplot.subplot(1,1,1)
    pyplot.loglog(x_data, y_data, '.')
    pyplot.xlabel('Gene pvals')
    pyplot.ylabel('Pleiotropic Alignment pvals')
    pyplot.title('Relationship Between Gene and Pleiotropic Significance')

    #pyplot.show()
    pyplot.savefig(output_file_name)

    # Clear it
    pyplot.clf()



def create_plot_GWAS_pdf(GWAS, database_cursor, bin_width, include_HLA_region, output_file_name, pdf_type):

    from es_toolbox import get_GWAS_pdf

    # Get the necessary data for the GWAS
    pdf_type = 'full'
    gene_list = []
    pval_pdf, thresh_pval, bins, counts = get_GWAS_pdf(gene_list, GWAS, database_cursor, bin_width, include_HLA_region, pdf_type)

    # First plot the counts
    pyplot.subplot(2,1,1)
    pyplot.bar(bins, counts, width=bin_width, log=True, align='edge')
    pyplot.xlabel('SNP pvals (-log10 units)')
    pyplot.ylabel('Number of SNPs')
    pyplot.title('SNP Counts and Probabilities for ' + GWAS)

    # Then plot the pdf
    pyplot.subplot(2,1,2)
    pyplot.bar(bins, pval_pdf, width=bin_width, log=True, align='edge')
    pyplot.xlabel('SNP pvals (-log10 units')
    pyplot.ylabel('Probability of SNPs')
    #pyplot.title('SNP p-value PDF')

    #pyplot.show()
    pyplot.savefig(output_file_name)

    # Clear it
    pyplot.clf()



def make_color_gene_list(gene_list, option='HLA'):
    color_gene_list = []

    if option == 'HLA':
        # Construct a list of HLA genes from file
        file = '/Users/chrisfuller/Documents/Code/vfa/data/HLA_gene_list.txt'
        HLA_genes_in_file = []
        FileReader = open(file, 'r')
        for line in FileReader:
            HLA_genes_in_file.append(line.strip())

        for item in gene_list:
            if item.name in HLA_genes_in_file:
                color_gene_list.append(item.name)

    if option == 'Top Genes':
        pval_threshold = float(1e-3)
        for item in gene_list:
            if item.pval < pval_threshold:
                 color_genes_list_of_genes.append(item.name)
        print ('Number of Genes Better than Threshold:' + str(len(Top_Genes_list_of_genes)))

    return color_gene_list


def PGC_pvals_plot(gene_list, file_name):

    # Import items
    from matplotlib import pyplot
    import numpy

    pvals = []
    PGC = []

    for gene in gene_list:
        for pair in gene.SNP_pairs:
            if pair.is_tag_SNP == 'Yes':
                pvals.append(float(pair.GWAS_SNP_pval))
                PGC.append(float(pair.PGC))

    # Convert to log pvals
    pvals = -numpy.log10(pvals)

    # Clip for plotting
    pvals = numpy.clip(pvals, 0, 6)
    PGC = numpy.clip(PGC, 0, 500)

    # Plot
    print (PGC)
    print (pvals)
    print (str(len(PGC)))
    print (str(len(pvals)))
    pyplot.plot(PGC, pvals, 'ro')

    # Now finalize the plot before we save it
    pyplot.title('Plot of SNP p-values versus PGC', fontsize=12)
    pyplot.xlabel('PGC (number of genes)', fontsize=12)
    pyplot.ylabel('SNP Significance (-log p-values)', fontsize=12)
    pyplot.ylim(0, 6)
    pyplot.xlim(0, 500)

    # Save it
    pyplot.savefig(file_name)

    # Clear it
    pyplot.clf()


def overlap_QQ_plots(list_of_gene_lists, output_file_name, GWAS):

    for gene_list in list_of_gene_lists:

        # Sort based on the p-values
        gene_list = sorted(gene_list, key=lambda x: x.pval, reverse=False)

        # Get the ordered list of -log10 p-values
        obs_pvals = []
        for item in gene_list:
            obs_pvals.append(-numpy.log10(item.pval))
        threshold = 6      # Display -log(p-val) greater than threshold as threshold
        obs_pvals = numpy.clip(obs_pvals, 0, threshold)

        # Create the Null Distribution
        N = len(obs_pvals)
        null_pvals = -numpy.log10(numpy.r_[1.:(N+1)]/N)


        # Plot
        pyplot.plot([0., max(null_pvals)],[0., max(null_pvals)], 'k')
        pyplot.plot(null_pvals,obs_pvals, '.b')
        pyplot.title('QQ Plots Summary for GWAS: ' + GWAS, fontsize=12)
        pyplot.xlabel('Null Distribution (-log p-values)', fontsize=12)
        pyplot.ylabel('Gene Significance (-log p-values)', fontsize=12)
        pyplot.ylim(0, threshold + 1)

    # Now compute confidence interval
    c975 = []
    c025 = []
    for i in range(1, N+1):
        #c975.append(-numpy.log10(qbeta(0.975,i,N-i+1)[0]))
        #c025.append(-numpy.log10(qbeta(0.025,i,N-i+1)[0]))
        c975.append(-numpy.log10(scipy.stats.beta.ppf(0.975,i,N-i+1)))
        c025.append(-numpy.log10(scipy.stats.beta.ppf(0.025,i,N-i+1)))

    # Add confidence intervals
    pyplot.plot(null_pvals,c975, 'r')
    pyplot.plot(null_pvals,c025, 'r')

    # Save it
    pyplot.savefig(output_file_name)

    # Clear it
    pyplot.clf()


