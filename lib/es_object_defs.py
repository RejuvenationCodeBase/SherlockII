# Objects used in the Empirical Sherlock algorithm (Gene, SNP, SNP_pair, etc.) are defined here

####################################################################################################################
###   Define Gene Object and annotation function
####################################################################################################################

def getHUGODB(database_cursor):
    database_cursor.execute("use hg19")
    SQL_string = "SELECT gd_app_sym, gd_app_name, gd_pub_ensembl_id FROM HUGO"
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    hugoDB = {}
    for r in qr:
        hugoDB[r[0]]=[r[1],r[2]]
    return hugoDB

def getGeneSymbolToID(database_cursor):
    SQL_string = "select geneSymbol, kgID from kgXref"
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    geneSymbolToID = {}
    for r in qr:
        if not (r[0] in geneSymbolToID):
            #only store the first encountered geneSymbol
            geneSymbolToID[r[0]]=r[1]
    return geneSymbolToID

def getGeneIDPos(database_cursor):
    geneIDPos = {}
    SQL_string = "select name, chrom, txStart, txEnd from knownGene"
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    for r in qr:
        geneIDPos[r[0]]=[r[1],r[2],r[3]]
    return geneIDPos

def getSNPsInLD(database_cursor, snpAnnotation):
    database_cursor.execute("use hg19")
    chrom_sort_dict = {'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22, 'chrX':23, 'chrY':24}
    chrSNPs = {}
    for snp in snpAnnotation:
        if snpAnnotation[snp][0] in chrom_sort_dict:
            chrSNPs.setdefault(snpAnnotation[snp][0],[])
            chrSNPs[snpAnnotation[snp][0]].append(snp)
    # store any snp with LD data
    snpInLD = {}
    # Identify the SNPs in the LD database
    queryTime = 0
    total = len(snpAnnotation)/10000
    for chrom in chrSNPs:
        names = []
        count = 0
        for snp in chrSNPs[chrom]:
            names.append(snp)
            count += 1
            if count>=10000:
                count=0
                querySNPInLD(database_cursor, chrom, names, snpInLD)
                queryTime += 1
                print ("Query SNP in LD for 10000 names:" +chrom+'\t'+ str(queryTime) + '/'+str(total))
                names=[]
        if names:
            querySNPInLD(database_cursor, chrom, names, snpInLD)
    return snpInLD

def querySNPInLD(database_cursor, chrom, names, snpInLD):
    nameTuple = tuple(names)
    hitted = {}
    SQL_stringA = "SELECT SNPA FROM LD_" + chrom + " WHERE SNPA in" + str(nameTuple)
    SQL_stringB = "SELECT SNPB FROM LD_" + chrom + " WHERE SNPB in" + str(nameTuple)
    database_cursor.execute(SQL_stringA)
    qr = database_cursor.fetchall()
    for r in qr:
        if not (r[0] in hitted):
            hitted[r[0]]=1
    #name2 need to be search if first attemp missed
    name2 = []
    for name in names:
        if not (name in hitted):
            name2.append(name)
    name2Tuple = tuple(name2)
    SQL_stringB = "SELECT SNPB FROM LD_" + chrom + " WHERE SNPB in" + str(name2Tuple)
    database_cursor.execute(SQL_stringB)
    qr = database_cursor.fetchall()
    for r in qr:
        if not (r[0] in hitted):
            hitted[r[0]] = 1
    for snp in hitted:
        snpInLD[snp]=1

def querySNPAnnotation(database_cursor, names, snpAnnotation):
    """Query annotations for a list of snps
    names: list of snp_ids
    snpAnnotation: dictionary to store the information for each snp
    """
    nameTuple = tuple(names)
    #only get information when the chromosome is in chrom_sort_dict
    #some snp_ids has multiple records (usually on a varient chromosome like
    #chr17-var-
    #make sure we only keep the records with the standard chromosome name
    #defined in chrom_sort_dict
    chrom_sort_dict = {'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22, 'chrX':23, 'chrY':24}

    database_cursor.execute("use hg19")
    SQL_string = "SELECT name, chrom, chromEnd, func, alleleFreqs FROM snp147 WHERE name in " + str(nameTuple)
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    for r in qr:
        if r[1] in chrom_sort_dict:
            snpAnnotation[r[0]]=[r[1],r[2],r[3]]
            try:
                foo = float(r[4].split(',')[0])
                boo = float(r[4].split(',')[1])
                try:
                    goo = float(r[4].split(',')[2])
                except:
                    goo = 'only two alleles'

                if goo != 'only two alleles':
                    MAF = sorted([foo, boo, goo])[1]
                else:
                    MAF = min([foo, boo])

            except:
                MAF = 'none'
            snpAnnotation[r[0]].append(MAF)


def getSNPAnnotation(database_cursor, allSNPs):
    snpAnnotation = {}
    names = []
    count = 0
    queryTime = 0
    total = len(allSNPs)/10000
    for snp in allSNPs:
        names.append(snp)
        count += 1
        if count>=10000:
            count=0
            querySNPAnnotation(database_cursor, names, snpAnnotation)
            queryTime += 1
            #print "Query Time for 10000 names:" + str(queryTime) + '/'+str(total)
            names=[]
    if names:
        querySNPAnnotation(database_cursor, names, snpAnnotation)
    return snpAnnotation


def annotate_Gene2(self, name, hugoDB, geneSymbolToID, geneIDPos):
    if name in hugoDB:
        self.name = name
        self.definition = hugoDB[name][0]
        self.ensembl = hugoDB[name][1]
    else:
        self.name = name
        self.definition = "Not Found in HUGO"

    chrom,start,end = ["Not Found"]*3
    if name in geneSymbolToID:
        kgID = geneSymbolToID[name]
        if kgID in geneIDPos:
            chrom = geneIDPos[kgID][0]
            start = geneIDPos[kgID][1]
            end = geneIDPos[kgID][2]

    self.chrom = chrom
    self.txStart = start
    self.txEnd = end

    return self


def annotate_Gene(self, name, database_cursor):
    """Used to load annotaiton information for a gene"""
    database_cursor.execute("use hg19")
    SQL_string = "SELECT gd_app_name, gd_pub_ensembl_id FROM HUGO WHERE gd_app_sym=" + "'" + name + "'"
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    if qr == ():
        self.name = name
        self.definition = "Not Found in HUGO"
    else:
        self.name = name
        self.definition = qr[0][0]
        self.ensembl = qr[0][1]

    # Get Gene location from UCSC known genes
    SQL_string = "select chrom, txStart, txEnd from knownGene where name=(select kgID from kgXref where geneSymbol = '" + name + "' limit 1)"
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    if qr == ():
        self.chrom = "Not Found"
        self.txStart = "Not Found"
        self.txEnd = "Not Found"
    else:
        self.chrom = qr[0][0]
        self.txStart = qr[0][1]
        self.txEnd = qr[0][2]

    return self

    # Need to get precise gene location from another database...


class Gene():
    """The Gene class contains all attributes of a given Gene"""
    def __init__(self, name='Not Loaded', rank='Not Loaded', cis_test='Not Loaded', trans_test_1='Not Loaded', trans_test_2='Not Loaded', trans_test_3='Not Loaded', test_used='Not Loaded', snps_used='Not Loaded', pleio_count='Not Loaded', p_a='Not Loaded', p_b='Not Loaded', p_b_given_a='Not Loaded', p_a_given_b='Not Loaded', pval='Not Loaded', qval='Not Loaded', eQTL='Not Loaded', GWAS='Not Loaded', definition='Not Loaded', score='Not Loaded', accession='Not Loaded', ensembl='Not Loaded', UCSCknown='Not Loaded', chrom='Not Loaded', txStart='Not Loaded', txEnd='Not Loaded', strand='Not Loaded', url='Not Loaded', SNP_pairs='Not Loaded', GWAS_catalog='Not Loaded', num_SNPs='Not Loaded', num_blocks_used='Not Loaded'):
        self.name = name
        self.rank = rank
        self.cis_test = cis_test
        self.trans_test_1 = trans_test_1
        self.trans_test_2 = trans_test_2
        self.trans_test_3 = trans_test_3
        self.test_used = test_used
        self.snps_used = snps_used
        self.pleio_count = pleio_count
        self.p_a = p_a
        self.p_b = p_b
        self.p_b_given_a = p_b_given_a
        self.p_a_given_b = p_a_given_b
        self.pval = pval
        self.qval = qval
        self.eQTL = eQTL
        self.GWAS = GWAS
        self.definition = definition
        self.score = score
        self.accession = accession
        self.ensembl = ensembl
        self.UCSCknown = UCSCknown
        self.chrom = chrom
        self.txStart = txStart
        self.txEnd = txEnd
        self.strand = strand
        self.url = url
        self.SNP_pairs = SNP_pairs
        self.GWAS_catalog = GWAS_catalog
        self.num_SNPs = num_SNPs
        self.num_blocks_used = num_blocks_used



    # Define a method for printing Gene objects
    def __str__(self):
        prt_str = '\n'
        prt_str = prt_str + '{:<}'.format("Gene: " + self.name + '\n')
        prt_str = prt_str + '{:<}'.format("Rank: " + str(self.rank) + '\n')
        if self.definition == 'Not Loaded':
            prt_str = prt_str + '{:<}'.format('\t' + "Gene definition not found in HUGO" + '\n')
        else:
            prt_str = prt_str + '{:<}'.format('\t' + self.definition + '\n')
        prt_str = prt_str + '{:<}'.format('\t' + "Has " + str(self.num_SNPs) + " in eQTL " + self.eQTL + '\n')
        if self.chrom == 'pseudo':
            # This is a pseudo gene, lacking a defined location
            prt_str = prt_str + '{:<}'.format('\t' + "Location: " + "Pseudo Gene" + '\n')
        else:
            prt_str = prt_str + '{:<}'.format('\t' + "Location: " + str(self.chrom) + ":" + str(self.txStart) + "-" + str(self.txEnd) + '\n')
        prt_str = prt_str + '{:<}'.format('\t' + "Accession: " + self.accession + '\n')
        prt_str = prt_str + '{:<}'.format('\t' + "Ensembl: " + self.ensembl + '\n')
        prt_str = prt_str + '{:<}'.format('\t' + "UCSCknown: " + self.UCSCknown + '\n')

        # Print GWAS catalog entries
        if self.GWAS_catalog == 'Not Loaded':
            prt_str = prt_str + '{:<}'.format('\t' + "This gene has no GWAS catalog entries:" + '\n')
        else:
            prt_str = prt_str + '{:<}'.format('\t' + "This gene has " + str(len(self.GWAS_catalog)) + " GWAS catalog entries:" + '\n')
            for entry in self.GWAS_catalog:
                prt_str = prt_str + '{:<}'.format('\t\t' + entry[0] + '\n')

        # Separate gene and GWAS alignment sections
        prt_str = prt_str + '\n'

        prt_str = prt_str + '{:<}'.format('\t' + "GWAS: " + self.GWAS + '\n')

        # Print the cis test results
        if self.cis_test != 'Not Loaded' and self.cis_test != 'None':
            prt_str = prt_str + '\t' + "Cis Test p-value: " + '{:02.2e}'.format(float(self.cis_test)) + '\n'
        else:
            prt_str = prt_str + '\t' + "Cis Test p-value: " + str(self.cis_test) + '\n'

        # Print the trans test 1 results
        if self.trans_test_1 != 'Not Loaded' and self.trans_test_1 != 'None':
            prt_str = prt_str + '\t' + "Trans Test 1 p-value: " + '{:02.2e}'.format(float(self.trans_test_1)) + '\n'
        else:
            prt_str = prt_str + '\t' + "Trans Test 1 p-value: " + str(self.trans_test_1) + '\n'

        # Print the trans test 2 results
        if self.trans_test_2 != 'Not Loaded' and self.trans_test_2 != 'None':
            prt_str = prt_str + '\t' + "Trans Test 2 p-value: " + '{:02.2e}'.format(float(self.trans_test_2)) + '\n'
        else:
            prt_str = prt_str + '\t' + "Trans Test 2 p-value: " + str(self.trans_test_2) + '\n'

        # Print the trans test 3 results
        if self.trans_test_3 != 'Not Loaded' and self.trans_test_3 != 'None' :
            prt_str = prt_str + '\t' + "Trans Test 3 p-value: " + '{:02.2e}'.format(float(self.trans_test_3)) + '\n'
        else:
            prt_str = prt_str + '\t' + "Trans Test 3 p-value: " + str(self.trans_test_3) + '\n'

        # Print the information for the test actually selected
        if self.test_used != 'Not Loaded' and self.num_blocks_used != 'Not Loaded':
            prt_str = prt_str + '\t' + "The " + str(self.test_used) + " test was used with " + str(self.num_blocks_used) + '\n'
        else:
            prt_str = prt_str + '\t' + "Test used Not Loaded" + '\n'

        # Print the pleiotropic count for this alignment
        prt_str = prt_str + '\t' + "Pleiotropic Count: " + str(self.pleio_count) + '\n'


        # Print the elements for the Bayesian Pleiotropic correction
        if self.p_a != 'Not Loaded':
            prt_str = prt_str + '\t' + "P(A) = " + '{:02.2e}'.format(float(self.p_a)) + '\n'
        else:
            prt_str = prt_str + '\t' + "P(A) = " + str(self.p_a) + '\n'
        if self.p_b != 'Not Loaded':
            prt_str = prt_str + '\t' + "P(B) = " + '{:02.2e}'.format(float(self.p_b)) + '\n'
        else:
            prt_str = prt_str + '\t' + "P(B) = " + str(self.p_b) + '\n'
        if self.p_b_given_a != 'Not Loaded':
            prt_str = prt_str + '\t' + "P(B|A) = " + '{:02.2e}'.format(float(self.p_b_given_a)) + '\n'
        else:
            prt_str = prt_str + '\t' + "P(B|A) = " + str(self.p_b_given_a) + '\n'
        if self.p_a_given_b != 'Not Loaded':
            prt_str = prt_str + '\t' + "P(A|B) = " + '{:02.2e}'.format(float(self.p_a_given_b)) + '\n'
        else:
            prt_str = prt_str + '\t' + "P(A|B) = " + str(self.p_a_given_b) + '\n'

        # Add a space
        prt_str = prt_str + '\n'

        # Print the p-value for this alignment
        if self.pval != 'Not Loaded':
            prt_str = prt_str + '\t' + "Final Gene p-value: " + '{:02.2e}'.format(float(self.pval)) + '\n'
        else:
            prt_str = prt_str + '\t' + "Final Gene p-value: " + str(self.pval) + '\n'


        # Print the q-value for this alignment
       	if self.qval != 'Not Loaded':
            prt_str = prt_str + '\t' + "With q-value: " + self.qval + '\n'
       	else:
            prt_str = prt_str + '\t' + "With q-value: " + str(self.qval) + '\n'

        # Separate alignment and SNP pairs section
        prt_str = prt_str + '\n'

        # Print header for the SNP pairs
        header = ''
        header = header + '\t'
        header = header + '{:<8}'.format('Block')
        header = header + '{:<8}'.format('Tag?')
        header = header + '{:<8}'.format('In_DB')
        header = header + '{:<12}'.format('Proximity')
        header = header + '{:<8}'.format('Chrom.')
        header = header + '{:<12}'.format('eQTL_SNP')
        header = header + '{:<8}'.format('PGC')
        header = header + '{:<12}'.format('Location')
        header = header + '{:<12}'.format('MAF')
        header = header + '{:<12}'.format('GWAS_SNP')
        header = header + '{:<12}'.format('Location')
        header = header + '{:<12}'.format('MAF')
        header = header + '{:<10}'.format('Delta')
        header = header + '{:<5}'.format('LD') + '   '
        header = header + '{:<8}'.format('eQTL_pval')
        header = header + '{:<8}'.format('GWAS_pval')
        prt_str = prt_str + header + '\n'

        # Print out the SNP pairs for each block below
        if self.SNP_pairs != 'Not Loaded':
            for pair in self.SNP_pairs:
                prt_str = prt_str + str(pair) + '\n'
        else:
            prt_str = prt_str + '\t' + "SNP Pair Blocks: None Loaded" + '\n'

        # Finally, return the full prt_str
        return prt_str

    def ucsc_url(self):
        # Create a URL to this location in the UCSC Genome Browser
        self.url = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=' + self.assembly + '&position=' + self.chrom + ':' + str(self.txStart) + '-' + str(self.txEnd)


####################################################################################################################
###   Define SNP Objects and annotation function
####################################################################################################################

def annotate_SNP2(self, name, snpAnnotation):
    """Load annotation from snpAnnotation"""
    if name in snpAnnotation:
        self.rsID = name
        self.chrom = snpAnnotation[name][0]
        self.location = snpAnnotation[name][1]
        self.function = snpAnnotation[name][2]
        self.MAF = snpAnnotation[name][3]
    else:
        self.rsID = name
        self.chrom = "missing"
        self.location = "missing"
        self.function = "missing"
        self.MAF = "missing"
    return self


def annotate_SNP(self, name, database_cursor):
    """Used to load annotation information for a SNP"""
    database_cursor.execute("use hg19")
    SQL_string = "SELECT chrom, chromEnd, func, alleleFreqs FROM snp135 WHERE name =" + "'" + name + "'"
    database_cursor.execute(SQL_string)
    qr = database_cursor.fetchall()
    if qr == ():
        self.rsID = name
        self.chrom = "missing"
        self.location = "missing"
        self.function = "missing"
        self.MAF = "missing"
    else:
        #sometimes there will be multiple records with the same snp name
        #usually coming from a varient chromosome like chr17-cat-hap1
        #use the first one (hopefully in form of chr17 not chr17-cat...)
        self.rsID = name
        self.chrom = qr[0][0]
        self.location = qr[0][1]
        self.function = qr[0][2]
        #print qr[0][3]
        try:
            foo = float(qr[0][3].split(',')[0])
            boo = float(qr[0][3].split(',')[1])

            try:
                goo = float(qr[0][3].split(',')[2])
            except:
                goo = 'only two alleles'

            if goo != 'only two alleles':
                MAF = sorted([foo, boo, goo])[1]
            else:
                MAF = min([foo, boo])

        except:
            MAF = 'none'
        self.MAF = MAF
    return self


class SNP():
    def __init__(self, rsID='Not Loaded', chrom='Not Loaded', location='Not Loaded', strand='Not Loaded', MAF='Not Loaded', function='Not Loaded', description='Not Loaded', pval='Not Loaded'):
        self.rsID = rsID
        self.chrom = chrom
        self.location = location
        self.strand = strand
        self.MAF = MAF
        self.function = function
        self.description = description
        self.pval = pval

    #Define a method for printing each SNP
    def __str__(self):
            prt_str = ''
            prt_str = str(self.rsID) + "\t" + str(self.chrom) + "\t" + str(self.location) + "\t" + str(self.MAF) + "\t" + str(self.function)
            if self.pval != 'Not Loaded':
                prt_str = prt_str + "\t" + str(self.pval)
            return prt_str



####################################################################################################################
###   Define SNP_pair Object
####################################################################################################################

# This class holds a pair of aligned SNPs:  one from GWAS, and one from eQTL

class SNP_pair:
    def __init__(self, block_number="Not Loaded", chrom="Not Loaded", eQTL_SNP='Not Loaded', eQTL_SNP_location='Not Loaded', eQTL_SNP_MAF='Not Loaded', PGC='NA', in_LD_DB='Not Loaded', is_cis_SNP='Not Loaded', is_tag_SNP='Not Loaded', GWAS_SNP='Not Loaded', GWAS_SNP_location='Not Loaded', GWAS_SNP_MAF='Not Loaded', separation="Not Loaded", LD_rsq='Not Loaded', eQTL_SNP_pval='Not Loaded', GWAS_SNP_pval="Not Loaded"):
        self.block_number = block_number
        self.chrom = chrom
        self.eQTL_SNP = eQTL_SNP
        self.eQTL_SNP_location = eQTL_SNP_location
        self.eQTL_SNP_MAF = eQTL_SNP_MAF
        self.PGC = PGC                  # This stands for Pleiotropic Gene Count (PGC)
        self.in_LD_DB = in_LD_DB
        self.is_cis_SNP = is_cis_SNP
        self.is_tag_SNP = is_tag_SNP
        self.GWAS_SNP = GWAS_SNP
        self.GWAS_SNP_location = GWAS_SNP_location
        self.GWAS_SNP_MAF = GWAS_SNP_MAF
        self.separation = separation
        self.LD_rsq = LD_rsq
        self.eQTL_SNP_pval = eQTL_SNP_pval
        self.GWAS_SNP_pval = GWAS_SNP_pval


    def map(self):
        """Find location information for SNPs in the pair."""

        try:
            separation = abs(self.GWAS_SNP_location - self.eQTL_SNP_location)
        except(ValueError, TypeError):
            separation = 'NA'
        else:
            separation = abs(self.GWAS_SNP_location - self.eQTL_SNP_location)

        return separation


    #Define a method for printing each SNP pair
    def __str__(self):
        self.separation = self.map()
        prt_str = ''
        prt_str = prt_str + '\t'
        prt_str = prt_str + '{:<8}'.format(str(self.block_number))
        prt_str = prt_str + '{:<8}'.format(self.is_tag_SNP)
        prt_str = prt_str + '{:<8}'.format(self.in_LD_DB)
        prt_str = prt_str + '{:<12}'.format(self.is_cis_SNP)
        prt_str = prt_str + '{:<8}'.format(self.chrom)
        prt_str = prt_str + '{:<12}'.format(self.eQTL_SNP)
        prt_str = prt_str + '{:<8}'.format(self.PGC)
        prt_str = prt_str + '{:<12}'.format(self.eQTL_SNP_location)
        prt_str = prt_str + '{:<12}'.format(self.eQTL_SNP_MAF)
        prt_str = prt_str + '{:<12}'.format(self.GWAS_SNP)
        prt_str = prt_str + '{:<12}'.format(self.GWAS_SNP_location)
        prt_str = prt_str + '{:<12}'.format(self.GWAS_SNP_MAF)
        prt_str = prt_str + '{:<10}'.format(self.separation)
        if self.LD_rsq != 'NA' and self.LD_rsq != 'Not Loaded' and self.LD_rsq != 'NF':
            prt_str = prt_str + '{:02.2f}'.format(float(self.LD_rsq)) + '   '
        else:
            prt_str = prt_str + '{:<5}'.format(self.LD_rsq) + '  '
        prt_str = prt_str + '{:10.4e}'.format(float(self.eQTL_SNP_pval)) + '   '
        if self.GWAS_SNP_pval != 'NA' and self.GWAS_SNP_pval != 'Not Loaded' and self.GWAS_SNP_pval != 'NF':
            prt_str = prt_str + '{:10.4e}'.format(float(self.GWAS_SNP_pval))
        else:
            prt_str = prt_str + '{:<8}'.format(self.GWAS_SNP_pval)
        return prt_str


