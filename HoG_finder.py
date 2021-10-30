import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import time
# usage: python horizontal_gene_finder.py -i file.txt -c 4 -d file2.txt -o file3.png
"""
starttime = time.time()

parser=argparse.ArgumentParser()

# define your arguments
parser.add_argument('-i', '--input_file', default="genomic_GC.txt")
parser.add_argument('-c', '--confidence', default=3)
parser.add_argument('-s', '--style_sheet', default='C:\\Users\\mattk\\Desktop\\BME163\\Assignments\\BME163')#'BME163')#
parser.add_argument('-o', '--output_file', default="horizontal_gene_visualized.png")
parser.add_argument('-d', '--output_details', default='horizontal_gene_details.txt')
parser.add_argument('-x1', '--xlim_low', default=None)
parser.add_argument('-x2', '--xlim_high', default=None)

# there are multiple ways you can access your arguments
# you can put all of the arguments into one variable like this:
args = parser.parse_args()
# which can then be accessed like this:
input_file = args.input_file
confidence = float( args.confidence )
style_sheet = args.style_sheet
output_file = args.output_file
details_file = args.output_details
input_xlim_low = args.xlim_low
input_xlim_high = args.xlim_high

try:
	plt.style.use(style_sheet)
except:
	print(f'No style sheet detected at {style_sheet}')

#######################################################################
#######################################################################

figureWidth=10# inches
figureHeight=5 # inches

plt.figure(figsize=(figureWidth,figureHeight))

panelWidth=8
panelHeight=3

relativePanelWidth=panelWidth/figureWidth
relativePanelHeight=panelHeight/figureHeight


##########################################################################
########### Panel creations ##############################################

#             left,bottom,width,height
panel=plt.axes([0.1,0.2,relativePanelWidth,relativePanelHeight])

##########################################################################
################## Data Reader ###########################################
def readData(input_file):
	readListX = []
	readListY = []
	with open(input_file, 'r') as inFile:
		for numLines, line in enumerate( inFile.readlines() ) :
			if line[0] != '#':
				data = line.replace(' ', '').rstrip().split('\t')
				genomicPosition = float( data[0] )
				GCcontent = float( data[1] )

				readListX.append(genomicPosition)
				readListY.append(GCcontent)

		return readListX, readListY


##########################################################################
################## Main ##################################################
##
# GC content
##

data = readData(input_file)
xData = data[0]
yData = data[1]


##########################################################################
########### Panel adjustments ##############################################
if input_xlim_low is None:
	xlim_low = min(xData)
else:
	xlim_low = float(input_xlim_low)

if input_xlim_high is None:
	xlim_high = max(xData)
else:
	xlim_high = float(input_xlim_high)


#xlim_low = min(xData) ; xlim_high = max(xData)
#xlim_low = 120000 ; xlim_high = 140000
ylim_low = 0.8*min(yData) ; ylim_high = 1.2*max(yData)

panel.set_xlim(xlim_low,xlim_high)
panel.set_ylim(ylim_low,ylim_high)

panel.tick_params(bottom=True, labelbottom=True,
	left=True, labelleft=True,
	right=False, labelright=False,
	top=False, labeltop=False)

panel.set_xticks(np.linspace(xlim_low, xlim_high, 20))

panel.set_xlabel('\nGenomic Location', fontsize = 12)
panel.set_ylabel('GC Content %\n', fontsize = 12)

panel.plot(xData,yData, marker='o',markerfacecolor=(0,0,0),
	markeredgecolor='black', markersize=0, markeredgewidth=0,
	linewidth=1, color = (0,0,0))
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Want to plot the median and std off of it

median = np.median(yData)
std = np.std(yData)
std = confidence*std
begToEnd = [xlim_low, xlim_high]

#median
panel.plot(begToEnd, 2*[median], lw=2, color=(0,0,0))
#std 1
panel.plot(begToEnd,2*[median+std], dashes=[1,2,2,2], lw=1.2, color=(0,0,0))
#std 2
panel.plot(begToEnd,2*[median-std], dashes=[1,2,2,2], lw=1.2, color=(0,0,0))
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Kind of want to highlight the regions outside of the confidence
listOfCandidates = []
for x,y in zip(xData, yData):
	if  y > median+std or median-std > y:
		listOfCandidates.append(x)
##########################################################################
########### Printing stuff ###############################################
with open(details_file, 'w') as file:
	header = '''# Summary file containing\n# locations that exceeded GC limits\n# High:{0} \n# Low:{1}\n'''.format(round( median+std, 2 ), round( median-std, 2))
	file.write(header)
	for cand in listOfCandidates:
		file.write(f'{cand}\n')
plt.savefig(output_file, dpi=600)
#plt.show()
"""
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

class Genome2Tetramer():
	'''
	Purpose:
		Given a genome, parse it, and transform it into a pandas df
		where columns are all possible tetramers.

	Usage:
		transformer = Genome2Tetramer()
		transformer.make_csv(inputfile, outputfile)
		'''
	import pandas as pd
	def __init__(self, inputGenome, outputCSV='tetraFrequency.csv', fragsize = 1000):
		print('hello, initialized')
		self.tetramerDicRepository = {key:[] for key in ['AAAA', 'AAAC', 'AAAG', 'AAAT', 'AACA', 'AACC', 'AACG', 'AACT', 'AAGA', 'AAGC', 'AAGG', 'AAGT', 'AATA', 'AATC', 'AATG', 'AATT', 'ACAA', 'ACAC', 'ACAG', 'ACAT', 'ACCA', 'ACCC', 'ACCG', 'ACCT', 'ACGA', 'ACGC', 'ACGG', 'ACGT', 'ACTA', 'ACTC', 'ACTG', 'ACTT', 'AGAA', 'AGAC', 'AGAG', 'AGAT', 'AGCA', 'AGCC', 'AGCG', 'AGCT', 'AGGA', 'AGGC', 'AGGG', 'AGGT', 'AGTA', 'AGTC', 'AGTG', 'AGTT', 'ATAA', 'ATAC', 'ATAG', 'ATAT', 'ATCA', 'ATCC', 'ATCG', 'ATCT', 'ATGA', 'ATGC', 'ATGG', 'ATGT', 'ATTA', 'ATTC', 'ATTG', 'ATTT', 'CAAA', 'CAAC', 'CAAG', 'CAAT', 'CACA', 'CACC', 'CACG', 'CACT', 'CAGA', 'CAGC', 'CAGG', 'CAGT', 'CATA', 'CATC', 'CATG', 'CATT', 'CCAA', 'CCAC', 'CCAG', 'CCAT', 'CCCA', 'CCCC', 'CCCG', 'CCCT', 'CCGA', 'CCGC', 'CCGG', 'CCGT', 'CCTA', 'CCTC', 'CCTG', 'CCTT', 'CGAA', 'CGAC', 'CGAG', 'CGAT', 'CGCA', 'CGCC', 'CGCG', 'CGCT', 'CGGA', 'CGGC', 'CGGG', 'CGGT', 'CGTA', 'CGTC', 'CGTG', 'CGTT', 'CTAA', 'CTAC', 'CTAG', 'CTAT', 'CTCA', 'CTCC', 'CTCG', 'CTCT', 'CTGA', 'CTGC', 'CTGG', 'CTGT', 'CTTA', 'CTTC', 'CTTG', 'CTTT', 'GAAA', 'GAAC', 'GAAG', 'GAAT', 'GACA', 'GACC', 'GACG', 'GACT', 'GAGA', 'GAGC', 'GAGG', 'GAGT', 'GATA', 'GATC', 'GATG', 'GATT', 'GCAA', 'GCAC', 'GCAG', 'GCAT', 'GCCA', 'GCCC', 'GCCG', 'GCCT', 'GCGA', 'GCGC', 'GCGG', 'GCGT', 'GCTA', 'GCTC', 'GCTG', 'GCTT', 'GGAA', 'GGAC', 'GGAG', 'GGAT', 'GGCA', 'GGCC', 'GGCG', 'GGCT', 'GGGA', 'GGGC', 'GGGG', 'GGGT', 'GGTA', 'GGTC', 'GGTG', 'GGTT', 'GTAA', 'GTAC', 'GTAG', 'GTAT', 'GTCA', 'GTCC', 'GTCG', 'GTCT', 'GTGA', 'GTGC', 'GTGG', 'GTGT', 'GTTA', 'GTTC', 'GTTG', 'GTTT', 'TAAA', 'TAAC', 'TAAG', 'TAAT', 'TACA', 'TACC', 'TACG', 'TACT', 'TAGA', 'TAGC', 'TAGG', 'TAGT', 'TATA', 'TATC', 'TATG', 'TATT', 'TCAA', 'TCAC', 'TCAG', 'TCAT', 'TCCA', 'TCCC', 'TCCG', 'TCCT', 'TCGA', 'TCGC', 'TCGG', 'TCGT', 'TCTA', 'TCTC', 'TCTG', 'TCTT', 'TGAA', 'TGAC', 'TGAG', 'TGAT', 'TGCA', 'TGCC', 'TGCG', 'TGCT', 'TGGA', 'TGGC', 'TGGG', 'TGGT', 'TGTA', 'TGTC', 'TGTG', 'TGTT', 'TTAA', 'TTAC', 'TTAG', 'TTAT', 'TTCA', 'TTCC', 'TTCG', 'TTCT', 'TTGA', 'TTGC', 'TTGG', 'TTGT', 'TTTA', 'TTTC', 'TTTG', 'TTTT']}
		self.tetramerDic = {key:0 for key in ['AAAA', 'AAAC', 'AAAG', 'AAAT', 'AACA', 'AACC', 'AACG', 'AACT', 'AAGA', 'AAGC', 'AAGG', 'AAGT', 'AATA', 'AATC', 'AATG', 'AATT', 'ACAA', 'ACAC', 'ACAG', 'ACAT', 'ACCA', 'ACCC', 'ACCG', 'ACCT', 'ACGA', 'ACGC', 'ACGG', 'ACGT', 'ACTA', 'ACTC', 'ACTG', 'ACTT', 'AGAA', 'AGAC', 'AGAG', 'AGAT', 'AGCA', 'AGCC', 'AGCG', 'AGCT', 'AGGA', 'AGGC', 'AGGG', 'AGGT', 'AGTA', 'AGTC', 'AGTG', 'AGTT', 'ATAA', 'ATAC', 'ATAG', 'ATAT', 'ATCA', 'ATCC', 'ATCG', 'ATCT', 'ATGA', 'ATGC', 'ATGG', 'ATGT', 'ATTA', 'ATTC', 'ATTG', 'ATTT', 'CAAA', 'CAAC', 'CAAG', 'CAAT', 'CACA', 'CACC', 'CACG', 'CACT', 'CAGA', 'CAGC', 'CAGG', 'CAGT', 'CATA', 'CATC', 'CATG', 'CATT', 'CCAA', 'CCAC', 'CCAG', 'CCAT', 'CCCA', 'CCCC', 'CCCG', 'CCCT', 'CCGA', 'CCGC', 'CCGG', 'CCGT', 'CCTA', 'CCTC', 'CCTG', 'CCTT', 'CGAA', 'CGAC', 'CGAG', 'CGAT', 'CGCA', 'CGCC', 'CGCG', 'CGCT', 'CGGA', 'CGGC', 'CGGG', 'CGGT', 'CGTA', 'CGTC', 'CGTG', 'CGTT', 'CTAA', 'CTAC', 'CTAG', 'CTAT', 'CTCA', 'CTCC', 'CTCG', 'CTCT', 'CTGA', 'CTGC', 'CTGG', 'CTGT', 'CTTA', 'CTTC', 'CTTG', 'CTTT', 'GAAA', 'GAAC', 'GAAG', 'GAAT', 'GACA', 'GACC', 'GACG', 'GACT', 'GAGA', 'GAGC', 'GAGG', 'GAGT', 'GATA', 'GATC', 'GATG', 'GATT', 'GCAA', 'GCAC', 'GCAG', 'GCAT', 'GCCA', 'GCCC', 'GCCG', 'GCCT', 'GCGA', 'GCGC', 'GCGG', 'GCGT', 'GCTA', 'GCTC', 'GCTG', 'GCTT', 'GGAA', 'GGAC', 'GGAG', 'GGAT', 'GGCA', 'GGCC', 'GGCG', 'GGCT', 'GGGA', 'GGGC', 'GGGG', 'GGGT', 'GGTA', 'GGTC', 'GGTG', 'GGTT', 'GTAA', 'GTAC', 'GTAG', 'GTAT', 'GTCA', 'GTCC', 'GTCG', 'GTCT', 'GTGA', 'GTGC', 'GTGG', 'GTGT', 'GTTA', 'GTTC', 'GTTG', 'GTTT', 'TAAA', 'TAAC', 'TAAG', 'TAAT', 'TACA', 'TACC', 'TACG', 'TACT', 'TAGA', 'TAGC', 'TAGG', 'TAGT', 'TATA', 'TATC', 'TATG', 'TATT', 'TCAA', 'TCAC', 'TCAG', 'TCAT', 'TCCA', 'TCCC', 'TCCG', 'TCCT', 'TCGA', 'TCGC', 'TCGG', 'TCGT', 'TCTA', 'TCTC', 'TCTG', 'TCTT', 'TGAA', 'TGAC', 'TGAG', 'TGAT', 'TGCA', 'TGCC', 'TGCG', 'TGCT', 'TGGA', 'TGGC', 'TGGG', 'TGGT', 'TGTA', 'TGTC', 'TGTG', 'TGTT', 'TTAA', 'TTAC', 'TTAG', 'TTAT', 'TTCA', 'TTCC', 'TTCG', 'TTCT', 'TTGA', 'TTGC', 'TTGG', 'TTGT', 'TTTA', 'TTTC', 'TTTG', 'TTTT']}

		self.fragsize = fragsize
		self.inputGenome = inputGenome
		self.outputCSV = outputCSV

		self.df = None # in order to store the created df

	def make_csv(self, write2file = True, freq=True):
		'''
		Main function to use. Given the input genome (fasta format),
		it will transform it into a tetranucleotide pandas file,
		where every line is frament and every column is a tetranucleotide.
		Make sure the input has only a single entry, it will not read past
		the first.
		'''
		genomeFragmentList = self._parse_genome()
		tetramerDictionary = self._fill_dict(genomeFragmentList)
		df = pd.DataFrame(data=tetramerDictionary)
		self.df = df
		if not freq:
			if write2file:
				df.to_csv(self.outputCSV, index=False)
			return df
		else:
			sums = df.sum()
			df_freq = df.loc[:,'AAAA':'TTTT'].div(sums, axis=1)
			if write2file:
				df_freq.to_csv(self.outputCSV, index=False)
			return df_freq
	def total_tetra(self):
		'''Returns a series of the overall genome freqs'''
		colSums = self.df.sum()
		total = colSums.sum()
		return colSums/total
	def give_GC(self):
		genomeGCList = []
		genomeFragmentList = self._parse_genome()
		for frag in genomeFragmentList:
			GCcount = 0
			for letter in frag:
				if letter == 'G' or letter == 'C':
					GCcount += 1
			genomeGCList.append(GCcount/self.fragsize)

		return pd.Series(genomeGCList)

	def _parse_genome(self):
		'''
		Parses the input genome into a list of fragments. Can only
		handle a single fasta sequence (1 header)
		'''
		lenFragment = self.fragsize
		fileName = self.inputGenome
		genomeFragmentList = []
		with open(fileName) as file:
			fragment = ""
			print(file.readline())
			for line in file.readlines():

				line = line.rstrip().upper()

				if len(fragment) < lenFragment:
					pass
				else:
					genomeFragmentList.append(fragment[:lenFragment])
					fragment = fragment[lenFragment:]

				try:
					if line[0] == '>':
						print(line)
						break
					else:
						fragment += line
				except IndexError:
					print(line)

		print(len(genomeFragmentList))
		return genomeFragmentList
	def _fill_dict(self, genomeFragmentList):
		'''
		Given the list of genome fragments, this transfoms it into
		a dictionary of counts where the keys are tetranucleotides
		'''
		begin = time.time()
		tetramerDic = self.tetramerDic
		tetramerDicRepository = self.tetramerDicRepository

		for n, fragment in enumerate(genomeFragmentList):
		    tetramerDic = {key:0 for key in ['AAAA', 'AAAC', 'AAAG', 'AAAT', 'AACA', 'AACC', 'AACG', 'AACT', 'AAGA', 'AAGC', 'AAGG', 'AAGT', 'AATA', 'AATC', 'AATG', 'AATT', 'ACAA', 'ACAC', 'ACAG', 'ACAT', 'ACCA', 'ACCC', 'ACCG', 'ACCT', 'ACGA', 'ACGC', 'ACGG', 'ACGT', 'ACTA', 'ACTC', 'ACTG', 'ACTT', 'AGAA', 'AGAC', 'AGAG', 'AGAT', 'AGCA', 'AGCC', 'AGCG', 'AGCT', 'AGGA', 'AGGC', 'AGGG', 'AGGT', 'AGTA', 'AGTC', 'AGTG', 'AGTT', 'ATAA', 'ATAC', 'ATAG', 'ATAT', 'ATCA', 'ATCC', 'ATCG', 'ATCT', 'ATGA', 'ATGC', 'ATGG', 'ATGT', 'ATTA', 'ATTC', 'ATTG', 'ATTT', 'CAAA', 'CAAC', 'CAAG', 'CAAT', 'CACA', 'CACC', 'CACG', 'CACT', 'CAGA', 'CAGC', 'CAGG', 'CAGT', 'CATA', 'CATC', 'CATG', 'CATT', 'CCAA', 'CCAC', 'CCAG', 'CCAT', 'CCCA', 'CCCC', 'CCCG', 'CCCT', 'CCGA', 'CCGC', 'CCGG', 'CCGT', 'CCTA', 'CCTC', 'CCTG', 'CCTT', 'CGAA', 'CGAC', 'CGAG', 'CGAT', 'CGCA', 'CGCC', 'CGCG', 'CGCT', 'CGGA', 'CGGC', 'CGGG', 'CGGT', 'CGTA', 'CGTC', 'CGTG', 'CGTT', 'CTAA', 'CTAC', 'CTAG', 'CTAT', 'CTCA', 'CTCC', 'CTCG', 'CTCT', 'CTGA', 'CTGC', 'CTGG', 'CTGT', 'CTTA', 'CTTC', 'CTTG', 'CTTT', 'GAAA', 'GAAC', 'GAAG', 'GAAT', 'GACA', 'GACC', 'GACG', 'GACT', 'GAGA', 'GAGC', 'GAGG', 'GAGT', 'GATA', 'GATC', 'GATG', 'GATT', 'GCAA', 'GCAC', 'GCAG', 'GCAT', 'GCCA', 'GCCC', 'GCCG', 'GCCT', 'GCGA', 'GCGC', 'GCGG', 'GCGT', 'GCTA', 'GCTC', 'GCTG', 'GCTT', 'GGAA', 'GGAC', 'GGAG', 'GGAT', 'GGCA', 'GGCC', 'GGCG', 'GGCT', 'GGGA', 'GGGC', 'GGGG', 'GGGT', 'GGTA', 'GGTC', 'GGTG', 'GGTT', 'GTAA', 'GTAC', 'GTAG', 'GTAT', 'GTCA', 'GTCC', 'GTCG', 'GTCT', 'GTGA', 'GTGC', 'GTGG', 'GTGT', 'GTTA', 'GTTC', 'GTTG', 'GTTT', 'TAAA', 'TAAC', 'TAAG', 'TAAT', 'TACA', 'TACC', 'TACG', 'TACT', 'TAGA', 'TAGC', 'TAGG', 'TAGT', 'TATA', 'TATC', 'TATG', 'TATT', 'TCAA', 'TCAC', 'TCAG', 'TCAT', 'TCCA', 'TCCC', 'TCCG', 'TCCT', 'TCGA', 'TCGC', 'TCGG', 'TCGT', 'TCTA', 'TCTC', 'TCTG', 'TCTT', 'TGAA', 'TGAC', 'TGAG', 'TGAT', 'TGCA', 'TGCC', 'TGCG', 'TGCT', 'TGGA', 'TGGC', 'TGGG', 'TGGT', 'TGTA', 'TGTC', 'TGTG', 'TGTT', 'TTAA', 'TTAC', 'TTAG', 'TTAT', 'TTCA', 'TTCC', 'TTCG', 'TTCT', 'TTGA', 'TTGC', 'TTGG', 'TTGT', 'TTTA', 'TTTC', 'TTTG', 'TTTT']}
		    for pos in range(0,len(fragment)-3):
		        tetramer = fragment[pos:pos+4]
		        try:
		            tetramerDic[tetramer] += 1
		        except KeyError:
		            pass

		    for key, value in tetramerDic.items():
		        tetramerDicRepository[key].append(value)

		    if n%500 == 0:
		        end = time.time()
		        print(f'Till {n} run took {end-begin}')
		        begin = time.time()
		print('done')
		return tetramerDicRepository

inputfile = "SaccCere/SaccCere_genomic.fna"#"Hodgkinia_cicadicola_tetund1_genome.fasta"#
outputfile = 'SaccCere/tetraFrequency.csv'
outFig = 'SaccCere/PCA.png'

transformer = Genome2Tetramer(inputfile, outputCSV=outputfile, fragsize=1000)
tetraFreq = transformer.make_csv(write2file=False)
GClist = transformer.give_GC()

####################################################################
####################################################################
################### Defining matplotlib canvas #####################
figureWidth = 10
figureHeight = 12
plt.figure(figsize=(figureWidth, figureHeight))

#>>>>>>>>>>>>> Data Panel
dataPanelHeight = 4 ; relativedataPanelHeight = dataPanelHeight/figureHeight
dataPanelWidth = 5 ; relativedataPanelWidth = dataPanelWidth/figureWidth

dataPanel = plt.axes([0.45, 0.15,relativedataPanelWidth,relativedataPanelHeight])
dataPanel.tick_params(bottom=False, labelbottom=False,
    left=False, labelleft=False,
    right=False, labelright=False,
    top=False, labeltop=False)

#>>>>>>>>>>>>> Var Panel
varPanelHeight = 4 ; relativeVarPanelHeight = varPanelHeight/figureHeight
varPanelWidth = 3.2 ; relativeVarPanelWidth = varPanelWidth/figureWidth

varPanel = plt.axes([0.05, 0.15,relativeVarPanelWidth,relativeVarPanelHeight])

#>>>>>>>>>>>>> Genome Panel
genomePanelHeight = 2 ; relativegenomePanelHeight = genomePanelHeight/figureHeight
genomePanelWidth = 9 ; relativegenomePanelWidth = genomePanelWidth/figureWidth

genomePanel = plt.axes([0.05, 0.6,relativegenomePanelWidth,relativegenomePanelHeight])
genomePanel.tick_params(bottom=True, labelbottom=True,
    left=True, labelleft=True,
    right=False, labelright=False,
    top=False, labeltop=False)
####################################################################
####################################################################
################### End of matplotlib canvas #######################


class TetraPCAPlotter():
	'''
	I want this to be given any df with tetramers
	and be able to PCA and plot it.


	'''
	import pandas as pd
	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib.patches as mplpatches
	from sklearn.decomposition import PCA
	from sklearn.preprocessing import StandardScaler
	def __init__(self, tetramerDF,
		dataPanel, variancePanel, genomePanel,
		fragsize=1000, n_components=6, x_axis_pca=1,
		y_axis_pca=2, outputFigure=None):

		self.tetramerDF = tetramerDF #raw tetramer dataframe
		self.n_components = n_components #how many pca components
		self.fragsize = fragsize #how large the fragments were

		self.dataPanel = dataPanel
		self.variancePanel = variancePanel
		self.genomePanel = genomePanel

		self.x_axis_pca = x_axis_pca
		self.y_axis_pca = y_axis_pca

		self.pca = None
		self.finalDF = None


	def make_plot(self):
		self._plotData()
		self._plotVariance()
		self._plotRainbowGenome()
	def _plotData(self):
		'''Plots the PCA data. Max and mins are a little weird'''
		finalDF, pca = self._curateData(self.tetramerDF, self.n_components)
		self.pca = pca
		colors = pd.Series(self._makeRainbowList(finalDF), name='color')
		finalDF = pd.concat([finalDF, colors], axis=1)

		x_pca_data = finalDF[f'pc{self.x_axis_pca}']
		y_pca_data = finalDF[f'pc{self.y_axis_pca}']

		self.dataPanel.scatter(x_pca_data,y_pca_data,c=finalDF['color'],s=7)

		xmin = min(list(x_pca_data))
		xmax = max(list(x_pca_data))
		self.dataPanel.set_xticks([x for x in np.linspace(xmin,xmax,20)])

		ymin = min(list(y_pca_data))
		ymax = max(list(y_pca_data))
		self.dataPanel.set_yticks([y for y in np.linspace(ymin,ymax,20)])

		self.dataPanel.set_xlabel(f'\nPC{self.x_axis_pca}', fontsize=18)
		self.dataPanel.set_ylabel(f'PC{self.y_axis_pca}\n', fontsize=18)
		self.dataPanel.grid(color='grey', linestyle='-', linewidth=1, alpha=0.1)

		xspan = xmax-xmin ; xbuffer = 0.1*xspan
		yspan = ymax-ymin ; ybuffer = 0.1*yspan
		self.dataPanel.set_xlim(xmin-xbuffer,xmax+xbuffer)
		self.dataPanel.set_ylim(ymin-ybuffer,ymax+ybuffer)
		self.finalDF = finalDF

	def boundries(self, x_boundry, y_boundry, verb_output=False):
		x_pca_data = self.finalDF[f'pc{self.x_axis_pca}']
		y_pca_data = self.finalDF[f'pc{self.y_axis_pca}']

		xmin = min(x_pca_data)
		xmax = max(x_pca_data)
		ymin = min(y_pca_data)
		ymax = max(y_pca_data)

		xspan = xmax-xmin ; xbuffer = 0.1*xspan
		yspan = ymax-ymin ; ybuffer = 0.1*yspan

		self.dataPanel.plot([xmin-xbuffer, xmax+xbuffer],[y_boundry,y_boundry],marker='o',markerfacecolor=(0,0,0),
			markersize=0, markeredgewidth=0, linewidth=1, color=(0,0,0) )

		self.dataPanel.plot([x_boundry,x_boundry], [ymin-ybuffer, ymax+ybuffer],marker='o',markerfacecolor=(0,0,0),
			markersize=0, markeredgewidth=0, linewidth=1, color=(0,0,0) )

		if verb_output:
			#with open('pca_interesting_points.txt', 'w') as file:
			firstPCA = self.fragsize * np.array( x_pca_data[ x_pca_data >= x_boundry ].index, dtype=int )
			secondPCA = self.fragsize * np.array( y_pca_data[ y_pca_data >= y_boundry ].index, dtype=int )
			print( firstPCA, secondPCA )

	def _plotRainbowGenome(self):
		n_chunk = self.finalDF.shape[0]
		for i in range(0,n_chunk):
			chunkSize =self.fragsize
			left = i*chunkSize
			rectangleHeight = 0.3
			rectangleWidth = chunkSize
			color = self.finalDF['color'].iloc[i]
			rectangle=mplpatches.Rectangle([left,0.5],
											rectangleWidth,
											rectangleHeight,
											facecolor=color,
											edgecolor='black',
											linewidth=0)
			self.genomePanel.add_patch(rectangle)

		self.genomePanel.set_ylim([0,1])
		self.genomePanel.set_xlim([0,chunkSize*n_chunk])

		xmin = 0
		xmax = self.fragsize*n_chunk
		self.genomePanel.set_xticks([x for x in np.linspace(xmin,xmax,30)])
		self.genomePanel.set_xticklabels(["{:.2e}".format( round(x,0) ) for x in np.linspace(xmin,xmax,30)], rotation=90)
		self.genomePanel.set_title('Colorful Genome Reference', fontdict={'fontsize':18})

	def plotGCGenome(self, GClist):
		x = np.array( GClist.index, dtype=int ) * self.fragsize
		y = GClist.values
		self.genomePanel.plot(x,y,marker='o',markerfacecolor=(0,0,0),
			markersize=0, markeredgewidth=0, linewidth=1, color=(0,0,0))
	def _plotGenomeFeatures(self):
		pass
	def _plotVariance(self):
		rectangleWidth = 0.7
		centeringShift = 1-rectangleWidth/2
		xTics = []
		xTicLabels = []
		for n, variance in enumerate(self.pca.explained_variance_ratio_):
		    rectangle=mplpatches.Rectangle( [n+centeringShift,0],rectangleWidth,variance,
		        facecolor=( 0.2,0.2,0.2 ),
		        edgecolor='black',
		        linewidth=0)
		    self.variancePanel.add_patch(rectangle)
		    xTics.append(n+1)
		    xTicLabels.append(f'PCA{n+1}')

		self.variancePanel.set_xlim(0.2,len(self.pca.explained_variance_ratio_)+0.8)
		self.variancePanel.set_ylim(0,max(self.pca.explained_variance_ratio_)*1.1)

		self.variancePanel.set_xticks(xTics)
		self.variancePanel.set_xticklabels(xTicLabels)

		self.variancePanel.set_ylabel('Explained Variance\n', fontsize=18)


	def _makeRainbowList(self, finalDf, colors=None):
	    '''Returns a list to make a rainbow'''
	    red = [x/255 for x in [255,0,0]]
	    orange = [x/255 for x in [245,182,66]]
	    yellow = [x/255 for x in [245,236,66]]
	    green = [x/255 for x in [66,245,84]]
	    blue = [x/255 for x in [66,69,245]]
	    violet = [x/255 for x in [245,66,242]]
	    if colors:
	        colors = colors
	    else:
	        colors = [red,orange,yellow,green,blue,violet]
	    transLen = int(finalDf.shape[0]/(len(colors)-1))
	    transMatrix = []

	    for n in range(0,len(colors)-1):
	        color1 = colors[n]
	        color2 = colors[n+1]
	        R = [x for x in np.linspace(color1[0], color2[0], transLen)]
	        G = [x for x in np.linspace(color1[1], color2[1], transLen)]
	        B = [x for x in np.linspace(color1[2], color2[2], transLen)]

	        for x in zip(R,G,B):
	            transMatrix.append(x)
	    while len(transMatrix) < finalDf.shape[0]:
	        transMatrix.append([0.5,0.5,0.5])

	    return transMatrix
	def _curateData(self, tetramerDF, n):
	    '''
	    Given a dataframe of 256 dimensions,
	    returns a df PCA of n dim and a pca for
	    explainedVariance
	    '''

	    features = tetramerDF.columns
	    x = tetramerDF.loc[:, features].values
	    x = StandardScaler().fit_transform(x)
	    pca = PCA(n_components=n)
	    columnNames = []
	    for column in range(0,n):
	        columnNames.append(f'pc{column+1}')
	    principalComponents = pca.fit_transform(x)
	    principalDf = pd.DataFrame(data=principalComponents, columns = columnNames)
	    return principalDf, pca

thing = TetraPCAPlotter(tetraFreq,dataPanel,varPanel,genomePanel,
	fragsize=transformer.fragsize,
	x_axis_pca=1, y_axis_pca=2)

thing.make_plot()
thing.plotGCGenome(GClist=GClist)
thing.boundries(8,8, verb_output=True)

plt.savefig(outFig, dpi=600)
plt.show()

