#This script is intended for use specifically with a degraded DNA protocol developed at the Centre for Biodiversity Genomics (Univeristy of Guelph).
#This script takes a PacBio CCS FASTA file and performs ths following:
#1) Extracts the file.
#2) Demultiplexes based on 96 forward UMIs and 96 reverse UMIs.
#3) Renames the demultiplexed FASTA files based on a list of sample IDs provided by the user.
#4) Identifies each read to order, then filters out reads that do not match expected order. This can be disabled when starting the script but we recommend
#	leaving it enabled to filter out background noise before de novo assembly.
#5) Reads are run through a custom R script for de novo assembly into full-length barcode sequences.

#User must provide the following before running this script:
#1) Compressed (zip or gz) FASTA file resulting from SMRTLink CCS analysis
#2) A file called samples.txt - this text file is a tab separated file with the following headers:
	#	Well - This can be anything, it's only used to make the file human-readable. We use plate records (A01-H12) to generate this file so well locators make it easy to visualize.
	#	ProcessID - This is the sample ID (cannot include spaces).
	#	Order - Ordinal-level taxonomy of the sample. Should be capitalized.
	#	Identification - Lowest taxonomy ID (if genus-species, use underscore between genus and species name). Can be "none" if no taxonomy available. These IDs are apprended to seqeunce headers.

#NOTE: The file samples.txt cannot contain spaces. Replace with underscores if necessary.
#NOTE: 'Samples.txt' must have sample IDs in same order their corresponding UMIs in '5-PacBio_Forward_1-96.txt' and '3-PacBio_Reverse_1-96.txt', as that order is hard coded below.
#NOTE: The script will check to make sure only one compressed FASTA file and one samples.txt exists before beginning.
#NOTE: All file paths, primer files, and reference library files are specific to our servers. Change paths to your own locations as necessary.
#NOTE: Curently can only handle 96 samples at a time.
#NOTE: The following third party programs must be installed on the server prior to use:
	#	1) Cutadapt
	#	2) fastx_barcode_splitter.pl
	#	3) fastx_reverse_complement
	#	4) assign_taxonomy.py (modified version)

#To start script, place compressed FASTA file and samples.txt into a working directory and type "bash NGSFT.bash" (include full file path of not in $PATH).

#!/bin/bash
clear
echo -e "**********NGSFT DE NOVO ASSEMBLER**********\n\n"
echo -e "You are about to start the de novo assembly of NGSFT Sequel data.
Please make sure the current directory contains only a single ".gz" or ".zip" file 
and "samples.txt".\n"
read -p 'Enter plate name: ' jobid
echo -e "\n"
read -p 'Type DISABLE the turn off the taxonomy filter, otherwise press any key: ' disable
echo -e "\n"
if [ "$disable" != "DISABLE" ]
then
	echo -e "Assembly will proceed with the taxonomy filter ENABLED."
else
	echo -e "You have DISABLED the taxonomy filter."
fi

read -p 'Type "go" to proceed or anything else to exit: ' proceed

if [ "$proceed" = "go" ] || [ "$proceed" = "Go" ] || [ "$proceed" = "GO" ]
then
	echo -e "\nANALYSIS STARTED\n"
else
	clear
	echo -e "\nAborted by user. Exiting script.\n\n"
	exit
fi

#Check current directory for a single ".gz" or ".zip" file and a single file named "samples.txt". Return error if not found.
gzcount=$(ls *.gz | wc -l)
zipcount=$(ls *.zip | wc -l)
txtcount=$(ls *.txt | wc -l)
txtpresent=$(ls *.txt)
fastacount=$(ls *.fasta | wc -l)

if [ "$gzcount" -lt 1 ] && [ "$zipcount" -lt 1 ] && [ "$fastacount" -lt 1 ]
then
	echo -e "\nERROR: Missing sequence data. Please correct and try again. Exiting script.\n\n"
	exit
fi

if [ "$gzcount" -gt 1 ]
then
	echo -e "\nERROR: More than one .gz file detected, please correct and try again. Exiting script. \n\n"
	exit
fi

if [ "$zipcount" -gt 1 ]
then
	echo -e "\nERROR: More than one .zip file detected, please correct and try again. Exiting script. \n\n"
	exit
fi

if [ "$fastacount" -gt 1 ]
then
	echo -e "\nERROR: More than one .fasta file detected, please correct and try again. Exiting script. \n\n"
	exit
fi

if [ "$gzcount" -eq 1 ] && [ "$zipcount" -eq 1 ] && [ "$fastacount" -eq 1 ]
then
	echo -e "\nERROR: More than one sequence file detected, please correct and try again. Exiting script. \n\n"
	exit
fi

if [ "$txtcount" -lt 1 ]
then
	echo -e "\nERROR: Missing samples.txt file. Please correct and try again. Exiting script.\n\n"
	exit
fi

if [ "$txtcount" -gt 1 ]
then
	echo -e "\nERROR: More than one text file detected. Please correct and try again. Exiting script.\n\n"
	exit
fi

if [ "$txtpresent" != "samples.txt" ]
then
	echo -e "\nERROR: Text file should be called 'samples.txt'. Please correct and try again. Exiting script.\n\n"
	exit
fi

if [ "$gzcount" -eq 1 ]
then
	fastqfile=$(ls *.gz)
	tar -xvzf $fastqfile
fi

if [ "$zipcount" -eq 1 ]
then
	fastqfile=$(ls *.zip)
	unzip -j $fastqfile 
	rm -r -f $fastqfile
fi

if [ "$fastacount" -eq 1 ]
then
	fastafile=$(ls *.fasta)
	rename $fastafile "ccs.fasta" $fastafile
	rm -r -f $fastafile
fi

listname=$(ls samples.txt)
#gunzip --stdout $file > ccs.fasta #turn this line on if using gz compressed FASTA

#This removes the pad sequence from each read
cutadapt -g ^TAG -o tmp1.fa ccs.fasta
cutadapt -g ^GTAG -o tmp2.fa tmp1.fa
cutadapt -g ^GGTAG -o tmp3.fa tmp2.fa
cutadapt -a CTA$ -o tmp4.fa tmp3.fa
cutadapt -a CTAC$ -o tmp5.fa tmp4.fa
cutadapt -a CTACC$ -m 220 -o tmp6.fa tmp5.fa

#This demultiplexes reads based on the forward UMI at 5' end, then forward at 3', then reverse at 3', then reverse at 5'
cat tmp6.fa | fastx_barcode_splitter.pl -bcfile /mnt/Files_1/sprosser/Primers/5-PacBio_Forward_1-96.txt --bol --mismatches 0 --prefix fwd5_ --suffix ".fa"
cat fwd5_unmatched.fa | fastx_barcode_splitter.pl -bcfile /mnt/Files_1/sprosser/Primers/3-PacBio_Forward_1-96.txt --eol --mismatches 0 --prefix fwd3_ --suffix ".faq"
cat fwd3_unmatched.faq | fastx_barcode_splitter.pl -bcfile /mnt/Files_1/sprosser/Primers/3-PacBio_Reverse_1-96.txt --eol --mismatches 0 --prefix rev3_ --suffix ".fa"
cat rev3_unmatched.fa | fastx_barcode_splitter.pl -bcfile /mnt/Files_1/sprosser/Primers/5-PacBio_Reverse_1-96.txt --bol --mismatches 0 --prefix rev5_ --suffix ".faq"
for f in *.faq; do fastx_reverse_complement -i $f -o $f.fa; done
rename ".faq" "" *.fa
rm *.faq tmp1.fa tmp2.fa tmp3.fa tmp4.fa tmp5.fa tmp6.fa ccs.fasta

#This combines all UMI splits from the same sample into a single file, and names the files based on the order in 'Samples.txt'
x=($(cat $listname))
cat fwd5_0001_AsF_lbc13.fa fwd3_0001_AsF_lbc13.fa rev5_0001_AsR_lbc76_rc.fa rev3_0001_AsR_lbc76_rc.fa > ${x[5]}"_READS.fasta"
cat fwd5_0002_AsF_lbc28.fa fwd3_0002_AsF_lbc28.fa rev5_0002_AsR_lbc260_rc.fa rev3_0002_AsR_lbc260_rc.fa > ${x[9]}"_READS.fasta"
cat fwd5_0003_AsF_lbc278.fa fwd3_0003_AsF_lbc278.fa rev5_0003_AsR_lbc301_rc.fa rev3_0003_AsR_lbc301_rc.fa > ${x[13]}"_READS.fasta"
cat fwd5_0004_AsF_lbc70.fa fwd3_0004_AsF_lbc70.fa rev5_0004_AsR_lbc157_rc.fa rev3_0004_AsR_lbc157_rc.fa > ${x[17]}"_READS.fasta"
cat fwd5_0005_AsF_lbc291.fa fwd3_0005_AsF_lbc291.fa rev5_0005_AsR_lbc53_rc.fa rev3_0005_AsR_lbc53_rc.fa > ${x[21]}"_READS.fasta"
cat fwd5_0006_AsF_lbc45.fa fwd3_0006_AsF_lbc45.fa rev5_0006_AsR_lbc248_rc.fa rev3_0006_AsR_lbc248_rc.fa > ${x[25]}"_READS.fasta"
cat fwd5_0007_AsF_lbc179.fa fwd3_0007_AsF_lbc179.fa rev5_0007_AsR_lbc255_rc.fa rev3_0007_AsR_lbc255_rc.fa > ${x[29]}"_READS.fasta"
cat fwd5_0008_AsF_lbc163.fa fwd3_0008_AsF_lbc163.fa rev5_0008_AsR_lbc313_rc.fa rev3_0008_AsR_lbc313_rc.fa > ${x[33]}"_READS.fasta"
cat fwd5_0009_AsF_lbc96.fa fwd3_0009_AsF_lbc96.fa rev5_0009_AsR_lbc240_rc.fa rev3_0009_AsR_lbc240_rc.fa > ${x[37]}"_READS.fasta"
cat fwd5_0010_AsF_lbc216.fa fwd3_0010_AsF_lbc216.fa rev5_0010_AsR_lbc138_rc.fa rev3_0010_AsR_lbc138_rc.fa > ${x[41]}"_READS.fasta"
cat fwd5_0011_AsF_lbc50.fa fwd3_0011_AsF_lbc50.fa rev5_0011_AsR_lbc196_rc.fa rev3_0011_AsR_lbc196_rc.fa > ${x[45]}"_READS.fasta"
cat fwd5_0012_AsF_lbc99.fa fwd3_0012_AsF_lbc99.fa rev5_0012_AsR_lbc234_rc.fa rev3_0012_AsR_lbc234_rc.fa > ${x[49]}"_READS.fasta"
cat fwd5_0013_AsF_lbc277.fa fwd3_0013_AsF_lbc277.fa rev5_0013_AsR_lbc21_rc.fa rev3_0013_AsR_lbc21_rc.fa > ${x[53]}"_READS.fasta"
cat fwd5_0014_AsF_lbc344.fa fwd3_0014_AsF_lbc344.fa rev5_0014_AsR_lbc116_rc.fa rev3_0014_AsR_lbc116_rc.fa > ${x[57]}"_READS.fasta"
cat fwd5_0015_AsF_lbc80.fa fwd3_0015_AsF_lbc80.fa rev5_0015_AsR_lbc246_rc.fa rev3_0015_AsR_lbc246_rc.fa > ${x[61]}"_READS.fasta"
cat fwd5_0016_AsF_lbc274.fa fwd3_0016_AsF_lbc274.fa rev5_0016_AsR_lbc322_rc.fa rev3_0016_AsR_lbc322_rc.fa > ${x[65]}"_READS.fasta"
cat fwd5_0017_AsF_lbc200.fa fwd3_0017_AsF_lbc200.fa rev5_0017_AsR_lbc19_rc.fa rev3_0017_AsR_lbc19_rc.fa > ${x[69]}"_READS.fasta"
cat fwd5_0018_AsF_lbc320.fa fwd3_0018_AsF_lbc320.fa rev5_0018_AsR_lbc215_rc.fa rev3_0018_AsR_lbc215_rc.fa > ${x[73]}"_READS.fasta"
cat fwd5_0019_AsF_lbc60.fa fwd3_0019_AsF_lbc60.fa rev5_0019_AsR_lbc187_rc.fa rev3_0019_AsR_lbc187_rc.fa > ${x[77]}"_READS.fasta"
cat fwd5_0020_AsF_lbc373.fa fwd3_0020_AsF_lbc373.fa rev5_0020_AsR_lbc376_rc.fa rev3_0020_AsR_lbc376_rc.fa > ${x[81]}"_READS.fasta"
cat fwd5_0021_AsF_lbc357.fa fwd3_0021_AsF_lbc357.fa rev5_0021_AsR_lbc298_rc.fa rev3_0021_AsR_lbc298_rc.fa > ${x[85]}"_READS.fasta"
cat fwd5_0022_AsF_lbc356.fa fwd3_0022_AsF_lbc356.fa rev5_0022_AsR_lbc51_rc.fa rev3_0022_AsR_lbc51_rc.fa > ${x[89]}"_READS.fasta"
cat fwd5_0023_AsF_lbc374.fa fwd3_0023_AsF_lbc374.fa rev5_0023_AsR_lbc62_rc.fa rev3_0023_AsR_lbc62_rc.fa > ${x[93]}"_READS.fasta"
cat fwd5_0024_AsF_lbc173.fa fwd3_0024_AsF_lbc173.fa rev5_0024_AsR_lbc128_rc.fa rev3_0024_AsR_lbc128_rc.fa > ${x[97]}"_READS.fasta"
cat fwd5_0025_AsF_lbc268.fa fwd3_0025_AsF_lbc268.fa rev5_0025_AsR_lbc247_rc.fa rev3_0025_AsR_lbc247_rc.fa > ${x[101]}"_READS.fasta"
cat fwd5_0026_AsF_lbc227.fa fwd3_0026_AsF_lbc227.fa rev5_0026_AsR_lbc328_rc.fa rev3_0026_AsR_lbc328_rc.fa > ${x[105]}"_READS.fasta"
cat fwd5_0027_AsF_lbc144.fa fwd3_0027_AsF_lbc144.fa rev5_0027_AsR_lbc231_rc.fa rev3_0027_AsR_lbc231_rc.fa > ${x[109]}"_READS.fasta"
cat fwd5_0028_AsF_lbc91.fa fwd3_0028_AsF_lbc91.fa rev5_0028_AsR_lbc210_rc.fa rev3_0028_AsR_lbc210_rc.fa > ${x[113]}"_READS.fasta"
cat fwd5_0029_AsF_lbc178.fa fwd3_0029_AsF_lbc178.fa rev5_0029_AsR_lbc146_rc.fa rev3_0029_AsR_lbc146_rc.fa > ${x[117]}"_READS.fasta"
cat fwd5_0030_AsF_lbc47.fa fwd3_0030_AsF_lbc47.fa rev5_0030_AsR_lbc306_rc.fa rev3_0030_AsR_lbc306_rc.fa > ${x[121]}"_READS.fasta"
cat fwd5_0031_AsF_lbc222.fa fwd3_0031_AsF_lbc222.fa rev5_0031_AsR_lbc93_rc.fa rev3_0031_AsR_lbc93_rc.fa > ${x[125]}"_READS.fasta"
cat fwd5_0032_AsF_lbc34.fa fwd3_0032_AsF_lbc34.fa rev5_0032_AsR_lbc309_rc.fa rev3_0032_AsR_lbc309_rc.fa > ${x[129]}"_READS.fasta"
cat fwd5_0033_AsF_lbc253.fa fwd3_0033_AsF_lbc253.fa rev5_0033_AsR_lbc134_rc.fa rev3_0033_AsR_lbc134_rc.fa > ${x[133]}"_READS.fasta"
cat fwd5_0034_AsF_lbc94.fa fwd3_0034_AsF_lbc94.fa rev5_0034_AsR_lbc84_rc.fa rev3_0034_AsR_lbc84_rc.fa > ${x[137]}"_READS.fasta"
cat fwd5_0035_AsF_lbc300.fa fwd3_0035_AsF_lbc300.fa rev5_0035_AsR_lbc252_rc.fa rev3_0035_AsR_lbc252_rc.fa > ${x[141]}"_READS.fasta"
cat fwd5_0036_AsF_lbc63.fa fwd3_0036_AsF_lbc63.fa rev5_0036_AsR_lbc176_rc.fa rev3_0036_AsR_lbc176_rc.fa > ${x[145]}"_READS.fasta"
cat fwd5_0037_AsF_lbc38.fa fwd3_0037_AsF_lbc38.fa rev5_0037_AsR_lbc382_rc.fa rev3_0037_AsR_lbc382_rc.fa > ${x[149]}"_READS.fasta"
cat fwd5_0038_AsF_lbc120.fa fwd3_0038_AsF_lbc120.fa rev5_0038_AsR_lbc132_rc.fa rev3_0038_AsR_lbc132_rc.fa > ${x[153]}"_READS.fasta"
cat fwd5_0039_AsF_lbc314.fa fwd3_0039_AsF_lbc314.fa rev5_0039_AsR_lbc131_rc.fa rev3_0039_AsR_lbc131_rc.fa > ${x[157]}"_READS.fasta"
cat fwd5_0040_AsF_lbc189.fa fwd3_0040_AsF_lbc189.fa rev5_0040_AsR_lbc161_rc.fa rev3_0040_AsR_lbc161_rc.fa > ${x[161]}"_READS.fasta"
cat fwd5_0041_AsF_lbc168.fa fwd3_0041_AsF_lbc168.fa rev5_0041_AsR_lbc335_rc.fa rev3_0041_AsR_lbc335_rc.fa > ${x[165]}"_READS.fasta"
cat fwd5_0042_AsF_lbc139.fa fwd3_0042_AsF_lbc139.fa rev5_0042_AsR_lbc347_rc.fa rev3_0042_AsR_lbc347_rc.fa > ${x[169]}"_READS.fasta"
cat fwd5_0043_AsF_lbc123.fa fwd3_0043_AsF_lbc123.fa rev5_0043_AsR_lbc22_rc.fa rev3_0043_AsR_lbc22_rc.fa > ${x[173]}"_READS.fasta"
cat fwd5_0044_AsF_lbc302.fa fwd3_0044_AsF_lbc302.fa rev5_0044_AsR_lbc43_rc.fa rev3_0044_AsR_lbc43_rc.fa > ${x[177]}"_READS.fasta"
cat fwd5_0045_AsF_lbc104.fa fwd3_0045_AsF_lbc104.fa rev5_0045_AsR_lbc345_rc.fa rev3_0045_AsR_lbc345_rc.fa > ${x[181]}"_READS.fasta"
cat fwd5_0046_AsF_lbc195.fa fwd3_0046_AsF_lbc195.fa rev5_0046_AsR_lbc95_rc.fa rev3_0046_AsR_lbc95_rc.fa > ${x[185]}"_READS.fasta"
cat fwd5_0047_AsF_lbc156.fa fwd3_0047_AsF_lbc156.fa rev5_0047_AsR_lbc183_rc.fa rev3_0047_AsR_lbc183_rc.fa > ${x[189]}"_READS.fasta"
cat fwd5_0048_AsF_lbc269.fa fwd3_0048_AsF_lbc269.fa rev5_0048_AsR_lbc11_rc.fa rev3_0048_AsR_lbc11_rc.fa > ${x[193]}"_READS.fasta"
cat fwd5_0049_AsF_lbc30.fa fwd3_0049_AsF_lbc30.fa rev5_0049_AsR_lbc66_rc.fa rev3_0049_AsR_lbc66_rc.fa > ${x[197]}"_READS.fasta"
cat fwd5_0050_AsF_lbc221.fa fwd3_0050_AsF_lbc221.fa rev5_0050_AsR_lbc24_rc.fa rev3_0050_AsR_lbc24_rc.fa > ${x[201]}"_READS.fasta"
cat fwd5_0051_AsF_lbc211.fa fwd3_0051_AsF_lbc211.fa rev5_0051_AsR_lbc224_rc.fa rev3_0051_AsR_lbc224_rc.fa > ${x[205]}"_READS.fasta"
cat fwd5_0052_AsF_lbc117.fa fwd3_0052_AsF_lbc117.fa rev5_0052_AsR_lbc136_rc.fa rev3_0052_AsR_lbc136_rc.fa > ${x[209]}"_READS.fasta"
cat fwd5_0053_AsF_lbc106.fa fwd3_0053_AsF_lbc106.fa rev5_0053_AsR_lbc130_rc.fa rev3_0053_AsR_lbc130_rc.fa > ${x[213]}"_READS.fasta"
cat fwd5_0054_AsF_lbc111.fa fwd3_0054_AsF_lbc111.fa rev5_0054_AsR_lbc167_rc.fa rev3_0054_AsR_lbc167_rc.fa > ${x[217]}"_READS.fasta"
cat fwd5_0055_AsF_lbc372.fa fwd3_0055_AsF_lbc372.fa rev5_0055_AsR_lbc2_rc.fa rev3_0055_AsR_lbc2_rc.fa > ${x[221]}"_READS.fasta"
cat fwd5_0056_AsF_lbc330.fa fwd3_0056_AsF_lbc330.fa rev5_0056_AsR_lbc133_rc.fa rev3_0056_AsR_lbc133_rc.fa > ${x[225]}"_READS.fasta"
cat fwd5_0057_AsF_lbc217.fa fwd3_0057_AsF_lbc217.fa rev5_0057_AsR_lbc188_rc.fa rev3_0057_AsR_lbc188_rc.fa > ${x[229]}"_READS.fasta"
cat fwd5_0058_AsF_lbc69.fa fwd3_0058_AsF_lbc69.fa rev5_0058_AsR_lbc232_rc.fa rev3_0058_AsR_lbc232_rc.fa > ${x[233]}"_READS.fasta"
cat fwd5_0059_AsF_lbc245.fa fwd3_0059_AsF_lbc245.fa rev5_0059_AsR_lbc98_rc.fa rev3_0059_AsR_lbc98_rc.fa > ${x[237]}"_READS.fasta"
cat fwd5_0060_AsF_lbc297.fa fwd3_0060_AsF_lbc297.fa rev5_0060_AsR_lbc166_rc.fa rev3_0060_AsR_lbc166_rc.fa > ${x[241]}"_READS.fasta"
cat fwd5_0061_AsF_lbc165.fa fwd3_0061_AsF_lbc165.fa rev5_0061_AsR_lbc290_rc.fa rev3_0061_AsR_lbc290_rc.fa > ${x[245]}"_READS.fasta"
cat fwd5_0062_AsF_lbc57.fa fwd3_0062_AsF_lbc57.fa rev5_0062_AsR_lbc316_rc.fa rev3_0062_AsR_lbc316_rc.fa > ${x[249]}"_READS.fasta"
cat fwd5_0063_AsF_lbc304.fa fwd3_0063_AsF_lbc304.fa rev5_0063_AsR_lbc59_rc.fa rev3_0063_AsR_lbc59_rc.fa > ${x[253]}"_READS.fasta"
cat fwd5_0064_AsF_lbc340.fa fwd3_0064_AsF_lbc340.fa rev5_0064_AsR_lbc279_rc.fa rev3_0064_AsR_lbc279_rc.fa > ${x[257]}"_READS.fasta"
cat fwd5_0065_AsF_lbc68.fa fwd3_0065_AsF_lbc68.fa rev5_0065_AsR_lbc44_rc.fa rev3_0065_AsR_lbc44_rc.fa > ${x[261]}"_READS.fasta"
cat fwd5_0066_AsF_lbc105.fa fwd3_0066_AsF_lbc105.fa rev5_0066_AsR_lbc257_rc.fa rev3_0066_AsR_lbc257_rc.fa > ${x[265]}"_READS.fasta"
cat fwd5_0067_AsF_lbc213.fa fwd3_0067_AsF_lbc213.fa rev5_0067_AsR_lbc118_rc.fa rev3_0067_AsR_lbc118_rc.fa > ${x[269]}"_READS.fasta"
cat fwd5_0068_AsF_lbc141.fa fwd3_0068_AsF_lbc141.fa rev5_0068_AsR_lbc377_rc.fa rev3_0068_AsR_lbc377_rc.fa > ${x[273]}"_READS.fasta"
cat fwd5_0069_AsF_lbc9.fa fwd3_0069_AsF_lbc9.fa rev5_0069_AsR_lbc229_rc.fa rev3_0069_AsR_lbc229_rc.fa > ${x[277]}"_READS.fasta"
cat fwd5_0070_AsF_lbc42.fa fwd3_0070_AsF_lbc42.fa rev5_0070_AsR_lbc326_rc.fa rev3_0070_AsR_lbc326_rc.fa > ${x[281]}"_READS.fasta"
cat fwd5_0071_AsF_lbc64.fa fwd3_0071_AsF_lbc64.fa rev5_0071_AsR_lbc154_rc.fa rev3_0071_AsR_lbc154_rc.fa > ${x[285]}"_READS.fasta"
cat fwd5_0072_AsF_lbc181.fa fwd3_0072_AsF_lbc181.fa rev5_0072_AsR_lbc71_rc.fa rev3_0072_AsR_lbc71_rc.fa > ${x[289]}"_READS.fasta"
cat fwd5_0073_AsF_lbc36.fa fwd3_0073_AsF_lbc36.fa rev5_0073_AsR_lbc74_rc.fa rev3_0073_AsR_lbc74_rc.fa > ${x[293]}"_READS.fasta"
cat fwd5_0074_AsF_lbc204.fa fwd3_0074_AsF_lbc204.fa rev5_0074_AsR_lbc236_rc.fa rev3_0074_AsR_lbc236_rc.fa > ${x[297]}"_READS.fasta"
cat fwd5_0075_AsF_lbc325.fa fwd3_0075_AsF_lbc325.fa rev5_0075_AsR_lbc276_rc.fa rev3_0075_AsR_lbc276_rc.fa > ${x[301]}"_READS.fasta"
cat fwd5_0076_AsF_lbc358.fa fwd3_0076_AsF_lbc358.fa rev5_0076_AsR_lbc56_rc.fa rev3_0076_AsR_lbc56_rc.fa > ${x[305]}"_READS.fasta"
cat fwd5_0077_AsF_lbc238.fa fwd3_0077_AsF_lbc238.fa rev5_0077_AsR_lbc124_rc.fa rev3_0077_AsR_lbc124_rc.fa > ${x[309]}"_READS.fasta"
cat fwd5_0078_AsF_lbc39.fa fwd3_0078_AsF_lbc39.fa rev5_0078_AsR_lbc265_rc.fa rev3_0078_AsR_lbc265_rc.fa > ${x[313]}"_READS.fasta"
cat fwd5_0079_AsF_lbc7.fa fwd3_0079_AsF_lbc7.fa rev5_0079_AsR_lbc379_rc.fa rev3_0079_AsR_lbc379_rc.fa > ${x[317]}"_READS.fasta"
cat fwd5_0080_AsF_lbc140.fa fwd3_0080_AsF_lbc140.fa rev5_0080_AsR_lbc35_rc.fa rev3_0080_AsR_lbc35_rc.fa > ${x[321]}"_READS.fasta"
cat fwd5_0081_AsF_lbc337.fa fwd3_0081_AsF_lbc337.fa rev5_0081_AsR_lbc152_rc.fa rev3_0081_AsR_lbc152_rc.fa > ${x[325]}"_READS.fasta"
cat fwd5_0082_AsF_lbc155.fa fwd3_0082_AsF_lbc155.fa rev5_0082_AsR_lbc126_rc.fa rev3_0082_AsR_lbc126_rc.fa > ${x[329]}"_READS.fasta"
cat fwd5_0083_AsF_lbc184.fa fwd3_0083_AsF_lbc184.fa rev5_0083_AsR_lbc310_rc.fa rev3_0083_AsR_lbc310_rc.fa > ${x[333]}"_READS.fasta"
cat fwd5_0084_AsF_lbc127.fa fwd3_0084_AsF_lbc127.fa rev5_0084_AsR_lbc305_rc.fa rev3_0084_AsR_lbc305_rc.fa > ${x[337]}"_READS.fasta"
cat fwd5_0085_AsF_lbc190.fa fwd3_0085_AsF_lbc190.fa rev5_0085_AsR_lbc261_rc.fa rev3_0085_AsR_lbc261_rc.fa > ${x[341]}"_READS.fasta"
cat fwd5_0086_AsF_lbc170.fa fwd3_0086_AsF_lbc170.fa rev5_0086_AsR_lbc1_rc.fa rev3_0086_AsR_lbc1_rc.fa > ${x[345]}"_READS.fasta"
cat fwd5_0087_AsF_lbc281.fa fwd3_0087_AsF_lbc281.fa rev5_0087_AsR_lbc8_rc.fa rev3_0087_AsR_lbc8_rc.fa > ${x[349]}"_READS.fasta"
cat fwd5_0088_AsF_lbc329.fa fwd3_0088_AsF_lbc329.fa rev5_0088_AsR_lbc266_rc.fa rev3_0088_AsR_lbc266_rc.fa > ${x[353]}"_READS.fasta"
cat fwd5_0089_AsF_lbc153.fa fwd3_0089_AsF_lbc153.fa rev5_0089_AsR_lbc283_rc.fa rev3_0089_AsR_lbc283_rc.fa > ${x[357]}"_READS.fasta"
cat fwd5_0090_AsF_lbc228.fa fwd3_0090_AsF_lbc228.fa rev5_0090_AsR_lbc122_rc.fa rev3_0090_AsR_lbc122_rc.fa > ${x[361]}"_READS.fasta"
cat fwd5_0091_AsF_lbc319.fa fwd3_0091_AsF_lbc319.fa rev5_0091_AsR_lbc354_rc.fa rev3_0091_AsR_lbc354_rc.fa > ${x[365]}"_READS.fasta"
cat fwd5_0092_AsF_lbc3.fa fwd3_0092_AsF_lbc3.fa rev5_0092_AsR_lbc311_rc.fa rev3_0092_AsR_lbc311_rc.fa > ${x[369]}"_READS.fasta"
cat fwd5_0093_AsF_lbc346.fa fwd3_0093_AsF_lbc346.fa rev5_0093_AsR_lbc180_rc.fa rev3_0093_AsR_lbc180_rc.fa > ${x[373]}"_READS.fasta"
cat fwd5_0094_AsF_lbc273.fa fwd3_0094_AsF_lbc273.fa rev5_0094_AsR_lbc151_rc.fa rev3_0094_AsR_lbc151_rc.fa > ${x[377]}"_READS.fasta"
cat fwd5_0095_AsF_lbc370.fa fwd3_0095_AsF_lbc370.fa rev5_0095_AsR_lbc172_rc.fa rev3_0095_AsR_lbc172_rc.fa > ${x[381]}"_READS.fasta"
cat fwd5_0096_AsF_lbc81.fa fwd3_0096_AsF_lbc81.fa rev5_0096_AsR_lbc239_rc.fa rev3_0096_AsR_lbc239_rc.fa > ${x[385]}"_READS.fasta"

rm *.fa

#This runs reads through taxonomy filter (which includes an R script that consolidates BLAST output) and/or de novo assembly R script
#Change file paths below to your own paths 

if [ "$disable" != "DISABLE" ]
then
	#BLAST reads and filter out those that don't match at the order level
	echo -e "Starting taxonomy filter"
	for f in *.fasta
	do
		assign_taxonomy.py -i $f -m blast -r /mnt/Files_1/sprosser/Reference_Libraries/BOLD_UniqueBINs_stripped.fasta -t /mnt/Files_1/sprosser/Reference_Libraries/BOLD_UniqueBINs_tax.txt
	done
	mv /mnt/Files_1/sprosser/ToNGSFT/samples.txt /mnt/Files_1/sprosser/R_Working_Directory
	mv /mnt/Files_1/sprosser/ToNGSFT/blast_assigned_taxonomy/*.txt /mnt/Files_1/sprosser/R_Working_Directory
	mv /mnt/Files_1/sprosser/ToNGSFT/*READS.fasta /mnt/Files_1/sprosser/R_Working_Directory
	echo -e "Starting de novo assembly with filtered reads"
	Rscript --verbose /mnt/Files_1/sprosser/Scripts/BLAST_Results_Analyzer_NGSFT.R
	mkdir /mnt/Files_1/sprosser/ToNGSFT/RawReads
	mkdir /mnt/Files_1/sprosser/ToNGSFT/NonTargetReads
	mv /mnt/Files_1/sprosser/R_Working_Directory/BLAST_Filtering_Summary.txt /mnt/Files_1/sprosser/ToNGSFT/
	mv /mnt/Files_1/sprosser/R_Working_Directory/*assignments.txt /mnt/Files_1/sprosser/ToNGSFT/blast_assigned_taxonomy
	mv /mnt/Files_1/sprosser/R_Working_Directory/*NONTARGETS.fasta /mnt/Files_1/sprosser/ToNGSFT/NonTargetReads
	mv /mnt/Files_1/sprosser/R_Working_Directory/*READS.fasta /mnt/Files_1/sprosser/ToNGSFT/RawReads
	
	#Run reads through de novo assembly script
	Rscript --verbose /mnt/Files_1/sprosser/Scripts/NGSFT.R
	mv /mnt/Files_1/sprosser/R_Working_Directory/samples.txt /mnt/Files_1/sprosser/ToNGSFT/
	mkdir /mnt/Files_1/sprosser/ToNGSFT/DeNovoAssembledReads
	mkdir /mnt/Files_1/sprosser/ToNGSFT/TargetReads
	mv /mnt/Files_1/sprosser/R_Working_Directory/*ASSEMBLED.fasta /mnt/Files_1/sprosser/ToNGSFT/DeNovoAssembledReads
	mv /mnt/Files_1/sprosser/R_Working_Directory/*TARGETS.fasta /mnt/Files_1/sprosser/ToNGSFT/TargetReads
	mv /mnt/Files_1/sprosser/R_Working_Directory/* /mnt/Files_1/sprosser/ToNGSFT/
else
	mv /mnt/Files_1/sprosser/ToNGSFT/*READS.fasta /mnt/Files_1/sprosser/R_Working_Directory
	mv /mnt/Files_1/sprosser/ToNGSFT/samples.txt /mnt/Files_1/sprosser/R_Working_Directory
	
	#Run reads through de novo assembly script
	echo -e "Starting de novo assembly with unfiltered reads"
	Rscript --verbose /mnt/Files_1/sprosser/Scripts/NGSFT.R
	mv /mnt/Files_1/sprosser/R_Working_Directory/samples.txt /mnt/Files_1/sprosser/ToNGSFT/
	mkdir /mnt/Files_1/sprosser/ToNGSFT/RawReads
	mkdir /mnt/Files_1/sprosser/ToNGSFT/DeNovoAssembledReads
	mv /mnt/Files_1/sprosser/R_Working_Directory/*ASSEMBLED.fasta /mnt/Files_1/sprosser/ToNGSFT/DeNovoAssembledReads
	mv /mnt/Files_1/sprosser/R_Working_Directory/*READS.fasta /mnt/Files_1/sprosser/ToNGSFT/RawReads
	mv /mnt/Files_1/sprosser/R_Working_Directory/* /mnt/Files_1/sprosser/ToNGSFT/
fi
rm error.log
mkdir /mnt/Files_1/sprosser/ToNGSFT/IntermediateFiles
mv /mnt/Files_1/sprosser/ToNGSFT/blast_assigned_taxonomy /mnt/Files_1/sprosser/ToNGSFT/IntermediateFiles
mv /mnt/Files_1/sprosser/ToNGSFT/BLAST_Filtering_Summary.txt /mnt/Files_1/sprosser/ToNGSFT/IntermediateFiles
mv '/mnt/Files_1/sprosser/ToNGSFT/Coverage Per Fragment.csv' /mnt/Files_1/sprosser/ToNGSFT/IntermediateFiles
mv '/mnt/Files_1/sprosser/ToNGSFT/Duplex Per Fragment.csv' /mnt/Files_1/sprosser/ToNGSFT/IntermediateFiles
mv /mnt/Files_1/sprosser/ToNGSFT/NonTargetReads /mnt/Files_1/sprosser/ToNGSFT/IntermediateFiles
mv /mnt/Files_1/sprosser/ToNGSFT/TargetReads /mnt/Files_1/sprosser/ToNGSFT/IntermediateFiles
mv /mnt/Files_1/sprosser/ToNGSFT/RawReads /mnt/Files_1/sprosser/ToNGSFT/IntermediateFiles
mv /mnt/Files_1/sprosser/ToNGSFT/*.gz /mnt/Files_1/sprosser/ToNGSFT/IntermediateFiles
mv /mnt/Files_1/sprosser/ToNGSFT/samples.txt /mnt/Files_1/sprosser/ToNGSFT/IntermediateFiles
cp '/mnt/Files_1/sprosser/ToNGSFT/Original Output.fas' /mnt/Files_1/sprosser/ToNGSFT/${jobid}_ToBOLD.fasta
mv '/mnt/Files_1/sprosser/ToNGSFT/Original Output.fas' /mnt/Files_1/sprosser/ToNGSFT/IntermediateFiles
mkdir /mnt/Files_1/sprosser/ToNGSFT/$jobid
mv /mnt/Files_1/sprosser/ToNGSFT/* /mnt/Files_1/sprosser/ToNGSFT/$jobid
zip -rmT $jobid{.zip,}
clear
echo -e "\n\n\nDE NOVO ASSEMBLY COMPLETE\n\n\n"