
import pandas as pd





# read in the samples.csv file
samples = pd.read_csv("samples.csv").set_index("sample")

# from samples, make a data frame of the distinct libraries so
# we can get the infile locations for them
libs = samples[["lib", "infile"]].drop_duplicates().set_index("lib")


################# This is all the stuff I did for dealing with the demulti-plexing problem ########### 

# convenient list of all the big files.
lib_list = list(libs.index)

# this is a function that returns all the extracted files associated
# with a library. It is NOT an input function
def extracted_files_from_lib(L):
	s = list(samples.loc[samples['lib'] == L].index)
	return expand("results/extracted/{S}.fq.gz", S = s)

# then, make a dictionary of those files for each library.
exf_dict = dict()
for L in lib_list:
	exf_dict[L] = extracted_files_from_lib(L)
# make a dictionary of infiles for each library too
inf_dict = dict()
for  L in lib_list:
	inf_dict[L] = libs.loc[L, "infile"]


# now, use python to write a separate rule for each library/big file
# in a new .smk file that we will include. 
rule_template = """

rule extract_{L}:
	input: inf_dict['{L}']
	output: exf_dict['{L}']
	shell:
		"for i in $(cat {{input}}); do echo $i > results/extracted/$i; done;"

"""

with  open('data-determined-demultiplexing-rules.smk','w') as myfile:
	for L in lib_list:
		myfile.write(rule_template.format(L = L))

# here we include that snakefile that has one rule for each
# big fastq that needs demultiplexing
include: "data-determined-demultiplexing-rules.smk"


############# DONE WITH DEALING WITH THE DEMULTIPLEXING PROBLEM #########################



def infile_from_lib(wildcards):
	return libs.loc[wildcards.lib, "infile"]


def flag_from_sample(wildcards):
	return r"results/flags/lib-flag-{lib}".format(lib = samples.loc[wildcards.sample, "lib"])


# rule all here is just requesting the variants of everyone
rule all:
	input:
		"results/vcf/everyone.vcf"







# then you can carry on with rules that request the extracted files as
# input.  When any such rule requests any extracted file, it will trigger
# evaluations of the appropriate instance of the exract_samples rule.
rule map_extracted_reads:
	input:
		fq = "results/extracted/{sample}.fq.gz",
	output:
		"results/mapped/{sample}.bam"
	shell:
		"echo {input.fq} > {output}"


# then, if there were a variant-calling step, it would request all of the
# bams as input
rule variants:
	input: 
		expand("results/mapped/{S}.bam", S = list(samples.index))
	output:
		"results/vcf/everyone.vcf"
	shell:
		"echo {input} > {output}"

