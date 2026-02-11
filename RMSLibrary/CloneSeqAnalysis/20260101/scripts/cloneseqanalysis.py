# A script that provides utilities to process Clone-seq data

# A typical set up
# CloneSeq/
#	References/
#	RawData/
#		FastQC_reports/
#	TrimmedData/
#		FastQC_reports/
#	Alignments/
#	Pileup/
#	Selections/

from rmsp import rmsutils
from rmsp.rmstemplate import *

	
	
def RMSTemplate_reference_indexing_20260101(
	rms,
	conda_env: str = "rms_cloneseq_analysis_20260101",
	ref_file: InputFileType = "test.fa",
	thread: int = 1,
	):
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	f = lambda **kwargs: [f'{kwargs["r"]}.{ext}' for ext in ["sa", "pac", "bwt", "ann", "amb"]]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"bwa index {r}",
		conda_env=conda_env,
		output_func=output_func,
		r=rms.file_from_path(ref_file)
	)

def RMSTemplate_data_preprocessing_20260101(
	rms, 
	conda_env: str = "rms_cloneseq_analysis_20260101",
	project_dir: str = "/home/happyuser/CloneSeq/",
	prefix: str="sample",
	# Inputs/Outputs
	ref_file: InputFileType = "test.fa",
	r1: InputFileType = "/home/happyuser/CloneSeq/RawData/example_read.fq.gz", 
	r2: InputFileType = None, 
	# General options
	thread: int = 1,
):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	RAW_DATA_DIR = f"{PROJECT_DIR}RawData/"
	TRIMMED_DATA_DIR = f"{PROJECT_DIR}TrimmedData/"	
	RAW_DATA_FASTQC_DIR = f"{RAW_DATA_DIR}FastQC_reports/"
	TRIMMED_DATA_FASTQC_DIR = f"{TRIMMED_DATA_DIR}FastQC_reports/"
	ALIGNMENT_DIR = f"{PROJECT_DIR}Alignments/"
	PILEUP_DIR = f"{PROJECT_DIR}Pileup/"
	for dname in [PROJECT_DIR, RAW_DATA_DIR, TRIMMED_DATA_DIR, RAW_DATA_FASTQC_DIR, TRIMMED_DATA_FASTQC_DIR, ALIGNMENT_DIR, PILEUP_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	

	alignment_file = f"{ALIGNMENT_DIR}{prefix}.bam"
	vcf_file = f"{PILEUP_DIR}{prefix}.vcf"

	ref_indice = {f"r{ext}":rms.file_from_path(f"{ref_file}.{ext}") for ext in ["sa", "pac", "bwt", "ann", "amb"]}
	
	if r2 is not None:
		# Data trimming and quality check
		t1 = f"{TRIMMED_DATA_DIR}{prefix}_1.fq.gz"
		t2 = f"{TRIMMED_DATA_DIR}{prefix}_2.fq.gz"
		ut1 = f"{TRIMMED_DATA_DIR}{prefix}_unpaired_1.fq.gz"
		ut2 = f"{TRIMMED_DATA_DIR}{prefix}_unpaired_2.fq.gz"
		f = lambda **kwargs: [kwargs[k] for k in ["o1", "o2", "u1", "u2"]]
		output_func = rms.register_pipe(f)
		run_template_bash(
			"trimmomatic PE -threads {thread} {i1} {i2} {o1} {u1} {o2} {u2} {parameter}",
			conda_env=conda_env,
			output_func=output_func,
			i1=rms.file_from_path(r1),
			i2=rms.file_from_path(r2),
			o1=t1,
			u1=ut1,
			o2=t2,
			u2=ut2,
			parameter="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30",
			thread=thread
		)
		
		def f(**kwargs):
			import os
			FASTQC_OUTDIR = os.path.abspath(kwargs['o'])
			output_files = []
			for i in [kwargs['i1'], kwargs['i2']]:
				n = os.path.basename(i)
				if n.endswith(".fastq.gz"):
					fn = n[:-9]
				elif n.endswith(".fq.gz"):
					fn = n[:-6]
				elif n.endswith(".fastq"):
					fn = n[:-6]
				elif n.endswith(".fq"):
					fn = n[:-3]
				else:
					fn = n
				output_files.append(FASTQC_OUTDIR + "/" + fn + "_fastqc.zip")
				output_files.append(FASTQC_OUTDIR + "/" + fn + "_fastqc.html")	
			return output_files
		output_func = rms.register_pipe(f)
		run_template_bash(
			"fastqc -q -t {thread} -o {o} {i1} {i2}",
			conda_env=conda_env,
			output_func=output_func,
			o=RAW_DATA_FASTQC_DIR, 
			i1=rms.file_from_path(r1),
			i2=rms.file_from_path(r2),
			thread=1, 
		)
		run_template_bash(
			"fastqc -q -t {thread} -o {o} {i1} {i2}",
			conda_env=conda_env,
			output_func=output_func,
			o=TRIMMED_DATA_FASTQC_DIR, 
			i1=rms.file_from_path(t1),
			i2=rms.file_from_path(t2),
			thread=1, 
		)
		# Alignment and pileup
		f = lambda **kwargs: [kwargs["o"], kwargs["o"] + ".bai"]
		output_func = rms.register_pipe(f)
		run_template_bash(
			"bwa mem -t {thread} {r} {i1} {i2} | samtools sort -@ {thread} -o {o} -",
			conda_env=conda_env,
			output_func=output_func,
			r=rms.file_from_path(ref_file),
			**ref_indice,
			i1=rms.file_from_path(t1),
			i2=rms.file_from_path(t2),
			o=alignment_file,
			thread=thread
		)
		
	else:
		t1 = f"{TRIMMED_DATA_DIR}{prefix}_filtered.fq.gz"
		f = lambda **kwargs: [kwargs[k] for k in ["o"] if k in kwargs]
		output_func = rms.register_pipe(f)
		run_template_bash(
			"trimmomatic SE -threads {thread} {i} {o} {parameter}",
			conda_env=conda_env,
			output_func=output_func,
			i=rms.file_from_path(r1),
			o=t1,
			parameter="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30",
			thread=thread
		)
		def f(**kwargs):
			import os
			FASTQC_OUTDIR = os.path.abspath(kwargs['o'])
			output_files = []
			i = kwargs['i']
			n = os.path.basename(i)
			if n.endswith(".fastq.gz"):
				fn = n[:-9]
			elif n.endswith(".fq.gz"):
				fn = n[:-6]
			elif n.endswith(".fastq"):
				fn = n[:-6]
			elif n.endswith(".fq"):
				fn = n[:-3]
			else:
				fn = n
			output_files.append(FASTQC_OUTDIR + "/" + fn + "_fastqc.zip")
			output_files.append(FASTQC_OUTDIR + "/" + fn + "_fastqc.html")	
			return output_files
		output_func = rms.register_pipe(f)
		run_template_bash(
			"fastqc -q -t {thread} -o {o} {i}",
			conda_env=conda_env,
			output_func=output_func,
			o=RAW_DATA_FASTQC_DIR, 
			i=rms.file_from_path(r1),
			thread=1, 
		)
		run_template_bash(
			"fastqc -q -t {thread} -o {o} {i}",
			conda_env=conda_env,
			output_func=output_func,
			o=TRIMMED_DATA_FASTQC_DIR, 
			i=rms.file_from_path(t1),
			thread=1, 
		)

		f = lambda **kwargs: [kwargs["o"], kwargs["o"] + ".bai"]
		output_func = rms.register_pipe(f)
		run_template_bash(
			"bwa mem -t {thread} {r} {i} | samtools sort -@ {thread} -o {o} -",
			conda_env=conda_env,
			output_func=output_func,
			r=rms.file_from_path(ref_file),
			**ref_indice,
			i=rms.file_from_path(t1),
			o=alignment_file,
			thread=thread
		)
		
	f = lambda **kwargs: [kwargs["o"]]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"bcftools mpileup --threads {thread} -B -a INFO/AD -d 100000 -L 100000 -f {r} {i} > {o}",
		conda_env=conda_env,
		output_func=output_func,
		r=rms.file_from_path(ref_file),
		i=rms.file_from_path(alignment_file),
		o=vcf_file,
		thread=thread
	)


def RMSTemplate_pool_selection_20260101(
	rms, 
	conda_env: str = "rms_cloneseq_analysis_20260101",
	project_dir: str = "/home/happyuser/CloneSeq/",
	
	# Inputs/Outputs
	ref_file: InputFileType = "test.fa",
	prefice: list[str]=[], 
	output_prefix: str="sample",
):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	PILEUP_DIR = f"{PROJECT_DIR}Pileup/"
	SELECTION_DIR = f"{PROJECT_DIR}Selections/"
	for dname in [PROJECT_DIR, PILEUP_DIR, SELECTION_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})

	selection_file = f"{SELECTION_DIR}{output_prefix}_match.txt"
	detailed_selection_file = f"{SELECTION_DIR}{output_prefix}_detailed_match.txt"
	
	f = lambda **kwargs: [kwargs[k] for k in ["o", "od"]]
	output_func = rms.register_pipe(f)
	
	ik = ["{i" + str(i) + "}" for i in range(len(prefice))]
	i_str = " ".join(ik)
	i_dict = {"i" + str(i): rms.file_from_path(f"{PILEUP_DIR}{p}.vcf") for i, p in enumerate(prefice)}
	pk = ["{p" + str(i) + "}" for i in range(len(prefice))]
	p_str = " ".join(pk)
	p_dict = {"p" + str(i): p for i, p in enumerate(prefice)}
	run_template_bash(
		"biodatatools generate_CloneSeq_pool_selection -r {r} -i " + i_str + " -pool_names " + p_str + " -o {o} -od {od}",
		conda_env=conda_env,
		output_func=output_func,
		r=rms.file_from_path(ref_file),
		**i_dict,
		**p_dict,
		o=selection_file,
		od=detailed_selection_file
	)

	
	
	