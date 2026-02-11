# Utilities to perform multiple steps in processing STARR-seq data

# A typical set up
# STARR/
#	References
#	RawData/
#		FastQC_reports/
#		fastp_reports/
#		LayoutInspection
#	TrimmedData/
#		FastQC_reports/
#	ProcessedData/
#	IDMap/
#	ElementActivity/
# 	QC/
#		IDMapQC/
# 		RankedLogFC/
# 		ReplicatesCorrelation/
#		OrientationCorrelation/

from rmsp import rmsutils
from rmsp.rmstemplate import *
from commonhelper import convert_to_bool

def RMSTemplate_plot_STARR_library_layouts_20241022(
	rms, 
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	i:str="InfoTable.csv", 
	output_prefix: str="batch",
	colnames:list[str]=["Read1_FileName", "Read2_FileName"],
	func:str="lambda row: row[\"Read1_FileName\"]"
):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	from biodatatools.utils.common import json_dump, OUTPUT_json_dump
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	json_dump = rms.register_pipe(json_dump, output_func=OUTPUT_json_dump)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	RAW_DATA_DIR = f"{PROJECT_DIR}RawData/"
	LAYOUT_INSPECTION_DIR = f"{RAW_DATA_DIR}LayoutInspection/"
	for dname in [PROJECT_DIR, RAW_DATA_DIR, LAYOUT_INSPECTION_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	json_file = f"{LAYOUT_INSPECTION_DIR}{output_prefix}.json.gz"
	if isinstance(func, str):
		func = eval(func, {})
	if isinstance(i, str):
		raise Exception("This function is not finalized yet!")
	else:
		df = i
	j = rmsutils.build_RMS_json_like_struct(
		rms, 
		{func(row):[rms.file_from_path(f"{RAW_DATA_DIR}{row[colname]}") if len(row[colname]) > 0 else None for colname in colnames] for _, row in df.iterrows()}
	)
	json_dump(json_file, j)
	
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"biodatatools plot_sequence_library_layouts -i {i} -o {o}",
		conda_env=conda_env,
		output_func=output_func,
		i=rms.file_from_path(json_file),
		o=f"{LAYOUT_INSPECTION_DIR}{output_prefix}.png", 
	)
	
	
	

def RMSTemplate_data_trimming_20240701(
	rms, 
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	prefix: str="sample",
	# Inputs/Outputs
	r1: InputFileType = "/home/happyuser/STARR/RawData/example_read1.fq.gz", 
	r2: InputFileType = "/home/happyuser/STARR/RawData/example_read2.fq.gz", 
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
	PROCESSED_READ_DIR = f"{PROJECT_DIR}ProcessedReads/"
	IDMAP_DIR = f"{PROJECT_DIR}IDMap/"
	
	RAW_DATA_FASTQC_DIR = f"{RAW_DATA_DIR}FastQC_reports/"
	TRIMMED_DATA_FASTQC_DIR = f"{TRIMMED_DATA_DIR}FastQC_reports/"
	FASTP_REPORT_DIR = f"{RAW_DATA_DIR}fastp_reports/"
	for dname in [PROJECT_DIR, RAW_DATA_DIR, TRIMMED_DATA_DIR, PROCESSED_READ_DIR, IDMAP_DIR, RAW_DATA_FASTQC_DIR, TRIMMED_DATA_FASTQC_DIR, FASTP_REPORT_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})

	fastp_html_report = f"{FASTP_REPORT_DIR}{prefix}.html"
	fastp_json_report = f"{FASTP_REPORT_DIR}{prefix}.json"

	

	if r2 is not None:
		t1 = f"{TRIMMED_DATA_DIR}{prefix}_1.fq.gz"
		t2 = f"{TRIMMED_DATA_DIR}{prefix}_2.fq.gz"
		f = lambda **kwargs: [kwargs[k] for k in ["o1", "o2", "html_report", "json_report"] if k in kwargs]
		output_func = rms.register_pipe(f)
		run_template_bash(
			"fastp -i {i1} -I {i2} -o {o1} -O {o2} -h {html_report} -j {json_report} --disable_adapter_trimming --trim_poly_g --cut_right --thread {thread}",
			conda_env=conda_env,
			output_func=output_func,
			i1=rms.file_from_path(r1),
			i2=rms.file_from_path(r2),
			o1=t1,
			o2=t2,
			html_report=fastp_html_report,
			json_report=fastp_json_report,
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
	else:
		t1 = f"{TRIMMED_DATA_DIR}{prefix}.fq.gz"
		f = lambda **kwargs: [kwargs[k] for k in ["o", "html_report", "json_report"] if k in kwargs]
		output_func = rms.register_pipe(f)
		run_template_bash(
			"fastp -i {i} -o {o} -h {html_report} -j {json_report} --disable_adapter_trimming --trim_poly_g --cut_right --thread {thread}",
			conda_env=conda_env,
			output_func=output_func,
			i=rms.file_from_path(r1),
			o=t1,
			html_report=fastp_html_report,
			json_report=fastp_json_report,
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


def RMSTemplate_process_fastq_to_delimited_20240701(
	rms, 
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	prefix: str="sample",
	output_prefix:str="sample_output",
	# Inputs/Outputs
	layout1:str=None,
	layout2:str=None,
	column_names:str=None,
	reverse_complements:str=None,
	# General options
	thread: int = 1,
):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	TRIMMED_DATA_DIR = f"{PROJECT_DIR}TrimmedData/"
	PROCESSED_READ_DIR = f"{PROJECT_DIR}ProcessedReads/"
	
	for dname in [PROJECT_DIR, TRIMMED_DATA_DIR, PROCESSED_READ_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})

	pfile = f"{PROCESSED_READ_DIR}{output_prefix}.tsv.gz"

	
	f = lambda **kwargs: [kwargs['o']]
	process_sequencing_reads_to_delimited_sequences_output_func = rms.register_pipe(f)

	params = {}
	if column_names is not None and len(column_names) > 0:
		params["column_names"] = column_names
	if reverse_complements is not None and len(reverse_complements) > 0:
		params["reverse_complements"] = reverse_complements
	param_str = "".join(["-"+k+" {" + k + "} " for k in params.keys()])
	
	if layout2 is not None: # paired
		t1 = f"{TRIMMED_DATA_DIR}{prefix}_1.fq.gz"
		t2 = f"{TRIMMED_DATA_DIR}{prefix}_2.fq.gz"
		run_template_bash(
			"biodatatools process_sequencing_reads_to_delimited_sequences -f1 {i1} -layout1 {l1} -f2 {i2} -layout2 {l2} -o {o} " + param_str + "-nthread {thread}",
			conda_env=conda_env,
			output_func=process_sequencing_reads_to_delimited_sequences_output_func,
			o=pfile, 
			i1=rms.file_from_path(t1),
			l1=layout1,
			i2=rms.file_from_path(t2),
			l2=layout2,
			**params,
			thread=thread, 
		)
	else:
		t1 = f"{TRIMMED_DATA_DIR}{prefix}.fq.gz"
		run_template_bash(
			"biodatatools process_sequencing_reads_to_delimited_sequences -f1 {i1} -layout1 {l1} -o {o} " + param_str + "-nthread {thread}",
			conda_env=conda_env,
			output_func=process_sequencing_reads_to_delimited_sequences_output_func,
			o=pfile, 
			i1=rms.file_from_path(t1),
			l1=layout1,
			**params,
			thread=thread, 
		)

def RMSTemplate_process_delimited_sequences_20240701(
	rms, 
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	prefix: str="sample",
	output_prefix:str="sample_output",
	refs:dict=None,
	process_params:dict={},
	# General options
	thread: int = 1,
):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	PROCESSED_READ_DIR = f"{PROJECT_DIR}ProcessedReads/"
	
	for dname in [PROJECT_DIR, PROCESSED_READ_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})

	ifile = f"{PROCESSED_READ_DIR}{prefix}.tsv.gz"
	ofile = f"{PROCESSED_READ_DIR}{output_prefix}.tsv.gz"

	
	f = lambda **kwargs: [kwargs['o']]
	process_sequencing_reads_to_delimited_sequences_output_func = rms.register_pipe(f)
	def f(**params):
		import csv
		from biodata.baseio import get_text_file_rootname
		def _quote_split(s, delimiter, quotechar='"'):
			return next(csv.reader([s],delimiter=delimiter))
		outputs = [params["o"]]
		if "cluster" in params:
			output_prefix = get_text_file_rootname(params["o"])
			t = params["cluster"]
			for s in _quote_split(t, " ", quotechar="'"):
				if s.startswith("'") and s.endswith("'"):
					s = s[1:-1]
				c = dict([_quote_split(e, "=") for e in _quote_split(s, ";")])
				outputs.append(output_prefix + "_cluster-" + c["column"] + ".tsv.gz")
				outputs.append(output_prefix + "_cluster-" + c["column"] + ".fa.gz")
		return sorted(outputs)
	process_delimited_sequences_output_func = rms.register_pipe(f)
	
	param_str = "".join(["-"+k+" {" + k + "} " for k in process_params.keys()])
	if refs is not None and len(refs) > 0:
		refs = {k:rms.file_from_path(v) for k, v in refs.items()}
		ref_param_str = "-ref " + "'" + ";".join(["{rk" + str(i) + "}={rv" + str(i) + "}" for i, k in enumerate(refs.keys())]) + "'" + " "
		ref_params = {**{"rk" + str(i):k for i, k in enumerate(refs.keys())}, **{"rv" + str(i):v for i, v in enumerate(refs.values())}}
	else:
		ref_param_str = ""
		ref_params = {}

	run_template_bash(
		"biodatatools process_delimited_sequences -i {i} -o {o} " + ref_param_str + param_str + "-nthread {thread}",
		conda_env=conda_env,
		output_func=process_delimited_sequences_output_func,
		i=rms.file_from_path(ifile),
		o=ofile, 
		**process_params,
		**ref_params,
		thread=thread, 
	)

def RMSTemplate_process_multi_delimited_sequences_20240701(
	rms,
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	reseq_prefixes:list[str]=[],
	output_prefix:str="sample_output",
	refs:dict=None,
	process_params:dict={},
	thread:int=1,
):
	reseq_prefixes = sorted(reseq_prefixes)
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	PROJECT_DIR = project_dir
	PROCESSED_READ_DIR = f"{PROJECT_DIR}ProcessedReads/"
	input_str = " ".join(["{i" + str(e) + "}" for e in range(len(reseq_prefixes))])
	reseq_files = {"i" + str(e): rms.file_from_path(f"{PROCESSED_READ_DIR}{s}.tsv.gz") for e, s in enumerate(reseq_prefixes)}
	param_str = "".join(["-"+k+" {" + k + "} " for k in process_params.keys()])
	if refs is not None and len(refs) > 0:
		refs = {k:rms.file_from_path(v) for k, v in refs.items()}
		ref_param_str = "-ref " + "'" + ";".join(["{rk" + str(i) + "}={rv" + str(i) + "}" for i, k in enumerate(refs.keys())]) + "'" + " "
		ref_params = {**{"rk" + str(i):k for i, k in enumerate(refs.keys())}, **{"rv" + str(i):v for i, v in enumerate(refs.values())}}
	else:
		ref_param_str = ""
		ref_params = {}
	def f(**params):
		import csv
		from biodata.baseio import get_text_file_rootname
		def _quote_split(s, delimiter, quotechar='"'):
			return next(csv.reader([s],delimiter=delimiter))
		outputs = [params["o"]]
		if "cluster" in params:
			output_prefix = get_text_file_rootname(params["o"])
			t = params["cluster"]
			for s in _quote_split(t, " ", quotechar="'"):
				if s.startswith("'") and s.endswith("'"):
					s = s[1:-1]
				c = dict([_quote_split(e, "=") for e in _quote_split(s, ";")])
				outputs.append(output_prefix + "_cluster-" + c["column"] + ".tsv.gz")
				outputs.append(output_prefix + "_cluster-" + c["column"] + ".fa.gz")
		return sorted(outputs)
	process_delimited_sequences_output_func = rms.register_pipe(f)
	run_template_bash(
		"biodatatools concat_delimited -i " + input_str + " -o - | biodatatools process_delimited_sequences -i - -o {o} " + ref_param_str + param_str + "-nthread {thread}",
		conda_env=conda_env,
		output_func=process_delimited_sequences_output_func,
		**reseq_files,
		o=f"{PROCESSED_READ_DIR}{output_prefix}.tsv.gz", 
		**process_params,
		**ref_params,
		thread=thread, 
	)
	
def RMSTemplate_cluster_delimited_sequences_20241024(
	rms,
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	prefixes:list[str]=[],
	columns:list[str]=[],
	output_prefix:str="sample_output",
	thread:int=1,
):
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	PROJECT_DIR = project_dir
	PROCESSED_READ_DIR = f"{PROJECT_DIR}ProcessedReads/"
	
	input_str = " ".join(["{i" + str(e) + "}" for e in range(len(prefixes))])
	files = {"i" + str(e): rms.file_from_path(f"{PROCESSED_READ_DIR}{s}.tsv.gz") for e, s in enumerate(prefixes)}
	f = lambda **kwargs: [kwargs['output_prefix'] + ".fa.gz", kwargs['output_prefix'] + ".tsv.gz"]
	cluster_delimited_sequences_output_func = rms.register_pipe(f)
	run_template_bash(
		"biodatatools cluster_delimited_sequences -i " + input_str + " -output_prefix {output_prefix} -column {column} -nthread {thread}",
		conda_env=conda_env,
		output_func=cluster_delimited_sequences_output_func,
		**files,
		output_prefix=f"{PROCESSED_READ_DIR}{output_prefix}", 
		column=" ".join(columns),
		thread=thread, 
	)
def RMSTemplate_generate_IDMap_20240701(
	rms,
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	prefix:str="sample",
	id_column:str="ID",
	target_column:str="",
	count_column:str="Count",
	params={}
):
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	PROJECT_DIR = project_dir
	PROCESSED_READ_DIR = f"{PROJECT_DIR}ProcessedReads/"

	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)

	# id_column_str = " ".join(["{ic" + str(i) + "}" for i, v in enumerate(id_column)])
	# id_column_dict = {"ic" + str(i):v for i, v in enumerate(id_column)}
	# target_column_str = " ".join(["{tc" + str(i) + "}" for i, v in enumerate(target_column)])
	# target_column_dict = {"tc" + str(i):v for i, v in enumerate(target_column)}
	name = "-".join(id_column.split(" ") + ([] if target_column is None or len(target_column) == 0 else target_column.split(" ")))
	param_str = "".join([" -"+k+" {" + k + "}" for k in params.keys()])
	if target_column is None or len(target_column) == 0:
		run_template_bash(
			"biodatatools generate_IDMap_from_delimited -i {i} -o {o} -id_column {id_column} -count_column {count_column}" + param_str,
			conda_env=conda_env,
			output_func=output_func,
			i=rms.file_from_path(f"{PROCESSED_READ_DIR}{prefix}.tsv.gz"),
			o=f"{PROCESSED_READ_DIR}{prefix}_{name}_map.tsv.gz", 
			id_column=id_column, 
			count_column=count_column,
			**params
		)
	else:
		run_template_bash(
			"biodatatools generate_IDMap_from_delimited -i {i} -o {o} -id_column {id_column} -target_column {target_column} -count_column {count_column}" + param_str,
			conda_env=conda_env,
			output_func=output_func,
			i=rms.file_from_path(f"{PROCESSED_READ_DIR}{prefix}.tsv.gz"),
			o=f"{PROCESSED_READ_DIR}{prefix}_{name}_map.tsv.gz", 
			id_column=id_column, 
			target_column=target_column,
			count_column=count_column,
			**params
		)
def RMSTemplate_join_IDMap_20241024(
	rms,
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	prefixes:list[str]=["sample1", "sample2"],
	output_prefix:str="output",
	id_column:str="ID",
	target_column:str=None,
):
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	PROJECT_DIR = project_dir
	PROCESSED_READ_DIR = f"{PROJECT_DIR}ProcessedReads/"

	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)

	input_str = " ".join(["{i" + str(e) + "}" for e in range(len(prefixes))])
	files = {"i" + str(e): rms.file_from_path(f"{PROCESSED_READ_DIR}{s}.tsv.gz") for e, s in enumerate(prefixes)}
	name = "-".join(id_column.split(" ") + ([] if target_column is None or len(target_column) == 0 else target_column.split(" ")))
	if target_column is None or len(target_column) == 0:
		run_template_bash(
			"biodatatools join_IDMap " + input_str + " -o {o} -id_column {id_column}",
			conda_env=conda_env,
			output_func=output_func,
			**files,
			o=f"{PROCESSED_READ_DIR}{output_prefix}_{name}_map.tsv.gz", 
			id_column=id_column, 
		)
	else:
		run_template_bash(
			"biodatatools join_IDMap -i " + input_str + " -o {o} -id_column {id_column} -target_column {target_column}",
			conda_env=conda_env,
			output_func=output_func,
			**files,
			o=f"{PROCESSED_READ_DIR}{output_prefix}_{name}_map.tsv.gz", 
			id_column=id_column, 
			target_column=target_column, 
		)
def RMSTemplate_STARR_element_call_20240701(
	rms,
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	dna_prefix:list[str]=[],
	rna_prefix:list[str]=[],
	output_prefix:str="batch_output",
	neg_ctrl:str="Some file",
	filter_func:str=None,
	element_column:list[str]=None, count_column:str=None
):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	PROCESSED_READ_DIR = f"{PROJECT_DIR}ProcessedReads/"
	ELEMENT_ACTIVITY_DIR = f"{PROJECT_DIR}ElementActivity/"
	
	for dname in [PROJECT_DIR, PROCESSED_READ_DIR, ELEMENT_ACTIVITY_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	element_counts_file = f"{ELEMENT_ACTIVITY_DIR}{output_prefix}_Element_Counts.tsv.gz"
	logFC_file = f"{ELEMENT_ACTIVITY_DIR}{output_prefix}_limma_logFCs.tsv.gz"
	activity_call_file = f"{ELEMENT_ACTIVITY_DIR}{output_prefix}_limma_calls.tsv.gz"
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	idnak = ["{idna" + str(i) + "}" for i in range(len(dna_prefix))]
	irnak = ["{irna" + str(i) + "}" for i in range(len(rna_prefix))]
	idna_str = " ".join(idnak)
	irna_str = " ".join(irnak)
	idna = {"idna" + str(i): rms.file_from_path(f"{PROCESSED_READ_DIR}{p}.tsv.gz") for i, p in enumerate(dna_prefix)}
	irna = {"irna" + str(i): rms.file_from_path(f"{PROCESSED_READ_DIR}{p}.tsv.gz") for i, p in enumerate(rna_prefix)}
	additional_func_str = ""
	params_dict = {}
	if filter_func is not None:
		additional_func_str += " -filter_func {filter_func}"
		params_dict["filter_func"] = filter_func
	if element_column is not None:
		additional_func_str += " -element_column {element_column}"
		params_dict["element_column"] = " ".join(element_column)
	if count_column is not None:
		additional_func_str += " -count_column {count_column}"
		params_dict["count_column"] = count_column
	
	run_template_bash(
		"biodatatools merge_STARR_element_counts -idna " + idna_str + " -irna " + irna_str + " -o {o}" + additional_func_str,
		conda_env=conda_env,
		output_func=output_func,
		**idna,
		**irna,
		**params_dict,
		o=element_counts_file, 
	)
	run_template_bash(
		"biodatatools generate_STARR_logFC_by_limma -i {i} -n {n} -o {o}",
		conda_env=conda_env,
		output_func=output_func,
		i=rms.file_from_path(element_counts_file),
		n=rms.file_from_path(neg_ctrl),
		o=logFC_file, 
	)
	run_template_bash(
		"biodatatools generate_STARR_activity_call -i {i} -n {n} -o {o}",
		conda_env=conda_env,
		output_func=output_func,
		i=rms.file_from_path(logFC_file),
		n=rms.file_from_path(neg_ctrl),
		o=activity_call_file, 
	)
def RMSTemplate_STARR_element_change_call_20240701(
	rms,
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	output_prefix:str="batch_output",
	neg_ctrl:str="negative_controls.txt",
	pair:str="pairs.txt",
):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	ELEMENT_ACTIVITY_DIR = f"{PROJECT_DIR}ElementActivity/"
	
	for dname in [PROJECT_DIR, ELEMENT_ACTIVITY_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	logFC_file = f"{ELEMENT_ACTIVITY_DIR}{output_prefix}_limma_logFCs.tsv.gz"
	change_call_file = f"{ELEMENT_ACTIVITY_DIR}{output_prefix}_limma_change_calls.tsv.gz"
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"biodatatools generate_STARR_activity_change_call -i {i} -n {n} -p {p} -o {o}",
		conda_env=conda_env,
		output_func=output_func,
		i=rms.file_from_path(logFC_file),
		n=rms.file_from_path(neg_ctrl),
		p=rms.file_from_path(pair),
		o=change_call_file, 
	)

def RMSTemplate_generate_QC_20240701(
	rms,
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	prefix:str="batch_output",
	):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	ELEMENT_ACTIVITY_DIR = f"{PROJECT_DIR}ElementActivity/"
	REPLICATES_CORRELATION_DIR = f"{PROJECT_DIR}QC/ReplicatesCorrelation/"
	RANKED_LOGFC_DIR = f"{PROJECT_DIR}QC/RankedLogFC/"
	ORIENTATION_CORRELATION_DIR = f"{PROJECT_DIR}QC/OrientationCorrelation/"
	for dname in [PROJECT_DIR, REPLICATES_CORRELATION_DIR, RANKED_LOGFC_DIR, ORIENTATION_CORRELATION_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	logFC_file = f"{ELEMENT_ACTIVITY_DIR}{prefix}_limma_logFCs.tsv.gz"
	activity_call_file = f"{ELEMENT_ACTIVITY_DIR}{prefix}_limma_calls.tsv.gz"
		
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f) # All output_funcs are the same
	
	run_template_bash(
		"biodatatools plot_STARR_replicates_correlation -i {i} -o {o}",
		conda_env=conda_env,
		output_func=output_func,
		i=rms.file_from_path(logFC_file),
		o=f"{REPLICATES_CORRELATION_DIR}{prefix}_replicates_correlation.png",
	)

	
	run_template_bash(
		"biodatatools plot_STARR_orientation_correlation -i {i} -o {o}",
		conda_env=conda_env,
		output_func=output_func,
		i=rms.file_from_path(activity_call_file),
		o=f"{ORIENTATION_CORRELATION_DIR}{prefix}_orientation_correlation.png",
	)

def RMSTemplate_generate_shared_ranked_logFCs_20240701(
	rms,
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	prefix:list[str]=[],
	names:list[str]=[],
	output_prefix:str="output",
	default_group_name:str=None,
	group:str=None,
	plot_kw:str=None,
	plot_kw_dict:str=None,
	
	):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	ELEMENT_ACTIVITY_DIR = f"{PROJECT_DIR}ElementActivity/"
	RANKED_LOGFC_DIR = f"{PROJECT_DIR}QC/RankedLogFC/"
	for dname in [PROJECT_DIR, RANKED_LOGFC_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	
	params = {}
	if default_group_name is not None:
		params["default_group_name"] = default_group_name
	if group is not None:
		params["group"] = group
	if plot_kw is not None:
		params["plot_kw"] = plot_kw
	if plot_kw_dict is not None:
		params["plot_kw_dict"] = plot_kw_dict
	param_str = "".join(["-"+k+" {" + k + "} " for k in params.keys()])
	
	input_files = [f"{ELEMENT_ACTIVITY_DIR}{i}_limma_calls.tsv.gz" for i in prefix]	
	dinputs = {"i" + str(i): f for i, f in enumerate(input_files)}
	dnames = {"n" + str(i): f for i, f in enumerate(names)}
	inputstr = " ".join(["{" + k + "}" for k in dinputs.keys()])
	namestr = " ".join(["{" + k + "}" for k in dnames.keys()])
	
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f) # All output_funcs are the same
	
	run_template_bash(
		"biodatatools plot_STARR_ranked_logFCs -i " + inputstr + " -n " + namestr + " " + param_str + "-o {o}",
		conda_env=conda_env,
		output_func=output_func,
		**dinputs,
		**dnames,
		**params,
		o=f"{RANKED_LOGFC_DIR}{output_prefix}_ranked_logFCs.png",
	)

def RMSTemplate_generate_IDMap_QC_20240701(
	rms,
	conda_env: str = "rms_starr_analysis_20240701",
	project_dir: str = "/home/happyuser/STARR/",
	prefix:str="sample",
	id_column:str="ID",
	target_column:str="", 
	refs:dict=None,
	default_group_name:str=None,
	group:str=None,
	plot_kw:str=None,
	plot_kw_dict:str=None,
	
	):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)

	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	PROCESSED_READ_DIR = f"{PROJECT_DIR}ProcessedReads/"
	IDMAP_QC_DIR = f"{PROJECT_DIR}QC/IDMapQC/"
	for dname in [PROJECT_DIR, IDMAP_QC_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	if refs is None or len(refs) == 0:
		raise Exception()
	params = {}
	if plot_kw is not None:
		params["plot_kw"] = plot_kw
	if plot_kw_dict is not None:
		params["plot_kw_dict"] = plot_kw_dict
	param_str = "".join([" -"+k+" {" + k + "}" for k in params.keys()])
	
	name = "-".join(id_column.split(" ") + ([] if target_column is None or len(target_column) == 0 else target_column.split(" ")))
	dinputs = {"r" + str(i): rms.file_from_path(f) for i, f in enumerate(refs.values())}
	inputstr = " ".join(["{" + k + "}" for k in dinputs.keys()])
	
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f) # All output_funcs are the same
	
	run_template_bash(
		"biodatatools plot_IDMap_target_ID_representations -i {i} -o {o} -target_column {target_column} -r " + inputstr + param_str,
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f'{PROCESSED_READ_DIR}{prefix}_{name}_map.tsv.gz'),
		target_column=target_column,
		**dinputs,
		**params,
		o=f"{IDMAP_QC_DIR}{prefix}_IDMap-representations.png",
	)
	
	run_template_bash(
		"biodatatools plot_ID_targets_pair_dominance -i {i} -o {o} -id_column {id_column} -target_column {target_column}" + param_str,
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f'{PROCESSED_READ_DIR}{prefix}.tsv.gz'),
		id_column=id_column,
		target_column=id_column + " " + target_column,
		**params,
		o=f"{IDMAP_QC_DIR}{prefix}_IDMap-dominance.png",
	)
