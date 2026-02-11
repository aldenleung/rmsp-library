# A preprocessing scripts that perform steps including the following
## STAR genome indexing
## Pre-processing 
### Data trimming
### Alignment
### Alignment filtering
### Alignment deduplication
## Peak calling
## Generate statistics and QC files

# A typical set up
# Resources/
#	Genome/
#		Human/
# 	Gencode/
#		Human/
# 	STAR_genome_index/
# PROcap/
#	RawData/
#		FastQC_reports/
#		fastp_reports/
#	TrimmedData/
#		FastQC_reports/
#	RawAlignments/ 			# This folder contains all bam files
#		bam_stat_reports/
#		umi_tools_reports/
#	Alignments/				# This folder contains all bigwig files and distance bedgraph files
#	Peaks
#		PINTS/			
#			Merged/
#		EnhancerNet/	
# 	QC/
#		Annotations/
# 		GeneBodyRatio/
# 		ReplicatesCorrelation/
#		RNALen/

from rmsp.rmstemplate import *
from commonhelper import convert_to_bool

def RMSTemplate_reference_download_20240501(
	rms,
	genome_link: str="https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz",
	annotation_link: str="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz",
	chrom_size_link: str="https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes",
	annotation_dir: str="Resources/Gencode/Human/",
	genome_dir: str="Resource/Genome/Human/",
	):
	requests.get()
	
	

def RMSTemplate_genome_indexing_20240501(
	rms,
	conda_env: str = "rms_procap_analysis_20240501",
	genome_file: InputFileType = "hg38.fa",
	annotation_file: InputFileType = "hg38.gtf",
	o: OutputFileType="STAR_genome_index/hg38",
	thread: int = 1,
	):
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	
	run_template_bash("STAR --runMode genomeGenerate --readFilesCommand zcat --runThreadN {thread} --genomeDir {genomeDir} --genomeFastaFiles {genome_file} --sjdbGTFfile {annotation_file}",
					 output_func=_OUTPUT_STAR_genome_generate_20210927,
					 genome_file=genome_file,
					 annotation_file=annotation_file,
					 genomeDir=genomeDir
	)
def RMSTemplate_generate_QC_gbratio_annotation_20240501(
	rms,
	conda_env: str = "rms_procap_analysis_20240501",
	genome_file: InputFileType = "hg38.fa",
	annotation_file: InputFileType = "hg38.gtf",
	o: OutputFileType="QC/Annotations/hg38_qc-gbratio.gtf",
	
	remove_overlapping_genes: convert_to_bool=True,
	min_transcript_size:int = 2000,
	discarded_chroms:list[str]=[],
	
	):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	QC_ANNO_DIR = os.path.dirname(o)
	for dname in [QC_ANNO_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	
	dc_dict = {f"dc{idx}":c for idx, c in enumerate(discarded_chroms)}
	dc_str = ",".join([f"'{{{k}}}'" for k in dc_dict.keys()]) 
	run_template_bash(
		"biodatatools filter_geneannotations -i {i} -o {o} -remove_overlapping_genes {remove_overlapping_genes} -filter_func \"lambda gff: gff.feature == 'transcript' and gff.attributes['gene_type'] == 'protein_coding' and len(gff.genomic_pos) >= {min_transcript_size} and gff.genomic_pos.name not in [" + dc_str + "]\"",
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(annotation_file),
		o = o,
		remove_overlapping_genes = remove_overlapping_genes,
		min_transcript_size = min_transcript_size,
		**dc_dict
	)		
def RMSTemplate_generate_QC_gbratio_annotation_20240901(
	rms,
	conda_env: str = "rms_procap_analysis_20240501",
	genome_file: InputFileType = "hg38.fa",
	annotation_file: InputFileType = "hg38.gtf",
	o: OutputFileType="QC/Annotations/hg38_qc-gbratio.gtf",
	
	remove_overlapping_genes: convert_to_bool=True,
	min_transcript_size:int = 2000,
	discarded_chroms:list[str]=[],
	
	):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	QC_ANNO_DIR = os.path.dirname(o)
	for dname in [QC_ANNO_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	
	dc_dict = {f"dc{idx}":c for idx, c in enumerate(discarded_chroms)}
	dc_str = ",".join([f"'{{{k}}}'" for k in dc_dict.keys()]) 
	run_template_bash(
		"biodatatools filter_geneannotations -i {i} -o {o} -remove_overlapping_genes {remove_overlapping_genes} -filter_func \"lambda gff: gff.feature == 'transcript' and ((gff.attributes['gene_type'] if 'gene_type' in gff.attributes else gff.attributes['gene_biotype']) == 'protein_coding') and len(gff.genomic_pos) >= {min_transcript_size} and gff.genomic_pos.name not in [" + dc_str + "]\"",
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(annotation_file),
		o = o,
		remove_overlapping_genes = remove_overlapping_genes,
		min_transcript_size = min_transcript_size,
		**dc_dict
	)		

	
def RMSTemplate_preprocessing_without_UMI_20240501(
	rms, 
	conda_env: str = "rms_procap_analysis_20240501",
	
	project_dir: str = "/home/happyuser/PROcap/",
	prefix: str="sample",
	
	# Inputs/Outputs
	genome_dir : InputFileType = "/home/happyuser/STAR_genome_index/hg38/", 
	chrom_size_file: InputFileType = "/home/happyuser/Genome/hg38.chrom.sizes", 
	r1: InputFileType = "/home/happyuser/PROcap/RawData/example_read1.fq.gz", 
	r2: InputFileType = "/home/happyuser/PROcap/RawData/example_read2.fq.gz", 
	
	# Experiment set up
	rna_strand: str = "forward",

	# General options
	thread: int = 1,
	# Other options
	remove_spliced_reads: bool = True
):
	'''
	
	'''	
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	RAW_DATA_DIR = f"{PROJECT_DIR}RawData/"
	TRIMMED_DATA_DIR = f"{PROJECT_DIR}TrimmedData/"
	RAW_ALIGNMENT_DIR = f"{PROJECT_DIR}RawAlignments/"
	ALIGNMENT_DIR = f"{PROJECT_DIR}Alignments/"
	
	RAW_DATA_FASTQC_DIR = f"{RAW_DATA_DIR}FastQC_reports/"
	TRIMMED_DATA_FASTQC_DIR = f"{TRIMMED_DATA_DIR}FastQC_reports/"
	FASTP_REPORT_DIR = f"{RAW_DATA_DIR}fastp_reports/"
	UMI_TOOLS_REPORT_DIR = f"{RAW_ALIGNMENT_DIR}umi_tools_reports/"
	BAM_STAT_DIR = f"{RAW_ALIGNMENT_DIR}bam_stat_reports/"
	for dname in [PROJECT_DIR, RAW_DATA_DIR, TRIMMED_DATA_DIR, RAW_ALIGNMENT_DIR, ALIGNMENT_DIR, RAW_DATA_FASTQC_DIR, TRIMMED_DATA_FASTQC_DIR, FASTP_REPORT_DIR, UMI_TOOLS_REPORT_DIR, BAM_STAT_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
		
	# Trimming
	t1 = f"{TRIMMED_DATA_DIR}{prefix}_1.fq.gz"
	t2 = f"{TRIMMED_DATA_DIR}{prefix}_2.fq.gz"
	fastp_html_report = f"{FASTP_REPORT_DIR}{prefix}.html"
	fastp_json_report = f"{FASTP_REPORT_DIR}{prefix}.json"
	f = lambda **kwargs: [kwargs[k] for k in ["o1", "o2", "html_report", "json_report"] if k in kwargs]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"fastp -i {i1} -I {i2} -o {o1} -O {o2} -h {html_report} -j {json_report} --thread {thread} --overlap_len_require {overlap_len_require} --length_required {length_required} --low_complexity_filter -g -p -c",
		conda_env=conda_env,
		output_func=output_func,
		i1=rms.file_from_path(r1),
		i2=rms.file_from_path(r2),
		o1=t1,
		o2=t2,
		html_report=fastp_html_report,
		json_report=fastp_json_report,
		thread=thread,
		overlap_len_require=18,
		length_required=18,
	)
	
	# Quality check for both raw and trimmed data
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
	
	# Alignments
	f = lambda **kwargs: [f"{kwargs['output_prefix']}Aligned.sortedByCoord.out.bam", f"{kwargs['output_prefix']}Log.final.out", f"{kwargs['output_prefix']}Log.out", f"{kwargs['output_prefix']}Log.progress.out", f"{kwargs['output_prefix']}ReadsPerGene.out.tab", f"{kwargs['output_prefix']}SJ.out.tab"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"STAR --readFilesCommand zcat --runThreadN {thread} --alignMatesGapMax 1000 --outFilterMultimapNmax 10 --outFilterMismatchNmax 1 --outFilterMultimapScoreRange 0 --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000 --quantMode GeneCounts --outSAMattributes All --genomeDir {genome_dir} --outFileNamePrefix {output_prefix} --readFilesIn {i1} {i2}", 
		conda_env=conda_env,
		output_func=output_func,
		genome_dir=rms.file_from_path(genome_dir),
		output_prefix=f"{RAW_ALIGNMENT_DIR}{prefix}_",
		i1=rms.file_from_path(t1), 
		i2=rms.file_from_path(t2),
		thread=thread,
	)
	
	# Filter uniquely mapped alignments
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	if remove_spliced_reads:
		run_template_bash(
			"TMPFILE=$(mktemp --suffix .bam); samtools view -hb -q 255 -f 3 -@ {thread} -o $TMPFILE {i}; biodatatools filter_bam_NCIGAR_reads -i $TMPFILE -o {o} -nthread {thread}; rm $TMPFILE",
			conda_env=conda_env,
			output_func=output_func,
			thread=thread,
			i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_Aligned.sortedByCoord.out.bam"),
			o = f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam"
			)
	else:
		run_template_bash(
			"samtools view -hb -q 255 -f 3 -@ {thread} -o {o} {i}", 
			conda_env=conda_env,
			output_func=output_func,
			thread=thread,
			i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_Aligned.sortedByCoord.out.bam"),
			o = f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam"
			)
	
	f = lambda **kwargs: [f"{kwargs['i']}.bai"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"samtools index {i}", 
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam"),
	)
	
	# bam statistics
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	for suffix in ["_Aligned.sortedByCoord.out", "_uniquely_mapped"]:
		run_template_bash(
			"samtools coverage {i} > {o}", 
			conda_env=conda_env,
			output_func=output_func,
			i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}{suffix}.bam"),
			o = f"{BAM_STAT_DIR}{prefix}{suffix}.bam.stat.txt",
		)
	
	# Conversion to bigwig and distance bed file
	f = lambda **kwargs: [kwargs['o'] + "_5pl.bw", kwargs['o'] + "_5mn.bw", kwargs['o'] + "_3pl.bw", kwargs['o'] + "_3mn.bw"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"biodatatools process_PROcap_bam_to_bigwig -i {i} -g {g} -o {o} -paired_end True -rna_strand {rna_strand}",
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam"),
		g = rms.file_from_path(chrom_size_file),
		o = f"{ALIGNMENT_DIR}{prefix}",
		rna_strand = rna_strand,
	)
	f = lambda **kwargs: [kwargs['o'] + "_dpl.bed.bgz", kwargs['o'] + "_dmn.bed.bgz"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"biodatatools process_PROcap_bam_to_TSS_RNA_len -i {i} -o {o} -paired_end True -rna_strand {rna_strand} -g {g}",
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam"),
		o = f"{ALIGNMENT_DIR}{prefix}",
		rna_strand = rna_strand,
		g = rms.file_from_path(chrom_size_file),
	)
	
	f = lambda **kwargs: [kwargs['i'] + ".tbi"]
	output_func = rms.register_pipe(f)	 
	for suffix in ["_dpl", "_dmn"]:
		run_template_bash(
			"tabix -s 1 -b 2 -e 3 -0 {i}",
			conda_env=conda_env,
			output_func=output_func,
			i = rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}{suffix}.bed.bgz")
		)

def RMSTemplate_preprocessing_with_UMI_20240501(
	rms, 
	conda_env: str = "rms_procap_analysis_20240501",
	
	project_dir: str = "/home/happyuser/PROcap/",
	prefix: str="sample",
	
	# Inputs/Outputs
	genome_dir : InputFileType = "/home/happyuser/STAR_genome_index/hg38/", 
	chrom_size_file: InputFileType = "/home/happyuser/Genome/hg38.chrom.sizes", 
	r1: InputFileType = "/home/happyuser/PROcap/RawData/example_read1.fq.gz", 
	r2: InputFileType = "/home/happyuser/PROcap/RawData/example_read2.fq.gz", 
	
	# Experiment set up
	rna_strand: str = "forward",
	umi_len: int = 6,
	umi_loc: str = "per_read",
	adapter_seq_r1: str = "ACGGAT", 
	adapter_seq_r2: str = "CGGCTT",

	# General options
	thread: int = 1,
	# Other options
	remove_spliced_reads: bool = True
):
	'''
	
	'''	
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	RAW_DATA_DIR = f"{PROJECT_DIR}RawData/"
	TRIMMED_DATA_DIR = f"{PROJECT_DIR}TrimmedData/"
	RAW_ALIGNMENT_DIR = f"{PROJECT_DIR}RawAlignments/"
	ALIGNMENT_DIR = f"{PROJECT_DIR}Alignments/"
	
	RAW_DATA_FASTQC_DIR = f"{RAW_DATA_DIR}FastQC_reports/"
	TRIMMED_DATA_FASTQC_DIR = f"{TRIMMED_DATA_DIR}FastQC_reports/"
	FASTP_REPORT_DIR = f"{RAW_DATA_DIR}fastp_reports/"
	UMI_TOOLS_REPORT_DIR = f"{RAW_ALIGNMENT_DIR}umi_tools_reports/"
	BAM_STAT_DIR = f"{RAW_ALIGNMENT_DIR}bam_stat_reports/"
	for dname in [PROJECT_DIR, RAW_DATA_DIR, TRIMMED_DATA_DIR, RAW_ALIGNMENT_DIR, ALIGNMENT_DIR, RAW_DATA_FASTQC_DIR, TRIMMED_DATA_FASTQC_DIR, FASTP_REPORT_DIR, UMI_TOOLS_REPORT_DIR, BAM_STAT_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
		
	# Trimming
	t1 = f"{TRIMMED_DATA_DIR}{prefix}_1.fq.gz"
	t2 = f"{TRIMMED_DATA_DIR}{prefix}_2.fq.gz"
	fastp_html_report = f"{FASTP_REPORT_DIR}{prefix}.html"
	fastp_json_report = f"{FASTP_REPORT_DIR}{prefix}.json"
	f = lambda **kwargs: [kwargs[k] for k in ["o1", "o2", "html_report", "json_report"] if k in kwargs]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"fastp -i {i1} -I {i2} -o {o1} -O {o2} -h {html_report} -j {json_report} --umi --umi_len={umi_len} --umi_loc={umi_loc} --thread {thread} --adapter_sequence {adapter_sequence} --adapter_sequence_r2 {adapter_sequence_r2} --overlap_len_require {overlap_len_require} --length_required {length_required} --low_complexity_filter -g -p -c",
		conda_env=conda_env,
		output_func=output_func,
		i1=rms.file_from_path(r1),
		i2=rms.file_from_path(r2),
		o1=t1,
		o2=t2,
		html_report=fastp_html_report,
		json_report=fastp_json_report,
		umi_len=umi_len,
		umi_loc=umi_loc,
		thread=thread,
		adapter_sequence=adapter_seq_r1,
		adapter_sequence_r2=adapter_seq_r2,
		overlap_len_require=18,
		length_required=18,
	)
	
	# Quality check for both raw and trimmed data
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
	
	# Alignments
	f = lambda **kwargs: [f"{kwargs['output_prefix']}Aligned.sortedByCoord.out.bam", f"{kwargs['output_prefix']}Log.final.out", f"{kwargs['output_prefix']}Log.out", f"{kwargs['output_prefix']}Log.progress.out", f"{kwargs['output_prefix']}ReadsPerGene.out.tab", f"{kwargs['output_prefix']}SJ.out.tab"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"STAR --readFilesCommand zcat --runThreadN {thread} --alignMatesGapMax 1000 --outFilterMultimapNmax 10 --outFilterMismatchNmax 1 --outFilterMultimapScoreRange 0 --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000 --quantMode GeneCounts --outSAMattributes All --genomeDir {genome_dir} --outFileNamePrefix {output_prefix} --readFilesIn {i1} {i2}", 
		conda_env=conda_env,
		output_func=output_func,
		genome_dir=rms.file_from_path(genome_dir),
		output_prefix=f"{RAW_ALIGNMENT_DIR}{prefix}_",
		i1=rms.file_from_path(t1), 
		i2=rms.file_from_path(t2),
		thread=thread,
	)
	
	# Filter uniquely mapped alignments
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	if remove_spliced_reads:
		run_template_bash(
			#"samtools view -h -F 4 {i} | awk '$0 ~ /^@/ || $6 !~ /N/' | samtools view -hbS -q 255 -f 3 -@ {thread} -o {o} -",
			"TMPFILE=$(mktemp --suffix .bam); samtools view -hb -q 255 -f 3 -@ {thread} -o $TMPFILE {i}; biodatatools filter_bam_NCIGAR_reads -i $TMPFILE -o {o} -nthread {thread}; rm $TMPFILE",
			conda_env=conda_env,
			output_func=output_func,
			thread=thread,
			i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_Aligned.sortedByCoord.out.bam"),
			o = f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam"
			)
	else:
		run_template_bash(
			"samtools view -hb -q 255 -f 3 -@ {thread} -o {o} {i}", 
			conda_env=conda_env,
			output_func=output_func,
			thread=thread,
			i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_Aligned.sortedByCoord.out.bam"),
			o = f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam"
			)
	
	f = lambda **kwargs: [f"{kwargs['i']}.bai"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"samtools index {i}", 
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam"),
	)
	
	# Reads deduplication based on UMIs
	f = lambda **kwargs: [kwargs['o'], kwargs['output_stats'] + "_per_umi_per_position.tsv", kwargs['output_stats'] + "_edit_distance.tsv", kwargs['output_stats'] + "_per_umi.tsv"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"umi_tools dedup --unpaired-reads=discard --umi-separator=: --paired -I {i} --output-stats={output_stats} -S {o}", 
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam"),
		bami = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam.bai"),
		o = f"{RAW_ALIGNMENT_DIR}{prefix}_dedup.bam",
		output_stats=f"{UMI_TOOLS_REPORT_DIR}{prefix}"
	)

	f = lambda **kwargs: [f"{kwargs['i']}.bai"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"samtools index {i}", 
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_dedup.bam"),
	)
	
	# bam statistics
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	for suffix in ["_Aligned.sortedByCoord.out", "_uniquely_mapped", "_dedup"]:
		run_template_bash(
			"samtools coverage {i} > {o}", 
			conda_env=conda_env,
			output_func=output_func,
			i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}{suffix}.bam"),
			o = f"{BAM_STAT_DIR}{prefix}{suffix}.bam.stat.txt",
		)
	
	# Conversion to bigwig and distance bed file
	f = lambda **kwargs: [kwargs['o'] + "_5pl.bw", kwargs['o'] + "_5mn.bw", kwargs['o'] + "_3pl.bw", kwargs['o'] + "_3mn.bw"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"biodatatools process_PROcap_bam_to_bigwig -i {i} -g {g} -o {o} -paired_end True -rna_strand {rna_strand}",
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_dedup.bam"),
		g = rms.file_from_path(chrom_size_file),
		o = f"{ALIGNMENT_DIR}{prefix}",
		rna_strand = rna_strand,
	)
	f = lambda **kwargs: [kwargs['o'] + "_dpl.bed.bgz", kwargs['o'] + "_dmn.bed.bgz"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"biodatatools process_PROcap_bam_to_TSS_RNA_len -i {i} -o {o} -paired_end True -rna_strand {rna_strand} -g {g}",
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_dedup.bam"),
		o = f"{ALIGNMENT_DIR}{prefix}",
		rna_strand = rna_strand,
		g = rms.file_from_path(chrom_size_file),
	)
	
	f = lambda **kwargs: [kwargs['i'] + ".tbi"]
	output_func = rms.register_pipe(f)	 
	for suffix in ["_dpl", "_dmn"]:
		run_template_bash(
			"tabix -s 1 -b 2 -e 3 -0 {i}",
			conda_env=conda_env,
			output_func=output_func,
			i = rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}{suffix}.bed.bgz")
		)
	
def RMSTemplate_merge_reseqs_20240501(
	rms,
	conda_env: str = "rms_procap_analysis_20240501",
	project_dir: str="/home/happyuser/PROcap/",
	prefixes: list[str] = ["sample1_rep1_reseq1", "sample1_rep1_reseq2"],
	chrom_size_file: str = "/home/happyuser/Resources/Genome/Human/hg38.chrom.sizes",
	rna_strand:str = "forward"
	):
	# Merge resequencing data. Since it is from the same replicate but sequenced multiple times, we need to start merging uniquely mapped alignments and redo deduplication. 
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	RAW_ALIGNMENT_DIR = f"{PROJECT_DIR}RawAlignments/"
	ALIGNMENT_DIR = f"{PROJECT_DIR}Alignments/"
	UMI_TOOLS_REPORT_DIR = f"{RAW_ALIGNMENT_DIR}umi_tools_reports/"
	BAM_STAT_DIR = f"{RAW_ALIGNMENT_DIR}bam_stat_reports/"
	for dname in [PROJECT_DIR, RAW_ALIGNMENT_DIR, ALIGNMENT_DIR, UMI_TOOLS_REPORT_DIR, BAM_STAT_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	
	KEY_MIDDLE = "_and_"
	KEY_BEFORE = "bsm_" # Begin reSeq Merge
	KEY_AFTER = "_esm" # End reSeq Merge
	
	prefixes = sorted(prefixes)
	output_prefix = KEY_BEFORE + KEY_MIDDLE.join(prefixes) + KEY_AFTER
	
	# Merging
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"samtools merge -o {o} " + " ".join([f"{{i{idx}}}" for idx in range(len(prefixes))]), 
		conda_env=conda_env,
		output_func=output_func,
		**{f"i{idx}": rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam") for idx, prefix in enumerate(prefixes)},
		o=f"{RAW_ALIGNMENT_DIR}{output_prefix}_uniquely_mapped.bam"
	)

	f = lambda **kwargs: [f"{kwargs['i']}.bai"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"samtools index {i}", 
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{output_prefix}_uniquely_mapped.bam"),
	)
	
	# The codes below will be the same as normal processing routine, except raw alignment stat is not generated
	prefix = output_prefix
	# Reads deduplication based on UMIs
	f = lambda **kwargs: [kwargs['o'], kwargs['output_stats'] + "_per_umi_per_position.tsv", kwargs['output_stats'] + "_edit_distance.tsv", kwargs['output_stats'] + "_per_umi.tsv"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"umi_tools dedup --unpaired-reads=discard --umi-separator=: --paired -I {i} --output-stats={output_stats} -S {o}", 
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam"),
		bami = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_uniquely_mapped.bam.bai"),
		o = f"{RAW_ALIGNMENT_DIR}{prefix}_dedup.bam",
		output_stats=f"{UMI_TOOLS_REPORT_DIR}{prefix}"
	)

	f = lambda **kwargs: [f"{kwargs['i']}.bai"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"samtools index {i}", 
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_dedup.bam"),
	)
		
	# bam statistics
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	for suffix in ["_uniquely_mapped", "_dedup"]: # No raw alignment stat
		run_template_bash(
			"samtools coverage {i} > {o}", 
			conda_env=conda_env,
			output_func=output_func,
			i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}{suffix}.bam"),
			o = f"{BAM_STAT_DIR}{prefix}{suffix}.bam.stat.txt",
		)
	
	# Conversion to bigwig and distance bed file
	f = lambda **kwargs: [kwargs['o'] + "_5pl.bw", kwargs['o'] + "_5mn.bw", kwargs['o'] + "_3pl.bw", kwargs['o'] + "_3mn.bw"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"biodatatools process_PROcap_bam_to_bigwig -i {i} -g {g} -o {o} -paired_end True -rna_strand {rna_strand}",
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_dedup.bam"),
		g = rms.file_from_path(chrom_size_file),
		o = f"{ALIGNMENT_DIR}{prefix}",
		rna_strand = rna_strand,
	)
	f = lambda **kwargs: [kwargs['o'] + "_dpl.bed.bgz", kwargs['o'] + "_dmn.bed.bgz"]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"biodatatools process_PROcap_bam_to_TSS_RNA_len -i {i} -o {o} -paired_end True -rna_strand {rna_strand} -g {g}",
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{RAW_ALIGNMENT_DIR}{prefix}_dedup.bam"),
		o = f"{ALIGNMENT_DIR}{prefix}",
		rna_strand = rna_strand,
		g=rms.file_from_path(chrom_size_file)
	)
	f = lambda **kwargs: [kwargs['i'] + ".tbi"]
	output_func = rms.register_pipe(f)	 
	for suffix in ["_dpl", "_dmn"]:
		run_template_bash(
			"tabix -s 1 -b 2 -e 3 -0 {i}",
			conda_env=conda_env,
			output_func=output_func,
			i = rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}{suffix}.bed.bgz")
		)


def RMSTemplate_merge_replicates_20240501(
	rms,
	conda_env: str = "rms_procap_analysis_20240501",
	project_dir: str="/home/happyuser/PROcap/",
	prefixes: list[str] = ["sample1_rep1", "sample1_rep2"],
	chrom_size_file: str = "/home/happyuser/Resources/Genome/Human/hg38.chrom.sizes",
	):
	# Merge bigwig and distance bed files
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	ALIGNMENT_DIR = f"{PROJECT_DIR}Alignments/"
	for dname in [PROJECT_DIR, ALIGNMENT_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	
	KEY_MIDDLE = "_and_"
	KEY_BEFORE = "brm_" # Begin Replicate Merge
	KEY_AFTER = "_erm" # End Replicate Merge
	
	prefixes = sorted(prefixes)
	output_prefix = KEY_BEFORE + KEY_MIDDLE.join(prefixes) + KEY_AFTER
	
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f) # All output_funcs are the same for the merging
	
	for suffix in "_5pl", "_3pl":
		input_files = {f"i{idx}":rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}{suffix}.bw") for idx, prefix in enumerate(prefixes)}
		output_file = f"{ALIGNMENT_DIR}{output_prefix}{suffix}.bw"
		run_template_bash(
			"biodatatools merge_bigwig -i " + " ".join(["{" + f"i{idx}" + "}" for idx in range(len(prefixes))]) + " -g {g} -o {o} -autosort {autosort} -filter_chr {filter_chr}",
			conda_env=conda_env,
			output_func=output_func,
			**input_files,
			g=rms.file_from_path(chrom_size_file),
			o=output_file,
			autosort=True,
			filter_chr=True,
		)
	for suffix in "_5mn", "_3mn":
		input_files = {f"i{idx}":rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}{suffix}.bw") for idx, prefix in enumerate(prefixes)}
		output_file = f"{ALIGNMENT_DIR}{output_prefix}{suffix}.bw"
		run_template_bash(
			"biodatatools merge_bigwig -i " + " ".join(["{" + f"i{idx}" + "}" for idx in range(len(prefixes))]) + " -g {g} -o {o} -threshold {threshold} -remove_zero {remove_zero} -autosort {autosort} -filter_chr {filter_chr}",
			conda_env=conda_env,
			output_func=output_func,
			**input_files,
			g=rms.file_from_path(chrom_size_file),
			o=output_file,
			threshold=-2147483648,
			remove_zero=True,
			autosort=True,
			filter_chr=True,
		)
	for suffix in "_dpl", "_dmn":
		input_files = {f"i{idx}":rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}{suffix}.bed.bgz") for idx, prefix in enumerate(prefixes)}
		output_file = f"{ALIGNMENT_DIR}{output_prefix}{suffix}.bed.bgz"
		run_template_bash(
			"biodatatools merge_PROcap_TSS_RNA_len -i " + " ".join(["{" + f"i{idx}" + "}" for idx in range(len(prefixes))]) + " -o {o}",
			conda_env=conda_env,
			output_func=output_func,
			**input_files,
			o=output_file,
		)
	f = lambda **kwargs: [kwargs['i'] + ".tbi"]
	output_func = rms.register_pipe(f)	 
	for suffix in ["_dpl", "_dmn"]:
		run_template_bash(
			"tabix -s 1 -b 2 -e 3 -0 {i}",
			conda_env=conda_env,
			output_func=output_func,
			i = rms.file_from_path(f"{ALIGNMENT_DIR}{output_prefix}{suffix}.bed.bgz")
		)

		
	

def RMSTemplate_peak_calling_pints_20240501(
	rms, 
	conda_env: str = "rms_procap_analysis_20240501",
	project_dir: str="/home/happyuser/PROcap/",
	prefix:str="sample1", 
	exp_type:str="PROcap", 
	thread:int=1
	):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	ALIGNMENT_DIR = f"{PROJECT_DIR}Alignments/"
	PEAK_PINTS_DIR = f"{PROJECT_DIR}Peaks/PINTS/"
	for dname in [PROJECT_DIR, ALIGNMENT_DIR, PEAK_PINTS_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})

	if exp_type == "PRO-cap" or exp_type == "PROcap":
		exp_type = "PROcap"
	elif exp_type == "PRO-seq" or exp_type == "PROseq":
		exp_type = "PROseq"

	f = lambda **kwargs: [f"{kwargs['output_dir']}/{kwargs['file_prefix']}_1_bidirectional_peaks.bed", f"{kwargs['output_dir']}/{kwargs['file_prefix']}_1_divergent_peaks.bed", f"{kwargs['output_dir']}/{kwargs['file_prefix']}_1_unidirectional_peaks.bed"]
	output_func = rms.register_pipe(f)

	if exp_type == "PROseq":
		bwpl = f"{ALIGNMENT_DIR}{prefix}_3pl.bw"
		bwmn = f"{ALIGNMENT_DIR}{prefix}_3mn.bw"
	else:
		bwpl = f"{ALIGNMENT_DIR}{prefix}_5pl.bw"
		bwmn = f"{ALIGNMENT_DIR}{prefix}_5mn.bw"
	run_template_bash(
		"pints_caller --save-to {output_dir} --file-prefix {file_prefix} --bw-pl {bwpl} --bw-mn {bwmn} --exp-type {exp_type} --min-lengths-opposite-peaks 5 --thread {thread}", 
		conda_env=conda_env, 
		output_func=output_func,
		output_dir=PEAK_PINTS_DIR,
		file_prefix=prefix, 
		bwpl=rms.file_from_path(bwpl),
		bwmn=rms.file_from_path(bwmn),
		exp_type=exp_type,
		thread=thread
	)
def RMSTemplate_generate_replicates_correlation_20240501(
	rms,
	conda_env: str = "rms_procap_analysis_20240501",
	project_dir: str="/home/happyuser/PROcap/",
	prefixes: list[str] = ["sample1_rep1", "sample1_rep2"],
	bin_size: int=50,
	chrom_size_file: str = "/home/happyuser/Resources/Genome/Human/hg38.chrom.sizes.selected.chrs",
	):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	ALIGNMENT_DIR = f"{PROJECT_DIR}Alignments/"
	REPLICATES_CORRELATION_DIR = f"{PROJECT_DIR}QC/ReplicatesCorrelation/"
	for dname in [PROJECT_DIR, ALIGNMENT_DIR, REPLICATES_CORRELATION_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
		
	KEY_MIDDLE = "_and_"
	KEY_BEFORE = "brm_" # Begin Replicate Merge
	KEY_AFTER = "_erm" # End Replicate Merge
	prefixes = sorted(prefixes)
	output_prefix = KEY_BEFORE + KEY_MIDDLE.join(prefixes) + KEY_AFTER
		
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f) # All output_funcs are the same
	
	for suffix in "_5pl", "_5mn":
		sample_names = {f"s{idx}":prefix for idx, prefix in enumerate(prefixes)}
		input_files = {f"i{idx}":rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}{suffix}.bw") for idx, prefix in enumerate(prefixes)}
		output_file = f"{REPLICATES_CORRELATION_DIR}{output_prefix}{suffix}_binned_count_table.txt.gz"
	
		run_template_bash(
			"biodatatools process_bigwigs_to_count_table -sample_names " + " ".join(["{" + f"s{idx}" + "}" for idx in range(len(prefixes))]) + " -i " + " ".join(["{" + f"i{idx}" + "}" for idx in range(len(prefixes))]) + " -o {o} -bin_size {bin_size} -g {g}",
			conda_env=conda_env,
			output_func=output_func,
			**sample_names,
			**input_files,
			o=output_file,
			bin_size=bin_size,
			g=rms.file_from_path(chrom_size_file)
	)
	run_template_bash(
		"biodatatools process_count_tables_to_correlation_table -i {i0} {i1} -o {o} -filter_func {filter_func} -value_func {value_func}",
		conda_env=conda_env,
		output_func=output_func,
		i0 = rms.file_from_path(f"{REPLICATES_CORRELATION_DIR}{output_prefix}_5pl_binned_count_table.txt.gz"),
		i1 = rms.file_from_path(f"{REPLICATES_CORRELATION_DIR}{output_prefix}_5mn_binned_count_table.txt.gz"),
		o = f"{REPLICATES_CORRELATION_DIR}{output_prefix}_correlation_table.txt.gz",
		filter_func = '"lambda x, y: abs(x) > 0 and abs(y) > 0"',
		value_func = '"lambda x: math.log10(abs(x))"'
	)
	run_template_bash(
		"biodatatools plot_count_tables_correlation -i {i0} {i1} -o {o} -filter_func {filter_func} -value_func {value_func} -fig_change_kw {figure_change_kw} -fig_save_kw {fig_save_kw}",
		conda_env=conda_env,
		output_func=output_func,
		i0 = rms.file_from_path(f"{REPLICATES_CORRELATION_DIR}{output_prefix}_5pl_binned_count_table.txt.gz"),
		i1 = rms.file_from_path(f"{REPLICATES_CORRELATION_DIR}{output_prefix}_5mn_binned_count_table.txt.gz"),
		o = f"{REPLICATES_CORRELATION_DIR}{output_prefix}_correlation.png",
		filter_func = '"lambda x, y: abs(x) > 0 and abs(y) > 0"',
		value_func = '"lambda x: math.log10(abs(x))"',
		figure_change_kw='\'{"fig_supxlabel_prop":{"text":"$log_{10}$ Reads"}, "fig_supylabel_prop":{"text":"$log_{10}$ Reads"}}\'',
		fig_save_kw='\'{"dpi":300, "bbox_inches":"tight"}\''
	)

	
def RMSTemplate_generate_PROcap_gbratio_20240501(
	rms,
	conda_env: str = "rms_procap_analysis_20240501",	
	project_dir: str="/home/happyuser/PROcap/",
	gbratio_annotation_file: InputFileType = "hg38_qc-gbratio.gtf",
	prefix:str="sample1_rep1"
	):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	ALIGNMENT_DIR = f"{PROJECT_DIR}Alignments/"
	GBRATIO_DIR = f"{PROJECT_DIR}QC/GeneBodyRatio/"
	for dname in [PROJECT_DIR, ALIGNMENT_DIR, GBRATIO_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})

	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	run_template_bash(
		"biodatatools generate_genebody_TSS_ratio_table -label {label} -ibwpl {ibwpl} -ibwmn {ibwmn} -iga {iga} -o {o}",
		conda_env=conda_env,
		output_func=output_func,
		label=prefix,
		ibwpl=rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}_5pl.bw"),
		ibwmn=rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}_5mn.bw"),
		iga=rms.file_from_path(gbratio_annotation_file),
		o=f"{GBRATIO_DIR}{prefix}_gbratio.txt"
	)	
	#input_labels = {f"l{idx}":prefix for idx, prefix in enumerate(prefixes)}
	#input_pl_files = {f"ibwpl{idx}":rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}_5pl.bw") for idx, prefix in enumerate(prefixes)}
	#input_mn_files = {f"ibwmn{idx}":rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}_5mn.bw") for idx, prefix in enumerate(prefixes)}
	#convert_to_input = lambda ts: " ".join(["{" + s + "}" for s in ts])
	
	# run_template_bash(
		# ("biodatatools generate_genebody_TSS_ratio_table" 
		# + " -label " + convert_to_input(input_labels.keys()) 
		# + " -ibwpl " + convert_to_input(input_pl_files.keys()) 
		# + " -ibwmn " + convert_to_input(input_mn_files.keys())
		# + " -iga {iga} -o {o}"),
		# conda_env=conda_env,
		# output_func=output_func,
		# **input_labels,
		# **input_pl_files,
		# **input_mn_files,
		# iga=rms.file_from_path(gbratio_annotation_file),
		# o=f"{GBRATIO_DIR}{output_name}"
	# )
	
def RMSTemplate_generate_PROcap_RNA_length_20240501(
	rms,
	conda_env: str = "rms_procap_analysis_20240501",	
	project_dir: str="/home/happyuser/PROcap/",
	prefix: str="sample1_rep1"
):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	ALIGNMENT_DIR = f"{PROJECT_DIR}Alignments/"
	RNALEN_DIR = f"{PROJECT_DIR}QC/RNALen/"
	for dname in [PROJECT_DIR, ALIGNMENT_DIR, RNALEN_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)

	run_template_bash(
		"biodatatools summarize_PROcap_TSS_RNA_len -i {i0} {i1} -o {o}",
		conda_env=conda_env,
		output_func=output_func,
		i0 = rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}_dpl.bed.bgz"),
		i1 = rms.file_from_path(f"{ALIGNMENT_DIR}{prefix}_dmn.bed.bgz"),
		o = f"{RNALEN_DIR}{prefix}_RNALen_stat.txt",
	)
	
	
def RMSTemplate_merge_PINTS_peaks_20240501(
	rms,
	conda_env: str = "rms_procap_analysis_20240501",	
	project_dir: str="/home/happyuser/PROcap/",
	prefixes: list[str] = ["sample1", "sample2"],
	output_prefix: str = "merged",
	min_overlap_len: int = 1,
	min_overlap_ratio: float = 0.0,
):
	import os
	from commandtemplate.conda import run_template_bash, OUTPUT_run_template_bash
	run_template_bash = rms.register_pipe(run_template_bash, output_func=OUTPUT_run_template_bash)
	if not project_dir.endswith("/"):
		project_dir = project_dir + "/"
	PROJECT_DIR = project_dir
	PEAK_PINTS_DIR = f"{PROJECT_DIR}Peaks/PINTS/"
	PEAK_PINTS_MERGED_DIR = f"{PROJECT_DIR}Peaks/PINTS/Merged/"

	for dname in [PROJECT_DIR, PEAK_PINTS_DIR, PEAK_PINTS_MERGED_DIR]:
		rms.os_operation("makedirs", [dname], {"exist_ok":True})
	prefixes = sorted(prefixes)
	f = lambda **kwargs: [kwargs['o']]
	output_func = rms.register_pipe(f)
	input_files = {f"i{idx}":rms.file_from_path(f"{PEAK_PINTS_DIR}{prefix}_1_divergent_peaks.bed") for idx, prefix in enumerate(prefixes)}
	run_template_bash(
		"biodatatools process_bed_overlapped_regions -i " + " ".join(f"{{i{n}}}" for n in range(len(prefixes))) + " -o {o} -min_overlap_len {min_overlap_len} -min_overlap_ratio {min_overlap_ratio}",
		conda_env=conda_env,
		output_func=output_func,
		**input_files,
		o = f"{PEAK_PINTS_MERGED_DIR}{output_prefix}_divergent_peaks.bed.bgz",
		min_overlap_len=min_overlap_len,
		min_overlap_ratio=min_overlap_ratio,
	)
	input_files = {f"i{idx}":rms.file_from_path(f"{PEAK_PINTS_DIR}{prefix}_1_unidirectional_peaks.bed") for idx, prefix in enumerate(prefixes)}
	run_template_bash(
		"biodatatools process_bed_overlapped_regions -i " + " ".join(f"{{i{n}}}" for n in range(len(prefixes))) + " -o {o} -min_overlap_len {min_overlap_len} -min_overlap_ratio {min_overlap_ratio}",
		conda_env=conda_env,
		output_func=output_func,
		**input_files,
		o = f"{PEAK_PINTS_MERGED_DIR}{output_prefix}_unidirectional_peaks.bed.bgz",
		min_overlap_len=min_overlap_len,
		min_overlap_ratio=min_overlap_ratio,
	)
	run_template_bash(
		"biodatatools filter_bed -i {i} -o {o} -non_overlap_regions {r}",
		conda_env=conda_env,
		output_func=output_func,
		i = rms.file_from_path(f"{PEAK_PINTS_MERGED_DIR}{output_prefix}_unidirectional_peaks.bed.bgz"),
		o = f"{PEAK_PINTS_MERGED_DIR}{output_prefix}_unidirectional-no-divergent_peaks.bed.bgz",
		r = rms.file_from_path(f"{PEAK_PINTS_MERGED_DIR}{output_prefix}_divergent_peaks.bed.bgz"),
	)
	
