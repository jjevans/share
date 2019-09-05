import sys
import jobs
import os
import os.path
import re
import subprocess
import shutil

#class for picard with gatk4
#requires java 1.8
#jje, ce 09282018
#bi r&d@qdx



#metrics files
#		detailmetrics = self.prefix + ".genotype_concordance_detail_metrics"
#		contingencymetrics = self.prefix + ".genotype_concordance_contingency_metrics"
#		summarymetrics = self.prefix + ".genotype_concordance_summary_metrics"
#		#concordvcf = self.prefix + ".genotype_concordance.vcf"#originally .gz prior to gunzip


class Picard():

	def __init__(self, reference=None, roi=None, picard="/media/bams/workspace/group_bi/util/picard.jar", java_exec=None, stdout=None, stderr=None, xms="8G", xmx="24G", do_print=True, nt=2):
	#def __init__(self, reference=None, dbsnp=None, known=list(), roi=None, gatk="gatk"):
		#dbsnp file, list of known indels (multiple files ok)
		#region of interest to limit to regions
		#if gatk path not provided assumes CLASSPATH env set for location of jar
		self.picard = picard
		self.java = java_exec if java_exec is not None else "java"
		self.roi = roi
		self.ref = reference
		self.stdout = stdout
		self.stderr = stderr
		self.xms = xms
		self.xmx = xmx#java max heap
		self.do_print = do_print
		self.nt = nt
		self.jobs = jobs.Process(nt=self.nt)

	def run(self, tool, *args, **kwargs):
		#run gatk command return file path to the output file
		#if gatk path not provided assumes gatk executable is in your PATH
		#input are the gatk walker name and the arguments for that walker
		# for instance, args run(HaplotypeCaller, R='reference.fasta', I='input.bam', O='output.vcf')
		# special key java_options, val=string of options equates to --java-options
		# gatk --java-options "-Xmx8G" HaplotypeCaller -R reference.fasta -I input.bam -O output.vcf
		#if argument (key) has dashes then assumed to be long args --reference (vs -R)
		#NOTE: if returned output of process is large then use 'stdout' and 'stderr' arguments 
		# and provide file obj or sys.stdout/sys.stderr
		#java_options_default = "-Xmx4G -Xms2G"

		#java heapsizes
		if "xms" in kwargs:#initial
			xms = kwargs["xms"]
			del kwargs["xms"]
		else:
			xms = self.xms

		if "xmx" in kwargs:#max
			xmx = kwargs["xmx"]
			del kwargs["xmx"]
		else:
			xmx = self.xmx

		java_wrap = self.java + " -Xms" + xms + " -Xmx" + xmx + " -jar"

		cmdlst = [java_wrap + " " + self.picard]
		
		if self.do_print:
			sys.stderr.write(" ".join(cmdlst)+"\n")

		
		#output file objects
		if 'stdout' in kwargs:
			stdout = kwargs['stdout']
			del kwargs['stdout']
		elif self.stdout:
			stdout = self.stdout
		else:
			stdout = subprocess.PIPE
		if 'stderr' in kwargs:
			stderr = kwargs['stderr']
			del kwargs['stderr']

		if 'do_wait' in kwargs:
			do_wait = kwargs['do_wait']
			del kwargs['do_wait']
		else:
			do_wait = True

		#tool walker
		cmdlst.append(tool)

		for k, v in kwargs.items():#add rest of args
			cmdlst += [k+"="+v]#I=file

		if self.do_print:
			sys.stderr.write("cmd: "+" ".join(cmdlst)+"\n\n")

		cmd = " ".join(cmdlst)

		
		#sp = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr, shell=True)
		#sp = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
		#return sp.communicate()
		
		if do_wait:
			sp = subprocess.Popen(cmd, stdout=stdout, stderr=stdout, shell=True)
			res = sp.communicate()

			if sp.returncode != 0:
				if res[1] is not None:#error
					for line in res[1].split("\n"):
						sys.stderr.write(line+"\n")

				message = "ERROR: command returned non-zero exit status: " + cmd + "\n"
				raise Exception(message)

		return res[0] if do_wait else cmd

class AlignUtil(Picard):
	#hsmetrics
	def __init__(self, picard="/media/bams/workspace/group_bi/util/picard.jar", java_exec=None, stdout=None, stderr=None, xms="2G", xmx="16G", do_print=True, do_format=True, do_mk_million=True, do_run=True, nt=1):
		#keep histogram and metrics lines in constructor
		self.metrics = None
		self.highqual = None
		self.unfiltered = None
		self.header = None

		self.picard = picard
		self.java = java_exec if java_exec is not None else "java"		
		self.stdout = stdout
		self.stderr = stderr
		self.xms = xms
		self.xmx = xmx
		self.do_print = do_print
		self.do_mk_million = do_mk_million
		self.do_format = do_format#make decimal a percent 0.99 to 99%, round to 3 decimal places
		self.num_digit = 3#number of decimal places to round if self.format
		self.do_run = do_run#do run picard hsmetrics or assume output already exists
		self.nt = nt
		self.intervallist = None
		self.jobs = jobs.Process(nt=self.nt)


	def many_hsmetrics(self, bams, bed, reference=None, do_clean=False, nt=None):
		#run hsmetrics in parallel.  bams is a list of filenames.  bed is a string for a single roi, 
		# outfile_prefix is a common prefix in output filename, nt indicates max number of jobs to run
		actives = list()
		numproc = 0

		if nt is None:
			nt = self.nt

		intervallist = None

		cmds = list()
		for bam in bams:

			##make intervallist file if the first bam, needs bam directory name
			#if self.intervallist is None:
			#	intervallist = os.path.join(os.path.dirname(os.path.abspath(bam)),  os.path.basename(bed) + ".intervallist")
			#	self.intervallist = intervallist

			#outfile = outfile_prefix + "." if outfile_prefix is not None else str()
			outfile = os.path.abspath(bam) + ".hsmetrics"

			cmds.append(self.hsmetrics(bam, bed, outfile=outfile, reference=None, do_clean=do_clean, do_wait=False, intervallist=intervallist))
			intervallist = self.intervallist

		#print "\n".join(cmds)
		statuses = self.jobs.run_many(cmds, nt=nt)

		#reset intervallist in constructor
		self.intervallist = None

		return statuses
		
	def hsmetrics(self, bam, bed=None, outfile=None, reference=None, do_clean=False, do_wait=True, intervallist=None):
		#CollectHsMetrics
		#runs and parses output
		#input is bam, bed intervals, output filename (or named bam + ".intervallist"), 
		# an optional referenence and whether to delete intermediat  files
		#cmds_only will not run the jobs.  instead it returns the command to run
		do_wait = do_wait if do_wait is not None else self.wait
		do_intervallist = True if intervallist is None else False
		 
		if bed is None and intervallist is None:
			message = "ERROR: either bed or intervallist must be provided at picard_util.hsmetrics for bam: " + bam + "\n"
			raise Exception(message)

		tool = "CollectHsMetrics"

		#intervallist filenaming.  roi name.intervallist in run directory with bams
		if intervallist is None:

			if outfile is not None:
				intervallist = os.path.join(os.path.dirname(os.path.abspath(outfile)), os.path.basename(bed) + ".intervallist")#created intervallist  from bed file
			else:
				intervallist = os.path.abspath(bam) + ".intervallist"

		#make intervallist
		if reference is not None:#find a sequence dictionary
			sd = reference
		else:
			sd = bam

		#make intervallist from bed file
		if do_intervallist:
			res = self.bed_to_intervallist(bed, sd, intervallist)
			self.intervallist = intervallist

		kwargs = dict()
		kwargs["do_wait"] = do_wait

		kwargs["I"] = bam
		kwargs["O"] = outfile
		kwargs["BI"] = intervallist
		kwargs["TI"] = intervallist
		kwargs["COVMAX"] = "20000"
		kwargs["MAX_RECORDS_IN_RAM"] = "1000000"

		if reference is not None:
			kwargs["R"] =  reference

		res = self.run(tool, **kwargs)

		if do_clean:
			os.remove(outfile)

		return res

	def qc_report_stat(self, bams, bed, do_clean=False, nt=None):
		#returns table of values desired for qc report for newborn/exome
		#takes bam, roi bed, reference fasta optional
		# creates file outfile, by befault bam + ".hsmetrics", 
		#do_clean will delete this file
		#returns two arrays, one for header column names, the 2nd for column vals(stats)
		#comment values returned (see header)
		colname_highq8 = "PCT_TARGETBASES_8X"
		colname_unfiltered20 = "PCT_TARGETBASES_20X"
		self.output_header = ["Read_Count", "Aligned_Reads", "Duplicate_Rate", "Coverage_8X", "Coverage_20X", "On-Target", "Mean_Target_Coverage", "Target_Pass_Rate", "Number_of_SNV", "Number_of_Indel", "Genotype_Variant_Concordance"]
		tbl = list()

		if nt is None:#number of threads/cores
			nt = self.nt

		#run picard collecthsmetrics
		if self.do_run:#assume or no that file output of hsmetrics already exists, don't run
			#don't define reference since allele dropout not evaluated in qc.
			#self.hsmetrics(bam, bed, outfile=outfile, reference=reference, do_clean=False)
			self.many_hsmetrics(bams, bed, nt=nt)
		
		#parse hsmetrics file
		for bam in bams:
			outfile = bam + ".hsmetrics"

			if os.path.exists(outfile):
				with open(outfile) as metrics_fh:
					metrics = self.collect_stats(metrics_fh)
			else:
				message = "ERROR: hsmetrics file not found.  to run picard set do_run_picard=True.  currently set to: " + str(self.do_run) + ", file: " + str(outfile) + "\n"
				raise Exception(message)

			highq = self.highq
			unfiltered = self.unfiltered

			#get desired columns
			colvals = list()
			colnames = list()

			#totalbases is unused
			#		totalbases = int(metrics[14])
			#totread=col[0], aligned reads=col[1], dupreads=col[2],ontarget=col3, meancov=col4, offtarget=col5
			##total reads #col0
			totalreads = int(metrics[5])

			#check total reads.  die if no reads
			if totalreads < 1:#col 0
				message = "ERROR: total reads found in the alignments metrics is less than 1. \n"
				raise Exception(message)

			elif self.do_mk_million:#report in unit millions of reads (divide by 1,000,000)
				totalreads = float(totalreads)/1000000
		
			if self.do_format:#round to 3 digits
				totalreads = round(totalreads, 3)

			colnames.append(self.header[5])
			colvals.append(totalreads)


			##percent aligned reads. #col1
			colnames.append(self.header[11])

			val = float(metrics[11])
			if self.do_format:
				val = round(val * 100, self.num_digit)
			colvals.append(str(val))


			##duplicate rate #col 2
			colnames.append(self.header[29])

			val = float(metrics[29])
			if self.do_format:
				val = round(val * 100, self.num_digit)
			colvals.append(str(val))


			##!done with bedtools in main script!
			#sum all reads over 8X
			#highq8 = sum(highq[8:])
			#done in coverage:##jje##colnames.append(colname_highq8)
			#done in coverage:##jje##colvals.append(str(highq8))#Coverage_8X

			##!done with bedtools in main script!
			#done in coverage:##jje####done in coverage:##jje##PCT_TARGETBASES_20X =  metrics[38]
			#done in coverage:##jje##colnames.append(self.header[38])
			#done in coverage:##jje##colvals.append(metrics[38])#Coverage_20X

			##!done with bedtools in main script! this is alternative to the metrics[38] immediately above!
			#sum all reads over 20X
			#unfiltered20 = sum(unfiltered[20:])
			#colnames.append(colname_unfiltered20)
			#colvals.append(str(unfiltered20))


			##On-Target percent#col 3
			colnames.append("PCT_ON_BAIT")

			val = 1-float(metrics[19])
			if self.do_format:
				val = round(val * 100, self.num_digit)
			colvals.append(str(val))

		
			##Mean_Target_Coverage #col 4
			colnames.append(self.header[22])
		
			val = float(metrics[22])
			if self.do_format:#only round, not a percent
				val = round(val, self.num_digit)
			colvals.append(str(val))

			##off-target #col 5
			colnames.append(self.header[19])

			val = float(metrics[19])
			if self.do_format:
				val = round(val * 100, self.num_digit)
			colvals.append(str(val))

			#outfh.write( totalreads + '\t' + "PhiX Placeholder" + '\t' + PF_UQ_READS_ALIGNED + '\t' + PCT_EXC_DUPE + "\t" + str(highq8) + '\t' +  PCT_TARGETBASES_20X+ '\t' + 	str(1-PCT_OFF_BAIT) + '\t' + "Mean Target Coverage Placeholder" + '\t' "Target Pass Rate Placeholder" + '\t' +  PCT_EXC_OFF_TARGET	+ '\n')

			tbl.append(colvals)

			#array order
			#totread=col[0], aligned reads=col[1], dupreads=col[2],ontarget=col3, meancov=col4, offtarget=col5

			self.colnames = colnames

		return tbl

	def qc_report_stats(self, bam, bed, outfile=None, reference=None, do_clean=False):
		#returns values desired for qc report for newborn/exome
		#written by witold w2
		#takes bam, roi bed, reference fasta optional
		# creates file outfile, by befault bam + ".hsmetrics", 
		#do_clean will delete this file
		#returns two arrays, one for header column names, the 2nd for column vals(stats)
		#comment values returned (see header)
		colname_highq8 = "PCT_TARGETBASES_8X"
		colname_unfiltered20 = "PCT_TARGETBASES_20X"
		self.output_header = ["Read_Count", "Aligned_Reads", "Duplicate_Rate", "Coverage_8X", "Coverage_20X", "On-Target", "Mean_Target_Coverage", "Target_Pass_Rate", "Number_of_SNV", "Number_of_Indel", "Genotype_Variant_Concordance"]

		if outfile is None:
			outfile = bam + ".hsmetrics"
		
		#run picard collecthsmetrics
		if self.do_run:#assume or no that file output of hsmetrics already exists, don't run
			#don't define reference since allele dropout not evaluated in qc.
			#self.hsmetrics(bam, bed, outfile=outfile, reference=reference, do_clean=False)
			self.hsmetrics(bam, bed, outfile=outfile, do_clean=do_clean)
		
		#parse hsmetrics file
		if os.path.exists(outfile):
			with open(outfile) as metrics_fh:
				metrics = self.collect_stats(metrics_fh)
		else:
			message = "ERROR: hsmetrics file not found.  to run picard set do_run_picard=True.  currently set to: " + str(self.do_run) + ", file: " + str(outfile) + "\n"
			raise Exception(message)

		highq = self.highq
		unfiltered = self.unfiltered

		#get desired columns
		colvals = list()
		colnames = list()

		#totalbases is unused
		#		totalbases = int(metrics[14])
		#totread=col[0], aligned reads=col[1], dupreads=col[2],ontarget=col3, meancov=col4, offtarget=col5
		##total reads #col0
		totalreads = int(metrics[5])

		#check total reads.  die if no reads
		if totalreads < 1:#col 0
			message = "ERROR: total reads found in the alignments metrics is less than 1. \n"
			raise Exception(message)

		elif self.do_mk_million:#report in unit millions of reads (divide by 1,000,000)
			totalreads = float(totalreads)/1000000
		
		if self.do_format:#round to 3 digits
			totalreads = round(totalreads, 3)

		colnames.append(self.header[5])
		colvals.append(totalreads)


		##percent aligned reads. #col1
		colnames.append(self.header[11])

		val = float(metrics[11])
		if self.do_format:
			val = round(val * 100, self.num_digit)
		colvals.append(str(val))


		##duplicate rate #col 2
		colnames.append(self.header[29])

		val = float(metrics[29])
		if self.do_format:
			val = round(val * 100, self.num_digit)
		colvals.append(str(val))


		##!done with bedtools in main script!
		#sum all reads over 8X
		#highq8 = sum(highq[8:])
		#done in coverage:##jje##colnames.append(colname_highq8)
		#done in coverage:##jje##colvals.append(str(highq8))#Coverage_8X

		##!done with bedtools in main script!
		#done in coverage:##jje####done in coverage:##jje##PCT_TARGETBASES_20X =  metrics[38]
		#done in coverage:##jje##colnames.append(self.header[38])
		#done in coverage:##jje##colvals.append(metrics[38])#Coverage_20X

		##!done with bedtools in main script! this is alternative to the metrics[38] immediately above!
		#sum all reads over 20X
		#unfiltered20 = sum(unfiltered[20:])
		#colnames.append(colname_unfiltered20)
		#colvals.append(str(unfiltered20))


		##On-Target percent#col 3
		colnames.append("PCT_ON_BAIT")

		val = 1-float(metrics[19])
		if self.do_format:
			val = round(val * 100, self.num_digit)
		colvals.append(str(val))

		
		##Mean_Target_Coverage #col 4
		colnames.append(self.header[22])
		
		val = float(metrics[22])
		if self.do_format:#only round, not a percent
			val = round(val, self.num_digit)
		colvals.append(str(val))

		##off-target #col 5
		colnames.append(self.header[19])

		val = float(metrics[19])
		if self.do_format:
			val = round(val * 100, self.num_digit)
		colvals.append(str(val))

		#outfh.write( totalreads + '\t' + "PhiX Placeholder" + '\t' + PF_UQ_READS_ALIGNED + '\t' + PCT_EXC_DUPE + "\t" + str(highq8) + '\t' +  PCT_TARGETBASES_20X+ '\t' +  str(1-PCT_OFF_BAIT) + '\t' + "Mean Target Coverage Placeholder" + '\t' "Target Pass Rate Placeholder" + '\t' +  PCT_EXC_OFF_TARGET	+ '\n')
		self.colnames = colnames

		#array order
		#totread=col[0], aligned reads=col[1], dupreads=col[2],ontarget=col3, meancov=col4, offtarget=col5

		return colvals

	def collect_stats(self, metrics_fh, do_clean=False):
		#w2's script, derived from 'hsmetrics_stats_v2.py'
		#gets values of metrics line and histogram
		#returns metrics (stats) line in an array of columns
		#sets to self the high quality and unfiltered histogram 
		# vals in a list (self.hist_highq, self.hist_unfiltered), 
		#and the header columns are in self.header
		metrics = None

		is_statline = False #indicates that line is one with stats
		is_histo = False

		#iterate line by line
		for line in metrics_fh.readlines():
			if is_statline:#line we want, set variables to column values

				#make array of columns of the tab delim statistics line
				metrics = line.rstrip().split("\t")

			#check if this line is header.  set switch to indicate next line is stats line
			if line.startswith("BAIT_SET"):
				is_statline = True

				self.header = line.rstrip().split("\t")#array of column names
			else:
				is_statline = False


		if metrics is None:#error if not found
			message = "ERROR: did not find the statistics line for hsmetrics\n"
			raise Exception(message)


		#histogram of coverage
		highq = list()
		unfiltered = list()

		metrics_fh.seek(0)#back to start until not use readlines

		for line in metrics_fh.readlines():

			if is_histo:#line we want, set variables to column values
				cols=line.rstrip().split("\t")
				#print str(len(cols))
				if len(cols)==3:
					highq.append(int(cols[1]))
					unfiltered.append(int(cols[2]))
				#make array of columns of the tab delim statistics line
	

			#check if this line is header.  set switch to indicate next line is histo line
			if line.startswith("coverage_or_base_quality"):
				is_histo = True

		if len(highq) == 0:#error collecting histogram
			message = "ERROR: did not accumulate the histogram of coverages in hsmetrics\n"
			raise Exception(message)

		self.highq = highq
		self.unfiltered = unfiltered

		return metrics

	#intervallists
	def bed_to_intervallist(self, bed, seqdict, outfile, do_unique=True):
		#run picard BedToIntervalList
		#input bed, a sequence dictionary (vcf, dict, bam) and the output filename
		#optional merge of overlapping region
		tool = "BedToIntervalList"

		if not os.path.exists(bed) or not os.path.exists(seqdict):
			message = "ERROR: bedfile or sequence dictionary file not found at bed_to_intervallist, " + str(bed) + ", " + str(seqdict) + "\n"
			raise Exception(message)

		kwargs = dict()
		kwargs["I"] = bed
		kwargs["SD"] = seqdict
		kwargs["O"] = outfile
		kwargs["UNIQUE"] = "true" if do_unique else "false"
		
		return self.run(tool, **kwargs)

class VcfConcordance():

	#def __init__(self, picard="/media/bams/workspace/group_bi/util/picard.jar", stdout=sys.stdout, stderr=sys.stderr, outprefix="vcfconcordance", do_clean=False, do_run=True, xms="2G", xmx="8G"):
	def __init__(self, picard="/media/bams/workspace/group_bi/util/picard.jar", stdout=None, stderr=None, outdir=None, outprefix="vcfconcordance", do_clean=False, do_run=True, do_print=False, xms="2G", xmx="8G", do_count=True, do_output_vcf=True):
		self.picard = picard
		self.stdout = stdout
		self.stderr = stderr
		self.prefix = outprefix
		self.do_clean = do_clean
		self.do_print = do_print
		self.do_output_vcf = do_output_vcf
		self.do_run = do_run#run picard GenotypeConcordance or assume output files already exist

		self.xms = xms
		self.xmx = xmx

		self.core = Picard(picard=picard, stdout=stdout, stderr=stderr, xms=xms, xmx=xmx, do_print=self.do_print)

		self.colnames = {"detail": list(), "summary": list(), "contingency": list()}#each colname except CALL_SAMPLE, TRUTH_SAMPLE
		self.names = None
		self.file_colnames = {"detail": list(), "summary": list(), "contingency": list()}#each file's header cols
		self.newcols = list()

		self.do_count = do_count#count num snp, indel, total overall
		self.numsnp = None
		self.numindel = None
		self.numtotal = None

		self.call = None
		self.truth = None
		
		self.outdir = outdir

	def genoconcord(self, call_vcf, truth_vcf, do_output_vcf=True, outdir=None, outprefix=None):
		#run and parse picard GenotypeConcordance

		if do_output_vcf:
			self.do_output_vcf = do_output_vcf

		if outprefix is not None:
			self.prefix = outprefix

		#provide the directory for the output, default is call vcf directory
		if outdir is not None:
			self.outdir = outdir
		else:
			self.outdir = os.path.dirname(os.path.abspath(call_vcf))


		if self.do_run:#call picard
			res = self.run_genoconcord(os.path.abspath(call_vcf), os.path.abspath(truth_vcf), do_output_vcf=self.do_output_vcf, outprefix=self.prefix, outdir=self.outdir)

		#parse files
		metrics = self.consume(outprefix=self.prefix, do_output_vcf=self.do_output_vcf)

		#remove files if desired
		if self.do_clean:
			self._cleanup()

		return metrics

	def run_genoconcord(self, call_vcf, truth_vcf, do_output_vcf=True, outprefix=None, outdir=None, ignore_filter=False, ignore_truth_hom=False):
		#execute picard GenotypeConcordance
		#creates files in the directory of the call_vcf
		tool = "GenotypeConcordance"
		
		if outprefix is not None:
			self.prefix = outprefix

		#location of output files, default call vcf dir
		if outdir is not None:
			self.outdir = outdir
		
		if self.outdir is None:#default to call vcf dir
			self.outdir = os.path.dirname(os.path.abspath(call_vcf))

		kwargs = dict()
		#kwargs["stdout"] = subprocess.PIPE if self.stdout is None else self.stdout
		#kwargs["stderr"] = subprocess.PIPE if self.stderr is None else self.stderr
		kwargs["CV"] = os.path.abspath(call_vcf)
		kwargs["TV"] = os.path.abspath(truth_vcf)
		kwargs["O"] = os.path.join(self.outdir, self.prefix)

		kwargs["OUTPUT_VCF"] = "true" if do_output_vcf else "false"
		kwargs["IGNORE_FILTER_STATUS"] = "true" if ignore_filter else "false"
		kwargs["MISSING_HOM"] = "true" if ignore_truth_hom else "false"

		res = None
		if self.do_run:
			res = self.core.run(tool, **kwargs)

		return res

	def consume(self, outprefix=None, outdir=None, do_output_vcf=False):
		#parse genoconcord files
		if do_output_vcf:
			self.do_output_vcf = do_output_vcf

		if outprefix is not None:
			self.prefix = outprefix

		if outdir is not None:
			self.outdir = outdir
		
		if not os.path.exists(self.outdir):
			message = "ERROR: no output directory provided at VcfConcordance.consume\n"
			raise Exception(message)

		#metrics files
		detailmetrics = os.path.join(self.outdir, self.prefix + ".genotype_concordance_detail_metrics")
		contingencymetrics = os.path.join(self.outdir, self.prefix + ".genotype_concordance_contingency_metrics")
		summarymetrics = os.path.join(self.outdir, self.prefix + ".genotype_concordance_summary_metrics")
		#concordvcf = self.prefix + ".genotype_concordance.vcf"#originally .gz prior to gunzip


		#detail metrics
		if os.path.exists(detailmetrics):
			with open(detailmetrics) as handle:
				detail = self.break_detail(handle)#self.detail_metrics(handle, do_pair=True)
				
				self.colnames["detail"] = self.names

				#strip CALL_SAMPLE, TRUTH_SAMPLE
		else:
			message = "ERROR: no detail metrics output file from picard GenotypeConcordance found: " + detailmetrics + "\n"
			raise Exception(message)

		#contingency metrics
		if os.path.exists(contingencymetrics):
			with open(contingencymetrics) as handle:
				contingency = self.break_metrics(handle)
				#self.out_c = self.pair_types(contingency)
				self.colnames["contingency"] = self.names

				if self.do_count:
					#get total snps, indels, and grand total
					#get total snps, indels, and grand total
					vals = contingency["SNP"][0]
					snpcoerce = [ int(val) for val in vals ]
					self.numsnp = sum(snpcoerce)

					vals = contingency["INDEL"][0]
					indelcoerce = [ int(val) for val in vals ]
					self.numindel = sum(indelcoerce)

					self.numtotal = self.numsnp + self.numindel

		else:
			message = "ERROR: no contingency metrics output file from picard GenotypeConcordance found: " + contingencymetrics + "\n"
			raise Exception(message)


		#summary metrics
		if os.path.exists(summarymetrics):
			with open(summarymetrics) as handle:
				summary = self.break_metrics(handle)
				self.colnames["summary"] = self.names

				#!!!here for adding, averaging GENOTYPE_CONCORDANCE from summary file

		else:
			message = "ERROR: no summary metrics output file from picard GenotypeConcordance found: " + summarymetrics + "\n"
			raise Exception(message)


		#parse variant by variant vcf
		'''if self.do_output_vcf:
			concordvcf = self.prefix + ".genotype_concordance.vcf"
			if os.path.exists(concordvcf):#gunzip vcf,  w/ clobber
				sp = subprocess.Popen("gunzip -f "+concordvcf+".gz", shell=True)
				sp.wait()

				#die if gunzip error
				if sp.returncode != 0:#error
					message = "ERROR: gunzip of genotype concordance vcf failed: " + concordvcf + ".gz\n"
					raise Exception(message)

			if os.path.exists(concordvcf):#already gunzipped
				with open(concordvcf) as handle:
					variants = self.concordance_vcf(handle)
			else:
				message = "ERROR: no concordance vcf output file from picard GenotypeConcordance found: " + concordvcf + "\n"
				raise Exception(message)'''

		return summary, detail, contingency


	def break_metrics(self, metric_fh, strip_sample=True):
		#parse 'summary_metrics' file and return dict of colnames and snp and indel vals
		#dict with colname and tuple of (SNP, INDEL) for each val
		#strip_sample remove CALL_SAMPLE, TRUTH_SAMPLE, puts into self always
		is_active = False
		concord = dict()

		for line in metric_fh.readlines():

			if is_active and line.rstrip() != "":
				cols = line.rstrip().split("\t")

				vtype = cols.pop(0)
				self.call = cols.pop(0) if strip_sample else None
				self.truth = cols.pop(0) if strip_sample else None

				if vtype not in concord:
					concord[vtype] = list()

				concord[vtype].append(cols)

				#for i, val in enumerate(self.names):

				#	if val not in concord:#init
				#		concord[val]= [None, None]

				#	if vtype == "SNP":
				#		concord[val][0] = cols[i]
				#	elif vtype == "INDEL":
				#		concord[val][1] = cols[i]
				#	else:
				#		is_active = False#done with getting stats

			elif line.startswith("VARIANT_TYPE"):#header
				cols = line.rstrip().split("\t")

				#to match the way the data rows pop
				cols.pop(0)

				if strip_sample:#remove CALL_SAMPLE, TRUTH_SAMPLE colnames
					del cols[0:2]

				self.names = cols

				is_active = True

		return concord

	def break_detail(self, metric_fh, strip_sample=True):#parse detail metrics into dict by variant type with tuple (SNP, INDEL) for each
		is_active = False
		concord = dict()

		for line in metric_fh.readlines():
			
			if is_active and line.rstrip() != "":
				cols = line.rstrip().split("\t")

				vtype = cols.pop(0)#variant type (snp, indel)
				
				if strip_sample:
					self.call = cols.pop(0)#call name
					self.truth = cols.pop(0)#truth name
				else:
					self.call = cols[0]
					self.truth = cols[1]

				if vtype not in concord:
					concord[vtype] = list()

				concord[vtype].append(cols)
				#for i, val in enumerate(self.names):
				#	if val not in concord:#init
				#		concord[val] = [None, None]

					#if vtype == "SNP":
					#	concord[val][0] = cols[i]
					#elif vtype == "INDEL":
					#	concord[val][1] = cols[i]
					#else:
					#	message = "WARNING: unfamiliar line in detail_metric file.  lines should start with SNP, INDEL or be empty\n"
					#	sys.stderr.write(message)
						
			elif line.startswith("VARIANT_TYPE"):#header
				cols = line.rstrip().split("\t")
				cols.pop(0)#remove type to match data

				#remove CALL_SAMPLE, TRUTH_SAMPLE
				if strip_sample:
					del cols[0:2]

				self.names = cols

				is_active = True

		return concord

	def form(self, dta):
		#make header and all data lines from summary_metric, contingency_metric, detail_metric dictionaries
		#only puts CALL_SAMPLE and TRUTH_SAMPLE in once
		outlines = list()

		outlines.append(self.mrg_head())

		outcols = list()
		for colname in self.file_colnames["summary"]:
			outcols += self.mrg_data(dta[colname])
		
		return

	def mrg_head(self):
		#make header based colnames for each of the details and summary files (self.file_colnames)
		#detail and summary file colnames prepended by SNP_ or INDEL_ 
		# adding each set of columns to array of new column names 
		#ignore CALL_SAMPLE and TRUTH_SAMPLE columns
		#!MUST RUN self.consume() first to make self.file_colnames
		header = ["RUNID", "FLOWCELL", "WELLNUM"]

		#summary_metric file
		cols = self.file_colnames["summary"]
		header += self.header_snp_indel(cols)

		#add contingency_metric variant types
		cols = self.file_colnames["contingency"]
		header += self.header_snp_indel(cols)

		#detail_metric file
		cols = [ col[2] for col in self.file_colnames["detail"]]
		header += self.header_snp_indel(cols)

		return header

	def mrg_data(self, dta):
		#put summary_metric, contingency_metric, and detail_metric data in order
		#skips CALL_SAMPLE and TRUTH_SAMPLE fields
		#must be output format from self.consume (tuples of snp, indel)
		snp = indel = list()
		
		for snpval, indelval in dta:
			snp.append(snpval)
			indel.append(indelval)

		return snp + indel
		
	def header_snp_indel(self, cols):
		#make header from colnames for both snps and indels
		#returns array of colnames prepended with "SNP_" followed by all colnames prepended by "INDEL_"
		snp = indel = list()

		for val in cols:
			if val != "CALL_SAMPLE" and val != "TRUTH_SAMPLE":
				snp.append("SNP_" + val)
				indel.append("INDEL_" + val)

		return snp + indel

	def zip_snp_indel(self, cols):
		#make data from colnames for both snps and indels
		#returns array of values of snps (elem 0) followed by indels (elem 1).  merged to same array
		snp = indel = list()

		for val in cols:
			if val != "CALL_SAMPLE" and val != "TRUTH_SAMPLE":
				snp.append(val[0])
				indel.append(val[1])

		return snp + indel

	def concordance_vcf(self, metric_fh):
		#get summary stats on variant by variant concordance vcf output (gunzipped filehandle)
		header = str()

		for line in metric_fh:
			
			if line.startswith("##"):#keep header
				header += line

			elif line.startswith("#CHROM"):#get sample names from header
				cols = line.rstrip().split("\t")
				
				#sample names compared
				call_smpl = cols[-2]
				self.call = call_smpl
				
				truth_smpl = cols[-1]
				self.truth = truth_smpl
			
				header += line

			else:
				info = self.eval_variant(line.rstrip())


		return

	#!!!!DO THIS OVER ALL VARIANTS.  TOTAL THIS IS FOR SINGLE VARIANT, RETURN STAT AND GET TOTALS FROM THERE!!!#
	def eval_variant(self, vcf_line):
		#determine if concordance between variant no call, correct call and disconconcordant
		#stats AD, GT, genoconcord_TP/FN/FP tuple
		#true negative (no call in either) returns None
		cols = vcf_line.rstrip().split("\t")
		allele_ref = cols[3]
		allele_alt = cols[4]

		concord_call = cols[7].split("=")[1]

		total_variants = 0
		total_snp = 0
		total_indel = 0

		total_empty = 0
		total_nocall = 0
		total_concord = 0
		total_concord_snp = 0
		total_concord_indel_5 = 0
		total_concord_indel_10 = 0

		total_discord = 0
		total_discord_snp = 0
		total_discord_indel_5 = 0
		total_discord_indel_10 = 0


		if concord_call == "EMPTY":
			total_empty += 1

		elif concord_call == "TN":
			total_nocall += 1

		elif concord_call == "TP" or concord_call == "TP,TN":#true positive
			total_concord += 1

		else:
			total_discordant -= 1
			raise Exception("ERROR: nonstandard call type genoconcord\n")

		#eval variant length
		lengths = [len(allele_ref), len(allele_alt)]
		sortlengths = sorted(lengths)
		variant_len = abs(sortlengths[1] - sortlengths[0])

		#is snp or indel
		if variant_len == 0 and len(allele_alt) == 1:
			is_snp = True
		else:
			is_snp = False

	def count_vtypes(self, vcf):
		#count number of indels and snps/delins and return tuple of the num snp, indel, total
		with open(vcf) as handle:
			variants = handle.readlines()

		#skip header
		for skipheader in variants:
			if skipheader.startswith("#CHROM"):
				break

		numtotal = 0
		numsnp = 0
		numindel = 0

		for variant in variants:
			cols = variant.rstrip().split("\t")
			ref = cols[3]
			alt = cols[4]

			lens = [len(ref), len(alt)]
			sortlens = sorted(lens)
			vlen = sortlens[1] - sortlens[0]

			if vlen == 0:
				numsnp += 1
			else:
				numindel += 1

			numtotal += 1

		return numsnp, numindel, numtotal

	def concord_col(self):
		#return array index number of GENOTYPE_CONCORDANCE in summary header (self.colnames["summary"])
		if "GENOTYPE_CONCORDANCE" in self.colnames["summary"]:

			try:
				return self.colnames["summary"].index("GENOTYPE_CONCORDANCE")
			except ValueError:
				pass

		return

	def get_summary_gconcord(self, summary_metrics):
		#return snp and indel tuple of GENOTYPE_CONCORDANCE from the summary metrics
		#takes output of self.break_metrics (array of arrays, header has column of gconcord vals)

		#find column number of GENOTYPE_CONCORDANCE
		colnum = self.concord_col()
		if colnum is None:
			message = "ERROR: cannot find column number of GENOTYPE_CONCORDANCE in summary metrics\n"
			raise Exception(message)

		#get values for snp and indel
		vals = [None, None]
		if "SNP" in summary_metrics and len(summary_metrics["SNP"]) > 0:
			vals[0] = summary_metrics["SNP"][0][colnum]
		else:
			message = "ERROR: SNP value not found in looking for GENOTYPE_CONCORDANCE metric from summary file.\n"
			raise Exception(message)

		if "INDEL" in summary_metrics and len(summary_metrics["INDEL"]) > 0:
			vals[1] = summary_metrics["INDEL"][0][colnum]
		else:
			message = "ERROR: INDEL value not found in looking for GENOTYPE_CONCORDANCE metric from summary file.\n"
			raise Exception(message)

		return vals
		
	def _cleanup(self):
		#remove 4 genoconcord files (details, summary, contingency, and vcf.gz and vcf.gz.tbi
		#returns filenames of removed files
		dfile = self.prefix + ".genotype_concordance_detail_metrics"
		sfile = self.prefix + ".genotype_concordance_summary_metrics"
		cfile = self.prefix + ".genotype_concordance_contingency_metrics"
		vcffile = self.prefix + ".genotype_concordance.vcf"
		vcfgz = vcffile + ".gz"
		vcfgztbi = vcfgz + ".tbi"

		removed = list()#files deleted
		
		res = self._rm_file(dfile)#detail_metrics
		if res is not None:
			removed.append(dfile) 

		res = self._rm_file(sfile)#summary_metrics
		if res is not None:
			removed.append(sfile) 

		res = self._rm_file(cfile)#contingency_metrics
		if res is not None:
			removed.append(cfile) 

		res = self._rm_file(vcffile)#vcf (gunzipped)
		if res is not None:
			removed.append(vcffile) 

		res = self._rm_file(vcfgz)#vcf.gz
		if res is not None:
			removed.append(vcfgz) 

		res = self._rm_file(vcfgztbi)#tabix vcf index
		if res is not None:
			removed.append(vcfgztbi) 

		return removed

	def _rm_file(self, filepath):
		#check for and remove file if exists, return None if not
		if os.path.exists(filepath):
			
			try:
				os.remove(filepath)
				return filepath
			except:
				message = "ERROR: removal of file in cleanup failed for: " + filepath + "\n"
				raise Exception(message)

		return

	'''def bed_to_intervallist(self, bed, seqdict, outfile, do_unique=True):
		#run picard BedToIntervalList
		#input bed, a sequence dictionary (vcf, dict, bam) and the output filename
		#optional merge of overlapping region
		tool = "BedToIntervalList"

		if not os.path.exists(bed) or not os.path.exists(seqdict):
			message = "ERROR: bedfile or sequence dictionary file not found at bed_to_intervallist, " + str(bed) + ", " + str(seqdict) + "\n"
			raise Exception(message)

			kwargs = dict()
			kwargs["I"] = bed
			kwargs["SD"] = seqdict
			kwargs["O"] = outfile
			kwargs["UNIQUE"] = "true" if do_unique else "false"
		
			return self.core.run(tool, **kwargs)'''

