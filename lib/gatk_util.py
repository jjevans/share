import os
import re
import subprocess
import sys

#variant calling (gatk) utilities

class Gatk():

	def __init__(self, reference=None, dbsnp=None, roi=None, gatk="gatk", stdout=None, stderr=None):
	#def __init__(self, reference=None, dbsnp=None, known=list(), roi=None, gatk="gatk"):
		#dbsnp file, list of known indels (multiple files ok)
		#region of interest to limit to regions
		#if gatk path not provided assumes CLASSPATH env set for location of jar
		self.gatk = gatk
		self.roi = roi
		self.ref = reference
		#self.known = known
		self.dbsnp = dbsnp
		self.stdout = stdout
		self.stderr = stderr

	def run(self, walker, *args, **kwargs):
		#run gatk command return file path to the output file
		#if gatk path not provided assumes gatk executable is in your PATH
	   #input are the gatk walker name and the arguments for that walker
		# for instance, args run(HaplotypeCaller, R='reference.fasta', I='input.bam', O='output.vcf')
		# special key java_options, val=string of options equates to --java-options
		# gatk --java-options "-Xmx8G" HaplotypeCaller -R reference.fasta -I input.bam -O output.vcf
		#if argument (key) has dashes then assumed to be long args --reference (vs -R)
		#NOTE: if returned output of process is large then use 'stdout' and 'stderr' arguments 
		# and provide file obj or sys.stdout/sys.stderr
		do_print = True
		java_options_default = "-Xmx4G -Xms2G"
		
		cmd = [self.gatk]

		if 'java_options' in kwargs:#set default options for java
			cmd += ["--java-options", kwargs["java_options"]]
			del kwargs['java_options']
		else:
			cmd += ["--java-options", java_options_default]
		
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
		elif self.stderr:
			stderr = self.stderr
		else:
			stderr = subprocess.PIPE

		#tool walker
		cmd.append(walker)

		for k, v in kwargs.items():#add rest of args
			cmd += [self.form_arg(k), str(v)]

		if do_print:
			sys.stderr.write("cmd: "+" ".join(cmd)+"\n\n")
		#sp = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		#sp = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
		sp = subprocess.Popen(cmd, stdout=self.stdout, stderr=self.stderr)
		return sp.communicate()

	def form_arg(self, arg):
		#put two dashes if all lower case otherwise one
		arg = arg.replace("_", "-")
		return "--" + arg if arg.islower() else "-" + arg

class Germline(Gatk):
	#germline variant calling with haplotypecaller
	
	def __init__(self, reference=None, roi=None, gatk="gatk", known=list(), stdout=None, stderr=None, **kwargs):
		#dbsnp file, list of known indels (multiple files ok)
		#region of interest to limit to regions
		#if gatk path not provided assumes PATH env set for executable
		self.gatk = gatk
		self.stdout = stdout
		self.stderr = stderr
		#self.dbsnp = dbsnp
		self.known = known
		self.roi = roi
		self.ref = reference
		self.maxalt = 2

		#haplotypecaller arguments
		self.hc = dict()


		self.hc['max_alternate_alleles'] = kwargs['max_alternate_alleles'] if 'max_alternate_alleles' in kwargs else self.maxalt
		self.hc['output_mode'] = 'EMIT_ALL_SITES'#emit all sites
		
		if 'dbsnp' in kwargs:
			self.dbsnp = kwargs['dbsnp']
			self.hc['dbsnp'] = self.dbsnp


	def haplotypecaller(self, bam, vcf, reference=None, roi=None, **kwargs):
		#call gatk4 haplotype caller by sys call
		#required args are a bam and a reference and 
		# a file path for vcf output
		java_options_default = "-Xms4g -Xmx12g"
		walker = "HaplotypeCaller"
		
		args = dict()
		args.update(self.hc)
		args.update(kwargs)#all args passed in

		args["I"] = bam
		args["O"] = vcf

		
		if reference:#reference seq file provided as args here or in constructor
			self.ref = reference
		elif not self.ref:
				message = "ERROR: haplotypecaller requires a reference sequence fasta"
				raise Exception(message)
		
		args["R"] = self.ref

		#java parameters for running gatk
		if 'java_options' in kwargs:
			args['java_options'] = kwargs['java_options']
		else:
			args['java_options'] = java_options_default

		if roi:#as arg or in constructor
			self.roi = roi
		if self.roi:#roi provided
			args["L"] = self.roi

		return self.run(walker, stdout=self.stdout, stderr=self.stderr, **args)

	def bqsr(self, bam, outfile=None, reference=None, roi=None, **kwargs):
		#BaseRecalibrator on a bam

		args = dict()
		args["R"] = self.ref

		args["I"] = bam

		if outfile is None:
			args["0"] = bam + ".bqsr"
		else:
			args["O"] = outfile


		if reference:#reference seq file provided as args here or in constructor
			self.ref = reference
		elif not self.ref:
				message = "ERROR: haplotypecaller requires a reference sequence fasta"
				raise Exception(message)
		
		args["R"] = self.ref

		#java parameters for running gatk
		if 'java_options' in kwargs:
			args['java_options'] = kwargs['java_options']
		else:
			args['java_options'] = java_options_default

		if roi:#as arg or in constructor
			self.roi = roi
		if self.roi:#roi provided
			args["L"] = self.roi

		return self.run(walker, stdout=self.stdout, stderr=self.stderr, **args)

class Util(Gatk):
	#gatk4 utilities
	def __init__(self, gatk="gatk", stdout=None, stderr=None):
		self.gatk = gatk
		self.stdout = stdout
		self.stderr = stderr

	def tribble_index(self, file):
		#input is bed or vcf file to index
		walker = "IndexFeatureFile"
		return self.run(walker, F=file, stdout=self.stdout, stderr=self.stderr)

	def min_represent(self, vcf, reference, do_trim=True):
		#do left align minimal representation on a vcf of variants and reference fasta
		#uses walker LeftAlignAndTrimVariants
		#default to trim alleles
		''' java -jar GenomeAnalysisTK.jar \
		-T LeftAlignAndTrimVariants \
		-R reference.fasta \
		--variant input.vcf \
		-o output.vcf \
		--splitMultiallelics \
		--dontTrimAlleles
		--keepOriginalAC'''
		walker = "LeftAlignAndTrimVariants"
		return self.run(walker, R=reference, variant=vcf, splitMultiallelics=True, dontTrimAlleles=do_trim, keepOriginalAC=True, stdout=self.stdout, stderr=self.stderr)

