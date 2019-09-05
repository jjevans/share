#!/usr/bin/env python
import gatk_util
import sys

#call gatk 4 haplotypecaller
#input bam, reference
#known indels are Mills and 1000G 
# from gatk resource bundle.  dbsnp
# all in cwd
#jje, ce 03022018
#bi r&d@qdx
dbsnp = "dbsnp_135.b37.vcf"
stdout = sys.stdout
stderr = sys.stderr

try:
    bam = sys.argv[1]
    ref = sys.argv[2]
    vcf = sys.argv[3]

except:
    print "usage: gatk_haplocaller.py align.bam reference.fa output.vcf"
    exit()

gatk_obj = gatk_util.Germline(reference=ref, dbsnp=dbsnp, stdout=stdout, stderr=stderr)

res = gatk_obj.haplotypecaller(bam, vcf)

print "\nstdout: " + str(res[0]) + "\nstderr: " + str(res[1]) + "\ndone.\n"

exit()
