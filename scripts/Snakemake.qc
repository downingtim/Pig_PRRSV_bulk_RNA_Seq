from pathlib import Path

configfile: "config.yaml"
SAMPLEDIR=Path(config["sampledir"])
SAMPLE_PARENT = str(Path(*SAMPLEDIR.parts[1:]))
if SAMPLE_PARENT != ".":
   SAMPLE_PARENT = SAMPLE_PARENT + "/"
else:
   SAMPLE_PARENT = ""
SAMPLEDIR = str(Path(SAMPLEDIR.parts[0]))
print(SAMPLE_PARENT)
print(SAMPLEDIR)

SAMPLES=config["samples"]
MINLEN = config.get("trimmer", {}).get("minlen", 120)
OUTDIR=config["outdir"]
QCTYPE=config["qctype"]
if QCTYPE == "": 
   QCTYPE = "trimmed"
print(QCTYPE)

rule all:
   input:
        fastqc=expand(OUTDIR + "/fastqc/"+SAMPLE_PARENT+"{ID}/{ID}_1_fastqc.zip",ID=SAMPLES),
        fastqcp=expand(OUTDIR + "/fastqc/"+SAMPLE_PARENT+"{ID}/1."+QCTYPE+"_fastqc.zip",ID=SAMPLES),
        QCTYPE1=expand(OUTDIR +"/" +QCTYPE+"/"+SAMPLE_PARENT+"{ID}/{ID}_1.fq.gz",ID=SAMPLES),
        QCTYPE2=expand(OUTDIR +"/" +QCTYPE+"/"+SAMPLE_PARENT+"{ID}/{ID}_2.fq.gz",ID=SAMPLES),
        OGlink1=expand(OUTDIR + "/fastqc/OG/"+SAMPLE_PARENT+"{ID}/1_fastqc.zip",ID=SAMPLES),
        OGlink2=expand(OUTDIR + "/fastqc/OG/"+SAMPLE_PARENT+"{ID}/2_fastqc.zip",ID=SAMPLES),
        QCTYPElink1=expand(OUTDIR + "/fastqc/"+QCTYPE+"/"+SAMPLE_PARENT+"{ID}/1_fastqc.zip",ID=SAMPLES),
        QCTYPElink2=expand(OUTDIR + "/fastqc/"+QCTYPE+"/"+SAMPLE_PARENT+"{ID}/2_fastqc.zip",ID=SAMPLES),
#        stat=expand("stat/"+SAMPLE_PARENT+"{ID}.stat",ID=SAMPLES),
   output:
        multiqc1=OUTDIR+"/fastqc/OG/"+SAMPLE_PARENT+"QC_report_multiqc_report.html",
        multiqc2=OUTDIR+"/fastqc/"+QCTYPE +"/"+SAMPLE_PARENT+"QC_report_multiqc_report.html"
   shell: "cd {OUTDIR}/fastqc/OG/{SAMPLE_PARENT};multiqc .  --interactive  --title \"QC_report\" -f;cd -;cd {OUTDIR}/fastqc/{QCTYPE}/{SAMPLE_PARENT};multiqc .  --interactive  --title \"QC_report\" -f;echo Done"

rule link:
   input:
        fastqc1=OUTDIR+"/fastqc/"+SAMPLE_PARENT+"{ID}/{ID}_1_fastqc.zip",
        fastqc2=OUTDIR+"/fastqc/"+SAMPLE_PARENT+"{ID}/{ID}_2_fastqc.zip",
        f1=OUTDIR+"/fastqc/"+SAMPLE_PARENT+"{ID}/1."+QCTYPE+"_fastqc.zip",
        f2=OUTDIR+"/fastqc/"+SAMPLE_PARENT+"{ID}/2."+QCTYPE+"_fastqc.zip"
   output:
        OGln1=OUTDIR+"/fastqc/OG/"+SAMPLE_PARENT+"{ID}/1_fastqc.zip",
        OGln2=OUTDIR+"/fastqc/OG/"+SAMPLE_PARENT+"{ID}/2_fastqc.zip",
        FPln1=OUTDIR+"/fastqc/"+QCTYPE+"/"+SAMPLE_PARENT+"{ID}/1_fastqc.zip",
        FPln2=OUTDIR+"/fastqc/"+QCTYPE+"/"+SAMPLE_PARENT+"{ID}/2_fastqc.zip",
   shell: "mkdir {OUTDIR}/fastqc/OG/{SAMPLE_PARENT}{wildcards.ID} -p;ln -sf $(readlink -f {input.fastqc1}) {output.OGln1};ln -sf $(readlink -f {input.fastqc2}) {output.OGln2};ln -sf $(readlink -f {input.f1}) {output.FPln1};ln -sf $(readlink -f {input.f2})  {output.FPln2};" 


rule fastqc:
   input:
        f1=SAMPLEDIR + "/"+SAMPLE_PARENT+"{ID}/{ID}_1.fq.gz",
        f2=SAMPLEDIR + "/"+SAMPLE_PARENT+"{ID}/{ID}_2.fq.gz",
   output:
        fastqc1=OUTDIR+"/fastqc/"+SAMPLE_PARENT+"{ID}/{ID}_1_fastqc.zip",
        fastqc2=OUTDIR+"/fastqc/"+SAMPLE_PARENT+"{ID}/{ID}_2_fastqc.zip"
   shell:
       """
        mkdir -p {OUTDIR}/fastqc/{SAMPLE_PARENT}{wildcards.ID} -p
        fastqc {input.f1} -o {OUTDIR}/fastqc/{SAMPLE_PARENT}{wildcards.ID}
        fastqc {input.f2} -o {OUTDIR}/fastqc/{SAMPLE_PARENT}{wildcards.ID}
       """

rule fastqcp:
   input:
        f1=OUTDIR +"/" +QCTYPE+"/"+SAMPLE_PARENT+"{ID}/{ID}_1.fq.gz",
        f2=OUTDIR +"/" +QCTYPE+"/"+SAMPLE_PARENT+"{ID}/{ID}_2.fq.gz"
   output:
        fastqc1=OUTDIR+"/fastqc/"+SAMPLE_PARENT+"{ID}/1."+QCTYPE+"_fastqc.zip",
        fastqc2=OUTDIR+"/fastqc/"+SAMPLE_PARENT+"{ID}/2."+QCTYPE+"_fastqc.zip"
   shell:
       """
        mkdir -p {OUTDIR}/fastqc/{SAMPLE_PARENT}{wildcards.ID}/{QCTYPE} -p
        fastqc {input.f1} -o {OUTDIR}/fastqc/{SAMPLE_PARENT}{wildcards.ID}/{QCTYPE}
        mv {OUTDIR}/fastqc/{SAMPLE_PARENT}{wildcards.ID}/{QCTYPE}/{wildcards.ID}_1_fastqc.zip {output.fastqc1}
        fastqc {input.f2} -o {OUTDIR}/fastqc/{SAMPLE_PARENT}{wildcards.ID}/{QCTYPE}
        mv {OUTDIR}/fastqc/{SAMPLE_PARENT}{wildcards.ID}/{QCTYPE}/{wildcards.ID}_2_fastqc.zip {output.fastqc2}
       """

rule fastp:
   input:
        f1=SAMPLEDIR + "/"+SAMPLE_PARENT+"{ID}/{ID}_1.fq.gz",
        f2=SAMPLEDIR + "/"+SAMPLE_PARENT+"{ID}/{ID}_2.fq.gz",
   output:
        fastqc1=OUTDIR+"/"+QCTYPE+"/"+SAMPLE_PARENT+"{ID}/{ID}_1.fq",
        fastqc2=OUTDIR+"/"+QCTYPE+"/"+SAMPLE_PARENT+"{ID}/{ID}_2.fq",
   params:
        adapters="adapters.fa",
   log:
        stdout=OUTDIR +"/"+QCTYPE+ "/"+SAMPLE_PARENT+"{ID}/log/stdout.txt",
        stderr=OUTDIR +"/"+QCTYPE+ "/"+SAMPLE_PARENT+"{ID}/log/stderr.txt",
   threads: 8
   shell:
       """
       mkdir -p {OUTDIR}/{QCTYPE}/{SAMPLE_PARENT}{wildcards.ID}/log

       if [ {QCTYPE} == "fastp" ];
       then
               if [ "{params.adapters}" != "" ];
               then
                       ADAPTER="--adapter_fasta {params.adapters}" 
               fi
#               fastp --overrepresentation_analysis $ADAPTER --html {OUTDIR}/{QCTYPE}/{SAMPLE_PARENT}{wildcards.ID}/1.{QCTYPE}.html  -i {input.f1} -o {output.fastqc1} > {log.stdout}  2>{log.stderr}
#               fastp --overrepresentation_analysis $ADAPTER --html {OUTDIR}/{QCTYPE}/{SAMPLE_PARENT}{wildcards.ID}/2.{QCTYPE}.html  -i {input.f2} -o {output.fastqc2} >> {log.stdout}  2>>{log.stderr}
               fastp --overrepresentation_analysis $ADAPTER --detect_adapter_for_pe --trim_poly_g \
                  --cut_mean_quality 20 --trim_poly_x --length_required 30 --thread {threads} \
                  --cut_window_size 4 --cut_front --cut_tail  --html {OUTDIR}/{QCTYPE}/{SAMPLE_PARENT}{wildcards.ID}/{QCTYPE}.html  --json {OUTDIR}/{QCTYPE}/{SAMPLE_PARENT}{wildcards.ID}/{QCTYPE}.json \
                  -i {input.f1}  -I {input.f2} -o {output.fastqc1} -O {output.fastqc2} > {log.stdout}  2>{log.stderr}

               cutadapt -a 'A{{10}}' -a 'T{{10}}' -A 'A{{10}}' -A 'T{{10}}' -o {OUTDIR}/{QCTYPE}/{SAMPLE_PARENT}{wildcards.ID}/fastp1c -p {OUTDIR}/{QCTYPE}/{SAMPLE_PARENT}{wildcards.ID}/fastp2c {output.fastqc1} {output.fastqc2}
               rm {output.fastqc1} {output.fastqc2}
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i {OUTDIR}/{QCTYPE}/{SAMPLE_PARENT}{wildcards.ID}/fastp1c -o {output.fastqc1} -v > {log.stdout} 2>{log.stderr}
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i {OUTDIR}/{QCTYPE}/{SAMPLE_PARENT}{wildcards.ID}/fastp2c -o {output.fastqc2} -v >> {log.stdout} 2>>{log.stderr}
       else
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i <(zcat {input.f1}) -o {output.fastqc1} -v > {log.stdout} 2>{log.stderr}
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i <(zcat {input.f2}) -o {output.fastqc2} -v >> {log.stdout} 2>>{log.stderr}
       fi

       """

rule gzip:
   input:
        fastqc1=OUTDIR + "/"+QCTYPE+"/"+SAMPLE_PARENT+"{ID}/{ID}_1.fq",
        fastqc2=OUTDIR + "/"+QCTYPE+"/"+SAMPLE_PARENT+"{ID}/{ID}_2.fq"
   output:
        gz1=OUTDIR + "/"+QCTYPE+"/"+SAMPLE_PARENT+"{ID}/{ID}_1.fq.gz",
        gz2=OUTDIR + "/"+QCTYPE+"/"+SAMPLE_PARENT+"{ID}/{ID}_2.fq.gz"
   shell:
       """
       gzip {input.fastqc1}
       gzip {input.fastqc2}
       """

