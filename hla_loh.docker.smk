##### hlatyping and loh of hla      #####
##### PolySolver + lohhla           #####
configfile: "config.yaml"


import os

WORKDIR=os.getcwd()
#### PolySolver ####
rule loh_type_polysolver:
    input:
        'bam/normal.sort.bam',
    output:
        'loh/polysolver/winners.hla.txt',
        
    params:
        soft=config['lohhla'],
        wkdir=WORKDIR,
    log:
        'logs/loh_type_polysolver.log'
    shell:
        """
docker load -i {params.soft}/images/polysolver.tar.gz 1>>{log} 2>&1
docker run -u $UID -v {params.wkdir}:{params.wkdir} -v {params.soft}:{params.soft} sachet/polysolver:v4 /bin/bash -c \
"cd {params.wkdir}
{params.soft}/scripts/shell_call_hla_type {input} Asian 1 hg19 STDFQ 0 `dirname {output}` 1>{log} 2>&1"
        """


rule loh_run_lohhla:
    input:
        normal='bam/normal.sort.bam',
        tumor='bam/tumor.sort.bam',
        tsv='hla/normal_result.tsv', #optitype result
        polysolver='loh/polysolver/winners.hla.txt',
    output:
        hla_type='loh/hla_optitype',
        solutions='loh/solutions',
        xls='loh/output/sample.10.DNA.HLAlossPrediction_CI.xls',
    params:
        soft=config['lohhla'],
        wkdir=WORKDIR,
    log:
        'logs/loh_run_lohhla.log'
    shell:
        """
{params.soft}/scripts/get_hla.py {input.tsv} `dirname {input.polysolver}` {output.hla_type} 1>{log} 2>&1
ln -sf {params.wkdir}/{input.normal}     loh/normal.bam     2>>{log}
ln -sf {params.wkdir}/{input.tumor}      loh/tumor.bam      2>>{log}
ln -sf {params.wkdir}/{input.normal}.bai loh/normal.bam.bai 2>>{log}
ln -sf {params.wkdir}/{input.tumor}.bai  loh/tumor.bam.bai  2>>{log}
echo -e "Ploidy\ttumorPurity\ttumorPloidy\ntumor\t2\t0.8\t2" 1>{output.solutions} 2>>{log}

if [ ! -s {output.hla_type} ]; then
  mkdir -p `dirname {output.xls}`
  touch {output.hla_type} {output.xls}
  echo 'hla_a, hla_b, hla_c are both HOM, skip loh_run_lohhla!' 1>>{log} 2>&1
else
  
  docker load -i {params.soft}/images/lohhla.tar.gz 1>>{log} 2>&1
  docker run -u $UID -v {params.wkdir}:{params.wkdir} -v {params.soft}:{params.soft} mleventhal/lohhla:latest /bin/bash -c \
"cd {params.wkdir}
export PATH=/bedtools2-master/bin/:$PATH
Rscript {params.soft}/scripts/LOHHLAscript.R \
--patientId sample \
--outputDir {params.wkdir}/loh/output \
--normalBAMfile {params.wkdir}/loh/normal.bam \
--BAMDir {params.wkdir}/loh/ \
--hlaPath {params.wkdir}/{output.hla_type} \
--HLAfastaLoc {params.soft}/data/abc_complete.fasta \
--CopyNumLoc {params.wkdir}/{output.solutions} \
--mappingStep TRUE \
--minCoverageFilter 10 \
--fishingStep TRUE \
--cleanUp FALSE \
--gatkDir /usr/picard.jar \
--novoDir /usr/ \
--HLAexonLoc {params.soft}/data/hla.dat 1>>{log} 2>&1 "
fi
        """

rule loh_get_result:
    input:
        rules.loh_run_lohhla.output.xls
    output:
        'loh/result.tsv'
    params:
        soft=config['lohhla'],
    log:
        'logs/loh_get_result.log'
    shell:
        "python {params.soft}/scripts/get_result.py {input} {output} 1>{log} 2>&1 "


rule loh_all:
    input:
        rules.loh_type_polysolver.output,
        rules.loh_get_result.output,
    output:
        'index/loh_all'
    shell:
        "echo -e 'Complete!\nResults are saved in <{input}>' > {output}"

