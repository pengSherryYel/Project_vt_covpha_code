#!/usr/bin/bash
##input
name=$1
fq1=$2
fq2=$3
filter_length=${4:-10000} ## this is used for filter scaffolds length above this

##fastpOutput
qcDir="qc/$name"
fastpfq1="$qcDir/$name.1.fastq.gz"
fastpfq2="$qcDir/$name.2.fastq.gz"

##spadesOutput
assembleDir="assemble_spades/$name"
scaffolds="$assembleDir/scaffolds.fasta"
scaffolds_limit_length="$assembleDir/scaffolds.gt${filter_length}.fasta"
scaffolds_limit_length2="$assembleDir/scaffolds.gt1000.fasta"

##mummerOutput
nucmerDir="map2ref_nucmer/$name"
ref="/home/viro/xue.peng/publicData/ncbi/genomes/refseq/GCF/002/116/925/GCF_002116925.1_ASM211692v1/GCF_002116925.1_ASM211692v1_genomic.fna.gz"
nucmerOpt="$nucmerDir/$name.nucmer"

##bowtie2Output

##QUASTOutput
quastOptDir="quast_assemble/$name"


##software
fastp="/home/viro/xue.peng/software_home/fastp/fastp"
spades="/home/viro/xue.peng/software_home/SPAdes-3.15.2-Linux/bin/spades.py"
bowtie2=""
nucmer="/home/viro/xue.peng/software_home/mummer-4.0.0/bin/nucmer"
quast="/home/viro/xue.peng/software_home/quast-5.0.2/quast.py"

## function
echo $name
function mkdirs(){
    dir_name=$1
    clean=${2:-1}
    if [[ -e $dir_name && $clean == 1 ]];then
        rm -rf $dir_name  && mkdir -p $dir_name
    elif [[ -e $dir_name && $clean == 0 ]];then
        echo "dir exist and not clean!!"
    else
        mkdir -p $dir_name
    fi
}

###qc
mkdir -p $qcDir
$fastp -i $fq1 -I $fq2 -o $fastpfq1 -O $fastpfq2 -q 20 -h $qcDir/$name.fastp.html -j $qcDir/$name.fastp.json -z 4 -n 10 -l 60 -5 -3 -W 4 -M 20 -c -g -x

###assemble
mkdirs $assembleDir
$spades --meta -1 $fastpfq1 -2 $fastpfq2 -o $assembleDir -t 20 -m 600
#

### filter length
## vt use both
seqkit seq -m $filter_length $scaffolds >$scaffolds_limit_length
seqkit seq -m 1000 $scaffolds >$scaffolds_limit_length2

### checkv
## vt use 1000
checkv_opt="./checkv/$name/proviruses_virus_all.level_2.fna"
checkv_opt_ori="./checkv/$name/proviruses_virus_all.fna"
###sh /home/viro/xue.peng/script/run_checkv.sh $name $scaffolds_limit_length "checkv"
sh /home/viro/xue.peng/script/run_checkv.sh $name $scaffolds_limit_length2 "checkv"


### vs1
## vt use checkv out ori
sh /home/viro/xue.peng/workplace_2023/vt/src/run_identify_virus_vs1.sh $checkv_opt_ori $name

### pharokka
#sh /home/viro/xue.peng/workplace_2023/vt/src/run_pharokka.sh $name $scaffolds_limit_length pharokka/$name
sh /home/viro/xue.peng/workplace_2023/vt/src/run_pharokka.sh $name $checkv_opt_ori pharokka/$name

### replidec
. /home/viro/xue.peng/software_home/miniconda3/etc/profile.d/conda.sh
conda activate replidec
#Replidec -i $scaffolds_limit_length -p multiSeqEachAsOne -w replidec/$name -c 1e-5 -m 1e-5 -b 1e-5
Replidec -i $checkv_opt_ori -p multiSeqEachAsOne -w replidec/$name -c 1e-5 -m 1e-5 -b 1e-5
conda deactivate

#### host
#. /home/viro/xue.peng/software_home/miniconda3/etc/profile.d/conda.sh
#conda activate iphop_env
#
#host_dir="./host_iphop/$name"
#mkdirs $host_dir
#
###vt use checkv ori
#iphop predict --out_dir $host_dir --db_dir /home/viro/xue.peng/software_home/iphop/iphop_db/Sept_2021_pub --num_threads 20 --fa_file $checkv_opt_ori
#conda deactivate

### bt alignment
#. /home/viro/xue.peng/script/bt2.sh
#btbuild $checkv_opt_ori bt2Index/$name "--threads 3"
#btalign bt2Index/$name $fastpfq1 $fastpfq2 $name "--sensitive-local -q -p 10"
#samtoolsStat ./btalign/$name/$name.bam

### run binning
#. /home/viro/xue.peng/software_home/miniconda3/etc/profile.d/conda.sh
#conda activate cocoNet
#mkdir -p binning_results/$name
#coconet run --fasta $checkv_opt_ori --bam ./btalign/$name/$name.sort.bam --output binning_results/$name --min-ctg-len 1000 --min-prevalence 2 --test-ratio 0.2 -t 10
#coconet run --fasta $checkv_opt_ori --bam ./btalign/$name/$name.sort.bam --output binning_results/$name --min-ctg-len 1000 -t 10 --features "composition" 


#. /home/viro/xue.peng/software_home/miniconda3/etc/profile.d/conda.sh
#conda activate vamb2
#rm -rf ./binning_results_vamb/$name 
#vamb --outdir ./binning_results_vamb/$name --fasta $checkv_opt_ori --bamfiles ./btalign/$name/$name.sort.bam -m 2000



### mmseqs tax
#mmseqs_tmp="mmseqs_tmp_$name"
#mkdir -p mmseqs/$name
#if [ -e $mmseqs_tmp ]; then rm -rf $mmseqs_tmp; fi
##mmseqs easy-taxonomy $checkv_opt /home/viro/xue.peng/software_home/mmseqs_db/swissprot_db/swissprot mmseqs/$name $mmseqs_tmp
##mmseqs easy-taxonomy $scaffolds_limit_length /home/viro/xue.peng/software_home/mmseqs_db/swissprot_db/swissprot mmseqs/$name $mmseqs_tmp
#mmseqs easy-taxonomy $scaffolds_limit_length2 /home/viro/xue.peng/software_home/mmseqs_db/swissprot_db/swissprot mmseqs/$name $mmseqs_tmp

###align2ref
#mkdir -p $nucmerDir
#refStat=`file $ref|grep gzip`
#if [ -n "$refStat" ];then
#    tmpRef="./$name.ref.fasta"
#    zcat $ref >$tmpRef
#    $nucmer -p $nucmerOpt $scaffolds $tmpRef && rm $tmpRef
#else
#    $nucmer -p $nucmerOpt $scaffolds $ref
#fi
#
##delta-filter -l 1000 -q m_harundinacea.delta > m_harundinacea_filter.delta
##show-coords -c -l -L 1000 -r -T m_harundinacea_filter.delta > m_harundinacea_filter_coords.txt

##assess quality
#mkdir -p $quastOptDir
#$quast $scaffolds -r $ref -1 $fastpfq1 -2 $fastpfq2 -o $quastOptDir




