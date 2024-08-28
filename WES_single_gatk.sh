#!/usr/bin/bash
#这个流程只适用与单样本，每一个样本只有一对illumina读数产生的PE fastq数据


#工具路径
bwa=/public/home/wulei/BioSoftware/bwa-0.7.17/bwa
samtools=/public/home/wulei/BioSoftware/samtools-1.17/samtools
gatk=/public/home/wulei/miniconda3/envs/gatk/share/gatk4-4.3.0.0-0/gatk

#reference
reference=/public2022/wulei/GRCh38/GRCh38.primary_assembly.genome.fa
GATK_bundle=/public2022/wulei/GRCh38/GATK
G=1000G_phase1.snps.high_confidence.hg38.vcf.gz
dbsnp=dbsnp_146.hg38.vcf.gz
mills=Mills_and_1000G_gold_standard.indels.hg38.vcf.gz


#从窗口中读取参数
fq1=$1 #绝对路径
fq2=$2 #绝对路径
RGID=$3 #read group 单样品就用lane1代替
PL=ILLUMINA #这个值不能乱写，只能是市面上有的，默认illumina
library=$4
sample=$5
PU=$6
outdir=$7 #所有中间文件的输出目录,注意末尾不要加 /

#创建一个输出目录
outdir=${outdir}/${sample}_GATK

#获取fastq文件的前缀名字，样本名字结尾都是.fastq.gz
#fq_file_name=`basename $fq1`   #这里的basename函数用于去除以/结尾的字符串的前缀，比如，basename $workfile ,其中workfile的内容是/home/wulei/project.sh 那么返回值就是 project.sh
#fq_file_name=${fq_file_name%%.1.fastq.gz}

#创建输出目录，这里输入的都是clean data，直接从比对开始进行
if [ ! -d $outdir/bwa ]
then mkdir -p $outdir/bwa
fi
if [ ! -d $outdir/gatk ]
then mkdir -p $outdir/gatk
fi

#使用bwa mem进行数据比对
time $bwa mem -t 8 -M -R "@RG\tID:${RGID}\tPL:${PL}\tLB:${library}\tSM:${sample}\tPU:$PU" ${reference} \
$fq1 $fq2 | $samtools view -Sb - > $outdir/bwa/${sample}.bam && echo "**BWA MEM done**"
#排序
wait
time $samtools sort -@ 4 -m 4G -O BAM -o $outdir/bwa/${sample}.sorted.bam $outdir/bwa/${sample}.bam && echo "**bam file was sorted**"
#给排完序的bam文件建立索引
wait
time $samtools index $outdir/bwa/${sample}.sorted.bam && echo "** ${sample}.sorted.bam index done **"

#标记重复序列
time $gatk MarkDuplicates \
-I $outdir/bwa/${sample}.sorted.bam \
-M $outdir/bwa/${sample}.markdup_metrics.txt \
-O $outdir/bwa/${sample}.sorted.markdup.bam && echo "** ${sample}.sorted.bam MarkDuplicates done **"
wait
#建索引
time $samtools index $outdir/bwa/${sample}.sorted.markdup.bam && echo "**${sample}.sorted.markdup.bam index done"
wait
#外显子组数据不需要进行SplitNCigarReads
#BQSR
time $gatk BaseRecalibrator \
-R $reference \
-I $outdir/bwa/${sample}.sorted.markdup.bam \
--known-sites ${GATK_bundle}/${G} \
--known-sites ${GATK_bundle}/${mills} \
--known-sites ${GATK_bundle}/${dbsnp} \
-O $outdir/bwa/${sample}.recal_data.table 
wait
time $gatk ApplyBQSR \
--bqsr-recal-file $outdir/bwa/${sample}.recal_data.table \
-R $reference \
-I $outdir/bwa/${sample}.sorted.markdup.bam \
-O $outdir/bwa/${sample}.sorted.markdup.BQSR.bam 
wait
echo "**BQSR DONE**"
#建索引
time $samtools index $outdir/bwa/${sample}.sorted.markdup.BQSR.bam && echo "**${sample}.sorted.markdup.BQSR.bam index done**"
wait
#chrom=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM )
#for i in ${chrom[@]};
#do
#time $gatk HaplotypeCaller \
#-R $reference \
#-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
#-L $i
#-O $outdir/gatk/${sample}.HC.${i}.vcf.gz && echo "**HaplotypeCaller done**"
#done && wait

index=0
for i in `ls /public2022/wulei/GRCh38/interval/new_interval/` ;
do
echo $i
let index+=1 
time $gatk HaplotypeCaller \
-R ${reference} \
-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
-L /public2022/wulei/GRCh38/interval/new_interval/${i}/scattered.interval_list \
-O ${outdir}/gatk/${sample}.HC.${index}.vcf.gz &
done && wait

#merge vcfs
#用生成的vcf文件名拼接输入字段
vcfs=""
for z in `ls $outdir/gatk/*.vcf.gz` ;
do
vcfs="-I $z $vcfs"
done && echo ${vcfs}
time $gatk MergeVcfs ${vcfs} -O ${outdir}/gatk/${sample}.HC.merge.vcf.gz && echo "** MergeVcfs done **"
wait
#解压
gunzip ${outdir}/gatk/${sample}.HC.merge.vcf.gz
#Hard Filter
mkdir ${outdir}/gatk/filter
# 使用SelectVariants，选出SNP
time $gatk SelectVariants \
-select-type SNP \
-V ${outdir}/gatk/${sample}.HC.merge.vcf \
-O ${outdir}/gatk/filter/${sample}.HC.snp.vcf
# 使用SelectVariants，选出Indel
time $gatk SelectVariants \
-select-type INDEL \
-V ${outdir}/gatk/${sample}.HC.merge.vcf \
-O  ${outdir}/gatk/filter/${sample}.HC.indel.vcf
wait
# 为SNP作硬过滤
time $gatk VariantFiltration \
-V ${outdir}/gatk/filter/${sample}.HC.snp.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O ${outdir}/gatk/filter/${sample}.HC.snp.filter.vcf
wait
# 为Indel作过滤
# 为Indel作过滤
time $gatk VariantFiltration \
--variant ${outdir}/gatk/filter/${sample}.HC.indel.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
--output ${outdir}/gatk/filter/${sample}.HC.indel.filter.vcf
wait
#time $gatk VariantFiltration --variant ${outdir}/gatk/filter/${sample}.HC.indel.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" --output ${outdir}/gatk/filter/${sample}.HC.indel.filter.vcf
wait
# 重新合并过滤后的SNP和Indel
time $gatk MergeVcfs \
-I ${outdir}/gatk/filter/${sample}.HC.snp.filter.vcf \
-I ${outdir}/gatk/filter/${sample}.HC.indel.filter.vcf \
-O ${outdir}/gatk/filter/${sample}.HC.filter.vcf
wait
#筛选出通过的行，写入PASS
awk '$7 == "PASS" { print }' ${outdir}/gatk/filter/${sample}.HC.filter.vcf > ${outdir}/gatk/filter/${sample}.HC.filter.PASS.vcf





#VQSR, 由于评价SNP和Indel质量高低的标准是不同的，因此，需要分SNP和Indel这两种不同的模式，分别进行质控

#先是SNP mode
#time $gatk VariantRecalibrator \
#-R $reference/Homo_sapiens_assembly38.fasta \
#-V $outdir/gatk/${sample}.HC.vcf.gz \
#-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf \
#-resource:omini,known=false,training=true,truth=false,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf \
#-resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf \
#-resource:dbsnp,known=true,training=false,truth=false,prior=5.0 $GATK_bundle/dbsnp_146.hg38.vcf \
#-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
#-mode SNP \
#-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
#-rscriptFile $outdir/gatk/${sample}.HC.snps.plots.R \
#--tranches-file $outdir/gatk/${sample}.HC.snps.tranches \
#-O $outdir/gatk/${sample}.HC.snps.recal && \

#time $gatk ApplyVQSR \
#-R $reference/Homo_sapiens_assembly38.fasta \
#-V $outdir/gatk/${sample}.HC.vcf.gz \
#--ts_filter_level 99.0 \
#--tranches-file $outdir/gatk/${sample}.HC.snps.tranches \
#-recalFile $outdir/gatk/${sample}.HC.snps.recal \
#-mode SNP \
#-O $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz && echo "** SNPs VQSR done **"

## 然后是Indel mode
#time $gatk VariantRecalibrator \
#$reference/Homo_sapiens_assembly38.fasta \
#-input $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz \
#resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
#-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
#-mode INDEL \
#--max-gaussians 6 \
#-rscriptFile $outdir/gatk/${sample}.HC.snps.indels.plots.R \
#--tranches-file $outdir/gatk/${sample}.HC.snps.indels.tranches \
#-O $outdir/gatk/${sample}.HC.snps.indels.recal && \
#time $gatk ApplyVQSR \
#-R $reference/Homo_sapiens_assembly38.fasta \
#-input $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz \
#--ts_filter_level 99.0 \
#--tranches-file $outdir/gatk/${sample}.HC.snps.indels.tranches \
#-recalFile $outdir/gatk/${sample}.HC.snps.indels.recal \
#-mode INDEL \
#-O $outdir/gatk/${sample}.HC.VQSR.vcf.gz && echo "** SNPs and Indels VQSR (${sample}.HC.VQSR.vcf.gz finish) done **

