gzip SRR15372941.all.txt
gzip ERR2213704.all.txt
# python3 NNpred.py -100 14 tax
# python3 NNpred.py -100 31 pfam
# python3 NNpred.py -100 47 kegg
# python3 NNpred.py -100 79 com
# python3 NNpred.py -100 46 kmc

for i in ../hg_bwa_ec/*fq; do j=$(basename $i); gzip $j; done
for i in *common.txt; do j=$(echo $i | sed s/.common.txt//); cp $i $j.txt; done
## Fiona for error correction ##################################################
PATH=/home/ben0522/toolkit/fiona-0.2.10-Linux-x86_64/bin/:$PATH
fiona -nt 16 -g 4639675 JZLC_04.read1.fq ../hg_bwa_ec/JZLC_04.read1.fq
fiona -nt 16 -g 4639675 JZLC_04.read2.fq ../hg_bwa_ec/JZLC_04.read2.fq

cd /home/ben0522/Data/wgs/JZLC/hg_bwa
for j in 'JZLC_08' 'JZLC_10' 'JZLC_15' 'JZLC_16' 'JZLC_19' 'JZLC_21' 'JZLC_22' 'JZLC_27' 'JZLC_30' 'JZLC_31' 'JZLC_32' 'JZLC_36' 'JZLC_37'
do
  fiona -nt 16 -g 4639675 ${j}.read1.fq ../hg_bwa_ec/${j}.read1.fq;
  fiona -nt 16 -g 4639675 ${j}.read2.fq ../hg_bwa_ec/${j}.read2.fq;
done

# for i in *read1.fq; do j=$(echo $i | sed s/.read1.fq//); fiona -nt 16 -g 4639675 ${j}.read1.fq ../hg_bwa_ec/${j}.read1.fq; fiona -nt 16 -g 4639675 ${j}.read2.fq ../hg_bwa_ec/${j}.read2.fq; done
cd /home/ben0522/Data/wgs/NSCLC/hg_bwa
for i in *fq; do j=$(echo $i | sed s/.fq//); fiona -nt 16 -g 4639675 $i ../hg_bwa_ec/$i; done

cd /home/ben0522/Data/wgs/LC/hg_bwa
for i in *fq; do j=$(echo $i | sed s/.fq//); fiona -nt 60 -g 4639675 $i ../hg_bwa_ec/$i; done

## pure test ###################################################################
# PF05221 PF02741 PF01913 PF02007
# PF01993 PF02240 PF02745 PF02662 PF02505 PF03201 PF04211 PF04210 PF05440 PF02249 PF06253 PF09505 PF09472
# "k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter_smithii"
# "k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Akkermansiaceae;g__Akkermansia;s__Akkermansia_muciniphila"
#  PF12953 PF04034 PF01950 PF04122 PF12837 PF01947 PF02249 PF07541 PF04675 PF02741 PF02240 PF05136 PF13856 PF13276 PF07991 PF05732 PF08542 PF01450 PF09001 PF04886 PF13787 PF14204 PF02241 PF02782 PF11858 PF04208 PF13542 PF12176 PF04010 PF13419 PF02359 PF11068 PF09195 PF02464 PF09991 PF01870 PF02505 PF00742
## naming
# 01993_02240_02745_02662_02505_03201_0421 methanogenesis
# 05221_02741_01913_02007 one_carbon
# 12953_04034_01950_04122_12837_01947_0224 PFS

# reproduce PFS predict response on DS12
python3 TestUsingGivenProfile1.py pfam pfs pfs1 -1 PF12953 PF04034 PF01950 PF04122 PF12837 PF01947 PF02249 PF07541 PF04675 PF02741 PF02240 PF05136 PF13856 PF13276 PF07991 PF05732 PF08542 PF01450 PF09001 PF04886 PF13787 PF14204 PF02241 PF02782 PF11858 PF04208 PF13542 PF12176 PF04010 PF13419 PF02359 PF11068 PF09195 PF02464 PF09991 PF01870 PF02505 PF00742
# test code
python3 TestUsingGivenProfile1.py pfam CB pfs1 -100 PF12953 PF04034 PF01950 PF04122 PF12837 PF01947 PF02249 PF07541 PF04675 PF02741 PF02240 PF05136 PF13856 PF13276 PF07991 PF05732 PF08542 PF01450 PF09001 PF04886 PF13787 PF14204 PF02241 PF02782 PF11858 PF04208 PF13542 PF12176 PF04010 PF13419 PF02359 PF11068 PF09195 PF02464 PF09991 PF01870 PF02505 PF00742
python3 TestUsingGivenProfile1.py pfam CB pfs2 -100 PF12953 PF04034 PF01950 PF04122 PF12837 PF01947 PF02249 PF07541 PF04675 PF02741 PF02240 PF05136 PF13856 PF13276 PF07991 PF05732 PF08542 PF01450 PF09001 PF04886 PF13787 PF14204 PF02241 PF02782 PF11858 PF04208 PF13542 PF12176 PF04010 PF13419 PF02359 PF11068 PF09195 PF02464 PF09991 PF01870 PF02505 PF00742
python3 TestUsingGivenProfile1.py pfam CB pfs1 @100 PF12953 PF04034 PF01950 PF04122 PF12837 PF01947 PF02249 PF07541 PF04675 PF02741 PF02240 PF05136 PF13856 PF13276 PF07991 PF05732 PF08542 PF01450 PF09001 PF04886 PF13787 PF14204 PF02241 PF02782 PF11858 PF04208 PF13542 PF12176 PF04010 PF13419 PF02359 PF11068 PF09195 PF02464 PF09991 PF01870 PF02505 PF00742
python3 TestUsingGivenProfile1.py pfam CB pfs2 @100 PF12953 PF04034 PF01950 PF04122 PF12837 PF01947 PF02249 PF07541 PF04675 PF02741 PF02240 PF05136 PF13856 PF13276 PF07991 PF05732 PF08542 PF01450 PF09001 PF04886 PF13787 PF14204 PF02241 PF02782 PF11858 PF04208 PF13542 PF12176 PF04010 PF13419 PF02359 PF11068 PF09195 PF02464 PF09991 PF01870 PF02505 PF00742

# pfs3 verify
python3 TestUsingGivenProfile1.py pfam CB pfs3 -100 PF12953 PF04034 PF01950 PF04122 PF12837 PF01947 PF02249 PF07541 PF04675 PF02741 PF02240 PF05136 PF13856 PF13276 PF07991 PF05732 PF08542 PF01450 PF09001 PF04886 PF13787 PF14204 PF02241 PF02782 PF11858 PF04208 PF13542 PF12176 PF04010 PF13419 PF02359 PF11068 PF09195 PF02464 PF09991 PF01870 PF02505 PF00742
python3 TestUsingGivenProfile1.py pfam CB pfs3 @100 PF12953 PF04034 PF01950 PF04122 PF12837 PF01947 PF02249 PF07541 PF04675 PF02741 PF02240 PF05136 PF13856 PF13276 PF07991 PF05732 PF08542 PF01450 PF09001 PF04886 PF13787 PF14204 PF02241 PF02782 PF11858 PF04208 PF13542 PF12176 PF04010 PF13419 PF02359 PF11068 PF09195 PF02464 PF09991 PF01870 PF02505 PF00742

# test kmer
for i in "ds2" "ds3" "ds12" "mix" "dds3"
do
python3 TestUsingGivenProfile1.py kmc CB ${i} -100 368404121285272800 1071198102319222352 77813908433010323 92101030321318200 132561842127585341 1045824668049546976 664337655841300785 129858417530129932 30562041644486551 742237302998555048 403264579655749200 389046521065649044 761102070482010181 198635915504996926 309348460840959304 625718481862414349 17237357211788346 122248166577946206 287829317361078646 311255633732041294 724534158377646977 4314370122243715 275726734608212782 287843020423446505 998268572611699968 501107951063285053 1031338988503573760 488992666311784825 351507614151509188 1228818274250464 652432849326687985 401542602153000604 257834747125893440 307683853259964324 346581801908411020 84472338756990243 929149815236608592 276119687823597792 534309776626258945 197030918235963445 240488189680386873 63331487973881882 232287453809152148 984317601898188804 253108951999189776 76920963314991081 253325951895527528 826027895456348480 139405177098967214 460136814016149824 883356503680324240 233405703026797106 62632047765384990 68949428847153384 267799525579805588 549686543164098488 288537580720274360 298646003603024141 69029921955899448 15026007787701204 512222576163816624 307460616980459514 423334264115694664 4915273097001857 478505893772214289 687625751457752320 304573493344656808 280791702686630910 60104031150804819 590978718166495497 4309339302947086 41662509805249591 3756501946925301 287839427271616696 739568964635095484 19661092388007428 166650019569360489 74603328882769145 360557618322279912 851698149268834624 1100949587861644544 547824663441960173 666600078277441956 71959856817904174 370392372810891241 421553632902513840 549425069723647380 58071863452288037 453248904005155441 557634292594797716 613831112438136409 105388408225628460 275797715388613536 85860039407933002 239562982127226498 360557303896073872
python3 TestUsingGivenProfile1.py kmc CB ${i} @100 368404121285272800 1071198102319222352 77813908433010323 92101030321318200 132561842127585341 1045824668049546976 664337655841300785 129858417530129932 30562041644486551 742237302998555048 403264579655749200 389046521065649044 761102070482010181 198635915504996926 309348460840959304 625718481862414349 17237357211788346 122248166577946206 287829317361078646 311255633732041294 724534158377646977 4314370122243715 275726734608212782 287843020423446505 998268572611699968 501107951063285053 1031338988503573760 488992666311784825 351507614151509188 1228818274250464 652432849326687985 401542602153000604 257834747125893440 307683853259964324 346581801908411020 84472338756990243 929149815236608592 276119687823597792 534309776626258945 197030918235963445 240488189680386873 63331487973881882 232287453809152148 984317601898188804 253108951999189776 76920963314991081 253325951895527528 826027895456348480 139405177098967214 460136814016149824 883356503680324240 233405703026797106 62632047765384990 68949428847153384 267799525579805588 549686543164098488 288537580720274360 298646003603024141 69029921955899448 15026007787701204 512222576163816624 307460616980459514 423334264115694664 4915273097001857 478505893772214289 687625751457752320 304573493344656808 280791702686630910 60104031150804819 590978718166495497 4309339302947086 41662509805249591 3756501946925301 287839427271616696 739568964635095484 19661092388007428 166650019569360489 74603328882769145 360557618322279912 851698149268834624 1100949587861644544 547824663441960173 666600078277441956 71959856817904174 370392372810891241 421553632902513840 549425069723647380 58071863452288037 453248904005155441 557634292594797716 613831112438136409 105388408225628460 275797715388613536 85860039407933002 239562982127226498 360557303896073872
done

# checked
for i in "ds2" "ds3" "ds12" "mix" "dds3"
do
  python3 TestUsingGivenProfile1.py pfam 'CB' ${i} -100  PF05221 PF02741 PF01913 PF02007
  python3 TestUsingGivenProfile1.py pfam 'CB' ${i} -100  PF01993 PF02240 PF02745 PF02662 PF02505 PF03201 PF04211 PF04210 PF05440 PF02249 PF06253 PF09505 PF09472
  python3 TestUsingGivenProfile1.py tax  'CB' ${i} -100   "k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter_smithii"
  python3 TestUsingGivenProfile1.py tax  'CB' ${i} -100   "k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Akkermansiaceae;g__Akkermansia;s__Akkermansia_muciniphila"
  python3 TestUsingGivenProfile1.py pfam 'CB' ${i} -100   PF12953 PF04034 PF01950 PF04122 PF12837 PF01947 PF02249 PF07541 PF04675 PF02741 PF02240 PF05136 PF13856 PF13276 PF07991 PF05732 PF08542 PF01450 PF09001 PF04886 PF13787 PF14204 PF02241 PF02782 PF11858 PF04208 PF13542 PF12176 PF04010 PF13419 PF02359 PF11068 PF09195 PF02464 PF09991 PF01870 PF02505 PF00742
done

for i in "ds2" "ds3" "ds12" "mix" "dds3"
do
  python3 TestUsingGivenProfile1.py pfam 'CB' ${i} @100  PF05221 PF02741 PF01913 PF02007
  python3 TestUsingGivenProfile1.py pfam 'CB' ${i} @100  PF01993 PF02240 PF02745 PF02662 PF02505 PF03201 PF04211 PF04210 PF05440 PF02249 PF06253 PF09505 PF09472
  python3 TestUsingGivenProfile1.py tax  'CB' ${i} @100   "k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter_smithii"
  python3 TestUsingGivenProfile1.py tax  'CB' ${i} @100   "k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Akkermansiaceae;g__Akkermansia;s__Akkermansia_muciniphila"
  python3 TestUsingGivenProfile1.py pfam 'CB' ${i} @100   PF12953 PF04034 PF01950 PF04122 PF12837 PF01947 PF02249 PF07541 PF04675 PF02741 PF02240 PF05136 PF13856 PF13276 PF07991 PF05732 PF08542 PF01450 PF09001 PF04886 PF13787 PF14204 PF02241 PF02782 PF11858 PF04208 PF13542 PF12176 PF04010 PF13419 PF02359 PF11068 PF09195 PF02464 PF09991 PF01870 PF02505 PF00742
done

## KMC #########################################################################
PATH=/home/ben0522/toolkit/KMC/bin/:$PATH
K=30
## count kmer
cd /home/ben0522/Data/wgs/JZLC/hg_bwa
for i in *read1.fq.gz; do j=$(echo $i | sed s/.read1.fq.gz//); kmc -k${K} -m258 -t8 -cs100000 @${j}.lst $j /home/ben0522/Data/wgs/JZLC/uproc_bwa; done
mv *kmc* /home/ben0522/Data/wgs/KMC/K${K}/raw
cd /home/ben0522/Data/wgs/NSCLC/hg_bwa
for i in *fq.gz; do j=$(echo $i | sed s/.fq.gz//); kmc -k${K} -m258 -t8 -cs100000 $i $j /home/ben0522/Data/wgs/NSCLC/uproc_bwa; done
mv *kmc* /home/ben0522/Data/wgs/KMC/K${K}/raw
cd /home/ben0522/Data/wgs/LC/hg_bwa
for i in *fq.gz; do j=$(echo $i | sed s/.fq.gz//); kmc -k${K} -m258 -t8 -cs100000 $i $j /home/ben0522/Data/wgs/LC/uproc_bwa; done
mv *kmc* /home/ben0522/Data/wgs/KMC/K${K}/raw

# count common kmer
cd /home/ben0522/Data/wgs/KMC/K${K}/raw
mv ../../*cp .
for i in {1..21}; do kmc_tools -t4 complex s${i}.cp; done
kmc_tools -t4 complex common.cp
mv *cp ../../
rm -f s*kmc*

# get kmer profile
for i in *pre
do
  j=$(echo $i | sed s/.kmc_pre//);
  kmc_tools simple $j result intersect $j.common -cs100000 -ocleft
  kmc_tools -t1 transform $j.common dump ${j}.txt
done

for i in {32..45};   do    bash runKMC.sh $i ; done

# metaphlan4 ###################################################################
source /home/ben0522/toolkit/metaphlan4/metaphlan4/bin/activate

ftr="bwa"
cd /home/ben0522/Data/wgs/JZLC/hg_${ftr}/
for i in *read1.fq.gz
do
  j=$(echo $i | sed s/.read1.fq.gz//);
  metaphlan --nproc 12 -t rel_ab_w_read_stats --input_type fastq --bowtie2out ../uproc_${ftr}/mtp/${j}.bowtie2.bz2 -o ../uproc_${ftr}/otu/${j}.txt ${j}.read1.fq.gz,${j}.read2.fq.gz;
done

cd /home/ben0522/Data/wgs/NSCLC/hg_${ftr}/
for i in *fq.gz
do
  j=$(echo $i | sed s/.fq.gz//);
  metaphlan --nproc 12 -t rel_ab_w_read_stats --input_type fastq --bowtie2out ../uproc_${ftr}/mtp/${j}.bowtie2.bz2 -o ../uproc_${ftr}/otu/${j}.txt ${j}.fq.gz;
done

cd /home/ben0522/Data/wgs/LC/hg_${ftr}/
for i in *fq.gz
do
  j=$(echo $i | sed s/.fq.gz//);
  metaphlan --nproc 12 -t rel_ab_w_read_stats --input_type fastq --bowtie2out ../uproc_${ftr}/mtp/${j}.bowtie2.bz2 -o ../uproc_${ftr}/otu/${j}.txt ${j}.fq.gz;
done

deactivate

# JZLC
ftr="star"
cd /home/ben0522/Data/wgs/JZLC/hg_${ftr}
PATH=/home/ben0522/toolkit/uproc/:$PATH
for j in 'JZLC_08' 'JZLC_21' 'JZLC_22' 'JZLC_32' 'JZLC_36'
do
  uproc-dna -t 12 -o ../uproc_${ftr}/pfam/$j.csv /home/ben0522/toolkit/uproc/pfam /home/ben0522/toolkit/uproc/model $j.read1.fq $j.read2.fq;
  uproc-dna -t 12 -o ../uproc_${ftr}/kegg/$j.csv /home/ben0522/toolkit/uproc/db   /home/ben0522/toolkit/uproc/model $j.read1.fq $j.read2.fq;
done

conda activate /home/ben0522/toolkit/metaphlan
for j in 'JZLC_08' 'JZLC_21' 'JZLC_22' 'JZLC_32' 'JZLC_36'
do
  metaphlan --nproc 12 -t rel_ab_w_read_stats --input_type fastq --bowtie2out ../uproc_${ftr}/mtp/$j.bowtie2.bz2 -o ../uproc_${ftr}/otu/$j.txt $j.read1.fq,$j.read2.fq;
done
conda deactivate


#### remove human contamination ################################################
## map using bwa, output ${sampleID}.sam in the same directory as raw reads
# pair end
for j in 'JZLC_08' 'JZLC_21' 'JZLC_22' 'JZLC_32' 'JZLC_36'
do
  python3 removeReadSAM.py ${j}.read1.fq.gz ${j}.sam ../hg_bwa/${j}.read1.fq
  python3 removeReadSAM.py ${j}.read2.fq.gz ${j}.sam ../hg_bwa/${j}.read2.fq
done

# single end
for j in 'ERR2213660' 'ERR2213665' 'ERR2213666' 'ERR2213675' 'ERR2213677'
do
  python3 removeReadSAM.py ${j}.fastq.gz ${j}.sam ../hg_bwa/${j}.fq
done

## map using star, output mapping/${sampleID}., check the perid!
# pair end
for j in 'JZLC_08' 'JZLC_21' 'JZLC_22' 'JZLC_32' 'JZLC_36'
do
  python3 removeReadSAM.py ${j}.read1.fq.gz mapping/${j}.Aligned.out.sam ../hg_star/${j}.read1.fq
  python3 removeReadSAM.py ${j}.read2.fq.gz mapping/${j}.Aligned.out.sam ../hg_star/${j}.read2.fq
done

# single end
for j in 'ERR2213660' 'ERR2213665' 'ERR2213666' 'ERR2213675' 'ERR2213677'
do
  python3 removeReadSAM.py ${j}.fastq.gz mapping/${j}.Aligned.out.sam ../hg_star/${j}.fq
done

#### deconseq ##################################################################
## input must be unzipped fastq !!!
## execute in the same directory
PATH=/home/ben0522/toolkit/deconseq/:$PATH
deconseq -f JZLC_04.read1.fastq -out_dir JZLC_04
deconseq -f JZLC_04.read2.fastq -out_dir JZLC_04

PATH=/home/ben0522/toolkit/deconseq/:$PATH
for j in 'JZLC_10' 'JZLC_15' 'JZLC_16' 'JZLC_19' 'JZLC_27' 'JZLC_30' 'JZLC_31' 'JZLC_37'
do
  deconseq -f $j.read1.fastq -out_dir $j
  deconseq -f $j.read2.fastq -out_dir $j
  rm -f $j/*clean*
  python3 removeReadFQ.py $j $j.read1.fastq $j.read1.fq
  python3 removeReadFQ.py $j $j.read2.fastq $j.read2.fq
  rm -f $j.read1.fastq $j.read2.fastq
done

#### samtools ##################################################################
PATH=/home/ben0522/test/DRAGoM_v2/lib/:$PATH
samtools view -F 4 JZLC_04.sam > JZLC_04.mapped.sam
samtools view -F -c 4 JZLC_04.sam # count mapping

#### bioNerDS ##################################################################
## this only works on bioNERDS directory
cd /home/ben0522/toolkit/bioNerDS
PATH=/home/ben0522/toolkit/Java/jdk1.7.0_80/bin/:$PATH

java -Xmx1200m -jar bioNerDS.jar --properties properties.conf --text PMC1482722.txt
java -Xmx1200m -jar bioNerDS.jar --properties properties.conf --text PMC1482722.1.txt

for i in outside*txt
do
  j=$(echo $i | sed s/\.txt//);
  java -Xmx1200m -jar bioNerDS.jar --properties properties.conf --text $i;
  mv out.normal.txt $j.tsv
done

for i in outside*tsv; do cat $i | tail -n +2 >> outside.tsv; done

# ssh download #################################################################
scp user@server:/path/to/remotefile.zip /Local/Target/Destination
scp user@server:/path/to/remotefile.zip /Local/Target/Destination

# qiime2 #######################################################################
conda activate /home/ben0522/toolkit/qiime2/qiime2-2020.2

conda deactivate

# AWS ##########################################################################
PATH=/home/ben0522/toolkit/aws/:$PATH

# list file
# https://www.ncbi.nlm.nih.gov/pmc/tools/pmcaws/
aws s3 ls s3://pmc-oa-opendata/ --no-sign-request

# download
aws s3api get-object --bucket DOC-EXAMPLE-BUCKET1 --key dir/my_images.tar.bz2 my_images.tar.bz2

# use sync to download
aws s3 sync s3://pmc-oa-opendata ./pmc-test/ --exclude "*" --include "/oa_comm/txt/all/"

# hmmsearch ####################################################################
PATH=/home/ben0522/toolkit/hmmer/binaries/:$PATH
hmmsearch --cpu 16 --tblout search.tblout ../28/top.hmm Akkermansia_muciniphila.pp.fna

# hmmgraspx ####################################################################
# watch out for directory
cd /home/ben0522/Data/wgs/JZLC/hmmgraspx
/usr/bin/time -v -o JZLC_04.log python batchHMMgraspx.py ../trim_2/JZLC_04.read1.fq.gz  ../trim_2/JZLC_04.read2.fq.gz
/usr/bin/time -v -o JZLC_08.log python batchHMMgraspx.py ../trim_1/trim_2.0/JZLC_08.read1.fq.gz  ../trim_1/trim_2.0/JZLC_08.read2.fq.gz -m 4000000

for i in 'JZLC_08' 'JZLC_10' 'JZLC_15' 'JZLC_16' 'JZLC_19' 'JZLC_21' 'JZLC_22' 'JZLC_27' 'JZLC_30' 'JZLC_31' 'JZLC_32' 'JZLC_36' 'JZLC_37'
do
  /usr/bin/time -v -o $i.log python batchHMMgraspx.py ../trim_2/$i.read1.fq.gz  ../trim_2/$i.read2.fq.gz

done

/usr/bin/time -v -o JZLC_04_4.log perl /home/ben0522/toolkit/HMMGRASPx/Scripts/RunHMMGRASPx.pl \
--hmm=/home/ben0522/Data/pfam/Pfam34-A.hmm \
--seq=/home/ben0522/Data/wgs/JZLC/hmmgraspx/JZLC_04.4000000.fa \
--out=JZLC_04_4 \
--home=/home/ben0522/toolkit/HMMGRASPx \
--index=JZLC_04_4 \
--param=/home/ben0522/toolkit/HMMGRASPx/Settings/param > screen_4


# uproc ########################################################################
PATH=/home/ben0522/toolkit/uproc/:$PATH
uproc-dna -t 16 -o SRR341583.csv /home/ben0522/toolkit/uproc/pfam /home/ben0522/toolkit/uproc/model SRR341583.forward.paired.fq  SRR341583.reverse.paired.fq

# for i in *.reverse.unpaired.fq.gz
for i in *.fastq
do
  j=$(echo $i | sed s/\.fastq//);
  uproc-dna -t 16 -o ../uproc/pfam/$j.csv /home/ben0522/toolkit/uproc/pfam /home/ben0522/toolkit/uproc/model $j.fastq;
  uproc-dna -t 16 -o ../uproc/kegg/$j.csv /home/ben0522/toolkit/uproc/db /home/ben0522/toolkit/uproc/model $j.fastq;
done

for i in *.read1.fq.gz
do
  j=$(echo $i | sed s/\.read1.fq.gz//);
  uproc-dna -t 24 -o ../uproc/pfam/$j.csv /home/ben0522/toolkit/uproc/pfam /home/ben0522/toolkit/uproc/model $i $j.read2.fq.gz;
  uproc-dna -t 24 -o ../uproc/kegg/$j.csv /home/ben0522/toolkit/uproc/db   /home/ben0522/toolkit/uproc/model $i $j.read2.fq.gz;
done

# metaphlan-env ################################################################
# conda create  --prefix /home/ben0522/toolkit/metaphlan
# conda install -c bioconda python=3.7 metaphlan

# conda create  --prefix /home/ben0522/toolkit/humanan
# conda install -c biobakery humann

# conda create  --prefix /home/ben0522/toolkit/picrust
# conda install -c bioconda -c conda-forge picrust
# ### picrust
# conda activate /home/ben0522/toolkit/picrust # works on 16s data
# conda deactivate

### metaphlan
# metaphlan doesn't care pair-end
# https://groups.google.com/g/metaphlan-users/c/FMApX7mQm2k?pli=1
conda activate /home/ben0522/toolkit/metaphlan
## raw read
for i in *.reverse.unpaired.fq.gz
do
  j=$(echo $i | sed s/\.reverse\.unpaired\.fq.gz//);
  metaphlan $j.read1.fq.gz,$j.read2.fq.gz,$j.forward.unpaired.fq.gz,$j.reverse.unpaired.fq.gz --bowtie2out ../mtp/$j.bowtie2.bz2 --nproc 16 --input_type fastq -o ../mtp/relative/$j.txt;
  metaphlan ../mtp/$j.bowtie2.bz2 -t rel_ab_w_read_stats --nproc 16 --input_type bowtie2out -o ../mtp/abs/$j.txt;
done
conda deactivate

for i in *
do
  j=$(echo $i | sed s/.fastq//);
  metaphlan $i --bowtie2out ../mtp/$j.bowtie2.bz2 --nproc 16 --input_type fastq -o ../mtp/$j.txt;
done

## bowtie output
for i in *
do
  j=$(echo $i | sed s/.bowtie2.bz2//);
  metaphlan $i -t rel_ab_w_read_stats --input_type bowtie2out --nproc 16 -o ../$j.txt
done
conda deactivate

### humanan
conda activate /home/ben0522/toolkit/humanan
for i in *fastq
do
  j=$(echo $i | sed s/.fastq//);
  humann --input $i --threads 16 --output ../hum \
  --protein-database /home/ben0522/toolkit/humann_dbs/uniref \
  --metaphlan-options "--bowtie2db /home/ben0522/toolkit/metaphlan/lib/python3.7/site-packages/metaphlan/metaphlan_databases"
done

for i in *.read1.fq.gz
do
  j=$(echo $i | sed s/\.read1\.fq.gz//);
  zcat $j.read1.fq.gz $j.read2.fq.gz > $j.fastq
  humann --input $j.fastq --threads 16 --output ../hum \
  --protein-database /home/ben0522/toolkit/humann_dbs/uniref \
  --metaphlan-options "--bowtie2db /home/ben0522/toolkit/metaphlan/lib/python3.7/site-packages/metaphlan/metaphlan_databases";
  rm -f $j.fastq;
done
conda deactivate

humann_join_tables  -i . -o total_genefamilies.tsv --file_name genefamilies
humann_join_tables  -i . -o total_pathabundance.tsv --file_name pathabundance

humann_renorm_table -i total_genefamilies.tsv -o total_genefamilies_cpm.tsv --units cpm
conda deactivate


# mOTUs_v2 #####################################################################
conda activate mOTUs_v2
# profile
# -c for absolute count
for i in *.reverse.unpaired.fq.gz
do
  j=$(echo $i | sed s/\.reverse\.unpaired\.fq.gz//);
  /home/ben0522/toolkit/mOTUs_v2/motus profile -t 16 -c -f $j.read1.fq.gz -r $j.read2.fq.gz -s $j.forward.unpaired.fq.gz,$j.reverse.unpaired.fq.gz > ../motu/abs/$j.txt
  /home/ben0522/toolkit/mOTUs_v2/motus profile -t 16    -f $j.read1.fq.gz -r $j.read2.fq.gz -s $j.forward.unpaired.fq.gz,$j.reverse.unpaired.fq.gz > ../motu/rel/$j.txt
done
conda deactivate

for i in *fastq
do
  j=$(echo $i | sed s/.fastq//);
  /home/ben0522/toolkit/mOTUs_v2/motus profile -t 16 -c -s $i > ../motu/abs/$j.txt;
  /home/ben0522/toolkit/mOTUs_v2/motus profile -t 16 -s $i > ../motu/rel/$j.txt;
done
conda deactivate

# merge
/home/ben0522/toolkit/mOTUs_v2/motus merge \
-i taxonomy_profile_1.txt,taxonomy_profile_2.txt > all_sample_profiles.txt
conda deactivate

# fastqc #######################################################################
/home/ben0522/toolkit/FastQC/fastqc -t 16 SRR341583.read1.fastq SRR341583.read2.fastq
# fastq-dump ###################################################################
/home/ben0522/toolkit/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump.2.9.0 ERR2213707
# time #########################################################################
/usr/bin/time -v -o log1 /home/ben0522/toolkit/sga/src/SGA/sga preprocess -o reads.pp.fastq --pe-mode 1 ../read1.180.fq ../read2.180.fq

# Trimmomatic ##################################################################
for i in *R1*;
do
  j=$(echo $i | sed s/_R1_001.fastq.gz//)_R2_001.fastq.gz;
  k=$(echo $i | sed s/_L007_R\._001.fastq.gz//);
  java -jar /home/ben0522/toolkit/Trimmomatic-0.38/trimmomatic-0.38.jar \
  PE $i $j \
  $k.read1.fq.gz $k.forward.unpaired.fq.gz \
  $k.read2.fq.gz $k.reverse.unpaired.fq.gz \
  ILLUMINACLIP:adapters.fa:2:30:10 \
  LEADING:30 \
  TRAILING:30 \
  SLIDINGWINDOW:4:30 \
  MINLEN:35 -threads 16;
done

for i in *;
do
  j=$(echo $i | sed s/.fastq.gz//);
  java -jar /home/ben0522/toolkit/Trimmomatic-0.38/trimmomatic-0.38.jar \
  SE $i ../trim/$j.fastq.gz \
  ILLUMINACLIP:/home/ben0522/toolkit/Trimmomatic-0.38/adapters/TruSeq2-SE.fa:2:30:10 \
  CROP:280 \
  LEADING:30 \
  TRAILING:30 \
  SLIDINGWINDOW:4:15 \
  MINLEN:35 -threads 16;
done

# AlienTrimmer #################################################################
java -jar /home/ben0522/toolkit/alienTrimmer/src/AlienTrimmer.jar \
...

# WGsim ########################################################################
# single-end
/home/ben0522/toolkit/wgsim/wgsim -N600000 -1100 -d0 -S9 -e0.001 -r0 \
refG.fna reads.fq /dev/null
# pair-end
/home/ben0522/toolkit/wgsim/wgsim -N2326963 -1100 -2100 -d180 -S7 -e0.01 -r0 \
stagger.fna read1.180.fq read2.180.fq

/home/ben0522/toolkit/wgsim/wgsim -N3300000 -1100 -2100 -d180 -S7 -e0.01 -r0 \
marine.fa marine.read1.fq marine.read2.fq
# blast ########################################################################
PATH=/home/ben0522/toolkit/blast/bin/:$PATH
PATH=/home/ben0522/toolkit/ncbi-magicblast-1.6.0/bin/:$PATH
# mk lib
makeblastdb \
-in /home/ben0522/Data/70M-genome/query/StringGraph2.fastq \
-dbtype nucl
# blast -outfmt 17 for sam
blastn -outfmt 6 \
-db /home/ben0522/Data/70M-genome/query/StringGraph1.fastq \
-query /home/ben0522/Data/70M-genome/query/reads2.fasta \
-out /home/ben0522/Data/70M-genome/BLAST/RVE2.txt

# magicblast
magicblast -db /home/ben0522/Data/RNAseq/hg38/hg38.fa -query JZLC_08.read1.fq.gz -infmt fastq -num_threads 5 > mapped.sam

# SPAdes #######################################################################
# pair-end mode
/home/ben0522/toolkit/SPAdes/bin/spades.py --only-assembler --meta -t 16 \
-1 ../refG/read1.180.fq \
-2 ../refG/read2.180.fq \
-o nec
/home/ben0522/toolkit/SPAdes/bin/spades.py --meta -t 8 \
-1 ../genome/read1.180.fq \
-2 ../genome/read2.180.fq \
-o ec

/home/ben0522/toolkit/SPAdes/bin/spades.py --meta -t 16 -m 1024 \
--12 ../reads/all_read.fq \
-o graph
# --only-assembler, to disable EC
# --meta, to indicate metagenomics

# single-end mode
/home/ben0522/toolkit/SPAdes/bin/spades.py --only-assembler -t 8 \
-s ../refG/SP.01.reads.fq \
-o single.01

# CM ###########################################################################
#cmpress
/home/ben0522/toolkit/infernal/cmpress /home/ben0522/Data/CM/41fam.cm

# CMsearch
/home/ben0522/toolkit/infernal/cmsearch --rfam --nohmmonly --cpu 16 --cut_ga \
--tblout read.tblout \
/home/ben0522/Data/CM/DLF148.cm \
SRR341588.reads.fa

# CMscan, hard-parallel
ssh watson.ittc.ku.edu
cd /home/ben0522/Data/CAMI/toy_medium/sga/nec/cg
/home/ben0522/toolkit/infernal/cmscan --fmt 1 --rfam --nohmmonly --cpu 1 -E 100 \
--tblout 1.tblout \
/home/ben0522/Data/CM/NLF148.cm \
1_og.fq

# CMscan, no parallel
/home/ben0522/toolkit/infernal/cmscan --fmt 1 --rfam --nohmmonly --cpu 4 -E 100 \
--tblout 1.tblout \
/home/ben0522/Data/CM/NLF148.cm \
1_og.fq

# bwa ##########################################################################
PATH=/home/ben0522/toolkit/bwa-master/:$PATH

# index
bwa index contigs.fasta

bwa mem -t 32 RF00001.fa 5SrRNA.fasta > test.sam

# single-end
bwa mem -t 8 -a -T 90 \
contigs.fasta /home/ben0522/Data/CAMI/toy_medium/reads/reads.fa > alldb.sam

# pair-end, interleaved
bwa mem -t 16 -p \
genome.fa ../reads/reads.nameback.fa > genome.pe.sam

# pair-end
bwa mem -t 8 \
NLF148.fna ../read/SRR341618.read1.fa ../read/SRR341618.read2.fa > NLF148.pe.sam

# cd_hit #######################################################################
/home/ben0522/toolkit/cdhit-master/cd-hit-est \
-M 0 -c 0.99 -T 16 \
-o 1.cd.fa -i 1_og.20.fa

# sga ##########################################################################
# single-end
# ec
/home/ben0522/toolkit/sga/src/SGA/sga preprocess -o reads.pp.fastq example.fq
/home/ben0522/toolkit/sga/src/SGA/sga index --no-reverse -t 16 reads.pp.fastq
/home/ben0522/toolkit/sga/src/SGA/sga correct -k 55 --learn -d 256 -t 16 reads.pp.fastq
/home/ben0522/toolkit/sga/src/SGA/sga index -t 16 reads.pp.ec.fa
/home/ben0522/toolkit/sga/src/SGA/sga rmdup -t 16 reads.pp.ec.fa
/home/ben0522/toolkit/sga/src/SGA/sga overlap -e 0 -t 16 -m 29 reads.pp.ec.rmdup.fa
/home/ben0522/toolkit/sga/src/SGA/sga assemble -m 33 --transitive-reduction reads.pp.ec.rmdup.asqg.gz

# nec
/home/ben0522/toolkit/sga/src/SGA/sga preprocess -o reads.pp.fastq example.fq
/home/ben0522/toolkit/sga/src/SGA/sga index -t 8 reads.pp.fastq
/home/ben0522/toolkit/sga/src/SGA/sga rmdup -t 8 reads.pp.fastq
/home/ben0522/toolkit/sga/src/SGA/sga overlap -e 0 -t 8 -m 29 reads.pp.rmdup.fa
/home/ben0522/toolkit/sga/src/SGA/sga assemble -m 33 --transitive-reduction reads.pp.rmdup.asqg.gz

# pair-end, nec
# read in 2 file
/home/ben0522/toolkit/sga/src/SGA/sga preprocess -o reads.pp.fastq --pe-mode 1 \
../read/pe/read.320.1.fq ../read/pe/read.320.2.fq
# read in 1 file
/home/ben0522/toolkit/sga/src/SGA/sga preprocess -o reads.pp.fastq --pe-mode 2 \
../read/pe/read.320.1.fq

/home/ben0522/toolkit/sga/src/SGA/sga index -t 16 reads.pp.fastq
/home/ben0522/toolkit/sga/src/SGA/sga rmdup -t 16 reads.pp.fastq
/home/ben0522/toolkit/sga/src/SGA/sga overlap -t 16 -m 30 reads.pp.rmdup.fa
/home/ben0522/toolkit/sga/src/SGA/sga assemble -m 33 reads.pp.rmdup.asqg.gz

# pair-end ec
# read in 2 file
/home/ben0522/toolkit/sga/src/SGA/sga preprocess -o reads.pp.fastq --pe-mode 1 \
../read/pe/read.320.1.fq ../read/pe/read.320.2.fq
# read in 1 file
/home/ben0522/toolkit/sga/src/SGA/sga preprocess -o reads.pp.fastq --pe-mode 2 \
../read/pe/read.320.1.fq

/home/ben0522/toolkit/sga/src/SGA/sga index --no-reverse -t 8 reads.pp.fastq
/home/ben0522/toolkit/sga/src/SGA/sga correct -k 55 --learn -d 256 -t 8 reads.pp.fastq
/home/ben0522/toolkit/sga/src/SGA/sga index -t 8 reads.pp.ec.fa
/home/ben0522/toolkit/sga/src/SGA/sga rmdup -t 16 reads.pp.ec.fa
/home/ben0522/toolkit/sga/src/SGA/sga overlap -t 16 -m 30 reads.pp.ec.rmdup.fa
/home/ben0522/toolkit/sga/src/SGA/sga assemble -m 33 reads.pp.ec.rmdup.asqg.gz

# complete pipeline
# 1: preprocess
/home/ben0522/toolkit/sga/src/SGA/sga preprocess -o reads.pp.fastq --pe-mode 1 ../read/pe/read.320.1.fq ../read/pe/read.320.2.fq
# 2: index raw
/home/ben0522/toolkit/sga/src/SGA/sga index --no-reverse -t 16 -d 5000000 reads.pp.fastq
# 3: error-correction
/home/ben0522/toolkit/sga/src/SGA/sga correct -k 55 --learn -d 256 -t 16 reads.pp.fastq
# compress file
# 4: index ec
/home/ben0522/toolkit/sga/src/SGA/sga index -t 16 -d 5000000 reads.pp.ec.fa
# 5: filter or rmdup
#/home/ben0522/toolkit/sga/src/SGA/sga filter -d 256 -x 2 -t 16 reads.pp.ec.fa
/home/ben0522/toolkit/sga/src/SGA/sga rmdup -t 16 reads.pp.ec.fa
# compress file
# clean pair-end read
# rename reads or wait after
# 6: index pass reads
/home/ben0522/toolkit/sga/src/SGA/sga index -t 16 reads.pp.ec.rmdup.clean.rename.fa
# 7: overlap
/home/ben0522/toolkit/sga/src/SGA/sga overlap -t 16 -m 33 reads.pp.ec.rmdup.clean.rename.fa
# 8: assemble
/home/ben0522/toolkit/sga/src/SGA/sga assemble -m 37 reads.pp.ec.rmdup.clean.rename.asqg.gz


# RNAstructure and RNAfold #####################################################
# Fold is used to predict the lowest free energy structure
# and a set of suboptimal structures, i.e. low free energy
# structures, using a variety of constraints.
export DATAPATH=/home/ben0522/toolkit/RNAstructure/data_tables/
# get .ct file
/home/ben0522/toolkit/RNAstructure/exe/Fold \
/home/ben0522/Data/seedSS/rsInput/AAGW02014762.1.fa \
/home/ben0522/Data/seedSS/rsOutput/AAGW02014762.1.ct -m 1 \
-c /home/ben0522/Data/seedSS/rsInput/AAGW02014762.1.con

# /home/ben0522/toolkit/RNAstructure/exe/Fold \
# /home/ben0522/Data/seedSS/rsInput/AAGW02014762.1.fa \
# /home/ben0522/Data/seedSS/rsRawOutput/AAGW02014762.1.ct -m 1

/home/ben0522/toolkit/RNAstructure/exe/ct2dot \
/home/ben0522/Data/seedSS/rsOutput/AAGW02014762.1.ct ALL \
/home/ben0522/Data/seedSS/rsOutput/AAGW02014762.1.txt
#
# # get .pfs file
# /home/ben0522/toolkit/RNAstructure/exe/partition \
# /home/ben0522/Data/SS/rf1.fa \
# /home/ben0522/Data/SS/fold1.pfs
#
# # draw SS
# /home/ben0522/toolkit/RNAstructure/exe/draw \
# /home/ben0522/Data/SS/fold1.ct \
# /home/ben0522/Data/SS/fold1.svg \
# --svg -n 1
#
#
# # RNAflod, can't change output dir
/home/ben0522/toolkit/ViennaRNA-2.4.8/exe/bin/RNAfold \
--noPS -C -i /home/ben0522/Data/seedSS/tt.txt >  /home/ben0522/Data/seedSS/tt.fold

/home/ben0522/toolkit/ViennaRNA-2.0.0/exe/bin/RNAfold \
--noPS -C >result_20_-10000.txt < testsequence_20_-10000.txt

# MEGAHIT ######################################################################
/home/ben0522/toolkit/megahit/build/megahit -1 refG/read1.180.fq -2 refG/read2.180.fq -o mega
