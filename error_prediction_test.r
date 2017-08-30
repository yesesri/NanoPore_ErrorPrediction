#setwd("C:/Users/cherukuri/Desktop/test/")####settheworkingdirectory####
options(warn=-1)
args = commandArgs(trailingOnly=TRUE)
require('stringr')
require('Biostrings')
require("glmnet")
library("Biostrings")
library("stringr")
library("glmnet")
print("Required libraries loaded")
y<-NULL
Z<-NULL
matrix<-NULL
p<-NULL
read_length<-NULL
sequence_datafile <- readDNAStringSet(args[1],"fasta")####inputdatasetinfastaformat####
print("Input dataset loaded")
for(f in seq(1,2,1)) #length(sequence_datafile)
{
each_read<-sequence_datafile[f]

for(ls in seq(1,100-21,1)) #nchar(each_read)
{
out<-capture.output(subseq(each_read,start=ls,end=ls+9))
y<-rbind(y,out)
pos<-ls
p<-rbind(p,pos+10)
r<-nchar(each_read)
read_length<-data.frame(rbind(read_length,r))
seq_mer<-data.frame(do.call('rbind',strsplit(as.character(y[,3]),' ',fixed=TRUE)))
left_mer<-data.frame(as.character(seq_mer$X6))
read_name<-data.frame(as.character(seq_mer[,ncol(seq_mer)-1]))
}
for(rs in seq(12,100-10,1)) #nchar(each_read)
{
out<-capture.output(subseq(each_read,start=rs,end=rs+9))
Z<-rbind(Z,out)
seq_rightmer<-data.frame(do.call('rbind',strsplit(as.character(y[,3]),'',fixed=TRUE)))
right_mer<-data.frame(as.character(seq_mer$X6))

}
mat_kmer<-cbind(left_mer,right_mer)
mat_kmer$x<-paste(mat_kmer$as.character.seq_mer.X6.,mat_kmer$as.character.seq_mer.X6.)
mat<-cbind(mat_kmer,read_name,p,read_length)
matrix<-rbind(matrix,mat)
}
print("calculating distance")
matrix$start_distance <- matrix[,5] - 1
print("start_distance calculation done")
matrix$end_distance<- matrix[,6] - matrix[,5]
print("end distance calculation done")
matrix$distance<-pmin(matrix$start_distance,matrix$end_distance)
print("extracted the shortest of the distances")
matrix<-matrix[,-(7:8)]
print("calculating ratio")
matrix$ratio <- round((matrix[,6]/matrix[,7]),2)
print("matrix generated")
names(matrix)[1]<- paste("left_10mer")
names(matrix)[2]<- paste("right_10mer")
names(matrix)[3]<- paste("10_mer")
names(matrix)[4]<- paste("read_name")
names(matrix)[5]<- paste("base_position")
names(matrix)[6]<- paste("read_length")
names(matrix)[7]<- paste("distance")
names(matrix)[8]<- paste("ratio")
print("frequency calculation started")
#write.table(matrix,file="kmer_matrix.txt",row.names=T,quote=F,sep="\t")#thematrixofextractedkmers#
A_final<-NULL
T_final<-NULL
G_final<-NULL
C_final<-NULL
AA_final<-NULL
TT_final<-NULL
GG_final<-NULL
CC_final<-NULL
AT_final<-NULL
AG_final<-NULL
AC_final<-NULL
TG_final<-NULL
TC_final<-NULL
TA_final<-NULL
GA_final<-NULL
GT_final<-NULL
GC_final<-NULL
CA_final<-NULL
CT_final<-NULL
CG_final<-NULL
ATT_final<-NULL
CTT_final<-NULL
GTT_final<-NULL
TTT_final<-NULL
ATG_final<-NULL
TGT_final<-NULL
GCT_final<-NULL
GGT_final<-NULL
CCT_final<-NULL
ACT_final<-NULL
TCT_final<-NULL
TAT_final<-NULL
TGG_final<-NULL
CAA_final<-NULL
AAT_final<-NULL
CAT_final<-NULL
GAA_final<-NULL
GAT_final<-NULL
AAA_final<-NULL
CGT_final<-NULL
TAA_final<-NULL
ATC_final<-NULL
CTC_final<-NULL
GTC_final<-NULL
TTC_final<-NULL
TGC_final<-NULL
GCC_final<-NULL
GGC_final<-NULL
CCC_final<-NULL
ACC_final<-NULL
TCC_final<-NULL
TAC_final<-NULL
CAG_final<-NULL
AAC_final<-NULL
CAC_final<-NULL
GAG_final<-NULL
GAC_final<-NULL
AAG_final<-NULL
CGC_final<-NULL
TAG_final<-NULL
ATA_final<-NULL
CTA_final<-NULL
GTA_final<-NULL
GCA_final<-NULL
GGA_final<-NULL
CCA_final<-NULL
ACA_final<-NULL
TCA_final<-NULL
CGA_final<-NULL
TGA_final<-NULL
CTG_final<-NULL
GTG_final<-NULL
GCG_final<-NULL
GGG_final<-NULL
CCG_final<-NULL
ACG_final<-NULL
TCG_final<-NULL
CGG_final<-NULL
TTA_final<-NULL
AGT_final<-NULL
AGA_final<-NULL
TTG_final<-NULL
AGC_final<-NULL
AGG_final<-NULL
for(d in 1:nrow(matrix))
{
kmer_seq<-str_c(matrix$`10_mer`[d],collapse="")
A<-str_count(kmer_seq,pattern='A')
T<-str_count(kmer_seq,pattern='T')
G<-str_count(kmer_seq,pattern='G')
C<-str_count(kmer_seq,pattern='c')
AA<-str_count(kmer_seq,pattern='AA')
TT<-str_count(kmer_seq,pattern='TT')
GG<-str_count(kmer_seq,pattern='GG')
CC<-str_count(kmer_seq,pattern='CC')
AT<-str_count(kmer_seq,pattern='AT')
AG<-str_count(kmer_seq,pattern='AG')
AC<-str_count(kmer_seq,pattern='AC')
TG<-str_count(kmer_seq,pattern='TG')
TC<-str_count(kmer_seq,pattern='TC')
TA<-str_count(kmer_seq,pattern='TA')
GA<-str_count(kmer_seq,pattern='GA')
GT<-str_count(kmer_seq,pattern='GT')
GC<-str_count(kmer_seq,pattern='GC')
CA<-str_count(kmer_seq,pattern='CA')
CT<-str_count(kmer_seq,pattern='CT')
CG<-str_count(kmer_seq,pattern='CG')
ATT<-str_count(kmer_seq,pattern="ATT")
CTT<-str_count(kmer_seq,pattern="CTT")
GTT<-str_count(kmer_seq,pattern="GTT")
TTT<-str_count(kmer_seq,pattern="TTT")
ATG<-str_count(kmer_seq,pattern="ATG")
TGT<-str_count(kmer_seq,pattern="TGT")
GCT<-str_count(kmer_seq,pattern="GCT")
GGT<-str_count(kmer_seq,pattern="GGT")
CCT<-str_count(kmer_seq,pattern="CCT")
ACT<-str_count(kmer_seq,pattern="ACT")
TCT<-str_count(kmer_seq,pattern="TCT")
TAT<-str_count(kmer_seq,pattern="TAT")
TGG<-str_count(kmer_seq,pattern="TGG")
CAA<-str_count(kmer_seq,pattern="CAA")
AAT<-str_count(kmer_seq,pattern="AAT")
CAT<-str_count(kmer_seq,pattern="CAT")
GAA<-str_count(kmer_seq,pattern="GAA")
GAT<-str_count(kmer_seq,pattern="GAT")
AAA<-str_count(kmer_seq,pattern="AAA")
CGT<-str_count(kmer_seq,pattern="CGT")
TAA<-str_count(kmer_seq,pattern="TAA")
ATC<-str_count(kmer_seq,pattern="ATC")
CTC<-str_count(kmer_seq,pattern="CTC")
GTC<-str_count(kmer_seq,pattern="GTC")
TTC<-str_count(kmer_seq,pattern="TTC")
TGC<-str_count(kmer_seq,pattern="TGC")
GCC<-str_count(kmer_seq,pattern="GCC")
GGC<-str_count(kmer_seq,pattern="GGC")
CCC<-str_count(kmer_seq,pattern="CCC")
ACC<-str_count(kmer_seq,pattern="ACC")
TCC<-str_count(kmer_seq,pattern="TCC")
TAC<-str_count(kmer_seq,pattern="TAC")
CAG<-str_count(kmer_seq,pattern="CAG")
AAC<-str_count(kmer_seq,pattern="AAC")
CAC<-str_count(kmer_seq,pattern="CAC")
GAG<-str_count(kmer_seq,pattern="GAG")
GAC<-str_count(kmer_seq,pattern="GAC")
AAG<-str_count(kmer_seq,pattern="AAG")
CGC<-str_count(kmer_seq,pattern="CGC")
TAG<-str_count(kmer_seq,pattern="TAG")
ATA<-str_count(kmer_seq,pattern="ATA")
CTA<-str_count(kmer_seq,pattern="CTA")
GTA<-str_count(kmer_seq,pattern="GTA")
GCA<-str_count(kmer_seq,pattern="GCA")
GGA<-str_count(kmer_seq,pattern="GGA")
CCA<-str_count(kmer_seq,pattern="CCA")
ACA<-str_count(kmer_seq,pattern="ACA")
TCA<-str_count(kmer_seq,pattern="TCA")
CGA<-str_count(kmer_seq,pattern="CGA")
TGA<-str_count(kmer_seq,pattern="TGA")
CTG<-str_count(kmer_seq,pattern="CTG")
GTG<-str_count(kmer_seq,pattern="GTG")
GCG<-str_count(kmer_seq,pattern="GCG")
GGG<-str_count(kmer_seq,pattern="GGG")
CCG<-str_count(kmer_seq,pattern="CCG")
ACG<-str_count(kmer_seq,pattern="ACG")
TCG<-str_count(kmer_seq,pattern="TCG")
CGG<-str_count(kmer_seq,pattern="CGG")
TTA<-str_count(kmer_seq,pattern="TTA")
AGT<-str_count(kmer_seq,pattern="AGT")
AGA<-str_count(kmer_seq,pattern="AGA")
TTG<-str_count(kmer_seq,pattern="TTG")
AGC<-str_count(kmer_seq,pattern="AGC")
AGG<-str_count(kmer_seq,pattern="AGG")
A_final<-rbind(A_final,A)
T_final<-rbind(T_final,T)
G_final<-rbind(G_final,G)
C_final<-rbind(C_final,C)
AA_final<-rbind(AA_final,AA)
TT_final<-rbind(TT_final,TT)
GG_final<-rbind(GG_final,GG)
CC_final<-rbind(CC_final,CC)
AT_final<-rbind(AT_final,AT)
AG_final<-rbind(AG_final,AG)
AC_final<-rbind(AC_final,AC)
TG_final<-rbind(TG_final,TG)
TC_final<-rbind(TC_final,TC)
TA_final<-rbind(TA_final,TA)
GA_final<-rbind(GA_final,GA)
GT_final<-rbind(GT_final,GT)
GC_final<-rbind(GC_final,GC)
CA_final<-rbind(CA_final,CA)
CT_final<-rbind(CT_final,CT)
CG_final<-rbind(CG_final,CG)
ATT_final<-rbind(ATT_final,ATT)
CTT_final<-rbind(CTT_final,CTT)
GTT_final<-rbind(GTT_final,GTT)
TTT_final<-rbind(TTT_final,TTT)
ATG_final<-rbind(ATG_final,ATG)
TGT_final<-rbind(TGT_final,TGT)
GCT_final<-rbind(GCT_final,GCT)
GGT_final<-rbind(GGT_final,GGT)
CCT_final<-rbind(CCT_final,CCT)
ACT_final<-rbind(ACT_final,ACT)
TCT_final<-rbind(TCT_final,TCT)
TAT_final<-rbind(TAT_final,TAT)
TGG_final<-rbind(TGG_final,TGG)
CAA_final<-rbind(CAA_final,CAA)
AAT_final<-rbind(AAT_final,AAT)
CAT_final<-rbind(CAT_final,CAT)
GAA_final<-rbind(GAA_final,GAA)
GAT_final<-rbind(GAT_final,GAT)
AAA_final<-rbind(AAA_final,AAA)
CGT_final<-rbind(CGT_final,CGT)
TAA_final<-rbind(TAA_final,TAA)
ATC_final<-rbind(ATC_final,ATC)
CTC_final<-rbind(CTC_final,CTC)
GTC_final<-rbind(GTC_final,GTC)
TTC_final<-rbind(TTC_final,TTC)
TGC_final<-rbind(TGC_final,TGC)
GCC_final<-rbind(GCC_final,GCC)
GGC_final<-rbind(GGC_final,GGC)
CCC_final<-rbind(CCC_final,CCC)
ACC_final<-rbind(ACC_final,ACC)
TCC_final<-rbind(TCC_final,TCC)
TAC_final<-rbind(TAC_final,TAC)
CAG_final<-rbind(CAG_final,CAG)
AAC_final<-rbind(AAC_final,AAC)
CAC_final<-rbind(CAC_final,CAC)
GAG_final<-rbind(GAG_final,GAG)
GAC_final<-rbind(GAC_final,GAC)
AAG_final<-rbind(AAG_final,AAG)
CGC_final<-rbind(CGC_final,CGC)
TAG_final<-rbind(TAG_final,TAG)
ATA_final<-rbind(ATA_final,ATA)
CTA_final<-rbind(CTA_final,CTA)
GTA_final<-rbind(GTA_final,GTA)
GCA_final<-rbind(GCA_final,GCA)
GGA_final<-rbind(GGA_final,GGA)
CCA_final<-rbind(CCA_final,CCA)
ACA_final<-rbind(ACA_final,ACA)
TCA_final<-rbind(TCA_final,TCA)
CGA_final<-rbind(CGA_final,CGA)
TGA_final<-rbind(TGA_final,TGA)
CTG_final<-rbind(CTG_final,CTG)
GTG_final<-rbind(GTG_final,GTG)
GCG_final<-rbind(GCG_final,GCG)
GGG_final<-rbind(GGG_final,GGG)
CCG_final<-rbind(CCG_final,CCG)
ACG_final<-rbind(ACG_final,ACG)
TCG_final<-rbind(TCG_final,TCG)
CGG_final<-rbind(CGG_final,CGG)
TTA_final<-rbind(TTA_final,TTA)
AGT_final<-rbind(AGT_final,AGT)
AGA_final<-rbind(AGA_final,AGA)
TTG_final<-rbind(TTG_final,TTG)
AGC_final<-rbind(AGC_final,AGC)
AGG_final<-rbind(AGG_final,AGG)
}
print("frequecy calculation done : preparing output")
single_nucleoyide_count<-data.frame(cbind(A_final,T_final,G_final,C_final))
di_nucleotide_count<- data.frame(cbind(AA_final,TT_final,GG_final,CC_final,AT_final,AG_final,AC_final,TG_final,TC_final,TA_final,GA_final,GT_final,GC_final,CA_final,CT_final,CG_final))
tri_nucleotide_count<- data.frame(cbind(ATT_final,CTT_final,GTT_final,TTT_final,ATG_final,TGT_final,GCT_final,GGT_final,CCT_final,ACT_final,TCT_final,TAT_final,TGG_final,CAA_final,AAT_final,CAT_final,GAA_final,GAT_final,AAA_final,CGT_final,TAA_final,ATC_final,CTC_final,GTC_final,TTC_final,TGC_final,GCC_final,GGC_final,CCC_final,ACC_final,TCC_final,TAC_final,CAG_final,AAC_final,CAC_final,GAG_final,GAC_final,AAG_final,CGC_final,TAG_final,ATA_final,CTA_final,GTA_final,GCA_final,GGA_final,CCA_final,ACA_final,TCA_final,CGA_final,TGA_final,CTG_final,GTG_final,GCG_final,GGG_final,CCG_final,ACG_final,TCG_final,CGG_final,TTA_final,AGT_final,AGA_final,TTG_final,AGC_final,AGG_final))
names(single_nucleoyide_count)[1]<-paste("A")
names(single_nucleoyide_count)[2]<-paste("T")
names(single_nucleoyide_count)[3]<-paste("G")
names(single_nucleoyide_count)[4]<-paste("C")
names(di_nucleotide_count)[1]<-paste("AA")
names(di_nucleotide_count)[2]<-paste("TT")
names(di_nucleotide_count)[3]<-paste("GG")
names(di_nucleotide_count)[4]<-paste("CC")
names(di_nucleotide_count)[5]<-paste("AT")
names(di_nucleotide_count)[6]<-paste("AG")
names(di_nucleotide_count)[7]<-paste("AC")
names(di_nucleotide_count)[8]<-paste("TG")
names(di_nucleotide_count)[9]<-paste("TC")
names(di_nucleotide_count)[10]<-paste("TA")
names(di_nucleotide_count)[11]<-paste("GA")
names(di_nucleotide_count)[12]<-paste("GT")
names(di_nucleotide_count)[13]<-paste("GC")
names(di_nucleotide_count)[14]<-paste("CA")
names(di_nucleotide_count)[15]<-paste("CT")
names(di_nucleotide_count)[16]<-paste("CG")
names(tri_nucleotide_count)[1]<-paste("ATT")
names(tri_nucleotide_count)[2]<-paste("CTT")
names(tri_nucleotide_count)[3]<-paste("GTT")
names(tri_nucleotide_count)[4]<-paste("TTT")
names(tri_nucleotide_count)[5]<-paste("ATG")
names(tri_nucleotide_count)[6]<-paste("TGT")
names(tri_nucleotide_count)[7]<-paste("GCT")
names(tri_nucleotide_count)[8]<-paste("GGT")
names(tri_nucleotide_count)[9]<-paste("CCT")
names(tri_nucleotide_count)[10]<-paste("ACT")
names(tri_nucleotide_count)[11]<-paste("TCT")
names(tri_nucleotide_count)[12]<-paste("TAT")
names(tri_nucleotide_count)[13]<-paste("TGG")
names(tri_nucleotide_count)[14]<-paste("CAA")
names(tri_nucleotide_count)[15]<-paste("AAT")
names(tri_nucleotide_count)[16]<-paste("CAT")
names(tri_nucleotide_count)[17]<-paste("GAA")
names(tri_nucleotide_count)[18]<-paste("GAT")
names(tri_nucleotide_count)[19]<-paste("AAA")
names(tri_nucleotide_count)[20]<-paste("CGT")
names(tri_nucleotide_count)[21]<-paste("TAA")
names(tri_nucleotide_count)[22]<-paste("ATC")
names(tri_nucleotide_count)[23]<-paste("CTC")
names(tri_nucleotide_count)[24]<-paste("GTC")
names(tri_nucleotide_count)[25]<-paste("TTC")
names(tri_nucleotide_count)[26]<-paste("TGC")
names(tri_nucleotide_count)[27]<-paste("CC")
names(tri_nucleotide_count)[28]<-paste("GGC")
names(tri_nucleotide_count)[29]<-paste("CCC")
names(tri_nucleotide_count)[30]<-paste("ACC")
names(tri_nucleotide_count)[31]<-paste("TCC")
names(tri_nucleotide_count)[32]<-paste("TAC")
names(tri_nucleotide_count)[33]<-paste("CAG")
names(tri_nucleotide_count)[34]<-paste("AAC")
names(tri_nucleotide_count)[35]<-paste("CAC")
names(tri_nucleotide_count)[36]<-paste("GAG")
names(tri_nucleotide_count)[37]<-paste("GAC")
names(tri_nucleotide_count)[38]<-paste("AAG")
names(tri_nucleotide_count)[39]<-paste("CGC")
names(tri_nucleotide_count)[40]<-paste("TAG")
names(tri_nucleotide_count)[41]<-paste("ATA")
names(tri_nucleotide_count)[42]<-paste("CTA")
names(tri_nucleotide_count)[43]<-paste("GTA")
names(tri_nucleotide_count)[44]<-paste("GCA")
names(tri_nucleotide_count)[45]<-paste("GGA")
names(tri_nucleotide_count)[46]<-paste("CCA")
names(tri_nucleotide_count)[47]<-paste("ACA")
names(tri_nucleotide_count)[48]<-paste("TCA")
names(tri_nucleotide_count)[49]<-paste("GA")
names(tri_nucleotide_count)[50]<-paste("TGA")
names(tri_nucleotide_count)[51]<-paste("CTG")
names(tri_nucleotide_count)[52]<-paste("GTG")
names(tri_nucleotide_count)[53]<-paste("GCG")
names(tri_nucleotide_count)[54]<-paste("GGG")
names(tri_nucleotide_count)[55]<-paste("CCG")
names(tri_nucleotide_count)[56]<-paste("ACG")
names(tri_nucleotide_count)[57]<-paste("TCG")
names(tri_nucleotide_count)[58]<-paste("CGG")
names(tri_nucleotide_count)[59]<-paste("TTA")
names(tri_nucleotide_count)[60]<-paste("AGT")
names(tri_nucleotide_count)[61]<-paste("AGA")
names(tri_nucleotide_count)[62]<-paste("TTG")
names(tri_nucleotide_count)[63]<-paste("AGC")
names(tri_nucleotide_count)[64]<-paste("AGG")
print("writing matrix to output file")
single_nucleoyide_prob<- round(single_nucleoyide_count/10,2)
di_nucleotide_prob<- round(di_nucleotide_count/5,2)
tri_nucleotide_prob<- round(tri_nucleotide_count/3,2)
frequency_mat <- cbind(single_nucleoyide_prob,di_nucleotide_prob,tri_nucleotide_prob)
out_matrix_final<- cbind(matrix,frequency_mat)
print("writing feature_matrix to output file")
write.table(out_matrix_final,file="feature_matrix.txt",row.names=F,quote=F,sep="\t")
print("matrix created")
print("Error prediction process started")
print("loading train dataset")
train<- read.delim(args[3]) 
x<- data.matrix(train[,-11])
x<- data.matrix(train[,8:94])
y<- as.factor(train[,7])
print("Modelling started")
cvfit <-  cv.glmnet(x, y, family = "binomial", type.measure = "class", alpha =0.5)
result<- coef( cvfit, s =  cvfit$lambda.min, exact = T)
result<- capture.output(print(result))
print("writing attribute importance to output file")
write.table(result, file= "attribute_importance_matrix.txt", quote = F, row.names = F)
print("loading input_dataset")
x<- data.matrix(out_matrix_final[,6:92])
print("Predicting the errors in the Input_dataset")
prediction <- predict(cvfit, newx = x, s =  cvfit$lambda.min, type = "class")
prediction<- gsub("0","error",prediction)
prediction<- gsub("1","noerror",prediction)
final_matrix<- cbind(out_matrix_final[,4:6],prediction) 
print("writing predictions to the output file")
names(final_matrix)[1]<-paste("read_name")
names(final_matrix)[2]<-paste("base_position")
names(final_matrix)[3]<-paste("read_length")
names(final_matrix)[4]<-paste("prediction")
write.table(final_matrix, file= args[2], quote = F, row.names = F,sep = "\t")
print(" Error Prediction Sucessful !!! ")