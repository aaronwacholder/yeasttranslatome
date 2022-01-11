####Code to analyze data and make figures for manuscript

##Load premade figures into R
##Figure 1A: Conceptual Schema
require("ggpubr")
require("grImport") 
PostScriptTrace("input_files/Schema_RiboInt6.eps",charpath=FALSE)
f1a_concept_arial<-readPicture("Schema_RiboInt6.eps.xml")

pdf("f1a_concept_arial.pdf")
grid.picture(f1a_concept_arial, gp=gpar(fontfamily="ArialMT"))
dev.off()
f1a<-pictureGrob(f1a_concept_arial) 

##Figure 1B: Methods Schema

require("ggpubr")
require("grImport") 
PostScriptTrace("input_files/Figure1_schema4_nolabel.eps",charpath=FALSE)
f1b_arial<-readPicture("Figure1_schema3_nolabel.eps.xml")

pdf("f1b_schema_arial.pdf")
grid.picture(f1b_arial, gp=gpar(fontfamily="ArialMT"))
dev.off()
f1b<-pictureGrob(f1b_arial) 

##Constructing ORF candidate set, data for Figure 1B
##To identify all ORFs in the yeast genome:
##Run command: YeastTranslatome -GetAllORFs
##output files: orfs_comp orf_overlap
##all coordinates are 0-indexed, so S. cer Chromosome I is contig 0

orfs<-read.csv("orfs_comp",sep=" ") ##all ATG-stop ORFs in yeast
orfs$length<-orfs$end-orfs$start+1

#identify set of ORFs "all_index" that are potential candidates for translation. 
all_index<-which(orfs$splice_gene=="X" & orfs$contig<16 & orfs$is_gene!="YDL185W") #ORFs overlapping spliced genes are excluded. "YDL185W" contains entein and is excluded. We also exclude ORFs on contig 16, which is mitochondria
overlaps<-read.csv("orf_overlap",sep=" ") #a file giving info on what annotation each ORF overlaps

#ORFs are candidates if they do not overlap an annotated gene on the same strand
candidate_index<-all_index[which((overlaps$sense_ver+overlaps$sense_unchar+overlaps$sense_te+overlaps$sense_blocked)[all_index]==0)]

#distinguishing antisense from intergene ORFs
antisense_index<-candidate_index[which((overlaps$anti_ver+overlaps$anti_unchar+overlaps$anti_te+overlaps$anti_blocked)[candidate_index]>0 )]
intergene_index<-candidate_index[which((overlaps$anti_ver+overlaps$anti_unchar+overlaps$anti_te+overlaps$anti_blocked)[candidate_index]==0 )]

#list of ORFs of various annotation categories
verified_all<-all_index[which(orfs$orf_class[all_index]=="Verified")]
unchar_all<-all_index[which(orfs$orf_class[all_index]=="Uncharacterized")]
dubious_all<-all_index[which(orfs$orf_class[all_index]=="Dubious")]
pseudogene_all<-all_index[which(orfs$orf_class[all_index]=="pseudogene")]
te_all<-all_index[which(orfs$orf_class[all_index]=="transposable_element_gene")]

#list of ORF candidates of various annotation categories
verified_index<-candidate_index[which(orfs$orf_class[candidate_index]=="Verified")]
unchar_index<-candidate_index[which(orfs$orf_class[candidate_index]=="Uncharacterized")]
dubious_index<-candidate_index[which(orfs$orf_class[candidate_index]=="Dubious")]
pseudogene_index<-candidate_index[which(orfs$orf_class[candidate_index]=="pseudogene")]
te_index<-candidate_index[which(orfs$orf_class[candidate_index]=="transposable_element_gene")]

unannotated_index<-candidate_index[which(orfs$orf_class[candidate_index]=="None"|orfs$orf_class[candidate_index]=="ARS")]

noncanonical_index<-c(dubious_index,unannotated_index,pseudogene_index)
noncanonical_intergene_index<-intersect(noncanonical_index,intergene_index) 
noncanonical_antisense_index<-intersect(noncanonical_index,antisense_index) 

canonical_index<-c(verified_index,unchar_index,te_index)
canonical_intergene_index<-intersect(canonical_index,intergene_index) 
canonical_antisense_index<-intersect(canonical_index,antisense_index) 

length(candidate_index) #number "candidates" for translation
length(noncanonical_index) #number noncanonical "candidates" for translation
length(canonical_index) #number canonical "candidates" for translation
length(c(canonical_index,noncanonical_index))

##Figure 1C: Power of accumulated studies: figures showing real ribosome profiling reads across ORFs 
#table of processed ribo-seq reads on + strand. 
mrf<-read.csv("input_files/mapped_ribseqs_combo_all_f",sep=" ")

experiments<-read.csv("input_files/riboseq_experiments_dataset.txt",sep="\t") #table with data describing experiments from which ribo-seq reads were retrieved
length(table(experiments$PMID)) #count studies
length(experiments[,1])  #count experiments
all_srr<-unique(gsub("_.*","",colnames(mrf)))
all_srr<-all_srr[-(1:3)]
length(unique(gsub("_.*","",colnames(mrf))))-3 #count experiments passing qc

#reorganizing table to facilitate plotting using ggplot
get_orf_reads_by_experiment<-function(orf_range)
{
	exp_names<-unique(gsub("_.*","",colnames(mrf)))
	reads_by_exp<-array(dim=c(length(orf_range),length(exp_names)))
	colnames(reads_by_exp)<-exp_names
	for(i in 1:length(colnames(reads_by_exp)))
	{
		colsel<-which(gsub("_.*","",colnames(mrf))==colnames(reads_by_exp)[i])
		if(length(colsel)>1)
		{
			reads_by_exp[,i]<-rowSums(mrf[orf_range,colsel])
		}
		else
		{
			reads_by_exp[,i]<-mrf[orf_range,colsel]
		}
	}
	return(reads_by_exp)
}

#assembling read data for an individual ORF to plot
library("tidyr")
orf_id<-9222

orf_start<-which(mrf$contig==orfs[orf_id,]$contig & mrf$pos==orfs[orf_id,]$start)
orf_end<-which(mrf$contig==orfs[orf_id,]$contig & mrf$pos==orfs[orf_id,]$end)
orf_mrf<-get_orf_reads_by_experiment(orf_start:orf_end)#mrf_by_exp[orf_start:orf_end,]
orf_mrf_gathered<-gather(as.data.frame(orf_mrf),"experiment","reads",4:length(colnames(orf_mrf)))

orf_length<-length(unique(orf_mrf_gathered$pos))
in_frame<-min(orf_mrf_gathered$pos)+(1:(orf_length/3))*3-3

orf_mrf_gathered<-subset(orf_mrf_gathered,reads>0)
orf_mrf_gathered$frame<-as.character((orf_mrf_gathered$pos-min(orf_mrf_gathered$pos))%%3+1)
orf_mrf_gathered$pos<-orf_mrf_gathered$pos+1


all_exp_names<-names(table(orf_mrf_gathered$experiment))
reads_counts<-array()
reads_max<-array()
for(i in 1:length(all_exp_names))
{
	reads_counts[i]<-sum(orf_mrf_gathered[orf_mrf_gathered$experiment==all_exp_names[i],]$reads)
	reads_max[i]<-max(orf_mrf_gathered[orf_mrf_gathered$experiment==all_exp_names[i],]$reads)
}

FIGNUM<-5 #Number of figures to show in plot
reads_max<-reads_max[rev(order(reads_counts))[1:FIGNUM]]
all_exp_names<-all_exp_names[rev(order(reads_counts))[1:FIGNUM]]
unique_srr<-length(unique(experiments$PMID[which(experiments$SRR %in% all_exp_names)]))

require("grid")
require("gridExtra")
require("ggplot2")
fig_layers<-list()

frameColors <-
  setNames( c('black', '#66FF33', 'orange')
            , c('1','2','3')  )

frameColors <-
  setNames( c('black', 'chocolate2', 'chartreuse4')
            , c('1','2','3')  )



ylabs<-c("","","","","")

top_margin<-c(25,5.5,5.5,5.5,5.5)
for(i in 1:FIGNUM)
{
	gdata<-orf_mrf_gathered[orf_mrf_gathered$experiment==all_exp_names[i],]
	fig_layers[[i]]<-
	ggplot() +
		theme_classic(base_size=16) +
		ggtitle(all_exp_names[i])+
		theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin = margin(top_margin[i], 5.5, 5.5, 5.5,unit="pt"),plot.title = element_text(size=10, hjust = 0.5)) +
		geom_bar(data=gdata,aes(x=pos,y=reads,fill=frame),position="stack", stat ="identity",width=1)+
		scale_x_continuous(limits=c(orfs[orf_id,]$start,orfs[orf_id,]$end))+
		scale_y_continuous(breaks=c(0,max(gdata$reads)))+
		labs(y=ylabs[i])+
		scale_fill_manual(values=frameColors) 
}

fig_arrow<-
ggplot()+
theme_classic()+
geom_segment(aes(x=604710,y=50, xend=604710,yend=10), arrow = arrow(length = unit(0.5, "cm")))+
scale_x_continuous(limits=c(orfs[orf_id,]$start,orfs[orf_id,]$end))+
annotate(geom = "text", x = 604725+6, y = 30, label = "+95 additional experiments", size = 5)+
theme(axis.title = element_blank(),axis.text = element_blank(), axis.ticks=element_blank(),panel.background = element_blank(),panel.grid = element_blank(), axis.line=element_blank())

png("arrow.png")
fig_arrow
dev.off()

fig_bottom<-
ggplot() +
	theme_classic(base_size=16) +
	ggtitle("Combined reads")+
	theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5,unit="pt"),plot.title = element_text(size=16, hjust = 0.5)) +
	geom_bar(data=orf_mrf_gathered,aes(x=pos,y=reads,fill=frame),position="stack", stat ="identity")+
	scale_x_continuous(limits=c(orfs[orf_id,]$start,orfs[orf_id,]$end),breaks=orfs[orf_id,]$start+c(1,15,30,45,60,75),labels=c("1\nStart codon",15,30,45,60,"75\nStop codon"))+
	labs(y="",x="ORF Position",fill="Position\nin codon")+
	scale_fill_manual(values=frameColors) #c("black","red","blue"))

require("cowplot")

legend <- get_legend(
  fig_bottom + theme(legend.box.margin = margin(0, 15, 0, 12), legend.title=element_text(size=12))
)

f1b_nolegend<-plot_grid(fig_layers[[1]],fig_layers[[2]],fig_layers[[3]],fig_layers[[4]],fig_layers[[5]], #fig_layers[[6]],fig_layers[[7]],fig_layers[[8]],fig_layers[[9]],fig_layers[[10]], 
	fig_arrow,fig_bottom+theme(legend.position = "none"), ncol = 1, align = "v",rel_heights=c((4+reads_max[1:FIGNUM]),10,40 ) )#=c(rep(1,length(fig_layers)),3))

y.grob <- textGrob("Read count", 
                   gp=gpar(fontsize=16), rot=90)

png("Figure1C_orf_reads_.png",height=1200,width=800)
f1c<-plot_grid(f1b_nolegend,legend, rel_widths = c(3, .4))
f1c<-grid.arrange(arrangeGrob(f1c,left=y.grob))
dev.off()

####Figure 1D-E: FDR calculations for translated ORF database

ribo<-read.csv("riboseq_orfs",sep=" ") #data on ribosome profiling reads on each ORF. File is generated using command: YeastTranslatome -IdentifyTranslatedORFs

num_hits_noncanonical<-array()
scrambled_hits_noncanonical<-array()
num_hits_canonical<-array()
scrambled_hits_canonical<-array()
#number of true hits and scrambled (negative control) hits at range of pval thresholds, used to calculate FDR
for(i in 1:200)
{
	num_hits_noncanonical[i]<-length(which(ribo$pval[noncanonical_index]<i/1000))
	scrambled_hits_noncanonical[i]<-length(which(ribo$scram_pval[noncanonical_index]<i/1000))
	num_hits_canonical[i]<-length(which(ribo$pval[canonical_index]<i/1000))
	scrambled_hits_canonical[i]<-length(which(ribo$scram_pval[canonical_index]<i/1000))
}
df_hits<-data.frame(pvals=(1:200)/1000,hits=c(num_hits_noncanonical,scrambled_hits_noncanonical),hit_type=c(rep("Actual",length(num_hits_noncanonical)),rep("Scrambled",length(scrambled_hits_noncanonical))))
require("scales")
png("fdr_tpm.png",height=480,width=480)
f1d<-ggplot(df_hits,aes(x=pvals,y=hits,col=hit_type))+
	geom_line(size=3)+
	scale_color_manual(values=c("blue","black"))+
	theme_classic(base_size=16)+
	theme(legend.title=element_blank(),legend.position=c(.8,.5),plot.margin = unit(c(5.5,5.5,5.5,20), "pt"))+
	labs(x="P-value threshold",y="Translated nORFs found")+
	geom_vline(xintercept=which.min(abs(scrambled_hits_noncanonical/num_hits_noncanonical-.05))/1000, linetype="dashed", color="black", size=1.0)+
	scale_y_continuous(labels=comma)
f1d
dev.off()

df_hits_canonical<-data.frame(pvals=(1:200)/1000,hits=c(num_hits_canonical,scrambled_hits_canonical),hit_type=c(rep("Actual",length(num_hits_canonical)),rep("Scrambled",length(scrambled_hits_canonical))))
require("scales")
png("fdr_tpm_canon.png",height=480,width=480)
f1e<-ggplot(df_hits_canonical,aes(x=pvals,y=hits,col=hit_type))+
	geom_line(size=3)+
	scale_color_manual(values=c("blue","black"))+
	theme_classic(base_size=16)+
	theme(legend.title=element_blank(),legend.position=c(.8,.5),plot.margin = unit(c(5.5,5.5,20.5,20), "pt"))+
	labs(x="P-value threshold",y="Translated cORFs found")+
	scale_y_continuous(labels=comma)
f1e
dev.off()

# accum<-list()
# for(i in 1:10)
# {
	# accum[[i]]<-read.csv(paste("/home/acwach/YeastTest/riboseq_orfs_accumulated_studies_",i-1,sep=""),sep=" ")

# }
# orfs_accum<-array()
# for(i in 1:length(accum))
# {
	# orfs_accum[i]<-length(which(accum[[i]]$pval[noncanonical_index]<.037))
# }


f1de<-plot_grid(f1d,f1e,nrow=2,labels=c('D','E'),label_size=30)

require("cowplot")
pdf("Figure1__top.pdf",width=1000*3/200,height=1300*3/200)
plot_grid(f1a,f1b, f1c, f1de, ncol = 2, nrow=2, rel_widths=c(1,1),rel_heights=c(1,1.2),labels = c('A', 'B', 'C', ''), label_size=30)
dev.off()

###Figure 2

#ORFs considered translated if they have p-value < .037, corresponding to 5% FDR as found in graph generated above
trans_index<-which(ribo$pval<.037)
nontrans_index<-which(ribo$pval>=.037)

#lists of translated ORFs
candidate_trans<-intersect(candidate_index,trans_index)
noncanonical_trans<-intersect(noncanonical_index,trans_index)
canonical_trans<-intersect(canonical_index,trans_index)
noncanonical_intergene_trans<-intersect(noncanonical_intergene_index,trans_index)
noncanonical_antisense_trans<-intersect(noncanonical_antisense_index,trans_index)

#rates at which canonical and noncanonical ORFs are identified as translated
length(canonical_trans)/length(canonical_index) #freq canonical translation
length(noncanonical_trans)/length(noncanonical_index) #freq noncanonical translation

##A.	Number ORFs identified as translated by annotation class, distinguishing intergenic vs. overlapping
verified_trans<-intersect(verified_index,candidate_trans) 
unchar_trans<-intersect(unchar_index,candidate_trans) 
dubious_trans<-intersect(dubious_index,candidate_trans) 
unannotated_trans<-intersect(unannotated_index,candidate_trans)
te_trans<-intersect(te_index,candidate_trans)
pseudogene_trans<-intersect(pseudogene_index,candidate_trans)

annotation_class_data<-data.frame(category=c("Verified","Uncharacterized","Transposable element","Dubious","Pseudogene","Unannotated"),
               proportion=c(
			   length(verified_trans)/length(verified_index),
			   length(unchar_trans)/length(unchar_index),
			   length(te_trans)/length(te_index),
			   length(dubious_trans)/length(dubious_index),
			   length(pseudogene_trans)/length(pseudogene_index),
			   length(unannotated_trans)/length(unannotated_index)),
               counts=c(length(verified_trans),length(unchar_trans),length(te_trans),length(dubious_trans),length(pseudogene_trans),length(unannotated_trans)),
			   coverage=c(sum(orfs$length[verified_trans]),sum(orfs$length[unchar_trans]),sum(orfs$length[te_trans]),sum(orfs$length[dubious_trans]),sum(orfs$length[pseudogene_trans]),
			   sum(orfs$length[unannotated_trans])),
			   canonical=c("Canonical","Canonical","Canonical","Noncanonical","Noncanonical","Noncanonical")
				)
annotation_class_data$category<-factor(annotation_class_data$category,levels=(c("Verified","Uncharacterized","Transposable element","Dubious","Pseudogene","Unannotated")))
annotation_class_data$percent<-paste(round(annotation_class_data$proportion*100),"%",sep="")
annotation_class_data$comma_count<-format(annotation_class_data$counts,big.mark=",",trim=T)

##Figures describing genomic context of translated ORFs 

font_size<-16
library("ggplot2")
require("grid")
require("gridExtra")
png("Figure2a_percent_calls.png",height=250,width = 400)
f2a<-ggplot(annotation_class_data)+
  theme_bw(base_size=font_size) +
  theme(legend.position="None", legend.title = element_blank(), legend.background = element_rect(fill="transparent"),plot.margin = unit(c(50,30,-15,15), "pt"),axis.text.x=element_text(angle = 45, hjust=1),axis.title=element_text(size=16))+
  geom_bar(aes(x=category,y=proportion*100,fill=canonical),stat="identity",width=.5)+
  scale_fill_manual(values=c("red","blue")) +
  labs(x="",y="Percent called translated")+
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))
f2a
dev.off()

png("Figure2b_counts_calls.png",height=250,width = 400)
f2b<-ggplot(annotation_class_data)+
  theme_bw(base_size=font_size) +
  theme(legend.position="None", legend.title = element_blank(),plot.margin = unit(c(50,30,-15,15), "pt"),axis.text.x=element_text(angle = 45, hjust=1))+
  geom_bar(aes(x=category,y=counts,fill=canonical),stat="identity",width=.5)+
    scale_fill_manual(values=c("red","blue")) +
  geom_text(stat='identity', aes(label=comma_count, x=category, y=counts), vjust=-.2) +
  labs(x="",y="Counts of translated ORFs")+
  scale_y_continuous(limits=c(0,22000),breaks=c(0,10000,20000),labels=c("0","10,000","20,000"))#+
f2b
dev.off()

#nORF vs cORF on length and translation rates

noncan_can<-c(noncanonical_trans,canonical_trans)

df_orflength<-data.frame(orflength=orfs$length[noncan_can],trans_rate=(ribo$frame0/orfs$length)[noncan_can],categ="Canonical",stringsAsFactors=F)
df_orflength$categ[1:length(noncanonical_trans)]<-"Noncanonical"
require("ggplot2")
png("Figure2c_orflengths.png")
f2c<-ggplot(df_orflength)+
  theme_classic(base_size=font_size) +
  theme(legend.title = element_blank(),legend.position=c(.25,.95))+
  geom_density(aes(orflength, fill=categ),alpha=.75)+
  labs(x="ORF length (bp)",y="Density")+
  scale_fill_manual(values=c("red","blue"))+
  scale_x_log10(breaks=c(10,100,1000,10000),labels=c("10","100","1,000","10,000"))+
  scale_y_continuous(limits=c(0,1.5))
f2c
dev.off()

require("ggplot2")
require("scales")
png("Figure2d_orf_rates.png")
f2d<-ggplot(df_orflength)+
  theme_classic(base_size=font_size) +
  theme(legend.title = element_blank(),legend.position=c(.25,.95))+
  geom_density(aes(trans_rate, fill=categ),alpha=.75)+
  labs(x="ORF translation rate (reads/length)",y="Density")+
  scale_fill_manual(values=c("red","blue"))+
  scale_x_log10(breaks=c(.01,1,100,10000),labels=c(".01","1","100","10,000"))
f2d
dev.off()

###Assembled Figure2 
require("cowplot")
png("Figure2_071221.png",width=2400/.75,height=2400/.75, res=300)
plot_grid(f2a,f2b,f2c,f2d,nrow=2,ncol=2,labels=c("A","B","C","D"),label_size=30, align="v")
dev.off()

###Supplementary Figure 1: translation patterns among ORFs highly correlated between studies

#read files giving ribo-seq information for each individual study
ribofiles<-list.files()
ribofiles<-ribofiles[grep("riboseq_orfs_",ribofiles)]#[1:42]
ribostudies<-list()
for(i in 1:length(ribofiles))
{
	ribostudies[[i]]<-read.csv(ribofiles[i],sep=" ")
}
 
filtered_ribostudies<-list()
filtered_ribofiles<-array()
ribothres<-0
n<-1
for(i in 1:length(ribostudies))
{
	if(sum(ribostudies[[i]]$frame0)>ribothres)
	{
		filtered_ribostudies[[n]]<-ribostudies[[i]]
		filtered_ribofiles[n]<-ribofiles[i]
		n<-n+1
	}
}

#correlation between riboseq reads per ORF betwen different studdies
ribocor<-array(dim=c(length(filtered_ribostudies),length(filtered_ribostudies)))
sumreads<-array()
avgcor<-array()
agree_count<-array(dim=c(length(filtered_ribostudies),length(filtered_ribostudies)))
for(i in 1:length(filtered_ribostudies))
{
	for(j in 1:length(filtered_ribostudies))
	{
		ribocor[i,j]<-cor(filtered_ribostudies[[i]]$frame0[candidate_index]/orfs$length[candidate_index],filtered_ribostudies[[j]]$frame0[candidate_index]/orfs$length[candidate_index])
		agree_count[i,j]<-length(intersect(which(filtered_ribostudies[[i]]$pval[candidate_index]<.02),which(filtered_ribostudies[[j]]$pval[candidate_index]<.02)))
	}
	sumreads[i]<-sum(filtered_ribostudies[[i]]$frame0)
	avgcor[i]<-mean(ribocor[i,])
}

colnames(ribocor)<-substr(filtered_ribofiles,14,200)
rownames(ribocor)<-substr(filtered_ribofiles,14,200)

require("heatmap3")
#cluster studies based on correlation matrix
h3<-heatmap3(ribocor,symm=TRUE,returnDistMatrix=T)

dfcor<-data.frame(rcor=avgcor, reads=sumreads/1000000, name=colnames(ribocor),clust_class=cutree( hclust(h3$DistMatrixR),8),id=1:length(avgcor))
dfcor[order(dfcor$clust_class),]

for(i in 1:length(dfcor$rcor))
{
	sel<-which(experiments$SRP==dfcor$name[i])[1]
	dfcor$chx[i]<-as.character(experiments$CHX[sel])
	dfcor$ypd[i]<-as.character(experiments$YPD[sel])
	dfcor$media[i]<-as.character(experiments$Growth.media[sel])
	dfcor$pmid[i]<-as.character(experiments$PMID[sel])
}

mean(ribocor[dfcor$clust_class==1,dfcor$clust_class==1])
length(which(dfcor$clust_class==1))

#get author from pmid file
pmid_map<-read.csv("input_files/pmid_to_author.csv",header=F,stringsAsFactors=F)

colnames(ribocor)<-pmid_map[match(dfcor$pmid,pmid_map[,1]),2]
rownames(ribocor)<-pmid_map[match(dfcor$pmid,pmid_map[,1]),2]

all_pmids<-unique(experiments$PMID)

experiments_count<-as.numeric(table(experiments$PMID)[as.character(all_pmids)])
experiments_passed<-as.numeric(table(experiments$PMID[which(experiments$SRR %in% all_srr)])[as.character(all_pmids)])

studies_table<-data.frame(pmids=all_pmids,experiments_count=experiments_count,experiments_passed=experiments_passed,ref=pmid_map[match(all_pmids,pmid_map[,1]),2])
studies_table[which(is.na(experiments_passed)),]$experiments_passed<-0

#Basis for SupplementaryTable2
write.table(studies_table,"studies_table",quote=F)

#SupplementaryFigure 1A: correlation in translation patterns between studies
ribocor_reordered<-ribocor[order(dfcor$clust_class),order(dfcor$clust_class)]
require("ggcorrplot")
png("SuppFigure1A_riboheat2.png",height=600,width=600)
s1a<-ggcorrplot(ribocor_reordered,type="full", tl.srt = 90)+
     scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1),name="Correlation")
s1a
dev.off()

df_cor<-data.frame(r18=filtered_ribostudies[[18]]$frame0/orfs$length,r16=filtered_ribostudies[[16]]$frame0/orfs$length)

require(scales)

png("SuppFigure1B_GerVsAlb.png")
s1b<-ggplot(df_cor)+
	geom_point(aes(x=r18,y=r16))+
	theme_classic(base_size=20)+
	theme(plot.margin = unit(c(5.5,5.5,5.5,20), "pt"))+
	scale_x_log10(limits=c(.001,10000),breaks=c(.01,1,100,10000),labels = c(".01","1","100","10,000"))+
	scale_y_log10(limits=c(.001,10000),breaks=c(.01,1,100,10000),labels = c(".01","1","100","10,000"))+
	labs(x="Reads per base, Gerashchenko et al 2014",y="Reads per base, Albert et al 2014")
s1b
dev.off()

require("cowplot")
pdf("SuppFigure1__.pdf",width=2000*2/200,height=700*3/200)
plot_grid(s1a,s1b, ncol = 2, nrow=1, rel_widths=c(1,1),labels = c('A', 'B'), label_size=30)
dev.off()


####Figure 3: replicability between studies and CHX analysis

#Figure 3A: replicability
pval_thres<-.037
rib16<-(which(filtered_ribostudies[[18]]$pval[noncanonical_index]<pval_thres))

hit_count<-array()
match_percent<-array()
for(i in 1:length(filtered_ribostudies))
{
	ribi<-(which(filtered_ribostudies[[i]]$pval[noncanonical_index]<pval_thres))
	match_percent[i]<-length(intersect(rib16,ribi))/length(ribi)
	hit_count[i]<-length(ribi)
}
sel<-which(dfcor$clust_class==1 & hit_count>400 &match_percent<1)

df_replic<-data.frame(agree=c(match_percent[sel],length(rib16)/length(noncanonical_index)),paper=c(colnames(ribocor)[sel],"Random expectation"))
df_replic$paper<-factor(df_replic$paper,levels=as.character(df_replic$paper))
df_replic<-df_replic[!is.na(df_replic$paper),]
png("agree.png")
f3a<-ggplot(df_replic)+
	theme_classic(base_size=16)+
	geom_bar(aes(x=paper,y=agree*100),stat="identity",fill="blue")+
	scale_y_continuous(limits=c(0,100))+
	coord_flip()+
	labs(x="Study",y="Percent nORFs replicated in Gerashchenko et al. 2014 data")+
	theme(axis.title.y = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
f3a
dev.off()

#CHX analysis

chx0<-read.csv("riboseq_orfs_chx_sampled0",sep=" ") #ribo-seq data no CHX experiments
chx1<-read.csv("riboseq_orfs_chx_sampled1",sep=" ") #ribo-seq data for CHX experiments

rel_chx_reads<-array()
rel_chx_reads[1]<-(sum(chx1$frame0[canonical_index])/sum(chx1$frame0))/(sum(chx0$frame0[canonical_index])/sum(chx0$frame0))
rel_chx_reads[2]<-(sum(chx1$frame0[noncanonical_index])/sum(chx1$frame0))/(sum(chx0$frame0[noncanonical_index])/sum(chx0$frame0))

genomic_context_data<-data.frame(category=c("cORF","nORF"),rel_chx=rel_chx_reads,canonical=c("Canonical","Noncanonical"))
genomic_context_data$category<-factor(genomic_context_data$category,levels=c("cORF","nORF"))

png("rel_chx.png",height=480,width = 400)
f3b<-ggplot(genomic_context_data)+
  theme_classic(base_size=16) +
  theme(legend.position="none",  axis.text.x=element_text(angle = 45, hjust=1), plot.margin = unit(c(45,5.5,5.5,40), "pt"))+
  geom_bar(aes(x=category,y=rel_chx,fill=canonical),stat="identity")+
  labs(x="",y="Read count enrichment CHX vs. no CHX")+
  scale_fill_manual(values=c("red","blue"))
f3b
dev.off()

chx0_index<-which(chx0$pval[candidate_index]<.03)
chx1_index<-which(chx1$pval[candidate_index]<.03)
chx1_excl_index<-candidate_index[which(chx1$pval[candidate_index]<.03 & chx0$pval[candidate_index]>.03)]

df<-data.frame(favored_position=c("1","2","3"), hits=c(sum(chx0$hits0[chx1_excl_index]),sum(chx0$hits1[chx1_excl_index]),sum(chx0$hits2[chx1_excl_index])))

png("trip_periodic.png",height=480,width = 400)
f3c<-ggplot(df)+
  theme_classic(base_size=16) +
  theme(legend.position="none",plot.margin = unit(c(5.5,5.5,5.5,40), "pt"))+
  geom_bar(aes(x=favored_position,y=hits),stat="identity")+
  labs(x="Majority position within codon",y="Count")+
  scale_y_continuous(breaks=c(5000,10000,15000),labels=c("5,000","10,000","15,000"))
f3c
dev.off()

png("Figure3_071221.png",width=1500*3.7,height=750*2.9,res=300)
plot_grid(f3a,f3b,f3c,ncol=3,nrow=1,label_size=30,labels=c("A","B","C"),rel_widths=c(1,.45,.45),align="h")
dev.off()


#####Figure 4: environmental and genomic context 

#Figure 4A: schema for genomic contexts
require("grImport") 

PostScriptTrace("genomic_context12.eps",charpath=FALSE)
f4a_arial<-readPicture("genomic_context12.eps.xml")

pdf("f4a_schema_arial.pdf")
grid.picture(f4a_arial, gp=gpar(fontfamily="ArialMT"))
dev.off()
f4a<-pictureGrob(f4a_arial,exp=.05)

#figure 4B-C: counts and percentages of translated ORFs in different contexts
shared<-read.csv("shared_transcripts",sep=" ")#info on annotations on same transcript as each ORF

uORF_index<-intersect(noncanonical_intergene_index,which(orfs$orf_id%in%shared$orf_id[shared$relation=="upstream" & shared$sgd_class=="gene"])) 
dORF_index<-intersect(noncanonical_intergene_index,which(orfs$orf_id%in%shared$orf_id[shared$relation=="downstream" & shared$sgd_class=="gene"]))
nORF_index<-intersect(noncanonical_intergene_index,which(!(orfs$orf_id%in%shared$orf_id[shared$sgd_class=="gene"])))#noncanonical_intergene_index[which(orfs$transcript_context[noncanonical_intergene_index]=="")]
wORF_index<-noncanonical_antisense_index[which(orfs$gene_overlap_length[noncanonical_antisense_index]==orfs$length[noncanonical_antisense_index]-1)]
oORF_index<-noncanonical_antisense_index[which(orfs$gene_overlap_length[noncanonical_antisense_index]<orfs$length[noncanonical_antisense_index]-1)]
tRNA_index<-intersect(noncanonical_intergene_index,which(orfs$orf_id%in%shared$orf_id[shared$sgd_class=="tRNA_gene"]))
snoRNA_index<-intersect(noncanonical_intergene_index,which(orfs$orf_id%in%shared$orf_id[shared$sgd_class=="snoRNA_gene"]))
ncRNA_index<-intersect(noncanonical_intergene_index,which(orfs$orf_id%in%shared$orf_id[shared$sgd_class=="ncRNA_gene"]))
tORF_index<-which(orfs$tseqs0+orfs$tseqs1>0)

uORF_trans<-intersect(uORF_index,candidate_trans) 
dORF_trans<-intersect(dORF_index,candidate_trans) 
nORF_trans<-intersect(nORF_index,candidate_trans) 
oORF_trans<-intersect(oORF_index,candidate_trans) 
wORF_trans<-intersect(wORF_index,candidate_trans) 
tRNA_trans<-intersect(tRNA_index,candidate_trans)
snoRNA_trans<-intersect(snoRNA_index,candidate_trans)
ncRNA_trans<-intersect(ncRNA_index,candidate_trans)

length(uORF_trans)/length(uORF_index) #uORF translation rate
length(dORF_trans)/length(dORF_index) #dORF translation rate
length(nORF_trans)/length(nORF_index) #nORF translation rate
length(oORF_trans)/length(oORF_index) #oORF translation rate
length(wORF_trans)/length(wORF_index) #wORF translation rate

genomic_context_data<-data.frame(category=c("uORF","dORF","Independent","Full overlap","Partial overlap","tRNA","snoRNA","ncRNA"),
counts=c(length(uORF_trans),length(dORF_trans),length(nORF_trans),length(wORF_trans),length(oORF_trans),length(tRNA_trans),length(snoRNA_trans),length(ncRNA_trans)),
prop=c(length(uORF_trans)/length(uORF_index),length(dORF_trans)/length(dORF_index),length(nORF_trans)/length(nORF_index),length(wORF_trans)/length(wORF_index),length(oORF_trans)/length(oORF_index),length(tRNA_trans)/length(tRNA_index),length(snoRNA_trans)/length(snoRNA_index),length(ncRNA_trans)/length(ncRNA_index)))
genomic_context_data$category<-factor(genomic_context_data$category,levels=c("uORF","dORF","Independent","Full overlap","Partial overlap","tRNA","snoRNA","ncRNA"))

genomic_context_data$comma_count<-format(genomic_context_data$counts,big.mark=",",trim=T)

dorf_vs_uorf<-rbind(c(length(uORF_trans),length(uORF_index)-length(uORF_trans)),c(length(dORF_trans),length(dORF_index)-length(dORF_trans)))
fisher.test(dorf_vs_uorf)

png("context_calls.png",height=380,width = 400)
f4b<-ggplot(genomic_context_data)+
  theme_bw(base_size=20) +
  theme(legend.position="none",plot.margin = unit(c(15,30,-15,15), "pt"),axis.text.x=element_text(angle = 45, hjust=1))+
  geom_bar(aes(x=category,y=counts),stat="identity",fill="blue",width=.5)+
  geom_text(stat='identity', aes(label=comma_count, x=category, y=counts), vjust=-.2) +
  labs(x="",y="Counts of translated nORFs")+
  scale_x_discrete(limits = (levels(genomic_context_data$category)))+
  scale_y_continuous(limits=c(0,6650),breaks=c(0,2000,4000,6000),labels=c("0","2,000","4,000","6,000"))#+
f4b
dev.off()

png("context_calls_prop.png",height=380,width = 400)
f4c<-ggplot(genomic_context_data)+
  theme_bw(base_size=20) +
  theme(legend.position="none",plot.margin = unit(c(15,30,-15,15), "pt"),axis.text.x=element_text(angle = 45, hjust=1))+
  geom_bar(aes(x=category,y=prop*100,fill=category),stat="identity",fill="blue",width=.5)+
  labs(x="",y="Percent called translated")+
  scale_x_discrete(limits = (levels(genomic_context_data$category)))+
  scale_y_continuous(limits=c(0,30))
f4c
dev.off()

####sd vs ypd media

ypd0<-list()
ypd1<-list()
ypd2<-list()
ypd3<-list()
ypd4<-list()

#experiments in different media, all subsampled to have same total number reads
for(i in 1:20)
{
	ypd0[[i]]<-read.csv(paste("riboseq_orfs_ypd0_",i-1,sep=""),sep=" ") #experiments with media other than YPD 
	ypd1[[i]]<-read.csv(paste("riboseq_orfs_ypd1_",i-1,sep=""),sep=" ") #experiments with YPD media 
	ypd2[[i]]<-read.csv(paste("riboseq_orfs_ypd2_",i-1,sep=""),sep=" ") #all experiments with known media
	ypd3[[i]]<-read.csv(paste("riboseq_orfs_ypd3_",i-1,sep=""),sep=" ") #all experiments, regardless of whether media is known
	ypd4[[i]]<-read.csv(paste("riboseq_orfs_ypd4_",i-1,sep=""),sep=" ") #experiment with SD
}

ypd_counts0<-array()
ypd_counts1<-array()
ypd_counts2<-array()
ypd_counts3<-array()
ypd_counts4<-array()

for(i in 1:20)
{
	ypd_counts0[i]<-length(which(ypd0[[i]]$pval[noncanonical_index]<.037))
	ypd_counts1[i]<-length(which(ypd1[[i]]$pval[noncanonical_index]<.037))
	ypd_counts2[i]<-length(which(ypd2[[i]]$pval[noncanonical_index]<.037))	
	ypd_counts3[i]<-length(which(ypd3[[i]]$pval[noncanonical_index]<.037))	
	ypd_counts4[i]<-length(which(ypd4[[i]]$pval[noncanonical_index]<.037))	
}

#ypd vs sd media comparison
require("ggplot2")
df3<-data.frame(reads=(1:17)*10, counts=c(ypd_counts1[1:17],ypd_counts4[1:17]), category=c(rep("Rich media",17),rep("Minimal media",17)))
png("ypd_growth0_gg.png")
f4d<-ggplot(df3)+
  theme_classic(base_size=20) +
  theme(legend.position=c(.7,.3),legend.title = element_blank(), axis.title.x=element_text(vjust=18))+
  geom_line(aes(x=reads,y=counts,col=category),size=1.5)+
  labs(x="Reads (millions)",y="Translated nORFs identified")+
  scale_color_manual(values=c("black","chocolate1"))+
  scale_y_continuous(limits=c(0,13000),breaks=c(0,5000,10000),labels=c("0","5,000","10,000"))
f4d
dev.off()

ypd_trans<-noncanonical_index[which(ypd1[[20]]$pval[noncanonical_index]<.037)]
sd_trans<-noncanonical_index[which(ypd4[[20]]$pval[noncanonical_index]<.037)]


ypd4_specific<-which(ypd4[[17]]$pval[noncanonical_index]<.001 & ypd1[[17]]$pval[noncanonical_index]>.037)
ypd1_specific<-which(ypd1[[17]]$pval[noncanonical_index]<.001 & ypd4[[17]]$pval[noncanonical_index]>.037)
both_conditions<-which(ypd1[[17]]$pval[noncanonical_index]<.001 & ypd4[[17]]$pval[noncanonical_index]<.001)

ypdsd_categs<-c("Rich (q<.001)\nnot minimal (q>.05) ","Minimal (q<.001)\nnot rich (q>.05)","Rich (q<.001)\nand minimal (q<.001)")
df_ypdsd<-data.frame(count=c(length(ypd1_specific),length(ypd4_specific),length(both_conditions)),name=ypdsd_categs)
df_ypdsd$name<-factor(df_ypdsd$name, levels=ypdsd_categs)

png("ypdsd.png")
f4e<-ggplot(df_ypdsd)+
	theme_classic(base_size=20)+
	geom_bar(aes(x=name,y=count),stat="identity",fill="blue")+
	labs(x="",y="Translated nORFs identified")+
	theme(axis.text.x=element_text(size=10))+
	scale_y_continuous(limits=c(0,3000),breaks=c(0,500,1000,1500,2000,2500,3000),labels=c("0","500","1,000","1,500","2,000","2,500","3,000"))
f4e
dev.off()

####amino acid composition

untrans_controls<-intersect(intersect(nontrans_index,intergene_index),noncanonical_index)
untrans_controls_antisense<-intersect(intersect(nontrans_index,antisense_index),noncanonical_index)

aaseqs<-list()
for(i in 1:length(orfs[,1]))
{
	aaseqs[[i]]<-strsplit(as.character(borfs[i,]$focal_aaseq),split="")[[1]]
}
aaseqs_nostart<-list()
for(i in 1:length(orfs[,1]))
{
	aaseqs_nostart[[i]]<-aaseqs[[i]][-1]
}

canonical_comp<-table(unlist(aaseqs_nostart[canonical_intergene_index]))
control_comp<-table(unlist(aaseqs_nostart[intersect(untrans_controls,which(ribo$pval>.5))]))
transient_comp<-table(unlist(aaseqs_nostart[noncanonical_intergene_trans]))

canonical_comp<-canonical_comp[-20]/sum(canonical_comp[-20])
control_comp<-control_comp[-20]/sum(control_comp[-20])
transient_comp<-transient_comp[-20]/sum(transient_comp[-20])

transient_comp_uorf<-table(unlist(aaseqs_nostart[uORF_trans]))
transient_comp_dorf<-table(unlist(aaseqs_nostart[intersect(dORF_index,noncanonical_intergene_trans)]))
transient_comp_norf<-table(unlist(aaseqs_nostart[intersect(nORF_index,noncanonical_intergene_trans)]))
transient_comp_worf<-table(unlist(aaseqs_nostart[intersect(wORF_index,noncanonical_antisense_trans)]))
transient_comp_oorf<-table(unlist(aaseqs_nostart[intersect(oORF_index,noncanonical_antisense_trans)]))

transient_comp_uorf<-transient_comp_uorf[-20]/sum(transient_comp_uorf[-20])
transient_comp_dorf<-transient_comp_dorf[-20]/sum(transient_comp_dorf[-20])
transient_comp_norf<-transient_comp_norf[-20]/sum(transient_comp_norf[-20])
transient_comp_worf<-transient_comp_worf[-20]/sum(transient_comp_worf[-20])
transient_comp_oorf<-transient_comp_oorf[-20]/sum(transient_comp_oorf[-20])


require("tidyr")
df_comp<-data.frame(aa= names(transient_comp),trans=as.numeric(transient_comp), trans_uorf=as.numeric(transient_comp_uorf), trans_dorf=as.numeric(transient_comp_dorf), trans_norf=as.numeric(transient_comp_norf), trans_worf=as.numeric(transient_comp_worf), trans_oorf=as.numeric(transient_comp_oorf), canon=as.numeric(canonical_comp), cont=as.numeric(control_comp))
df_comp_gathered<-gather(df_comp, "orf_type","comp",2:8)
df_comp_gathered$aa<-factor(df_comp_gathered$aa,levels=names(canonical_comp)[order(canonical_comp/control_comp)])
df_comp_gathered$orf_type<-factor(df_comp_gathered$orf_type,levels=c("canon","trans","trans_uorf","trans_dorf","trans_norf","trans_oorf","trans_worf"))


library("ggplot2")
png("transient_composition_context.png")
f4f<-ggplot(df_comp_gathered,aes(x=aa, y=comp/cont, color=orf_type))+
	geom_point()+
	geom_line(aes(group = orf_type))+
	scale_color_manual(values=c("red","blue","green","brown","yellow","orange","purple"),labels = c("cORFs", "nORFs","uORFs","dORFs", "Independent nORFs", "Partial overlap antisense nORFs", "Full overlap antisense nORFs"  ))+
	theme_classic(base_size=20)+
	theme(legend.position=c(.4,.8),legend.title=element_blank(), axis.title.x=element_text(vjust=16))+
	labs(x="Encoded amino acid",y="Frequency relative to\nuntranslated controls")+
	scale_y_continuous(breaks=(1:40)/4)+
	geom_hline(yintercept=1.0, linetype="dashed", color="grey", size=1.0)
f4f
dev.off()

get_controls<-function(test_set,control_set)
{
	matched_test<-array()
	matched_control<-array()
	for(i in 1:length(test_set))
	{
		sel<-which(orfs$length[control_set]==orfs$length[test_set[i]])
		if(length(sel)==1)
		{
			matched_test[i]<-test_set[i]
			matched_control[i]<-control_set[sel]
		}
		if(length(sel)>1)
		{
			matched_test[i]<-test_set[i]
			matched_control[i]<-control_set[sample(sel,1)]
		}
	}
	matched_test<-na.omit(matched_test)
	matched_control<-na.omit(matched_control)
	return_set<-data.frame(test=matched_test,cont=matched_control)
	return(return_set)
}

uORF_matched<-get_controls(uORF_trans,intersect(uORF_index,untrans_controls))
transient_comp_uorf<-table(unlist(aaseqs_nostart[uORF_matched$test]))
transient_comp_uorf_control<-table(unlist(aaseqs_nostart[uORF_matched$cont]))

dORF_matched<-get_controls(dORF_trans,intersect(dORF_index,untrans_controls))
transient_comp_dorf<-table(unlist(aaseqs_nostart[dORF_matched$test]))
transient_comp_dorf_control<-table(unlist(aaseqs_nostart[dORF_matched$cont]))

nORF_matched<-get_controls(nORF_trans,intersect(nORF_index,untrans_controls))
transient_comp_norf<-table(unlist(aaseqs_nostart[nORF_matched$test]))
transient_comp_norf_control<-table(unlist(aaseqs_nostart[nORF_matched$cont]))

wORF_matched<-get_controls(wORF_trans,intersect(wORF_index,untrans_controls_antisense))
transient_comp_worf<-table(unlist(aaseqs_nostart[wORF_matched$test]))
transient_comp_worf_control<-table(unlist(aaseqs_nostart[wORF_matched$cont]))

oORF_matched<-get_controls(oORF_trans,intersect(oORF_index,untrans_controls_antisense))
transient_comp_oorf<-table(unlist(aaseqs_nostart[oORF_matched$test]))
transient_comp_oorf_control<-table(unlist(aaseqs_nostart[oORF_matched$cont]))

chisq.test(rbind(transient_comp_uorf,transient_comp_uorf_control))[3]
chisq.test(rbind(transient_comp_dorf,transient_comp_dorf_control))[3]
chisq.test(rbind(transient_comp_norf,transient_comp_norf_control))[3]
chisq.test(rbind(transient_comp_worf,transient_comp_worf_control))[3]
chisq.test(rbind(transient_comp_oorf,transient_comp_oorf_control))[3]


transient_comp_uorf<-transient_comp_uorf[-20]/sum(transient_comp_uorf[-20])
transient_comp_dorf<-transient_comp_dorf[-20]/sum(transient_comp_dorf[-20])
transient_comp_norf<-transient_comp_norf[-20]/sum(transient_comp_norf[-20])
transient_comp_worf<-transient_comp_worf[-20]/sum(transient_comp_worf[-20])
transient_comp_oorf<-transient_comp_oorf[-20]/sum(transient_comp_oorf[-20])

transient_comp_uorf_control<-transient_comp_uorf_control[-20]/sum(transient_comp_uorf_control[-20])
transient_comp_dorf_control<-transient_comp_dorf_control[-20]/sum(transient_comp_dorf_control[-20])
transient_comp_norf_control<-transient_comp_norf_control[-20]/sum(transient_comp_norf_control[-20])
transient_comp_worf_control<-transient_comp_worf_control[-20]/sum(transient_comp_worf_control[-20])
transient_comp_oorf_control<-transient_comp_oorf_control[-20]/sum(transient_comp_oorf_control[-20])

df_comp<-data.frame(aa= names(transient_comp),trans=as.numeric(transient_comp)/as.numeric(control_comp), trans_uorf=as.numeric(transient_comp_uorf)/as.numeric(transient_comp_uorf_control), trans_dorf=as.numeric(transient_comp_dorf)/as.numeric(transient_comp_dorf_control), trans_norf=as.numeric(transient_comp_norf)/as.numeric(transient_comp_norf_control), trans_worf=as.numeric(transient_comp_worf)/as.numeric(transient_comp_worf_control), trans_oorf=as.numeric(transient_comp_oorf)/as.numeric(transient_comp_oorf_control), canon=as.numeric(canonical_comp)/as.numeric(control_comp), cont=as.numeric(control_comp))
df_comp_gathered<-gather(df_comp, "orf_type","comp",2:8)
df_comp_gathered$aa<-factor(df_comp_gathered$aa,levels=names(canonical_comp)[order(canonical_comp/control_comp)])
df_comp_gathered$orf_type<-factor(df_comp_gathered$orf_type,levels=c("canon","trans","trans_uorf","trans_dorf","trans_norf","trans_oorf","trans_worf"))

library("ggplot2")
png("transient_composition_context_control_context.png")
f4g<-ggplot(df_comp_gathered,aes(x=aa, y=comp, color=orf_type))+
	geom_point()+
	geom_line(aes(group = orf_type))+
	scale_color_manual(values=c("red","blue","green","brown","yellow","orange","purple"),labels = c("cORFs", "nORFs","uORFs","dORFs", "Independent nORFs", "Partial overlap antisense nORFs", "Full overlap antisense nORFs"  ))+
	theme_classic(base_size=20)+
	theme(legend.position=c(.4,.8),legend.title=element_blank(), axis.title.x=element_text(vjust=16))+
	labs(x="Encoded amino acid",y="Frequency relative to\ncontext-matched controls")+
	scale_y_continuous(breaks=(1:40)/4)+
	geom_hline(yintercept=1.0, linetype="dashed", color="grey", size=1.0)
f4g
dev.off()

require("ggpubr")
require("grImport") 
require("cowplot")
f4bcde<-plot_grid(ncol=2,nrow=3,f4b,f4c,f4f,f4g,f4d,f4e,labels=c("B","C","D","E","F","G"),label_size=30,align="vh")

###Assembling 4
pdf("Figure4_101221.pdf",width=2900/200,height=4400/200)
plot_grid(f4a,f4bcde,ncol=1,nrow=2,labels=c("A",""),label_size=30, rel_heights=c(.4,2))
dev.off()

##############                                              #########################
######################### EVOLUTIONARY ANALYSIS BEGINS ##############################
##############                                              #########################

##schema flowhart

require("ggpubr")
require("grImport") 
PostScriptTrace("input_files/flowchart_evolve14.eps",charpath=FALSE)
f5a_flowchart_arial<-readPicture("flowchart_evolve14.eps.xml")

pdf("f45_flowchart_arial14.pdf")
grid.picture(f5a_flowchart_arial, gp=gpar(fontfamily="ArialMT"))
dev.off()
f5a_flowchart<-pictureGrob(f5a_flowchart_arial) 

###requirements to analyze ORFs individually
###good homology

borfs<-read.csv("orfs_scer_blast_matched_validated",sep=" ",stringsAsFactors=F) #ALL ATG-stop ORFs in yeast, including homology data using synteny or blast

#an alignment between S. cereisiae and another saccharomyces yeast is considered good homology if it is in synteny or the sequence itself is determined to be homologous at p-value <.01
blast_good_homology<-array(dim=c(length(orfs[,1]),7))
for(i in 1:7)
{
	blast_good_homology[,i]<-(borfs[,118+i]<.01 & (borfs[,178+i]==0|borfs[,186+i]==1))
}

#reading frame conservation calcuated for each ORF among all species for which it has good homology
spp_included<-1:7
blast_frame_conservation<-array()
for(i in 1:length(orfs[,1]))
{
	blast_frame_conservation[i]<-mean(as.numeric(borfs[i,(23:29)[spp_included]][which(blast_good_homology[i,spp_included])]))/(orfs$length[i])
}

homology_pass_index<-which(rowSums(blast_good_homology)>3) #index of ORFs with at least 4 species with good homology among saccharomyces, required for being in the high information group

homology_intergene_pass<-intersect(homology_pass_index,intergene_index)
homology_antisense_pass<-intersect(homology_pass_index,antisense_index)

homology_intergene_trans<-(intersect(homology_intergene_pass,candidate_trans))
homology_antisense_trans<-(intersect(homology_antisense_pass,candidate_trans))
length(homology_intergene_trans)

homology_canonical_trans<-intersect(canonical_trans,homology_pass_index)
homology_noncanonical_trans<-intersect(noncanonical_trans,homology_pass_index)

homology_canonical_intergene_trans<-intersect(homology_canonical_trans,intergene_index)
homology_noncanonical_intergene_trans<-intersect(homology_noncanonical_trans,intergene_index)

intergene_trans<-intersect(intergene_index,trans_index)
bad_homology_trans<-intersect(intersect(which(rowSums(blast_good_homology)<4),candidate_trans),intergene_trans)

mult<-read.csv("orfs_mult_align",sep=" ",stringsAsFactors = F) #multiple alignment of homologs for each ORF

#here we calculate nucleotide differences from the ORF among its homologs. Too few differences means insufficient power for using variation to infer selection
seqs<-list() 
for(i in 1:8)
{
  seqs[[i]]<-list()
}

for(i in 1:length(mult[,1]))
{
  for(j in 1:8)
  {
    seqs[[j]][[i]]<-strsplit(mult[i,j+1],split="")[[1]]
  }
}

diverse<-rep(0,length(seqs[[1]]))
alldiffs<-array()
meddiffs<-array()
for(i in 1:length(seqs[[1]]))
{
	diffs<-0
	pos<-0
	numdiffs<-array()
	j<-1
	for(k in 2:length(seqs))
	{
		if(blast_good_homology[i,k-1])
		{
			diffs<-diffs+length(which(seqs[[j]][[i]]!=seqs[[k]][[i]] & seqs[[j]][[i]]!="-" &  seqs[[k]][[i]]!="-"))
			pos<-pos+length(which(seqs[[j]][[i]]!="-" & seqs[[k]][[i]]!="-"))
			numdiffs[k]<-length(which(seqs[[j]][[i]]!=seqs[[k]][[i]] & seqs[[j]][[i]]!="-" &  seqs[[k]][[i]]!="-"))
		}
	}
	diverse[i]<-diffs/pos
	alldiffs[i]<-diffs
	meddiffs[i]<-median(numdiffs,na.rm=T)
}

diverse_index<-which(meddiffs>20) #need more than 20 median differences to be considered sufficiently diverse for high info group
diverse_homology_intergene_pass<-intersect(diverse_index,homology_intergene_pass)
diverse_intergene_homology_trans<-intersect(homology_intergene_trans,diverse_index)
diverse_antisense_homology_trans<-intersect(homology_antisense_trans,diverse_index)

nondiverse_index<-which(meddiffs<=20)

###SET of informative ORFs to study using RFC, having both sufficient homology and informative variation
length(diverse_intergene_homology_trans) #intergenic informative ORFs
length(diverse_antisense_homology_trans) #antisense informative ORFs
informative_orfs<-c(diverse_intergene_homology_trans,diverse_antisense_homology_trans)
##information describing plot
frame_conservation<-blast_frame_conservation
unconserved_index<-which(frame_conservation<=.6)
intermediate_index<-which(frame_conservation>.6 & frame_conservation<=.8)
conserved_index<-which(frame_conservation>.8)

class_unconserved<-intersect(diverse_intergene_homology_trans,unconserved_index)
class_intermediate<-intersect(diverse_intergene_homology_trans,intermediate_index)
class_conserved<-intersect(diverse_intergene_homology_trans,conserved_index)

class_unconserved_antisense<-intersect(diverse_antisense_homology_trans,unconserved_index)

length(class_unconserved)
length(class_intermediate)
length(class_conserved)

length(class_conserved)/length(diverse_intergene_homology_trans)
length(class_intermediate)/length(diverse_intergene_homology_trans)
length(class_unconserved)/length(diverse_intergene_homology_trans)

class_conserved_canonical<-intersect(canonical_index,class_conserved)
class_conserved_noncanonical<-intersect(noncanonical_index,class_conserved)

length(class_conserved_canonical)/length(intersect(canonical_index,diverse_intergene_homology_trans))
length(intersect(noncanonical_index,class_unconserved))/length(intersect(noncanonical_index,diverse_intergene_homology_trans))


length(class_unconserved_antisense)/length(diverse_antisense_homology_trans)
##the plot itself

canonical_diverse_intergene_homology_trans<-intersect(diverse_intergene_homology_trans,canonical_index)
noncanonical_diverse_intergene_homology_trans<-intersect(diverse_intergene_homology_trans,noncanonical_index)

canonical_diverse_antisense_homology_trans<-intersect(diverse_antisense_homology_trans,canonical_index)
noncanonical_diverse_antisense_homology_trans<-intersect(diverse_antisense_homology_trans,noncanonical_index)

allcats<-array()
allcats[1:length(canonical_diverse_intergene_homology_trans)]<-"Canonical"
allcats[length(allcats)+1:length(noncanonical_diverse_intergene_homology_trans)]<-"Noncanonical"

all_in_order<-c(canonical_diverse_intergene_homology_trans,noncanonical_diverse_intergene_homology_trans)
evol_data<-data.frame(cons=frame_conservation[all_in_order],category=allcats,diffs=meddiffs[all_in_order])
evol_data$category<-factor(evol_data$category,levels=c("Canonical","Noncanonical"))

dummy<-evol_data[which(evol_data$cons>.99 & evol_data$category=="Canonical"),]
full_evol_data<-evol_data
evol_data<-evol_data[which(evol_data$cons<.99|evol_data$category!="Canonical"),]
evol_data<-rbind(evol_data,dummy[1:700,])

require("ggplot2")
png("conservation_dichotomy_bottom.png",height=480,width=600)
f5b_bottom<-ggplot(evol_data[evol_data$cons<1.95,], aes(cons, fill = category)) +
  scale_fill_manual(values=c("red","blue"))+
  geom_histogram(binwidth = .05)+
  theme_classic(base_size=30) +
  theme(axis.title.y=element_blank(),legend.title=element_blank(),strip.background = element_blank(),strip.text.y=element_blank(),legend.position = c(0.3, 0.9))+
  labs(x=expression(paste("Frame conservation of ",italic("S. cerevisiae")," ORF with ",italic("sensu stricto"))),y="")+
  coord_cartesian(ylim=c(0,1000))+
  geom_vline(xintercept=.6, linetype="dashed", color="black", size=1.0)+
  geom_vline(xintercept=.8, linetype="dashed", color="black", size=1.0)
f5b_bottom
dev.off()
        
png("conservation_dichotomy_top.png",height=480,width=600)
f5b_top<-ggplot(full_evol_data[full_evol_data$cons<1.95,], aes(cons, fill = category)) +
  scale_fill_manual(values=c("red","blue"))+
  geom_histogram(binwidth = .05)+
  theme_classic(base_size=30) +
  guides(fill = guide_legend(override.aes = list(fill = c("white","white"))))+ 
  theme(axis.title.y=element_blank(),legend.title=element_blank(),legend.text = element_text(color = "white"), strip.background = element_blank(),strip.text.y=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),legend.position = c(0.3, 0.9))+ 
  labs(x=expression(paste("Frame conservation of ",italic("S. cer")," ORF with ",italic("sensu stricto"))),y="")+
  coord_cartesian(ylim=c(3600,4200))+
  scale_y_continuous(breaks=seq(3500,4000,250))+
  geom_vline(xintercept=.6, linetype="dashed", color="black", size=1.0)+
  geom_vline(xintercept=.8, linetype="dashed", color="black", size=1.0)
f5b_top
dev.off()

require("ggpubr")
png("conservation_dichotomy_combo.png",height=600,width=600)
f5b<-ggarrange(f5b_top, f5b_bottom, 
                        nrow = 2, ncol = 1,heights=c(1,2))
f5b<-annotate_figure(f5b, left = text_grob("Count", rot = 90,size=30))
f5b
dev.off()

###nORF new genes
##newly identified conserved genes analyis

class_conserved_noncanonical<-class_conserved_noncanonical[which(overlaps$ltr[class_conserved_noncanonical]==0)]
newgenes<-intersect(class_conserved_noncanonical,which(orfs$is_gene!="YAR073W" & orfs$orf_class!="pseudogene")) #which(gene_sim$best_evalue>.05 & orfs$orf_class!="pseudogene"))

tree<-read.csv("tree_dnds",sep=" ")

dnds_pval<-array()
dnds_omega<-(tree[newgenes,]$nonsyn_obs/tree[newgenes,]$nonsyn_exp)/(tree[newgenes,]$syn_obs/tree[newgenes,]$syn_exp)

for(index in 1:length(newgenes))
{
	if(tree[newgenes[index],]$syn_obs+tree[newgenes[index],]$nonsyn_obs>0)
	{
		dnds_pval[index]<-binom.test(tree[newgenes[index],]$nonsyn_obs, tree[newgenes[index],]$syn_obs+tree[newgenes[index],]$nonsyn_obs,tree[newgenes[index],]$nonsyn_exp/(tree[newgenes[index],]$nonsyn_exp+tree[newgenes[index],]$syn_exp))[[3]]
	}
}

#BLASTP analyis
pres<-read.csv("phylo_pres",sep=" ")
pres<-pres[,2:length(pres[1,])]
min_pres<-apply(pres,1,min)
num_pres<-apply(pres,1,function(x) length(which(x<10^-4)))

#tBLASTN analyis
prestb<-read.csv("phylo_pres_tblastn_scrambled",sep=" ")
prestb<-prestb[,2:length(prestb[1,])]
min_prestb<-apply(prestb,1,min)
num_prestb<-apply(prestb,1,function(x) length(which(x<10^-4)))

#S cer self-blast
gene_sim<-read.csv("gene_sim",sep=" ")

#pnps analysis 
porfs<-read.csv("orfs_pop",sep=" ")

bin_pvals<-array()
purify<-array()
omega<-array()
for(i in 1:length(porfs[,1]))
{
  if(porfs$syn_vars[i]+porfs$nonsyn_vars[i]>0)
  {
    bin_pvals[i]<-binom.test(porfs$nonsyn_vars[i],porfs$nonsyn_vars[i]+porfs$syn_vars[i],porfs$nonsyn_exp[i]/(porfs$nonsyn_exp[i]+porfs$syn_exp[i]))[[3]]
  }
  purify[i]<-porfs$nonsyn_vars[i]/(porfs$nonsyn_vars[i]+porfs$syn_vars[i])>porfs$nonsyn_exp[i]/(porfs$nonsyn_exp[i]+porfs$syn_exp[i])
  omega[i]<-(porfs$nonsyn_vars[i]/porfs$nonsyn_exp[i])/(porfs$syn_vars[i]/porfs$syn_exp[i])
}

df_newgenes<-data.frame(index=newgenes,chr=orfs[newgenes,]$contig+1,start=orfs[newgenes,]$start+1,end=orfs[newgenes,]$end+1,strand=orfs[newgenes,]$strand,blastp=min_pres[newgenes],tblastn=min_prestb[newgenes],orf_class=orfs$orf_class[newgenes],rfc=frame_conservation[newgenes],nuclenth=orfs$length[newgenes],gene_sim=gene_sim$best_evalue[newgenes])
df_newgenes$coord_str<-paste("chr",as.roman(df_newgenes$chr),":",df_newgenes$start,"-",df_newgenes$end,sep="")
translation_percentile<-array()

for(i in 1:length(newgenes))
{
	translation_percentile[i]<-length(which(ribo$frame0[noncanonical_trans]<ribo$frame0[newgenes[i]]))/length(noncanonical_trans)
}
df_newgenes$tran_per<-translation_percentile
df_newgenes$dnds<-dnds_omega
df_newgenes$dnds_pval<-dnds_pval
df_newgenes$pnps<-omega[newgenes]
df_newgenes$pnps_pval<-bin_pvals[newgenes]

write.table(df_newgenes,"newgenes_table",row.names=F,quote=F,sep="\t")

#### ANALYZING LOW-INFORMATION ORFS ####

#####PAIRWISE DN/DS ANALYSIS AGAINST S. PARADOXUS #####


##low info analysis groups 
low_info<-trans_index[which(!(trans_index %in% diverse_intergene_homology_trans))]
blast_nonoverlap_index<-all_index[which((overlaps$sense_ver+overlaps$sense_unchar+overlaps$sense_te+overlaps$anti_ver+overlaps$anti_unchar+overlaps$anti_te+overlaps$ars + overlaps$tRNA + overlaps$sno + overlaps$sn + overlaps$ltr+overlaps$anti_blocked +overlaps$sense_blocked )[all_index]==0)]
low_info_intergene<-intersect(low_info,intergene_index)
low_info_antisense<-intersect(low_info,antisense_index)
low_info_antisense1<-low_info_antisense[which(orfs$gene_rel_frame[low_info_antisense]==-1)]
low_info_antisense2<-low_info_antisense[which(orfs$gene_rel_frame[low_info_antisense]==-2)]
low_info_antisense3<-low_info_antisense[which(orfs$gene_rel_frame[low_info_antisense]==-3)]

no_ancient_homologs<-which(min_prestb>10^-4)
low_info_intergene_young<-intersect(low_info_intergene,no_ancient_homologs) #we analyze these in this section
low_info_can<-intersect(low_info_intergene_young,canonical_index)
low_info_noncan<-intersect(low_info_intergene_young,noncanonical_index)
low_info_noncan_blast<-intersect(low_info_noncan,blast_nonoverlap_index)

low_info_uorf_blast<-intersect(low_info_noncan_blast,uORF_index)
low_info_dorf_blast<-intersect(low_info_noncan_blast,dORF_index)
low_info_norf_blast<-intersect(low_info_noncan_blast,nORF_index)

low_info_can<-intersect(low_info_intergene_young,canonical_index)
low_info_noncan<-intersect(low_info_intergene_young,noncanonical_index)

low_info_can_blast<-intersect(low_info_can,blast_nonoverlap_index)
low_info_noncan_blast<-intersect(low_info_noncan,blast_nonoverlap_index)
can_orphans<-intersect(canonical_intergene_index,no_ancient_homologs)


#scrambled controls
scrambledtb<-read.csv("phylo_pres_tblastn_scrambled",sep=" ")
num_pres_scrambled<-apply(scrambledtb,1,function(x) length(which(x>10^-4 & x<10^-2)))
num_prestb_weak<-apply(prestb,1,function(x) length(which(x>10^-4 & x<10^-2)))

####functions

spar_dnds<-read.csv("spar_dnds_all_",sep=" ")

get_dnds_omega3<-function(x,s0=spar_dnds)
{
	orf_subset<-x
	s0<-s0[which(s0$orf_id %in% orfs$orf_id[orf_subset]),]
	T<-table(s0$context[which(s0$syn==1)])
	nucs<-c("A","C","G","T")
	syn_exp<-0
	nonsyn_exp<-0
	syn_obs<-0
	nonsyn_obs<-0
	sum_syn_base<-0
	sum_nonsyn_base<-0
	syn_diff<-array()
	for(i in 1:length(T))
	{
		for(j in 1:4)
		{
			nonsyn_base<-length(which(s0$syn==0 & s0$context==names(T)[i] & s0$sub==nucs[j]))
			nonsyn_sub<-length(which(s0$syn==0 & s0$context==names(T)[i] & s0$sub==nucs[j] & s0$actual==1))
			syn_base<-length(which(s0$syn==1 & s0$context==names(T)[i] & s0$sub==nucs[j]))
			syn_sub<-length(which(s0$syn==1 & s0$context==names(T)[i] & s0$sub==nucs[j] & s0$actual==1))
			if(syn_base>0 & nonsyn_base>0)
			{
				sum_syn_base<-sum_syn_base+syn_base
				sum_nonsyn_base<-sum_nonsyn_base+nonsyn_base
				rate<-(syn_sub+nonsyn_sub)/(syn_base+nonsyn_base)
				syn_exp<-syn_exp+rate*syn_base
				nonsyn_exp<-nonsyn_exp+rate*nonsyn_base
				syn_obs<-syn_obs+syn_sub
				nonsyn_obs<-nonsyn_obs+nonsyn_sub
				syn_diff[i]<-abs(rate*syn_base-syn_sub)+abs(rate*nonsyn_base-nonsyn_sub)
			}
		}
	}
	sample_syn_count<-rbinom(1000,sum_syn_base,syn_obs/sum_syn_base)
	sample_nonsyn_count<-rbinom(1000,sum_nonsyn_base,nonsyn_obs/sum_nonsyn_base)
	sample_omega<-(sample_nonsyn_count/nonsyn_exp)/(sample_syn_count/syn_exp)
	
	df_dnds<-data.frame(omega=(nonsyn_obs/nonsyn_exp)/(syn_obs/syn_exp),
		pval=binom.test(nonsyn_obs,nonsyn_obs+syn_obs,nonsyn_exp/(syn_exp+nonsyn_exp))[[3]],
		syn_exp=syn_exp,nonsyn_exp=nonsyn_exp,syn_obs=syn_obs,nonsyn_obs=nonsyn_obs,
		sum_syn_base=sum_syn_base,sum_nonsyn_base=sum_nonsyn_base,se=sd(sample_omega),ci975=quantile(sample_omega,.975),ci025=quantile(sample_omega,.025))
	return(df_dnds)
}

get_dnds_omega4<-function(x,s0=spar_dnds)
{
	orf_subset<-x
	s0<-s0[which(s0$orf_id %in% orfs$orf_id[orf_subset]),]
	T<-table(s0$context[which(s0$syn==1)])
	nucs<-c("A","C","G","T")
	syn_exp<-0
	nonsyn_exp<-0
	syn_obs<-0
	nonsyn_obs<-0
	sum_syn_base<-0
	sum_nonsyn_base<-0
	syn_diff<-array()
		#print(i)
	for(i in 1:4)
	{
		for(j in 1:4)
		{
			nonsyn_base<-length(which(s0$syn==0  & s0$sub==nucs[j] & s0$base==nucs[i]))
			nonsyn_sub<-length(which(s0$syn==0  & s0$sub==nucs[j] & s0$base==nucs[i] & s0$actual==1))
			syn_base<-length(which(s0$syn==1 & s0$sub==nucs[j] & s0$base==nucs[i]))
			syn_sub<-length(which(s0$syn==1  & s0$sub==nucs[j] &  s0$base==nucs[i] & s0$actual==1))
			if(syn_base>0 & nonsyn_base>0)
			{
				sum_syn_base<-sum_syn_base+syn_base
				sum_nonsyn_base<-sum_nonsyn_base+nonsyn_base
				rate<-(syn_sub+nonsyn_sub)/(syn_base+nonsyn_base)
				syn_exp<-syn_exp+rate*syn_base
				nonsyn_exp<-nonsyn_exp+rate*nonsyn_base
				syn_obs<-syn_obs+syn_sub
				nonsyn_obs<-nonsyn_obs+nonsyn_sub
				syn_diff[i]<-abs(rate*syn_base-syn_sub)+abs(rate*nonsyn_base-nonsyn_sub)
			}
		}
	}
	sample_syn_count<-rbinom(1000,sum_syn_base,syn_obs/sum_syn_base)
	sample_nonsyn_count<-rbinom(1000,sum_nonsyn_base,nonsyn_obs/sum_nonsyn_base)
	sample_omega<-(sample_nonsyn_count/nonsyn_exp)/(sample_syn_count/syn_exp)
	
	df_dnds<-data.frame(omega=(nonsyn_obs/nonsyn_exp)/(syn_obs/syn_exp),
		pval=binom.test(nonsyn_obs,nonsyn_obs+syn_obs,nonsyn_exp/(syn_exp+nonsyn_exp))[[3]],
		syn_exp=syn_exp,nonsyn_exp=nonsyn_exp,syn_obs=syn_obs,nonsyn_obs=nonsyn_obs,
		sum_syn_base=sum_syn_base,sum_nonsyn_base=sum_nonsyn_base,se=sd(sample_omega))
	return(df_dnds)
}

get_omega<-function(x) {
return( (sum(porfs$nonsyn_vars[x])/sum(porfs$nonsyn_exp[x]))/(sum(porfs$syn_vars[x])/sum(porfs$syn_exp[x])))
}
get_omega_pval<-function(x){
return(binom.test(sum(porfs$nonsyn_vars[x]),sum(porfs$nonsyn_vars[x])+sum(porfs$syn_vars[x]),sum(porfs$nonsyn_exp[x])/(sum(porfs$nonsyn_exp[x])+sum(porfs$syn_exp[x])))[[3]])
}

get_omega_se<-function(x,n=1000)
{
	sample_omegas<-array()
	for(i in 1:n)
	{
		x_sample<-sample(x,replace=TRUE)
		sample_omegas[i]<-get_omega(x_sample)
	}
	return(sd(sample_omegas))
}

combine_omega_ests<-function(x,y)
{
	combo_omega<-data.frame(omega=1,se=0)
	combo_omega$omega<-((x$nonsyn_obs+y$nonsyn_obs)/(x$nonsyn_exp+y$nonsyn_exp))/((x$syn_obs+y$syn_obs)/(x$syn_exp+y$syn_exp))

	sum_syn_base<-x$sum_syn_base+y$sum_syn_base
	syn_obs<-x$syn_obs+y$syn_obs
	sum_nonsyn_base<-x$sum_nonsyn_base+y$sum_nonsyn_base
	nonsyn_obs<-x$nonsyn_obs+y$nonsyn_obs
	nonsyn_exp<-x$nonsyn_exp+y$nonsyn_exp
	syn_exp<-x$syn_exp+y$syn_exp

	sample_syn_count<-rbinom(1000,sum_syn_base,syn_obs/sum_syn_base)
	sample_nonsyn_count<-rbinom(1000,sum_nonsyn_base,nonsyn_obs/sum_nonsyn_base)
	sample_omega<-(sample_nonsyn_count/nonsyn_exp)/(sample_syn_count/syn_exp)

	combo_omega$se<- sd(sample_omega)
	return(combo_omega)

}

####Analysis of selection in groups of ORFs grouped by decile of various properties 

####translation rate deciles

trates<-ribo$frame0[low_info_noncan]/orfs$length[low_info_noncan]
trates_blast<-ribo$frame0[low_info_noncan_blast]/orfs$length[low_info_noncan_blast]
trate_quantiles<-quantile(trates,seq(0,1,.1))
trate_groups<-list()
trate_groups_blast<-list()
for(i in 1:(length(trate_quantiles)-1))
{
	trate_groups[[i]]<-low_info_noncan[which(trates>=trate_quantiles[i] & trates<trate_quantiles[i+1])]
	trate_groups_blast[[i]]<-low_info_noncan_blast[which(trates_blast>=trate_quantiles[i] & trates_blast<trate_quantiles[i+1])]
}
trate_omegas<-array()
trate_omegas_se<-array()
trate_omegas_pvals<-array()
trate_weak<-array()
trate_weak_control<-array()
trate_counts<-array()
trate_dnds<-array()
trate_dnds_pvals<-array()
trate_dnds_se<-array()

for(i in 1:(length(trate_quantiles)-1))
{
	dnds_obj<-get_dnds_omega3(trate_groups[[i]])
	trate_dnds[i]<-dnds_obj$omega
	trate_dnds_pvals[i]<-dnds_obj$pval
	trate_dnds_se[i]<-dnds_obj$se
	trate_omegas[i]<-get_omega(trate_groups[[i]])
	trate_omegas_se[i]<-get_omega_se(trate_groups[[i]])
	trate_omegas_pvals[i]<-get_omega_pval(trate_groups[[i]])
	trate_weak[i]<-length(which(num_prestb_weak[trate_groups_blast[[i]]]>1))
	trate_weak_control[i]<-length(which(num_pres_scrambled[trate_groups_blast[[i]]]>1))	
	trate_counts[i]<-length(trate_groups_blast[[i]])
}

###coding scores groups

cs<-read.csv("coding_scores_",sep= " ")#coding scores for ORFs
cscores<-cs$score[low_info_noncan]
cscores_blast<-cs$score[low_info_noncan_blast]

cscores_quantiles<-quantile(cscores,seq(0,1,.1))
cscores_groups<-list()
cscores_groups_blast<-list()

for(i in 1:(length(cscores_quantiles)-1))
{
	cscores_groups[[i]]<-low_info_noncan[which(cscores>=cscores_quantiles[i] & cscores<cscores_quantiles[i+1])]
	cscores_groups_blast[[i]]<-low_info_noncan_blast[which(cscores_blast>=cscores_quantiles[i] & cscores_blast<cscores_quantiles[i+1])]
}
cscores_omegas<-array()
cscores_omegas_se<-array()
cscores_omegas_pvals<-array()

cscores_weak<-array()
cscores_weak_control<-array()
cscores_counts<-array()

cscores_dnds<-array()
cscores_dnds_pvals<-array()
cscores_dnds_se<-array()
for(i in 1:(length(cscores_quantiles)-1))
{
	dnds_obj<-get_dnds_omega3(cscores_groups[[i]])
	cscores_omegas[i]<-get_omega(cscores_groups[[i]])
	cscores_omegas_se[i]<-get_omega_se(cscores_groups[[i]])
	cscores_omegas_pvals[i]<-get_omega_pval(cscores_groups[[i]])
	cscores_weak[i]<-length(which(num_prestb_weak[cscores_groups_blast[[i]]]>1))
	cscores_weak_control[i]<-length(which(num_pres_scrambled[cscores_groups_blast[[i]]]>1))	
	cscores_counts[i]<-length(cscores_groups_blast[[i]])
	cscores_dnds[i]<-dnds_obj$omega #get_dnds_omega2(cscores_groups[[i]])
	cscores_dnds_se[i]<-dnds_obj$se
	cscores_dnds_pvals[i]<-dnds_obj$pval#get_dnds_omega2_pval(cscores_groups[[i]])
}

####length groups
lens<-orfs$length[low_info_noncan]

lens_blast<-orfs$length[low_info_noncan_blast]

lens_quantiles<-quantile(lens,seq(0,1,.1))
lens_groups<-list()
lens_groups_blast<-list()

for(i in 1:(length(lens_quantiles)-1))
{
	lens_groups[[i]]<-low_info_noncan[which(lens>=lens_quantiles[i] & lens<lens_quantiles[i+1])]
	lens_groups_blast[[i]]<-low_info_noncan_blast[which(lens_blast>=lens_quantiles[i] & lens_blast<lens_quantiles[i+1])]
}

lens_omegas<-array()
lens_omegas_se<-array()
lens_omegas_pvals<-array()
lens_weak<-array()
lens_weak_control<-array()
lens_counts<-array()
lens_dnds<-array()
lens_dnds_pvals<-array()
lens_dnds_se<-array()

for(i in 1:(length(lens_quantiles)-1))
{
	dnds_obj<-get_dnds_omega3(lens_groups[[i]])
	lens_omegas[i]<-get_omega(lens_groups[[i]])
	lens_omegas_se[i]<-get_omega_se(lens_groups[[i]])
	lens_omegas_pvals[i]<-get_omega_pval(lens_groups[[i]])
	lens_weak[i]<-length(which(num_prestb_weak[lens_groups_blast[[i]]]>1))
	lens_weak_control[i]<-length(which(num_pres_scrambled[lens_groups_blast[[i]]]>1))	
	lens_counts[i]<-length(lens_groups_blast[[i]])
	lens_dnds[i]<-dnds_obj$omega#get_dnds_omega2(lens_groups[[i]])
	lens_dnds_pvals[i]<-dnds_obj$pval#get_dnds_omega2_pval(lens_groups[[i]])
	lens_dnds_se[i]<-dnds_obj$se
}


low_info_uorf<-intersect(low_info_noncan,uORF_index)
low_info_dorf<-intersect(low_info_noncan,dORF_index)
low_info_norf<-intersect(low_info_noncan,nORF_index)


context_dnds<-array()
context_dnds[1]<-get_dnds_omega3(low_info_uorf)$omega
context_dnds[2]<-get_dnds_omega3(low_info_dorf)$omega
context_dnds[3]<-get_dnds_omega3(low_info_norf)$omega

context_dnds_pvals<-array()
context_dnds_pvals[1]<-get_dnds_omega3(low_info_uorf)$pval
context_dnds_pvals[2]<-get_dnds_omega3(low_info_dorf)$pval
context_dnds_pvals[3]<-get_dnds_omega3(low_info_norf)$pval

context_dnds_se<-array()
context_dnds_se[1]<-get_dnds_omega3(low_info_uorf)$se
context_dnds_se[2]<-get_dnds_omega3(low_info_dorf)$se
context_dnds_se[3]<-get_dnds_omega3(low_info_norf)$se

context_omegas<-array()
context_omegas_se<-array()
context_omegas_pvals<-array()

context_omegas[1]<-get_omega(low_info_uorf)
context_omegas[2]<-get_omega(low_info_dorf)
context_omegas[3]<-get_omega(low_info_norf)

context_omegas_se[1]<-get_omega_se(low_info_uorf)
context_omegas_se[2]<-get_omega_se(low_info_dorf)
context_omegas_se[3]<-get_omega_se(low_info_norf)

context_omegas_pvals[1]<-get_omega_pval(low_info_uorf)
context_omegas_pvals[2]<-get_omega_pval(low_info_dorf)
context_omegas_pvals[3]<-get_omega_pval(low_info_norf)

context_weak<-array()
context_weak[1]<-length(which(num_prestb_weak[low_info_uorf_blast]>1))
context_weak[2]<-length(which(num_prestb_weak[low_info_dorf_blast]>1))
context_weak[3]<-length(which(num_prestb_weak[low_info_norf_blast]>1))

context_weak_control<-array()
context_weak_control[1]<-length(which(num_pres_scrambled[low_info_uorf_blast]>1))	
context_weak_control[2]<-length(which(num_pres_scrambled[low_info_dorf_blast]>1))	
context_weak_control[3]<-length(which(num_pres_scrambled[low_info_norf_blast]>1))	

context_counts<-array()
context_counts[1]<-length(low_info_uorf_blast)
context_counts[2]<-length(low_info_dorf_blast)
context_counts[3]<-length(low_info_norf_blast)

canonical_dnds<-array()
canonical_dnds[1]<-get_dnds_omega3(low_info_can)$omega
canonical_dnds[2]<-get_dnds_omega3(canonical_trans)$omega

canonical_dnds_pvals<-array()
canonical_dnds_pvals[1]<-get_dnds_omega3(low_info_can)$pval
canonical_dnds_pvals[2]<-get_dnds_omega3(canonical_trans)$pval

canonical_dnds_se<-array()
canonical_dnds_se[1]<-get_dnds_omega3(low_info_can)$se
canonical_dnds_se[2]<-get_dnds_omega3(canonical_trans)$se

canonical_omega_se<-array()
canonical_omega<-array()
canonical_omega[2]<-get_omega(canonical_trans)
canonical_omega_se[2]<-get_omega_se(canonical_trans)
canonical_omega[1]<-get_omega(low_info_can)
canonical_omega_se[1]<-get_omega_se(low_info_can)

low_info_corf_weak<-length(which(num_prestb_weak[low_info_can_blast]>0))
low_info_corf_weak_control<-length(which(num_pres_scrambled[low_info_can_blast]>0))
low_info_corf_counts<-length(low_info_can_blast)

can_orphan_weak<-length(which(num_prestb_weak[can_orphans]>1))
can_orphan_weak_control<-length(which(num_pres_scrambled[can_orphans]>1))
can_orphan_counts<-length(can_orphans)

noncan_weak<-length(which(num_prestb_weak[low_info_noncan_blast]>1))
noncan_weak_control<-length(which(num_pres_scrambled[low_info_noncan_blast]>1))
noncan_weak_counts<-length(low_info_noncan_blast)


antisyn1<-read.csv("spar_antisyn_all_1",sep=" ")
antisyn2<-read.csv("spar_antisyn_all_2",sep=" ")
antisyn3<-read.csv("spar_antisyn_all_3",sep=" ")

noncanonical_intergene_trans<-intersect(noncanonical_intergene_index,trans_index)
noncanonical_antisense_trans<-intersect(noncanonical_antisense_index,trans_index)

anti1<-noncanonical_antisense_trans[which(orfs$gene_rel_frame[noncanonical_antisense_trans]==-1)]
anti2<-noncanonical_antisense_trans[which(orfs$gene_rel_frame[noncanonical_antisense_trans]==-2)]
anti3<-noncanonical_antisense_trans[which(orfs$gene_rel_frame[noncanonical_antisense_trans]==-3)]

can_antisense1<-intersect(canonical_trans,antisense_index)[which(orfs$gene_rel_frame[intersect(canonical_trans,antisense_index)]==-1)]
can_antisense2<-intersect(canonical_trans,antisense_index)[which(orfs$gene_rel_frame[intersect(canonical_trans,antisense_index)]==-2)]
can_antisense3<-intersect(canonical_trans,antisense_index)[which(orfs$gene_rel_frame[intersect(canonical_trans,antisense_index)]==-3)]

anti_dnds<-array()
anti_dnds_pvals<-array()
anti_dnds_se<-array()

anti2_dnds<-get_dnds_omega3(anti2,antisyn2)
anti3_dnds<-get_dnds_omega3(anti3,antisyn3)

anti_dnds[1]<-combine_omega_ests(anti2_dnds,anti3_dnds)$omegaanti_dnds_se[1]<- combine_omega_ests(anti2_dnds,anti3_dnds)$se

can_anti2_dnds<-get_dnds_omega3(can_antisense2,antisyn2)
can_anti3_dnds<-get_dnds_omega3(can_antisense3,antisyn3)
anti_dnds[2]<-combine_omega_ests(can_anti2_dnds,can_anti3_dnds)$omegaanti_dnds_se[2]<- combine_omega_ests(can_anti2_dnds,can_anti3_dnds)$se

antipop1<-read.csv("anti_pnps_pop_1",sep=" ")
antipop2<-read.csv("anti_pnps_pop_2",sep=" ")
antipop3<-read.csv("anti_pnps_pop_3",sep=" ")

anti_pnps<-array()
anti_pnps_pvals<-array()
anti_pnps_se<-array()

anti2_pnps<-get_dnds_omega4(anti2,antipop2)
anti3_pnps<-get_dnds_omega4(anti3,antipop3)
anti_pnps[1]<-combine_omega_ests(anti2_pnps,anti3_pnps)$omega
anti_pnps_se[1]<-combine_omega_ests(anti2_pnps,anti3_pnps)$se

can_anti2_pnps<-get_dnds_omega4(can_antisense2,antipop2)
can_anti3_pnps<-get_dnds_omega4(can_antisense3,antipop3)

anti_pnps[2]<-combine_omega_ests(can_anti2_pnps,can_anti3_pnps)$omega
anti_pnps_se[2]<-combine_omega_ests(can_anti2_pnps,can_anti3_pnps)$se


# anti_dnds[1]<-get_dnds_omega3(anti2,antisyn2)$omega
# anti_dnds[2]<-get_dnds_omega3(anti3,antisyn3)$omega
# anti_dnds[3]<-get_dnds_omega3(can_antisense2,antisyn2)$omega
# anti_dnds[4]<-get_dnds_omega3(can_antisense3,antisyn3)$omega

# anti_dnds_pvals[1]<-get_dnds_omega3(anti2,antisyn2)$pval
# anti_dnds_pvals[2]<-get_dnds_omega3(anti3,antisyn3)$pval
# anti_dnds_pvals[3]<-get_dnds_omega3(can_antisense2,antisyn2)$pval
# anti_dnds_pvals[4]<-get_dnds_omega3(can_antisense3,antisyn3)$pval

# anti_dnds_se[1]<-get_dnds_omega3(anti2,antisyn2)$se
# anti_dnds_se[2]<-get_dnds_omega3(anti3,antisyn3)$se
# anti_dnds_se[3]<-get_dnds_omega3(can_antisense2,antisyn2)$se
# anti_dnds_se[4]<-get_dnds_omega3(can_antisense3,antisyn3)$se




trate_labels<-paste("t",1:10,sep="")#as.character(round(trate_mean,2))
cscores_labels<-paste("c",1:10,sep="")#as.character(round(cscores_mean,2))
lens_labels<-paste("l",1:10,sep="")#as.character(round(lens_mean,2))
context_labels<-c("uORF","dORF","Independent")


all_labels<-(c(trate_labels,"x1",cscores_labels,"x2",lens_labels,"x3",context_labels,"x4","Low info cORF","All nonoverlap cORF","x5","Antisense nORF", "Antisense cORF"))


grouptype<-array()
grouptype[1:11]<-"translation"
grouptype[12:22]<-"cscore"
grouptype[23:33]<-"length"
grouptype[34:37]<-"context"
grouptype[38:40]<-"conserved"
grouptype[41:42]<-"antisense"

df_trate_dnds<-data.frame(dnds=c(trate_dnds,0,cscores_dnds,0,lens_dnds,0,context_dnds,0,canonical_dnds,0,anti_dnds),
	labs=all_labels,
	group=grouptype,
	dnds_se=c(trate_dnds_se,0,cscores_dnds_se,0,lens_dnds_se,0,context_dnds_se,0,canonical_dnds_se,0,anti_dnds_se))
df_trate_dnds$labs<-factor(df_trate_dnds$labs,levels=all_labels)

#figure 4b: dn/ds
require("ggplot2")
png("dndsall.png",height=880,width=800)
f5d<-ggplot(df_trate_dnds, aes(x=labs,y=dnds)) +
  geom_bar(stat="identity",aes(fill=group)) +
  theme_classic(base_size=30) +
  theme(legend.position = "none",plot.margin=margin(c(5.5,5.5,5.5,50), "pt"),axis.title.y = element_blank()) +
  geom_errorbar(aes(x=labs,ymin=dnds-dnds_se,ymax=dnds+dnds_se)) +
  labs(x="",y="dN/dS")+
  annotate(geom = "text", x = 32+7, y = -.24, label = "Translation rate\ndecile", size = 8)+
  annotate(geom = "text", x = 21+7, y = -.24, label = "Coding score\ndecile", size = 8)+
  annotate(geom = "text", x = 10+7, y = -.24, label = "ORF length\ndecile", size = 8)+
  coord_flip(ylim=c(0,1.2),clip="off",expand = FALSE)+
  scale_y_continuous(breaks=c(.25,.5,.75,1))+
  scale_x_discrete(limits = rev(levels(df_trate_dnds$labs)), breaks = all_labels[-which(all_labels%in%c("x1","x2","x3","x4","x5") )] , labels=c(1:10,1:10,1:10,c("uORF","dORF","Independent","Low info cORF","All nonoverlap cORF","Antisense nORF", "Antisense cORF")) )+
  geom_hline(yintercept=1, linetype="dashed", color="black", size=1.0)
f5d
dev.off()


#######pnps


df_trate_pnps<-data.frame(pnps=c(trate_omegas,0,cscores_omegas,0,lens_omegas,0,context_omegas,0,canonical_omega,0,anti_pnps),
	labs=all_labels,
	group=grouptype,
	pnps_se=c(trate_omegas_se,0,cscores_omegas_se,0,lens_omegas_se,0,context_omegas_se,0,canonical_omega_se,0,anti_pnps_se))
df_trate_pnps$labs<-factor(df_trate_pnps$labs,levels=all_labels)


#figure 4c: pn/ps
require("ggplot2")
png("pnpsall.png",height=880,width=800)
f5c<-ggplot(df_trate_pnps, aes(x=labs,y=pnps)) +
  geom_bar(stat="identity",aes(fill=group)) +
  theme_classic(base_size=30) +
  theme(legend.position = "none",plot.margin=margin(c(5.5,5.5,5.5,50), "pt"),axis.title.y = element_blank()) +
  geom_errorbar(aes(x=labs,ymin=pnps-pnps_se,ymax=pnps+pnps_se)) +
  labs(x="",y="pN/pS")+
  annotate(geom = "text", x = 32+7, y = -.24, label = "Translation rate\ndecile", size = 8)+
  annotate(geom = "text", x = 21+7, y = -.24, label = "Coding score\ndecile", size = 8)+
  annotate(geom = "text", x = 10+7, y = -.24, label = "ORF length\ndecile", size = 8)+
  coord_flip(ylim=c(0,1.2),clip="off",expand = FALSE)+
  scale_y_continuous(breaks=c(.25,.5,.75,1))+
  scale_x_discrete(limits = rev(levels(df_trate_pnps$labs)), breaks = all_labels[-which(all_labels%in%c("x1","x2","x3","x4","x5") )] , labels=c(1:10,1:10,1:10,c("uORF","dORF","Independent","Low info cORF","All nonoverlap cORF","Antisense nORF", "Antisense cORF")) )+
  geom_hline(yintercept=1, linetype="dashed", color="black", size=1.0)
f5c
dev.off()

##figure 4d: weak

blast_grouptype<-array()
blast_grouptype[1:22]<-"translation"
blast_grouptype[23:44]<-"cscore"
blast_grouptype[45:66]<-"length"
blast_grouptype[67:74]<-"context"

is_control<-c(rep("no",length(trate_weak)),rep("yes",length(trate_weak_control)),"no","yes",rep("no",length(trate_weak)),rep("yes",length(trate_weak_control)),"no","yes",rep("no",length(trate_weak)),rep("yes",length(trate_weak_control)),"no","yes",rep("no",length(context_weak)),rep("yes",length(context_weak_control)),"no","yes")
blast_labs<-c(trate_labels,trate_labels,"x0","x0",cscores_labels,cscores_labels,"x1","x1",lens_labels,lens_labels,"x2","x2",context_labels,context_labels,"x3","x3")
blast_weak<-c(trate_weak_control,trate_weak,0,0,cscores_weak_control,cscores_weak,0,0,lens_weak_control,lens_weak,0,0,context_weak_control,context_weak,0,0)
blast_counts<-c(trate_counts,trate_counts,0,0,cscores_counts,cscores_counts,0,0,lens_counts,lens_counts,0,0,context_counts,context_counts,0,0)
df_blast<-data.frame(weak=blast_weak,is_control=is_control,labs=blast_labs,grouptype=blast_grouptype,counts=blast_counts)
df_blast$labs<-factor(df_blast$labs,levels=unique(blast_labs))

df_blast$weak_se<-sqrt((df_blast$weak/df_blast$counts)*(1-(df_blast$weak/df_blast$counts))/df_blast$counts)



require("ggplot2")
png("blast_trate.png",height=880,width=800)
f5e<-ggplot(df_blast, aes(x=labs,y=weak/counts, group=is_control , fill=blast_grouptype)) +
  geom_bar(stat="identity",width=0.6,position=position_dodge(width=.7),aes(alpha=is_control) ) +
  scale_alpha_discrete(range=c(.5,1),labels=c("Scrambled sequence","ORF sequence"))+
  theme_classic(base_size=30) +
  theme(legend.position = c(.6,.8),plot.margin=margin(c(5.5,5.5,5.5,50), "pt"),axis.title.y = element_blank()) +
  geom_errorbar(aes(x=labs,ymin=weak/counts-weak_se,ymax=weak/counts+weak_se),width=0.6,position=position_dodge(width=.7)) +
  labs(x="",y="Frequency of weak TBLASTN matches")+
  annotate(geom = "text", x = 32, y = -.01, label = "Translation rate\ndecile", size = 8)+
  annotate(geom = "text", x = 21, y = -.01, label = "Coding score\ndecile", size = 8)+
  annotate(geom = "text", x = 10, y = -.01, label = "ORF length\ndecile", size = 8)+
  coord_flip(ylim=c(0,.062),clip="off",expand = FALSE)+
  scale_x_discrete(limits = rev(levels(df_blast$labs)), breaks = unique(blast_labs[-which(blast_labs%in%c("x0","x1","x2","x3","x4") )]) , labels=c(1:10,1:10,1:10,c("uORF","dORF","Independent"))) +
  guides(fill=FALSE, alpha=guide_legend(title=NULL,reverse = TRUE,override.aes=list(fill="#C77CFF")))
f5e
dev.off()


df_blast_overall<-data.frame(weak=c(can_orphan_weak,can_orphan_weak_control),counts=c(can_orphan_counts,can_orphan_counts),labs=c("cORFs, no\nstrong homologs","Scrambled controls"))
df_blast_overall$weak_se<-sqrt((df_blast_overall$weak/df_blast_overall$counts)*(1-(df_blast_overall$weak/df_blast_overall$counts))/df_blast_overall$counts)

require("ggplot2")
png("blast_trate1.png",height=800,width=800)
f4f<-ggplot(df_blast_overall,aes(x=labs,y=weak/counts))+
  geom_bar(stat="identity",width=0.6,position=position_dodge(width=.7)) +
  theme_classic(base_size=30) +
  geom_errorbar(aes(x=labs,ymin=weak/counts-weak_se,ymax=weak/counts+weak_se),width=0.6,position=position_dodge(width=.7)) +
  labs(x="",y="Frequency of weak TBLASTN matches")+
  theme(axis.text=element_text(size=30),plot.margin=margin(c(5.5,5.5,-14,5.5), "pt"))#+
 # coord_flip()
f4f
dev.off()


###custom figure to show categories of annotated ORFs


transient_set<-c(intersect(class_unconserved,no_ancient_homologs),class_unconserved_antisense,intersect(low_info_intergene_young,noncanonical_index))

require("plotrix")
start_x<-0
start_y<-90
side<-(2)
gap<-0
between_gap<-65
between_gap2<-95
legend_x<-62
text_size<-2.4
legend_y<-60

png("FigCustom2.png",height=500,width=500*13/10)
fig_cust<-~{plot(1:130,(1:130)*10/13,type="n",xlab="",ylab="",xaxt='n',yaxt='n',bty="n")
text(0,99,"Transient",cex=text_size,adj=0)
text(0,99-6,"18,755 ORFs",cex=text_size,adj=0)

text(0+between_gap,99,"Conserved",cex=text_size,adj=0)
text(0+between_gap,99-6,"5,076 ORFs",cex=text_size,adj=0)

text(0+between_gap2,99,"Undetermined",cex=text_size,adj=0)
text(0+between_gap2,99-6,"486 ORFs",cex=text_size,adj=0)

rect(45+legend_x,legend_y,45+legend_x+side,legend_y-side,col="pink")			
text(60+legend_x,legend_y-1,"15 cORFs",cex=text_size)

rect(45+legend_x,legend_y-5,45+legend_x+side,legend_y-side-5,col="cyan")			
text(60+legend_x,legend_y-1-5,"15 nORFs",cex=text_size)


undetermined_set<-candidate_trans[which(!(candidate_trans %in% c(transient_set,class_conserved)))]
length(undetermined_set)

#transient
for(i in 1:45)
{
	for(j in 1:31)
	{	
		if(j==1&i<=7)
		{
			rect(start_x+(j-1)*side,start_y-(i-1)*side,start_x+j*side-gap,start_y-i*side+gap,col="pink")			
		}
		else if(j<31|i<=37)
		{
			rect(start_x+(j-1)*side,start_y-(i-1)*side,start_x+j*side-gap,start_y-i*side+gap,col="cyan")
		}
	}
}
#conserved
for(i in 1:45)
{
	for(j in 1:8)
	{	
		if(j==1&i<=1)
		{
			rect(between_gap+start_x+(j-1)*side,start_y-(i-1)*side,between_gap+start_x+j*side-gap,start_y-i*side+gap,col="cyan")			
		}
		else if(j<8|i<=23)
		{
			rect(between_gap+start_x+(j-1)*side,start_y-(i-1)*side,between_gap+start_x+j*side-gap,start_y-i*side+gap,col="pink")
		}
	}
}

#undetermined
for(i in 1:45)
{
	for(j in 1:1)
	{	
		if(j==1&i<=1)
		{
			rect(between_gap2+start_x+(j-1)*side,start_y-(i-1)*side,between_gap2+start_x+j*side-gap,start_y-i*side+gap,col="pink")			
		}
		else if(i<=32)
		{
			rect(between_gap2+start_x+(j-1)*side,start_y-(i-1)*side,between_gap2+start_x+j*side-gap,start_y-i*side+gap,col="cyan")
		}
	}
}
}
dev.off()





f5left<-plot_grid(f5a_flowchart,f5d,nrow=2,ncol=1,labels=c("A","D"),label_size=30)
f5mid<-plot_grid(f5b,f5e,nrow=2,ncol=1,rel_heights=c(1,1),labels=c("B","E"),label_size=30)
f5right<-plot_grid(f5c,fig_cust,nrow=2,ncol=1,rel_heights=c(1,1),labels=c("C","F"),label_size=30)

require("cowplot")
png("Figure5___.png",width=3000,height=2000)
plot_grid(f5left,f5mid,f5right,ncol=3,nrow=1)
dev.off()


####Supplementary table 6: ORF data

df_orf_data<-data.frame(orf_id=orfs$id, chromosome=orfs$contig+1, first_coord=orfs$start+1, last_coord=orfs$end+1, strand=orfs$strand) 
df_orf_data$strand[which(orfs$strand==0)]<-"+"
df_orf_data$strand[which(orfs$strand==1)]<-"-"
df_orf_data$gene_sysname<-orfs$is_gene
df_orf_data$orf_class<-orfs$orf_class
df_orf_data$length<-orfs$length
df_orf_data$rfc<-frame_conservation
df_orf_data$rfc[which(is.na(df_orf_data$rfc))]<-(-1)
df_orf_data$cs<-cs$score
df_orf_data$scer_gene_eval<-gene_sim$best_evalue
df_orf_data$tblasn_eval<-min_prestb

df_orf_data$is_intergenic<-0
df_orf_data$is_intergenic[intergene_index]<-1

df_orf_data$is_transient<-0
df_orf_data$is_transient[transient_set]<-1

df_orf_data$up_extend_gene<-orfs$up_extend_gene
df_orf_data$is_candidate<-rep(0,length(orfs[,1]))
df_orf_data$is_candidate[candidate_index]<-1

qvals<-array()
for(i in 1:length(orfs[,1]))
{
	qvals[i]<-length(which(ribo$scram_pval[noncanonical_index]<=ribo$pval[i]))/length(which(ribo$pval[noncanonical_index]<=ribo$pval[i]))
	if(!is.na(ribo$pval[i]) & length(which(ribo$scram_pval[noncanonical_index]<=ribo$pval[i]))==0 & length(which(ribo$pval[noncanonical_index]<=ribo$pval[i]))==0)
	{
		qvals[i]<-0
	}
}

df_orf_data$pval<-ribo$pval
df_orf_data$qvalue<-qvals
df_orf_data$reads<-ribo$frame0


# write.table(df_orf_data,"SupplementaryTable6_orf_data",row.names=F,quote=F,sep="\t")


ypd_end<-read.csv("riboseq_orfsypd_endpoint1",sep=" ")
ypd_nochx<-read.csv("riboseq_orfsypd_nochx_endpoint1",sep=" ")

df_orf_data$chx_plus_pval<-chx1$pval
df_orf_data$chx_minus_pval<-chx0$pval
df_orf_data$ypd_media_pval<-ypd_end$pval
df_orf_data$sd_media_pval<-ypd4[[20]]$pval
df_orf_data$ypd_chx_minus_pval<-ypd_nochx$pval

df_orf_data$chx_plus_reads<-chx1$frame0
df_orf_data$chx_minus_reads<-chx0$frame0
df_orf_data$ypd_media_reads<-ypd_end$frame0
df_orf_data$sd_media_reads<-ypd4[[20]]$frame0
df_orf_data$ypd_chx_minus_reads<-ypd_nochx$frame0

df_orf_data<-cbind(df_orf_data,overlaps[,8:25])

neigh<-read.csv("orf_neighbors",sep=" ",stringsAsFactors=F)
neigh[which(neigh$up_neighbor==""),]$up_neighbor<-"X"
neigh[which(neigh$down_neighbor==""),]$down_neighbor<-"X"
df_orf_data<-cbind(df_orf_data,neigh[,2:5])

df_orf_data$start_codon_overlap<-neigh$start_codon_overlap
df_orf_data$stop_codon_overlap<-neigh$stop_codon_overlap


#df_orf_data$up_dist<-df_orf_data$up_neighbor_dist
#df_orf_data$down_dist<-df_orf_data$down_neighbor_dist
#df_orf_data$up_dist[which(df_orf_data$strand=="-")]<-neigh$down_neighbor_dist[which(df_orf_data$strand=="-")]
#df_orf_data$down_dist[which(df_orf_data$strand=="-")]<-neigh$up_neighbor_dist[which(df_orf_data$strand=="-")]
df_orf_data$orf_coord_id<-paste("orf_chr_",df_orf_data$chr,"_",df_orf_data$first_coord,"_",df_orf_data$last_coord,"_",df_orf_data$strand,sep="")

df_orf_data<-cbind(df_orf_data,borfs[,39:45])

df_orf_data$informative<-0
df_orf_data$informative[informative_orfs]<-1

write.table(df_orf_data,"orf_data",row.names=F,quote=F,sep="\t",col.names=T)

#write.table(df_orf_data,"orf_data",row.names=F,quote=F,sep="\t",col.names=F)

##########################
#####PHENOTYPES###########
##########################
##########################

##### localization and protein detection

###Cyclops dataset
cyclops1<-read.csv("CYCLoPs Table3-1_WT1_LOCscore.csv",skip=3,stringsAsFactors=F)
cyclops2<-read.csv("CYCLoPs Table3-2_WT2_LOCscore.csv",skip=3,stringsAsFactors=F)
cyclops3<-read.csv("CYCLoPs Table3-3_WT3_LOCscore.csv",skip=3,stringsAsFactors=F)
cyclops_all<-read.csv("CYCLoPs Table1_GFPIntensityMatrix.csv",skip=4,stringsAsFactors=F)

named_transient_orfs<-transient_set[which(orfs$is_gene[transient_set]!="X")]
named_transient_orfs_intergenic<-intersect(named_transient_orfs,intergene_index)

cyclops_localization<-array()
for(i in 1:length(named_transient_orfs))
{
	sel<-which(cyclops1$ORF==orfs$is_gene[named_transient_orfs[i]])
	if(length(sel)>0)
	{
		cyclops_localization[i]<-names(which.max(cyclops1[sel,3:18]))
	}
	sel<-which(cyclops2$ORF==orfs$is_gene[named_transient_orfs[i]])
	if(length(sel)>0)
	{
		cyclops_localization[i]<-names(which.max(cyclops2[sel,3:18]))
	}
	sel<-which(cyclops3$ORF==orfs$is_gene[named_transient_orfs[i]])
	if(length(sel)>0)
	{
		cyclops_localization[i]<-names(which.max(cyclops3[sel,3:18]))
	}
}

if(length(cyclops_localization)<length(named_transient_orfs))
{
	cyclops_localization[length(named_transient_orfs)]<-NA
}

cswat_localization1<-array()
cswat_localization2<-array()
cswat_localization<-array()

localize<-read.csv("local.csv",stringsAsFactors=F)
for(i in 1:length(named_transient_orfs))
{
	sel<-which(localize$Systematic.Name==orfs$is_gene[named_transient_orfs[i]])
	if(length(sel)>0)
	{
		cswat_localization1[i]<-localize[sel,]$mNG.I
		cswat_localization2[i]<-localize[sel,]$mNG.II
		if(cswat_localization1[i] %in% c("cyto","mito"))
		{
			cswat_localization[i]<-cswat_localization1[i]
		}
		else
		{
			cswat_localization[i]<-cswat_localization2[i]
		}
	}
}

localization<-cyclops_localization
localization[which(cswat_localization=="cyto")]<-"Cytoplasm"
localization[which(cswat_localization=="mito")]<-"Mitochondria"


###detection using star database

star<-read.csv("star_results.csv",skip=2,stringsAsFactors=F)
starcount<-rep(0,length(named_transient_orfs))
for(i in 1:length(named_transient_orfs))
{
	sel<-which(star$Systematic.Name==orfs$is_gene[named_transient_orfs[i]])
	if(length(sel)>0)
	{
		starcount[i]<-sum(star[sel,3:24],na.rm=T)
	}
}

###detection using cswap databse from meurer

cswatm<-read.csv("/home/acwach/translatome2/cswap_meurer.csv")
cswatm0<-which(cswatm$Systematic.Name %in% orfs$is_gene[named_transient_orfs])# orfs$is_gene[transient_set_full])
cswat_found<-cswatm0[which(cswatm[cswatm0,]$fold_bkg.1>1.2|cswatm[cswatm0,]$fold_bkg.2>1.2|cswatm[cswatm0,]$fold_bkg>1.2)]
table(orfs$orf_class[which(orfs$is_gene %in% cswatm$Systematic.Name[cswatm0[cswat_found]])])
length(which(cswatm[cswatm0,]$fold_bkg.1>1.2|cswatm[cswatm0,]$fold_bkg.2>1.2|cswatm[cswatm0,]$fold_bkg>1.2))/length(cswatm0)
om<-which(orfs$is_gene%in%cswatm$Systematic.Name)
#tom<-intersect(om,transient_set_full)


# cswat<-read.csv("screen_C_SWAT_GFP-stats.csv",sep=" ")
# orfs$cswat<-cswat[match(orfs$is_gene,cswat$ORF),]$avgGFP
# orfs$cswatm1<-cswatm[match(orfs$is_gene,cswatm$Systematic.Name),]$intensity
# orfs$cswatm2<-cswatm[match(orfs$is_gene,cswatm$Systematic.Name),]$intensity.1
# orfs$cswatm3<-cswatm[match(orfs$is_gene,cswatm$Systematic.Name),]$intensity.2

# cswat_transient<-transient_set_full[(which(orfs$is_gene[transient_set_full] %in% cswat$ORF))]

#denominator:
possible_hits<- named_transient_orfs[which(orfs$is_gene[named_transient_orfs]%in%cyclops_all$ORF | orfs$is_gene[named_transient_orfs]%in% cswatm$Systematic.Name)]
observed_hits<-unique(named_transient_orfs[c(which(orfs$is_gene[named_transient_orfs]%in% cswatm[cswat_found,]$Systematic.Name),which(!is.na(cyclops_localization)))])

#discovered<-cswatm$Systematic.Name[which(cswatm$fold_bkg.1>1.2|cswatm$fold_bkg.2>1.2|cswatm$fold_bkg>1.2)]
#undiscovered<-cswatm$Systematic.Name[-which(cswatm$fold_bkg.1>1.2|cswatm$fold_bkg.2>1.2|cswatm$fold_bkg>1.2)]

#detected_protein<-named_transient_orfs[which(!is.na(cswat_localization)|starcount>0|!is.na(cyclops_localization))]

#named_transient_orfs[which(!is.na(cswat_localization)|starcount>0|!is.na(cyclops_localization))]

df_localize<-data.frame(categ=c("Total detected","Detected, unknown\nlocalization","Detected in\ncytoplasm","Detected in\nmitochondria","Detected in\nnucleus"),counts=c(length(detected_protein),length(detected_protein)-sum(table(localization)),table(localization)["Cytoplasm"],table(localization)["Mitochondria"],table(localization)["Nucleus"] ))

df_localize$categ<-factor(df_localize$categ,levels=rev(c("Total detected","Detected, unknown\nlocalization","Detected in\ncytoplasm","Detected in\nmitochondria","Detected in\nnucleus")))

require("ggplot2")
png("localize.png",height=3000,width=2000, res=300)
f6a<-ggplot(df_localize)+
	theme_classic(base_size=20) +
	#theme(axis.text.x=element_text(angle = 45, hjust=1))+
	geom_bar(aes(x=categ,y=counts),position="stack",stat="identity",width=.5)+
	labs(y="Annotated transient ORF\nproteins detected",x="")+
	coord_flip()
f6a
dev.off()

	
#annalysis of phenotypes in guo et al. 

#no transient ORF is essential
guo1<-read.csv("rep1_screens_guo2018.csv",stringsAsFactors=F)
guo1[which(guo1$Essentiality.screen.Target %in% orfs$is_gene[named_transient_orfs]),]$P.Values

#transient ORFs found in screens
allpheno<-vector()
sel<-which(guo1$Heat.screen.Target %in% orfs$is_gene[named_transient_orfs])
guo1$Heat.screen.Target[sel[which(guo1[sel,12]==T)]]
allpheno<-c(allpheno,guo1$Heat.screen.Target[sel[which(guo1[sel,12]==T)]])

sel<-which(guo1$HU.screen.Target %in% orfs$is_gene[named_transient_orfs])
guo1$HU.screen.Target[sel[which(guo1[sel,18]==T)]]
allpheno<-c(allpheno,guo1$HU.screen.Target[sel[which(guo1[sel,18]==T)]])

sel<-which(guo1$FLU.screen.Target %in% orfs$is_gene[named_transient_orfs])
guo1$FLU.screen.Target[sel[which(guo1[sel,24]==T)]]
allpheno<-c(allpheno,guo1$HU.screen.Target[sel[which(guo1[sel,18]==T)]],guo1$FLU.screen.Target[sel[which(guo1[sel,24]==T)]])

allguo<-unique(c(guo1$Essentiality.screen.Target,guo1$Heat.screen.Target,guo1$HU.screen.Target,guo1$FLU.screen.Target))
table(orfs$orf_class[which(orfs$is_gene %in% intersect(allguo,orfs$is_gene[named_transient_orfs]))])

####genetic interaction analysis 

####genetic interaction data

#interaction and data downloaded from thecellmap.org
NXN_int<-read.csv("SGA_NxN.txt",sep = "\t",stringsAsFactors=F)
EXE_int<-read.csv("SGA_ExE.txt",sep = "\t",stringsAsFactors=F)
NXE_int<-read.csv("SGA_ExN_NxE.txt",sep = "\t",stringsAsFactors=F)
strains<-read.csv("strain_ids_and_single_mutant_fitness.csv")

##add column with systematic name to all costanzo datasets
NXN_int$Query.sysname<-substr(NXN_int$Query.Strain.ID,1,regexpr('_',NXN_int$Query.Strain.ID)-1)
NXN_int$Array.sysname<-substr(NXN_int$Array.Strain.ID,1,regexpr('_',NXN_int$Array.Strain.ID)-1)
NXE_int$Query.sysname<-substr(NXE_int$Query.Strain.ID,1,regexpr('_',NXE_int$Query.Strain.ID)-1)
NXE_int$Array.sysname<-substr(NXE_int$Array.Strain.ID,1,regexpr('_',NXE_int$Array.Strain.ID)-1)
EXE_int$Query.sysname<-substr(EXE_int$Query.Strain.ID,1,regexpr('_',EXE_int$Query.Strain.ID)-1)
EXE_int$Array.sysname<-substr(EXE_int$Array.Strain.ID,1,regexpr('_',EXE_int$Array.Strain.ID)-1)

#find all transient orfss in costanzo by looking at all categories (503 matches)
all_names_in_cos<-unique(c(NXN_int$Query.sysname,NXN_int$Array.sysname,NXE_int$Query.sysname,NXE_int$Array.sysname,EXE_int$Query.sysname,EXE_int$Array.sysname))

allmatches<-list()
nematches<-list()
nnmatches<-list()
sys_nematches<-list()
sys_nnmatches<-list()

thres<-(-.2)
nnavail<-array()
neavail<-array()
stdevs<-list()
eps<-list()
sysmatches<-list()
sdthres<-10
anns<-orfs$is_gene[intersect(named_transient_orfs,intergene_index)]
for(i in 1:length(anns))
{
  nnavail[i]<-length(c(which(NXN_int$Query.sysname==anns[i]),which(NXN_int$Array.sysname==anns[i])))
  neavail[i]<-length(c(which(NXE_int$Query.sysname==anns[i]),which(NXE_int$Array.sysname==anns[i])))
  sel0<-which(NXN_int$Query.sysname==anns[i] & NXN_int$Genetic<(thres) & NXN_int$P.value<.05 & NXN_int$Double.mutant.fitness.standard.deviation<sdthres)
  sel1<-which(NXN_int$Array.sysname==anns[i] & NXN_int$Genetic<(thres) & NXN_int$P.value<.05 & NXN_int$Double.mutant.fitness.standard.deviation<sdthres)
  sel2<-which(NXE_int$Query.sysname==anns[i] & NXE_int$Genetic<(thres) & NXE_int$P.value<.05 & NXE_int$Double.mutant.fitness.standard.deviation<sdthres)
  sel3<-which(NXE_int$Array.sysname==anns[i] & NXE_int$Genetic<(thres) & NXE_int$P.value<.05 & NXE_int$Double.mutant.fitness.standard.deviation<sdthres)
  nnmatches[[i]]<-c(NXN_int[sel0,]$Array.allele.name,NXN_int[sel1,]$Query.allele.name)
  nematches[[i]]<-c(NXE_int[sel2,]$Array.allele.name,NXE_int[sel3,]$Query.allele.name)
  sys_nnmatches[[i]]<-c(NXN_int[sel0,]$Array.sysname,NXN_int[sel1,]$Query.sysname)
  sys_nematches[[i]]<-c(NXE_int[sel2,]$Array.sysname,NXE_int[sel3,]$Query.sysname)
 
  allmatches[[i]]<-c(nnmatches[[i]],nematches[[i]])
  sysmatches[[i]]<-c(sys_nnmatches[[i]],sys_nematches[[i]])  
  stdevs[[i]]<-c(NXN_int[sel0,]$Double.mutant.fitness.standard.deviation,NXN_int[sel1,]$Double.mutant.fitness.standard.deviation,NXE_int[sel2,]$Double.mutant.fitness.standard.deviation,NXE_int[sel3,]$Double.mutant.fitness.standard.deviation)
  eps[[i]]<-c(NXN_int[sel0,]$Genetic.interaction,NXN_int[sel1,]$Genetic.interaction,NXE_int[sel2,]$Genetic.interaction,NXE_int[sel3,]$Genetic.interaction)
}

count_matches<-array()
#write interaction files
for(i in 1:length(sysmatches))
{
  write.table(sysmatches[[i]],paste("inter",i,sep="",collapse=""),quote=F,col.names=F,row.names = F)
  count_matches[i]<-length(sysmatches[[i]])
}


which(NXN_int$Query.sysname%in%anns & NXN_int$P.value<.05 & NXN_int$Genetic.interaction.score..ε.<(-.2))
write.table(all_names_in_cos,"costanzo_genes",quote=F,col.names=F,row.names=F)

##Run Ontologizer on each list of genes
#java -jar Ontologizer.jar -g go.obo -a gene_association.sgd -p costanzo_genes -m Benjamini-Hochberg -c Term-For-Term -s inter



#read GO files
GO<-list()
for(i in 1:length(sysmatches))
{
  filename<-paste("table-inter",i,"-Term-For-Term-Benjamini-Hochberg.txt",sep="",collapse = "")
  if(file.exists(filename))
  {
    if(length(readLines(filename))>0)
    {
      GO[[i]]<-read.csv(filename,sep="\t")  
    }
  }
}

numgo<-rep(0,length(GO))
gohits<-array(0)
for(i in 1:length(GO))
{
  if(length(GO[[i]])>0)
  {
    numgo[i]<-length(which(GO[[i]]$p.adjusted<.05 & GO[[i]]$Study.term>1))
    if(numgo[i]>0)
    {
      gohits<-rbind(gohits,GO[[i]][1,])
    }
  }  
}

gohits<-gohits[-1,]
gohits<-cbind(anns[which(numgo>0)],gohits)

write.csv(gohits,"gohits0",quote=F,row.names=F)


#Genetic Interation Network Code ##Omer Acar
library(tidyverse)
library(glue)
library(igraph)
library(ggsignif)
library(cowplot)
library(readxl)
# Read Input files -----------------
# Read costanzo 2016 dataset downloaded on Mon, Nov  5 2018 from http://thecellmap.org/costanzo2016/
# SGA_NxN for nonessential-nonessential network
# SGA_ExN_NxE for essential-nonessential & nonessential-essential network
# SGA_ExE for essential-essential network

filename_SGA_NxN <- "SGA_NxN.txt"
filename_SGA_ExN_NxE <- "SGA_ExN_NxE.txt"
filename_SGA_ExE <- "SGA_ExE.txt"

SGA_NxN <- read_delim(filename_SGA_NxN, delim = "\t") %>% mutate(data_source = 'NxN')
SGA_ExN_NxE <- read_delim(filename_SGA_ExN_NxE, delim = "\t") %>% mutate(data_source = 'ExN_NxE')
SGA_ExE <- read_delim(filename_SGA_ExE, delim = "\t") %>% mutate(data_source = 'ExE')


# Mutate systematic names---------
# Original input files has extensions on SGD systematic names (like 'YAL002_sn273')
# This makes finding genes in the future analysis difficult, so I remove those parts after '_' character
# for both query and array id columns

SGA_NxN$`Array Strain ID` <- sapply(strsplit(SGA_NxN$`Array Strain ID`, "_"), "[", 1)
SGA_NxN$`Query Strain ID` <- sapply(strsplit(SGA_NxN$`Query Strain ID`, "_"), "[", 1)
SGA_ExN_NxE$`Array Strain ID` <- sapply(strsplit(SGA_ExN_NxE$`Array Strain ID`, "_"), "[", 1)
SGA_ExN_NxE$`Query Strain ID` <- sapply(strsplit(SGA_ExN_NxE$`Query Strain ID`, "_"), "[", 1)
SGA_ExE$`Array Strain ID` <- sapply(strsplit(SGA_ExE$`Array Strain ID`, "_"), "[", 1)
SGA_ExE$`Query Strain ID` <- sapply(strsplit(SGA_ExE$`Query Strain ID`, "_"), "[", 1)

# Get the data frame from edited input data
net.df <- bind_rows(SGA_NxN, SGA_ExN_NxE, SGA_ExE)
#following alleles has 2 temperature reported while supplementary methods of Costanzo 2016 says they used particular temperatures and reported only 1 according to quality. Thus I stick with their explanation and remove 30 degrees data
supplist<- c('vac8-supp1','mtg1-supp1','tfb1-6-supp1','mob2-11-supp1')
net.df.filtered=net.df %>% filter(!(`Query allele name`=='med6-ts'&`Arraytype/Temp`=='DMA30')&!(`Query allele name`%in%supplist&`Arraytype/Temp`=='DMA30')&!(`Query allele name`%in%supplist&`Arraytype/Temp`=='TSA30'))

# # this script reads SGA_data_combined data frame which contains costanzo 2016 data
# # and creates a data frame containing SGD systematic gene names and the Allele names used in the experiments

# Take query data and array data as separate data frames
# give meaningful column names to both
# combine two data frames, remove suppressor mutations, extract unique rows and save
q.data <- net.df[, c(1, 2)] %>% distinct()
a.data <- net.df[, c(3, 4)] %>% distinct()
colnames(q.data) <- c("Systematic gene name", "Allele Gene name")
colnames(a.data) <- c("Systematic gene name", "Allele Gene name")

strain_ids <- bind_rows(q.data, a.data) %>%
  distinct() %>%
  filter(grepl("supp", `Allele Gene name`) == F)


essetial.query <- net.df %>%
  filter(data_source=='ExE') %>%
  select(`Query allele name`) %>%
  pull()
essetial.array <- net.df %>%
  filter(data_source=='ExE') %>%
  select(`Array allele name`) %>%
  pull()
exn.query <- net.df %>%
  filter(data_source == "ExN_NxE" & (`Arraytype/Temp` == "TSA26" | `Arraytype/Temp` == "TSA30")==F) %>%
  select(`Query allele name`) %>%
  distinct() %>%
  pull()
essential.alleles <- unique(c(essetial.query,essetial.array,exn.query))

exp.number.data <- strain_ids %>%
  mutate(maincat = ifelse(`Allele Gene name`%in%essential.alleles,'essential','nonessential'))

write_csv(exp.number.data, "strain_ids_with_experiment_count_all.csv")

## Inputs -------------



#pgs <- readr::read_csv('transient_annotated_intergenic_041421')
overlappingorfs <- read_csv("overlappingorfs.csv")
#net.df <- readRDS("analysis/data/derived_data/SGA_data_combined.rds.gz")
net_df_significant <- filter(net.df, `P-value` <= 0.05 & `Query Strain ID` %in% overlappingorfs$orf_name == FALSE & `Array Strain ID` %in% overlappingorfs$orf_name == FALSE)
net_df_significant_sl <- filter(net_df_significant, `Genetic interaction score (ε)` <= -0.2)#& `Double mutant fitness standard deviation`<=0.1)
exp.data <- read_csv("strain_ids_with_experiment_count_all.csv")%>% mutate(group=ifelse(`Systematic gene name`%in%anns,'proto-gene',maincat))
strain_ids_smf <- read_excel("strain_ids_and_single_mutant_fitness.xlsx")

## Calculate strong (eps<-0.2) and lethal (eps<-0.35) interactions -----
strong_orf_list <- net.df %>% filter(`P-value`<=0.05,`Genetic interaction score (ε)` <= -0.2) %>% select(`Query Strain ID`, `Array Strain ID`) %>% as.list() %>% unlist() %>% unique()
lethal_orf_list <- net.df %>% filter(`P-value`<=0.05,`Genetic interaction score (ε)` <= -0.35) %>% select(`Query Strain ID`, `Array Strain ID`) %>% as.list() %>% unlist() %>% unique()
exp.data <- exp.data %>% mutate(strong_interaction=ifelse(`Systematic gene name`%in%strong_orf_list,TRUE,FALSE),
                          lethal_interaction=ifelse(`Systematic gene name`%in%lethal_orf_list,TRUE,FALSE))

exp_data_nones=exp.data %>% filter(group!='essential')

## Compare interaction ratios of transient orfs and nonessential orfs with fisher.test ----
strong_interaction_pval <- table(exp_data_nones$group,exp_data_nones$strong_interaction) %>% fisher.test() #0.04652
lethal_interaction_pval <- table(exp_data_nones$group,exp_data_nones$lethal_interaction) %>% fisher.test() #0.02636

percent_plot_data <- exp_data_nones %>%
  select(group,strong_interaction) %>%
  group_by(group,strong_interaction) %>%
  summarise(str_count=n()) %>%
  mutate(freq = str_count / sum(str_count)) %>%
  filter(strong_interaction==TRUE) %>% mutate(cat='strong') %>% #%>% gather()
  bind_rows(
    exp_data_nones %>% select(group,lethal_interaction) %>% group_by(group,lethal_interaction) %>% summarise(leth_count=n()) %>% mutate(freq = leth_count / sum(leth_count)) %>% filter(lethal_interaction==TRUE) %>% mutate(cat='lethal')
  ) %>% mutate(cat=factor(cat,levels=c('strong','lethal')))%>%
  ungroup(group) %>%
  mutate(group = factor(group, levels = c("proto-gene", "nonessential")))


percent_plot <-percent_plot_data %>%
  ggplot(aes(x=cat,y=freq,fill=group))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.95)) + theme_classic() + theme(plot.margin = margin(45,5.5,5.5,35, "pt"),
    axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),axis.title.y=element_text(size=20),
    legend.text = element_text(size = 15), legend.position = "right", legend.title = element_blank(), legend.margin = margin(t=-0.25,unit='in'),
    legend.spacing.x = unit(0.1, 'in')
  ) +
  coord_cartesian(clip="off")+
  #scale_fill_manual(labels = c("Transient ORFs", "Nonessential\nGenes"), values = c("#1CBDC2", "#EF3D23")) +
  scale_fill_manual(labels = c("Transient ORFs", "Nonessential\nGenes"), values = c("grey", "black")) +
  scale_x_discrete(labels = c(expression(epsilon * "< -0.2"), expression(epsilon * "<-0.35"))) +
  ylab("Percent with at least one\ninteraction at given threshold") + scale_y_continuous(labels = scales::percent) + xlab("") +
  geom_signif(
    y_position = c(1.05, 0.85), xmin = c(0.8, 1.8), xmax = c(1.2, 2.2),
    annotation = c(glue("p={format.pval(strong_interaction_pval$p.value,2)}"),glue("p={format.pval(lethal_interaction_pval$p.value,2)}")), tip_length = .1, textsize=6
  )
ggsave(plot = percent_plot,filename = 'Figure7B.pdf',width = 4,height = 4)


#compare smf to wild type fitness of 1 -----
strain_ids_smf$smf <- strain_ids_smf[, 4] %>%
  pull() %>%
  as.numeric()
strain_ids_joined <- exp.data %>%
  inner_join(strain_ids_smf) %>%
  pivot_wider(names_from = group, values_from = smf)

wt_pvalue <- strain_ids_joined$`proto-gene`[is.na(strain_ids_joined$`proto-gene`) == F] %>% t.test(mu = 1) #p=0.05663

#compare nonessential smf to transient smf distributions ---------
figure7a <- strain_ids_joined %>%
  # filter(subcat != "essential") %>%
  ggplot() +
  geom_histogram(aes(x = `proto-gene`), fill = "#21BDC2", color = "#21BDC2", binwidth=0.01) +
  geom_density(aes(x = nonessential), color = "#EF4024") +
  theme_classic() +
  scale_x_continuous(name = "Single Mutant Fitness", limits = c(0.5, NA)) +
  scale_y_continuous(name = "Counts", limits=c(0,20)) +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20), plot.margin = margin(5.5, 5.5, 5.5, 38.5,unit="pt"))
ggsave("Figure7A.pdf", plot=figure7a, width = 2, height = 2)

##Number of interactions histogram -----
colname <- colnames(exp.data)[2]
pgs_allele <- exp.data %>%
  filter(group=='proto-gene') %>%
  select(`Allele Gene name`) %>%
  pull()
allnet <- graph_from_data_frame(net_df_significant_sl[,c(2,4,6)],directed=FALSE)
allnet <- simplify(allnet)
allnet.deg<-degree(allnet)



#if(plot_histogram){
allnet.deg.pgs <- allnet.deg[names(allnet.deg)%in%pgs_allele]
f6d<-ggplot()+geom_histogram(binwidth=1,aes(x=allnet.deg.pgs),fill="#21BDC2")+theme_classic()+xlab('Number of interactions at  < -0.2')+
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20), 	plot.margin = margin(5.5,5.5,5.5,130, "pt"))+ylab('Count')
ggsave('Figure7D.pdf',width=2,height = 2)


## get yer175w-a ego graph to csv -----
induced_subgraph(allnet,ego(allnet,order=1,nodes='yer175w-a')[[1]])%>%get.edgelist()%>%as.data.frame()%>%
  write_csv('yer175w-a_ego.csv', col_names = FALSE)




## Interaction Density -----
thr <- -0.2
allnet <- graph_from_data_frame(net_df_significant[,c(2,4,6)],directed=FALSE)
adj <- as_adjacency_matrix(allnet,type='both',attr="Genetic interaction score (ε)") %>% as.matrix()
adj_nonw <- as_adjacency_matrix(allnet,type='both') %>% as.matrix()
colname <- colnames(exp.data)[2]

exp.data$essint <- NA
exp.data$nonesint <- NA
exp.data$esscount <- NA
exp.data$nonescount <- NA
for( i in seq_along(exp.data$`Allele Gene name`)){
  name <- exp.data$`Allele Gene name`[i]
  if(name%in%colnames(adj)){
    interactions <- adj[name,]
    int_nonw <- adj_nonw[name,]
    lethal_int_name <- names(interactions[interactions<=thr])
    exp.data$essint[i] <- sum(lethal_int_name%in%exp.data$`Allele Gene name`[exp.data$maincat=='essential'])
    exp.data$nonesint[i] <- sum(lethal_int_name%in%exp.data$`Allele Gene name`[exp.data$maincat=='nonessential'])
    #names(int_nonw[int_nonw==1])
    exp.data$esscount[i] <- sum(names(int_nonw[int_nonw==1])%in%exp.data$`Allele Gene name`[exp.data$maincat=='essential'])
    exp.data$nonescount[i] <- sum(names(int_nonw[int_nonw==1])%in%exp.data$`Allele Gene name`[exp.data$maincat=='nonessential'])
  }

}
library("coin")

exp.data <- exp.data %>%
  dplyr::mutate(essint_density = essint / esscount, nonesint_density = nonesint / nonescount)

pe12 = coin::independence_test(group~essint,exp.data %>% mutate(group = as.factor(group)) %>% filter(group!='essential'))
pn12 = coin::independence_test(group~nonesint,exp.data %>% mutate(group = as.factor(group)) %>% filter(group!='essential'))
#1.028e-05 nones~pgs for nones ints
#0.0165 nones~pgs for ess ints
#p-value = 0.00124 nones~pgs for nones ints
#p-value = 0.01571 nones~pgs for ess ints
exp.data$group <- factor(exp.data$group,levels = c('proto-gene','nonessential','essential'))
plt <- exp.data %>%
  select(group, essint_density, nonesint_density) %>%
  gather(type, int_density, -group) %>%
  ggplot(aes(x = type, y = int_density, fill = group)) + # geom_boxplot()
  stat_summary(fun.y = mean, geom = "bar", position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 1), width = 0.5) +
  ylab("Interaction density") +
  xlab("") +
  scale_x_discrete(labels = c("Interactions with\nEssential genes", "Interactions with\nNonssential genes")) +
  theme_classic() + theme(
	plot.margin = margin(5.5,5.5,5.5,35, "pt"),
    axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 15), legend.position = "bottom", legend.title = element_blank(), legend.margin = margin(t=-0.25,unit='in'),
    legend.spacing.x = unit(0.1, 'in')
  ) +
  scale_fill_manual(labels = c("Transient ORFs", "Nonessential\nGenes", "Essential\nGenes"), values = c("#1CBDC2", "#EF3D23", "#FAA51A")) +
  geom_signif(
    annotations = c(paste0("p= ", formatC(coin::pvalue(pe12), digits = 2)), paste0("p= ", formatC(coin::pvalue(pn12), digits = 1))),
    y_position = c(0.07, 0.05), xmin = c(0.65, 1.65), xmax = c(1, 2), tip_length = 0.005, textsize = 6
  )
ggsave('Figure7C.pdf',plot= plt , width = 4,height = 3)




#overall phenotype analysis

#strong genetic interaction 
#null mutant phenotype
#overexpression phenotype

phenotypes<-read.csv("SupplementaryTable5_Phenotypes.txt",sep="\t")
null_phenotype<-rep(0,length(named_transient_orfs))
oe_phenotype<-rep(0,length(named_transient_orfs))
gi_phenotype<-rep(0,length(named_transient_orfs))

null_phenotype[which(orfs$is_gene[named_transient_orfs] %in% phenotypes$Systematic.name[which(phenotypes$Evidence.type%in%c("Null mutant","Null mutant, overexpression mutant")  ) ])]<-1
oe_phenotype[which(orfs$is_gene[named_transient_orfs] %in% phenotypes$Systematic.name[which(phenotypes$Evidence.type%in%c("Overexpression mutant","Null mutant, overexpression mutant")  ) ])]<-1
gi_phenotype[which(orfs$is_gene[named_transient_orfs] %in% as.character(anns[which(count_matches>0)]))]<-1


#n1<-as.character(orfs[named_transient_orfs,]$is_gene)
#n2<-as.character(anns)
#cbind(n1[order(n1)],n2[order(n2)])

# library("eulerr")

# pheno_venn <-  c(A = sum(null_phenotype), B = sum(oe_phenotype), C = sum(gi_phenotype), 
                    # "A&B" = length(which(null_phenotype+oe_phenotype==2)),
					# "A&C" = length(which(null_phenotype+gi_phenotype==2)),
                    # "B&C" = length(which(gi_phenotype+oe_phenotype==2)),
					# "A&B&C" =length(which(gi_phenotype+oe_phenotype+null_phenotype==3))
					# )

# png("pheno_venn.png")					
# pheno_fit <- euler(pheno_venn, shape = "ellipse", input="union")
# f6f<-plot(pheno_fit,quantities = list(cex=1.5),adjust_labels=TRUE,legend = list(labels=c("Null mutant phenotype","Overexpression phenotype","Genetic interaction"),side = "right",cex=1.3))
# f6f
# dev.off()
# library("png")
# library("grid")

yer_png<-readPNG("G_cropped.png")
yer_fig <- rasterGrob(yer_png, interpolate=TRUE, width=.75)


se<-function(x) {return(sd(x)/sqrt(length(x)))}
se_pro<-function(p,n){return(sqrt(p*(1-p)/n))}

library("ggplot2")
df_micro<-data.frame(categ=c("Transient cORF","Transient nORF"),percent_detected=c(.81,.72))
df_micro$se<-c(se_pro(.81,88),se_pro(.72,36))
png("microscopy.png",height=500,width=500)
f6a<-ggplot(df_micro)+
	geom_bar(stat="identity",aes(x=categ,y=percent_detected,fill=categ))+
	geom_errorbar(aes(x=categ,ymin=percent_detected-se,ymax=percent_detected+se))+
	theme_classic(base_size=20)+
	theme(legend.position = "none",plot.margin = margin(5.5, 5.5, 5.5, 35.5,unit="pt"))+
	scale_fill_manual(values=c("red","blue"))+
	labs(x="",y="Proportion of tested proteins detected")
f6a
dev.off()

###phenotype heatmap

del_col<-read.table("strain_heterozygous_diploid.txt",sep="\t",header=1,stringsAsFactors=F)
del_col$ORF_name<-trimws(del_col$ORF_name)

unique(phenotypes$Systematic.name[which(phenotypes$Type.of.screen=="Yeast deletion collection")])

length(intersect(del_col$ORF_name,orfs$is_gene[named_transient_orfs]))
length(unique(phenotypes$Systematic.name[which(phenotypes$Type.of.screen=="Yeast deletion collection")]))


df_pheno<-data.frame(id=orfs$is_gene[named_transient_orfs])
df_pheno$orf_class<-as.character(orfs$orf_class[named_transient_orfs])
df_pheno$orf_class[which(df_pheno$orf_class!="Dubious")]<-(-2)#"Canonical"
df_pheno$orf_class[which(df_pheno$orf_class=="Dubious")]<-(-1)#"Noncanonical"
df_pheno$orf_class<-as.numeric(df_pheno$orf_class)

df_pheno$detected<-0
df_pheno[which(named_transient_orfs %in% possible_hits),]$detected<-1
df_pheno[which(named_transient_orfs %in% observed_hits),]$detected<-2

#df_pheno$localized<-0
#df_pheno[which(named_transient_orfs %in% possible_hits),]$localized<-1
#df_pheno[which(!is.na(localization)),]$localized<-2

df_pheno$del_col<-0
df_pheno[which(orfs$is_gene[named_transient_orfs] %in% del_col$ORF_name),]$del_col<-1
df_pheno[which(orfs$is_gene[named_transient_orfs] %in% phenotypes$Systematic.name[which(phenotypes$Type.of.screen=="Yeast deletion collection")]),]$del_col<-2

df_pheno$null_pheno<-1
df_pheno$null_pheno[which(orfs$is_gene[named_transient_orfs] %in% phenotypes$Systematic.name[which(phenotypes$Evidence.type=="Null mutant"|phenotypes$Evidence.type=="Null mutant, overexpression mutant")])]<-2

df_pheno$oe_pheno<-1
df_pheno$oe_pheno[which(orfs$is_gene[named_transient_orfs] %in% phenotypes$Systematic.name[which(phenotypes$Evidence.type=="Overexpression mutant"|phenotypes$Evidence.type=="Null mutant, overexpression mutant")])]<-2

df_pheno$gi_pheno<-0
df_pheno$gi_pheno[which(orfs$is_gene[named_transient_orfs] %in% intersect(all_names_in_cos, orfs$is_gene[intergene_index]))]<-1
df_pheno$gi_pheno[which(orfs$is_gene[named_transient_orfs] %in% as.character(anns[which(count_matches>0)]))]<-2

df_pheno$go_pheno<-0
df_pheno$go_pheno[which(orfs$is_gene[named_transient_orfs] %in% intersect(all_names_in_cos, orfs$is_gene[intergene_index]))]<-1
df_pheno$go_pheno[which(orfs$is_gene[named_transient_orfs] %in% anns[which(numgo>0)])]<-2

#df_pheno<-df_pheno[order(df_pheno$orf_class,df_pheno$detected,df_pheno$del_col,df_pheno$null_pheno,df_pheno$oe_pheno,df_pheno$gi_pheno,df_pheno$go_pheno),]
dfpx<-df_pheno[,3:8]
dfpx[dfpx==1]<-0
df_pheno<-df_pheno[rev(order(rowSums(dfpx))),]


library("ggplot2")
library("tidyr")
dfp<-df_pheno %>% pivot_longer(2:length(colnames(df_pheno)), names_to = "pheno", values_to = "status")

dfp$status<-factor(dfp$status,levels=c("-2","-1","4","0","1","2")) #as.character(dfp$status)
dfp$pheno<-factor(dfp$pheno,levels=rev(c("orf_class","detected","localized","del_col","null_pheno","oe_pheno","gi_pheno","go_pheno")))
dfp$id<-factor(dfp$id,levels=as.character(unique(dfp$id)))
dfp$height<-1
dfp$height[which(dfp$pheno=="orf_class")]<-.5

png("pheno_tile.png",width=2000,height=1000)
f6f<-ggplot(dfp)+
geom_tile(aes(y=pheno,x=id,fill=status,height=height))+
scale_fill_manual(labels = c("Canonical","Noncanonical","\n\n\n","Phenotype not tested", "Phenotype not found", "Phenotype shown"),values=c("red","blue","white","grey81","grey52","#21BDC2"),drop=FALSE)+
scale_y_discrete(labels=rev(c("","Protein\ndetected","Deletion\ncollection\nscreen","Null mutant\nscreen","Overexpresion\nmutant screen","Genetic\ninteraction","GO-associated\ninteractors")))+
xlab("Annotated transient ORFs")+
ylab("")+
theme_classic(base_size=20) +
theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.line=element_blank(),legend.title = element_blank(),plot.margin = margin(5.5, 5.5, 5.5, 35.5,unit="pt")) 
f6f
dev.off()



require("cowplot")

# abcd<-plot_grid(f6a,figure7a, percent_plot, plt,  ncol = 4, nrow=1, rel_widths=c(1,.8,1,1),labels = c('A', 'B', 'C', 'D'), label_size=30, align="vh",axis="b")
# efg<-plot_grid(f6d, yer_fig, f6f,  ncol = 3, nrow=1, rel_widths=c(1,.8,2),labels = c('E', 'F', 'G'), label_size=30, align="vh",axis="b")

# pdf("Figure6.pdf",width=1800*3/200,height=700*3/200)
# plot_grid(abcd,efg,  ncol = 1, nrow=2, align="vh",axis="b")
# dev.off()


# pdf("Figure6.pdf",width=1800*3/200,height=700*3/200)
# plot_grid(f6a,figure7a, percent_plot, plt, f6d, yer_fig, f6f,  ncol = 4, nrow=2, rel_widths=c(1,.8,1,1),labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G'), label_size=30, align="vh",axis="b")
# dev.off()


abc<-plot_grid(f6a,figure7a, percent_plot,ncol=3,nrow=1,rel_widths=c(.6,1,1),labels=c('A','B','C'),label_size=30)
def<-plot_grid(yer_fig,f6f,rel_widths=c(.5,1),labels=c('D','E'),label_size=30)

pdf("Figure6.pdf",width=1800*3/200,height=1000*3/200)
plot_grid(abc, def, ncol = 1, nrow=2)
dev.off()

png("Figure6.png",width=1800*3/200,height=1000*3/200,units='in',res=600)
plot_grid(abc, def, ncol = 1, nrow=2)
dev.off()


# pdf("Figure6.pdf",width=1800*3/200,height=700*3/200)
# plot_grid(f6a,figure7a, percent_plot, yer_fig, f6f,  ncol = 3, nrow=2, rel_widths=c(1,.8,1),labels = c('A', 'B', 'C', 'D', 'E'), label_size=30, align="vh",axis="b")
# dev.off()

# require("ggplot2")
# png("fig7a_pnps.png",height=480,width=800)
# f7a<-ggplot(df7, aes(x=category,y=omega)) +
  # geom_bar(stat="identity",aes(fill=category)) +
  # theme_classic(base_size=30) +
  # theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust=1)) +
  # geom_errorbar(aes(x=category,ymin=omega-se,ymax=omega+se)) +
  # labs(x="",y="pN/pS")+
  # #scale_fill_brewer(palette="pal8")+
  # scale_fill_manual(values=fig7_colors)+
  # geom_text(data=data.frame(), aes(x=levels(df7$cat), y=1.2, label=is_sigf7 ), col='black', size=10)+
  # geom_hline(yintercept=1, linetype="dashed", color="black", size=1.0)
# f7a
# dev.off()





#d

#named_transient_intergenic_orfs<-intersect(named_transient_orfs,named_transient_orfs_intergenic

####last figure

untrans_controls<-intersect(intersect(nontrans_index,intergene_index),noncanonical_index)
untrans_controls_antisense<-intersect(intersect(nontrans_index,antisense_index),noncanonical_index)

transient_dubious<-intersect(intergene_index,intersect(named_transient_orfs,which(orfs$orf_class=="Dubious")))
transient_can<-intersect(intergene_index,intersect(named_transient_orfs,which(orfs$orf_class!="Dubious")))

transient_unann<-intersect(transient_set,intergene_index) 
gene_control<-intersect(canonical_trans,conserved_index)

cat7<-c("Untranslated controls","Transient dubious nORFs", "Transient cORFs" , "Transient unannotated nORFs","Conserved cORFs")

omega7<-array()
omega7[1]<-get_omega(untrans_controls)
omega7[4]<-get_omega(transient_unann)
omega7[3]<-get_omega(transient_can)
omega7[2]<-get_omega(transient_dubious)
omega7[5]<-get_omega(gene_control)

omega7p<-array()
omega7p[1]<-get_omega_pval(untrans_controls)
omega7p[2]<-get_omega_pval(transient_dubious)
omega7p[3]<-get_omega_pval(transient_can)
omega7p[4]<-get_omega_pval(transient_unann)
omega7p[5]<-get_omega_pval(gene_control)

omega7se<-array()
omega7se[1]<-get_omega_se(untrans_controls)
omega7se[2]<-get_omega_se(transient_dubious)
omega7se[3]<-get_omega_se(transient_can)
omega7se[4]<-get_omega_se(transient_unann)
omega7se[5]<-get_omega_se(gene_control)

orfs_pi<-read.csv("nuc_diverse",sep=" ")
pi7<-array()
pi7[1]<-mean(orfs_pi$pi[untrans_controls])
pi7[2]<-mean(orfs_pi$pi[transient_dubious])
pi7[3]<-mean(orfs_pi$pi[transient_can])
pi7[4]<-mean(orfs_pi$pi[transient_unann])
pi7[5]<-mean(orfs_pi$pi[gene_control])


pi7se<-array()
pi7se[1]<-se(orfs_pi$pi[untrans_controls])
pi7se[2]<-se(orfs_pi$pi[transient_dubious])
pi7se[3]<-se(orfs_pi$pi[transient_can])
pi7se[4]<-se(orfs_pi$pi[transient_unann])
pi7se[5]<-se(orfs_pi$pi[gene_control])

spar<-array()
spar[1]<-mean((borfs$overlap1/orfs$length)[untrans_controls])
spar[2]<-mean((borfs$overlap1/orfs$length)[transient_dubious])
spar[3]<-mean((borfs$overlap1/orfs$length)[transient_can])
spar[4]<-mean((borfs$overlap1/orfs$length)[transient_unann])
spar[5]<-mean((borfs$overlap1/orfs$length)[gene_control])

sparse<-array()
sparse[1]<-se((borfs$overlap1/orfs$length)[untrans_controls])
sparse[2]<-se((borfs$overlap1/orfs$length)[transient_dubious])
sparse[3]<-se((borfs$overlap1/orfs$length)[transient_can])
sparse[4]<-se((borfs$overlap1/orfs$length)[transient_unann])
sparse[5]<-se((borfs$overlap1/orfs$length)[gene_control])


can_match<-array()
con_match<-array()
for(i in 1:length(transient_can))
{
	sel<-transient_unann[which(orfs$length[transient_unann]==orfs$length[transient_can[i]])]
	if(length(sel)==1)
	{
		can_match[i]<-sel
	}
	else if(length(sel)>1)
	{
		can_match[i]<-sample(sel,1)
	}
	sel<-gene_control[which(orfs$length[gene_control]==orfs$length[transient_can[i]])]
	a<-4
	while(length(sel)==0)
	{
		sel<-gene_control[which(abs(orfs$length[gene_control]-orfs$length[transient_can[i]])<a)]
		a<-a+3
	}
	if(length(sel)==1)
	{
		con_match[i]<-sel
	}
	else if(length(sel)>1)
	{
		con_match[i]<-sample(sel,1)
	}
}


mean(log((ribo$frame0/orfs$length)[transient_can]))
mean(log((ribo$frame0/orfs$length)[can_match]))
mean(log((ribo$frame0/orfs$length)[con_match]))


trans7<-array()
#trans7[1]<-mean((ribo$frame0/orfs$length)[untrans_controls])
trans7[2]<-mean(log((ribo$frame0/orfs$length)[transient_dubious]))
trans7[3]<-mean(log((ribo$frame0/orfs$length)[transient_can]))
trans7[4]<-mean(log((ribo$frame0/orfs$length)[can_match]))
trans7[5]<-mean(log((ribo$frame0/orfs$length)[con_match]))

trans7se<-array()
trans7se[2]<-se(log((ribo$frame0/orfs$length)[transient_dubious]))
trans7se[3]<-se(log((ribo$frame0/orfs$length)[transient_can]))
trans7se[4]<-se(log((ribo$frame0/orfs$length)[transient_can]))
trans7se[5]<-se(log((ribo$frame0/orfs$length)[con_match]))


cs7<-array()
cs7[1]<-mean(cs$score[untrans_controls])
cs7[2]<-mean(cs$score[transient_dubious])
cs7[3]<-mean(cs$score[transient_can])
cs7[4]<-mean(cs$score[transient_unann])
cs7[5]<-mean(cs$score[gene_control])

cs7se<-array()
cs7se[1]<-se(cs$score[untrans_controls])
cs7se[2]<-se(cs$score[transient_dubious])
cs7se[3]<-se(cs$score[transient_can])
cs7se[4]<-se(cs$score[transient_unann])
cs7se[5]<-se(cs$score[gene_control])


len7<-array()
len7[1]<-mean(orfs$length[untrans_controls])
len7[2]<-mean(orfs$length[transient_dubious])
len7[3]<-mean(orfs$length[transient_can])
len7[4]<-mean(orfs$length[transient_unann])
len7[5]<-mean(orfs$length[gene_control])

len7se<-array()
len7se[1]<-se(orfs$length[untrans_controls])
len7se[2]<-se(orfs$length[transient_dubious])
len7se[3]<-se(orfs$length[transient_can])
len7se[4]<-se(orfs$length[transient_unann])
len7se[5]<-se(orfs$length[gene_control])

df7<-data.frame(category=cat7, omega=omega7, se=omega7se, pi7=pi7, pi7se=pi7se, spar=spar,sparse=sparse, trans=trans7, transse=trans7se, cs=cs7, csse=cs7se, len=len7, lense=len7se)
df7$category<-factor(df7$category,levels=cat7)


fig7_colors<-c("grey","red","blue","red","blue")
is_sigf7<-c("","","","","***")
require("ggplot2")
png("fig7a_pnps.png",height=480,width=800)
f7a<-ggplot(df7, aes(x=category,y=omega)) +
  geom_bar(stat="identity",aes(fill=category)) +
  theme_classic(base_size=30) +
  theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust=1)) +
  geom_errorbar(aes(x=category,ymin=omega-se,ymax=omega+se)) +
  labs(x="",y="pN/pS")+
  #scale_fill_brewer(palette="pal8")+
  scale_fill_manual(values=fig7_colors)+
  geom_text(data=data.frame(), aes(x=levels(df7$cat), y=1.2, label=is_sigf7 ), col='black', size=10)+
  geom_hline(yintercept=1, linetype="dashed", color="black", size=1.0)
f7a
dev.off()

png("fig7b_pi.png",height=480,width=800)
f7b<-ggplot(df7, aes(x=category,y=pi7)) +
  geom_bar(stat="identity",aes(fill=category)) +
  theme_classic(base_size=30) +
  theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust=1)) +
  geom_errorbar(aes(x=category,ymin=pi7-pi7se,ymax=pi7+pi7se)) +
  labs(x="",y="Nucleotide diversity")+
  #scale_fill_brewer(palette="pal8")+
  scale_fill_manual(values=fig7_colors)+
  geom_text(data=data.frame(), aes(x=levels(df7$cat), y=.009, label=is_sigf7 ), col='black', size=10)
f7b
dev.off()

png("fig7c_spar.png",height=480,width=800)
f7c<-ggplot(df7, aes(x=category,y=spar)) +
  geom_bar(stat="identity",aes(fill=category)) +
  theme_classic(base_size=30) +
  theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust=1)) +
  geom_errorbar(aes(x=category,ymin=spar-sparse,ymax=spar+sparse)) +
  labs(x="",y="RFC with S. paradoxus")+
  #scale_fill_brewer(palette="pal8")+
  scale_fill_manual(values=fig7_colors)+
  geom_text(data=data.frame(), aes(x=levels(df7$cat), y=1, label=is_sigf7 ), col='black', size=10)
f7c
dev.off()

png("fig7d_trans.png",height=480,width=800)
f7d<-ggplot(df7[2:5,], aes(x=category,y=trans)) +
  geom_bar(stat="identity",aes(fill=category)) +
  #scale_y_log10()+
  theme_classic(base_size=30) +
  theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust=1)) +
  geom_errorbar(aes(x=category,ymin=trans-transse,ymax=trans+transse)) +
  labs(x="",y="Log ribo-seq reads per base")+
  #scale_fill_brewer(palette="pal8")+
  scale_fill_manual(values=fig7_colors[2:5])
f7d
dev.off()

png("fig7f_cs.png",height=480,width=800)
f7f<-ggplot(df7, aes(x=category,y=cs)) +
  geom_bar(stat="identity",aes(fill=category)) +
  #scale_y_log10()+
  theme_classic(base_size=30) +
  theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust=1)) +
  geom_errorbar(aes(x=category,ymin=cs-csse,ymax=cs+csse)) +
  labs(x="",y="Coding score")+
  #scale_fill_brewer(palette="pal8")+
  scale_fill_manual(values=fig7_colors)
f7f
dev.off()

png("fig7g_len.png",height=480,width=800)
f7g<-ggplot(df7[2:4,], aes(x=category,y=len)) +
  geom_bar(stat="identity",aes(fill=category)) +
  #scale_y_log10()+
  theme_classic(base_size=30) +
  theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust=1)) +
  geom_errorbar(aes(x=category,ymin=len-lense,ymax=len+lense)) +
  labs(x="",y="Length (bp)")+
  #scale_fill_brewer(palette="pal8")+
  scale_fill_manual(values=fig7_colors[2:4])
f7g
dev.off()

require("cowplot")
png("Figure7.png",width=2000,height=1800)
plot_grid(f7a,f7b, f7c, f7d, f7f, f7g, ncol = 3, nrow=2, rel_widths=c(1,1,1),rel_heights=c(1,1),labels = c('A', 'B', 'C', 'D', 'E', 'F'), label_size=30,align="hv")
dev.off()