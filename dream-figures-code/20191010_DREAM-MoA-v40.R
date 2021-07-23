



####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########      ########      ########     ########             #########           ##       ############      ##
##     ##     ##     ##     ##     ##           ##     ##            ##      ##         ####           ##          ####
##     ##     ##     ##    ##      ##           ##     ##            ##       ##       ##  ##          ##         ##  ##
##     ########      #######       ########     ########             ##        ##     ##    ##         ##        ##    ##
##     ##            ##    ##      ##           ##                   ##       ##     ##########        ##       ##########
##     ##            ##     ##     ##           ##                   ##      ##     ##        ##       ##      ##        ##
##     ##            ##      ##    ########     ##                   #########     ##          ##      ##     ##          ##
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#1 SLICE RNAseq & DOSE-RESPONSE BY HAVING KINOME DATA:
#2 SLICE THAT DATA BY HAVING DRUGS IN EVERY LINE (except TCCSUP & BT20) and DOSE-RESPONSE
#3 FINAL PRODUCT:  71 kinase inhibitors with dose-response & RNAseq in 11 lines

setwd("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA")


############################################################################################################################
#1 SLICE KINOME_KD section of PANACEA (dose-response, RNAseq, annotation)
############################################################################################################################

{

###############################
#CUMC RAW DATA:
###############################
#LOAD entrez2hugo converter
#load("/Volumes/ifs-1/archive/shares/n-of-1-plateseq-data/n1screen_paper/n1database_org/context_annotations/entrez2hugo_LIST.RData")
#LOAD truSeq:
#load("/Volumes/ifs-1/archive/shares/n-of-1-plateseq-data/truseq-clines_dset.rda")

#LOAD DOSE-RESPONSES:
load("/Volumes/ifs-1/archive/shares/n-of-1-plateseq-data/n1screen_paper/n1database_org/HTS_doseresponse_final.RData")

#LOAD RAW RNAseq:   MOCK (steady-state), DMSO, DMSO+drug (perturbation)
load("/Volumes/ifs-1/archive/shares/n-of-1-plateseq-data/n1screen_paper/n1database_org/data/master_n1screen_expCOMBAT.RData")

#load("/Volumes/ifs-1/archive/shares/n-of-1-plateseq-data/n1screen_paper/n1database_org/data/individual lines/indiv_n1screen_expCOMBAT.RData")


###############################
#GOLD STD:  kinome bining vs DrugBank
###############################

#LOAD Kinome-binding
load("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-05-24 - N1screen - drug-clean-ness/OVERLAP_n1-Kinome-Screens.RData")
rm(kinome_selectivity_scores_OVERLAP)

#write.csv(n1druglibrary_annotation_CTEP_OVERLAP,file="n1druglibrary_annotation_CTEP_OVERLAP.csv")

#SIMPLIFIED FURTHER BY TOP 84 drugs:
#load("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-09-20-PANGEA-paper/simple_PANGEA-HALLMARKS_KinomeBinding_overlap.RData")



############################################################################################################################
#SLICE OUT KINASE-INHIBITORS RNAseq (which have Kd's)
############################################################################################################################
#EXTRACT KINASE-INHIBITORS HAVE Kds for:
KINASE_DRUGS_CAPS<-c(rownames(kinome_inhibitor_matrix_pKa_OVERLAP),"DMSO","MOCK","UNTREATED")

#EXTRACT PANACEA DRUG-NAMES IN ORDER
panacea_master_drug_list<-sapply(strsplit(colnames(master_n1screen_expCOMBAT),"_"),function(x) x[1])
panacea_master_drug_list_CAPS<-gsub("-","",gsub("_","",gsub(" ","",toupper(panacea_master_drug_list))))


#SLICE PANACEA Based on Inhibitors: 3,211 / 10,172
kinome_inhibitors_idx<-which(panacea_master_drug_list_CAPS %in% KINASE_DRUGS_CAPS)

master_n1screen_expCOMBAT_DREAMslice<-master_n1screen_expCOMBAT[,kinome_inhibitors_idx]



save(master_n1screen_expCOMBAT_DREAMslice,file="master_n1screen_expCOMBAT_DREAMslice.RData")


############################################################################################################################
#SLICE OUT KINASE-INHIBITORS DOSE RESPONSES
############################################################################################################################
unique_drugs<-unique(sapply(strsplit(colnames(master_n1screen_expCOMBAT_DREAMslice),"_"),function(x) x[1]))

HTS_dose_response_DREAMslice<-list()

#line="ASPC1"
for (line in names(hts_dose_response_data_fixed)){
  #ID columns with kinase inhibitors w/n DREAM DRUGS (unique_drugs above)
  avg_matrix<-hts_dose_response_data_fixed[[line]]$'avg'
  drug_names_order<-colnames(avg_matrix)
  overlap_idx<-which(drug_names_order %in% unique_drugs)
  
  #SLICE/STORE dose-response slice:
  HTS_dose_response_DREAMslice[[line]]<-avg_matrix[,overlap_idx]
    
  
}



save(HTS_dose_response_DREAMslice,file="HTS_dose_response_DREAMslice.RData")

}

############################################################################################################################
#2  FURTHER SLICE #1 by having dose-response, RNAseq & annotation across ALL LINES (except tccsup and bt20 as have slightly different libraryies...)
############################################################################################################################

{
load("HTS_dose_response_DREAMslice.RData")

load("master_n1screen_expCOMBAT_DREAMslice.RData")



############################################################################################################################
#RNAseq Census
############################################################################################################################

#DRUG CENSUS:  Remove everything <9 lines (Erlotinib)
sort(table(sapply(strsplit(colnames(master_n1screen_expCOMBAT_DREAMslice),"_"),function(x) x[1])),decreasing=T)/2
#DMSO    UNTREATED     Afatinib     Axitinib      AZD2014    Bosutinib Cabozantinib    Ceritinib   Crizotinib   Dabrafenib    Dasatinib    Gefitinib    Ibrutinib     Imatinib 
#282.0        282.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0 
#Lapatinib   Lenvatinib    Nilotinib  Osimertinib    Pazopanib    Ponatinib  Regorafenib  Ruxolitinib  Saracatinib  Selumetinib    Sorafenib    Sunitinib  Tofacitinib   Trametinib 
#13.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0         13.0 
#Vandetanib  Vemurafenib       AEE788    Alisertib       AMG900       AT9283      AZD1480      AZD5363    Bafetinib       BI2536     Brivanib    Cediranib     CP724714   Crenolanib 
#13.0         13.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0 
#Dacomitinib   Dinaciclib    Dovitinib  Enzastaurin    Foretinib    GSK461364     Icotinib       KW2449    Linifanib   Linsitinib      MGCD265       MK1775       MK2206    Motesanib 
#12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0 
#Neratinib   Nintedanib   PF04691502   Pimasertib  Quizartinib       RAF265   Rigosertib       TAK733    Telatinib   Tivantinib    Tivozanib   Varlitinib Galunisertib      MLN2480 
#12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         12.0         11.5         11.0 
#Momelotinib    Vatalanib   Volasertib       AMG208  Baricitinib  Binimetinib   Canertinib   Danusertib     ENMD2076 Gilteritinib       P27600   Pacritinib  Refametinib   Tandutinib 
#11.0         11.0         11.0         10.0         10.0         10.0         10.0         10.0         10.0         10.0         10.0         10.0         10.0         10.0 
#UCN01    Erlotinib Sotrastaurin   Fedratinib      AZD7762  Cobimetinib  Palbociclib   GSK1059615    GSK690693     IMD 0354      OSI-930 
#10.0          9.0          8.0          7.0          6.0          6.0          3.0          2.0          2.0          2.0          1.0 



#CELL-LINE CENSUS:  remove TCCSUP (too small a library) remove BT20 (to different a library)
sort(table(sapply(strsplit(colnames(master_n1screen_expCOMBAT_DREAMslice),"_"),function(x) x[length(x)])),decreasing=T)
#  aspc  hf2597   lncap hcc1143   du145   efo21   panc1     u87   h1793    hsts    krj1    bt20     tcc 
#296     296     294     292     272     270     266     234     232     218     218     197     126 


############################################################################################################################
#RNAseq Census
############################################################################################################################


sort(sapply(HTS_dose_response_DREAMslice, function(x) ncol(x)),decreasing=TRUE)
#  H1793     U87   DU145   EFO21 HCC1143   ASPC1  HF2597    BT20   LNCAP    HSTS    KRJ1   PANC1  TCCSUP 
#92      91      89      89      89      88      88      87      87      86      86      86      86 





############################################################################################################################
#FINAL LIST OF DRUGS:   In every cell-line EXCEPT TCCSUP and BT20 (too different drug-library)
############################################################################################################################

RNAseq_minusTCC<-master_n1screen_expCOMBAT_DREAMslice[,-which(sapply(strsplit(colnames(master_n1screen_expCOMBAT_DREAMslice),"_"),function(x) x[length(x)]) %in% c("tcc","bt20"))]

sort(table(sapply(strsplit(colnames(RNAseq_minusTCC),"_"),function(x) x[1])),decreasing=T)/2

#DMSO    UNTREATED       AEE788     Afatinib    Alisertib       AMG900       AT9283     Axitinib      AZD1480      AZD2014      AZD5363    Bafetinib       BI2536    Bosutinib 
#256          256           11           11           11           11           11           11           11           11           11           11           11           11 
#Brivanib Cabozantinib    Cediranib    Ceritinib     CP724714   Crenolanib   Crizotinib   Dabrafenib  Dacomitinib    Dasatinib   Dinaciclib    Dovitinib  Enzastaurin    Foretinib 
#11           11           11           11           11           11           11           11           11           11           11           11           11           11 
#Galunisertib    Gefitinib    GSK461364    Ibrutinib     Icotinib     Imatinib       KW2449    Lapatinib   Lenvatinib    Linifanib   Linsitinib      MGCD265       MK1775       MK2206 
#11           11           11           11           11           11           11           11           11           11           11           11           11           11 
#MLN2480  Momelotinib    Motesanib    Neratinib    Nilotinib   Nintedanib  Osimertinib    Pazopanib   PF04691502   Pimasertib    Ponatinib  Quizartinib       RAF265  Regorafenib 
#11           11           11           11           11           11           11           11           11           11           11           11           11           11 
#Rigosertib  Ruxolitinib  Saracatinib  Selumetinib    Sorafenib    Sunitinib       TAK733    Telatinib   Tivantinib    Tivozanib  Tofacitinib   Trametinib   Vandetanib   Varlitinib 
#11           11           11           11           11           11           11           11           11           11           11           11           11           11 
#Vatalanib  Vemurafenib   Volasertib       AMG208  Baricitinib  Binimetinib   Canertinib   Danusertib     ENMD2076 Gilteritinib       P27600   Pacritinib  Refametinib   Tandutinib 
#11           11           11            9            9            9            9            9            9            9            9            9            9            9 
#UCN01    Erlotinib Sotrastaurin   Fedratinib      AZD7762  Cobimetinib   GSK1059615    GSK690693     IMD 0354  Palbociclib      OSI-930 
#9            8            8            7            6            5            2            2            2            2            1 



drugs_in_everyline<-names(which(table(sapply(strsplit(colnames(RNAseq_minusTCC),"_"),function(x) x[1]))==22))

#MAKE VECTOR TO CHECK AGAINST EVERY DOSE-RESPONSE DATA:  good all 71 drugs have dose-response data for all
drug_rna_dose_intersect<-drugs_in_everyline

line_names<-names(HTS_dose_response_DREAMslice)[-which(names(HTS_dose_response_DREAMslice) %in% c("TCCSUP","BT20"))]
for (line in line_names){
  drug_rna_dose_intersect<-intersect(colnames(HTS_dose_response_DREAMslice[[line]]),drug_rna_dose_intersect)
}

annotations_71drugs<-which(rownames(n1druglibrary_annotation_CTEP) %in% drugs_in_everyline)


###################################################################################
#MAKE FINAL DATA-SETS:   annotation, RNAseq + dose-responses
###################################################################################

#SLICE ANNOTATIONS:
DREAM_drugLibrary_annotation<-n1druglibrary_annotation_CTEP[drugs_in_everyline,]

#SLICE DOSE-RESPONSES:
DREAM_Dose_Responses<-list()
line_names<-c("ASPC1","DU145","EFO21","H1793","HCC1143","HF2597","HSTS","KRJ1","LNCAP","PANC1","U87")
#line="ASPC1"
for (line in line_names){

  #SLICE/STORE dose-response slice:
  DREAM_Dose_Responses[[line]]<-HTS_dose_response_DREAMslice[[line]][,drugs_in_everyline]
  
  
}


#SLICE RNAseq:
drugset_idx<-which(sapply(strsplit(colnames(RNAseq_minusTCC),"_"),function(x) x[1]) %in% c(drugs_in_everyline,"MOCK","DMSO","UNTREATED"))

DREAM_RNAseq<-RNAseq_minusTCC[,drugset_idx]


#SANITY CHECK on RNAseq:
sort(table(sapply(strsplit(colnames(DREAM_RNAseq),"_"),function(x) x[1])),decreasing=T)/2
#DMSO    UNTREATED       AEE788     Afatinib    Alisertib       AMG900       AT9283     Axitinib      AZD1480      AZD2014      AZD5363    Bafetinib       BI2536    Bosutinib     Brivanib Cabozantinib    Cediranib 
#256          256           11           11           11           11           11           11           11           11           11           11           11           11           11           11           11 
#Ceritinib     CP724714   Crenolanib   Crizotinib   Dabrafenib  Dacomitinib    Dasatinib   Dinaciclib    Dovitinib  Enzastaurin    Foretinib Galunisertib    Gefitinib    GSK461364    Ibrutinib     Icotinib     Imatinib 
#11           11           11           11           11           11           11           11           11           11           11           11           11           11           11           11           11 
#KW2449    Lapatinib   Lenvatinib    Linifanib   Linsitinib      MGCD265       MK1775       MK2206      MLN2480  Momelotinib    Motesanib    Neratinib    Nilotinib   Nintedanib  Osimertinib    Pazopanib   PF04691502 
#11           11           11           11           11           11           11           11           11           11           11           11           11           11           11           11           11 
#Pimasertib    Ponatinib  Quizartinib       RAF265  Regorafenib   Rigosertib  Ruxolitinib  Saracatinib  Selumetinib    Sorafenib    Sunitinib       TAK733    Telatinib   Tivantinib    Tivozanib  Tofacitinib   Trametinib 
#11           11           11           11           11           11           11           11           11           11           11           11           11           11           11           11           11 
#Vandetanib   Varlitinib    Vatalanib  Vemurafenib   Volasertib 
#11           11           11           11           11 




save(list=c("DREAM_drugLibrary_annotation","DREAM_Dose_Responses","DREAM_RNAseq"),file="DREAM_kinaseInhibitor_ALL-DATA.RData")


}


############################################################################################################################
#3  PLOT DRUG-OFF-TARGET-NESS
############################################################################################################################




load("/Volumes/ifs-1/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-05-24 - N1screen - drug-clean-ness/drugClassAvg_offtargets_analysis.RData")

library("gplots")
library("RColorBrewer")


hcluster<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
  #distance matrix:
  dist_columns <- dist(t(matrix))
  dist_rows <- dist(matrix)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "NES",key.xlab = "NES")
}





###################################################
#CLUSTERING RESULTS:  Assess GLOBAL-KINOME EFFECTS
###################################################
top_ranked_kinase_var2<-names(sort(apply(n1drugClass_vs_kinome,2,var),decreasing=TRUE))
hcluster(n1drugClass_vs_kinome[,top_ranked_kinase_var2[1:100]],col_size=0.5,row_size = 1,title="Canonical-Drug-Class-Avg Off-Targets")

###################################################
#UNCLUSTERED matching kinase order: ASSESS ON-TARGET
###################################################
overlapping_kinases<-intersect(rownames(n1drugClass_vs_kinome),colnames(n1drugClass_vs_kinome))

matrix<-n1drugClass_vs_kinome[overlapping_kinases,overlapping_kinases]
dist_columns <- dist(t(matrix))
dist_rows <- dist(matrix)

#hclustering:
hclust_columns<-hclust(dist_columns, method="average")
hclust_rows<-hclust(dist_rows, method="average")


heatmap.2(matrix,margin=c(7,10), 
          Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_rows), cexCol =1.5,cexRow = 1.5,
          col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",xlab="kinases-binding",ylab="kinase-class",
          main="Kinase-Inhibitors On-Target Matrix",keysize=1,key.title = "-log10(Kd)",key.xlab = "-log10(Kd)")






####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########      ########      ########     ########    #########    ####     ##  ##########                  #########        ##      ##########       ##                                    
##     ##     ##     ##     ##     ##          ##           ##           ## ##    ##      ##                      ##      ##      ####         ##          ####                              
##     ##     ##     ##    ##      ##          ##           ##           ##  ##   ##      ##                      ##       ##    ##  ##        ##         ##  ##                             
##     ########      #######       ########     #######     ########     ##   ##  ##      ##                      ##       ##   ##    ##       ##        ##    ##                                
##     ##            ##    ##      ##                 ##    ##           ##    ## ##      ##                      ##       ##  ##########      ##       ##########                                                                                             
##     ##            ##     ##     ##                 ##    ##           ##     ####      ##                      ##      ##   ##      ##      ##       ##      ##                           
##     ##            ##      ##    ########     #######     #########    ##      ###      ##                      #########    ##      ##      ##       ##      ##                                      
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#1   PRESENT DIFFERENT 30 SUBSETS / 70 DRUG DATA to present to DREAM PEOPLE NEXT WEEK
#2   FINAL PREP OF DATA TO GIVE THEM



#setwd("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA")
setwd("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/")

load("DREAM_kinaseInhibitor_ALL-DATA.RData")


#write.csv(DREAM_drugLibrary_annotation,file="DREAM_drugLibrary_annotation.csv")


#LOAD Kinome-binding
load("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-05-24 - N1screen - drug-clean-ness/OVERLAP_n1-Kinome-Screens.RData")
rm(kinome_selectivity_scores_OVERLAP)

##################################################################################################################################
#PLOTTING FUNCTIONS:
##################################################################################################################################

library("gplots", lib.loc="~/Library/R/3.3/library")
library("RColorBrewer", lib.loc="~/Library/R/3.3/library")
library("devtools")

heatmap.3 <- function(x,Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,distfun = dist,hclustfun = hclust,dendrogram = c("both","row", "column", "none"),symm = FALSE,scale = c("none","row", "column"),na.rm = TRUE,revC = identical(Colv,"Rowv"),add.expr,
                      breaks,symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",col = "heat.colors",colsep,rowsep,sepcolor = "white",sepwidth = c(0.05, 0.05),cellnote,notecex = 1,notecol = "cyan",na.color = par("bg"),
                      trace = c("none", "column","row", "both"),tracecol = "cyan",hline = median(breaks),vline = median(breaks),linecol = tracecol,margins = c(5,5),ColSideColors,RowSideColors,side.height.fraction=0.3,cexRow = 0.2 + 1/log10(nr),cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,labCol = NULL,key = TRUE,keysize = 1.5,density.info = c("none", "histogram", "density"),denscol = tracecol,symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,main = NULL,xlab = NULL,ylab = NULL,lmat = NULL,lhei = NULL,lwid = NULL,ColSideColorsSize = 1,RowSideColorsSize = 1,KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}



hcluster_corr_n1screen_DrugTargets_hallmarks<-function(matrix,database,target_vector,target_colors,main_title){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix_idx<-which(rownames(matrix) %in% rownames(database))
  
  matrix<-matrix[matrix_idx,]
  
  
  
  #FILTER MATRIX BY ONLY DRUGS ANNOTATED: (reduce noise of drugs that do nothing...)
  master_indx_filter<-c()
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #2 extract class-indices in rlab
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(rownames(matrix),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 make master index to filter matrix
    master_indx_filter<-c(master_indx_filter,rlab_class_idx)
  }
  master_indx_filter_unique<-unique(master_indx_filter)
  
  matrix<-matrix[master_indx_filter_unique,]
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(target_vector),ncol=length(rownames(matrix)))
  rownames(rlab)<-target_vector
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #1 extract drugs-vector from n1screen_drugDB
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(colnames(rlab),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-target_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  top_var<-names(sort(apply(matrix,2,var),decreasing=T))[1:20]
  matrix_topvar<-matrix[,top_var]
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  cor_columns <- cor(matrix_topvar,method="spearman")
  cor_rows <- cor(t(matrix_topvar),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(length(target_vector)>1){
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(20,12),
              RowSideColors=rlab, symbreaks=TRUE, key=TRUE, symkey=TRUE,
              density.info="none", trace="none", main=main_title, cexRow=1.5,cexCol=1.2,col=colorRampPalette(c("blue","blue", "white", "red", "red"))(n = 20),
              RowSideColorsSize=5, KeyValueName="NES",keysize = 0.8,xlab="Hallmarks",ylab="Drug-Averages (across cell-lines)")
  }else{
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(20,12),
              symbreaks=TRUE, key=TRUE, symkey=TRUE,
              density.info="none", trace="none", main=main_title, cexRow=1.5,cexCol=1.2,col=colorRampPalette(c("blue","blue", "white", "red", "red"))(n = 20),
              KeyValueName="NES",keysize = 0.8,xlab="Hallmarks",ylab="Drug-Averages (across cell-lines)")
  }
  
  
}



hcluster_euclid_n1screen_DrugTargets_kinasesDREAM<-function(matrix,database,target_vector,target_colors,main_title,num_var=30){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix_idx<-which(rownames(matrix) %in% rownames(database))
  
  matrix<-matrix[matrix_idx,]
  
  
  
  #FILTER MATRIX BY ONLY DRUGS ANNOTATED: (reduce noise of drugs that do nothing...)
  master_indx_filter<-c()
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #2 extract class-indices in rlab
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(rownames(matrix),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 make master index to filter matrix
    master_indx_filter<-c(master_indx_filter,rlab_class_idx)
  }
  master_indx_filter_unique<-unique(master_indx_filter)
  
  matrix<-matrix[master_indx_filter_unique,]
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(target_vector),ncol=length(rownames(matrix)))
  rownames(rlab)<-target_vector
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #1 extract drugs-vector from n1screen_drugDB
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(colnames(rlab),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-target_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  top_var<-names(sort(apply(matrix,2,var),decreasing=T))[1:num_var]
  matrix_topvar<-matrix[,top_var]
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  
  dist_columns <- dist(t(matrix_topvar))
  dist_rows <- dist(matrix_topvar)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(length(target_vector)>1){
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(12,12),
              RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
              density.info="none", trace="none", main=main_title, cexRow=1.5,cexCol=1.2,col=colorRampPalette(c("blue","white", "red"))(n = 20),
              RowSideColorsSize=5, KeyValueName="-log10(Kd)",keysize = 0.8,xlab="Kinases",ylab=paste0(as.character(nrow(matrix_topvar))," Drugs"))
  }else{
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(12,12),
              symbreaks=FALSE, key=TRUE, symkey=FALSE,
              density.info="none", trace="none", main=main_title, cexRow=1.5,cexCol=1.2,col=colorRampPalette(c("blue","white", "red"))(n = 20),
              KeyValueName="-log10(Kd)",keysize = 0.8,xlab="Kinases",ylab=paste0(as.character(nrow(matrix_topvar))," Drugs"))
  }
  
  
}

##################################################################################################################################
#PLOTTING FUNCTIONS:
##################################################################################################################################

DREAM_kinome_Kds<-kinome_inhibitor_matrix_pKa_OVERLAP
rownames(DREAM_drugLibrary_annotation)<-toupper(rownames(DREAM_drugLibrary_annotation))


#MIX OF INHIBITORS:  some dissagreement with Drug-Banks
hcluster_euclid_n1screen_DrugTargets_kinasesDREAM(matrix=DREAM_kinome_Kds,
                                             database=DREAM_drugLibrary_annotation,
                                             target_vector = c("EGFR","BRAF","PLK1","FGFR1","ALK","KIT"),
                                             target_colors = c("darkred","red","blue3","purple","darkorange","steelblue"),
                                             main_title = "KINASE INHIBITORS:  Kinome Database & Drug-Bank disagree")




#3 classes: All aggree with Drug-Bank + variable off-targets
hcluster_euclid_n1screen_DrugTargets_kinasesDREAM(matrix=DREAM_kinome_Kds,
                                               database=DREAM_drugLibrary_annotation,
                                               target_vector = c("EGFR","FLT3","ABL1","RET","MET","AKT2"),
                                               target_colors = c("darkred","darkorange","steelblue","blue","darkgreen","red"),
                                               main_title = "KINASE INHIBITORS: Kinome Database & Drug-Bank Agreee")


##################################################################################################################################
#EXPORT DATA
##################################################################################################################################


####################################
#DEFINE DRUG-LIST
####################################
target_vector = c("EGFR","FLT3","ABL1","RET","MET","AKT2")
database=DREAM_drugLibrary_annotation
final_drugs<-c()
for (row in 1:length(target_vector)){
  target<-target_vector[row]
  
  #2 extract class-indices in rlab
  targets_list<-strsplit(database$target,",")
  db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
  class_names_tmp<-rownames(database[db_class_idx,])
  
  
  #3 make master index to filter matrix
  final_drugs<-c(final_drugs,class_names_tmp)
}
final_drugs_unique<-unique(final_drugs)



#HARD-WIRE DRUG LIST FROM ABOVE:
drug_list<-c("AEE788","AFATINIB","DACOMITINIB","GEFITINIB","ICOTINIB","LAPATINIB","NERATINIB","OSIMERTINIB","VANDETANIB","VARLITINIB","CABOZANTINIB","CEDIRANIB","CRENOLANIB","DOVITINIB","KW2449","LINIFANIB","PONATINIB","QUIZARTINIB","SORAFENIB","SUNITINIB","BAFETINIB","BOSUTINIB","DASATINIB","IMATINIB","NILOTINIB","REGORAFENIB","CRIZOTINIB","FORETINIB","MGCD265","TIVANTINIB","AZD5363","MK2206")




####################################
#FORMATTING RNAseq Matrix:
####################################
sort(table(sapply(strsplit(colnames(DREAM_RNAseq),"_"),function(x) x[length(x)])),decreasing=T)#ONLY DIFFERENCE SHOULD BE # DMSO & UNTREATED
#aspc hcc1143  hf2597   lncap   du145   efo21   panc1    hsts    krj1   h1793     u87 
#262     262     262     262     238     238     238     214     214     198     198 


sort(table(sapply(strsplit(colnames(DREAM_RNAseq),"_"),function(x) x[1])),decreasing=T)
#DMSO    UNTREATED       AEE788     Afatinib    Alisertib       AMG900       AT9283     Axitinib      AZD1480      AZD2014      AZD5363    Bafetinib       BI2536    Bosutinib     Brivanib 
#512          512           22           22           22           22           22           22           22           22           22           22           22           22           22 
#Cabozantinib    Cediranib    Ceritinib     CP724714   Crenolanib   Crizotinib   Dabrafenib  Dacomitinib    Dasatinib   Dinaciclib    Dovitinib  Enzastaurin    Foretinib Galunisertib    Gefitinib 
#22           22           22           22           22           22           22           22           22           22           22           22           22           22           22 
#GSK461364    Ibrutinib     Icotinib     Imatinib       KW2449    Lapatinib   Lenvatinib    Linifanib   Linsitinib      MGCD265       MK1775       MK2206      MLN2480  Momelotinib    Motesanib 
#22           22           22           22           22           22           22           22           22           22           22           22           22           22           22 
#Neratinib    Nilotinib   Nintedanib  Osimertinib    Pazopanib   PF04691502   Pimasertib    Ponatinib  Quizartinib       RAF265  Regorafenib   Rigosertib  Ruxolitinib  Saracatinib  Selumetinib 
#22           22           22           22           22           22           22           22           22           22           22           22           22           22           22 
#Sorafenib    Sunitinib       TAK733    Telatinib   Tivantinib    Tivozanib  Tofacitinib   Trametinib   Vandetanib   Varlitinib    Vatalanib  Vemurafenib   Volasertib 
#22           22           22           22           22           22           22           22           22           22           22           22           22



####################################
#EXPORTING DATA:  ANNOTATION
####################################
drug_list<-c("AEE788","AFATINIB","DACOMITINIB","GEFITINIB","ICOTINIB","LAPATINIB","NERATINIB","OSIMERTINIB","VANDETANIB","VARLITINIB","CABOZANTINIB","CEDIRANIB","CRENOLANIB","DOVITINIB","KW2449","LINIFANIB","PONATINIB","QUIZARTINIB","SORAFENIB","SUNITINIB","BAFETINIB","BOSUTINIB","DASATINIB","IMATINIB","NILOTINIB","REGORAFENIB","CRIZOTINIB","FORETINIB","MGCD265","TIVANTINIB","AZD5363","MK2206")


FINAL_DREAM_drugLibrary_annotation<-DREAM_drugLibrary_annotation[drug_list,-1]

write.csv(FINAL_DREAM_drugLibrary_annotation,file = "FINAL_DREAM_drugLibrary_annotation.csv")

####################################
#EXPORTING DATA:  ANNOTATION
####################################
drug_list<-c("AEE788","AFATINIB","DACOMITINIB","GEFITINIB","ICOTINIB","LAPATINIB","NERATINIB","OSIMERTINIB","VANDETANIB","VARLITINIB","CABOZANTINIB","CEDIRANIB","CRENOLANIB","DOVITINIB","KW2449","LINIFANIB","PONATINIB","QUIZARTINIB","SORAFENIB","SUNITINIB","BAFETINIB","BOSUTINIB","DASATINIB","IMATINIB","NILOTINIB","REGORAFENIB","CRIZOTINIB","FORETINIB","MGCD265","TIVANTINIB","AZD5363","MK2206")


FINAL_DREAM_kinome_Kds<-kinome_inhibitor_matrix_pKa_OVERLAP[drug_list,]

write.csv(FINAL_DREAM_kinome_Kds,file = "FINAL_DREAM_kinome_Kds.csv")



####################################
#EXPORTING DATA:  RNAseq csv + RData-list
####################################

drug_list<-c("AEE788","AFATINIB","DACOMITINIB","GEFITINIB","ICOTINIB","LAPATINIB","NERATINIB","OSIMERTINIB","VANDETANIB","VARLITINIB","CABOZANTINIB","CEDIRANIB","CRENOLANIB","DOVITINIB","KW2449","LINIFANIB","PONATINIB","QUIZARTINIB","SORAFENIB","SUNITINIB","BAFETINIB","BOSUTINIB","DASATINIB","IMATINIB","NILOTINIB","REGORAFENIB","CRIZOTINIB","FORETINIB","MGCD265","TIVANTINIB","AZD5363","MK2206")


DREAM_RNAseq_upper<-DREAM_RNAseq
colnames(DREAM_RNAseq_upper)<-gsub("ASPC","ASPC1",gsub("__","_0_",toupper(colnames(DREAM_RNAseq))))




drugs_in_order<-sapply(strsplit(colnames(DREAM_RNAseq_upper),"_"),function(x) x[1])
lines_in_order<-sapply(strsplit(colnames(DREAM_RNAseq_upper),"_"),function(x) x[4])
unique_lines<-names(table(sapply(strsplit(colnames(DREAM_RNAseq_upper),"_"),function(x) x[4])))
#names(table(sapply(strsplit(colnames(DREAM_RNAseq_upper),"_"),function(x) x[length(x)])))

#SLICE DOSE-RESPONSES:
FINAL_DREAM_RNAseq<-list()
#line="ASPC1"
for (line in unique_lines){
  drug_idx<-which(drugs_in_order %in% c(drug_list,"UNTREATED","DMSO"))
  line_idx<-which(lines_in_order == line)
  
  line_specific_drugs_idx<-intersect(drug_idx,line_idx)
  
  #SLICE/STORE dose-response slice:
  FINAL_DREAM_RNAseq[[line]]<-DREAM_RNAseq_upper[,line_specific_drugs_idx]
  
  
  tmp_matrix<-FINAL_DREAM_RNAseq[[line]]
  write.csv(tmp_matrix,file=paste0(line,"-RNAseq-Perturbations.csv"))
}


####################################
#EXPORTING DATA:  DOSE-RESPONSE csv + RData-list
####################################
drug_list<-c("AEE788","AFATINIB","DACOMITINIB","GEFITINIB","ICOTINIB","LAPATINIB","NERATINIB","OSIMERTINIB","VANDETANIB","VARLITINIB","CABOZANTINIB","CEDIRANIB","CRENOLANIB","DOVITINIB","KW2449","LINIFANIB","PONATINIB","QUIZARTINIB","SORAFENIB","SUNITINIB","BAFETINIB","BOSUTINIB","DASATINIB","IMATINIB","NILOTINIB","REGORAFENIB","CRIZOTINIB","FORETINIB","MGCD265","TIVANTINIB","AZD5363","MK2206")



#SLICE DOSE-RESPONSES:
FINAL_DREAM_Dose_Responses<-list()
line_names<-c("ASPC1","DU145","EFO21","H1793","HCC1143","HF2597","HSTS","KRJ1","LNCAP","PANC1","U87")
#line="ASPC1"
for (line in line_names){
  tmp_matrix<-DREAM_Dose_Responses[[line]]
  colnames(tmp_matrix)<-toupper(colnames(tmp_matrix))
  
  all_na_idx<-which(apply(tmp_matrix,1,function(x) sum(is.na(x)))==ncol(tmp_matrix))
  
  tmp_matrix_filter_na<-tmp_matrix[-all_na_idx,]
  
  tmp_matrix_FINAL<-tmp_matrix_filter_na[,drug_list]
    
    
  #SLICE/STORE dose-response slice:
  FINAL_DREAM_Dose_Responses[[line]]<-tmp_matrix_FINAL
  
  write.csv(tmp_matrix_FINAL,file=paste0(line,"-Dose-Responses.csv"))
}


save(list=c("drug_list","FINAL_DREAM_drugLibrary_annotation","FINAL_DREAM_kinome_Kds","FINAL_DREAM_RNAseq","FINAL_DREAM_Dose_Responses"),file="DREAM_kinaseInhibitor_FINAL-DATA.RData")





####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########      ########      ########     ########    #########    ####     ##  ##########                  #########        ##      ##########       ##                                    
##     ##     ##     ##     ##     ##          ##           ##           ## ##    ##      ##                      ##      ##      ####         ##          ####                              
##     ##     ##     ##    ##      ##          ##           ##           ##  ##   ##      ##                      ##       ##    ##  ##        ##         ##  ##                             
##     ########      #######       ########     #######     ########     ##   ##  ##      ##                      ##       ##   ##    ##       ##        ##    ##                                
##     ##            ##    ##      ##                 ##    ##           ##    ## ##      ##                      ##       ##  ##########      ##       ##########                                                                                             
##     ##            ##     ##     ##                 ##    ##           ##     ####      ##                      ##      ##   ##      ##      ##       ##      ##                           
##     ##            ##      ##    ########     #######     #########    ##      ###      ##                      #########    ##      ##      ##       ##      ##                                      
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#FOLLOW UP ON FINALIZED DATA:  10-16-2019



#setwd("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA")
setwd("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/")



load("DREAM_kinaseInhibitor_FINAL-DATA.RData")
###load("DREAM_kinaseInhibitor_ALL-DATA.RData")


load("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-07-15 - REDO PANGEA ANALYSIS/HALLMARKS_drugAVERAGE.RData")


##################################################################################################################################
#PLOTTING FUNCTIONS:
##################################################################################################################################

library("gplots")
library("RColorBrewer")
library("devtools")

heatmap.3 <- function(x,Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,distfun = dist,hclustfun = hclust,dendrogram = c("both","row", "column", "none"),symm = FALSE,scale = c("none","row", "column"),na.rm = TRUE,revC = identical(Colv,"Rowv"),add.expr,
                      breaks,symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",col = "heat.colors",colsep,rowsep,sepcolor = "white",sepwidth = c(0.05, 0.05),cellnote,notecex = 1,notecol = "cyan",na.color = par("bg"),
                      trace = c("none", "column","row", "both"),tracecol = "cyan",hline = median(breaks),vline = median(breaks),linecol = tracecol,margins = c(5,5),ColSideColors,RowSideColors,side.height.fraction=0.3,cexRow = 0.2 + 1/log10(nr),cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,labCol = NULL,key = TRUE,keysize = 1.5,density.info = c("none", "histogram", "density"),denscol = tracecol,symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,main = NULL,xlab = NULL,ylab = NULL,lmat = NULL,lhei = NULL,lwid = NULL,ColSideColorsSize = 1,RowSideColorsSize = 1,KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}



hcluster_corr_n1screen_DrugTargets_hallmarks<-function(matrix,database,target_vector,target_colors,main_title){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix_idx<-which(rownames(matrix) %in% rownames(database))
  
  matrix<-matrix[matrix_idx,]
  
  
  
  #FILTER MATRIX BY ONLY DRUGS ANNOTATED: (reduce noise of drugs that do nothing...)
  master_indx_filter<-c()
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #2 extract class-indices in rlab
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(rownames(matrix),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 make master index to filter matrix
    master_indx_filter<-c(master_indx_filter,rlab_class_idx)
  }
  master_indx_filter_unique<-unique(master_indx_filter)
  
  matrix<-matrix[master_indx_filter_unique,]
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(target_vector),ncol=length(rownames(matrix)))
  rownames(rlab)<-target_vector
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #1 extract drugs-vector from n1screen_drugDB
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(colnames(rlab),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-target_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  top_var<-names(sort(apply(matrix,2,var),decreasing=T))[1:20]
  matrix_topvar<-matrix[,top_var]
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  cor_columns <- cor(matrix_topvar,method="spearman")
  cor_rows <- cor(t(matrix_topvar),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(length(target_vector)>1){
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(20,12),
              RowSideColors=rlab, symbreaks=TRUE, key=TRUE, symkey=TRUE,
              density.info="none", trace="none", main=main_title, cexRow=1.5,cexCol=1.2,col=colorRampPalette(c("blue","blue", "white", "red", "red"))(n = 20),
              RowSideColorsSize=5, KeyValueName="NES",keysize = 0.8,xlab="Hallmarks",ylab="Drug-Averages (across cell-lines)")
  }else{
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(20,12),
              symbreaks=TRUE, key=TRUE, symkey=TRUE,
              density.info="none", trace="none", main=main_title, cexRow=1.5,cexCol=1.2,col=colorRampPalette(c("blue","blue", "white", "red", "red"))(n = 20),
              KeyValueName="NES",keysize = 0.8,xlab="Hallmarks",ylab="Drug-Averages (across cell-lines)")
  }
  
  
}



hcluster_euclid_n1screen_DrugTargets_kinasesDREAM<-function(matrix,database,target_vector,target_colors,main_title,num_var=30){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix_idx<-which(rownames(matrix) %in% rownames(database))
  
  matrix<-matrix[matrix_idx,]
  
  
  
  #FILTER MATRIX BY ONLY DRUGS ANNOTATED: (reduce noise of drugs that do nothing...)
  master_indx_filter<-c()
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #2 extract class-indices in rlab
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(rownames(matrix),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 make master index to filter matrix
    master_indx_filter<-c(master_indx_filter,rlab_class_idx)
  }
  master_indx_filter_unique<-unique(master_indx_filter)
  
  matrix<-matrix[master_indx_filter_unique,]
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(target_vector),ncol=length(rownames(matrix)))
  rownames(rlab)<-target_vector
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #1 extract drugs-vector from n1screen_drugDB
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(colnames(rlab),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-target_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  top_var<-names(sort(apply(matrix,2,var),decreasing=T))[1:num_var]
  matrix_topvar<-matrix[,top_var]
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  
  dist_columns <- dist(t(matrix_topvar))
  dist_rows <- dist(matrix_topvar)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(length(target_vector)>1){
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(12,12),
              RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
              density.info="none", trace="none", main=main_title, cexRow=1.5,cexCol=1.2,col=colorRampPalette(c("blue","white", "red"))(n = 20),
              RowSideColorsSize=5, KeyValueName="-log10(Kd)",keysize = 0.8,xlab="Kinases",ylab=paste0(as.character(nrow(matrix_topvar))," Drugs"))
  }else{
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(12,12),
              symbreaks=FALSE, key=TRUE, symkey=FALSE,
              density.info="none", trace="none", main=main_title, cexRow=1.5,cexCol=1.2,col=colorRampPalette(c("blue","white", "red"))(n = 20),
              KeyValueName="-log10(Kd)",keysize = 0.8,xlab="Kinases",ylab=paste0(as.character(nrow(matrix_topvar))," Drugs"))
  }
  
  
}

##################################################################################################################################
#PLOTTING FUNCTIONS:
##################################################################################################################################

#MIX OF INHIBITORS:  some dissagreement with Drug-Banks
hcluster_euclid_n1screen_DrugTargets_kinasesDREAM(matrix=FINAL_DREAM_kinome_Kds,
                                                  database=FINAL_DREAM_drugLibrary_annotation,
                                                target_vector = c("EGFR","FLT3","ABL1","RET","MET","AKT2"),
                                                target_colors = c("darkred","darkorange","steelblue","blue","darkgreen","red"),
                                                main_title = "KINASE INHIBITORS: Kinome Database & Drug-Bank Agreee")








####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##    ########     ########    ##     ##     ########              ###############     ##       ########         #########    #########   ###########      #######                                                                                                                        
##    ##     ##    ##     ##   ##     ##    ##                           ##           ####      ##     ##       ##            ##               ##         ##                                   
##    ##      ##   ##     ##   ##     ##   ##                            ##          ##  ##     ##     ##      ##             ##               ##         ##                                                      
##    ##      ##   ########    ##     ##   ##    ####    ######          ##         ##    ##    ########       ##     ####    #######          ##          #######                             
##    ##      ##   ##    ##    ##     ##   ##       ##                   ##        ##########   ##    ##       ##        ##   ##               ##                ##                                              
##    ##     ##    ##     ##   ##     ##    ##      ##                   ##        ##      ##   ##     ##       ##      ##    ##               ##                ##                                                 
##    ########     ##      ##   #######      ########                    ##        ##      ##   ##      ##       ########     ##########       ##          #######                                                                            
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#GOAL:  Define:  "all drugable targets" according to drug-bank + 600 kinome

setwd("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA")

#ALL DRUG-BANK
load("/Volumes/ifs/archive/c2b2/ac_lab/efd2115/00 HOME ARCHIVE/02 Research/02 Coding-Data/2017/2017-06-21 - N1-benchmarks - drug-MR-mut/03 DrugBank - N1 database/DrugBank_df.RData")
##load("/Volumes/ifs/archive/c2b2/ac_lab/efd2115/00 HOME ARCHIVE/02 Research/02 Coding-Data/2017/2017-06-21 - N1-benchmarks - drug-MR-mut/03 DrugBank - N1 database/DrugBank_AClab.RData")

#ENTREZ2HUGO
load("/Volumes/ifs/archive/shares/n-of-1-plateseq-data/n1screen_paper/n1database_org/context_annotations/entrez2hugo_LIST.RData")

#LOAD KINOME CLASSIFICATIONS:
load("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-06-03 KINASE_targets_effectors/KINOME-classifications.RData")


####################
#1 clean up DRUG-BANK "target column"
####################

#remove text surrounding HUGO id:
targets_pruned<-gsub(".*GenAtlas","",DrugBank_df[,"targets"])
targets_pruned<-gsub("GenBank.*","",targets_pruned)

#remove text surrounding "DNA" target info:
targets_pruned<-gsub("Human.*","",targets_pruned)
targets_pruned<-gsub("BE.......","",targets_pruned)

#UNIQUE
targets_pruned_unique<-unique(targets_pruned)

#OVERLAP w/ HUGO-names:
targets_pruned_unique_hugo<-targets_pruned_unique[which(targets_pruned_unique %in% names(hugo2entrez))]

#UNIQUE:

####################
#2 ADD KINOME-GENES TO DRUG-BANK TARGETS:
####################
kinome_gene_names<-rownames(kinome_classification_hugo_filtered)

kinome_gene_names_hugo<-kinome_gene_names[which(kinome_gene_names %in% names(hugo2entrez))]

FINAL_TARGETLIST_DREAM<-sort(unique(c(kinome_gene_names_hugo,targets_pruned_unique_hugo)))

save(FINAL_TARGETLIST_DREAM,file="FINAL_TARGETLIST_DREAM.RData")


write.csv(FINAL_TARGETLIST_DREAM,file="drugBank_kinome_list.csv")























####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################  
##      ########     #########    ##########        ##                                                               
##     ##           ##            ##               ####                                                         
##    ##            ##            ##              ##  ##                                                         
##    ##    ####     ########     ########       ##    ##                                                            
##    ##       ##           ##    ##            ##########                                                         
##     ##      ##           ##    ##           ##        ##                                                         
##      ########     ########     ##########  ##          ##          
####################################################################################################################################################################################################################################################################  
####################################################################################################################################################################################################################################################################  

setwd("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA")

#INPUT: matrix of (drug x genes) or (drugs x lines)

library("viper")
library("gplots")
library("RColorBrewer", lib.loc="~/Library/R/3.3/library")
library("devtools")
heatmap.3 <- function(x,Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,distfun = dist,hclustfun = hclust,dendrogram = c("both","row", "column", "none"),symm = FALSE,scale = c("none","row", "column"),na.rm = TRUE,revC = identical(Colv,"Rowv"),add.expr,
                      breaks,symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",col = "heat.colors",colsep,rowsep,sepcolor = "white",sepwidth = c(0.05, 0.05),cellnote,notecex = 1,notecol = "cyan",na.color = par("bg"),
                      trace = c("none", "column","row", "both"),tracecol = "cyan",hline = median(breaks),vline = median(breaks),linecol = tracecol,margins = c(5,5),ColSideColors,RowSideColors,side.height.fraction=0.3,cexRow = 0.2 + 1/log10(nr),cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,labCol = NULL,key = TRUE,keysize = 1.5,density.info = c("none", "histogram", "density"),denscol = tracecol,symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,main = NULL,xlab = NULL,ylab = NULL,lmat = NULL,lhei = NULL,lwid = NULL,ColSideColorsSize = 1,RowSideColorsSize = 1,KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}


#Asymmetric Key: log10(Kd)
hcluster_corr_MULTI_sidebar_G2M_mod2_nonsym<-function(matrix,database,DB_col,DB_class,class_names,class_colors,main_title,G2M_plus=TRUE,x_lab="Genes"){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  #0 RENAME MATRIX to names compatible with database 
  rownames(database)<-make.names(tolower(rownames(database)),unique = TRUE)
  
  rownames(matrix)<-make.names(tolower(rownames(matrix)),unique = TRUE)
  
  
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix<-matrix[rownames(database),]
  
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(class_names),ncol=length(rownames(matrix)))
  rownames(rlab)<-class_names
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(class_names)){
    #1 extract 1st rows DB-entry
    db_col_tmp<-DB_col[row]
    db_class_tmp<-DB_class[row]
    
    #2 extract class-indices in rlab
    db_class_idx<-which(database[,db_col_tmp] == db_class_tmp)
    class_names_tmp<-rownames(database[db_class_idx,])
    rlab_class_idx<-which(colnames(rlab) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-class_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  
  
  #########################################
  #XX G2M TOP-bar
  #########################################
  
  ###1 extract relevant G2M genes:
  G2M_chckpt<-c("TOP2A","CENPF","CHEK1","FOXM1","CENPI","BLM","BRIP1","MCM6","MCM4","E2F1","E2F2","E2F3","E2F8","BUB1B","MYBL2","MCM7","FEN1","AURKA","AURKB","CENPU","FANCA","TYMS","RAD51","HMGA1","DNMT3B","E2F7","BRCA1","BRCA2","PRKDC","CENPK","MCM3","DMNT1","FANCD2","RAN","RANBP1","MCM2","CENPZ","RACGAP1","PCNA","H2AFX","MCM5")
  G2M_chckpt_plus<-unique(c("TOP2A","CENPF","CHEK1","FOXM1","CENPI","BLM","BRIP1","MCM6","MCM4","E2F1","E2F2","E2F3","E2F8","BUB1B","MYBL2","MCM7","FEN1","AURKA","AURKB","CENPU","FANCA","TYMS","RAD51","HMGA1","DNMT3B","E2F7","BRCA1","BRCA2","PRKDC","CENPK","MCM3","DMNT1","FANCD2","RAN","RANBP1","MCM2","CENPZ","RACGAP1","PCNA","H2AFX","MCM5","TRIP13","ASF1B","UHRF1","ATAD2","ECT2","TTK","WHSC1","ARHGAP11A","ZNF367","TCF19","HELLS","RACGAP1","MKI67","SPAG5","UBE2C","TONSL","RUVBLT","ZWINT","HMGB2","DEPDC1","CBX3","RFC4","CHAF1A","GTSE1","HMGN1","PSMC3IP","HDGF","CKS1B","TCEB1","CKS2","CCNA2","ENY2","GGCT","ILE2","HSPE1","PSRC1","NDC80","TMPO","CDK2","GMPS","EIF4EBP1","PA2G4","GMNN","ZC3H15","PTTG1","ZC3H1","CYCS","GTF2E2","HLTF","TCF19","WHSC1","TTF2","SHOX2","EIF2AK2","ARNTL2","UHRF1","GTSE1","DEPDC1","SUV39H2","ZWINT","YEASTS2","LRPPRC","NAA15","POLD1","RFC4","RBL1","MCM3","ILF3","CDK7","CIT","GLRA2","PHF19","CDCA7","ZNF695","EZH2","PCNA","HNRNPAB","ZCRB1","SUB1","FH","HSPE1","RAB1A","CCDC47","COPS3","H2AFX","DAXX","HDGF","APH1A","PUF60","POLR3K","HSBP1","IDH2","VPS72","PRMT1","UBE2L3","PPP2CA","HSPA5","MCM2","HDAC2","CEBPZ","TIMELESS","UBEC2","PTTG1","TRAIP","MNX1","ARHGEF39","CHEK2","IQGAP3","HSPD1"))
  if (G2M_plus == TRUE){
    G2M_chckpt_filtered<-G2M_chckpt_plus[which(G2M_chckpt_plus %in% colnames(matrix))]
  }
  else {
    G2M_chckpt_filtered<-G2M_chckpt[which(G2M_chckpt %in% colnames(matrix))]
  }
  
  ###2 MAKE column color-key:
  clab<-matrix("white",nrow=length(colnames(matrix)),ncol=1)
  
  matrix_G2M_idx<-which(colnames(matrix) %in% G2M_chckpt_filtered)
  
  clab[matrix_G2M_idx,1]<-"darkblue"
  
  colnames(clab)=c("G2M")
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  cor_columns <- cor(matrix,method="spearman")
  cor_rows <- cor(t(matrix),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  heatmap.3(matrix, 
            Rowv=as.dendrogram(hclust_rows),
            Colv=as.dendrogram(hclust_columns), 
            na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,6),
            ColSideColors=clab, RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
            density.info="none", trace="none", main=main_title, cexRow=0.6,cexCol=0.3,col=colorRampPalette(c("blue", "white", "red","red"))(n = 20),
            ColSideColorsSize=1, RowSideColorsSize=5, KeyValueName="spearman",keysize = 0.8,xlab=x_lab,ylab="Drugs")
  
}



GSEA_kinase_predictions<-function(kinase_inhibitors_BINDING2,kinase_prediction_matrix_na,threshold=6,title_main="True-Targets w/n Predictions"){
  
  #############################################################################################
  #1 FORMAT DATA SO CAN DO aREA ENRICHEMENT EASILY:
  #############################################################################################
  
  
  ########################
  #BINARY-IZE AUC RESPONSE: i.e. define "GOLD STANDARDS"
  ########################
  
  stat_sig_drug_idx<-which(kinase_inhibitors_BINDING2 > threshold)
  
  kinase_inhibitors_BINDING2_binary<-kinase_inhibitors_BINDING2
  kinase_inhibitors_BINDING2_binary[stat_sig_drug_idx]<-1
  kinase_inhibitors_BINDING2_binary[-stat_sig_drug_idx]<-0
  
  
  percent_sensitive_drugs=length(stat_sig_drug_idx)/(nrow(kinase_inhibitors_BINDING2_binary)*ncol(kinase_inhibitors_BINDING2_binary))
  #4.7%
  
  
  
  
  ########################
  #MAKE MATRIX OF RANKED DRUGS / CELL-LINE:  rank  x cell-line
  ########################
  kinase_PREDICTIONS_formatted<-t(kinase_prediction_matrix_na)
  
  
  
  
  ##############
  #NOTE: more positive score = better binding (i.e. predictions convert VIPER into -log10(pKd))
  ##############
  #sort(kinase_PREDICTIONS_formatted[,"ERLOTINIB"],decreasing=T)[1:10]
  #MAP2K2   MAP2K1     PFKP     BRAF      NLK    TESK1      KIT    IGF1R     INSR     ARAF 
  #260.4869 231.2391 194.9972 166.3905 143.4486 143.4486 143.4486 129.8185 129.8185 125.8992 
  
  #MAKE RANKED MATRIX:
  num_kinases<-length(rownames(kinase_PREDICTIONS_formatted))
  
  kinaseRanks_per_drug<-matrix(nrow=num_kinases,ncol=0)
  
  rownames(kinaseRanks_per_drug)<-1:num_kinases
  
  col_names<-c()
  for (column in colnames(kinase_PREDICTIONS_formatted)){
    #1 store cell-line name:
    col_names<-append(col_names,column)
    #2 organize kinases by rank
    temp_ranked_kinases<-names(sort(kinase_PREDICTIONS_formatted[,column],decreasing=TRUE))
    #3 append ranked kinase-list matrix:
    kinaseRanks_per_drug<-cbind(kinaseRanks_per_drug,temp_ranked_kinases)
    #4 name columns:
    colnames(kinaseRanks_per_drug)<-col_names
  }
  
  
  
  
  ######################################################################################################################################################################################
  #PLOTTING AND ANALYSIS   (original function)
  ######################################################################################################################################################################################
  gene_set_matrix=kinase_inhibitors_BINDING2_binary#BINARY GOLD STANDARD DEFINITIONS = define regulon
  gene_all_matrix=kinaseRanks_per_drug #KINASE-RANKING/LINE = for enrichment visualization
  gene_all_data_matrix=kinase_PREDICTIONS_formatted#KINASE-SCORE/DRUG 
  
  
  
  ###################
  #1 MAKE PLOT
  ###################
  blank_enrichment_plot_matrix<-matrix(1,ncol=ncol(gene_all_matrix),nrow=nrow(gene_all_matrix))
  colnames(blank_enrichment_plot_matrix)<-colnames(gene_all_matrix)
  
  #column="AFATINIB"
  for (column in colnames(blank_enrichment_plot_matrix)){
    #1 extract gene-set w/ 1 in matrix:
    column_gene_set<-names(which(gene_set_matrix[,column]==1))
    #[1] "EGFR"     "MAPKAPK2"
    
    #3 make plottable matrix:
    gene_set_idx_wn_plot_column<-which(gene_all_matrix[,column] %in% column_gene_set)
    blank_enrichment_plot_matrix[gene_set_idx_wn_plot_column,column]<-2
  }
  
  visualize_plot <- blank_enrichment_plot_matrix[ nrow(blank_enrichment_plot_matrix):1, ]
  
  ###################
  #2 PLOT ENRICMENT:
  ###################
  image(x=1:ncol(visualize_plot),y=1:nrow(visualize_plot),z=t(visualize_plot),axes=FALSE,xlab="kinase inhibitors",ylab="ranked Kinases",col=c("white","red"))
  par(cex.lab=1,cex.axis=1)
  axis(2,at=seq(from=228, to=0, by = -10), labels = TRUE)
  #axis(2,at=seq(3, nrow(visualize_plot)-10, by = 10), labels = rev(as.character(seq(10, 90, by = 10))))#LOOKS LIKE I NEED TO MANUALLY ADJUST AXIS TO MATCH VISUALIZATION...
  par(cex.lab=1,cex.axis=min(1,10/ncol(visualize_plot)))
  axis(1,at=1:ncol(visualize_plot), labels = colnames(visualize_plot),srt=10,tick=FALSE)
  
  box(lwd=2)
  abline(v=1:ncol(visualize_plot)-0.5,lwd=1)
  abline(h=floor(nrow(visualize_plot)/2),lwd=1,lty=1)
  
  #title(main =title, font.main = 8)
  
  
  ####################
  #3 CALC aREA ENRICHMENT
  ####################
  #FIGURE OUT WHICH LINES HAVE >1 drug that they are sensitive too...
  num_sel_drugs_per_line<-sort(colSums(gene_set_matrix),decreasing=TRUE)
  regulon_sufficient_lines<-names(which(num_sel_drugs_per_line>1))
  
  
  #MAKE REGULON: for ONLY cell-lines w/ >1 drug
  library(viper)
  gene_set_regulon<-list()
  for (column in regulon_sufficient_lines){
    column_gene_set<-names(which(gene_set_matrix[,column]==1))
    num_genes_wn_set<-length(column_gene_set)
    
    gene_set_regulon[[column]]<-list()
    gene_set_regulon[[column]]$tfmode<-rep(1,num_genes_wn_set)
    gene_set_regulon[[column]]$likelihood<-rep(1,num_genes_wn_set)
    names(gene_set_regulon[[column]]$tfmode)<-column_gene_set
  }
  
  #gene_set_regulon$AFATINIB
  #$tfmode
  #EGFR MAPKAPK2 
  #1        1 
  #$likelihood
  #[1] 1 1
  
  #RUN aREA (use ONLY >1 drug-line-NETWORK and >1 drug-line-DATA)
  aREA_all<-aREA(gene_all_data_matrix[,regulon_sufficient_lines],gene_set_regulon,minsize=2)
  
  #extract aREA on self-predictions:
  #View(aREA_all$nes)
  aREA_all_diag<-c()
  for (i in 1:nrow(aREA_all$nes)){
    aREA_all_diag<-append(aREA_all_diag,aREA_all$nes[i,i])
  }
  names(aREA_all_diag)<-rownames(aREA_all$nes)
  
  #sort aREA:
  sorted_aREA_vector<-sort(aREA_all_diag,decreasing=TRUE)
  sorted_aREA_names<-names(sorted_aREA_vector)
  
  #stouffer integration:
  integrated_p<-sum(aREA_all_diag)/sqrt(length(aREA_all_diag))
  
  sorted_visualize_plot<-visualize_plot[,sorted_aREA_names]
  
  ##################
  #PLOT aREA vs RANDOM-null:
  ##################
  ##NULL distribution: (enrichment for random combinations)
  hist(aREA_all$nes, col=rgb(0,0,0,0.2),main=paste0("NES for true Kinase Targets w/n Predictions (stouffer nes: ",as.character(round(integrated_p,3)),")"),xlab="NES(aREA)",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,freq = FALSE,breaks=100,ylim=c(0,0.5))
  #paired distribution: (top achilles wn VIPER)
  hist(sorted_aREA_vector, col=rgb(1,0,0,0.5),freq = FALSE, add=T,breaks=100,ylim=c(0,0.5))
  box()
  legend("topright",legend=c("NULL","Drug-Predictions"),fill=c(rgb(0,0,0,0.2), rgb(1,0,0,0.5)),cex=1.5)
  
  
  
  
  
  ###################
  #3 PLOT ENRICMENT:  sorted by aREA
  ###################
  image(x=1:ncol(sorted_visualize_plot),y=1:nrow(sorted_visualize_plot),z=t(sorted_visualize_plot),axes=FALSE,xlab="kinase inhibitors",ylab="ranked Kinases",col=c("white","red"))
  par(cex.lab=1,cex.axis=1)
  #axis(2,at=seq(from=93, to=3, by = -10), labels = TRUE)
  axis(2)#true rank inverted here...
  par(cex.lab=1,cex.axis=min(1,10/ncol(sorted_visualize_plot)))
  axis(1,at=1:ncol(sorted_visualize_plot), labels = colnames(sorted_visualize_plot),srt=10,tick=FALSE)
  
  box(lwd=2)
  abline(v=1:ncol(sorted_visualize_plot)-0.5,lwd=1)
  abline(h=floor(nrow(sorted_visualize_plot)/2),lwd=1,lty=1)
  
  title(main =paste0(title_main," (stouffer nes: ",as.character(round(integrated_p,3)),")"), font.main = 8)
  
  return(sorted_aREA_names)
}




################################################################################################
#LOAD DATA:
################################################################################################

load("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-05-24 - N1screen - drug-clean-ness/OVERLAP_n1-Kinome-Screens.RData")
rm(kinome_selectivity_scores_OVERLAP)


load("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-09-01 - REGULON-improvement/PREDICTIONS_kinome-binding-GES-LASSO.RData")

sorted_aREA_names<-GSEA_kinase_predictions(kinase_inhibitors_BINDING2,kinase_prediction_gesLASSO,threshold=6)

sorted_aREA_names<-GSEA_kinase_predictions(kinase_inhibitors_BINDING2,t(kinase_inhibitors_BINDING2),threshold=6)



####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################  
##      ##           #############                    ########       ########       ##       ###     ##   ###     ##    ########  ########                                                        
##      ##                ##                         ##             ##             ####      ####    ##   ####    ##    ##        ##     ##                   
##      ##                ##                         ##            ##             ##  ##     ## ##   ##   ## ##   ##    ##        ##     ##                  
##      ##                ##          ########        ########     ##            ##    ##    ##  ##  ##   ##  ##  ##    #######   ########                                          
##      ##                ##                                 ##    ##           ##########   ##   ## ##   ##   ## ##    ##        ##    ##                       
##      ##                ##                                 ##     ##          ##      ##   ##    ####   ##    ####    ##        ##     ##              
##      ##########        ##                          ########       ########   ##      ##   ##     ###   ##     ###    ######### ##      ##
####################################################################################################################################################################################################################################################################  
####################################################################################################################################################################################################################################################################  


setwd("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/")

#LOAD PARED LT SCANNER DATA:  parsed in #0 below
load("LTscanner_parsed-data.RData")

load("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-05-24 - N1screen - drug-clean-ness/OVERLAP_n1-Kinome-Screens.RData")
rm(kinome_selectivity_scores_OVERLAP)


load("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-09-01 - REGULON-improvement/PREDICTIONS_kinome-binding-GES-LASSO.RData")

#############################################################################################################
#PLOTTING FUNCTIONS
#############################################################################################################

library("viper")
library("gplots")
library("RColorBrewer", lib.loc="~/Library/R/3.3/library")
library("devtools")
heatmap.3 <- function(x,Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,distfun = dist,hclustfun = hclust,dendrogram = c("both","row", "column", "none"),symm = FALSE,scale = c("none","row", "column"),na.rm = TRUE,revC = identical(Colv,"Rowv"),add.expr,
                      breaks,symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",col = "heat.colors",colsep,rowsep,sepcolor = "white",sepwidth = c(0.05, 0.05),cellnote,notecex = 1,notecol = "cyan",na.color = par("bg"),
                      trace = c("none", "column","row", "both"),tracecol = "cyan",hline = median(breaks),vline = median(breaks),linecol = tracecol,margins = c(5,5),ColSideColors,RowSideColors,side.height.fraction=0.3,cexRow = 0.2 + 1/log10(nr),cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,labCol = NULL,key = TRUE,keysize = 1.5,density.info = c("none", "histogram", "density"),denscol = tracecol,symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,main = NULL,xlab = NULL,ylab = NULL,lmat = NULL,lhei = NULL,lwid = NULL,ColSideColorsSize = 1,RowSideColorsSize = 1,KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}


#Asymmetric Key: log10(Kd)
hcluster_corr_MULTI_sidebar_G2M_mod2_nonsym<-function(matrix,database,DB_col,DB_class,class_names,class_colors,main_title,G2M_plus=TRUE,x_lab="Genes"){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  #0 RENAME MATRIX to names compatible with database 
  rownames(database)<-make.names(tolower(rownames(database)),unique = TRUE)
  
  rownames(matrix)<-make.names(tolower(rownames(matrix)),unique = TRUE)
  
  
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix<-matrix[rownames(database),]
  
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(class_names),ncol=length(rownames(matrix)))
  rownames(rlab)<-class_names
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(class_names)){
    #1 extract 1st rows DB-entry
    db_col_tmp<-DB_col[row]
    db_class_tmp<-DB_class[row]
    
    #2 extract class-indices in rlab
    db_class_idx<-which(database[,db_col_tmp] == db_class_tmp)
    class_names_tmp<-rownames(database[db_class_idx,])
    rlab_class_idx<-which(colnames(rlab) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-class_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  
  
  #########################################
  #XX G2M TOP-bar
  #########################################
  
  ###1 extract relevant G2M genes:
  G2M_chckpt<-c("TOP2A","CENPF","CHEK1","FOXM1","CENPI","BLM","BRIP1","MCM6","MCM4","E2F1","E2F2","E2F3","E2F8","BUB1B","MYBL2","MCM7","FEN1","AURKA","AURKB","CENPU","FANCA","TYMS","RAD51","HMGA1","DNMT3B","E2F7","BRCA1","BRCA2","PRKDC","CENPK","MCM3","DMNT1","FANCD2","RAN","RANBP1","MCM2","CENPZ","RACGAP1","PCNA","H2AFX","MCM5")
  G2M_chckpt_plus<-unique(c("TOP2A","CENPF","CHEK1","FOXM1","CENPI","BLM","BRIP1","MCM6","MCM4","E2F1","E2F2","E2F3","E2F8","BUB1B","MYBL2","MCM7","FEN1","AURKA","AURKB","CENPU","FANCA","TYMS","RAD51","HMGA1","DNMT3B","E2F7","BRCA1","BRCA2","PRKDC","CENPK","MCM3","DMNT1","FANCD2","RAN","RANBP1","MCM2","CENPZ","RACGAP1","PCNA","H2AFX","MCM5","TRIP13","ASF1B","UHRF1","ATAD2","ECT2","TTK","WHSC1","ARHGAP11A","ZNF367","TCF19","HELLS","RACGAP1","MKI67","SPAG5","UBE2C","TONSL","RUVBLT","ZWINT","HMGB2","DEPDC1","CBX3","RFC4","CHAF1A","GTSE1","HMGN1","PSMC3IP","HDGF","CKS1B","TCEB1","CKS2","CCNA2","ENY2","GGCT","ILE2","HSPE1","PSRC1","NDC80","TMPO","CDK2","GMPS","EIF4EBP1","PA2G4","GMNN","ZC3H15","PTTG1","ZC3H1","CYCS","GTF2E2","HLTF","TCF19","WHSC1","TTF2","SHOX2","EIF2AK2","ARNTL2","UHRF1","GTSE1","DEPDC1","SUV39H2","ZWINT","YEASTS2","LRPPRC","NAA15","POLD1","RFC4","RBL1","MCM3","ILF3","CDK7","CIT","GLRA2","PHF19","CDCA7","ZNF695","EZH2","PCNA","HNRNPAB","ZCRB1","SUB1","FH","HSPE1","RAB1A","CCDC47","COPS3","H2AFX","DAXX","HDGF","APH1A","PUF60","POLR3K","HSBP1","IDH2","VPS72","PRMT1","UBE2L3","PPP2CA","HSPA5","MCM2","HDAC2","CEBPZ","TIMELESS","UBEC2","PTTG1","TRAIP","MNX1","ARHGEF39","CHEK2","IQGAP3","HSPD1"))
  if (G2M_plus == TRUE){
    G2M_chckpt_filtered<-G2M_chckpt_plus[which(G2M_chckpt_plus %in% colnames(matrix))]
  }
  else {
    G2M_chckpt_filtered<-G2M_chckpt[which(G2M_chckpt %in% colnames(matrix))]
  }
  
  ###2 MAKE column color-key:
  clab<-matrix("white",nrow=length(colnames(matrix)),ncol=1)
  
  matrix_G2M_idx<-which(colnames(matrix) %in% G2M_chckpt_filtered)
  
  clab[matrix_G2M_idx,1]<-"darkblue"
  
  colnames(clab)=c("G2M")
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  cor_columns <- cor(matrix,method="spearman")
  cor_rows <- cor(t(matrix),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  heatmap.3(matrix, 
            Rowv=as.dendrogram(hclust_rows),
            Colv=as.dendrogram(hclust_columns), 
            na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,6),
            ColSideColors=clab, RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
            density.info="none", trace="none", main=main_title, cexRow=0.6,cexCol=0.3,col=colorRampPalette(c("blue", "white", "red","red"))(n = 20),
            ColSideColorsSize=1, RowSideColorsSize=5, KeyValueName="spearman",keysize = 0.8,xlab=x_lab,ylab="Drugs")
  
}



GSEA_kinase_predictions<-function(kinase_inhibitors_BINDING2,kinase_prediction_matrix_na,threshold=6,title_main="True-Targets w/n Predictions"){
  
  #############################################################################################
  #1 FORMAT DATA SO CAN DO aREA ENRICHEMENT EASILY:
  #############################################################################################
  
  
  ########################
  #BINARY-IZE AUC RESPONSE: i.e. define "GOLD STANDARDS"
  ########################
  
  stat_sig_drug_idx<-which(kinase_inhibitors_BINDING2 > threshold)
  
  kinase_inhibitors_BINDING2_binary<-kinase_inhibitors_BINDING2
  kinase_inhibitors_BINDING2_binary[stat_sig_drug_idx]<-1
  kinase_inhibitors_BINDING2_binary[-stat_sig_drug_idx]<-0
  
  
  percent_sensitive_drugs=length(stat_sig_drug_idx)/(nrow(kinase_inhibitors_BINDING2_binary)*ncol(kinase_inhibitors_BINDING2_binary))
  #4.7%
  
  
  
  
  ########################
  #MAKE MATRIX OF RANKED DRUGS / CELL-LINE:  rank  x cell-line
  ########################
  kinase_PREDICTIONS_formatted<-t(kinase_prediction_matrix_na)
  
  
  
  
  ##############
  #NOTE: more positive score = better binding (i.e. predictions convert VIPER into -log10(pKd))
  ##############
  #sort(kinase_PREDICTIONS_formatted[,"ERLOTINIB"],decreasing=T)[1:10]
  #MAP2K2   MAP2K1     PFKP     BRAF      NLK    TESK1      KIT    IGF1R     INSR     ARAF 
  #260.4869 231.2391 194.9972 166.3905 143.4486 143.4486 143.4486 129.8185 129.8185 125.8992 
  
  #MAKE RANKED MATRIX:
  num_kinases<-length(rownames(kinase_PREDICTIONS_formatted))
  
  kinaseRanks_per_drug<-matrix(nrow=num_kinases,ncol=0)
  
  rownames(kinaseRanks_per_drug)<-1:num_kinases
  
  col_names<-c()
  for (column in colnames(kinase_PREDICTIONS_formatted)){
    #1 store cell-line name:
    col_names<-append(col_names,column)
    #2 organize kinases by rank
    temp_ranked_kinases<-names(sort(kinase_PREDICTIONS_formatted[,column],decreasing=TRUE))
    #3 append ranked kinase-list matrix:
    kinaseRanks_per_drug<-cbind(kinaseRanks_per_drug,temp_ranked_kinases)
    #4 name columns:
    colnames(kinaseRanks_per_drug)<-col_names
  }
  
  
  
  
  ######################################################################################################################################################################################
  #PLOTTING AND ANALYSIS   (original function)
  ######################################################################################################################################################################################
  gene_set_matrix=kinase_inhibitors_BINDING2_binary#BINARY GOLD STANDARD DEFINITIONS = define regulon
  gene_all_matrix=kinaseRanks_per_drug #KINASE-RANKING/LINE = for enrichment visualization
  gene_all_data_matrix=kinase_PREDICTIONS_formatted#KINASE-SCORE/DRUG 
  
  
  
  ###################
  #1 MAKE PLOT
  ###################
  blank_enrichment_plot_matrix<-matrix(1,ncol=ncol(gene_all_matrix),nrow=nrow(gene_all_matrix))
  colnames(blank_enrichment_plot_matrix)<-colnames(gene_all_matrix)
  
  #column="AFATINIB"
  for (column in colnames(blank_enrichment_plot_matrix)){
    #1 extract gene-set w/ 1 in matrix:
    column_gene_set<-names(which(gene_set_matrix[,column]==1))
    #[1] "EGFR"     "MAPKAPK2"
    
    #3 make plottable matrix:
    gene_set_idx_wn_plot_column<-which(gene_all_matrix[,column] %in% column_gene_set)
    blank_enrichment_plot_matrix[gene_set_idx_wn_plot_column,column]<-2
  }
  
  visualize_plot <- blank_enrichment_plot_matrix[ nrow(blank_enrichment_plot_matrix):1, ]
  
  ###################
  #2 PLOT ENRICMENT:
  ###################
  image(x=1:ncol(visualize_plot),y=1:nrow(visualize_plot),z=t(visualize_plot),axes=FALSE,xlab="kinase inhibitors",ylab="ranked Kinases",col=c("white","red"))
  par(cex.lab=1,cex.axis=1)
  axis(2,at=seq(from=228, to=0, by = -10), labels = TRUE)
  #axis(2,at=seq(3, nrow(visualize_plot)-10, by = 10), labels = rev(as.character(seq(10, 90, by = 10))))#LOOKS LIKE I NEED TO MANUALLY ADJUST AXIS TO MATCH VISUALIZATION...
  par(cex.lab=1,cex.axis=min(1,10/ncol(visualize_plot)))
  axis(1,at=1:ncol(visualize_plot), labels = colnames(visualize_plot),srt=10,tick=FALSE)
  
  box(lwd=2)
  abline(v=1:ncol(visualize_plot)-0.5,lwd=1)
  abline(h=floor(nrow(visualize_plot)/2),lwd=1,lty=1)
  
  #title(main =title, font.main = 8)
  
  
  ####################
  #3 CALC aREA ENRICHMENT
  ####################
  #FIGURE OUT WHICH LINES HAVE >1 drug that they are sensitive too...
  num_sel_drugs_per_line<-sort(colSums(gene_set_matrix),decreasing=TRUE)
  regulon_sufficient_lines<-names(which(num_sel_drugs_per_line>1))
  
  
  #MAKE REGULON: for ONLY cell-lines w/ >1 drug
  library(viper)
  gene_set_regulon<-list()
  for (column in regulon_sufficient_lines){
    column_gene_set<-names(which(gene_set_matrix[,column]==1))
    num_genes_wn_set<-length(column_gene_set)
    
    gene_set_regulon[[column]]<-list()
    gene_set_regulon[[column]]$tfmode<-rep(1,num_genes_wn_set)
    gene_set_regulon[[column]]$likelihood<-rep(1,num_genes_wn_set)
    names(gene_set_regulon[[column]]$tfmode)<-column_gene_set
  }
  
  #gene_set_regulon$AFATINIB
  #$tfmode
  #EGFR MAPKAPK2 
  #1        1 
  #$likelihood
  #[1] 1 1
  
  #RUN aREA (use ONLY >1 drug-line-NETWORK and >1 drug-line-DATA)
  aREA_all<-aREA(gene_all_data_matrix[,regulon_sufficient_lines],gene_set_regulon,minsize=2)
  
  #extract aREA on self-predictions:
  #View(aREA_all$nes)
  aREA_all_diag<-c()
  for (i in 1:nrow(aREA_all$nes)){
    aREA_all_diag<-append(aREA_all_diag,aREA_all$nes[i,i])
  }
  names(aREA_all_diag)<-rownames(aREA_all$nes)
  
  #sort aREA:
  sorted_aREA_vector<-sort(aREA_all_diag,decreasing=TRUE)
  sorted_aREA_names<-names(sorted_aREA_vector)
  
  #stouffer integration:
  integrated_p<-sum(aREA_all_diag)/sqrt(length(aREA_all_diag))
  
  sorted_visualize_plot<-visualize_plot[,sorted_aREA_names]
  
  ##################
  #PLOT aREA vs RANDOM-null:
  ##################
  ##NULL distribution: (enrichment for random combinations)
  hist(aREA_all$nes, col=rgb(0,0,0,0.2),main=paste0("NES for true Kinase Targets w/n Predictions (stouffer nes: ",as.character(round(integrated_p,3)),")"),xlab="NES(aREA)",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,freq = FALSE,breaks=100,ylim=c(0,0.5))
  #paired distribution: (top achilles wn VIPER)
  hist(sorted_aREA_vector, col=rgb(1,0,0,0.5),freq = FALSE, add=T,breaks=100,ylim=c(0,0.5))
  box()
  legend("topright",legend=c("NULL","Drug-Predictions"),fill=c(rgb(0,0,0,0.2), rgb(1,0,0,0.5)),cex=1.5)
  
  
  
  
  
  ###################
  #3 PLOT ENRICMENT:  sorted by aREA
  ###################
  image(x=1:ncol(sorted_visualize_plot),y=1:nrow(sorted_visualize_plot),z=t(sorted_visualize_plot),axes=FALSE,xlab="kinase inhibitors",ylab="ranked Kinases",col=c("white","red"))
  par(cex.lab=1,cex.axis=1)
  #axis(2,at=seq(from=93, to=3, by = -10), labels = TRUE)
  axis(2)#true rank inverted here...
  par(cex.lab=1,cex.axis=min(1,10/ncol(sorted_visualize_plot)))
  axis(1,at=1:ncol(sorted_visualize_plot), labels = colnames(sorted_visualize_plot),srt=10,tick=FALSE)
  
  box(lwd=2)
  abline(v=1:ncol(sorted_visualize_plot)-0.5,lwd=1)
  abline(h=floor(nrow(sorted_visualize_plot)/2),lwd=1,lty=1)
  
  title(main =paste0(title_main," (stouffer nes: ",as.character(round(integrated_p,3)),")"), font.main = 8)
  
  return(sorted_aREA_vector)
}



#############################################################################################################
#0 PROCESSING RAW LT-SCANNER DATA:
#############################################################################################################
{

##################################################
#UNIPROT to HUGO converter:
##################################################
hugo_uniprot_df<-read.table("hugo_uniprot.txt",sep="\t",header=T,as.is=T)
na_idx<-which(hugo_uniprot_df[,2]=="")

unipro2hugo<-hugo_uniprot_df[-na_idx,1]
names(unipro2hugo)<-hugo_uniprot_df[-na_idx,2]




##################################################
#PARSE LT-Scanner Data:
##################################################
#IMPORT RAW DATA:
raw_LTscanner<-read.csv("PrePCI_predictions_49of94_drugs.csv",as.is=T,header=T)

#APPEND HUGO NAMES TO RAW LT-SCANNER:
raw_LTscanner_hugo<-cbind(raw_LTscanner[,],unipro2hugo[raw_LTscanner[,3]],stringsAsFactors=F)

na_indx2<-which(is.na(raw_LTscanner_hugo[,5])==TRUE)
raw_LTscanner_hugo_na<-raw_LTscanner_hugo[-na_indx2,]

#MAKE MATRIX:
gene_names<-unique(raw_LTscanner_hugo_na[,5])
cmpd_names<-unique(raw_LTscanner_hugo_na[,1])

LTscanner_matrix<-matrix(nrow=length(cmpd_names),ncol=length(gene_names))
rownames(LTscanner_matrix)<-cmpd_names
colnames(LTscanner_matrix)<-gene_names

#FILL MATRIX
for (indx in 1:nrow(raw_LTscanner_hugo_na)){
  gene_name<-raw_LTscanner_hugo_na[indx,5]
  cmpd_name<-raw_LTscanner_hugo_na[indx,1]
    
  LTscanner_matrix[cmpd_name,gene_name]<-raw_LTscanner_hugo_na[indx,4] 
}


rowSums(LTscanner_matrix)


##################################################
#PROCESS LT-Scanner Data:  NA's, kinase overlap
##################################################
LTscanner_matrix_na<-LTscanner_matrix
LTscanner_matrix_na[is.na(LTscanner_matrix)]<-0

#195/216 kinases overlap
overlapping_kinases<-colnames(LTscanner_matrix_na)[which(colnames(LTscanner_matrix_na) %in% rownames(kinase_inhibitors_BINDING2))]

LTscanner_matrix_na_kinases<-LTscanner_matrix_na[,overlapping_kinases]
  
  
save(list=c("LTscanner_matrix","LTscanner_matrix_na","LTscanner_matrix_na_kinases"),file="LTscanner_parsed-data.RData")  

}



#############################################################################################################
#1 COMPARE LT-Scanner to LASSO MODEL:
#############################################################################################################
overlapping_drugs<-rownames(LTscanner_matrix_na_kinases)[which(rownames(LTscanner_matrix_na_kinases) %in% colnames(kinase_inhibitors_BINDING2))]
overlapping_kinases<-colnames(LTscanner_matrix_na_kinases)[which(colnames(LTscanner_matrix_na_kinases) %in% rownames(kinase_inhibitors_BINDING2))]


#UPPER LIMIT:
UPPERLIMIT_aREA<-GSEA_kinase_predictions(kinase_inhibitors_BINDING2[overlapping_kinases,overlapping_drugs],t(kinase_inhibitors_BINDING2[overlapping_kinases,overlapping_drugs]),threshold=6)

#LASSO:
LASSO_aREA<-GSEA_kinase_predictions(kinase_inhibitors_BINDING2[overlapping_kinases,overlapping_drugs],kinase_prediction_gesLASSO[overlapping_drugs,overlapping_kinases],threshold=6)


#LT-Scanner:
LTscanner_aREA<-GSEA_kinase_predictions(kinase_inhibitors_BINDING2[overlapping_kinases,overlapping_drugs],
                                        LTscanner_matrix_na_kinases[overlapping_drugs,overlapping_kinases],threshold=6)

plot(LASSO_aREA,LTscanner_aREA[names(LASSO_aREA)],cex.axis=1,cex.lab=1.5,
     col="white",xlab="mRNA LASSO Model Performance (z-score)",ylab="LT-scanner Model Performance (z-score)")
abline(h=0,lty=2)
abline(v=0,lty=2)
text(LASSO_aREA,LTscanner_aREA[names(LASSO_aREA)],labels=names(LASSO_aREA),cex=1.2)

#############################################################################################################
#2 KINASE RANKING
#############################################################################################################
#UPPER LIMIT:
UPPERLIMIT_aREA_kinases<-GSEA_kinase_predictions(t(kinase_inhibitors_BINDING2[overlapping_kinases,overlapping_drugs]),
                                                 kinase_inhibitors_BINDING2[overlapping_kinases,overlapping_drugs],threshold=6)

#LASSO:   doesn't work because this is a LOOCV procedure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##########THEREFORE:   each drug sees a different model and so rankings can only be assessed w/n drug column
##########                                                                                   NOT w/n a kinase column
#OK SO POORLY STAT-SIGNIFICANT, overall, because Kd' are apples to oranges between inhibitors (only apples to apples w/n 1 inhibitor)
#THAT SAID:   still have predicting binding constants and can use relative ranking to get an idea of performance per kinases
LASSO_aREA_kinases<-GSEA_kinase_predictions(t(kinase_inhibitors_BINDING2[overlapping_kinases,overlapping_drugs]),
                                            t(kinase_prediction_gesLASSO[overlapping_drugs,overlapping_kinases]),threshold=6)




#LT-Scanner:
LTscanner_aREA_kinases<-GSEA_kinase_predictions(t(kinase_inhibitors_BINDING2[overlapping_kinases,overlapping_drugs]),
                                        t(LTscanner_matrix_na_kinases[overlapping_drugs,overlapping_kinases]),threshold=6)



plot(LASSO_aREA_kinases,LTscanner_aREA_kinases[names(LASSO_aREA_kinases)],cex.axis=1,cex.lab=1.5,
     col="white",xlab="mRNA LASSO Model Performance (z-score)",ylab="LT-scanner Model Performance (z-score)")
abline(h=0,lty=2)
abline(v=0,lty=2)
text(LASSO_aREA_kinases,LTscanner_aREA_kinases[names(LASSO_aREA_kinases)],labels=names(LASSO_aREA_kinases),cex=1.2)



write.csv(cbind('mRNA-lasso'=LASSO_aREA_kinases,"LTscanner"=LTscanner_aREA_kinases[names(LASSO_aREA_kinases)]),file="KINASES_LTscanner_rnaLasso_BENCHMARK-by-Kds.csv")


write.csv(cbind('mRNA-lasso'=LASSO_aREA,"LTscanner"=LTscanner_aREA[names(LASSO_aREA)]),file="DRUGS_LTscanner_rnaLasso_BENCHMARK-by-Kds.csv")










####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########      ########      ########     ########    #########    ####     ##  ##########                  #########        ##      ##########       ##                                    
##     ##     ##     ##     ##     ##          ##           ##           ## ##    ##      ##                      ##      ##      ####         ##          ####                              
##     ##     ##     ##    ##      ##          ##           ##           ##  ##   ##      ##                      ##       ##    ##  ##        ##         ##  ##                             
##     ########      #######       ########     #######     ########     ##   ##  ##      ##                      ##       ##   ##    ##       ##        ##    ##                                
##     ##            ##    ##      ##                 ##    ##           ##    ## ##      ##                      ##       ##  ##########      ##       ##########                                                                                             
##     ##            ##     ##     ##                 ##    ##           ##     ####      ##                      ##      ##   ##      ##      ##       ##      ##                           
##     ##            ##      ##    ########     #######     #########    ##      ###      ##                      #########    ##      ##      ##       ##      ##                                      
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#DATE: 4-1-2020
#START MAKING FIGURES FOR PAPER:
#  (1)  INTRO FIGURES TO DATA-SET
#  (2)  ANALYSIS of results


#setwd("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA")
setwd("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/")



load("DREAM_kinaseInhibitor_FINAL-DATA.RData")
###load("DREAM_kinaseInhibitor_ALL-DATA.RData")
rm(FINAL_DREAM_Dose_Responses)
rm(FINAL_DREAM_RNAseq)

load("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-07-15 - REDO PANGEA ANALYSIS/HALLMARKS_drugAVERAGE.RData")
rm(PANGEA_ges_HALLMARKS_aREA)
rm(PANGEA_ges_HALLMARKS_targetAVG)
rm(PANGEA_vpr_HALLMARKS_drugAVG)
rm(PANGEA_vpr_HALLMARKS_target)
PANACEA_ges_drugAvg_Hallmarks<-PANGEA_ges_HALLMARKS_drugAVG$nes
rm(PANGEA_ges_HALLMARKS_drugAVG)


colnames(PANACEA_ges_drugAvg_Hallmarks)<-gsub("-","",gsub(" ","",gsub("_","",toupper(colnames(PANACEA_ges_drugAvg_Hallmarks)))))

OVERLAP_ges_drugAvg_Hallmarks<-PANACEA_ges_drugAvg_Hallmarks[,rownames(FINAL_DREAM_kinome_Kds)]
##################################################################################################################################
#PLOTTING FUNCTIONS:
##################################################################################################################################

library("gplots")
library("RColorBrewer")
library("devtools")

heatmap.3 <- function(x,Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,distfun = dist,hclustfun = hclust,dendrogram = c("both","row", "column", "none"),symm = FALSE,scale = c("none","row", "column"),na.rm = TRUE,revC = identical(Colv,"Rowv"),add.expr,
                      breaks,symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",col = "heat.colors",colsep,rowsep,sepcolor = "white",sepwidth = c(0.05, 0.05),cellnote,notecex = 1,notecol = "cyan",na.color = par("bg"),
                      trace = c("none", "column","row", "both"),tracecol = "cyan",hline = median(breaks),vline = median(breaks),linecol = tracecol,margins = c(5,5),ColSideColors,RowSideColors,side.height.fraction=0.3,cexRow = 0.2 + 1/log10(nr),cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,labCol = NULL,key = TRUE,keysize = 1.5,density.info = c("none", "histogram", "density"),denscol = tracecol,symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,main = NULL,xlab = NULL,ylab = NULL,lmat = NULL,lhei = NULL,lwid = NULL,ColSideColorsSize = 1,RowSideColorsSize = 1,KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}



hcluster_corr_n1screen_DrugTargets_hallmarks<-function(matrix,database,target_vector,target_colors,main_title,num_var=30){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix_idx<-which(rownames(matrix) %in% rownames(database))
  
  matrix<-matrix[matrix_idx,]
  
  
  
  #FILTER MATRIX BY ONLY DRUGS ANNOTATED: (reduce noise of drugs that do nothing...)
  master_indx_filter<-c()
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #2 extract class-indices in rlab
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(rownames(matrix),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 make master index to filter matrix
    master_indx_filter<-c(master_indx_filter,rlab_class_idx)
  }
  master_indx_filter_unique<-unique(master_indx_filter)
  
  matrix<-matrix[master_indx_filter_unique,]
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(target_vector),ncol=length(rownames(matrix)))
  rownames(rlab)<-target_vector
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #1 extract drugs-vector from n1screen_drugDB
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(colnames(rlab),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-target_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  top_var<-names(sort(apply(matrix,2,var),decreasing=T))[1:num_var]
  matrix_topvar<-matrix[,top_var]
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  cor_columns <- cor(matrix_topvar,method="spearman")
  cor_rows <- cor(t(matrix_topvar),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(length(target_vector)>1){
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(20,12),
              RowSideColors=rlab, symbreaks=TRUE, key=TRUE, symkey=TRUE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1,col=colorRampPalette(c("blue","blue", "white", "red", "red"))(n = 20),
              RowSideColorsSize=5, KeyValueName="NES",keysize = 0.8,xlab="Hallmarks",ylab="Drug-Averages (across cell-lines)")
  }else{
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(20,12),
              symbreaks=TRUE, key=TRUE, symkey=TRUE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1,col=colorRampPalette(c("blue","blue", "white", "red", "red"))(n = 20),
              KeyValueName="NES",keysize = 0.8,xlab="Hallmarks",ylab="Drug-Averages (across cell-lines)")
  }
  
  return(rownames(matrix_topvar[hclust_rows$order,]))
}



hcluster_euclid_n1screen_DrugTargets_kinasesDREAM<-function(matrix,database,target_vector,target_colors,main_title,num_var=30){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix_idx<-which(rownames(matrix) %in% rownames(database))
  
  matrix<-matrix[matrix_idx,]
  
  
  
  #FILTER MATRIX BY ONLY DRUGS ANNOTATED: (reduce noise of drugs that do nothing...)
  master_indx_filter<-c()
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #2 extract class-indices in rlab
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(rownames(matrix),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 make master index to filter matrix
    master_indx_filter<-c(master_indx_filter,rlab_class_idx)
  }
  master_indx_filter_unique<-unique(master_indx_filter)
  
  matrix<-matrix[master_indx_filter_unique,]
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(target_vector),ncol=length(rownames(matrix)))
  rownames(rlab)<-target_vector
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #1 extract drugs-vector from n1screen_drugDB
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(colnames(rlab),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-target_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  top_var<-names(sort(apply(matrix,2,var),decreasing=T))[1:num_var]
  matrix_topvar<-matrix[,top_var]
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  
  dist_columns <- dist(t(matrix_topvar))
  dist_rows <- dist(matrix_topvar)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(length(target_vector)>1){
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(12,12),
              RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1.2,col=colorRampPalette(c("blue","white", "red"))(n = 20),
              RowSideColorsSize=5, KeyValueName="-log10(Kd)",keysize = 0.8,xlab="Kinases",ylab=paste0(as.character(nrow(matrix_topvar))," Drugs"))
  }else{
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(12,12),
              symbreaks=FALSE, key=TRUE, symkey=FALSE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1.2,col=colorRampPalette(c("blue","white", "red"))(n = 20),
              KeyValueName="-log10(Kd)",keysize = 0.8,xlab="Kinases",ylab=paste0(as.character(nrow(matrix_topvar))," Drugs"))
  }
  
  return(rownames(matrix_topvar[hclust_rows$order,]))
}



hcluster_corr_n1screen_DrugTargets_hallmarks_kinomeCLUSTER<-function(matrix,database,target_vector,target_colors,main_title,num_var=30,
                                                                     rownames_order){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix_idx<-which(rownames(matrix) %in% rownames(database))
  
  matrix<-matrix[matrix_idx,]
  
  
  
  #FILTER MATRIX BY ONLY DRUGS ANNOTATED: (reduce noise of drugs that do nothing...)

  
  matrix<-matrix[kinome_row_order,]
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(target_vector),ncol=length(rownames(matrix)))
  rownames(rlab)<-target_vector
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #1 extract drugs-vector from n1screen_drugDB
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(colnames(rlab),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-target_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  top_var<-names(sort(apply(matrix,2,var),decreasing=T))[1:num_var]
  matrix_topvar<-matrix[,top_var]
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  cor_columns <- cor(matrix_topvar,method="spearman")
  cor_rows <- cor(t(matrix_topvar),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(length(target_vector)>1){
    heatmap.3(matrix_topvar, 
              Rowv=NA,
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="column", margins=c(20,12),
              RowSideColors=rlab, symbreaks=TRUE, key=TRUE, symkey=TRUE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1,col=colorRampPalette(c("blue","blue", "white", "red", "red"))(n = 20),
              RowSideColorsSize=5, KeyValueName="NES",keysize = 0.8,xlab="Hallmarks",ylab="Drug-Averages (across cell-lines)")
  }else{
    heatmap.3(matrix_topvar, 
              Rowv=NA,
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="column", margins=c(20,12),
              symbreaks=TRUE, key=TRUE, symkey=TRUE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1,col=colorRampPalette(c("blue","blue", "white", "red", "red"))(n = 20),
              KeyValueName="NES",keysize = 0.8,xlab="Hallmarks",ylab="Drug-Averages (across cell-lines)")
  }
  
  
}

hcluster_euclid_n1screen_DrugTargets_kinasesDREAM_panaceaCLUSTER<-function(matrix,database,target_vector,target_colors,main_title,num_var=30,
                                                                           rownames_order){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix_idx<-which(rownames(matrix) %in% rownames(database))
  
  matrix<-matrix[matrix_idx,]
  
  
  
  #FILTER MATRIX BY ONLY DRUGS ANNOTATED: (reduce noise of drugs that do nothing...)

  
  matrix<-matrix[rownames_order,]
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(target_vector),ncol=length(rownames(matrix)))
  rownames(rlab)<-target_vector
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #1 extract drugs-vector from n1screen_drugDB
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(colnames(rlab),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-target_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  top_var<-names(sort(apply(matrix,2,var),decreasing=T))[1:num_var]
  matrix_topvar<-matrix[,top_var]
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  
  dist_columns <- dist(t(matrix_topvar))
  dist_rows <- dist(matrix_topvar)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(length(target_vector)>1){
    heatmap.3(matrix_topvar, 
              Rowv=NA,
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="column", margins=c(12,12),
              RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1.2,col=colorRampPalette(c("blue","white", "red"))(n = 20),
              RowSideColorsSize=5, KeyValueName="-log10(Kd)",keysize = 0.8,xlab="Kinases",ylab=paste0(as.character(nrow(matrix_topvar))," Drugs"))
  }else{
    heatmap.3(matrix_topvar, 
              Rowv=NA,
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="column", margins=c(12,12),
              symbreaks=FALSE, key=TRUE, symkey=FALSE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1.2,col=colorRampPalette(c("blue","white", "red"))(n = 20),
              KeyValueName="-log10(Kd)",keysize = 0.8,xlab="Kinases",ylab=paste0(as.character(nrow(matrix_topvar))," Drugs"))
  }
  
  return(rownames(matrix_topvar[hclust_rows$order,]))
}



###################################################################################################################################################
#1 INTRODUCTION FIGURES
###################################################################################################################################################

#MIX OF INHIBITORS:  some dissagreement with Drug-Banks
kinome_row_order<-hcluster_euclid_n1screen_DrugTargets_kinasesDREAM(matrix=FINAL_DREAM_kinome_Kds,
                                                  database=FINAL_DREAM_drugLibrary_annotation,
                                                  target_vector = c("EGFR","FLT3","ABL1","RET","MET","AKT2"),
                                                  target_colors = c("darkred","darkorange","steelblue","blue","darkgreen","red"),
                                                  main_title = "KINOME-BINDING DATABASE")



panacea_row_order<-hcluster_corr_n1screen_DrugTargets_hallmarks(matrix=t(OVERLAP_ges_drugAvg_Hallmarks),
                                                  database=FINAL_DREAM_drugLibrary_annotation,
                                                  target_vector = c("EGFR","FLT3","ABL1","RET","MET","AKT2"),
                                                  target_colors = c("darkred","darkorange","steelblue","blue","darkgreen","red"),
                                                  main_title = "Differential-mRNA Database")


#############################
#RE-ORDER CLUSTERS TO MATCH OTHER-DATA-SET DRUG-order
##############################
hcluster_euclid_n1screen_DrugTargets_kinasesDREAM_panaceaCLUSTER(matrix=FINAL_DREAM_kinome_Kds,
                                                                    database=FINAL_DREAM_drugLibrary_annotation,
                                                                    target_vector = c("EGFR","FLT3","ABL1","RET","MET","AKT2"),
                                                                    target_colors = c("darkred","darkorange","steelblue","blue","darkgreen","red"),
                                                                    main_title = "KINOME-BINDING DATABASE",
                                                                    rownames_order=panacea_row_order[length(panacea_row_order):1])

hcluster_corr_n1screen_DrugTargets_hallmarks_kinomeCLUSTER(matrix=t(OVERLAP_ges_drugAvg_Hallmarks),
                                             database=FINAL_DREAM_drugLibrary_annotation,
                                             target_vector = c("EGFR","FLT3","ABL1","RET","MET","AKT2"),
                                             target_colors = c("darkred","darkorange","steelblue","blue","darkgreen","red"),
                                             main_title = "Differential-mRNA Database",
                                             rownames_order=kinome_row_order)




###################################################################################################################################################
#2 PLOTTING LEADER-BOARD SCORES:
###################################################################################################################################################
LB_SC1<-read.csv("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/raw_scores/LB_sc1.csv",header=F)$"V1"
LB_SC2<-read.csv("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/raw_scores/LB_sc2.csv",header=F)$"V1"


hist(LB_SC2[1:25],breaks=20,col="yellowgreen",xlim=c(0,60),ylim=c(0,35),cex.main=1.5,
     main="Subchallenge 2 (SC2) Leaderboard Results",xlab="-log10(p-value)",cex=1.5,cex.axis=1.5,cex.lab=1.5)
hist(LB_SC2[26:53],breaks=6,col="gold2",add=T)
hist(LB_SC2[54:86],breaks=1,col="firebrick",add=T)





hist(LB_SC1[1:14],breaks=6,col="yellowgreen",xlim=c(0,17),ylim=c(0,10),cex.main=1.5,
     main="Subchallenge 1 (SC1) Leaderboard Results",xlab="-log10(p-value)",cex=1.5,cex.axis=1.5,cex.lab=1.5)
hist(LB_SC1[14:39],breaks=20,col="gold2",add=T)
hist(LB_SC1[40:86],breaks=3,col="firebrick",add=T)




























































####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################  
##     ########   ##   ###     ##      ##        ##                ########   ##      ##    ########       #######                                         
##     ##         ##   ####    ##     ####       ##               ##          ##      ##    ##     ##     ##
##     ##         ##   ## ##   ##    ##  ##      ##               ##          ##      ##    ##     ##     ##
##     ########   ##   ##  ##  ##   ##    ##     ##                #######    ##      ##    ########       #######                                      
##     ##         ##   ##   ## ##  ##########    ##                      ##   ##      ##    ##     ##            ##
##     ##         ##   ##    ####  ##      ##    ##                      ##    ##    ##     ##     ##            ##
##     ##         ##   ##     ###  ##      ##    ########          #######      ######      ########       #######          
####################################################################################################################################################################################################################################################################  
####################################################################################################################################################################################################################################################################  


setwd("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/")


#GOLD STANDARD:
load("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-05-24 - N1screen - drug-clean-ness/OVERLAP_n1-Kinome-Screens.RData")
rm(kinome_selectivity_scores_OVERLAP)


#LOAD PARSED FINAL SUBMISSIONS:  #0 below
load("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/final_results_processed.RData")


#LOAD LOOCV LASSO-model
library("viper")
library("gplots")
library("RColorBrewer", lib.loc="~/Library/R/3.3/library")
library("devtools")

GSEA_kinase_predictions_FINAL_plot<-function(kinase_inhibitors_BINDING2,kinase_prediction_matrix_na,threshold=6,title_main="True-Targets w/n Predictions"){
  
  #############################################################################################
  #1 FORMAT DATA SO CAN DO aREA ENRICHEMENT EASILY:
  #############################################################################################
  
  
  ########################
  #BINARY-IZE AUC RESPONSE: i.e. define "GOLD STANDARDS"
  ########################
  
  stat_sig_drug_idx<-which(kinase_inhibitors_BINDING2 > threshold)
  
  kinase_inhibitors_BINDING2_binary<-kinase_inhibitors_BINDING2
  kinase_inhibitors_BINDING2_binary[stat_sig_drug_idx]<-1
  kinase_inhibitors_BINDING2_binary[-stat_sig_drug_idx]<-0
  
  
  percent_sensitive_drugs=length(stat_sig_drug_idx)/(nrow(kinase_inhibitors_BINDING2_binary)*ncol(kinase_inhibitors_BINDING2_binary))
  #4.7%
  
  
  
  
  ########################
  #MAKE MATRIX OF RANKED DRUGS / CELL-LINE:  rank  x cell-line
  ########################
  kinase_PREDICTIONS_formatted<-t(kinase_prediction_matrix_na)
  
  
  
  
  ##############
  #NOTE: more positive score = better binding (i.e. predictions convert VIPER into -log10(pKd))
  ##############
  #sort(kinase_PREDICTIONS_formatted[,"ERLOTINIB"],decreasing=T)[1:10]
  #MAP2K2   MAP2K1     PFKP     BRAF      NLK    TESK1      KIT    IGF1R     INSR     ARAF 
  #260.4869 231.2391 194.9972 166.3905 143.4486 143.4486 143.4486 129.8185 129.8185 125.8992 
  
  #MAKE RANKED MATRIX:
  num_kinases<-length(rownames(kinase_PREDICTIONS_formatted))
  
  kinaseRanks_per_drug<-matrix(nrow=num_kinases,ncol=0)
  
  rownames(kinaseRanks_per_drug)<-1:num_kinases
  
  col_names<-c()
  for (column in colnames(kinase_PREDICTIONS_formatted)){
    #1 store cell-line name:
    col_names<-append(col_names,column)
    #2 organize kinases by rank
    temp_ranked_kinases<-names(sort(kinase_PREDICTIONS_formatted[,column],decreasing=TRUE))
    #3 append ranked kinase-list matrix:
    kinaseRanks_per_drug<-cbind(kinaseRanks_per_drug,temp_ranked_kinases)
    #4 name columns:
    colnames(kinaseRanks_per_drug)<-col_names
  }
  
  
  
  
  ######################################################################################################################################################################################
  #PLOTTING AND ANALYSIS   (original function)
  ######################################################################################################################################################################################
  gene_set_matrix=kinase_inhibitors_BINDING2_binary#BINARY GOLD STANDARD DEFINITIONS = define regulon
  gene_all_matrix=kinaseRanks_per_drug #KINASE-RANKING/LINE = for enrichment visualization
  gene_all_data_matrix=kinase_PREDICTIONS_formatted#KINASE-SCORE/DRUG 
  
  
  
  ###################
  #1 MAKE PLOT
  ###################
  blank_enrichment_plot_matrix<-matrix(1,ncol=ncol(gene_all_matrix),nrow=nrow(gene_all_matrix))
  colnames(blank_enrichment_plot_matrix)<-colnames(gene_all_matrix)
  
  #column="AFATINIB"
  for (column in colnames(blank_enrichment_plot_matrix)){
    #1 extract gene-set w/ 1 in matrix:
    column_gene_set<-names(which(gene_set_matrix[,column]==1))
    #[1] "EGFR"     "MAPKAPK2"
    
    #3 make plottable matrix:
    gene_set_idx_wn_plot_column<-which(gene_all_matrix[,column] %in% column_gene_set)
    blank_enrichment_plot_matrix[gene_set_idx_wn_plot_column,column]<-2
  }
  
  visualize_plot <- blank_enrichment_plot_matrix[ nrow(blank_enrichment_plot_matrix):1, ]
  

  
  
  ####################
  #3 CALC aREA ENRICHMENT
  ####################
  #FIGURE OUT WHICH LINES HAVE >1 drug that they are sensitive too...
  num_sel_drugs_per_line<-sort(colSums(gene_set_matrix),decreasing=TRUE)
  regulon_sufficient_lines<-names(which(num_sel_drugs_per_line>0))
  
  
  #MAKE REGULON: for ONLY cell-lines w/ >1 drug
  library(viper)
  gene_set_regulon<-list()
  for (column in regulon_sufficient_lines){
    column_gene_set<-names(which(gene_set_matrix[,column]==1))
    num_genes_wn_set<-length(column_gene_set)
    
    gene_set_regulon[[column]]<-list()
    gene_set_regulon[[column]]$tfmode<-rep(1,num_genes_wn_set)
    gene_set_regulon[[column]]$likelihood<-rep(1,num_genes_wn_set)
    names(gene_set_regulon[[column]]$tfmode)<-column_gene_set
  }
  
  #gene_set_regulon$AFATINIB
  #$tfmode
  #EGFR MAPKAPK2 
  #1        1 
  #$likelihood
  #[1] 1 1
  
  #RUN aREA (use ONLY >1 drug-line-NETWORK and >1 drug-line-DATA)
  aREA_all<-aREA(gene_all_data_matrix[,regulon_sufficient_lines],gene_set_regulon,minsize=1)
  
  #extract aREA on self-predictions:
  #View(aREA_all$nes)
  aREA_all_diag<-c()
  for (i in 1:nrow(aREA_all$nes)){
    aREA_all_diag<-append(aREA_all_diag,aREA_all$nes[i,i])
  }
  names(aREA_all_diag)<-rownames(aREA_all$nes)
  
  #sort aREA:
  sorted_aREA_vector<-sort(aREA_all_diag,decreasing=TRUE)
  sorted_aREA_names<-names(sorted_aREA_vector)
  
  #stouffer integration:
  integrated_p<-sum(aREA_all_diag)/sqrt(length(aREA_all_diag))
  
  sorted_visualize_plot<-visualize_plot[,sorted_aREA_names]
  

  
  
  
  ###################
  #3 PLOT ENRICMENT:  sorted by aREA
  ###################
  par(mar=c(10, 4, 4, 2))
  image(x=1:ncol(sorted_visualize_plot),y=1:nrow(sorted_visualize_plot),z=t(sorted_visualize_plot),axes=FALSE,xlab=NA,ylab="ranked Kinases",col=c("white","red"))
  par(cex.lab=1,cex.axis=1,mar=c(10, 4, 4, 2))
  #axis(2,at=seq(from=93, to=3, by = -10), labels = TRUE)
  axis(2)#true rank inverted here...
  par(cex.lab=1,cex.axis=1,las=3,mar=c(10, 4, 4, 2))
  axis(1,at=1:ncol(sorted_visualize_plot), labels = paste0(names(sorted_aREA_vector)," (",as.character(round(sorted_aREA_vector,1)),")"),srt=10,tick=FALSE)
  
  box(lwd=2)
  abline(v=1:ncol(sorted_visualize_plot)-0.5,lwd=1)
  abline(h=floor(nrow(sorted_visualize_plot)/2),lwd=1,lty=1)
  
  title(main =paste0(title_main," (stouffer: ",as.character(round(integrated_p,1)),")"), font.main = 8,cex.main=1.7)
  
  return(sorted_aREA_vector)
}


GSEA_kinase_predictions_FINAL_plot_kinases<-function(kinase_inhibitors_BINDING2,kinase_prediction_matrix_na,threshold=6,title_main="True-Targets w/n Predictions"){
  
  #############################################################################################
  #1 FORMAT DATA SO CAN DO aREA ENRICHEMENT EASILY:
  #############################################################################################
  
  
  ########################
  #BINARY-IZE AUC RESPONSE: i.e. define "GOLD STANDARDS"
  ########################
  
  stat_sig_drug_idx<-which(kinase_inhibitors_BINDING2 > threshold)
  
  kinase_inhibitors_BINDING2_binary<-kinase_inhibitors_BINDING2
  kinase_inhibitors_BINDING2_binary[stat_sig_drug_idx]<-1
  kinase_inhibitors_BINDING2_binary[-stat_sig_drug_idx]<-0
  
  
  percent_sensitive_drugs=length(stat_sig_drug_idx)/(nrow(kinase_inhibitors_BINDING2_binary)*ncol(kinase_inhibitors_BINDING2_binary))
  #4.7%
  
  
  
  
  ########################
  #MAKE MATRIX OF RANKED DRUGS / CELL-LINE:  rank  x cell-line
  ########################
  kinase_PREDICTIONS_formatted<-t(kinase_prediction_matrix_na)
  
  
  
  
  ##############
  #NOTE: more positive score = better binding (i.e. predictions convert VIPER into -log10(pKd))
  ##############
  #sort(kinase_PREDICTIONS_formatted[,"ERLOTINIB"],decreasing=T)[1:10]
  #MAP2K2   MAP2K1     PFKP     BRAF      NLK    TESK1      KIT    IGF1R     INSR     ARAF 
  #260.4869 231.2391 194.9972 166.3905 143.4486 143.4486 143.4486 129.8185 129.8185 125.8992 
  
  #MAKE RANKED MATRIX:
  num_kinases<-length(rownames(kinase_PREDICTIONS_formatted))
  
  kinaseRanks_per_drug<-matrix(nrow=num_kinases,ncol=0)
  
  rownames(kinaseRanks_per_drug)<-1:num_kinases
  
  col_names<-c()
  for (column in colnames(kinase_PREDICTIONS_formatted)){
    #1 store cell-line name:
    col_names<-append(col_names,column)
    #2 organize kinases by rank
    temp_ranked_kinases<-names(sort(kinase_PREDICTIONS_formatted[,column],decreasing=TRUE))
    #3 append ranked kinase-list matrix:
    kinaseRanks_per_drug<-cbind(kinaseRanks_per_drug,temp_ranked_kinases)
    #4 name columns:
    colnames(kinaseRanks_per_drug)<-col_names
  }
  
  
  
  
  ######################################################################################################################################################################################
  #PLOTTING AND ANALYSIS   (original function)
  ######################################################################################################################################################################################
  gene_set_matrix=kinase_inhibitors_BINDING2_binary#BINARY GOLD STANDARD DEFINITIONS = define regulon
  gene_all_matrix=kinaseRanks_per_drug #KINASE-RANKING/LINE = for enrichment visualization
  gene_all_data_matrix=kinase_PREDICTIONS_formatted#KINASE-SCORE/DRUG 
  
  
  
  ###################
  #1 MAKE PLOT
  ###################
  blank_enrichment_plot_matrix<-matrix(1,ncol=ncol(gene_all_matrix),nrow=nrow(gene_all_matrix))
  colnames(blank_enrichment_plot_matrix)<-colnames(gene_all_matrix)
  
  #column="AFATINIB"
  for (column in colnames(blank_enrichment_plot_matrix)){
    #1 extract gene-set w/ 1 in matrix:
    column_gene_set<-names(which(gene_set_matrix[,column]==1))
    #[1] "EGFR"     "MAPKAPK2"
    
    #3 make plottable matrix:
    gene_set_idx_wn_plot_column<-which(gene_all_matrix[,column] %in% column_gene_set)
    blank_enrichment_plot_matrix[gene_set_idx_wn_plot_column,column]<-2
  }
  
  visualize_plot <- blank_enrichment_plot_matrix[ nrow(blank_enrichment_plot_matrix):1, ]
  
  
  
  
  ####################
  #3 CALC aREA ENRICHMENT
  ####################
  #FIGURE OUT WHICH LINES HAVE >1 drug that they are sensitive too...
  num_sel_drugs_per_line<-sort(colSums(gene_set_matrix),decreasing=TRUE)
  regulon_sufficient_lines<-names(which(num_sel_drugs_per_line>0))
  
  
  #MAKE REGULON: for ONLY cell-lines w/ >1 drug
  library(viper)
  gene_set_regulon<-list()
  for (column in regulon_sufficient_lines){
    column_gene_set<-names(which(gene_set_matrix[,column]==1))
    num_genes_wn_set<-length(column_gene_set)
    
    gene_set_regulon[[column]]<-list()
    gene_set_regulon[[column]]$tfmode<-rep(1,num_genes_wn_set)
    gene_set_regulon[[column]]$likelihood<-rep(1,num_genes_wn_set)
    names(gene_set_regulon[[column]]$tfmode)<-column_gene_set
  }
  
  #gene_set_regulon$AFATINIB
  #$tfmode
  #EGFR MAPKAPK2 
  #1        1 
  #$likelihood
  #[1] 1 1
  
  #RUN aREA (use ONLY >1 drug-line-NETWORK and >1 drug-line-DATA)
  aREA_all<-aREA(gene_all_data_matrix[,regulon_sufficient_lines],gene_set_regulon,minsize=1)
  
  #extract aREA on self-predictions:
  #View(aREA_all$nes)
  aREA_all_diag<-c()
  for (i in 1:nrow(aREA_all$nes)){
    aREA_all_diag<-append(aREA_all_diag,aREA_all$nes[i,i])
  }
  names(aREA_all_diag)<-rownames(aREA_all$nes)
  
  #sort aREA:
  sorted_aREA_vector<-sort(aREA_all_diag,decreasing=TRUE)
  sorted_aREA_names<-names(sorted_aREA_vector)
  
  #stouffer integration:
  integrated_p<-sum(aREA_all_diag)/sqrt(length(aREA_all_diag))
  
  sorted_visualize_plot<-visualize_plot[,sorted_aREA_names[1:50]]
  
  
  
  
  
  ###################
  #3 PLOT ENRICMENT:  sorted by aREA
  ###################
  par(mar=c(10, 4, 4, 2))
  image(x=1:ncol(sorted_visualize_plot),y=1:nrow(sorted_visualize_plot),z=t(sorted_visualize_plot),axes=FALSE,xlab=NA,ylab="ranked Drugs",col=c("white","red"))
  par(cex.lab=1,cex.axis=1,mar=c(10, 4, 4, 2))
  #axis(2,at=seq(from=93, to=3, by = -10), labels = TRUE)
  axis(2)#true rank inverted here...
  par(cex.lab=1,cex.axis=0.7,las=3,mar=c(10, 4, 4, 2))
  axis(1,at=1:ncol(sorted_visualize_plot), labels = paste0(names(sorted_aREA_vector[1:50])," (",as.character(round(sorted_aREA_vector[1:50],1)),")"),srt=10,tick=FALSE)
  
  box(lwd=2)
  abline(v=1:ncol(sorted_visualize_plot)-0.5,lwd=1)
  abline(h=floor(nrow(sorted_visualize_plot)/2),lwd=1,lty=1)
  
  title(main =paste0(title_main," (stouffer: ",as.character(round(integrated_p,1)),")"), font.main = 8,cex.main=1.7)
  
  return(sorted_aREA_vector)
}




FISHER_kinase_predictions_FINAL_plot<-function(kinase_inhibitors_BINDING2,kinase_prediction_matrix_na,threshold=6,title_main="True-Targets w/n Predictions"){
  
  #############################################################################################
  #1 FORMAT DATA SO CAN DO aREA ENRICHEMENT EASILY:
  #############################################################################################
  
  
  ########################
  #BINARY-IZE AUC RESPONSE: i.e. define "GOLD STANDARDS"
  ########################
  
  stat_sig_drug_idx<-which(kinase_inhibitors_BINDING2 > threshold)
  
  kinase_inhibitors_BINDING2_binary<-kinase_inhibitors_BINDING2
  kinase_inhibitors_BINDING2_binary[stat_sig_drug_idx]<-1
  kinase_inhibitors_BINDING2_binary[-stat_sig_drug_idx]<-0
  
  
  percent_sensitive_drugs=length(stat_sig_drug_idx)/(nrow(kinase_inhibitors_BINDING2_binary)*ncol(kinase_inhibitors_BINDING2_binary))
  #4.7%
  
  
  
  
  ########################
  #MAKE MATRIX OF RANKED DRUGS / CELL-LINE:  rank  x cell-line
  ########################
  kinase_PREDICTIONS_formatted<-t(kinase_prediction_matrix_na)
  
  
  
  
  ##############
  #NOTE: more positive score = better binding (i.e. predictions convert VIPER into -log10(pKd))
  ##############
  #sort(kinase_PREDICTIONS_formatted[,"ERLOTINIB"],decreasing=T)[1:10]
  #MAP2K2   MAP2K1     PFKP     BRAF      NLK    TESK1      KIT    IGF1R     INSR     ARAF 
  #260.4869 231.2391 194.9972 166.3905 143.4486 143.4486 143.4486 129.8185 129.8185 125.8992 
  
  #MAKE RANKED MATRIX:
  num_kinases<-length(rownames(kinase_PREDICTIONS_formatted))
  
  kinaseRanks_per_drug<-matrix(nrow=num_kinases,ncol=0)
  
  rownames(kinaseRanks_per_drug)<-1:num_kinases
  
  col_names<-c()
  for (column in colnames(kinase_PREDICTIONS_formatted)){
    #1 store cell-line name:
    col_names<-append(col_names,column)
    #2 organize kinases by rank
    temp_ranked_kinases<-names(sort(kinase_PREDICTIONS_formatted[,column],decreasing=TRUE))
    #3 append ranked kinase-list matrix:
    kinaseRanks_per_drug<-cbind(kinaseRanks_per_drug,temp_ranked_kinases)
    #4 name columns:
    colnames(kinaseRanks_per_drug)<-col_names
  }
  
  
  
  
  ######################################################################################################################################################################################
  #PLOTTING AND ANALYSIS   (original function)
  ######################################################################################################################################################################################
  gene_set_matrix=kinase_inhibitors_BINDING2_binary#BINARY GOLD STANDARD DEFINITIONS = define regulon
  gene_all_matrix=kinaseRanks_per_drug #KINASE-RANKING/LINE = for enrichment visualization
  gene_all_data_matrix=kinase_PREDICTIONS_formatted#KINASE-SCORE/DRUG 
  
  
  
  ###################
  #1 MAKE PLOT
  ###################
  blank_enrichment_plot_matrix<-matrix(1,ncol=ncol(gene_all_matrix),nrow=nrow(gene_all_matrix))
  colnames(blank_enrichment_plot_matrix)<-colnames(gene_all_matrix)
  
  #column="AFATINIB"
  for (column in colnames(blank_enrichment_plot_matrix)){
    #1 extract gene-set w/ 1 in matrix:
    column_gene_set<-names(which(gene_set_matrix[,column]==1))
    #[1] "EGFR"     "MAPKAPK2"
    
    #3 make plottable matrix:
    gene_set_idx_wn_plot_column<-which(gene_all_matrix[,column] %in% column_gene_set)
    blank_enrichment_plot_matrix[gene_set_idx_wn_plot_column,column]<-2
  }
  
  visualize_plot <- blank_enrichment_plot_matrix[ nrow(blank_enrichment_plot_matrix):1, ]
  
  
  
  
  ####################
  #3 CALC FISHER EXACT TEST RANKING OF COLUMNS (drugs)
  ####################
  #FIGURE OUT WHICH LINES HAVE >1 drug that they are sensitive too...
  num_sel_drugs_per_line<-sort(colSums(gene_set_matrix),decreasing=TRUE)
  regulon_sufficient_lines<-names(which(num_sel_drugs_per_line>0))
  
  aREA_vector<-c()
  #FOR EACH COLUMN CALCULATE & STORE FISHER EXACT TEST FOR HITS w/n TOP-10
  #column=regulon_sufficient_lines[1]
  for (column in regulon_sufficient_lines){
    
    
    
    sub_uM_hits<-names(which(gene_set_matrix[,column]==1))
    hits_total<-length(sub_uM_hits)
    
    top_10_predicted_genes<-names(sort(gene_all_data_matrix[,column],decreasing=T)[1:10])
    
    hits_top10<-length(which(sub_uM_hits %in% top_10_predicted_genes))
    
    hits_bottom243<-hits_total-hits_top10
    
    tmp_test<-fisher.test(rbind(c(hits_top10,hits_bottom243),
                      c((10-hits_top10),(243-hits_bottom243))))
    
    tmp_log10_negP<-(-log10(tmp_test$p.value))
    
    aREA_vector<-c(aREA_vector,tmp_log10_negP)
  }
  names(aREA_vector)<-regulon_sufficient_lines

  
  #sort aREA:
  sorted_aREA_vector<-sort(aREA_vector,decreasing=TRUE)
  sorted_aREA_names<-names(sorted_aREA_vector)
  
  #stouffer integration:
  integrated_p<-sum(sorted_aREA_vector)/length(sorted_aREA_vector)
  
  sorted_visualize_plot<-visualize_plot[,sorted_aREA_names]
  
  
  
  
  
  ###################
  #3 PLOT ENRICMENT:  sorted by aREA
  ###################
  par(mar=c(10, 4, 4, 2))
  image(x=1:ncol(sorted_visualize_plot),y=1:nrow(sorted_visualize_plot),z=t(sorted_visualize_plot),axes=FALSE,xlab=NA,ylab="ranked Kinases",col=c("white","red"))
  par(cex.lab=1,cex.axis=1,mar=c(10, 4, 4, 2))
  #axis(2,at=seq(from=93, to=3, by = -10), labels = TRUE)
  axis(2)#true rank inverted here...
  par(cex.lab=1,cex.axis=1,las=3,mar=c(10, 4, 4, 2))
  axis(1,at=1:ncol(sorted_visualize_plot), labels = paste0(names(sorted_aREA_vector)," (",as.character(round(sorted_aREA_vector,1)),")"),srt=10,tick=FALSE)
  
  box(lwd=2)
  abline(v=1:ncol(sorted_visualize_plot)-0.5,lwd=1)
  abline(h=243,lwd=1.5,lty=1)
  
  title(main =paste0(title_main," (avg -Log10(p): ",as.character(round(integrated_p,1)),")"), font.main = 8,cex.main=1.7)
  
  return(sorted_aREA_vector)
}

###################################################################################################################################################
#0 PROCESS FINAL SUBMISSIONS:   convert drug-names, slice out kinases
###################################################################################################################################################
{
kinases_final<-colnames(kinome_inhibitor_matrix_pKa_OVERLAP)
cmpd_id_key<-read.csv("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/random_cmpd_map.csv",as.is=T,header=T)

key2cmpd<-cmpd_id_key$cmpd
names(key2cmpd)<-cmpd_id_key$cmpd_id


setwd("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/07 FINAL Results/")

#MASTER LIST TO STORE SUBMISSIONS:
final_results_processed<-list()
#LIST RAW SUBMISSION FILES:
file_list<-list.files(pattern=".csv")
#file=file_list[1]
for (file in file_list){
  team_name<-strsplit(file,split=".csv")[[1]]
  
  raw_data<-data.matrix(read.csv(file,as.is=T,header=T,row.names = 1))
  
  raw_data_kinases<-raw_data[which(rownames(raw_data) %in% kinases_final),]
  
  final_data<-raw_data_kinases
  colnames(final_data)<-key2cmpd[colnames(raw_data_kinases)]
  
  
  final_results_processed[[team_name]]<-t(final_data)
}


save(final_results_processed,file="/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/final_results_processed.RData")

}


###################################################################################################################################################
#1 DRUG-WISE GSEA validation:
###################################################################################################################################################

overlapping_drugs<-rownames(final_results_processed$AMbeRland)
overlapping_kinases<-colnames(final_results_processed$AMbeRland)


#UPPER LIMIT:
UPPERLIMIT_aREA<-GSEA_kinase_predictions_FINAL_plot(t(kinome_inhibitor_matrix_pKa_OVERLAP[overlapping_drugs,overlapping_kinases]),
                                         kinome_inhibitor_matrix_pKa_OVERLAP[overlapping_drugs,overlapping_kinases])
#TEST PARTICIPANT:
netphar_test<-GSEA_kinase_predictions_FINAL_plot(t(kinome_inhibitor_matrix_pKa_OVERLAP[overlapping_drugs,overlapping_kinases]),
                                        final_results_processed$netphar,title_main = "True-Targets w/n Netphar Predictions")




finalScores_drugs<-matrix(nrow=length(overlapping_drugs),ncol=0)
rownames(finalScores_drugs)<-overlapping_drugs

#team="netphar"
for (team in names(final_results_processed)){
  #extract result vector
  tmp_result_vector<-GSEA_kinase_predictions_FINAL_plot(t(kinome_inhibitor_matrix_pKa_OVERLAP[overlapping_drugs,overlapping_kinases]),
                                     final_results_processed[[team]],title_main = paste0("True-Targets w/n ",team," Predictions"))
  
  #append to master-matrix:
  finalScores_drugs<-cbind(finalScores_drugs[overlapping_drugs,],tmp_result_vector[overlapping_drugs])
  
  jpeg(filename = paste0(team,"-drug-Scores2.jpg"), width = 800, height = 600)
  GSEA_kinase_predictions_FINAL_plot(t(kinome_inhibitor_matrix_pKa_OVERLAP[overlapping_drugs,overlapping_kinases]),
                                     final_results_processed[[team]],title_main = paste0("True-Targets w/n ",team," Predictions"))
  
  dev.off()
}

colnames(finalScores_drugs)<-names(final_results_processed)

finalScores_drugs_ceiling<-cbind(finalScores_drugs[overlapping_drugs,],"UPPER"=UPPERLIMIT_aREA[overlapping_drugs])

save(finalScores_drugs_ceiling,file="finalScores_drugs.RData")

###################################################################################################################################################
#2 DRUG-WISE FISHER-top10 validation:
###################################################################################################################################################

overlapping_drugs<-rownames(final_results_processed$AMbeRland)
overlapping_kinases<-colnames(final_results_processed$AMbeRland)


#UPPER LIMIT:
UPPERLIMIT_aREA<-FISHER_kinase_predictions_FINAL_plot(t(kinome_inhibitor_matrix_pKa_OVERLAP[overlapping_drugs,overlapping_kinases]),
                                                    kinome_inhibitor_matrix_pKa_OVERLAP[overlapping_drugs,overlapping_kinases])
#TEST PARTICIPANT:
netphar_test<-FISHER_kinase_predictions_FINAL_plot(t(kinome_inhibitor_matrix_pKa_OVERLAP[overlapping_drugs,overlapping_kinases]),
                                                 final_results_processed$netphar,title_main = "True-Targets w/n Netphar Predictions")





finalScores_drugs_top10<-matrix(nrow=length(overlapping_drugs),ncol=0)
rownames(finalScores_drugs_top10)<-overlapping_drugs

#team="netphar"
for (team in names(final_results_processed)){
  #extract result vector
  tmp_result_vector<-FISHER_kinase_predictions_FINAL_plot(t(kinome_inhibitor_matrix_pKa_OVERLAP[overlapping_drugs,overlapping_kinases]),
                                                        final_results_processed[[team]],title_main = paste0("True-Targets w/n ",team," Predictions"))
  
  #append to master-matrix:
  finalScores_drugs_top10<-cbind(finalScores_drugs_top10[overlapping_drugs,],tmp_result_vector[overlapping_drugs])
  
  jpeg(filename = paste0("SC1 00 - -",team,"-Fisher.jpg"), width = 800, height = 600)
  FISHER_kinase_predictions_FINAL_plot(t(kinome_inhibitor_matrix_pKa_OVERLAP[overlapping_drugs,overlapping_kinases]),
                                     final_results_processed[[team]],title_main = paste0("True-Targets w/n ",team," Predictions"))
  
  dev.off()
}

colnames(finalScores_drugs_top10)<-names(final_results_processed)

finalScores_drugs_top10_ceiling<-cbind(finalScores_drugs_top10[overlapping_drugs,],"UPPER"=UPPERLIMIT_aREA[overlapping_drugs])

save(finalScores_drugs_top10_ceiling,file="finalScores_drugs_FISHER.RData")














####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      ########   ##         ##     ##     ########   ##########   #########   ########     ########    ###     ##      #########                                                      
##     ##          ##         ##     ##    ##              ##       ##          ##     ##       ##       ####    ##     ## 
##    ##           ##         ##     ##    ##              ##       ##          ##     ##       ##       ## ##   ##    ##                 
##    ##           ##         ##     ##     #######        ##       #######     ########        ##       ##  ##  ##    ##    ####
##    ##           ##         ##     ##           ##       ##       ##          ##    ##        ##       ##   ## ##    ##       ##
##     ##          ##         ##     ##           ##       ##       ##          ##     ##       ##       ##    ####     ##      ##       
##      ########   ########    #######     ########        ##       #########   ##      ##   ########    ##     ###      ########     
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################

setwd("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/")

load("finalScores_drugs.RData")

load("finalScores_drugs_FISHER.RData")


library("gplots", lib.loc="~/Library/R/3.3/library")
library("RColorBrewer", lib.loc="~/Library/R/3.3/library")


hcluster<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
  #distance matrix:
  dist_columns <- dist(t(matrix))
  dist_rows <- dist(matrix)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "NES",key.xlab = "NES")
}

hcluster_corr<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
  #distance matrix:
  cor_columns <- cor(matrix,method="spearman")
  cor_rows <- cor(t(matrix),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(10,12), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "NES",key.xlab = "NES")
}


#GSEA- all predictions

hcluster_corr(t(finalScores_drugs_ceiling[,-12]),title="Drug-wise Scores per Team",col_size = 0.8,row_size=1.5)

hcluster_corr(t(apply(finalScores_drugs_ceiling[,-12],2,rank)),title="Drug-wise Scores per Team",col_size = 1.3,row_size=1.5)


#FISHER - top 10 predictions

hcluster_corr(t(finalScores_drugs_top10_ceiling[,-12]),title="Drug-wise Scores per Team",col_size = 0.8,row_size=1.5)

hcluster_corr(t(apply(finalScores_drugs_top10_ceiling[,-12],2,rank)),title="Drug-wise Scores per Team",col_size = 1.3,row_size=1.5)



####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##    ########     ########    ##     ##     ########              ########       ##       ###     ##   ##     ##                                       
##    ##     ##    ##     ##   ##     ##    ##                     ##     ##     ####      ####    ##   ##    ##                     
##    ##      ##   ##     ##   ##     ##   ##                      ##     ##    ##  ##     ## ##   ##   ##   ##                         
##    ##      ##   ########    ##     ##   ##    ####    ######    ########    ##    ##    ##  ##  ##   ######                    
##    ##      ##   ##    ##    ##     ##   ##       ##             ##     ##  ##########   ##   ## ##   ##   ##                      
##    ##     ##    ##     ##   ##     ##    ##      ##             ##     ##  ##      ##   ##    ####   ##    ##          
##    ########     ##      ##   #######      ########              ########   ##      ##   ##     ###   ##     ##                                                 
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##      ########     #########    ##########        ##              
##     ##           ##            ##               ####                                        
##    ##            ##            ##              ##  ##                              
##    ##    ####     ########     ########       ##    ##                                                                                                                                  
##    ##       ##           ##    ##            ##########                                                         
##     ##      ##           ##    ##           ##        ##                                                 
##      ########     ########     ##########  ##          ##                                                                                
####################################################################################################################################################################################################################################################################
##    ########     ########    ##     ##     ########             ##############     ##        ########      ########    ########   ##########                  ##       ##      ##     #########                                                                  
##    ##     ##    ##     ##   ##     ##    ##                          ##          ####       ##     ##    ##           ##             ##                     ####      ##      ##    ##                                              
##    ##      ##   ##     ##   ##     ##   ##                           ##         ##  ##      ##     ##   ##            ##             ##                    ##  ##     ##      ##   ##                                               
##    ##      ##   ########    ##     ##   ##    ####    ######         ##        ##    ##     ########    ##     ####   ########       ##      ######       ##    ##     ##    ##    ##     #####                                 
##    ##      ##   ##    ##    ##     ##   ##       ##                  ##       ##########    ##    ##    ##        ##  ##             ##                  ##########     ##  ##     ##         ##                                      
##    ##     ##    ##     ##   ##     ##    ##      ##                  ##      ##        ##   ##     ##    ##      ##   ##             ##                 ##        ##     ####       ##       ##                                        
##    ########     ##      ##   #######      ########                   ##     ##          ##  ##      ##    ########    ########       ##                ##          ##     ##         #########                                                                         
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################


setwd("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/")

load("/Volumes/archive_shares/n-of-1-plateseq-data/n1screen_paper/n1database_v2/n1database_org/n1druglibrary_annotation_CTEP.RData")
n1druglibrary_annotation_CAPS<-n1druglibrary_annotation_CTEP
rownames(n1druglibrary_annotation_CAPS)<-make.names(gsub("_","",gsub("-","",gsub(" ","",toupper(rownames(n1druglibrary_annotation_CAPS))))),unique=T)

load("finalScores_drugs.RData")
load("finalScores_drugs_FISHER.RData")



########################
#SC1 & SC2 tightly correlated
########################
plot(rowMeans(apply(finalScores_drugs_ceiling[,],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),],2,rank)),col="white")
text(rowMeans(apply(finalScores_drugs_ceiling,2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),],2,rank)),labels=rownames(finalScores_drugs_ceiling))

########################
#Netphar vs Atom
########################

plot(rowMeans(apply(finalScores_drugs_ceiling[,c("netphar","netphar")],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),c("Atom","Atom")],2,rank)),col="white")
text(rowMeans(apply(finalScores_drugs_ceiling[,c("netphar","netphar")],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),c("Atom","Atom")],2,rank)),labels=rownames(finalScores_drugs_ceiling))

plot(rowMeans(apply(finalScores_drugs_ceiling[,c("netphar","netphar")],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),c("SBNB","SBNB")],2,rank)),col="white")
text(rowMeans(apply(finalScores_drugs_ceiling[,c("netphar","netphar")],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),c("SBNB","SBNB")],2,rank)),labels=rownames(finalScores_drugs_ceiling))

plot(rowMeans(apply(finalScores_drugs_ceiling[,c("Atom","Atom")],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),c("SBNB","SBNB")],2,rank)),col="white")
text(rowMeans(apply(finalScores_drugs_ceiling[,c("Atom","Atom")],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),c("SBNB","SBNB")],2,rank)),labels=rownames(finalScores_drugs_ceiling))


library("viper")
library("gplots")
library("RColorBrewer", lib.loc="~/Library/R/3.3/library")
library("devtools")
heatmap.3_sidebarLabel <- function(x,drug_label_size=1,Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,distfun = dist,hclustfun = hclust,dendrogram = c("both","row", "column", "none"),symm = FALSE,scale = c("none","row", "column"),na.rm = TRUE,revC = identical(Colv,"Rowv"),add.expr,
                                   breaks,symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",col = "heat.colors",colsep,rowsep,sepcolor = "white",sepwidth = c(0.05, 0.05),cellnote,notecex = 1,notecol = "cyan",na.color = par("bg"),
                                   trace = c("none", "column","row", "both"),tracecol = "cyan",hline = median(breaks),vline = median(breaks),linecol = tracecol,margins = c(5,5),ColSideColors,RowSideColors,side.height.fraction=0.3,cexRow = 0.2 + 1/log10(nr),cexCol = 0.2 + 1/log10(nc),
                                   labRow = NULL,labCol = NULL,key = TRUE,keysize = 1.5,density.info = c("none", "histogram", "density"),denscol = tracecol,symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                                   densadj = 0.25,main = NULL,xlab = NULL,ylab = NULL,lmat = NULL,lhei = NULL,lwid = NULL,ColSideColorsSize = 1,RowSideColorsSize = 1,KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE,cex.axis=drug_label_size)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

final_targets<-unique(unlist(strsplit(n1druglibrary_annotation_CAPS$target,",")))


DREAM_aREA_drugTargets<-function(matrix,database,DB_targets,top_drug_num=25,stat_sig=c("all_sig","inhibit_sig","activate_sig","top_num")){

  trim25_pangea_matrix<-matrix
  
  #######################################
  #RANKED DRUGS & LINES PLOTTING:
  #######################################
  RANKED_DRUGS_all<-names(sort(rowMeans(trim25_pangea_matrix,na.rm=T),decreasing=FALSE))
  RANKED_LINES_all<-names(sort(colMeans(trim25_pangea_matrix[RANKED_DRUGS_all[1:20],],na.rm=T),decreasing=FALSE))
  
  if (stat_sig=="activate_sig"){
    RANKED_LINES_all<-names(sort(colMeans(trim25_pangea_matrix[names(sort(rowMeans(trim25_pangea_matrix,na.rm=T),decreasing=T))[1:20],],na.rm=T),decreasing=T))
    
  }
  
  matrix<-trim25_pangea_matrix[RANKED_DRUGS_all,RANKED_LINES_all]
  
  
  
  avg_drug_rank<-sort(rowMeans(trim25_pangea_matrix,na.rm=T),decreasing=T)
  
  
  
  
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #0 RENAME MATRIX to names compatible with database 
  rownames(database)<-toupper(rownames(database))
  rownames(matrix)<-toupper(rownames(matrix))
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(DB_targets),ncol=length(rownames(matrix)))
  rownames(rlab)<-DB_targets
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  #target="MAP2K1"
  for (target in DB_targets){
    
    #1 extract drug-set given target
    targets_list<-strsplit(database$target,",")
    target_drug_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    target_drug_names<-rownames(database[target_drug_idx,])
    rlab_class_idx<-which(colnames(rlab) %in% target_drug_names)
    
    #3 replace "white" w/ color assigned to each class:
    rlab[target,rlab_class_idx]<-"black"
  }
  
  
  
  
  #######################################
  #RANKED DRUG-CLASSES & LINES PLOTTING:
  #######################################
  #extract average NES-vector to run aREA on:
  avg_nes_gene_vector<-as.matrix(sort(rowMeans(matrix,na.rm=T),decreasing=FALSE))
  
  drug_class_regulon<-list()
  
  #MAKE drug-class regulon:
  #target="MAP2K1"
  for (target in DB_targets){
    
    
    targets_list<-strsplit(database$target,",")
    target_drug_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    target_drug_names<-rownames(database[target_drug_idx,])
    target_drug_names_filtered<-target_drug_names[which(target_drug_names %in% rownames(matrix))]
    num_drugs<-length(target_drug_names_filtered)
    
    #3 define regulon:
    drug_class_regulon[[target]]<-list()
    drug_class_regulon[[target]]$tfmode<-rep(1,num_drugs)
    drug_class_regulon[[target]]$likelihood<-rep(1,num_drugs)
    names(drug_class_regulon[[target]]$tfmode)<-target_drug_names_filtered
    
  }
  
  
  ################
  #USE: STAT-SIG or RANKING
  ################
  #RUN aREA to get ranking:
  aREA_object<-aREA(avg_nes_gene_vector,drug_class_regulon,minsize = 1)
  ranked_drug_classes_values<-sort(aREA_object$nes[,1],decreasing=F)
  ranked_drug_classes<-names(ranked_drug_classes_values)
  
  #PLOT ALL STATISTICALLY SIGNIFICANT DRUGS:
  if (stat_sig=="all_sig"){
    rlab_sorted<-rlab[ranked_drug_classes[which(abs(ranked_drug_classes_values)>1.96)],]
  }
  #PLOT TOP # OF DRUGS
  if(stat_sig=="top_num"){
    rlab_sorted<-rlab[ranked_drug_classes[1:top_drug_num],]
  }
  
  #PLOT ONLY INHIBITING STAT-SIG DRUGS:
  if (stat_sig=="inhibit_sig"){
    rlab_sorted<-rlab[ranked_drug_classes[which(ranked_drug_classes_values<(-1.96))],]
  }
  
  #PLOT ONLY ACTIVATING STAT-SIG DRUGS:
  if (stat_sig=="activate_sig"){
    rlab_sorted<-rlab[ranked_drug_classes[which(ranked_drug_classes_values>-10)],]
  }
  
  
  label_size=min(max(6/nrow(rlab_sorted),1),2.75)
  
  
  
  
  
  
  
  #########################################
  #COLUMN TOP-BAR: cell-lines:
  #########################################
  
  
  
  
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(stat_sig=="activate_sig"){
    heatmap.3_sidebarLabel(matrix[(nrow(matrix):1),], drug_label_size=label_size,
                           Rowv=NA,
                           Colv=NA, 
                           na.rm = TRUE, scale="none", dendrogram="none", margins=c(13,8),
                           RowSideColors=rlab_sorted[nrow(rlab_sorted):1,ncol(rlab_sorted):1], symbreaks=FALSE, key=TRUE, symkey=FALSE,
                           density.info="none", trace="none", main=NA, 
                           cexRow=1,cexCol=1.5,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
                           RowSideColorsSize=45, KeyValueName="NES",keysize = 0.8,xlab=NA,ylab=NA)
  }else{
    heatmap.3_sidebarLabel(matrix, drug_label_size=label_size,
                           Rowv=NA,
                           Colv=NA, 
                           na.rm = TRUE, scale="none", dendrogram="none", margins=c(10,6),
                           RowSideColors=rlab_sorted, symbreaks=TRUE, key=TRUE, symkey=TRUE,
                           density.info="none", trace="none", main=NA, 
                           cexRow=0.2,cexCol=1.5,col=colorRampPalette(c("blue", "blue","white", "red", "red"))(n = 20),
                           RowSideColorsSize=15, KeyValueName="NES",keysize = 0.8,xlab="cell-lines",ylab="Drugs")
  }
  
  
  
  return(avg_drug_rank)
}

####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################

##################################
#GSEA metric = SC2
##################################
#final targets include representative Dr

#RAW SCORE:
nes_avg<-DREAM_aREA_drugTargets(matrix=finalScores_drugs_ceiling[,c("Atom","netphar","SBNB")],
                       database=n1druglibrary_annotation_CTEP,
                       DB_targets = final_targets,
                       stat_sig = "activate_sig")

#RANK-Transformed Score:
rank_avg<-DREAM_aREA_drugTargets(matrix=apply(finalScores_drugs_ceiling[,c("Atom","netphar","SBNB")],2,rank),
                                 database=n1druglibrary_annotation_CTEP,
                                 DB_targets = final_targets,
                                 stat_sig = "activate_sig")


n1druglibrary_annotation_CAPS[names(rank_avg),"DrugBank_targets"]


rank_transform<-apply(finalScores_drugs_ceiling[,-12],2,rank)

write.csv(rank_transform,file="rank-tansformed_SC2.csv")
##################################
#FISHER metric = SC1
##################################


#RAW SCORE:
nes_avg<-DREAM_aREA_drugTargets(matrix=finalScores_drugs_top10_ceiling[,c("Atom","netphar","SBNB")],
                                database=n1druglibrary_annotation_CTEP,
                                DB_targets = final_targets,
                                stat_sig = "activate_sig")

#RANK-Transformed Score:
rank_avg<-DREAM_aREA_drugTargets(matrix=apply(finalScores_drugs_top10_ceiling[,c("Atom","netphar","SBNB")],2,rank),
                                 database=n1druglibrary_annotation_CTEP,
                                 DB_targets = final_targets,
                                 stat_sig = "activate_sig")


n1druglibrary_annotation_CAPS[names(rank_avg),"DrugBank_targets"]



####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##      ########    #######     ########      #######    ##              ##      ########### ########    ########   ####    ##                                                                            
##     ##          ##     ##    ##     ##     ##    ##   ##             ####         ##         ##      ##      ##  ## ##   ##                
##    ##           ##     ##    ##     ##     ##    ##   ##            ##  ##        ##         ##      ##      ##  ##  ##  ##                         
##    ##           ##     ##    ########      #######    ##           ##    ##       ##         ##      ##      ##  ##   ## ##                 
##    ##           ##     ##    ##    ##      ##   ##    ##          ##########      ##         ##      ##      ##  ##    ####                    
##     ##          ##     ##    ##     ##     ##    ##   ##         ##        ##     ##         ##      ##      ##  ##     ###                  
##      ########    #######     ##      ##    ##     ##  ########  ##          ##    ##      ########    ########   ##      ##                             
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#UNDERSTAND CORRELATION STRUCTURE: kinome-binding vs MR-perturbations (across overlapping 88 drugs)




load("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-06-03 KINASE_targets_effectors/simple_PANGEA-HALLMARKS_KinomeBinding_overlap.RData")

CORRmatrix_kinomeBIND_effectorACTIVITY<-cor(kinase_inhibitors_Binding,kinase_inhibitors_avgHALLMARKS,method="spearman")

plot(CORRmatrix_kinomeBIND_effectorACTIVITY["MAP2K2",],CORRmatrix_kinomeBIND_effectorACTIVITY["BRAF",])



##################################################################
#PLOTTING FUNCTIONS:
##################################################################
library("gplots")
library("RColorBrewer")
library("devtools")

hcluster<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
  #distance matrix:
  dist_columns <- dist(t(matrix))
  dist_rows <- dist(matrix)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "-log10(Kd)",key.xlab = "-log10(Kd)")
}
hcluster_corr<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
  #distance matrix:
  cor_columns <- cor(matrix,method="spearman")
  cor_rows <- cor(t(matrix),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(15,15), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "NES",key.xlab = "NES")
}

na_idx<-which(is.na(rowMeans(CORRmatrix_kinomeBIND_effectorACTIVITY))==TRUE)
CORRmatrix_kinomeBIND_effectorACTIVITY_na<-CORRmatrix_kinomeBIND_effectorACTIVITY[-na_idx,]


hcluster_corr(CORRmatrix_kinomeBIND_effectorACTIVITY_na,col_size = 0.5,row_size=0.4)

ranked_genes<-names(sort(apply(CORRmatrix_kinomeBIND_effectorACTIVITY_na,2,var),decreasing=TRUE))
ranked_kinases<-names(sort(apply(CORRmatrix_kinomeBIND_effectorACTIVITY_na,1,var),decreasing=TRUE))


hcluster_corr(t(CORRmatrix_kinomeBIND_effectorACTIVITY_na[ranked_kinases[1:200],ranked_genes[1:40]]),col_size = 0.6,row_size=1,title="Kinome-N1-Correlation:  top-2000-var genes")



##################################################################
#GENOME & KINOME ANNOTATION:
##################################################################
library("gplots")
library("RColorBrewer")
library("devtools")
heatmap.3 <- function(x,Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,distfun = dist,hclustfun = hclust,dendrogram = c("both","row", "column", "none"),symm = FALSE,scale = c("none","row", "column"),na.rm = TRUE,revC = identical(Colv,"Rowv"),add.expr,
                      breaks,symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",col = "heat.colors",colsep,rowsep,sepcolor = "white",sepwidth = c(0.05, 0.05),cellnote,notecex = 1,notecol = "cyan",na.color = par("bg"),
                      trace = c("none", "column","row", "both"),tracecol = "cyan",hline = median(breaks),vline = median(breaks),linecol = tracecol,margins = c(5,5),ColSideColors,RowSideColors,side.height.fraction=0.3,cexRow = 0.2 + 1/log10(nr),cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,labCol = NULL,key = TRUE,keysize = 1.5,density.info = c("none", "histogram", "density"),denscol = tracecol,symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,main = NULL,xlab = NULL,ylab = NULL,lmat = NULL,lhei = NULL,lwid = NULL,ColSideColorsSize = 1,RowSideColorsSize = 1,KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

load("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-06-03 KINASE_targets_effectors/KINOME-classifications.RData")

#Symetric Key:  NES's
hcluster_corr_MULTI_sidebar_G2M_mod2<-function(matrix,database,DB_col,DB_class,class_names,class_colors,main_title,G2M_plus=TRUE,x_lab="Genes"){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  #0 RENAME MATRIX to names compatible with database 
  rownames(database)<-make.names(toupper(rownames(database)),unique = TRUE)
  
  rownames(matrix)<-make.names(toupper(rownames(matrix)),unique = TRUE)
  
  
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix<-matrix[rownames(database),]
  
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(class_names),ncol=length(rownames(matrix)))
  rownames(rlab)<-class_names
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(class_names)){
    #1 extract 1st rows DB-entry
    db_col_tmp<-DB_col[row]
    db_class_tmp<-DB_class[row]
    
    #2 extract class-indices in rlab
    db_class_idx<-which(database[,db_col_tmp] == db_class_tmp)
    class_names_tmp<-rownames(database[db_class_idx,])
    rlab_class_idx<-which(colnames(rlab) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-class_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  
  
  #########################################
  #XX G2M TOP-bar
  #########################################
  
  ###1 extract relevant G2M genes:
  G2M_chckpt<-c("TOP2A","CENPF","CHEK1","FOXM1","CENPI","BLM","BRIP1","MCM6","MCM4","E2F1","E2F2","E2F3","E2F8","BUB1B","MYBL2","MCM7","FEN1","AURKA","AURKB","CENPU","FANCA","TYMS","RAD51","HMGA1","DNMT3B","E2F7","BRCA1","BRCA2","PRKDC","CENPK","MCM3","DMNT1","FANCD2","RAN","RANBP1","MCM2","CENPZ","RACGAP1","PCNA","H2AFX","MCM5")
  G2M_chckpt_plus<-unique(c("TOP2A","CENPF","CHEK1","FOXM1","CENPI","BLM","BRIP1","MCM6","MCM4","E2F1","E2F2","E2F3","E2F8","BUB1B","MYBL2","MCM7","FEN1","AURKA","AURKB","CENPU","FANCA","TYMS","RAD51","HMGA1","DNMT3B","E2F7","BRCA1","BRCA2","PRKDC","CENPK","MCM3","DMNT1","FANCD2","RAN","RANBP1","MCM2","CENPZ","RACGAP1","PCNA","H2AFX","MCM5","TRIP13","ASF1B","UHRF1","ATAD2","ECT2","TTK","WHSC1","ARHGAP11A","ZNF367","TCF19","HELLS","RACGAP1","MKI67","SPAG5","UBE2C","TONSL","RUVBLT","ZWINT","HMGB2","DEPDC1","CBX3","RFC4","CHAF1A","GTSE1","HMGN1","PSMC3IP","HDGF","CKS1B","TCEB1","CKS2","CCNA2","ENY2","GGCT","ILE2","HSPE1","PSRC1","NDC80","TMPO","CDK2","GMPS","EIF4EBP1","PA2G4","GMNN","ZC3H15","PTTG1","ZC3H1","CYCS","GTF2E2","HLTF","TCF19","WHSC1","TTF2","SHOX2","EIF2AK2","ARNTL2","UHRF1","GTSE1","DEPDC1","SUV39H2","ZWINT","YEASTS2","LRPPRC","NAA15","POLD1","RFC4","RBL1","MCM3","ILF3","CDK7","CIT","GLRA2","PHF19","CDCA7","ZNF695","EZH2","PCNA","HNRNPAB","ZCRB1","SUB1","FH","HSPE1","RAB1A","CCDC47","COPS3","H2AFX","DAXX","HDGF","APH1A","PUF60","POLR3K","HSBP1","IDH2","VPS72","PRMT1","UBE2L3","PPP2CA","HSPA5","MCM2","HDAC2","CEBPZ","TIMELESS","UBEC2","PTTG1","TRAIP","MNX1","ARHGEF39","CHEK2","IQGAP3","HSPD1"))
  if (G2M_plus == TRUE){
    G2M_chckpt_filtered<-G2M_chckpt_plus[which(G2M_chckpt_plus %in% colnames(matrix))]
  }
  else {
    G2M_chckpt_filtered<-G2M_chckpt[which(G2M_chckpt %in% colnames(matrix))]
  }
  
  ###2 MAKE column color-key:
  clab<-matrix("white",nrow=length(colnames(matrix)),ncol=1)
  
  matrix_G2M_idx<-which(colnames(matrix) %in% G2M_chckpt_filtered)
  
  clab[matrix_G2M_idx,1]<-"darkblue"
  
  colnames(clab)=c("G2M")
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  cor_columns <- cor(matrix,method="spearman")
  cor_rows <- cor(t(matrix),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  heatmap.3(matrix, 
            Rowv=as.dendrogram(hclust_rows),
            Colv=as.dendrogram(hclust_columns), 
            na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,6),
            ColSideColors=clab, RowSideColors=rlab, symbreaks=TRUE, key=TRUE, symkey=TRUE,
            density.info="none", trace="none", main=main_title, cexRow=0.6,cexCol=0.1,col=colorRampPalette(c("blue",  "white",  "red"))(n = 20),
            ColSideColorsSize=1, RowSideColorsSize=5, KeyValueName="spearman",keysize = 0.8,xlab=x_lab,ylab="kinases")
  
}

hcluster_corr_MULTI_sidebar_G2M_mod2_transpose<-function(matrix,database,DB_col,DB_class,class_names,class_colors,main_title,G2M_plus=TRUE,x_lab="Genes"){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  #0 RENAME MATRIX to names compatible with database 
  rownames(database)<-make.names(toupper(rownames(database)),unique = TRUE)
  
  rownames(matrix)<-make.names(toupper(rownames(matrix)),unique = TRUE)
  
  
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix<-matrix[rownames(database),]
  
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(class_names),ncol=length(rownames(matrix)))
  rownames(rlab)<-class_names
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(class_names)){
    #1 extract 1st rows DB-entry
    db_col_tmp<-DB_col[row]
    db_class_tmp<-DB_class[row]
    
    #2 extract class-indices in rlab
    db_class_idx<-which(database[,db_col_tmp] == db_class_tmp)
    class_names_tmp<-rownames(database[db_class_idx,])
    rlab_class_idx<-which(colnames(rlab) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-class_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  clab<-t(rlab)
  matrix<-t(matrix)
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  cor_columns <- cor(matrix,method="spearman")
  cor_rows <- cor(t(matrix),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  heatmap.3(matrix, 
            Rowv=as.dendrogram(hclust_rows),
            Colv=as.dendrogram(hclust_columns), 
            na.rm = TRUE, scale="none", dendrogram="both", margins=c(3,20),
            ColSideColors=clab, symbreaks=TRUE, key=TRUE, symkey=TRUE,
            density.info="none", trace="none", main=main_title, cexRow=0.9,cexCol=0.3,col=colorRampPalette(c("blue",  "white",  "red"))(n = 20),
            ColSideColorsSize=5, KeyValueName="spearman",keysize = 0.8,xlab="kinases",)
  
}




hcluster_corr_MULTI_sidebar_G2M_mod2(CORRmatrix_kinomeBIND_effectorACTIVITY_na,kinome_classification_hugo_filtered,
                                     DB_col = c(3,3,3,3,3,3,3),
                                     DB_class = c("AGC","CAMK","CK1","CMGC","STE","TK","TKL"),
                                     class_names =c("AGC","CAMK","CK1","CMGC","STE","TK","TKL"),
                                     class_colors = c("darkred","red","darkorange","darkgreen","cyan4","darkblue","blueviolet"),
                                     main_title = "KINOME-binding/N1-screen Correlation Matrix:")



hcluster_corr_MULTI_sidebar_G2M_mod2_transpose(CORRmatrix_kinomeBIND_effectorACTIVITY_na[,ranked_genes[1:40]],kinome_classification_hugo_filtered,
                                     DB_col = c(3,3,3,3,3,3,3),
                                     DB_class = c("AGC","CAMK","CK1","CMGC","STE","TK","TKL"),
                                     class_names =c("AGC","CAMK","CK1","CMGC","STE","TK","TKL"),
                                     class_colors = c("darkred","red","darkorange","darkgreen","cyan4","darkblue","blueviolet"),
                                     main_title = "KINOME-binding/N1-screen Correlation Matrix:")





####################################################################################################################################
#TRANSFORM KINASES INTO KEGG PATHWAYS
####################################################################################################################################

library(viper)
load("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/11 Kinome-Explanation/kegg_pathways_regulon.RData")

#MAKE REGULON BELOW:
{
###########################
#1  IMPORT KEGG
###########################


library("GSA")


data <-  GSA.read.gmt("c2.cp.kegg.v7.1.symbols.gmt")

data$genesets

data$geneset.names

###########################
#2  MAKE HALLMARK REGULON:
###########################



kegg_pathways_regulon<-list()

#idx=1
for (idx in 1:length(data$geneset.names)){
  hallmark<-gsub("KEGG_","",data$geneset.names[idx])
  hallmark_geneset<-data$genesets[[idx]]
  
  
  
  #2 MAKE REGULON:
  column_gene_set<-hallmark_geneset
  num_genes_wn_set<-length(hallmark_geneset)
  
  kegg_pathways_regulon[[hallmark]]<-list()
  kegg_pathways_regulon[[hallmark]]$tfmode<-rep(1,num_genes_wn_set)
  kegg_pathways_regulon[[hallmark]]$likelihood<-rep(1,num_genes_wn_set)
  names(kegg_pathways_regulon[[hallmark]]$tfmode)<-column_gene_set
}


save(kegg_pathways_regulon,file="kegg_pathways_regulon.RData")


}


KinasePathway_mRNAHallmark_corr<-aREA(CORRmatrix_kinomeBIND_effectorACTIVITY_na,kegg_pathways_regulon,minsize = 2)


hcluster_corr(KinasePathway_mRNAHallmark_corr$nes,col_size = 0.6,row_size=0.5,title="Kinome-N1-Correlation:  top-2000-var genes")


sort(rownames(KinasePathway_mRNAHallmark_corr$nes))

KinasePathway_mRNAHallmark_corr$nes



                        
kegg_relevant<-c("APOPTOSIS","CALCIUM_SIGNALING_PATHWAY","CELL_CYCLE","CHEMOKINE_SIGNALING_PATHWAY","CYTOSOLIC_DNA_SENSING_PATHWAY",
                 "ERBB_SIGNALING_PATHWAY","INSULIN_SIGNALING_PATHWAY","JAK_STAT_SIGNALING_PATHWAY","MAPK_SIGNALING_PATHWAY",
                 "MTOR_SIGNALING_PATHWAY","NUCLEOTIDE_EXCISION_REPAIR","P53_SIGNALING_PATHWAY","PHOSPHATIDYLINOSITOL_SIGNALING_SYSTEM",
                 "REGULATION_OF_AUTOPHAGY","TGF_BETA_SIGNALING_PATHWAY","TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","VEGF_SIGNALING_PATHWAY","WNT_SIGNALING_PATHWAY")

hallmark_relevant<-c("APOPTOSIS","DNA_REPAIR","E2F_TARGETS","FATTY_ACID_METABOLISM","G2M_CHECKPOINT","GLYCOLYSIS","HEDGEHOG_SIGNALING",
  "HYPOXIA","IL2_STAT5_SIGNALING","IL6_JAK_STAT3_SIGNALING","INFLAMMATORY_RESPONSE",
  "INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","KRAS_SIGNALING_DN","KRAS_SIGNALING_UP","MITOTIC_SPINDLE",
  "MTORC1_SIGNALING","MYC_TARGETS_V1","MYC_TARGETS_V2","NOTCH_SIGNALING",
  "OXIDATIVE_PHOSPHORYLATION","P53_PATHWAY","PANCREAS_BETA_CELLS","PEROXISOME","PI3K_AKT_MTOR_SIGNALING",
  "PROTEIN_SECRETION","REACTIVE_OXIGEN_SPECIES_PATHWAY","TGF_BETA_SIGNALING","TNFA_SIGNALING_VIA_NFKB",
  "UNFOLDED_PROTEIN_RESPONSE","WNT_BETA_CATENIN_SIGNALING")

hcluster_corr(KinasePathway_mRNAHallmark_corr$nes[kegg_relevant,hallmark_relevant],col_size = 1,row_size=1,title="KEGG-Pathways vs mRNA-Programs")

sort(KinasePathway_mRNAHallmark_corr$nes[kegg_relevant,"MYC_TARGETS_V1"])
sort(KinasePathway_mRNAHallmark_corr$nes[kegg_relevant,"APOPTOSIS"])


kegg_relevant2<-c("APOPTOSIS","CALCIUM_SIGNALING_PATHWAY","CELL_CYCLE",
                 "ERBB_SIGNALING_PATHWAY","INSULIN_SIGNALING_PATHWAY","JAK_STAT_SIGNALING_PATHWAY","MAPK_SIGNALING_PATHWAY",
                 "MTOR_SIGNALING_PATHWAY","P53_SIGNALING_PATHWAY","PHOSPHATIDYLINOSITOL_SIGNALING_SYSTEM",
                 "TGF_BETA_SIGNALING_PATHWAY","VEGF_SIGNALING_PATHWAY")

kegg_relevant3<-c("CALCIUM_SIGNALING_PATHWAY","CELL_CYCLE",
                  "ERBB_SIGNALING_PATHWAY","INSULIN_SIGNALING_PATHWAY","JAK_STAT_SIGNALING_PATHWAY",
                  "P53_SIGNALING_PATHWAY",
                  "TGF_BETA_SIGNALING_PATHWAY","VEGF_SIGNALING_PATHWAY")

hcluster(KinasePathway_mRNAHallmark_corr$nes[kegg_relevant3,c("MYC_TARGETS_V1","MYC_TARGETS_V2","E2F_TARGETS","G2M_CHECKPOINT","TGF_BETA_SIGNALING","APOPTOSIS")],col_size = 1,row_size=1,title="KEGG-Pathways vs mRNA-Programs")













####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##    ########     ########    ##     ##     ########             ########       ##      ###     ##   ##     ## 
##    ##     ##    ##     ##   ##     ##    ##                    ##     ##     ####     ####    ##   ##    ##               
##    ##      ##   ##     ##   ##     ##   ##                     ##     ##    ##  ##    ## ##   ##   ##   ##                 
##    ##      ##   ########    ##     ##   ##    ####    ######   ########    ##    ##   ##  ##  ##   ######         
##    ##      ##   ##    ##    ##     ##   ##       ##            ##     ##  ##########  ##   ## ##   ##   ##                                  
##    ##     ##    ##     ##   ##     ##    ##      ##            ##     ##  ##      ##  ##    ####   ##    ##
##    ########     ##      ##   #######      ########             ########   ##      ##  ##     ###   ##     ##                                    
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      ########    #######      ####    ##    ########  ##     ##   #######    ##    #######    ###     ##                                                
##     ##          ##     ##     ## ##   ##    ##        ##     ##  ##          ##   ##     ##   ####    ##                                                  
##    ##           ##     ##     ##  ##  ##    ##        ##     ##  ##          ##   ##     ##   ## ##   ##                                                           
##    ##           ##     ##     ##   ## ##    #######   ##     ##   #######    ##   ##     ##   ##  ##  ##                                                            
##    ##           ##     ##     ##    ####    ##        ##     ##         ##   ##   ##     ##   ##   ## ##                                                     
##     ##          ##     ##     ##     ###    ##        ##     ##         ##   ##   ##     ##   ##    ####                                                   
##      ########    #######      ##      ##    ##         #######    #######    ##    #######    ##     ###                                                                         
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#NOT ENOUGH TEAMS TO DO DEEP DIVE IN METHODS OR DATA BASES AS:
#         -> data sources quite diverse (each person had a unique data-set and so hard to deconvolue NORM & Data biases)
#GOAL:  #1 analyze drugbank vs kinome gold-standards
#       #2 see if non-drug-bank targets picked up in submissions...




setwd("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/")


#GOLD STANDARD:
load("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-05-24 - N1screen - drug-clean-ness/OVERLAP_n1-Kinome-Screens.RData")
rm(kinome_selectivity_scores_OVERLAP)


load("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/final_results_processed.RData")


####################################################################################################################################
#1 RECOVERY OF DRUG-BANK TARGETS
#####################################################################################################################################
kinome_GoldStandard<-kinome_inhibitor_matrix_pKa_OVERLAP[rownames(final_results_processed$SBNB),which(colnames(kinome_inhibitor_matrix_pKa_OVERLAP) %in% colnames(final_results_processed$SBNB))]

max(kinome_GoldStandard)

threshold_matrix<-matrix(ncol=length(seq(4.6,8.5,0.01)),nrow=length(rownames(final_results_processed$SBNB)))
rownames(threshold_matrix)<-rownames(final_results_processed$SBNB)
colnames(threshold_matrix)<-seq(4.6,8.5,0.01)

#drug="SORAFENIB"
for (drug in rownames(threshold_matrix)){
  cutoff="8.5"
  for (cutoff in colnames(threshold_matrix)){
    #ID DrugBank Kinases (for this drug):
    targets_list<-strsplit(n1druglibrary_annotation_CTEP_OVERLAP[drug,"target"],",")[[1]]
    KINASES_drugbank<-targets_list[which(targets_list %in% colnames(kinome_GoldStandard))]
    
    #ID KINOME Drugs (for this kinase)
    KINASES_kinome<-names(which(kinome_GoldStandard[drug,]>as.numeric(cutoff)))
    
    #FRACTION DRUG-BANK KINASES Recovered
    drugbank_recovered<-length(intersect(KINASES_kinome,KINASES_drugbank))/length(KINASES_drugbank)
    
    #STORE:
    threshold_matrix[drug,cutoff]<-drugbank_recovered
  }
}

plot(as.numeric(colnames(threshold_matrix)),colMeans(threshold_matrix),type="n",main="Recovered DrugBank Targets Based on Kd Threshold",
     ylab="Fraction Drug-Bank Kinase Targets",xlab="Affinity Threshold Used (-log10(Kd))",cex.main=1.5,cex.lab=1.5,cex.axis=1.5,xlim=c(8.5,4.5))
lines(as.numeric(colnames(threshold_matrix)),colMeans(threshold_matrix),lwd=3)
abline(v=6,col="red",lty=2)

####################################################################################################################################
#  Kinases-Drug <uM events NOT in DrugBank:
#####################################################################################################################################
kinome_inhibitor_matrix_pKa_OVERLAP_filtered<-kinome_inhibitor_matrix_pKa_OVERLAP[,which(colnames(kinome_inhibitor_matrix_pKa_OVERLAP) %in% colnames(final_results_processed$SBNB))]

drug_unique_num<-c()
drug_overlap_num<-c()

Bootstrapped_targets<-matrix(nrow=nrow(final_results_processed$SBNB),ncol=6)
rownames(Bootstrapped_targets)<-rownames(final_results_processed$SBNB)
colnames(Bootstrapped_targets)<-c("Atom_overlap","Atom_kinome","netphar_overlap","netphar_kinome","SBNB_overlap","SBNB_kinome")

#drug="DASATINIB"
for (drug in rownames(final_results_processed$SBNB)){
  #################################
  #1 Characterize Gold-Standard Overlap:
  #################################
  
  #ID DrugBank Kinases (for this drug):
  targets_list<-strsplit(n1druglibrary_annotation_CTEP_OVERLAP[drug,"target"],",")[[1]]
  KINASES_drugbank<-targets_list[which(targets_list %in% colnames(kinome_inhibitor_matrix_pKa_OVERLAP_filtered))]

  #ID KINOME Drugs (for this kinase)
  KINASES_kinome<-names(which(kinome_inhibitor_matrix_pKa_OVERLAP_filtered[drug,]>6))
  
  #PERCENT OF DrugBank drugs w/n kinome:

  tmp_num<-(length(KINASES_kinome)-length(which(KINASES_drugbank %in% KINASES_kinome)))
  drug_unique_num<-append(drug_unique_num,tmp_num)
  
  tmp_overlap<-length(which(KINASES_drugbank %in% KINASES_kinome))
  drug_overlap_num<-append(drug_overlap_num,tmp_overlap)
  
  
  #################################
  #2 POPULATE RANK-MATRIX OF PREDICTIONS
  #################################
  
  Bootstrapped_targets[drug,"Atom_overlap"]<-mean(final_results_processed$Atom[drug,KINASES_drugbank])/mean(final_results_processed$Atom[drug,])
  Bootstrapped_targets[drug,"Atom_kinome"]<-mean(final_results_processed$Atom[drug,KINASES_kinome])/mean(final_results_processed$Atom[drug,])
  
  
  Bootstrapped_targets[drug,"SBNB_overlap"]<-mean(final_results_processed$SBNB[drug,KINASES_drugbank])/mean(final_results_processed$SBNB[drug,])
  Bootstrapped_targets[drug,"SBNB_kinome"]<-mean(final_results_processed$SBNB[drug,KINASES_kinome])/mean(final_results_processed$SBNB[drug,])
  
  
  Bootstrapped_targets[drug,"netphar_overlap"]<-mean(final_results_processed$netphar[drug,KINASES_drugbank])/mean(final_results_processed$netphar[drug,])
  Bootstrapped_targets[drug,"netphar_kinome"]<-mean(final_results_processed$netphar[drug,KINASES_kinome])/mean(final_results_processed$netphar[drug,])
  
}

##########################
#TARGET-DATA OVERLAP: bar-plot
##########################
names(drug_unique_num)<-rownames(final_results_processed$SBNB)
names(drug_overlap_num)<-rownames(final_results_processed$SBNB)


DrugBank_Kinome_compare_drugs<-cbind(drug_overlap_num,drug_unique_num)
names_order<-names(sort(rowSums(DrugBank_Kinome_compare_drugs),decreasing=T))

names_order2<-names(sort(drug_overlap_num,decreasing=T))


par(las=2) # make label text perpendicular to axis
par(mar=c(9,5,4,5),cex=1.2) # increase y-axis margin.
#PLOT:
barplot(t(DrugBank_Kinome_compare_drugs[names_order,]),horiz=F,cex.lab=2,
        ylab="# of Kinase-Targets",xlab="",main="Gold-Standard Comparisons:",
        cex.axis = 1.5,col=c("black","red"),cex.main=2)
legend("topright",c("DrugBank-Kinome Overlap","Unique Kinome-Targets"),
       fill=c("black","red"),text.col=c("black","red"),cex=1.7)

save(list=c("DrugBank_Kinome_compare_drugs","Bootstrapped_targets"),file="KINOME-DrugBank_library-compare.RData")

##########################
#PREDICTIONS COMPARISONS:
##########################

plot(Bootstrapped_targets[,"Atom_overlap"],Bootstrapped_targets[,"Atom_kinome"], col="white"
     ,xlab="DrugBank-Targets Average-Score (average)",ylab="Unique Kinome-Targets Avg-Score")
text(Bootstrapped_targets[,"Atom_overlap"],Bootstrapped_targets[,"Atom_kinome"], labels=rownames(Bootstrapped_targets))


plot(Bootstrapped_targets[,"netphar_overlap"],Bootstrapped_targets[,"netphar_kinome"], col="white"
     ,xlab="DrugBank-Targets Average-Score (average)",ylab="Unique Kinome-Targets Avg-Score")
text(Bootstrapped_targets[,"netphar_overlap"],Bootstrapped_targets[,"netphar_kinome"], labels=rownames(Bootstrapped_targets))

plot(Bootstrapped_targets[,"SBNB_overlap"],Bootstrapped_targets[,"SBNB_kinome"], col="white"
     ,xlab="DrugBank-Targets Average-Score (average)",ylab="Unique Kinome-Targets Avg-Score")
text(Bootstrapped_targets[,"SBNB_overlap"],Bootstrapped_targets[,"SBNB_kinome"], labels=rownames(Bootstrapped_targets))

#CONCENSUS:
plot(rowMeans(apply(Bootstrapped_targets[,grep("overlap",colnames(Bootstrapped_targets))],2,rank)),
     rowMeans(apply(Bootstrapped_targets[,grep("kinome",colnames(Bootstrapped_targets))],2,rank)), 
     col="white",xlab="DrugBank-Targets Average-Score (Rank)",ylab="Unique Kinome-Targets Avg-Score (Rank)")
text(rowMeans(apply(Bootstrapped_targets[,grep("overlap",colnames(Bootstrapped_targets))],2,rank)),
     rowMeans(apply(Bootstrapped_targets[,grep("kinome",colnames(Bootstrapped_targets))],2,rank)),
     labels=rownames(Bootstrapped_targets))











####################################################################################################################################
#  Kinases-Drug <uM events NOT in DrugBank:
#####################################################################################################################################

kinome_unique_fract<-c()
kinome_unique_num<-c()
overlap_num<-c()


#target="EGFR"
for (target in colnames(kinome_inhibitor_matrix_pKa_OVERLAP)){
  #ID DrugBank Drugs (for this kinase):
  targets_list<-strsplit(n1druglibrary_annotation_CTEP_OVERLAP$target,",")
  target_drug_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
  DRUGS_drugbank<-rownames(n1druglibrary_annotation_CTEP_OVERLAP[target_drug_idx,])
  
  #ID KINOME Drugs (for this kinase)
  DRUGS_kinome<-names(which(kinome_inhibitor_matrix_pKa_OVERLAP[,target]>6))
  
  #PERCENT OF DrugBank drugs w/n kinome:
  tmp_fract<-(length(DRUGS_kinome)-length(which(DRUGS_drugbank %in% DRUGS_kinome)))/length(DRUGS_kinome)
  kinome_unique_fract<-append(kinome_unique_fract,tmp_fract)
  
  tmp_num<-(length(DRUGS_kinome)-length(which(DRUGS_drugbank %in% DRUGS_kinome)))
  kinome_unique_num<-append(kinome_unique_num,tmp_num)
  
  tmp_overlap<-length(which(DRUGS_drugbank %in% DRUGS_kinome))
  overlap_num<-append(overlap_num,tmp_overlap)
  
}
names(kinome_unique_num)<-colnames(kinome_inhibitor_matrix_pKa_OVERLAP)
names(kinome_unique_fract)<-colnames(kinome_inhibitor_matrix_pKa_OVERLAP)
names(overlap_num)<-colnames(kinome_inhibitor_matrix_pKa_OVERLAP)


DrugBank_Kinome_compare<-cbind(overlap_num,kinome_unique_num)
names_order<-names(sort(rowSums(DrugBank_Kinome_compare),decreasing=T))



par(las=2) # make label text perpendicular to axis
par(mar=c(5,5,4,5),cex=1) # increase y-axis margin.
#PLOT:
barplot(t(DrugBank_Kinome_compare[names_order[1:40],]),horiz=F,cex.main=1.5,cex.lab=1.5,
        ylab="Correlation",xlab="135 drugs",main="PANACEA-L1000 Correlation:  Drug-Signatures",
        cex.axis = 1.5,col=c("black","red"))


legend("topright",legend_info,cex=1,fill="red",text.col="darkred")












