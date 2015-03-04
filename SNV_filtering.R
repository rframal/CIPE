#Filtrando SNVs

exomas<-read.table("all_fs_qd_toR.csv", header=T, sep="\t")

#Filtrar SNVs com cobertura total >= 10,  )

subexomas<-subset(exomas,
                  (as.numeric(exomas$MIC204S_COVv) + as.numeric(exomas$MIC204S_COVr)) >=10 &
                    (as.numeric(exomas$MIC204T_COVv) + as.numeric(exomas$MIC204T_COVr)) >=10 &
                    (as.numeric(exomas$M168S_COVv) + as.numeric(exomas$M168S_COVr)) >=10 &
                    (as.numeric(exomas$M168T_COVv) + as.numeric(exomas$M168T_COVr)) >=10 &
                    (as.numeric(exomas$MIC192S_COVv) + as.numeric(exomas$MIC192S_COVr)) >=10 &
                    (as.numeric(exomas$MIC192T_COVv) + as.numeric(exomas$MIC192T_COVr)) >=10 &  
                    (as.numeric(exomas$M27S_COVv) + as.numeric(exomas$M27S_COVr)) >=10 &
                    (as.numeric(exomas$M27T_COVv) + as.numeric(exomas$M27T_COVr)) >=10 &
                    (as.numeric(exomas$M56S_COVv) + as.numeric(exomas$M56S_COVr)) >=10 &
                    (as.numeric(exomas$M56T_COVv) + as.numeric(exomas$M56T_COVr)) >=10 &
                    (as.numeric(exomas$M155S_COVv) + as.numeric(exomas$M155S_COVr)) >=10 &
                    (as.numeric(exomas$M155T_COVv) + as.numeric(exomas$M155T_COVr)) >=10 &
                                        
                    exomas$Variant_class2=="MISSENSE" & 
                    exomas$Type_SNV_or_INDEL=="SNV")


subexomas<-rbind(subexomas,
       subset(exomas,
       (as.numeric(exomas$MIC204S_COVv) + as.numeric(exomas$MIC204S_COVr)) >=10 &
         (as.numeric(exomas$MIC204T_COVv) + as.numeric(exomas$MIC204T_COVr)) >=10 &
         (as.numeric(exomas$M168S_COVv) + as.numeric(exomas$M168S_COVr)) >=10 &
         (as.numeric(exomas$M168T_COVv) + as.numeric(exomas$M168T_COVr)) >=10 &
         (as.numeric(exomas$MIC192S_COVv) + as.numeric(exomas$MIC192S_COVr)) >=10 &
         (as.numeric(exomas$MIC192T_COVv) + as.numeric(exomas$MIC192T_COVr)) >=10 &  
         (as.numeric(exomas$M27S_COVv) + as.numeric(exomas$M27S_COVr)) >=10 &
         (as.numeric(exomas$M27T_COVv) + as.numeric(exomas$M27T_COVr)) >=10 &
         (as.numeric(exomas$M56S_COVv) + as.numeric(exomas$M56S_COVr)) >=10 &
         (as.numeric(exomas$M56T_COVv) + as.numeric(exomas$M56T_COVr)) >=10 &
         (as.numeric(exomas$M155S_COVv) + as.numeric(exomas$M155S_COVr)) >=10 &
         (as.numeric(exomas$M155T_COVv) + as.numeric(exomas$M155T_COVr)) >=10 &
         
         exomas$Variant_class2=="NONSENSE" & 
         exomas$Type_SNV_or_INDEL=="SNV")
          ) 


#Trabalhando apenas com o paciente MIC204 

#1) Filtrar SNVs do sangue com freq variante >= 20%  )

MIC204S<-subset(subexomas, subexomas$MIC204S_FreqV >= 20 )[,c(1,2,46)]

#3) Filtrar SNVs do tumor com freq variante >= 20%  )

MIC204T<-subset(subexomas, subexomas$MIC204T_FreqV >= 20 )[,c(1,2,3,46,59)]

#5) Fazendo o merge entre os SNVs de sangue e tumor 

MIC204T_merge<-merge(MIC204S, MIC204T, by=c("chr", "pos"), all.y=T)

colnames(MIC204T_merge)<-c("chr","pos","dbSNPID.sangue","gene","dbSNPID.tumor","cosmic")

#6) #As linhas com NA na coluna dbSNPID.sangue são mutações somáticas

MIC204T_merge2<-MIC204T_merge[is.na(MIC204T_merge$dbSNPID.sangue),]

#8) Obtendo a lista de genes

which(colnames(subexomas)=="MIC204T_COVr")
which(colnames(subexomas)=="MIC204T_COVv")
which(colnames(subexomas)=="MIC204T_FreqV")

MIC204T_list<-merge(subexomas, MIC204T_merge2, by=c("chr", "pos"))[,c(1,2,3,31,32,44,46,56,57,63,67:70)]

#9) Escrevendo os resultados

write.table(file="somatic_mutations_new.txt", MIC204T_list, row.names=F, sep="\t")

##########################################################################################
#Trabalhando apenas com o paciente M168

#1) Filtrar SNVs do sangue com freq variante >= 20%  )

M168S<-subset(subexomas, subexomas$M168S_FreqV >= 20 )[,c(1,2,46)]

#3) Filtrar SNVs do tumor com freq variante >= 20%  )

M168T<-subset(subexomas, subexomas$M168T_FreqV >= 20 )[,c(1,2,3,46,59)]

#5) Fazendo o merge entre os SNVs de sangue e tumor 

M168T_merge<-merge(M168S, M168T, by=c("chr", "pos"), all.y=T)

colnames(M168T_merge)<-c("chr","pos","dbSNPID.sangue","gene","dbSNPID.tumor","cosmic")

#6) #As linhas com NA na coluna dbSNPID.sangue são mutações somáticas

M168T_merge2<-M168T_merge[is.na(M168T_merge$dbSNPID.sangue),]

#8) Obtendo a lista de genes

which(colnames(subexomas)=="M168T_COVr")
which(colnames(subexomas)=="M168T_COVv")
which(colnames(subexomas)=="M168T_FreqV")

M168T_list<-merge(subexomas, M168T_merge2, by=c("chr", "pos"))[,c(1,2,3,15,16,36,46,56,57,63,67:70)]

#9) Escrevendo os resultados

write.table(file="somatic_mutations_new.txt", M168T_list, row.names=F, sep="\t", append=T)

#################################################################################################
#Trabalhando apenas com o paciente MIC192

#1) Filtrar SNVs do sangue com freq variante >= 20%  )

MIC192S<-subset(subexomas, subexomas$MIC192S_FreqV >= 20 )[,c(1,2,46)]

#3) Filtrar SNVs do tumor com freq variante >= 20%  )

MIC192T<-subset(subexomas, subexomas$MIC192T_FreqV >= 20 )[,c(1,2,3,46,59)]

#5) Fazendo o merge entre os SNVs de sangue e tumor 

MIC192T_merge<-merge(MIC192S, MIC192T, by=c("chr", "pos"), all.y=T)

colnames(MIC192T_merge)<-c("chr","pos","dbSNPID.sangue","gene","dbSNPID.tumor","cosmic")

#6) #As linhas com NA na coluna dbSNPID.sangue são mutações somáticas

MIC192T_merge2<-MIC192T_merge[is.na(MIC192T_merge$dbSNPID.sangue),]

#8) Obtendo a lista de genes

which(colnames(subexomas)=="MIC192T_COVr")
which(colnames(subexomas)=="MIC192T_COVv")
which(colnames(subexomas)=="MIC192T_FreqV")

MIC192T_list<-merge(subexomas, MIC192T_merge2, by=c("chr", "pos"))[,c(1,2,3,27,28,42,46,56,57,63,67:70)]

#9) Escrevendo os resultados

write.table(file="somatic_mutations_new.txt", MIC192T_list, row.names=F, sep="\t", append=T)

###########################################################################################################

#Trabalhando apenas com o paciente M27

#1) Filtrar SNVs do sangue com freq variante >= 20%  )

M27S<-subset(subexomas, subexomas$M27S_FreqV >= 20 )[,c(1,2,46)]

#3) Filtrar SNVs do tumor com freq variante >= 20%  )

M27T<-subset(subexomas, subexomas$M27T_FreqV >= 20 )[,c(1,2,3,46,59)]

#5) Fazendo o merge entre os SNVs de sangue e tumor 

M27T_merge<-merge(M27S, M27T, by=c("chr", "pos"), all.y=T)

colnames(M27T_merge)<-c("chr","pos","dbSNPID.sangue","gene","dbSNPID.tumor","cosmic")

#6) #As linhas com NA na coluna dbSNPID.sangue são mutações somáticas

M27T_merge2<-M27T_merge[is.na(M27T_merge$dbSNPID.sangue),]

#8) Obtendo a lista de genes

which(colnames(subexomas)=="M27T_COVr")
which(colnames(subexomas)=="M27T_COVv")
which(colnames(subexomas)=="M27T_FreqV")


M27T_list<-merge(subexomas, M27T_merge2, by=c("chr", "pos"))[,c(1,2,3,19,20,38,46,56,57,63,67:70)]

#9) Escrevendo os resultados

write.table(file="somatic_mutations_new.txt", M27T_list, row.names=F, sep="\t", append=T)

######################################################################################################

#Trabalhando apenas com o paciente M56

#1) Filtrar SNVs do sangue com freq variante >= 20%  )

M56S<-subset(subexomas, subexomas$M56S_FreqV >= 20 )[,c(1,2,46)]

#3) Filtrar SNVs do tumor com freq variante >= 20%  )

M56T<-subset(subexomas, subexomas$M56T_FreqV >= 20 )[,c(1,2,3,46,59)]

#5) Fazendo o merge entre os SNVs de sangue e tumor 

M56T_merge<-merge(M56S, M56T, by=c("chr", "pos"), all.y=T)

colnames(M56T_merge)<-c("chr","pos","dbSNPID.sangue","gene","dbSNPID.tumor","cosmic")

#6) #As linhas com NA na coluna dbSNPID.sangue são mutações somáticas

M56T_merge2<-M56T_merge[is.na(M56T_merge$dbSNPID.sangue),]

#8) Obtendo a lista de genes

which(colnames(subexomas)=="M56T_COVr")
which(colnames(subexomas)=="M56T_COVv")
which(colnames(subexomas)=="M56T_FreqV")

M56T_list<-merge(subexomas, M56T_merge2, by=c("chr", "pos"))[,c(1,2,3,23,24,40,46,56,57,63,67:70)]

#9) Escrevendo os resultados

write.table(file="somatic_mutations_new.txt", M56T_list, row.names=F, sep="\t", append=T)

###############################################################################################

#Trabalhando apenas com o paciente M155

#1) Filtrar SNVs do sangue com freq variante >= 20%  )

M155S<-subset(subexomas, subexomas$M155S_FreqV >= 20 )[,c(1,2,46)]

#3) Filtrar SNVs do tumor com freq variante >= 20%  )

M155T<-subset(subexomas, subexomas$M155T_FreqV >= 20 )[,c(1,2,3,46,59)]

#5) Fazendo o merge entre os SNVs de sangue e tumor 

M155T_merge<-merge(M155S, M155T, by=c("chr", "pos"), all.y=T)

colnames(M155T_merge)<-c("chr","pos","dbSNPID.sangue","gene","dbSNPID.tumor","cosmic")

#6) #As linhas com NA na coluna dbSNPID.sangue são mutações somáticas

M155T_merge2<-M155T_merge[is.na(M155T_merge$dbSNPID.sangue),]

#8) Obtendo a lista de genes

which(colnames(subexomas)=="M155T_COVr")
which(colnames(subexomas)=="M155T_COVv")
which(colnames(subexomas)=="M155T_FreqV")

M155T_list<-merge(subexomas, M155T_merge2, by=c("chr", "pos"))[,c(1,2,3,11,12,34,46,56,57,63,67:70)]

#9) Escrevendo os resultados

write.table(file="somatic_mutations_new.txt", M155T_list, row.names=F, sep="\t", append=T)

################################### Barplots cobertura ####################################################
05/02/15

MIC204T<-read.table("MIC204T.Coverage_full.txt", sep="\t", header=T, skip=9)
MIC204S<-read.table("MIC204S.Coverage_full.txt", sep="\t", header=T, skip=9)

MIC192T<-read.table("MIC192T.Coverage_full.txt", sep="\t", header=T, skip=9)
MIC192S<-read.table("MIC192S.Coverage_full.txt", sep="\t", header=T, skip=9)

M168T<-read.table("M168T.Coverage_full.txt", sep="\t", header=T, skip=9)
M168S<-read.table("M168S.Coverage_full.txt", sep="\t", header=T, skip=9)

M155T<-read.table("M155T.Coverage_full.txt", sep="\t", header=T, skip=9)
M155S<-read.table("M155S.Coverage_full.txt", sep="\t", header=T, skip=9)

M56T<-read.table("M56T.Coverage_full.txt", sep="\t", header=T, skip=9)
M56S<-read.table("M56S.Coverage_full.txt", sep="\t", header=T, skip=9)

M27T<-read.table("M27T.Coverage_full.txt", sep="\t", header=T, skip=9)
M27S<-read.table("M27S.Coverage_full.txt", sep="\t", header=T, skip=9)

pdf(file="exoma_Andrea_morehalf.pdf")
par(mfrow=c(2,3))

#MIC204
barplot(matrix(
  c(
    length(unique(MIC204T[MIC204T$X......2 >= 50,1]))/1000,
    length(unique(MIC204S[MIC204S$X......2 >= 50,1]))/1000,
    length(unique(MIC204T[MIC204T$X......3 >= 50,1]))/1000,
    length(unique(MIC204S[MIC204S$X......3 >= 50,1]))/1000,
    length(unique(MIC204T[MIC204T$X......4 >= 50,1]))/1000,
    length(unique(MIC204S[MIC204S$X......4 >= 50,1]))/1000
    
  ),
  ncol=3        ), 
  beside=T, 
  ylim=c(0,70), 
  names=c(">10x", ">20x",">100x"),
  main=c("MIC204 - IDC"),
  sub="Only genes with horizontal coverage > 50%",
  col.main="blue",
  ylab="locus counts (in thousands)"
)
legend(6.5,69,legend=c("tumor","sangue"), fill=c("dark grey", "light grey"), cex=0.8)


#MIC192
barplot(matrix(
  c(
    length(unique(MIC192T[MIC192T$X......2 >= 50,1]))/1000,
    length(unique(MIC192S[MIC192S$X......2 >= 50,1]))/1000,
    length(unique(MIC192T[MIC192T$X......3 >= 50,1]))/1000,
    length(unique(MIC192S[MIC192S$X......3 >= 50,1]))/1000,
    length(unique(MIC192T[MIC192T$X......4 >= 50,1]))/1000,
    length(unique(MIC192S[MIC192S$X......4 >= 50,1]))/1000
    
  ),
  ncol=3        ), 
  beside=T, 
  ylim=c(0,100), 
  names=c(">10x",">20x",">100x"),
  main="MIC192 - IDC",
  sub="Only genes with horizontal coverage > 50%",
  col.main="blue",
  ylab="locus counts (in thousands)"
)

#M168

barplot(matrix(
  c(
    length(unique(M168T[M168T$X......2 >= 50,1]))/1000,
    length(unique(M168S[M168S$X......2 >= 50,1]))/1000,
    length(unique(M168T[M168T$X......3 >= 50,1]))/1000,
    length(unique(M168S[M168S$X......3 >= 50,1]))/1000,
    length(unique(M168T[M168T$X......4 >= 50,1]))/1000,
    length(unique(M168S[M168S$X......4 >= 50,1]))/1000
    
  ),
  ncol=3        ), 
  beside=T, 
  ylim=c(0,100), 
  names=c(">10x",">20x",">100x"),
  main="M168 - IDC",
  sub="Only genes with horizontal coverage > 50%",
  col.main="blue",
  ylab="locus counts (in thousands)"
)

#M155T

barplot(matrix(
  c(
    length(unique(M155T[M155T$X......2 >= 50,1]))/1000,
    length(unique(M155S[M155S$X......2 >= 50,1]))/1000,
    length(unique(M155T[M155T$X......3 >= 50,1]))/1000,
    length(unique(M155S[M155S$X......3 >= 50,1]))/1000,
    length(unique(M155T[M155T$X......4 >= 50,1]))/1000,
    length(unique(M155S[M155S$X......4 >= 50,1]))/1000
    
  ),
  ncol=3        ), 
  beside=T, 
  ylim=c(0,100), 
  names=c(">10x",">20",">100"),
  main="M155 - DCIS",
  sub="Only genes with horizontal coverage > 50%",
  ylab="locus counts (in thousands)"
)

#M56

barplot(matrix(
  c(
    length(unique(M56T[M56T$X......2 >= 50,1]))/1000,
    length(unique(M56S[M56S$X......2 >= 50,1]))/1000,
    length(unique(M56T[M56T$X......3 >= 50,1]))/1000,
    length(unique(M56S[M56S$X......3 >= 50,1]))/1000,
    length(unique(M56T[M56T$X......4 >= 50,1]))/1000,
    length(unique(M56S[M56S$X......4 >= 50,1]))/1000
    
  ),
  ncol=3        ), 
  beside=T, 
  ylim=c(0,100), 
  names=c("10x",">20x",">100x"),
  main="M56 - DCIS",
  sub="Only genes with horizontal coverage > 50%",
  ylab="locus counts (in thousands)"
  
)

#M27

barplot(matrix(
  c(
    length(unique(M27T[M27T$X......2 >= 50,1]))/1000,
    length(unique(M27S[M27S$X......2 >= 50,1]))/1000,
    length(unique(M27T[M27T$X......3 >= 50,1]))/1000,
    length(unique(M27S[M27S$X......3 >= 50,1]))/1000,
    length(unique(M27T[M27T$X......4 >= 50,1]))/1000,
    length(unique(M27S[M27S$X......4 >= 50,1]))/1000
    
  ),
  ncol=3        ), 
  beside=T, 
  ylim=c(0,100), 
  names=c(">10x",">20x",">100x"),
  main="M27 - DCIS",
  sub="Only genes with horizontal coverage > 50%",
  ylab="locus counts (in thousands)"

)


dev.off()

