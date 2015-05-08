mincover_vr<-function(v){
  #Função para somar coberturas variante  e referencia a cada amostra. Retorna a cobertura mínima entre as amostras
  k<-length(v)
  covers<-v[seq(1,k-1,2)]+v[seq(2,k,2)] #reference + variante
  #covers<-v[seq(1,k/2)]+v[seq((k/2)+1,k)] #minha alteração
  return(min(covers))
}
Variantes<-function(Variant_file,sufixos=c("T","S"),
                    mutations=c("MISSENSE","NONSENSE"),type=c("SNV"),location="somatic", sample="",
                    cover_limit=20,freq_lim_neg=20,freq_lim_pos=20,
                    barplot=TRUE, barplot_file="Exoma_barplot.pdf",
                    outfile=paste(sample, "Somatic_mutations.txt",sep="_")
                     )
                    #outfile="Somatic_mutations_new.txt"
                     {
  
  #Filtra SNVs com cobertura suficiente e seleciona SNVs somáticas
  #Filtrando SNVs
  exomas<-read.table(Variant_file, header=T, sep="\t")

  amostras <- grep(sample, colnames(exomas), value=T)     # nome da amostra escolhida

  amostras <- grep("_COV", amostras, value=T)            # selecionando so as colunas do coverage bruto (sangue e tumor)  , removendo frequencia da variante
  
  #verifica se o coverage minimo por paciente (sangue ou tumor)  é superior ao escolhido , retorna True , False
  cond_cover <- apply(exomas[,amostras],1,mincover_vr)>=cover_limit 
  #Filtrar mutações dos tipos contidos em "mutations"
  cond_class <- exomas$Variant_class2 %in% mutations
  #Filtrar mutações dos tipos contidos em "type"
  cond_type <- exomas$Type_SNV_or_INDEL %in% type
  
  #Gerando o subset do dataset original com as condições selecionadas
  subexomas <- exomas[cond_cover & cond_class & cond_type,]
  
  amostras <- sample
  #amostras <- grep("_COVr",amostras,value=TRUE)
  #amostras <- grep(paste(sufixos[1],"_",sep=""),amostras,value=TRUE)
  #amostras <- sapply(amostras,function(str){strsplit(str,sufixos[1],fixed=TRUE)[[1]][1]},USE.NAMES=FALSE)
  
  #Obtêm lista de amostras do arquivo 'Variant_file'
  inicio=1
  for (amostra in amostras){
    #A cada amostra, filtra mutações de interesse (somaticas ou germinativas) e as salva em arquivo output.
    col_num_FreqV_1<-which(colnames(subexomas)==paste(amostra,sufixos[1],"_FreqV",sep=""))
    col_num_FreqV_2<-which(colnames(subexomas)==paste(amostra,sufixos[2],"_FreqV",sep=""))
    col_num_COVr <- which(colnames(subexomas)==paste(amostra,sufixos[1],"_COVr",sep=""))
    col_num_COVv <- which(colnames(subexomas)==paste(amostra,sufixos[1],"_COVv",sep=""))
    col_num_COVr_sangue <- which(colnames(subexomas)==paste(amostra,sufixos[2],"_COVr",sep="")) #  alteração
    col_num_COVv_sangue <- which(colnames(subexomas)==paste(amostra,sufixos[2],"_COVv",sep="")) # alteração
    cond_sufixo1 <- subexomas[,col_num_FreqV_1] >= freq_lim_pos #presente no tumor
    cond_sufixo1[is.na(cond_sufixo1)]<-FALSE
    if (location=="somatic"){
      cond_sufixo2 <- subexomas[,col_num_FreqV_2] < freq_lim_neg #ausente no sangue
    }else{
      cond_sufixo2 <- subexomas[,col_num_FreqV_2] >= freq_lim_pos #presente no sangue
    }
    cond_sufixo2[is.na(cond_sufixo2)]<-FALSE
    This_mutations<-subexomas[cond_sufixo1 & cond_sufixo2,]
    #De acordo com 'location', filtra mutações somáticas, ausentes no sangue, ou germinativas, presentes em sangue e tumor.
    This_data<-data.frame(Pac=rep(amostra,NROW(This_mutations)),This_mutations[,c(1:3,col_num_COVr_sangue,col_num_COVv_sangue,col_num_COVr,col_num_COVv,col_num_FreqV_1,45:71)])
    colnames(This_data)<-c("Sample",colnames(subexomas)[1:3],
                           paste(sufixos[2],"COVr",sep="_"), paste(sufixos[2],"COVv",sep="_"),
                           paste(sufixos[1],"COVr",sep="_"),paste(sufixos[1],"COVv",sep="_"),
                           paste(sufixos[1],"FreqV",sep="_"),
                           colnames(subexomas)[c(45:71)])
    if(inicio==1){
      Final_data<-This_data
      inicio<-0
    }else{
      Final_data<-rbind(Final_data,This_data)
    }
  }
  write.table(file=outfile, Final_data, row.names=FALSE, col.names=TRUE, sep="\t", append=FALSE)
  if(barplot){
    #Se opção 'barplot' for TRUE, gera pdf com gráficos.
    nsamp<-length(amostras)
    nr<-floor(nsamp^0.5)
    nc<-ceiling(nsamp/nr)
    pdf(file=barplot_file)
    par(mfrow=c(nr,nc))
    for (amostra in amostras){
      files<-paste(amostra,sufixos,".Coverage_full.txt",sep="")
      F1<-read.table(files[1],sep="\t", header=T, skip=9)
      F2<-read.table(files[2],sep="\t", header=T, skip=9)
      Genecounts<-matrix(c(sum(F1$X......2 >= 50)/1000,sum(F2$X......2 >= 50)/1000,
                           sum(F1$X......3 >= 50)/1000,sum(F2$X......3 >= 50)/1000,
                           sum(F1$X......4 >= 50)/1000,sum(F2$X......4 >= 50)/1000),2,3)
      barplot(Genecounts, beside=T, ylim=c(0,70), names=c(">10x", ">20x",">100x"),
              main=c(amostra), sub="Only genes with horizontal coverage > 50%",
              col.main="blue", ylab="loci counts (in thousands)")
    if (amostra=="M155")   {legend(6.5,69,legend=c("tumor","sangue"), fill=c("dark grey", "light grey"), cex=0.8)}
      
    }
    dev.off()
  }
  #output: dataframe com mutações filtradas.
  return(Final_data)
}
