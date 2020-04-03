setwd("/home/albert/Bioinfo/MachineLearning/TareasFinales")
###Carga librerias####
library(CoExpNets)
library(CoExpROSMAP)
library(WGCNA)
source("EnrichmentAnnotationFunctions.R")#Carga las funciones
CoExpROSMAP::initDb()
#########DEJAR EN TRUE###############
loadnetwork<-TRUE
###Consulta de datos disponibles####
for(data.family in getNetworkCategories()){
  cat(paste0("Family of gene expression dataset ",data.family," available\n"))#Familia de datos, solo hay una
  tissues = getAvailableNetworks(data.family)
  for(tissue in tissues){#Datos en concreto. ad="Alzheimer Disease"
    cat(paste0("Gene expression profiling dataset for tissue ",tissue,"\n is ",
               getExprDataFromTissue(which.one=data.family,tissue=tissue,only.file = T)," available\n"))
  }
}
###Obtencion de covariables####
nets<-getAvailableNetworks("CoExpROSMAP")
#Visualizacion de la distribucion del atributo braaksc 
par(mfrow=c(2,2))
devnull<-lapply(nets, function(x){
  covs<-CoExpROSMAP::getCovariates(tissue=x, which.one = "CoExpROSMAP")
  barplot(table(covs$braaksc),main=paste("Braak stage in", x, sep = " "))
  return(covs)
})
par(mfrow=c(1,1))
###Extrae los datos de expresion con los que se va a trabajar####
expr_data<-getExprDataFromTissue(tissue = "allsamples", which.one = "CoExpROSMAP")
dim(expr_data)
###Extrae los nombres de muestras a retener: Aquellas con braaksc = 5 y 6
metadata<-CoExpROSMAP::getCovariates(tissue="allsamples", which.one="CoExpROSMAP")
rowstoinclude<-rownames(metadata)[c(which(metadata$braaksc==5), which(metadata$braaksc==6))]
###Filtra el conjunto de datos
expr_data<-expr_data[rownames(expr_data) %in% rowstoinclude,]
###Crea la red####
if(loadnetwork==FALSE){#Evita ejecutar otra vez si ya se ha creado
  system("mkdir CoexpressionResults")
  expr_network<-CoExpNets::getDownstreamNetwork(tissue = "MyRosMap",
                                              n.iterations = 20,
                                              net.type = "signed",
                                              debug = FALSE,
                                              expr.data = expr_data,
                                              job.path = "./CoexpressionResults")
}else{
  expr_network<-readRDS("./CoexpressionResults/netMyRosMap.13.it.20.rds")
}
networkname<-"./CoexpressionResults/netMyRosMap.13.it.20.rds"
#Algunos genes tienen un ID x.algo. Formatea los ID a x sin el .algo
names(expr_network$moduleColors)<-gsub(pattern = "\\.\\d+", replacement = "", 
                                       perl = TRUE, names(expr_network$moduleColors))
###Crea anotacion de la red####
#Las anotaciones GO y de tipo celular que genera getDownstreamNetwork() no son identicas a las que se generan aquí
go_data<-CoExpNets::getGProfilerOnNet(net.file=networkname, exclude.iea=FALSE, out.file=paste(
  networkname, "gprofII.csv", sep = "_"))
##Se puede trabajar directamente con go_data, es un data frame
CoExpNets::initDb()
cells_data<-CoExpNets::cellTypeByModule(return.processed = FALSE,
                                        tissue="MyRosMap",
                                        which.one = "new",
                                        plot.file=paste(networkname,".celltype.pdf", sep = "_"),
                                        net.in=networkname,
                                        legend="ROS/MAP cell type signals")

cells_data<-as.data.frame(cells_data)#Idem con cells_data
write.csv(cells_data, "./CoexpressionResults/CelltypeCorrelation.csv", quote = FALSE)
###Relacion con covariables####
#Estan ya guardadas en metadata
#Ordenalos de la misma forma
filt_metadata<-metadata[rowstoinclude, ]
m<-match(rownames(expr_data), rownames(filt_metadata))
filt_metadata<-filt_metadata[m,]
stopifnot(identical(rownames(expr_data),rownames(filt_metadata)))#Asegurate
covarCorr_data<-CoExpNets::corWithCatTraits(which.one="new", tissue=networkname,
                                            covs=filt_metadata,
                                            covlist = colnames(metadata))
write.csv(covarCorr_data, file = "./CoexpressionResults/CovariableCorrelation.csv", quote = FALSE)



###Tamaño de los modulos####
CoExpNets::plotModSizes(tissue=networkname,which.one="new")
###Clustering jerarquico de los modulos####
plotEGClustering(which.one="new",tissue=networkname)
###Correlacion entre modulos####
library(corrplot)
egs<-getNetworkEigengenes(which.one="new",tissue=networkname)
corrplot(cor(egs), type = "upper", tl.cex = 0.5, order='hclust')

##############################
####Seleccion de modulos####
#Basados en correlacion con covariables
selectModules<-as.character(sapply(1:nrow(covarCorr_data), function(x){
  if(length(which((covarCorr_data[x,]>1.5))) > 0){ #Coge los modulos que tengan corr > 1.5 en CUALQUIER covariable
    return(rownames(covarCorr_data)[x])
  }
}))
selectModules<-selectModules[selectModules!="NULL"]
#Basados en numero de genes
modulesize<-as.numeric(sapply(unique(expr_network$moduleColors), function(x){
  l<-length(which(expr_network$moduleColors == x))
  return(l)
}))
modsize_df<-data.frame(ModName=unique(expr_network$moduleColors), ModSize=modulesize)
modsize_df<-modsize_df[order(modsize_df$ModSize, decreasing = TRUE),]
selectModulesII<-as.character(modsize_df$ModName[1:4])
#############################
####Analisis de los modulos#### 
totalmod<-c(selectModules, selectModulesII)
annotations <- c("CC","MF","BP","keg","rea", "celltype")
#Bucle iterando por cada módulo todas las anotaciones posibles.
annot_list <- lapply(totalmod, function(x){
  aux_list <- lapply(annotations,function(y){
    if(y != "celltype"){
      return(functional_annotator(data=go_data, modulename = x, annot_type = y))
    }else{
      return(celltype_annotator(data = cells_data, modulename = x))
    }
  })
  names(aux_list) <- annotations
  return(aux_list)
  })
names(annot_list) <- totalmod
###Extraer genes con MM más alto de cada modulo###
mm<-getMM(tissue = networkname, which.one = "new", genes = colnames(expr_data),
          expr.data.file = expr_data) #Extrae el mm de todos los genes
TopGenes_list<-lapply(totalmod, function(x){
  df<-topmm_genextractor(mmdata = mm, modulename = x, cutoff = 5)
  return(df)
})
names(TopGenes_list)<-totalmod







