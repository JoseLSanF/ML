celltype_annotator<-function(data, modulename){
  w<-which(data[, modulename]<1)
  annotations<-data.frame(AnnotName=rownames(data)[w], AnnotScore=data[w,modulename])
  return(annotations)
}
functional_annotator <- function(data, modulename, annot_type) {
  module <-  data[data$query.number == modulename & data$domain == annot_type,]
  module <- module[order(module$p.value)[1:10],c("term.name","domain","p.value")]
  return(module[,c(1,3)])
}
topmm_genextractor<-function(mmdata, modulename, cutoff){
  genes<-mmdata[mmdata$module ,]
  genes<-genes[order(genes$mm, decreasing = TRUE),]
  return(genes[1:cutoff,])
}
