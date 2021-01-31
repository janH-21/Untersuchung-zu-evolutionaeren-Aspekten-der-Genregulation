
#setwd("")
rm(list=ls())
require("seewave")
require("tools")


################################ Check Ergebnisdateien #############################
listdir <- dir(pattern="*.txt")

if("result.txt" %in% listdir){
    
  # Überspringe bereits bearbetiete Dateien
  results <- read.table("result.txt", header = T)
  listdir <- listdir[listdir != "result.txt"]
  nfiles <- length(listdir)
    
  for(element in listdir){
    
    head_elements <- strsplit(file_path_sans_ext(element),"_")
    el_frequency <- strsplit(head_elements[[1]][2],"f-")[[1]][2]
    el_d0 <- strsplit(head_elements[[1]][3],"d0-")[[1]][2]
    el_k1 <- strsplit(head_elements[[1]][4],"k1-")[[1]][2]
    #el_off <- strsplit(head_elements[[1]][5],"off-")[[1]][2]
    #el_bas <- strsplit(head_elements[[1]][5],"bas-")[[1]][2]
    
    for(iter in 1:nrow(results)){
      
      entry <- results[iter,]
      
      if(as.double(el_frequency) == entry$frequency 
         && as.double(el_d0) == entry$d0 
         && as.double(el_k1) == entry$k1){
         #&& as.double(el_off) == entry$off){
         #&& as.double(el_bas) == entry$bas){
        listdir <- listdir[listdir != element]
      }
      
    }
  }
  
  nfiles <- nfiles - length(listdir)
  print(paste(as.character(nfiles), "files skipped."))
  
} else {
  
  # Wenn noch keine result.txt Datei vorhanden
  header <- c("mab","mac","mad","energy","start_reading","stop_reading", "max_target","frequency","d0","k1")
  #header <- c("mab","mac","mad","energy","start_reading","stop_reading", "max_target","frequency","d0","k1","off")
  #header <- c("mab","mac","mad","energy","start_reading","stop_reading", "max_target","frequency","d0","k1","bas")
  write(header,"result.txt", ncolumns = 11, sep="\t", append=TRUE)
  
} 


############################ Berechne Energie und Phase ###########################
# Browse Dateiliste
iterTotal <- length(listdir)
iterCount <- 0

for(j in listdir){
  
  # f1 header:  t \t k0 \t mRNA \t ptf \t target
  f1 <- read.table(j, header=T) 
  maxi <- length(f1$k0)
  los <- maxi*0.8 # los muss je nach Frequenz geändert werden -> immer manuell checken
  
  # Berechne instantane Fequenz und Verstärkung
  k0 <- f1$k0[los:maxi] - mean(f1$k0[los:maxi])
  mRNA <- f1$mRNA[los:maxi] - mean(f1$mRNA[los:maxi])
  ptf <- f1$ptf[los:maxi] - mean(f1$ptf[los:maxi])
  target <- f1$target[los:maxi] - mean(f1$target[los:maxi])
  a <- ifreq(k0, 1000,plot=FALSE)
  b <- ifreq(mRNA, 1000,plot=FALSE)
  c <- ifreq(ptf, 1000,plot=FALSE)
  d <- ifreq(target, 1000,plot=FALSE)
  ab <- median((unwrap(a$p[,2])-unwrap(b$p[,2])))
  ac <- median((unwrap(a$p[,2])-unwrap(c$p[,2])))
  ad <- median((unwrap(a$p[,2])-unwrap(d$p[,2])))
  
  mab <- -(ab%%(2*pi))
  mac <- -(ac%%(2*pi))
  mad <- -(ad%%(2*pi))
  
  energy <- sd(f1$target[los:maxi])/sd(f1$k0[los:maxi])
  
  # Wähle Parameter für Augabedatei
  head_elements <- strsplit(file_path_sans_ext(j),"_")
  frequency <- strsplit(head_elements[[1]][2],"f-")[[1]][2]
  d0 <- strsplit(head_elements[[1]][3],"d0-")[[1]][2]
  k1 <- strsplit(head_elements[[1]][4],"k1-")[[1]][2]
  #off <- strsplit(head_elements[[1]][5],"off-")[[1]][2]
  #bas <- strsplit(head_elements[[1]][5],"bas-")[[1]][2]
  
  # speichere Ergebnisse
  toFile<- c(mab,mac,mad,energy,los,maxi, max(f1$target[los:maxi]),frequency, d0, k1, off)
  write(toFile,"result.txt", ncolumns = 11, sep="\t", append=TRUE)
  
  iterCount <- iterCount + 1
  print(paste("File", as.character(iterCount), "of", as.character(iterTotal), "processed."))
  
}

