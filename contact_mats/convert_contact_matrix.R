library('fields')


# convert C to 10 year age-groups
convert_bins_5to10 <- function(C_by5, pop) {
  
  l <- dim(C_by5)[1]
  C_by10 <- matrix(nrow = 9, ncol = 9)
  
  col_count <- 1
  for (col in seq(1,l-1,by = 2)){
    row_count <- 1
    for (row in seq(1,l-1, by = 2)){
      p1 <- pop[row]
      p2 <- pop[row + 1]
      C_by10[row_count, col_count] <- ((C_by5[row, col] + C_by5[row, col + 1])*p1 +
                                       (C_by5[row + 1, col] + C_by5[row + 1, col + 1])*p2)/(p1+p2)
      row_count <- row_count + 1
    }
    col_count <- col_count + 1
  }
  C_by10
}

getContactMatrices <- function() {
  # Read matrix data from file
  data_dir <- getwd();
  household <- read.csv(paste(data_dir,"/household_18.csv",sep=""), header = FALSE);
  community <- read.csv(paste(data_dir,"/community_18.csv",sep=""), header = FALSE);
  school <- read.csv(paste(data_dir,"/school_18.csv",sep=""), header = FALSE);
  work <- read.csv(paste(data_dir,"/work_18.csv",sep=""), header = FALSE);
  return(list("household"=household,"work"=work,"school"=school,"community"=community));
}


pop_5y <- c(0.1010,0.0904,0.0857,0.0757,0.0739,0.0712,0.0729,0.0698,0.0622,
            0.0578,0.0527,0.0490,0.0410,0.0291,0.0245,0.0182,0.0109,0.0140)

# age demographics by 10 year age bin 
pop_10y     <- c(0.1914,0.1614,0.1451,0.1427,0.1200,0.1017,0.0701,0.0427,0.0249)
pop_10y_new <- c(0.1947,0.1649,0.1394,0.1287,0.1177,0.0919,0.0795,0.0532,0.0300)

nag5 <- 18 # number of age groups
ag_labs5 <- c("0 to 4" , "5 to 9", "10 to 14" ,"15 to 19", "20 to 24", "25 to 29", 
              "30 to 34", "35 to 39", "40 to 44", "45 to 49", "50 to 54", "55 to 59",
              "60 to 64", "65 to 69", "70 to 74", "75-79", "80 to 84", "85+")

nag10 <- 9
ag_labs10 <- c("0 to 9", "10 to 19", "20 to 29", "30 to 39", "40 to 49", "50 to 59", "60 to 69", "70 to 79" ,"80+")


Cmats <- getContactMatrices()

settings <- c('household','work','school','community')

for(i in 1:length(settings)) {
  
  setting <- settings[i]
  
  C_by5 <- as.matrix(Cmats[[setting]])
  rownames(C_by5) <- as.list(ag_labs5)
  colnames(C_by5) <- as.list(ag_labs5)
  
  x11()
  image.plot(C_by5, axes=F, useRaster=T, xlab='Age of Individual',ylab='Age of Contact')
  mtext(text=ag_labs5, side=2, line=0.3, at=seq(0,1,length.out=nag5), las=1, cex=0.8)
  mtext(text=ag_labs5, side=1, line=0.3, at=seq(0,1,length.out=nag5), las=2, cex=0.8)
  
  C_by10 <- convert_bins_5to10(C_by5,pop_5y)
  
  C_by10 <- sapply(1:nag10, function(k) C_by10[,k]*pop_10y_new[k]/pop_10y[k])
  
  rownames(C_by10) <- as.list(ag_labs10)
  colnames(C_by10) <- as.list(ag_labs10)

  x11()
  image.plot(C_by10, axes=F, useRaster=T, xlab='Age of Individual',ylab='Age of Contact',main=settings[i])
  mtext(text=ag_labs10, side=2, line=0.3, at=seq(0,1,length.out=nag10), las=1, cex=0.8)
  mtext(text=ag_labs10, side=1, line=0.3, at=seq(0,1,length.out=nag10), las=2, cex=0.8)
  
  saveRDS(C_by10, paste0("mat_", setting, ".RDS"))
}
