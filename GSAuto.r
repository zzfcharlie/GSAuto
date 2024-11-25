# library packages
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(illuminaio))








GSAuto <- function(idat_dir,platform){
  # Construct the file path based on the platform
  platform_paths <- c("27k", "450k", "850k", "935k", "MSA")
  if (platform %in% platform_paths) {
    base_path <- paste0("data/", platform, "/")
    controlIdx <- readRDS(paste0(base_path, "negControl.rds"))
    manifest <- readRDS(paste0(base_path, "manifest.rds"))
  } else {
    stop("Invalid platform specified.")
  }
  
  
  
  ### Sample information
  files <- list.files(idat_dir)
  files <- files[grepl(".idat", files)]
  files <- gsub("_Grn.idat", '', files)
  files <- gsub("_Red.idat", '', files)
  files <- unique(files)
  file <- paste(idat_dir,files,sep = "/")
  
  cat('loading', length(file), 'idat from', idat_dir, '\n')
  
  ####Initialization####
  detect_rate <- c()
  qc_res <- c()
  beta_matrix <- matrix(rep(0),nrow(manifest),length(files))
  pval_matrix <- matrix(rep(0),nrow(manifest),length(files))
  U_matrix <- matrix(rep(0),nrow(manifest),length(files))
  M_matrix <- matrix(rep(0),nrow(manifest),length(files))
  Intensity_matrix <- matrix(rep(0),nrow(manifest),length(files))

  
  cat('\nStart processing...\n')
  
  
  for(i in 1:length(files)){
    # Read IDAT data
    grn.name <- paste0(file[[i]],"_Grn.idat")
    red.name <- paste0(file[[i]],"_Red.idat")

    ida.grn <- suppressWarnings(illuminaio::readIDAT(grn.name))
    ida.red <- suppressWarnings(illuminaio::readIDAT(red.name))
    #Mean intensity (cy3 represents green intensity, cy5 represents red intenstiy)
    ida.Data <- cbind(cy3=ida.grn$Quants[,"Mean"],
                      cy5=ida.red$Quants[,"Mean"])
    colnames(ida.Data) <- c('G', 'R') # 1105209 × 2


    ### Classify probes based on different Probe Design

    ## Type I green channel
    ## IordG <- manifest[((manifest$DESIGN=='I')&(manifest$col=='G')),]
    IordG <- manifest[(!is.na(manifest$col))&(manifest$col=='G'),]
    ## 2-channel for green probes' U allele
    IuG2ch <- ida.Data[match(IordG$U, rownames(ida.Data)),]
    ## 2-channel for green probes' M allele
    ImG2ch <- ida.Data[match(IordG$M, rownames(ida.Data)),]
    IG.sset <- as.matrix(data.frame(M = ImG2ch[,'G'], U = IuG2ch[,'G'],
                                    row.names = IordG$Probe_ID))
    IG.sset[rowSums(is.na(IG.sset)) > 0] <- 0




    ## Type I red channel
    IordR <- manifest[(!is.na(manifest$col))&(manifest$col=='R'),]
    ## 2-channel for red probes' U allele
    IuR2ch <- ida.Data[match(IordR$U, rownames(ida.Data)),]
    ## 2-channel for red probes' M allele
    ImR2ch <- ida.Data[match(IordR$M, rownames(ida.Data)),]
    IR.sset <- as.matrix(data.frame(M = ImR2ch[,'R'], U = IuR2ch[,'R'],
                                    row.names = IordR$Probe_ID))
    IR.sset[rowSums(is.na(IR.sset)) > 0] <- 0

    ## Type II
    IIord <- manifest[is.na(manifest$col),]
    signal.II <- ida.Data[match(IIord$U, rownames(ida.Data)),c('G','R')]
    colnames(signal.II) <- c('M', 'U')
    rownames(signal.II) <- IIord$Probe_ID
    II.sset <- signal.II
    sset <- rbind(IR.sset,IG.sset,II.sset)


    M_matrix[,i] <- as.vector(sset[,1])
    U_matrix[,i] <- as.vector(sset[,2])





    ## Total intensities
    IR.sset.new <- as.data.frame(rowSums(IR.sset))
    IG.sset.new <- as.data.frame(rowSums(IG.sset))
    II.sset.new <- as.data.frame(rowSums(II.sset))


    colnames(IR.sset.new)[1] <- "Intensities"
    colnames(IG.sset.new)[1] <- "Intensities"
    colnames(II.sset.new)[1] <- "Intensities"

    new.sset <- rbind(IR.sset.new,IG.sset.new,II.sset.new)
    Intensity_matrix[,i] <- new.sset$Intensities



    probe_ind <- data.frame(index = rownames(ida.Data))


    negctls <- ida.Data[as.character(controlIdx),]

    negctls_sum <- rowSums(negctls)
    IR.pval <- 1 - (findInterval(IR.sset.new$Intensities, sort(2 * negctls[, 2]),left.open = TRUE))/nrow(negctls)
    IG.pval <- 1 - (findInterval(IG.sset.new$Intensities, sort(2 * negctls[, 1]),left.open = TRUE))/nrow(negctls)
    II.pval <- 1 - (findInterval(II.sset.new$Intensities, sort(negctls_sum),left.open = TRUE))/nrow(negctls)




    pval <- c(IR.pval,IG.pval,II.pval)


    ## Beta
    betas <- c(
      pmax(IR.sset[,'M'],0) / (abs(IR.sset[,'M'])+abs(IR.sset[,'U'])+100),
      pmax(IG.sset[,'M'],0) / (abs(IG.sset[,'M'])+abs(IG.sset[,'U'])+100),
      pmax(II.sset[,'M'],0) / (abs(II.sset[,'M'])+abs(II.sset[,'U'])+100))

    final.sset <- data.frame(betas = as.vector(betas), pval = pval)
    rownames(final.sset) <- rownames(new.sset)


    detect_rate[i] <- 1-sum(final.sset$pval > 0.05, na.rm = T)/nrow(manifest)
    qc_res[i] <- ifelse(detect_rate[i] > 0.96, 'yes', 'no')
    beta_matrix[,i] <- as.vector(betas)
    pval_matrix[,i] <- pval


  }

  detect_re <- data.frame(detection_rate = detect_rate,qc = qc_res)
  rownames(detect_re) <- files
  rownames(beta_matrix) <- rownames(final.sset)
  colnames(beta_matrix) <- files
  rownames(pval_matrix) <- rownames(final.sset)
  colnames(pval_matrix) <- files
  rownames(M_matrix) <- rownames(final.sset)
  colnames(M_matrix) <- files
  rownames(U_matrix) <- rownames(final.sset)
  colnames(U_matrix) <- files
  rownames(Intensity_matrix) <- rownames(final.sset)
  colnames(Intensity_matrix) <- files

  result_list <- list(
    detect_re = detect_re,
    beta_matrix = beta_matrix,
    pval_matrix = pval_matrix,
    M_matrix = M_matrix,
    U_matrix = U_matrix,
    Intensity_matrix = Intensity_matrix
  )
  return(result_list)
}














