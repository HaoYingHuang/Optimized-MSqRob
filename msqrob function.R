convertMenuItem <- function(mi,tabName) {
  mi$children[[1]]$attribs['data-toggle']="tab"
  mi$children[[1]]$attribs['data-value'] = tabName
  if(length(mi$attribs$class)>0 && mi$attribs$class=="treeview"){
    mi$attribs$class=NULL
  }
  mi
}


my_topFeatures <- function(models, contrast, adjust.method, sort = TRUE, alpha = 1) {
  if (is(contrast, "matrix")) {
    if (ncol(contrast) > 1) {
      stop("Argument contrast is matrix with more than one column, only one contrast is allowed")
    }
  }
  logFC <- vapply(models,
                  getContrast,
                  numeric(1),
                  L = contrast
  )
  se <- sqrt(vapply(models,
                    varContrast,
                    numeric(1),
                    L = contrast
  ))
  df <- vapply(models, getDfPosterior, numeric(1))
  
  ## 去除缺失值 ##
  logFC <- logFC[!is.na(logFC)]
  se <- se[!is.na(se)]
  df <- df[!is.na(df)]
  
  t <- logFC / se
  pval <- pt(-abs(t), df) * 2
  
  if (adjust.method=="fdrtool"){
    adjPval <- fdrtool::fdrtool(pval,statistic = "pvalue",plot = F)
    adjPval <- adjPval$lfdr
  }else{
    adjPval <- p.adjust(pval, method = adjust.method)
  }
  
  out <- data.frame(logFC, se, df, t, pval, adjPval)
  if (sort) {
    if (alpha < 1) {
      ids <- adjPval < alpha
      out <- na.exclude(out[ids, ])
    }
    return(out[order(out$pval), ])
  } else {
    return(out)
  }
}

## filter function ##-------------------
my_filter_missval <- function(assay,thr,cond,condition){
  bin_data <- assay
  idx <- is.na(assay)
  bin_data[!idx] <- 1
  bin_data[idx] <- 0
  colnames(bin_data) <- cond
  for_keep_frame <- data.frame(label=colnames(assay),ID=cond,condition=condition)
  rownames(for_keep_frame) <- cond
  
  # Filter se on the maximum allowed number of
  # missing values per condition (defined by thr)
  keep <- bin_data %>%data.frame() %>%rownames_to_column(var = "Protein.IDs") 
  keep <- gather(keep,ID, value, -Protein.IDs)%>%left_join(., for_keep_frame, by = "ID") 
  # keep <- group_by(keep,Protein.IDs, condition)%>%dplyr::summarize(miss_val = table(condition)[1]-sum(value)) 
  keep <- group_by(keep,Protein.IDs, condition) %>%summarize(miss_val = n() - sum(value))
  #keep <- filter(keep,miss_val <= 0)
  keep <- spread(keep,condition, miss_val)  ## 可能要改 ##
  filtered <- keep 
  for (i in 1:nrow(filtered)){
    if(thr%in%filtered[i,]){
      filtered$filter[i] <-"yes" 
    }else{
      filtered$filter[i] <-"no" 
    }
  }
  return(filtered)
} 
# my_filter_missval <- function(assay,thr,ecols,cond,condition) {
#   
#   # Make assay values binary (1 = valid value)
#   bin_data <- assay
#   idx <- is.na(assay)
#   bin_data[!idx] <- 1
#   bin_data[idx] <- 0
#   colnames(bin_data) <- cond
#   for_keep_frame <- data.frame(label=colnames(assay)[ecols],ID=cond,condition=condition)
#   rownames(for_keep_frame) <- cond
#   
#   # Filter se on the maximum allowed number of
#   # missing values per condition (defined by thr)
#   keep <- bin_data %>%data.frame() %>%rownames_to_column() 
#   keep <- gather(keep,ID, value, -rowname)%>%left_join(., for_keep_frame, by = "ID") 
#   keep <- group_by(keep,rowname, condition) %>%summarize(miss_val = n() - sum(value)) 
#   keep <- filter(keep,miss_val <= thr) %>%spread(condition, miss_val)  ## 可能要改 ##
#   
#   assay <- assay[keep$rowname, ]
#   return(assay)
# }

## logTransform normalize function ##-----------
logTransform_and_normalize <- function(objectRaw,method=c("quantiles","vsn"),type){
  if(type=="peptide"){
    Raw <- "peptideRaw"
    Log <- "peptideLog"
    Norm <- "peptideNorm"
  }else if(type=="protein"){
    Raw <- "proteinRaw"
    Log <- "proteinLog"
    Norm <- "protein"
  }
  
  if(method=="vsn"){
    objectRaw <- QFeatures::normalize(objectRaw, i = Raw, method = method, name = Norm)
    objectRaw <- QFeatures::logTransform(objectRaw, base = 2, i = Raw, name = Log)
  }else{
    objectRaw <- QFeatures::logTransform(objectRaw, base = 2, i = Raw, name = Log)
    objectRaw <- QFeatures::normalize(objectRaw, i = Log, method = method, name = Norm)
  }
  return(objectRaw)
}

## ctrl_ecols function ##--------------------------
## ctrl_ecols== "ATG14" "con" ##
get_ctrl_ecols <- function(ecols,file){
  ctrl_ecols <- colnames(file)[ecols]%>%strsplit("\\.")   ##ctrl备选
  for (i in 1:length(ctrl_ecols)){
    ctrl_ecols[i] <- strsplit(ctrl_ecols[[i]][length(ctrl_ecols[[i]])],"_")[[1]][1]
  }
  ctrl_ecols <- unlist(ctrl_ecols)[!duplicated(unlist(ctrl_ecols))]
  return(ctrl_ecols)
}


## cond_vector ##---------------------
## cond_vector== "con_1"   "con_2"   "con_3"   "ATG14_1" "ATG14_2" "ATG14_3" ##
get_cond_vector <- function(corrected_file,ecols){
  cond_vector <- vector()
  for (i in 1:length(strsplit(colnames(corrected_file)[ecols],"\\."))){
    sp <- strsplit(colnames(corrected_file)[ecols],"\\.")
    cond_vector[[i]] <- sp[[i]][length(sp[[i]])]
  }
  return(cond_vector)
}

## condition_vector ##---------------------
## == con con con ATG14 ATG14 ATG14 named=peptideRaw ##
get_condition_vector <- function(cond_vector){
  condition_vector <- vector()
  for (i in 1:length(cond_vector)){
    condition_vector[i] <- strsplit(cond_vector,"_")[[i]][1]
  }
  names(condition_vector) <- c(rep("peptideRaw",length(condition_vector)))
  condition_vector <- as.factor(condition_vector)
  return(condition_vector)
}

## ecols function ##-------------------
get_ecols <- function(file,intensity_colmun){
  ecols <- grep(paste0("^",intensity_colmun,"."),colnames(file))
  return(ecols)
}

## corrected_file ##----------------------
## 按照Control将原文件的intensity列重新排列 ##
get_corrected_file <- function(Control,file,ecols){
  x <- grep(Control,colnames(file)[ecols])
  df <- file[ecols][,x]
  ctrl_column <- colnames(df)
  
  corrected_file <- relocate(file, match(ctrl_column,colnames(file)), .before = ecols[1])
  return(corrected_file)
}

## corrected_ecols ##-------------
## 得到corrected之后所需数据列
get_corrected_ecols <- function (corrected_file,Control,Experiment,intensity_colmun){
  if (intensity_colmun=="LFQ"){
    con <- paste0("LFQ.","intensity.",Control)
    exp <- paste0("LFQ.","intensity.",Experiment)
    corrected_ecols <- c(grep(con,colnames(corrected_file)),grep(exp,colnames(corrected_file)))
  }else{
    con <- paste0("Intensity.",Control)
    exp <- paste0("Intensity.",Experiment)
    corrected_ecols <- c(grep(con,colnames(corrected_file)),grep(exp,colnames(corrected_file)))
  }
  return(corrected_ecols)
}

## condition_matched_frame ##-----
## 得到替换条件的frame 
get_condition_matched_frame <- function(file,ecols){
  corrected_condition <- get_ctrl_ecols(file=file,ecols=ecols)
  final_condition <- vector()
  for (i in 2:length(corrected_condition)){
    final_condition[1] <- corrected_condition[1]
    final_condition[i] <- paste0(corrected_condition[i],"_vs_",corrected_condition[1])
  }
  condition_matched_frame <- data.frame(origin_condition=corrected_condition,
                                        matched_condition=LETTERS[1:length(corrected_condition)],
                                        replicate=rep(length(ecols)/length(corrected_condition),length(corrected_condition)),
                                        condition=final_condition)
  return(condition_matched_frame)
}

##得到替换条件frame中的condition
get_condition_from_frame <- function(condition_matched_frame){
  condition <- condition_matched_frame$condition[-1]
  return(condition)
}

## 得到符合condition的coldata
get_col_vector <- function(condition_matched_frame,type=c("peptide","protein")){
  col_vector <- vector()
  for(i in 1:nrow(condition_matched_frame)){
    col_vector <- sort(rep(condition_matched_frame$matched_condition,
                           condition_matched_frame$replicate[1]))
  }
  
  if(type=="peptide"){
    names(col_vector) <- rep("peptideRaw",length(col_vector))
  }else if(type=="protein"){
    names(col_vector) <- rep("proteinRaw",length(col_vector))
  }
  
  return(col_vector)
}

## reactive analysis function ##----------------------
my_test <- function(pe_imp,condition_matched_frame,file_type){
  if(file_type=="peptides"){
    pe <- suppressWarnings(msqrob(object = pe_imp,i = "protein",formula = ~condition))
    
    L <- my_makeContrast(condition_matched_frame)
    pe <- hypothesisTest(pe, i = "protein", L)
    
    return(pe)
  }else if(file_type=="protein groups"){
    pe_result <- suppressWarnings(msqrob(object = pe_imp,i = "protein",formula = ~condition))
    L <- my_makeContrast(condition_matched_frame)
    pe_result <- hypothesisTest(pe_result, i = "protein", L)
    
    return(pe_result)
  }
}

my_imp <- function(pe,filter_threshold,cond_vector,condition_vector,ID_colmun,
                   imputation_type,file_type){
  if(file_type=="peptides"){
    df_missval <- assay(pe[["protein"]])
    
    filtered <- my_filter_missval(assay(pe[["protein"]]),thr = filter_threshold,cond = cond_vector,condition = condition_vector)
    colnames(filtered)[1] <- ID_colmun
    rowData(pe[["protein"]]) <- left_join(as.data.frame(rowData(pe[["protein"]])),
                                          filtered,by = ID_colmun)
    pe[["protein"]] <- pe[["protein"]][rowData(pe[["protein"]])$filter== "yes", ]
    
    pe_imp <- pe
    if (imputation_type=="RF"){
      assay_protein <- assay(pe_imp[["protein"]])
      assay_protein <- missForest(assay_protein)
      assay(pe_imp[["protein"]]) <- assay_protein$ximp
    }else if(imputation_type=="No imputation"){
      assay_protein <- assay(pe_imp[["protein"]])
      assay(pe_imp[["protein"]]) <- assay_protein
    }else{
      assay_protein <- assay(pe_imp[["protein"]])
      assay(pe_imp[["protein"]]) <- impute_matrix(assay_protein, method = imputation_type)
    }
    
    pe_result <- pe_imp
    
    return(pe_result)
  }else if(file_type=="protein groups"){
    pe <- zeroIsNA(pe, "proteinRaw")
    filtered <- my_filter_missval(assay(pe[["proteinRaw"]]),thr = filter_threshold,cond = cond_vector,condition = condition_vector)
    
    rowData(pe[["proteinRaw"]]) <- left_join(as.data.frame(rowData(pe[["proteinRaw"]])),
                                             filtered,by ="Protein.IDs")
    try(pe[["proteinRaw"]] <- pe[["proteinRaw"]][rowData(pe[["proteinRaw"]])$Reverse != "+", ], silent = T)
    try(pe[["proteinRaw"]] <- pe[["proteinRaw"]][rowData(pe[["proteinRaw"]])$Potential.contaminant != "+", ], silent = T)
    pe[["proteinRaw"]] <- pe[["proteinRaw"]][rowData(pe[["proteinRaw"]])$filter== "yes", ]
    
    pe_log <- pe
    pe_log <- logTransform_and_normalize(pe_log,method = "quantiles",type = "protein")
    
    df_missval <- assay(pe_log[["protein"]])
    ## imputation ##
    if (imputation_type=="RF"){
      assay_protein <- missForest(df_missval)
      # assay(test[["protein"]]) <- reactive({assay_protein$ximp})
      assay(pe_log[["protein"]]) <- assay_protein$ximp
    }else if(imputation_type=="No imputation"){
      assay(pe_log[["protein"]]) <- df_missval
    }else{
      assay(pe_log[["protein"]]) <- impute_matrix(df_missval, method = imputation_type)
    }
    
    return(pe_log)
  }
}

read2agg <- function(corrected_file,ecols,condition_matched_frame,ID_colmun,
                     Genename_colmun,file_type){
  if(file_type=="peptides"){
    pe <- readQFeatures(table = corrected_file, fnames = 1, ecol = ecols,name = "peptideRaw", sep = "\t")
    
    col_vector <- get_col_vector(condition_matched_frame,type = "peptide")
    colData(pe)$condition <- as.factor(col_vector)
    
    rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)##
    pe <- zeroIsNA(pe, "peptideRaw")
    #filter
    pe[["peptideRaw"]] <- pe[["peptideRaw"]][rowData(pe[["peptideRaw"]])$Proteins%in% smallestUniqueGroups(rowData(pe[["peptideRaw"]])$Proteins), ]
    pe[["peptideRaw"]] <- pe[["peptideRaw"]][rowData(pe[["peptideRaw"]])$Reverse != "+", ]
    pe[["peptideRaw"]] <- pe[["peptideRaw"]][rowData(pe[["peptideRaw"]])$Potential.contaminant != "+", ]
    pe[["peptideRaw"]] <- pe[["peptideRaw"]][rowData(pe[["peptideRaw"]])$nNonZero >= 1, ]## 要改 ##
    
    pe_log <- pe
    pe_log <- logTransform_and_normalize(pe_log,method = "quantiles",type = "peptide")
    pe_log <- suppressWarnings(aggregateFeatures(pe_log,i = "peptideNorm",fcol = ID_colmun,na.rm = TRUE,name = "protein"))
    rowdata <- as.data.frame(rowData(pe_log[["protein"]]))
    
    if(!any(colnames(rowdata)==Genename_colmun)){
      for(i in 1:nrow(rowdata)){
        rowdata$Gene.names[i]<- corrected_file[corrected_file[,which(colnames(corrected_file)==ID_colmun)]%in%rowdata[,which(colnames(rowdata)==ID_colmun)][i],
                                               which(colnames(corrected_file)==Genename_colmun)]
      }
    }
    rowData(pe_log[["protein"]]) <- rowdata
    
    return(pe_log)
  }else if(file_type=="protein groups"){
    data_unique <- make_unique(corrected_file, Genename_colmun, ID_colmun, delim = ";")
    pro <-readQFeatures(table = data_unique,fnames = 1,ecol = ecols,name = "proteinRaw", sep="\t")
    
    # condition_matched_frame <- get_condition_matched_frame(file = data_unique,ecols = ecols)
    col_vector <- get_col_vector(condition_matched_frame,type = "protein")
    colData(pro)$condition <- as.factor(col_vector)
    
    return(pro)
  }
}

## get Filtered result ##------------------
## only remove the row contain NA

Add_Filtered_result <- function(object,LogFC,alpha,addedname,condition_matched_frame,ID_colmun){
  # Result_row <- as.data.frame(rowData(object[["protein"]]))
  # Result_row <- Result_row[!is.na(Result_row$logFC),]%>%rownames_to_column(var = "Proteins")
  Result_row <- as.data.frame(rowData(object[["protein"]]))
  Result_row <- Result_row[!is.na(Result_row$conditionB.logFC),]
  
  if(!any(colnames(as.data.frame(rowData(object[["protein"]])))=="Proteins")){
    Result_row <- rownames_to_column(Result_row,var = "Proteins")
  }else{
    Result_row <- Result_row
  }
  
  if(length(which(colnames(Result_row)=="Reverse"))==1){
    Result_row <- Result_row[!Result_row$Reverse=="+",]
    Result_row <- Result_row[,-grep("Reverse",colnames(Result_row))]
  }else{
    Result_row <- Result_row
  }
  
  Result_row <- Result_row[,c(grep("Gene.names",colnames(Result_row)),grep("Leading.razor.protein",colnames(Result_row)),
                              grep("Protein.group.IDs",colnames(Result_row)),grep("Proteins",colnames(Result_row)),
                              grep("Majority.protein.IDs",colnames(Result_row)),grep("Protein.IDs",colnames(Result_row)),
                              grep("logFC",colnames(Result_row)),grep("pval",colnames(Result_row)),grep("adjPval",colnames(Result_row)),
                              grep("se$",colnames(Result_row)),grep("df$",colnames(Result_row)),grep("t$",colnames(Result_row)))]
  sig_frame <- list()
  for (i in 1:(length(condition_matched_frame$condition)-1)){
    sig_frame[[i]] <- 1:nrow(Result_row)
  }
  sig_frame <- as.data.frame(sig_frame)
  colnames(sig_frame) <- paste0("condition",condition_matched_frame$matched_condition[2:nrow(condition_matched_frame)]
                                ,".significant")
  
  for (j in 1:ncol(sig_frame)){
    condition <- strsplit(colnames(sig_frame),"\\.")[[j]][1]
    for (i in 1:nrow(Result_row)){
      if(abs(Result_row[,which(colnames(Result_row)==paste0(condition,".logFC"))][i])>= LogFC
         & Result_row[,which(colnames(Result_row)==paste0(condition,".adjPval"))][i]<= alpha){
        sig_frame[i,j] <- "TRUE"
      }else{
        sig_frame[i,j] <- "FALSE"
      }
    }
  }
  
  Result_row <- cbind(Result_row,sig_frame)
  
  Result_assay <- as.data.frame(assay(object[["protein"]]))
  # Result_assay <- Result_assay[rownames(Result_assay)%in%Result_row$Proteins,]
  Result_assay <- Result_assay[rownames(Result_assay)%in%(Result_row[,which(colnames(Result_row)==ID_colmun)]),]
  se <- SummarizedExperiment(assays = Result_assay,
                             rowData = Result_row,
                             colData = colData(object))
  # se <- se[rowData(se)$significant >= "true", ]
  object <- addAssay(object,se,name = addedname)
  
  
  return(object)
}
# 使 用 例 pe <- Add_Filtered_result(pe,type = "peptides",LogFC = 1,alpha = 0.05,addedname = "try")


## get msqrob option table ##-------------------
get_msqrob_table <- function(object){
  temp <- as.data.frame(rowData(object[["FilteredResult"]]))%>%filter(significant=="TRUE")
  temp <- temp[,-c(grep("Leading.razor.protein",colnames(temp)),
                   grep("Taxonomy.IDs",colnames(temp)),
                   grep("Reverse",colnames(temp)),
                   grep("Potential.contaminant",colnames(temp)),
                   grep("msqrobModels",colnames(temp)))]
  rownames(temp) <- 1:nrow(temp)
  
  return(temp)
}
# get_msqrob_table(object = pe())

## output result ##
## add SummarizedExperiment for result ##
output_result <- function(object,LogFC,Padj){
  Norm_assay <- assay(object[["proteinNorm"]])
  Result_row <- rowData(object[["proteinNorm"]])$conditionB
  Result_row <- dplyr::mutate(as.data.frame(Result_row),significant=logFC)
  for (i in 1:nrow(Result_row)){
    if(abs(Result_row$logFC[i])>= LogFC & Result_row$adjPval[i]<= Padj){
      Result_row$significant[i] <- "true"
    }else{
      Result_row$significant[i] <- "false"
    }
  }
  se <- SummarizedExperiment(assays = Norm_assay,
                             rowData = Result_row,
                             colData = colData(object))
  se <- se[rowData(se)$significant >= "true", ]
  object <- addAssay(object,se,name = "FilteredResult")
  
  return(object)
}

## result from object
result_from_object <- function(object,condition_matched_frame){
  temp <- as.data.frame(rowData(object[["FilteredResult"]]))
  temp <- temp[,c(grep("Gene.names",colnames(temp)),grep("Proteins",colnames(temp)),
                  grep("protein",colnames(temp)),grep("Leading.razor.protein",colnames(temp)),
                  grep("logFC",colnames(temp)),grep("pval",colnames(temp)),
                  grep("adjPval",colnames(temp)),grep("significant",colnames(temp)))]
  
  colname <- colnames(temp)
  colname[grep("logFC",colname)] <- paste0(condition_matched_frame$condition[2:nrow(condition_matched_frame)],".logFC")
  colname[grep("pval",colname)] <- paste0(condition_matched_frame$condition[2:nrow(condition_matched_frame)],".pval")
  colname[grep("adjPval",colname)] <- paste0(condition_matched_frame$condition[2:nrow(condition_matched_frame)],".adjPval")
  colname[grep("significant",colname)] <- paste0(condition_matched_frame$condition[2:nrow(condition_matched_frame)],".significant")
  colnames(temp) <- colname
  if(!"Majority.protein.IDs"%in%colnames(temp)){
    temp <- rownames_to_column(temp,var = "Majority.protein.IDs")
  }
  
  for (i in 1:nrow(temp)){
    temp$Majority.protein.IDs[i] <- strsplit(temp$Majority.protein.IDs[i],";")[[1]][1]
    if(temp$Gene.names[i]==""){
      temp$Gene.names[i] <- temp$Majority.protein.IDs[i]
    }
    temp$Gene.names[i] <- strsplit(temp$Gene.names[i],";")[[1]][1]
  }
  
  return(temp)
  
}

## my makeContrast
my_makeContrast <- function(condition_matched_frame){
  temp1 <- vector()
  temp2 <- vector()
  for (i in 2:nrow(condition_matched_frame)){
    temp1[i-1] <- paste0("condition",condition_matched_frame$matched_condition[i],"=0")
    temp2[i-1] <- paste0("condition",condition_matched_frame$matched_condition[i])
  }
  
  L <- msqrob2::makeContrast(temp1,parameterNames = temp2)
  return(L)
}


## assign msqrob result##----------
assign_result <- function(name,object){
  assign(name,object,envir = get_msqrob_result_Env())
}

## plot MDS
plot_MDS <- function(object,type=c("peptide","protein"),top = top){
  if(type=="peptide"){
    assay <- assay(object[["peptideNorm"]])
  }else{
    assay <- assay(object[["protein"]])
  }
  assay_col <- vector()
  for(i in 1:length(colnames(assay))){
    assay_col[i] <- strsplit(colnames(assay)[i],"\\.")[[1]][length(strsplit(colnames(assay)[i],"\\.")[[1]])]
  }
  colnames(assay) <- assay_col
  limma::plotMDS(assay, col = as.numeric(colData(object)$condition),top=top)
}

## plot cor
plot_cor <- function (object, significant = TRUE, lower = -1, upper = 1, pal = "PRGn", condition_matched_frame,condition,
                      pal_rev = FALSE, indicate = NULL, font_size = 12, plot = TRUE,  add_values = FALSE, value_size = 10, digits = 2,
                      ...) 
{
  assertthat::assert_that(inherits(object, "QFeatures"), 
                          is.logical(significant), length(significant) == 1, is.numeric(lower), 
                          length(lower) == 1, is.numeric(upper), length(upper) == 
                            1, is.character(pal), length(pal) == 1, is.logical(pal_rev), 
                          length(pal_rev) == 1, is.numeric(font_size), length(font_size) == 
                            1, is.logical(plot), length(plot) == 1)
  if (!(lower >= -1 & upper >= -1 & lower <= 1 & upper <= 1)) {
    stop("'lower' and/or 'upper' arguments are not valid\n         Run plot_pca() with 'lower' and 'upper' between -1 and 1", 
         call. = FALSE)
  }
  pals <- RColorBrewer::brewer.pal.info %>% rownames_to_column() %>% 
    filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel", 
         " (qualitative panels also not allowed)\n", 
         "Run plot_pca() with one of the following 'pal' options: ", 
         paste(pals$rowname, collapse = "', '"), "'", 
         call. = FALSE)
  }
  if (any(is.na(assay(object[["FilteredResult"]])))) {
    stop("Missing values in '", deparse(substitute(object)), 
         "'. Use plot_dist() instead")
  }
  if (!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))
    col_data <- colData(object) %>% as.data.frame()
    columns <- colnames(col_data)
    if (any(!indicate %in% columns)) {
      stop("'", paste0(indicate, collapse = "' and/or '"), 
           "' column(s) is/are not present in ", deparse(substitute(object)), 
           ".\nValid columns are: '", paste(columns, 
                                            collapse = "', '"), "'.", call. = FALSE)
    }
    anno <- colData(object) %>% data.frame() %>% select(indicate)
    names <- colnames(anno)
    anno_col <- vector(mode = "list", length = length(names))
    names(anno_col) <- names
    for (i in names) {
      var = anno[[i]] %>% unique() %>% sort()
      if (length(var) == 1) 
        cols <- c("black")
      if (length(var) == 2) 
        cols <- c("orangered", "cornflowerblue")
      if (length(var) < 7 & length(var) > 2) 
        cols <- RColorBrewer::brewer.pal(length(var), 
                                         "Pastel1")
      if (length(var) >= 7) 
        cols <- RColorBrewer::brewer.pal(length(var), 
                                         "Set3")
      names(cols) <- var
      anno_col[[i]] <- cols
    }
    ha1 = HeatmapAnnotation(df = anno, col = anno_col, show_annotation_name = TRUE)
  }
  else {
    ha1 <- NULL
  }
  
  matched_cond <- paste0("condition",
                         condition_matched_frame$matched_condition[grep(condition,condition_matched_frame$condition)])
  if (significant) {
    if (!paste0(matched_cond,".significant") %in% colnames(rowData(object[["FilteredResult"]]))) { ## co↑ co↓
      stop("'significant' column is not present in '", 
           deparse(substitute(object)), "'\nRun add_rejections() to obtain the required column", 
           call. = FALSE)
    }
    temp <- object[as.data.frame(rowData(object[["FilteredResult"]]))[,paste0(matched_cond,".significant")]=="TRUE",]## co↑ co↓
  }
  cor_mat <- cor(assay(temp[["FilteredResult"]]))
  ht1 = Heatmap(cor_mat, col = circlize::colorRamp2(seq(lower, 
                                                        upper, ((upper - lower)/7)), if (pal_rev) {
                                                          rev(RColorBrewer::brewer.pal(8, pal))
                                                        }
                                                    else {
                                                      RColorBrewer::brewer.pal(8, pal)
                                                    }), heatmap_legend_param = list(color_bar = "continuous", 
                                                                                    legend_direction = "horizontal", legend_width = unit(5, 
                                                                                                                                         "cm"), title_position = "topcenter"), 
                name = "Pearson correlation", column_names_gp = gpar(fontsize = font_size), 
                row_names_gp = gpar(fontsize = font_size), top_annotation = ha1, cell_fun = if(add_values) {function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf(paste("%.", digits, "f", sep = ""), cor_mat[i, j]), x, y, gp = gpar(fontsize = value_size))}} else {NULL}, 
                ...)
  if (plot) {
    draw(ht1, heatmap_legend_side = "top")
  }
  else {
    df <- as.data.frame(cor_mat)
    return(df)
  }
}

## plot dis
plot_dist <- function (object, condition,condition_matched_frame ,significant = TRUE, pal = "YlOrRd", pal_rev = TRUE, 
                       indicate = NULL, font_size = 12, plot = TRUE, add_values = FALSE, value_size = 10, digits = 2, ...) 
{
  assertthat::assert_that(inherits(object, "QFeatures"), 
                          is.logical(significant), length(significant) == 1, is.character(pal), 
                          length(pal) == 1, is.logical(pal_rev), length(pal_rev) == 
                            1, is.numeric(font_size), length(font_size) == 1, 
                          is.logical(plot), length(plot) == 1)
  pals <- RColorBrewer::brewer.pal.info %>% rownames_to_column() %>% 
    filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel", 
         " (qualitative panels also not allowed)\n", 
         "Run plot_pca() with one of the following 'pal' options: ", 
         paste(pals$rowname, collapse = "', '"), "'", 
         call. = FALSE)
  }
  if (!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))
    col_data <- colData(object) %>% as.data.frame()
    columns <- colnames(col_data)
    if (any(!indicate %in% columns)) {
      stop("'", paste0(indicate, collapse = "' and/or '"), 
           "' column(s) is/are not present in ", deparse(substitute(object)), 
           ".\nValid columns are: '", paste(columns, 
                                            collapse = "', '"), "'.", call. = FALSE)
    }
    anno <- colData(object) %>% data.frame() %>% select(indicate)
    names <- colnames(anno)
    anno_col <- vector(mode = "list", length = length(names))
    names(anno_col) <- names
    for (i in names) {
      var = anno[[i]] %>% unique() %>% sort()
      if (length(var) == 1) 
        cols <- c("black")
      if (length(var) == 2) 
        cols <- c("orangered", "cornflowerblue")
      if (length(var) < 7 & length(var) > 2) 
        cols <- RColorBrewer::brewer.pal(length(var), 
                                         "Pastel1")
      if (length(var) >= 7) 
        cols <- RColorBrewer::brewer.pal(length(var), 
                                         "Set3")
      names(cols) <- var
      anno_col[[i]] <- cols
    }
    ha1 = HeatmapAnnotation(df = anno, col = anno_col, show_annotation_name = TRUE)
  }else {
    ha1 <- NULL
  }
  
  matched_cond <- paste0("condition",
                         condition_matched_frame$matched_condition[grep(condition,condition_matched_frame$condition)])
  
  if (significant) {
    if (!paste0(matched_cond,".significant") %in% colnames(rowData(object[["FilteredResult"]]))) { ## co↑ co↓
      stop("'significant' column is not present in '", 
           deparse(substitute(object)), "'\nRun add_rejections() to obtain the required column", 
           call. = FALSE)
    }
    temp <- object[as.data.frame(rowData(object[["FilteredResult"]]))[,paste0(matched_cond,".significant")]=="TRUE",]## co↑ co↓
  }
  dist_mat <- cluster::daisy(t(assay(temp[["FilteredResult"]])), metric = "gower") %>% 
    as.matrix()
  max <- max(dist_mat)
  ht1 = Heatmap(dist_mat, col = circlize::colorRamp2(seq(0, 
                                                         max, ((max)/7)), if (pal_rev) {
                                                           rev(RColorBrewer::brewer.pal(8, pal))
                                                         }
                                                     else {
                                                       RColorBrewer::brewer.pal(8, pal)
                                                     }), heatmap_legend_param = list(color_bar = "continuous", 
                                                                                     legend_direction = "horizontal", legend_width = unit(5, 
                                                                                                                                          "cm"), title_position = "topcenter"), 
                name = "Gower's distance", column_names_gp = gpar(fontsize = font_size), 
                row_names_gp = gpar(fontsize = font_size), top_annotation = ha1, cell_fun = if(add_values) {function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf(paste("%.", digits, "f", sep = ""), dist_mat[i, j]), x, y, gp = gpar(fontsize = value_size))}} else {NULL},
                ...)
  if (plot) {
    draw(ht1, heatmap_legend_side = "top")
  }
  else {
    df <- as.data.frame(dist_mat)
    return(df)
  }
}

## plot cv
plot_cvs<-function(object) {
  
  ## backtransform data
  untransformed_intensity<- 2^(assay(object[["FilteredResult"]]))
  exp_design <-as.data.frame(colData(object))
  exp_design <- tibble::rownames_to_column(exp_design,var = "ID")
  condition <- vector()
  for (i in 1:nrow(exp_design)){
    condition[i] <- strsplit(exp_design$ID,"\\.")[[i]][length(strsplit(exp_design$ID,"\\.")[[i]])]
    condition[i] <- strsplit(condition[i],split = "_")[[1]][1]
  }
  exp_design <- data.frame(ID=exp_design$ID,condition=condition)
  
  ### merge untransformed to exp design and calculate cvs
  
  cvs_group<- untransformed_intensity %>% data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::gather("ID", "Intensity", -rowname) %>%
    dplyr::left_join(.,data.frame(exp_design), by="ID") %>%
    dplyr::group_by(rowname,condition) %>%
    dplyr::summarise(cvs=coef_variation(Intensity)) %>%
    dplyr::group_by(condition)%>%
    dplyr::mutate(condition_median=median(cvs))
  
  p1 <-  ggplot(cvs_group, aes(cvs, color=condition, fill=condition)) +
    geom_histogram(alpha=.5, bins= 20, show.legend = FALSE) +
    facet_wrap(~condition) +
    geom_vline(aes(xintercept=condition_median, group=condition),color='grey40',
               linetype="dashed") +
    labs(title= 'Sample Coefficient of Variation', x="Coefficient of Variation", y="Count") +
    theme_DEP2() +
    theme(plot.title = element_text(hjust = 0.5,face = "bold")) 
  
  p1 +geom_text(aes(x=max(cvs_group$cvs)-0.6,
                    y=max(ggplot_build(p1)$data[[1]]$ymax*1.1), 
                    label=paste0("Median =",round(condition_median,2)*100,"%",by="")),
                show.legend = FALSE, size=4)
  
}


## my plot volcano ##---------------------
plot_volcano <- function (row_data,condition, condition_matched_frame,label_size = 3, add_names = TRUE, adjusted = FALSE, 
                          plot = TRUE, same_width = TRUE, my_breaks = FALSE, mybreaks = NULL) 
{
  if (is.integer(label_size)) 
    label_size <- as.numeric(label_size)
  
  # row_data <- get("frame",envir = get_msqrob_result_Env())
  if (any(!c("Gene.names", "Majority.protein.IDs") %in% colnames(row_data))) {               ## co↑co↓
    stop(paste0("'Gene.names' and/or 'Majority.protein.IDs' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun make_unique() to obtain required columns."), 
         call. = FALSE)
  }
  if (length(grep("logFC", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."), 
         call. = FALSE)
  }
  if (length(grep("significant", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_significant' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required columns."), 
         call. = FALSE)
  }
  
  
  if (length(grep(paste(condition, ".logFC", sep = ""), 
                  colnames(row_data))) == 0) {
    valid_cntrsts <- row_data %>% data.frame() %>% dplyr::select(ends_with("logFC")) %>% 
      colnames(.) %>% gsub(".logFC", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '", 
                                paste0(valid_cntrsts, collapse = "', '"), "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n", 
         valid_cntrsts_msg, call. = FALSE)
  }
  #when name have Disadvantages eg: when contrast = c("H2A_vs_Biotin","uH2A_vs_Biotin"), when set contrast = c("H2A_vs_Biotin"), can find two cols,and then error, the same problem with p_values and signif
  #diff <- grep(paste(contrast, "_diff", sep = ""), 
  #    colnames(row_data))
  diff <- match(paste(condition, ".logFC", sep = ""), 
                colnames(row_data))
  if (adjusted) {
    p_values <- match(paste(condition, ".adjPval", sep = ""), 
                      colnames(row_data))
  }
  else {
    p_values <- match(paste(condition, ".pval", sep = ""), 
                      colnames(row_data))
  }
  signif <- match(paste(condition, ".significant", sep = ""), 
                  colnames(row_data))
  
  df <- data.frame(x = row_data[, diff],y = -log10(row_data[,p_values]), 
                   significant = as.logical(row_data[, signif]), name = row_data$Gene.names,ID=row_data$Majority.protein.IDs) %>% 
    dplyr::filter(!is.na(significant)) %>% arrange(significant)
  for(i in 1:nrow(df)){
    if(df[i,4]==""){
      df[i,4] <- df[i,5]
    }
  }
  
  name1 <- gsub("_vs_.*", "", condition)
  name2 <- gsub(".*_vs_", "", condition)
  xlimit <- ceiling(max(c(abs(min(df$x)), abs(max(df$x)))))
  p <- ggplot(df, aes(x, y)) + geom_vline(xintercept = 0) + 
    geom_point(aes(col = significant)) + 
    geom_text(data = data.frame(), aes(x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1,0), 
                                       vjust = c(-1, -1), label = c(name1, name2), size = 5, fontface = "bold")) + 
    #scale_x_continuous(limits = c(-xlimit, xlimit), breaks = seq(-15, 15, by = 5)) +
    labs(title = condition, 
         x = expression(log[2] ~ "Fold change")) + theme_DEP1() + 
    theme(legend.position = "none") + scale_color_manual(values = c(`TRUE` = "black", 
                                                                    `FALSE` = "grey"))
  
  if (add_names) {
    p <- p + ggrepel::geom_text_repel(data = dplyr::filter(df, significant), 
                                      aes(label = name), size = label_size,
                                      box.padding = unit(0.1,"lines"), point.padding = unit(0.1, "lines"), 
                                      segment.size = 0.5,
                                      max.overlaps = 10) ##默认是10 需改大
  }
  if (adjusted) {
    p <- p + labs(y = expression(-log[10] ~ "Adjusted p-value"))
  }
  else {
    p <- p + labs(y = expression(-log[10] ~ "P-value"))
    
  }
  if(same_width){
    p <- p + scale_x_continuous(limits = c(-xlimit, xlimit))
  }
  if(same_width & my_breaks){
    p <- p + scale_x_continuous(limits = c(-xlimit, xlimit), breaks = mybreaks)
  }
  if (plot) {
    return(p)
  }
  else {
    df <- df %>% dplyr::select(name, x, y, significant) %>% arrange(desc(x))
    colnames(df)[c(1, 2, 3)] <- c("protein", "log2_fold_change", 
                                  "p_value_-log10")
    if (adjusted) {
      colnames(df)[3] <- "adjusted_p_value_-log10"
    }
    return(df)
  }
}
## 使 用 例 p <- my_plot_volcano(pe,"B1_vs_A1")


## my plot pca##------------
# object <- get("pe",envir = get_msqrob_result_Env())
# x = 1
# y = 2
# # indicate = c("condition","replicate")
# indicate = "condition"
# label = FALSE
# n = 500
# point_size = 4
# label_size = 3
# plot = TRUE
# if_square = FALSE
plot_pca <- function (object, x = 1, y = 2,condition_matched_frame,
                      label = FALSE, n, point_size = 4,label_size = 3, plot = TRUE, if_square = FALSE) 
{
  assay_data <- assay(object[["protein"]])
  if(any(is.na(assay_data))){
    stop("There is na in the data, PCA cannot be performed")
  }
  # assay_col <- vector()
  # for(i in 1:length(colnames(assay_data))){
  #   assay_col[i] <- strsplit(colnames(assay_data)[i],"\\.")[[1]][length(strsplit(colnames(assay_data)[i],"\\.")[[1]])]
  # }
  # colnames(assay_data) <- assay_col
  
  var <- apply(assay_data, 1, sd)                       
  
  df <- assay_data[order(var, decreasing = TRUE)[seq_len(n)],]
  pca <- prcomp(t(df), scale = FALSE)
  pca_df <- pca$x %>% data.frame() %>% rownames_to_column(var = "label") 
  experimental_design <- as.data.frame(colData(object))%>%rownames_to_column(var = "label")
  for(i in 1:nrow(experimental_design)){
    experimental_design$replicate[i] <-strsplit(experimental_design$label[i],"_")[[1]][length(strsplit(experimental_design$label[3],"_")[[1]])]
    for(j in 1:nrow(condition_matched_frame)){
      if(experimental_design$condition[i]==condition_matched_frame$matched_condition[j]){
        experimental_design$origin_condition[i] <- condition_matched_frame$origin_condition[j]
      }
    }
  }
  
  plot_df <- left_join(pca_df,experimental_design, by = "label")
  plot_df <- plot_df[,-which(colnames(plot_df)=="condition")]
  colnames(plot_df)[which(colnames(plot_df)=="origin_condition")] <- "condition"
  percent <- round(100 * pca$sdev^2/sum(pca$sdev^2), 1)
  for (feat in c("replicate","condition")) {
    plot_df[[feat]] <- as.factor(plot_df[[feat]])
  }
  
  indicate <- c("condition", 
                "replicate")
  point_size = 4
  
  rownames(plot_df) <- plot_df$label
  limit = unlist(plot_df[,c(paste0("PC", 1), paste0("PC", 2))])
  p <- ggplot(plot_df, aes(get(paste0("PC", 1)), get(paste0("PC",2)))) +
    labs(title = paste0("PCA plot - top ",as.character(n), " variable proteins"), x = paste0("PC",1, ": ", percent[1], "%"), y = paste0("PC",2, ": ", percent[2], "%")) + 
    coord_fixed() + geom_point(aes(col = .data[[indicate[1]]], 
                                   shape = .data[[indicate[2]]]), size = point_size) + 
    labs(col = indicate[1], shape = indicate[2])
  # + theme_DEP1()
  if(if_square){
    p <- p +
      scale_x_continuous(limits = range(limit)) +
      scale_y_continuous(limits = range(limit))
  }
  
  
  if (length(indicate) == 0) {
    p <- p + geom_point(size = point_size)
  }
  if (length(indicate) == 1) {
    p <- p + geom_point(aes(col = .data[[indicate[1]]]), 
                        size = point_size) + labs(col = indicate[1])
  }
  if (length(indicate) == 2) {
    p <- p + geom_point(aes(col = .data[[indicate[1]]], 
                            shape = .data[[indicate[2]]]), size = point_size) + 
      labs(col = indicate[1], shape = indicate[2])
  }
  if (length(indicate) == 3) {
    p <- p + geom_point(aes(col = .data[[indicate[1]]], 
                            shape = .data[[indicate[2]]]), size = point_size) + 
      facet_wrap(~.data[[indicate[3]]]) +
      labs(col = indicate[1], shape = indicate[2])
  }
  if (label) {
    p <- p + geom_text(aes(label = rowname), size = label_size)
  }
  if (plot) {
    return(p)
  }else {
    df <- plot_df %>% dplyr::select(rowname, paste0("PC", c(x, 
                                                            y)), match(indicate, colnames(pca_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }
}

## my plot heatmap ##---------------
# msqrob <- object
# # type = c("contrast", "centered")
# type="centered"
# manual = FALSE
# manual_name = NULL
# same_trend = FALSE
# kmeans = FALSE
# k = 6
# color = "RdBu"
# col_limit = 6
# indicate = NULL
# # clustering_distance = c("euclidean", "maximum", "manhattan", "canberra",
# #                         "binary", "minkowski", "pearson", "spearman", 
# #                         "kendall", "gower")
# clustering_distance = "euclidean"
# row_font_size = 6
# col_font_size = 10
# plot = TRUE
# if_mysplit = FALSE
# mysplit = NULL
# if_rowname_color = FALSE
plot_heatmap <- function (object,row_data,select,condition_matched_frame,type = c("contrast", "centered"), manual = FALSE, manual_name = NULL, same_trend = FALSE,
                          kmeans = FALSE, k = 6, color, col_limit, indicate = NULL,
                          clustering_distance = c("euclidean", "maximum", "manhattan", "canberra",
                                                  "binary", "minkowski", "pearson", "spearman", 
                                                  "kendall", "gower"), row_font_size = 6, col_font_size = 10,
                          plot = TRUE, if_mysplit = FALSE, mysplit = NULL, if_rowname_color = FALSE, ...) 
{
  library(ComplexHeatmap)
  library(tidyr)
  library(tibble)
  if (is.integer(k)) 
    k <- as.numeric(k)
  if (is.integer(col_limit)) 
    col_limit <- as.numeric(col_limit)
  if (is.integer(row_font_size)) 
    row_font_size <- as.numeric(row_font_size)
  if (is.integer(col_font_size)) 
    col_font_size <- as.numeric(col_font_size)
  assertthat::assert_that(inherits(object, "QFeatures"), 
                          is.character(type), is.logical(kmeans), is.numeric(k), 
                          length(k) == 1, is.numeric(col_limit), length(col_limit) == 
                            1, is.numeric(row_font_size), length(row_font_size) == 
                            1, is.numeric(col_font_size), length(col_font_size) == 
                            1, is.logical(plot), length(plot) == 1)
  
  # set rowname of heatmap color, which Peptides == 1 to purple, ==2 to blue
  # if(if_rowname_color) {
  #   rowData(object)$row_name_color = "black"
  #   rowData(object)$row_name_color[which(rowData(object)$Peptides == 1)] = "purple"
  #   rowData(object)$row_name_color[which(rowData(object)$Peptides == 2)] = "blue"
  # }
  
  type <- match.arg(type,c("contrast", "centered"))
  clustering_distance <- match.arg(clustering_distance,c("euclidean", "maximum", "manhattan", "canberra",
                                                         "binary", "minkowski", "pearson", "spearman", 
                                                         "kendall", "gower"))
  
  # row_data <- get("frame",envir = get_msqrob_result_Env())%>%table_select(select = select)
  assay <- as.data.frame(assay(object[["protein"]]))
  row_data <- table_select(row_data,select = select)
  df<- assay[rownames(assay)%in%row_data$Proteins,]
  df1 <- df[,grep(condition_matched_frame[1,1],colnames(df))]
  
  origin <- vector()
  for(i in 1:length(select)){
    origin[i] <- condition_matched_frame$origin_condition[grep(select[i],condition_matched_frame$condition)]
    colnames(df) <- gsub(origin[i],select[i],colnames(df))
  }
  
  df2 <- df[,grep("vs",colnames(df))]
  for(i in 1:length(select)){
    origin[i] <- condition_matched_frame$origin_condition[grep(select[i],condition_matched_frame$condition)]
    colnames(df2) <- gsub(select[i],origin[i],colnames(df2))
  }
  rm(origin)
  df <- cbind(df1,df2)
  
  df_mean<-apply(df,1,mean)
  df <- data.frame(df,mean=df_mean)
  df <- df[,1:(ncol(df)-1)] - df$mean
  df_colname <-vector() 
  for(i in 1:ncol(df)){
    df_colname[i] <- strsplit(colnames(df)[i],"\\.")[[1]][length(strsplit(colnames(df)[1],"\\.")[[1]])]
  }
  colnames(df) <- df_colname
  
  dfrowname <- vector()
  for(i in 1:nrow(df)){
    dfrowname[i] <- row_data$Gene.names[row_data$Proteins%in%rownames(df)[i]]
  }
  dfrowname <- make.names(dfrowname,unique = T)
  rownames(df) <- dfrowname
  
  # col_data <- rownames_to_column(as.data.frame(colData(object)),var = "label")
  # for (i in 1:nrow(col_data)){
  #   col_data$replicate[i] <- strsplit(col_data$label,"_")[[i]][2]
  # }
  # rownames(col_data) <- col_data$label
  
  # if (any(!c("label", "condition", "replicate") %in% 
  #         colnames(col_data))) {
  #   stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '", 
  #               deparse(substitute(object)), "'"), call. = FALSE)
  # }
  if (length(grep("logFC", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' columns are not present in '", 
                deparse(substitute(object)), "'.\nRun test_diff() to obtain the required columns."), 
         call. = FALSE)
  }
  # if (!"significant" %in% colnames(row_data)) {
  #   stop(paste0("'significant' column is not present in '", 
  #               deparse(substitute(object)), "'.\nRun add_rejections() to obtain the required column."), 
  #        call. = FALSE)
  # }
  if (!is.null(indicate) & type == "contrast") {
    warning("Heatmap annotation only applicable for type = 'centered'", 
            call. = FALSE)
  }
  
  if (!is.null(indicate) & type == "centered") {   ##有待探究##
    ha1 <- get_annotation(object, indicate)
  } else {
    ha1 <- NULL
  }
  #filtered <- res[row_data$significant, ]
  ## Filter for significant proteins only
  library(stringr)
  
  # filtered <- object[["FilteredResult"]]
  # filtered <- as.data.frame(rowData(filtered))
  # filtered <- filtered[rowData(filtered)$significant =="TRUE", ]
  # pe[["peptideLog"]] <- pe[["peptideLog"]][rowData(pe[["peptideLog"]])$nNonZero >= 2, ]
  
  # res <- get("frame",envir = get_msqrob_result_Env())
  # cols <- grep("significant", colnames(res))
  # names <- colnames(res)[cols]
  # names <- gsub(".significant", "", names)
  # selectizeInput("select",
  #                "Select direct comparisons",
  #                choices=names,
  #                multiple = TRUE)
  
  
  # # } else{
  if(any(is.na(assay))) {
    warning("Missing values in '", deparse(substitute(object)), "'. ",
            "Using clustering_distance = 'gower'",
            call. = FALSE)
    clustering_distance <- "gower"
    obs_NA <- TRUE
  } else {
    obs_NA <- FALSE
  }
  
  
  # if(type == "centered") {
  #   if(manual){
  #     rowData(filtered)$mean <- rowMeans(assay(filtered)[ , ind], na.rm = TRUE)
  #     df <- assay(filtered)[ , ind] - rowData(filtered, use.names = FALSE)$mean
  #   } else{
  #     rowData(filtered)$mean <- rowMeans(assay(filtered), na.rm = TRUE)
  #     df <- assay(filtered) - rowData(filtered, use.names = FALSE)$mean
  #   }
  # }
  # Get contrast fold changes ('contrast')   ## 这种type暂时没用 ## 
  # if(type == "contrast") {
  #   df <- rowData(filtered, use.names = FALSE) %>%
  #     data.frame() %>%
  #     column_to_rownames(var = "name") %>%
  #     dplyr::select(ends_with("_diff"))
  #   colnames(df) <-
  #     gsub("_diff", "", colnames(df)) %>%
  #     gsub("_vs_", " vs ", .)
  #   df <- as.matrix(df) 
  #   if(manual){
  #     i = gsub("_vs_", " vs ", manual_name)
  #     ii = as.data.frame(df[ , i])
  #     colnames(ii) = i
  #     df = as.matrix(ii)
  #   }
  # }
  
  # Facultative kmeans clustering
  if(kmeans & obs_NA) {
    warning("Cannot perform kmeans clustering with missing values",
            call. = FALSE)
    kmeans <- FALSE
  }
  if(kmeans & !obs_NA) {
    set.seed(1)
    df_kmeans <- kmeans(df, k)
    if(type == "centered") {
      # Order the k-means clusters according to the maximum fold change
      # in all samples averaged over the proteins in the cluster
      order <- data.frame(df) %>%
        cbind(., cluster = df_kmeans$cluster) %>%
        mutate(row = apply(.[, seq_len(ncol(.) - 1)], 1, function(x) max(x))) %>%
        group_by(cluster) %>%
        summarize(index = sum(row)/n()) %>%
        arrange(desc(index)) %>%
        pull(cluster) %>%
        match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
    if(type == "contrast") {
      # Order the k-means clusters according to their average fold change
      order <- data.frame(df) %>%
        cbind(df, cluster = df_kmeans$cluster) %>%
        gather(condition, diff, -cluster) %>%
        group_by(cluster) %>%
        summarize(row = mean(diff)) %>%
        arrange(desc(row)) %>%
        pull(cluster) %>%
        match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
  }
  
  
  if(clustering_distance == "gower") {
    clustering_distance <- function(x) {
      dist <- cluster::daisy(x, metric = "gower")
      dist[is.na(dist)] <- max(dist, na.rm = TRUE)
      return(dist)
    }
  }
  
  # Legend info
  legend <- ifelse(type == "contrast",
                   "log2 Fold change",
                   "log2 Centered intensity")
  
  if(if_mysplit) {
    cluster_row_slices = FALSE
    cluster_column_slices = FALSE
    split = if(kmeans) {factor(df_kmeans$cluster, levels = mysplit)} else {NULL}
  } else {
    cluster_row_slices = TRUE
    cluster_column_slices = TRUE
    split = if(kmeans) {df_kmeans$cluster} else {NULL}
  }
  
  #set rowname color
  row_name_color = if(if_rowname_color) {rowData(filtered)$row_name_color} else {"black"}
  
  # Heatmap
  df <- as.matrix(df)
  ht1 = ComplexHeatmap::Heatmap(df,
                                col = circlize::colorRamp2(
                                  seq(-col_limit, col_limit, (col_limit/5)),
                                  rev(RColorBrewer::brewer.pal(11, color))),
                                split = mysplit,
                                cluster_row_slices = cluster_row_slices,
                                cluster_column_slices = cluster_column_slices,
                                # cluster_rows = col_clust,
                                # cluster_columns = row_clust,
                                row_names_side = "left",
                                column_names_side = "top",
                                clustering_distance_rows = clustering_distance,
                                clustering_distance_columns = clustering_distance,
                                heatmap_legend_param = list(color_bar = "continuous",
                                                            legend_direction = "horizontal",
                                                            legend_width = unit(5, "cm"),
                                                            title_position = "lefttop"),
                                name = legend,
                                row_names_gp = gpar(fontsize = row_font_size, col = row_name_color),
                                column_names_gp = gpar(fontsize = col_font_size),
                                top_annotation = ha1)
  if(plot) {
    # Plot
    draw(ht1, heatmap_legend_side = "top")
  } else {
    # Return data.frame
    colnames(df) <- gsub(" ", "_", colnames(df))
    df <- df[, unlist(column_order(ht1))]
    if(kmeans) {
      df <- cbind(df, k = df_kmeans$cluster)
    }
    return <- df[unlist(row_order(ht1)),]
    data.frame(protein = row.names(return), return) %>%
      mutate(order = row_number())
  }
}


coef_variation<-function(x){
  coef=sd(x)/mean(x)
}

table_select <- function(frame,select){
  temp <- frame
  selected_frame <- list()
  for(i in 1:length(select)){
    selected_frame[[i]] <- temp[,grep(select[i],colnames(temp))]
  }
  temp2 <- do.call(cbind,selected_frame)
  temp1 <- temp[,c(1:2)]
  if(!any(colnames(temp1)=="Majority.protein.IDs")){
    temp1$Majority.protein.IDs <- temp$Majority.protein.IDs
  }
  if(!any(colnames(temp1)=="Proteins")){
    temp1$Proteins <- temp$Proteins
  }
  temp <- cbind(temp1,temp2)
  sig_col <- grep("significant",colnames(temp))
  
  allsig_col <- vector()
  for(i in 1:nrow(temp)){
    allsig_col[i] <- any(as.data.frame(temp[,sig_col])[i,]=="TRUE")
  }
  temp <- temp[allsig_col,]
  rownames(temp) <- 1:nrow(temp)
  
  
  return(temp)
}


########################
## annoation function ##
########################
rm_digit_end <- function(x){
  gsub("\\.\\d*$", "", x)
}

## the1stname ##--------------
the1stname <- function(gene.name){
  names <- strsplit(as.character(gene.name),";")[[1]][1]
  return(names)
}

## to entrezid function ##-------------
my_to_entrezid <- function(orgDB = org.Hs.eg.db, gene) {
  ids1 <- try(mapIds(x = orgDB, keys = gene, keytype = "SYMBOL", column = "ENTREZID"), silent = TRUE)
  ids2 <- try(mapIds(x = orgDB, keys = gene, keytype = "ENSEMBL", column = "ENTREZID"), silent = TRUE)
  ids3 <- try(mapIds(x = orgDB, keys = gene, keytype = "UNIPROT", column = "ENTREZID"), silent = TRUE)
  ids4 <- try(mapIds(x = orgDB, keys = gene, keytype = "ALIAS", column = "ENTREZID"), silent = TRUE)
  
  ids_lis <- list(ids1 = ids1, ids2 = ids2, ids3 = ids3, ids4 = ids4)
  
  ids_ind <- c(class(ids1) != "try-error", class(ids2) != "try-error", class(ids3) != "try-error", class(ids4) != "try-error")
  
  my_ids <- ids_lis[ids_ind]
  my_ids1 <- as.data.frame(do.call(cbind, my_ids), stringsAsFactors = F)
  
  my_ids1$id <- apply(my_ids1, 1, function(i){
    if(all(is.na(i))) {id = NA
    } else {
      id =  i[which(!is.na(i))[1]]
    }
    return(id)
  })
  return(my_ids1)
}

## Gene annotate function (return bgo)##-----------
Geneannotate <- function(geneid,keytype="ENTREZID",genedb){
  GOdata <- GO.db
  ##ago <- NA
  #ago<- AnnotationDbi::select(genedb, keys=geneid, columns = c("GO","ENTREZID","ENSEMBL", "SYMBOL", "PFAM", "UNIPROT", "GENENAME", "PATH"), keytype = keytype)
  ago<- AnnotationDbi::select(genedb, keys=geneid, columns = c("GO","ENTREZID", "SYMBOL", "PFAM", "GENENAME", "PATH"), keytype = keytype)
  ago <- ago %>% mutate(numb=c(1:nrow(ago)))
  # try(temgo<-AnnotationDbi::select(GOdata,key=ago$GO,columns = c("DEFINITION","TERM"),keytype = "GOID"))
  try(temgo<-AnnotationDbi::select(GOdata,key=ago$GO,columns = c("TERM"),keytype = "GOID"))
  temgo <- temgo %>% mutate(numb=c(1:nrow(ago)))
  bgo <- merge(ago,temgo,by="numb")
  
  kegg <- ko2name(ko = paste("ko","",bgo$PATH, sep = ""))
  names(kegg)[2] = "kegg"
  names(kegg)[1] = "PATH"
  kegg[ , 1] = gsub("ko", "", kegg[ , 1])
  kegg = kegg[!duplicated(kegg$kegg), ]
  
  bgo <- dplyr::left_join(x = bgo, y = kegg, by = "PATH")
  
  react = AnnotationDbi::select(reactome.db, keys = geneid, columns = c("REACTOMEID", "PATHNAME"), keytypes="ENTREZID") 
  names(react)[3] = "react"
  bgo <- dplyr::left_join(x = bgo, y = react, by = "ENTREZID")
  #################### test
  pfamAC = as.data.frame(PFAMID)
  names(pfamAC)[1:2] = c("PFAM", "PFAM_ID")
  
  bgo = dplyr::left_join(x = bgo, y = pfamAC, by = "PFAM")
  
  IDS <<- unique(bgo$ENTREZID)
  return(bgo)
}

## mergego function (return ) ## -----------
mergego <- function(x,bgo){
  library(dplyr)
  subgo <- bgo %>% dplyr::filter(ENTREZID==x) 
  subgo <- subgo[!duplicated(subgo$GO),] 
  GOCC <- subgo%>%dplyr::filter(ONTOLOGY=="CC") %>%dplyr::select("TERM")%>%unlist()%>% paste0(collapse = " ; ")
  GOBP <- subgo%>%dplyr::filter(ONTOLOGY=="BP") %>%dplyr::select("TERM")%>%unlist()%>% paste0(collapse = " ;  ")
  GOMF <- subgo%>%dplyr::filter(ONTOLOGY=="MF") %>%dplyr::select("TERM")%>%unlist()%>% paste0(collapse = " ;  ")
  
  #DECC <- subgo%>%dplyr::filter(ONTOLOGY=="CC") %>%dplyr::select("DEFINITION")%>%unlist()%>% paste0(collapse = " ;  ")
  #DEBP <- subgo%>%dplyr::filter(ONTOLOGY=="BP") %>%dplyr::select("DEFINITION")%>%unlist()%>% paste0(collapse = " ;  ")
  #DEMF <- subgo%>%dplyr::filter(ONTOLOGY=="MF") %>%dplyr::select("DEFINITION")%>%unlist()%>% paste0(collapse = " ;  ")
  #GOs <- c(GOCC,DECC,GOBP,DEBP,GOMF,DEMF)
  GOs <- c(GOCC,GOBP,GOMF)
  
  subpfam <- bgo %>% dplyr::filter(ENTREZID==x) 
  subpfam <- subpfam[!duplicated(subpfam$PFAM_ID),]
  pfam <- subpfam %>% dplyr::select("PFAM_ID") %>% unlist() %>% paste0(collapse = " ; ")
  
  
  sub_gene.descri <- bgo %>% dplyr::filter(ENTREZID==x)
  sub_gene.descri <- sub_gene.descri[!duplicated(sub_gene.descri$GENENAME),]
  gene.descri <- sub_gene.descri %>% dplyr::select("GENENAME") %>% unlist() %>% paste0(collapse = " ; ")
  
  subkegg <- bgo %>% dplyr::filter(ENTREZID==x) 
  subkegg <- subkegg[!duplicated(subkegg$kegg),]
  kegg <- subkegg %>% dplyr::select("kegg") %>% unlist() %>% paste0(collapse = " ; ")
  
  subreact <- bgo %>% dplyr::filter(ENTREZID==x) 
  subreact <- subreact[!duplicated(subreact$react),]
  react <- subreact %>% dplyr::select("react") %>% unlist() %>% paste0(collapse = " ; ")
  
  res = c(GOs, pfam, gene.descri, kegg, react)
  return(res)
}

## for gene list tool ##------------

get_msqrob_genelist_Env <- function () {
  if (!exists(".msqrob_genelist_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".msqrob_genelist_Env", new.env(), envir=envir)
  }
  get(".msqrob_genelist_Env", envir = .GlobalEnv)
}

get_imported_genelist_Env <- function () {
  if (!exists(".genelist_imported_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".genelist_imported_Env", new.env(), envir=envir)
  }
  get(".genelist_imported_Env", envir = .GlobalEnv)
}

get_all_genelist_Env <- function () {
  if (!exists(".genelist_all_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".genelist_all_Env", new.env(), envir=envir)
  }
  get(".genelist_all_Env", envir = .GlobalEnv)
}

# assign DEP-LFQ and DEG-RNAseq all contrasts with log2 fc to this envir in order to import to GSEA options
get_gsea_genelist_Env <- function () {
  if (!exists(".genelist_gsea_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".genelist_gsea_Env", new.env(), envir=envir)
  }
  get(".genelist_gsea_Env", envir = .GlobalEnv)
} 

get_msqrob_result_Env <- function () {
  if (!exists(".msqrob_result_Env ", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".msqrob_result_Env ", new.env(), envir=envir)
  }
  get(".msqrob_result_Env ", envir = .GlobalEnv)
}

## get all lists##-------
get_all_lists <- function(){
  ls(envir = get_all_genelist_Env())
}

## all msqrob siglist##-------
all_msqrob_siglist <- function(frame,alpha,lfc)
{
  if (is.integer(alpha))
    alpha <- as.numeric(alpha)
  if (is.integer(lfc))
    lfc <- as.numeric(lfc)
  assertthat::assert_that(inherits(frame, "data.frame"),
                          is.numeric(alpha), length(alpha) == 1, is.numeric(lfc),
                          length(lfc) == 1)
  # row_data <- rowData(dep, use.names = FALSE) %>% as.data.frame()
  # res <- as.data.frame(rowData(object[["FilteredResult"]]))
  res <- frame
  if(!length(grep("_",res$Gene.names))==0){
    message("Your input file version does not match the requirements, please rename the file as required, 
            otherwise the analysis will most likely fail.")
  }else{
    res <- res[!res$Majority.protein.IDs==res$Gene.names,]
  }
  
  
  
  if (any(!c("Gene.names", "Majority.protein.IDs") %in% colnames(res))) {
    stop("'Gene.names' and/or 'ID' columns are not present in '",
         deparse(substitute(frame)), "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (length(grep("adjPval|logFC", colnames(res))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(frame)), "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }
  cols_p <- grep("adjPval", colnames(res))
  cols_diff <- grep("logFC", colnames(res))
  cols_sig <- grep("significant", colnames(res))
  
  genelist_env = get_msqrob_genelist_Env()
  genelist_all_env = get_all_genelist_Env()
  rm(list = ls(envir = genelist_env),envir = genelist_env)
  
  
  for(i in 1:length(cols_p)){
    tem_all = res[,cols_sig[i]]=="TRUE"
    temp_all = res[which(tem_all),c(which(colnames(res)=="Gene.names"),which(colnames(res)=="Majority.protein.IDs")
                            ,cols_diff[i],cols_p[i])]
    if(nrow(temp_all)>0){
      colnames(temp_all) = c("symbol","ID","L2FC","p.adj")
      temp_all$significant = TRUE
      temp_all = temp_all[ , c(1, 3, 2, 4, 5)]
      assign(value = temp_all, x = paste(colnames(res)[cols_p[i]] %>% gsub("adjPval",paste(alpha,lfc,"Sig",sep = "_"),.), sep = "") ,
             envir = genelist_env)
      assign(value = temp_all, x = paste(colnames(res)[cols_p[i]] %>% gsub("adjPval",paste(alpha,lfc,"Sig",sep = "_"),.), sep = "") ,
             envir = genelist_all_env)
    }
    
    tem_up = (res[,cols_sig[i]]=="TRUE") & (res[,cols_diff[i]] > lfc)
    temp_up = res[which(tem_up),c(which(colnames(res)=="Gene.names"),which(colnames(res)=="Majority.protein.IDs")
                                  ,cols_diff[i],cols_p[i])]
    if(nrow(temp_up)>0){
      colnames(temp_up) = c("symbol","ID","L2FC","p.adj")
      temp_up$significant = TRUE
      temp_up = temp_up[ , c(1, 3, 2, 4, 5)]
      assign(value = temp_up, x = paste(colnames(res)[cols_p[i]] %>% gsub("adjPval",paste(alpha,lfc,"Up",sep = "_"),.), sep = "") ,
             envir = genelist_env)
      assign(value = temp_up, x = paste(colnames(res)[cols_p[i]] %>% gsub("adjPval",paste(alpha,lfc,"Up",sep = "_"),.), sep = "") ,
             envir = genelist_all_env)
    }
    
    tem_dn = (res[,cols_sig[i]]=="TRUE") & (res[,cols_diff[i]] < -lfc)
    temp_dn = res[which(tem_dn),c(which(colnames(res)=="Gene.names"),which(colnames(res)=="Majority.protein.IDs")
                                  ,cols_diff[i],cols_p[i])]
    if(nrow(temp_dn)>0){
      colnames(temp_dn) = c("symbol","ID","L2FC","p.adj")
      temp_dn$significant = TRUE
      temp_dn = temp_dn[ , c(1, 3, 2, 4, 5)]
      assign(value = temp_dn, x = paste(colnames(res)[cols_p[i]] %>% gsub("adjPval",paste(alpha,lfc,"Down",sep = "_"),.), sep = "") ,
             envir = genelist_env)
      assign(value = temp_dn, x = paste(colnames(res)[cols_p[i]] %>% gsub("adjPval",paste(alpha,lfc,"Down",sep = "_"),.), sep = "") ,
             envir = genelist_all_env)
    }
  }
  return(ls(envir = get_msqrob_genelist_Env()))
}
# 使 用 例 all_msqrob_siglist(msqrob_result,lfc = 1,alpha = 0.05)


## module function ##-------------
Genelist_UImodule <- function(id, label = "Counter"){
  
  ns <- NS(id)
  tagList(
    # navbarPage("Genelist tool",
    fluidRow(
      column(width = 12,
             tabBox(title = "Genelist tool", width = 12,
                    tabPanel("Generate gene list by cutoff",
                             column(width = 6,
                                    h3("Significant gene lists by threshold"),
                                    # fluidRow(
                                    #   conditionalPanel(condition = paste("input['", id, "-genelistFrom']", " ==  'DEG' || input.threshold_method == 'intersect'" , sep = ""), #ns = ns, #input.threshold_method == 'intersect' ||  paste("input['", id, "-genelistFrom']", " ==  'DEG'" , sep = "")
                                    #                    box(numericInput(ns("p_genelist"),
                                    #                                     "adjusted P cutoff",
                                    #                                     min = 0.0001, max = 0.1, value = 0.05),
                                    #                        width = 6),
                                    #                    box(numericInput(ns("lfc_genelist"),
                                    #                                     "log2 fold change cutoff",
                                    #                                     min = 0, max = 10, value = 1),
                                    #                        width = 6)  
                                    #   )
                                    #   
                                    #   # box(numericInput(ns("p_genelist"),
                                    #   #                       "adjusted P cutoff",
                                    #   #                       min = 0.0001, max = 0.1, value = 0.05),
                                    #   #          width = 6),
                                    #   # box(numericInput(ns("lfc_genelist"),
                                    #   #                       "log2 fold change cutoff",
                                    #   #                       min = 0, max = 10, value = 1),
                                    #   #          width = 6)                            
                                    # ),
                                    fluidRow(
                                      box(selectizeInput(ns("genelistFrom"),
                                                         "Generate gene list from",
                                                         choices=c("MSqRob"," "), selected = " ",
                                                         multiple = F), width = 6),
                                      # box(selectizeInput(ns("selected_compare"),
                                      #                    "Select groups",
                                      #                    choices=NULL,
                                      #                    multiple = F),width = 4),
                                      box(uiOutput(ns("selected_compare_ui")),width = 6)
                                    ),
                                    fluidRow(
                                      
                                      box(selectizeInput(ns("type_genelist"),
                                                         "Significant type",
                                                         choices=c("Significant(Up&Down)","Upregulated(Up)","Downregulated(Down)"),
                                                         multiple = F), width = 6)                            
                                    ),
                                    fluidRow(
                                      column(width = 12,
                                             actionBttn(inputId = ns("Generate_genelist"),
                                                        label = "Generate",
                                                        style = "jelly",
                                                        color = "primary")                              
                                      )                             
                                    )
                             ),
                             column(width = 5,
                                    br(),
                                    br(),
                                    br(),
                                    h4("existing gene lists"),
                                    htmlOutput(ns("existed_lists"))
                             )
                    ),
                    tabPanel("Import gene list",
                             h3("Import gene list"),
                             fluidRow(
                               column(h4("Import by pasting gene list"),
                                      textAreaInput(inputId = ns("text_input_genelist"),
                                                    label = "Please paste your gene list",
                                                    placeholder = "TP53\nMAGED2\n",
                                                    rows = 8, width = "350px"),
                                      textInput(ns("name_genelist"),"the name of imported gene list","Imported list 1", width = "350px"),
                                      actionBttn(
                                        inputId = ns("Import_genelist"),
                                        label = "Import",
                                        style = "jelly",
                                        color = "primary"
                                      ), 
                                      br(),
                                      br(),
                                      textOutput(ns("input_genelist")),
                                      tableOutput(ns("dataframe"))                                                               
                                      ,width = 5),
                               column(h4("existing imported gene lists"),
                                      htmlOutput(ns("existed_imported_lists")),
                                      width = 7)
                             ),
                             # tags$br(),
                             # textInput(ns("name_genelist"),"the name of import Genelist","Imported list 1", width = "350px"),
                             
                             # br(),
                             # br(),
                             # textOutput(ns("input_genelist")),
                             # tableOutput(ns("dataframe"))
                    )         
                    
             )
      )
    )
    
    # )
    
  )
  
}


Genelist_server_module <- function(id,alpha,object,lfc,condition_matched_frame) {
  moduleServer(
    id,
    function(input, output, session,id=id) {
      
      # rm all objects of the envir that we build, in order to make each run app rm objects last run
      rm(list = ls(envir = get_msqrob_genelist_Env()),envir = get_msqrob_genelist_Env())
      rm(list = ls(envir = get_imported_genelist_Env()),envir = get_imported_genelist_Env())
      rm(list = ls(envir = get_all_genelist_Env()),envir = get_all_genelist_Env())
      rm(list = ls(envir = get_gsea_genelist_Env()),envir = get_gsea_genelist_Env())
      
      msqrob_frame <- reactive({
        msqrob_frame <- result_from_object(object=object(),condition_matched_frame())
        # test_frame <<- msqrob_frame
      })
      
      
      lists_msqrob <- reactive({
        lists_msqrob <- all_msqrob_siglist(frame = msqrob_frame(),alpha = alpha(),lfc = lfc())
        # testlist <<- lists_msqrob
        })
      
      
      output$existed_lists = renderUI({
        validate(need(!is.null(msqrob_frame()), "Please go to the MSqRob options, and do the differential analysis first"))
        if(!is.null(msqrob_frame())) {
          lis <- unique(c(lists_msqrob()))
          # lis <- unique(ls(envir = get_msqrob_genelist_Env()))
        } 
        HTML(paste(lis, collapse = "</br>"))
      })
      
      # 
      # output$existed_lists <- renderUI({
      #   HTML(paste(ls(envir = get_msqrob_genelist_Env()), collapse = "</br>"))
      # })
      
      
      output$input_genelist = renderText({
        thename <- input$name_genelist
      })
      
      #
      #       output$Vennplot = renderPlot({NULL})
      
      
      cat("bb")
      
      output$selected_compare_ui <- renderUI({
        froms <- input$genelistFrom
        if(froms == "MSqRob"){
          validate(need(!is.null(msqrob_frame()), "Please go to the MSqRob options, and do the differential analysis first"))
          frame <- msqrob_frame()
          cols <- grep(".significant", colnames(frame))
          names <- colnames(frame)[cols]
          names <- gsub(".significant", "", names)
          selectizeInput(session$ns("selected_compare"),"Select groups",choices=names, multiple = F)#getDefaultReactiveDomain()
        }
      })

      
      ## Event generate  genelist
      observeEvent(input$Generate_genelist,{
        # if(length(input$selected_compare)==0){
        #   sendSweetAlert(
        #     session = session,
        #     title = "Warning !",
        #     text = "please choose compare",
        #     type = "warning"
        #   )
        # }else{
          # aaa <- all_msqrob_siglist(object = msqrob_object,alpha = alpha,lfc = lfc)
          froms = input$genelistFrom
          list_type = input$type_genelist
          # test_type <<- list_type
          type <- switch(list_type,
                              `Significant(Up&Down)` = "Sig",
                              `Upregulated(Up)`= "Up",
                              `Downregulated(Down)` ="Down")
          
          
          # if(froms == "Msqrob"){
          thelist <- get(paste0(input$selected_compare,".",alpha(),"_",lfc(),
                                "_",type),envir = get_msqrob_genelist_Env())
          # testthelist <<- paste0(input$selected_compare,".",alpha,"_",lfc,"_",type) 
          the_name <- paste0(input$selected_compare,".",alpha(),"_",lfc(),
                             "_",type)
          
          if(nrow(thelist)==0){
            sendSweetAlert(
              session = session,
              title = "Warning !",
              text = "no expected protein under this threshold\nPlease back to Msqrob option to do hypothesis test",
              type = "warning"
            )
          }
          # thelist <- get_msqrob_siglist(msqrob_object = msqrob(),
          #                                 alpha = input$p_genelist , lfc = input$lfc_genelist,
          #                                 compare = input$selected_compare,type = list_type) 
          else{
            sendSweetAlert(
              session = session,
              title = "Success",
              text = paste("the list saved as: ", the_name,
                           "  (","total ",nrow(thelist)," proteins",")",
                           sep = ""),
              type = "success"
            )
          }
          
      })
      
      ## Event import  genelist
      observeEvent(input$Import_genelist, {
        texts <- input$text_input_genelist
        ## warning if
        if(nchar(texts) == 0){
          sendSweetAlert(
            session = session,
            title = "Warning !",
            text = "your input text is empty",
            type = "warning"
          )
          import_table = NULL
        }else{
          import_table = try(read.table(text = texts,header = F,sep = "\t",stringsAsFactors = F),silent = T)
          
          ## read in talbe
          if(nrow(import_table) ==1 ){
            ## warning if only have one gene
            sendSweetAlert(
              session = session,
              title = "Warning !",
              text = "please input list with multiple genes",
              type = "warning"
            )
            import_table = NULL
          }else{
            if(ncol(import_table)==1){
              colnames(import_table) = "symbol"
            }else if(ncol(import_table)==2){
              colnames(import_table) = c("symbol","FC")
            }else if(ncol(import_table)>2){
              import_table = import_table[,1:2]
              colnames(import_table) = c("symbol","FC")
            }
            ## assign table
            assign(x = input$name_genelist,value = import_table,envir = get_imported_genelist_Env())
            assign(x = input$name_genelist,value = import_table,envir = get_all_genelist_Env())
            
            sendSweetAlert(
              session = session,
              title = "Success",
              text = paste("list import success: (total ",nrow(import_table),")",sep = ""),
              type = "success"
            )
          }
        }
        output$dataframe <- renderTable({
          import_table
        })
        
        output$existed_imported_lists <- renderUI({
          HTML(paste(ls(envir = get_imported_genelist_Env()), collapse = "</br>"))
        })
        
      })
    }
  )
}


## ORA function ##---------
the1stname <- function(gene.name){
  names <- strsplit(as.character(gene.name),";")[[1]][1]
  return(names)
}

## to entrezid function ##-------------
my_to_entrezid <- function(orgDB = org.Hs.eg.db, gene) {
  ids1 <- try(mapIds(x = orgDB, keys = gene, keytype = "SYMBOL", column = "ENTREZID"), silent = TRUE)
  ids2 <- try(mapIds(x = orgDB, keys = gene, keytype = "ENSEMBL", column = "ENTREZID"), silent = TRUE)
  ids3 <- try(mapIds(x = orgDB, keys = gene, keytype = "UNIPROT", column = "ENTREZID"), silent = TRUE)
  ids4 <- try(mapIds(x = orgDB, keys = gene, keytype = "ALIAS", column = "ENTREZID"), silent = TRUE)
  
  ids_lis <- list(ids1 = ids1, ids2 = ids2, ids3 = ids3, ids4 = ids4)
  
  ids_ind <- c(class(ids1) != "try-error", class(ids2) != "try-error", class(ids3) != "try-error", class(ids4) != "try-error")
  
  my_ids <- ids_lis[ids_ind]
  my_ids1 <- as.data.frame(do.call(cbind, my_ids), stringsAsFactors = F)
  
  my_ids1$id <- apply(my_ids1, 1, function(i){
    if(all(is.na(i))) {id = NA
    } else {
      id =  i[which(!is.na(i))[1]]
    }
    return(id)
  })
  return(my_ids1)
}

## Gene annotate function (return bgo)##-----------
Geneannotate <- function(geneid,keytype="ENTREZID",genedb){
  GOdata <- GO.db
  ##ago <- NA
  #ago<- AnnotationDbi::select(genedb, keys=geneid, columns = c("GO","ENTREZID","ENSEMBL", "SYMBOL", "PFAM", "UNIPROT", "GENENAME", "PATH"), keytype = keytype)
  ago<- AnnotationDbi::select(genedb, keys=geneid, columns = c("GO","ENTREZID", "SYMBOL", "PFAM", "GENENAME", "PATH"), keytype = keytype)
  ago <- ago %>% mutate(numb=c(1:nrow(ago)))
  # try(temgo<-AnnotationDbi::select(GOdata,key=ago$GO,columns = c("DEFINITION","TERM"),keytype = "GOID"))
  try(temgo<-AnnotationDbi::select(GOdata,key=ago$GO,columns = c("TERM"),keytype = "GOID"))
  temgo <- temgo %>% mutate(numb=c(1:nrow(ago)))
  bgo <- merge(ago,temgo,by="numb")
  
  kegg <- ko2name(ko = paste("ko","",bgo$PATH, sep = ""))
  names(kegg)[2] = "kegg"
  names(kegg)[1] = "PATH"
  kegg[ , 1] = gsub("ko", "", kegg[ , 1])
  kegg = kegg[!duplicated(kegg$kegg), ]
  
  bgo <- dplyr::left_join(x = bgo, y = kegg, by = "PATH")
  
  react = AnnotationDbi::select(reactome.db, keys = geneid, columns = c("REACTOMEID", "PATHNAME"), keytypes="ENTREZID") 
  names(react)[3] = "react"
  bgo <- dplyr::left_join(x = bgo, y = react, by = "ENTREZID")
  #################### test
  pfamAC = as.data.frame(PFAMID)
  names(pfamAC)[1:2] = c("PFAM", "PFAM_ID")
  
  bgo = dplyr::left_join(x = bgo, y = pfamAC, by = "PFAM")
  
  IDS <<- unique(bgo$ENTREZID)
  return(bgo)
}

## mergego function (return ) ## -----------
mergego <- function(x,bgo){
  library(dplyr)
  subgo <- bgo %>% dplyr::filter(ENTREZID==x) 
  subgo <- subgo[!duplicated(subgo$GO),] 
  GOCC <- subgo%>%dplyr::filter(ONTOLOGY=="CC") %>%dplyr::select("TERM")%>%unlist()%>% paste0(collapse = " ; ")
  GOBP <- subgo%>%dplyr::filter(ONTOLOGY=="BP") %>%dplyr::select("TERM")%>%unlist()%>% paste0(collapse = " ;  ")
  GOMF <- subgo%>%dplyr::filter(ONTOLOGY=="MF") %>%dplyr::select("TERM")%>%unlist()%>% paste0(collapse = " ;  ")
  
  #DECC <- subgo%>%dplyr::filter(ONTOLOGY=="CC") %>%dplyr::select("DEFINITION")%>%unlist()%>% paste0(collapse = " ;  ")
  #DEBP <- subgo%>%dplyr::filter(ONTOLOGY=="BP") %>%dplyr::select("DEFINITION")%>%unlist()%>% paste0(collapse = " ;  ")
  #DEMF <- subgo%>%dplyr::filter(ONTOLOGY=="MF") %>%dplyr::select("DEFINITION")%>%unlist()%>% paste0(collapse = " ;  ")
  #GOs <- c(GOCC,DECC,GOBP,DEBP,GOMF,DEMF)
  GOs <- c(GOCC,GOBP,GOMF)
  
  subpfam <- bgo %>% dplyr::filter(ENTREZID==x) 
  subpfam <- subpfam[!duplicated(subpfam$PFAM_ID),]
  pfam <- subpfam %>% dplyr::select("PFAM_ID") %>% unlist() %>% paste0(collapse = " ; ")
  
  
  sub_gene.descri <- bgo %>% dplyr::filter(ENTREZID==x)
  sub_gene.descri <- sub_gene.descri[!duplicated(sub_gene.descri$GENENAME),]
  gene.descri <- sub_gene.descri %>% dplyr::select("GENENAME") %>% unlist() %>% paste0(collapse = " ; ")
  
  subkegg <- bgo %>% dplyr::filter(ENTREZID==x) 
  subkegg <- subkegg[!duplicated(subkegg$kegg),]
  kegg <- subkegg %>% dplyr::select("kegg") %>% unlist() %>% paste0(collapse = " ; ")
  
  subreact <- bgo %>% dplyr::filter(ENTREZID==x) 
  subreact <- subreact[!duplicated(subreact$react),]
  react <- subreact %>% dplyr::select("react") %>% unlist() %>% paste0(collapse = " ; ")
  
  res = c(GOs, pfam, gene.descri, kegg, react)
  return(res)
}


## goAnalysis ##-------------
goAnalysis <- function(df, df_with_lg2fc = FALSE, organism="Human", species_df){
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)
  
  # if(organism == "human"){
  #   library(org.Hs.eg.db)
  #   orgDB <<- org.Hs.eg.db
  #   #kegg_organism <- "hsa"
  # }
  # if(organism == "mouse"){
  #   library(org.Mm.eg.db)
  #   orgDB <<- org.Mm.eg.db
  #   #kegg_organism <- "mmu"
  # }
  # if(organism == "rat"){
  #   library(org.Rn.eg.db)
  #   orgDB <<- org.Rn.eg.db
  #   #kegg_organism <- "rno"
  # }
  
  
  # ids <- bitr(df$name, fromType = "SYMBOL",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB)
  # ids1 <- try(bitr(setdiff(df$name, ids$SYMBOL), fromType = "UNIPROT",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB))
  # if(class(ids1) != "try-error") {
  #   names(ids1)[1] = "SYMBOL"
  #   ids = rbind(ids, ids1)
  # }
  
  # names(ids)[1] = "name"
  #  ids <- inner_join(ids, df, by = "name")
  
  ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)
  
  ids <- inner_join(ids2, df, by = "name")
  
  if(df_with_lg2fc){
    ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
    de = ids$fc
    names(de) = unlist(ids$ENTREZID)
    de = sort(de, decreasing = T)
  } else {
    ids <- ids[!is.na(ids$ENTREZID), ]
    de = unlist(ids$ENTREZID)
    names(de) = de
  }
  
  reat_ALL <- clusterProfiler::enrichGO(gene = names(de), OrgDb = orgDB, ont = "ALL", 
                                        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  
  reat_BP <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "BP")
  reat_BP@ontology = "BP"
  
  reat_CC <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "CC")
  reat_CC@ontology = "CC"
  
  reat_MF <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "MF")
  reat_MF@ontology = "MF"
  
  reat <- list(ALL = reat_ALL, BP = reat_BP, CC = reat_CC, MF = reat_MF, de = de)
  return(reat)
}

## giveGO_res_and_table ##--------
giveGO_res_and_table <- function(reat, ont = "BP", pCutoff = 0.05, p.adj.cutoff = 0.05, q.cutoff = 0.2, simplify = FALSE){
  res <- reat[[ont]]
  sig_res <- clusterProfiler.dplyr::filter(res, pvalue < pCutoff, p.adjust < p.adj.cutoff, qvalue < q.cutoff)
  
  if(simplify) {
    if(ont == "ALL") {
      sig_res@ontology = "BP"
      sig_res <- clusterProfiler::simplify(sig_res, cutoff=0.7, by="p.adjust", select_fun=min)
    } else {
      sig_res <- clusterProfiler::simplify(sig_res, cutoff=0.7, by="p.adjust", select_fun=min)
    }
  }
  
  res_table <- as.data.frame(res)
  sig_res_table <- as.data.frame(sig_res)
  
  sig_res@result = sig_res@result[ , c(2:ncol(sig_res@result), 1)]
  
  table <- list(all_table = res_table, sig_table = sig_res_table, sig_res = sig_res, de = reat[["de"]])
  
  return(table)
}

## Heatplot ##------------------
Heatplot <-  function (x, showCategory = 30, foldChange = NULL) {
  n <- update_n(x, showCategory)
  geneSets <- extract_geneSets(x, n)
  foldChange <- fc_readable(x, foldChange)
  d <- list2df(geneSets)
  if (!is.null(foldChange)) {
    d$foldChange <- foldChange[as.character(d[, 2])]
    p <- ggplot(d, aes_(~Gene, ~categoryID)) + geom_tile(aes_(fill = ~foldChange), 
                                                         color = "white") + scale_fill_continuous(low = "blue", 
                                                                                                  high = "red", name = "log2 fold change")
  }
  else {
    p <- ggplot(d, aes_(~Gene, ~categoryID)) + geom_tile(color = "white")
  }
  p + xlab(NULL) + ylab(NULL) + theme_minimal() + theme(panel.grid.major = element_blank(), 
                                                        axis.text.x = element_text(angle = 60, hjust = 1))
}     
environment(Heatplot) = asNamespace("enrichplot")

## Cnetplot ##----------------
Cnetplot <- function (x, showCategory = 5, foldChange = NULL, layout = "kk", 
                      colorEdge = FALSE, circular = FALSE, node_label = "all", 
                      ...) {
  node_label <- match.arg(node_label, c("category", "gene", 
                                        "all", "none"))
  if (circular) {
    layout <- "linear"
    geom_edge <- geom_edge_arc
  }
  else {
    geom_edge <- geom_edge_link
  }
  geneSets <- extract_geneSets(x, showCategory)
  g <- list2graph(geneSets)
  foldChange <- fc_readable(x, foldChange)
  size <- sapply(geneSets, length)
  V(g)$size <- min(size)/2
  n <- length(geneSets)
  V(g)$size[1:n] <- size
  if (colorEdge) {
    E(g)$category <- rep(names(geneSets), sapply(geneSets, 
                                                 length))
    edge_layer <- geom_edge(aes_(color = ~category), alpha = 0.8)
  }
  else {
    edge_layer <- geom_edge(alpha = 0.8, colour = "darkgrey")
  }
  if (!is.null(foldChange)) {
    fc <- foldChange[V(g)$name[(n + 1):length(V(g))]]
    V(g)$color <- NA
    V(g)$color[(n + 1):length(V(g))] <- fc
    palette <- fc_palette(fc)
    p <- ggraph(g, layout = layout, circular = circular) + 
      edge_layer + geom_node_point(aes_(color = ~as.numeric(as.character(color)), 
                                        size = ~size)) + scale_color_gradientn(name = "log2 fold change", 
                                                                               colors = palette, na.value = "#E5C494")
  }
  else {
    V(g)$color <- "#B3B3B3"
    V(g)$color[1:n] <- "#E5C494"
    p <- ggraph(g, layout = layout, circular = circular) + 
      edge_layer + geom_node_point(aes_(color = ~I(color), 
                                        size = ~size))
  }
  p <- p + scale_size(range = c(3, 10), breaks = unique(round(seq(min(size), 
                                                                  max(size), length.out = 4)))) + theme_void()
  if (node_label == "category") {
    p <- p + geom_node_text(aes_(label = ~name), data = p$data[1:n, 
    ], repel = TRUE)
  }
  else if (node_label == "gene") {
    p <- p + geom_node_text(aes_(label = ~name), data = p$data[-c(1:n), 
    ], repel = TRUE)
  }
  else if (node_label == "all") {
    p <- p + 
      geom_node_text(aes_(label = ~name), data = p$data, repel = TRUE, max.overlaps = Inf)
    # geom_node_text(aes_(label = ~name), data = p$data[-c(1:n), 
    # ], repel = TRUE) +
    # geom_node_text(aes_(label = ~name), data = p$data[1:n, 
    # ], repel = TRUE)
  }
  return(p)
}

fc_palette <- function (fc) 
{
  if (all(fc > 0, na.rm = TRUE)) {
    palette <- color_palette(c("blue", "red"))
  }
  else if (all(fc < 0, na.rm = TRUE)) {
    palette <- color_palette(c("green", "blue"))
  }
  else {
    palette <- color_palette(c("darkgreen", "#0AFF34", 
                               "#B3B3B3", "#FF6347", "red"))
  }
  return(palette)
}

environment(Cnetplot) = asNamespace("enrichplot")
environment(fc_palette) = asNamespace("enrichplot")

## Emapplot ##---------------
Emapplot <- function (x, showCategory = 30, color = "p.adjust", layout = "kk", 
                      ...) 
{
  n <- update_n(x, showCategory)
  geneSets <- geneInCategory(x)
  y <- as.data.frame(x)
  if (is.numeric(n)) {
    y <- y[1:n, ]
  }
  else {
    y <- y[match(n, y$Description), ]
    n <- length(n)
  }
  if (n == 0) {
    stop("no enriched term found...")
  }
  else if (n == 1) {
    g <- graph.empty(0, directed = FALSE)
    g <- add_vertices(g, nv = 1)
    V(g)$name <- y$Description
    V(g)$color <- "red"
    return(ggraph(g) + geom_node_point(color = "red", 
                                       size = 5) + geom_node_text(aes_(label = ~name)))
  }
  else {
    id <- y[, 1]
    geneSets <- geneSets[id]
    n <- nrow(y)
    w <- matrix(NA, nrow = n, ncol = n)
    colnames(w) <- rownames(w) <- y$Description
    for (i in 1:n) {
      for (j in i:n) {
        w[i, j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
      }
    }
    wd <- melt(w)
    wd <- wd[wd[, 1] != wd[, 2], ]
    wd <- wd[!is.na(wd[, 3]), ]
    g <- graph.data.frame(wd[, -3], directed = FALSE)
    E(g)$width = sqrt(wd[, 3] * 5)
    g <- delete.edges(g, E(g)[wd[, 3] < 0.2])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == 
                                                        y$Description)))
    cnt <- sapply(geneSets[idx], length)
    V(g)$size <- cnt
    colVar <- y[idx, color]
    V(g)$color <- colVar
  }
  p <- ggraph(g, layout = layout)
  if (length(E(g)$width) > 0) {
    p <- p + geom_edge_link(alpha = 0.8, aes_(width = ~I(width)), 
                            colour = "darkgrey")
  }
  p + geom_node_point(aes_(color = ~color, size = ~size)) + 
    geom_node_text(aes_(label = ~name), repel = TRUE) + theme_void() + 
    scale_color_continuous(low = "red", high = "blue", 
                           name = color, guide = guide_colorbar(reverse = TRUE)) + 
    scale_size(range = c(3, 8))
}
environment(Emapplot) = asNamespace("enrichplot")

## my_barplot ##-------------
my_barplot <- function(res, ShowCategory = 20, color = "p.adjust", ont = "BP", Split = FALSE) {
  library(clusterProfiler.dplyr)
  library(ggplot2)
  if(ont == "ALL") {
    if(!Split){
      sig_res_for_ALL <- res$sig_res %>% clusterProfiler.dplyr::arrange(p.adjust)
      barplot(sig_res_for_ALL, showCategory = ShowCategory, color = color) + facet_grid(ONTOLOGY~., scales = "free", space = "free")
    } else {
      barplot(res$sig_res, showCategory = ShowCategory, color = color, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free", space = "free")
    }
  } else {
    barplot(res$sig_res, showCategory = ShowCategory, color = color)
  }
}

## my_dotplot ##----------
my_dotplot <- function(res, ShowCategory = 20, color = "p.adjust", ont = "BP", Split = FALSE) {
  library(clusterProfiler.dplyr)
  library(ggplot2)
  if(ont == "ALL") {
    if(!Split){
      sig_res_for_ALL <- res$sig_res %>% clusterProfiler.dplyr::arrange(p.adjust)
      dotplot(sig_res_for_ALL, showCategory = ShowCategory, color = color) + facet_grid(ONTOLOGY~., scales = "free", space = "free")
    } else {
      dotplot(res$sig_res, showCategory = ShowCategory, color = color, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free", space = "free")
    }
  } else {
    dotplot(res$sig_res, showCategory = ShowCategory, color = color)
  }
}

## my_dotplot_opt ##-------------
my_dotplot_opt <- function(res, color = "p.adjust", size = "Count", title = "", decreasing = TRUE, ShowCategory = 20){
  df = as.data.frame(res$sig_res)
  
  if(nrow(df) > ShowCategory) {
    df = df[1:ShowCategory, ]
  }
  df$x = ""
  orderBy = color
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels = unique(df$Description[idx]))
  
  p_up = as.numeric(format(max(df[ , color]), digits = 2))
  p_low = as.numeric(format(min(df[ , color]), digits = 2))
  
  if(color == "p.adjust"){
    col.name = "FDR"
  }
  if(color == "pvalue"){
    col.name = "Pvalue"
  }
  
  
  a =  ggplot(df, aes_string(x = "x", y = "Description", size = size, color = color)) + 
    geom_point() + 
    geom_point(shape = 21, color = "black", stroke = 0.6, position = "identity")+
    scale_colour_gradient2(low = "#df293f", mid = ifelse(!max(df[ , color]) == min(df[ , color]), "#f9f9f9", "#f19c9b"), high = "#959cc8", midpoint = min(df[ , color]) + (max(df[ , color])-min(df[ , color]))/2, name = col.name, breaks = unique(c(min(df[ , color]), max(df[ , color]))), limits = c(min(df[ , color]), max(df[ , color])), labels = unique(c(format(p_low, scientific = TRUE), format(p_up, scientific = TRUE)))) + 
    ylab(NULL) + 
    ggtitle(title) + 
    scale_size(name = "Number of proteins found", range = c(3, 8), breaks = c(min(df[ , "Count"]), max(df[ , "Count"])), limits = c(min(df[ , "Count"]), max(df[ , "Count"]))) +
    theme(panel.background=element_blank(),
          axis.ticks=element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key=element_blank(),
          legend.spacing = unit(3, "lines"),
          legend.justification="bottom", legend.position=c(0.15,0),#0.15
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          #plot.margin = margin(t = 0, r = 0, b = 0.1, l = 15, unit = "cm"),
          title = element_text(size = 10)#all title elements: plot, axes, legends
    )+
    labs(x = NULL) +
    guides(
      color = guide_colorbar(order = 1, frame.colour = "black", frame.linetype = 1, frame.linewidth = 1, draw.ulim = TRUE, draw.llim = TRUE, raster = TRUE, ticks = FALSE, label = TRUE, title.position = "left", title.theme = element_text(angle = 90,  hjust = 0.5, size = 10), label.theme = element_text(angle = 90,  hjust = 0.5, size = 10), reverse = TRUE, barwidth = unit(0.5, units = "cm"), barheight = unit(2.5, units = "cm")),
      size = guide_legend(order = 0, title.position = "right", title.theme = element_text(angle = 90,  hjust = 0.5, size = 10), label.theme = element_text(size = 10), override.aes = list(alpha = 1, bg = "grey", color = "black"), reverse = TRUE, keyheight = 4)
    ) +
    scale_x_discrete(expand = expand_scale(mult = 0, add = c(0.008, 0.25)))+#0.005, 0.25
    scale_y_discrete()
  return(a)
}

##  my_emaplot ##-------
my_emaplot <- function(res, ShowCategory = 30, color = "p.adjust", layout = "kk", ont = "BP") {
  if(ont == "ALL") {
    sig_res_for_ALL <- res$sig_res %>% clusterProfiler.dplyr::arrange(p.adjust)
    try(Emapplot(sig_res_for_ALL, showCategory = ShowCategory, color = color, layout = layout), silent = T)
  } else {
    try(Emapplot(res$sig_res, showCategory = ShowCategory, color = color, layout = layout), silent = T)
  }
}

## my_cnetplot##------------
my_cnetplot <- function(res, ShowCategory = 5, circular = TRUE, colorEdge = TRUE, df_with_lg2fc = FALSE, ont = "BP") {
  if(ont == "ALL") {
    sig_res_for_ALL <- res$sig_res %>% clusterProfiler.dplyr::arrange(p.adjust)
    if(df_with_lg2fc) {
      try(Cnetplot(x = sig_res_for_ALL, showCategory = ShowCategory, foldChange = res$de, circular = circular, colorEdge = colorEdge), silent = T)
    } else {
      try(Cnetplot(x = sig_res_for_ALL, showCategory = ShowCategory, foldChange=NULL, circular = circular, colorEdge = colorEdge), silent = T)
    }
    
  } else {
    if(df_with_lg2fc) {
      try(Cnetplot(x = res$sig_res, showCategory = ShowCategory, foldChange = res$de, circular = circular, colorEdge = colorEdge), silent = T)
    } else {
      try(Cnetplot(x = res$sig_res, showCategory = ShowCategory, foldChange=NULL, circular = circular, colorEdge = colorEdge), silent = T)
    }
  }
}

## my_goplot ##----------
my_goplot <- function(res, ShowCategory = 10, color = "p.adjust", ont = "BP", Layout = "kk",circular = TRUE) {
  library(ggplotify)
  if(! ont == "ALL"){
    a = try(goplot(res$sig_res, showCategory = ShowCategory,color = color, layout = Layout,geom = "text", circular = circular), silent = T)
    if(!class(a) == "try-error") {
      print(a)
    } else {
      print("can not plot")
    }
  } else {
    print("Ontology ALL: can not plot for goplot")
  }
}

## my_plotGOgraph ##---------
my_plotGOgraph <- function(res, firstSigNodes = 10, ont = "BP") {
  library(ggplotify)
  if(! ont == "ALL"){
    a = function() {try(plotGOgraph(res$sig_res, firstSigNodes = firstSigNodes), silent = T)}
    a = try(as.ggplot(a), silent = T)
    if(!class(a) == "try-error") {
      print(a)
    } else {
      print("can not plot")
    }
  } else {
    print("Ontology ALL: can not plot for plotGOgraph")
  }
}

## my_heatplot ##----------
my_heatplot <- function(res, ShowCategory = 30, df_with_lg2fc = FALSE, ont = "BP") {
  if(ont == "ALL") {
    sig_res_for_ALL <- res$sig_res %>% clusterProfiler.dplyr::arrange(p.adjust)
    if(df_with_lg2fc) {
      try(Heatplot(sig_res_for_ALL, showCategory = ShowCategory, foldChange = res$de), silent = T)
    } else {
      try(Heatplot(sig_res_for_ALL, showCategory = ShowCategory, foldChange = NULL), silent = T)
    }
  } else {
    if(df_with_lg2fc) {
      try(Heatplot(res$sig_res, showCategory = ShowCategory, foldChange = res$de), silent = T)
    } else {
      try(Heatplot(res$sig_res, showCategory = ShowCategory, foldChange = NULL), silent = T)
    }
  }
}


## dataframe for annotation, ORA, GSEA, PPI, DEG-RNAseq select organism----------
annoSpecies_df <<-
  data.frame(
    species = c(
      "", "Anopheles", "Arabidopsis", "Bovine", "Worm",
      "Canine", "Fly", "Zebrafish", "E coli strain K12",
      "E coli strain Sakai", "Chicken", "Human", "Mouse",
      "Rhesus", "Malaria", "Chimp", "Rat",
      "Yeast", "Pig",
      "Xenopus"
    ),
    pkg = c(
      "", "org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db",
      "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db",
      "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db",
      "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db", "org.Rn.eg.db",
      "org.Sc.sgd.db",  "org.Ss.eg.db", 
      "org.Xl.eg.db"
    ),
    stringsAsFactors = FALSE
  )
annoSpecies_df$organism = c("", "aga", "ath", "bta", "cel", "cfa", "dme", "dre", "eco", "ecs", "gga", "hsa", "mmu", "mcc", "pfa", "ptr", "rno", "sce", "ssc", "xla")
annoSpecies_df_for_reactome <- annoSpecies_df[c(12,17,13,5,18,8,7), ]
annoSpecies_df_for_reactome$organism = c("human", "rat", "mouse", "worm", "yeast", "zebrafish", "fly")
annoSpecies_df <- annoSpecies_df[order(annoSpecies_df$species), ]
rownames(annoSpecies_df) <- annoSpecies_df$species # easier to access afterwards


## givekegg_reat_res_and_table ##-----------
givekegg_reat_res_and_table <- function(reat, pCutoff = 0.05, p.adj.cutoff = 0.05, q.cutoff = 0.2){
  res <- reat$res
  sig_res <- clusterProfiler.dplyr::filter(res, pvalue < pCutoff, p.adjust < p.adj.cutoff, qvalue < q.cutoff)
  
  res_table <- as.data.frame(res)
  sig_res_table <- as.data.frame(sig_res)
  
  table <- list(all_table = res_table, sig_table = sig_res_table, sig_res = sig_res, de = reat[["de"]])
  
  return(table)
}

keggAnalysis <- function(df, organism="Human", df_with_lg2fc = FALSE, species_df){
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)
  
  organism <- species_df$organism[species_df$species == organism]
  
  # if(organism == "hsa"){
  #   library(org.Hs.eg.db)
  #   orgDB <<- org.Hs.eg.db
  #   #kegg_organism <- "hsa"
  # }
  # if(organism == "mmu"){
  #   library(org.Mm.eg.db)
  #   orgDB <<- org.Mm.eg.db
  #   #kegg_organism <- "mmu"
  # }
  # if(organism == "rno"){
  #   library(org.Rn.eg.db)
  #   orgDB <<- org.Rn.eg.db
  #   #kegg_organism <- "rno"
  # }
  
  
  # ids <- bitr(df$name, fromType = "SYMBOL",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB)
  # ids1 <- try(bitr(setdiff(df$name, ids$SYMBOL), fromType = "UNIPROT",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB))
  # if(class(ids1) != "try-error") {
  #   names(ids1)[1] = "SYMBOL"
  #   ids = rbind(ids, ids1)
  # }
  
  # names(ids)[1] = "name"
  #  ids <- inner_join(ids, df, by = "name")
  
  ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)
  
  ids <- inner_join(ids2, df, by = "name")
  
  if(df_with_lg2fc){
    ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
    de = ids$fc
    names(de) = unlist(ids$ENTREZID)
    de = sort(de, decreasing = T)
  } else {
    ids <- ids[!is.na(ids$ENTREZID), ]
    de = unlist(ids$ENTREZID)
    names(de) = de
  }
  
  reat <- enrichKEGG(gene = names(de), organism = organism,
                     pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  reat <- setReadable(reat, OrgDb = orgDB, keyType="ENTREZID")
  
  reat <- list(res = reat, de = de)
  return(reat)
}

#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism:  one of "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"
#must be form reat <- keggAnalysis(), because following use name reat
#species_df: data frame which have columns species, pkg and organism
reactAnalysis <- function(df, organism="Human", df_with_lg2fc = FALSE, species_df){
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)
  
  organism <- species_df$organism[species_df$species == organism]
  
  # if(organism == "human"){
  #   library(org.Hs.eg.db)
  #   orgDB <<- org.Hs.eg.db
  #   #kegg_organism <- "hsa"
  # }
  # if(organism == "mouse"){
  #   library(org.Mm.eg.db)
  #   orgDB <<- org.Mm.eg.db
  #   #kegg_organism <- "mmu"
  # }
  # if(organism == "rat"){
  #   library(org.Rn.eg.db)
  #   orgDB <<- org.Rn.eg.db
  #   #kegg_organism <- "rno"
  # }
  
  
  # ids <- bitr(df$name, fromType = "SYMBOL",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB)
  # ids1 <- try(bitr(setdiff(df$name, ids$SYMBOL), fromType = "UNIPROT",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB))
  # if(class(ids1) != "try-error") {
  #   names(ids1)[1] = "SYMBOL"
  #   ids = rbind(ids, ids1)
  # }
  
  # names(ids)[1] = "name"
  # ids <- inner_join(ids, df, by = "name")
  
  ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)
  
  ids <- inner_join(ids2, df, by = "name")
  
  if(df_with_lg2fc){
    ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
    de = ids$fc
    names(de) = unlist(ids$ENTREZID)
    de = sort(de, decreasing = T)
  } else {
    ids <- ids[!is.na(ids$ENTREZID), ]
    de = unlist(ids$ENTREZID)
    names(de) = de
  }
  
  reat <- ReactomePA::enrichPathway(gene = names(de), organism = organism,
                                    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable=TRUE)
  
  reat <- list(res = reat, de = de)
  return(reat)
}

# reat: keggAnalysis or reactAnalysis funtion return value
# pCutoff: pvalue cutoff
# p.adj.cutoff: padj cutoff
# q.cutoff: qvalue cutoff
givekegg_reat_res_and_table <- function(reat, pCutoff = 0.05, p.adj.cutoff = 0.05, q.cutoff = 0.2){
  res <- reat$res
  sig_res <- clusterProfiler.dplyr::filter(res, pvalue < pCutoff, p.adjust < p.adj.cutoff, qvalue < q.cutoff)
  
  res_table <- as.data.frame(res)
  sig_res_table <- as.data.frame(sig_res)
  
  table <- list(all_table = res_table, sig_table = sig_res_table, sig_res = sig_res, de = reat[["de"]])
  
  return(table)
}

#plot_type: barplot, dotplot,emaplot, cnetplot(based on your significant limit)
# res_kegg or res_react: givekegg_reat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# barplot(res_kegg$sig_res, showCategory = ShowCategory, color = color)
# barplot(res_react$sig_res, showCategory = ShowCategory, color = color)

# res_kegg or res_react: givekegg_reat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# dotplot(res_kegg$sig_res, showCategory = ShowCategory, color = color)
# dotplot(res_react$sig_res, showCategory = ShowCategory, color = color)

# res: givekegg_reat_res_and_table or giveGO_res_and_table funtion return value
# color: p.adjust" or "pvalue"
# ShowCategory: number of go terms of plot
my_dotplot_opt <- function(res, color = "p.adjust", size = "Count", title = "", decreasing = TRUE, ShowCategory = 20){
  df = as.data.frame(res$sig_res)
  
  if(nrow(df) > ShowCategory) {
    df = df[1:ShowCategory, ]
  }
  df$x = ""
  orderBy = color
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels = unique(df$Description[idx]))
  
  p_up = as.numeric(format(max(df[ , color]), digits = 2))
  p_low = as.numeric(format(min(df[ , color]), digits = 2))
  
  if(color == "p.adjust"){
    col.name = "FDR"
  }
  if(color == "pvalue"){
    col.name = "Pvalue"
  }
  
  
  a =  ggplot(df, aes_string(x = "x", y = "Description", size = size, color = color)) + 
    geom_point() + 
    geom_point(shape = 21, color = "black", stroke = 0.6, position = "identity")+
    scale_colour_gradient2(low = "#df293f", mid = ifelse(!max(df[ , color]) == min(df[ , color]), "#f9f9f9", "#f19c9b"), high = "#959cc8", midpoint = min(df[ , color]) + (max(df[ , color])-min(df[ , color]))/2, name = col.name, breaks = unique(c(min(df[ , color]), max(df[ , color]))), limits = c(min(df[ , color]), max(df[ , color])), labels = unique(c(format(p_low, scientific = TRUE), format(p_up, scientific = TRUE)))) + 
    ylab(NULL) + 
    ggtitle(title) + 
    scale_size(name = "Number of proteins found", range = c(3, 8), breaks = c(min(df[ , "Count"]), max(df[ , "Count"])), limits = c(min(df[ , "Count"]), max(df[ , "Count"]))) +
    theme(panel.background=element_blank(),
          axis.ticks=element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key=element_blank(),
          legend.spacing = unit(3, "lines"),
          legend.justification="bottom", legend.position=c(0.15,0),#0.15
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          #plot.margin = margin(t = 0, r = 0, b = 0.1, l = 15, unit = "cm"),
          title = element_text(size = 10)#all title elements: plot, axes, legends
    )+
    labs(x = NULL) +
    guides(
      color = guide_colorbar(order = 1, frame.colour = "black", frame.linetype = 1, frame.linewidth = 1, draw.ulim = TRUE, draw.llim = TRUE, raster = TRUE, ticks = FALSE, label = TRUE, title.position = "left", title.theme = element_text(angle = 90,  hjust = 0.5, size = 10), label.theme = element_text(angle = 90,  hjust = 0.5, size = 10), reverse = TRUE, barwidth = unit(0.5, units = "cm"), barheight = unit(2.5, units = "cm")),
      size = guide_legend(order = 0, title.position = "right", title.theme = element_text(angle = 90,  hjust = 0.5, size = 10), label.theme = element_text(size = 10), override.aes = list(alpha = 1, bg = "grey", color = "black"), reverse = TRUE, keyheight = 4)
    ) +
    scale_x_discrete(expand = expand_scale(mult = 0, add = c(0.008, 0.25)))+#0.005, 0.25
    scale_y_discrete()
  return(a)
}

# res: givekegg_reat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
emaplot_for_react_kegg <- function(res, ShowCategory = 30, color = "p.adjust", layout = "kk") {
  try(Emapplot(res$sig_res, showCategory = ShowCategory, color = color, layout = layout), silent = T) 
}

# res: givekegg_reat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# circular: if circular layout
# colorEdge: if colorEdge
# df_with_lg2fc: if df_with_lg2fc , the same as function goAnalysis
cnetplot_for_react_kegg <- function(res, ShowCategory = 5, circular = TRUE, colorEdge = TRUE, df_with_lg2fc = FALSE) {
  if(df_with_lg2fc) {
    try(Cnetplot(x = res$sig_res, showCategory = ShowCategory, foldChange = res$de, circular = circular, colorEdge = colorEdge), silent = T)
  } else {
    try(Cnetplot(x = res$sig_res, showCategory = ShowCategory, foldChange=NULL, circular = circular, colorEdge = colorEdge), silent = T)
  }
  
}

# res: givekegg_reat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# df_with_lg2fc: if df_with_lg2fc , the same as function goAnalysis
heatplot_for_react_kegg <- function(res, ShowCategory = 30, df_with_lg2fc = FALSE) {
  if(df_with_lg2fc) {
    try(Heatplot(res$sig_res, showCategory = ShowCategory, foldChange = res$de), silent = T)
  } else {
    try(Heatplot(res$sig_res, showCategory = ShowCategory, foldChange = NULL), silent = T)
  }   
}


## gsea part ##---------
gsegoAnalysis <- function(df, organism="Human", species_df){
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)
  
  
  # if(organism == "human"){
  #   library(org.Hs.eg.db)
  #   orgDB <<- org.Hs.eg.db
  #   #kegg_organism <- "hsa"
  # }
  # if(organism == "mouse"){
  #   library(org.Mm.eg.db)
  #   orgDB <<- org.Mm.eg.db
  #   #kegg_organism <- "mmu"
  # }
  # if(organism == "rat"){
  #   library(org.Rn.eg.db)
  #   orgDB <<- org.Rn.eg.db
  #   #kegg_organism <- "rno"
  # }
  
  
  # ids <- bitr(df$name, fromType = "SYMBOL",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB)
  # ids1 <- try(bitr(setdiff(df$name, ids$SYMBOL), fromType = "UNIPROT",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB))
  # if(class(ids1) != "try-error") {
  #   names(ids1)[1] = "SYMBOL"
  #   ids = rbind(ids, ids1)
  # }
  
  # names(ids)[1] = "name"
  # ids <- inner_join(ids, df, by = "name")
  
  ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)
  
  ids <- inner_join(ids2, df, by = "name")
  
  ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
  de = ids$fc
  names(de) = unlist(ids$ENTREZID)
  de = sort(de, decreasing = T)
  
  set.seed(1234)
  reat_ALL <- try(gseGO(gene = de, OrgDb = orgDB, ont = "ALL", 
                        pAdjustMethod = "BH", nPerm = 1000, pvalueCutoff = 1, verbose = FALSE, seed = FALSE), silent = TRUE)
  reat_ALL <- setReadable(reat_ALL, OrgDb = orgDB, keyType="ENTREZID")
  
  reat_BP <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "BP")
  reat_BP@setType = "BP"
  
  reat_CC <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "CC")
  reat_CC@setType = "CC"
  
  reat_MF <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "MF")
  reat_MF@setType = "MF"
  
  reat <- list(ALL = reat_ALL, BP = reat_BP, CC = reat_CC, MF = reat_MF, de = de)
  return(reat)
}

# reat: gsegoAnalysis funtion return value
# ont: one of ALL, BP, CC, MF
# pCutoff: pvalue cutoff
# p.adj.cutoff: padj cutoff
# NES.cutoff: NES cutoff,it means abs(NES) > NES.cutoff
# simplify: if remove redundancy of enriched GO terms
# Phenotype: the Phenotype you want to show, one or both of c("activated", "suppressed"), default: c("activated", "suppressed")
givegseGO_res_and_table <- function(reat, ont = "BP", pCutoff = 0.05, p.adj.cutoff = 0.25, NES.cutoff = 1, simplify = FALSE, Phenotype = c("activated", "suppressed")){
  res <- reat[[ont]]
  res@result$phenotype = ifelse(res@result$NES > 0, "activated", "suppressed")
  res@result$group = paste(res@result$ONTOLOGY, res@result$phenotype, sep = "_")
  
  if(length(Phenotype) == 1) {
    res <- clusterProfiler.dplyr::filter(res, phenotype == Phenotype)#
  }
  
  sig_res <- clusterProfiler.dplyr::filter(res, pvalue < pCutoff, p.adjust < p.adj.cutoff, abs(NES) > NES.cutoff)
  
  if(simplify) {
    if(ont == "ALL") {
      sig_res@setType = "BP"
      sig_res <- clusterProfiler::simplify(sig_res, cutoff=0.7, by="p.adjust", select_fun=min)
    } else {
      sig_res <- clusterProfiler::simplify(sig_res, cutoff=0.7, by="p.adjust", select_fun=min)
    }
  }
  
  res_table <- as.data.frame(res)
  sig_res_table <- as.data.frame(sig_res)
  
  sig_res@result = sig_res@result[ , c(2:ncol(sig_res@result), 1)]
  
  table <- list(all_table = res_table, sig_table = sig_res_table, sig_res = sig_res, de = reat[["de"]])
  
  return(table)
}

#plot_type: barplot, dotplot,emaplot, cnetplot, heatplot, gseaplot2

# res: givegseGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# Split: act when ont == "ALL", if split by ontology , TRUE or FALSE
gse_dotplot <- function(res, ShowCategory = 20, color = "p.adjust", ont = "BP", Split = FALSE) {
  library(clusterProfiler.dplyr)
  library(ggplot2)
  if(ont == "ALL") {
    if(Split) {
      enrichplot:::dotplot(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split = "group") + facet_grid(ONTOLOGY ~ ., scale="free", space = "free")
    } else {
      enrichplot:::dotplot(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split = "phenotype") + facet_grid(ONTOLOGY ~ ., scale="free", space = "free")
    }
  } else {
    enrichplot:::dotplot(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split="group")
  }
}

# res: givegseGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# Split: act when ont == "ALL", if split by ontology , TRUE or FALSE
gse_barplot <- function(res, ShowCategory = 20, color = "p.adjust", ont = "BP", Split = FALSE) {
  library(clusterProfiler.dplyr)
  library(ggplot2)
  if(ont == "ALL") {
    if(Split) {
      enrichplot:::barplot.enrichResult(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split = "group") + facet_grid(ONTOLOGY ~ ., scale="free", space = "free") + labs(y = "NES")
    } else {
      enrichplot:::barplot.enrichResult(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split = "phenotype") + facet_grid(ONTOLOGY ~ ., scale="free", space = "free") + labs(y = "NES")
    }
  } else {
    enrichplot:::barplot.enrichResult(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split="group") + labs(y = "NES")
  }
}

# res: giveGO_res_and_table or givegseGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# try(emapplot(res$sig_res, showCategory = ShowCategory, color = color, layout = "kk"), silent = T)


# res: giveGO_res_and_table or givegseGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# foldChange: res$de
# circular: TRUE or FALSE
# try(cnetplot(res$sig_res, showCategory = ShowCategory, color = color, foldChange = res$de, circular = circular, colorEdge = TRUE), silent = T)

# ShowCategory: number of go terms of plot
# foldChange: res$de
# try(heatplot(res$sig_res, showCategory = ShowCategory, foldChange = res$de), silent = T)

#enrichplot:::tableGrob2 修改，将ptable的color去掉,使description的颜色由绿色变为黑色，仅适用于ES有一条线
mytableGrob2 <- function (d, p = NULL) 
{
  d <- d[order(rownames(d)), ]
  tp <- gridExtra::tableGrob(d)#, theme = gridExtra::ttheme_default(base_size = 9),theme = gridExtra::ttheme_minimal()
  if (is.null(p)) {
    return(tp)
  }
  pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
  j <- which(tp$layout$name == "rowhead-fg")
  for (i in seq_along(pcol)) {
    tp$grobs[j][[i + 1]][["gp"]] = grid::gpar(col = "black", cex = 0.8)
  }
  return(tp)
}


my_gseaplot2 <- function (x, geneSetID, title = "", color = "green", base_size = 11, 
                          rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
                          ES_geom = "line") 
{
  library(ggplot2)
  library(dplyr)
  # library(openxlsx)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if (length(geneSetID) == 1) {
    gsdata <- enrichplot:::gsInfo(x, geneSetID)
  }
  else {
    gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
  }
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0))
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
                          size = 1)
  }
  else {
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
                           size = 1, data = subset(gsdata, position == 1))
  }
  p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8), 
                                legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
  p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(), 
                                                            axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                                                            plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
                                                                                 unit = "cm"))
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == 
                   term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
                                                           ymax = ~ymax, color = ~Description)) + xlab(NULL) + 
    ylab(NULL) + theme_classic(base_size) + theme(legend.position = "none", 
                                                  plot.margin = margin(t = -0.1, b = 0, unit = "cm"), 
                                                  axis.ticks = element_blank(), axis.text = element_blank(), 
                                                  axis.line.x = element_blank()) + scale_x_continuous(expand = c(0, 
                                                                                                                 0)) + scale_y_continuous(expand = c(0, 0))
  if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) 
      inv <- inv + 1
    col = c(rev(RColorBrewer::brewer.pal(5, "Blues")), RColorBrewer::brewer.pal(5, "Reds"))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[unique(inv)])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                    xmax = xmax, col = col[unique(inv)])
    p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
                              ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
                         alpha = 0.9, inherit.aes = FALSE)
  }
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
                                             y = ~y, yend = 0), color = "grey")
  p.pos <- p.pos + ylab("Ranked list metric") + xlab("Rank in Ordered Dataset") + 
    theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, 
                               l = 0.2, unit = "cm"))
  if (!is.null(title) && !is.na(title) && title != "") 
    p.res <- p.res + ggtitle(title)
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust", "NES")]
    pd <- pd[order(pd[, 1], decreasing = FALSE), ]
    rownames(pd) <-  "" 
    
    pd <- pd[, -1]
    pd <- round(pd, 4)
    tp <- mytableGrob2(pd, p.res) 
    p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
                                                                         xmin = quantile(p.res$data$x, 0.7), xmax = quantile(p.res$data$x, 
                                                                                                                             0.95), ymin = quantile(p.res$data$runningScore, 
                                                                                                                                                    0.75), ymax = quantile(p.res$data$runningScore, 
                                                                                                                                                                           0.9))
  }
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                         axis.ticks.x = element_line(), axis.text.x = element_text())
  if (length(subplots) == 1) 
    return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
                                                      r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
  if (length(rel_heights) > length(subplots)) 
    rel_heights <- rel_heights[subplots]
  PloT <- plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
  PloT + theme(plot.margin = margin(t = 1, 
                                    r = 1, b = 0.2, l = 0.2, unit = "cm"))
}

gsekeggAnalysis <- function(df, organism="Human", species_df){
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)
  
  organism <- species_df$organism[species_df$species == organism]
  
  
  # if(organism == "hsa"){
  #   library(org.Hs.eg.db)
  #   orgDB <<- org.Hs.eg.db
  #   #kegg_organism <- "hsa"
  # }
  # if(organism == "mmu"){
  #   library(org.Mm.eg.db)
  #   orgDB <<- org.Mm.eg.db
  #   #kegg_organism <- "mmu"
  # }
  # if(organism == "rno"){
  #   library(org.Rn.eg.db)
  #   orgDB <<- org.Rn.eg.db
  #   #kegg_organism <- "rno"
  # }
  
  
  # ids <- bitr(df$name, fromType = "SYMBOL",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB)
  # ids1 <- try(bitr(setdiff(df$name, ids$SYMBOL), fromType = "UNIPROT",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB))
  # if(class(ids1) != "try-error") {
  #   names(ids1)[1] = "SYMBOL"
  #   ids = rbind(ids, ids1)
  # }
  
  # names(ids)[1] = "name"
  # ids <- inner_join(ids, df, by = "name")
  
  #   de = ids$fc
  #   names(de) = ids$ENTREZID
  #   de = sort(de, decreasing = T)
  
  ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)
  
  ids <- inner_join(ids2, df, by = "name")
  
  ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
  de = ids$fc
  names(de) = unlist(ids$ENTREZID)
  de = sort(de, decreasing = T)    
  
  set.seed(1234)
  reat <- try(gseKEGG(gene = de, organism = organism,
                      pAdjustMethod = "BH", nPerm = 1000, pvalueCutoff = 1, verbose = FALSE, seed = FALSE), silent = TRUE)
  reat <- setReadable(reat, OrgDb = orgDB, keyType="ENTREZID")
  
  reat <- list(res = reat, de = de)
  return(reat)
}

#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism:  one of "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"
#must be form reat <- gsereactAnalysis(), because following use name reat
#species_df: data frame which have columns species, pkg and organism
gsereactAnalysis <- function(df, organism="Human", species_df){
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)
  
  organism <- species_df$organism[species_df$species == organism]
  
  # if(organism == "human"){
  #   library(org.Hs.eg.db)
  #   orgDB <<- org.Hs.eg.db
  #   #kegg_organism <- "hsa"
  # }
  # if(organism == "mouse"){
  #   library(org.Mm.eg.db)
  #   orgDB <<- org.Mm.eg.db
  #   #kegg_organism <- "mmu"
  # }
  # if(organism == "rat"){
  #   library(org.Rn.eg.db)
  #   orgDB <<- org.Rn.eg.db
  #   #kegg_organism <- "rno"
  # }
  
  
  # ids <- bitr(df$name, fromType = "SYMBOL",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB)
  # ids1 <- try(bitr(setdiff(df$name, ids$SYMBOL), fromType = "UNIPROT",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB))
  # if(class(ids1) != "try-error") {
  #   names(ids1)[1] = "SYMBOL"
  #   ids = rbind(ids, ids1)
  # }
  
  # names(ids)[1] = "name"
  # ids <- inner_join(ids, df, by = "name")
  
  #   de = ids$fc
  #   names(de) = ids$ENTREZID
  #   de = sort(de, decreasing = T)
  
  ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)
  
  ids <- inner_join(ids2, df, by = "name")
  
  ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
  de = ids$fc
  names(de) = unlist(ids$ENTREZID)
  de = sort(de, decreasing = T)    
  
  set.seed(1234)
  reat <- try(ReactomePA::gsePathway(gene = de, organism = organism,
                                     pAdjustMethod = "BH", nPerm = 1000, pvalueCutoff = 1, verbose = FALSE, seed = FALSE), silent = TRUE)
  reat <- setReadable(reat, OrgDb = orgDB, keyType="ENTREZID")
  
  reat <- list(res = reat, de = de)
  return(reat)
}


# reat: gsekeggAnalysis or gsereactAnalysis funtion return value
# pCutoff: pvalue cutoff
# p.adj.cutoff: padj cutoff
# NES.cutoff: NES cutoff,it means abs(NES) > NES.cutoff
# Phenotype: the Phenotype you want to show, one or both of c("activated", "suppressed"), default: c("activated", "suppressed")
give_gsekegg_gsereat_res_and_table <- function(reat, pCutoff = 0.05, p.adj.cutoff = 0.25, NES.cutoff = 1, Phenotype = c("activated", "suppressed")){
  res <- reat$res
  res@result$phenotype = ifelse(res@result$NES > 0, "activated", "suppressed")
  
  if(length(Phenotype) == 1) {
    res <- clusterProfiler.dplyr::filter(res, phenotype == Phenotype)#
  }
  
  sig_res <- clusterProfiler.dplyr::filter(res, pvalue < pCutoff, p.adjust < p.adj.cutoff, abs(NES) > NES.cutoff)
  
  res_table <- as.data.frame(res)
  sig_res_table <- as.data.frame(sig_res)
  
  table <- list(all_table = res_table, sig_table = sig_res_table, sig_res = sig_res, de = reat[["de"]])
  
  return(table)
}