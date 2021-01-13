
#版本信息：v2019-12-27
#修复了rep后文件无法读取的bug

##FCS_Edit是可以用来删除FCS的通道或者更改通道名称，三个参数第一个为原始FCS文件所在的位置，
##第二个为输出FCS文件的位置，默认值是在原始FCS目录下建立Output文件夹
##第三个为markers_edit_csv，包含删除和更名的信息。




FCS_Edit<-function( raw_fcs_dir,output_fcs_dir=NULL,markers_edit=NULL){


  markers_edit[,c(4,5,6)][is.na(markers_edit[,c(4,5,6)])]<-""
  

  if(is.null(output_fcs_dir)){

    output_fcs_dir<-paste0(raw_fcs_dir,"/output")
  }


  files <- list.files(path=raw_fcs_dir,pattern='.fcs$', full=TRUE)

  raw_flowset=suppressWarnings(read.flowSet(files,transformation = FALSE))

  if(sum(is.na(markers_edit[,"original"]))>0){
  marker_original=filter(markers_edit,original==1)} else
  marker_original=markers_edit
  
  
  if(nrow(parameters(raw_flowset[[1]]))!=nrow(marker_original)){

    cat("Error:FCS parameters is not compatible with csv file, please select the right csv. \n")
    return(FALSE)
  }



  rename_para_id<-which(marker_original$rename!="")

  del_para_id<-which(markers_edit$markers_to_del==1)
  del_para_id_ori<-which(marker_original$markers_to_del==1)
  add_para_id<-which(markers_edit$values_to_add!="")
  rep_para_id<-which(markers_edit$values_to_replace!="")


for(i in 1:length(raw_flowset)){

    raw_file_name<-raw_flowset@phenoData$name[i]


    raw_flowframe <- raw_flowset[[i]]

    exprs <- raw_flowframe@exprs
    params <- parameters(raw_flowframe)
    pd <- pData(params)
    keyval <- keyword(raw_flowframe)

   # str(exprs)

    #####change marker name

    fcs_desc<-pd$desc
    fcs_names<-pd$name

    for(m in rename_para_id){
      #m=rename_para_id[2]
      fcs_desc[[m]]<-as.character(marker_original[m,5])
      if(fcs_names[m]==as.character(marker_original[m,2])){
        fcs_names[m]<-as.character(marker_original[m,5])
        }
    }
    pd$desc<-fcs_desc



    for(m in rename_para_id){
#      m=2
      key_rename_id<-NULL
      for(n in 1:length(keyval)){

      if(keyval[[n]]==as.character(marker_original[m,2])) {key_rename_id<-c(key_rename_id,n)}

      }

      if(!is.null(key_rename_id)){
      for(o in key_rename_id){
      keyval[[o]]<-as.character(marker_original[m,5])}
      }
     }




    ####add channel  and del channel

    exprs_add<-exprs
    Exprs<-as.data.frame(exprs)

    #增加通道数据
    
    if (!is.na(add_para_id[1])){
    
    for(i_add in c(1:length(add_para_id))){
    #i_add=1
    
      indata<-eval(parse(text = as.character(markers_edit[add_para_id[i_add],4])))
      add_data<-matrix(data=indata,ncol=1)
    colnames(add_data)<-markers_edit[add_para_id[i_add],1]
    
    
    exprs_add<-cbind(exprs_add[,c(1:(add_para_id[i_add]-1)),drop=F],
                     add_data,
                     exprs_add[,c(add_para_id[i_add]:ncol(exprs_add)),drop=F]
                     )
                     
    }
    }
    

    Exprs<-as.data.frame(exprs_add)
    #replace通道数据
    if (!is.na(rep_para_id[1])){
      for(i_rep in c(1:length(rep_para_id))){
        #i_rep=2
        indata<-eval(parse(text = as.character(markers_edit[rep_para_id[i_rep],6])))
        add_data<-matrix(data=indata,ncol=1)
        exprs_add[,rep_para_id[i_rep]]<-indata
      }
    }
    

    
    if (!is.na(del_para_id[1])){
    exprs_del<-exprs_add[,-1*del_para_id]}else
    exprs_del<-exprs_add

    
    keyval_add<-keyval
    
    pd_add<-pd
    
    if (!is.na(add_para_id[1])){
    
    
    for (j_add in 1:length(add_para_id)) {
      #j_add=1
      channel_name <- as.character(markers_edit[add_para_id[j_add],1])
      channel_desc <- as.character(markers_edit[add_para_id[j_add],2])

      minRange <- ceiling(min(exprs_add[,add_para_id[j_add]]))
      maxRange <- ceiling(max(exprs_add[,add_para_id[j_add]]))
      channel_range <- maxRange - minRange+1  #   修改一个计算错误

      plist <- matrix(c(channel_name, channel_desc, channel_range, 
                        minRange, maxRange))
      
      
      rownames(plist) <- c("name", "desc", "range", "minRange", 
                           "maxRange")
      colnames(plist) <- paste0("$Padd_",j_add)

      pd_add <- rbind(pd_add[c(1:(add_para_id[j_add]-1)),], 
                      t(plist),
                      pd_add[c(add_para_id[j_add]:nrow(pd_add)),])
      
      keyval_add[[paste("$P", "add_",j_add, "B", sep = "")]] <- "32"
      keyval_add[[paste("$P", "add_",j_add, "R", sep = "")]] <- toString(channel_range)
      keyval_add[[paste("$P", "add_",j_add, "E", sep = "")]] <- "0,0"
      keyval_add[[paste("$P", "add_",j_add, "N", sep = "")]] <- channel_name
      keyval_add[[paste("$P", "add_",j_add, "S", sep = "")]] <- channel_desc  ###added to show shortname of parameters

    }
      
    }
    
    #修改replace通道的metadata
    if (!is.na(rep_para_id[1])){
      
      
      for (j_rep in 1:length(rep_para_id)) {

        minRange <- ceiling(min(exprs_add[,rep_para_id[j_rep]]))
        maxRange <- ceiling(max(exprs_add[,rep_para_id[j_rep]]))
        channel_range <- maxRange - minRange+1  #   修改一个计算错误

        
        channel_name <- as.character(markers_edit[rep_para_id[j_rep],1])   
        channel_desc <- as.character(markers_edit[rep_para_id[j_rep],2])
        
        plist <- matrix(c(channel_name, channel_desc, channel_range, 
                          minRange, maxRange))
        
        
        pd_add$minRange<-minRange
        pd_add$maxRange<-maxRange
        pd_add$range<-channel_range
        
        keyval_add[[paste(rownames(pd_add[rep_para_id[j_rep],]), "R", sep = "")]] <- toString(channel_range)   #修复一处bug

        
      }
      
    }
    
    
    
    if (!is.na(del_para_id[1])) 
      {pd_del<-pd_add[-1*del_para_id,]} else
       pd_del<-pd_add
    
    pData(params)<- pd_del



    pattern_del<-paste0("\\$P",del_para_id_ori,"[A-Z]")
    key_del_id=NULL
    for(j in 1:length(pattern_del)){
      key_del_id <- c(key_del_id,grep(pattern = pattern_del[j],names(keyval)))
    }

    keyval_del<-keyval_add[-1*key_del_id]



    pd_del_rownames<-row.names(pd_del)
    row.names(pd_del)<-paste0("$P",c(1:length(pd_del_rownames)))

    pd_del_rownames<-cbind(pd_del_rownames,paste0(row.names(pd_del),"L"),row.names(pd_del))
    
    names(colnames(exprs_del)) <- paste0(row.names(pd_del),"N")

    key_del_name<-names(keyval_del)
    params_del<-params
    pData(params_del)<- pd_del

    for(k in 1:nrow(pd_del_rownames)){
      #k=11
      pattern_rename<-paste0("\\",pd_del_rownames[k,1],"[A-Z]","$")
      key_rename_id<-grep(pattern = pattern_rename,key_del_name)
      key_del_name[key_rename_id]<-sub(pd_del_rownames[k,1],pd_del_rownames[k,2],key_del_name[key_rename_id],fixed = T)

    }
    for(k in 1:nrow(pd_del_rownames)){
      #k=11
      pattern_rename<-paste0("\\",pd_del_rownames[k,2],"[A-Z]","$")
      key_rename_id<-grep(pattern = pattern_rename,key_del_name)
      key_del_name[key_rename_id]<-sub(pd_del_rownames[k,2],pd_del_rownames[k,3],key_del_name[key_rename_id],fixed = T)
      
    }
    
    names(keyval_del)<-key_del_name

    #####Write fcs files:

    if (!dir.exists(output_fcs_dir)) {
      dir.create(output_fcs_dir)
    }


    cat("Save to file:",paste0(output_fcs_dir, "/",raw_file_name, "\n"))

    
     
    
 out_frame <- flowFrame(exprs = exprs_del, parameters = params_del,
                        description = keyval_del) 

 
 head(exprs_del)
 


suppressWarnings(write.FCS(out_frame, paste0(output_fcs_dir, "/",raw_file_name)))


  }

  return(TRUE)

}
