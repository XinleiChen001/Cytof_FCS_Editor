
###############################################################################
#                                                                             #
#         FCS文件通道修改，   版本号：2.001     发布时间：2019-12-27          #    
#                                                                             #
###############################################################################


library(flowCore)
library(dplyr)


#选择目标文件目录

wdir <-"D:/Customer Data/2019_12_24_Channel_Del/FCS_Edit_2019_9_22/1_FCS_Edit"   #工作目录
rawdir<-"D:/Customer Data/2019_12_24_Channel_Del/FCS_Edit_2019_9_22/0_raw_data"   #原始文件目录


setwd(wdir)
source("./backup/FCS_Edit_backup.R")  #载入FCS_Edit.R



#获取工作目录中所有FCS文件名称，注意文件夹中的FCS文件应该是完全相同的Panel
files <- list.files(path=rawdir,pattern='.fcs$', full=TRUE)

files

#FCS文件读取为flowset数据类型
raw_flowset=suppressWarnings(read.flowSet(files[1],transformation = FALSE))


#输出marker列表

regen_file = TRUE
#regen_file = FALSE

if (regen_file) {
markers_edit_raw<-raw_flowset[[1]]@parameters@data[,c(1,2)]
markers_edit_raw<-cbind(markers_edit_raw,matrix(rep("",5*nrow(markers_edit_raw)),ncol = 5))
colnames(markers_edit_raw)<-c("name","desc","original","values_to_add","rename","values_to_replace","markers_to_del")
markers_edit_raw$original=1

write.csv(markers_edit_raw,"markers_edit.csv",row.names = F)

}
#regen..后面的这段代码仅需运行一次，用来输出marker列表，运行后regen_file要设置为FALSE，以免覆盖。

#用Excel打开markers_edit.csv,按照说明设置各列内容，保存退出

markers_edit<-read.csv("markers_edit.csv",stringsAsFactors = F)

markers_edit
FCS_Edit(raw_fcs_dir = rawdir,
         output_fcs_dir = paste0(wdir,"/output") ,
         markers_edit=markers_edit)



