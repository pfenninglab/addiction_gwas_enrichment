ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

# make ldsc regression ldsc file
dir = '/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/annotations/Fullard_CNNscored_EUR2'
pd = read.table("Fullard_cnn_scored_peaks.txt")

files = list.files(path = dir, full.names = T)
files = unique(ss(files,'\\.'))
bgFile = grep('HoneyBadger_Fullard_CNN_scored',files, value = T)
lab = gsub('_cnnV3_scored_Fullard_peaks.bed.gz|Ctx\\.|Str\\.', '',basename(pd$V2))
lab = gsub('cutOff','',lab)
files = file.path(dir, lab)
Celltypes = basename(files)
files = paste0(files,'.,',bgFile,'.')

ldctsFile = 'Fullard_cnn_scored_peaks.ldcts'
write.table(data.frame(Celltypes,files), file = ldctsFile, sep = '\t',quote = F, row.names = F, col.names = F)

# 
# ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
# options(stringsAsFactors = F)
# 
# # make ldsc regression ldsc file
# dir = '/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/annotations/Fullard_SVMscored_EUR'
# pd = read.table("Fullard_svm_scored_peaks.txt")
# 
# files = list.files(path = dir, full.names = T)
# files = unique(ss(files,'\\.'))
# bgFile = grep('HoneyBadger_Fullard_SVM_scored',files, value = T)
# files = files[basename(files) %in% pd$V1]
# Celltypes = basename(files)
# files = paste0(files,'.,',bgFile,'.')
# 
# ldctsFile = 'Fullard_svm_scored_peaks.ldcts'
# write.table(data.frame(Celltypes,files), file = ldctsFile, sep = '\t',quote = F, row.names = F, col.names = F)
# 
# 
# 
# 
# ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
# options(stringsAsFactors = F)
# 
# # make ldsc regression ldsc file
# dir = '/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/annotations/Fullard_BOTHscored_EUR'
# pd = read.table("Fullard_both_scored_peaks.txt")
# 
# files = list.files(path = dir, full.names = T)
# files = unique(ss(files,'\\.'))
# bgFile = grep('HoneyBadger_Fullard_both_scored',files, value = T)
# files = files[basename(files) %in% pd$V1]
# Celltypes = basename(files)
# files = paste0(files,'.,',bgFile,'.')
# 
# ldctsFile = 'Fullard_both_scored_peaks.ldcts'
# write.table(data.frame(Celltypes,files), file = ldctsFile, sep = '\t',quote = F, row.names = F, col.names = F)
# 
# 
