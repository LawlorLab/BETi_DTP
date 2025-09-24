
samplesheet <- read.table("/active/taylor_s/people/nkatiy/CP-E_sarcoma_BETi/lawlorlab_shireen_2023.07.12_bulk_rnaseq_count_nf/resources/paired_end_samplesheet_A673_3days.csv", sep=",", header=T)
#print(samplesheet)

samplesheet <- samplesheet[samplesheet$cell_line=="A673",]
#print(samplesheet)

files_len <- nrow(samplesheet)
#print(files_len)
#files = list.files(path = "/active/taylor_s/people/nkatiy/CP-CDK8_RNA-seq/evansii_m_2023.01.23_CDK8_bulk_RNAseq/CDK8_bulk_RNAseq_count_nf/paired_end_results/star/", pattern = "ReadsPerGene", full.names=TRUE)

filepath = "/active/taylor_s/people/nkatiy/CP-E_sarcoma_BETi/lawlorlab_shireen_2023.07.12_bulk_rnaseq_count_nf/paired_end_results/star/"
file <- paste(filepath,samplesheet$id[1],".ReadsPerGene.out.tab", sep="")

count_mat <- read.table(file, skip = 4, sep="\t", header=F, row.names=1)
#print(head(count_mat))
count_mat1 <- count_mat[,1]

for(i in 2:files_len){
	file <- paste(filepath,samplesheet$id[i],".ReadsPerGene.out.tab", sep="")
	#print(file)
	count_dat <- read.table(file, skip = 4, sep="\t", header=F, row.names=1)
	count_mat1 <- cbind(count_mat1, count_dat[,1])
	rownames(count_mat1) <- row.names(count_dat)
	}

colnames(count_mat1) <- samplesheet$id

write.table(count_mat1, "additional_results/Read_count_matrix_cell_line_A673_3days.txt", sep="\t", quote=F)

