
################  Using CorShrink per gene correlation matrix shrinkage  #########################

library(data.table)
data <- data.frame(fread("../data/GTEX_V6/cis_gene_expression.txt"))
matdata <- data[,-c(1,2)]
voom_matdata <- t(limma::voom(matdata))
tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]
samples_id <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,1]
Numgenes <- 100;


person_label=read.table("../data/GTEX_V6/person_identifier_labels_with_numbers.txt");
person_label_num=as.numeric(person_label[,2]);


U=unique(tissue_labels);

cor_tissue1_tissue2_genes=array(0,c(length(U),length(U),Numgenes));

num_samples_common=matrix(0,length(U),length(U));

for(l in 1:(length(U)-1))
{
  for(k in (l+1):length(U))
  {
    index1=which(!is.na(match(tissue_labels, U[l]))); ## which samples come from tissue U[l]
    index2=which(!is.na(match(tissue_labels, U[k]))); ## which samples come from tissue U[k]
    person_label1=person_label_num[index1]; ## person labels for these samples in U[l]
    person_label2=person_label_num[index2]; ## person labels for these samples in U[k]
    person_label_pooled <- intersect(person_label1, person_label2);

    num_samples_common[l,k]=length(person_label_pooled);
    num_samples_common[k,l]=num_samples_common[l,k];

    if(length(person_label_pooled)==0)
    {
      cor_tissue1_tissue2_genes[l,k,]=rep(NA,Numgenes);
      cor_tissue1_tissue2_genes[k,l,]=rep(NA,Numgenes);
    }
    if(length(person_label_pooled)!=0)
    {
      index_1_main=index1[which(!is.na(match(person_label1,person_label_pooled)))]
      index_2_main=index2[which(!is.na(match(person_label2,person_label_pooled)))]

      temp1=matrix(as.numeric(as.character(as.matrix(as.vector(voom_matdata[index_1_main, 1:Numgenes])))),nrow=length(index_1_main))
      temp2=matrix(as.numeric(as.character(as.matrix(as.vector(voom_matdata[index_2_main, 1:Numgenes])))),nrow=length(index_2_main))

      cor_result=diag(cor(temp1,temp2,method="spearman"))
      cor_tissue1_tissue2_genes[l,k,]=cor_result;
      cor_tissue1_tissue2_genes[k,l,]=cor_result;
    }
  }
  cat("we are at tissue ", l, "\n");
}

save(cor_tissue1_tissue2_genes, file = "../rda/cor_tissues_non_ash_voom.rda")
save(num_samples_common, file = "../rda/common_samples_voom.rda")

