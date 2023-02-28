library(BiocManager)
library(tximport)
library(readr)
library(DESeq2)
library(edgeR)
#BiocManager::install("edgeR")
samples <- read_tsv("H:/Downloads/06c4bb74-2223-4107-b34b-92478f1422fa.kallisto.abundance.tsv.gz")

#  "H:/Downloads/06c4bb74-2223-4107-b34b-92478f1422fa.kallisto.abundance.tsv.gz", header = TRUE)
##########################################################
# Below is code to make tx2gene function
#
##########################################################

library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

edb <- EnsDb.Hsapiens.v86
txdb <- transcripts(edb, 
                    columns = c(listColumns(edb, "tx"), c(listColumns(edb, "gene"))),
                    return.type = "DataFrame")
tx2gene <- tibble("tx_id" = txdb@listData[["tx_id"]], 
                  "gene_id" = txdb@listData[["gene_name"]])


txi_merged = tximport(files = c("/Users/joshin2/Downloads/f59ac261-5f1b-474c-848a-e788e923fbdf.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/f0cb8ab8-17cf-44c7-abf9-bb9dcd72f734.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/e4fecc0d-5af1-4e6b-92f7-6e699af08305.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/ccbdb608-4bfb-45d8-92fb-e5dba5dc93f0.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/c31620a9-e5c3-4dfd-8b32-1244ddf08c39.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/ac764709-65d5-4a9f-9fb6-382bc3b6c60e.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/7b929701-3f57-4d6b-aaac-b5de3d2bf3a8.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/7a96d448-ab36-4097-b1c9-afd1dec42d22.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/753650c1-34bc-4b78-afae-9494692d7b0b.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/729b560b-1bc7-4127-b099-60b5fb58760c.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/5abc4934-1d10-4eb1-94cf-2b7e08d6f3a6.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/69791b1f-b3eb-4775-b416-32f9960be617.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/49d2fd08-b8dd-4f0d-8447-6dbd1bfc3f2b.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/62795d35-9e9f-408e-8a1d-05e2861d95c3.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/45b3cb32-d0be-4e7e-8fdb-b48881d53488.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/33716491-6ea9-4ae3-bbb5-5aec9c822cff.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/2df28126-3e66-4ba6-8b96-a13d8607f839.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/2d66286a-bb58-4bc0-b89c-d865e83c693a.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/2bf7928c-4fab-48c5-8e44-0b59baaad7e0.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/293ee9ad-bece-49c0-9561-2dbd33e44fb8.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/25af94d7-c03f-4506-a7ed-9aee0025fca1.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/226cac03-0544-471a-8808-5b4c43f0727d.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/1d892c75-1fa2-431c-b7e2-a95d1c839f37.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/1c35a73f-23fb-4212-b5c7-1f59487bf0d6.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/1b42aa58-6e2c-4bf6-b885-551a5652b281.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/139a5e81-bd60-45f5-84d2-75f622e449b7.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/0afff214-e585-4cec-9501-1456d3836ac7.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/09764b99-fab0-4079-a780-081710698034.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/ffb9694e-6a53-47ab-8204-75ca08c93f7a.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/fc0c9795-03e4-45fe-a69f-46d20cd6e10f.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/e6535989-48cd-40e0-978c-b910d72b82ff.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/e58e2312-353c-45c1-a946-7d8cc6500efc.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/d4fa5be8-26e9-4a0a-a665-6bfd12c92b11.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/cfde8c98-f5dc-426b-8e4d-7c8387fcadb6.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/c4a48436-5548-4d54-9f86-677636f4a930.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/c3b96c55-f6da-4032-9c46-9b8f96aae754.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/b4365569-74e2-412b-91a4-186c4ef73c7e.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/b07bcad0-0523-4507-8345-93ba17eee5e0.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/afd583f5-6b9d-4166-b7aa-6eaa4460da0a.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/8d19e7e3-dd0a-4974-ae0d-ef8a6e68d4c9.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/88cb44da-d0d0-4ec0-a767-bbb297cb8621.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/8540d8ad-cb4d-4878-9311-ba472366cf64.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/75769edc-f168-455f-ac5e-7a9f8090757d.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/636274ce-d529-4eb8-81e6-8d262f0bc30a.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/62783acd-e7aa-4d6f-b018-71240622f3c0.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/41b304fd-1a41-44d5-bcca-07df95861995.kallisto.abundance.tsv.gz", 
                                "/Users/joshin2/Downloads/4021abc2-b0dc-412c-8416-db5937f67f27.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/16c71be5-b6d8-4382-b3d2-bffe61441de4.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/1500a8bf-4a78-421f-a1d0-c4445473a16f.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/1052dec8-e735-4bfc-8c72-8abd4d742121.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/0b8f07e3-747e-48d3-b172-629e2be36b77.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/07a1f9b2-5988-4d0b-a903-42c1f8d652da.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/063395ae-879e-43c1-b700-752aee34f093.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/e10a2795-ed4d-4a9a-aa0c-66f28282f1d7.kallisto.abundance.tsv.gz", #recur here start!
                                "/Users/joshin2/Downloads/d9c9477a-431e-4e41-9788-563b83048ba0.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/d64c0bf9-934b-4287-930e-1d7ca6c8ceb9.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/9f5cf8a1-257e-48de-a062-2f6e60a66bf6.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/6bfbfc4a-d280-476f-b0bc-059ef308943d.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/3c837573-37d1-42f5-b3f2-96911536145a.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/279f08fc-9521-409c-af63-bc39316703df.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/171e5ed7-9e88-46bb-bd61-a6fe49caa80e.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/0db0af2a-2792-4944-a636-a329190fbf8a.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/0a4cec33-95b8-4f47-875f-8b4674c7f1ff.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/080d0f69-5133-4cd1-9e58-8b7a7d4d36b3.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/f9dc9c3b-2593-4728-a05a-2f4fdd454249.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/dcdc79e8-e061-498b-81af-16eb1d1d4d09.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/cf80fd30-0f27-4621-a150-425055fc89c2.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/c6c464fc-901c-4199-bd39-2b7b13612fbd.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/9fb8cfc1-cb1b-4ffe-8fe4-5e1c808e3b1b.kallisto.abundance.tsv.gz",
                                "/Users/joshin2/Downloads/50b1a57c-d95e-463e-9146-97921d8f0846.kallisto.abundance.tsv.gz"),
                      type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
library(tidyverse)
sample_id <- tibble(sample = c(1:70), condition = factor(c("Initial", "Initial", "Initial", 
                                                           "Initial", "Initial", "Initial", 
                                                           "Initial", "Initial", "Initial", 
                                                           "Initial", "Initial", "Initial",
                                                           "Initial", "Initial", "Initial", 
                                                           "Initial", "Initial", "Initial", 
                                                           "Initial", "Initial", "Initial", 
                                                           "Initial", "Initial", "Initial",
                                                           "Initial", "Initial", "Initial",
                                                           "Initial", "Initial", "Initial", 
                                                           "Initial", "Initial", "Initial", 
                                                           "Initial", "Initial", "Initial", 
                                                           "Initial", "Initial", "Initial",
                                                           "Initial", "Initial", "Initial", 
                                                           "Initial", "Initial", "Initial", 
                                                           "Initial", "Initial", "Initial", 
                                                           "Initial", "Initial", "Initial",
                                                           "Initial", "Initial",
                                                           "Rec/Prog", 
                                                           "Rec/Prog", "Rec/Prog", "Rec/Prog",
                                                           "Rec/Prog", "Rec/Prog", "Rec/Prog", 
                                                           "Rec/Prog", "Rec/Prog", "Rec/Prog",
                                                           "Rec/Prog", "Rec/Prog", "Rec/Prog", 
                                                           "Rec/Prog", "Rec/Prog", "Rec/Prog",
                                                           "Rec/Prog")))
txi2deseq_merged <- DESeqDataSetFromTximport(txi_merged,
                                             colData = sample_id,
                                             design = ~ condition)
deseq_cranio <- DESeq(txi2deseq_merged)
deseq_results <- results(deseq_cranio)


library(fgsea)
library(data.table)
library(apeglm)
set.seed(2)

resultsNames(deseq_cranio)
lfc_cranio_shrink = lfcShrink(deseq_cranio, 
                              coef="condition_Rec.Prog_vs_Initial", type="apeglm")
res_ordered <-  deseq_results[order(deseq_results$pvalue),]
summary(res_ordered)
sum(res_ordered$padj < 0.1, na.rm=TRUE)

plotMA(res_ordered, ylim = c(-2,2), xlim = c(0, 20))
library(org.Hs.eg.db)
symbol2ensbl <- AnnotationDbi::select(org.Hs.eg.db,
                                      key=rownames(lfc_cranio_shrink), 
                                      columns="ENSEMBL",
                                      keytype="SYMBOL")
symbol2ensbl <- as_tibble(symbol2ensbl)
symbol2ensbl

cranio_ordered <- deseq_results[order(deseq_results$stat),]

cranio_xx <- as.data.frame(cranio_ordered)

library(tidyverse)
cranio_necessities <- cranio_xx %>% 
  rownames_to_column("SYMBOL") %>%
  dplyr::select(SYMBOL, stat) %>% 
  as.data.frame() %>%
  na.omit(.) %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(cranio_necessities)

library(msigdbr)
genesets_HALL = msigdbr(species = "Homo sapiens", category = "H")
msigdbr_list_HALL = split(x = genesets_HALL$gene_symbol, f = genesets_HALL$gs_name)
fgsea_HALL = fgsea(pathways= msigdbr_list_HALL, ranks, minSize = 15, 
                   maxSize = 500) %>% 
  as_tibble() %>% 
  arrange(padj)
fgseaResTidy <- fgsea_HALL %>%
  as_tibble() %>%
  arrange(desc(NES))

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

genesets_reactome = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME" )
msigdbr_list_reactome = split(x = genesets_reactome$gene_symbol, f = genesets_reactome$gs_name)
fgsea_reactome = fgsea(pathways= msigdbr_list_reactome, ranks, minSize = 15, 
                       maxSize = 500) %>% 
  as_tibble() %>% 
  arrange(padj)
fgsea_reactome_tidy <- fgsea_reactome %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  rownames_to_column()

genesets_kegg = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG" )
msigdbr_list_kegg = split(x = genesets_kegg$gene_symbol, f = genesets_kegg$gs_name)
fgsea_kegg = fgsea(pathways= msigdbr_list_kegg, ranks, minSize = 15, 
                   maxSize = 500) %>% 
  as_tibble() %>% 
  arrange(padj)
fgsea_kegg_tidy <- fgsea_kegg %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  rownames_to_column()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

dat.gct_normalhypo <- read.delim(file="H:/Desktop/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
                                 skip=2)

