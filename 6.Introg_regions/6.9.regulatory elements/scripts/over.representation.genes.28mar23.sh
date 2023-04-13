# Running over representation analyses - for genes regulated by sEs in MG introgressed regions & MG adaptive introgressed regions

[http://www.webgestalt.org/results/1679990339/#](http://www.webgestalt.org/results/1679990339/#)

#1) candidate genes: mg adaptive introg sE genes

**Job summary**

- **Enrichment method:** ORA
- **Organism:** hsapiens
- **Enrichment Categories:** geneontology_Biological_Process
- **Interesting list:** mg.adap.se.genes_1679990339.txt. **ID type:** ensembl_gene_id
- The interesting list contains **45** user IDs in which **45** user IDs are unambiguously mapped to **45** unique entrezgene IDs and **0** user IDs can not be mapped to any entrezgene ID.
- The GO Slim summary are based upon the **45** unique entrezgene IDs.
- Among **45** unique entrezgene IDs, **42** IDs are annotated to the selected functional categories and also in the reference list, which are used for the enrichment analysis.
- **Reference list:** uploads/background.se.genes_1679990339.txt **ID type:** ensembl_gene_id
- The reference list can be mapped to **3027** entrezgene IDs and **2780** IDs are annotated to the selected functional categories that are used as the reference for the enrichment analysis.

**Parameters for the enrichment analysis:**

- **Minimum number of IDs in the category:** 5
- **Maximum number of IDs in the category:** 2000
- **FDR Method:** BH
- **Significance Level:** Top 10

Based on the above parameters, **10** categories are identified as enriched categories and all are shown in this report.

#but all have fdr>0.05

![Screenshot 2023-03-28 at 11.20.36.png](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/2655030a-e50a-4812-af3c-56da1838f344/Screenshot_2023-03-28_at_11.20.36.png)

ERROR: No significant results are found under threshold fdr 0.05.

[http://www.webgestalt.org/process.php](http://www.webgestalt.org/process.php)

[http://www.webgestalt.org/results/1679995309/#](http://www.webgestalt.org/results/1679995309/#)

#2) candidate genes: mg introg sE genes

- **Enrichment method:** ORA
- **Organism:** hsapiens
- **Enrichment Categories:** geneontology_Biological_Process
- **Interesting list:** mg.se.genes_1679995309.txt. **ID type:** ensembl_gene_id
- The interesting list contains **237** user IDs in which **237** user IDs are unambiguously mapped to **237** unique entrezgene IDs and **0** user IDs can not be mapped to any entrezgene ID.
- The GO Slim summary are based upon the **237** unique entrezgene IDs.
- Among **237** unique entrezgene IDs, **216** IDs are annotated to the selected functional categories and also in the reference list, which are used for the enrichment analysis.
- **Reference list:** uploads/background.se.genes_1679995309.txt **ID type:** ensembl_gene_id
- The reference list can be mapped to **3027** entrezgene IDs and **2780** IDs are annotated to the selected functional categories that are used as the reference for the enrichment analysis.

**Parameters for the enrichment analysis:**

- **Minimum number of IDs in the category:** 5
- **Maximum number of IDs in the category:** 2000
- **FDR Method:** BH
- **Significance Level:** Top 10

Based on the above parameters, **10** categories are identified as enriched categories and all are shown in this report.

![Screenshot 2023-03-28 at 11.22.54.png](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/cac535df-70ce-44d1-91db-df252a61f445/Screenshot_2023-03-28_at_11.22.54.png)

ERROR: No significant results are found under threshold fdr 0.05.

#————————————————-

#me to MK

I performed over-representation analyses in gene ontology terms for two sets of candidate genes, taking the genes which are 1:1 orthologs across primates (human/chimp/gorilla/orang/macaque). using the webgestalt tool, also used for such an ORA analysis in the LCLs paper

Candidate sets -

1) genes regulated by sEs in MG introgressed regions (235 genes)

2) genes regulated by sEs in MG adaptive introgressed regions (45 genes)

Background set - all genes regulated by sEs in gorillas (3010 genes)

Results

1) [http://www.webgestalt.org/results/1679995309/#](http://www.webgestalt.org/results/1679995309/#)

2) [http://www.webgestalt.org/results/1679990339/#](http://www.webgestalt.org/results/1679990339/#)

No category reaches significance at FDR <0.05 with Benjamin-Hochberg correction.

The top categories for each list are related to the LCL cell type (immune functions, brain functions which LCLs also seem to pick up).