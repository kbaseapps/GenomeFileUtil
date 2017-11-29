The files listed here are the truncated contents of a Poplar genome from Phytozome.
The naming convention of the files is typical of a Phytozome directory, though it is common for there to be other files in the same directory (at the bottom is a couple of example `find`s)

The assembly will be found in the assembly directory:
Test_v1.0.fa.gz

The GFF file will be found in the annotation folder:
Test_v1.0.gene.gff3.gz

Note the different use of identifiers and names, in Phytozome, the true original identifier is the Name (for gene and mRNA) and we use mRNA.Name+".CDS" for the CDS identifier
In other words, we ignore the use of the gene model version number in the ID field
$ gunzip -c Test_v1.0.gene.gff3.gz | head
Chr01   phytozomev10    gene    1660    2502    .       -       .       ID=Potri.001G000100.v3.0;Name=Potri.001G000100
Chr01   phytozomev10    mRNA    1660    2502    .       -       .       ID=Potri.001G000100.1.v3.0;Name=Potri.001G000100.1;pacid=27043735;longest=1;Parent=Potri.001G000100.v3.0
Chr01   phytozomev10    CDS     1660    2502    .       -       0       ID=Potri.001G000100.1.v3.0.CDS.1;Parent=Potri.001G000100.1.v3.0;pacid=27043735
Chr01   phytozomev10    gene    2906    6646    .       -       .       ID=Potri.001G000200.v3.0;Name=Potri.001G000200
Chr01   phytozomev10    mRNA    2906    6646    .       -       .       ID=Potri.001G000200.1.v3.0;Name=Potri.001G000200.1;pacid=27045395;longest=1;Parent=Potri.001G000200.v3.0
Chr01   phytozomev10    CDS     6501    6644    .       -       0       ID=Potri.001G000200.1.v3.0.CDS.1;Parent=Potri.001G000200.1.v3.0;pacid=27045395
Chr01   phytozomev10    five_prime_UTR  6645    6646    .       -       .       ID=Potri.001G000200.1.v3.0.five_prime_UTR.1;Parent=Potri.001G000200.1.v3.0;pacid=27045395
Chr01   phytozomev10    CDS     3506    3928    .       -       0       ID=Potri.001G000200.1.v3.0.CDS.2;Parent=Potri.001G000200.1.v3.0;pacid=27045395
Chr01   phytozomev10    CDS     2906    3475    .       -       0       ID=Potri.001G000200.1.v3.0.CDS.3;Parent=Potri.001G000200.1.v3.0;pacid=27045395
Chr01   phytozomev10    gene    8391    8860    .       +       .       ID=Potri.001G000300.v3.0;Name=Potri.001G000300

The protein file is redundant since we translate the protein from the sequence extracted using the GFF coordinates:
Test_v1.0.protein.fa.gz

We use it for two purposes, to test the translations, and, in the case where the genome was assembled from an external lab, to replace the translations.

There are several other files that contain functional annotation:

The 'synonym' file contains aliases for the features within:
Test_v1.0.synonym.txt.gz

The identifiers used are those of the mRNA biotype, but can/should also apply to both the gene and CDS biotype:
$ gunzip -c Test_v1.0.synonym.txt.gz | head
Potri.001G000200.1      POPTR_0001s03740
Potri.001G000300.1      POPTR_0001s03750
Potri.001G000500.1      POPTR_0001s03760        fgenesh4_pg.C_LG_I000006
Potri.001G000600.1      POPTR_0001s03770        eugene3.00010005
Potri.001G000700.1      POPTR_0001s03780        Pt-ALPHA.5      gw1.I.6338.1    POPTR_0001s03790
Potri.001G000700.2      POPTR_0001s03780        Pt-ALPHA.5      gw1.I.6338.1    POPTR_0001s03790
Potri.001G000800.1      POPTR_0001s03800
Potri.001G000900.1      POPTR_0001s03810        fgenesh4_pg.C_LG_I000010
Potri.001G000900.2      POPTR_0001s03810        fgenesh4_pg.C_LG_I000010
Potri.001G000900.3      POPTR_0001s03810        fgenesh4_pg.C_LG_I000010

The 'defline' are the functional annotations that were assigned to each feature by the external lab that assembled the genome
Test_v1.0.defline.txt.gz

The identifiers used are those of the mRNA biotype, but can/should also apply to both the gene and CDS biotype:
$ gunzip -c Test_v1.0.defline.txt.gz | head
Potri.001G000500.1      multi-copper oxidase type 1 family protein; similar to pollen-specific BP10 protein (SP|Q00624|Brassica)(napus ; [ co-ortholog (3of4) of At1g55570, At3g13390, ]
Potri.001G000600.1      multi-copper oxidase type 1 family protein; similar to pollen-specific BP10 protein (SP|Q00624|Brassica)(napus ; [ co-ortholog (2of4) of At1g55570, At3g13390, ]
Potri.001G000700.1      similar to PUR alpha-1 protein; similar to PUR alpha-1 GI:5081612 from (Arabidopsis thaliana); similar to PUR alpha-1 protein; similar to PUR alpha-1 GI:5081612 from (Arabidopsis thaliana); [ co-ortholog (1of2) of At2g32080, ]
Potri.001G000700.2      similar to PUR alpha-1 protein; similar to PUR alpha-1 GI:5081612 from (Arabidopsis thaliana); similar to PUR alpha-1 protein; similar to PUR alpha-1 GI:5081612 from (Arabidopsis thaliana); [ co-ortholog (1of2) of At2g32080, ]
Potri.001G000900.1      similar to expressed protein in Arabidopsis thaliana; [ ortholog of At3g13410,]
Potri.001G000900.2      similar to expressed protein in Arabidopsis thaliana; [ ortholog of At3g13410,]
Potri.001G000900.3      similar to expressed protein in Arabidopsis thaliana; [ ortholog of At3g13410,]
Potri.001G000900.4      similar to expressed protein in Arabidopsis thaliana; [ ortholog of At3g13410,]
Potri.001G000900.5      similar to expressed protein in Arabidopsis thaliana; [ ortholog of At3g13410,]
Potri.001G000900.6      similar to expressed protein in Arabidopsis thaliana; [ ortholog of At3g13410,]

Finally, the 'annotation_info' file contains all the functional annotations that were performed by Phytozome itself:
Test_v1.0.annotation_info.txt.gz

The first few columns contain the feature ids, note that there are typically multiple rows for the locusName(s). The headers are mostly self-explanatory, the 'arabi-' keyword corresponds to orthologous projections from Arabidopsis. Some Phytozome files have 'rice-' keywords in their annotation_info file which correspond to orthologous projects from Rice.
$ gunzip -c Test_v1.0.annotation_info.txt.gz | head
#pacId  locusName       transcriptName  peptideName     Pfam    Panther KOG     KEGG/ec KO      GO      Best-hit-arabi-name     arabi-symbol    arabi-defline
27043735        Potri.001G000100        Potri.001G000100.1      Potri.001G000100.1      PF14111 PTHR31286,PTHR31286:SF2                                 AT2G01050.1             zinc ion binding;nucleic acid binding
27045395        Potri.001G000200        Potri.001G000200.1      Potri.001G000200.1      PF14111 PTHR31286,PTHR31286:SF2         3.5.1.28,3.2.1.96                       AT2G01050.1             zinc ion binding;nucleic acid binding
27045442        Potri.001G000300        Potri.001G000300.1      Potri.001G000300.1
27046906        Potri.001G000400        Potri.001G000400.1      Potri.001G000400.1      PF11779 PTHR33727                                       AT2G30942.1             Protein of unknown function (DUF3317)
27046909        Potri.001G000400        Potri.001G000400.2      Potri.001G000400.2      PF11779 PTHR33727                                       AT2G30942.1             Protein of unknown function (DUF3317)
27046907        Potri.001G000400        Potri.001G000400.3      Potri.001G000400.3      PF11779 PTHR33727                                       AT2G30942.1             Protein of unknown function (DUF3317)
27046908        Potri.001G000400        Potri.001G000400.4      Potri.001G000400.4      PF11779 PTHR33727                                       AT2G30942.1             Protein of unknown function (DUF3317)
27042524        Potri.001G000500        Potri.001G000500.1      Potri.001G000500.1      PF07731,PF07732,PF00394 PTHR11709,PTHR11709:SF27                1.10.3.3        K00423  GO:GO:0055114,GO:GO:0016491,GO:GO:0005507 AT1G55570.1     sks12   SKU5  similar 12
27043996        Potri.001G000600        Potri.001G000600.1      Potri.001G000600.1      PF07731,PF07732,PF00394 PTHR11709,PTHR11709:SF27                1.10.3.3        K00423  GO:GO:0055114,GO:GO:0016491,GO:GO:0005507 AT1G55570.1     sks12   SKU5  similar 12

$ find Ptrichocarpa/v2.2/
Ptrichocarpa/v2.2/
Ptrichocarpa/v2.2/related_data
Ptrichocarpa/v2.2/related_data/Ptrichocarpa_156_repeats.gff3.gz
Ptrichocarpa/v2.2/assembly
Ptrichocarpa/v2.2/assembly/Ptrichocarpa_156_softmasked.fa.gz
Ptrichocarpa/v2.2/assembly/Ptrichocarpa_156_hardmasked.fa.gz
Ptrichocarpa/v2.2/assembly/Ptrichocarpa_156_RM_16mer.fa.gz
Ptrichocarpa/v2.2/assembly/Ptrichocarpa_156_RM.fa.gz
Ptrichocarpa/v2.2/assembly/Ptrichocarpa_156.fa.gz
Ptrichocarpa/v2.2/annotation
Ptrichocarpa/v2.2/annotation/Ptrichocarpa_156_transcript.fa.gz
Ptrichocarpa/v2.2/annotation/Ptrichocarpa_156_synonym.txt.gz
Ptrichocarpa/v2.2/annotation/Ptrichocarpa_156_peptide.fa.gz
Ptrichocarpa/v2.2/annotation/Ptrichocarpa_156_gene_exons.gff3.gz
Ptrichocarpa/v2.2/annotation/Ptrichocarpa_156_gene.gff3.gz
Ptrichocarpa/v2.2/annotation/Ptrichocarpa_156_defline.txt.gz
Ptrichocarpa/v2.2/annotation/Ptrichocarpa_156_cds.fa.gz
Ptrichocarpa/v2.2/annotation/Ptrichocarpa_156_annotation_info.txt.gz
Ptrichocarpa/v2.2/annotation/Ptrichocarpa_156_edited_gene.gff3.gz
Ptrichocarpa/v2.2/Ptrichocarpa_156_readme.txt

$ find Creinhardtii/v4.0/
Creinhardtii/v4.0/
Creinhardtii/v4.0/assembly
Creinhardtii/v4.0/assembly/Chlre4_masked_genomic_scaffolds.fasta.gz
Creinhardtii/v4.0/assembly/Chlre4_genomic_scaffolds.fasta.gz
Creinhardtii/v4.0/annotation
Creinhardtii/v4.0/annotation/FrozenGeneCatalog_20080828_transcripts.fasta.gz
Creinhardtii/v4.0/annotation/FrozenGeneCatalog_20080828_proteins.fasta.gz
Creinhardtii/v4.0/annotation/FrozenGeneCatalog_20080828_genes.gff.gz
Creinhardtii/v4.0/annotation/FrozenGeneCatalog_20080828_IPR.tab.gz
Creinhardtii/v4.0/annotation/Chlre4_best_transcripts_2.fasta.gz
Creinhardtii/v4.0/annotation/Chlre4_best_proteins.fasta.gz
Creinhardtii/v4.0/annotation/Chlre4_best_genes.gff.gz
Creinhardtii/v4.0/annotation/Chlre4_all_transcripts_2.fasta.gz
Creinhardtii/v4.0/annotation/Chlre4_all_proteins.fasta.gz
Creinhardtii/v4.0/annotation/Chlre4_RepeatMasker.gff
Creinhardtii/v4.0/annotation/Chlre4_KOG.tab.gz
Creinhardtii/v4.0/annotation/Chlre4_KEGG.tab.gz
Creinhardtii/v4.0/annotation/Chlre4_IPR.tab.gz
Creinhardtii/v4.0/annotation/Chlre4_GO.tab.gz
Creinhardtii/v4.0/annotation/Chlre4_Augustus5_transcripts.fasta.gz
Creinhardtii/v4.0/annotation/Chlre4_Augustus5_proteins.fasta.gz
Creinhardtii/v4.0/annotation/Chlre4_Augustus5_genes.gff.gz
Creinhardtii/v4.0/ESTs
Creinhardtii/v4.0/ESTs/Chlre4_ESTs.fasta.gz
