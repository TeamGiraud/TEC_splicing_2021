TEC_splicing_2021

Annotations for Alternative splicing events (ASEs) were generated from the hg38 or mm10 GFT files (igenome_hg38.gtf, igenome_mm10.gtf) using MATS 3.0.8.
We considered the MATS-generated “fromGTF.XXX.txt” files for the ASEs: A3SS, A5SS, RI and SE.
Human ASE annotations are in the folder: “Human hg38-based ASEs”
Mouse ASE annotations are in the folder: “Mouse mm10-based ASEs”
Alternatively, custom files based on the same file structure as the “fromGTF.XXX.txt” files can be considered for downstream analysis.

“ASE_label.R” parses the hg39 or mm10 GFT files and the MATS-generated ASE annotation files to tag each transcript isoform with splice-in or -out ASEs.
A single transcript isoform can feature multiple ASEs
1 is for a spliced-in SE
2 a spliced-in A3SS
3 a spliced-in A5SS
4 a spliced-in RI
9 is for a spliced-out ASE
0 is for an unrelated ASE
 The output file is named “ASE_label_output.csv”. An example is shown in the “Sample” folder.
“ASE_label_V2.R” is used by default. A more recent V3 version fixing some “warning issues” is provided.

“PSI_ASE.R” parses the “ASE_label_output.csv” file and the frequency files (for genes and isoforms) obtained from Tophat2, as exemplified for a mTEChi sample in the “Sample” folder.
“PSI_ASE.R” calculates the percent splicing inclusion (PSI) of each ASE listed in the “ASE_label_output.csv” file.
The output file is names “PSI_ASE_output.csv”. The output PSI file for the exemplified mTEChi sample is provided in the “Sample” folder.

