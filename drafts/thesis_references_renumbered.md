# Thesis references — renumbered (Sept 2026)

Changes from the current list (38 entries):

1. **Duplicate removed.** Old [26] and old [37] are the same paper (Petukhov 2022). Both now cite [31].
2. **One new reference.** Haage et al. 2019 replaces Klemm [8] on the microglia/macrophage
   marker-overlap claim (Introduction and Significance). Klemm reports that cell-type-specific
   programs are maintained, so it does not support that claim.
3. **Four existing references now first cited in the Significance section.** Luecken (old 27),
   Aran/SingleR (old 31), Andreatta (old 29) and Zhong (old 30) move forward, because
   numbering follows order of first appearance.
4. **Two new references for the SingleR brain references** (Stage 1, the two "(Ref)"
   placeholders). Both references come from Tabula Muris Senis (FACS; `annotate.R`):
   the atlas paper [33] and the Bioconductor data package used to load it [34].

Assumptions: the current list is in order of first appearance, and Research Objectives and
Significance sits at the end of the Introduction, after Rao and before Eisenbach.

---

## Old → new mapping

| Old | New | Short ID |
|---|---|---|
| 1–10 | 1–10 | unchanged |
| — | **11** | Haage 2019 (new) |
| 11 | 12 | Kierdorf 2019 |
| 12 | 13 | Gonzalez 2022 |
| 13 | 14 | Rao 2021 |
| 27 | 15 | Luecken & Theis 2019 |
| 31 | 16 | Aran 2019 (SingleR) |
| 29 | 17 | Andreatta 2025 |
| 30 | 18 | Zhong 2024 |
| 14 | 19 | Eisenbach 1983 |
| 15 | 20 | Chen 2025 |
| 16 | 21 | Liu 2019 |
| 17 | 22 | Peng 2020 |
| 18 | 23 | Ratzabi 2026 |
| 19 | 24 | Bennett 2016 |
| 20 | 25 | Keren-Shaul 2017 |
| 21 | 26 | Qian 2011 |
| 22 | 27 | Kazanietz 2019 |
| 23 | 28 | Yang & Wang 2015 |
| 24 | 29 | Stringer 2021 (Cellpose) |
| 25 | 30 | NanoString 2023 |
| 26 | 31 | Petukhov 2022 |
| 28 | 32 | Crowell 2025 |
| — | **33** | Tabula Muris Consortium 2020 (new) |
| — | **34** | Soneson 2025, TabulaMurisSenisData (new) |
| 32 | 35 | Kumar 2018 (GSE103548) |
| 33 | 36 | He 2022 (CosMx) |
| 34 | 37 | Danaher 2024 |
| 35 | 38 | Currie 1995 |
| 36 | 39 | Yang 2020 (DecontX) |
| 37 | 31 | Petukhov 2022 (duplicate of old 26) |
| 38 | 40 | Wu 2025 (FastReseg) |

**Replacing numbers in Word:** many numbers swap with each other (e.g. 27→15 while 15→20), so a
plain find-and-replace will overwrite numbers you've already changed. Either first replace every
old number with a temporary marker (e.g. [27] → [#15#]) and strip the markers at the end, or move
the references into a reference manager (Zotero/Mendeley/EndNote), which renumbers automatically.

## In-text changes besides renumbering

- Introduction: "For example, microglia can acquire transcriptional features that resemble those of
  infiltrating macrophages **[8]**" → **[11]**.
- Significance, paragraph 1: [6, 8] · [6, 10] · [11] · [9] · [14] · [6]
- Significance, paragraph 2: SingleR/clustering sentence [15, 16] · copy-number sentence [17] ·
  misclassification sentence [17, 18]
- Tumor identification, Stage 1: SingleR [16] · brain structural reference [33, 34] · brain immune
  reference [33, 34] · GSE103548 [35] · D122 subclone [19].
- Data-quality section: old [37] → [31].
- Tumor identification chapter: check the in-text numbers before mapping. The markdown draft cites
  [30, 31] / [32] SingleR / [33] LLC, one higher than this list ([29, 30] / [31] / [32]).

---

## References

[1] Achrol AS, Rennert RC, Anders C, Soffietti R, Ahluwalia MS, Nayak L, et al. Brain metastases. Nature Reviews Disease Primers. 2019;5(1):5. doi:10.1038/s41572-018-0055-y.

[2] Paisana E, Cascão R, Alvoeiro M, Félix F, Martins G, Guerreiro C, et al. Immunotherapy in lung cancer brain metastases. npj Precision Oncology. 2025;9:130. doi:10.1038/s41698-025-00901-0.

[3] Quail DF, Joyce JA. The microenvironmental landscape of brain tumors. Cancer Cell. 2017;31(3):326-341. doi:10.1016/j.ccell.2017.02.009.

[4] Hanahan D, Weinberg RA. Hallmarks of cancer: the next generation. Cell. 2011;144(5):646-674. doi:10.1016/j.cell.2011.02.013.

[5] Ginhoux F, Jung S. Monocytes and macrophages: developmental pathways and tissue homeostasis. Nature Reviews Immunology. 2014;14(6):392-404. doi:10.1038/nri3671.

[6] Bowman RL, Klemm F, Akkari L, Pyonteck SM, Sevenich L, Quail DF, et al. Macrophage ontogeny underlies differences in tumor-specific education in brain malignancies. Cell Reports. 2016;17(9):2445-2459. doi:10.1016/j.celrep.2016.10.052.

[7] Feng Y, Hu X, Zhang Y, Wang Y. The role of microglia in brain metastases: mechanisms and strategies. Aging and Disease. 2024;15(1):169-185. doi:10.14336/AD.2023.0514.

[8] Klemm F, Maas RR, Bowman RL, Kornete M, Soukup K, Nassiri S, et al. Interrogation of the microenvironmental landscape in brain tumors reveals disease-specific alterations of immune cells. Cell. 2020;181(7):1643-1660. doi:10.1016/j.cell.2020.05.007.

[9] Van Hove H, Martens L, Scheyltjens I, De Vlaminck K, Pombo Antunes AR, De Prijck S, et al. A single-cell atlas of mouse brain macrophages reveals unique transcriptional identities shaped by ontogeny and tissue environment. Nature Neuroscience. 2019;22(6):1021-1035. doi:10.1038/s41593-019-0393-4.

[10] Sun R, Jiang H. Border-associated macrophages in the central nervous system. Journal of Neuroinflammation. 2024;21(1):67. doi:10.1186/s12974-024-03059-x.

[11] Haage V, Semtner M, Vidal RO, Hernandez DP, Pong WW, Chen Z, et al. Comprehensive gene expression meta-analysis identifies signature genes that distinguish microglia from peripheral monocytes/macrophages in health and glioma. Acta Neuropathologica Communications. 2019;7(1):20. doi:10.1186/s40478-019-0665-y.

[12] Kierdorf K, Masuda T, Jordão MJC, Prinz M. Macrophages at CNS interfaces: ontogeny and function in health and disease. Nature Reviews Neuroscience. 2019;20(9):547-562. doi:10.1038/s41583-019-0201-x.

[13] Gonzalez H, Mei W, Robles I, Hagerling C, Allen BM, Nanjaraj A, et al. Cellular architecture of human brain metastases. Cell. 2022;185(4):729-745.e20. doi:10.1016/j.cell.2021.12.043.

[14] Rao A, Barkley D, França GS, Yanai I. Exploring tissue architecture using spatial transcriptomics. Nature. 2021;596(7871):211-220. doi:10.1038/s41586-021-03634-9.

[15] Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. Molecular Systems Biology. 2019;15(6):e8746. doi:10.15252/msb.20188746.

[16] Aran D, Looney AP, Liu L, Wu E, Fong V, Hsu A, et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology. 2019;20(2):163-172. doi:10.1038/s41590-018-0276-y.

[17] Andreatta M, Garnica J, Carmona SJ. Identification of malignant cells in single-cell transcriptomics data. Communications Biology. 2025;8(1):1264. doi:10.1038/s42003-025-08695-4.

[18] Zhong Z, Hou J, Yao Z, Dong L, Liu F, Yue J, et al. Domain generalization enables general cancer cell annotation in single-cell and spatial transcriptomics. Nature Communications. 2024;15(1):1929. doi:10.1038/s41467-024-46413-6.

[19] Eisenbach L, Segal S, Feldman M. MHC imbalance and metastatic spread in Lewis lung carcinoma clones. International Journal of Cancer. 1983;32(1):113-120. doi:10.1002/ijc.2910320118.

[20] Chen Y, Zhang A, Wang J, Pan H, Liu L, Li R. Refining lung cancer brain metastasis models for spatiotemporal dynamic research and personalized therapy. Cancers. 2025;17(9):1588. doi:10.3390/cancers17091588.

[21] Liu Z, Gu Y, Chakarov S, Bleriot C, Kwok I, Chen X, et al. Fate mapping via Ms4a3-expression history traces monocyte-derived cells. Cell. 2019;178(6):1509-1525.e19. doi:10.1016/j.cell.2019.08.009.

[22] Peng L, Wang Y, Fei S, Wei C, Tong F, Wu G, et al. The effect of combining Endostar with radiotherapy on blood vessels, tumor-associated macrophages, and T cells in brain metastases of Lewis lung cancer. Translational Lung Cancer Research. 2020;9(3):745-760. doi:10.21037/tlcr-20-500.

[23] Ratzabi A, Caspit IM, Telechi I, Kim JS, Vaknine H, Blinder P, Jung S, Stein R. Brain metastases exhibit distinct spatial patterns of resident and infiltrating macrophages. Cell Death Discovery. 2026;12(1):211. doi:10.1038/s41420-026-03084-0.

[24] Bennett ML, Bennett FC, Liddelow SA, Ajami B, Zamanian JL, Fernhoff NB, et al. New tools for studying microglia in the mouse and human CNS. Proceedings of the National Academy of Sciences. 2016;113(12):E1738-E1746. doi:10.1073/pnas.1525528113.

[25] Keren-Shaul H, Spinrad A, Weiner A, Matcovitch-Natan O, Dvir-Szternfeld R, Ulland TK, et al. A unique microglia type associated with restricting development of Alzheimer's disease. Cell. 2017;169(7):1276-1290.e17. doi:10.1016/j.cell.2017.05.018.

[26] Qian BZ, Li J, Zhang H, Kitamura T, Zhang J, Campion LR, et al. CCL2 recruits inflammatory monocytes to facilitate breast-tumour metastasis. Nature. 2011;475(7355):222-225. doi:10.1038/nature10138.

[27] Kazanietz MG, Durando M, Cooke M. CXCL13 and its receptor CXCR5 in cancer: inflammation, immune response, and beyond. Frontiers in Endocrinology. 2019;10:471. doi:10.3389/fendo.2019.00471.

[28] Yang Z, Wang KKW. Glial fibrillary acidic protein: from intermediate filament assembly and gliosis to neurobiomarker. Trends in Neurosciences. 2015;38(6):364-374. doi:10.1016/j.tins.2015.04.003.

[29] Stringer C, Wang T, Michaelos M, Pachitariu M. Cellpose: a generalist algorithm for cellular segmentation. Nature Methods. 2021;18(1):100-106. doi:10.1038/s41592-020-01018-x.

[30] NanoString Technologies. Evaluating the technical performance of single-cell spatial biology with CosMx Spatial Molecular Imager. Seattle (WA): NanoString Technologies; 2023. White paper.

[31] Petukhov V, Xu RJ, Soldatov RA, Cadinu P, Khodosevich K, Moffitt JR, et al. Cell segmentation in imaging-based spatial transcriptomics. Nature Biotechnology. 2022;40(3):345-354. doi:10.1038/s41587-021-01044-w.

[32] Crowell HL, Dong Y, Billato I, Cai P, Emons M, Gunz S, et al. Orchestrating spatial transcriptomics analysis with Bioconductor. bioRxiv. 2025. doi:10.1101/2025.11.20.688607.

[33] Tabula Muris Consortium. A single-cell transcriptomic atlas characterizes ageing tissues in the mouse. Nature. 2020;583(7817):590-595. doi:10.1038/s41586-020-2496-1.

[34] Soneson C, Machlab D, Marini F, Astrologo S. TabulaMurisSenisData: bulk and single-cell RNA-seq data from the Tabula Muris Senis project. R package version 1.16.0. Bioconductor; 2025. doi:10.18129/B9.bioc.TabulaMurisSenisData.

[35] Kumar R. Analysis of mRNA expression in Lewis lung carcinoma (LLC1) cells and MLE 12 cells [dataset]. Gene Expression Omnibus; 2018. Accession no. GSE103548. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103548. Associated publication: Ryan ZC, Craig TA, Wang X, Delmotte P, Salisbury JL, Lanza IR, et al. 1α,25-dihydroxyvitamin D3 mitigates cancer cell-mediated mitochondrial dysfunction in human skeletal muscle cells. Biochemical and Biophysical Research Communications. 2018;496(2):746-752. doi:10.1016/j.bbrc.2018.01.092.

[36] He S, Bhatt R, Brown C, Brown EA, Buhr DL, Chantranuvatana K, et al. High-plex imaging of RNA and proteins at subcellular resolution in fixed tissue by spatial molecular imaging. Nature Biotechnology. 2022;40(12):1794-1806. doi:10.1038/s41587-022-01483-z.

[37] Danaher P. How does background impact CosMx data, and when does it matter? CosMx Analysis Scratch Space, Bruker Spatial Biology; 2024 Aug 21. https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/background/ (accessed [date]).

[38] Currie LA. Nomenclature in evaluation of analytical methods including detection and quantification capabilities (IUPAC Recommendations 1995). Pure and Applied Chemistry. 1995;67(10):1699-1723. doi:10.1351/pac199567101699.

[39] Yang S, Corbett SE, Koga Y, Wang Z, Johnson WE, Yajima M, et al. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. Genome Biology. 2020;21(1):57. doi:10.1186/s13059-020-1950-6.

[40] Wu L, Beechem JM, Danaher P. Using transcripts to refine image based cell segmentation with FastReseg. Scientific Reports. 2025;15:30508. doi:10.1038/s41598-025-08733-5.
