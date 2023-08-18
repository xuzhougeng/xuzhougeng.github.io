---
title: 「文献」多倍体植物基因组测序组装当前策略
date: 2019-10-09 20:46:39.461
updated: 2019-10-09 20:46:39.461
url: /archives/Current-Strategies-of-Polyploid-Plant-Genome-Sequence-Assembly
categories: 文献阅读
tags: 
---

# 「文献」多倍体植物基因组测序组装当前策略

文献地址: [Current Strategies of Polyploid Plant Genome Sequence Assembly](https://www.frontiersin.org/articles/10.3389/fpls.2018.01660/full)

基因组多倍化主要发生在被子植物中。很多多倍体植物都在农业生产上有重大的价值，例如小麦(*Triticum aestivum*)，花生(*Arachis hypogaea*)，十字花科，马铃薯(*Solanum tuberosum*)，燕麦(*Avena sativa*)，香蕉(*Musa sp.*)，草莓(*Fragaria ananassa*)，咖啡( *Coffea arabica*)等。

多倍体分为两种类型，来自于全基因组加倍的同源多倍体(Autopolyploidy)和物种间/物种内杂交后染色体加倍的异源多倍体(allopolyploidy). 同源多倍体通常会有育性上的问题，而异源多倍体则可能出现杂交优势(heterosis or hybrid vigor).多倍体在表型和基因型上的关系更加复杂，例如它们需要比较复杂的调控才能保证同源基因相互间的表达一致。

在基因组组装上，同源多倍体相对异源多倍体更加困难。这是因为全基因组加倍事件之后通常还会跟着基因组重拍(genome rearrangement), 非典型重组(atypical recombination), 可移动因子启动(transposable element activation)，减数分裂/有丝分裂缺陷(meiotic/mitotic defects)，以及内含子扩张(intron expansions)与DNA缺失。因此组装基因组一大挑战就是不能错误组装了两个亚基因组中的相似片段。 

作者在[NCBI](https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/)查询并总结了到2018年为止已发表的多倍体物种，我更新了草莓(Fragaria × ananassa)，香蕉(*Musa balbisiana*)和甘蔗(*Saccharum spontaneum* L.)


| ID   | Organism name                   | Genome size (Mb) | Current status | 1st Release date in NCBI | Ploidy level                                 | References/center                                            |
| ---- | ------------------------------- | ---------------- | -------------- | ------------------------ | -------------------------------------------- | ------------------------------------------------------------ |
| 1    | Arabidopsis lyrata subsp lyrata | 206.823          | Scaffold       | 2009/11/30               | Tetraploid                                   | Hu et al., 2011                                              |
| 2    | Glycine max                     | 978.972          | Chromosome     | 2010/1/5                 | Allotetraploid                               | Schmutz et al., 2010                                         |
| 3    | Triticum aestivum               | 15344.7          | Chromosome 3B  | 2010/7/15                | Allohexaploid                                | Choulet et al., 2010                                         |
| 4    | Solanum tuberosum               | 705.934          | Scaffold       | 2011/5/24                | Autotetraploid                               | Potato Genome Sequencing Consortium, 2011                    |
| 5    | Actinidia chinensis             | 604.217          | Contig         | 2013/9/16                | Tetraploid                                   | Huang et al., 2013                                           |
| 6    | Fragaria orientalis             | 214.356          | Scaffold       | 2013/11/27               | Tetraploid                                   | Hirakawa et al., 2014                                        |
| 7    | Fragaria x ananassa             | 805.488          | Chromosome     | 2019/2/25                | Allooctaploid                                | Edger et al., 2019                                           |
| 8    | Beta vulgaris                   | 566.55           | Chromosome     | 2013/12/18               | 2n, 4n (Beyaz et al., 2013)                  | Dohm et al., 2014                                            |
| 9    | Oryza minuta                    | 45.1659          | Chromosome     | 2014/4/16                | Tetraploid                                   | Oryza Chr3 Short Arm Comparative Sequencing Project          |
| 10   | Camelina sativa                 | 641.356          | Chromosome     | 2014/4/17                | Hexaploid                                    | Kagale et al., 2014                                          |
| 11   | Brassica napus                  | 976.191          | Chromosome     | 2014/5/5                 | Allotetraploid                               | Chalhoub et al., 2014                                        |
| 12   | Brassica oleracea var. oleracea | 488.954          | Chromosome     | 2014/5/22                | Hexaploid                                    | NCBI                                                         |
| 13   | Nicotiana tabacum               | 3643.47          | Scaffold       | 2014/5/29                | Allotetraploid                               | Sierro et al., 2014                                          |
| 14   | Eragrostis tef                  | 607.318          | Scaffold       | 2015/4/8                 | Allotetraploid                               | Cannarozzi et al., 2014                                      |
| 15   | Gossypium hirsutum              | 2189.14          | Chromosome     | 2015/4/29                | Allotetraploid                               | Li et al., 2015                                              |
| 16   | Zoysia japonica                 | 334.384          | Scaffold       | 2016/3/15                | Tetraploid                                   | Tanaka et al., 2016                                          |
| 17   | Zoysia matrella                 | 563.439          | Scaffold       | 2016/3/15                | Allotetraploid                               | Tanaka et al., 2016                                          |
| 18   | Zoysia pacifica                 | 397.01           | Scaffold       | 2016/3/15                | Allotetraploid                               | Tanaka et al., 2016                                          |
| 19   | Musa itinerans                  | 455.349          | Scaffold       | 2016/5/21                | 2n, 3n hybrids (Wu et al.,2016)              | South China Botanic Garden, CAS                              |
| 20   | Rosa x damascena                | 711.72           | Scaffold       | 2016/6/13                | Tetraploid                                   | BIO-FD & C CO., LTD                                          |
| 21   | Chenopodium quinoa              | 1333.55          | Scaffold       | 2016/7/11                | Tetraploid                                   | Jarvis et al., 2017                                          |
| 22   | Brassica juncea var. tumida     | 954.861          | Chromosome     | 2016/7/19                | Allotetraploid                               | Zhejiang University                                          |
| 23   | Hibiscus syriacus               | 1748.25          | Scaffold       | 2016/7/29                | 2n, 3n, 4n (Van   Huylenbroeck et al., 2000) | Korea Research Institute of Science and Biotechnology (Kim et   al., 2017) |
| 24   | Gossypium barbadense            | 2566.74          | Scaffold       | 2016/10/28               | Tetraploid                                   | Huazhong Agricultural University                             |
| 25   | Momordica charantia             | 285.614          | Scaffold       | 2016/12/27               | 2n to 6n (Kausar et al., 2015)               | Urasaki et al., 2016                                         |
| 26   | Drosera capensis                | 263.788          | Scaffold       | 2016/12/30               | Tetraploid (Rothfels and Heimburger, 1968)   | Butts et al., 2016                                           |
| 27   | Capsella bursa-pastoris         | 268.431          | Scaffold       | 2017/1/29                | Tetraploid                                   | Lomonosov Moscow State University                            |
| 28   | Saccharum hybrid cultivar       | 1169.95          | Contig         | 2017/3/3                 | It varies (D’Hont, 2005)                     | Riaño-Pachón and Mattiello, 2017                             |
| 29   | Xerophyta viscosa               | 295.462          | Scaffold       | 2017/3/31                | Hexaploid                                    | Costa et al., 2017                                           |
| 30   | Triticum dicoccoides            | 10495            | Chromosome     | 2017/5/18                | Tetraploid                                   | WEWseq consortium                                            |
| 31   | Utricularia gibba               | 100.689          | Chromosome     | 2017/5/31                | 16-ploid                                     | Lan et al., 2017                                             |
| 32   | Eleusine coracana               | 1195.99          | Scaffold       | 2017/6/8                 | Allotetraploid                               | Hittalmani et al., 2017                                      |
| 33   | Dioscorea rotundata             | 456.675          | Chromosome     | 2017/7/28                | Tetraploid                                   | Iwate Biotechnology Research Center                          |
| 34   | Ipomoea batatas                 | 837.013          | Contig         | 2017/8/26                | Autohexaploid                                | Yang et al., 2017                                            |
| 35   | Echinochloa crus-galli          | 1486.61          | Scaffold       | 2017/10/23               | Hexaploid                                    | Zhejiang University                                          |
| 36   | Pachycereus pringlei            | 629.656          | Scaffold       | 2017/10/31               | Autotetraploid                               | Zhou et al., 2017                                            |
| 37   | Olea europaea                   | 1141.15          | Chromosome     | 2017/11/1                | 2n, 4n, 6n (Besnard et al.,      2007)       | Unver et al., 2017                                           |
| 38   | Monotropa hypopitys             | 2197.49          | Contig         | 2018/1/3                 | Hexaploid                                    | Institute of Bioengineering, RAS                             |
| 39   | Dactylis glomerata              | 839.915          | Scaffold       | 2018/1/19                | Autotetraploid                               | Sichuan Agricultural University                              |
| 40   | Panicum miliaceum               | 848.309          | Scaffold       | 2018/1/23                | Allotetraploid                               | China Agricultural University                                |
| 41   | Euphorbia esula                 | 1124.89          | Scaffold       | 2018/2/6                 | Hexaploid                                    | USDA-ARS                                                     |
| 42   | Santalum album                  | 220.961          | Scaffold       | 2018/2/12                | 2n, 4n etc (Xin-Hua et al.,      2010)       | Center for Cellular and Molecular Platforms                  |
| 43   | Avena sativa                    | 67.3266          | Contig         | 2018/2/26                | Hexaploid                                    | The Sainsbury Laboratory                                     |
| 44   | Panicum miliaceum               | 850.677          | Chromosome     | 2018/4/9                 | Tetraploid                                   | Shanghai Center for Plant Stress Biology                     |
| 45   | Arachis monticola               | 2618.65          | Chromosome     | 2018/4/23                | Tetraploid                                   | Henan Agricultural University                                |
| 46   | Arachis hypogaea                | 2538.28          | Chromosome     | 2018/5/2                 | Allotetraploid                               | International Peanut Genome Initiative                       |
| 47   | Artemisia annua                 | 1792.86          | Scaffold       | 2018/5/8                 | Tetraploid                                   | Shen et al., 2018                                            |
| 48   | Saccharum spontaneum L.         | 2.9 G            | Chromosome     | 2018/09/10               | octoploid                                    | Zhang et al., 2018                                           |
| 49   | Musa balbisiana                 | 430              | Chromosome     | 2019/7/15                | Tetraploid                                   | Wang et al., 2019                                            |

在倍性预测上，有两种方法可以使用

- 使用流式细胞仪测量C-值, 参考网站[Plant DNA C-values Database](https://cvalues.science.kew.org/)
- 基于NGS的生信策略: ploidyNG ConPADE 

而在单倍型组装上，作者列了如下工具，当然最靠谱的肯定是最新的，也就是HapCUT2

- HapCompass (Aguiar and Istrail, 2012) 
- HaploSim (Bastiaansen et al., 2012) 
- HapCut (Bansal and Bafna, 2008) 
- HapCUT2 (Edge et al., 2017) 

在解决多倍体问题上，作者给出了两种策略

- 基因组上: 尽量挑选单倍型，或者先测二倍体祖先
- 分析流程上:  三代测序, BioNano, HiC

最终，作者总结了目前植物可用的资源网站


| DB name            | Resources                                        | Plants                                                       | URL                                                          |
| ------------------ | ------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Genbank            | Genomic                                          | Various plant species                                        | https://www.ncbi.nlm.nih.gov/genbank/                        |
| EMBL               | Genomic                                          | Various plant species                                        | https://www.ebi.ac.uk/                                       |
| DDBJ               | Genomic                                          | Various plant species                                        | http://www.ddbj.nig.ac.jp/                                   |
| UniProt            | Protein and functional                           | Various plant species                                        | http://www.uniprot.org/                                      |
| NCBI               | Genomic                                          | Various plant species                                        | https://www.ncbi.nlm.nih.gov/                                |
| GOLD               | Genomic, metagenomics, transcriptomic            | Various plant species                                        | https://gold.jgi.doe.gov/cgi-bin/GOLD/bin/gold.cgi           |
| Phytozom           | Genomic                                          | 92 assembled and annotated plant species                     | https://phytozome.jgi.doe.gov/pz/portal.html                 |
| Plantgdb           | Genomic, transcriptomic                          | 27 assembled and annotated plant species                     | http://www.plantgdb.org/                                     |
| Sol                | Genomic                                          | 11 Solanaceae species                                        | https://solgenomics.net/                                     |
| Gramene            | Genomic, genetic markers, QTLs                   | 53 plant species                                             | http://www.gramene.org/                                      |
| MaizeGCB           | Genomic, annotations, tool host                  | Zea mays                                                     | https://www.maizegdb.org/                                    |
| Tair               | Genetic and molecular biology data               | Arabidopsis thaliana                                         | https://www.arabidopsis.org/                                 |
| CottonGE           | Genomic, Genetic and breeding resources          | 49 Gossypium species                                         | https://www.arabidopsis.org/                                 |
| PLEXdb             | Gene expression                                  | 14 plant species                                             | http://www.plexdb.org/                                       |
| RicePro            | Gene expression                                  | Oryza sativa                                                 | http://ricexpro.dna.affrc.go.jp/                             |
| CerealsDB          | Genetic markers                                  | Triticum aestivum                                            | http://www.cerealsdb.uk.net/cerealgenomics/CerealsDB/indexNEW.php |
| PeanutBa           | Genome, MAS, QTLs, Germplasm                     | Arachis hypogaea                                             | https://peanutbase.org/                                      |
| SoyKb              | Genetic markers, genomic resources               | Glycine max                                                  | http://soykb.org/                                            |
| SoyBase            | Genetic markers, QTLs, genomic resources         | G. max                                                       | https://soybase.org/                                         |
| PGDBj              | Genetic markers, QTLs, genomic resources         | 80 plant species                                             | http://pgdbj.jp/                                             |
| SNP-Seek           | Genotype, Phenotype and Variety information      | O. sativa                                                    | http://snp-seek.irri.org/                                    |
| GrainGene          | Genome, Genetic markers, QTLs, genomic resources | T. aestivum, Hordeum vulgare, Secale cereale, Avena sativa etc | https://wheat.pw.usda.gov/GG3/                               |
| ASRP               | small RNA                                        | A. thaliana                                                  | http://asrp.danforthcenter.org/                              |
| CSRDB              | small RNA                                        | Z. mays                                                      | http://sundarlab.ucdavis.edu/smrnas/                         |
| BrassicaIn         | Genomic                                          | 7 Brassica species                                           | http://brassica.info/                                        |
| BRAD               | Genomics, Genetic Markers and Maps               | Brassica                                                     | http://brassicadb.org/brad/                                  |
| Ensembl Plants     | Genomic                                          | 45 plant species                                             | http://plants.ensembl.org/index.html                         |
| Ipomoea Genome Hub | Genomic, EST                                     | Ipomoea batatas                                              | https://ipomoea-genome.org/                                  |
| PGSC               | Genomic, annotation                              | S. tuberosum, S.chacoense                                    | http://solanaceae.plantbiology.msu.edu/pgsc_download.shtml   |
| GDR                | Genomics, Genetics, breeding                     | Rosaceae                                                     | https://www.rosaceae.org/analysis/266                        |
| HWG                | Genomics, Transcriptomics, Genetic Markers       | Forest trees and woody plants                                | https://www.hardwoodgenomics.org/                            |
