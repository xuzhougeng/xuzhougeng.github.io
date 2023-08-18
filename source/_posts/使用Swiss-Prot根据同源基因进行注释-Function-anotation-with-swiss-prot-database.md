---
title: 使用Swiss-Prot根据同源基因进行注释
date: 2019-09-02 10:22:20.408
updated: 2019-09-03 16:12:26.285
url: /archives/Function-anotation-with-swiss-prot-database
categories: 生信软件工具箱
tags: 注释
---

第一步: 在[uniprot](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/)下载UniProt 上植物dat格式的注释文件。

```bash
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_plants.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_plants.dat.gz
```

将两个dat合并到成一个文件

```bash
zcat uniprot_sprot_plants.dat.gz uniprot_trembl_plants.dat.gz > uniprot_plants.dat
```

第二步: 从dat中提取fasta序列

```bash
dat=uniprot_plants.dat
awk '{if (/^ /) {gsub(/ /, ""); print} else if (/^AC/) print ">" $2}' $1 > ${1%%.dat}_AC.fasta
```

第三步: 建立DIAMOND或NCBI BLAST+索引

```bash
diamond makedb --in uniprot_plants_AC.fasta -d uniprot_plants_AC
```

第四步: 使用DIAMOND或NCBI BLAST+进行比对

```bash
diamond blastp -d /data/database/UniProt-Plant/uniprot_plants_AC.dmnd -q proteins.fasta --evalue 1e-5 > blastp.outfmt6
```

第五步: 从DIMAMOND或NCBI BLAST+的比对结果中筛选每个query的最佳subject

```bash
python -m jcvi.formats.blast best -n 1 blastp.outfmt6
```

第六步: 使用add_annotation_from_dat.py(代码在[GitHub](https://github.com/xuzhougeng/myscripts)上)根据blastp输出从dat中提取GO/KEGG/同源基因。运行在Python2/3环境中，需要安装BioPython

```bash
python ~/myscripts/add_annotation_from_dat.py blastp.outfmt6.best /data/database/UniProt-Plant/uniprot_plants.dat
```

之后会输出swiss_annotation.tsv， 输出信息包括如下几列

- gene id
- uniprot accession
- identity
- homology species
- EnsemblPlants
- GO ID
- GO component, CC/MF/BP
- evidence

