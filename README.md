# Projet ReproHackaton

## Ressources
- Dépôt github (contient le code du projet) 
- Google drive (contient les fichiers binaires (articles, documentations, rapport)) 
- Données: issues de Harbour et al. (Nat. Genet. 2013) disponibles au SRA, numéros [SRA062369](https://www.ncbi.nlm.nih.gov/sra?term=SRA062369) et [SRA062359](https://www.ncbi.nlm.nih.gov/sra?term=SRA062359). Attention : *Exome sequences and RNA-seq data are available at the NCBI Sequence Read Archive (SRA) under accessions SRA062369 and SRA062359, respectively).* A voir.

## Objectifs
Dans un article de 2013, [Harbour et al.](https://drive.google.com/file/d/1mR2oxIx7IG2UqzZr1kt1vVCcWvMr6b8B/view?usp=sharing) publient des données de RNAseq de patients atteints d'un mélanome, et ayant ou non le gène SF3B1 muté ; alors que cette protéine est impliquée dans l'épissage, ils ne mettent pas en évidence de différence entre les deux groupes de patients. Un second article publié par [Furney et al.](https://drive.google.com/file/d/1MSxQ1XNcuXBHLKFrOiXP3Xhky4Q00pmb/view?usp=sharing) met en évidence des différences d'épissage entre les deux groupes en ré-analysant le même jeu de données.

Qu'en est-il réellement? Nous allons analyser ces mêmes données et tenter de trouver des différences d'epression chez certains gènes entre les deux groupes de patients. L'enjeu principal est de produire une analyse qui soit complètement reproductible.

## Questions pour la réunion du 16/10/2020

### Concernant les articles
- Différences entre les deux jeux de données ? RNAseq vs. WES ? Deux cohortes de patients ?

### Techniques
- Formuler plus précisément la problématique: gènes différentiellement exprimés ? Gènes différentiellement épissés ? Se focaliser sur quoi ?
- Utiliser dockerhub?
- Importance des résultats obtenus dans l'évaluation ? Si oui, comment estimer si les différences qu'on trouve entre échantillons sont significatives ?
- Rapport en anglais?