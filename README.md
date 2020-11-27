# Projet ReproHackaton

## Organisation du dépôt

### Branche main

Cette branche contient l'ensemble du code nécessaire à l'exécution du pipeline. Le fichier [main.nf](https://github.com/hugovaysset/ReproHackaton/blob/main/main.nf) contient l'ensemble des process nécessaire à l'éxécution du pipeline. Le fichier de configuration [nextflow.config](https://github.com/hugovaysset/ReproHackaton/blob/main/nextflow.config) spécifie la configuration utlisée par le pipeline et en particulier les images docker utilisées par chacun des process. Le dossier [bin](https://github.com/hugovaysset/ReproHackaton/tree/main/bin) contient les scripts R utilisées par les process [statAnalysis](https://github.com/hugovaysset/ReproHackaton/blob/main/main.nf#L205) et [statAnalysisSplicing](https://github.com/hugovaysset/ReproHackaton/blob/main/main.nf#L225) pour l'analyse de l'expresion différentielle et l'épissage atlernatif respectivement. 
Le dossier [Docker](https://github.com/hugovaysset/ReproHackaton/tree/main/Docker) contient les Dockerfiles utilisés pour la création des images docker nécessaire à l'exécution du pipeline nextflow. 

### Banche results_DE_genes 
Cette branche contient les résultats d'une première version du pipeline ne permettant que la recherche de gènes différentiellement exprimés. 
Le dossier [reports](https://github.com/hugovaysset/ReproHackaton/tree/results_DE_genes/reports) contient les rapports d'exécution du pipeline généré par nextflow.
Le dossier [results](https://github.com/hugovaysset/ReproHackaton/tree/results_DE_genes/results/DE_genes) contient la table de gènes différentiellement exprimés détectés par `DESeq2`. 

### Branche results_splicing
Cette branche contient les résultat de la dernière version du pipeline avec le contrôle qualité des reads, la recherche de gènes différentiellement exprimés et la recherche d'épissage alternatif. 
Le dossier [reports](https://github.com/hugovaysset/ReproHackaton/tree/results_splicing/reports) contient les rapports d'exécution du pipeline généré par nextflow.
Le dossier [results](https://github.com/hugovaysset/ReproHackaton/tree/results_splicing/results) contient des résultats de l'exécution du pipeline : 
- Contrôle qualité : [FastQC](https://github.com/hugovaysset/ReproHackaton/tree/results_splicing/results/fastqc_results) et [Fastq Screen](https://github.com/hugovaysset/ReproHackaton/tree/results_splicing/results/fqscreen_results)
- [Épissage alternatif](https://github.com/hugovaysset/ReproHackaton/tree/results_splicing/results/DE_splicing)
- [Expression différentielle](https://github.com/hugovaysset/ReproHackaton/tree/results_splicing/results/DE_genes)
- [Rapport MultiQC](https://github.com/hugovaysset/ReproHackaton/tree/results_splicing/results/multiqc_results)

## Ressources
- Dépôt github (contient le code du projet)
- Google drive (contient les fichiers binaires (articles, documentations, rapport))
- Données: issues de Harbour et al. (Nat. Genet. 2013) disponibles au [SRA](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP017413&o=acc_s%3Aa).

## Objectifs
Dans un article de 2013, [Harbour et al.](https://drive.google.com/file/d/1mR2oxIx7IG2UqzZr1kt1vVCcWvMr6b8B/view?usp=sharing) ([Identifiant SRA](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP017413)) publient des données de RNAseq de patients atteints d'un mélanome, et ayant ou non le gène *SF3B1* muté ; alors que cette protéine est impliquée dans l'épissage, ils ne mettent pas en évidence de différence entre les deux groupes de patients. Un second article publié par [Furney et al.](https://drive.google.com/file/d/1MSxQ1XNcuXBHLKFrOiXP3Xhky4Q00pmb/view?usp=sharing) met en évidence des différences d'épissage entre les deux groupes en ré-analysant le même jeu de données.

Nous nous intéressons dans un premier temps à identifier quels sont les gènes qui sont différentiellement exprimés entre les deux cohortes de patients. Nous portons une attention particulière à la reproductibilité de notre pipeline d'analyse bioinformatique.

## Méthodes

Pour garantir la reproductibilité de notre pipeline d'analyse, nous utilisons le gestionnaire de workflow *nextflow*. Les étapes du workflow sont les suivantes:
1. Récupération des données de séquençage à l'aide de [SRA toolkit](https://github.com/ncbi/sra-tools).
2. Récupération du génome de référence sur le site du [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39) et les [annotations](ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz) de ce génome.
3. Contrôle qualité des reads avec [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) et [Fastq Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
4. Alignement des reads sur le génome de référence en utilisant l'outil [STAR](https://github.com/alexdobin/STAR).
5. Attribution des reads à chacun des gènes à partir du résultat de l'alginement (fichier BAM) via l'outil [featureCounts](http://bioinf.wehi.edu.au/featureCounts/).
6. Analyse statistique avec le package R [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
7. Analyse de l'épissage alternatif avec le package R [DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html)
8. Aggrégation des résultats avec [MultiQC](https://multiqc.info/)

## Notes

### Lancer le pipeline sur la VM
Pour pouvoir lancer le pipeline nextflow an arrière plan et donc pouvoir se déconnecter de sa session ssh: 
- Lancer le processus en arrière plan : `nextflow -bg -q run main.nf`
- L'option `-bg` permet l'exécution du pipeline en arrière plan
- L'option `-q` évite l'envoie les messages de progression du pipeline sur la console
- `nextflow log` permet de connaitre le statut d'éxécution du pipeline.
