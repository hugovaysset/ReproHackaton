# Projet ReproHackaton

## Ressources
- Dépôt github (contient le code du projet)
- Google drive (contient les fichiers binaires (articles, documentations, rapport))
- Données: issues de Harbour et al. (Nat. Genet. 2013) disponibles au [SRA](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP017413&o=acc_s%3Aa).

## Objectifs
Dans un article de 2013, [Harbour et al.](https://drive.google.com/file/d/1mR2oxIx7IG2UqzZr1kt1vVCcWvMr6b8B/view?usp=sharing) ([Identifiant SRA](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP017413)) publient des données de RNAseq de patients atteints d'un mélanome, et ayant ou non le gène *SF3B1* muté ; alors que cette protéine est impliquée dans l'épissage, ils ne mettent pas en évidence de différence entre les deux groupes de patients. Un second article publié par [Furney et al.](https://drive.google.com/file/d/1MSxQ1XNcuXBHLKFrOiXP3Xhky4Q00pmb/view?usp=sharing) met en évidence des différences d'épissage entre les deux groupes en ré-analysant le même jeu de données.

Nous nous intéressons dans un premier temps à identifier quels sont les gènes qui sont différentiellement exprimés entre les deux cohortes de patients. Nous portons une attention particulière à la reproductibilité de notre pipeline d'analyse bioinformatique.

## Méthodes

Pour garantir la reproductibilité de notre pipeline d'analyse, nous utilisons le gestionnaire de workflow *nextflow*. Les étapes du workflow sont les suivantes:
1. Récupération des données de séquençage.
2. Récupération du génome de référence sur le site du [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39) et les [annotations](ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz) de ce génome.
3. Alignement des reads sur le génome de référence en utilisant l'outil [STAR](https://github.com/alexdobin/STAR).
4. Attribution des reads à chacun des gènes à partir du résultat de l'alginement (fichier BAM) via l'outil.
5. Analyse statistique avec le package R [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

## Notes

### Lancer le pipeline sur la VM
Pour pouvoir lancer le pipeline nextflow an arrière plan et donc pouvoir se déconnecter de sa session ssh: 
- Lancer le processus en arrière plan : `nextflow -bg -q run main.nf`
- L'option `-bg` permet l'exécution du pipeline en arrière plan
- L'option `-q` évite l'envoie les messages de progression du pipeline sur la console
- `nextflow log` permet de connaitre le statut d'éxécution du pipeline.
