# Dataset construction
We curated a set of marker proteins which we expected to be present in all eCIS systems of the families afp11, afp15 and afp1-5 from the dbeCIS database. We then searched for these protein against NCBI's databases env_nt and nt_prok (version 2022-06-14) using tblastn v2.13.0 (min identity 25%, coverage 50%) and extracted all DNA loci within 50 kbp of at least one hit from each of the tree marker protein families. We then searched for homologs to our experimentally validated toxins in these loci using tblasn with the same criteria.
We then extracted the 20 N-terminal amino acids, and homology reduced these with a 90% identity threshold. The search yielded 8 new potential signal peptides.

# Model building
We now had 17 prospective signal peptides in our positive dataset: 6 experimentally validated to pack, 3 additional sequences obtained by mutagenesis, and the 8 sequences found by searching NCBI as described above. We created a negative dataset comprised of 3 experimentally verified negative sequences, 4 mutated signal peptided experimentally verified to not pack, 21 toxins from types III, VI, VII secretion systems which we assume to be too distantly related to the eCIS in question to be able to pack, as well as 48 eCIS antitoxins.

Visual inspection of the positive dataset showed a high number of polar and charged amino acids, as well as conspicious pattern of spaced lysine residues; however, a concrete motif was not clear [TODO: Refer to LOGOs figure from paper?]. To search more systematically for a motif, we used our positive and negative datasets to fit a variety of models to predict if a sequence was positive (i.e. led to packing of the cargo protein) or negative, based on our positive and negative datasets. We used three-fold cross validation to measure the model accuracy.
A simple logistic model achieved 96.4% accuracy, compared to a 80% accuracy from a baseline model (always predicting negative).
The logistic model predicted "MPYSSNSKKNETHSKKNERD" to be the optimal 20 amino acid signal peptide.

We note that a homology reduction threshold of 90% is still quite high, meaning the sequences in our dataset were still homologous, such that the models suffered from data leakage between the training and validation splits. For that reason, we have low confidence in the assesment of our model's accuracy. Indeed, the models presumably just learns to recognize any sequence that looks like SE18NT20 [TODO: REAL NAME?], variants of which comprise much of the dataset, hence why the "optimal" sequence shares 15/20 aa with this sequence. Setting a stricter homology reduction threshold of 50% reduces the size of our positive dataset to just 4 sequences, too low for validating a statistical model.

# Code availability
The code used in the bioinformatic analysis, including virtual environments with all software version used can be found at https://github.com/jakobnissen/ecis_search
