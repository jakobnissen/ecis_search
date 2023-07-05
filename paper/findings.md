# Main findings

### Data search
* HMM searching the SKESKESKE is not useful: Too many hits and no way of determining if it's meaningful or not
* Searching in NCBI's prokaryote nucleotide database yielded 8 new potential toxins.
    - Two of these have low (< 70%) amino acid identity to any of the already tested toxins
    - Searching in a large collection of environmental contigs (NCBI's nt_env) yielded no further toxins.

The final dataset constists of:
    6 experimentally verified sequences:
        - SE18NT20
        - YR17NT20
        - CyaA
        - EPTox20
        - YPTox20
        - SETox20
    3 experimentally verified mutations of Afp18N20
    8 positives found by homology search in NCBI

    3 experimentally verified negatives
        - ANV21619.1 AfpX17
        - SEAfp17-AAT48354.1-KHA73_24135
        - PAU_RS10120-RHS

    4 experimentally verified mutational negatives
    21 toxins from types II, VI, VII assumed to be negative

    The positives were homology reduced at 90% identity ONLY on the first 20aa

### Homology reduction:
Positives at 70% reduction:
SE18NT20                      MPYSSESKEKETHSKETERD
YR17NT20                      MPYFNKSKKNEIRPEKSKEE
CyaA                          MPRYSNSQRTPTQSTKNTRR
EPTox20                       MPYFNELNEKETRSKETESG
CP025698.1_1685283:1688714    MPYSRESKEKDTHAKGSKQD
FM162591.1_3913240:3914247    MPRYANYQINPKQNIKNSHG

Positives at 90% reduction:
SE18NT20                      MPYSSESKEKETHSKETERD
YR17NT20                      MPYFNKSKKNEIRPEKSKEE
CyaA                          MPRYSNSQRTPTQSTKNTRR
EPTox20                       MPYFNELNEKETRSKETESG
YPTox20                       MLYSSESKEKKTHSKETERD
SETox20                       MPYSSESKLKDTHLKEAESD
Afp18N20KtAaa                 MPYSSESAEAETHSAETERD
Afp18N20KTtAaa                MPYSSESAEAEAHSAEAERD
Afp18N20EtAaa                 MPYSSASKAKATHSKATARD
CP020335.1_1949579:1950583    MPRYSNSQRIPTQNSKNSRR
CP025698.1_1685283:1688714    MPYSRESKEKDTHAKGSKQD
MT039196.1_91004:95170        MPYSRESKEKETHPKETKQD
CP020335.1_2002383:2003390    MPRYSNSQKIPTQGTKNTRR
FM162591.1_3913240:3914247    MPRYANYQINPKQNIKNSHG
CP009539.1_1748617:1754976    MPYSNKSKKNEIRSEKSNEE

### Pattern finding
* Simply always predicting "negative" gives 83.4% accuracy
* Logistic regression on the number of each amino acid (disregarding position): 96.5% accuracy
* Autocorrelation, best combination: 90.8 %
* Logistic, regression on each AA on each position: 98.0% accuracy

Neither of these can be particular trusted: With 90% homology, the sequences are still to similar.
With a reasonably homology reduction (70%), there are only 6 positives and the models cannot be run.

The predicted optimal sequence is below, compared to the original AfpN18:
```
Optimal MPYYSNSKEKETHSKKNERD
        ||| | |||||||||  |||
AfpN18  MPYSSESKEKETHSKETERD
```
