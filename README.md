# bTBtreereconstruction
## Brief description

R scripts used to simulate *Mycobacterium bovis* sequence data, then reconstruct and evaluate bTB transmission trees.
R scripts are divided into 9 categories.

1. Sequence simulation
2. Transmission tree reconstruction
3. Convergence check
4. Description of genetic data
5. Accuracy
6. Outbreak size
7. Host-species contribution
8. Presence of super-spreaders
9. Host-species of the index case

## Coded variables: R *vs.* article

| Value in R | Name in article |
| ----------- | ----------- |
| B1 | Reference transmission scenario |
| A1 | Dead-end host scenario |
| B2 | Badger index case scenario |
| S1 | Single-host system scenario |
| S4 | High mutation rate scenario |
| 1 | Reference sampling scheme |
| 2 | Temporal bias or T |
| 3 | Wild boar bias or SW |
| 4 | Badger bias or SB |
| 5 | T + SW |
| 6 | T + SB |

## Recurrent R objects

| Name in R | Object |
| ----------- | ----------- |
| samp | Name of the transmission scenario |
| Ttree | Reference transmission tree |
| j | Number of the reference transmission tree |
| sim |  ID of the reference transmission tree or samp and j |
| nb_scheme | Number of sampling schemes to consider, 1 or 6 |
