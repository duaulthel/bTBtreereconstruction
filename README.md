# bTBtreereconstruction
## Brief description

R scripts used to simulate *Mycobacterium bovis* sequence data, then reconstruct and evaluate bTB transmission trees.
R scripts (in 0.Code) are divided into 9 categories.

1. Sequence simulation

| Input | Output |
| ----------- | ----------- |
| Reference trees (Ttree files) | Simulated sequences (seq files) |

2. Transmission tree reconstruction

| Input | Output |
| ----------- | ----------- |
| Reference trees | Parameters for *outbreaker2* and *TransPhylo* (generation/sampling time files) |
| Simulated sequences + Parameters | Transmission trees (seqtrack_tree or res_tree or (WIW)ttree files) |

3. Convergence check

| Input | Output |
| ----------- | ----------- |
| Output from 2. (res or mcmc files) | outbreaker/transphylo_mcmc files |

4. Description of genetic data

| Input | Output |
| ----------- | ----------- |
| Simulated sequences and reference trees | trans_div files and S1 Fig |

5. Accuracy

| Input | Output |
| ----------- | ----------- |
| Reference trees | Add transmission pairs that can be reconstructed (Ttree_det files) |
| Simulated sequences, reference and reconstructed trees | outbreaker/seqTrack/transphylo_acc files |
| outbreaker/seqTrack/transphylo_acc files | GLM(M) Binomial |

6. Outbreak size

| Input | Output |
| ----------- | ----------- |
| Reference trees | Subtree induced by the sampling scheme (Ttree_det_sup files) |
| Simulated sequences, reference and reconstructed trees | outbreaker/transphylo_size files |
| outbreaker/transphylo_size files | GLM(M) Negative Binomial or Figs 3, 5 and 7 |

7. Host-species contribution

| Input | Output |
| ----------- | ----------- |
| Simulated sequences, reference and reconstructed trees | outbreaker/transphylo_reff files |
| outbreaker/transphylo_reff files | GLM(M) Negative Binomial or Figs 4, 6 and 8 |

8. Presence of super-spreaders

| Input | Output |
| ----------- | ----------- |
| Simulated sequences, reference and reconstructed trees | outbreaker/seqTrack/transphylo_spp files |

9. Host-species of the index case

| Input | Output |
| ----------- | ----------- |
| Simulated sequences and reconstructed trees | outbreaker/seqTrack/transphylo_index files |
| outbreaker/seqTrack/transphylo_index files | GLM Binomial |

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
