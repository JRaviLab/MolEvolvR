# MolEvolvR diagram for mapping input types and advanced options to backend processes

## How to export as PNG/SVG?

The simplest way is to copy and paste the `mermaid` code block into mermaidchart.com and export.

Mermaid CLI is not used because I/we have not installed the [node module](https://github.com/mermaid-js/mermaid-cli/tree/master#mermaid-cli) 
on the server hosting the MolEvolvR application.

```mermaid
---
title: MolEvolvR proceses to run for each input and advanced option
---
flowchart LR

    SELECT_INPUT{Select input type}
    subgraph "Input types"
        subgraph "Sequence-like inputs"
            FASTA(FASTA)
            MSA(MSA)
        end
        ACCNUM(ACCNUM)

        subgraph "Tabular inputs"
            BLAST_OUTPUT(BLAST output)
            INTERPROSCAN_OUTPUT(Interproscan output)
        end
    end

        SELECT_INPUT --> FASTA
        SELECT_INPUT --> ACCNUM
        SELECT_INPUT --> MSA
        SELECT_INPUT --> BLAST_OUTPUT
        SELECT_INPUT --> INTERPROSCAN_OUTPUT
    
    subgraph "Advanced options "
        PHYLO_AND_DOMARCH(Phylogenetic analysis & Domain architecture)
        HOMOLOGY_AND_DOMARCH(Homlogy search & Domain architecture)
        DOMARCH(Domain architecture ONLY)
        HOMOLOGY(Homology search ONLY)
    end

        %% Input type to options
        FASTA --> PHYLO_AND_DOMARCH
        FASTA --> HOMOLOGY_AND_DOMARCH
        FASTA --> DOMARCH
        FASTA --> HOMOLOGY

        ACCNUM --> PHYLO_AND_DOMARCH
        ACCNUM --> HOMOLOGY_AND_DOMARCH
        ACCNUM --> DOMARCH
        ACCNUM --> HOMOLOGY

        MSA --> PHYLO_AND_DOMARCH
        MSA --> HOMOLOGY_AND_DOMARCH
        MSA --> DOMARCH
        MSA --> HOMOLOGY
        
        BLAST_OUTPUT --> PHYLO_AND_DOMARCH
        
        INTERPROSCAN_OUTPUT --> HOMOLOGY_AND_DOMARCH
        INTERPROSCAN_OUTPUT --> DOMARCH
        
    subgraph "Processes "
        DELTABLAST(Deltablast)
        DELTABLAST_CLEANUP(Deltablast cleanup)
        INTERPROSCAN(Interproscan)
        IPR2LIN(ipr2lin)
        IPR2DA(ipr2da)
        BLASTCLUST(blastclust)
        CLUST2TABLE(clust2table)
    end

        %% Options to processes
        HOMOLOGY_AND_DOMARCH --> DELTABLAST
        HOMOLOGY_AND_DOMARCH --> DELTABLAST_CLEANUP
        HOMOLOGY_AND_DOMARCH --> BLASTCLUST
        HOMOLOGY_AND_DOMARCH --> CLUST2TABLE

        HOMOLOGY --> DELTABLAST
        HOMOLOGY --> DELTABLAST_CLEANUP
        HOMOLOGY --> BLASTCLUST
        HOMOLOGY --> CLUST2TABLE

        PHYLO_AND_DOMARCH --> INTERPROSCAN
        PHYLO_AND_DOMARCH --> IPR2LIN
        PHYLO_AND_DOMARCH --> IPR2DA
        PHYLO_AND_DOMARCH --> BLASTCLUST
        PHYLO_AND_DOMARCH --> CLUST2TABLE

        DOMARCH --> INTERPROSCAN
        DOMARCH --> IPR2LIN
        DOMARCH --> IPR2DA
        DOMARCH --> BLASTCLUST
        DOMARCH --> CLUST2TABLE

```
