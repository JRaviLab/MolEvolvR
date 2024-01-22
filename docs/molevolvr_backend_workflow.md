``` mermaid
---
title: MolEvolvR workflow
---

flowchart TD
    %% Input types
    A(User input)
    A --> Z
    Z{Select input type}
    subgraph Input types
        B(FASTA)
        C(ACCNUM)
        D(MSA)
        E(BLAST output)
        F(Interproscan output)
    end
        Z --> B
        Z --> C
        Z --> D
        Z --> E
        Z --> F
    %% Options for each input type
        subgraph Advanced options
            G(Phylogenetic analysis & Domain architecture)
            H(Homlogy search & Domain architecture)
            I(Domain architecture ONLY)
            J(Homology search ONLY)
        end
        %%Y{Advanced options}
        %% FASTA
        B --> G
        B --> H
        B --> I
        B --> J
        %% ACCNUM
        C --> G
        C --> H
        C --> I
        C --> J
        %% MSA
        D --> G
        D --> H
        D --> I
        D --> J
        %% Blast output
        E --> G
        %% Interproscan output
        F --> H
        F --> I
        %% Processes
        subgraph Processes
            K(Deltablast)
            L(Deltablast cleanup)
            M(Interproscan)
            N(ipr2lin)
            O(ipr2da)
            P(blastclust)
            Q(clust2table)
        end
        %% Options' processes
            %% Homology Search & Domain Architecture
            H --> K
            H --> L
            H --> P
            H --> Q
            %% Homology search only
            J --> K
            J --> L
            J --> P
            J --> Q
            %% Phylogenetic analysis & Domain Architecture
            G --> M
            G --> N
            G --> O
            G --> P
            G --> Q
            %% Domain Architecture ONLY
            J --> M
            J --> N
            J --> O
            J --> P
            J --> Q

```
