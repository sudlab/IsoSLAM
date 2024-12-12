# Workflow

This page gives an overview of the workflow undertaken by IsoSLAM.

**WORK IN PROGRESS** - This is very much a work in progress and is not yet complete. Contributions are welcome.

```{mermaid}
%%{init: {'theme': 'base', 'gitGraph': {'rotateCommitLabel': true}
         }
}%%
graph TD;

  subgraph input
  A1([STAR Aligned Reads]) --> B1
  A2([.FASTA file]) --> B1
  A3([Contigs]) --> B1
  A4([.BED file]) --> B1
  A5([.GTF file]) --> B1
  end

  subgraph Python

  B1[pipeline.py] --> B2


  style A1 fill:#648FFF,stroke:#000000
  style A2 fill:#648FFF,stroke:#000000
  style A3 fill:#648FFF,stroke:#000000
  style A4 fill:#648FFF,stroke:#000000
  style A5 fill:#648FFF,stroke:#000000
  end

  B2([Pipeline_slam_3UIs.py]) --> statistics


  subgraph statistics
  D1([get_pair_pvalues])
  D2([get_pair_half_lives])
  D3([get_over_time_pvalues])
  D4([get_interaction_over_time_pvals])

  end
  style B1 fill:#00FF90,stroke:#000000
  style B2 fill:#00FF90,stroke:#000000





  style D1 fill:#FE6100,stroke:#000000
  style D2 fill:#FE6100,stroke:#000000
  style D3 fill:#FE6100,stroke:#000000
  style D4 fill:#FE6100,stroke:#000000

  style E1 fill:#FEE100,stroke:#000000

```

## IsoSLAM

A number of pre-processing steps are undertaken prior to IsoSLAM work being done. The following is work in progress as
the code is refactored.

1. Iterate over `.bam` file and pair segments. If two or more `AlignedSegments` with the same `query_name` are found
   then `n > 1` segments are dropped.
2. Pairs of segments (individual `AlignedSegments`) are then assessed and if they are `Assigned` the `start`, `end`,
   `length`, `status` (i.e. `Assigned`), `transcript_id`, `block_start` and `block_end` are extracted.

## Updating

The above diagram is written in [Mermaid][mermaid]. You can view the source code in the IsoSLAM repository and
develop/modify it using the [Mermaid Live Editor][mermaid-live] and make pull-requests to update this documentation.

[mermaid]: https://mermaid.js.org/
[mermaid-live]: https://mermaid.live
