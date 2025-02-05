# Workflow

This page gives an overview of the workflow undertaken by IsoSLAM, it is a **WORK IN PROGRESS** as the code base is
underoing refactoring.

- This is very much a work in progress and is not yet complete. Contributions are welcome.

```mermaid
flowchart TB
    subgraph Input
        BAM[(BAM Files)]
        VCF[("VCF Files")]
        GTF[("GTF Files")]
        Config[("Config Files")]
    end

    subgraph Core["Core Processing"]
        IO["Input/Output Handler"]:::core
        Process["Processing Engine"]:::core
        Pipeline["SLAM Pipeline"]:::core
        Summary["Summary Generator"]:::core
    end

    subgraph Support["Support Components"]
        Logger["Logging System"]:::support
        Utils["Utility Functions"]:::support
        DefaultConfig["Configuration Manager"]:::support
    end

    subgraph Testing["Testing & Documentation"]
        TestInfra["Testing Infrastructure"]:::test
        TestRes["Test Resources"]:::test
        APIDoc["API Documentation"]:::doc
    end

    subgraph Integration["External Integration"]
        RScript["R Script Integration"]:::integration
        CICD["CI/CD Pipeline"]:::integration
        PipeConfig["Pipeline Configuration"]:::integration
    end

    %% Main Data Flow
    BAM --> IO
    VCF --> IO
    GTF --> IO
    Config --> DefaultConfig

    IO --> Process
    Process --> Pipeline
    Pipeline --> Summary

    %% Support Flow
    DefaultConfig -.-> IO
    DefaultConfig -.-> Process
    DefaultConfig -.-> Pipeline

    Logger -.-> IO
    Logger -.-> Process
    Logger -.-> Pipeline
    Logger -.-> Summary

    Utils --> IO
    Utils --> Process
    Utils --> Pipeline

    %% Integration Flow
    Pipeline --> RScript
    PipeConfig -.-> Pipeline
    CICD -.-> TestInfra

    %% Testing Flow
    TestRes --> TestInfra
    TestInfra -.-> Core

    %% Click Events
    click IO "https://github.com/sudlab/IsoSLAM/blob/main/isoslam/io.py"
    click Process "https://github.com/sudlab/IsoSLAM/blob/main/isoslam/processing.py"
    click Pipeline "https://github.com/sudlab/IsoSLAM/tree/main/isoslam/pipeline_slam_3UIs/"
    click Summary "https://github.com/sudlab/IsoSLAM/blob/main/isoslam/summary.py"
    click Logger "https://github.com/sudlab/IsoSLAM/blob/main/isoslam/logging.py"
    click Utils "https://github.com/sudlab/IsoSLAM/blob/main/isoslam/utils.py"
    click DefaultConfig "https://github.com/sudlab/IsoSLAM/blob/main/isoslam/default_config.yaml"
    click TestInfra "https://github.com/sudlab/IsoSLAM/tree/main/tests/"
    click APIDoc "https://github.com/sudlab/IsoSLAM/tree/main/docs/api/"
    click CICD "https://github.com/sudlab/IsoSLAM/tree/main/.github/workflows/"
    click PipeConfig "https://github.com/sudlab/IsoSLAM/blob/main/isoslam/pipeline_slam_3UIs/pipeline.yml"
    click RScript "https://github.com/sudlab/IsoSLAM/blob/main/isoslam/pipeline_slam_3UIs/summarize_counts.R"
    click TestRes "https://github.com/sudlab/IsoSLAM/tree/main/tests/resources/"

    %% Styling
    classDef core fill:#90EE90,stroke:#333,stroke-width:2px
    classDef support fill:#FFE4B5,stroke:#333,stroke-width:2px
    classDef test fill:#DDA0DD,stroke:#333,stroke-width:2px
    classDef doc fill:#87CEEB,stroke:#333,stroke-width:2px
    classDef integration fill:#F08080,stroke:#333,stroke-width:2px

    %% Legend
    subgraph Legend
        L1["Core Components"]:::core
        L2["Support Components"]:::support
        L3["Testing Components"]:::test
        L4["Documentation"]:::doc
        L5["Integration"]:::integration
    end
```

The above diagram is written in [Mermaid][mermaid] and generated using [GitDiagram][gitdiagram]. You can view the source
code in the IsoSLAM repository and develop/modify it using the [Mermaid Live Editor][mermaid-live] and make
pull-requests to update this documentation.

## IsoSLAM

A number of pre-processing steps are undertaken prior to IsoSLAM work being done. The following is work in progress as
the code is refactored.

1. Iterate over `.bam` file and pair segments. If two or more `AlignedSegments` with the same `query_name` are found
   then `n > 1` segments are dropped.
2. Pairs of segments (individual `AlignedSegments`) are then assessed and if they are `Assigned` the `start`, `end`,
   `length`, `status` (i.e. `Assigned`), `transcript_id`, `block_start` and `block_end` are extracted.

## Descriptive Workflow

<!-- markdownlint-disable MD033 -->
<style>
body {
    counter-reset: h3counter;
}
h3 {
    counter-increment: h3counter;
}
h3:before {
    content: counter(h3counter) ". ";
}
</style>

### Read alignments are loaded

### The gene transcript these are within are identified

### Introns within these genes are identified

### Retain reads that are overlap with introns or splice ends

```text
                                                 Genome ----------------------->

Read Alignment Blocks                     |>>>>|                     |>>>>>>>>>>>>>>>|
Transcript 1:                      |===========|---------------------|=========|-------------------|==========|
Transcript 2:              |====|------------------------------------|=========|-------------------|==========|

Introns[0]                                      ---------------------
Introns[1]                                                                      -------------------
Introns[2]                       ------------------------------------
Introns[3]                                                                      -------------------

Exon : |========|   Read alignment block:  |>>>>>|
```

[gitdiagram]: https://gitdiagram.com/sudlab/IsoSLAM
[mermaid]: https://mermaid.js.org/
[mermaid-live]: https://mermaid.live
