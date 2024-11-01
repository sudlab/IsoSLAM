# slam_3UIs

This pipeline is meant to document the development of our SLAMseq read analysis.

Therefore it should maintain the outputs as things are more refined.

For the final release we will just use a refined pipeline (i.e. we will remove bits dealing with unstranded / non SNP filtered)

## Installation

You need a virtual environment with (at least) the following packages installed...

```bash
pip install cgatcore gevent sqlalchemy apsw
```

## Usage

The `pipeline_slam_3UIs/pipeline.yaml` file needs to be edited to reflect the location of key files. Once this reflects
your local file system you can run the following.

```bash
python -m pipeline_slam_3UIs --local make full -v5
```
