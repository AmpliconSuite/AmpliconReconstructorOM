# AmpliconReconstructor (AR)
Reconstructs complex variation using BioNano optical mapping data and an NGS-derived breakpoint graph. A preprint of the publication associated with this work is available on bioRxiv:

["AmpliconReconstructor: Integrated analysis of NGS and optical mapping resolves the complex structures of focal amplifications in cancer"](https://www.biorxiv.org/content/10.1101/2020.01.22.916031v1)

## Contents:
1. [Dependencies](#dependencies)
2. [Installation](#installation)
3. [Inputs & Outputs](#inputs-&-outputs)
4. [Running AR](#usage)
5. [SegAligner documentation](#segaligner-documentation)

## Dependencies
AR uses Python 2.7 and C++11 (with g++ as the compiler) and a Unix-based OS. Both of these are standard in any modern version of Linux. AR has been tested on Ubuntu 16.04 and Ubuntu 18.04.

A version of Python including commonly used scientific libraries (such as Anaconda) will be sufficient as the only Python dependencies are
 - Matplotlib 2.0.0 (or higher)
    * To ensure you meet the version dependency for Matplotlib, do
`pip install --upgrade matplotlib`.
 - numpy

AR can produce optional visualizations of the reconstructed amplicons, this requires the 
[CycleViz](https://github.com/jluebeck/CycleViz). Instructions for installing CycleViz are included below, but you may want to satisfy the dependencies now.
 
- Python `intervaltree` library
    * `pip install intervaltree`     
- ttf-mscorefonts library (optional, but recommened)
    * `sudo apt-get install ttf-mscorefonts-installer`

## Installation
To install AR, we add some variables to the .bashrc file (located in your home directory). We provide some bash commands to automate this process. Typically this installation process can be completed in 5 or less minutes.

1. Download AR.

    `git clone https://github.com/jluebeck/AmpliconReconstructor && cd AmpliconReconstructor`

2. Add AmpliconReconstructor and SegAligner to path.

    `echo "export AR_SRC=$PWD" >> ~/.bashrc`
    
    `echo "export SA_SRC=$PWD/SegAligner" >> ~/.bashrc`

3. A SegAligner binary compatible with a Linux x86 64-bit architecture is included. If you want to compile the SegAligner binary yourself, run the following.

    ```
    cd SegAligner
    make 
    cd ..
    ```

4. Add AR python libs to `$PYTHONPATH` variable.

    `echo "export PYTHONPATH=$PYTHONPATH:$AR_SRC >> ~/.bashrc"`

5. *(Optional, but highly recommended)* Install [CycleViz](https://github.com/jluebeck/CycleViz).
   
    ```
   # before running, make sure dependencies for CycleViz (listed in "Dependencies" section) are satisfied
   # now, install CycleViz
   cd .. #or wherever you want to install CycleViz
   git clone https://github.com/jluebeck/CycleViz
   echo "export CV_SRC=$PWD/CycleViz" >> ~/.bashrc
   ```
   
6. Make the changes to `.bashrc` live for this session.

    `source ~/.bashrc`


## Inputs & Outputs
### Inputs:

At a high level, AR accepts as inputs assembled Bionano contigs (`.cmap`) and an AA-formatted breakpoint graph file (`*_graph.txt`). 

##### Generating an *in silico* digested graph
Once AA has been run on your NGS data, please convert the resulting `*_graph.txt` file to CMAP form, using the script `generate_cmap.py`.

An example command is

`$AR_SRC/generate_cmap.py -r [path_to_reference.fasta] -e [enzyme] -g [path/to/your_graph.txt]`

##### The YAML file

A sample .yaml template is included in the AR source directory.

The .yaml file should specify the following properties for each entry (`sample_name`)


```
sample_name:
	path: /some/path/to/your/samples/      - Path prefix which will be applied to all other input filenames

	graph: sample_graph.txt                - AA-formatted breakpoint graph file. Assumes "graph" is located under "path".

	contigs: sample_EXP_REFINEFINAL1.cmap  - assembled OM contigs. Assumes "contigs" is located under "path".

	cmap: sample_graph_DLE1.cmap           - in-silico CMAP generated from AA-formatted graph file

	instrument: [Irys/Saphyr]

	enzyme: [BspQI/DLE1]

	min_map_len: ~                         - [Optional] set a custom minimum number of labels for SegAligner alignment (advanced option)

	xmaps: ~                               - [Optional] provide a .xmap-formatted file of alignments between "cmap" and reference genome. Disables SegAligner alignment (advanced option)
```

### Outputs:
AR outputs a collection of possible reconstruction paths, ordered by total alignment score. 
The output file format is the same as the [AA cycles file format](https://github.com/jluebeck/AmpliconArchitect/blob/master/README.md#file-formats).

The files created by AR are placed into three folders

1. `alignments/`
    * Contains the individual graph segment alignments to contigs, produced by SegAligner
2. `reconstructions/`
    * AR reconstructions and reconstruction alignments
3. `visualizations/`
    * Visualizations of AR reconstructions
    
The `reconstructions/` folder contains a number of files.

a) `[sample_name]_paths_cycles.txt` - AA cycles-formatted output describing the reconstructed genomic paths. This is the primary output used for interpreting reconstructions.

b) `[sample_name]_path_[N]_aln.txt` - SegAligner-formatted OM alignment of entire reconstruction path (segments to scaffolds).

c) `[sample_name]_scaffold_paths.txt` - AA cycles-formatted output describing the heaviest weight paths for each individual OM contig.

d) `[sample_name]_scaffold_path_[N]_aln.txt` - SegAligner-formatted OM alignment of graph segments to individual OM contig.

e) `[sample_name]_data.json` - A `json` representation of the unresolved reconstruction graph. Can be uploaded to [ScaffoldGraphViewer](https://jluebeck.github.io/ScaffoldGraphViewer/).

f)  `[sample_name]_run.log` - Verbose output about reconstruction process.

  

## Usage
AR requires a number of inputs. To simplify running AR, we use a wrapper script `AmpliconReconstructorOM.py`, and we use a YAML file to specify the sample-specific files and paths AR will use. Information about creating the YAML file with your inputs (including a template) is located in the section "Preparing your files".

An example invokation of AR is as follows

`python AmpliconReconstructorOM.py -i samples.yaml --outdir /path/to/output/directory --run_name [name for this run] --threads 24 `

A descripton of other command line arguments can produced by running AR with the `--help` flag.

### Running an AR test
You can test AR on using previously published data. The GBM39 cell line has been previously characterized by AR in [Wu, et al., <em>Nature</em>](https://www.nature.com/articles/s41586-019-1763-5). We have made WGS data and OM data publicly available:

- [GBM39 WGS data](https://www.ncbi.nlm.nih.gov/sra/SRX2006441[accn])

- [GBM39 OM data](https://submit.ncbi.nlm.nih.gov/subs/supfiles/SUB6698144/overview) (.cmap file)

You may either generate an AA breakpoint graph from the WGS data yourself (see [PrepareAA](https://github.com/jluebeck/PrepareAA) for details) or we also provide the pre-generated GBM39 AA breakpoint graph file in the AR data repo (`AmpliconReconstructor/test_files/`). 

- **If starting from BAM file:** 

If you wish to generate the AA breakpoint graph from scratch, you can use PrepareAA on the downloaded BAM file from SRA. An example command is below. This may take 1-2 hours on a standard desktop.

`/path/to/PrepareAA/PrepareAA.py -s GBM39  -t 8 --cnvkit_dir /path/to/cnvkit.py --rscript_path /path/to/Rscript --sorted_bam FF18.cs.bam --run_AA`

This will output a file `GBM39_amplicon1_graph.txt` which can be used in the next step.

- **If starting from OM data + pre-generated graph:** 

Starting from the pre-generated input data should run on a standard desktop in less than 10-15 minutes.

An example .yaml file is provided in the test_files directory. By default AR will not produce CycleViz images unless either the `CV_SRC` bash variable is set (CycleViz installation) or `--CV_path /path/to/CycleViz/` is specified.

To run the test, do 

`cd test_files`

`python $AR_SRC/AmpliconReconstructorOM.py -i gbm39_test.yaml --outdir GBM39_AR_output --run_name GBM39_test --nthreads [num threads]`

In the resulting folder you will see AR output files, including reconstructed cycles, and (if CycleViz specified) visualization plots. Furthermore, the SegAligner and combined alignment for the entire amplicon will be present. The SegAligner files end with `*_aln.txt`. 

## SegAligner
SegAligner is a multi-threaded C++ aligner for BioNano optical map contigs and in-silico digested genomic reference segments. It additionally supports alignment of contigs to the full reference genome to identify the location of candidate regions of the contigs belonging to regions of the reference. 


SegAligner is wrapped inside AmpliconReconstructorOM, but can be invoked on its own. For installation please see the instructions for AR installation above.


Command line arguments for SegAligner are as follows - the first two arguments are positional and are paths to (1) the in-silico reference genome CMAP (full or collection of extracted segments) and (2) the BioNano contig CMAP. The most simple command would be:

`SegAligner your_reference_segs.CMAP EXP_REFINEFINAL1.cmap`

However, there are other arguments you will probably want to consider: 

- `-nthreads=1` Sets the number of threads to use (default 1, HIGHLY recommend at least 8 for typical datasets).

- `-prefix=SA_output` Sets the filename output prefix (default "SA_output")).

- `-min_labs=10` Sets the minimum number of labels on a cmap entry in order to attempt alignment (default 10. For Saphyr DLE1 data, recommend 12).

- `-n_detect_scores=500` Set the number of in the E-value distribution (default 500).

- `-gen=2` Set the "generation" of BioNano instrument. "1" = Irys, "2" = Saphyr (default 2).

- `-local` Perform local alignment (default semi-global).

- `-fitting_aln` Performs fitting alignment (conflicts with `-local`). Not recommended for large datasets.

- `-no_tip_aln` Turns off sensitive search for overlapping alignments. 

- `-nl` Turns off banded alignment. Not recommended for large datasets (memory + speed issues).

- `-detection` Used for detection of contig locations in reference genome. Sets `-local_aln`, `-no_tip_aln`, and other internal variables.


