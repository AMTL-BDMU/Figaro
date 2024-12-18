### Running Amplicon-based ONT

```
nextflow run Figaro \
	--ontAmplicon \
	--inDir <input directory> \
	--outDir <output directory> \
	--reference <reference.fasta> \
	--primer <primer bed file> \
	--medakaModel <medaka model>
```
| Parameter | Description |
|-------| -----|
| inDir | directory containing sub-directories with non-concatenated fastq files or directory containing concatenated fastq files |
| reference | custom reference, e.g. whole genome or pol gene only |
| primer | bed file used for primer trimming. See sample `assets/primers/primerpair_pol.bed`. If primers are not within the reference fasta, just use `0` for the genomic range in the bed file. |
| medakaModel | see [list](https://github.com/nanoporetech/medaka/blob/master/medaka/options.py) |

<br> <br>

### Running Shotgun short-reads

```
nextflow run Figaro \
	--illuminaShotgun \
	--inDir <directory containing forward and reverse fastq> \
	--outDir <output directory>
```

| Parameter | Description |
|-------| -----|
| inDir | directory containing the forward and reverse fastq files |

<br> <br>

### Updating HIVDB Algorithm
Download xml files [here](https://hivdb.stanford.edu/page/algorithm-updates) then revise `--sierraXML` in the `nextflow.config`.
