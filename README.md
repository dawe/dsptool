# DSPtool

# Table of Contents

# Usage

#### Required Extensions:
 There is a list of required packages to run the DSPtool, 
 It is suggested that use a virtual environment to prevent any possible conflict between packages. It is possible to download the packages separately or use the "dsptool.yml" file with the following command.  
> conda env create -f dsptool.yml --name dsptool

#### Example Files
BigWig file   [ENCFF409PGL.bw](https://drive.google.com/open?id=1JYTv_Zj-M6xtNzed5Mk3saG6LMbDUu70)

Bed file   [ENCFF992MBC.bed](https://github.com/dawe/dsptool/blob/master/ENCFF992MBC.bed)

#### Basic usage is as follows:
python3.7  filter-module.py [-h] -i path/inputfile.bw -o path/output.bw [-f FILTER] [-s SIZE] [-r REGION] [-l INTERVAL] [-S STEP] [-e ENTIRE] [--list-filters] [-V VERSION]

```
>>> python3.7  filter-module.py -i ENCFF409PGL.bw -o MyResult.bw -l ENCFF992MBC.bed
```
