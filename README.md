# AIControl

AIControl simulates/imputes control experiments and identify peaks for your ChIPseq data.
![alt text](images/concept.png)

## Required libraries for AIControl

Use ```pkg.add()``` and ```pkg.clone()``` to install libraries.
- DataFrames
- JLD
- Distributions
- PureSeq (https://github.com/slundberg/PureSeq.jl)

## Usage
1. Download binned control data from the link on our project website (http://suinlee.cs.washington.edu/projects/aicontrol/). Alternatively, you can bin your own data using ```WriteBinned``` function in the PureSeq package. We are working on the better compression of binned data.  
2. Set parameters appropriately in ```AIControl.config```. Particularly, ```ctrldir```, ```xtxdir```, and ```metafile``` is important. 
3. Call ```julia aicontrol.jl _your_bamfile_```

For questions, please e-mail hiranumn at cs dot washington dot edu
