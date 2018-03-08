# AIControl

Simulates control background experiments with machine learning from ENCODE control data.

## Required libraries for AIControl

Use ```pkg.add()``` and ```pkg.clone()``` to install libraries.
- DataFrames
- JLD
- Distributions
- PureSeq (https://github.com/slundberg/PureSeq.jl)

## Usage
1a. Download binned control data from our project website.  
1b. Alternatively, you can bin additional control data, if you have your own control data or want finer resolution bins.  
2. Set parameters appropriately in ```AIControl.config```

julia aicontrol.jl _bamfile_

## Parameters
You can edit aicontrol.config to change parameters of peak calling
