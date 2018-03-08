using DataFrames
using PureSeq
using JLD
using Distributions

# Check if the controls are available as binned data
function checkControlAvailability(controlFiles, basedir, suffix)
    ret = String[]
    for c in controlFiles
        if isfile("$(basedir)$(c).f$(suffix)") && isfile("$(basedir)$(c).r$(suffix)")
            push!(ret, c) 
        end
    end
    ret
end

# Load control files from specified directory. 
function load_controls(controlFiles, basedir, direction, suffix)
    readers = BinnedReader[]
    control_names = Any[]
    for i in 1:length(controlFiles)
        filename = "$(basedir)$(controlFiles[i])$(direction)$(suffix)"
        push!(readers, BinnedReader(filename))
        push!(control_names, filename)
    end
    readers, control_names
end

function compute_xtx(controlFiles, ifForward, num_chroms, basedir; verbose=false, blockSize=10000, blockLimit=Inf)
    
    ##############
    # PARAMETERS #
    ##############
    if num_chroms > length(ReferenceContigs_hg38.sizes) || num_chroms < 1
        num_chroms = length(ReferenceContigs_hg38.sizes)
    end
    
    binsize = 100 # number of bp in each bin.
    training_limit = Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize)) # of bins used for training.
    suffix = "bin100" # suffix for the binned file.
    if ifForward == true direction = ".f" else direction = ".r" end #.r or .f for file suffix.
    
    ####################
    # BUILDING READERS # 
    ####################  
    #Build control readers.
    readers, control_names = load_controls(controlFiles, basedir, direction, suffix)
    
    #######################
    # PREPARING MATRICIES #
    #######################
    P = length(readers) + 1 # +1 is coming from the bias term
    XtX1 = zeros(P,P)
    XtX2 = zeros(P,P)
    
    #####################
    # COMPUTING WEIGHTS #
    #####################
    count = 0
    for control in denseblocks(readers, blockSize, constantColumn=true, loop=true)
        count += 1
        #checking the termination condition
        if count > blockLimit break end
        if count*blockSize > training_limit break end
                
        #updating our estimates
        if count % 2 == 0
            XtX1 .+= control'*control
        else
            XtX2 .+= control'*control
        end
            
        if count % 1000 == 0 && verbose
            println("Processed $(count*blockSize) out of $training_limit bins...")
        end
    end

    writedlm("xtx/xtx$(direction)e.dlm",  XtX1)
    writedlm("xtx/xtx$(direction)o.dlm",  XtX2)
    
    ################
    # CLEANING IOs #
    ################
    for i in length(readers) close(readers[i]) end
end

# Load parameters
file = readlines(open("aicontrol.config"))
parameters = filter(x->!(length(x)==0||x[1]=='#'), [strip(i) for i in file])
configs = Dict()
for p in parameters
    line = [strip(i) for i in split(p, '=')]
    configs[line[1]] = line[2]
end

binsize = configs["binsize"]
chrom = parse(Int, configs["chrom"])

#metadata = readtable(configs["metafile"])
#basectrls = collect(metadata[(metadata[:IFCTRL].== true), :ID])
basectrls = collect(Set([split(i, ".")[1] for i in readdir(configs["ctrldir"])]))
temp = ["$(i).bam" for i in basectrls]
controlfiles = checkControlAvailability(temp, configs["ctrldir"], "bin$(binsize)")
writedlm("xtx/ctrllist.dlm", controlfiles)
println("Computing XtX matrices for $(length(controlfiles)) control files ...")

#compute_xtx(controlfiles, true, -1, configs["ctrldir"], verbose=true)
#compute_xtx(controlfiles, false, -1, configs["ctrldir"], verbose=true)