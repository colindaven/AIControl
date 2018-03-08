####################
# Import libraries #
####################

using PureSeq
using Distributions
using JLD
using DataFrames

######################
# Defining functions #
######################

function checkControlAvailability(controlFiles, basedir, suffix, excludes=[])
    ret = String[]
    mask = Bool[]
    for c in controlFiles
        if isfile("$(basedir)$(c).f$(suffix)") && isfile("$(basedir)$(c).r$(suffix)")
            if !(c[1:end-4] in excludes)
                push!(ret, c)
                push!(mask, true)
            else
                push!(mask, false)
            end
        end
    end
    push!(mask, true)
    ret, mask
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

function filter2d(matrix, ymask, xmask)
    temp = filter(x->x[2], [(matrix[i,1:end], ymask[i]) for i in 1:size(matrix)[1]])
    temp = [i[1] for i in temp]
    temp = hcat(temp...)
    
    temp = filter(x->x[2], [(temp[i,1:end], xmask[i]) for i in 1:size(temp)[1]])
    temp = [i[1] for i in temp]
    temp = hcat(temp...)
end

function smooth(a, width)
    if width > length(a)
        width = length(a)
    end
    
    ret = copy(a)
    counts = ones(length(a))
    
    #Aggregating the values
    for i in 1:width
        for j in 1:length(a)-i
            ret[j] += a[j+i]
            ret[end-j+1] += a[end-j-i+1]
            counts[j] += 1
            counts[end-j+1] += 1
        end 
    end
    ret./counts
end

type PeakWriter
    #a output stream object to write to
    Outstream
    contigs::ReferenceContigs
    cur_ref::Int64
end

function PeakWriter(output_stream, contigs)
    sw = PeakWriter(output_stream, contigs, 1)
    sw
end

function WritePeak(sw::PeakWriter, binSize::Int64, binPosStart::Int64, binPosEnd::Int64, pval::Float64, fold)
    startPos = binSize*(binPosStart-1)+1
    endPos = binSize*binPosEnd
    
    while startPos > sw.contigs.offsets[sw.cur_ref]+sw.contigs.sizes[sw.cur_ref]
        sw.cur_ref += 1
    end
    
    RNAME = sw.contigs.names[sw.cur_ref]
    startPos = startPos-sw.contigs.offsets[sw.cur_ref]
    endPos = endPos-sw.contigs.offsets[sw.cur_ref]
    println(sw.Outstream, RNAME,"\t",startPos,"\t",endPos,"\t-\t-\t-\t", fold,"\t", pval,"\t")
end

function compute_weights(targetFile, controlFiles, ifForward, num_chroms, mask, basedir, binsize, xtxdir; verbose=true, blockSize=10000, blockLimit=Inf)
    
    ##############
    # PARAMETERS #
    ##############
    if num_chroms > length(ReferenceContigs_hg38.sizes) || num_chroms < 1
        num_chroms = length(ReferenceContigs_hg38.sizes)
    end
    
    training_limit = Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize)) # of bins used for training.
    suffix = "bin$((binsize))" # suffix for the binned file.
    if ifForward == true direction = ".f" else direction = ".r" end #.r or .f for file suffix.
    
    ####################
    # BUILDING READERS # 
    ####################  
    #Build control readers.
    readers, control_names = load_controls(controlFiles, basedir, direction, suffix)
    #Build target control reader.
    target_reader = BinnedReader(targetFile)
    
    #######################
    # PREPARING MATRICIES #
    #######################
    P = length(readers) + 1 # +1 is coming from the constant term
    XtX1 = readdlm("$(xtxdir)/xtx$(direction)e.dlm")
    XtX2 = readdlm("$(xtxdir)/xtx$(direction)o.dlm")
    
    XtX1 = filter2d(XtX1, mask, mask)
    XtX2 = filter2d(XtX2, mask, mask)
    
    Xty1 = zeros((P,1))
    Xty2 = zeros((P,1))
    
    assert(P==size(XtX1)[1])
    assert(P==size(XtX2)[2])
    
    #####################
    # COMPUTING WEIGHTS #
    #####################
    count = 0
    for (target,control) in zip(denseblocks([target_reader], blockSize), denseblocks(readers, blockSize, constantColumn=true, loop=true))
        count += 1
        #checking the termination condition
        if count > blockLimit break end
        if count*blockSize > training_limit break end
        
        #updating our estimates
        if count % 2 == 0
            Xty1 .+= control'*target
        else
            Xty2 .+= control'*target
        end
            
        if count % 1000 == 0 && verbose
            println("Processed $(count*blockSize) out of $training_limit bins... ($targetFile)")
        end
    end
    
    weights1 = inv(XtX1 + 0.00001*eye(P))*Xty1
    weights2 = inv(XtX2 + 0.00001*eye(P))*Xty2
    
    ################
    # CLEANING IOs #
    ################
    for i in length(readers) close(readers[i]) end
    close(target_reader)
    weights1, weights2, control_names 
end

function compute_fits(targetFile, controlFiles, weight1, weight2, ifForward, num_chroms, basedir, binsize; verbose=true, blockSize=100000, blockLimit=Inf)

    ##############
    # PARAMETERS #
    ##############
    if num_chroms > length(ReferenceContigs_hg38.sizes) || num_chroms < 1
        num_chroms = length(ReferenceContigs_hg38.sizes)
    end
    
    training_limit = Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize)) # of bins used for training.
    suffix = "bin$((binsize))" # suffix for the binned file.
    if ifForward == true direction = ".f" else direction = ".r" end #.r or .f for file suffix.
    
    ####################
    # BUILDING READERS # 
    ####################  
    #Build control readers.
    readers, control_names = load_controls(controlFiles, basedir, direction, suffix)
    target_reader = BinnedReader(targetFile)
    
    #########################
    # PREPAREING DATA ARRAY #
    #########################
    regression_fit = zeros(Float16, training_limit)

    #################
    # PRODUCING FIT #
    #################
    # use denseblocks to iterate over blocks of the target and control data
    count = 0
    binpos = 0
    for (target,control) in zip(denseblocks([target_reader], blockSize), denseblocks(readers, blockSize, constantColumn=true, loop=true))
        count += 1
        #checking the termination condition
        if count > blockLimit break end
        if count*blockSize > training_limit break end
        
        #producing regression fit. 
        if count % 2 == 0 pred = control*weight2 else pred = control*weight1 end
        
        for j in 1:length(pred)
            binpos +=1
            try 
                regression_fit[binpos] = pred[j]
            end
        end   
        
        if count % 100 == 0 && verbose
            println("Processed $(count*blockSize) out of $training_limit bins... ($targetFile)")
        end
    end
    
    regression_fit
end

function callpeaks(targetFile, regfit, ifForward, num_chroms, binsize; base=1, smoothing=false)
       
    if ifForward == true direction = ".f" else direction = ".r" end #.r or .f for file suffix.
    if num_chroms > length(ReferenceContigs_hg38.sizes) || num_chroms < 1
        num_chroms = length(ReferenceContigs_hg38.sizes)
    end
    
    #loading signals
    b = BinnedReader(targetFile)
    target = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize)))
    while !eof(b) && position(b) < length(target)
        target[position(b)] = value(b)
        advance!(b)
    end
    close(b)

    #smoothing ?
    if smoothing
        smooth1000 = smooth(regfit, 10)
        smooth5000 = smooth(regfit, 50)
        smooth10000 = smooth(regfit, 100)
    end
    m = mean(target)
    
    pvals = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize)))
    lambdas = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize)))
    folds = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes[1:num_chroms])/binsize)))

    for i in 1:length(regfit)
        if smoothing
            lambda = maximum([regfit[i], smooth1000[i], smooth5000[i], smooth10000[i], base])
        else
            lambda = maximum([regfit[i], base])
        end
        pval = -1*log(10, 1-cdf(Poisson(lambda), target[i]))
        fold = target[i]/lambda
        if pval == Inf pval=1000.0 end

        #This version of the code will assign pval to individual bins.
        pvals[i] = pval
        lambdas[i] = lambda
        folds[i] = fold
    end
    
    pvals, folds, lambdas, target
end

function generate_peakfile(forward_p, reverse_p, folds_f, folds_r, name, th, offset, binsize)
    pvals = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes)/binsize)))
    folds = zeros(Float32, Int(ceil(sum(ReferenceContigs_hg38.sizes)/binsize)))
    
    for i in 1:length(forward_p)
        try
            pvals[i+Int(floor(offset/binsize))] = forward_p[i]
            folds[i+Int(floor(offset/binsize))] = folds_f[i]
        end
    end
    
    for i in 1:length(reverse_p)
        try
            if reverse_p[i] < pvals[i-Int(floor(offset/binsize))] 
                pvals[i-Int(floor(offset/binsize))] = reverse_p[i]
            end
            if folds_r[i] < folds[i-Int(floor(offset/binsize))] 
                folds[i-Int(floor(offset/binsize))] = folds_r[i]
            end
        end
    end
    
    fout = open("$(name).peaks","w")
    pw = PeakWriter(fout, ReferenceContigs_hg38)
    for i in 1:length(pvals)
        #This version of the code will assign pval to individual bins.
        pval = pvals[i]
        if pval >= th
            WritePeak(pw, binsize, i, i, Int(ceil(pval*1000))/1000.0, folds[i]) 
        end
    end
    close(fout)
    pvals
end

###################
# Load parameters #
###################

# Checking target file existance
targetfilename = ARGS[1]
if isfile(targetfilename) 
    println("Working on: ", targetfilename)
else
    println(targetfilename, " does not exist.")
    return -1
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
binname_f = "$(targetfilename).fbin$(binsize)"
binname_r = "$(targetfilename).rbin$(binsize)"
prefix = split(targetfilename, ".")[1]
chrom = parse(Int, configs["chrom"])

##################
# Peak Generator #
##################

println("Binning target file...")
if !isfile(binname_f)
    println("Binning $(targetfilename) (Forward)")
    write_binned(targetfilename, parse(Int, configs["binsize"]), :forward, skipDup=true)
end
if !isfile(binname_r)
    println("Binning $(targetfilename) (Reverse)")
    write_binned(targetfilename, parse(Int, configs["binsize"]), :reverse, skipDup=true)
end

metadata = readtable(configs["metafile"])
basectrls = collect(metadata[(metadata[:IFCTRL].== true), :ID])
temp = ["$(i).bam" for i in basectrls]

########################################################
#ToDo: Add check to see if you need to recompute xtxs. #  
########################################################
controlfiles, mask = checkControlAvailability(temp, configs["ctrldir"], "bin$(binsize)")

println("Computing weights...")
if !isfile("$(prefix)_wfe.dlm")
    wfe, wfo = compute_weights(binname_f, controlfiles, true, chrom, mask, configs["ctrldir"], parse(Int, binsize), configs["xtxdir"])
    writedlm("$(prefix)_wfe.dlm", wfe)
    writedlm("$(prefix)_wfo.dlm", wfo)
else
    wfe = readdlm("$(prefix)_wfe.dlm")
    wfo = readdlm("$(prefix)_wfo.dlm")
end

if !isfile("$(prefix)_wre.dlm")
    wre, wro = compute_weights(binname_r, controlfiles, false, chrom, mask, configs["ctrldir"], parse(Int, binsize), configs["xtxdir"])
    writedlm("$(prefix)_wre.dlm", wre)
    writedlm("$(prefix)_wro.dlm", wro)
else
    wre = readdlm("$(prefix)_wre.dlm")
    wro = readdlm("$(prefix)_wro.dlm")
end

println("Learning background noise...")
if !isfile("$(prefix)_fitf.dlm")
    fitf = compute_fits(binname_f, controlfiles, wfe, wfo, true, chrom, configs["ctrldir"], parse(Int, binsize))
    writedlm("$(prefix)_fitf.dlm", fitf)
else
    fitf = readdlm("$(prefix)_fitf.dlm")
end

if !isfile("$(prefix)_fitr.dlm")
    fitr = compute_fits(binname_r, controlfiles, wre, wro, false, chrom, configs["ctrldir"], parse(Int, binsize))
    writedlm("$(prefix)_fitr.dlm", fitr)
else
    fitr = readdlm("$(prefix)_fitr.dlm")
end

println("Identifiying peaks...")
pvals_f, folds_f, dmy, dmy = callpeaks(binname_f, fitf, true, chrom, parse(Int, binsize))
pvals_r, folds_r, dmy, dmy = callpeaks(binname_r, fitr, false, chrom, parse(Int, binsize))

generate_peakfile(pvals_f, pvals_r, folds_f, folds_r, "$(prefix)",
                    parse(Float64, configs["th"]),
                    parse(Int, configs["shift"]),
                    parse(Int, binsize))

println("Created peak file at $(prefix).peaks ...")
println("Done...")