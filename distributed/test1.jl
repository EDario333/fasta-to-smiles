# 100 seqs, 8 threads, with try... finally

using Distributed

addprocs([("edario@192.168.0.11", :auto), ("edario@192.168.0.200", :auto)], exename="/home/edario/julia-1.6.3/bin/julia")

@everywhere using FASTX,BioFetch,PyCall
@everywhere @pyimport rdkit.Chem as Chem
@everywhere @pyimport rdkit.Chem.Draw as Draw

using InformationMeasures,CPUTime

function fetchseqs()
    seqs=fetchseq("QIS30335.1 QMI90807.1 QJS39579.1 QJS39567.1 QMS52716.1 QIS30425.1 QOU93974.1 QMI94525.1 QJC19455.1 QKU28906.1 YP_009724390.1 BCN86353.1 QJF75467.1 QJX45031.1 QJR85953.1 QII57278.1 QRN64146.1 QNN86157.1 7KRQ_A QJF77846.1 QRN78371.1 QMS52716.1 QIZ16509.1 QKU32813.1 QIZ97039.1 QJQ84843.1 QKS90791.1 QQP45825.1 QJG65956.1 QMJ01317.1 QOH28993.1 QPN00919.1 QOS49673.1 QNN86445.1 QKV39923.1 QNA39510.1 QIU81585.1 QMX86773.1 QJT73034.1 QIZ16559.1 QOF10625.1 QWC76579.1 QSE99577.1 QOF08429.1 QOE75686.1 QQ076328.1 QJD20632.1 QIS60546.1 QTI95854.1 QOS50952.1 QXE71677.1 QII57268.1 QIT07011.1 QWC79193.1 QMX85049.1 QIO04367.1 QJW70147.1 QKU28894.1 QXR88761.1 QOU94298.1 QJS53398.1 QOU91202.1 QQP47157.1 QLC93524.1 QWM87761.1 QXR89097.1 QOF12353.1 QOU89966.1 QQK89979.1 QHU79173.2 QII87830.1 QLC48052.1 QJE38426.1 QIS60489.1 QLI51781.1 QKV41471.1 QKY78084.1 QRA18933.1 QJR84873.1 QWC77885.1 QRN63738.1 QMT96172.1 QJQ39524.1 QIU80973.1 QOF15989.1 QRN58518.1 QJQ84676.1 QKU33257.1 QJR84837.1 QWC76135.1 QOU86714.1 QIS60906.1 QKU31481.1 QQP44733.1 QJF76007.1 QWJ89677.1 QKV35819.1 QPF48560.1 QIA98583.1 QIS30165.1")

    return seqs
end

@everywhere function compute(coll)
    fe = open("errors", "w")
    #Threads.@threads for seq in seqs
    master_task = @distributed for x in coll
        fs = FASTX.FASTA.sequence(global seqs[x])
        s = FASTX.FASTA.sequence(global seqs[x], 330:583)
        id = FASTX.FASTA.identifier(global seqs[x])
        desc = FASTX.FASTA.description(global seqs[x])
        println(fs)
        println(s)
        println(id)
        println(desc)
    end
    println("Beginning the process...")
    wait(master_task)
    println("All workers done!")
end

#@time @CPUtime 
#@time begin
#    seqs=fetchseqs()
#end
#@time @CPUtime 

pmap(compute, 1:length(global seqs))

#=
@time begin
    compute(seqs)
end
=#

#@everywhere 
rmprocs()
