# 100 seqs, single thread, with try... finally

using FASTX,BioFetch,PyCall
@pyimport rdkit.Chem as Chem
@pyimport rdkit.Chem.Draw as Draw

using InformationMeasures,CPUTime

function fetchseqs()
    seqs=fetchseq("QIS30335.1 QMI90807.1 QJS39579.1 QJS39567.1 QMS52716.1 QIS30425.1 QOU93974.1 QMI94525.1 QJC19455.1 QKU28906.1 YP_009724390.1 BCN86353.1 QJF75467.1 QJX45031.1 QJR85953.1 QII57278.1 QRN64146.1 QNN86157.1 7KRQ_A QJF77846.1 QRN78371.1 QMS52716.1 QIZ16509.1 QKU32813.1 QIZ97039.1 QJQ84843.1 QKS90791.1 QQP45825.1 QJG65956.1 QMJ01317.1 QOH28993.1 QPN00919.1 QOS49673.1 QNN86445.1 QKV39923.1 QNA39510.1 QIU81585.1 QMX86773.1 QJT73034.1 QIZ16559.1 QOF10625.1 QWC76579.1 QSE99577.1 QOF08429.1 QOE75686.1 QQ076328.1 QJD20632.1 QIS60546.1 QTI95854.1 QOS50952.1 QXE71677.1 QII57268.1 QIT07011.1 QWC79193.1 QMX85049.1 QIO04367.1 QJW70147.1 QKU28894.1 QXR88761.1 QOU94298.1 QJS53398.1 QOU91202.1 QQP47157.1 QLC93524.1 QWM87761.1 QXR89097.1 QOF12353.1 QOU89966.1 QQK89979.1 QHU79173.2 QII87830.1 QLC48052.1 QJE38426.1 QIS60489.1 QLI51781.1 QKV41471.1 QKY78084.1 QRA18933.1 QJR84873.1 QWC77885.1 QRN63738.1 QMT96172.1 QJQ39524.1 QIU80973.1 QOF15989.1 QRN58518.1 QJQ84676.1 QKU33257.1 QJR84837.1 QWC76135.1 QOU86714.1 QIS60906.1 QKU31481.1 QQP44733.1 QJF76007.1 QWJ89677.1 QKV35819.1 QPF48560.1 QIA98583.1 QIS30165.1")

    return seqs
end

function compute(seqs)
    fe = open("errors", "w")
    for seq in seqs
        fs = FASTX.FASTA.sequence(seq)
        s = FASTX.FASTA.sequence(seq, 330:583)
        id = FASTX.FASTA.identifier(seq)
        desc = FASTX.FASTA.description(seq)

        # FASTA molecule
        mol_fasta = Chem.MolFromFASTA("$s")  #mol_fasta = pycall(Chem.MolFromFASTA, Any, "$s")

        if mol_fasta==nothing
            println("Error with sequence ID: $id")
            println("Full sequence is: $fs")
            println("Sequence from 330-583: $s")
            println("")
            write(fe, string(id, '\n'))
        else
            # Write FASTA file (FULL sequence)
            rec = FASTX.FASTA.Record("$id $desc", fs)
            f = FASTX.FASTA.Writer(open("fasta/full/$id.fasta", "w"))
            write(f, rec)
            close(f)

            # Write FASTA file (sequence from 330:583)
            rec = FASTX.FASTA.Record("$id:330-583 $desc", s)
            f = FASTX.FASTA.Writer(open("fasta/$id.fasta", "w"))
            write(f, rec)
            close(f)

            # Molecule img from FASTA
            Draw.MolToFile(mol_fasta,"fasta/img/$id.png")

            # For 3D (not working yet!)
            #m2=Chem.AddHs(mol_fasta)
            #Draw.MolToImage(m2)
            #Draw.MolToFile(m2,"fasta/img/$id-3D.png")
            
            # FASTA to SMILES
            smiles = Chem.MolToSmiles(mol_fasta)

            # Write SMI file
            f = open("smiles/$id.smi", "w")
            write(f, smiles)
            close(f)

            # Molecule img from SMILES
            mol_smile = Chem.MolFromSmiles(smiles)
            Draw.MolToFile(mol_smile,"smiles/img/$id.png")
        end
    end
    close(fe)
end

@time @CPUtime seqs=fetchseqs()
@time @CPUtime compute(seqs)
