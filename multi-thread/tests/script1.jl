# 2 seqs, 8 threads, not try...finally

using FASTX,BioFetch,PyCall
@pyimport rdkit.Chem as Chem
@pyimport rdkit.Chem.Draw as Draw

using InformationMeasures,CPUTime

function fetchseqs()
    seqs=fetchseq("QIS30335.1 QMI90807.1")

    return seqs
end

function compute(seqs)
    fe = open("errors", "w")
    l= ReentrantLock()
    Threads.@threads for seq in seqs
        fs = FASTX.FASTA.sequence(seq)
        s = FASTX.FASTA.sequence(seq, 330:583)
        id = FASTX.FASTA.identifier(seq)
        desc = FASTX.FASTA.description(seq)

        # FASTA molecule
        #=
        try
            while !thing_we_are_waiting_for
                wait(c)
            end
        finally
            unlock(c)
        end
        =#
        lock(l)
        mol_fasta = Chem.MolFromFASTA("$s")  #mol_fasta = pycall(Chem.MolFromFASTA, Any, "$s")
        unlock(l)
        #wait(mol_fasta)

        if mol_fasta==nothing
            println("Error with sequence ID: $id")
            println("Full sequence is: $fs")
            println("Sequence from 330-583: $s")
            println("")
            write(fe, string(id, '\n'))
        else
            # Write FASTA file (FULL sequence)
            lock(l)
            rec = FASTX.FASTA.Record("$id $desc", fs)
            f = FASTX.FASTA.Writer(open("fasta/full/$id.fasta", "w"))
            write(f, rec)
            close(f)
            unlock(l)

            # Write FASTA file (sequence from 330:583)
            lock(l)
            rec = FASTX.FASTA.Record("$id:330-583 $desc", s)
            f = FASTX.FASTA.Writer(open("fasta/$id.fasta", "w"))
            write(f, rec)
            close(f)
            unlock(l)

            # Molecule img from FASTA
            lock(l)
            Draw.MolToFile(mol_fasta,"fasta/img/$id.png")
            unlock(l)

            # For 3D (not working yet!)
            #m2=Chem.AddHs(mol_fasta)
            #Draw.MolToImage(m2)
            #Draw.MolToFile(m2,"fasta/img/$id-3D.png")
            
            # FASTA to SMILES
            lock(l)
            smiles = Chem.MolToSmiles(mol_fasta)
            unlock(l)

            # Write SMI file
            lock(l)
            f = open("smiles/$id.smi", "w")
            write(f, smiles)
            close(f)
            unlock(l)

            # Molecule img from SMILES
            lock(l)
            mol_smile = Chem.MolFromSmiles(smiles)
            Draw.MolToFile(mol_smile,"smiles/img/$id.png")
            unlock(l)
        end
    end
    close(fe)
end

@time @CPUtime seqs=fetchseqs()
@time @CPUtime compute(seqs)
