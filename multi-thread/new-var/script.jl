# 1 seq, 8 threads, with try... finally
# run: julia --thread 8 script.jl

using FASTX,BioFetch,PyCall
@pyimport rdkit.Chem as Chem
@pyimport rdkit.Chem.Draw as Draw

using InformationMeasures,CPUTime

function fetchseqs()
    seqs=fetchseq("B.1.1.529")

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
        mol_fasta = nothing
        try
            lock(l)
            mol_fasta = Chem.MolFromFASTA("$s")  #mol_fasta = pycall(Chem.MolFromFASTA, Any, "$s")
        finally
            unlock(l)
        end

        if mol_fasta==nothing
            println("Error with sequence ID: $id")
            println("Full sequence is: $fs")
            println("Sequence from 330-583: $s")
            println("")
            try
                lock(l)
                write(fe, string(id, '\n'))
            finally
                unlock(l)
            end
        else
            # Write FASTA file (FULL sequence)
            try
                lock(l)
                rec = FASTX.FASTA.Record("$id $desc", fs)
                f = FASTX.FASTA.Writer(open("fasta/full/$id.fasta", "w"))
                write(f, rec)
                close(f)
            finally
                unlock(l)
            end

            # Write FASTA file (sequence from 330:583)
            try
                lock(l)
                rec = FASTX.FASTA.Record("$id:330-583 $desc", s)
                f = FASTX.FASTA.Writer(open("fasta/$id.fasta", "w"))
                write(f, rec)
                close(f)
            finally
                unlock(l)
            end

            # Molecule img from FASTA
            try
                lock(l)
                Draw.MolToFile(mol_fasta,"fasta/img/$id.png")
            finally
                unlock(l)
            end

            # For 3D (not working yet!)
            #m2=Chem.AddHs(mol_fasta)
            #Draw.MolToImage(m2)
            #Draw.MolToFile(m2,"fasta/img/$id-3D.png")
            
            # FASTA to SMILES
            smiles = nothing
            try
                lock(l)
                smiles = Chem.MolToSmiles(mol_fasta)
            finally
                unlock(l)
            end

            # Write SMI file
            try
                lock(l)
                f = open("smiles/$id.smi", "w")
                write(f, smiles)
                close(f)
            finally
                unlock(l)
            end

            # Molecule img from SMILES
            try
                lock(l)
                mol_smile = Chem.MolFromSmiles(smiles)
                Draw.MolToFile(mol_smile,"smiles/img/$id.png")
            finally
                unlock(l)
            end
        end
    end
    close(fe)
end

@time @CPUtime seqs=fetchseqs()
@time @CPUtime compute(seqs)
