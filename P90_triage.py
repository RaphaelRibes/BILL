from Bio import SeqIO

#This code is designed to extract only P90 FASTA sequences and organize them into two separate files based on their variant type.
#P90 sequences with the same length as the ORF99 reference sequence are classified as SNP variants and stored in the first file.
#P90 sequences with a different length than the ORF99 reference are classified as structural variants (SV) and stored in the second file.

def orf_99_sorting(fasta_file, snp_file, sv_file):

    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    #snp_file = "/home/farah/Bureau/BILL/SNP_0RF99.fasta"
    #sv_file = "/home/farah/Bureau/BILL/SV_ORF99.fasta"

    ref = sequences[0]
    ref_len = (len(ref.seq))

    with open(snp_file, "a") as snp_out, open(sv_file, "a") as sv_out:
        for seq in sequences:
            if "P90" in seq.name:
                if len(seq.seq) == ref_len:
                    SeqIO.write (seq, snp_out, "fasta")
                    print(f"{seq.id} saved in {snp_file}")
                else:
                    SeqIO.write(seq, sv_out, "fasta")
                    print(f"{seq.id} saved in {sv_file}")

    print("done")

#snp_file = "/home/farah/Bureau/BILL/SNP_0RF99.fasta" -> destinantion file for snp variants 
#sv_file = "/home/farah/Bureau/BILL/SV_ORF99.fasta" -> destination file for SV

orf_99_sorting("/home/farah/Bureau/BILL/ORF99_aa.fasta","/home/farah/Bureau/BILL/SNP_0RF99.fasta","/home/farah/Bureau/BILL/SV_ORF99.fasta")