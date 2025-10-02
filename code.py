# Installing Biopython 
!pip install biopython

from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Blast import NCBIWWW, NCBIXML

# Set your email for NCBI Entrez (mandatory)
Entrez.email = "examplmail@email.com"  # Replace with a valid email ID

# Fetching nucleotide sequence by accession number
accession = "NM_001301717"
handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

print("Accession:", record.id)
print("Organism:", record.annotations['organism'])

# Recording and Displaying DNA sequence
dna_seq = record.seq
print("\n DNA sequence:\n", dna_seq)

# Transcription of  mRNA sequence
mrna_seq = dna_seq.transcribe()
print("\n mRNA sequence:\n", mrna_seq)

# Creating a complementary strand 
complement_seq = dna_seq.complement()
print("\nComplementary strand:\n", complement_seq)

# Translation of sequence (trim to multiple of 3)
protein_seq = dna_seq.translate()
print("\nTranslated protein sequence:\n", protein_seq)

# GC Content calculation
gc_percent = round(gc_fraction(dna_seq) * 100, 2)
print("\nGC Content (%):", gc_percent)

# Performing BLASTn
print("\nRunning BLASTn on sequence. **************Please wait, this may take some time *************")
result_handle = NCBIWWW.qblast("blastn", "nt", str(dna_seq))

with open("blastn_result.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
print("BLASTn results saved to blastn_result.xml")

# Parsing BLASTn XML result and display top 3 alignments
def print_top_alignments(blast_xml_file, top_n=3):
    with open(blast_xml_file) as result_handle:
        blast_record = NCBIXML.read(result_handle)
        print(f"\nTop {top_n} Alignments for BLASTn:\n")
        
        count = 0
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print("Sequence:", alignment.title)
                print("Length:", alignment.length)
                print(f"Score: {hsp.score}, E-value: {hsp.expect}")
                print("Query:", hsp.query)
                print("Match:", hsp.match)
                print("Subject:", hsp.sbjct)
                print()
                count += 1
                if count == top_n:
                    return

print_top_alignments("blastn_result.xml", top_n=3)
