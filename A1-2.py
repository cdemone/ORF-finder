# imports needed components of Biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable

# store bacterial translation table as a variable
Bacterial_Table = CodonTable.unambiguous_dna_by_name["Bacterial"]

# read sequence
for record in SeqIO.parse("BT_4072.fasta", "fasta"):

    # print descriptor line
    print(record.id)

    # store sequence as string
    fwd_sequence = str(record.seq)

    # print sequence length
    print(str(len(fwd_sequence)) + "bp")

    # convert sequence to unambiguous IUPAC notation
    formatted_seq = Seq(fwd_sequence, IUPAC.unambiguous_dna)

    # get reverse complement of sequence and store as string
    rev_comp_sequence = str(formatted_seq.reverse_complement())

    # find and translate ORF using given sequence and frame
    def ORF_finder(frame, sequence):

        # indicates if start codon has been located
        start_found = False

        while frame <= len(sequence):

            # look at sequence one codon at a time, beginning at specified reading frame
            codon = sequence[frame:frame + 3]

            # find start codon and indicate it has been located
            if codon in Bacterial_Table.start_codons and start_found is False:
                start_codon = frame
                start_found = True

            # after start codon has been located, find stop codon
            if codon in Bacterial_Table.stop_codons and start_found is True:
                stop_codon = frame

                # store and return ORF
                ORF = Seq(sequence[start_codon:stop_codon + 3], IUPAC.unambiguous_dna)
                return ORF

            # return message if no ORF is found
            if frame == len(sequence) + 1 and start_found is False:
                return 'No ORF found at this sequence position'

            # go to next codon
            else:
                frame += 3

    frame = 0
    ORF_list = []

    for frame in range(0, 3, 1):

        # tell user if there is no ORF at specified reading frame
        if ORF_finder(frame, fwd_sequence) == 'No ORF found at this sequence position':
            print("No Forward ORF " + str(frame))

        if ORF_finder(frame, rev_comp_sequence) == 'No ORF found at this sequence position':
            print("No Reverse ORF " + str(frame))

        # print DNA and amino acid sequence lengths for ORFs found
        # make a list of translated ORFs
        else:
            print("Forward ORF " + str(frame + 1) + ": " + str(len(ORF_finder(frame, fwd_sequence))) + "bp, "
                  + str(len(ORF_finder(frame, fwd_sequence).translate(table="Bacterial", to_stop=True))) + "aa")

            print("Reverse ORF " + str(frame + 1) + ": " + str(len(ORF_finder(frame, rev_comp_sequence))) + "bp, "
                  + str(len(ORF_finder(frame, rev_comp_sequence).translate(table="Bacterial", to_stop=True))) + "aa")

            ORF_list.append(str(ORF_finder(frame, fwd_sequence).translate(table="Bacterial", to_stop=True)))

            ORF_list.append(str(ORF_finder(frame, rev_comp_sequence).translate(table="Bacterial", to_stop=True)))

    # save longest amino acid sequence as a text file
    file = open('translatedORF.txt', 'w')
    file.write(max(ORF_list, key=len))
