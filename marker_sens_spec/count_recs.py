from Bio import SeqIO

def count_records(filename):
    count = 0
    for record in SeqIO.parse(filename, "fasta"):
        count += 1
    print(count)
count_records(r"C:\Users\abel\Documents\hydrogenases\marker_sens_spec\[FeFe]_Group_A1\0.8_2_test.fasta")
count_records(r"C:\Users\abel\Documents\hydrogenases\marker_sens_spec\[FeFe]_Group_A1_negative.fasta")

count_records(r"C:\Users\abel\Documents\hydrogenases\marker_sens_spec\[NiFe]_Group_1e\0.8_2_test.fasta")
count_records(r"C:\Users\abel\Documents\hydrogenases\marker_sens_spec\[NiFe]_Group_1e_negative.fasta")
count_records(r"C:\Users\abel\Documents\hydrogenases\marker_sens_spec\negative.fasta")