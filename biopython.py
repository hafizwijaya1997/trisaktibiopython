import streamlit as st
import matplotlib.pyplot as plt 
import matplotlib
matplotlib.use("Agg")
from Bio.Seq import Seq 
from Bio import SeqIO
from collections import Counter
import neatbio.sequtils as utils
import numpy as np 
from PIL import Image 
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def delta(x,y):
    return 0 if x == y else 1


def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))


def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]

def plotMatrix(M,t, seq1, seq2, nonblank = chr(0x25A0), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1,M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)


def dotplot(seq1,seq2,k = 1,t = 1):
    M = makeMatrix(seq1,seq2,k)
    plotMatrix(M, t, seq1,seq2) #experiment with character choice


# Convert to Fxn
def dotplotx(seq1,seq2):
    plt.imshow(np.array(makeMatrix(seq1,seq2,1)))
    # on x-axis list all sequences of seq 2
    xt=plt.xticks(np.arange(len(list(seq2))),list(seq2))
    # on y-axis list all sequences of seq 1
    yt=plt.yticks(np.arange(len(list(seq1))),list(seq1))
    plt.show()


def gc_content(seq):
    result = float(str(seq).count('G') + str(seq).count('C'))/len(seq) * 100
    return result

def at_content(seq):
    result = float(str(seq).count('A') + str(seq).count('T'))/len(seq) * 100
    return result


def main():
    image = Image.open('b2.png')
    st.image(image,width=200)
    st.title("Trisakti Bioinformatics Application")
    st.title("Powered by Python")
    
    menu = ["Intro", "DNA Sequence", "DotPlot", "About"]
    choice = st.sidebar.selectbox("Select Activity",menu)
    
    if choice == "Intro":
        st.subheader("Intro to BioInformatics")
        image = Image.open('dna.png')
        st.image(image,width=800)
        
        st.subheader("Bioinformatics")
        "Bioinformatika merupakan cabang ilmu dari biologi yang mengkombinasikan penggunaan komputerisasi dengan karakterisik molekuler biologi. Kombinasi ini disebabkan karena perkembangan teknologi informasi yang sangat pesat, sehingga memudahkan untuk dilakukannuya penelitian, serta dapat memberikan informasi yang akurat berdasarkan pengelolaan data.Bioinformatika mempelajari interpretasi biologis data, terutama data asam nukleat dan asam amino, dan mempelajari aturan molekuler dan sistem yang mengatur atau memengaruhi struktur/ sequen, fungsi, dan evolusi berbagai bentuk kehidupan dari pendekatan komputasi. Kata komputasi tidak hanya berarti dengan komputer tetapi mengacu pada analisis data dengan matematika, metode statistik, dan algoritmik, yang sebagian besar perlu diimplementasikan program komputer."
        
        st.subheader("DNA")
        "DNA adalah singkatan dari asam deoksiribonukleat, yang merupakan molekul yang menyimpan informasi genetik utama dalam sel. Nukleotida terdiri dari tiga bagian: gugus fosfat, gula pentosa (gula ribosa), dan basa. Basisnya terdiri dari empat jenis: adenin (A), guanin (G), sitosin (C), dan timin (T). A dan G adalah purin dengan dua cincin menyatu. C dan T adalah pirimidin dengan satu cincin tunggal. Selain DNA, ada jenis nukleotida lain yang disebut RNA atau asam ribonukleat. Untuk RNA juga terdiri dari empat jenis yang sama dengan DNA kecuali timin (T) yang digantikan oleh uracil (U) di RNA. DNA biasanya terdiri dari dua untai yang berjalan berlawanan arah. Tulang belakang dari setiap untai adalah serangkaian kelompok pentosa dan fosfat. Ikatan hydrogen antara purin dan pirimidin memegang dua untai DNA bersamaan, membentuk heliks ganda yang terkenal." 
        
        "Dalam ikatan hidrogen, basa A selalu berpasangan dengan basa T didudukan dengan yang lainnya dan G selalu dengan C. Mekanisme ini disebut pemasangan pasangan. RNA biasanya satu untai. Ketika untai RNA berpasangan dengan untai DNA, maka aturan pasangan-pasangan menjadi A-U, T-A, G-C, dan C-G. Gula ribosa disebut gula pentosa karena mengandung lima karbon, masing-masing bernomor 10 –50. Definisi arah untai DNA atau RNA juga didasarkan pada penomoran ini, sehingga kedua ujung untai DNA atau RNA disebut ujung 50 dan ujung 30. Serangkaian pangkalan di sepanjang untai disebut urutan DNA atau RNA dan dapat dilihat sebagai string karakter yang dikomposisikan alfabet A, C, G, dan T (U untuk RNA). Pada heliks DNA ganda, kedua untai berjalan berlawanan. adalah contoh dari segmen urutan DNA untai ganda."
        
    elif choice == "DNA Sequence":
        st.subheader("DNA Sequence Analysis")
        
        seq_file = st.file_uploader("Upload FASTA File", type=["fasta","fa","txt"])
        
        if seq_file is not None:
            dna_record = SeqIO.read(seq_file,"fasta")
            #st.write(dna_record)
            dna_seq = dna_record.seq
            
            details = st.radio("Details",("Description","Sequence"))
            if details == "Description":
                st.write(dna_record.description)
            elif details == "Sequence":
                st.write(dna_record.seq)
                
            # Frekuensi Nucleotide
            st.subheader("Nucleotide Frequency")
            dna_freq = Counter(dna_seq)
            st.write(dna_freq)
            adenine_color = st.beta_color_picker("Adenine Color")
            thymine_color = st.beta_color_picker("Thymine Color")
            guanine_color = st.beta_color_picker("Guanine Color")
            cytosil_color = st.beta_color_picker("Cytosil Color")
            
            if st.button("Plot Frequency"):
                barlist = plt.bar(dna_freq.keys(),dna_freq.values())
                barlist[2].set_color(adenine_color)
                barlist[3].set_color(thymine_color)
                barlist[1].set_color(guanine_color)
                barlist[0].set_color(cytosil_color)
                
                
                st.pyplot()
                
            st.subheader("DNA Composition")
            gc_score = utils.gc_content(str(dna_seq))
            at_score = utils.at_content(str(dna_seq))
            st.json({"GC Content":gc_score,"AT Content":at_score})
            
            # Count Nucleotide
            nt_count = st.text_input("Enter Nucleotide","Type Nucleotide Alphabet")
            st.write("Number of {} Nucleotide is ::{}".format((nt_count),str(dna_seq).count(nt_count)))
            
            # Protein Synthesis
            st.subheader("Protein Synthesis")
            p1 = dna_seq.translate()
            aa_freq = Counter(str(p1))

            if st.checkbox("Transcription"):
                st.write(dna_seq.transcribe())

            elif st.checkbox("Translation"):
                st.write(dna_seq.translate())

            elif st.checkbox("Complement"):
                st.write(dna_seq.complement())

            elif st.checkbox("AA Frequency"):
                st.write(aa_freq)

            elif st.checkbox("Plot AA Frequency"):
                aa_color = st.beta_color_picker("Pick An Amino Acid Color")
                #barlist = plt.bar(aa_freq.keys(),aa_freq.values())
                #barlist[2].set_color(aa_color)
                plt.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
                st.pyplot()

            elif st.checkbox("Full Amino Acid Name"):
                aa_name = str(p1).replace("*","")
                aa3 = utils.convert_1to3(aa_name)
                st.write(aa_name)
                st.write("=========================")
                st.write(aa3)
                
                st.write("=========================")
                st.write(utils.get_acid_name(aa3))
                
        
    elif choice == "DotPlot":
        st.subheader("Generate Dot Plot For Two Sequences")
        seq_file1 = st.file_uploader("Upload 1st FASTA File",type=["fasta","fa"])
        seq_file2 = st.file_uploader("Upload 2nd FASTA File",type=["fasta","fa"])

        if seq_file1 and seq_file2 is not None:
            dna_record1 = SeqIO.read(seq_file1,"fasta")
            dna_record2 = SeqIO.read(seq_file2,"fasta")
            # st.write(dna_record)
            dna_seq1 = dna_record1.seq
            dna_seq2 = dna_record2.seq
            details = st.radio("Details",("Description","Sequence"))
            if details == "Description":
                st.write(dna_record1.description)
                st.write("=====================")
                st.write(dna_record2.description)
            elif details == "Sequence":
                st.write(dna_record1.seq)
                st.write("=====================")
                st.write(dna_record2.seq)
            cus_limit = st.number_input("Select Max number of Nucleotide (Minimum 100)",100,40000,10000)
            if st.button("Dot Plot"):
                    st.write("Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
                    dotplotx(dna_seq1[0:cus_limit],dna_seq2[0:cus_limit])

                    st.pyplot()
            elif st.button("Similarity"):
                    st.write("Similarity of Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
                    r = pairwise2.align.globalxx(dna_seq1[0:cus_limit],dna_seq2[0:cus_limit],one_alignment_only=True,score_only=True)
                    r/len(dna_seq1[0:cus_limit]) * 100
        
    elif choice == "About":
        st.subheader("Website ini dalam tahap pengembangan & digunakan untuk project penelitian.") 
        st.subheader("contact : hafiz065001600009.trisakti.ac.id")
        
        
if  __name__ == '__main__':
    main()