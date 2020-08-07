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
    
    menu = ["Intro", "About", "DNA Sequence", "DotPlot"]
    choice = st.sidebar.selectbox("Select Activity",menu)
    
    if choice == "Intro":
        st.subheader("Intro to BioInformatics")
        image = Image.open('dna.png')
        st.image(image,width=800)
        
        st.subheader("Bioinformatics")
        "Bioinformatika merupakan cabang ilmu dari biologi yang mengkombinasikan penggunaan komputerisasi dengan karakterisik molekuler biologi. Kombinasi ini disebabkan karena perkembangan teknologi informasi yang sangat pesat, sehingga memudahkan untuk dilakukannya penelitian, serta dapat memberikan informasi yang akurat berdasarkan pengelolaan data. Bioinformatika mempelajari interpretasi biologis data, dan evolusi berbagai bentuk kehidupan dari pendekatan komputasi."
        
        st.subheader("DNA")
        "DNA adalah singkatan dari asam deoksiribonukleat, yang merupakan molekul yang menyimpan informasi genetik utama dalam sel. Nukleotida terdiri dari tiga bagian: gugus fosfat, gula pentosa (gula ribosa), dan basa. Basisnya terdiri dari empat jenis: adenin (A), guanin (G), sitosin (C), dan timin (T). A dan G adalah purin dengan dua cincin menyatu. C dan T adalah pirimidin dengan satu cincin tunggal. Selain DNA, ada jenis nukleotida lain yang disebut RNA atau asam ribonukleat."
        
        st.subheader("Protein/Amino Acid")
        "Asam amino adalah senyawa organik yang memiliki gugus fungsional karboksil (-COOH) dan amina (biasanya -NH2). Dalam biokimia sering kali pengertiannya dipersempit: keduanya terikat pada satu atom karbon (C) yang sama (disebut atom C alfa atau α). Gugus karboksil memberikan sifat asam dan gugus amina memberikan sifat basa. Dalam bentuk larutan, asam amino bersifat amfoterik: cenderung menjadi asam pada larutan basa dan menjadi basa pada larutan asam."
        
        st.subheader("Biopython")
        "Biopython adalah seperangkat alat yang tersedia secara gratis untuk komputasi biologis yang ditulis dengan Python oleh tim pengembang internasional. Aplikasi ini dibuat dan dikembangkan dengan bahasa pemrograman python yang mana menggunakan library biopython untuk proses eksplorasi dan ekstraksi data. Dalam eksplorasi dan ekstraksi data, biopython dapat melakukan proses comparing sequences DNA dan transtlation Nucleotide DNA ke Amino Acid penyusun protein. Berikut ini merupakan tabel Codon transtlation dari Nucleotide ke Amino Acid."
        image = Image.open('protein.png')
        st.image(image,width=800)
        
    elif choice == "About":
        
        st.subheader("Sejarah Awal Coronavirus")
        "Coronavirus pertama kali ditemukan pada pertengahan tahun 1960 dengan jenis HCoV-229E. Virus ini bermutasi selama 56 tahun sampai pada tahun 2020 tercatat ada tujuh dari banyaknya jenis spesies virus corona yang menginfeksi manusia muali dari Alpha Coronavirus, Beta Coronavirus, SARS, dan juga MERS. Evolusi dari jenis spesies virus corona dapat terlihat pada gambar di bawah ini."
        image = Image.open('mutasi.PNG')
        st.image(image,width=800)
        
        st.subheader("Corona Virus Disease 2019")
        "Pandemi koronavirus 2019 (bahasa Inggris: coronavirus disease 2019, disingkat COVID-19) adalah penyakit menular yang disebabkan oleh SARS-CoV-2, salah satu jenis koronavirus. Penyakit ini mengakibatkan pandemi koronavirus 2019–2020.Penderita COVID-19 dapat mengalami demam, batuk kering, dan kesulitan bernapas.Sakit tenggorokan, pilek, atau bersin-bersin lebih jarang ditemukan.Pada penderita yang paling rentan, penyakit ini dapat berujung pada pneumonia dan kegagalan multiorgan.Infeksi menyebar dari satu orang ke orang lain melalui percikan (droplet) dari saluran pernapasan yang sering dihasilkan saat batuk atau bersin. Waktu dari paparan virus hingga timbulnya gejala klinis berkisar antara 1–14 hari dengan rata-rata 5 hari. Metode standar diagnosis adalah uji reaksi berantai polimerase transkripsi-balik (rRT-PCR) dari usap nasofaring atau sampel dahak dengan hasil dalam beberapa jam hingga 2 hari. Pemeriksaan antibodi dari sampel serum darah juga dapat digunakan dengan hasil dalam beberapa hari. Infeksi juga dapat didiagnosis dari kombinasi gejala, faktor risiko, dan pemindaian tomografi terkomputasi pada dada yang menunjukkan gejala pneumonia."
        
        st.subheader("Tujuan Pembuatan Aplikasi")
        "Tujuan pembuatan dan pengembangan aplikasi ini adalah agar dapat membantu para peneliti dalam menganalisis informasi yang ada pada data DNA dengan bentuk visualisasi diagram plot. Juga hasil informasi dari ekstraksi data DNA menghasilkan pola/patern transtlation protein dari sample DNA Coronavirus. Aplikasi ini dapat memberi hasil persentase similaritas sequencing alignment DNA, sehingga terlihat adanya mutasi genetik yang terjadi."
        
        st.subheader("Website ini dalam tahap pengembangan & digunakan untuk project penelitian.") 
        st.subheader("contact : hafiz065001600009.trisakti.ac.id")
        
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

            elif st.checkbox("Amino Acid Frequency"):
                st.write(aa_freq)

            elif st.checkbox("Plot Amino Acid Frequency"):
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
            dna_freq1 = Counter(dna_seq1)
            dna_freq2 = Counter(dna_seq2)
            p1 = dna_seq1.translate()
            aa_freq1 = Counter(str(p1))
            p2 = dna_seq2.translate()
            aa_freq2 = Counter(str(p2))
            details = st.radio("Details",("Description","Sequence","Nucleotide Frequency","Nucleotide Plot Frequency","Amino Acid Frequency","Amino Acid Plot Frequency"))
            if details == "Description":
                st.write(dna_record1.description)
                st.write("=====================")
                st.write(dna_record2.description)
                st.subheader("DNA Composition")
                gc_score1 = utils.gc_content(str(dna_seq1))
                at_score1 = utils.at_content(str(dna_seq1))
                st.json({"GC Content":gc_score1,"AT Content":at_score1})
                gc_score2 = utils.gc_content(str(dna_seq2))
                at_score2 = utils.at_content(str(dna_seq2))
                st.json({"GC Content":gc_score2,"AT Content":at_score2})
            elif details == "Sequence":
                st.write(dna_record1.seq)
                st.write("=========================================================================")
                st.write(dna_record2.seq)
            elif details == "Nucleotide Frequency":
                st.write(dna_freq1)
                st.write("=====================")
                st.write(dna_freq2)
            elif details == "Nucleotide Plot Frequency":
                barlist = plt.bar(dna_freq1.keys(),dna_freq1.values())
                st.pyplot()
                st.write("==========================================================================")
                barlist = plt.bar(dna_freq2.keys(),dna_freq2.values())
                st.pyplot()
            elif details == "Amino Acid Frequency":
                st.write(aa_freq1)
                st.write("=====================")
                st.write(aa_freq2)
            elif details == "Amino Acid Plot Frequency":
                plt.bar(aa_freq1.keys(),aa_freq1.values())
                st.pyplot()
                st.write("==========================================================================")
                plt.bar(aa_freq2.keys(),aa_freq2.values())
                st.pyplot()
            cus_limit = st.number_input("Select Max number of Nucleotide (Minimum 100)",100,40000,10000)
            if st.button("Dot Plot"):
                    st.write("Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
                    dotplotx(dna_seq1[0:cus_limit],dna_seq2[0:cus_limit])

                    st.pyplot()
            elif st.button("Similarity"):
                    st.write("Similarity of Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
                    r = pairwise2.align.globalxx(dna_seq1[0:cus_limit],dna_seq2[0:cus_limit],one_alignment_only=True,score_only=True)
                    r/len(dna_seq1[0:cus_limit]) * 100
        
        
if  __name__ == '__main__':
    main()
