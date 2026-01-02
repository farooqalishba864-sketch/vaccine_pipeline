import streamlit as st
from mhcflurry import Class1PresentationPredictor
from Bio import SeqIO
from io import StringIO
import pandas as pd

st.title("MHC-I Epitope Prediction (All Alleles, User Threshold, Proper Allele)")

uploaded = st.file_uploader("Upload FASTA file", type=["fasta", "fa"])
threshold = st.number_input(
    "Set IC50 threshold (nM) to define strong binders", 
    min_value=1, max_value=10000, value=500, step=50
)

ALLELES = [
    "HLA-A*01:01","HLA-A*02:01","HLA-A*02:03","HLA-A*02:06",
    "HLA-A*03:01","HLA-A*11:01","HLA-A*23:01","HLA-A*24:02",
    "HLA-A*26:01","HLA-A*30:01","HLA-A*30:02","HLA-A*31:01",
    "HLA-A*32:01","HLA-A*33:01","HLA-A*68:01","HLA-A*68:02",
    "HLA-B*07:02","HLA-B*08:01","HLA-B*15:01","HLA-B*35:01",
    "HLA-B*40:01","HLA-B*44:02","HLA-B*44:03","HLA-B*51:01",
    "HLA-B*53:01","HLA-B*57:01","HLA-B*58:01"
]

PEPTIDE_LENGTH = 9

def chunk_list(lst, size=6):
    for i in range(0, len(lst), size):
        yield lst[i:i+size]

if uploaded:
    st.info("⏳ Prediction running… this may take some time for multiple proteins and alleles.")
    predictor = Class1PresentationPredictor.load()
    fasta = StringIO(uploaded.getvalue().decode())
    results = []

    for record in SeqIO.parse(fasta, "fasta"):
        sequence = str(record.seq)
        peptides = [sequence[i:i+PEPTIDE_LENGTH] for i in range(len(sequence)-PEPTIDE_LENGTH+1)]

        for allele_chunk in chunk_list(ALLELES, 6):
            df = predictor.predict(peptides=peptides, alleles=allele_chunk)
            df["Protein"] = record.id
            df["Peptide_Length"] = PEPTIDE_LENGTH

            # 🔹 Ensure allele column exists
            if "allele" not in df.columns and "mhc_allele" in df.columns:
                df.rename(columns={"mhc_allele":"allele"}, inplace=True)
            if "allele" not in df.columns:
                df["allele"] = allele_chunk[0]  # fallback, at least first allele

            df["Strong_Binder"] = df["affinity"] <= threshold

            results.append(df)

    if results:
        final = pd.concat(results, ignore_index=True)
        final.sort_values(by="Strong_Binder", ascending=False, inplace=True)

        st.subheader("Prediction Results (Preview)")
        st.dataframe(final[["Protein","peptide","allele","affinity","Strong_Binder"]].head(50))

        st.download_button(
            "Download Full CSV",
            final.to_csv(index=False),
            "mhc1_full_predictions.csv"
        )
    else:
        st.warning("No peptides found in uploaded FASTA.")
