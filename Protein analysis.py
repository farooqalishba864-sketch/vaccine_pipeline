import streamlit as st
import pandas as pd
import numpy as np
from Bio import SeqIO
from io import StringIO

# ================= CONFIG ================= #
st.set_page_config(page_title="Vaccine Candidate Screening", layout="wide")
VAXIJEN_THRESHOLD = 0.4

# ================= PHYSICOCHEMICAL PROPERTIES ================= #
AA_PROP = {
    'A': 0.62, 'C': 0.29, 'D': -0.90, 'E': -0.74, 'F': 1.19,
    'G': 0.48, 'H': -0.40, 'I': 1.38, 'K': -1.50, 'L': 1.06,
    'M': 0.64, 'N': -0.78, 'P': 0.12, 'Q': -0.85, 'R': -2.53,
    'S': -0.18, 'T': -0.05, 'V': 1.08, 'W': 0.81, 'Y': 0.26
}

# ================= FUNCTIONS ================= #
def auto_cross_covariance(signal, lag=5):
    mean = np.mean(signal)
    return [
        np.mean([(signal[i] - mean) * (signal[i + l] - mean)
                 for i in range(len(signal) - l)])
        for l in range(1, lag + 1)
    ]

def vaxijen_score(seq):
    weights = {'A':0.5,'C':0.3,'D':0.2,'E':0.2,'F':0.8,'G':0.4,
               'I':0.6,'L':0.7,'M':0.5,'V':0.6,'W':0.9,'Y':0.7}
    score = sum(weights.get(a, 0.4) for a in seq) / len(seq)
    return round(score, 3), score >= VAXIJEN_THRESHOLD

def allertop(seq):
    if seq.count('C') / len(seq) > 0.12:
        return "Allergen"
    for motif in ['KYSL', 'PFYI', 'NGG']:
        if motif in seq:
            return "Allergen"
    return "Non-Allergen"

def toxinpred(seq):
    kr_fraction = (seq.count('K') + seq.count('R')) / len(seq)
    if kr_fraction > 0.25 and len(seq) <= 30:
        return "Toxic"
    return "Non-Toxic"

def immunogenicity_like(sequence):
    signal = []
    for aa in sequence:
        if aa in AA_PROP:
            val = AA_PROP[aa]
            if aa in "WFY": val *= 1.3
            if aa in "KR": val *= 0.8
            signal.append(val)
    if len(signal) < 9:
        return None, "Insufficient length"
    acc = auto_cross_covariance(signal)
    score = sum(w * f for w, f in zip([0.4,0.3,0.2,0.1,0.05], acc))
    return round(score, 3), "Immunogenic" if score > 0 else "Non-Immunogenic"

# ================= UI ================= #
st.title("🧬 Physicochemical Vaccine Candidate Screening Tool")
st.markdown("""
**Offline | FASTA-based | Rule-driven Reverse Vaccinology Pipeline**
""")

uploaded_file = st.file_uploader("Upload Protein FASTA File", type=["fasta", "fa", "faa", "txt"])

if uploaded_file:
    fasta_text = uploaded_file.read().decode("utf-8")
    records = list(SeqIO.parse(StringIO(fasta_text), "fasta"))
    st.success(f"✅ {len(records)} proteins loaded")

    results = []

    for r in records:
        seq = str(r.seq)
        vax_score, antigen_flag = vaxijen_score(seq)
        allergen = allertop(seq)
        toxic = toxinpred(seq)
        immuno_score, immuno_label = immunogenicity_like(seq)

        results.append({
            "Protein_ID": r.id,
            "Sequence": seq,
            "Length": len(seq),
            "VaxiJen_Score": vax_score,
            "Antigenicity": "Antigen" if antigen_flag else "Non-Antigen",
            "Allergenicity": allergen,
            "Toxicity": toxic,
            "Immunogenicity_Score": immuno_score,
            "Immunogenicity": immuno_label
        })

    df = pd.DataFrame(results)

    # ================= FULL RESULTS ================= #
    st.subheader("🔬 Full Protein Analysis")
    st.dataframe(df.drop(columns=["Sequence"]), use_container_width=True)

    st.download_button(
        "⬇ Download Full Results (CSV)",
        df.drop(columns=["Sequence"]).to_csv(index=False),
        "protein_analysis_full.csv",
        "text/csv"
    )

    # ================= IDEAL FILTER ================= #
    ideal_df = df[
        (df["Antigenicity"] == "Antigen") &
        (df["Allergenicity"] == "Non-Allergen") &
        (df["Toxicity"] == "Non-Toxic") &
        (df["Immunogenicity"] == "Immunogenic")
    ]

    # ================= STATISTICS ================= #
    st.subheader("📊 Protein Analysis Statistics")

    stats = {
        "Total Proteins": len(df),
        "Antigenic": (df["Antigenicity"] == "Antigen").sum(),
        "Non-Allergen": (df["Allergenicity"] == "Non-Allergen").sum(),
        "Non-Toxic": (df["Toxicity"] == "Non-Toxic").sum(),
        "Immunogenic": (df["Immunogenicity"] == "Immunogenic").sum(),
        "Ideal Vaccine Candidates": len(ideal_df)
    }

    st.table(pd.DataFrame(stats.items(), columns=["Category", "Count"]))

    # ================= FILTERED PREVIEW ================= #
    st.subheader("🔍 Ideal Vaccine Candidate Proteins")

    if not ideal_df.empty:
        st.dataframe(
            ideal_df.drop(columns=["Sequence"]),
            use_container_width=True
        )

        # CSV download
        st.download_button(
            "⬇ Download Filtered Candidates (CSV)",
            ideal_df.drop(columns=["Sequence"]).to_csv(index=False),
            "vaccine_candidates_filtered.csv",
            "text/csv"
        )

        # FASTA generation
        fasta_out = ""
        for _, row in ideal_df.iterrows():
            fasta_out += f">{row['Protein_ID']}\n{row['Sequence']}\n"

        st.download_button(
            "⬇ Download Ideal Proteins (FASTA)",
            fasta_out,
            "ideal_vaccine_candidates.fasta",
            "text/plain"
        )
    else:
        st.warning("No proteins satisfied all vaccine candidate criteria.")
