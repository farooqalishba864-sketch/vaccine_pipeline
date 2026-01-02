import streamlit as st
import pandas as pd
from Bio import SeqIO

# ==============================
# BepiPred-style propensity scales
# ==============================

PARKER_SCALE = {
    'A': -0.5, 'C': -1.0, 'D': 3.0, 'E': 3.0, 'F': -2.5,
    'G': 0.0, 'H': -0.5, 'I': -1.8, 'K': 3.0, 'L': -1.8,
    'M': -1.3, 'N': 0.2, 'P': 0.0, 'Q': 0.2, 'R': 3.0,
    'S': 0.3, 'T': -0.4, 'V': -1.5, 'W': -3.4, 'Y': -2.3
}

FLEXIBILITY_SCALE = {
    'A': 1.0, 'C': 0.77, 'D': 1.04, 'E': 1.09, 'F': 0.77,
    'G': 1.19, 'H': 0.81, 'I': 0.88, 'K': 1.04, 'L': 0.81,
    'M': 0.94, 'N': 1.15, 'P': 1.32, 'Q': 1.10, 'R': 1.05,
    'S': 1.08, 'T': 1.09, 'V': 0.87, 'W': 0.81, 'Y': 0.76
}

EMINI_SCALE = {
    'A': 0.87, 'C': 1.52, 'D': 1.46, 'E': 1.46, 'F': 1.52,
    'G': 0.84, 'H': 1.22, 'I': 1.02, 'K': 1.11, 'L': 1.02,
    'M': 1.09, 'N': 1.29, 'P': 0.96, 'Q': 1.29, 'R': 1.11,
    'S': 1.20, 'T': 1.20, 'V': 0.96, 'W': 1.14, 'Y': 1.25
}

# ==============================
# Helper functions
# ==============================

def sliding_window(values, window):
    return [
        sum(values[i:i + window]) / window
        for i in range(len(values) - window + 1)
    ]


def normalize(values):
    min_v, max_v = min(values), max(values)
    if max_v == min_v:
        return [0.0] * len(values)
    return [(v - min_v) / (max_v - min_v) for v in values]


def bepipred_bcell_epitopes(fasta_file, window_size, threshold):
    results = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)

        hydro = [PARKER_SCALE.get(a, 0) for a in seq]
        flex = [FLEXIBILITY_SCALE.get(a, 0) for a in seq]
        access = [EMINI_SCALE.get(a, 0) for a in seq]

        hydro_w = normalize(sliding_window(hydro, window_size))
        flex_w = normalize(sliding_window(flex, window_size))
        access_w = normalize(sliding_window(access, window_size))

        combined = [(h + f + a) / 3 for h, f, a in zip(hydro_w, flex_w, access_w)]

        for i, score in enumerate(combined):
            if score >= threshold:
                results.append({
                    "Protein_ID": record.id,
                    "Epitope": seq[i:i + window_size],
                    "Start": i + 1,
                    "End": i + window_size,
                    "BepiPred_Score": round(score, 3)
                })

    df = pd.DataFrame(results)

    # 🔥 SORT HIGH → LOW SCORE
    if not df.empty:
        df = df.sort_values(
            by="BepiPred_Score",
            ascending=False
        ).reset_index(drop=True)

    return df

# ==============================
# Streamlit UI
# ==============================

st.set_page_config(page_title="B-cell Epitope Predictor", layout="wide")

st.title("🧬 B-cell Epitope Prediction (BepiPred-style)")
st.markdown(
    "Predict **linear B-cell epitopes** and rank them by **BepiPred score**."
)

uploaded_file = st.file_uploader(
    "Upload Protein FASTA file",
    type=["fasta", "fa", "faa"]
)

col1, col2 = st.columns(2)

with col1:
    window_size = st.slider("Sliding Window Size", 5, 15, 9)

with col2:
    threshold = st.slider("Prediction Threshold", 0.3, 0.9, 0.55, 0.01)

if uploaded_file:
    with st.spinner("Predicting B-cell epitopes..."):
        with open("temp.fasta", "wb") as f:
            f.write(uploaded_file.read())

        df = bepipred_bcell_epitopes(
            fasta_file="temp.fasta",
            window_size=window_size,
            threshold=threshold
        )

    st.success(f"✅ {len(df)} B-cell epitopes identified (sorted by score)")

    if not df.empty:
        st.dataframe(df, use_container_width=True)

        st.download_button(
            label="⬇ Download CSV",
            data=df.to_csv(index=False),
            file_name="bepipred_bcell_epitopes_sorted.csv",
            mime="text/csv"
        )
    else:
        st.warning("No epitopes found at the selected threshold.")
