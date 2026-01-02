import streamlit as st
from Bio import SeqIO
import io
import pandas as pd
import re
import requests
from bs4 import BeautifulSoup
import textwrap
import time

# ===============================
# Robust ProtParam function
# ===============================
def get_protparam(seq, max_retries=3, wait=0.5):
    """
    Fetch ProtParam data from Expasy and parse Molecular weight,
    Instability index, Instability label, Aliphatic index, GRAVY score.
    Returns a dict with None for missing fields.
    """
    url = 'https://web.expasy.org/cgi-bin/protparam/protparam'
    data = {'prot_id':'', 'mandatory':'', None:'', 'sequence':seq}
    
    retries = 0
    while retries < max_retries:
        try:
            res = requests.post(url, data=data, timeout=20)
            if res.status_code != 200:
                raise Exception(f"HTTP {res.status_code}")
            text = BeautifulSoup(res.text, 'html.parser').get_text()
            break
        except (requests.Timeout, requests.RequestException):
            retries += 1
            time.sleep(wait)
    else:
        raise RuntimeError(f"Failed to get ProtParam data after {max_retries} attempts")
    
    def parse_value(pattern, default=None):
        m = re.search(pattern, text)
        if m:
            try:
                return float(re.findall(r'-?\d+\.?\d*', m.group(0))[0])
            except:
                return default
        return default
    
    mol = parse_value(r'Molecular weight: -?[0-9]+\.[0-9]+')
    inst_score = parse_value(r'computed to be -?[0-9]+\.[0-9]+')
    inst_label = 'unstable' if 'unstable' in text.lower() else 'Stable'
    alip = parse_value(r'Aliphatic index: -?[0-9]+\.[0-9]+')
    gravy_score = parse_value(
    r'Grand average of hydropathicity \(GRAVY\):\s*-?\d+\.?\d*')
    
    return {
        'Molecular weight': mol,
        'Instability index': inst_score,
        'Instability': inst_label,
        'Aliphatic index': alip,
        'GRAVY score': gravy_score
    }

# ===============================
# Streamlit UI
# ===============================
st.set_page_config(page_title="Protein ProtParam Analyzer", layout="wide")
st.title("🧪 ProtParam Protein Parameter Analyzer")

st.markdown("""
Upload a **FASTA file** and compute **ProtParam parameters** (Molecular Weight, Instability, Aliphatic Index, GRAVY).  
Optionally, split sequences into **30-residue chunks** for detailed analysis.
""")

uploaded_file = st.file_uploader("Upload FASTA File", type=["fasta", "fa"])
chunk_option = st.checkbox("Split sequences into 30-residue chunks", value=False)

if uploaded_file:
    fasta_text = io.TextIOWrapper(uploaded_file, encoding="utf-8")
    records = list(SeqIO.parse(fasta_text, "fasta"))
    seqs = [str(r.seq).upper() for r in records]

    all_results = []

    with st.spinner("Computing ProtParam parameters..."):
        if chunk_option:
            # Split into 30-residue chunks
            for rec in records:
                protein_id = rec.id
                chunks = textwrap.wrap(str(rec.seq).upper(), 30)
                for idx, chunk in enumerate(chunks, 1):
                    res = get_protparam(chunk)
                    res['Protein_ID'] = protein_id
                    res['Chunk'] = idx
                    res['Sequence'] = chunk
                    all_results.append(res)
        else:
            for rec in records:
                protein_id = rec.id
                res = get_protparam(str(rec.seq).upper())
                res['Protein_ID'] = protein_id
                res['Sequence'] = str(rec.seq).upper()
                res['Chunk'] = 1
                all_results.append(res)

    df = pd.DataFrame(all_results)
    df = df[['Protein_ID', 'Chunk', 'Sequence', 'Molecular weight', 'Instability index', 'Instability', 'Aliphatic index', 'GRAVY score']]

    st.success(f"✅ Computed ProtParam for {len(df)} sequences/chunks")
    st.dataframe(df, use_container_width=True)

    st.download_button(
        label="⬇ Download Results as CSV",
        data=df.to_csv(index=False),
        file_name="protparam_results.csv",
        mime="text/csv"
    )
