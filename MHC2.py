import streamlit as st
from Bio import SeqIO
import io
import os
import subprocess
import pandas as pd

# --------- Streamlit UI ---------
st.title("MHC-II Epitope Prediction with MixMHC2pred")
st.write("""
Upload a FASTA file and predict 12-aa peptides against 10 DR, 10 DQ, and 10 DP alleles.
""")

fasta_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa","faa"])

if fasta_file is not None:
    st.success("FASTA file uploaded!")

    # --------- Convert uploaded file to text mode ---------
    fasta_text = io.TextIOWrapper(fasta_file, encoding='utf-8')

    # --------- Function: make 12-aa peptides ---------
    def make_12mer_peptides(fasta_handle):
        peptides = []
        for record in SeqIO.parse(fasta_handle, "fasta"):
            seq = str(record.seq)
            for i in range(len(seq) - 11):
                peptide = seq[i:i+12]
                peptides.append(peptide)
        return peptides

    peptides = make_12mer_peptides(fasta_text)

    # --------- Create MixMHC2pred input file ---------
    input_txt = "mixmhc_input.txt"
    with open(input_txt, "w") as f:
        for pep in peptides:
            f.write(f"{pep}\n")  # Using --no_context
    st.write(f"Created MixMHC2pred input file: `{input_txt}` with {len(peptides)} peptides")

    # --------- Allele lists ---------
    DR_alleles = [
        "DRB1_01_01","DRB1_03_01","DRB1_04_01","DRB1_07_01","DRB1_09_01",
        "DRB1_11_01","DRB1_13_01","DRB1_15_01","DRB3_01_01","DRB5_01_01"
    ]
    DQ_alleles = [
        "DQA1_01_01__DQB1_05_01","DQA1_01_02__DQB1_05_02","DQA1_03_01__DQB1_03_02",
        "DQA1_04_01__DQB1_04_02","DQA1_05_01__DQB1_02_01","DQA1_05_01__DQB1_03_01",
        "DQA1_02_01__DQB1_02_02","DQA1_01_01__DQB1_06_02","DQA1_03_03__DQB1_03_01",
        "DQA1_02_01__DQB1_03_03"
    ]
    DP_alleles = [
        "DPA1_01_03__DPB1_01_01","DPA1_02_01__DPB1_01_01","DPA1_01_03__DPB1_02_01",
        "DPA1_02_01__DPB1_02_01","DPA1_01_03__DPB1_04_01","DPA1_02_01__DPB1_04_01",
        "DPA1_01_03__DPB1_05_01","DPA1_02_01__DPB1_05_01","DPA1_01_03__DPB1_09_01",
        "DPA1_02_01__DPB1_09_01"
    ]

    # --------- Function to run MixMHC2pred ---------
    def run_mixmhcpred(input_file, output_file, alleles):
        allele_str = " ".join(alleles)
        cmd = f"MixMHC2pred.exe -i {input_file} -o {output_file} -a {allele_str} --no_context"
        st.write(f"Running: `{cmd}`")
        subprocess.run(cmd, shell=True)
        st.success(f"Results saved to `{output_file}`")
        return output_file

    # --------- Function to read MixMHC2pred output safely ---------
    def read_mixmhc_file(file):
        with open(file, 'r') as f:
            lines = f.readlines()
        header_line_index = 0
        for i, line in enumerate(lines):
            if line.startswith("Peptide"):
                header_line_index = i
                break
        df = pd.read_csv(file, sep="\t", header=header_line_index)
        if 'Context' in df.columns:
            df = df.drop(columns=['Context'])
        return df

    # --------- Run MixMHC2pred on button click ---------
    if st.button("Run MixMHC2pred"):
        # Run DR, DQ, DP
        out_dr = "results_DR.txt"
        out_dq = "results_DQ.txt"
        out_dp = "results_DP.txt"
        run_mixmhcpred(input_txt, out_dr, DR_alleles)
        run_mixmhcpred(input_txt, out_dq, DQ_alleles)
        run_mixmhcpred(input_txt, out_dp, DP_alleles)

        st.success("All predictions completed!")

        # Read all results
        df_dr = read_mixmhc_file(out_dr)
        df_dq = read_mixmhc_file(out_dq)
        df_dp = read_mixmhc_file(out_dp)

        # Merge all results into one table
        df_all = pd.concat([df_dr, df_dq, df_dp], ignore_index=True)
        st.write("Merged results from DR/DQ/DP alleles:")
        st.dataframe(df_all.head(20))

        # --------- Download button ---------
        csv_all = df_all.to_csv(index=False)
        st.download_button(
            label="Download all results as CSV",
            data=csv_all,
            file_name="all_MHCII_results.csv",
            mime="text/csv"
        )

        # --------- Top 15 epitopes button ---------
        if st.button("Show Top 15 Epitopes"):
            # Assuming %Rank columns exist for all alleles
            rank_cols = [c for c in df_all.columns if "%Rank" in c]
            df_all['min_rank'] = df_all[rank_cols].min(axis=1)
            df_top15 = df_all.nsmallest(15, 'min_rank')
            st.write("Top 15 predicted epitopes:")
            st.dataframe(df_top15.drop(columns=['min_rank']))
