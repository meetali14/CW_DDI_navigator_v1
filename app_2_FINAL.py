# Updated DDI Streamlit app
import os
import pandas as pd
import streamlit as st

# Constants
DATA_FILE_PATH = "tool.xlsx"
DATA_FILE_PATH2 = "ddiTable.xlsx"

# Load Data
def load_data(file_path):
    if not os.path.exists(file_path):
        st.error(f"The data file '{file_path}' does not exist.")
        return None
    try:
        return pd.read_excel(file_path)
    except Exception as e:
        st.error(f"An error occurred while loading the data: {e}")
        return None

# Search Drug Interaction
def search_drug_interaction(df, drug1, drug2):
    drug1, drug2 = drug1.strip().lower(), drug2.strip().lower()
    try:
        filtered_data = df[
            (df["Object compound"].str.lower() == drug1)
            & (df["Precipitant compound"].str.lower() == drug2)
        ]
        return filtered_data
    except KeyError as e:
        st.error(f"The expected column '{e}' is missing from the data file.")
        return pd.DataFrame()
    except Exception as e:
        st.error(f"An error occurred during the search: {e}")
        return pd.DataFrame()

# Drug-Drug Interaction Analysis
def ddi_analysis_for_drugs(file_path, drug1, drug2):
    try:
        df = pd.read_excel(file_path)
    except FileNotFoundError:
        return f"Error: The file '{file_path}' was not found.", []
    except Exception as e:
        return f"Error: {str(e)}", []

    try:
        drug1, drug2 = drug1.strip().lower(), drug2.strip().lower()
        # make safe lowercasing only if column exists
        if "Compound name" in df.columns:
            df["Compound name"] = df["Compound name"].astype(str).str.lower()

        drug1_genes = df[df["Compound name"] == drug1]["Interacting gene name"].unique() if "Compound name" in df.columns else []
        drug2_genes = df[df["Compound name"] == drug2]["Interacting gene name"].unique() if "Compound name" in df.columns else []

        if len(drug1_genes) == 0:
            return f"Warning: The drug '{drug1}' was not found in the dataset.", []
        if len(drug2_genes) == 0:
            return f"Warning: The drug '{drug2}' was not found in the dataset.", []
    except KeyError as e:
        return f"Error: Column '{e}' not found in the dataframe.", []
    except Exception as e:
        return f"Error: {str(e)}", []

    common_genes = set([str(g).upper() for g in drug1_genes if pd.notna(g)]) & set([str(g).upper() for g in drug2_genes if pd.notna(g)])
    result = []
    victim_perpetrator_results = []
    if common_genes:
        result = [
            f"Common genes found between '{drug1}' and '{drug2}': {', '.join(sorted(common_genes))}",
            "Drug-Drug Interaction (DDI) is possible based on these common genes.",
        ]

        # Build interactions list of tuples: (gene, substrate(if any), inducer(if any), inhibitor(if any))
        interactions = []
        for gene in common_genes:
            # default "-"
            sub_val = "-"
            ind_val = "-"
            inh_val = "-"
            try:
                # substrate check for drug1
                sub_rows = df[(df["Compound name"] == drug1) & (df["Interacting gene name"].astype(str).str.upper() == gene)]
                if "Compound Gene Relation" in df.columns:
                    try:
                        if any(sub_rows["Compound Gene Relation"].astype(str).str.lower().str.contains("substrate", na=False)):
                            sub_val = drug1
                    except Exception:
                        # fallback: ignore relation parsing errors
                        pass
                else:
                    # if no relation column, still mark substrate if the row exists
                    if not sub_rows.empty:
                        sub_val = drug1
            except Exception:
                pass

            try:
                # inducer check for drug2
                ind_rows = df[(df["Compound name"] == drug2) & (df["Interacting gene name"].astype(str).str.upper() == gene)]
                if "Compound Gene Relation" in df.columns:
                    try:
                        if any(ind_rows["Compound Gene Relation"].astype(str).str.lower().str.contains("inducer", na=False)):
                            ind_val = drug2
                    except Exception:
                        pass
                else:
                    if not ind_rows.empty:
                        # heuristic: if row exists, keep "-" unless relation indicates inducer
                        pass
            except Exception:
                pass

            try:
                # inhibitor check for drug2
                inh_rows = df[(df["Compound name"] == drug2) & (df["Interacting gene name"].astype(str).str.upper() == gene)]
                if "Compound Gene Relation" in df.columns:
                    try:
                        if any(inh_rows["Compound Gene Relation"].astype(str).str.lower().str.contains("inhibitor", na=False)):
                            inh_val = drug2
                    except Exception:
                        pass
            except Exception:
                pass

            interactions.append((gene, sub_val, ind_val, inh_val))

        if interactions:
            interaction_df = pd.DataFrame(interactions, columns=["Gene", "Substrate", "Inducer", "Inhibitor"])
            result.append(interaction_df)

            for _, row in interaction_df.iterrows():
                # victim-perpetrator detection (substrate==drug1 and (inducer==drug2 or inhibitor==drug2))
                if str(row["Substrate"]).strip().lower() == drug1.lower() and (
                    str(row["Inducer"]).strip().lower() == drug2.lower()
                    or str(row["Inhibitor"]).strip().lower() == drug2.lower()
                ):
                    victim_perpetrator_results.append(
                        f"'{row['Gene']}', '{drug1}' is likely to be a victim drug."
                    )
                else:
                    victim_perpetrator_results.append(
                        f"'{row['Gene']}', no sufficient information found for '{drug1}' as victim drug and '{drug2}' as perpetrator drug."
                    )

    else:
        result = [
            f"No common genes found between '{drug1}' and '{drug2}'. Cannot report possibility of Drug-Drug Interaction (DDI)."
        ]

    return result, victim_perpetrator_results

# ---------- Classification ----------
def classify_ddi_rule(common_genes_found: bool, victim_perpetrator_found: bool, clinical_evidence: bool):
    c = bool(common_genes_found)
    v = bool(victim_perpetrator_found)
    e = bool(clinical_evidence)

    if c and v and e:
        return "DDI (Yes)", 1, "Common genes + Victim-Perpetrator evidence + Clinical evidence"
    elif (not c) and (not v) and e:
        return "DDI (Yes)", 2, "No common gene & No victim-perpetrator, but Clinical evidence exists"
    elif c and (not v) and e:
        return "DDI (Yes)", 3, "Common genes found + Clinical evidence, but no victim-perpetrator identified"
    elif c and v and (not e):
        return "DDI (Yes)", 4, "Common genes + Victim-Perpetrator evidence but NO clinical evidence"
    elif (not c) and (not v) and (not e):
        return "DDI (No)", 5, "No common genes, no victim-perpetrator, and no clinical evidence"
    elif c and (not v) and (not e):
        return "DDI (No)", 6, "Common genes found but no victim-perpetrator and no clinical evidence"
    else:
        return "DDI (Undetermined)", 0, "Combination not covered by rules; manual review recommended"

# ---------- DDI strength mapping (per your reference table) ----------
def classify_ddi_strength(sub_imp, perp_imp):
    """
    Returns (ddi_strength_label, expected_auc_text, clinical_note)
    Accept expected values: 'Strong', 'Moderate', 'Weak' or '-'
    """
    s = str(sub_imp).capitalize() if sub_imp is not None else "-"
    p = str(perp_imp).capitalize() if perp_imp is not None else "-"

    table = {
        ("Strong", "Strong"): ("Strong", "Inhibitor: ≥5× ↑AUC (≥400%); Inducer: ≥5× ↓AUC (≥80%)", "Avoid co-administration or major dose adjustment; close monitoring; treat as critical for NTI drugs."),
        ("Strong", "Moderate"): ("Moderate", "Inhibitor: ~2–5× ↑AUC (100–400%); Inducer: ~2–5× ↓AUC (50–80%)", "Likely dose adjustment or close monitoring; consider alternatives; stronger action for NTI."),
        ("Strong", "Weak"): ("Weak", "Inhibitor: <2× ↑AUC (<100%); Inducer: <2× ↓AUC (<50%)", "Usually minor — monitor clinically; escalate if NTI or other risks."),
        ("Moderate", "Strong"): ("Strong / Moderate", "Inhibitor: often 2–5× ↑AUC; Inducer: often 2–5× ↓AUC", "Often clinically meaningful — treat as Moderate, escalate to Strong if other risk factors present."),
        ("Moderate", "Moderate"): ("Moderate", "Inhibitor: ~2× ↑AUC (20–200%); Inducer: ~2× ↓AUC (20–50%)", "Monitor, possible dose change; label warning."),
        ("Moderate", "Weak"): ("Weak", "Inhibitor: <2× ↑AUC; Inducer: <2× ↓AUC", "Usually no major action; monitor if multiple perpetrators or NTI substrate."),
        ("Weak", "Strong"): ("Moderate", "Inhibitor: commonly ≤2–5× ↑AUC; Inducer: ≤2–5× ↓AUC", "Sometimes clinically relevant — monitor; consider context."),
        ("Weak", "Moderate"): ("Weak / Moderate", "Inhibitor: <2× ↑AUC; Inducer: <2× ↓AUC", "Typically weak; monitor when other risks co-exist."),
        ("Weak", "Weak"): ("Weak", "Trivial AUC change", "No action for most drugs."),
    }
    return table.get((s, p), ("Undetermined", "No mapping available.", "Manual review recommended."))

# Streamlit UI
st.title("Drug-Drug Interaction Tool")
st.write("This tool allows you to search for interactions between two drugs. Select drugs to see detailed analysis.")

drug_data = load_data(DATA_FILE_PATH2)
interaction_result = pd.DataFrame()
ggi_text = []
victim_perpetrator_results = []

if drug_data is not None:
    # Build selectable drug list - handle case-insensitive and missing column safely
    try:
        custom_drug_list = sorted(drug_data["Compound name"].dropna().unique().tolist())
    except Exception:
        # fallback if column missing
        st.error("Column 'Compound name' not found in ddiTable.xlsx.")
        custom_drug_list = []

    drug1 = st.selectbox("Select the name of the first drug:", custom_drug_list)
    filtered_drug_list = [drug for drug in custom_drug_list if drug != drug1]
    drug2 = st.selectbox("Select the name of the second drug:", filtered_drug_list)

    if st.button("Analyze"):
        df = load_data(DATA_FILE_PATH)
        if df is not None:
            interaction_result = search_drug_interaction(df, drug1, drug2)

            SELECTED_COLUMNS = [
                "Object compound", "Precipitant compound", "Gene Name", "Gene type",
                "Impact on Object concentration", "AUC change (%)", "AUC fold change",
                "Evidence type", "Reference"
            ]
            SELECTED_COLUMNS_RENAME_MAP = {
                "Object compound": "Victim Drug",
                "Precipitant compound": "Precipitant Drug",
                "Gene Name": "Gene Name",
                "Gene type": "Gene type",
                "Impact on Object concentration": "Impact on Victim Drug concentration",
                "AUC change (%)": "AUC change (%)",
                "AUC fold change": "AUC fold change",
                "Evidence type": "Evidence type",
                "Reference": "Reference"
            }

            # DDI Analysis
            st.subheader("DDI Analysis:")
            ggi_text, victim_perpetrator_results = ddi_analysis_for_drugs(DATA_FILE_PATH2, drug1, drug2)

            # Show messages (strings in ggi_text)
            for line in ggi_text:
                if isinstance(line, str):
                    st.text(line)

            # Show gene-level table and capture it for impacts
            st.subheader("Potential Drug-Drug Interactions Based on their Gene:")
            interaction_df_for_impacts = None
            for line in ggi_text:
                if isinstance(line, pd.DataFrame):
                    st.table(line)
                    interaction_df_for_impacts = line.copy()

            # Victim-perpetrator text
            st.subheader("Victim-Perpetrator Analysis:")
            if victim_perpetrator_results:
                for vp_result in victim_perpetrator_results:
                    if isinstance(vp_result, str):
                        st.text(vp_result)

            # Clinical evidence table
            st.subheader("Evidence of Documented Drug Interactions:")
            if not interaction_result.empty:
                try:
                    interaction_result = interaction_result[SELECTED_COLUMNS].rename(columns=SELECTED_COLUMNS_RENAME_MAP)
                    st.table(interaction_result)
                except Exception:
                    st.warning("Evidence table missing expected columns; showing raw data.")
                    st.dataframe(interaction_result)
            else:
                st.write(f"No clinical evidence found between {drug1} and {drug2} to confirm the possibility of DDI.")

            # Classification (run after analysis)
            common_genes_found = any(isinstance(l, str) and "Common genes found" in l for l in ggi_text) or any(isinstance(l, pd.DataFrame) and not l.empty for l in ggi_text)
            victim_perpetrator_found = any(isinstance(vp, str) and "is likely to be a victim drug" in vp.lower() for vp in victim_perpetrator_results)
            clinical_evidence = not interaction_result.empty
            classification, rule_number, explanation = classify_ddi_rule(common_genes_found, victim_perpetrator_found, clinical_evidence)

            st.subheader("Automated Classification")
            col1, col2, col3 = st.columns([2, 1, 3])
            with col1:
                if classification == "DDI (Yes)":
                    st.success(classification)
                elif classification == "DDI (No)":
                    st.info(classification)
                else:
                    st.warning(classification)
            with col2:
                st.metric("Rule", f"{rule_number if rule_number else 'N/A'}")
            with col3:
                st.write(f"**Reason:** {explanation}")

            # -------- Gene Impact Sub-categorization --------
            # Will be used for per-gene DDI Strength Conclusions and Rule-based recommendation
            rule_based_label = "Undetermined"

            if interaction_df_for_impacts is None or interaction_df_for_impacts.empty:
                st.info("No gene-level interaction table available to extract impact strengths.")
            else:
                try:
                    ddi_table = pd.read_excel(DATA_FILE_PATH2)
                except Exception as e:
                    st.error(f"Cannot load ddi table to fetch Impact values: {e}")
                    ddi_table = pd.DataFrame()

                if ddi_table.empty:
                    st.warning("DDI source table is empty — cannot fetch impact strengths.")
                else:
                    # Normalize expected columns if present
                    # Make lowercase for compound and uppercase for gene matching
                    if "Compound name" in ddi_table.columns:
                        ddi_table["Compound name"] = ddi_table["Compound name"].astype(str).str.lower()
                    if "Interacting gene name" in ddi_table.columns:
                        ddi_table["Interacting gene name"] = ddi_table["Interacting gene name"].astype(str).str.upper()

                    # Prepare columns in the interaction DataFrame for impacts
                    interaction_df_for_impacts["Substrate Impact"] = "-"
                    interaction_df_for_impacts["Inducer Impact"] = "-"
                    interaction_df_for_impacts["Inhibitor Impact"] = "-"
                    # New: Reference column (collect references from ddi_table for that gene & compounds)
                    interaction_df_for_impacts["Reference"] = "-"

                    def normalize_impact(val):
                        if pd.isna(val):
                            return "-"
                        v = str(val).lower()
                        if "strong" in v:
                            return "Strong"
                        if "moderate" in v:
                            return "Moderate"
                        if "weak" in v:
                            return "Weak"
                        # numeric fallback
                        try:
                            num = float(val)
                            if abs(num) >= 2.0:
                                return "Strong"
                            elif abs(num) >= 1.0:
                                return "Moderate"
                            else:
                                return "Weak"
                        except Exception:
                            return str(val)

                    # Lookup impacts per gene-role and gather references
                    for idx, row in interaction_df_for_impacts.iterrows():
                        gene = str(row["Gene"]).upper()

                        # Substrate impact lookup (drug1)
                        try:
                            mask_sub = (ddi_table["Compound name"] == drug1.lower()) & (ddi_table["Interacting gene name"] == gene)
                            sub_hits = ddi_table[mask_sub]
                            if not sub_hits.empty and "Impact" in sub_hits.columns:
                                interaction_df_for_impacts.at[idx, "Substrate Impact"] = normalize_impact(sub_hits["Impact"].dropna().iloc[0])
                        except Exception:
                            pass

                        # Inducer impact lookup (drug2)
                        try:
                            mask_ind = (ddi_table["Compound name"] == drug2.lower()) & (ddi_table["Interacting gene name"] == gene)
                            ind_hits = ddi_table[mask_ind]
                            if not ind_hits.empty and "Compound Gene Relation" in ind_hits.columns and "Impact" in ind_hits.columns:
                                ind_hits = ind_hits[ind_hits["Compound Gene Relation"].astype(str).str.lower().str.contains("inducer", na=False)]
                                if not ind_hits.empty:
                                    interaction_df_for_impacts.at[idx, "Inducer Impact"] = normalize_impact(ind_hits["Impact"].dropna().iloc[0])
                        except Exception:
                            pass

                        # Inhibitor impact lookup (drug2)
                        try:
                            mask_inh = (ddi_table["Compound name"] == drug2.lower()) & (ddi_table["Interacting gene name"] == gene)
                            inh_hits = ddi_table[mask_inh]
                            if not inh_hits.empty and "Compound Gene Relation" in inh_hits.columns and "Impact" in inh_hits.columns:
                                inh_hits = inh_hits[inh_hits["Compound Gene Relation"].astype(str).str.lower().str.contains("inhibitor", na=False)]
                                if not inh_hits.empty:
                                    interaction_df_for_impacts.at[idx, "Inhibitor Impact"] = normalize_impact(inh_hits["Impact"].dropna().iloc[0])
                        except Exception:
                            pass

                        # --- Collect Reference(s) for this gene from ddi_table for both compounds ---
                        try:
                            refs = []
                            # references for drug1 + gene
                            mask_ref1 = (ddi_table["Compound name"] == drug1.lower()) & (ddi_table["Interacting gene name"] == gene)
                            hits_ref1 = ddi_table[mask_ref1]
                            if not hits_ref1.empty and "Reference" in hits_ref1.columns:
                                refs.extend([str(x).strip() for x in hits_ref1["Reference"].dropna().unique()])

                            # references for drug2 + gene
                            mask_ref2 = (ddi_table["Compound name"] == drug2.lower()) & (ddi_table["Interacting gene name"] == gene)
                            hits_ref2 = ddi_table[mask_ref2]
                            if not hits_ref2.empty and "Reference" in hits_ref2.columns:
                                refs.extend([str(x).strip() for x in hits_ref2["Reference"].dropna().unique()])

                            if refs:
                                refs_unique = sorted(set(refs))
                                interaction_df_for_impacts.at[idx, "Reference"] = "; ".join(refs_unique)
                        except Exception:
                            pass

                    # Present the Gene Impact Strengths table including Reference column
                    present_cols = ["Gene","Substrate","Substrate Impact","Inducer","Inducer Impact","Inhibitor","Inhibitor Impact","Reference"]
                    present_cols = [c for c in present_cols if c in interaction_df_for_impacts.columns]
                    st.subheader("Gene Impact Strengths")
                    st.table(interaction_df_for_impacts[present_cols])

                    # -------- DDI Strength Conclusion (per gene + rule-based overall label) --------
                    st.subheader("DDI Strength Conclusions (Rule-based)")
                    severity_rank = {"Strong": 4, "Strong / Moderate": 4, "Moderate": 3, "Weak / Moderate": 2.5, "Weak": 2, "Undetermined": 1}
                    # ranking for perp impact choice
                    rank_simple = {"Strong": 3, "Moderate": 2, "Weak": 1, "-": 0}

                    overall_scores = []  # collect severity ranks to compute rule-based overall

                    # Build a list of conclusion dicts to tabulate
                    conclusions = []

                    for idx, row in interaction_df_for_impacts.iterrows():
                        gene = row.get("Gene", "")
                        substrate_imp = row.get("Substrate Impact", "-")
                        # choose strongest perpetrator impact among Inducer and Inhibitor
                        ind_imp = row.get("Inducer Impact", "-")
                        inh_imp = row.get("Inhibitor Impact", "-")
                        # pick stronger impact by rank_simple
                        perp_imp = ind_imp if rank_simple.get(str(ind_imp).capitalize(),0) >= rank_simple.get(str(inh_imp).capitalize(),0) else inh_imp

                        ddi_strength, auc_change, notes = classify_ddi_strength(substrate_imp, perp_imp)

                        # map ddi_strength to a numeric score for overall (higher = worse)
                        if ddi_strength == "Strong":
                            score = 4
                        elif ddi_strength in ["Strong / Moderate"]:
                            score = 3.5
                        elif ddi_strength == "Moderate":
                            score = 3
                        elif ddi_strength in ["Weak / Moderate"]:
                            score = 2.5
                        elif ddi_strength == "Weak":
                            score = 2
                        elif ddi_strength == "Undetermined":
                            score = 1
                        else:
                            score = 1

                        overall_scores.append(score)

                        # Collect Reference
                        ref = row.get("Reference", "-") if "Reference" in row else "-"

                        # Append to conclusions list
                        conclusions.append({
                            "Gene": gene,
                            "Predicted DDI Strength": ddi_strength,
                            "Expected AUC Change": auc_change,
                            "Clinical Note": notes,
                            "Reference": ref
                        })

                    # Display conclusions as a table
                    if conclusions:
                        conclusions_df = pd.DataFrame(conclusions)
                        # Order columns for clarity
                        cols_order = ["Gene", "Predicted DDI Strength", "Expected AUC Change", "Clinical Note", "Reference"]
                        cols_order = [c for c in cols_order if c in conclusions_df.columns]
                        st.table(conclusions_df[cols_order])
                    else:
                        st.info("No per-gene conclusions could be computed.")

                    # Derive rule-based overall strength label (for RECOMMENDATION section)
                    if overall_scores:
                        max_score = max(overall_scores)
                        if max_score >= 3.75:
                            rule_based_label = "Strong"
                        elif max_score >= 3.0:
                            rule_based_label = "Moderate"
                        elif max_score >= 2.5:
                            rule_based_label = "Weak / Moderate"
                        elif max_score >= 2.0:
                            rule_based_label = "Weak"
                        else:
                            rule_based_label = "Undetermined"
                    else:
                        rule_based_label = "Undetermined"

            # -------- DDI Strength Conclusions (Clinical evidence-based) --------
            st.subheader("DDI Strength Conclusions (Clinical evidence-based)")

            clinical_conclusions = []
            clinical_label_overall = "Undetermined"
            clinical_scores = []

            if not interaction_result.empty:

                clinical_df = interaction_result.copy()

                # Ensure necessary columns exist (for safety)
                for col in [
                    "Victim Drug",
                    "Precipitant Drug",
                    "Gene Name",
                    "AUC change (%)",
                    "AUC fold change",
                    "Impact on Victim Drug concentration",
                ]:
                    if col not in clinical_df.columns:
                        clinical_df[col] = None

                def classify_fold(fold, impact):
                    # fold is AUC fold change; impact tells increase/decrease
                    try:
                        if fold is None or pd.isna(fold):
                            return "Undetermined", "No AUC fold change reported."
                        f = float(fold)
                    except Exception:
                        return "Undetermined", f"AUC fold change '{fold}' could not be interpreted numerically."

                    direction = str(impact).lower() if impact is not None else ""

                    # Magnitude-based classification (fold-wise)
                    if abs(f) >= 5:
                        strength = "Strong"
                        mag_text = "≥ 5-fold change → Strong DDI."
                    elif abs(f) >= 2:
                        strength = "Moderate"
                        mag_text = "≈ 2–5-fold change → Moderate DDI."
                    elif abs(f) >= 1:
                        strength = "Weak / Moderate"
                        mag_text = "≈ 1–2-fold change → Weak to Moderate DDI."
                    elif abs(f) > 0:
                        strength = "Weak"
                        mag_text = "< 1-fold change → Weak DDI."
                    else:
                        strength = "Undetermined"
                        mag_text = "No meaningful fold change reported."

                    # Directional explanation
                    if "increase" in direction:
                        dir_text = f"The precipitant appears to increase exposure of the victim drug by ~{f:.2f}-fold."
                    elif "decrease" in direction:
                        dir_text = f"The precipitant appears to decrease exposure of the victim drug by ~{f:.2f}-fold."
                    else:
                        dir_text = f"Observed AUC fold change of ~{f:.2f} (direction not clearly specified)."

                    explanation = f"{dir_text} {mag_text}"
                    return strength, explanation

                for _, row in clinical_df.iterrows():
                    fold = row.get("AUC fold change")
                    pct = row.get("AUC change (%)")
                    impact = row.get("Impact on Victim Drug concentration")

                    strength, explanation = classify_fold(fold, impact)

                    # score for clinical overall label
                    if strength == "Strong":
                        score = 4
                    elif strength == "Moderate":
                        score = 3
                    elif strength == "Weak / Moderate":
                        score = 2.5
                    elif strength == "Weak":
                        score = 2
                    else:
                        score = 1
                    clinical_scores.append(score)

                    clinical_conclusions.append({
                        "Victim Drug": row.get("Victim Drug"),
                        "Precipitant Drug": row.get("Precipitant Drug"),
                        "Gene Name": row.get("Gene Name"),
                        "AUC change (%)": pct,
                        "AUC fold change": fold,
                        "Clinical DDI Strength": strength,
                        "Interpretation": explanation
                    })

                if clinical_conclusions:
                    clinical_conclusions_df = pd.DataFrame(clinical_conclusions)
                    st.table(clinical_conclusions_df)

                if clinical_scores:
                    max_score_clin = max(clinical_scores)
                    if max_score_clin >= 3.75:
                        clinical_label_overall = "Strong"
                    elif max_score_clin >= 3.0:
                        clinical_label_overall = "Moderate"
                    elif max_score_clin >= 2.5:
                        clinical_label_overall = "Weak / Moderate"
                    elif max_score_clin >= 2.0:
                        clinical_label_overall = "Weak"
                    else:
                        clinical_label_overall = "Undetermined"
                else:
                    clinical_label_overall = "Undetermined"
            else:
                st.info("No AUC data available in the clinical evidence table to derive clinical DDI strength.")
                clinical_label_overall = "Undetermined"

            # -------- Unified RECOMMENDATION section --------
            st.subheader("RECOMMENDATION:")

            recommendation_texts = {
                "Strong": "High-risk DDI; avoid co-administration or apply major dose changes; close monitoring required.",
                "Moderate": "Clinically meaningful DDI; consider dose adjustment or close monitoring.",
                "Weak / Moderate": "Usually minor but monitor if other risk factors present.",
                "Weak": "Likely minimal clinical consequence.",
                "Undetermined": "Manual review recommended."
            }

            rule_text = recommendation_texts.get(rule_based_label or "Undetermined")
            clin_text = recommendation_texts.get(clinical_label_overall or "Undetermined")

            st.markdown(f"**(A) Rule-based:** {rule_based_label or 'Undetermined'} — {rule_text}")
            st.markdown(f"**(B) Clinical evidence-based:** {clinical_label_overall or 'Undetermined'} — {clin_text}")

    if st.button("Reset"):
        st.session_state.clear()
