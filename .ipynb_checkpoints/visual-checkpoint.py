import json, sys, re
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt

# -------- helpers --------
def make_res_label(item):
    """Return a readable label like 'GLU40_A' or 'GLU40' from a PLIP residue dict or string."""
    if not isinstance(item, dict):
        s = str(item)
        m = re.match(r'([A-Z]{2,4})\s*(-?\d+)\s*([A-Za-z]?)', s)
        if m:
            return f"{m.group(1)}{m.group(2)}_{m.group(3)}" if m.group(3) else f"{m.group(1)}{m.group(2)}"
        return s
    name = item.get("resname") or item.get("residue_name") or item.get("residue") or item.get("name")
    num = item.get("resnr") or item.get("residue_number") or item.get("residue_id") or item.get("resSeq") or item.get("number")
    chain = item.get("chain") or item.get("chain_id") or item.get("chainId") or ""
    if name is None and isinstance(item.get("label"), str):
        # sometimes PLIP puts full label in a text field
        txt = item.get("label")
        m = re.match(r'([A-Z]{2,4})\s*(-?\d+)\s*([A-Za-z]?)', txt)
        if m:
            name, num, chain = m.group(1), m.group(2), m.group(3) or chain
    name = name or "UNK"
    num = num or "?"
    try:
        num = int(re.search(r"-?\d+", str(num)).group())
    except:
        num = str(num)
    return f"{name}{num}_{chain}" if chain else f"{name}{num}"

def normalize_json(plip_json):
    """
    Map many possible PLIP JSON keys to canonical bins:
      'H-bonds', 'Hydrophobic', 'Ionic', 'Water bridges'
    Each value will be a list of residue-like dicts or strings.
    """
    mapping = {
        'hydrogen': 'H-bonds',
        'hbond': 'H-bonds',
        'hb': 'H-bonds',
        'hydrophobic': 'Hydrophobic',
        'hydro': 'Hydrophobic',
        'hyd': 'Hydrophobic',
        'ionic': 'Ionic',
        'salt': 'Ionic',
        'water': 'Water bridges',
        'bridge': 'Water bridges',
        'waters': 'Water bridges'
    }
    out = {v: [] for v in set(mapping.values())}

    # If top-level canonical keys already present, use them
    lower_keys = {k.lower(): k for k in plip_json.keys()}
    for canonical in ["hydrogen_bonds","hydrophobic_contacts","ionic_interactions","water_bridges"]:
        if canonical in lower_keys:
            out_map = {
                'hydrogen_bonds': 'H-bonds',
                'hydrophobic_contacts': 'Hydrophobic',
                'ionic_interactions': 'Ionic',
                'water_bridges': 'Water bridges'
            }[canonical]
            out[out_map].extend(plip_json[lower_keys[canonical]])

    # Otherwise iterate keys heuristically
    for k, val in plip_json.items():
        if not isinstance(val, list):
            continue
        lk = k.lower()
        chosen = None
        for kw, mapped in mapping.items():
            if kw in lk:
                chosen = mapped
                break
        if chosen:
            out[chosen].extend(val)
            continue
        # try inspecting first list item for 'type' hints
        if val and isinstance(val[0], dict):
            t = (val[0].get('type') or val[0].get('interaction') or "")
            for kw, mapped in mapping.items():
                if isinstance(t, str) and kw in t.lower():
                    chosen = mapped
                    break
        if chosen:
            out[chosen].extend(val)
        else:
            # fallback: if we can't classify, put into Hydrophobic (most common)
            out['Hydrophobic'].extend(val)
    return out

def build_counts(normalized):
    counts = defaultdict(lambda: {'H-bonds':0, 'Hydrophobic':0, 'Ionic':0, 'Water bridges':0})
    total = 0
    for itype, items in normalized.items():
        for it in items:
            label = make_res_label(it)
            counts[label][itype] += 1
            total += 1
    if not counts:
        return pd.DataFrame(columns=['H-bonds','Hydrophobic','Ionic','Water bridges']), total
    df = pd.DataFrame.from_dict(counts, orient='index').fillna(0).astype(int)
    # sort by residue number when possible
    def sort_key(lbl):
        m = re.search(r"(-?\d+)", lbl)
        return int(m.group(1)) if m else 10**9
    df = df.reindex(sorted(df.index, key=sort_key))
    return df, total

# -------- plotting --------
def plot_stacked_fraction(df, total, out_png, title="Protein-Ligand Contacts (PLIP)"):
    # compute fractions (per total interactions)
    frac = df.divide(total) if total>0 else df.astype(float)
    cols = ['H-bonds','Hydrophobic','Ionic','Water bridges']
    frac = frac[cols]

    # color mapping chosen to match typical PLIP legend (customize if you want)
    colors = {
        'H-bonds': '#4CAF50',        # green
        'Hydrophobic': '#9B59B6',    # purple
        'Ionic': '#FF1493',          # magenta/pink
        'Water bridges': '#1E88E5'   # blue
    }
    color_list = [colors[c] for c in frac.columns]

    fig, ax = plt.subplots(figsize=(12,5))
    frac.plot(kind='bar', stacked=True, ax=ax, color=color_list, width=0.8)
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_ylabel("Interactions Fraction")
    ax.set_xlabel("Residue")
    ax.legend(title="Interaction type", loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=4)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Saved plot to: {out_png}")

# -------- main CLI --------
def main(json_path, out_png, out_csv):
    with open(json_path) as f:
        plip_json = json.load(f)
    normalized = normalize_json(plip_json)
    df, total = build_counts(normalized)
    df.to_csv(out_csv)
    print(f"Wrote counts CSV to {out_csv}. Total interactions = {total}")
    plot_stacked_fraction(df, total, out_png)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python visualize_plip_json.py report.json out_plot.png out_counts.csv")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])