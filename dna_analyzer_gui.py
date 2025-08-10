import tkinter as tk
from tkinter import messagebox, scrolledtext, filedialog
import matplotlib.pyplot as plt
import csv
import os
import json

# Constants
NUCLEOTIDE_WEIGHTS = {'A':313.21, 'T':304.2, 'G':329.21, 'C':289.18}

# Genetic code dictionary for translation (Standard)
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def count_nucleotides(dna_seq):
    dna_seq = dna_seq.upper()
    return {nuc: dna_seq.count(nuc) for nuc in "ATGC"}

def gc_content(dna_seq):
    dna_seq = dna_seq.upper()
    g = dna_seq.count('G')
    c = dna_seq.count('C')
    return round(((g + c) / len(dna_seq)) * 100, 2) if len(dna_seq) > 0 else 0

def complementary_strand(dna_seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    dna_seq = dna_seq.upper()
    return ''.join(complement.get(base, '') for base in dna_seq)

def reverse_complement(dna_seq):
    return complementary_strand(dna_seq)[::-1]

def find_repeats(dna_seq, length):
    repeats = {}
    dna_seq = dna_seq.upper()
    for i in range(len(dna_seq) - length + 1):
        subseq = dna_seq[i:i+length]
        repeats[subseq] = repeats.get(subseq, 0) + 1
    return {seq: count for seq, count in repeats.items() if count > 1}

def molecular_weight(dna_seq):
    counts = count_nucleotides(dna_seq)
    return round(sum(NUCLEOTIDE_WEIGHTS[nuc]*count for nuc, count in counts.items()), 2)

def nucleotide_percentage(dna_seq):
    counts = count_nucleotides(dna_seq)
    length = len(dna_seq)
    if length == 0:
        return {nuc: 0 for nuc in counts}
    return {nuc: round((count / length) * 100, 2) for nuc, count in counts.items()}

def plot_bar_chart(counts):
    nucleotides = list(counts.keys())
    values = list(counts.values())
    plt.bar(nucleotides, values, color=['blue','orange','green','red'])
    plt.title('Nucleotide Counts')
    plt.xlabel('Nucleotide')
    plt.ylabel('Count')
    plt.show()

def translate_dna(dna_seq):
    dna_seq = dna_seq.upper()
    protein = []
    for i in range(0, len(dna_seq)-2, 3):
        codon = dna_seq[i:i+3]
        protein.append(CODON_TABLE.get(codon, 'X'))  # X for unknown codon
    return ''.join(protein)

def find_orfs(dna_seq, min_length=30):
    seq = dna_seq.upper()
    orfs = []
    for frame in range(3):
        start = None
        for i in range(frame, len(seq)-2, 3):
            codon = seq[i:i+3]
            if codon == "ATG" and start is None:
                start = i
            if codon in ("TAA", "TAG", "TGA") and start is not None:
                length = i + 3 - start
                if length >= min_length:
                    orfs.append((start+1, i+3, seq[start:i+3], translate_dna(seq[start:i+3])))
                start = None
        # If no stop codon found till end
        if start is not None and (len(seq) - start) >= min_length:
            orfs.append((start+1, len(seq), seq[start:], translate_dna(seq[start:])))
    return orfs

def reverse_sequence(dna_seq):
    return dna_seq[::-1]

# UI Functions

def analyze():
    dna_seq = entry.get().upper().strip()
    if len(dna_seq) == 0:
        set_status("Please enter a DNA sequence.", "red")
        return
    if any(nuc not in "ATGC" for nuc in dna_seq):
        set_status("Invalid sequence! Use only A, T, G, C.", "red")
        return

    set_status("Sequence valid. Analyzing...", "green")
    counts = count_nucleotides(dna_seq)
    gc = gc_content(dna_seq)
    comp = complementary_strand(dna_seq)
    rev_comp = reverse_complement(dna_seq)
    repeats_2 = find_repeats(dna_seq, 2)
    repeats_3 = find_repeats(dna_seq, 3)
    weight = molecular_weight(dna_seq)
    percentages = nucleotide_percentage(dna_seq)
    protein = translate_dna(dna_seq)
    orfs = find_orfs(dna_seq)

    output_text.config(state=tk.NORMAL)
    output_text.delete(1.0, tk.END)
    output_text.insert(tk.END, f"Nucleotide Counts:\n")
    for nuc, count in counts.items():
        output_text.insert(tk.END, f"  {nuc}: {count}\n")

    output_text.insert(tk.END, f"\nNucleotide Percentages:\n")
    for nuc, perc in percentages.items():
        output_text.insert(tk.END, f"  {nuc}: {perc}%\n")

    output_text.insert(tk.END, f"\nGC Content: {gc}%\n")
    output_text.insert(tk.END, f"Molecular Weight (approx): {weight} Da\n")
    output_text.insert(tk.END, f"Complementary Strand: {comp}\n")
    output_text.insert(tk.END, f"Reverse Complement Strand: {rev_comp}\n")

    output_text.insert(tk.END, "\nRepeated sequences (length 2):\n")
    if repeats_2:
        for seq, count in repeats_2.items():
            output_text.insert(tk.END, f"  {seq}: {count} times\n")
    else:
        output_text.insert(tk.END, "  None found\n")

    output_text.insert(tk.END, "\nRepeated sequences (length 3):\n")
    if repeats_3:
        for seq, count in repeats_3.items():
            output_text.insert(tk.END, f"  {seq}: {count} times\n")
    else:
        output_text.insert(tk.END, "  None found\n")

    output_text.insert(tk.END, "\nProtein Translation:\n")
    output_text.insert(tk.END, f"{protein}\n")

    output_text.insert(tk.END, "\nOpen Reading Frames (ORFs) >= 30 nt:\n")
    if orfs:
        for start, end, seq, prot in orfs:
            output_text.insert(tk.END, f"  Positions {start}-{end}, Protein: {prot}\n")
    else:
        output_text.insert(tk.END, "  None found\n")

    output_text.insert(tk.END, "\nNucleotide Counts Visualization:\n")
    for nuc, count in counts.items():
        output_text.insert(tk.END, f"{nuc}: {'#' * count}\n")

    output_text.config(state=tk.DISABLED)

    plot_bar_chart(counts)

def save_results_csv():
    content = output_text.get(1.0, tk.END).strip()
    if not content:
        messagebox.showwarning("Warning", "No analysis results to save.")
        return

    file_path = filedialog.asksaveasfilename(
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        title="Save Analysis Results"
    )
    if not file_path:
        return

    try:
        with open(file_path, mode='w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            lines = content.splitlines()
            for line in lines:
                if ':' in line:
                    parts = [p.strip() for p in line.split(':', 1)]
                    writer.writerow(parts)
                else:
                    writer.writerow([line])
        messagebox.showinfo("Success", f"Results saved as {os.path.basename(file_path)}")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to save file:\n{e}")

def clear_all():
    entry.delete(0, tk.END)
    output_text.config(state=tk.NORMAL)
    output_text.delete(1.0, tk.END)
    output_text.config(state=tk.DISABLED)
    set_status("", "black")

def set_status(msg, color):
    status_label.config(text=msg, fg=color)

def show_help():
    help_text = (
        "DNA Analyzer Help\n\n"
        "1. Enter a DNA sequence (A, T, G, C only) or load from file.\n"
        "2. Click Analyze to see nucleotide counts, GC content, molecular weight, complements, repeats, translation, ORFs, and visualization.\n"
        "3. Use buttons to save results, copy results or sequences, reverse input, save/load project.\n\n"
        "Created by YourName"
    )
    messagebox.showinfo("Help / About", help_text)

def load_file():
    file_path = filedialog.askopenfilename(
        title="Open DNA Sequence File",
        filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
    )
    if not file_path:
        return
    try:
        with open(file_path, 'r') as f:
            content = f.read().strip().replace('\n','').replace(' ','')
        entry.delete(0, tk.END)
        entry.insert(0, content)
        set_status(f"Loaded sequence from {os.path.basename(file_path)}", "green")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load file:\n{e}")

def copy_results():
    content = output_text.get(1.0, tk.END).strip()
    if not content:
        messagebox.showwarning("Warning", "No analysis results to copy.")
        return
    root.clipboard_clear()
    root.clipboard_append(content)
    set_status("Results copied to clipboard.", "green")

def copy_complement():
    comp = complementary_strand(entry.get().strip())
    if not comp:
        messagebox.showwarning("Warning", "No sequence to copy complement.")
        return
    root.clipboard_clear()
    root.clipboard_append(comp)
    set_status("Complementary strand copied to clipboard.", "green")

def copy_reverse_complement():
    rev_comp = reverse_complement(entry.get().strip())
    if not rev_comp:
        messagebox.showwarning("Warning", "No sequence to copy reverse complement.")
        return
    root.clipboard_clear()
    root.clipboard_append(rev_comp)
    set_status("Reverse complementary strand copied to clipboard.", "green")

def reverse_input():
    seq = entry.get().strip()
    if not seq:
        set_status("No sequence to reverse.", "red")
        return
    entry.delete(0, tk.END)
    entry.insert(0, reverse_sequence(seq))
    set_status("Sequence reversed.", "green")

def save_project():
    project_data = {
        "sequence": entry.get(),
        "results": output_text.get(1.0, tk.END)
    }
    file_path = filedialog.asksaveasfilename(
        defaultextension=".json",
        filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        title="Save Project"
    )
    if not file_path:
        return
    try:
        with open(file_path, 'w') as f:
            json.dump(project_data, f)
        set_status(f"Project saved to {os.path.basename(file_path)}", "green")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to save project:\n{e}")

def load_project():
    file_path = filedialog.askopenfilename(
        title="Open Project",
        filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
    )
    if not file_path:
        return
    try:
        with open(file_path, 'r') as f:
            project_data = json.load(f)
        entry.delete(0, tk.END)
        entry.insert(0, project_data.get("sequence", ""))
        output_text.config(state=tk.NORMAL)
        output_text.delete(1.0, tk.END)
        output_text.insert(tk.END, project_data.get("results", ""))
        output_text.config(state=tk.DISABLED)
        set_status(f"Project loaded from {os.path.basename(file_path)}", "green")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load project:\n{e}")

def on_entry_change(event):
    set_status("", "black")

# --- GUI setup ---

root = tk.Tk()
root.title("DNA Sequence Analyzer")

label_font = ("Arial", 11)
button_font = ("Arial", 10)

frame_top = tk.Frame(root)
frame_top.pack(padx=10, pady=10)

tk.Label(frame_top, text="Enter DNA Sequence (A, T, G, C only):", font=label_font).pack(anchor='w')
entry = tk.Entry(frame_top, width=60, font=label_font)
entry.pack(pady=5)
entry.bind("<KeyRelease>", on_entry_change)

btn_frame = tk.Frame(frame_top)
btn_frame.pack(pady=5)

analyze_btn = tk.Button(btn_frame, text="Analyze", command=analyze, font=button_font, width=12)
analyze_btn.grid(row=0, column=0, padx=4, pady=4)

save_btn = tk.Button(btn_frame, text="Save Results", command=save_results_csv, font=button_font, width=12)
save_btn.grid(row=0, column=1, padx=4, pady=4)

copy_btn = tk.Button(btn_frame, text="Copy Results", command=copy_results, font=button_font, width=12)
copy_btn.grid(row=0, column=2, padx=4, pady=4)

load_btn = tk.Button(btn_frame, text="Load File", command=load_file, font=button_font, width=12)
load_btn.grid(row=0, column=3, padx=4, pady=4)

clear_btn = tk.Button(btn_frame, text="Clear", command=clear_all, font=button_font, width=12)
clear_btn.grid(row=0, column=4, padx=4, pady=4)

help_btn = tk.Button(btn_frame, text="Help / About", command=show_help, font=button_font, width=12)
help_btn.grid(row=0, column=5, padx=4, pady=4)

# New row for more features
btn_frame2 = tk.Frame(frame_top)
btn_frame2.pack(pady=5)

copy_comp_btn = tk.Button(btn_frame2, text="Copy Complement", command=copy_complement, font=button_font, width=15)
copy_comp_btn.grid(row=0, column=0, padx=4, pady=4)

copy_revcomp_btn = tk.Button(btn_frame2, text="Copy Rev. Complement", command=copy_reverse_complement, font=button_font, width=18)
copy_revcomp_btn.grid(row=0, column=1, padx=4, pady=4)

reverse_seq_btn = tk.Button(btn_frame2, text="Reverse Input", command=reverse_input, font=button_font, width=12)
reverse_seq_btn.grid(row=0, column=2, padx=4, pady=4)

save_proj_btn = tk.Button(btn_frame2, text="Save Project", command=save_project, font=button_font, width=12)
save_proj_btn.grid(row=0, column=3, padx=4, pady=4)

load_proj_btn = tk.Button(btn_frame2, text="Load Project", command=load_project, font=button_font, width=12)
load_proj_btn.grid(row=0, column=4, padx=4, pady=4)

output_text = scrolledtext.ScrolledText(root, width=70, height=20, font=("Courier New", 10))
output_text.pack(padx=10, pady=10)
output_text.config(state=tk.DISABLED)

status_label = tk.Label(root, text="", font=label_font)
status_label.pack(pady=5)

root.mainloop()
