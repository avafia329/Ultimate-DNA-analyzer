# app.py — Final Ultimate DNA Analyzer (no login)
import os, io, json, time, re
from datetime import datetime
from collections import Counter

from flask import Flask, render_template, request, jsonify, send_file, url_for
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Restriction import RestrictionBatch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader

# --- App init ---
app = Flask(__name__)
app.config['SECRET_KEY'] = 'replace-with-any-string'

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, 'static', 'outputs')
PROJECTS_DIR = os.path.join(BASE_DIR, 'projects')
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(PROJECTS_DIR, exist_ok=True)

# --- Constants ---
STANDARD_CODON_TABLE = CodonTable.unambiguous_dna_by_name["Standard"]
COMMON_ENZYMES = ['EcoRI', 'BamHI', 'HindIII', 'NotI', 'XhoI']
MW = {'A': 313.21, 'T': 304.2, 'G': 329.21, 'C': 289.18}

# keep last analysis in memory for quick PDF export
latest_results = {}

# ------------------ Helpers / Bio functions ------------------

def clean_seq(raw: str) -> str:
    if not raw: return ''
    s = re.sub('[^ATGCatgc]', '', raw)
    return s.upper()

def nucleotide_counts(seq: str):
    counts = {n: seq.count(n) for n in 'ATGC'}
    total = len(seq)
    percentages = {n: round((counts[n]/total)*100,2) if total else 0 for n in counts}
    return counts, percentages

def gc_content(seq: str):
    return round((seq.count('G') + seq.count('C')) / len(seq) * 100, 2) if seq else 0

def complement(seq: str):
    comp = {'A':'T','T':'A','G':'C','C':'G'}
    return ''.join(comp.get(b, b) for b in seq)

def reverse_complement(seq: str):
    return complement(seq)[::-1]

def molecular_weight(seq: str):
    return round(sum(MW.get(n,0) for n in seq), 2)

def find_repeats(seq: str, length: int):
    d = {}
    for i in range(len(seq)-length+1):
        sub = seq[i:i+length]
        d[sub] = d.get(sub,0)+1
    return {k:v for k,v in d.items() if v>1}

def translate_protein(seq: str):
    prot = ''
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        aa = STANDARD_CODON_TABLE.forward_table.get(codon, '')
        if aa == '': continue
        prot += aa
    return prot

def find_orfs(seq: str, min_prot_len=5):
    orfs = []
    L = len(seq)
    stops = {'TAA','TAG','TGA'}
    for frame in range(3):
        i = frame
        while i <= L-3:
            if seq[i:i+3] == 'ATG':
                j = i
                prot = ''
                while j <= L-3:
                    cod = seq[j:j+3]
                    if cod in stops:
                        break
                    aa = STANDARD_CODON_TABLE.forward_table.get(cod,'')
                    if aa == '':
                        break
                    prot += aa
                    j += 3
                if prot.startswith('M') and len(prot) >= min_prot_len:
                    orfs.append({'start': i+1, 'end': j+3, 'protein': prot})
                i = j+3
            else:
                i += 3
    return orfs

def gc_skew(seq: str, window=10):
    if len(seq) < window: return []
    out = []
    for i in range(len(seq)-window+1):
        w = seq[i:i+window]
        g = w.count('G'); c = w.count('C')
        out.append(round((g-c)/(g+c) if (g+c)!=0 else 0, 3))
    return out

def codon_usage(seq: str):
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    cnt = Counter([c for c in codons if len(c)==3])
    tot = sum(cnt.values()) or 1
    return {c: round((f/tot)*100,2) for c,f in cnt.items()}

def find_restriction_sites(seq: str):
    # RestrictionBatch accepts Seq or string; we use Seq for safety
    sb = Seq(seq)
    rb = RestrictionBatch(COMMON_ENZYMES)
    analysis = rb.search(sb)
    return {enzyme: positions for enzyme, positions in analysis.items() if positions}

def melting_temp(seq: str):
    try:
        return round(mt.Tm_Wallace(Seq(seq)), 2)
    except Exception:
        return None

def cpg_islands(seq: str, window=200):
    islands = []
    if len(seq) < window: return islands
    for i in range(0, len(seq)-window+1, window//2):
        w = seq[i:i+window]
        gc = (w.count('G')+w.count('C'))/window
        cpg = w.count('CG')
        expected = (w.count('C') * w.count('G'))/window if window else 0
        obs = (cpg/expected) if expected else 0
        if gc > 0.5 and obs > 0.6:
            islands.append((i+1, i+window))
    return islands

def primer_design(seq: str):
    primers = []
    for L in range(18,23):
        for i in range(0, len(seq)-L+1):
            frag = seq[i:i+L]
            if frag[0] not in ('A','T'): continue
            tm = melting_temp(frag)
            if tm and 55 <= tm <= 65:
                primers.append({'start': i+1, 'seq': frag, 'tm': tm})
                if len(primers) >= 5:
                    return primers
    return primers

# ---------------- plotting helpers (returning file paths) ----------------

def save_nuc_plot(counts, fname):
    labels = ['A','T','G','C']
    values = [counts.get(l,0) for l in labels]
    fig, ax = plt.subplots(figsize=(6,2.4))
    bars = ax.bar(labels, values, color=['#4caf50','#f44336','#2196f3','#ffeb3b'])
    ax.set_ylabel('Count'); ax.set_title('Nucleotide counts')
    for bar, val in zip(bars, values):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.1, str(val), ha='center', fontsize=9)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, fname)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path

def save_gc_plot(skew_vals, fname):
    fig, ax = plt.subplots(figsize=(6,2.4))
    if skew_vals:
        ax.plot(range(1,len(skew_vals)+1), skew_vals, color='purple')
        ax.set_title('GC Skew (window=10)')
        ax.set_xlabel('Window index')
        ax.set_ylabel('(G-C)/(G+C)')
    else:
        ax.text(0.5,0.5,'No GC skew (sequence too short)', ha='center', va='center')
        ax.axis('off')
    plt.tight_layout()
    path = os.path.join(OUT_DIR, fname)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path

def save_protein_preview(protein_seq, fname):
    # simple spiral-like 3d-ish image using colors by hydrophobicity
    fig = plt.figure(figsize=(5,4))
    try:
        ax = fig.add_subplot(111, projection='3d')
        n = len(protein_seq)
        if n == 0:
            ax.text(0.5,0.5,0.5,'No protein', ha='center', va='center')
        else:
            t = np.linspace(0, 4*np.pi, n)
            x = np.cos(t); y = np.sin(t); z = np.linspace(0,4,n)
            hydrophobic = set('AILMVFWY')
            colors = ['#d73027' if aa in hydrophobic else '#4575b4' for aa in protein_seq]
            ax.plot(x,y,z, color='gray', lw=1)
            ax.scatter(x,y,z, c=colors, s=25)
            ax.set_axis_off()
    except Exception:
        plt.text(0.5,0.5,'3D preview not available', ha='center', va='center')
        plt.axis('off')
    plt.tight_layout()
    path = os.path.join(OUT_DIR, fname)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path

# ---------------- project functions ----------------
def safe_name(name):
    s = re.sub(r'[^A-Za-z0-9_\- ]','_', str(name)).strip()
    return s or f"project_{int(time.time())}"

def save_project(name, data):
    fname = safe_name(name)+'.json'
    path = os.path.join(PROJECTS_DIR, fname)
    with open(path, 'w', encoding='utf8') as f:
        json.dump(data, f, indent=2)
    return fname

def list_projects():
    files = [f for f in os.listdir(PROJECTS_DIR) if f.endswith('.json')]
    files.sort(key=lambda x: os.path.getmtime(os.path.join(PROJECTS_DIR,x)), reverse=True)
    return files

def load_project(fname):
    path = os.path.join(PROJECTS_DIR, fname)
    if not os.path.exists(path):
        return None
    with open(path, 'r', encoding='utf8') as f:
        return json.load(f)

# ---------------- PDF generation ----------------

def generate_pdf(results):
    # create images in memory and embed them
    nuc_name = f"nuc_{int(time.time()*1000)}.png"
    gc_name = f"gc_{int(time.time()*1000)}.png"
    prot_name = f"prot_{int(time.time()*1000)}.png"
    nuc_path = save_nuc_plot(results.get('counts',{}), nuc_name)
    gc_path = save_gc_plot(results.get('gc_skew',[]), gc_name)
    prot_path = save_protein_preview(results.get('protein',''), prot_name)

    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=(595,842))
    c.setFont("Helvetica-Bold", 16)
    c.drawString(40,800,"DNA Analysis Report")
    c.setFont("Helvetica",10)
    c.drawString(40,782, f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    c.drawString(40,764, f"Sequence length: {len(results.get('sequence',''))}")
    y = 740
    c.setFont("Helvetica-Bold",12); c.drawString(40,y, "Nucleotide counts:")
    y -= 16; c.setFont("Helvetica",10)
    for b in ['A','T','G','C']:
        c.drawString(50,y, f"{b}: {results['counts'].get(b,0)} ({results['percentages'].get(b,0)}%)")
        y -= 14
    y -= 6
    c.drawString(40,y, f"GC content: {results.get('gc_content',0)}%")
    y -= 20
    try:
        img1 = ImageReader(nuc_path); c.drawImage(img1,40,y-140, width=260, height=120)
        img2 = ImageReader(gc_path); c.drawImage(img2,320,y-140, width=220, height=120)
        y -= 160
    except Exception:
        pass
    try:
        img3 = ImageReader(prot_path); c.drawImage(img3, 40, y-210, width=260, height=200); y -= 210
    except Exception:
        pass
    c.setFont("Helvetica-Bold",12); c.drawString(40,y, "ORFs:")
    y -= 16; c.setFont("Helvetica",10)
    for orf in results.get('orfs', [])[:10]:
        c.drawString(50,y, f"Start: {orf['start']}, End: {orf['end']}, Protein len: {len(orf['protein'])}")
        y -= 12
        if y < 80:
            c.showPage(); y = 760
    c.showPage(); c.save()
    buf.seek(0)
    return buf

# ---------------- Routes ----------------

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html', results=None)

@app.route('/analyze', methods=['POST'])
def analyze():
    global latest_results
    # accept JSON or form
    seq = ''
    if request.is_json:
        payload = request.get_json()
        seq = payload.get('sequence','') or payload.get('dna_seq','')
    else:
        seq = request.form.get('dna_seq','')

    seq = clean_seq(seq)
    if not seq:
        return jsonify({'error': 'No valid DNA sequence (A,T,G,C) provided.'}), 400

    # ensure Seq object where needed
    # perform analyses
    counts, percentages = nucleotide_counts(seq)
    gc = gc_content(seq)
    comp = complement(seq)
    rev_comp = reverse_complement(seq)
    mw = molecular_weight(seq)
    repeats2 = find_repeats(seq,2)
    repeats3 = find_repeats(seq,3)
    protein = translate_protein(seq)
    orfs = find_orfs(seq)
    gcsk = gc_skew(seq)
    codon_stats = codon_usage(seq)
    restriction_sites = find_restriction_sites(seq)
    mt_val = melting_temp(seq)
    cpg = cpg_islands(seq)
    primers = primer_design(seq)

    # produce images and store relative paths
    timestamp = int(time.time()*1000)
    nuc_img = f"nuc_{timestamp}.png"; gc_img = f"gc_{timestamp}.png"; prot_img = f"prot_{timestamp}.png"
    nuc_path = save_nuc_plot(counts, nuc_img)
    gc_path = save_gc_plot(gcsk, gc_img)
    prot_path = save_protein_preview(protein, prot_img)

    results = {
        'sequence': seq,
        'counts': counts,
        'percentages': percentages,
        'gc_content': gc,
        'complement': comp,
        'reverse_complement': rev_comp,
        'molecular_weight': mw,
        'repeats2': repeats2,
        'repeats3': repeats3,
        'protein': protein,
        'orfs': orfs,
        'codon_usage': codon_stats,
        'restriction_sites': restriction_sites,
        'melting_temp': mt_val,
        'cpg_islands': cpg,
        'primers': primers,
        'gc_skew': gcsk,
        'nuc_png': os.path.relpath(nuc_path, BASE_DIR).replace('\\','/'),
        'gc_png': os.path.relpath(gc_path, BASE_DIR).replace('\\','/'),
        'prot_png': os.path.relpath(prot_path, BASE_DIR).replace('\\','/'),
    }

    latest_results = results
    # return JSON (frontend expects JSON in AJAX mode)
    summary = f"Length: {len(seq)} | GC%: {gc} | ORFs: {len(orfs)}"
    return jsonify({'text': summary, 'counts': counts, 'percentages': percentages, 'gc_skew': gcsk, 'full': results})

@app.route('/feature/<name>', methods=['POST'])
def feature(name):
    payload = request.get_json(force=True)
    seq = clean_seq(payload.get('sequence',''))
    if not seq:
        return jsonify({'result': 'No valid sequence provided.'}), 400
    name = name.lower()
    if name == 'codon': return jsonify({'result': codon_usage(seq)})
    if name == 'restriction': return jsonify({'result': find_restriction_sites(seq)})
    if name == 'melting': return jsonify({'result': melting_temp(seq)})
    if name == 'cpg': return jsonify({'result': cpg_islands(seq)})
    if name == 'primer': return jsonify({'result': primer_design(seq)})
    if name == 'orf': return jsonify({'result': find_orfs(seq)})
    if name == 'protein': return jsonify({'result': translate_protein(seq)})
    return jsonify({'result': f'Unknown feature: {name}'}), 400

@app.route('/save_project', methods=['POST'])
def save_project_route():
    payload = request.get_json(force=True)
    name = payload.get('name') or f'project_{int(time.time())}'
    results = payload.get('results') or {}
    if not results:
        return jsonify({'ok': False, 'error': 'No results to save'})
    fname = save_project(name, results)
    return jsonify({'ok': True, 'filename': fname})

@app.route('/list_projects', methods=['GET'])
def list_projects_route():
    return jsonify({'projects': list_projects()})

@app.route('/load_project', methods=['GET'])
def load_project_route():
    name = request.args.get('name')
    if not name: return jsonify({'ok': False, 'error': 'Missing name'}), 400
    data = load_project(name)
    if not data: return jsonify({'ok': False, 'error': 'Not found'}), 404
    return jsonify({'ok': True, 'data': data})

@app.route('/download_pdf', methods=['POST'])
def download_pdf():
    # Accept JSON results or use latest_results
    try:
        payload = request.get_json(force=True)
        results = payload.get('results') or latest_results
    except Exception:
        results = latest_results
    if not results:
        return jsonify({'error': 'No analysis available to export.'}), 400
    pdf_buf = generate_pdf(results)
    return send_file(pdf_buf, as_attachment=True, download_name='dna_report.pdf', mimetype='application/pdf')

@app.route('/static/outputs/<path:fname>')
def serve_output(fname):
    path = os.path.join(OUT_DIR, fname)
    if not os.path.exists(path):
        return "Not found", 404
    return send_file(path, mimetype='image/png')

# ----- run -----
if __name__ == '__main__':
    app.run(debug=True)
