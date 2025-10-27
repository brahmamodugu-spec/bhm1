# app.py
# 🧬 Sir Brahmam Labs — DNA QC & Reflective Symmetry (BDNT) Pro+
# Elegant UI, live theme, STATIC photo (NO DNA overlay), QR code, mirror map,
# JSON report always; PDF report only if ReportLab is available.
# FASTA streaming with strict 200MB cap; optional FASTQ quality; simple login.

import io, os, re, json, math, gzip, mmap, tempfile, base64, datetime
from typing import Optional, Tuple, Iterable, List

import streamlit as st
from PIL import Image
import qrcode
import matplotlib.pyplot as plt

# ---------- Try ReportLab; fall back cleanly if missing ----------
REPORTLAB_OK = True
try:
    from reportlab.pdfgen import canvas
    from reportlab.lib.utils import ImageReader
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
except Exception:
    REPORTLAB_OK = False

# ======================
# Limits & Defaults
# ======================
MAX_FILE_MB = 200
MAX_FILE_BYTES = MAX_FILE_MB * 1024 * 1024

DEFAULT_BRAND   = "Sir Brahmam Labs"
DEFAULT_TAGLINE = "Mathematical DNA Insight — Educational, Elegant, Empowering"

# Login (optional). Keep empty to disable login screen.
ADMIN_USER = "admin"
ADMIN_PASS = "brahmam"   # change this in production

# ======================
# Helpers
# ======================
def img_to_base64(path: str) -> Optional[str]:
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")

def inject_css(primary="#5B7FFF", accent="#22c55e", warn="#f59e0b", danger="#ef4444"):
    BG_GRADIENT = f"linear-gradient(135deg, {primary} 0%, {accent} 40%, {danger} 100%)"
    st.markdown(f"""
    <style>
    .header-band {{
      background: {BG_GRADIENT};
      padding: 26px 20px;
      border-radius: 18px;
      color: white !important;
    }}
    .header-title {{
      font-size: 38px;
      font-weight: 800;
      margin: 0;
      line-height: 1.05;
    }}
    .header-sub {{
      font-size: 16px;
      opacity: 0.95;
    }}
    .card {{
      background: rgba(255,255,255,0.65);
      backdrop-filter: blur(6px);
      border: 1px solid rgba(255,255,255,0.65);
      border-radius: 16px;
      padding: 18px 18px 12px;
      box-shadow: 0 6px 32px rgba(0,0,0,0.08);
      margin-bottom: 8px;
    }}
    .metric-pill {{
      display: inline-block;
      padding: 6px 12px;
      border-radius: 999px;
      font-weight: 700;
      margin-right: 8px;
    }}
    .pill-good {{ background: #dcfce7; color: #166534; }}
    .pill-warn {{ background: #fef9c3; color: #92400e; }}
    .pill-bad  {{ background: #fee2e2; color: #991b1b; }}
    .cta-btn {{
      display:inline-block;padding:10px 16px;border-radius:12px;
      background:{accent};color:white;font-weight:700;text-decoration:none;
    }}
    .small-note {{ opacity: 0.8; font-size: 13px; }}
    h3 span.emoji {{ margin-right:6px; }}
    hr.soft {{ border: none; border-top: 1px dashed rgba(0,0,0,0.2); margin: 12px 0; }}

    /* ---------- Avatar (NO DNA overlay) ---------- */
    .avatar-orbit {{
      position: relative;
      width: 220px;
      height: 220px;
      margin: 8px auto 14px auto;
      filter: drop-shadow(0 8px 24px rgba(0,0,0,0.15));
    }}
    .avatar-orbit .avatar {{
      position: absolute;
      top: 30px; left: 30px;
      width: 160px; height: 160px;
      border-radius: 50%;
      object-fit: cover;
      border: 4px solid rgba(255,255,255,0.9);
      box-shadow: 0 6px 28px rgba(0,0,0,0.2);
      background: #fff;
    }}
    .avatar-orbit .orbit-ring {{
      position: absolute;
      top: 0; left: 0;
      width: 220px; height: 220px;
      border-radius: 50%;
      border: 1.5px dashed rgba(255,255,255,0.8);
    }}
    </style>
    """, unsafe_allow_html=True)

# ======================
# BDNT Core (streaming & safe)
# ======================
DIGIT = {"A": 0, "T": 1, "C": 2, "G": 3}
COMP  = {"A": "T", "T": "A", "C": "G", "G": "C"}

def clean_line(line: str) -> str:
    return re.sub(r"[^ATCGN]", "", line.upper())

def open_possibly_gzip(file_like) -> Iterable[str]:
    head = file_like.read(2); file_like.seek(0)
    if head == b"\x1f\x8b":  # gzip
        with gzip.GzipFile(fileobj=file_like) as gz:
            for raw in gz:
                yield raw.decode("utf-8", errors="ignore")
    else:
        for raw in file_like:
            yield raw.decode("utf-8", errors="ignore") if isinstance(raw, bytes) else raw

def parse_fasta_to_tempfile(upload, max_bytes: int = MAX_FILE_BYTES) -> Tuple[str, int, int, int, int, int, int]:
    total = 0
    written = 0
    counts = {"A":0, "T":0, "C":0, "G":0, "N":0}
    longest_run = 0
    run_char = None
    run_len  = 0
    tf = tempfile.NamedTemporaryFile(delete=False, suffix=".dna", mode="wb")
    temp_path = tf.name
    try:
        for line in open_possibly_gzip(upload):
            if line.startswith(">"):
                continue
            s = clean_line(line)
            if not s:
                continue
            b = s.encode("ascii")
            if max_bytes is not None and written + len(b) > max_bytes:
                remain = max(0, max_bytes - written)
                if remain > 0:
                    tf.write(b[:remain])
                    for ch in s[:remain]:
                        total += 1
                        if ch in counts: counts[ch] += 1
                        if run_char == ch:
                            run_len += 1
                        else:
                            run_char = ch; run_len = 1
                        longest_run = max(longest_run, run_len)
                raise ValueError(f"Processed sequence truncated at {MAX_FILE_MB} MB limit.")
            tf.write(b)
            written += len(b)
            for ch in s:
                total += 1
                if ch in counts: counts[ch] += 1
                if run_char == ch:
                    run_len += 1
                else:
                    run_char = ch
                    run_len = 1
                longest_run = max(longest_run, run_len)
    finally:
        tf.close()
    return temp_path, total, counts["A"], counts["T"], counts["C"], counts["G"], counts["N"], longest_run

def compute_entropy(a, t, c, g):
    tot = a+t+c+g
    if tot == 0: return 0.0
    H = 0.0
    for x in (a,t,c,g):
        if x == 0: continue
        p = x / tot
        H -= p * math.log2(p)
    return H  # bits (max 2)

def reverse_complement_base(b: str) -> str:
    return COMP.get(b, "N")

def compute_RDI_from_tempfile(temp_path: str, length: int) -> float:
    if length == 0: return 1.0
    matches = 0
    with open(temp_path, "rb") as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        lo, hi = 0, length-1
        while lo <= hi:
            b_lo = chr(mm[lo]); b_hi = chr(mm[hi])
            if b_lo == reverse_complement_base(b_hi):
                matches += 1
            lo += 1; hi -= 1
        mm.close()
    ratio = matches / length
    return 1.0 - ratio

def base4_eval_small(seq: str) -> int:
    val = 0
    for ch in seq:
        val = val*4 + DIGIT.get(ch,0)
    return val

def is_palindrome_in_base(n: int, base: int) -> bool:
    if n == 0: return True
    digits = []
    x = n
    while x > 0:
        digits.append(x % base)
        x //= base
    return digits == digits[::-1]

def reflective_equilibrium_small(temp_path: str, length: int) -> Tuple[Optional[bool], Optional[bool]]:
    MAX_RE_LEN = 100_000
    if length == 0: return False, False
    if length > MAX_RE_LEN: return None, None
    with open(temp_path,"rb") as f:
        seq = f.read().decode("ascii")
    rc = "".join(reverse_complement_base(ch) for ch in reversed(seq))
    Ex = base4_eval_small(seq)
    Ey = base4_eval_small(rc)
    R  = Ex * Ey
    pal10 = is_palindrome_in_base(R,10)
    pal4  = is_palindrome_in_base(R,4)
    return pal10, pal4

def parse_fastq_quality(upload) -> Tuple[Optional[float], Optional[float]]:
    upload.seek(0)
    state=0; seq_len=0
    total_q=0; total_bases=0; total_q30=0
    for line in open_possibly_gzip(upload):
        line=line.rstrip("\n")
        if state==0:
            if not line.startswith("@"): continue
            state=1
        elif state==1:
            seq_len=len(line.strip()); state=2
        elif state==2:
            state=3
        elif state==3:
            qual=line.strip()
            if len(qual)==seq_len:
                for ch in qual:
                    q=ord(ch)-33
                    total_q+=q; total_bases+=1
                    if q>=30: total_q30+=1
            state=0
    if total_bases==0: return None, None
    return total_q/total_bases, total_q30/total_bases

# ======================
# Mirror Map (heatline)
# ======================
def mirror_match_series(temp_path: str, length: int, sample_points: int = 5000) -> List[int]:
    if length == 0: return []
    idxs = list(range(length))
    if length > sample_points:
        step = length / sample_points
        idxs = [int(i*step) for i in range(sample_points)]
    with open(temp_path, "rb") as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        series = []
        for i in idxs:
            j = length-1-i
            b_i = chr(mm[i]); b_j = chr(mm[j])
            series.append(1 if b_i == reverse_complement_base(b_j) else 0)
        mm.close()
    return series

def plot_mirror_heatline(series: List[int], title: str = "Mirror Map (1=match, 0=mismatch)") -> None:
    plt.figure(figsize=(8, 2.5))
    plt.plot(series, linewidth=1.0)
    plt.title(title)
    plt.xlabel("Position (sampled)")
    plt.ylabel("Match")
    plt.ylim(-0.1, 1.1)
    st.pyplot(plt.gcf())
    plt.close()

# ======================
# PDF Report (only if ReportLab is available)
# ======================
if REPORTLAB_OK:
    def make_pdf_bytes(analysis: dict, img: Optional[Image.Image], primary_hex: str) -> bytes:
        buf = io.BytesIO()
        c = canvas.Canvas(buf, pagesize=A4)
        W,H = A4
        margin=40

        rgb = tuple(int(primary_hex.lstrip("#")[i:i+2],16)/255 for i in (0,2,4))
        c.setFillColorRGB(*rgb)
        c.rect(0, H-80, W, 80, stroke=0, fill=1)
        c.setFillColor(colors.white)
        c.setFont("Helvetica-Bold", 20)
        c.drawString(margin, H-50, f"{analysis.get('brand', DEFAULT_BRAND)} — DNA QC & BDNT Report")
        c.setFont("Helvetica", 11)
        c.drawString(margin, H-66, f"{analysis.get('tagline', DEFAULT_TAGLINE)}")

        if img is not None:
            img2 = img.copy(); img2.thumbnail((100,100))
            c.drawImage(ImageReader(img2), W-margin-100, H-110, width=90, height=90, mask='auto')

        y = H - 110 - margin
        c.setFillColorRGB(0,0,0)

        def line(text, size=12, bold=False, color=None, gap=16):
            nonlocal y
            if color is not None:
                c.setFillColor(color)
            else:
                c.setFillColorRGB(0,0,0)
            c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
            c.drawString(margin, y, text)
            y -= gap

        line(f"Institute/Brand: {analysis.get('org','')}", 11)
        line(f"Lead: {analysis.get('person','')}", 11)
        line(f"Contact: {analysis.get('phone','')} | {analysis.get('email','')}", 11)
        line(f"Website: {analysis.get('website','')}", 11)
        line(f"Report Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}", 11, gap=22)

        line("Results", 14, bold=True, color=colors.HexColor(primary_hex), gap=18)
        line(f"Length (bp): {analysis['length_bp']}")
        line(f"GC%: {analysis['GC_percent']}")
        line(f"N%: {analysis['N_percent']}")
        line(f"Shannon entropy (bits): {analysis['entropy_bits']}")
        line(f"Longest homopolymer: {analysis['longest_homopolymer']}")
        line(f"Reflective Equilibrium base-10: {analysis['reflective_equilibrium_base10']}")
        line(f"Reflective Equilibrium base-4: {analysis['reflective_equilibrium_base4']}")
        line(f"RDI (0 good → 1 bad): {analysis['RDI']}")
        if analysis.get("mean_phred_Q") is not None:
            line(f"Mean Phred Q (FASTQ): {analysis['mean_phred_Q']}")
        if analysis.get("frac_Q30") is not None:
            line(f"% bases ≥ Q30 (FASTQ): {round(analysis['frac_Q30']*100,2)}%")

        y -= 10
        c.setFillColorRGB(0.27,0.27,0.27)
        text = ("This report is for education/research demonstration. "
                "It does not diagnose disease. For clinical decisions, use validated "
                "variant-calling pipelines and ACMG interpretation.")
        c.setFont("Helvetica-Oblique", 10)
        for linewrap in [text[i:i+95] for i in range(0,len(text),95)]:
            c.drawString(margin, y, linewrap); y -= 12

        c.showPage(); c.save()
        buf.seek(0)
        return buf.getvalue()
else:
    def make_pdf_bytes(*_, **__):
        raise RuntimeError("ReportLab not installed")

# ======================
# App Init & Theme
# ======================
st.set_page_config(page_title="🧬 Sir Brahmam Labs — DNA QC & BDNT", page_icon="🧬", layout="wide")

with st.sidebar:
    st.markdown("### 🎨 Theme")
    PRIMARY = st.color_picker("Primary", "#5B7FFF")
    ACCENT  = st.color_picker("Accent",  "#22c55e")
    WARN    = "#f59e0b"
    DANGER  = "#ef4444"
inject_css(PRIMARY, ACCENT, WARN, DANGER)

# Optional login gate
if ADMIN_USER and ADMIN_PASS:
    if "auth" not in st.session_state:
        st.session_state.auth = False
    if not st.session_state.auth:
        with st.form("login"):
            st.markdown('<div class="card"><h3>🔐 Login</h3>', unsafe_allow_html=True)
            u = st.text_input("Username", "")
            p = st.text_input("Password", "", type="password")
            submitted = st.form_submit_button("Enter")
            st.markdown("</div>", unsafe_allow_html=True)
            if submitted:
                if u == ADMIN_USER and p == ADMIN_PASS:
                    st.session_state.auth = True
                else:
                    st.error("Invalid credentials")
        if not st.session_state.auth:
            st.stop()

# Header
st.markdown(f"""
<div class="header-band">
  <div class="header-title">🧬 {DEFAULT_BRAND}</div>
  <div class="header-sub">{DEFAULT_TAGLINE}</div>
</div>
""", unsafe_allow_html=True)
st.write("")

# ======================
# Sidebar — Static Photo (NO DNA overlay) + Business + QR
# ======================
with st.sidebar:
    st.markdown("### 👨‍🏫 Your Profile")
    photo64 = img_to_base64("brand_photo.jpg")
    if photo64:
        st.markdown(
            f"""
<div class="avatar-orbit">
  <div class="orbit-ring"></div>
  <img class="avatar" src="data:image/jpeg;base64,{photo64}" alt="Profile"/>
</div>
""",
            unsafe_allow_html=True,
        )
    else:
        st.warning("⚠️ Place your static photo as **brand_photo.jpg** in the app folder.")
        up_photo = st.file_uploader("Upload your photo (JPG/PNG)", type=["jpg","jpeg","png"])
        if up_photo is not None:
            img = Image.open(up_photo).convert("RGB")
            img.save("brand_photo.jpg")
            st.success("Saved as brand_photo.jpg. Reload the page to see your photo.")

    st.markdown("### 🏷️ Business Info")
    brand  = st.text_input("Brand", DEFAULT_BRAND)
    tagline= st.text_input("Tagline", DEFAULT_TAGLINE)
    org    = st.text_input("Institute / Brand", "New Era High School, Panchgani")
    person = st.text_input("Owner / Lead", "Mr. Brahmam Odugu")
    phone  = st.text_input("Contact Number", "+91-XXXXXXXXXX")
    email  = st.text_input("Email", "info@example.com")
    website= st.text_input("Website", "https://example.com")

    st.markdown("### 🔗 QR Code")
    qr_link = st.text_input("QR link (WhatsApp/Website)", "https://wa.me/91XXXXXXXXXX")
    if st.button("Generate QR"):
        qr_img = qrcode.make(qr_link)
        st.image(qr_img, caption="Scan me", use_container_width=True)
        st.session_state.qr_img = qr_img

    st.markdown("""<hr class="soft" />""", unsafe_allow_html=True)
    st.caption("Educational tool — not for medical diagnosis.")

# ======================
# Tabs
# ======================
tab_home, tab_analyze, tab_report, tab_about = st.tabs(
    ["🏠 Home", "🧪 DNA Analysis", "📄 Report Generator", "💼 Commercial"]
)

# HOME
with tab_home:
    st.markdown(
        f"""
<div class="card">
  <h3><span class="emoji">✨</span> Welcome</h3>
  <p><b>{brand}</b> blends mathematics and genomics to offer an educational view of DNA integrity and reflective symmetry (BDNT).</p>
  <ul>
    <li>Plain-English results everyone can understand</li>
    <li>Handles FASTA via streaming (≤ {MAX_FILE_MB} MB cap)</li>
    <li>Mirror map visualization (mismatch heatline)</li>
    <li>JSON report download; PDF export if ReportLab available</li>
    <li>Your branding, colors & photo</li>
  </ul>
  <a class="cta-btn" href="#dna-analysis">Start DNA Analysis</a>
</div>
""",
        unsafe_allow_html=True,
    )

# ANALYSIS
with tab_analyze:
    st.markdown('<a name="dna-analysis"></a>', unsafe_allow_html=True)
    st.markdown(f"""<div class="card"><h3><span class="emoji">🧪</span> Upload & Run</h3>""", unsafe_allow_html=True)

    colU1, colU2 = st.columns([2,2])
    with colU1:
        fasta = st.file_uploader(
            f"👉 Upload DNA (FASTA/plain; .gz supported) — Max {MAX_FILE_MB} MB",
            type=["fa","fna","fasta","txt","gz"]
        )
    with colU2:
        fastq = st.file_uploader(
            f"Optional FASTQ for base quality (.gz supported) — Max {MAX_FILE_MB} MB",
            type=["fastq","fq","gz"]
        )

    # Size check
    def _reject_if_too_big(f, label):
        if f is None:
            return False
        size = getattr(f, "size", None)
        if size is not None and size > MAX_FILE_BYTES:
            st.error(f"{label} is {size/1024/1024:.1f} MB, exceeding the {MAX_FILE_MB} MB limit.")
            return True
        return False

    too_big = _reject_if_too_big(fasta, "FASTA") or _reject_if_too_big(fastq, "FASTQ")

    st.markdown("**Mirror Map Options**")
    colM1, colM2 = st.columns([2,2])
    with colM1:
        sample_pts = st.slider("Sample points for Mirror Map", min_value=500, max_value=20000, value=5000, step=500)
    with colM2:
        show_mirror_map = st.checkbox("Show Mirror Map (heatline)", value=True)

    run = st.button("🚀 Run Analysis")
    st.markdown("</div>", unsafe_allow_html=True)

    if run:
        if fasta is None:
            st.error("Please upload a FASTA/plain DNA file.")
            st.stop()
        if too_big:
            st.stop()

        st.info(f"Reading & analyzing DNA (streaming, hard limit {MAX_FILE_MB} MB)…")
        fasta.seek(0)
        try:
            temp_path, length, A, T, C, G, N, longest_hp = parse_fasta_to_tempfile(fasta, max_bytes=MAX_FILE_BYTES)
        except ValueError as e:
            st.error(str(e))
            st.stop()

        if length==0:
            st.error("No valid A/T/C/G/N bases found.")
            st.stop()

        gc = (G+C)/length*100.0
        n_pct = (N/length*100.0)
        H = compute_entropy(A,T,C,G)
        rdi = compute_RDI_from_tempfile(temp_path, length)
        pal10, pal4 = reflective_equilibrium_small(temp_path, length)

        mean_q = frac_q30 = None
        if fastq is not None:
            if getattr(fastq, "size", 0) > MAX_FILE_BYTES:
                st.warning("FASTQ size exceeds limit; skipping quality metrics.")
            else:
                fastq.seek(0)
                mean_q, frac_q30 = parse_fastq_quality(fastq)

        # Friendly outputs
        st.markdown(f"""
<div class="card">
  <h3><span class="emoji">🔬</span> 1️⃣ Length (bp): {length}</h3>
  <p>Your DNA sequence is <b>{length}</b> base pairs long — suitable for QC and symmetry analysis within the cap.</p>

  <h3><span class="emoji">🧬</span> 2️⃣ GC% = {gc:.2f}</h3>
  <p><b>G + C = {gc:.2f}%</b>. Normal biological DNA often lies in <b>30–70%</b>.</p>
  <p>{"✅ Your GC content is within this normal zone." if 30<=gc<=70 else '<span class="metric-pill pill-warn">⚠️ Outside typical range</span>'}</p>

  <h3><span class="emoji">🧫</span> 3️⃣ N% = {n_pct:.2f}</h3>
  <p>‘N’ are unknown bases.</p>
  <p>{"✅ No ambiguous bases detected." if n_pct==0 else '<span class="metric-pill pill-warn">⚠️ Contains ambiguous bases</span>'}</p>

  <h3><span class="emoji">📊</span> 4️⃣ Shannon entropy = {H:.3f} bits</h3>
  <p>Maximum is <b>2.000 bits</b> for perfectly even A/T/C/G.</p>
  <p>{"✅ Near-ideal diversity." if abs(H-2.0)<=0.05 else '<span class="metric-pill pill-warn">⚠️ Low diversity</span>'}</p>

  <h3><span class="emoji">🧩</span> 5️⃣ Longest homopolymer = {longest_hp}</h3>
  <p>A homopolymer is a run like “AAAAA”.</p>
  <p>{"✅ Within normal range (≤ 6 typical)." if longest_hp<=6 else ('<span class="metric-pill pill-warn">⚠️ Long</span>' if longest_hp<10 else '<span class="metric-pill pill-bad">❌ Very long</span>')}</p>
</div>
""", unsafe_allow_html=True)

        # RE + RDI
        st.markdown("<div class='card'>", unsafe_allow_html=True)
        if pal10 is None or pal4 is None:
            st.markdown(f"""
  <h3><span class="emoji">💡</span> 6️⃣ Reflective Equilibrium tests</h3>
  <p>For very large sequences, the exact reflective product is too large to compute reliably in-browser.</p>
  <p><b>We skip palindromicity test</b> and rely on RDI (below) for symmetry deviation.</p>
""", unsafe_allow_html=True)
        else:
            st.markdown(f"""
  <h3><span class="emoji">💡</span> 6️⃣ Reflective Equilibrium (base-10): {"<span class='metric-pill pill-good'>PASS ✅</span>" if pal10 else "<span class='metric-pill pill-bad'>FAIL ❌</span>"}</h3>
  <h3><span class="emoji">💡</span> 7️⃣ Reflective Equilibrium (base-4): {"<span class='metric-pill pill-good'>PASS ✅</span>" if pal4 else "<span class='metric-pill pill-bad'>FAIL ❌</span>"}</h3>
  <p class="small-note">Computed via base-4 encoding × mirror-complement and palindromicity check.</p>
""", unsafe_allow_html=True)

        rdi_text = "✅ Very symmetric." if rdi<0.1 else ("🟡 Mild deviation." if rdi<0.3 else "⚠️ Strong asymmetry.")
        pill_cls = "pill-good" if rdi<0.3 else "pill-warn"
        st.markdown(f"""
  <h3><span class="emoji">⚖️</span> 8️⃣ Reflective Deviation Index (RDI) = {rdi:.3f}</h3>
  <p>RDI ranges 0 → 1 (0 = perfect mirror complement symmetry).</p>
  <p><span class="metric-pill {pill_cls}">{rdi_text}</span></p>
</div>
""", unsafe_allow_html=True)

        # Mirror Map
        if show_mirror_map:
            st.markdown("""<div class="card"><h3><span class="emoji">🗺️</span> Mirror Map</h3>""", unsafe_allow_html=True)
            series = mirror_match_series(temp_path, length, sample_points=sample_pts)
            plot_mirror_heatline(series, title=f"Mirror Map (sampled {len(series)} points)")
            st.caption("1 = mirror-complement match, 0 = mismatch (sampled positions across the whole sequence).")
            st.markdown("</div>", unsafe_allow_html=True)

        # Summary table
        status_gc = "✅" if 30 <= gc <= 70 else "⚠️"
        status_n  = "✅" if n_pct == 0 else "⚠️"
        status_H  = "✅" if abs(H - 2.0) <= 0.05 else "⚠️"
        status_hp = "✅" if longest_hp <= 6 else ("⚠️" if longest_hp < 10 else "❌")
        status_re = "—" if pal10 is None else ("✅" if (pal10 or pal4) else "❌")
        status_rdi= "✅" if rdi < 0.3 else "⚠️"

        st.markdown(f"""
<div class="card">
  <h3><span class="emoji">🧠</span> Summary</h3>
</div>
""", unsafe_allow_html=True)
        st.table({
            "Metric": ["GC %", "N %", "Entropy", "Homopolymer", "Reflective Eq.", "RDI"],
            "Interpretation": [
                f"{gc:.1f}%",
                f"{n_pct:.2f}%",
                f"{H:.3f} bits",
                f"{longest_hp}",
                ("(skipped for very large)" if pal10 is None else ("PASS" if (pal10 or pal4) else "FAIL")),
                f"{rdi:.3f}"
            ],
            "Status": [status_gc, status_n, status_H, status_hp, status_re, status_rdi]
        })

        # Educational alert (non-diagnostic)
        reasons = []
        if not(30<=gc<=70): reasons.append("Unusual GC%")
        if n_pct>5: reasons.append("High N%")
        if longest_hp>=10: reasons.append("Very long homopolymer")
        if rdi>=0.30: reasons.append("Strong reflective asymmetry (high RDI)")
        if reasons:
            st.error("Educational alert (non-diagnostic): " + "; ".join(reasons))
        else:
            st.success("No QC red flags detected. (Still not medical clearance.)")

        # Save for report tab
        st.session_state.analysis = {
            "brand": brand, "tagline": tagline,
            "length_bp": int(length), "GC_percent": round(gc,3), "N_percent": round(n_pct,3),
            "entropy_bits": round(H,3), "longest_homopolymer": int(longest_hp),
            "reflective_equilibrium_base10": None if pal10 is None else bool(pal10),
            "reflective_equilibrium_base4": None if pal4 is None else bool(pal4),
            "RDI": round(rdi,3),
            "mean_phred_Q": None if mean_q is None else round(float(mean_q),3),
            "frac_Q30": None if frac_q30 is None else round(float(frac_q30),4),
            "org": org, "person": person, "phone": phone, "email": email, "website": website
        }

        # Cleanup temp
        try: os.remove(temp_path)
        except Exception: pass

# REPORT
with tab_report:
    st.markdown(f"""<div class="card"><h3><span class="emoji">📄</span> Report</h3>""", unsafe_allow_html=True)
    st.caption("Generate downloads after you run an analysis.")
    if "analysis" not in st.session_state:
        st.warning("Run an analysis first (see the DNA Analysis tab).")
    else:
        data = st.session_state.analysis

        # JSON
        st.download_button(
            "⬇️ Download JSON Report",
            data=json.dumps(data, indent=2),
            file_name="dna_report.json",
            mime="application/json"
        )

        # PDF (only if ReportLab present)
        if REPORTLAB_OK:
            img_for_pdf = None
            if os.path.exists("brand_photo.jpg"):
                img_for_pdf = Image.open("brand_photo.jpg")
            if st.button("🧾 Generate PDF"):
                pdf_bytes = make_pdf_bytes(data, img_for_pdf, PRIMARY)
                st.download_button("⬇️ Download PDF Report", data=pdf_bytes,
                                   file_name="dna_report.pdf", mime="application/pdf")
        else:
            st.info("PDF export is unavailable because ReportLab is not installed in this environment. "
                    "Use the JSON download above, or install ReportLab to enable PDF.")

    st.markdown("</div>", unsafe_allow_html=True)

# COMMERCIAL
with tab_about:
    st.markdown(f"""
<div class="card">
  <h3><span class="emoji">💼</span> {brand} — Commercial</h3>
  <p><b>Services</b></p>
  <ul>
    <li>DNA QC & Reflective Symmetry (BDNT) educational reports</li>
    <li>School workshops & live demonstrations (biology + mathematics)</li>
    <li>Custom curriculum content & student competitions</li>
    <li>Institutional branding on scientific reports</li>
  </ul>
  <hr class="soft" />
  <p><b>Packages</b></p>
  <ul>
    <li><b>Starter</b> — 10 reports/month, standard branding</li>
    <li><b>Pro</b> — 50 reports/month, custom branding + workshop</li>
    <li><b>Institution</b> — Unlimited reports, annual training & events</li>
  </ul>
  <hr class="soft" />
  <p><b>Contact</b><br/>
     {org}<br/>
     {person}<br/>
     {phone} | {email}<br/>
     {website}
  </p>
  <p><b>Scan to connect:</b></p>
</div>
""", unsafe_allow_html=True)

    if st.session_state.get("qr_img") is not None:
        st.image(st.session_state.qr_img, width=180)
    else:
        st.caption("Use the sidebar to generate a QR code to your WhatsApp or website.")

# Footer
st.caption("© {} {} • For education and demonstrations; not a clinical diagnostic."
           .format(datetime.datetime.now().year, brand))
