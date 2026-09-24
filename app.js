'use strict';

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fmt = n => n.toLocaleString('en-US');
const err = msg => `<p class="err">${esc(msg)}</p>`;

const COMP = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N', R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K', B: 'V', V: 'B', D: 'H', H: 'D' };
const revcomp = s => s.split('').reverse().map(c => COMP[c] || 'N').join('');

// Accepts raw sequence or FASTA; returns {seq} or {error}
function cleanDNA(text) {
  const body = text.split('\n').filter(l => !l.trim().startsWith('>')).join('');
  const s = body.replace(/[\s\d]/g, '').toUpperCase().replace(/U/g, 'T');
  if (!s) return { error: 'Please enter a DNA sequence.' };
  const bad = [...new Set(s.replace(/[ACGTN]/g, ''))];
  if (bad.length) return { error: `Invalid characters in the sequence: ${bad.join(' ')} (only A, C, G, T, N are allowed).` };
  return { seq: s };
}

function cleanProtein(text) {
  const body = text.split('\n').filter(l => !l.trim().startsWith('>')).join('');
  let s = body.replace(/[\s\d]/g, '').toUpperCase();
  const hadStop = /[*]/.test(s.replace(/[*]+$/, ''));
  s = s.replace(/[*\-]/g, '');
  if (!s) return { error: 'Please enter an amino acid sequence.' };
  const bad = [...new Set(s.replace(/[A-Z]/g, ''))];
  if (bad.length) return { error: `Invalid characters in the sequence: ${bad.join(' ')}` };
  return { seq: s, hadStop };
}

// Standard genetic code, codons ordered TCAG
const CODE = (() => {
  const b = 'TCAG', aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', t = {};
  let k = 0;
  for (const x of b) for (const y of b) for (const z of b) t[x + y + z] = aa[k++];
  return t;
})();
const translate = s => { let p = ''; for (let i = 0; i + 3 <= s.length; i += 3) p += CODE[s.slice(i, i + 3)] || 'X'; return p; };

const FRAMES = { f1: 'Forward frame 1', f2: 'Forward frame 2', f3: 'Forward frame 3', r1: 'Reverse frame 1', r2: 'Reverse frame 2', r3: 'Reverse frame 3' };
function frameSeq(seq, f) {
  const k = +f[1] - 1;
  return (f[0] === 'f' ? seq : revcomp(seq)).slice(k);
}

// Glocal alignment with affine gaps (Gotoh): `b` is aligned end-to-end,
// `a` (the reference) may have free unaligned ends.
function glocal(a, b, score, gapOpen, gapExt) {
  const n = a.length, m = b.length, W = m + 1, NEG = -1e9;
  const M = new Float64Array((n + 1) * W).fill(NEG), X = new Float64Array((n + 1) * W).fill(NEG), Y = new Float64Array((n + 1) * W).fill(NEG);
  const tM = new Uint8Array((n + 1) * W), tX = new Uint8Array((n + 1) * W), tY = new Uint8Array((n + 1) * W);
  for (let i = 0; i <= n; i++) M[i * W] = 0;
  for (let i = 0; i <= n; i++) {
    for (let j = 1; j <= m; j++) {
      const c = i * W + j;
      // Y: gap in a (consume b[j-1])
      const yo = M[c - 1] + gapOpen, ye = Y[c - 1] + gapExt;
      if (yo >= ye) { Y[c] = yo; tY[c] = 0; } else { Y[c] = ye; tY[c] = 2; }
      if (i === 0) continue;
      const d = c - W - 1, s = score(a[i - 1], b[j - 1]);
      let best = M[d], from = 0;
      if (X[d] > best) { best = X[d]; from = 1; }
      if (Y[d] > best) { best = Y[d]; from = 2; }
      M[c] = best + s; tM[c] = from;
      // X: gap in b (consume a[i-1])
      const xo = M[c - W] + gapOpen, xe = X[c - W] + gapExt;
      if (xo >= xe) { X[c] = xo; tX[c] = 0; } else { X[c] = xe; tX[c] = 1; }
    }
  }
  let bi = 0, bs = NEG, st = 0;
  for (let i = 0; i <= n; i++) {
    const c = i * W + m;
    if (M[c] > bs) { bs = M[c]; bi = i; st = 0; }
    if (Y[c] > bs) { bs = Y[c]; bi = i; st = 2; }
  }
  let i = bi, j = m, ra = [], rb = [];
  while (j > 0) {
    const c = i * W + j;
    if (st === 0) { ra.push(a[i - 1]); rb.push(b[j - 1]); st = tM[c]; i--; j--; }
    else if (st === 1) { ra.push(a[i - 1]); rb.push('-'); st = tX[c]; i--; }
    else { ra.push('-'); rb.push(b[j - 1]); st = tY[c]; j--; }
  }
  return { aStart: i, aEnd: bi, alnA: ra.reverse().join(''), alnB: rb.reverse().join(''), score: bs };
}

function copyButton(text) {
  const b = document.createElement('button');
  b.className = 'ghost'; b.textContent = 'Copy';
  b.onclick = () => navigator.clipboard.writeText(text).then(() => { b.textContent = 'Copied'; setTimeout(() => (b.textContent = 'Copy'), 1200); });
  return b;
}

// ---------------------------------------------------------------------------
// Tool 1 – sequence annotation (k-mer seeding + glocal alignment against the locus)
// ---------------------------------------------------------------------------
let GENOME = null, KMERS = null;
const K = 11;

async function loadGenome() {
  if (GENOME) return GENOME;
  GENOME = await (await fetch('genome.json')).json();
  KMERS = new Map();
  const g = GENOME.seq;
  for (let i = 0; i + K <= g.length; i++) {
    const k = g.slice(i, i + K);
    if (k.includes('N')) continue;
    const l = KMERS.get(k); l ? l.push(i) : KMERS.set(k, [i]);
  }
  return GENOME;
}

function bestDiagonal(q) {
  const votes = new Map();
  for (let i = 0; i + K <= q.length; i++) {
    const hits = KMERS.get(q.slice(i, i + K));
    if (hits) for (const p of hits) votes.set(p - i, (votes.get(p - i) || 0) + 1);
  }
  let d = null, v = 0;
  for (const [k, c] of votes) if (c > v) { v = c; d = k; }
  return { d, v };
}

async function runTool1() {
  const out = $('t1-out');
  const c = cleanDNA($('t1-in').value);
  if (c.error) return (out.innerHTML = err(c.error));
  const query = c.seq;
  if (query.length < 20) return (out.innerHTML = err('Sequence too short: the minimum length is 20 bp.'));
  if (query.length > 5000) return (out.innerHTML = err('Sequence too long: the maximum length is 5,000 bp.'));
  out.innerHTML = '<p class="ok">Searching…</p>';
  const G = await loadGenome();

  const fw = bestDiagonal(query), rv = bestDiagonal(revcomp(query));
  const strand = fw.v >= rv.v ? '+' : '-';
  const hit = strand === '+' ? fw : rv;
  if (hit.v < 3) return (out.innerHTML = err('No match found in the indexed region (HFE locus, GRCh38).'));
  const q = strand === '+' ? query : revcomp(query);
  const ws = Math.max(0, hit.d - 60), we = Math.min(G.seq.length, hit.d + q.length + 60);
  const aln = glocal(G.seq.slice(ws, we), q, (x, y) => (x === y && x !== 'N' ? 2 : -3), -6, -2);

  // Map each query base to a genome coordinate (1-based)
  const coord = new Array(q.length).fill(null);
  let gi = ws + aln.aStart, qi = 0, matches = 0, mism = 0, gaps = 0;
  for (let k = 0; k < aln.alnA.length; k++) {
    const x = aln.alnA[k], y = aln.alnB[k];
    if (x !== '-' && y !== '-') { coord[qi] = G.start + gi; x === y ? matches++ : mism++; gi++; qi++; }
    else if (x === '-') { coord[qi] = G.start + gi - 1; gaps++; qi++; }
    else { gi++; gaps++; }
  }
  const identity = matches / aln.alnA.length;
  if (identity < 0.9) return (out.innerHTML = err(`No good match found (best identity ${(identity * 100).toFixed(1)}%).`));
  const gStart = G.start + ws + aln.aStart, gEnd = G.start + ws + aln.aEnd - 1;

  // Pick the gene whose MANE Select exons overlap most of the amplicon
  const inExon = (ex, p) => ex.some(([s, e]) => p >= s && p <= e);
  let best = null, bestN = 0;
  for (const gene of G.genes) {
    const n = coord.filter(p => p !== null && inExon(gene.mane.exons, p)).length;
    if (n > bestN) { bestN = n; best = gene; }
  }
  const overlapping = G.genes.filter(g => g.start <= gEnd && g.end >= gStart);

  let html = `<dl class="kv">
    <dt>Location</dt><dd>${G.chrom}:${fmt(gStart)}-${fmt(gEnd)} (${strand} strand)</dd>
    <dt>Identity</dt><dd>${(identity * 100).toFixed(2)}% (${matches}/${aln.alnA.length} bp${mism ? `, ${mism} mismatches` : ''}${gaps ? `, ${gaps} gap positions` : ''})</dd>`;
  if (!best) {
    html += `<dt>Gene</dt><dd>${overlapping.length ? `${overlapping.map(g => g.name).join(', ')} (intronic only, no MANE Select exon overlap)` : 'Intergenic'}</dd></dl>`;
    return (out.innerHTML = html);
  }
  const ex = best.mane.exons, nEx = ex.length;
  const exNum = idx => (best.strand === '+' ? idx + 1 : nEx - idx);
  const hitEx = ex.map((e, i) => ({ e, n: exNum(i) })).filter(({ e }) => e[0] <= gEnd && e[1] >= gStart);
  html += `<dt>Gene</dt><dd><b>${esc(best.name)}</b> · ${esc(best.id)} · ${esc(best.type.replace(/_/g, ' '))}</dd>
    <dt>Gene location</dt><dd>${G.chrom}:${fmt(best.start)}-${fmt(best.end)} (${best.strand})</dd>
    <dt>MANE Select</dt><dd>${esc(best.mane.id)} (${esc(best.mane.name)})</dd>
    <dt>Exons covered</dt><dd>${hitEx.map(({ e, n }) => `exon ${n} of ${nEx} (${G.chrom}:${fmt(e[0])}-${fmt(e[1])})`).join('; ')}</dd></dl>`;

  // Format in the orientation the user pasted it
  const orient = strand === '+' ? i => i : i => q.length - 1 - i;
  const isEx = new Array(query.length);
  coord.forEach((p, i) => (isEx[orient(i)] = p !== null && inExon(ex, p)));
  let plain = '', rich = '';
  for (let i = 0; i < query.length; i += 60) {
    const num = `<span class="ln" data-n="${String(i + 1).padStart(5)}  "></span>`;
    let pl = '', ri = '';
    for (let j = i; j < Math.min(i + 60, query.length); j++) {
      const b = isEx[j] ? query[j] : query[j].toLowerCase();
      pl += b; ri += isEx[j] ? `<span class="ex">${b}</span>` : `<span class="in">${b}</span>`;
    }
    plain += pl; rich += num + ri + '\n';
  }
  html += `<p><b>Formatted sequence</b> (UPPERCASE = exon, lowercase = intron)</p><pre>${rich}</pre>`;
  out.innerHTML = html;
  out.appendChild(copyButton(plain));
}

// ---------------------------------------------------------------------------
// Tool 2 – six-frame translation
// ---------------------------------------------------------------------------
function runTool2() {
  const out = $('t2-out');
  const c = cleanDNA($('t2-in').value);
  if (c.error) return (out.innerHTML = err(c.error));
  let html = '';
  for (const f of Object.keys(FRAMES)) {
    const p = translate(frameSeq(c.seq, f));
    let lines = '';
    for (let i = 0; i < p.length; i += 60) lines += p.slice(i, i + 60).replace(/\*/g, '<span class="stop">*</span>') + '\n';
    html += `<div class="frame"><b>${FRAMES[f]}</b><pre>${lines || '(too short)'}</pre></div>`;
  }
  out.innerHTML = html;
}

// ---------------------------------------------------------------------------
// Tool 3 – align peptide to UniProt reference
// ---------------------------------------------------------------------------
const uniprotCache = {};
async function fetchUniprot(acc) {
  acc = acc.trim().toUpperCase();
  if (!acc) throw new Error('Please enter a UniProt accession code.');
  if (!/^[A-Z0-9]+(-\d+)?$/.test(acc)) throw new Error(`"${acc}" does not look like a UniProt accession code.`);
  if (uniprotCache[acc]) return uniprotCache[acc];
  const r = await fetch(`https://rest.uniprot.org/uniprotkb/${acc}.fasta`);
  const t = r.ok ? await r.text() : '';
  if (!t.startsWith('>')) throw new Error(`Could not retrieve UniProt entry ${acc}.`);
  const lines = t.trim().split('\n');
  return (uniprotCache[acc] = { acc, header: lines[0].slice(1), seq: lines.slice(1).join('').trim() });
}

const protScore = (x, y) => (x === y ? 2 : -1);

async function runTool3() {
  const out = $('t3-out');
  const c = cleanProtein($('t3-in').value);
  if (c.error) return (out.innerHTML = err(c.error));
  let ref;
  try { ref = await fetchUniprot($('uniprot').value); } catch (e) { return (out.innerHTML = err(e.message)); }
  const aln = glocal(ref.seq, c.seq, protScore, -5, -0.5);
  const hl = parseInt($('t3-pos').value, 10);

  // Extend alignment with the unaligned reference ends so the full protein is shown
  const A = ref.seq.slice(0, aln.aStart) + aln.alnA + ref.seq.slice(aln.aEnd);
  const B = '-'.repeat(aln.aStart) + aln.alnB + '-'.repeat(ref.seq.length - aln.aEnd);
  let ident = 0, cols = 0;
  for (let k = 0; k < aln.alnA.length; k++) { cols++; if (aln.alnA[k] === aln.alnB[k]) ident++; }

  const L = 50;
  let txt = '', refPos = 0;
  for (let i = 0; i < A.length; i += L) {
    const first = refPos + 1;
    let ra = '', mm = '', rb = '';
    for (let k = i; k < Math.min(i + L, A.length); k++) {
      const x = A[k], y = B[k];
      if (x !== '-') refPos++;
      const on = x !== '-' && refPos === hl;
      ra += on ? `<mark>${x}</mark>` : x;
      rb += on ? `<mark>${y}</mark>` : y;
      mm += x === y && x !== '-' ? '|' : (x !== '-' && y !== '-' ? '.' : ' ');
    }
    txt += `<span class="ln" data-n="Reference [${String(first).padStart(4)}] "></span>${ra}\n<span class="ln" data-n="${' '.repeat(17)}"></span>${mm}\n<span class="ln" data-n="Yours${' '.repeat(12)}"></span>${rb}\n\n`;
  }
  out.innerHTML = `<dl class="kv">
      <dt>Reference</dt><dd>${esc(ref.header)} (${ref.seq.length} aa)</dd>
      <dt>Your sequence</dt><dd>aligns to reference residues ${aln.aStart + 1}–${aln.aEnd}; identity ${ident}/${cols} (${(100 * ident / cols).toFixed(1)}%)</dd>
    </dl>
    ${c.hadStop ? '<p class="err">Your sequence contained stop codons (*). Are you sure you translated the right frame / region?</p>' : ''}
    ${hl && hl > ref.seq.length ? `<p class="err">Position ${hl} is beyond the end of the protein (${ref.seq.length} aa).</p>` : ''}
    <pre>${txt}</pre>
    <p class="hint"><code>|</code> identical · <code>.</code> different · <code>-</code> gap</p>`;
}

// keep the two UniProt fields in sync
['uniprot', 'uniprot2'].forEach((id, k, ids) => $(id).addEventListener('input', e => ($(ids[1 - k]).value = e.target.value)));

// ---------------------------------------------------------------------------
// Tools 4 & 5 – structure
// ---------------------------------------------------------------------------
const AA3 = { ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M',
  PHE: 'F', PRO: 'P', SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V', MSE: 'M', SEC: 'U', PYL: 'O' };
const PALETTE = ['#4e79a7', '#f28e2b', '#59a14f', '#b07aa1', '#76b7b2', '#edc948', '#9c755f', '#bab0ac', '#ff9da7'];
let S = null; // loaded structure state

async function loadStructure() {
  const msg = $('t4-msg');
  const pdb = $('pdb').value.trim().toUpperCase();
  if (!/^[0-9][A-Z0-9]{3}$/.test(pdb)) return (msg.innerHTML = err('Please enter a valid 4-character PDB code.'));
  let ref;
  try { ref = await fetchUniprot($('uniprot2').value); } catch (e) { return (msg.innerHTML = err(e.message)); }
  if (S && S.pdb === pdb && S.ref.acc === ref.acc) return S;
  msg.innerHTML = `<p class="ok">Downloading ${pdb}…</p>`;
  const r = await fetch(`https://files.rcsb.org/download/${pdb}.cif`);
  if (!r.ok) return (msg.innerHTML = err(`Could not download PDB entry ${pdb}.`));
  const cif = await r.text();

  // entity names (optional, for the legend)
  const entityOf = {};
  try {
    const qry = `{entry(entry_id:"${pdb}"){polymer_entities{rcsb_polymer_entity{pdbx_description}rcsb_polymer_entity_container_identifiers{auth_asym_ids}}}}`;
    const j = await (await fetch('https://data.rcsb.org/graphql?query=' + encodeURIComponent(qry))).json();
    j.data.entry.polymer_entities.forEach((e, k) => e.rcsb_polymer_entity_container_identifiers.auth_asym_ids.forEach(ch =>
      (entityOf[ch] = { k, name: e.rcsb_polymer_entity.pdbx_description })));
  } catch (e) { /* legend falls back to chain ids */ }

  $('t4-view').innerHTML = '';
  $('t4-view').classList.remove('hidden');
  const viewer = $3Dmol.createViewer($('t4-view'), { backgroundColor: 'white' });
  const model = viewer.addModel(cif, 'cif');
  const atoms = model.selectedAtoms({});

  // residues per chain (protein only, first alt-loc)
  const chains = {};
  for (const a of atoms) {
    if (a.hetflag && a.resn !== 'MSE') continue;
    const aa = AA3[a.resn];
    if (!aa) continue;
    const ch = chains[a.chain] || (chains[a.chain] = { id: a.chain, res: [], byResi: {} });
    let res = ch.byResi[a.resi];
    if (!res) { res = ch.byResi[a.resi] = { resi: a.resi, aa, atoms: {} }; ch.res.push(res); }
    if (!res.atoms[a.atom]) res.atoms[a.atom] = a;
  }
  // map every chain to UniProt numbering
  const mapped = [];
  for (const ch of Object.values(chains)) {
    ch.res = ch.res.filter(r => r.atoms.CA);
    ch.seq = ch.res.map(r => r.aa).join('');
    ch.entity = entityOf[ch.id];
    if (ch.seq.length < 20) continue;
    const aln = glocal(ref.seq, ch.seq, protScore, -5, -0.5);
    let ui = aln.aStart, ri = 0, same = 0, pairs = 0;
    const map = {};
    for (let k = 0; k < aln.alnA.length; k++) {
      const x = aln.alnA[k], y = aln.alnB[k];
      if (x !== '-') ui++;
      if (y !== '-') {
        if (x !== '-') { pairs++; if (x === y) same++; ch.res[ri].uni = ui; map[ui] = ch.res[ri]; }
        ri++;
      }
    }
    if (pairs >= 20 && same / pairs >= 0.9) { ch.uniMap = map; mapped.push(ch); }
    else ch.res.forEach(r => delete r.uni);
  }
  if (!mapped.length) return (msg.innerHTML = err(`None of the chains in ${pdb} match UniProt ${ref.acc}. Check both codes.`));

  // colours by entity (or chain)
  const ids = Object.keys(chains).sort();
  const colorKey = ch => (chains[ch].entity ? 'e' + chains[ch].entity.k : ch);
  const keys = [...new Set(ids.map(colorKey))];
  const colorOf = {};
  ids.forEach(ch => (colorOf[ch] = PALETTE[keys.indexOf(colorKey(ch)) % PALETTE.length]));
  const groups = {};
  ids.forEach(ch => { const k = colorKey(ch); (groups[k] = groups[k] || { name: chains[ch].entity ? chains[ch].entity.name : 'Chain', chains: [] }).chains.push(ch); });
  $('t4-legend').innerHTML = Object.entries(groups).map(([k, g]) =>
    `<span style="--c:${colorOf[g.chains[0]]}">${esc(g.name.toLowerCase().replace(/^\w/, c => c.toUpperCase()))} (chain ${g.chains.join(', ')})</span>`).join('') +
    '<span style="--c:red">your subsequence</span>';

  S = { pdb, ref, cif, viewer, chains, mapped, colorOf };
  resetStyle();
  viewer.zoomTo(); viewer.render();
  const m0 = mapped[0], r0 = m0.res.find(r => r.uni);
  msg.innerHTML = `<p class="ok">Loaded ${pdb}: ${ids.length} protein chains. Chain${mapped.length > 1 ? 's' : ''} ${mapped.map(c => c.id).join(', ')} match ${ref.acc}
    (UniProt residues ${m0.res.find(r => r.uni).uni}–${[...m0.res].reverse().find(r => r.uni).uni} are resolved in chain ${m0.id}).
    Numbering: UniProt ${r0.uni} = PDB ${r0.resi} in chain ${m0.id}.</p>`;
  buildDisulfideTable();
  return S;
}

function resetStyle() {
  S.viewer.removeAllLabels(); S.viewer.removeAllShapes();
  S.viewer.setStyle({}, {});
  for (const [ch, col] of Object.entries(S.colorOf)) S.viewer.setStyle({ chain: ch }, { cartoon: { color: col } });
}

async function runTool4() {
  const msg = $('t4-msg');
  if (!(await loadStructure())) return;
  const sub = $('t4-sub').value.replace(/[\s\d]/g, '').toUpperCase();
  if (!sub) return;
  if (!/^[A-Z]+$/.test(sub)) return (msg.innerHTML += err('The subsequence must contain amino acid letters only.'));

  // look in mapped chains first (UniProt numbering), then any chain; fall back to the UniProt sequence
  let hits = [];
  const order = [...S.mapped, ...Object.values(S.chains).filter(c => !c.uniMap)];
  for (const ch of order) {
    let from = ch.seq.indexOf(sub);
    while (from !== -1) { hits.push({ ch, res: ch.res.slice(from, from + sub.length) }); from = ch.seq.indexOf(sub, from + 1); }
  }
  if (!hits.length) {
    const u = S.ref.seq.indexOf(sub);
    if (u !== -1) {
      const ch = S.mapped[0], res = [];
      for (let n = u + 1; n <= u + sub.length; n++) if (ch.uniMap[n]) res.push(ch.uniMap[n]);
      if (!res.length) return (msg.innerHTML += err(`"${sub}" is in the UniProt sequence (residues ${u + 1}–${u + sub.length}) but those residues are not resolved in the structure.`));
      hits.push({ ch, res, partial: res.length < sub.length });
    }
  }
  resetStyle();
  if (!hits.length) { S.viewer.render(); return (msg.innerHTML = err(`"${sub}" was not found in ${S.pdb} or in ${S.ref.acc}.`)); }

  const { ch, res, partial } = hits[0];
  const resis = res.map(r => r.resi);
  S.viewer.setStyle({ chain: ch.id, resi: resis }, { cartoon: { color: 'red' }, stick: {} });
  for (const r of res) {
    S.viewer.addLabel(`${r.aa}${r.uni ?? r.resi}`, { fontColor: 'black', backgroundColor: '#ffe08a', backgroundOpacity: 0.9, fontSize: 11, inFront: true },
      { chain: ch.id, resi: r.resi, atom: 'CA' });
  }
  S.viewer.zoomTo({ chain: ch.id, resi: resis }); S.viewer.zoom(0.8); S.viewer.render();
  const numb = ch.uniMap ? `UniProt residues ${res[0].uni}–${res[res.length - 1].uni}` : `PDB residues ${res[0].resi}–${res[res.length - 1].resi} (this chain is not ${S.ref.acc})`;
  const others = hits.slice(1).map(h => h.ch.id);
  msg.innerHTML = `<p class="ok">Found in chain ${ch.id}: ${numb} (PDB numbering ${res[0].resi}–${res[res.length - 1].resi})${partial ? ' – only partly resolved in the structure' : ''}.
    ${others.length ? `Also present in chain${others.length > 1 ? 's' : ''} ${[...new Set(others)].join(', ')}.` : ''}</p>`;
}

// ---- Tool 5: disulfide bridges ----
let SSviewer = null;
function buildDisulfideTable() {
  const out = $('t5-out');
  const sg = [];
  for (const ch of Object.values(S.chains)) for (const r of ch.res) if (r.aa === 'C' && r.atoms.SG) sg.push({ ch, r, a: r.atoms.SG });
  const dist = (p, q) => Math.hypot(p.x - q.x, p.y - q.y, p.z - q.z);
  const label = (ch, r) => `C${r.uni ?? r.resi}${ch.uniMap ? '' : ' (PDB no.)'} chain ${ch.id}`;

  const ch = S.mapped[0];
  const rows = [];
  [...S.ref.seq].forEach((aa, i) => {
    if (aa !== 'C') return;
    const n = i + 1, r = ch.uniMap[n];
    if (!r) return rows.push(`<tr><td>C${n}</td><td class="ok">not resolved in the structure</td><td></td><td></td><td></td></tr>`);
    if (!r.atoms.SG) return rows.push(`<tr><td>C${n}</td><td>${r.resi}</td><td class="ok">no SG atom</td><td></td><td></td></tr>`);
    let best = null;
    for (const o of sg) if (o.r !== r) { const d = dist(r.atoms.SG, o.a); if (d <= 2.5 && (!best || d < best.d)) best = { ...o, d }; }
    rows.push(`<tr><td>C${n}</td><td>${r.resi}</td>
      <td>${best ? label(best.ch, best.r) : '<span class="ok">none (free thiol)</span>'}</td>
      <td>${best ? best.d.toFixed(2) + ' Å' : ''}</td>
      <td>${best ? `<button class="ghost" data-a="${ch.id}:${r.resi}" data-b="${best.ch.id}:${best.r.resi}">Show</button>` : ''}</td></tr>`);
  });
  out.innerHTML = `<p class="ok">Cysteines of ${S.ref.acc} in chain ${ch.id} of ${S.pdb} (UniProt numbering).</p>
    <table><thead><tr><th>Cysteine</th><th>PDB no.</th><th>Bonded to</th><th>S–S distance</th><th></th></tr></thead><tbody>${rows.join('')}</tbody></table>`;
  out.querySelectorAll('button[data-a]').forEach(b => (b.onclick = () => showBridge(b.dataset.a, b.dataset.b)));
  $('t5-view').classList.add('hidden');
}

function showBridge(a, b) {
  const [ca, ra] = a.split(':'), [cb, rb] = b.split(':');
  const A = S.chains[ca].byResi[ra], B = S.chains[cb].byResi[rb];
  const el = $('t5-view');
  el.classList.remove('hidden');
  if (!SSviewer || SSviewer.pdb !== S.pdb) {
    el.innerHTML = '';
    SSviewer = $3Dmol.createViewer(el, { backgroundColor: 'white' });
    SSviewer.addModel(S.cif, 'cif'); SSviewer.pdb = S.pdb;
  }
  const v = SSviewer;
  v.removeAllLabels(); v.removeAllShapes();
  v.setStyle({}, { cartoon: { color: '#9a9a9a', opacity: 0.55 } });
  const sel = [{ chain: ca, resi: A.resi }, { chain: cb, resi: B.resi }];
  sel.forEach(s => v.setStyle(s, { cartoon: { color: '#9a9a9a', opacity: 0.55 }, stick: { colorscheme: 'yellowCarbon', radius: 0.25 } }));
  [[ca, A], [cb, B]].forEach(([c, r]) => v.addLabel(`C${r.uni ?? r.resi} (chain ${c})`,
    { fontColor: 'black', backgroundColor: 'white', backgroundOpacity: 0.85, fontSize: 13, inFront: true }, { chain: c, resi: r.resi, atom: 'CA' }));
  const p = A.atoms.SG, q = B.atoms.SG;
  v.addCylinder({ start: { x: p.x, y: p.y, z: p.z }, end: { x: q.x, y: q.y, z: q.z }, radius: 0.18, color: 'yellow', fromCap: 1, toCap: 1 });
  v.zoomTo({ or: sel }); v.zoom(0.8); v.render();
  el.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
}

// ---------------------------------------------------------------------------
// Tool 6 – codons to amino acids
// ---------------------------------------------------------------------------
function runTool6() {
  const out = $('t6-out');
  const c = cleanDNA($('t6-in').value);
  if (c.error) return (out.innerHTML = err(c.error));
  const f = $('t6-frame').value, s = frameSeq(c.seq, f), k = +f[1] - 1, L = c.seq.length;
  let html = '';
  for (let i = 0; i + 3 <= s.length; i += 3) {
    const cod = s.slice(i, i + 3), aa = CODE[cod] || 'X';
    // positions of the codon in the sequence as pasted (1-based)
    const pos = f[0] === 'f' ? `${k + i + 1}–${k + i + 3}` : `${L - (k + i + 2)}–${L - (k + i)}`;
    html += `<div class="cd${aa === '*' ? ' s' : ''}"><span>${cod}</span><span class="aa">${aa}</span><span class="pos">${pos}</span></div>`;
  }
  out.innerHTML = `<p class="hint">${FRAMES[f]}. Numbers under each codon: its nucleotide positions in the sequence you pasted${f[0] === 'r' ? ' (codons are read on the reverse-complement strand)' : ''}.</p><div class="codons">${html}</div>`;
}

// ---------------------------------------------------------------------------
// Tool 7 – chromatogram viewer (ABIF parser + canvas)
// ---------------------------------------------------------------------------
function parseABIF(buf) {
  const dv = new DataView(buf);
  const str = (o, n) => String.fromCharCode(...new Uint8Array(buf, o, n));
  if (str(0, 4) !== 'ABIF') throw new Error('This is not an ABIF (.ab1) file.');
  const n = dv.getInt32(18), off = dv.getInt32(26), tags = {};
  for (let i = 0; i < n; i++) {
    const e = off + i * 28;
    const name = str(e, 4) + dv.getInt32(e + 4);
    const size = dv.getInt32(e + 16), cnt = dv.getInt32(e + 12);
    tags[name] = { type: dv.getInt16(e + 8), cnt, size, off: size <= 4 ? e + 20 : dv.getInt32(e + 20) };
  }
  const chars = t => (t ? str(t.off, t.cnt) : '');
  const shorts = t => { const a = new Int16Array(t.cnt); for (let i = 0; i < t.cnt; i++) a[i] = dv.getInt16(t.off + 2 * i); return a; };
  const order = chars(tags.FWO_1) || 'GATC';
  const traces = {};
  [...order].forEach((b, i) => (traces[b] = shorts(tags['DATA' + (9 + i)])));
  const bases = chars(tags.PBAS2 || tags.PBAS1).toUpperCase();
  const ploc = shorts(tags.PLOC2 || tags.PLOC1);
  let name = '';
  if (tags.SMPL1) { const t = tags.SMPL1; name = str(t.off + 1, new Uint8Array(buf, t.off, 1)[0]); }
  return { traces, bases, ploc, name, len: traces[order[0]].length };
}

const IUPAC = { A: 'A', C: 'C', G: 'G', T: 'T', R: 'AG', Y: 'CT', S: 'CG', W: 'AT', K: 'GT', M: 'AC', B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT' };
const compatible = (call, q) => (IUPAC[call] || 'ACGT').includes(q);

let CH = null;
function drawChromo() {
  if (!CH) return;
  const box = $('t7-scroll'), cv = $('t7-canvas'), dpr = window.devicePixelRatio || 1;
  const w = box.clientWidth, h = 260, zoom = +$('t7-zoom').value, gain = +$('t7-gain').value;
  cv.style.width = w + 'px'; cv.width = w * dpr; cv.height = h * dpr;
  const ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const css = getComputedStyle(document.documentElement), col = b => css.getPropertyValue('--' + b).trim();
  ctx.clearRect(0, 0, w, h);
  const x0 = box.scrollLeft / zoom, x1 = x0 + w / zoom, top = 38, base = h - 10, yMax = CH.yMax / gain;
  const X = i => (i - x0) * zoom;

  if (CH.hit) {
    const a = CH.ploc[CH.hit.start], b = CH.ploc[CH.hit.end - 1];
    ctx.fillStyle = col('mark'); ctx.globalAlpha = 0.6; ctx.fillRect(X(a) - 6 * zoom, 0, X(b) - X(a) + 12 * zoom, h); ctx.globalAlpha = 1;
  }
  for (const b of 'ACGT') {
    const t = CH.traces[b]; if (!t) continue;
    ctx.strokeStyle = col(b); ctx.lineWidth = 1.3; ctx.beginPath();
    const s = Math.max(0, Math.floor(x0) - 1), e = Math.min(t.length - 1, Math.ceil(x1) + 1);
    for (let i = s; i <= e; i++) {
      const y = base - Math.min(1, t[i] / yMax) * (base - top);
      i === s ? ctx.moveTo(X(i), y) : ctx.lineTo(X(i), y);
    }
    ctx.stroke();
  }
  ctx.textAlign = 'center'; ctx.font = 'bold 13px ui-monospace, Menlo, monospace';
  for (let k = 0; k < CH.bases.length; k++) {
    const p = CH.ploc[k]; if (p < x0 - 10 || p > x1 + 10) continue;
    const c = CH.bases[k];
    ctx.fillStyle = 'ACGT'.includes(c) ? col(c) : col('fg');
    ctx.fillText(c, X(p), 30);
    if ((k + 1) % 10 === 0) { ctx.fillStyle = col('muted'); ctx.font = '10px ui-monospace, Menlo, monospace'; ctx.fillText(k + 1, X(p), 12); ctx.font = 'bold 13px ui-monospace, Menlo, monospace'; }
  }
}

function setSpacer() { $('t7-spacer').style.width = CH.len * +$('t7-zoom').value + 'px'; }

async function loadChromo(file) {
  const msg = $('t7-msg');
  try { CH = parseABIF(await file.arrayBuffer()); } catch (e) { CH = null; $('t7-scroll').classList.add('hidden'); return (msg.innerHTML = err(e.message)); }
  const all = Object.values(CH.traces).flatMap(t => Array.from(t.subarray(CH.ploc[0] || 0, CH.ploc[CH.ploc.length - 1] || t.length))).sort((a, b) => a - b);
  CH.yMax = all[Math.floor(all.length * 0.995)] || 1;
  CH.hit = null;
  $('t7-scroll').classList.remove('hidden');
  setSpacer(); $('t7-scroll').scrollLeft = 0; drawChromo();
  msg.innerHTML = `<p class="ok">${esc(CH.name || file.name)} · ${CH.bases.length} bases called.</p>`;
  if ($('t7-q').value.trim()) findInChromo();
}

function findInChromo() {
  const msg = $('t7-msg');
  if (!CH) return (msg.innerHTML = err('Open an .ab1 file first.'));
  const c = cleanDNA($('t7-q').value);
  if (c.error) return (msg.innerHTML = err(c.error));
  const q = c.seq;
  if (q.length < 8) return (msg.innerHTML = err('Use at least 8 bases so the search is specific.'));
  // best window allowing mismatches; ambiguous calls (e.g. R = A/G) count as matches
  let best = { mm: Infinity };
  for (const [s, rc] of [[q, false], [revcomp(q), true]]) {
    for (let i = 0; i + s.length <= CH.bases.length; i++) {
      let mm = 0;
      for (let j = 0; j < s.length && mm < best.mm; j++) if (!compatible(CH.bases[i + j], s[j])) mm++;
      if (mm < best.mm) best = { mm, start: i, end: i + s.length, rc };
    }
  }
  const maxMM = Math.max(1, Math.floor(q.length / 10));
  if (best.mm > maxMM) { CH.hit = null; drawChromo(); return (msg.innerHTML = err(`Sequence not found in this read (best match has ${best.mm} mismatches).`)); }
  CH.hit = best;
  const zoom = +$('t7-zoom').value, box = $('t7-scroll');
  box.scrollLeft = (CH.ploc[best.start] + CH.ploc[best.end - 1]) / 2 * zoom - box.clientWidth / 2;
  drawChromo();
  const called = CH.bases.slice(best.start, best.end);
  const amb = [...called].map((b, j) => ('ACGT'.includes(b) ? null : `${b} at base ${best.start + j + 1}`)).filter(Boolean);
  msg.innerHTML = `<p class="ok">Found at bases ${best.start + 1}–${best.end}${best.rc ? ' as the <b>reverse complement</b> (this read is on the other strand)' : ''}${best.mm ? `, ${best.mm} mismatch${best.mm > 1 ? 'es' : ''}` : ''}.
    Called: <code>${esc(called)}</code>${amb.length ? ` · ambiguous calls: ${esc(amb.join(', '))}` : ''}</p>`;
}

$('t7-file').addEventListener('change', e => e.target.files[0] && loadChromo(e.target.files[0]));
$('t7-scroll').addEventListener('scroll', () => requestAnimationFrame(drawChromo));
$('t7-zoom').addEventListener('input', () => {
  if (!CH) return;
  const box = $('t7-scroll'), old = box.scrollWidth, mid = (box.scrollLeft + box.clientWidth / 2) / old;
  setSpacer(); box.scrollLeft = mid * box.scrollWidth - box.clientWidth / 2; drawChromo();
});
$('t7-gain').addEventListener('input', drawChromo);
window.addEventListener('resize', drawChromo);
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', drawChromo);

// ---------------------------------------------------------------------------
// Wiring
// ---------------------------------------------------------------------------
const guard = (fn, outId) => async () => {
  try { await fn(); } catch (e) { console.error(e); $(outId).innerHTML = err('Something went wrong: ' + e.message); }
};
$('t1-go').onclick = guard(runTool1, 't1-out');
$('t2-go').onclick = guard(runTool2, 't2-out');
$('t3-go').onclick = guard(runTool3, 't3-out');
$('t3-pos').addEventListener('change', () => $('t3-out').innerHTML && guard(runTool3, 't3-out')());
$('t4-load').onclick = guard(loadStructure, 't4-msg');
$('t4-go').onclick = guard(runTool4, 't4-msg');
$('t4-sub').addEventListener('keydown', e => e.key === 'Enter' && $('t4-go').click());
$('t6-go').onclick = guard(runTool6, 't6-out');
$('t7-go').onclick = guard(async () => findInChromo(), 't7-msg');
$('t7-q').addEventListener('keydown', e => e.key === 'Enter' && $('t7-go').click());
loadGenome().catch(() => {});
