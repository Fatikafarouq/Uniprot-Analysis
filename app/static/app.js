const $ = (id) => document.getElementById(id);
const statusBox = $("status");
const choices = $("choices");
const results = $("results");

function escapeHtml(value) {
  return String(value ?? "").replace(/[&<>'"]/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;',"'":'&#39;','"':'&quot;'}[ch]));
}

function showStatus(message, warning=false) {
  statusBox.hidden = !message;
  statusBox.className = `status ${warning ? 'warning' : ''}`;
  statusBox.textContent = message || "";
}

function renderWarnings(warnings=[]) {
  if (!warnings.length) return "";
  return `<div class="status warning"><strong>Source note:</strong><ul>${warnings.map(w => `<li>${escapeHtml(w.message)}</li>`).join('')}</ul></div>`;
}

async function post(url, body) {
  const response = await fetch(url, {method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify(body)});
  if (!response.ok) throw new Error(`Request failed (${response.status})`);
  return response.json();
}

function resetMain() {
  choices.hidden = true; choices.innerHTML = "";
  results.hidden = true; results.innerHTML = "";
}

async function runLookup(query, extra={}) {
  resetMain();
  showStatus("Checking UniProt…");
  try {
    const data = await post('/api/lookup', {query, ...extra});
    showStatus("");
    renderLookup(data, query);
  } catch (err) {
    showStatus(err.message, true);
  }
}

function renderLookup(data, originalQuery) {
  if (data.status === 'ready') return renderReady(data);

  if (data.status === 'needs_organism_choice') {
    choices.hidden = false;
    choices.innerHTML = `<h2>Which organism did you mean?</h2>${renderWarnings(data.warnings)}` +
      data.options.map(opt => `<div class="choice"><div><strong>${escapeHtml(opt.label)}</strong><div class="small">${escapeHtml(opt.rank)} · UniProt taxonomy ${escapeHtml(opt.taxon_id)}</div></div><button data-taxon="${opt.taxon_id}" data-name="${escapeHtml(opt.label)}" data-phrase="${escapeHtml(opt.input_phrase || '')}">Use this organism</button></div>`).join('') +
      `<div class="choice"><div><strong>None of these</strong><div class="small">Continue without choosing an organism.</div></div><button data-broad="1">Continue broadly</button></div>`;
    choices.querySelectorAll('[data-taxon]').forEach(btn => btn.addEventListener('click', () => runLookup(originalQuery, {taxon_id:Number(btn.dataset.taxon), organism_name:btn.dataset.name, organism_phrase:btn.dataset.phrase, resolve_organism:false})));
    choices.querySelector('[data-broad]').addEventListener('click', () => runLookup(originalQuery, {resolve_organism:false}));
    return;
  }

  if (data.status === 'needs_protein_choice') {
    choices.hidden = false;
    choices.innerHTML = `<h2>Choose the direct match you meant</h2><p>${escapeHtml(data.message)}</p>${renderWarnings(data.warnings)}` +
      data.options.map(opt => `<div class="choice"><div><strong>${escapeHtml(opt.gene)}</strong><div>${escapeHtml(opt.organism)}${opt.scientific_name && opt.scientific_name !== opt.organism ? ` · <i>${escapeHtml(opt.scientific_name)}</i>` : ''}</div><div class="small">${escapeHtml(opt.record_count_in_search)} matching record(s) in this search</div></div><button data-gene="${escapeHtml(opt.gene)}" data-taxon="${opt.taxon_id}" data-species="${escapeHtml(opt.organism)}">Explain these records</button></div>`).join('');
    choices.querySelectorAll('[data-gene]').forEach(btn => btn.addEventListener('click', async () => {
      showStatus('Fetching the full UniProt record set…');
      const payload = await post('/api/explain', {gene:btn.dataset.gene, taxon_id:Number(btn.dataset.taxon), species_name:btn.dataset.species});
      showStatus(''); renderReady(payload);
    }));
    return;
  }

  if (data.status === 'no_direct_match') {
    results.hidden = false;
    results.innerHTML = `<div class="summary"><h2>No direct identity match</h2><p>${escapeHtml(data.message)}</p><p>You can use the discovery results below, but they are deliberately kept separate from direct gene/protein identity.</p></div>${renderWarnings(data.warnings)}${renderDiscoveryCards(data.discovery || [])}`;
    return;
  }

  showStatus(data.message || 'No result was returned.', data.status !== 'ready');
}

function renderReady(data) {
  results.hidden = false;
  const cards = (data.records || []).map(item => {
    const r = item.record;
    const review = r.review || {};
    const evidence = r.function_evidence || {};
    const links = r.links || {};
    const isoforms = r.isoforms || [];
    const pdbLinks = (links.pdb || []).map(x => `<a href="${x.url}" target="_blank" rel="noopener">PDB ${escapeHtml(x.id)}</a>`).join(' · ');
    return `<article class="card">
      <div><span class="badge ${escapeHtml(review.status)}">${escapeHtml(review.label)}</span></div>
      <h3>${escapeHtml(item.accession)}</h3>
      <div class="meta">${escapeHtml(r.name)} · ${escapeHtml(item.length ?? 'unknown')} aa</div>
      <ul>${(item.sentences || []).map(x => `<li>${escapeHtml(x)}</li>`).join('')}</ul>
      <div class="evidence-box"><strong>Function evidence:</strong> ${escapeHtml(evidence.summary || 'Not available')}</div>
      ${isoforms.length ? `<details><summary>Isoforms described in this UniProt entry (${isoforms.length})</summary><ul>${isoforms.map(i => `<li>${escapeHtml(i.ids?.join(', ') || i.name || 'Unnamed isoform')} — ${escapeHtml(i.sequence_status || 'status unavailable')}</li>`).join('')}</ul></details>` : ''}
      <details><summary>Database details</summary>
        <p><strong>Protein existence:</strong> ${escapeHtml(r.existence)}</p>
        <p><strong>Ensembl transcript:</strong> ${escapeHtml(r.transcript || 'No cross-reference returned by UniProt')}</p>
        ${r.appris ? `<p><strong>APPRIS:</strong> ${escapeHtml(r.appris)}</p>` : ''}
        ${r.canonical ? `<p><strong>Ensembl canonical:</strong> Yes</p>` : ''}
        ${r.gene_centric ? `<p><strong>UniProt gene-centric representative:</strong> Yes</p>` : ''}
      </details>
      <div class="links">
        ${links.uniprot ? `<a href="${links.uniprot}" target="_blank" rel="noopener">UniProt</a>` : ''}
        ${links.alphafold ? `<a href="${links.alphafold}" target="_blank" rel="noopener">AlphaFold DB</a>` : ''}
        ${pdbLinks || (links.pdb_search ? `<a href="${links.pdb_search}" target="_blank" rel="noopener">Search PDB</a>` : '')}
      </div>
    </article>`;
  }).join('');

  results.innerHTML = `${renderWarnings(data.warnings)}<section class="summary"><h2>UniProt has ${escapeHtml(data.record_count)} record${data.record_count === 1 ? '' : 's'} for ${escapeHtml(data.gene)} in ${escapeHtml(data.species)}.</h2><p>${escapeHtml(data.review_summary)}</p><p>${escapeHtml(data.review_explanation)}</p>${data.isoform_explanation ? `<p><strong>Isoforms:</strong> ${escapeHtml(data.isoform_explanation)}</p>` : ''}<p class="small">The cards below are comparisons, not rankings or recommendations.</p></section><div class="grid">${cards}</div>`;
}

function renderDiscoveryCards(items) {
  if (!items.length) return `<div class="summary"><p>No traceable discovery result was found.</p></div>`;
  return `<div class="grid">${items.map(item => `<article class="card"><h3>${escapeHtml(item.protein_name || item.gene || item.accession)}</h3><div class="meta">${escapeHtml(item.gene || 'No gene label')} · ${escapeHtml(item.organism || 'Unknown organism')} · ${escapeHtml(item.accession)}</div><strong>Why this appeared</strong>${(item.why || []).map(ctx => `<div class="evidence-box"><div class="small">${escapeHtml(ctx.source)}</div>${escapeHtml(ctx.text)}</div>`).join('')}</article>`).join('')}</div>`;
}

$('lookup-form').addEventListener('submit', (event) => {
  event.preventDefault();
  const query = $('query').value.trim();
  if (query) runLookup(query);
});

$('show-discovery').addEventListener('click', () => {
  const panel = $('discovery-panel');
  panel.hidden = !panel.hidden;
});

$('discovery-form').addEventListener('submit', async (event) => {
  event.preventDefault();
  const query = $('discovery-query').value.trim();
  if (!query) return;
  $('discovery-results').innerHTML = '<p>Checking UniProt…</p>';
  try {
    const data = await post('/api/discover', {query});
    $('discovery-results').innerHTML = renderWarnings(data.warnings) + renderDiscoveryCards(data.results || []);
  } catch (err) {
    $('discovery-results').innerHTML = `<div class="status warning">${escapeHtml(err.message)}</div>`;
  }
});
