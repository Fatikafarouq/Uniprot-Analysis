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
  let payload = null;
  try { payload = await response.json(); } catch (_) {}
  if (!response.ok) throw new Error(payload?.message || `Request failed (${response.status})`);
  return payload;
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

  if (data.status === 'needs_spelling_confirmation') {
    const suggestion = data.suggestion?.phrase || '';
    choices.hidden = false;
    choices.innerHTML = `<h2>Possible spelling correction</h2><p>${escapeHtml(data.message)}</p><p class="small">The suggestion comes from wording found in UniProt records for the resolved organism. It will not be used unless you confirm it.</p>${renderWarnings(data.warnings)}<div class="choice"><div><strong>${escapeHtml(suggestion)}</strong></div><div><button id="accept-spelling">Use this spelling</button> <button id="reject-spelling" class="secondary">Keep my wording</button></div></div>`;
    $('accept-spelling').addEventListener('click', () => runLookup(originalQuery, {confirmed_spelling:suggestion}));
    $('reject-spelling').addEventListener('click', () => runLookup(originalQuery, {confirmed_spelling:''}));
    return;
  }

  if (data.status === 'needs_protein_choice') {
    choices.hidden = false;
    choices.innerHTML = `<h2>Which protein did you mean?</h2><p>${escapeHtml(data.message)}</p>${renderWarnings(data.warnings)}` +
      data.options.map(opt => `<div class="choice"><div><strong>${escapeHtml(opt.gene)}</strong><div>${escapeHtml(opt.organism)}${opt.scientific_name && opt.scientific_name !== opt.organism ? ` · <i>${escapeHtml(opt.scientific_name)}</i>` : ''}</div>${(opt.reasons || []).map(reason => `<div class="small">${escapeHtml(reason)}</div>`).join('')}<div class="small">${escapeHtml(opt.record_count_in_search)} UniProt search record(s) carried the direct name evidence.</div></div><button data-gene="${escapeHtml(opt.gene)}" data-taxon="${opt.taxon_id}" data-species="${escapeHtml(opt.organism)}">Explain all records</button></div>`).join('');
    choices.querySelectorAll('[data-gene]').forEach(btn => btn.addEventListener('click', async () => {
      showStatus('Fetching the full UniProt record set…');
      try {
        const payload = await post('/api/explain', {gene:btn.dataset.gene, taxon_id:Number(btn.dataset.taxon), species_name:btn.dataset.species});
        showStatus(''); renderReady(payload);
      } catch (err) { showStatus(err.message, true); }
    }));
    return;
  }

  if (data.status === 'no_direct_match') {
    results.hidden = false;
    const organismNote = data.organism?.name
      ? `<p class="small"><strong>Organism resolved:</strong> ${escapeHtml(data.organism.input_phrase || data.organism.name)} → ${escapeHtml(data.organism.name)}. Discovery is not restricted to that organism; other UniProt source organisms can appear when the evidence is traceable.</p>`
      : '';
    const relatedNote = data.search_mode === 'related' && data.related_from
      ? `<p class="small"><strong>Related search:</strong> These results use “${escapeHtml(data.query)}” because the full concept “${escapeHtml(data.related_from)}” produced no explainable result.</p>`
      : '';
    results.innerHTML = `<div class="summary"><h2>No direct identity match</h2><p>${escapeHtml(data.message)}</p>${organismNote}${relatedNote}<p>The cards below are discovery options. Each one shows the exact evidence that caused it to appear.</p></div>${renderWarnings(data.warnings)}${renderDiscoveryCards(data.discovery || [])}`;
    bindDiscoveryButtons(results);
    return;
  }

  showStatus(data.message || 'No result was returned.', data.status !== 'ready');
}

function renderReady(data) {
  results.hidden = false;
  choices.hidden = true;
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

  const directReasons = (data.search_context?.direct_reasons || []);
  const contextHtml = directReasons.length ? `<div class="evidence-box"><strong>Why this was treated as a direct identity match:</strong><ul>${directReasons.map(x => `<li>${escapeHtml(x)}</li>`).join('')}</ul></div>` : '';

  results.innerHTML = `${renderWarnings(data.warnings)}<section class="summary"><h2>UniProt has ${escapeHtml(data.record_count)} record${data.record_count === 1 ? '' : 's'} for ${escapeHtml(data.gene)} in ${escapeHtml(data.species)}.</h2>${contextHtml}<p>${escapeHtml(data.review_summary)}</p><p>${escapeHtml(data.review_explanation)}</p>${data.isoform_explanation ? `<p><strong>Isoforms:</strong> ${escapeHtml(data.isoform_explanation)}</p>` : ''}<p class="small">The cards below are comparisons, not rankings or recommendations.</p></section><div class="grid">${cards}</div>`;
}

function discoveryCard(item) {
  const inspect = item.inspect || {};
  const inspectAttrs = inspect.mode === 'gene' && inspect.gene && inspect.taxon_id
    ? `data-inspect-gene="${escapeHtml(inspect.gene)}" data-inspect-taxon="${escapeHtml(inspect.taxon_id)}" data-inspect-species="${escapeHtml(inspect.species_name || item.organism || '')}"`
    : (inspect.accession ? `data-inspect-accession="${escapeHtml(inspect.accession)}" data-inspect-species="${escapeHtml(inspect.species_name || item.organism || '')}"` : '');
  const buttonLabel = inspect.mode === 'gene' && inspect.gene ? `Explain all ${escapeHtml(inspect.gene)} records` : 'Inspect this UniProt entry';
  const evidenceCount = item.evidence_record_count > 1 ? `<p class="small">${escapeHtml(item.evidence_record_count)} records for this gene carried qualifying evidence in the discovery search.</p>` : '';

  return `<article class="card">
    <h3>${escapeHtml(item.protein_name || item.gene || item.accession)}</h3>
    <div class="meta">${escapeHtml(item.gene || 'No gene label')} · ${escapeHtml(item.organism || 'Unknown organism')} · ${escapeHtml(item.accession || '')}</div>
    ${item.relationship ? `<p class="small"><strong>${escapeHtml(item.relationship)}</strong></p>` : ''}
    ${item.connection_explanation ? `<p>${escapeHtml(item.connection_explanation)}</p>` : ''}
    ${item.function_note && item.evidence_type !== 'function' ? `<p><strong>What UniProt says it does:</strong> ${escapeHtml(item.function_note)}</p>` : ''}
    ${evidenceCount}
    <strong>Evidence shown by UniProt</strong>
    ${(item.why || []).map(ctx => `<div class="evidence-box"><div class="small">${escapeHtml(ctx.source)}</div>${escapeHtml(ctx.text)}</div>`).join('')}
    ${inspectAttrs ? `<div class="links"><button class="inspect-discovery" ${inspectAttrs}>${buttonLabel}</button></div>` : ''}
  </article>`;
}

function renderDiscoveryCards(items) {
  if (!items.length) return `<div class="summary"><p>No traceable discovery result was found.</p></div>`;

  const groups = [
    ['mentioned_organism', 'Proteins from the organism you mentioned'],
    ['virus_host', 'Viral proteins linked to the mentioned host'],
    ['global', 'Other traceable UniProt connections'],
  ];

  const sections = groups.map(([key, label]) => {
    const rows = items.filter(item => (item.bucket || 'global') === key);
    if (!rows.length) return '';
    return `<section><h2>${escapeHtml(label)}</h2><div class="grid">${rows.map(discoveryCard).join('')}</div></section>`;
  }).join('');

  return sections || `<div class="grid">${items.map(discoveryCard).join('')}</div>`;
}

function bindDiscoveryButtons(root) {
  root.querySelectorAll('[data-inspect-gene]').forEach(btn => btn.addEventListener('click', async () => {
    showStatus('Fetching the full UniProt record set…');
    try {
      const payload = await post('/api/explain', {
        gene: btn.dataset.inspectGene,
        taxon_id: Number(btn.dataset.inspectTaxon),
        species_name: btn.dataset.inspectSpecies,
      });
      showStatus(''); renderReady(payload);
      window.scrollTo({top: 0, behavior: 'smooth'});
    } catch (err) { showStatus(err.message, true); }
  }));

  root.querySelectorAll('[data-inspect-accession]').forEach(btn => btn.addEventListener('click', async () => {
    showStatus('Fetching the UniProt entry…');
    try {
      const payload = await post('/api/entry', {
        accession: btn.dataset.inspectAccession,
        species_name: btn.dataset.inspectSpecies,
      });
      showStatus(''); renderReady(payload);
      window.scrollTo({top: 0, behavior: 'smooth'});
    } catch (err) { showStatus(err.message, true); }
  }));
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
    bindDiscoveryButtons($('discovery-results'));
  } catch (err) {
    $('discovery-results').innerHTML = `<div class="status warning">${escapeHtml(err.message)}</div>`;
  }
});
