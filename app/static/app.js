const $ = (id) => document.getElementById(id);
const statusBox = $("status");
const choices = $("choices");
const results = $("results");

const responseCache = new Map();
let activeController = null;
let currentReady = null;
let currentFilter = "all";
let currentRecordSearch = "";
let visibleLimit = 12;
let selectedAccessions = new Set();

function escapeHtml(value) {
  return String(value ?? "").replace(/[&<>'"]/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;',"'":'&#39;','"':'&quot;'}[ch]));
}

function showStatus(message, warning=false, loading=false) {
  statusBox.hidden = !message;
  statusBox.className = `status ${warning ? 'warning' : ''}`;
  statusBox.innerHTML = message ? `${loading ? '<div class="loading"><span class="spinner"></span><span>' : ''}${escapeHtml(message)}${loading ? '</span></div>' : ''}` : "";
}

function renderWarnings(warnings=[]) {
  if (!warnings.length) return "";
  return `<div class="status warning"><strong>Source note</strong><ul>${warnings.map(w => `<li>${escapeHtml(w.message)}</li>`).join('')}</ul></div>`;
}


function normalizeTextBlock(value) {
  return String(value ?? '').replace(/\s+/g, ' ').trim();
}

function renderMechanismBlock(mechanism, compact=false) {
  if (!mechanism) return '';
  const sentences = (mechanism.summary || []).map(normalizeTextBlock).filter(Boolean);
  if (!sentences.length) return '';
  const sourceRows = (mechanism.sources || []).map(src => {
    const field = escapeHtml(src.field || 'UniProt annotation');
    const text = escapeHtml(normalizeTextBlock(src.text || ''));
    return `<li><strong>${field}</strong><div>${text}</div></li>`;
  }).join('');
  return `<section class="biology-block mechanism-block">
    <div class="biology-kicker">How this protein works</div>
    <div class="mechanism-summary">${sentences.map((sentence, index) => `<p class="${index === 0 ? 'biology-lead' : 'mechanism-step'}">${escapeHtml(sentence)}</p>`).join('')}</div>
    ${mechanism.evidence_note ? `<p class="evidence-note">${escapeHtml(mechanism.evidence_note)}</p>` : ''}
    ${sourceRows && !compact ? `<details class="mechanism-provenance"><summary>How UniProt supports this explanation</summary><ul>${sourceRows}</ul></details>` : ''}
  </section>`;
}

function renderDiseaseItems(diseases=[], limit=2) {
  if (!diseases.length) return '';
  const visible = diseases.slice(0, limit);
  const hidden = diseases.slice(limit);
  const item = d => `<li><strong>${escapeHtml(d.name || d.acronym || 'Disease annotation')}</strong>${d.acronym && d.name ? ` <span class="muted-inline">(${escapeHtml(d.acronym)})</span>` : ''}${d.description ? `<div>${escapeHtml(normalizeTextBlock(d.description))}</div>` : ''}${d.note ? `<div class="disease-note">${escapeHtml(normalizeTextBlock(d.note))}</div>` : ''}${d.evidence_summary ? `<div class="evidence-note">${escapeHtml(d.evidence_summary)}</div>` : ''}</li>`;
  return `<section class="biology-block disease-block">
    <div class="biology-kicker">Disease relevance</div>
    <p class="small biology-source-note">UniProt exposes ${diseases.length} disease annotation${diseases.length === 1 ? '' : 's'} for this entry.</p>
    <ul class="disease-list">${visible.map(item).join('')}</ul>
    ${hidden.length ? `<details><summary>Show ${hidden.length} more disease annotation${hidden.length === 1 ? '' : 's'}</summary><ul class="disease-list">${hidden.map(item).join('')}</ul></details>` : ''}
  </section>`;
}

function renderBiologicalOverview(overview) {
  if (!overview) return '';
  const hasMechanism = (overview.mechanism?.summary || []).length;
  const hasDiseases = (overview.diseases || []).length;
  if (!hasMechanism && !hasDiseases) return '';
  const sourceLabel = overview.source_review?.label || 'UniProtKB';
  return `<section class="biological-overview">
    <div class="overview-heading"><div><p class="eyebrow">Biological explanation</p><h3>How ${escapeHtml(overview.protein_name || 'this protein')} works</h3></div><div class="source-pill">Source: ${escapeHtml(overview.source_accession || '')} · ${escapeHtml(sourceLabel)}</div></div>
    ${renderMechanismBlock(overview.mechanism, false)}
    ${renderDiseaseItems(overview.diseases || [], 2)}
  </section>`;
}

async function post(url, body, signal=null) {
  const response = await fetch(url, {
    method: 'POST',
    headers: {'Content-Type':'application/json'},
    body: JSON.stringify(body),
    signal,
  });
  let payload = null;
  try { payload = await response.json(); } catch (_) {}
  if (!response.ok) throw new Error(payload?.detail || payload?.message || `Request failed (${response.status})`);
  return payload;
}

function lookupCacheKey(body) {
  return JSON.stringify(['lookup', body]);
}

async function cachedPost(url, body, {signal=null, cache=true}={}) {
  const key = JSON.stringify([url, body]);
  if (cache && responseCache.has(key)) return responseCache.get(key);
  const payload = await post(url, body, signal);
  if (cache) responseCache.set(key, payload);
  return payload;
}

function geneDownloadUrl(gene, taxonId, format) {
  return `/api/download/gene?gene=${encodeURIComponent(gene)}&taxon_id=${encodeURIComponent(taxonId)}&format=${encodeURIComponent(format)}`;
}

function entryDownloadUrl(accession, format) {
  return `/api/download/entry/${encodeURIComponent(accession)}?format=${encodeURIComponent(format)}`;
}

function formatDownloadLinks(urlBuilder, label='Download') {
  return `<div class="download-group"><span class="download-label">${escapeHtml(label)}:</span>${['fasta','tsv','json','xml','txt'].map(fmt => `<a class="download-link" href="${urlBuilder(fmt)}">${fmt === 'txt' ? 'UniProt TXT' : fmt.toUpperCase()}</a>`).join('')}</div>`;
}

function compactDownloadMenu(urlBuilder, label='Download') {
  return `<details class="download-menu"><summary>${escapeHtml(label)}</summary><div class="download-menu-body">${['fasta','tsv','json','xml','txt'].map(fmt => `<a class="download-link" href="${urlBuilder(fmt)}">${fmt === 'txt' ? 'UniProt TXT' : fmt.toUpperCase()}</a>`).join('')}</div></details>`;
}

async function downloadSelected(format) {
  const ids = [...selectedAccessions];
  if (!ids.length) {
    showStatus('Select at least one record first.', true);
    return;
  }
  showStatus(`Preparing ${ids.length} selected record${ids.length === 1 ? '' : 's'} as ${format.toUpperCase()}…`, false, true);
  try {
    const response = await fetch('/api/download/selected', {
      method: 'POST',
      headers: {'Content-Type':'application/json'},
      body: JSON.stringify({accessions: ids, format}),
    });
    if (!response.ok) {
      let msg = `Download failed (${response.status})`;
      try { msg = (await response.json()).detail || msg; } catch (_) {}
      throw new Error(msg);
    }
    const blob = await response.blob();
    const disposition = response.headers.get('Content-Disposition') || '';
    const match = disposition.match(/filename="?([^";]+)"?/i);
    const filename = match?.[1] || `uniprot_selected_${ids.length}.${format}`;
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
    showStatus('');
  } catch (err) {
    showStatus(err.message, true);
  }
}

async function runLookup(query, overrides={}) {
  if (activeController) activeController.abort();
  activeController = new AbortController();
  choices.hidden = true;
  results.hidden = true;
  currentReady = null;
  selectedAccessions.clear();
  const body = {query, ...overrides};
  const started = performance.now();
  showStatus('Searching UniProt…', false, true);
  const slowTimer = setTimeout(() => {
    showStatus('Still working — checking UniProt records and traceable evidence…', false, true);
  }, 8000);
  const hardTimer = setTimeout(() => {
    activeController?.abort();
  }, 45000);

  try {
    const data = await cachedPost('/api/lookup', body, {signal: activeController.signal});
    clearTimeout(slowTimer);
    clearTimeout(hardTimer);
    const elapsed = Math.max(0.1, (performance.now() - started) / 1000).toFixed(1);
    showStatus(`Search completed in ${elapsed}s.`);
    setTimeout(() => { if (!statusBox.classList.contains('warning')) showStatus(''); }, 1800);
    handleLookupResponse(data, query);
  } catch (err) {
    clearTimeout(slowTimer);
    clearTimeout(hardTimer);
    if (err.name === 'AbortError') {
      showStatus('This search took longer than 45 seconds, so it was stopped. Try again; cached UniProt responses may make the next attempt faster.', true);
      return;
    }
    showStatus(err.message, true);
  }
}
function handleLookupResponse(data, originalQuery) {
  if (data.status === 'ready') {
    renderReady(data);
    return;
  }

  if (data.status === 'needs_organism_choice') {
    choices.hidden = false;
    choices.innerHTML = `<div class="summary"><h2>Which organism did you mean?</h2><p>${escapeHtml(data.message)}</p></div>${renderWarnings(data.warnings)}` +
      data.options.map(opt => `<div class="choice"><div><strong>${escapeHtml(opt.label)}</strong><div class="small">${escapeHtml(opt.rank || '')} · Taxonomy ID ${escapeHtml(opt.taxon_id)}</div></div><button data-organism-taxon="${escapeHtml(opt.taxon_id)}" data-organism-name="${escapeHtml(opt.label)}" data-organism-phrase="${escapeHtml(opt.input_phrase)}">Use this organism</button></div>`).join('') +
      `<div class="choice"><div><strong>Continue without choosing an organism</strong><div class="small">The search will stay broad and label the source organism for each result.</div></div><button class="secondary" data-no-organism>Continue broadly</button></div>`;

    choices.querySelectorAll('[data-organism-taxon]').forEach(btn => btn.addEventListener('click', () => {
      runLookup(originalQuery, {
        taxon_id: Number(btn.dataset.organismTaxon),
        organism_name: btn.dataset.organismName,
        organism_phrase: btn.dataset.organismPhrase,
        resolve_organism: false,
      });
    }));
    choices.querySelector('[data-no-organism]')?.addEventListener('click', () => runLookup(originalQuery, {resolve_organism:false}));
    return;
  }

  if (data.status === 'needs_protein_choice') {
    choices.hidden = false;
    choices.innerHTML = `<div class="summary"><h2>Choose the matching protein or gene</h2><p>${escapeHtml(data.message)}</p></div>${renderWarnings(data.warnings)}` +
      data.options.map(opt => {
        const dl = opt.gene && opt.taxon_id ? compactDownloadMenu(fmt => geneDownloadUrl(opt.gene, opt.taxon_id, fmt), 'Download record set') : '';
        return `<div class="choice"><div><strong>${escapeHtml(opt.gene)}</strong><div>${escapeHtml(opt.organism)}${opt.scientific_name && opt.scientific_name !== opt.organism ? ` · <i>${escapeHtml(opt.scientific_name)}</i>` : ''}</div>${(opt.reasons || []).map(reason => `<div class="small">${escapeHtml(reason)}</div>`).join('')}<div class="small">${escapeHtml(opt.record_count_in_search)} search record(s) carried direct-name evidence.</div><div class="links">${dl}</div></div><div class="choice-actions"><button data-gene="${escapeHtml(opt.gene)}" data-taxon="${escapeHtml(opt.taxon_id)}" data-species="${escapeHtml(opt.organism)}">Open record set</button></div></div>`;
      }).join('');
    choices.querySelectorAll('[data-gene]').forEach(btn => btn.addEventListener('click', async () => {
      await openGeneSet(btn.dataset.gene, Number(btn.dataset.taxon), btn.dataset.species);
    }));
    return;
  }

  if (data.status === 'no_direct_match') {
    results.hidden = false;
    const organismNote = data.organism?.name
      ? `<p class="small"><strong>Organism resolved:</strong> ${escapeHtml(data.organism.input_phrase || data.organism.name)} → ${escapeHtml(data.organism.name)}. Other source organisms can still appear when the evidence is traceable.</p>`
      : '';
    const relatedNote = data.search_mode === 'related' && data.related_from
      ? `<p class="small"><strong>Related search:</strong> These results use “${escapeHtml(data.query)}” because the full concept “${escapeHtml(data.related_from)}” produced no explainable result.</p>`
      : '';
    const spelling = data.can_check_spelling && data.organism?.taxon_id
      ? `<button class="secondary" id="check-spelling" data-phrase="${escapeHtml(data.query)}" data-taxon="${escapeHtml(data.organism.taxon_id)}">Check spelling suggestions</button>`
      : '';
    results.innerHTML = `<div class="summary"><h2>${(data.discovery || []).length ? 'Traceable UniProt connections' : 'No traceable UniProt match'}</h2><p>${escapeHtml(data.message)}</p>${organismNote}${relatedNote}${spelling ? `<div class="links">${spelling}</div>` : ''}</div>${renderWarnings(data.warnings)}${renderDiscoveryCards(data.discovery || [])}`;
    bindDiscoveryButtons(results);
    results.querySelector('#check-spelling')?.addEventListener('click', async (event) => {
      const button = event.currentTarget;
      button.disabled = true;
      button.textContent = 'Checking…';
      try {
        const spellingData = await post('/api/spelling', {phrase:button.dataset.phrase, taxon_id:Number(button.dataset.taxon)});
        if (spellingData.suggestion?.phrase) {
          button.outerHTML = `<button id="use-spelling" class="secondary">Search “${escapeHtml(spellingData.suggestion.phrase)}”</button>`;
          results.querySelector('#use-spelling')?.addEventListener('click', () => runLookup(originalQuery, {confirmed_spelling: spellingData.suggestion.phrase}));
        } else {
          button.textContent = 'No confident spelling suggestion found';
        }
      } catch (err) {
        button.disabled = false;
        button.textContent = 'Check spelling suggestions';
        showStatus(err.message, true);
      }
    });
    return;
  }

  showStatus(data.message || 'No result was returned.', data.status !== 'ready');
}

async function openGeneSet(gene, taxonId, speciesName) {
  const body = {gene, taxon_id: taxonId, species_name: speciesName};
  showStatus('Fetching the full UniProt record set…', false, true);
  try {
    const payload = await cachedPost('/api/explain', body);
    showStatus('');
    renderReady(payload);
    window.scrollTo({top: 0, behavior: 'smooth'});
  } catch (err) {
    showStatus(err.message, true);
  }
}

async function openAccession(accession, speciesName) {
  showStatus('Fetching the UniProt entry…', false, true);
  try {
    const payload = await cachedPost('/api/entry', {accession, species_name:speciesName});
    showStatus('');
    renderReady(payload);
    window.scrollTo({top: 0, behavior: 'smooth'});
  } catch (err) {
    showStatus(err.message, true);
  }
}

function recordMatchesFilter(item) {
  const review = item.record?.review?.status || 'unknown';
  if (currentFilter !== 'all' && review !== currentFilter) return false;
  if (!currentRecordSearch) return true;
  const haystack = [item.accession, item.record?.name, item.record?.gene, item.record?.transcript].filter(Boolean).join(' ').toLowerCase();
  return haystack.includes(currentRecordSearch.toLowerCase());
}

function recordCard(item) {
  const r = item.record || {};
  const review = r.review || {};
  const mechanism = r.mechanism || {};
  const links = r.links || {};
  const isoforms = r.isoforms || [];
  const accession = item.accession;
  const pdb = links.pdb || [];
  const structureLink = pdb.length
    ? `<a href="${links.pdb_search}" target="_blank" rel="noopener">PDB structures (${pdb.length})</a>`
    : (links.pdb_search ? `<a href="${links.pdb_search}" target="_blank" rel="noopener">Search PDB</a>` : '');
  const explanation = (item.sentences || []).length
    ? `<ul class="difference-list">${item.sentences.map(x => `<li>${escapeHtml(x)}</li>`).join('')}</ul>`
    : `<p class="small">No major differences were exposed by the comparison fields used in this prototype.</p>`;

  return `<article class="card record-card" data-record-card data-review="${escapeHtml(review.status || 'unknown')}" data-accession="${escapeHtml(accession)}">
    <div class="card-head">
      <div><span class="badge ${escapeHtml(review.status)}">${escapeHtml(review.label)}</span><h3>${escapeHtml(accession)}</h3></div>
      <label class="select-control"><input class="card-select" type="checkbox" aria-label="Select ${escapeHtml(accession)}" data-select-accession="${escapeHtml(accession)}" ${selectedAccessions.has(accession) ? 'checked' : ''}><span>Select</span></label>
    </div>
    <div class="meta">${escapeHtml(r.name)} · ${escapeHtml(item.length ?? 'unknown')} aa${r.gene ? ` · ${escapeHtml(r.gene)}` : ''}</div>
    ${renderMechanismBlock(mechanism)}
    ${renderDiseaseItems(r.disease_evidence || [], 2)}
    <div class="card-section"><div class="card-section-title">How this record differs from the others</div>${explanation}</div>
    ${isoforms.length ? `<details><summary>Isoforms described in this entry (${isoforms.length})</summary><ul class="isoform-list">${isoforms.map(i => { const ids = Array.isArray(i.ids) ? i.ids.filter(Boolean).join(', ') : ''; const label = ids || i.name || 'Unnamed isoform'; const status = i.sequence_status || 'status unavailable'; return `<li><strong>${escapeHtml(label)}</strong>${i.name && ids && i.name !== ids ? ` · ${escapeHtml(i.name)}` : ''} <span class="muted-inline">— ${escapeHtml(status)}</span></li>`; }).join('')}</ul></details>` : ''}
    <details><summary>Database details</summary><div class="detail-grid"><div><span>Protein existence</span><strong>${escapeHtml(r.existence)}</strong></div><div><span>Ensembl transcript</span><strong>${escapeHtml(r.transcript || 'Not returned')}</strong></div></div></details>
    <div class="record-actions">
      ${links.uniprot ? `<a class="button-link secondary" href="${links.uniprot}" target="_blank" rel="noopener">UniProt</a>` : ''}
      ${links.alphafold ? `<a class="button-link secondary" href="${links.alphafold}" target="_blank" rel="noopener">AlphaFold</a>` : ''}
      ${structureLink ? `<span class="structure-link">${structureLink}</span>` : ''}
      ${compactDownloadMenu(fmt => entryDownloadUrl(accession, fmt), 'Download record')}
    </div>
  </article>`;
}
function renderReady(data) {
  currentReady = data;
  currentFilter = 'all';
  currentRecordSearch = '';
  visibleLimit = 12;
  selectedAccessions.clear();
  choices.hidden = true;
  results.hidden = false;

  const directReasons = data.search_context?.direct_reasons || [];
  const contextHtml = directReasons.length ? `<div class="identity-note"><strong>Why this matched directly</strong><ul>${directReasons.map(x => `<li>${escapeHtml(x)}</li>`).join('')}</ul></div>` : '';
  const differenceSummary = (data.difference_summary || []).length
    ? `<div class="difference-summary"><h3>What differs across these records</h3><ul>${data.difference_summary.map(x => `<li>${escapeHtml(x)}</li>`).join('')}</ul></div>`
    : '';
  const allDownloads = data.download_scope?.mode === 'gene'
    ? formatDownloadLinks(fmt => geneDownloadUrl(data.download_scope.gene, data.download_scope.taxon_id, fmt), 'Download full record set')
    : formatDownloadLinks(fmt => entryDownloadUrl(data.accessions?.[0], fmt), 'Download record');

  results.innerHTML = `${renderWarnings(data.warnings)}
    <section class="summary result-summary">
      <div class="summary-topline"><div><p class="eyebrow">Record explanation</p><h2>${escapeHtml(data.gene)} <span>in ${escapeHtml(data.species)}</span></h2></div><div class="record-count"><strong>${escapeHtml(data.record_count)}</strong><span>UniProtKB record${data.record_count === 1 ? '' : 's'}</span></div></div>
      <div class="review-strip"><span class="badge reviewed">${escapeHtml(data.reviewed_count ?? 0)} reviewed</span><span class="badge unreviewed">${escapeHtml(data.unreviewed_count ?? 0)} unreviewed</span></div>
      ${contextHtml}
      ${renderBiologicalOverview(data.biological_overview)}
      ${differenceSummary}
      <div class="summary-actions">${allDownloads}</div>
      <div class="learn-row">
        <details><summary>What does reviewed vs unreviewed mean?</summary><p>${escapeHtml(data.review_explanation)}</p></details>
        ${data.isoform_explanation ? `<details><summary>What is a UniProt isoform?</summary><p>${escapeHtml(data.isoform_explanation)}</p></details>` : ''}
      </div>
      ${data.download_scope?.mode === 'gene' ? `<div class="secondary-context"><button class="tertiary" id="load-external">Load Ensembl / APPRIS context</button><div id="external-panel" class="external-panel"></div></div>` : ''}
    </section>
    <section class="toolbar">
      <div class="toolbar-row">
        <div class="segmented" role="group" aria-label="Filter records">
          <button class="filter-button active" data-filter="all">All ${escapeHtml(data.record_count)}</button>
          <button class="filter-button secondary" data-filter="reviewed">Reviewed ${escapeHtml(data.reviewed_count ?? 0)}</button>
          <button class="filter-button secondary" data-filter="unreviewed">Unreviewed ${escapeHtml(data.unreviewed_count ?? 0)}</button>
        </div>
        <input id="record-search" placeholder="Filter records by accession, name, gene or transcript" />
      </div>
      <div class="toolbar-row selection-row">
        <button class="secondary" id="select-visible">Select visible</button>
        <button class="tertiary" id="clear-selection">Clear</button>
        <span class="small" id="selection-count">0 selected</span>
        <div class="selected-actions"><span class="download-label">Download selected</span>${['fasta','tsv','json','xml','txt'].map(fmt => `<button class="secondary selected-download" data-format="${fmt}">${fmt === 'txt' ? 'TXT' : fmt.toUpperCase()}</button>`).join('')}</div>
      </div>
    </section>
    <div id="record-grid" class="grid"></div>
    <button id="show-more" class="secondary show-more" hidden>Show more records</button>
    <div id="empty-records" class="summary empty" hidden>No records match this filter.</div>`;

  bindReadyControls();
  renderRecordGrid();
}
function bindReadyControls() {
  results.querySelectorAll('[data-filter]').forEach(btn => btn.addEventListener('click', () => {
    currentFilter = btn.dataset.filter;
    visibleLimit = 12;
    results.querySelectorAll('[data-filter]').forEach(x => x.classList.toggle('active', x === btn));
    renderRecordGrid();
  }));

  results.querySelector('#record-search')?.addEventListener('input', (event) => {
    currentRecordSearch = event.target.value.trim();
    visibleLimit = 12;
    renderRecordGrid();
  });

  results.querySelector('#show-more')?.addEventListener('click', () => {
    visibleLimit += 24;
    renderRecordGrid();
  });

  results.querySelector('#select-visible')?.addEventListener('click', () => {
    getVisibleRecordItems().forEach(item => selectedAccessions.add(item.accession));
    renderRecordGrid();
  });

  results.querySelector('#clear-selection')?.addEventListener('click', () => {
    selectedAccessions.clear();
    renderRecordGrid();
  });

  results.querySelectorAll('.selected-download').forEach(btn => btn.addEventListener('click', () => downloadSelected(btn.dataset.format)));

  results.querySelector('#load-external')?.addEventListener('click', async (event) => {
    const button = event.currentTarget;
    button.disabled = true;
    button.textContent = 'Loading…';
    const panel = results.querySelector('#external-panel');
    try {
      const data = await cachedPost('/api/external-context', {gene:currentReady.gene, taxon_id:currentReady.taxon_id});
      const e = data.external_annotations || {};
      panel.innerHTML = `${renderWarnings(data.warnings)}<div class="evidence-box"><strong>Linked annotation context</strong><p>Ensembl gene: ${escapeHtml(e.ensembl_gene || 'Not available')}</p><p>Canonical transcript: ${escapeHtml(e.ensembl_canonical || 'Not available')}</p><p>UniProt gene-centric representative: ${escapeHtml(e.gene_centric_accession || 'Not available')}</p><p>APPRIS annotations loaded: ${escapeHtml(Object.keys(e.appris || {}).length)}</p></div>`;
      button.remove();
    } catch (err) {
      button.disabled = false;
      button.textContent = 'Load Ensembl / APPRIS context';
      panel.innerHTML = `<div class="status warning">${escapeHtml(err.message)}</div>`;
    }
  });
}

function getFilteredRecordItems() {
  if (!currentReady) return [];
  const items = (currentReady.records || []).filter(recordMatchesFilter);
  return items.sort((a, b) => {
    const rank = {reviewed:0, unreviewed:1, unknown:2};
    const ar = rank[a.record?.review?.status] ?? 3;
    const br = rank[b.record?.review?.status] ?? 3;
    return ar - br || String(a.accession).localeCompare(String(b.accession));
  });
}

function getVisibleRecordItems() {
  return getFilteredRecordItems().slice(0, visibleLimit);
}

function updateSelectionCount() {
  const node = results.querySelector('#selection-count');
  if (node) node.textContent = `${selectedAccessions.size} selected`;
}

function renderRecordGrid() {
  const grid = results.querySelector('#record-grid');
  if (!grid) return;
  const all = getFilteredRecordItems();
  const visible = all.slice(0, visibleLimit);
  grid.innerHTML = visible.map(recordCard).join('');
  results.querySelector('#empty-records').hidden = all.length !== 0;
  const more = results.querySelector('#show-more');
  more.hidden = visible.length >= all.length;
  if (!more.hidden) more.textContent = `Show more (${all.length - visible.length} remaining)`;

  grid.querySelectorAll('[data-select-accession]').forEach(box => box.addEventListener('change', () => {
    if (box.checked) selectedAccessions.add(box.dataset.selectAccession);
    else selectedAccessions.delete(box.dataset.selectAccession);
    updateSelectionCount();
  }));
  updateSelectionCount();
}

function discoveryCard(item) {
  const inspect = item.inspect || {};
  const isGene = inspect.mode === 'gene' && inspect.gene && inspect.taxon_id;
  const inspectAttrs = isGene
    ? `data-inspect-gene="${escapeHtml(inspect.gene)}" data-inspect-taxon="${escapeHtml(inspect.taxon_id)}" data-inspect-species="${escapeHtml(inspect.species_name || item.organism || '')}"`
    : (inspect.accession ? `data-inspect-accession="${escapeHtml(inspect.accession)}" data-inspect-species="${escapeHtml(inspect.species_name || item.organism || '')}"` : '');
  const buttonLabel = isGene ? `Open ${escapeHtml(inspect.gene)} record set` : 'Open UniProt entry';
  const evidenceCount = item.evidence_record_count > 1 ? `<p class="small">${escapeHtml(item.evidence_record_count)} search records for this gene carried qualifying evidence.</p>` : '';
  const dl = isGene
    ? compactDownloadMenu(fmt => geneDownloadUrl(inspect.gene, inspect.taxon_id, fmt), 'Download record set')
    : (inspect.accession ? compactDownloadMenu(fmt => entryDownloadUrl(inspect.accession, fmt), 'Download entry') : '');
  const uniLink = item.accession ? `<a href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(item.accession)}/entry" target="_blank" rel="noopener">Evidence entry</a>` : '';
  const contexts = item.why || [];
  const visibleEvidence = contexts.slice(0, 2).map(ctx => `<div class="evidence-box"><div class="small">${escapeHtml(ctx.source)}</div>${escapeHtml(ctx.text)}</div>`).join('');
  const moreEvidence = contexts.length > 2 ? `<details><summary>More evidence (${contexts.length - 2})</summary>${contexts.slice(2).map(ctx => `<div class="evidence-box"><div class="small">${escapeHtml(ctx.source)}</div>${escapeHtml(ctx.text)}</div>`).join('')}</details>` : '';

  return `<article class="card discovery-card">
    <div class="card-head"><div><h3>${escapeHtml(item.protein_name || item.gene || item.accession)}</h3><div class="meta">${escapeHtml(item.gene || 'No gene label')} · ${escapeHtml(item.organism || 'Unknown organism')} · ${escapeHtml(item.accession || '')}</div></div></div>
    ${item.relationship ? `<span class="relationship-chip">${escapeHtml(item.relationship)}</span>` : ''}
    ${item.connection_explanation ? `<p>${escapeHtml(item.connection_explanation)}</p>` : ''}
    ${item.function_note && item.evidence_type !== 'function' ? `<p class="small"><strong>UniProt function:</strong> ${escapeHtml(item.function_note)}</p>` : ''}
    ${evidenceCount}
    <div class="card-section-title">Why this appeared</div>
    ${visibleEvidence}${moreEvidence}
    <div class="record-actions">${inspectAttrs ? `<button class="inspect-discovery" ${inspectAttrs}>${buttonLabel}</button>` : ''}${uniLink}${dl}</div>
  </article>`;
}
function renderDiscoveryCards(items) {
  if (!items.length) return `<div class="summary empty"><p>No traceable discovery result was found.</p></div>`;
  const groups = [
    ['mentioned_organism', 'Proteins from the organism you mentioned'],
    ['virus_host', 'Viral proteins linked to the mentioned host'],
    ['global', 'Other traceable UniProt connections'],
  ];
  const sections = groups.map(([key, label]) => {
    const rows = items.filter(item => (item.bucket || 'global') === key);
    if (!rows.length) return '';
    return `<section><h2 class="section-title">${escapeHtml(label)}</h2><div class="grid">${rows.map(discoveryCard).join('')}</div></section>`;
  }).join('');
  return sections || `<div class="grid">${items.map(discoveryCard).join('')}</div>`;
}

function bindDiscoveryButtons(root) {
  root.querySelectorAll('[data-inspect-gene]').forEach(btn => btn.addEventListener('click', () => {
    openGeneSet(btn.dataset.inspectGene, Number(btn.dataset.inspectTaxon), btn.dataset.inspectSpecies);
  }));
  root.querySelectorAll('[data-inspect-accession]').forEach(btn => btn.addEventListener('click', () => {
    openAccession(btn.dataset.inspectAccession, btn.dataset.inspectSpecies);
  }));
}

$('lookup-form').addEventListener('submit', (event) => {
  event.preventDefault();
  const query = $('query').value.trim();
  if (query) runLookup(query);
});
