/**
 * Cross-run comparison page (#/compare)
 */
async function renderCompare() {
  const el = document.getElementById('app-content');
  el.innerHTML = `
    <h1 style="margin-bottom:1.25rem">Compare Pipeline Runs</h1>
    <div class="card">
      <div class="card-header">Select Completed Jobs to Compare</div>
      <div class="card-body">
        <div id="compare-jobs-list" class="compare-jobs-list">
          <div class="spinner"></div>
        </div>
        <div class="compare-controls" style="margin-top:1rem">
          <div class="form-group" style="margin-bottom:0;flex:0 0 auto">
            <label for="compare-topn">Top N</label>
            <input type="number" id="compare-topn" class="form-input" value="20" min="1" max="100" style="width:80px">
          </div>
          <button class="btn btn-primary" id="compare-btn" onclick="runComparison()">Compare Selected</button>
        </div>
      </div>
    </div>
    <div id="compare-results" style="margin-top:1.5rem"></div>
  `;
  loadCompareJobs();
}

async function loadCompareJobs() {
  const container = document.getElementById('compare-jobs-list');
  try {
    const data = await api.get('/api/jobs/?status=completed&page_size=100');
    if (!data.jobs.length) {
      container.innerHTML = '<p style="color:var(--text-secondary)">No completed jobs available for comparison.</p>';
      return;
    }
    container.innerHTML = data.jobs.map(j => {
      const created = new Date(j.created_at).toLocaleDateString(undefined, { month: 'short', day: 'numeric' });
      return `
        <label class="compare-job-check">
          <input type="checkbox" value="${j.id}" class="compare-checkbox">
          <strong>${escapeHTML(j.target)}</strong>
          <span class="badge badge-sm-mol" style="font-size:.65rem">${j.modality}</span>
          <span style="color:var(--text-secondary);font-size:.75rem">${created}</span>
        </label>
      `;
    }).join('');
  } catch (err) {
    container.innerHTML = `<p style="color:var(--danger)">${err.message}</p>`;
  }
}

async function runComparison() {
  const checked = Array.from(document.querySelectorAll('.compare-checkbox:checked')).map(cb => cb.value);
  if (checked.length < 2) {
    showToast('Select at least 2 jobs to compare', 'error');
    return;
  }
  const topN = parseInt(document.getElementById('compare-topn').value, 10) || 20;
  const btn = document.getElementById('compare-btn');
  const results = document.getElementById('compare-results');

  btn.disabled = true;
  btn.textContent = 'Comparing...';
  results.innerHTML = '<div style="text-align:center;padding:2rem"><div class="spinner spinner-lg"></div></div>';

  try {
    const data = await api.post('/api/compare', { job_ids: checked, top_n: topN });
    let html = `
      <div class="card">
        <div class="card-header">Comparison Results (${data.total_compared} candidates compared)</div>
        <div class="card-body">
    `;

    // Unique counts
    html += `<div class="kpi-grid" style="grid-template-columns:repeat(${Object.keys(data.unique_per_job).length},1fr);margin-bottom:1rem">`;
    for (const [jid, count] of Object.entries(data.unique_per_job)) {
      html += `
        <div class="kpi-card">
          <div class="kpi-label">Unique to ${jid.slice(0, 8)}</div>
          <div class="kpi-value">${count}</div>
        </div>
      `;
    }
    html += '</div>';

    // Shared candidates table
    if (data.shared_candidates.length) {
      html += `
        <h3 style="margin-bottom:.75rem">Shared Candidates (${data.shared_candidates.length})</h3>
        <div class="table-wrap"><table>
          <thead><tr>
            <th>Identifier</th>
            ${data.job_ids.map(jid => `<th>${jid.slice(0, 8)}...</th>`).join('')}
          </tr></thead>
          <tbody>
            ${data.shared_candidates.map(c => `
              <tr>
                <td>${escapeHTML(c.identifier)}</td>
                ${data.job_ids.map(jid => `<td>${c.scores[jid] != null ? c.scores[jid].toFixed(3) : '—'}</td>`).join('')}
              </tr>
            `).join('')}
          </tbody>
        </table></div>
      `;
    } else {
      html += '<p style="color:var(--text-secondary)">No shared candidates found between selected runs.</p>';
    }

    html += '</div></div>';
    results.innerHTML = html;
  } catch (err) {
    results.innerHTML = `<div class="card"><div class="card-body" style="color:var(--danger)">${err.message}</div></div>`;
  } finally {
    btn.disabled = false;
    btn.textContent = 'Compare Selected';
  }
}
