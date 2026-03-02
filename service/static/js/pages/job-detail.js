/**
 * Job detail page (#/jobs/:id)
 * Shows pipeline stepper, live log streaming, and results summary.
 */
const MODULE_STEPS = [
  { tag: 'M1', name: 'Target Prep', step: 1 },
  { tag: 'M2', name: 'Library Screen', step: 2 },
  { tag: 'M3', name: 'De Novo Design', step: 3 },
  { tag: 'M4', name: 'Structure Pred', step: 4 },
  { tag: 'M4.5', name: 'Docking', step: 5 },
  { tag: 'M4.6', name: 'Perturbation', step: 6 },
  { tag: 'M7', name: 'ADMET', step: 7 },
  { tag: 'M5', name: 'Scoring', step: 8 },
  { tag: 'M8', name: 'Delivery', step: 9 },
  { tag: 'M9', name: 'Report', step: 10 },
];

let _logAbort = null;

function renderJobDetail(params) {
  stopLogStream();
  const el = document.getElementById('app-content');
  el.innerHTML = '<div style="text-align:center;padding:3rem"><div class="spinner spinner-lg"></div></div>';
  loadJobDetail(params.id);
}

async function loadJobDetail(jobId) {
  const el = document.getElementById('app-content');
  try {
    const job = await api.get(`/api/jobs/${jobId}`);
    const fmtDate = (d) => d ? new Date(d).toLocaleString() : '—';

    // Action buttons
    let actions = '';
    if (job.status === 'running' || job.status === 'pending') {
      actions += `<button class="btn btn-danger btn-sm" onclick="cancelJob('${job.id}')">Cancel</button>`;
    }
    if (job.status === 'completed') {
      actions += `<a href="#/jobs/${job.id}/dashboard" class="btn btn-primary btn-sm">View Dashboard</a>`;
      actions += `<button class="btn btn-sm" onclick="rerunJob('${job.id}')">Rerun</button>`;
    }

    el.innerHTML = `
      <div class="detail-header">
        <h1>${escapeHTML(job.target)} <span class="badge badge-${job.status}">${job.status}</span></h1>
        <div class="detail-actions">${actions}</div>
      </div>
      <div class="detail-timestamps">
        <span>Created: ${fmtDate(job.created_at)}</span>
        <span>Started: ${fmtDate(job.started_at)}</span>
        <span>Completed: ${fmtDate(job.completed_at)}</span>
      </div>
      <div class="card" style="margin-bottom:1.5rem">
        <div class="card-header">Pipeline Progress</div>
        <div class="card-body">
          <div class="pipeline-stepper" id="stepper" style="padding-bottom:2rem">
            ${buildStepperHTML(job.progress_step)}
          </div>
        </div>
      </div>
      <div id="log-section"></div>
      <div id="results-section"></div>
    `;

    if (job.status === 'running') {
      renderLogStream(jobId);
    } else if (job.status === 'completed') {
      renderResultsSummary(jobId);
    } else if (job.status === 'failed' && job.error_message) {
      document.getElementById('log-section').innerHTML = `
        <div class="card"><div class="card-header">Error</div>
        <div class="card-body" style="color:var(--danger)">${escapeHTML(job.error_message)}</div></div>
      `;
    }
  } catch (err) {
    el.innerHTML = `<div class="empty-state"><h3>Job not found</h3><p>${err.message}</p><a href="#/jobs" class="btn">Back to Jobs</a></div>`;
  }
}

function buildStepperHTML(currentStep) {
  return MODULE_STEPS.map((m, i) => {
    let cls = '';
    if (m.step < currentStep) cls = 'completed';
    else if (m.step === currentStep) cls = 'active';
    const line = i < MODULE_STEPS.length - 1 ? `<div class="step-line ${m.step < currentStep ? '' : ''}"></div>` : '';
    return `
      <div class="step ${cls}" id="step-${m.step}">
        <div class="step-circle">
          ${cls === 'completed' ? '✓' : m.tag}
          <div class="step-label">${m.name}</div>
        </div>
      </div>
      ${line}
    `;
  }).join('');
}

function updateStepper(step) {
  MODULE_STEPS.forEach(m => {
    const stepEl = document.getElementById(`step-${m.step}`);
    if (!stepEl) return;
    stepEl.className = 'step';
    if (m.step < step) stepEl.className = 'step completed';
    else if (m.step === step) stepEl.className = 'step active';
    const circle = stepEl.querySelector('.step-circle');
    if (circle) {
      circle.childNodes[0].textContent = m.step < step ? '✓' : m.tag;
    }
  });
}

async function renderLogStream(jobId) {
  const section = document.getElementById('log-section');
  section.innerHTML = `
    <div class="card">
      <div class="card-header">Live Log Stream <span class="spinner" style="width:16px;height:16px;border-width:2px"></span></div>
      <div class="card-body" style="padding:0">
        <div class="log-stream" id="log-output"></div>
      </div>
    </div>
  `;
  const logEl = document.getElementById('log-output');
  try {
    for await (const event of api.streamSSE(`/api/results/${jobId}/logs`)) {
      if (!document.getElementById('log-output')) break; // navigated away
      if (event.type === 'progress') {
        updateStepper(event.step);
        appendLog(logEl, event.message, 'progress-line');
      } else if (event.type === 'log') {
        appendLog(logEl, event.message);
      } else if (event.type === 'complete') {
        appendLog(logEl, '✓ Pipeline complete!', 'progress-line');
        showToast('Pipeline completed!', 'success');
        setTimeout(() => loadJobDetail(jobId), 1000);
        break;
      } else if (event.type === 'error' || event.type === 'timeout') {
        appendLog(logEl, event.message, 'error-line');
        break;
      }
    }
  } catch (e) {
    if (document.getElementById('log-output')) {
      appendLog(logEl, 'Log stream disconnected', 'error-line');
    }
  }
}

function appendLog(container, text, cls) {
  const line = document.createElement('div');
  line.className = 'log-line' + (cls ? ` ${cls}` : '');
  line.textContent = text;
  container.appendChild(line);
  container.scrollTop = container.scrollHeight;
}

function stopLogStream() {
  // SSE streams auto-stop when navigating away since the element check fails
}

async function renderResultsSummary(jobId) {
  const section = document.getElementById('results-section');
  try {
    const summary = await api.get(`/api/results/${jobId}`);
    const candidates = await api.get(`/api/results/${jobId}/candidates?page_size=5&sort_by=composite_score&descending=true`);

    section.innerHTML = `
      <div class="results-section">
        <h2>Results Summary</h2>
        <div class="kpi-grid" style="grid-template-columns:repeat(3,1fr)">
          <div class="kpi-card green">
            <div class="kpi-label">Total Candidates</div>
            <div class="kpi-value">${summary.total_candidates ?? 0}</div>
          </div>
          <div class="kpi-card">
            <div class="kpi-label">Top Composite Score</div>
            <div class="kpi-value">${summary.top_score != null ? summary.top_score.toFixed(3) : '—'}</div>
          </div>
          <div class="kpi-card orange">
            <div class="kpi-label">Output Files</div>
            <div class="kpi-value">${summary.files?.length ?? 0}</div>
          </div>
        </div>
        ${candidates.candidates.length ? `
          <h3 style="margin-bottom:.75rem">Top 5 Candidates</h3>
          <div class="table-wrap">
            <table>
              <thead><tr>
                <th>Rank</th><th>ID</th><th>Type</th><th>Binding</th><th>ADMET</th><th>Composite</th>
              </tr></thead>
              <tbody>
                ${candidates.candidates.map(c => `
                  <tr>
                    <td>${c.rank}</td>
                    <td>${escapeHTML(c.candidate_id)}</td>
                    <td>${escapeHTML(c.candidate_type || c.modality)}</td>
                    <td>${c.binding_score.toFixed(3)}</td>
                    <td>${c.admet_score.toFixed(3)}</td>
                    <td><strong>${c.composite_score.toFixed(3)}</strong></td>
                  </tr>
                `).join('')}
              </tbody>
            </table>
          </div>
        ` : ''}
        <div style="margin-top:1rem">
          <a href="#/jobs/${jobId}/dashboard" class="btn btn-primary">View Full Dashboard</a>
        </div>
      </div>
    `;
  } catch (e) {
    section.innerHTML = `<p style="color:var(--text-secondary)">Results not available yet.</p>`;
  }
}

async function cancelJob(jobId) {
  if (!confirm('Cancel this pipeline run?')) return;
  try {
    await api.del(`/api/jobs/${jobId}`);
    showToast('Job cancelled', 'info');
    loadJobDetail(jobId);
  } catch (err) {
    showToast(err.message, 'error');
  }
}

async function rerunJob(jobId) {
  try {
    const job = await api.post(`/api/jobs/${jobId}/rerun`);
    showToast(`Rerun submitted as ${job.id}`, 'success');
    router.navigate(`#/jobs/${job.id}`);
  } catch (err) {
    showToast(err.message, 'error');
  }
}
