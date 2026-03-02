/**
 * Jobs gallery page (#/jobs)
 */
let _jobsPollingTimers = [];

function renderJobsGallery() {
  clearJobsPolling();
  const el = document.getElementById('app-content');
  el.innerHTML = `
    <div class="jobs-header">
      <h1>Your Pipeline Jobs</h1>
      <div class="jobs-filters">
        <select class="form-select" id="status-filter" style="width:auto">
          <option value="">All Status</option>
          <option value="pending">Pending</option>
          <option value="running">Running</option>
          <option value="completed">Completed</option>
          <option value="failed">Failed</option>
          <option value="cancelled">Cancelled</option>
        </select>
        <a href="#/jobs/new" class="btn btn-primary">+ New Job</a>
      </div>
    </div>
    <div id="jobs-container"></div>
    <div id="jobs-pagination"></div>
  `;
  document.getElementById('status-filter').addEventListener('change', () => loadJobs(1));
  loadJobs(1);
}

async function loadJobs(page) {
  clearJobsPolling();
  const container = document.getElementById('jobs-container');
  const pagEl = document.getElementById('jobs-pagination');
  if (!container) return;
  const statusFilter = document.getElementById('status-filter')?.value || '';
  const params = `page=${page}&page_size=12${statusFilter ? '&status=' + statusFilter : ''}`;

  container.innerHTML = '<div style="text-align:center;padding:2rem"><div class="spinner spinner-lg"></div></div>';
  try {
    const data = await api.get(`/api/jobs/?${params}`);
    if (!data.jobs.length) {
      container.innerHTML = `
        <div class="empty-state">
          <div class="empty-icon">🧪</div>
          <h3>No jobs yet</h3>
          <p>Launch your first drug discovery pipeline to get started.</p>
          <a href="#/jobs/new" class="btn btn-primary">+ New Job</a>
        </div>
      `;
      pagEl.innerHTML = '';
      return;
    }
    container.innerHTML = `<div class="jobs-grid">${data.jobs.map(j => jobCardHTML(j)).join('')}</div>`;
    // Poll progress for running jobs
    data.jobs.filter(j => j.status === 'running').forEach(j => {
      const timer = setInterval(() => pollJobProgress(j.id), 3000);
      _jobsPollingTimers.push(timer);
    });
    // Pagination
    const totalPages = Math.ceil(data.total / data.page_size);
    if (totalPages > 1) {
      pagEl.innerHTML = `
        <div class="pagination">
          <button class="btn btn-sm" ${page <= 1 ? 'disabled' : ''} onclick="loadJobs(${page - 1})">Prev</button>
          <span class="page-info">Page ${page} of ${totalPages}</span>
          <button class="btn btn-sm" ${page >= totalPages ? 'disabled' : ''} onclick="loadJobs(${page + 1})">Next</button>
        </div>
      `;
    } else {
      pagEl.innerHTML = '';
    }
  } catch (err) {
    container.innerHTML = `<div class="empty-state"><h3>Error loading jobs</h3><p>${err.message}</p></div>`;
  }
}

function jobCardHTML(j) {
  const modalityBadge = j.modality === 'peptide'
    ? '<span class="badge badge-peptide">Peptide</span>'
    : '<span class="badge badge-sm-mol">Small Mol</span>';
  let extraContent = '';
  if (j.status === 'running') {
    const pct = j.progress_total > 0 ? Math.round((j.progress_step / j.progress_total) * 100) : 0;
    extraContent = `
      <div class="job-card-progress">
        <div class="progress-text"><span>Step ${j.progress_step}/${j.progress_total}</span><span>${j.current_module || ''}</span></div>
        <div class="progress-bar"><div class="progress-bar-fill" style="width:${pct}%" id="prog-${j.id}"></div></div>
      </div>
    `;
  } else if (j.status === 'completed') {
    extraContent = `
      <div class="job-card-score">
        <span>Candidates: <strong>${j.total_candidates ?? '—'}</strong></span>
        <span>Top Score: <strong>${j.top_score != null ? j.top_score.toFixed(3) : '—'}</strong></span>
      </div>
    `;
  } else if (j.status === 'failed' && j.error_message) {
    extraContent = `<div class="job-card-error" title="${escapeHTML(j.error_message)}">${escapeHTML(j.error_message)}</div>`;
  }
  const created = new Date(j.created_at).toLocaleDateString(undefined, { month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit' });
  return `
    <a class="job-card status-${j.status}" href="#/jobs/${j.id}" id="jcard-${j.id}">
      <div class="job-card-target">${escapeHTML(j.target)}</div>
      <div class="job-card-badges">
        ${modalityBadge}
        <span class="badge badge-${j.status}">${j.status}</span>
        <span class="badge" style="background:var(--bg-hover);color:var(--text-secondary)">${j.mode}</span>
      </div>
      ${extraContent}
      <div class="job-card-meta">${created}</div>
    </a>
  `;
}

async function pollJobProgress(jobId) {
  try {
    const p = await api.get(`/api/jobs/${jobId}/progress`);
    const bar = document.getElementById(`prog-${jobId}`);
    if (bar) {
      const pct = p.total > 0 ? Math.round((p.step / p.total) * 100) : 0;
      bar.style.width = pct + '%';
      const textEl = bar.closest('.job-card-progress')?.querySelector('.progress-text');
      if (textEl) textEl.innerHTML = `<span>Step ${p.step}/${p.total}</span><span>${p.module || ''}</span>`;
    }
    if (p.status === 'completed' || p.status === 'failed') {
      loadJobs(1); // refresh
    }
  } catch (e) {
    // ignore polling errors
  }
}

function clearJobsPolling() {
  _jobsPollingTimers.forEach(t => clearInterval(t));
  _jobsPollingTimers = [];
}

function escapeHTML(str) {
  if (!str) return '';
  const div = document.createElement('div');
  div.textContent = str;
  return div.innerHTML;
}
