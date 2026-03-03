/**
 * Job dashboard page (#/jobs/:id/dashboard)
 * Embeds the per-run HTML report in an iframe.
 */
async function renderJobDashboard(params) {
  const el = document.getElementById('app-content');
  el.innerHTML = '<div style="text-align:center;padding:3rem"><div class="spinner spinner-lg"></div></div>';
  try {
    const job = await api.get(`/api/jobs/${params.id}`);
    if (job.status !== 'completed' || !job.result_dir) {
      el.innerHTML = `
        <div class="empty-state">
          <div class="empty-icon">📊</div>
          <h3>Dashboard not available</h3>
          <p>The pipeline must complete before the dashboard is available.</p>
          <a href="#/jobs/${params.id}" class="btn">Back to Job</a>
        </div>
      `;
      return;
    }

    // result_dir is like "job_results/UUID/timestamp_target_..."
    // /reports mount maps to job_results/, so strip that prefix
    // Normalize backslashes (Windows) to forward slashes for URLs
    const resultPath = job.result_dir.replace(/^job_results[\\/]/, '').replace(/\\/g, '/');
    const reportURL = `/reports/${resultPath}/report/final_report.html`;

    el.innerHTML = `
      <div class="dashboard-page">
        <div class="back-bar">
          <a href="#/jobs/${params.id}" class="btn btn-sm">← Back to Job</a>
          <span style="margin-left:.75rem;font-weight:600">${escapeHTML(job.target)} — Interactive Dashboard</span>
        </div>
        <iframe src="${reportURL}" id="dashboard-iframe" frameborder="0" allowfullscreen></iframe>
      </div>
    `;

    // Sync theme with iframe after load
    const iframe = document.getElementById('dashboard-iframe');
    iframe.addEventListener('load', () => {
      try {
        const theme = document.documentElement.getAttribute('data-theme') || 'light';
        iframe.contentDocument.documentElement.setAttribute('data-theme', theme);
      } catch (e) {
        // cross-origin will silently fail, which is fine for same-origin
      }
    });
  } catch (err) {
    el.innerHTML = `<div class="empty-state"><h3>Error</h3><p>${err.message}</p><a href="#/jobs/${params.id}" class="btn">Back</a></div>`;
  }
}
