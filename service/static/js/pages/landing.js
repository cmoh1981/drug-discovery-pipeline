/**
 * Landing / dashboard page (#/)
 */
function renderLanding() {
  const el = document.getElementById('app-content');
  el.innerHTML = `
    <div class="landing-hero">
      <div class="health-banner" id="health-banner"></div>
      <h1>Drug Discovery <span class="accent">Pipeline</span></h1>
      <p class="tagline">AI-powered drug candidate discovery — from target gene to ranked, delivery-ready candidates in a single run.</p>
      <div class="feature-grid">
        <div class="feature-card">
          <div class="feature-icon">🧬</div>
          <h3>10 Modules</h3>
          <p>Target prep, screening, de novo design, docking, ADMET, scoring, and more</p>
        </div>
        <div class="feature-card">
          <div class="feature-icon">🗄️</div>
          <h3>6+ Databases</h3>
          <p>ChEMBL, BindingDB, UniProt, PDB, DrugBank, PubChem integration</p>
        </div>
        <div class="feature-card">
          <div class="feature-icon">💊</div>
          <h3>Small Molecules & Peptides</h3>
          <p>Dual modality support with modality-aware ADMET and delivery systems</p>
        </div>
        <div class="feature-card">
          <div class="feature-icon">📊</div>
          <h3>Interactive Reports</h3>
          <p>Rich dashboards with charts, tables, filtering, and dark/light theme</p>
        </div>
      </div>
      <div class="hero-cta">
        <a href="#/jobs/new" class="btn btn-primary btn-lg">Launch New Pipeline</a>
      </div>
    </div>
    <div class="recent-jobs-section" id="recent-jobs"></div>
  `;
  loadHealth();
  if (api.isLoggedIn()) loadRecentJobs();
}

async function loadHealth() {
  try {
    const h = await api.get('/api/health');
    const banner = document.getElementById('health-banner');
    if (!banner) return;
    const cls = h.status === 'ok' ? 'ok' : 'degraded';
    const dot = h.status === 'ok' ? '●' : '◐';
    banner.className = `health-banner ${cls}`;
    banner.innerHTML = `${dot} System ${h.status} &middot; DB ${h.database} &middot; ${h.active_jobs} active job${h.active_jobs !== 1 ? 's' : ''}`;
  } catch (e) {
    // health check is best-effort
  }
}

async function loadRecentJobs() {
  try {
    const data = await api.get('/api/jobs?page_size=3');
    const container = document.getElementById('recent-jobs');
    if (!container || !data.jobs.length) return;
    container.innerHTML = `
      <h2>Recent Jobs</h2>
      <div class="jobs-grid" style="grid-template-columns:repeat(3,1fr)">
        ${data.jobs.map(j => jobCardHTML(j)).join('')}
      </div>
    `;
  } catch (e) {
    // not logged in or no jobs
  }
}
