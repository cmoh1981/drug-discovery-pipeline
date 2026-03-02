/**
 * Landing / dashboard page (#/)
 */
function renderLanding() {
  const el = document.getElementById('app-content');
  el.innerHTML = `
    <div class="landing-hero">
      <div class="health-banner" id="health-banner"></div>
      <div class="company-name">Brown Biotech Inc.</div>
      <h1>Drug Discovery <span class="accent">Pipeline</span></h1>
      <p class="tagline">AI-powered computational drug discovery — from target gene to ranked, delivery-ready candidates in a single automated run.</p>
      <div class="hero-cta">
        <a href="#/jobs/new" class="btn btn-primary btn-lg">Launch Pipeline</a>
        <a href="#/jobs" class="btn btn-outline btn-lg">View Results</a>
      </div>
      <div class="metrics-bar">
        <div class="metric">
          <div class="metric-value">10</div>
          <div class="metric-label">Pipeline Modules</div>
        </div>
        <div class="metric">
          <div class="metric-value">6+</div>
          <div class="metric-label">Integrated Databases</div>
        </div>
        <div class="metric">
          <div class="metric-value">2</div>
          <div class="metric-label">Modalities</div>
        </div>
        <div class="metric">
          <div class="metric-value">500+</div>
          <div class="metric-label">Candidates / Run</div>
        </div>
      </div>
    </div>

    <div class="pricing-notice">
      <div class="pricing-icon">&#x1F9EA;</div>
      <div class="pricing-text">
        <strong>First 3 pipeline runs are free</strong>
        <span>For additional runs and enterprise access, contact us at <a href="mailto:brownbio.ocm@gmail.com">brownbio.ocm@gmail.com</a></span>
      </div>
    </div>

    <div class="section-header">
      <div class="section-tag">Platform Capabilities</div>
      <h2>End-to-End Computational Drug Discovery</h2>
    </div>
    <div class="feature-grid">
      <div class="feature-card">
        <div class="feature-icon blue">&#x1F9EC;</div>
        <h3>Target Preparation</h3>
        <p>UniProt sequence retrieval, PDB structure resolution, and binding-site identification via automated workflows</p>
      </div>
      <div class="feature-card">
        <div class="feature-icon teal">&#x1F5C4;</div>
        <h3>Multi-Database Screening</h3>
        <p>ChEMBL, BindingDB, UniProt, PDB, DrugBank, and PubChem integration with circuit-breaker resilience</p>
      </div>
      <div class="feature-card">
        <div class="feature-icon purple">&#x1F48A;</div>
        <h3>Dual Modality</h3>
        <p>Small molecule and peptide discovery with modality-aware ADMET profiling and delivery system selection</p>
      </div>
      <div class="feature-card">
        <div class="feature-icon green">&#x1F4CA;</div>
        <h3>Interactive Dashboards</h3>
        <p>Per-run reports with 5-tab Chart.js dashboards, candidate ranking, and cross-run comparison tools</p>
      </div>
    </div>

    <div class="pipeline-overview">
      <h3>Pipeline Architecture</h3>
      <div class="pipeline-stages">
        <div class="pipeline-stage">
          <div class="stage-num">1</div>
          <div class="stage-name">Target Prep</div>
          <div class="stage-tag">UniProt/PDB</div>
        </div>
        <div class="pipeline-stage">
          <div class="stage-num">2</div>
          <div class="stage-name">Ligand Screen</div>
          <div class="stage-tag">ChEMBL</div>
        </div>
        <div class="pipeline-stage">
          <div class="stage-num">3</div>
          <div class="stage-name">Compound Enrich</div>
          <div class="stage-tag">PubChem</div>
        </div>
        <div class="pipeline-stage">
          <div class="stage-num">4</div>
          <div class="stage-name">De Novo Design</div>
          <div class="stage-tag">BRICS</div>
        </div>
        <div class="pipeline-stage">
          <div class="stage-num">5</div>
          <div class="stage-name">BindingDB</div>
          <div class="stage-tag">Cross-ref</div>
        </div>
        <div class="pipeline-stage">
          <div class="stage-num">6</div>
          <div class="stage-name">DrugBank</div>
          <div class="stage-tag">Safety</div>
        </div>
        <div class="pipeline-stage">
          <div class="stage-num">7</div>
          <div class="stage-name">Docking</div>
          <div class="stage-tag">Vina/GPU</div>
        </div>
        <div class="pipeline-stage">
          <div class="stage-num">8</div>
          <div class="stage-name">ADMET</div>
          <div class="stage-tag">Profiling</div>
        </div>
        <div class="pipeline-stage">
          <div class="stage-num">9</div>
          <div class="stage-name">Scoring</div>
          <div class="stage-tag">Ranking</div>
        </div>
        <div class="pipeline-stage">
          <div class="stage-num">10</div>
          <div class="stage-name">Reporting</div>
          <div class="stage-tag">Dashboard</div>
        </div>
      </div>
    </div>

    <div class="contact-bar">
      <div class="contact-info">
        <strong>Brown Biotech Inc.</strong> &mdash; <a href="mailto:brownbio.ocm@gmail.com">brownbio.ocm@gmail.com</a>
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
    const dot = h.status === 'ok' ? '\u25CF' : '\u25D0';
    banner.className = `health-banner ${cls}`;
    banner.innerHTML = `${dot} System ${h.status} &middot; DB ${h.database} &middot; ${h.active_jobs} active job${h.active_jobs !== 1 ? 's' : ''}`;
  } catch (e) {
    // health check is best-effort
  }
}

async function loadRecentJobs() {
  try {
    const data = await api.get('/api/jobs/?page_size=3');
    const container = document.getElementById('recent-jobs');
    if (!container || !data.jobs.length) return;
    container.innerHTML = `
      <div class="section-header">
        <div class="section-tag">Recent Activity</div>
        <h2>Your Latest Pipeline Runs</h2>
      </div>
      <div class="jobs-grid" style="grid-template-columns:repeat(3,1fr)">
        ${data.jobs.map(j => jobCardHTML(j)).join('')}
      </div>
    `;
  } catch (e) {
    // not logged in or no jobs
  }
}
