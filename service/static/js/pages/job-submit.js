/**
 * Job submission form (#/jobs/new)
 */
function renderJobSubmit() {
  const el = document.getElementById('app-content');
  el.innerHTML = `
    <div class="submit-page">
      <div class="card">
        <div class="card-header"><h2>Launch New Pipeline</h2></div>
        <div class="card-body">
          <form id="submit-form">
            <div class="form-group">
              <label for="sub-target">Target Gene *</label>
              <input type="text" id="sub-target" class="form-input" placeholder="e.g. EGFR, GLP1R, ESRRG" required>
              <div class="form-hint">UniProt gene symbol for the drug target</div>
            </div>

            <div class="form-row">
              <div class="form-group">
                <label>Modality *</label>
                <div class="radio-group">
                  <label class="radio-option">
                    <input type="radio" name="modality" value="small_molecule" checked> Small Molecule
                  </label>
                  <label class="radio-option">
                    <input type="radio" name="modality" value="peptide"> Peptide
                  </label>
                </div>
              </div>
              <div class="form-group">
                <label>Mode *</label>
                <div class="radio-group">
                  <label class="radio-option">
                    <input type="radio" name="mode" value="antagonist" checked> Antagonist
                  </label>
                  <label class="radio-option">
                    <input type="radio" name="mode" value="agonist"> Agonist
                  </label>
                </div>
              </div>
            </div>

            <div class="form-group">
              <label for="sub-tissue">Tissue Context</label>
              <input type="text" id="sub-tissue" class="form-input" placeholder="e.g. lung, liver, brain (optional)">
              <div class="form-hint">Tissue type for perturbation biology context</div>
            </div>

            <div class="collapsible-header" onclick="toggleAdvanced(this)">
              <span class="collapse-icon">▶</span> Advanced Options
            </div>
            <div class="collapsible-body" id="advanced-opts">
              <div class="form-group">
                <label for="sub-topn">Top N Candidates</label>
                <div class="range-group">
                  <input type="range" id="sub-topn" min="1" max="500" value="20" oninput="document.getElementById('topn-val').textContent=this.value">
                  <span class="range-value" id="topn-val">20</span>
                </div>
              </div>

              <div class="form-group">
                <label class="checkbox-option">
                  <input type="checkbox" id="sub-runpod"> Use RunPod GPU (faster structure prediction)
                </label>
              </div>

              <div class="form-group" id="template-group">
                <label for="sub-template">Config Template</label>
                <select id="sub-template" class="form-select">
                  <option value="">Default Configuration</option>
                </select>
              </div>
            </div>

            <div style="margin-top:1.5rem">
              <button type="submit" class="btn btn-primary btn-lg" style="width:100%;justify-content:center" id="submit-btn">
                🚀 Launch Pipeline
              </button>
            </div>
          </form>
        </div>
      </div>
    </div>
  `;
  loadTemplates();
  document.getElementById('submit-form').addEventListener('submit', handleSubmit);
}

function toggleAdvanced(header) {
  header.classList.toggle('open');
  document.getElementById('advanced-opts').classList.toggle('open');
}

async function loadTemplates() {
  try {
    const templates = await api.get('/api/jobs/templates');
    const select = document.getElementById('sub-template');
    if (!select) return;
    templates.forEach(t => {
      const opt = document.createElement('option');
      opt.value = t.config_json;
      opt.textContent = t.name + (t.description ? ` — ${t.description}` : '');
      select.appendChild(opt);
    });
  } catch (e) {
    // templates are optional
  }
}

async function handleSubmit(e) {
  e.preventDefault();
  const btn = document.getElementById('submit-btn');
  const target = document.getElementById('sub-target').value.trim();
  if (!target) { showToast('Target gene is required', 'error'); return; }

  const modality = document.querySelector('input[name="modality"]:checked').value;
  const mode = document.querySelector('input[name="mode"]:checked').value;
  const tissue = document.getElementById('sub-tissue').value.trim();
  const top_n = parseInt(document.getElementById('sub-topn').value, 10);
  const use_runpod = document.getElementById('sub-runpod').checked;

  let config_overrides = null;
  const templateVal = document.getElementById('sub-template').value;
  if (templateVal) {
    try { config_overrides = JSON.parse(templateVal); } catch (e) { /* ignore */ }
  }

  btn.disabled = true;
  btn.textContent = 'Submitting...';
  try {
    const job = await api.post('/api/jobs/', { target, modality, mode, tissue, top_n, use_runpod, config_overrides });
    showToast(`Job submitted for ${target}`, 'success');
    router.navigate(`#/jobs/${job.id}`);
  } catch (err) {
    showToast(err.message || 'Submission failed', 'error');
  } finally {
    btn.disabled = false;
    btn.textContent = '🚀 Launch Pipeline';
  }
}
