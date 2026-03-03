/**
 * Pricing page with subscription tiers and PortOne SDK integration (#/pricing)
 * Accessible without login — subscribe buttons require login.
 */

const PLANS = [
  {
    tier: 'free',
    name: 'Free',
    price: '₩0',
    period: '',
    desc: 'Get started with basic pipeline runs',
    features: [
      { text: '3 pipeline runs total', has: true },
      { text: 'Standard docking', has: true },
      { text: 'Basic results dashboard', has: true },
      { text: 'GPU acceleration', has: false },
      { text: 'Priority queue', has: false },
      { text: 'Dedicated support', has: false },
    ],
    btnClass: 'pricing-btn-outline',
    btnText: 'Get Started Free',
  },
  {
    tier: 'pro',
    name: 'Pro',
    price: '₩39,000',
    period: '/월',
    desc: 'For active researchers needing GPU power',
    featured: true,
    features: [
      { text: '50 pipeline runs/month', has: true },
      { text: 'Standard docking', has: true },
      { text: 'Full results dashboard', has: true },
      { text: 'GPU acceleration (RunPod)', has: true },
      { text: 'Priority queue', has: true },
      { text: 'Dedicated support', has: false },
    ],
    btnClass: 'pricing-btn-primary',
    btnText: 'Subscribe',
  },
  {
    tier: 'enterprise',
    name: 'Enterprise',
    price: '₩129,000',
    period: '/월',
    desc: 'Unlimited runs with premium support',
    features: [
      { text: 'Unlimited pipeline runs', has: true },
      { text: 'Standard docking', has: true },
      { text: 'Full results dashboard', has: true },
      { text: 'GPU acceleration (RunPod)', has: true },
      { text: 'Priority queue', has: true },
      { text: 'Dedicated support', has: true },
    ],
    btnClass: 'pricing-btn-primary',
    btnText: 'Subscribe',
  },
];

const SERVICE_FEATURES = [
  {
    icon: '🧬',
    title: 'Target Preparation',
    desc: 'Automated UniProt sequence retrieval, PDB structure resolution, binding pocket identification, and anti-target profiling.',
  },
  {
    icon: '🗄️',
    title: 'Multi-Database Screening',
    desc: 'Integrated search across ChEMBL, BindingDB, PubChem, DrugBank, UniProt, and PDB with circuit-breaker resilience.',
  },
  {
    icon: '💊',
    title: 'Dual Modality',
    desc: 'Support for both small molecule and peptide discovery with modality-aware ADMET profiling and scoring.',
  },
  {
    icon: '🔬',
    title: 'De Novo Design',
    desc: 'AI-powered candidate generation using BRICS fragmentation (small molecules) and PepMLM (peptides).',
  },
  {
    icon: '⚡',
    title: 'GPU-Accelerated Docking',
    desc: 'AutoDock Vina molecular docking with optional RunPod GPU dispatch for faster structure prediction.',
  },
  {
    icon: '📊',
    title: 'Interactive Dashboards',
    desc: 'Per-run HTML reports with 5-tab Chart.js dashboards, candidate ranking, and cross-run comparison.',
  },
];

async function renderPricing() {
  const el = document.getElementById('app-content');
  const loggedIn = api.isLoggedIn();

  let quota = { tier: 'free', runs_used: 0, runs_limit: 3 };
  let payments = [];

  if (loggedIn) {
    try {
      [quota, payments] = await Promise.all([
        api.get('/api/subscription/'),
        api.get('/api/subscription/payments').catch(() => []),
      ]);
    } catch (e) {
      // fallback to defaults
    }
  }

  // ── Service overview section ──
  const serviceHtml = SERVICE_FEATURES.map(f => `
    <div class="service-feature-card">
      <div class="service-feature-icon">${f.icon}</div>
      <h4>${f.title}</h4>
      <p>${f.desc}</p>
    </div>
  `).join('');

  // ── Plan cards ──
  const cardsHtml = PLANS.map(plan => {
    const isCurrent = loggedIn && plan.tier === quota.tier;
    const featuresHtml = plan.features.map(f =>
      `<li><span class="${f.has ? 'check' : 'cross'}">${f.has ? '✓' : '✕'}</span> ${f.text}</li>`
    ).join('');

    let btnHtml;
    if (!loggedIn) {
      if (plan.tier === 'free') {
        btnHtml = `<a href="#/register" class="pricing-btn pricing-btn-primary" style="text-decoration:none;display:block;text-align:center">Sign Up Free</a>`;
      } else {
        btnHtml = `<a href="#/register" class="pricing-btn ${plan.btnClass}" style="text-decoration:none;display:block;text-align:center">Sign Up to Subscribe</a>`;
      }
    } else if (isCurrent) {
      btnHtml = `<button class="pricing-btn pricing-btn-outline" disabled>Current Plan</button>`;
    } else if (plan.tier === 'free') {
      btnHtml = `<button class="pricing-btn pricing-btn-outline" disabled>Free Tier</button>`;
    } else {
      btnHtml = `<button class="pricing-btn ${plan.btnClass}" onclick="subscribeToPlan('${plan.tier}')">
        ${plan.btnText}
      </button>`;
    }

    return `
      <div class="pricing-card${plan.featured ? ' featured' : ''}${isCurrent ? ' current' : ''}">
        ${isCurrent ? '<span class="pricing-current-badge">Current Plan</span>' : ''}
        <div class="pricing-tier-name">${plan.name}</div>
        <div class="pricing-price">${plan.price}<span class="period">${plan.period}</span></div>
        <div class="pricing-desc">${plan.desc}</div>
        <ul class="pricing-features">${featuresHtml}</ul>
        ${btnHtml}
      </div>
    `;
  }).join('');

  // ── Payments section (logged in only) ──
  let paymentsHtml = '';
  if (loggedIn && payments.length > 0) {
    const rows = payments.map(p => `
      <tr>
        <td>${new Date(p.created_at).toLocaleDateString()}</td>
        <td>${p.tier || '—'}</td>
        <td>₩${p.amount.toLocaleString()}</td>
        <td>${p.status}</td>
      </tr>
    `).join('');
    paymentsHtml = `
      <div class="payments-section">
        <h3>Payment History</h3>
        <table class="payments-table">
          <thead><tr><th>Date</th><th>Plan</th><th>Amount</th><th>Status</th></tr></thead>
          <tbody>${rows}</tbody>
        </table>
      </div>
    `;
  }

  let cancelHtml = '';
  if (loggedIn && quota.tier !== 'free') {
    cancelHtml = `
      <div style="text-align:center;margin-top:1.5rem">
        <button class="pricing-btn pricing-btn-outline" style="max-width:300px" onclick="cancelSubscription()">
          Cancel Subscription
        </button>
      </div>
    `;
  }

  // ── Pipeline overview ──
  const pipelineSteps = [
    { num: '1', name: 'Target Prep', tag: 'UniProt/PDB' },
    { num: '2', name: 'Library Screen', tag: 'ChEMBL' },
    { num: '3', name: 'De Novo Design', tag: 'AI Generation' },
    { num: '4', name: 'Structure Prediction', tag: 'ESMFold' },
    { num: '5', name: 'Molecular Docking', tag: 'AutoDock Vina' },
    { num: '6', name: 'Perturbation Biology', tag: 'CMap' },
    { num: '7', name: 'ADMET Profiling', tag: 'Safety' },
    { num: '8', name: 'Composite Scoring', tag: 'Ranking' },
    { num: '9', name: 'Delivery Design', tag: 'Tissue-specific' },
    { num: '10', name: 'Final Report', tag: 'Dashboard' },
  ];
  const pipelineHtml = pipelineSteps.map(s => `
    <div class="pipeline-stage">
      <div class="stage-num">${s.num}</div>
      <div class="stage-name">${s.name}</div>
      <div class="stage-tag">${s.tag}</div>
    </div>
  `).join('');

  el.innerHTML = `
    <div class="pricing-page">
      <div class="pricing-header">
        <h2>AI-Powered Drug Discovery Platform</h2>
        <p>From target gene to ranked, delivery-ready candidates in a single automated run</p>
      </div>

      <div class="pricing-service-section">
        <h3>What You Get</h3>
        <div class="service-features-grid">${serviceHtml}</div>
      </div>

      <div class="pricing-pipeline-section">
        <h3>10-Module Pipeline</h3>
        <div class="pipeline-stages">${pipelineHtml}</div>
      </div>

      <div class="pricing-header" style="margin-top:2rem">
        <h2>Choose Your Plan</h2>
        <p>Scale your drug discovery pipeline with GPU-accelerated runs</p>
      </div>
      <div class="pricing-grid">${cardsHtml}</div>
      ${cancelHtml}
      ${paymentsHtml}

      <div class="pricing-faq">
        <h3>Frequently Asked Questions</h3>
        <div class="faq-grid">
          <div class="faq-item">
            <h4>What is a pipeline run?</h4>
            <p>Each run takes a target gene and produces a ranked list of drug candidates through all 10 modules — target analysis, screening, docking, ADMET, scoring, and reporting.</p>
          </div>
          <div class="faq-item">
            <h4>What does GPU acceleration do?</h4>
            <p>GPU-accelerated runs use RunPod cloud GPUs for faster molecular docking and structure prediction, reducing run times significantly.</p>
          </div>
          <div class="faq-item">
            <h4>Can I cancel anytime?</h4>
            <p>Yes. Cancel at any time and keep access until the end of your billing period. No long-term contracts.</p>
          </div>
          <div class="faq-item">
            <h4>What modalities are supported?</h4>
            <p>Both small molecules and peptides. The pipeline adapts its screening, design, ADMET, and delivery modules based on the chosen modality.</p>
          </div>
        </div>
      </div>

      <div class="pricing-contact">
        <strong>Questions?</strong> Contact us at <a href="mailto:brownbio.ocm@gmail.com">brownbio.ocm@gmail.com</a>
      </div>
    </div>
  `;
}

async function subscribeToPlan(tier) {
  if (!api.isLoggedIn()) {
    showToast('Please sign in first', 'error');
    router.navigate('#/login');
    return;
  }

  if (typeof PortOne === 'undefined') {
    showToast('Payment SDK not loaded. Please refresh the page.', 'error');
    return;
  }

  try {
    // 1. Request billing key via PortOne popup
    const response = await PortOne.requestIssueBillingKey({
      storeId: window.__PORTONE_STORE_ID || '',
      channelKey: window.__PORTONE_CHANNEL_KEY || '',
      billingKeyMethod: 'CARD',
    });

    if (response.code) {
      showToast(response.message || 'Card registration failed', 'error');
      return;
    }

    // 2. Send billing key to backend
    await api.post('/api/subscription/register-billing-key', {
      billing_key: response.billingKey,
    });

    // 3. Subscribe to the selected tier
    await api.post('/api/subscription/subscribe', { tier });

    showToast(`Subscribed to ${tier.charAt(0).toUpperCase() + tier.slice(1)} plan!`, 'success');
    router.navigate('#/pricing');
  } catch (err) {
    showToast(err.message || 'Subscription failed', 'error');
  }
}

async function cancelSubscription() {
  if (!confirm('Are you sure you want to cancel? You will keep access until the end of your billing period.')) {
    return;
  }
  try {
    const result = await api.post('/api/subscription/cancel');
    showToast(result.message || 'Subscription cancelled', 'info');
    router.navigate('#/pricing');
  } catch (err) {
    showToast(err.message || 'Cancel failed', 'error');
  }
}
