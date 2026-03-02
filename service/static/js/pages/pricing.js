/**
 * Pricing page with subscription tiers and PortOne SDK integration (#/pricing)
 */

const PLANS = [
  {
    tier: 'free',
    name: 'Free',
    price: '$0',
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
    btnText: 'Current Plan',
  },
  {
    tier: 'pro',
    name: 'Pro',
    price: '$29',
    period: '/mo',
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
    price: '$99',
    period: '/mo',
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

async function renderPricing() {
  const el = document.getElementById('app-content');
  el.innerHTML = '<div class="pricing-page"><p>Loading...</p></div>';

  let quota = { tier: 'free', runs_used: 0, runs_limit: 3 };
  let payments = [];
  try {
    [quota, payments] = await Promise.all([
      api.get('/api/subscription/'),
      api.get('/api/subscription/payments').catch(() => []),
    ]);
  } catch (e) {
    // Not logged in or error — show default
  }

  const cardsHtml = PLANS.map(plan => {
    const isCurrent = plan.tier === quota.tier;
    const featuresHtml = plan.features.map(f =>
      `<li><span class="${f.has ? 'check' : 'cross'}">${f.has ? '✓' : '✕'}</span> ${f.text}</li>`
    ).join('');

    let btnHtml;
    if (isCurrent) {
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

  let paymentsHtml = '';
  if (payments.length > 0) {
    const rows = payments.map(p => `
      <tr>
        <td>${new Date(p.created_at).toLocaleDateString()}</td>
        <td>${p.tier || '—'}</td>
        <td>$${(p.amount / 100).toFixed(2)} ${p.currency}</td>
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
  if (quota.tier !== 'free') {
    cancelHtml = `
      <div style="text-align:center;margin-top:1.5rem">
        <button class="pricing-btn pricing-btn-outline" style="max-width:300px" onclick="cancelSubscription()">
          Cancel Subscription
        </button>
      </div>
    `;
  }

  el.innerHTML = `
    <div class="pricing-page">
      <div class="pricing-header">
        <h2>Choose Your Plan</h2>
        <p>Scale your drug discovery pipeline with GPU-accelerated runs</p>
      </div>
      <div class="pricing-grid">${cardsHtml}</div>
      ${cancelHtml}
      ${paymentsHtml}
    </div>
  `;
}

async function subscribeToPlan(tier) {
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
