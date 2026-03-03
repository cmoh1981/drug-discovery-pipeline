/**
 * App initialization: route registration, theme toggle, auth guard, sidebar.
 */
(function () {
  // ── Theme ────────────────────────────────────────────
  const THEME_KEY = 'dd_theme';

  function getTheme() {
    return localStorage.getItem(THEME_KEY) || 'light';
  }

  function setTheme(theme) {
    localStorage.setItem(THEME_KEY, theme);
    document.documentElement.setAttribute('data-theme', theme);
    const btn = document.getElementById('theme-toggle');
    if (btn) btn.textContent = theme === 'dark' ? '☀️ Light' : '🌙 Dark';
    // Sync iframe if present
    const iframe = document.getElementById('dashboard-iframe');
    if (iframe) {
      try { iframe.contentDocument.documentElement.setAttribute('data-theme', theme); } catch (e) {}
    }
  }

  function toggleTheme() {
    setTheme(getTheme() === 'dark' ? 'light' : 'dark');
  }

  // ── Auth guard ───────────────────────────────────────
  function requireAuth(handler) {
    return function (params) {
      if (!api.isLoggedIn()) {
        router.navigate('#/login');
        return;
      }
      handler(params);
    };
  }

  // ── Sidebar active state ─────────────────────────────
  function updateSidebar() {
    const hash = window.location.hash || '#/';
    document.querySelectorAll('#sidebar nav a').forEach(a => {
      const href = a.getAttribute('href');
      if (href === hash || (href !== '#/' && hash.startsWith(href))) {
        a.classList.add('active');
      } else {
        a.classList.remove('active');
      }
    });
    // Update page title in top bar
    const titles = {
      '#/': 'Dashboard',
      '#/login': 'Sign In',
      '#/register': 'Create Account',
      '#/jobs': 'Pipeline Jobs',
      '#/jobs/new': 'New Job',
      '#/compare': 'Compare Runs',
      '#/pricing': 'Pricing',
    };
    const titleEl = document.querySelector('.page-title');
    if (titleEl) {
      titleEl.textContent = titles[hash] || 'Drug Discovery Pipeline';
    }
  }

  // ── User menu ────────────────────────────────────────
  window.updateUserMenu = function () {
    const menuEl = document.getElementById('user-menu');
    if (!menuEl) return;
    if (api.isLoggedIn() && window.__ddUser) {
      const initials = (window.__ddUser.full_name || window.__ddUser.email || '?')
        .split(' ').map(w => w[0]).join('').toUpperCase().slice(0, 2);
      menuEl.innerHTML = `
        <div class="user-avatar">${initials}</div>
        <span>${escapeHTML(window.__ddUser.email)}</span>
        <button id="logout-btn" onclick="doLogout()">Logout</button>
      `;
    } else {
      menuEl.innerHTML = `<a href="#/login" class="btn btn-sm btn-primary">Sign In</a>`;
    }
  };

  window.doLogout = function () {
    api.clearToken();
    window.__ddUser = null;
    updateUserMenu();
    showToast('Logged out', 'info');
    router.navigate('#/login');
  };

  // ── Toast helper ─────────────────────────────────────
  window.showToast = function (msg, type) {
    type = type || 'info';
    const container = document.getElementById('toast-container');
    if (!container) return;
    const toast = document.createElement('div');
    toast.className = `toast toast-${type}`;
    toast.textContent = msg;
    container.appendChild(toast);
    setTimeout(() => toast.remove(), 4000);
  };

  // ── Mobile sidebar toggle ────────────────────────────
  window.toggleSidebar = function () {
    document.getElementById('sidebar').classList.toggle('open');
    document.getElementById('sidebar-overlay').classList.toggle('open');
  };

  // ── Route registration ───────────────────────────────
  router.add('#/', renderLanding);
  router.add('#/login', renderLogin);
  router.add('#/register', renderRegister);
  router.add('#/jobs', requireAuth(renderJobsGallery));
  router.add('#/jobs/new', requireAuth(renderJobSubmit));
  router.add('#/jobs/:id/dashboard', requireAuth(renderJobDashboard));
  router.add('#/jobs/:id', requireAuth(renderJobDetail));
  router.add('#/compare', requireAuth(renderCompare));
  router.add('#/pricing', renderPricing);

  // ── Boot ─────────────────────────────────────────────
  window.addEventListener('DOMContentLoaded', async () => {
    setTheme(getTheme());
    document.getElementById('theme-toggle')?.addEventListener('click', toggleTheme);

    // Restore user session
    if (api.isLoggedIn()) {
      try {
        window.__ddUser = await api.get('/api/auth/me');
      } catch (e) {
        api.clearToken();
      }
    }
    updateUserMenu();

    // Load PortOne config for frontend SDK
    try {
      const cfg = await api.get('/api/subscription/config');
      window.__PORTONE_STORE_ID = cfg.store_id || '';
      window.__PORTONE_CHANNEL_KEY = cfg.channel_key || '';
    } catch (e) {
      // PortOne not configured — pricing page will show error
    }

    window.addEventListener('hashchange', () => {
      updateSidebar();
      clearJobsPolling();
    });

    router.init();
    updateSidebar();
  });
})();
