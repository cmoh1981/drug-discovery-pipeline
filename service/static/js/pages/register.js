/**
 * Register page (#/register)
 */
function renderRegister() {
  const el = document.getElementById('app-content');
  el.innerHTML = `
    <div class="auth-page">
      <div class="auth-card card">
        <div class="card-body">
          <h2>Create Account</h2>
          <p class="auth-subtitle">Start discovering drug candidates</p>
          <form id="register-form">
            <div class="form-group">
              <label for="reg-email">Email</label>
              <input type="email" id="reg-email" class="form-input" placeholder="you@example.com" required>
            </div>
            <div class="form-group">
              <label for="reg-password">Password</label>
              <input type="password" id="reg-password" class="form-input" placeholder="Min 8 characters" required minlength="8">
            </div>
            <div class="form-row">
              <div class="form-group">
                <label for="reg-name">Full Name</label>
                <input type="text" id="reg-name" class="form-input" placeholder="Optional">
              </div>
              <div class="form-group">
                <label for="reg-org">Organization</label>
                <input type="text" id="reg-org" class="form-input" placeholder="Optional">
              </div>
            </div>
            <button type="submit" class="btn btn-primary" style="width:100%;justify-content:center" id="reg-btn">Create Account</button>
          </form>
          <div class="auth-footer">
            Already have an account? <a href="#/login">Sign in</a>
          </div>
        </div>
      </div>
    </div>
  `;
  document.getElementById('register-form').addEventListener('submit', handleRegister);
}

async function handleRegister(e) {
  e.preventDefault();
  const btn = document.getElementById('reg-btn');
  const email = document.getElementById('reg-email').value.trim();
  const password = document.getElementById('reg-password').value;
  const full_name = document.getElementById('reg-name').value.trim();
  const organization = document.getElementById('reg-org').value.trim();

  if (!email || password.length < 8) {
    showToast('Password must be at least 8 characters', 'error');
    return;
  }

  btn.disabled = true;
  btn.textContent = 'Creating account...';
  try {
    await api.post('/api/auth/register', { email, password, full_name, organization });
    // Auto-login
    const data = await api.post('/api/auth/login', { email, password });
    api.setToken(data.access_token);
    const user = await api.get('/api/auth/me');
    window.__ddUser = user;
    updateUserMenu();
    showToast('Account created! Welcome aboard.', 'success');
    router.navigate('#/jobs');
  } catch (err) {
    showToast(err.message || 'Registration failed', 'error');
  } finally {
    btn.disabled = false;
    btn.textContent = 'Create Account';
  }
}
