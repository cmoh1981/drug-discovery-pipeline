/**
 * Login page (#/login)
 */
function renderLogin() {
  const el = document.getElementById('app-content');
  el.innerHTML = `
    <div class="auth-page">
      <div class="auth-card card">
        <div class="card-body">
          <h2>Welcome Back</h2>
          <p class="auth-subtitle">Sign in to your pipeline account</p>
          <form id="login-form">
            <div class="form-group">
              <label for="login-email">Email</label>
              <input type="email" id="login-email" class="form-input" placeholder="you@example.com" required>
            </div>
            <div class="form-group">
              <label for="login-password">Password</label>
              <input type="password" id="login-password" class="form-input" placeholder="Your password" required>
            </div>
            <button type="submit" class="btn btn-primary" style="width:100%;justify-content:center" id="login-btn">Sign In</button>
          </form>
          <div class="auth-footer">
            Don't have an account? <a href="#/register">Create one</a>
          </div>
        </div>
      </div>
    </div>
  `;
  document.getElementById('login-form').addEventListener('submit', handleLogin);
}

async function handleLogin(e) {
  e.preventDefault();
  const btn = document.getElementById('login-btn');
  const email = document.getElementById('login-email').value.trim();
  const password = document.getElementById('login-password').value;
  if (!email || !password) return;

  btn.disabled = true;
  btn.textContent = 'Signing in...';
  try {
    const data = await api.post('/api/auth/login', { email, password });
    api.setToken(data.access_token);
    // Fetch user info
    const user = await api.get('/api/auth/me');
    window.__ddUser = user;
    updateUserMenu();
    showToast('Logged in successfully', 'success');
    router.navigate('#/jobs');
  } catch (err) {
    showToast(err.message || 'Invalid credentials', 'error');
  } finally {
    btn.disabled = false;
    btn.textContent = 'Sign In';
  }
}
