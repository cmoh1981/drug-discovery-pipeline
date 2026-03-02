/**
 * Auth-aware fetch wrapper + SSE streaming.
 * All API calls go through this module.
 */
const api = (() => {
  const TOKEN_KEY = 'dd_token';

  function getToken() {
    return localStorage.getItem(TOKEN_KEY);
  }

  function setToken(token) {
    localStorage.setItem(TOKEN_KEY, token);
  }

  function clearToken() {
    localStorage.removeItem(TOKEN_KEY);
  }

  function isLoggedIn() {
    return !!getToken();
  }

  async function request(method, path, body) {
    const headers = { 'Content-Type': 'application/json' };
    const token = getToken();
    if (token) {
      headers['Authorization'] = `Bearer ${token}`;
    }
    const opts = { method, headers };
    if (body !== undefined) {
      opts.body = JSON.stringify(body);
    }
    const res = await fetch(path, opts);
    if (res.status === 401) {
      clearToken();
      window.location.hash = '#/login';
      throw new Error('Session expired');
    }
    if (res.status === 204) return null;
    if (!res.ok) {
      const err = await res.json().catch(() => ({ detail: res.statusText }));
      throw new Error(err.detail || `Request failed: ${res.status}`);
    }
    return res.json();
  }

  function get(path) { return request('GET', path); }
  function post(path, body) { return request('POST', path, body); }
  function del(path) { return request('DELETE', path); }

  /**
   * SSE streaming via fetch + ReadableStream.
   * Returns an async generator yielding parsed JSON events.
   * Uses fetch instead of EventSource to support Authorization header.
   */
  async function* streamSSE(path) {
    const token = getToken();
    const headers = {};
    if (token) {
      headers['Authorization'] = `Bearer ${token}`;
    }
    const res = await fetch(path, { headers });
    if (!res.ok) {
      throw new Error(`SSE connection failed: ${res.status}`);
    }
    const reader = res.body.getReader();
    const decoder = new TextDecoder();
    let buffer = '';

    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      buffer += decoder.decode(value, { stream: true });
      const lines = buffer.split('\n');
      buffer = lines.pop() || '';
      for (const line of lines) {
        const trimmed = line.trim();
        if (trimmed.startsWith('data: ')) {
          try {
            yield JSON.parse(trimmed.slice(6));
          } catch (e) {
            // skip malformed events
          }
        }
      }
    }
  }

  return { getToken, setToken, clearToken, isLoggedIn, get, post, del, streamSSE };
})();
