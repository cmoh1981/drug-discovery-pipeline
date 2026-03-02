/**
 * Minimal hash-based router with :param support.
 */
const router = (() => {
  const routes = [];

  function add(pattern, handler) {
    // Convert pattern like '#/jobs/:id' to regex
    const parts = pattern.replace(/:[^/]+/g, '([^/]+)');
    const paramNames = (pattern.match(/:[^/]+/g) || []).map(p => p.slice(1));
    const regex = new RegExp(`^${parts}$`);
    routes.push({ regex, paramNames, handler, pattern });
  }

  function resolve() {
    const hash = window.location.hash || '#/';
    for (const route of routes) {
      const match = hash.match(route.regex);
      if (match) {
        const params = {};
        route.paramNames.forEach((name, i) => {
          params[name] = decodeURIComponent(match[i + 1]);
        });
        route.handler(params);
        return;
      }
    }
    // Fallback to landing
    window.location.hash = '#/';
  }

  function navigate(path) {
    window.location.hash = path;
  }

  function init() {
    window.addEventListener('hashchange', resolve);
    resolve();
  }

  return { add, resolve, navigate, init };
})();
