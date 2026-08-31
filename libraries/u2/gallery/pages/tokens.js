function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

export async function render(main) {
  main.append(el('h1', null, 'Design tokens'));
  const intro = el('p');
  intro.innerHTML = 'L0 of u2: <code>--dg-*</code> custom properties codifying current Datagrok ' +
    'computed styles. Single source of truth for old CSS and u2 components.';
  const toggle = el('label', 'u2-density-toggle');
  const cb = el('input');
  cb.type = 'checkbox';
  cb.checked = document.documentElement.hasAttribute('data-dg-density');
  cb.addEventListener('change', () => cb.checked ?
    document.documentElement.setAttribute('data-dg-density', 'compact') :
    document.documentElement.removeAttribute('data-dg-density'));
  toggle.append(cb, ' compact density');
  const sheet = el('div');
  main.append(intro, toggle, sheet);

  const css = await (await fetch('../css/tokens.css')).text();
  const seen = new Set();
  const tokens = [...css.matchAll(/(--dg-[\w-]+):\s*([^;]+);/g)]
    .map(([, name, value]) => ({name, value: value.trim()}))
    .filter((t) => !seen.has(t.name) && seen.add(t.name));
  const style = getComputedStyle(document.documentElement);
  const resolved = (name) => style.getPropertyValue(name).trim();
  const section = (title) => sheet.append(el('h2', null, title));
  const isColor = (t) => (t.value.startsWith('#') || t.value.startsWith('var(')) && !/solid|px/.test(t.value);

  section('Colors');
  const grid = el('div', 'u2-token-grid');
  for (const t of tokens.filter(isColor)) {
    const sw = el('div', 'u2-token-swatch');
    const chip = el('div', 'chip');
    chip.style.background = `var(${t.name})`;
    const info = el('div');
    info.append(el('div', 'name', t.name), el('div', 'value', `${t.value} → ${resolved(t.name)}`));
    sw.append(chip, info);
    grid.append(sw);
  }
  sheet.append(grid);

  section('Typography');
  for (const t of tokens.filter((t) => t.name.includes('font') || t.name.includes('line-height'))) {
    const row = el('div', 'u2-token-row');
    const demo = el('div', null, 'Grokking data 13px — the quick brown fox');
    if (t.name.includes('family')) demo.style.fontFamily = `var(${t.name})`;
    else if (t.name.includes('size')) demo.style.fontSize = `var(${t.name})`;
    else if (t.name.includes('weight')) demo.style.fontWeight = `var(${t.name})`;
    else demo.textContent = resolved(t.name);
    row.append(el('div', 'name', t.name), demo);
    sheet.append(row);
  }

  section('Spacing');
  for (const t of tokens.filter((t) => t.name.startsWith('--dg-space'))) {
    const row = el('div', 'u2-token-row');
    const bar = el('div', 'u2-token-bar');
    bar.style.width = `var(${t.name})`;
    row.append(el('div', 'name', t.name), bar, el('div', 'value', resolved(t.name)));
    sheet.append(row);
  }

  section('Control metrics (toggle density above)');
  for (const t of tokens.filter((t) => t.name.includes('height') && !t.name.includes('line'))) {
    const row = el('div', 'u2-token-row');
    row.style.alignItems = 'center';
    const demo = el('div', 'u2-metric-demo');
    demo.style.height = `var(${t.name})`;
    demo.textContent = resolved(t.name);
    row.append(el('div', 'name', t.name), demo);
    sheet.append(row);
  }

  section('Borders, shadows, layers');
  for (const t of tokens.filter((t) => /border|shadow|radius|z-/.test(t.name) && !isColor(t))) {
    const row = el('div', 'u2-token-row');
    row.append(el('div', 'name', t.name), el('div', 'value', t.value));
    sheet.append(row);
  }
}
