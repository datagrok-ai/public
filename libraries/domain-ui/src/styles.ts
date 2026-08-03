/**
 * The `domain-ui-*` stylesheet.
 *
 * It is injected once at runtime rather than imported as a `.css` file: this
 * library compiles with `tsc` alone and is consumed by packages whose webpack
 * config carries no css-loader, so a stylesheet import would break their build.
 * Colours are design tokens only (`core/client/xamgle/web/datagrok.css` `:root`)
 * — the in-grid CELL colours are canvas ints and live in `domain-grid.ts`.
 *
 * @module styles
 */

const STYLE_ID = 'domain-ui-styles';

// The containers themselves are already flex boxes (`d4-flex-col` / `d4-flex-row`
// from ui.divV / ui.divH) — everything below is additive.
const DOMAIN_UI_CSS = `
.domain-ui-grid {
  height: 100%;
}

.domain-ui-grid-toolbar {
  flex: 0 0 auto;
  align-items: center;
  gap: 4px;
  padding: 4px 8px;
  background-color: var(--grey-1);
  border-bottom: 1px solid var(--grey-2);
}

.domain-ui-save-bar {
  align-items: center;
  gap: 4px;
  margin-left: auto;
}

.domain-ui-edit-count {
  color: var(--text-color-light);
  font-size: 12px;
  white-space: nowrap;
  padding-right: 4px;
}

.domain-ui-grid-toolbar button[disabled] {
  color: var(--grey-3);
  cursor: default;
}
`;

/** Adds the stylesheet to the document once — safe to call per component. */
export function applyDomainUiStyles(): void {
  if (typeof document === 'undefined' || document.getElementById(STYLE_ID) != null)
    return;
  const style = document.createElement('style');
  style.id = STYLE_ID;
  style.textContent = DOMAIN_UI_CSS;
  document.head.appendChild(style);
}
