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

/* The grid host inside it: a canvas grid has no intrinsic size, so it needs the
   space the toolbar leaves — in BOTH directions. */
.domain-ui-grid > .ui-box {
  flex-grow: 1;
  min-height: 0;
  width: 100%;
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

.domain-ui-entity-list,
.domain-ui-entity-page {
  height: 100%;
}

/* The list's content area: whatever the toolbar leaves. min-height:0 is what
   lets a canvas grid inside it shrink instead of overflowing the page. */
.domain-ui-list-body {
  display: flex;
  flex-direction: column;
  flex-grow: 1;
  min-height: 0;
}

/* The app's status-bar slot: the pending-change count and its actions, on the
   one line the platform gives a view at the bottom of the window. */
.domain-ui-status-bar {
  display: flex;
  align-items: center;
  gap: 8px;
  color: var(--text-color-light);
  font-size: 12px;
  white-space: nowrap;
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
