/** The portal's own skin, over the u2 tokens so it follows the platform theme. */
const CSS = `
.pp-root { display: flex; flex-direction: row; align-items: stretch; justify-content: flex-start; gap: 12px; flex: 1 1 auto; align-self: stretch; height: 100%; width: 100%; min-width: 0; min-height: 0; padding: 8px; box-sizing: border-box; overflow: hidden; }
.pp-root > * { align-self: stretch; }
.pp-portfolios { flex: 0 0 190px; justify-content: flex-start; gap: 6px; overflow-y: auto; }
.pp-board { flex: 1 1 0; min-width: 0; gap: 8px; overflow: hidden; justify-content: flex-start; }
.pp-detail { flex: 0 0 300px; justify-content: flex-start; gap: 6px; overflow-y: auto; border-left: 1px solid var(--u2-border, #e0e0e0); padding-left: 12px; }
.pp-heading { font-weight: 600; margin-top: 8px; color: var(--u2-text-muted, #666); text-transform: uppercase; font-size: 11px; }
.pp-title { font-size: 16px; font-weight: 600; }
.pp-muted { color: var(--u2-text-muted, #888); font-size: 12px; }
.pp-portfolio { padding: 8px; border-radius: 6px; border: 1px solid var(--u2-border, #e0e0e0); cursor: pointer; gap: 2px; }
.pp-portfolio:hover, .pp-card:hover { background: var(--u2-hover, rgba(0,0,0,.04)); }
.pp-portfolio.pp-selected, .pp-card.pp-selected { border-color: var(--u2-accent, #2083d5); box-shadow: 0 0 0 1px var(--u2-accent, #2083d5); }
.pp-portfolio-name { font-weight: 600; }
.pp-board-header { gap: 8px; align-items: center; flex-wrap: wrap; flex: 0 0 auto; }
.pp-board-header .pp-title { white-space: nowrap; }
.pp-columns { gap: 8px; overflow-x: auto; flex: 1 1 auto; align-items: flex-start; padding-bottom: 8px; }
.pp-column { flex: 0 0 160px; justify-content: flex-start; gap: 6px; background: var(--u2-surface-2, rgba(0,0,0,.03)); border-radius: 6px; padding: 6px; max-height: 100%; overflow-y: auto; }
.pp-column-header { gap: 0; min-height: 44px; }
.pp-column-empty { flex: 0 0 56px; }
.pp-column-empty .pp-column-title { white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.pp-column-empty .pp-column-header .pp-muted { display: none; }
.pp-column-title { font-weight: 600; font-size: 12px; }
.pp-count { font-size: 11px; color: var(--u2-text-muted, #888); }
.pp-card { background: var(--u2-surface, #fff); border: 1px solid var(--u2-border, #e0e0e0); border-radius: 5px; padding: 6px; cursor: pointer; gap: 4px; }
.pp-card-title { font-size: 12px; font-weight: 500; }
.pp-card-meta { gap: 4px; align-items: center; flex-wrap: wrap; }
.pp-track { gap: 2px; flex-wrap: wrap; }
.pp-step { font-size: 11px; padding: 2px 6px; border-radius: 3px; background: var(--u2-surface-2, rgba(0,0,0,.05)); }
.pp-step-done { background: var(--u2-success-bg, #d8f0de); }
.pp-step-current { background: var(--u2-accent, #2083d5); color: #fff; font-weight: 600; }
.pp-step { cursor: pointer; }
.pp-step-viewed { outline: 2px solid var(--u2-accent, #2083d5); outline-offset: 1px; }
.pp-readiness { gap: 6px; align-items: center; font-size: 12px; }
.pp-criterion-block { gap: 2px; }
.pp-criterion-detail { padding-left: 8px; gap: 1px; }
.pp-gap { font-size: 11px; color: var(--u2-text-muted, #888); }
.pp-units { gap: 6px; }
.pp-unit-group { gap: 2px; }
.pp-unit-kind { font-size: 12px; font-weight: 600; }
.pp-unit { gap: 4px; align-items: center; flex-wrap: wrap; font-size: 12px; }
.pp-criterion { gap: 6px; align-items: center; flex-wrap: wrap; font-size: 12px; }
.pp-reason { color: var(--u2-error, #c0392b); font-size: 12px; }
.pp-actions { gap: 6px; }
.pp-history-row, .pp-lineage { font-size: 12px; }
.pp-detail textarea { width: 100%; box-sizing: border-box; }
`;

let installed = false;
export function installStyles(): void {
  if (installed)
    return;
  installed = true;
  const style = document.createElement('style');
  style.dataset.owner = 'project-portal-core';
  style.textContent = CSS;
  document.head.append(style);
}
