/** The portal's own skin, over the u2 tokens so it follows the platform theme. */
const CSS = `
.pp-msg-nosend .u2-msg-send { display: none; }
.pp-root { display: flex; flex-direction: row; align-items: stretch; justify-content: flex-start; gap: 12px; flex: 1 1 auto; align-self: stretch; height: 100%; width: 100%; min-width: 0; min-height: 0; padding: 8px; box-sizing: border-box; overflow: hidden; }
.pp-root > * { align-self: stretch; }
.pp-portfolios { flex: 1 1 0; max-width: 720px; justify-content: flex-start; gap: 6px; overflow-y: auto; }
.pp-portfolio-cards { flex-wrap: wrap; gap: 8px; }
.pp-portfolio-cards > .pp-portfolio { flex: 0 0 220px; }
.pp-board { flex: 1 1 0; min-width: 0; gap: 8px; overflow: hidden; justify-content: flex-start; }
.pp-detail { justify-content: flex-start; gap: 6px; padding: 4px 8px; box-sizing: border-box; }
.pp-heading { font-weight: 600; margin-top: 8px; color: var(--u2-text-muted, #666); text-transform: uppercase; font-size: 11px; }
.pp-title { font-size: 16px; font-weight: 600; }
.pp-icon.grok-icon { color: var(--u2-text-muted, #888); cursor: default; margin-right: 6px; font-size: inherit; }
.pp-icon-line { align-items: center; gap: 0; }
.pp-icon-line .pp-icon { flex: 0 0 auto; }
.pp-muted { color: var(--u2-text-muted, #888); font-size: 12px; }
.pp-question { font-style: italic; font-size: 12px; }
.pp-portfolio { padding: 8px; border-radius: 6px; border: 1px solid var(--u2-border, #e0e0e0); cursor: pointer; gap: 2px; }
.pp-portfolio:hover, .pp-card:hover { background: var(--u2-hover, rgba(0,0,0,.04)); }
.pp-portfolio.pp-selected, .pp-card.pp-selected { border-color: var(--u2-accent, #2083d5); box-shadow: 0 0 0 1px var(--u2-accent, #2083d5); }
.pp-portfolio-name { font-weight: 600; }
.pp-search-box { gap: 4px; margin-bottom: 6px; }
.pp-search { flex: 0 0 auto; width: 100%; box-sizing: border-box; }
.pp-search-results, .pp-recent { gap: 3px; }
.pp-search-hit { padding: 4px 6px; border-radius: 4px; border: 1px solid var(--u2-border, #e0e0e0); cursor: pointer; gap: 1px; }
.pp-search-hit:hover { background: var(--u2-hover, rgba(0,0,0,.04)); }
.pp-search-hit .pp-card-title { font-size: 12px; }
.pp-search-hit .pp-muted { font-size: 11px; }
.pp-board-header { gap: 8px; align-items: center; flex-wrap: wrap; flex: 0 0 auto; }
.pp-access-notice { gap: 8px; align-items: center; flex: 0 0 auto; padding: 6px 8px; margin: 4px 0;
  border: 1px solid var(--u2-border, #ddd); border-radius: 4px; }
.pp-access-text { font-size: 12px; }
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
.pp-attention-item { gap: 1px; font-size: 12px; }
.pp-criterion { gap: 6px; align-items: center; flex-wrap: wrap; font-size: 12px; }
.pp-reason { color: var(--u2-error, #c0392b); font-size: 12px; }
.pp-post-body { white-space: pre-wrap; font-size: 12px; }
.pp-actions { gap: 6px; }
.pp-data-actions { flex-wrap: wrap; }
.pp-history-row, .pp-lineage { font-size: 12px; }
.pp-detail textarea { width: 100%; box-sizing: border-box; }
/* a long gate panel scrolls; its rows must not shrink (the rationale box collapsed to a line) */
.pp-detail > * { flex-shrink: 0; }
.pp-runs { gap: 4px; }
.pp-run { gap: 1px; font-size: 12px; padding: 4px 0; border-bottom: 1px solid var(--u2-border, #eee); }
.pp-run-meta { gap: 4px; align-items: center; flex-wrap: wrap; }
.pp-link { color: var(--u2-accent, #2083d5); cursor: pointer; }
.pp-link:hover { text-decoration: underline; }
.pp-model-card { gap: 4px; font-size: 12px; max-width: 520px; }
.pp-model-card > * { flex-shrink: 0; }
.pp-card-text { white-space: pre-wrap; }
.pp-card-table { border-collapse: collapse; font-size: 12px; }
.pp-table.u2-table td { white-space: normal; height: auto; padding-top: 3px; padding-bottom: 3px; }
.pp-card-table th { text-align: left; font-weight: 600; color: var(--u2-text-muted, #666); padding: 2px 10px 2px 0; }
.pp-card-table td { padding: 2px 10px 2px 0; border-top: 1px solid var(--u2-border, #eee); }
.pp-card-facts { gap: 2px; margin-top: 6px; }
.pp-card-fact { gap: 6px; }
.pp-card-fact > .pp-muted { flex: 0 0 110px; }
.pp-save-run { gap: 8px; min-width: 400px; }
/* u2 pins its editors to the platform's fixed widths; in the package's forms and filter rows they take the room there is */
.pp-save-run > .u2-input-root:not([data-u2='bool-input']) { width: 100%; }
.pp-save-run > .u2-input-root:not([data-u2='bool-input']) > .u2-input-box,
.pp-fill > .u2-input-root > .u2-input-box { flex: 1 1 auto; min-width: 0; }
/* as specific as u2's own text-box pin (.u2-input-root[data-u2='text-input'] > .u2-input-box > input.u2-input-editor) */
.pp-save-run > .u2-input-root[data-u2]:not([data-u2='bool-input']) > .u2-input-box > .u2-input-editor,
.pp-fill > .u2-input-root[data-u2] > .u2-input-box > .u2-input-editor { width: 100%; box-sizing: border-box; }
.pp-fill > .u2-input-root { flex: 1 1 0; min-width: 0; }
[data-u2-name="decision-results"] > .pp-unit { align-items: flex-start; }
.pp-field { gap: 8px; align-items: center; }
.pp-field-label { flex: 0 0 110px; font-size: 12px; color: var(--u2-text-muted, #666); }
.pp-select, .pp-text { flex: 1 1 auto; min-width: 0; height: 26px; font-size: 13px; }
.pp-warnings { gap: 2px; }
.pp-choose-run { gap: 8px; min-width: 560px; max-height: 70vh; overflow-y: auto; }
.pp-choose-run > * { flex-shrink: 0; }
.pp-choose-run textarea { width: 100%; box-sizing: border-box; }
.pp-candidate { display: flex; gap: 8px; align-items: flex-start; cursor: pointer; }
.pp-candidate input { margin-top: 6px; }
.pp-compare-box { overflow-x: auto; }
.pp-compare tr.pp-differs td { font-weight: 600; }
.pp-criterion-run { gap: 4px; align-items: center; flex-wrap: wrap; font-size: 12px; }
.pp-preferences { gap: 8px; align-items: center; margin-left: auto; }
.pp-preference { gap: 4px; align-items: center; }
.pp-preference .u2-input-root { flex-shrink: 1; min-width: 0; padding: 0; }
.pp-preference select.u2-input-editor { width: auto; max-width: 180px; }
.pp-select-small { height: 22px; font-size: 12px; flex: 0 0 auto; }
.pp-quick-cards { gap: 6px; flex-wrap: wrap; }
.pp-quick { flex: 0 0 132px; gap: 2px; padding: 6px; border: 1px solid var(--u2-border, #e0e0e0); border-radius: 5px; cursor: pointer; font-size: 12px; }
.pp-quick:hover { background: var(--u2-hover, rgba(0,0,0,.04)); }
.pp-quick-label { font-weight: 600; color: var(--u2-accent, #2083d5); }
.pp-quick-off { cursor: default; opacity: .5; }
.pp-quick-off:hover { background: none; }
.pp-role-form { align-items: center; flex-wrap: wrap; }
.pp-role-form .pp-text { flex: 1 1 120px; height: 22px; font-size: 12px; }
.pp-evidence-picker { gap: 2px; font-size: 12px; }
.pp-evidence-pick { display: flex; align-items: center; gap: 4px; cursor: pointer; }
.pp-evidence { padding-left: 12px; gap: 1px; }
.pp-provenance { gap: 4px; min-width: 520px; max-height: 70vh; overflow-y: auto; }
.pp-provenance > * { flex-shrink: 0; }
.pp-provenance-tree { gap: 2px; }
.pp-document-card { gap: 2px; font-size: 12px; padding: 4px 0; border-bottom: 1px solid var(--u2-border, #eee); }
.pp-my-programs { width: 100%; height: 100%; overflow-y: auto; }
.pp-my-programs .pp-tiles { display: grid; grid-template-columns: repeat(auto-fill, minmax(210px, 1fr)); gap: 8px; align-items: start; padding: 2px; }
.pp-tile { gap: 3px; }
.pp-tile-stage { font-size: 12px; }
.pp-tile-track { gap: 2px; flex-wrap: nowrap; }
.pp-tile-seg { flex: 1 1 0; height: 4px; border-radius: 2px; background: var(--u2-surface-2, rgba(0,0,0,.08)); }
.pp-tile-seg-done { background: var(--u2-success, #5cb85c); }
.pp-tile-seg-current { background: var(--u2-accent, #2083d5); }
.pp-pinned-mark { color: var(--u2-accent, #2083d5); font-size: 11px; }
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
