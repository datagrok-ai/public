import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {convertIvpToLatex, ConvertOptions} from 'diff-grok';

import {LINK} from '../ui-constants';

type Format = 'latex' | 'markdown';

/* Mirrors diff-grok/src/latex-export/output/tex-formatter.ts (not re-exported at top level).
   Packages match what buildLatexDocument emits — align (amsmath), symbols (amssymb),
   booktabs rules (booktabs). Re-check on diff-grok version bumps. */
function wrapTex(content: string): string {
  return [
    '\\documentclass{article}',
    '\\usepackage{amsmath}',
    '\\usepackage{amssymb}',
    '\\usepackage{booktabs}',
    '',
    '\\begin{document}',
    '',
    content,
    '',
    '\\end{document}',
  ].join('\n');
}

function isoDate(): string {
  return new Date().toISOString().slice(0, 10);
}

function extForFormat(f: Format): string {
  return f === 'latex' ? 'tex' : 'md';
}

function titleForFormat(f: Format): string {
  return f === 'latex' ? 'Export as LaTeX' : 'Export as Markdown';
}

const MULT_CDOT = 'a · b';
const MULT_JUXT = 'ab';

export function showExportDialog(
  ivpText: string,
  modelName: string,
  initialFormat: Format,
): void {
  const formatInput = ui.input.choice<string>('Format', {
    value: initialFormat,
    items: ['markdown', 'latex'],
    tooltipText: 'Output format: LaTeX (.tex) or Markdown with LaTeX math (.md).',
  });

  const metaInput = ui.input.bool('Title & description', {
    value: true,
    tooltipText: 'Include metadata sections: model name, description, and comment.',
  });
  const initsInput = ui.input.bool('Initial conditions', {
    value: true,
    tooltipText: 'Include the initial conditions table (variable values at the start of the simulation).',
  });
  const paramsInput = ui.input.bool('Parameters', {
    value: true,
    tooltipText: 'Include the parameters table (values that generate UI controls for model exploration).',
  });
  const constsInput = ui.input.bool('Constants', {
    value: true,
    tooltipText: 'Include the constants table (fixed values used in equations).',
  });

  const compactInput = ui.input.bool('Compact layout', {
    value: true,
    tooltipText: 'Compact mode: no section headings, inline initial conditions. ' +
      'Suitable for small models (up to three equations).',
  });
  const cdotInput = ui.input.choice<string>('Multiplication', {
    value: MULT_CDOT,
    items: [MULT_CDOT, MULT_JUXT],
    tooltipText: 'Symbol used for multiplication. "a · b" renders \\cdot; ' +
      '"ab" uses juxtaposition, common in physics papers.',
  });
  const standaloneInput = ui.input.bool('Standalone document', {
    value: true,
    tooltipText: 'Wrap the output with \\documentclass{article} and the required \\usepackage lines ' +
      'so the file compiles with pdflatex out of the box. ' +
      'Turn off to get a body fragment for pasting into an existing LaTeX document. Ignored for Markdown.',
  });

  const groupHeader = (text: string): HTMLElement => {
    const h = ui.element('div');
    h.textContent = text;
    h.style.cssText =
      'margin: 12px 0 4px 0; font-size: 11px; font-weight: 600; ' +
      'text-transform: uppercase; letter-spacing: 0.4px; color: var(--grey-5);';
    return h;
  };

  const optionsPanel = ui.divV([
    groupHeader('Format'),
    formatInput.root,
    groupHeader('Include'),
    metaInput.root,
    initsInput.root,
    paramsInput.root,
    constsInput.root,
    groupHeader('Style'),
    compactInput.root,
    cdotInput.root,
    standaloneInput.root,
  ]);
  optionsPanel.style.cssText =
    'flex: 0 0 260px; width: 260px; box-sizing: border-box; ' +
    'padding: 4px 16px 4px 4px; overflow-y: auto; overflow-x: hidden; ' +
    'border-right: 1px solid var(--grey-2);';

  const previewEl = ui.element('pre');
  previewEl.style.cssText =
    'font-family: var(--grok-font-family-monospace, monospace); ' +
    'font-size: 12px; margin: 0; padding: 12px 16px; ' +
    'white-space: pre-wrap; word-break: break-word; ' +
    'box-sizing: border-box; height: 100%; width: 100%; overflow: auto; ' +
    'background: var(--grey-1); border-radius: 4px;';

  const copyIcon = ui.icons.copy(async () => {
    if (!lastResult) return;
    try {
      await navigator.clipboard.writeText(lastResult);
      grok.shell.info('Copied to clipboard');
    } catch {
      grok.shell.error('Clipboard access denied');
    }
  }, 'Copy the generated source to the clipboard.');
  copyIcon.style.cssText =
    'position: absolute; top: 10px; right: 14px; z-index: 2; ' +
    'padding: 4px 6px; cursor: pointer; ' +
    'color: var(--grey-5); background: var(--grey-1); border-radius: 3px;';

  const previewPanel = ui.div([previewEl, copyIcon]);
  previewPanel.style.cssText =
    'flex: 1 1 auto; min-width: 0; height: 100%; ' +
    'display: flex; padding-left: 12px; box-sizing: border-box; ' +
    'position: relative;';

  const fileNameInput = ui.input.string('File name', {
    value: `${sanitizeBase(modelName)}-${isoDate()}.${extForFormat(initialFormat)}`,
    nullable: false,
    tooltipText: 'Name of the file created by Download. ' +
      'The extension updates automatically with the selected format. Required.',
  });

  fileNameInput.root.style.flex = '1 1 auto';
  const actionsRow = ui.divH([fileNameInput.root], {style: {
    alignItems: 'center', padding: '10px 4px 0 4px',
  }});

  const split = ui.divH([optionsPanel, previewPanel]);
  split.style.cssText = 'height: 480px; width: 100%; box-sizing: border-box;';

  const body = ui.divV([split, actionsRow]);
  body.style.cssText = 'width: 100%; box-sizing: border-box;';

  const dialog = ui.dialog({
    title: titleForFormat(initialFormat),
    helpUrl: LINK.LOAD_SAVE,
  }).add(body);

  let lastResult: string | null = null;

  const collectOptions = (): Partial<ConvertOptions> => ({
    format: (formatInput.value as Format) ?? 'latex',
    includeMetadata: metaInput.value,
    includeInits: initsInput.value,
    includeParameters: paramsInput.value,
    includeConstants: constsInput.value,
    compact: compactInput.value,
    useCdot: cdotInput.value === MULT_CDOT,
  });

  let hasResult = false;
  const hasFileName = (): boolean => (fileNameInput.value ?? '').trim().length > 0;

  const updateButtons = (): void => {
    copyIcon.style.display = hasResult ? '' : 'none';
    const dl = dialog.getButton('Download');
    if (dl) dl.disabled = !(hasResult && hasFileName());
  };

  const regenerate = (): void => {
    try {
      const opts = collectOptions();
      const converted = convertIvpToLatex(ivpText, opts);
      const isLatex = opts.format === 'latex';
      const result = (isLatex && standaloneInput.value) ? wrapTex(converted) : converted;
      previewEl.textContent = result;
      if (!result.trim()) {
        lastResult = null;
        hasResult = false;
      } else {
        lastResult = result;
        hasResult = true;
      }
    } catch (e) {
      previewEl.textContent = `// Error: ${(e as Error).message}`;
      lastResult = null;
      hasResult = false;
    }
    updateButtons();
  };

  fileNameInput.onChanged.subscribe(() => updateButtons());

  formatInput.onChanged.subscribe(() => regenerate());
  metaInput.onChanged.subscribe(() => regenerate());
  initsInput.onChanged.subscribe(() => regenerate());
  paramsInput.onChanged.subscribe(() => regenerate());
  constsInput.onChanged.subscribe(() => regenerate());
  compactInput.onChanged.subscribe(() => regenerate());
  cdotInput.onChanged.subscribe(() => regenerate());
  standaloneInput.onChanged.subscribe(() => regenerate());

  formatInput.onChanged.subscribe(() => {
    const fmt = (formatInput.value as Format) ?? 'latex';
    const newExt = extForFormat(fmt);
    fileNameInput.value = fileNameInput.value.replace(/\.(md|tex)$/i, `.${newExt}`);
    standaloneInput.root.style.display = (fmt === 'latex') ? '' : 'none';
    dialog.title = titleForFormat(fmt);
  });

  // hide Standalone up-front if Markdown was the initial format
  if (initialFormat === 'markdown')
    standaloneInput.root.style.display = 'none';

  dialog.addButton('Download', () => {
    if (!lastResult) return;
    DG.Utils.download(fileNameInput.value, lastResult);
  });

  regenerate();
  dialog.show({modal: true, resizable: true, center: true, width: 960, height: 640});
}

function sanitizeBase(name: string): string {
  const cleaned = (name || '').trim().replace(/[^A-Za-z0-9._-]+/g, '-').replace(/^-+|-+$/g, '');
  return cleaned || 'diff-export';
}

