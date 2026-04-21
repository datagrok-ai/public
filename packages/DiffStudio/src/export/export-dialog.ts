import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import type {EditorView} from 'codemirror';
import type {EditorState, Compartment, Extension} from '@codemirror/state';

import {convertIvpToLatex, ConvertOptions} from 'diff-grok';

import {LINK} from '../ui-constants';

//@ts-ignore
import '../../css/app-styles.css';

type Format = 'latex' | 'markdown';

/** Cached CodeMirror 6 modules; populated on first call to `getCM()` and reused for later dialogs. */
let _cm: {
  basicSetup: typeof import('codemirror')['basicSetup'];
  EditorView: typeof EditorView;
  EditorState: typeof EditorState;
  Compartment: typeof Compartment;
  markdown: typeof import('@codemirror/lang-markdown')['markdown'];
  StreamLanguage: typeof import('@codemirror/language')['StreamLanguage'];
  stex: typeof import('@codemirror/legacy-modes/mode/stex')['stex'];
} | null = null;

/** Lazily imports CodeMirror 6 plus the language modes used by the export preview
 * (Markdown, and stex via StreamLanguage for LaTeX). Results are cached. */
async function getCM() {
  if (!_cm) {
    const [{basicSetup, EditorView}, {EditorState, Compartment}, {markdown}, {StreamLanguage}, {stex}] =
      await Promise.all([
        import('codemirror'),
        import('@codemirror/state'),
        import('@codemirror/lang-markdown'),
        import('@codemirror/language'),
        import('@codemirror/legacy-modes/mode/stex'),
      ]);
    _cm = {basicSetup, EditorView, EditorState, Compartment, markdown, StreamLanguage, stex};
  }
  return _cm;
} // getCM

/** Wraps a LaTeX body fragment into a standalone compilable document with the preamble
 * used by `diff-grok`'s `buildLatexDocument` (amsmath, amssymb, booktabs). */
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

/** Returns today's date as an ISO-8601 `YYYY-MM-DD` string, used in default export file names. */
function isoDate(): string {
  return new Date().toISOString().slice(0, 10);
}

/** Returns the file-name extension that matches the chosen export format (`tex` or `md`). */
function extForFormat(f: Format): string {
  return f === 'latex' ? 'tex' : 'md';
}

/** Returns the human-readable dialog title for the chosen export format. */
function titleForFormat(f: Format): string {
  return f === 'latex' ? 'Export as LaTeX' : 'Export as Markdown';
}

const MULT_CDOT = 'a · b';
const MULT_JUXT = 'ab';

/** Opens the interactive export dialog: lets the user tune LaTeX/Markdown conversion options,
 * live-preview the result, copy it to the clipboard, or download it as a file. */
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

  const compactInput = ui.input.bool('Compact', {
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
  const standaloneInput = ui.input.bool('Standalone', {
    value: true,
    tooltipText: 'Wrap the output with \\documentclass{article} and the required \\usepackage lines ' +
      'so the file compiles with pdflatex out of the box. ' +
      'Turn off to get a body fragment for pasting into an existing LaTeX document. Ignored for Markdown.',
  });

  const extraEls: HTMLElement[] = [
    metaInput.root,
    initsInput.root,
    paramsInput.root,
    constsInput.root,
    cdotInput.root,
  ];

  let extrasVisible = false;
  /** Shows or hides the secondary export options toggled by the settings icon. */
  const applyExtrasVisibility = (): void => {
    for (const el of extraEls)
      el.style.display = extrasVisible ? '' : 'none';
  };

  const settingsIcon = ui.iconFA('cog', () => {
    extrasVisible = !extrasVisible;
    applyExtrasVisibility();
  }, 'Customize export');
  settingsIcon.classList.add('diff-studio-export-settings-icon');

  const formEl = ui.form([
    formatInput,
    compactInput,
    standaloneInput,
    metaInput,
    initsInput,
    paramsInput,
    constsInput,
    cdotInput,
  ]);
  formEl.classList.add('diff-studio-export-form');

  const optionsPanel = ui.divH([settingsIcon, formEl]);
  optionsPanel.classList.add('diff-studio-export-options');

  applyExtrasVisibility();

  const previewEl = ui.div();
  previewEl.classList.add('diff-studio-export-preview');

  let editorView: EditorView | null = null;
  let langCompartment: Compartment | null = null;
  let pendingDoc = '';

  /** Returns the CodeMirror language extension matching the currently selected format:
   * `stex` (via `StreamLanguage`) for LaTeX, `markdown()` for Markdown. Falls back to
   * an empty extension array when CM modules haven't finished loading yet. */
  const getLanguageExt = (): Extension => {
    if (!_cm) return [];
    const fmt = (formatInput.value as Format) ?? 'latex';
    return fmt === 'markdown' ? _cm.markdown() : _cm.StreamLanguage.define(_cm.stex);
  };

  /** Replaces the preview editor's document with `text`. If the editor is still being
   * loaded, stashes the value in `pendingDoc` so it is applied on first mount. */
  const setPreviewText = (text: string): void => {
    if (editorView) {
      editorView.dispatch({
        changes: {from: 0, to: editorView.state.doc.length, insert: text},
      });
    } else
      pendingDoc = text;
  };

  /** Reconfigures the language compartment so syntax highlighting follows the format
   * input (LaTeX vs. Markdown). No-op until the editor has been mounted. */
  const syncLanguage = (): void => {
    if (!editorView || !langCompartment) return;
    editorView.dispatch({effects: langCompartment.reconfigure(getLanguageExt())});
  };

  /** Mounts the CodeMirror preview into `previewEl` as a read-only editor with a
   * reconfigurable language compartment, then flushes any text queued while loading. */
  const initEditor = async (): Promise<void> => {
    const cm = await getCM();
    langCompartment = new cm.Compartment();
    editorView = new cm.EditorView({
      doc: pendingDoc,
      extensions: [
        cm.basicSetup,
        langCompartment.of(getLanguageExt()),
        cm.EditorView.editable.of(false),
        cm.EditorState.readOnly.of(true),
      ],
      parent: previewEl,
    });
    pendingDoc = '';
  };
  initEditor();

  const copyFeedback = ui.element('span');
  copyFeedback.classList.add('diff-studio-export-copy-feedback');

  let feedbackTimer: number | null = null;
  /** Briefly shows an inline status badge next to the copy icon. Used instead of
   * `grok.shell.info/error` because the modal overlay dims shell-level balloons. */
  const showCopyFeedback = (text: string, isError: boolean): void => {
    copyFeedback.textContent = text;
    copyFeedback.classList.toggle('diff-studio-export-copy-feedback--error', isError);
    copyFeedback.classList.add('diff-studio-export-copy-feedback--visible');
    if (feedbackTimer !== null) window.clearTimeout(feedbackTimer);
    feedbackTimer = window.setTimeout(() => {
      copyFeedback.classList.remove('diff-studio-export-copy-feedback--visible');
      feedbackTimer = null;
    }, 1500);
  };

  const copyIcon = ui.icons.copy(async () => {
    if (!lastResult) return;
    try {
      await navigator.clipboard.writeText(lastResult);
      showCopyFeedback('Copied', false);
    } catch {
      showCopyFeedback('Copy failed', true);
    }
  }, 'Copy the generated source to the clipboard.');
  copyIcon.classList.add('diff-studio-export-copy-icon');

  const previewPanel = ui.div([previewEl, copyFeedback, copyIcon]);
  previewPanel.classList.add('diff-studio-export-preview-host');

  const fileNameInput = ui.input.string('File name', {
    value: `${sanitizeBase(modelName)}-${isoDate()}.${extForFormat(initialFormat)}`,
    nullable: false,
    tooltipText: 'Name of the file created by Download. ' +
      'The extension updates automatically with the selected format.',
  });
  fileNameInput.root.classList.add('diff-studio-export-filename');

  const actionsRow = ui.divH([fileNameInput.root]);
  actionsRow.classList.add('diff-studio-export-actions-row');

  const split = ui.divH([optionsPanel, previewPanel]);
  split.classList.add('diff-studio-export-split');

  const body = ui.divV([split, actionsRow]);
  body.classList.add('diff-studio-export-body');

  const dialog = ui.dialog({
    title: titleForFormat(initialFormat),
    helpUrl: LINK.DIF_STUDIO_REL,
  }).add(body);

  let lastResult: string | null = null;

  /** Collects the current dialog inputs into a `ConvertOptions` object for `convertIvpToLatex`. */
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
  /** Returns `true` when the file-name input contains a non-blank value. */
  const hasFileName = (): boolean => (fileNameInput.value ?? '').trim().length > 0;

  /** Syncs button/icon states with the current preview and file-name validity:
   * hides the copy icon when there is nothing to copy and disables Download accordingly. */
  const updateButtons = (): void => {
    copyIcon.style.display = hasResult ? '' : 'none';
    const dl = dialog.getButton('Download');
    if (dl) dl.disabled = !(hasResult && hasFileName());
  };

  /** Regenerates the preview from the current options; wraps LaTeX into a standalone document
   * when requested, and reports conversion errors inline in the preview pane. */
  const regenerate = (): void => {
    try {
      const opts = collectOptions();
      const converted = convertIvpToLatex(ivpText, opts);
      const isLatex = opts.format === 'latex';
      const result = (isLatex && standaloneInput.value) ? wrapTex(converted) : converted;
      setPreviewText(result);
      if (!result.trim()) {
        lastResult = null;
        hasResult = false;
      } else {
        lastResult = result;
        hasResult = true;
      }
    } catch (e) {
      setPreviewText(`// Error: ${(e as Error).message}`);
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
    syncLanguage();
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
} // showExportDialog

/** Sanitizes a model name for use as a file-name base: strips unsafe characters, trims
 * stray separators, and falls back to `'diff-export'` when nothing usable remains. */
function sanitizeBase(name: string): string {
  const cleaned = (name || '').trim().replace(/[^A-Za-z0-9._-]+/g, '-').replace(/^-+|-+$/g, '');
  return cleaned || 'diff-export';
}

