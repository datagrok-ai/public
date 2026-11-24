import {praseAsync, renderDocument} from 'docx-preview';

const OPTIONS = {ignoreWidth: true, ignoreHeight: true, useBase64URL: true};

export const isSafeUrl = (url: string | null | undefined) => /^(https?:|mailto:)/i.test(url ?? '');

/** docx-preview writes a symbol's character code into innerHTML; anything but a hex code is replaced before rendering. */
function sanitizeSymbols(node: any, seen = new Set<object>()): void {
  if (node === null || typeof node !== 'object' || seen.has(node))
    return;
  seen.add(node);
  if (node.type === 'symbol' && !/^[0-9a-f]{1,6}$/i.test(node.char ?? ''))
    node.char = '25A1';
  for (const [key, value] of Object.entries(node)) {
    if (!key.startsWith('_'))
      sanitizeSymbols(value, seen);
  }
}

/** Renders into a shadow root, so the document's CSS (built from document values) stays inside the preview. */
export async function renderDocx(bytes: Uint8Array, host: HTMLElement): Promise<void> {
  const doc = await praseAsync(bytes, OPTIONS);
  sanitizeSymbols(doc);
  const container = document.createElement('div');
  container.className = 'fe-docx';
  host.append(container);
  const root = container.attachShadow({mode: 'open'});
  renderDocument(doc, root as unknown as HTMLElement, undefined, OPTIONS);
  for (const a of Array.from(root.querySelectorAll('a[href]'))) {
    const href = a.getAttribute('href')!;
    if (!href.startsWith('#') && !isSafeUrl(href))
      a.removeAttribute('href');
  }
  Object.assign(root.querySelector<HTMLElement>('.docx-wrapper')!.style, {background: 'var(--grey-1)', padding: '16px'});
}
