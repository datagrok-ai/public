/**
 * Registry for the "Workspace preview" host — the area below the Spotlight widget
 * in the Welcome view where the Workspace tab renders a live preview of the
 * selected pinned entity (for a query / function, the result dataframe).
 *
 * The Welcome view owns the DOM element and registers it here; the Workspace tab
 * (nested inside the Spotlight widget inside the Welcome view) pulls it out and
 * writes into it. At most one host is active at a time — the latest registration
 * wins.
 */

let host: HTMLElement | null = null;
let customizeEl: HTMLElement | null = null;

export function setWorkspacePreviewHost(el: HTMLElement | null): void {
  host = el;
}

export function setWorkspaceCustomizeEl(el: HTMLElement | null): void {
  customizeEl = el;
}

export function getWorkspacePreviewHost(): HTMLElement | null {
  return host;
}

/** On the panel while a preview has the home page to itself: the other widgets step aside, so that the preview is
 * right under the Spotlight and takes the rest of the page. */
const FOCUS_CLASS = 'power-pack-preview-focus';

export function clearWorkspacePreview(): void {
  if (!host) return;
  host.innerHTML = '';
  host.style.display = 'none';
  host.parentElement?.classList.remove(FOCUS_CLASS);
  if (customizeEl) customizeEl.style.display = '';
}

export function showWorkspacePreview(content: HTMLElement, title?: string, onOpen?: () => void): void {
  if (!host)
    return;
  host.innerHTML = '';
  host.parentElement?.classList.remove(FOCUS_CLASS);
  if (title || onOpen) {
    const header = document.createElement('div');
    header.className = 'pp-workspace-preview-header';
    if (title) {
      const titleEl = document.createElement('div');
      titleEl.className = 'pp-workspace-preview-title';
      titleEl.textContent = title;
      header.appendChild(titleEl);
    }
    if (onOpen) {
      const btn = document.createElement('button');
      btn.className = 'pp-workspace-preview-open-btn';
      btn.textContent = 'Open';
      btn.addEventListener('click', (e) => {
        e.stopPropagation();
        onOpen();
      });
      header.appendChild(btn);
    }
    host.appendChild(header);
  }
  host.appendChild(content);
  host.style.display = '';
  if (customizeEl) customizeEl.style.display = 'none';
}

/** A package Spotlight tab's results: the other widgets step aside until the Spotlight's tab changes. */
export function showFocusedPreview(content: HTMLElement, title?: string): void {
  showWorkspacePreview(content, title);
  host?.parentElement?.classList.add(FOCUS_CLASS);
}
