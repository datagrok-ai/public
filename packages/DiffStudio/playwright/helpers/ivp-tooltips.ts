import { Page } from '@playwright/test';

/**
 * Read an IVP file from `System:AppData/DiffStudio/library/<file>` and parse its
 * `[tooltip in brackets]` annotations into a `caption → tooltip` map.
 *
 * IVP line format (relevant parts):
 *   <var> = <value> {caption: <visible name>; category: ...; ...} [<tooltip text>]
 *   <var> = <value> {category: ...; ...} [<tooltip text>]   // variable name doubles as caption
 *
 * Tooltips that are bracketed at the end of the line become the UI tooltip when the user
 * hovers the input label. Inputs without a `[<tooltip>]` get the platform default.
 *
 * Reading the file uses `grok.dapi.files.readAsText` — a single API call to access platform
 * file state. The corresponding UI gesture (clicking the file in Browse) opens a preview view
 * and would destroy any model state in the active view, so the API read is the appropriate
 * fallback per the spec.
 */
export async function readIvpTooltips(page: Page, ivpPath: string): Promise<Map<string, string>> {
  const text = await page.evaluate(async (p) => {
    try {
      return await (window as any).grok.dapi.files.readAsText(p);
    } catch {
      return null;
    }
  }, ivpPath);
  const out = new Map<string, string>();
  if (!text) return out;
  for (const raw of (text as string).split(/\r?\n/)) {
    const line = raw.trim();
    if (!line || line.startsWith('#')) continue;
    // Match: <name> = <value> [optional {annotations}] [tooltip]
    const m = line.match(/^([A-Za-z_][\w-]*)\s*=\s*[^{]*?(\{[^}]*\})?\s*\[([^\]]+)\]\s*$/);
    if (!m) continue;
    const varName = m[1];
    const annots = m[2] ?? '';
    const tooltip = m[3].trim();
    // Caption: from {caption: ...} if present, otherwise the variable name itself.
    const capMatch = annots.match(/caption\s*:\s*([^;}]+?)\s*(?:;|}|$)/);
    const caption = capMatch ? capMatch[1].trim() : varName;
    out.set(caption, tooltip);
  }
  return out;
}

/** Map a caption to the safeName used in `input-host-<safeName>` (`[:_; *\[\]{}|]` → `-`). */
export function captionToSafeName(caption: string): string {
  return caption.replace(/[:_; *\\\[\]{}|]/g, '-');
}
