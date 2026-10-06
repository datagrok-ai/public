import { Page } from '@playwright/test';

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

    const m = line.match(/^([A-Za-z_][\w-]*)\s*=\s*[^{]*?(\{[^}]*\})?\s*\[([^\]]+)\]\s*$/);
    if (!m) continue;
    const varName = m[1];
    const annots = m[2] ?? '';
    const tooltip = m[3].trim();

    const capMatch = annots.match(/caption\s*:\s*([^;}]+?)\s*(?:;|}|$)/);
    const caption = capMatch ? capMatch[1].trim() : varName;
    out.set(caption, tooltip);
  }
  return out;
}

export function captionToSafeName(caption: string): string {
  return caption.replace(/[:_; *\\\[\]{}|]/g, '-');
}
