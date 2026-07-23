import * as fs from 'node:fs';
import * as path from 'node:path';

// Generates workspace/help/INDEX.md — one line per published doc page, built from frontmatter
// (title / description / keywords, with first-H1 fallback). Grounding a help answer then costs
// one Read of the index + one Read of the right page instead of an unbounded grep chain
// (measured at 8-33 model round-trips per help turn before this existed — see docs/LATENCY.md).
// Regenerated at startup and after each workspace git sync; INDEX.md is untracked, so
// `git reset --hard` leaves it alone.

const SKIP_DIRS = new Set(['_internal', '_access', '_stashed', '_deprecated', '_stash', 'uploads', 'img']);

interface PageEntry {
  rel: string;
  title: string;
  description: string;
  keywords: string[];
}

function parseFrontmatter(text: string): {title?: string; description?: string; keywords: string[]} {
  const out: {title?: string; description?: string; keywords: string[]} = {keywords: []};
  const lines = text.split(/\r?\n/);
  if (!/^---\s*$/.test(lines[0] ?? ''))
    return out;
  let inKeywords = false;
  for (let i = 1; i < lines.length; i++) {
    const line = lines[i];
    if (/^---\s*$/.test(line))
      break;
    const kw = line.match(/^\s+-\s+(.+)$/);
    if (inKeywords && kw) {
      out.keywords.push(kw[1].trim().replace(/^["']|["']$/g, ''));
      continue;
    }
    inKeywords = false;
    const field = line.match(/^(\w[\w-]*):\s*(.*)$/);
    if (!field)
      continue;
    const [, key, value] = field;
    if (key === 'keywords')
      inKeywords = value.trim() === '' || value.trim() === '|';
    else if (key === 'title')
      out.title = value.trim().replace(/^["']|["']$/g, '');
    else if (key === 'description')
      out.description = value.trim().replace(/^["']|["']$/g, '');
  }
  return out;
}

function firstH1(text: string): string | undefined {
  return text.match(/^#\s+(.+)$/m)?.[1]?.trim();
}

function collectPages(helpDir: string): PageEntry[] {
  const entries: PageEntry[] = [];
  const walk = (dir: string): void => {
    for (const name of fs.readdirSync(dir, {withFileTypes: true})) {
      const full = path.join(dir, name.name);
      if (name.isDirectory()) {
        if (!SKIP_DIRS.has(name.name))
          walk(full);
        continue;
      }
      // -test.md files are internal QA scenarios — indexing them would draw help answers
      // to test checklists instead of the real page (they stay greppable on disk).
      if (!/\.(md|mdx)$/.test(name.name) || name.name.startsWith('_') ||
          /-test\.md$/.test(name.name) || name.name === 'INDEX.md' || name.name === 'CLAUDE.md')
        continue;
      let text: string;
      try {
        text = fs.readFileSync(full, 'utf8');
      } catch {
        continue;
      }
      const fm = parseFrontmatter(text);
      const title = fm.title ?? firstH1(text) ?? name.name.replace(/\.(md|mdx)$/, '');
      entries.push({
        rel: path.relative(helpDir, full).split(path.sep).join('/'),
        title,
        description: fm.description ?? '',
        keywords: fm.keywords,
      });
    }
  };
  walk(helpDir);
  return entries.sort((a, b) => a.rel.localeCompare(b.rel));
}

export function buildHelpIndex(workspaceDir: string): void {
  const helpDir = path.join(workspaceDir, 'help');
  if (!fs.existsSync(helpDir))
    return;
  const pages = collectPages(helpDir);
  const lines: string[] = [
    '# Help documentation index',
    '',
    'One line per page: `path — Title — description [keywords]`. To answer a platform question:',
    'find the relevant page below (titles + keywords), then Read `workspace/help/<path>`.',
    'If no page covers the topic, the docs do not cover it — say so instead of searching further.',
    '',
  ];
  let section = '';
  for (const p of pages) {
    const top = p.rel.split('/')[0];
    if (top !== section) {
      section = top;
      lines.push(`## ${top}`, '');
    }
    const desc = p.description ? ` — ${p.description}` : '';
    const kw = p.keywords.length ? ` [${p.keywords.join(', ')}]` : '';
    lines.push(`- ${p.rel} — ${p.title}${desc}${kw}`);
  }
  lines.push('');
  fs.writeFileSync(path.join(helpDir, 'INDEX.md'), lines.join('\n'));
  console.log(`help-index: ${pages.length} pages indexed`);
}
