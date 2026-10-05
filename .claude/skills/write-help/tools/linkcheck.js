// Checks relative links and #anchors in help pages.
// Usage: node linkcheck.js [page.md ...]   (paths relative to the current directory)
// Without arguments: the .md/.mdx files changed on this branch since it left origin/master, plus new untracked ones.
const fs = require('fs');
const path = require('path');
const {execSync} = require('child_process');

const slug = (h) => {
  const id = h.match(/\{#([^}]+)\}\s*$/);
  if (id) return id[1];
  return h.replace(/<[^>]+>/g, '').replace(/`/g, '').replace(/\*\*/g, '').replace(/\[([^\]]*)\]\([^)]*\)/g, '$1')
    .trim().toLowerCase().replace(/[^\w\s-]/g, '').replace(/\s/g, '-');
};
// Docusaurus suffixes repeated heading ids with -1, -2, ...; explicit id="..." attributes count too.
const anchors = (file) => {
  const text = fs.readFileSync(file, 'utf8').replace(/```[\s\S]*?```/g, '');
  const result = new Set();
  const seen = {};
  for (const line of text.split(/\r?\n/).filter((l) => /^#{2,6}\s/.test(l))) {
    const s = slug(line.replace(/^#+\s*/, ''));
    result.add(seen[s] ? `${s}-${seen[s]}` : s);
    seen[s] = (seen[s] ?? 0) + 1;
  }
  for (const m of text.matchAll(/\bid=["']([^"']+)["']/g)) result.add(m[1]);
  return result;
};
// Docusaurus resolves `page`, `page.md`, `page.mdx`, and `dir/` (its index or same-named page).
const resolveDoc = (p) => {
  for (const c of [p, `${p}.md`, `${p}.mdx`, path.join(p, 'index.md'), path.join(p, `${path.basename(p)}.md`), path.join(p, `${path.basename(p)}.mdx`)])
    if (fs.existsSync(c) && fs.statSync(c).isFile()) return c;
  return fs.existsSync(p) ? p : null;
};
const git = (cmd) => execSync(cmd).toString().split(/\r?\n/).filter(Boolean);

let files = process.argv.slice(2);
if (files.length === 0)
  files = [...new Set([
    ...git('git diff --name-only --relative --merge-base origin/master -- "*.md" "*.mdx"'),
    ...git('git ls-files --others --exclude-standard -- "*.md" "*.mdx"'),
  ])].filter((f) => fs.existsSync(f));

let bad = 0;
for (const f of files.filter((f) => !fs.existsSync(f))) { console.log(`${f}: file not found`); bad++; }
files = files.filter((f) => fs.existsSync(f));
for (const file of files) {
  const text = fs.readFileSync(file, 'utf8').replace(/```[\s\S]*?```/g, '').replace(/<!--[\s\S]*?-->/g, '');
  const re = /\]\(([^)\s]+)(?:\s+"[^"]*")?\)|require\(['"]([^'"]+)['"]\)|\bsrc=["']([^"'{]+)["']/g;
  let m;
  while ((m = re.exec(text))) {
    const link = m[1] || m[2] || m[3];
    if (/^([a-z]+:|\/\/)/i.test(link) || link.startsWith('/')) continue;
    const [p, a] = link.split('#');
    const target = p ? resolveDoc(path.resolve(path.dirname(file), decodeURI(p))) : file;
    if (!target) { console.log(`${file}: missing ${link}`); bad++; continue; }
    if (a && /\.mdx?$/.test(target) && !anchors(target).has(a)) { console.log(`${file}: no anchor ${link}`); bad++; }
  }
}
console.log(bad ? `${bad} problems in ${files.length} files` : `all links ok in ${files.length} files`);
process.exitCode = bad ? 1 : 0;
