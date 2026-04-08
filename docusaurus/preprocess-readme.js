// Pre-processes the JS API readme to convert relative links to absolute URLs.
// TypeDoc 0.28 rewrites relative links to _media/ paths, breaking Docusaurus link resolution.
// This script creates a processed copy with absolute URLs that TypeDoc can safely consume.

const fs = require('fs');
const path = require('path');

const src = path.resolve(__dirname, '../help/develop/packages/js-api.md');
const dest = path.resolve(__dirname, 'js-api-readme.md');

let content = fs.readFileSync(src, 'utf8');

// Convert relative markdown/mdx links to absolute URLs
// ../../foo.md(x) -> https://datagrok.ai/help/foo
content = content.replace(/\(\.\.\/\.\.\/([^)]*?)\.(mdx?)(#[^)]*)?\)/g, '(https://datagrok.ai/help/$1$3)');
// ../foo.md(x) -> https://datagrok.ai/help/develop/foo (negative lookahead to avoid matching ../../)
content = content.replace(/\(\.\.\/(?!\.\.)([^)]*?)\.(mdx?)(#[^)]*)?\)/g, '(https://datagrok.ai/help/develop/$1$3)');
// ./foo.md(x) -> https://datagrok.ai/help/develop/packages/foo
content = content.replace(/\(\.\/([^)]*?)\.(mdx?)(#[^)]*)?\)/g, '(https://datagrok.ai/help/develop/packages/$1$3)');

// Convert relative image links to absolute URLs (avoids _media rewriting by TypeDoc 0.28)
// Handle optional title attribute in markdown image links, e.g. (path.png "title")
// ../../foo.png -> https://datagrok.ai/help/foo.png
content = content.replace(/\(\.\.\/\.\.\/([^)"]*?\.(png|gif|jpg|jpeg))(\s+"[^"]*")?\)/gi, '(https://datagrok.ai/help/$1$3)');
// ../foo.png -> https://datagrok.ai/help/develop/foo.png (negative lookahead)
content = content.replace(/\(\.\.\/(?!\.\.)([^)"]*?\.(png|gif|jpg|jpeg))(\s+"[^"]*")?\)/gi, '(https://datagrok.ai/help/develop/$1$3)');
// ./foo.png -> https://datagrok.ai/help/develop/packages/foo.png
content = content.replace(/\(\.\/([^)"]*?\.(png|gif|jpg|jpeg))(\s+"[^"]*")?\)/gi, '(https://datagrok.ai/help/develop/packages/$1$3)');

fs.writeFileSync(dest, content);
console.log('Preprocessed JS API readme -> ' + dest);
