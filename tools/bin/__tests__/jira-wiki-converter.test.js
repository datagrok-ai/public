"use strict";

var _vitest = require("vitest");
var _report = require("../commands/report");
(0, _vitest.describe)('markdownToJiraWiki — basic rules', () => {
  (0, _vitest.it)('converts H1', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('# Title')).toBe('h1. Title');
  });
  (0, _vitest.it)('converts H2 / H3 / H6', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('## Sub')).toBe('h2. Sub');
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('### Sub-sub')).toBe('h3. Sub-sub');
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('###### Deep')).toBe('h6. Deep');
  });
  (0, _vitest.it)('converts bold', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('hello **world** foo')).toBe('hello *world* foo');
  });
  (0, _vitest.it)('converts italic with single asterisks', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('hello *world* foo')).toBe('hello _world_ foo');
  });
  (0, _vitest.it)('converts strikethrough', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('~~gone~~')).toBe('-gone-');
  });
  (0, _vitest.it)('converts links', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('see [docs](https://x.example/y)')).toBe('see [docs|https://x.example/y]');
  });
  (0, _vitest.it)('converts blockquote', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('> quoted line')).toBe('bq. quoted line');
  });
  (0, _vitest.it)('converts unordered list', () => {
    const md = '- one\n- two\n- three';
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)(md)).toBe('* one\n* two\n* three');
  });
  (0, _vitest.it)('converts one level of nested unordered list', () => {
    const md = '- top\n  - sub\n- back';
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)(md)).toBe('* top\n** sub\n* back');
  });
  (0, _vitest.it)('converts ordered list', () => {
    const md = '1. one\n2. two\n3. three';
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)(md)).toBe('# one\n# two\n# three');
  });
  (0, _vitest.it)('converts inline code', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('use `foo()` here')).toBe('use {{foo()}} here');
  });
  (0, _vitest.it)('converts plain code fence to noformat', () => {
    const md = '```\nraw code\nmore\n```';
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)(md)).toBe('{noformat}\nraw code\nmore\n{noformat}');
  });
  (0, _vitest.it)('converts code fence with language tag to {code:lang}', () => {
    const md = '```js\nconst x = 1;\n```';
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)(md)).toBe('{code:js}\nconst x = 1;\n{code}');
  });
  (0, _vitest.it)('converts HTML entities &nbsp; / &amp; / &lt; / &gt;', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('a&nbsp;b')).toBe('a b');
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('a&amp;b')).toBe('a&b');
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('a&lt;b&gt;c')).toBe('a<b>c');
  });
});
(0, _vitest.describe)('markdownToJiraWiki — content-protection inside code fences', () => {
  (0, _vitest.it)('does not transform `{{plates}}` inside a fenced block (run-#4 false-positive case)', () => {
    const md = ['Before', '```', 'context: {{plates}} should stay literal', '# not a heading', '- not a list', '```', 'After'].join('\n');
    const out = (0, _report.markdownToJiraWiki)(md);
    (0, _vitest.expect)(out).toContain('{noformat}\ncontext: {{plates}} should stay literal');
    (0, _vitest.expect)(out).toContain('# not a heading');
    (0, _vitest.expect)(out).toContain('- not a list');
    (0, _vitest.expect)(out).not.toContain('h1.');
  });
  (0, _vitest.it)('does not transform headings/links/lists inside an inline code span', () => {
    const md = 'use `# not a heading` and `[ref](u)` here';
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)(md)).toBe('use {{# not a heading}} and {{[ref](u)}} here');
  });
});
(0, _vitest.describe)('markdownToJiraWiki — edge cases', () => {
  (0, _vitest.it)('handles bold containing italic: **bold *inside* bold**', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('**bold *inside* bold**')).toBe('*bold _inside_ bold*');
  });
  (0, _vitest.it)('preserves bold and italic on the same line', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('**a** and *b*')).toBe('*a* and _b_');
  });
  (0, _vitest.it)('handles a heading whose text contains backticks', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('## use `foo()` for x')).toBe('h2. use {{foo()}} for x');
  });
  (0, _vitest.it)('passes plain text through unchanged', () => {
    const plain = 'Just a normal sentence with no markdown.';
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)(plain)).toBe(plain);
  });
  (0, _vitest.it)('returns empty string unchanged', () => {
    (0, _vitest.expect)((0, _report.markdownToJiraWiki)('')).toBe('');
  });
  (0, _vitest.it)('handles a multi-block document end-to-end', () => {
    const md = ['# Handoff', '', 'Some **bold** intro and a [link](https://x.example).', '', '## Findings', '', '- first', '- second', '', '> note: this matters', '', '```js', 'const x = 1;', '```', '', 'Done&nbsp;here.'].join('\n');
    const out = (0, _report.markdownToJiraWiki)(md);
    (0, _vitest.expect)(out).toContain('h1. Handoff');
    (0, _vitest.expect)(out).toContain('h2. Findings');
    (0, _vitest.expect)(out).toContain('Some *bold* intro and a [link|https://x.example].');
    (0, _vitest.expect)(out).toContain('* first\n* second');
    (0, _vitest.expect)(out).toContain('bq. note: this matters');
    (0, _vitest.expect)(out).toContain('{code:js}\nconst x = 1;\n{code}');
    (0, _vitest.expect)(out).toContain('Done here.');
    (0, _vitest.expect)(out).not.toContain('&nbsp;');
  });
});