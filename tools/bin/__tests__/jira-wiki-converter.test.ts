import {describe, it, expect} from 'vitest';
import {markdownToJiraWiki} from '../commands/report';

describe('markdownToJiraWiki — basic rules', () => {
  it('converts H1', () => {
    expect(markdownToJiraWiki('# Title')).toBe('h1. Title');
  });

  it('converts H2 / H3 / H6', () => {
    expect(markdownToJiraWiki('## Sub')).toBe('h2. Sub');
    expect(markdownToJiraWiki('### Sub-sub')).toBe('h3. Sub-sub');
    expect(markdownToJiraWiki('###### Deep')).toBe('h6. Deep');
  });

  it('converts bold', () => {
    expect(markdownToJiraWiki('hello **world** foo')).toBe('hello *world* foo');
  });

  it('converts italic with single asterisks', () => {
    expect(markdownToJiraWiki('hello *world* foo')).toBe('hello _world_ foo');
  });

  it('converts strikethrough', () => {
    expect(markdownToJiraWiki('~~gone~~')).toBe('-gone-');
  });

  it('converts links', () => {
    expect(markdownToJiraWiki('see [docs](https://x.example/y)'))
      .toBe('see [docs|https://x.example/y]');
  });

  it('converts blockquote', () => {
    expect(markdownToJiraWiki('> quoted line')).toBe('bq. quoted line');
  });

  it('converts unordered list', () => {
    const md = '- one\n- two\n- three';
    expect(markdownToJiraWiki(md)).toBe('* one\n* two\n* three');
  });

  it('converts one level of nested unordered list', () => {
    const md = '- top\n  - sub\n- back';
    expect(markdownToJiraWiki(md)).toBe('* top\n** sub\n* back');
  });

  it('converts ordered list', () => {
    const md = '1. one\n2. two\n3. three';
    expect(markdownToJiraWiki(md)).toBe('# one\n# two\n# three');
  });

  it('converts inline code', () => {
    expect(markdownToJiraWiki('use `foo()` here')).toBe('use {{foo()}} here');
  });

  it('converts plain code fence to noformat', () => {
    const md = '```\nraw code\nmore\n```';
    expect(markdownToJiraWiki(md)).toBe('{noformat}\nraw code\nmore\n{noformat}');
  });

  it('converts code fence with language tag to {code:lang}', () => {
    const md = '```js\nconst x = 1;\n```';
    expect(markdownToJiraWiki(md)).toBe('{code:js}\nconst x = 1;\n{code}');
  });

  it('converts HTML entities &nbsp; / &amp; / &lt; / &gt;', () => {
    expect(markdownToJiraWiki('a&nbsp;b')).toBe('a b');
    expect(markdownToJiraWiki('a&amp;b')).toBe('a&b');
    expect(markdownToJiraWiki('a&lt;b&gt;c')).toBe('a<b>c');
  });
});

describe('markdownToJiraWiki — content-protection inside code fences', () => {
  it('does not transform `{{plates}}` inside a fenced block (run-#4 false-positive case)', () => {
    const md = [
      'Before',
      '```',
      'context: {{plates}} should stay literal',
      '# not a heading',
      '- not a list',
      '```',
      'After',
    ].join('\n');
    const out = markdownToJiraWiki(md);
    expect(out).toContain('{noformat}\ncontext: {{plates}} should stay literal');
    expect(out).toContain('# not a heading');
    expect(out).toContain('- not a list');
    expect(out).not.toContain('h1.');
  });

  it('does not transform headings/links/lists inside an inline code span', () => {
    const md = 'use `# not a heading` and `[ref](u)` here';
    expect(markdownToJiraWiki(md))
      .toBe('use {{# not a heading}} and {{[ref](u)}} here');
  });
});

describe('markdownToJiraWiki — edge cases', () => {
  it('handles bold containing italic: **bold *inside* bold**', () => {
    expect(markdownToJiraWiki('**bold *inside* bold**'))
      .toBe('*bold _inside_ bold*');
  });

  it('preserves bold and italic on the same line', () => {
    expect(markdownToJiraWiki('**a** and *b*')).toBe('*a* and _b_');
  });

  it('handles a heading whose text contains backticks', () => {
    expect(markdownToJiraWiki('## use `foo()` for x'))
      .toBe('h2. use {{foo()}} for x');
  });

  it('passes plain text through unchanged', () => {
    const plain = 'Just a normal sentence with no markdown.';
    expect(markdownToJiraWiki(plain)).toBe(plain);
  });

  it('returns empty string unchanged', () => {
    expect(markdownToJiraWiki('')).toBe('');
  });

  it('handles a multi-block document end-to-end', () => {
    const md = [
      '# Handoff',
      '',
      'Some **bold** intro and a [link](https://x.example).',
      '',
      '## Findings',
      '',
      '- first',
      '- second',
      '',
      '> note: this matters',
      '',
      '```js',
      'const x = 1;',
      '```',
      '',
      'Done&nbsp;here.',
    ].join('\n');
    const out = markdownToJiraWiki(md);
    expect(out).toContain('h1. Handoff');
    expect(out).toContain('h2. Findings');
    expect(out).toContain('Some *bold* intro and a [link|https://x.example].');
    expect(out).toContain('* first\n* second');
    expect(out).toContain('bq. note: this matters');
    expect(out).toContain('{code:js}\nconst x = 1;\n{code}');
    expect(out).toContain('Done here.');
    expect(out).not.toContain('&nbsp;');
  });
});
