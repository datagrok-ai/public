import {buttons, clickOk, css, expect, field, flush, input, mount, need, sources, type} from '../../lib/harness.mjs';

async function build(mod) {
  const calls = {created: [], cancelled: 0};
  const root = mount(need(mod, 'templateEditor')((v) => void calls.created.push(v), () => void calls.cancelled++));
  await flush();
  const cancel = buttons(root, /^cancel$/)[0];
  const create = buttons(root, /^create$/)[0];
  expect(cancel !== undefined && create !== undefined, 'no CANCEL and CREATE buttons');
  return {root, calls, cancel, create};
}

const INTERNAL = /\.(u2|d4|ui|grok)-|\[data-u2/;
const INLINE_STYLE = /\.style\.\w+\s*=(?!=)|\.style\.(setProperty|cssText)\b|setAttribute\(\s*['"]style['"]/;
const editorOf = (label, root) => field(input(label, root));
const widest = (el) => Math.min(...['width', 'max-width'].map((p) => css.px(css.computed(el, p)) ?? Infinity));

export const checks = {
  'fields': async (mod) => {
    const {root} = await build(mod);
    for (const label of ['Name', 'Description']) {
      const editor = editorOf(label, root);
      const box = editor.closest('.u2-input-box') ?? editor.parentElement;
      expect(css.grow(box) >= 1 && (css.grow(editor) >= 1 || css.computed(editor, 'width') === '100%'),
        `the ${label} editor does not take the width of its row`);
    }
    const key = editorOf('Key', root);
    expect(widest(key) <= 200, `the Key editor is not kept short (${css.computed(key, 'width') ?? 'no width'})`);
  },
  'width-cap': async (mod) => {
    const {root, create} = await build(mod);
    const name = input('Name', root).root;
    const capped = [root, ...root.querySelectorAll('*')].filter((el) => el.contains(name) && el.contains(create))
      .some((el) => widest(el) <= 640);
    expect(capped, 'nothing keeps the editor within 640px');
  },
  'buttons-row': async (mod) => {
    const {cancel, create} = await build(mod);
    const row = cancel.parentElement;
    expect(row.contains(create), 'CANCEL and CREATE are not in one row');
    expect(row.children.indexOf(cancel) < row.children.indexOf(create), 'CREATE comes before CANCEL');
    expect(['flex', 'inline-flex'].includes(css.computed(row, 'display')) &&
      !String(css.computed(row, 'flex-direction')).startsWith('column'), 'the buttons are not laid out in a row');
    expect(['flex-end', 'end', 'right'].includes(css.computed(row, 'justify-content')),
      'the buttons are not right-aligned');
    const gap = css.px(css.computed(row, 'column-gap'));
    expect(gap === 8, `the buttons are ${gap ?? 'not'}px apart, not 8px`);
  },
  'actions': async (mod) => {
    const {root, calls, cancel, create} = await build(mod);
    await type(input('Name', root), 'Kinase panel');
    await type(input('Key', root), 'KIN');
    await type(input('Description', root), 'Kinase hits');
    await clickOk(create, 'CREATE');
    const wanted = [{name: 'Kinase panel', key: 'KIN', description: 'Kinase hits'}];
    expect(JSON.stringify(calls.created) === JSON.stringify(wanted), `onCreate got ${JSON.stringify(calls.created)}`);
    await clickOk(cancel, 'CANCEL');
    expect(calls.cancelled === 1, `onCancel ran ${calls.cancelled} times`);
  },
  'package-css-clean': async (mod) => {
    await build(mod);
    for (const rule of css.parse(css.packageCss())) {
      for (const part of rule.parts)
        expect(!INTERNAL.test(part.selector), `package CSS styles u2's insides: ${part.selector}`);
      const forced = rule.decls.find((d) => d.important);
      expect(forced === undefined, `package CSS forces ${forced?.prop} with !important`);
    }
  },
  'no-inline-styles': async (mod, ctx) => {
    const hit = INLINE_STYLE.exec(sources(ctx.file));
    expect(hit === null, `the solution writes inline styles (${hit?.[0]})`);
  },
};
