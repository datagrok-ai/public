import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// DG.ClickUtils: getFullPath, getElementLoggingName, getClickElementDescription,
// sanitizeCssAttrValue. Pure DOM/string logging helpers.
category('AI: Utils: ClickUtils', () => {
  function withElement(body: (el: HTMLElement) => void): void {
    const parent = document.createElement('div');
    parent.id = 'ai-click-parent';
    const el = document.createElement('button');
    el.setAttribute('name', 'ai-click-button');
    el.className = 'ai-click-target';
    el.textContent = 'Click me';
    parent.appendChild(el);
    document.body.appendChild(parent);
    try {
      body(el);
    } finally {
      parent.remove();
    }
  }

  test('getFullPath returns a non-empty path for an attached element', async () => {
    withElement((el) => {
      const path = DG.ClickUtils.getFullPath(el);
      expect(typeof path, 'string');
      expect(path.length > 0, true);
    });
  });

  test('getElementLoggingName returns a string', async () => {
    withElement((el) => {
      expect(typeof DG.ClickUtils.getElementLoggingName(el), 'string');
    });
  });

  test('getClickElementDescription returns a string', async () => {
    withElement((el) => {
      expect(typeof DG.ClickUtils.getClickElementDescription(el), 'string');
    });
  });

  test('sanitizeCssAttrValue preserves a plain value and stays a string', async () => {
    const plain = DG.ClickUtils.sanitizeCssAttrValue('simple-value');
    expect(typeof plain, 'string');
    expect(plain.indexOf('simple-value') >= 0, true);
  });

  test('sanitizeCssAttrValue handles special characters without throwing', async () => {
    for (const v of ['a"b', 'a\'b', 'a\nb', 'a\\b', ''])
      expect(typeof DG.ClickUtils.sanitizeCssAttrValue(v), 'string');
  });
}, {owner: 'agolovko@datagrok.ai'});
