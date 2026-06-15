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

  test('getFullPath reconstructs the path from the element name attribute', async () => {
    withElement((el) => {
      // The button carries name=ai-click-button; the parent div has only an id, which
      // contributes no logging name, so the path resolves to the button's name.
      const path = DG.ClickUtils.getFullPath(el);
      expect(path, 'ai-click-button');
    });
  });

  test('getElementLoggingName returns the element name attribute', async () => {
    withElement((el) => {
      expect(DG.ClickUtils.getElementLoggingName(el), 'ai-click-button');
    });
  });

  test('getClickElementDescription returns the element name attribute', async () => {
    withElement((el) => {
      expect(DG.ClickUtils.getClickElementDescription(el), 'ai-click-button');
    });
  });

  test('sanitizeCssAttrValue passes a plain value through unchanged', async () => {
    expect(DG.ClickUtils.sanitizeCssAttrValue('simple-value'), 'simple-value');
  });

  test('sanitizeCssAttrValue strips the known role prefixes', async () => {
    // button- prefix is dropped and hyphens become spaces.
    expect(DG.ClickUtils.sanitizeCssAttrValue('button-Save'), 'Save');
    // dialog- prefix becomes a "Dialog: " label with hyphens turned into spaces.
    expect(DG.ClickUtils.sanitizeCssAttrValue('dialog-edit-user'), 'Dialog: edit user');
    // input-host- prefix becomes an "Input: " label.
    expect(DG.ClickUtils.sanitizeCssAttrValue('input-host-age'), 'Input: age');
  });
}, {owner: 'agolovko@datagrok.ai'});
