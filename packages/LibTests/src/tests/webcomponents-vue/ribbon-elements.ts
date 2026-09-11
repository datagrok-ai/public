import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {createMenuItem, createPanelItem, defaultDebounce, updatePanelItem}
  from '@datagrok-libraries/webcomponents-vue/src/ViewService/ribbon-elements';
import {identityDebounce, panelItem, menuItem} from './ribbon-test-utils';

const deps = (warned: string[] = []) =>
  ({debounce: identityDebounce, warn: (message: string) => warned.push(message)});

category('WebComponentsVue: Ribbon elements', () => {
  test('renders fa icon element', async () => {
    const state = createPanelItem(panelItem({icon: 'bug'}), 'fa:bug', deps());
    expect(state.el.classList.contains('grok-icon'));
    expect(state.el.classList.contains('fal'));
    expect(state.el.classList.contains('fa-bug'));
  });

  test('renders image icon element with size', async () => {
    const item = panelItem({icon: {path: '/icons/i.svg', width: 24, height: 24}});
    const state = createPanelItem(item, 'img:/icons/i.svg', deps());
    expect(state.el.classList.contains('image-icon'));
    expect(state.el.style.backgroundImage.includes('/icons/i.svg'));
    expectDeepEqual(state.el.style.width, '24px');
    expectDeepEqual(state.el.style.height, '24px');
  });

  test('toggles active background', async () => {
    const state = createPanelItem(panelItem(), 'fa:save', deps());
    updatePanelItem(state, panelItem({active: true}));
    expectDeepEqual(state.el.style.backgroundColor, 'var(--grey-1)');
    updatePanelItem(state, panelItem());
    expectDeepEqual(state.el.style.backgroundColor, '');
  });

  test('applies and clears disabled style', async () => {
    const state = createPanelItem(panelItem(), 'fa:save', deps());
    updatePanelItem(state, panelItem({disabled: true}));
    expectDeepEqual(state.el.style.opacity, '0.4');
    expectDeepEqual(state.el.style.filter, 'grayscale(1)');
    updatePanelItem(state, panelItem());
    expectDeepEqual(state.el.style.opacity, '');
    expectDeepEqual(state.el.style.filter, '');
  });

  test('skips disabled style when style is none', async () => {
    const state = createPanelItem(panelItem(), 'fa:save', deps());
    updatePanelItem(state, panelItem({disabled: true, disabledStyle: 'none'}));
    expectDeepEqual(state.el.style.opacity, '');
  });

  test('gates clicks and warns in popup mode', async () => {
    let ran = 0;
    const warned: string[] = [];
    const state = createPanelItem(panelItem({onClick: () => ran++}), 'fa:save', deps(warned));
    state.el.dispatchEvent(new MouseEvent('click'));
    expectDeepEqual(ran, 1);
    updatePanelItem(state, panelItem({
      onClick: () => ran++,
      disabled: true, disabledStyle: 'none', disabledReason: 'busy', disabledReasonMode: 'popup',
    }));
    state.el.dispatchEvent(new MouseEvent('click'));
    expectDeepEqual(ran, 1);
    expectDeepEqual(warned, ['busy']);
  });

  test('debounces disabled style with real timers', async () => {
    const warned: string[] = [];
    const state = createPanelItem(
      panelItem({disabledStyleDebounce: 20}), 'fa:save',
      {debounce: defaultDebounce, warn: (message: string) => warned.push(message)});
    updatePanelItem(state, panelItem({disabled: true, disabledStyleDebounce: 20}));
    expectDeepEqual(state.el.style.opacity, '');
    await new Promise((resolve) => setTimeout(resolve, 60));
    expectDeepEqual(state.el.style.opacity, '0.4');
  });

  test('renders menu item with optional icon', async () => {
    const plain = createMenuItem(menuItem({text: 'Report'}), 'Report|');
    expectDeepEqual(plain.el.textContent, 'Report');
    expectDeepEqual(plain.el.querySelector('i'), null);
    const withIcon = createMenuItem(menuItem({text: 'Run', icon: 'play'}), 'Run|play');
    expect(withIcon.el.textContent?.includes('Run'));
    expect(withIcon.el.querySelector('i')?.classList.contains('fa-play'));
  });
});
