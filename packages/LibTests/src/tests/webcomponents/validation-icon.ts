import * as ui from 'datagrok-api/ui';
import {category, test, expect, waitUntil, delay} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import type {
  ValidationIcon, ValidationIconInput,
} from '@datagrok-libraries/webcomponents/src/ValidationIcon/ValidationIcon';
import {awaitWebComponents} from '../../test-utils';
import {ADVICE, consistency, findLink, iconOf, openPopover} from './fixtures';

const iconName = (root: HTMLElement) => iconOf(root)?.className.match(/fa-([a-z-]+)/)?.[1];

async function mountIcon(setup?: (icon: ValidationIcon) => void) {
  await awaitWebComponents();
  const host = ui.div();
  document.body.append(host);
  const icon = document.createElement('dg-validation-icon') as ValidationIcon;
  setup?.(icon);
  host.append(icon);
  const events: [string, any][] = [];
  for (const name of ['consistency-reset', 'action-request', 'show-dataframe-diff'])
    icon.addEventListener(name, (ev: any) => events.push([name, ev.additionalParams ?? ev.detail ?? null]));
  return {host, icon, events, dispose: () => host.remove()};
}

category('WebComponents: ValidationIcon', () => {
  test('the most severe status picks the icon and its color', async () => {
    const {host, icon, dispose} = await mountIcon();
    try {
      const cases: [ValidationIconInput, string, string][] = [
        [{validation: {errors: ['e'], warnings: ['w']}}, 'exclamation-circle', 'var(--red-3)'],
        [{validation: {warnings: ['w'], notifications: ['n']}}, 'exclamation-circle', 'var(--orange-2)'],
        [{validation: {notifications: ['n']}}, 'info-circle', 'var(--blue-1)'],
        [{validation: {errors: ['e']}, consistency: consistency('restricted', true)},
          'handshake-slash', 'var(--red-3)'],
        [{consistency: consistency('info', true)}, 'handshake-slash', 'var(--blue-1)'],
      ];
      for (const [status, name, color] of cases) {
        icon.validationStatus = status;
        await waitUntil(() => iconName(host) === name && iconOf(host)?.style.getPropertyValue('color') === color,
          () => `${JSON.stringify(status)}: expected ${name} ${color}, got ${iconName(host)} ` +
            `${iconOf(host)?.style.getPropertyValue('color')}`, 3000);
      }
    } finally {
      dispose();
    }
  });

  test('shows a link for workflow-set values', async () => {
    const {host, icon, dispose} = await mountIcon();
    try {
      for (const [status, expected] of [
        [consistency('restricted'), 'link'],
        [consistency('disabled'), 'link'],
        [consistency('restricted', true), 'handshake-slash'],
        [consistency('info'), undefined],
      ] as const) {
        icon.validationStatus = {consistency: status};
        await waitUntil(() => iconName(host) === expected,
          `${status.restriction}: expected ${expected}, got ${iconName(host)}`, 3000);
      }
    } finally {
      dispose();
    }
  });

  test('clearing the status hides the icon', async () => {
    const {host, icon, dispose} = await mountIcon();
    try {
      icon.validationStatus = {validation: {errors: ['e']}};
      await waitUntil(() => icon.style.display === 'flex', 'icon not shown', 3000);
      icon.validationStatus = undefined;
      await waitUntil(() => icon.style.display === 'none', 'icon not hidden', 3000);
      expect(iconOf(host), null);
    } finally {
      dispose();
    }
  });

  test('an advice link sends its additional params with the action request', async () => {
    const {host, icon, events, dispose} = await mountIcon();
    try {
      icon.validationStatus = ADVICE;
      const popover = await openPopover(host);
      findLink(popover, 'Set to 10')!.click();
      expectDeepEqual(events, [['action-request', {val: 10}]]);
    } finally {
      dispose();
    }
  });

  test('an inconsistent value offers a reset, and a diff for dataframes only', async () => {
    for (const isDataFrame of [false, true]) {
      const {host, icon, events, dispose} = await mountIcon((icon) => icon.isDataFrame = isDataFrame);
      try {
        icon.validationStatus = {consistency: consistency('restricted', true, 'computed')};
        const popover = await openPopover(host);
        expect(!!findLink(popover, 'Show dataframes diff'), isDataFrame);
        findLink(popover, 'Reset to consistent value')!.click();
        findLink(popover, 'Show dataframes diff')?.click();
        expectDeepEqual(events.map(([name]) => name),
          isDataFrame ? ['consistency-reset', 'show-dataframe-diff'] : ['consistency-reset']);
        if (isDataFrame)
          expect(events[1][1], 'computed');
      } finally {
        dispose();
      }
    }
  });

  test('hover mode shows a tooltip instead of a popover', async () => {
    const {host, icon, dispose} = await mountIcon((icon) => icon.useHover = true);
    try {
      icon.validationStatus = {validation: {errors: ['e']}};
      await waitUntil(() => iconName(host) === 'exclamation-circle', 'icon not rendered', 3000);
      expect(host.querySelector('.dg-validation-icon-popover'), null);
    } finally {
      dispose();
    }
  });

  test('a click inside the popover keeps its content', async () => {
    const {host, icon, dispose} = await mountIcon();
    try {
      icon.validationStatus = {validation: {errors: ['A long explanation']}};
      const popover = await openPopover(host);
      const content = popover.firstElementChild as HTMLElement;
      (content.querySelector('div') ?? content).click();
      await delay(50);
      expect(popover.firstElementChild === content, true, 'the popover re-rendered on a click inside it');
    } finally {
      dispose();
    }
  });
});
