/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/data/pzfx-files.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/curves.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {doubleClickOn, isExpanded} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelOpen, packageInstalled, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A GraphPad Prism file opens from Browse", () => {
  const session = feature(test, "features/data/pzfx-files.feature", import.meta.url);
  test("A double-click previews the file in a view of its name", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And simple mode is off", () => simpleModeOff(page));
    await session.step(10, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(11, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(12, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(13, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(14, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(15, "And Files---App-Data---Curves tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Curves tree node inside browse tree")));
    await session.step(18, "When user double-clicks on Files---App-Data---Curves---dose-response-4pl.pzfx tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Curves---dose-response-4pl.pzfx tree node inside browse tree")));
    await session.step(19, "Then the \"dose-response-4pl.pzfx\" view should be current", () => viewIsCurrent(page, "dose-response-4pl.pzfx"));
    await session.step(20, "And no errors should have been logged", () => noErrors(page));
  });
  test("Open from the context menu opens the file", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And simple mode is off", () => simpleModeOff(page));
    await session.step(10, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(11, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(12, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(13, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(14, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(15, "And Files---App-Data---Curves tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Curves tree node inside browse tree")));
    await session.step(23, "When user picks \"Open\" from the context menu of Files---App-Data---Curves---dose-response-4pl.pzfx tree node inside browse tree", () => pickFromContextMenu(page, "Open", el("Files---App-Data---Curves---dose-response-4pl.pzfx tree node inside browse tree")));
    await session.step(24, "Then the \"dose-response-4pl.pzfx\" view should be current", () => viewIsCurrent(page, "dose-response-4pl.pzfx"));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
  });
});
