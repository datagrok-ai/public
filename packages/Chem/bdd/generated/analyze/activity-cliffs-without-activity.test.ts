/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/activity-cliffs-without-activity.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [GROK-20915]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clickAtOnce, okDisabledThroughout, quietWindow, watchNextOk} from '../../bindings/dialogs.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {noNewColumn, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Activity Cliffs on a table with no numeric column", () => {
  const session = feature(test, "features/analyze/activity-cliffs-without-activity.feature", import.meta.url);
  test("OK clicked at once runs nothing and the dialog stays open", {tag: ["@regression", "@realizes:GROK-20915"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And user opens smiles-only dataset", () => openDataset(page, ds("smiles-only")));
    await session.step(21, "Given user watches the OK button of the next dialog", () => watchNextOk(page));
    await session.step(22, "When user picks \"Chem > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Activity Cliffs..."));
    await session.step(23, "And user clicks on OK button in \"Activity Cliffs\" dialog at once", () => clickAtOnce(page, el("OK button in \"Activity Cliffs\" dialog")));
    await session.step(24, "Then the OK button should have been disabled when it appeared and when it was clicked", () => okDisabledThroughout(page));
    await session.step(25, "And no error or warning balloon and no error should appear for 2 seconds", () => quietWindow(page, 2));
    await session.step(26, "And \"Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Activity Cliffs\" dialog"), "visible"));
    await session.step(27, "And OK button in \"Activity Cliffs\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Activity Cliffs\" dialog"), "disabled"));
    await session.step(28, "And no new column should have been added", () => noNewColumn(page));
    await session.step(29, "And the table should have 1000 rows", () => rowCount(page, 1000));
  });
});
