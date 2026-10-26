/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/activity-cliffs-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.activity-cliffs]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {newColumnMatching, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openDataset, switchTableView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Activity Cliffs regressions from closed issues", () => {
  const session = feature(test, "features/analyze/activity-cliffs-regressions.feature", import.meta.url);
  test("Activity Cliffs regressions from closed issues", {tag: ["@journey", "@realizes:chem.cp.activity-cliffs"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(12, "And user opens activity-cliffs dataset", () => openDataset(page, ds("activity-cliffs")));
    await run.scenario("The Activity cliffs panel stays with the view of its analysis (public#1512)", async () => {
      await session.step(15, "When user picks \"Chem > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Activity Cliffs..."));
      await session.step(16, "And user clicks on OK button in \"Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Activity Cliffs\" dialog")));
      await session.step(17, "Then a new column matching \"^sali_\" should have been added", () => newColumnMatching(page, "^sali_"));
      await session.step(18, "And the \"cliffs\" reading of scatter plot viewer should be 2", () => readingIs(page, "cliffs", el("scatter plot viewer"), 2));
      await session.step(19, "When user clicks on \"cliffs\" button", () => clickOn(page, el("\"cliffs\" button")));
      await session.step(20, "Then \"Activity cliffs\" dock panel should be visible", () => shouldBe(page, el("\"Activity cliffs\" dock panel"), "visible"));
      await session.step(21, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(22, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
      await session.step(23, "And \"Activity cliffs\" dock panel should not be visible", () => shouldNotBe(page, el("\"Activity cliffs\" dock panel"), "visible"));
      await session.step(24, "Given user switches to the \"activity_cliffs_test\" table view", () => switchTableView(page, "activity_cliffs_test"));
      await session.step(25, "Then the \"activity_cliffs_test\" view should be current", () => viewIsCurrent(page, "activity_cliffs_test"));
      await session.step(26, "And \"Activity cliffs\" dock panel should be visible", () => shouldBe(page, el("\"Activity cliffs\" dock panel"), "visible"));
      await session.step(27, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
