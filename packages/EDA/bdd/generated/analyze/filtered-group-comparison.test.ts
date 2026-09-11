/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/filtered-group-comparison.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [ml.menu.analyze.group-comparison.control-comparisons]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, filterTo} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Group comparison of a filtered table", () => {
  const session = feature(test, "features/analyze/filtered-group-comparison.feature", import.meta.url);
  test("Group comparison of a filtered table", {tag: ["@journey", "@eda", "@realizes:ml.menu.analyze.group-comparison.control-comparisons", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await run.scenario("Control comparisons count only the women when the table is filtered to them", async () => {
      await session.step(17, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(18, "Then 3243 rows should pass the filter", () => filterPasses(page, 3243));
      await session.step(19, "When user picks \"ML > Analyze > Group Comparison > Control Comparisons...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Group Comparison > Control Comparisons..."));
      await session.step(20, "Then 3243 rows should pass the filter", () => filterPasses(page, 3243));
      await session.step(21, "When user clicks on Run button in \"Control comparisons\" dialog", () => clickOn(page, el("Run button in \"Control comparisons\" dialog")));
      await session.step(22, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(23, "And 3243 rows should pass the filter", () => filterPasses(page, 3243));
      await session.step(24, "And second grid viewer should be bound to table \"Control comparisons result\"", () => boundTable(page, el("second grid viewer"), "Control comparisons result"));
      await session.step(25, "And the \"text of cell 1 of Group\" reading of second grid viewer should be \"Black\"", () => readingReads(page, "text of cell 1 of Group", el("second grid viewer"), "Black"));
      await session.step(26, "And the \"text of cell 1 of n\" reading of second grid viewer should be \"104\"", () => readingReads(page, "text of cell 1 of n", el("second grid viewer"), "104"));
      await session.step(27, "And the \"text of cell 2 of n\" reading of second grid viewer should be \"2823\"", () => readingReads(page, "text of cell 2 of n", el("second grid viewer"), "2823"));
      await session.step(28, "And the \"text of cell 3 of n\" reading of second grid viewer should be \"279\"", () => readingReads(page, "text of cell 3 of n", el("second grid viewer"), "279"));
    }, {knownFailure: true});
    run.finish();
  });
});
