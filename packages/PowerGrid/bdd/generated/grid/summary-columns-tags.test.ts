/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/grid/summary-columns-tags.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColor, areaColors, noErrors, pickFromAreaContextMenu, readingReads, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid Tags column", () => {
  const session = feature(test, "features/grid/summary-columns-tags.feature", import.meta.url);
  test("Grid Tags column", {tag: ["@journey", "@viewers", "@realizes:viewers.grid", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await run.scenario("A Tags column marks the rows that carry the flag and leaves the rest of the grid alone (GROK-20888)", async () => {
      await session.step(23, "When user picks \"Add > Summary Columns > Sparklines\" from the context menu of the \"cell 2 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Sparklines", "cell 2 of USUBJID", el("grid")));
      await session.step(24, "And user picks \"Add > Summary Columns > Tags\" from the context menu of the \"cell 2 of Sparklines\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Tags", "cell 2 of Sparklines", el("grid")));
      await session.step(25, "Then the \"cell type of Tags\" reading of grid should be \"tags\"", () => readingReads(page, "cell type of Tags", el("grid"), "tags"));
      await session.step(26, "And the \"cell 1 of Tags\" area of grid should contain the color \"#FFFFFF\"", () => areaColor(page, "cell 1 of Tags", el("grid"), "#FFFFFF"));
      await session.step(27, "And the \"cell 3 of Tags\" area of grid should be painted in at least 2 colors", () => areaColors(page, "cell 3 of Tags", el("grid"), 2));
      await session.step(28, "And the \"cell 1 of Sparklines\" area of grid should contain the color \"#FFFFFF\"", () => areaColor(page, "cell 1 of Sparklines", el("grid"), "#FFFFFF"));
      await session.step(29, "And the \"cell 1 of SEVERITY\" area of grid should contain the color \"#FFFFFF\"", () => areaColor(page, "cell 1 of SEVERITY", el("grid"), "#FFFFFF"));
      await session.step(30, "And the \"text of cell 1 of SEVERITY\" reading of grid should be \"High\"", () => readingReads(page, "text of cell 1 of SEVERITY", el("grid"), "High"));
      await session.step(31, "And the \"text of cell 1 of CONTROL\" reading of grid should be \"false\"", () => readingReads(page, "text of cell 1 of CONTROL", el("grid"), "false"));
      await session.step(32, "And the \"text of cell 3 of CONTROL\" reading of grid should be \"true\"", () => readingReads(page, "text of cell 3 of CONTROL", el("grid"), "true"));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true});
    run.finish();
  });
});
