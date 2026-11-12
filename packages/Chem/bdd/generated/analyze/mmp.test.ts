/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/mmp.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.mmp-analysis, public#4435]
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
import {mmpReady} from '../../bindings/scaffold-tree.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterBetween, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clickArea, noErrors, readingAtLeast, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Matched Molecular Pairs over molecules and their activities", () => {
  const session = feature(test, "features/analyze/mmp.feature", import.meta.url);
  test("Matched Molecular Pairs over molecules and their activities", {tag: ["@journey", "@realizes:chem.cp.mmp-analysis", "@GROK-18517", "@realizes:public#4435"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens sar-small dataset", () => openDataset(page, ds("sar-small")));
    await run.scenario("The dialog opens on the molecules column and the activity", async () => {
      await session.step(18, "When user picks \"Chem > Analyze > Matched Molecular Pairs...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Matched Molecular Pairs..."));
      await session.step(19, "Then \"Matched Molecular Pairs\" dialog should be visible", () => shouldBe(page, el("\"Matched Molecular Pairs\" dialog"), "visible"));
      await session.step(20, "And Column input in \"Matched Molecular Pairs\" dialog should contain text \"smiles\"", () => shouldContainText(page, el("Column input in \"Matched Molecular Pairs\" dialog"), "smiles"));
      await session.step(21, "And Activities input in \"Matched Molecular Pairs\" dialog should contain text \"Activities(0)\"", () => shouldContainText(page, el("Activities input in \"Matched Molecular Pairs\" dialog"), "Activities(0)"));
      await session.step(22, "When user clicks on editor of Activities input in \"Matched Molecular Pairs\" dialog", () => clickOn(page, el("editor of Activities input in \"Matched Molecular Pairs\" dialog")));
      await session.step(23, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(24, "When user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(25, "And user clicks on the \"cell 6 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 6 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(26, "Then the \"text of cell 6 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"LD(50)\"", () => readingReads(page, "text of cell 6 of __name", el("grid viewer in \"Select columns...\" dialog"), "LD(50)"));
      await session.step(27, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(28, "Then Activities input in \"Matched Molecular Pairs\" dialog should contain text \"Activities(1)\"", () => shouldContainText(page, el("Activities input in \"Matched Molecular Pairs\" dialog"), "Activities(1)"));
      await session.step(29, "And Scaling input in \"Matched Molecular Pairs\" dialog should be visible", () => shouldBe(page, el("Scaling input in \"Matched Molecular Pairs\" dialog"), "visible"));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The run builds the viewer with its substitutions and pairs", async () => {
      await session.step(33, "When user clicks on OK button in \"Matched Molecular Pairs\" dialog", () => clickOn(page, el("OK button in \"Matched Molecular Pairs\" dialog")));
      await session.step(34, "Then Matched Molecular Pairs Analysis viewer should be visible", () => shouldBe(page, el("Matched Molecular Pairs Analysis viewer"), "visible"));
      await session.step(35, "And Matched Molecular Pairs Analysis viewer should have finished its analysis", () => mmpReady(page, el("Matched Molecular Pairs Analysis viewer")));
      await session.step(36, "And the \"activities\" reading of Matched Molecular Pairs Analysis viewer should be \"LD(50)\"", () => readingReads(page, "activities", el("Matched Molecular Pairs Analysis viewer"), "LD(50)"));
      await session.step(37, "And the \"molecules column\" reading of Matched Molecular Pairs Analysis viewer should be \"smiles\"", () => readingReads(page, "molecules column", el("Matched Molecular Pairs Analysis viewer"), "smiles"));
      await session.step(38, "And the \"substitutions\" reading of Matched Molecular Pairs Analysis viewer should be at least 1", () => readingAtLeast(page, "substitutions", el("Matched Molecular Pairs Analysis viewer"), 1));
      await session.step(39, "And the \"pairs\" reading of Matched Molecular Pairs Analysis viewer should be at least 1", () => readingAtLeast(page, "pairs", el("Matched Molecular Pairs Analysis viewer"), 1));
      await session.step(40, "And the table should have 200 rows", () => rowCount(page, 200));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The tabs of the viewer are its four analyses", async () => {
      await session.step(44, "Then Substitutions tab in Matched Molecular Pairs Analysis viewer should be visible", () => shouldBe(page, el("Substitutions tab in Matched Molecular Pairs Analysis viewer"), "visible"));
      await session.step(45, "And Fragments tab in Matched Molecular Pairs Analysis viewer should be visible", () => shouldBe(page, el("Fragments tab in Matched Molecular Pairs Analysis viewer"), "visible"));
      await session.step(46, "And Cliffs tab in Matched Molecular Pairs Analysis viewer should be visible", () => shouldBe(page, el("Cliffs tab in Matched Molecular Pairs Analysis viewer"), "visible"));
      await session.step(47, "And Generation tab in Matched Molecular Pairs Analysis viewer should be visible", () => shouldBe(page, el("Generation tab in Matched Molecular Pairs Analysis viewer"), "visible"));
      await session.step(48, "And the \"tab\" reading of Matched Molecular Pairs Analysis viewer should be \"Substitutions\"", () => readingReads(page, "tab", el("Matched Molecular Pairs Analysis viewer"), "Substitutions"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Generation tab fills its grid", async () => {
      await session.step(52, "When user clicks on Generation tab in Matched Molecular Pairs Analysis viewer", () => clickOn(page, el("Generation tab in Matched Molecular Pairs Analysis viewer")));
      await session.step(53, "Then the \"tab\" reading of Matched Molecular Pairs Analysis viewer should be \"Generation\"", () => readingReads(page, "tab", el("Matched Molecular Pairs Analysis viewer"), "Generation"));
      await session.step(54, "And the \"generated\" reading of Matched Molecular Pairs Analysis viewer should be at least 1", () => readingAtLeast(page, "generated", el("Matched Molecular Pairs Analysis viewer"), 1));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A run with two activities names both and gives each its mean difference", async () => {
      await session.step(59, "When user clicks on close icon of Matched Molecular Pairs Analysis viewer", () => clickOn(page, el("close icon of Matched Molecular Pairs Analysis viewer")));
      await session.step(60, "Then Matched Molecular Pairs Analysis viewer should be absent", () => shouldBe(page, el("Matched Molecular Pairs Analysis viewer"), "absent"));
      await session.step(61, "When user picks \"Chem > Analyze > Matched Molecular Pairs...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Matched Molecular Pairs..."));
      await session.step(62, "Then \"Matched Molecular Pairs\" dialog should be visible", () => shouldBe(page, el("\"Matched Molecular Pairs\" dialog"), "visible"));
      await session.step(63, "When user clicks on editor of Activities input in \"Matched Molecular Pairs\" dialog", () => clickOn(page, el("editor of Activities input in \"Matched Molecular Pairs\" dialog")));
      await session.step(64, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(65, "When user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(66, "And user clicks on the \"cell 6 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 6 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(67, "And user clicks on the \"cell 7 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 7 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(68, "Then the \"text of cell 6 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"LD(50)\"", () => readingReads(page, "text of cell 6 of __name", el("grid viewer in \"Select columns...\" dialog"), "LD(50)"));
      await session.step(69, "And the \"text of cell 7 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"In-vivo Activity\"", () => readingReads(page, "text of cell 7 of __name", el("grid viewer in \"Select columns...\" dialog"), "In-vivo Activity"));
      await session.step(70, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(71, "Then Activities input in \"Matched Molecular Pairs\" dialog should contain text \"Activities(2)\"", () => shouldContainText(page, el("Activities input in \"Matched Molecular Pairs\" dialog"), "Activities(2)"));
      await session.step(72, "When user clicks on OK button in \"Matched Molecular Pairs\" dialog", () => clickOn(page, el("OK button in \"Matched Molecular Pairs\" dialog")));
      await session.step(73, "Then Matched Molecular Pairs Analysis viewer should be visible", () => shouldBe(page, el("Matched Molecular Pairs Analysis viewer"), "visible"));
      await session.step(74, "And Matched Molecular Pairs Analysis viewer should have finished its analysis", () => mmpReady(page, el("Matched Molecular Pairs Analysis viewer")));
      await session.step(75, "And the \"activities\" reading of Matched Molecular Pairs Analysis viewer should be \"LD(50), In-vivo Activity\"", () => readingReads(page, "activities", el("Matched Molecular Pairs Analysis viewer"), "LD(50), In-vivo Activity"));
      await session.step(76, "And the \"fragment columns\" reading of Matched Molecular Pairs Analysis viewer should contain \"Δ LD(50)\"", () => readingContains(page, "fragment columns", el("Matched Molecular Pairs Analysis viewer"), "Δ LD(50)"));
      await session.step(77, "And the \"fragment columns\" reading of Matched Molecular Pairs Analysis viewer should contain \"Δ In-vivo Activity\"", () => readingContains(page, "fragment columns", el("Matched Molecular Pairs Analysis viewer"), "Δ In-vivo Activity"));
      await session.step(78, "And the \"pairs\" reading of Matched Molecular Pairs Analysis viewer should be at least 1", () => readingAtLeast(page, "pairs", el("Matched Molecular Pairs Analysis viewer"), 1));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A run on filtered data opens one view of the rows that pass", async () => {
      await session.step(83, "When user clicks on close icon of Matched Molecular Pairs Analysis viewer", () => clickOn(page, el("close icon of Matched Molecular Pairs Analysis viewer")));
      await session.step(84, "And user filters rows where \"Mol Wt.\" is between 300 and 1000", () => filterBetween(page, "Mol Wt.", 300, 1000));
      await session.step(85, "And user picks \"Chem > Analyze > Matched Molecular Pairs...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Matched Molecular Pairs..."));
      await session.step(86, "And user clicks on editor of Activities input in \"Matched Molecular Pairs\" dialog", () => clickOn(page, el("editor of Activities input in \"Matched Molecular Pairs\" dialog")));
      await session.step(87, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(88, "And user clicks on the \"cell 6 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 6 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(89, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(90, "And user clicks on OK button in \"Matched Molecular Pairs\" dialog", () => clickOn(page, el("OK button in \"Matched Molecular Pairs\" dialog")));
      await session.step(91, "Then Matched Molecular Pairs Analysis viewer should have finished its analysis", () => mmpReady(page, el("Matched Molecular Pairs Analysis viewer")));
      await session.step(92, "And the open table views should be exactly \"sar_small, sar_small (Filtered)\"", () => openTableViewsExactly(page, "sar_small, sar_small (Filtered)"));
    });
    run.finish();
  });
});
