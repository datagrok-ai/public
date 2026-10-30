/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sar/peptide-generation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {peptidesInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, selectIn, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly, everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Peptide Generation viewer designs new peptides and opens them as a table", () => {
  const session = feature(test, "features/sar/peptide-generation.feature", import.meta.url);
  test("The Peptide Generation viewer designs new peptides and opens them as a table", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And simple mode is off", () => simpleModeOff(page));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(24, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await run.scenario("The Add Viewer dialog adds the viewer, which generates 50 peptides for the highest activity", async () => {
      await session.step(27, "Then \"Add viewer\" icon in toolbar should be visible", () => shouldBe(page, el("\"Add viewer\" icon in toolbar"), "visible"));
      await session.step(28, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
      await session.step(29, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
      await session.step(30, "When user types \"peptide\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "peptide", el("viewer gallery search in \"Add Viewer\" dialog")));
      await session.step(31, "And user clicks on first \"Peptide Generation\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Peptide Generation\" card in \"Add Viewer\" dialog")));
      await session.step(32, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
      await session.step(33, "And the open tableview should have 1 Peptide Generation viewer", () => viewerCount(page, 1, "Peptide Generation"));
      await session.step(34, "And Peptide Generation viewer should contain text \"50 peptides generated for highest activity\"", () => shouldContainText(page, el("Peptide Generation viewer"), "50 peptides generated for highest activity"));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Activity, Activity target and Peptide count in the properties change what is generated", async () => {
      await session.step(38, "When user picks \"Properties...\" from the viewer menu of Peptide Generation viewer", () => pickFromViewerMenu(page, "Properties...", el("Peptide Generation viewer")));
      await session.step(39, "And user selects \"IC50\" in \"Activity\" property in context panel", () => selectIn(page, "IC50", el("\"Activity\" property in context panel")));
      await session.step(40, "And user expands \"Misc\" category in context panel", () => expand(page, el("\"Misc\" category in context panel")));
      await session.step(41, "And user selects \"Low\" in \"Activity Target\" property in context panel", () => selectIn(page, "Low", el("\"Activity Target\" property in context panel")));
      await session.step(42, "Then Peptide Generation viewer should contain text \"50 peptides generated for lowest activity\"", () => shouldContainText(page, el("Peptide Generation viewer"), "50 peptides generated for lowest activity"));
      await session.step(43, "When user enters \"10\" into \"Peptide Count\" property in context panel", () => enterInto(page, "10", el("\"Peptide Count\" property in context panel")));
      await session.step(44, "Then Peptide Generation viewer should contain text \"10 peptides generated for lowest activity\"", () => shouldContainText(page, el("Peptide Generation viewer"), "10 peptides generated for lowest activity"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The plus icon opens the generated peptides as a new table", async () => {
      await session.step(48, "When user clicks on \"Open the generated peptides as a new table\" icon in Peptide Generation viewer", () => clickOn(page, el("\"Open the generated peptides as a new table\" icon in Peptide Generation viewer")));
      await session.step(49, "Then the \"Generated Peptides\" table view should open with 10 rows", () => tableViewOpened(page, "Generated Peptides", 10));
      await session.step(50, "And the table should have the columns \"Peptide, Predicted activity, Confidence, Min support, Mean p-value, Significant positions, Novel, Basis, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17\"", () => columnsExactly(page, "Peptide, Predicted activity, Confidence, Min support, Mean p-value, Significant positions, Novel, Basis, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17"));
      await session.step(51, "And every value of \"Novel\" column should match \"^true$\"", () => everyValueMatches(page, "Novel", "^true$"));
      await session.step(52, "And table \"peptides\" should have 647 rows", () => tableRows(page, "peptides", 647));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
