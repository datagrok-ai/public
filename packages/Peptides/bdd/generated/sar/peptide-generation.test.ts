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
import {clickOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {addViewer, noErrors, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Peptide Generation viewer designs new peptides and opens them as a table", () => {
  const session = feature(test, "features/sar/peptide-generation.feature", import.meta.url);
  test("The viewer generates peptides for the target and count set, and its plus icon opens them as a table", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(14, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(15, "And user adds Peptide Generation viewer", () => addViewer(page, "Peptide Generation"));
    await session.step(16, "Then Peptide Generation viewer should contain text \"50 peptides generated for highest activity\"", () => shouldContainText(page, el("Peptide Generation viewer"), "50 peptides generated for highest activity"));
    await session.step(17, "When user sets \"activityColumnName\" property of Peptide Generation viewer to \"IC50\"", () => setProperty(page, "activityColumnName", el("Peptide Generation viewer"), "IC50"));
    await session.step(18, "And user sets \"activityTarget\" property of Peptide Generation viewer to \"Low\"", () => setProperty(page, "activityTarget", el("Peptide Generation viewer"), "Low"));
    await session.step(19, "And user sets \"peptideCount\" property of Peptide Generation viewer to \"10\"", () => setProperty(page, "peptideCount", el("Peptide Generation viewer"), "10"));
    await session.step(20, "Then Peptide Generation viewer should contain text \"10 peptides generated for lowest activity\"", () => shouldContainText(page, el("Peptide Generation viewer"), "10 peptides generated for lowest activity"));
    await session.step(21, "When user clicks on \"Open the generated peptides as a new table\" icon in Peptide Generation viewer", () => clickOn(page, el("\"Open the generated peptides as a new table\" icon in Peptide Generation viewer")));
    await session.step(22, "Then the \"Generated Peptides\" table view should open with 10 rows", () => tableViewOpened(page, "Generated Peptides", 10));
    await session.step(23, "And no errors should have been logged", () => noErrors(page));
  });
});
