/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sar/mutation-cliffs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {peptidesInitialized, sarReady} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColors, areaNoColor, areaPainted, areaRepainted, noBalloons, noErrors, painted, readingIs, readingReads, setProperty, takeSnapshot, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Compute and visualize peptide mutation cliffs", () => {
  const session = feature(test, "features/sar/mutation-cliffs.feature", import.meta.url);
  test("Compute and visualize peptide mutation cliffs", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(18, "And user opens peptides dataset keeping the first 200 rows", () => openDatasetRows(page, ds("peptides"), 200));
    await session.step(19, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(20, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(21, "When user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(22, "And user enters \"93\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "93", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(23, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(24, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(25, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(26, "Then no errors should have been logged", () => noErrors(page));
    await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("Switching from counts to mutation cliffs repaints the populated cell", async () => {
      await session.step(30, "When user clicks on \"Invariant Map\" checkbox in Sequence Variability Map viewer", () => clickOn(page, el("\"Invariant Map\" checkbox in Sequence Variability Map viewer")));
      await session.step(31, "Then the \"mode\" reading of Sequence Variability Map viewer should be \"Invariant Map\"", () => readingReads(page, "mode", el("Sequence Variability Map viewer"), "Invariant Map"));
      await session.step(32, "And the \"cell A at 2\" area of Sequence Variability Map viewer should be painted", () => areaPainted(page, "cell A at 2", el("Sequence Variability Map viewer")));
      await session.step(33, "And the \"cell A at 3\" area of Sequence Variability Map viewer should be painted in at least 1 colors", () => areaColors(page, "cell A at 3", el("Sequence Variability Map viewer"), 1));
      await session.step(34, "When user takes a snapshot of Sequence Variability Map viewer", () => takeSnapshot(page, el("Sequence Variability Map viewer")));
      await session.step(35, "When user clicks on \"Mutation Cliffs\" checkbox in Sequence Variability Map viewer", () => clickOn(page, el("\"Mutation Cliffs\" checkbox in Sequence Variability Map viewer")));
      await session.step(36, "Then the \"mode\" reading of Sequence Variability Map viewer should be \"Mutation Cliffs\"", () => readingReads(page, "mode", el("Sequence Variability Map viewer"), "Mutation Cliffs"));
      await session.step(37, "And the \"cell A at 2\" area of Sequence Variability Map viewer should have repainted", () => areaRepainted(page, "cell A at 2", el("Sequence Variability Map viewer")));
      await session.step(38, "And the \"cell A at 2\" area of Sequence Variability Map viewer should be painted in at least 1 colors", () => areaColors(page, "cell A at 2", el("Sequence Variability Map viewer"), 1));
      await session.step(39, "And the \"cell A at 3\" area of Sequence Variability Map viewer should be painted in no color", () => areaNoColor(page, "cell A at 3", el("Sequence Variability Map viewer")));
      await session.step(40, "Then no errors should have been logged", () => noErrors(page));
      await session.step(41, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The position chart draws the peptides of position-two cliffs and none at position one", async () => {
      await session.step(44, "Given user adds a Sequence Mutation Cliffs viewer with:", () => addViewerWith(page, "Sequence Mutation Cliffs", [["sequenceColumnName","AlignedSequence"],["activityColumnName","IC50"],["position","2"]]), [["sequenceColumnName","AlignedSequence"],["activityColumnName","IC50"],["position","2"]]);
      await session.step(48, "Then Sequence Mutation Cliffs viewer should be added to the open tableview", () => viewerAdded(page, "Sequence Mutation Cliffs"));
      await session.step(49, "And the \"position\" reading of Sequence Mutation Cliffs viewer should be 2", () => readingIs(page, "position", el("Sequence Mutation Cliffs viewer"), 2));
      await session.step(50, "And line chart viewer in Sequence Mutation Cliffs viewer should be painted", () => painted(page, el("line chart viewer in Sequence Mutation Cliffs viewer")));
      await session.step(51, "When user sets \"position\" property of Sequence Mutation Cliffs viewer to \"1\"", () => setProperty(page, "position", el("Sequence Mutation Cliffs viewer"), "1"));
      await session.step(52, "Then the \"position\" reading of Sequence Mutation Cliffs viewer should be 1", () => readingIs(page, "position", el("Sequence Mutation Cliffs viewer"), 1));
      await session.step(53, "And the \"message\" reading of Sequence Mutation Cliffs viewer should be \"No mutation cliffs found for the selected position.\"", () => readingReads(page, "message", el("Sequence Mutation Cliffs viewer"), "No mutation cliffs found for the selected position."));
      await session.step(54, "When user sets \"position\" property of Sequence Mutation Cliffs viewer to \"2\"", () => setProperty(page, "position", el("Sequence Mutation Cliffs viewer"), "2"));
      await session.step(55, "Then the \"position\" reading of Sequence Mutation Cliffs viewer should be 2", () => readingIs(page, "position", el("Sequence Mutation Cliffs viewer"), 2));
      await session.step(56, "And line chart viewer in Sequence Mutation Cliffs viewer should be painted", () => painted(page, el("line chart viewer in Sequence Mutation Cliffs viewer")));
      await session.step(57, "Then no errors should have been logged", () => noErrors(page));
      await session.step(58, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The cluster summary is built on the MCL clusters", async () => {
      await session.step(61, "Then MCL viewer should be added to the open tableview", () => viewerAdded(page, "MCL"));
      await session.step(62, "And Logo Summary Table viewer should be added to the open tableview", () => viewerAdded(page, "Logo Summary Table"));
      await session.step(63, "And the table should have a column \"Cluster (MCL)\"", () => hasColumn(page, "Cluster (MCL)"));
      await session.step(64, "And the \"clusters column\" reading of Logo Summary Table viewer should be \"Cluster (MCL)\"", () => readingReads(page, "clusters column", el("Logo Summary Table viewer"), "Cluster (MCL)"));
      await session.step(65, "And Logo Summary Table viewer should be painted", () => painted(page, el("Logo Summary Table viewer")));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
