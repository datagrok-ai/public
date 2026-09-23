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
import {cliffParticipants} from '../../bindings/cliffs.js';
import {mutationPairs} from '../../bindings/exports.js';
import {peptidesInitialized, sarReady} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, tableColumns} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {closeCurrentView, openDatasetRows, switchTableView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColors, areaNoColor, areaPainted, areaRepainted, noBalloons, noErrors, painted, pickFromContextMenu, readingAtLeast, readingBetween, readingFinite, readingIs, readingReads, setProperty, showsRows, takeSnapshot, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Compute and visualize peptide mutation cliffs", () => {
  const session = feature(test, "features/sar/mutation-cliffs.feature", import.meta.url);
  test("Compute and visualize peptide mutation cliffs", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(13, "And user opens peptides dataset keeping the first 200 rows", () => openDatasetRows(page, ds("peptides"), 200));
    await session.step(14, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(15, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(16, "When user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(17, "And user enters \"93\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "93", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(18, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(19, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(20, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(21, "Then no errors should have been logged", () => noErrors(page));
    await session.step(22, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("The variability map reports computed mutation pairs and populated statistics", async () => {
      await session.step(25, "Then the table should have 200 rows", () => rowCount(page, 200));
      await session.step(26, "And the \"positions\" reading of Sequence Variability Map viewer should be 17", () => readingIs(page, "positions", el("Sequence Variability Map viewer"), 17));
      await session.step(27, "And the \"cliff cells\" reading of Sequence Variability Map viewer should be at least 1", () => readingAtLeast(page, "cliff cells", el("Sequence Variability Map viewer"), 1));
      await session.step(28, "And the \"cliff pairs\" reading of Sequence Variability Map viewer should be 4484", () => readingIs(page, "cliff pairs", el("Sequence Variability Map viewer"), 4484));
      await session.step(29, "And the \"unique cliff pairs\" reading of Sequence Variability Map viewer should be 2242", () => readingIs(page, "unique cliff pairs", el("Sequence Variability Map viewer"), 2242));
      await session.step(30, "And the \"count of cell A at 2\" reading of Sequence Variability Map viewer should be 59", () => readingIs(page, "count of cell A at 2", el("Sequence Variability Map viewer"), 59));
      await session.step(31, "And the \"mean difference of cell A at 2\" reading of Sequence Variability Map viewer should be a finite number", () => readingFinite(page, "mean difference of cell A at 2", el("Sequence Variability Map viewer")));
      await session.step(32, "And the \"p-value of cell A at 2\" reading of Sequence Variability Map viewer should be between 0 and 1", () => readingBetween(page, "p-value of cell A at 2", el("Sequence Variability Map viewer"), 0, 1));
      await session.step(33, "And the \"cliffs of cell A at 2\" reading of Sequence Variability Map viewer should be 100", () => readingIs(page, "cliffs of cell A at 2", el("Sequence Variability Map viewer"), 100));
      await session.step(34, "And \"2\" column should have semantic type \"Monomer\"", () => columnSemType(page, "2", "Monomer"));
      await session.step(35, "And the table should have a column \"17\"", () => hasColumn(page, "17"));
      await session.step(36, "Then no errors should have been logged", () => noErrors(page));
      await session.step(37, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Switching from counts to mutation cliffs repaints the populated cell", async () => {
      await session.step(40, "When user clicks on \"Invariant Map\" checkbox in Sequence Variability Map viewer", () => clickOn(page, el("\"Invariant Map\" checkbox in Sequence Variability Map viewer")));
      await session.step(41, "Then the \"mode\" reading of Sequence Variability Map viewer should be \"Invariant Map\"", () => readingReads(page, "mode", el("Sequence Variability Map viewer"), "Invariant Map"));
      await session.step(42, "And the \"cell A at 2\" area of Sequence Variability Map viewer should be painted", () => areaPainted(page, "cell A at 2", el("Sequence Variability Map viewer")));
      await session.step(43, "And the \"cell A at 3\" area of Sequence Variability Map viewer should be painted in at least 1 colors", () => areaColors(page, "cell A at 3", el("Sequence Variability Map viewer"), 1));
      await session.step(44, "When user takes a snapshot of Sequence Variability Map viewer", () => takeSnapshot(page, el("Sequence Variability Map viewer")));
      await session.step(45, "When user clicks on \"Mutation Cliffs\" checkbox in Sequence Variability Map viewer", () => clickOn(page, el("\"Mutation Cliffs\" checkbox in Sequence Variability Map viewer")));
      await session.step(46, "Then the \"mode\" reading of Sequence Variability Map viewer should be \"Mutation Cliffs\"", () => readingReads(page, "mode", el("Sequence Variability Map viewer"), "Mutation Cliffs"));
      await session.step(47, "And the \"cell A at 2\" area of Sequence Variability Map viewer should have repainted", () => areaRepainted(page, "cell A at 2", el("Sequence Variability Map viewer")));
      await session.step(48, "And the \"cliffs of cell A at 2\" reading of Sequence Variability Map viewer should be 100", () => readingIs(page, "cliffs of cell A at 2", el("Sequence Variability Map viewer"), 100));
      await session.step(49, "And the \"cell A at 2\" area of Sequence Variability Map viewer should be painted in at least 1 colors", () => areaColors(page, "cell A at 2", el("Sequence Variability Map viewer"), 1));
      await session.step(50, "And the \"cliffs of cell A at 3\" reading of Sequence Variability Map viewer should be 0", () => readingIs(page, "cliffs of cell A at 3", el("Sequence Variability Map viewer"), 0));
      await session.step(51, "And the \"cell A at 3\" area of Sequence Variability Map viewer should be painted in no color", () => areaNoColor(page, "cell A at 3", el("Sequence Variability Map viewer")));
      await session.step(52, "Then no errors should have been logged", () => noErrors(page));
      await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The position chart draws precisely the peptides participating in position-two cliffs", async () => {
      await session.step(56, "Given user adds a Sequence Mutation Cliffs viewer with:", () => addViewerWith(page, "Sequence Mutation Cliffs", [["sequenceColumnName","AlignedSequence"],["activityColumnName","IC50"],["position","2"]]), [["sequenceColumnName","AlignedSequence"],["activityColumnName","IC50"],["position","2"]]);
      await session.step(60, "Then Sequence Mutation Cliffs viewer should be added to the open tableview", () => viewerAdded(page, "Sequence Mutation Cliffs"));
      await session.step(61, "And the \"position\" reading of Sequence Mutation Cliffs viewer should be 2", () => readingIs(page, "position", el("Sequence Mutation Cliffs viewer"), 2));
      await session.step(62, "And the \"cliff rows\" reading of Sequence Mutation Cliffs viewer should be 115", () => readingIs(page, "cliff rows", el("Sequence Mutation Cliffs viewer"), 115));
      await session.step(63, "And line chart viewer in Sequence Mutation Cliffs viewer should show 115 rows", () => showsRows(page, el("line chart viewer in Sequence Mutation Cliffs viewer"), 115));
      await session.step(64, "And the position 2 cliff chart should contain exactly its participating peptides", () => cliffParticipants(page, 2));
      await session.step(65, "And line chart viewer in Sequence Mutation Cliffs viewer should be painted", () => painted(page, el("line chart viewer in Sequence Mutation Cliffs viewer")));
      await session.step(66, "When user sets \"position\" property of Sequence Mutation Cliffs viewer to \"1\"", () => setProperty(page, "position", el("Sequence Mutation Cliffs viewer"), "1"));
      await session.step(67, "Then the \"cliff rows\" reading of Sequence Mutation Cliffs viewer should be 0", () => readingIs(page, "cliff rows", el("Sequence Mutation Cliffs viewer"), 0));
      await session.step(68, "And the \"message\" reading of Sequence Mutation Cliffs viewer should be \"No mutation cliffs found for the selected position.\"", () => readingReads(page, "message", el("Sequence Mutation Cliffs viewer"), "No mutation cliffs found for the selected position."));
      await session.step(69, "When user sets \"position\" property of Sequence Mutation Cliffs viewer to \"2\"", () => setProperty(page, "position", el("Sequence Mutation Cliffs viewer"), "2"));
      await session.step(70, "Then the \"cliff rows\" reading of Sequence Mutation Cliffs viewer should be 115", () => readingIs(page, "cliff rows", el("Sequence Mutation Cliffs viewer"), 115));
      await session.step(71, "And line chart viewer in Sequence Mutation Cliffs viewer should be painted", () => painted(page, el("line chart viewer in Sequence Mutation Cliffs viewer")));
      await session.step(72, "And the position 2 cliff chart should contain exactly its participating peptides", () => cliffParticipants(page, 2));
      await session.step(73, "Then no errors should have been logged", () => noErrors(page));
      await session.step(74, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Export contains every unique computed pair with its source sequences and activities", async () => {
      await session.step(77, "When user picks \"Export > Export Mutation Cliffs...\" from the context menu of Sequence Variability Map viewer", () => pickFromContextMenu(page, "Export > Export Mutation Cliffs...", el("Sequence Variability Map viewer")));
      await session.step(78, "Then \"Export Mutation Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Export Mutation Cliffs\" dialog"), "visible"));
      await session.step(79, "When user clicks on OK button in \"Export Mutation Cliffs\" dialog", () => clickOn(page, el("OK button in \"Export Mutation Cliffs\" dialog")));
      await session.step(80, "Then the \"Mutation Cliffs\" view should be current", () => viewIsCurrent(page, "Mutation Cliffs"));
      await session.step(81, "And table \"Mutation Cliffs\" should have columns \"Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta\"", () => tableColumns(page, "Mutation Cliffs", "Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta"));
      await session.step(82, "And the table should have 2242 rows", () => rowCount(page, 2242));
      await session.step(83, "And the mutation-cliff export should contain every single-mutation pair from table \"peptides\"", () => mutationPairs(page, "peptides"));
      await session.step(84, "When user closes the current view", () => closeCurrentView(page));
      await session.step(85, "And user switches to the \"peptides\" table view", () => switchTableView(page, "peptides"));
      await session.step(86, "Then no errors should have been logged", () => noErrors(page));
      await session.step(87, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The cluster summary accounts for all source peptides", async () => {
      await session.step(90, "Then MCL viewer should be added to the open tableview", () => viewerAdded(page, "MCL"));
      await session.step(91, "And Logo Summary Table viewer should be added to the open tableview", () => viewerAdded(page, "Logo Summary Table"));
      await session.step(92, "And the table should have a column \"Cluster (MCL)\"", () => hasColumn(page, "Cluster (MCL)"));
      await session.step(93, "And the \"clusters column\" reading of Logo Summary Table viewer should be \"Cluster (MCL)\"", () => readingReads(page, "clusters column", el("Logo Summary Table viewer"), "Cluster (MCL)"));
      await session.step(94, "And the \"members total\" reading of Logo Summary Table viewer should be 200", () => readingIs(page, "members total", el("Logo Summary Table viewer"), 200));
      await session.step(95, "And Logo Summary Table viewer should be painted", () => painted(page, el("Logo Summary Table viewer")));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
      await session.step(97, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
