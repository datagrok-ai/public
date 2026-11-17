/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-scatter-labels.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot, viewers.scatter-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {noneSelected, onlyOfSelected, selectNoRows, selectWhereBetween, selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, clickArea, eventFired, hoverArea, listenFor, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {innerPropertyShouldBe, noSuchReading, setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot labels in trellis cells", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-scatter-labels.feature", import.meta.url);
  test("Scatter plot labels in trellis cells", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(23, "And user sets \"labelColumnNames\" inner property of trellis plot viewer to \"USUBJID\"", () => setInnerProperty(page, "labelColumnNames", el("trellis plot viewer"), "USUBJID"));
    await session.step(24, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await session.step(25, "And \"showLabelsFor\" inner property of trellis plot viewer should be \"FilteredSelected\"", () => innerPropertyShouldBe(page, "showLabelsFor", el("trellis plot viewer"), "FilteredSelected"));
    await session.step(26, "And the \"cell labels shown F | Caucasian\" reading of trellis plot viewer should be 0", () => readingIs(page, "cell labels shown F | Caucasian", el("trellis plot viewer"), 0));
    await session.step(27, "And the \"cell labels shown F | Other\" reading of trellis plot viewer should be 0", () => readingIs(page, "cell labels shown F | Other", el("trellis plot viewer"), 0));
    await run.scenario("A cell labels the selected rows it draws, whatever Show Labels For says", async () => {
      await session.step(30, "When user selects rows where \"AGE\" is between 20 and 25", () => selectWhereBetween(page, "AGE", 20, 25));
      await session.step(31, "Then the \"cell labels shown F | Black\" reading of trellis plot viewer should be 1", () => readingIs(page, "cell labels shown F | Black", el("trellis plot viewer"), 1));
      await session.step(32, "And the \"cell labels shown F | Other\" reading of trellis plot viewer should be 3", () => readingIs(page, "cell labels shown F | Other", el("trellis plot viewer"), 3));
      await session.step(33, "And the \"cell labels shown M | Black\" reading of trellis plot viewer should be 0", () => readingIs(page, "cell labels shown M | Black", el("trellis plot viewer"), 0));
      await session.step(34, "And the \"cell labels shown F | Caucasian\" reading of trellis plot viewer should be 20", () => readingIs(page, "cell labels shown F | Caucasian", el("trellis plot viewer"), 20));
      await session.step(35, "When user sets \"showLabelsFor\" inner property of trellis plot viewer to \"Filtered\"", () => setInnerProperty(page, "showLabelsFor", el("trellis plot viewer"), "Filtered"));
      await session.step(36, "Then the \"cell labels shown F | Other\" reading of trellis plot viewer should be 3", () => readingIs(page, "cell labels shown F | Other", el("trellis plot viewer"), 3));
      await session.step(37, "And the \"cell labels shown M | Black\" reading of trellis plot viewer should be 0", () => readingIs(page, "cell labels shown M | Black", el("trellis plot viewer"), 0));
      await session.step(38, "When user sets \"showLabelsFor\" inner property of trellis plot viewer to \"All\"", () => setInnerProperty(page, "showLabelsFor", el("trellis plot viewer"), "All"));
      await session.step(39, "Then the \"cell labels shown F | Other\" reading of trellis plot viewer should be 3", () => readingIs(page, "cell labels shown F | Other", el("trellis plot viewer"), 3));
      await session.step(40, "When user selects no rows", () => selectNoRows(page));
      await session.step(41, "Then the \"cell labels shown F | Other\" reading of trellis plot viewer should be 0", () => readingIs(page, "cell labels shown F | Other", el("trellis plot viewer"), 0));
      await session.step(42, "When user sets \"showLabelsFor\" inner property of trellis plot viewer to \"FilteredSelected\"", () => setInnerProperty(page, "showLabelsFor", el("trellis plot viewer"), "FilteredSelected"));
      await session.step(43, "Then the \"cell labels shown F | Black\" reading of trellis plot viewer should be 0", () => readingIs(page, "cell labels shown F | Black", el("trellis plot viewer"), 0));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cell shows at most 20 labels, with Display Labels Auto or Always, and none with Never", async () => {
      await session.step(47, "When user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
      await session.step(48, "Then the \"cell labels shown F | Caucasian\" reading of trellis plot viewer should be 20", () => readingIs(page, "cell labels shown F | Caucasian", el("trellis plot viewer"), 20));
      await session.step(49, "And the \"cell labels shown F | Black\" reading of trellis plot viewer should be 17", () => readingIs(page, "cell labels shown F | Black", el("trellis plot viewer"), 17));
      await session.step(50, "And the \"cell labels shown M | Caucasian\" reading of trellis plot viewer should be 0", () => readingIs(page, "cell labels shown M | Caucasian", el("trellis plot viewer"), 0));
      await session.step(51, "When user sets \"displayLabels\" inner property of trellis plot viewer to \"Always\"", () => setInnerProperty(page, "displayLabels", el("trellis plot viewer"), "Always"));
      await session.step(52, "Then the \"cell labels shown F | Caucasian\" reading of trellis plot viewer should be 20", () => readingIs(page, "cell labels shown F | Caucasian", el("trellis plot viewer"), 20));
      await session.step(53, "And the \"cell labels shown F | Black\" reading of trellis plot viewer should be 17", () => readingIs(page, "cell labels shown F | Black", el("trellis plot viewer"), 17));
      await session.step(54, "When user sets \"displayLabels\" inner property of trellis plot viewer to \"Never\"", () => setInnerProperty(page, "displayLabels", el("trellis plot viewer"), "Never"));
      await session.step(55, "Then the \"cell labels shown F | Caucasian\" reading of trellis plot viewer should be 0", () => readingIs(page, "cell labels shown F | Caucasian", el("trellis plot viewer"), 0));
      await session.step(56, "And the \"cell labels shown F | Black\" reading of trellis plot viewer should be 0", () => readingIs(page, "cell labels shown F | Black", el("trellis plot viewer"), 0));
      await session.step(57, "When user sets \"displayLabels\" inner property of trellis plot viewer to \"Auto\"", () => setInnerProperty(page, "displayLabels", el("trellis plot viewer"), "Auto"));
      await session.step(58, "Then the \"cell labels shown F | Caucasian\" reading of trellis plot viewer should be 20", () => readingIs(page, "cell labels shown F | Caucasian", el("trellis plot viewer"), 20));
      await session.step(59, "When user selects no rows", () => selectNoRows(page));
      await session.step(60, "Then the \"cell labels shown F | Caucasian\" reading of trellis plot viewer should be 0", () => readingIs(page, "cell labels shown F | Caucasian", el("trellis plot viewer"), 0));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The cell opened full screen labels by its own Show Labels For, with no cap, and keeps the selection", async () => {
      await session.step(64, "When user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
      await session.step(65, "And user sets \"showLabelsFor\" inner property of trellis plot viewer to \"Selected\"", () => setInnerProperty(page, "showLabelsFor", el("trellis plot viewer"), "Selected"));
      await session.step(66, "Then the \"cell labels shown F | Caucasian\" reading of trellis plot viewer should be 20", () => readingIs(page, "cell labels shown F | Caucasian", el("trellis plot viewer"), 20));
      await session.step(67, "When user hovers over the \"cell body F | Caucasian\" area of trellis plot viewer", () => hoverArea(page, "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(68, "And user clicks on the \"full screen icon\" area of trellis plot viewer", () => clickArea(page, "full screen icon", el("trellis plot viewer")));
      await session.step(69, "Then \"SEX: F, RACE: Caucasian\" dialog should be visible", () => shouldBe(page, el("\"SEX: F, RACE: Caucasian\" dialog"), "visible"));
      await session.step(70, "And the \"full screen labels shown\" reading of trellis plot viewer should be 429", () => readingIs(page, "full screen labels shown", el("trellis plot viewer"), 429));
      await session.step(71, "When user listens for \"d4-trellis-plot-inner-viewer-clicked\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-inner-viewer-clicked", el("trellis plot viewer")));
      await session.step(72, "And user clicks on the \"top left corner of full screen view\" area of trellis plot viewer", () => clickArea(page, "top left corner of full screen view", el("trellis plot viewer")));
      await session.step(73, "Then \"d4-trellis-plot-inner-viewer-clicked\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-inner-viewer-clicked", el("trellis plot viewer")));
      await session.step(74, "And only rows where \"SEX\" is \"F\" should be selected", () => onlyOfSelected(page, "SEX", "F"));
      await session.step(75, "And the \"full screen labels shown\" reading of trellis plot viewer should be 429", () => readingIs(page, "full screen labels shown", el("trellis plot viewer"), 429));
      await session.step(76, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(77, "Then \"SEX: F, RACE: Caucasian\" dialog should be absent", () => shouldBe(page, el("\"SEX: F, RACE: Caucasian\" dialog"), "absent"));
      await session.step(78, "And trellis plot viewer should not report a \"full screen labels shown\" reading", () => noSuchReading(page, el("trellis plot viewer"), "full screen labels shown"));
      await session.step(79, "When user selects no rows", () => selectNoRows(page));
      await session.step(80, "And user sets \"showLabelsFor\" inner property of trellis plot viewer to \"FilteredSelected\"", () => setInnerProperty(page, "showLabelsFor", el("trellis plot viewer"), "FilteredSelected"));
      await session.step(81, "Then the \"cell labels shown F | Caucasian\" reading of trellis plot viewer should be 0", () => readingIs(page, "cell labels shown F | Caucasian", el("trellis plot viewer"), 0));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on empty space in a cell keeps the selection a standalone scatter plot clears", async () => {
      await session.step(85, "Given user listens for \"d4-trellis-plot-inner-viewer-clicked\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-inner-viewer-clicked", el("trellis plot viewer")));
      await session.step(86, "When user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
      await session.step(87, "And user clicks on the \"top left corner of cell body F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "top left corner of cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(88, "Then \"d4-trellis-plot-inner-viewer-clicked\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-inner-viewer-clicked", el("trellis plot viewer")));
      await session.step(89, "And only rows where \"SEX\" is \"F\" should be selected", () => onlyOfSelected(page, "SEX", "F"));
      await session.step(90, "And the \"cell labels shown F | Caucasian\" reading of trellis plot viewer should be 20", () => readingIs(page, "cell labels shown F | Caucasian", el("trellis plot viewer"), 20));
      await session.step(91, "When user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(92, "And user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(93, "Then no rows should be selected", () => noneSelected(page));
      await session.step(94, "And the \"cell labels shown F | Caucasian\" reading of trellis plot viewer should be 0", () => readingIs(page, "cell labels shown F | Caucasian", el("trellis plot viewer"), 0));
      await session.step(95, "When user clicks on close icon of scatter plot viewer", () => clickOn(page, el("close icon of scatter plot viewer")));
      await session.step(96, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
