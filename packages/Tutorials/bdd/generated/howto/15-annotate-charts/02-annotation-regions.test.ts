/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/15-annotate-charts/02-annotation-regions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {clickOnAxes, dragOnAxes, pointAtOnAxes} from '../../../bindings/howto-annotations.js';
import {viewerShare} from '../../../bindings/howto-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, shouldBe, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {noneOfSelected, selectedRowCount, someOfSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, noErrors, pickFromAreaContextMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Annotations", () => {
  const session = feature(test, "features/howto/15-annotate-charts/02-annotation-regions.feature", import.meta.url);
  test("Mark a group of points with an annotation region", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens iris dataset", () => openDataset(page, ds("iris")));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","Petal.Length"],["Y","Petal.Width"],["Color","Species"]]), [["X","Petal.Length"],["Y","Petal.Width"],["Color","Species"]]);
    await session.step(25, "And scatter plot viewer takes 65% of the space it shares with the grid", () => viewerShare(page, el("scatter plot viewer"), 65));
    await session.step(26, "When the video says \"An annotation region marks and names a part of a chart\"", () => videoSays(page, "An annotation region marks and names a part of a chart"));
    await session.step(28, "And user picks \"Tools > Draw Annotation Region\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tools > Draw Annotation Region", "view", el("scatter plot viewer")), undefined, "Right-click the plot: Tools > Draw Annotation Region");
    await session.step(29, "Then the \"region drawing mode\" reading of scatter plot viewer should be \"true\"", () => readingReads(page, "region drawing mode", el("scatter plot viewer"), "true"));
    await session.step(31, "When user drags across scatter plot viewer from 4.45, 1.85 to 5.15, 1.35 on its axes", () => dragOnAxes(page, el("scatter plot viewer"), 4.45, 1.85, 5.15, 1.35), undefined, "Drag a box around the points to mark");
    await session.step(32, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(34, "When user enters \"Overlap\" into \"Title\" text input in \"Formula Lines\" dialog", () => enterInto(page, "Overlap", el("\"Title\" text input in \"Formula Lines\" dialog")), undefined, "Give it a name and a note");
    await session.step(35, "And user enters \"Versicolor and virginica mix here\" into \"Description\" text area in \"Formula Lines\" dialog", () => enterInto(page, "Versicolor and virginica mix here", el("\"Description\" text area in \"Formula Lines\" dialog")));
    await session.step(36, "And user clicks on OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(37, "Then scatter plot viewer should have a 'region \"Overlap\"' area", () => hasArea(page, el("scatter plot viewer"), "region \"Overlap\""));
    await session.step(38, "And the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
    await session.step(40, "When user points at 4.82, 1.45 inside the 'region \"Overlap\"' area of scatter plot viewer", () => pointAtOnAxes(page, 4.82, 1.45, "region \"Overlap\"", el("scatter plot viewer")), undefined, "Point at it: the name, the note and how many rows it holds");
    await session.step(41, "Then tooltip should contain text \"26 rows\"", () => shouldContainText(page, el("tooltip"), "26 rows"));
    await session.step(43, "When user clicks at 4.82, 1.45 inside the 'region \"Overlap\"' area of scatter plot viewer", () => clickOnAxes(page, 4.82, 1.45, "region \"Overlap\"", el("scatter plot viewer")), undefined, "A click selects the rows inside");
    await session.step(44, "Then 26 rows should be selected", () => selectedRowCount(page, 26));
    await session.step(45, "And some rows where \"Species\" is \"versicolor\" should be selected", () => someOfSelected(page, "Species", "versicolor"));
    await session.step(46, "And some rows where \"Species\" is \"virginica\" should be selected", () => someOfSelected(page, "Species", "virginica"));
    await session.step(47, "And no rows where \"Species\" is \"setosa\" should be selected", () => noneOfSelected(page, "Species", "setosa"));
    await session.step(49, "When user hovers over \"Selected: 26\" text in status bar", () => hoverOver(page, el("\"Selected: 26\" text in status bar")), undefined, "The status bar counts them: 26 selected");
    await session.step(50, "Then no errors should have been logged", () => noErrors(page));
  });
});
