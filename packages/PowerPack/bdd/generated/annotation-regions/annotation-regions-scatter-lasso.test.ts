/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/annotation-regions/annotation-regions-scatter-lasso.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, powerpack.dialogs.formula-lines]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hasNoArea, noErrors, pickFromAreaContextMenu, pickFromContextMenu, propertyShouldBe, propertyShouldContain, readingIs, readingReads, resizeTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragLasso} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A lasso drawn on the scatter plot becomes an area region", () => {
  const session = feature(test, "features/annotation-regions/annotation-regions-scatter-lasso.feature", import.meta.url);
  test("The lasso becomes a titled area region, removed again from the dialog", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["lassoTool","true"],["markerDefaultSize","2"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["lassoTool","true"],["markerDefaultSize","2"]]);
    await session.step(19, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(22, "Then the \"viewer regions\" reading of scatter plot viewer should be 0", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 0));
    await session.step(23, "When user picks \"Tools > Draw Annotation Region\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tools > Draw Annotation Region", el("scatter plot viewer")));
    await session.step(24, "Then the \"region drawing mode\" reading of scatter plot viewer should be \"true\"", () => readingReads(page, "region drawing mode", el("scatter plot viewer"), "true"));
    await session.step(25, "When user drags a lasso over the \"view\" area of scatter plot viewer", () => dragLasso(page, "view", el("scatter plot viewer")));
    await session.step(26, "Then the \"region drawing mode\" reading of scatter plot viewer should be \"false\"", () => readingReads(page, "region drawing mode", el("scatter plot viewer"), "false"));
    await session.step(27, "And the \"viewer regions\" reading of scatter plot viewer should be 1", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 1));
    await session.step(28, "And scatter plot viewer should have a \"region 1\" area", () => hasArea(page, el("scatter plot viewer"), "region 1"));
    await session.step(29, "And \"annotationRegions\" property of scatter plot viewer should contain \"\\\"type\\\":\\\"area\\\"\"", () => propertyShouldContain(page, "annotationRegions", el("scatter plot viewer"), "\"type\":\"area\""));
    await session.step(30, "And \"annotationRegions\" property of scatter plot viewer should contain \"\\\"x\\\":\\\"WEIGHT\\\",\\\"y\\\":\\\"HEIGHT\\\"\"", () => propertyShouldContain(page, "annotationRegions", el("scatter plot viewer"), "\"x\":\"WEIGHT\",\"y\":\"HEIGHT\""));
    await session.step(31, "And \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(32, "And editor of \"X column\" input in \"Formula Lines\" dialog should have text \"WEIGHT\"", () => shouldHaveText(page, el("editor of \"X column\" input in \"Formula Lines\" dialog"), "WEIGHT"));
    await session.step(33, "And editor of \"Y column\" input in \"Formula Lines\" dialog should have text \"HEIGHT\"", () => shouldHaveText(page, el("editor of \"Y column\" input in \"Formula Lines\" dialog"), "HEIGHT"));
    await session.step(34, "When user enters \"Lasso A\" into Title input in \"Formula Lines\" dialog", () => enterInto(page, "Lasso A", el("Title input in \"Formula Lines\" dialog")));
    await session.step(35, "And user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(36, "Then the \"Formula Lines\" dialog should close", () => dialogCloses(page, "Formula Lines"));
    await session.step(37, "And scatter plot viewer should have a \"region Lasso A\" area", () => hasArea(page, el("scatter plot viewer"), "region Lasso A"));
    await session.step(38, "And the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
    await session.step(39, "And \"Lasso Tool\" property of scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "Lasso Tool", el("scatter plot viewer"), "true"));
    await session.step(40, "When user picks \"Tools > Formula Lines...\" from the context menu of the \"bottom edge of view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tools > Formula Lines...", "bottom edge of view", el("scatter plot viewer")));
    await session.step(41, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(42, "When user clicks on Delete button in \"Formula Lines\" dialog", () => clickOn(page, el("Delete button in \"Formula Lines\" dialog")));
    await session.step(43, "And user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(44, "Then the \"viewer regions\" reading of scatter plot viewer should be 0", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 0));
    await session.step(45, "And scatter plot viewer should not have a \"region Lasso A\" area", () => hasNoArea(page, el("scatter plot viewer"), "region Lasso A"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
});
