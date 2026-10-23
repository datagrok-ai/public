/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/formula-lines/formula-lines-axis-annotations-box-plot.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot, powerpack.dialogs.formula-lines]
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
import {clickOn, shouldBe, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areasSameSize, closeContextMenu, hasArea, hasNoArea, menuDoesNotList, menuLists, noErrors, pickFromAreaContextMenu, propertyShouldContain, readingIs, resizeTo, rightClickArea, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Annotations group on the box plot's value axis", () => {
  const session = feature(test, "features/formula-lines/formula-lines-axis-annotations-box-plot.feature", import.meta.url);
  test("The value axis adds a line, a band and a region on AGE", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["categoryColumnNames","RACE"],["valueColumnName","AGE"]]), [["categoryColumnNames","RACE"],["valueColumnName","AGE"]]);
    await session.step(20, "And user resizes box plot viewer to 800 by 500", () => resizeTo(page, el("box plot viewer"), 800, 500));
    await session.step(23, "When user picks \"Annotations > Add Line\" from the context menu of the \"bottom edge of y axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Line", "bottom edge of y axis", el("box plot viewer")));
    await session.step(24, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(25, "And editor of Column input in \"Formula Lines\" dialog should have text \"AGE\"", () => shouldHaveText(page, el("editor of Column input in \"Formula Lines\" dialog"), "AGE"));
    await session.step(26, "And Value input in \"Formula Lines\" dialog should have value \"45\"", () => shouldHaveValue(page, el("Value input in \"Formula Lines\" dialog"), "45"));
    await session.step(27, "When user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(28, "Then the \"Formula Lines\" dialog should close", () => dialogCloses(page, "Formula Lines"));
    await session.step(29, "And \"formulaLines\" property of box plot viewer should contain \"${AGE} = 45.0\"", () => propertyShouldContain(page, "formulaLines", el("box plot viewer"), "${AGE} = 45.0"));
    await session.step(30, "When user sets \"formulaLines\" property of box plot viewer to \"\"", () => setProperty(page, "formulaLines", el("box plot viewer"), ""));
    await session.step(31, "And user picks \"Annotations > Add Band\" from the context menu of the \"bottom edge of y axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Band", "bottom edge of y axis", el("box plot viewer")));
    await session.step(32, "And user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(33, "Then \"formulaLines\" property of box plot viewer should contain \"${AGE} in (36.0, 56.0)\"", () => propertyShouldContain(page, "formulaLines", el("box plot viewer"), "${AGE} in (36.0, 56.0)"));
    await session.step(34, "When user sets \"formulaLines\" property of box plot viewer to \"\"", () => setProperty(page, "formulaLines", el("box plot viewer"), ""));
    await session.step(35, "Then the \"viewer regions\" reading of box plot viewer should be 0", () => readingIs(page, "viewer regions", el("box plot viewer"), 0));
    await session.step(36, "When user picks \"Annotations > Add Region\" from the context menu of the \"bottom edge of y axis\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Region", "bottom edge of y axis", el("box plot viewer")));
    await session.step(37, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(38, "When user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(39, "Then \"annotationRegions\" property of box plot viewer should contain \"${AGE} = 36.0\"", () => propertyShouldContain(page, "annotationRegions", el("box plot viewer"), "${AGE} = 36.0"));
    await session.step(40, "And \"annotationRegions\" property of box plot viewer should contain \"${AGE} = 56.0\"", () => propertyShouldContain(page, "annotationRegions", el("box plot viewer"), "${AGE} = 56.0"));
    await session.step(41, "And the \"viewer regions\" reading of box plot viewer should be 1", () => readingIs(page, "viewer regions", el("box plot viewer"), 1));
    await session.step(42, "And box plot viewer should have a \"region 1\" area", () => hasArea(page, el("box plot viewer"), "region 1"));
    await session.step(43, "And the \"region 1\" and \"view\" areas of box plot viewer should be the same width", () => areasSameSize(page, "region 1", "view", el("box plot viewer"), "width"));
    await session.step(44, "When user sets \"annotationRegions\" property of box plot viewer to \"[]\"", () => setProperty(page, "annotationRegions", el("box plot viewer"), "[]"));
    await session.step(45, "Then box plot viewer should not have a \"region 1\" area", () => hasNoArea(page, el("box plot viewer"), "region 1"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
  test("The category axis has no Annotations group", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:powerpack.dialogs.formula-lines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["categoryColumnNames","RACE"],["valueColumnName","AGE"]]), [["categoryColumnNames","RACE"],["valueColumnName","AGE"]]);
    await session.step(20, "And user resizes box plot viewer to 800 by 500", () => resizeTo(page, el("box plot viewer"), 800, 500));
    await session.step(49, "When user right-clicks on the \"x axis\" area of box plot viewer", () => rightClickArea(page, "x axis", el("box plot viewer")));
    await session.step(50, "Then the open menu should list \"Show Category Axis\"", () => menuLists(page, "Show Category Axis"));
    await session.step(51, "And the open menu should not list \"Annotations\"", () => menuDoesNotList(page, "Annotations"));
    await session.step(52, "And the open menu should not list \"Add Line\"", () => menuDoesNotList(page, "Add Line"));
    await session.step(53, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(54, "And user right-clicks on the \"y axis\" area of box plot viewer", () => rightClickArea(page, "y axis", el("box plot viewer")));
    await session.step(55, "Then the open menu should list \"Annotations\"", () => menuLists(page, "Annotations"));
    await session.step(56, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(57, "Then no errors should have been logged", () => noErrors(page));
  });
});
