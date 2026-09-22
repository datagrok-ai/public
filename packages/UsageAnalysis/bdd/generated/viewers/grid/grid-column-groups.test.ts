/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-column-groups.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, isExpanded, pressKey, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {columnsSelected, noColumnsSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, contextPanelOpen, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColor, areaNotColor, clickArea, clickAreaHolding, hasArea, hasNoArea, loadLayout, noErrors, saveLayoutToServer, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid column groups", () => {
  const session = feature(test, "features/viewers/grid/grid-column-groups.feature", import.meta.url);
  test("Grid column groups", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(55, "Given user is logged in", () => loggedIn(page));
    await session.step(56, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(57, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(58, "And grid should not have a \"group Person\" area", () => hasNoArea(page, el("grid"), "group Person"));
    await run.scenario("Group columns... from the Context Panel draws a band in the group's colour", async () => {
      await session.step(61, "When user clicks on the \"header AGE\" area of grid holding Control", () => clickAreaHolding(page, "header AGE", el("grid"), "Control"));
      await session.step(62, "And user clicks on the \"header SEX\" area of grid holding Control", () => clickAreaHolding(page, "header SEX", el("grid"), "Control"));
      await session.step(63, "Then columns \"AGE, SEX\" should be selected", () => columnsSelected(page, "AGE, SEX"));
      await session.step(64, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(65, "And Items accordion header in context panel is expanded", () => isExpanded(page, el("Items accordion header in context panel")));
      await session.step(66, "Then context panel should contain text \"AGE\"", () => shouldContainText(page, el("context panel"), "AGE"));
      await session.step(67, "And context panel should contain text \"SEX\"", () => shouldContainText(page, el("context panel"), "SEX"));
      await session.step(68, "Given Actions accordion header in context panel is expanded", () => isExpanded(page, el("Actions accordion header in context panel")));
      await session.step(69, "When user clicks on \"Group columns...\" text in context panel", () => clickOn(page, el("\"Group columns...\" text in context panel")));
      await session.step(70, "Then Group input in dialog should be visible", () => shouldBe(page, el("Group input in dialog"), "visible"));
      await session.step(71, "When user types \"Person\" into Group input in dialog", () => typeInto(page, "Person", el("Group input in dialog")));
      await session.step(72, "And user types \"#FF0000\" into Color input in dialog", () => typeInto(page, "#FF0000", el("Color input in dialog")));
      await session.step(73, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(74, "Then dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(75, "And grid should have a \"group Person\" area", () => hasArea(page, el("grid"), "group Person"));
      await session.step(76, "And \"AGE\" column should have tag \"group\" equal to \"Person\"", () => columnTag(page, "AGE", "group", "Person"));
      await session.step(77, "And \"SEX\" column should have tag \"group\" equal to \"Person\"", () => columnTag(page, "SEX", "group", "Person"));
      await session.step(78, "And the \"group Person\" area of grid should contain the color \"#FF0000\"", () => areaColor(page, "group Person", el("grid"), "#FF0000"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shift+click across the group, a click on its band and Escape raise no error", async () => {
      await session.step(82, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
      await session.step(83, "And user clicks on the \"header SEX\" area of grid holding Shift", () => clickAreaHolding(page, "header SEX", el("grid"), "Shift"));
      await session.step(84, "Then no errors should have been logged", () => noErrors(page));
      await session.step(85, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(86, "Then no columns should be selected", () => noColumnsSelected(page));
      await session.step(87, "When user clicks on the \"group Person\" area of grid", () => clickArea(page, "group Person", el("grid")));
      await session.step(88, "Then columns \"AGE, SEX\" should be selected", () => columnsSelected(page, "AGE, SEX"));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
      await session.step(90, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(91, "Then no columns should be selected", () => noColumnsSelected(page));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout saved to the server brings the group back over a fresh view of the table", async () => {
      await session.step(95, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(96, "And user closes all views", () => closeAllViews(page));
      await session.step(97, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(98, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
      await session.step(99, "And grid should not have a \"group Person\" area", () => hasNoArea(page, el("grid"), "group Person"));
      await session.step(100, "And \"AGE\" column should have tag \"group\" equal to \"\"", () => columnTag(page, "AGE", "group", ""));
      await session.step(101, "When user loads the saved layout", () => loadLayout(page));
      await session.step(102, "Then grid should have a \"group Person\" area", () => hasArea(page, el("grid"), "group Person"));
      await session.step(103, "And \"AGE\" column should have tag \"group\" equal to \"Person\"", () => columnTag(page, "AGE", "group", "Person"));
      await session.step(104, "And \"SEX\" column should have tag \"group\" equal to \"Person\"", () => columnTag(page, "SEX", "group", "Person"));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second group gets its own band and colour; the first one's band selected and Escape raise no error", async () => {
      await session.step(108, "When user clicks on the \"header HEIGHT\" area of grid holding Control", () => clickAreaHolding(page, "header HEIGHT", el("grid"), "Control"));
      await session.step(109, "And user clicks on the \"header WEIGHT\" area of grid holding Control", () => clickAreaHolding(page, "header WEIGHT", el("grid"), "Control"));
      await session.step(110, "Then columns \"HEIGHT, WEIGHT\" should be selected", () => columnsSelected(page, "HEIGHT, WEIGHT"));
      await session.step(111, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(112, "Then context panel should contain text \"HEIGHT\"", () => shouldContainText(page, el("context panel"), "HEIGHT"));
      await session.step(113, "And context panel should contain text \"WEIGHT\"", () => shouldContainText(page, el("context panel"), "WEIGHT"));
      await session.step(114, "Given Actions accordion header in context panel is expanded", () => isExpanded(page, el("Actions accordion header in context panel")));
      await session.step(115, "When user clicks on \"Group columns...\" text in context panel", () => clickOn(page, el("\"Group columns...\" text in context panel")));
      await session.step(116, "And user types \"Body\" into Group input in dialog", () => typeInto(page, "Body", el("Group input in dialog")));
      await session.step(117, "And user types \"#0000FF\" into Color input in dialog", () => typeInto(page, "#0000FF", el("Color input in dialog")));
      await session.step(118, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(119, "Then dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(120, "And grid should have a \"group Body\" area", () => hasArea(page, el("grid"), "group Body"));
      await session.step(121, "And grid should have a \"group Person\" area", () => hasArea(page, el("grid"), "group Person"));
      await session.step(122, "And \"HEIGHT\" column should have tag \"group\" equal to \"Body\"", () => columnTag(page, "HEIGHT", "group", "Body"));
      await session.step(123, "And \"AGE\" column should have tag \"group\" equal to \"Person\"", () => columnTag(page, "AGE", "group", "Person"));
      await session.step(124, "And the \"group Body\" area of grid should contain the color \"#0000FF\"", () => areaColor(page, "group Body", el("grid"), "#0000FF"));
      await session.step(125, "And the \"group Body\" area of grid should not contain the color \"#FF0000\"", () => areaNotColor(page, "group Body", el("grid"), "#FF0000"));
      await session.step(126, "When user clicks on the \"group Person\" area of grid", () => clickArea(page, "group Person", el("grid")));
      await session.step(127, "Then columns \"AGE, SEX\" should be selected", () => columnsSelected(page, "AGE, SEX"));
      await session.step(128, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(129, "Then no columns should be selected", () => noColumnsSelected(page));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout saved to the server restores the groups, their columns and their colours", async () => {
      await session.step(133, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(134, "And user closes all views", () => closeAllViews(page));
      await session.step(135, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(136, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
      await session.step(137, "And grid should not have a \"group Person\" area", () => hasNoArea(page, el("grid"), "group Person"));
      await session.step(138, "And grid should not have a \"group Body\" area", () => hasNoArea(page, el("grid"), "group Body"));
      await session.step(139, "And \"AGE\" column should have tag \"group\" equal to \"\"", () => columnTag(page, "AGE", "group", ""));
      await session.step(140, "When user loads the saved layout", () => loadLayout(page));
      await session.step(141, "Then grid should have a \"group Person\" area", () => hasArea(page, el("grid"), "group Person"));
      await session.step(142, "And grid should have a \"group Body\" area", () => hasArea(page, el("grid"), "group Body"));
      await session.step(143, "And \"AGE\" column should have tag \"group\" equal to \"Person\"", () => columnTag(page, "AGE", "group", "Person"));
      await session.step(144, "And \"SEX\" column should have tag \"group\" equal to \"Person\"", () => columnTag(page, "SEX", "group", "Person"));
      await session.step(145, "And \"HEIGHT\" column should have tag \"group\" equal to \"Body\"", () => columnTag(page, "HEIGHT", "group", "Body"));
      await session.step(146, "And \"WEIGHT\" column should have tag \"group\" equal to \"Body\"", () => columnTag(page, "WEIGHT", "group", "Body"));
      await session.step(147, "And the \"group Person\" area of grid should contain the color \"#FF0000\"", () => areaColor(page, "group Person", el("grid"), "#FF0000"));
      await session.step(148, "And the \"group Body\" area of grid should contain the color \"#0000FF\"", () => areaColor(page, "group Body", el("grid"), "#0000FF"));
      await session.step(149, "And the \"group Person\" area of grid should not contain the color \"#0000FF\"", () => areaNotColor(page, "group Person", el("grid"), "#0000FF"));
      await session.step(150, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Both groups and their columns come back from a project", async () => {
      await session.step(153, "When user saves the current view as project \"zz-grid-column-groups-2\"", () => saveAsProject(page, "zz-grid-column-groups-2"));
      await session.step(154, "And user closes all views", () => closeAllViews(page));
      await session.step(155, "And user opens the \"zz-grid-column-groups-2\" project", () => openProject(page, "zz-grid-column-groups-2"));
      await session.step(156, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
      await session.step(157, "And grid should have a \"group Person\" area", () => hasArea(page, el("grid"), "group Person"));
      await session.step(158, "And grid should have a \"group Body\" area", () => hasArea(page, el("grid"), "group Body"));
      await session.step(159, "And \"AGE\" column should have tag \"group\" equal to \"Person\"", () => columnTag(page, "AGE", "group", "Person"));
      await session.step(160, "And \"WEIGHT\" column should have tag \"group\" equal to \"Body\"", () => columnTag(page, "WEIGHT", "group", "Body"));
      await session.step(161, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Grouping the same columns again gives the band a new name and colour", async () => {
      await session.step(164, "When user clicks on the \"group Person\" area of grid", () => clickArea(page, "group Person", el("grid")));
      await session.step(165, "Then columns \"AGE, SEX\" should be selected", () => columnsSelected(page, "AGE, SEX"));
      await session.step(166, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(167, "Then context panel should contain text \"AGE\"", () => shouldContainText(page, el("context panel"), "AGE"));
      await session.step(168, "And context panel should contain text \"SEX\"", () => shouldContainText(page, el("context panel"), "SEX"));
      await session.step(169, "Given Actions accordion header in context panel is expanded", () => isExpanded(page, el("Actions accordion header in context panel")));
      await session.step(170, "When user clicks on \"Group columns...\" text in context panel", () => clickOn(page, el("\"Group columns...\" text in context panel")));
      await session.step(171, "And user types \"Measures\" into Group input in dialog", () => typeInto(page, "Measures", el("Group input in dialog")));
      await session.step(172, "And user types \"#d62728\" into Color input in dialog", () => typeInto(page, "#d62728", el("Color input in dialog")));
      await session.step(173, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(174, "Then dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(175, "And grid should have a \"group Measures\" area", () => hasArea(page, el("grid"), "group Measures"));
      await session.step(176, "And grid should not have a \"group Person\" area", () => hasNoArea(page, el("grid"), "group Person"));
      await session.step(177, "And the \"group Measures\" area of grid should contain the color \"#d62728\"", () => areaColor(page, "group Measures", el("grid"), "#d62728"));
      await session.step(178, "And \"SEX\" column should have tag \"group\" equal to \"Measures\"", () => columnTag(page, "SEX", "group", "Measures"));
      await session.step(179, "And \"AGE\" column should have tag \"group\" equal to \"Measures\"", () => columnTag(page, "AGE", "group", "Measures"));
      await session.step(180, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(181, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Grouping with an empty name ungroups the columns", async () => {
      await session.step(184, "When user clicks on the \"group Measures\" area of grid", () => clickArea(page, "group Measures", el("grid")));
      await session.step(185, "Then columns \"AGE, SEX\" should be selected", () => columnsSelected(page, "AGE, SEX"));
      await session.step(186, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(187, "Then context panel should contain text \"AGE\"", () => shouldContainText(page, el("context panel"), "AGE"));
      await session.step(188, "Given Actions accordion header in context panel is expanded", () => isExpanded(page, el("Actions accordion header in context panel")));
      await session.step(189, "When user clicks on \"Group columns...\" text in context panel", () => clickOn(page, el("\"Group columns...\" text in context panel")));
      await session.step(190, "And user clears Group input in dialog", () => clearField(page, el("Group input in dialog")));
      await session.step(191, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(192, "Then dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(193, "And grid should not have a \"group Measures\" area", () => hasNoArea(page, el("grid"), "group Measures"));
      await session.step(194, "And grid should have a \"group Body\" area", () => hasArea(page, el("grid"), "group Body"));
      await session.step(195, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(196, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
