/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewer/properties.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [dendrogram.cp.viewer-from-newick-prop]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, clickOn, enterInto, expand, pressKeyIn, selectIn, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {setTableTag} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, noErrors, pickColorSwatch, pickFromContextMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Dendrogram viewer's properties from the property panel", () => {
  const session = feature(test, "features/viewer/properties.feature", import.meta.url);
  test("The Dendrogram viewer's properties from the property panel", {tag: ["@journey", "@realizes:dendrogram.cp.viewer-from-newick-prop", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens a table \"leaves4\" with:", () => openTableOf(page, "leaves4", [["leaf","value"],["a","1"],["b","2"],["c","3"],["d","4"]]), [["leaf","value"],["a","1"],["b","2"],["c","3"],["d","4"]]);
    await session.step(24, "And user sets tag \".newick\" of the table to \"((a:1,b:1):1,(c:1,d:1):1);\"", () => setTableTag(page, ".newick", "((a:1,b:1):1,(c:1,d:1):1);"));
    await session.step(25, "And user sets tag \".newick-alt\" of the table to \"(((a:1,b:1):1,c:1):1,d:1);\"", () => setTableTag(page, ".newick-alt", "(((a:1,b:1):1,c:1):1,d:1);"));
    await session.step(26, "And user adds Dendrogram viewer with:", () => addViewerWith(page, "Dendrogram", [["newick","((a:1,b:1):1,(c:1,d:1):1);"],["nodeColumnName","leaf"]]), [["newick","((a:1,b:1):1,(c:1,d:1):1);"],["nodeColumnName","leaf"]]);
    await session.step(29, "Then the \"newick\" reading of Dendrogram viewer should be \"((a:1,b:1):1,(c:1,d:1):1);\"", () => readingReads(page, "newick", el("Dendrogram viewer"), "((a:1,b:1):1,(c:1,d:1):1);"));
    await session.step(30, "And the \"leaves\" reading of Dendrogram viewer should be \"a, b, c, d\"", () => readingReads(page, "leaves", el("Dendrogram viewer"), "a, b, c, d"));
    await session.step(31, "When user picks \"Properties...\" from the context menu of Dendrogram viewer", () => pickFromContextMenu(page, "Properties...", el("Dendrogram viewer")));
    await run.scenario("The property panel shows the viewer's categories", async () => {
      await session.step(34, "Then Data category should be visible", () => shouldBe(page, el("Data category"), "visible"));
      await session.step(35, "And Style category should be visible", () => shouldBe(page, el("Style category"), "visible"));
      await session.step(36, "And Behavior category should be visible", () => shouldBe(page, el("Behavior category"), "visible"));
      await session.step(37, "And \"Newick Tag\" property should be visible", () => shouldBe(page, el("\"Newick Tag\" property"), "visible"));
      await session.step(38, "When user expands Style category", () => expand(page, el("Style category")));
      await session.step(39, "And user expands Behavior category", () => expand(page, el("Behavior category")));
      await session.step(40, "Then \"Line Width\" property should be visible", () => shouldBe(page, el("\"Line Width\" property"), "visible"));
      await session.step(41, "And \"Show Tooltip\" property should be visible", () => shouldBe(page, el("\"Show Tooltip\" property"), "visible"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Newick wins over Newick Tag; with Newick cleared the named tag is drawn", async () => {
      await session.step(45, "When user selects \".newick-alt\" in \"Newick Tag\" property", () => selectIn(page, ".newick-alt", el("\"Newick Tag\" property")));
      await session.step(46, "Then the \"newick\" reading of Dendrogram viewer should be \"((a:1,b:1):1,(c:1,d:1):1);\"", () => readingReads(page, "newick", el("Dendrogram viewer"), "((a:1,b:1):1,(c:1,d:1):1);"));
      await session.step(47, "When user clears Newick property", () => clearField(page, el("Newick property")));
      await session.step(48, "And user presses Enter in Newick property", () => pressKeyIn(page, "Enter", el("Newick property")));
      await session.step(49, "Then the \"newick\" reading of Dendrogram viewer should be \"(((a:1,b:1):1,c:1):1,d:1);\"", () => readingReads(page, "newick", el("Dendrogram viewer"), "(((a:1,b:1):1,c:1):1,d:1);"));
      await session.step(50, "And the \"leaves\" reading of Dendrogram viewer should be \"a, b, c, d\"", () => readingReads(page, "leaves", el("Dendrogram viewer"), "a, b, c, d"));
      await session.step(51, "When user selects \"\" in \"Newick Tag\" property", () => selectIn(page, "", el("\"Newick Tag\" property")));
      await session.step(52, "Then the \"newick\" reading of Dendrogram viewer should be \"((a:1,b:1):1,(c:1,d:1):1);\"", () => readingReads(page, "newick", el("Dendrogram viewer"), "((a:1,b:1):1,(c:1,d:1):1);"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The node column binds the grid's current row to a node", async () => {
      await session.step(56, "When user clicks on the \"cell 2 of leaf\" area of grid", () => clickArea(page, "cell 2 of leaf", el("grid")));
      await session.step(57, "Then the \"current node\" reading of Dendrogram viewer should be \"b\"", () => readingReads(page, "current node", el("Dendrogram viewer"), "b"));
      await session.step(58, "When user selects \"value\" in Node property", () => selectIn(page, "value", el("Node property")));
      await session.step(59, "And user clicks on the \"cell 3 of leaf\" area of grid", () => clickArea(page, "cell 3 of leaf", el("grid")));
      await session.step(60, "Then the \"current node\" reading of Dendrogram viewer should be \"\"", () => readingReads(page, "current node", el("Dendrogram viewer"), ""));
      await session.step(61, "When user selects \"leaf\" in Node property", () => selectIn(page, "leaf", el("Node property")));
      await session.step(62, "And user clicks on the \"cell 4 of leaf\" area of grid", () => clickArea(page, "cell 4 of leaf", el("grid")));
      await session.step(63, "Then the \"current node\" reading of Dendrogram viewer should be \"d\"", () => readingReads(page, "current node", el("Dendrogram viewer"), "d"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color and Color Aggr Type color the tree by each aggregation", async () => {
      await session.step(67, "When user selects \"value\" in Color property", () => selectIn(page, "value", el("Color property")));
      await session.step(68, "And user selects \"avg\" in \"Color Aggr Type\" property", () => selectIn(page, "avg", el("\"Color Aggr Type\" property")));
      await session.step(69, "Then the \"color coding\" reading of Dendrogram viewer should be \"avg of value\"", () => readingReads(page, "color coding", el("Dendrogram viewer"), "avg of value"));
      await session.step(70, "When user selects \"min\" in \"Color Aggr Type\" property", () => selectIn(page, "min", el("\"Color Aggr Type\" property")));
      await session.step(71, "Then the \"color coding\" reading of Dendrogram viewer should be \"min of value\"", () => readingReads(page, "color coding", el("Dendrogram viewer"), "min of value"));
      await session.step(72, "When user selects \"max\" in \"Color Aggr Type\" property", () => selectIn(page, "max", el("\"Color Aggr Type\" property")));
      await session.step(73, "Then the \"color coding\" reading of Dendrogram viewer should be \"max of value\"", () => readingReads(page, "color coding", el("Dendrogram viewer"), "max of value"));
      await session.step(74, "When user selects \"med\" in \"Color Aggr Type\" property", () => selectIn(page, "med", el("\"Color Aggr Type\" property")));
      await session.step(75, "Then the \"color coding\" reading of Dendrogram viewer should be \"med of value\"", () => readingReads(page, "color coding", el("Dendrogram viewer"), "med of value"));
      await session.step(76, "When user selects \"count\" in \"Color Aggr Type\" property", () => selectIn(page, "count", el("\"Color Aggr Type\" property")));
      await session.step(77, "Then the \"color coding\" reading of Dendrogram viewer should be \"count of value\"", () => readingReads(page, "color coding", el("Dendrogram viewer"), "count of value"));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Line Width, Node Size and Show Grid reach the main style", async () => {
      await session.step(81, "When user enters \"0\" into \"Line Width\" property", () => enterInto(page, "0", el("\"Line Width\" property")));
      await session.step(82, "Then the \"line width\" reading of Dendrogram viewer should be 0", () => readingIs(page, "line width", el("Dendrogram viewer"), 0));
      await session.step(83, "When user enters \"16\" into \"Line Width\" property", () => enterInto(page, "16", el("\"Line Width\" property")));
      await session.step(84, "Then the \"line width\" reading of Dendrogram viewer should be 16", () => readingIs(page, "line width", el("Dendrogram viewer"), 16));
      await session.step(85, "When user enters \"2.5\" into \"Line Width\" property", () => enterInto(page, "2.5", el("\"Line Width\" property")));
      await session.step(86, "Then the \"line width\" reading of Dendrogram viewer should be 2.5", () => readingIs(page, "line width", el("Dendrogram viewer"), 2.5));
      await session.step(87, "When user enters \"0\" into \"Node Size\" property", () => enterInto(page, "0", el("\"Node Size\" property")));
      await session.step(88, "Then the \"node size\" reading of Dendrogram viewer should be 0", () => readingIs(page, "node size", el("Dendrogram viewer"), 0));
      await session.step(89, "When user enters \"16\" into \"Node Size\" property", () => enterInto(page, "16", el("\"Node Size\" property")));
      await session.step(90, "Then the \"node size\" reading of Dendrogram viewer should be 16", () => readingIs(page, "node size", el("Dendrogram viewer"), 16));
      await session.step(91, "When user enters \"4\" into \"Node Size\" property", () => enterInto(page, "4", el("\"Node Size\" property")));
      await session.step(92, "Then the \"node size\" reading of Dendrogram viewer should be 4", () => readingIs(page, "node size", el("Dendrogram viewer"), 4));
      await session.step(93, "When user checks \"Show Grid\" property", () => check(page, el("\"Show Grid\" property")));
      await session.step(94, "Then the \"show grid\" reading of Dendrogram viewer should be \"true\"", () => readingReads(page, "show grid", el("Dendrogram viewer"), "true"));
      await session.step(95, "When user unchecks \"Show Grid\" property", () => uncheck(page, el("\"Show Grid\" property")));
      await session.step(96, "Then the \"show grid\" reading of Dendrogram viewer should be \"false\"", () => readingReads(page, "show grid", el("Dendrogram viewer"), "false"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The style colors reach their styles", async () => {
      await session.step(100, "When user clicks on value of \"Main Color\" property", () => clickOn(page, el("value of \"Main Color\" property")));
      await session.step(101, "And user picks the color \"#d62728\" in the color picker dialog", () => pickColorSwatch(page, "#d62728"));
      await session.step(102, "And user clicks on label of \"Step Zoom\" property", () => clickOn(page, el("label of \"Step Zoom\" property")));
      await session.step(103, "Then the \"main color\" reading of Dendrogram viewer should be \"#d62728\"", () => readingReads(page, "main color", el("Dendrogram viewer"), "#d62728"));
      await session.step(104, "When user clicks on value of \"Light Color\" property", () => clickOn(page, el("value of \"Light Color\" property")));
      await session.step(105, "And user picks the color \"#d62728\" in the color picker dialog", () => pickColorSwatch(page, "#d62728"));
      await session.step(106, "And user clicks on label of \"Step Zoom\" property", () => clickOn(page, el("label of \"Step Zoom\" property")));
      await session.step(107, "Then the \"light color\" reading of Dendrogram viewer should be \"#d62728\"", () => readingReads(page, "light color", el("Dendrogram viewer"), "#d62728"));
      await session.step(108, "When user clicks on value of \"Current Color\" property", () => clickOn(page, el("value of \"Current Color\" property")));
      await session.step(109, "And user picks the color \"#d62728\" in the color picker dialog", () => pickColorSwatch(page, "#d62728"));
      await session.step(110, "And user clicks on label of \"Step Zoom\" property", () => clickOn(page, el("label of \"Step Zoom\" property")));
      await session.step(111, "Then the \"current color\" reading of Dendrogram viewer should be \"#d62728\"", () => readingReads(page, "current color", el("Dendrogram viewer"), "#d62728"));
      await session.step(112, "When user clicks on value of \"Selections Color\" property", () => clickOn(page, el("value of \"Selections Color\" property")));
      await session.step(113, "And user picks the color \"#d62728\" in the color picker dialog", () => pickColorSwatch(page, "#d62728"));
      await session.step(114, "And user clicks on label of \"Step Zoom\" property", () => clickOn(page, el("label of \"Step Zoom\" property")));
      await session.step(115, "Then the \"selections color\" reading of Dendrogram viewer should be \"#d62728\"", () => readingReads(page, "selections color", el("Dendrogram viewer"), "#d62728"));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Mouse Over Color reaches its style", async () => {
      await session.step(120, "When user clicks on value of \"Mouse Over Color\" property", () => clickOn(page, el("value of \"Mouse Over Color\" property")));
      await session.step(121, "And user picks the color \"#d62728\" in the color picker dialog", () => pickColorSwatch(page, "#d62728"));
      await session.step(122, "And user clicks on label of \"Step Zoom\" property", () => clickOn(page, el("label of \"Step Zoom\" property")));
      await session.step(123, "Then the \"mouse over color\" reading of Dendrogram viewer should be \"#d62728\"", () => readingReads(page, "mouse over color", el("Dendrogram viewer"), "#d62728"));
    }, {knownFailure: true});
    await run.scenario("Show Labels makes the tree draw its labels", async () => {
      await session.step(127, "When user checks \"Show Labels\" property", () => check(page, el("\"Show Labels\" property")));
      await session.step(128, "Then the \"labels drawn\" reading of Dendrogram viewer should be \"true\"", () => readingReads(page, "labels drawn", el("Dendrogram viewer"), "true"));
    }, {knownFailure: true});
    await run.scenario("Font reaches the labels", async () => {
      await session.step(132, "When user enters \"12pt monospace\" into Font property", () => enterInto(page, "12pt monospace", el("Font property")));
      await session.step(133, "Then the \"label font\" reading of Dendrogram viewer should be \"12pt monospace\"", () => readingReads(page, "label font", el("Dendrogram viewer"), "12pt monospace"));
    }, {knownFailure: true});
    await run.scenario("Step sets the leaf row spacing", async () => {
      await session.step(137, "When user enters \"40\" into Step property", () => enterInto(page, "40", el("Step property")));
      await session.step(138, "Then the \"row step\" reading of Dendrogram viewer should be \"40\"", () => readingReads(page, "row step", el("Dendrogram viewer"), "40"));
    }, {knownFailure: true});
    await run.scenario("Step Zoom sets the zoom step", async () => {
      await session.step(142, "When user enters \"2\" into \"Step Zoom\" property", () => enterInto(page, "2", el("\"Step Zoom\" property")));
      await session.step(143, "Then the \"zoom step\" reading of Dendrogram viewer should be \"2\"", () => readingReads(page, "zoom step", el("Dendrogram viewer"), "2"));
    }, {knownFailure: true});
    run.finish();
  });
});
