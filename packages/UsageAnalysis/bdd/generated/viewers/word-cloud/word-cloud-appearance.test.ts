/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/word-cloud/word-cloud-appearance.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.word-cloud, public#4424]
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
import {clickOn, enterInto, isExpanded, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaNarrower, areaShorter, areaTaller, areasSameSize, closeContextMenu, hasArea, menuLists, noErrors, openContextMenu, painted, pickFromContextMenu, propertyShouldBe, readingHigherThanRemembered, readingIs, readingLowerThanRemembered, readingReads, rememberReading, repainted, reportsNoError, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaBiggerThanArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Word cloud text size, rotation, shape and font", () => {
  const session = feature(test, "features/viewers/word-cloud/word-cloud-appearance.feature", import.meta.url);
  test("Word cloud text size, rotation, shape and font", {tag: ["@journey", "@viewers", "@realizes:viewers.word-cloud", "@realizes:public#4424"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(28, "And user adds a word cloud viewer with:", () => addViewerWith(page, "word cloud", [["wordColumnName","RACE"],["minRotationDegree","0"],["maxRotationDegree","0"]]), [["wordColumnName","RACE"],["minRotationDegree","0"],["maxRotationDegree","0"]]);
    await session.step(32, "Then word cloud viewer should be visible", () => shouldBe(page, el("word cloud viewer"), "visible"));
    await session.step(33, "And the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
    await session.step(34, "And the \"rows of word \\\"Caucasian\\\"\" reading of word cloud viewer should be 896", () => readingIs(page, "rows of word \"Caucasian\"", el("word cloud viewer"), 896));
    await session.step(35, "And word cloud viewer should be painted", () => painted(page, el("word cloud viewer")));
    await run.scenario("Text size spreads the boxes by count, and equal min and max flattens them", async () => {
      await session.step(38, "Then the \"word \\\"Caucasian\\\"\" area of word cloud viewer should be taller than the \"word \\\"Other\\\"\" area", () => areaBiggerThanArea(page, "word \"Caucasian\"", el("word cloud viewer"), "taller", "word \"Other\""));
      await session.step(39, "And the \"word \\\"Other\\\"\" area of word cloud viewer should be taller than the \"word \\\"Black\\\"\" area", () => areaBiggerThanArea(page, "word \"Other\"", el("word cloud viewer"), "taller", "word \"Black\""));
      await session.step(40, "When user sets properties of word cloud viewer:", () => setProperties(page, el("word cloud viewer"), [["minTextSize","20"],["maxTextSize","20"]]), [["minTextSize","20"],["maxTextSize","20"]]);
      await session.step(43, "Then the \"word \\\"Caucasian\\\"\" and \"word \\\"Black\\\"\" areas of word cloud viewer should be the same height", () => areasSameSize(page, "word \"Caucasian\"", "word \"Black\"", el("word cloud viewer"), "height"));
      await session.step(44, "And the \"word \\\"Caucasian\\\"\" and \"word \\\"Other\\\"\" areas of word cloud viewer should be the same height", () => areasSameSize(page, "word \"Caucasian\"", "word \"Other\"", el("word cloud viewer"), "height"));
      await session.step(45, "And the \"word \\\"Caucasian\\\"\" and \"word \\\"Asian\\\"\" areas of word cloud viewer should be the same height", () => areasSameSize(page, "word \"Caucasian\"", "word \"Asian\"", el("word cloud viewer"), "height"));
      await session.step(46, "And the \"word \\\"Caucasian\\\"\" area of word cloud viewer should be shorter than before", () => areaShorter(page, "word \"Caucasian\"", el("word cloud viewer")));
      await session.step(47, "And the \"word \\\"Caucasian\\\"\" area of word cloud viewer should be wider than the \"word \\\"Black\\\"\" area", () => areaBiggerThanArea(page, "word \"Caucasian\"", el("word cloud viewer"), "wider", "word \"Black\""));
      await session.step(48, "When user sets properties of word cloud viewer:", () => setProperties(page, el("word cloud viewer"), [["minTextSize","14"],["maxTextSize","100"]]), [["minTextSize","14"],["maxTextSize","100"]]);
      await session.step(51, "Then the \"word \\\"Caucasian\\\"\" area of word cloud viewer should be taller than before", () => areaTaller(page, "word \"Caucasian\"", el("word cloud viewer")));
      await session.step(52, "And the \"word \\\"Caucasian\\\"\" area of word cloud viewer should be taller than the \"word \\\"Black\\\"\" area", () => areaBiggerThanArea(page, "word \"Caucasian\"", el("word cloud viewer"), "taller", "word \"Black\""));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A rotation of 90 degrees turns every box on its side", async () => {
      await session.step(57, "When user sets properties of word cloud viewer:", () => setProperties(page, el("word cloud viewer"), [["minTextSize","20"],["maxTextSize","20"]]), [["minTextSize","20"],["maxTextSize","20"]]);
      await session.step(60, "Then the \"word \\\"Caucasian\\\"\" area of word cloud viewer should be wider than the \"word \\\"Black\\\"\" area", () => areaBiggerThanArea(page, "word \"Caucasian\"", el("word cloud viewer"), "wider", "word \"Black\""));
      await session.step(61, "When user sets properties of word cloud viewer:", () => setProperties(page, el("word cloud viewer"), [["minRotationDegree","90"],["maxRotationDegree","90"]]), [["minRotationDegree","90"],["maxRotationDegree","90"]]);
      await session.step(64, "Then the \"word \\\"Caucasian\\\"\" area of word cloud viewer should be narrower than before", () => areaNarrower(page, "word \"Caucasian\"", el("word cloud viewer")));
      await session.step(65, "And the \"word \\\"Caucasian\\\"\" area of word cloud viewer should be taller than before", () => areaTaller(page, "word \"Caucasian\"", el("word cloud viewer")));
      await session.step(66, "And the \"word \\\"Caucasian\\\"\" and \"word \\\"Black\\\"\" areas of word cloud viewer should be the same width", () => areasSameSize(page, "word \"Caucasian\"", "word \"Black\"", el("word cloud viewer"), "width"));
      await session.step(67, "And the \"word \\\"Caucasian\\\"\" area of word cloud viewer should be taller than the \"word \\\"Black\\\"\" area", () => areaBiggerThanArea(page, "word \"Caucasian\"", el("word cloud viewer"), "taller", "word \"Black\""));
      await session.step(68, "When user sets properties of word cloud viewer:", () => setProperties(page, el("word cloud viewer"), [["minRotationDegree","0"],["maxRotationDegree","0"],["minTextSize","14"],["maxTextSize","100"]]), [["minRotationDegree","0"],["maxRotationDegree","0"],["minTextSize","14"],["maxTextSize","100"]]);
      await session.step(73, "Then the \"word \\\"Caucasian\\\"\" area of word cloud viewer should be wider than the \"word \\\"Black\\\"\" area", () => areaBiggerThanArea(page, "word \"Caucasian\"", el("word cloud viewer"), "wider", "word \"Black\""));
      await session.step(74, "When user clicks on settings icon of word cloud viewer", () => clickOn(page, el("settings icon of word cloud viewer")));
      await session.step(75, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(76, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(77, "When user enters \"0\" into \"Rotation Step\" property in context panel", () => enterInto(page, "0", el("\"Rotation Step\" property in context panel")));
      await session.step(78, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(79, "When user enters \"5\" into \"Rotation Step\" property in context panel", () => enterInto(page, "5", el("\"Rotation Step\" property in context panel")));
      await session.step(80, "Then \"Rotation Step\" property of word cloud viewer should be \"5\"", () => propertyShouldBe(page, "Rotation Step", el("word cloud viewer"), "5"));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
      await session.step(82, "When user presses F4", () => pressKey(page, "F4"));
      await session.step(83, "Then context panel should be hidden", () => shouldBe(page, el("context panel"), "hidden"));
    });
    await run.scenario("Every shape keeps all four words and their counts", async () => {
      await session.step(86, "Then \"shape\" property of word cloud viewer should be \"circle\"", () => propertyShouldBe(page, "shape", el("word cloud viewer"), "circle"));
      await session.step(87, "When user sets \"shape\" property of word cloud viewer to \"diamond\"", () => setProperty(page, "shape", el("word cloud viewer"), "diamond"));
      await session.step(88, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(89, "And the \"rows of word \\\"Caucasian\\\"\" reading of word cloud viewer should be 896", () => readingIs(page, "rows of word \"Caucasian\"", el("word cloud viewer"), 896));
      await session.step(90, "And word cloud viewer should have a \"word \\\"Asian\\\"\" area", () => hasArea(page, el("word cloud viewer"), "word \"Asian\""));
      await session.step(91, "And word cloud viewer should have repainted", () => repainted(page, el("word cloud viewer")));
      await session.step(92, "When user sets \"shape\" property of word cloud viewer to \"star\"", () => setProperty(page, "shape", el("word cloud viewer"), "star"));
      await session.step(93, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(94, "And the \"rows of word \\\"Asian\\\"\" reading of word cloud viewer should be 15", () => readingIs(page, "rows of word \"Asian\"", el("word cloud viewer"), 15));
      await session.step(95, "And word cloud viewer should have a \"word \\\"Black\\\"\" area", () => hasArea(page, el("word cloud viewer"), "word \"Black\""));
      await session.step(96, "When user sets \"shape\" property of word cloud viewer to \"pentagon\"", () => setProperty(page, "shape", el("word cloud viewer"), "pentagon"));
      await session.step(97, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(98, "And word cloud viewer should have a \"word \\\"Other\\\"\" area", () => hasArea(page, el("word cloud viewer"), "word \"Other\""));
      await session.step(99, "When user sets \"shape\" property of word cloud viewer to \"triangle\"", () => setProperty(page, "shape", el("word cloud viewer"), "triangle"));
      await session.step(100, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(101, "And word cloud viewer should have a \"word \\\"Caucasian\\\"\" area", () => hasArea(page, el("word cloud viewer"), "word \"Caucasian\""));
      await session.step(102, "When user sets \"shape\" property of word cloud viewer to \"circle\"", () => setProperty(page, "shape", el("word cloud viewer"), "circle"));
      await session.step(103, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(104, "And word cloud viewer should report no error", () => reportsNoError(page, el("word cloud viewer")));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The font reading is the font the words were drawn with, not the property", async () => {
      await session.step(108, "Then the \"font\" reading of word cloud viewer should be \"bold sans-serif\"", () => readingReads(page, "font", el("word cloud viewer"), "bold sans-serif"));
      await session.step(109, "When user sets \"bold\" property of word cloud viewer to \"false\"", () => setProperty(page, "bold", el("word cloud viewer"), "false"));
      await session.step(110, "Then the \"font\" reading of word cloud viewer should be \"normal sans-serif\"", () => readingReads(page, "font", el("word cloud viewer"), "normal sans-serif"));
      await session.step(111, "When user sets \"fontFamily\" property of word cloud viewer to \"monospace\"", () => setProperty(page, "fontFamily", el("word cloud viewer"), "monospace"));
      await session.step(112, "Then the \"font\" reading of word cloud viewer should be \"normal monospace\"", () => readingReads(page, "font", el("word cloud viewer"), "normal monospace"));
      await session.step(113, "And the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(114, "When user sets properties of word cloud viewer:", () => setProperties(page, el("word cloud viewer"), [["bold","true"],["fontFamily","sans-serif"]]), [["bold","true"],["fontFamily","sans-serif"]]);
      await session.step(117, "Then the \"font\" reading of word cloud viewer should be \"bold sans-serif\"", () => readingReads(page, "font", el("word cloud viewer"), "bold sans-serif"));
      await session.step(118, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Grid Size spaces the words, so a coarser grid fits fewer of them", async () => {
      await session.step(121, "When user sets properties of word cloud viewer:", () => setProperties(page, el("word cloud viewer"), [["wordColumnName","DEMOG"],["drawOutOfBound","false"]]), [["wordColumnName","DEMOG"],["drawOutOfBound","false"]]);
      await session.step(124, "Then \"gridSize\" property of word cloud viewer should be \"8\"", () => propertyShouldBe(page, "gridSize", el("word cloud viewer"), "8"));
      await session.step(125, "When user remembers the \"words\" reading of word cloud viewer", () => rememberReading(page, "words", el("word cloud viewer")));
      await session.step(126, "And user sets \"gridSize\" property of word cloud viewer to \"2\"", () => setProperty(page, "gridSize", el("word cloud viewer"), "2"));
      await session.step(127, "Then the \"words\" reading of word cloud viewer should be higher than remembered", () => readingHigherThanRemembered(page, "words", el("word cloud viewer")));
      await session.step(128, "When user sets \"gridSize\" property of word cloud viewer to \"40\"", () => setProperty(page, "gridSize", el("word cloud viewer"), "40"));
      await session.step(129, "Then the \"words\" reading of word cloud viewer should be lower than remembered", () => readingLowerThanRemembered(page, "words", el("word cloud viewer")));
      await session.step(130, "When user sets properties of word cloud viewer:", () => setProperties(page, el("word cloud viewer"), [["gridSize","8"],["drawOutOfBound","true"],["wordColumnName","RACE"]]), [["gridSize","8"],["drawOutOfBound","true"],["wordColumnName","RACE"]]);
      await session.step(134, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(135, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The context menu offers the viewer's commands and opens its properties", async () => {
      await session.step(138, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(139, "When user opens the context menu of word cloud viewer", () => openContextMenu(page, el("word cloud viewer")));
      await session.step(140, "Then the open menu should list \"Properties...\"", () => menuLists(page, "Properties..."));
      await session.step(141, "And the open menu should list \"General\"", () => menuLists(page, "General"));
      await session.step(142, "And the open menu should list \"Pick Up / Apply\"", () => menuLists(page, "Pick Up / Apply"));
      await session.step(143, "And the open menu should list \"To Script\"", () => menuLists(page, "To Script"));
      await session.step(144, "And the open menu should list \"General > Save as PNG\"", () => menuLists(page, "General > Save as PNG"));
      await session.step(145, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(146, "And user picks \"Properties...\" from the context menu of word cloud viewer", () => pickFromContextMenu(page, "Properties...", el("word cloud viewer")));
      await session.step(147, "Then the context panel should show \"Word cloud\"", () => contextPanelShows(page, "Word cloud"));
      await session.step(148, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
