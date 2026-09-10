/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/word-cloud/word-cloud.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.word-cloud]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, clearSelection, filterPasses, noneSelected, onlyOfSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaAtLeastTall, clickArea, hasArea, hasNoArea, hoverArea, noErrors, oneTooltip, painted, pointerAway, readingIs, readingReads, repainted, reportsError, reportsNoError, setProperties, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaBiggerThanArea, noSuchReading, readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Word cloud counts, the category gate and the viewer filter", () => {
  const session = feature(test, "features/viewers/word-cloud/word-cloud.feature", import.meta.url);
  test("Word cloud counts, the category gate and the viewer filter", {tag: ["@journey", "@viewers", "@realizes:viewers.word-cloud"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "And user adds a word cloud viewer", () => addViewer(page, "word cloud"));
    await session.step(23, "Then word cloud viewer should be visible", () => shouldBe(page, el("word cloud viewer"), "visible"));
    await session.step(24, "And 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(25, "And the \"column\" reading of word cloud viewer should be \"SEX\"", () => readingReads(page, "column", el("word cloud viewer"), "SEX"));
    await session.step(26, "And word cloud viewer should be painted", () => painted(page, el("word cloud viewer")));
    await run.scenario("The cloud draws one word per distinct value with that value's row count", async () => {
      await session.step(29, "Then the \"words\" reading of word cloud viewer should be 2", () => readingIs(page, "words", el("word cloud viewer"), 2));
      await session.step(30, "And the \"word names\" reading of word cloud viewer should be \"F, M\"", () => readingReads(page, "word names", el("word cloud viewer"), "F, M"));
      await session.step(31, "And the \"rows of word \\\"F\\\"\" reading of word cloud viewer should be 553", () => readingIs(page, "rows of word \"F\"", el("word cloud viewer"), 553));
      await session.step(32, "And the \"rows of word \\\"M\\\"\" reading of word cloud viewer should be 447", () => readingIs(page, "rows of word \"M\"", el("word cloud viewer"), 447));
      await session.step(33, "And the \"rows shown\" reading of word cloud viewer should be 1000", () => readingIs(page, "rows shown", el("word cloud viewer"), 1000));
      await session.step(34, "And word cloud viewer should have a \"word \\\"F\\\"\" area", () => hasArea(page, el("word cloud viewer"), "word \"F\""));
      await session.step(35, "And word cloud viewer should have a \"word \\\"M\\\"\" area", () => hasArea(page, el("word cloud viewer"), "word \"M\""));
      await session.step(36, "And word cloud viewer should report no error", () => reportsNoError(page, el("word cloud viewer")));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Column RACE draws one word per race with its row count", async () => {
      await session.step(40, "When user sets \"wordColumnName\" property of word cloud viewer to \"RACE\"", () => setProperty(page, "wordColumnName", el("word cloud viewer"), "RACE"));
      await session.step(41, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(42, "And the \"word names\" reading of word cloud viewer should contain \"Caucasian\"", () => readingContains(page, "word names", el("word cloud viewer"), "Caucasian"));
      await session.step(43, "And the \"word names\" reading of word cloud viewer should contain \"Asian\"", () => readingContains(page, "word names", el("word cloud viewer"), "Asian"));
      await session.step(44, "And the \"rows of word \\\"Caucasian\\\"\" reading of word cloud viewer should be 896", () => readingIs(page, "rows of word \"Caucasian\"", el("word cloud viewer"), 896));
      await session.step(45, "And the \"rows of word \\\"Other\\\"\" reading of word cloud viewer should be 62", () => readingIs(page, "rows of word \"Other\"", el("word cloud viewer"), 62));
      await session.step(46, "And the \"rows of word \\\"Black\\\"\" reading of word cloud viewer should be 27", () => readingIs(page, "rows of word \"Black\"", el("word cloud viewer"), 27));
      await session.step(47, "And the \"rows of word \\\"Asian\\\"\" reading of word cloud viewer should be 15", () => readingIs(page, "rows of word \"Asian\"", el("word cloud viewer"), 15));
      await session.step(48, "And word cloud viewer should have repainted", () => repainted(page, el("word cloud viewer")));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The word of the bigger count gets the bigger box", async () => {
      await session.step(55, "Given user sets properties of word cloud viewer:", () => setProperties(page, el("word cloud viewer"), [["minRotationDegree","0"],["maxRotationDegree","0"]]));
      await session.step(58, "Then the \"word \\\"Caucasian\\\"\" area of word cloud viewer should be taller than the \"word \\\"Other\\\"\" area", () => areaBiggerThanArea(page, "word \"Caucasian\"", el("word cloud viewer"), "taller", "word \"Other\""));
      await session.step(59, "And the \"word \\\"Other\\\"\" area of word cloud viewer should be taller than the \"word \\\"Asian\\\"\" area", () => areaBiggerThanArea(page, "word \"Other\"", el("word cloud viewer"), "taller", "word \"Asian\""));
      await session.step(60, "And the \"word \\\"Caucasian\\\"\" area of word cloud viewer should be at least 60 pixels tall", () => areaAtLeastTall(page, "word \"Caucasian\"", el("word cloud viewer"), 60));
      await session.step(61, "When user sets properties of word cloud viewer:", () => setProperties(page, el("word cloud viewer"), [["minRotationDegree","-30"],["maxRotationDegree","30"]]));
      await session.step(64, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column with more than 500 categories reports the gate and shows no cloud at all", async () => {
      await session.step(68, "When user sets \"wordColumnName\" property of word cloud viewer to \"USUBJID\"", () => setProperty(page, "wordColumnName", el("word cloud viewer"), "USUBJID"));
      await session.step(69, "Then word cloud viewer should report the error \"The Word cloud viewer requires categorical column with 500 or fewer unique categories\"", () => reportsError(page, el("word cloud viewer"), "The Word cloud viewer requires categorical column with 500 or fewer unique categories"));
      await session.step(70, "And word cloud viewer should not have a \"word \\\"Caucasian\\\"\" area", () => hasNoArea(page, el("word cloud viewer"), "word \"Caucasian\""));
      await session.step(71, "And word cloud viewer should not have a \"view\" area", () => hasNoArea(page, el("word cloud viewer"), "view"));
      await session.step(72, "And word cloud viewer should not report a \"words\" reading", () => noSuchReading(page, el("word cloud viewer"), "words"));
      await session.step(73, "And word cloud viewer should not report a \"word names\" reading", () => noSuchReading(page, el("word cloud viewer"), "word names"));
      await session.step(74, "And word cloud viewer should not report a \"rows of word \\\"Caucasian\\\"\" reading", () => noSuchReading(page, el("word cloud viewer"), "rows of word \"Caucasian\""));
      await session.step(75, "And the \"column\" reading of word cloud viewer should be \"USUBJID\"", () => readingReads(page, "column", el("word cloud viewer"), "USUBJID"));
      await session.step(76, "And the \"rows shown\" reading of word cloud viewer should be 1000", () => readingIs(page, "rows shown", el("word cloud viewer"), 1000));
      await session.step(77, "When user sets \"wordColumnName\" property of word cloud viewer to \"RACE\"", () => setProperty(page, "wordColumnName", el("word cloud viewer"), "RACE"));
      await session.step(78, "Then word cloud viewer should report no error", () => reportsNoError(page, el("word cloud viewer")));
      await session.step(79, "And the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(80, "And word cloud viewer should have a \"word \\\"Caucasian\\\"\" area", () => hasArea(page, el("word cloud viewer"), "word \"Caucasian\""));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hovering a word shows how many rows carry it", async () => {
      await session.step(87, "Given user sets properties of word cloud viewer:", () => setProperties(page, el("word cloud viewer"), [["minRotationDegree","0"],["maxRotationDegree","0"]]));
      await session.step(90, "When user hovers over the \"word \\\"Caucasian\\\"\" area of word cloud viewer", () => hoverArea(page, "word \"Caucasian\"", el("word cloud viewer")));
      await session.step(91, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(92, "And tooltip should contain text \"896 rows\"", () => shouldContainText(page, el("tooltip"), "896 rows"));
      await session.step(93, "When user moves the pointer away from word cloud viewer", () => pointerAway(page, el("word cloud viewer")));
      await session.step(94, "And user hovers over the \"word \\\"Other\\\"\" area of word cloud viewer", () => hoverArea(page, "word \"Other\"", el("word cloud viewer")));
      await session.step(95, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(96, "And tooltip should contain text \"62 rows\"", () => shouldContainText(page, el("tooltip"), "62 rows"));
      await session.step(97, "When user moves the pointer away from word cloud viewer", () => pointerAway(page, el("word cloud viewer")));
      await session.step(98, "And user sets properties of word cloud viewer:", () => setProperties(page, el("word cloud viewer"), [["minRotationDegree","-30"],["maxRotationDegree","30"]]));
      await session.step(101, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(102, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clicking a word selects exactly that word's rows", async () => {
      await session.step(105, "Given user clears the row selection", () => clearSelection(page));
      await session.step(106, "When user clicks on the \"word \\\"Other\\\"\" area of word cloud viewer", () => clickArea(page, "word \"Other\"", el("word cloud viewer")));
      await session.step(107, "Then 62 rows should be selected", () => selectedRowCount(page, 62));
      await session.step(108, "And only rows where \"RACE\" is \"Other\" should be selected", () => onlyOfSelected(page, "RACE", "Other"));
      await session.step(109, "When user clicks on the \"word \\\"Asian\\\"\" area of word cloud viewer", () => clickArea(page, "word \"Asian\"", el("word cloud viewer")));
      await session.step(110, "Then 15 rows should be selected", () => selectedRowCount(page, 15));
      await session.step(111, "And only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(112, "When user clears the row selection", () => clearSelection(page));
      await session.step(113, "Then no rows should be selected", () => noneSelected(page));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A table filter moves the counts and leaves the names alone", async () => {
      await session.step(117, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(118, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(119, "And the \"rows shown\" reading of word cloud viewer should be 447", () => readingIs(page, "rows shown", el("word cloud viewer"), 447));
      await session.step(120, "And the \"rows of word \\\"Caucasian\\\"\" reading of word cloud viewer should be 416", () => readingIs(page, "rows of word \"Caucasian\"", el("word cloud viewer"), 416));
      await session.step(121, "And the \"rows of word \\\"Other\\\"\" reading of word cloud viewer should be 14", () => readingIs(page, "rows of word \"Other\"", el("word cloud viewer"), 14));
      await session.step(122, "And the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
      await session.step(123, "And the \"word names\" reading of word cloud viewer should contain \"Caucasian\"", () => readingContains(page, "word names", el("word cloud viewer"), "Caucasian"));
      await session.step(124, "And word cloud viewer should report no error", () => reportsNoError(page, el("word cloud viewer")));
      await session.step(125, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(126, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(127, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(128, "And the \"rows of word \\\"Caucasian\\\"\" reading of word cloud viewer should be 896", () => readingIs(page, "rows of word \"Caucasian\"", el("word cloud viewer"), 896));
      await session.step(129, "And the \"rows shown\" reading of word cloud viewer should be 1000", () => readingIs(page, "rows shown", el("word cloud viewer"), 1000));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title bar closes the cloud", async () => {
      await session.step(133, "When user clicks on close icon of word cloud viewer", () => clickOn(page, el("close icon of word cloud viewer")));
      await session.step(134, "Then word cloud viewer should be absent", () => shouldBe(page, el("word cloud viewer"), "absent"));
      await session.step(135, "And the open tableview should have 0 word cloud viewers", () => viewerCount(page, 0, "word cloud"));
      await session.step(136, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
