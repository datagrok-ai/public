/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/viewer-rules-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [public#4423, public#4425, core#1174, core#1200, public#4424]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterBetween, filterPasses, filterPassesAll, filterPassesFewer, filterTo, removeColumn, renameColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, clickArea, loadLayout, noErrors, propertyShouldBe, propertyShouldNotBe, readingIs, reportsNoError, saveLayout, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Charts viewer edge cases shared by every viewer: no rows, a column removed or renamed, a layout with a missing column, a filter left by a closed viewer, values out of bounds, empty categories, a category limit", () => {
  const session = feature(test, "features/viewers/viewer-rules-edge-cases.feature", import.meta.url);
  test("With every row filtered out a sunburst shows no row and no message, and the rows come back [viewer=sunburst, property=hierarchyColumnNames, value=SEX, RACE]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(30, "When user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(32, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(33, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(34, "And the \"rows shown\" reading of sunburst viewer should be 0", () => readingIs(page, "rows shown", el("sunburst viewer"), 0));
    await session.step(35, "And sunburst viewer should report no error", () => reportsNoError(page, el("sunburst viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of sunburst viewer should be 1000", () => readingIs(page, "rows shown", el("sunburst viewer"), 1000));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a tree shows no row and no message, and the rows come back [viewer=tree, property=hierarchyColumnNames, value=CONTROL, SEX, RACE]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(30, "When user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(32, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(33, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(34, "And the \"rows shown\" reading of tree viewer should be 0", () => readingIs(page, "rows shown", el("tree viewer"), 0));
    await session.step(35, "And tree viewer should report no error", () => reportsNoError(page, el("tree viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of tree viewer should be 1000", () => readingIs(page, "rows shown", el("tree viewer"), 1000));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a word cloud shows no row and no message, and the rows come back [viewer=word cloud, property=wordColumnName, value=RACE]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(30, "When user adds a word cloud viewer with:", () => addViewerWith(page, "word cloud", [["wordColumnName","RACE"]]), [["wordColumnName","RACE"]]);
    await session.step(32, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(33, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(34, "And the \"rows shown\" reading of word cloud viewer should be 0", () => readingIs(page, "rows shown", el("word cloud viewer"), 0));
    await session.step(35, "And word cloud viewer should report no error", () => reportsNoError(page, el("word cloud viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of word cloud viewer should be 1000", () => readingIs(page, "rows shown", el("word cloud viewer"), 1000));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("The sizeColumnName column removed from the table leaves the settings of the tree [viewer=tree, property=sizeColumnName, column=HEIGHT]", {tag: ["@viewers", "@edge-cases", "@realizes:public#4423"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(47, "When user adds a tree viewer with:", () => addViewerWith(page, "tree", [["sizeColumnName","HEIGHT"]]), [["sizeColumnName","HEIGHT"]]);
    await session.step(49, "And user removes \"HEIGHT\" column", () => removeColumn(page, "HEIGHT"));
    await session.step(50, "Then no errors should have been logged", () => noErrors(page));
    await session.step(51, "And \"sizeColumnName\" property of tree viewer should not be \"HEIGHT\"", () => propertyShouldNotBe(page, "sizeColumnName", el("tree viewer"), "HEIGHT"));
    await session.step(52, "And tree viewer should report no error", () => reportsNoError(page, el("tree viewer")));
  });
  test("The colorColumnName column removed from the table leaves the settings of the tree [viewer=tree, property=colorColumnName, column=AGE]", {tag: ["@viewers", "@edge-cases", "@realizes:public#4423"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(47, "When user adds a tree viewer with:", () => addViewerWith(page, "tree", [["colorColumnName","AGE"]]), [["colorColumnName","AGE"]]);
    await session.step(49, "And user removes \"AGE\" column", () => removeColumn(page, "AGE"));
    await session.step(50, "Then no errors should have been logged", () => noErrors(page));
    await session.step(51, "And \"colorColumnName\" property of tree viewer should not be \"AGE\"", () => propertyShouldNotBe(page, "colorColumnName", el("tree viewer"), "AGE"));
    await session.step(52, "And tree viewer should report no error", () => reportsNoError(page, el("tree viewer")));
  });
  test("The RACE column renamed is followed by the hierarchyColumnNames setting of the tree [viewer=tree, property=hierarchyColumnNames, value=CONTROL, SEX, RACE, column=RACE, new name=ORIGIN, renamed=CONTROL, SEX, ORIGIN]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(61, "When user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(63, "And user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(64, "Then \"hierarchyColumnNames\" property of tree viewer should be \"CONTROL, SEX, ORIGIN\"", () => propertyShouldBe(page, "hierarchyColumnNames", el("tree viewer"), "CONTROL, SEX, ORIGIN"));
    await session.step(65, "And tree viewer should report no error", () => reportsNoError(page, el("tree viewer")));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
  test("The RACE column renamed is followed by the wordColumnName setting of the word cloud [viewer=word cloud, property=wordColumnName, value=RACE, column=RACE, new name=ORIGIN, renamed=ORIGIN]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(61, "When user adds a word cloud viewer with:", () => addViewerWith(page, "word cloud", [["wordColumnName","RACE"]]), [["wordColumnName","RACE"]]);
    await session.step(63, "And user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(64, "Then \"wordColumnName\" property of word cloud viewer should be \"ORIGIN\"", () => propertyShouldBe(page, "wordColumnName", el("word cloud viewer"), "ORIGIN"));
    await session.step(65, "And word cloud viewer should report no error", () => reportsNoError(page, el("word cloud viewer")));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
  test("The RACE column renamed is followed by the hierarchyColumnNames setting of the sunburst [viewer=sunburst, property=hierarchyColumnNames, value=SEX, RACE, column=RACE, new name=ORIGIN, renamed=SEX, ORIGIN]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(61, "When user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(63, "And user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(64, "Then \"hierarchyColumnNames\" property of sunburst viewer should be \"SEX, ORIGIN\"", () => propertyShouldBe(page, "hierarchyColumnNames", el("sunburst viewer"), "SEX, ORIGIN"));
    await session.step(65, "And sunburst viewer should report no error", () => reportsNoError(page, el("sunburst viewer")));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout naming a SEX column the table no longer has is applied to the sunburst, which draws without it [viewer=sunburst, property=hierarchyColumnNames, value=SEX, RACE, column=SEX, left=SEX, RACE]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(75, "When user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(77, "And user saves the layout of the current table view", () => saveLayout(page));
    await session.step(78, "And user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(79, "And user loads the saved layout", () => loadLayout(page));
    await session.step(80, "Then sunburst viewer should be visible", () => shouldBe(page, el("sunburst viewer"), "visible"));
    await session.step(81, "And \"hierarchyColumnNames\" property of sunburst viewer should be \"SEX, RACE\"", () => propertyShouldBe(page, "hierarchyColumnNames", el("sunburst viewer"), "SEX, RACE"));
    await session.step(82, "And sunburst viewer should report no error", () => reportsNoError(page, el("sunburst viewer")));
    await session.step(83, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout naming a RACE column the table no longer has is applied to the tree, which draws without it [viewer=tree, property=hierarchyColumnNames, value=CONTROL, SEX, RACE, column=RACE, left=CONTROL, SEX, RACE]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(75, "When user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(77, "And user saves the layout of the current table view", () => saveLayout(page));
    await session.step(78, "And user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(79, "And user loads the saved layout", () => loadLayout(page));
    await session.step(80, "Then tree viewer should be visible", () => shouldBe(page, el("tree viewer"), "visible"));
    await session.step(81, "And \"hierarchyColumnNames\" property of tree viewer should be \"CONTROL, SEX, RACE\"", () => propertyShouldBe(page, "hierarchyColumnNames", el("tree viewer"), "CONTROL, SEX, RACE"));
    await session.step(82, "And tree viewer should report no error", () => reportsNoError(page, el("tree viewer")));
    await session.step(83, "And no errors should have been logged", () => noErrors(page));
  });
  test("Closing a sunburst gives back the rows its click filtered out [viewer=sunburst, property=hierarchyColumnNames, value=SEX, RACE, area=segment F]", {tag: ["@viewers", "@edge-cases", "@realizes:public#4425"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(91, "When user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
    await session.step(93, "And user sets \"onClick\" property of sunburst viewer to \"Filter\"", () => setProperty(page, "onClick", el("sunburst viewer"), "Filter"));
    await session.step(94, "And user clicks on the \"segment F\" area of sunburst viewer", () => clickArea(page, "segment F", el("sunburst viewer")));
    await session.step(95, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
    await session.step(96, "When user clicks on close icon of sunburst viewer", () => clickOn(page, el("close icon of sunburst viewer")));
    await session.step(97, "Then the open tableview should have 0 sunburst viewers", () => viewerCount(page, 0, "sunburst"));
    await session.step(98, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(99, "And no errors should have been logged", () => noErrors(page));
  });
  test("Closing a tree gives back the rows its click filtered out [viewer=tree, property=hierarchyColumnNames, value=CONTROL, SEX, RACE, area=branch All | false | F | Asian]", {tag: ["@viewers", "@edge-cases", "@realizes:public#4425"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(91, "When user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
    await session.step(93, "And user sets \"onClick\" property of tree viewer to \"Filter\"", () => setProperty(page, "onClick", el("tree viewer"), "Filter"));
    await session.step(94, "And user clicks on the \"branch All | false | F | Asian\" area of tree viewer", () => clickArea(page, "branch All | false | F | Asian", el("tree viewer")));
    await session.step(95, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
    await session.step(96, "When user clicks on close icon of tree viewer", () => clickOn(page, el("close icon of tree viewer")));
    await session.step(97, "Then the open tableview should have 0 tree viewers", () => viewerCount(page, 0, "tree"));
    await session.step(98, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(99, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the column of the category a tree click filtered by gives the rows back [viewer=tree, property=hierarchyColumnNames, value=CONTROL, SEX, RACE, area=branch All | false | F | Asian, column=RACE]", {tag: ["@viewers", "@edge-cases", "@known-failure:114", "@realizes:core#1174"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(108, "When user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","CONTROL, SEX, RACE"]]), [["hierarchyColumnNames","CONTROL, SEX, RACE"]]);
      await session.step(110, "And user sets \"onClick\" property of tree viewer to \"Filter\"", () => setProperty(page, "onClick", el("tree viewer"), "Filter"));
      await session.step(111, "And user clicks on the \"branch All | false | F | Asian\" area of tree viewer", () => clickArea(page, "branch All | false | F | Asian", el("tree viewer")));
      await session.step(112, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(113, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
      await session.step(114, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
    }, 114);
  });
  test("Removing the column of the category a sunburst click filtered by gives the rows back [viewer=sunburst, property=hierarchyColumnNames, value=SEX, RACE, area=segment F, column=SEX]", {tag: ["@viewers", "@edge-cases", "@known-failure:114", "@realizes:core#1174"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(108, "When user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","SEX, RACE"]]), [["hierarchyColumnNames","SEX, RACE"]]);
      await session.step(110, "And user sets \"onClick\" property of sunburst viewer to \"Filter\"", () => setProperty(page, "onClick", el("sunburst viewer"), "Filter"));
      await session.step(111, "And user clicks on the \"segment F\" area of sunburst viewer", () => clickArea(page, "segment F", el("sunburst viewer")));
      await session.step(112, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(113, "When user removes \"SEX\" column", () => removeColumn(page, "SEX"));
      await session.step(114, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
    }, 114);
  });
  test("Initial Tree Depth typed in the context panel outside its declared bounds is not taken by the tree [viewer=tree, category=Style, caption=Initial Tree Depth, typed=6, kept=3]", {tag: ["@viewers", "@edge-cases", "@known-failure:129", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(124, "When user adds a tree viewer", () => addViewer(page, "tree"));
      await session.step(125, "And user clicks on settings icon of tree viewer", () => clickOn(page, el("settings icon of tree viewer")));
      await session.step(126, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(127, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(128, "When user enters \"6\" into \"Initial Tree Depth\" property in context panel", () => enterInto(page, "6", el("\"Initial Tree Depth\" property in context panel")));
      await session.step(129, "Then \"Initial Tree Depth\" property of tree viewer should be \"3\"", () => propertyShouldBe(page, "Initial Tree Depth", el("tree viewer"), "3"));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
    }, 129);
  });
  test("Font Size typed in the context panel outside its declared bounds is not taken by the tree [viewer=tree, category=Style, caption=Font Size, typed=31, kept=12]", {tag: ["@viewers", "@edge-cases", "@known-failure:129", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(124, "When user adds a tree viewer", () => addViewer(page, "tree"));
      await session.step(125, "And user clicks on settings icon of tree viewer", () => clickOn(page, el("settings icon of tree viewer")));
      await session.step(126, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(127, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(128, "When user enters \"31\" into \"Font Size\" property in context panel", () => enterInto(page, "31", el("\"Font Size\" property in context panel")));
      await session.step(129, "Then \"Font Size\" property of tree viewer should be \"12\"", () => propertyShouldBe(page, "Font Size", el("tree viewer"), "12"));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
    }, 129);
  });
  test("Rotation Step typed in the context panel outside its declared bounds is not taken by the word cloud [viewer=word cloud, category=Misc, caption=Rotation Step, typed=0, kept=5]", {tag: ["@viewers", "@edge-cases", "@known-failure:129", "@realizes:public#4424"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(124, "When user adds a word cloud viewer", () => addViewer(page, "word cloud"));
      await session.step(125, "And user clicks on settings icon of word cloud viewer", () => clickOn(page, el("settings icon of word cloud viewer")));
      await session.step(126, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(127, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(128, "When user enters \"0\" into \"Rotation Step\" property in context panel", () => enterInto(page, "0", el("\"Rotation Step\" property in context panel")));
      await session.step(129, "Then \"Rotation Step\" property of word cloud viewer should be \"5\"", () => propertyShouldBe(page, "Rotation Step", el("word cloud viewer"), "5"));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
    }, 129);
  });
  test("initialTreeDepth set through the API outside its declared bounds is not taken by the tree [viewer=tree, property=initialTreeDepth, set=6, kept=3]", {tag: ["@viewers", "@edge-cases", "@known-failure:146", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(144, "When user adds a tree viewer", () => addViewer(page, "tree"));
      await session.step(145, "And user sets \"initialTreeDepth\" property of tree viewer to \"6\"", () => setProperty(page, "initialTreeDepth", el("tree viewer"), "6"));
      await session.step(146, "Then \"initialTreeDepth\" property of tree viewer should be \"3\"", () => propertyShouldBe(page, "initialTreeDepth", el("tree viewer"), "3"));
      await session.step(147, "And no errors should have been logged", () => noErrors(page));
    }, 146);
  });
  test("fontSize set through the API outside its declared bounds is not taken by the tree [viewer=tree, property=fontSize, set=31, kept=12]", {tag: ["@viewers", "@edge-cases", "@known-failure:146", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(144, "When user adds a tree viewer", () => addViewer(page, "tree"));
      await session.step(145, "And user sets \"fontSize\" property of tree viewer to \"31\"", () => setProperty(page, "fontSize", el("tree viewer"), "31"));
      await session.step(146, "Then \"fontSize\" property of tree viewer should be \"12\"", () => propertyShouldBe(page, "fontSize", el("tree viewer"), "12"));
      await session.step(147, "And no errors should have been logged", () => noErrors(page));
    }, 146);
  });
  test("The word cloud counts an empty value as a word with no text", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(156, "Given user opens a table \"with-empty\" with:", () => openTableOf(page, "with-empty", [["site","arm"],["A","x"],["A",""],["B","y"]]), [["site","arm"],["A","x"],["A",""],["B","y"]]);
    await session.step(161, "When user adds a word cloud viewer with:", () => addViewerWith(page, "word cloud", [["wordColumnName","arm"]]), [["wordColumnName","arm"]]);
    await session.step(163, "Then no errors should have been logged", () => noErrors(page));
    await session.step(164, "And word cloud viewer should report no error", () => reportsNoError(page, el("word cloud viewer")));
    await session.step(165, "And the \"rows shown\" reading of word cloud viewer should be 3", () => readingIs(page, "rows shown", el("word cloud viewer"), 3));
    await session.step(166, "And the \"words\" reading of word cloud viewer should be 3", () => readingIs(page, "words", el("word cloud viewer"), 3));
  });
  test("A tree over a column of more than 500 categories draws once a filter passes fewer than 500 rows [viewer=tree, message=The Tree viewer requires at least one categorical column with fewer than 500 unique categories]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(169, "When user adds a tree viewer with:", () => addViewerWith(page, "tree", [["hierarchyColumnNames","USUBJID"]]), [["hierarchyColumnNames","USUBJID"]]);
    await session.step(171, "Then \"The Tree viewer requires at least one categorical column with fewer than 500 unique categories\" text in tree viewer should be visible", () => shouldBe(page, el("\"The Tree viewer requires at least one categorical column with fewer than 500 unique categories\" text in tree viewer"), "visible"));
    await session.step(172, "When user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
    await session.step(173, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
    await session.step(174, "And \"The Tree viewer requires at least one categorical column with fewer than 500 unique categories\" text in tree viewer should be absent", () => shouldBe(page, el("\"The Tree viewer requires at least one categorical column with fewer than 500 unique categories\" text in tree viewer"), "absent"));
    await session.step(175, "And no errors should have been logged", () => noErrors(page));
  });
  test("A sunburst over a column of more than 500 categories draws once a filter passes fewer than 500 rows [viewer=sunburst, message=The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(169, "When user adds a sunburst viewer with:", () => addViewerWith(page, "sunburst", [["hierarchyColumnNames","USUBJID"]]), [["hierarchyColumnNames","USUBJID"]]);
    await session.step(171, "Then \"The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories\" text in sunburst viewer should be visible", () => shouldBe(page, el("\"The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories\" text in sunburst viewer"), "visible"));
    await session.step(172, "When user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
    await session.step(173, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
    await session.step(174, "And \"The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories\" text in sunburst viewer should be absent", () => shouldBe(page, el("\"The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories\" text in sunburst viewer"), "absent"));
    await session.step(175, "And no errors should have been logged", () => noErrors(page));
  });
});
