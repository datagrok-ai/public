/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-chrome-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.chrome, public#667, public#719, public#1745, public#2925, public#338]
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
import {clickOn, doubleClickOn, enterInto, isExpanded, pressKey, selectIn, shouldBe, shouldHaveText, typeAtCaret} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, contextPanelOpen, contextPanelShows, openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, loadLayout, noBalloons, noErrors, propertyShouldBe, saveLayoutToServer, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu, stackViewer, viewersTabbed} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Viewer chrome regressions — the context panel, the title and layouts", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-chrome-regressions.feature", import.meta.url);
  test("The gear of a viewer on the first table shows its properties after a viewer on the second table was edited (public#667)", {tag: ["@viewers", "@realizes:viewers.chrome", "@realizes:public#667"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(25, "Given user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(26, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(27, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(28, "When user clicks on settings icon of histogram viewer", () => clickOn(page, el("settings icon of histogram viewer")));
    await session.step(29, "Then the context panel should show \"Histogram\"", () => contextPanelShows(page, "Histogram"));
    await session.step(30, "Given \"X Axis\" category in context panel is expanded", () => isExpanded(page, el("\"X Axis\" category in context panel")));
    await session.step(31, "When user enters \"30\" in \"Bins\" property in context panel", () => enterInto(page, "30", el("\"Bins\" property in context panel")));
    await session.step(32, "Then \"Bins\" property of histogram viewer should be \"30\"", () => propertyShouldBe(page, "Bins", el("histogram viewer"), "30"));
    await session.step(33, "Given user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
    await session.step(34, "When user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
    await session.step(35, "Then the context panel should show \"Scatter plot\"", () => contextPanelShows(page, "Scatter plot"));
    await session.step(36, "And \"Bins\" property in context panel should be absent", () => shouldBe(page, el("\"Bins\" property in context panel"), "absent"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("After a property is edited in the context panel, the gear of another viewer switches the panel to it (public#719)", {tag: ["@viewers", "@realizes:viewers.chrome", "@realizes:public#719"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(42, "Given user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(43, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(44, "When user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
    await session.step(45, "Then the context panel should show \"Scatter plot\"", () => contextPanelShows(page, "Scatter plot"));
    await session.step(46, "Given \"Color\" category in context panel is expanded", () => isExpanded(page, el("\"Color\" category in context panel")));
    await session.step(47, "When user selects \"RACE\" in \"Color\" property in context panel", () => selectIn(page, "RACE", el("\"Color\" property in context panel")));
    await session.step(48, "Then \"Color\" property of scatter plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Color", el("scatter plot viewer"), "RACE"));
    await session.step(49, "When user clicks on settings icon of histogram viewer", () => clickOn(page, el("settings icon of histogram viewer")));
    await session.step(50, "Then the context panel should show \"Histogram\"", () => contextPanelShows(page, "Histogram"));
    await session.step(51, "And \"Marker Type\" property in context panel should be absent", () => shouldBe(page, el("\"Marker Type\" property in context panel"), "absent"));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
  });
  test("A viewer stacked over the grid and closed does not come back as a floating window from the layout (public#1745)", {tag: ["@viewers", "@realizes:viewers.chrome", "@realizes:public#1745"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(57, "Given user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(58, "When user stacks scatter plot viewer onto grid as a tab", () => stackViewer(page, el("scatter plot viewer"), el("grid")));
    await session.step(59, "Then the viewers \"demog-1000, Scatter plot\" should be the tabs of one panel", () => viewersTabbed(page, "demog-1000, Scatter plot"));
    await session.step(60, "When user picks \"General > Close\" from the viewer menu of scatter plot viewer", () => pickFromViewerMenu(page, "General > Close", el("scatter plot viewer")));
    await session.step(61, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
    await session.step(62, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(63, "And user loads the saved layout", () => loadLayout(page));
    await session.step(64, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
    await session.step(65, "And grid should be visible", () => shouldBe(page, el("grid"), "visible"));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
  test("A title typed into the viewer's header becomes its Title and comes back from a layout (public#2925)", {tag: ["@viewers", "@realizes:viewers.chrome", "@realizes:public#2925"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(71, "Given user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(72, "And user sets \"Title\" property of scatter plot viewer to \"Plot\"", () => setProperty(page, "Title", el("scatter plot viewer"), "Plot"));
    await session.step(73, "Then title of scatter plot viewer should have text \"Plot\"", () => shouldHaveText(page, el("title of scatter plot viewer"), "Plot"));
    await session.step(74, "When user double-clicks on title of scatter plot viewer", () => doubleClickOn(page, el("title of scatter plot viewer")));
    await session.step(75, "And user types \"Heights\" at the caret", () => typeAtCaret(page, "Heights"));
    await session.step(76, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(77, "Then title of scatter plot viewer should have text \"Heights\"", () => shouldHaveText(page, el("title of scatter plot viewer"), "Heights"));
    await session.step(78, "And \"Title\" property of scatter plot viewer should be \"Heights\"", () => propertyShouldBe(page, "Title", el("scatter plot viewer"), "Heights"));
    await session.step(79, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(80, "And user closes all views", () => closeAllViews(page));
    await session.step(81, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(82, "And user loads the saved layout", () => loadLayout(page));
    await session.step(83, "Then \"Title\" property of scatter plot viewer should be \"Heights\"", () => propertyShouldBe(page, "Title", el("scatter plot viewer"), "Heights"));
    await session.step(84, "And title of scatter plot viewer should have text \"Heights\"", () => shouldHaveText(page, el("title of scatter plot viewer"), "Heights"));
    await session.step(85, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout with a viewer bound to another table applies over the tables opened anew with no warning (public#338)", {tag: ["@viewers", "@realizes:viewers.chrome", "@realizes:public#338"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(90, "Given user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await session.step(91, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(92, "When user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
    await session.step(93, "And user selects \"demog-1000\" in \"Table\" property in context panel", () => selectIn(page, "demog-1000", el("\"Table\" property in context panel")));
    await session.step(94, "Then scatter plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("scatter plot viewer"), "demog-1000"));
    await session.step(95, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(96, "And user closes all views", () => closeAllViews(page));
    await session.step(97, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(98, "And user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await session.step(99, "And user loads the saved layout", () => loadLayout(page));
    await session.step(100, "Then scatter plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("scatter plot viewer"), "demog-1000"));
    await session.step(101, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(102, "And no errors should have been logged", () => noErrors(page));
  });
});
