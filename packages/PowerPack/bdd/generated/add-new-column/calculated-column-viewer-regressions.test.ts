/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/calculated-column-viewer-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.dialogs.add-new-column, public#434]
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
import {clickOn, holdsCode, isExpanded, replaceCode, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, columnIsCurrentObject} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, readingAsRemembered, readingHigherThanRemembered, readingNotAsRemembered, rememberReading, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Calculated column regressions — a viewer follows a formula changed in Edit Column Formula", () => {
  const session = feature(test, "features/add-new-column/calculated-column-viewer-regressions.feature", import.meta.url);
  test("A viewer on a calculated column redraws after the formula is changed, and keeps the new values (public#434)", {tag: ["@realizes:powerpack.dialogs.add-new-column", "@realizes:public#434"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(18, "Given user adds a calculated column \"W2\" with formula \"${WEIGHT} / 2\"", () => addCalculated(page, "W2", "${WEIGHT} / 2"));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","HEIGHT"],["Y","W2"]]), [["X","HEIGHT"],["Y","W2"]]);
    await session.step(22, "When user remembers the \"y axis max\" reading of scatter plot viewer", () => rememberReading(page, "y axis max", el("scatter plot viewer")));
    await session.step(23, "Given the \"W2\" column is the current object", () => columnIsCurrentObject(page, "W2"));
    await session.step(24, "Then the context panel should show \"W2\"", () => contextPanelShows(page, "W2"));
    await session.step(25, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
    await session.step(26, "Then formula pane editor in Formula pane in context panel should hold the code \"${WEIGHT} / 2\"", () => holdsCode(page, el("formula pane editor in Formula pane in context panel"), "${WEIGHT} / 2"));
    await session.step(27, "When user clicks on \"Edit in dialog\" button in context panel", () => clickOn(page, el("\"Edit in dialog\" button in context panel")));
    await session.step(28, "Then \"Edit Column Formula\" dialog should be visible", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "visible"));
    await session.step(29, "When user replaces the code of formula editor with \"${WEIGHT} * 2\"", () => replaceCode(page, el("formula editor"), "${WEIGHT} * 2"));
    await session.step(30, "And user clicks on OK button in \"Edit Column Formula\" dialog", () => clickOn(page, el("OK button in \"Edit Column Formula\" dialog")));
    await session.step(31, "Then \"Edit Column Formula\" dialog should be hidden", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "hidden"));
    await session.step(32, "And the \"y axis max\" reading of scatter plot viewer should be higher than remembered", () => readingHigherThanRemembered(page, "y axis max", el("scatter plot viewer")));
    await session.step(33, "When user remembers the \"y axis max\" reading of scatter plot viewer", () => rememberReading(page, "y axis max", el("scatter plot viewer")));
    await session.step(34, "And user sets \"Y\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Y", el("scatter plot viewer"), "AGE"));
    await session.step(35, "Then the \"y axis max\" reading of scatter plot viewer should not be as remembered", () => readingNotAsRemembered(page, "y axis max", el("scatter plot viewer")));
    await session.step(36, "When user sets \"Y\" property of scatter plot viewer to \"W2\"", () => setProperty(page, "Y", el("scatter plot viewer"), "W2"));
    await session.step(37, "Then the \"y axis max\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "y axis max", el("scatter plot viewer")));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
});
