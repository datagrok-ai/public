/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/containers/layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.accordion, u2.breadcrumbs, u2.toolbar, u2.splitter]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, collapse, dragTo, expand, shouldBe, shouldContainText, shouldHaveText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Layout containers", () => {
  const session = feature(test, "features/demo/containers/layout.feature", import.meta.url);
  test("Accordion panes expand and collapse", {tag: ["@demo", "@realizes:u2.accordion", "@realizes:u2.breadcrumbs", "@realizes:u2.toolbar", "@realizes:u2.splitter"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user opens the \"Layout\" demo page", () => openDemoPage(page, "Layout"));
    enter(page, "U2 Demo");
    await session.step(10, "Then General pane should be expanded", () => shouldBe(page, el("General pane"), "expanded"));
    await session.step(11, "And Advanced pane should be collapsed", () => shouldBe(page, el("Advanced pane"), "collapsed"));
    await session.step(12, "When user expands Advanced pane", () => expand(page, el("Advanced pane")));
    await session.step(13, "Then Advanced pane should be expanded", () => shouldBe(page, el("Advanced pane"), "expanded"));
    await session.step(14, "And Advanced pane should contain text \"Lazily built on first expand\"", () => shouldContainText(page, el("Advanced pane"), "Lazily built on first expand"));
    await session.step(15, "When user collapses Advanced pane", () => collapse(page, el("Advanced pane")));
    await session.step(16, "Then Advanced pane should be collapsed", () => shouldBe(page, el("Advanced pane"), "collapsed"));
  });
  test("Breadcrumbs truncate the path", {tag: ["@demo", "@realizes:u2.accordion", "@realizes:u2.breadcrumbs", "@realizes:u2.toolbar", "@realizes:u2.splitter"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user opens the \"Layout\" demo page", () => openDemoPage(page, "Layout"));
    enter(page, "U2 Demo");
    await session.step(19, "Then value of path readout should have text \"Home > Projects > Demo > Tables > demog\"", () => shouldHaveText(page, el("value of path readout"), "Home > Projects > Demo > Tables > demog"));
    await session.step(20, "When user clicks on Projects breadcrumb", () => clickOn(page, el("Projects breadcrumb")));
    await session.step(21, "Then value of path readout should have text \"Home > Projects\"", () => shouldHaveText(page, el("value of path readout"), "Home > Projects"));
    await session.step(22, "And breadcrumbs should not contain text \"demog\"", () => shouldNotContainText(page, el("breadcrumbs"), "demog"));
  });
  test("A panel-local toolbar with buttons, a toggle and a menu", {tag: ["@demo", "@realizes:u2.accordion", "@realizes:u2.breadcrumbs", "@realizes:u2.toolbar", "@realizes:u2.splitter"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user opens the \"Layout\" demo page", () => openDemoPage(page, "Layout"));
    enter(page, "U2 Demo");
    await session.step(25, "When user clicks on Run button in toolbar", () => clickOn(page, el("Run button in toolbar")));
    await session.step(26, "Then value of toolbar readout should contain text \"Run clicked\"", () => shouldContainText(page, el("value of toolbar readout"), "Run clicked"));
    await session.step(27, "When user clicks on Save button in toolbar", () => clickOn(page, el("Save button in toolbar")));
    await session.step(28, "Then value of toolbar readout should contain text \"Saved\"", () => shouldContainText(page, el("value of toolbar readout"), "Saved"));
    await session.step(29, "When user clicks on Compact button in toolbar", () => clickOn(page, el("Compact button in toolbar")));
    await session.step(30, "Then Compact button in toolbar should be selected", () => shouldBe(page, el("Compact button in toolbar"), "selected"));
    await session.step(31, "And value of toolbar readout should contain text \"compact true\"", () => shouldContainText(page, el("value of toolbar readout"), "compact true"));
    await session.step(32, "When user clicks on Help button in toolbar", () => clickOn(page, el("Help button in toolbar")));
    await session.step(33, "And user clicks on \"About u2\" menu item", () => clickOn(page, el("\"About u2\" menu item")));
    await session.step(34, "Then value of toolbar readout should contain text \"u2: next-gen Datagrok UI library\"", () => shouldContainText(page, el("value of toolbar readout"), "u2: next-gen Datagrok UI library"));
  });
  test("The splitter reports its sizes", {tag: ["@demo", "@realizes:u2.accordion", "@realizes:u2.breadcrumbs", "@realizes:u2.toolbar", "@realizes:u2.splitter"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user opens the \"Layout\" demo page", () => openDemoPage(page, "Layout"));
    enter(page, "U2 Demo");
    await session.step(37, "Then value of sizes readout should have text \"0.3 / 0.7\"", () => shouldHaveText(page, el("value of sizes readout"), "0.3 / 0.7"));
    await session.step(38, "When user drags sash in splitter to Content heading", () => dragTo(page, el("sash in splitter"), el("Content heading")));
    await session.step(39, "Then value of sizes readout should not contain text \"0.3 / 0.7\"", () => shouldNotContainText(page, el("value of sizes readout"), "0.3 / 0.7"));
  });
});
