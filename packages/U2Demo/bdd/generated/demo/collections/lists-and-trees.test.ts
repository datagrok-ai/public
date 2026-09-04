/* ---
generated: features/demo/collections/lists-and-trees.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.list, u2.tree]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, expand, pressKeyIn, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Virtual lists and trees", () => {
  const session = feature(test);
  test("Selecting in a virtual list of 100,000 rows", {tag: ["@demo", "@realizes:u2.list", "@realizes:u2.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Lists\" demo page", () => openDemoPage(page, "Lists"));
    enter(page, "U2 Demo");
    await test.step("Then value of selectedIndex readout should have text \"-1\"", () => shouldHaveText(page, el("value of selectedIndex readout"), "-1"));
    await test.step("When user clicks on \"item 5\" item in list", () => clickOn(page, el("\"item 5\" item in list")));
    await test.step("Then \"item 5\" item in list should be selected", () => shouldBe(page, el("\"item 5\" item in list"), "selected"));
    await test.step("And value of selectedIndex readout should have text \"5\"", () => shouldHaveText(page, el("value of selectedIndex readout"), "5"));
    await test.step("When user presses ArrowDown in list", () => pressKeyIn(page, "ArrowDown", el("list")));
    await test.step("Then value of selectedIndex readout should have text \"6\"", () => shouldHaveText(page, el("value of selectedIndex readout"), "6"));
    await test.step("And \"item 6\" item in list should be selected", () => shouldBe(page, el("\"item 6\" item in list"), "selected"));
  });
  test("Expanding branches, lazy ones included", {tag: ["@demo", "@realizes:u2.list", "@realizes:u2.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Trees\" demo page", () => openDemoPage(page, "Trees"));
    enter(page, "U2 Demo");
    await test.step("Then \"Demo project\" tree node should be collapsed", () => shouldBe(page, el("\"Demo project\" tree node"), "collapsed"));
    await test.step("When user expands \"Demo project\" tree node", () => expand(page, el("\"Demo project\" tree node")));
    await test.step("Then Tables tree node should be visible", () => shouldBe(page, el("Tables tree node"), "visible"));
    await test.step("When user expands Tables tree node", () => expand(page, el("Tables tree node")));
    await test.step("And user clicks on demog tree node", () => clickOn(page, el("demog tree node")));
    await test.step("Then value of selectedNode readout should have text \"demog\"", () => shouldHaveText(page, el("value of selectedNode readout"), "demog"));
    await test.step("And demog tree node should be selected", () => shouldBe(page, el("demog tree node"), "selected"));
    await test.step("When user expands \"Server (lazy)\" tree node", () => expand(page, el("\"Server (lazy)\" tree node")));
    await test.step("Then dataset-0.csv tree node should be visible", () => shouldBe(page, el("dataset-0.csv tree node"), "visible"));
  });
  test("Revealing a path expands its ancestors", {tag: ["@demo", "@realizes:u2.list", "@realizes:u2.tree"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Trees\" demo page", () => openDemoPage(page, "Trees"));
    enter(page, "U2 Demo");
    await test.step("When user clicks on \"Reveal demog\" button", () => clickOn(page, el("\"Reveal demog\" button")));
    await test.step("Then demog tree node should be selected", () => shouldBe(page, el("demog tree node"), "selected"));
    await test.step("And value of selectedNode readout should have text \"demog\"", () => shouldHaveText(page, el("value of selectedNode readout"), "demog"));
  });
});
