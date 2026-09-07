/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/forms/property-grid.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.property-grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {collapse, enterInto, expand, selectIn, shouldBe, shouldContainText, shouldHaveText, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The property grid", () => {
  const session = feature(test, "features/demo/forms/property-grid.feature", import.meta.url);
  test("Editing a property replaces the value record", {tag: ["@demo", "@realizes:u2.property-grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(6, "Given user opens the \"Property grid\" demo page", () => openDemoPage(page, "Property grid"));
    enter(page, "U2 Demo");
    await session.step(9, "Then value of \"last change\" readout should have text \"(none)\"", () => shouldHaveText(page, el("value of \"last change\" readout"), "(none)"));
    await session.step(10, "When user types \"Heatmap\" into title property", () => typeInto(page, "Heatmap", el("title property")));
    await session.step(11, "Then value of \"last change\" readout should have text \"title → Heatmap\"", () => shouldHaveText(page, el("value of \"last change\" readout"), "title → Heatmap"));
    await session.step(12, "And value of values readout should contain text \"\\\"title\\\":\\\"Heatmap\\\"\"", () => shouldContainText(page, el("value of values readout"), "\"title\":\"Heatmap\""));
    await session.step(13, "When user unchecks showLegend property", () => uncheck(page, el("showLegend property")));
    await session.step(14, "Then value of \"last change\" readout should have text \"showLegend → false\"", () => shouldHaveText(page, el("value of \"last change\" readout"), "showLegend → false"));
    await session.step(15, "When user selects \"top\" in position property", () => selectIn(page, "top", el("position property")));
    await session.step(16, "Then value of \"last change\" readout should have text \"position → top\"", () => shouldHaveText(page, el("value of \"last change\" readout"), "position → top"));
    await session.step(17, "When user enters \"0.5\" into opacity property", () => enterInto(page, "0.5", el("opacity property")));
    await session.step(18, "Then value of \"last change\" readout should have text \"opacity → 0.5\"", () => shouldHaveText(page, el("value of \"last change\" readout"), "opacity → 0.5"));
  });
  test("Categories collapse", {tag: ["@demo", "@realizes:u2.property-grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(6, "Given user opens the \"Property grid\" demo page", () => openDemoPage(page, "Property grid"));
    enter(page, "U2 Demo");
    await session.step(21, "Then Layout category should be expanded", () => shouldBe(page, el("Layout category"), "expanded"));
    await session.step(22, "And width property should be visible", () => shouldBe(page, el("width property"), "visible"));
    await session.step(23, "When user collapses Layout category", () => collapse(page, el("Layout category")));
    await session.step(24, "Then Layout category should be collapsed", () => shouldBe(page, el("Layout category"), "collapsed"));
    await session.step(25, "And width property should be hidden", () => shouldBe(page, el("width property"), "hidden"));
    await session.step(26, "When user expands Layout category", () => expand(page, el("Layout category")));
    await session.step(27, "Then width property should be visible", () => shouldBe(page, el("width property"), "visible"));
  });
});
