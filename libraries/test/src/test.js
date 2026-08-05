var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
import { testData } from './dataframe-utils';
import { changeOptionsSaveLayout, filterAsync, loadLayout, selectFilterChangeCurrent, testViewerInternal } from './test-viewer-utils';
const STANDART_TIMEOUT = 30000;
const BENCHMARK_TIMEOUT = 10800000;
const stdLog = console.log.bind(console);
const stdInfo = console.info.bind(console);
const stdWarn = console.warn.bind(console);
const stdError = console.error.bind(console);
export const tests = {};
const autoTestsCatName = 'Auto Tests';
const demoCatName = 'Demo';
const detectorsCatName = 'Detectors';
const coreCatName = 'Core';
const wasRegistered = {};
export let currentCategory;
export var assure;
(function (assure) {
    function notNull(value, name) {
        if (value == null)
            throw new Error(`${name == null ? 'Value' : name} not defined`);
    }
    assure.notNull = notNull;
})(assure || (assure = {}));
export class TestContext {
    constructor(catchUnhandled, report, returnOnFail) {
        this.catchUnhandled = true;
        this.report = false;
        this.returnOnFail = false;
        this.closeAll = true;
        if (catchUnhandled !== undefined)
            this.catchUnhandled = catchUnhandled;
        if (report !== undefined)
            this.report = report;
        if (returnOnFail !== undefined)
            this.returnOnFail = returnOnFail;
    }
    ;
}
export class Test {
    constructor(category, name, test, options) {
        var _a;
        this.category = category;
        this.name = name;
        options !== null && options !== void 0 ? options : (options = {});
        (_a = options.timeout) !== null && _a !== void 0 ? _a : (options.timeout = STANDART_TIMEOUT);
        this.options = options;
        this.test = () => __awaiter(this, void 0, void 0, function* () {
            return new Promise((resolve, reject) => __awaiter(this, void 0, void 0, function* () {
                var _b;
                let result = '';
                try {
                    if (DG.Test.isInDebug)
                        debugger;
                    let res = yield test();
                    try {
                        result = (_b = res === null || res === void 0 ? void 0 : res.toString()) !== null && _b !== void 0 ? _b : '';
                    }
                    catch (e) {
                        result = 'Can\'t convert test\'s result to string';
                        console.error(`Can\'t convert test\'s result to string in the ${this.category}:${this.name} test`);
                    }
                }
                catch (e) {
                    reject(e);
                }
                resolve(result);
            }));
        });
    }
}
export class Category {
}
export class NodeTestExecutionOptions {
}
export class TestExecutionOptions {
}
export function testEvent(event, handler, trigger, ms = 0, reason = `timeout`) {
    return __awaiter(this, void 0, void 0, function* () {
        return new Promise((resolve, reject) => {
            const sub = event.subscribe((args) => {
                try {
                    handler(args);
                    resolve('OK');
                }
                catch (e) {
                    reject(e);
                }
                finally {
                    sub.unsubscribe();
                    clearTimeout(timeout);
                }
            });
            const timeout = setTimeout(() => {
                sub.unsubscribe();
                // eslint-disable-next-line prefer-promise-reject-errors
                reject(reason);
            }, ms);
            trigger();
        });
    });
}
export function testEventAsync(event, handler, trigger, ms = 0, reason = `timeout`) {
    return __awaiter(this, void 0, void 0, function* () {
        return new Promise((resolve, reject) => {
            const sub = event.subscribe((args) => {
                handler(args).then(() => {
                    resolve('OK');
                }).catch((e) => {
                    reject(e);
                }).finally(() => {
                    sub.unsubscribe();
                    clearTimeout(timeout);
                });
            });
            const timeout = setTimeout(() => {
                sub.unsubscribe();
                // eslint-disable-next-line prefer-promise-reject-errors
                reject(reason);
            }, ms);
            trigger();
        });
    });
}
export function test(name, test, options) {
    if (tests[currentCategory] == undefined)
        tests[currentCategory] = {};
    if (tests[currentCategory].tests == undefined)
        tests[currentCategory].tests = [];
    tests[currentCategory].tests.push(new Test(currentCategory, name, test, options));
}
/* Tests two objects for equality, throws an exception if they are not equal. */
export function expect(actual, expected = true, error) {
    if (error)
        error = `${error}, `;
    else
        error = '';
    if (actual !== expected)
        throw new Error(`${error}Expected "${expected}", got "${actual}"`);
}
export function expectFloat(actual, expected, tolerance = 0.001, error) {
    if ((actual === Number.POSITIVE_INFINITY && expected === Number.POSITIVE_INFINITY) ||
        (actual === Number.NEGATIVE_INFINITY && expected === Number.NEGATIVE_INFINITY) ||
        (actual === Number.NaN && expected === Number.NaN) || (isNaN(actual) && isNaN(expected)))
        return;
    const areEqual = Math.abs(actual - expected) < tolerance;
    expect(areEqual, true, `${error !== null && error !== void 0 ? error : ''} (tolerance = ${tolerance}; a = ${actual}, e = ${expected})`);
    if (!areEqual)
        throw new Error(`Expected ${expected}, got ${actual} (tolerance = ${tolerance})`);
}
export function expectTable(actual, expected, error) {
    const expectedRowCount = expected.rowCount;
    const actualRowCount = actual.rowCount;
    expect(actualRowCount, expectedRowCount, `${error !== null && error !== void 0 ? error : ''}, row count`);
    for (const column of expected.columns) {
        const actualColumn = actual.columns.byName(column.name);
        if (actualColumn == null)
            throw new Error(`Column ${column.name} not found`);
        if (actualColumn.type != column.type)
            throw new Error(`Column ${column.name} type expected ${column.type} got ${actualColumn.type}`);
        for (let i = 0; i < expectedRowCount; i++) {
            const value = column.get(i);
            const actualValue = actualColumn.get(i);
            if (column.type == DG.TYPE.FLOAT)
                expectFloat(actualValue, value, 0.0001, error);
            else if (column.type == DG.TYPE.DATE_TIME)
                expect(actualValue.isSame(value), true, error);
            else
                expect(actualValue, value, error);
        }
    }
}
export function expectObject(actual, expected) {
    for (const [expectedKey, expectedValue] of Object.entries(expected)) {
        if (!actual.hasOwnProperty(expectedKey))
            throw new Error(`Expected property "${expectedKey}" not found`);
        const actualValue = actual[expectedKey];
        if (actualValue instanceof Array && expectedValue instanceof Array)
            expectArray(actualValue, expectedValue);
        else if (actualValue instanceof Object && expectedValue instanceof Object)
            expectObject(actualValue, expectedValue);
        else if (Number.isFinite(actualValue) && Number.isFinite(expectedValue))
            expectFloat(actualValue, expectedValue);
        else if (actualValue != expectedValue)
            throw new Error(`Expected (${expectedValue}) for key '${expectedKey}', got (${actualValue})`);
    }
}
export function expectArray(actual, expected) {
    const actualLength = actual.length;
    const expectedLength = expected.length;
    if (actualLength != expectedLength) {
        throw new Error(`Arrays are of different length: actual array length is ${actualLength} ` +
            `and expected array length is ${expectedLength}`);
    }
    for (let i = 0; i < actualLength; i++) {
        if (actual[i] instanceof Array && expected[i] instanceof Array)
            expectArray(actual[i], expected[i]);
        else if (actual[i] instanceof Object && expected[i] instanceof Object)
            expectObject(actual[i], expected[i]);
        else if (actual[i] != expected[i])
            throw new Error(`Expected ${expected[i]} at position ${i}, got ${actual[i]}`);
    }
}
/* Defines a test suite. */
export function category(category, tests_, options) {
    var _a;
    currentCategory = category;
    tests_();
    if (tests[currentCategory]) {
        tests[currentCategory].clear = (_a = options === null || options === void 0 ? void 0 : options.clear) !== null && _a !== void 0 ? _a : true;
        tests[currentCategory].timeout = options === null || options === void 0 ? void 0 : options.timeout;
        tests[currentCategory].benchmarks = options === null || options === void 0 ? void 0 : options.benchmarks;
        tests[currentCategory].stressTests = options === null || options === void 0 ? void 0 : options.stressTests;
        tests[currentCategory].owner = options === null || options === void 0 ? void 0 : options.owner;
        tests[currentCategory].node = options === null || options === void 0 ? void 0 : options.node;
    }
}
/* Defines a function to be executed before the tests in this category are executed. */
export function before(before) {
    if (tests[currentCategory] == undefined)
        tests[currentCategory] = {};
    tests[currentCategory].before = before;
}
/* Defines a function to be executed after the tests in this category are executed. */
export function after(after) {
    if (tests[currentCategory] == undefined)
        tests[currentCategory] = {};
    tests[currentCategory].after = after;
}
function addNamespace(s, f) {
    return s.replace(new RegExp(f.name, 'gi'), f.nqName);
}
/** Whether a test matches the node/browser split requested by the run options. */
function matchesNodeTarget(t, cat, options) {
    var _a, _b, _c;
    const isNode = (_c = (_b = (_a = t.options) === null || _a === void 0 ? void 0 : _a.node) !== null && _b !== void 0 ? _b : cat.node) !== null && _c !== void 0 ? _c : false;
    return options.nodeOnly ? isNode : options.excludeNodeTests ? !isNode : true;
}
export function initAutoTests(package_, module) {
    var _a, _b, _c, _d, _e, _f, _g, _h, _j, _k;
    return __awaiter(this, void 0, void 0, function* () {
        const packageId = package_.id;
        if (wasRegistered[packageId])
            return;
        const moduleTests = module ? module.tests : tests;
        if (package_.name === 'DevTools' || (!!module && module._package.name === 'DevTools')) {
            for (const f of window.dartTests) {
                const arr = f.name.split(/\s*\|\s*!/g);
                let name = (_a = arr.pop()) !== null && _a !== void 0 ? _a : f.name;
                let cat = arr.length ? coreCatName + ': ' + arr.join(': ') : coreCatName;
                let fullName = name.split(' | ');
                name = fullName[fullName.length - 1];
                fullName.unshift(cat);
                fullName.pop();
                cat = fullName.join(': ');
                if (moduleTests[cat] === undefined)
                    moduleTests[cat] = { tests: [], clear: true };
                moduleTests[cat].tests.push(new Test(cat, name, f.test, { isAggregated: false, timeout: (_c = (_b = f.options) === null || _b === void 0 ? void 0 : _b.timeout) !== null && _c !== void 0 ? _c : STANDART_TIMEOUT, skipReason: (_d = f.options) === null || _d === void 0 ? void 0 : _d.skipReason, owner: (_e = f.options) === null || _e === void 0 ? void 0 : _e.owner, benchmark: (_g = (_f = f.options) === null || _f === void 0 ? void 0 : _f.benchmark) !== null && _g !== void 0 ? _g : false }));
            }
        }
        const moduleAutoTests = [];
        const moduleDemo = [];
        const moduleDetectors = [];
        const packFunctions = yield grok.dapi.functions.filter(`package.id = "${packageId}"`).list();
        const reg = new RegExp(/skip:\s*([^,\s]+)|wait:\s*(\d+)|cat:\s*([^,\s]+)|timeout:\s*(\d+)/g);
        for (const f of packFunctions) {
            const tests = f.options['test'];
            const demo = f.options['demoPath'];
            if ((tests && Array.isArray(tests) && tests.length)) {
                for (let i = 0; i < tests.length; i++) {
                    const res = tests[i].matchAll(reg);
                    const map = {};
                    Array.from(res).forEach((arr) => {
                        if (arr[0].startsWith('skip'))
                            map['skip'] = arr[1];
                        else if (arr[0].startsWith('wait'))
                            map['wait'] = parseInt(arr[2]);
                        else if (arr[0].startsWith('cat'))
                            map['cat'] = arr[3];
                        else if (arr[0].startsWith('timeout'))
                            map['timeout'] = parseInt(arr[4]);
                    });
                    const test = new Test((_h = map.cat) !== null && _h !== void 0 ? _h : autoTestsCatName, tests.length === 1 ? f.name : `${f.name} ${i + 1}`, () => __awaiter(this, void 0, void 0, function* () {
                        const res = yield grok.functions.eval(addNamespace(tests[i], f));
                        if (map.wait)
                            yield delay(map.wait);
                        // eslint-disable-next-line no-throw-literal
                        if (typeof res === 'boolean' && !res)
                            throw `Failed: ${tests[i]}, expected true, got ${res}`;
                    }), { skipReason: map.skip, timeout: DG.Test.isInBenchmark ? (_j = map.benchmarkTimeout) !== null && _j !== void 0 ? _j : BENCHMARK_TIMEOUT : (_k = map.timeout) !== null && _k !== void 0 ? _k : STANDART_TIMEOUT,
                        // Query tests evaluate server-side (expressions like ExpectTable(Query(), OpenFile(...))
                        // resolve through core functions available in the Node runtime), so they can run headless;
                        // tests of client JS functions need the browser-loaded package.
                        node: f instanceof DG.DataQuery });
                    if (map.cat) {
                        const cat = map.cat;
                        if (moduleTests[cat] === undefined)
                            moduleTests[cat] = { tests: [], clear: true };
                        // only before/after can be defined in ts files tests under the category
                        if (!moduleTests[cat].tests)
                            moduleTests[cat].tests = [];
                        moduleTests[cat].tests.push(test);
                    }
                    else
                        moduleAutoTests.push(test);
                }
            }
            if (demo) {
                const wait = f.options['demoWait'] ? parseInt(f.options['demoWait']) : undefined;
                const test = new Test(demoCatName, f.friendlyName, () => __awaiter(this, void 0, void 0, function* () {
                    yield delay(300);
                    grok.shell.clearLastError();
                    yield f.apply();
                    yield delay(wait ? wait : 2000);
                    const unhandled = yield grok.shell.lastError;
                    if (unhandled)
                        throw new Error(unhandled);
                }), { skipReason: f.options['demoSkip'] });
                moduleDemo.push(test);
            }
            if (f.hasTag('semTypeDetector')) {
                let detectorsTestData = testData;
                if (f.options['testData']) {
                    detectorsTestData = yield grok.data.files.openTable(`System:AppData/${package_.nqName}/${f.options['testData']}`);
                }
                const test = new Test(detectorsCatName, f.friendlyName, () => __awaiter(this, void 0, void 0, function* () {
                    const arr = [];
                    console.log(`System:AppData/${package_.nqName}/${f.options['testData']}`);
                    for (const col of detectorsTestData.clone().columns) {
                        const res = yield f.apply([col]);
                        arr.push(res || col.semType);
                    }
                    const resArr = arr.filter((i) => i);
                    expect(resArr.length, 1);
                    if (f.options['testDataColumnName'])
                        expect(resArr[0], f.options['testDataColumnName']);
                }), { skipReason: f.options['skipTest'] });
                moduleDetectors.push(test);
            }
        }
        wasRegistered[packageId] = true;
        if (moduleAutoTests.length > 0)
            moduleTests[autoTestsCatName] = { tests: moduleAutoTests, clear: true };
        if (moduleDemo.length > 0)
            moduleTests[demoCatName] = { tests: moduleDemo, clear: true };
        if (moduleDetectors.length > 0)
            moduleTests[detectorsCatName] = { tests: moduleDetectors, clear: false };
    });
}
function redefineConsole() {
    const logs = [];
    console.log = (...args) => {
        logs.push(...args);
        stdLog(...args);
    };
    console.info = (...args) => {
        logs.push(...args);
        stdInfo(...args);
    };
    console.warn = (...args) => {
        logs.push(...args);
        stdWarn(...args);
    };
    console.error = (...args) => {
        logs.push(...args);
        stdError(...args);
    };
    return logs;
}
function resetConsole() {
    console.log = stdLog;
    console.info = stdInfo;
    console.warn = stdWarn;
    console.error = stdError;
}
export function runTests(options) {
    var _a, _b;
    var _c;
    return __awaiter(this, void 0, void 0, function* () {
        const package_ = (options === null || options === void 0 ? void 0 : options.nodeOptions) ? options.nodeOptions.package : grok.functions.getCurrentCall().func.package;
        if (!package_)
            throw new Error('Can\'t run tests outside of the package');
        const match = (_a = package_.packageOwner) === null || _a === void 0 ? void 0 : _a.match(/<([^>]*)>/);
        const packageOwner = match ? match[1] : '';
        if (package_ != undefined)
            yield initAutoTests(package_);
        const results = [];
        console.log(`Running tests...`);
        console.log(options);
        options !== null && options !== void 0 ? options : (options = {});
        (_b = (_c = options).testContext) !== null && _b !== void 0 ? _b : (_c.testContext = new TestContext());
        grok.shell.clearLastError();
        const logs = redefineConsole();
        yield invokeTests(tests, options);
        for (let r of results) {
            r.result = r.result.toString().replace(/"/g, '\'');
            if (r.logs != undefined)
                r.logs = r.logs.toString().replace(/"/g, '\'');
        }
        return results;
        function invokeCategoryMethod(method, category) {
            return __awaiter(this, void 0, void 0, function* () {
                let invokationResult = undefined;
                try {
                    if (method !== undefined) {
                        yield timeout(() => __awaiter(this, void 0, void 0, function* () {
                            yield method();
                        }), 100000, `before ${category}: timeout error`);
                    }
                }
                catch (x) {
                    invokationResult = yield getResult(x);
                }
                return invokationResult;
            });
        }
        function invokeTestsInCategory(category, options, isTargetCategory) {
            var _a, _b, _c, _d, _e, _f, _g, _h, _j, _k, _l, _m, _o, _p, _q, _r, _s, _t, _u;
            return __awaiter(this, void 0, void 0, function* () {
                let t = ((_a = category.tests) !== null && _a !== void 0 ? _a : []).filter((e) => matchesNodeTarget(e, category, options));
                const res = [];
                // let memoryUsageBefore = (window?.performance as any)?.memory?.usedJSHeapSize;
                const widgetsBefore = getWidgetsCountSafe();
                if (category.clear) {
                    let skippingTests = isTargetCategory && options.skipToTest != undefined;
                    for (let i = 0; i < t.length; i++) {
                        if (t[i].options) {
                            if (((_b = t[i].options) === null || _b === void 0 ? void 0 : _b.benchmark) === undefined) {
                                if (!t[i].options)
                                    t[i].options = {};
                                t[i].options.benchmark = (_c = category.benchmarks) !== null && _c !== void 0 ? _c : false;
                            }
                        }
                        let test = t[i];
                        if (options.test)
                            if (options.test.toLowerCase() !== test.name.toLowerCase())
                                continue;
                        if (skippingTests) {
                            if ((options === null || options === void 0 ? void 0 : options.skipToTest) != undefined && test.name.toLowerCase().trim() === (options === null || options === void 0 ? void 0 : options.skipToTest.toLowerCase().trim())) {
                                // Found the target test, stop skipping after this one
                                skippingTests = false;
                            }
                            else
                                continue;
                        }
                        if (test === null || test === void 0 ? void 0 : test.options) {
                            test.options.owner = (_g = (_f = (_e = (_d = t[i].options) === null || _d === void 0 ? void 0 : _d.owner) !== null && _e !== void 0 ? _e : category === null || category === void 0 ? void 0 : category.owner) !== null && _f !== void 0 ? _f : packageOwner) !== null && _g !== void 0 ? _g : '';
                        }
                        // let isGBEnable = (window as any).gc && test.options?.skipReason == undefined;
                        // console.log(`********${isGBEnable}`);
                        // if (isGBEnable)
                        //   await (window as any).gc();
                        // memoryUsageBefore = (window?.performance as any)?.memory?.usedJSHeapSize;
                        let testRun = yield execTest(test, options === null || options === void 0 ? void 0 : options.test, logs, DG.Test.isInBenchmark ? (_j = (_h = t[i].options) === null || _h === void 0 ? void 0 : _h.benchmarkTimeout) !== null && _j !== void 0 ? _j : BENCHMARK_TIMEOUT : (_l = (_k = t[i].options) === null || _k === void 0 ? void 0 : _k.timeout) !== null && _l !== void 0 ? _l : STANDART_TIMEOUT, package_.name, options.verbose);
                        // if (isGBEnable)
                        //   await (window as any).gc();
                        if (testRun) {
                            res.push(Object.assign(Object.assign({}, testRun), { widgetsDifference: getWidgetsCountSafe() - widgetsBefore }));
                            // Return early if returnOnFail is set and test failed (but ignore failure for the skipToTest test itself)
                            if (options.returnOnFail && options.skipToTest !== test.name && !testRun.success && !testRun.skipped)
                                return res;
                        }
                        // res.push({ ...testRun, memoryDelta: (window?.performance as any)?.memory?.usedJSHeapSize - memoryUsageBefore, widgetsDelta: getWidgetsCountSafe() - widgetsBefore });
                        if (!options.nodeOptions && ((_m = options.testContext) === null || _m === void 0 ? void 0 : _m.closeAll) !== false) {
                            grok.shell.closeAll();
                            DG.Balloon.closeAll();
                        }
                    }
                }
                else {
                    let skippingTests = isTargetCategory && options.skipToTest != undefined;
                    for (let i = 0; i < t.length; i++) {
                        let test = t[i];
                        if (options.test)
                            if (options.test.toLowerCase() !== test.name.toLowerCase())
                                continue;
                        if (skippingTests) {
                            if ((options === null || options === void 0 ? void 0 : options.skipToTest) != undefined && test.name.toLowerCase().trim() === (options === null || options === void 0 ? void 0 : options.skipToTest.toLowerCase().trim())) {
                                // Found the target test, stop skipping after this one
                                skippingTests = false;
                            }
                            continue; // Skip this test (including the target)
                        }
                        if (test === null || test === void 0 ? void 0 : test.options) {
                            test.options.owner = (_r = (_q = (_p = (_o = t[i].options) === null || _o === void 0 ? void 0 : _o.owner) !== null && _p !== void 0 ? _p : category === null || category === void 0 ? void 0 : category.owner) !== null && _q !== void 0 ? _q : packageOwner) !== null && _r !== void 0 ? _r : '';
                        }
                        // let isGBEnable = (window as any).gc && test.options?.skipReason == undefined;
                        // console.log(`********${isGBEnable}`);
                        // if (isGBEnable)
                        //   await (window as any).gc();
                        // memoryUsageBefore = (window?.performance as any)?.memory?.usedJSHeapSize;
                        let testRun = yield execTest(test, options === null || options === void 0 ? void 0 : options.test, logs, DG.Test.isInBenchmark ? (_t = (_s = t[i].options) === null || _s === void 0 ? void 0 : _s.benchmarkTimeout) !== null && _t !== void 0 ? _t : BENCHMARK_TIMEOUT : (_u = t[i].options) === null || _u === void 0 ? void 0 : _u.timeout, package_.name, options.verbose);
                        // if (isGBEnable)
                        //   await (window as any).gc();
                        if (testRun) {
                            res.push(Object.assign(Object.assign({}, testRun), { widgetsDifference: getWidgetsCountSafe() - widgetsBefore }));
                            // Return early if returnOnFail is set and test failed (but ignore failure for the skipToTest test itself)
                            if (options.returnOnFail && options.skipToTest !== test.name && !testRun.success && !testRun.skipped)
                                return res;
                        }
                        // res.push({ ...testRun, memoryDelta: (window?.performance as any)?.memory?.usedJSHeapSize - memoryUsageBefore, widgetsDifference: getWidgetsCountSafe() - widgetsBefore });
                    }
                }
                return res;
            });
        }
        function getWidgetsCountSafe() {
            var _a;
            if (typeof process !== 'undefined')
                return 0;
            let length = -1;
            try {
                length = DG.Widget.getAll().length;
            }
            catch (e) {
                console.warn((_a = e.message) !== null && _a !== void 0 ? _a : e);
            }
            return length;
        }
        function invokeTests(categoriesToInvoke, options) {
            var _a, _b, _c, _d;
            return __awaiter(this, void 0, void 0, function* () {
                try {
                    let skippingCategories = (options === null || options === void 0 ? void 0 : options.skipToCategory) != undefined;
                    let isTargetCategory = false;
                    for (const [key, value] of Object.entries(categoriesToInvoke)) {
                        if ((_a = options.exclude) === null || _a === void 0 ? void 0 : _a.some((c) => key.startsWith(c)))
                            continue;
                        if ((options === null || options === void 0 ? void 0 : options.category) != null && !key.toLowerCase().startsWith(options === null || options === void 0 ? void 0 : options.category.toLowerCase().trim()))
                            continue;
                        if (skippingCategories) {
                            if (isTargetCategory)
                                skippingCategories = false;
                            else {
                                if ((options === null || options === void 0 ? void 0 : options.skipToCategory) != null && key.toLowerCase().trim() === (options === null || options === void 0 ? void 0 : options.skipToCategory.toLowerCase().trim())) {
                                    isTargetCategory = true;
                                }
                                else {
                                    // Haven't found the target category yet, keep skipping
                                    continue;
                                }
                            }
                        }
                        let t = ((_b = value.tests) !== null && _b !== void 0 ? _b : []).filter((e) => matchesNodeTarget(e, value, options));
                        // Node/browser split: if no tests of this category belong to the requested target,
                        // skip it entirely — including before()/after().
                        if ((options.nodeOnly || options.excludeNodeTests) && t.length === 0)
                            continue;
                        const skipped = value.tests == undefined ? undefined : t.every((e) => {
                            var _a;
                            return ((_a = e.options) === null || _a === void 0 ? void 0 : _a.skipReason)
                                || ((options === null || options === void 0 ? void 0 : options.test) != null && options.test.toLowerCase() !== e.name.toLowerCase());
                        });
                        if (!skipped) {
                            const skippedCount = t.filter((e) => { var _a; return ((_a = e.options) === null || _a === void 0 ? void 0 : _a.skipReason) || ((options === null || options === void 0 ? void 0 : options.test) != null && options.test.toLowerCase() !== e.name.toLowerCase()); }).length;
                            stdLog(`Package testing: Started {{${key}}}${skippedCount > 0 ? ` skipped {{${skippedCount}}}` : ''}`);
                            value.beforeStatus = yield invokeCategoryMethod(value.before, key);
                        }
                        if (options.stressTest) {
                            t = t.filter((e) => { var _a; return (_a = e.options) === null || _a === void 0 ? void 0 : _a.stressTest; });
                            t = shuffle(t);
                        }
                        if (((_d = (_c = options.tags) === null || _c === void 0 ? void 0 : _c.length) !== null && _d !== void 0 ? _d : 0) > 0) {
                            t = t.filter((e) => { var _a, _b; return (_b = (_a = e.options) === null || _a === void 0 ? void 0 : _a.tags) === null || _b === void 0 ? void 0 : _b.some(tag => { var _a; return ((_a = options === null || options === void 0 ? void 0 : options.tags) !== null && _a !== void 0 ? _a : []).includes(tag); }); });
                        }
                        let res;
                        if (value.beforeStatus) {
                            res = Array.from(t.map((testElem) => {
                                return {
                                    date: new Date().toISOString(),
                                    category: key,
                                    name: testElem.name,
                                    success: false,
                                    result: 'before() failed',
                                    ms: 0,
                                    skipped: false,
                                    logs: '',
                                    owner: packageOwner,
                                    package: package_.name,
                                    widgetsDifference: 0,
                                    flaking: DG.Test.isReproducing
                                };
                            }));
                            res.forEach((test) => __awaiter(this, void 0, void 0, function* () { return yield grok.shell.reportTest('package', test); }));
                        }
                        else
                            res = yield invokeTestsInCategory(value, options, skippingCategories);
                        const data = res.filter((d) => d.result != 'skipped');
                        if (!skipped)
                            value.afterStatus = yield invokeCategoryMethod(value.after, key);
                        // Clear after category
                        // grok.shell.closeAll();
                        // DG.Balloon.closeAll();
                        if (value.afterStatus) {
                            stdLog(`Package testing: Category after() {{${key}}} failed`);
                            stdLog(`Package testing: Result for {{${key}}} after: ${value.afterStatus}`);
                            data.push({
                                date: new Date().toISOString(),
                                category: key,
                                name: 'after',
                                success: false,
                                result: value.afterStatus,
                                ms: 0,
                                skipped: false,
                                logs: '',
                                owner: packageOwner,
                                package: package_.name,
                                widgetsDifference: 0,
                                flaking: DG.Test.isReproducing
                            });
                        }
                        if (value.beforeStatus) {
                            stdLog(`Package testing: Category before() {{${key}}} failed`);
                            stdLog(`Package testing: Result for {{${key}}} before: ${value.beforeStatus}`);
                            data.push({
                                date: new Date().toISOString(),
                                category: key,
                                name: 'before',
                                success: false,
                                result: value.beforeStatus,
                                ms: 0,
                                skipped: false,
                                logs: '',
                                owner: packageOwner,
                                package: package_.name,
                                widgetsDifference: 0,
                                flaking: DG.Test.isReproducing
                            });
                        }
                        results.push(...data);
                        // If returnOnFail is set and a test failed (other than skipToTest), stop processing more categories
                        if (options.returnOnFail && data.some((d) => !d.success && !d.skipped && d.name !== options.skipToTest))
                            break;
                    }
                }
                finally {
                    resetConsole();
                }
                if (options.testContext.catchUnhandled && (!DG.Test.isInBenchmark)) {
                    yield delay(1000);
                    const error = yield grok.shell.lastError;
                    if (error != undefined) {
                        const params = {
                            logs: '',
                            date: new Date().toISOString(),
                            category: 'Unhandled exceptions',
                            name: 'Exception',
                            result: error !== null && error !== void 0 ? error : '',
                            success: !error,
                            ms: 0,
                            skipped: false,
                            owner: packageOwner !== null && packageOwner !== void 0 ? packageOwner : '',
                            'package': package_.name,
                            widgetsDifference: 0
                        };
                        stdLog(`Package testing: Unhandled Exception: ${error}`);
                        results.push(Object.assign(Object.assign({}, params), { 'flaking': DG.Test.isReproducing && !error }));
                        params.package = package_.name;
                        yield grok.shell.reportTest('package', params);
                    }
                }
            });
        }
    });
}
function getResult(x) {
    return __awaiter(this, void 0, void 0, function* () {
        return `${x.toString()}\n${x.stack ? (yield DG.Logger.translateStackTrace(x.stack)) : ''}`;
    });
}
function execTest(t, predicate, logs, testTimeout, packageName, verbose) {
    var _a, _b, _c, _d, _e, _f, _g, _h, _j, _k, _l, _m, _o, _p;
    return __awaiter(this, void 0, void 0, function* () {
        logs.length = 0;
        let r;
        let type = 'package';
        const filter = predicate != undefined && (t.name.toLowerCase() !== predicate.toLowerCase());
        let skip = ((_a = t.options) === null || _a === void 0 ? void 0 : _a.skipReason) || filter;
        let skipReason = filter ? 'skipped' : (_b = t.options) === null || _b === void 0 ? void 0 : _b.skipReason;
        if (DG.Test.isInBenchmark && !((_c = t.options) === null || _c === void 0 ? void 0 : _c.benchmark)) {
            stdLog(`Package testing: Skipped {{${t.category}}} {{${t.name}}} doesnt available in benchmark mode`);
            return undefined;
        }
        if (skip && !DG.Test.isInBenchmark)
            stdLog(`Package testing: Skipped {{${t.category}}} {{${t.name}}}`);
        if (!skip)
            stdLog(`Package testing: Started {{${t.category}}} {{${t.name}}}`);
        const start = Date.now();
        const startDate = new Date(start).toISOString();
        try {
            if (skip)
                r = { name: t.name, owner: (_e = (_d = t.options) === null || _d === void 0 ? void 0 : _d.owner) !== null && _e !== void 0 ? _e : '', category: t.category, logs: '', date: startDate, success: true, result: skipReason, ms: 0, skipped: true, package: packageName !== null && packageName !== void 0 ? packageName : '', flaking: DG.Test.isReproducing };
            else {
                let timeout_ = testTimeout !== null && testTimeout !== void 0 ? testTimeout : STANDART_TIMEOUT;
                if (DG.Test.isProfiling)
                    console.profile(`${t.category}: ${t.name}`);
                r = { name: t.name, owner: (_g = (_f = t.options) === null || _f === void 0 ? void 0 : _f.owner) !== null && _g !== void 0 ? _g : '', category: t.category, logs: '', date: startDate, success: true, result: (_h = (yield timeout(t.test, timeout_)).toString()) !== null && _h !== void 0 ? _h : 'OK', ms: 0, skipped: false, package: packageName !== null && packageName !== void 0 ? packageName : '', flaking: DG.Test.isReproducing };
                if (DG.Test.isProfiling) {
                    console.profileEnd(`${t.category}: ${t.name}`);
                    grok.shell.info(`Profiling of ${t.category}: ${t.name} finished \n Please ensure that you have opened DevTools (F12) / Performance panel before test starts.`);
                }
            }
        }
        catch (x) {
            stdError(x);
            r = { name: t.name, owner: (_k = (_j = t.options) === null || _j === void 0 ? void 0 : _j.owner) !== null && _k !== void 0 ? _k : '', category: t.category, logs: '', date: startDate, success: false, result: yield getResult(x), ms: 0, skipped: false, package: packageName !== null && packageName !== void 0 ? packageName : '', flaking: false };
        }
        if (((_l = t.options) === null || _l === void 0 ? void 0 : _l.isAggregated) && r.result.constructor === DG.DataFrame) {
            const col = r.result.col('success');
            if (col)
                r.success = col.stats.sum === col.length;
            if (!verbose) {
                const df = r.result;
                df.columns.remove('stack');
                df.rows.removeWhere((r) => r.get('success'));
                r.result = df;
            }
            r.result = r.result.toCsv();
        }
        r.logs = logs.join('\n');
        r.ms = Date.now() - start;
        if (!skip)
            stdLog(`Package testing: Finished {{${t.category}}} {{${t.name}}} with {{${r.success ? 'success' : 'error'}}} for ${r.ms} ms`);
        if (!r.success) {
            stdLog(`Package testing: Result for {{${t.category}}} {{${t.name}}}: ${r.result}`);
        }
        r.category = t.category;
        r.name = t.name;
        r.owner = (_o = (_m = t.options) === null || _m === void 0 ? void 0 : _m.owner) !== null && _o !== void 0 ? _o : '';
        if (!filter) {
            let params = {
                'success': r.success, 'result': r.result, 'ms': r.ms, 'date': r.date,
                'skipped': r.skipped, 'category': t.category, 'name': t.name, 'logs': r.logs, 'owner': r.owner,
                'flaking': DG.Test.isReproducing && r.success,
                'package': r.package
            };
            if (r.result.constructor == Object) {
                const res = Object.keys(r.result).reduce((acc, k) => (Object.assign(Object.assign({}, acc), { ['result.' + k]: r.result[k] })), {});
                params = Object.assign(Object.assign({}, params), res);
            }
            if (params.result instanceof DG.DataFrame)
                params.result = JSON.stringify((_p = params.result) === null || _p === void 0 ? void 0 : _p.toJson()) || '';
            yield grok.shell.reportTest(type, params);
        }
        return r;
    });
}
export function shuffle(array) {
    const newArr = array.slice();
    newArr.sort(() => Math.random() - 0.5);
    return newArr;
}
/* Waits [ms] milliseconds */
export function delay(ms) {
    return __awaiter(this, void 0, void 0, function* () {
        yield new Promise((r) => setTimeout(r, ms));
    });
}
export function awaitCheck(checkHandler, error = 'Timeout exceeded', wait = 500, interval = 50) {
    return __awaiter(this, void 0, void 0, function* () {
        return new Promise((resolve, reject) => {
            setTimeout(() => {
                clearInterval(intervalId);
                reject(new Error(error));
            }, wait);
            // @ts-ignore
            const intervalId = setInterval(() => {
                if (checkHandler()) {
                    clearInterval(intervalId);
                    resolve(null);
                }
            }, interval);
        });
    });
}
// Returns test execution result or an error in case of timeout
export function timeout(func, testTimeout, timeoutReason = 'EXECUTION TIMEOUT') {
    return __awaiter(this, void 0, void 0, function* () {
        let timeout = null;
        const timeoutPromise = new Promise((_, reject) => {
            timeout = setTimeout(() => {
                // eslint-disable-next-line prefer-promise-reject-errors
                reject(timeoutReason);
            }, testTimeout);
        });
        try {
            return yield Promise.race([func(), timeoutPromise]);
        }
        finally {
            if (timeout)
                clearTimeout(timeout);
        }
    });
}
export function isDialogPresent(dialogTitle) {
    const dialogs = DG.Dialog.getOpenDialogs();
    for (let i = 0; i < dialogs.length; i++) {
        if (dialogs[i].title == dialogTitle)
            return true;
    }
    return false;
}
/** Expects an asynchronous {@link action} to throw an exception. Use {@link check} to perform
 * deeper inspection of the exception if necessary.
 * @param  {function(): Promise<void>} action
 * @param  {function(any): boolean} check
 * @return {Promise<void>}
 */
export function expectExceptionAsync(action, check) {
    return __awaiter(this, void 0, void 0, function* () {
        let caught = false;
        let checked = false;
        try {
            yield action();
        }
        catch (e) {
            caught = true;
            checked = !check || check(e);
        }
        finally {
            if (!caught)
                throw new Error('An exception is expected but not thrown');
            if (!checked)
                throw new Error('An expected exception is thrown, but it does not satisfy the condition');
        }
    });
}
const catDF = DG.DataFrame.fromColumns([DG.Column.fromStrings('col', ['val1', 'val2', 'val3'])]);
/**
 * Universal test for viewers. It search viewers in DOM by tags: canvas, svg, img, input, h1, a
 * @param  {string} v Viewer name
 * @param  {_DG.DataFrame} df Dataframe to use. Should have at least 3 rows
 * @param  {boolean} options.detectSemanticTypes Specify whether to detect semantic types or not
 * @param  {boolean} options.readOnly If set to true, the dataframe will not be modified during the test
 * @param  {boolean} options.arbitraryDfTest If set to false, test on arbitrary dataframe
 * (one categorical column) will not be performed
 * @param  {object} options List of options (optional)
 * @return {Promise<void>} The test is considered successful if it completes without errors
 */
export function testViewer(v, df, options) {
    var _a;
    return __awaiter(this, void 0, void 0, function* () {
        const packageName = (_a = options === null || options === void 0 ? void 0 : options.packageName) !== null && _a !== void 0 ? _a : '';
        if (options === null || options === void 0 ? void 0 : options.detectSemanticTypes)
            yield grok.data.detectSemanticTypes(df);
        const tv = grok.shell.addTableView(df);
        try {
            //1. Open, do nothing and close
            yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded);
            //in case viewer with async rendering - wait for render to complete
            if (options === null || options === void 0 ? void 0 : options.awaitViewer)
                yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, undefined, options.awaitViewer);
            //2. Open viewer, run selection, filter, etc. and close
            if (!(options === null || options === void 0 ? void 0 : options.readOnly)) {
                yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, selectFilterChangeCurrent);
                if (options === null || options === void 0 ? void 0 : options.awaitViewer)
                    yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, selectFilterChangeCurrent, options.awaitViewer);
            }
            //2. Open viewer, change options, save layout and close
            let propsAndLayout = null;
            propsAndLayout = yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, changeOptionsSaveLayout);
            if (options === null || options === void 0 ? void 0 : options.awaitViewer)
                propsAndLayout = yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, changeOptionsSaveLayout, options.awaitViewer);
            //3. Load layout
            yield testViewerInternal(tv, v, packageName, grok.events.onViewLayoutApplied, loadLayout, undefined, propsAndLayout === null || propsAndLayout === void 0 ? void 0 : propsAndLayout.layout, { savedProps: propsAndLayout === null || propsAndLayout === void 0 ? void 0 : propsAndLayout.savedProps });
            if (options === null || options === void 0 ? void 0 : options.awaitViewer)
                yield testViewerInternal(tv, v, packageName, grok.events.onViewLayoutApplied, loadLayout, options.awaitViewer, propsAndLayout === null || propsAndLayout === void 0 ? void 0 : propsAndLayout.layout, { savedProps: propsAndLayout === null || propsAndLayout === void 0 ? void 0 : propsAndLayout.savedProps });
            //4. Open viewer on arbitary dataset
            if ((options === null || options === void 0 ? void 0 : options.arbitraryDfTest) !== false) {
                tv.dataFrame = catDF;
                yield delay(50);
                yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded);
                if (options === null || options === void 0 ? void 0 : options.awaitViewer)
                    yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, undefined, options.awaitViewer);
            }
            //5. Call postponed filtering
            yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, filterAsync);
            if (options === null || options === void 0 ? void 0 : options.awaitViewer)
                yield testViewerInternal(tv, v, packageName, grok.events.onViewerAdded, filterAsync, options.awaitViewer);
        }
        finally {
            // closeAll() is handling by common test workflow
            // grok.shell.closeAll();
            // DG.Balloon.closeAll();
        }
    });
}
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoidGVzdC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzIjpbInRlc3QudHMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6Ijs7Ozs7Ozs7O0FBS0EsT0FBTyxFQUFFLFFBQVEsRUFBRSxNQUFNLG1CQUFtQixDQUFDO0FBRTdDLE9BQU8sRUFBRSx1QkFBdUIsRUFBRSxXQUFXLEVBQUUsVUFBVSxFQUFFLHlCQUF5QixFQUFFLGtCQUFrQixFQUFFLE1BQU0scUJBQXFCLENBQUM7QUFFdEksTUFBTSxnQkFBZ0IsR0FBRyxLQUFLLENBQUM7QUFDL0IsTUFBTSxpQkFBaUIsR0FBRyxRQUFRLENBQUM7QUFFbkMsTUFBTSxNQUFNLEdBQUcsT0FBTyxDQUFDLEdBQUcsQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDekMsTUFBTSxPQUFPLEdBQUcsT0FBTyxDQUFDLElBQUksQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDM0MsTUFBTSxPQUFPLEdBQUcsT0FBTyxDQUFDLElBQUksQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDM0MsTUFBTSxRQUFRLEdBQUcsT0FBTyxDQUFDLEtBQUssQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFFN0MsTUFBTSxDQUFDLE1BQU0sS0FBSyxHQUVkLEVBQUUsQ0FBQztBQUVQLE1BQU0sZ0JBQWdCLEdBQUcsWUFBWSxDQUFDO0FBQ3RDLE1BQU0sV0FBVyxHQUFHLE1BQU0sQ0FBQztBQUMzQixNQUFNLGdCQUFnQixHQUFHLFdBQVcsQ0FBQztBQUNyQyxNQUFNLFdBQVcsR0FBRyxNQUFNLENBQUM7QUFDM0IsTUFBTSxhQUFhLEdBQStCLEVBQUUsQ0FBQztBQUNyRCxNQUFNLENBQUMsSUFBSSxlQUF1QixDQUFDO0FBRW5DLE1BQU0sS0FBVyxNQUFNLENBS3RCO0FBTEQsV0FBaUIsTUFBTTtJQUNyQixTQUFnQixPQUFPLENBQUMsS0FBVSxFQUFFLElBQWE7UUFDL0MsSUFBSSxLQUFLLElBQUksSUFBSTtZ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QVcsQ0FBQyxDQUFDO1lBQ3JGLElBQUksT0FBTyxhQUFQLE9BQU8sdUJBQVAsT0FBTyxDQUFFLFdBQVc7Z0JBQ3RCLE1BQU0sa0JBQWtCLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRSxXQUFXLEVBQUUsSUFBSSxDQUFDLE1BQU0sQ0FBQyxhQUFhLEVBQUUsV0FBVyxFQUFFLE9BQVEsQ0FBQyxXQUFXLENBQUMsQ0FBQztTQUU5RztnQkFBUztZQUNSLGlEQUFpRDtZQUNqRCx5QkFBeUI7WUFDekIseUJBQXlCO1NBQzFCOztDQUNGIiwic291cmNlc0NvbnRlbnQiOlsiaW1wb3J0IHR5cGUgKiBhcyBfZ3JvayBmcm9tICdkYXRhZ3Jvay1hcGkvZ3Jvayc7XHJcbmltcG9ydCB0eXBlICogYXMgX0RHIGZyb20gJ2RhdGFncm9rLWFwaS9kZyc7XHJcbmRlY2xhcmUgbGV0IGdyb2s6IHR5cGVvZiBfZ3JvaywgREc6IHR5cGVvZiBfREc7XHJcblxyXG5pbXBvcnQgeyBPYnNlcnZhYmxlIH0gZnJvbSAncnhqcyc7XHJcbmltcG9ydCB7IHRlc3REYXRhIH0gZnJvbSAnLi9kYXRhZnJhbWUtdXRpbHMnO1xyXG5pbXBvcnQgVGltZW91dCA9IE5vZGVKUy5UaW1lb3V0O1xyXG5pbXBvcnQgeyBjaGFuZ2VPcHRpb25zU2F2ZUxheW91dCwgZmlsdGVyQXN5bmMsIGxvYWRMYXlvdXQsIHNlbGVjdEZpbHRlckNoYW5nZUN1cnJlbnQsIHRlc3RWaWV3ZXJJbnRlcm5hbCB9IGZyb20gJy4vdGVzdC12aWV3ZXItdXRpbHMnO1xyXG5cclxuY29uc3QgU1RBTkRBUlRfVElNRU9VVCA9IDMwMDAwO1xyXG5jb25zdCBCRU5DSE1BUktfVElNRU9VVCA9IDEwODAwMDAwO1xyXG5cclxuY29uc3Qgc3RkTG9nID0gY29uc29sZS5sb2cuYmluZChjb25zb2xlKTtcclxuY29uc3Qgc3RkSW5mbyA9IGNvbnNvbGUuaW5mby5iaW5kKGNvbnNvbGUpO1xyXG5jb25zdCBzdGRXYXJuID0gY29uc29sZS53YXJuLmJpbmQoY29uc29sZSk7XHJcbmNvbnN0IHN0ZEVycm9yID0gY29uc29sZS5lcnJvci5iaW5kKGNvbnNvbGUpO1xyXG5cclxuZXhwb3J0IGNvbnN0IHRlc3RzOiB7XHJcbiAgW2tleTogc3RyaW5nXTogQ2F0ZWdvcnlcclxufSA9IHt9O1xyXG5cclxuY29uc3QgYXV0b1Rlc3RzQ2F0TmFtZSA9ICdBdXRvIFRlc3RzJztcclxuY29uc3QgZGVtb0NhdE5hbWUgPSAnRGVtbyc7XHJcbmNvbnN0IGRldGVjdG9yc0NhdE5hbWUgPSAnRGV0ZWN0b3JzJztcclxuY29uc3QgY29yZUNhdE5hbWUgPSAnQ29yZSc7XHJcbmNvbnN0IHdhc1JlZ2lzdGVyZWQ6IHsgW2tleTogc3RyaW5nXTogYm9vbGVhbiB9ID0ge307XHJcbmV4cG9ydCBsZXQgY3VycmVudENhdGVnb3J5OiBzdHJpbmc7XHJcblxyXG5leHBvcnQgbmFtZXNwYWNlIGFzc3VyZSB7XHJcbiAgZXhwb3J0IGZ1bmN0aW9uIG5vdE51bGwodmFsdWU6IGFueSwgbmFtZT86IHN0cmluZykge1xyXG4gICAgaWYgKHZhbHVlID09IG51bGwpXHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihgJHtuYW1lID09IG51bGwgPyAnVmFsdWUnIDogbmFtZX0gbm90IGRlZmluZWRgKTtcclxuICB9XHJcbn1cclxuXHJcbmV4cG9ydCBpbnRlcmZhY2UgVGVzdE9wdGlvbnMge1xyXG4gIHRpbWVvdXQ/OiBudW1iZXI7XHJcbiAgYmVuY2htYXJrV2FyblRpbWVvdXQ/OiBudW1iZXI7XHJcbiAgYmVuY2htYXJrVGltZW91dD86IG51bWJlcjtcclxuICB1bmhhbmRsZWRFeGNlcHRpb25UaW1lb3V0PzogbnVtYmVyO1xyXG4gIHNraXBSZWFzb24/OiBzdHJpbmc7XHJcbiAgaXNBZ2dyZWdhdGVkPzogYm9vbGVhbjtcclxuICBiZW5jaG1hcms/OiBib29sZWFuO1xyXG4gIHN0cmVzc1Rlc3Q/OiBib29sZWFuO1xyXG4gIG93bmVyPzogc3RyaW5nO1xyXG4gIHRhZ3M/OiBzdHJpbmdbXTtcclxuICAvKiogVGVzdCBuZWVkcyBubyBicm93c2VyOiBpdCB1c2VzIG9ubHkgTm9kZS1hdmFpbGFibGUgQVBJIGFuZCBjYW4gcnVuIGhlYWRsZXNzIHVuZGVyIHRoZSBqcy1hcGkgTm9kZSBydW50aW1lLiAqL1xyXG4gIG5vZGU/OiBib29sZWFuO1xyXG59XHJcblxyXG5leHBvcnQgaW50ZXJmYWNlIFRlc3RSZXN1bHQge1xyXG4gIGRhdGU6IHN0cmluZztcclxuICBjYXRlZ29yeTogc3RyaW5nO1xyXG4gIG5hbWU6IHN0cmluZztcclxuICBzdWNjZXNzOiBib29sZWFuO1xyXG4gIHJlc3VsdDogYW55O1xyXG4gIG1zOiBudW1iZXI7XHJcbiAgc2tpcHBlZDogYm9vbGVhbjtcclxuICBsb2dzOiBzdHJpbmc7XHJcbiAgb3duZXI6IHN0cmluZztcclxuICBwYWNrYWdlOiBzdHJpbmc7XHJcbiAgZmxha2luZzogYm9vbGVhbjtcclxufVxyXG5cclxuXHJcbmV4cG9ydCBpbnRlcmZhY2UgVGVzdFJlc3VsdEV4dGVuZGVkIGV4dGVuZHMgVGVzdFJlc3VsdHtcclxuICB3aWRnZXRzRGlmZmVyZW5jZTogbnVtYmVyO1xyXG59XHJcblxyXG5leHBvcnQgaW50ZXJmYWNlIENhdGVnb3J5T3B0aW9ucyB7XHJcbiAgY2xlYXI/OiBib29sZWFuO1xyXG4gIHRpbWVvdXQ/OiBudW1iZXI7XHJcbiAgYmVuY2htYXJrcz86IGJvb2xlYW47XHJcbiAgc3RyZXNzVGVzdHM/OiBib29sZWFuO1xyXG4gIG93bmVyPzogc3RyaW5nO1xyXG4gIC8qKiBEZWZhdWx0IGZvciBhbGwgdGVzdHMgaW4gdGhlIGNhdGVnb3J5IHRoYXQgZG9uJ3Qgc2V0IHRoZWlyIG93biBgbm9kZWAgb3B0aW9uLiAqL1xyXG4gIG5vZGU/OiBib29sZWFuO1xyXG59XHJcblxyXG5leHBvcnQgY2xhc3MgVGVzdENvbnRleHQge1xyXG4gIHN0cmVzc1Rlc3Q/OiBib29sZWFuO1xyXG4gIGNhdGNoVW5oYW5kbGVkID0gdHJ1ZTtcclxuICByZXBvcnQgPSBmYWxzZTtcclxuICByZXR1cm5PbkZhaWwgPSBmYWxzZTtcclxuICBjbG9zZUFsbCA9IHRydWU7XHJcblxyXG4gIGNvbnN0cnVjdG9yKGNhdGNoVW5oYW5kbGVkPzogYm9vbGVhbiwgcmVwb3J0PzogYm9vbGVhbiwgcmV0dXJuT25GYWlsPzogYm9vbGVhbikge1xyXG4gICAgaWYgKGNhdGNoVW5oYW5kbGVkICE9PSB1bmRlZmluZWQpIHRoaXMuY2F0Y2hVbmhhbmRsZWQgPSBjYXRjaFVuaGFuZGxlZDtcclxuICAgIGlmIChyZXBvcnQgIT09IHVuZGVmaW5lZCkgdGhpcy5yZXBvcnQgPSByZXBvcnQ7XHJcbiAgICBpZiAocmV0dXJuT25GYWlsICE9PSB1bmRlZmluZWQpIHRoaXMucmV0dXJuT25GYWlsID0gcmV0dXJuT25GYWlsO1xyXG4gIH07XHJcbn1cclxuXHJcbmV4cG9ydCBjbGFzcyBUZXN0IHtcclxuICB0ZXN0OiAoKSA9PiBQcm9taXNlPGFueT47XHJcbiAgbmFtZTogc3RyaW5nO1xyXG4gIGNhdGVnb3J5OiBzdHJpbmc7XHJcbiAgb3B0aW9ucz86IFRlc3RPcHRpb25zO1xyXG5cclxuICBjb25zdHJ1Y3RvcihjYXRlZ29yeTogc3RyaW5nLCBuYW1lOiBzdHJpbmcsIHRlc3Q6ICgpID0+IFByb21pc2U8YW55Piwgb3B0aW9ucz86IFRlc3RPcHRpb25zKSB7XHJcbiAgICB0aGlzLmNhdGVnb3J5ID0gY2F0ZWdvcnk7XHJcbiAgICB0aGlzLm5hbWUgPSBuYW1lO1xyXG4gICAgb3B0aW9ucyA/Pz0ge307XHJcbiAgICBvcHRpb25zLnRpbWVvdXQgPz89IFNUQU5EQVJUX1RJTUVPVVQ7XHJcbiAgICB0aGlzLm9wdGlvbnMgPSBvcHRpb25zO1xyXG4gICAgdGhpcy50ZXN0ID0gYXN5bmMgKCk6IFByb21pc2U8YW55PiA9PiB7XHJcbiAgICAgIHJldHVybiBuZXcgUHJvbWlzZShhc3luYyAocmVzb2x2ZSwgcmVqZWN0KSA9PiB7XHJcbiAgICAgICAgbGV0IHJlc3VsdCA9ICcnO1xyXG4gICAgICAgIHRyeSB7XHJcbiAgICAgICAgICBpZiAoREcuVGVzdC5pc0luRGVidWcpXHJcbiAgICAgICAgICAgIGRlYnVnZ2VyO1xyXG5cclxuICAgICAgICAgIGxldCByZXMgPSBhd2FpdCB0ZXN0KCk7XHJcbiAgICAgICAgICB0cnkge1xyXG4gICAgICAgICAgICByZXN1bHQgPSByZXM/LnRvU3RyaW5nKCkgPz8gJyc7XHJcbiAgICAgICAgICB9XHJcbiAgICAgICAgICBjYXRjaCAoZSkge1xyXG4gICAgICAgICAgICByZXN1bHQgPSAnQ2FuXFwndCBjb252ZXJ0IHRlc3RcXCdzIHJlc3VsdCB0byBzdHJpbmcnO1xyXG4gICAgICAgICAgICBjb25zb2xlLmVycm9yKGBDYW5cXCd0IGNvbnZlcnQgdGVzdFxcJ3MgcmVzdWx0IHRvIHN0cmluZyBpbiB0aGUgJHt0aGlzLmNhdGVnb3J5fToke3RoaXMubmFtZX0gdGVzdGApO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgIH0gY2F0Y2ggKGU6IGFueSkge1xyXG4gICAgICAgICAgcmVqZWN0KGUpO1xyXG4gICAgICAgIH1cclxuICAgICAgICByZXNvbHZlKHJlc3VsdCk7XHJcbiAgICAgIH0pO1xyXG4gICAgfTtcclxuICB9XHJcbn1cclxuXHJcbmV4cG9ydCBjbGFzcyBDYXRlZ29yeSB7XHJcbiAgdGVzdHM/OiBUZXN0W107XHJcbiAgYmVmb3JlPzogKCkgPT4gUHJvbWlzZTx2b2lkPjtcclxuICBhZnRlcj86ICgpID0+IFByb21pc2U8dm9pZD47XHJcblxyXG4gIGJlZm9yZVN0YXR1cz86IHN0cmluZztcclxuICBhZnRlclN0YXR1cz86IHN0cmluZztcclxuICBjbGVhcj86IGJvb2xlYW47XHJcbiAgdGltZW91dD86IG51bWJlcjtcclxuICBiZW5jaG1hcmtzPzogYm9vbGVhbjtcclxuICBiZW5jaG1hcmtUaW1lb3V0PzogbnVtYmVyO1xyXG4gIHN0cmVzc1Rlc3RzPzogYm9vbGVhbjtcclxuICBvd25lcj86IHN0cmluZztcclxuICBub2RlPzogYm9vbGVhbjtcclxufVxyXG5cclxuZXhwb3J0IGNsYXNzIE5vZGVUZXN0RXhlY3V0aW9uT3B0aW9ucyB7XHJcbiAgcGFja2FnZSE6IF9ERy5QYWNrYWdlO1xyXG59XHJcblxyXG5leHBvcnQgY2xhc3MgVGVzdEV4ZWN1dGlvbk9wdGlvbnMge1xyXG4gIGNhdGVnb3J5Pzogc3RyaW5nO1xyXG4gIHRlc3Q/OiBzdHJpbmc7XHJcbiAgdGVzdENvbnRleHQ/OiBUZXN0Q29udGV4dDtcclxuICBleGNsdWRlPzogc3RyaW5nW107XHJcbiAgdmVyYm9zZT86IGJvb2xlYW47XHJcbiAgc3RyZXNzVGVzdD86IGJvb2xlYW47XHJcbiAgdGFncz86IHN0cmluZ1tdO1xyXG4gIG5vZGVPcHRpb25zPzogTm9kZVRlc3RFeGVjdXRpb25PcHRpb25zO1xyXG4gIHNraXBUb0NhdGVnb3J5Pzogc3RyaW5nO1xyXG4gIHNraXBUb1Rlc3Q/OiBzdHJpbmc7XHJcbiAgcmV0dXJuT25GYWlsPzogYm9vbGVhbjtcclxuICAvKiogUnVuIG9ubHkgdGVzdHMgbWFya2VkIGBub2RlOiB0cnVlYCAodGhlIGhlYWRsZXNzIE5vZGUgcGFzcyBvZiBgZ3JvayB0ZXN0YCkuICovXHJcbiAgbm9kZU9ubHk/OiBib29sZWFuO1xyXG4gIC8qKiBTa2lwIHRlc3RzIG1hcmtlZCBgbm9kZTogdHJ1ZWAgKHRoZSBicm93c2VyIHBhc3Mgb2YgYGdyb2sgdGVzdGAgYWZ0ZXIgYSBOb2RlIHBhc3MgYWxyZWFkeSByYW4gdGhlbSkuICovXHJcbiAgZXhjbHVkZU5vZGVUZXN0cz86IGJvb2xlYW47XHJcbn1cclxuXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiB0ZXN0RXZlbnQ8VD4oZXZlbnQ6IE9ic2VydmFibGU8VD4sXHJcbiAgaGFuZGxlcjogKGFyZ3M6IFQpID0+IHZvaWQsIHRyaWdnZXI6ICgpID0+IHZvaWQsIG1zOiBudW1iZXIgPSAwLCByZWFzb246IHN0cmluZyA9IGB0aW1lb3V0YFxyXG4pOiBQcm9taXNlPGFueT4ge1xyXG4gIHJldHVybiBuZXcgUHJvbWlzZSgocmVzb2x2ZSwgcmVqZWN0KSA9PiB7XHJcbiAgICBjb25zdCBzdWIgPSBldmVudC5zdWJzY3JpYmUoKGFyZ3M6IFQpID0+IHtcclxuICAgICAgdHJ5IHtcclxuICAgICAgICBoYW5kbGVyKGFyZ3MpO1xyXG4gICAgICAgIHJlc29sdmUoJ09LJyk7XHJcbiAgICAgIH0gY2F0Y2ggKGUpIHtcclxuICAgICAgICByZWplY3QoZSk7XHJcbiAgICAgIH0gZmluYWxseSB7XHJcbiAgICAgICAgc3ViLnVuc3Vic2NyaWJlKCk7XHJcbiAgICAgICAgY2xlYXJUaW1lb3V0KHRpbWVvdXQpO1xyXG4gICAgICB9XHJcbiAgICB9KTtcclxuICAgIGNvbnN0IHRpbWVvdXQgPSBzZXRUaW1lb3V0KCgpID0+IHtcclxuICAgICAgc3ViLnVuc3Vic2NyaWJlKCk7XHJcbiAgICAgIC8vIGVzbGludC1kaXNhYmxlLW5leHQtbGluZSBwcmVmZXItcHJvbWlzZS1yZWplY3QtZXJyb3JzXHJcbiAgICAgIHJlamVjdChyZWFzb24pO1xyXG4gICAgfSwgbXMpO1xyXG4gICAgdHJpZ2dlcigpO1xyXG4gIH0pO1xyXG59XHJcblxyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gdGVzdEV2ZW50QXN5bmM8VD4oZXZlbnQ6IE9ic2VydmFibGU8VD4sXHJcbiAgaGFuZGxlcjogKGFyZ3M6IFQpID0+IFByb21pc2U8dm9pZD4sIHRyaWdnZXI6ICgpID0+IHZvaWQsIG1zOiBudW1iZXIgPSAwLCByZWFzb246IHN0cmluZyA9IGB0aW1lb3V0YFxyXG4pOiBQcm9taXNlPGFueT4ge1xyXG4gIHJldHVybiBuZXcgUHJvbWlzZSgocmVzb2x2ZSwgcmVqZWN0KSA9PiB7XHJcbiAgICBjb25zdCBzdWIgPSBldmVudC5zdWJzY3JpYmUoKGFyZ3M6IFQpID0+IHtcclxuICAgICAgaGFuZGxlcihhcmdzKS50aGVuKCgpID0+IHtcclxuICAgICAgICByZXNvbHZlKCdPSycpO1xyXG4gICAgICB9KS5jYXRjaCgoZSkgPT4ge1xyXG4gICAgICAgIHJlamVjdChlKTtcclxuICAgICAgfSkuZmluYWxseSgoKSA9PiB7XHJcbiAgICAgICAgc3ViLnVuc3Vic2NyaWJlKCk7XHJcbiAgICAgICAgY2xlYXJUaW1lb3V0KHRpbWVvdXQpO1xyXG4gICAgICB9KTtcclxuICAgIH0pO1xyXG4gICAgY29uc3QgdGltZW91dCA9IHNldFRpbWVvdXQoKCkgPT4ge1xyXG4gICAgICBzdWIudW5zdWJzY3JpYmUoKTtcclxuICAgICAgLy8gZXNsaW50LWRpc2FibGUtbmV4dC1saW5lIHByZWZlci1wcm9taXNlLXJlamVjdC1lcnJvcnNcclxuICAgICAgcmVqZWN0KHJlYXNvbik7XHJcbiAgICB9LCBtcyk7XHJcbiAgICB0cmlnZ2VyKCk7XHJcbiAgfSk7XHJcbn1cclxuXHJcbmV4cG9ydCBmdW5jdGlvbiB0ZXN0KG5hbWU6IHN0cmluZywgdGVzdDogKCkgPT4gUHJvbWlzZTxhbnk+LCBvcHRpb25zPzogVGVzdE9wdGlvbnMpOiB2b2lkIHtcclxuICBpZiAodGVzdHNbY3VycmVudENhdGVnb3J5XSA9PSB1bmRlZmluZWQpXHJcbiAgICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldID0ge307XHJcbiAgaWYgKHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0udGVzdHMgPT0gdW5kZWZpbmVkKVxyXG4gICAgdGVzdHNbY3VycmVudENhdGVnb3J5XS50ZXN0cyA9IFtdO1xyXG4gIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0udGVzdHMhLnB1c2gobmV3IFRlc3QoY3VycmVudENhdGVnb3J5LCBuYW1lLCB0ZXN0LCBvcHRpb25zKSk7XHJcbn1cclxuXHJcbi8qIFRlc3RzIHR3byBvYmplY3RzIGZvciBlcXVhbGl0eSwgdGhyb3dzIGFuIGV4Y2VwdGlvbiBpZiB0aGV5IGFyZSBub3QgZXF1YWwuICovXHJcbmV4cG9ydCBmdW5jdGlvbiBleHBlY3QoYWN0dWFsOiBhbnksIGV4cGVjdGVkOiBhbnkgPSB0cnVlLCBlcnJvcj86IHN0cmluZyk6IHZvaWQge1xyXG4gIGlmIChlcnJvcilcclxuICAgIGVycm9yID0gYCR7ZXJyb3J9LCBgO1xyXG4gIGVsc2UgZXJyb3IgPSAnJztcclxuICBpZiAoYWN0dWFsICE9PSBleHBlY3RlZClcclxuICAgIHRocm93IG5ldyBFcnJvcihgJHtlcnJvcn1FeHBlY3RlZCBcIiR7ZXhwZWN0ZWR9XCIsIGdvdCBcIiR7YWN0dWFsfVwiYCk7XHJcbn1cclxuXHJcbmV4cG9ydCBmdW5jdGlvbiBleHBlY3RGbG9hdChhY3R1YWw6IG51bWJlciwgZXhwZWN0ZWQ6IG51bWJlciwgdG9sZXJhbmNlID0gMC4wMDEsIGVycm9yPzogc3RyaW5nKTogdm9pZCB7XHJcbiAgaWYgKChhY3R1YWwgPT09IE51bWJlci5QT1NJVElWRV9JTkZJTklUWSAmJiBleHBlY3RlZCA9PT0gTnVtYmVyLlBPU0lUSVZFX0lORklOSVRZKSB8fFxyXG4gICAgKGFjdHVhbCA9PT0gTnVtYmVyLk5FR0FUSVZFX0lORklOSVRZICYmIGV4cGVjdGVkID09PSBOdW1iZXIuTkVHQVRJVkVfSU5GSU5JVFkpIHx8XHJcbiAgICAoYWN0dWFsID09PSBOdW1iZXIuTmFOICYmIGV4cGVjdGVkID09PSBOdW1iZXIuTmFOKSB8fCAoaXNOYU4oYWN0dWFsKSAmJiBpc05hTihleHBlY3RlZCkpKVxyXG4gICAgcmV0dXJuO1xyXG4gIGNvbnN0IGFyZUVxdWFsID0gTWF0aC5hYnMoYWN0dWFsIC0gZXhwZWN0ZWQpIDwgdG9sZXJhbmNlO1xyXG4gIGV4cGVjdChhcmVFcXVhbCwgdHJ1ZSwgYCR7ZXJyb3IgPz8gJyd9ICh0b2xlcmFuY2UgPSAke3RvbGVyYW5jZX07IGEgPSAke2FjdHVhbH0sIGUgPSAke2V4cGVjdGVkfSlgKTtcclxuICBpZiAoIWFyZUVxdWFsKVxyXG4gICAgdGhyb3cgbmV3IEVycm9yKGBFeHBlY3RlZCAke2V4cGVjdGVkfSwgZ290ICR7YWN0dWFsfSAodG9sZXJhbmNlID0gJHt0b2xlcmFuY2V9KWApO1xyXG59XHJcblxyXG5leHBvcnQgZnVuY3Rpb24gZXhwZWN0VGFibGUoYWN0dWFsOiBfREcuRGF0YUZyYW1lLCBleHBlY3RlZDogX0RHLkRhdGFGcmFtZSwgZXJyb3I/OiBzdHJpbmcpOiB2b2lkIHtcclxuICBjb25zdCBleHBlY3RlZFJvd0NvdW50ID0gZXhwZWN0ZWQucm93Q291bnQ7XHJcbiAgY29uc3QgYWN0dWFsUm93Q291bnQgPSBhY3R1YWwucm93Q291bnQ7XHJcbiAgZXhwZWN0KGFjdHVhbFJvd0NvdW50LCBleHBlY3RlZFJvd0NvdW50LCBgJHtlcnJvciA/PyAnJ30sIHJvdyBjb3VudGApO1xyXG5cclxuICBmb3IgKGNvbnN0IGNvbHVtbiBvZiBleHBlY3RlZC5jb2x1bW5zKSB7XHJcbiAgICBjb25zdCBhY3R1YWxDb2x1bW4gPSBhY3R1YWwuY29sdW1ucy5ieU5hbWUoY29sdW1uLm5hbWUpO1xyXG4gICAgaWYgKGFjdHVhbENvbHVtbiA9PSBudWxsKVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYENvbHVtbiAke2NvbHVtbi5uYW1lfSBub3QgZm91bmRgKTtcclxuICAgIGlmIChhY3R1YWxDb2x1bW4udHlwZSAhPSBjb2x1bW4udHlwZSlcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKGBDb2x1bW4gJHtjb2x1bW4ubmFtZX0gdHlwZSBleHBlY3RlZCAke2NvbHVtbi50eXBlfSBnb3QgJHthY3R1YWxDb2x1bW4udHlwZX1gKTtcclxuICAgIGZvciAobGV0IGkgPSAwOyBpIDwgZXhwZWN0ZWRSb3dDb3VudDsgaSsrKSB7XHJcbiAgICAgIGNvbnN0IHZhbHVlID0gY29sdW1uLmdldChpKTtcclxuICAgICAgY29uc3QgYWN0dWFsVmFsdWUgPSBhY3R1YWxDb2x1bW4uZ2V0KGkpO1xyXG4gICAgICBpZiAoY29sdW1uLnR5cGUgPT0gREcuVFlQRS5GTE9BVClcclxuICAgICAgICBleHBlY3RGbG9hdChhY3R1YWxWYWx1ZSwgdmFsdWUsIDAuMDAwMSwgZXJyb3IpO1xyXG4gICAgICBlbHNlIGlmIChjb2x1bW4udHlwZSA9PSBERy5UWVBFLkRBVEVfVElNRSlcclxuICAgICAgICBleHBlY3QoYWN0dWFsVmFsdWUuaXNTYW1lKHZhbHVlKSwgdHJ1ZSwgZXJyb3IpO1xyXG4gICAgICBlbHNlXHJcbiAgICAgICAgZXhwZWN0KGFjdHVhbFZhbHVlLCB2YWx1ZSwgZXJyb3IpO1xyXG4gICAgfVxyXG4gIH1cclxufVxyXG5cclxuZXhwb3J0IGZ1bmN0aW9uIGV4cGVjdE9iamVjdChhY3R1YWw6IHsgW2tleTogc3RyaW5nXTogYW55IH0sIGV4cGVjdGVkOiB7IFtrZXk6IHN0cmluZ106IGFueSB9KSB7XHJcbiAgZm9yIChjb25zdCBbZXhwZWN0ZWRLZXksIGV4cGVjdGVkVmFsdWVdIG9mIE9iamVjdC5lbnRyaWVzKGV4cGVjdGVkKSkge1xyXG4gICAgaWYgKCFhY3R1YWwuaGFzT3duUHJvcGVydHkoZXhwZWN0ZWRLZXkpKVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYEV4cGVjdGVkIHByb3BlcnR5IFwiJHtleHBlY3RlZEtleX1cIiBub3QgZm91bmRgKTtcclxuXHJcbiAgICBjb25zdCBhY3R1YWxWYWx1ZSA9IGFjdHVhbFtleHBlY3RlZEtleV07XHJcbiAgICBpZiAoYWN0dWFsVmFsdWUgaW5zdGFuY2VvZiBBcnJheSAmJiBleHBlY3RlZFZhbHVlIGluc3RhbmNlb2YgQXJyYXkpXHJcbiAgICAgIGV4cGVjdEFycmF5KGFjdHVhbFZhbHVlLCBleHBlY3RlZFZhbHVlKTtcclxuICAgIGVsc2UgaWYgKGFjdHVhbFZhbHVlIGluc3RhbmNlb2YgT2JqZWN0ICYmIGV4cGVjdGVkVmFsdWUgaW5zdGFuY2VvZiBPYmplY3QpXHJcbiAgICAgIGV4cGVjdE9iamVjdChhY3R1YWxWYWx1ZSwgZXhwZWN0ZWRWYWx1ZSk7XHJcbiAgICBlbHNlIGlmIChOdW1iZXIuaXNGaW5pdGUoYWN0dWFsVmFsdWUpICYmIE51bWJlci5pc0Zpbml0ZShleHBlY3RlZFZhbHVlKSlcclxuICAgICAgZXhwZWN0RmxvYXQoYWN0dWFsVmFsdWUsIGV4cGVjdGVkVmFsdWUpO1xyXG4gICAgZWxzZSBpZiAoYWN0dWFsVmFsdWUgIT0gZXhwZWN0ZWRWYWx1ZSlcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKGBFeHBlY3RlZCAoJHtleHBlY3RlZFZhbHVlfSkgZm9yIGtleSAnJHtleHBlY3RlZEtleX0nLCBnb3QgKCR7YWN0dWFsVmFsdWV9KWApO1xyXG4gIH1cclxufVxyXG5cclxuZXhwb3J0IGZ1bmN0aW9uIGV4cGVjdEFycmF5KGFjdHVhbDogQXJyYXlMaWtlPGFueT4sIGV4cGVjdGVkOiBBcnJheUxpa2U8YW55Pikge1xyXG4gIGNvbnN0IGFjdHVhbExlbmd0aCA9IGFjdHVhbC5sZW5ndGg7XHJcbiAgY29uc3QgZXhwZWN0ZWRMZW5ndGggPSBleHBlY3RlZC5sZW5ndGg7XHJcblxyXG4gIGlmIChhY3R1YWxMZW5ndGggIT0gZXhwZWN0ZWRMZW5ndGgpIHtcclxuICAgIHRocm93IG5ldyBFcnJvcihgQXJyYXlzIGFyZSBvZiBkaWZmZXJlbnQgbGVuZ3RoOiBhY3R1YWwgYXJyYXkgbGVuZ3RoIGlzICR7YWN0dWFsTGVuZ3RofSBgICtcclxuICAgICAgYGFuZCBleHBlY3RlZCBhcnJheSBsZW5ndGggaXMgJHtleHBlY3RlZExlbmd0aH1gKTtcclxuICB9XHJcblxyXG4gIGZvciAobGV0IGkgPSAwOyBpIDwgYWN0dWFsTGVuZ3RoOyBpKyspIHtcclxuICAgIGlmIChhY3R1YWxbaV0gaW5zdGFuY2VvZiBBcnJheSAmJiBleHBlY3RlZFtpXSBpbnN0YW5jZW9mIEFycmF5KVxyXG4gICAgICBleHBlY3RBcnJheShhY3R1YWxbaV0sIGV4cGVjdGVkW2ldKTtcclxuICAgIGVsc2UgaWYgKGFjdHVhbFtpXSBpbnN0YW5jZW9mIE9iamVjdCAmJiBleHBlY3RlZFtpXSBpbnN0YW5jZW9mIE9iamVjdClcclxuICAgICAgZXhwZWN0T2JqZWN0KGFjdHVhbFtpXSwgZXhwZWN0ZWRbaV0pO1xyXG4gICAgZWxzZSBpZiAoYWN0dWFsW2ldICE9IGV4cGVjdGVkW2ldKVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYEV4cGVjdGVkICR7ZXhwZWN0ZWRbaV19IGF0IHBvc2l0aW9uICR7aX0sIGdvdCAke2FjdHVhbFtpXX1gKTtcclxuICB9XHJcbn1cclxuXHJcbi8qIERlZmluZXMgYSB0ZXN0IHN1aXRlLiAqL1xyXG5leHBvcnQgZnVuY3Rpb24gY2F0ZWdvcnkoY2F0ZWdvcnk6IHN0cmluZywgdGVzdHNfOiAoKSA9PiB2b2lkLCBvcHRpb25zPzogQ2F0ZWdvcnlPcHRpb25zKTogdm9pZCB7XHJcbiAgY3VycmVudENhdGVnb3J5ID0gY2F0ZWdvcnk7XHJcbiAgdGVzdHNfKCk7XHJcbiAgaWYgKHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0pIHtcclxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0uY2xlYXIgPSBvcHRpb25zPy5jbGVhciA/PyB0cnVlO1xyXG4gICAgdGVzdHNbY3VycmVudENhdGVnb3J5XS50aW1lb3V0ID0gb3B0aW9ucz8udGltZW91dDtcclxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0uYmVuY2htYXJrcyA9IG9wdGlvbnM/LmJlbmNobWFya3M7XHJcbiAgICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldLnN0cmVzc1Rlc3RzID0gb3B0aW9ucz8uc3RyZXNzVGVzdHM7XHJcbiAgICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldLm93bmVyID0gb3B0aW9ucz8ub3duZXI7XHJcbiAgICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldLm5vZGUgPSBvcHRpb25zPy5ub2RlO1xyXG4gIH1cclxufVxyXG5cclxuLyogRGVmaW5lcyBhIGZ1bmN0aW9uIHRvIGJlIGV4ZWN1dGVkIGJlZm9yZSB0aGUgdGVzdHMgaW4gdGhpcyBjYXRlZ29yeSBhcmUgZXhlY3V0ZWQuICovXHJcbmV4cG9ydCBmdW5jdGlvbiBiZWZvcmUoYmVmb3JlOiAoKSA9PiBQcm9taXNlPHZvaWQ+KTogdm9pZCB7XHJcbiAgaWYgKHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0gPT0gdW5kZWZpbmVkKVxyXG4gICAgdGVzdHNbY3VycmVudENhdGVnb3J5XSA9IHt9O1xyXG4gIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0uYmVmb3JlID0gYmVmb3JlO1xyXG59XHJcblxyXG4vKiBEZWZpbmVzIGEgZnVuY3Rpb24gdG8gYmUgZXhlY3V0ZWQgYWZ0ZXIgdGhlIHRlc3RzIGluIHRoaXMgY2F0ZWdvcnkgYXJlIGV4ZWN1dGVkLiAqL1xyXG5leHBvcnQgZnVuY3Rpb24gYWZ0ZXIoYWZ0ZXI6ICgpID0+IFByb21pc2U8dm9pZD4pOiB2b2lkIHtcclxuICBpZiAodGVzdHNbY3VycmVudENhdGVnb3J5XSA9PSB1bmRlZmluZWQpXHJcbiAgICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldID0ge307XHJcbiAgdGVzdHNbY3VycmVudENhdGVnb3J5XS5hZnRlciA9IGFmdGVyO1xyXG59XHJcblxyXG5mdW5jdGlvbiBhZGROYW1lc3BhY2Uoczogc3RyaW5nLCBmOiBfREcuRnVuYyk6IHN0cmluZyB7XHJcbiAgcmV0dXJuIHMucmVwbGFjZShuZXcgUmVnRXhwKGYubmFtZSwgJ2dpJyksIGYubnFOYW1lKTtcclxufVxyXG5cclxuLyoqIFdoZXRoZXIgYSB0ZXN0IG1hdGNoZXMgdGhlIG5vZGUvYnJvd3NlciBzcGxpdCByZXF1ZXN0ZWQgYnkgdGhlIHJ1biBvcHRpb25zLiAqL1xyXG5mdW5jdGlvbiBtYXRjaGVzTm9kZVRhcmdldCh0OiBUZXN0LCBjYXQ6IENhdGVnb3J5LCBvcHRpb25zOiBUZXN0RXhlY3V0aW9uT3B0aW9ucyk6IGJvb2xlYW4ge1xyXG4gIGNvbnN0IGlzTm9kZSA9IHQub3B0aW9ucz8ubm9kZSA/PyBjYXQubm9kZSA/PyBmYWxzZTtcclxuICByZXR1cm4gb3B0aW9ucy5ub2RlT25seSA/IGlzTm9kZSA6IG9wdGlvbnMuZXhjbHVkZU5vZGVUZXN0cyA/ICFpc05vZGUgOiB0cnVlO1xyXG59XHJcblxyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gaW5pdEF1dG9UZXN0cyhwYWNrYWdlXzogX0RHLlBhY2thZ2UsIG1vZHVsZT86IGFueSkge1xyXG4gIGNvbnN0IHBhY2thZ2VJZCA9IHBhY2thZ2VfLmlkO1xyXG4gIGlmICh3YXNSZWdpc3RlcmVkW3BhY2thZ2VJZF0pIHJldHVybjtcclxuICBjb25zdCBtb2R1bGVUZXN0cyA9IG1vZHVsZSA/IG1vZHVsZS50ZXN0cyA6IHRlc3RzO1xyXG4gIGlmIChwYWNrYWdlXy5uYW1lID09PSAnRGV2VG9vbHMnIHx8ICghIW1vZHVsZSAmJiBtb2R1bGUuX3BhY2thZ2UubmFtZSA9PT0gJ0RldlRvb2xzJykpIHtcclxuICAgIGZvciAoY29uc3QgZiBvZiAoPGFueT53aW5kb3cpLmRhcnRUZXN0cykge1xyXG4gICAgICBjb25zdCBhcnIgPSBmLm5hbWUuc3BsaXQoL1xccypcXHxcXHMqIS9nKTtcclxuICAgICAgbGV0IG5hbWUgPSBhcnIucG9wKCkgPz8gZi5uYW1lO1xyXG4gICAgICBsZXQgY2F0ID0gYXJyLmxlbmd0aCA/IGNvcmVDYXROYW1lICsgJzogJyArIGFyci5qb2luKCc6ICcpIDogY29yZUNhdE5hbWU7XHJcbiAgICAgIGxldCBmdWxsTmFtZTogc3RyaW5nW10gPSBuYW1lLnNwbGl0KCcgfCAnKTtcclxuICAgICAgbmFtZSA9IGZ1bGxOYW1lW2Z1bGxOYW1lLmxlbmd0aCAtIDFdO1xyXG4gICAgICBmdWxsTmFtZS51bnNoaWZ0KGNhdCk7XHJcbiAgICAgIGZ1bGxOYW1lLnBvcCgpO1xyXG4gICAgICBjYXQgPSBmdWxsTmFtZS5qb2luKCc6ICcpO1xyXG4gICAgICBpZiAobW9kdWxlVGVzdHNbY2F0XSA9PT0gdW5kZWZpbmVkKVxyXG4gICAgICAgIG1vZHVsZVRlc3RzW2NhdF0gPSB7IHRlc3RzOiBbXSwgY2xlYXI6IHRydWUgfTtcclxuICAgICAgbW9kdWxlVGVzdHNbY2F0XS50ZXN0cy5wdXNoKG5ldyBUZXN0KGNhdCwgbmFtZSwgZi50ZXN0LCB7IGlzQWdncmVnYXRlZDogZmFsc2UsIHRpbWVvdXQ6IGYub3B0aW9ucz8udGltZW91dCA/PyBTVEFOREFSVF9USU1FT1VULCBza2lwUmVhc29uOiBmLm9wdGlvbnM/LnNraXBSZWFzb24sIG93bmVyOiBmLm9wdGlvbnM/Lm93bmVyLCBiZW5jaG1hcms6IGYub3B0aW9ucz8uYmVuY2htYXJrID8/IGZhbHNlIH0pKTtcclxuICAgIH1cclxuICB9XHJcbiAgY29uc3QgbW9kdWxlQXV0b1Rlc3RzID0gW107XHJcbiAgY29uc3QgbW9kdWxlRGVtbyA9IFtdO1xyXG4gIGNvbnN0IG1vZHVsZURldGVjdG9ycyA9IFtdO1xyXG4gIGNvbnN0IHBhY2tGdW5jdGlvbnMgPSBhd2FpdCBncm9rLmRhcGkuZnVuY3Rpb25zLmZpbHRlcihgcGFja2FnZS5pZCA9IFwiJHtwYWNrYWdlSWR9XCJgKS5saXN0KCk7XHJcbiAgY29uc3QgcmVnID0gbmV3IFJlZ0V4cCgvc2tpcDpcXHMqKFteLFxcc10rKXx3YWl0OlxccyooXFxkKyl8Y2F0OlxccyooW14sXFxzXSspfHRpbWVvdXQ6XFxzKihcXGQrKS9nKTtcclxuICBmb3IgKGNvbnN0IGYgb2YgcGFja0Z1bmN0aW9ucykge1xyXG4gICAgY29uc3QgdGVzdHMgPSBmLm9wdGlvbnNbJ3Rlc3QnXTtcclxuICAgIGNvbnN0IGRlbW8gPSBmLm9wdGlvbnNbJ2RlbW9QYXRoJ107XHJcbiAgICBpZiAoKHRlc3RzICYmIEFycmF5LmlzQXJyYXkodGVzdHMpICYmIHRlc3RzLmxlbmd0aCkpIHtcclxuICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCB0ZXN0cy5sZW5ndGg7IGkrKykge1xyXG4gICAgICAgIGNvbnN0IHJlcyA9ICh0ZXN0c1tpXSBhcyBzdHJpbmcpLm1hdGNoQWxsKHJlZyk7XHJcbiAgICAgICAgY29uc3QgbWFwOiB7IHNraXA/OiBzdHJpbmcsIHdhaXQ/OiBudW1iZXIsIGNhdD86IHN0cmluZywgdGltZW91dD86IG51bWJlciwgYmVuY2htYXJrVGltZW91dD86IG51bWJlciB9ID0ge307XHJcbiAgICAgICAgQXJyYXkuZnJvbShyZXMpLmZvckVhY2goKGFycikgPT4ge1xyXG4gICAgICAgICAgaWYgKGFyclswXS5zdGFydHNXaXRoKCdza2lwJykpIG1hcFsnc2tpcCddID0gYXJyWzFdO1xyXG4gICAgICAgICAgZWxzZSBpZiAoYXJyWzBdLnN0YXJ0c1dpdGgoJ3dhaXQnKSkgbWFwWyd3YWl0J10gPSBwYXJzZUludChhcnJbMl0pO1xyXG4gICAgICAgICAgZWxzZSBpZiAoYXJyWzBdLnN0YXJ0c1dpdGgoJ2NhdCcpKSBtYXBbJ2NhdCddID0gYXJyWzNdO1xyXG4gICAgICAgICAgZWxzZSBpZiAoYXJyWzBdLnN0YXJ0c1dpdGgoJ3RpbWVvdXQnKSkgbWFwWyd0aW1lb3V0J10gPSBwYXJzZUludChhcnJbNF0pO1xyXG4gICAgICAgIH0pO1xyXG4gICAgICAgIGNvbnN0IHRlc3QgPSBuZXcgVGVzdChtYXAuY2F0ID8/IGF1dG9UZXN0c0NhdE5hbWUsIHRlc3RzLmxlbmd0aCA9PT0gMSA/IGYubmFtZSA6IGAke2YubmFtZX0gJHtpICsgMX1gLCBhc3luYyAoKSA9PiB7XHJcbiAgICAgICAgICBjb25zdCByZXMgPSBhd2FpdCBncm9rLmZ1bmN0aW9ucy5ldmFsKGFkZE5hbWVzcGFjZSh0ZXN0c1tpXSwgZikpO1xyXG4gICAgICAgICAgaWYgKG1hcC53YWl0KSBhd2FpdCBkZWxheShtYXAud2FpdCk7XHJcbiAgICAgICAgICAvLyBlc2xpbnQtZGlzYWJsZS1uZXh0LWxpbmUgbm8tdGhyb3ctbGl0ZXJhbFxyXG4gICAgICAgICAgaWYgKHR5cGVvZiByZXMgPT09ICdib29sZWFuJyAmJiAhcmVzKSB0aHJvdyBgRmFpbGVkOiAke3Rlc3RzW2ldfSwgZXhwZWN0ZWQgdHJ1ZSwgZ290ICR7cmVzfWA7XHJcbiAgICAgICAgfSwgeyBza2lwUmVhc29uOiBtYXAuc2tpcCwgdGltZW91dDogREcuVGVzdC5pc0luQmVuY2htYXJrID8gbWFwLmJlbmNobWFya1RpbWVvdXQgPz8gQkVOQ0hNQVJLX1RJTUVPVVQgOiBtYXAudGltZW91dCA/PyBTVEFOREFSVF9USU1FT1VULFxyXG4gICAgICAgICAgLy8gUXVlcnkgdGVzdHMgZXZhbHVhdGUgc2VydmVyLXNpZGUgKGV4cHJlc3Npb25zIGxpa2UgRXhwZWN0VGFibGUoUXVlcnkoKSwgT3BlbkZpbGUoLi4uKSlcclxuICAgICAgICAgIC8vIHJlc29sdmUgdGhyb3VnaCBjb3JlIGZ1bmN0aW9ucyBhdmFpbGFibGUgaW4gdGhlIE5vZGUgcnVudGltZSksIHNvIHRoZXkgY2FuIHJ1biBoZWFkbGVzcztcclxuICAgICAgICAgIC8vIHRlc3RzIG9mIGNsaWVudCBKUyBmdW5jdGlvbnMgbmVlZCB0aGUgYnJvd3Nlci1sb2FkZWQgcGFja2FnZS5cclxuICAgICAgICAgIG5vZGU6IGYgaW5zdGFuY2VvZiBERy5EYXRhUXVlcnkgfSk7XHJcbiAgICAgICAgaWYgKG1hcC5jYXQpIHtcclxuICAgICAgICAgIGNvbnN0IGNhdDogc3RyaW5nID0gbWFwLmNhdDtcclxuICAgICAgICAgIGlmIChtb2R1bGVUZXN0c1tjYXRdID09PSB1bmRlZmluZWQpXHJcbiAgICAgICAgICAgIG1vZHVsZVRlc3RzW2NhdF0gPSB7IHRlc3RzOiBbXSwgY2xlYXI6IHRydWUgfTtcclxuXHJcbiAgICAgICAgICAvLyBvbmx5IGJlZm9yZS9hZnRlciBjYW4gYmUgZGVmaW5lZCBpbiB0cyBmaWxlcyB0ZXN0cyB1bmRlciB0aGUgY2F0ZWdvcnlcclxuICAgICAgICAgIGlmICghbW9kdWxlVGVzdHNbY2F0XS50ZXN0cylcclxuICAgICAgICAgICAgbW9kdWxlVGVzdHNbY2F0XS50ZXN0cyA9IFtdO1xyXG4gICAgICAgICAgbW9kdWxlVGVzdHNbY2F0XS50ZXN0cy5wdXNoKHRlc3QpO1xyXG4gICAgICAgIH1cclxuICAgICAgICBlbHNlXHJcbiAgICAgICAgICBtb2R1bGVBdXRvVGVzdHMucHVzaCh0ZXN0KTtcclxuICAgICAgfVxyXG4gICAgfVxyXG4gICAgaWYgKGRlbW8pIHtcclxuICAgICAgY29uc3Qgd2FpdCA9IGYub3B0aW9uc1snZGVtb1dhaXQnXSA/IHBhcnNlSW50KGYub3B0aW9uc1snZGVtb1dhaXQnXSkgOiB1bmRlZmluZWQ7XHJcbiAgICAgIGNvbnN0IHRlc3QgPSBuZXcgVGVzdChkZW1vQ2F0TmFtZSwgZi5mcmllbmRseU5hbWUsIGFzeW5jICgpID0+IHtcclxuICAgICAgICBhd2FpdCBkZWxheSgzMDApO1xyXG4gICAgICAgIGdyb2suc2hlbGwuY2xlYXJMYXN0RXJyb3IoKTtcclxuICAgICAgICBhd2FpdCBmLmFwcGx5KCk7XHJcbiAgICAgICAgYXdhaXQgZGVsYXkod2FpdCA/IHdhaXQgOiAyMDAwKTtcclxuICAgICAgICBjb25zdCB1bmhhbmRsZWQgPSBhd2FpdCBncm9rLnNoZWxsLmxhc3RFcnJvcjtcclxuICAgICAgICBpZiAodW5oYW5kbGVkKVxyXG4gICAgICAgICAgdGhyb3cgbmV3IEVycm9yKHVuaGFuZGxlZCk7XHJcbiAgICAgIH0sIHsgc2tpcFJlYXNvbjogZi5vcHRpb25zWydkZW1vU2tpcCddIH0pO1xyXG4gICAgICBtb2R1bGVEZW1vLnB1c2godGVzdCk7XHJcbiAgICB9XHJcbiAgICBpZiAoZi5oYXNUYWcoJ3NlbVR5cGVEZXRlY3RvcicpKSB7XHJcbiAgICAgIGxldCBkZXRlY3RvcnNUZXN0RGF0YSA9IHRlc3REYXRhO1xyXG4gICAgICBpZiAoZi5vcHRpb25zWyd0ZXN0RGF0YSddKSB7XHJcbiAgICAgICAgZGV0ZWN0b3JzVGVzdERhdGEgPSBhd2FpdCBncm9rLmRhdGEuZmlsZXMub3BlblRhYmxlKGBTeXN0ZW06QXBwRGF0YS8ke3BhY2thZ2VfLm5xTmFtZX0vJHtmLm9wdGlvbnNbJ3Rlc3REYXRhJ119YCk7XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIGNvbnN0IHRlc3QgPSBuZXcgVGVzdChkZXRlY3RvcnNDYXROYW1lLCBmLmZyaWVuZGx5TmFtZSwgYXN5bmMgKCkgPT4ge1xyXG4gICAgICAgIGNvbnN0IGFyciA9IFtdO1xyXG4gICAgICAgIGNvbnNvbGUubG9nKGBTeXN0ZW06QXBwRGF0YS8ke3BhY2thZ2VfLm5xTmFtZX0vJHtmLm9wdGlvbnNbJ3Rlc3REYXRhJ119YCk7XHJcblxyXG4gICAgICAgIGZvciAoY29uc3QgY29sIG9mIGRldGVjdG9yc1Rlc3REYXRhLmNsb25lKCkuY29sdW1ucykge1xyXG4gICAgICAgICAgY29uc3QgcmVzID0gYXdhaXQgZi5hcHBseShbY29sXSk7XHJcbiAgICAgICAgICBhcnIucHVzaChyZXMgfHwgY29sLnNlbVR5cGUpO1xyXG4gICAgICAgIH1cclxuICAgICAgICBjb25zdCByZXNBcnIgPSBhcnIuZmlsdGVyKChpKSA9PiBpKTtcclxuICAgICAgICBleHBlY3QocmVzQXJyLmxlbmd0aCwgMSk7XHJcblxyXG4gICAgICAgIGlmIChmLm9wdGlvbnNbJ3Rlc3REYXRhQ29sdW1uTmFtZSddKVxyXG4gICAgICAgICAgZXhwZWN0KHJlc0FyclswXSwgZi5vcHRpb25zWyd0ZXN0RGF0YUNvbHVtbk5hbWUnXSk7XHJcblxyXG4gICAgICB9LCB7IHNraXBSZWFzb246IGYub3B0aW9uc1snc2tpcFRlc3QnXSB9KTtcclxuICAgICAgbW9kdWxlRGV0ZWN0b3JzLnB1c2godGVzdCk7XHJcbiAgICB9XHJcbiAgfVxyXG4gIHdhc1JlZ2lzdGVyZWRbcGFja2FnZUlkXSA9IHRydWU7XHJcbiAgaWYgKG1vZHVsZUF1dG9UZXN0cy5sZW5ndGggPiAwKVxyXG4gICAgbW9kdWxlVGVzdHNbYXV0b1Rlc3RzQ2F0TmFtZV0gPSB7IHRlc3RzOiBtb2R1bGVBdXRvVGVzdHMsIGNsZWFyOiB0cnVlIH07XHJcbiAgaWYgKG1vZHVsZURlbW8ubGVuZ3RoID4gMClcclxuICAgIG1vZHVsZVRlc3RzW2RlbW9DYXROYW1lXSA9IHsgdGVzdHM6IG1vZHVsZURlbW8sIGNsZWFyOiB0cnVlIH07XHJcbiAgaWYgKG1vZHVsZURldGVjdG9ycy5sZW5ndGggPiAwKVxyXG4gICAgbW9kdWxlVGVzdHNbZGV0ZWN0b3JzQ2F0TmFtZV0gPSB7IHRlc3RzOiBtb2R1bGVEZXRlY3RvcnMsIGNsZWFyOiBmYWxzZSB9O1xyXG59XHJcblxyXG5mdW5jdGlvbiByZWRlZmluZUNvbnNvbGUoKTogYW55W10ge1xyXG4gIGNvbnN0IGxvZ3M6IGFueVtdID0gW107XHJcbiAgY29uc29sZS5sb2cgPSAoLi4uYXJncykgPT4ge1xyXG4gICAgbG9ncy5wdXNoKC4uLmFyZ3MpO1xyXG4gICAgc3RkTG9nKC4uLmFyZ3MpO1xyXG4gIH07XHJcbiAgY29uc29sZS5pbmZvID0gKC4uLmFyZ3MpID0+IHtcclxuICAgIGxvZ3MucHVzaCguLi5hcmdzKTtcclxuICAgIHN0ZEluZm8oLi4uYXJncyk7XHJcbiAgfTtcclxuICBjb25zb2xlLndhcm4gPSAoLi4uYXJncykgPT4ge1xyXG4gICAgbG9ncy5wdXNoKC4uLmFyZ3MpO1xyXG4gICAgc3RkV2FybiguLi5hcmdzKTtcclxuICB9O1xyXG4gIGNvbnNvbGUuZXJyb3IgPSAoLi4uYXJncykgPT4ge1xyXG4gICAgbG9ncy5wdXNoKC4uLmFyZ3MpO1xyXG4gICAgc3RkRXJyb3IoLi4uYXJncyk7XHJcbiAgfTtcclxuICByZXR1cm4gbG9ncztcclxufVxyXG5cclxuZnVuY3Rpb24gcmVzZXRDb25zb2xlKCk6IHZvaWQge1xyXG4gIGNvbnNvbGUubG9nID0gc3RkTG9nO1xyXG4gIGNvbnNvbGUuaW5mbyA9IHN0ZEluZm87XHJcbiAgY29uc29sZS53YXJuID0gc3RkV2FybjtcclxuICBjb25zb2xlLmVycm9yID0gc3RkRXJyb3I7XHJcbn1cclxuXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBydW5UZXN0cyhvcHRpb25zPzogVGVzdEV4ZWN1dGlvbk9wdGlvbnMpIDogUHJvbWlzZTxUZXN0UmVzdWx0RXh0ZW5kZWRbXT57XHJcblxyXG4gIGNvbnN0IHBhY2thZ2VfOiBfREcuUGFja2FnZSA9IG9wdGlvbnM/Lm5vZGVPcHRpb25zID8gb3B0aW9ucy5ub2RlT3B0aW9ucy5wYWNrYWdlIDogZ3Jvay5mdW5jdGlvbnMuZ2V0Q3VycmVudENhbGwoKS5mdW5jLnBhY2thZ2U7XHJcbiAgaWYgKCFwYWNrYWdlXylcclxuICAgIHRocm93IG5ldyBFcnJvcignQ2FuXFwndCBydW4gdGVzdHMgb3V0c2lkZSBvZiB0aGUgcGFja2FnZScpO1xyXG4gIGNvbnN0IG1hdGNoID0gcGFja2FnZV8ucGFja2FnZU93bmVyPy5tYXRjaCgvPChbXj5dKik+Lyk7XHJcbiAgY29uc3QgcGFja2FnZU93bmVyID0gbWF0Y2ggPyBtYXRjaFsxXSA6ICcnO1xyXG4gIGlmIChwYWNrYWdlXyAhPSB1bmRlZmluZWQpXHJcbiAgICBhd2FpdCBpbml0QXV0b1Rlc3RzKHBhY2thZ2VfKTtcclxuICBjb25zdCByZXN1bHRzOlRlc3RSZXN1bHRFeHRlbmRlZFtdID0gW107XHJcbiAgY29uc29sZS5sb2coYFJ1bm5pbmcgdGVzdHMuLi5gKTtcclxuICBjb25zb2xlLmxvZyhvcHRpb25zKTtcclxuICBvcHRpb25zID8/PSB7fTtcclxuICBvcHRpb25zIS50ZXN0Q29udGV4dCA/Pz0gbmV3IFRlc3RDb250ZXh0KCk7XHJcbiAgZ3Jvay5zaGVsbC5jbGVhckxhc3RFcnJvcigpO1xyXG4gIGNvbnN0IGxvZ3MgPSByZWRlZmluZUNvbnNvbGUoKTtcclxuXHJcbiAgYXdhaXQgaW52b2tlVGVzdHModGVzdHMsIG9wdGlvbnMpO1xyXG5cclxuICBmb3IgKGxldCByIG9mIHJlc3VsdHMpIHtcclxuICAgIHIucmVzdWx0ID0gci5yZXN1bHQudG9TdHJpbmcoKS5yZXBsYWNlKC9cIi9nLCAnXFwnJyk7XHJcbiAgICBpZiAoci5sb2dzICE9IHVuZGVmaW5lZClcclxuICAgICAgci5sb2dzID0gci5sb2dzIS50b1N0cmluZygpLnJlcGxhY2UoL1wiL2csICdcXCcnKTtcclxuICB9XHJcbiAgcmV0dXJuIHJlc3VsdHM7XHJcblxyXG4gIGFzeW5jIGZ1bmN0aW9uIGludm9rZUNhdGVnb3J5TWV0aG9kKG1ldGhvZDogKCgpID0+IFByb21pc2U8dm9pZD4pIHwgdW5kZWZpbmVkLCBjYXRlZ29yeTogc3RyaW5nKTogUHJvbWlzZTxzdHJpbmcgfCB1bmRlZmluZWQ+IHtcclxuICAgIGxldCBpbnZva2F0aW9uUmVzdWx0ID0gdW5kZWZpbmVkO1xyXG4gICAgdHJ5IHtcclxuICAgICAgaWYgKG1ldGhvZCAhPT0gdW5kZWZpbmVkKSB7XHJcbiAgICAgICAgYXdhaXQgdGltZW91dChhc3luYyAoKSA9PiB7XHJcbiAgICAgICAgICBhd2FpdCBtZXRob2QoKTtcclxuICAgICAgICB9LCAxMDAwMDAsIGBiZWZvcmUgJHtjYXRlZ29yeX06IHRpbWVvdXQgZXJyb3JgKTtcclxuICAgICAgfVxyXG4gICAgfSBjYXRjaCAoeDogYW55KSB7XHJcbiAgICAgIGludm9rYXRpb25SZXN1bHQgPSBhd2FpdCBnZXRSZXN1bHQoeCk7XHJcbiAgICB9XHJcbiAgICByZXR1cm4gaW52b2thdGlvblJlc3VsdFxyXG4gIH1cclxuXHJcbiAgYXN5bmMgZnVuY3Rpb24gaW52b2tlVGVzdHNJbkNhdGVnb3J5KGNhdGVnb3J5OiBDYXRlZ29yeSwgb3B0aW9uczogVGVzdEV4ZWN1dGlvbk9wdGlvbnMsIGlzVGFyZ2V0Q2F0ZWdvcnk6IGJvb2xlYW4pOiBQcm9taXNlPFRlc3RSZXN1bHRFeHRlbmRlZFtdPiB7XHJcbiAgICBsZXQgdCA9IChjYXRlZ29yeS50ZXN0cyA/PyBbXSkuZmlsdGVyKChlKSA9PiBtYXRjaGVzTm9kZVRhcmdldChlLCBjYXRlZ29yeSwgb3B0aW9ucykpO1xyXG4gICAgY29uc3QgcmVzIDogVGVzdFJlc3VsdEV4dGVuZGVkW10gPSBbXTtcclxuICAgIC8vIGxldCBtZW1vcnlVc2FnZUJlZm9yZSA9ICh3aW5kb3c/LnBlcmZvcm1hbmNlIGFzIGFueSk/Lm1lbW9yeT8udXNlZEpTSGVhcFNpemU7XHJcbiAgICBjb25zdCB3aWRnZXRzQmVmb3JlID0gZ2V0V2lkZ2V0c0NvdW50U2FmZSgpO1xyXG5cclxuICAgIGlmIChjYXRlZ29yeS5jbGVhcikge1xyXG4gICAgICAgIGxldCBza2lwcGluZ1Rlc3RzID0gaXNUYXJnZXRDYXRlZ29yeSAmJiBvcHRpb25zLnNraXBUb1Rlc3QgIT0gdW5kZWZpbmVkO1xyXG4gICAgICBmb3IgKGxldCBpID0gMDsgaSA8IHQubGVuZ3RoOyBpKyspIHtcclxuXHJcbiAgICAgICAgaWYgKHRbaV0ub3B0aW9ucykge1xyXG4gICAgICAgICAgaWYgKHRbaV0ub3B0aW9ucz8uYmVuY2htYXJrID09PSB1bmRlZmluZWQpIHtcclxuICAgICAgICAgICAgaWYgKCF0W2ldLm9wdGlvbnMpXHJcbiAgICAgICAgICAgICAgdFtpXS5vcHRpb25zID0ge31cclxuICAgICAgICAgICAgdFtpXS5vcHRpb25zIS5iZW5jaG1hcmsgPSBjYXRlZ29yeS5iZW5jaG1hcmtzID8/IGZhbHNlO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgIH1cclxuICAgICAgICBsZXQgdGVzdCA9IHRbaV07XHJcbiAgICAgICAgaWYgKG9wdGlvbnMudGVzdClcclxuICAgICAgICAgIGlmIChvcHRpb25zLnRlc3QudG9Mb3dlckNhc2UoKSAhPT0gdGVzdC5uYW1lLnRvTG93ZXJDYXNlKCkpXHJcbiAgICAgICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICAgIGlmIChza2lwcGluZ1Rlc3RzKSB7XHJcbiAgICAgICAgICBpZiAob3B0aW9ucz8uc2tpcFRvVGVzdCAhPSB1bmRlZmluZWQgJiYgdGVzdC5uYW1lLnRvTG93ZXJDYXNlKCkudHJpbSgpID09PSBvcHRpb25zPy5za2lwVG9UZXN0LnRvTG93ZXJDYXNlKCkudHJpbSgpKSB7XHJcbiAgICAgICAgICAgIC8vIEZvdW5kIHRoZSB0YXJnZXQgdGVzdCwgc3RvcCBza2lwcGluZyBhZnRlciB0aGlzIG9uZVxyXG4gICAgICAgICAgICBza2lwcGluZ1Rlc3RzID0gZmFsc2U7XHJcbiAgICAgICAgICB9IGVsc2VcclxuICAgICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICAgIH1cclxuICAgICAgICBpZiAodGVzdD8ub3B0aW9ucykge1xyXG4gICAgICAgICAgdGVzdC5vcHRpb25zLm93bmVyID0gdFtpXS5vcHRpb25zPy5vd25lciA/PyBjYXRlZ29yeT8ub3duZXIgPz8gcGFja2FnZU93bmVyID8/ICcnO1xyXG4gICAgICAgIH1cclxuICAgICAgICAvLyBsZXQgaXNHQkVuYWJsZSA9ICh3aW5kb3cgYXMgYW55KS5nYyAmJiB0ZXN0Lm9wdGlvbnM/LnNraXBSZWFzb24gPT0gdW5kZWZpbmVkO1xyXG4gICAgICAgIC8vIGNvbnNvbGUubG9nKGAqKioqKioqKiR7aXNHQkVuYWJsZX1gKTtcclxuICAgICAgICAvLyBpZiAoaXNHQkVuYWJsZSlcclxuICAgICAgICAvLyAgIGF3YWl0ICh3aW5kb3cgYXMgYW55KS5nYygpO1xyXG4gICAgICAgIC8vIG1lbW9yeVVzYWdlQmVmb3JlID0gKHdpbmRvdz8ucGVyZm9ybWFuY2UgYXMgYW55KT8ubWVtb3J5Py51c2VkSlNIZWFwU2l6ZTtcclxuICAgICAgICBsZXQgdGVzdFJ1biA9IGF3YWl0IGV4ZWNUZXN0KFxyXG4gICAgICAgICAgICB0ZXN0LFxyXG4gICAgICAgICAgICBvcHRpb25zPy50ZXN0LFxyXG4gICAgICAgICAgICBsb2dzLCBERy5UZXN0LmlzSW5CZW5jaG1hcmsgPyB0W2ldLm9wdGlvbnM/LmJlbmNobWFya1RpbWVvdXQgPz8gQkVOQ0hNQVJLX1RJTUVPVVQgOiB0W2ldLm9wdGlvbnM/LnRpbWVvdXQgPz8gU1RBTkRBUlRfVElNRU9VVCxcclxuICAgICAgICAgICAgcGFja2FnZV8ubmFtZSxcclxuICAgICAgICAgICAgb3B0aW9ucy52ZXJib3NlXHJcbiAgICAgICAgKTtcclxuXHJcbiAgICAgICAgLy8gaWYgKGlzR0JFbmFibGUpXHJcbiAgICAgICAgLy8gICBhd2FpdCAod2luZG93IGFzIGFueSkuZ2MoKTtcclxuICAgICAgICBpZiAodGVzdFJ1bikge1xyXG4gICAgICAgICAgcmVzLnB1c2goeyAuLi50ZXN0UnVuLCAgd2lkZ2V0c0RpZmZlcmVuY2U6IGdldFdpZGdldHNDb3VudFNhZmUoKSAtIHdpZGdldHNCZWZvcmUgfSk7XHJcbiAgICAgICAgICAvLyBSZXR1cm4gZWFybHkgaWYgcmV0dXJuT25GYWlsIGlzIHNldCBhbmQgdGVzdCBmYWlsZWQgKGJ1dCBpZ25vcmUgZmFpbHVyZSBmb3IgdGhlIHNraXBUb1Rlc3QgdGVzdCBpdHNlbGYpXHJcbiAgICAgICAgICBpZiAob3B0aW9ucy5yZXR1cm5PbkZhaWwgJiYgb3B0aW9ucy5za2lwVG9UZXN0ICE9PSB0ZXN0Lm5hbWUgJiYgIXRlc3RSdW4uc3VjY2VzcyAmJiAhdGVzdFJ1bi5za2lwcGVkKVxyXG4gICAgICAgICAgICByZXR1cm4gcmVzO1xyXG4gICAgICAgIH1cclxuICAgICAgICAvLyByZXMucHVzaCh7IC4uLnRlc3RSdW4sIG1lbW9yeURlbHRhOiAod2luZG93Py5wZXJmb3JtYW5jZSBhcyBhbnkpPy5tZW1vcnk/LnVzZWRKU0hlYXBTaXplIC0gbWVtb3J5VXNhZ2VCZWZvcmUsIHdpZGdldHNEZWx0YTogZ2V0V2lkZ2V0c0NvdW50U2FmZSgpIC0gd2lkZ2V0c0JlZm9yZSB9KTtcclxuXHJcbiAgICAgICAgaWYgKCFvcHRpb25zLm5vZGVPcHRpb25zICYmIG9wdGlvbnMudGVzdENvbnRleHQ/LmNsb3NlQWxsICE9PSBmYWxzZSkge1xyXG4gICAgICAgICAgZ3Jvay5zaGVsbC5jbG9zZUFsbCgpO1xyXG4gICAgICAgICAgREcuQmFsbG9vbi5jbG9zZUFsbCgpO1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG4gICAgfSBlbHNlIHtcclxuICAgICAgbGV0IHNraXBwaW5nVGVzdHMgPSBpc1RhcmdldENhdGVnb3J5ICYmIG9wdGlvbnMuc2tpcFRvVGVzdCAhPSB1bmRlZmluZWQ7XHJcbiAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgdC5sZW5ndGg7IGkrKykge1xyXG4gICAgICAgIGxldCB0ZXN0ID0gdFtpXTtcclxuICAgICAgICBpZiAob3B0aW9ucy50ZXN0KVxyXG4gICAgICAgICAgaWYgKG9wdGlvbnMudGVzdC50b0xvd2VyQ2FzZSgpICE9PSB0ZXN0Lm5hbWUudG9Mb3dlckNhc2UoKSlcclxuICAgICAgICAgICAgY29udGludWU7XHJcbiAgICAgICAgaWYgKHNraXBwaW5nVGVzdHMpIHtcclxuICAgICAgICAgIGlmIChvcHRpb25zPy5za2lwVG9UZXN0ICE9IHVuZGVmaW5lZCAmJiB0ZXN0Lm5hbWUudG9Mb3dlckNhc2UoKS50cmltKCkgPT09IG9wdGlvbnM/LnNraXBUb1Rlc3QudG9Mb3dlckNhc2UoKS50cmltKCkpIHtcclxuICAgICAgICAgICAgLy8gRm91bmQgdGhlIHRhcmdldCB0ZXN0LCBzdG9wIHNraXBwaW5nIGFmdGVyIHRoaXMgb25lXHJcbiAgICAgICAgICAgIHNraXBwaW5nVGVzdHMgPSBmYWxzZTtcclxuICAgICAgICAgIH1cclxuICAgICAgICAgIGNvbnRpbnVlOyAgLy8gU2tpcCB0aGlzIHRlc3QgKGluY2x1ZGluZyB0aGUgdGFyZ2V0KVxyXG4gICAgICAgIH1cclxuXHJcbiAgICAgICAgaWYgKHRlc3Q/Lm9wdGlvbnMpIHtcclxuICAgICAgICAgIHRlc3Qub3B0aW9ucy5vd25lciA9IHRbaV0ub3B0aW9ucz8ub3duZXIgPz8gY2F0ZWdvcnk/Lm93bmVyID8/IHBhY2thZ2VPd25lciA/PyAnJztcclxuICAgICAgICB9XHJcbiAgICAgICAgLy8gbGV0IGlzR0JFbmFibGUgPSAod2luZG93IGFzIGFueSkuZ2MgJiYgdGVzdC5vcHRpb25zPy5za2lwUmVhc29uID09IHVuZGVmaW5lZDtcclxuICAgICAgICAvLyBjb25zb2xlLmxvZyhgKioqKioqKioke2lzR0JFbmFibGV9YCk7XHJcbiAgICAgICAgLy8gaWYgKGlzR0JFbmFibGUpXHJcbiAgICAgICAgLy8gICBhd2FpdCAod2luZG93IGFzIGFueSkuZ2MoKTtcclxuICAgICAgICAvLyBtZW1vcnlVc2FnZUJlZm9yZSA9ICh3aW5kb3c/LnBlcmZvcm1hbmNlIGFzIGFueSk/Lm1lbW9yeT8udXNlZEpTSGVhcFNpemU7XHJcbiAgICAgICAgbGV0IHRlc3RSdW4gPSBhd2FpdCBleGVjVGVzdChcclxuICAgICAgICAgICAgdGVzdCxcclxuICAgICAgICAgICAgb3B0aW9ucz8udGVzdCxcclxuICAgICAgICAgICAgbG9ncyxcclxuICAgICAgICAgICAgREcuVGVzdC5pc0luQmVuY2htYXJrID8gdFtpXS5vcHRpb25zPy5iZW5jaG1hcmtUaW1lb3V0ID8/IEJFTkNITUFSS19USU1FT1VUIDogdFtpXS5vcHRpb25zPy50aW1lb3V0LFxyXG4gICAgICAgICAgICBwYWNrYWdlXy5uYW1lLFxyXG4gICAgICAgICAgICBvcHRpb25zLnZlcmJvc2VcclxuICAgICAgICApO1xyXG5cclxuICAgICAgICAvLyBpZiAoaXNHQkVuYWJsZSlcclxuICAgICAgICAvLyAgIGF3YWl0ICh3aW5kb3cgYXMgYW55KS5nYygpO1xyXG5cclxuICAgICAgICBpZiAodGVzdFJ1bikge1xyXG4gICAgICAgICAgcmVzLnB1c2goeyAuLi50ZXN0UnVuLCB3aWRnZXRzRGlmZmVyZW5jZTogZ2V0V2lkZ2V0c0NvdW50U2FmZSgpIC0gd2lkZ2V0c0JlZm9yZSB9KTtcclxuICAgICAgICAgIC8vIFJldHVybiBlYXJseSBpZiByZXR1cm5PbkZhaWwgaXMgc2V0IGFuZCB0ZXN0IGZhaWxlZCAoYnV0IGlnbm9yZSBmYWlsdXJlIGZvciB0aGUgc2tpcFRvVGVzdCB0ZXN0IGl0c2VsZilcclxuICAgICAgICAgIGlmIChvcHRpb25zLnJldHVybk9uRmFpbCAmJiBvcHRpb25zLnNraXBUb1Rlc3QgIT09IHRlc3QubmFtZSAmJiAhdGVzdFJ1bi5zdWNjZXNzICYmICF0ZXN0UnVuLnNraXBwZWQpXHJcbiAgICAgICAgICAgIHJldHVybiByZXM7XHJcbiAgICAgICAgfVxyXG4gICAgICAgIC8vIHJlcy5wdXNoKHsgLi4udGVzdFJ1biwgbWVtb3J5RGVsdGE6ICh3aW5kb3c/LnBlcmZvcm1hbmNlIGFzIGFueSk/Lm1lbW9yeT8udXNlZEpTSGVhcFNpemUgLSBtZW1vcnlVc2FnZUJlZm9yZSwgd2lkZ2V0c0RpZmZlcmVuY2U6IGdldFdpZGdldHNDb3VudFNhZmUoKSAtIHdpZGdldHNCZWZvcmUgfSk7XHJcblxyXG4gICAgICB9XHJcbiAgICB9XHJcbiAgICByZXR1cm4gcmVzO1xyXG4gIH1cclxuXHJcbiAgZnVuY3Rpb24gZ2V0V2lkZ2V0c0NvdW50U2FmZSgpIHtcclxuICAgIGlmICh0eXBlb2YgcHJvY2VzcyAhPT0gJ3VuZGVmaW5lZCcpXHJcbiAgICAgIHJldHVybiAwO1xyXG4gICAgbGV0IGxlbmd0aCA9IC0xO1xyXG4gICAgdHJ5IHtcclxuICAgICAgbGVuZ3RoID0gREcuV2lkZ2V0LmdldEFsbCgpLmxlbmd0aDtcclxuICAgIH0gY2F0Y2ggKGU6IGFueSkge1xyXG4gICAgICBjb25zb2xlLndhcm4oZS5tZXNzYWdlID8/IGUpO1xyXG4gICAgfVxyXG4gICAgcmV0dXJuIGxlbmd0aDtcclxuICB9XHJcblxyXG4gIGFzeW5jIGZ1bmN0aW9uIGludm9rZVRlc3RzKGNhdGVnb3JpZXNUb0ludm9rZTogeyBba2V5OiBzdHJpbmddOiBDYXRlZ29yeSB9LCBvcHRpb25zOiBUZXN0RXhlY3V0aW9uT3B0aW9ucykge1xyXG4gICAgdHJ5IHtcclxuICAgICAgbGV0IHNraXBwaW5nQ2F0ZWdvcmllcyA9IG9wdGlvbnM/LnNraXBUb0NhdGVnb3J5ICE9IHVuZGVmaW5lZDtcclxuICAgICAgbGV0IGlzVGFyZ2V0Q2F0ZWdvcnkgPSBmYWxzZTtcclxuICAgICAgZm9yIChjb25zdCBba2V5LCB2YWx1ZV0gb2YgT2JqZWN0LmVudHJpZXMoY2F0ZWdvcmllc1RvSW52b2tlKSkge1xyXG4gICAgICAgICAgaWYgKG9wdGlvbnMuZXhjbHVkZT8uc29tZSgoYykgPT4ga2V5LnN0YXJ0c1dpdGgoYykpKVxyXG4gICAgICAgICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICAgICAgaWYgKG9wdGlvbnM/LmNhdGVnb3J5ICE9IG51bGwgJiYgIWtleS50b0xvd2VyQ2FzZSgpLnN0YXJ0c1dpdGgob3B0aW9ucz8uY2F0ZWdvcnkudG9Mb3dlckNhc2UoKS50cmltKCkpKVxyXG4gICAgICAgICAgICAgIGNvbnRpbnVlO1xyXG5cclxuICAgICAgICAgIGlmIChza2lwcGluZ0NhdGVnb3JpZXMpIHtcclxuICAgICAgICAgICAgICBpZiAoaXNUYXJnZXRDYXRlZ29yeSlcclxuICAgICAgICAgICAgICAgICAgc2tpcHBpbmdDYXRlZ29yaWVzID0gZmFsc2U7XHJcbiAgICAgICAgICAgICAgZWxzZSB7XHJcbiAgICAgICAgICAgICAgICAgIGlmIChvcHRpb25zPy5za2lwVG9DYXRlZ29yeSAhPSBudWxsICYmIGtleS50b0xvd2VyQ2FzZSgpLnRyaW0oKSA9PT0gb3B0aW9ucz8uc2tpcFRvQ2F0ZWdvcnkudG9Mb3dlckNhc2UoKS50cmltKCkpIHtcclxuICAgICAgICAgICAgICAgICAgICAgIGlzVGFyZ2V0Q2F0ZWdvcnkgPSB0cnVlO1xyXG4gICAgICAgICAgICAgICAgICB9IGVsc2Uge1xyXG4gICAgICAgICAgICAgICAgICAgICAgLy8gSGF2ZW4ndCBmb3VuZCB0aGUgdGFyZ2V0IGNhdGVnb3J5IHlldCwga2VlcCBza2lwcGluZ1xyXG4gICAgICAgICAgICAgICAgICAgICAgY29udGludWU7XHJcbiAgICAgICAgICAgICAgICAgIH1cclxuICAgICAgICAgICAgICB9XHJcbiAgICAgICAgICB9XHJcbiAgICAgICAgICBsZXQgdCA9ICh2YWx1ZS50ZXN0cyA/PyBbXSkuZmlsdGVyKChlKSA9PiBtYXRjaGVzTm9kZVRhcmdldChlLCB2YWx1ZSwgb3B0aW9ucykpO1xyXG4gICAgICAgICAgLy8gTm9kZS9icm93c2VyIHNwbGl0OiBpZiBubyB0ZXN0cyBvZiB0aGlzIGNhdGVnb3J5IGJlbG9uZyB0byB0aGUgcmVxdWVzdGVkIHRhcmdldCxcclxuICAgICAgICAgIC8vIHNraXAgaXQgZW50aXJlbHkg4oCUIGluY2x1ZGluZyBiZWZvcmUoKS9hZnRlcigpLlxyXG4gICAgICAgICAgaWYgKChvcHRpb25zLm5vZGVPbmx5IHx8IG9wdGlvbnMuZXhjbHVkZU5vZGVUZXN0cykgJiYgdC5sZW5ndGggPT09IDApXHJcbiAgICAgICAgICAgICAgY29udGludWU7XHJcblxyXG4gICAgICAgICAgY29uc3Qgc2tpcHBlZCA9IHZhbHVlLnRlc3RzID09IHVuZGVmaW5lZCA/IHVuZGVmaW5lZCA6IHQuZXZlcnkoKGU6IFRlc3QpID0+IGUub3B0aW9ucz8uc2tpcFJlYXNvblxyXG4gICAgICAgICAgICAgIHx8IChvcHRpb25zPy50ZXN0ICE9IG51bGwgJiYgb3B0aW9ucy50ZXN0LnRvTG93ZXJDYXNlKCkgIT09IGUubmFtZS50b0xvd2VyQ2FzZSgpKSk7XHJcblxyXG4gICAgICAgICAgaWYgKCFza2lwcGVkKSB7XHJcbiAgICAgICAgICAgICAgY29uc3Qgc2tpcHBlZENvdW50ID0gdC5maWx0ZXIoKGU6IFRlc3QpID0+XHJcbiAgICAgICAgICAgICAgICBlLm9wdGlvbnM/LnNraXBSZWFzb24gfHwgKG9wdGlvbnM/LnRlc3QgIT0gbnVsbCAmJiBvcHRpb25zLnRlc3QudG9Mb3dlckNhc2UoKSAhPT0gZS5uYW1lLnRvTG93ZXJDYXNlKCkpXHJcbiAgICAgICAgICAgICAgKS5sZW5ndGg7XHJcbiAgICAgICAgICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFN0YXJ0ZWQge3ske2tleX19fSR7c2tpcHBlZENvdW50ID4gMCA/IGAgc2tpcHBlZCB7eyR7c2tpcHBlZENvdW50fX19YCA6ICcnfWApO1xyXG4gICAgICAgICAgICAgIHZhbHVlLmJlZm9yZVN0YXR1cyA9IGF3YWl0IGludm9rZUNhdGVnb3J5TWV0aG9kKHZhbHVlLmJlZm9yZSwga2V5KTtcclxuICAgICAgICAgIH1cclxuXHJcbiAgICAgICAgICBpZiAob3B0aW9ucy5zdHJlc3NUZXN0KSB7XHJcbiAgICAgICAgICAgICAgdCA9IHQuZmlsdGVyKChlKSA9PiBlLm9wdGlvbnM/LnN0cmVzc1Rlc3QpO1xyXG4gICAgICAgICAgICAgIHQgPSBzaHVmZmxlKHQpO1xyXG4gICAgICAgICAgfVxyXG5cclxuICAgICAgICAgIGlmICgob3B0aW9ucy50YWdzPy5sZW5ndGggPz8gMCkgPiAwKSB7XHJcbiAgICAgICAgICAgICAgdCA9IHQuZmlsdGVyKChlKSA9PlxyXG4gICAgICAgICAgICAgICAgICBlLm9wdGlvbnM/LnRhZ3M/LnNvbWUodGFnID0+IChvcHRpb25zPy50YWdzID8/IFtdKS5pbmNsdWRlcyh0YWcpKVxyXG4gICAgICAgICAgICAgICk7XHJcbiAgICAgICAgICB9XHJcblxyXG4gICAgICAgICAgbGV0IHJlczogVGVzdFJlc3VsdEV4dGVuZGVkW107XHJcbiAgICAgICAgICBpZiAodmFsdWUuYmVmb3JlU3RhdHVzKSB7XHJcbiAgICAgICAgICAgICAgcmVzID0gQXJyYXkuZnJvbSh0Lm1hcCgodGVzdEVsZW0pID0+IHtcclxuICAgICAgICAgICAgICAgICAgcmV0dXJuIHtcclxuICAgICAgICAgICAgICAgICAgICAgIGRhdGU6IG5ldyBEYXRlKCkudG9JU09TdHJpbmcoKSxcclxuICAgICAgICAgICAgICAgICAgICAgIGNhdGVnb3J5OiBrZXksXHJcbiAgICAgICAgICAgICAgICAgICAgICBuYW1lOiB0ZXN0RWxlbS5uYW1lLFxyXG4gICAgICAgICAgICAgICAgICAgICAgc3VjY2VzczogZmFsc2UsXHJcbiAgICAgICAgICAgICAgICAgICAgICByZXN1bHQ6ICdiZWZvcmUoKSBmYWlsZWQnLFxyXG4gICAgICAgICAgICAgICAgICAgICAgbXM6IDAsXHJcbiAgICAgICAgICAgICAgICAgICAgICBza2lwcGVkOiBmYWxzZSxcclxuICAgICAgICAgICAgICAgICAgICAgIGxvZ3M6ICcnLFxyXG4gICAgICAgICAgICAgICAgICAgICAgb3duZXI6IHBhY2thZ2VPd25lcixcclxuICAgICAgICAgICAgICAgICAgICAgIHBhY2thZ2U6IHBhY2thZ2VfLm5hbWUsXHJcbiAgICAgICAgICAgICAgICAgICAgICB3aWRnZXRzRGlmZmVyZW5jZTogMCxcclxuICAgICAgICAgICAgICAgICAgICAgIGZsYWtpbmc6IERHLlRlc3QuaXNSZXByb2R1Y2luZ1xyXG4gICAgICAgICAgICAgICAgICB9O1xyXG4gICAgICAgICAgICAgIH0pKTtcclxuICAgICAgICAgICAgICByZXMuZm9yRWFjaChhc3luYyAodGVzdCkgPT4gYXdhaXQgZ3Jvay5zaGVsbC5yZXBvcnRUZXN0KCdwYWNrYWdlJywgdGVzdCkpO1xyXG4gICAgICAgICAgfSBlbHNlXHJcbiAgICAgICAgICAgICAgcmVzID0gYXdhaXQgaW52b2tlVGVzdHNJbkNhdGVnb3J5KHZhbHVlLCBvcHRpb25zLCBza2lwcGluZ0NhdGVnb3JpZXMpO1xyXG4gICAgICAgICAgY29uc3QgZGF0YTogVGVzdFJlc3VsdEV4dGVuZGVkW10gPSByZXMuZmlsdGVyKChkKSA9PiBkLnJlc3VsdCAhPSAnc2tpcHBlZCcpO1xyXG5cclxuICAgICAgICAgIGlmICghc2tpcHBlZClcclxuICAgICAgICAgICAgICB2YWx1ZS5hZnRlclN0YXR1cyA9IGF3YWl0IGludm9rZUNhdGVnb3J5TWV0aG9kKHZhbHVlLmFmdGVyLCBrZXkpO1xyXG5cclxuICAgICAgICAgIC8vIENsZWFyIGFmdGVyIGNhdGVnb3J5XHJcbiAgICAgICAgICAvLyBncm9rLnNoZWxsLmNsb3NlQWxsKCk7XHJcbiAgICAgICAgICAvLyBERy5CYWxsb29uLmNsb3NlQWxsKCk7XHJcbiAgICAgICAgICBpZiAodmFsdWUuYWZ0ZXJTdGF0dXMpIHtcclxuICAgICAgICAgICAgICBzdGRMb2coYFBhY2thZ2UgdGVzdGluZzogQ2F0ZWdvcnkgYWZ0ZXIoKSB7eyR7a2V5fX19IGZhaWxlZGApO1xyXG4gICAgICAgICAgICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBSZXN1bHQgZm9yIHt7JHtrZXl9fX0gYWZ0ZXI6ICR7dmFsdWUuYWZ0ZXJTdGF0dXN9YCk7XHJcbiAgICAgICAgICAgICAgZGF0YS5wdXNoKHtcclxuICAgICAgICAgICAgICAgICAgZGF0ZTogbmV3IERhdGUoKS50b0lTT1N0cmluZygpLFxyXG4gICAgICAgICAgICAgICAgICBjYXRlZ29yeToga2V5LFxyXG4gICAgICAgICAgICAgICAgICBuYW1lOiAnYWZ0ZXInLFxyXG4gICAgICAgICAgICAgICAgICBzdWNjZXNzOiBmYWxzZSxcclxuICAgICAgICAgICAgICAgICAgcmVzdWx0OiB2YWx1ZS5hZnRlclN0YXR1cyxcclxuICAgICAgICAgICAgICAgICAgbXM6IDAsXHJcbiAgICAgICAgICAgICAgICAgIHNraXBwZWQ6IGZhbHNlLFxyXG4gICAgICAgICAgICAgICAgICBsb2dzOiAnJyxcclxuICAgICAgICAgICAgICAgICAgb3duZXI6IHBhY2thZ2VPd25lcixcclxuICAgICAgICAgICAgICAgICAgcGFja2FnZTogcGFja2FnZV8ubmFtZSxcclxuICAgICAgICAgICAgICAgICAgd2lkZ2V0c0RpZmZlcmVuY2U6IDAsXHJcbiAgICAgICAgICAgICAgICAgIGZsYWtpbmc6IERHLlRlc3QuaXNSZXByb2R1Y2luZ1xyXG4gICAgICAgICAgICAgIH0pO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgICAgaWYgKHZhbHVlLmJlZm9yZVN0YXR1cykge1xyXG4gICAgICAgICAgICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBDYXRlZ29yeSBiZWZvcmUoKSB7eyR7a2V5fX19IGZhaWxlZGApO1xyXG4gICAgICAgICAgICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBSZXN1bHQgZm9yIHt7JHtrZXl9fX0gYmVmb3JlOiAke3ZhbHVlLmJlZm9yZVN0YXR1c31gKTtcclxuICAgICAgICAgICAgICBkYXRhLnB1c2goe1xyXG4gICAgICAgICAgICAgICAgICBkYXRlOiBuZXcgRGF0ZSgpLnRvSVNPU3RyaW5nKCksXHJcbiAgICAgICAgICAgICAgICAgIGNhdGVnb3J5OiBrZXksXHJcbiAgICAgICAgICAgICAgICAgIG5hbWU6ICdiZWZvcmUnLFxyXG4gICAgICAgICAgICAgICAgICBzdWNjZXNzOiBmYWxzZSxcclxuICAgICAgICAgICAgICAgICAgcmVzdWx0OiB2YWx1ZS5iZWZvcmVTdGF0dXMsXHJcbiAgICAgICAgICAgICAgICAgIG1zOiAwLFxyXG4gICAgICAgICAgICAgICAgICBza2lwcGVkOiBmYWxzZSxcclxuICAgICAgICAgICAgICAgICAgbG9nczogJycsXHJcbiAgICAgICAgICAgICAgICAgIG93bmVyOiBwYWNrYWdlT3duZXIsXHJcbiAgICAgICAgICAgICAgICAgIHBhY2thZ2U6IHBhY2thZ2VfLm5hbWUsXHJcbiAgICAgICAgICAgICAgICAgIHdpZGdldHNEaWZmZXJlbmNlOiAwLFxyXG4gICAgICAgICAgICAgICAgICBmbGFraW5nOiBERy5UZXN0LmlzUmVwcm9kdWNpbmdcclxuICAgICAgICAgICAgICB9KTtcclxuICAgICAgICAgIH1cclxuICAgICAgICAgIHJlc3VsdHMucHVzaCguLi5kYXRhKTtcclxuXHJcbiAgICAgICAgICAvLyBJZiByZXR1cm5PbkZhaWwgaXMgc2V0IGFuZCBhIHRlc3QgZmFpbGVkIChvdGhlciB0aGFuIHNraXBUb1Rlc3QpLCBzdG9wIHByb2Nlc3NpbmcgbW9yZSBjYXRlZ29yaWVzXHJcbiAgICAgICAgICBpZiAob3B0aW9ucy5yZXR1cm5PbkZhaWwgJiYgZGF0YS5zb21lKChkKSA9PiAhZC5zdWNjZXNzICYmICFkLnNraXBwZWQgJiYgZC5uYW1lICE9PSBvcHRpb25zLnNraXBUb1Rlc3QpKVxyXG4gICAgICAgICAgICAgIGJyZWFrO1xyXG4gICAgICB9XHJcbiAgICB9IGZpbmFsbHkge1xyXG4gICAgICByZXNldENvbnNvbGUoKTtcclxuICAgIH1cclxuICAgIGlmIChvcHRpb25zLnRlc3RDb250ZXh0IS5jYXRjaFVuaGFuZGxlZCAmJiAoIURHLlRlc3QuaXNJbkJlbmNobWFyaykpIHtcclxuICAgICAgYXdhaXQgZGVsYXkoMTAwMCk7XHJcbiAgICAgIGNvbnN0IGVycm9yID0gYXdhaXQgZ3Jvay5zaGVsbC5sYXN0RXJyb3I7XHJcbiAgICAgIGlmIChlcnJvciAhPSB1bmRlZmluZWQpIHtcclxuICAgICAgICAgIGNvbnN0IHBhcmFtczogYW55ID0ge1xyXG4gICAgICAgICAgICAgIGxvZ3M6ICcnLFxyXG4gICAgICAgICAgICAgIGRhdGU6IG5ldyBEYXRlKCkudG9JU09TdHJpbmcoKSxcclxuICAgICAgICAgICAgICBjYXRlZ29yeTogJ1VuaGFuZGxlZCBleGNlcHRpb25zJyxcclxuICAgICAgICAgICAgICBuYW1lOiAnRXhjZXB0aW9uJyxcclxuICAgICAgICAgICAgICByZXN1bHQ6IGVycm9yID8/ICcnLFxyXG4gICAgICAgICAgICAgIHN1Y2Nlc3M6ICFlcnJvcixcclxuICAgICAgICAgICAgICBtczogMCxcclxuICAgICAgICAgICAgICBza2lwcGVkOiBmYWxzZSxcclxuICAgICAgICAgICAgICBvd25lcjogcGFja2FnZU93bmVyID8/ICcnLFxyXG4gICAgICAgICAgICAgICdwYWNrYWdlJzogcGFja2FnZV8ubmFtZSxcclxuICAgICAgICAgICAgICB3aWRnZXRzRGlmZmVyZW5jZTogMFxyXG4gICAgICAgICAgfTtcclxuICAgICAgICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBVbmhhbmRsZWQgRXhjZXB0aW9uOiAke2Vycm9yfWApO1xyXG5cclxuICAgICAgICAgIHJlc3VsdHMucHVzaCh7Li4ucGFyYW1zLCAnZmxha2luZyc6IERHLlRlc3QuaXNSZXByb2R1Y2luZyAmJiAhZXJyb3J9KTtcclxuICAgICAgICAgICg8YW55PnBhcmFtcykucGFja2FnZSA9IHBhY2thZ2VfLm5hbWU7XHJcbiAgICAgICAgICBhd2FpdCBncm9rLnNoZWxsLnJlcG9ydFRlc3QoJ3BhY2thZ2UnLCBwYXJhbXMpO1xyXG4gICAgICB9XHJcbiAgICB9XHJcbiAgfVxyXG59XHJcblxyXG5hc3luYyBmdW5jdGlvbiBnZXRSZXN1bHQoeDogYW55KTogUHJvbWlzZTxzdHJpbmc+IHtcclxuICByZXR1cm4gYCR7eC50b1N0cmluZygpfVxcbiR7eC5zdGFjayA/IChhd2FpdCBERy5Mb2dnZXIudHJhbnNsYXRlU3RhY2tUcmFjZSh4LnN0YWNrKSkgOiAnJ31gO1xyXG59XHJcblxyXG5hc3luYyBmdW5jdGlvbiBleGVjVGVzdCh0OiBUZXN0LCBwcmVkaWNhdGU6IHN0cmluZyB8IHVuZGVmaW5lZCwgbG9nczogYW55W10sXHJcbiAgdGVzdFRpbWVvdXQ/OiBudW1iZXIsIHBhY2thZ2VOYW1lPzogc3RyaW5nLCB2ZXJib3NlPzogYm9vbGVhblxyXG4pOiBQcm9taXNlPFRlc3RSZXN1bHQgfCB1bmRlZmluZWQ+IHtcclxuICBsb2dzLmxlbmd0aCA9IDA7XHJcbiAgbGV0IHI6IFRlc3RSZXN1bHQ7XHJcbiAgbGV0IHR5cGU6IHN0cmluZyA9ICdwYWNrYWdlJztcclxuICBjb25zdCBmaWx0ZXIgPSBwcmVkaWNhdGUgIT0gdW5kZWZpbmVkICYmICh0Lm5hbWUudG9Mb3dlckNhc2UoKSAhPT0gcHJlZGljYXRlLnRvTG93ZXJDYXNlKCkpO1xyXG4gIGxldCBza2lwID0gdC5vcHRpb25zPy5za2lwUmVhc29uIHx8IGZpbHRlcjtcclxuICBsZXQgc2tpcFJlYXNvbiA9IGZpbHRlciA/ICdza2lwcGVkJyA6IHQub3B0aW9ucz8uc2tpcFJlYXNvbjtcclxuXHJcbiAgaWYgKERHLlRlc3QuaXNJbkJlbmNobWFyayAmJiAhdC5vcHRpb25zPy5iZW5jaG1hcmspIHtcclxuICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBTa2lwcGVkIHt7JHt0LmNhdGVnb3J5fX19IHt7JHt0Lm5hbWV9fX0gZG9lc250IGF2YWlsYWJsZSBpbiBiZW5jaG1hcmsgbW9kZWApO1xyXG4gICAgcmV0dXJuIHVuZGVmaW5lZDtcclxuICB9XHJcblxyXG4gIGlmIChza2lwICYmICFERy5UZXN0LmlzSW5CZW5jaG1hcmspXHJcbiAgICBzdGRMb2coYFBhY2thZ2UgdGVzdGluZzogU2tpcHBlZCB7eyR7dC5jYXRlZ29yeX19fSB7eyR7dC5uYW1lfX19YCk7XHJcbiAgaWYgKCFza2lwKVxyXG4gICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFN0YXJ0ZWQge3ske3QuY2F0ZWdvcnl9fX0ge3ske3QubmFtZX19fWApO1xyXG4gIGNvbnN0IHN0YXJ0ID0gRGF0ZS5ub3coKTtcclxuICBjb25zdCBzdGFydERhdGUgPSBuZXcgRGF0ZShzdGFydCkudG9JU09TdHJpbmcoKTtcclxuICB0cnkge1xyXG4gICAgaWYgKHNraXApXHJcbiAgICAgIHIgPSB7IG5hbWU6IHQubmFtZSwgb3duZXI6dC5vcHRpb25zPy5vd25lciA/PyAnJywgY2F0ZWdvcnk6IHQuY2F0ZWdvcnksIGxvZ3M6ICcnLCBkYXRlOiBzdGFydERhdGUsIHN1Y2Nlc3M6IHRydWUsIHJlc3VsdDogc2tpcFJlYXNvbiEsIG1zOiAwLCBza2lwcGVkOiB0cnVlLCBwYWNrYWdlOiBwYWNrYWdlTmFtZSA/PyAnJywgZmxha2luZzogREcuVGVzdC5pc1JlcHJvZHVjaW5nfTtcclxuICAgIGVsc2Uge1xyXG4gICAgICBsZXQgdGltZW91dF8gPSB0ZXN0VGltZW91dCA/PyBTVEFOREFSVF9USU1FT1VUO1xyXG5cclxuICAgICAgaWYgKERHLlRlc3QuaXNQcm9maWxpbmcpXHJcbiAgICAgICAgY29uc29sZS5wcm9maWxlKGAke3QuY2F0ZWdvcnl9OiAke3QubmFtZX1gKTtcclxuXHJcbiAgICAgIHIgPSB7IG5hbWU6IHQubmFtZSwgb3duZXI6dC5vcHRpb25zPy5vd25lciA/PyAnJywgY2F0ZWdvcnk6IHQuY2F0ZWdvcnksIGxvZ3M6ICcnLCBkYXRlOiBzdGFydERhdGUsIHN1Y2Nlc3M6IHRydWUsIHJlc3VsdDogKGF3YWl0IHRpbWVvdXQodC50ZXN0LCB0aW1lb3V0XykpLnRvU3RyaW5nKCkgPz8gJ09LJywgbXM6IDAsIHNraXBwZWQ6IGZhbHNlICwgcGFja2FnZTogcGFja2FnZU5hbWUgPz8gJycsIGZsYWtpbmc6IERHLlRlc3QuaXNSZXByb2R1Y2luZ307XHJcblxyXG4gICAgICBpZiAoREcuVGVzdC5pc1Byb2ZpbGluZykge1xyXG4gICAgICAgIGNvbnNvbGUucHJvZmlsZUVuZChgJHt0LmNhdGVnb3J5fTogJHt0Lm5hbWV9YCk7XHJcbiAgICAgICAgZ3Jvay5zaGVsbC5pbmZvKGBQcm9maWxpbmcgb2YgJHt0LmNhdGVnb3J5fTogJHt0Lm5hbWV9IGZpbmlzaGVkIFxcbiBQbGVhc2UgZW5zdXJlIHRoYXQgeW91IGhhdmUgb3BlbmVkIERldlRvb2xzIChGMTIpIC8gUGVyZm9ybWFuY2UgcGFuZWwgYmVmb3JlIHRlc3Qgc3RhcnRzLmApO1xyXG4gICAgICB9XHJcbiAgICB9XHJcbiAgfSBjYXRjaCAoeDogYW55KSB7XHJcbiAgICBzdGRFcnJvcih4KTtcclxuICAgIHIgPSB7IG5hbWU6IHQubmFtZSwgb3duZXI6dC5vcHRpb25zPy5vd25lciA/PyAnJywgY2F0ZWdvcnk6IHQuY2F0ZWdvcnksIGxvZ3M6ICcnLCBkYXRlOiBzdGFydERhdGUsIHN1Y2Nlc3M6IGZhbHNlLCByZXN1bHQ6IGF3YWl0IGdldFJlc3VsdCh4KSwgbXM6IDAsIHNraXBwZWQ6IGZhbHNlLCBwYWNrYWdlOiBwYWNrYWdlTmFtZSA/PyAnJywgZmxha2luZzogZmFsc2V9O1xyXG4gIH1cclxuICBpZiAodC5vcHRpb25zPy5pc0FnZ3JlZ2F0ZWQgJiYgci5yZXN1bHQuY29uc3RydWN0b3IgPT09IERHLkRhdGFGcmFtZSkge1xyXG4gICAgY29uc3QgY29sID0gci5yZXN1bHQuY29sKCdzdWNjZXNzJyk7XHJcbiAgICBpZiAoY29sKVxyXG4gICAgICByLnN1Y2Nlc3MgPSBjb2wuc3RhdHMuc3VtID09PSBjb2wubGVuZ3RoO1xyXG4gICAgaWYgKCF2ZXJib3NlKSB7XHJcbiAgICAgIGNvbnN0IGRmID0gci5yZXN1bHQ7XHJcbiAgICAgIGRmLmNvbHVtbnMucmVtb3ZlKCdzdGFjaycpO1xyXG4gICAgICBkZi5yb3dzLnJlbW92ZVdoZXJlKChyKSA9PiByLmdldCgnc3VjY2VzcycpKTtcclxuICAgICAgci5yZXN1bHQgPSBkZjtcclxuICAgIH1cclxuICAgIHIucmVzdWx0ID0gci5yZXN1bHQudG9Dc3YoKTtcclxuICB9XHJcbiAgci5sb2dzID0gbG9ncy5qb2luKCdcXG4nKTtcclxuICByLm1zID0gRGF0ZS5ub3coKSAtIHN0YXJ0O1xyXG4gIGlmICghc2tpcClcclxuICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBGaW5pc2hlZCB7eyR7dC5jYXRlZ29yeX19fSB7eyR7dC5uYW1lfX19IHdpdGgge3ske3Iuc3VjY2VzcyA/ICdzdWNjZXNzJyA6ICdlcnJvcid9fX0gZm9yICR7ci5tc30gbXNgKTtcclxuICBpZiAoIXIuc3VjY2Vzcykge1xyXG4gICAgICBzdGRMb2coYFBhY2thZ2UgdGVzdGluZzogUmVzdWx0IGZvciB7eyR7dC5jYXRlZ29yeX19fSB7eyR7dC5uYW1lfX19OiAke3IucmVzdWx0fWApO1xyXG4gIH1cclxuICByLmNhdGVnb3J5ID0gdC5jYXRlZ29yeTtcclxuICByLm5hbWUgPSB0Lm5hbWU7XHJcbiAgci5vd25lciA9IHQub3B0aW9ucz8ub3duZXIgPz8gJyc7XHJcbiAgaWYgKCFmaWx0ZXIpIHtcclxuICAgIGxldCBwYXJhbXMgPSB7XHJcbiAgICAgICdzdWNjZXNzJzogci5zdWNjZXNzLCAncmVzdWx0Jzogci5yZXN1bHQsICdtcyc6IHIubXMsICdkYXRlJzogci5kYXRlLFxyXG4gICAgICAnc2tpcHBlZCc6IHIuc2tpcHBlZCwgJ2NhdGVnb3J5JzogdC5jYXRlZ29yeSwgJ25hbWUnOiB0Lm5hbWUsICdsb2dzJzogci5sb2dzLCAnb3duZXInOiByLm93bmVyLFxyXG4gICAgICAnZmxha2luZyc6IERHLlRlc3QuaXNSZXByb2R1Y2luZyAmJiByLnN1Y2Nlc3MsXHJcbiAgICAgICdwYWNrYWdlJzogci5wYWNrYWdlXHJcbiAgICB9O1xyXG4gICAgaWYgKHIucmVzdWx0LmNvbnN0cnVjdG9yID09IE9iamVjdCkge1xyXG4gICAgICBjb25zdCByZXMgPSBPYmplY3Qua2V5cyhyLnJlc3VsdCkucmVkdWNlKChhY2MsIGspID0+ICh7IC4uLmFjYywgWydyZXN1bHQuJyArIGtdOiByLnJlc3VsdFtrXSB9KSwge30pO1xyXG4gICAgICBwYXJhbXMgPSB7IC4uLnBhcmFtcywgLi4ucmVzIH07XHJcbiAgICB9XHJcblxyXG4gICAgaWYgKHBhcmFtcy5yZXN1bHQgaW5zdGFuY2VvZiBERy5EYXRhRnJhbWUpXHJcbiAgICAgIHBhcmFtcy5yZXN1bHQgPSBKU09OLnN0cmluZ2lmeShwYXJhbXMucmVzdWx0Py50b0pzb24oKSkgfHwgJyc7XHJcbiAgICBhd2FpdCBncm9rLnNoZWxsLnJlcG9ydFRlc3QodHlwZSwgcGFyYW1zKTtcclxuICB9XHJcbiAgcmV0dXJuIHI7XHJcbn1cclxuXHJcbmV4cG9ydCBmdW5jdGlvbiBzaHVmZmxlKGFycmF5OiBhbnlbXSk6IGFueVtdIHtcclxuICBjb25zdCBuZXdBcnIgPSBhcnJheS5zbGljZSgpO1xyXG4gIG5ld0Fyci5zb3J0KCgpID0+IE1hdGgucmFuZG9tKCkgLSAwLjUpO1xyXG4gIHJldHVybiBuZXdBcnI7XHJcbn1cclxuXHJcbi8qIFdhaXRzIFttc10gbWlsbGlzZWNvbmRzICovXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBkZWxheShtczogbnVtYmVyKSB7XHJcbiAgYXdhaXQgbmV3IFByb21pc2UoKHIpID0+IHNldFRpbWVvdXQociwgbXMpKTtcclxufVxyXG5cclxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIGF3YWl0Q2hlY2soY2hlY2tIYW5kbGVyOiAoKSA9PiBib29sZWFuLFxyXG4gIGVycm9yOiBzdHJpbmcgPSAnVGltZW91dCBleGNlZWRlZCcsIHdhaXQ6IG51bWJlciA9IDUwMCwgaW50ZXJ2YWw6IG51bWJlciA9IDUwKTogUHJvbWlzZTxhbnk+IHtcclxuICByZXR1cm4gbmV3IFByb21pc2UoKHJlc29sdmUsIHJlamVjdCkgPT4ge1xyXG4gICAgc2V0VGltZW91dCgoKSA9PiB7XHJcbiAgICAgIGNsZWFySW50ZXJ2YWwoaW50ZXJ2YWxJZCk7XHJcbiAgICAgIHJlamVjdChuZXcgRXJyb3IoZXJyb3IpKTtcclxuICAgIH0sIHdhaXQpO1xyXG4gICAgLy8gQHRzLWlnbm9yZVxyXG4gICAgY29uc3QgaW50ZXJ2YWxJZDogVGltZW91dCA9IHNldEludGVydmFsKCgpID0+IHtcclxuICAgICAgaWYgKGNoZWNrSGFuZGxlcigpKSB7XHJcbiAgICAgICAgY2xlYXJJbnRlcnZhbChpbnRlcnZhbElkKTtcclxuICAgICAgICByZXNvbHZlKG51bGwpO1xyXG4gICAgICB9XHJcbiAgICB9LCBpbnRlcnZhbCk7XHJcbiAgfSk7XHJcbn1cclxuXHJcbi8vIFJldHVybnMgdGVzdCBleGVjdXRpb24gcmVzdWx0IG9yIGFuIGVycm9yIGluIGNhc2Ugb2YgdGltZW91dFxyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gdGltZW91dChmdW5jOiAoKSA9PiBQcm9taXNlPGFueT4sIHRlc3RUaW1lb3V0OiBudW1iZXIsIHRpbWVvdXRSZWFzb246IHN0cmluZyA9ICdFWEVDVVRJT04gVElNRU9VVCcpOiBQcm9taXNlPGFueT4ge1xyXG4gIGxldCB0aW1lb3V0OiBhbnkgPSBudWxsO1xyXG4gIGNvbnN0IHRpbWVvdXRQcm9taXNlID0gbmV3IFByb21pc2U8YW55PigoXywgcmVqZWN0KSA9PiB7XHJcbiAgICB0aW1lb3V0ID0gc2V0VGltZW91dCgoKSA9PiB7XHJcbiAgICAgIC8vIGVzbGludC1kaXNhYmxlLW5leHQtbGluZSBwcmVmZXItcHJvbWlzZS1yZWplY3QtZXJyb3JzXHJcbiAgICAgIHJlamVjdCh0aW1lb3V0UmVhc29uKTtcclxuICAgIH0sIHRlc3RUaW1lb3V0KTtcclxuICB9KTtcclxuICB0cnkge1xyXG4gICAgcmV0dXJuIGF3YWl0IFByb21pc2UucmFjZShbZnVuYygpLCB0aW1lb3V0UHJvbWlzZV0pO1xyXG4gIH0gZmluYWxseSB7XHJcbiAgICBpZiAodGltZW91dClcclxuICAgICAgY2xlYXJUaW1lb3V0KHRpbWVvdXQpO1xyXG4gIH1cclxufVxyXG5cclxuZXhwb3J0IGZ1bmN0aW9uIGlzRGlhbG9nUHJlc2VudChkaWFsb2dUaXRsZTogc3RyaW5nKTogYm9vbGVhbiB7XHJcbiAgY29uc3QgZGlhbG9ncyA9IERHLkRpYWxvZy5nZXRPcGVuRGlhbG9ncygpO1xyXG4gIGZvciAobGV0IGkgPSAwOyBpIDwgZGlhbG9ncy5sZW5ndGg7IGkrKykge1xyXG4gICAgaWYgKGRpYWxvZ3NbaV0udGl0bGUgPT0gZGlhbG9nVGl0bGUpXHJcbiAgICAgIHJldHVybiB0cnVlO1xyXG4gIH1cclxuICByZXR1cm4gZmFsc2U7XHJcbn1cclxuXHJcbi8qKiBFeHBlY3RzIGFuIGFzeW5jaHJvbm91cyB7QGxpbmsgYWN0aW9ufSB0byB0aHJvdyBhbiBleGNlcHRpb24uIFVzZSB7QGxpbmsgY2hlY2t9IHRvIHBlcmZvcm1cclxuICogZGVlcGVyIGluc3BlY3Rpb24gb2YgdGhlIGV4Y2VwdGlvbiBpZiBuZWNlc3NhcnkuXHJcbiAqIEBwYXJhbSAge2Z1bmN0aW9uKCk6IFByb21pc2U8dm9pZD59IGFjdGlvblxyXG4gKiBAcGFyYW0gIHtmdW5jdGlvbihhbnkpOiBib29sZWFufSBjaGVja1xyXG4gKiBAcmV0dXJuIHtQcm9taXNlPHZvaWQ+fVxyXG4gKi9cclxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIGV4cGVjdEV4Y2VwdGlvbkFzeW5jKGFjdGlvbjogKCkgPT4gUHJvbWlzZTx2b2lkPixcclxuICBjaGVjaz86IChleGNlcHRpb246IGFueSkgPT4gYm9vbGVhbik6IFByb21pc2U8dm9pZD4ge1xyXG4gIGxldCBjYXVnaHQ6IGJvb2xlYW4gPSBmYWxzZTtcclxuICBsZXQgY2hlY2tlZDogYm9vbGVhbiA9IGZhbHNlO1xyXG4gIHRyeSB7XHJcbiAgICBhd2FpdCBhY3Rpb24oKTtcclxuICB9IGNhdGNoIChlKSB7XHJcbiAgICBjYXVnaHQgPSB0cnVlO1xyXG4gICAgY2hlY2tlZCA9ICFjaGVjayB8fCBjaGVjayhlKTtcclxuICB9IGZpbmFsbHkge1xyXG4gICAgaWYgKCFjYXVnaHQpXHJcbiAgICAgIHRocm93IG5ldyBFcnJvcignQW4gZXhjZXB0aW9uIGlzIGV4cGVjdGVkIGJ1dCBub3QgdGhyb3duJyk7XHJcbiAgICBpZiAoIWNoZWNrZWQpXHJcbiAgICAgIHRocm93IG5ldyBFcnJvcignQW4gZXhwZWN0ZWQgZXhjZXB0aW9uIGlzIHRocm93biwgYnV0IGl0IGRvZXMgbm90IHNhdGlzZnkgdGhlIGNvbmRpdGlvbicpO1xyXG4gIH1cclxufVxyXG5cclxuY29uc3QgY2F0REYgPSBERy5EYXRhRnJhbWUuZnJvbUNvbHVtbnMoW0RHLkNvbHVtbi5mcm9tU3RyaW5ncygnY29sJywgWyd2YWwxJywgJ3ZhbDInLCAndmFsMyddKV0pO1xyXG5cclxuLyoqXHJcbiAqIFVuaXZlcnNhbCB0ZXN0IGZvciB2aWV3ZXJzLiBJdCBzZWFyY2ggdmlld2VycyBpbiBET00gYnkgdGFnczogY2FudmFzLCBzdmcsIGltZywgaW5wdXQsIGgxLCBhXHJcbiAqIEBwYXJhbSAge3N0cmluZ30gdiBWaWV3ZXIgbmFtZVxyXG4gKiBAcGFyYW0gIHtfREcuRGF0YUZyYW1lfSBkZiBEYXRhZnJhbWUgdG8gdXNlLiBTaG91bGQgaGF2ZSBhdCBsZWFzdCAzIHJvd3NcclxuICogQHBhcmFtICB7Ym9vbGVhbn0gb3B0aW9ucy5kZXRlY3RTZW1hbnRpY1R5cGVzIFNwZWNpZnkgd2hldGhlciB0byBkZXRlY3Qgc2VtYW50aWMgdHlwZXMgb3Igbm90XHJcbiAqIEBwYXJhbSAge2Jvb2xlYW59IG9wdGlvbnMucmVhZE9ubHkgSWYgc2V0IHRvIHRydWUsIHRoZSBkYXRhZnJhbWUgd2lsbCBub3QgYmUgbW9kaWZpZWQgZHVyaW5nIHRoZSB0ZXN0XHJcbiAqIEBwYXJhbSAge2Jvb2xlYW59IG9wdGlvbnMuYXJiaXRyYXJ5RGZUZXN0IElmIHNldCB0byBmYWxzZSwgdGVzdCBvbiBhcmJpdHJhcnkgZGF0YWZyYW1lXHJcbiAqIChvbmUgY2F0ZWdvcmljYWwgY29sdW1uKSB3aWxsIG5vdCBiZSBwZXJmb3JtZWRcclxuICogQHBhcmFtICB7b2JqZWN0fSBvcHRpb25zIExpc3Qgb2Ygb3B0aW9ucyAob3B0aW9uYWwpXHJcbiAqIEByZXR1cm4ge1Byb21pc2U8dm9pZD59IFRoZSB0ZXN0IGlzIGNvbnNpZGVyZWQgc3VjY2Vzc2Z1bCBpZiBpdCBjb21wbGV0ZXMgd2l0aG91dCBlcnJvcnNcclxuICovXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiB0ZXN0Vmlld2VyKHY6IHN0cmluZywgZGY6IF9ERy5EYXRhRnJhbWUsIG9wdGlvbnM/OiB7XHJcbiAgZGV0ZWN0U2VtYW50aWNUeXBlcz86IGJvb2xlYW4sIHJlYWRPbmx5PzogYm9vbGVhbiwgYXJiaXRyYXJ5RGZUZXN0PzogYm9vbGVhbixcclxuICBwYWNrYWdlTmFtZT86IHN0cmluZywgYXdhaXRWaWV3ZXI/OiAodmlld2VyOiBfREcuVmlld2VyKSA9PiBQcm9taXNlPHZvaWQ+XHJcbn0pOiBQcm9taXNlPHZvaWQ+IHtcclxuICBjb25zdCBwYWNrYWdlTmFtZSA9IG9wdGlvbnM/LnBhY2thZ2VOYW1lID8/ICcnO1xyXG4gIGlmIChvcHRpb25zPy5kZXRlY3RTZW1hbnRpY1R5cGVzKVxyXG4gICAgYXdhaXQgZ3Jvay5kYXRhLmRldGVjdFNlbWFudGljVHlwZXMoZGYpO1xyXG4gIGNvbnN0IHR2ID0gZ3Jvay5zaGVsbC5hZGRUYWJsZVZpZXcoZGYpO1xyXG5cclxuICB0cnkge1xyXG4gICAgLy8xLiBPcGVuLCBkbyBub3RoaW5nIGFuZCBjbG9zZVxyXG4gICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCk7XHJcbiAgICAvL2luIGNhc2Ugdmlld2VyIHdpdGggYXN5bmMgcmVuZGVyaW5nIC0gd2FpdCBmb3IgcmVuZGVyIHRvIGNvbXBsZXRlXHJcbiAgICBpZiAob3B0aW9ucz8uYXdhaXRWaWV3ZXIpXHJcbiAgICAgIGF3YWl0IHRlc3RWaWV3ZXJJbnRlcm5hbCh0diwgdiwgcGFja2FnZU5hbWUsIGdyb2suZXZlbnRzLm9uVmlld2VyQWRkZWQsIHVuZGVmaW5lZCwgb3B0aW9ucyEuYXdhaXRWaWV3ZXIpO1xyXG5cclxuICAgIC8vMi4gT3BlbiB2aWV3ZXIsIHJ1biBzZWxlY3Rpb24sIGZpbHRlciwgZXRjLiBhbmQgY2xvc2VcclxuICAgIGlmICghb3B0aW9ucz8ucmVhZE9ubHkpIHtcclxuICAgICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCwgc2VsZWN0RmlsdGVyQ2hhbmdlQ3VycmVudCk7XHJcbiAgICAgIGlmIChvcHRpb25zPy5hd2FpdFZpZXdlcilcclxuICAgICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCBzZWxlY3RGaWx0ZXJDaGFuZ2VDdXJyZW50LCBvcHRpb25zIS5hd2FpdFZpZXdlcik7XHJcbiAgICB9XHJcblxyXG4gICAgLy8yLiBPcGVuIHZpZXdlciwgY2hhbmdlIG9wdGlvbnMsIHNhdmUgbGF5b3V0IGFuZCBjbG9zZVxyXG4gICAgbGV0IHByb3BzQW5kTGF5b3V0OiB7IGxheW91dDogYW55LCBzYXZlZFByb3BzOiBhbnkgfSB8IG51bGwgPSBudWxsO1xyXG4gICAgcHJvcHNBbmRMYXlvdXQgPSBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCBjaGFuZ2VPcHRpb25zU2F2ZUxheW91dCk7XHJcbiAgICBpZiAob3B0aW9ucz8uYXdhaXRWaWV3ZXIpXHJcbiAgICAgIHByb3BzQW5kTGF5b3V0ID0gYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCxcclxuICAgICAgICBjaGFuZ2VPcHRpb25zU2F2ZUxheW91dCwgb3B0aW9ucyEuYXdhaXRWaWV3ZXIpXHJcblxyXG4gICAgLy8zLiBMb2FkIGxheW91dFxyXG4gICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3TGF5b3V0QXBwbGllZCwgbG9hZExheW91dCwgdW5kZWZpbmVkLCBwcm9wc0FuZExheW91dD8ubGF5b3V0LFxyXG4gICAgICB7IHNhdmVkUHJvcHM6IHByb3BzQW5kTGF5b3V0Py5zYXZlZFByb3BzIH0pO1xyXG4gICAgaWYgKG9wdGlvbnM/LmF3YWl0Vmlld2VyKVxyXG4gICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdMYXlvdXRBcHBsaWVkLCBsb2FkTGF5b3V0LCBvcHRpb25zIS5hd2FpdFZpZXdlcixcclxuICAgICAgICBwcm9wc0FuZExheW91dD8ubGF5b3V0LCB7IHNhdmVkUHJvcHM6IHByb3BzQW5kTGF5b3V0Py5zYXZlZFByb3BzIH0pO1xyXG5cclxuICAgIC8vNC4gT3BlbiB2aWV3ZXIgb24gYXJiaXRhcnkgZGF0YXNldFxyXG4gICAgaWYgKG9wdGlvbnM/LmFyYml0cmFyeURmVGVzdCAhPT0gZmFsc2UpIHtcclxuICAgICAgdHYuZGF0YUZyYW1lID0gY2F0REY7XHJcbiAgICAgIGF3YWl0IGRlbGF5KDUwKTtcclxuICAgICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCk7XHJcbiAgICAgIGlmIChvcHRpb25zPy5hd2FpdFZpZXdlcilcclxuICAgICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCB1bmRlZmluZWQsIG9wdGlvbnMhLmF3YWl0Vmlld2VyKTtcclxuICAgIH1cclxuXHJcbiAgICAvLzUuIENhbGwgcG9zdHBvbmVkIGZpbHRlcmluZ1xyXG4gICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCwgZmlsdGVyQXN5bmMpO1xyXG4gICAgaWYgKG9wdGlvbnM/LmF3YWl0Vmlld2VyKVxyXG4gICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCBmaWx0ZXJBc3luYywgb3B0aW9ucyEuYXdhaXRWaWV3ZXIpO1xyXG5cclxuICB9IGZpbmFsbHkge1xyXG4gICAgLy8gY2xvc2VBbGwoKSBpcyBoYW5kbGluZyBieSBjb21tb24gdGVzdCB3b3JrZmxvd1xyXG4gICAgLy8gZ3Jvay5zaGVsbC5jbG9zZUFsbCgpO1xyXG4gICAgLy8gREcuQmFsbG9vbi5jbG9zZUFsbCgpO1xyXG4gIH1cclxufVxyXG4iXX0=