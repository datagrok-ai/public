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
                    }), { skipReason: map.skip, timeout: DG.Test.isInBenchmark ? (_j = map.benchmarkTimeout) !== null && _j !== void 0 ? _j : BENCHMARK_TIMEOUT : (_k = map.timeout) !== null && _k !== void 0 ? _k : STANDART_TIMEOUT });
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
            var _a, _b, _c, _d, _e, _f, _g, _h, _j, _k, _l, _m, _o, _p, _q, _r, _s, _t;
            return __awaiter(this, void 0, void 0, function* () {
                let t = (_a = category.tests) !== null && _a !== void 0 ? _a : [];
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
                        if (!options.nodeOptions) {
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
                            test.options.owner = (_q = (_p = (_o = (_m = t[i].options) === null || _m === void 0 ? void 0 : _m.owner) !== null && _o !== void 0 ? _o : category === null || category === void 0 ? void 0 : category.owner) !== null && _p !== void 0 ? _p : packageOwner) !== null && _q !== void 0 ? _q : '';
                        }
                        // let isGBEnable = (window as any).gc && test.options?.skipReason == undefined;
                        // console.log(`********${isGBEnable}`);
                        // if (isGBEnable)
                        //   await (window as any).gc();
                        // memoryUsageBefore = (window?.performance as any)?.memory?.usedJSHeapSize;
                        let testRun = yield execTest(test, options === null || options === void 0 ? void 0 : options.test, logs, DG.Test.isInBenchmark ? (_s = (_r = t[i].options) === null || _r === void 0 ? void 0 : _r.benchmarkTimeout) !== null && _s !== void 0 ? _s : BENCHMARK_TIMEOUT : (_t = t[i].options) === null || _t === void 0 ? void 0 : _t.timeout, package_.name, options.verbose);
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
            var _a, _b, _c, _d, _e, _f;
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
                        //@ts-ignore
                        const skipped = (_b = value.tests) === null || _b === void 0 ? void 0 : _b.every((t) => {
                            var _a;
                            return ((_a = t.options) === null || _a === void 0 ? void 0 : _a.skipReason)
                                || ((options === null || options === void 0 ? void 0 : options.test) != null && options.test.toLowerCase() !== t.name.toLowerCase());
                        });
                        if (!skipped) {
                            //@ts-ignore
                            const skippedCount = ((_c = value.tests) !== null && _c !== void 0 ? _c : []).filter((t) => { var _a; return ((_a = t.options) === null || _a === void 0 ? void 0 : _a.skipReason) || ((options === null || options === void 0 ? void 0 : options.test) != null && options.test.toLowerCase() !== t.name.toLowerCase()); }).length;
                            stdLog(`Package testing: Started {{${key}}}${skippedCount > 0 ? ` skipped {{${skippedCount}}}` : ''}`);
                            value.beforeStatus = yield invokeCategoryMethod(value.before, key);
                        }
                        let t = (_d = value.tests) !== null && _d !== void 0 ? _d : [];
                        if (options.stressTest) {
                            t = t.filter((e) => { var _a; return (_a = e.options) === null || _a === void 0 ? void 0 : _a.stressTest; });
                            t = shuffle(t);
                        }
                        if (((_f = (_e = options.tags) === null || _e === void 0 ? void 0 : _e.length) !== null && _f !== void 0 ? _f : 0) > 0) {
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
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoidGVzdC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzIjpbInRlc3QudHMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6Ijs7Ozs7Ozs7O0FBS0EsT0FBTyxFQUFFLFFBQVEsRUFBRSxNQUFNLG1CQUFtQixDQUFDO0FBRTdDLE9BQU8sRUFBRSx1QkFBdUIsRUFBRSxXQUFXLEVBQUUsVUFBVSxFQUFFLHlCQUF5QixFQUFFLGtCQUFrQixFQUFFLE1BQU0scUJBQXFCLENBQUM7QUFFdEksTUFBTSxnQkFBZ0IsR0FBRyxLQUFLLENBQUM7QUFDL0IsTUFBTSxpQkFBaUIsR0FBRyxRQUFRLENBQUM7QUFFbkMsTUFBTSxNQUFNLEdBQUcsT0FBTyxDQUFDLEdBQUcsQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDekMsTUFBTSxPQUFPLEdBQUcsT0FBTyxDQUFDLElBQUksQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDM0MsTUFBTSxPQUFPLEdBQUcsT0FBTyxDQUFDLElBQUksQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDM0MsTUFBTSxRQUFRLEdBQUcsT0FBTyxDQUFDLEtBQUssQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFFN0MsTUFBTSxDQUFDLE1BQU0sS0FBSyxHQUVkLEVBQUUsQ0FBQztBQUVQLE1BQU0sZ0JBQWdCLEdBQUcsWUFBWSxDQUFDO0FBQ3RDLE1BQU0sV0FBVyxHQUFHLE1BQU0sQ0FBQztBQUMzQixNQUFNLGdCQUFnQixHQUFHLFdBQVcsQ0FBQztBQUNyQyxNQUFNLFdBQVcsR0FBRyxNQUFNLENBQUM7QUFDM0IsTUFBTSxhQUFhLEdBQStCLEVBQUUsQ0FBQztBQUNyRCxNQUFNLENBQUMsSUFBSSxlQUF1QixDQUFDO0FBRW5DLE1BQU0sS0FBVyxNQUFNLENBS3RCO0FBTEQsV0FBaUIsTUFBTTtJQUNyQixTQUFnQixPQUFPLENBQUMsS0FBVSxFQUFFLElBQWE7UUFDL0MsSUFBSSxLQUFLLElBQUksSUFBSTtZ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QyxFQUFFLFdBQVcsRUFBRSxJQUFJLENBQUMsTUFBTSxDQUFDLGFBQWEsRUFBRSxXQUFXLENBQUMsQ0FBQztZQUNyRixJQUFJLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxXQUFXO2dCQUN0QixNQUFNLGtCQUFrQixDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUUsV0FBVyxFQUFFLElBQUksQ0FBQyxNQUFNLENBQUMsYUFBYSxFQUFFLFdBQVcsRUFBRSxPQUFRLENBQUMsV0FBVyxDQUFDLENBQUM7U0FFOUc7Z0JBQVM7WUFDUixpREFBaUQ7WUFDakQseUJBQXlCO1lBQ3pCLHlCQUF5QjtTQUMxQjs7Q0FDRiIsInNvdXJjZXNDb250ZW50IjpbImltcG9ydCB0eXBlICogYXMgX2dyb2sgZnJvbSAnZGF0YWdyb2stYXBpL2dyb2snO1xuaW1wb3J0IHR5cGUgKiBhcyBfREcgZnJvbSAnZGF0YWdyb2stYXBpL2RnJztcbmRlY2xhcmUgbGV0IGdyb2s6IHR5cGVvZiBfZ3JvaywgREc6IHR5cGVvZiBfREc7XG5cbmltcG9ydCB7IE9ic2VydmFibGUgfSBmcm9tICdyeGpzJztcbmltcG9ydCB7IHRlc3REYXRhIH0gZnJvbSAnLi9kYXRhZnJhbWUtdXRpbHMnO1xuaW1wb3J0IFRpbWVvdXQgPSBOb2RlSlMuVGltZW91dDtcbmltcG9ydCB7IGNoYW5nZU9wdGlvbnNTYXZlTGF5b3V0LCBmaWx0ZXJBc3luYywgbG9hZExheW91dCwgc2VsZWN0RmlsdGVyQ2hhbmdlQ3VycmVudCwgdGVzdFZpZXdlckludGVybmFsIH0gZnJvbSAnLi90ZXN0LXZpZXdlci11dGlscyc7XG5cbmNvbnN0IFNUQU5EQVJUX1RJTUVPVVQgPSAzMDAwMDtcbmNvbnN0IEJFTkNITUFSS19USU1FT1VUID0gMTA4MDAwMDA7XG5cbmNvbnN0IHN0ZExvZyA9IGNvbnNvbGUubG9nLmJpbmQoY29uc29sZSk7XG5jb25zdCBzdGRJbmZvID0gY29uc29sZS5pbmZvLmJpbmQoY29uc29sZSk7XG5jb25zdCBzdGRXYXJuID0gY29uc29sZS53YXJuLmJpbmQoY29uc29sZSk7XG5jb25zdCBzdGRFcnJvciA9IGNvbnNvbGUuZXJyb3IuYmluZChjb25zb2xlKTtcblxuZXhwb3J0IGNvbnN0IHRlc3RzOiB7XG4gIFtrZXk6IHN0cmluZ106IENhdGVnb3J5XG59ID0ge307XG5cbmNvbnN0IGF1dG9UZXN0c0NhdE5hbWUgPSAnQXV0byBUZXN0cyc7XG5jb25zdCBkZW1vQ2F0TmFtZSA9ICdEZW1vJztcbmNvbnN0IGRldGVjdG9yc0NhdE5hbWUgPSAnRGV0ZWN0b3JzJztcbmNvbnN0IGNvcmVDYXROYW1lID0gJ0NvcmUnO1xuY29uc3Qgd2FzUmVnaXN0ZXJlZDogeyBba2V5OiBzdHJpbmddOiBib29sZWFuIH0gPSB7fTtcbmV4cG9ydCBsZXQgY3VycmVudENhdGVnb3J5OiBzdHJpbmc7XG5cbmV4cG9ydCBuYW1lc3BhY2UgYXNzdXJlIHtcbiAgZXhwb3J0IGZ1bmN0aW9uIG5vdE51bGwodmFsdWU6IGFueSwgbmFtZT86IHN0cmluZykge1xuICAgIGlmICh2YWx1ZSA9PSBudWxsKVxuICAgICAgdGhyb3cgbmV3IEVycm9yKGAke25hbWUgPT0gbnVsbCA/ICdWYWx1ZScgOiBuYW1lfSBub3QgZGVmaW5lZGApO1xuICB9XG59XG5cbmV4cG9ydCBpbnRlcmZhY2UgVGVzdE9wdGlvbnMge1xuICB0aW1lb3V0PzogbnVtYmVyO1xuICBiZW5jaG1hcmtXYXJuVGltZW91dD86IG51bWJlcjtcbiAgYmVuY2htYXJrVGltZW91dD86IG51bWJlcjtcbiAgdW5oYW5kbGVkRXhjZXB0aW9uVGltZW91dD86IG51bWJlcjtcbiAgc2tpcFJlYXNvbj86IHN0cmluZztcbiAgaXNBZ2dyZWdhdGVkPzogYm9vbGVhbjtcbiAgYmVuY2htYXJrPzogYm9vbGVhbjtcbiAgc3RyZXNzVGVzdD86IGJvb2xlYW47XG4gIG93bmVyPzogc3RyaW5nO1xuICB0YWdzPzogc3RyaW5nW107XG59XG5cbmV4cG9ydCBpbnRlcmZhY2UgVGVzdFJlc3VsdCB7XG4gIGRhdGU6IHN0cmluZztcbiAgY2F0ZWdvcnk6IHN0cmluZztcbiAgbmFtZTogc3RyaW5nO1xuICBzdWNjZXNzOiBib29sZWFuO1xuICByZXN1bHQ6IGFueTtcbiAgbXM6IG51bWJlcjtcbiAgc2tpcHBlZDogYm9vbGVhbjtcbiAgbG9nczogc3RyaW5nO1xuICBvd25lcjogc3RyaW5nO1xuICBwYWNrYWdlOiBzdHJpbmc7XG4gIGZsYWtpbmc6IGJvb2xlYW47XG59XG5cblxuZXhwb3J0IGludGVyZmFjZSBUZXN0UmVzdWx0RXh0ZW5kZWQgZXh0ZW5kcyBUZXN0UmVzdWx0e1xuICB3aWRnZXRzRGlmZmVyZW5jZTogbnVtYmVyO1xufVxuXG5leHBvcnQgaW50ZXJmYWNlIENhdGVnb3J5T3B0aW9ucyB7XG4gIGNsZWFyPzogYm9vbGVhbjtcbiAgdGltZW91dD86IG51bWJlcjtcbiAgYmVuY2htYXJrcz86IGJvb2xlYW47XG4gIHN0cmVzc1Rlc3RzPzogYm9vbGVhbjtcbiAgb3duZXI/OiBzdHJpbmc7XG59XG5cbmV4cG9ydCBjbGFzcyBUZXN0Q29udGV4dCB7XG4gIHN0cmVzc1Rlc3Q/OiBib29sZWFuO1xuICBjYXRjaFVuaGFuZGxlZCA9IHRydWU7XG4gIHJlcG9ydCA9IGZhbHNlO1xuICByZXR1cm5PbkZhaWwgPSBmYWxzZTtcblxuICBjb25zdHJ1Y3RvcihjYXRjaFVuaGFuZGxlZD86IGJvb2xlYW4sIHJlcG9ydD86IGJvb2xlYW4sIHJldHVybk9uRmFpbD86IGJvb2xlYW4pIHtcbiAgICBpZiAoY2F0Y2hVbmhhbmRsZWQgIT09IHVuZGVmaW5lZCkgdGhpcy5jYXRjaFVuaGFuZGxlZCA9IGNhdGNoVW5oYW5kbGVkO1xuICAgIGlmIChyZXBvcnQgIT09IHVuZGVmaW5lZCkgdGhpcy5yZXBvcnQgPSByZXBvcnQ7XG4gICAgaWYgKHJldHVybk9uRmFpbCAhPT0gdW5kZWZpbmVkKSB0aGlzLnJldHVybk9uRmFpbCA9IHJldHVybk9uRmFpbDtcbiAgfTtcbn1cblxuZXhwb3J0IGNsYXNzIFRlc3Qge1xuICB0ZXN0OiAoKSA9PiBQcm9taXNlPGFueT47XG4gIG5hbWU6IHN0cmluZztcbiAgY2F0ZWdvcnk6IHN0cmluZztcbiAgb3B0aW9ucz86IFRlc3RPcHRpb25zO1xuXG4gIGNvbnN0cnVjdG9yKGNhdGVnb3J5OiBzdHJpbmcsIG5hbWU6IHN0cmluZywgdGVzdDogKCkgPT4gUHJvbWlzZTxhbnk+LCBvcHRpb25zPzogVGVzdE9wdGlvbnMpIHtcbiAgICB0aGlzLmNhdGVnb3J5ID0gY2F0ZWdvcnk7XG4gICAgdGhpcy5uYW1lID0gbmFtZTtcbiAgICBvcHRpb25zID8/PSB7fTtcbiAgICBvcHRpb25zLnRpbWVvdXQgPz89IFNUQU5EQVJUX1RJTUVPVVQ7XG4gICAgdGhpcy5vcHRpb25zID0gb3B0aW9ucztcbiAgICB0aGlzLnRlc3QgPSBhc3luYyAoKTogUHJvbWlzZTxhbnk+ID0+IHtcbiAgICAgIHJldHVybiBuZXcgUHJvbWlzZShhc3luYyAocmVzb2x2ZSwgcmVqZWN0KSA9PiB7XG4gICAgICAgIGxldCByZXN1bHQgPSAnJztcbiAgICAgICAgdHJ5IHtcbiAgICAgICAgICBpZiAoREcuVGVzdC5pc0luRGVidWcpXG4gICAgICAgICAgICBkZWJ1Z2dlcjtcblxuICAgICAgICAgIGxldCByZXMgPSBhd2FpdCB0ZXN0KCk7XG4gICAgICAgICAgdHJ5IHtcbiAgICAgICAgICAgIHJlc3VsdCA9IHJlcz8udG9TdHJpbmcoKSA/PyAnJztcbiAgICAgICAgICB9XG4gICAgICAgICAgY2F0Y2ggKGUpIHtcbiAgICAgICAgICAgIHJlc3VsdCA9ICdDYW5cXCd0IGNvbnZlcnQgdGVzdFxcJ3MgcmVzdWx0IHRvIHN0cmluZyc7XG4gICAgICAgICAgICBjb25zb2xlLmVycm9yKGBDYW5cXCd0IGNvbnZlcnQgdGVzdFxcJ3MgcmVzdWx0IHRvIHN0cmluZyBpbiB0aGUgJHt0aGlzLmNhdGVnb3J5fToke3RoaXMubmFtZX0gdGVzdGApO1xuICAgICAgICAgIH1cbiAgICAgICAgfSBjYXRjaCAoZTogYW55KSB7XG4gICAgICAgICAgcmVqZWN0KGUpO1xuICAgICAgICB9XG4gICAgICAgIHJlc29sdmUocmVzdWx0KTtcbiAgICAgIH0pO1xuICAgIH07XG4gIH1cbn1cblxuZXhwb3J0IGNsYXNzIENhdGVnb3J5IHtcbiAgdGVzdHM/OiBUZXN0W107XG4gIGJlZm9yZT86ICgpID0+IFByb21pc2U8dm9pZD47XG4gIGFmdGVyPzogKCkgPT4gUHJvbWlzZTx2b2lkPjtcblxuICBiZWZvcmVTdGF0dXM/OiBzdHJpbmc7XG4gIGFmdGVyU3RhdHVzPzogc3RyaW5nO1xuICBjbGVhcj86IGJvb2xlYW47XG4gIHRpbWVvdXQ/OiBudW1iZXI7XG4gIGJlbmNobWFya3M/OiBib29sZWFuO1xuICBiZW5jaG1hcmtUaW1lb3V0PzogbnVtYmVyO1xuICBzdHJlc3NUZXN0cz86IGJvb2xlYW47XG4gIG93bmVyPzogc3RyaW5nO1xufVxuXG5leHBvcnQgY2xhc3MgTm9kZVRlc3RFeGVjdXRpb25PcHRpb25zIHtcbiAgcGFja2FnZSE6IF9ERy5QYWNrYWdlO1xufVxuXG5leHBvcnQgY2xhc3MgVGVzdEV4ZWN1dGlvbk9wdGlvbnMge1xuICBjYXRlZ29yeT86IHN0cmluZztcbiAgdGVzdD86IHN0cmluZztcbiAgdGVzdENvbnRleHQ/OiBUZXN0Q29udGV4dDtcbiAgZXhjbHVkZT86IHN0cmluZ1tdO1xuICB2ZXJib3NlPzogYm9vbGVhbjtcbiAgc3RyZXNzVGVzdD86IGJvb2xlYW47XG4gIHRhZ3M/OiBzdHJpbmdbXTtcbiAgbm9kZU9wdGlvbnM/OiBOb2RlVGVzdEV4ZWN1dGlvbk9wdGlvbnM7XG4gIHNraXBUb0NhdGVnb3J5Pzogc3RyaW5nO1xuICBza2lwVG9UZXN0Pzogc3RyaW5nO1xuICByZXR1cm5PbkZhaWw/OiBib29sZWFuO1xufVxuXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gdGVzdEV2ZW50PFQ+KGV2ZW50OiBPYnNlcnZhYmxlPFQ+LFxuICBoYW5kbGVyOiAoYXJnczogVCkgPT4gdm9pZCwgdHJpZ2dlcjogKCkgPT4gdm9pZCwgbXM6IG51bWJlciA9IDAsIHJlYXNvbjogc3RyaW5nID0gYHRpbWVvdXRgXG4pOiBQcm9taXNlPGFueT4ge1xuICByZXR1cm4gbmV3IFByb21pc2UoKHJlc29sdmUsIHJlamVjdCkgPT4ge1xuICAgIGNvbnN0IHN1YiA9IGV2ZW50LnN1YnNjcmliZSgoYXJnczogVCkgPT4ge1xuICAgICAgdHJ5IHtcbiAgICAgICAgaGFuZGxlcihhcmdzKTtcbiAgICAgICAgcmVzb2x2ZSgnT0snKTtcbiAgICAgIH0gY2F0Y2ggKGUpIHtcbiAgICAgICAgcmVqZWN0KGUpO1xuICAgICAgfSBmaW5hbGx5IHtcbiAgICAgICAgc3ViLnVuc3Vic2NyaWJlKCk7XG4gICAgICAgIGNsZWFyVGltZW91dCh0aW1lb3V0KTtcbiAgICAgIH1cbiAgICB9KTtcbiAgICBjb25zdCB0aW1lb3V0ID0gc2V0VGltZW91dCgoKSA9PiB7XG4gICAgICBzdWIudW5zdWJzY3JpYmUoKTtcbiAgICAgIC8vIGVzbGludC1kaXNhYmxlLW5leHQtbGluZSBwcmVmZXItcHJvbWlzZS1yZWplY3QtZXJyb3JzXG4gICAgICByZWplY3QocmVhc29uKTtcbiAgICB9LCBtcyk7XG4gICAgdHJpZ2dlcigpO1xuICB9KTtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHRlc3RFdmVudEFzeW5jPFQ+KGV2ZW50OiBPYnNlcnZhYmxlPFQ+LFxuICBoYW5kbGVyOiAoYXJnczogVCkgPT4gUHJvbWlzZTx2b2lkPiwgdHJpZ2dlcjogKCkgPT4gdm9pZCwgbXM6IG51bWJlciA9IDAsIHJlYXNvbjogc3RyaW5nID0gYHRpbWVvdXRgXG4pOiBQcm9taXNlPGFueT4ge1xuICByZXR1cm4gbmV3IFByb21pc2UoKHJlc29sdmUsIHJlamVjdCkgPT4ge1xuICAgIGNvbnN0IHN1YiA9IGV2ZW50LnN1YnNjcmliZSgoYXJnczogVCkgPT4ge1xuICAgICAgaGFuZGxlcihhcmdzKS50aGVuKCgpID0+IHtcbiAgICAgICAgcmVzb2x2ZSgnT0snKTtcbiAgICAgIH0pLmNhdGNoKChlKSA9PiB7XG4gICAgICAgIHJlamVjdChlKTtcbiAgICAgIH0pLmZpbmFsbHkoKCkgPT4ge1xuICAgICAgICBzdWIudW5zdWJzY3JpYmUoKTtcbiAgICAgICAgY2xlYXJUaW1lb3V0KHRpbWVvdXQpO1xuICAgICAgfSk7XG4gICAgfSk7XG4gICAgY29uc3QgdGltZW91dCA9IHNldFRpbWVvdXQoKCkgPT4ge1xuICAgICAgc3ViLnVuc3Vic2NyaWJlKCk7XG4gICAgICAvLyBlc2xpbnQtZGlzYWJsZS1uZXh0LWxpbmUgcHJlZmVyLXByb21pc2UtcmVqZWN0LWVycm9yc1xuICAgICAgcmVqZWN0KHJlYXNvbik7XG4gICAgfSwgbXMpO1xuICAgIHRyaWdnZXIoKTtcbiAgfSk7XG59XG5cbmV4cG9ydCBmdW5jdGlvbiB0ZXN0KG5hbWU6IHN0cmluZywgdGVzdDogKCkgPT4gUHJvbWlzZTxhbnk+LCBvcHRpb25zPzogVGVzdE9wdGlvbnMpOiB2b2lkIHtcbiAgaWYgKHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0gPT0gdW5kZWZpbmVkKVxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0gPSB7fTtcbiAgaWYgKHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0udGVzdHMgPT0gdW5kZWZpbmVkKVxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0udGVzdHMgPSBbXTtcbiAgdGVzdHNbY3VycmVudENhdGVnb3J5XS50ZXN0cyEucHVzaChuZXcgVGVzdChjdXJyZW50Q2F0ZWdvcnksIG5hbWUsIHRlc3QsIG9wdGlvbnMpKTtcbn1cblxuLyogVGVzdHMgdHdvIG9iamVjdHMgZm9yIGVxdWFsaXR5LCB0aHJvd3MgYW4gZXhjZXB0aW9uIGlmIHRoZXkgYXJlIG5vdCBlcXVhbC4gKi9cbmV4cG9ydCBmdW5jdGlvbiBleHBlY3QoYWN0dWFsOiBhbnksIGV4cGVjdGVkOiBhbnkgPSB0cnVlLCBlcnJvcj86IHN0cmluZyk6IHZvaWQge1xuICBpZiAoZXJyb3IpXG4gICAgZXJyb3IgPSBgJHtlcnJvcn0sIGA7XG4gIGVsc2UgZXJyb3IgPSAnJztcbiAgaWYgKGFjdHVhbCAhPT0gZXhwZWN0ZWQpXG4gICAgdGhyb3cgbmV3IEVycm9yKGAke2Vycm9yfUV4cGVjdGVkIFwiJHtleHBlY3RlZH1cIiwgZ290IFwiJHthY3R1YWx9XCJgKTtcbn1cblxuZXhwb3J0IGZ1bmN0aW9uIGV4cGVjdEZsb2F0KGFjdHVhbDogbnVtYmVyLCBleHBlY3RlZDogbnVtYmVyLCB0b2xlcmFuY2UgPSAwLjAwMSwgZXJyb3I/OiBzdHJpbmcpOiB2b2lkIHtcbiAgaWYgKChhY3R1YWwgPT09IE51bWJlci5QT1NJVElWRV9JTkZJTklUWSAmJiBleHBlY3RlZCA9PT0gTnVtYmVyLlBPU0lUSVZFX0lORklOSVRZKSB8fFxuICAgIChhY3R1YWwgPT09IE51bWJlci5ORUdBVElWRV9JTkZJTklUWSAmJiBleHBlY3RlZCA9PT0gTnVtYmVyLk5FR0FUSVZFX0lORklOSVRZKSB8fFxuICAgIChhY3R1YWwgPT09IE51bWJlci5OYU4gJiYgZXhwZWN0ZWQgPT09IE51bWJlci5OYU4pIHx8IChpc05hTihhY3R1YWwpICYmIGlzTmFOKGV4cGVjdGVkKSkpXG4gICAgcmV0dXJuO1xuICBjb25zdCBhcmVFcXVhbCA9IE1hdGguYWJzKGFjdHVhbCAtIGV4cGVjdGVkKSA8IHRvbGVyYW5jZTtcbiAgZXhwZWN0KGFyZUVxdWFsLCB0cnVlLCBgJHtlcnJvciA/PyAnJ30gKHRvbGVyYW5jZSA9ICR7dG9sZXJhbmNlfTsgYSA9ICR7YWN0dWFsfSwgZSA9ICR7ZXhwZWN0ZWR9KWApO1xuICBpZiAoIWFyZUVxdWFsKVxuICAgIHRocm93IG5ldyBFcnJvcihgRXhwZWN0ZWQgJHtleHBlY3RlZH0sIGdvdCAke2FjdHVhbH0gKHRvbGVyYW5jZSA9ICR7dG9sZXJhbmNlfSlgKTtcbn1cblxuZXhwb3J0IGZ1bmN0aW9uIGV4cGVjdFRhYmxlKGFjdHVhbDogX0RHLkRhdGFGcmFtZSwgZXhwZWN0ZWQ6IF9ERy5EYXRhRnJhbWUsIGVycm9yPzogc3RyaW5nKTogdm9pZCB7XG4gIGNvbnN0IGV4cGVjdGVkUm93Q291bnQgPSBleHBlY3RlZC5yb3dDb3VudDtcbiAgY29uc3QgYWN0dWFsUm93Q291bnQgPSBhY3R1YWwucm93Q291bnQ7XG4gIGV4cGVjdChhY3R1YWxSb3dDb3VudCwgZXhwZWN0ZWRSb3dDb3VudCwgYCR7ZXJyb3IgPz8gJyd9LCByb3cgY291bnRgKTtcblxuICBmb3IgKGNvbnN0IGNvbHVtbiBvZiBleHBlY3RlZC5jb2x1bW5zKSB7XG4gICAgY29uc3QgYWN0dWFsQ29sdW1uID0gYWN0dWFsLmNvbHVtbnMuYnlOYW1lKGNvbHVtbi5uYW1lKTtcbiAgICBpZiAoYWN0dWFsQ29sdW1uID09IG51bGwpXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYENvbHVtbiAke2NvbHVtbi5uYW1lfSBub3QgZm91bmRgKTtcbiAgICBpZiAoYWN0dWFsQ29sdW1uLnR5cGUgIT0gY29sdW1uLnR5cGUpXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYENvbHVtbiAke2NvbHVtbi5uYW1lfSB0eXBlIGV4cGVjdGVkICR7Y29sdW1uLnR5cGV9IGdvdCAke2FjdHVhbENvbHVtbi50eXBlfWApO1xuICAgIGZvciAobGV0IGkgPSAwOyBpIDwgZXhwZWN0ZWRSb3dDb3VudDsgaSsrKSB7XG4gICAgICBjb25zdCB2YWx1ZSA9IGNvbHVtbi5nZXQoaSk7XG4gICAgICBjb25zdCBhY3R1YWxWYWx1ZSA9IGFjdHVhbENvbHVtbi5nZXQoaSk7XG4gICAgICBpZiAoY29sdW1uLnR5cGUgPT0gREcuVFlQRS5GTE9BVClcbiAgICAgICAgZXhwZWN0RmxvYXQoYWN0dWFsVmFsdWUsIHZhbHVlLCAwLjAwMDEsIGVycm9yKTtcbiAgICAgIGVsc2UgaWYgKGNvbHVtbi50eXBlID09IERHLlRZUEUuREFURV9USU1FKVxuICAgICAgICBleHBlY3QoYWN0dWFsVmFsdWUuaXNTYW1lKHZhbHVlKSwgdHJ1ZSwgZXJyb3IpO1xuICAgICAgZWxzZVxuICAgICAgICBleHBlY3QoYWN0dWFsVmFsdWUsIHZhbHVlLCBlcnJvcik7XG4gICAgfVxuICB9XG59XG5cbmV4cG9ydCBmdW5jdGlvbiBleHBlY3RPYmplY3QoYWN0dWFsOiB7IFtrZXk6IHN0cmluZ106IGFueSB9LCBleHBlY3RlZDogeyBba2V5OiBzdHJpbmddOiBhbnkgfSkge1xuICBmb3IgKGNvbnN0IFtleHBlY3RlZEtleSwgZXhwZWN0ZWRWYWx1ZV0gb2YgT2JqZWN0LmVudHJpZXMoZXhwZWN0ZWQpKSB7XG4gICAgaWYgKCFhY3R1YWwuaGFzT3duUHJvcGVydHkoZXhwZWN0ZWRLZXkpKVxuICAgICAgdGhyb3cgbmV3IEVycm9yKGBFeHBlY3RlZCBwcm9wZXJ0eSBcIiR7ZXhwZWN0ZWRLZXl9XCIgbm90IGZvdW5kYCk7XG5cbiAgICBjb25zdCBhY3R1YWxWYWx1ZSA9IGFjdHVhbFtleHBlY3RlZEtleV07XG4gICAgaWYgKGFjdHVhbFZhbHVlIGluc3RhbmNlb2YgQXJyYXkgJiYgZXhwZWN0ZWRWYWx1ZSBpbnN0YW5jZW9mIEFycmF5KVxuICAgICAgZXhwZWN0QXJyYXkoYWN0dWFsVmFsdWUsIGV4cGVjdGVkVmFsdWUpO1xuICAgIGVsc2UgaWYgKGFjdHVhbFZhbHVlIGluc3RhbmNlb2YgT2JqZWN0ICYmIGV4cGVjdGVkVmFsdWUgaW5zdGFuY2VvZiBPYmplY3QpXG4gICAgICBleHBlY3RPYmplY3QoYWN0dWFsVmFsdWUsIGV4cGVjdGVkVmFsdWUpO1xuICAgIGVsc2UgaWYgKE51bWJlci5pc0Zpbml0ZShhY3R1YWxWYWx1ZSkgJiYgTnVtYmVyLmlzRmluaXRlKGV4cGVjdGVkVmFsdWUpKVxuICAgICAgZXhwZWN0RmxvYXQoYWN0dWFsVmFsdWUsIGV4cGVjdGVkVmFsdWUpO1xuICAgIGVsc2UgaWYgKGFjdHVhbFZhbHVlICE9IGV4cGVjdGVkVmFsdWUpXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYEV4cGVjdGVkICgke2V4cGVjdGVkVmFsdWV9KSBmb3Iga2V5ICcke2V4cGVjdGVkS2V5fScsIGdvdCAoJHthY3R1YWxWYWx1ZX0pYCk7XG4gIH1cbn1cblxuZXhwb3J0IGZ1bmN0aW9uIGV4cGVjdEFycmF5KGFjdHVhbDogQXJyYXlMaWtlPGFueT4sIGV4cGVjdGVkOiBBcnJheUxpa2U8YW55Pikge1xuICBjb25zdCBhY3R1YWxMZW5ndGggPSBhY3R1YWwubGVuZ3RoO1xuICBjb25zdCBleHBlY3RlZExlbmd0aCA9IGV4cGVjdGVkLmxlbmd0aDtcblxuICBpZiAoYWN0dWFsTGVuZ3RoICE9IGV4cGVjdGVkTGVuZ3RoKSB7XG4gICAgdGhyb3cgbmV3IEVycm9yKGBBcnJheXMgYXJlIG9mIGRpZmZlcmVudCBsZW5ndGg6IGFjdHVhbCBhcnJheSBsZW5ndGggaXMgJHthY3R1YWxMZW5ndGh9IGAgK1xuICAgICAgYGFuZCBleHBlY3RlZCBhcnJheSBsZW5ndGggaXMgJHtleHBlY3RlZExlbmd0aH1gKTtcbiAgfVxuXG4gIGZvciAobGV0IGkgPSAwOyBpIDwgYWN0dWFsTGVuZ3RoOyBpKyspIHtcbiAgICBpZiAoYWN0dWFsW2ldIGluc3RhbmNlb2YgQXJyYXkgJiYgZXhwZWN0ZWRbaV0gaW5zdGFuY2VvZiBBcnJheSlcbiAgICAgIGV4cGVjdEFycmF5KGFjdHVhbFtpXSwgZXhwZWN0ZWRbaV0pO1xuICAgIGVsc2UgaWYgKGFjdHVhbFtpXSBpbnN0YW5jZW9mIE9iamVjdCAmJiBleHBlY3RlZFtpXSBpbnN0YW5jZW9mIE9iamVjdClcbiAgICAgIGV4cGVjdE9iamVjdChhY3R1YWxbaV0sIGV4cGVjdGVkW2ldKTtcbiAgICBlbHNlIGlmIChhY3R1YWxbaV0gIT0gZXhwZWN0ZWRbaV0pXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYEV4cGVjdGVkICR7ZXhwZWN0ZWRbaV19IGF0IHBvc2l0aW9uICR7aX0sIGdvdCAke2FjdHVhbFtpXX1gKTtcbiAgfVxufVxuXG4vKiBEZWZpbmVzIGEgdGVzdCBzdWl0ZS4gKi9cbmV4cG9ydCBmdW5jdGlvbiBjYXRlZ29yeShjYXRlZ29yeTogc3RyaW5nLCB0ZXN0c186ICgpID0+IHZvaWQsIG9wdGlvbnM/OiBDYXRlZ29yeU9wdGlvbnMpOiB2b2lkIHtcbiAgY3VycmVudENhdGVnb3J5ID0gY2F0ZWdvcnk7XG4gIHRlc3RzXygpO1xuICBpZiAodGVzdHNbY3VycmVudENhdGVnb3J5XSkge1xuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0uY2xlYXIgPSBvcHRpb25zPy5jbGVhciA/PyB0cnVlO1xuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0udGltZW91dCA9IG9wdGlvbnM/LnRpbWVvdXQ7XG4gICAgdGVzdHNbY3VycmVudENhdGVnb3J5XS5iZW5jaG1hcmtzID0gb3B0aW9ucz8uYmVuY2htYXJrcztcbiAgICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldLnN0cmVzc1Rlc3RzID0gb3B0aW9ucz8uc3RyZXNzVGVzdHM7XG4gICAgdGVzdHNbY3VycmVudENhdGVnb3J5XS5vd25lciA9IG9wdGlvbnM/Lm93bmVyO1xuICB9XG59XG5cbi8qIERlZmluZXMgYSBmdW5jdGlvbiB0byBiZSBleGVjdXRlZCBiZWZvcmUgdGhlIHRlc3RzIGluIHRoaXMgY2F0ZWdvcnkgYXJlIGV4ZWN1dGVkLiAqL1xuZXhwb3J0IGZ1bmN0aW9uIGJlZm9yZShiZWZvcmU6ICgpID0+IFByb21pc2U8dm9pZD4pOiB2b2lkIHtcbiAgaWYgKHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0gPT0gdW5kZWZpbmVkKVxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0gPSB7fTtcbiAgdGVzdHNbY3VycmVudENhdGVnb3J5XS5iZWZvcmUgPSBiZWZvcmU7XG59XG5cbi8qIERlZmluZXMgYSBmdW5jdGlvbiB0byBiZSBleGVjdXRlZCBhZnRlciB0aGUgdGVzdHMgaW4gdGhpcyBjYXRlZ29yeSBhcmUgZXhlY3V0ZWQuICovXG5leHBvcnQgZnVuY3Rpb24gYWZ0ZXIoYWZ0ZXI6ICgpID0+IFByb21pc2U8dm9pZD4pOiB2b2lkIHtcbiAgaWYgKHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0gPT0gdW5kZWZpbmVkKVxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0gPSB7fTtcbiAgdGVzdHNbY3VycmVudENhdGVnb3J5XS5hZnRlciA9IGFmdGVyO1xufVxuXG5mdW5jdGlvbiBhZGROYW1lc3BhY2Uoczogc3RyaW5nLCBmOiBfREcuRnVuYyk6IHN0cmluZyB7XG4gIHJldHVybiBzLnJlcGxhY2UobmV3IFJlZ0V4cChmLm5hbWUsICdnaScpLCBmLm5xTmFtZSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBpbml0QXV0b1Rlc3RzKHBhY2thZ2VfOiBfREcuUGFja2FnZSwgbW9kdWxlPzogYW55KSB7XG4gIGNvbnN0IHBhY2thZ2VJZCA9IHBhY2thZ2VfLmlkO1xuICBpZiAod2FzUmVnaXN0ZXJlZFtwYWNrYWdlSWRdKSByZXR1cm47XG4gIGNvbnN0IG1vZHVsZVRlc3RzID0gbW9kdWxlID8gbW9kdWxlLnRlc3RzIDogdGVzdHM7XG4gIGlmIChwYWNrYWdlXy5uYW1lID09PSAnRGV2VG9vbHMnIHx8ICghIW1vZHVsZSAmJiBtb2R1bGUuX3BhY2thZ2UubmFtZSA9PT0gJ0RldlRvb2xzJykpIHtcbiAgICBmb3IgKGNvbnN0IGYgb2YgKDxhbnk+d2luZG93KS5kYXJ0VGVzdHMpIHtcbiAgICAgIGNvbnN0IGFyciA9IGYubmFtZS5zcGxpdCgvXFxzKlxcfFxccyohL2cpO1xuICAgICAgbGV0IG5hbWUgPSBhcnIucG9wKCkgPz8gZi5uYW1lO1xuICAgICAgbGV0IGNhdCA9IGFyci5sZW5ndGggPyBjb3JlQ2F0TmFtZSArICc6ICcgKyBhcnIuam9pbignOiAnKSA6IGNvcmVDYXROYW1lO1xuICAgICAgbGV0IGZ1bGxOYW1lOiBzdHJpbmdbXSA9IG5hbWUuc3BsaXQoJyB8ICcpO1xuICAgICAgbmFtZSA9IGZ1bGxOYW1lW2Z1bGxOYW1lLmxlbmd0aCAtIDFdO1xuICAgICAgZnVsbE5hbWUudW5zaGlmdChjYXQpO1xuICAgICAgZnVsbE5hbWUucG9wKCk7XG4gICAgICBjYXQgPSBmdWxsTmFtZS5qb2luKCc6ICcpO1xuICAgICAgaWYgKG1vZHVsZVRlc3RzW2NhdF0gPT09IHVuZGVmaW5lZClcbiAgICAgICAgbW9kdWxlVGVzdHNbY2F0XSA9IHsgdGVzdHM6IFtdLCBjbGVhcjogdHJ1ZSB9O1xuICAgICAgbW9kdWxlVGVzdHNbY2F0XS50ZXN0cy5wdXNoKG5ldyBUZXN0KGNhdCwgbmFtZSwgZi50ZXN0LCB7IGlzQWdncmVnYXRlZDogZmFsc2UsIHRpbWVvdXQ6IGYub3B0aW9ucz8udGltZW91dCA/PyBTVEFOREFSVF9USU1FT1VULCBza2lwUmVhc29uOiBmLm9wdGlvbnM/LnNraXBSZWFzb24sIG93bmVyOiBmLm9wdGlvbnM/Lm93bmVyLCBiZW5jaG1hcms6IGYub3B0aW9ucz8uYmVuY2htYXJrID8/IGZhbHNlIH0pKTtcbiAgICB9XG4gIH1cbiAgY29uc3QgbW9kdWxlQXV0b1Rlc3RzID0gW107XG4gIGNvbnN0IG1vZHVsZURlbW8gPSBbXTtcbiAgY29uc3QgbW9kdWxlRGV0ZWN0b3JzID0gW107XG4gIGNvbnN0IHBhY2tGdW5jdGlvbnMgPSBhd2FpdCBncm9rLmRhcGkuZnVuY3Rpb25zLmZpbHRlcihgcGFja2FnZS5pZCA9IFwiJHtwYWNrYWdlSWR9XCJgKS5saXN0KCk7XG4gIGNvbnN0IHJlZyA9IG5ldyBSZWdFeHAoL3NraXA6XFxzKihbXixcXHNdKyl8d2FpdDpcXHMqKFxcZCspfGNhdDpcXHMqKFteLFxcc10rKXx0aW1lb3V0OlxccyooXFxkKykvZyk7XG4gIGZvciAoY29uc3QgZiBvZiBwYWNrRnVuY3Rpb25zKSB7XG4gICAgY29uc3QgdGVzdHMgPSBmLm9wdGlvbnNbJ3Rlc3QnXTtcbiAgICBjb25zdCBkZW1vID0gZi5vcHRpb25zWydkZW1vUGF0aCddO1xuICAgIGlmICgodGVzdHMgJiYgQXJyYXkuaXNBcnJheSh0ZXN0cykgJiYgdGVzdHMubGVuZ3RoKSkge1xuICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCB0ZXN0cy5sZW5ndGg7IGkrKykge1xuICAgICAgICBjb25zdCByZXMgPSAodGVzdHNbaV0gYXMgc3RyaW5nKS5tYXRjaEFsbChyZWcpO1xuICAgICAgICBjb25zdCBtYXA6IHsgc2tpcD86IHN0cmluZywgd2FpdD86IG51bWJlciwgY2F0Pzogc3RyaW5nLCB0aW1lb3V0PzogbnVtYmVyLCBiZW5jaG1hcmtUaW1lb3V0PzogbnVtYmVyIH0gPSB7fTtcbiAgICAgICAgQXJyYXkuZnJvbShyZXMpLmZvckVhY2goKGFycikgPT4ge1xuICAgICAgICAgIGlmIChhcnJbMF0uc3RhcnRzV2l0aCgnc2tpcCcpKSBtYXBbJ3NraXAnXSA9IGFyclsxXTtcbiAgICAgICAgICBlbHNlIGlmIChhcnJbMF0uc3RhcnRzV2l0aCgnd2FpdCcpKSBtYXBbJ3dhaXQnXSA9IHBhcnNlSW50KGFyclsyXSk7XG4gICAgICAgICAgZWxzZSBpZiAoYXJyWzBdLnN0YXJ0c1dpdGgoJ2NhdCcpKSBtYXBbJ2NhdCddID0gYXJyWzNdO1xuICAgICAgICAgIGVsc2UgaWYgKGFyclswXS5zdGFydHNXaXRoKCd0aW1lb3V0JykpIG1hcFsndGltZW91dCddID0gcGFyc2VJbnQoYXJyWzRdKTtcbiAgICAgICAgfSk7XG4gICAgICAgIGNvbnN0IHRlc3QgPSBuZXcgVGVzdChtYXAuY2F0ID8/IGF1dG9UZXN0c0NhdE5hbWUsIHRlc3RzLmxlbmd0aCA9PT0gMSA/IGYubmFtZSA6IGAke2YubmFtZX0gJHtpICsgMX1gLCBhc3luYyAoKSA9PiB7XG4gICAgICAgICAgY29uc3QgcmVzID0gYXdhaXQgZ3Jvay5mdW5jdGlvbnMuZXZhbChhZGROYW1lc3BhY2UodGVzdHNbaV0sIGYpKTtcbiAgICAgICAgICBpZiAobWFwLndhaXQpIGF3YWl0IGRlbGF5KG1hcC53YWl0KTtcbiAgICAgICAgICAvLyBlc2xpbnQtZGlzYWJsZS1uZXh0LWxpbmUgbm8tdGhyb3ctbGl0ZXJhbFxuICAgICAgICAgIGlmICh0eXBlb2YgcmVzID09PSAnYm9vbGVhbicgJiYgIXJlcykgdGhyb3cgYEZhaWxlZDogJHt0ZXN0c1tpXX0sIGV4cGVjdGVkIHRydWUsIGdvdCAke3Jlc31gO1xuICAgICAgICB9LCB7IHNraXBSZWFzb246IG1hcC5za2lwLCB0aW1lb3V0OiBERy5UZXN0LmlzSW5CZW5jaG1hcmsgPyBtYXAuYmVuY2htYXJrVGltZW91dCA/PyBCRU5DSE1BUktfVElNRU9VVCA6IG1hcC50aW1lb3V0ID8/IFNUQU5EQVJUX1RJTUVPVVQgfSk7XG4gICAgICAgIGlmIChtYXAuY2F0KSB7XG4gICAgICAgICAgY29uc3QgY2F0OiBzdHJpbmcgPSBtYXAuY2F0O1xuICAgICAgICAgIGlmIChtb2R1bGVUZXN0c1tjYXRdID09PSB1bmRlZmluZWQpXG4gICAgICAgICAgICBtb2R1bGVUZXN0c1tjYXRdID0geyB0ZXN0czogW10sIGNsZWFyOiB0cnVlIH07XG5cbiAgICAgICAgICAvLyBvbmx5IGJlZm9yZS9hZnRlciBjYW4gYmUgZGVmaW5lZCBpbiB0cyBmaWxlcyB0ZXN0cyB1bmRlciB0aGUgY2F0ZWdvcnlcbiAgICAgICAgICBpZiAoIW1vZHVsZVRlc3RzW2NhdF0udGVzdHMpXG4gICAgICAgICAgICBtb2R1bGVUZXN0c1tjYXRdLnRlc3RzID0gW107XG4gICAgICAgICAgbW9kdWxlVGVzdHNbY2F0XS50ZXN0cy5wdXNoKHRlc3QpO1xuICAgICAgICB9XG4gICAgICAgIGVsc2VcbiAgICAgICAgICBtb2R1bGVBdXRvVGVzdHMucHVzaCh0ZXN0KTtcbiAgICAgIH1cbiAgICB9XG4gICAgaWYgKGRlbW8pIHtcbiAgICAgIGNvbnN0IHdhaXQgPSBmLm9wdGlvbnNbJ2RlbW9XYWl0J10gPyBwYXJzZUludChmLm9wdGlvbnNbJ2RlbW9XYWl0J10pIDogdW5kZWZpbmVkO1xuICAgICAgY29uc3QgdGVzdCA9IG5ldyBUZXN0KGRlbW9DYXROYW1lLCBmLmZyaWVuZGx5TmFtZSwgYXN5bmMgKCkgPT4ge1xuICAgICAgICBhd2FpdCBkZWxheSgzMDApO1xuICAgICAgICBncm9rLnNoZWxsLmNsZWFyTGFzdEVycm9yKCk7XG4gICAgICAgIGF3YWl0IGYuYXBwbHkoKTtcbiAgICAgICAgYXdhaXQgZGVsYXkod2FpdCA/IHdhaXQgOiAyMDAwKTtcbiAgICAgICAgY29uc3QgdW5oYW5kbGVkID0gYXdhaXQgZ3Jvay5zaGVsbC5sYXN0RXJyb3I7XG4gICAgICAgIGlmICh1bmhhbmRsZWQpXG4gICAgICAgICAgdGhyb3cgbmV3IEVycm9yKHVuaGFuZGxlZCk7XG4gICAgICB9LCB7IHNraXBSZWFzb246IGYub3B0aW9uc1snZGVtb1NraXAnXSB9KTtcbiAgICAgIG1vZHVsZURlbW8ucHVzaCh0ZXN0KTtcbiAgICB9XG4gICAgaWYgKGYuaGFzVGFnKCdzZW1UeXBlRGV0ZWN0b3InKSkge1xuICAgICAgbGV0IGRldGVjdG9yc1Rlc3REYXRhID0gdGVzdERhdGE7XG4gICAgICBpZiAoZi5vcHRpb25zWyd0ZXN0RGF0YSddKSB7XG4gICAgICAgIGRldGVjdG9yc1Rlc3REYXRhID0gYXdhaXQgZ3Jvay5kYXRhLmZpbGVzLm9wZW5UYWJsZShgU3lzdGVtOkFwcERhdGEvJHtwYWNrYWdlXy5ucU5hbWV9LyR7Zi5vcHRpb25zWyd0ZXN0RGF0YSddfWApO1xuICAgICAgfVxuXG4gICAgICBjb25zdCB0ZXN0ID0gbmV3IFRlc3QoZGV0ZWN0b3JzQ2F0TmFtZSwgZi5mcmllbmRseU5hbWUsIGFzeW5jICgpID0+IHtcbiAgICAgICAgY29uc3QgYXJyID0gW107XG4gICAgICAgIGNvbnNvbGUubG9nKGBTeXN0ZW06QXBwRGF0YS8ke3BhY2thZ2VfLm5xTmFtZX0vJHtmLm9wdGlvbnNbJ3Rlc3REYXRhJ119YCk7XG5cbiAgICAgICAgZm9yIChjb25zdCBjb2wgb2YgZGV0ZWN0b3JzVGVzdERhdGEuY2xvbmUoKS5jb2x1bW5zKSB7XG4gICAgICAgICAgY29uc3QgcmVzID0gYXdhaXQgZi5hcHBseShbY29sXSk7XG4gICAgICAgICAgYXJyLnB1c2gocmVzIHx8IGNvbC5zZW1UeXBlKTtcbiAgICAgICAgfVxuICAgICAgICBjb25zdCByZXNBcnIgPSBhcnIuZmlsdGVyKChpKSA9PiBpKTtcbiAgICAgICAgZXhwZWN0KHJlc0Fyci5sZW5ndGgsIDEpO1xuXG4gICAgICAgIGlmIChmLm9wdGlvbnNbJ3Rlc3REYXRhQ29sdW1uTmFtZSddKVxuICAgICAgICAgIGV4cGVjdChyZXNBcnJbMF0sIGYub3B0aW9uc1sndGVzdERhdGFDb2x1bW5OYW1lJ10pO1xuXG4gICAgICB9LCB7IHNraXBSZWFzb246IGYub3B0aW9uc1snc2tpcFRlc3QnXSB9KTtcbiAgICAgIG1vZHVsZURldGVjdG9ycy5wdXNoKHRlc3QpO1xuICAgIH1cbiAgfVxuICB3YXNSZWdpc3RlcmVkW3BhY2thZ2VJZF0gPSB0cnVlO1xuICBpZiAobW9kdWxlQXV0b1Rlc3RzLmxlbmd0aCA+IDApXG4gICAgbW9kdWxlVGVzdHNbYXV0b1Rlc3RzQ2F0TmFtZV0gPSB7IHRlc3RzOiBtb2R1bGVBdXRvVGVzdHMsIGNsZWFyOiB0cnVlIH07XG4gIGlmIChtb2R1bGVEZW1vLmxlbmd0aCA+IDApXG4gICAgbW9kdWxlVGVzdHNbZGVtb0NhdE5hbWVdID0geyB0ZXN0czogbW9kdWxlRGVtbywgY2xlYXI6IHRydWUgfTtcbiAgaWYgKG1vZHVsZURldGVjdG9ycy5sZW5ndGggPiAwKVxuICAgIG1vZHVsZVRlc3RzW2RldGVjdG9yc0NhdE5hbWVdID0geyB0ZXN0czogbW9kdWxlRGV0ZWN0b3JzLCBjbGVhcjogZmFsc2UgfTtcbn1cblxuZnVuY3Rpb24gcmVkZWZpbmVDb25zb2xlKCk6IGFueVtdIHtcbiAgY29uc3QgbG9nczogYW55W10gPSBbXTtcbiAgY29uc29sZS5sb2cgPSAoLi4uYXJncykgPT4ge1xuICAgIGxvZ3MucHVzaCguLi5hcmdzKTtcbiAgICBzdGRMb2coLi4uYXJncyk7XG4gIH07XG4gIGNvbnNvbGUuaW5mbyA9ICguLi5hcmdzKSA9PiB7XG4gICAgbG9ncy5wdXNoKC4uLmFyZ3MpO1xuICAgIHN0ZEluZm8oLi4uYXJncyk7XG4gIH07XG4gIGNvbnNvbGUud2FybiA9ICguLi5hcmdzKSA9PiB7XG4gICAgbG9ncy5wdXNoKC4uLmFyZ3MpO1xuICAgIHN0ZFdhcm4oLi4uYXJncyk7XG4gIH07XG4gIGNvbnNvbGUuZXJyb3IgPSAoLi4uYXJncykgPT4ge1xuICAgIGxvZ3MucHVzaCguLi5hcmdzKTtcbiAgICBzdGRFcnJvciguLi5hcmdzKTtcbiAgfTtcbiAgcmV0dXJuIGxvZ3M7XG59XG5cbmZ1bmN0aW9uIHJlc2V0Q29uc29sZSgpOiB2b2lkIHtcbiAgY29uc29sZS5sb2cgPSBzdGRMb2c7XG4gIGNvbnNvbGUuaW5mbyA9IHN0ZEluZm87XG4gIGNvbnNvbGUud2FybiA9IHN0ZFdhcm47XG4gIGNvbnNvbGUuZXJyb3IgPSBzdGRFcnJvcjtcbn1cblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHJ1blRlc3RzKG9wdGlvbnM/OiBUZXN0RXhlY3V0aW9uT3B0aW9ucykgOiBQcm9taXNlPFRlc3RSZXN1bHRFeHRlbmRlZFtdPntcblxuICBjb25zdCBwYWNrYWdlXzogX0RHLlBhY2thZ2UgPSBvcHRpb25zPy5ub2RlT3B0aW9ucyA/IG9wdGlvbnMubm9kZU9wdGlvbnMucGFja2FnZSA6IGdyb2suZnVuY3Rpb25zLmdldEN1cnJlbnRDYWxsKCkuZnVuYy5wYWNrYWdlO1xuICBpZiAoIXBhY2thZ2VfKVxuICAgIHRocm93IG5ldyBFcnJvcignQ2FuXFwndCBydW4gdGVzdHMgb3V0c2lkZSBvZiB0aGUgcGFja2FnZScpO1xuICBjb25zdCBtYXRjaCA9IHBhY2thZ2VfLnBhY2thZ2VPd25lcj8ubWF0Y2goLzwoW14+XSopPi8pO1xuICBjb25zdCBwYWNrYWdlT3duZXIgPSBtYXRjaCA/IG1hdGNoWzFdIDogJyc7XG4gIGlmIChwYWNrYWdlXyAhPSB1bmRlZmluZWQpXG4gICAgYXdhaXQgaW5pdEF1dG9UZXN0cyhwYWNrYWdlXyk7XG4gIGNvbnN0IHJlc3VsdHM6VGVzdFJlc3VsdEV4dGVuZGVkW10gPSBbXTtcbiAgY29uc29sZS5sb2coYFJ1bm5pbmcgdGVzdHMuLi5gKTtcbiAgY29uc29sZS5sb2cob3B0aW9ucyk7XG4gIG9wdGlvbnMgPz89IHt9O1xuICBvcHRpb25zIS50ZXN0Q29udGV4dCA/Pz0gbmV3IFRlc3RDb250ZXh0KCk7XG4gIGdyb2suc2hlbGwuY2xlYXJMYXN0RXJyb3IoKTtcbiAgY29uc3QgbG9ncyA9IHJlZGVmaW5lQ29uc29sZSgpO1xuXG4gIGF3YWl0IGludm9rZVRlc3RzKHRlc3RzLCBvcHRpb25zKTtcblxuICBmb3IgKGxldCByIG9mIHJlc3VsdHMpIHtcbiAgICByLnJlc3VsdCA9IHIucmVzdWx0LnRvU3RyaW5nKCkucmVwbGFjZSgvXCIvZywgJ1xcJycpO1xuICAgIGlmIChyLmxvZ3MgIT0gdW5kZWZpbmVkKVxuICAgICAgci5sb2dzID0gci5sb2dzIS50b1N0cmluZygpLnJlcGxhY2UoL1wiL2csICdcXCcnKTtcbiAgfVxuICByZXR1cm4gcmVzdWx0cztcblxuICBhc3luYyBmdW5jdGlvbiBpbnZva2VDYXRlZ29yeU1ldGhvZChtZXRob2Q6ICgoKSA9PiBQcm9taXNlPHZvaWQ+KSB8IHVuZGVmaW5lZCwgY2F0ZWdvcnk6IHN0cmluZyk6IFByb21pc2U8c3RyaW5nIHwgdW5kZWZpbmVkPiB7XG4gICAgbGV0IGludm9rYXRpb25SZXN1bHQgPSB1bmRlZmluZWQ7XG4gICAgdHJ5IHtcbiAgICAgIGlmIChtZXRob2QgIT09IHVuZGVmaW5lZCkge1xuICAgICAgICBhd2FpdCB0aW1lb3V0KGFzeW5jICgpID0+IHtcbiAgICAgICAgICBhd2FpdCBtZXRob2QoKTtcbiAgICAgICAgfSwgMTAwMDAwLCBgYmVmb3JlICR7Y2F0ZWdvcnl9OiB0aW1lb3V0IGVycm9yYCk7XG4gICAgICB9XG4gICAgfSBjYXRjaCAoeDogYW55KSB7XG4gICAgICBpbnZva2F0aW9uUmVzdWx0ID0gYXdhaXQgZ2V0UmVzdWx0KHgpO1xuICAgIH1cbiAgICByZXR1cm4gaW52b2thdGlvblJlc3VsdFxuICB9XG5cbiAgYXN5bmMgZnVuY3Rpb24gaW52b2tlVGVzdHNJbkNhdGVnb3J5KGNhdGVnb3J5OiBDYXRlZ29yeSwgb3B0aW9uczogVGVzdEV4ZWN1dGlvbk9wdGlvbnMsIGlzVGFyZ2V0Q2F0ZWdvcnk6IGJvb2xlYW4pOiBQcm9taXNlPFRlc3RSZXN1bHRFeHRlbmRlZFtdPiB7XG4gICAgbGV0IHQgPSBjYXRlZ29yeS50ZXN0cyA/PyBbXTtcbiAgICBjb25zdCByZXMgOiBUZXN0UmVzdWx0RXh0ZW5kZWRbXSA9IFtdO1xuICAgIC8vIGxldCBtZW1vcnlVc2FnZUJlZm9yZSA9ICh3aW5kb3c/LnBlcmZvcm1hbmNlIGFzIGFueSk/Lm1lbW9yeT8udXNlZEpTSGVhcFNpemU7XG4gICAgY29uc3Qgd2lkZ2V0c0JlZm9yZSA9IGdldFdpZGdldHNDb3VudFNhZmUoKTtcblxuICAgIGlmIChjYXRlZ29yeS5jbGVhcikge1xuICAgICAgICBsZXQgc2tpcHBpbmdUZXN0cyA9IGlzVGFyZ2V0Q2F0ZWdvcnkgJiYgb3B0aW9ucy5za2lwVG9UZXN0ICE9IHVuZGVmaW5lZDtcbiAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgdC5sZW5ndGg7IGkrKykge1xuXG4gICAgICAgIGlmICh0W2ldLm9wdGlvbnMpIHtcbiAgICAgICAgICBpZiAodFtpXS5vcHRpb25zPy5iZW5jaG1hcmsgPT09IHVuZGVmaW5lZCkge1xuICAgICAgICAgICAgaWYgKCF0W2ldLm9wdGlvbnMpXG4gICAgICAgICAgICAgIHRbaV0ub3B0aW9ucyA9IHt9XG4gICAgICAgICAgICB0W2ldLm9wdGlvbnMhLmJlbmNobWFyayA9IGNhdGVnb3J5LmJlbmNobWFya3MgPz8gZmFsc2U7XG4gICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIGxldCB0ZXN0ID0gdFtpXTtcbiAgICAgICAgaWYgKG9wdGlvbnMudGVzdClcbiAgICAgICAgICBpZiAob3B0aW9ucy50ZXN0LnRvTG93ZXJDYXNlKCkgIT09IHRlc3QubmFtZS50b0xvd2VyQ2FzZSgpKVxuICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgIGlmIChza2lwcGluZ1Rlc3RzKSB7XG4gICAgICAgICAgaWYgKG9wdGlvbnM/LnNraXBUb1Rlc3QgIT0gdW5kZWZpbmVkICYmIHRlc3QubmFtZS50b0xvd2VyQ2FzZSgpLnRyaW0oKSA9PT0gb3B0aW9ucz8uc2tpcFRvVGVzdC50b0xvd2VyQ2FzZSgpLnRyaW0oKSkge1xuICAgICAgICAgICAgLy8gRm91bmQgdGhlIHRhcmdldCB0ZXN0LCBzdG9wIHNraXBwaW5nIGFmdGVyIHRoaXMgb25lXG4gICAgICAgICAgICBza2lwcGluZ1Rlc3RzID0gZmFsc2U7XG4gICAgICAgICAgfSBlbHNlXG4gICAgICAgICAgY29udGludWU7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKHRlc3Q/Lm9wdGlvbnMpIHtcbiAgICAgICAgICB0ZXN0Lm9wdGlvbnMub3duZXIgPSB0W2ldLm9wdGlvbnM/Lm93bmVyID8/IGNhdGVnb3J5Py5vd25lciA/PyBwYWNrYWdlT3duZXIgPz8gJyc7XG4gICAgICAgIH1cbiAgICAgICAgLy8gbGV0IGlzR0JFbmFibGUgPSAod2luZG93IGFzIGFueSkuZ2MgJiYgdGVzdC5vcHRpb25zPy5za2lwUmVhc29uID09IHVuZGVmaW5lZDtcbiAgICAgICAgLy8gY29uc29sZS5sb2coYCoqKioqKioqJHtpc0dCRW5hYmxlfWApO1xuICAgICAgICAvLyBpZiAoaXNHQkVuYWJsZSlcbiAgICAgICAgLy8gICBhd2FpdCAod2luZG93IGFzIGFueSkuZ2MoKTtcbiAgICAgICAgLy8gbWVtb3J5VXNhZ2VCZWZvcmUgPSAod2luZG93Py5wZXJmb3JtYW5jZSBhcyBhbnkpPy5tZW1vcnk/LnVzZWRKU0hlYXBTaXplO1xuICAgICAgICBsZXQgdGVzdFJ1biA9IGF3YWl0IGV4ZWNUZXN0KFxuICAgICAgICAgICAgdGVzdCxcbiAgICAgICAgICAgIG9wdGlvbnM/LnRlc3QsXG4gICAgICAgICAgICBsb2dzLCBERy5UZXN0LmlzSW5CZW5jaG1hcmsgPyB0W2ldLm9wdGlvbnM/LmJlbmNobWFya1RpbWVvdXQgPz8gQkVOQ0hNQVJLX1RJTUVPVVQgOiB0W2ldLm9wdGlvbnM/LnRpbWVvdXQgPz8gU1RBTkRBUlRfVElNRU9VVCxcbiAgICAgICAgICAgIHBhY2thZ2VfLm5hbWUsXG4gICAgICAgICAgICBvcHRpb25zLnZlcmJvc2VcbiAgICAgICAgKTtcblxuICAgICAgICAvLyBpZiAoaXNHQkVuYWJsZSlcbiAgICAgICAgLy8gICBhd2FpdCAod2luZG93IGFzIGFueSkuZ2MoKTtcbiAgICAgICAgaWYgKHRlc3RSdW4pIHtcbiAgICAgICAgICByZXMucHVzaCh7IC4uLnRlc3RSdW4sICB3aWRnZXRzRGlmZmVyZW5jZTogZ2V0V2lkZ2V0c0NvdW50U2FmZSgpIC0gd2lkZ2V0c0JlZm9yZSB9KTtcbiAgICAgICAgICAvLyBSZXR1cm4gZWFybHkgaWYgcmV0dXJuT25GYWlsIGlzIHNldCBhbmQgdGVzdCBmYWlsZWQgKGJ1dCBpZ25vcmUgZmFpbHVyZSBmb3IgdGhlIHNraXBUb1Rlc3QgdGVzdCBpdHNlbGYpXG4gICAgICAgICAgaWYgKG9wdGlvbnMucmV0dXJuT25GYWlsICYmIG9wdGlvbnMuc2tpcFRvVGVzdCAhPT0gdGVzdC5uYW1lICYmICF0ZXN0UnVuLnN1Y2Nlc3MgJiYgIXRlc3RSdW4uc2tpcHBlZClcbiAgICAgICAgICAgIHJldHVybiByZXM7XG4gICAgICAgIH1cbiAgICAgICAgLy8gcmVzLnB1c2goeyAuLi50ZXN0UnVuLCBtZW1vcnlEZWx0YTogKHdpbmRvdz8ucGVyZm9ybWFuY2UgYXMgYW55KT8ubWVtb3J5Py51c2VkSlNIZWFwU2l6ZSAtIG1lbW9yeVVzYWdlQmVmb3JlLCB3aWRnZXRzRGVsdGE6IGdldFdpZGdldHNDb3VudFNhZmUoKSAtIHdpZGdldHNCZWZvcmUgfSk7XG5cbiAgICAgICAgaWYgKCFvcHRpb25zLm5vZGVPcHRpb25zKSB7XG4gICAgICAgICAgZ3Jvay5zaGVsbC5jbG9zZUFsbCgpO1xuICAgICAgICAgIERHLkJhbGxvb24uY2xvc2VBbGwoKTtcbiAgICAgICAgfVxuICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICBsZXQgc2tpcHBpbmdUZXN0cyA9IGlzVGFyZ2V0Q2F0ZWdvcnkgJiYgb3B0aW9ucy5za2lwVG9UZXN0ICE9IHVuZGVmaW5lZDtcbiAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgdC5sZW5ndGg7IGkrKykge1xuICAgICAgICBsZXQgdGVzdCA9IHRbaV07XG4gICAgICAgIGlmIChvcHRpb25zLnRlc3QpXG4gICAgICAgICAgaWYgKG9wdGlvbnMudGVzdC50b0xvd2VyQ2FzZSgpICE9PSB0ZXN0Lm5hbWUudG9Mb3dlckNhc2UoKSlcbiAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICBpZiAoc2tpcHBpbmdUZXN0cykge1xuICAgICAgICAgIGlmIChvcHRpb25zPy5za2lwVG9UZXN0ICE9IHVuZGVmaW5lZCAmJiB0ZXN0Lm5hbWUudG9Mb3dlckNhc2UoKS50cmltKCkgPT09IG9wdGlvbnM/LnNraXBUb1Rlc3QudG9Mb3dlckNhc2UoKS50cmltKCkpIHtcbiAgICAgICAgICAgIC8vIEZvdW5kIHRoZSB0YXJnZXQgdGVzdCwgc3RvcCBza2lwcGluZyBhZnRlciB0aGlzIG9uZVxuICAgICAgICAgICAgc2tpcHBpbmdUZXN0cyA9IGZhbHNlO1xuICAgICAgICAgIH1cbiAgICAgICAgICBjb250aW51ZTsgIC8vIFNraXAgdGhpcyB0ZXN0IChpbmNsdWRpbmcgdGhlIHRhcmdldClcbiAgICAgICAgfVxuXG4gICAgICAgIGlmICh0ZXN0Py5vcHRpb25zKSB7XG4gICAgICAgICAgdGVzdC5vcHRpb25zLm93bmVyID0gdFtpXS5vcHRpb25zPy5vd25lciA/PyBjYXRlZ29yeT8ub3duZXIgPz8gcGFja2FnZU93bmVyID8/ICcnO1xuICAgICAgICB9XG4gICAgICAgIC8vIGxldCBpc0dCRW5hYmxlID0gKHdpbmRvdyBhcyBhbnkpLmdjICYmIHRlc3Qub3B0aW9ucz8uc2tpcFJlYXNvbiA9PSB1bmRlZmluZWQ7XG4gICAgICAgIC8vIGNvbnNvbGUubG9nKGAqKioqKioqKiR7aXNHQkVuYWJsZX1gKTtcbiAgICAgICAgLy8gaWYgKGlzR0JFbmFibGUpXG4gICAgICAgIC8vICAgYXdhaXQgKHdpbmRvdyBhcyBhbnkpLmdjKCk7XG4gICAgICAgIC8vIG1lbW9yeVVzYWdlQmVmb3JlID0gKHdpbmRvdz8ucGVyZm9ybWFuY2UgYXMgYW55KT8ubWVtb3J5Py51c2VkSlNIZWFwU2l6ZTtcbiAgICAgICAgbGV0IHRlc3RSdW4gPSBhd2FpdCBleGVjVGVzdChcbiAgICAgICAgICAgIHRlc3QsXG4gICAgICAgICAgICBvcHRpb25zPy50ZXN0LFxuICAgICAgICAgICAgbG9ncyxcbiAgICAgICAgICAgIERHLlRlc3QuaXNJbkJlbmNobWFyayA/IHRbaV0ub3B0aW9ucz8uYmVuY2htYXJrVGltZW91dCA/PyBCRU5DSE1BUktfVElNRU9VVCA6IHRbaV0ub3B0aW9ucz8udGltZW91dCxcbiAgICAgICAgICAgIHBhY2thZ2VfLm5hbWUsXG4gICAgICAgICAgICBvcHRpb25zLnZlcmJvc2VcbiAgICAgICAgKTtcblxuICAgICAgICAvLyBpZiAoaXNHQkVuYWJsZSlcbiAgICAgICAgLy8gICBhd2FpdCAod2luZG93IGFzIGFueSkuZ2MoKTtcblxuICAgICAgICBpZiAodGVzdFJ1bikge1xuICAgICAgICAgIHJlcy5wdXNoKHsgLi4udGVzdFJ1biwgd2lkZ2V0c0RpZmZlcmVuY2U6IGdldFdpZGdldHNDb3VudFNhZmUoKSAtIHdpZGdldHNCZWZvcmUgfSk7XG4gICAgICAgICAgLy8gUmV0dXJuIGVhcmx5IGlmIHJldHVybk9uRmFpbCBpcyBzZXQgYW5kIHRlc3QgZmFpbGVkIChidXQgaWdub3JlIGZhaWx1cmUgZm9yIHRoZSBza2lwVG9UZXN0IHRlc3QgaXRzZWxmKVxuICAgICAgICAgIGlmIChvcHRpb25zLnJldHVybk9uRmFpbCAmJiBvcHRpb25zLnNraXBUb1Rlc3QgIT09IHRlc3QubmFtZSAmJiAhdGVzdFJ1bi5zdWNjZXNzICYmICF0ZXN0UnVuLnNraXBwZWQpXG4gICAgICAgICAgICByZXR1cm4gcmVzO1xuICAgICAgICB9XG4gICAgICAgIC8vIHJlcy5wdXNoKHsgLi4udGVzdFJ1biwgbWVtb3J5RGVsdGE6ICh3aW5kb3c/LnBlcmZvcm1hbmNlIGFzIGFueSk/Lm1lbW9yeT8udXNlZEpTSGVhcFNpemUgLSBtZW1vcnlVc2FnZUJlZm9yZSwgd2lkZ2V0c0RpZmZlcmVuY2U6IGdldFdpZGdldHNDb3VudFNhZmUoKSAtIHdpZGdldHNCZWZvcmUgfSk7XG5cbiAgICAgIH1cbiAgICB9XG4gICAgcmV0dXJuIHJlcztcbiAgfVxuXG4gIGZ1bmN0aW9uIGdldFdpZGdldHNDb3VudFNhZmUoKSB7XG4gICAgaWYgKHR5cGVvZiBwcm9jZXNzICE9PSAndW5kZWZpbmVkJylcbiAgICAgIHJldHVybiAwO1xuICAgIGxldCBsZW5ndGggPSAtMTtcbiAgICB0cnkge1xuICAgICAgbGVuZ3RoID0gREcuV2lkZ2V0LmdldEFsbCgpLmxlbmd0aDtcbiAgICB9IGNhdGNoIChlOiBhbnkpIHtcbiAgICAgIGNvbnNvbGUud2FybihlLm1lc3NhZ2UgPz8gZSk7XG4gICAgfVxuICAgIHJldHVybiBsZW5ndGg7XG4gIH1cblxuICBhc3luYyBmdW5jdGlvbiBpbnZva2VUZXN0cyhjYXRlZ29yaWVzVG9JbnZva2U6IHsgW2tleTogc3RyaW5nXTogQ2F0ZWdvcnkgfSwgb3B0aW9uczogVGVzdEV4ZWN1dGlvbk9wdGlvbnMpIHtcbiAgICB0cnkge1xuICAgICAgbGV0IHNraXBwaW5nQ2F0ZWdvcmllcyA9IG9wdGlvbnM/LnNraXBUb0NhdGVnb3J5ICE9IHVuZGVmaW5lZDtcbiAgICAgIGxldCBpc1RhcmdldENhdGVnb3J5ID0gZmFsc2U7XG4gICAgICBmb3IgKGNvbnN0IFtrZXksIHZhbHVlXSBvZiBPYmplY3QuZW50cmllcyhjYXRlZ29yaWVzVG9JbnZva2UpKSB7XG4gICAgICAgICAgaWYgKG9wdGlvbnMuZXhjbHVkZT8uc29tZSgoYykgPT4ga2V5LnN0YXJ0c1dpdGgoYykpKVxuICAgICAgICAgICAgICBjb250aW51ZTtcbiAgICAgICAgICBpZiAob3B0aW9ucz8uY2F0ZWdvcnkgIT0gbnVsbCAmJiAha2V5LnRvTG93ZXJDYXNlKCkuc3RhcnRzV2l0aChvcHRpb25zPy5jYXRlZ29yeS50b0xvd2VyQ2FzZSgpLnRyaW0oKSkpXG4gICAgICAgICAgICAgIGNvbnRpbnVlO1xuXG4gICAgICAgICAgaWYgKHNraXBwaW5nQ2F0ZWdvcmllcykge1xuICAgICAgICAgICAgICBpZiAoaXNUYXJnZXRDYXRlZ29yeSlcbiAgICAgICAgICAgICAgICAgIHNraXBwaW5nQ2F0ZWdvcmllcyA9IGZhbHNlO1xuICAgICAgICAgICAgICBlbHNlIHtcbiAgICAgICAgICAgICAgICAgIGlmIChvcHRpb25zPy5za2lwVG9DYXRlZ29yeSAhPSBudWxsICYmIGtleS50b0xvd2VyQ2FzZSgpLnRyaW0oKSA9PT0gb3B0aW9ucz8uc2tpcFRvQ2F0ZWdvcnkudG9Mb3dlckNhc2UoKS50cmltKCkpIHtcbiAgICAgICAgICAgICAgICAgICAgICBpc1RhcmdldENhdGVnb3J5ID0gdHJ1ZTtcbiAgICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgICAgLy8gSGF2ZW4ndCBmb3VuZCB0aGUgdGFyZ2V0IGNhdGVnb3J5IHlldCwga2VlcCBza2lwcGluZ1xuICAgICAgICAgICAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICB9XG4gICAgICAgICAgfVxuICAgICAgICAgIC8vQHRzLWlnbm9yZVxuICAgICAgICAgIGNvbnN0IHNraXBwZWQgPSB2YWx1ZS50ZXN0cz8uZXZlcnkoKHQ6IFRlc3QpID0+IHQub3B0aW9ucz8uc2tpcFJlYXNvblxuICAgICAgICAgICAgICB8fCAob3B0aW9ucz8udGVzdCAhPSBudWxsICYmIG9wdGlvbnMudGVzdC50b0xvd2VyQ2FzZSgpICE9PSB0Lm5hbWUudG9Mb3dlckNhc2UoKSkpO1xuXG4gICAgICAgICAgaWYgKCFza2lwcGVkKSB7XG4gICAgICAgICAgICAgIC8vQHRzLWlnbm9yZVxuICAgICAgICAgICAgICBjb25zdCBza2lwcGVkQ291bnQgPSAodmFsdWUudGVzdHMgPz8gW10pLmZpbHRlcigodDogVGVzdCkgPT5cbiAgICAgICAgICAgICAgICB0Lm9wdGlvbnM/LnNraXBSZWFzb24gfHwgKG9wdGlvbnM/LnRlc3QgIT0gbnVsbCAmJiBvcHRpb25zLnRlc3QudG9Mb3dlckNhc2UoKSAhPT0gdC5uYW1lLnRvTG93ZXJDYXNlKCkpXG4gICAgICAgICAgICAgICkubGVuZ3RoO1xuICAgICAgICAgICAgICBzdGRMb2coYFBhY2thZ2UgdGVzdGluZzogU3RhcnRlZCB7eyR7a2V5fX19JHtza2lwcGVkQ291bnQgPiAwID8gYCBza2lwcGVkIHt7JHtza2lwcGVkQ291bnR9fX1gIDogJyd9YCk7XG4gICAgICAgICAgICAgIHZhbHVlLmJlZm9yZVN0YXR1cyA9IGF3YWl0IGludm9rZUNhdGVnb3J5TWV0aG9kKHZhbHVlLmJlZm9yZSwga2V5KTtcbiAgICAgICAgICB9XG4gICAgICAgICAgbGV0IHQgPSB2YWx1ZS50ZXN0cyA/PyBbXTtcblxuICAgICAgICAgIGlmIChvcHRpb25zLnN0cmVzc1Rlc3QpIHtcbiAgICAgICAgICAgICAgdCA9IHQuZmlsdGVyKChlKSA9PiBlLm9wdGlvbnM/LnN0cmVzc1Rlc3QpO1xuICAgICAgICAgICAgICB0ID0gc2h1ZmZsZSh0KTtcbiAgICAgICAgICB9XG5cbiAgICAgICAgICBpZiAoKG9wdGlvbnMudGFncz8ubGVuZ3RoID8/IDApID4gMCkge1xuICAgICAgICAgICAgICB0ID0gdC5maWx0ZXIoKGUpID0+XG4gICAgICAgICAgICAgICAgICBlLm9wdGlvbnM/LnRhZ3M/LnNvbWUodGFnID0+IChvcHRpb25zPy50YWdzID8/IFtdKS5pbmNsdWRlcyh0YWcpKVxuICAgICAgICAgICAgICApO1xuICAgICAgICAgIH1cblxuICAgICAgICAgIGxldCByZXM6IFRlc3RSZXN1bHRFeHRlbmRlZFtdO1xuICAgICAgICAgIGlmICh2YWx1ZS5iZWZvcmVTdGF0dXMpIHtcbiAgICAgICAgICAgICAgcmVzID0gQXJyYXkuZnJvbSh0Lm1hcCgodGVzdEVsZW0pID0+IHtcbiAgICAgICAgICAgICAgICAgIHJldHVybiB7XG4gICAgICAgICAgICAgICAgICAgICAgZGF0ZTogbmV3IERhdGUoKS50b0lTT1N0cmluZygpLFxuICAgICAgICAgICAgICAgICAgICAgIGNhdGVnb3J5OiBrZXksXG4gICAgICAgICAgICAgICAgICAgICAgbmFtZTogdGVzdEVsZW0ubmFtZSxcbiAgICAgICAgICAgICAgICAgICAgICBzdWNjZXNzOiBmYWxzZSxcbiAgICAgICAgICAgICAgICAgICAgICByZXN1bHQ6ICdiZWZvcmUoKSBmYWlsZWQnLFxuICAgICAgICAgICAgICAgICAgICAgIG1zOiAwLFxuICAgICAgICAgICAgICAgICAgICAgIHNraXBwZWQ6IGZhbHNlLFxuICAgICAgICAgICAgICAgICAgICAgIGxvZ3M6ICcnLFxuICAgICAgICAgICAgICAgICAgICAgIG93bmVyOiBwYWNrYWdlT3duZXIsXG4gICAgICAgICAgICAgICAgICAgICAgcGFja2FnZTogcGFja2FnZV8ubmFtZSxcbiAgICAgICAgICAgICAgICAgICAgICB3aWRnZXRzRGlmZmVyZW5jZTogMCxcbiAgICAgICAgICAgICAgICAgICAgICBmbGFraW5nOiBERy5UZXN0LmlzUmVwcm9kdWNpbmdcbiAgICAgICAgICAgICAgICAgIH07XG4gICAgICAgICAgICAgIH0pKTtcbiAgICAgICAgICAgICAgcmVzLmZvckVhY2goYXN5bmMgKHRlc3QpID0+IGF3YWl0IGdyb2suc2hlbGwucmVwb3J0VGVzdCgncGFja2FnZScsIHRlc3QpKTtcbiAgICAgICAgICB9IGVsc2VcbiAgICAgICAgICAgICAgcmVzID0gYXdhaXQgaW52b2tlVGVzdHNJbkNhdGVnb3J5KHZhbHVlLCBvcHRpb25zLCBza2lwcGluZ0NhdGVnb3JpZXMpO1xuICAgICAgICAgIGNvbnN0IGRhdGE6IFRlc3RSZXN1bHRFeHRlbmRlZFtdID0gcmVzLmZpbHRlcigoZCkgPT4gZC5yZXN1bHQgIT0gJ3NraXBwZWQnKTtcblxuICAgICAgICAgIGlmICghc2tpcHBlZClcbiAgICAgICAgICAgICAgdmFsdWUuYWZ0ZXJTdGF0dXMgPSBhd2FpdCBpbnZva2VDYXRlZ29yeU1ldGhvZCh2YWx1ZS5hZnRlciwga2V5KTtcblxuICAgICAgICAgIC8vIENsZWFyIGFmdGVyIGNhdGVnb3J5XG4gICAgICAgICAgLy8gZ3Jvay5zaGVsbC5jbG9zZUFsbCgpO1xuICAgICAgICAgIC8vIERHLkJhbGxvb24uY2xvc2VBbGwoKTtcbiAgICAgICAgICBpZiAodmFsdWUuYWZ0ZXJTdGF0dXMpIHtcbiAgICAgICAgICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IENhdGVnb3J5IGFmdGVyKCkge3ske2tleX19fSBmYWlsZWRgKTtcbiAgICAgICAgICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFJlc3VsdCBmb3Ige3ske2tleX19fSBhZnRlcjogJHt2YWx1ZS5hZnRlclN0YXR1c31gKTtcbiAgICAgICAgICAgICAgZGF0YS5wdXNoKHtcbiAgICAgICAgICAgICAgICAgIGRhdGU6IG5ldyBEYXRlKCkudG9JU09TdHJpbmcoKSxcbiAgICAgICAgICAgICAgICAgIGNhdGVnb3J5OiBrZXksXG4gICAgICAgICAgICAgICAgICBuYW1lOiAnYWZ0ZXInLFxuICAgICAgICAgICAgICAgICAgc3VjY2VzczogZmFsc2UsXG4gICAgICAgICAgICAgICAgICByZXN1bHQ6IHZhbHVlLmFmdGVyU3RhdHVzLFxuICAgICAgICAgICAgICAgICAgbXM6IDAsXG4gICAgICAgICAgICAgICAgICBza2lwcGVkOiBmYWxzZSxcbiAgICAgICAgICAgICAgICAgIGxvZ3M6ICcnLFxuICAgICAgICAgICAgICAgICAgb3duZXI6IHBhY2thZ2VPd25lcixcbiAgICAgICAgICAgICAgICAgIHBhY2thZ2U6IHBhY2thZ2VfLm5hbWUsXG4gICAgICAgICAgICAgICAgICB3aWRnZXRzRGlmZmVyZW5jZTogMCxcbiAgICAgICAgICAgICAgICAgIGZsYWtpbmc6IERHLlRlc3QuaXNSZXByb2R1Y2luZ1xuICAgICAgICAgICAgICB9KTtcbiAgICAgICAgICB9XG4gICAgICAgICAgaWYgKHZhbHVlLmJlZm9yZVN0YXR1cykge1xuICAgICAgICAgICAgICBzdGRMb2coYFBhY2thZ2UgdGVzdGluZzogQ2F0ZWdvcnkgYmVmb3JlKCkge3ske2tleX19fSBmYWlsZWRgKTtcbiAgICAgICAgICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFJlc3VsdCBmb3Ige3ske2tleX19fSBiZWZvcmU6ICR7dmFsdWUuYmVmb3JlU3RhdHVzfWApO1xuICAgICAgICAgICAgICBkYXRhLnB1c2goe1xuICAgICAgICAgICAgICAgICAgZGF0ZTogbmV3IERhdGUoKS50b0lTT1N0cmluZygpLFxuICAgICAgICAgICAgICAgICAgY2F0ZWdvcnk6IGtleSxcbiAgICAgICAgICAgICAgICAgIG5hbWU6ICdiZWZvcmUnLFxuICAgICAgICAgICAgICAgICAgc3VjY2VzczogZmFsc2UsXG4gICAgICAgICAgICAgICAgICByZXN1bHQ6IHZhbHVlLmJlZm9yZVN0YXR1cyxcbiAgICAgICAgICAgICAgICAgIG1zOiAwLFxuICAgICAgICAgICAgICAgICAgc2tpcHBlZDogZmFsc2UsXG4gICAgICAgICAgICAgICAgICBsb2dzOiAnJyxcbiAgICAgICAgICAgICAgICAgIG93bmVyOiBwYWNrYWdlT3duZXIsXG4gICAgICAgICAgICAgICAgICBwYWNrYWdlOiBwYWNrYWdlXy5uYW1lLFxuICAgICAgICAgICAgICAgICAgd2lkZ2V0c0RpZmZlcmVuY2U6IDAsXG4gICAgICAgICAgICAgICAgICBmbGFraW5nOiBERy5UZXN0LmlzUmVwcm9kdWNpbmdcbiAgICAgICAgICAgICAgfSk7XG4gICAgICAgICAgfVxuICAgICAgICAgIHJlc3VsdHMucHVzaCguLi5kYXRhKTtcblxuICAgICAgICAgIC8vIElmIHJldHVybk9uRmFpbCBpcyBzZXQgYW5kIGEgdGVzdCBmYWlsZWQgKG90aGVyIHRoYW4gc2tpcFRvVGVzdCksIHN0b3AgcHJvY2Vzc2luZyBtb3JlIGNhdGVnb3JpZXNcbiAgICAgICAgICBpZiAob3B0aW9ucy5yZXR1cm5PbkZhaWwgJiYgZGF0YS5zb21lKChkKSA9PiAhZC5zdWNjZXNzICYmICFkLnNraXBwZWQgJiYgZC5uYW1lICE9PSBvcHRpb25zLnNraXBUb1Rlc3QpKVxuICAgICAgICAgICAgICBicmVhaztcbiAgICAgIH1cbiAgICB9IGZpbmFsbHkge1xuICAgICAgcmVzZXRDb25zb2xlKCk7XG4gICAgfVxuICAgIGlmIChvcHRpb25zLnRlc3RDb250ZXh0IS5jYXRjaFVuaGFuZGxlZCAmJiAoIURHLlRlc3QuaXNJbkJlbmNobWFyaykpIHtcbiAgICAgIGF3YWl0IGRlbGF5KDEwMDApO1xuICAgICAgY29uc3QgZXJyb3IgPSBhd2FpdCBncm9rLnNoZWxsLmxhc3RFcnJvcjtcbiAgICAgIGlmIChlcnJvciAhPSB1bmRlZmluZWQpIHtcbiAgICAgICAgICBjb25zdCBwYXJhbXM6IGFueSA9IHtcbiAgICAgICAgICAgICAgbG9nczogJycsXG4gICAgICAgICAgICAgIGRhdGU6IG5ldyBEYXRlKCkudG9JU09TdHJpbmcoKSxcbiAgICAgICAgICAgICAgY2F0ZWdvcnk6ICdVbmhhbmRsZWQgZXhjZXB0aW9ucycsXG4gICAgICAgICAgICAgIG5hbWU6ICdFeGNlcHRpb24nLFxuICAgICAgICAgICAgICByZXN1bHQ6IGVycm9yID8/ICcnLFxuICAgICAgICAgICAgICBzdWNjZXNzOiAhZXJyb3IsXG4gICAgICAgICAgICAgIG1zOiAwLFxuICAgICAgICAgICAgICBza2lwcGVkOiBmYWxzZSxcbiAgICAgICAgICAgICAgb3duZXI6IHBhY2thZ2VPd25lciA/PyAnJyxcbiAgICAgICAgICAgICAgJ3BhY2thZ2UnOiBwYWNrYWdlXy5uYW1lLFxuICAgICAgICAgICAgICB3aWRnZXRzRGlmZmVyZW5jZTogMFxuICAgICAgICAgIH07XG4gICAgICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFVuaGFuZGxlZCBFeGNlcHRpb246ICR7ZXJyb3J9YCk7XG5cbiAgICAgICAgICByZXN1bHRzLnB1c2goey4uLnBhcmFtcywgJ2ZsYWtpbmcnOiBERy5UZXN0LmlzUmVwcm9kdWNpbmcgJiYgIWVycm9yfSk7XG4gICAgICAgICAgKDxhbnk+cGFyYW1zKS5wYWNrYWdlID0gcGFja2FnZV8ubmFtZTtcbiAgICAgICAgICBhd2FpdCBncm9rLnNoZWxsLnJlcG9ydFRlc3QoJ3BhY2thZ2UnLCBwYXJhbXMpO1xuICAgICAgfVxuICAgIH1cbiAgfVxufVxuXG5hc3luYyBmdW5jdGlvbiBnZXRSZXN1bHQoeDogYW55KTogUHJvbWlzZTxzdHJpbmc+IHtcbiAgcmV0dXJuIGAke3gudG9TdHJpbmcoKX1cXG4ke3guc3RhY2sgPyAoYXdhaXQgREcuTG9nZ2VyLnRyYW5zbGF0ZVN0YWNrVHJhY2UoeC5zdGFjaykpIDogJyd9YDtcbn1cblxuYXN5bmMgZnVuY3Rpb24gZXhlY1Rlc3QodDogVGVzdCwgcHJlZGljYXRlOiBzdHJpbmcgfCB1bmRlZmluZWQsIGxvZ3M6IGFueVtdLFxuICB0ZXN0VGltZW91dD86IG51bWJlciwgcGFja2FnZU5hbWU/OiBzdHJpbmcsIHZlcmJvc2U/OiBib29sZWFuXG4pOiBQcm9taXNlPFRlc3RSZXN1bHQgfCB1bmRlZmluZWQ+IHtcbiAgbG9ncy5sZW5ndGggPSAwO1xuICBsZXQgcjogVGVzdFJlc3VsdDtcbiAgbGV0IHR5cGU6IHN0cmluZyA9ICdwYWNrYWdlJztcbiAgY29uc3QgZmlsdGVyID0gcHJlZGljYXRlICE9IHVuZGVmaW5lZCAmJiAodC5uYW1lLnRvTG93ZXJDYXNlKCkgIT09IHByZWRpY2F0ZS50b0xvd2VyQ2FzZSgpKTtcbiAgbGV0IHNraXAgPSB0Lm9wdGlvbnM/LnNraXBSZWFzb24gfHwgZmlsdGVyO1xuICBsZXQgc2tpcFJlYXNvbiA9IGZpbHRlciA/ICdza2lwcGVkJyA6IHQub3B0aW9ucz8uc2tpcFJlYXNvbjtcblxuICBpZiAoREcuVGVzdC5pc0luQmVuY2htYXJrICYmICF0Lm9wdGlvbnM/LmJlbmNobWFyaykge1xuICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBTa2lwcGVkIHt7JHt0LmNhdGVnb3J5fX19IHt7JHt0Lm5hbWV9fX0gZG9lc250IGF2YWlsYWJsZSBpbiBiZW5jaG1hcmsgbW9kZWApO1xuICAgIHJldHVybiB1bmRlZmluZWQ7XG4gIH1cblxuICBpZiAoc2tpcCAmJiAhREcuVGVzdC5pc0luQmVuY2htYXJrKVxuICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBTa2lwcGVkIHt7JHt0LmNhdGVnb3J5fX19IHt7JHt0Lm5hbWV9fX1gKTtcbiAgaWYgKCFza2lwKVxuICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBTdGFydGVkIHt7JHt0LmNhdGVnb3J5fX19IHt7JHt0Lm5hbWV9fX1gKTtcbiAgY29uc3Qgc3RhcnQgPSBEYXRlLm5vdygpO1xuICBjb25zdCBzdGFydERhdGUgPSBuZXcgRGF0ZShzdGFydCkudG9JU09TdHJpbmcoKTtcbiAgdHJ5IHtcbiAgICBpZiAoc2tpcClcbiAgICAgIHIgPSB7IG5hbWU6IHQubmFtZSwgb3duZXI6dC5vcHRpb25zPy5vd25lciA/PyAnJywgY2F0ZWdvcnk6IHQuY2F0ZWdvcnksIGxvZ3M6ICcnLCBkYXRlOiBzdGFydERhdGUsIHN1Y2Nlc3M6IHRydWUsIHJlc3VsdDogc2tpcFJlYXNvbiEsIG1zOiAwLCBza2lwcGVkOiB0cnVlLCBwYWNrYWdlOiBwYWNrYWdlTmFtZSA/PyAnJywgZmxha2luZzogREcuVGVzdC5pc1JlcHJvZHVjaW5nfTtcbiAgICBlbHNlIHtcbiAgICAgIGxldCB0aW1lb3V0XyA9IHRlc3RUaW1lb3V0ID8/IFNUQU5EQVJUX1RJTUVPVVQ7XG5cbiAgICAgIGlmIChERy5UZXN0LmlzUHJvZmlsaW5nKVxuICAgICAgICBjb25zb2xlLnByb2ZpbGUoYCR7dC5jYXRlZ29yeX06ICR7dC5uYW1lfWApO1xuXG4gICAgICByID0geyBuYW1lOiB0Lm5hbWUsIG93bmVyOnQub3B0aW9ucz8ub3duZXIgPz8gJycsIGNhdGVnb3J5OiB0LmNhdGVnb3J5LCBsb2dzOiAnJywgZGF0ZTogc3RhcnREYXRlLCBzdWNjZXNzOiB0cnVlLCByZXN1bHQ6IChhd2FpdCB0aW1lb3V0KHQudGVzdCwgdGltZW91dF8pKS50b1N0cmluZygpID8/ICdPSycsIG1zOiAwLCBza2lwcGVkOiBmYWxzZSAsIHBhY2thZ2U6IHBhY2thZ2VOYW1lID8/ICcnLCBmbGFraW5nOiBERy5UZXN0LmlzUmVwcm9kdWNpbmd9O1xuXG4gICAgICBpZiAoREcuVGVzdC5pc1Byb2ZpbGluZykge1xuICAgICAgICBjb25zb2xlLnByb2ZpbGVFbmQoYCR7dC5jYXRlZ29yeX06ICR7dC5uYW1lfWApO1xuICAgICAgICBncm9rLnNoZWxsLmluZm8oYFByb2ZpbGluZyBvZiAke3QuY2F0ZWdvcnl9OiAke3QubmFtZX0gZmluaXNoZWQgXFxuIFBsZWFzZSBlbnN1cmUgdGhhdCB5b3UgaGF2ZSBvcGVuZWQgRGV2VG9vbHMgKEYxMikgLyBQZXJmb3JtYW5jZSBwYW5lbCBiZWZvcmUgdGVzdCBzdGFydHMuYCk7XG4gICAgICB9XG4gICAgfVxuICB9IGNhdGNoICh4OiBhbnkpIHtcbiAgICBzdGRFcnJvcih4KTtcbiAgICByID0geyBuYW1lOiB0Lm5hbWUsIG93bmVyOnQub3B0aW9ucz8ub3duZXIgPz8gJycsIGNhdGVnb3J5OiB0LmNhdGVnb3J5LCBsb2dzOiAnJywgZGF0ZTogc3RhcnREYXRlLCBzdWNjZXNzOiBmYWxzZSwgcmVzdWx0OiBhd2FpdCBnZXRSZXN1bHQoeCksIG1zOiAwLCBza2lwcGVkOiBmYWxzZSwgcGFja2FnZTogcGFja2FnZU5hbWUgPz8gJycsIGZsYWtpbmc6IGZhbHNlfTtcbiAgfVxuICBpZiAodC5vcHRpb25zPy5pc0FnZ3JlZ2F0ZWQgJiYgci5yZXN1bHQuY29uc3RydWN0b3IgPT09IERHLkRhdGFGcmFtZSkge1xuICAgIGNvbnN0IGNvbCA9IHIucmVzdWx0LmNvbCgnc3VjY2VzcycpO1xuICAgIGlmIChjb2wpXG4gICAgICByLnN1Y2Nlc3MgPSBjb2wuc3RhdHMuc3VtID09PSBjb2wubGVuZ3RoO1xuICAgIGlmICghdmVyYm9zZSkge1xuICAgICAgY29uc3QgZGYgPSByLnJlc3VsdDtcbiAgICAgIGRmLmNvbHVtbnMucmVtb3ZlKCdzdGFjaycpO1xuICAgICAgZGYucm93cy5yZW1vdmVXaGVyZSgocikgPT4gci5nZXQoJ3N1Y2Nlc3MnKSk7XG4gICAgICByLnJlc3VsdCA9IGRmO1xuICAgIH1cbiAgICByLnJlc3VsdCA9IHIucmVzdWx0LnRvQ3N2KCk7XG4gIH1cbiAgci5sb2dzID0gbG9ncy5qb2luKCdcXG4nKTtcbiAgci5tcyA9IERhdGUubm93KCkgLSBzdGFydDtcbiAgaWYgKCFza2lwKVxuICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBGaW5pc2hlZCB7eyR7dC5jYXRlZ29yeX19fSB7eyR7dC5uYW1lfX19IHdpdGgge3ske3Iuc3VjY2VzcyA/ICdzdWNjZXNzJyA6ICdlcnJvcid9fX0gZm9yICR7ci5tc30gbXNgKTtcbiAgaWYgKCFyLnN1Y2Nlc3MpIHtcbiAgICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBSZXN1bHQgZm9yIHt7JHt0LmNhdGVnb3J5fX19IHt7JHt0Lm5hbWV9fX06ICR7ci5yZXN1bHR9YCk7XG4gIH1cbiAgci5jYXRlZ29yeSA9IHQuY2F0ZWdvcnk7XG4gIHIubmFtZSA9IHQubmFtZTtcbiAgci5vd25lciA9IHQub3B0aW9ucz8ub3duZXIgPz8gJyc7XG4gIGlmICghZmlsdGVyKSB7XG4gICAgbGV0IHBhcmFtcyA9IHtcbiAgICAgICdzdWNjZXNzJzogci5zdWNjZXNzLCAncmVzdWx0Jzogci5yZXN1bHQsICdtcyc6IHIubXMsICdkYXRlJzogci5kYXRlLFxuICAgICAgJ3NraXBwZWQnOiByLnNraXBwZWQsICdjYXRlZ29yeSc6IHQuY2F0ZWdvcnksICduYW1lJzogdC5uYW1lLCAnbG9ncyc6IHIubG9ncywgJ293bmVyJzogci5vd25lcixcbiAgICAgICdmbGFraW5nJzogREcuVGVzdC5pc1JlcHJvZHVjaW5nICYmIHIuc3VjY2VzcyxcbiAgICAgICdwYWNrYWdlJzogci5wYWNrYWdlXG4gICAgfTtcbiAgICBpZiAoci5yZXN1bHQuY29uc3RydWN0b3IgPT0gT2JqZWN0KSB7XG4gICAgICBjb25zdCByZXMgPSBPYmplY3Qua2V5cyhyLnJlc3VsdCkucmVkdWNlKChhY2MsIGspID0+ICh7IC4uLmFjYywgWydyZXN1bHQuJyArIGtdOiByLnJlc3VsdFtrXSB9KSwge30pO1xuICAgICAgcGFyYW1zID0geyAuLi5wYXJhbXMsIC4uLnJlcyB9O1xuICAgIH1cblxuICAgIGlmIChwYXJhbXMucmVzdWx0IGluc3RhbmNlb2YgREcuRGF0YUZyYW1lKVxuICAgICAgcGFyYW1zLnJlc3VsdCA9IEpTT04uc3RyaW5naWZ5KHBhcmFtcy5yZXN1bHQ/LnRvSnNvbigpKSB8fCAnJztcbiAgICBhd2FpdCBncm9rLnNoZWxsLnJlcG9ydFRlc3QodHlwZSwgcGFyYW1zKTtcbiAgfVxuICByZXR1cm4gcjtcbn1cblxuZXhwb3J0IGZ1bmN0aW9uIHNodWZmbGUoYXJyYXk6IGFueVtdKTogYW55W10ge1xuICBjb25zdCBuZXdBcnIgPSBhcnJheS5zbGljZSgpO1xuICBuZXdBcnIuc29ydCgoKSA9PiBNYXRoLnJhbmRvbSgpIC0gMC41KTtcbiAgcmV0dXJuIG5ld0Fycjtcbn1cblxuLyogV2FpdHMgW21zXSBtaWxsaXNlY29uZHMgKi9cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBkZWxheShtczogbnVtYmVyKSB7XG4gIGF3YWl0IG5ldyBQcm9taXNlKChyKSA9PiBzZXRUaW1lb3V0KHIsIG1zKSk7XG59XG5cbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBhd2FpdENoZWNrKGNoZWNrSGFuZGxlcjogKCkgPT4gYm9vbGVhbixcbiAgZXJyb3I6IHN0cmluZyA9ICdUaW1lb3V0IGV4Y2VlZGVkJywgd2FpdDogbnVtYmVyID0gNTAwLCBpbnRlcnZhbDogbnVtYmVyID0gNTApOiBQcm9taXNlPGFueT4ge1xuICByZXR1cm4gbmV3IFByb21pc2UoKHJlc29sdmUsIHJlamVjdCkgPT4ge1xuICAgIHNldFRpbWVvdXQoKCkgPT4ge1xuICAgICAgY2xlYXJJbnRlcnZhbChpbnRlcnZhbElkKTtcbiAgICAgIHJlamVjdChuZXcgRXJyb3IoZXJyb3IpKTtcbiAgICB9LCB3YWl0KTtcbiAgICAvLyBAdHMtaWdub3JlXG4gICAgY29uc3QgaW50ZXJ2YWxJZDogVGltZW91dCA9IHNldEludGVydmFsKCgpID0+IHtcbiAgICAgIGlmIChjaGVja0hhbmRsZXIoKSkge1xuICAgICAgICBjbGVhckludGVydmFsKGludGVydmFsSWQpO1xuICAgICAgICByZXNvbHZlKG51bGwpO1xuICAgICAgfVxuICAgIH0sIGludGVydmFsKTtcbiAgfSk7XG59XG5cbi8vIFJldHVybnMgdGVzdCBleGVjdXRpb24gcmVzdWx0IG9yIGFuIGVycm9yIGluIGNhc2Ugb2YgdGltZW91dFxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHRpbWVvdXQoZnVuYzogKCkgPT4gUHJvbWlzZTxhbnk+LCB0ZXN0VGltZW91dDogbnVtYmVyLCB0aW1lb3V0UmVhc29uOiBzdHJpbmcgPSAnRVhFQ1VUSU9OIFRJTUVPVVQnKTogUHJvbWlzZTxhbnk+IHtcbiAgbGV0IHRpbWVvdXQ6IGFueSA9IG51bGw7XG4gIGNvbnN0IHRpbWVvdXRQcm9taXNlID0gbmV3IFByb21pc2U8YW55PigoXywgcmVqZWN0KSA9PiB7XG4gICAgdGltZW91dCA9IHNldFRpbWVvdXQoKCkgPT4ge1xuICAgICAgLy8gZXNsaW50LWRpc2FibGUtbmV4dC1saW5lIHByZWZlci1wcm9taXNlLXJlamVjdC1lcnJvcnNcbiAgICAgIHJlamVjdCh0aW1lb3V0UmVhc29uKTtcbiAgICB9LCB0ZXN0VGltZW91dCk7XG4gIH0pO1xuICB0cnkge1xuICAgIHJldHVybiBhd2FpdCBQcm9taXNlLnJhY2UoW2Z1bmMoKSwgdGltZW91dFByb21pc2VdKTtcbiAgfSBmaW5hbGx5IHtcbiAgICBpZiAodGltZW91dClcbiAgICAgIGNsZWFyVGltZW91dCh0aW1lb3V0KTtcbiAgfVxufVxuXG5leHBvcnQgZnVuY3Rpb24gaXNEaWFsb2dQcmVzZW50KGRpYWxvZ1RpdGxlOiBzdHJpbmcpOiBib29sZWFuIHtcbiAgY29uc3QgZGlhbG9ncyA9IERHLkRpYWxvZy5nZXRPcGVuRGlhbG9ncygpO1xuICBmb3IgKGxldCBpID0gMDsgaSA8IGRpYWxvZ3MubGVuZ3RoOyBpKyspIHtcbiAgICBpZiAoZGlhbG9nc1tpXS50aXRsZSA9PSBkaWFsb2dUaXRsZSlcbiAgICAgIHJldHVybiB0cnVlO1xuICB9XG4gIHJldHVybiBmYWxzZTtcbn1cblxuLyoqIEV4cGVjdHMgYW4gYXN5bmNocm9ub3VzIHtAbGluayBhY3Rpb259IHRvIHRocm93IGFuIGV4Y2VwdGlvbi4gVXNlIHtAbGluayBjaGVja30gdG8gcGVyZm9ybVxuICogZGVlcGVyIGluc3BlY3Rpb24gb2YgdGhlIGV4Y2VwdGlvbiBpZiBuZWNlc3NhcnkuXG4gKiBAcGFyYW0gIHtmdW5jdGlvbigpOiBQcm9taXNlPHZvaWQ+fSBhY3Rpb25cbiAqIEBwYXJhbSAge2Z1bmN0aW9uKGFueSk6IGJvb2xlYW59IGNoZWNrXG4gKiBAcmV0dXJuIHtQcm9taXNlPHZvaWQ+fVxuICovXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gZXhwZWN0RXhjZXB0aW9uQXN5bmMoYWN0aW9uOiAoKSA9PiBQcm9taXNlPHZvaWQ+LFxuICBjaGVjaz86IChleGNlcHRpb246IGFueSkgPT4gYm9vbGVhbik6IFByb21pc2U8dm9pZD4ge1xuICBsZXQgY2F1Z2h0OiBib29sZWFuID0gZmFsc2U7XG4gIGxldCBjaGVja2VkOiBib29sZWFuID0gZmFsc2U7XG4gIHRyeSB7XG4gICAgYXdhaXQgYWN0aW9uKCk7XG4gIH0gY2F0Y2ggKGUpIHtcbiAgICBjYXVnaHQgPSB0cnVlO1xuICAgIGNoZWNrZWQgPSAhY2hlY2sgfHwgY2hlY2soZSk7XG4gIH0gZmluYWxseSB7XG4gICAgaWYgKCFjYXVnaHQpXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoJ0FuIGV4Y2VwdGlvbiBpcyBleHBlY3RlZCBidXQgbm90IHRocm93bicpO1xuICAgIGlmICghY2hlY2tlZClcbiAgICAgIHRocm93IG5ldyBFcnJvcignQW4gZXhwZWN0ZWQgZXhjZXB0aW9uIGlzIHRocm93biwgYnV0IGl0IGRvZXMgbm90IHNhdGlzZnkgdGhlIGNvbmRpdGlvbicpO1xuICB9XG59XG5cbmNvbnN0IGNhdERGID0gREcuRGF0YUZyYW1lLmZyb21Db2x1bW5zKFtERy5Db2x1bW4uZnJvbVN0cmluZ3MoJ2NvbCcsIFsndmFsMScsICd2YWwyJywgJ3ZhbDMnXSldKTtcblxuLyoqXG4gKiBVbml2ZXJzYWwgdGVzdCBmb3Igdmlld2Vycy4gSXQgc2VhcmNoIHZpZXdlcnMgaW4gRE9NIGJ5IHRhZ3M6IGNhbnZhcywgc3ZnLCBpbWcsIGlucHV0LCBoMSwgYVxuICogQHBhcmFtICB7c3RyaW5nfSB2IFZpZXdlciBuYW1lXG4gKiBAcGFyYW0gIHtfREcuRGF0YUZyYW1lfSBkZiBEYXRhZnJhbWUgdG8gdXNlLiBTaG91bGQgaGF2ZSBhdCBsZWFzdCAzIHJvd3NcbiAqIEBwYXJhbSAge2Jvb2xlYW59IG9wdGlvbnMuZGV0ZWN0U2VtYW50aWNUeXBlcyBTcGVjaWZ5IHdoZXRoZXIgdG8gZGV0ZWN0IHNlbWFudGljIHR5cGVzIG9yIG5vdFxuICogQHBhcmFtICB7Ym9vbGVhbn0gb3B0aW9ucy5yZWFkT25seSBJZiBzZXQgdG8gdHJ1ZSwgdGhlIGRhdGFmcmFtZSB3aWxsIG5vdCBiZSBtb2RpZmllZCBkdXJpbmcgdGhlIHRlc3RcbiAqIEBwYXJhbSAge2Jvb2xlYW59IG9wdGlvbnMuYXJiaXRyYXJ5RGZUZXN0IElmIHNldCB0byBmYWxzZSwgdGVzdCBvbiBhcmJpdHJhcnkgZGF0YWZyYW1lXG4gKiAob25lIGNhdGVnb3JpY2FsIGNvbHVtbikgd2lsbCBub3QgYmUgcGVyZm9ybWVkXG4gKiBAcGFyYW0gIHtvYmplY3R9IG9wdGlvbnMgTGlzdCBvZiBvcHRpb25zIChvcHRpb25hbClcbiAqIEByZXR1cm4ge1Byb21pc2U8dm9pZD59IFRoZSB0ZXN0IGlzIGNvbnNpZGVyZWQgc3VjY2Vzc2Z1bCBpZiBpdCBjb21wbGV0ZXMgd2l0aG91dCBlcnJvcnNcbiAqL1xuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHRlc3RWaWV3ZXIodjogc3RyaW5nLCBkZjogX0RHLkRhdGFGcmFtZSwgb3B0aW9ucz86IHtcbiAgZGV0ZWN0U2VtYW50aWNUeXBlcz86IGJvb2xlYW4sIHJlYWRPbmx5PzogYm9vbGVhbiwgYXJiaXRyYXJ5RGZUZXN0PzogYm9vbGVhbixcbiAgcGFja2FnZU5hbWU/OiBzdHJpbmcsIGF3YWl0Vmlld2VyPzogKHZpZXdlcjogX0RHLlZpZXdlcikgPT4gUHJvbWlzZTx2b2lkPlxufSk6IFByb21pc2U8dm9pZD4ge1xuICBjb25zdCBwYWNrYWdlTmFtZSA9IG9wdGlvbnM/LnBhY2thZ2VOYW1lID8/ICcnO1xuICBpZiAob3B0aW9ucz8uZGV0ZWN0U2VtYW50aWNUeXBlcylcbiAgICBhd2FpdCBncm9rLmRhdGEuZGV0ZWN0U2VtYW50aWNUeXBlcyhkZik7XG4gIGNvbnN0IHR2ID0gZ3Jvay5zaGVsbC5hZGRUYWJsZVZpZXcoZGYpO1xuXG4gIHRyeSB7XG4gICAgLy8xLiBPcGVuLCBkbyBub3RoaW5nIGFuZCBjbG9zZVxuICAgIGF3YWl0IHRlc3RWaWV3ZXJJbnRlcm5hbCh0diwgdiwgcGFja2FnZU5hbWUsIGdyb2suZXZlbnRzLm9uVmlld2VyQWRkZWQpO1xuICAgIC8vaW4gY2FzZSB2aWV3ZXIgd2l0aCBhc3luYyByZW5kZXJpbmcgLSB3YWl0IGZvciByZW5kZXIgdG8gY29tcGxldGVcbiAgICBpZiAob3B0aW9ucz8uYXdhaXRWaWV3ZXIpXG4gICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCB1bmRlZmluZWQsIG9wdGlvbnMhLmF3YWl0Vmlld2VyKTtcblxuICAgIC8vMi4gT3BlbiB2aWV3ZXIsIHJ1biBzZWxlY3Rpb24sIGZpbHRlciwgZXRjLiBhbmQgY2xvc2VcbiAgICBpZiAoIW9wdGlvbnM/LnJlYWRPbmx5KSB7XG4gICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCBzZWxlY3RGaWx0ZXJDaGFuZ2VDdXJyZW50KTtcbiAgICAgIGlmIChvcHRpb25zPy5hd2FpdFZpZXdlcilcbiAgICAgICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCwgc2VsZWN0RmlsdGVyQ2hhbmdlQ3VycmVudCwgb3B0aW9ucyEuYXdhaXRWaWV3ZXIpO1xuICAgIH1cblxuICAgIC8vMi4gT3BlbiB2aWV3ZXIsIGNoYW5nZSBvcHRpb25zLCBzYXZlIGxheW91dCBhbmQgY2xvc2VcbiAgICBsZXQgcHJvcHNBbmRMYXlvdXQ6IHsgbGF5b3V0OiBhbnksIHNhdmVkUHJvcHM6IGFueSB9IHwgbnVsbCA9IG51bGw7XG4gICAgcHJvcHNBbmRMYXlvdXQgPSBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCBjaGFuZ2VPcHRpb25zU2F2ZUxheW91dCk7XG4gICAgaWYgKG9wdGlvbnM/LmF3YWl0Vmlld2VyKVxuICAgICAgcHJvcHNBbmRMYXlvdXQgPSBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLFxuICAgICAgICBjaGFuZ2VPcHRpb25zU2F2ZUxheW91dCwgb3B0aW9ucyEuYXdhaXRWaWV3ZXIpXG5cbiAgICAvLzMuIExvYWQgbGF5b3V0XG4gICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3TGF5b3V0QXBwbGllZCwgbG9hZExheW91dCwgdW5kZWZpbmVkLCBwcm9wc0FuZExheW91dD8ubGF5b3V0LFxuICAgICAgeyBzYXZlZFByb3BzOiBwcm9wc0FuZExheW91dD8uc2F2ZWRQcm9wcyB9KTtcbiAgICBpZiAob3B0aW9ucz8uYXdhaXRWaWV3ZXIpXG4gICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdMYXlvdXRBcHBsaWVkLCBsb2FkTGF5b3V0LCBvcHRpb25zIS5hd2FpdFZpZXdlcixcbiAgICAgICAgcHJvcHNBbmRMYXlvdXQ/LmxheW91dCwgeyBzYXZlZFByb3BzOiBwcm9wc0FuZExheW91dD8uc2F2ZWRQcm9wcyB9KTtcblxuICAgIC8vNC4gT3BlbiB2aWV3ZXIgb24gYXJiaXRhcnkgZGF0YXNldFxuICAgIGlmIChvcHRpb25zPy5hcmJpdHJhcnlEZlRlc3QgIT09IGZhbHNlKSB7XG4gICAgICB0di5kYXRhRnJhbWUgPSBjYXRERjtcbiAgICAgIGF3YWl0IGRlbGF5KDUwKTtcbiAgICAgIGF3YWl0IHRlc3RWaWV3ZXJJbnRlcm5hbCh0diwgdiwgcGFja2FnZU5hbWUsIGdyb2suZXZlbnRzLm9uVmlld2VyQWRkZWQpO1xuICAgICAgaWYgKG9wdGlvbnM/LmF3YWl0Vmlld2VyKVxuICAgICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCB1bmRlZmluZWQsIG9wdGlvbnMhLmF3YWl0Vmlld2VyKTtcbiAgICB9XG5cbiAgICAvLzUuIENhbGwgcG9zdHBvbmVkIGZpbHRlcmluZ1xuICAgIGF3YWl0IHRlc3RWaWV3ZXJJbnRlcm5hbCh0diwgdiwgcGFja2FnZU5hbWUsIGdyb2suZXZlbnRzLm9uVmlld2VyQWRkZWQsIGZpbHRlckFzeW5jKTtcbiAgICBpZiAob3B0aW9ucz8uYXdhaXRWaWV3ZXIpXG4gICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCBmaWx0ZXJBc3luYywgb3B0aW9ucyEuYXdhaXRWaWV3ZXIpO1xuXG4gIH0gZmluYWxseSB7XG4gICAgLy8gY2xvc2VBbGwoKSBpcyBoYW5kbGluZyBieSBjb21tb24gdGVzdCB3b3JrZmxvd1xuICAgIC8vIGdyb2suc2hlbGwuY2xvc2VBbGwoKTtcbiAgICAvLyBERy5CYWxsb29uLmNsb3NlQWxsKCk7XG4gIH1cbn1cbiJdfQ==