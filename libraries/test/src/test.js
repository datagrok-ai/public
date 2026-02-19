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
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoidGVzdC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzIjpbInRlc3QudHMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6Ijs7Ozs7Ozs7O0FBS0EsT0FBTyxFQUFFLFFBQVEsRUFBRSxNQUFNLG1CQUFtQixDQUFDO0FBRTdDLE9BQU8sRUFBRSx1QkFBdUIsRUFBRSxXQUFXLEVBQUUsVUFBVSxFQUFFLHlCQUF5QixFQUFFLGtCQUFrQixFQUFFLE1BQU0scUJBQXFCLENBQUM7QUFFdEksTUFBTSxnQkFBZ0IsR0FBRyxLQUFLLENBQUM7QUFDL0IsTUFBTSxpQkFBaUIsR0FBRyxRQUFRLENBQUM7QUFFbkMsTUFBTSxNQUFNLEdBQUcsT0FBTyxDQUFDLEdBQUcsQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDekMsTUFBTSxPQUFPLEdBQUcsT0FBTyxDQUFDLElBQUksQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDM0MsTUFBTSxPQUFPLEdBQUcsT0FBTyxDQUFDLElBQUksQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDM0MsTUFBTSxRQUFRLEdBQUcsT0FBTyxDQUFDLEtBQUssQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFFN0MsTUFBTSxDQUFDLE1BQU0sS0FBSyxHQUVkLEVBQUUsQ0FBQztBQUVQLE1BQU0sZ0JBQWdCLEdBQUcsWUFBWSxDQUFDO0FBQ3RDLE1BQU0sV0FBVyxHQUFHLE1BQU0sQ0FBQztBQUMzQixNQUFNLGdCQUFnQixHQUFHLFdBQVcsQ0FBQztBQUNyQyxNQUFNLFdBQVcsR0FBRyxNQUFNLENBQUM7QUFDM0IsTUFBTSxhQUFhLEdBQStCLEVBQUUsQ0FBQztBQUNyRCxNQUFNLENBQUMsSUFBSSxlQUF1QixDQUFDO0FBRW5DLE1BQU0sS0FBVyxNQUFNLENBS3RCO0FBTEQsV0FBaUIsTUFBTTtJQUNyQixTQUFnQixPQUFPLENBQUMsS0FBVSxFQUFFLElBQWE7UUFDL0MsSUFBSSxLQUFLLElBQUksSUFBSTtZQUNmLE1BQU0sSUFBSSxLQUFLLENBQUMsR0FBRyxJQUFJLElBQUksSUFBSSxDQUFDLENBQUMsQ0FBQyxPQUFPLENBQUMsQ0FBQyxDQUFDLElBQUksY0FBYyxDQUFDLENBQUM7SUFDcEUsQ0FBQztJQUhlLGNBQU8sVUFHdEIsQ0FBQTtBQUNILENBQUMsRUFMZ0IsTUFBTSxLQUFOLE1BQU0sUUFLdEI7QUEwQ0QsTUFBTSxPQUFPLFdBQVc7SUFNdEIsWUFBWSxjQUF3QixFQUFFLE1BQWdCLEVBQUUsWUFBc0I7UUFKOUUsbUJBQWMsR0FBRyxJQUFJLENBQUM7UUFDdEIsV0FBTSxHQUFHLEtBQUssQ0FBQztRQUNmLGlCQUFZLEdBQUcsS0FBSyxDQUFDO1FBR25CLElBQUksY0FBYyxLQUFLLFNBQVM7WUFBRSxJQUFJLENBQUMsY0FBYyxHQUFHLGNBQWMsQ0FBQztRQUN2RSxJQUFJLE1BQU0sS0FBSyxTQUFTO1lBQUUsSUFBSSxDQUFDLE1BQU0sR0FBRyxNQUFNLENBQUM7UUFDL0MsSUFBSSxZQUFZLEtBQUssU0FBUztZQUFFLElBQUksQ0FBQyxZQUFZLEdBQUcsWUFBWSxDQUFDO0lBQ25FLENBQUM7SUFBQSxDQUFDO0NBQ0g7QUFFRCxNQUFNLE9BQU8sSUFBSTtJQU1mLFlBQVksUUFBZ0IsRUFBRSxJQUFZLEVBQUUsSUFBd0IsRUFBRSxPQUFxQjs7UUFDekYsSUFBSSxDQUFDLFFBQVEsR0FBRyxRQUFRLENBQUM7UUFDekIsSUFBSSxDQUFDLElBQUksR0FBRyxJQUFJLENBQUM7UUFDakIsT0FBTyxhQUFQLE9BQU8sY0FBUCxPQUFPLElBQVAsT0FBTyxHQUFLLEVBQUUsRUFBQztRQUNmLE1BQUEsT0FBTyxDQUFDLE9BQU8sb0NBQWYsT0FBTyxDQUFDLE9BQU8sR0FBSyxnQkFBZ0IsRUFBQztRQUNyQyxJQUFJLENBQUMsT0FBTyxHQUFHLE9BQU8sQ0FBQztRQUN2QixJQUFJLENBQUMsSUFBSSxHQUFHLEdBQXVCLEVBQUU7WUFDbkMsT0FBTyxJQUFJLE9BQU8sQ0FBQyxDQUFPLE9BQU8sRUFBRSxNQUFNLEVBQUUsRUFBRTs7Z0JBQzNDLElBQUksTUFBTSxHQUFHLEVBQUUsQ0FBQztnQkFDaEIsSUFBSTtvQkFDRixJQUFJLEVBQUUsQ0FBQyxJQUFJLENBQUMsU0FBUzt3QkFDbkIsUUFBUSxDQUFDO29CQUVYLElBQUksR0FBRyxHQUFHLE1BQU0sSUFBSSxFQUFFLENBQUM7b0JBQ3ZCLElBQUk7d0JBQ0YsTUFBTSxHQUFHLE1BQUEsR0FBRyxhQUFILEdBQUcsdUJBQUgsR0FBRyxDQUFFLFFBQVEsRUFBRSxtQ0FBSSxFQUFFLENBQUM7cUJBQ2hDO29CQUNELE9BQU8sQ0FBQyxFQUFFO3dCQUNSLE1BQU0sR0FBRyx5Q0FBeUMsQ0FBQzt3QkFDbkQsT0FBTyxDQUFDLEtBQUssQ0FBQyxrREFBa0QsSUFBSSxDQUFDLFFBQVEsSUFBSSxJQUFJLENBQUMsSUFBSSxPQUFPLENBQUMsQ0FBQztxQkFDcEc7aUJBQ0Y7Z0JBQUMsT0FBTyxDQUFNLEVBQUU7b0JBQ2YsTUFBTSxDQUFDLENBQUMsQ0FBQyxDQUFDO2lCQUNYO2dCQUNELE9BQU8sQ0FBQyxNQUFNLENBQUMsQ0FBQztZQUNsQixDQUFDLENBQUEsQ0FBQyxDQUFDO1FBQ0wsQ0FBQyxDQUFBLENBQUM7SUFDSixDQUFDO0NBQ0Y7QUFFRCxNQUFNLE9BQU8sUUFBUTtDQWFwQjtBQUVELE1BQU0sT0FBTyx3QkFBd0I7Q0FFcEM7QUFFRCxNQUFNLE9BQU8sb0JBQW9CO0NBWWhDO0FBRUQsTUFBTSxVQUFnQixTQUFTLENBQUksS0FBb0IsRUFDckQsT0FBMEIsRUFBRSxPQUFtQixFQUFFLEtBQWEsQ0FBQyxFQUFFLFNBQWlCLFNBQVM7O1FBRTNGLE9BQU8sSUFBSSxPQUFPLENBQUMsQ0FBQyxPQUFPLEVBQUUsTUFBTSxFQUFFLEVBQUU7WUFDckMsTUFBTSxHQUFHLEdBQUcsS0FBSyxDQUFDLFNBQVMsQ0FBQyxDQUFDLElBQU8sRUFBRSxFQUFFO2dCQUN0QyxJQUFJO29CQUNGLE9BQU8sQ0FBQyxJQUFJLENBQUMsQ0FBQztvQkFDZCxPQUFPLENBQUMsSUFBSSxDQUFDLENBQUM7aUJBQ2Y7Z0JBQUMsT0FBTyxDQUFDLEVBQUU7b0JBQ1YsTUFBTSxDQUFDLENBQUMsQ0FBQyxDQUFDO2lCQUNYO3dCQUFTO29CQUNSLEdBQUcsQ0FBQyxXQUFXLEVBQUUsQ0FBQztvQkFDbEIsWUFBWSxDQUFDLE9BQU8sQ0FBQyxDQUFDO2lCQUN2QjtZQUNILENBQUMsQ0FBQyxDQUFDO1lBQ0gsTUFBTSxPQUFPLEdBQUcsVUFBVSxDQUFDLEdBQUcsRUFBRTtnQkFDOUIsR0FBRyxDQUFDLFdBQVcsRUFBRSxDQUFDO2dCQUNsQix3REFBd0Q7Z0JBQ3hELE1BQU0sQ0FBQyxNQUFNLENBQUMsQ0FBQztZQUNqQixDQUFDLEVBQUUsRUFBRSxDQUFDLENBQUM7WUFDUCxPQUFPLEVBQUUsQ0FBQztRQUNaLENBQUMsQ0FBQyxDQUFDO0lBQ0wsQ0FBQztDQUFBO0FBRUQsTUFBTSxVQUFnQixjQUFjLENBQUksS0FBb0IsRUFDMUQsT0FBbUMsRUFBRSxPQUFtQixFQUFFLEtBQWEsQ0FBQyxFQUFFLFNBQWlCLFNBQVM7O1FBRXBHLE9BQU8sSUFBSSxPQUFPLENBQUMsQ0FBQyxPQUFPLEVBQUUsTUFBTSxFQUFFLEVBQUU7WUFDckMsTUFBTSxHQUFHLEdBQUcsS0FBSyxDQUFDLFNBQVMsQ0FBQyxDQUFDLElBQU8sRUFBRSxFQUFFO2dCQUN0QyxPQUFPLENBQUMsSUFBSSxDQUFDLENBQUMsSUFBSSxDQUFDLEdBQUcsRUFBRTtvQkFDdEIsT0FBTyxDQUFDLElBQUksQ0FBQyxDQUFDO2dCQUNoQixDQUFDLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxDQUFDLEVBQUUsRUFBRTtvQkFDYixNQUFNLENBQUMsQ0FBQyxDQUFDLENBQUM7Z0JBQ1osQ0FBQyxDQUFDLENBQUMsT0FBTyxDQUFDLEdBQUcsRUFBRTtvQkFDZCxHQUFHLENBQUMsV0FBVyxFQUFFLENBQUM7b0JBQ2xCLFlBQVksQ0FBQyxPQUFPLENBQUMsQ0FBQztnQkFDeEIsQ0FBQyxDQUFDLENBQUM7WUFDTCxDQUFDLENBQUMsQ0FBQztZQUNILE1BQU0sT0FBTyxHQUFHLFVBQVUsQ0FBQyxHQUFHLEVBQUU7Z0JBQzlCLEdBQUcsQ0FBQyxXQUFXLEVBQUUsQ0FBQztnQkFDbEIsd0RBQXdEO2dCQUN4RCxNQUFNLENBQUMsTUFBTSxDQUFDLENBQUM7WUFDakIsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxDQUFDO1lBQ1AsT0FBTyxFQUFFLENBQUM7UUFDWixDQUFDLENBQUMsQ0FBQztJQUNMLENBQUM7Q0FBQTtBQUVELE1BQU0sVUFBVSxJQUFJLENBQUMsSUFBWSxFQUFFLElBQXdCLEVBQUUsT0FBcUI7SUFDaEYsSUFBSSxLQUFLLENBQUMsZUFBZSxDQUFDLElBQUksU0FBUztRQUNyQyxLQUFLLENBQUMsZUFBZSxDQUFDLEdBQUcsRUFBRSxDQUFDO0lBQzlCLElBQUksS0FBSyxDQUFDLGVBQWUsQ0FBQyxDQUFDLEtBQUssSUFBSSxTQUFTO1FBQzNDLEtBQUssQ0FBQyxlQUFlLENBQUMsQ0FBQyxLQUFLLEdBQUcsRUFBRSxDQUFDO0lBQ3BDLEtBQUssQ0FBQyxlQUFlLENBQUMsQ0FBQyxLQUFNLENBQUMsSUFBSSxDQUFDLElBQUksSUFBSSxDQUFDLGVBQWUsRUFBRSxJQUFJLEVBQUUsSUFBSSxFQUFFLE9BQU8sQ0FBQyxDQUFDLENBQUM7QUFDckYsQ0FBQztBQUVELGdGQUFnRjtBQUNoRixNQUFNLFVBQVUsTUFBTSxDQUFDLE1BQVcsRUFBRSxXQUFnQixJQUFJLEVBQUUsS0FBYztJQUN0RSxJQUFJLEtBQUs7UUFDUCxLQUFLLEdBQUcsR0FBRyxLQUFLLElBQUksQ0FBQzs7UUFDbEIsS0FBSyxHQUFHLEVBQUUsQ0FBQztJQUNoQixJQUFJLE1BQU0sS0FBSyxRQUFRO1FBQ3JCLE1BQU0sSUFBSSxLQUFLLENBQUMsR0FBRyxLQUFLLGFBQWEsUUFBUSxXQUFXLE1BQU0sR0FBRyxDQUFDLENBQUM7QUFDdkUsQ0FBQztBQUVELE1BQU0sVUFBVSxXQUFXLENBQUMsTUFBYyxFQUFFLFFBQWdCLEVBQUUsU0FBUyxHQUFHLEtBQUssRUFBRSxLQUFjO0lBQzdGLElBQUksQ0FBQyxNQUFNLEtBQUssTUFBTSxDQUFDLGlCQUFpQixJQUFJLFFBQVEsS0FBSyxNQUFNLENBQUMsaUJBQWlCLENBQUM7UUFDaEYsQ0FBQyxNQUFNLEtBQUssTUFBTSxDQUFDLGlCQUFpQixJQUFJLFFBQVEsS0FBSyxNQUFNLENBQUMsaUJBQWlCLENBQUM7UUFDOUUsQ0FBQyxNQUFNLEtBQUssTUFBTSxDQUFDLEdBQUcsSUFBSSxRQUFRLEtBQUssTUFBTSxDQUFDLEdBQUcsQ0FBQyxJQUFJLENBQUMsS0FBSyxDQUFDLE1BQU0sQ0FBQyxJQUFJLEtBQUssQ0FBQyxRQUFRLENBQUMsQ0FBQztRQUN4RixPQUFPO0lBQ1QsTUFBTSxRQUFRLEdBQUcsSUFBSSxDQUFDLEdBQUcsQ0FBQyxNQUFNLEdBQUcsUUFBUSxDQUFDLEdBQUcsU0FBUyxDQUFDO0lBQ3pELE1BQU0sQ0FBQyxRQUFRLEVBQUUsSUFBSSxFQUFFLEdBQUcsS0FBSyxhQUFMLEtBQUssY0FBTCxLQUFLLEdBQUksRUFBRSxpQkFBaUIsU0FBUyxTQUFTLE1BQU0sU0FBUyxRQUFRLEdBQUcsQ0FBQyxDQUFDO0lBQ3BHLElBQUksQ0FBQyxRQUFRO1FBQ1gsTUFBTSxJQUFJLEtBQUssQ0FBQyxZQUFZLFFBQVEsU0FBUyxNQUFNLGlCQUFpQixTQUFTLEdBQUcsQ0FBQyxDQUFDO0FBQ3RGLENBQUM7QUFFRCxNQUFNLFVBQVUsV0FBVyxDQUFDLE1BQXFCLEVBQUUsUUFBdUIsRUFBRSxLQUFjO0lBQ3hGLE1BQU0sZ0JBQWdCLEdBQUcsUUFBUSxDQUFDLFFBQVEsQ0FBQztJQUMzQyxNQUFNLGNBQWMsR0FBRyxNQUFNLENBQUMsUUFBUSxDQUFDO0lBQ3ZDLE1BQU0sQ0FBQyxjQUFjLEVBQUUsZ0JBQWdCLEVBQUUsR0FBRyxLQUFLLGFBQUwsS0FBSyxjQUFMLEtBQUssR0FBSSxFQUFFLGFBQWEsQ0FBQyxDQUFDO0lBRXRFLEtBQUssTUFBTSxNQUFNLElBQUksUUFBUSxDQUFDLE9BQU8sRUFBRTtRQUNyQyxNQUFNLFlBQVksR0FBRyxNQUFNLENBQUMsT0FBTyxDQUFDLE1BQU0sQ0FBQyxNQUFNLENBQUMsSUFBSSxDQUFDLENBQUM7UUFDeEQsSUFBSSxZQUFZLElBQUksSUFBSTtZQUN0QixNQUFNLElBQUksS0FBSyxDQUFDLFVBQVUsTUFBTSxDQUFDLElBQUksWUFBWSxDQUFDLENBQUM7UUFDckQsSUFBSSxZQUFZLENBQUMsSUFBSSxJQUFJLE1BQU0sQ0FBQyxJQUFJO1lBQ2xDLE1BQU0sSUFBSSxLQUFLLENBQUMsVUFBVSxNQUFNLENBQUMsSUFBSSxrQkFBa0IsTUFBTSxDQUFDLElBQUksUUFBUSxZQUFZLENBQUMsSUFBSSxFQUFFLENBQUMsQ0FBQztRQUNqRyxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsZ0JBQWdCLEVBQUUsQ0FBQyxFQUFFLEVBQUU7WUFDekMsTUFBTSxLQUFLLEdBQUcsTUFBTSxDQUFDLEdBQUcsQ0FBQyxDQUFDLENBQUMsQ0FBQztZQUM1QixNQUFNLFdBQVcsR0FBRyxZQUFZLENBQUMsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDO1lBQ3hDLElBQUksTUFBTSxDQUFDLElBQUksSUFBSSxFQUFFLENBQUMsSUFBSSxDQUFDLEtBQUs7Z0JBQzlCLFdBQVcsQ0FBQyxXQUFXLEVBQUUsS0FBSyxFQUFFLE1BQU0sRUFBRSxLQUFLLENBQUMsQ0FBQztpQkFDNUMsSUFBSSxNQUFNLENBQUMsSUFBSSxJQUFJLEVBQUUsQ0FBQyxJQUFJLENBQUMsU0FBUztnQkFDdkMsTUFBTSxDQUFDLFdBQVcsQ0FBQyxNQUFNLENBQUMsS0FBSyxDQUFDLEVBQUUsSUFBSSxFQUFFLEtBQUssQ0FBQyxDQUFDOztnQkFFL0MsTUFBTSxDQUFDLFdBQVcsRUFBRSxLQUFLLEVBQUUsS0FBSyxDQUFDLENBQUM7U0FDckM7S0FDRjtBQUNILENBQUM7QUFFRCxNQUFNLFVBQVUsWUFBWSxDQUFDLE1BQThCLEVBQUUsUUFBZ0M7SUFDM0YsS0FBSyxNQUFNLENBQUMsV0FBVyxFQUFFLGFBQWEsQ0FBQyxJQUFJLE1BQU0sQ0FBQyxPQUFPLENBQUMsUUFBUSxDQUFDLEVBQUU7UUFDbkUsSUFBSSxDQUFDLE1BQU0sQ0FBQyxjQUFjLENBQUMsV0FBVyxDQUFDO1lBQ3JDLE1BQU0sSUFBSSxLQUFLLENBQUMsc0JBQXNCLFdBQVcsYUFBYSxDQUFDLENBQUM7UUFFbEUsTUFBTSxXQUFXLEdBQUcsTUFBTSxDQUFDLFdBQVcsQ0FBQyxDQUFDO1FBQ3hDLElBQUksV0FBVyxZQUFZLEtBQUssSUFBSSxhQUFhLFlBQVksS0FBSztZQUNoRSxXQUFXLENBQUMsV0FBVyxFQUFFLGFBQWEsQ0FBQyxDQUFDO2FBQ3JDLElBQUksV0FBVyxZQUFZLE1BQU0sSUFBSSxhQUFhLFlBQVksTUFBTTtZQUN2RSxZQUFZLENBQUMsV0FBVyxFQUFFLGFBQWEsQ0FBQyxDQUFDO2FBQ3RDLElBQUksTUFBTSxDQUFDLFFBQVEsQ0FBQyxXQUFXLENBQUMsSUFBSSxNQUFNLENBQUMsUUFBUSxDQUFDLGFBQWEsQ0FBQztZQUNyRSxXQUFXLENBQUMsV0FBVyxFQUFFLGFBQWEsQ0FBQyxDQUFDO2FBQ3JDLElBQUksV0FBVyxJQUFJLGFBQWE7WUFDbkMsTUFBTSxJQUFJLEtBQUssQ0FBQyxhQUFhLGFBQWEsY0FBYyxXQUFXLFdBQVcsV0FBVyxHQUFHLENBQUMsQ0FBQztLQUNqRztBQUNILENBQUM7QUFFRCxNQUFNLFVBQVUsV0FBVyxDQUFDLE1BQXNCLEVBQUUsUUFBd0I7SUFDMUUsTUFBTSxZQUFZLEdBQUcsTUFBTSxDQUFDLE1BQU0sQ0FBQztJQUNuQyxNQUFNLGNBQWMsR0FBRyxRQUFRLENBQUMsTUFBTSxDQUFDO0lBRXZDLElBQUksWUFBWSxJQUFJLGNBQWMsRUFBRTtRQUNsQyxNQUFNLElBQUksS0FBSyxDQUFDLDBEQUEwRCxZQUFZLEdBQUc7WUFDdkYsZ0NBQWdDLGNBQWMsRUFBRSxDQUFDLENBQUM7S0FDckQ7SUFFRCxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsWUFBWSxFQUFFLENBQUMsRUFBRSxFQUFFO1FBQ3JDLElBQUksTUFBTSxDQUFDLENBQUMsQ0FBQyxZQUFZLEtBQUssSUFBSSxRQUFRLENBQUMsQ0FBQyxDQUFDLFlBQVksS0FBSztZQUM1RCxXQUFXLENBQUMsTUFBTSxDQUFDLENBQUMsQ0FBQyxFQUFFLFFBQVEsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO2FBQ2pDLElBQUksTUFBTSxDQUFDLENBQUMsQ0FBQyxZQUFZLE1BQU0sSUFBSSxRQUFRLENBQUMsQ0FBQyxDQUFDLFlBQVksTUFBTTtZQUNuRSxZQUFZLENBQUMsTUFBTSxDQUFDLENBQUMsQ0FBQyxFQUFFLFFBQVEsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO2FBQ2xDLElBQUksTUFBTSxDQUFDLENBQUMsQ0FBQyxJQUFJLFFBQVEsQ0FBQyxDQUFDLENBQUM7WUFDL0IsTUFBTSxJQUFJLEtBQUssQ0FBQyxZQUFZLFFBQVEsQ0FBQyxDQUFDLENBQUMsZ0JBQWdCLENBQUMsU0FBUyxNQUFNLENBQUMsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDO0tBQ2pGO0FBQ0gsQ0FBQztBQUVELDJCQUEyQjtBQUMzQixNQUFNLFVBQVUsUUFBUSxDQUFDLFFBQWdCLEVBQUUsTUFBa0IsRUFBRSxPQUF5Qjs7SUFDdEYsZUFBZSxHQUFHLFFBQVEsQ0FBQztJQUMzQixNQUFNLEVBQUUsQ0FBQztJQUNULElBQUksS0FBSyxDQUFDLGVBQWUsQ0FBQyxFQUFFO1FBQzFCLEtBQUssQ0FBQyxlQUFlLENBQUMsQ0FBQyxLQUFLLEdBQUcsTUFBQSxPQUFPLGFBQVAsT0FBTyx1QkFBUCxPQUFPLENBQUUsS0FBSyxtQ0FBSSxJQUFJLENBQUM7UUFDdEQsS0FBSyxDQUFDLGVBQWUsQ0FBQyxDQUFDLE9BQU8sR0FBRyxPQUFPLGFBQVAsT0FBTyx1QkFBUCxPQUFPLENBQUUsT0FBTyxDQUFDO1FBQ2xELEtBQUssQ0FBQyxlQUFlLENBQUMsQ0FBQyxVQUFVLEdBQUcsT0FBTyxhQUFQLE9BQU8sdUJBQVAsT0FBTyxDQUFFLFVBQVUsQ0FBQztRQUN4RCxLQUFLLENBQUMsZUFBZSxDQUFDLENBQUMsV0FBVyxHQUFHLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxXQUFXLENBQUM7UUFDMUQsS0FBSyxDQUFDLGVBQWUsQ0FBQyxDQUFDLEtBQUssR0FBRyxPQUFPLGFBQVAsT0FBTyx1QkFBUCxPQUFPLENBQUUsS0FBSyxDQUFDO0tBQy9DO0FBQ0gsQ0FBQztBQUVELHVGQUF1RjtBQUN2RixNQUFNLFVBQVUsTUFBTSxDQUFDLE1BQTJCO0lBQ2hELElBQUksS0FBSyxDQUFDLGVBQWUsQ0FBQyxJQUFJLFNBQVM7UUFDckMsS0FBSyxDQUFDLGVBQWUsQ0FBQyxHQUFHLEVBQUUsQ0FBQztJQUM5QixLQUFLLENBQUMsZUFBZSxDQUFDLENBQUMsTUFBTSxHQUFHLE1BQU0sQ0FBQztBQUN6QyxDQUFDO0FBRUQsc0ZBQXNGO0FBQ3RGLE1BQU0sVUFBVSxLQUFLLENBQUMsS0FBMEI7SUFDOUMsSUFBSSxLQUFLLENBQUMsZUFBZSxDQUFDLElBQUksU0FBUztRQUNyQyxLQUFLLENBQUMsZUFBZSxDQUFDLEdBQUcsRUFBRSxDQUFDO0lBQzlCLEtBQUssQ0FBQyxlQUFlLENBQUMsQ0FBQyxLQUFLLEdBQUcsS0FBSyxDQUFDO0FBQ3ZDLENBQUM7QUFFRCxTQUFTLFlBQVksQ0FBQyxDQUFTLEVBQUUsQ0FBVztJQUMxQyxPQUFPLENBQUMsQ0FBQyxPQUFPLENBQUMsSUFBSSxNQUFNLENBQUMsQ0FBQyxDQUFDLElBQUksRUFBRSxJQUFJLENBQUMsRUFBRSxDQUFDLENBQUMsTUFBTSxDQUFDLENBQUM7QUFDdkQsQ0FBQztBQUVELE1BQU0sVUFBZ0IsYUFBYSxDQUFDLFFBQXFCLEVBQUUsTUFBWTs7O1FBQ3JFLE1BQU0sU0FBUyxHQUFHLFFBQVEsQ0FBQyxFQUFFLENBQUM7UUFDOUIsSUFBSSxhQUFhLENBQUMsU0FBUyxDQUFDO1lBQUUsT0FBTztRQUNyQyxNQUFNLFdBQVcsR0FBRyxNQUFNLENBQUMsQ0FBQyxDQUFDLE1BQU0sQ0FBQyxLQUFLLENBQUMsQ0FBQyxDQUFDLEtBQUssQ0FBQztRQUNsRCxJQUFJLFFBQVEsQ0FBQyxJQUFJLEtBQUssVUFBVSxJQUFJLENBQUMsQ0FBQyxDQUFDLE1BQU0sSUFBSSxNQUFNLENBQUMsUUFBUSxDQUFDLElBQUksS0FBSyxVQUFVLENBQUMsRUFBRTtZQUNyRixLQUFLLE1BQU0sQ0FBQyxJQUFVLE1BQU8sQ0FBQyxTQUFTLEVBQUU7Z0JBQ3ZDLE1BQU0sR0FBRyxHQUFHLENBQUMsQ0FBQyxJQUFJLENBQUMsS0FBSyxDQUFDLFlBQVksQ0FBQyxDQUFDO2dCQUN2QyxJQUFJLElBQUksR0FBRyxNQUFBLEdBQUcsQ0FBQyxHQUFHLEVBQUUsbUNBQUksQ0FBQyxDQUFDLElBQUksQ0FBQztnQkFDL0IsSUFBSSxHQUFHLEdBQUcsR0FBRyxDQUFDLE1BQU0sQ0FBQyxDQUFDLENBQUMsV0FBVyxHQUFHLElBQUksR0FBRyxHQUFHLENBQUMsSUFBSSxDQUFDLElBQUksQ0FBQyxDQUFDLENBQUMsQ0FBQyxXQUFXLENBQUM7Z0JBQ3pFLElBQUksUUFBUSxHQUFhLElBQUksQ0FBQyxLQUFLLENBQUMsS0FBSyxDQUFDLENBQUM7Z0JBQzNDLElBQUksR0FBRyxRQUFRLENBQUMsUUFBUSxDQUFDLE1BQU0sR0FBRyxDQUFDLENBQUMsQ0FBQztnQkFDckMsUUFBUSxDQUFDLE9BQU8sQ0FBQyxHQUFHLENBQUMsQ0FBQztnQkFDdEIsUUFBUSxDQUFDLEdBQUcsRUFBRSxDQUFDO2dCQUNmLEdBQUcsR0FBRyxRQUFRLENBQUMsSUFBSSxDQUFDLElBQUksQ0FBQyxDQUFDO2dCQUMxQixJQUFJLFdBQVcsQ0FBQyxHQUFHLENBQUMsS0FBSyxTQUFTO29CQUNoQyxXQUFXLENBQUMsR0FBRyxDQUFDLEdBQUcsRUFBRSxLQUFLLEVBQUUsRUFBRSxFQUFFLEtBQUssRUFBRSxJQUFJLEVBQUUsQ0FBQztnQkFDaEQsV0FBVyxDQUFDLEdBQUcsQ0FBQyxDQUFDLEtBQUssQ0FBQyxJQUFJLENBQUMsSUFBSSxJQUFJLENBQUMsR0FBRyxFQUFFLElBQUksRUFBRSxDQUFDLENBQUMsSUFBSSxFQUFFLEVBQUUsWUFBWSxFQUFFLEtBQUssRUFBRSxPQUFPLEVBQUUsTUFBQSxNQUFBLENBQUMsQ0FBQyxPQUFPLDBDQUFFLE9BQU8sbUNBQUksZ0JBQWdCLEVBQUUsVUFBVSxFQUFFLE1BQUEsQ0FBQyxDQUFDLE9BQU8sMENBQUUsVUFBVSxFQUFFLEtBQUssRUFBRSxNQUFBLENBQUMsQ0FBQyxPQUFPLDBDQUFFLEtBQUssRUFBRSxTQUFTLEVBQUUsTUFBQSxNQUFBLENBQUMsQ0FBQyxPQUFPLDBDQUFFLFNBQVMsbUNBQUksS0FBSyxFQUFFLENBQUMsQ0FBQyxDQUFDO2FBQzFPO1NBQ0Y7UUFDRCxNQUFNLGVBQWUsR0FBRyxFQUFFLENBQUM7UUFDM0IsTUFBTSxVQUFVLEdBQUcsRUFBRSxDQUFDO1FBQ3RCLE1BQU0sZUFBZSxHQUFHLEVBQUUsQ0FBQztRQUMzQixNQUFNLGFBQWEsR0FBRyxNQUFNLElBQUksQ0FBQyxJQUFJLENBQUMsU0FBUyxDQUFDLE1BQU0sQ0FBQyxpQkFBaUIsU0FBUyxHQUFHLENBQUMsQ0FBQyxJQUFJLEVBQUUsQ0FBQztRQUM3RixNQUFNLEdBQUcsR0FBRyxJQUFJLE1BQU0sQ0FBQyxvRUFBb0UsQ0FBQyxDQUFDO1FBQzdGLEtBQUssTUFBTSxDQUFDLElBQUksYUFBYSxFQUFFO1lBQzdCLE1BQU0sS0FBSyxHQUFHLENBQUMsQ0FBQyxPQUFPLENBQUMsTUFBTSxDQUFDLENBQUM7WUFDaEMsTUFBTSxJQUFJLEdBQUcsQ0FBQyxDQUFDLE9BQU8sQ0FBQyxVQUFVLENBQUMsQ0FBQztZQUNuQyxJQUFJLENBQUMsS0FBSyxJQUFJLEtBQUssQ0FBQyxPQUFPLENBQUMsS0FBSyxDQUFDLElBQUksS0FBSyxDQUFDLE1BQU0sQ0FBQyxFQUFFO2dCQUNuRCxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsS0FBSyxDQUFDLE1BQU0sRUFBRSxDQUFDLEVBQUUsRUFBRTtvQkFDckMsTUFBTSxHQUFHLEdBQUksS0FBSyxDQUFDLENBQUMsQ0FBWSxDQUFDLFFBQVEsQ0FBQyxHQUFHLENBQUMsQ0FBQztvQkFDL0MsTUFBTSxHQUFHLEdBQWdHLEVBQUUsQ0FBQztvQkFDNUcsS0FBSyxDQUFDLElBQUksQ0FBQyxHQUFHLENBQUMsQ0FBQyxPQUFPLENBQUMsQ0FBQyxHQUFHLEVBQUUsRUFBRTt3QkFDOUIsSUFBSSxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUMsVUFBVSxDQUFDLE1BQU0sQ0FBQzs0QkFBRSxHQUFHLENBQUMsTUFBTSxDQUFDLEdBQUcsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDOzZCQUMvQyxJQUFJLEdBQUcsQ0FBQyxDQUFDLENBQUMsQ0FBQyxVQUFVLENBQUMsTUFBTSxDQUFDOzRCQUFFLEdBQUcsQ0FBQyxNQUFNLENBQUMsR0FBRyxRQUFRLENBQUMsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUM7NkJBQzlELElBQUksR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDLFVBQVUsQ0FBQyxLQUFLLENBQUM7NEJBQUUsR0FBRyxDQUFDLEtBQUssQ0FBQyxHQUFHLEdBQUcsQ0FBQyxDQUFDLENBQUMsQ0FBQzs2QkFDbEQsSUFBSSxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUMsVUFBVSxDQUFDLFNBQVMsQ0FBQzs0QkFBRSxHQUFHLENBQUMsU0FBUyxDQUFDLEdBQUcsUUFBUSxDQUFDLEdBQUcsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO29CQUMzRSxDQUFDLENBQUMsQ0FBQztvQkFDSCxNQUFNLElBQUksR0FBRyxJQUFJLElBQUksQ0FBQyxNQUFBLEdBQUcsQ0FBQyxHQUFHLG1DQUFJLGdCQUFnQixFQUFFLEtBQUssQ0FBQyxNQUFNLEtBQUssQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsSUFBSSxDQUFDLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxJQUFJLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxFQUFFLEdBQVMsRUFBRTt3QkFDaEgsTUFBTSxHQUFHLEdBQUcsTUFBTSxJQUFJLENBQUMsU0FBUyxDQUFDLElBQUksQ0FBQyxZQUFZLENBQUMsS0FBSyxDQUFDLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxDQUFDLENBQUM7d0JBQ2pFLElBQUksR0FBRyxDQUFDLElBQUk7NEJBQUUsTUFBTSxLQUFLLENBQUMsR0FBRyxDQUFDLElBQUksQ0FBQyxDQUFDO3dCQUNwQyw0Q0FBNEM7d0JBQzVDLElBQUksT0FBTyxHQUFHLEtBQUssU0FBUyxJQUFJLENBQUMsR0FBRzs0QkFBRSxNQUFNLFdBQVcsS0FBSyxDQUFDLENBQUMsQ0FBQyx3QkFBd0IsR0FBRyxFQUFFLENBQUM7b0JBQy9GLENBQUMsQ0FBQSxFQUFFLEVBQUUsVUFBVSxFQUFFLEdBQUcsQ0FBQyxJQUFJLEVBQUUsT0FBTyxFQUFFLEVBQUUsQ0FBQyxJQUFJLENBQUMsYUFBYSxDQUFDLENBQUMsQ0FBQyxNQUFBLEdBQUcsQ0FBQyxnQkFBZ0IsbUNBQUksaUJBQWlCLENBQUMsQ0FBQyxDQUFDLE1BQUEsR0FBRyxDQUFDLE9BQU8sbUNBQUksZ0JBQWdCLEVBQUUsQ0FBQyxDQUFDO29CQUMzSSxJQUFJLEdBQUcsQ0FBQyxHQUFHLEVBQUU7d0JBQ1gsTUFBTSxHQUFHLEdBQVcsR0FBRyxDQUFDLEdBQUcsQ0FBQzt3QkFDNUIsSUFBSSxXQUFXLENBQUMsR0FBRyxDQUFDLEtBQUssU0FBUzs0QkFDaEMsV0FBVyxDQUFDLEdBQUcsQ0FBQyxHQUFHLEVBQUUsS0FBSyxFQUFFLEVBQUUsRUFBRSxLQUFLLEVBQUUsSUFBSSxFQUFFLENBQUM7d0JBRWhELHdFQUF3RTt3QkFDeEUsSUFBSSxDQUFDLFdBQVcsQ0FBQyxHQUFHLENBQUMsQ0FBQyxLQUFLOzRCQUN6QixXQUFXLENBQUMsR0FBRyxDQUFDLENBQUMsS0FBSyxHQUFHLEVBQUUsQ0FBQzt3QkFDOUIsV0FBVyxDQUFDLEdBQUcsQ0FBQyxDQUFDLEtBQUssQ0FBQyxJQUFJLENBQUMsSUFBSSxDQUFDLENBQUM7cUJBQ25DOzt3QkFFQyxlQUFlLENBQUMsSUFBSSxDQUFDLElBQUksQ0FBQyxDQUFDO2lCQUM5QjthQUNGO1lBQ0QsSUFBSSxJQUFJLEVBQUU7Z0JBQ1IsTUFBTSxJQUFJLEdBQUcsQ0FBQyxDQUFDLE9BQU8sQ0FBQyxVQUFVLENBQUMsQ0FBQyxDQUFDLENBQUMsUUFBUSxDQUFDLENBQUMsQ0FBQyxPQUFPLENBQUMsVUFBVSxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsU0FBUyxDQUFDO2dCQUNqRixNQUFNLElBQUksR0FBRyxJQUFJLElBQUksQ0FBQyxXQUFXLEVBQUUsQ0FBQyxDQUFDLFlBQVksRUFBRSxHQUFTLEVBQUU7b0JBQzVELE1BQU0sS0FBSyxDQUFDLEdBQUcsQ0FBQyxDQUFDO29CQUNqQixJQUFJLENBQUMsS0FBSyxDQUFDLGNBQWMsRUFBRSxDQUFDO29CQUM1QixNQUFNLENBQUMsQ0FBQyxLQUFLLEVBQUUsQ0FBQztvQkFDaEIsTUFBTSxLQUFLLENBQUMsSUFBSSxDQUFDLENBQUMsQ0FBQyxJQUFJLENBQUMsQ0FBQyxDQUFDLElBQUksQ0FBQyxDQUFDO29CQUNoQyxNQUFNLFNBQVMsR0FBRyxNQUFNLElBQUksQ0FBQyxLQUFLLENBQUMsU0FBUyxDQUFDO29CQUM3QyxJQUFJLFNBQVM7d0JBQ1gsTUFBTSxJQUFJLEtBQUssQ0FBQyxTQUFTLENBQUMsQ0FBQztnQkFDL0IsQ0FBQyxDQUFBLEVBQUUsRUFBRSxVQUFVLEVBQUUsQ0FBQyxDQUFDLE9BQU8sQ0FBQyxVQUFVLENBQUMsRUFBRSxDQUFDLENBQUM7Z0JBQzFDLFVBQVUsQ0FBQyxJQUFJLENBQUMsSUFBSSxDQUFDLENBQUM7YUFDdkI7WUFDRCxJQUFJLENBQUMsQ0FBQyxNQUFNLENBQUMsaUJBQWlCLENBQUMsRUFBRTtnQkFDL0IsSUFBSSxpQkFBaUIsR0FBRyxRQUFRLENBQUM7Z0JBQ2pDLElBQUksQ0FBQyxDQUFDLE9BQU8sQ0FBQyxVQUFVLENBQUMsRUFBRTtvQkFDekIsaUJBQWlCLEdBQUcsTUFBTSxJQUFJLENBQUMsSUFBSSxDQUFDLEtBQUssQ0FBQyxTQUFTLENBQUMsa0JBQWtCLFFBQVEsQ0FBQyxNQUFNLElBQUksQ0FBQyxDQUFDLE9BQU8sQ0FBQyxVQUFVLENBQUMsRUFBRSxDQUFDLENBQUM7aUJBQ25IO2dCQUVELE1BQU0sSUFBSSxHQUFHLElBQUksSUFBSSxDQUFDLGdCQUFnQixFQUFFLENBQUMsQ0FBQyxZQUFZLEVBQUUsR0FBUyxFQUFFO29CQUNqRSxNQUFNLEdBQUcsR0FBRyxFQUFFLENBQUM7b0JBQ2YsT0FBTyxDQUFDLEdBQUcsQ0FBQyxrQkFBa0IsUUFBUSxDQUFDLE1BQU0sSUFBSSxDQUFDLENBQUMsT0FBTyxDQUFDLFVBQVUsQ0FBQyxFQUFFLENBQUMsQ0FBQztvQkFFMUUsS0FBSyxNQUFNLEdBQUcsSUFBSSxpQkFBaUIsQ0FBQyxLQUFLLEVBQUUsQ0FBQyxPQUFPLEVBQUU7d0JBQ25ELE1BQU0sR0FBRyxHQUFHLE1BQU0sQ0FBQyxDQUFDLEtBQUssQ0FBQyxDQUFDLEdBQUcsQ0FBQyxDQUFDLENBQUM7d0JBQ2pDLEdBQUcsQ0FBQyxJQUFJLENBQUMsR0FBRyxJQUFJLEdBQUcsQ0FBQyxPQUFPLENBQUMsQ0FBQztxQkFDOUI7b0JBQ0QsTUFBTSxNQUFNLEdBQUcsR0FBRyxDQUFDLE1BQU0sQ0FBQyxDQUFDLENBQUMsRUFBRSxFQUFFLENBQUMsQ0FBQyxDQUFDLENBQUM7b0JBQ3BDLE1BQU0sQ0FBQyxNQUFNLENBQUMsTUFBTSxFQUFFLENBQUMsQ0FBQyxDQUFDO29CQUV6QixJQUFJLENBQUMsQ0FBQyxPQUFPLENBQUMsb0JBQW9CLENBQUM7d0JBQ2pDLE1BQU0sQ0FBQyxNQUFNLENBQUMsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDLE9BQU8sQ0FBQyxvQkFBb0IsQ0FBQyxDQUFDLENBQUM7Z0JBRXZELENBQUMsQ0FBQSxFQUFFLEVBQUUsVUFBVSxFQUFFLENBQUMsQ0FBQyxPQUFPLENBQUMsVUFBVSxDQUFDLEVBQUUsQ0FBQyxDQUFDO2dCQUMxQyxlQUFlLENBQUMsSUFBSSxDQUFDLElBQUksQ0FBQyxDQUFDO2FBQzVCO1NBQ0Y7UUFDRCxhQUFhLENBQUMsU0FBUyxDQUFDLEdBQUcsSUFBSSxDQUFDO1FBQ2hDLElBQUksZUFBZSxDQUFDLE1BQU0sR0FBRyxDQUFDO1lBQzVCLFdBQVcsQ0FBQyxnQkFBZ0IsQ0FBQyxHQUFHLEVBQUUsS0FBSyxFQUFFLGVBQWUsRUFBRSxLQUFLLEVBQUUsSUFBSSxFQUFFLENBQUM7UUFDMUUsSUFBSSxVQUFVLENBQUMsTUFBTSxHQUFHLENBQUM7WUFDdkIsV0FBVyxDQUFDLFdBQVcsQ0FBQyxHQUFHLEVBQUUsS0FBSyxFQUFFLFVBQVUsRUFBRSxLQUFLLEVBQUUsSUFBSSxFQUFFLENBQUM7UUFDaEUsSUFBSSxlQUFlLENBQUMsTUFBTSxHQUFHLENBQUM7WUFDNUIsV0FBVyxDQUFDLGdCQUFnQixDQUFDLEdBQUcsRUFBRSxLQUFLLEVBQUUsZUFBZSxFQUFFLEtBQUssRUFBRSxLQUFLLEVBQUUsQ0FBQzs7Q0FDNUU7QUFFRCxTQUFTLGVBQWU7SUFDdEIsTUFBTSxJQUFJLEdBQVUsRUFBRSxDQUFDO0lBQ3ZCLE9BQU8sQ0FBQyxHQUFHLEdBQUcsQ0FBQyxHQUFHLElBQUksRUFBRSxFQUFFO1FBQ3hCLElBQUksQ0FBQyxJQUFJLENBQUMsR0FBRyxJQUFJLENBQUMsQ0FBQztRQUNuQixNQUFNLENBQUMsR0FBRyxJQUFJLENBQUMsQ0FBQztJQUNsQixDQUFDLENBQUM7SUFDRixPQUFPLENBQUMsSUFBSSxHQUFHLENBQUMsR0FBRyxJQUFJLEVBQUUsRUFBRTtRQUN6QixJQUFJLENBQUMsSUFBSSxDQUFDLEdBQUcsSUFBSSxDQUFDLENBQUM7UUFDbkIsT0FBTyxDQUFDLEdBQUcsSUFBSSxDQUFDLENBQUM7SUFDbkIsQ0FBQyxDQUFDO0lBQ0YsT0FBTyxDQUFDLElBQUksR0FBRyxDQUFDLEdBQUcsSUFBSSxFQUFFLEVBQUU7UUFDekIsSUFBSSxDQUFDLElBQUksQ0FBQyxHQUFHLElBQUksQ0FBQyxDQUFDO1FBQ25CLE9BQU8sQ0FBQyxHQUFHLElBQUksQ0FBQyxDQUFDO0lBQ25CLENBQUMsQ0FBQztJQUNGLE9BQU8sQ0FBQyxLQUFLLEdBQUcsQ0FBQyxHQUFHLElBQUksRUFBRSxFQUFFO1FBQzFCLElBQUksQ0FBQyxJQUFJLENBQUMsR0FBRyxJQUFJLENBQUMsQ0FBQztRQUNuQixRQUFRLENBQUMsR0FBRyxJQUFJLENBQUMsQ0FBQztJQUNwQixDQUFDLENBQUM7SUFDRixPQUFPLElBQUksQ0FBQztBQUNkLENBQUM7QUFFRCxTQUFTLFlBQVk7SUFDbkIsT0FBTyxDQUFDLEdBQUcsR0FBRyxNQUFNLENBQUM7SUFDckIsT0FBTyxDQUFDLElBQUksR0FBRyxPQUFPLENBQUM7SUFDdkIsT0FBTyxDQUFDLElBQUksR0FBRyxPQUFPLENBQUM7SUFDdkIsT0FBTyxDQUFDLEtBQUssR0FBRyxRQUFRLENBQUM7QUFDM0IsQ0FBQztBQUVELE1BQU0sVUFBZ0IsUUFBUSxDQUFDLE9BQThCOzs7O1FBRTNELE1BQU0sUUFBUSxHQUFnQixDQUFBLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxXQUFXLEVBQUMsQ0FBQyxDQUFDLE9BQU8sQ0FBQyxXQUFXLENBQUMsT0FBTyxDQUFDLENBQUMsQ0FBQyxJQUFJLENBQUMsU0FBUyxDQUFDLGNBQWMsRUFBRSxDQUFDLElBQUksQ0FBQyxPQUFPLENBQUM7UUFDaEksSUFBSSxDQUFDLFFBQVE7WUFDWCxNQUFNLElBQUksS0FBSyxDQUFDLHlDQUF5QyxDQUFDLENBQUM7UUFDN0QsTUFBTSxLQUFLLEdBQUcsTUFBQSxRQUFRLENBQUMsWUFBWSwwQ0FBRSxLQUFLLENBQUMsV0FBVyxDQUFDLENBQUM7UUFDeEQsTUFBTSxZQUFZLEdBQUcsS0FBSyxDQUFDLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLEVBQUUsQ0FBQztRQUMzQyxJQUFJLFFBQVEsSUFBSSxTQUFTO1lBQ3ZCLE1BQU0sYUFBYSxDQUFDLFFBQVEsQ0FBQyxDQUFDO1FBQ2hDLE1BQU0sT0FBTyxHQUF3QixFQUFFLENBQUM7UUFDeEMsT0FBTyxDQUFDLEdBQUcsQ0FBQyxrQkFBa0IsQ0FBQyxDQUFDO1FBQ2hDLE9BQU8sQ0FBQyxHQUFHLENBQUMsT0FBTyxDQUFDLENBQUM7UUFDckIsT0FBTyxhQUFQLE9BQU8sY0FBUCxPQUFPLElBQVAsT0FBTyxHQUFLLEVBQUUsRUFBQztRQUNmLFlBQUEsT0FBUSxFQUFDLFdBQVcsdUNBQVgsV0FBVyxHQUFLLElBQUksV0FBVyxFQUFFLEVBQUM7UUFDM0MsSUFBSSxDQUFDLEtBQUssQ0FBQyxjQUFjLEVBQUUsQ0FBQztRQUM1QixNQUFNLElBQUksR0FBRyxlQUFlLEVBQUUsQ0FBQztRQUUvQixNQUFNLFdBQVcsQ0FBQyxLQUFLLEVBQUUsT0FBTyxDQUFDLENBQUM7UUFFbEMsS0FBSyxJQUFJLENBQUMsSUFBSSxPQUFPLEVBQUU7WUFDckIsQ0FBQyxDQUFDLE1BQU0sR0FBRyxDQUFDLENBQUMsTUFBTSxDQUFDLFFBQVEsRUFBRSxDQUFDLE9BQU8sQ0FBQyxJQUFJLEVBQUUsSUFBSSxDQUFDLENBQUM7WUFDbkQsSUFBSSxDQUFDLENBQUMsSUFBSSxJQUFJLFNBQVM7Z0JBQ3JCLENBQUMsQ0FBQyxJQUFJLEdBQUcsQ0FBQyxDQUFDLElBQUssQ0FBQyxRQUFRLEVBQUUsQ0FBQyxPQUFPLENBQUMsSUFBSSxFQUFFLElBQUksQ0FBQyxDQUFDO1NBQ25EO1FBQ0QsT0FBTyxPQUFPLENBQUM7UUFFZixTQUFlLG9CQUFvQixDQUFDLE1BQXlDLEVBQUUsUUFBZ0I7O2dCQUM3RixJQUFJLGdCQUFnQixHQUFHLFNBQVMsQ0FBQztnQkFDakMsSUFBSTtvQkFDRixJQUFJLE1BQU0sS0FBSyxTQUFTLEVBQUU7d0JBQ3hCLE1BQU0sT0FBTyxDQUFDLEdBQVMsRUFBRTs0QkFDdkIsTUFBTSxNQUFNLEVBQUUsQ0FBQzt3QkFDakIsQ0FBQyxDQUFBLEVBQUUsTUFBTSxFQUFFLFVBQVUsUUFBUSxpQkFBaUIsQ0FBQyxDQUFDO3FCQUNqRDtpQkFDRjtnQkFBQyxPQUFPLENBQU0sRUFBRTtvQkFDZixnQkFBZ0IsR0FBRyxNQUFNLFNBQVMsQ0FBQyxDQUFDLENBQUMsQ0FBQztpQkFDdkM7Z0JBQ0QsT0FBTyxnQkFBZ0IsQ0FBQTtZQUN6QixDQUFDO1NBQUE7UUFFRCxTQUFlLHFCQUFxQixDQUFDLFFBQWtCLEVBQUUsT0FBNkIsRUFBRSxnQkFBeUI7OztnQkFDL0csSUFBSSxDQUFDLEdBQUcsTUFBQSxRQUFRLENBQUMsS0FBSyxtQ0FBSSxFQUFFLENBQUM7Z0JBQzdCLE1BQU0sR0FBRyxHQUEwQixFQUFFLENBQUM7Z0JBQ3RDLGdGQUFnRjtnQkFDaEYsTUFBTSxhQUFhLEdBQUcsbUJBQW1CLEVBQUUsQ0FBQztnQkFFNUMsSUFBSSxRQUFRLENBQUMsS0FBSyxFQUFFO29CQUNoQixJQUFJLGFBQWEsR0FBRyxnQkFBZ0IsSUFBSSxPQUFPLENBQUMsVUFBVSxJQUFJLFNBQVMsQ0FBQztvQkFDMUUsS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLENBQUMsQ0FBQyxNQUFNLEVBQUUsQ0FBQyxFQUFFLEVBQUU7d0JBRWpDLElBQUksQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLE9BQU8sRUFBRTs0QkFDaEIsSUFBSSxDQUFBLE1BQUEsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLE9BQU8sMENBQUUsU0FBUyxNQUFLLFNBQVMsRUFBRTtnQ0FDekMsSUFBSSxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQyxPQUFPO29DQUNmLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQyxPQUFPLEdBQUcsRUFBRSxDQUFBO2dDQUNuQixDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsT0FBUSxDQUFDLFNBQVMsR0FBRyxNQUFBLFFBQVEsQ0FBQyxVQUFVLG1DQUFJLEtBQUssQ0FBQzs2QkFDeEQ7eUJBQ0Y7d0JBQ0QsSUFBSSxJQUFJLEdBQUcsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO3dCQUNoQixJQUFJLE9BQU8sQ0FBQyxJQUFJOzRCQUNkLElBQUksT0FBTyxDQUFDLElBQUksQ0FBQyxXQUFXLEVBQUUsS0FBSyxJQUFJLENBQUMsSUFBSSxDQUFDLFdBQVcsRUFBRTtnQ0FDeEQsU0FBUzt3QkFDYixJQUFJLGFBQWEsRUFBRTs0QkFDakIsSUFBSSxDQUFBLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxVQUFVLEtBQUksU0FBUyxJQUFJLElBQUksQ0FBQyxJQUFJLENBQUMsV0FBVyxFQUFFLENBQUMsSUFBSSxFQUFFLE1BQUssT0FBTyxhQUFQLE9BQU8sdUJBQVAsT0FBTyxDQUFFLFVBQVUsQ0FBQyxXQUFXLEdBQUcsSUFBSSxFQUFFLENBQUEsRUFBRTtnQ0FDbkgsc0RBQXNEO2dDQUN0RCxhQUFhLEdBQUcsS0FBSyxDQUFDOzZCQUN2Qjs7Z0NBQ0QsU0FBUzt5QkFDVjt3QkFDRCxJQUFJLElBQUksYUFBSixJQUFJLHVCQUFKLElBQUksQ0FBRSxPQUFPLEVBQUU7NEJBQ2pCLElBQUksQ0FBQyxPQUFPLENBQUMsS0FBSyxHQUFHLE1BQUEsTUFBQSxNQUFBLE1BQUEsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLE9BQU8sMENBQUUsS0FBSyxtQ0FBSSxRQUFRLGFBQVIsUUFBUSx1QkFBUixRQUFRLENBQUUsS0FBSyxtQ0FBSSxZQUFZLG1DQUFJLEVBQUUsQ0FBQzt5QkFDbkY7d0JBQ0QsZ0ZBQWdGO3dCQUNoRix3Q0FBd0M7d0JBQ3hDLGtCQUFrQjt3QkFDbEIsZ0NBQWdDO3dCQUNoQyw0RUFBNEU7d0JBQzVFLElBQUksT0FBTyxHQUFHLE1BQU0sUUFBUSxDQUN4QixJQUFJLEVBQ0osT0FBTyxhQUFQLE9BQU8sdUJBQVAsT0FBTyxDQUFFLElBQUksRUFDYixJQUFJLEVBQUUsRUFBRSxDQUFDLElBQUksQ0FBQyxhQUFhLENBQUMsQ0FBQyxDQUFDLE1BQUEsTUFBQSxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsT0FBTywwQ0FBRSxnQkFBZ0IsbUNBQUksaUJBQWlCLENBQUMsQ0FBQyxDQUFDLE1BQUEsTUFBQSxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsT0FBTywwQ0FBRSxPQUFPLG1DQUFJLGdCQUFnQixFQUM3SCxRQUFRLENBQUMsSUFBSSxFQUNiLE9BQU8sQ0FBQyxPQUFPLENBQ2xCLENBQUM7d0JBRUYsa0JBQWtCO3dCQUNsQixnQ0FBZ0M7d0JBQ2hDLElBQUksT0FBTyxFQUFFOzRCQUNYLEdBQUcsQ0FBQyxJQUFJLGlDQUFNLE9BQU8sS0FBRyxpQkFBaUIsRUFBRSxtQkFBbUIsRUFBRSxHQUFHLGFBQWEsSUFBRyxDQUFDOzRCQUNwRiwwR0FBMEc7NEJBQzFHLElBQUksT0FBTyxDQUFDLFlBQVksSUFBSSxPQUFPLENBQUMsVUFBVSxLQUFLLElBQUksQ0FBQyxJQUFJLElBQUksQ0FBQyxPQUFPLENBQUMsT0FBTyxJQUFJLENBQUMsT0FBTyxDQUFDLE9BQU87Z0NBQ2xHLE9BQU8sR0FBRyxDQUFDO3lCQUNkO3dCQUNELHdLQUF3Szt3QkFFeEssSUFBSSxDQUFDLE9BQU8sQ0FBQyxXQUFXLEVBQUU7NEJBQ3hCLElBQUksQ0FBQyxLQUFLLENBQUMsUUFBUSxFQUFFLENBQUM7NEJBQ3RCLEVBQUUsQ0FBQyxPQUFPLENBQUMsUUFBUSxFQUFFLENBQUM7eUJBQ3ZCO3FCQUNGO2lCQUNGO3FCQUFNO29CQUNMLElBQUksYUFBYSxHQUFHLGdCQUFnQixJQUFJLE9BQU8sQ0FBQyxVQUFVLElBQUksU0FBUyxDQUFDO29CQUN4RSxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxDQUFDLE1BQU0sRUFBRSxDQUFDLEVBQUUsRUFBRTt3QkFDakMsSUFBSSxJQUFJLEdBQUcsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO3dCQUNoQixJQUFJLE9BQU8sQ0FBQyxJQUFJOzRCQUNkLElBQUksT0FBTyxDQUFDLElBQUksQ0FBQyxXQUFXLEVBQUUsS0FBSyxJQUFJLENBQUMsSUFBSSxDQUFDLFdBQVcsRUFBRTtnQ0FDeEQsU0FBUzt3QkFDYixJQUFJLGFBQWEsRUFBRTs0QkFDakIsSUFBSSxDQUFBLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxVQUFVLEtBQUksU0FBUyxJQUFJLElBQUksQ0FBQyxJQUFJLENBQUMsV0FBVyxFQUFFLENBQUMsSUFBSSxFQUFFLE1BQUssT0FBTyxhQUFQLE9BQU8sdUJBQVAsT0FBTyxDQUFFLFVBQVUsQ0FBQyxXQUFXLEdBQUcsSUFBSSxFQUFFLENBQUEsRUFBRTtnQ0FDbkgsc0RBQXNEO2dDQUN0RCxhQUFhLEdBQUcsS0FBSyxDQUFDOzZCQUN2Qjs0QkFDRCxTQUFTLENBQUUsd0NBQXdDO3lCQUNwRDt3QkFFRCxJQUFJLElBQUksYUFBSixJQUFJLHVCQUFKLElBQUksQ0FBRSxPQUFPLEVBQUU7NEJBQ2pCLElBQUksQ0FBQyxPQUFPLENBQUMsS0FBSyxHQUFHLE1BQUEsTUFBQSxNQUFBLE1BQUEsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLE9BQU8sMENBQUUsS0FBSyxtQ0FBSSxRQUFRLGFBQVIsUUFBUSx1QkFBUixRQUFRLENBQUUsS0FBSyxtQ0FBSSxZQUFZLG1DQUFJLEVBQUUsQ0FBQzt5QkFDbkY7d0JBQ0QsZ0ZBQWdGO3dCQUNoRix3Q0FBd0M7d0JBQ3hDLGtCQUFrQjt3QkFDbEIsZ0NBQWdDO3dCQUNoQyw0RUFBNEU7d0JBQzVFLElBQUksT0FBTyxHQUFHLE1BQU0sUUFBUSxDQUN4QixJQUFJLEVBQ0osT0FBTyxhQUFQLE9BQU8sdUJBQVAsT0FBTyxDQUFFLElBQUksRUFDYixJQUFJLEVBQ0osRUFBRSxDQUFDLElBQUksQ0FBQyxhQUFhLENBQUMsQ0FBQyxDQUFDLE1BQUEsTUFBQSxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsT0FBTywwQ0FBRSxnQkFBZ0IsbUNBQUksaUJBQWlCLENBQUMsQ0FBQyxDQUFDLE1BQUEsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLE9BQU8sMENBQUUsT0FBTyxFQUNuRyxRQUFRLENBQUMsSUFBSSxFQUNiLE9BQU8sQ0FBQyxPQUFPLENBQ2xCLENBQUM7d0JBRUYsa0JBQWtCO3dCQUNsQixnQ0FBZ0M7d0JBRWhDLElBQUksT0FBTyxFQUFFOzRCQUNYLEdBQUcsQ0FBQyxJQUFJLGlDQUFNLE9BQU8sS0FBRSxpQkFBaUIsRUFBRSxtQkFBbUIsRUFBRSxHQUFHLGFBQWEsSUFBRyxDQUFDOzRCQUNuRiwwR0FBMEc7NEJBQzFHLElBQUksT0FBTyxDQUFDLFlBQVksSUFBSSxPQUFPLENBQUMsVUFBVSxLQUFLLElBQUksQ0FBQyxJQUFJLElBQUksQ0FBQyxPQUFPLENBQUMsT0FBTyxJQUFJLENBQUMsT0FBTyxDQUFDLE9BQU87Z0NBQ2xHLE9BQU8sR0FBRyxDQUFDO3lCQUNkO3dCQUNELDZLQUE2SztxQkFFOUs7aUJBQ0Y7Z0JBQ0QsT0FBTyxHQUFHLENBQUM7O1NBQ1o7UUFFRCxTQUFTLG1CQUFtQjs7WUFDMUIsSUFBSSxPQUFPLE9BQU8sS0FBSyxXQUFXO2dCQUNoQyxPQUFPLENBQUMsQ0FBQztZQUNYLElBQUksTUFBTSxHQUFHLENBQUMsQ0FBQyxDQUFDO1lBQ2hCLElBQUk7Z0JBQ0YsTUFBTSxHQUFHLEVBQUUsQ0FBQyxNQUFNLENBQUMsTUFBTSxFQUFFLENBQUMsTUFBTSxDQUFDO2FBQ3BDO1lBQUMsT0FBTyxDQUFNLEVBQUU7Z0JBQ2YsT0FBTyxDQUFDLElBQUksQ0FBQyxNQUFBLENBQUMsQ0FBQyxPQUFPLG1DQUFJLENBQUMsQ0FBQyxDQUFDO2FBQzlCO1lBQ0QsT0FBTyxNQUFNLENBQUM7UUFDaEIsQ0FBQztRQUVELFNBQWUsV0FBVyxDQUFDLGtCQUErQyxFQUFFLE9BQTZCOzs7Z0JBQ3ZHLElBQUk7b0JBQ0YsSUFBSSxrQkFBa0IsR0FBRyxDQUFBLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxjQUFjLEtBQUksU0FBUyxDQUFDO29CQUM5RCxJQUFJLGdCQUFnQixHQUFHLEtBQUssQ0FBQztvQkFDN0IsS0FBSyxNQUFNLENBQUMsR0FBRyxFQUFFLEtBQUssQ0FBQyxJQUFJLE1BQU0sQ0FBQyxPQUFPLENBQUMsa0JBQWtCLENBQUMsRUFBRTt3QkFDM0QsSUFBSSxNQUFBLE9BQU8sQ0FBQyxPQUFPLDBDQUFFLElBQUksQ0FBQyxDQUFDLENBQUMsRUFBRSxFQUFFLENBQUMsR0FBRyxDQUFDLFVBQVUsQ0FBQyxDQUFDLENBQUMsQ0FBQzs0QkFDL0MsU0FBUzt3QkFDYixJQUFJLENBQUEsT0FBTyxhQUFQLE9BQU8sdUJBQVAsT0FBTyxDQUFFLFFBQVEsS0FBSSxJQUFJLElBQUksQ0FBQyxHQUFHLENBQUMsV0FBVyxFQUFFLENBQUMsVUFBVSxDQUFDLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxRQUFRLENBQUMsV0FBVyxHQUFHLElBQUksRUFBRSxDQUFDOzRCQUNsRyxTQUFTO3dCQUViLElBQUksa0JBQWtCLEVBQUU7NEJBQ3BCLElBQUksZ0JBQWdCO2dDQUNoQixrQkFBa0IsR0FBRyxLQUFLLENBQUM7aUNBQzFCO2dDQUNELElBQUksQ0FBQSxPQUFPLGFBQVAsT0FBTyx1QkFBUCxPQUFPLENBQUUsY0FBYyxLQUFJLElBQUksSUFBSSxHQUFHLENBQUMsV0FBVyxFQUFFLENBQUMsSUFBSSxFQUFFLE1BQUssT0FBTyxhQUFQLE9BQU8sdUJBQVAsT0FBTyxDQUFFLGNBQWMsQ0FBQyxXQUFXLEdBQUcsSUFBSSxFQUFFLENBQUEsRUFBRTtvQ0FDOUcsZ0JBQWdCLEdBQUcsSUFBSSxDQUFDO2lDQUMzQjtxQ0FBTTtvQ0FDSCx1REFBdUQ7b0NBQ3ZELFNBQVM7aUNBQ1o7NkJBQ0o7eUJBQ0o7d0JBQ0QsWUFBWTt3QkFDWixNQUFNLE9BQU8sR0FBRyxNQUFBLEtBQUssQ0FBQyxLQUFLLDBDQUFFLEtBQUssQ0FBQyxDQUFDLENBQU8sRUFBRSxFQUFFOzs0QkFBQyxPQUFBLENBQUEsTUFBQSxDQUFDLENBQUMsT0FBTywwQ0FBRSxVQUFVO21DQUM5RCxDQUFDLENBQUEsT0FBTyxhQUFQLE9BQU8sdUJBQVAsT0FBTyxDQUFFLElBQUksS0FBSSxJQUFJLElBQUksT0FBTyxDQUFDLElBQUksQ0FBQyxXQUFXLEVBQUUsS0FBSyxDQUFDLENBQUMsSUFBSSxDQUFDLFdBQVcsRUFBRSxDQUFDLENBQUE7eUJBQUEsQ0FBQyxDQUFDO3dCQUV2RixJQUFJLENBQUMsT0FBTyxFQUFFOzRCQUNWLFlBQVk7NEJBQ1osTUFBTSxZQUFZLEdBQUcsQ0FBQyxNQUFBLEtBQUssQ0FBQyxLQUFLLG1DQUFJLEVBQUUsQ0FBQyxDQUFDLE1BQU0sQ0FBQyxDQUFDLENBQU8sRUFBRSxFQUFFLFdBQzFELE9BQUEsQ0FBQSxNQUFBLENBQUMsQ0FBQyxPQUFPLDBDQUFFLFVBQVUsS0FBSSxDQUFDLENBQUEsT0FBTyxhQUFQLE9BQU8sdUJBQVAsT0FBTyxDQUFFLElBQUksS0FBSSxJQUFJLElBQUksT0FBTyxDQUFDLElBQUksQ0FBQyxXQUFXLEVBQUUsS0FBSyxDQUFDLENBQUMsSUFBSSxDQUFDLFdBQVcsRUFBRSxDQUFDLENBQUEsRUFBQSxDQUN4RyxDQUFDLE1BQU0sQ0FBQzs0QkFDVCxNQUFNLENBQUMsOEJBQThCLEdBQUcsS0FBSyxZQUFZLEdBQUcsQ0FBQyxDQUFDLENBQUMsQ0FBQyxjQUFjLFlBQVksSUFBSSxDQUFDLENBQUMsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxDQUFDOzRCQUN2RyxLQUFLLENBQUMsWUFBWSxHQUFHLE1BQU0sb0JBQW9CLENBQUMsS0FBSyxDQUFDLE1BQU0sRUFBRSxHQUFHLENBQUMsQ0FBQzt5QkFDdEU7d0JBQ0QsSUFBSSxDQUFDLEdBQUcsTUFBQSxLQUFLLENBQUMsS0FBSyxtQ0FBSSxFQUFFLENBQUM7d0JBRTFCLElBQUksT0FBTyxDQUFDLFVBQVUsRUFBRTs0QkFDcEIsQ0FBQyxHQUFHLENBQUMsQ0FBQyxNQUFNLENBQUMsQ0FBQyxDQUFDLEVBQUUsRUFBRSxXQUFDLE9BQUEsTUFBQSxDQUFDLENBQUMsT0FBTywwQ0FBRSxVQUFVLENBQUEsRUFBQSxDQUFDLENBQUM7NEJBQzNDLENBQUMsR0FBRyxPQUFPLENBQUMsQ0FBQyxDQUFDLENBQUM7eUJBQ2xCO3dCQUVELElBQUksQ0FBQyxNQUFBLE1BQUEsT0FBTyxDQUFDLElBQUksMENBQUUsTUFBTSxtQ0FBSSxDQUFDLENBQUMsR0FBRyxDQUFDLEVBQUU7NEJBQ2pDLENBQUMsR0FBRyxDQUFDLENBQUMsTUFBTSxDQUFDLENBQUMsQ0FBQyxFQUFFLEVBQUUsZUFDZixPQUFBLE1BQUEsTUFBQSxDQUFDLENBQUMsT0FBTywwQ0FBRSxJQUFJLDBDQUFFLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxXQUFDLE9BQUEsQ0FBQyxNQUFBLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxJQUFJLG1DQUFJLEVBQUUsQ0FBQyxDQUFDLFFBQVEsQ0FBQyxHQUFHLENBQUMsQ0FBQSxFQUFBLENBQUMsQ0FBQSxFQUFBLENBQ3BFLENBQUM7eUJBQ0w7d0JBRUQsSUFBSSxHQUF5QixDQUFDO3dCQUM5QixJQUFJLEtBQUssQ0FBQyxZQUFZLEVBQUU7NEJBQ3BCLEdBQUcsR0FBRyxLQUFLLENBQUMsSUFBSSxDQUFDLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxRQUFRLEVBQUUsRUFBRTtnQ0FDaEMsT0FBTztvQ0FDSCxJQUFJLEVBQUUsSUFBSSxJQUFJLEVBQUUsQ0FBQyxXQUFXLEVBQUU7b0NBQzlCLFFBQVEsRUFBRSxHQUFHO29DQUNiLElBQUksRUFBRSxRQUFRLENBQUMsSUFBSTtvQ0FDbkIsT0FBTyxFQUFFLEtBQUs7b0NBQ2QsTUFBTSxFQUFFLGlCQUFpQjtvQ0FDekIsRUFBRSxFQUFFLENBQUM7b0NBQ0wsT0FBTyxFQUFFLEtBQUs7b0NBQ2QsSUFBSSxFQUFFLEVBQUU7b0NBQ1IsS0FBSyxFQUFFLFlBQVk7b0NBQ25CLE9BQU8sRUFBRSxRQUFRLENBQUMsSUFBSTtvQ0FDdEIsaUJBQWlCLEVBQUUsQ0FBQztvQ0FDcEIsT0FBTyxFQUFFLEVBQUUsQ0FBQyxJQUFJLENBQUMsYUFBYTtpQ0FDakMsQ0FBQzs0QkFDTixDQUFDLENBQUMsQ0FBQyxDQUFDOzRCQUNKLEdBQUcsQ0FBQyxPQUFPLENBQUMsQ0FBTyxJQUFJLEVBQUUsRUFBRSxnREFBQyxPQUFBLE1BQU0sSUFBSSxDQUFDLEtBQUssQ0FBQyxVQUFVLENBQUMsU0FBUyxFQUFFLElBQUksQ0FBQyxDQUFBLEdBQUEsQ0FBQyxDQUFDO3lCQUM3RTs7NEJBQ0csR0FBRyxHQUFHLE1BQU0scUJBQXFCLENBQUMsS0FBSyxFQUFFLE9BQU8sRUFBRSxrQkFBa0IsQ0FBQyxDQUFDO3dCQUMxRSxNQUFNLElBQUksR0FBeUIsR0FBRyxDQUFDLE1BQU0sQ0FBQyxDQUFDLENBQUMsRUFBRSxFQUFFLENBQUMsQ0FBQyxDQUFDLE1BQU0sSUFBSSxTQUFTLENBQUMsQ0FBQzt3QkFFNUUsSUFBSSxDQUFDLE9BQU87NEJBQ1IsS0FBSyxDQUFDLFdBQVcsR0FBRyxNQUFNLG9CQUFvQixDQUFDLEtBQUssQ0FBQyxLQUFLLEVBQUUsR0FBRyxDQUFDLENBQUM7d0JBRXJFLHVCQUF1Qjt3QkFDdkIseUJBQXlCO3dCQUN6Qix5QkFBeUI7d0JBQ3pCLElBQUksS0FBSyxDQUFDLFdBQVcsRUFBRTs0QkFDbkIsTUFBTSxDQUFDLHVDQUF1QyxHQUFHLFdBQVcsQ0FBQyxDQUFDOzRCQUM5RCxNQUFNLENBQUMsaUNBQWlDLEdBQUcsYUFBYSxLQUFLLENBQUMsV0FBVyxFQUFFLENBQUMsQ0FBQzs0QkFDN0UsSUFBSSxDQUFDLElBQUksQ0FBQztnQ0FDTixJQUFJLEVBQUUsSUFBSSxJQUFJLEVBQUUsQ0FBQyxXQUFXLEVBQUU7Z0NBQzlCLFFBQVEsRUFBRSxHQUFHO2dDQUNiLElBQUksRUFBRSxPQUFPO2dDQUNiLE9BQU8sRUFBRSxLQUFLO2dDQUNkLE1BQU0sRUFBRSxLQUFLLENBQUMsV0FBVztnQ0FDekIsRUFBRSxFQUFFLENBQUM7Z0NBQ0wsT0FBTyxFQUFFLEtBQUs7Z0NBQ2QsSUFBSSxFQUFFLEVBQUU7Z0NBQ1IsS0FBSyxFQUFFLFlBQVk7Z0NBQ25CLE9BQU8sRUFBRSxRQUFRLENBQUMsSUFBSTtnQ0FDdEIsaUJBQWlCLEVBQUUsQ0FBQztnQ0FDcEIsT0FBTyxFQUFFLEVBQUUsQ0FBQyxJQUFJLENBQUMsYUFBYTs2QkFDakMsQ0FBQyxDQUFDO3lCQUNOO3dCQUNELElBQUksS0FBSyxDQUFDLFlBQVksRUFBRTs0QkFDcEIsTUFBTSxDQUFDLHdDQUF3QyxHQUFHLFdBQVcsQ0FBQyxDQUFDOzRCQUMvRCxNQUFNLENBQUMsaUNBQWlDLEdBQUcsY0FBYyxLQUFLLENBQUMsWUFBWSxFQUFFLENBQUMsQ0FBQzs0QkFDL0UsSUFBSSxDQUFDLElBQUksQ0FBQztnQ0FDTixJQUFJLEVBQUUsSUFBSSxJQUFJLEVBQUUsQ0FBQyxXQUFXLEVBQUU7Z0NBQzlCLFFBQVEsRUFBRSxHQUFHO2dDQUNiLElBQUksRUFBRSxRQUFRO2dDQUNkLE9BQU8sRUFBRSxLQUFLO2dDQUNkLE1BQU0sRUFBRSxLQUFLLENBQUMsWUFBWTtnQ0FDMUIsRUFBRSxFQUFFLENBQUM7Z0NBQ0wsT0FBTyxFQUFFLEtBQUs7Z0NBQ2QsSUFBSSxFQUFFLEVBQUU7Z0NBQ1IsS0FBSyxFQUFFLFlBQVk7Z0NBQ25CLE9BQU8sRUFBRSxRQUFRLENBQUMsSUFBSTtnQ0FDdEIsaUJBQWlCLEVBQUUsQ0FBQztnQ0FDcEIsT0FBTyxFQUFFLEVBQUUsQ0FBQyxJQUFJLENBQUMsYUFBYTs2QkFDakMsQ0FBQyxDQUFDO3lCQUNOO3dCQUNELE9BQU8sQ0FBQyxJQUFJLENBQUMsR0FBRyxJQUFJLENBQUMsQ0FBQzt3QkFFdEIsb0dBQW9HO3dCQUNwRyxJQUFJLE9BQU8sQ0FBQyxZQUFZLElBQUksSUFBSSxDQUFDLElBQUksQ0FBQyxDQUFDLENBQUMsRUFBRSxFQUFFLENBQUMsQ0FBQyxDQUFDLENBQUMsT0FBTyxJQUFJLENBQUMsQ0FBQyxDQUFDLE9BQU8sSUFBSSxDQUFDLENBQUMsSUFBSSxLQUFLLE9BQU8sQ0FBQyxVQUFVLENBQUM7NEJBQ25HLE1BQU07cUJBQ2I7aUJBQ0Y7d0JBQVM7b0JBQ1IsWUFBWSxFQUFFLENBQUM7aUJBQ2hCO2dCQUNELElBQUksT0FBTyxDQUFDLFdBQVksQ0FBQyxjQUFjLElBQUksQ0FBQyxDQUFDLEVBQUUsQ0FBQyxJQUFJLENBQUMsYUFBYSxDQUFDLEVBQUU7b0JBQ25FLE1BQU0sS0FBSyxDQUFDLElBQUksQ0FBQyxDQUFDO29CQUNsQixNQUFNLEtBQUssR0FBRyxNQUFNLElBQUksQ0FBQyxLQUFLLENBQUMsU0FBUyxDQUFDO29CQUN6QyxJQUFJLEtBQUssSUFBSSxTQUFTLEVBQUU7d0JBQ3BCLE1BQU0sTUFBTSxHQUFROzRCQUNoQixJQUFJLEVBQUUsRUFBRTs0QkFDUixJQUFJLEVBQUUsSUFBSSxJQUFJLEVBQUUsQ0FBQyxXQUFXLEVBQUU7NEJBQzlCLFFBQVEsRUFBRSxzQkFBc0I7NEJBQ2hDLElBQUksRUFBRSxXQUFXOzRCQUNqQixNQUFNLEVBQUUsS0FBSyxhQUFMLEtBQUssY0FBTCxLQUFLLEdBQUksRUFBRTs0QkFDbkIsT0FBTyxFQUFFLENBQUMsS0FBSzs0QkFDZixFQUFFLEVBQUUsQ0FBQzs0QkFDTCxPQUFPLEVBQUUsS0FBSzs0QkFDZCxLQUFLLEVBQUUsWUFBWSxhQUFaLFlBQVksY0FBWixZQUFZLEdBQUksRUFBRTs0QkFDekIsU0FBUyxFQUFFLFFBQVEsQ0FBQyxJQUFJOzRCQUN4QixpQkFBaUIsRUFBRSxDQUFDO3lCQUN2QixDQUFDO3dCQUNGLE1BQU0sQ0FBQyx5Q0FBeUMsS0FBSyxFQUFFLENBQUMsQ0FBQzt3QkFFekQsT0FBTyxDQUFDLElBQUksaUNBQUssTUFBTSxLQUFFLFNBQVMsRUFBRSxFQUFFLENBQUMsSUFBSSxDQUFDLGFBQWEsSUFBSSxDQUFDLEtBQUssSUFBRSxDQUFDO3dCQUNoRSxNQUFPLENBQUMsT0FBTyxHQUFHLFFBQVEsQ0FBQyxJQUFJLENBQUM7d0JBQ3RDLE1BQU0sSUFBSSxDQUFDLEtBQUssQ0FBQyxVQUFVLENBQUMsU0FBUyxFQUFFLE1BQU0sQ0FBQyxDQUFDO3FCQUNsRDtpQkFDRjs7U0FDRjs7Q0FDRjtBQUVELFNBQWUsU0FBUyxDQUFDLENBQU07O1FBQzdCLE9BQU8sR0FBRyxDQUFDLENBQUMsUUFBUSxFQUFFLEtBQUssQ0FBQyxDQUFDLEtBQUssQ0FBQyxDQUFDLENBQUMsQ0FBQyxNQUFNLEVBQUUsQ0FBQyxNQUFNLENBQUMsbUJBQW1CLENBQUMsQ0FBQyxDQUFDLEtBQUssQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLEVBQUUsRUFBRSxDQUFDO0lBQzdGLENBQUM7Q0FBQTtBQUVELFNBQWUsUUFBUSxDQUFDLENBQU8sRUFBRSxTQUE2QixFQUFFLElBQVcsRUFDekUsV0FBb0IsRUFBRSxXQUFvQixFQUFFLE9BQWlCOzs7UUFFN0QsSUFBSSxDQUFDLE1BQU0sR0FBRyxDQUFDLENBQUM7UUFDaEIsSUFBSSxDQUFhLENBQUM7UUFDbEIsSUFBSSxJQUFJLEdBQVcsU0FBUyxDQUFDO1FBQzdCLE1BQU0sTUFBTSxHQUFHLFNBQVMsSUFBSSxTQUFTLElBQUksQ0FBQyxDQUFDLENBQUMsSUFBSSxDQUFDLFdBQVcsRUFBRSxLQUFLLFNBQVMsQ0FBQyxXQUFXLEVBQUUsQ0FBQyxDQUFDO1FBQzVGLElBQUksSUFBSSxHQUFHLENBQUEsTUFBQSxDQUFDLENBQUMsT0FBTywwQ0FBRSxVQUFVLEtBQUksTUFBTSxDQUFDO1FBQzNDLElBQUksVUFBVSxHQUFHLE1BQU0sQ0FBQyxDQUFDLENBQUMsU0FBUyxDQUFDLENBQUMsQ0FBQyxNQUFBLENBQUMsQ0FBQyxPQUFPLDBDQUFFLFVBQVUsQ0FBQztRQUU1RCxJQUFJLEVBQUUsQ0FBQyxJQUFJLENBQUMsYUFBYSxJQUFJLENBQUMsQ0FBQSxNQUFBLENBQUMsQ0FBQyxPQUFPLDBDQUFFLFNBQVMsQ0FBQSxFQUFFO1lBQ2xELE1BQU0sQ0FBQyw4QkFBOEIsQ0FBQyxDQUFDLFFBQVEsUUFBUSxDQUFDLENBQUMsSUFBSSx1Q0FBdUMsQ0FBQyxDQUFDO1lBQ3RHLE9BQU8sU0FBUyxDQUFDO1NBQ2xCO1FBRUQsSUFBSSxJQUFJLElBQUksQ0FBQyxFQUFFLENBQUMsSUFBSSxDQUFDLGFBQWE7WUFDaEMsTUFBTSxDQUFDLDhCQUE4QixDQUFDLENBQUMsUUFBUSxRQUFRLENBQUMsQ0FBQyxJQUFJLElBQUksQ0FBQyxDQUFDO1FBQ3JFLElBQUksQ0FBQyxJQUFJO1lBQ1AsTUFBTSxDQUFDLDhCQUE4QixDQUFDLENBQUMsUUFBUSxRQUFRLENBQUMsQ0FBQyxJQUFJLElBQUksQ0FBQyxDQUFDO1FBQ3JFLE1BQU0sS0FBSyxHQUFHLElBQUksQ0FBQyxHQUFHLEVBQUUsQ0FBQztRQUN6QixNQUFNLFNBQVMsR0FBRyxJQUFJLElBQUksQ0FBQyxLQUFLLENBQUMsQ0FBQyxXQUFXLEVBQUUsQ0FBQztRQUNoRCxJQUFJO1lBQ0YsSUFBSSxJQUFJO2dCQUNOLENBQUMsR0FBRyxFQUFFLElBQUksRUFBRSxDQUFDLENBQUMsSUFBSSxFQUFFLEtBQUssRUFBQyxNQUFBLE1BQUEsQ0FBQyxDQUFDLE9BQU8sMENBQUUsS0FBSyxtQ0FBSSxFQUFFLEVBQUUsUUFBUSxFQUFFLENBQUMsQ0FBQyxRQUFRLEVBQUUsSUFBSSxFQUFFLEVBQUUsRUFBRSxJQUFJLEVBQUUsU0FBUyxFQUFFLE9BQU8sRUFBRSxJQUFJLEVBQUUsTUFBTSxFQUFFLFVBQVcsRUFBRSxFQUFFLEVBQUUsQ0FBQyxFQUFFLE9BQU8sRUFBRSxJQUFJLEVBQUUsT0FBTyxFQUFFLFdBQVcsYUFBWCxXQUFXLGNBQVgsV0FBVyxHQUFJLEVBQUUsRUFBRSxPQUFPLEVBQUUsRUFBRSxDQUFDLElBQUksQ0FBQyxhQUFhLEVBQUMsQ0FBQztpQkFDdE47Z0JBQ0gsSUFBSSxRQUFRLEdBQUcsV0FBVyxhQUFYLFdBQVcsY0FBWCxXQUFXLEdBQUksZ0JBQWdCLENBQUM7Z0JBRS9DLElBQUksRUFBRSxDQUFDLElBQUksQ0FBQyxXQUFXO29CQUNyQixPQUFPLENBQUMsT0FBTyxDQUFDLEdBQUcsQ0FBQyxDQUFDLFFBQVEsS0FBSyxDQUFDLENBQUMsSUFBSSxFQUFFLENBQUMsQ0FBQztnQkFFOUMsQ0FBQyxHQUFHLEVBQUUsSUFBSSxFQUFFLENBQUMsQ0FBQyxJQUFJLEVBQUUsS0FBSyxFQUFDLE1BQUEsTUFBQSxDQUFDLENBQUMsT0FBTywwQ0FBRSxLQUFLLG1DQUFJLEVBQUUsRUFBRSxRQUFRLEVBQUUsQ0FBQyxDQUFDLFFBQVEsRUFBRSxJQUFJLEVBQUUsRUFBRSxFQUFFLElBQUksRUFBRSxTQUFTLEVBQUUsT0FBTyxFQUFFLElBQUksRUFBRSxNQUFNLEVBQUUsTUFBQSxDQUFDLE1BQU0sT0FBTyxDQUFDLENBQUMsQ0FBQyxJQUFJLEVBQUUsUUFBUSxDQUFDLENBQUMsQ0FBQyxRQUFRLEVBQUUsbUNBQUksSUFBSSxFQUFFLEVBQUUsRUFBRSxDQUFDLEVBQUUsT0FBTyxFQUFFLEtBQUssRUFBRyxPQUFPLEVBQUUsV0FBVyxhQUFYLFdBQVcsY0FBWCxXQUFXLEdBQUksRUFBRSxFQUFFLE9BQU8sRUFBRSxFQUFFLENBQUMsSUFBSSxDQUFDLGFBQWEsRUFBQyxDQUFDO2dCQUVwUSxJQUFJLEVBQUUsQ0FBQyxJQUFJLENBQUMsV0FBVyxFQUFFO29CQUN2QixPQUFPLENBQUMsVUFBVSxDQUFDLEdBQUcsQ0FBQyxDQUFDLFFBQVEsS0FBSyxDQUFDLENBQUMsSUFBSSxFQUFFLENBQUMsQ0FBQztvQkFDL0MsSUFBSSxDQUFDLEtBQUssQ0FBQyxJQUFJLENBQUMsZ0JBQWdCLENBQUMsQ0FBQyxRQUFRLEtBQUssQ0FBQyxDQUFDLElBQUksd0dBQXdHLENBQUMsQ0FBQztpQkFDaEs7YUFDRjtTQUNGO1FBQUMsT0FBTyxDQUFNLEVBQUU7WUFDZixRQUFRLENBQUMsQ0FBQyxDQUFDLENBQUM7WUFDWixDQUFDLEdBQUcsRUFBRSxJQUFJLEVBQUUsQ0FBQyxDQUFDLElBQUksRUFBRSxLQUFLLEVBQUMsTUFBQSxNQUFBLENBQUMsQ0FBQyxPQUFPLDBDQUFFLEtBQUssbUNBQUksRUFBRSxFQUFFLFFBQVEsRUFBRSxDQUFDLENBQUMsUUFBUSxFQUFFLElBQUksRUFBRSxFQUFFLEVBQUUsSUFBSSxFQUFFLFNBQVMsRUFBRSxPQUFPLEVBQUUsS0FBSyxFQUFFLE1BQU0sRUFBRSxNQUFNLFNBQVMsQ0FBQyxDQUFDLENBQUMsRUFBRSxFQUFFLEVBQUUsQ0FBQyxFQUFFLE9BQU8sRUFBRSxLQUFLLEVBQUUsT0FBTyxFQUFFLFdBQVcsYUFBWCxXQUFXLGNBQVgsV0FBVyxHQUFJLEVBQUUsRUFBRSxPQUFPLEVBQUUsS0FBSyxFQUFDLENBQUM7U0FDbk47UUFDRCxJQUFJLENBQUEsTUFBQSxDQUFDLENBQUMsT0FBTywwQ0FBRSxZQUFZLEtBQUksQ0FBQyxDQUFDLE1BQU0sQ0FBQyxXQUFXLEtBQUssRUFBRSxDQUFDLFNBQVMsRUFBRTtZQUNwRSxNQUFNLEdBQUcsR0FBRyxDQUFDLENBQUMsTUFBTSxDQUFDLEdBQUcsQ0FBQyxTQUFTLENBQUMsQ0FBQztZQUNwQyxJQUFJLEdBQUc7Z0JBQ0wsQ0FBQyxDQUFDLE9BQU8sR0FBRyxHQUFHLENBQUMsS0FBSyxDQUFDLEdBQUcsS0FBSyxHQUFHLENBQUMsTUFBTSxDQUFDO1lBQzNDLElBQUksQ0FBQyxPQUFPLEVBQUU7Z0JBQ1osTUFBTSxFQUFFLEdBQUcsQ0FBQyxDQUFDLE1BQU0sQ0FBQztnQkFDcEIsRUFBRSxDQUFDLE9BQU8sQ0FBQyxNQUFNLENBQUMsT0FBTyxDQUFDLENBQUM7Z0JBQzNCLEVBQUUsQ0FBQyxJQUFJLENBQUMsV0FBVyxDQUFDLENBQUMsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxDQUFDLENBQUMsR0FBRyxDQUFDLFNBQVMsQ0FBQyxDQUFDLENBQUM7Z0JBQzdDLENBQUMsQ0FBQyxNQUFNLEdBQUcsRUFBRSxDQUFDO2FBQ2Y7WUFDRCxDQUFDLENBQUMsTUFBTSxHQUFHLENBQUMsQ0FBQyxNQUFNLENBQUMsS0FBSyxFQUFFLENBQUM7U0FDN0I7UUFDRCxDQUFDLENBQUMsSUFBSSxHQUFHLElBQUksQ0FBQyxJQUFJLENBQUMsSUFBSSxDQUFDLENBQUM7UUFDekIsQ0FBQyxDQUFDLEVBQUUsR0FBRyxJQUFJLENBQUMsR0FBRyxFQUFFLEdBQUcsS0FBSyxDQUFDO1FBQzFCLElBQUksQ0FBQyxJQUFJO1lBQ1AsTUFBTSxDQUFDLCtCQUErQixDQUFDLENBQUMsUUFBUSxRQUFRLENBQUMsQ0FBQyxJQUFJLGFBQWEsQ0FBQyxDQUFDLE9BQU8sQ0FBQyxDQUFDLENBQUMsU0FBUyxDQUFDLENBQUMsQ0FBQyxPQUFPLFVBQVUsQ0FBQyxDQUFDLEVBQUUsS0FBSyxDQUFDLENBQUM7UUFDakksSUFBSSxDQUFDLENBQUMsQ0FBQyxPQUFPLEVBQUU7WUFDWixNQUFNLENBQUMsaUNBQWlDLENBQUMsQ0FBQyxRQUFRLFFBQVEsQ0FBQyxDQUFDLElBQUksT0FBTyxDQUFDLENBQUMsTUFBTSxFQUFFLENBQUMsQ0FBQztTQUN0RjtRQUNELENBQUMsQ0FBQyxRQUFRLEdBQUcsQ0FBQyxDQUFDLFFBQVEsQ0FBQztRQUN4QixDQUFDLENBQUMsSUFBSSxHQUFHLENBQUMsQ0FBQyxJQUFJLENBQUM7UUFDaEIsQ0FBQyxDQUFDLEtBQUssR0FBRyxNQUFBLE1BQUEsQ0FBQyxDQUFDLE9BQU8sMENBQUUsS0FBSyxtQ0FBSSxFQUFFLENBQUM7UUFDakMsSUFBSSxDQUFDLE1BQU0sRUFBRTtZQUNYLElBQUksTUFBTSxHQUFHO2dCQUNYLFNBQVMsRUFBRSxDQUFDLENBQUMsT0FBTyxFQUFFLFFBQVEsRUFBRSxDQUFDLENBQUMsTUFBTSxFQUFFLElBQUksRUFBRSxDQUFDLENBQUMsRUFBRSxFQUFFLE1BQU0sRUFBRSxDQUFDLENBQUMsSUFBSTtnQkFDcEUsU0FBUyxFQUFFLENBQUMsQ0FBQyxPQUFPLEVBQUUsVUFBVSxFQUFFLENBQUMsQ0FBQyxRQUFRLEVBQUUsTUFBTSxFQUFFLENBQUMsQ0FBQyxJQUFJLEVBQUUsTUFBTSxFQUFFLENBQUMsQ0FBQyxJQUFJLEVBQUUsT0FBTyxFQUFFLENBQUMsQ0FBQyxLQUFLO2dCQUM5RixTQUFTLEVBQUUsRUFBRSxDQUFDLElBQUksQ0FBQyxhQUFhLElBQUksQ0FBQyxDQUFDLE9BQU87Z0JBQzdDLFNBQVMsRUFBRSxDQUFDLENBQUMsT0FBTzthQUNyQixDQUFDO1lBQ0YsSUFBSSxDQUFDLENBQUMsTUFBTSxDQUFDLFdBQVcsSUFBSSxNQUFNLEVBQUU7Z0JBQ2xDLE1BQU0sR0FBRyxHQUFHLE1BQU0sQ0FBQyxJQUFJLENBQUMsQ0FBQyxDQUFDLE1BQU0sQ0FBQyxDQUFDLE1BQU0sQ0FBQyxDQUFDLEdBQUcsRUFBRSxDQUFDLEVBQUUsRUFBRSxDQUFDLGlDQUFNLEdBQUcsS0FBRSxDQUFDLFNBQVMsR0FBRyxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsTUFBTSxDQUFDLENBQUMsQ0FBQyxJQUFHLEVBQUUsRUFBRSxDQUFDLENBQUM7Z0JBQ3JHLE1BQU0sbUNBQVEsTUFBTSxHQUFLLEdBQUcsQ0FBRSxDQUFDO2FBQ2hDO1lBRUQsSUFBSSxNQUFNLENBQUMsTUFBTSxZQUFZLEVBQUUsQ0FBQyxTQUFTO2dCQUN2QyxNQUFNLENBQUMsTUFBTSxHQUFHLElBQUksQ0FBQyxTQUFTLENBQUMsTUFBQSxNQUFNLENBQUMsTUFBTSwwQ0FBRSxNQUFNLEVBQUUsQ0FBQyxJQUFJLEVBQUUsQ0FBQztZQUNoRSxNQUFNLElBQUksQ0FBQyxLQUFLLENBQUMsVUFBVSxDQUFDLElBQUksRUFBRSxNQUFNLENBQUMsQ0FBQztTQUMzQztRQUNELE9BQU8sQ0FBQyxDQUFDOztDQUNWO0FBRUQsTUFBTSxVQUFVLE9BQU8sQ0FBQyxLQUFZO0lBQ2xDLE1BQU0sTUFBTSxHQUFHLEtBQUssQ0FBQyxLQUFLLEVBQUUsQ0FBQztJQUM3QixNQUFNLENBQUMsSUFBSSxDQUFDLEdBQUcsRUFBRSxDQUFDLElBQUksQ0FBQyxNQUFNLEVBQUUsR0FBRyxHQUFHLENBQUMsQ0FBQztJQUN2QyxPQUFPLE1BQU0sQ0FBQztBQUNoQixDQUFDO0FBRUQsNkJBQTZCO0FBQzdCLE1BQU0sVUFBZ0IsS0FBSyxDQUFDLEVBQVU7O1FBQ3BDLE1BQU0sSUFBSSxPQUFPLENBQUMsQ0FBQyxDQUFDLEVBQUUsRUFBRSxDQUFDLFVBQVUsQ0FBQyxDQUFDLEVBQUUsRUFBRSxDQUFDLENBQUMsQ0FBQztJQUM5QyxDQUFDO0NBQUE7QUFFRCxNQUFNLFVBQWdCLFVBQVUsQ0FBQyxZQUEyQixFQUMxRCxRQUFnQixrQkFBa0IsRUFBRSxPQUFlLEdBQUcsRUFBRSxXQUFtQixFQUFFOztRQUM3RSxPQUFPLElBQUksT0FBTyxDQUFDLENBQUMsT0FBTyxFQUFFLE1BQU0sRUFBRSxFQUFFO1lBQ3JDLFVBQVUsQ0FBQyxHQUFHLEVBQUU7Z0JBQ2QsYUFBYSxDQUFDLFVBQVUsQ0FBQyxDQUFDO2dCQUMxQixNQUFNLENBQUMsSUFBSSxLQUFLLENBQUMsS0FBSyxDQUFDLENBQUMsQ0FBQztZQUMzQixDQUFDLEVBQUUsSUFBSSxDQUFDLENBQUM7WUFDVCxhQUFhO1lBQ2IsTUFBTSxVQUFVLEdBQVksV0FBVyxDQUFDLEdBQUcsRUFBRTtnQkFDM0MsSUFBSSxZQUFZLEVBQUUsRUFBRTtvQkFDbEIsYUFBYSxDQUFDLFVBQVUsQ0FBQyxDQUFDO29CQUMxQixPQUFPLENBQUMsSUFBSSxDQUFDLENBQUM7aUJBQ2Y7WUFDSCxDQUFDLEVBQUUsUUFBUSxDQUFDLENBQUM7UUFDZixDQUFDLENBQUMsQ0FBQztJQUNMLENBQUM7Q0FBQTtBQUVELCtEQUErRDtBQUMvRCxNQUFNLFVBQWdCLE9BQU8sQ0FBQyxJQUF3QixFQUFFLFdBQW1CLEVBQUUsZ0JBQXdCLG1CQUFtQjs7UUFDdEgsSUFBSSxPQUFPLEdBQVEsSUFBSSxDQUFDO1FBQ3hCLE1BQU0sY0FBYyxHQUFHLElBQUksT0FBTyxDQUFNLENBQUMsQ0FBQyxFQUFFLE1BQU0sRUFBRSxFQUFFO1lBQ3BELE9BQU8sR0FBRyxVQUFVLENBQUMsR0FBRyxFQUFFO2dCQUN4Qix3REFBd0Q7Z0JBQ3hELE1BQU0sQ0FBQyxhQUFhLENBQUMsQ0FBQztZQUN4QixDQUFDLEVBQUUsV0FBVyxDQUFDLENBQUM7UUFDbEIsQ0FBQyxDQUFDLENBQUM7UUFDSCxJQUFJO1lBQ0YsT0FBTyxNQUFNLE9BQU8sQ0FBQyxJQUFJLENBQUMsQ0FBQyxJQUFJLEVBQUUsRUFBRSxjQUFjLENBQUMsQ0FBQyxDQUFDO1NBQ3JEO2dCQUFTO1lBQ1IsSUFBSSxPQUFPO2dCQUNULFlBQVksQ0FBQyxPQUFPLENBQUMsQ0FBQztTQUN6QjtJQUNILENBQUM7Q0FBQTtBQUVELE1BQU0sVUFBVSxlQUFlLENBQUMsV0FBbUI7SUFDakQsTUFBTSxPQUFPLEdBQUcsRUFBRSxDQUFDLE1BQU0sQ0FBQyxjQUFjLEVBQUUsQ0FBQztJQUMzQyxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsT0FBTyxDQUFDLE1BQU0sRUFBRSxDQUFDLEVBQUUsRUFBRTtRQUN2QyxJQUFJLE9BQU8sQ0FBQyxDQUFDLENBQUMsQ0FBQyxLQUFLLElBQUksV0FBVztZQUNqQyxPQUFPLElBQUksQ0FBQztLQUNmO0lBQ0QsT0FBTyxLQUFLLENBQUM7QUFDZixDQUFDO0FBRUQ7Ozs7O0dBS0c7QUFDSCxNQUFNLFVBQWdCLG9CQUFvQixDQUFDLE1BQTJCLEVBQ3BFLEtBQW1DOztRQUNuQyxJQUFJLE1BQU0sR0FBWSxLQUFLLENBQUM7UUFDNUIsSUFBSSxPQUFPLEdBQVksS0FBSyxDQUFDO1FBQzdCLElBQUk7WUFDRixNQUFNLE1BQU0sRUFBRSxDQUFDO1NBQ2hCO1FBQUMsT0FBTyxDQUFDLEVBQUU7WUFDVixNQUFNLEdBQUcsSUFBSSxDQUFDO1lBQ2QsT0FBTyxHQUFHLENBQUMsS0FBSyxJQUFJLEtBQUssQ0FBQyxDQUFDLENBQUMsQ0FBQztTQUM5QjtnQkFBUztZQUNSLElBQUksQ0FBQyxNQUFNO2dCQUNULE1BQU0sSUFBSSxLQUFLLENBQUMseUNBQXlDLENBQUMsQ0FBQztZQUM3RCxJQUFJLENBQUMsT0FBTztnQkFDVixNQUFNLElBQUksS0FBSyxDQUFDLHdFQUF3RSxDQUFDLENBQUM7U0FDN0Y7SUFDSCxDQUFDO0NBQUE7QUFFRCxNQUFNLEtBQUssR0FBRyxFQUFFLENBQUMsU0FBUyxDQUFDLFdBQVcsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxNQUFNLENBQUMsV0FBVyxDQUFDLEtBQUssRUFBRSxDQUFDLE1BQU0sRUFBRSxNQUFNLEVBQUUsTUFBTSxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFFakc7Ozs7Ozs7Ozs7R0FVRztBQUNILE1BQU0sVUFBZ0IsVUFBVSxDQUFDLENBQVMsRUFBRSxFQUFpQixFQUFFLE9BRzlEOzs7UUFDQyxNQUFNLFdBQVcsR0FBRyxNQUFBLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxXQUFXLG1DQUFJLEVBQUUsQ0FBQztRQUMvQyxJQUFJLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxtQkFBbUI7WUFDOUIsTUFBTSxJQUFJLENBQUMsSUFBSSxDQUFDLG1CQUFtQixDQUFDLEVBQUUsQ0FBQyxDQUFDO1FBQzFDLE1BQU0sRUFBRSxHQUFHLElBQUksQ0FBQyxLQUFLLENBQUMsWUFBWSxDQUFDLEVBQUUsQ0FBQyxDQUFDO1FBRXZDLElBQUk7WUFDRiwrQkFBK0I7WUFDL0IsTUFBTSxrQkFBa0IsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFLFdBQVcsRUFBRSxJQUFJLENBQUMsTUFBTSxDQUFDLGFBQWEsQ0FBQyxDQUFDO1lBQ3hFLG1FQUFtRTtZQUNuRSxJQUFJLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxXQUFXO2dCQUN0QixNQUFNLGtCQUFrQixDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUUsV0FBVyxFQUFFLElBQUksQ0FBQyxNQUFNLENBQUMsYUFBYSxFQUFFLFNBQVMsRUFBRSxPQUFRLENBQUMsV0FBVyxDQUFDLENBQUM7WUFFM0csdURBQXVEO1lBQ3ZELElBQUksQ0FBQyxDQUFBLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxRQUFRLENBQUEsRUFBRTtnQkFDdEIsTUFBTSxrQkFBa0IsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFLFdBQVcsRUFBRSxJQUFJLENBQUMsTUFBTSxDQUFDLGFBQWEsRUFBRSx5QkFBeUIsQ0FBQyxDQUFDO2dCQUNuRyxJQUFJLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxXQUFXO29CQUN0QixNQUFNLGtCQUFrQixDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUUsV0FBVyxFQUFFLElBQUksQ0FBQyxNQUFNLENBQUMsYUFBYSxFQUFFLHlCQUF5QixFQUFFLE9BQVEsQ0FBQyxXQUFXLENBQUMsQ0FBQzthQUM1SDtZQUVELHVEQUF1RDtZQUN2RCxJQUFJLGNBQWMsR0FBNEMsSUFBSSxDQUFDO1lBQ25FLGNBQWMsR0FBRyxNQUFNLGtCQUFrQixDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUUsV0FBVyxFQUFFLElBQUksQ0FBQyxNQUFNLENBQUMsYUFBYSxFQUFFLHVCQUF1QixDQUFDLENBQUM7WUFDbEgsSUFBSSxPQUFPLGFBQVAsT0FBTyx1QkFBUCxPQUFPLENBQUUsV0FBVztnQkFDdEIsY0FBYyxHQUFHLE1BQU0sa0JBQWtCLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRSxXQUFXLEVBQUUsSUFBSSxDQUFDLE1BQU0sQ0FBQyxhQUFhLEVBQ3JGLHVCQUF1QixFQUFFLE9BQVEsQ0FBQyxXQUFXLENBQUMsQ0FBQTtZQUVsRCxnQkFBZ0I7WUFDaEIsTUFBTSxrQkFBa0IsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFLFdBQVcsRUFBRSxJQUFJLENBQUMsTUFBTSxDQUFDLG1CQUFtQixFQUFFLFVBQVUsRUFBRSxTQUFTLEVBQUUsY0FBYyxhQUFkLGNBQWMsdUJBQWQsY0FBYyxDQUFFLE1BQU0sRUFDekgsRUFBRSxVQUFVLEVBQUUsY0FBYyxhQUFkLGNBQWMsdUJBQWQsY0FBYyxDQUFFLFVBQVUsRUFBRSxDQUFDLENBQUM7WUFDOUMsSUFBSSxPQUFPLGFBQVAsT0FBTyx1QkFBUCxPQUFPLENBQUUsV0FBVztnQkFDdEIsTUFBTSxrQkFBa0IsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFLFdBQVcsRUFBRSxJQUFJLENBQUMsTUFBTSxDQUFDLG1CQUFtQixFQUFFLFVBQVUsRUFBRSxPQUFRLENBQUMsV0FBVyxFQUM1RyxjQUFjLGFBQWQsY0FBYyx1QkFBZCxjQUFjLENBQUUsTUFBTSxFQUFFLEVBQUUsVUFBVSxFQUFFLGNBQWMsYUFBZCxjQUFjLHVCQUFkLGNBQWMsQ0FBRSxVQUFVLEVBQUUsQ0FBQyxDQUFDO1lBRXhFLG9DQUFvQztZQUNwQyxJQUFJLENBQUEsT0FBTyxhQUFQLE9BQU8sdUJBQVAsT0FBTyxDQUFFLGVBQWUsTUFBSyxLQUFLLEVBQUU7Z0JBQ3RDLEVBQUUsQ0FBQyxTQUFTLEdBQUcsS0FBSyxDQUFDO2dCQUNyQixNQUFNLEtBQUssQ0FBQyxFQUFFLENBQUMsQ0FBQztnQkFDaEIsTUFBTSxrQkFBa0IsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFLFdBQVcsRUFBRSxJQUFJLENBQUMsTUFBTSxDQUFDLGFBQWEsQ0FBQyxDQUFDO2dCQUN4RSxJQUFJLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxXQUFXO29CQUN0QixNQUFNLGtCQUFrQixDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUUsV0FBVyxFQUFFLElBQUksQ0FBQyxNQUFNLENBQUMsYUFBYSxFQUFFLFNBQVMsRUFBRSxPQUFRLENBQUMsV0FBVyxDQUFDLENBQUM7YUFDNUc7WUFFRCw2QkFBNkI7WUFDN0IsTUFBTSxrQkFBa0IsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFLFdBQVcsRUFBRSxJQUFJLENBQUMsTUFBTSxDQUFDLGFBQWEsRUFBRSxXQUFXLENBQUMsQ0FBQztZQUNyRixJQUFJLE9BQU8sYUFBUCxPQUFPLHVCQUFQLE9BQU8sQ0FBRSxXQUFXO2dCQUN0QixNQUFNLGtCQUFrQixDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUUsV0FBVyxFQUFFLElBQUksQ0FBQyxNQUFNLENBQUMsYUFBYSxFQUFFLFdBQVcsRUFBRSxPQUFRLENBQUMsV0FBVyxDQUFDLENBQUM7U0FFOUc7Z0JBQVM7WUFDUixpREFBaUQ7WUFDakQseUJBQXlCO1lBQ3pCLHlCQUF5QjtTQUMxQjs7Q0FDRiIsInNvdXJjZXNDb250ZW50IjpbImltcG9ydCB0eXBlICogYXMgX2dyb2sgZnJvbSAnZGF0YWdyb2stYXBpL2dyb2snO1xyXG5pbXBvcnQgdHlwZSAqIGFzIF9ERyBmcm9tICdkYXRhZ3Jvay1hcGkvZGcnO1xyXG5kZWNsYXJlIGxldCBncm9rOiB0eXBlb2YgX2dyb2ssIERHOiB0eXBlb2YgX0RHO1xyXG5cclxuaW1wb3J0IHsgT2JzZXJ2YWJsZSB9IGZyb20gJ3J4anMnO1xyXG5pbXBvcnQgeyB0ZXN0RGF0YSB9IGZyb20gJy4vZGF0YWZyYW1lLXV0aWxzJztcclxuaW1wb3J0IFRpbWVvdXQgPSBOb2RlSlMuVGltZW91dDtcclxuaW1wb3J0IHsgY2hhbmdlT3B0aW9uc1NhdmVMYXlvdXQsIGZpbHRlckFzeW5jLCBsb2FkTGF5b3V0LCBzZWxlY3RGaWx0ZXJDaGFuZ2VDdXJyZW50LCB0ZXN0Vmlld2VySW50ZXJuYWwgfSBmcm9tICcuL3Rlc3Qtdmlld2VyLXV0aWxzJztcclxuXHJcbmNvbnN0IFNUQU5EQVJUX1RJTUVPVVQgPSAzMDAwMDtcclxuY29uc3QgQkVOQ0hNQVJLX1RJTUVPVVQgPSAxMDgwMDAwMDtcclxuXHJcbmNvbnN0IHN0ZExvZyA9IGNvbnNvbGUubG9nLmJpbmQoY29uc29sZSk7XHJcbmNvbnN0IHN0ZEluZm8gPSBjb25zb2xlLmluZm8uYmluZChjb25zb2xlKTtcclxuY29uc3Qgc3RkV2FybiA9IGNvbnNvbGUud2Fybi5iaW5kKGNvbnNvbGUpO1xyXG5jb25zdCBzdGRFcnJvciA9IGNvbnNvbGUuZXJyb3IuYmluZChjb25zb2xlKTtcclxuXHJcbmV4cG9ydCBjb25zdCB0ZXN0czoge1xyXG4gIFtrZXk6IHN0cmluZ106IENhdGVnb3J5XHJcbn0gPSB7fTtcclxuXHJcbmNvbnN0IGF1dG9UZXN0c0NhdE5hbWUgPSAnQXV0byBUZXN0cyc7XHJcbmNvbnN0IGRlbW9DYXROYW1lID0gJ0RlbW8nO1xyXG5jb25zdCBkZXRlY3RvcnNDYXROYW1lID0gJ0RldGVjdG9ycyc7XHJcbmNvbnN0IGNvcmVDYXROYW1lID0gJ0NvcmUnO1xyXG5jb25zdCB3YXNSZWdpc3RlcmVkOiB7IFtrZXk6IHN0cmluZ106IGJvb2xlYW4gfSA9IHt9O1xyXG5leHBvcnQgbGV0IGN1cnJlbnRDYXRlZ29yeTogc3RyaW5nO1xyXG5cclxuZXhwb3J0IG5hbWVzcGFjZSBhc3N1cmUge1xyXG4gIGV4cG9ydCBmdW5jdGlvbiBub3ROdWxsKHZhbHVlOiBhbnksIG5hbWU/OiBzdHJpbmcpIHtcclxuICAgIGlmICh2YWx1ZSA9PSBudWxsKVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYCR7bmFtZSA9PSBudWxsID8gJ1ZhbHVlJyA6IG5hbWV9IG5vdCBkZWZpbmVkYCk7XHJcbiAgfVxyXG59XHJcblxyXG5leHBvcnQgaW50ZXJmYWNlIFRlc3RPcHRpb25zIHtcclxuICB0aW1lb3V0PzogbnVtYmVyO1xyXG4gIGJlbmNobWFya1dhcm5UaW1lb3V0PzogbnVtYmVyO1xyXG4gIGJlbmNobWFya1RpbWVvdXQ/OiBudW1iZXI7XHJcbiAgdW5oYW5kbGVkRXhjZXB0aW9uVGltZW91dD86IG51bWJlcjtcclxuICBza2lwUmVhc29uPzogc3RyaW5nO1xyXG4gIGlzQWdncmVnYXRlZD86IGJvb2xlYW47XHJcbiAgYmVuY2htYXJrPzogYm9vbGVhbjtcclxuICBzdHJlc3NUZXN0PzogYm9vbGVhbjtcclxuICBvd25lcj86IHN0cmluZztcclxuICB0YWdzPzogc3RyaW5nW107XHJcbn1cclxuXHJcbmV4cG9ydCBpbnRlcmZhY2UgVGVzdFJlc3VsdCB7XHJcbiAgZGF0ZTogc3RyaW5nO1xyXG4gIGNhdGVnb3J5OiBzdHJpbmc7XHJcbiAgbmFtZTogc3RyaW5nO1xyXG4gIHN1Y2Nlc3M6IGJvb2xlYW47XHJcbiAgcmVzdWx0OiBhbnk7XHJcbiAgbXM6IG51bWJlcjtcclxuICBza2lwcGVkOiBib29sZWFuO1xyXG4gIGxvZ3M6IHN0cmluZztcclxuICBvd25lcjogc3RyaW5nO1xyXG4gIHBhY2thZ2U6IHN0cmluZztcclxuICBmbGFraW5nOiBib29sZWFuO1xyXG59XHJcblxyXG5cclxuZXhwb3J0IGludGVyZmFjZSBUZXN0UmVzdWx0RXh0ZW5kZWQgZXh0ZW5kcyBUZXN0UmVzdWx0e1xyXG4gIHdpZGdldHNEaWZmZXJlbmNlOiBudW1iZXI7XHJcbn1cclxuXHJcbmV4cG9ydCBpbnRlcmZhY2UgQ2F0ZWdvcnlPcHRpb25zIHtcclxuICBjbGVhcj86IGJvb2xlYW47XHJcbiAgdGltZW91dD86IG51bWJlcjtcclxuICBiZW5jaG1hcmtzPzogYm9vbGVhbjtcclxuICBzdHJlc3NUZXN0cz86IGJvb2xlYW47XHJcbiAgb3duZXI/OiBzdHJpbmc7XHJcbn1cclxuXHJcbmV4cG9ydCBjbGFzcyBUZXN0Q29udGV4dCB7XHJcbiAgc3RyZXNzVGVzdD86IGJvb2xlYW47XHJcbiAgY2F0Y2hVbmhhbmRsZWQgPSB0cnVlO1xyXG4gIHJlcG9ydCA9IGZhbHNlO1xyXG4gIHJldHVybk9uRmFpbCA9IGZhbHNlO1xyXG5cclxuICBjb25zdHJ1Y3RvcihjYXRjaFVuaGFuZGxlZD86IGJvb2xlYW4sIHJlcG9ydD86IGJvb2xlYW4sIHJldHVybk9uRmFpbD86IGJvb2xlYW4pIHtcclxuICAgIGlmIChjYXRjaFVuaGFuZGxlZCAhPT0gdW5kZWZpbmVkKSB0aGlzLmNhdGNoVW5oYW5kbGVkID0gY2F0Y2hVbmhhbmRsZWQ7XHJcbiAgICBpZiAocmVwb3J0ICE9PSB1bmRlZmluZWQpIHRoaXMucmVwb3J0ID0gcmVwb3J0O1xyXG4gICAgaWYgKHJldHVybk9uRmFpbCAhPT0gdW5kZWZpbmVkKSB0aGlzLnJldHVybk9uRmFpbCA9IHJldHVybk9uRmFpbDtcclxuICB9O1xyXG59XHJcblxyXG5leHBvcnQgY2xhc3MgVGVzdCB7XHJcbiAgdGVzdDogKCkgPT4gUHJvbWlzZTxhbnk+O1xyXG4gIG5hbWU6IHN0cmluZztcclxuICBjYXRlZ29yeTogc3RyaW5nO1xyXG4gIG9wdGlvbnM/OiBUZXN0T3B0aW9ucztcclxuXHJcbiAgY29uc3RydWN0b3IoY2F0ZWdvcnk6IHN0cmluZywgbmFtZTogc3RyaW5nLCB0ZXN0OiAoKSA9PiBQcm9taXNlPGFueT4sIG9wdGlvbnM/OiBUZXN0T3B0aW9ucykge1xyXG4gICAgdGhpcy5jYXRlZ29yeSA9IGNhdGVnb3J5O1xyXG4gICAgdGhpcy5uYW1lID0gbmFtZTtcclxuICAgIG9wdGlvbnMgPz89IHt9O1xyXG4gICAgb3B0aW9ucy50aW1lb3V0ID8/PSBTVEFOREFSVF9USU1FT1VUO1xyXG4gICAgdGhpcy5vcHRpb25zID0gb3B0aW9ucztcclxuICAgIHRoaXMudGVzdCA9IGFzeW5jICgpOiBQcm9taXNlPGFueT4gPT4ge1xyXG4gICAgICByZXR1cm4gbmV3IFByb21pc2UoYXN5bmMgKHJlc29sdmUsIHJlamVjdCkgPT4ge1xyXG4gICAgICAgIGxldCByZXN1bHQgPSAnJztcclxuICAgICAgICB0cnkge1xyXG4gICAgICAgICAgaWYgKERHLlRlc3QuaXNJbkRlYnVnKVxyXG4gICAgICAgICAgICBkZWJ1Z2dlcjtcclxuXHJcbiAgICAgICAgICBsZXQgcmVzID0gYXdhaXQgdGVzdCgpO1xyXG4gICAgICAgICAgdHJ5IHtcclxuICAgICAgICAgICAgcmVzdWx0ID0gcmVzPy50b1N0cmluZygpID8/ICcnO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgICAgY2F0Y2ggKGUpIHtcclxuICAgICAgICAgICAgcmVzdWx0ID0gJ0NhblxcJ3QgY29udmVydCB0ZXN0XFwncyByZXN1bHQgdG8gc3RyaW5nJztcclxuICAgICAgICAgICAgY29uc29sZS5lcnJvcihgQ2FuXFwndCBjb252ZXJ0IHRlc3RcXCdzIHJlc3VsdCB0byBzdHJpbmcgaW4gdGhlICR7dGhpcy5jYXRlZ29yeX06JHt0aGlzLm5hbWV9IHRlc3RgKTtcclxuICAgICAgICAgIH1cclxuICAgICAgICB9IGNhdGNoIChlOiBhbnkpIHtcclxuICAgICAgICAgIHJlamVjdChlKTtcclxuICAgICAgICB9XHJcbiAgICAgICAgcmVzb2x2ZShyZXN1bHQpO1xyXG4gICAgICB9KTtcclxuICAgIH07XHJcbiAgfVxyXG59XHJcblxyXG5leHBvcnQgY2xhc3MgQ2F0ZWdvcnkge1xyXG4gIHRlc3RzPzogVGVzdFtdO1xyXG4gIGJlZm9yZT86ICgpID0+IFByb21pc2U8dm9pZD47XHJcbiAgYWZ0ZXI/OiAoKSA9PiBQcm9taXNlPHZvaWQ+O1xyXG5cclxuICBiZWZvcmVTdGF0dXM/OiBzdHJpbmc7XHJcbiAgYWZ0ZXJTdGF0dXM/OiBzdHJpbmc7XHJcbiAgY2xlYXI/OiBib29sZWFuO1xyXG4gIHRpbWVvdXQ/OiBudW1iZXI7XHJcbiAgYmVuY2htYXJrcz86IGJvb2xlYW47XHJcbiAgYmVuY2htYXJrVGltZW91dD86IG51bWJlcjtcclxuICBzdHJlc3NUZXN0cz86IGJvb2xlYW47XHJcbiAgb3duZXI/OiBzdHJpbmc7XHJcbn1cclxuXHJcbmV4cG9ydCBjbGFzcyBOb2RlVGVzdEV4ZWN1dGlvbk9wdGlvbnMge1xyXG4gIHBhY2thZ2UhOiBfREcuUGFja2FnZTtcclxufVxyXG5cclxuZXhwb3J0IGNsYXNzIFRlc3RFeGVjdXRpb25PcHRpb25zIHtcclxuICBjYXRlZ29yeT86IHN0cmluZztcclxuICB0ZXN0Pzogc3RyaW5nO1xyXG4gIHRlc3RDb250ZXh0PzogVGVzdENvbnRleHQ7XHJcbiAgZXhjbHVkZT86IHN0cmluZ1tdO1xyXG4gIHZlcmJvc2U/OiBib29sZWFuO1xyXG4gIHN0cmVzc1Rlc3Q/OiBib29sZWFuO1xyXG4gIHRhZ3M/OiBzdHJpbmdbXTtcclxuICBub2RlT3B0aW9ucz86IE5vZGVUZXN0RXhlY3V0aW9uT3B0aW9ucztcclxuICBza2lwVG9DYXRlZ29yeT86IHN0cmluZztcclxuICBza2lwVG9UZXN0Pzogc3RyaW5nO1xyXG4gIHJldHVybk9uRmFpbD86IGJvb2xlYW47XHJcbn1cclxuXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiB0ZXN0RXZlbnQ8VD4oZXZlbnQ6IE9ic2VydmFibGU8VD4sXHJcbiAgaGFuZGxlcjogKGFyZ3M6IFQpID0+IHZvaWQsIHRyaWdnZXI6ICgpID0+IHZvaWQsIG1zOiBudW1iZXIgPSAwLCByZWFzb246IHN0cmluZyA9IGB0aW1lb3V0YFxyXG4pOiBQcm9taXNlPGFueT4ge1xyXG4gIHJldHVybiBuZXcgUHJvbWlzZSgocmVzb2x2ZSwgcmVqZWN0KSA9PiB7XHJcbiAgICBjb25zdCBzdWIgPSBldmVudC5zdWJzY3JpYmUoKGFyZ3M6IFQpID0+IHtcclxuICAgICAgdHJ5IHtcclxuICAgICAgICBoYW5kbGVyKGFyZ3MpO1xyXG4gICAgICAgIHJlc29sdmUoJ09LJyk7XHJcbiAgICAgIH0gY2F0Y2ggKGUpIHtcclxuICAgICAgICByZWplY3QoZSk7XHJcbiAgICAgIH0gZmluYWxseSB7XHJcbiAgICAgICAgc3ViLnVuc3Vic2NyaWJlKCk7XHJcbiAgICAgICAgY2xlYXJUaW1lb3V0KHRpbWVvdXQpO1xyXG4gICAgICB9XHJcbiAgICB9KTtcclxuICAgIGNvbnN0IHRpbWVvdXQgPSBzZXRUaW1lb3V0KCgpID0+IHtcclxuICAgICAgc3ViLnVuc3Vic2NyaWJlKCk7XHJcbiAgICAgIC8vIGVzbGludC1kaXNhYmxlLW5leHQtbGluZSBwcmVmZXItcHJvbWlzZS1yZWplY3QtZXJyb3JzXHJcbiAgICAgIHJlamVjdChyZWFzb24pO1xyXG4gICAgfSwgbXMpO1xyXG4gICAgdHJpZ2dlcigpO1xyXG4gIH0pO1xyXG59XHJcblxyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gdGVzdEV2ZW50QXN5bmM8VD4oZXZlbnQ6IE9ic2VydmFibGU8VD4sXHJcbiAgaGFuZGxlcjogKGFyZ3M6IFQpID0+IFByb21pc2U8dm9pZD4sIHRyaWdnZXI6ICgpID0+IHZvaWQsIG1zOiBudW1iZXIgPSAwLCByZWFzb246IHN0cmluZyA9IGB0aW1lb3V0YFxyXG4pOiBQcm9taXNlPGFueT4ge1xyXG4gIHJldHVybiBuZXcgUHJvbWlzZSgocmVzb2x2ZSwgcmVqZWN0KSA9PiB7XHJcbiAgICBjb25zdCBzdWIgPSBldmVudC5zdWJzY3JpYmUoKGFyZ3M6IFQpID0+IHtcclxuICAgICAgaGFuZGxlcihhcmdzKS50aGVuKCgpID0+IHtcclxuICAgICAgICByZXNvbHZlKCdPSycpO1xyXG4gICAgICB9KS5jYXRjaCgoZSkgPT4ge1xyXG4gICAgICAgIHJlamVjdChlKTtcclxuICAgICAgfSkuZmluYWxseSgoKSA9PiB7XHJcbiAgICAgICAgc3ViLnVuc3Vic2NyaWJlKCk7XHJcbiAgICAgICAgY2xlYXJUaW1lb3V0KHRpbWVvdXQpO1xyXG4gICAgICB9KTtcclxuICAgIH0pO1xyXG4gICAgY29uc3QgdGltZW91dCA9IHNldFRpbWVvdXQoKCkgPT4ge1xyXG4gICAgICBzdWIudW5zdWJzY3JpYmUoKTtcclxuICAgICAgLy8gZXNsaW50LWRpc2FibGUtbmV4dC1saW5lIHByZWZlci1wcm9taXNlLXJlamVjdC1lcnJvcnNcclxuICAgICAgcmVqZWN0KHJlYXNvbik7XHJcbiAgICB9LCBtcyk7XHJcbiAgICB0cmlnZ2VyKCk7XHJcbiAgfSk7XHJcbn1cclxuXHJcbmV4cG9ydCBmdW5jdGlvbiB0ZXN0KG5hbWU6IHN0cmluZywgdGVzdDogKCkgPT4gUHJvbWlzZTxhbnk+LCBvcHRpb25zPzogVGVzdE9wdGlvbnMpOiB2b2lkIHtcclxuICBpZiAodGVzdHNbY3VycmVudENhdGVnb3J5XSA9PSB1bmRlZmluZWQpXHJcbiAgICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldID0ge307XHJcbiAgaWYgKHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0udGVzdHMgPT0gdW5kZWZpbmVkKVxyXG4gICAgdGVzdHNbY3VycmVudENhdGVnb3J5XS50ZXN0cyA9IFtdO1xyXG4gIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0udGVzdHMhLnB1c2gobmV3IFRlc3QoY3VycmVudENhdGVnb3J5LCBuYW1lLCB0ZXN0LCBvcHRpb25zKSk7XHJcbn1cclxuXHJcbi8qIFRlc3RzIHR3byBvYmplY3RzIGZvciBlcXVhbGl0eSwgdGhyb3dzIGFuIGV4Y2VwdGlvbiBpZiB0aGV5IGFyZSBub3QgZXF1YWwuICovXHJcbmV4cG9ydCBmdW5jdGlvbiBleHBlY3QoYWN0dWFsOiBhbnksIGV4cGVjdGVkOiBhbnkgPSB0cnVlLCBlcnJvcj86IHN0cmluZyk6IHZvaWQge1xyXG4gIGlmIChlcnJvcilcclxuICAgIGVycm9yID0gYCR7ZXJyb3J9LCBgO1xyXG4gIGVsc2UgZXJyb3IgPSAnJztcclxuICBpZiAoYWN0dWFsICE9PSBleHBlY3RlZClcclxuICAgIHRocm93IG5ldyBFcnJvcihgJHtlcnJvcn1FeHBlY3RlZCBcIiR7ZXhwZWN0ZWR9XCIsIGdvdCBcIiR7YWN0dWFsfVwiYCk7XHJcbn1cclxuXHJcbmV4cG9ydCBmdW5jdGlvbiBleHBlY3RGbG9hdChhY3R1YWw6IG51bWJlciwgZXhwZWN0ZWQ6IG51bWJlciwgdG9sZXJhbmNlID0gMC4wMDEsIGVycm9yPzogc3RyaW5nKTogdm9pZCB7XHJcbiAgaWYgKChhY3R1YWwgPT09IE51bWJlci5QT1NJVElWRV9JTkZJTklUWSAmJiBleHBlY3RlZCA9PT0gTnVtYmVyLlBPU0lUSVZFX0lORklOSVRZKSB8fFxyXG4gICAgKGFjdHVhbCA9PT0gTnVtYmVyLk5FR0FUSVZFX0lORklOSVRZICYmIGV4cGVjdGVkID09PSBOdW1iZXIuTkVHQVRJVkVfSU5GSU5JVFkpIHx8XHJcbiAgICAoYWN0dWFsID09PSBOdW1iZXIuTmFOICYmIGV4cGVjdGVkID09PSBOdW1iZXIuTmFOKSB8fCAoaXNOYU4oYWN0dWFsKSAmJiBpc05hTihleHBlY3RlZCkpKVxyXG4gICAgcmV0dXJuO1xyXG4gIGNvbnN0IGFyZUVxdWFsID0gTWF0aC5hYnMoYWN0dWFsIC0gZXhwZWN0ZWQpIDwgdG9sZXJhbmNlO1xyXG4gIGV4cGVjdChhcmVFcXVhbCwgdHJ1ZSwgYCR7ZXJyb3IgPz8gJyd9ICh0b2xlcmFuY2UgPSAke3RvbGVyYW5jZX07IGEgPSAke2FjdHVhbH0sIGUgPSAke2V4cGVjdGVkfSlgKTtcclxuICBpZiAoIWFyZUVxdWFsKVxyXG4gICAgdGhyb3cgbmV3IEVycm9yKGBFeHBlY3RlZCAke2V4cGVjdGVkfSwgZ290ICR7YWN0dWFsfSAodG9sZXJhbmNlID0gJHt0b2xlcmFuY2V9KWApO1xyXG59XHJcblxyXG5leHBvcnQgZnVuY3Rpb24gZXhwZWN0VGFibGUoYWN0dWFsOiBfREcuRGF0YUZyYW1lLCBleHBlY3RlZDogX0RHLkRhdGFGcmFtZSwgZXJyb3I/OiBzdHJpbmcpOiB2b2lkIHtcclxuICBjb25zdCBleHBlY3RlZFJvd0NvdW50ID0gZXhwZWN0ZWQucm93Q291bnQ7XHJcbiAgY29uc3QgYWN0dWFsUm93Q291bnQgPSBhY3R1YWwucm93Q291bnQ7XHJcbiAgZXhwZWN0KGFjdHVhbFJvd0NvdW50LCBleHBlY3RlZFJvd0NvdW50LCBgJHtlcnJvciA/PyAnJ30sIHJvdyBjb3VudGApO1xyXG5cclxuICBmb3IgKGNvbnN0IGNvbHVtbiBvZiBleHBlY3RlZC5jb2x1bW5zKSB7XHJcbiAgICBjb25zdCBhY3R1YWxDb2x1bW4gPSBhY3R1YWwuY29sdW1ucy5ieU5hbWUoY29sdW1uLm5hbWUpO1xyXG4gICAgaWYgKGFjdHVhbENvbHVtbiA9PSBudWxsKVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYENvbHVtbiAke2NvbHVtbi5uYW1lfSBub3QgZm91bmRgKTtcclxuICAgIGlmIChhY3R1YWxDb2x1bW4udHlwZSAhPSBjb2x1bW4udHlwZSlcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKGBDb2x1bW4gJHtjb2x1bW4ubmFtZX0gdHlwZSBleHBlY3RlZCAke2NvbHVtbi50eXBlfSBnb3QgJHthY3R1YWxDb2x1bW4udHlwZX1gKTtcclxuICAgIGZvciAobGV0IGkgPSAwOyBpIDwgZXhwZWN0ZWRSb3dDb3VudDsgaSsrKSB7XHJcbiAgICAgIGNvbnN0IHZhbHVlID0gY29sdW1uLmdldChpKTtcclxuICAgICAgY29uc3QgYWN0dWFsVmFsdWUgPSBhY3R1YWxDb2x1bW4uZ2V0KGkpO1xyXG4gICAgICBpZiAoY29sdW1uLnR5cGUgPT0gREcuVFlQRS5GTE9BVClcclxuICAgICAgICBleHBlY3RGbG9hdChhY3R1YWxWYWx1ZSwgdmFsdWUsIDAuMDAwMSwgZXJyb3IpO1xyXG4gICAgICBlbHNlIGlmIChjb2x1bW4udHlwZSA9PSBERy5UWVBFLkRBVEVfVElNRSlcclxuICAgICAgICBleHBlY3QoYWN0dWFsVmFsdWUuaXNTYW1lKHZhbHVlKSwgdHJ1ZSwgZXJyb3IpO1xyXG4gICAgICBlbHNlXHJcbiAgICAgICAgZXhwZWN0KGFjdHVhbFZhbHVlLCB2YWx1ZSwgZXJyb3IpO1xyXG4gICAgfVxyXG4gIH1cclxufVxyXG5cclxuZXhwb3J0IGZ1bmN0aW9uIGV4cGVjdE9iamVjdChhY3R1YWw6IHsgW2tleTogc3RyaW5nXTogYW55IH0sIGV4cGVjdGVkOiB7IFtrZXk6IHN0cmluZ106IGFueSB9KSB7XHJcbiAgZm9yIChjb25zdCBbZXhwZWN0ZWRLZXksIGV4cGVjdGVkVmFsdWVdIG9mIE9iamVjdC5lbnRyaWVzKGV4cGVjdGVkKSkge1xyXG4gICAgaWYgKCFhY3R1YWwuaGFzT3duUHJvcGVydHkoZXhwZWN0ZWRLZXkpKVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYEV4cGVjdGVkIHByb3BlcnR5IFwiJHtleHBlY3RlZEtleX1cIiBub3QgZm91bmRgKTtcclxuXHJcbiAgICBjb25zdCBhY3R1YWxWYWx1ZSA9IGFjdHVhbFtleHBlY3RlZEtleV07XHJcbiAgICBpZiAoYWN0dWFsVmFsdWUgaW5zdGFuY2VvZiBBcnJheSAmJiBleHBlY3RlZFZhbHVlIGluc3RhbmNlb2YgQXJyYXkpXHJcbiAgICAgIGV4cGVjdEFycmF5KGFjdHVhbFZhbHVlLCBleHBlY3RlZFZhbHVlKTtcclxuICAgIGVsc2UgaWYgKGFjdHVhbFZhbHVlIGluc3RhbmNlb2YgT2JqZWN0ICYmIGV4cGVjdGVkVmFsdWUgaW5zdGFuY2VvZiBPYmplY3QpXHJcbiAgICAgIGV4cGVjdE9iamVjdChhY3R1YWxWYWx1ZSwgZXhwZWN0ZWRWYWx1ZSk7XHJcbiAgICBlbHNlIGlmIChOdW1iZXIuaXNGaW5pdGUoYWN0dWFsVmFsdWUpICYmIE51bWJlci5pc0Zpbml0ZShleHBlY3RlZFZhbHVlKSlcclxuICAgICAgZXhwZWN0RmxvYXQoYWN0dWFsVmFsdWUsIGV4cGVjdGVkVmFsdWUpO1xyXG4gICAgZWxzZSBpZiAoYWN0dWFsVmFsdWUgIT0gZXhwZWN0ZWRWYWx1ZSlcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKGBFeHBlY3RlZCAoJHtleHBlY3RlZFZhbHVlfSkgZm9yIGtleSAnJHtleHBlY3RlZEtleX0nLCBnb3QgKCR7YWN0dWFsVmFsdWV9KWApO1xyXG4gIH1cclxufVxyXG5cclxuZXhwb3J0IGZ1bmN0aW9uIGV4cGVjdEFycmF5KGFjdHVhbDogQXJyYXlMaWtlPGFueT4sIGV4cGVjdGVkOiBBcnJheUxpa2U8YW55Pikge1xyXG4gIGNvbnN0IGFjdHVhbExlbmd0aCA9IGFjdHVhbC5sZW5ndGg7XHJcbiAgY29uc3QgZXhwZWN0ZWRMZW5ndGggPSBleHBlY3RlZC5sZW5ndGg7XHJcblxyXG4gIGlmIChhY3R1YWxMZW5ndGggIT0gZXhwZWN0ZWRMZW5ndGgpIHtcclxuICAgIHRocm93IG5ldyBFcnJvcihgQXJyYXlzIGFyZSBvZiBkaWZmZXJlbnQgbGVuZ3RoOiBhY3R1YWwgYXJyYXkgbGVuZ3RoIGlzICR7YWN0dWFsTGVuZ3RofSBgICtcclxuICAgICAgYGFuZCBleHBlY3RlZCBhcnJheSBsZW5ndGggaXMgJHtleHBlY3RlZExlbmd0aH1gKTtcclxuICB9XHJcblxyXG4gIGZvciAobGV0IGkgPSAwOyBpIDwgYWN0dWFsTGVuZ3RoOyBpKyspIHtcclxuICAgIGlmIChhY3R1YWxbaV0gaW5zdGFuY2VvZiBBcnJheSAmJiBleHBlY3RlZFtpXSBpbnN0YW5jZW9mIEFycmF5KVxyXG4gICAgICBleHBlY3RBcnJheShhY3R1YWxbaV0sIGV4cGVjdGVkW2ldKTtcclxuICAgIGVsc2UgaWYgKGFjdHVhbFtpXSBpbnN0YW5jZW9mIE9iamVjdCAmJiBleHBlY3RlZFtpXSBpbnN0YW5jZW9mIE9iamVjdClcclxuICAgICAgZXhwZWN0T2JqZWN0KGFjdHVhbFtpXSwgZXhwZWN0ZWRbaV0pO1xyXG4gICAgZWxzZSBpZiAoYWN0dWFsW2ldICE9IGV4cGVjdGVkW2ldKVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoYEV4cGVjdGVkICR7ZXhwZWN0ZWRbaV19IGF0IHBvc2l0aW9uICR7aX0sIGdvdCAke2FjdHVhbFtpXX1gKTtcclxuICB9XHJcbn1cclxuXHJcbi8qIERlZmluZXMgYSB0ZXN0IHN1aXRlLiAqL1xyXG5leHBvcnQgZnVuY3Rpb24gY2F0ZWdvcnkoY2F0ZWdvcnk6IHN0cmluZywgdGVzdHNfOiAoKSA9PiB2b2lkLCBvcHRpb25zPzogQ2F0ZWdvcnlPcHRpb25zKTogdm9pZCB7XHJcbiAgY3VycmVudENhdGVnb3J5ID0gY2F0ZWdvcnk7XHJcbiAgdGVzdHNfKCk7XHJcbiAgaWYgKHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0pIHtcclxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0uY2xlYXIgPSBvcHRpb25zPy5jbGVhciA/PyB0cnVlO1xyXG4gICAgdGVzdHNbY3VycmVudENhdGVnb3J5XS50aW1lb3V0ID0gb3B0aW9ucz8udGltZW91dDtcclxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0uYmVuY2htYXJrcyA9IG9wdGlvbnM/LmJlbmNobWFya3M7XHJcbiAgICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldLnN0cmVzc1Rlc3RzID0gb3B0aW9ucz8uc3RyZXNzVGVzdHM7XHJcbiAgICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldLm93bmVyID0gb3B0aW9ucz8ub3duZXI7XHJcbiAgfVxyXG59XHJcblxyXG4vKiBEZWZpbmVzIGEgZnVuY3Rpb24gdG8gYmUgZXhlY3V0ZWQgYmVmb3JlIHRoZSB0ZXN0cyBpbiB0aGlzIGNhdGVnb3J5IGFyZSBleGVjdXRlZC4gKi9cclxuZXhwb3J0IGZ1bmN0aW9uIGJlZm9yZShiZWZvcmU6ICgpID0+IFByb21pc2U8dm9pZD4pOiB2b2lkIHtcclxuICBpZiAodGVzdHNbY3VycmVudENhdGVnb3J5XSA9PSB1bmRlZmluZWQpXHJcbiAgICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldID0ge307XHJcbiAgdGVzdHNbY3VycmVudENhdGVnb3J5XS5iZWZvcmUgPSBiZWZvcmU7XHJcbn1cclxuXHJcbi8qIERlZmluZXMgYSBmdW5jdGlvbiB0byBiZSBleGVjdXRlZCBhZnRlciB0aGUgdGVzdHMgaW4gdGhpcyBjYXRlZ29yeSBhcmUgZXhlY3V0ZWQuICovXHJcbmV4cG9ydCBmdW5jdGlvbiBhZnRlcihhZnRlcjogKCkgPT4gUHJvbWlzZTx2b2lkPik6IHZvaWQge1xyXG4gIGlmICh0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldID09IHVuZGVmaW5lZClcclxuICAgIHRlc3RzW2N1cnJlbnRDYXRlZ29yeV0gPSB7fTtcclxuICB0ZXN0c1tjdXJyZW50Q2F0ZWdvcnldLmFmdGVyID0gYWZ0ZXI7XHJcbn1cclxuXHJcbmZ1bmN0aW9uIGFkZE5hbWVzcGFjZShzOiBzdHJpbmcsIGY6IF9ERy5GdW5jKTogc3RyaW5nIHtcclxuICByZXR1cm4gcy5yZXBsYWNlKG5ldyBSZWdFeHAoZi5uYW1lLCAnZ2knKSwgZi5ucU5hbWUpO1xyXG59XHJcblxyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gaW5pdEF1dG9UZXN0cyhwYWNrYWdlXzogX0RHLlBhY2thZ2UsIG1vZHVsZT86IGFueSkge1xyXG4gIGNvbnN0IHBhY2thZ2VJZCA9IHBhY2thZ2VfLmlkO1xyXG4gIGlmICh3YXNSZWdpc3RlcmVkW3BhY2thZ2VJZF0pIHJldHVybjtcclxuICBjb25zdCBtb2R1bGVUZXN0cyA9IG1vZHVsZSA/IG1vZHVsZS50ZXN0cyA6IHRlc3RzO1xyXG4gIGlmIChwYWNrYWdlXy5uYW1lID09PSAnRGV2VG9vbHMnIHx8ICghIW1vZHVsZSAmJiBtb2R1bGUuX3BhY2thZ2UubmFtZSA9PT0gJ0RldlRvb2xzJykpIHtcclxuICAgIGZvciAoY29uc3QgZiBvZiAoPGFueT53aW5kb3cpLmRhcnRUZXN0cykge1xyXG4gICAgICBjb25zdCBhcnIgPSBmLm5hbWUuc3BsaXQoL1xccypcXHxcXHMqIS9nKTtcclxuICAgICAgbGV0IG5hbWUgPSBhcnIucG9wKCkgPz8gZi5uYW1lO1xyXG4gICAgICBsZXQgY2F0ID0gYXJyLmxlbmd0aCA/IGNvcmVDYXROYW1lICsgJzogJyArIGFyci5qb2luKCc6ICcpIDogY29yZUNhdE5hbWU7XHJcbiAgICAgIGxldCBmdWxsTmFtZTogc3RyaW5nW10gPSBuYW1lLnNwbGl0KCcgfCAnKTtcclxuICAgICAgbmFtZSA9IGZ1bGxOYW1lW2Z1bGxOYW1lLmxlbmd0aCAtIDFdO1xyXG4gICAgICBmdWxsTmFtZS51bnNoaWZ0KGNhdCk7XHJcbiAgICAgIGZ1bGxOYW1lLnBvcCgpO1xyXG4gICAgICBjYXQgPSBmdWxsTmFtZS5qb2luKCc6ICcpO1xyXG4gICAgICBpZiAobW9kdWxlVGVzdHNbY2F0XSA9PT0gdW5kZWZpbmVkKVxyXG4gICAgICAgIG1vZHVsZVRlc3RzW2NhdF0gPSB7IHRlc3RzOiBbXSwgY2xlYXI6IHRydWUgfTtcclxuICAgICAgbW9kdWxlVGVzdHNbY2F0XS50ZXN0cy5wdXNoKG5ldyBUZXN0KGNhdCwgbmFtZSwgZi50ZXN0LCB7IGlzQWdncmVnYXRlZDogZmFsc2UsIHRpbWVvdXQ6IGYub3B0aW9ucz8udGltZW91dCA/PyBTVEFOREFSVF9USU1FT1VULCBza2lwUmVhc29uOiBmLm9wdGlvbnM/LnNraXBSZWFzb24sIG93bmVyOiBmLm9wdGlvbnM/Lm93bmVyLCBiZW5jaG1hcms6IGYub3B0aW9ucz8uYmVuY2htYXJrID8/IGZhbHNlIH0pKTtcclxuICAgIH1cclxuICB9XHJcbiAgY29uc3QgbW9kdWxlQXV0b1Rlc3RzID0gW107XHJcbiAgY29uc3QgbW9kdWxlRGVtbyA9IFtdO1xyXG4gIGNvbnN0IG1vZHVsZURldGVjdG9ycyA9IFtdO1xyXG4gIGNvbnN0IHBhY2tGdW5jdGlvbnMgPSBhd2FpdCBncm9rLmRhcGkuZnVuY3Rpb25zLmZpbHRlcihgcGFja2FnZS5pZCA9IFwiJHtwYWNrYWdlSWR9XCJgKS5saXN0KCk7XHJcbiAgY29uc3QgcmVnID0gbmV3IFJlZ0V4cCgvc2tpcDpcXHMqKFteLFxcc10rKXx3YWl0OlxccyooXFxkKyl8Y2F0OlxccyooW14sXFxzXSspfHRpbWVvdXQ6XFxzKihcXGQrKS9nKTtcclxuICBmb3IgKGNvbnN0IGYgb2YgcGFja0Z1bmN0aW9ucykge1xyXG4gICAgY29uc3QgdGVzdHMgPSBmLm9wdGlvbnNbJ3Rlc3QnXTtcclxuICAgIGNvbnN0IGRlbW8gPSBmLm9wdGlvbnNbJ2RlbW9QYXRoJ107XHJcbiAgICBpZiAoKHRlc3RzICYmIEFycmF5LmlzQXJyYXkodGVzdHMpICYmIHRlc3RzLmxlbmd0aCkpIHtcclxuICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCB0ZXN0cy5sZW5ndGg7IGkrKykge1xyXG4gICAgICAgIGNvbnN0IHJlcyA9ICh0ZXN0c1tpXSBhcyBzdHJpbmcpLm1hdGNoQWxsKHJlZyk7XHJcbiAgICAgICAgY29uc3QgbWFwOiB7IHNraXA/OiBzdHJpbmcsIHdhaXQ/OiBudW1iZXIsIGNhdD86IHN0cmluZywgdGltZW91dD86IG51bWJlciwgYmVuY2htYXJrVGltZW91dD86IG51bWJlciB9ID0ge307XHJcbiAgICAgICAgQXJyYXkuZnJvbShyZXMpLmZvckVhY2goKGFycikgPT4ge1xyXG4gICAgICAgICAgaWYgKGFyclswXS5zdGFydHNXaXRoKCdza2lwJykpIG1hcFsnc2tpcCddID0gYXJyWzFdO1xyXG4gICAgICAgICAgZWxzZSBpZiAoYXJyWzBdLnN0YXJ0c1dpdGgoJ3dhaXQnKSkgbWFwWyd3YWl0J10gPSBwYXJzZUludChhcnJbMl0pO1xyXG4gICAgICAgICAgZWxzZSBpZiAoYXJyWzBdLnN0YXJ0c1dpdGgoJ2NhdCcpKSBtYXBbJ2NhdCddID0gYXJyWzNdO1xyXG4gICAgICAgICAgZWxzZSBpZiAoYXJyWzBdLnN0YXJ0c1dpdGgoJ3RpbWVvdXQnKSkgbWFwWyd0aW1lb3V0J10gPSBwYXJzZUludChhcnJbNF0pO1xyXG4gICAgICAgIH0pO1xyXG4gICAgICAgIGNvbnN0IHRlc3QgPSBuZXcgVGVzdChtYXAuY2F0ID8/IGF1dG9UZXN0c0NhdE5hbWUsIHRlc3RzLmxlbmd0aCA9PT0gMSA/IGYubmFtZSA6IGAke2YubmFtZX0gJHtpICsgMX1gLCBhc3luYyAoKSA9PiB7XHJcbiAgICAgICAgICBjb25zdCByZXMgPSBhd2FpdCBncm9rLmZ1bmN0aW9ucy5ldmFsKGFkZE5hbWVzcGFjZSh0ZXN0c1tpXSwgZikpO1xyXG4gICAgICAgICAgaWYgKG1hcC53YWl0KSBhd2FpdCBkZWxheShtYXAud2FpdCk7XHJcbiAgICAgICAgICAvLyBlc2xpbnQtZGlzYWJsZS1uZXh0LWxpbmUgbm8tdGhyb3ctbGl0ZXJhbFxyXG4gICAgICAgICAgaWYgKHR5cGVvZiByZXMgPT09ICdib29sZWFuJyAmJiAhcmVzKSB0aHJvdyBgRmFpbGVkOiAke3Rlc3RzW2ldfSwgZXhwZWN0ZWQgdHJ1ZSwgZ290ICR7cmVzfWA7XHJcbiAgICAgICAgfSwgeyBza2lwUmVhc29uOiBtYXAuc2tpcCwgdGltZW91dDogREcuVGVzdC5pc0luQmVuY2htYXJrID8gbWFwLmJlbmNobWFya1RpbWVvdXQgPz8gQkVOQ0hNQVJLX1RJTUVPVVQgOiBtYXAudGltZW91dCA/PyBTVEFOREFSVF9USU1FT1VUIH0pO1xyXG4gICAgICAgIGlmIChtYXAuY2F0KSB7XHJcbiAgICAgICAgICBjb25zdCBjYXQ6IHN0cmluZyA9IG1hcC5jYXQ7XHJcbiAgICAgICAgICBpZiAobW9kdWxlVGVzdHNbY2F0XSA9PT0gdW5kZWZpbmVkKVxyXG4gICAgICAgICAgICBtb2R1bGVUZXN0c1tjYXRdID0geyB0ZXN0czogW10sIGNsZWFyOiB0cnVlIH07XHJcblxyXG4gICAgICAgICAgLy8gb25seSBiZWZvcmUvYWZ0ZXIgY2FuIGJlIGRlZmluZWQgaW4gdHMgZmlsZXMgdGVzdHMgdW5kZXIgdGhlIGNhdGVnb3J5XHJcbiAgICAgICAgICBpZiAoIW1vZHVsZVRlc3RzW2NhdF0udGVzdHMpXHJcbiAgICAgICAgICAgIG1vZHVsZVRlc3RzW2NhdF0udGVzdHMgPSBbXTtcclxuICAgICAgICAgIG1vZHVsZVRlc3RzW2NhdF0udGVzdHMucHVzaCh0ZXN0KTtcclxuICAgICAgICB9XHJcbiAgICAgICAgZWxzZVxyXG4gICAgICAgICAgbW9kdWxlQXV0b1Rlc3RzLnB1c2godGVzdCk7XHJcbiAgICAgIH1cclxuICAgIH1cclxuICAgIGlmIChkZW1vKSB7XHJcbiAgICAgIGNvbnN0IHdhaXQgPSBmLm9wdGlvbnNbJ2RlbW9XYWl0J10gPyBwYXJzZUludChmLm9wdGlvbnNbJ2RlbW9XYWl0J10pIDogdW5kZWZpbmVkO1xyXG4gICAgICBjb25zdCB0ZXN0ID0gbmV3IFRlc3QoZGVtb0NhdE5hbWUsIGYuZnJpZW5kbHlOYW1lLCBhc3luYyAoKSA9PiB7XHJcbiAgICAgICAgYXdhaXQgZGVsYXkoMzAwKTtcclxuICAgICAgICBncm9rLnNoZWxsLmNsZWFyTGFzdEVycm9yKCk7XHJcbiAgICAgICAgYXdhaXQgZi5hcHBseSgpO1xyXG4gICAgICAgIGF3YWl0IGRlbGF5KHdhaXQgPyB3YWl0IDogMjAwMCk7XHJcbiAgICAgICAgY29uc3QgdW5oYW5kbGVkID0gYXdhaXQgZ3Jvay5zaGVsbC5sYXN0RXJyb3I7XHJcbiAgICAgICAgaWYgKHVuaGFuZGxlZClcclxuICAgICAgICAgIHRocm93IG5ldyBFcnJvcih1bmhhbmRsZWQpO1xyXG4gICAgICB9LCB7IHNraXBSZWFzb246IGYub3B0aW9uc1snZGVtb1NraXAnXSB9KTtcclxuICAgICAgbW9kdWxlRGVtby5wdXNoKHRlc3QpO1xyXG4gICAgfVxyXG4gICAgaWYgKGYuaGFzVGFnKCdzZW1UeXBlRGV0ZWN0b3InKSkge1xyXG4gICAgICBsZXQgZGV0ZWN0b3JzVGVzdERhdGEgPSB0ZXN0RGF0YTtcclxuICAgICAgaWYgKGYub3B0aW9uc1sndGVzdERhdGEnXSkge1xyXG4gICAgICAgIGRldGVjdG9yc1Rlc3REYXRhID0gYXdhaXQgZ3Jvay5kYXRhLmZpbGVzLm9wZW5UYWJsZShgU3lzdGVtOkFwcERhdGEvJHtwYWNrYWdlXy5ucU5hbWV9LyR7Zi5vcHRpb25zWyd0ZXN0RGF0YSddfWApO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBjb25zdCB0ZXN0ID0gbmV3IFRlc3QoZGV0ZWN0b3JzQ2F0TmFtZSwgZi5mcmllbmRseU5hbWUsIGFzeW5jICgpID0+IHtcclxuICAgICAgICBjb25zdCBhcnIgPSBbXTtcclxuICAgICAgICBjb25zb2xlLmxvZyhgU3lzdGVtOkFwcERhdGEvJHtwYWNrYWdlXy5ucU5hbWV9LyR7Zi5vcHRpb25zWyd0ZXN0RGF0YSddfWApO1xyXG5cclxuICAgICAgICBmb3IgKGNvbnN0IGNvbCBvZiBkZXRlY3RvcnNUZXN0RGF0YS5jbG9uZSgpLmNvbHVtbnMpIHtcclxuICAgICAgICAgIGNvbnN0IHJlcyA9IGF3YWl0IGYuYXBwbHkoW2NvbF0pO1xyXG4gICAgICAgICAgYXJyLnB1c2gocmVzIHx8IGNvbC5zZW1UeXBlKTtcclxuICAgICAgICB9XHJcbiAgICAgICAgY29uc3QgcmVzQXJyID0gYXJyLmZpbHRlcigoaSkgPT4gaSk7XHJcbiAgICAgICAgZXhwZWN0KHJlc0Fyci5sZW5ndGgsIDEpO1xyXG5cclxuICAgICAgICBpZiAoZi5vcHRpb25zWyd0ZXN0RGF0YUNvbHVtbk5hbWUnXSlcclxuICAgICAgICAgIGV4cGVjdChyZXNBcnJbMF0sIGYub3B0aW9uc1sndGVzdERhdGFDb2x1bW5OYW1lJ10pO1xyXG5cclxuICAgICAgfSwgeyBza2lwUmVhc29uOiBmLm9wdGlvbnNbJ3NraXBUZXN0J10gfSk7XHJcbiAgICAgIG1vZHVsZURldGVjdG9ycy5wdXNoKHRlc3QpO1xyXG4gICAgfVxyXG4gIH1cclxuICB3YXNSZWdpc3RlcmVkW3BhY2thZ2VJZF0gPSB0cnVlO1xyXG4gIGlmIChtb2R1bGVBdXRvVGVzdHMubGVuZ3RoID4gMClcclxuICAgIG1vZHVsZVRlc3RzW2F1dG9UZXN0c0NhdE5hbWVdID0geyB0ZXN0czogbW9kdWxlQXV0b1Rlc3RzLCBjbGVhcjogdHJ1ZSB9O1xyXG4gIGlmIChtb2R1bGVEZW1vLmxlbmd0aCA+IDApXHJcbiAgICBtb2R1bGVUZXN0c1tkZW1vQ2F0TmFtZV0gPSB7IHRlc3RzOiBtb2R1bGVEZW1vLCBjbGVhcjogdHJ1ZSB9O1xyXG4gIGlmIChtb2R1bGVEZXRlY3RvcnMubGVuZ3RoID4gMClcclxuICAgIG1vZHVsZVRlc3RzW2RldGVjdG9yc0NhdE5hbWVdID0geyB0ZXN0czogbW9kdWxlRGV0ZWN0b3JzLCBjbGVhcjogZmFsc2UgfTtcclxufVxyXG5cclxuZnVuY3Rpb24gcmVkZWZpbmVDb25zb2xlKCk6IGFueVtdIHtcclxuICBjb25zdCBsb2dzOiBhbnlbXSA9IFtdO1xyXG4gIGNvbnNvbGUubG9nID0gKC4uLmFyZ3MpID0+IHtcclxuICAgIGxvZ3MucHVzaCguLi5hcmdzKTtcclxuICAgIHN0ZExvZyguLi5hcmdzKTtcclxuICB9O1xyXG4gIGNvbnNvbGUuaW5mbyA9ICguLi5hcmdzKSA9PiB7XHJcbiAgICBsb2dzLnB1c2goLi4uYXJncyk7XHJcbiAgICBzdGRJbmZvKC4uLmFyZ3MpO1xyXG4gIH07XHJcbiAgY29uc29sZS53YXJuID0gKC4uLmFyZ3MpID0+IHtcclxuICAgIGxvZ3MucHVzaCguLi5hcmdzKTtcclxuICAgIHN0ZFdhcm4oLi4uYXJncyk7XHJcbiAgfTtcclxuICBjb25zb2xlLmVycm9yID0gKC4uLmFyZ3MpID0+IHtcclxuICAgIGxvZ3MucHVzaCguLi5hcmdzKTtcclxuICAgIHN0ZEVycm9yKC4uLmFyZ3MpO1xyXG4gIH07XHJcbiAgcmV0dXJuIGxvZ3M7XHJcbn1cclxuXHJcbmZ1bmN0aW9uIHJlc2V0Q29uc29sZSgpOiB2b2lkIHtcclxuICBjb25zb2xlLmxvZyA9IHN0ZExvZztcclxuICBjb25zb2xlLmluZm8gPSBzdGRJbmZvO1xyXG4gIGNvbnNvbGUud2FybiA9IHN0ZFdhcm47XHJcbiAgY29uc29sZS5lcnJvciA9IHN0ZEVycm9yO1xyXG59XHJcblxyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gcnVuVGVzdHMob3B0aW9ucz86IFRlc3RFeGVjdXRpb25PcHRpb25zKSA6IFByb21pc2U8VGVzdFJlc3VsdEV4dGVuZGVkW10+e1xyXG5cclxuICBjb25zdCBwYWNrYWdlXzogX0RHLlBhY2thZ2UgPSBvcHRpb25zPy5ub2RlT3B0aW9ucyA/IG9wdGlvbnMubm9kZU9wdGlvbnMucGFja2FnZSA6IGdyb2suZnVuY3Rpb25zLmdldEN1cnJlbnRDYWxsKCkuZnVuYy5wYWNrYWdlO1xyXG4gIGlmICghcGFja2FnZV8pXHJcbiAgICB0aHJvdyBuZXcgRXJyb3IoJ0NhblxcJ3QgcnVuIHRlc3RzIG91dHNpZGUgb2YgdGhlIHBhY2thZ2UnKTtcclxuICBjb25zdCBtYXRjaCA9IHBhY2thZ2VfLnBhY2thZ2VPd25lcj8ubWF0Y2goLzwoW14+XSopPi8pO1xyXG4gIGNvbnN0IHBhY2thZ2VPd25lciA9IG1hdGNoID8gbWF0Y2hbMV0gOiAnJztcclxuICBpZiAocGFja2FnZV8gIT0gdW5kZWZpbmVkKVxyXG4gICAgYXdhaXQgaW5pdEF1dG9UZXN0cyhwYWNrYWdlXyk7XHJcbiAgY29uc3QgcmVzdWx0czpUZXN0UmVzdWx0RXh0ZW5kZWRbXSA9IFtdO1xyXG4gIGNvbnNvbGUubG9nKGBSdW5uaW5nIHRlc3RzLi4uYCk7XHJcbiAgY29uc29sZS5sb2cob3B0aW9ucyk7XHJcbiAgb3B0aW9ucyA/Pz0ge307XHJcbiAgb3B0aW9ucyEudGVzdENvbnRleHQgPz89IG5ldyBUZXN0Q29udGV4dCgpO1xyXG4gIGdyb2suc2hlbGwuY2xlYXJMYXN0RXJyb3IoKTtcclxuICBjb25zdCBsb2dzID0gcmVkZWZpbmVDb25zb2xlKCk7XHJcblxyXG4gIGF3YWl0IGludm9rZVRlc3RzKHRlc3RzLCBvcHRpb25zKTtcclxuXHJcbiAgZm9yIChsZXQgciBvZiByZXN1bHRzKSB7XHJcbiAgICByLnJlc3VsdCA9IHIucmVzdWx0LnRvU3RyaW5nKCkucmVwbGFjZSgvXCIvZywgJ1xcJycpO1xyXG4gICAgaWYgKHIubG9ncyAhPSB1bmRlZmluZWQpXHJcbiAgICAgIHIubG9ncyA9IHIubG9ncyEudG9TdHJpbmcoKS5yZXBsYWNlKC9cIi9nLCAnXFwnJyk7XHJcbiAgfVxyXG4gIHJldHVybiByZXN1bHRzO1xyXG5cclxuICBhc3luYyBmdW5jdGlvbiBpbnZva2VDYXRlZ29yeU1ldGhvZChtZXRob2Q6ICgoKSA9PiBQcm9taXNlPHZvaWQ+KSB8IHVuZGVmaW5lZCwgY2F0ZWdvcnk6IHN0cmluZyk6IFByb21pc2U8c3RyaW5nIHwgdW5kZWZpbmVkPiB7XHJcbiAgICBsZXQgaW52b2thdGlvblJlc3VsdCA9IHVuZGVmaW5lZDtcclxuICAgIHRyeSB7XHJcbiAgICAgIGlmIChtZXRob2QgIT09IHVuZGVmaW5lZCkge1xyXG4gICAgICAgIGF3YWl0IHRpbWVvdXQoYXN5bmMgKCkgPT4ge1xyXG4gICAgICAgICAgYXdhaXQgbWV0aG9kKCk7XHJcbiAgICAgICAgfSwgMTAwMDAwLCBgYmVmb3JlICR7Y2F0ZWdvcnl9OiB0aW1lb3V0IGVycm9yYCk7XHJcbiAgICAgIH1cclxuICAgIH0gY2F0Y2ggKHg6IGFueSkge1xyXG4gICAgICBpbnZva2F0aW9uUmVzdWx0ID0gYXdhaXQgZ2V0UmVzdWx0KHgpO1xyXG4gICAgfVxyXG4gICAgcmV0dXJuIGludm9rYXRpb25SZXN1bHRcclxuICB9XHJcblxyXG4gIGFzeW5jIGZ1bmN0aW9uIGludm9rZVRlc3RzSW5DYXRlZ29yeShjYXRlZ29yeTogQ2F0ZWdvcnksIG9wdGlvbnM6IFRlc3RFeGVjdXRpb25PcHRpb25zLCBpc1RhcmdldENhdGVnb3J5OiBib29sZWFuKTogUHJvbWlzZTxUZXN0UmVzdWx0RXh0ZW5kZWRbXT4ge1xyXG4gICAgbGV0IHQgPSBjYXRlZ29yeS50ZXN0cyA/PyBbXTtcclxuICAgIGNvbnN0IHJlcyA6IFRlc3RSZXN1bHRFeHRlbmRlZFtdID0gW107XHJcbiAgICAvLyBsZXQgbWVtb3J5VXNhZ2VCZWZvcmUgPSAod2luZG93Py5wZXJmb3JtYW5jZSBhcyBhbnkpPy5tZW1vcnk/LnVzZWRKU0hlYXBTaXplO1xyXG4gICAgY29uc3Qgd2lkZ2V0c0JlZm9yZSA9IGdldFdpZGdldHNDb3VudFNhZmUoKTtcclxuXHJcbiAgICBpZiAoY2F0ZWdvcnkuY2xlYXIpIHtcclxuICAgICAgICBsZXQgc2tpcHBpbmdUZXN0cyA9IGlzVGFyZ2V0Q2F0ZWdvcnkgJiYgb3B0aW9ucy5za2lwVG9UZXN0ICE9IHVuZGVmaW5lZDtcclxuICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCB0Lmxlbmd0aDsgaSsrKSB7XHJcblxyXG4gICAgICAgIGlmICh0W2ldLm9wdGlvbnMpIHtcclxuICAgICAgICAgIGlmICh0W2ldLm9wdGlvbnM/LmJlbmNobWFyayA9PT0gdW5kZWZpbmVkKSB7XHJcbiAgICAgICAgICAgIGlmICghdFtpXS5vcHRpb25zKVxyXG4gICAgICAgICAgICAgIHRbaV0ub3B0aW9ucyA9IHt9XHJcbiAgICAgICAgICAgIHRbaV0ub3B0aW9ucyEuYmVuY2htYXJrID0gY2F0ZWdvcnkuYmVuY2htYXJrcyA/PyBmYWxzZTtcclxuICAgICAgICAgIH1cclxuICAgICAgICB9XHJcbiAgICAgICAgbGV0IHRlc3QgPSB0W2ldO1xyXG4gICAgICAgIGlmIChvcHRpb25zLnRlc3QpXHJcbiAgICAgICAgICBpZiAob3B0aW9ucy50ZXN0LnRvTG93ZXJDYXNlKCkgIT09IHRlc3QubmFtZS50b0xvd2VyQ2FzZSgpKVxyXG4gICAgICAgICAgICBjb250aW51ZTtcclxuICAgICAgICBpZiAoc2tpcHBpbmdUZXN0cykge1xyXG4gICAgICAgICAgaWYgKG9wdGlvbnM/LnNraXBUb1Rlc3QgIT0gdW5kZWZpbmVkICYmIHRlc3QubmFtZS50b0xvd2VyQ2FzZSgpLnRyaW0oKSA9PT0gb3B0aW9ucz8uc2tpcFRvVGVzdC50b0xvd2VyQ2FzZSgpLnRyaW0oKSkge1xyXG4gICAgICAgICAgICAvLyBGb3VuZCB0aGUgdGFyZ2V0IHRlc3QsIHN0b3Agc2tpcHBpbmcgYWZ0ZXIgdGhpcyBvbmVcclxuICAgICAgICAgICAgc2tpcHBpbmdUZXN0cyA9IGZhbHNlO1xyXG4gICAgICAgICAgfSBlbHNlXHJcbiAgICAgICAgICBjb250aW51ZTtcclxuICAgICAgICB9XHJcbiAgICAgICAgaWYgKHRlc3Q/Lm9wdGlvbnMpIHtcclxuICAgICAgICAgIHRlc3Qub3B0aW9ucy5vd25lciA9IHRbaV0ub3B0aW9ucz8ub3duZXIgPz8gY2F0ZWdvcnk/Lm93bmVyID8/IHBhY2thZ2VPd25lciA/PyAnJztcclxuICAgICAgICB9XHJcbiAgICAgICAgLy8gbGV0IGlzR0JFbmFibGUgPSAod2luZG93IGFzIGFueSkuZ2MgJiYgdGVzdC5vcHRpb25zPy5za2lwUmVhc29uID09IHVuZGVmaW5lZDtcclxuICAgICAgICAvLyBjb25zb2xlLmxvZyhgKioqKioqKioke2lzR0JFbmFibGV9YCk7XHJcbiAgICAgICAgLy8gaWYgKGlzR0JFbmFibGUpXHJcbiAgICAgICAgLy8gICBhd2FpdCAod2luZG93IGFzIGFueSkuZ2MoKTtcclxuICAgICAgICAvLyBtZW1vcnlVc2FnZUJlZm9yZSA9ICh3aW5kb3c/LnBlcmZvcm1hbmNlIGFzIGFueSk/Lm1lbW9yeT8udXNlZEpTSGVhcFNpemU7XHJcbiAgICAgICAgbGV0IHRlc3RSdW4gPSBhd2FpdCBleGVjVGVzdChcclxuICAgICAgICAgICAgdGVzdCxcclxuICAgICAgICAgICAgb3B0aW9ucz8udGVzdCxcclxuICAgICAgICAgICAgbG9ncywgREcuVGVzdC5pc0luQmVuY2htYXJrID8gdFtpXS5vcHRpb25zPy5iZW5jaG1hcmtUaW1lb3V0ID8/IEJFTkNITUFSS19USU1FT1VUIDogdFtpXS5vcHRpb25zPy50aW1lb3V0ID8/IFNUQU5EQVJUX1RJTUVPVVQsXHJcbiAgICAgICAgICAgIHBhY2thZ2VfLm5hbWUsXHJcbiAgICAgICAgICAgIG9wdGlvbnMudmVyYm9zZVxyXG4gICAgICAgICk7XHJcblxyXG4gICAgICAgIC8vIGlmIChpc0dCRW5hYmxlKVxyXG4gICAgICAgIC8vICAgYXdhaXQgKHdpbmRvdyBhcyBhbnkpLmdjKCk7XHJcbiAgICAgICAgaWYgKHRlc3RSdW4pIHtcclxuICAgICAgICAgIHJlcy5wdXNoKHsgLi4udGVzdFJ1biwgIHdpZGdldHNEaWZmZXJlbmNlOiBnZXRXaWRnZXRzQ291bnRTYWZlKCkgLSB3aWRnZXRzQmVmb3JlIH0pO1xyXG4gICAgICAgICAgLy8gUmV0dXJuIGVhcmx5IGlmIHJldHVybk9uRmFpbCBpcyBzZXQgYW5kIHRlc3QgZmFpbGVkIChidXQgaWdub3JlIGZhaWx1cmUgZm9yIHRoZSBza2lwVG9UZXN0IHRlc3QgaXRzZWxmKVxyXG4gICAgICAgICAgaWYgKG9wdGlvbnMucmV0dXJuT25GYWlsICYmIG9wdGlvbnMuc2tpcFRvVGVzdCAhPT0gdGVzdC5uYW1lICYmICF0ZXN0UnVuLnN1Y2Nlc3MgJiYgIXRlc3RSdW4uc2tpcHBlZClcclxuICAgICAgICAgICAgcmV0dXJuIHJlcztcclxuICAgICAgICB9XHJcbiAgICAgICAgLy8gcmVzLnB1c2goeyAuLi50ZXN0UnVuLCBtZW1vcnlEZWx0YTogKHdpbmRvdz8ucGVyZm9ybWFuY2UgYXMgYW55KT8ubWVtb3J5Py51c2VkSlNIZWFwU2l6ZSAtIG1lbW9yeVVzYWdlQmVmb3JlLCB3aWRnZXRzRGVsdGE6IGdldFdpZGdldHNDb3VudFNhZmUoKSAtIHdpZGdldHNCZWZvcmUgfSk7XHJcblxyXG4gICAgICAgIGlmICghb3B0aW9ucy5ub2RlT3B0aW9ucykge1xyXG4gICAgICAgICAgZ3Jvay5zaGVsbC5jbG9zZUFsbCgpO1xyXG4gICAgICAgICAgREcuQmFsbG9vbi5jbG9zZUFsbCgpO1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG4gICAgfSBlbHNlIHtcclxuICAgICAgbGV0IHNraXBwaW5nVGVzdHMgPSBpc1RhcmdldENhdGVnb3J5ICYmIG9wdGlvbnMuc2tpcFRvVGVzdCAhPSB1bmRlZmluZWQ7XHJcbiAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgdC5sZW5ndGg7IGkrKykge1xyXG4gICAgICAgIGxldCB0ZXN0ID0gdFtpXTtcclxuICAgICAgICBpZiAob3B0aW9ucy50ZXN0KVxyXG4gICAgICAgICAgaWYgKG9wdGlvbnMudGVzdC50b0xvd2VyQ2FzZSgpICE9PSB0ZXN0Lm5hbWUudG9Mb3dlckNhc2UoKSlcclxuICAgICAgICAgICAgY29udGludWU7XHJcbiAgICAgICAgaWYgKHNraXBwaW5nVGVzdHMpIHtcclxuICAgICAgICAgIGlmIChvcHRpb25zPy5za2lwVG9UZXN0ICE9IHVuZGVmaW5lZCAmJiB0ZXN0Lm5hbWUudG9Mb3dlckNhc2UoKS50cmltKCkgPT09IG9wdGlvbnM/LnNraXBUb1Rlc3QudG9Mb3dlckNhc2UoKS50cmltKCkpIHtcclxuICAgICAgICAgICAgLy8gRm91bmQgdGhlIHRhcmdldCB0ZXN0LCBzdG9wIHNraXBwaW5nIGFmdGVyIHRoaXMgb25lXHJcbiAgICAgICAgICAgIHNraXBwaW5nVGVzdHMgPSBmYWxzZTtcclxuICAgICAgICAgIH1cclxuICAgICAgICAgIGNvbnRpbnVlOyAgLy8gU2tpcCB0aGlzIHRlc3QgKGluY2x1ZGluZyB0aGUgdGFyZ2V0KVxyXG4gICAgICAgIH1cclxuXHJcbiAgICAgICAgaWYgKHRlc3Q/Lm9wdGlvbnMpIHtcclxuICAgICAgICAgIHRlc3Qub3B0aW9ucy5vd25lciA9IHRbaV0ub3B0aW9ucz8ub3duZXIgPz8gY2F0ZWdvcnk/Lm93bmVyID8/IHBhY2thZ2VPd25lciA/PyAnJztcclxuICAgICAgICB9XHJcbiAgICAgICAgLy8gbGV0IGlzR0JFbmFibGUgPSAod2luZG93IGFzIGFueSkuZ2MgJiYgdGVzdC5vcHRpb25zPy5za2lwUmVhc29uID09IHVuZGVmaW5lZDtcclxuICAgICAgICAvLyBjb25zb2xlLmxvZyhgKioqKioqKioke2lzR0JFbmFibGV9YCk7XHJcbiAgICAgICAgLy8gaWYgKGlzR0JFbmFibGUpXHJcbiAgICAgICAgLy8gICBhd2FpdCAod2luZG93IGFzIGFueSkuZ2MoKTtcclxuICAgICAgICAvLyBtZW1vcnlVc2FnZUJlZm9yZSA9ICh3aW5kb3c/LnBlcmZvcm1hbmNlIGFzIGFueSk/Lm1lbW9yeT8udXNlZEpTSGVhcFNpemU7XHJcbiAgICAgICAgbGV0IHRlc3RSdW4gPSBhd2FpdCBleGVjVGVzdChcclxuICAgICAgICAgICAgdGVzdCxcclxuICAgICAgICAgICAgb3B0aW9ucz8udGVzdCxcclxuICAgICAgICAgICAgbG9ncyxcclxuICAgICAgICAgICAgREcuVGVzdC5pc0luQmVuY2htYXJrID8gdFtpXS5vcHRpb25zPy5iZW5jaG1hcmtUaW1lb3V0ID8/IEJFTkNITUFSS19USU1FT1VUIDogdFtpXS5vcHRpb25zPy50aW1lb3V0LFxyXG4gICAgICAgICAgICBwYWNrYWdlXy5uYW1lLFxyXG4gICAgICAgICAgICBvcHRpb25zLnZlcmJvc2VcclxuICAgICAgICApO1xyXG5cclxuICAgICAgICAvLyBpZiAoaXNHQkVuYWJsZSlcclxuICAgICAgICAvLyAgIGF3YWl0ICh3aW5kb3cgYXMgYW55KS5nYygpO1xyXG5cclxuICAgICAgICBpZiAodGVzdFJ1bikge1xyXG4gICAgICAgICAgcmVzLnB1c2goeyAuLi50ZXN0UnVuLCB3aWRnZXRzRGlmZmVyZW5jZTogZ2V0V2lkZ2V0c0NvdW50U2FmZSgpIC0gd2lkZ2V0c0JlZm9yZSB9KTtcclxuICAgICAgICAgIC8vIFJldHVybiBlYXJseSBpZiByZXR1cm5PbkZhaWwgaXMgc2V0IGFuZCB0ZXN0IGZhaWxlZCAoYnV0IGlnbm9yZSBmYWlsdXJlIGZvciB0aGUgc2tpcFRvVGVzdCB0ZXN0IGl0c2VsZilcclxuICAgICAgICAgIGlmIChvcHRpb25zLnJldHVybk9uRmFpbCAmJiBvcHRpb25zLnNraXBUb1Rlc3QgIT09IHRlc3QubmFtZSAmJiAhdGVzdFJ1bi5zdWNjZXNzICYmICF0ZXN0UnVuLnNraXBwZWQpXHJcbiAgICAgICAgICAgIHJldHVybiByZXM7XHJcbiAgICAgICAgfVxyXG4gICAgICAgIC8vIHJlcy5wdXNoKHsgLi4udGVzdFJ1biwgbWVtb3J5RGVsdGE6ICh3aW5kb3c/LnBlcmZvcm1hbmNlIGFzIGFueSk/Lm1lbW9yeT8udXNlZEpTSGVhcFNpemUgLSBtZW1vcnlVc2FnZUJlZm9yZSwgd2lkZ2V0c0RpZmZlcmVuY2U6IGdldFdpZGdldHNDb3VudFNhZmUoKSAtIHdpZGdldHNCZWZvcmUgfSk7XHJcblxyXG4gICAgICB9XHJcbiAgICB9XHJcbiAgICByZXR1cm4gcmVzO1xyXG4gIH1cclxuXHJcbiAgZnVuY3Rpb24gZ2V0V2lkZ2V0c0NvdW50U2FmZSgpIHtcclxuICAgIGlmICh0eXBlb2YgcHJvY2VzcyAhPT0gJ3VuZGVmaW5lZCcpXHJcbiAgICAgIHJldHVybiAwO1xyXG4gICAgbGV0IGxlbmd0aCA9IC0xO1xyXG4gICAgdHJ5IHtcclxuICAgICAgbGVuZ3RoID0gREcuV2lkZ2V0LmdldEFsbCgpLmxlbmd0aDtcclxuICAgIH0gY2F0Y2ggKGU6IGFueSkge1xyXG4gICAgICBjb25zb2xlLndhcm4oZS5tZXNzYWdlID8/IGUpO1xyXG4gICAgfVxyXG4gICAgcmV0dXJuIGxlbmd0aDtcclxuICB9XHJcblxyXG4gIGFzeW5jIGZ1bmN0aW9uIGludm9rZVRlc3RzKGNhdGVnb3JpZXNUb0ludm9rZTogeyBba2V5OiBzdHJpbmddOiBDYXRlZ29yeSB9LCBvcHRpb25zOiBUZXN0RXhlY3V0aW9uT3B0aW9ucykge1xyXG4gICAgdHJ5IHtcclxuICAgICAgbGV0IHNraXBwaW5nQ2F0ZWdvcmllcyA9IG9wdGlvbnM/LnNraXBUb0NhdGVnb3J5ICE9IHVuZGVmaW5lZDtcclxuICAgICAgbGV0IGlzVGFyZ2V0Q2F0ZWdvcnkgPSBmYWxzZTtcclxuICAgICAgZm9yIChjb25zdCBba2V5LCB2YWx1ZV0gb2YgT2JqZWN0LmVudHJpZXMoY2F0ZWdvcmllc1RvSW52b2tlKSkge1xyXG4gICAgICAgICAgaWYgKG9wdGlvbnMuZXhjbHVkZT8uc29tZSgoYykgPT4ga2V5LnN0YXJ0c1dpdGgoYykpKVxyXG4gICAgICAgICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICAgICAgaWYgKG9wdGlvbnM/LmNhdGVnb3J5ICE9IG51bGwgJiYgIWtleS50b0xvd2VyQ2FzZSgpLnN0YXJ0c1dpdGgob3B0aW9ucz8uY2F0ZWdvcnkudG9Mb3dlckNhc2UoKS50cmltKCkpKVxyXG4gICAgICAgICAgICAgIGNvbnRpbnVlO1xyXG5cclxuICAgICAgICAgIGlmIChza2lwcGluZ0NhdGVnb3JpZXMpIHtcclxuICAgICAgICAgICAgICBpZiAoaXNUYXJnZXRDYXRlZ29yeSlcclxuICAgICAgICAgICAgICAgICAgc2tpcHBpbmdDYXRlZ29yaWVzID0gZmFsc2U7XHJcbiAgICAgICAgICAgICAgZWxzZSB7XHJcbiAgICAgICAgICAgICAgICAgIGlmIChvcHRpb25zPy5za2lwVG9DYXRlZ29yeSAhPSBudWxsICYmIGtleS50b0xvd2VyQ2FzZSgpLnRyaW0oKSA9PT0gb3B0aW9ucz8uc2tpcFRvQ2F0ZWdvcnkudG9Mb3dlckNhc2UoKS50cmltKCkpIHtcclxuICAgICAgICAgICAgICAgICAgICAgIGlzVGFyZ2V0Q2F0ZWdvcnkgPSB0cnVlO1xyXG4gICAgICAgICAgICAgICAgICB9IGVsc2Uge1xyXG4gICAgICAgICAgICAgICAgICAgICAgLy8gSGF2ZW4ndCBmb3VuZCB0aGUgdGFyZ2V0IGNhdGVnb3J5IHlldCwga2VlcCBza2lwcGluZ1xyXG4gICAgICAgICAgICAgICAgICAgICAgY29udGludWU7XHJcbiAgICAgICAgICAgICAgICAgIH1cclxuICAgICAgICAgICAgICB9XHJcbiAgICAgICAgICB9XHJcbiAgICAgICAgICAvL0B0cy1pZ25vcmVcclxuICAgICAgICAgIGNvbnN0IHNraXBwZWQgPSB2YWx1ZS50ZXN0cz8uZXZlcnkoKHQ6IFRlc3QpID0+IHQub3B0aW9ucz8uc2tpcFJlYXNvblxyXG4gICAgICAgICAgICAgIHx8IChvcHRpb25zPy50ZXN0ICE9IG51bGwgJiYgb3B0aW9ucy50ZXN0LnRvTG93ZXJDYXNlKCkgIT09IHQubmFtZS50b0xvd2VyQ2FzZSgpKSk7XHJcblxyXG4gICAgICAgICAgaWYgKCFza2lwcGVkKSB7XHJcbiAgICAgICAgICAgICAgLy9AdHMtaWdub3JlXHJcbiAgICAgICAgICAgICAgY29uc3Qgc2tpcHBlZENvdW50ID0gKHZhbHVlLnRlc3RzID8/IFtdKS5maWx0ZXIoKHQ6IFRlc3QpID0+XHJcbiAgICAgICAgICAgICAgICB0Lm9wdGlvbnM/LnNraXBSZWFzb24gfHwgKG9wdGlvbnM/LnRlc3QgIT0gbnVsbCAmJiBvcHRpb25zLnRlc3QudG9Mb3dlckNhc2UoKSAhPT0gdC5uYW1lLnRvTG93ZXJDYXNlKCkpXHJcbiAgICAgICAgICAgICAgKS5sZW5ndGg7XHJcbiAgICAgICAgICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFN0YXJ0ZWQge3ske2tleX19fSR7c2tpcHBlZENvdW50ID4gMCA/IGAgc2tpcHBlZCB7eyR7c2tpcHBlZENvdW50fX19YCA6ICcnfWApO1xyXG4gICAgICAgICAgICAgIHZhbHVlLmJlZm9yZVN0YXR1cyA9IGF3YWl0IGludm9rZUNhdGVnb3J5TWV0aG9kKHZhbHVlLmJlZm9yZSwga2V5KTtcclxuICAgICAgICAgIH1cclxuICAgICAgICAgIGxldCB0ID0gdmFsdWUudGVzdHMgPz8gW107XHJcblxyXG4gICAgICAgICAgaWYgKG9wdGlvbnMuc3RyZXNzVGVzdCkge1xyXG4gICAgICAgICAgICAgIHQgPSB0LmZpbHRlcigoZSkgPT4gZS5vcHRpb25zPy5zdHJlc3NUZXN0KTtcclxuICAgICAgICAgICAgICB0ID0gc2h1ZmZsZSh0KTtcclxuICAgICAgICAgIH1cclxuXHJcbiAgICAgICAgICBpZiAoKG9wdGlvbnMudGFncz8ubGVuZ3RoID8/IDApID4gMCkge1xyXG4gICAgICAgICAgICAgIHQgPSB0LmZpbHRlcigoZSkgPT5cclxuICAgICAgICAgICAgICAgICAgZS5vcHRpb25zPy50YWdzPy5zb21lKHRhZyA9PiAob3B0aW9ucz8udGFncyA/PyBbXSkuaW5jbHVkZXModGFnKSlcclxuICAgICAgICAgICAgICApO1xyXG4gICAgICAgICAgfVxyXG5cclxuICAgICAgICAgIGxldCByZXM6IFRlc3RSZXN1bHRFeHRlbmRlZFtdO1xyXG4gICAgICAgICAgaWYgKHZhbHVlLmJlZm9yZVN0YXR1cykge1xyXG4gICAgICAgICAgICAgIHJlcyA9IEFycmF5LmZyb20odC5tYXAoKHRlc3RFbGVtKSA9PiB7XHJcbiAgICAgICAgICAgICAgICAgIHJldHVybiB7XHJcbiAgICAgICAgICAgICAgICAgICAgICBkYXRlOiBuZXcgRGF0ZSgpLnRvSVNPU3RyaW5nKCksXHJcbiAgICAgICAgICAgICAgICAgICAgICBjYXRlZ29yeToga2V5LFxyXG4gICAgICAgICAgICAgICAgICAgICAgbmFtZTogdGVzdEVsZW0ubmFtZSxcclxuICAgICAgICAgICAgICAgICAgICAgIHN1Y2Nlc3M6IGZhbHNlLFxyXG4gICAgICAgICAgICAgICAgICAgICAgcmVzdWx0OiAnYmVmb3JlKCkgZmFpbGVkJyxcclxuICAgICAgICAgICAgICAgICAgICAgIG1zOiAwLFxyXG4gICAgICAgICAgICAgICAgICAgICAgc2tpcHBlZDogZmFsc2UsXHJcbiAgICAgICAgICAgICAgICAgICAgICBsb2dzOiAnJyxcclxuICAgICAgICAgICAgICAgICAgICAgIG93bmVyOiBwYWNrYWdlT3duZXIsXHJcbiAgICAgICAgICAgICAgICAgICAgICBwYWNrYWdlOiBwYWNrYWdlXy5uYW1lLFxyXG4gICAgICAgICAgICAgICAgICAgICAgd2lkZ2V0c0RpZmZlcmVuY2U6IDAsXHJcbiAgICAgICAgICAgICAgICAgICAgICBmbGFraW5nOiBERy5UZXN0LmlzUmVwcm9kdWNpbmdcclxuICAgICAgICAgICAgICAgICAgfTtcclxuICAgICAgICAgICAgICB9KSk7XHJcbiAgICAgICAgICAgICAgcmVzLmZvckVhY2goYXN5bmMgKHRlc3QpID0+IGF3YWl0IGdyb2suc2hlbGwucmVwb3J0VGVzdCgncGFja2FnZScsIHRlc3QpKTtcclxuICAgICAgICAgIH0gZWxzZVxyXG4gICAgICAgICAgICAgIHJlcyA9IGF3YWl0IGludm9rZVRlc3RzSW5DYXRlZ29yeSh2YWx1ZSwgb3B0aW9ucywgc2tpcHBpbmdDYXRlZ29yaWVzKTtcclxuICAgICAgICAgIGNvbnN0IGRhdGE6IFRlc3RSZXN1bHRFeHRlbmRlZFtdID0gcmVzLmZpbHRlcigoZCkgPT4gZC5yZXN1bHQgIT0gJ3NraXBwZWQnKTtcclxuXHJcbiAgICAgICAgICBpZiAoIXNraXBwZWQpXHJcbiAgICAgICAgICAgICAgdmFsdWUuYWZ0ZXJTdGF0dXMgPSBhd2FpdCBpbnZva2VDYXRlZ29yeU1ldGhvZCh2YWx1ZS5hZnRlciwga2V5KTtcclxuXHJcbiAgICAgICAgICAvLyBDbGVhciBhZnRlciBjYXRlZ29yeVxyXG4gICAgICAgICAgLy8gZ3Jvay5zaGVsbC5jbG9zZUFsbCgpO1xyXG4gICAgICAgICAgLy8gREcuQmFsbG9vbi5jbG9zZUFsbCgpO1xyXG4gICAgICAgICAgaWYgKHZhbHVlLmFmdGVyU3RhdHVzKSB7XHJcbiAgICAgICAgICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IENhdGVnb3J5IGFmdGVyKCkge3ske2tleX19fSBmYWlsZWRgKTtcclxuICAgICAgICAgICAgICBzdGRMb2coYFBhY2thZ2UgdGVzdGluZzogUmVzdWx0IGZvciB7eyR7a2V5fX19IGFmdGVyOiAke3ZhbHVlLmFmdGVyU3RhdHVzfWApO1xyXG4gICAgICAgICAgICAgIGRhdGEucHVzaCh7XHJcbiAgICAgICAgICAgICAgICAgIGRhdGU6IG5ldyBEYXRlKCkudG9JU09TdHJpbmcoKSxcclxuICAgICAgICAgICAgICAgICAgY2F0ZWdvcnk6IGtleSxcclxuICAgICAgICAgICAgICAgICAgbmFtZTogJ2FmdGVyJyxcclxuICAgICAgICAgICAgICAgICAgc3VjY2VzczogZmFsc2UsXHJcbiAgICAgICAgICAgICAgICAgIHJlc3VsdDogdmFsdWUuYWZ0ZXJTdGF0dXMsXHJcbiAgICAgICAgICAgICAgICAgIG1zOiAwLFxyXG4gICAgICAgICAgICAgICAgICBza2lwcGVkOiBmYWxzZSxcclxuICAgICAgICAgICAgICAgICAgbG9nczogJycsXHJcbiAgICAgICAgICAgICAgICAgIG93bmVyOiBwYWNrYWdlT3duZXIsXHJcbiAgICAgICAgICAgICAgICAgIHBhY2thZ2U6IHBhY2thZ2VfLm5hbWUsXHJcbiAgICAgICAgICAgICAgICAgIHdpZGdldHNEaWZmZXJlbmNlOiAwLFxyXG4gICAgICAgICAgICAgICAgICBmbGFraW5nOiBERy5UZXN0LmlzUmVwcm9kdWNpbmdcclxuICAgICAgICAgICAgICB9KTtcclxuICAgICAgICAgIH1cclxuICAgICAgICAgIGlmICh2YWx1ZS5iZWZvcmVTdGF0dXMpIHtcclxuICAgICAgICAgICAgICBzdGRMb2coYFBhY2thZ2UgdGVzdGluZzogQ2F0ZWdvcnkgYmVmb3JlKCkge3ske2tleX19fSBmYWlsZWRgKTtcclxuICAgICAgICAgICAgICBzdGRMb2coYFBhY2thZ2UgdGVzdGluZzogUmVzdWx0IGZvciB7eyR7a2V5fX19IGJlZm9yZTogJHt2YWx1ZS5iZWZvcmVTdGF0dXN9YCk7XHJcbiAgICAgICAgICAgICAgZGF0YS5wdXNoKHtcclxuICAgICAgICAgICAgICAgICAgZGF0ZTogbmV3IERhdGUoKS50b0lTT1N0cmluZygpLFxyXG4gICAgICAgICAgICAgICAgICBjYXRlZ29yeToga2V5LFxyXG4gICAgICAgICAgICAgICAgICBuYW1lOiAnYmVmb3JlJyxcclxuICAgICAgICAgICAgICAgICAgc3VjY2VzczogZmFsc2UsXHJcbiAgICAgICAgICAgICAgICAgIHJlc3VsdDogdmFsdWUuYmVmb3JlU3RhdHVzLFxyXG4gICAgICAgICAgICAgICAgICBtczogMCxcclxuICAgICAgICAgICAgICAgICAgc2tpcHBlZDogZmFsc2UsXHJcbiAgICAgICAgICAgICAgICAgIGxvZ3M6ICcnLFxyXG4gICAgICAgICAgICAgICAgICBvd25lcjogcGFja2FnZU93bmVyLFxyXG4gICAgICAgICAgICAgICAgICBwYWNrYWdlOiBwYWNrYWdlXy5uYW1lLFxyXG4gICAgICAgICAgICAgICAgICB3aWRnZXRzRGlmZmVyZW5jZTogMCxcclxuICAgICAgICAgICAgICAgICAgZmxha2luZzogREcuVGVzdC5pc1JlcHJvZHVjaW5nXHJcbiAgICAgICAgICAgICAgfSk7XHJcbiAgICAgICAgICB9XHJcbiAgICAgICAgICByZXN1bHRzLnB1c2goLi4uZGF0YSk7XHJcblxyXG4gICAgICAgICAgLy8gSWYgcmV0dXJuT25GYWlsIGlzIHNldCBhbmQgYSB0ZXN0IGZhaWxlZCAob3RoZXIgdGhhbiBza2lwVG9UZXN0KSwgc3RvcCBwcm9jZXNzaW5nIG1vcmUgY2F0ZWdvcmllc1xyXG4gICAgICAgICAgaWYgKG9wdGlvbnMucmV0dXJuT25GYWlsICYmIGRhdGEuc29tZSgoZCkgPT4gIWQuc3VjY2VzcyAmJiAhZC5za2lwcGVkICYmIGQubmFtZSAhPT0gb3B0aW9ucy5za2lwVG9UZXN0KSlcclxuICAgICAgICAgICAgICBicmVhaztcclxuICAgICAgfVxyXG4gICAgfSBmaW5hbGx5IHtcclxuICAgICAgcmVzZXRDb25zb2xlKCk7XHJcbiAgICB9XHJcbiAgICBpZiAob3B0aW9ucy50ZXN0Q29udGV4dCEuY2F0Y2hVbmhhbmRsZWQgJiYgKCFERy5UZXN0LmlzSW5CZW5jaG1hcmspKSB7XHJcbiAgICAgIGF3YWl0IGRlbGF5KDEwMDApO1xyXG4gICAgICBjb25zdCBlcnJvciA9IGF3YWl0IGdyb2suc2hlbGwubGFzdEVycm9yO1xyXG4gICAgICBpZiAoZXJyb3IgIT0gdW5kZWZpbmVkKSB7XHJcbiAgICAgICAgICBjb25zdCBwYXJhbXM6IGFueSA9IHtcclxuICAgICAgICAgICAgICBsb2dzOiAnJyxcclxuICAgICAgICAgICAgICBkYXRlOiBuZXcgRGF0ZSgpLnRvSVNPU3RyaW5nKCksXHJcbiAgICAgICAgICAgICAgY2F0ZWdvcnk6ICdVbmhhbmRsZWQgZXhjZXB0aW9ucycsXHJcbiAgICAgICAgICAgICAgbmFtZTogJ0V4Y2VwdGlvbicsXHJcbiAgICAgICAgICAgICAgcmVzdWx0OiBlcnJvciA/PyAnJyxcclxuICAgICAgICAgICAgICBzdWNjZXNzOiAhZXJyb3IsXHJcbiAgICAgICAgICAgICAgbXM6IDAsXHJcbiAgICAgICAgICAgICAgc2tpcHBlZDogZmFsc2UsXHJcbiAgICAgICAgICAgICAgb3duZXI6IHBhY2thZ2VPd25lciA/PyAnJyxcclxuICAgICAgICAgICAgICAncGFja2FnZSc6IHBhY2thZ2VfLm5hbWUsXHJcbiAgICAgICAgICAgICAgd2lkZ2V0c0RpZmZlcmVuY2U6IDBcclxuICAgICAgICAgIH07XHJcbiAgICAgICAgICBzdGRMb2coYFBhY2thZ2UgdGVzdGluZzogVW5oYW5kbGVkIEV4Y2VwdGlvbjogJHtlcnJvcn1gKTtcclxuXHJcbiAgICAgICAgICByZXN1bHRzLnB1c2goey4uLnBhcmFtcywgJ2ZsYWtpbmcnOiBERy5UZXN0LmlzUmVwcm9kdWNpbmcgJiYgIWVycm9yfSk7XHJcbiAgICAgICAgICAoPGFueT5wYXJhbXMpLnBhY2thZ2UgPSBwYWNrYWdlXy5uYW1lO1xyXG4gICAgICAgICAgYXdhaXQgZ3Jvay5zaGVsbC5yZXBvcnRUZXN0KCdwYWNrYWdlJywgcGFyYW1zKTtcclxuICAgICAgfVxyXG4gICAgfVxyXG4gIH1cclxufVxyXG5cclxuYXN5bmMgZnVuY3Rpb24gZ2V0UmVzdWx0KHg6IGFueSk6IFByb21pc2U8c3RyaW5nPiB7XHJcbiAgcmV0dXJuIGAke3gudG9TdHJpbmcoKX1cXG4ke3guc3RhY2sgPyAoYXdhaXQgREcuTG9nZ2VyLnRyYW5zbGF0ZVN0YWNrVHJhY2UoeC5zdGFjaykpIDogJyd9YDtcclxufVxyXG5cclxuYXN5bmMgZnVuY3Rpb24gZXhlY1Rlc3QodDogVGVzdCwgcHJlZGljYXRlOiBzdHJpbmcgfCB1bmRlZmluZWQsIGxvZ3M6IGFueVtdLFxyXG4gIHRlc3RUaW1lb3V0PzogbnVtYmVyLCBwYWNrYWdlTmFtZT86IHN0cmluZywgdmVyYm9zZT86IGJvb2xlYW5cclxuKTogUHJvbWlzZTxUZXN0UmVzdWx0IHwgdW5kZWZpbmVkPiB7XHJcbiAgbG9ncy5sZW5ndGggPSAwO1xyXG4gIGxldCByOiBUZXN0UmVzdWx0O1xyXG4gIGxldCB0eXBlOiBzdHJpbmcgPSAncGFja2FnZSc7XHJcbiAgY29uc3QgZmlsdGVyID0gcHJlZGljYXRlICE9IHVuZGVmaW5lZCAmJiAodC5uYW1lLnRvTG93ZXJDYXNlKCkgIT09IHByZWRpY2F0ZS50b0xvd2VyQ2FzZSgpKTtcclxuICBsZXQgc2tpcCA9IHQub3B0aW9ucz8uc2tpcFJlYXNvbiB8fCBmaWx0ZXI7XHJcbiAgbGV0IHNraXBSZWFzb24gPSBmaWx0ZXIgPyAnc2tpcHBlZCcgOiB0Lm9wdGlvbnM/LnNraXBSZWFzb247XHJcblxyXG4gIGlmIChERy5UZXN0LmlzSW5CZW5jaG1hcmsgJiYgIXQub3B0aW9ucz8uYmVuY2htYXJrKSB7XHJcbiAgICBzdGRMb2coYFBhY2thZ2UgdGVzdGluZzogU2tpcHBlZCB7eyR7dC5jYXRlZ29yeX19fSB7eyR7dC5uYW1lfX19IGRvZXNudCBhdmFpbGFibGUgaW4gYmVuY2htYXJrIG1vZGVgKTtcclxuICAgIHJldHVybiB1bmRlZmluZWQ7XHJcbiAgfVxyXG5cclxuICBpZiAoc2tpcCAmJiAhREcuVGVzdC5pc0luQmVuY2htYXJrKVxyXG4gICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFNraXBwZWQge3ske3QuY2F0ZWdvcnl9fX0ge3ske3QubmFtZX19fWApO1xyXG4gIGlmICghc2tpcClcclxuICAgIHN0ZExvZyhgUGFja2FnZSB0ZXN0aW5nOiBTdGFydGVkIHt7JHt0LmNhdGVnb3J5fX19IHt7JHt0Lm5hbWV9fX1gKTtcclxuICBjb25zdCBzdGFydCA9IERhdGUubm93KCk7XHJcbiAgY29uc3Qgc3RhcnREYXRlID0gbmV3IERhdGUoc3RhcnQpLnRvSVNPU3RyaW5nKCk7XHJcbiAgdHJ5IHtcclxuICAgIGlmIChza2lwKVxyXG4gICAgICByID0geyBuYW1lOiB0Lm5hbWUsIG93bmVyOnQub3B0aW9ucz8ub3duZXIgPz8gJycsIGNhdGVnb3J5OiB0LmNhdGVnb3J5LCBsb2dzOiAnJywgZGF0ZTogc3RhcnREYXRlLCBzdWNjZXNzOiB0cnVlLCByZXN1bHQ6IHNraXBSZWFzb24hLCBtczogMCwgc2tpcHBlZDogdHJ1ZSwgcGFja2FnZTogcGFja2FnZU5hbWUgPz8gJycsIGZsYWtpbmc6IERHLlRlc3QuaXNSZXByb2R1Y2luZ307XHJcbiAgICBlbHNlIHtcclxuICAgICAgbGV0IHRpbWVvdXRfID0gdGVzdFRpbWVvdXQgPz8gU1RBTkRBUlRfVElNRU9VVDtcclxuXHJcbiAgICAgIGlmIChERy5UZXN0LmlzUHJvZmlsaW5nKVxyXG4gICAgICAgIGNvbnNvbGUucHJvZmlsZShgJHt0LmNhdGVnb3J5fTogJHt0Lm5hbWV9YCk7XHJcblxyXG4gICAgICByID0geyBuYW1lOiB0Lm5hbWUsIG93bmVyOnQub3B0aW9ucz8ub3duZXIgPz8gJycsIGNhdGVnb3J5OiB0LmNhdGVnb3J5LCBsb2dzOiAnJywgZGF0ZTogc3RhcnREYXRlLCBzdWNjZXNzOiB0cnVlLCByZXN1bHQ6IChhd2FpdCB0aW1lb3V0KHQudGVzdCwgdGltZW91dF8pKS50b1N0cmluZygpID8/ICdPSycsIG1zOiAwLCBza2lwcGVkOiBmYWxzZSAsIHBhY2thZ2U6IHBhY2thZ2VOYW1lID8/ICcnLCBmbGFraW5nOiBERy5UZXN0LmlzUmVwcm9kdWNpbmd9O1xyXG5cclxuICAgICAgaWYgKERHLlRlc3QuaXNQcm9maWxpbmcpIHtcclxuICAgICAgICBjb25zb2xlLnByb2ZpbGVFbmQoYCR7dC5jYXRlZ29yeX06ICR7dC5uYW1lfWApO1xyXG4gICAgICAgIGdyb2suc2hlbGwuaW5mbyhgUHJvZmlsaW5nIG9mICR7dC5jYXRlZ29yeX06ICR7dC5uYW1lfSBmaW5pc2hlZCBcXG4gUGxlYXNlIGVuc3VyZSB0aGF0IHlvdSBoYXZlIG9wZW5lZCBEZXZUb29scyAoRjEyKSAvIFBlcmZvcm1hbmNlIHBhbmVsIGJlZm9yZSB0ZXN0IHN0YXJ0cy5gKTtcclxuICAgICAgfVxyXG4gICAgfVxyXG4gIH0gY2F0Y2ggKHg6IGFueSkge1xyXG4gICAgc3RkRXJyb3IoeCk7XHJcbiAgICByID0geyBuYW1lOiB0Lm5hbWUsIG93bmVyOnQub3B0aW9ucz8ub3duZXIgPz8gJycsIGNhdGVnb3J5OiB0LmNhdGVnb3J5LCBsb2dzOiAnJywgZGF0ZTogc3RhcnREYXRlLCBzdWNjZXNzOiBmYWxzZSwgcmVzdWx0OiBhd2FpdCBnZXRSZXN1bHQoeCksIG1zOiAwLCBza2lwcGVkOiBmYWxzZSwgcGFja2FnZTogcGFja2FnZU5hbWUgPz8gJycsIGZsYWtpbmc6IGZhbHNlfTtcclxuICB9XHJcbiAgaWYgKHQub3B0aW9ucz8uaXNBZ2dyZWdhdGVkICYmIHIucmVzdWx0LmNvbnN0cnVjdG9yID09PSBERy5EYXRhRnJhbWUpIHtcclxuICAgIGNvbnN0IGNvbCA9IHIucmVzdWx0LmNvbCgnc3VjY2VzcycpO1xyXG4gICAgaWYgKGNvbClcclxuICAgICAgci5zdWNjZXNzID0gY29sLnN0YXRzLnN1bSA9PT0gY29sLmxlbmd0aDtcclxuICAgIGlmICghdmVyYm9zZSkge1xyXG4gICAgICBjb25zdCBkZiA9IHIucmVzdWx0O1xyXG4gICAgICBkZi5jb2x1bW5zLnJlbW92ZSgnc3RhY2snKTtcclxuICAgICAgZGYucm93cy5yZW1vdmVXaGVyZSgocikgPT4gci5nZXQoJ3N1Y2Nlc3MnKSk7XHJcbiAgICAgIHIucmVzdWx0ID0gZGY7XHJcbiAgICB9XHJcbiAgICByLnJlc3VsdCA9IHIucmVzdWx0LnRvQ3N2KCk7XHJcbiAgfVxyXG4gIHIubG9ncyA9IGxvZ3Muam9pbignXFxuJyk7XHJcbiAgci5tcyA9IERhdGUubm93KCkgLSBzdGFydDtcclxuICBpZiAoIXNraXApXHJcbiAgICBzdGRMb2coYFBhY2thZ2UgdGVzdGluZzogRmluaXNoZWQge3ske3QuY2F0ZWdvcnl9fX0ge3ske3QubmFtZX19fSB3aXRoIHt7JHtyLnN1Y2Nlc3MgPyAnc3VjY2VzcycgOiAnZXJyb3InfX19IGZvciAke3IubXN9IG1zYCk7XHJcbiAgaWYgKCFyLnN1Y2Nlc3MpIHtcclxuICAgICAgc3RkTG9nKGBQYWNrYWdlIHRlc3Rpbmc6IFJlc3VsdCBmb3Ige3ske3QuY2F0ZWdvcnl9fX0ge3ske3QubmFtZX19fTogJHtyLnJlc3VsdH1gKTtcclxuICB9XHJcbiAgci5jYXRlZ29yeSA9IHQuY2F0ZWdvcnk7XHJcbiAgci5uYW1lID0gdC5uYW1lO1xyXG4gIHIub3duZXIgPSB0Lm9wdGlvbnM/Lm93bmVyID8/ICcnO1xyXG4gIGlmICghZmlsdGVyKSB7XHJcbiAgICBsZXQgcGFyYW1zID0ge1xyXG4gICAgICAnc3VjY2Vzcyc6IHIuc3VjY2VzcywgJ3Jlc3VsdCc6IHIucmVzdWx0LCAnbXMnOiByLm1zLCAnZGF0ZSc6IHIuZGF0ZSxcclxuICAgICAgJ3NraXBwZWQnOiByLnNraXBwZWQsICdjYXRlZ29yeSc6IHQuY2F0ZWdvcnksICduYW1lJzogdC5uYW1lLCAnbG9ncyc6IHIubG9ncywgJ293bmVyJzogci5vd25lcixcclxuICAgICAgJ2ZsYWtpbmcnOiBERy5UZXN0LmlzUmVwcm9kdWNpbmcgJiYgci5zdWNjZXNzLFxyXG4gICAgICAncGFja2FnZSc6IHIucGFja2FnZVxyXG4gICAgfTtcclxuICAgIGlmIChyLnJlc3VsdC5jb25zdHJ1Y3RvciA9PSBPYmplY3QpIHtcclxuICAgICAgY29uc3QgcmVzID0gT2JqZWN0LmtleXMoci5yZXN1bHQpLnJlZHVjZSgoYWNjLCBrKSA9PiAoeyAuLi5hY2MsIFsncmVzdWx0LicgKyBrXTogci5yZXN1bHRba10gfSksIHt9KTtcclxuICAgICAgcGFyYW1zID0geyAuLi5wYXJhbXMsIC4uLnJlcyB9O1xyXG4gICAgfVxyXG5cclxuICAgIGlmIChwYXJhbXMucmVzdWx0IGluc3RhbmNlb2YgREcuRGF0YUZyYW1lKVxyXG4gICAgICBwYXJhbXMucmVzdWx0ID0gSlNPTi5zdHJpbmdpZnkocGFyYW1zLnJlc3VsdD8udG9Kc29uKCkpIHx8ICcnO1xyXG4gICAgYXdhaXQgZ3Jvay5zaGVsbC5yZXBvcnRUZXN0KHR5cGUsIHBhcmFtcyk7XHJcbiAgfVxyXG4gIHJldHVybiByO1xyXG59XHJcblxyXG5leHBvcnQgZnVuY3Rpb24gc2h1ZmZsZShhcnJheTogYW55W10pOiBhbnlbXSB7XHJcbiAgY29uc3QgbmV3QXJyID0gYXJyYXkuc2xpY2UoKTtcclxuICBuZXdBcnIuc29ydCgoKSA9PiBNYXRoLnJhbmRvbSgpIC0gMC41KTtcclxuICByZXR1cm4gbmV3QXJyO1xyXG59XHJcblxyXG4vKiBXYWl0cyBbbXNdIG1pbGxpc2Vjb25kcyAqL1xyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gZGVsYXkobXM6IG51bWJlcikge1xyXG4gIGF3YWl0IG5ldyBQcm9taXNlKChyKSA9PiBzZXRUaW1lb3V0KHIsIG1zKSk7XHJcbn1cclxuXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBhd2FpdENoZWNrKGNoZWNrSGFuZGxlcjogKCkgPT4gYm9vbGVhbixcclxuICBlcnJvcjogc3RyaW5nID0gJ1RpbWVvdXQgZXhjZWVkZWQnLCB3YWl0OiBudW1iZXIgPSA1MDAsIGludGVydmFsOiBudW1iZXIgPSA1MCk6IFByb21pc2U8YW55PiB7XHJcbiAgcmV0dXJuIG5ldyBQcm9taXNlKChyZXNvbHZlLCByZWplY3QpID0+IHtcclxuICAgIHNldFRpbWVvdXQoKCkgPT4ge1xyXG4gICAgICBjbGVhckludGVydmFsKGludGVydmFsSWQpO1xyXG4gICAgICByZWplY3QobmV3IEVycm9yKGVycm9yKSk7XHJcbiAgICB9LCB3YWl0KTtcclxuICAgIC8vIEB0cy1pZ25vcmVcclxuICAgIGNvbnN0IGludGVydmFsSWQ6IFRpbWVvdXQgPSBzZXRJbnRlcnZhbCgoKSA9PiB7XHJcbiAgICAgIGlmIChjaGVja0hhbmRsZXIoKSkge1xyXG4gICAgICAgIGNsZWFySW50ZXJ2YWwoaW50ZXJ2YWxJZCk7XHJcbiAgICAgICAgcmVzb2x2ZShudWxsKTtcclxuICAgICAgfVxyXG4gICAgfSwgaW50ZXJ2YWwpO1xyXG4gIH0pO1xyXG59XHJcblxyXG4vLyBSZXR1cm5zIHRlc3QgZXhlY3V0aW9uIHJlc3VsdCBvciBhbiBlcnJvciBpbiBjYXNlIG9mIHRpbWVvdXRcclxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHRpbWVvdXQoZnVuYzogKCkgPT4gUHJvbWlzZTxhbnk+LCB0ZXN0VGltZW91dDogbnVtYmVyLCB0aW1lb3V0UmVhc29uOiBzdHJpbmcgPSAnRVhFQ1VUSU9OIFRJTUVPVVQnKTogUHJvbWlzZTxhbnk+IHtcclxuICBsZXQgdGltZW91dDogYW55ID0gbnVsbDtcclxuICBjb25zdCB0aW1lb3V0UHJvbWlzZSA9IG5ldyBQcm9taXNlPGFueT4oKF8sIHJlamVjdCkgPT4ge1xyXG4gICAgdGltZW91dCA9IHNldFRpbWVvdXQoKCkgPT4ge1xyXG4gICAgICAvLyBlc2xpbnQtZGlzYWJsZS1uZXh0LWxpbmUgcHJlZmVyLXByb21pc2UtcmVqZWN0LWVycm9yc1xyXG4gICAgICByZWplY3QodGltZW91dFJlYXNvbik7XHJcbiAgICB9LCB0ZXN0VGltZW91dCk7XHJcbiAgfSk7XHJcbiAgdHJ5IHtcclxuICAgIHJldHVybiBhd2FpdCBQcm9taXNlLnJhY2UoW2Z1bmMoKSwgdGltZW91dFByb21pc2VdKTtcclxuICB9IGZpbmFsbHkge1xyXG4gICAgaWYgKHRpbWVvdXQpXHJcbiAgICAgIGNsZWFyVGltZW91dCh0aW1lb3V0KTtcclxuICB9XHJcbn1cclxuXHJcbmV4cG9ydCBmdW5jdGlvbiBpc0RpYWxvZ1ByZXNlbnQoZGlhbG9nVGl0bGU6IHN0cmluZyk6IGJvb2xlYW4ge1xyXG4gIGNvbnN0IGRpYWxvZ3MgPSBERy5EaWFsb2cuZ2V0T3BlbkRpYWxvZ3MoKTtcclxuICBmb3IgKGxldCBpID0gMDsgaSA8IGRpYWxvZ3MubGVuZ3RoOyBpKyspIHtcclxuICAgIGlmIChkaWFsb2dzW2ldLnRpdGxlID09IGRpYWxvZ1RpdGxlKVxyXG4gICAgICByZXR1cm4gdHJ1ZTtcclxuICB9XHJcbiAgcmV0dXJuIGZhbHNlO1xyXG59XHJcblxyXG4vKiogRXhwZWN0cyBhbiBhc3luY2hyb25vdXMge0BsaW5rIGFjdGlvbn0gdG8gdGhyb3cgYW4gZXhjZXB0aW9uLiBVc2Uge0BsaW5rIGNoZWNrfSB0byBwZXJmb3JtXHJcbiAqIGRlZXBlciBpbnNwZWN0aW9uIG9mIHRoZSBleGNlcHRpb24gaWYgbmVjZXNzYXJ5LlxyXG4gKiBAcGFyYW0gIHtmdW5jdGlvbigpOiBQcm9taXNlPHZvaWQ+fSBhY3Rpb25cclxuICogQHBhcmFtICB7ZnVuY3Rpb24oYW55KTogYm9vbGVhbn0gY2hlY2tcclxuICogQHJldHVybiB7UHJvbWlzZTx2b2lkPn1cclxuICovXHJcbmV4cG9ydCBhc3luYyBmdW5jdGlvbiBleHBlY3RFeGNlcHRpb25Bc3luYyhhY3Rpb246ICgpID0+IFByb21pc2U8dm9pZD4sXHJcbiAgY2hlY2s/OiAoZXhjZXB0aW9uOiBhbnkpID0+IGJvb2xlYW4pOiBQcm9taXNlPHZvaWQ+IHtcclxuICBsZXQgY2F1Z2h0OiBib29sZWFuID0gZmFsc2U7XHJcbiAgbGV0IGNoZWNrZWQ6IGJvb2xlYW4gPSBmYWxzZTtcclxuICB0cnkge1xyXG4gICAgYXdhaXQgYWN0aW9uKCk7XHJcbiAgfSBjYXRjaCAoZSkge1xyXG4gICAgY2F1Z2h0ID0gdHJ1ZTtcclxuICAgIGNoZWNrZWQgPSAhY2hlY2sgfHwgY2hlY2soZSk7XHJcbiAgfSBmaW5hbGx5IHtcclxuICAgIGlmICghY2F1Z2h0KVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoJ0FuIGV4Y2VwdGlvbiBpcyBleHBlY3RlZCBidXQgbm90IHRocm93bicpO1xyXG4gICAgaWYgKCFjaGVja2VkKVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoJ0FuIGV4cGVjdGVkIGV4Y2VwdGlvbiBpcyB0aHJvd24sIGJ1dCBpdCBkb2VzIG5vdCBzYXRpc2Z5IHRoZSBjb25kaXRpb24nKTtcclxuICB9XHJcbn1cclxuXHJcbmNvbnN0IGNhdERGID0gREcuRGF0YUZyYW1lLmZyb21Db2x1bW5zKFtERy5Db2x1bW4uZnJvbVN0cmluZ3MoJ2NvbCcsIFsndmFsMScsICd2YWwyJywgJ3ZhbDMnXSldKTtcclxuXHJcbi8qKlxyXG4gKiBVbml2ZXJzYWwgdGVzdCBmb3Igdmlld2Vycy4gSXQgc2VhcmNoIHZpZXdlcnMgaW4gRE9NIGJ5IHRhZ3M6IGNhbnZhcywgc3ZnLCBpbWcsIGlucHV0LCBoMSwgYVxyXG4gKiBAcGFyYW0gIHtzdHJpbmd9IHYgVmlld2VyIG5hbWVcclxuICogQHBhcmFtICB7X0RHLkRhdGFGcmFtZX0gZGYgRGF0YWZyYW1lIHRvIHVzZS4gU2hvdWxkIGhhdmUgYXQgbGVhc3QgMyByb3dzXHJcbiAqIEBwYXJhbSAge2Jvb2xlYW59IG9wdGlvbnMuZGV0ZWN0U2VtYW50aWNUeXBlcyBTcGVjaWZ5IHdoZXRoZXIgdG8gZGV0ZWN0IHNlbWFudGljIHR5cGVzIG9yIG5vdFxyXG4gKiBAcGFyYW0gIHtib29sZWFufSBvcHRpb25zLnJlYWRPbmx5IElmIHNldCB0byB0cnVlLCB0aGUgZGF0YWZyYW1lIHdpbGwgbm90IGJlIG1vZGlmaWVkIGR1cmluZyB0aGUgdGVzdFxyXG4gKiBAcGFyYW0gIHtib29sZWFufSBvcHRpb25zLmFyYml0cmFyeURmVGVzdCBJZiBzZXQgdG8gZmFsc2UsIHRlc3Qgb24gYXJiaXRyYXJ5IGRhdGFmcmFtZVxyXG4gKiAob25lIGNhdGVnb3JpY2FsIGNvbHVtbikgd2lsbCBub3QgYmUgcGVyZm9ybWVkXHJcbiAqIEBwYXJhbSAge29iamVjdH0gb3B0aW9ucyBMaXN0IG9mIG9wdGlvbnMgKG9wdGlvbmFsKVxyXG4gKiBAcmV0dXJuIHtQcm9taXNlPHZvaWQ+fSBUaGUgdGVzdCBpcyBjb25zaWRlcmVkIHN1Y2Nlc3NmdWwgaWYgaXQgY29tcGxldGVzIHdpdGhvdXQgZXJyb3JzXHJcbiAqL1xyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gdGVzdFZpZXdlcih2OiBzdHJpbmcsIGRmOiBfREcuRGF0YUZyYW1lLCBvcHRpb25zPzoge1xyXG4gIGRldGVjdFNlbWFudGljVHlwZXM/OiBib29sZWFuLCByZWFkT25seT86IGJvb2xlYW4sIGFyYml0cmFyeURmVGVzdD86IGJvb2xlYW4sXHJcbiAgcGFja2FnZU5hbWU/OiBzdHJpbmcsIGF3YWl0Vmlld2VyPzogKHZpZXdlcjogX0RHLlZpZXdlcikgPT4gUHJvbWlzZTx2b2lkPlxyXG59KTogUHJvbWlzZTx2b2lkPiB7XHJcbiAgY29uc3QgcGFja2FnZU5hbWUgPSBvcHRpb25zPy5wYWNrYWdlTmFtZSA/PyAnJztcclxuICBpZiAob3B0aW9ucz8uZGV0ZWN0U2VtYW50aWNUeXBlcylcclxuICAgIGF3YWl0IGdyb2suZGF0YS5kZXRlY3RTZW1hbnRpY1R5cGVzKGRmKTtcclxuICBjb25zdCB0diA9IGdyb2suc2hlbGwuYWRkVGFibGVWaWV3KGRmKTtcclxuXHJcbiAgdHJ5IHtcclxuICAgIC8vMS4gT3BlbiwgZG8gbm90aGluZyBhbmQgY2xvc2VcclxuICAgIGF3YWl0IHRlc3RWaWV3ZXJJbnRlcm5hbCh0diwgdiwgcGFja2FnZU5hbWUsIGdyb2suZXZlbnRzLm9uVmlld2VyQWRkZWQpO1xyXG4gICAgLy9pbiBjYXNlIHZpZXdlciB3aXRoIGFzeW5jIHJlbmRlcmluZyAtIHdhaXQgZm9yIHJlbmRlciB0byBjb21wbGV0ZVxyXG4gICAgaWYgKG9wdGlvbnM/LmF3YWl0Vmlld2VyKVxyXG4gICAgICBhd2FpdCB0ZXN0Vmlld2VySW50ZXJuYWwodHYsIHYsIHBhY2thZ2VOYW1lLCBncm9rLmV2ZW50cy5vblZpZXdlckFkZGVkLCB1bmRlZmluZWQsIG9wdGlvbnMhLmF3YWl0Vmlld2VyKTtcclxuXHJcbiAgICAvLzIuIE9wZW4gdmlld2VyLCBydW4gc2VsZWN0aW9uLCBmaWx0ZXIsIGV0Yy4gYW5kIGNsb3NlXHJcbiAgICBpZiAoIW9wdGlvbnM/LnJlYWRPbmx5KSB7XHJcbiAgICAgIGF3YWl0IHRlc3RWaWV3ZXJJbnRlcm5hbCh0diwgdiwgcGFja2FnZU5hbWUsIGdyb2suZXZlbnRzLm9uVmlld2VyQWRkZWQsIHNlbGVjdEZpbHRlckNoYW5nZUN1cnJlbnQpO1xyXG4gICAgICBpZiAob3B0aW9ucz8uYXdhaXRWaWV3ZXIpXHJcbiAgICAgICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCwgc2VsZWN0RmlsdGVyQ2hhbmdlQ3VycmVudCwgb3B0aW9ucyEuYXdhaXRWaWV3ZXIpO1xyXG4gICAgfVxyXG5cclxuICAgIC8vMi4gT3BlbiB2aWV3ZXIsIGNoYW5nZSBvcHRpb25zLCBzYXZlIGxheW91dCBhbmQgY2xvc2VcclxuICAgIGxldCBwcm9wc0FuZExheW91dDogeyBsYXlvdXQ6IGFueSwgc2F2ZWRQcm9wczogYW55IH0gfCBudWxsID0gbnVsbDtcclxuICAgIHByb3BzQW5kTGF5b3V0ID0gYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCwgY2hhbmdlT3B0aW9uc1NhdmVMYXlvdXQpO1xyXG4gICAgaWYgKG9wdGlvbnM/LmF3YWl0Vmlld2VyKVxyXG4gICAgICBwcm9wc0FuZExheW91dCA9IGF3YWl0IHRlc3RWaWV3ZXJJbnRlcm5hbCh0diwgdiwgcGFja2FnZU5hbWUsIGdyb2suZXZlbnRzLm9uVmlld2VyQWRkZWQsXHJcbiAgICAgICAgY2hhbmdlT3B0aW9uc1NhdmVMYXlvdXQsIG9wdGlvbnMhLmF3YWl0Vmlld2VyKVxyXG5cclxuICAgIC8vMy4gTG9hZCBsYXlvdXRcclxuICAgIGF3YWl0IHRlc3RWaWV3ZXJJbnRlcm5hbCh0diwgdiwgcGFja2FnZU5hbWUsIGdyb2suZXZlbnRzLm9uVmlld0xheW91dEFwcGxpZWQsIGxvYWRMYXlvdXQsIHVuZGVmaW5lZCwgcHJvcHNBbmRMYXlvdXQ/LmxheW91dCxcclxuICAgICAgeyBzYXZlZFByb3BzOiBwcm9wc0FuZExheW91dD8uc2F2ZWRQcm9wcyB9KTtcclxuICAgIGlmIChvcHRpb25zPy5hd2FpdFZpZXdlcilcclxuICAgICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3TGF5b3V0QXBwbGllZCwgbG9hZExheW91dCwgb3B0aW9ucyEuYXdhaXRWaWV3ZXIsXHJcbiAgICAgICAgcHJvcHNBbmRMYXlvdXQ/LmxheW91dCwgeyBzYXZlZFByb3BzOiBwcm9wc0FuZExheW91dD8uc2F2ZWRQcm9wcyB9KTtcclxuXHJcbiAgICAvLzQuIE9wZW4gdmlld2VyIG9uIGFyYml0YXJ5IGRhdGFzZXRcclxuICAgIGlmIChvcHRpb25zPy5hcmJpdHJhcnlEZlRlc3QgIT09IGZhbHNlKSB7XHJcbiAgICAgIHR2LmRhdGFGcmFtZSA9IGNhdERGO1xyXG4gICAgICBhd2FpdCBkZWxheSg1MCk7XHJcbiAgICAgIGF3YWl0IHRlc3RWaWV3ZXJJbnRlcm5hbCh0diwgdiwgcGFja2FnZU5hbWUsIGdyb2suZXZlbnRzLm9uVmlld2VyQWRkZWQpO1xyXG4gICAgICBpZiAob3B0aW9ucz8uYXdhaXRWaWV3ZXIpXHJcbiAgICAgICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCwgdW5kZWZpbmVkLCBvcHRpb25zIS5hd2FpdFZpZXdlcik7XHJcbiAgICB9XHJcblxyXG4gICAgLy81LiBDYWxsIHBvc3Rwb25lZCBmaWx0ZXJpbmdcclxuICAgIGF3YWl0IHRlc3RWaWV3ZXJJbnRlcm5hbCh0diwgdiwgcGFja2FnZU5hbWUsIGdyb2suZXZlbnRzLm9uVmlld2VyQWRkZWQsIGZpbHRlckFzeW5jKTtcclxuICAgIGlmIChvcHRpb25zPy5hd2FpdFZpZXdlcilcclxuICAgICAgYXdhaXQgdGVzdFZpZXdlckludGVybmFsKHR2LCB2LCBwYWNrYWdlTmFtZSwgZ3Jvay5ldmVudHMub25WaWV3ZXJBZGRlZCwgZmlsdGVyQXN5bmMsIG9wdGlvbnMhLmF3YWl0Vmlld2VyKTtcclxuXHJcbiAgfSBmaW5hbGx5IHtcclxuICAgIC8vIGNsb3NlQWxsKCkgaXMgaGFuZGxpbmcgYnkgY29tbW9uIHRlc3Qgd29ya2Zsb3dcclxuICAgIC8vIGdyb2suc2hlbGwuY2xvc2VBbGwoKTtcclxuICAgIC8vIERHLkJhbGxvb24uY2xvc2VBbGwoKTtcclxuICB9XHJcbn1cclxuIl19