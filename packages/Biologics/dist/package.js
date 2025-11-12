var biologics;
/******/ (() => { // webpackBootstrap
/******/ 	"use strict";
/******/ 	var __webpack_modules__ = ({

/***/ "../../libraries/db-explorer/src/db-explorer.ts":
/*!******************************************************!*\
  !*** ../../libraries/db-explorer/src/db-explorer.ts ***!
  \******************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   DBExplorer: () => (/* binding */ DBExplorer)
/* harmony export */ });
/* harmony import */ var datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! datagrok-api/grok */ "datagrok-api/grok");
/* harmony import */ var datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var datagrok_api_dg__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! datagrok-api/dg */ "datagrok-api/dg");
/* harmony import */ var datagrok_api_dg__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _object_handlers__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./object-handlers */ "../../libraries/db-explorer/src/object-handlers.ts");
/* eslint-disable max-len */



class DBExplorer {
    connectionName;
    schemaName;
    nqName;
    dataSourceName;
    schemasLoaded = false;
    connection = null;
    references = {};
    referencedBy = {};
    _dbLoadPromise;
    objHandlers = [];
    loadingFailed = false;
    constructor(connectionName, schemaName, nqName, dataSourceName) {
        this.connectionName = connectionName;
        this.schemaName = schemaName;
        this.nqName = nqName;
        this.dataSourceName = dataSourceName;
        const handler = new _object_handlers__WEBPACK_IMPORTED_MODULE_2__.DBExplorerObjectHandler({ valueConverter: (a) => a,
            joinOptions: [] }, async () => {
            return await this.dbSchema;
        }, this.schemaName);
        datagrok_api_dg__WEBPACK_IMPORTED_MODULE_1__.ObjectHandler.register(handler);
        this.objHandlers.push(handler);
    }
    get dbSchema() {
        if (this.loadingFailed)
            return Promise.resolve(null);
        this._dbLoadPromise ??= this.loadDbSchema();
        return this._dbLoadPromise
            .then(() => ({ schema: { references: this.references, referencedBy: this.referencedBy }, connection: this.connection }))
            .catch((e) => {
            console.error(e);
            return null;
        });
    }
    async loadDbSchema() {
        try {
            const connections = await datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.dapi.connections.filter(`name="${this.connectionName}"`).list();
            this.connection = connections.find((c) => (!this.nqName || c.nqName?.toLowerCase() === this.nqName.toLowerCase()) && (!this.dataSourceName || c.dataSource?.toLowerCase() === this.dataSourceName.toLowerCase())) ?? null;
            if (this.connection == null) {
                console.warn(`Connection ${this.connectionName} not found, Object handlers not registered`);
                this.loadingFailed = true;
                return;
            }
            const tables = await datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.dapi.connections.getSchema(this.connection, this.schemaName);
            if (!tables)
                throw new Error(`Schema ${this.schemaName} not found`);
            tables.forEach((table) => {
                const tableName = table.friendlyName ?? table.name;
                this.references[tableName] = {};
                const t = this.references[tableName];
                table.columns.forEach((column) => {
                    const ref = column.referenceInfo;
                    if (ref && ref.table && ref.column) {
                        t[column.name] = { refTable: ref.table, refColumn: ref.column };
                        if (!this.referencedBy[ref.table])
                            this.referencedBy[ref.table] = {};
                        if (!this.referencedBy[ref.table][ref.column])
                            this.referencedBy[ref.table][ref.column] = [];
                        this.referencedBy[ref.table][ref.column].push({ refTable: tableName, refColumn: column.name });
                    }
                });
            });
            this.schemasLoaded = true;
        }
        catch (_e) {
            this.loadingFailed = true;
            console.warn('Failed to load DB schema, Object handlers not registered');
            console.error(_e);
        }
    }
    async addCustomRelation(tableName, columnName, refTable, refColumn) {
        if (!this.schemasLoaded) {
            try {
                await this._dbLoadPromise;
            }
            catch (_e) {
                console.error(_e);
            }
        }
        if (!this.schemasLoaded || !this.connection || !this.references || !this.referencedBy) {
            console.warn('Failed to add custom relation, DB schema not loaded');
            return;
        }
        if (!this.references[tableName])
            this.references[tableName] = {};
        this.references[tableName][columnName] = { refTable, refColumn };
        if (!this.referencedBy[refTable])
            this.referencedBy[refTable] = {};
        this.referencedBy[refTable][refColumn] ??= [];
        this.referencedBy[refTable][refColumn].push({ refTable: tableName, refColumn: columnName });
    }
    addEntryPoint(semanticType, tableName, columnName, options) {
        const schemaPromise = () => this.dbSchema;
        const fullOpts = { ...{ valueConverter: (a) => a, joinOptions: [] }, ...(options ?? {}) };
        const handler = new _object_handlers__WEBPACK_IMPORTED_MODULE_2__.SemValueObjectHandler(semanticType, tableName, columnName, fullOpts, schemaPromise, this.schemaName);
        datagrok_api_dg__WEBPACK_IMPORTED_MODULE_1__.ObjectHandler.register(handler);
        this.objHandlers.push(handler);
        return this;
    }
    addJoinOptions(joinOpts) {
        this.objHandlers.forEach((handler) => handler.options.joinOptions.push(...joinOpts));
        return this;
    }
    addCustomRenderer(check, renderer) {
        this.objHandlers.forEach((handler) => handler.addCustomRenderer(check, renderer));
        return this;
    }
    addHeaderReplacers(replacers) {
        this.objHandlers.forEach((handler) => handler.addHeaderReplacers(replacers));
        return this;
    }
    addEntryPointValueConverter(func) {
        this.objHandlers.forEach((handler) => handler.options.valueConverter = func);
        return this;
    }
    addDefaultHeaderReplacerColumns(columns) {
        this.objHandlers.forEach((handler) => handler.addDefaultHeaderReplacerColumns(columns));
        return this;
    }
    addUniqueColumns(columns) {
        this.objHandlers.forEach((handler) => handler.addUniqueColumns(columns));
        return this;
    }
    addCustomSelectedColumns(columns) {
        this.objHandlers.forEach((handler) => handler.addCustomSelectedColumns(columns));
        return this;
    }
    static initFromConfig(config) {
        const exp = new DBExplorer(config.connectionName, config.schemaName, config.nqName, config.dataSourceName);
        for (const [semType, entry] of Object.entries(config.entryPoints))
            exp.addEntryPoint(semType, entry.table, entry.column, { regexpExample: entry.regexpExample });
        if (config.joinOptions)
            exp.addJoinOptions(config.joinOptions);
        if (config.headerNames)
            exp.addHeaderReplacers(config.headerNames);
        if (config.uniqueColumns)
            exp.addUniqueColumns(config.uniqueColumns);
        if (config.customSelectedColumns)
            exp.addCustomSelectedColumns(config.customSelectedColumns);
        return exp;
    }
    static async initFromConfigPath(_package, configPath = 'db-explorer/db-explorer-config.json') {
        try {
            const config = await _package.files.readAsText(configPath);
            return DBExplorer.initFromConfig(JSON.parse(config));
        }
        catch (e) {
            datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.shell.error(`Failed to load db-explorer config from ${configPath}`);
            console.error(e);
        }
        return null;
    }
}


/***/ }),

/***/ "../../libraries/db-explorer/src/object-handlers.ts":
/*!**********************************************************!*\
  !*** ../../libraries/db-explorer/src/object-handlers.ts ***!
  \**********************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   DBExplorerObjectHandler: () => (/* binding */ DBExplorerObjectHandler),
/* harmony export */   SemValueObjectHandler: () => (/* binding */ SemValueObjectHandler)
/* harmony export */ });
/* harmony import */ var datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! datagrok-api/grok */ "datagrok-api/grok");
/* harmony import */ var datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! datagrok-api/ui */ "datagrok-api/ui");
/* harmony import */ var datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! datagrok-api/dg */ "datagrok-api/dg");
/* harmony import */ var datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./types */ "../../libraries/db-explorer/src/types.ts");
/* harmony import */ var _renderer__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./renderer */ "../../libraries/db-explorer/src/renderer.ts");
/* eslint-disable max-len */





class DBExplorerObjectHandler extends datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__.ObjectHandler {
    options;
    schemaInfoPromise;
    get type() {
        return 'db-explorer-value';
    }
    renderer;
    isApplicable(x) {
        return x instanceof _types__WEBPACK_IMPORTED_MODULE_3__.DBValueObject;
    }
    constructor(options, schemaInfoPromise, schemaName) {
        super();
        this.options = options;
        this.schemaInfoPromise = schemaInfoPromise;
        this.renderer = new _renderer__WEBPACK_IMPORTED_MODULE_4__.DBExplorerRenderer(schemaInfoPromise, schemaName, options);
    }
    addCustomRenderer(check, renderer) {
        this.renderer.addCustomRenderer(check, renderer);
    }
    addDefaultHeaderReplacerColumns(columns) {
        this.renderer.addDefaultHeaderReplacerColumns(columns);
        return this;
    }
    addCustomSelectedColumns(columns) {
        this.renderer.addCustomSelectedColumns(columns);
    }
    addUniqueColumns(columns) {
        this.renderer.addUniqueColNames(columns);
    }
    addHeaderReplacers(replacers) {
        this.renderer.addHeaderReplacers(replacers);
    }
    renderInnerCard(x) {
        return datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.wait(async () => {
            return await this.renderer.renderTable(x.table, x.column, this.options.valueConverter(x.value ?? ''), this.options.joinOptions);
        });
    }
    renderCard(x, _context) {
        const c = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.card(this.renderInnerCard(x));
        c.style.width = 'unset';
        c.addEventListener('click', () => (datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.shell.o = x));
        return c;
    }
    renderTooltip(x, _context) {
        return this.renderInnerCard(x);
    }
    renderInnerProperties(tableName, x, paneName) {
        const acc = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.accordion(tableName);
        acc.addPane(paneName, () => this.renderInnerCard(x), true);
        this.renderer
            .getTable(tableName, x.column, this.options.valueConverter(x.value ?? ''), this.options.joinOptions)
            .then((df) => {
            this.renderer.renderAssociations(acc, this.schemaInfoPromise, tableName, df);
        });
        return acc;
    }
    renderProperties(x, _context) {
        const acc = this.renderInnerProperties(x.table, x, this.options.valueConverter(x.value ?? '').toString());
        if (x.semValue) {
            const origAcc = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.panels.infoPanel(x.semValue);
            origAcc.context = x.semValue;
            if (x.semValue.semType)
                x.semValue.tags[datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__.Tags.Quality] = x.semValue.semType;
            origAcc?.end();
            return datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.divV([acc.root, origAcc.root]);
        }
        return acc.root;
    }
}
class SemValueObjectHandler extends DBExplorerObjectHandler {
    semanticType;
    tableName;
    columnName;
    get type() {
        return this.semanticType;
    }
    _rgExample = null;
    get regexpExample() {
        return this._rgExample;
    }
    isApplicable(x) {
        return x instanceof datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__.SemanticValue && x.semType == this.semanticType;
    }
    constructor(semanticType, tableName, columnName, options, schemaInfoPromise, schemaName) {
        super(options, schemaInfoPromise, schemaName);
        this.semanticType = semanticType;
        this.tableName = tableName;
        this.columnName = columnName;
        if (options && options.regexpExample && options.regexpExample.nonVariablePart && options.regexpExample.regexpMarkup)
            this._rgExample = options.regexpExample;
    }
    renderCard(x, context) {
        return super.renderCard(new _types__WEBPACK_IMPORTED_MODULE_3__.DBValueObject(this.tableName, this.columnName, this.options.valueConverter(x.value ?? ''), x), context);
    }
    renderTooltip(x, context) {
        return super.renderTooltip(new _types__WEBPACK_IMPORTED_MODULE_3__.DBValueObject(this.tableName, this.columnName, this.options.valueConverter(x.value ?? ''), x), context);
    }
    renderProperties(x, context) {
        return super.renderProperties(new _types__WEBPACK_IMPORTED_MODULE_3__.DBValueObject(this.tableName, this.columnName, this.options.valueConverter(x.value ?? ''), x), context);
    }
}


/***/ }),

/***/ "../../libraries/db-explorer/src/query.ts":
/*!************************************************!*\
  !*** ../../libraries/db-explorer/src/query.ts ***!
  \************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   queryDB: () => (/* binding */ queryDB),
/* harmony export */   queryDBMultiple: () => (/* binding */ queryDBMultiple)
/* harmony export */ });
/* harmony import */ var datagrok_api_dg__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! datagrok-api/dg */ "datagrok-api/dg");
/* harmony import */ var datagrok_api_dg__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_0__);

async function queryDB(connection, tableName, match, matchValue, schemaName, joinOptions = []) {
    if (connection == null)
        return datagrok_api_dg__WEBPACK_IMPORTED_MODULE_0__.DataFrame.create(0);
    const matchValueStr = typeof matchValue === 'string' ? `'${matchValue}'` : matchValue;
    const startCharCode = 98; // ascii b
    const applicableJoins = joinOptions.filter((opt) => opt.fromTable === tableName);
    const tableAliases = applicableJoins.map((_, i) => String.fromCharCode(startCharCode + i));
    const otherCols = applicableJoins
        .map((opt, i) => {
        const alias = tableAliases[i];
        return opt.select.map((col) => `${alias}.${col}`).join(', ');
    })
        .join(', ') + (applicableJoins.length > 0 ? ', ' : '');
    const joinStr = applicableJoins
        .map((opt, i) => {
        const alias = tableAliases[i];
        return `left join "${schemaName}".${opt.tableName} ${alias} on a.${opt.columnName} = ${alias}.${opt.onColumn}`;
    })
        .join(' ');
    const q = connection.query('getDBValueInfo', `
        --name: getDBValueInfo
        --output: dataframe result
        select ${otherCols} a.* 
        from "${schemaName}".${tableName} a
        ${joinStr}
        where a.${match} = ${matchValueStr}
    `);
    return (await q.apply({})) ?? datagrok_api_dg__WEBPACK_IMPORTED_MODULE_0__.DataFrame.create(0);
}
async function queryDBMultiple(connection, tableName, match, matchValues) {
    if (connection == null)
        return datagrok_api_dg__WEBPACK_IMPORTED_MODULE_0__.DataFrame.create(0);
    if ((matchValues?.length ?? 0) < 1)
        return datagrok_api_dg__WEBPACK_IMPORTED_MODULE_0__.DataFrame.create(0);
    // matchValues can be more than 1000, which is the limit for oracle and others
    // we need to split it in chunks of 500, to be on the safe side
    const chunkSize = 500;
    let outDataFrame = null;
    for (let i = 0; i < matchValues.length; i += chunkSize) {
        const chunk = matchValues.slice(i, i + chunkSize);
        const q = connection.query('getDBValueInfoMult', `
            --name: getDBValueInfoMult
            --output: dataframe result
            select * from ${tableName} where ${match} in (${chunk.map((v) => `'${v}'`).join(',')})
        `);
        const res = await q.apply({});
        if (res) {
            if (!outDataFrame)
                outDataFrame = res;
            else
                outDataFrame.append(res, true);
        }
    }
    return outDataFrame ?? datagrok_api_dg__WEBPACK_IMPORTED_MODULE_0__.DataFrame.create(0);
}


/***/ }),

/***/ "../../libraries/db-explorer/src/renderer.ts":
/*!***************************************************!*\
  !*** ../../libraries/db-explorer/src/renderer.ts ***!
  \***************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   DBExplorerRenderer: () => (/* binding */ DBExplorerRenderer),
/* harmony export */   MAX_MULTIROW_VALUES: () => (/* binding */ MAX_MULTIROW_VALUES),
/* harmony export */   getLoaderDiv: () => (/* binding */ getLoaderDiv),
/* harmony export */   imageRenderer: () => (/* binding */ imageRenderer),
/* harmony export */   moleculeRenderer: () => (/* binding */ moleculeRenderer),
/* harmony export */   ownIdRenderer: () => (/* binding */ ownIdRenderer),
/* harmony export */   rawImageRenderer: () => (/* binding */ rawImageRenderer),
/* harmony export */   textRenderer: () => (/* binding */ textRenderer)
/* harmony export */ });
/* harmony import */ var datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! datagrok-api/grok */ "datagrok-api/grok");
/* harmony import */ var datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! datagrok-api/ui */ "datagrok-api/ui");
/* harmony import */ var datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! datagrok-api/dg */ "datagrok-api/dg");
/* harmony import */ var datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./types */ "../../libraries/db-explorer/src/types.ts");
/* harmony import */ var _query__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./query */ "../../libraries/db-explorer/src/query.ts");
/* eslint-disable max-len */





const MAX_MULTIROW_VALUES = 10;
/** normal renderer supporting copy and elipsis */
function textRenderer(value, withTooltip = true) {
    const nameHost = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.div(value.toString());
    nameHost.style.maxWidth = '200px';
    nameHost.style.overflow = 'hidden';
    nameHost.style.whiteSpace = 'nowrap';
    nameHost.style.textOverflow = 'ellipsis';
    // ability to copy
    const menu = datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__.Menu.popup();
    menu.item('Copy', () => {
        navigator?.clipboard?.writeText(value.toString());
    });
    nameHost.addEventListener('contextmenu', (e) => {
        e.preventDefault();
        e.stopImmediatePropagation();
        setTimeout(() => menu.show());
    });
    // approximately what will fit in 150 px
    if (value.toString().length > 20 && withTooltip)
        datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.tooltip.bind(nameHost, value.toString());
    return nameHost;
}
function moleculeRenderer(value) {
    try {
        const molDiv = datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.chem.drawMolecule(value, 200, 200);
        return molDiv;
    }
    catch (e) {
        console.error(e);
    }
    return textRenderer(value);
}
function imageRenderer(fullUrl, useProxy = true) {
    const nameHost = textRenderer(fullUrl);
    const loaderDiv = getLoaderDiv();
    const host = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.divH([nameHost, loaderDiv]);
    datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.tooltip.bind(nameHost, 'Unable to get image');
    async function replaceWithImage() {
        try {
            const qRes = useProxy ? await datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.dapi.fetchProxy(fullUrl, {}) : await fetch(fullUrl);
            if (!qRes)
                throw new Error('');
            const blob = await qRes.blob();
            if (!blob)
                throw new Error('');
            const canvas = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.canvas(200, 200);
            const ctx = canvas.getContext('2d');
            if (!ctx)
                throw new Error('');
            const image = new Image();
            image.onload = () => {
                ctx.drawImage(image, 0, 0, 200, 200);
            };
            image.src = URL.createObjectURL(blob);
            nameHost.remove();
            host.appendChild(canvas);
        }
        catch (e) {
            console.error(e);
        }
        finally {
            loaderDiv.remove();
        }
    }
    replaceWithImage();
    return host;
}
function rawImageRenderer(rawImage) {
    const nameHost = textRenderer('image');
    const loaderDiv = getLoaderDiv();
    const host = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.divH([nameHost, loaderDiv]);
    datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.tooltip.bind(nameHost, 'Unable to render image');
    async function replaceWithImage() {
        try {
            if (!rawImage)
                throw new Error('Empty image string');
            const canvas = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.canvas(200, 200);
            const ctx = canvas.getContext('2d');
            if (!ctx)
                throw new Error('');
            const image = new Image();
            image.onload = () => {
                ctx.drawImage(image, 0, 0, 200, 200);
            };
            image.src = 'data:image/png;base64,' + rawImage;
            nameHost.remove();
            host.appendChild(canvas);
        }
        catch (e) {
            console.error(e);
        }
        finally {
            loaderDiv.remove();
        }
    }
    replaceWithImage();
    return host;
}
function ownIdRenderer(id, tableName, colName) {
    const nameHost = textRenderer(id.toString(), false);
    nameHost.style.color = 'var(--blue-1)';
    nameHost.style.cursor = 'pointer';
    datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.tooltip.bind(nameHost, 'Click to explore');
    nameHost.addEventListener('click', (e) => {
        e.stopImmediatePropagation();
        datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.shell.o = new _types__WEBPACK_IMPORTED_MODULE_3__.DBValueObject(tableName, colName, id);
    });
    return nameHost;
}
function removeEmptyCols(df) {
    df.columns.names().forEach((colName) => {
        if (!df.col(colName) || df.col(colName).isNone(0))
            df.columns.remove(colName);
    });
    return df;
}
function getLoaderDiv() {
    const div = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.div([], { style: { width: '50px', height: '24px', position: 'relative' } });
    div.innerHTML = `<div class="grok-loader"><div></div><div></div><div></div><div></div></div>`;
    return div;
}
class DBExplorerRenderer {
    schemaInfoPromise;
    schemaName;
    entryPointOptions;
    valueReplacers = [];
    customRenderers = [];
    headerReplacers = {};
    uniqueColNames = {};
    customSelectedColumns = {};
    defaultHeaderReplacerColumns = ['name'];
    constructor(schemaInfoPromise, schemaName, entryPointOptions) {
        this.schemaInfoPromise = schemaInfoPromise;
        this.schemaName = schemaName;
        this.entryPointOptions = entryPointOptions;
    }
    addHeaderReplacers(replacers) {
        this.headerReplacers = { ...this.headerReplacers, ...replacers };
    }
    addDefaultHeaderReplacerColumns(columns) {
        this.defaultHeaderReplacerColumns = [...this.defaultHeaderReplacerColumns, ...columns];
    }
    addCustomSelectedColumns(columns) {
        const columnSets = {};
        Object.entries(columns).forEach(([tableName, cols]) => {
            columnSets[tableName] = new Set(cols);
        });
        this.customSelectedColumns = { ...this.customSelectedColumns, ...columnSets };
    }
    addUniqueColNames(cols) {
        this.uniqueColNames = { ...this.uniqueColNames, ...cols };
    }
    addValueReplacer(check, replacer) {
        this.valueReplacers.push({ check, replacer });
    }
    addCustomRenderer(check, renderer) {
        this.customRenderers.push({ check, renderer });
    }
    refIdRenderer(connection, id, tableName, colName) {
        const loaderDiv = getLoaderDiv();
        const nameHost = textRenderer(id.toString(), false);
        const host = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.divH([nameHost, loaderDiv]);
        function unableToRetrieve(reason = `Unable to retrieve ${tableName} information`) {
            loaderDiv.remove();
            datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.tooltip.bind(nameHost, reason);
        }
        try {
            (0,_query__WEBPACK_IMPORTED_MODULE_4__.queryDB)(connection, tableName, colName, id, this.schemaName, this.entryPointOptions.joinOptions)
                .then((df) => {
                if (df.rowCount == 0) {
                    unableToRetrieve();
                    return;
                }
                const clearedDF = removeEmptyCols(df);
                const replacer = this.valueReplacers.find((replacer) => replacer.check(tableName, colName, id));
                if (replacer)
                    nameHost.textContent = replacer.replacer(clearedDF, id);
                datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.tooltip.bind(nameHost, () => datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.wait(async () => await this.renderDataFrame(clearedDF, tableName)));
                nameHost.addEventListener('click', (e) => {
                    e.stopImmediatePropagation();
                    datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.shell.o = new _types__WEBPACK_IMPORTED_MODULE_3__.DBValueObject(tableName, colName, id);
                });
                nameHost.style.color = 'var(--blue-1)';
                nameHost.style.cursor = 'pointer';
            })
                .finally(() => {
                loaderDiv.remove();
            });
        }
        catch (e) {
            console.error(e);
            unableToRetrieve();
        }
        return host;
    }
    async getTable(tableName, match, matchValue, joinOptions = []) {
        const schemaAndConnection = await this.schemaInfoPromise();
        const res = await (0,_query__WEBPACK_IMPORTED_MODULE_4__.queryDB)(schemaAndConnection?.connection ?? null, tableName, match, matchValue, this.schemaName, joinOptions);
        return res;
    }
    async renderMultiRowTable(tableName, match, matchValue, joinOptions = []) {
        const df = await this.getTable(tableName, match, matchValue, joinOptions);
        if (df.rowCount < 1)
            return { root: datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.divText('ID not found'), rowCount: 0 };
        if (df.rowCount === 1)
            return { root: this.renderDataFrame(df, tableName), rowCount: 1 };
        const acc = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.accordion(`Multiple rows for ${tableName}`);
        const dfMaxCount = Math.min(df.rowCount, 10);
        const replaceColName = this.headerReplacers[tableName];
        for (let i = 0; i < dfMaxCount; i++) {
            let paneName = `Row ${i + 1}`;
            if (replaceColName && df.col(replaceColName)?.get(i)) {
                paneName = df.col(replaceColName).get(i).toString();
            }
            else {
                const f = this.defaultHeaderReplacerColumns.find((col) => df.col(col) && df.col(col).get(i));
                if (f)
                    paneName = df.col(f).get(i).toString();
            }
            acc.addPane(paneName, () => {
                const rowBitset = datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__.BitSet.create(df.rowCount);
                rowBitset.set(i, true);
                const rowDf = df.clone(rowBitset);
                return datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.wait(async () => await this.renderDataFrame(rowDf, tableName));
            });
        }
        return { root: acc.root, rowCount: df.rowCount };
    }
    async renderTable(tableName, match, matchValue, joinOptions = []) {
        const df = await this.getTable(tableName, match, matchValue, joinOptions);
        if (df.rowCount < 1)
            return datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.divText('ID not found');
        return this.renderDataFrame(df, tableName);
    }
    async renderDataFrame(df, tableName) {
        const schemaAndConnection = await this.schemaInfoPromise();
        if (!schemaAndConnection)
            return datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.divText('Schema information is not available');
        const clearedDF = removeEmptyCols(df);
        let entries = Object.entries(clearedDF.toJson()[0]);
        if (this.customSelectedColumns[tableName]) {
            entries = entries.filter(([key]) => this.customSelectedColumns[tableName].has(key));
            // reorder entries according to customSelectedColumns if present
            const cols = Array.from(this.customSelectedColumns[tableName]);
            entries.sort((a, b) => {
                return cols.indexOf(a[0]) - cols.indexOf(b[0]);
            });
        }
        const mainTable = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.table(entries, (entry) => {
            const key = entry[0];
            const value = entry[1];
            const cutomRenderer = this.customRenderers.find((renderer) => renderer.check(tableName, key, value));
            if (cutomRenderer)
                return [textRenderer(key), cutomRenderer.renderer(value, schemaAndConnection.connection)];
            const refInfo = schemaAndConnection.schema.references?.[tableName]?.[key];
            if (refInfo)
                return [textRenderer(key), this.refIdRenderer(schemaAndConnection.connection, value, refInfo.refTable, refInfo.refColumn)];
            const isUnqueCol = this.uniqueColNames[tableName] === key;
            if (isUnqueCol)
                return [textRenderer(key), ownIdRenderer(value, tableName, key)];
            return [textRenderer(key), textRenderer(value)];
        });
        return datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.divV([mainTable]);
    }
    async renderAssociations(acc, schemaPromise, curTable, curDf) {
        const schemaAndConnection = await schemaPromise();
        if (!schemaAndConnection)
            return;
        const addAllAssociated = async (tableName, colName, value) => {
            const pi = datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__.TaskBarProgressIndicator.create('Opening all associated entries');
            try {
                const assocDf = await this.getTable(tableName, colName, value, this.entryPointOptions.joinOptions);
                if (assocDf) {
                    if (assocDf.rowCount === 0) {
                        datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.shell.info(`No associated ${tableName} entries found`);
                        return;
                    }
                    assocDf.name = `${tableName}`;
                    datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.shell.addTableView(assocDf);
                }
            }
            catch (e) {
                console.error(e);
                datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.shell.error(`Failed to open associated ${tableName} entries`);
            }
            finally {
                pi.close();
            }
        };
        const addIconToPane = (pane, tableName, colName, value) => {
            const icon = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.icons.add(() => { }, `Add all associated ${tableName} entries to workspace`);
            // need separate event to avoid click on accordion pane
            icon.addEventListener('click', (e) => {
                e.stopImmediatePropagation();
                e.preventDefault();
                addAllAssociated(tableName, colName, value);
            });
            pane.root.getElementsByClassName('d4-accordion-pane-header')?.[0]?.appendChild(icon);
        };
        const attachOpenInWorkspaceMenu = (tableName, colName, value, pane) => {
            const menu = datagrok_api_dg__WEBPACK_IMPORTED_MODULE_2__.Menu.popup();
            menu.item(`Add all associated ${tableName} entries to workspace`, async () => {
                addAllAssociated(tableName, colName, value);
            });
            pane.root.addEventListener('contextmenu', (e) => {
                e.preventDefault();
                setTimeout(() => menu.show());
            });
            addIconToPane(pane, tableName, colName, value);
        };
        const refTable = schemaAndConnection.schema.referencedBy[curTable];
        if (!refTable || Object.keys(refTable).length == 0)
            return;
        const _linksPane = acc.addPane('Links', () => {
            const linksAcc = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.accordion(`Links to ${curTable}`);
            if ((curDf?.rowCount ?? 0) === 0)
                return datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.divText('No data');
            Object.entries(refTable).forEach(([refedColumn, refInfo]) => {
                const val = curDf.col(refedColumn)?.get(0);
                if (!val)
                    return;
                const _pane = linksAcc.addPane(refedColumn, () => {
                    const colAcc = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.accordion(`Links to ${curTable}_${refedColumn}`);
                    // there can be cases when same column is referneced by two or more columns in other table.
                    // example is chembl table molecule_hiearchy where all 3 columns are referencing to molregno
                    // need to account for such cases
                    const tableRefs = new Map();
                    refInfo.forEach((ref) => {
                        if (!tableRefs.has(ref.refTable))
                            tableRefs.set(ref.refTable, []);
                        tableRefs.get(ref.refTable).push(ref);
                    });
                    tableRefs.forEach((refInfos, refTable) => {
                        const singleColPane = colAcc.addPane(refTable, () => {
                            if (refInfos.length === 1) {
                                return datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.wait(async () => {
                                    const res = await this.renderMultiRowTable(refInfos[0].refTable, refInfos[0].refColumn, val, this.entryPointOptions.joinOptions);
                                    if (res.rowCount > MAX_MULTIROW_VALUES) {
                                        singleColPane.name = `${singleColPane.name} (${MAX_MULTIROW_VALUES} / ${res.rowCount})`;
                                        addIconToPane(singleColPane, refInfos[0].refTable, refInfos[0].refColumn, val);
                                    }
                                    return res.root;
                                });
                            }
                            const multiTableAcc = datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.accordion(`Multiple links to ${refTable}`);
                            refInfos.forEach((ref) => {
                                const colPane = multiTableAcc.addPane(ref.refColumn, () => {
                                    return datagrok_api_ui__WEBPACK_IMPORTED_MODULE_1__.wait(async () => {
                                        const res = await this.renderMultiRowTable(ref.refTable, ref.refColumn, val, this.entryPointOptions.joinOptions);
                                        if (res.rowCount > MAX_MULTIROW_VALUES) {
                                            colPane.name = `${singleColPane.name} (${MAX_MULTIROW_VALUES} / ${res.rowCount})`;
                                            addIconToPane(colPane, ref.refTable, ref.refColumn, val);
                                        }
                                        return res.root;
                                    });
                                });
                                attachOpenInWorkspaceMenu(ref.refTable, ref.refColumn, val, colPane);
                            });
                            return multiTableAcc.root;
                        });
                        if (refInfos.length === 1)
                            attachOpenInWorkspaceMenu(refInfos[0].refTable, refInfos[0].refColumn, val, singleColPane);
                    });
                    return colAcc.root;
                });
            });
            return linksAcc.root;
        });
    }
}


/***/ }),

/***/ "../../libraries/db-explorer/src/types.ts":
/*!************************************************!*\
  !*** ../../libraries/db-explorer/src/types.ts ***!
  \************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   DBValueObject: () => (/* binding */ DBValueObject)
/* harmony export */ });
class DBValueObject {
    table;
    column;
    value;
    semValue;
    constructor(table, column, value, semValue) {
        this.table = table;
        this.column = column;
        this.value = value;
        this.semValue = semValue;
    }
}


/***/ }),

/***/ "./src/config.ts":
/*!***********************!*\
  !*** ./src/config.ts ***!
  \***********************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   biologicsConfig: () => (/* binding */ biologicsConfig)
/* harmony export */ });
const biologicsConfig = {
    'connectionName': 'Biologics',
    'schemaName': 'biologics',
    'dataSourceName': 'postgresDart',
    'entryPoints': {
        'DG_BIOLOGICS_DRUG_ID': {
            'table': 'drugs',
            'column': 'identifier',
            'regexpExample': {
                'example': 'GROKMOL-000001',
                'nonVariablePart': 'GROKMOL-',
                'regexpMarkup': 'GROKMOL-{d6}'
            }
        },
        'DG_BIOLOGICS_SEQUENCE_ID': {
            'table': 'sequences',
            'column': 'identifier',
            'regexpExample': {
                'example': 'GROKSEQ-000001',
                'nonVariablePart': 'GROKSEQ-',
                'regexpMarkup': 'GROKSEQ-{d6}'
            }
        },
        'DG_BIOLOGICS_LINKER_ID': {
            'table': 'linkers',
            'column': 'identifier',
            'regexpExample': {
                'example': 'GROKLINKER-000001',
                'nonVariablePart': 'GROKLINKER-',
                'regexpMarkup': 'GROKLINKER-{d6}'
            }
        },
        'DG_BIOLOGICS_ADC_ID': {
            'table': 'adc',
            'column': 'identifier',
            'regexpExample': {
                'example': 'GROKADC-000001',
                'nonVariablePart': 'GROKADC-',
                'regexpMarkup': 'GROKADC-{d6}'
            }
        },
        'DG_BIOLOGICS_ORGANISM_ID': {
            'table': 'target_organisms',
            'column': 'identifier',
            'regexpExample': {
                'example': 'GROKORG-000001',
                'nonVariablePart': 'GROKORG-',
                'regexpMarkup': 'GROKORG-{d6}'
            }
        },
        'DG_BIOLOGICS_PURIFICATION_ID': {
            'table': 'purification_batches',
            'column': 'identifier',
            'regexpExample': {
                'example': 'GROKPUR-000001',
                'nonVariablePart': 'GROKPUR-',
                'regexpMarkup': 'GROKPUR-{d6}'
            }
        },
        'DG_BIOLOGICS_EXPRESSION_ID': {
            'table': 'expression_batches',
            'column': 'identifier',
            'regexpExample': {
                'example': 'GROKEXP-000001',
                'nonVariablePart': 'GROKEXP-',
                'regexpMarkup': 'GROKEXP-{d6}'
            }
        }
    },
    'joinOptions': [
        {
            'fromTable': 'adc',
            'columnName': 'drug_id',
            'tableName': 'drugs',
            'onColumn': 'id',
            'select': ['smiles as compound_structure', 'name as compound_name']
        },
        {
            'fromTable': 'adc',
            'columnName': 'antibody_id',
            'tableName': 'sequences',
            'onColumn': 'id',
            'select': ['sequence as antibody_sequence', 'name as sequence_name']
        },
        {
            'fromTable': 'assay_results',
            'columnName': 'assay_id',
            'tableName': 'assay_types',
            'onColumn': 'id',
            'select': ['name']
        }
    ],
    'headerNames': {
        'smiles': 'Compound'
    },
    'uniqueColumns': {
        'adc': 'identifier',
        'drugs': 'identifier',
        'sequences': 'identifier',
        'linkers': 'identifier',
        'target_organisms': 'identifier',
        'purification_batches': 'identifier',
        'expression_batches': 'identifier'
    },
    'customSelectedColumns': {
        'adc': [
            'name',
            'identifier',
            'antibody_sequence',
            'antibody_id',
            'drug_id',
            'linker_id',
            'compound_structure',
            'glyph'
        ]
    }
};


/***/ }),

/***/ "./src/glyphs/glyphs.ts":
/*!******************************!*\
  !*** ./src/glyphs/glyphs.ts ***!
  \******************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   glyphPool: () => (/* binding */ glyphPool)
/* harmony export */ });
// Auto-generated file, do not edit
const glyphPool = [
    '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'
];


/***/ }),

/***/ "./src/randsmiles.ts":
/*!***************************!*\
  !*** ./src/randsmiles.ts ***!
  \***************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   smi: () => (/* binding */ smi)
/* harmony export */ });
const smi = ["CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3", "COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCCCc3cccnc3", "CC(C(=O)NCCS)c1cccc(c1)C(=O)c2ccccc2", "FC(F)(F)c1ccc(OC2CCNCC2)cc1", "CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCCc3ccccc3", "COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5", "CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCO[N+](=O)[O-]", "CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCO", "CN1CCC(CC1)Oc2ccc(cc2)C(F)(F)F", "COc1cc(C)ccc1OC(=O)C(C)c2ccc(CC(C)C)cc2", "CC(C)Cc1ccc(cc1)C(C)C(=O)OCCCc2cccnc2", "COc1ccc(\\\\C=N\\\\NC(=N)N)c(Cl)c1OC", "Nc1ncnc2c1c(Br)cn2[C@@H]3OC[C@@H](O)[C@H]3O", "CNc1ncnc2c1c(I)cn2[C@@H]3O[C@H](C)[C@@H](O)[C@H]3O", "CN1CCC(O)(CC1)c2ccccc2", "OC(COc1ccccc1)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F", "OC(COc1ccc(Cl)cc1)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F", "OC(COc1ccc(Br)cc1)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F", "COC(=O)c1ccccc1OCC(O)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F", "CCC1(CC)CC(CCNC(=O)c2ccc(OC)cc2)OC1=O", "COc1ccc(cc1)N(Cc2ccccc2)Cc3ccc(OC)c(O)c3", "COc1ccc(cc1)N(Cc2ccc(Br)cc2)Cc3ccc(OC)c(O)c3", "COc1ccc(cc1)N(Cc2ccc(F)cc2)Cc3ccc(OC)c(O)c3", "COc1ccc(cc1)N(Cc2ccc(Cl)cc2)Cc3ccc(OC)c(O)c3", "CCCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCCC)C", "COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccc(Cl)cc3Cl)cc1O", "COc1ccc(CN(Cc2ccccc2)Cc3ccc(Br)cc3)cc1O", "COc1ccc(CN(CCc2ccccc2)Cc3ccccc3)cc1O", "COc1ccc(CCN(Cc2ccccc2)Cc3ccc(OC)c(O)c3)cc1", "COc1ccc2cc(ccc2c1)C(C)C(=O)N3CCCC3C(=O)OCCCc4cccnc4", "COc1ccc(CN(Cc2ccccc2)c3ccc(Br)cc3)cc1O", "COc1ccc(CN(Cc2ccccc2)c3cccc(Cl)c3)cc1O", "COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccccc3)cc1O", "COc1ccc(CN(CCc2ccc(OC)c(OC)c2)Cc3ccccc3)cc1O", "CC(C)OC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)C)C", "COc1ccc(CCN(Cc2ccc(Br)cc2)Cc3ccc(OC)c(O)c3)cc1", "COc1ccc(CN(Cc2ccc(Br)cc2)c3ccc(Br)cc3)cc1O", "COc1ccc(CN(Cc2ccc(Br)cc2)c3cccc(Cl)c3)cc1O", "COc1ccc(CN(CCc2ccc(OC)c(OC)c2)Cc3ccc(Br)cc3)cc1O", "CCCN1OC(=CC1=O)C", "COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccc(F)cc3)cc1O", "COc1ccc(CN(CCc2ccc(Br)cc2)Cc3cccc(Cl)c3)cc1O", "COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccc(Cl)cc3)cc1O", "COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccc(Br)cc3)cc1O", "CC(C)(NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CCCNC(=N)N)NC(=O)CN)C(=O)N[C@@H](CC(=O)N)C(=O)O", "COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccc(C)cc3)cc1O", "COc1ccc(CN(CCc2ccc(Br)cc2)Cc3c(F)c(F)c(F)c(F)c3F)cc1O", "COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccc(Cl)c(Cl)c3)cc1O", "COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccc(cc3)C(C)(C)C)cc1O", "CC(C)COC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC(C)C)C", "[I-].C[N+](C)(C)C[C@H]1CO[C@@H](O1)C(c2ccccc2)c3ccccc3", "CCC1(CC)CC(CNC(=O)c2ccc(Br)cc2)OC1=O", "CCC1(CC)CC(CO)OC1=O", "CCC1(CC)CC(CCNC(=O)c2ccc(Br)cc2)OC1=O", "CCC1(CC)CC(CCNC(=O)c2ccc(C)cc2)OC1=O", "CC1=C(C(C(=C(C)N1)C(=O)OC(C)(C)C)c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)(C)C", "CCCCCCCC(=O)NCCC1CC(CC)(CC)C(=O)O1", "CC(=O)C(=O)c1cc2c(cn1)[nH]c3ccccc23", "COC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC)C", "CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC)C", "CC1=C(C(C(=C(C)N1)C(=O)OC2CCCC2)c3csc(n3)c4ccc(Cl)cc4)C(=O)OC5CCCC5", "CC1=C(C(C(=C(C)N1)C(=O)OC2CCCCC2)c3csc(n3)c4ccc(Cl)cc4)C(=O)OC5CCCCC5", "CC1=C(C(C(=C(C)N1)C(=O)OCc2ccccc2)c3csc(n3)c4ccc(Cl)cc4)C(=O)OCc5ccccc5", "CC1=C(C(C(=C(C)N1)C(=O)OCCc2ccccc2)c3csc(n3)c4ccc(Cl)cc4)C(=O)OCCc5ccccc5", "CC(OC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)c4ccccc4)C)c5ccccc5", "CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC)C", "COC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC4CCCCC4)C", "CCCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC)C", "CCCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC)C", "COC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)C)C", "CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)C)C", "COC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC(C)C)C", "CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC(C)C)C", "COC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)(C)C)C", "CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)(C)C)C", "COC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC4CCCC4)C", "CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC4CCCC4)C", "COC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCc4ccccc4)C", "CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCc4ccccc4)C", "COC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCCc4ccccc4)C", "CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCCc4ccccc4)C", "COC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)c4ccccc4)C", "CCCCCCCCCCCCCCCCC(CCCCCCCCCCCCCC)C(=O)OC[C@H]1O[C@H](O[C@H]2O[C@H](COC(=O)C(CCCCCCCCCCCCCC)CCCCCCCCCCCCCCCC)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", "CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)c4ccccc4)C", "COC(=O)C1=C[C@@H](O)[C@@H](O)[C@H](C)C1", "Oc1ccc(cc1)[C@H]2CC(=O)c3c(O)cc(Oc4c(O)cc5O[C@H](CC(=O)c5c4O)c6ccc(O)cc6)cc3O2", "O=C(N\\\\N=C/1\\\\C(=O)Nc2ccccc12)NN3C(=O)c4ccccc4N=C3c5ccccc5", "COc1cc(\\\\C=N\\\\NC(=O)NN2C(=O)c3ccccc3N=C2c4ccccc4)ccc1O", "CC1(C)C2CCC1(C)\\\\C(=N/NC(=O)NN3C(=O)c4ccccc4N=C3c5ccccc5)\\\\C2", "O=C(NN=C(c1ccccc1)c2ccccc2)NN3C(=O)c4ccccc4N=C3c5ccccc5", "CCCCCCCCCCCCC(CCCCCCCCCCCC)C(=O)OC[C@H]1O[C@H](O[C@H]2O[C@H](COC(=O)C(CCCCCCCCCCCC)CCCCCCCCCCCC)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", "CN(C)c1ccc(\\\\C=N\\\\NC(=O)NN2C(=O)c3ccccc3N=C2c4ccccc4)cc1", "COc1ccc2cc(\\\\C=N\\\\NC(=O)C3=C(O)c4ccccc4S(=O)(=O)N3C)c(Cl)nc2c1", "Fc1ccc(F)c(c1)C(=O)OCCCN2C(=O)c3ccccc3C2=O", "[O-][N+](=O)c1cc(ccc1Cl)C(=O)OCCCN2C(=O)c3ccccc3C2=O", "Fc1ccc(F)c(c1)C(=O)OCCN2C(=O)[C@@H]3CC=CC[C@@H]3C2=O", "[O-][N+](=O)c1cc(ccc1Cl)C(=O)OCCN2C(=O)[C@@H]3CC=CC[C@@H]3C2=O", "OC(COc1cccc(c1)C(F)(F)F)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F", "Fc1ccc(F)c(c1)C(=O)OCCN2C(=O)[C@@H]3CCCC[C@@H]3C2=O", "[O-][N+](=O)c1cc(ccc1Cl)C(=O)OCCN2C(=O)[C@@H]3CCCC[C@@H]3C2=O", "[O-][N+](=O)c1cc(ccc1Cl)C(=O)OCCOc2ccc(Br)cc2", "[O-][N+](=O)c1ccc(OCCOC(=O)c2ccc(Cl)c(c2)[N+](=O)[O-])cc1", "Clc1ccc(cc1Cl)C(=O)OCCS(=O)(=O)c2ccccc2", "[O-][N+](=O)c1cc(ccc1Cl)C(=O)OCCc2c[nH]c3ccccc23", "Clc1ccc(cc1Cl)C(=O)OC2Cc3cccc4cccc2c34", "[O-][N+](=O)c1cc(ccc1Cl)C(=O)OC2Cc3cccc4cccc2c34", "Cc1cc(ccc1Br)C(=O)OC2Cc3cccc4cccc2c34", "Clc1ccc(cn1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "Cc1oc(C)c(c1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "CCC(C(=O)OCCN1C(=O)c2ccccc2C1=O)c3ccccc3", "O=C(OCCN1C(=O)c2ccccc2C1=O)c3ccc4ccccc4c3", "Fc1ccc(F)c(c1)C(=O)OCN2C(=O)c3ccccc3C2=O", "[O-][N+](=O)c1cc(ccc1Cl)C(=O)OCN2C(=O)c3ccccc3C2=O", "CC(COC(=O)c1cc(F)ccc1F)N2C(=O)c3ccccc3C2=O", "CC(CN1C(=O)c2ccccc2C1=O)OC(=O)c3cc(F)ccc3F", "CC(CN1C(=O)c2ccccc2C1=O)OC(=O)c3ccc(Cl)c(c3)[N+](=O)[O-]", "Cc1ccc2C(=O)N(CCOC(=O)c3cc(F)ccc3F)C(=O)c2c1", "Oc1ccc(Cl)cc1C(=O)\\\\C=C\\\\c2ccccc2Cl", "Cc1ccc2C(=O)N(CCOC(=O)c3ccc(Cl)c(c3)[N+](=O)[O-])C(=O)c2c1", "CC(C)(C)c1ccc2C(=O)N(CCOC(=O)c3cc(F)ccc3F)C(=O)c2c1", "Fc1cc(F)cc(c1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "Clc1ccc(cc1Cl)C(=O)OCCN2C(=O)c3ccccc3C2=O", "COc1cc(C)ccc1OC(=O)C2CCCN2C(=O)C(C)c3cccc(c3)C(=O)c4ccccc4", "Clc1cc(Cl)cc(c1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "[O-][N+](=O)c1ccc(C(=O)OCCN2C(=O)c3ccccc3C2=O)c(Cl)c1", "Cc1ccc(F)cc1C(=O)OCCN2C(=O)c3ccccc3C2=O", "Fc1ccc(cc1C(=O)OCCN2C(=O)c3ccccc3C2=O)C(F)(F)F", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H](C(C)C)C(=O)O", "Fc1ccc(Cl)cc1C(=O)OCCN2C(=O)c3ccccc3C2=O", "[O-][N+](=O)c1ccc(F)c(c1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "Fc1ccc(Cl)c(c1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "[O-][N+](=O)c1ccc(Cl)c(c1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H](CS)C(=O)O", "Fc1ccc(Br)c(c1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "[O-][N+](=O)c1cc(ccc1Cl)C(=O)OCCN2C(=O)c3ccccc3C2=O", "Cc1cc(ccc1Br)C(=O)OCCN2C(=O)c3ccccc3C2=O", "Clc1ncccc1C(=O)OCCN2C(=O)c3ccccc3C2=O", "Cc1cccc(c1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "CC1=Nc2ccccc2C(=O)N1NC(=O)Nc3ccc(Cl)cc3", "CC1=Nc2ccccc2C(=O)N1NC(=O)Nc3ccc(Br)cc3", "CC1=Nc2ccccc2C(=O)N1NC(=O)Nc3ccc(I)cc3", "CC1=Nc2ccccc2C(=O)N1NC(=O)Nc3ccc(cc3)[N+](=O)[O-]", "CC1=Nc2ccccc2C(=O)N1NC(=O)Nc3ccccc3C", "CC1=Nc2ccccc2C(=O)N1NC(=O)Nc3cccc(C)c3", "CC1=Nc2ccccc2C(=O)N1NC(=O)Nc3ccc(C)cc3", "CCc1ccc(NC(=O)NN2C(=Nc3ccccc3C2=O)C)cc1", "COc1ccc(NC(=O)NN2C(=Nc3ccccc3C2=O)C)cc1", "CCOc1ccc(NC(=O)NN2C(=Nc3ccccc3C2=O)C)cc1", "CC1=Nc2ccccc2C(=O)N1NC(=O)Nc3ccc(F)cc3", "COc1ccc2c(c1)c(CC(=O)Oc3cc(O)c4C(=O)C[C@H](Oc4c3)c5ccc(O)cc5)c(C)n2C(=O)c6ccc(Cl)cc6", "COc1ccc2c(c1)c(CC(=O)Oc3cc(O)c4C(=O)C[C@H](Oc4c3)c5ccc(OC)c(O)c5)c(C)n2C(=O)c6ccc(Cl)cc6", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\CC(=O)O", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H](C)C(=O)O", "COc1ccc(cc1)C2=CS\\\\C(=N/NC(=O)c3ccc(O)cc3)\\\\N2c4ccc(OC)cc4OC", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H](CC(C)C)C(=O)O", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H]([C@@H](C)O)C(=O)O", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H](Cc3ccccc3)C(=O)O", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H](Cc3ccc(O)cc3)C(=O)O", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H](Cc3cnc[nH]3)C(=O)O", "CCCCCCCCCCCCCCC(=O)N[C@@H](CCC(=O)O)C(=O)N[C@H]1CNC(=O)[C@@H]2CCCN2C(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CC(=O)O)NC(=O)CNC(=O)[C@H](CC(=O)O)NC(=O)CNC(=O)[C@@H]3CCCCN3C1=O)[C@@H](C)O)[C@H](C)CC", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H](Cc3c[nH]c4ccccc34)C(=O)O", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H](CC(=O)O)C(=O)O", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H](CCC(=O)O)C(=O)O", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H](CCCN)C(=O)O", "COc1ccc2cc(ccc2c1)C(C)C(=O)OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)C(OC(=O)C)\\\\C=N\\\\[C@@H](Cc3ccc(O)c(O)c3)C(=O)O", "CCCCCCCCCCCCCC(CCCCCCCCCCCCC)C(=O)OC[C@H]1O[C@H](O[C@H]2O[C@H](COC(=O)C(CCCCCCCCCCCCC)CCCCCCCCCCCCC)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", "Cc1ccc(cc1)C2=CS\\\\C(=N/NC(=O)c3ccc(O)cc3)\\\\N2c4c(C)cccc4C", "COc1ccc(N2\\\\C(=N\\\\NC(=O)c3ccc(O)cc3)\\\\SC=C2c4ccc(Cl)cc4)c(OC)c1", "Oc1ccc(cc1)C(=O)N\\\\N=C\\\\2/SC=C(N2c3ccc(Cl)cc3Cl)c4ccc(Cl)cc4", "COc1ccc(cc1)C2=NN(C(=O)c3ccc(O)cc3)C(=S)S2", "CCOC(=O)C1=C(C)N(CC(O)COc2ccc(\\\\C=N\\\\C(=S)Nc3ccc(OC)cc3)cc2)C(=S)NC1c4cccc(c4)[N+](=O)[O-]", "Oc1ccc(cc1)C(=O)N2N=C(SC2=S)c3ccc(O)cc3", "Oc1ccc(cc1)C(=O)N2N=C(SC2=S)c3c(Cl)cccc3Cl", "CCOC(=O)C1=C(C)N(CC(O)COc2ccc(\\\\C=N\\\\C(=S)Nc3ccccc3OC)cc2)C(=S)NC1c4ccccc4[N+](=O)[O-]", "CCOC(=O)C1=C(C)N(CC(O)COc2ccc(\\\\C=N\\\\C(=S)Nc3ccccc3OC)cc2)C(=S)NC1c4cccc(c4)[N+](=O)[O-]", "CCOC(=O)C1=C(C)N(CC(O)COc2ccc(\\\\C=N\\\\C(=S)Nc3ccc(OC)cc3)cc2)C(=S)NC1c4ccccc4[N+](=O)[O-]", "COc1ccc2cc(\\\\C=N\\\\NC(=O)C3=C(O)c4ccccc4S(=O)(=O)N3)c(Cl)nc2c1", "Fc1ccc(cc1)N2CCN(CCCCN3C(=O)C4=C(SCCS4)C3=O)CC2", "O=C1N(CCCCN2CCN(CC2)c3ccccn3)C(=O)C4=C1SCCS4", "O=C1N(CCCN2CCN(CC2)c3ccccn3)C(=O)C4=C1SCCS4", "Fc1ccc(cc1)N2CCN(CCCN3C(=O)C4=C(SCCS4)C3=O)CC2", "COc1ccc2cc(ccc2c1)C(C)C(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC", "O=C1N(CCCCN2CCN(CC2)c3ccccc3)C(=O)C4=C1SCCS4", "O=C1OC(=O)C2=C1SCCS2", "O=C1NC(=O)C2=C1SCCS2", "CN(C)CC(O)CN1C(=O)C2=C(SCCS2)C1=O", "CCN(CC)CC(O)CN1C(=O)C2=C(SCCS2)C1=O", "CC(C)(C)NCC(O)CN1C(=O)C2=C(SCCS2)C1=O", "C\\\\C=C(/C)\\\\O[C@H]1[C@H](O)[C@]2(COC(=O)C)[C@H](O)C[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CCC(O[C@@H]6O[C@@H]([C@@H](O)[C@H](O[C@@H]7OC[C@H](O)[C@H](O)[C@H]7O[C@@H]8OC[C@@H](O)[C@H](O)[C@H]8O)[C@H]6O[C@]9(C)O[C@H](CO)[C@H](O)[C@H](O)[C@H]9O)C(=O)O)[C@@](C)(C=O)[C@@H]5CC[C@@]34C)[C@@H]2CC1(C)C", "OC(CN1CCCCC1)CN2C(=O)C3=C(SCCS3)C2=O", "CC1CCN(CC(O)CN2C(=O)C3=C(SCCS3)C2=O)CC1", "OC(CN1CCN(CC1)c2ccccc2)CN3C(=O)C4=C(SCCS4)C3=O", "OC(CN1CCN(CC1)c2ccc(F)cc2)CN3C(=O)C4=C(SCCS4)C3=O", "OC1=C(C(C2=C(O)c3ccccc3OC2=O)c4ccc5CC=CCc5c4)C(=O)Oc6ccccc16", "CCCCCCCCCCCCCCC(CCCCCCCCCCCCCC)C(=O)OC[C@H]1O[C@H](O[C@H]2O[C@H](COC(=O)C(CCCCCCCCCCCCCC)CCCCCCCCCCCCCC)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", "OC(CN1CCN(CC1)c2ccccn2)CN3C(=O)C4=C(SCCS4)C3=O", "COc1ccccc1N2CCN(CC(O)CN3C(=O)C4=C(SCCS4)C3=O)CC2", "OC1=C(C(C2=C(O)c3ccccc3OC2=O)c4ccc(OCc5ccccc5)c(OCc6ccccc6)c4)C(=O)Oc7ccccc17", "OC1=C(C(C2=C(O)c3ccc(O)cc3OC2=O)c4ccc5CC=CCc5c4)C(=O)Oc6cc(O)ccc16", "OC1=C(C(C2=C(O)c3ccc(O)cc3OC2=O)c4ccc(CCc5ccccc5)cc4)C(=O)Oc6cc(O)ccc16", "CCCCCCCCCCCCCCCC(CCCCCCCCCCCCCCC)C(=O)OC[C@H]1O[C@H](O[C@H]2O[C@H](COC(=O)C(CCCCCCCCCCCCCCC)CCCCCCCCCCCCCCC)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", "OC1=C(C(=O)Oc2ccccc12)c3ccccc3C4=C(O)c5ccccc5OC4=O", "OC1=C(C(=O)Oc2ccccc12)c3cccc(c3)C4=C(O)c5ccccc5OC4=O", "CCn1cc(C(=O)\\\\C=C(/O)\\\\C(=O)O)c2ccccc12", "CCn1cc(C(=O)\\\\C=C(/O)\\\\C(=O)O)c2cc3OCOc3cc12", "OC1=C(C(C2=C(O)c3ccccc3OC2=O)c4cnc5CC=CCc5c4)C(=O)Oc6ccccc16", "O=C(NN=C(c1ccccc1)c2ccccc2)c3ccncc3", "Brc1ccc(cc1)\\\\C(=N\\\\NC(=O)c2ccncc2)\\\\c3ccccc3", "Cc1ccc2cc(\\\\C=N\\\\NC(=O)C3=C(O)c4ccccc4S(=O)(=O)N3)c(Cl)nc2c1", "Cc1cccc2cc(\\\\C=N\\\\NC(=O)C3=C(O)c4ccccc4S(=O)(=O)N3)c(Cl)nc12", "COc1ccc2nc(Cl)c(\\\\C=N\\\\NC(=O)C3=C(O)c4ccccc4S(=O)(=O)N3)cc2c1", "Cc1ccc2nc(Cl)c(\\\\C=N\\\\NC(=O)C3=C(O)c4ccccc4S(=O)(=O)N3)cc2c1", "COC(=O)[C@]1(C)CCC[C@@]2(C)[C@H]1CC(=O)c3cc(c(Nc4ccc(OC)cc4)cc23)[N+](=O)[O-]", "CCCCCCCCCCCCCCCCC(CCCCCCCCCCCCCCCC)C(=O)OC[C@H]1O[C@H](O[C@H]2O[C@H](COC(=O)C(CCCCCCCCCCCCCCCC)CCCCCCCCCCCCCCCC)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", "COC(=O)[C@]1(C)CCC[C@@]2(C)[C@H]1CC(=O)c3cc(c(Nc4ccc(C)cc4)cc23)[N+](=O)[O-]", "COC(=O)[C@]1(C)CCC[C@@]2(C)[C@H]1CC(=O)c3cc(c(Nc4ccccc4)cc23)[N+](=O)[O-]", "COC(=O)[C@]1(C)CCC[C@@]2(C)[C@H]1CC(=O)c3cc(c(NCCN)cc23)[N+](=O)[O-]", "COC(=O)[C@]1(C)CCC[C@@]2(C)[C@H]1CC(=O)c3cc(c(NCCO)cc23)[N+](=O)[O-]", "COC(=O)[C@]1(C)CCC[C@@]2(C)[C@H]1CC(=O)c3cc(c(NC(C)(C)C)cc23)[N+](=O)[O-]", "CCCCNc1cc2c(cc1[N+](=O)[O-])C(=O)C[C@H]3[C@@](C)(CCC[C@]23C)C(=O)OC", "COC(=O)[C@]1(C)CCC[C@@]2(C)[C@H]1CC(=O)c3cc(c(NC(C)C)cc23)[N+](=O)[O-]", "CCCNc1cc2c(cc1[N+](=O)[O-])C(=O)C[C@H]3[C@@](C)(CCC[C@]23C)C(=O)OC", "CCNc1cc2c(cc1[N+](=O)[O-])C(=O)C[C@H]3[C@@](C)(CCC[C@]23C)C(=O)OC", "CNc1cc2c(cc1[N+](=O)[O-])C(=O)C[C@H]3[C@@](C)(CCC[C@]23C)C(=O)OC", "COc1cc2CCN(C)[C@H]3Cc4ccc(O)c(c4)c5cc(C[C@H]6N(C)CCc7cc(OC)c(Oc(c1O)c23)cc67)ccc5OC", "OC(COc1ccc(cc1)C(F)(F)F)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F", "COc1ccc2C[C@@H]3N(C)CCc4cc5Oc6c(OC)cc7CCN=C(Cc8ccc(Oc1c2)cc8)c7c6Oc5cc34", "C[C@H]1O[C@H](OC[C@H]2O[C@@H](Oc3cc(O)c4C(=O)C(=C(Oc4c3)c5ccc(O)cc5)O)[C@H](O[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O", "C\\\\C=C\\\\1/CN2CC[C@@]34[C@@H]2C[C@@H]1[C@H](CO)[C@@H]3N(C(=O)C)c5ccccc45", "COC(=O)c1ccc(OCC(O)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F)cc1", "COc1cccc(OCC(O)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F)c1", "COc1ccc(OCC(O)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F)cc1", "CC(=O)c1ccc(OCC(O)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F)cc1", "CCC(=O)c1ccc(OCC(O)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F)cc1", "CC(=O)Nc1ccc(OCC(O)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F)cc1", "OC(COc1ccc(cc1)C#N)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F", "CCCNCC#CC(=O)Nc1ccc2ncc(C#N)c(Nc3cccc(Br)c3)c2c1", "O=C1NC(=O)C(=CN2CCOCC2)C(=O)N1", "CC(C)NCC#CC(=O)Nc1ccc2ncc(C#N)c(Nc3cccc(Br)c3)c2c1", "CCC1(CC)CC(CCNC(=O)c2ccccc2)OC1=O", "CCC1(CC)CC(CCOC(=O)c2ccc(Br)cc2)OC1=O", "CCC(=O)NCCC1CC(CC)(CC)C(=O)O1", "CCC1(CC)CC(CCOC(=O)c2ccc(OC)cc2)OC1=O", "CCCCCCCCCCCCCCC[C@@H](O)[C@@H](CCCCCCCCCCCCCC)C(=O)OC[C@H]1O[C@H](O[C@H]2O[C@H](COC(=O)[C@H](CCCCCCCCCCCCCC)[C@H](O)CCCCCCCCCCCCCCC)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", "CCC1(CC)CC(CCOC(=O)c2ccc(C)cc2)OC1=O", "CCC1(CC)CC(CCOC(=O)c2ccccc2)OC1=O", "CCC(=O)OCCC1CC(CC)(CC)C(=O)O1", "CCCCCCCC(=O)OCCC1CC(CC)(CC)C(=O)O1", "CCC1(CC)CC(COC(=O)c2ccc(Br)cc2)OC1=O", "CC(C)Cc1ccc(cc1)C(C)C2=NNC(=S)N2c3ccc(cc3)C(=O)NNC(=O)CO[N+](=O)[O-]", "CCC1(CC)CC(COC(=O)c2ccc(C)cc2)OC1=O", "CCC1(CC)CC(COC(=O)c2ccc(OC)cc2)OC1=O", "COc1ccc(cc1OC)C2OC(=O)C(=C2I)O", "CCC1(CC)CC(COC(=O)c2ccccc2)OC1=O", "CCCCCCCC(=O)OCC1CC(CC)(CC)C(=O)O1", "CCC(=O)OCC1CC(CC)(CC)C(=O)O1", "Nc1nc(Cl)c2c(ncn2Cc3cccc(CCl)c3)n1", "CCOC(=O)N1C(CC(=O)C)N(C(=O)OCC)c2ccccc12", "ClCc1ccc(Cn2cnc3ncnc(Cl)c23)cc1", "CCOC(=O)N1C(CC(=O)C)N(C(=O)OCC)c2cc(C)c(C)cc12", "CCOC(=O)N1C(CC(=O)c2ccccc2)N(C(=O)OCC)c3ccccc13", "CCOC(=O)N1C(CC(=O)c2ccccc2)N(C(=O)OCC)c3cc(C)c(C)cc13", "COc1ccc(cc1)C2C\\\\C(=N\\\\OCc3ccccc3)\\\\C(C)C(N2C)c4ccc(OC)cc4", "ClCc1ccc(Cn2cnc3nc(Cl)nc(Cl)c23)cc1", "CCOC(=O)N1C(CC(=O)c2ccc(O)cc2)N(C(=O)OCC)c3ccccc13", "CCOC(=O)N1C(CC(=O)c2ccc(OC)cc2)N(C(=O)OCC)c3cc(C)c(C)cc13", "CCOC(=O)N1C(CC(=O)Cc2ccccc2)N(C(=O)OCC)c3ccccc13", "CCOC(=O)N1C(CC(=O)\\\\C=C\\\\c2ccccc2)N(C(=O)OCC)c3cc(C)c(C)cc13", "ClCc1cccc(Cn2cnc3nc(Cl)nc(Cl)c23)c1", "C[C@H]1O[C@H](OC2=C(Oc3cc(O)cc(O)c3C2=O)c4ccc(O)c(O)c4)[C@@H](O)[C@@H](O)[C@@H]1O", "Nc1nc(Cl)c2c(ncn2Cc3ccc(CCl)cc3)n1", "ClCc1ccccc1Cn2cnc3c(Cl)ncnc23", "ClCc1ccccc1Cn2cnc3c(Cl)nc(Cl)nc23", "ClCc1cccc(Cn2cnc3c(Cl)ncnc23)c1", "ClCc1cccc(Cn2cnc3c(Cl)nc(Cl)nc23)c1", "ClCc1ccc(Cn2cnc3c(Cl)ncnc23)cc1", "ClCc1ccc(Cn2cnc3c(Cl)nc(Cl)nc23)cc1", "ClCc1ccccc1Cn2cnc3ncnc(Cl)c23", "ClCc1ccccc1Cn2cnc3nc(Cl)nc(Cl)c23", "Nc1nc(Cl)c2c(ncn2Cc3ccccc3CCl)n1", "ClCc1cccc(Cn2cnc3ncnc(Cl)c23)c1", "Cc1ccc(cc1)c2cc(nc(N)n2)c3ccccc3O", "Nc1nc(cc(n1)c2ccccc2O)c3ccc(F)cc3", "COc1ccc(cc1)c2cc(nc(N)n2)c3ccccc3O", "Nc1nc(cc(n1)c2ccccc2O)c3ccc(Cl)cc3", "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@H](CC[C@]34C)OC(=O)Cc5ccccc5)[C@@H]1CC/C/2=N\\\\O", "Cc1cccc(c1)c2cc(nc(N)n2)c3ccccc3O", "CCn1c(cc2ccccc12)C(=O)\\\\C=C(/O)\\\\C(=O)O", "COc1cccc(c1)c2cc(nc(N)n2)c3ccccc3O", "Nc1nc(cc(n1)c2ccccc2O)c3cccc(Cl)c3", "COc1ccc(cc1OC)c2cc(nc(N)n2)c3ccccc3O", "Nc1nc(cc(n1)c2ccccc2O)c3occc3", "Cc1noc(n1)c2ccccn2", "Nc1nc(cc(n1)c2ccccc2O)c3cccs3", "CC1C(N(C)C(C/C/1=N/OCc2ccccc2)c3ccccc3)c4ccccc4", "CCC1C(N(C)C(C/C/1=N/OCc2ccccc2)c3ccccc3)c4ccccc4", "CC(C)C1C(N(C)C(C/C/1=N/OCc2ccccc2)c3ccccc3)c4ccccc4", "o1ncnc1c2cccnc2", "CC(C)CNCC(C)(C)N1C(=O)c2ccccc2C1=O", "CN1C(CC(=NOCc2ccccc2)CC1c3ccc(Cl)cc3)c4ccc(Cl)cc4", "CC1C(N(C)C(C/C/1=N/OCc2ccccc2)c3ccc(Cl)cc3)c4ccc(Cl)cc4", "CN1C(CC(=NOCc2ccccc2)CC1c3ccc(C)cc3)c4ccc(C)cc4", "CC1C(N(C)C(C/C/1=N/OCc2ccccc2)c3ccc(C)cc3)c4ccc(C)cc4", "COc1ccc(cc1)C2CC(=NOCc3ccccc3)CC(N2C)c4ccc(OC)cc4", "COc1ccc(cc1)C2CC(=NN2c3ccccc3)c4ccc(O)c(C)c4", "Cc1cc(ccc1O)C2=NN(C(C2)c3ccc(Cl)cc3)c4ccccc4", "CN(C)c1ccc(cc1)C2CC(=NN2c3ccccc3)c4ccc(O)c(C)c4", "Cc1cc(ccc1O)C2=NN(C(C2)c3ccccc3)c4ccccc4", "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@H](CC[C@]34C)OC(=O)COc5ccccc5)[C@@H]1CC/C/2=N\\\\O", "COc1ccc(cc1OC)C2CC(=NN2c3ccccc3)c4ccc(O)c(C)c4", "COc1cc(cc(OC)c1OC)C2CC(=NN2c3ccccc3)c4ccc(O)c(C)c4", "Cc1cc(ccc1O)C2=NN(C(C2)c3ccc(F)cc3)c4ccccc4", "Cc1cc(ccc1O)C2=NN(C(C2)c3ccccc3Cl)c4ccccc4", "C\\\\C=C(/C)\\\\O[C@H]1[C@H](CO)[C@]2(C)[C@H](O)C[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CCC(O[C@@H]6O[C@@H]([C@@H](O)[C@H](O[C@@H]7OC[C@H](O)[C@H](O)[C@H]7O[C@@H]8OC[C@@H](O)[C@H](O)[C@H]8O)[C@H]6O[C@]9(C)O[C@H](CO)[C@H](O)[C@H](O)[C@H]9O)C(=O)O)[C@@](C)(C=O)[C@@H]5CC[C@@]34C)[C@@H]2CC1(C)C", "Cc1cc(ccc1O)C2=NN(C(C2)c3c(Cl)cccc3Cl)c4ccccc4", "Cc1cc(ccc1O)C2=NN(C(C2)c3cccc(c3)[N+](=O)[O-])c4ccccc4", "Cc1cc(ccc1O)C2=NN(C(C2)c3occc3)c4ccccc4", "Clc1ccc(cc1)c2cc(C3=CC(=O)C=CC3=O)n(n2)c4ccccc4", "CC(C)(CO)NC=C1C(=O)NC(=O)NC1=O", "Cc1ccc(cc1)c2cc(C3=CC(=O)C=CC3=O)n(n2)c4ccc(cc4)S(=O)(=O)N", "NS(=O)(=O)c1ccc(cc1)n2nc(cc2C3=CC(=O)C=CC3=O)c4ccc(Cl)cc4", "CC(CO)(CO)Nc1nc2ccccc2[nH]1", "OCCN(Cc1nc2ccccc2[nH]1)c3nc4ccccc4[nH]3", "CC(CO)(CO)NCc1nc2ccccc2[nH]1", "CCOC(=O)N1CCN(CC1)C=C2C(=O)NC(=O)NC2=O", "OC(CNC1CCCCC1)CN2C(=O)c3ccccc3C2=O", "CC(CO)(CO)NCCCN1C(=O)c2ccccc2C1=O", "CC(C)(CO)N1C(=O)c2ccccc2C1=O", "CC(C)(CO)NCC(C)(C)N1C(=O)c2ccccc2C1=O", "COc1ccc(cc1)C(=O)O[C@H]2CC[C@]3(C)[C@H]4CC[C@@]5(C)[C@@H](CC/C/5=N\\\\O)[C@@H]4CC=C3C2", "Cc1c(CCc2ccccc2)oc3cccc(OCCNCc4cccnc4)c13", "Cc1c(oc2cccc(OCCNCc3cccnc3)c12)C(=O)Nc4ccccc4", "Cc1c(CSc2ccccc2)oc3cccc(OCCNCc4cccnc4)c13", "OC1=C(Br)C(OC1=O)c2ccccc2", "CCOC(=O)C1=C(C)N(CC(O)COc2ccc(\\\\C=N\\\\C(=S)Nc3ccc(Cl)cc3)cc2)C(=S)NC1c4ccccc4[N+](=O)[O-]", "CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4", "CCOC(=O)c1oc2cccc(OCCNCc3ccccc3)c2c1C", "CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C", "CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1CC", "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@H](CC[C@]34C)OC(=O)c5ccc(Cl)cc5)[C@@H]1CC/C/2=N\\\\O", "CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C(C)C", "CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1", "CC(=O)c1nccn1c2oc3cccc(OCCNCc4cccnc4)c3c2C", "COCc1ccc2oc(cc2c1)c3oc4cccc(OCCNCc5cccnc5)c4c3C", "CCCCc1nc(Cl)c(C=O)n1CC(=O)c2ccccc2", "Cc1c(COc2ccccc2C#N)oc3cccc(OCCNCc4cccnc4)c13", "Cc1c(COc2ccccc2)oc3cccc(OCCNCc4cccnc4)c13", "Cc1c(COc2cccc(F)c2)oc3cccc(OCCNCc4cccnc4)c13", "Cc1c(COc2cccc(F)c2F)oc3cccc(OCCNCc4cccnc4)c13", "Cc1c(COc2ccccc2F)oc3cccc(OCCNCc4cccnc4)c13", "Cc1c(COc2ccc(F)c(F)c2F)oc3cccc(OCCNCc4cccnc4)c13", "Cc1c(COc2ccc(F)cc2F)oc3cccc(OCCNCc4cccnc4)c13", "Cc1c(COc2ccc(F)cc2)oc3cccc(OCCNCc4cccnc4)c13", "CCCCCCCCCCCCCCCC(=O)OCCCCCCOC(=O)CCCCCCCCCCCCCCC", "Cc1c(COc2ccc(Cl)cc2)oc3cccc(OCCNCc4cccnc4)c13", "Cc1c(COc2ccc(cc2)C#N)oc3cccc(OCCNCc4cccnc4)c13", "CCCCc1nc(Cl)c(C=O)n1Cc2c(C)onc2C", "CCCCc1nc(Cl)c(C=O)n1CC(=O)c2ccc(Cl)cc2", "CC1=C(C(C(=C(C)N1)C(=O)Nc2ccc(Cl)cc2)c3ccncc3)C(=O)Nc4ccc(Cl)cc4", "CCCCc1nc(Cl)c(C=O)n1Cc2ccc(cc2)c3ccccc3C(=O)OC", "CCCCc1nc(Cl)c(C=O)n1Cc2ccc(cc2)c3ccccc3C#N", "CCCCc1nc(Cl)c(C=O)n1Cc2ccc(cc2)[N+](=O)[O-]", "CCCCc1nc(Cl)c(C=O)n1Cc2cc3OCOc3cc2C", "CCCCc1nc(Cl)c(C=O)n1Cc2ccccc2", "CCCCc1nc(Cl)c(C=O)n1Cc2cc(OC)c(OC)cc2Br", "Fc1cccc(c1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "Clc1cccc(c1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "[O-][N+](=O)c1cccc(c1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "CC(C)(C)c1ccc(cc1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "CCCCCc1ccc(cc1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "Brc1ccc(cc1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "Fc1cccc(C(=O)OCCN2C(=O)c3ccccc3C2=O)c1F", "Fc1ccc(C(=O)OCCN2C(=O)c3ccccc3C2=O)c(F)c1", "Fc1ccc(F)c(c1)C(=O)OCCN2C(=O)c3ccccc3C2=O", "Fc1cccc(F)c1C(=O)OCCN2C(=O)c3ccccc3C2=O", "Fc1ccc(cc1F)C(=O)OCCN2C(=O)c3ccccc3C2=O", "Clc1ccccc1C(=O)OCCN2C(=O)c3ccccc3C2=O", "[O-][N+](=O)c1ccccc1C(=O)OCCN2C(=O)c3ccccc3C2=O", "CCCCCCCCCCCCCCCC(=O)OCC(CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC", "COc1ccc(cc1)C2=NOC(C2CN3CCOCC3)c4c[nH]c5ccccc45", "Cc1ccc(cc1)C2=NOC(C2CN3CCOCC3)c4c[nH]c5ccccc45", "Clc1ccc(cc1)C2=NOC(C2CN3CCOCC3)c4c[nH]c5ccccc45", "C(C1C(ON=C1c2ccccc2)c3c[nH]c4ccccc34)N5CCOCC5", "CCCCCCCCCCCCNCC(=O)N[C@@H](CC(=O)N)C(=O)N[C@H]1CNC(=O)[C@@H]2CCCN2C(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CC(=O)O)NC(=O)CNC(=O)[C@H](CC(=O)O)NC(=O)CNC(=O)[C@@H]3CCCCN3C1=O)[C@@H](C)O)[C@H](C)CC", "COc1ccc(cc1)C2=NOC(C2CN3CCCCC3)c4c[nH]c5ccccc45", "Cc1ccc(cc1)C2=NOC(C2CN3CCCCC3)c4c[nH]c5ccccc45", "Clc1ccc(cc1)C2=NOC(C2CN3CCCCC3)c4c[nH]c5ccccc45", "C(C1C(ON=C1c2ccccc2)c3c[nH]c4ccccc34)N5CCCCC5", "CCCCCCCCCCCCCCC[C@@H](OC)[C@@H](CCCCCCCCCCCCCC)C(=O)OC[C@H]1O[C@H](O[C@H]2O[C@H](COC(=O)[C@H](CCCCCCCCCCCCCC)[C@@H](CCCCCCCCCCCCCCC)OC)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", "COc1ccc(NC(=O)C2=C(C)NC(=O)NC2c3ccccc3)cc1", "COc1ccc(NC(=O)C2=C(C)NC(=S)NC2c3ccccc3)cc1", "COc1ccc(NC(=O)C2=C(C)NC(=O)NC2c3ccc(OC)cc3)cc1", "COc1ccc(NC(=O)C2=C(C)NC(=S)NC2c3ccc(OC)cc3)cc1", "Cc1ccc(cc1)C(=O)O[C@H]2CC[C@]3(C)[C@H]4CC[C@@]5(C)[C@@H](CC/C/5=N\\\\O)[C@@H]4CC=C3C2", "COc1ccc(NC(=O)C2=C(C)NC(=O)NC2c3occc3)cc1", "COc1ccc(NC(=O)C2=C(C)NC(=S)NC2c3occc3)cc1", "COc1ccc(NC(=O)C2=C(C)NC(=O)NC2c3ccccc3O)cc1", "COc1ccc(NC(=O)C2=C(C)NC(=S)NC2c3ccccc3O)cc1", "CCOC(=O)C1=C(C)N(CC(O)COc2ccc(\\\\C=N\\\\C(=S)Nc3ccc(cc3)C(=O)O)cc2)C(=S)NC1c4cccc(c4)[N+](=O)[O-]", "COc1ccc(NC(=O)C2=C(C)NC(=O)NC2C\\\\C=C\\\\c3ccccc3)cc1", "COc1ccc(NC(=O)C2=C(C)NC(=O)NC2c3cccc(c3)[N+](=O)[O-])cc1", "COc1ccc(NC(=O)C2=C(C)NC(=S)NC2c3cccc(c3)[N+](=O)[O-])cc1", "COc1ccc(NC(=O)C2=C(C)NC(=O)NC2c3ccc(O)c(OC)c3)cc1", "CCCCCCCCCCCCCCCCNC(=O)N[C@@H](CC(=O)N)C(=O)N[C@H]1CNC(=O)[C@@H]2CCCN2C(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CC(=O)O)NC(=O)CNC(=O)[C@H](CC(=O)O)NC(=O)CNC(=O)[C@@H]3CCCCN3C1=O)[C@@H](C)O)[C@H](C)CC", "COc1ccc(NC(=O)C2=C(C)NC(=O)NC2c3cccc(Cl)c3)cc1", "COc1ccc(NC(=O)C2=C(C)NC(=S)NC2c3cccc(Cl)c3)cc1", "COc1ccc(NC(=O)C2=C(C)NC(=O)NC2c3ccc(cc3)N(C)C)cc1", "COc1ccc(NC(=O)C2=C(C)NC(=S)NC2c3ccc(cc3)N(C)C)cc1", "CN1c2c(C)n(nc2c3ccccc3S1(=O)=O)c4ccc(cc4)C(=O)\\\\C=C\\\\c5ccc(Cl)cc5", "COc1ccc(NC(=O)C2=C(C)NC(=O)NC2c3ccncc3)cc1", "COc1ccc(NC(=O)C2=C(C)NC(=O)NC2)cc1", "CC1=C(C(C(=C(C)N1)C(=O)Nc2ccc(Cl)cc2)c3ccccc3)C(=O)Nc4ccc(Cl)cc4", "COc1ccc(cc1)C2C(=C(C)NC(=C2C(=O)Nc3ccc(Cl)cc3)C)C(=O)Nc4ccc(Cl)cc4", "CC1=C(C(C(=C(C)N1)C(=O)Nc2ccc(Cl)cc2)c3occc3)C(=O)Nc4ccc(Cl)cc4", "COc1cc(ccc1O)C2C(=C(C)NC(=C2C(=O)Nc3ccc(Cl)cc3)C)C(=O)Nc4ccc(Cl)cc4", "CC[C@@]1(C[C@H]2CN(CCc3c([nH]c4ccccc34)[C@@](C2)(C(=O)OC)c5cc6c(cc5OC)N(C)[C@H]7[C@](O)([C@H](OC(=O)C)[C@]8(CC)C=CCN9CC[C@]67[C@H]89)C(=O)OC)C1)NC(=S)Nc%10ccccc%10", "CCOC(=O)C1=C(C)N(CC(O)COc2ccc(\\\\C=N\\\\C(=S)Nc3ccc(Cl)cc3)cc2)C(=S)NC1c4cccc(c4)[N+](=O)[O-]", "Cc1ccc2[nH]c3c(CCc4c3nc5ccc(Cl)cc5c4c6ccccc6)c2c1", "Cc1cccc2c3CCc4c(nc5ccc(Cl)cc5c4c6ccccc6)c3[nH]c12", "CC(C)Cc1ccc(cc1)C(C)C2=NNC(=S)N2c3ccc(cc3)C(=O)NNC(=O)CCO[N+](=O)[O-]", "OC1=C(Oc2c(O)ccc(O)c2C1=O)c3ccc(O)c(O)c3", "Cc1ccc2c3CCc4c(nc5ccc(Cl)cc5c4c6ccccc6)c3[nH]c2c1", "CC(C)Cc1ccc(cc1)C(C)C2=NNC(=S)N2c3ccc(C(=O)NNC(=O)CCO[N+](=O)[O-])c(c3)C(=O)NNC(=O)CCO[N+](=O)[O-]", "CC(C)Cc1ccc(cc1)C(C)C2=NNC(=S)N2c3ccccc3CC(=O)NNC(=O)CCO[N+](=O)[O-]", "CC(C)Cc1ccc(cc1)C(C)C2=NNC(=S)N2c3ccccc3CC(=O)NNC(=O)CO[N+](=O)[O-]", "CC(C)Cc1ccc(cc1)C(C)C2=NNC(=S)N2c3ccc(C(=O)NNC(=O)CO[N+](=O)[O-])c(c3)C(=O)NNC(=O)CO[N+](=O)[O-]", "CC[C@@]1(C[C@H]2CN(CCc3c([nH]c4ccccc34)[C@@](C2)(C(=O)OC)c5cc6c(cc5OC)N(C)[C@H]7[C@](O)([C@H](OC(=O)C)[C@]8(CC)C=CCN9CC[C@]67[C@H]89)C(=O)OC)C1)NC(=S)NCCc%10ccc(F)cc%10", "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@H](CC[C@]34C)OC(=O)CCl)[C@@H]1CC/C/2=N\\\\O", "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@H](CC[C@]34C)OC(=O)c5ccccc5)[C@@H]1CC/C/2=N\\\\O", "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@H](CC[C@]34C)OC(=O)c5ccc(cc5)[N+](=O)[O-])[C@@H]1CC/C/2=N\\\\O", "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@H](CC[C@]34C)OC(=O)c5ccc(N)cc5)[C@@H]1CC/C/2=N\\\\O", "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@H](CC[C@]34C)OC(=O)c5ccc(O)cc5)[C@@H]1CC/C/2=N\\\\O", "C[C@]12CC[C@H](O)CC1=CC[C@H]3[C@@H]4CC\\\\C(=N/O)\\\\[C@@]4(C)CC[C@H]23", "COc1ccc(cc1)c2nn(cc2C3OC(=O)C(=C3Br)O)c4ccccc4", "OC1=C(Br)C(OC1=O)c2cn(nc2c3ccc(Br)cc3)c4ccccc4", "Cc1nn(c(C)c1C2OC(=O)C(=C2Cl)O)c3ccccc3", "OC1=C(I)C(OC1=O)c2ccccc2", "OC1=C(I)C(OC1=O)c2ccc(Cl)cc2", "Oc1ccc(cn1)c2ccc3N=C(NCCN4CCOCC4)C(=O)N(CC5CCCCC5)c3n2", "Cc1nn(c(C)c1C2OC(=O)C(=C2I)O)c3ccccc3", "BrC1=C(OCC(=O)NCc2ccccc2)C(=O)OC1c3ccccc3", "Cl\\\\C(=C/c1ccc(Cl)cc1)\\\\C2=Nc3ccccc3NC2=O", "I\\\\C(=C/c1ccccc1)\\\\C2=Nc3ccccc3NC2=O", "Clc1ccc(\\\\C=C(/I)\\\\C2=Nc3ccccc3NC2=O)cc1", "Cc1nn(c(C)c1\\\\C=C(/I)\\\\C2=Nc3ccccc3NC2=O)c4ccccc4", "OC(=O)\\\\C(=C\\\\c1ccc(Cl)cc1)\\\\Br", "Cc1nn(c(C)c1\\\\C=C(/Br)\\\\C(=O)O)c2ccccc2", "OC(=O)\\\\C(=C\\\\c1cn(nc1c2ccc(Br)cc2)c3ccccc3)\\\\Br", "CN(CCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2)CC(O)COc3ccccc3C", "COc1ccccc1OCC(O)CN(C)CCC(Oc2ccc(cc2)C(F)(F)F)c3ccccc3", "OC(=O)C(=O)\\\\C=C\\\\c1cn(nc1c2ccccc2)c3ccccc3", "COc1ccc(cc1)c2nn(cc2\\\\C=C\\\\C(=O)C(=O)O)c3ccccc3", "OC(=O)C(=O)\\\\C=C\\\\c1cn(nc1c2ccc(Br)cc2)c3ccccc3", "Cn1c(cc2ccccc12)C(=O)\\\\C=C(/O)\\\\C(=O)O", "Cn1c(cc2cc3OCOc3cc12)C(=O)\\\\C=C(/O)\\\\C(=O)O", "CCn1c(cc2cc3OCOc3cc12)C(=O)\\\\C=C(/O)\\\\C(=O)O", "OC(=O)\\\\C(=C\\\\C(=O)c1cc2ccccc2n1Cc3ccccc3)\\\\O", "OC(=O)\\\\C(=C\\\\C(=O)c1cc2cc3OCOc3cc2n1Cc4ccccc4)\\\\O", "Cn1cc(C(=O)\\\\C=C(/O)\\\\C(=O)O)c2ccccc12", "Cn1cc(C(=O)\\\\C=C(/O)\\\\C(=O)O)c2cc3OCOc3cc12", "COc1cccc2cc(C3SC(=NN3C(=O)C)NC(=O)C)c(Cl)nc12", "COc1ccc2cc(C3SC(=NN3C(=O)C)NC(=O)C)c(Cl)nc2c1", "COc1cc2cc(C3SC(=NN3C(=O)C)NC(=O)C)c(Cl)nc2cc1OC", "CC(=O)NC1=NN(C(S1)c2cc3ccc(Cl)cc3nc2Cl)C(=O)C", "CCC(=O)c1ccc(OCC(O)CN(C)CCC(Oc2ccc(cc2)C(F)(F)F)c3ccccc3)cc1", "CC(=O)NC1=NN(C(S1)c2cc3cc(Br)ccc3nc2Cl)C(=O)C", "CC(=O)NC1=NN(C(S1)c2cc3cc(C)ccc3nc2Cl)C(=O)C", "COc1ccc2nc(Cl)c(cc2c1)C3SC(=NN3C(=O)C)NC(=O)C", "CC(=O)NC1=NN(C(S1)c2cc3ccc(C)cc3nc2Cl)C(=O)C", "Oc1ncc(cn1)c2ccc3N=C(NCC4CCOCC4)C(=O)N(CC5CCCCC5)c3n2", "C\\\\C=C\\\\1/CN2CC[C@@]34[C@@H]2C[C@@H]1[C@H]5CO[C@H]([C@@H]6[C@H]7C[C@@H]8N(CC[C@]89[C@H]6N(C(=O)C)c%10ccccc9%10)C/C/7=C\\\\C)N([C@H]35)c%11ccccc4%11", "CCc1ccc2c3[nH]c4ccccc4c3CC[n+]2c1", "CN1CCc2c([nH]c3ccccc23)[C@@H]1C[C@@H]4C[C@@H]5N(CCc6c5[nH]c7ccccc67)C[C@@H]4C=C.OC(=O)C(=O)O", "COc1cc2CCN(C)[C@@H]3Cc4ccc(O)c(c4)c5cc(C[C@@H]6N(C)CCc7cc8Oc1c(Oc8cc67)c23)ccc5OC", "C\\\\C=C/1\\\\CN2CCc3c([C@@H]2C[C@@H]1Cc4nccc5c6ccccc6[nH]c45)n(C)c7ccccc37", "CN(CCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2)CC(O)COc3ccc(cc3)C(F)(F)F", "O[C@@]12CCO[C@@]1(Nc3ccccc23)[C@H]4C[C@@H]5CCN4C[C@@H]5C=C", "O[C@H]1Cc2c(O)cc3O[C@@]4(Oc5cc(O)cc(O)c5[C@@H]([C@H]4O)c3c2O[C@@H]1c6cc(O)c(O)c(O)c6)c7ccc(O)cc7", "C\\\\C=C/1\\\\CN(C)CC[C@@]23[C@@H]4[C@H](CO[C@H](O)C(=O)N4c5ccccc25)[C@H]1CC3=O", "C\\\\C=C/1\\\\CN2CC[C@@]34[C@@H]2C[C@@H]1[C@@H](C=O)[C@@H]3N(C(=O)C)c5ccccc45", "COC(=O)C1=CO[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H]3[C@@H]1C=C[C@]34OC(=O)C(=C4)C(O)c5ccc(O)c(OC)c5", "CC[C@@H]1CN2CC[C@@]34[C@@H]2C[C@@H]1[C@H]5CO[C@H]([C@@H]6[C@H]7C[C@@H]8N(CC[C@]89[C@H]6N(C(=O)C)c%10ccccc9%10)C/C/7=C/C)N([C@H]35)c%11ccccc4%11", "COc1cccc2[nH]c3[C@H](C[C@H]4C[C@@H]5N(CCc6c5[nH]c7ccccc67)C[C@@H]4C=C)NCCc3c12", "C[C@@H]1NC[C@]23CC[C@H]4[C@@H](CC=C5C[C@H](CC[C@]45C)N(C)C)[C@@H]2CC[C@H]13", "C[C@@H]1N=C[C@]23CC[C@H]4[C@@H](CC=C5C[C@@H](N)CC[C@]45C)[C@@H]2CC[C@H]13", "COCCc1c([nH]c2ccccc12)[C@@H]3C[C@@H]4CCN3C[C@@H]4C=C", "COc1ccc2nccc(C(=O)[C@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1", "C\\\\C=C/1\\\\CN2CC[C@@]34[C@@H]2C[C@@H]1[C@H](C=O)[C@@H]3N(C(=O)C)c5ccccc45", "CN1CCc2c([nH]c3ccccc23)[C@@H]1C[C@H]4C[C@@H]5N(CC[C@@]56C(=O)Nc7cc(O)ccc67)C[C@@H]4C=C", "COc1cc2CCN(C)[C@H]3Cc4ccc(O)c(Oc5ccc(C[C@@H]6N(C)CCc7c(O)c(OC)c(OC)c(Oc1cc23)c67)cc5)c4", "COC(=O)C(CO)C1CC2c3[nH]c4ccccc4c3CC[N+]2(C)C/C/1=C/C", "O[C@H]1Cc2c(O)cc3O[C@@]4(Oc5c([C@@H]([C@H]4O)c3c2O[C@@H]1c6ccc(O)cc6)c(O)cc7O[C@@]8(Oc9cc(O)cc(O)c9[C@@H]([C@@H]8O)c57)c%10ccc(O)cc%10)c%11ccc(O)cc%11", "C[C@@H]1N=C[C@]23C[C@@H](O)[C@H]4[C@@H](CCC5=CC(=O)C=C[C@]45C)[C@@H]2CC[C@H]13", "CO[C@H]1OC[C@@H]2[C@H]3CC(=O)[C@]4(CCN(C)C/C/3=C\\\\C)[C@H]2N(C1=O)c5ccccc45", "CC[C@@H]1CN2CC[C@@]34[C@@H]2C[C@@H]1[C@H]5CO[C@H]([C@@H]6[C@H]7C[C@@H]8N(CC[C@]89[C@H]6N(C(=O)C)c%10c(O)cccc9%10)C[C@H]7CC)N([C@H]35)c%11ccccc4%11", "COc1cc2CCN(C)[C@@H]3Cc4ccc5Oc6c(Oc5c4)c(OC)cc7CCN(C)[C@@H](Cc8ccc(Oc(c1OC)c23)cc8)c67", "CN(CCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2)CC(O)COc3ccc(cc3)C(=O)C", "CN(CCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2)CC(O)COc3cccc(c3)C(F)(F)F", "CN(CCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2)CC(O)COc3ccc(NC(=O)C)cc3", "CN(CCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2)CC(O)COc3ccc(Br)cc3", "COc1ccc(OCC(O)CN(C)CCC(Oc2ccc(cc2)C(F)(F)F)c3ccccc3)cc1", "CN(CCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2)CC(O)COc3ccc(C=O)cc3", "CN(CCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2)CC(O)COc3ccc(Cl)cc3", "CCNc1ccc(OCC(O)CN(C)CCC(Oc2ccc(cc2)C(F)(F)F)c3ccccc3)cc1", "CN(CCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2)CC(O)COc3ccc4ccccc4c3", "CN(CCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2)CC(O)COc3ccc(cc3)C#N", "CN(CCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2)CC(O)COc3ccccc3", "CC(=O)NS(=O)(=O)c1ccc(cc1)n2nc(cc2C3=CC(=O)C=CC3=O)c4ccccc4", "Cc1nn(c(C)c1CC(=O)c2cc(O)ccc2O)c3ccc(Br)cc3", "COc1ccc(OC)c(c1)c2cc(nn2c3ccc(cc3)S(=O)(=O)N)c4ccc(Cl)cc4", "CC1=NN(C(=O)C1CC(=O)c2cc(O)ccc2O)c3ccc(Br)cc3", "Oc1ccc(O)c(c1)C(=O)CC2C(=O)NN(C2=O)c3ccc(Br)cc3", "COc1ccc(OC)c(c1)c2cc(nn2c3ccccc3)c4ccc(Cl)cc4", "COc1ccc(OC)c(c1)c2cc(nn2c3ccc(cc3)S(=O)(=O)N)c4ccccc4", "COc1ccc(OC)c(c1)c2cc(nn2c3ccc(cc3)S(=O)(=O)NC(=O)C)c4ccccc4", "COc1ccc(OC)c(c1)c2cc(nn2c3ccc(cc3)S(=O)(=O)NC(=O)C)c4ccc(Cl)cc4", "Oc1ccc(O)c(c1)c2cc(nn2c3ccccc3)c4ccccc4", "OC1(C(=O)Nc2ccc(Cl)cc12)C(F)(F)F", "Oc1ccc(O)c(c1)c2cc(nn2c3ccccc3)c4ccc(Cl)cc4", "CC(=O)N1C(=O)C(O)(c2cc(Cl)ccc12)C(F)(F)F", "NS(=O)(=O)c1ccc(cc1)n2nc(cc2c3cc(O)ccc3O)c4ccccc4", "NS(=O)(=O)c1ccc(cc1)n2nc(cc2c3cc(O)ccc3O)c4ccc(Cl)cc4", "CC(=O)NS(=O)(=O)c1ccc(cc1)n2nc(cc2c3cc(O)ccc3O)c4ccccc4", "CC(=O)NS(=O)(=O)c1ccc(cc1)n2nc(cc2c3cc(O)ccc3O)c4ccc(Cl)cc4", "CN1C(=O)C(O)(c2cc(Br)ccc12)C(F)(F)F", "CC(=O)N1C(=O)C(O)(c2cc(C)ccc12)C(F)(F)F", "[O-][N+](=O)c1cc2oc(nc2cc1Cl)C3CCCCC3", "CCOc1ccc(Cc2oc3ccc(cc3n2)[N+](=O)[O-])cc1", "Nc1ccc(cc1)c2oc3ccc(Cl)cc3n2", "Cc1ccc2[nH]c(Cc3ccc(Br)cc3)nc2c1", "Clc1ccc2[nH]c(CC3CCCCC3)nc2c1", "Fc1ccc(cc1)C(=O)Nc2ccc3oc(nc3c2)c4ccccc4", "CCc1ccc(cc1)c2oc3ccc(NC(=O)Cc4ccc(Br)cc4)cc3n2", "OC[C@@H](O)[C@@H](O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@@H](O)C(=O)N2CCN(CC2)c3ccc(cc3)c4c5C=Cc(n5)c(c6ccccc6)c7ccc([nH]7)c(c8ccccc8)c9ccc(n9)c(c%10C=Cc4[nH]%10)c%11ccccc%11", "Clc1ccc(CSc2nc(Cl)c(\\\\C=C\\\\3/SC(=O)N(Cc4ccc(Cl)cc4Cl)C3=O)s2)cc1", "Clc1nc(sc1\\\\C=C\\\\2/SC(=O)N(Cc3ccc(Br)cc3)C2=O)N4CCCCC4", "[O-][N+](=O)c1ccc(CN2C(=O)S\\\\C(=C/c3sc(nc3Cl)N4CCCCC4)\\\\C2=O)cc1", "Clc1ccc(CN2C(=O)S\\\\C(=C/c3sc(nc3Cl)N4CCCCC4)\\\\C2=O)c(Cl)c1", "Clc1ccc(CSc2nc(Cl)c(\\\\C=C\\\\3/SC(=O)N(Cc4ccccc4)C3=O)s2)cc1", "Fc1ccc(CN2C(=O)S\\\\C(=C/c3sc(SCc4ccc(Cl)cc4)nc3Cl)\\\\C2=O)cc1", "Clc1ccc(CSc2nc(Cl)c(\\\\C=C\\\\3/SC(=O)N(Cc4ccc(Cl)cc4)C3=O)s2)cc1", "Clc1ccc(CSc2nc(Cl)c(\\\\C=C\\\\3/SC(=O)N(Cc4ccc(Br)cc4)C3=O)s2)cc1", "[O-][N+](=O)c1ccc(CN2C(=O)S\\\\C(=C/c3sc(SCc4ccc(Cl)cc4)nc3Cl)\\\\C2=O)cc1", "Cc1noc(n1)c2cccnc2", "Clc1nc(sc1\\\\C=C\\\\2/SC(=O)N(Cc3ccccc3)C2=O)N4CCCCC4", "Fc1ccc(CN2C(=O)S\\\\C(=C/c3sc(nc3Cl)N4CCCCC4)\\\\C2=O)cc1", "Clc1ccc(CN2C(=O)S\\\\C(=C/c3sc(nc3Cl)N4CCCCC4)\\\\C2=O)cc1", "COc1cc(cc(OC)c1OC)[C@H]2[C@@H]3[C@H](COC3=O)[C@@H](OC(=O)c4ccncc4)c5cc6OCOc6cc25", "COc1cc(cc(OC)c1OC)[C@H]2[C@@H]3[C@H](COC3=O)[C@@H](OC(=O)c4cccnc4)c5cc6OCOc6cc25", "COc1cc(cc(OC)c1OC)[C@H]2[C@@H]3[C@H](COC3=O)[C@@H](OC(=O)c4ccccn4)c5cc6OCOc6cc25", "COc1cc(cc(OC)c1OC)[C@H]2[C@@H]3[C@H](COC3=O)[C@@H](OC(=O)c4cccnc4Cl)c5cc6OCOc6cc25", "COc1cc(cc(OC)c1OC)[C@H]2[C@@H]3[C@H](COC3=O)[C@@H](OC(=O)c4ccc(Cl)nc4)c5cc6OCOc6cc25", "COc1cc(cc(OC)c1OC)[C@H]2[C@@H]3[C@H](COC3=O)[C@@H](OC(=O)c4cncc(Br)c4)c5cc6OCOc6cc25", "COc1cc(cc(OC)c1O)[C@H]2[C@@H]3[C@H](COC3=O)[C@H](OC(=O)c4ccncc4)c5cc6OCOc6cc25", "COc1cc(cc(OC)c1O)[C@H]2[C@@H]3[C@H](COC3=O)[C@H](OC(=O)c4cccnc4)c5cc6OCOc6cc25", "COc1cc(cc(OC)c1O)[C@H]2[C@@H]3[C@H](COC3=O)[C@H](OC(=O)c4ccccn4)c5cc6OCOc6cc25", "COc1cc(cc(OC)c1O)[C@H]2[C@@H]3[C@H](COC3=O)[C@H](OC(=O)c4cccnc4Cl)c5cc6OCOc6cc25", "COc1cc(cc(OC)c1O)[C@H]2[C@@H]3[C@H](COC3=O)[C@H](OC(=O)c4ccc(Cl)nc4)c5cc6OCOc6cc25", "COc1cc(cc(OC)c1O)[C@H]2[C@@H]3[C@H](COC3=O)[C@H](OC(=O)c4cncc(Br)c4)c5cc6OCOc6cc25", "C[C@H]1OCC\\\\C(=C\\\\C(=O)OCC23CCC(=C[C@H]2O[C@@H]4C[C@@H](OC(=O)\\\\C=C/C=C/[C@H]1O)[C@@]3(C)[C@]45CO5)C)\\\\C", "C[C@@H]1OCC\\\\C(=C\\\\C(=O)OCC23CCC(=C[C@H]2O[C@@H]4C[C@@H](OC(=O)\\\\C=C/C=C/[C@H]1O)[C@@]3(C)[C@]45CO5)C)\\\\C", "C\\\\C=C(\\\\C)/O[C@H]1[C@H](O)[C@]2(CO)[C@@H](C[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CCC(O[C@@H]6O[C@@H]([C@@H](O)[C@H](O[C@@H]7OC[C@H](O)[C@H](O)[C@H]7O[C@@H]8OC[C@@H](O)[C@H](O)[C@H]8O)[C@H]6O[C@]9(C)O[C@H](CO)[C@H](O)[C@H](O)[C@H]9O)C(=O)O)[C@@](C)(C=O)[C@@H]5CC[C@@]34C)[C@@H]2CC1(C)C)OC(=O)C", "CCN(CC)CCOC(=O)Cc1c(C)n(C(=O)c2ccc(Cl)cc2)c3ccc(OC)cc13", "COc1ccc2c(c1)c(CC(=O)OCCN3CCCC3)c(C)n2C(=O)c4ccc(Cl)cc4", "COc1ccc2c(c1)c(CC(=O)OCCN(C)C)c(C)n2C(=O)c3ccc(Cl)cc3", "Clc1ccc(CN2CCOCC2)cc1", "C(N1CCOCC1)c2ccc(cc2)c3ccccc3", "Cc1cccc(c1)C(=S)N2CCCCC2", "NCc1cnn(O)c1Br", "Brc1ccc(cc1)C(=S)N2CCOCC2", "S=C(NCCc1ccccc1)c2ccccc2", "Oc1ccc(cc1O)C(=S)NCCc2ccccc2", "COC(=O)c1ccc(cc1)C(=S)NCCc2ccccc2", "Brc1ccc(cc1)C(=S)NCCc2ccccc2", "Clc1ccc(cc1)C(=S)NCCc2ccccc2", "Clc1ccccc1C(=S)NCCc2ccccc2", "Clc1ccc(C(=S)NCCc2ccccc2)c(Cl)c1", "Clc1ccc(cc1Cl)C(=S)NCCc2ccccc2", "Fc1ccc(cc1)C(=S)NCCc2ccccc2", "S=C(Cc1ccc(cc1)c2ccccc2)NCCc3ccccc3", "Clc1ccccc1C(=O)N2CCOCC2", "Clc1ccc(cc1)C(=O)N2CCOCC2", "O=C(N1CCOCC1)c2ccc(cc2)c3ccccc3", "Clc1ccccc1CN2CCOCC2", "S=C(N1CCOCC1)c2ccccc2", "Clc1ccccc1C(=S)N2CCOCC2", "Clc1cccc(c1)C(=S)N2CCOCC2", "Clc1ccc(C(=S)N2CCOCC2)c(Cl)c1", "Clc1ccc(cc1Cl)C(=S)N2CCOCC2", "S=C(N1CCOCC1)c2ccc(cc2)C3CCCCC3", "[O-][N+](=O)c1cccc(c1)C(=S)N2CCOCC2", "Clc1ccccc1CC(=S)N2CCOCC2", "S=C(N1CCCCC1)c2ccccc2", "Cc1cc(C)c2nc(sc2c1)N3C(=O)c4ccccc4N=C3c5ccccc5", "C\\\\C(=N/NC(=O)c1ccc(F)cc1)\\\\c2ccncc2", "CCOc1cccc2sc(nc12)N3C(=O)c4ccccc4N=C3c5ccccc5", "Brc1cc(Br)c2N=C(N(C(=O)c2c1)c3nc4ccccc4s3)c5ccccc5", "COc1cccc2sc(nc12)N3C(=O)c4cc(Br)cc(Br)c4N=C3c5ccccc5", "COc1ccc2sc(nc2c1)N3C(=O)c4cc(Br)cc(Br)c4N=C3c5ccccc5", "COc1ccc2nc(sc2c1)N3C(=O)c4cc(Br)cc(Br)c4N=C3c5ccccc5", "Cc1cccc2sc(nc12)N3C(=O)c4cc(Br)cc(Br)c4N=C3c5ccccc5", "Cc1ccc2sc(nc2c1)N3C(=O)c4cc(Br)cc(Br)c4N=C3c5ccccc5", "Clc1ccc2sc(nc2c1)N3C(=O)c4cc(Br)cc(Br)c4N=C3c5ccccc5", "CCOc1ccc2nc(sc2c1)N3C(=O)c4cc(Br)cc(Br)c4N=C3c5ccccc5", "CCOc1ccc2sc(nc2c1)N3C(=O)c4cc(Br)cc(Br)c4N=C3c5ccccc5", "Brc1ccc2nc(sc2c1)N3C(=O)c4cc(Br)cc(Br)c4N=C3c5ccccc5", "[O-][N+](=O)c1ccc2nc(sc2c1)N3C(=O)c4cc(Br)cc(Br)c4N=C3c5ccccc5", "Cc1ccc(cc1)C2NC(C3CCCC2C3=NNC(=O)c4ccc(N)cc4)c5ccc(C)cc5", "Cc1cc(C)c2nc(sc2c1)N3C(=O)c4cc(Br)cc(Br)c4N=C3c5ccccc5", "CCOc1cccc2sc(nc12)N3C(=O)c4cc(Br)cc(Br)c4N=C3c5ccccc5", "COc1ccc2nc(sc2c1)N3C(=O)c4ccccc4N=C3c5ccccc5", "CC1=C(C(C(=C(C)N1)C(=O)Nc2ccc(Cl)cc2)c3ccccc3O)C(=O)Nc4ccc(Cl)cc4", "CCOc1ccc2nc(sc2c1)N3C(=O)c4ccccc4N=C3c5ccccc5", "CCOc1ccc2sc(nc2c1)N3C(=O)c4ccccc4N=C3c5ccccc5", "[O-][N+](=O)c1ccc2nc(sc2c1)N3C(=O)c4ccccc4N=C3c5ccccc5", "C\\\\C(=N/NC(=O)c1ccc(F)cc1)\\\\c2occc2", "CN(C)c1ccc(cc1)C2C(=C(C)NC(=C2C(=O)Nc3ccc(Cl)cc3)C)C(=O)Nc4ccc(Cl)cc4", "C\\\\C(=N/NC(=O)c1ccc(F)cc1)\\\\c2cccs2", "Fc1ccc(cc1)C(=O)N\\\\N=C\\\\c2c(Cl)cccc2Cl", "CC1=C(C(C(=C(C)N1)C(=O)Nc2ccc(Cl)cc2)c3cccc(Cl)c3)C(=O)Nc4ccc(Cl)cc4", "CC1=C(C(C\\\\C=C\\\\c2ccccc2)C(=C(C)N1)C(=O)Nc3ccc(Cl)cc3)C(=O)Nc4ccc(Cl)cc4", "CC1=C(CC(=C(C)N1)C(=O)Nc2ccc(Cl)cc2)C(=O)Nc3ccc(Cl)cc3", "OC1=C(Br)C(OC1=O)c2ccc(Cl)cc2", "Cc1nn(c(C)c1C2OC(=O)C(=C2Br)O)c3ccccc3", "OC1=C(Br)C(OC1=O)c2cn(nc2c3ccccc3)c4ccccc4", "CC(C)N(CCOc1cc(O)c2C(=O)C=C(Oc2c1)c3ccccc3)C(C)C", "Oc1cc(OCCNC2CCCCC2)cc3OC(=CC(=O)c13)c4ccccc4", "CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](Cc6ccccc6)C(C)(C)[C@@H]5CC[C@@]34C)[C@@H]2C1)C(=O)O", "Oc1cc(OCCN2CCOCC2)cc3OC(=CC(=O)c13)c4ccccc4", "Oc1cc(OCCN2CCN(CC2)c3ccccc3)cc4OC(=CC(=O)c14)c5ccccc5", "Oc1cc(OCCn2ccnc2)cc3OC(=CC(=O)c13)c4ccccc4", "Cc1nccn1CCOc2cc(O)c3C(=O)C=C(Oc3c2)c4ccccc4", "OCCNCCOc1cc(O)c2C(=O)C=C(Oc2c1)c3ccccc3", "OCCN(CCO)CCOc1cc(O)c2C(=O)C=C(Oc2c1)c3ccccc3", "Oc1cc(OCCNc2ccccc2)cc3OC(=CC(=O)c13)c4ccccc4", "CNCCOc1cc(O)c2C(=O)C=C(Oc2c1)c3ccccc3", "CC(C)NCCOc1cc(O)c2C(=O)C=C(Oc2c1)c3ccccc3", "Nc1ccc(cc1)C(=O)NN=C2C3CCCC2C(NC3c4ccccc4Cl)c5ccccc5Cl", "CCCCNCCOc1cc(O)c2C(=O)C=C(Oc2c1)c3ccccc3", "CC(C)(C)NCCOc1cc(O)c2C(=O)C=C(Oc2c1)c3ccccc3", "CN(C)CCOc1cc(O)c2C(=O)C=C(Oc2c1)c3ccccc3", "O=C1NC=NC2=C1C(c3ccccc3)c4ccc5cccnc5c4O2", "COc1ccc(cc1)C2NC(C3CCCC2C3=NNC(=O)c4ccc(N)cc4)c5ccc(OC)cc5", "NC1=NC(=O)NC2=C1C(c3ccccc3)c4ccccc4O2", "NC1=NC(=O)NC2=C1C(c3ccccc3)c4ccc(O)cc4O2", "NC1=NC(=O)NC2=C1C(c3ccccc3)c4ccc5ccccc5c4O2", "NC1=NC(=O)NC2=C1C(c3ccccc3)c4c(O2)ccc5ccccc45", "Oc1ccc2C(C3=C(Oc2c1)N=CNC3=O)c4ccccc4", "NC1=NC(=O)NC2=C1C(c3ccccc3)c4ccc5cccnc5c4O2", "NC1=NC(=S)NC2=C1C(c3ccccc3)c4ccccc4O2", "NC1=NC(=S)NC2=C1C(c3ccccc3)c4ccc(O)cc4O2", "NC1=NC(=S)NC2=C1C(c3ccccc3)c4ccc5ccccc5c4O2", "Nc1ccc(cc1)C(=O)NN=C2C3CCCC2C(NC3c4ccc(Br)cc4)c5ccc(Br)cc5", "NC1=NC(=S)NC2=C1C(c3ccccc3)c4c(O2)ccc5ccccc45", "NC1=NC(=S)NC2=C1C(c3ccccc3)c4ccc5cccnc5c4O2", "S=C1NC2=C(C(c3ccccc3)c4ccccc4O2)C(=S)N1", "Oc1ccc2C(C3=C(NC(=S)NC3=S)Oc2c1)c4ccccc4", "Nc1ccc(cc1)C(=O)NN=C2C3CCCC2C(NC3c4ccc(F)cc4)c5ccc(F)cc5", "S=C1NC2=C(C(c3ccccc3)c4ccc5ccccc5c4O2)C(=S)N1", "S=C1NC2=C(C(c3ccccc3)c4c(O2)ccc5ccccc45)C(=S)N1", "S=C1NC2=C(C(c3ccccc3)c4ccc5cccnc5c4O2)C(=S)N1", "O=C1NC=NC2=C1C(c3ccccc3)c4ccccc4O2", "Nc1ccc(cc1)C(=O)NN=C2C3CCCC2C(NC3c4ccc(Cl)cc4)c5ccc(Cl)cc5", "O=C1NC=NC2=C1C(c3ccccc3)c4ccc5ccccc5c4O2", "O=C1NC=NC2=C1C(c3ccccc3)c4c(O2)ccc5ccccc45", "CC(=O)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@H](CC[C@]4(C)[C@H]3CC[C@]12C)OC(=O)CCC(=O)O[C@@]5(CC[C@H]6[C@@H]7CCC8=Cc9oncc9C[C@]8(C)[C@H]7CC[C@]56C)C#C", "CC(=O)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@H](CC[C@]4(C)[C@H]3CC[C@]12C)OC(=O)CCC(=O)O[C@@]5(CC[C@H]6[C@@H]7CCC8=C(CNCCN)c9oncc9C[C@]8(C)[C@H]7CC[C@]56C)C#C", "Nc1ccc(cc1)C(=O)NN=C2C3CCCC2C(NC3c4ccccc4)c5ccccc5", "CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](Cc6ccccc6C(=O)O)C(C)(C)[C@@H]5CC[C@@]34C)[C@@H]2C1)C(=O)O", "CC1(C)CC[C@]2(CCCCC(=O)NC(Cc3ccccc3)C(=O)O)CC[C@]4(C)C(=CC[C@@H]5[C@@]6(C)CC[C@H](OCc7ccc(cc7)C(=O)O)C(C)(C)[C@@H]6CC[C@@]45C)[C@@H]2C1", "COc1ccc(CC(NC(=O)C[C@]23CCC(C)(C)C[C@H]2C4=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]5(C)[C@]4(C)CC3)C(=O)O)cc1OC", "COc1cccc2[nH]c3[C@@H](C[C@H]4C[C@@H]5N(CCc6c5[nH]c7ccccc67)C[C@@H]4C=C)NCCc3c12", "COc1cc(CC(NC(=O)C[C@]23CCC(C)(C)C[C@H]2C4=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]5(C)[C@]4(C)CC3)C(=O)O)cc(OC)c1", "Cc1ccccc1CC(NC(=O)C[C@]23CCC(C)(C)C[C@H]2C4=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]5(C)[C@]4(C)CC3)C(=O)O", "Cc1ccc(CC(NC(=O)C[C@]23CCC(C)(C)C[C@H]2C4=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]5(C)[C@]4(C)CC3)C(=O)O)cc1", "COc1ccccc1CC(NC(=O)C[C@]23CCC(C)(C)C[C@H]2C4=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]5(C)[C@]4(C)CC3)C(=O)O", "COc1cccc(CC(NC(=O)C[C@]23CCC(C)(C)C[C@H]2C4=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]5(C)[C@]4(C)CC3)C(=O)O)c1", "Cc1c(CN)cnn1O", "COc1ccc(CC(NC(=O)C[C@]23CCC(C)(C)C[C@H]2C4=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]5(C)[C@]4(C)CC3)C(=O)O)cc1", "CC1(C)CC[C@]2(CC(=O)NC(Cc3ccc4OCOc4c3)C(=O)O)CC[C@]5(C)C(=CC[C@@H]6[C@@]7(C)CC[C@H](O)C(C)(C)[C@@H]7CC[C@@]56C)[C@@H]2C1", "COc1cccc(CC(NC(=O)C[C@]23CCC(C)(C)C[C@H]2C4=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]5(C)[C@]4(C)CC3)C(=O)O)c1OC", "CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)C(C)(C)[C@@H]5CC[C@@]34C)[C@@H]2C1)C(=O)OCCCCCCCCCCCCC(=O)O", "CC[C@@H]1CN2CCc3c([nH]c4ccccc34)[C@H]2C[C@@H]1CCO", "C\\\\C=C/1\\\\CN2CC[C@@]34[C@@H]2C[C@@H]1[C@@H](CO)C3N(C(=O)C)c5ccccc45", "C\\\\C=C/1\\\\CN2CC[C@@]34[C@@H]2C[C@@H]1[C@@H](COC(=O)C)[C@@H]3N(C(=O)C)c5ccccc45", "C\\\\C=C(\\\\C)/C(=O)O[C@H]1[C@H]2[C@@H](C[C@@H](C)[C@@H]3C=CC(=O)[C@@]13C)OC(=O)[C@H]2C", "C[C@H]1N[C@@H](CCCCCCCCCCC(=O)C)CC[C@H]1O", "CC(C)(C)c1ccc(cc1)S(=O)(=O)N2CCN(CC2)C(c3ccccc3)c4ccccc4", "COc1cc2CCN(C)[C@H]3Cc4ccc(Oc5cc(C[C@@H]6NCCc7cc8Oc1c(Oc8cc67)c23)ccc5O)cc4", "CI.COc1cc2CCN(C)[C@H]3Cc4ccc(Oc5cc(C[C@@H]6N(C)CCc7cc8Oc1c(Oc8cc67)c23)ccc5O)cc4", "CN[C@H]1CC[C@]2(C)[C@H]3CC[C@]45C=N[C@@H](C)[C@H]4CC[C@H]5[C@@H]3CC=C2C1", "COc1cc2CCN(C)[C@H]3Cc4ccc(Oc5cc(C[C@H]6NCCc7cc(OC)c(Oc(c1O)c23)cc67)ccc5O)cc4", "CC[C@@H]1CN2CC[C@@]34[C@@H]2C[C@@H]1[C@H]5CO[C@H]([C@@H]6[C@H]7C[C@@H]8N(CC[C@]89[C@H]6N(C(=O)C)c%10ccccc9%10)C[C@H]7CC)N([C@H]35)c%11ccccc4%11", "CCCN(CC=C)C(=O)C(Cl)Cl", "O=C(CSc1nc2ccccc2s1)N\\\\N=C\\\\3/C(=O)Nc4ccccc34", "N(\\\\N=C\\\\c1ccc(Oc2ccccc2)cc1)c3nc4ccccc4s3", "[O-][N+](=O)c1ccc(Oc2ccc(\\\\C=N\\\\Nc3nc4ccccc4s3)cc2)cc1", "Cc1ccc(Oc2ccc(\\\\C=N\\\\Nc3nc4ccccc4s3)cc2)cc1", "CO[C@@H]1[C@H]2OC(C)(C)O[C@H]2O[C@@H]1[C@@H]3CC(=O)N(C(=O)N3c4occc4)c5ccc(cc5)C#N", "Clc1ccc(Oc2ccc(\\\\C=N\\\\Nc3nc4ccccc4s3)cc2)cc1", "Brc1ccc(Oc2ccc(\\\\C=N\\\\Nc3nc4ccccc4s3)cc2)cc1", "Fc1ccc(Oc2ccc(\\\\C=N\\\\Nc3nc4ccccc4s3)cc2)cc1", "Cc1cc(Oc2ccc(\\\\C=N\\\\Nc3nc4ccccc4s3)cc2)ccc1Cl", "Fc1ccc(CN(Cc2cccc(c2)C3=CC(=O)c4ccccc4O3)C(=O)\\\\C=C\\\\c5ccccc5)cc1", "N(\\\\N=C\\\\c1ccc(Oc2ccc3ccccc3c2)cc1)c4nc5ccccc5s4", "C1Oc2ccc(Oc3ccc(\\\\C=N\\\\Nc4nc5ccccc5s4)cc3)cc2O1", "O=C1Nc2ccccc2/C/1=N/Nc3nc4ccccc4s3", "O=C(CSc1nc2ccccc2s1)N\\\\N=C\\\\c3ccc(Oc4ccccc4)cc3", "[O-][N+](=O)c1ccc(Oc2ccc(\\\\C=N\\\\NC(=O)CSc3nc4ccccc4s3)cc2)cc1", "Cc1ccc(Oc2ccc(\\\\C=N\\\\NC(=O)CSc3nc4ccccc4s3)cc2)cc1", "Clc1ccc(Oc2ccc(\\\\C=N\\\\NC(=O)CSc3nc4ccccc4s3)cc2)cc1", "Brc1ccc(Oc2ccc(\\\\C=N\\\\NC(=O)CSc3nc4ccccc4s3)cc2)cc1", "Fc1ccc(Oc2ccc(\\\\C=N\\\\NC(=O)CSc3nc4ccccc4s3)cc2)cc1", "Cc1cc(Oc2ccc(\\\\C=N\\\\NC(=O)CSc3nc4ccccc4s3)cc2)ccc1Cl", "O=C(CSc1nc2ccccc2s1)N\\\\N=C\\\\c3ccc(Oc4ccc5ccccc5c4)cc3", "O=C(CSc1nc2ccccc2s1)N\\\\N=C\\\\c3ccc(Oc4ccc5OCOc5c4)cc3", "NNC(=O)COc1ccc2ccccc2c1", "O=C(COc1cccc2ccccc12)N\\\\N=C\\\\C=C\\\\c3ccccc3", "CO[C@@H]1[C@H]2OC(C)(C)O[C@H]2O[C@@H]1[C@@H]3CC(=O)N(C(=O)N3Cc4ccc(OC)cc4)c5cccc(c5)C(=O)C", "Oc1ccc2ccccc2c1\\\\C=N\\\\NC(=O)COc3cccc4ccccc34", "O=C(COc1cccc2ccccc12)N\\\\N=C\\\\c3c[nH]c4ccccc34", "C\\\\C(=N/NC(=O)COc1cccc2ccccc12)\\\\c3ccccc3", "C\\\\C(=N/NC(=O)COc1cccc2ccccc12)\\\\c3ccc(Cl)cc3", "CO[C@@H]1[C@H]2OC(C)(C)O[C@H]2O[C@@H]1[C@@H]3CC(=O)N(C(=O)N3Cc4ccccc4)c5cc(C)cc(C)c5", "C\\\\C(=N/NC(=O)COc1cccc2ccccc12)\\\\c3ccccc3Br", "C\\\\C(=N/NC(=O)COc1cccc2ccccc12)\\\\c3ccc(Br)cc3", "C\\\\C(=N/NC(=O)COc1cccc2ccccc12)\\\\c3ccc(cc3Br)C#N", "C\\\\C(=N/NC(=O)COc1cccc2ccccc12)\\\\c3cccc(c3)[N+](=O)[O-]", "CO[C@@H]1[C@H]2OC(C)(C)O[C@H]2O[C@@H]1[C@@H]3CC(=O)N(Cc4ccccc4)C(=O)N3Cc5ccccc5", "C\\\\C(=N/NC(=O)COc1cccc2ccccc12)\\\\c3ccc(cc3)[N+](=O)[O-]", "C\\\\C(=N/NC(=O)COc1cccc2ccccc12)\\\\c3ccccc3O", "C\\\\C(=N/NC(=O)COc1cccc2ccccc12)\\\\c3ccc(O)cc3", "C\\\\C(=N/NC(=O)COc1cccc2ccccc12)\\\\c3ccc(O)cc3O", "CO[C@@H]1[C@H]2OC(C)(C)O[C@H]2O[C@@H]1[C@@H]3CC(=O)N(C(=O)N3Cc4ccccc4)c5ccc(C)c(Cl)c5", "CC(C)OC(Cc1ccc(OCCc2noc(n2)C3CCCCC3)cc1)C(=O)O", "CC(C)Cc1ccc(cc1)c2onc(COc3ccc(CC(OC(C)C)C(=O)O)cc3)n2", "CC(C)C[C@@H]1N(C[C@H]2O[C@@H]3OC(C)(C)O[C@@H]3[C@H]2OCc4ccccc4)C(=O)N(C1=O)c5cccc(c5)C(=O)C", "CCCCc1ccc(cc1)N2C(=O)C[C@@H]([C@H]3O[C@@H]4OC(C)(C)O[C@@H]4[C@H]3OC)N(C2=O)c5occc5", "COc1ccc(CN2[C@@H](CC(=O)N(C2=O)c3ccc(cc3)C(=O)C)[C@H]4O[C@@H]5OC(C)(C)O[C@@H]5[C@H]4OCc6ccccc6)cc1", "CO[C@@H]1[C@H]2OC(C)(C)O[C@H]2O[C@@H]1[C@@H]3CC(=O)N(C(=O)N3c4occc4)c5cccc(c5)C(=O)C", "CO[C@@H]1[C@H]2OC(C)(C)O[C@H]2O[C@@H]1[C@@H]3CC(=O)N(C(=O)N3Cc4ccccc4O)c5cccc(c5)C(=O)C", "CO[C@@H]1[C@H]2OC(C)(C)O[C@H]2O[C@@H]1[C@@H]3CC(=O)N(Cc4ccccc4)C(=O)N3Cc5ccccc5O", "CO[C@@H]1[C@H]2OC(C)(C)O[C@H]2O[C@@H]1[C@@H]3CC(=O)N(C(=O)N3Cc4ccccc4O)c5ccc(C)c(Cl)c5", "CO[C@@H]1[C@H]2OC(C)(C)O[C@H]2O[C@@H]1[C@@H]3CC(=O)N(C(=O)N3Cc4ccccc4O)c5ccc(cc5)C#N", "CC1(C)O[C@H]2O[C@H]([C@@H]3CC(=O)N(C(=O)N3)c4ccc(Cl)cc4)[C@H](OCc5ccccc5)[C@H]2O1", "CO[C@H]1[C@H](O[C@@H]2OC(C)(C)O[C@H]12)[C@@H]3CC(=O)N(C(=O)N3)c4ccc(Cl)cc4", "CC(=O)c1ccc(cc1)N2C(=O)C[C@H](NC2=O)[C@H]3O[C@@H]4OC(C)(C)O[C@@H]4[C@H]3OCc5ccccc5", "CCCCCCCCCCCCCCCCNC(CC(=O)NO)[C@H]1O[C@@H]2OC(C)(C)O[C@@H]2[C@H]1OCc3ccccc3", "CC1(C)O[C@H]2O[C@H](C(CC(=O)NO)N3CCCCC3)[C@H](OCc4ccccc4)[C@H]2O1", "CCOC(=O)CC([C@H]1O[C@@H]2OC(C)(C)O[C@@H]2[C@H]1OC)N3CCC[C@H]3C(=O)OCC4c5ccccc5c6ccccc46", "CCOC(=O)CC([C@H]1O[C@@H]2OC(C)(C)O[C@@H]2[C@H]1OC)N3CCC[C@H]3C(=O)O", "Clc1ccc(CN(Cc2cccc(c2)C3=CC(=O)c4ccccc4O3)C(=O)\\\\C=C\\\\c5ccccc5)cc1", "CCOC(=O)CC(N[C@@H](C(C)C)C(=O)O)[C@H]1O[C@@H]2OC(C)(C)O[C@@H]2[C@H]1OC", "CCOC(=O)CC([C@H]1O[C@@H]2OC(C)(C)O[C@@H]2[C@H]1OC)N(Cc3ccccc3O)C(=O)Nc4cccc(c4)C(=O)C", "CCOC(=O)CC([C@H]1O[C@@H]2OC(C)(C)O[C@@H]2[C@H]1OC)N(Cc3ccccc3O)C(=O)NCc4ccccc4", "CCOC(=O)CC(NC(=O)Nc1ccccc1)[C@H]2O[C@@H]3OC(C)(C)O[C@@H]3[C@H]2OCc4ccccc4", "CCOC(=O)CC([C@H]1O[C@@H]2OC(C)(C)O[C@@H]2[C@H]1OCc3ccccc3)N(Cc4ccccc4)C(=O)NCc5ccccc5", "CN1c2c(C)n(nc2c3ccccc3S1(=O)=O)c4ccc(cc4)C(=O)\\\\C=C\\\\c5cccc(Cl)c5", "CCOC(=O)CC([C@H]1O[C@@H]2OC(C)(C)O[C@@H]2[C@H]1OC)N(CCC(C)CC(C)C)C(=O)NCc3ccccc3", "CCOC(=O)CC([C@H]1O[C@@H]2OC(C)(C)O[C@@H]2[C@H]1OCc3ccccc3)N(C(=O)Nc4ccc(Cl)cc4)c5occc5", "CCCCN(C(CC(=O)OCC)[C@H]1O[C@@H]2OC(C)(C)O[C@@H]2[C@H]1OC)C(=S)Nc3cccc(Cl)c3", "CCCCN(C(CC(=O)OCC)[C@H]1O[C@@H]2OC(C)(C)O[C@@H]2[C@H]1OC)C(=S)Nc3ccccc3Cl", "CO[C@H]1[C@H](O[C@@H]2OC(C)(C)O[C@H]12)[C@H](CC(=O)N)N(C3O[C@H]4OC(C)(C)O[C@H]4[C@@H]3OC)C(=O)Nc5ccc(C)c(Cl)c5", "CO[C@H]1[C@H](O[C@@H]2OC(C)(C)O[C@H]12)[C@H](CC(=O)N)N(C3O[C@H]4OC(C)(C)O[C@H]4[C@@H]3OC)C(=O)Nc5cccc(c5)C(=O)C", "CO[C@H]1[C@H](O[C@@H]2OC(C)(C)O[C@H]12)[C@H](CC(=O)N)N(C3O[C@H]4OC(C)(C)O[C@H]4[C@@H]3OC)C(=O)NCc5ccccc5", "Fc1ccc(cc1)c2cc3ccccc3[nH]2", "Clc1ccc(cc1Cl)c2cc3ccccc3[nH]2", "[O-][N+](=O)c1cc(c2[nH]c(cc2c1)c3ccccc3)[N+](=O)[O-]", "[O-][N+](=O)c1cc(ccc1Cl)c2cc3ccccc3[nH]2", "Fc1ccc(cc1)c2[nH]c3ccccc3c2C=O", "Clc1ccc(cc1)c2[nH]c3ccccc3c2C=O", "Nc1ccc(cc1)c2[nH]c3ccccc3c2C=O", "[O-][N+](=O)c1cc(ccc1Cl)c2[nH]c3ccccc3c2C=O", "Cc1ccc(cc1)c2cc(N)cc(c2)c3ccccc3", "Cc1ccc(cc1)c2cc(N)cc(c2)c3ccc(C)cc3", "Nc1cc(cc(c1)c2ccc(Cl)cc2)c3ccccc3", "Nc1cc(cc(c1)c2ccc(Cl)cc2)c3ccc(Cl)cc3", "CCOc1ccc(cc1)c2cc(N)cc(c2)c3ccccc3", "CCOc1ccc(cc1)c2cc(N)cc(c2)c3ccc(OCC)cc3", "COc1ccc(cc1)c2cc(N)cc(c2)c3ccccc3", "CSc1ccc(cc1)c2cc(N)cc(c2)c3ccccc3", "COC(=O)c1cc(Cl)nn1C", "CS(=O)(=O)c1ccc(cc1)c2cc(N)cc(c2)c3ccccc3", "Cc1c(CC(=O)O)cc(cc1c2ccccc2)c3ccccc3", "COc1ccccc1C(=O)N2C=C(C)N(C2=S)c3ccccc3", "COc1ccccc1C(=O)N2C=C(C)N(C2=S)c3ccccc3Cl", "COc1ccccc1C(=O)N2C=C(C)N(C2=S)c3cccc(Cl)c3", "COc1ccccc1C(=O)N2C=C(C)N(C2=S)c3ccc(Cl)cc3", "COc1ccccc1C(=O)N2C=C(C)N(C2=S)c3ccc(Br)cc3", "COc1ccccc1C(=O)N2C=C(C)N(C2=S)c3ccccc3C", "COc1ccccc1C(=O)N2C=C(C)N(C2=S)c3cccc(C)c3", "COc1ccccc1N2C(=CN(C(=O)c3ccccc3OC)C2=S)C", "COc1cccc(c1)N2C(=CN(C(=O)c3ccccc3OC)C2=S)C", "CC1=CN(C(=O)c2cccc(C)c2)C(=S)N1c3ccccc3Cl", "CC1=CN(C(=O)c2cccc(C)c2)C(=S)N1c3cccc(Cl)c3", "CC1=CN(C(=O)c2cccc(C)c2)C(=S)N1c3ccc(Cl)cc3", "CC1=CN(C(=O)c2cccc(C)c2)C(=S)N1c3ccc(Br)cc3", "CC1=CN(C(=O)c2cccc(C)c2)C(=S)N1c3ccccc3C", "CC1=CN(C(=O)c2cccc(C)c2)C(=S)N1c3cccc(C)c3", "COc1ccccc1N2C(=CN(C(=O)c3cccc(C)c3)C2=S)C", "COc1cccc(c1)N2C(=CN(C(=O)c3cccc(C)c3)C2=S)C", "CC1=CN(C(=O)c2ccc(C)cc2)C(=S)N1c3ccccc3", "CC1=CN(C(=O)c2ccc(C)cc2)C(=S)N1c3ccccc3Cl", "CC1=CN(C(=O)c2ccc(C)cc2)C(=S)N1c3cccc(Cl)c3", "CC1=CN(C(=O)c2ccc(C)cc2)C(=S)N1c3ccc(Cl)cc3", "CC1=CN(C(=O)c2ccc(C)cc2)C(=S)N1c3ccc(Br)cc3", "CC1=CN(C(=O)c2ccc(C)cc2)C(=S)N1c3ccccc3C", "CC1=CN(C(=O)c2ccc(C)cc2)C(=S)N1c3cccc(C)c3", "COc1ccccc1N2C(=CN(C(=O)c3ccc(C)cc3)C2=S)C", "COc1cccc(c1)N2C(=CN(C(=O)c3ccc(C)cc3)C2=S)C", "Clc1ccc(CNCc2cccc(c2)C3=CC(=O)c4ccccc4O3)cc1Cl", "O=C(\\\\C=C\\\\c1ccccc1)N(Cc2ccccc2)Cc3cccc(c3)C4=CC(=O)c5ccccc5O4", "Clc1ccc(CN(Cc2cccc(c2)C3=CC(=O)c4ccccc4O3)C(=O)\\\\C=C\\\\c5ccccc5)c(Cl)c1", "Clc1ccc(CN(Cc2cccc(c2)C3=CC(=O)c4ccccc4O3)C(=O)\\\\C=C\\\\c5ccccc5)cc1Cl", "O=C1C=C(Oc2ccccc12)c3cccc(CNCc4ccccc4)c3", "Fc1ccc(CNCc2cccc(c2)C3=CC(=O)c4ccccc4O3)cc1", "Clc1ccc(CNCc2cccc(c2)C3=CC(=O)c4ccccc4O3)cc1", "Clc1ccc(CNCc2cccc(c2)C3=CC(=O)c4ccccc4O3)c(Cl)c1", "CS(=O)(=O)N1CCN(CC1)C(c2ccccc2)c3ccccc3", "Cc1ccc(cc1)S(=O)(=O)N2CCN(CC2)C(c3ccccc3)c4ccccc4", "Clc1ccc(cc1)S(=O)(=O)N2CCN(CC2)C(c3ccccc3)c4ccccc4", "Sc1oc(nn1)c2cccc(Cl)c2", "Cc1onc(C)c1S(=O)(=O)N2CCN(CC2)C(c3ccccc3)c4ccccc4", "CCC1C(N(C(CC1=O)c2ccc(Cl)cc2)C(=O)CCl)c3ccc(Cl)cc3", "CC1C(N(C(C(C)C1=O)c2ccc(Cl)cc2)C(=O)CCl)c3ccc(Cl)cc3", "COc1ccc(cc1)C2CC(=O)C(C)C(N2C(=O)CCl)c3ccc(OC)cc3", "NC(Cc1occc1)C(=O)O", "COc1ccc(cc1)C2C(C)C(=O)C(C)C(N2C(=O)CCl)c3ccc(OC)cc3", "CC1C(N(C(CC1=O)c2ccc(C)cc2)C(=O)CCl)c3ccc(C)cc3", "CC1C(N(C(C(C)C1=O)c2ccc(C)cc2)C(=O)CCl)c3ccc(C)cc3", "ClCC(=O)N1C(CC(=O)CC1c2ccccc2)c3ccccc3", "[O-][N+](=O)c1ccccc1\\\\C=N\\\\N=C\\\\2/C[C@H](S[C@H](C2)c3ccccc3)c4ccccc4", "CC1C(N(C(CC1=O)c2ccccc2)C(=O)CCl)c3ccccc3", "CC(C)C1C(N(C(CC1=O)c2ccccc2)C(=O)CCl)c3ccccc3", "Fc1ccc(cc1)C(=O)N\\\\N=C\\\\c2cccc(Oc3ccccc3)c2", "Oc1c(Cl)cc(Cl)cc1\\\\C=N\\\\NC(=O)c2ccc(F)cc2", "Fc1ccc(cc1)C(=O)N\\\\N=C\\\\c2cc(Cl)cc(Cl)c2", "Fc1ccc(cc1)C(=O)N\\\\N=C\\\\c2ccc(Cl)cc2Cl", "CN(C)CCCOc1ccc(\\\\C=N\\\\NC(=O)c2ccc(F)cc2)cc1", "Fc1ccc(cc1)C(=O)N\\\\N=C\\\\c2cccc(F)c2", "C\\\\C(=N/NC(=O)c1ccc(F)cc1)\\\\c2ccc[nH]2", "C\\\\C(=N/NC(=O)c1ccc(F)cc1)\\\\c2ccc(Cl)s2", "COc1cc(ccc1O)C2CC(=NN2)C(=NNc3ccc(cc3)[N+](=O)[O-])C4=NNC(C4)c5ccc(O)c(OC)c5", "C\\\\C(=N/NC(=O)c1ccc(F)cc1)\\\\c2ccccc2", "C\\\\C(=N/NC(=O)c1ccc(F)cc1)\\\\c2ccc(Br)s2", "C\\\\C(=N/NC(=O)c1ccc(F)cc1)\\\\C=C", "CC\\\\C(=N/NC(=O)c1ccc(F)cc1)\\\\c2ccccc2", "Brc1ccccc1NC(=S)Nc2nc3c(ccc4ccccc34)s2", "COc1cc(\\\\C=C\\\\C(=O)\\\\C(=N/Nc2ccc(cc2)[N+](=O)[O-])\\\\C3=NC(=S)NC(C3)c4ccc(O)c(OC)c4)ccc1O", "Brc1ccc(NC(=S)Nc2nc3c(ccc4ccccc34)s2)cc1", "Clc1ccccc1NC(=S)Nc2nc3c(ccc4ccccc34)s2", "Clc1ccc(NC(=S)Nc2nc3c(ccc4ccccc34)s2)cc1", "[O-][N+](=O)c1ccccc1NC(=S)Nc2nc3c(ccc4ccccc34)s2", "COc1cc(ccc1O)C2ONC(=C2)C(=NNc3ccc(cc3)[N+](=O)[O-])C4=CC(ON4)c5ccc(O)c(OC)c5", "[O-][N+](=O)c1ccc(NC(=S)Nc2nc3c(ccc4ccccc34)s2)cc1", "Cc1ccccc1NC(=S)Nc2nc3c(ccc4ccccc34)s2", "Cc1ccc(NC(=S)Nc2nc3c(ccc4ccccc34)s2)cc1", "O=C(Nc1ccccc1)Nc2nc3c(ccc4ccccc34)s2", "S=C(Nc1ccccc1)Nc2nc3c(ccc4ccccc34)s2", "COc1ccc(NC(=O)Nc2nc3c(ccc4ccccc34)s2)cc1", "COc1ccccc1NC(=O)Nc2nc3c(ccc4ccccc34)s2", "COc1ccc(NC(=O)Nc2nc3c(ccc4ccccc34)s2)c(OC)c1", "COc1cc(NC(=O)Nc2nc3c(ccc4ccccc34)s2)cc(OC)c1OC", "Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5", "O=C(Nc1ccc(Oc2ccccc2)cc1)Nc3nc4c(ccc5ccccc45)s3", "Fc1ccccc1NC(=O)Nc2nc3c(ccc4ccccc34)s2", "Fc1ccc(NC(=O)Nc2nc3c(ccc4ccccc34)s2)cc1", "Brc1ccccc1NC(=O)Nc2nc3c(ccc4ccccc34)s2", "COc1cc(\\\\C=C\\\\C(=O)C(=NNc2ccc(Cl)cc2)C(=O)\\\\C=C\\\\c3ccc(O)c(OC)c3)ccc1O", "Brc1ccc(NC(=O)Nc2nc3c(ccc4ccccc34)s2)cc1", "Clc1ccccc1NC(=O)Nc2nc3c(ccc4ccccc34)s2", "Clc1ccc(NC(=O)Nc2nc3c(ccc4ccccc34)s2)cc1", "FC(F)(F)c1ccc(Cl)c(NC(=O)Nc2nc3c(ccc4ccccc34)s2)c1", "COc1cc(\\\\C=C\\\\C(=O)C(=NNc2ccc(cc2)[N+](=O)[O-])C(=O)\\\\C=C\\\\c3ccc(O)c(OC)c3)ccc1O", "[O-][N+](=O)c1ccccc1NC(=O)Nc2nc3c(ccc4ccccc34)s2", "[O-][N+](=O)c1ccc(NC(=O)Nc2nc3c(ccc4ccccc34)s2)cc1", "Cc1ccccc1NC(=O)Nc2nc3c(ccc4ccccc34)s2", "Cc1ccc(NC(=O)Nc2nc3c(ccc4ccccc34)s2)cc1", "COc1cc(\\\\C=C\\\\C(=O)\\\\C(=N/Nc2ccc(cc2)[N+](=O)[O-])\\\\C(=O)C(Br)C(Br)c3ccc(O)c(OC)c3)ccc1O", "COc1ccc(NC(=S)Nc2nc3c(ccc4ccccc34)s2)cc1", "COc1ccccc1NC(=S)Nc2nc3c(ccc4ccccc34)s2", "Fc1ccccc1NC(=S)Nc2nc3c(ccc4ccccc34)s2", "Fc1ccc(NC(=S)Nc2nc3c(ccc4ccccc34)s2)cc1", "CSc1ccc(\\\\C=C\\\\C(=O)c2ccc(F)cc2)cc1", "CSc1ccc(\\\\C=C\\\\C(=O)c2ccc(cc2)c3ccccc3)cc1", "CSc1ccc(\\\\C=C\\\\C(=O)c2cccs2)cc1", "CSc1ccc(\\\\C=C\\\\C(=O)c2ccc(Br)s2)cc1", "OC(=O)c1cc(nc2cc(ccc12)N3CCNCC3)c4ccccc4", "NS(=O)(=O)c1ccc(Nc2ccc3c(cc(nc3c2)c4ccccc4)C(=O)O)cc1", "NC(=O)Nc1ccc2c(cc(nc2c1)c3ccccc3)C(=O)O", "NC(=S)Nc1ccc2c(cc(nc2c1)c3ccccc3)C(=O)O", "OC(=O)c1cc(nc2cc(Nc3nccs3)ccc12)c4ccccc4", "COc1ccc(Nc2ccc3c(cc(nc3c2)c4ccccc4)C(=O)O)cc1", "COc1cc(ccc1O)C2CC(=NC(=S)N2)C(=NNc3ccc(cc3)[N+](=O)[O-])C4=NC(=S)NC(C4)c5ccc(O)c(OC)c5", "OC(=O)CNc1ccc2c(cc(nc2c1)c3ccccc3)C(=O)O", "OC(=O)c1cc(nc2cc(Nn3cnnc3)ccc12)c4ccccc4", "OC(=O)c1cc(nc2cc(Nc3ccncc3)ccc12)c4ccccc4", "OC(=O)c1cc(nc2cc(Cl)ccc12)c3ccccc3", "Cc1ccc(Nc2ccc3c(cc(nc3c2)c4ccccc4)C(=O)O)cc1", "OC(=O)c1cc(nc2cc(Nc3ccc(cc3)[N+](=O)[O-])ccc12)c4ccccc4", "OC(=O)c1cc(nc2cc(Nc3ccccc3)ccc12)c4ccccc4", "COc1ccc(NN=C(C(=O)\\\\C=C\\\\c2ccc(O)c(OC)c2)C(=O)\\\\C=C\\\\c3ccc(O)c(OC)c3)cc1", "COc1cc(\\\\C=C\\\\C(=O)C(=NNc2ccc(C)cc2)C(=O)\\\\C=C\\\\c3ccc(O)c(OC)c3)ccc1O", "COc1cc(\\\\C=C\\\\C(=O)C(=NNc2ccccc2)C(=O)\\\\C=C\\\\c3ccc(O)c(OC)c3)ccc1O", "CCN1C=C(C(=O)c2cc(F)c(cc12)N3CCNCC3)c4oc(nn4)c5ccccc5", "CCOC(=O)C(C(CC(=O)C(=NNc1ccc(Cl)cc1)C(=O)CC(C(C(=O)C)C(=O)OCC)c2ccc(O)c(OC)c2)c3ccc(O)c(OC)c3)C(=O)C", "COc1cc(\\\\C=C\\\\c2cc(\\\\C=C\\\\c3ccc(O)c(OC)c3)n(NN)n2)ccc1O", "COc1cc(\\\\C=C\\\\c2cc(\\\\C=C\\\\c3ccc(O)c(OC)c3)n(N=O)n2)ccc1O", "COc1cc(\\\\C=C\\\\c2cc(\\\\C=C\\\\c3ccc(O)c(OC)c3)n(n2)C(=O)CCl)ccc1O", "CC[C@H]1OC(=O)[C@H](C)[C@@H](O)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)N(C)C)[C@](C)(O)C[C@@H](C)CN(CCCNC(=O)[C@H]3[C@H](C)C[C@H]4[C@@H]5CCC6=CC(=O)C=C[C@]6(C)[C@@]5(F)[C@@H](O)C[C@]34C)[C@H](C)[C@@H](O)[C@]1(C)O", "COc1cc(\\\\C=C\\\\c2cc(\\\\C=C\\\\c3cc(OC)c(O)c(c3)N=Nc4ccc(Cl)cc4)[nH]n2)cc(N=Nc5ccc(Cl)cc5)c1O", "COc1cc(\\\\C=C\\\\c2cc(\\\\C=C\\\\c3cc(OC)c(O)c(c3)N=Nc4ccc(Cl)cc4)n(C=O)n2)cc(N=Nc5ccc(Cl)cc5)c1O", "COc1cc(\\\\C=C\\\\c2cc(\\\\C=C\\\\c3cc(OC)c(O)c(c3)N=Nc4ccc(Cl)cc4)n(CO)n2)cc(N=Nc5ccc(Cl)cc5)c1O", "CN(C)c1ccc(\\\\C=C\\\\C(=O)c2ccc(O)cc2O)cc1", "Oc1ccc(C(=O)\\\\C=C\\\\c2ccccc2[N+](=O)[O-])c(O)c1", "CCN(CC)C(=O)c1ccccc1N", "Oc1ccc(C(=O)\\\\C=C\\\\c2cccc(c2)[N+](=O)[O-])c(O)c1", "Cc1ccsc1c2nc3ccc4C(=O)c5ccccc5C(=O)c4c3[nH]2", "Cc1ccc(s1)c2nc3ccc4C(=O)c5ccccc5C(=O)c4c3[nH]2", "Cn1cccc1c2nc3ccc4C(=O)c5ccccc5C(=O)c4c3[nH]2", "Cn1cccc1c2nc3cc4nc5ccccc5nc4cc3[nH]2", "O=C1c2ccccc2C(=O)c3c1ccc4nc([nH]c34)c5ccc6ccccc6c5", "Oc1ccc2ccccc2c1c3nc4ccc5C(=O)c6ccccc6C(=O)c5c4[nH]3", "Cc1ccsc1c2nc3cc4nc5ccccc5nc4cc3[nH]2", "Cc1ccc(s1)c2nc3cc4nc5ccccc5nc4cc3[nH]2", "CCCCNC(=O)c1ccccc1N", "c1ccc2c(cccc2c1)c3nc4cc5nc6ccccc6nc5cc4[nH]3", "c1ccc2cc(ccc2c1)c3nc4cc5nc6ccccc6nc5cc4[nH]3", "Oc1ccc2ccccc2c1c3nc4cc5nc6ccccc6nc5cc4[nH]3", "Cc1ccsc1\\\\C=N/C(=N)Nc2nc3ccccc3[nH]2", "CO[C@H]1[C@H](O[C@@H]2OC(C)(C)O[C@H]12)[C@H](CC(=O)N)N(Cc3ccc(OC)cc3)C(=O)NCc4ccccc4", "Cc1ccc(\\\\C=N/C(=N)Nc2nc3ccccc3[nH]2)s1", "Cn1cccc1\\\\C=N\\\\C(=N)Nc2nc3ccccc3[nH]2", "N=C(Nc1nc2ccccc2[nH]1)\\\\N=C\\\\c3cccc4ccccc34", "N=C(Nc1nc2ccccc2[nH]1)\\\\N=C/c3ccc4ccccc4c3", "Oc1ccc2ccccc2c1\\\\C=N\\\\C(=N)Nc3nc4ccccc4[nH]3", "[Br-].[Br-].Clc1ccccc1N2CCN(CC2)C(=O)c3ccc[n+](Cc4ccc(C[n+]5cccc(c5)C(=O)N6CCN(CC6)c7ccccc7Cl)cc4)c3", "Cl.Clc1ccccc1N2CCN(CC2)C(=O)c3cccnc3", "Cl.O=C(N1CCN(CC1)c2ncccn2)c3cccnc3", "[Cl-].Clc1ccccc1N2CCN(CC2)C(=O)c3ccc[n+](Cc4ccccc4)c3", "[Br-].[Br-].O=C(N1CCN(CC1)c2ncccn2)c3ccc[n+](Cc4ccc(C[n+]5cccc(c5)C(=O)N6CCN(CC6)c7ncccn7)cc4)c3", "[Cl-].Fc1ccccc1C[n+]2cccc(c2)C(=O)N3CCN(CC3)c4ccccc4Cl", "[Cl-].Fc1ccc(C[n+]2cccc(c2)C(=O)N3CCN(CC3)c4ccccc4Cl)cc1", "[Br-].Clc1ccccc1N2CCN(CC2)C(=O)c3ccc[n+](Cc4ccc(Br)cc4)c3", "[Cl-].COc1cc(cc(OC)c1OC)C(=O)[n+]2cccc(c2)C(=O)N3CCN(CC3)c4ccccc4Cl", "CN(C)C(=O)c1ccccc1N", "[Br-].Brc1ccc(C[n+]2cccc(c2)C(=O)N3CCN(CC3)c4ncccn4)cc1", "Cl.Clc1ccccc1N2CCN(CC2)C(=O)C3CCCN(Cc4ccccc4)C3", "Cl.Fc1ccccc1CN2CCCC(C2)C(=O)N3CCN(CC3)c4ccccc4Cl", "Cl.Fc1ccc(CN2CCCC(C2)C(=O)N3CCN(CC3)c4ccccc4Cl)cc1", "Br.Clc1ccccc1N2CCN(CC2)C(=O)C3CCCN(Cc4ccc(Br)cc4)C3", "Br.[O-][N+](=O)c1ccc(CN2CCCC(C2)C(=O)N3CCN(CC3)c4ccccc4Cl)cc1", "Cl.Clc1ccccc1N2CCN(CC2)C(=O)C3CCCN(C3)C(=O)c4ccc(Br)cc4", "Br.Brc1ccc(CN2CCCC(C2)C(=O)N3CCN(CC3)c4ccccn4)cc1", "Cl.COc1cc(cc(OC)c1OC)C(=O)N2CCCC(C2)C(=O)N3CCN(CC3)c4ncccn4", "Br.Clc1ccccc1N2CCN(CC2)C(=O)C3CCCN(Cc4ccc(CN5CCCC(C5)C(=O)N6CCN(CC6)c7ccccc7Cl)cc4)C3", "Cl.Cl.O=C(C1CCCN(Cc2ccc(CN3CCCC(C3)C(=O)N4CCN(CC4)c5ncccn5)cc2)C1)N6CCN(CC6)c7ncccn7", "[Cl-].[Cl-].O=C(N1CCN(CC1)C(=O)c2ccc[n+](Cc3ccccc3)c2)c4ccc[n+](Cc5ccccc5)c4", "Cl.Cl.O=C(C1CCCN(Cc2ccccc2)C1)N3CCN(CC3)C(=O)C4CCCN(Cc5ccccc5)C4", "Clc1ccccc1C2=CC(=O)c3ccccc3O2", "Clc1cc(I)c2OC(=CC(=O)c2c1)c3ccccc3", "Cc1c(Cl)cc2C(=O)C=C(Oc2c1I)c3ccccc3", "Cc1c(Cl)cc2C(=O)C=C(Oc2c1I)c3ccc(Cl)cc3Cl", "Clc1ccc(C2=CC(=O)c3cc(Cl)cc(I)c3O2)c(Cl)c1", "CC1(C)CC[C@]2(CC(=O)NC(Cc3ccccc3)C(=O)O)CC[C@]4(C)C(=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]45C)[C@@H]2C1", "CC1(C)CC[C@]2(CC(=O)NC(Cc3cnc[nH]3)C(=O)O)CC[C@]4(C)C(=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]45C)[C@@H]2C1", "Fc1ccc(Cl)cc1I", "CC1(C)CC[C@]2(CC(=O)NC(CC3=CCc4ccccc34)C(=O)O)CC[C@]5(C)C(=CC[C@@H]6[C@@]7(C)CC[C@H](O)C(C)(C)[C@@H]7CC[C@@]56C)[C@@H]2C1", "CSCCC(NC(=O)C[C@]12CCC(C)(C)C[C@H]1C3=CC[C@@H]4[C@@]5(C)CC[C@H](O)C(C)(C)[C@@H]5CC[C@@]4(C)[C@]3(C)CC2)C(=O)O", "CC1(C)CC[C@]2(CC(=O)NC(Cc3ccccc3Cl)C(=O)O)CC[C@]4(C)C(=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]45C)[C@@H]2C1", "CC1(C)CC[C@]2(CC(=O)NC(Cc3cccc(Cl)c3)C(=O)O)CC[C@]4(C)C(=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]45C)[C@@H]2C1", "CC1(C)CC[C@]2(CC(=O)NC(Cc3ccc(Cl)cc3)C(=O)O)CC[C@]4(C)C(=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]45C)[C@@H]2C1", "Fc1cc(Cl)ccc1Cl", "CC1(C)CC[C@]2(CC(=O)NC(Cc3ccc(F)cc3)C(=O)O)CC[C@]4(C)C(=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]45C)[C@@H]2C1", "CC1(C)CC[C@]2(CC(=O)NC(Cc3ccc(cc3)[N+](=O)[O-])C(=O)O)CC[C@]4(C)C(=CC[C@@H]5[C@@]6(C)CC[C@H](O)C(C)(C)[C@@H]6CC[C@@]45C)[C@@H]2C1", "CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)C(C)(C)[C@@H]5CC[C@@]34C)[C@@H]2C1)C(=O)OCCC(=O)O", "CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)C(C)(C)[C@@H]5CC[C@@]34C)[C@@H]2C1)C(=O)OCCCCC(=O)O", "CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)C(C)(C)[C@@H]5CC[C@@]34C)[C@@H]2C1)C(=O)OCCCCCCC(=O)O", "COc1ccccc1NC(=O)c2ccccc2N", "CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)C(C)(C)[C@@H]5CC[C@@]34C)[C@@H]2C1)C(=O)OCCCCCCCCC(=O)O", "CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)C(C)(C)[C@@H]5CC[C@@]34C)[C@@H]2C1)C(=O)OCCCCCCCCCCC(=O)O", "CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)C(C)(C)[C@@H]5CC[C@@]34C)[C@@H]2C1)C(=O)NC(Cc6ccccc6)C(=O)O", "CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)C(C)(C)[C@@H]5CC[C@@]34C)[C@@H]2C1)C(=O)NC(CCC(=O)O)C(=O)O", "CC1=C(C(NC(=S)N1)c2cccc(c2)[N+](=O)[O-])C(=O)N(CCCl)CCCl", "Fc1cc(Cl)ccc1Br", "CN(C)c1ccc(cc1)C2NC(=O)NC(=C2C(=O)Nc3ccc(NC(=O)C)cc3)C", "CN(C)c1ccc(cc1)C2NC(=O)NC(=C2C(=O)NNC(=O)c3ccccc3O)C", "CN(C)c1ccc(cc1)C2NC(=S)NC(=C2C(=O)NNC(=O)c3ccccc3O)C", "CC1=C(C(NC(=O)N1)c2cccc(c2)[N+](=O)[O-])C(=O)NNC(=O)c3ccccc3O", "CC1=C(C(NC(=S)N1)c2cccc(c2)[N+](=O)[O-])C(=O)NNC(=O)c3ccccc3O", "COc1cc(cc(OC)c1OC)C2NC(=S)NC(=C2C(=O)NNC(=O)c3ccccc3O)C", "C\\\\C=C\\\\C(=O)NNC(=O)C1=C(C)NC(=S)NC1c2cccc(c2)[N+](=O)[O-]", "CC1=C(C(NC(=S)N1)c2cccc(c2)[N+](=O)[O-])C(=O)N(CCO)CCO", "CCCCN(CCCC)C(=O)c1ccccc1N", "Nc1ccccc1C(=O)Nc2ccc3ccccc3c2", "Nc1ccccc1C(=O)Oc2cccc3cccnc23", "CC(C)COC(=O)c1ccccc1N", "CCCCCCCOC(=O)c1ccccc1N", "CC(C)CCOC(=O)c1ccccc1N"];


/***/ }),

/***/ "datagrok-api/dg":
/*!*********************!*\
  !*** external "DG" ***!
  \*********************/
/***/ ((module) => {

module.exports = DG;

/***/ }),

/***/ "datagrok-api/grok":
/*!***********************!*\
  !*** external "grok" ***!
  \***********************/
/***/ ((module) => {

module.exports = grok;

/***/ }),

/***/ "datagrok-api/ui":
/*!*********************!*\
  !*** external "ui" ***!
  \*********************/
/***/ ((module) => {

module.exports = ui;

/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			// no module.id needed
/******/ 			// no module.loaded needed
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId](module, module.exports, __webpack_require__);
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/************************************************************************/
/******/ 	/* webpack/runtime/compat get default export */
/******/ 	(() => {
/******/ 		// getDefaultExport function for compatibility with non-harmony modules
/******/ 		__webpack_require__.n = (module) => {
/******/ 			var getter = module && module.__esModule ?
/******/ 				() => (module['default']) :
/******/ 				() => (module);
/******/ 			__webpack_require__.d(getter, { a: getter });
/******/ 			return getter;
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/define property getters */
/******/ 	(() => {
/******/ 		// define getter functions for harmony exports
/******/ 		__webpack_require__.d = (exports, definition) => {
/******/ 			for(var key in definition) {
/******/ 				if(__webpack_require__.o(definition, key) && !__webpack_require__.o(exports, key)) {
/******/ 					Object.defineProperty(exports, key, { enumerable: true, get: definition[key] });
/******/ 				}
/******/ 			}
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/hasOwnProperty shorthand */
/******/ 	(() => {
/******/ 		__webpack_require__.o = (obj, prop) => (Object.prototype.hasOwnProperty.call(obj, prop))
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/make namespace object */
/******/ 	(() => {
/******/ 		// define __esModule on exports
/******/ 		__webpack_require__.r = (exports) => {
/******/ 			if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
/******/ 				Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
/******/ 			}
/******/ 			Object.defineProperty(exports, '__esModule', { value: true });
/******/ 		};
/******/ 	})();
/******/ 	
/************************************************************************/
var __webpack_exports__ = {};
// This entry needs to be wrapped in an IIFE because it needs to be isolated against other modules in the chunk.
(() => {
/*!************************!*\
  !*** ./src/package.ts ***!
  \************************/
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   _package: () => (/* binding */ _package),
/* harmony export */   autostartBiologics: () => (/* binding */ autostartBiologics),
/* harmony export */   createDemoBiologicsData: () => (/* binding */ createDemoBiologicsData),
/* harmony export */   info: () => (/* binding */ info),
/* harmony export */   populateAdcGlyphs: () => (/* binding */ populateAdcGlyphs)
/* harmony export */ });
/* harmony import */ var datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! datagrok-api/grok */ "datagrok-api/grok");
/* harmony import */ var datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var datagrok_api_dg__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! datagrok-api/dg */ "datagrok-api/dg");
/* harmony import */ var datagrok_api_dg__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(datagrok_api_dg__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _randsmiles__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./randsmiles */ "./src/randsmiles.ts");
/* harmony import */ var _glyphs_glyphs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./glyphs/glyphs */ "./src/glyphs/glyphs.ts");
/* harmony import */ var _datagrok_libraries_db_explorer_src_db_explorer__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! @datagrok-libraries/db-explorer/src/db-explorer */ "../../libraries/db-explorer/src/db-explorer.ts");
/* harmony import */ var _datagrok_libraries_db_explorer_src_renderer__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! @datagrok-libraries/db-explorer/src/renderer */ "../../libraries/db-explorer/src/renderer.ts");
/* harmony import */ var _config__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./config */ "./src/config.ts");
/* eslint-disable max-lines-per-function */
/* eslint-disable max-len */
/* eslint-disable max-lines */
/* Do not change these import lines to match external modules in webpack configuration */



// Placeholder imports for glyph PNGs (ensure your bundler loads them as base64 strings)




const _package = new datagrok_api_dg__WEBPACK_IMPORTED_MODULE_1__.Package();
//name: info
function info() {
    datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.shell.info(_package.webRoot);
}
//name: createDemoBiologicsData
async function createDemoBiologicsData() {
    const cn = 'Biologics:biologics'; // connection name
    const AA = 'ACDEFGHIKLMNPQRSTVWYGGGGGGAAAADEE';
    const randInt = (min, max) => Math.floor(min + Math.random() * (max - min + 1));
    const randPick = (arr) => arr[randInt(0, arr.length - 1)];
    const escape = (s) => s.replace(/'/g, '\'\'').replaceAll('@', '');
    const exec = async (sql) => await datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.data.db.query(cn, sql);
    // 0. Clear existing demo data (keep assay types & target organisms)
    await exec(`
    DELETE FROM biologics.assay_results;
    DELETE FROM biologics.adc;
    DELETE FROM biologics.expression_batches;
    DELETE FROM biologics.purification_batches;
    DELETE FROM biologics.linkers;
    DELETE FROM biologics.drugs;
    DELETE FROM biologics.sequences;
  `); // cascade should handle the rest
    const getConservedRegion = () => {
        const conservedRegionSeq = 'CYSASNITLYCVHQRGGGGAAAAGGGGGGGSSSVTVS';
        let actCons = conservedRegionSeq.substring(0, Math.max(Math.floor(Math.random() * (conservedRegionSeq.length - 1)), Math.floor(conservedRegionSeq.length / 2)));
        const substitutions = randInt(1, 4);
        for (let i = 0; i < substitutions; i++) {
            const pos = randInt(0, actCons.length - 1);
            const newAA = randPick(AA);
            actCons = actCons.substring(0, pos) + newAA + actCons.substring(pos + 1);
        }
        return actCons;
    };
    // 1. Generate 500 protein sequences (~1000 length)
    const seqCount = 200;
    const sequences = [];
    for (let i = 1; i <= seqCount; i++) {
        const len = 950 + randInt(0, 100); // 950-1050
        let s = '';
        for (let j = 0; j < len; j++) {
            s += AA[randInt(0, AA.length - 1)];
            if (j % 150 == 1) {
                s += getConservedRegion();
                j = s.length + 1;
            }
        }
        sequences.push({ name: `Sequence_${i}`, seq: s });
    }
    async function insertSequences(list) {
        const chunkSize = 20;
        const ids = [];
        for (let i = 0; i < list.length; i += chunkSize) {
            const chunk = list.slice(i, i + chunkSize);
            const values = chunk
                .map((c) => `('PROTEIN','${escape(c.seq)}','${escape(c.name)}')`)
                .join(',');
            const sql = `
        INSERT INTO biologics.sequences(sequence_type, sequence, name)
        VALUES ${values}
        RETURNING id, name
      `;
            const df = await exec(sql);
            for (let r = 0; r < df.rowCount; r++)
                ids.push(df.get('id', r));
        }
        return ids;
    }
    // 2. Drugs (empty smiles list for user to fill later)
    const drugSmiles = [
    // { name: 'Drug_1', smiles: 'CCO...' },
    ];
    for (let i = 1; i <= 200; i++)
        drugSmiles.push({ name: `Drug_${i}`, smiles: _randsmiles__WEBPACK_IMPORTED_MODULE_2__.smi[i] });
    async function insertDrugs(list) {
        if (list.length === 0)
            return [];
        const values = list.map((d) => `('${escape(d.name)}','${escape(d.smiles)}')`).join(',');
        const sql = `
      INSERT INTO biologics.drugs(name, smiles)
      VALUES ${values}
      RETURNING id
    `;
        const df = await exec(sql);
        const ids = [];
        for (let r = 0; r < df.rowCount; r++)
            ids.push(df.get('id', r));
        return ids;
    }
    // 3. Linkers (mix of SMALL & PROTEIN)
    const linkerProteinSeqs = ['GGGGS', 'GGSGGS', 'GGGSGGGS', 'GSGSG', 'GPGPG'];
    const linkerSmallSmiles = ['CCOCCOCCOCCOCCOCCOCCOCCO', 'CCNCCNCCNCCNCCNCCNCCN', 'CNC(=O)CCNC(=O)CCNC(=O)CCNC(=O)C', 'NCC(=O)NCC(=O)NCC(=O)NCC=O', 'COCOCOCOCOCOCOCO'];
    const linkers = [];
    for (let i = 0; i < 5; i++)
        linkers.push({ linker_type: 'PROTEIN', linker_sequence: linkerProteinSeqs[i] });
    for (let i = 0; i < 5; i++)
        linkers.push({ linker_type: 'SMALL', linker_molecule_smiles: linkerSmallSmiles[i] });
    async function insertLinkers(list) {
        if (!list.length)
            return [];
        const values = list.map((l) => {
            if (l.linker_type === 'PROTEIN')
                return `('PROTEIN', NULL, '${escape(l.linker_sequence)}')`;
            else
                return `('SMALL', '${escape(l.linker_molecule_smiles)}', NULL)`;
        }).join(',');
        const sql = `
      INSERT INTO biologics.linkers(linker_type, linker_molecule_smiles, linker_sequence)
      VALUES ${values}
      RETURNING id
    `;
        const df = await exec(sql);
        const ids = [];
        for (let r = 0; r < df.rowCount; r++)
            ids.push(df.get('id', r));
        return ids;
    }
    // 4. Fetch existing assay_types & target organisms
    const assayTypesDf = await exec('SELECT id, name FROM biologics.assay_types');
    const assayTypes = [];
    for (let r = 0; r < assayTypesDf.rowCount; r++)
        assayTypes.push({ id: assayTypesDf.get('id', r), name: assayTypesDf.get('name', r) });
    const orgDf = await exec('SELECT id, name FROM biologics.target_organisms');
    const organisms = [];
    for (let r = 0; r < orgDf.rowCount; r++)
        organisms.push(orgDf.get('id', r));
    // 5. Insert everything
    const sequenceIds = await insertSequences(sequences);
    const drugIds = await insertDrugs(drugSmiles);
    const linkerIds = await insertLinkers(linkers);
    // 6. Purification batches (random subset)
    const purBatchesValues = [];
    new Array(sequenceIds.length * 3).fill(null).map((_, i) => randPick(sequenceIds)).forEach((id) => {
        purBatchesValues.push(`(${id}, 'PurBatch_${id}_${Math.floor(Math.random() * 1000)}', 'Auto-generated purification batch')`);
    });
    if (purBatchesValues.length) {
        await exec(`
      INSERT INTO biologics.purification_batches(sequence_id, name, notes)
      VALUES ${purBatchesValues.join(',')}
    `);
    }
    // 7. Expression batches (random subset)
    const exprSystems = ['CHO', 'HEK293', 'E.coli'];
    const exprValues = [];
    new Array(sequenceIds.length * 3).fill(null).map((_, i) => randPick(sequenceIds)).forEach((id) => {
        exprValues.push(`(${id}, '${escape(randPick(exprSystems))}', ${(Math.random() * 150).toFixed(2)}, 'ExprBatch_${id}_${(Math.floor(Math.random() * 1000))}', 'Auto-generated expression batch')`);
    });
    if (exprValues.length) {
        await exec(`
      INSERT INTO biologics.expression_batches(sequence_id, expression_system, yield_mg, name, notes)
      VALUES ${exprValues.join(',')}
    `);
    }
    // 8. ADCs (only if we have drugs)
    let adcCount = 0;
    const adcIds = [];
    if (drugIds.length && linkerIds.length) {
        const adcValues = [];
        const adcToMake = Math.min(drugIds.length * linkerIds.length, sequenceIds.length * 3);
        for (let i = 0; i < adcToMake; i++) {
            const antibodyId = randPick(sequenceIds);
            const drugId = randPick(drugIds);
            const linkerId = randPick(linkerIds);
            const name = `ADC_${i + 1}`;
            const glyph = ''; // placeholder (base64 PNG string)
            adcValues.push(`('${escape(name)}', ${antibodyId}, ${linkerId}, ${drugId}, ${glyph === '' ? 'NULL' : `'${escape(glyph)}'`})`);
        }
        if (adcValues.length) {
            const df = await exec(`
        INSERT INTO biologics.adc(name, antibody_id, linker_id, drug_id, glyph)
        VALUES ${adcValues.join(',')}
        RETURNING id
      `);
            for (let r = 0; r < df.rowCount; r++)
                adcIds.push(df.get('id', r));
            adcCount = adcValues.length;
        }
    }
    // 9. Assay results (random)
    function randomAssayValue(name) {
        if (/IC50/i.test(name))
            return { val: +(Math.random() * 500).toFixed(2), units: 'nM' };
        if (/Caspase/i.test(name))
            return { val: +(Math.random() * 10).toFixed(2), units: 'RFU' };
        if (/half-life/i.test(name))
            return { val: +(Math.random() * 240).toFixed(1), units: 'min' };
        if (/Binding affinity/i.test(name))
            return { val: +(Math.random() * 50).toFixed(2), units: 'nM' };
        if (/Cell binding/i.test(name))
            return { val: +(Math.random() * 100).toFixed(2), units: 'nM' };
        if (/DAR/i.test(name))
            return { val: +(2 + Math.random() * 6).toFixed(2), units: 'ratio' };
        if (/Cmax/i.test(name))
            return { val: +(Math.random() * 100).toFixed(2), units: 'µg/mL' };
        if (/Tmax/i.test(name))
            return { val: +(Math.random() * 48).toFixed(2), units: 'h' };
        if (/AUC/i.test(name))
            return { val: +(Math.random() * 5000).toFixed(1), units: 'µg·h/mL' };
        return { val: +(Math.random() * 100).toFixed(2), units: '' };
    }
    const assayResultChunks = [];
    const maxResults = 5000; // cap
    let inserted = 0;
    for (const seqId of sequenceIds) {
        const perSeqAssays = randInt(3, 26);
        for (let i = 0; i < perSeqAssays; i++) {
            const at = randPick(assayTypes);
            const org = randPick(organisms);
            const av = randomAssayValue(at.name);
            assayResultChunks.push(`(${at.id}, ${randPick(adcIds)}, ${org}, ${av.val}, '${escape(av.units)}')`);
            inserted++;
            if (inserted >= maxResults)
                break;
        }
        if (inserted >= maxResults)
            break;
    }
    if (assayResultChunks.length) {
        // Insert in chunks to avoid overly large statements
        const chunkSize = 200;
        for (let i = 0; i < assayResultChunks.length; i += chunkSize) {
            const part = assayResultChunks.slice(i, i + chunkSize);
            await exec(`
        INSERT INTO biologics.assay_results(assay_id, adc_id, target_organism_id, result_value, units)
        VALUES ${part.join(',')}
      `);
        }
    }
    await populateAdcGlyphs(adcCount + 1);
    return {
        sequences: sequenceIds.length,
        drugs: drugIds.length,
        linkers: linkerIds.length,
        assay_results: inserted,
        purification_batches: purBatchesValues.length,
        expression_batches: exprValues.length,
        adcs: adcCount
    };
}
//name: populateAdcGlyphs
//description: Populates missing ADC glyphs with random PNG (base64) strings
async function populateAdcGlyphs(limit = 50) {
    const cn = 'Biologics:biologics';
    const exec = async (sql) => await datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.data.db.query(cn, sql);
    const escape = (s) => s.replace(/'/g, '\'\'');
    // Use real imported glyph images if available, else fall back
    if (!_glyphs_glyphs__WEBPACK_IMPORTED_MODULE_3__.glyphPool.length)
        throw new Error('No glyph images available');
    // Fetch ADC ids without glyphs (NULL or empty)
    const df = await exec(`SELECT id FROM biologics.adc WHERE (glyph IS NULL OR glyph='') LIMIT ${limit}`);
    if (df.rowCount === 0)
        return { updated: 0 };
    const updates = [];
    for (let i = 0; i < df.rowCount; i++) {
        const id = df.get('id', i);
        const glyph = _glyphs_glyphs__WEBPACK_IMPORTED_MODULE_3__.glyphPool[Math.floor(Math.random() * _glyphs_glyphs__WEBPACK_IMPORTED_MODULE_3__.glyphPool.length)];
        updates.push(`UPDATE biologics.adc SET glyph='${escape(glyph)}' WHERE id=${id}`);
    }
    const batchSize = 20;
    for (let i = 0; i < updates.length; i += batchSize)
        await exec(updates.slice(i, (i + batchSize)).join(';') + ';');
    return { updated: updates.length };
}
//name: autostartbiologics
//tags: autostart
async function autostartBiologics() {
    const exp = _datagrok_libraries_db_explorer_src_db_explorer__WEBPACK_IMPORTED_MODULE_4__.DBExplorer.initFromConfig(_config__WEBPACK_IMPORTED_MODULE_6__.biologicsConfig);
    if (!exp) {
        datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.shell.error('Failed to initialize Biologics DB Explorer');
        return;
    }
    exp.addCustomRenderer((_, colName, value) => {
        const lc = colName?.toLowerCase() || '';
        return (lc === 'structure' || lc.includes('smiles') || lc.includes('compound_structure')) && typeof value === 'string' && datagrok_api_grok__WEBPACK_IMPORTED_MODULE_0__.chem.checkSmiles(value);
    }, (value) => (0,_datagrok_libraries_db_explorer_src_renderer__WEBPACK_IMPORTED_MODULE_5__.moleculeRenderer)(value));
    exp.addCustomRenderer((_, colName, value) => {
        const lc = colName?.toLowerCase() || '';
        return (lc === 'glyph' || lc === 'image' || lc === 'png' || lc === 'thumbnail') && typeof value === 'string' && value.startsWith('iVBORw0KGgo');
    }, (value) => (0,_datagrok_libraries_db_explorer_src_renderer__WEBPACK_IMPORTED_MODULE_5__.rawImageRenderer)(value));
    // exp.addDefaultHeaderReplacerColumns(['units']);
    console.log('Biologics object handlers registered');
}

})();

biologics = __webpack_exports__;
/******/ })()
;
//# sourceMappingURL=package.js.map