// @ts-check
// Sidebar for docusaurus-debug.config.js, whose typedoc output is ./js-api (typedoc-sidebar.js reads ./api).
/** @type {import('@docusaurus/plugin-content-docs').SidebarsConfig} */
module.exports = {
    api: [
        {
            type: 'category',
            label: 'JavaScript API',
            link: { type: 'doc', id: 'index' },
            items: require('./js-api/typedoc-sidebar.cjs'),
        },
    ],
};
