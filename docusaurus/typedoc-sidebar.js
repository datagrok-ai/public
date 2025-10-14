// @ts-check
/** @type {import('@docusaurus/plugin-content-docs').SidebarsConfig} */
const sidebars = {
    api: [
        {
            type: 'doc',
            id: 'index',
            label: 'API Overview',
        },
        {
            type: 'category',
            label: 'JavaScript API',
            link: { type: "doc", id: "js/index" },
            items: require('./api/js/typedoc-sidebar'),
        },
        {
            type: 'category',
            label: 'Python API',
            link: { type: "doc", id: "py/index" },
            items: require('./api/py/typedoc-sidebar'),
        },
    ],
};

module.exports = sidebars;