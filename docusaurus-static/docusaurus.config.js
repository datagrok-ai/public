// @ts-check
// Note: type annotations allow type checking and IDEs autocompletion

const router = process.env.DOCUSAURUS_ROUTER;

/** @type {import('@docusaurus/types').Config} */
const config = {
    title: 'Datagrok',
    tagline: 'Datagrok: Swiss Army Knife for Data',
    url: 'http://localhost:8080',
    baseUrl: '/',
    onBrokenLinks: 'throw',
    onBrokenMarkdownLinks: 'throw',
    onBrokenAnchors: 'throw',
    onDuplicateRoutes: 'throw',
    favicon: 'favicon/favicon.ico',
    staticDirectories: ['static'],

    // Even if you don't use internalization, you can use this field to set useful
    // metadata like html lang. For example, if your site is Chinese, you may want
    // to replace "en" with "zh-Hans".
    i18n: {
        defaultLocale: 'en', locales: ['en']
    },

    future: {
        experimental_router: router,
    },

    presets: [['classic', /** @type {import('@docusaurus/preset-classic').Options} */
        ({
            docs: {
                sidebarPath: require.resolve('./sidebar-empty.js'),
                path: '../help',
                routeBasePath: 'help',
                exclude: ['**/_*/**', '_*/**', '**/_*', '**/*-test.md']
            }
        }),],],

    themeConfig: {
        colorMode: {
            defaultMode: 'light', disableSwitch: true
        },
        customCss: require.resolve('./src/css/custom.css'),
    }
};

module.exports = config;
