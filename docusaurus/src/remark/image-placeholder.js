const visit = require('unist-util-visit');
const path = require('path');
const fs = require('fs')
const logger = require('@docusaurus/logger')

const plugin = () => {
    const isUrl = urlString => {
        try {
            return Boolean(new URL(urlString));
        } catch (e) {
            return false;
        }
    }
    const transformer = (root, file) => {
        visit(root, 'image', (node) => {
            if (!node.url || (!fs.existsSync(path.join(path.dirname(file.history[0]), node.url)) && !isUrl(node.url))) {
                logger.report(
                    'warn',
                    )`Image couldn't be resolved: (url=${node.url}) in path=${file.history[0]}`;
                node.url = "/img/image_placeholder.png";
            }
        });
    };
    return transformer;
};

module.exports = plugin;
