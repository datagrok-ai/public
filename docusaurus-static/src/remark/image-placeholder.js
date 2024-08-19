import {visit} from 'unist-util-visit';
import {join, dirname} from 'path';
import {existsSync} from 'fs';
import logger from '@docusaurus/logger';

const plugin = (options) => {
    const isUrl = urlString => {
        try {
            return Boolean(new URL(urlString));
        } catch (e) {
            return false;
        }
    }
    const transformer = async (tree, file) => {
        visit(tree, 'image', (node) => {
            if (!node.url || (!existsSync(join(dirname(file.history[0]), node.url)) && !isUrl(node.url))) {
                logger.warn`Image couldn't be resolved: (url=${node.url}) in path=${file.history[0]}`;
                node.url = "/img/image_placeholder.png";
            }
        });
    };
    return transformer;
};

export default plugin;
