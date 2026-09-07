"use strict";
/**
 * Loaded by the platform at startup without the main package bundle. Besides semantic type
 * detectors (none yet), it hosts the autostart that registers the ComputeRun object handler,
 * so saved Compute2 runs (e.g. linked into a space) render and open through their own editor
 * even before any Compute2 function has run. The heavy work (preview, open) goes through
 * DG.Func.byName('Compute2:...'), which loads the main bundle on demand.
 */
class Compute2PackageDetectors extends DG.Package {
    //name: autostart
    //meta.role: autostart
    //description: Registers the Compute2 run handler at platform startup
    autostart() {
        if (DG.ObjectHandler.list().some((h) => h.type === 'ComputeRun'))
            return;
        DG.ObjectHandler.register(new ComputeRunHandler());
    }
}
function openComputeRun(id) {
    return DG.Func.byName('Compute2:OpenWorkflowRun').prepare({ id: id }).call();
}
function computeRunEditorOf(call) {
    const editor = call.func && call.func.options ? call.func.options['editor'] : null;
    return typeof editor === 'string' && editor.startsWith('Compute2:') ? editor : null;
}
/** Same contract as the ModelHandler one entity level up: makes saved Compute2 runs render
 * and open through their editor instead of the platform's generic funccall view. */
class ComputeRunHandler extends DG.ObjectHandler {
    get type() {
        return 'ComputeRun';
    }
    isApplicable(x) {
        try {
            const call = DG.toJs(x);
            return call instanceof DG.FuncCall && computeRunEditorOf(call) != null;
        }
        catch (_) {
            return false;
        }
    }
    async getById(id) {
        return grok.dapi.functions.calls.allPackageVersions().include('session.user,func.package').find(id);
    }
    async refresh(x) {
        return this.getById(DG.toJs(x).id);
    }
    // The run name is its title; a run saved without one falls back to the model name.
    getCaption(x) {
        const call = DG.toJs(x);
        return call.options['title'] || this.modelName(call);
    }
    modelName(call) {
        return (call.func && (call.func.friendlyName || call.func.name)) || 'Run';
    }
    startedText(call) {
        try {
            const started = call.started;
            return started && started.isValid() ? started.format('YYYY-MM-DD HH:mm') : '';
        }
        catch (_) {
            return '';
        }
    }
    subtitleText(call, author) {
        const tags = (call.options['tags'] || []).map((t) => '#' + t).join(' ');
        return [this.modelName(call), this.startedText(call), author, tags].filter((s) => s).join(' · ');
    }
    async authorOf(call) {
        if (call.author && call.author.friendlyName)
            return call.author.friendlyName;
        try {
            const full = await grok.dapi.functions.calls.allPackageVersions()
                .include('session.user,func.package').find(call.id);
            return (full && full.author && full.author.friendlyName) || '';
        }
        catch (_) {
            return '';
        }
    }
    renderIcon(x) {
        const call = DG.toJs(x);
        const icon = ui.iconFA('play-circle');
        // The gallery and the Browse tree consult only renderIcon on JS handlers (caption/markup/card
        // go through the unbridged Dart meta, whose fallback caption is the raw call signature), so
        // the item is reworked from inside the returned icon: run title as the name, plus
        // model/date/author/tags as the gallery subtitle. Unrecognized hosts are left untouched.
        void this.patchHostItem(icon, call);
        return icon;
    }
    // The icon is returned detached and inserted later, so poll until it lands in a known host;
    // re-assert once afterwards in case the host rewrites the caption on a late render pass.
    async patchHostItem(icon, call) {
        const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
        for (let attempt = 0; attempt < 20; attempt++) {
            await sleep(attempt === 0 ? 0 : 100);
            if (!icon.isConnected)
                continue;
            if (this.tryPatchOnce(icon, call))
                setTimeout(() => this.tryPatchOnce(icon, call), 1000);
            return;
        }
    }
    tryPatchOnce(icon, call) {
        const treeItem = icon.closest('.d4-tree-view-item');
        if (treeItem) {
            const label = treeItem.querySelector('.d4-tree-view-item-label');
            if (!label)
                return false;
            for (const n of Array.from(label.childNodes)) {
                if (n.nodeType === Node.TEXT_NODE)
                    n.remove();
            }
            const started = this.startedText(call);
            label.appendChild(document.createTextNode(this.getCaption(call) + (started ? ' (' + started + ')' : '')));
            return true;
        }
        const linkLabel = icon.closest('.d4-link-label');
        if (!linkLabel || !linkLabel.closest('.grok-gallery-grid'))
            return false;
        const label = linkLabel.querySelector('label');
        if (label)
            label.textContent = this.getCaption(call);
        const subtitle = Array.from(linkLabel.children)
            .find((el) => el.tagName === 'SPAN' && el !== icon.parentElement);
        if (subtitle) {
            subtitle.textContent = ' ' + this.subtitleText(call, (call.author && call.author.friendlyName) || '');
            void this.authorOf(call).then((author) => {
                subtitle.textContent = ' ' + this.subtitleText(call, author);
            });
        }
        return true;
    }
    renderMarkup(x) {
        const call = DG.toJs(x);
        const markup = ui.span([ui.iconFA('play-circle'), ui.label(this.getCaption(call))]);
        markup.ondblclick = () => void openComputeRun(call.id);
        return markup;
    }
    renderTooltip(x) {
        const call = DG.toJs(x);
        return ui.divV([
            ui.divText(this.getCaption(call)),
            ui.divText(this.subtitleText(call, (call.author && call.author.friendlyName) || '')),
        ]);
    }
    renderCard(x) {
        const call = DG.toJs(x);
        const subtitle = ui.divText(this.subtitleText(call, (call.author && call.author.friendlyName) || ''), { style: { color: 'var(--grey-4)' } });
        void this.authorOf(call).then((author) => {
            subtitle.textContent = this.subtitleText(call, author);
        });
        const card = ui.divV([
            ui.divH([ui.iconFA('play-circle'), ui.divText(this.getCaption(call), { style: { fontWeight: 'bold' } })], { style: { alignItems: 'center', gap: '6px' } }),
            subtitle,
        ], 'd4-gallery-item');
        card.ondblclick = () => void openComputeRun(call.id);
        return card;
    }
    renderProperties(x) {
        const call = DG.toJs(x);
        const info = ui.div([]);
        const fill = (author) => {
            ui.empty(info);
            const map = { 'Model': this.modelName(call), 'Started': this.startedText(call) };
            if (author)
                map['Author'] = author;
            const tags = (call.options['tags'] || []).join(', ');
            if (tags)
                map['Tags'] = tags;
            info.appendChild(ui.tableFromMap(map));
        };
        fill((call.author && call.author.friendlyName) || '');
        void this.authorOf(call).then(fill);
        return ui.divV([
            ui.h2(this.getCaption(call)),
            info,
            ui.div([ui.bigButton('Open run', () => void openComputeRun(call.id))], { style: { marginTop: '8px' } }),
        ]);
    }
    async renderPreview(x) {
        const call = DG.toJs(x);
        const editorCall = call.func.prepare({});
        editorCall.aux.initialRunId = call.id;
        // Return the editor view directly so the docked preview is the view RFV/TreeWizard controls
        // (wrapping it makes a separate view the tracked preview — see compute-utils ModelHandler).
        const res = await DG.Func.byName(computeRunEditorOf(call)).prepare({ call: editorCall }).call();
        return res.getOutputParamValue();
    }
    init() {
        this.registerParamFunc('Open run', (x) => void openComputeRun(DG.toJs(x).id));
    }
}
