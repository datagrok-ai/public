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
    // The run name: the legacy title when set, otherwise the friendly name, else the model name.
    getCaption(x) {
        const call = DG.toJs(x);
        return call.options['title'] || call.friendlyName || this.modelName(call);
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
    renderIcon(_x) {
        return ui.iconFA('play-circle');
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
