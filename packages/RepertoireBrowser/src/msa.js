import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export class MsaMethods {

    init(view, inputs) {
        this.view = view;
        this.table = view.table;
        this.msaInput = inputs.msaContentChoice;
        this.msaHostL = inputs.msa_host_L;
        this.msaHostH = inputs.msa_host_H;

        this.msaOpts = {
            el: this.msaHostL,
            vis: {
                conserv: false,
                overviewbox: false,
                consensus: true,
                seqlogo: true,
                scaleslider: false,
            },
            conf: {
                dropImport: true
            },
            bootstrapMenu: false,
            zoomer: {
                boxRectHeight: 1,
                boxRectWidth: 1,
                labelNameLength: 110,
                labelFontsize: 10,
                labelIdLength: 20
            }
        };
        this.msaL = new msa.msa(this.msaOpts);
        this.msaOpts.el = this.msaHostH;
        this.msaH = new msa.msa(this.msaOpts);
        this.gffParser = msa.io.gff;

        this.seqHeavyCol = this.table.col('sequence_alignment_heavy');
        this.seqLightCol = this.table.col('sequence_alignment_light');
        this.germHeavyCol = this.table.col('germline_alignment_heavy');
        this.germLightCol = this.table.col('germline_alignment_light');

        this.seqAlignHeavyAACol = this.table.col('sequence_alignment_aa_heavy');
        this.seqAlignLightAACol = this.table.col('sequence_alignment_aa_light');
        this.germAlignHeavyAACol = this.table.col('germline_alignment_aa_heavy');
        this.germAlignLightAACol = this.table.col('germline_alignment_aa_light');

        this.vStartHeavy = this.table.col('v_alignment_start_heavy');
        this.dStartHeavy = this.table.col('d_alignment_start_heavy');
        this.jStartHeavy = this.table.col('j_alignment_start_heavy');
        this.vEndHeavy = this.table.col('v_alignment_end_heavy');
        this.dEndHeavy = this.table.col('d_alignment_end_heavy');
        this.jEndHeavy = this.table.col('j_alignment_end_heavy');

        this.vStartLight = this.table.col('v_alignment_start_light');
        this.dStartLight = this.table.col('d_alignment_start_light');
        this.jStartLight = this.table.col('j_alignment_start_light');
        this.vEndLight = this.table.col('v_alignment_end_light');
        this.dEndLight = this.table.col('d_alignment_end_light');
        this.jEndLight = this.table.col('j_alignment_end_light');

        this.drawAlignments();
    }

    msaRender(m, msaFasta, gffAnnotations = {}) {
        const seqs = msa.io.fasta.parse(msaFasta);
        m.seqs.features = {};
        m.seqs.reset(seqs);
        m.seqs.removeAllFeatures();
        for (let annotation of Object.values(gffAnnotations)) {
            const features = this.gffParser.parseSeqs(annotation);
            m.seqs.addFeatures(features);
        }
        m.render();
    }

    gffAnnotator(seqid, source, type, start, end, score, strand, phase, attributes) {
        return `${seqid}\t${source}\t${type}\t${start}\t${end}\t${score}\t${strand}\t${phase}\t${attributes}\n`;
    }

    annotateRegions(startCol, endCol, i, chain, name, color) {
        let an = '';
        if (startCol && endCol) {
            an += this.gffAnnotator(`germline_align_${chain}`, '.', 'gene',
            startCol.get(i), endCol.get(i), '.', '+', '.', `Name=${name};Color=${color}`);
        }
        return an;
    }

    /** Translates an AA sequence into GFF formatted annotations. */
    getAA_gff(molId, aaStr) {
        let gff = "##gff-version 3\n";
        let pos = 2;
        for (i = 0; i < aaStr.length; i++) {
            let aaBase = aaStr[i];
            let line = `${molId} . p	${pos} 	${pos}	.	.  +  Name=${aaBase};Color=gray\n`;
            gff = gff + line;
            pos = pos + 3;
        }
        return gff;
    }

    /** Formats AA column values into GFF annotation. */
    makeColGFF(col, rowIdx, seqId) {
        let gff = null;
        if (col) {
            const seq = col.get(rowIdx);
            gff = this.getAA_gff(seqId, seq);
        }
        return gff;
    }

    msaOnCols(ids, cols, idx, m, opts = {}) {
        if (cols.some(col => col == null) || ids.length != cols.length) return;
        const toFasta = (id, seq) => `>${id}\n${seq}\n`;
        let seqs = '';
        for (let i = 0; i < ids.length; i++) {
            seqs += toFasta(ids[i], cols[i].get(idx));
        }
        this.msaRender(m, seqs, opts);
    }

    drawAlignments() {
        const mode = this.msaInput.value;
        const idx = this.table.currentRow.idx;

        if (mode === 'AA MSA') {
            this.msaOnCols(['seq_align_heavy_aa', 'germline_align_heavy_aa'], [this.seqAlignHeavyAACol, this.germAlignHeavyAACol], idx, this.msaH);
            this.msaOnCols(['seq_align_light_aa', 'germline_align_light_aa'], [this.seqAlignLightAACol, this.germAlignLightAACol], idx, this.msaL);
        } else if (mode === 'DNA MSA') {
            this.msaOnCols(['seq_align_heavy', 'germline_align_heavy'], [this.seqHeavyCol, this.germHeavyCol], idx, this.msaH);
            this.msaOnCols(['seq_align_light', 'germline_align_light'], [this.seqLightCol, this.germLightCol], idx, this.msaL);
        } else if (mode === 'Hybrid') {
            let annotsH = {};
            let headerLen = 16;

            let gffGermH = '##gff-version 3\n';
            gffGermH += this.annotateRegions(this.vStartHeavy, this.vEndHeavy, idx, 'heavy', 'V region', 'violet');
            gffGermH += this.annotateRegions(this.dStartHeavy, this.dEndHeavy, idx, 'heavy', 'D region', 'green');
            gffGermH += this.annotateRegions(this.jStartHeavy, this.jEndHeavy, idx, 'heavy', 'J region', 'blue');
            if (gffGermH.length > headerLen) annotsH['gffGerm'] = gffGermH;

            let gffAASeqH = this.makeColGFF(this.seqAlignHeavyAACol, idx, "seq_align_heavy");
            let gffAAGermH = this.makeColGFF(this.germAlignHeavyAACol, idx, "germline_align_heavy");
            if (gffAASeqH.length > headerLen) annotsH['gffAASeq'] = gffAASeqH;
            if (gffAAGermH.length > headerLen) annotsH['gffAAGerm'] = gffAAGermH;

            this.msaOnCols(['seq_align_heavy', 'germline_align_heavy'], [this.seqHeavyCol, this.germHeavyCol], idx, this.msaH, annotsH);

            let annotsL = {};
            let gffGermL = '##gff-version 3\n';
            gffGermL += this.annotateRegions(this.vStartLight, this.vEndLight, idx, 'light', 'V region', 'violet');
            gffGermL += this.annotateRegions(this.dStartLight, this.dEndLight, idx, 'light', 'D region', 'green');
            gffGermL += this.annotateRegions(this.jStartLight, this.jEndLight, idx, 'light', 'J region', 'blue');
            if (gffGermL.length > headerLen) annotsL['gffGerm'] = gffGermL;

            let gffAASeqL = this.makeColGFF(this.seqAlignLightAACol, idx, "seq_align_light");
            let gffAAGermL = this.makeColGFF(this.germAlignLightAACol, idx, "germline_align_light");
            if (gffAASeqL.length > headerLen) annotsL['gffAASeq'] = gffAASeqL;
            if (gffAAGermL.length > headerLen) annotsL['gffAAGerm'] = gffAAGermL;

            this.msaOnCols(['seq_align_light', 'germline_align_light'], [this.seqLightCol, this.germLightCol], idx, this.msaL, annotsL);
        }
    }

}