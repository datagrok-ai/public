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

    msaRender(m, msa_fasta, {
        gff_annots = null,
        gff_aa_annots1 = null,
        gff_aa_annots2 = null } = {}) {
        const seqs = msa.io.fasta.parse(msa_fasta);
        m.seqs.reset(seqs);
        // m.seqs.removeAllFeatures();
        // if (gff_aa_annots1) {
        //     const features = gffParser.parseSeqs(gff_aa_annots1);
        //     m.seqs.addFeatures(features);
        // }
        // if (gff_aa_annots2) {
        //     const features = gffParser.parseSeqs(gff_aa_annots2);
        //     m.seqs.addFeatures(features);
        // }
        // if (gff_annots) {
        //     const features = gffParser.parseSeqs(gff_annots);
        //     m.seqs.addFeatures(features);
        // }
        m.render();
    }

    gffAnnotator(seqid, source, type, start, end, score, strand, phase, attributes) {
        return `${seqid}\t${source}\t${type}\t${start}\t${end}\t${score}\t${strand}\t${phase}\t${attributes}\n`;
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
            if (this.seqHeavyCol && this.germHeavyCol) {
                const seqsH = `>seq_align_heavy\n${this.seqHeavyCol.get(idx)}\n>germline_align_heavy\n${this.germHeavyCol.get(idx)}\n`;
    
                const gffAnnotsH = '##gff-version 3\n' + (
                    (this.vStartHeavy && this.vEndHeavy) ? this.gffAnnotator('germline_align_heavy', '.', 'gene',
                    this.vStartHeavy.get(idx), this.vEndHeavy.get(idx), '.', '+', '.', 'Name=V region;Color=violet') : '') + (
                    (this.dStartHeavy && this.dEndHeavy) ? this.gffAnnotator('germline_align_heavy', '.', 'gene',
                    this.dStartHeavy.get(idx), this.dEndHeavy.get(idx), '.', '+', '.', 'Name=D region;Color=green') : '') + (
                    (this.jStartHeavy && this.jEndHeavy) ? this.gffAnnotator('germline_align_heavy', '.', 'gene',
                    this.jStartHeavy.get(idx), this.jEndHeavy.get(idx), '.', '+', '.', 'Name=J region;Color=blue') : '');
    
                let gffAASeq = this.makeColGFF(this.seqAlignHeavyAACol, idx, "seq_align_heavy");
                let gffAAGerm = this.makeColGFF(this.germAlignHeavyAACol, idx, "germline_align_heavy");
    
                this.msaRender(this.msaH, seqsH, {
                    gffAnnotsH: gffAnnotsH.length > 16 ? gffAnnotsH : null,
                    gffAASeq: gffAASeq.length > 16   ? gffAASeq : null,
                    gffAAGerm: gffAAGerm.length > 16  ? gffAAGerm : null
                });
            }
            if (this.seqLightCol && this.germLightCol) {
                const seqsL = `>seq_align_light\n${this.seqLightCol.get(idx)}\n>germline_align_light\n${this.germLightCol.get(idx)}\n`;
    
                const gffAnnotsL = '##gff-version 3\n' + (
                    (this.vStartLight && this.vEndLight) ? this.gffAnnotator('germline_align_light', '.', 'gene',
                    this.vStartLight.get(idx), this.vEndLight.get(idx), '.', '+', '.', 'Name=V region;Color=violet') : '') + (
                    (this.dStartLight && this.dEndLight) ? this.gffAnnotator('germline_align_light', '.', 'gene',
                    this.dStartLight.get(idx), this.dEndLight.get(idx), '.', '+', '.', 'Name=D region;Color=green') : '') + (
                    (this.jStartLight && this.jEndLight) ? this.gffAnnotator('germline_align_light', '.', 'gene',
                    this.jStartLight.get(idx), this.jEndLight.get(idx), '.', '+', '.', 'Name=J region;Color=blue') : '');
    
                let gffAASeq = this.makeColGFF(this.seqAlignLightAACol, idx, "seq_align_light");
                let gffAAGerm = this.makeColGFF(this.germAlignLightAACol, idx, "germline_align_light");
    
                this.msaRender(this.msaL, seqsL, {
                    gffAnnotsL: gffAnnotsL.length > 16 ? gffAnnotsL : null,
                    gffAASeq: gffAASeq.length > 16   ? gffAASeq : null,
                    gffAAGerm: gffAAGerm.length > 16  ? gffAAGerm : null
                });
            }
        }
    }

}