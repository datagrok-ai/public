import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export default class MsaMethods {

    // ---- MSA ----
    static msaRender(m, msa_fasta, gff_annots = null,
                       gff_aa_annots1 = null,
                       gff_aa_annots2 = null) {
        const seqs = msa.io.fasta.parse(msa_fasta);
        m.seqs.reset(seqs);
        m.seqs.removeAllFeatures();
        if (gff_aa_annots1) {
            const features = gffParser.parseSeqs(gff_aa_annots1);
            m.seqs.addFeatures(features);
        }
        if (gff_aa_annots2) {
            const features = gffParser.parseSeqs(gff_aa_annots2);
            m.seqs.addFeatures(features);
        }
        if (gff_annots) {
            const features = gffParser.parseSeqs(gff_annots);
            m.seqs.addFeatures(features);
        }
        m.render();
    }

    static gffAnnotator(seqid, source, type, start, end, score, strand, phase, attributes) {
        return `${seqid}\t${source}\t${type}\t${start}\t${end}\t${score}\t${strand}\t${phase}\t${attributes}\n`;
    }

    // Translate AA sequence into a GFF formatted annotations
    static getAA_gff(molId, aaStr) {
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

    // Format AA column values into GFF annotation
    static makeColGFF(col, rowIdx, seqId) {
        let gff = null;
        if (col) {
            const seq = col.get(rowIdx);
            gff = getAA_gff(seqId, seq);
        }
        return gff;
    }

    static drawAlignments() {
        const idx = table.currentRow.idx;
        if (seqHeavyCol && germHeavyCol) {
            const seqsH = `>seq_align_heavy\n${seqHeavyCol.get(idx)}\n>germline_align_heavy\n${germHeavyCol.get(idx)}\n`;

            const gffAnnotsH = '##gff-version 3\n' + (
                (vStartHeavy && vEndHeavy) ? gffAnnotator('germline_align_heavy', '.', 'gene',
                    vStartHeavy.get(idx), vEndHeavy.get(idx), '.', '+', '.', 'Name=V region;Color=violet') : '') + (
                (dStartHeavy && dEndHeavy) ? gffAnnotator('germline_align_heavy', '.', 'gene',
                    dStartHeavy.get(idx), dEndHeavy.get(idx), '.', '+', '.', 'Name=D region;Color=green') : '') + (
                (jStartHeavy && jEndHeavy) ? gffAnnotator('germline_align_heavy', '.', 'gene',
                    jStartHeavy.get(idx), jEndHeavy.get(idx), '.', '+', '.', 'Name=J region;Color=blue') : '');

            let gffAASeq = makeColGFF(seqAlignHeavyAACol, idx, "seq_align_heavy");
            let gffAAGerm = makeColGFF(germAlignHeavyAACol, idx, "germline_align_heavy");

            msaRender(msaH, seqsH,
                gffAnnotsH.length > 16 ? gffAnnotsH : null,
                gffAASeq.length > 16   ? gffAASeq : null,
                gffAAGerm.length > 16  ? gffAAGerm : null
            );
        }
        if (seqLightCol && germLightCol) {
            const seqsL = `>seq_align_light\n${seqLightCol.get(idx)}\n>germline_align_light\n${germLightCol.get(idx)}\n`;

            const gffAnnotsL = '##gff-version 3\n' + (
                (vStartLight && vEndLight) ? gffAnnotator('germline_align_light', '.', 'gene',
                    vStartLight.get(idx), vEndLight.get(idx), '.', '+', '.', 'Name=V region;Color=violet') : '') + (
                (dStartLight && dEndLight) ? gffAnnotator('germline_align_light', '.', 'gene',
                    dStartLight.get(idx), dEndLight.get(idx), '.', '+', '.', 'Name=D region;Color=green') : '') + (
                (jStartLight && jEndLight) ? gffAnnotator('germline_align_light', '.', 'gene',
                    jStartLight.get(idx), jEndLight.get(idx), '.', '+', '.', 'Name=J region;Color=blue') : '');

            let gffAASeq = makeColGFF(seqAlignLightAACol, idx, "seq_align_light");
            let gffAAGerm = makeColGFF(germAlignLightAACol, idx, "germline_align_light");

            msaRender(msaL, seqsL,
                gffAnnotsL.length > 16 ? gffAnnotsL : null,
                gffAASeq.length > 16   ? gffAASeq : null,
                gffAAGerm.length > 16  ? gffAAGerm : null
            );
        }
    }

}