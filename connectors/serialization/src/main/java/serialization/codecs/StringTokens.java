package serialization.codecs;

import java.util.ArrayList;
import java.util.List;

import serialization.BufferAccessor;
import serialization.IntColumn;
import serialization.StringColumn;

// Port of StringTokens (ddt/lib/src/serialization/string_tokens.dart, GROK-20761); same
// structure, names and traversal order, so the two implementations stay diffable.
// Splits sorted string categories that share one template into literal / int / string parts,
// so that each part can be encoded by the codec that fits it.
// See core/docs/features/column-encoders-v2/DESIGN.md, section 5.
public class StringTokens {
    public static final int VERSION = 1;
    public static final int MAX_INT_DIGITS = 9;
    public static final int MAX_LONG_DIGITS = 15;
    public static final int MAX_HEX_DIGITS = 7;
    public static final long LONG_SPLIT = 1000000000L;

    static final int _DIGIT = 0;
    static final int _LETTER = 1;
    static final int _SEPARATOR = 2;

    public List<StringTokenPart> parts = new ArrayList<>();
    public int emptyCatIdx = -1;
    public int catCount = 0;
    public int valueCount = 0;

    private List<TokenSegment> segments;
    private List<FieldStats> fields;
    private List<List<StringTokenPart>> fieldParts;
    private int[] runBounds;

    static int classOf(int c) {
        if (c >= 48 && c <= 57)
            return _DIGIT;
        if (c < 128 && !(c >= 65 && c <= 90) && !(c >= 97 && c <= 122))
            return _SEPARATOR;
        return _LETTER;
    }

    public String template() {
        StringBuilder sb = new StringBuilder();
        for (StringTokenPart part : parts)
            sb.append(part);
        return sb.toString();
    }

    // Returns null when [categories] do not share one template, or when no int part would result.
    public static StringTokens analyze(List<String> categories) {
        StringTokens t = new StringTokens();
        t.catCount = categories.size();
        t.emptyCatIdx = categories.indexOf("");
        t.valueCount = t.catCount - (t.emptyCatIdx >= 0 ? 1 : 0);
        if (t.valueCount < 2)
            return null;

        t.buildTemplate(categories.get(t.emptyCatIdx == 0 ? 1 : 0));
        if (t.fields.isEmpty())
            return null;

        for (int i = 0; i < categories.size(); i++)
            if (i != t.emptyCatIdx && !t.scan(categories.get(i)))
                return null;

        t.classify();
        boolean hasIntPart = false;
        for (StringTokenPart p : t.parts)
            hasIntPart |= p.kind == StringTokenPart.INT || p.kind == StringTokenPart.LONG_INT;
        if (!hasIntPart)
            return null;

        t.extract(categories);
        return t;
    }

    private void buildTemplate(String first) {
        segments = new ArrayList<>();
        fields = new ArrayList<>();
        int pos = 0;
        while (pos < first.length()) {
            boolean sep = classOf(first.charAt(pos)) == _SEPARATOR;
            int end = pos;
            while (end < first.length() && (classOf(first.charAt(end)) == _SEPARATOR) == sep)
                end++;
            if (sep)
                segments.add(TokenSegment.sep(first.substring(pos, end)));
            else {
                segments.add(TokenSegment.ofField(fields.size()));
                fields.add(new FieldStats());
            }
            pos = end;
        }
    }

    private boolean scan(String s) {
        int pos = 0;
        int len = s.length();
        for (TokenSegment seg : segments) {
            if (seg.isLiteral()) {
                String lit = seg.literal;
                if (pos + lit.length() > len)
                    return false;
                for (int k = 0; k < lit.length(); k++)
                    if (s.charAt(pos + k) != lit.charAt(k))
                        return false;
                pos += lit.length();
            }
            else {
                int end = fields.get(seg.field).addAndScan(s, pos);
                if (end == pos)
                    return false;
                pos = end;
            }
        }
        return pos == len;
    }

    private void classify() {
        fieldParts = new ArrayList<>();
        for (int i = 0; i < fields.size(); i++)
            fieldParts.add(new ArrayList<>());
        int maxRuns = 1;
        for (TokenSegment seg : segments) {
            if (seg.isLiteral()) {
                StringTokenPart p = new StringTokenPart(StringTokenPart.LITERAL);
                p.text = seg.literal;
                parts.add(p);
                continue;
            }
            FieldStats st = fields.get(seg.field);
            if (st.consistent) {
                maxRuns = Math.max(maxRuns, st.sig.size());
                for (int r = 0; r < st.sig.size(); r++) {
                    if (st.sig.get(r) == _DIGIT)
                        addIntParts(seg.field, r, st.runWidthMin.get(r), st.runWidthMax.get(r), st.runLeadingZero.get(r));
                    else {
                        StringTokenPart p = new StringTokenPart(StringTokenPart.STRING);
                        p.run = r;
                        addPart(seg.field, p);
                    }
                }
            }
            else if ((st.hexUpper || st.hexLower) && st.widthMin == st.widthMax && st.widthMax <= MAX_HEX_DIGITS) {
                StringTokenPart p = new StringTokenPart(StringTokenPart.INT);
                p.flags = StringTokenPart.FLAG_PADDED
                        | (st.hexUpper ? StringTokenPart.FLAG_HEX_UPPER : StringTokenPart.FLAG_HEX_LOWER);
                p.width = st.widthMax;
                addPart(seg.field, p);
            }
            else
                addPart(seg.field, new StringTokenPart(StringTokenPart.STRING));
        }
        runBounds = new int[maxRuns + 1];
    }

    private void addPart(int field, StringTokenPart part) {
        parts.add(part);
        fieldParts.get(field).add(part);
    }

    private void addIntParts(int field, int run, int wMin, int wMax, boolean leadingZero) {
        if (wMin == wMax && wMax > MAX_INT_DIGITS && wMax <= MAX_LONG_DIGITS) {
            StringTokenPart p = new StringTokenPart(StringTokenPart.LONG_INT);
            p.flags = StringTokenPart.FLAG_PADDED;
            p.width = wMax;
            p.run = run;
            addPart(field, p);
            return;
        }
        if (wMin == wMax) {
            int chunks = (wMax + MAX_INT_DIGITS - 1) / MAX_INT_DIGITS;
            int offset = 0;
            for (int c = 0; c < chunks; c++) {
                int w = c == 0 ? wMax - MAX_INT_DIGITS * (chunks - 1) : MAX_INT_DIGITS;
                StringTokenPart p = new StringTokenPart(StringTokenPart.INT);
                p.flags = StringTokenPart.FLAG_PADDED;
                p.width = w;
                p.run = run;
                p.chunkOffset = offset;
                addPart(field, p);
                offset += w;
            }
        }
        else if (!leadingZero && wMax <= MAX_INT_DIGITS) {
            StringTokenPart p = new StringTokenPart(StringTokenPart.INT);
            p.run = run;
            addPart(field, p);
        }
        else {
            StringTokenPart p = new StringTokenPart(StringTokenPart.STRING);
            p.run = run;
            addPart(field, p);
        }
    }

    private void extract(List<String> categories) {
        for (StringTokenPart part : parts) {
            if (part.kind == StringTokenPart.INT)
                part.ints = new int[valueCount];
            else if (part.kind == StringTokenPart.LONG_INT)
                part.longs = new long[valueCount];
            else if (part.kind == StringTokenPart.STRING)
                part.strings = new String[valueCount];
        }

        int j = 0;
        for (int i = 0; i < categories.size(); i++) {
            if (i == emptyCatIdx)
                continue;
            String s = categories.get(i);
            int pos = 0;
            for (TokenSegment seg : segments) {
                if (seg.isLiteral()) {
                    pos += seg.literal.length();
                    continue;
                }
                List<StringTokenPart> fp = fieldParts.get(seg.field);
                int end = pos;
                if (fields.get(seg.field).consistent) {
                    int runs = 0;
                    runBounds[0] = pos;
                    int prev = -1;
                    for (; end < s.length(); end++) {
                        int c = classOf(s.charAt(end));
                        if (c == _SEPARATOR)
                            break;
                        if (prev >= 0 && c != prev)
                            runBounds[++runs] = end;
                        prev = c;
                    }
                    runBounds[runs + 1] = end;
                    for (StringTokenPart part : fp) {
                        int rs = runBounds[part.run];
                        int re = runBounds[part.run + 1];
                        if (part.kind == StringTokenPart.STRING)
                            part.strings[j] = s.substring(rs, re);
                        else if (part.kind == StringTokenPart.LONG_INT)
                            part.longs[j] = parseLong(s, rs, re);
                        else if (part.chunkOffset >= 0)
                            part.ints[j] = (int) parseLong(s, rs + part.chunkOffset, rs + part.chunkOffset + part.width);
                        else
                            part.ints[j] = (int) parseLong(s, rs, re);
                    }
                }
                else {
                    while (end < s.length() && classOf(s.charAt(end)) != _SEPARATOR)
                        end++;
                    StringTokenPart part = fp.get(0);
                    if (part.kind == StringTokenPart.STRING)
                        part.strings[j] = s.substring(pos, end);
                    else
                        part.ints[j] = parseHex(s, pos, end);
                }
                pos = end;
            }
            j++;
        }

        // A string part whose values are all identical collapses into a literal (Dart does this
        // when building the nested part columns).
        for (StringTokenPart part : parts) {
            if (part.kind != StringTokenPart.STRING)
                continue;
            boolean constant = true;
            for (int i = 1; i < part.strings.length && constant; i++)
                constant = part.strings[i].equals(part.strings[0]);
            if (constant) {
                part.kind = StringTokenPart.LITERAL;
                part.text = part.strings[0];
                part.strings = null;
            }
        }
    }

    private static long parseLong(String s, int start, int end) {
        long v = 0;
        for (int i = start; i < end; i++)
            v = v * 10 + (s.charAt(i) - 48);
        return v;
    }

    private static int parseHex(String s, int start, int end) {
        int v = 0;
        for (int i = start; i < end; i++) {
            int c = s.charAt(i);
            v = v * 16 + (c <= 57 ? c - 48 : c <= 70 ? c - 55 : c - 87);
        }
        return v;
    }

    // ---- writer side ----

    // Payload size (excluding the shared row->category index column, which costs the same for
    // every string encoder and cancels out of the id-4 vs id-0 comparison). Builds the nested
    // part columns; encode() reuses them.
    public int estimate() {
        int size = 13;
        for (StringTokenPart part : parts) {
            size += 1;
            if (part.kind == StringTokenPart.LITERAL)
                size += 10 + part.text.length();
            else if (part.kind == StringTokenPart.LONG_INT)
                size += 2 + estimateLongPart(part);
            else if (part.kind == StringTokenPart.INT) {
                part.intColumn = new IntColumn("", part.ints);
                size += 6 + 4 + IntColumn.Encoding.choose(part.ints, part.ints.length).size;
            }
            else {
                part.stringColumn = new StringColumn("", part.strings);
                size += 4 + part.stringColumn.estimateCategoriesPayload();
            }
        }
        return size;
    }

    // Estimates the cheaper of the two long-part payload modes and stores the choice on [part].
    // Mode 0: Int64 first, Int64 deltaBase, BitPacking of (delta - deltaBase) over consecutive
    // values in category order. Mode 1: two nested IntColumns (v / LONG_SPLIT, v % LONG_SPLIT).
    static int estimateLongPart(StringTokenPart part) {
        long[] longs = part.longs;
        int n = longs.length;
        long dMin = 0, dMax = 0;
        boolean any = false;
        for (int i = 1; i < n; i++) {
            long d = longs[i] - longs[i - 1];
            dMin = any ? Math.min(dMin, d) : d;
            dMax = any ? Math.max(dMax, d) : d;
            any = true;
        }
        long[] range = new long[]{dMin, dMax, 0};
        int bits = BitPacking.bitsFor(range);
        int mode0 = bits < 0 ? Integer.MAX_VALUE : 17 + BitPacking.sizeInBytes(n - 1, bits);

        part.hi = new int[n];
        part.lo = new int[n];
        for (int i = 0; i < n; i++) {
            part.hi[i] = (int) (longs[i] / LONG_SPLIT);
            part.lo[i] = (int) (longs[i] % LONG_SPLIT);
        }
        int mode1 = 9 + IntColumn.Encoding.choose(part.hi, n).size
                + IntColumn.Encoding.choose(part.lo, n).size;

        part.longMode = mode0 <= mode1 ? 0 : 1;
        part.longDeltaRange = range;
        return part.longMode == 0 ? mode0 : mode1;
    }

    static void encodeLongPart(BufferAccessor buf, StringTokenPart part) {
        buf.writeInt8((byte) part.longMode);
        if (part.longMode == 0) {
            long[] longs = part.longs;
            long dBase = part.longDeltaRange[0];
            buf.writeInt64(longs[0]);
            buf.writeInt64(dBase);
            int[] rel = new int[longs.length - 1];
            for (int i = 1; i < longs.length; i++)
                rel[i - 1] = (int) (longs[i] - longs[i - 1] - dBase);
            BitPacking.write(buf, rel, 0, rel.length, new long[]{0, part.longDeltaRange[1] - dBase, 0});
        }
        else {
            new IntColumn("", part.hi).encode(buf);
            new IntColumn("", part.lo).encode(buf);
        }
    }

    static void decodeLongPart(BufferAccessor buf, StringTokenPart part) {
        int mode = buf.readInt8();
        if (mode == 0) {
            long acc = buf.readInt64();
            long dBase = buf.readInt64();
            int[] rel = BitPacking.read(buf);
            long[] longs = new long[rel.length + 1];
            longs[0] = acc;
            for (int i = 0; i < rel.length; i++)
                longs[i + 1] = acc += dBase + rel[i];
            part.longs = longs;
        }
        else {
            IntColumn hi = new IntColumn("", 0);
            hi.decode(buf);
            IntColumn lo = new IntColumn("", 0);
            lo.decode(buf);
            int[] h = (int[]) hi.toArray();
            int[] l = (int[]) lo.toArray();
            long[] longs = new long[h.length];
            for (int i = 0; i < longs.length; i++)
                longs[i] = h[i] * LONG_SPLIT + l[i];
            part.longs = longs;
        }
    }

    // Writes everything after the 4-byte encoder id, up to (not including) the row->category
    // index column (written by the caller, as for every string encoder).
    public void encode(BufferAccessor buf) {
        buf.writeInt8((byte) VERSION);
        buf.writeInt32(emptyCatIdx);
        buf.writeInt32(catCount);
        buf.writeInt32(parts.size());
        for (StringTokenPart part : parts) {
            buf.writeInt8((byte) part.kind);
            if (part.kind == StringTokenPart.LITERAL) {
                buf.writeString(part.text);
                continue;
            }
            if (part.kind == StringTokenPart.INT || part.kind == StringTokenPart.LONG_INT) {
                buf.writeInt8((byte) part.flags);
                buf.writeInt8((byte) part.width);
            }
            if (part.kind == StringTokenPart.LONG_INT)
                encodeLongPart(buf, part);
            else if (part.kind == StringTokenPart.INT)
                part.intColumn.encode(buf);
            else
                part.stringColumn.encode(buf);
        }
    }

    // ---- reader side ----

    // Ports StringTokensEncoder.decode (string_column_encoders.dart) up to, but not including,
    // the trailing indices int column (read by the caller).
    public static String[] decode(BufferAccessor buf, String colName) {
        int version = buf.readInt8();
        if (version != VERSION)
            throw new RuntimeException("decoding " + colName + ": string:tokens version " + version + " not supported");
        StringTokens tokens = new StringTokens();
        tokens.emptyCatIdx = buf.readInt32();
        tokens.catCount = buf.readInt32();
        int partCount = buf.readInt32();
        for (int i = 0; i < partCount; i++) {
            StringTokenPart part = new StringTokenPart(buf.readInt8());
            if (part.kind == StringTokenPart.LITERAL)
                part.text = buf.readString();
            else if (part.kind == StringTokenPart.INT) {
                part.flags = buf.readInt8();
                part.width = buf.readInt8();
                IntColumn c = new IntColumn("", 0);
                c.decode(buf);
                part.ints = (int[]) c.toArray();
            }
            else if (part.kind == StringTokenPart.LONG_INT) {
                part.flags = buf.readInt8();
                part.width = buf.readInt8();
                decodeLongPart(buf, part);
            }
            else if (part.kind == StringTokenPart.STRING) {
                StringColumn c = new StringColumn("");
                c.decode(buf);
                part.stringColumn = c;
            }
            else
                throw new RuntimeException("decoding " + colName + ": string:tokens part kind " + part.kind + " not supported");
            tokens.parts.add(part);
        }
        return tokens.buildCategories();
    }

    // Rebuilds the category list from decoded [parts]. Mirrors the Dart digit-table writer:
    // one char buffer reused across categories, one String allocation per category.
    public String[] buildCategories() {
        String[] cats = new String[catCount];
        Builder b = new Builder();
        int j = 0;
        for (int i = 0; i < catCount; i++) {
            if (i == emptyCatIdx) {
                cats[i] = "";
                continue;
            }
            b.pos = 0;
            for (StringTokenPart part : parts) {
                if (part.kind == StringTokenPart.LITERAL)
                    b.writeText(part.text);
                else if (part.kind == StringTokenPart.INT)
                    b.writeInt(part.ints[j],
                            (part.flags & StringTokenPart.FLAG_PADDED) != 0 ? part.width : 0,
                            part.isHex() ? 16 : 10,
                            (part.flags & StringTokenPart.FLAG_HEX_LOWER) != 0 ? 87 : 55);
                else if (part.kind == StringTokenPart.LONG_INT)
                    b.writeInt(part.longs[j], part.width, 10, 0);
                else
                    b.writeText(part.stringColumn.get(j));
            }
            cats[i] = new String(b.buf, 0, b.pos);
            j++;
        }
        return cats;
    }

    private static class Builder {
        char[] buf = new char[64];
        int pos = 0;

        void ensure(int extra) {
            if (pos + extra > buf.length) {
                char[] grown = new char[Math.max(buf.length * 2, pos + extra)];
                System.arraycopy(buf, 0, grown, 0, pos);
                buf = grown;
            }
        }

        void writeText(String text) {
            ensure(text.length());
            for (int k = 0; k < text.length(); k++)
                buf[pos + k] = text.charAt(k);
            pos += text.length();
        }

        // Writes [v] right-aligned in [width] digits (0 = natural width) without intermediate strings.
        void writeInt(long v, int width, int base, int letterBase) {
            if (width == 0) {
                width = 1;
                for (long t = v; t >= base; t /= base)
                    width++;
            }
            ensure(width);
            int end = pos + width;
            long t = v;
            for (int k = end - 1; k >= pos; k--) {
                int d = (int) (t % base);
                buf[k] = (char) (d < 10 ? 48 + d : letterBase + d);
                t /= base;
            }
            pos = end;
        }
    }

    private static class TokenSegment {
        final String literal;
        final int field;

        private TokenSegment(String literal, int field) {
            this.literal = literal;
            this.field = field;
        }

        static TokenSegment sep(String literal) { return new TokenSegment(literal, -1); }
        static TokenSegment ofField(int field) { return new TokenSegment(null, field); }
        boolean isLiteral() { return literal != null; }
    }

    // Per-field statistics accumulated over all categories in one pass.
    private static class FieldStats {
        List<Integer> sig;
        List<Integer> runWidthMin = new ArrayList<>();
        List<Integer> runWidthMax = new ArrayList<>();
        List<Boolean> runLeadingZero = new ArrayList<>();
        boolean consistent = true;
        Integer widthMin;
        Integer widthMax;
        boolean hexUpper = true;
        boolean hexLower = true;

        // Scans one field starting at [start] in a single pass: classifies characters,
        // accumulates run and hex stats, and returns the field end.
        int addAndScan(String s, int start) {
            int len = s.length();
            boolean building = sig == null;
            if (building)
                sig = new ArrayList<>();
            int runIdx = 0;
            int runStart = start;
            int cls = -1;
            int end = start;
            for (; end < len; end++) {
                int cu = s.charAt(end);
                int c;
                if (cu >= 48 && cu <= 57)
                    c = _DIGIT;
                else if (cu < 128 && !(cu >= 65 && cu <= 90) && !(cu >= 97 && cu <= 122))
                    break;
                else {
                    c = _LETTER;
                    hexUpper = hexUpper && cu >= 65 && cu <= 70;
                    hexLower = hexLower && cu >= 97 && cu <= 102;
                }
                if (c != cls) {
                    if (cls >= 0)
                        endRun(building, cls, runIdx++, s, runStart, end);
                    runStart = end;
                    cls = c;
                }
            }
            if (cls >= 0)
                endRun(building, cls, runIdx++, s, runStart, end);
            if (!building && consistent && runIdx != sig.size())
                consistent = false;
            int width = end - start;
            if (widthMin == null || width < widthMin) widthMin = width;
            if (widthMax == null || width > widthMax) widthMax = width;
            return end;
        }

        void endRun(boolean building, int cls, int runIdx, String s, int runStart, int runEnd) {
            int w = runEnd - runStart;
            boolean lz = cls == _DIGIT && w > 1 && s.charAt(runStart) == 48;
            if (building) {
                sig.add(cls);
                runWidthMin.add(w);
                runWidthMax.add(w);
                runLeadingZero.add(lz);
            }
            else if (consistent) {
                if (runIdx >= sig.size() || sig.get(runIdx) != cls)
                    consistent = false;
                else {
                    if (w < runWidthMin.get(runIdx)) runWidthMin.set(runIdx, w);
                    if (w > runWidthMax.get(runIdx)) runWidthMax.set(runIdx, w);
                    if (lz) runLeadingZero.set(runIdx, true);
                }
            }
        }
    }
}
