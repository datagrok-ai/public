/*

The MIT License (MIT)

Copyright (c) 2016 Tom Zoehner
Copyright (c) 2018 Thomas Bluemel

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

export class EMFJSError extends Error {
    constructor(message: string) {
        super(message); // 'Error' breaks prototype chain here
        Object.setPrototypeOf(this, new.target.prototype); // restore prototype chain
    }
}

let isLoggingEnabled = true;

export function loggingEnabled(enabled: boolean): void {
    isLoggingEnabled = enabled;
}

export class Helper {
    public static GDI = {
        FormatSignature: {
            ENHMETA_SIGNATURE: 0x464D4520,
            EPS_SIGNATURE: 0x46535045,
        },
        BITMAPINFOHEADER_SIZE: 40,
        BITMAPCOREHEADER_SIZE: 12,
        RecordType: {
            EMR_POLYBEZIER: 0x00000002,
            EMR_POLYGON: 0x00000003,
            EMR_POLYLINE: 0x00000004,
            EMR_POLYBEZIERTO: 0x00000005,
            EMR_POLYLINETO: 0x00000006,
            EMR_POLYPOLYLINE: 0x00000007,
            EMR_POLYPOLYGON: 0x00000008,
            EMR_SETWINDOWEXTEX: 0x00000009,
            EMR_SETWINDOWORGEX: 0x0000000A,
            EMR_SETVIEWPORTEXTEX: 0x0000000B,
            EMR_SETVIEWPORTORGEX: 0x0000000C,
            EMR_SETBRUSHORGEX: 0x0000000D,
            EMR_EOF: 0x0000000E,
            EMR_SETPIXELV: 0x0000000F,
            EMR_SETMAPPERFLAGS: 0x00000010,
            EMR_SETMAPMODE: 0x00000011,
            EMR_SETBKMODE: 0x00000012,
            EMR_SETPOLYFILLMODE: 0x00000013,
            EMR_SETROP2: 0x00000014,
            EMR_SETSTRETCHBLTMODE: 0x00000015,
            EMR_SETTEXTALIGN: 0x00000016,
            EMR_SETCOLORADJUSTMENT: 0x00000017,
            EMR_SETTEXTCOLOR: 0x00000018,
            EMR_SETBKCOLOR: 0x00000019,
            EMR_OFFSETCLIPRGN: 0x0000001A,
            EMR_MOVETOEX: 0x0000001B,
            EMR_SETMETARGN: 0x0000001C,
            EMR_EXCLUDECLIPRECT: 0x0000001D,
            EMR_INTERSECTCLIPRECT: 0x0000001E,
            EMR_SCALEVIEWPORTEXTEX: 0x0000001F,
            EMR_SCALEWINDOWEXTEX: 0x00000020,
            EMR_SAVEDC: 0x00000021,
            EMR_RESTOREDC: 0x00000022,
            EMR_SETWORLDTRANSFORM: 0x00000023,
            EMR_MODIFYWORLDTRANSFORM: 0x00000024,
            EMR_SELECTOBJECT: 0x00000025,
            EMR_CREATEPEN: 0x00000026,
            EMR_CREATEBRUSHINDIRECT: 0x00000027,
            EMR_DELETEOBJECT: 0x00000028,
            EMR_ANGLEARC: 0x00000029,
            EMR_ELLIPSE: 0x0000002A,
            EMR_RECTANGLE: 0x0000002B,
            EMR_ROUNDRECT: 0x0000002C,
            EMR_ARC: 0x0000002D,
            EMR_CHORD: 0x0000002E,
            EMR_PIE: 0x0000002F,
            EMR_SELECTPALETTE: 0x00000030,
            EMR_CREATEPALETTE: 0x00000031,
            EMR_SETPALETTEENTRIES: 0x00000032,
            EMR_RESIZEPALETTE: 0x00000033,
            EMR_REALIZEPALETTE: 0x00000034,
            EMR_EXTFLOODFILL: 0x00000035,
            EMR_LINETO: 0x00000036,
            EMR_ARCTO: 0x00000037,
            EMR_POLYDRAW: 0x00000038,
            EMR_SETARCDIRECTION: 0x00000039,
            EMR_SETMITERLIMIT: 0x0000003A,
            EMR_BEGINPATH: 0x0000003B,
            EMR_ENDPATH: 0x0000003C,
            EMR_CLOSEFIGURE: 0x0000003D,
            EMR_FILLPATH: 0x0000003E,
            EMR_STROKEANDFILLPATH: 0x0000003F,
            EMR_STROKEPATH: 0x00000040,
            EMR_FLATTENPATH: 0x00000041,
            EMR_WIDENPATH: 0x00000042,
            EMR_SELECTCLIPPATH: 0x00000043,
            EMR_ABORTPATH: 0x00000044,
            EMR_COMMENT: 0x00000046,
            EMR_FILLRGN: 0x00000047,
            EMR_FRAMERGN: 0x00000048,
            EMR_INVERTRGN: 0x00000049,
            EMR_PAINTRGN: 0x0000004A,
            EMR_EXTSELECTCLIPRGN: 0x0000004B,
            EMR_BITBLT: 0x0000004C,
            EMR_STRETCHBLT: 0x0000004D,
            EMR_MASKBLT: 0x0000004E,
            EMR_PLGBLT: 0x0000004F,
            EMR_SETDIBITSTODEVICE: 0x00000050,
            EMR_STRETCHDIBITS: 0x00000051,
            EMR_EXTCREATEFONTINDIRECTW: 0x00000052,
            EMR_EXTTEXTOUTA: 0x00000053,
            EMR_EXTTEXTOUTW: 0x00000054,
            EMR_POLYBEZIER16: 0x00000055,
            EMR_POLYGON16: 0x00000056,
            EMR_POLYLINE16: 0x00000057,
            EMR_POLYBEZIERTO16: 0x00000058,
            EMR_POLYLINETO16: 0x00000059,
            EMR_POLYPOLYLINE16: 0x0000005A,
            EMR_POLYPOLYGON16: 0x0000005B,
            EMR_POLYDRAW16: 0x0000005C,
            EMR_CREATEMONOBRUSH: 0x0000005D,
            EMR_CREATEDIBPATTERNBRUSHPT: 0x0000005E,
            EMR_EXTCREATEPEN: 0x0000005F,
            EMR_POLYTEXTOUTA: 0x00000060,
            EMR_POLYTEXTOUTW: 0x00000061,
            EMR_SETICMMODE: 0x00000062,
            EMR_CREATECOLORSPACE: 0x00000063,
            EMR_SETCOLORSPACE: 0x00000064,
            EMR_DELETECOLORSPACE: 0x00000065,
            EMR_GLSRECORD: 0x00000066,
            EMR_GLSBOUNDEDRECORD: 0x00000067,
            EMR_PIXELFORMAT: 0x00000068,
            EMR_DRAWESCAPE: 0x00000069,
            EMR_EXTESCAPE: 0x0000006A,
            EMR_SMALLTEXTOUT: 0x0000006C,
            EMR_FORCEUFIMAPPING: 0x0000006D,
            EMR_NAMEDESCAPE: 0x0000006E,
            EMR_COLORCORRECTPALETTE: 0x0000006F,
            EMR_SETICMPROFILEA: 0x00000070,
            EMR_SETICMPROFILEW: 0x00000071,
            EMR_ALPHABLEND: 0x00000072,
            EMR_SETLAYOUT: 0x00000073,
            EMR_TRANSPARENTBLT: 0x00000074,
            EMR_GRADIENTFILL: 0x00000076,
            EMR_SETLINKEDUFIS: 0x00000077,
            EMR_SETTEXTJUSTIFICATION: 0x00000078,
            EMR_COLORMATCHTOTARGETW: 0x00000079,
            EMR_CREATECOLORSPACEW: 0x0000007A,
        },
        MetafileEscapes: {
            NEWFRAME: 0x0001,
            ABORTDOC: 0x0002,
            NEXTBAND: 0x0003,
            SETCOLORTABLE: 0x0004,
            GETCOLORTABLE: 0x0005,
            FLUSHOUT: 0x0006,
            DRAFTMODE: 0x0007,
            QUERYESCSUPPORT: 0x0008,
            SETABORTPROC: 0x0009,
            STARTDOC: 0x000a,
            ENDDOC: 0x000b,
            GETPHYSPAGESIZE: 0x000c,
            GETPRINTINGOFFSET: 0x000d,
            GETSCALINGFACTOR: 0x000e,
            META_ESCAPE_ENHANCED_METAFILE: 0x000f,
            SETPENWIDTH: 0x0010,
            SETCOPYCOUNT: 0x0011,
            SETPAPERSOURCE: 0x0012,
            PASSTHROUGH: 0x0013,
            GETTECHNOLOGY: 0x0014,
            SETLINECAP: 0x0015,
            SETLINEJOIN: 0x0016,
            SETMITERLIMIT: 0x0017,
            BANDINFO: 0x0018,
            DRAWPATTERNRECT: 0x0019,
            GETVECTORPENSIZE: 0x001a,
            GETVECTORBRUSHSIZE: 0x001b,
            ENABLEDUPLEX: 0x001c,
            GETSETPAPERBINS: 0x001d,
            GETSETPRINTORIENT: 0x001e,
            ENUMPAPERBINS: 0x001f,
            SETDIBSCALING: 0x0020,
            EPSPRINTING: 0x0021,
            ENUMPAPERMETRICS: 0x0022,
            GETSETPAPERMETRICS: 0x0023,
            POSTSCRIPT_DATA: 0x0025,
            POSTSCRIPT_IGNORE: 0x0026,
            GETDEVICEUNITS: 0x002a,
            GETEXTENDEDTEXTMETRICS: 0x0100,
            GETPAIRKERNTABLE: 0x0102,
            EXTTEXTOUT: 0x0200,
            GETFACENAME: 0x0201,
            DOWNLOADFACE: 0x0202,
            METAFILE_DRIVER: 0x0801,
            QUERYDIBSUPPORT: 0x0c01,
            BEGIN_PATH: 0x1000,
            CLIP_TO_PATH: 0x1001,
            END_PATH: 0x1002,
            OPEN_CHANNEL: 0x100e,
            DOWNLOADHEADER: 0x100f,
            CLOSE_CHANNEL: 0x1010,
            POSTSCRIPT_PASSTHROUGH: 0x1013,
            ENCAPSULATED_POSTSCRIPT: 0x1014,
            POSTSCRIPT_IDENTIFY: 0x1015,
            POSTSCRIPT_INJECTION: 0x1016,
            CHECKJPEGFORMAT: 0x1017,
            CHECKPNGFORMAT: 0x1018,
            GET_PS_FEATURESETTING: 0x1019,
            MXDC_ESCAPE: 0x101a,
            SPCLPASSTHROUGH2: 0x11d8,
        },
        MapMode: {
            MM_TEXT: 1,
            MM_LOMETRIC: 2,
            MM_HIMETRIC: 3,
            MM_LOENGLISH: 4,
            MM_HIENGLISH: 5,
            MM_TWIPS: 6,
            MM_ISOTROPIC: 7,
            MM_ANISOTROPIC: 8,
        },
        StretchMode: {
            BLACKONWHITE: 1,
            WHITEONBLACK: 2,
            COLORONCOLOR: 3,
            HALFTONE: 4,
        },
        MixMode: {
            TRANSPARENT: 1,
            OPAQUE: 2,
        },
        BrushStyle: {
            BS_SOLID: 0,
            BS_NULL: 1,
            BS_HATCHED: 2,
            BS_PATTERN: 3,
            BS_INDEXED: 4,
            BS_DIBPATTERN: 5,
            BS_DIBPATTERNPT: 6,
            BS_PATTERN8X8: 7,
            BS_DIBPATTERN8X8: 8,
            BS_MONOPATTERN: 9,
        },
        PenStyle: {
            PS_COSMETIC: 0x00000000,
            PS_ENDCAP_ROUND: 0x00000000,
            PS_JOIN_ROUND: 0x00000000,
            PS_SOLID: 0x00000000,
            PS_DASH: 0x00000001,
            PS_DOT: 0x00000002,
            PS_DASHDOT: 0x00000003,
            PS_DASHDOTDOT: 0x00000004,
            PS_NULL: 0x00000005,
            PS_INSIDEFRAME: 0x00000006,
            PS_USERSTYLE: 0x00000007,
            PS_ALTERNATE: 0x00000008,
            PS_ENDCAP_SQUARE: 0x00000100,
            PS_ENDCAP_FLAT: 0x00000200,
            PS_JOIN_BEVEL: 0x00001000,
            PS_JOIN_MITER: 0x00002000,
            PS_GEOMETRIC: 0x00010000,
        },
        PolygonFillMode: {
            ALTERNATE: 1,
            WINDING: 2,
        },
        BitmapCompression: {
            BI_RGB: 0,
            BI_RLE8: 1,
            BI_RLE4: 2,
            BI_BITFIELDS: 3,
            BI_JPEG: 4,
            BI_PNG: 5,
        },
        RegionMode: {
            RGN_AND: 1,
            RGN_OR: 2,
            RGN_XOR: 3,
            RGN_DIFF: 4,
            RGN_COPY: 5,
        },
        StockObject: {
            WHITE_BRUSH: 0x80000000,
            LTGRAY_BRUSH: 0x80000001,
            GRAY_BRUSH: 0x80000002,
            DKGRAY_BRUSH: 0x80000003,
            BLACK_BRUSH: 0x80000004,
            NULL_BRUSH: 0x80000005,
            WHITE_PEN: 0x80000006,
            BLACK_PEN: 0x80000007,
            NULL_PEN: 0x80000008,
            OEM_FIXED_FONT: 0x8000000A,
            ANSI_FIXED_FONT: 0x8000000B,
            ANSI_VAR_FONT: 0x8000000C,
            SYSTEM_FONT: 0x8000000D,
            DEVICE_DEFAULT_FONT: 0x8000000E,
            DEFAULT_PALETTE: 0x8000000F,
            SYSTEM_FIXED_FONT: 0x80000010,
            DEFAULT_GUI_FONT: 0x80000011,
            DC_BRUSH: 0x80000012,
            DC_PEN: 0x80000013,
        },
    };

    public static _uniqueId = 0;

    public static log(message: string): void {
        if (isLoggingEnabled) {
            console.log(message);
        }
    }

    public static _makeUniqueId(prefix: string): string {
        return "EMFJS_" + prefix + (this._uniqueId++);
    }

    public static _writeUint32Val(uint8arr: Uint8Array, pos: number, val: number): void {
        uint8arr[pos++] = val & 0xff;
        uint8arr[pos++] = (val >>> 8) & 0xff;
        uint8arr[pos++] = (val >>> 16) & 0xff;
        uint8arr[pos++] = (val >>> 24) & 0xff;
    }

    public static _blobToBinary(blob: Uint8Array): string {
        let ret = "";
        const len = blob.length;
        for (let i = 0; i < len; i++) {
            ret += String.fromCharCode(blob[i]);
        }
        return ret;
    }
}
