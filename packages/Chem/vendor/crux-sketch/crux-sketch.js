import{$ as e,$n as t,$r as n,$t as r,A as i,Ai as a,An as o,At as s,B as c,Bn as l,Br as u,Bt as d,C as f,Ci as p,Cn as m,Cr as h,Ct as g,D as _,Dn as v,Dr as y,Dt as b,E as ee,Ei as te,En as ne,Er as re,Et as ie,F as ae,Fi as x,Fn as oe,Fr as se,Ft as ce,G as le,Gn as ue,Gr as de,Gt as fe,H as pe,Hn as S,Hr as me,Ht as he,I as ge,Ii as _e,In as ve,Ir as ye,It as be,J as xe,Jn as Se,Jr as Ce,Jt as we,K as Te,Kn as Ee,Kr as De,Kt as Oe,L as ke,Li as Ae,Ln as C,Lr as je,Lt as Me,M as Ne,Mi as w,Mn as Pe,Mr as Fe,Mt as Ie,N as Le,Ni as T,Nn as Re,Nr as ze,Nt as E,O as Be,Oi as Ve,On as He,Or as Ue,Ot as We,P as Ge,Pi as Ke,Pn as qe,Pr as Je,Pt as Ye,Q as Xe,Qn as Ze,Qr as Qe,Qt as $e,R as et,Ri as tt,Rn as nt,Rr as rt,Rt as it,S as at,Si as ot,Sn as st,Sr as ct,St as lt,T as ut,Ti as dt,Tn as ft,Tr as pt,Tt as mt,U as ht,Un as gt,Ur as _t,Ut as vt,V as yt,Vn as bt,Vr as xt,Vt as St,W as Ct,Wn as wt,Wr as Tt,Wt as Et,X as Dt,Xn as Ot,Xr as kt,Xt as At,Y as jt,Yn as Mt,Yr as Nt,Yt as Pt,Z as Ft,Zn as It,Zr as Lt,Zt as Rt,_ as zt,_i as Bt,_n as Vt,_r as Ht,_t as Ut,a as Wt,ai as Gt,an as Kt,ar as qt,at as Jt,b as Yt,bi as Xt,bn as Zt,br as Qt,bt as $t,c as en,ci as tn,cr as nn,ct as rn,d as an,di as D,dn as on,dr as sn,dt as cn,ei as ln,en as un,er as dn,et as fn,f as pn,fi as mn,fn as O,fr as hn,ft as gn,g as _n,gi as vn,gn as k,gr as yn,gt as bn,h as xn,hi as Sn,hn as Cn,hr as wn,ht as Tn,i as En,ii as Dn,in as On,ir as kn,it as An,j as jn,ji as Mn,jn as Nn,jr as Pn,jt as Fn,k as In,ki as Ln,kn as Rn,kr as zn,kt as Bn,l as Vn,li as Hn,ln as Un,lr as Wn,lt as Gn,m as Kn,mi as qn,mn as Jn,mr as Yn,mt as Xn,n as Zn,ni as Qn,nn as $n,nr as er,nt as tr,o as nr,oi as rr,on as ir,or as ar,ot as or,p as sr,pi as cr,pn as A,pr as lr,pt as ur,q as dr,qn as fr,qr as pr,qt as mr,r as hr,ri as gr,rn as _r,rt as vr,s as yr,si as br,sn as xr,sr as Sr,st as Cr,t as wr,ti as Tr,tn as Er,tr as Dr,tt as Or,u as kr,ui as Ar,un as jr,ur as Mr,ut as Nr,v as Pr,vi as Fr,vn as Ir,vr as Lr,vt as Rr,w as zr,wi as Br,wn as Vr,wr as Hr,wt as Ur,x as Wr,xi as Gr,xn as Kr,xr as qr,xt as Jr,y as Yr,yi as Xr,yn as Zr,yr as Qr,yt as $r,z as ei,zn as ti,zr as ni,zt as ri}from"./shared-BKc5K24C.js";var ii=class e{static __wrap(t){let n=Object.create(e.prototype);return n.__wbg_ptr=t,di.register(n,n.__wbg_ptr,n),n}__destroy_into_raw(){let e=this.__wbg_ptr;return this.__wbg_ptr=0,di.unregister(this),e}free(){let e=this.__destroy_into_raw();L.__wbg_sketchmol_free(e,0)}static abbreviationsJson(){let e,t;try{let n=L.sketchmol_abbreviationsJson();return e=n[0],t=n[1],j(n[0],n[1])}finally{L.__wbindgen_free(e,t,1)}}addAtom(e,t,n){let r=P(e,L.__wbindgen_malloc,L.__wbindgen_realloc),i=I,a=L.sketchmol_addAtom(this.__wbg_ptr,r,i,t,n);if(a[2])throw F(a[1]);return a[0]>>>0}addBond(e,t,n,r){let i=P(n,L.__wbindgen_malloc,L.__wbindgen_realloc),a=I,o=P(r,L.__wbindgen_malloc,L.__wbindgen_realloc),s=I,c=L.sketchmol_addBond(this.__wbg_ptr,e,t,i,a,o,s);if(c[2])throw F(c[1]);return c[0]>>>0}addHs(e,t){var n=M(e)?0:N(e,L.__wbindgen_malloc),r=I;let i=L.sketchmol_addHs(this.__wbg_ptr,n,r,!M(t),M(t)?0:t);if(i[3])throw F(i[2]);var a=mi(i[0],i[1]).slice();return L.__wbindgen_free(i[0],i[1]*4,4),a}addSGroup(e,t,n,r){let i=P(e,L.__wbindgen_malloc,L.__wbindgen_realloc),a=I,o=P(r,L.__wbindgen_malloc,L.__wbindgen_realloc),s=I,c=L.sketchmol_addSGroup(this.__wbg_ptr,i,a,t,n,o,s);if(c[2])throw F(c[1]);return c[0]>>>0}aromatize(e){let t,n;try{var r=M(e)?0:N(e,L.__wbindgen_malloc),i=I;let s=L.sketchmol_aromatize(this.__wbg_ptr,r,i);var a=s[0],o=s[1];if(s[3])throw a=0,o=0,F(s[2]);return t=a,n=o,j(a,o)}finally{L.__wbindgen_free(t,n,1)}}atomAttachmentPoints(e){let t=L.sketchmol_atomAttachmentPoints(this.__wbg_ptr,e);if(t[2])throw F(t[1]);return t[0]>>>0}atomJson(e){let t,n;try{let a=L.sketchmol_atomJson(this.__wbg_ptr,e);var r=a[0],i=a[1];if(a[3])throw r=0,i=0,F(a[2]);return t=r,n=i,j(r,i)}finally{L.__wbindgen_free(t,n,1)}}atomQueryJson(e){let t,n;try{let a=L.sketchmol_atomQueryJson(this.__wbg_ptr,e);var r=a[0],i=a[1];if(a[3])throw r=0,i=0,F(a[2]);return t=r,n=i,j(r,i)}finally{L.__wbindgen_free(t,n,1)}}atomRNumbers(e){let t=L.sketchmol_atomRNumbers(this.__wbg_ptr,e);if(t[2])throw F(t[1]);return F(t[0])}bondQueryJson(e){let t,n;try{let a=L.sketchmol_bondQueryJson(this.__wbg_ptr,e);var r=a[0],i=a[1];if(a[3])throw r=0,i=0,F(a[2]);return t=r,n=i,j(r,i)}finally{L.__wbindgen_free(t,n,1)}}chiralFlag(){return L.sketchmol_chiralFlag(this.__wbg_ptr)!==0}clearAtomQuery(e){let t=L.sketchmol_clearAtomQuery(this.__wbg_ptr,e);if(t[2])throw F(t[1]);return t[0]!==0}clearBondQuery(e){let t=L.sketchmol_clearBondQuery(this.__wbg_ptr,e);if(t[2])throw F(t[1]);return t[0]!==0}clear(){L.sketchmol_clear(this.__wbg_ptr)}clone(){let t=L.sketchmol_clone(this.__wbg_ptr);return e.__wrap(t)}contractAbbreviation(e){let t=L.sketchmol_contractAbbreviation(this.__wbg_ptr,e);if(t[2])throw F(t[1]);return t[0]!==0}contractToAbbreviation(e){let t=L.sketchmol_contractToAbbreviation(this.__wbg_ptr,e);if(t[2])throw F(t[1]);return F(t[0])}dearomatize(e){let t,n;try{var r=M(e)?0:N(e,L.__wbindgen_malloc),i=I;let s=L.sketchmol_dearomatize(this.__wbg_ptr,r,i);var a=s[0],o=s[1];if(s[3])throw a=0,o=0,F(s[2]);return t=a,n=o,j(a,o)}finally{L.__wbindgen_free(t,n,1)}}descriptors(e){let t,n;try{var r=M(e)?0:N(e,L.__wbindgen_malloc),i=I;let s=L.sketchmol_descriptors(this.__wbg_ptr,r,i);var a=s[0],o=s[1];if(s[3])throw a=0,o=0,F(s[2]);return t=a,n=o,j(a,o)}finally{L.__wbindgen_free(t,n,1)}}doubleBondPosition(e){let t,n;try{let a=L.sketchmol_doubleBondPosition(this.__wbg_ptr,e);var r=a[0],i=a[1];if(a[3])throw r=0,i=0,F(a[2]);return t=r,n=i,j(r,i)}finally{L.__wbindgen_free(t,n,1)}}drawnBondOrders(e){let t=P(e,L.__wbindgen_malloc,L.__wbindgen_realloc),n=I,r=L.sketchmol_drawnBondOrders(this.__wbg_ptr,t,n);if(r[2])throw F(r[1]);return F(r[0])}static empty(){let t=L.sketchmol_empty();return e.__wrap(t)}expandAbbreviation(e){let t=L.sketchmol_expandAbbreviation(this.__wbg_ptr,e);if(t[2])throw F(t[1]);return t[0]!==0}extract(t){let n=L.sketchmol_extract(this.__wbg_ptr,t);if(n[2])throw F(n[1]);return e.__wrap(n[0])}geometryJson(){let e,t;try{let n=L.sketchmol_geometryJson(this.__wbg_ptr);return e=n[0],t=n[1],j(n[0],n[1])}finally{L.__wbindgen_free(e,t,1)}}hasQueryIn(e,t){let n=N(e,L.__wbindgen_malloc),r=I,i=N(t,L.__wbindgen_malloc),a=I,o=L.sketchmol_hasQueryIn(this.__wbg_ptr,n,r,i,a);if(o[2])throw F(o[1]);return o[0]!==0}hasQuery(){return L.sketchmol_hasQuery(this.__wbg_ptr)!==0}insertAbbreviation(e,t,n,r){let i=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),a=I,o=L.sketchmol_insertAbbreviation(this.__wbg_ptr,M(e)?2**53-1:e>>>0,i,a,n,r);if(o[2])throw F(o[1]);return o[0]>>>0}insert(t){return pi(t,e),L.sketchmol_insert(this.__wbg_ptr,t.__wbg_ptr)>>>0}invertChirality(e){let t=L.sketchmol_invertChirality(this.__wbg_ptr,e);if(t[1])throw F(t[0])}isEdited(){return L.sketchmol_isEdited(this.__wbg_ptr)!==0}isEmpty(){return L.sketchmol_isEmpty(this.__wbg_ptr)!==0}layoutLike(t,n){let r,i;try{pi(t,e);let s=P(n,L.__wbindgen_malloc,L.__wbindgen_realloc),c=I,l=L.sketchmol_layoutLike(this.__wbg_ptr,t.__wbg_ptr,s,c);var a=l[0],o=l[1];if(l[3])throw a=0,o=0,F(l[2]);return r=a,i=o,j(a,o)}finally{L.__wbindgen_free(r,i,1)}}layout(e,t,n){var r=M(e)?0:N(e,L.__wbindgen_malloc),i=I,a=M(n)?0:P(n,L.__wbindgen_malloc,L.__wbindgen_realloc),o=I;let s=L.sketchmol_layout(this.__wbg_ptr,r,i,!M(t),M(t)?0:t,a,o);if(s[1])throw F(s[0])}massInfo(e){let t,n;try{let a=L.sketchmol_massInfo(this.__wbg_ptr,M(e)?0:fi(e));var r=a[0],i=a[1];if(a[3])throw r=0,i=0,F(a[2]);return t=r,n=i,j(r,i)}finally{L.__wbindgen_free(t,n,1)}}meanBondLength(){return L.sketchmol_meanBondLength(this.__wbg_ptr)}needsV3000(){return L.sketchmol_needsV3000(this.__wbg_ptr)!==0}normalizeDepiction(e,t){let n=L.sketchmol_normalizeDepiction(this.__wbg_ptr,M(e)?2**53-1:e>>0,!M(t),M(t)?0:t);if(n[2])throw F(n[1]);return n[0]}numAtoms(){return L.sketchmol_numAtoms(this.__wbg_ptr)>>>0}numBonds(){return L.sketchmol_numBonds(this.__wbg_ptr)>>>0}static parse(t,n,r,i,a){let o=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),s=I;var c=M(n)?0:P(n,L.__wbindgen_malloc,L.__wbindgen_realloc),l=I,u=M(i)?0:P(i,L.__wbindgen_malloc,L.__wbindgen_realloc),d=I;let f=L.sketchmol_parse(o,s,c,l,M(r)?16777215:+!!r,u,d,M(a)?16777215:+!!a);if(f[2])throw F(f[1]);return e.__wrap(f[0])}potentialStereoCentres(){let e=L.sketchmol_potentialStereoCentres(this.__wbg_ptr);var t=mi(e[0],e[1]).slice();return L.__wbindgen_free(e[0],e[1]*4,4),t}problems(){return L.sketchmol_problems(this.__wbg_ptr)}readWarnings(){return L.sketchmol_readWarnings(this.__wbg_ptr)}recomputeWedges(e){var t=M(e)?0:N(e,L.__wbindgen_malloc),n=I;let r=L.sketchmol_recomputeWedges(this.__wbg_ptr,t,n);if(r[2])throw F(r[1]);return r[0]!==0}removeAtom(e){let t=L.sketchmol_removeAtom(this.__wbg_ptr,e);if(t[1])throw F(t[0])}removeAtoms(e){let t=N(e,L.__wbindgen_malloc),n=I,r=L.sketchmol_removeAtoms(this.__wbg_ptr,t,n);if(r[1])throw F(r[0])}removeBond(e){let t=L.sketchmol_removeBond(this.__wbg_ptr,e);if(t[1])throw F(t[0])}removeHs(e){var t=M(e)?0:N(e,L.__wbindgen_malloc),n=I;let r=L.sketchmol_removeHs(this.__wbg_ptr,t,n);if(r[3])throw F(r[2]);var i=mi(r[0],r[1]).slice();return L.__wbindgen_free(r[0],r[1]*4,4),i}removeSGroup(e){let t=L.sketchmol_removeSGroup(this.__wbg_ptr,e);if(t[1])throw F(t[0])}renderSvg(e){let t,n;try{let a=P(e,L.__wbindgen_malloc,L.__wbindgen_realloc),o=I,s=L.sketchmol_renderSvg(this.__wbg_ptr,a,o);var r=s[0],i=s[1];if(s[3])throw r=0,i=0,F(s[2]);return t=r,n=i,j(r,i)}finally{L.__wbindgen_free(t,n,1)}}replaceSGroup(e,t){let n=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I,i=L.sketchmol_replaceSGroup(this.__wbg_ptr,e,n,r);if(i[2])throw F(i[1]);return F(i[0])}reverseBond(e){let t=L.sketchmol_reverseBond(this.__wbg_ptr,e);if(t[1])throw F(t[0])}ringBonds(){let e=L.sketchmol_ringBonds(this.__wbg_ptr);if(e[2])throw F(e[1]);return F(e[0])}setAtomAlias(e,t){var n=M(t)?0:P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I;let i=L.sketchmol_setAtomAlias(this.__wbg_ptr,e,n,r);if(i[1])throw F(i[0])}setAtomAttachmentPoints(e,t){let n=L.sketchmol_setAtomAttachmentPoints(this.__wbg_ptr,e,t);if(n[2])throw F(n[1]);return n[0]!==0}setAtomCharge(e,t){let n=L.sketchmol_setAtomCharge(this.__wbg_ptr,e,t);if(n[1])throw F(n[0])}setAtomElement(e,t){let n=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I,i=L.sketchmol_setAtomElement(this.__wbg_ptr,e,n,r);if(i[2])throw F(i[1]);return i[0]!==0}setAtomElements(e,t){let n=N(e,L.__wbindgen_malloc),r=I,i=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),a=I,o=L.sketchmol_setAtomElements(this.__wbg_ptr,n,r,i,a);if(o[2])throw F(o[1]);return o[0]!==0}setAtomIsotope(e,t){let n=L.sketchmol_setAtomIsotope(this.__wbg_ptr,e,t);if(n[1])throw F(n[0])}setAtomLabel(e,t){let n=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I,i=L.sketchmol_setAtomLabel(this.__wbg_ptr,e,n,r);if(i[2])throw F(i[1]);return F(i[0])}setAtomPosition(e,t,n){let r=L.sketchmol_setAtomPosition(this.__wbg_ptr,e,t,n);if(r[1])throw F(r[0])}setAtomQueries(e,t){let n=N(e,L.__wbindgen_malloc),r=I,i=L.sketchmol_setAtomQueries(this.__wbg_ptr,n,r,t);if(i[2])throw F(i[1]);return i[0]!==0}setAtomQuerySmarts(e,t){let n=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I,i=L.sketchmol_setAtomQuerySmarts(this.__wbg_ptr,e,n,r);if(i[2])throw F(i[1]);return F(i[0])}setAtomQuery(e,t){let n=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I,i=L.sketchmol_setAtomQuery(this.__wbg_ptr,e,n,r);if(i[2])throw F(i[1]);return i[0]!==0}setAtomRNumbers(e,t){let n=L.sketchmol_setAtomRNumbers(this.__wbg_ptr,e,t);if(n[2])throw F(n[1]);return n[0]!==0}setAtomRadical(e,t){let n=L.sketchmol_setAtomRadical(this.__wbg_ptr,e,t);if(n[1])throw F(n[0])}setAtomValence(e,t){let n=L.sketchmol_setAtomValence(this.__wbg_ptr,e,!M(t),M(t)?0:t);if(n[1])throw F(n[0])}setAtropisomer(e,t){var n=M(t)?0:P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I;let i=L.sketchmol_setAtropisomer(this.__wbg_ptr,e,n,r);if(i[2])throw F(i[1]);return i[0]!==0}setBondQueries(e,t){let n=N(e,L.__wbindgen_malloc),r=I,i=L.sketchmol_setBondQueries(this.__wbg_ptr,n,r,t);if(i[2])throw F(i[1]);return i[0]!==0}setBondQuerySmarts(e,t){let n=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I,i=L.sketchmol_setBondQuerySmarts(this.__wbg_ptr,e,n,r);if(i[2])throw F(i[1]);return F(i[0])}setBondQuery(e,t){let n=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I,i=L.sketchmol_setBondQuery(this.__wbg_ptr,e,n,r);if(i[2])throw F(i[1]);return i[0]!==0}setBondStereo(e,t,n){let r=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),i=I,a=L.sketchmol_setBondStereo(this.__wbg_ptr,e,r,i,M(n)?2**53-1:n>>>0);if(a[1])throw F(a[0])}setBondType(e,t){let n=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I,i=L.sketchmol_setBondType(this.__wbg_ptr,e,n,r);if(i[1])throw F(i[0])}setBondTypes(e,t){let n=N(e,L.__wbindgen_malloc),r=I,i=L.sketchmol_setBondTypes(this.__wbg_ptr,n,r,t);if(i[1])throw F(i[0])}setChiralFlag(e){L.sketchmol_setChiralFlag(this.__wbg_ptr,e)}setDoubleBondPosition(e,t){let n=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I,i=L.sketchmol_setDoubleBondPosition(this.__wbg_ptr,e,n,r);if(i[1])throw F(i[0])}setSGroupProps(e,t){let n=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I,i=L.sketchmol_setSGroupProps(this.__wbg_ptr,e,n,r);if(i[1])throw F(i[0])}setStereoCentre(e,t){var n=M(t)?0:P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I;let i=L.sketchmol_setStereoCentre(this.__wbg_ptr,e,n,r);if(i[2])throw F(i[1]);return i[0]!==0}setStereoGroup(e,t,n,r){let i=P(t,L.__wbindgen_malloc,L.__wbindgen_realloc),a=I,o=L.sketchmol_setStereoGroup(this.__wbg_ptr,e,i,a,!M(n),M(n)?0:n,M(r)?0:fi(r));if(o[2])throw F(o[1]);return o[0]>>>0}sgroupsJson(){let e,t;try{let n=L.sketchmol_sgroupsJson(this.__wbg_ptr);return e=n[0],t=n[1],j(n[0],n[1])}finally{L.__wbindgen_free(e,t,1)}}stereoGroups(){let e,t;try{let n=L.sketchmol_stereoGroups(this.__wbg_ptr);return e=n[0],t=n[1],j(n[0],n[1])}finally{L.__wbindgen_free(e,t,1)}}stereoInfo(e){let t,n;try{let r=L.sketchmol_stereoInfo(this.__wbg_ptr,M(e)?16777215:+!!e);return t=r[0],n=r[1],j(r[0],r[1])}finally{L.__wbindgen_free(t,n,1)}}stereoisomerCount(){return L.sketchmol_stereoisomerCount(this.__wbg_ptr)}straightenDepiction(e){L.sketchmol_straightenDepiction(this.__wbg_ptr,M(e)?16777215:+!!e)}straightenTurn(e,t){let n=N(e,L.__wbindgen_malloc),r=I,i=L.sketchmol_straightenTurn(this.__wbg_ptr,n,r,t);if(i[2])throw F(i[1]);return i[0]}toCxsmilesCanonical(){let e,t;try{let n=L.sketchmol_toCxsmilesCanonical(this.__wbg_ptr);return e=n[0],t=n[1],j(n[0],n[1])}finally{L.__wbindgen_free(e,t,1)}}toCxsmiles(e){let t,n;try{let r=L.sketchmol_toCxsmiles(this.__wbg_ptr,M(e)?16777215:+!!e);return t=r[0],n=r[1],j(r[0],r[1])}finally{L.__wbindgen_free(t,n,1)}}toMolblockV3000(){let e,t;try{let i=L.sketchmol_toMolblockV3000(this.__wbg_ptr);var n=i[0],r=i[1];if(i[3])throw n=0,r=0,F(i[2]);return e=n,t=r,j(n,r)}finally{L.__wbindgen_free(e,t,1)}}toMolblock(){let e,t;try{let i=L.sketchmol_toMolblock(this.__wbg_ptr);var n=i[0],r=i[1];if(i[3])throw n=0,r=0,F(i[2]);return e=n,t=r,j(n,r)}finally{L.__wbindgen_free(e,t,1)}}toSmarts(){let e,t;try{let i=L.sketchmol_toSmarts(this.__wbg_ptr);var n=i[0],r=i[1];if(i[3])throw n=0,r=0,F(i[2]);return e=n,t=r,j(n,r)}finally{L.__wbindgen_free(e,t,1)}}toSmilesCanonical(){let e,t;try{let n=L.sketchmol_toSmilesCanonical(this.__wbg_ptr);return e=n[0],t=n[1],j(n[0],n[1])}finally{L.__wbindgen_free(e,t,1)}}toSmiles(){let e,t;try{let n=L.sketchmol_toSmiles(this.__wbg_ptr);return e=n[0],t=n[1],j(n[0],n[1])}finally{L.__wbindgen_free(e,t,1)}}writeWarnings(e,t){let n=P(e,L.__wbindgen_malloc,L.__wbindgen_realloc),r=I,i=L.sketchmol_writeWarnings(this.__wbg_ptr,n,r,M(t)?16777215:+!!t);if(i[2])throw F(i[1]);return F(i[0])}};Symbol.dispose&&(ii.prototype[Symbol.dispose]=ii.prototype.free);function ai(e){let t,n;try{let r=P(e,L.__wbindgen_malloc,L.__wbindgen_realloc),i=I,a=L.detectFormat(r,i);return t=a[0],n=a[1],j(a[0],a[1])}finally{L.__wbindgen_free(t,n,1)}}function oi(){let e,t;try{let n=L.engineRev();return e=n[0],t=n[1],j(n[0],n[1])}finally{L.__wbindgen_free(e,t,1)}}function si(){return L.liveMolCount()>>>0}function ci(e){return L.useLegacyStereoPerception(e)!==0}function li(){let e,t;try{let n=L.version();return e=n[0],t=n[1],j(n[0],n[1])}finally{L.__wbindgen_free(e,t,1)}}function ui(){return{__proto__:null,"./sketch_wasm_bg.js":{__proto__:null,__wbg___wbindgen_number_get_2e0e7dee9f701a71:function(e,t){let n=t,r=typeof n==`number`?n:void 0;gi().setFloat64(e+8,M(r)?0:r,!0),gi().setInt32(e+0,!M(r),!0)},__wbg___wbindgen_string_get_0380ccaa2f57f0d9:function(e,t){let n=t,r=typeof n==`string`?n:void 0;var i=M(r)?0:P(r,L.__wbindgen_malloc,L.__wbindgen_realloc),a=I;gi().setInt32(e+4,a,!0),gi().setInt32(e+0,i,!0)},__wbg___wbindgen_throw_41e9ee4f547fc59a:function(e,t){throw Error(j(e,t))},__wbg_error_757e9472f8410341:function(e,t){let n,r;try{n=e,r=t,console.error(j(e,t))}finally{L.__wbindgen_free(n,r,1)}},__wbg_get_unchecked_288889d017702237:function(e,t){return e[t>>>0]},__wbg_length_d4bdea10311bd9cf:function(e){return e.length},__wbg_new_227d7c05414eb861:function(){return Error()},__wbg_new_343a093a3c2ffb4e:function(e,t){return Error(j(e,t))},__wbg_new_617a8cdb8bb1130e:function(){return{}},__wbg_new_ee2291f50781bf1d:function(){return[]},__wbg_new_with_length_469fcc27bd71672e:function(e){return Array(e>>>0)},__wbg_push_2baf45db356cf468:function(e,t){return e.push(t)},__wbg_set_145a351398b48c65:function(){return xi(function(e,t,n){return Reflect.set(e,t,n)},arguments)},__wbg_set_bea140a88be9b277:function(e,t,n){e[t>>>0]=n},__wbg_stack_3b0d974bbf31e44f:function(e,t){let n=t.stack,r=P(n,L.__wbindgen_malloc,L.__wbindgen_realloc),i=I;gi().setInt32(e+4,i,!0),gi().setInt32(e+0,r,!0)},__wbindgen_generic_0000000000000001:function(e){return e},__wbindgen_generic_0000000000000002:function(e,t){return j(e,t)},__wbindgen_init_externref_table:function(){let e=L.__wbindgen_externrefs,t=e.grow(4);e.set(0,void 0),e.set(t+0,void 0),e.set(t+1,null),e.set(t+2,!0),e.set(t+3,!1)}}}}var di=typeof FinalizationRegistry>`u`?{register:()=>{},unregister:()=>{}}:new FinalizationRegistry(e=>L.__wbg_sketchmol_free(e,1));function fi(e){let t=L.__externref_table_alloc();return L.__wbindgen_externrefs.set(t,e),t}function pi(e,t){if(!(e instanceof t))throw Error(`expected instance of ${t.name}`)}function mi(e,t){return e>>>=0,vi().subarray(e/4,e/4+t)}var hi=null;function gi(){return(hi===null||hi.buffer.detached===!0||hi.buffer.detached===void 0&&hi.buffer!==L.memory.buffer)&&(hi=new DataView(L.memory.buffer)),hi}function j(e,t){return Ti(e>>>0,t)}var _i=null;function vi(){return(_i===null||_i.byteLength===0)&&(_i=new Uint32Array(L.memory.buffer)),_i}var yi=null;function bi(){return(yi===null||yi.byteLength===0)&&(yi=new Uint8Array(L.memory.buffer)),yi}function xi(e,t){try{return e.apply(this,t)}catch(e){let t=fi(e);L.__wbindgen_exn_store(t)}}function M(e){return e==null}function N(e,t){let n=t(e.length*4,4)>>>0;return vi().set(e,n/4),I=e.length,n}function P(e,t,n){if(n===void 0){let n=Ei.encode(e),r=t(n.length,1)>>>0;return bi().subarray(r,r+n.length).set(n),I=n.length,r}let r=e.length,i=t(r,1)>>>0,a=bi(),o=0;for(;o<r;o++){let t=e.charCodeAt(o);if(t>127)break;a[i+o]=t}if(o!==r){o!==0&&(e=e.slice(o)),i=n(i,r,r=o+e.length*3,1)>>>0;let t=bi().subarray(i+o,i+r),a=Ei.encodeInto(e,t);o+=a.written,i=n(i,r,o,1)>>>0}return I=o,i}function F(e){let t=L.__wbindgen_externrefs.get(e);return L.__externref_table_dealloc(e),t}var Si=new TextDecoder(`utf-8`,{ignoreBOM:!0,fatal:!0});Si.decode();var Ci=2146435072,wi=0;function Ti(e,t){return wi+=t,wi>=Ci&&(Si=new TextDecoder(`utf-8`,{ignoreBOM:!0,fatal:!0}),Si.decode(),wi=t),Si.decode(bi().subarray(e,e+t))}var Ei=new TextEncoder;`encodeInto`in Ei||(Ei.encodeInto=function(e,t){let n=Ei.encode(e);return t.set(n),{read:e.length,written:n.length}});var I=0,L;function Di(e,t){return L=e.exports,hi=null,_i=null,yi=null,L.__wbindgen_start(),L}async function Oi(e,t){if(typeof Response==`function`&&e instanceof Response){if(!e.ok)throw Error(`failed to fetch Wasm: ${e.status} ${e.statusText} fetching '${e.url}'`);if(typeof WebAssembly.instantiateStreaming==`function`)try{return await WebAssembly.instantiateStreaming(e,t)}catch(t){if(n(e.type)&&e.headers.get(`Content-Type`)!==`application/wasm`)console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n",t);else throw t}let r=await e.arrayBuffer();return await WebAssembly.instantiate(r,t)}{let n=await WebAssembly.instantiate(e,t);return n instanceof WebAssembly.Instance?{instance:n,module:e}:n}function n(e){switch(e){case`basic`:case`cors`:case`default`:return!0}return!1}}async function ki(e){if(L!==void 0)return L;e!==void 0&&(Object.getPrototypeOf(e)===Object.prototype?{module_or_path:e}=e:console.warn(`using deprecated parameters for the initialization function; pass a single object instead`)),e===void 0&&(e=new URL(`sketch_wasm_bg-BvHC48hh.wasm`,import.meta.url));let t=ui();(typeof e==`string`||typeof Request==`function`&&e instanceof Request||typeof URL==`function`&&e instanceof URL)&&(e=fetch(e));let{instance:n,module:r}=await Oi(await e,t);return Di(n,r)}var Ai=e=>typeof e==`object`&&!!e&&!(e instanceof Error)&&typeof e.message==`string`,ji=Object.freeze({version:li,engineRev:oi,detectFormat:ai,liveMolCount:si,parse:(e,t,n,r,i)=>ii.parse(e,t,n,r,i),useLegacyStereoPerception:ci,empty:()=>ii.empty(),abbreviationsJson:()=>ii.abbreviationsJson(),cxsmilesOf:e=>{let t=null;try{return t=ii.parse(e,`smiles`,!0),t.toCxsmiles()}catch(e){if(Ai(e))return null;throw e}finally{try{t?.free()}catch{}}}}),Mi;function Ni(){return Mi??=ki({module_or_path:new URL(`sketch_wasm_bg-BvHC48hh.wasm`,import.meta.url)}).then(()=>ji,e=>{throw Mi=void 0,e}),Mi}function Pi(e){let t,n=null;return{load:()=>t??=e().then(e=>e.start()).then(e=>n=e,e=>{throw t=void 0,e}),loaded:()=>n}}var Fi=Pi(()=>import(`./cdxml-C0oob9kG.js`)),Ii=Pi(()=>import(`./inchi-CzWcatSB.js`));function Li(){return Fi.load()}function Ri(){return Fi.loaded()}function zi(){return Ii.load()}function Bi(){return Ii.loaded()}var Vi={rest:`#ffffff`,lit:`#f4faff`,hoverRing:`#6b7a90`,hoverRingPx:1,pressedRing:`#1c65c0`,pressedRingPx:2},R={fill:`#ffffff`,ring:`#1a5fb4`,ringPx:1},Hi=Object.freeze({"alkali-metal":{fill:`#ffe1d6`,stripe:`#e8956c`},"alkaline-earth-metal":{fill:`#ffecc7`,stripe:`#e5b84f`},"transition-metal":{fill:`#fdf3c4`,stripe:`#d9c23a`},"post-transition-metal":{fill:`#e3eed8`,stripe:`#9cc17a`},metalloid:{fill:`#d7efe5`,stripe:`#62b894`},nonmetal:{fill:`#ffffff`,stripe:`#5b9bd5`},halogen:{fill:`#ffffff`,stripe:`#3fb8c4`},"noble-gas":{fill:`#e8e1f7`,stripe:`#a58fe0`},lanthanide:{fill:`#f7dfef`,stripe:`#de8ec2`},actinide:{fill:`#f3dde0`,stripe:`#d98a97`}}),Ui=Object.freeze({top:[{label:p.groupHistory,tools:[`undo`,`redo`]},{label:p.groupCanvas,tools:[`clear`]},{label:p.groupZoom,tools:[`zoom`]},{label:p.groupStructure,tools:[`clean`,`structure.stereo-panel`,`hydrogens`]},{label:p.groupClipboard,tools:[`cut`,`copy`,`paste`]},{label:p.groupSketcher,tools:[`options`]}],left:[{label:p.groupSelection,tools:[{palette:`select`,label:p.selectionTools,tools:[`select.rect`,`select.lasso`,`select.structure`]}]},{label:p.groupHandErase,tools:[`hand`,`erase`]},{label:p.groupBonds,tools:[`bond.single`,`bond.double`,`bond.triple`,`bond.aromatic`,{palette:`bond.query`,label:p.queryBonds,tools:[`bond.any`,`bond.single-or-double`,`bond.single-or-aromatic`,`bond.double-or-aromatic`]},{palette:`bond.stereo`,label:p.stereoBonds,tools:[`bond.wedge`,`bond.hash`,`bond.wavy`,`bond.crossed`]},`bond.dative`,`chain`]},{label:p.groupRGroups,tools:[{palette:`rgroup`,label:p.rgroupTools,tools:[`rgroup.label`,`rgroup.attachment`]}]},{label:p.groupSGroups,tools:[`sgroup`]},{label:p.groupCharges,tools:[`charge.plus`,`charge.minus`]}],right:[{label:p.groupElements,tools:[`element.h`,`element.c`,`element.n`,`element.o`,`element.s`,`element.p`,`element.f`,`element.cl`,`element.br`,`element.i`,`element.table`]}],bottom:[{label:p.groupRings,tools:[...Hn.map(Gr),Ar]}]}),Wi=/^[a-z][a-z0-9-]*(\.[a-z][a-z0-9-]*)*$/,Gi=e=>typeof e==`object`&&!!e&&!Array.isArray(e);function Ki(e,t){if(!Array.isArray(e))return[];let n=[];for(let r of e){if(!Gi(r)||!Array.isArray(r.tools))continue;let e=[];for(let n of r.tools){let r=n===`structure.enhanced-stereo`?qi:n;if(typeof r==`string`)Sn(r)&&r!==`element.query`&&!t.has(r)&&(t.add(r),r===`zoom`&&t.add(`zoom.level`),r===`zoom.level`&&t.add(`zoom`),e.push(Bt(r)?{kind:`readout`,id:r}:{kind:`button`,id:r}));else if(Gi(r)&&typeof r.palette==`string`&&Wi.test(r.palette)&&!Sn(r.palette)&&!t.has(r.palette)&&Array.isArray(r.tools)){let n=r.tools.filter(e=>typeof e==`string`&&Sn(e)&&e!==`element.query`&&!qn(e)&&!Bt(e)&&!vn(e)&&!Xt(e)&&!t.has(e));if(!n.length)continue;t.add(r.palette),n.forEach(e=>t.add(e)),e.push({kind:`palette`,id:r.palette,label:typeof r.label==`string`&&r.label.trim()?r.label:D(n[0]),members:n})}}e.length&&n.push({label:typeof r.label==`string`&&r.label.trim()?r.label:``,items:e})}return n}var qi=`structure.stereo-panel`,Ji=`structure.enhanced-stereo`;function Yi(e){let t=Gi(e)&&Array.isArray(e.hidden)?e.hidden:[];return new Set(t.filter(e=>typeof e==`string`&&(Sn(e)||rr.has(e))).flatMap(e=>e===`structure.enhanced-stereo`||e===`structure.stereo-panel`?[qi,Ji]:[e]))}var Xi=[`top`,`left`,`right`,`bottom`];function Zi(e,t=new Set){let n=Yi(e),r=Gi(e)&&Xi.some(t=>t in e),i=Gi(e)&&(r||!n.size)?e:Ui,a=new Set([...n,...t]);return{top:Ki(i.top,a),left:Ki(i.left,a),right:Ki(i.right,a),bottom:Ki(i.bottom,a)}}function Qi(e){return[...e.top,...e.left,...e.right,...e.bottom].flatMap(e=>e.items.flatMap(e=>e.kind===`palette`?e.members:[e.id]))}function $i(e){return new Set(Qi(e).map(mn).filter(e=>e!==null))}var ea=new Map([`groupHistory`,`groupCanvas`,`groupZoom`,`groupStructure`,`groupClipboard`,`groupSketcher`,`groupSelection`,`selectionTools`,`groupHandErase`,`groupBonds`,`stereoBonds`,`queryBonds`,`groupRGroups`,`rgroupTools`,`groupSGroups`,`groupCharges`,`groupElements`,`groupRings`].map(e=>[p[e],e]));function ta(e,t){let n=ea.get(e);return n?t[n]:e}function na(e,t){let n=e.members[0];return n!==void 0&&e.label===D(n)?D(n,t):ta(e.label,t)}function ra(e,t){return!(e>0)||!(t>0)?0:Math.max(1,Math.ceil(e/t-1e-9))}var ia=(e,t,n)=>n?A(e,t):Ir(A(e,t));function aa(e,t,n,r,i=1){let a=r*Math.cos(30*Math.PI/180),o=r*Math.sin(30*Math.PI/180);return Array.from({length:n},(n,r)=>{let s=r+1,c=O(e,t,s*a);return s%2?O(c,t+i*90,o):c})}function oa(e,t,n,r,i){let a=ia(e,t,r);return{direction:a,points:aa(e,a,ra(Jn(e,t),n),n,i?-1:1)}}function sa(e,t){let n=Math.max(12,.6*t);return{width:Math.max(n*1.2,n*(.62*e.length+.6)),height:n*1.35}}var ca=.97,la=2.4,ua=.5,da=3;function fa(e,t,n){return Math.hypot(Math.max(Math.abs(e.x-t.x)-n.width/2,0),Math.max(Math.abs(e.y-t.y)-n.height/2,0))}function pa(e){let[t,n]=[e.label.width,e.label.height],r=r=>({x:Math.min(Math.max(r.x,t/2),Math.max(t/2,e.area.width-t/2)),y:Math.min(Math.max(r.y,n/2),Math.max(n/2,e.area.height-n/2))});if(e.touch){let i=e.pointer.y-(30+n/2);if(i-n/2<0){let n=e.pointer.x<e.area.width/2?1:-1;return r({x:e.pointer.x+n*(30+t/2),y:e.pointer.y})}let a=[e.pointer.x,(e.pointer.x+e.end.x)/2,e.end.x].map(e=>r({x:e,y:i}));return a.find(t=>Math.hypot(t.x-e.end.x,t.y-e.end.y)<=la*e.bondPx)??a[a.length-1]}let i=e=>e*Math.PI/180,a={x:Math.cos(i(e.direction)),y:-Math.sin(i(e.direction))},o=e.flip?-1:1,s={x:-Math.sin(i(e.direction))*o,y:-Math.cos(i(e.direction))*o},c=t=>(t.x-e.origin.x)*a.x+(t.y-e.origin.y)*a.y,l=t=>(t.x-e.origin.x)*s.x+(t.y-e.origin.y)*s.y,u=e.bondPx,[d,f]=[Math.abs(a.x)*(t/2)+Math.abs(a.y)*(n/2),Math.abs(s.x)*(t/2)+Math.abs(s.y)*(n/2)],[p,m,h]=[c(e.end),ua*u/2,ua*u/2+da],g={x:e.origin.x+p*a.x,y:e.origin.y+p*a.y},_=[];for(let t of[-.9,1.05,-1.05,1.2,-1.2,1.4,-1.4])for(let n of[.3,0,.6,-.3,-.6]){let i=r({x:g.x-n*u*a.x+t*u*s.x,y:g.y-n*u*a.y+t*u*s.y}),o=Math.abs(l(i)-m)<f+h&&c(i)-d<p+da&&c(i)+d>-3,v=fa(e.pointer,i,e.label);_.push({at:i,near:Math.hypot(i.x-e.end.x,i.y-e.end.y),gap:v,clear:v>=8&&!o})}let v=e=>e.sort((e,t)=>e.near-t.near)[0]?.at??null;return v(_.filter(e=>e.clear&&e.near<=ca*u))??v(_.filter(e=>e.clear))??_.sort((e,t)=>t.gap-e.gap)[0].at}var ma=Object.freeze({points:Object.freeze([]),reach:0,labels:Object.freeze([])}),ha=(e,t,n)=>n<t?t:Math.min(n,Math.max(t,e));function ga(e,t,n=[]){if(!e.length&&!n.length)return null;let r=[...e.flatMap(e=>[e.x-t,e.x+t]),...n.flatMap(e=>[e.x,e.x+e.width])],i=[...e.flatMap(e=>[e.y-t,e.y+t]),...n.flatMap(e=>[e.y,e.y+e.height])];return{left:Math.min(...r),top:Math.min(...i),right:Math.max(...r),bottom:Math.max(...i)}}function _a(e,t){let n=1/0;for(let r of t.points)n=Math.min(n,Math.hypot(e.x-r.x,e.y-r.y)-t.reach);for(let r of t.labels){let t=Math.max(r.x-e.x,0,e.x-(r.x+r.width)),i=Math.max(r.y-e.y,0,e.y-(r.y+r.height));n=Math.min(n,Math.hypot(t,i))}return n}function va(e,t){let n=e.reach+12;return{points:e.points.filter(e=>Math.abs(e.y-t)<n),reach:e.reach,labels:e.labels.filter(e=>e.y-12<t&&t<e.y+e.height+12)}}function ya(e,t,n){let r=null;for(let i of e){let e=va(n,i.y),a=ha(t,i.from,i.to);for(let t=0;;t++){let n=(t===0?[a]:[a+t*2,a-t*2]).filter(e=>e>=i.from&&e<=i.to);if(!n.length)break;for(let t of n){let n={x:t,y:i.y},a=_a(n,e);if(a>=12)return{at:n,below:i.below};(!r||a>r.clear+.5)&&(r={at:n,below:i.below,clear:a})}}}return r??{at:{x:ha(t,e[0].from,e[0].to),y:e[0].y},below:e[0].below}}function ba(e,t,n=ma){let{width:r,height:i}=t,a=(e.left+e.right)/2,[o,s]=[16,r-4-12],c={from:ha(e.left,o,s),to:ha(e.right,o,s)},[l,u]=[16,i-4-12],d={y:e.top-24,below:!1,...c},f={y:e.bottom+24,below:!0,...c},p=d.y>=l&&d.y<=u,m=f.y>=l&&f.y<=u,h=[];p&&e.top-76>=4&&h.push(d),m&&e.bottom+76<=i-4&&h.push(f),p&&!h.includes(d)&&h.push(d),m&&!h.includes(f)&&h.push(f),h.push({y:ha(d.y,l,u),below:!1,from:o,to:s},{y:ha(f.y,l,u),below:!0,from:o,to:s});let{at:g,below:_}=ya(h,a,n),v=_?e.bottom:e.top,y=_?g.y-12:g.y+12,b=_?y-v:v-y,ee=g.x>=e.left&&g.x<=e.right&&b>0&&b<=12+1e-9?[{x:g.x,y:v},{x:g.x,y}]:null,te=_?g.y+12+8:g.y-12-8-32,ne=ha(g.x-48,4,r-4-96);return te>=4&&te+32<=i-4?{handle:g,stem:ee,actions:{x:ne,y:te},below:_,beside:!1}:{handle:g,stem:ee,actions:xa(g,t,_),below:_,beside:!0}}function xa(e,t,n){let{width:r,height:i}=t,a=ha(e.y-16,4,i-4-32),o=e.x+12+8;if(o+96<=r-4)return{x:o,y:a};let s=e.x-12-8-96;if(s>=4)return{x:s,y:a};let c=ha(e.x-48,4,r-4-96),l=e.y-12-8-32,u=e.y+12+8,d=l>=4,f=u+32<=i-4;return(n||!d)&&f?{x:c,y:u}:{x:c,y:ha(l,4,i-4-32)}}var Sa=(e,t,n)=>Math.hypot(e.x-t.x,e.y-t.y)<=(n?22:14);function Ca(e,t,n,r){let i=Math.hypot(e.x-t.x,e.y-t.y);return i>(n?14:9)?!1:r.every(t=>Math.hypot(e.x-t.x,e.y-t.y)>=i-1)}function wa(e,t,n){return!e.onAtom&&Sa(e.at,t.handle,e.touch)?`rotate`:e.alt&&Ca(e.at,t.centre,e.touch,n)?`centre`:null}var Ta=class{#e=null;#t=null;get lastLabel(){return this.#t}clicked(e){this.#e=e}tapCount(e,t){let n=this.#e;return!n||!n.touch?1:t-n.time<=300&&Math.hypot(e.x-n.at.x,e.y-n.at.y)<=10?2:1}double(e){let t=this.#e;return this.#e=null,t&&t.atomId===e?t.token:null}committed(e){let t=e.trim();t&&(this.#t=t)}forget(){this.#e=null}};function Ea(e,t){return Object.freeze({kind:`format`,format:e,...t.code?{code:t.code}:{},message:t.message,atoms:Object.freeze([...t.atoms]),bonds:Object.freeze([...t.bonds??[]])})}function Da(e){return(typeof e.readWarnings==`function`?e.readWarnings():[]).map(e=>Object.freeze({kind:`format`,format:e.format,code:e.code,message:e.message,atoms:Object.freeze([]),bonds:Object.freeze([])}))}function Oa(e){return typeof e==`string`?e.trim().toLowerCase()!==`false`:e===!0}var ka=e=>/^﻿?\s*</.test(e),Aa=Object.freeze({format:`reaction`,message:`Reactions are not supported`});function ja(e){return typeof WebAssembly<`u`&&e instanceof WebAssembly.RuntimeError?Na(e):Ma(e)&&e.kind===`reaction`?Aa:e}function Ma(e){return typeof e==`object`&&!!e&&!(e instanceof Error)&&typeof e.kind==`string`}function Na(e){let t=e instanceof Error?`${e.name}: ${e.message}`:String(e);return Object.freeze({format:`cdxml`,message:`Crux Sketch: the CDXML module failed (${t})`})}function Pa(e){let t=e instanceof Error?`${e.name}: ${e.message}`:String(e);return Object.freeze({format:`cdxml`,message:`Crux Sketch: the CDXML module could not be loaded (${t})`,code:`unavailable`})}var Fa=Object.freeze({format:`cdxml`,message:`Crux Sketch: the CDXML module is loading; await loadFormat("cdxml") before reading the cdxml getter`,code:`loading`}),Ia=Object.freeze({format:`cdxml`,message:`Crux Sketch: this sketcher was destroyed before the CDXML module loaded, so its drawing was never written as CDXML`});function La(e){return Object.freeze({kind:`format`,format:`cdxml`,code:e.kind,message:e.message,atoms:Object.freeze([]),bonds:Object.freeze([])})}function Ra(e,t){let n=Array.from({length:e},(e,t)=>t),r=e=>{for(;n[e]!==e;)e=n[e]=n[n[e]];return e},i=e;for(let a of t){if(a.a>=e||a.b>=e)continue;let[t,o]=[r(a.a),r(a.b)];t!==o&&(n[t]=o,--i)}return i}function za(e,t,n,r){let i;try{i=t.cdxmlToMol(n)}catch(e){throw ja(e)}let a=e.parse(i.molblock,`molV3000`,r);try{for(let e of i.aliases)a.setAtomAlias(e.atom,e.alias);let e=Yn(a.geometryJson()),t=Da(a).map(e=>Object.freeze({...e,format:`cdxml`}));return{mol:a,drops:[...i.readWarnings.map(La),...t],structures:Ra(e.atoms.length,e.bonds),left:i.readWarnings.map(e=>Object.freeze({kind:e.kind,count:e.count}))}}catch(e){try{a.free()}catch{}throw e}}function Ba(e){return{molblock:e.toMolblockV3000(),aliases:[...new Set(e.writeWarnings(`molV3000`).filter(e=>e.code===`alias`).flatMap(e=>e.atoms))].sort((e,t)=>e-t).flatMap(t=>{let n=JSON.parse(e.atomJson(t)).alias;return typeof n==`string`&&n?[{atom:t,alias:n}]:[]})}}function Va(e){return Ea(`cdxml`,e)}function Ha(e,t){let n=e.molToCdxml(t.molblock,t.aliases.map(e=>({atom:e.atom,alias:e.alias})));return Object.freeze({text:n.cdxml,warnings:Object.freeze(n.writeWarnings.map(Va))})}var Ua=[`structure`,`text`,`graphic`,`arrow`,`picture`,`bracket`];function Wa(e,t,n,r,i=[]){let a=e=>r.filter(t=>t.kind===e).reduce((e,t)=>e+t.count,0),o=[...Ua.filter(e=>a(e)>0).map(t=>e.leftOut(t,a(t))),...[...new Set(i)].map(t=>e.readLost(t))],s=a(`alias`);return!o.length&&!s?null:e.readReport(t,n,o,s)}function Ga(e,t=[]){let n=new Set(t.map(e=>e.kind));return[...new Set(e.flatMap(e=>e.code&&!n.has(e.code)?[e.code]:[]))]}var Ka=e=>/^\uFEFF?\s*InChI=/.test(e),qa=/^\uFEFF?\s*(?:InChIKey=)?[A-Z]{14}-[A-Z]{10}-[A-Z]\s*$/,Ja=e=>qa.test(e),Ya=Object.freeze({format:`inchi`,message:`An InChIKey cannot be drawn: it is a hash of an InChI, and no structure can be made of it`}),Xa=Object.freeze({format:`inchi`,message:`Crux Sketch: the InChI module is loading; await loadFormat("inchi") before reading the inchi or inchiKey getter`,code:`loading`});function Za(e){let t=e instanceof Error?`${e.name}: ${e.message}`:String(e);return Object.freeze({format:`inchi`,message:`Crux Sketch: the InChI module could not be loaded (${t})`,code:`unavailable`})}var Qa=Object.freeze({format:`inchi`,message:`Crux Sketch: this sketcher was destroyed before the InChI module loaded, so its drawing was never written as InChI`});function $a(e){let t=e instanceof Error?`${e.name}: ${e.message}`:String(e);return Object.freeze({format:`inchi`,message:`Crux Sketch: the InChI module failed (${t})`})}function eo(e){if(typeof WebAssembly<`u`&&e instanceof WebAssembly.RuntimeError)return $a(e);let t=e&&typeof e==`object`&&typeof e.message==`string`?e.message:String(e);return Object.freeze({format:`inchi`,message:t})}var to=new Set([`undefined-stereo`,`charges-rearranged`,`unusual-valence`,`ambiguous-stereo`,`library`]);function no(e){let t=e.replace(/^\uFEFF/,``);return t.split(/\r?\n/).map(e=>e.trim()).find(e=>e.startsWith(`InChI=`))??t.trim()}function ro(e,t,n){if(Ja(t))throw Ya;let r=no(t),i;try{i=e.inchiToMol(r).smiles}catch(e){throw eo(e)}let a=n(i),o=Da(a).map(e=>Object.freeze({...e,format:`inchi`})),s=io.exec(r);if(s){let e=s[1]===`3`?`racemic (/s3)`:`relative (/s2)`;o.push(Ea(`inchi`,{code:`stereo-groups`,message:`the InChI's ${e} stereo is read as absolute centres`,atoms:[]}))}return{mol:a,drops:o}}var io=/\/s([23])(?=\/|$)/;function ao(e){return Ea(`inchi`,{message:e.message,atoms:e.atoms,bonds:e.bonds,code:e.code})}function oo(e,t){let n;try{n=e.molToInchi(t)}catch(e){throw eo(e)}if(n.code===2||n.code===3)throw Object.freeze({format:`inchi`,message:`the InChI library could not write this drawing (${n.message||`return code ${n.code}`})`});let r=``;if(n.inchi)try{r=e.inchiToInchiKey(n.inchi)}catch(e){throw eo(e)}return Object.freeze({inchi:n.inchi,key:r,warnings:Object.freeze(n.writeWarnings.map(ao))})}function so(e){return[...new Set(e.warnings.flatMap(e=>e.code&&!to.has(e.code)?[e.code]:[]))]}var co={frame(e){let t=requestAnimationFrame(()=>e());return()=>cancelAnimationFrame(t)},painted(e){let t,n=requestAnimationFrame(()=>{t=setTimeout(e,0)});return()=>{cancelAnimationFrame(n),t!==void 0&&clearTimeout(t)}},after(e,t){let n=setTimeout(t,e);return()=>clearTimeout(n)}},lo=Object.freeze({atoms:Object.freeze([]),bonds:Object.freeze([]),bondLength:0});function uo(e,t,n=null){if(!e||!e.atoms.length)return lo;let r=n?.k??1,i=(e,i)=>Object.freeze({x:r*e+(n?.x??0)+t.x,y:r*i+(n?.y??0)+t.y}),a=e.atoms.map(e=>e.hidden?null:i(e.px.x,e.px.y)),o=e.bonds.map(t=>{let[n,r]=[e.atoms[t.ends[0]],e.atoms[t.ends[1]]];return t.hidden||!n||!r||n.hidden||r.hidden||t.ends[0]===t.ends[1]?null:i((n.px.x+r.px.x)/2,(n.px.y+r.px.y)/2)});return Object.freeze({atoms:Object.freeze(a),bonds:Object.freeze(o),bondLength:r*e.bondPx})}function fo(e){return e.ctrlKey||e.metaKey||e.altKey?null:/^[A-Za-z0-9]$/.test(e.key)?e.key:null}var po=class{#e=null;press(e,t,n=!0){let r=this.#e;return n&&r&&t-r.time<=500&&r.char.toLowerCase()!==e.toLowerCase()?(this.#e=null,r.char+e):(this.#e={char:e,time:t},null)}reset(){this.#e=null}},mo=`abbreviation.place`;function ho(e,t,n){let r=e.labelGroupOf(t);return r?o(r,n):S(`atom.label`,[{op:`setAtomLabel`,atom:t,text:n}])}var go=(e,t,n=null)=>({op:`insertAbbreviation`,atom:n,label:e,x:t.x,y:t.y});function _o(e,t,n,r){return t?t.kind===`atom`?ho(e,t.index,r):null:S(mo,[go(r,n)])}function vo(e,t,n,r){if(`point`in n)return S(mo,[go(r,n.point,t)]);if(n.atom===t)return null;let i=e.bondBetween(t,n.atom)>=0?[]:[{op:`addBond`,a:t,b:n.atom,type:`single`,stereo:`none`}];return S(mo,[...i,...ho(e,n.atom,r)?.ops??[]])}var yo=(e,t,n,r)=>t>=0?ho(e,t,r):S(mo,[go(r,n)]);function bo(e,t,n,r,i=!1){let a=Oe(n.atoms);if(!a)return null;let o=Er(n.atoms,n.bonds,e.docLen),s={left:o*a.left,bottom:o*a.bottom,right:o*a.right,top:o*a.top},c=i?_r(s,ir(e.atoms.map(e=>e.at),e.bonds),r,e.docLen):{x:r.x-(s.left+s.right)/2,y:r.y-(s.bottom+s.top)/2};return{op:`paste`,text:t.molblock,format:`molV2000`,kekule:n.kekule,scale:o,dx:c.x,dy:c.y}}function xo(e,t){let n=new Set([t]),r=[t];for(;r.length;){let t=r.pop();for(let i of e.bondsOf(t)){let a=e.other(i,t);n.has(a)||(n.add(a),r.push(a))}}return[...n].sort((e,t)=>e-t)}var So=1.25;function Co(e,t,n,r){let i=-1/0;for(let n of xo(e,t)){let t=e.atoms[n];i=Math.max(i,t.px.x,t.labelBox?t.labelBox.x+t.labelBox.width:t.px.x)}let a=e.toMolecule({x:i,y:0}).x,o=Oe(n.atoms),s=o?Er(n.atoms,n.bonds,e.docLen)*(o.right-o.left)/2:0;return{x:a+So*e.docLen+s,y:r}}function wo(e,t){let n=e.contractedGroupOf(t);return n?n.bonds.length===0?[t]:null:e.groupOf(t)?null:e.bondsOf(t).length===0?[t]:null}function To(e,t,n,r,i){if(!t)return Eo(e,n,r,i);let a=t.kind===`atom`?t.index:e.bonds[t.index]?.a;if(a===void 0||!e.atoms[a])return null;let o=t.kind===`atom`?wo(e,a):null;if(o){let t=bo(e,r,i,e.atoms[a].at);return t&&S(`salt.place`,[t,...gn(e,o,[])])}let s=bo(e,r,i,Co(e,a,i,n.y),!0);return s&&S(`salt.place`,[s])}function Eo(e,t,n,r){let i=bo(e,n,r,t);return i&&S(`salt.place`,[i])}var Do=.35,Oo=3,ko=3/19.2,Ao=.25,jo=1;function Mo(e,t,n){let[r,i]=[n.x-t.x,n.y-t.y],a=r*r+i*i,o=a?Math.max(0,Math.min(1,((e.x-t.x)*r+(e.y-t.y)*i)/a)):0;return Math.hypot(e.x-(t.x+o*r),e.y-(t.y+o*i))}function No(e,t){let n=Math.max(t.x-e.x,0,e.x-(t.x+t.width)),r=Math.max(t.y-e.y,0,e.y-(t.y+t.height));return Math.hypot(n,r)}var Po=Object.freeze({atom:12,bond:8});function Fo(e,t,n,r,i={}){let a=Math.max(12,Do*r)+(i.slop?.atom??0),o=jo+(i.slop?.atom??0),s=12+(i.slop?.bond??0),c=Math.min(Math.max(Oo,ko*r),Ao*r),l=null,u=1/0;return e.forEach((e,t)=>{if(i.excludeAtoms?.has(t)||e.hidden)return;let r=Math.hypot(n.x-e.x,n.y-e.y),s=1/0;for(let t of e.glyphs)s=Math.min(s,No(n,t));let d=s<=jo;if(r>a&&s>o)return;let f=d?-c:r<=a?r-c:s;f<u&&(u=f,l={kind:`atom`,index:t})}),i.atomsOnly||t.forEach((t,r)=>{let[i,a]=[e[t.a],e[t.b]];if(!i||!a||i.hidden||a.hidden)return;let o=Mo(n,i,a);o>s||o<u&&(u=o,l={kind:`bond`,index:r})}),l}function Io(e,t){let n=Math.max(e.x-(t.x+t.width),0,t.x-(e.x+e.width)),r=Math.max(e.y-(t.y+t.height),0,t.y-(e.y+e.height));return Math.hypot(n,r)}function Lo(e,t,n){let[r,i]=[t.x-e.x,t.y-e.y],[a,o]=[0,1],s=[[-r,e.x-n.x],[r,n.x+n.width-e.x],[-i,e.y-n.y],[i,n.y+n.height-e.y]];for(let[e,t]of s){if(e===0){if(t<0)return!1;continue}let n=t/e;if(e<0){if(n>o)return!1;a=Math.max(a,n)}else{if(n<a)return!1;o=Math.min(o,n)}}return!0}function Ro(e,t,n){if(Lo(e,t,n))return 0;let r=[{x:n.x,y:n.y},{x:n.x+n.width,y:n.y},{x:n.x,y:n.y+n.height},{x:n.x+n.width,y:n.y+n.height}];return Math.min(No(e,n),No(t,n),...r.map(n=>Mo(n,e,t)))}function zo(e,t,n,r,i={}){let a=Math.max(12,Do*r)+(i.slop?.atom??0),o=jo+(i.slop?.atom??0),s=12+(i.slop?.bond??0);for(let t of e)if(!t.hidden){if(No(t,n)<=a)return!0;for(let e of t.glyphs)if(Io(e,n)<=o)return!0}for(let r of t){let[t,i]=[e[r.a],e[r.b]];if(t&&i&&!t.hidden&&!i.hidden&&Ro(t,i,n)<=s)return!0}return!1}var z=(e,t,n=!1)=>({key:e,scope:`atom`,action:`atom.element`,element:t,...n?{alias:n}:{}}),Bo=e=>({key:e,scope:`atom`,action:`atom.sprout`,sprout:e}),Vo=(e,t)=>({key:e,scope:`bond`,action:`bond.fuse`,ring:t}),Ho=(e,t,n)=>[{key:e,scope:`selection`,action:`selection.nudge`,nudge:{x:t*1,y:n*1}},{key:e,shift:!0,scope:`selection`,action:`selection.nudge`,nudge:{x:t*10,y:n*10}}],Uo=(e,t)=>[...[`atom`,`bond`].flatMap(n=>[{key:e,scope:n,action:`hotspot.walk`,arrow:t},{key:e,shift:!0,scope:n,action:`hotspot.jump`,arrow:t},{key:e,shift:!0,alt:!0,scope:n,action:`hotspot.hop`,arrow:t}]),{key:e,scope:`generic`,action:`hotspot.start`,arrow:t}],Wo=(e,t)=>({key:e,alt:!0,scope:`any`,action:`selection.rotate`,turn:t}),Go=Object.freeze([{key:`1`,scope:`bond`,action:`bond.single`},{key:`2`,scope:`bond`,action:`bond.double`},{key:`3`,scope:`bond`,action:`bond.triple`},{key:`w`,scope:`bond`,action:`bond.wedge`},{key:`h`,scope:`bond`,action:`bond.hash`},{key:`W`,scope:`bond`,action:`bond.hash`},{key:`y`,scope:`bond`,action:`bond.wavy`},{key:`l`,scope:`bond`,action:`bond.left`},{key:`c`,scope:`bond`,action:`bond.centre`},{key:`r`,scope:`bond`,action:`bond.right`},Vo(`a`,`benzene`),Vo(`4`,`cyclobutane`),Vo(`5`,`cyclopentane`),Vo(`6`,`cyclohexane`),Vo(`7`,`cycloheptane`),Vo(`8`,`cyclooctane`),{key:`x`,scope:`generic`,action:`tool.bond.single`},{key:`j`,scope:`generic`,action:`tool.ring.benzene`},{key:`J`,scope:`generic`,action:`tool.ring.cyclopentadiene`},{key:`X`,scope:`generic`,action:`tool.chain`},...[`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`,`0`,`a`].map(Bo),...[`K`,`v`,`u`,`z`].map(Bo),...[`Z`,`M`,`k`].map(Bo),z(`c`,`C`),z(`n`,`N`),z(`w`,`N`,!0),z(`o`,`O`),z(`q`,`O`,!0),z(`s`,`S`),z(`p`,`P`),z(`f`,`F`),z(`h`,`H`),z(`i`,`I`),{key:`d`,scope:`atom`,action:`atom.deuterium`},z(`l`,`Cl`),z(`C`,`Cl`,!0),z(`b`,`Br`),z(`B`,`B`),z(`S`,`Si`),z(`L`,`Li`),{key:`Backspace`,scope:`atom`,action:`atom.unlabel`},{key:`Delete`,scope:`atom`,action:`atom.unlabel`},{key:`/`,scope:`atom`,action:`atom.properties`},{key:`?`,scope:`atom`,action:`atom.properties`},{key:`Enter`,scope:`atom`,action:`atom.edit`},{key:`r`,scope:`atom`,action:`atom.rgroup`},{key:`+`,scope:`atom`,action:`atom.charge.plus`},{key:`-`,scope:`atom`,action:`atom.charge.minus`},{key:`Delete`,ctrl:!0,scope:`any`,action:`clear`},{key:`Backspace`,ctrl:!0,scope:`any`,action:`clear`},{key:`z`,ctrl:!0,scope:`any`,action:`undo`},{key:`z`,ctrl:!0,shift:!0,scope:`any`,action:`redo`},{key:`y`,ctrl:!0,scope:`any`,action:`redo`,alias:!0},{key:`=`,ctrl:!0,scope:`any`,action:`zoom.in`},{key:`+`,ctrl:!0,scope:`any`,action:`zoom.in`,alias:!0},{key:`-`,ctrl:!0,scope:`any`,action:`zoom.out`},{key:`0`,ctrl:!0,scope:`any`,action:`zoom.reset`},{key:`F5`,scope:`any`,action:`zoom.reset`,canvas:!0},{key:`F6`,scope:`any`,action:`zoom.fit`,canvas:!0},{key:`F7`,scope:`any`,action:`zoom.in`,canvas:!0},{key:`F8`,scope:`any`,action:`zoom.out`,canvas:!0},{key:`h`,ctrl:!0,alt:!0,scope:`any`,action:`tool.hand`},{key:`k`,ctrl:!0,shift:!0,scope:`any`,action:`clean`},{key:`a`,ctrl:!0,scope:`any`,action:`select.all`},{key:`a`,ctrl:!0,shift:!0,scope:`any`,action:`select.none`},{key:` `,scope:`generic`,action:`select.last`},{key:`Delete`,scope:`selection`,action:`selection.delete`},{key:`Backspace`,scope:`selection`,action:`selection.delete`},...Ho(`ArrowUp`,0,1),...Ho(`ArrowDown`,0,-1),...Ho(`ArrowLeft`,-1,0),...Ho(`ArrowRight`,1,0),{key:`j`,ctrl:!0,scope:`any`,action:`select.join`},Wo(`ArrowUp`,15),Wo(`ArrowDown`,-15),Wo(`ArrowLeft`,1),Wo(`ArrowRight`,-1),{key:`h`,alt:!0,scope:`any`,action:`flip.horizontal`},{key:`v`,alt:!0,scope:`any`,action:`flip.vertical`},{key:`h`,ctrl:!0,shift:!0,scope:`any`,action:`rotate180.horizontal`},{key:`v`,ctrl:!0,shift:!0,scope:`any`,action:`rotate180.vertical`},{key:`x`,ctrl:!0,scope:`any`,action:`cut`},{key:`c`,ctrl:!0,scope:`any`,action:`copy`},{key:`v`,ctrl:!0,scope:`any`,action:`paste`},{key:`c`,ctrl:!0,alt:!0,scope:`any`,action:`copy.as`,copy:`smiles`},{key:`o`,ctrl:!0,alt:!0,scope:`any`,action:`copy.as`,copy:`molV2000`},{key:`o`,ctrl:!0,alt:!0,shift:!0,scope:`any`,action:`copy.as`,copy:`molV3000`},{key:`d`,ctrl:!0,scope:`any`,action:`copy.as`,copy:`cdxml`},{key:`p`,ctrl:!0,alt:!0,scope:`any`,action:`paste.as`,paste:`smiles`},{key:`p`,ctrl:!0,alt:!0,shift:!0,scope:`any`,action:`paste.as`,paste:`mol`},{key:`/`,scope:`bond`,action:`bond.properties`},{key:`?`,scope:`bond`,action:`bond.properties`},{key:`F10`,shift:!0,scope:`any`,action:`menu.open`,canvas:!0},{key:`ContextMenu`,scope:`any`,action:`menu.open`,canvas:!0},...Uo(`ArrowUp`,`up`),...Uo(`ArrowDown`,`down`),...Uo(`ArrowLeft`,`left`),...Uo(`ArrowRight`,`right`),{key:`Enter`,scope:`selection`,action:`selection.end`,idle:!0},{key:`Enter`,scope:`generic`,action:`hotspot.start`},{key:`Enter`,shift:!0,scope:`any`,action:`structure.new`,canvas:!0},{key:` `,scope:`atom`,action:`select.hotspot`},{key:` `,scope:`bond`,action:`select.hotspot`},{key:`.`,scope:`atom`,action:`atom.attachment`},{key:`?`,scope:`generic`,action:`keys.open`},...Kr.map(([e,t])=>({key:e,scope:`atom`,action:`atom.abbreviation`,abbreviation:t})),{key:`=`,scope:`atom`,action:`atom.nickname`},{key:`a`,alt:!0,scope:`any`,action:`structure.aromatize`},{key:`a`,ctrl:!0,alt:!0,scope:`any`,action:`structure.dearomatize`}]),Ko=(e,t)=>({key:e,scope:`generic`,action:`tool.choose`,tool:t}),qo=(e,t)=>({key:e,scope:`generic`,action:`tool.cycle`,tools:t}),Jo=[[`h`,`H`],[`c`,`C`],[`n`,`N`],[`o`,`O`],[`s`,`S`],[`p`,`P`],[`f`,`F`],[`l`,`Cl`],[`b`,`Br`],[`i`,`I`],[`S`,`Si`],[`N`,`Na`],[`B`,`B`]],Yo=[`single`,`double`,`triple`,`aromatic`],Xo=Object.freeze([...Go.filter(e=>!(e.scope===`atom`&&e.action===`atom.sprout`&&/^[0-9]$/.test(e.key))),...Yo.map((e,t)=>({key:String(t+1),scope:`atom`,action:`atom.bond`,bondType:e})),...Jo.map(([e,t])=>Ko(e,cr(t))),qo(`1`,[`bond.single`,`bond.wedge`,`bond.hash`,`bond.wavy`]),qo(`2`,[`bond.double`,`bond.crossed`]),Ko(`3`,`bond.triple`),Ko(`4`,`bond.aromatic`),qo(`t`,Hn.map(Gr)),{key:`T`,scope:`generic`,action:`library.open`},Ko(`Delete`,`erase`),Ko(`Backspace`,`erase`),{key:`l`,ctrl:!0,scope:`any`,action:`layout`,canvas:!0},{key:`l`,ctrl:!0,shift:!0,scope:`any`,action:`clean`,canvas:!0},{key:`m`,ctrl:!0,shift:!0,scope:`any`,action:`copy.as`,copy:`molV2000`}]),Zo=Object.freeze([{key:`x`,scope:`atom`,action:`atom.generic`,generic:`X`}]),Qo=Object.freeze({chemdraw:Go,ketcher:Xo}),$o=e=>typeof e==`string`&&e.trim().toLowerCase()===`ketcher`?`ketcher`:`chemdraw`,es=e=>!!e.ctrl&&/^[a-z]$/.test(e.key),ts=e=>!!e.alt&&!e.ctrl&&/^[a-z]$/.test(e.key);function ns(e,t){return e.key.toLowerCase()===t||!/^[a-z]$/i.test(e.key)&&e.code===`Key${t.toUpperCase()}`}function rs(e,t,n,r=!1){let i=!!(t.ctrlKey||t.metaKey),a=!!t.altKey;for(let o of e)if(!!o.ctrl===i&&!!o.alt===a&&(!o.canvas||r)){if(es(o)){if(!(o.alt?ns(t,o.key):o.key===t.key.toLowerCase())||!!o.shift!=!!t.shiftKey)continue}else if(ts(o)){if(!ns(t,o.key)||!!o.shift!=!!t.shiftKey)continue}else if(o.key!==t.key||o.key.length>1&&!!t.shiftKey!=!!o.shift)continue;if(o.scope===`any`||o.scope===n)return o}return null}function is(e,t,n,r,i=!1){let a=rs(e,t,n??`generic`,i);return a||!n||r?a:rs(e,t,`generic`,i)}function as(e,t,n=!1,r=!1){let i=rs(e,t,`selection`,n);return i?.scope===`selection`&&!(i.idle&&r)?i:null}var os={"bond.single":`single`,"bond.double":`double`,"bond.triple":`triple`,"bond.wedge":`wedge`,"bond.hash":`hash`};function ss(e,t,n,r){let i=rs(e,t,`bond`),a=i?.scope===`bond`?os[i.action]:void 0;if(a&&r)return{kind:`bond`,type:a};let o=rs(e,t,`atom`);return o?.scope===`atom`&&o.action===`atom.element`&&o.element?n?{kind:`element`,symbol:o.element}:{kind:`none`}:o?.scope===`atom`&&(o.action===`atom.sprout`||o.action===`atom.bond`||o.action===`atom.abbreviation`)||i?.scope===`bond`&&(a||i.action===`bond.fuse`)?{kind:`none`}:null}var cs=Object.freeze({Backspace:`⌫`,Delete:`⌦`,Enter:`↩`,Escape:`⎋`,Tab:`⇥`,ArrowUp:`↑`,ArrowDown:`↓`,ArrowLeft:`←`,ArrowRight:`→`}),ls=e=>!e.ctrl&&!e.alt&&/^[A-Z]$/.test(e.key);function us(e,t=`Ctrl`){let n=!!e.shift||ls(e),r=e.key===` `?`Space`:e.key,i=!!(e.ctrl||e.alt||n),a=i&&r.length===1?r.toUpperCase():r;return t===`⌘`?`${e.alt?`⌥`:``}${n?`⇧`:``}${e.ctrl?`⌘`:``}${cs[e.key]??a}`:i?[...e.ctrl?[t]:[],...e.alt?[`Alt`]:[],...n?[`Shift`]:[],a].join(`+`):r}function ds(e,t,n=`Ctrl`,r=` or `){let i=e.filter(e=>e.action===t&&!e.alias).map(e=>us(e,n));return i.length?i.join(r):void 0}var fs=(e,t,n)=>e.filter(e=>e.action===t&&(t===`copy.as`?e.copy:e.paste)===n&&!e.alias);function ps(e,t,n,r=`Ctrl`,i=` or `){let a=fs(e,t,n).map(e=>us(e,r));return a.length?a.join(i):void 0}function ms(e,t,n,r=!1){let i=fs(e,t,n).map(e=>gs(e,r));return i.length?i.join(` `):void 0}function hs(e,t,n=!1){let r=e.filter(e=>e.action===t&&!e.alias).map(e=>gs(e,n));return r.length?r.join(` `):void 0}var gs=(e,t=!1)=>[...e.ctrl?[t?`Meta`:`Control`]:[],...e.alt?[`Alt`]:[],...e.shift||ls(e)?[`Shift`]:[],e.key===`+`?`Plus`:e.key===` `?`Space`:e.key.length===1?e.key.toUpperCase():e.key].join(`+`);function _s(e,t,n,r=`Ctrl`,i=!1){let a=e.find(e=>e.action===t&&!e.alias&&(n===void 0||(t===`copy.as`?e.copy:e.paste)===n));return a?{label:us(a,r),aria:gs(a,i)}:void 0}function vs(e,t,n=`Ctrl`){let r=e.find(e=>e.action===`atom.element`&&e.element===t&&!e.alias);return r?us(r,n):void 0}function ys(e){switch(e.action){case`tool.bond.single`:return[`bond.single`];case`tool.ring.benzene`:return[`ring.benzene`];case`tool.ring.cyclopentadiene`:return[`ring.cyclopentadiene`];case`tool.chain`:return[`chain`];case`tool.hand`:return[`hand`];case`tool.choose`:return e.tool?[e.tool]:[];case`tool.cycle`:return e.tools??[];default:return[]}}function bs(e,t,n=`Ctrl`,r=` or `,i){let a=e.filter(e=>!e.alias&&(e.action===`tool.cycle`?i?(e.tools??[]).includes(t):e.tools?.[0]===t:ys(e).includes(t))).map(e=>i&&e.action===`tool.cycle`?i(us(e,n),e.tools?.[0]===t):us(e,n));return a.length?[...new Set(a)].join(r):void 0}var xs=Object.freeze({clear:`clear`,undo:`undo`,redo:`redo`,"zoom.in":`zoom.in`,"zoom.out":`zoom.out`,"zoom.reset":`zoom.reset`,"zoom.fit":`zoom.fit`,layout:`layout`,clean:`clean`,cut:`cut`,copy:`copy`,paste:`paste`,"keys.open":`keys`,"settings.open":`settings`,"library.open":`template.library`});function Ss(e,t){return t.size?e.flatMap(e=>{let n=xs[e.action];if(n&&t.has(n))return[];if(e.action===`tool.cycle`){let n=(e.tools??[]).filter(e=>!t.has(e));return n.length?[n.length===e.tools?.length?e:{...e,tools:n}]:[]}let r=ys(e);return r.length&&r.every(e=>t.has(e))?[]:[e]}):e}var Cs=new Set([`clear`,`undo`,`redo`,`zoom.in`,`zoom.out`,`zoom.reset`,`zoom.fit`,`layout`,`clean`,`select.all`,`select.none`,`select.join`,`flip.horizontal`,`flip.vertical`,`rotate180.horizontal`,`rotate180.vertical`,`keys.open`,`settings.open`,`library.open`]),ws=new Set([`atom.properties`,`atom.edit`,`atom.rgroup`,`atom.charge.plus`,`atom.charge.minus`,`atom.attachment`,`atom.unlabel`,`atom.deuterium`]),Ts=new Set([`bond.single`,`bond.double`,`bond.triple`,`bond.wedge`,`bond.hash`,`bond.wavy`,`bond.properties`]),Es=[`Backspace`,`Delete`,`Enter`,`Escape`,`Tab`,`ArrowUp`,`ArrowDown`,`ArrowLeft`,`ArrowRight`,`Home`,`End`,`PageUp`,`PageDown`,`Insert`,`ContextMenu`,...Array.from({length:12},(e,t)=>`F${t+1}`)];function Ds(e){if(typeof e!=`string`||!e.trim())return null;let t=e.trim(),n=t===`+`?[`+`]:t.endsWith(`++`)?[...t.slice(0,-2).split(`+`),`+`]:t.split(`+`),r=n.pop()??``,[i,a,o]=[!1,!1,!1];for(let e of n){let t=e.trim().toLowerCase();if([`ctrl`,`control`,`cmd`,`command`,`meta`,`mod`,`⌘`].includes(t))i=!0;else if([`alt`,`option`,`⌥`].includes(t))a=!0;else if([`shift`,`⇧`].includes(t))o=!0;else return null}let s=r===`Space`||r===`space`?` `:r===`Plus`||r===`plus`?`+`:r;if(s.length===1)/^[a-z]$/i.test(s)?i||a?s=s.toLowerCase():o&&([s,o]=[s.toUpperCase(),!1]):!i&&!a&&(o=!1);else{let e=Es.find(e=>e.toLowerCase()===s.toLowerCase());if(!e)return null;s=e}return{key:s,...i?{ctrl:i}:{},...a?{alt:a}:{},...o?{shift:o}:{}}}function Os(e){let t=!!(e.ctrl||e.alt)&&/^[a-z]$/i.test(e.key),n=t?e.key.toLowerCase():e.key,r=(t||e.key.length>1)&&!!e.shift;return`${e.ctrl?`C`:``}${e.alt?`A`:``}${r?`S`:``}:${n}`}function ks(e,t){if(Cs.has(t))return{...e,scope:`any`,action:t};if(ws.has(t))return{...e,scope:`atom`,action:t};if(Ts.has(t))return{...e,scope:`bond`,action:t};if(t===`selection.delete`)return{...e,scope:`selection`,action:`selection.delete`};let n=/^tool\.(.+)$/.exec(t)?.[1];if(n&&Xr(n)&&n!==`template.library`)return{...e,scope:`generic`,action:`tool.choose`,tool:n};let r=/^element\.([A-Z][a-z]{0,2})$/.exec(t)?.[1];return r&&Xr(cr(r))?{...e,scope:`atom`,action:`atom.element`,element:r}:null}function As(e){if(typeof e!=`object`||!e||Array.isArray(e))return Object.freeze({});let t={};for(let[n,r]of Object.entries(e)){let e=Ds(n);e&&(r===null?t[n]=null:typeof r==`string`&&ks(e,r)&&(t[n]=r))}return Object.freeze(t)}function js(e,t){let n=new Set,r=[];for(let[e,i]of Object.entries(t)){let t=Ds(e);if(!t)continue;n.add(Os(t));let a=i===null?null:ks(t,i);a&&r.push(a)}return n.size?[...r,...e.filter(e=>!n.has(Os(e)))]:e}function Ms(e,t){return e===`single`&&(t===`up`||t===`wedge`)?`wedge`:e===`single`&&(t===`down`||t===`hash`)?`hash`:e===`double`&&(t===`either`||t===`crossed`)?`crossed`:t===`either`||t===`wavy`?`wavy`:e===`single`||e===`double`||e===`triple`||e===`aromatic`||e===`dative`?e:`other`}function Ns(e,t){switch(t){case`single`:return e.typeSingle;case`double`:return e.typeDouble;case`triple`:return e.typeTriple;case`aromatic`:return e.typeAromatic;case`wedge`:return e.typeWedge;case`hash`:return e.typeHash;case`wavy`:return e.typeWavy;case`crossed`:return e.typeCrossed;case`dative`:return e.typeDative;default:return e.typeOther}}function Ps(e,t,n){if(!Number.isInteger(n)||n===0)return t;let r=Math.abs(n),i=RegExp(`${r===1?``:r}[+−-]$`),a=i.test(t)?t.replace(i,``):t;return e.format(n>0?`chargeSaidPlus`:`chargeSaidMinus`,{label:a||t,count:r})}var Fs=(e,t,n,r)=>e.format(`hotAtom`,{label:t,number:n+1,bonds:e.bonds(r)}),Is=(e,t,n,r,i)=>e.format(`hotBond`,{number:t+1,type:Ns(e,n),first:Math.min(r,i)+1,second:Math.max(r,i)+1});function Ls(e,t,n){return t||n?e.format(`selected`,{what:e.atomsAndBonds(t,n)}):e.selectionCleared}function Rs(e,t){let n=t??``,r=n.split(`.`)[0];return n===`clipboard.paste`?e.paste:n===`clipboard.drop`||n===`api.add`?e.stepAdd:n===`selection.delete`||n===`atom.erase`||n===`bond.erase`?e.stepDelete:n===`clear`?e.clearCanvas:n===`layout`?e.layout:n===`clean`?e.cleanUp:n===`straighten`?e.straighten:n===`bond.ez`?e.ezGroup:n===`stereo.chiral-flag`?e.chiralFlag:n===`stereo.group`?e.enhancedStereo:n===`stereo.centre`?e.stereochemistry:n===`stereo.atropisomer`?e.atropisomerGroup:n===`stereo.remove`?e.removeStereo:n===`stereo.rewedge`?e.recomputeWedges:n===`structure.aromatize`?e.aromatize:n===`structure.dearomatize`?e.dearomatize:n===`bond.position`?e.doubleBondPosition:n===`bond.reverse`?e.changeDirection:n===`bond.change`?e.stepBond:n===`atom.element`||n===`atom.label`?e.stepLabel:n===`atom.properties`?e.atomProperties:n===`atom.rgroup`?e.rgroupDialog:n===`atom.attachment-points`?e.attachmentPoint:n===`atom.charge`?e.charge:n===`atom.join`?e.stepJoin:n===`abbreviation.expand`?e.expandAbbreviation:n===`abbreviation.contract`?e.contractAbbreviation:n===`abbreviation.remove`?e.removeAbbreviation:n===`abbreviation.contract-label`?e.contractLabel:n===`abbreviation.replace`?e.stepLabel:n===`sgroup.add`||n===`sgroup.edit`?e.stepSGroup:n===`sgroup.remove`?e.removeSGroup:n===`selection.rotate`||n.endsWith(`.rotate180`)?e.stepRotate:n.endsWith(`.flip`)?e.stepFlip:n.endsWith(`.move`)||n.endsWith(`.merge`)||n===`selection.nudge`?e.stepMove:n.startsWith(`selection.copy`)?e.stepDuplicate:r===`ring`?e.stepRing:r===`template`?e.stepTemplate:r===`chain`?e.chain:r===`bond`||n===`atom.place`||n===`atom.grow`||n===`atom.sprout`?e.stepAdd:e.stepEdit}var zs=.5;function Bs(e,t){let n=new Map(e.atoms.map(e=>[e.id,e])),r=new Map(e.bonds.map(e=>[e.id,e])),i=new Set(t.atoms.map(e=>e.id)),a=new Set(t.bonds.map(e=>e.id)),o={added:[],removed:0,addedBonds:0,removedBonds:0,atoms:[],bonds:[],moved:0};return t.atoms.forEach((e,t)=>{let r=n.get(e.id);r?(r.label!==e.label&&o.atoms.push(t),Math.hypot(r.x-e.x,r.y-e.y)>zs&&(o.moved+=1)):o.added.push(e)}),t.bonds.forEach((e,t)=>{let n=r.get(e.id);n?n.type!==e.type&&o.bonds.push(t):o.addedBonds+=1}),o.removed=e.atoms.filter(e=>!i.has(e.id)).length,o.removedBonds=e.bonds.filter(e=>!a.has(e.id)).length,o}var Vs=(e,t,n)=>e.atomsAndBonds(t,n);function Hs(e,t,n,r){if(t.kind===`undo`)return e.format(`undone`,{name:Rs(e,t.label)});if(t.kind===`redo`)return e.format(`redone`,{name:Rs(e,t.label)});if(t.kind===`write`||!n)return e.format(`newDrawing`,{what:r.atoms.length?Vs(e,r.atoms.length,r.bonds.length):e.descriptionEmpty});let i=t.said??null;return i?.kind===`tool`?e.format(`toolThen`,{tool:e.format(`toolChosen`,{tool:i.tool}),what:Us(e,{kind:`edit`,label:t.label},n,r)}):Us(e,t,n,r)}function Us(e,t,n,r){let i=t.label??``,a=t.said??null;if(a?.kind===`ez`)return e.format(`bondNow`,{number:a.bond+1,type:a.choice===`unknown`?e.typeCrossed:a.choice});if(a?.kind===`flag`)return e.format(a.on?`switchedOn`:`switchedOff`,{name:e.chiralFlag});if(a?.kind===`centre`)return e.format(`atomNow`,{number:a.atom+1,label:a.config??e.stereoUnspecified});if(a?.kind===`axis`)return e.format(`bondNow`,{number:a.bond+1,type:a.label});if(a?.kind===`stereo-removed`||a?.kind===`rewedged`)return Ws(e,a);if(a?.kind===`aromaticity`)return Mt(e,a.command,a.rings,a.refused,a.reason);if(a?.kind===`cleaned`)return e.format(`cleanedAtoms`,{what:e.atoms(a.atoms)});let o=Bs(n,r);if(i===`clear`)return e.canvasCleared;if(i===`layout`)return e.laidOut;if(i===`clean`)return e.cleanedUp;if(i===`straighten`)return e.straightened;if(i===`clipboard.paste`)return e.format(`pasted`,{what:e.atoms(o.added.length)});if(i===`clipboard.drop`||i===`api.add`)return e.format(`added`,{what:e.atoms(o.added.length)});if(i===`selection.rotate`||i.endsWith(`.rotate180`))return e.format(`rotated`,{what:e.atoms(o.moved)});if(i.endsWith(`.flip`))return e.format(`flipped`,{what:e.atoms(o.moved)});if(i===`atom.join`)return e.format(`joined`,{what:e.atoms(o.removed+1)});if(i===`abbreviation.expand`||i===`abbreviation.contract`||i===`abbreviation.remove`||i===`abbreviation.contract-label`)return Rs(e,i);if(o.added.length===1&&o.addedBonds<=1&&!o.removed)return e.format(`added`,{what:o.added[0].label});if(o.added.length||o.addedBonds)return e.format(`added`,{what:Vs(e,o.added.length,o.addedBonds)});if(o.removed||o.removedBonds)return e.format(`deleted`,{what:Vs(e,o.removed,o.removedBonds)});if(o.bonds.length&&!o.atoms.length){let t=new Set(o.bonds.map(e=>r.bonds[e].type));return o.bonds.length===1?e.format(`bondNow`,{number:o.bonds[0]+1,type:Ns(e,r.bonds[o.bonds[0]].type)}):t.size===1?e.format(`bondsNow`,{bonds:e.bonds(o.bonds.length),type:Ns(e,r.bonds[o.bonds[0]].type)}):e.format(`changed`,{what:e.bonds(o.bonds.length)})}if(o.atoms.length&&!o.bonds.length){let t=new Set(o.atoms.map(e=>r.atoms[e].label));return o.atoms.length===1?e.format(`atomNow`,{number:o.atoms[0]+1,label:r.atoms[o.atoms[0]].label}):t.size===1?e.format(`atomsNow`,{atoms:e.atoms(o.atoms.length),label:r.atoms[o.atoms[0]].label}):e.format(`changed`,{what:e.atoms(o.atoms.length)})}return o.atoms.length||o.bonds.length?e.format(`changed`,{what:Vs(e,o.atoms.length,o.bonds.length)}):o.moved?e.format(`moved`,{what:e.atoms(o.moved)}):Rs(e,i)}function Ws(e,t){let n=[...t.centres?[e.count(`centreCount`,t.centres)]:[],...t.axes?[e.count(`axisCount`,t.axes)]:[],...t.kind===`stereo-removed`&&t.marks&&!t.centres&&!t.axes?[e.count(`stereoBondCount`,t.marks)]:[]],r=n.length>1?e.format(`both`,{first:n[0],second:n[1]}):n[0]??``;return e.format(t.kind===`stereo-removed`?`stereoRemoved`:`wedgesPlaced`,{what:r})}function Gs(e,t){if(!t.atoms)return e.descriptionEmpty;let n=t.smiles.length>200?`${t.smiles.slice(0,200)}…`:t.smiles,r=Vs(e,t.atoms,t.bonds);return e.format(`description`,{counts:t.formula?e.format(`withFormula`,{formula:t.formula,counts:r}):r,smiles:n})}var Ks={up:{x:0,y:-1},down:{x:0,y:1},left:{x:-1,y:0},right:{x:1,y:0}};function qs(e,t,n){let[r,i]=[t.x-e.x,t.y-e.y],a=Math.hypot(r,i);if(a<1e-9)return null;let o=Ks[n],s=Math.max(-1,Math.min(1,(r*o.x+i*o.y)/a));return Math.acos(s)*180/Math.PI}var Js=1e-6;function Ys(e){let t=e.atoms.map(()=>[]);return e.bonds.forEach((e,n)=>{t[e.a]?.push(n),e.b!==e.a&&t[e.b]?.push(n)}),t}var Xs=(e,t,n)=>e.bonds[t].a===n?e.bonds[t].b:e.bonds[t].a;function Zs(e,t,n,r){let i=-1,a=1/0;for(let o of r){let r=qs(e.atoms[t],e.atoms[Xs(e,o,t)],n);r===null||r>75.000001||(r<a-Js||Math.abs(r-a)<=Js&&o<i)&&(i=o,a=r)}return i}var Qs=(e,t)=>{let{a:n,b:r}=e.bonds[t];return{x:(e.atoms[n].x+e.atoms[r].x)/2,y:(e.atoms[n].y+e.atoms[r].y)/2}};function $s(e,t,n){let r=Qs(e,t),{a:i,b:a}=e.bonds[t],[o,s]=[qs(r,e.atoms[i],n),qs(r,e.atoms[a],n)],c=e=>e!==null&&e<89.999999;return c(o)&&(!c(s)||o<=s)?i:c(s)?a:-1}function ec(e,t,n){if(t.kind===`atom`){if(!e.atoms[t.index])return null;let r=Zs(e,t.index,n,Ys(e)[t.index]??[]);return r>=0?{kind:`bond`,index:r}:null}if(!e.bonds[t.index])return null;let r=$s(e,t.index,n);return r>=0?{kind:`atom`,index:r}:null}function tc(e,t,n){let r=Ys(e),i,a=0;if(t.kind===`bond`){if(!e.bonds[t.index]||(i=$s(e,t.index,n),i<0))return null;a=1}else{if(!e.atoms[t.index])return null;i=t.index}let o=new Set([i]);for(;a<3;){let t=Zs(e,i,n,(r[i]??[]).filter(t=>!o.has(Xs(e,t,i))));if(t<0)break;i=Xs(e,t,i),o.add(i),a+=1}return a?{kind:`atom`,index:i}:null}function nc(e){let t=e.atoms.map((e,t)=>t),n=e=>{for(;t[e]!==e;)t[e]=t[t[e]],e=t[e];return e};for(let r of e.bonds)e.atoms[r.a]&&e.atoms[r.b]&&(t[n(r.a)]=n(r.b));let r=new Map;return e.atoms.map((e,t)=>{let i=n(t);return r.has(i)||r.set(i,r.size),r.get(i)})}function rc(e,t,n){let r=t.kind===`atom`?e.atoms[t.index]:e.bonds[t.index]?Qs(e,t.index):void 0;if(!r)return null;let i=nc(e),a=i[t.kind===`atom`?t.index:e.bonds[t.index].a],o=-1,s=1/0;return e.atoms.forEach((e,t)=>{if(i[t]===a)return;let c=qs(r,e,n);if(c===null||c>75.000001)return;let l=Math.hypot(e.x-r.x,e.y-r.y);l<s-1e-9&&(o=t,s=l)}),o>=0?{kind:`atom`,index:o}:null}function ic(e,t,n,r){let i=new Set(e.atoms.filter(e=>t.atoms[e]));for(let n of e.bonds)t.bonds[n]&&(i.add(t.bonds[n].a),i.add(t.bonds[n].b));if(n?.kind===`atom`&&i.has(n.index))return n.index;let a=[...i].sort((e,t)=>t-e);return a.find(r)??a[0]??-1}function ac(e,t,n){let r=Ys(e),i=t.kind===`atom`?[t.index]:e.bonds[t.index]?[e.bonds[t.index].a,e.bonds[t.index].b]:[],a=i.map(e=>n[e]??-1).filter(e=>e>=0);if(t.kind===`bond`&&a.length)return Math.min(...a);let o=new Set(i),s=i;for(;s.length;){let t=[];for(let n of s)for(let i of r[n]??[]){let r=Xs(e,i,n);o.has(r)||(o.add(r),t.push(r))}let i=t.map(e=>n[e]??-1).filter(e=>e>=0);if(i.length)return Math.min(...i);s=t}return-1}var oc={anyBond:{viewBox:`0 0 48 48`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M15.9424 32.0782C16.7178 32.8648 16.7087 34.1311 15.9221 34.9065L12.2761 38.5005C11.4895 39.2759 10.2232 39.2668 9.44774 38.4802C8.67232 37.6936 8.68141 36.4273 9.46805 35.6518L13.114 32.0578C13.9007 31.2824 15.167 31.2915 15.9424 32.0782Z" fill="currentColor"/><path fill-rule="evenodd" clip-rule="evenodd" d="M23.4574 24.5383C24.2326 25.3251 24.2232 26.5914 23.4363 27.3666L19.7903 30.9586C19.0035 31.7338 17.7372 31.7244 16.962 30.9375C16.1868 30.1507 16.1962 28.8844 16.9831 28.1092L20.6291 24.5172C21.416 23.742 22.6822 23.7514 23.4574 24.5383Z" fill="currentColor"/><path fill-rule="evenodd" clip-rule="evenodd" d="M31.0131 17.0363C31.7883 17.8232 31.7789 19.0895 30.992 19.8647L27.346 23.4567C26.5592 24.2319 25.2929 24.2224 24.5177 23.4356C23.7425 22.6487 23.7519 21.3824 24.5388 20.6072L28.1848 17.0152C28.9716 16.24 30.2379 16.2495 31.0131 17.0363Z" fill="currentColor"/><path fill-rule="evenodd" clip-rule="evenodd" d="M38.5528 9.52005C39.3282 10.3067 39.3191 11.573 38.5324 12.3484L34.8864 15.9424C34.0998 16.7178 32.8335 16.7087 32.0581 15.9221C31.2827 15.1355 31.2918 13.8692 32.0784 13.0937L35.7244 9.49974C36.511 8.72432 37.7773 8.73341 38.5528 9.52005Z" fill="currentColor"/>`},singleOrDouble:{viewBox:`0 0 48 48`,body:`<path d="M19.0775 31.6002L15.5335 35.0882C15.3498 35.2871 15.1275 35.4466 14.8803 35.5571C14.633 35.6675 14.3659 35.7266 14.0952 35.7308C13.8244 35.735 13.5556 35.6841 13.3051 35.5813C13.0546 35.4785 12.8276 35.3258 12.6378 35.1327C12.448 34.9395 12.2994 34.7098 12.2011 34.4575C12.1027 34.2052 12.0566 33.9356 12.0655 33.6649C12.0745 33.3943 12.1383 33.1283 12.2531 32.883C12.368 32.6378 12.5314 32.4184 12.7335 32.2382L16.2775 28.7442C16.4613 28.5453 16.6835 28.3857 16.9308 28.2753C17.178 28.1648 17.4451 28.1057 17.7159 28.1016C17.9866 28.0974 18.2554 28.1483 18.5059 28.2511C18.7565 28.3539 18.9835 28.5065 19.1733 28.6997C19.363 28.8928 19.5116 29.1225 19.61 29.3748C19.7083 29.6271 19.7544 29.8968 19.7455 30.1674C19.7365 30.4381 19.6727 30.7041 19.5579 30.9493C19.4431 31.1946 19.2797 31.414 19.0775 31.5942V31.6002Z" fill="currentColor"/><path d="M30.2007 19.8321C29.8044 19.8312 29.4172 19.7126 29.0883 19.4913C28.7595 19.27 28.5039 18.9559 28.3538 18.589C28.2038 18.2221 28.1662 17.8189 28.2458 17.4306C28.3254 17.0423 28.5185 16.6864 28.8007 16.4081L32.3447 12.9141C32.5285 12.7152 32.7507 12.5557 32.998 12.4452C33.2452 12.3347 33.5123 12.2756 33.7831 12.2715C34.0538 12.2673 34.3226 12.3182 34.5731 12.421C34.8237 12.5238 35.0507 12.6764 35.2405 12.8696C35.4302 13.0628 35.5788 13.2925 35.6772 13.5448C35.7756 13.7971 35.8217 14.0667 35.8127 14.3374C35.8037 14.608 35.7399 14.874 35.6251 15.1193C35.5103 15.3645 35.3469 15.5839 35.1447 15.7641L31.6007 19.2581C31.2272 19.6254 30.7246 19.8315 30.2007 19.8321Z" fill="currentColor"/><path d="M22.1598 27.7479C21.763 27.7475 21.3753 27.629 21.046 27.4076C20.7168 27.1862 20.4608 26.8718 20.3107 26.5045C20.1606 26.1372 20.1232 25.7335 20.2031 25.3449C20.2831 24.9562 20.4769 24.6002 20.7598 24.3219L24.3038 20.8279C24.6872 20.4861 25.1874 20.3049 25.7008 20.3219C26.2142 20.3389 26.7014 20.5527 27.0614 20.9192C27.4214 21.2856 27.6266 21.7765 27.6345 22.2901C27.6424 22.8037 27.4523 23.3006 27.1038 23.6779L23.5598 27.1719C23.1866 27.54 22.6839 27.7468 22.1598 27.7479Z" fill="currentColor"/><path d="M22.3999 16.352C22.0036 16.3512 21.6164 16.2325 21.2876 16.0112C20.9587 15.7899 20.7031 15.4759 20.5531 15.109C20.4031 14.7421 20.3654 14.3389 20.445 13.9506C20.5246 13.5623 20.7177 13.2064 20.9999 12.928L24.5439 9.43403C24.727 9.23294 24.9492 9.07131 25.1969 8.959C25.4446 8.8467 25.7126 8.78609 25.9845 8.78088C26.2564 8.77568 26.5265 8.82597 26.7783 8.92871C27.0301 9.03144 27.2583 9.18445 27.449 9.37837C27.6396 9.5723 27.7887 9.80308 27.8872 10.0566C27.9856 10.3101 28.0313 10.581 28.0215 10.8528C28.0116 11.1246 27.9465 11.3915 27.83 11.6372C27.7135 11.883 27.5481 12.1024 27.3439 12.282L23.7999 15.776C23.4268 16.1441 22.9241 16.3509 22.3999 16.352Z" fill="currentColor"/><path d="M14.3716 24.2659C13.9752 24.265 13.5881 24.1464 13.2592 23.9251C12.9304 23.7037 12.6748 23.3897 12.5247 23.0228C12.3747 22.6559 12.3371 22.2527 12.4167 21.8644C12.4962 21.4761 12.6894 21.1202 12.9716 20.8419L16.5116 17.3479C16.6989 17.1635 16.9206 17.0178 17.1642 16.9191C17.4078 16.8205 17.6684 16.7707 17.9312 16.7728C18.194 16.7748 18.4538 16.8286 18.6958 16.931C18.9378 17.0335 19.1572 17.1826 19.3416 17.3699C19.526 17.5571 19.6717 17.7789 19.7704 18.0225C19.869 18.266 19.9188 18.5266 19.9167 18.7894C19.9147 19.0522 19.8609 19.312 19.7585 19.554C19.656 19.796 19.5069 20.0155 19.3196 20.1999L15.7756 23.6899C15.4015 24.0589 14.8971 24.2658 14.3716 24.2659Z" fill="currentColor"/><path d="M29.9458 31.2283C29.5495 31.2274 29.1623 31.1088 28.8335 30.8874C28.5046 30.6661 28.249 30.3521 28.099 29.9852C27.9489 29.6183 27.9113 29.2151 27.9909 28.8268C28.0705 28.4385 28.2636 28.0826 28.5458 27.8043L32.0918 24.3103C32.2749 24.1092 32.4971 23.9475 32.7448 23.8352C32.9925 23.7229 33.2605 23.6623 33.5324 23.6571C33.8043 23.6519 34.0744 23.7022 34.3262 23.8049C34.578 23.9077 34.8062 24.0607 34.9969 24.2546C35.1875 24.4485 35.3366 24.6793 35.4351 24.9328C35.5335 25.1863 35.5792 25.4573 35.5694 25.729C35.5595 26.0008 35.4944 26.2677 35.3779 26.5135C35.2614 26.7592 35.096 26.9786 34.8918 27.1582L31.3458 30.6522C30.9727 31.0203 30.47 31.2271 29.9458 31.2283Z" fill="currentColor"/><path d="M21.9165 39.1421C21.5202 39.1412 21.133 39.0226 20.8042 38.8013C20.4753 38.5799 20.2197 38.2659 20.0697 37.899C19.9197 37.5321 19.882 37.1289 19.9616 36.7406C20.0412 36.3523 20.2343 35.9964 20.5165 35.7181L24.0625 32.2241C24.2456 32.023 24.4678 31.8613 24.7155 31.749C24.9632 31.6367 25.2312 31.5761 25.5031 31.5709C25.775 31.5657 26.0451 31.616 26.2969 31.7187C26.5487 31.8215 26.7769 31.9745 26.9676 32.1684C27.1582 32.3623 27.3073 32.5931 27.4058 32.8466C27.5042 33.1001 27.5499 33.3711 27.5401 33.6428C27.5302 33.9146 27.4651 34.1815 27.3486 34.4273C27.2321 34.673 27.0667 34.8924 26.8625 35.0721L23.3205 38.5661C22.9464 38.9351 22.442 39.142 21.9165 39.1421Z" fill="currentColor"/>`},singleOrAromatic:{viewBox:`0 0 48 48`,body:`<path d="M15.8863 40.04C15.4899 40.0392 15.1027 39.9205 14.7739 39.6992C14.4451 39.4779 14.1894 39.1639 14.0394 38.797C13.8894 38.4301 13.8518 38.0269 13.9313 37.6386C14.0109 37.2503 14.2041 36.8944 14.4863 36.616L34.0863 17.292C34.27 17.0931 34.4923 16.9336 34.7395 16.8231C34.9867 16.7127 35.2538 16.6536 35.5246 16.6494C35.7954 16.6452 36.0642 16.6961 36.3147 16.7989C36.5652 16.9017 36.7922 17.0544 36.982 17.2475C37.1718 17.4407 37.3204 17.6704 37.4187 17.9227C37.5171 18.175 37.5632 18.4446 37.5543 18.7153C37.5453 18.9859 37.4815 19.2519 37.3667 19.4972C37.2518 19.7424 37.0884 19.9618 36.8863 20.142L17.2863 39.464C16.9131 39.8321 16.4104 40.0389 15.8863 40.04Z" fill="currentColor"/><path d="M17.5844 27.342L14.0404 30.836C13.8566 31.0349 13.6344 31.1944 13.3871 31.3049C13.1399 31.4153 12.8728 31.4744 12.602 31.4786C12.3313 31.4828 12.0625 31.4319 11.8119 31.3291C11.5614 31.2263 11.3344 31.0736 11.1446 30.8805C10.9548 30.6873 10.8063 30.4576 10.7079 30.2053C10.6095 29.953 10.5634 29.6834 10.5724 29.4127C10.5813 29.1421 10.6452 28.8761 10.76 28.6308C10.8748 28.3856 11.0382 28.1662 11.2404 27.986L14.7844 24.492C14.9681 24.2931 15.1903 24.1335 15.4376 24.0231C15.6848 23.9126 15.9519 23.8535 16.2227 23.8494C16.4935 23.8452 16.7622 23.8961 17.0128 23.9989C17.2633 24.1017 17.4903 24.2543 17.6801 24.4475C17.8699 24.6406 18.0185 24.8703 18.1168 25.1226C18.2152 25.3749 18.2613 25.6446 18.2523 25.9152C18.2434 26.1859 18.1796 26.4519 18.0647 26.6971C17.9499 26.9424 17.7865 27.1618 17.5844 27.342Z" fill="currentColor"/><path d="M24.441 19.6001C24.0446 19.5993 23.6574 19.4806 23.3286 19.2593C22.9997 19.038 22.7441 18.724 22.5941 18.3571C22.4441 17.9902 22.4065 17.587 22.486 17.1987C22.5656 16.8104 22.7587 16.4545 23.041 16.1761L26.581 12.6901C26.764 12.489 26.9862 12.3274 27.2339 12.2151C27.4816 12.1028 27.7496 12.0422 28.0215 12.037C28.2934 12.0318 28.5635 12.0821 28.8153 12.1848C29.0671 12.2875 29.2953 12.4405 29.486 12.6345C29.6766 12.8284 29.8257 13.0592 29.9242 13.3127C30.0226 13.5662 30.0683 13.8371 30.0585 14.1089C30.0487 14.3807 29.9835 14.6476 29.867 14.8933C29.7505 15.1391 29.5851 15.3585 29.381 15.5381L25.845 19.0321C25.4695 19.3982 24.9653 19.6022 24.441 19.6001Z" fill="currentColor"/><path d="M32.1133 11.96C33.2179 11.96 34.1133 11.0645 34.1133 9.95996C34.1133 8.85539 33.2179 7.95996 32.1133 7.95996C31.0087 7.95996 30.1133 8.85539 30.1133 9.95996C30.1133 11.0645 31.0087 11.96 32.1133 11.96Z" fill="currentColor"/><path d="M20.3105 23.762C21.4151 23.762 22.3105 22.8665 22.3105 21.762C22.3105 20.6574 21.4151 19.762 20.3105 19.762C19.206 19.762 18.3105 20.6574 18.3105 21.762C18.3105 22.8665 19.206 23.762 20.3105 23.762Z" fill="currentColor"/>`},doubleOrAromatic:{viewBox:`0 0 48 48`,body:`<path d="M17.5717 27.342L14.0277 30.836C13.8439 31.0349 13.6217 31.1944 13.3744 31.3049C13.1272 31.4153 12.8601 31.4744 12.5893 31.4786C12.3186 31.4828 12.0498 31.4319 11.7993 31.3291C11.5487 31.2263 11.3217 31.0736 11.1319 30.8805C10.9421 30.6873 10.7936 30.4576 10.6952 30.2053C10.5968 29.953 10.5507 29.6834 10.5597 29.4127C10.5686 29.1421 10.6325 28.8761 10.7473 28.6308C10.8621 28.3856 11.0255 28.1662 11.2277 27.986L14.7717 24.492C14.9554 24.2931 15.1777 24.1335 15.4249 24.0231C15.6721 23.9126 15.9392 23.8535 16.21 23.8494C16.4808 23.8452 16.7496 23.8961 17.0001 23.9989C17.2506 24.1017 17.4776 24.2543 17.6674 24.4475C17.8572 24.6406 18.0058 24.8703 18.1041 25.1226C18.2025 25.3749 18.2486 25.6446 18.2396 25.9152C18.2307 26.1859 18.1669 26.4519 18.052 26.6971C17.9372 26.9424 17.7738 27.1618 17.5717 27.342Z" fill="currentColor"/><path d="M24.4273 19.5998C24.0309 19.599 23.6437 19.4804 23.3149 19.259C22.9861 19.0377 22.7304 18.7237 22.5804 18.3568C22.4304 17.9899 22.3928 17.5867 22.4723 17.1984C22.5519 16.8101 22.7451 16.4542 23.0273 16.1758L26.5713 12.6818C26.9551 12.3441 27.4538 12.1662 27.9647 12.1847C28.4757 12.2031 28.9602 12.4166 29.3186 12.7812C29.6771 13.1458 29.8823 13.6338 29.8921 14.145C29.9019 14.6562 29.7155 15.1518 29.3713 15.5298L25.8313 19.0318C25.4558 19.3979 24.9517 19.6019 24.4273 19.5998Z" fill="currentColor"/><path d="M32.1016 11.96C33.2061 11.96 34.1016 11.0645 34.1016 9.95996C34.1016 8.85539 33.2061 7.95996 32.1016 7.95996C30.997 7.95996 30.1016 8.85539 30.1016 9.95996C30.1016 11.0645 30.997 11.96 32.1016 11.96Z" fill="currentColor"/><path d="M20.2979 23.762C21.4024 23.762 22.2979 22.8665 22.2979 21.762C22.2979 20.6574 21.4024 19.762 20.2979 19.762C19.1933 19.762 18.2979 20.6574 18.2979 21.762C18.2979 22.8665 19.1933 23.762 20.2979 23.762Z" fill="currentColor"/><path d="M20.8881 36.0361L17.3441 39.5301C17.1603 39.729 16.9381 39.8885 16.6908 39.999C16.4436 40.1094 16.1765 40.1685 15.9057 40.1727C15.635 40.1768 15.3662 40.126 15.1157 40.0232C14.8651 39.9204 14.6381 39.7677 14.4483 39.5746C14.2585 39.3814 14.11 39.1517 14.0116 38.8994C13.9132 38.6471 13.8671 38.3775 13.8761 38.1068C13.885 37.8362 13.9489 37.5701 14.0637 37.3249C14.1785 37.0797 14.3419 36.8603 14.5441 36.6801L18.0801 33.1861C18.2638 32.9872 18.4861 32.8276 18.7333 32.7172C18.9805 32.6067 19.2476 32.5476 19.5184 32.5434C19.7892 32.5393 20.058 32.5902 20.3085 32.693C20.559 32.7958 20.786 32.9484 20.9758 33.1416C21.1656 33.3347 21.3142 33.5644 21.4125 33.8167C21.5109 34.069 21.557 34.3387 21.548 34.6093C21.5391 34.88 21.4753 35.146 21.3604 35.3912C21.2456 35.6365 21.0822 35.8559 20.8801 36.0361H20.8881Z" fill="currentColor"/><path d="M27.7437 28.302C27.3473 28.3011 26.9601 28.1825 26.6313 27.9612C26.3025 27.7398 26.0468 27.4258 25.8968 27.0589C25.7468 26.692 25.7092 26.2888 25.7888 25.9005C25.8683 25.5122 26.0615 25.1563 26.3437 24.878L29.8877 21.384C30.0708 21.1829 30.293 21.0213 30.5406 20.909C30.7883 20.7967 31.0563 20.736 31.3282 20.7308C31.6001 20.7256 31.8702 20.7759 32.122 20.8787C32.3738 20.9814 32.6021 21.1344 32.7927 21.3283C32.9834 21.5223 33.1325 21.753 33.2309 22.0065C33.3294 22.2601 33.3751 22.531 33.3652 22.8028C33.3554 23.0745 33.2902 23.3415 33.1737 23.5872C33.0573 23.8329 32.8919 24.0523 32.6877 24.232L29.1437 27.726C28.7705 28.094 28.2678 28.3008 27.7437 28.302Z" fill="currentColor"/><path d="M35.418 20.6541C36.5225 20.6541 37.418 19.7586 37.418 18.6541C37.418 17.5495 36.5225 16.6541 35.418 16.6541C34.3134 16.6541 33.418 17.5495 33.418 18.6541C33.418 19.7586 34.3134 20.6541 35.418 20.6541Z" fill="currentColor"/><path d="M23.6143 32.4561C24.7188 32.4561 25.6143 31.5606 25.6143 30.4561C25.6143 29.3515 24.7188 28.4561 23.6143 28.4561C22.5097 28.4561 21.6143 29.3515 21.6143 30.4561C21.6143 31.5606 22.5097 32.4561 23.6143 32.4561Z" fill="currentColor"/>`},single:{viewBox:`0 0 48 48`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M35.2266 12.9341C36.002 13.7208 35.9929 14.9871 35.2062 15.7625L15.6022 35.0865C14.8156 35.8619 13.5493 35.8528 12.7739 35.0661C11.9985 34.2795 12.0076 33.0132 12.7942 32.2378L32.3982 12.9138C33.1849 12.1384 34.4512 12.1475 35.2266 12.9341Z" fill="currentColor"/>`},double:{viewBox:`0 0 48 48`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M36.9162 17.3142C37.6916 18.1009 37.6823 19.3672 36.8956 20.1425L17.2896 39.4645C16.5029 40.2398 15.2366 40.2306 14.4612 39.4439C13.6859 38.6571 13.6952 37.3909 14.4819 36.6155L34.0879 17.2935C34.8746 16.5182 36.1409 16.5274 36.9162 17.3142Z" fill="currentColor"/> <path fill-rule="evenodd" clip-rule="evenodd" d="M33.5383 8.5561C34.3136 9.34283 34.3044 10.6091 33.5177 11.3845L13.9117 30.7065C13.1249 31.4818 11.8586 31.4725 11.0833 30.6858C10.308 29.8991 10.3172 28.6328 11.104 27.8575L30.71 8.53547C31.4967 7.76014 32.763 7.76938 33.5383 8.5561Z" fill="currentColor"/>`},triple:{viewBox:`0 0 48 48`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M35.1584 12.9282C35.9338 13.7149 35.9245 14.9812 35.1378 15.7565L15.5318 35.0785C14.7451 35.8539 13.4788 35.8446 12.7034 35.0579C11.9281 34.2712 11.9373 33.0049 12.7241 32.2295L32.3301 12.9075C33.1168 12.1322 34.3831 12.1414 35.1584 12.9282Z" fill="currentColor"/> <path fill-rule="evenodd" clip-rule="evenodd" d="M31.3059 5.49018C32.0812 6.27691 32.072 7.54321 31.2852 8.31854L11.6792 27.6405C10.8925 28.4159 9.62623 28.4066 8.8509 27.6199C8.07557 26.8332 8.08481 25.5669 8.87153 24.7915L28.4775 5.46955C29.2643 4.69422 30.5306 4.70346 31.3059 5.49018Z" fill="currentColor"/> <path fill-rule="evenodd" clip-rule="evenodd" d="M39.0129 20.3662C39.7882 21.1529 39.779 22.4192 38.9923 23.1945L19.3863 42.5165C18.5996 43.2918 17.3333 43.2826 16.5579 42.4959C15.7826 41.7092 15.7918 40.4429 16.5786 39.6675L36.1846 20.3455C36.9713 19.5702 38.2376 19.5794 39.0129 20.3662Z" fill="currentColor"/>`},aromatic:{viewBox:`0 0 48 48`,body:`<path d="M15.8863 40.0398C15.4899 40.0389 15.1027 39.9203 14.7739 39.699C14.4451 39.4777 14.1894 39.1636 14.0394 38.7967C13.8894 38.4299 13.8518 38.0267 13.9313 37.6384C14.0109 37.25 14.2041 36.8941 14.4863 36.6158L34.0863 17.2918C34.4697 16.95 34.9699 16.7688 35.4833 16.7858C35.9967 16.8027 36.4839 17.0166 36.8439 17.383C37.2039 17.7494 37.4091 18.2403 37.417 18.7539C37.4249 19.2675 37.2348 19.7645 36.8863 20.1418L17.2863 39.4638C16.9131 39.8318 16.4104 40.0387 15.8863 40.0398Z" fill="currentColor"/> <path d="M17.4555 27.2141L13.9115 30.708C13.7277 30.907 13.5055 31.0665 13.2582 31.1769C13.011 31.2874 12.7439 31.3465 12.4731 31.3507C12.2023 31.3548 11.9336 31.304 11.683 31.2011C11.4325 31.0983 11.2055 30.9457 11.0157 30.7525C10.8259 30.5594 10.6774 30.3297 10.579 30.0774C10.4806 29.8251 10.4345 29.5554 10.4435 29.2848C10.4524 29.0141 10.5162 28.7481 10.6311 28.5029C10.7459 28.2576 10.9093 28.0383 11.1115 27.858L14.6555 24.364C14.8392 24.1651 15.0614 24.0056 15.3087 23.8952C15.5559 23.7847 15.823 23.7256 16.0938 23.7214C16.3646 23.7173 16.6333 23.7681 16.8839 23.8709C17.1344 23.9738 17.3614 24.1264 17.5512 24.3195C17.741 24.5127 17.8895 24.7424 17.9879 24.9947C18.0863 25.247 18.1324 25.5167 18.1234 25.7873C18.1145 26.058 18.0507 26.324 17.9358 26.5692C17.821 26.8145 17.6576 27.0338 17.4555 27.2141Z" fill="currentColor"/> <path d="M28.5679 15.4541C28.1715 15.4532 27.7844 15.3346 27.4555 15.1133C27.1267 14.892 26.871 14.5779 26.721 14.211C26.571 13.8441 26.5334 13.4409 26.613 13.0526C26.6925 12.6643 26.8857 12.3084 27.1679 12.0301L30.7139 8.53609C31.0978 8.19834 31.5964 8.02043 32.1074 8.03891C32.6183 8.05739 33.1028 8.27086 33.4613 8.63545C33.8197 9.00005 34.0249 9.48809 34.0347 9.99929C34.0445 10.5105 33.8581 11.006 33.5139 11.3841L29.9719 14.8781C29.5978 15.2471 29.0934 15.454 28.5679 15.4541Z" fill="currentColor"/> <path d="M20.5386 23.3681C20.1422 23.3673 19.7551 23.2487 19.4262 23.0273C19.0974 22.806 18.8418 22.492 18.6917 22.1251C18.5417 21.7582 18.5041 21.355 18.5837 20.9667C18.6632 20.5784 18.8564 20.2225 19.1386 19.9441L22.6806 16.4501C22.8637 16.2491 23.0859 16.0874 23.3336 15.9751C23.5813 15.8628 23.8492 15.8022 24.1211 15.797C24.3931 15.7918 24.6632 15.8421 24.915 15.9448C25.1668 16.0476 25.395 16.2006 25.5856 16.3945C25.7763 16.5884 25.9254 16.8192 26.0238 17.0727C26.1223 17.3262 26.168 17.5971 26.1582 17.8689C26.1483 18.1407 26.0832 18.4076 25.9667 18.6534C25.8502 18.8991 25.6848 19.1185 25.4806 19.2981L21.9426 22.8001C21.5672 23.1662 21.063 23.3702 20.5386 23.3681Z" fill="currentColor"/>`},wedge:{viewBox:`0 0 24 24`,body:`<path d="M18.7631 10.8294L4.70145 18.9724C4.61237 19.0237 4.50889 19.0442 4.40699 19.0308C4.30509 19.0174 4.21043 18.9708 4.13764 18.8983C4.06484 18.8257 4.01795 18.7312 4.0042 18.6294C3.99046 18.5275 4.01063 18.424 4.06159 18.3347L12.1679 4.23411C12.2037 4.17195 12.2534 4.11887 12.3131 4.07897C12.3727 4.03908 12.4408 4.01344 12.5119 4.00403C12.5831 3.99463 12.6554 4.00172 12.7234 4.02474C12.7914 4.04777 12.8531 4.08612 12.9039 4.13683L18.8604 10.0901C18.912 10.1408 18.951 10.2027 18.9746 10.2711C18.9981 10.3394 19.0054 10.4123 18.996 10.484C18.9866 10.5556 18.9606 10.6241 18.9202 10.6841C18.8798 10.744 18.826 10.7938 18.7631 10.8294Z" fill="currentColor"/>`},hash:{viewBox:`0 0 24 24`,body:`<path d="M17.3106 11.9541L18.9033 11.0172C19.1852 10.8514 19.2346 10.464 19.0033 10.2327L12.9171 4.14649C12.6859 3.91521 12.2985 3.96461 12.1326 4.24654L11.1957 5.83928L17.3106 11.9541Z" fill="currentColor"/> <path d="M9.34387 8.98742L10.3914 7.20656L15.9433 12.7584L14.1624 13.806L9.34387 8.98742Z" fill="currentColor"/> <path d="M7.49201 12.1356L8.53958 10.3547L12.7951 14.6102L11.0143 15.6578L7.49201 12.1356Z" fill="currentColor"/> <path d="M5.64016 15.2837L6.68773 13.5029L9.64697 16.4621L7.86611 17.5097L5.64016 15.2837Z" fill="currentColor"/> <path d="M6.49882 18.314L5.28615 19.0273C4.53026 19.4719 3.6779 18.6196 4.12254 17.8637L4.83588 16.651L6.49882 18.314Z" fill="currentColor"/>`},wavy:{viewBox:`0 0 24 24`,body:`<path d="M5 18.5L7 15L11 16.5L9.5 11.5L14.5 14L11.5 8L19 11.5L13 5" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>`},crossed:{viewBox:`0 0 48 48`,body:`<path d="M39.5635 18.7178L25.4655 25.5178L18.2155 39.3838C18.0462 39.7078 17.7914 39.9793 17.4786 40.1686C17.1658 40.358 16.8071 40.458 16.4415 40.4578C16.1189 40.4586 15.8009 40.3803 15.5155 40.2298C15.0457 39.984 14.6926 39.5618 14.534 39.0559C14.3753 38.55 14.4241 38.0018 14.6695 37.5318L19.4315 28.4218L10.1715 32.8858C9.93373 33.0125 9.67279 33.09 9.40434 33.1134C9.13589 33.1368 8.86547 33.1058 8.60932 33.0221C8.35318 32.9384 8.11659 32.8038 7.91376 32.6264C7.71092 32.449 7.54603 32.2324 7.42899 31.9897C7.31194 31.747 7.24515 31.4831 7.23263 31.2139C7.22011 30.9447 7.26211 30.6758 7.35612 30.4233C7.45013 30.1707 7.5942 29.9398 7.77968 29.7443C7.96517 29.5488 8.18824 29.3929 8.43551 29.2858L22.5335 22.4858L29.7835 8.61576C30.0291 8.14553 30.4514 7.79212 30.9576 7.63328C31.4638 7.47443 32.0123 7.52317 32.4825 7.76876C32.9527 8.01435 33.3062 8.43668 33.465 8.94284C33.6238 9.449 33.5751 9.99753 33.3295 10.4678L28.5675 19.5778L37.8275 15.1138C38.0653 14.987 38.3262 14.9095 38.5947 14.8861C38.8631 14.8627 39.1336 14.8937 39.3897 14.9774C39.6458 15.0611 39.8824 15.1957 40.0853 15.3731C40.2881 15.5505 40.453 15.7671 40.57 16.0098C40.6871 16.2525 40.7539 16.5164 40.7664 16.7856C40.7789 17.0548 40.7369 17.3237 40.6429 17.5762C40.5489 17.8288 40.4048 18.0597 40.2193 18.2552C40.0339 18.4507 39.8108 18.6066 39.5635 18.7138V18.7178Z" fill="currentColor"/>`},erase:{viewBox:`0 0 48 48`,body:`<path d="M40.2844 18.4577L30.3704 8.54566C29.9954 8.17071 29.4868 7.96008 28.9564 7.96008C28.4261 7.96008 27.9175 8.17071 27.5424 8.54566L7.71643 28.3717C7.34149 28.7467 7.13086 29.2553 7.13086 29.7857C7.13086 30.316 7.34149 30.8246 7.71643 31.1997L14.8004 38.2817C15.3576 38.8392 16.0192 39.2815 16.7474 39.5833C17.4756 39.8851 18.2562 40.0404 19.0444 40.0404C19.8327 40.0404 20.6132 39.8851 21.3414 39.5833C22.0696 39.2815 22.7312 38.8392 23.2884 38.2817L40.2884 21.2817C40.6634 20.9066 40.874 20.398 40.874 19.8677C40.874 19.3373 40.6634 18.8287 40.2884 18.4537L40.2844 18.4577ZM20.4604 35.4577C20.2747 35.6439 20.0539 35.7917 19.811 35.8925C19.568 35.9933 19.3075 36.0452 19.0444 36.0452C18.7814 36.0452 18.5209 35.9933 18.2779 35.8925C18.0349 35.7917 17.8142 35.6439 17.6284 35.4577L11.9584 29.7857L19.7824 21.9617L26.8664 29.0477L20.4604 35.4577Z" fill="currentColor"/>`},clear:{viewBox:`0 0 24 24`,body:`<path d="M18.6554 7.52861L15.0194 3.89261C14.8319 3.70506 14.5776 3.59967 14.3124 3.59961H6.39844C6.13322 3.59961 5.87887 3.70497 5.69133 3.8925C5.50379 4.08004 5.39844 4.33439 5.39844 4.59961V19.1477C5.39844 19.279 5.4243 19.409 5.47456 19.5303C5.52481 19.6517 5.59847 19.7619 5.69133 19.8548C5.78419 19.9476 5.89443 20.0213 6.01575 20.0715C6.13708 20.1218 6.26712 20.1477 6.39844 20.1477H17.9484C18.0798 20.1477 18.2098 20.1218 18.3311 20.0715C18.4524 20.0213 18.5627 19.9476 18.6555 19.8548C18.7484 19.7619 18.8221 19.6517 18.8723 19.5303C18.9226 19.409 18.9484 19.279 18.9484 19.1477V8.23561C18.9484 7.97041 18.843 7.7161 18.6554 7.52861ZM16.9484 18.1437H7.39844V5.59561H13.3504V9.19561H16.9504L16.9484 18.1437Z" fill="currentColor"/> <path d="M10.9219 14.5498H8.92188V16.5498H10.9219V14.5498Z" fill="currentColor"/> <path d="M10.9219 7.91992H8.92188V13.4219H10.9219V7.91992Z" fill="currentColor"/>`},undo:{viewBox:`0 0 24 24`,body:`<path d="M12.375 9C10.3875 9 8.5875 9.75 7.2 10.95L4.5 8.25V15H11.25L8.535 12.285C9.5775 11.415 10.905 10.875 12.375 10.875C15.03 10.875 17.2875 12.6075 18.075 15L19.8525 14.415C18.81 11.2725 15.8625 9 12.375 9Z" fill="currentColor"/>`},redo:{viewBox:`0 0 24 24`,body:`<path d="M16.8012 10.95C15.4137 9.75 13.6137 9 11.6262 9C8.13875 9 5.19125 11.2725 4.15625 14.415L5.92625 15C6.71375 12.6075 8.96375 10.875 11.6262 10.875C13.0887 10.875 14.4237 11.415 15.4662 12.285L12.7512 15H19.5012V8.25L16.8012 10.95Z" fill="currentColor"/>`},open:{viewBox:`0 0 48 48`,body:`<path d="M24.0002 33.2303L9.78223 19.0123L12.6102 16.1843L24.0002 27.5743L35.3882 16.1843L38.2182 19.0123L24.0002 33.2303Z" fill="currentColor"/>`},table:{viewBox:`0 0 24 24`,body:`<path d="M7.44474 15.5H4.14632C3.78937 15.5 3.5 15.7686 3.5 16.1V19.162C3.5 19.4934 3.78937 19.762 4.14632 19.762H7.44474C7.80169 19.762 8.09106 19.4934 8.09106 19.162V16.1C8.09106 15.7686 7.80169 15.5 7.44474 15.5Z" fill="currentColor"/> <path d="M13.6518 15.5H10.3534C9.9964 15.5 9.70703 15.7686 9.70703 16.1V19.162C9.70703 19.4934 9.9964 19.762 10.3534 19.762H13.6518C14.0087 19.762 14.2981 19.4934 14.2981 19.162V16.1C14.2981 15.7686 14.0087 15.5 13.6518 15.5Z" fill="currentColor"/> <path d="M19.8549 15.5H16.5565C16.1995 15.5 15.9102 15.7686 15.9102 16.1V19.162C15.9102 19.4934 16.1995 19.762 16.5565 19.762H19.8549C20.2119 19.762 20.5012 19.4934 20.5012 19.162V16.1C20.5012 15.7686 20.2119 15.5 19.8549 15.5Z" fill="currentColor"/> <path d="M5 13.5C4.72386 13.5 4.5 13.2761 4.5 13V5C4.5 4.72386 4.72386 4.5 5 4.5H8.12468C8.78163 4.5 9.34861 4.63184 9.82562 4.89551C10.3055 5.15625 10.6754 5.521 10.9353 5.98975C11.1952 6.45557 11.3252 6.99756 11.3252 7.61572C11.3252 8.23682 11.1924 8.78027 10.9267 9.24609C10.664 9.70898 10.2883 10.0679 9.79991 10.3228C9.31148 10.5776 8.73165 10.7051 8.06041 10.7051H6.32391C6.04777 10.7051 5.82391 10.4812 5.82391 10.2051V9.49121C5.82391 9.21507 6.04777 8.99121 6.32391 8.99121H7.66624C7.98615 8.99121 8.25321 8.93408 8.46744 8.81982C8.68452 8.70557 8.84876 8.5459 8.96015 8.34082C9.07155 8.13281 9.12725 7.89111 9.12725 7.61572C9.12725 7.3374 9.07155 7.09717 8.96015 6.89502C8.84876 6.68994 8.68452 6.53174 8.46744 6.42041C8.25036 6.30908 7.98329 6.25342 7.66624 6.25342H7.12082C6.84468 6.25342 6.62082 6.47728 6.62082 6.75342V13C6.62082 13.2761 6.39697 13.5 6.12082 13.5H5Z" fill="currentColor"/> <path d="M12.5835 6.2666C12.3074 6.2666 12.0835 6.04274 12.0835 5.7666V5C12.0835 4.72386 12.3074 4.5 12.5835 4.5H19C19.2761 4.5 19.5 4.72386 19.5 5V5.7666C19.5 6.04274 19.2761 6.2666 19 6.2666H17.3393C17.0632 6.2666 16.8393 6.49046 16.8393 6.7666V13C16.8393 13.2761 16.6155 13.5 16.3393 13.5H15.2485C14.9724 13.5 14.7485 13.2761 14.7485 13V6.7666C14.7485 6.49046 14.5246 6.2666 14.2485 6.2666H12.5835Z" fill="currentColor"/>`},benzene:{viewBox:`0 0 32 32`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M15.3329 4.53606C15.7455 4.29788 16.2538 4.29788 16.6663 4.53606L25.5943 9.69072C26.0068 9.9289 26.261 10.3691 26.261 10.8454V21.1548C26.261 21.6311 26.0068 22.0713 25.5943 22.3095L16.6663 27.4641C16.2538 27.7023 15.7455 27.7023 15.3329 27.4641L6.40494 22.3095C5.99241 22.0713 5.73828 21.6311 5.73828 21.1548V10.8454C5.73828 10.3691 5.99241 9.9289 6.40494 9.69072L15.3329 4.53606ZM8.40495 11.6152V20.385L15.9996 24.7698L23.5943 20.385V11.6152L15.9996 7.23036L8.40495 11.6152Z" fill="currentColor"/> <path d="M21.0837 13.0666V18.9332L22.4171 19.7039V12.2932L21.0837 13.0666ZM15.9997 8.58789L9.58105 12.2932L10.9144 13.0626L15.9997 10.1332V8.58922V8.58789ZM9.58105 19.7039L15.9997 23.4092V21.8666L10.9144 18.9332L9.58105 19.7039Z" fill="currentColor"/>`},cyclopentadiene:{viewBox:`0 0 32 32`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M15.2157 5.41598C15.683 5.07643 16.3159 5.07643 16.7832 5.41598L25.9365 12.0666C26.4038 12.4062 26.5994 13.008 26.4209 13.5573L22.9249 24.3173C22.7464 24.8667 22.2344 25.2386 21.6568 25.2386H10.3421C9.76448 25.2386 9.25253 24.8667 9.07404 24.3173L5.57804 13.5573C5.39955 13.008 5.59508 12.4062 6.06238 12.0666L15.2157 5.41598ZM8.41353 13.6546L11.3109 22.572H20.6881L23.5854 13.6546L15.9995 8.14277L8.41353 13.6546Z" fill="currentColor"/> <path fill-rule="evenodd" clip-rule="evenodd" d="M22.3492 14.6101L20.1665 21.3274L18.8984 20.9154L21.0811 14.1981L22.3492 14.6101Z" fill="currentColor"/> <path fill-rule="evenodd" clip-rule="evenodd" d="M11.834 21.3274L9.65137 14.6101L10.9194 14.1981L13.1021 20.9154L11.834 21.3274Z" fill="currentColor"/>`},cyclohexane:{viewBox:`0 0 32 32`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M15.3329 4.53606C15.7455 4.29788 16.2538 4.29788 16.6663 4.53606L25.5943 9.69072C26.0068 9.9289 26.261 10.3691 26.261 10.8454V21.1548C26.261 21.6311 26.0068 22.0713 25.5943 22.3095L16.6663 27.4641C16.2538 27.7023 15.7455 27.7023 15.3329 27.4641L6.40494 22.3095C5.99241 22.0713 5.73828 21.6311 5.73828 21.1548V10.8454C5.73828 10.3691 5.99241 9.9289 6.40494 9.69072L15.3329 4.53606ZM8.40495 11.6152V20.385L15.9996 24.7698L23.5943 20.385V11.6152L15.9996 7.23036L8.40495 11.6152Z" fill="currentColor"/>`},cyclopentane:{viewBox:`0 0 32 32`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M15.2157 6.21597C15.683 5.87642 16.3159 5.87642 16.7832 6.21597L25.9365 12.8666C26.4038 13.2062 26.5994 13.808 26.4209 14.3573L22.9249 25.1173C22.7464 25.6667 22.2344 26.0386 21.6568 26.0386H10.3421C9.76448 26.0386 9.25253 25.6667 9.07404 25.1173L5.57804 14.3573C5.39955 13.808 5.59508 13.2062 6.06238 12.8666L15.2157 6.21597ZM8.41353 14.4546L11.3109 23.372H20.6881L23.5854 14.4546L15.9995 8.94276L8.41353 14.4546Z" fill="currentColor"/>`},cyclopropane:{viewBox:`0 0 32 32`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M15.9995 6.71075C16.4759 6.71075 16.9161 6.9649 17.1543 7.37745L25.8156 22.3801C26.0538 22.7927 26.0538 23.3009 25.8156 23.7134C25.5774 24.126 25.1372 24.3801 24.6609 24.3801H7.33822C6.86187 24.3801 6.4217 24.126 6.18352 23.7134C5.94534 23.3009 5.94534 22.7927 6.1835 22.3801L14.8448 7.37745C15.083 6.9649 15.5232 6.71075 15.9995 6.71075ZM15.9995 10.7109L9.64756 21.7134H22.3515L15.9995 10.7109Z" fill="currentColor"/>`},cyclobutane:{viewBox:`0 0 32 32`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M5.32812 6.66128C5.32812 5.9249 5.92508 5.32794 6.66146 5.32794H25.3401C26.0765 5.32794 26.6735 5.9249 26.6735 6.66128V25.3399C26.6735 26.0763 26.0765 26.6733 25.3401 26.6733H6.66146C5.92508 26.6733 5.32812 26.0763 5.32812 25.3399V6.66128ZM7.99479 7.99461V24.0066H24.0068V7.99461H7.99479Z" fill="currentColor"/>`},cycloheptane:{viewBox:`0 0 32 32`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M15.4219 5.62533C15.7874 5.44935 16.2131 5.44935 16.5786 5.62533L24.0146 9.20533C24.3802 9.38135 24.6458 9.71431 24.7361 10.1099L26.5735 18.1579C26.6638 18.5535 26.569 18.9688 26.316 19.286L21.1693 25.7394C20.9163 26.0566 20.5327 26.2413 20.1269 26.2413H11.8736C11.4678 26.2413 11.0841 26.0566 10.8311 25.7394L5.68448 19.286C5.43146 18.9688 5.33669 18.5535 5.42701 18.1579L7.26434 10.1099C7.35466 9.71431 7.62023 9.38135 7.98585 9.20533L15.4219 5.62533ZM9.7212 11.3295L8.16972 18.1254L12.5156 23.5747H19.4848L23.8308 18.1254L22.2793 11.3295L16.0002 8.30649L9.7212 11.3295Z" fill="currentColor"/>`},cyclooctane:{viewBox:`0 0 32 32`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M11.3569 6.12527C11.6069 5.87522 11.9461 5.73474 12.2997 5.73474H19.6997C20.0533 5.73474 20.3925 5.87522 20.6425 6.12527L25.8745 11.3573C26.1246 11.6073 26.265 11.9465 26.265 12.3001V19.7001C26.265 20.0537 26.1246 20.3928 25.8745 20.6429L20.6425 25.8749C20.3925 26.1249 20.0533 26.2654 19.6997 26.2654H12.2997C11.9461 26.2654 11.6069 26.1249 11.3569 25.8749L6.1249 20.6429C5.87485 20.3928 5.73438 20.0537 5.73438 19.7001V12.3001C5.73438 11.9465 5.87485 11.6073 6.1249 11.3573L11.3569 6.12527ZM12.852 8.40141L8.40104 12.8524V19.1478L12.852 23.5987H19.1474L23.5984 19.1478V12.8524L19.1474 8.40141H12.852Z" fill="currentColor"/>`},chain:{viewBox:`0 0 32 32`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M26.3232 10.6201C26.84 11.1446 26.8339 11.9888 26.3094 12.5057L20.8947 17.8417C20.5652 18.1664 20.0908 18.2971 19.6414 18.187L12.6445 16.4722L7.65205 21.3923C7.12757 21.9092 6.28337 21.903 5.76648 21.3786C5.2496 20.8541 5.25576 20.0099 5.78024 19.493L11.3029 14.0503C11.6325 13.7256 12.1068 13.5948 12.5562 13.705L19.5532 15.4198L24.4376 10.6063C24.9621 10.0894 25.8063 10.0956 26.3232 10.6201Z" fill="currentColor"/>`},library:{viewBox:`0 0 24 24`,body:`<path d="M7.44474 15.5H4.14632C3.78937 15.5 3.5 15.7686 3.5 16.1V19.162C3.5 19.4934 3.78937 19.762 4.14632 19.762H7.44474C7.80169 19.762 8.09106 19.4934 8.09106 19.162V16.1C8.09106 15.7686 7.80169 15.5 7.44474 15.5Z" fill="currentColor"/><path d="M13.6518 15.5H10.3534C9.9964 15.5 9.70703 15.7686 9.70703 16.1V19.162C9.70703 19.4934 9.9964 19.762 10.3534 19.762H13.6518C14.0087 19.762 14.2981 19.4934 14.2981 19.162V16.1C14.2981 15.7686 14.0087 15.5 13.6518 15.5Z" fill="currentColor"/><path d="M19.8549 15.5H16.5565C16.1995 15.5 15.9102 15.7686 15.9102 16.1V19.162C15.9102 19.4934 16.1995 19.762 16.5565 19.762H19.8549C20.2119 19.762 20.5012 19.4934 20.5012 19.162V16.1C20.5012 15.7686 20.2119 15.5 19.8549 15.5Z" fill="currentColor"/><path d="M10.8048 7.00003C10.5441 7.00003 10.3455 6.78614 10.2145 6.56077C10.1647 6.47523 10.1011 6.40006 10.0236 6.33526C9.83609 6.17617 9.55058 6.09662 9.16705 6.09662C8.92273 6.09662 8.72245 6.12645 8.5662 6.18611C8.41279 6.24293 8.29916 6.32106 8.22529 6.42049C8.15143 6.51992 8.11308 6.63356 8.11023 6.7614C8.10455 6.86651 8.12302 6.96168 8.16563 7.04691C8.21109 7.12929 8.28211 7.20458 8.3787 7.27276C8.47529 7.3381 8.59887 7.39776 8.74944 7.45174C8.90001 7.50572 9.07899 7.55401 9.28637 7.59662L10.0023 7.75003C10.4852 7.85231 10.8986 7.98725 11.2423 8.15486C11.5861 8.32248 11.8673 8.51992 12.0861 8.74719C12.3048 8.97162 12.4653 9.22447 12.5676 9.50572C12.6727 9.78697 12.7267 10.0938 12.7296 10.4262C12.7267 11 12.5832 11.4858 12.2992 11.8836C12.0151 12.2813 11.6088 12.5838 11.0804 12.7912C10.5548 12.9986 9.92274 13.1023 9.1841 13.1023C8.42558 13.1023 7.76364 12.9901 7.1983 12.7657C6.6358 12.5412 6.1983 12.1961 5.8858 11.7301C5.65013 11.3734 5.50332 10.941 5.44537 10.4328C5.41399 10.1577 5.6421 9.93185 5.919 9.93185H7.19294C7.45876 9.93185 7.66417 10.1468 7.7599 10.3948C7.78906 10.4703 7.82487 10.5404 7.86734 10.6051C7.98666 10.787 8.15427 10.9247 8.37018 11.0185C8.58893 11.1122 8.84887 11.1591 9.15001 11.1591C9.40285 11.1591 9.6145 11.1279 9.78495 11.0654C9.95541 11.0029 10.0847 10.9162 10.1727 10.8054C10.2608 10.6946 10.3063 10.5682 10.3091 10.4262C10.3063 10.2926 10.2622 10.1762 10.177 10.0767C10.0946 9.97447 9.95825 9.88356 9.76791 9.80401C9.57756 9.72162 9.32046 9.64492 8.9966 9.5739L8.12728 9.3864C7.35455 9.21878 6.74518 8.93895 6.29916 8.54691C5.85597 8.15202 5.6358 7.61367 5.63864 6.93185C5.6358 6.37787 5.78353 5.8935 6.08183 5.47873C6.38296 5.06111 6.79916 4.73583 7.33041 4.50287C7.8645 4.26992 8.47671 4.15344 9.16705 4.15344C9.8716 4.15344 10.481 4.27134 10.9952 4.50714C11.5094 4.74293 11.9057 5.07532 12.1841 5.5043C12.3799 5.80099 12.5083 6.1328 12.5691 6.49972C12.6145 6.7731 12.3849 7.00003 12.1077 7.00003H10.8048Z" fill="currentColor"/><path d="M14.1308 13C13.8547 13 13.6308 12.7762 13.6308 12.5V4.77276C13.6308 4.49662 13.8547 4.27276 14.1308 4.27276H15.5001C15.7763 4.27276 16.0001 4.49662 16.0001 4.77276V10.5909C16.0001 10.8671 16.224 11.0909 16.5001 11.0909H19.0286C19.3047 11.0909 19.5286 11.3148 19.5286 11.5909V12.5C19.5286 12.7762 19.3047 13 19.0286 13H14.1308Z" fill="currentColor"/>`},zoomIn:{viewBox:`0 0 16 16`,body:`<path d="M12.6654 8.66683H8.66536V12.6668H7.33203V8.66683H3.33203V7.3335H7.33203V3.3335H8.66536V7.3335H12.6654V8.66683Z" fill="currentColor"/>`},zoomOut:{viewBox:`0 0 16 16`,body:`<path d="M12.6654 8.66683H3.33203V7.3335H12.6654V8.66683Z" fill="currentColor"/>`},zoomReset:{viewBox:`0 0 1024 1024`,body:`<path d="M921.728 654.4a428.352 428.352 0 1 0-647.04 203.264c3.2 2.24 7.488 1.856 10.24-0.896l43.2-43.2a8 8 0 0 0-1.28-12.352 352.256 352.256 0 1 1 462.592-67.968l-70.4-69.76a8 8 0 0 0-13.632 4.608l-27.904 219.776a8 8 0 0 0 9.024 8.96l219.52-29.952a8 8 0 0 0 4.48-13.632l-67.008-66.368a431.36 431.36 0 0 0 78.208-132.48z" fill="currentColor"/>`},hand:{viewBox:`0 0 24 24`,body:`<path d="M20.4338 7.67261C20.3338 7.07261 19.7208 5.94161 18.1138 5.94161C17.9877 5.9414 17.8618 5.94908 17.7368 5.96461V5.17461C17.737 5.07291 17.7218 4.97176 17.6918 4.87461C17.5148 4.31061 16.8158 3.24361 15.4298 3.24361C15.2016 3.24249 14.9743 3.27141 14.7538 3.32961C14.7095 3.25021 14.6608 3.17342 14.6078 3.09961C14.1548 2.45461 13.4218 2.09961 12.5448 2.09961C11.6678 2.09961 10.9608 2.46061 10.5198 3.11661C10.3468 3.37494 10.2243 3.66369 10.1588 3.96761C10.0237 3.94716 9.88733 3.93713 9.75076 3.93761C8.27676 3.93761 7.62676 5.08761 7.48876 5.69361C7.47277 5.76653 7.46472 5.84096 7.46476 5.91561V11.0996L7.33576 10.9996C6.74676 10.4996 5.30076 9.87861 3.94376 11.0286C2.58676 12.1786 2.92876 13.6926 3.51276 14.4806C3.96576 15.2536 8.02476 21.9056 13.7358 21.9056C16.9868 21.9056 18.6458 20.3936 19.4648 19.1246C20.0103 18.2829 20.3467 17.3228 20.4458 16.3246C20.4458 16.3036 20.4458 16.2826 20.4458 16.2616V7.83761C20.4462 7.78238 20.4422 7.72719 20.4338 7.67261ZM18.4478 16.2146C18.3608 16.8886 18.1214 17.534 17.7478 18.1016C16.9428 19.3016 15.5948 19.9016 13.7398 19.9016C8.97176 19.9016 5.25576 13.4896 5.21976 13.4256C5.19504 13.3803 5.16518 13.338 5.13076 13.2996C4.83576 12.8916 5.13076 12.6446 5.23676 12.5516C5.50376 12.3246 5.77676 12.3146 6.04776 12.5206L7.81676 13.9726C7.96236 14.0952 8.13991 14.1736 8.32855 14.1988C8.5172 14.224 8.70911 14.1948 8.88174 14.1147C9.05437 14.0345 9.20056 13.9068 9.30312 13.7465C9.40569 13.5862 9.46037 13.3999 9.46076 13.2096V6.08061C9.48827 6.02928 9.53075 5.98754 9.58255 5.96093C9.63435 5.93432 9.69301 5.92409 9.75076 5.93161C9.82113 5.92666 9.89166 5.9387 9.9564 5.96672C10.0211 5.99474 10.0782 6.03792 10.1228 6.09261V10.7926C10.1228 11.0578 10.2281 11.3122 10.4157 11.4997C10.6032 11.6873 10.8575 11.7926 11.1228 11.7926C11.388 11.7926 11.6423 11.6873 11.8299 11.4997C12.0174 11.3122 12.1228 11.0578 12.1228 10.7926V4.34161C12.1908 4.09361 12.4228 4.09361 12.5448 4.09361C12.7448 4.09361 12.9588 4.12561 13.0268 4.34861V10.7896C13.0268 11.0548 13.1321 11.3092 13.3197 11.4967C13.5072 11.6843 13.7615 11.7896 14.0268 11.7896C14.292 11.7896 14.5463 11.6843 14.7339 11.4967C14.9214 11.3092 15.0268 11.0548 15.0268 10.7896V5.58961C15.0427 5.52323 15.0676 5.45932 15.1008 5.39961C15.1498 5.31761 15.2088 5.24961 15.4258 5.24961C15.4855 5.24639 15.545 5.25848 15.5987 5.28473C15.6525 5.31098 15.6986 5.35053 15.7328 5.39961V10.7896C15.7328 11.0548 15.8381 11.3092 16.0257 11.4967C16.2132 11.6843 16.4675 11.7896 16.7328 11.7896C16.998 11.7896 17.2523 11.6843 17.4399 11.4967C17.6274 11.3092 17.7328 11.0548 17.7328 10.7896V8.04561C17.8432 7.96998 17.9763 7.93467 18.1098 7.94561C18.1748 7.94561 18.3678 7.95161 18.4438 8.02361L18.4478 16.2146Z" fill="currentColor"/>`},fit:{viewBox:`0 0 24 24`,body:`<path d="M4 9V4h5M15 4h5v5M20 15v5h-5M9 20H4v-5" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>`},layout:{viewBox:`0 0 24 24`,body:`<path d="M20.0546 15.0566C20.0608 15.3433 20.0097 15.6284 19.9042 15.8951C19.7988 16.1617 19.6411 16.4047 19.4406 16.6096C19.24 16.8146 19.0005 16.9774 18.7361 17.0886C18.4718 17.1998 18.1879 17.257 17.9011 17.257C17.6144 17.257 17.3305 17.1998 17.0661 17.0886C16.8018 16.9774 16.5623 16.8146 16.3617 16.6096C16.1611 16.4047 16.0035 16.1617 15.898 15.8951C15.7926 15.6284 15.7414 15.3433 15.7476 15.0566C15.7474 14.9424 15.7567 14.8283 15.7756 14.7156L14.0196 13.6996L13.1356 13.1886L12.2526 13.6986L10.7616 14.5596C10.8008 14.7213 10.8202 14.8872 10.8196 15.0536C10.8205 15.5239 10.6673 15.9816 10.3833 16.3566C10.0994 16.7316 9.70051 17.0032 9.24758 17.13C8.79466 17.2567 8.31267 17.2316 7.87537 17.0584C7.43808 16.8853 7.06955 16.5736 6.82616 16.1712C6.58278 15.7687 6.47795 15.2976 6.5277 14.8299C6.57746 14.3622 6.77907 13.9237 7.10168 13.5815C7.42429 13.2392 7.85014 13.0121 8.31408 12.9348C8.77802 12.8576 9.2545 12.9344 9.67062 13.1536L11.3706 12.1736L12.2546 11.6636V8.73761C11.8052 8.53548 11.4386 8.18534 11.2161 7.74562C10.9937 7.3059 10.9287 6.80318 11.032 6.32134C11.1354 5.8395 11.4009 5.40767 11.7841 5.0979C12.1674 4.78812 12.6453 4.61914 13.1381 4.61914C13.6309 4.61914 14.1088 4.78812 14.4921 5.0979C14.8754 5.40767 15.1408 5.8395 15.2442 6.32134C15.3476 6.80318 15.2826 7.3059 15.0601 7.74562C14.8376 8.18534 14.4711 8.53548 14.0216 8.73761V11.6626L14.9046 12.1716L16.7516 13.2386C17.077 13.0334 17.4515 12.9191 17.836 12.9078C18.2205 12.8964 18.6011 12.9883 18.938 13.1739C19.275 13.3596 19.556 13.6321 19.7519 13.9632C19.9477 14.2944 20.0512 14.6719 20.0516 15.0566H20.0546Z" fill="currentColor"/> <path fill-rule="evenodd" clip-rule="evenodd" d="M4 3C4.55228 3 5 3.44772 5 4V19H20C20.5523 19 21 19.4477 21 20C21 20.5523 20.5523 21 20 21H4C3.44772 21 3 20.5523 3 20V4C3 3.44772 3.44772 3 4 3Z" fill="currentColor"/>`},clean:{viewBox:`0 0 24 24`,body:`<path d="M19.8484 16.3124C19.8484 16.9587 19.5917 17.5785 19.1347 18.0356C18.6776 18.4926 18.0578 18.7494 17.4114 18.7494C16.7651 18.7494 16.1452 18.4926 15.6882 18.0356C15.2312 17.5785 14.9744 16.9587 14.9744 16.3124C14.9741 16.1831 14.9845 16.0539 15.0054 15.9264L13.0184 14.7794L12.0184 14.2014L11.0184 14.7784L9.33144 15.7524C9.37582 15.9357 9.39798 16.1237 9.39744 16.3124C9.39733 16.8448 9.22284 17.3626 8.90064 17.7865C8.57845 18.2104 8.12628 18.5171 7.61328 18.6598C7.10028 18.8024 6.55468 18.7731 6.0599 18.5764C5.56513 18.3796 5.14841 18.0262 4.87347 17.5702C4.59852 17.1143 4.48049 16.5808 4.5374 16.0514C4.59432 15.5219 4.82306 15.0257 5.18864 14.6386C5.55423 14.2515 6.03654 13.9948 6.56183 13.9077C7.08712 13.8206 7.62648 13.9079 8.09744 14.1564L10.0194 13.0474L11.0194 12.4694V9.16036C10.511 8.93159 10.0964 8.53541 9.84479 8.03792C9.59316 7.54042 9.51971 6.97169 9.6367 6.42659C9.75369 5.88149 10.054 5.39298 10.4877 5.04256C10.9213 4.69213 11.4619 4.50098 12.0194 4.50098C12.577 4.50098 13.1176 4.69213 13.5512 5.04256C13.9848 5.39298 14.2852 5.88149 14.4022 6.42659C14.5192 6.97169 14.4457 7.54042 14.1941 8.03792C13.9425 8.53541 13.5279 8.93159 13.0194 9.16036V12.4704L14.0194 13.0474L16.1094 14.2544C16.4779 14.0206 16.9023 13.8901 17.3384 13.8765C17.7745 13.8629 18.2063 13.9667 18.5886 14.177C18.9709 14.3873 19.2897 14.6964 19.5116 15.072C19.7336 15.4477 19.8506 15.876 19.8504 16.3124H19.8484Z" fill="currentColor"/>`},rgroupLabel:{viewBox:`0 0 32 32`,body:`<path d="M27.9272 22.052H25.2605V11.2813H22.8525V8.61465H27.9272V22.052Z" fill="currentColor"/> <path d="M19.3963 16.176C20.8536 15.5093 21.7576 14.2533 21.7576 12.6947C21.7576 11.228 21.0616 10.0173 19.8003 9.30401C18.7003 8.68268 17.4989 8.61601 17.0349 8.61601H11.7803V22.052H14.4469V16.6933C15.1896 16.6933 16.0816 16.6933 16.6909 16.6933L19.2243 22.0547H22.1736L19.3963 16.176ZM14.4469 14.0267V11.2827H17.0349C17.5393 11.2795 18.0372 11.397 18.4869 11.6253C18.9163 11.868 19.0909 12.1773 19.0909 12.692C19.0909 13.6947 17.7229 13.996 16.9989 14.0147C16.4749 14.032 15.3496 14.0307 14.4469 14.0267Z" fill="currentColor"/> <path d="M6.3492 19.1373L4.46387 17.2507L6.41187 15.304L4.46387 13.356L6.3492 11.4693L10.1825 15.304L6.3492 19.1373Z" fill="currentColor"/>`},attachmentPoint:{viewBox:`0 0 32 32`,body:`<path d="M24.9626 13.3347H23.3626V9.07468H22.5986V7.47334H24.9626V13.3347Z" fill="currentColor"/> <path d="M21.4667 10.9787C21.7432 10.8234 21.9736 10.5976 22.1343 10.3242C22.295 10.0509 22.3802 9.73977 22.3814 9.42268C22.3814 8.48934 21.6867 7.47334 20.1641 7.47334H17.7334V13.34H19.3334V11.3333C19.5187 11.3333 19.7107 11.3333 19.8761 11.3333L20.8201 13.3333H22.5894L21.4667 10.9787ZM19.3334 9.73334V9.06668H20.1734C20.3174 9.06668 20.7907 9.09201 20.7907 9.41601C20.7907 9.61868 20.4174 9.71868 20.1521 9.72668C19.9734 9.73334 19.6427 9.73334 19.3334 9.73334Z" fill="currentColor"/> <path d="M18.2586 17.1187H16.5652C16.1368 17.118 15.7262 16.9475 15.4233 16.6446C15.1204 16.3417 14.9499 15.9311 14.9492 15.5027V5.95871C14.9499 5.53046 15.1204 5.11997 15.4234 4.81727C15.7263 4.51458 16.137 4.34439 16.5652 4.34404H18.2586V6.34404H16.9492V15.1187H18.2586V17.1187Z" fill="currentColor"/> <path d="M26.1963 17.1187H24.5029V15.1187H25.8123V6.34404H24.5029V4.34404H26.1963C26.6245 4.34439 27.0351 4.51458 27.3381 4.81727C27.641 5.11997 27.8116 5.53046 27.8123 5.95871V15.5027C27.8116 15.9311 27.6411 16.3417 27.3382 16.6446C27.0353 16.9475 26.6246 17.118 26.1963 17.1187Z" fill="currentColor"/> <path d="M13.88 16.7557L7.85547 22.7803L9.26968 24.1945L15.2942 18.1699L13.88 16.7557Z" fill="currentColor"/> <path d="M11.4794 25.0853L5.31543 26.7373L6.9661 20.572L11.4794 25.0853Z" fill="currentColor"/>`},selectRect:{viewBox:`0 0 48 48`,body:`<path d="M36.2109 42.2441C36.1045 42.2985 35.9877 42.3299 35.8683 42.3361C35.7489 42.3424 35.6295 42.3234 35.518 42.2803C35.4064 42.2373 35.3052 42.1712 35.2209 42.0864C35.1367 42.0016 35.0712 41.8999 35.0289 41.7881L33.0909 37.1881L30.8609 39.1501C30.7105 39.2885 30.5133 39.365 30.3089 39.3641C30.1902 39.3632 30.073 39.3375 29.9648 39.2886C29.8567 39.2397 29.76 39.1687 29.6809 39.0801C29.512 38.8909 29.4194 38.6457 29.4209 38.3921V23.8121C29.4198 23.5585 29.5124 23.3135 29.6809 23.1241C29.76 23.0355 29.8567 22.9645 29.9648 22.9156C30.073 22.8667 30.1902 22.8409 30.3089 22.8401C30.5199 22.8428 30.7225 22.9225 30.8789 23.0641L41.0909 32.4301C41.2744 32.6048 41.3876 32.8406 41.4089 33.0931C41.4303 33.3455 41.3584 33.597 41.2069 33.8001C41.0831 33.967 40.9036 34.0841 40.7009 34.1301L37.9009 34.7301L39.8569 39.3301C39.9649 39.5575 39.9809 39.8178 39.9016 40.0567C39.8222 40.2956 39.6536 40.4946 39.4309 40.6121L36.2109 42.2441Z" fill="currentColor"/> <path d="M8.56055 20.1022V25.5802H12.5605V20.1022H8.56055ZM19.9405 32.9622V36.9622H25.4205V32.9622H19.9405ZM31.5265 8.72021H29.4205V12.7202H31.5265C31.6943 12.7205 31.8603 12.7538 32.0152 12.8183C32.17 12.8828 32.3106 12.9772 32.4289 13.0961C32.5473 13.215 32.641 13.356 32.7048 13.5111C32.7685 13.6663 32.8011 13.8325 32.8005 14.0002V16.1022H36.8005V14.0002C36.8 12.6011 36.2444 11.2593 35.2556 10.2694C34.2668 9.27949 32.9257 8.72233 31.5265 8.72021ZM13.8365 8.72021C12.4371 8.7218 11.0954 9.27873 10.1062 10.2687C9.11698 11.2587 8.56108 12.6007 8.56055 14.0002V16.1022H12.5605V14.0002C12.5605 13.6618 12.695 13.3372 12.9343 13.0979C13.1736 12.8586 13.4981 12.7242 13.8365 12.7242H15.9405V8.72421L13.8365 8.72021ZM19.9405 8.72021V12.7202H25.4205V8.72021H19.9405ZM32.8005 20.1022V21.8002L36.8005 25.4682V20.1022H32.8005ZM13.8365 32.9622C13.4981 32.9622 13.1736 32.8278 12.9343 32.5885C12.695 32.3492 12.5605 32.0246 12.5605 31.6862V29.5802H8.56055V31.6862C8.56213 33.085 9.11851 34.4261 10.1076 35.4152C11.0967 36.4043 12.4378 36.9606 13.8365 36.9622H15.9405V32.9622H13.8365Z" fill="currentColor"/>`},selectLasso:{viewBox:`0 0 48 48`,body:`<path d="M36.209 42.2438C36.1025 42.2983 35.9858 42.3296 35.8664 42.3359C35.747 42.3421 35.6276 42.3231 35.516 42.2801C35.4045 42.2371 35.3032 42.171 35.219 42.0861C35.1347 42.0013 35.0693 41.8997 35.027 41.7878L33.089 37.1878L30.859 39.1498C30.7086 39.2883 30.5114 39.3647 30.307 39.3638C30.1883 39.363 30.071 39.3373 29.9629 39.2884C29.8547 39.2394 29.758 39.1684 29.679 39.0798C29.5101 38.8906 29.4174 38.6455 29.419 38.3918V23.8118C29.4178 23.5583 29.5104 23.3133 29.679 23.1238C29.758 23.0353 29.8547 22.9643 29.9629 22.9153C30.071 22.8664 30.1883 22.8407 30.307 22.8398C30.5179 22.8426 30.7206 22.9222 30.877 23.0638V23.0638L41.089 32.4298C41.2725 32.6045 41.3856 32.8404 41.407 33.0928C41.4284 33.3453 41.3565 33.5968 41.205 33.7998C41.0811 33.9668 40.9016 34.0838 40.699 34.1298L37.899 34.7298L39.855 39.3298C39.963 39.5572 39.979 39.8176 39.8996 40.0565C39.8202 40.2954 39.6516 40.4943 39.429 40.6118L36.209 42.2438Z" fill="currentColor"/> <path d="M30.6734 8.50393C30.4277 8.03388 30.0054 7.68066 29.4994 7.52193C28.264 7.13666 26.9971 6.86093 25.7134 6.69793C25.3334 6.64993 24.9374 6.60993 24.5394 6.57993C24.0113 6.5446 23.4906 6.71875 23.09 7.06465C22.6895 7.41056 22.4413 7.90036 22.3994 8.42793C22.3994 8.47793 22.3994 8.52993 22.3994 8.58193C22.4011 9.08594 22.593 9.57073 22.9367 9.93933C23.2805 10.3079 23.7507 10.5331 24.2534 10.5699C24.5694 10.5919 24.8934 10.6259 25.2154 10.6659C26.2653 10.7987 27.3013 11.0243 28.3114 11.3399C28.5056 11.4005 28.7079 11.4315 28.9114 11.4319C29.2555 11.432 29.5938 11.3432 29.8937 11.1743C30.1935 11.0053 30.4446 10.7619 30.6228 10.4674C30.8009 10.173 30.9001 9.83762 30.9108 9.49364C30.9214 9.14967 30.8431 8.80878 30.6834 8.50393H30.6734ZM18.5934 8.90393C18.4528 8.51883 18.1972 8.18625 17.8613 7.95122C17.5254 7.71619 17.1254 7.59007 16.7154 7.58993C16.4812 7.5902 16.2489 7.63219 16.0294 7.71393C14.3808 8.30305 12.8383 9.15469 11.4614 10.2359C11.137 10.4945 10.901 10.8476 10.7861 11.2462C10.6713 11.6449 10.6832 12.0694 10.8203 12.4609C10.9574 12.8525 11.2129 13.1917 11.5514 13.4316C11.8898 13.6715 12.2945 13.8002 12.7094 13.7999C13.1618 13.8016 13.601 13.6477 13.9534 13.3639C14.9922 12.5492 16.1559 11.9076 17.3994 11.4639C17.6469 11.3746 17.8744 11.2373 18.0687 11.0598C18.2631 10.8824 18.4205 10.6684 18.532 10.43C18.6434 10.1916 18.7068 9.93358 18.7184 9.67068C18.7299 9.40778 18.6895 9.14518 18.5994 8.89793L18.5934 8.90393ZM40.1034 16.1559C39.3434 14.5377 38.2857 13.0769 36.9854 11.8499C36.7976 11.6605 36.5738 11.5105 36.3271 11.4088C36.0805 11.3072 35.816 11.2559 35.5492 11.258C35.2824 11.2601 35.0188 11.3155 34.7738 11.4211C34.5288 11.5266 34.3073 11.6801 34.1225 11.8725C33.9377 12.0649 33.7933 12.2923 33.6977 12.5414C33.6021 12.7905 33.5574 13.0561 33.566 13.3228C33.5747 13.5894 33.6366 13.8516 33.7481 14.094C33.8596 14.3363 34.0185 14.5539 34.2154 14.7339C35.1537 15.6151 35.9183 16.6646 36.4694 17.8279C36.5792 18.0665 36.7349 18.2812 36.9276 18.4596C37.1204 18.6381 37.3464 18.7768 37.5927 18.8679C37.8391 18.959 38.101 19.0007 38.3634 18.9905C38.6259 18.9804 38.8838 18.9187 39.1224 18.8089C39.361 18.6991 39.5756 18.5434 39.7541 18.3507C39.9325 18.158 40.0712 17.932 40.1623 17.6856C40.2534 17.4393 40.2951 17.1774 40.285 16.9149C40.2749 16.6524 40.2132 16.3945 40.1034 16.1559V16.1559ZM11.3634 22.1199C11.0077 21.2472 10.8224 20.3144 10.8174 19.3719C10.817 19.1086 10.8317 18.8455 10.8614 18.5839C10.8915 18.323 10.8699 18.0587 10.7978 17.8061C10.7257 17.5535 10.6046 17.3176 10.4414 17.1118C10.2781 16.906 10.0759 16.7344 9.84637 16.6068C9.6168 16.4791 9.36433 16.398 9.10338 16.3679C9.02631 16.3599 8.94887 16.3559 8.87138 16.3559C8.38129 16.3571 7.9087 16.5383 7.54333 16.8649C7.17796 17.1916 6.94525 17.641 6.88938 18.1279C6.84125 18.5408 6.81721 18.9562 6.81738 19.3719C6.82223 20.8217 7.10535 22.257 7.65138 23.5999C7.85612 24.0814 8.2414 24.4636 8.7246 24.6643C9.20779 24.865 9.75041 24.8684 10.2361 24.6737C10.7217 24.4791 11.1118 24.1018 11.3225 23.6229C11.5333 23.144 11.5479 22.6015 11.3634 22.1119V22.1199ZM39.2294 22.7199C38.9774 22.6087 38.7049 22.5515 38.4294 22.5519C38.0401 22.5517 37.6592 22.6652 37.3334 22.8783C37.0076 23.0915 36.7511 23.3951 36.5954 23.7519C36.3513 24.3081 36.0418 24.8331 35.6734 25.3159L38.6294 28.0279C39.2913 27.2134 39.8403 26.3133 40.2614 25.3519C40.3726 25.0999 40.4298 24.8274 40.4294 24.5519C40.4292 24.163 40.3156 23.7825 40.1024 23.4571C39.8893 23.1317 39.5859 22.8755 39.2294 22.7199ZM17.4734 28.1199C16.2551 27.5431 15.124 26.7978 14.1134 25.9059C13.9189 25.7229 13.6898 25.5807 13.4395 25.4877C13.1893 25.3946 12.9229 25.3527 12.6561 25.3642C12.3894 25.3758 12.1276 25.4406 11.8863 25.555C11.645 25.6693 11.4291 25.8308 11.2512 26.0299C11.0733 26.229 10.9371 26.4618 10.8505 26.7144C10.764 26.967 10.729 27.2344 10.7475 27.5007C10.7659 27.7671 10.8376 28.0271 10.9581 28.2653C11.0787 28.5036 11.2457 28.7153 11.4494 28.8879C12.7398 30.0306 14.1846 30.986 15.7414 31.7259C16.0104 31.8575 16.3059 31.9259 16.6054 31.9259C16.9815 31.9242 17.3495 31.8171 17.6677 31.6166C17.9859 31.4161 18.2415 31.1304 18.4054 30.7919C18.537 30.5215 18.6054 30.2247 18.6054 29.9239C18.6045 29.5471 18.498 29.1782 18.2979 28.8589C18.0978 28.5397 17.8121 28.2829 17.4734 28.1179V28.1199Z" fill="currentColor"/> <path d="M27.6912 29.6101V36.8541C26.3612 38.0881 24.1972 38.8801 22.0012 38.7421C20.9132 38.6741 19.7092 38.2841 18.4372 37.8721C17.1652 37.4601 15.7112 36.9861 14.8652 37.0541C13.3652 37.1701 11.8112 38.4201 11.3412 38.8901C10.9662 39.2654 10.4574 39.4763 9.92689 39.4765C9.39635 39.4767 8.88747 39.2661 8.51219 38.8911C8.1369 38.5161 7.92597 38.0074 7.92578 37.4768C7.92559 36.9463 8.13617 36.4374 8.51119 36.0621C8.79319 35.7821 11.3452 33.3141 14.5572 33.0621C16.1912 32.9361 17.9572 33.5081 19.6692 34.0621C20.6512 34.3801 21.6692 34.7101 22.2512 34.7461C23.6232 34.8321 24.6772 34.2361 25.0132 33.8921C25.0132 33.8661 25.0132 33.8381 24.9992 33.8081C24.6812 33.8081 24.3552 33.7981 24.0312 33.7821C23.5016 33.7567 23.0038 33.5221 22.6471 33.1298C22.2905 32.7375 22.1042 32.2197 22.1292 31.6901C22.1534 31.1761 22.3747 30.6913 22.7472 30.3363C23.1197 29.9813 23.6146 29.7835 24.1292 29.7841H24.2252C24.5132 29.7961 24.7992 29.8041 25.0752 29.8041C25.9509 29.8072 26.8255 29.7424 27.6912 29.6101V29.6101Z" fill="currentColor"/>`},selectStructure:{viewBox:`0 0 48 48`,body:`<path d="M36.085 41.8959C35.9784 41.9503 35.8615 41.9816 35.742 41.9877C35.6225 41.9938 35.5031 41.9746 35.3915 41.9314C35.2799 41.8882 35.1787 41.8219 35.0946 41.7368C35.0104 41.6518 34.9451 41.5499 34.903 41.4379L32.965 36.8379L30.733 38.7999C30.5836 38.9388 30.3869 39.0153 30.183 39.0139C30.0643 39.0131 29.9471 38.9874 29.8389 38.9385C29.7308 38.8895 29.634 38.8185 29.555 38.7299C29.3851 38.5412 29.2916 38.2959 29.293 38.0419V23.4639C29.2912 23.2099 29.3847 22.9645 29.555 22.7759C29.634 22.6874 29.7308 22.6163 29.8389 22.5674C29.9471 22.5185 30.0643 22.4928 30.183 22.4919C30.3934 22.4946 30.5954 22.5743 30.751 22.7159H30.761L40.973 32.0819C41.1561 32.2568 41.269 32.4927 41.2903 32.745C41.3117 32.9973 41.2401 33.2487 41.089 33.4519C40.9648 33.6194 40.7845 33.7365 40.581 33.7819L37.781 34.3819L39.737 38.9819C39.845 39.2099 39.8607 39.4707 39.781 39.7099C39.7012 39.9492 39.5321 40.1484 39.309 40.2659L36.085 41.8959Z" fill="currentColor"/> <path d="M34.3705 14.8379L21.0665 8.23792C20.6874 8.05115 20.2589 7.98925 19.8424 8.06111C19.426 8.13298 19.043 8.3349 18.7485 8.63792L8.34845 19.2519C8.05316 19.5539 7.86103 19.9415 7.7996 20.3594C7.73816 20.7772 7.81057 21.2038 8.00645 21.5779L14.8885 34.7419C15.0602 35.0636 15.3166 35.3322 15.63 35.5187C15.9434 35.7052 16.3018 35.8024 16.6665 35.7999C16.7784 35.8005 16.8902 35.7911 17.0005 35.7719L27.2905 34.0299V29.9759L17.7665 31.5879L12.2305 20.9999L20.5885 12.4639L31.2905 17.7759L30.9485 20.0879L34.5105 23.3539L35.4605 16.9219C35.5222 16.5039 35.45 16.0772 35.2541 15.7028C35.0583 15.3284 34.749 15.0256 34.3705 14.8379Z" fill="currentColor"/>`},rotate:{viewBox:`0 0 48 48`,body:`<path d="M30.4892 35.5656L22.5518 42.2926L21.9292 38.8543C19.5112 38.5171 17.2121 37.5949 15.2315 36.1676C13.2509 34.7403 11.6485 32.8512 10.5636 30.6642C9.47864 28.4772 8.94391 26.0584 9.0058 23.6179C9.06769 21.1773 9.72433 18.7888 10.9187 16.6595L13.95 18.3639C13.0814 19.9136 12.5829 21.6432 12.4934 23.4174C12.4038 25.1917 12.7258 26.9627 13.434 28.5919C14.1422 30.2212 15.2175 31.6647 16.576 32.8095C17.9344 33.9544 19.5391 34.7697 21.2648 35.1917L20.6961 32.0613L30.4892 35.5656Z" fill="currentColor"/> <path d="M37.0822 31.3413L34.0509 29.637C34.9195 28.0873 35.418 26.3577 35.5075 24.5834C35.597 22.8091 35.2751 21.0382 34.5669 19.4089C33.8587 17.7797 32.7834 16.3362 31.4249 15.1913C30.0665 14.0464 28.4618 13.2311 26.7361 12.8091L27.3048 15.9483L17.5117 12.4352L25.4491 5.70825L26.0717 9.14652C28.4897 9.48372 30.7887 10.406 32.7694 11.8333C34.75 13.2605 36.3524 15.1497 37.4373 17.3367C38.5222 19.5237 39.057 21.9425 38.9951 24.383C38.9332 26.8236 38.2766 29.2121 37.0822 31.3413Z" fill="currentColor"/>`},flipH:{viewBox:`0 0 48 48`,body:`<path d="M22.001 6L22.001 42H26.001V6H22.001Z" fill="currentColor"/> <path d="M18.271 40H11.601C10.6462 40 9.73052 39.6207 9.05539 38.9456C8.38026 38.2705 8.00098 37.3548 8.00098 36.4V11.6C8.00098 10.6452 8.38026 9.72955 9.05539 9.05442C9.73052 8.37928 10.6462 8 11.601 8H18.271V12H12.001V36H18.271V40Z" fill="currentColor"/> <path d="M36.001 8H36.401C37.3558 8 38.2714 8.37928 38.9466 9.05442C39.6217 9.72955 40.001 10.6452 40.001 11.6V12H36.001V8Z" fill="currentColor"/> <path d="M40.001 36V36.4C40.001 37.3548 39.6217 38.2705 38.9466 38.9456C38.2714 39.6207 37.3558 40 36.401 40H36.001V36H40.001Z" fill="currentColor"/> <path d="M34.535 36H30.001V40H34.535V36Z" fill="currentColor"/> <path d="M34.535 8H30.001V12H34.535V8Z" fill="currentColor"/> <path d="M36.001 14.5981V19.1321H40.001V14.5981H36.001Z" fill="currentColor"/> <path d="M36.001 21.732V26.266H40.001V21.732H36.001Z" fill="currentColor"/> <path d="M36.001 28.866V33.4H40.001V28.866H36.001Z" fill="currentColor"/>`},flipV:{viewBox:`0 0 48 48`,body:`<path d="M42 21.9996H6V25.9996H42V21.9996Z" fill="currentColor"/> <path d="M40 18.2696H36V11.9996H12V18.2696H8V11.5996C8 10.6449 8.37928 9.72918 9.05442 9.05405C9.72955 8.37892 10.6452 7.99963 11.6 7.99963H36.4C37.3548 7.99963 38.2705 8.37892 38.9456 9.05405C39.6207 9.72918 40 10.6449 40 11.5996V18.2696Z" fill="currentColor"/> <path d="M40 36.0156V36.4156C40 37.3704 39.6207 38.2861 38.9456 38.9612C38.2705 39.6363 37.3548 40.0156 36.4 40.0156H36V36.0156H40Z" fill="currentColor"/> <path d="M12 40.0157H11.6C10.6452 40.0157 9.72955 39.6365 9.05442 38.9613C8.37928 38.2862 8 37.3705 8 36.4157L8 36.0157H12V40.0157Z" fill="currentColor"/> <path d="M12 34.5496V30.0156H8V34.5496H12Z" fill="currentColor"/> <path d="M40 34.5496V30.0156H36V34.5496H40Z" fill="currentColor"/> <path d="M33.3992 36.0156H28.8652V40.0156H33.3992V36.0156Z" fill="currentColor"/> <path d="M26.2674 36.0156H21.7334V40.0156H26.2674V36.0156Z" fill="currentColor"/> <path d="M19.1346 36.0156H14.6006V40.0156H19.1346V36.0156Z" fill="currentColor"/>`},delete:{viewBox:`0 0 24 24`,body:`<path d="M15.8128 20.4834H8.18281C7.70542 20.4834 7.24759 20.2938 6.91002 19.9562C6.57245 19.6186 6.38281 19.1608 6.38281 18.6834V7.35742H17.6128V18.6834C17.6128 19.1608 17.4232 19.6186 17.0856 19.9562C16.748 20.2938 16.2902 20.4834 15.8128 20.4834ZM8.38281 18.4834H15.6128V9.35742H8.38281V18.4834Z" fill="currentColor"/> <path d="M19.3053 4.61291V6.61291H4.69531V4.61291H8.12831L9.16231 3.48291H14.8383L15.8723 4.61291H19.3053Z" fill="currentColor"/>`},cut:{viewBox:`0 0 24 24`,body:`<path d="M17.25 5.25L12.75 9.75L14.25 11.25L19.5 6V5.25H17.25ZM12 12.375C11.9005 12.375 11.8052 12.3355 11.7348 12.2652C11.6645 12.1948 11.625 12.0995 11.625 12C11.625 11.9005 11.6645 11.8052 11.7348 11.7348C11.8052 11.6645 11.9005 11.625 12 11.625C12.0995 11.625 12.1948 11.6645 12.2652 11.7348C12.3355 11.8052 12.375 11.9005 12.375 12C12.375 12.0995 12.3355 12.1948 12.2652 12.2652C12.1948 12.3355 12.0995 12.375 12 12.375ZM7.5 18C7.10217 18 6.72064 17.842 6.43934 17.5607C6.15803 17.2794 6 16.8978 6 16.5C6 15.6675 6.675 15 7.5 15C7.89782 15 8.27935 15.158 8.56066 15.4393C8.84196 15.7206 9 16.1022 9 16.5C9 17.3325 8.325 18 7.5 18ZM7.5 9C7.10217 9 6.72064 8.84196 6.43934 8.56066C6.15803 8.27935 6 7.89782 6 7.5C6 6.6675 6.675 6 7.5 6C7.89782 6 8.27935 6.15803 8.56066 6.43934C8.84196 6.72064 9 7.10217 9 7.5C9 8.3325 8.325 9 7.5 9ZM10.23 8.73C10.4025 8.355 10.5 7.9425 10.5 7.5C10.5 6.70435 10.1839 5.94129 9.62132 5.37868C9.05871 4.81607 8.29565 4.5 7.5 4.5C6.70435 4.5 5.94129 4.81607 5.37868 5.37868C4.81607 5.94129 4.5 6.70435 4.5 7.5C4.5 8.29565 4.81607 9.05871 5.37868 9.62132C5.94129 10.1839 6.70435 10.5 7.5 10.5C7.9425 10.5 8.355 10.4025 8.73 10.23L10.5 12L8.73 13.77C8.355 13.5975 7.9425 13.5 7.5 13.5C6.70435 13.5 5.94129 13.8161 5.37868 14.3787C4.81607 14.9413 4.5 15.7043 4.5 16.5C4.5 17.2956 4.81607 18.0587 5.37868 18.6213C5.94129 19.1839 6.70435 19.5 7.5 19.5C8.29565 19.5 9.05871 19.1839 9.62132 18.6213C10.1839 18.0587 10.5 17.2956 10.5 16.5C10.5 16.0575 10.4025 15.645 10.23 15.27L12 13.5L17.25 18.75H19.5V18L10.23 8.73Z" fill="currentColor"/>`},copy:{viewBox:`0 0 24 24`,body:`<path d="M17.25 18.75H9V8.25H17.25V18.75ZM17.25 6.75H9C8.60217 6.75 8.22064 6.90803 7.93934 7.18934C7.65803 7.47064 7.5 7.85217 7.5 8.25V18.75C7.5 19.1478 7.65803 19.5294 7.93934 19.8107C8.22064 20.092 8.60217 20.25 9 20.25H17.25C17.6478 20.25 18.0294 20.092 18.3107 19.8107C18.592 19.5294 18.75 19.1478 18.75 18.75V8.25C18.75 7.85217 18.592 7.47064 18.3107 7.18934C18.0294 6.90803 17.6478 6.75 17.25 6.75ZM15 3.75H6C5.60217 3.75 5.22064 3.90804 4.93934 4.18934C4.65804 4.47064 4.5 4.85217 4.5 5.25V15.75H6V5.25H15V3.75Z" fill="currentColor"/>`},paste:{viewBox:`0 0 24 24`,body:`<path d="M17.25 18H6.75V6H8.25V8.25H15.75V6H17.25V18ZM12 4.5C12.1989 4.5 12.3897 4.57902 12.5303 4.71967C12.671 4.86032 12.75 5.05109 12.75 5.25C12.75 5.44891 12.671 5.63968 12.5303 5.78033C12.3897 5.92098 12.1989 6 12 6C11.8011 6 11.6103 5.92098 11.4697 5.78033C11.329 5.63968 11.25 5.44891 11.25 5.25C11.25 5.05109 11.329 4.86032 11.4697 4.71967C11.6103 4.57902 11.8011 4.5 12 4.5ZM17.25 4.5H14.115C13.8 3.63 12.975 3 12 3C11.025 3 10.2 3.63 9.885 4.5H6.75C6.35217 4.5 5.97064 4.65803 5.68934 4.93934C5.40804 5.22064 5.25 5.60217 5.25 6V18C5.25 18.3978 5.40804 18.7794 5.68934 19.0607C5.97064 19.342 6.35217 19.5 6.75 19.5H17.25C17.6478 19.5 18.0294 19.342 18.3107 19.0607C18.592 18.7794 18.75 18.3978 18.75 18V6C18.75 5.60217 18.592 5.22064 18.3107 4.93934C18.0294 4.65803 17.6478 4.5 17.25 4.5Z" fill="currentColor"/>`},chargePlus:{viewBox:`0 0 24 24`,body:`<path d="M21.5 6.499V8.499H17.5V12.5H15.5V8.499H11.5V6.499H15.5V2.5H17.5V6.499H21.5Z" fill="currentColor"/> <path d="M5.81933 19.3202C5.68199 19.7265 5.30087 20 4.87199 20H4.42109C3.72897 20 3.24613 19.3138 3.47981 18.6624L6.70812 9.66236C6.85064 9.26504 7.22729 9 7.6494 9H9.35585C9.77812 9 10.1549 9.26525 10.2973 9.66279L13.521 18.6628C13.7543 19.3142 13.2715 20 12.5796 20H12.1285C11.6994 20 11.3181 19.7262 11.1809 19.3197L8.56171 11.556C8.55318 11.5307 8.52947 11.5137 8.50278 11.5137C8.4761 11.5137 8.45239 11.5307 8.44385 11.556L5.81933 19.3202ZM5.4277 16.584C5.4277 16.0827 5.8341 15.6763 6.33542 15.6763H10.6367C11.138 15.6763 11.5444 16.0827 11.5444 16.584C11.5444 17.0853 11.138 17.4917 10.6367 17.4917H6.33542C5.8341 17.4917 5.4277 17.0853 5.4277 16.584Z" fill="currentColor"/>`},chargeMinus:{viewBox:`0 0 24 24`,body:`<path d="M20 4H12V6H20V4Z" fill="currentColor"/> <path d="M5.81933 19.3202C5.68199 19.7265 5.30087 20 4.87199 20H4.42109C3.72897 20 3.24613 19.3138 3.47981 18.6624L6.70812 9.66236C6.85064 9.26504 7.22729 9 7.6494 9H9.35585C9.77812 9 10.1549 9.26525 10.2973 9.66279L13.521 18.6628C13.7543 19.3142 13.2715 20 12.5796 20H12.1285C11.6994 20 11.3181 19.7262 11.1809 19.3197L8.56171 11.556C8.55318 11.5307 8.52947 11.5137 8.50278 11.5137C8.4761 11.5137 8.45239 11.5307 8.44385 11.556L5.81933 19.3202ZM5.4277 16.584C5.4277 16.0827 5.8341 15.6763 6.33542 15.6763H10.6367C11.138 15.6763 11.5444 16.0827 11.5444 16.584C11.5444 17.0853 11.138 17.4917 10.6367 17.4917H6.33542C5.8341 17.4917 5.4277 17.0853 5.4277 16.584Z" fill="currentColor"/>`},enhancedStereo:{viewBox:`0 0 24 24`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M13.9998 4.5C13.9998 5.63772 13.2398 6.59794 12.1998 6.90068V9.86285C12.8744 10.0661 13.4313 10.5556 13.7392 11.1958L15.5061 10.4381C15.5019 10.3926 15.4998 10.3466 15.4998 10.3C15.4998 9.47157 16.1714 8.8 16.9998 8.8C17.8282 8.8 18.4998 9.47157 18.4998 10.3C18.4998 11.1284 17.8282 11.8 16.9998 11.8C16.6545 11.8 16.3364 11.6833 16.0829 11.4872L13.9994 12.3808C13.9957 12.6821 13.9421 12.971 13.8469 13.2391L15.3232 14.512C15.8019 14.1888 16.3789 14 17 14C18.6569 14 20 15.3431 20 17C20 18.6569 18.6569 20 17 20C15.3431 20 14 18.6569 14 17C14 16.5023 14.1212 16.0329 14.3357 15.6197L12.968 14.4404C12.5558 14.7505 12.0484 14.9334 11.4997 14.9334C11.0699 14.9334 10.6654 14.8212 10.3122 14.6235L7.96239 16.0656C7.98711 16.2067 8 16.3518 8 16.5C8 17.8807 6.88071 19 5.5 19C4.11929 19 3 17.8807 3 16.5C3 15.1193 4.11929 14 5.5 14C6.1963 14 6.82611 14.2847 7.27941 14.7439L9.27167 13.5213C9.09775 13.1688 8.99966 12.7696 8.99966 12.3467C8.99966 11.2076 9.71135 10.2404 10.6998 9.89527V6.86926C9.71142 6.53565 8.99978 5.60093 8.99978 4.5C8.99978 3.11929 10.1191 2 11.4998 2C12.8805 2 13.9998 3.11929 13.9998 4.5ZM10.3654 4.5C10.3654 5.12651 10.8733 5.6344 11.4998 5.6344C12.1263 5.6344 12.6342 5.12651 12.6342 4.5C12.6342 3.87349 12.1263 3.3656 11.4998 3.3656C10.8733 3.3656 10.3654 3.87349 10.3654 4.5ZM11.4997 13.4549C10.9081 13.4549 10.4286 12.9587 10.4286 12.3467C10.4286 11.7347 10.9081 11.2385 11.4997 11.2385C12.0912 11.2385 12.5707 11.7347 12.5707 12.3467C12.5707 12.9587 12.0912 13.4549 11.4997 13.4549ZM4.29848 16.5C4.29848 17.1636 4.83642 17.7015 5.5 17.7015C6.16358 17.7015 6.70152 17.1636 6.70152 16.5C6.70152 15.8364 6.16358 15.2985 5.5 15.2985C4.83642 15.2985 4.29848 15.8364 4.29848 16.5ZM17 18.5134C16.1642 18.5134 15.4866 17.8358 15.4866 17C15.4866 16.1642 16.1642 15.4866 17 15.4866C17.8358 15.4866 18.5134 16.1642 18.5134 17C18.5134 17.8358 17.8358 18.5134 17 18.5134ZM16.9998 10.9993C16.6136 10.9993 16.3005 10.6862 16.3005 10.3C16.3005 9.9138 16.6136 9.60072 16.9998 9.60072C17.386 9.60072 17.6991 9.9138 17.6991 10.3C17.6991 10.6862 17.386 10.9993 16.9998 10.9993Z" fill="currentColor"/>`},sgroup:{viewBox:`0 0 24 24`,body:`<path d="M7 4H4v16h3M17 4h3v16h-3" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/><path d="M8 14l2.5-4 3 4 2.5-4" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/>`},keys:{viewBox:`0 0 24 24`,body:`<path d="M11.25 16.5H12.75V15H11.25V16.5ZM12 4.5C11.0151 4.5 10.0398 4.69399 9.12987 5.0709C8.21993 5.44781 7.39314 6.00026 6.6967 6.6967C5.29018 8.10322 4.5 10.0109 4.5 12C4.5 13.9891 5.29018 15.8968 6.6967 17.3033C7.39314 17.9997 8.21993 18.5522 9.12987 18.9291C10.0398 19.306 11.0151 19.5 12 19.5C13.9891 19.5 15.8968 18.7098 17.3033 17.3033C18.7098 15.8968 19.5 13.9891 19.5 12C19.5 11.0151 19.306 10.0398 18.9291 9.12987C18.5522 8.21993 17.9997 7.39314 17.3033 6.6967C16.6069 6.00026 15.7801 5.44781 14.8701 5.0709C13.9602 4.69399 12.9849 4.5 12 4.5ZM12 18C8.6925 18 6 15.3075 6 12C6 8.6925 8.6925 6 12 6C15.3075 6 18 8.6925 18 12C18 15.3075 15.3075 18 12 18ZM12 7.5C11.2043 7.5 10.4413 7.81607 9.87868 8.37868C9.31607 8.94129 9 9.70435 9 10.5H10.5C10.5 10.1022 10.658 9.72064 10.9393 9.43934C11.2206 9.15803 11.6022 9 12 9C12.3978 9 12.7794 9.15803 13.0607 9.43934C13.342 9.72064 13.5 10.1022 13.5 10.5C13.5 12 11.25 11.8125 11.25 14.25H12.75C12.75 12.5625 15 12.375 15 10.5C15 9.70435 14.6839 8.94129 14.1213 8.37868C13.5587 7.81607 12.7956 7.5 12 7.5Z" fill="currentColor"/>`},settings:{viewBox:`0 0 24 24`,body:`<path d="M12.0001 9C12.7958 9 13.5588 9.31607 14.1214 9.87868C14.684 10.4413 15.0001 11.2043 15.0001 12C15.0001 12.7956 14.684 13.5587 14.1214 14.1213C13.5588 14.6839 12.7958 15 12.0001 15C11.2045 15 10.4414 14.6839 9.8788 14.1213C9.31619 13.5587 9.00012 12.7956 9.00012 12C9.00012 11.2043 9.31619 10.4413 9.8788 9.87868C10.4414 9.31607 11.2045 9 12.0001 9ZM12.0001 10.5C11.6023 10.5 11.2208 10.658 10.9395 10.9393C10.6582 11.2206 10.5001 11.6022 10.5001 12C10.5001 12.3978 10.6582 12.7794 10.9395 13.0607C11.2208 13.342 11.6023 13.5 12.0001 13.5C12.3979 13.5 12.7795 13.342 13.0608 13.0607C13.3421 12.7794 13.5001 12.3978 13.5001 12C13.5001 11.6022 13.3421 11.2206 13.0608 10.9393C12.7795 10.658 12.3979 10.5 12.0001 10.5ZM10.5001 19.5C10.3126 19.5 10.1551 19.365 10.1251 19.185L9.84762 17.1975C9.37512 17.01 8.97012 16.755 8.58012 16.455L6.71262 17.2125C6.54762 17.2725 6.34512 17.2125 6.25512 17.0475L4.75512 14.4525C4.65762 14.2875 4.70262 14.085 4.84512 13.9725L6.42762 12.7275L6.37512 12L6.42762 11.25L4.84512 10.0275C4.70262 9.915 4.65762 9.7125 4.75512 9.5475L6.25512 6.9525C6.34512 6.7875 6.54762 6.72 6.71262 6.7875L8.58012 7.5375C8.97012 7.245 9.37512 6.99 9.84762 6.8025L10.1251 4.815C10.1551 4.635 10.3126 4.5 10.5001 4.5H13.5001C13.6876 4.5 13.8451 4.635 13.8751 4.815L14.1526 6.8025C14.6251 6.99 15.0301 7.245 15.4201 7.5375L17.2876 6.7875C17.4526 6.72 17.6551 6.7875 17.7451 6.9525L19.2451 9.5475C19.3426 9.7125 19.2976 9.915 19.1551 10.0275L17.5726 11.25L17.6251 12L17.5726 12.75L19.1551 13.9725C19.2976 14.085 19.3426 14.2875 19.2451 14.4525L17.7451 17.0475C17.6551 17.2125 17.4526 17.28 17.2876 17.2125L15.4201 16.4625C15.0301 16.755 14.6251 17.01 14.1526 17.1975L13.8751 19.185C13.8451 19.365 13.6876 19.5 13.5001 19.5H10.5001ZM11.4376 6L11.1601 7.9575C10.2601 8.145 9.46512 8.625 8.88762 9.2925L7.08012 8.5125L6.51762 9.4875L8.10012 10.65C7.80012 11.5275 7.80012 12.48 8.10012 13.35L6.51012 14.52L7.07262 15.495L8.89512 14.715C9.47262 15.375 10.2601 15.855 11.1526 16.035L11.4301 18H12.5701L12.8476 16.0425C13.7401 15.855 14.5276 15.375 15.1051 14.715L16.9276 15.495L17.4901 14.52L15.9001 13.3575C16.2001 12.48 16.2001 11.5275 15.9001 10.65L17.4826 9.4875L16.9201 8.5125L15.1126 9.2925C14.5351 8.625 13.7401 8.145 12.8401 7.965L12.5626 6H11.4376Z" fill="currentColor"/>`},scrollBack:{viewBox:`0 0 24 24`,body:`<path d="M15 5l-7 7 7 7" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/>`},scrollForward:{viewBox:`0 0 24 24`,body:`<path d="M9 5l7 7-7 7" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/>`},scrollUp:{viewBox:`0 0 24 24`,body:`<path d="M5 15l7-7 7 7" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/>`},scrollDown:{viewBox:`0 0 24 24`,body:`<path d="M5 9l7 7 7-7" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/>`},dative:{viewBox:`0 0 48 48`,body:`<path d="M38.2004 11.2782L32.6464 23.0382C32.4124 23.5322 31.7324 23.5962 31.3004 23.1642L29.6504 21.5142L14.3724 36.5682C14.1886 36.7671 13.9664 36.9266 13.7192 37.0371C13.4719 37.1476 13.2048 37.2067 12.9341 37.2108C12.6633 37.215 12.3945 37.1641 12.144 37.0613C11.8935 36.9585 11.6664 36.8059 11.4766 36.6127C11.2869 36.4195 11.1383 36.1898 11.0399 35.9375C10.9416 35.6852 10.8955 35.4156 10.9044 35.1449C10.9134 34.8743 10.9772 34.6083 11.092 34.363C11.2068 34.1178 11.3703 33.8984 11.5724 33.7182L26.8244 18.6862L25.0964 16.9582C24.6644 16.5242 24.7284 15.8462 25.2224 15.6122L36.9824 10.0582C37.6904 9.7202 38.5324 10.5622 38.2004 11.2782Z" fill="currentColor"/>`},hydrogens:{viewBox:`0 0 24 24`,body:`<path d="M1.9794 6V1.63636H2.63778V3.53267H4.81321V1.63636H5.47372V6H4.81321V4.0973H2.63778V6H1.9794Z" fill="currentColor" fill-rule="evenodd" clip-rule="evenodd"/><path d="M9.9794 6V1.63636H10.6378V3.53267H12.8132V1.63636H13.4737V6H12.8132V4.0973H10.6378V6H9.9794Z" fill="currentColor" fill-rule="evenodd" clip-rule="evenodd"/><path d="M1.9794 14V9.63636H2.63778V11.5327H4.81321V9.63636H5.47372V14H4.81321V12.0973H2.63778V14H1.9794Z" fill="currentColor" fill-rule="evenodd" clip-rule="evenodd"/><path d="M9.9794 14V9.63636H10.6378V11.5327H12.8132V9.63636H13.4737V14H12.8132V12.0973H10.6378V14H9.9794Z" fill="currentColor" fill-rule="evenodd" clip-rule="evenodd"/><path d="M15.4855 6.8855L19 10.4V3.37097L15.4855 6.8855Z" fill="currentColor" fill-rule="evenodd" clip-rule="evenodd"/><path d="M15.4855 17.9855L19 21.5V14.4709L15.4855 17.9855Z" fill="currentColor" fill-rule="evenodd" clip-rule="evenodd"/><circle cx="8" cy="20" r="2" fill="currentColor"/><circle cx="8" cy="8" r="2" fill="currentColor"/><path d="M19 18.7C21.2091 18.7 23 16.9091 23 14.7V10C23 7.79086 21.2091 6 19 6V7.5C20.3807 7.5 21.5 8.61929 21.5 10V14.7C21.5 16.0807 20.3807 17.2 19 17.2V18.7Z" fill="currentColor" fill-rule="evenodd" clip-rule="evenodd"/>`},corner:{viewBox:`0 0 8 8`,body:`<path d="M8 0V8H0Z" fill="currentColor"/>`},menuMark:{viewBox:`0 0 7 4`,body:`<path d="M1 0.8L3.5 3.2L6 0.8" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round" stroke-linejoin="round"/>`},aromatize:{viewBox:`0 0 24 24`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M8.61309 7.54726L4.39649 3.33066L3.33583 4.39132L7.55855 8.61403L4.58825 8.61403L4.58825 10.114L10.1131 10.114L10.1131 4.58919L8.61309 4.58919L8.61309 7.54726ZM8.04762 16.5941V11.7589H10.0476V15.4394L14.0104 17.7273L17.9731 15.4394V10.8636L14.0104 8.57575L11.7432 9.88472V7.57532L14.0104 6.26635L19.9731 9.70894V16.5941L14.0104 20.0367L8.04762 16.5941ZM12.9201 13.1516C12.9201 12.5494 13.4082 12.0613 14.0104 12.0613C14.6126 12.0613 15.1007 12.5494 15.1007 13.1516C15.1007 13.7537 14.6126 14.2419 14.0104 14.2419C13.4082 14.2419 12.9201 13.7537 12.9201 13.1516ZM14.0104 10.5613C12.5798 10.5613 11.4201 11.721 11.4201 13.1516C11.4201 14.5821 12.5798 15.7419 14.0104 15.7419C15.441 15.7419 16.6007 14.5821 16.6007 13.1516C16.6007 11.721 15.441 10.5613 14.0104 10.5613Z" fill="currentColor"/>`},dearomatize:{viewBox:`0 0 24 24`,body:`<path fill-rule="evenodd" clip-rule="evenodd" d="M12.2568 14.1153L9.98962 15.4242L6.02687 13.1363V8.56055L9.98962 6.27266L13.9524 8.56055V12.2411H15.9524V7.40585L9.98962 3.96326L4.02687 7.40585V14.291L9.98962 17.7336L12.2568 16.4247V14.1153ZM13.8869 13.886L19.4117 13.886V15.386L16.4415 15.386L20.6642 19.6087L19.6035 20.6693L15.3869 16.4527L15.3869 19.4108H13.8869L13.8869 13.886ZM10.2385 8.47162V7.22958L9.85857 7.46082C8.99387 7.98705 8.12848 8.51438 7.2624 9.0428L6.91243 9.25633L7.93143 9.87739L10.2385 8.47162ZM12.8365 12.8766V8.81261L11.7967 9.44518V12.2441L12.8365 12.8766ZM10.2385 13.2259L7.9307 11.8074L6.91614 12.4381L10.2385 14.4683V13.2259Z" fill="currentColor"/>`},removeStereo:{viewBox:`0 0 24 24`,body:`<path d="M3 21L11.5 7.5L15.5 11.5Z" fill="currentColor" stroke="currentColor" stroke-width="1" stroke-linejoin="round"/><path d="M15 15l5.5 5.5M20.5 15l-5.5 5.5" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"/>`},recomputeWedges:{viewBox:`0 0 24 24`,body:`<path d="M3 21L10 12L13 15Z" fill="currentColor" stroke="currentColor" stroke-width="1" stroke-linejoin="round"/><path d="M21 7a4 4 0 1 1-2-3.46" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"/><path d="M18.8 1.3L19.5 3.8L17 4.5" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/>`}},sc=new Set([`path`,`line`,`polyline`,`polygon`,`rect`,`circle`,`ellipse`,`text`]),cc=[`d`,`points`,`x`,`y`,`x1`,`y1`,`x2`,`y2`,`cx`,`cy`,`r`,`rx`,`ry`,`width`,`height`,`fill`,`stroke`,`stroke-width`,`opacity`,`style`,`transform`];function lc(e){let t=e.getAttribute(`pointer-events`);return t!==null&&t!==`none`||/\b(atom|bond)-selector\b/.test(e.getAttribute(`class`)??``)}var uc=e=>/\bCIP_Code\b/.test(e.getAttribute(`class`)??``);function dc(e){let t=[];for(let n of Array.from(e.querySelectorAll(`*`))){let e=n.localName;!sc.has(e)||lc(n)||uc(n)||t.push([e,...cc.map(e=>n.getAttribute(e)??``),Ht(n),n.textContent??``].join(`|`))}return t.join(`
`)}function fc(e){let t=e.getAttribute(`data-id`);return t===null?null:e.hasAttribute(`data-atom`)?`atom:${t}`:e.hasAttribute(`data-bond`)?`bond:${t}`:null}function pc(e,t){return[...e].filter(([e,n])=>t.get(e)!==n).map(([e])=>e)}function mc(e){let t=new Map;if(!e)return t;for(let n of Array.from(e.querySelectorAll(`g[data-atom], g[data-bond]`))){let e=fc(n);e&&t.set(e,{el:n,sig:dc(n)})}return t}var hc=new WeakMap;function gc(e){if(!e)return{groups:new Map,elements:new Map,added:new Map};let t=hc.get(e);if(!t){let n=mc(e);t={groups:new Map([...n].map(([e,t])=>[e,t.sig])),elements:new Map([...n].map(([e,t])=>[e,t.el])),added:new Map(Array.from(e.querySelectorAll(`[data-added-for]`)).map(e=>[dc(e),e]))},hc.set(e,t)}return t}var _c=/-?\d+(?:\.\d+)?(?:e[-+]?\d+)?/gi;function vc(e){if(/[^MLml\d\s.,eE+-]/.test(e)||/[ml]/.test(e.replace(/e[-+]?\d/gi,``)))return[];let t=[],n=null;for(let r of e.match(/[ML][^ML]*/g)??[]){let e=r.match(_c)?.map(Number)??[],i=[];for(let t=0;t+1<e.length;t+=2)i.push({x:e[t],y:e[t+1]});let a=r[0]===`M`?null:n;for(let e of i)a&&t.push([a,e]),a=e;n=a}return t}function yc(e){let t=Qr(e),n=/(?:^|;)\s*stroke\s*:\s*([^;]+)/.exec(t)?.[1]?.trim()??e.getAttribute(`stroke`)??`none`;return n!==`none`&&n!==`transparent`}function bc(e,t){let n=[];for(let r of Array.from(e.querySelectorAll(`path`)))if(!lc(r)&&yc(r))for(let[e,i]of vc(r.getAttribute(`d`)??``))Math.hypot(i.x-e.x,i.y-e.y)>.05&&n.push({a:e,b:i,group:t});return n}function xc(e,t){let n=Math.hypot(e.b.x-e.a.x,e.b.y-e.a.y)||1,[r,i]=[(e.b.x-e.a.x)/n,(e.b.y-e.a.y)/n],[a,o]=[t.x-e.a.x,t.y-e.a.y];return{along:a*r+o*i,off:Math.abs(a*i-o*r)}}function Sc(e,t,n=.5){let[r,i]=[Math.atan2(e.b.y-e.a.y,e.b.x-e.a.x),Math.atan2(t.b.y-t.a.y,t.b.x-t.a.x)],a=Math.abs((r-i)*180/Math.PI%180+180)%180;if(Math.min(a,180-a)>2)return!1;let[o,s]=[xc(e,t.a),xc(e,t.b)];if(o.off>n||s.off>n)return!1;let c=Math.hypot(e.b.x-e.a.x,e.b.y-e.a.y),[l,u]=[Math.min(o.along,s.along),Math.max(o.along,s.along)];return l<=c+n&&u>=-n}var Cc=1,wc=32;function Tc(e){let[t,n]=[Math.floor((Math.min(e.a.x,e.b.x)-Cc)/wc),Math.floor((Math.max(e.a.x,e.b.x)+Cc)/wc)],[r,i]=[Math.floor((Math.min(e.a.y,e.b.y)-Cc)/wc),Math.floor((Math.max(e.a.y,e.b.y)+Cc)/wc)],a=[];for(let e=t;e<=n;e++)for(let t=r;t<=i;t++)a.push(`${e},${t}`);return a}function Ec(e,t){let n=new Set(t),r=new Map,i=new Map;for(let t of e){let e=r.get(t.group);if(e?e.push(t):r.set(t.group,[t]),!n.has(t.group))for(let e of Tc(t))(i.get(e)??i.set(e,[]).get(e)).push(t)}let a=e.filter(e=>n.has(e.group));for(;a.length;){let e=a.pop();for(let t of Tc(e))for(let o of i.get(t)??[])!n.has(o.group)&&Sc(e,o)&&(n.add(o.group),a.push(...r.get(o.group)??[]))}return n}var Dc=/<g\b[^>]*\bdata-(?:atom|bond|added-for)=[^>]*>[\s\S]*?<\/g>/g,Oc=/^<g\b[^>]*>/,kc=/\bdata-id=(?:'([^']*)'|"([^"]*)")/,Ac=/<path\b[^>]*\bCIP_Code\b[^>]*>(?:\s*<\/path>)?/g;function jc(e){let t=[],n=[];for(let r of e.matchAll(Dc)){let e=r[0],i=Oc.exec(e)?.[0]??``,a=e.replace(Ac,``);if(/\bdata-added-for=/.test(i)){n.push({key:null,added:!0,text:e,same:a});continue}let o=kc.exec(i);if(!o)continue;let s=/\bdata-atom=/.test(i)?`atom`:`bond`;t.push({key:`${s}:${o[1]??o[2]}`,added:!1,text:e,same:a})}return[...t,...n]}function Mc(e,t){let n=/<svg\b[^>]*>/.exec(e)?.[0]??``,r=RegExp(`\\s${t}=(?:'([^']*)'|"([^"]*)")`).exec(n);return r?r[1]??r[2]??null:null}var Nc=new WeakMap;function Pc(e,t){let n=Nc.get(e);if(!n){let r=jc(t),i=Array.from(e.querySelectorAll(`[data-added-for]`)),a=r.filter(e=>e.added);n={keyed:new Map(r.filter(e=>!e.added).map(e=>[e.key,e.same])),added:new Map(a.length===i.length?a.map((e,t)=>[e.same,i[t]]):[])},Nc.set(e,n)}return n}var Fc=new WeakMap;function Ic(e,t){let n=Fc.get(e);return n||(n=bc(e,0).map(e=>[e.a,e.b]),Fc.set(e,n)),n.map(([e,n])=>({a:e,b:n,group:t}))}function Lc(e,t,n,r){let i=jc(e),a=gc(t),o=Pc(t,n);if(o.keyed.size!==a.groups.size)return null;let s=new Map,c=e=>{let t=r(e.map(e=>i[e].text));return t.length===e.length&&(e.forEach((e,n)=>s.set(e,t[n])),!0)},l=i.map((e,t)=>e.added?o.added.has(e.same)?-1:t:o.keyed.get(e.key)===e.same?-1:t).filter(e=>e>=0);if(!c(l))return null;let u=new Set,d=[],f=new Set;for(let e of l){let t=i[e],n=dc(s.get(e));t.added?(f.add(n),a.added.has(n)||u.add(e)):a.groups.get(t.key)!==n&&(u.add(e),d.push(t.key))}let p=[];i.forEach((e,t)=>{let n=s.get(t);if(n){p.push(...bc(n,t));return}let r=e.added?o.added.get(e.same):a.elements.get(e.key);r&&(p.push(...Ic(r,t)),e.added&&f.add(dc(r)))});let m=Ec(p,u);if(!c([...m].filter(e=>!s.has(e))))return null;let h=new Set(i.filter(e=>!e.added).map(e=>e.key)),g=[...a.groups.keys()].filter(e=>!h.has(e)),_=[...d,...g].map(e=>a.elements.get(e)).filter(e=>e!==void 0);for(let[e,t]of a.added)f.has(e)||_.push(t);return{parts:[...m].sort((e,t)=>e-t).map(e=>s.get(e)),replaced:_}}function Rc(e,t){let n=mc(e),r=gc(t),i=pc(new Map([...n].map(([e,t])=>[e,t.sig])),r.groups),a=new Set(i),o=[...n.values()].map(e=>e.el),s=Array.from(e.querySelectorAll(`[data-added-for]`)),c=new Set;[...n.keys()].forEach((e,t)=>a.has(e)&&c.add(t));let l=new Set(s.map(dc));s.forEach((e,t)=>!r.added.has(dc(e))&&c.add(o.length+t));let u=[...o,...s],d=Ec(u.flatMap((e,t)=>bc(e,t)),c),f=[...r.groups.keys()].filter(e=>!n.has(e)),p=[...i,...f].map(e=>r.elements.get(e)).filter(e=>e!==void 0);for(let[e,t]of r.added)l.has(e)||p.push(t);return{parts:u.filter((e,t)=>d.has(t)),replaced:p}}function zc(e){return[...Array.from(e.querySelectorAll(`g[data-atom], g[data-bond]`)).filter(e=>fc(e)!==null),...Array.from(e.querySelectorAll(`[data-added-for]`))]}var Bc=`http://www.w3.org/2000/svg`,Vc=.6,Hc=`#1c65c0`,Uc=`#8cc4ff`,Wc=`#1c65c0`,Gc=`#ffc54d`,Kc=`#a15c00`,qc=`#d1242f`,Jc=`#7ddc8c`,Yc=`#1b7a3a`,Xc=`0.3`,Zc=`#5c6470`,Qc=`#1c65c0`,$c=(10+2*ct)/20,el=`#1c65c0`,tl=.6,nl=e=>Math.max(2,ct/20*e);function rl(e,t,n){let r=nl(t),i=r+n/2,[a,o]=[e.width+2*i,e.height+2*i];return{x:e.x-i,y:e.y-i,width:a,height:o,rx:Math.min(2*r+n/2,a/2,o/2)}}function B(e,t,n={}){let r=e.createElementNS(Bc,t);for(let[e,t]of Object.entries(n))r.setAttribute(e,String(t));return r}function il(e,t,n){let r=t.bondPx,{fill:i,stroke:a,outline:o}=n,s=[];t.bonds.forEach(([c,l],u)=>{let d=Math.hypot(l.x-c.x,l.y-c.y),f=.26083333333333336*r+o/2,[p,m]=[Math.max(.85*d,6),2*f],h=Math.min(m/3,p/2),g=`translate(${(c.x+l.x)/2} ${(c.y+l.y)/2}) rotate(${Math.atan2(l.y-c.y,l.x-c.x)*180/Math.PI})`;s.push(B(e,`rect`,{x:-p/2,y:-m/2,width:p,height:m,rx:h,transform:g,"data-band":``,...i}));let _=p/2-h;if(_>0&&!t.short?.[u]&&!n.bandsFilled){let t=`M ${-_} ${-f} L ${_} ${-f} M ${-_} ${f} L ${_} ${f}`;s.push(B(e,`path`,{d:t,fill:`none`,"stroke-linecap":`butt`,transform:g,...a}))}});for(let c of t.atoms){let t=n.grow??0;if(c.label){let n=rl(c.label,r,o);s.push(B(e,`rect`,{x:n.x-t,y:n.y-t,width:n.width+2*t,height:n.height+2*t,rx:n.rx+t,...i,...a}))}else s.push(B(e,`circle`,{cx:c.x,cy:c.y,r:Math.max(.4*r,5)+t,...i,...a}))}return s}var al=`#1c65c0`,ol=4,sl=6;function cl(e,t,n,r=1){let i=Math.min(Math.max(r,.5),220/Math.max(t.width,t.height,1)),[a,o]=[t.width*i+12,t.height*i+12],s=e.y-sl-o,c=e.y+e.height+sl,l=s>=ol||c+o>n.height-ol?Math.max(ol,s):c;return{x:Math.max(ol,Math.min(e.x+e.width/2-a/2,n.width-ol-a)),y:l,width:a,height:o}}var ll=16;function ul(e,t,n=`; `){let r=[],i=``;for(let a of e)i&&[...i].length+n.length+[...a].length>t?(r.push(i+n),i=a):i=i?i+n+a:a;return i&&r.push(i),r}var dl=e=>Math.max(12,Math.min(36,Math.floor((e-8-10)/hl)));function fl(e,t,n){let r=10+hl*Math.max(0,...t.map(e=>[...e.trimEnd()].length)),i=ll*Math.max(1,t.length)+4,a=e.y-sl-i,o=e.y+e.height+sl,s=a>=ol||o+i>n.height-ol?Math.max(ol,a):o;return{x:Math.max(ol,Math.min(e.x+e.width/2-r/2,n.width-ol-r)),y:s,width:r,height:i}}var pl=`#1f2328`,ml=12,hl=7,gl=class{element;chrome;#e;#t;#n;#r;#i;#a;#o;#s;#c;#l;#u;constructor(e){this.element=B(e,`svg`,{class:`overlay`,"aria-hidden":`true`,focusable:`false`}),this.#e=B(e,`g`,{"data-testid":`overlay.problems`,display:`none`}),this.#t=B(e,`g`,{"data-testid":`overlay.selection`,display:`none`}),this.#n=B(e,`g`,{"data-testid":`overlay.area`,display:`none`}),this.#r=B(e,`g`,{"data-testid":`overlay.hover`,display:`none`}),this.#i=B(e,`g`,{"data-testid":`overlay.preview`,display:`none`}),this.#a=B(e,`g`,{"data-testid":`overlay.chain-count`,display:`none`}),this.element.append(this.#e,this.#t,this.#n,this.#r,this.#i,this.#a),this.chrome=B(e,`svg`,{class:`chrome`,"aria-hidden":`true`,focusable:`false`}),this.#s=B(e,`g`,{"data-testid":`overlay.rotate-centre`,display:`none`}),this.#o=B(e,`g`,{"data-testid":`overlay.rotate-handle`,display:`none`}),this.#c=B(e,`g`,{"data-testid":`overlay.rotate-angle`,display:`none`}),this.#l=B(e,`g`,{"data-testid":`overlay.abbreviation-preview`,display:`none`}),this.#u=B(e,`g`,{"data-testid":`overlay.query-tooltip`,display:`none`}),this.chrome.append(this.#s,this.#o,this.#c,this.#l,this.#u)}showQueryTip(e,t,n,r=`; `){let i=this.chrome.ownerDocument,a=ul(t,dl(n.width),r),o=fl(e,a,n),s=B(i,`text`,{x:o.x+5,y:o.y+2,fill:pl,"font-size":ml,"font-family":`system-ui, sans-serif`,"text-anchor":`start`,"dominant-baseline":`central`});a.forEach((e,t)=>s.append(Object.assign(B(i,`tspan`,{x:o.x+5,y:o.y+2+ll*(t+.5)}),{textContent:e}))),this.#u.replaceChildren(B(i,`rect`,{x:o.x,y:o.y,width:o.width,height:o.height,rx:3,fill:`#ffffff`,stroke:al,"stroke-width":1}),s),this.#u.setAttribute(`display`,`inline`)}hideQueryTip(){this.#d(this.#u)}get queryTip(){return this.#u.getAttribute(`display`)===`none`?null:this.#u.textContent??null}showAbbreviation(e,t,n,r=1){let i=this.chrome.ownerDocument,a=i.createElement(`template`),o=t.indexOf(`<svg`);a.innerHTML=Qt(o>=0?t.slice(o):t);let s=a.content.firstElementChild;if(!s||s.localName!==`svg`){this.hideAbbreviation();return}let c=i.adoptNode(s);Lr(i);let l=cl(e,{width:Number.parseFloat(c.getAttribute(`width`)??``)||0,height:Number.parseFloat(c.getAttribute(`height`)??``)||0},n,r);c.setAttribute(`x`,String(l.x+6)),c.setAttribute(`y`,String(l.y+6)),c.setAttribute(`width`,String(l.width-12)),c.setAttribute(`height`,String(l.height-12)),c.setAttribute(`overflow`,`visible`),this.#l.replaceChildren(B(i,`rect`,{x:l.x,y:l.y,width:l.width,height:l.height,rx:4,fill:`#ffffff`,stroke:al,"stroke-width":1.5}),c),this.#l.setAttribute(`display`,`inline`)}hideAbbreviation(){this.#d(this.#l)}get abbreviationShown(){return this.#l.getAttribute(`display`)!==`none`}showChrome(e){let t=this.chrome.ownerDocument;if(e.handle){let{at:n,stem:r}=e.handle,i=oc.rotate,[,,a]=i.viewBox.split(` `).map(Number),o=1.4*12,s=[];r&&s.push(B(t,`line`,{x1:r[0].x,y1:r[0].y,x2:r[1].x,y2:r[1].y,stroke:al,"stroke-width":1.5})),s.push(B(t,`circle`,{cx:n.x,cy:n.y,r:11,fill:`#ffffff`,stroke:al,"stroke-width":2}));let c=B(t,`g`,{transform:`translate(${n.x-o/2} ${n.y-o/2}) scale(${o/a})`,fill:al,color:al});c.innerHTML=i.body,s.push(c),this.#o.replaceChildren(...s),this.#o.setAttribute(`display`,`inline`)}else this.#d(this.#o);if(e.centre){let{x:n,y:r}=e.centre;this.#s.replaceChildren(B(t,`path`,{d:`M ${n-10} ${r} L ${n+10} ${r} M ${n} ${r-10} L ${n} ${r+10}`,stroke:al,"stroke-width":1.5,fill:`none`}),B(t,`circle`,{cx:n,cy:r,r:5,fill:`#ffffff`,"fill-opacity":.6,stroke:al,"stroke-width":1.5})),this.#s.setAttribute(`display`,`inline`)}else this.#d(this.#s);if(e.angle){let{at:n,text:r,leftward:i}=e.angle,a=8+hl*r.length,o=i?n.x-a:n.x;this.#c.replaceChildren(B(t,`rect`,{x:o,y:n.y-10,width:a,height:20,rx:3,fill:`#ffffff`,stroke:al,"stroke-width":1}),Object.assign(B(t,`text`,{x:o+a/2,y:n.y,fill:pl,"font-size":ml,"font-family":`system-ui, sans-serif`,"text-anchor":`middle`,"dominant-baseline":`central`}),{textContent:r})),this.#c.setAttribute(`display`,`inline`)}else this.#d(this.#c)}hideChrome(){this.#d(this.#o),this.#d(this.#s),this.#d(this.#c)}#d(e){(e.getAttribute(`display`)!==`none`||e.firstChild)&&(e.replaceChildren(),e.setAttribute(`display`,`none`))}showHover(e){this.#f+=1;let t=this.element.ownerDocument,n=!!e.selected&&!e.warn,r=e.quiet&&!e.warn?1:Math.max(1,.06*e.bondPx),i=il(t,e,{fill:e.warn?{fill:Gc,"fill-opacity":`0.55`}:n?{fill:Jc,"fill-opacity":e.quiet?Xc:`0.6`}:{fill:Uc,"fill-opacity":e.quiet?`0.2`:`0.45`},stroke:e.warn?{stroke:Kc,"stroke-width":r}:{stroke:Wc,"stroke-width":r,...e.quiet?{"stroke-dasharray":`2 2`}:{}},outline:r});this.#r.replaceChildren(...i),this.#r.setAttribute(`display`,i.length?`inline`:`none`),this.#r.setAttribute(`data-quiet`,String(!!e.quiet)),this.#r.setAttribute(`data-warn`,String(!!e.warn)),this.#r.setAttribute(`data-selected`,String(n))}get hoverVersion(){return this.#f}#f=0;hideHover(){this.#f+=1,this.#r.replaceChildren(),this.#r.setAttribute(`display`,`none`),this.#r.removeAttribute(`data-quiet`),this.#r.removeAttribute(`data-warn`),this.#r.removeAttribute(`data-selected`)}showSelection(e){let t=Math.max(1,.06*e.bondPx),n=il(this.element.ownerDocument,e,{fill:{fill:Jc,"fill-opacity":`0.45`},stroke:{stroke:Yc,"stroke-width":t},outline:t,bandsFilled:!0});this.#t.replaceChildren(...n),this.#t.setAttribute(`display`,n.length?`inline`:`none`),this.#t.removeAttribute(`transform`)}hideSelection(){this.#t.replaceChildren(),this.#t.setAttribute(`display`,`none`),this.#t.removeAttribute(`transform`)}showProblems(e){let t=Math.max(1.5,.06*e.bondPx),n=il(this.element.ownerDocument,{...e,bonds:[]},{fill:{fill:`none`},stroke:{stroke:qc,"stroke-width":t,"stroke-dasharray":`${2*t} ${1.5*t}`},outline:t,grow:3});this.#e.replaceChildren(...n),this.#e.setAttribute(`display`,n.length?`inline`:`none`)}hideProblems(){(this.#e.getAttribute(`display`)!==`none`||this.#e.firstChild)&&(this.#e.replaceChildren(),this.#e.setAttribute(`display`,`none`))}moveSelection(e,t){_l(this.#t,e,t)}showArea(e){let t=this.element.ownerDocument,n={fill:Qc,"fill-opacity":`0.06`,stroke:Zc,"stroke-width":1,"stroke-dasharray":`4 3`},r;if(e.kind===`rect`){let[i,a]=[Math.min(e.from.x,e.to.x),Math.min(e.from.y,e.to.y)];r=B(t,`rect`,{x:i,y:a,width:Math.abs(e.to.x-e.from.x),height:Math.abs(e.to.y-e.from.y),...n})}else r=B(t,`path`,{d:`${e.path.map((e,t)=>`${t?`L`:`M`} ${e.x} ${e.y}`).join(` `)} Z`,"stroke-linejoin":`round`,...n});this.#n.replaceChildren(r),this.#n.setAttribute(`display`,`inline`)}hideArea(){this.#n.replaceChildren(),this.#n.setAttribute(`display`,`none`)}showPreview(e,t,n,r){let i=this.element.ownerDocument,a=Math.max(1,n/12);this.#i.replaceChildren(B(i,`line`,{x1:e.x,y1:e.y,x2:t.x,y2:t.y,stroke:el,"stroke-width":a,"stroke-linecap":`round`}),r?this.#p(i,t,r,n):B(i,`circle`,{cx:t.x,cy:t.y,r:Math.max(2,n/8),fill:el})),this.#i.setAttribute(`display`,`inline`),this.#i.removeAttribute(`transform`)}showLabel(e,t,n){this.#i.replaceChildren(this.#p(this.element.ownerDocument,e,t,n)),this.#i.setAttribute(`display`,`inline`),this.#i.removeAttribute(`transform`)}#p(e,t,n,r){let i=Math.max(11,.7*r),a=Math.max(9,.36*i*Math.max(1.6,n.symbol.length+.6)),o=B(e,`g`,{opacity:tl});return o.append(B(e,`circle`,{cx:t.x,cy:t.y,r:a,fill:`#ffffff`,stroke:n.colour,"stroke-width":Math.max(1,r/16)}),Object.assign(B(e,`text`,{x:t.x,y:t.y,fill:n.colour,"font-size":i,"font-weight":600,"font-family":`system-ui, sans-serif`,"text-anchor":`middle`,"dominant-baseline":`central`}),{textContent:n.symbol})),o}moveWith(e){for(let t of[this.element,this.chrome])e?t.setAttribute(`viewBox`,e.join(` `)):t.removeAttribute(`viewBox`)}hidePreview(){this.#i.replaceChildren(),this.#i.setAttribute(`display`,`none`),this.#i.removeAttribute(`transform`)}showDrawn(e,t,n=Vc){let r=this.element.ownerDocument,i=B(r,`svg`,{x:0,y:0,width:`100%`,height:`100%`,overflow:`visible`});t.viewBox&&i.setAttribute(`viewBox`,t.viewBox),t.preserveAspectRatio&&i.setAttribute(`preserveAspectRatio`,t.preserveAspectRatio);let a=B(r,`g`,{opacity:n});for(let t of e)a.append(vl(r,t));i.append(a),this.#i.replaceChildren(i),this.#i.setAttribute(`display`,`inline`),this.#i.removeAttribute(`transform`)}moveDrawn(e,t){_l(this.#i,e,t)}showCount(e,t,n){let r=this.element.ownerDocument,i=Math.max(12,.6*n),{width:a,height:o}=sa(t,n);this.#a.replaceChildren(B(r,`rect`,{x:e.x-a/2,y:e.y-o/2,width:a,height:o,rx:o/4,fill:`#ffffff`,stroke:Hc,"stroke-width":1}),Object.assign(B(r,`text`,{x:e.x,y:e.y,fill:Hc,"font-size":i,"font-weight":700,"font-family":`system-ui, sans-serif`,"text-anchor":`middle`,"dominant-baseline":`central`}),{textContent:t})),this.#a.setAttribute(`display`,`inline`)}hideCount(){this.#a.replaceChildren(),this.#a.setAttribute(`display`,`none`)}};function _l(e,t,n){t===0&&n===0?e.removeAttribute(`transform`):e.setAttribute(`transform`,`translate(${t} ${n})`)}function vl(e,t){let n=e.importNode(t,!0);for(let e of[n,...Array.from(n.querySelectorAll(`*`))]){if(e!==n&&lc(e)){e.remove();continue}let t=Ht(e);for(let t of e.getAttributeNames())(t.startsWith(`data-`)||t===`class`||t===`id`)&&e.removeAttribute(t);t&&e.setAttribute(`class`,t)}return n}function yl(e,t){let n=t??e.atoms.map(e=>e.index);return[...new Set(n)].filter(t=>e.atoms[t]!==void 0&&!e.hiddenAtoms.has(t)).sort((e,t)=>e-t)}var bl=(e,t)=>{let n=e.atoms[t];return n!==void 0&&!e.hiddenAtoms.has(t)&&(n.element||e.info(t)?.symbol)===`H`},xl=(e,t)=>t.some(t=>(e.info(t)?.hCount??0)>0),Sl=(e,t)=>t.some(t=>bl(e,t)||e.bondsOf(t).some(n=>bl(e,e.other(n,t))));function Cl(e,t,n=!1){let r=yl(e,t);return{add:!n&&xl(e,r),remove:Sl(e,r)}}function wl(e,t,n,r=!1){let i=yl(e,n);return t===`add`?!r&&xl(e,i)?S(`hydrogens.add`,[{op:`addHs`,atoms:i,bondLength:e.docLen}]):null:Sl(e,i)?S(`hydrogens.remove`,[{op:`removeHs`,atoms:i}]):null}var Tl=`single`,El=`double`,Dl=(e,t=90)=>({size:e,vertex:t,orders:Array(e).fill(Tl),onSingle:0,onDouble:0}),Ol=Object.freeze({benzene:{size:6,vertex:90,orders:[El,Tl,El,Tl,El,Tl],onSingle:1,onDouble:0},cyclopentadiene:{size:5,vertex:90,orders:[Tl,El,Tl,El,Tl],onSingle:0,onDouble:1},cyclohexane:Dl(6),cyclopentane:Dl(5),cyclopropane:Dl(3),cyclobutane:Dl(4,45),cycloheptane:Dl(7),cyclooctane:Dl(8,22.5)}),kl=.4,Al=(e,t)=>e/(2*Math.sin(Math.PI/t)),jl=(e,t)=>e/(2*Math.tan(Math.PI/t));function Ml(e,t,n,r){let i=Al(r,t);return Array.from({length:t},(r,a)=>O(e,n+a*360/t,i))}var Nl=e=>e.map((t,n)=>e[(n-1+e.length)%e.length]);function Pl(e,t,n,r){return Ml(t,e.size,e.vertex+n,r)}function Fl(e,t,n,r){let i=Jn(e,t),a={x:(e.x+t.x)/2,y:(e.y+t.y)/2},o=O(a,A(e,t)+r*90,jl(i,n)),s=Al(i,n),c=A(o,e),l=k(A(o,t)-c);return Array.from({length:n},(n,r)=>r===0?{...e}:r===1?{...t}:O(o,c+r*l,s))}function Il(e,t,n,r){let i=Ml(O(e,n,Al(r,t)),t,n+180,r);return i[0]={...e},i}function Ll(e,t,n,r){return Il(O(e,n,r),t,n,r)}function Rl(e,t){return t.length?t.length===1?k(A(t[0],e)):Zt(t.map(t=>A(e,t))):0}function zl(e,t,n){let r=n?Vt(n,e,t):0,i=Jn(e,t)||1;return Math.abs(r)/i<1e-6?1:r>0?-1:1}function Bl(e,t,n){let r=Vt(n,e,t),i=Jn(e,t)||1;return Math.abs(r)/i<1e-9?null:r>0?1:-1}function Vl(e,t,n,r=[],i=[]){let a=e.map((e,t)=>r[t]??-1),o=new Set([...r.filter(e=>e>=0),...i]),s=[];e.forEach((e,r)=>{a[r]>=0||t.forEach((t,i)=>{if(o.has(i))return;let a=Jn(t,e);a<=n&&s.push({k:r,atom:i,d:a})})}),s.sort((e,t)=>e.d-t.d);let c=new Set;for(let e of s)a[e.k]>=0||c.has(e.atom)||(a[e.k]=e.atom,c.add(e.atom));return a}function Hl(e,t,n,r){let i=null;n.forEach((n,a)=>{if(n<0)return;let[o,s]=[t[n].x-e[a].x,t[n].y-e[a].y],c=Math.hypot(o,s);c<=r&&(!i||c<i.d)&&(i={dx:o,dy:s,d:c})});let a=i;return a?e.map(e=>({x:e.x+a.dx,y:e.y+a.dy})):[...e]}function Ul(e,t,n=!1){let r=e.length;return e.map((i,a)=>e[((n?t-a:t+a)%r+r)%r])}var Wl=(e,t)=>e.length===t.length&&e.every((e,n)=>e===t[n]);function Gl(e,t){let n=[];for(let r of[!1,!0])for(let i=0;i<e.length;i++){let a=Ul(e,i,r);n.some(e=>Wl(e.orders,a))||n.push({orders:a,rank:t(a)})}return n}function Kl(e,t){let n=t!==`single`,r=n?e.onDouble:e.onSingle,i=[Ul(e.orders,r),Ul(e.orders,r,!0)];return Gl(e.orders,e=>i.some(t=>Wl(t,e))?0:e[0]===El===n?1:2)}function ql(e){let t=[Ul(e.orders,0),Ul(e.orders,e.size-1,!0)];return Gl(e.orders,e=>+!t.some(t=>Wl(t,e)))}var Jl=(e,t)=>{let[n,r]=[e[t],e[(t+1)%e.length]];return!n.unsaturated&&!r.unsaturated&&n.room>=1&&r.room>=1};function Yl(e,t,n){let r=t.length,i=e=>e.map((e,r)=>!n[r]&&e===El&&Jl(t,r)?El:Tl),a=null;for(let t of e){if(t.orders.length!==r)continue;let e=i(t.orders),n=e.filter(e=>e===El).length;(!a||n>a.count||n===a.count&&t.rank<a.rank)&&(a={orders:e,count:n,rank:t.rank})}return a?.orders??Array(r).fill(Tl)}var Xl=.5,Zl=1e-9;function Ql(e,t){return t?Math.abs(e.x)>=Math.abs(e.y)?{x:e.x,y:0}:{x:0,y:e.y}:{x:e.x,y:e.y}}var V={single:1,aromatic:1.5,double:2,triple:3},$l=e=>Object.prototype.hasOwnProperty.call(V,e),eu=e=>e===`none`||e===`wedge`||e===`hash`,tu=(e,t)=>e<t?`${e} ${t}`:`${t} ${e}`;function nu(e,t,n){let r=e=>t.get(e)??e,i=new Set(t.values()),a=new Map;for(let n of e.bonds)!t.has(n.a)&&!t.has(n.b)&&a.set(tu(n.a,n.b),{units:$l(n.order)?V[n.order]:1,type:n.order,stereo:n.stereo,bond:n.index,add:null});let o=[],s=new Set,c=[];for(let i of[...t.keys()].sort((e,t)=>e-t))for(let t of e.bondsOf(i)){if(s.has(t))continue;s.add(t);let i=e.bonds[t];if(!$l(i.order))return null;let l=V[i.order],[u,d]=[r(i.a),r(i.b)];if(u===d)continue;let f=a.get(tu(u,d));if(f){if(!$l(f.type))return null;l>f.units&&(f.add?(f.add.type=i.order,i.order!==`single`&&(f.add.stereo=`none`)):(o.push({op:`setBondType`,bond:f.bond,type:i.order}),f.stereo!==`none`&&i.order!==`single`&&o.push({op:`setBondStereo`,bond:f.bond,stereo:`none`})),f.units=l,f.type=i.order,i.order!==`single`&&(f.stereo=`none`));continue}if(!eu(i.stereo))return null;let p=i.order===`single`?n?n(t):e.wedgeOf(t):null,m=p&&(p.from===i.a||p.from===i.b)?p:null,[h,g]=m&&m.from===i.b?[d,u]:[u,d],_={op:`addBond`,a:h,b:g,type:i.order,stereo:m?.stereo??`none`};c.push({bond:t,add:_}),a.set(tu(u,d),{units:l,type:i.order,stereo:i.stereo,bond:-1,add:_})}let l=[];for(let{bond:e,add:t}of c)o.push({op:`removeBond`,bond:e-l.filter(t=>t<e).length},t),l.push(e);for(let e of[...t.keys()].sort((e,t)=>t-e))o.push({op:`removeAtom`,atom:e});let u=new Map,d=(e,t)=>{i.has(e)&&u.set(e,(u.get(e)??0)+t)};for(let[e,t]of a)for(let n of e.split(` `).map(Number))d(n,t.units);for(let t of e.bonds)for(let e of[t.a,t.b])d(e,-($l(t.order)?V[t.order]:1));let f=[...i].sort((e,t)=>e-t).find(t=>!Ye(e.info(t),u.get(t)??0))??-1;return{ops:o,room:f<0,blocker:f}}var ru=Object.freeze({atoms:Object.freeze([]),bonds:Object.freeze([])}),iu=new Map,au=Object.freeze({tx:null,targets:ru,refused:!1,blocker:-1,crowded:!1,merged:iu}),ou=(e,t=!1)=>Object.freeze({tx:null,targets:ru,refused:!0,blocker:e,crowded:t,merged:iu}),su=(e,t)=>Math.hypot(e.x-t.x,e.y-t.y),cu=(e,t)=>({x:(e.x+t.x)/2,y:(e.y+t.y)/2});function lu(e,t,n,r){let i=e.docLen,a=new Map,o=new Set,s=[],c=e.bonds.filter(t=>r(t.a)&&r(t.b)&&e.atoms[t.a]&&e.atoms[t.b]),l=[],u=Math.cos(10*Math.PI/180);for(let r of n){let[n,a]=[t.get(r.a),t.get(r.b)];if(!n||!a)continue;let o=su(n,a);if(!(o>0))continue;let s=cu(n,a);for(let t of c){let[c,d]=[e.atoms[t.a].at,e.atoms[t.b].at],f=su(s,cu(c,d));if(f>.32*i)continue;let p=su(c,d);if(!(p>0))continue;let m=o/p;if(m<.8||m>=1.2)continue;let h=((a.x-n.x)*(d.x-c.x)+(a.y-n.y)*(d.y-c.y))/(o*p);Math.abs(h)<u||l.push({a:r.a,b:r.b,onto:t.index,ends:h>=0?[t.a,t.b]:[t.b,t.a],d:f})}}l.sort((e,t)=>e.d-t.d||e.onto-t.onto);for(let e of l)a.has(e.a)||a.has(e.b)||o.has(e.ends[0])||o.has(e.ends[1])||(a.set(e.a,e.ends[0]),a.set(e.b,e.ends[1]),o.add(e.ends[0]).add(e.ends[1]),s.push(e.onto));let d=kl*i,f=[];for(let[n,i]of t){if(a.has(n))continue;let t=-1,o=1/0;for(let n of e.atoms){if(!r(n.index))continue;let e=su(i,n.at);e<=d&&e<o&&([t,o]=[n.index,e])}t>=0&&f.push({key:n,atom:t,d:o})}f.sort((e,t)=>e.d-t.d||e.key-t.key);let p=-1;for(let e of f){if(o.has(e.atom)){p<0&&(p=e.atom);continue}a.set(e.key,e.atom),o.add(e.atom)}return{atoms:a,bonds:s.sort((e,t)=>e-t),crowded:p}}var uu=e=>({atoms:[...new Set(e.atoms.values())].sort((e,t)=>e-t),bonds:e.bonds});function du(e,t,n,r=`selection`,i){if(Math.hypot(n.x,n.y)<Zl||t.some(t=>!e.atoms[t]))return au;let a=e.wholeGroups(t);if(!a.length)return au;let o=new Set(a),s=new Map(a.map(t=>[t,{x:e.atoms[t].at.x+n.x,y:e.atoms[t].at.y+n.y}])),c=new Map([...s].filter(([t])=>!e.contractedGroupOf(t))),l=e.bonds.filter(e=>c.has(e.a)&&c.has(e.b)),u=t=>!o.has(t)&&!e.atoms[t]?.hidden,d=lu(e,c,l,u);if(d.crowded>=0)return ou(d.crowded,!0);let f=null;if(d.atoms.size){if(f=nu(e,d.atoms,i),!f)return ou([...d.atoms.values()][0]);if(!f.room)return ou(f.blocker)}let p=Xl*e.docLen;for(let t of c.keys()){if(d.atoms.has(t))continue;let n=s.get(t),r=e.atoms.filter(e=>u(e.index)&&su(n,e.at)<=p).sort((e,t)=>su(n,e.at)-su(n,t.at)).find(n=>!nu(e,new Map([[t,n.index]]))?.room);if(r)return ou(r.index)}let m=a.filter(e=>!d.atoms.has(e)).map(e=>({op:`setAtomPosition`,atom:e,x:s.get(e).x,y:s.get(e).y}));return f&&m.push(...f.ops),{tx:S(`${r}.${f?`merge`:`move`}`,m),targets:uu(d),refused:!1,blocker:-1,crowded:!1,merged:d.atoms}}function fu(e,t,n,r){let i=t.kind===`bond`?e.bonds[t.index]:void 0,a=t.kind===`atom`?[t.index]:i?[i.a,i.b]:[];return a.length?du(e,a,n,t.kind,r):au}function pu(e,t,n){if(Math.hypot(n.x,n.y)<Zl)return null;let r=[...new Set(t)].filter(t=>e.atoms[t]).sort((e,t)=>e-t);return S(`selection.nudge`,r.map(t=>({op:`setAtomPosition`,atom:t,x:e.atoms[t].at.x+n.x,y:e.atoms[t].at.y+n.y})))}function mu(e,t,n){let r=[...new Set(t)].filter(t=>e.atoms[t]);if(r.length<2)return au;let[i,...a]=r,o=new Map(a.map(e=>[e,i])),s=nu(e,o,n);return s?s.room?{tx:S(`atom.join`,s.ops),targets:{atoms:[i],bonds:[]},refused:!1,blocker:-1,crowded:!1,merged:o}:ou(s.blocker):ou(i)}var hu={doubleBondPosition:()=>`auto`},gu=(e,t)=>e<t?`${e} ${t}`:`${t} ${e}`,_u=(e,t)=>Math.hypot(e.x-t.x,e.y-t.y);function vu(e,t){let n=new Set(t.atoms.filter(t=>e.atoms[t]));for(let r of t.bonds)e.bonds[r]&&n.add(e.bonds[r].a).add(e.bonds[r].b);for(let t of e.wholeGroups([...n]))n.add(t);for(let t of e.withCopies([...n]))n.add(t);let r=e.bonds.filter(e=>n.has(e.a)&&n.has(e.b)).map(e=>e.index);return{atoms:[...n].sort((e,t)=>e-t),bonds:r}}var yu=e=>$l(e)?V[e]:1;function bu(e,t,n,r=hu){if(Math.hypot(n.x,n.y)<1e-9)return au;let{atoms:i,bonds:a}=vu(e,t);if(!i.length)return au;let o=new Map(i.map(t=>[t,{x:e.atoms[t].at.x+n.x,y:e.atoms[t].at.y+n.y}])),s=new Map([...o].filter(([t])=>!e.contractedGroupOf(t)&&!e.hiddenAtoms.has(t))),c=lu(e,s,a.map(t=>e.bonds[t]).filter(e=>s.has(e.a)&&s.has(e.b)),t=>!e.atoms[t]?.hidden);if(c.crowded>=0)return ou(c.crowded,!0);let l=a.filter(t=>c.atoms.has(e.bonds[t].a)||c.atoms.has(e.bonds[t].b));for(let t of l)if(!$l(e.bonds[t].order))return au;let u=t=>a.reduce((n,r)=>n+(e.bonds[r].a===t||e.bonds[r].b===t?yu(e.bonds[r].order):0),0),d=Xl*e.docLen;for(let t of s.keys()){if(c.atoms.has(t))continue;let n=o.get(t),r=e.atoms.filter(e=>!e.hidden&&_u(n,e.at)<=d).sort((e,t)=>_u(n,e.at)-_u(n,t.at)).find(n=>!Ye(e.info(n.index),u(t)));if(r)return ou(r.index)}let f=i.filter(e=>!c.atoms.has(e)),p=f.length?[{op:`duplicate`,atoms:f,dx:n.x,dy:n.y}]:[],m=new Map(f.map((e,t)=>[e,{added:t}]));for(let[e,t]of c.atoms)m.set(e,t);let h=new Set(f),g=e.bonds.filter(e=>h.has(e.a)&&h.has(e.b)).length,_=new Map(e.bonds.map(e=>[gu(e.a,e.b),e])),v=new Map,y=(e,t)=>{typeof e==`number`&&v.set(e,(v.get(e)??0)+t)},b=[],ee=0;for(let t of l){let n=e.bonds[t],i=n.order,[a,o]=[m.get(n.a),m.get(n.b)];if(typeof a==`number`&&typeof o==`number`){if(a===o)continue;let e=_.get(gu(a,o));if(e){let t=$l(e.order)?V[e.order]:1;V[i]>t&&(y(a,V[i]-t),y(o,V[i]-t),p.push({op:`setBondType`,bond:e.index,type:i}),e.stereo!==`none`&&i!==`single`&&p.push({op:`setBondStereo`,bond:e.index,stereo:`none`}));continue}}let s=i===`single`?r.wedge?r.wedge(t):e.wedgeOf(t):null,c=s&&(s.from===n.a||s.from===n.b)?s:null,[l,u]=c&&c.from===n.b?[o,a]:[a,o];p.push({op:`addBond`,a:l,b:u,type:i,stereo:c?.stereo??`none`}),y(a,V[i]),y(o,V[i]);let d=i===`double`?r.doubleBondPosition(t):`auto`;d!==`auto`&&b.push({op:`setDoubleBondPosition`,bond:e.bonds.length+g+ee,pos:d}),ee+=1}p.push(...b);let te=[...v.entries()].sort((e,t)=>e[0]-t[0]).find(([t,n])=>!Ye(e.info(t),n))?.[0];if(te!==void 0)return ou(te);let ne={atoms:[...new Set(c.atoms.values())].sort((e,t)=>e-t),bonds:c.bonds};return{tx:S(c.atoms.size?`selection.copy-merge`:`selection.copy`,p),targets:ne,refused:!1,blocker:-1,crowded:!1,merged:c.atoms}}var xu=Object.freeze([90,-90,120,-120,180]),Su=1e-9,Cu=1e-9;function H(e){let t=e%360;return t<=-180&&(t+=360),t>180&&(t-=360),t}var wu=(e,t)=>Math.atan2(t.y-e.y,t.x-e.x)*180/Math.PI;function Tu(e,t,n){let r=n*Math.PI/180,[i,a]=[Math.cos(r),Math.sin(r)],[o,s]=[e.x-t.x,e.y-t.y];return{x:t.x+o*i-s*a,y:t.y+o*a+s*i}}function Eu(e){if(!e.length)return{x:0,y:0};let[t,n]=[e.map(e=>e.x),e.map(e=>e.y)];return{x:(Math.min(...t)+Math.max(...t))/2,y:(Math.min(...n)+Math.max(...n))/2}}function Du(e,t){return e.bonds.filter(e=>t.has(e.a)!==t.has(e.b)).map(e=>e.index)}function Ou(e,t){let n=[...new Set(t)].filter(t=>e.atoms[t]).sort((e,t)=>e-t);if(!n.length)return null;let r=new Set(n),i=Du(e,r),a=-1;if(i.length===1){let t=e.bonds[i[0]];a=r.has(t.a)?t.b:t.a}else if(i.length>1){let[t,...n]=i.map(t=>e.bonds[t]);a=[t.a,t.b].find(e=>n.every(t=>t.a===e||t.b===e))??-1}return{atoms:n,attachments:i,centre:a>=0?{...e.atoms[a].at}:Eu(n.map(t=>e.atoms[t].at)),atom:a}}function ku(e,t,n){let r=e.atoms[n];if(!r)return null;let i=new Set(t),a=[],o=[];for(let t of e.bondsOf(n)){let s=e.other(t,n);if(s===n||!e.atoms[s])continue;let c=wu(r.at,e.atoms[s].at);i.has(s)?a.push(c):o.push(c)}if(a.length!==1||!o.length)return null;let s=[];if(o.length===1)for(let e of xu)s.push(H(o[0]+e));else for(let e=0;e<o.length;e++)for(let t=e+1;t<o.length;t++){let n=((o[t]-o[e])%360+360)%360,r=o[e]+n/2;Math.abs(n-180)<1e-6?s.push(H(r),H(r+180)):s.push(H(n>180?r:r+180))}return{from:H(a[0]),to:s}}function Au(e,t,n=null){let r=H(e);if(t)return H(Math.round(r/1)*1);if(n){let e=null,t=1/0;for(let i of n.to){let a=Math.abs(H(n.from+r-i));a<=10.000000001&&a<t&&([e,t]=[H(i-n.from),a])}if(e!==null)return e}return H(Math.round(r/15)*15)}function ju(e,t,n,r){let i=H(r);if(Math.abs(i)<Su)return null;let a=[];for(let r of[...new Set(t)].filter(t=>e.atoms[t]).sort((e,t)=>e-t)){let t=e.atoms[r].at,o=Tu(t,n,i);Math.hypot(o.x-t.x,o.y-t.y)<Cu||a.push({op:`setAtomPosition`,atom:r,x:o.x,y:o.y})}return S(`selection.rotate`,a)}var Mu=1e-9;function Nu(e,t,n){let[r,i]=[n.x-t.x,n.y-t.y],a=r*r+i*i;if(!(a>0))return{x:e.x,y:e.y};let o=((e.x-t.x)*r+(e.y-t.y)*i)/a,s={x:t.x+o*r,y:t.y+o*i};return{x:2*s.x-e.x,y:2*s.y-e.y}}var Pu=e=>e===`wedge`?`hash`:`wedge`;function Fu(e,t,n,r,i,a=!1,o=()=>!1){let s=[...new Set(t)].filter(t=>e.atoms[t]).sort((e,t)=>e-t);if(!s.length)return{tx:null,refused:!1,reason:null,attachments:[]};let c=new Set(s),l=Du(e,c);if(l.length>1)return{tx:null,refused:!0,reason:`bonds`,attachments:l};if(r===`turn`&&l.length===1&&o(l[0]))return{tx:null,refused:!0,reason:`ez`,attachments:l};let u;if(l.length===1){let t=e.bonds[l[0]],[n,r]=[e.atoms[t.a].at,e.atoms[t.b].at];u=e=>Nu(e,n,r)}else{let t=Eu(s.map(t=>e.atoms[t].at));u=n===`horizontal`?e=>({x:2*t.x-e.x,y:e.y}):e=>({x:e.x,y:2*t.y-e.y})}let d=[];for(let t of s){let n=e.atoms[t].at,r=u(n);Math.hypot(r.x-n.x,r.y-n.y)<Mu||d.push({op:`setAtomPosition`,atom:t,x:r.x,y:r.y})}for(let t of e.bonds){let e=i(t.index);if(!e||!c.has(e.from)||e.from!==t.a&&e.from!==t.b)continue;let n=r===`mirror`?e.stereo:Pu(e.stereo),a=t.a!==e.from;(a||t.stereo!==n)&&(a&&d.push({op:`reverseBond`,bond:t.index}),d.push({op:`setBondStereo`,bond:t.index,stereo:n}))}return r===`mirror`&&d.length&&d.push({op:`invertChirality`,atoms:s}),{tx:S(`${a?`drawing`:`selection`}.${r===`mirror`?`flip`:`rotate180`}`,d),refused:!1,reason:null,attachments:l}}var Iu=Object.freeze([`E`,`Z`,`unknown`]);function Lu(e,t){let n=e.bonds[t];if(!n)return null;let r=(n,r)=>{let i=new Set([n]),a=[n];for(;a.length;){let n=a.pop();for(let r of e.bondsOf(n)){if(r===t)continue;let o=e.other(r,n);i.has(o)||(i.add(o),a.push(o))}}return i.has(r)?new Set:i},i=r(n.a,n.b);if(!i.size)return null;let a=r(n.b,n.a);return{a:[...i].sort((e,t)=>e-t),b:[...a].sort((e,t)=>e-t)}}var Ru=e=>e===`E`||e===`Z`?e:null,zu=(e,t)=>e.bonds[t.bond]?.stereo===`crossed`||t.status===`unknown`,Bu=(e,t)=>Ru(zu(e,t)?t.drawn:t.cip);function Vu(e,t,n){let r=ei(t,n),i=e.bonds[n];if(!r||!i||i.order!==`double`)return null;let a=zu(e,r)?`unknown`:Ru(r.cip),o=Lu(e,n)===null;if(!o)return{offered:Iu,checked:a,ring:o};let s=Bu(e,r);return{offered:s?[s,`unknown`]:[`unknown`],checked:a,ring:o}}function Hu(e,t,n,r,i){let a=Vu(e,t,n),o=ei(t,n);if(!a||!o||!a.offered.includes(r)||a.checked===r)return null;let s=zu(e,o),c={kind:`ez`,bond:n,choice:r};if(r===`unknown`)return s?null:gt(S(`bond.ez`,[{op:`setBondStereo`,bond:n,stereo:`crossed`}]),c);let l=Bu(e,o);if(!l)return null;let u=s?[{op:`setBondStereo`,bond:n,stereo:`none`}]:[];if(l!==r){let t=Lu(e,n);if(!t)return null;let r=Fu(e,t.a.length<=t.b.length?t.a:t.b,`horizontal`,`turn`,i);if(!r.tx)return null;u.push(...r.tx.ops)}return gt(S(`bond.ez`,u),c)}var Uu=Object.freeze({abs:`stereoMeaningAbs`,or:`stereoMeaningOr`,and:`stereoMeaningAnd`,known:`stereoMeaningKnown`,relativeOr:`stereoMeaningRelativeOr`,relativeAnd:`stereoMeaningRelativeAnd`,unspecified:`stereoMeaningUnspecified`,ezKnown:`stereoMeaningEzKnown`,ezUnknown:`stereoMeaningEzUnknown`,axisKnown:`stereoMeaningAxisKnown`,axisUnspecified:`stereoMeaningAxisUnspecified`}),Wu=Object.freeze(Object.keys(Uu)),Gu=(e,t)=>t[Uu[e]].trim(),Ku=e=>e?.type===`or`?`or`:e?.type===`and`?`and`:`abs`,qu=(e,t)=>e?.type===`or`?`relativeOr`:e?.type===`and`?`relativeAnd`:t,Ju=(e,t)=>`${e.elements[t]||`?`}${t+1}`,Yu=(e,t)=>{let n=e.bonds[t];return n?n.a<=n.b?[n.a,n.b]:[n.b,n.a]:[-1,-1]},Xu=(e,t)=>e.type===`abs`?`abs`:e.type===`and`?t.andName(e.id):t.orName(e.id),Zu=(e,t)=>e===null?t===null:t!==null&&e.type===t.type&&(e.type===`abs`||e.id===t.id);function Qu(e,t,n){if(t===n.setRorSFirst)return{group:null,face:`—`,name:n.format(`stereoGroupChip`,{group:n.setRorSFirst.toLowerCase()}),disabled:t,meaning:null};let r=e?Xu(e,n):n.stereoGroupNone;return{group:e,face:r,name:e?n.format(`stereoGroupChip`,{group:r}):n.format(`stereoGroupChip`,{group:`${n.stereoGroupNone}, ${n.readAsAbsolute}`}),disabled:t,meaning:t===null?Ku(e):null}}function $u(e,t,n,r){let i=t===`r`||t===`s`,a=e===`specified`?t?.toUpperCase():void 0,o=i?`r`:r.stereoR,s=i?`s`:r.stereoS,c=e===`unknown`,l=qu(n,`known`);return[{id:`r`,face:o,name:o,pressed:a===`R`,unset:!1,meaning:l},{id:`s`,face:s,name:s,pressed:a===`S`,unset:!1,meaning:l},{id:`unspecified`,face:`?`,name:c?r.ezUnknown:r.stereoUnspecifiedItem,pressed:e!==`specified`,unset:!0,meaning:`unspecified`}]}var ed={E:`e`,Z:`z`,unknown:`unknown`};function td(e,t){return(e?.offered??[`E`,`Z`,`unknown`]).map(n=>({id:ed[n],face:n===`unknown`?`?`:n===`E`?t.ezE:t.ezZ,name:n===`unknown`?t.ezUnknown:n===`E`?t.ezE:t.ezZ,pressed:e?.checked===n,unset:n===`unknown`,meaning:n===`unknown`?`ezUnknown`:`ezKnown`}))}function nd(e,t,n){let r=e?.toUpperCase(),i=qu(t,`axisKnown`);return[{id:`p`,face:n.atropP,name:n.atropP,pressed:r===`P`,unset:!1,meaning:i},{id:`m`,face:n.atropM,name:n.atropM,pressed:r===`M`,unset:!1,meaning:i}]}function rd(e,t,n){let r=t.find(e=>e.pressed),i=r?r.meaning:e===`double`?`ezUnknown`:e===`axis`?`axisUnspecified`:`unspecified`;return r&&!r.unset&&(n?.type===`or`||n?.type===`and`)?[i,Ku(n)]:[i]}function id(e,t){let n=e.find(e=>e.pressed);return n?n.unset?n.name.toLowerCase():n.name:t.stereoUnspecified}function ad(e,t,n){let r=[];for(let i of[...e.atoms].sort((e,t)=>e.atom-t.atom)){if(i.atom>=t.elements.length)continue;let e=i.status===`specified`?i.group:null,a=$u(i.status,i.cip,e,n),o=i.status===`specified`?Qu(i.group,null,n):Qu(null,n.setRorSFirst,n),s=Ju(t,i.atom),c=id(a,n),l=i.status===`specified`?`, ${n.format(`stereoGroupSaid`,{group:o.face})}`:``;r.push({key:`atom.${i.atom}`,kind:`centre`,atom:i.atom,bond:-1,atoms:[i.atom],name:s,kindWord:n.stereoCentreWord,status:i.status,cip:i.cip,segments:a,chip:o,label:`${n.format(`atomNumber`,{number:i.atom+1})}, ${n.stereoCentreWord}, ${c}${l}`,meanings:rd(`centre`,a,e),problem:i.status===`unspecified`})}for(let i of[...e.bonds].sort((e,t)=>e.bond-t.bond)){if(!t.bonds[i.bond])continue;let[e,a]=Yu(t,i.bond),o=td(t.ez(i.bond),n);r.push({key:`bond.${i.bond}`,kind:`double`,atom:-1,bond:i.bond,atoms:[e,a],name:`${Ju(t,e)}=${Ju(t,a)}`,kindWord:n.stereoDoubleWord,status:i.status,cip:i.cip,segments:o,chip:null,label:`${n.format(`bondNumber`,{number:i.bond+1})}, ${n.stereoDoubleWord}, ${id(o,n)}`,meanings:rd(`double`,o,null),problem:i.status===`unspecified`})}for(let i of[...e.axes].sort((e,t)=>e.bond-t.bond)){if(!t.bonds[i.bond])continue;let[e,a]=Yu(t,i.bond),o=t.axisGroups&&i.status===`specified`?i.group:null,s=nd(i.cip,o,n),c=t.axisGroups?Qu(i.group,null,n):Qu(i.group,n.axesNotGrouped,n);r.push({key:`bond.${i.bond}`,kind:`axis`,atom:-1,bond:i.bond,atoms:[e,a],name:`${Ju(t,e)}–${Ju(t,a)}`,kindWord:n.stereoAxisWord,status:i.status,cip:i.cip,segments:s,chip:c,label:`${n.format(`bondNumber`,{number:i.bond+1})}, ${n.stereoAxisWord}, ${id(s,n)}${t.axisGroups?`, ${n.format(`stereoGroupSaid`,{group:c.face})}`:``}`,meanings:rd(`axis`,s,o),problem:!1})}return r}function od(e,t){let n=new Set(t.atoms),r=new Set(t.bonds);return new Set(e.filter(e=>e.kind===`centre`?n.has(e.atom):r.has(e.bond)).map(e=>e.key))}function sd(e,t){let n=[],r=[];for(let i of e)t.has(i.key)&&(i.kind===`centre`?n.push(i.atom):r.push(i.bond));return{atoms:n,bonds:r}}function cd(e,t,n){if(!n)return new Set([t]);let r=new Set(e);return r.has(t)?r.delete(t):r.add(t),r}var ld=e=>e.chip!==null&&e.chip.disabled===null;function ud(e,t,n){let r=e.find(e=>e.key===n);return!r||!ld(r)?[]:t.has(n)&&t.size>1?e.filter(e=>t.has(e.key)&&ld(e)):[r]}var dd=(e,t)=>Math.max(0,...e.groups.filter(e=>e.type===t).map(e=>e.id))+1;function fd(e,t,n,r,i){let a=t.find(e=>e.key===r),o=ud(t,n,r),s=o.length>1?i.format(`stereoGroupOfSelected`,{count:o.length}):i.format(`stereoGroupOf`,{name:a?.name??``}),c=o.length&&o.every(e=>Zu(e.chip?.group??null,o[0].chip?.group??null))?o[0].chip?.group??null:void 0,l=e=>c!==void 0&&Zu(c,e),u=e=>{let n=t.filter(t=>t.chip&&Zu(t.chip.group,e));return i.count(n.every(e=>e.kind===`centre`)?`groupSizeCentres`:`groupSizeElements`,n.length)},d=t=>[...new Set(e.groups.filter(e=>e.type===t).map(e=>e.id))].sort((e,t)=>e-t),f=[{id:`abs`,label:Gu(`abs`,i)?Xu({type:`abs`,id:0},i):i.stereoAbs,hint:null,checked:l({type:`abs`,id:0}),choice:{type:`abs`},meaning:`abs`}];for(let e of d(`and`))f.push({id:`and.${e}`,label:i.andName(e),hint:u({type:`and`,id:e}),checked:l({type:`and`,id:e}),choice:{type:`and`,id:e},meaning:`and`});for(let e of d(`or`))f.push({id:`or.${e}`,label:i.orName(e),hint:u({type:`or`,id:e}),checked:l({type:`or`,id:e}),choice:{type:`or`,id:e},meaning:`or`});f.push({id:`and.new`,label:i.format(`newAndGroup`,{name:i.andName(dd(e,`and`))}),hint:null,checked:!1,choice:{type:`and`,id:`new`},meaning:`and`}),f.push({id:`or.new`,label:i.format(`newOrGroup`,{name:i.orName(dd(e,`or`))}),hint:null,checked:!1,choice:{type:`or`,id:`new`},meaning:`or`}),f.push({id:`none`,label:i.stereoNone,hint:i.readAsAbsolute,checked:l(null),choice:{type:`none`},meaning:null});let p=a?.chip?.group??null,m=t.filter(e=>ld(e)&&Zu(e.chip?.group??null,p)).map(e=>e.key);return f.push({id:`select-group`,label:p?i.format(`selectAllOfGroup`,{group:Xu(p,i)}):i.selectAllUngrouped,hint:null,checked:null,choice:null,meaning:null}),{heading:s,items:f,selectKeys:m}}function pd(e,t){let n=[...e.atoms.filter(e=>e.status===`specified`).map(e=>e.group),...e.axes.filter(e=>e.status===`specified`).map(e=>e.group)];if(!n.length)return!1;let r=n[0];return!!(r&&r.type===t&&n.every(e=>e&&e.type===t&&e.id===r.id))}function md(e,t,n){if(t===null||Number.isNaN(t)||t<1)return null;if(t===1/0)return n.standsForMany;if(t===1)return n.count(`standsFor`,1);let r=e.atoms.some(e=>e.status!==`specified`)||e.bonds.some(e=>e.status===`unspecified`);return t===2&&!r&&pd(e,`and`)?n.standsForRacemate:t===2&&!r&&pd(e,`or`)?n.standsForOneOf:n.count(`standsFor`,t)}function hd(e,t){let n=/^(.*\S)\s+(\S+)$/s.exec(e.trim()),r=/^(\S+)\s+(.*\S)$/s.exec(t.trim());return{head:n?n[1]:``,tail:n?n[2]:e.trim(),first:r?r[1]:t.trim(),rest:r?r[2]:``}}function gd(e,t,n,r){let i=t.filter(e=>e.problem).length,a=yt(e),o=a?r.format(`flagSaid`,{flag:a===`absolute`?r.flagAbsolute:a===`racemic`?r.flagRacemic:r.flagRelative}):null;return{elements:r.count(`stereoElements`,t.length),unspecified:i?r.count(`unspecifiedCount`,i):null,count:md(e,n,r),flag:o}}function _d(e,t){let n=t.toLowerCase(),r=n===`u`?e.kind===`centre`?`unspecified`:e.kind===`double`?`unknown`:null:n===`r`||n===`s`||n===`e`||n===`z`||n===`p`||n===`m`?n:null;return r&&e.segments.some(e=>e.id===r)?r:null}function vd(e,t,n){return t?n===`Home`?0:n===`End`?t-1:n===`ArrowDown`?Math.min(t-1,e+1):n===`ArrowUp`?Math.max(0,e-1):e:-1}var yd=Object.freeze({regular:288,compact:216}),bd=Object.freeze({share:.42,min:172,max:300});function xd(e,t,n){let r=t?yd.compact:yd.regular,i=Math.round(Math.min(bd.max,Math.max(bd.min,e.height*bd.share)));return e.width-r>=300+(n===`bottom`||n===`cover`?24:0)?{place:`right`,width:r,height:e.height}:e.width>=yd.compact&&e.height-i>=140+(n===`cover`?24:0)?{place:`bottom`,width:e.width,height:i}:{place:`cover`,width:e.width,height:e.height}}var Sd=e=>{let t=typeof e==`string`?e.trim().toLowerCase():``;return t===`open`||t===`hidden`?t:`closed`},Cd=Object.freeze({elements:Object.freeze([]),count:null,flag:null});function wd(e,t,n){let r=t.map(t=>Object.freeze({kind:t.kind===`centre`?`centre`:t.kind===`double`?`double-bond`:`axis`,atoms:Object.freeze([...t.atoms]),bond:t.kind===`centre`?null:t.bond,label:t.cip,status:t.status,group:t.kind===`centre`?Td(e,t.atom):t.kind===`axis`?e.axes.find(e=>e.bond===t.bond)?.group??null:null}));return Object.freeze({elements:Object.freeze(r),count:n!==null&&n>0?n:null,flag:yt(e)})}var Td=(e,t)=>{let n=e.atoms.find(e=>e.atom===t)?.group??null;return n?Object.freeze({type:n.type,id:n.type===`abs`?0:n.id}):null},Ed={kind:`bond`,type:`single`,stereo:`none`},Dd=()=>`auto`,Od=e=>e!==`single`;function kd(e,t){let n=t.length;return t.map((r,i)=>{let a=t[(i+1)%n];return r>=0&&a>=0&&e.bondBetween(r,a)>=0})}function Ad(e,t=[]){let n=e.length;return e.map((r,i)=>+!e[i]+ +!e[(i-1+n)%n]+(t[i]??0))}function jd(e,t,n,r=[]){return t.filter((t,i)=>t>=0&&!r.includes(i)&&!Fn(e,t,n[i])).length===0}function Md(e,t,n,r=[]){let i=Ad(n,r);return t.map((t,n)=>{if(t<0)return{unsaturated:!1,room:1/0};let r=e.info(t)?.hCount;return{unsaturated:e.bondsOf(t).some(t=>Od(e.bonds[t].order)),room:r===void 0?1/0:r-i[n]}})}function Nd(e,t,n,r){let i=e.length,a=[],o=[],s=0;e.forEach((e,n)=>{t[n]>=0?o.push(t[n]):(a.push({op:`addAtom`,symbol:`C`,x:e.x,y:e.y}),o.push({added:s++}))});for(let e=0;e<i;e++)r[e]||a.push({op:`addBond`,a:o[e],b:o[(e+1)%i],type:n[e],stereo:`none`});return{ops:a,refs:o}}var Pd=e=>kl*e.docLen,Fd=e=>e.atoms.map(e=>e.hidden?lr:e.at);function Id(e,t,n,r,i){let a=Ol[e],o=Pl(a,n,r,t.docLen),s=Hl(o,Fd(t),Vl(o,Fd(t),Pd(t)),1*t.docLen/(t.bondPx||1)),c=Vl(s,Fd(t),Pd(t)),l=kd(t,c),u=l.indexOf(!0);if(u>=0){let n=t.bondBetween(c[u],c[(u+1)%a.size]),{a:r,b:i}=t.bonds[n];return Vd(e,t,n,Bl(t.atoms[r].at,t.atoms[i].at,zd(s))??Bd(t,n))}if(!jd(t,c,Ad(l)))return null;let d=i?[Nl(a.orders),a.orders]:[a.orders,Nl(a.orders)],f=Md(t,c,l),p=Yl(d.map((e,t)=>({orders:e,rank:t})),f,l),m=a.orders.filter(e=>e===`double`).length;return p.filter(e=>e===`double`).length<m?(c.forEach((e,n)=>{e>=0&&(f[n].unsaturated||f[n].room<1)&&Ie({atom:e,symbol:t.atoms[e]?.element||t.info(e)?.symbol||``})}),null):S(`ring.place`,Nd(s,c,p,l).ops)}function Ld(e,t){let{a:n,b:r}=e.bonds[t],i=new Map([[n,n]]),a=[n];for(;a.length;){let n=a.shift();if(n===r)break;for(let r of e.bondsOf(n)){if(r===t)continue;let o=e.other(r,n);i.has(o)||(i.set(o,n),a.push(o))}}if(!i.has(r))return null;let o=[r];for(;o[o.length-1]!==n;)o.push(i.get(o[o.length-1]));return o}function Rd(e,t){let n=new Set([t]),r=[t];for(;r.length;){let t=r.shift();for(let i of e.bondsOf(t)){let a=e.other(i,t);n.has(a)||(n.add(a),r.push(a))}}return[...n]}var zd=e=>({x:e.reduce((e,t)=>e+t.x,0)/e.length,y:e.reduce((e,t)=>e+t.y,0)/e.length});function Bd(e,t){let{a:n,b:r}=e.bonds[t],i=Ld(e,t)??Rd(e,n);return zl(e.atoms[n].at,e.atoms[r].at,zd(i.map(t=>e.atoms[t].at)))}function Vd(e,t,n,r){let i=Ol[e],{a,b:o,order:s}=t.bonds[n],c=Fl(t.atoms[a].at,t.atoms[o].at,i.size,r??Bd(t,n)),l=Vl(c,Fd(t),Pd(t),[a,o]),u=kd(t,l);if(!jd(t,l,Ad(u),[0,1]))return null;let d=Yl(Kl(i,s),Md(t,l,u),u);return S(`ring.fuse`,Nd(c,l,d,u).ops)}function Hd(e,t){let n=e.info(t)?.hCount;return n===void 0||n>=2}function Ud(e,t,n){let r=t.atoms[n].at,i=t.bondsOf(n).map(e=>t.atoms[t.other(e,n)].at),a=Il(r,e.size,Rl(r,i),t.docLen),o=Vl(a,Fd(t),Pd(t),[n]),s=kd(t,o);if(!jd(t,o,Ad(s),[0]))return null;let{ops:c,refs:l}=Nd(a,o,Yl(ql(e),Md(t,o,s),s),s),u=S(`ring.spiro`,c);return u?{tx:u,refs:l}:null}function Wd(e,t,n){let r=t.atoms[n].at,i=Zr(r,Ut(t,n),`single`),a=Ll(r,e.size,i,t.docLen),o=Vl(a,Fd(t),Pd(t),[],[n]),s=kd(t,o),c=o[0]>=0&&t.bondBetween(n,o[0])>=0,l=[+!c];if(!jd(t,o,Ad(s,l)))return null;let{ops:u,refs:d}=Nd(a,o,Yl([{orders:e.orders,rank:0}],Md(t,o,s,l),s),s),f=S(`ring.sprout`,c?u:[...u,{op:`addBond`,a:n,b:d[0],type:`single`,stereo:`none`}]);return f?{tx:f,refs:d}:null}function Gd(e,t,n,r){let i=Ol[e];return(!r&&Hd(t,n)?Ud(i,t,n):Wd(i,t,n))?.tx??null}var Kd=e=>({onCanvas:(t,n,r,i)=>Id(e,t,n,r,i),onAtom:(t,n,r)=>Gd(e,t,n,r),onBond:(t,n,r)=>Vd(e,t,n,r)});function qd(e,t,n,r,i,a){let o=`atom`in t?t.atom:-1,s=o>=0?e.atoms[o].at:t.point,{direction:c,points:l}=oa(s,n,e.docLen,r,i),u=a?l.slice(0,Yd(e,l,a)):l;if(!u.length)return null;let d=u.length-1,f=Vl([u[d]],Fd(e),Pd(e),[],o>=0?[o]:[])[0],p=u.map((e,t)=>t===d?f:-1),m=[],h=0,g=o>=0?o:{added:h++};o<0&&m.push({op:`addAtom`,symbol:`C`,x:s.x,y:s.y});let _=u.map((e,t)=>p[t]>=0?p[t]:(m.push({op:`addAtom`,symbol:`C`,x:e.x,y:e.y}),{added:h++})),v=new Set(p.filter(e=>e>=0)),y=new Map,b=g;for(let t of _){if(!(typeof b==`number`&&typeof t==`number`&&(b===t||e.bondBetween(b,t)>=0))){m.push({op:`addBond`,a:b,b:t,type:`single`,stereo:`none`});for(let e of[b,t])typeof e==`number`&&v.has(e)&&y.set(e,(y.get(e)??0)+1)}b=t}return[...y].filter(([t,n])=>!Fn(e,t,n)).length?null:{tx:S(`chain.draw`,m),count:u.length+(o>=0?0:1),origin:s,end:f>=0?e.atoms[f].at:u[d],direction:c}}function Jd(e,t,n){if(!Number.isInteger(n)||n<1||n>50)return null;let r=`atom`in t?t.atom:-1;if(`atom`in t&&(r<0||r>=e.atoms.length))return null;let i=r>=0?e.atoms[r].at:t.point,a=r>=0?n:n-1;if(a===0)return{tx:S(`chain.draw`,[{op:`addAtom`,symbol:`C`,x:i.x,y:i.y}]),count:1,origin:i,end:i,direction:0};let o=r>=0?k(Zr(i,Ut(e,r),`single`)-30):0;return qd(e,t,O(i,o,(a-.5)*e.docLen),!0,!1)}function Yd(e,t,n){let r=t.findIndex(t=>{let r=e.toCanvas(t);return r.x<0||r.y<0||r.x>n.width||r.y>n.height});return r<0?t.length:r}Object.freeze({a:`benzene`,4:`cyclobutane`,5:`cyclopentane`,6:`cyclohexane`,7:`cycloheptane`,8:`cyclooctane`});var Xd=1.75;function Zd(e,t){let n=e.atoms[t].at,r=Zr(n,Ut(e,t),`single`),i=e.bondsOf(t).map(r=>A(n,e.atoms[e.other(r,t)].at));if(i.length===0)return[r,k(r+120)];if(i.length===1){let e=i[0],[t,n]=[k(e+120),k(e-120)];return Math.abs(k(t-r))<Math.abs(k(n-r))?[t,n]:[n,t]}let a=Zt(i),o=[...i].sort((e,t)=>e-t),s=0;return o.forEach((e,t)=>{let n=t+1<o.length?o[t+1]:o[0]+360;s=Math.max(s,n-e)}),[k(a-s/6),k(a+s/6)]}var U=(e,t,n,r,i=`single`,a=`none`)=>[{op:`addAtom`,symbol:e,x:t.x,y:t.y},{op:`addBond`,a:n,b:{added:r},type:i,stereo:a}];function Qd(e,t,n,r,i,a,o){let s=Zr(e,[{at:t,order:`single`,beyond:[...n]}],`single`),c=A(e,t),[l,u]=[k(c+120),k(c-120)],d=Math.abs(k(l-s))<Math.abs(k(u-s))?u:l;return[...U(a[0].symbol,O(e,s,o),r,i,a[0].type,a[0].stereo),...U(a[1].symbol,O(e,d,o),r,i+1,a[1].type,a[1].stereo)]}function $d(e,t){return(e.atoms[t]?.element||e.info(t)?.symbol)===`C`&&e.bondsOf(t).some(t=>e.bonds[t]?.type===`aromatic`||e.perceivedAromatic(t))}var ef=e=>e===`double`?2:e===`triple`?3:1;function tf(e,t,n){let r=Vl(t.filter(e=>e.op===`addAtom`).map(e=>({x:e.x,y:e.y})),Fd(e),Pd(e),[],n);if(r.every(e=>e<0))return{ops:[...t],ref:e=>e};let i=0,a=r.map(e=>e<0?i++:-1),o=e=>typeof e==`number`?e:r[e.added]>=0?r[e.added]:{added:a[e.added]},s=new Set(r.filter(e=>e>=0)),c=new Map,l=[],u=0;for(let n of t)if(n.op===`addAtom`)r[u++]<0&&l.push(n);else if(n.op===`addBond`){let[t,r]=[o(n.a),o(n.b)];if(typeof t==`number`&&typeof r==`number`&&(t===r||e.bondBetween(t,r)>=0))continue;l.push({...n,a:t,b:r});for(let e of[t,r])typeof e==`number`&&s.has(e)&&c.set(e,(c.get(e)??0)+ef(n.type))}else l.push(n);return[...c].filter(([t,n])=>!Fn(e,t,n)).length?null:{ops:l,ref:o}}function nf(e,t,n,r=Dd){let i=t.atoms[n].at,a=t.docLen,o=t.info(n)?.hCount??0,s=t.bondsOf(n).length,c=Zr(i,Ut(t,n),`single`),l=t.bondsOf(n).map(e=>t.atoms[t.other(e,n)].at),u=(e,r,i)=>{let a=tf(t,r,[n]),o=a?S(e,a.ops):null;return o&&a?{tx:o,hotspot:i===null?null:a.ref(i)}:null},d=(e,r)=>{let i=Ol[e],a=r?Ud(i,t,n):Wd(i,t,n);return a?{tx:{...a.tx,label:`atom.sprout`},hotspot:a.refs[+!!r]}:null};switch(e){case`1`:{let e=lt(Ed,t,n,r);if(!e)return null;let i=[...e.ops].reverse().find(e=>e.op===`addBond`);return{tx:{...e,label:`atom.sprout`},hotspot:i&&i.op===`addBond`?i.b:null}}case`2`:{if(o>=3){let[e,r]=Zd(t,n);return u(`atom.sprout`,[...U(`O`,O(i,e,a),n,0,`double`),...U(`C`,O(i,r,a),n,1)],{added:0})}if(o===2)return u(`atom.sprout`,U(`O`,O(i,Zr(i,Ut(t,n),`double`),a),n,0,`double`),{added:0});let e=O(i,c,a);return u(`atom.sprout`,[...U(`C`,e,n,0),...Qd(e,i,l,{added:0},1,[{symbol:`O`,type:`double`},{symbol:`C`,type:`single`}],a)],{added:0})}case`3`:case`a`:return d(`benzene`,o>=3);case`4`:case`5`:{let r=e===`4`?`wedge`:`hash`;if(o>=3){let[e,o]=Zd(t,n);return u(`atom.sprout`,[...U(`C`,O(i,e,a),n,0,`single`,r),...U(`C`,O(i,o,a),n,1)],{added:0})}if($d(t,n)){let e=O(i,c,a);return u(`atom.sprout`,[...U(`C`,e,n,0),...Qd(e,i,l,{added:0},1,[{symbol:`C`,type:`single`,stereo:r},{symbol:`C`,type:`single`}],a)],{added:0})}return u(`atom.sprout`,U(`C`,O(i,c,a),n,0,`single`,r),{added:0})}case`6`:return d(`cyclohexane`,Hd(t,n));case`7`:return d(`cyclopentane`,Hd(t,n));case`8`:return o<2?null:u(`atom.sprout`,U(`C`,O(i,Zr(i,Ut(t,n),`double`),a),n,0,`double`),{added:0});case`9`:{if(o>=2){let[e,r]=Zd(t,n);return u(`atom.sprout`,[...U(`C`,O(i,e,a),n,0),...U(`C`,O(i,r,a),n,1)],{added:0})}let e=O(i,c,a);return u(`atom.sprout`,[...U(`C`,e,n,0),...Qd(e,i,l,{added:0},1,[{symbol:`C`,type:`single`},{symbol:`C`,type:`single`}],a)],{added:0})}case`0`:if(s===1){let e=k(2*A(i,l[0])-c);return u(`atom.sprout`,U(`C`,O(i,e,a),n,0),{added:0})}return s<2?null:u(`atom.sprout`,U(`C`,O(i,c,Xd*a),n,0),{added:0});case`K`:{let e=o>=3,t=e?i:O(i,c,a),r=e?n:{added:0},s=+!e,d=e&&l.length?A(i,l[0]):k(c+180),f=[180,90,-90].flatMap((e,n)=>U(`C`,O(t,k(d+e),a),r,s+n));return u(`atom.sprout`,e?f:[...U(`C`,t,n,0),...f],{added:s})}case`v`:return d(`cyclopropane`,Hd(t,n));case`u`:return d(`cyclobutane`,Hd(t,n));case`Z`:case`M`:case`k`:{let t=o>=3,r=e===`Z`?`N`:e===`M`?`Mg`:`S`,s=t?n:{added:0},d=t?i:O(i,c,a),f=t?[{op:`setAtomElement`,atom:n,symbol:r}]:U(r,d,n,0),p=+!t,m=t&&l.length?k(A(i,l[0])+180):c,h=t?l.length?l[0]:O(i,k(c+180),a):i;if(e===`k`)return f.push(...Qd(d,h,[],s,p,[{symbol:`O`,type:`double`},{symbol:`O`,type:`double`}],a)),u(`atom.sprout`,f,s);if(e===`M`)return f.push(...U(`Br`,O(d,m,a),s,p)),u(`atom.sprout`,f,{added:p});let g=k(A(h,d)+60),_=O(d,g,a),v=O(_,g,a);return f.push(...U(`N`,_,s,p,`double`),...U(`N`,v,{added:p},p+1,`double`)),f.push({op:`setAtomCharge`,atom:{added:p},charge:1},{op:`setAtomCharge`,atom:{added:p+1},charge:-1}),u(`atom.sprout`,f,{added:p+1})}case`z`:{if(o>=3){let e=l.length?k(A(i,l[0])+180):c,t=O(i,e,a);return u(`atom.sprout`,[{op:`addAtom`,symbol:`C`,x:t.x,y:t.y},{op:`addBond`,a:n,b:{added:0},type:`triple`,stereo:`none`}],{added:0})}let e=O(i,c,a),t=O(e,c,a);return u(`atom.sprout`,[...U(`C`,e,n,0),{op:`addAtom`,symbol:`C`,x:t.x,y:t.y},{op:`addBond`,a:{added:0},b:{added:1},type:`triple`,stereo:`none`}],{added:1})}}}var rf=`sgroup.add`,af=`sgroup.edit`,of=`sgroup.remove`,sf=[`SUP`,`DAT`,`MUL`,`SRU`],cf=[`atom`,`bond`,`fragment`,`group`,`multifragment`],lf=[`attached`,`absolute`,`relative`],uf=[`HT`,`HH`,`EU`],df=e=>sf.includes(e);function ff(e,t=`DAT`){return{type:t,name:``,context:e,fieldName:``,value:``,placement:`attached`,count:`1`,label:`n`,connect:`HT`}}function pf(e){let t=typeof e==`string`?e:``;return t.charAt(24)===`A`?`attached`:t.charAt(25)===`R`?`relative`:t.length>=26?`absolute`:`attached`}var W=e=>typeof e==`string`?e:e==null?``:String(e);function mf(e,t){let n=String(e.type??``);if(!df(n))return null;let r=ff(t,n),i=e.props??{};switch(n){case`SUP`:return{...r,name:W(e.label)};case`DAT`:return{...r,fieldName:W(i.FIELDNAME),value:(e.fields??[]).join(`
`),placement:pf(i.FIELDDISP)};case`MUL`:return{...r,count:W(e.label)||`1`};case`SRU`:{let t=W(i.CONNECT).toUpperCase();return{...r,label:W(e.label)||`n`,connect:uf.includes(t)?t:`HT`}}}}function hf(e,t){let n=Array.from({length:e},(e,t)=>t),r=e=>{for(;n[e]!==e;)n[e]=n[n[e]],e=n[e];return e};for(let i of t){if(i.a<0||i.b<0||i.a>=e||i.b>=e)continue;let[t,a]=[r(i.a),r(i.b)];t!==a&&(n[Math.max(t,a)]=Math.min(t,a))}let i=new Map;return n.map((e,t)=>{let n=r(t);return i.has(n)||i.set(n,i.size),i.get(n)})}function gf(e,t){let n=hf(e.atoms,e.bonds),r=new Set(t.filter(t=>t>=0&&t<e.atoms).map(e=>n[e]));return n.flatMap((e,t)=>r.has(e)?[t]:[])}function _f(e,t,n){let r=[...new Set(t)].sort((e,t)=>e-t);if(r.length===1&&!n.length)return`atom`;if(n.length===1&&r.length===2){let t=e.bonds[n[0]];if(t&&[t.a,t.b].sort((e,t)=>e-t).join(` `)===r.join(` `))return`bond`}if(!r.length)return`group`;let i=hf(e.atoms,e.bonds),a=new Set(r.map(e=>i[e]));return i.every((e,t)=>!a.has(e)||r.includes(t))?a.size===1?`fragment`:`multifragment`:`group`}function vf(e,t,n){let r=[...new Set(t.atoms)].sort((e,t)=>e-t);if(n.type!==`DAT`)return{atoms:r,bonds:[]};switch(n.context){case`bond`:{let n=t.bonds.length===1?t.bonds[0]:r.length===2?e.bonds.findIndex(e=>e.a===r[0]&&e.b===r[1]||e.a===r[1]&&e.b===r[0]):-1;return n>=0?{atoms:[],bonds:[n]}:{atoms:r,bonds:[]}}case`fragment`:case`multifragment`:return{atoms:gf(e,r),bonds:[]};default:return{atoms:r,bonds:[]}}}function yf(e,t){switch(e.type){case`SUP`:{let n=e.name.trim();return n?n.length>15?t.sgroupNameTooLong(15):null:t.sgroupNeedsName}case`DAT`:return!e.fieldName.trim()||!e.value.trim()?t.sgroupNeedsField:e.fieldName.trim().length>30?t.sgroupFieldNameTooLong(30):e.value.split(`
`).some(e=>e.length>200)?t.sgroupValueTooLong(200):null;case`MUL`:return bf(e.count)===null?t.sgroupCountRange(200):null;case`SRU`:{let n=e.label.trim();return!n||n.includes(`"`)||n.includes(`
`)?t.sgroupNeedsLabel:n.length>69?t.sgroupLabelTooLong(69):null}}}function bf(e){let t=e.trim();if(!/^\d{1,3}$/.test(t))return null;let n=Number(t);return n>=1&&n<=200?n:null}function xf(e,t){if(!e.length)return{x:0,y:0};let n=e.map(e=>e.x),r=e.map(e=>e.y);return{x:Math.max(...n)+.5*t,y:(Math.min(...r)+Math.max(...r))/2}}var Sf=e=>Math.round(e*1e4)/1e4+0;function Cf(e,t,n,r){if(e.placement===`attached`)return{placement:`attached`};let i=xf(t.map(n),r);if(e.placement===`absolute`)return{placement:`absolute`,x:Sf(i.x),y:Sf(i.y)};let a=n(Math.min(...t));return{placement:`relative`,x:Sf(i.x-a.x),y:Sf(i.y-a.y)}}function wf(e,t,n,r,i){switch(e.type){case`SUP`:return{label:e.name.trim()};case`DAT`:return{fieldName:e.fieldName.trim(),value:e.value,...Cf(e,t,n,r),...i?{context:e.context}:{}};case`MUL`:return{count:bf(e.count)??1};case`SRU`:return{label:e.label.trim(),connect:e.connect}}}function Tf(e,t,n,r){if(yf(n,r))return null;let{atoms:i,bonds:a}=vf(e,t,n),o=[...new Set([...i,...a.flatMap(t=>e.bonds[t]?[e.bonds[t].a,e.bonds[t].b]:[])])];if(!o.length)return null;let s=wf(n,o,t=>e.position(t),e.docLen,!0);return S(rf,[{op:`addSGroup`,type:n.type,atoms:i,bonds:a,props:JSON.stringify(s)}])}function Ef(e,t,n,r){if(yf(n,r)||n.type!==t.type)return null;let i=mf(t,n.context);if(!i)return null;let a=wf(n,[...t.atoms??[]],t=>e.position(t),e.docLen,!1),o={};switch(n.type){case`SUP`:n.name.trim()!==i.name&&(o.label=a.label);break;case`DAT`:n.fieldName.trim()!==i.fieldName&&(o.fieldName=a.fieldName),n.value!==i.value&&(o.value=a.value),n.placement!==i.placement&&Object.assign(o,{placement:a.placement,...a.x===void 0?{}:{x:a.x,y:a.y}});break;case`MUL`:bf(n.count)!==bf(i.count)&&(o.count=a.count);break;case`SRU`:n.label.trim()!==i.label&&(o.label=a.label),n.connect!==i.connect&&(o.connect=a.connect)}return Object.keys(o).length?S(af,[{op:`setSGroupProps`,k:t.index,props:JSON.stringify(o)}]):null}function Df(e){return e>=0?S(of,[{op:`removeSGroup`,k:e}]):null}function Of(e,t,n){let r=e.trim();return/partial s-group overlapping/i.test(r)?n.sgroupPartialOverlap:/of context \w+ holds/i.test(r)?n.sgroupContextMismatch(n.sgroupContextName(t.context)):/multiple group needs no crossing bond or two/i.test(r)?n.sgroupMulBonds:n.sgroupRefused}function kf(e,t){return t.length?e.filter(e=>e.type!==`SUP`&&t.every(t=>(e.atoms??[]).includes(t))).sort((e,t)=>(e.atoms?.length??0)-(t.atoms?.length??0)):[]}function Af(e,t){let n=e.props??{};switch(e.type){case`DAT`:{let r=W(n.FIELDNAME),i=(e.fields??[]).join(` `).trim();return t.sgroupTipData(r,i)}case`SRU`:return t.sgroupTipSru(W(e.label)||`n`,t.sgroupConnectName(uf.includes(W(n.CONNECT).toUpperCase())?W(n.CONNECT).toUpperCase():`HT`));case`MUL`:return t.sgroupTipMul(W(e.label));case`SUP`:return t.sgroupTipSup(W(e.label));default:return t.sgroupTipOther(String(e.type??``))}}var jf=e=>({x:e.reduce((e,t)=>e+t.x,0)/e.length,y:e.reduce((e,t)=>e+t.y,0)/e.length});function Mf(e,t,n,r,i){let a=i*Math.PI/180,[o,s]=[Math.cos(a),Math.sin(a)];return e.atoms.map(e=>{let[i,a]=[e.x-t.x,e.y-t.y];return{x:n.x+r*(i*o-a*s),y:n.y+r*(i*s+a*o)}})}var Nf=(e,t)=>t/Un(e);function Pf(e,t,n){if(!n.length)return null;let r=jf(n.map(t=>e.atoms[t])),i={x:r.x-e.atoms[t].x,y:r.y-e.atoms[t].y};return Math.hypot(i.x,i.y)>1e-9?i:null}function Ff(e,t,n,r){return Mf(e,jf(e.atoms),t,Nf(e,r),n)}function If(e,t,n,r,i,a){let o=Pf(e,t,a),s=o?r-Math.atan2(o.y,o.x)*180/Math.PI:0;return Mf(e,e.atoms[t],n,Nf(e,i),s)}function Lf(e,t,n,r,i){let[a,o]=[e.atoms[t],e.atoms[n]],s=Math.hypot(o.x-a.x,o.y-a.y)||1,c=A(r,i)-Math.atan2(o.y-a.y,o.x-a.x)*180/Math.PI;return Mf(e,a,r,Math.hypot(i.x-r.x,i.y-r.y)/s,c)}function Rf(e,t,n){let r=(n.x-t.x)*(e.y-t.y)-(n.y-t.y)*(e.x-t.x),i=Math.hypot(n.x-t.x,n.y-t.y)||1;return Math.abs(r)/i<1e-9?0:Math.sign(r)}var zf=e=>kl*e.docLen,Bf=e=>e.atoms.map(e=>e.hidden?lr:e.at),Vf=e=>e===`double`?1:e===`triple`?2:0,Hf=[`single`,`double`,`triple`],Uf=8;function Wf(e,t){if(!t.size)return[];let n=t=>{let n=e.bonds[t];return n.stereo===`none`?n.order===`double`?`D`:n.order===`single`?`S`:null:null},r=e.atoms.map(()=>[]);e.bonds.forEach((e,t)=>{r[e.a].push({to:e.b,bond:t}),r[e.b].push({to:e.a,bond:t})});let i=new Set,a=[],o=t=>{let n=[...t].sort((e,t)=>e-t).join(` `);if(i.has(n))return;i.add(n);let r=new Set(t);a.push(e.bonds.map((e,t)=>r.has(t)?{...e,order:e.order===`double`?`single`:`double`}:e))};for(let e of t){let t=[],i=[e],a=s=>{if(!(t.length>=Uf))for(let{to:c,bond:l}of r[s]){let r=n(l);if(!(!r||t.includes(l)||t.length&&n(t[t.length-1])===r)){if(c===e){t.length>=3&&(t.length+1)%2==0&&n(t[0])!==r&&o([...t,l]);continue}i.includes(c)||(t.push(l),i.push(c),a(c),t.pop(),i.pop())}}};a(e)}return a}function Gf(e,t,n,r){if(t.unplaceable)return null;let i=t.atoms.map((e,t)=>n.own.get(t)??-1),a=Vl(n.points,Bf(e),zf(e),i,n.link?[n.link.atom]:[]),o=new Set(a.map((e,t)=>e>=0?t:-1).filter(e=>e>=0)),s=qf(e,t,n,r,a);for(let i of Wf(t,o)){if(s.lost===0)break;let o=qf(e,{...t,bonds:i},n,r,a);o.lost<s.lost&&(o.tx||!s.tx)&&(s=o)}return s.tx}var Kf=e=>e===`double`||e===`aromatic`?1:e===`triple`?2:0;function qf(e,t,n,r,i){let a=new Set(n.own.values()),o=[],s=[],c=0;n.points.forEach((e,n)=>{i[n]>=0?s.push(i[n]):(o.push({op:`addAtom`,symbol:t.atoms[n].symbol,x:e.x,y:e.y}),s.push({added:c++}))});let l=(t,n)=>typeof t==`number`&&typeof n==`number`&&(t===n||e.bondBetween(t,n)>=0),u=t.bonds.filter(e=>!l(s[e.a],s[e.b])),d=n.link&&!l(n.link.atom,s[n.link.to])?n.link:void 0,f=new Map,p=(e,t)=>typeof e==`number`&&f.set(e,(f.get(e)??0)+t);for(let e of u)e.order!==`dative`&&(p(s[e.a],1),p(s[e.b],1));d&&(p(d.atom,1),p(s[d.to],1));let m=t=>e.bondsOf(t).some(t=>e.bonds[t].order!==`single`&&e.bonds[t].order!==`dative`),h=t=>{let n=e.info(t)?.hCount;return(n===void 0?1/0:n)-(f.get(t)??0)},g=u.map(e=>{if(e.order===`single`||e.order===`aromatic`||e.order===`dative`)return e.order;let t=[s[e.a],s[e.b]].filter(e=>typeof e==`number`),n=Vf(e.order);for(;n>0&&!t.every(e=>!m(e)&&h(e)>=n);)n--;for(let e of t)f.set(e,(f.get(e)??0)+n);return Hf[n]}),_=0;u.forEach((e,t)=>_+=Math.max(0,Vf(e.order)-Vf(g[t])));for(let n of t.bonds){let[t,r]=[s[n.a],s[n.b]];if(typeof t!=`number`||typeof r!=`number`||t===r)continue;let i=e.bondBetween(t,r);i>=0&&(_+=Math.max(0,Vf(n.order)-Kf(e.bonds[i].order)))}return[...f].filter(([t,n])=>!a.has(t)&&t!==d?.atom&&!Fn(e,t,n)).length?{tx:null,lost:_}:(u.forEach((e,t)=>o.push({op:`addBond`,a:s[e.a],b:s[e.b],type:g[t],stereo:g[t]===`single`?e.stereo:`none`})),d&&o.push({op:`addBond`,a:d.atom,b:s[d.to],type:`single`,stereo:`none`}),o.push(...Jf(e,t,s)),t.atoms.forEach((e,t)=>{let n=s[t];typeof n!=`number`&&(e.charge&&o.push({op:`setAtomCharge`,atom:n,charge:e.charge}),e.isotope&&o.push({op:`setAtomIsotope`,atom:n,isotope:e.isotope}),e.radical&&o.push({op:`setAtomRadical`,atom:n,radical:e.radical}),e.valence!==null&&o.push({op:`setAtomValence`,atom:n,valence:e.valence}))}),{tx:S(r,o),lost:_})}function Jf(e,t,n){let r=[];for(let i of t.bonds){let[t,a]=[n[i.a],n[i.b]];if(i.stereo!==`wedge`&&i.stereo!==`hash`||i.order!==`single`||typeof t!=`number`||typeof a!=`number`||t===a)continue;let o=e.bondBetween(t,a),s=o>=0?e.bonds[o]:void 0;s&&s.order===`single`&&(s.stereo===`none`?(s.a!==t&&r.push({op:`reverseBond`,bond:o}),r.push({op:`setBondStereo`,bond:o,stereo:i.stereo})):(s.stereo===`wedge`||s.stereo===`hash`)&&s.a===a&&r.push({op:`setBondStereo`,bond:o,stereo:i.stereo,from:t}))}return r}function Yf(e,t,n,r){let i=Ff(e,n,r,t.docLen),a=Hl(i,Bf(t),Vl(i,Bf(t),zf(t)),1*t.docLen/(t.bondPx||1)),o=Vl(a,Bf(t),zf(t)),s=e.bonds.find(e=>o[e.a]>=0&&o[e.b]>=0&&o[e.a]!==o[e.b]&&t.bondBetween(o[e.a],o[e.b])>=0);if(s){let[n,r]=[o[s.a],o[s.b]];return Gf(t,e,{points:Lf(e,s.a,s.b,t.atoms[n].at,t.atoms[r].at),own:new Map([[s.a,n],[s.b,r]])},`template.fuse`)}return Gf(t,e,{points:a,own:new Map},`template.place`)}var Xf=(e,t)=>e.bonds.filter(e=>e.a===t||e.b===t).length,Zf=(e,t)=>Rl(e.atoms[t].at,e.bondsOf(t).map(n=>e.atoms[e.other(n,t)].at));function Qf(e,t,n,r){let i=e.attachAtom,a=e.atoms[i]?.symbol;if(a===void 0)return null;let o=t.info(n)?.hCount,s=o===void 0||o>=Xf(e,i),c=t.atoms[n].at,l=Zf(t,n),u=e.atoms.map((e,t)=>t).filter(e=>e!==i);return a===`H`||!r&&a===`C`&&s?s?Gf(t,e,{points:If(e,i,c,l,t.docLen,u),own:new Map([[i,n]])},`template.through`):(Ie({atom:n,symbol:t.atoms[n]?.element||t.info(n)?.symbol||``}),null):e.atoms[i].hCount<1?(Ie({atom:-1,symbol:a}),null):Gf(t,e,{points:If(e,i,O(c,l,t.docLen),l,t.docLen,e.atoms.map((e,t)=>t)),own:new Map,link:{atom:n,to:i}},`template.sprout`)}function $f(e,t,n,r){let i=e.bonds[e.attachBond];if(!i)return null;let{a,b:o}=t.bonds[n],[s,c]=[t.atoms[a].at,t.atoms[o].at],l=r??Bd(t,n),u=e.atoms.map((e,t)=>t).filter(e=>e!==i.a&&e!==i.b),d=[{points:Lf(e,i.a,i.b,s,c),own:new Map([[i.a,a],[i.b,o]])},{points:Lf(e,i.a,i.b,c,s),own:new Map([[i.a,o],[i.b,a]])}],f=e=>u.length?Rf(jf(u.map(t=>e[t])),s,c):0;return Gf(t,e,d.find(e=>f(e.points)===l)??d.find(e=>f(e.points)===0)??d[0],`template.fuse`)}var ep=e=>({onCanvas:(t,n,r)=>Yf(e,t,n,r),onAtom:(t,n,r)=>Qf(e,t,n,r),onBond:(t,n,r)=>$f(e,t,n,r)}),tp=e=>typeof e==`string`&&e.trim().toLowerCase()===`hidden`?`hidden`:`shown`,np=e=>e===`AtomValenceException`?`valence`:`aromaticity`,rp={valence:0,aromaticity:1,label:2,abbreviation:3,stereo:4,contact:5,linear:6},ip=e=>!(!e||e.symbol!==`C`||e.charge!==0||e.isotope!==0||e.alias!==null||e.shown);function ap(e){let t=0;for(let n of e)if(n===`(`)t+=1;else if(n===`)`&&--t<0)return!1;return t===0}var op=e=>e&&e.atomicNumber===0&&e.symbol===`*`&&e.alias!==null?e.alias:null;function sp(e,t,n=p){return e===`aromaticity`?n.problemAromatic:e===`stereo`?n.problemStereo:e===`abbreviation`?n.problemAbbreviation:e===`contact`?n.problemContact:e===`linear`?n.problemLinear:e===`label`?ap(op(t)??``)?n.problemLabel:n.problemParentheses:ip(t)?n.problemCarbon:n.problemValence(t?.symbol||`?`)}function cp(e){let t=e=>e.atoms.length?Math.min(...e.atoms):1/0;return e.map((e,t)=>({p:e,k:t})).sort((e,n)=>t(e.p)-t(n.p)||rp[e.p.kind]-rp[n.p.kind]||e.k-n.k).map(e=>e.p)}var lp=new Set([`wedge`,`hash`,`wavy`]);function up(e,t,n){let r=e.filter(e=>e.status===`unspecified`&&Number.isInteger(e.atom)&&e.atom>=0&&e.atom<n);if(!r.length)return[];let i=Array.from({length:n},(e,t)=>t),a=e=>{for(;i[e]!==e;)e=i[e]=i[i[e]];return e},o=e=>Number.isInteger(e)&&e>=0&&e<n;for(let e of t)o(e.a)&&o(e.b)&&(i[a(e.a)]=a(e.b));let s=new Set(t.filter(e=>lp.has(e.stereo)&&o(e.a)&&o(e.b)).map(e=>a(e.a)));return[...new Set(r.map(e=>e.atom).filter(e=>s.has(a(e))))].sort((e,t)=>e-t)}function dp(e,t,n){let r=e=>Number.isInteger(e)&&e>=0&&e<n,i=new Map(e.atoms.filter(e=>r(e.atom)).map(e=>[e.atom,e])),a=new Set(e.axes.flatMap(e=>t[e.bond]?[t[e.bond].a,t[e.bond].b]:[])),o=[],s=new Map;for(let e of t)lp.has(e.stereo)&&r(e.a)&&r(e.b)&&(s.set(e.a,[...s.get(e.a)??[],e.stereo]),e.stereo!==`wavy`&&(!i.has(e.a)&&!a.has(e.a)?o.push({check:`no-centre`,atoms:[e.a]}):i.has(e.a)&&i.has(e.b)&&o.push({check:`between`,atoms:[Math.min(e.a,e.b),Math.max(e.a,e.b)]})));for(let[e,t]of s){let n=i.get(e);n&&t.length>1&&t.some(e=>e!==`wavy`)&&n.status!==`specified`&&o.push({check:`conflict`,atoms:[e]})}let c=new Set(e.groups.filter(e=>e.type===`and`||e.type===`or`).flatMap(e=>e.atoms));for(let e of i.values())e.status!==`specified`&&(c.has(e.atom)||e.group&&e.group.type!==`abs`)&&o.push({check:`group`,atoms:[e.atom]});let l=new Set;return o.filter(e=>{let t=`${e.check}:${e.atoms.join(`,`)}`;return!l.has(t)&&(l.add(t),!0)}).sort((e,t)=>e.atoms[0]-t.atoms[0])}function fp(e,t=p){return e===`no-centre`?t.problemStereoNoCentre:e===`between`?t.problemStereoBetween:e===`conflict`?t.problemStereoConflict:t.problemStereoGroup}function pp(e){let t;try{t=JSON.parse(e||`[]`)}catch{return[]}if(!Array.isArray(t))return[];let n=new Set;for(let e of t)e.type===`SUP`&&e.contracted&&Number.isInteger(e.labelAtom)&&(e.bonds?.length??0)>(e.attachmentPoints?.length??0)&&n.add(e.labelAtom);return[...n].sort((e,t)=>e-t)}var mp=.05;function hp(e,t,n,r=mp){let i=r*(Number.isFinite(n)&&n>0?n:1.5),a=Math.max(i,Number.isFinite(n)&&n>0?n:1.5),o=t=>!!e[t]&&!e[t].hidden&&Number.isFinite(e[t].x)&&Number.isFinite(e[t].y),s=(e,t)=>`${e} ${t}`,c=new Map;e.forEach((e,t)=>{if(!o(t))return;let n=s(Math.floor(e.x/a),Math.floor(e.y/a)),r=c.get(n);r?r.push(t):c.set(n,[t])});let l=(e,t,n,r)=>{let o=[];for(let l=Math.floor((e-i)/a);l<=Math.floor((n+i)/a);l++)for(let e=Math.floor((t-i)/a);e<=Math.floor((r+i)/a);e++)o.push(...c.get(s(l,e))??[]);return o},u=new Set(t.map(e=>e.a<e.b?`${e.a} ${e.b}`:`${e.b} ${e.a}`)),d=[],f=new Set;e.forEach((t,n)=>{if(o(n))for(let r of l(t.x,t.y,t.x,t.y))r<=n||u.has(`${n} ${r}`)||Math.hypot(e[r].x-t.x,e[r].y-t.y)<i&&(d.push([n,r]),f.add(n).add(r))});let p=new Set;for(let n of t){if(!o(n.a)||!o(n.b))continue;let[t,r]=[e[n.a],e[n.b]];for(let a of l(Math.min(t.x,r.x),Math.min(t.y,r.y),Math.max(t.x,r.x),Math.max(t.y,r.y)))a===n.a||a===n.b||f.has(a)||p.has(a)||_p(e[a],t,r)<i&&p.add(a)}for(let e of p)d.push([e]);return d.sort((e,t)=>e[0]-t[0]||e.length-t.length)}function gp(e,t,n){let r=new Map;t.forEach((e,t)=>{for(let n of[e.a,e.b])r.set(n,[...r.get(n)??[],t])});let i=[];for(let[a,o]of r){if(o.length!==2||!e[a]||e[a].hidden)continue;let[r,s]=o.map(e=>t[e]);if(r.type!==s.type||r.stereo!==s.stereo)continue;let[c,l]=[r,s].map(t=>e[t.a===a?t.b:t.a]);if(!c||!l||c.hidden||l.hidden)continue;let[u,d,f,p]=[c.x-e[a].x,c.y-e[a].y,l.x-e[a].x,l.y-e[a].y],m=Math.hypot(u,d),h=Math.hypot(f,p);m>0&&h>0&&((u*f+d*p)/(m*h)>-Math.cos(2*Math.PI/180)||n(a)&&i.push(a))}return i.sort((e,t)=>e-t)}function _p(e,t,n){let[r,i]=[n.x-t.x,n.y-t.y],a=r||i?Math.max(0,Math.min(1,((e.x-t.x)*r+(e.y-t.y)*i)/(r*r+i*i))):0;return Math.hypot(e.x-(t.x+a*r),e.y-(t.y+a*i))}function vp(e,t,n,r=[],i=[],a=[],o=[],s=[]){let c=[];for(let r of e){let e=[...new Set(r.atoms.filter(e=>Number.isInteger(e)&&e>=0&&e<t))].sort((e,t)=>e-t);if(!e.length)continue;let i=np(r.type);c.push({kind:i,atoms:e,message:sp(i,n(e[0]))})}for(let e=0;e<t;e++){let t=n(e);op(t)!==null&&c.push({kind:`label`,atoms:[e],message:sp(`label`,t)})}let l=new Set(s.filter(e=>e.check===`conflict`||e.check===`group`).flatMap(e=>e.atoms));for(let e of new Set(r))Number.isInteger(e)&&e>=0&&e<t&&!l.has(e)&&c.push({kind:`stereo`,atoms:[e],message:sp(`stereo`,null)});for(let e of s){let n=[...new Set(e.atoms.filter(e=>Number.isInteger(e)&&e>=0&&e<t))].sort((e,t)=>e-t);n.length&&c.push({kind:`stereo`,atoms:n,message:fp(e.check),check:e.check})}for(let e of new Set(i))Number.isInteger(e)&&e>=0&&e<t&&c.push({kind:`abbreviation`,atoms:[e],message:sp(`abbreviation`,null)});for(let e of a){let n=[...new Set(e.filter(e=>Number.isInteger(e)&&e>=0&&e<t))].sort((e,t)=>e-t);n.length&&c.push({kind:`contact`,atoms:n,message:sp(`contact`,null)})}for(let e of new Set(o))Number.isInteger(e)&&e>=0&&e<t&&c.push({kind:`linear`,atoms:[e],message:sp(`linear`,null)});return cp(c)}function yp(e){let t=new Map;for(let n of e)for(let e of n.atoms){let r=t.get(e);(!r||rp[n.kind]<rp[r.kind])&&t.set(e,n)}return t}var bp=(e,t)=>`${e.kind}:${e.atoms.map(t).sort((e,t)=>e-t).join(`,`)}`;function xp(e,t,n){let r=me(e.atomJson(t));if(!r||r.radical===0)return null;let i=n(e.toMolblock(),`molV2000`);try{return(me(i.atomJson(t))?.radical??0)<r.radical?r.symbol:null}finally{try{i.free()}catch{}}}var Sp=Object.freeze(Array.from({length:32},(e,t)=>t+1));function Cp(e){let t=new Set(e),n=1;for(;t.has(n);)n+=1;return n}function wp(e){return{opened:[...e],pressed:new Set(e.filter(e=>e>=1&&e<=32))}}function Tp(e,t){if(!Number.isInteger(t)||t<1||t>32)return e;let n=new Set(e.pressed);return n.has(t)?n.delete(t):n.add(t),{opened:e.opened,pressed:n}}function Ep(e){let t=e.opened.filter(e=>e>32);return[...new Set([...e.pressed,...t])].sort((e,t)=>e-t)}function Dp(e,t){let n={ArrowRight:1,ArrowLeft:-1,ArrowDown:8,ArrowUp:-8};if(t===`Home`)return 1;if(t===`End`)return 32;let r=e+(n[t]??0);return r>=1&&r<=32?r:e}var Op=e=>[`ArrowRight`,`ArrowLeft`,`ArrowDown`,`ArrowUp`,`Home`,`End`].includes(e),kp=(e,t)=>e.length===t.length&&e.every((e,n)=>e===t[n]),Ap=e=>[...new Set(e)].sort((e,t)=>e-t);function jp(e,t,n){let r=Ap(t),i=Ap(n);return kp(r,i)?null:S(`atom.rgroup`,[{op:`setAtomRNumbers`,atom:e,numbers:i}])}function Mp(e,t){let n=Ap(t);return n.length?S(`atom.rgroup`,[{op:`addAtom`,symbol:`*`,x:e.x,y:e.y},{op:`setAtomRNumbers`,atom:{added:0},numbers:n}]):null}function Np(e,t,n){let r=Ap(n);return r.length?`point`in t?S(`atom.rgroup`,[{op:`addAtom`,symbol:`*`,x:t.point.x,y:t.point.y},{op:`addBond`,a:e,b:{added:0},type:`single`,stereo:`none`},{op:`setAtomRNumbers`,atom:{added:0},numbers:r}]):t.atom===e?null:S(`atom.rgroup`,[...t.bonded?[]:[{op:`addBond`,a:e,b:t.atom,type:`single`,stereo:`none`}],...kp(Ap(t.was),r)?[]:[{op:`setAtomRNumbers`,atom:t.atom,numbers:r}]]):null}function Pp(e,t,n){return jp(e,t,[Cp(n)])}function Fp(e,t,n){return Np(e,t,[Cp(n)])}var Ip=e=>({primary:!!(e&1),secondary:!!(e&2)}),Lp=(e,t)=>!!e|(t?2:0);function Rp(e,t,n){return(t&3)==(n&3)?null:S(`atom.attachment-points`,[{op:`setAtomAttachmentPoints`,atom:e,mask:n&3}])}var zp=Object.freeze({atoms:Object.freeze([]),bonds:Object.freeze([])});function Bp(e,t,n){return n.x>=Math.min(e.x,t.x)&&n.x<=Math.max(e.x,t.x)&&n.y>=Math.min(e.y,t.y)&&n.y<=Math.max(e.y,t.y)}function Vp(e,t){if(e.length<3)return!1;let n=!1;for(let r=0,i=e.length-1;r<e.length;i=r++){let[a,o]=[e[r],e[i]];a.y>t.y!=o.y>t.y&&t.x<(o.x-a.x)*(t.y-a.y)/(o.y-a.y)+a.x&&(n=!n)}return n}var Hp=(e,t,n)=>(t.x-e.x)*(n.y-e.y)-(t.y-e.y)*(n.x-e.x);function Up(e,t,n,r){let[i,a,o,s]=[Hp(n,r,e),Hp(n,r,t),Hp(e,t,n),Hp(e,t,r)];return(i>0&&a<0||i<0&&a>0)&&(o>0&&s<0||o<0&&s>0)}function Wp(e,t,n){if(!Vp(e,t)||!Vp(e,n))return!1;for(let r=0,i=e.length-1;r<e.length;i=r++)if(Up(t,n,e[i],e[r]))return!1;return!0}function Gp(e,t){let n=t.kind===`rect`?e=>Bp(t.from,t.to,e):e=>Vp(t.path,e),r=e.atoms.flatMap((e,t)=>n(e)?[t]:[]),i=new Set(r);return{atoms:r,bonds:e.bonds.flatMap((n,r)=>!i.has(n.a)||!i.has(n.b)||t.kind===`lasso`&&!Wp(t.path,e.atoms[n.a],e.atoms[n.b])?[]:[r])}}function Kp(e,t){let n=e.atoms.map(()=>[]);e.bonds.forEach((e,t)=>{n[e.a]?.push(t),n[e.b]?.push(t)});let r=new Set,i=[...t.atoms,...t.bonds.flatMap(t=>e.bonds[t]?[e.bonds[t].a,e.bonds[t].b]:[])].filter(t=>e.atoms[t]);for(let e of i)r.add(e);for(;i.length;){let t=i.pop();for(let a of n[t]??[]){let{a:n,b:o}=e.bonds[a],s=n===t?o:n;r.has(s)||(r.add(s),i.push(s))}}return{atoms:[...r].sort((e,t)=>e-t),bonds:e.bonds.flatMap((e,t)=>r.has(e.a)&&r.has(e.b)?[t]:[])}}function qp(e){let t=-1;return e.atoms.forEach((n,r)=>{(t<0||n.id>e.atoms[t].id)&&(t=r)}),t<0?zp:Kp(e,{atoms:[t],bonds:[]})}var Jp=e=>({atoms:e.atoms.map((e,t)=>t),bonds:e.bonds.map((e,t)=>t)}),Yp=class{#e=new Set;#t=new Set;get empty(){return this.#e.size===0&&this.#t.size===0}ids(){return{atoms:[...this.#e].sort((e,t)=>e-t),bonds:[...this.#t].sort((e,t)=>e-t)}}atomsInOrder(e){let t=new Map(e.atoms.map((e,t)=>[e.id,t]));return[...this.#e].flatMap(e=>t.has(e)?[t.get(e)]:[])}setIds(e,t){this.clear();for(let t of e)this.#e.add(t);for(let e of t)this.#t.add(e)}moved(e){let{atoms:t,bonds:n}=this.indices(e),r=new Set(t);for(let t of n)r.add(e.bonds[t].a).add(e.bonds[t].b);return[...r].sort((e,t)=>e-t)}clear(){let e=!this.empty;return this.#e.clear(),this.#t.clear(),e}replace(e,t){this.clear(),this.add(e,t)}add(e,t){for(let n of t.atoms)e.atoms[n]&&this.#e.add(e.atoms[n].id);for(let n of t.bonds)e.bonds[n]&&this.#t.add(e.bonds[n].id)}remove(e,t){for(let n of t.atoms)e.atoms[n]&&this.#e.delete(e.atoms[n].id);for(let n of t.bonds)e.bonds[n]&&this.#t.delete(e.bonds[n].id)}under(e,t){if(this.empty)return!1;let n=new Set(this.moved(e));if(t.kind===`atom`)return n.has(t.index);let r=e.bonds[t.index];return!!r&&(this.#t.has(r.id)||n.has(r.a)&&n.has(r.b))}holds(e,t){let n=t.atoms.filter(t=>e.atoms[t]),r=t.bonds.filter(t=>e.bonds[t]);return n.length+r.length>0&&n.every(t=>this.#e.has(e.atoms[t].id))&&r.every(t=>this.#t.has(e.bonds[t].id))}toggle(e,t){this.holds(e,t)?this.remove(e,t):this.add(e,t)}indices(e){return{atoms:e.atoms.flatMap((e,t)=>this.#e.has(e.id)&&!e.hidden?[t]:[]),bonds:e.bonds.flatMap((e,t)=>this.#t.has(e.id)&&!e.hidden?[t]:[])}}prune(e){let t=new Set(e.atoms.map(e=>e.id)),n=new Set(e.bonds.map(e=>e.id));for(let e of[...this.#e])t.has(e)||this.#e.delete(e);for(let e of[...this.#t])n.has(e)||this.#t.delete(e)}},Xp=e=>({shift:e.shiftKey,alt:e.altKey,ctrl:e.ctrlKey||e.metaKey}),Zp=4e3,Qp=1e4,$p=.4,em=(e,t,n)=>n>=0?e.thisAtom(t||`?`):t||`?`,tm=(e,t)=>e.atoms[t]?.element||e.info(t)?.symbol||`?`;function nm(e,t,n){let r=n*Math.PI/180,[i,a]=[Math.cos(r),Math.sin(r)],[o,s]=[e.x-t.x,e.y-t.y];return{x:t.x+o*i+s*a,y:t.y-o*a+s*i}}var rm=(e,t)=>e===null||t===null?e===t:Math.abs(e.x-t.x)<.01&&Math.abs(e.y-t.y)<.01,im=(e,t)=>Math.atan2(-(t.y-e.y),t.x-e.x)*180/Math.PI;function am(e,t,n){let r=ye(e);return[r,0,0,-r,t/2-r*e.cx,n/2+r*e.cy]}var om=(e,t)=>!!(e&&t&&e.kind===t.kind&&e.id===t.id),sm=e=>e?.kind===`bond`||e?.kind===`element`,cm=(e,t)=>!t.grip&&!t.still&&(e?.kind===`hand`||e?.kind===`select`&&t.hit!==null),lm=new Set([`zoom.in`,`zoom.out`,`zoom.reset`,`zoom.fit`,`select.all`,`select.none`,`select.last`,`select.hotspot`,`selection.end`,`copy`,`copy.as`,`hotspot.walk`,`hotspot.jump`,`hotspot.hop`,`hotspot.start`,`keys.open`,`settings.open`]),um=(e,t)=>!t.grip&&t.hit===null&&(t.group??-1)<0&&(e.kind===`select`&&e.mode!==`structure`||e.kind===`erase`||e.kind===`sgroup`),dm=new Set,fm=(e,t)=>({atoms:[...new Set([...e.atoms,...t.atoms])].sort((e,t)=>e-t),bonds:[...new Set([...e.bonds,...t.bonds])].sort((e,t)=>e-t)});function pm(e,t){return{...e,atoms:e.atoms.map(e=>{let n=t(e),r=e.label?{...e.label,x:e.label.x+n.x-e.x,y:e.label.y+n.y-e.y}:e.label;return{...e,...n,label:r}}),bonds:e.bonds.map(([e,n])=>[t(e),t(n)])}}var mm=e=>({symbol:`R${e}`,colour:`#000000`}),hm=(e,t)=>e.atoms[t]?.element===`R#`,gm=e=>Math.max(12,10/20*e),_m=e=>({atoms:e.atoms.map(e=>e.hidden?lr:e.px),bonds:e.hiddenAtoms.size?e.bonds.map(e=>e.hidden?vm:{a:e.ends[0],b:e.ends[1]}):e.bonds}),vm=Object.freeze({a:-1,b:-1}),ym=Object.freeze({}),bm=Object.freeze({slop:Po}),xm=(e,t)=>(e.info(t)?.hCount??0)>0,Sm=class{#e;#t;#n;#r;#i=Go;#a=p;get strings(){return this.#a}set strings(e){this.#a=e}get keymap(){return this.#i}set keymap(e){this.#i=e}eraseWhole=!1;hoverKeys=!1;set pointerKind(e){this.#o=e===`mouse`?ym:bm}#o=ym;#s=Gt;#c=null;#l=null;#u=null;#d=null;#f=null;#p=null;#m=new Set;#h=null;#g=new po;#_=null;#v=-1;#y=null;#b=null;#x=null;#S=null;#C=[];#w=new Ta;#T=null;#E=null;#D;#O=!1;#k={};#A=!1;#j=!1;#M=null;#N=null;#P=null;#F=null;#I=null;#L=[];#R=null;#z;#B=null;#V=null;#H=null;#U=-1/0;#W=null;#G=null;#K=null;#q={message:null,pointer:!1};#J=!0;#Y=null;#X=[];#Z=null;#Q=new Yp;#$=Dn;#ee=[];#te=null;#ne=null;#re=null;#ie=null;#ae=!1;#oe=null;#se=null;#ce=0;#le=null;#ue=null;#de=null;#fe=0;#pe=!1;#me=JSON.stringify({atoms:[],bonds:[]});#he=null;#ge=null;#_e=null;#ve=null;#ye=!1;constructor(e,t,n,r,i=()=>performance.now(),a=co){this.#e=e,this.#t=t,this.#n=n,this.#r=r,this.#D=i,this.#z=a}get tool(){return this.#s}get element(){let e=ot(this.#s);return e?.kind===`element`?e.symbol:null}get viewOnly(){return this.#be}set viewOnly(e){if(e===this.#be||this.#ye)return;this.#be=e,e&&(this.cancelGesture(),this.#zt());let t=this.#xr();t&&(this.#Xn(t),this.#Hr(t))}#be=!1;#xe(){let e=ot(this.#s);return this.#be?e?.kind===`select`?e:ot(`select.rect`):e?.kind===`bond`?this.#Se(e):e?.kind===`query-atom`?this.#p?{kind:`element`,symbol:`*`,query:this.#p}:null:e}#Se(e){return this.#Ce&&e.type===`aromatic`&&!e.query?{...e,query:`aromatic`}:e}get#Ce(){return this.#e.mode===`query`}get queryAtom(){return this.#p}setQueryAtom(e){if(this.#Ce){if(this.#p=e,this.#s===`element.query`){this.#Tn(),this.#Yt(),this.#Kt(),this.#r.toolChanged(br);return}this.#De(br)}}set queryOff(e){this.#m=e,tn.has(this.#s)&&this.#we(this.#s)&&this.#De(Gt)}#we(e){let t=ot(e);if(t?.kind===`bond`)return!!(t.query&&this.#m.has(t.query));let n=e===`element.query`&&this.#p?Ce(this.#p):null;return n!==null&&this.#m.has(n)}#Te(e,t){return e!==null&&this.#m.has(e)?this.#a.format(`queryElementOff`,{name:t}):null}modeChanged(){!this.#Ce&&tn.has(this.#s)&&this.#De(Gt)}get template(){return this.#c}setTemplate(e){this.#c=e,this.#l=ep(e),this.#Ee(`template`),e.unplaceable&&this.notify(this.#mt(e))}setGroup(e){this.#d=e,this.#Ee(`group`)}setSalt(e){this.#f=e,this.#Ee(`salt`)}get picked(){return this.#u}#Ee(e){if(this.#u=e,this.#s===`template.library`){this.#Tn(),this.#Qt(),this.#Yt(),this.#Kt(),this.#r.toolChanged(Ar);return}this.#De(Ar)}setTool(e){if(!Xr(e))return;let t=ot(e);if(t?.kind===`element`&&this.#Fn(t.symbol)||e===this.#s||e===`template.library`&&!this.#u||tn.has(e)&&(!this.#Ce||e===`element.query`&&!this.#p||this.#we(e)))return;if(e===`erase`&&!this.#Q.empty&&!this.#S){this.deleteSelection();return}let n=e===`sgroup`&&!this.#Q.empty&&!this.#S?this.#xr():null,r=n?this.#Q.indices(n):null;this.#De(e),n&&r&&(r.atoms.length||r.bonds.length)&&this.#jr(n,r.atoms,r.bonds)}#De(e){this.#s=e,Fr(e)&&(this.#$=e),this.#Tn(),this.#w.forget(),this.#St(),e!==`hand`&&this.#In(),this.#Qt(),this.#Yt();let t=this.#xr();t&&this.#Hr(t),this.#Kt(),this.#r.toolChanged(e)}afterRender(){this.#it+=1,this.#_=null,this.#Br(),this.#Dr(),this.#Lr(),this.#C=[],this.#ee=[],this.#X=[],this.#yr=null,this.#ne=null,this.#L=[],this.#St(!1),this.#Qt(),this.#F!==null&&this.#Yt(),this.#I=null,this.#e.generation!==this.#v&&(this.#v=this.#e.generation,this.#x=null,this.#b=null,this.#N=null,this.#P=null,this.#Q.clear(),this.#Tn(),this.#w.forget(),this.#zt());let e=this.#xr();if(!e){this.#n.hideHover(),this.#n.hideSelection(),this.#n.hideProblems(),this.#pr();return}this.#Q.prune(e),this.#fr(e),this.#wt(e);let t=this.#P===null?void 0:e.atoms[this.#P];this.#P=null,t&&(this.#x={kind:`atom`,id:t.id},this.#N=this.#y?{...this.#y}:null),this.#N&&!this.#Sr(e)&&(this.#N=null),this.#y&&!this.#S?(this.#Cr(e,this.#y),this.#sn()&&this.#en()):this.#Hr(e),this.#y&&!this.#S&&this.#sn()||this.#xt(),this.#Jt(e);let n=this.#nt;this.#nt=null,n&&this.#rt(n)}pointerDown(e,t,n=!1,r={}){if(this.#S||this.#T)return!1;let i=this.#xr();if(!i)return!1;this.#Br(),this.#Lr(),this.#W&&!this.#W.pointer&&this.#bt(null),this.#y=e,this.#k=r,r.alt&&(this.#A=!0),this.#O=!1;let a=i.toMolecule(e),o=this.#oe?Fo(i.hitAtoms(),i.hitBonds(),e,i.bondPx,this.#o):null,s=this.#Qn(i,e,n,!!r.alt,o);if(s){this.#w.forget(),this.#b=null,this.#Vr();let r=s===`rotate`?this.#$n(i,e):null;return this.#M=null,this.#S={pointerId:t,start:e,hit:o,dragging:!1,cancelled:!1,count:1,touch:n,path:[e],startMol:a,pathMol:[a],moving:null,whole:!1,copy:!1,pressSelected:!1,grip:s,turn:r},this.#ne=null,!0}let c=this.#Cr(i,e);this.#Dr();let l=this.#wr(i,e,c,n);l>=0&&(c=null);let u=n?this.#w.tapCount(e,this.#D()):1,d={pointerId:t,start:e,hit:c,dragging:!1,cancelled:!1,count:u,touch:n,path:[e],startMol:a,pathMol:[a],moving:null,whole:!1,copy:!1,pressSelected:!1,grip:null,turn:null,group:l,still:this.#be},f=this.#xe();if(c&&(f?.kind===`select`||f?.kind===`hand`)){let e=this.#Q.under(i,c);!e&&f.kind===`select`&&(this.#zn(i,c,!!r.shift,f.mode===`structure`,!1),d.pressSelected=!0,e=this.#Q.under(i,c)),d.whole=e,d.moving=i.wholeGroups(e?this.#Q.moved(i):c.kind===`atom`?[c.index]:[i.bonds[c.index].a,i.bonds[c.index].b]),d.copy=e&&!!r.ctrl}return this.#S=d,this.#M=null,this.#ne=null,!0}clickCount(e){let t=this.#S;t&&!t.touch&&Number.isInteger(e)&&e>0&&(t.count=e)}pointerMove(e,t,n={}){let r=this.#S;if(r&&t!==void 0&&t!==r.pointerId)return;this.#y&&(this.#y.x!==e.x||this.#y.y!==e.y)&&(this.#nt=null),this.#y=e,this.#k=n,r&&n.alt&&(this.#A=!0),this.#U=this.#D(),this.#O=!1;let i=this.#xr();if(!i)return;if(!r){this.#j=!1,this.#W&&!this.#W.pointer&&!this.#dt()&&this.#bt(null),this.#Cr(i,e),this.#Kt();return}if(r.cancelled){this.#Cr(i,e);return}if(r.path.push(e),r.pathMol.push(i.toMolecule(e)),!r.dragging&&Math.hypot(e.x-r.start.x,e.y-r.start.y)>=4&&(r.dragging=!0),!r.dragging)return;if(r.grip){this.#$t();return}let a=this.#xe();if(cm(a,r)){this.#$t(),this.#ir();return}if(a&&um(a,r)){this.#lr(i,r,e),this.#ir();return}if(a?.kind===`select`)return;if(a?.kind===`erase`||a?.kind===`charge`){this.#Vr();return}if(this.#sn()||a?.kind===`chain`){this.#$t(),this.#Hr(i);return}if(a?.kind===`rgroup`&&r.hit?.kind!==`atom`){this.#Yt(),this.#Vr();return}let o=a?.kind===`rgroup`?null:this.#Gt();this.#F=null;let s=this.#cn();if(o&&(r.hit?.kind!==`atom`&&(a?.kind===`element`||s?.kind===`group`)||s?.kind===`salt`)){this.#qt(e,o,i.bondPx),this.#Hr(i);return}let c=this.#Jr(i,r,e);this.#R=null;let l=c&&a?this.#_t(`drag:${JSON.stringify([c.from,c.to])}`,()=>this.#gt(i,a,c)):null;this.#bt(l);let u=a?.kind===`rgroup`&&c?mm(Cp(this.#Nt(i,r,`atom`in c.to?c.to.atom:-1))):o;c&&!l?this.#n.showPreview(c.fromPx,c.toPx,i.bondPx,u??void 0):this.#Yt();let d=c&&`atom`in c.to?c.to.atom:-1;this.#Hr(i,d>=0?{kind:`atom`,index:d}:void 0)}pointerUp(e,t,n={}){let r=this.#S;if(!r||t!==void 0&&t!==r.pointerId)return;this.#S=null,this.#M=null,this.#or(),this.#Qt(),this.#Yt(),this.#n.hideArea(),this.#te=null,this.#ne=null,r.dragging&&this.#St(),this.#k=n,n.alt&&(this.#A=!0),r.path.push(e),r.touch&&(this.#j=!0),this.#y=r.touch||this.#Sn(e)?e:null,this.#y||(this.#b=null);let i=this.#xr();if(!i)return;r.pathMol.push(i.toMolecule(e));let a=!1;if(r.grip&&!r.dragging&&!r.cancelled&&r.hit&&(r.grip=null,r.turn=null,this.#de=null,this.#fe=0),r.grip)this.#w.forget(),!r.cancelled&&r.dragging&&r.grip===`rotate`&&r.turn?a=this.#st(ju(i,r.turn.atoms,r.turn.centre,this.#er(i,r,e)),i,[],`bond`,!0,!0):!r.cancelled&&r.dragging&&r.grip===`centre`&&this.#rr(i,e),this.#de=null,this.#fe=0;else if(!r.cancelled){let t=this.#xe();if(!r.dragging&&r.count>=2&&r.hit?.kind===`atom`&&sm(t)){this.#At(i,i.atoms[r.hit.index].id);return}if(!r.dragging&&(r.group??-1)>=0&&!this.#be){this.#Or(i,r,t);return}let n={plan:null,noRoom:[]},o=`bond`,s=this.#sn(),c=null;if(t?.kind===`select`&&!(r.dragging&&r.hit))r.dragging?um(t,r)&&this.#ur(i,r,e):((!r.pressSelected||r.count>=2)&&this.#Rn(i,r),r.count>=2&&!r.hit&&this.#r.zoom(`fit`));else if(!r.still){if(cm(t,r)){if(r.dragging){let{plan:t}=this.#fn(i,r,e);n.plan=t.tx,c=r.whole?this.#Bn(i,r,t):null}}else if(t?.kind===`erase`&&r.dragging&&um(t,r))n.plan=this.#dr(i,r,e);else if(this.#cn())n.plan=this.#un(i,r,e,this.#cn());else if(s)n=E(()=>this.#dn(s,i,r,e)),o=t?.kind===`template`?`template`:`ring`;else if(t?.kind===`chain`){if(r.dragging){let t=E(()=>this.#yn(i,r,e));n={plan:t.plan?.tx??null,noRoom:t.noRoom},o=`chain`}else if(r.hit?.kind===`bond`){let e=r.hit.index;n=E(()=>g(t,i,e,this.#at))}else{let e=r.hit?.kind===`atom`?{atomId:i.atoms[r.hit.index].id,last:this.#vn}:{point:i.toMolecule(r.start),last:this.#vn};if(this.#k.alt&&this.#vn!==null){let e=Jd(i,r.hit?.kind===`atom`?{atom:r.hit.index}:{point:i.toMolecule(r.start)},this.#vn);n={plan:e?.tx??null,noRoom:!e&&r.hit?.kind===`atom`?[{atom:r.hit.index,symbol:tm(i,r.hit.index)}]:[]},o=`chain`}else this.#r.openAddChain?.(e)}}else if(t?.kind===`rgroup`||t?.kind===`attachment`)r.dragging?t.kind===`rgroup`&&(n.plan=this.#Mt(i,r,e)):this.#jt(t.kind,i,r);else if(t?.kind===`sgroup`)r.dragging?um(t,r)&&this.#Ar(i,r,e):this.#kr(i,r);else if(t&&r.dragging){if(t.kind===`element`&&r.hit?.kind!==`atom`)n.plan=b(t.symbol,i.toMolecule(e),this.#wn(i,e),i,t.query);else{let a=this.#Jr(i,r,e);a&&(n=E(()=>ie(t,i,a.from,a.to,this.#at)))}}else if(t){let e=r.hit;n=E(()=>e?e.kind===`atom`?lt(t,i,e.index,this.#at):g(t,i,e.index,this.#at,this.eraseWhole):Ur(t,i,i.toMolecule(r.start)))}}if(a=this.#st(n.plan,i,n.noRoom,o,!0),a&&c&&this.#Vn(c),r.dragging)this.#w.forget();else{let e=r.hit?.kind===`atom`?i.atoms[r.hit.index].id:null;this.#w.clicked({atomId:e,token:a?this.#e.stateToken:null,at:r.start,time:this.#D(),touch:r.touch})}}a||(this.#y?this.#Cr(i,this.#y):(this.#b=null,this.#Hr(i)),this.#fr(i),this.#sn()&&this.#$t())}get busy(){return this.#S!==null}positions(e,t=null){return uo(this.#ye?null:this.#xr(),e,t)}grabs(e,t=!1,n=!1){if(this.#s!==`hand`||this.#T)return!1;let r=this.#xr();if(!r)return!1;let i=Fo(r.hitAtoms(),r.hitBonds(),e,r.bondPx,this.#o);return i!==null||this.#Qn(r,e,t,n,i)!==null}navigating(){this.#Tn(),this.#y=null,this.#b=null,this.#Br(),this.#Dr(),this.#O=!1,this.#Lr(),this.#St(),this.#Qt(),this.#Yt();let e=this.#xr();e&&this.#Hr(e)}viewMoving(){this.#Br(),this.#Dr(),this.#Lr(),this.#Qt(),this.#Yt()}pointerCancel(){this.#Tn(),this.#St();let e=this.#xr();e&&this.#Hr(e)}forgetRefusal(){if(this.#ye||this.#S||!this.#W)return;this.#St();let e=this.#xr();e&&this.#Hr(e)}pointerLeave(){if(this.#S)return;this.#y=null,this.#b=null,this.#Br(),this.#Dr(),this.#W?.notice||this.#dt()?this.#K=null:this.#St(),this.#Lr(),this.#Qt(),this.#Yt();let e=this.#xr();e&&this.#Hr(e)}dragElement(e,t){let n=this.#xr();if(!t||!n){let e=this.#E;this.#E=null,this.#Yt(),e&&this.#St(),e&&n&&this.#Hr(n);return}this.#E={symbol:e,at:t},this.#F=null,this.#R=null,this.#n.showLabel(t,this.#Wt(e),n.bondPx);let r=this.#wn(n,t);r>=0?this.#Hr(n,{kind:`atom`,index:r}):this.#Vr()}dropElement(e,t){this.#E=null,this.#Yt();let n=this.#xr();if(!n)return!1;this.#w.forget(),this.#St();let r=this.#st(b(e,n.toMolecule(t),this.#wn(n,t),n),n);return r||this.#Hr(n),r}cancelDrop(){if(!this.#E)return;this.#E=null,this.#Yt(),this.#St();let e=this.#xr();e&&this.#Hr(e)}keyDown(e,t=!1){this.#nt=null;let n=this.#S;if(e.key===`Escape`)return this.#g.reset(),this.cancelGesture();if(n&&e.key===`Alt`)return this.#A=!0,this.#Oe({...this.#k,alt:!0}),!0;if(n)return!1;let r=this.#xr(),i=r&&(!this.hoverKeys||this.#ot()||this.#O)?this.#Sr(r):null,a=this.#b!==null;if(r&&!this.#Q.empty&&!a&&!this.#be){let t=this.#Q.indices(r),n=ss(this.#i,e,t.atoms.length>0,t.bonds.length>0);if(n)return this.#w.forget(),this.#W&&!this.#W.pointer&&this.#bt(null),this.#Pn(r,n,t),!0}let o=(r&&!this.#Q.empty?as(this.#i,e,t,a):null)??is(this.#i,e,i?i.kind:null,this.#ot(),t);if(!this.#be&&this.#ke(r,e,o,a))return!0;if(!o||this.#be&&!lm.has(o.action))return!1;this.#w.forget(),this.#W&&!this.#W.pointer&&this.#bt(null);let s=!!(r&&i&&o.scope===i.kind);s&&(this.#O=!0);let c=this.#s,l=this.#Ae(o,r,i);return s&&r&&!this.#S&&this.#xr()===r&&this.#Hr(r),this.#s!==c&&!this.#ye&&this.#r.announce?.(this.#a.format(`toolChosen`,{tool:D(this.#s,this.#a)}),`tool`),l}keyUp(e){return e.key!==`Alt`||!this.#A?!1:(this.#A=!1,this.#S&&this.#Oe({...this.#k,alt:!1}),!0)}#Oe(e){let t=this.#S;if(!t||!this.#y){this.#k=e;return}this.pointerMove(this.#y,t.pointerId,e)}#ke(e,t,n,r){let i=fo(t);if(i===null||e===null||!this.#Q.empty||r||this.#T||n&&n.scope!==`generic`||!this.#r.openLookup)return this.#g.reset(),!1;let a=this.#g.press(i,this.#D(),!n);return a!==null&&(this.#w.forget(),this.#r.openLookup({mode:`lookup`,text:a,at:this.#y?{...this.#y}:null}),!0)}applyNickname(e,t){let n=this.#xr(),r=n?.atomById(e)??-1;return!n||r<0||this.#S||this.#ye?!1:(this.#w.forget(),this.#st(ho(n,r,t),n))}#Ae(e,t,n){switch(e.action){case`clear`:return this.clear(),!0;case`undo`:return this.#e.undo(`user`),!0;case`redo`:return this.#e.redo(`user`),!0;case`tool.bond.single`:return this.setTool(`bond.single`),!0;case`tool.ring.benzene`:return this.setTool(`ring.benzene`),!0;case`tool.ring.cyclopentadiene`:return this.setTool(`ring.cyclopentadiene`),!0;case`tool.chain`:return this.setTool(`chain`),!0;case`tool.hand`:return this.setTool(`hand`),!0;case`tool.choose`:return e.tool&&this.setTool(e.tool),!0;case`tool.cycle`:{let t=e.tools??[],n=t.indexOf(this.#s);return t.length&&this.setTool(t[(n+1)%t.length]),!0}case`keys.open`:return this.#r.openKeys?.(),!0;case`settings.open`:return this.#r.openSettings?.(),!0;case`library.open`:return this.#r.openLibrary?.(),!0;case`atom.bond`:if(t&&n?.kind===`atom`){let r=E(()=>lt(this.#Se({kind:`bond`,type:e.bondType??`single`,stereo:`none`}),t,n.index,this.#at));this.#st(r.plan,t,r.noRoom,`bond`)}return!0;case`select.all`:return this.selectAll(),!0;case`select.none`:return this.#In(),!0;case`select.last`:return this.#Ln(),!0;case`selection.delete`:return this.deleteSelection(),!0;case`selection.nudge`:return t&&e.nudge&&this.#Hn(t,e.nudge),!0;case`select.join`:return t&&this.#Wn(t),!0;case`selection.rotate`:return t&&e.turn&&this.#Gn(t,e.turn),!0;case`flip.horizontal`:return this.flip(`horizontal`,`mirror`),!0;case`structure.aromatize`:case`structure.dearomatize`:return this.aromaticity(e.action===`structure.aromatize`?`aromatize`:`dearomatize`),!0;case`flip.vertical`:return this.flip(`vertical`,`mirror`),!0;case`rotate180.horizontal`:return this.flip(`horizontal`,`turn`),!0;case`rotate180.vertical`:return this.flip(`vertical`,`turn`),!0;case`cut`:case`copy`:case`paste`:return!1;case`copy.as`:return e.copy&&this.#r.copyAs?.(e.copy),!0;case`paste.as`:return e.paste&&this.#r.pasteAs?.(e.paste),!0;case`zoom.in`:return this.#r.zoom(`in`),!0;case`zoom.out`:return this.#r.zoom(`out`),!0;case`zoom.reset`:return this.#r.zoom(`reset`),!0;case`zoom.fit`:return this.#r.zoom(`fit`),!0;case`layout`:case`clean`:return this.#r.relayout?this.#r.relayout(e.action):this.relayout(e.action),!0;case`atom.nickname`:if(t&&n?.kind===`atom`&&this.#r.openLookup){let e=t.atoms[n.index];this.#r.openLookup({mode:`nickname`,text:``,at:{...e.px},atomId:e.id})}return!0;case`atom.abbreviation`:return t&&n?.kind===`atom`&&e.abbreviation&&this.#st(this.#Vt(t,n.index,e.abbreviation),t),!0;case`atom.properties`:if(t&&n?.kind===`atom`&&t.contractedGroupOf(n.index))return!0;if(t&&n?.kind===`atom`){let e=t.atoms[n.index].id;hm(t,n.index)&&this.#r.openRGroup?this.#r.openRGroup({atomId:e}):this.#r.openAtomProperties(e)}return!0;case`atom.rgroup`:return t&&n?.kind===`atom`&&!t.contractedGroupOf(n.index)&&this.#Pt(t,n.index),!0;case`atom.generic`:{if(!this.#Ce||!e.generic)return!1;let r=this.#Te(e.generic,e.generic);return r?this.notify(r):t&&n?.kind===`atom`&&this.#st(ri(t,n.index,{kind:`generic`,generic:e.generic}),t),!0}case`bond.properties`:return t&&n?.kind===`bond`&&this.#r.openBondProperties?.([t.bonds[n.index].id]),!0;case`menu.open`:return this.#r.openMenu?.(),!0;case`atom.edit`:return t&&n?.kind===`atom`&&this.#Rt(t,n.index,!1),!0;case`hotspot.walk`:case`hotspot.jump`:case`hotspot.hop`:if(t&&n&&e.arrow){let r={kind:n.kind,index:n.index},i=_m(t),a=e.action===`hotspot.walk`?ec(i,r,e.arrow):e.action===`hotspot.jump`?tc(i,r,e.arrow):rc(i,r,e.arrow);a?this.#En(t,a):this.#r.announce?.(this.#a.nothingThatWay,`hotspot`)}return!0;case`hotspot.start`:return t?this.#On(t,!!e.arrow):!1;case`structure.new`:return t?this.#kn(t):!1;case`selection.end`:return t&&this.#jn(t),!0;case`select.hotspot`:return!t||!n||this.#ot()?!1:(this.#Mn(t,n),!0);case`atom.attachment`:return t&&n?.kind===`atom`&&!t.contractedGroupOf(n.index)&&this.#Nn(t,n.index),!0;case`atom.element`:case`atom.deuterium`:case`atom.unlabel`:return t&&n?.kind===`atom`&&this.#st(Rr(e,t,n.index),t),!0;case`atom.sprout`:return t&&n?.kind===`atom`&&e.sprout&&this.#bn(e.sprout,t,n.index),!0;case`atom.charge.plus`:case`atom.charge.minus`:return t&&n?.kind===`atom`&&this.#st($r(t,n.index,e.action===`atom.charge.plus`?1:-1),t),!0;case`bond.fuse`:if(t&&n?.kind===`bond`&&e.ring){let r=e.ring,i=E(()=>Vd(r,t,n.index));this.#st(i.plan,t,i.noRoom,`ring`)}return!0;default:if(t&&n?.kind===`bond`){let r=()=>this.#st(Bn(e.action,t,n.index,this.#at),t);e.action===`bond.left`||e.action===`bond.centre`||e.action===`bond.right`?this.#tt(t,[n.index],r):r()}return!0}}cancelGesture(){let e=this.#S;if(!e||e.cancelled)return!1;e.cancelled=!0,this.#or(),this.#Qt(),this.#Yt(),this.#St(),this.#n.hideArea(),this.#te=null,this.#ne=null,e.dragging&&(this.#b=null);let t=this.#xr();return t&&(this.#Hr(t),this.#fr(t)),!0}historyMoved(){this.cancelGesture(),this.#Q.clear(),this.#te=null,this.#ne=null,this.#ue=null}keyboardFocus(){if(this.#O||this.#ye)return;this.#O=!0;let e=this.#xr();e&&!this.#S&&this.#Hr(e)}keyboardGone(){if(!this.#O||this.#ye)return;this.#O=!1;let e=this.#xr();e&&!this.#S&&this.#Hr(e)}get keyboard(){return this.#O}clearSticky(){if(!this.#x||this.#S||this.#ot())return!1;this.#x=null;let e=this.#xr();return e?this.#Hr(e):this.#n.hideHover(),!0}escape(){if(this.#S||this.#ye)return!1;if(!this.#Q.empty)return this.#In(),!0;let e=!!this.#x&&!this.#ot();if(e&&this.#O)return this.clearSticky();let t=e&&this.clearSticky();return Fr(this.#s)?t:(this.setTool(this.#$),!0)}get lastSelectTool(){return this.#$}selectionCentre(){let e=this.#xr();if(!e||this.#Q.empty)return null;let t=this.#Q.moved(e).map(t=>e.atoms[t].at);if(!t.length)return null;let[n,r]=[t.map(e=>e.x),t.map(e=>e.y)];return{x:(Math.min(...n)+Math.max(...n))/2,y:(Math.min(...r)+Math.max(...r))/2}}get selected(){let e=this.#xr();return e?this.#Q.indices(e):zp}get selectedAtoms(){let e=this.#xr();if(!e||this.#Q.empty)return[];let t=this.#Q.indices(e),n=new Set(t.atoms);for(let r of t.bonds){let t=e.bonds[r];t&&n.add(t.a).add(t.b)}return e.wholeGroups([...n])}selectAll(){if(this.#S||this.#ye)return;Fr(this.#s)||this.setTool(this.#$);let e=this.#xr();e&&(this.#Q.replace(e,Jp(e)),this.#fr(e))}deleteSelection(){let e=this.#xr();if(!e||this.#S||this.#Q.empty)return!1;let t=this.#Q.indices(e),n=this.#st(mt(e,t.atoms,t.bonds),e);return n&&this.#In(),n}flip(e,t){let n=this.#xr();if(!n||this.#S||this.#T||this.#ye)return!1;let r=this.#Q.empty,i=r?n.atoms.map(e=>e.index):n.wholeGroups(this.#Q.moved(n));if(!i.length)return!1;let a=Fu(n,i,e,t,this.#pn().wedge??(e=>n.wedgeOf(e)),r,e=>t===`turn`&&Cm(this.#e.stereo(),n,e));if(a.refused){let e=a.reason===`ez`?this.#a.cannotTurnEz:this.#a.cannotFlip(a.attachments.length,t);return this.#lt({atoms:[],bonds:a.attachments,message:e,pointer:!1,subject:this.#yt(n)},n),!1}return this.#st(a.tx,n)}get hasSelection(){return!this.#Q.empty}copySelection(){let e=this.#xr();if(!e||this.#S||this.#T||this.#ye||this.#Q.empty)return null;let{atoms:t}=vu(e,this.#Q.indices(e));return t.length?this.#e.copyText(t,this.#e.edited?[]:bn(e)):null}copyPicture(){let e=this.#xr();if(!e||this.#S||this.#T||this.#ye||this.#Q.empty)return``;let{atoms:t}=vu(e,this.#Q.indices(e));return t.length?this.#e.copySvg(t):``}paste(e){let t=this.#Pe();if(!t)return!1;if(Ja(e))return this.notify(this.#a.inchiKeyRefused),!1;let n=this.#e.formatOf(e),r=$n(n,this.#Ce);return r===`read`?this.#je(n,()=>this.paste(e))&&this.#Fe(t,e,n,null,this.#a.pasteNothing,`paste`):(this.notify(r===`query-only`?this.#a.format(`pasteQueryOnly`,{format:this.#a.formatName(n)}):r===`not-yet`?this.#a.pasteNotYet(n):r===`reaction`?this.#a.reactionsRefused:this.#a.pasteNothing),!1)}#je(e,t){if(this.#e.formatLoaded(e))return!0;let n=this.#Me===null;this.#Me=t;let r=e===`inchi`?this.#a.inchiLoading:this.#a.cdxmlLoading,i=e===`inchi`?this.#a.inchiNotLoaded:this.#a.cdxmlNotLoaded;return this.notify(r,!1),this.#r.announce?.(r,`tool`),n&&this.#e.loadFormat(e).then(()=>{let e=this.#Me;if(this.#Me=null,!this.#ye){if(this.#W?.notice&&this.#W.message===r){this.#bt(null);let e=this.#xr();e&&!this.#S&&this.#Hr(e)}e?.()}},()=>{this.#Me=null,this.#ye||this.notify(i)}),!1}#Me=null;reportRead(e){let t=this.#e.lastRead;if(!t)return;let n=Wa(this.#a,e,t.structures,t.left,t.lost);n&&this.notify(n,!1)}copyAs(e,t){let n=this.#xr();if(!n||!n.atoms.length||this.#S||this.#T||this.#ye||!$e(e))return null;let r=this.#Q.empty?n.atoms.map((e,t)=>t):vu(n,this.#Q.indices(n)).atoms;if(!r.length)return null;let i=this.#e.edited?[]:bn(n);return this.#Ne(r,i,e,t)}#Ne(e,t,n,r){return n===`cdxml`?this.#e.copyCdxml(e,t,r):n===`inchi`||n===`inchikey`?this.#e.copyInchi(e,t,n,r):this.#e.copyText(e,t,n,r)}pasteAs(e,t){let n=this.#Pe();if(!n)return!1;let i=this.#e.formatOf(e);if(i===`reaction`)return this.notify(this.#a.reactionsRefused),!1;let a=On(t,i);return a&&!r(t,this.#Ce)?(this.notify(this.#a.format(`pasteQueryOnly`,{format:this.#a.formatName(t)})),!1):a?a===`inchi`&&Ja(e)?(this.notify(this.#a.inchiKeyRefused),!1):this.#je(a,()=>this.pasteAs(e,t))&&this.#Fe(n,e,a,null,this.#a.pasteNotAs(t),`paste`):(this.notify(this.#a.pasteNotYet(t)),!1)}canPasteAs(e,t){let n=On(t,this.#e.formatOf(e));return!n||!r(t,this.#Ce)?!1:this.#e.formatLoaded(n)?!!this.#e.inspect(e,n,Kt)?.atoms.length:this.#e.formatOf(e)===n}pasteAt(e,t,n,r){let i=this.#Pe();return i?this.#je(t,()=>this.pasteAt(e,t,n,r))&&this.#Fe(i,e,t,i.toMolecule(n),r,`drop`):!1}#Pe(){let e=this.#xr();return!e||!this.#e.view||this.#S||this.#T||this.#ye?null:e}#Fe(e,t,n,r,i,a){let o=this.#e.view,s=o?this.#e.readOnce(t,n,Kt,i=>this.#Le(e,o,i,t,n,r,a)):null;return s===null?(this.notify(this.#e.inspectError?.format===`reaction`?n===`cdxml`?this.#a.cdxmlReactionRefused:this.#a.reactionsRefused:i),!1):s}#Ie(e,t,n,r){let i=Oe(n.atoms);if(!i)return null;let a=Er(n.atoms,n.bonds,e.docLen),o={left:a*i.left,bottom:a*i.bottom,right:a*i.right,top:a*i.top};return{scale:a,to:r?{x:r.x-(o.left+o.right)/2,y:r.y-(o.bottom+o.top)/2}:_r(o,ir(e.atoms.map(e=>e.at),e.bonds),{x:t.cx,y:t.cy},e.docLen)}}addMolecule(e,t){let n=this.#xr(),r=this.#e.view,i=t??this.#e.formatOf(e),a=e=>Object.freeze({format:i,message:e});return!n||!r||this.#ye?a(this.#a.pasteNothing):this.#S||this.#T?a(`The sketcher is busy: a press is held or the label editor is open`):Ja(e)?a(this.#a.inchiKeyRefused):i===`reaction`?a(this.#a.reactionsRefused):this.#e.formatLoaded(i)?this.#e.readOnce(e,i,Kt,t=>{let a=this.#Ie(n,r,t,null);if(!a)return null;let o=S(`api.add`,[{op:`paste`,text:e,format:i,kekule:t.kekule,scale:a.scale,dx:a.to.x,dy:a.to.y}]);return this.#w.forget(),o?this.#st({...o,source:`api`},n):!1})===null?this.#e.inspectError?.format===`reaction`?a(this.#a.reactionsRefused):this.#e.inspectError??a(this.#a.pasteNothing):null:(this.#e.loadFormat(i).then(()=>{let n=this.#ye?null:this.addMolecule(e,t);n&&this.#e.refuseRead(n)},()=>void 0),null)}#Le(e,t,n,r,i,a,o){let s=this.#Ie(e,t,n,a);if(!s)return null;let{scale:c,to:l}=s,u=S(o===`drop`?`clipboard.drop`:`clipboard.paste`,[{op:`paste`,text:r,format:i,kekule:n.kekule,scale:c,dx:l.x,dy:l.y}]),d={atoms:new Set(e.atoms.map(e=>e.id)),bonds:new Set(e.bonds.map(e=>e.id))};return this.#w.forget(),this.#st(u,e)?(this.#Vn({atoms:[],pairs:[],before:d}),this.#Un(1,1),this.reportRead(o),!0):!1}notify(e,t=!0){if(!this.#ye){if(this.#r.noticed?.(e,`notice`),this.#r.notice){this.#r.notice(e,t);return}this.#Re(e)}}#Re(e){if(this.#ye)return;this.#q={message:null,pointer:!1},this.#r.refusal?.(null,!1),this.#bt({atoms:[],bonds:[],message:e,pointer:!1,subject:``,notice:!0});let t=this.#xr();t&&!this.#S&&this.#Hr(t)}get editingLabel(){return this.#T!==null}commitLabel(e,t=!1){let n=this.#T;if(!n)return;let r=this.#r.labelText().trim();if(r.length>80){if(!t){this.#r.refuseLabel(this.#a.aliasLimit);return}this.#zt(),e&&this.#r.focusCanvas();return}let i=r&&r!==n.opened?this.#Ut(r)??this.#ze(n.atomId,r):null;if(i){if(!t){this.#r.refuseLabel(i),this.notify(i);return}this.#zt(),this.notify(i),e&&this.#r.focusCanvas();return}this.#zt(),r&&r!==n.opened&&(this.#w.committed(r),this.#Bt(n.atomId,r)),e&&this.#r.focusCanvas()}#ze(e,t){let n=this.#xr(),r=n?.atomById(e)??-1;if(!n||r<0)return null;let i=n.labelGroupOf(r),a=i?r-i.atoms.filter(e=>e<r).length:r;return this.radicalRefusalOf(this.#Vt(n,r,t),[a])}radicalRefusalOf(e,t){let n=this.#xr(),r=e&&n&&(this.#e.edited?e:s(e,n));if(!r)return null;for(let e of t){let t=this.#e.radicalDropped(r,e);if(t!==null)return this.#a.noRoomForRadical(t)}return null}cancelLabel(){this.#T&&(this.#zt(),this.#r.focusCanvas())}labelPress(e){let t=this.#T,n=this.#w.lastLabel;if(!t||!t.byClicks||e<3||!n)return!1;this.#zt();let r=this.#ze(t.atomId,n);return r?this.notify(r):this.#Bt(t.atomId,n),this.#r.focusCanvas(),!0}clear(){if(this.#be)return;let e=Jr(this.#e.isEmpty);e&&this.#e.edit(e)&&this.#In()}relayout(e){if(this.#S||this.#e.isEmpty||this.#ye)return!1;let t=this.#xr();if(!t)return!1;let n=this.#Q.empty?void 0:t.wholeGroups(this.#Q.moved(t));if(n&&!n.length)return!1;let r=this.#e.lastError,i=this.#st(fr(e,t.docLen,n),t),a=this.#e.lastError;return!i&&n&&a&&a!==r&&this.notify(/keeps the stereo/.test(a.message)?this.#a.layoutLosesStereo:this.#a.layoutRefused),i}hydrogens(e){if(this.#S||this.#e.isEmpty||this.#ye)return!1;let t=this.#xr();return t?this.#st(wl(t,e,this.#Q.empty?null:this.#Q.moved(t),this.#Ce),t):!1}hydrogensAvailable(){let e=this.#xr();return!e||this.#ye||this.#e.isEmpty?{add:!1,remove:!1}:Cl(e,this.#Q.empty?null:this.#Q.moved(e),this.#Ce)}aromaticity(e){if(this.#S||this.#e.isEmpty||this.#ye)return!1;let t=this.#xr();return t?this.#Be(t,e,this.#Q.empty?null:this.#Q.moved(t)):!1}aromaticityState(e){let t=this.#xr();return!t||this.#ye||this.#e.isEmpty?`none`:dn(this.#e.aromaticity(e,this.#Q.empty?null:this.#Q.moved(t),this.#Ve(t)))}#Be(e,n,r){if(this.#be)return!1;let i=this.#e.aromaticity(n,r,this.#Ve(e));if(!i)return!1;let a=t(this.#a,n,i),o=Ze(n,r,i);return o?this.#st(o,e)?(a&&this.notify(a,!1),!0):!1:(a&&this.notify(a),!1)}#Ve(e){return this.#e.edited?[]:bn(e)}#He(e,t){if(t.kind===`canvas`)return null;let n=this.#$e(e,t);return[...n.atoms,...n.bonds.flatMap(t=>e.bonds[t]?[e.bonds[t].a,e.bonds[t].b]:[])]}menuAromaticity(e,t){let n=this.#xr();return!n||this.#S||n.isEmpty||this.#ye?!1:this.#Be(n,t,this.#He(n,e))}#Ue(e,t){if(!e||e.isEmpty)return;let n=t.kind===`canvas`?null:[...t.atoms,...t.bonds.flatMap(t=>e.bonds[t]?[e.bonds[t].a,e.bonds[t].b]:[])],r=this.#e.aromaticity(`aromatize`,n,this.#Ve(e)),i=this.#e.aromaticity(`dearomatize`,n,this.#Ve(e));if(r||i)return{aromatize:dn(r),dearomatize:dn(i)}}menuTarget(e){let t=this.#xr();if(!t||this.#ye)return null;let n=this.#Q.indices(t),r,i;if(e){let a=Fo(t.hitAtoms(),t.hitBonds(),e,t.bondPx,this.#o),o=!a&&this.#oe!==null&&this.#Qn(t,e,!1,!1,null)===`rotate`;r=zr(a,o||a!==null&&this.#Q.under(t,a),n),i=e}else r=_n(n,this.#Sr(t)),i=this.#We(t,r);return{target:{...r,atomIds:r.atoms.map(e=>t.atoms[e].id),bondIds:r.bonds.map(e=>t.bonds[e].id)},anchor:i}}#We(e,t){let n=[...t.atoms.map(t=>e.atoms[t].px),...t.bonds.flatMap(t=>e.bonds[t].ends.map(t=>e.atoms[t].px))];if(!n.length){let e=this.#t.size;return e?{x:e.width/2,y:e.height/2}:{x:0,y:0}}let[r,i]=[n.map(e=>e.x),n.map(e=>e.y)];return{x:(Math.min(...r)+Math.max(...r))/2,y:(Math.min(...i)+Math.max(...i))/2}}menuFacts(e){let t=this.#xr();return{target:e,rLabel:e=>!!(t&&hm(t,e)),choice:e=>t?or(t,e):null,wedged:e=>!!t?.wedgeOf(e),double:e=>t?.bonds[e]?.order===`double`,position:e=>{let n=this.#e.doubleBondPosition(e);return n===`auto`?t?.bonds[e]?.crossedShown?`centre`:t?.placement(e)??`auto`:n},empty:!t||t.atoms.length===0,stereo:this.#Ye(t),abbreviations:this.#qe(t,e),sgroups:this.#Ke(t,e),query:this.#Ge(t,e),hydrogens:t?{...Cl(t,e.kind===`canvas`?null:e.atoms,this.#Ce),...this.#Ce?{queryMode:!0}:{}}:void 0,aromaticity:this.#Ue(t,e),structures:t&&e.kind===`selection`?Ge(t,[...e.atoms,...e.bonds.flatMap(e=>t.bonds[e]?[t.bonds[e].a,t.bonds[e].b]:[])]).length:0}}#Ge(e,t){let n=e?t.atoms.filter(t=>e.atoms[t]):[],r=e?t.bonds.filter(t=>e.bonds[t]):[];return{on:this.#Ce,atomsTake:e!==null&&n.length>0&&ce(e,n),bondsTake:e!==null&&r.length>0&&be(e,r),field:t=>e?vt(e,n,t):void 0,bondType:()=>e?he(e,r,`type`):void 0,topology:()=>e?he(e,r,`topology`):void 0,off:e=>this.#m.has(e)}}menuQueryField(e,t,n){let r=this.#xr();if(!r||this.#S||!this.#Ce)return!1;let i=this.#$e(r,e).atoms,a=t===`unsaturated`?vt(r,i,`unsaturated`)!==!0:n;return this.#st(d(r,i,t,a),r)}tryQuery(e,t){let n=this.#xr(),r=e&&n&&!this.#e.edited?s(e,n):e,i=this.#e.tryQuery(r,t);return i.ok?{ok:!0,smarts:Nt(i.query)?.smarts??null}:i}menuClearQuery(e){let t=this.#xr();return!t||this.#S||!this.#Ce?!1:this.#st(it(t,this.#$e(t,e).atoms),t)}menuBondQuery(e,t){let n=this.#xr();return!n||this.#S||!this.#Ce||this.#m.has(t)?!1:this.#st(Me(n,this.#$e(n,e).bonds,t),n)}menuTopology(e,t){let n=this.#xr();return!n||this.#S||!this.#Ce?!1:this.#st(St(n,this.#$e(n,e).bonds,t),n)}#Ke(e,t){if(!e||t.kind===`canvas`||!e.sgroups.length)return{edit:!1,remove:!1};let n=this.#Pr(e,t);return n===null?{edit:!1,remove:!1}:{edit:df(e.sgroups[n].type),remove:!0}}#qe(e,t){let n={alone:!1,expand:!1,contract:!1,remove:!1,contractLabel:!1};if(!e||t.kind!==`atom`&&t.kind!==`selection`)return n;let r=this.#Je(e,t),i=ft(e,r),a=t.kind===`atom`&&t.atoms.length===1&&e.labelGroupOf(t.atoms[0])!==null,o=t.kind===`selection`&&!i.remove.length&&m(e,r)&&this.#e.contractible(r,this.#e.edited?[]:bn(e));return{alone:a,expand:i.expand.length>0,contract:i.contract.length>0,remove:i.remove.length>0,contractLabel:o}}#Je(e,t){let n=new Set(t.atoms.filter(t=>e.atoms[t]));for(let r of t.bonds)e.bonds[r]&&n.add(e.bonds[r].a).add(e.bonds[r].b);return[...n].sort((e,t)=>e-t)}menuAbbreviation(e,t){let n=this.#xr();if(!n||this.#S)return!1;let r=this.#Je(n,this.#$e(n,e));if(t===`contract-label`)return this.#st(v(r),n);let i=ft(n,r);return this.#st(t===`expand`?He(i.expand):t===`contract`?ne(i.contract):Rn(i.remove),n)}dissolveAnswered(e){let t=this.#ge;if(this.#ge=null,!t||!e||this.#ye||this.#S)return!1;let n=this.#xr();return!n||this.#e.stateToken!==t.token?!1:this.#st(Pe(t.tx,t.groups),n,[],`bond`,!1,t.keepCentre,!0)}#Ye(t){let n=this.#e,r=()=>n.stereo();return{centre:e=>{let n=t?ke(r(),e):null;return n?n.status===`specified`?`specified`:`centre`:null},specified:e=>pe(r(),e),potential:e=>{let t=e===null?null:new Set(e);return r().atoms.some(e=>t===null||t.has(e.atom))},choice:e=>t?c(n.stereo(!0),e):void 0,ez:e=>t?Vu(t,n.stereo(!0),e):null,noEz:e=>!!(t&&e.every(e=>t.bonds[e]?.order===`double`)&&n.noEz(e)),axis:e=>t?et(n.stereo(!0),e):null,removable:e=>!!(t&&ht(Or(r(),t.bonds,e))),rewedgeable:n=>!!(t&&Ct(e(r(),t.bonds,n))),labels:n.stereoShown.labels,chiralFlag:r().chiralFlag,axisGroups:n.axisGroups,axisBonds:()=>r().axes.filter(e=>e.status===`specified`).map(e=>e.bond)}}menuEz(e,t){let n=this.#xr();if(!n||this.#S)return!1;let r=this.#$e(n,e).bonds[0];if(r===void 0)return!1;let i=this.#pn();return this.#st(Hu(n,this.#e.stereo(!0),r,t,i.wedge??(e=>n.wedgeOf(e))),n)}menuStereoCentre(e,t){let n=this.#xr();if(!n||this.#S)return!1;let r=this.#$e(n,e).atoms[0];if(r===void 0)return!1;let i=this.#e.lastError,a=this.#st(Ft(this.#e.stereo(!0),r,t),n),o=this.#e.lastError;return!a&&t&&o&&o!==i&&this.notify(this.#a.format(`noCentreConfig`,{config:t,symbol:tm(n,r)})),a}menuAtropisomer(e,t){let n=this.#xr();if(!n||this.#S)return!1;let r=this.#$e(n,e).bonds[0];if(r===void 0)return!1;let i=this.#e.lastError,a=this.#st(dr(this.#e.stereo(!0),r,t),n),o=this.#e.lastError;return!a&&o&&o!==i&&this.notify(this.#a.format(`noAxisConfig`,{config:t})),a}menuStereoCommand(e,t){let n=this.#xr();if(!n||this.#S||n.isEmpty)return!1;let r=null;e.kind!==`canvas`&&(r=this.#$e(n,e));let i=this.#e.stereo(),a=t===`remove`?Dt(i,n.bonds,r):jt(i,n.bonds,r);return!a||!this.#st(a,n)?!1:((a.said?.kind===`stereo-removed`||a.said?.kind===`rewedged`)&&this.notify(Ws(this.#a,a.said),!1),!0)}menuChiralFlag(){let e=this.#xr();if(!e||this.#S||e.isEmpty)return!1;let t=this.#e.stereo();return this.#st(xe(t,!t.chiralFlag),e)}stereoTargetsOf(e){let t=this.#xr();if(!t||this.#ye)return null;let n=this.#e.stereo(),r,i;if(e){let n=e.kind===`canvas`?null:this.#$e(t,e);[r,i]=[n?.atoms??null,n?.bonds??null]}else[r,i]=this.#Q.empty?[null,null]:[this.#Q.indices(t).atoms,this.#Q.indices(t).bonds];let a=tr(n,r),o=this.#e.axisGroups?n.axes.filter(e=>e.status===`specified`&&(i===null||i.includes(e.bond))).map(e=>e.bond):[];return{atomIds:a.map(e=>t.atoms[e].id),bondIds:o.flatMap(e=>t.bonds[e]?[t.bonds[e].id]:[]),opening:le(n,a)}}applyStereoGroup(e,t){let n=this.#xr();if(!n||this.#S)return!1;let r=e.map(e=>n.atomById(e)).filter(e=>e>=0);return this.#st(Xe(this.#e.stereo(),r,t),n)}stereoRows(){let e=this.#xr();if(!e||this.#ye)return null;let t=this.#e.stereo(!0);return{info:t,rows:ad(t,{elements:e.atoms.map(e=>e.element),bonds:e.bonds.map(e=>({a:e.a,b:e.b})),ez:n=>Vu(e,t,n),axisGroups:this.#e.axisGroups},this.#a)}}selectItems(e){if(this.#S||this.#ye)return;!Fr(this.#s)&&this.#s!==`hand`&&this.setTool(this.#$);let t=this.#xr();if(!t)return;let n=e.atoms.filter(e=>t.atoms[e]),r=e.bonds.filter(e=>t.bonds[e]);!n.length&&!r.length?this.#In():this.#Q.replace(t,{atoms:n,bonds:r}),this.#fr(t)}markItems(e){let t=this.#xr();!t||this.#ye||this.#S||(e&&(e.atoms.length||e.bonds.length)?this.#br(t,e):(this.#Vr(),this.#Wr=null,this.#Hr(t)))}panelConfig(e,t){let n=this.#xr();if(!n||this.#S||this.#ye)return!1;let r=this.#e.lastError;if(e.kind===`centre`){if(!n.atoms[e.atom]||t!==`r`&&t!==`s`&&t!==`unspecified`)return!1;let i=t===`r`?`R`:t===`s`?`S`:null,a=this.#st(Ft(this.#e.stereo(!0),e.atom,i),n),o=this.#e.lastError;return!a&&i&&o&&o!==r&&this.notify(this.#a.format(`noCentreConfig`,{config:i,symbol:tm(n,e.atom)})),a}if(!n.bonds[e.bond])return!1;if(e.kind===`double`){if(t!==`e`&&t!==`z`&&t!==`unknown`)return!1;let r=this.#pn();return this.#st(Hu(n,this.#e.stereo(!0),e.bond,t===`unknown`?`unknown`:t===`e`?`E`:`Z`,r.wedge??(e=>n.wedgeOf(e))),n)}if(t!==`p`&&t!==`m`)return!1;let i=t===`p`?`P`:`M`,a=this.#st(dr(this.#e.stereo(!0),e.bond,i),n),o=this.#e.lastError;return!a&&o&&o!==r&&this.notify(this.#a.format(`noAxisConfig`,{config:i})),a}panelGroup(e,t,n){let r=this.#xr();if(!r||this.#S||this.#ye)return!1;let i=e.filter(e=>r.atoms[e]),a=this.#e.axisGroups?t.filter(e=>r.bonds[e]):[],o=n.type===`abs`||n.type===`none`?{type:n.type}:{type:n.type,group:n.id??`new`};return this.#st(Xe(this.#e.stereo(),i,o,a),r)}panelStereoCommand(e){let t=this.#xr();if(!t||this.#S||t.isEmpty||this.#ye)return!1;let n=this.#Q.empty?null:this.#Q.indices(t),r=this.#e.stereo(),i=e===`remove`?Dt(r,t.bonds,n):jt(r,t.bonds,n);return!i||!this.#st(i,t)?!1:((i.said?.kind===`stereo-removed`||i.said?.kind===`rewedged`)&&this.notify(Ws(this.#a,i.said),!1),!0)}panelCommandsAct(){let t=this.#xr();if(!t||t.isEmpty)return{remove:!1,rewedge:!1};let n=this.#Q.empty?null:this.#Q.indices(t),r=this.#e.stereo();return{remove:ht(Or(r,t.bonds,n)),rewedge:Ct(e(r,t.bonds,n))}}atomPx(e){let t=this.#xr();return t?(e??t.atoms.map(e=>e.index)).flatMap(e=>t.atoms[e]&&!t.atoms[e].hidden?[t.atoms[e].px]:[]):[]}#Xe=``;#Ze(e){let t=e&&!this.#he&&this.#ot()&&!this.#S?this.#Sr(e):null,n=t?`${t.kind} ${t.index}`:``;n!==this.#Xe&&(this.#Xe=n,this.#r.hotChanged?.(t))}menuOpened(e){if(this.#ye)return;this.#Tn(),this.#Br(),this.#Dr(),this.#he=e,this.#Lr(),this.#y=null,this.#b=null,this.#St(),this.#Qt(),this.#Yt();let t=this.#xr(),n=t?this.#Qe(t):null;t&&n&&(this.#x={kind:n.kind,id:n.kind===`atom`?t.atoms[n.index].id:t.bonds[n.index].id}),t&&this.#Hr(t)}menuClosed(){if(!this.#he)return;this.#he=null;let e=this.#xr();e&&!this.#ye&&this.#Hr(e)}#Qe(e){let t=this.#he;if(!t||t.kind!==`atom`&&t.kind!==`bond`)return null;let n=t.kind===`atom`?e.atomById(t.atomIds[0]??-1):e.bondById(t.bondIds[0]??-1);return n>=0?{kind:t.kind,index:n}:null}#$e(e,t){return{atoms:t.atomIds.map(t=>e.atomById(t)).filter(e=>e>=0),bonds:t.bondIds.map(t=>e.bondById(t)).filter(e=>e>=0)}}#et(e,t){let n=this.#$e(e,t);return[...n.atoms,...n.bonds.flatMap(t=>e.bonds[t]?[e.bonds[t].a,e.bonds[t].b]:[])]}menuAlign(e,t){let n=this.#xr();return!n||this.#S?!1:this.#st(Ne(n,this.#et(n,e),t),n)}menuDistribute(e,t){let n=this.#xr();return!n||this.#S?!1:this.#st(Le(n,this.#et(n,e),t),n)}menuCharge(e,t){let n=this.#xr();if(!n||this.#S)return!1;let{tx:r,past:i}=$t(n,this.#$e(n,e).atoms,t);return r?this.#st(r,n):(i>=0&&this.notify(this.#a.chargeLimit(n.info(i)?.symbol||`?`,t)),!1)}menuHydrogens(e,t){let n=this.#xr();return!n||this.#S?!1:this.#st(wl(n,t,e.kind===`canvas`?null:this.#$e(n,e).atoms,this.#Ce),n)}menuBondChoice(e,t){let n=this.#xr();return!n||this.#S?!1:this.#st(Cr(n,this.#$e(n,e).bonds,t,this.#at,this.#Ce),n)}applyBondChoice(e,t){let n=this.#xr();if(!n||this.#S)return!1;let r=e.map(e=>n.bondById(e)).filter(e=>e>=0);return this.#st(Cr(n,r,t,this.#at,this.#Ce),n)}bondFormOf(e){let t=this.#xr(),n=t?e.map(e=>t.bondById(e)).filter(e=>e>=0):[];return!t||!n.length?null:{form:An(t,n),shows:this.#Ce||n.some(e=>Jt(t,e))}}applyBondForm(e,t,n){let r=this.#xr();if(!r||this.#S)return{ok:!0};let i=e.map(e=>r.bondById(e)).filter(e=>e>=0),a=rn(r,i,t,n,this.#at,this.#Ce);if(a&&n.custom===!0){let e=this.#e.tryQuery(this.#e.edited?a:s(a,r));if(!e.ok)return e}return this.#st(a,r),{ok:!0}}bondSmarts(e,t,n){let r=this.#xr(),i=r?e.map(e=>r.bondById(e)).filter(e=>e>=0):[];if(!r||!i.length)return``;let a=rn(r,i,t,{...n,custom:!1,filled:void 0},this.#at,this.#Ce),o=this.#e.tryQuery(a&&!this.#e.edited?s(a,r):a,null,i[0]);return o.ok?kt(o.query)?.smarts??``:``}bondChoiceOf(e){let t=this.#xr(),n=t?e.map(e=>t.bondById(e)).filter(e=>e>=0):[];if(!t||!n.length)return null;let r=n.map(e=>or(t,e));return{shared:r.every(e=>e===r[0])?r[0]:null,count:n.length}}menuReverse(e){let t=this.#xr();return!t||this.#S?!1:this.#st(Nr(t,this.#$e(t,e).bonds),t)}menuPosition(e,t){let n=this.#xr();if(!n||this.#S)return!1;let r=this.#$e(n,e).bonds;return this.#tt(n,r,()=>this.#st(Gn(n,r,t,this.#at),n))}#tt(e,t,n){let r=this.#it,i=n();return t.some(t=>e.bonds[t]?.order===`double`&&e.bonds[t]?.crossedShown)?(this.#rt(this.#a.crossedPositionKept),this.#nt=i&&this.#it===r?this.#a.crossedPositionKept:null,i):i}#nt=null;#rt(e){if(this.#ye)return;this.#q={message:null,pointer:!1},this.#r.refusal?.(null,!1);let t={atoms:[],bonds:[],message:e,pointer:!1,subject:``,notice:!0};this.#bt(t),this.#ut(t,Zp);let n=this.#xr();n&&!this.#S&&this.#Hr(n)}#it=0;menuDelete(e){let t=this.#xr();if(!t||this.#S)return!1;if(e.fromSelection)return this.deleteSelection();let n=this.#$e(t,e);return this.#st(mt(t,n.atoms,n.bonds),t)}menuCopy(e,t,n){let r=this.#xr();if(!r||this.#S||this.#T||this.#ye||t!==void 0&&!$e(t))return null;let{atoms:i}=vu(r,this.#$e(r,e));if(!i.length)return null;let a=this.#e.edited?[]:bn(r);return t?this.#Ne(i,a,t,n):this.#e.copyText(i,a)}menuSelectStructure(e){if(this.#S||this.#ye)return;Fr(this.#s)||this.setTool(this.#$);let t=this.#xr();t&&(this.#Q.replace(t,Kp(t,this.#$e(t,e))),this.#fr(t))}menuAtomIds(e){let t=this.#xr();return t?e.atomIds.filter(e=>{let n=t.atomById(e);return n>=0&&!hm(t,n)&&!t.contractedGroupOf(n)}):[]}menuEditLabel(e){let t=this.#xr(),n=t?t.atomById(e.atomIds[0]??-1):-1;if(!(!t||n<0||this.#S||this.#T)){if(hm(t,n)){this.#r.openRGroup?.({atomId:t.atoms[n].id});return}this.#Rt(t,n,!1)}}apply(e){let t=this.#xr();return t?this.#st(e,t):!1}atomIndex(e){return this.#xr()?.atomById(e)??-1}bondIndex(e){return this.#xr()?.bondById(e)??-1}destroy(){this.#ye=!0,this.#or(),this.#Qt(),this.#S=null,this.#_=null,this.#Wr=null,this.#C=[],this.#L=[],this.#x=null,this.#b=null,this.#T=null,this.#E=null,this.#N=null,this.#P=null,this.#F=null,this.#I=null,this.#R=null,this.#W=null,this.#G?.cancel(),this.#G=null,this.#K=null,this.#q={message:null,pointer:!1},this.#Z=null,this.#Y=null,this.#X=[],this.#Q.clear(),this.#ee=[],this.#te=null,this.#ne=null,this.#oe=null,this.#ue=null,this.#de=null,this.#le=null,this.#he=null,this.#ge=null,this.#_e?.cancel?.(),this.#_e=null,this.#ve?.cancel?.(),this.#ve=null}#at=e=>this.#e.doubleBondPosition(e);#ot(){return om(this.#x,this.#b)}#st(e,t,n=[],r=`bond`,i=!1,a=!1,o=!1){if(this.#be&&e?.source!==`api`)return!1;if(!e){let e=this.#pt(t,n,r,i);return e&&this.#lt(e,t),!1}let c=e,l=o?[]:Vr(t,e);if(l.length){if(this.#r.confirmDissolve)return this.#ge={tx:e,groups:l,token:this.#e.stateToken,keepCentre:a},this.#r.confirmDissolve(l.map(e=>e.label??``)),!1;c=Pe(e,l)}let u=this.#e.edited?c:s(c,t);if(!u)return!1;let d=this.#Tt(t),f=u.ops.some(e=>e.op===`removeAtom`||e.op===`removeBond`)?this.#Sr(t):null,p=f?{hot:f,graph:_m(t),ids:t.atoms.map(e=>e.id)}:null;if(!this.#e.edit(u)){let e=this.#e.noEzRefused;if(e!==null&&t.bonds[e]){let n=this.#ht(t,e,i);this.#ye||this.#r.noticed?.(n.message,`refusal`),this.#bt(n),this.#ut(n),this.#xr()===t&&!this.#S&&this.#Hr(t)}return!1}if(!a&&this.#ue){this.#ue=null;let e=this.#xr();e&&this.#Xn(e)}return p&&this.#ct(p),this.#Dt(d),!0}#ct(e){let t=this.#xr();if(!t||this.#Sr(t))return;let n=ac(e.graph,e.hot,e.ids.map(e=>t.atomById(e)));n<0||(this.#x={kind:`atom`,id:t.atoms[n].id},this.#N=this.#y?{...this.#y}:null,this.#S||this.#Hr(t))}#lt(e,t){this.#ye||this.#r.noticed?.(e.message,`refusal`);let n=this.#W?.pointer&&this.#W.subject===e.subject&&!e.pointer?{...e,pointer:!0}:e;this.#bt(n),this.#ut(n,Qp,!n.pointer),this.#xr()===t&&!this.#S&&this.#Hr(t)}#ut(e,t=Qp,n=!0){let r=this.#z.after(t,()=>{if(this.#ye||this.#G?.refusal!==e||(this.#G=null,this.#W!==e))return;this.#K=null,this.#bt(null);let t=this.#xr();t&&!this.#S&&this.#Hr(t)});this.#G={refusal:e,cancel:r,held:n}}#dt(){return this.#W!==null&&this.#G?.refusal===this.#W&&this.#G.held}#ft(e,t,n,r,i){return t?null:this.#pt(e,n,r,i)}#pt(e,t,n,r){let i=this.#c;return n===`template`&&this.#u===`template`&&i?.unplaceable?{atoms:[],bonds:[],message:this.#mt(i),pointer:r,subject:this.#yt(e)}:this.#vt(e,t,n,r)}#mt(e){return e.cannotMake?.includes(`radical`)?this.#a.templateRadicals:this.#a.templateBondType}#ht(e,t,n){return{atoms:[],bonds:[t],message:this.#a.noEzLinearCentre,pointer:n,subject:this.#yt(e)}}#gt(e,t,n){let{plan:r,noRoom:i}=E(()=>ie(t,e,n.from,n.to,this.#at));return this.#ft(e,r,i,`bond`,!0)}#_t(e,t){if(this.#K?.key===e)return this.#K.refusal;let n=t();return this.#K={key:e,refusal:n},n}#vt(e,t,n,r){if(!t.length)return null;let i=t[0];return{atoms:t.filter(t=>t.atom>=0&&e.atoms[t.atom]).map(e=>e.atom),bonds:[],message:this.#a.noRoom(em(this.#a,i.symbol,i.atom),n),pointer:r,subject:this.#yt(e)}}#yt(e){let t=this.#b;return t?`${t.kind}:${t.id}`:``}#bt(e,t=!0){this.#G&&this.#G.refusal!==e&&(this.#G.cancel(),this.#G=null),this.#W=e,t&&this.#xt()}#xt(){if(this.#ye)return;let e=this.#W?.message??this.#Z,t=!!this.#W?.pointer;(e!==this.#q.message||t!==this.#q.pointer)&&(this.#q={message:e,pointer:t},this.#r.refusal?.(e,t))}#St(e=!0){this.#K=null,this.#bt(null,e)}get problemMarks(){return this.#J}set problemMarks(e){if(e===this.#J||this.#ye)return;this.#J=e;let t=this.#xr();t&&this.#wt(t)}#Ct(e){return this.#Y?.model!==e&&(this.#Y={model:e,byAtom:yp(this.#e.problems)}),this.#Y.byAtom}#wt(e){for(let e of this.#X)e.removeAttribute(`data-problem`);this.#X=[];let t=this.#Ct(e),n=[...t.keys()].filter(t=>e.atoms[t]).sort((e,t)=>e-t);if(n.length){let e=this.#vr();for(let r of n){let n=e.atoms.get(r);n&&(n.setAttribute(`data-problem`,t.get(r).kind),this.#X.push(n))}}!n.length||!this.#J?this.#n.hideProblems():this.#n.showProblems(this.#Kr(e,n,[]))}#Tt(e){return new Set(this.#e.problems.map(t=>bp(t,t=>e.atoms[t]?.id??-1)))}#Et(e,t){return t?t.check?fp(t.check,this.#a):sp(t.kind,e.info(t.atoms[0]??-1),this.#a):null}#Dt(e){let t=this.#xr();if(!t)return;let n=this.#e.problems.find(n=>!e.has(bp(n,e=>t.atoms[e]?.id??-1)));n&&this.#Re(this.#Et(t,n)??n.message)}#Ot(e,t,n){let r=n&&t?.kind===`atom`&&!this.#S&&!this.#E?this.#Et(e,this.#Ct(e).get(t.index)):null;r!==this.#Z&&(this.#Z=r,this.#xt())}#kt(e){if(this.#S||this.#E||this.#T)return;let t=this.#yt(e),n=this.#xe();if(this.#W&&this.#W.subject!==t&&!this.#dt()&&this.#bt(null),this.#sn())return;let r=this.#b,i=r?r.kind===`atom`?e.atomById(r.id):e.bondById(r.id):-1;if(!r||i<0||!n||!(n.kind===`bond`||n.kind===`element`||n.kind===`chain`&&r.kind===`bond`)){this.#W?.pointer&&this.#bt(null);return}let a=this.#_t(`hover:${this.#s}:${t}`,()=>{if(r.kind===`bond`&&n.kind===`bond`&&(n.stereo===`crossed`||n.stereo===`wavy`)&&e.bonds[i]?.order===`double`&&this.#e.noEz([i]))return this.#ht(e,i,!0);let{plan:t,noRoom:a}=E(()=>r.kind===`atom`?lt(n,e,i,this.#at):g(n,e,i,this.#at));return this.#ft(e,t,a,`bond`,!0)}),o=this.#W;!(a&&o&&this.#dt()&&o.message===a.message&&`${o.atoms}|${o.bonds}`==`${a.atoms}|${a.bonds}`)&&a!==this.#W&&(a||this.#W?.pointer)&&this.#bt(a)}#At(e,t){let n=this.#w.double(t),r=n!==null&&this.#e.takeBack(n)?this.#xr():e,i=r?.atomById(t)??-1;if(!(!r||i<0)){if(hm(r,i)&&this.#r.openRGroup){this.#Tn(),this.#Yt(),this.#r.openRGroup({atomId:t});return}this.#Rt(r,i,!0)}}#jt(e,t,n){if(n.hit?.kind===`bond`)return;let r=n.hit?.kind===`atom`?t.atoms[n.hit.index].id:null;if(e===`attachment`){r!==null&&this.#r.openAttachmentPoints?.(r);return}this.#r.openRGroup?.(r===null?{at:t.toMolecule(n.start)}:{atomId:r})}#Mt(e,t,n){if(t.hit?.kind!==`atom`)return null;let r=this.#Jr(e,t,n);if(!r||!(`atom`in r.from))return null;let i=r.from.atom;if(!(`atom`in r.to))return Fp(i,{point:r.to.point},this.#Nt(e,t,-1));let a=r.to.atom;return Fp(i,{atom:a,bonded:e.bondBetween(i,a)>=0,was:this.#Ft(e,a)},this.#Nt(e,t,a))}#Nt(e,t,n){return t.rLabels??=this.#Lt(e),this.#It(e,n,t.rLabels)}#Pt(e,t){this.#st(Pp(t,this.#Ft(e,t),this.#It(e,t)),e)}#Ft(e,t){return hm(e,t)?this.#e.atomRNumbers(t)??[]:[]}#It(e,t,n=this.#Lt(e)){let r=[];for(let[e,i]of n)e!==t&&r.push(...i);return r}#Lt(e){let t=new Map;return e.atoms.forEach((e,n)=>{e.element===`R#`&&t.set(n,this.#e.atomRNumbers(n)??[])}),t}#Rt(e,t,n){let r=e.atoms[t],i=e.labelGroupOf(t)?.label??r.label??(r.element&&r.element!==`C`?r.element:``);this.#Tn(),this.#Yt(),this.#Br(),this.#T={atomId:r.id,opened:i,byClicks:n},this.#Lr(),this.#r.showLabelEditor(r.px,i,gm(e.bondPx)),this.#Xn(e)}#zt(){if(!this.#T)return;this.#T=null,this.#r.hideLabelEditor();let e=this.#xr();e&&this.#Xn(e)}#Bt(e,t){let n=this.#xr(),r=n?.atomById(e)??-1;!n||r<0||this.#st(this.#Vt(n,r,t),n)}#Vt(e,t,n){let r=this.#Ht(n);return r?ri(e,t,r):ho(e,t,n)}#Ht(e){if(!this.#Ce)return null;let t=Lt(e,e=>Ve(e)!==void 0);return t?ln(t):null}#Ut(e){if(!this.#Ce)return null;let t=Lt(e,e=>Ve(e)!==void 0);if(!t)return null;if(t.kind===`refused`)return t.why===`unknown-element`?this.#a.notAnElement(t.element??e):t.why===`too-many`?this.#a.listTooLong:this.#a.listEmpty;let n=ln(t);return n?this.#Te(Ce(n),e.trim()):null}#Wt(e){return{symbol:e,colour:Ln(e,this.#t.colours)}}#Gt(){let e=this.element;if(e)return this.#Wt(e);let t=this.#s===`element.query`?this.#p:null;if(t)return{symbol:n(t),colour:`#000000`};let r=this.#cn();return r?{symbol:r.kind===`group`?r.label:r.salt.abbreviation,colour:`#000000`}:null}#Kt(){if(this.#S||this.#E)return;if(this.#sn()){this.#$t();return}let e=this.#y&&!this.#T?this.#xr():null,t=this.#Gt();this.#F=null,t&&this.#y&&e?this.#qt(this.#y,t,e.bondPx):this.#Yt()}#qt(e,t,n){this.#R={at:{...e},label:t},this.#n.showLabel(e,t,n)}#Jt(e){let t=this.#T;if(t){let n=e.atomById(t.atomId);n>=0&&this.#r.moveLabelEditor(e.atoms[n].px,gm(e.bondPx))}let n=this.#E;if(n){this.dragElement(n.symbol,n.at);return}let r=this.#R;r&&!this.#S&&(this.#R={at:r.at,label:this.#Wt(r.label.symbol)},this.#n.showLabel(r.at,this.#R.label,e.bondPx))}#Yt(){this.#I=null,this.#se=null,this.#Xt()}#Xt(){this.#F=null,this.#se&&(this.#se.shown=!1),this.#R=null,this.#Zt([]),this.#n.hidePreview(),this.#n.hideCount()}#Zt(e){for(let t of this.#L)e.includes(t)||t.removeAttribute(`visibility`);for(let t of e)t.setAttribute(`visibility`,`hidden`);this.#L=[...e]}#Qt(){let e=[this.#B,this.#V,this.#H];this.#B=this.#V=this.#H=null;for(let t of e)t?.()}#$t(){this.#ye||this.#B||this.#H||(this.#B=this.#z.frame(()=>{this.#B=null,this.#tn()}))}#en(){this.#ye||(this.#H?.(),this.#H=this.#z.painted(()=>{this.#H=null,this.#tn()}))}#tn(){if(this.#ye||this.#E||this.#T||!this.#y)return;let e=this.#xr(),t=this.#xe();if(!e)return;let n=this.#S,r=this.#sn();if(n){if(n.cancelled||!n.dragging)return;let i=this.#y;if(n.grip===`rotate`)this.#tr(e,n,i);else if(n.grip===`centre`)this.#nr(e,i);else if(cm(t,n))this.#mn(e,n,i);else if(r){let t=E(()=>this.#dn(r,e,n,i));this.#in(t.plan,e,{noRoom:t.noRoom,why:this.#rn()})}else if(t?.kind===`chain`){let{plan:t,noRoom:r}=E(()=>this.#yn(e,n,i));t?.tx&&this.#in(t.tx,e,{why:`chain`})?this.#on(e,t,n.touch,i):(this.#Yt(),t?.tx||this.#_n(e,this.#vt(e,r,`chain`,!0)))}return}r&&!this.#j&&this.#nn(r,e,this.#y)}#nn(e,t,n){let r=Fo(t.hitAtoms(),t.hitBonds(),n,t.bondPx,this.#o);if(!r){this.#V?.(),this.#V=null,this.#I=``;let r=!!this.#k.shift,i=E(()=>e.onCanvas(t,t.toMolecule(n),0,r));this.#in(i.plan,t,{noRoom:i.noRoom,why:this.#rn()});return}let i=r.kind===`atom`?`atom:${t.atoms[r.index].id}`:`bond:${t.bonds[r.index].id}`,a=150-(this.#D()-this.#U);if(this.#V?.(),this.#V=null,i!==this.#I&&a>0){this.#Yt(),this.#W&&this.#W.subject!==this.#yt(t)&&(this.#bt(null),this.#Hr(t)),this.#V=this.#z.after(a,()=>{this.#V=null,this.#$t()});return}this.#I=i;let o=!!this.#k.alt,s=E(()=>r.kind===`atom`?e.onAtom(t,r.index,o):e.onBond(t,r.index));this.#in(s.plan,t,{noRoom:s.noRoom,why:this.#rn()})}#rn(){return this.#xe()?.kind===`template`?`template`:`ring`}#in(e,t,n={}){let r=n.opacity??.6,i=n.judge??!0;if(this.#se&&(this.#se.shown=!1),!e)return this.#Xt(),i&&this.#_n(t,this.#pt(t,n.noRoom??[],n.why??`ring`,!0)),!1;let a=bt(e),o=a||this.#e.edited?e:s(e,t);if(!o)return this.#Xt(),!1;let c=`${a?`alone `:``}${r} ${JSON.stringify(o.ops)}`;if(c===this.#F)return!0;i&&this.#_n(t,null);let l=a?void 0:this.#an(o);if(l===null)return this.#Xt(),!1;let u,d;if(l)({shown:u,frame:d}=l);else{let e=this.#e.preview(o,a),t=e?ar(e):null;if(!t)return this.#Xt(),!1;u=a?{parts:zc(t),replaced:[]}:Rc(t,this.#t.svg),d={viewBox:t.getAttribute(`viewBox`),preserveAspectRatio:t.getAttribute(`preserveAspectRatio`)}}return this.#R=null,this.#n.showDrawn(u.parts,d,r),this.#Zt(u.replaced),this.#F=c,!0}#an(e){let t=this.#t.svg,n=this.#t.text;if(!t||!n)return;let r=this.#e.previewText(e);if(r==null)return r;let i=t.ownerDocument,a=Lc(r,t,n,e=>qt(i,e));return a?{shown:a,frame:{viewBox:Mc(r,`viewBox`),preserveAspectRatio:Mc(r,`preserveAspectRatio`)}}:void 0}#on(e,t,n,r){let i=String(t.count),a=pa({origin:e.toCanvas(t.origin),end:e.toCanvas(t.end),direction:t.direction,flip:!!this.#k.ctrl,bondPx:e.bondPx,pointer:r,touch:n,area:this.#xn()??{width:1/0,height:1/0},label:sa(i,e.bondPx)});this.#n.showCount(a,i,e.bondPx)}#sn(){let e=this.#xe();return e?.kind===`ring`?Kd(e.ring):e?.kind===`template`&&this.#u===`template`?this.#l:null}#cn(){return this.#s===`template.library`?this.#u===`group`&&this.#d?{kind:`group`,label:this.#d}:this.#u===`salt`&&this.#f?{kind:`salt`,salt:this.#f}:null:null}#ln(e){return this.#h?.index!==e.index&&(this.#h={index:e.index,part:this.#e.inspect(e.molblock,`molV2000`,Kt)}),this.#h.part}#un(e,t,n,r){let i=t.hit?{kind:t.hit.kind,index:t.hit.index}:null;if(r.kind===`salt`){let a=this.#ln(r.salt);return a?t.dragging?Eo(e,e.toMolecule(n),r.salt,a):To(e,i,e.toMolecule(t.start),r.salt,a):null}if(!t.dragging)return _o(e,i,e.toMolecule(t.start),r.label);if(t.hit?.kind===`atom`){let i=this.#Jr(e,t,n);return i&&`atom`in i.from?vo(e,i.from.atom,i.to,r.label):null}return t.hit?null:yo(e,this.#wn(e,n),e.toMolecule(n),r.label)}#dn(e,t,n,r){let i=!!this.#k.alt;if(n.hit?.kind===`atom`)return e.onAtom(t,n.hit.index,i);if(n.hit?.kind===`bond`){let i=t.bonds[n.hit.index],a=n.dragging?Bl(t.atoms[i.a].at,t.atoms[i.b].at,t.toMolecule(r)):null;return e.onBond(t,n.hit.index,a??void 0)}let a=n.dragging?Ir(im(n.start,r)):0;return e.onCanvas(t,t.toMolecule(n.start),a,!!this.#k.shift)}#fn(e,t,n){if(!t.hit||!t.moving)return{plan:au,delta:{x:0,y:0}};let r=e.toMolecule(n),i=Ql({x:r.x-t.startMol.x,y:r.y-t.startMol.y},!!this.#k.shift);if(this.#M={kind:`move`,length:Math.hypot(i.x,i.y),angle:Math.atan2(i.y,i.x)*180/Math.PI},t.copy)return{plan:bu(e,this.#Q.indices(e),i,this.#pn()),delta:i};let a=this.#pn().wedge;return t.whole?{plan:du(e,t.moving,i,`selection`,a),delta:i}:{plan:fu(e,t.hit,i,a),delta:i}}#pn(){return{doubleBondPosition:e=>this.#e.doubleBondPosition(e)}}#mn(e,t,n){let{plan:r,delta:i}=this.#fn(e,t,n),a=r.refused&&r.blocker>=0?e.atoms[r.blocker]:void 0,o=a?em(this.#a,tm(e,a.index),a.index):``;if(this.#bt(a?{atoms:[a.index],bonds:[],message:r.crowded?this.#a.tooClose(o):this.#a.noRoom(o,`merge`),pointer:!0,subject:``}:null),r.targets.atoms.length||r.targets.bonds.length?this.#br(e,r.targets):a?this.#Hr(e,{kind:`atom`,index:a.index}):this.#Vr(),t.copy&&t.whole&&r.tx&&!r.merged.size&&this.#gn(e,t,r.tx,i))return;let s=r.refused?this.#hn(e,t,i):null;if(!(t.copy&&t.whole&&s&&this.#gn(e,t,s,i,.4))){if(this.#in(s??r.tx,e,{opacity:s?$p:1,judge:!1}),t.whole){let n=e.toCanvas(t.startMol),a=e.toCanvas({x:t.startMol.x+i.x,y:t.startMol.y+i.y});this.#fr(e,dm,r.tx||s?{x:a.x-n.x,y:a.y-n.y}:null)}else t.moving&&!this.#Q.empty&&this.#fr(e,new Set(t.moving))}}#hn(e,t,n){if(t.copy&&t.whole){let t=vu(e,this.#Q.indices(e)).atoms;return t.length?S(`selection.copy`,[{op:`duplicate`,atoms:t,dx:n.x,dy:n.y}]):null}return t.moving?.length?pu(e,t.moving,n):null}#gn(e,t,n,r,i=1){let a=n.ops[0];if(a?.op!==`duplicate`)return!1;let o=this.#se;if(!o||o.gesture!==t||o.model!==e){let r=this.#e.edited?n:s(n,e),c=r&&this.#e.partDrawing(a.atoms,r.ops.slice(0,r.prepared??0));if(!c)return!1;let l=ar(c);o={gesture:t,model:e,parts:l?zc(l):[],frame:{viewBox:l?.getAttribute(`viewBox`)??null,preserveAspectRatio:l?.getAttribute(`preserveAspectRatio`)??null},shown:!1,opacity:i,marks:null},this.#se=o}let c=e.toCanvas(t.startMol),l=e.toCanvas({x:t.startMol.x+r.x,y:t.startMol.y+r.y}),u={x:l.x-c.x,y:l.y-c.y};return(!o.shown||o.opacity!==i)&&(this.#n.showDrawn(o.parts,o.frame,i),o.opacity=i,this.#Zt([]),this.#F=null,this.#R=null,o.shown=!0),this.#n.moveDrawn(u.x,u.y),o.marks&&o.marks.version===this.#ce?this.#n.moveSelection(u.x-o.marks.at.x,u.y-o.marks.at.y):(this.#fr(e,dm,u),o.marks={at:u,version:this.#ce}),!0}#_n(e,t){if(!this.#W&&!t){this.#xt();return}this.#bt(t),(!this.#S||this.#sn()||this.#xe()?.kind===`chain`)&&this.#Hr(e)}addChain(e,t){let n=this.#xr();if(!n||this.#S||this.#ye)return!1;let r=`atomId`in e?n.atomById(e.atomId):-1;if(`atomId`in e&&r<0)return!1;this.#vn=t;let i=Jd(n,r>=0?{atom:r}:{point:e.point},t);return this.#st(i?.tx??null,n,!i&&r>=0?[{atom:r,symbol:tm(n,r)}]:[],`chain`)}#vn=null;#yn(e,t,n){return t.hit?.kind===`bond`?null:qd(e,t.hit?.kind===`atom`?{atom:t.hit.index}:{point:e.toMolecule(t.start)},e.toMolecule(n),!!this.#k.alt,!!this.#k.ctrl,this.#xn())}#bn(e,t,n){let{plan:r,noRoom:i}=E(()=>nf(e,t,n,this.#at));if(!r){this.#st(null,t,i,e===`6`||e===`7`||e===`3`||e===`a`?`ring`:`bond`);return}let a=r.hotspot;this.#P=a===null?null:typeof a==`number`?a:t.atoms.length+a.added,this.#st(r.tx,t)||(this.#P=null)}#xn(){let e=this.#t.svg,t=Number.parseFloat(e?.getAttribute(`width`)??``)||0,n=Number.parseFloat(e?.getAttribute(`height`)??``)||0;return t>0&&n>0?{width:t,height:n}:void 0}#Sn(e){let t=this.#t.svg,n=Number.parseFloat(t?.getAttribute(`width`)??``)||0,r=Number.parseFloat(t?.getAttribute(`height`)??``)||0;return e.x>=0&&e.y>=0&&e.x<=n&&e.y<=r}#Cn(e){let t=this.#N;return t?Math.hypot(e.x-t.x,e.y-t.y)<1||(this.#N=null,!1):!1}#wn(e,t){let n=Fo(e.hitAtoms(),[],t,e.bondPx,{atomsOnly:!0});return n?.kind===`atom`?n.index:-1}#Tn(){if(!this.#S)return;let e=this.#S.dragging;if(this.#S=null,this.#de=null,this.#or(),this.#Yt(),this.#n.hideArea(),this.#ne=null,this.#te||e){this.#te=null;let e=this.#xr();e&&this.#fr(e)}}#En(e,t){let n=t.kind===`atom`?e.atoms[t.index]:e.bonds[t.index];n&&(this.#x={kind:t.kind,id:n.id},this.#N=this.#y?{...this.#y}:null,this.#O=!0,this.#St(),this.#Hr(e),this.#Dn(e))}#Dn(e){let t=this.#Sr(e);if(!(!t||this.#ye)){if(t.kind===`atom`){let n=e.atoms[t.index],r=e.info(t.index),i=r&&!e.labelGroupOf(t.index)?ee(r.query??null,r.symbol,this.#a):null,a=Ps(this.#a,i?i.join(`; `):n.label||tm(e,t.index),r?.charge??0);this.#r.announce?.(Fs(this.#a,a,t.index,e.bondsOf(t.index).length),`hotspot`)}else{let n=e.bonds[t.index],r=e.wedgeOf(t.index);this.#r.announce?.(Is(this.#a,t.index,Ms(n.order,r?r.stereo:n.stereo),n.ends[0],n.ends[1]),`hotspot`)}}}#On(e,t){if(this.#Sr(e)||!this.#Q.empty||this.#T)return!1;if(e.atoms.length){if(!t)return!1;let n=ic({atoms:e.atoms.filter(e=>!e.hidden).map(e=>e.index),bonds:[]},_m(e),null,t=>xm(e,t));return n<0?!1:(this.#En(e,{kind:`atom`,index:n}),!0)}let n=this.#xn();if(!n)return!1;let r=e.toMolecule({x:n.width/2,y:n.height/2}),i=this.#An(e,r);if(!i){this.#De(Gt);let t=this.#An(e,r);i=gt(t,{kind:`tool`,tool:D(Gt,this.#a)})}if(!i)return!1;if(this.#w.forget(),!this.#st(i,e))return!0;let a=this.#xr();return a?.atoms.length&&this.#En(a,{kind:`atom`,index:a.atoms.length-1}),!0}#kn(e){let t=this.#e.view;if(!t||this.#T||this.#S)return!1;let n=e.docLen,r={x:Math.max(0,...e.atoms.map(e=>e.at.x))+10*n,y:0},i=this.#An(e,r);if(!i)return!1;let a=i.ops.flatMap(e=>e.op===`addAtom`?[{x:e.x-r.x,y:e.y-r.y}]:[]),o=Oe(a)??{left:-n,bottom:-n,right:n,top:n},s=e.atoms.filter(e=>!e.hidden),c=_r(o,s.length?ir(e.atoms.map(e=>e.at),e.bonds):[],{x:t.cx,y:t.cy},n),l=this.#An(e,c);if(!l)return!1;if(this.#w.forget(),!this.#st(l,e))return!0;let u=this.#xr();return u?.atoms.length&&this.#En(u,{kind:`atom`,index:u.atoms.length-1}),!0}#An(e,t){let n=this.#xe();if(!n)return null;let r=this.#sn();if(r)return r.onCanvas(e,t,0,!1);let i=this.#cn();if(i?.kind===`group`)return _o(e,null,t,i.label);if(i?.kind===`salt`){let n=this.#ln(i.salt);return n&&To(e,null,t,i.salt,n)}return Ur(n,e,t)}#jn(e){let t=this.#Sr(e),n=ic(this.#Q.indices(e),_m(e),t,t=>xm(e,t));this.#In(),n>=0&&this.#En(e,{kind:`atom`,index:n})}#Mn(e,t){let n=t.kind===`atom`?e.atoms[t.index].id:e.bonds[t.index].id;this.setTool(`select.rect`);let r=this.#xr();if(!r)return;let i=t.kind===`atom`?r.atomById(n):r.bondById(n);i<0||(this.#Q.replace(r,Kp(r,t.kind===`atom`?{atoms:[i],bonds:[]}:{atoms:[],bonds:[i]})),this.#fr(r))}#Nn(e,t){let n=this.#e.atomAttachmentPoints(t)??0,r=n&1?n&2?-1:n|2:n|1;if(r<0){let n=this.#a.noFreeAttachment(em(this.#a,tm(e,t),t));this.#lt({atoms:[t],bonds:[],message:n,pointer:!1,subject:this.#yt(e)},e);return}this.#st(Rp(t,n,r),e)}#Pn(e,t,n){if(t.kind===`element`){let r=n.atoms.flatMap(t=>e.labelGroupOf(t)??[]);if(!r.length)this.#st(We(n.atoms,t.symbol),e);else{let i=n.atoms.filter(t=>!e.labelGroupOf(t)),a=[...We(i,t.symbol)?.ops??[],...Nn(r,t.symbol)];this.#st(S(`atom.element`,a),e)}}else t.kind===`bond`&&this.#st(Cr(e,n.bonds,t.type,this.#at),e)}#Fn(e){if(this.#Q.empty||this.#S)return!1;let t=this.#xr();if(!t)return!1;let n=this.#Q.indices(t);return n.atoms.length?(this.#w.forget(),this.#W&&!this.#W.pointer&&this.#bt(null),this.#Pn(t,{kind:`element`,symbol:e},n),!0):!1}#In(){let e=this.#Q.clear()||this.#te!==null;if(this.#te=null,this.#ne=null,!e)return;let t=this.#xr();t?this.#fr(t):this.#n.hideSelection()}#Ln(){this.setTool(`select.rect`);let e=this.#xr();e&&(this.#Q.replace(e,qp(e)),this.#fr(e))}#Rn(e,t){let n=this.#xe(),r=!!this.#k.shift,i=t.hit;if(!i){r||this.#Q.clear(),this.#fr(e);return}let a=t.count>=2;this.#zn(e,i,r,a||n?.kind===`select`&&n.mode===`structure`,a)}#zn(e,t,n,r,i){let a=t.kind===`atom`?{atoms:[t.index],bonds:[]}:{atoms:[],bonds:[t.index]},o=r?Kp(e,a):a;n?i?this.#Q.add(e,o):this.#Q.toggle(e,o):this.#Q.replace(e,o),this.#fr(e)}#Bn(e,t,n){if(!n.tx)return null;let r=t=>e.atoms[n.merged.get(t)??t]?.id??-1,i=this.#Q.atomsInOrder(e),a=this.#Q.indices(e).bonds.map(t=>e.bonds[t]),o=e=>[r(e.a),r(e.b)];return t.copy?{atoms:i.filter(e=>n.merged.has(e)).map(r),pairs:a.filter(e=>n.merged.has(e.a)&&n.merged.has(e.b)).map(o),before:{atoms:new Set(e.atoms.map(e=>e.id)),bonds:new Set(e.bonds.map(e=>e.id))}}:{atoms:i.map(r),pairs:a.map(o),before:null}}#Vn(e){let t=this.#xr();if(!t)return;let n=e.before,r=new Set(t.atoms.map(e=>e.id)),i=[...e.atoms,...n?t.atoms.filter(e=>!n.atoms.has(e.id)).map(e=>e.id):[]].filter(e=>r.has(e)),a=(e,t)=>e<t?`${e} ${t}`:`${t} ${e}`,o=new Set(e.pairs.map(([e,t])=>a(e,t))),s=t.bonds.filter(e=>n&&!n.bonds.has(e.id)||o.has(a(t.atoms[e.a].id,t.atoms[e.b].id))).map(e=>e.id);this.#Q.setIds([...new Set(i)],s),this.#fr(t)}#Hn(e,t){let n=e.wholeGroups(this.#Q.moved(e)),r=e.docLen/20;n.length&&this.#st(pu(e,n,{x:t.x*r,y:t.y*r}),e)&&this.#Un(t.x,-t.y)}#Un(e,t){let n=this.#xr(),r=this.#xn();if(!n||!r||!this.#r.pan)return;let i=this.#Q.moved(n).map(e=>n.atoms[e].px);if(!i.length)return;let[a,o]=[i.map(e=>e.x),i.map(e=>e.y)],s={left:Math.min(...a),top:Math.min(...o),right:Math.max(...a),bottom:Math.max(...o)},c=Fe(s,r,e,t);(c.dx||c.dy)&&this.#r.pan(c.dx,c.dy)}#Wn(e){let t=mu(e,this.#Q.atomsInOrder(e),this.#pn().wedge);t.tx?this.#st(t.tx,e):t.refused&&t.blocker>=0&&this.#st(null,e,[{atom:t.blocker,symbol:tm(e,t.blocker)}],`merge`)}rotateSelection(e,t){let n=this.#xr();if(!n||this.#S||this.#ye||this.#Q.empty)return!1;let r=n.wholeGroups(this.#Q.moved(n));return r.length?this.#st(ju(n,r,this.#qn(n,r).centre,t===`cw`?-e:e),n,[],`bond`,!1,!0):!1}#Gn(e,t){let n=e.wholeGroups(this.#Q.moved(e));n.length&&this.#st(ju(e,n,this.#qn(e,n).centre,t),e,[],`bond`,!1,!0)}#Kn(){return JSON.stringify(this.#Q.ids())}#qn(e,t){let n=this.#ue;if(n&&n.key===this.#Kn()){if(n.atomId===null)return{centre:{...n.at},atom:-1};let t=e.atomById(n.atomId);if(t>=0)return{centre:{...e.atoms[t].at},atom:t}}let r=Ou(e,t);return r?{centre:r.centre,atom:r.atom}:{centre:{x:0,y:0},atom:-1}}#Jn(e){if(this.#T||this.#ye)return null;let t=this.#Q.moved(e);if(t.length<2)return null;let n=this.#xn();if(!n)return null;let r=Math.max(10,.4*e.bondPx+2),i=t.flatMap(t=>e.atoms[t].labelBox?[e.atoms[t].labelBox]:[]),a=ga(t.filter(t=>!e.atoms[t].hidden).map(t=>e.atoms[t].px),r,i);return a?ba(a,n,this.#Yn(e)):null}#Yn(e){return{points:e.atoms.filter(e=>!e.hidden).map(e=>e.px),reach:Math.max(12,Do*e.bondPx),labels:e.atoms.flatMap(e=>e.labelBox?[e.labelBox]:[])}}#Xn(e){if(this.#ye)return;let t=this.#S,n=!!(t&&t.dragging&&!t.cancelled),r=n?t?.grip??null:null,i=r===`rotate`?t?.turn??null:null,a=this.#be||n&&!r?null:i?i.layout:this.#Jn(e);if(!a){this.#oe=null,this.#n.hideChrome(),this.#Zn(null);return}let o=i?i.centre:this.#qn(e,this.#Q.moved(e)).centre,s=e.toCanvas(o),c={at:a.handle,stem:a.stem},l=null;if(i){let e=this.#fe,t=nm(a.handle,s,e);c={at:t,stem:a.stem&&[nm(a.stem[0],s,e),nm(a.stem[1],s,e)]};let n=this.#xn()?.width??1/0,r=this.#a.degrees(Math.round(e)||0),i=t.x+12+6+8+7*r.length>n-4;l={at:{x:i?t.x-12-6:t.x+12+6,y:t.y},text:r,leftward:i}}r===`centre`&&this.#de&&(s=this.#de.atom>=0?e.atoms[this.#de.atom].px:e.toCanvas(this.#de.at)),this.#n.showChrome({handle:c,centre:s,angle:l}),this.#oe=n?null:{layout:a,centre:s},this.#Zn(n?null:a.actions)}#Zn(e){rm(e,this.#le)||(this.#le=e?{x:e.x,y:e.y}:null,this.#r.selectionActions?.(this.#le))}#Qn(e,t,n,r,i){let a=this.#oe;return a?wa({at:t,touch:n,alt:r,onAtom:i?.kind===`atom`},{handle:a.layout.handle,centre:a.centre},e.atoms.map(e=>e.hidden?lr:e.px)):null}#$n(e,t){let n=this.#oe,r=e.wholeGroups(this.#Q.moved(e));if(!n||!r.length)return null;this.#fe=0;let{centre:i,atom:a}=this.#qn(e,r);return{atoms:r,centre:i,atom:a,from:wu(i,e.toMolecule(t)),snaps:a>=0?ku(e,r,a):null,layout:n.layout}}#er(e,t,n){let r=t.turn;return r?Au(wu(r.centre,e.toMolecule(n))-r.from,!!(this.#k.alt||this.#k.ctrl),r.snaps):0}#tr(e,t,n){let r=t.turn;if(!r)return;let i=this.#er(e,t,n);this.#fe=i;let a=ju(e,r.atoms,r.centre,i);a?this.#in(a,e,{opacity:1,judge:!1}):this.#Xt();let o=e.toCanvas(r.centre);this.#fr(e,dm,i?e=>nm(e,o,i):null)}#nr(e,t){let n=this.#wn(e,t);this.#de=n>=0?{atom:n,at:{...e.atoms[n].at}}:{atom:-1,at:e.toMolecule(t)},this.#Xn(e)}#rr(e,t){let n=this.#wn(e,t);this.#ue={key:this.#Kn(),atomId:n>=0?e.atoms[n].id:null,at:n>=0?{...e.atoms[n].at}:e.toMolecule(t)},this.#de=null}#ir(){if(!this.#r.pan||this.#ye)return;let e=this.#y,t=this.#xn();if(!e||!t||!Pn(e,t)){this.#or();return}if(this.#ae)return;let n=this.#ie;n&&this.#re&&Math.hypot(e.x-n.x,e.y-n.y)<=6||(this.#re?.(),this.#ie={...e},this.#re=this.#z.after(300,()=>{this.#re=null,this.#ie=null,this.#ae=!0,this.#ar()}))}#ar(){this.#re=null;let e=this.#S,t=this.#y,n=this.#xe(),r=this.#xn(),i=!!(n&&e&&um(n,e)),a=e&&!e.cancelled&&e.dragging&&t&&r&&this.#r.pan&&!this.#ye&&(i||cm(n,e))?Pn(t,r):null;if(!e||!t||!a){this.#or();return}this.#r.pan?.(a.dx,a.dy);let o=this.#xr();if(!o||this.#S!==e){this.#or();return}i?this.#lr(o,e,t):this.#tn(),this.#re=this.#z.frame(()=>this.#ar())}#or(){this.#re?.(),this.#re=null,this.#ie=null,this.#ae=!1}#sr(e,t,n){let r=this.#xe();return r?.kind===`select`&&r.mode===`lasso`!=!!this.#k.alt?{kind:`lasso`,path:t.pathMol.map(t=>e.toCanvas(t))}:{kind:`rect`,from:e.toCanvas(t.startMol),to:n}}#cr(e,t){return Gp({atoms:e.atoms.map(e=>e.hidden?lr:e.px),bonds:e.bonds.map(e=>({a:e.ends[0],b:e.ends[1]}))},t)}#lr(e,t,n){let r=this.#sr(e,t,n);this.#n.showArea(r);let i=this.#cr(e,r),a=this.#s===`erase`,o=a?ur(e,i.atoms,i.bonds):this.#k.shift?fm(this.#Q.indices(e),i):i,s=JSON.stringify(o);if(s!==this.#ne){if(this.#ne=s,a){this.#br(e,o);return}this.#te=o,this.#fr(e)}}#ur(e,t,n){let r=this.#cr(e,this.#sr(e,t,n));this.#k.shift?this.#Q.add(e,r):this.#Q.replace(e,r)}#dr(e,t,n){let r=this.#cr(e,this.#sr(e,t,n));return mt(e,r.atoms,r.bonds)}#fr(e,t=dm,n=null){this.#gr(e,t,n),!this.#S&&this.#x&&this.#Hr(e),this.#Xn(e),this.#pr()}#pr(){if(this.#ye)return;let e=this.#mr();e!==this.#hr&&(this.#hr=e,this.#r.stereoAvailable?.(e));let t=JSON.stringify(this.#Q.ids());t!==this.#me&&(this.#me=t,this.#r.selectionUpdated?.());let n=!this.#Q.empty;n!==this.#pe&&(this.#pe=n,this.#r.selectionChanged?.(n))}#mr(){let e=this.#xr();return e?pe(this.#e.stereo(),this.#Q.empty?null:this.#Q.indices(e).atoms):!1}#hr=!1;#gr(e,t,n){this.#ce+=1;for(let e of this.#ee)e.removeAttribute(`data-selected`);this.#ee=[];let r=this.#te??this.#Q.indices(e),i=r.atoms.filter(n=>e.atoms[n]&&!t.has(n)),a=r.bonds.filter(n=>e.bonds[n]&&!t.has(e.bonds[n].a)&&!t.has(e.bonds[n].b));if(!i.length&&!a.length){this.#n.hideSelection();return}let o=this.#vr();for(let e of i)this.#_r(o.atoms.get(e));for(let e of a)this.#_r(o.bonds.get(e));let s=this.#Kr(e,i,a);n?typeof n==`function`?this.#n.showSelection(pm(s,n)):this.#n.showSelection(pm(s,e=>({x:e.x+n.x,y:e.y+n.y}))):this.#n.showSelection(s)}#_r(e){e&&(e.setAttribute(`data-selected`,`true`),this.#ee.push(e))}#vr(){if(this.#yr)return this.#yr;let e=new Map,t=new Map;for(let n of Array.from(this.#t.svg?.querySelectorAll(`g[data-atom], g[data-bond]`)??[])){let r=n.getAttribute(`data-atom`);r===null?t.set(Number(n.getAttribute(`data-bond`)),n):e.set(Number(r),n)}return this.#yr={atoms:e,bonds:t},this.#yr}#yr=null;#br(e,t){for(let e of this.#C)e.removeAttribute(`data-hovered`);if(this.#C=[],!t.atoms.length&&!t.bonds.length){this.#n.hideHover();return}let n=this.#vr();for(let e of t.atoms)this.#qr(n.atoms.get(e));for(let e of t.bonds)this.#qr(n.bonds.get(e));this.#n.showHover({...this.#Kr(e,t.atoms,t.bonds),quiet:!1})}#xr(){if(this.#_)return this.#_;let e=this.#t.svg,t=this.#e.view,n=e&&t?this.#e.geometry():null;if(!e||!t||n===null)return null;let r=am(t,Number.parseFloat(e.getAttribute(`width`)??``)||0,Number.parseFloat(e.getAttribute(`height`)??``)||0);return this.#_=new hn(Yn(n),yn(e),t.docLen,r,e=>me(this.#e.atomJson(e)),wn(this.#e.sgroups()),e=>kt(this.#e.bondQueryJson(e)),()=>({perceived:this.#e.aromaticRingBonds(),asRead:!this.#e.edited})),this.#_}#Sr(e){let t=this.#x;if(!t)return null;let n=t.kind===`atom`?e.atomById(t.id):e.bondById(t.id);return n>=0?{kind:t.kind,index:n}:null}#Cr(e,t){let n=Fo(e.hitAtoms(),e.hitBonds(),t,e.bondPx,this.#o);return this.#b=n?{kind:n.kind,id:n.kind===`atom`?e.atoms[n.index].id:e.bonds[n.index].id}:null,this.#b&&!this.#Cn(t)&&(this.#x=this.#b),this.#kt(e),this.#Hr(e),this.#Rr(e,n),this.#Fr(e,n),this.#Er(e,t,n),n}#wr(e,t,n,r=!1){if(!e.sgroups.length)return-1;let i=e.groupAt(t,r?Po.bond:0);return i<0||n?.kind===`atom`&&e.sgroups[i]?.type!==`DAT`?-1:i}#Tr=null;#Er(e,t,n){let r=!this.#S&&!this.#he&&!this.#T&&!this.#E&&!this.#j?this.#wr(e,t,n):-1,i=r>=0?e.sgroups[r]:null;if(!i){this.#Dr();return}let a=Af(i,this.#a);(this.#Tr?.k!==r||this.#Tr.text!==a)&&(this.#Tr={k:r,text:a},this.#r.sgroupTip?.(a,t))}#Dr(){this.#Tr&&(this.#Tr=null,this.#r.sgroupTip?.(null,null))}#Or(e,t,n){let r=t.group??-1,i=e.sgroups[r];if(i){if(this.#w.clicked({atomId:null,token:null,at:t.start,time:this.#D(),touch:t.touch}),n?.kind===`erase`){this.#st(Df(r),e,[],`bond`,!0);return}if(n?.kind===`sgroup`||t.count>=2){df(i.type)&&this.#r.openSGroup?.({kind:`edit`,group:r,token:this.#e.stateToken});return}if(this.#y&&this.#Er(e,this.#y,null),n?.kind===`select`){let t=e.groupItems(r);this.#k.shift?this.#Q.add(e,t):this.#Q.replace(e,t),this.#fr(e)}}}#kr(e,t){let n=t.hit;if(n){if(n.kind===`atom`){let t=e.labelGroupOf(n.index);if(t){this.#Dr(),this.#r.openSGroup?.({kind:`edit`,group:t.index,token:this.#e.stateToken});return}this.#jr(e,[n.index],[])}else this.#jr(e,[],[n.index])}}#Ar(e,t,n){let r=this.#cr(e,this.#sr(e,t,n));(r.atoms.length||r.bonds.length)&&this.#jr(e,r.atoms,r.bonds)}#jr(e,t,n){this.#Dr(),this.#r.openSGroup?.({kind:`new`,atomIds:t.filter(t=>e.atoms[t]).map(t=>e.atoms[t].id),bondIds:n.filter(t=>e.bonds[t]).map(t=>e.bonds[t].id),token:this.#e.stateToken})}#Mr(e,t){let n=t.atomIds.map(t=>e.atomById(t)),r=t.bondIds.map(t=>e.bondById(t));if(n.some(e=>e<0)||r.some(e=>e<0))return null;let i=new Set(n);for(let t of r)i.add(e.bonds[t].a).add(e.bonds[t].b);let a=e.wholeGroups([...i]);return a.length?{atoms:a,bonds:r}:null}#Nr(e){return{atoms:e.atoms.length,bonds:e.bonds,position:t=>e.atoms[t]?.at??{x:0,y:0},docLen:e.docLen}}sgroupOpening(e){let t=this.#xr();if(!t||e.token!==this.#e.stateToken)return null;if(e.kind===`edit`){let n=t.sgroups[e.group];if(!n)return null;let r=mf(n,_f(this.#Nr(t),n.atoms??[],n.type===`DAT`?n.bonds??[]:[]));return r?{form:r,editing:!0}:null}let n=this.#Mr(t,e);if(!n)return null;let r=e.bondIds.length===1&&n.atoms.length===2||n.bonds.length?n.bonds:[];return{form:ff(_f(this.#Nr(t),n.atoms,r)),editing:!1}}applySGroup(e,t){let n=this.#xr();if(!n||this.#S||e.token!==this.#e.stateToken)return null;let r=this.#Nr(n),i;if(e.kind===`edit`){let a=n.sgroups[e.group];if(!a)return null;i=Ef(r,a,t,this.#a)}else{let a=this.#Mr(n,e);if(!a)return null;i=Tf(r,a,t,this.#a)}if(!i)return null;let a=this.#e.sgroupRefusal(i);return a===null?(this.#st(i,n),null):Of(a,t,this.#a)}#Pr(e,t){let n=this.#Je(e,t).filter(t=>!e.hiddenAtoms.has(t)),r=kf(e.sgroups,n)[0];return r?r.index:null}menuSGroup(e,t){let n=this.#xr();if(!n||this.#S)return!1;let r=this.#Pr(n,this.#$e(n,e));return r===null?!1:t===`remove`?this.#st(Df(r),n):(df(n.sgroups[r].type)&&this.#r.openSGroup?.({kind:`edit`,group:r,token:this.#e.stateToken}),!1)}#Fr(e,t){let n=t?.kind===`atom`&&!this.#S&&!this.#he&&!this.#T&&!this.#E&&!this.#j&&!e.labelGroupOf(t.index)?t.index:-1,r=n>=0?e.info(n):null,i=r?ee(r.query??null,r.symbol,this.#a):null,a=i?e.atoms[n].id:-1;if(this.#ve&&this.#ve.atomId===a||(this.#Lr(),!i))return;let o={atomId:a,cancel:null,shown:!1};this.#ve=o,o.cancel=this.#z.after(300,()=>{o.cancel=null,this.#ve===o&&!this.#S&&this.#Ir(o)})}#Ir(e){let t=this.#xr(),n=t?.atomById(e.atomId)??-1,r=this.#xn();if(!t||n<0||!r||this.#ye)return;let i=t.atoms[n],a=t.info(n),o=a?ee(a.query??null,a.symbol,this.#a):null;o&&(this.#n.showQueryTip(i.labelBox??{x:i.px.x,y:i.px.y,width:0,height:0},o,r,`; `),e.shown=!0)}#Lr(){let e=this.#ve;e&&(this.#ve=null,e.cancel?.(),e.shown&&this.#n.hideQueryTip())}#Rr(e,t){let n=t?.kind===`atom`&&!this.#S&&!this.#he&&!this.#T&&!this.#E&&!this.#j?e.labelGroupOf(t.index):null,r=n?e.atoms[n.labelAtom].id:-1;if(this.#_e&&this.#_e.atomId===r||(this.#Br(),!n))return;let i={atomId:r,cancel:null,shown:!1};this.#_e=i,i.cancel=this.#z.after(300,()=>{i.cancel=null,this.#_e===i&&!this.#S&&this.#zr(i)})}#zr(e){let t=this.#xr(),n=t?.atomById(e.atomId)??-1,r=t&&n>=0?t.labelGroupOf(n):null,i=this.#xn(),a=t&&n>=0?t.atoms[n].labelBox??{x:t.atoms[n].px.x,y:t.atoms[n].px.y,width:0,height:0}:null;if(!t||!r||!i||!a||this.#ye)return;let o=this.#e.groupPreview(r.atoms,st(r.points),t.docLen);o&&(this.#n.showAbbreviation(a,o,i,t.bondPx/qr),e.shown=!0)}#Br(){let e=this.#_e;e&&(this.#_e=null,e.cancel?.(),e.shown&&this.#n.hideAbbreviation())}#Vr(){for(let e of this.#C)e.removeAttribute(`data-hovered`);this.#C=[],this.#n.hideHover()}#Hr(e,t){this.#Ur(e,t),this.#r.hotChanged&&this.#Ze(e)}#Ur(e,t){for(let e of this.#C)e.removeAttribute(`data-hovered`);if(this.#C=[],this.#he&&!t){let t=this.#Qe(e);if(!t){this.#n.hideHover();return}this.#qr(this.#t.svg?.querySelector(`g[data-${t.kind}="${t.index}"]`));let[n,r]=[t.kind===`atom`?[t.index]:[],t.kind===`bond`?[t.index]:[]];this.#n.showHover({...this.#Kr(e,n,r),quiet:!1,selected:this.#Gr(e,n,r)});return}let n=t??this.#Sr(e);!n&&this.#x&&!t&&(this.#x=null);let r=!!t||this.#ot();this.#Ot(e,n,r||this.#O);let i=this.#W;if(i&&!i.notice){let t=n&&(r||this.#O)?n:null,a=t?.kind===`bond`&&e.bonds[t.index]?[e.bonds[t.index].a,e.bonds[t.index].b]:[],o=i.atoms.filter(e=>!a.includes(e)),s=[...new Set([...t?.kind===`atom`?[t.index]:[],...o])].filter(t=>e.atoms[t]),c=[...new Set([...t?.kind===`bond`?[t.index]:[],...i.bonds])].filter(t=>e.bonds[t]);if(s.length||c.length){let n=this.#t.svg;t&&this.#qr(n?.querySelector(`g[data-${t.kind}="${t.index}"]`)),this.#n.showHover({...this.#Kr(e,s,c),warn:!0});return}}if(!n){this.#n.hideHover();return}if(!r&&!this.#O){this.#n.hideHover();return}let a=n.kind===`atom`?[n.index]:[],o=n.kind===`bond`?[n.index]:[];if(this.#s===`erase`&&!t&&(n.kind===`atom`?{atoms:a,bonds:o}=Xn(e,n.index):(this.eraseWhole||!cn[e.bonds[n.index].order])&&({atoms:a,bonds:o}=Tn(e,n.index))),this.#s===`select.structure`&&!t&&r){let t=Kp(e,{atoms:a,bonds:o});[a,o]=[[...t.atoms],[...t.bonds]]}let s=this.#t.svg;for(let e of a)this.#qr(s?.querySelector(`g[data-atom="${e}"]`));for(let e of o)this.#qr(s?.querySelector(`g[data-bond="${e}"]`));let c=this.#Gr(e,a,o),l=`${a.join(`,`)}|${o.join(`,`)}|${r}|${c}`,u=this.#Wr;u&&u.model===e&&u.key===l&&Number.isFinite(u.version)&&u.version===this.#n.hoverVersion||(this.#n.showHover({...this.#Kr(e,a,o),quiet:!r,selected:c}),this.#Wr={model:e,key:l,version:this.#n.hoverVersion})}#Wr=null;#Gr(e,t,n){if(!t.length&&!n.length)return!1;let r=this.#te;if(!r)return this.#Q.holds(e,{atoms:t,bonds:n});let[i,a]=[new Set(r.atoms),new Set(r.bonds)];return t.every(e=>i.has(e))&&n.every(e=>a.has(e))}#Kr(e,t,n){let r=t.filter(t=>!e.atoms[t]?.hidden),i=n.filter(t=>!e.bonds[t]?.hidden);return{atoms:r.map(t=>({...e.atoms[t].px,label:e.atoms[t].labelBox})),bonds:i.map(t=>[e.atoms[e.bonds[t].ends[0]].px,e.atoms[e.bonds[t].ends[1]].px]),bondPx:e.bondPx,short:i.map(t=>e.drawnLength(t)<$c*e.bondPx)}}#qr(e){e&&(e.setAttribute(`data-hovered`,`true`),this.#C.push(e))}#Jr(e,t,n){if(t.hit?.kind===`bond`)return null;let r=t.hit?.kind===`atom`?t.hit.index:-1,i=new Set(r>=0?[r]:[]),a=r>=0?{atom:r}:{point:e.toMolecule(t.start)},o=r>=0?e.atoms[r].px:t.start,s=`atom`in a?e.atoms[a.atom].at:a.point,c=Fo(e.hitAtoms(),[],n,e.bondPx,{excludeAtoms:i,atomsOnly:!0});if(c)return this.#Yr({from:a,to:{atom:c.index},fromPx:o,toPx:e.atoms[c.index].px},s,e.atoms[c.index].at);let l=!!this.#k.alt,u=l?e.toMolecule(n):O(s,Ir(im(t.start,n)),e.docLen);if(l&&Math.hypot(u.x-s.x,u.y-s.y)<.3*e.docLen)return null;let d=Cn(e.atoms.map(e=>e.at),u,on*e.docLen,e.hiddenAtoms.size?new Set([...i,...e.hiddenAtoms]):i);return d>=0?this.#Yr({from:a,to:{atom:d},fromPx:o,toPx:e.atoms[d].px},s,e.atoms[d].at):this.#Yr({from:a,to:{point:u},fromPx:o,toPx:e.toCanvas(u)},s,u)}#Yr(e,t,n){return this.#M={kind:`draw`,length:Math.hypot(n.x-t.x,n.y-t.y),angle:Math.atan2(n.y-t.y,n.x-t.x)*180/Math.PI},e}drawingReaches(e){let t=this.#xr();return t?zo(t.hitAtoms(),t.hitBonds(),e,t.bondPx,this.#o):!1}pointerFacts(){let e=this.#S;if(e)return e.dragging&&!e.cancelled?this.#M:null;let t=this.#b;if(!t||!this.#y||this.#j||this.#T)return null;let n=this.#xr();if(!n)return null;if(t.kind===`atom`){let e=n.atomById(t.id),r=e>=0?n.info(e):null;return r?{kind:`atom`,number:e+1,element:r.alias??r.symbol,charge:r.charge,h:r.hCount,valence:r.totalValence??null}:null}let r=n.bondById(t.id),i=r>=0?n.bonds[r]:null;if(!i)return null;let[a,o]=[n.atoms[i.a].at,n.atoms[i.b].at];return{kind:`bond`,number:r+1,type:Ns(this.#a,Ms(i.order,i.stereo)),length:Math.hypot(o.x-a.x,o.y-a.y)}}};function Cm(e,t,n){let r=ei(e,n);return!!(r&&r.status===`specified`&&t.bonds[n]?.order===`double`&&t.bonds[n]?.stereo!==`crossed`)}var wm=Symbol(`Crux Sketch: destroyed`),Tm=()=>{},Em=class{#e;#t=!1;#n=new Set;#r=new Set;#i=new Set;constructor(e){this.#e=e}get ended(){return this.#t}get pending(){return{jobs:this.#n.size,attached:this.#r.size,waits:this.#i.size}}get busy(){if(this.#i.size)return!0;for(let e of this.#n)if(!e.quiet)return!0;return!1}frame(e){return this.#a(e=>this.#e.frame(e),e)}painted(e,t=!1){return this.#a(e=>this.#e.painted(e),e,t)}after(e,t){return this.#a(t=>this.#e.after(e,t),t)}expire(e,t){return this.#a(t=>this.#e.after(e,t),t,!0)}attach(e){if(this.#t)return e(),Tm;let t=!1,n=()=>{t||(t=!0,this.#r.delete(n),e())};return this.#r.add(n),n}alive(e){return this.#t?(e.catch(Tm),Promise.reject(wm)):new Promise((t,n)=>{this.#i.add(n),e.then(e=>{this.#i.delete(n),t(e)},e=>{this.#i.delete(n),n(e)})})}end(){if(this.#t)return;this.#t=!0;let e=[...this.#n];this.#n.clear();for(let t of e)t.cancel();let t=[...this.#r];for(let e of t)try{e()}catch{}this.#r.clear();let n=[...this.#i];this.#i.clear();for(let e of n)e(wm)}#a(e,t,n=!1){if(this.#t)return Tm;let r={cancel:Tm,quiet:n};return this.#n.add(r),r.cancel=e(()=>{!this.#n.delete(r)||this.#t||t()}),()=>{this.#n.delete(r)&&r.cancel()}}},Dm=[`any`,`generic`,`atom`,`bond`,`selection`];function Om(e,t){return e===`any`?t.sheetAnywhere:e===`generic`?t.sheetNothingHovered:e===`atom`?t.sheetOverAtom:e===`bond`?t.sheetOverBond:t.sheetSelection}var km={1:`sprout1`,2:`sprout2`,3:`sprout3`,4:`sprout4`,5:`sprout5`,6:`sprout6`,7:`sprout7`,8:`sprout8`,9:`sprout9`,0:`sprout0`,a:`sprout3`,K:`sproutTBu`,v:`sproutCyclopropane`,u:`sproutCyclobutane`,z:`sproutAlkyne`,Z:`sproutAzide`,M:`sproutMgBr`,k:`sproutSulfonyl`},Am={"bond.single":`singleBond`,"bond.double":`doubleBond`,"bond.triple":`tripleBond`,"bond.wedge":`wedgeBond`,"bond.hash":`hashedWedgeBond`,"bond.wavy":`wavyBond`,"bond.left":`sheetPositionLeft`,"bond.centre":`sheetPositionCentre`,"bond.right":`sheetPositionRight`,clear:`clearCanvas`,undo:`undo`,redo:`redo`,"atom.deuterium":`sheetDeuterium`,"atom.unlabel":`sheetUnlabel`,"atom.properties":`atomProperties`,"atom.edit":`editLabel`,"atom.rgroup":`rgroupLabel`,"atom.nickname":`chooseNickname`,"atom.generic":`sheetGenericX`,"atom.charge.plus":`chargePlus`,"atom.charge.minus":`chargeMinus`,"zoom.in":`zoomIn`,"zoom.out":`zoomOut`,"zoom.reset":`zoomReset`,"zoom.fit":`zoomFit`,layout:`layout`,clean:`cleanUp`,"select.all":`selectAll`,"select.none":`sheetSelectNone`,"select.last":`sheetSelectLast`,"selection.delete":`deleteSelection`,"select.join":`stepJoin`,"flip.horizontal":`flipHorizontal`,"flip.vertical":`flipVertical`,"rotate180.horizontal":`sheetTurnHorizontal`,"rotate180.vertical":`sheetTurnVertical`,"structure.aromatize":`aromatize`,"structure.dearomatize":`dearomatize`,cut:`cut`,copy:`copy`,paste:`paste`,"bond.properties":`bondProperties`,"menu.open":`sheetMenu`,"hotspot.walk":`sheetWalk`,"hotspot.jump":`sheetJump`,"hotspot.hop":`sheetHop`,"hotspot.start":`sheetStart`,"structure.new":`sheetNewStructure`,"selection.end":`sheetEndSelection`,"select.hotspot":`selectStructure`,"atom.attachment":`attachmentPoint`,"keys.open":`keyboardShortcuts`,"settings.open":`settings`,"library.open":`structureLibrary`};function jm(e,t){switch(e.action){case`tool.bond.single`:case`tool.ring.benzene`:case`tool.ring.cyclopentadiene`:case`tool.chain`:case`tool.hand`:case`tool.choose`:{let n=e.action===`tool.choose`?e.tool:{"tool.bond.single":`bond.single`,"tool.ring.benzene":`ring.benzene`,"tool.ring.cyclopentadiene":`ring.cyclopentadiene`,"tool.chain":`chain`,"tool.hand":`hand`}[e.action];return n?t.format(`toolChosen`,{tool:D(n,t)}):e.action}case`tool.cycle`:return t.format(`sheetInTurn`,{tools:(e.tools??[]).map(e=>D(e,t)).join(`, `)});case`atom.element`:return e.element?t.elementName(e.element):e.action;case`atom.sprout`:return t[km[e.sprout??``]??`sprout1`];case`atom.bond`:return D(`bond.${e.bondType??`single`}`,t);case`bond.fuse`:return t.format(`sheetFuse`,{ring:e.ring?D(`ring.${e.ring}`,t):``});case`selection.nudge`:return t.format(`sheetNudge`,{n:Math.hypot(e.nudge?.x??0,e.nudge?.y??0)});case`selection.rotate`:return t.format(`sheetRotate`,{n:e.turn??0});case`copy.as`:return t.format(`sheetCopyAs`,{format:t.formatName(e.copy??``)});case`paste.as`:return t.format(`sheetPasteAs`,{format:t.formatName(e.paste??``)});case`atom.abbreviation`:return t.format(`sheetAbbreviation`,{label:e.abbreviation??``});default:{let n=Am[e.action];return n?t[n]:e.action}}}function Mm(e,t,n=`Ctrl`){let r=new Map(Dm.map(e=>[e,new Map]));for(let i of e){let e=r.get(i.scope);if(!e)continue;let a=jm(i,t),o=e.get(a)??[],s=us(i,n);o.includes(s)||o.push(s),e.set(a,o)}return Dm.flatMap(e=>{let n=[...(r.get(e)??new Map).entries()].map(([e,t])=>({what:e,keys:t}));return n.length?[{scope:e,title:Om(e,t),rows:n}]:[]})}var Nm=e=>{let t=typeof e==`string`?e.trim().toLowerCase():``;return t===`collapsed`||t===`hidden`?t:`shown`},Pm=e=>typeof e==`number`&&Number.isFinite(e)?e:0,Fm=e=>Array.isArray(e)?e.filter(e=>Number.isInteger(e)&&e>=0):[];function Im(e){if(!e)return null;let t;try{t=JSON.parse(e)}catch{return null}if(!t||typeof t!=`object`)return null;if(t.query===!0)return{query:!0,formula:``,averageMass:0,exactMass:0,fragments:[]};if(typeof t.formula!=`string`)return null;let n=(Array.isArray(t.fragments)?t.fragments:[]).map(e=>({atoms:Fm(e?.atoms),formula:typeof e?.formula==`string`?e.formula:``,averageMass:Pm(e?.averageMass),exactMass:Pm(e?.exactMass),charge:Pm(e?.charge),radicalElectrons:Pm(e?.radicalElectrons),rSites:Pm(e?.rSites),dummies:Pm(e?.dummies)}));return{query:!1,formula:t.formula,averageMass:Pm(t.averageMass),exactMass:Pm(t.exactMass),fragments:n}}function Lm(e){return!e||e.query?null:Object.freeze({formula:e.formula,molecularWeight:e.averageMass,exactMass:e.exactMass,fragments:Object.freeze(e.fragments.map(e=>Object.freeze({atoms:Object.freeze([...e.atoms]),formula:e.formula,molecularWeight:e.averageMass,exactMass:e.exactMass})))})}function Rm(e){let t=[],n=(e,n)=>{let r=t.at(-1);r&&r.kind===n?t[t.length-1]={text:r.text+e,kind:n}:t.push({text:e,kind:n})},r=0;for(;r<e.length;){let t=e.slice(r),i=/^\[(\d+)([A-Z][a-z]?)\]/.exec(t);if(i){n(i[1],`sup`),n(i[2],`text`),r+=i[0].length;continue}let a=/^([+-])(\d*)$/.exec(t);if(a){n(`${a[2]&&a[2]!==`1`?a[2]:``}${a[1]===`-`?`−`:`+`}`,`sup`);break}let o=/^\d+/.exec(t);if(o){n(o[0],`sub`),r+=o[0].length;continue}let s=e[r];n(s,s===`•`?`sup`:`text`),r+=1}return t}function zm(e){let t=[];return e.split(`.`).forEach((e,n)=>{n&&t.push({text:`.`,kind:`text`}),t.push(...Rm(e))}),t.reduce((e,t)=>{let n=e.at(-1);return n&&n.kind===t.kind?e[e.length-1]={text:n.text+t.text,kind:t.kind}:e.push(t),e},[])}var Bm=e=>e.toFixed(2),Vm=e=>e.toFixed(4),Hm=e=>e-12,Um=(e,t)=>e+1<=Hm(t);function Wm(e,t){return!t||!(t.width>0)||!(t.height>0)?!0:t.height<200?!1:Um(e,t.width)}function Gm(e,t,n,r=!1){return e===`hidden`||r?`hidden`:e===`collapsed`?`folded`:t||n?`open`:`folded`}function Km(e,t,n){let r=Gm(e,t,n);return r===`hidden`?{state:e,opened:n}:r===`open`?t?{state:`collapsed`,opened:!1}:{state:e,opened:!1}:{state:`shown`,opened:!t}}var qm=(e,t,n)=>e&&t===n,Jm={0:`⁰`,1:`¹`,2:`²`,3:`³`,4:`⁴`,5:`⁵`,6:`⁶`,7:`⁷`,8:`⁸`,9:`⁹`,"+":`⁺`,"−":`⁻`,"-":`⁻`},Ym=e=>e.map(e=>e.kind===`sup`?[...e.text].map(e=>Jm[e]??e).join(``):e.text).join(``),Xm=(e,t)=>Object.freeze({id:e,items:Object.freeze([...t])}),Zm=Object.freeze([Xm(`formula-and-mass`,[`formula`,`molecular-weight`,`exact-mass`]),Xm(`notations`,[`smiles`,`cxsmiles`,`smarts`,`inchikey`]),Xm(`counts`,[`atoms`,`heavy-atoms`,`bonds`,`rings`,`fragments`,`stereocentres`,`charge`]),Xm(`properties`,[`clogp`,`tpsa`,`hbd`,`hba`,`rotatable-bonds`,`lipinski-violations`,`qed`]),Xm(`pointer`,[`under-pointer`,`drawing`])]),Qm=Object.freeze(Zm.flatMap(e=>[...e.items])),$m=Object.freeze([`formula`,`molecular-weight`,`exact-mass`]),eh=new Set(Qm),th=e=>typeof e==`string`&&eh.has(e);function nh(e){let t;if(Array.isArray(e))t=e;else if(typeof e==`string`)t=e.split(/[\s,]+/).map(e=>e.trim().toLowerCase());else return $m;let n=new Set(t.filter(th));return Object.freeze(Qm.filter(e=>n.has(e)))}var rh=(e,t)=>e.length===t.length&&e.every((e,n)=>e===t[n]);function ih(e,t,n){let r=new Set(e);return n?r.add(t):r.delete(t),Object.freeze(Qm.filter(e=>r.has(e)))}var ah=[`CrippenClogP`,`tpsa`,`NumHBD`,`NumHBA`,`lipinskiHBD`,`lipinskiHBA`,`NumRotatableBonds`,`NumHeavyAtoms`,`NumAtoms`,`NumRings`,`NumAtomStereoCenters`,`NumUnspecifiedAtomStereoCenters`,`amw`,`qed`];function oh(e){if(!e)return null;let t;try{t=JSON.parse(e)}catch{return null}if(!t||typeof t!=`object`||Array.isArray(t))return null;let n={query:t.query===!0};for(let e of ah){let r=t[e];n[e]=typeof r==`number`&&Number.isFinite(r)?r:null}return Object.freeze(n)}function sh(e){return[e.amw,e.CrippenClogP,e.lipinskiHBD,e.lipinskiHBA].some(e=>e===null)?null:Number(e.amw>500)+Number(e.CrippenClogP>5)+Number(e.lipinskiHBD>5)+Number(e.lipinskiHBA>10)}var ch=Object.freeze({kind:`empty`}),lh=Ym,uh=e=>e>0?`+${e}`:e<0?`−${-e}`:`0`;function dh(e,t){let n=e.toFixed(t);return/^-0\.?0*$/.test(n)?n.slice(1):n.replace(`-`,`−`)}function fh(e){let t=Math.round(e)%360;return t>180&&(t-=360),t<=-180&&(t+=360),String(t===0?0:t).replace(`-`,`−`)}function ph(e){let t=t=>({name:e.formatName(t),label:e.formatName(t)}),n=e=>({name:e,label:e});return{formula:{name:e.readoutFormula,label:``},"molecular-weight":{name:e.infoWeightName,label:e.infoWeight,labelName:e.infoWeightName},"exact-mass":n(e.infoExactMass),smiles:t(`smiles`),cxsmiles:t(`cxsmiles`),smarts:t(`smarts`),inchikey:t(`inchikey`),atoms:n(e.readoutAtoms),"heavy-atoms":n(e.readoutHeavyAtoms),bonds:n(e.readoutBonds),rings:n(e.readoutRings),fragments:n(e.readoutFragments),stereocentres:n(e.readoutStereocentres),charge:n(e.readoutCharge),clogp:n(e.readoutClogp),tpsa:n(e.readoutTpsa),hbd:{name:e.readoutHbdName,label:e.readoutHbd,labelName:e.readoutHbdName},hba:{name:e.readoutHbaName,label:e.readoutHba,labelName:e.readoutHbaName},"rotatable-bonds":n(e.readoutRotatable),"lipinski-violations":n(e.readoutLipinski),qed:n(e.readoutQed),"under-pointer":{name:e.readoutUnderPointer,label:``},drawing:{name:e.readoutDrawing,label:``}}}function mh(e,t){switch(e){case`formula-and-mass`:return t.infoPanel;case`notations`:return t.readoutNotations;case`counts`:return t.readoutCounts;case`properties`:return t.readoutProperties;default:return t.readoutPointer}}var hh=new Set([`formula`,`molecular-weight`,`exact-mass`]),gh=new Set([`atoms`,`heavy-atoms`,`rings`,`stereocentres`]),_h=new Set([`clogp`,`tpsa`,`hbd`,`hba`,`rotatable-bonds`,`lipinski-violations`,`qed`]);function vh(e,t,n,r,i){if(r||!e.length)return ch;let a=ph(i),o=new Set(e),s=[],c=(e,t,n={})=>({id:e,...a[e],value:t,copy:t,...n}),l,u=()=>l===void 0?l=t.mass():l,d,f=()=>d===void 0?d=t.descriptors():d,p=e.some(e=>hh.has(e)||e===`fragments`||e===`charge`||e===`smiles`||e===`cxsmiles`||e===`inchikey`||gh.has(e)||e===`bonds`||_h.has(e))?!!u()?.query:!1,m=!1,h=e=>{if(m)return;let n=t.text(`smarts`);n&&(m=!0,s.push(c(`smarts`,n,{title:e??(n.length>24?n:void 0)})))};for(let e of Qm)if(o.has(e)){if(hh.has(e)){if(p){h(i.infoQuery);continue}let t=u();if(!t||!t.formula)continue;if(e===`formula`){let n=zm(t.formula),r=lh(n);s.push(c(e,r,{runs:n}))}else s.push(c(e,e===`molecular-weight`?Bm(t.averageMass):Vm(t.exactMass)));continue}switch(e){case`smarts`:h();break;case`smiles`:case`cxsmiles`:case`inchikey`:{if(p)break;let n=t.text(e);n===null?s.push(c(e,`…`,{copy:``,pending:!0})):n&&s.push(c(e,n,n.length>24?{title:n}:{}));break}case`bonds`:p||s.push(c(e,String(t.bonds())));break;case`fragments`:case`charge`:{let t=u();if(p||!t||!t.fragments.length)break;e===`fragments`?s.push(c(e,String(t.fragments.length))):s.push(c(e,uh(t.fragments.reduce((e,t)=>e+t.charge,0))));break}case`atoms`:case`heavy-atoms`:case`rings`:case`stereocentres`:{let t=p?null:f();if(t===`pending`){s.push(c(e,`…`,{copy:``,pending:!0}));break}if(!t||t.query)break;let n=e===`atoms`?t.NumAtoms:e===`heavy-atoms`?t.NumHeavyAtoms:e===`rings`?t.NumRings:t.NumAtomStereoCenters;if(n===null)break;let r=e===`stereocentres`?t.NumUnspecifiedAtomStereoCenters??0:0;s.push(c(e,r>0?i.format(`readoutUnspecified`,{count:n,unspecified:r}):String(n)));break}case`clogp`:case`tpsa`:case`hbd`:case`hba`:case`rotatable-bonds`:case`lipinski-violations`:case`qed`:{let t=p?null:f();if(t===`pending`){s.push(c(e,`…`,{copy:``,pending:!0}));break}if(!t||t.query)break;let n=e===`clogp`?t.CrippenClogP:e===`tpsa`?t.tpsa:e===`hbd`?t.NumHBD:e===`hba`?t.NumHBA:e===`rotatable-bonds`?t.NumRotatableBonds:e===`lipinski-violations`?sh(t):t.qed;if(n===null)break;s.push(c(e,e===`clogp`||e===`tpsa`||e===`qed`?dh(n,2):String(n)));break}case`under-pointer`:case`drawing`:{let n=t.pointer();if(!n)break;e===`under-pointer`&&(n.kind===`atom`||n.kind===`bond`)&&s.push(c(e,yh(n,i),{pointer:!0})),e===`drawing`&&(n.kind===`draw`||n.kind===`move`)&&s.push(c(e,yh(n,i),{pointer:!0}));break}}}let g=l??null;if(g&&!g.query&&g.fragments.some(e=>e.rSites>0)){let e=s.map(e=>e.id).lastIndexOf(o.has(`exact-mass`)?`exact-mass`:`molecular-weight`);e>=0&&(s[e]={...s[e],note:i.infoRNotCounted})}return s.length?Object.freeze({kind:`items`,selection:n,items:Object.freeze(s)}):ch}function yh(e,t){if(e.kind===`atom`){let n=t.format(`readoutAtom`,{n:e.number,element:e.element,charge:uh(e.charge),h:e.h});return e.valence===null?n:t.format(`readoutAtomValence`,{atom:n,valence:e.valence})}return e.kind===`bond`?t.format(`readoutBond`,{n:e.number,type:e.type,length:dh(e.length,2)}):t.format(e.kind===`draw`?`readoutDrawn`:`readoutMoved`,{length:dh(e.length,2),angle:fh(e.angle)})}function bh(e,t){if(e.kind!==`items`||t.kind!==`items`||e.selection!==t.selection||!e.items.some(e=>e.pending))return e;let n=e.items.map(e=>{let n=e.pending?t.items.find(t=>t.id===e.id&&t.value!==`…`):void 0;return n?{...e,value:n.value,title:n.title,note:n.note}:e});return Object.freeze({...e,items:Object.freeze(n)})}var xh=e=>[e.label,e.value,e.note??``].filter(Boolean).join(` `);function Sh(e,t){if(e.kind===`empty`)return``;let n=e.items.filter(e=>!(e.pending&&e.value===`…`));return[...e.selection?[t.infoSelection]:[],...n.map(xh)].join(` `)}function Ch(e,t,n,r){if(e.reduce((e,t)=>e+t,0)<=n)return e.map((e,t)=>t);let i=n-r,a=new Set,o=[...e.keys()].sort((e,n)=>Number(!!t[n])-Number(!!t[e])||e-n);for(let n of o){if(e[n]>i){if(!t[n])break;continue}a.add(n),i-=e[n]}return[...a].sort((e,t)=>e-t)}var wh=Object.freeze({atomColours:`element`,warningMarks:`shown`,analysisPanel:`shown`,stereoLabels:`hidden`,stereoLabelStyle:`standard`,stereoFlag:`shown`,unspecifiedDoubleBonds:`crossed`,methylLabels:`hidden`,atomNumbers:`hidden`,keymap:`chemdraw`,readoutItems:$m,stereoPanel:`closed`}),Th=Object.freeze([`atomColours`,`warningMarks`,`analysisPanel`,`stereoLabels`,`stereoLabelStyle`,`stereoFlag`,`unspecifiedDoubleBonds`,`methylLabels`,`atomNumbers`,`keymap`,`readoutItems`,`stereoPanel`]),Eh=e=>typeof e==`string`?e.trim().toLowerCase():``,Dh=(e,t,n)=>Eh(e)===t?t:n;function Oh(e){let t=typeof e==`object`&&e?e:{};return Object.freeze({atomColours:a(t.atomColours),warningMarks:tp(t.warningMarks),analysisPanel:Nm(t.analysisPanel),stereoLabels:Dh(t.stereoLabels,`shown`,`hidden`),stereoLabelStyle:Dh(t.stereoLabelStyle,`registration`,`standard`),stereoFlag:Dh(t.stereoFlag,`hidden`,`shown`),unspecifiedDoubleBonds:Dh(t.unspecifiedDoubleBonds,`plain`,`crossed`),methylLabels:Eh(t.methylLabels)===`all`?`all`:Dh(t.methylLabels,`shown`,`hidden`),atomNumbers:Dh(t.atomNumbers,`shown`,`hidden`),keymap:$o(t.keymap),readoutItems:nh(t.readoutItems),stereoPanel:Sd(t.stereoPanel)})}var kh=(e,t)=>Th.every(n=>Ah(n,e[n],t[n]));function Ah(e,t,n){return e===`readoutItems`?Array.isArray(t)&&Array.isArray(n)&&rh(t,n):t===n}var jh=Object.freeze(Th.filter(e=>e!==`stereoPanel`)),Mh=`crux-sketch:settings`;function Nh(e){if(typeof e!=`object`||!e)return{};let{version:t,settings:n}=e;if(t!==1||typeof n!=`object`||!n||Array.isArray(n))return{};let r=n,i=Oh(r),a={};for(let e of jh)Object.prototype.hasOwnProperty.call(r,e)&&(a[e]=i[e]);return a}var Ph=e=>({version:1,settings:e});function Fh(e,t){let n={};for(let r of jh)Ah(r,e[r],t[r])||(n[r]=t[r]);return n}function Ih(e,t,n){let r={...e};for(let e of jh)t[e]!==void 0&&!n.has(e)&&(r[e]=t[e]);return Oh(r)}function Lh(e){try{return e?.localStorage??null}catch{return null}}function Rh(e){try{let t=Lh(e)?.getItem(Mh);return t?Nh(JSON.parse(t)):{}}catch{return{}}}function zh(e,t){if(Object.keys(t).length)try{let n=Lh(e);if(!n)return;let r=n.getItem(Mh),i={};try{i=r?Nh(JSON.parse(r)):{}}catch{}n.setItem(Mh,JSON.stringify(Ph({...i,...t})))}catch{}}var Bh=e=>!(e===!1||typeof e==`string`&&e.trim().toLowerCase()===`false`),Vh=`compact`,Hh=Object.freeze({width:720,height:600});function Uh(e){let t=typeof e==`string`?e.trim().toLowerCase():``;return t===`regular`||t===`compact`||t===`auto`?t:Vh}function Wh(e,t){return e===`auto`?!t||!(t.width>0)||!(t.height>0)?`regular`:t.width<Hh.width||t.height<Hh.height?`compact`:`regular`:e}var Gh=`dialog.settings`,Kh=[{name:`atomColours`,segment:`atom-colours`,label:`settingAtomColours`,check:[`element`,`black`]},{name:`warningMarks`,segment:`warning-marks`,label:`settingWarningMarks`,check:[`shown`,`hidden`]},{name:`analysisPanel`,segment:`readout`,label:`settingReadout`,options:[[`shown`,`readoutShown`],[`collapsed`,`readoutCollapsed`],[`hidden`,`readoutHidden`]]},{name:`stereoLabels`,segment:`stereo-labels`,label:`settingStereoLabels`,check:[`shown`,`hidden`]},{name:`stereoLabelStyle`,segment:`stereo-label-style`,label:`settingStereoLabelStyle`,options:[[`standard`,`stereoLabelStyleStandard`],[`registration`,`stereoLabelStyleRegistration`]]},{name:`stereoFlag`,segment:`stereo-flag`,label:`settingStereoFlag`,check:[`shown`,`hidden`]},{name:`unspecifiedDoubleBonds`,segment:`unspecified-double-bonds`,label:`settingUnspecified`,check:[`crossed`,`plain`]},{name:`methylLabels`,segment:`methyl-labels`,label:`settingMethyl`,options:[[`hidden`,`carbonLabelsNone`],[`shown`,`carbonLabelsTerminal`],[`all`,`carbonLabelsAll`]]},{name:`atomNumbers`,segment:`atom-numbers`,label:`settingAtomNumbers`,check:[`shown`,`hidden`]},{name:`keymap`,segment:`keymap`,label:`settingKeymap`,options:[[`chemdraw`,`keymapChemDraw`],[`ketcher`,`keymapKetcher`]]},{name:`readoutItems`,segment:`readout-items`,label:`readoutItemsHeading`,items:!0}],qh=[{id:`drawing`,heading:`settingsDrawing`,fields:[`atomColours`,`warningMarks`,`methylLabels`,`atomNumbers`]},{id:`stereo`,heading:`stereochemistry`,fields:[`stereoLabels`,`stereoLabelStyle`,`stereoFlag`,`unspecifiedDoubleBonds`]},{id:`keyboard`,heading:`settingsKeyboard`,fields:[`keymap`]},{id:`readout`,heading:`settingsReadout`,fields:[`analysisPanel`,`readoutItems`]}],Jh=qh.flatMap(e=>e.fields.map(e=>Kh.find(t=>t.name===e)));if([...Kh.map(e=>e.name)].sort().join()!==Th.filter(e=>e!==`stereoPanel`).slice().sort().join()||Jh.some(e=>!e)||[...Jh.map(e=>e.name)].sort().join()!==[...Kh.map(e=>e.name)].sort().join())throw Error(`Crux Sketch: the Settings dialog and the settings model disagree`);var Yh=class{#e;#t;#n=null;#r=wh;constructor(e,t){this.#e=e,this.#t=t}get open(){return this.#n!==null}show(e){if(this.#n!==null||!this.#t.isConnected)return;this.#n=e??{kind:`canvas`},this.#r=this.#e.settings(),this.#h();let t=this.#u();t&&(t.showModal(),this.#t.querySelector(`[data-testid="${Gh}.${Jh[0].segment}"]`)?.focus())}cancel(){let e=this.#n;e!==null&&(this.#l(),this.#e.restoreFocus(e))}relabel(){this.#n!==null&&this.#h()}#i(){let e=this.#n;if(e===null)return;let t=this.#r;this.#l(),this.#e.apply(t),this.#e.restoreFocus(e)}#a(){let e=this.#o();kh(this.#r,e)||(this.#r=e,this.#h())}#o(){let e=new Set(Kh.map(e=>e.name)),t=Object.fromEntries(Th.filter(t=>!e.has(t)).map(e=>[e,this.#r[e]]));return Oh({...wh,...t})}#s(e,t){this.#r=Oh({...this.#r,readoutItems:ih(this.#r.readoutItems,e,t)}),this.#h()}#c(e,t){this.#r=Oh({...this.#r,[e]:t}),this.#h()}#l(){this.#u()?.close(),this.#n=null,this.#h()}#u(){return this.#t.querySelector(`dialog[data-testid="${Gh}"]`)}#d(e){if(e.key===`Escape`){e.preventDefault(),this.cancel();return}if(e.key===`Tab`){kr(e);return}let t=e.composedPath()[0];e.key===`Enter`&&t instanceof HTMLInputElement&&t.type===`checkbox`&&(e.preventDefault(),this.#i())}#f(e,t){let n=`st-${e.segment}`,r=`${Gh}.${e.segment}`;if(e.items)return this.#p(e,t);let i=this.#r[e.name];if(e.check){let[a,o]=e.check;return x`<div class="check-row">
        <input
          type="checkbox"
          id=${n}
          data-testid=${r}
          .checked=${i===a}
          @change=${t=>this.#c(e.name,t.target.checked?a:o)}
        /><label for=${n}>${t[e.label]}</label>
      </div>`}return x`<div class="select-row">
      <label for=${n}>${t[e.label]}</label>
      <select id=${n} data-testid=${r} .value=${i} @change=${t=>this.#c(e.name,t.target.value)}>
        ${(e.options??[]).map(([e,n])=>x`<option value=${e} ?selected=${e===i}>${t[n]}</option>`)}
      </select>
    </div>`}#p(e,t){let n=ph(t),r=new Set(this.#r.readoutItems);return x`<div class="readout-items" role="group" aria-labelledby="st-${e.segment}">
      <div class="readout-items-heading" id="st-${e.segment}">${t[e.label]}</div>
      ${Zm.map(i=>x`<div class="readout-group" role="group" aria-labelledby="st-${e.segment}-${i.id}">
          <div class="readout-group-name" id="st-${e.segment}-${i.id}">${mh(i.id,t)}</div>
          <div class="readout-group-items">
            ${i.items.map(t=>{let i=`st-${e.segment}-${t}`;return x`<div class="check-row">
                <input
                  type="checkbox"
                  id=${i}
                  data-testid="${Gh}.${e.segment}.${t}"
                  .checked=${r.has(t)}
                  @change=${e=>this.#s(t,e.target.checked)}
                /><label for=${i}>${n[t].name}</label>
              </div>`})}
          </div>
        </div>`)}
    </div>`}#m(e,t){let n=`st-sec-${e.id}`;return x`<div class="st-section st-${e.id}" role="group" aria-labelledby=${n} data-section=${e.id}>
      <h3 class="st-heading" id=${n}>${t[e.heading]}</h3>
      ${e.fields.map(e=>this.#f(Kh.find(t=>t.name===e),t))}
    </div>`}#h(){let e=this.#e.strings(),t=kh(this.#r,this.#o());T(this.#n===null?w:x`<dialog
            class="dialog settings"
            data-testid=${Gh}
            aria-labelledby="st-title"
            @keydown=${e=>this.#d(e)}
            @cancel=${e=>{e.preventDefault(),this.cancel()}}
          >
            <div class="dialog-head">
              <h2 id="st-title">${e.settings}</h2>
              <button type="button" class="close" data-testid="${Gh}.close" aria-label=${e.close} @click=${()=>this.cancel()}>×</button>
            </div>
            <div class="settings-fields">
              <div class="st-column">${qh.filter(e=>e.id!==`readout`).map(t=>this.#m(t,e))}</div>
              <div class="st-column">${qh.filter(e=>e.id===`readout`).map(t=>this.#m(t,e))}</div>
            </div>
            <div class="dialog-foot">
              <button type="button" class="reset" data-testid="${Gh}.reset" aria-disabled=${t?`true`:w} @click=${()=>this.#a()}>${e.reset}</button>
              <button type="button" class="primary" data-testid="${Gh}.apply" @click=${()=>this.#i()}>${e.apply}</button>
              <button type="button" data-testid="${Gh}.cancel" @click=${()=>this.cancel()}>${e.cancel}</button>
            </div>
          </dialog>`,this.#t)}},Xh=`dialog.keys`,Zh=class{#e;#t;#n=null;constructor(e,t){this.#e=e,this.#t=t}get open(){return this.#n!==null}show(e){if(this.#n!==null||!this.#t.isConnected)return;this.#n=e??{kind:`canvas`},this.#i();let t=this.#t.querySelector(`dialog[data-testid="${Xh}"]`);t&&(t.showModal(),this.#t.querySelector(`[data-testid="${Xh}.close"]`)?.focus())}cancel(){let e=this.#n;e!==null&&(this.#t.querySelector(`dialog[data-testid="${Xh}"]`)?.close(),this.#n=null,this.#i(),this.#e.restoreFocus(e))}refresh(){this.#n!==null&&this.#i()}#r(e){e.key===`Escape`?(e.preventDefault(),this.cancel()):e.key===`Tab`&&kr(e)}#i(){let e=this.#e.strings();T(this.#n===null?w:x`<dialog
            class="dialog keys"
            data-testid=${Xh}
            aria-labelledby="ks-title"
            @keydown=${e=>this.#r(e)}
            @cancel=${e=>{e.preventDefault(),this.cancel()}}
          >
            <div class="dialog-head">
              <h2 id="ks-title">${e.keyboardShortcuts}</h2>
              <button type="button" class="close" data-testid="${Xh}.close" aria-label=${e.close} @click=${()=>this.cancel()}>×</button>
            </div>
            <div class="sheet" tabindex="0" role="region" aria-labelledby="ks-title" data-testid="${Xh}.sheet">
              ${this.#e.sheet().map(t=>x`<table data-testid="${Xh}.${t.scope}">
                  <caption>${t.title}</caption>
                  <tbody>
                    <tr>
                      <th scope="col">${e.sheetWhat}</th>
                      <th scope="col">${e.sheetKeys}</th>
                    </tr>
                    ${t.rows.map(t=>x`<tr>
                        <th scope="row">${t.what}</th>
                        <td>${t.keys.map((t,n)=>x`${n?x`<span class="or">${e.keysOr}</span>`:w}<kbd>${t}</kbd>`)}</td>
                      </tr>`)}
                  </tbody>
                </table>`)}
            </div>
          </dialog>`,this.#t)}},Qh=`dialog.about`,$h=class{#e;#t;#n=null;constructor(e,t){this.#e=e,this.#t=t}get open(){return this.#n!==null}show(e){if(this.#n!==null||!this.#t.isConnected)return;this.#n=e??{kind:`canvas`},this.#i();let t=this.#t.querySelector(`dialog[data-testid="${Qh}"]`);t&&(t.showModal(),this.#t.querySelector(`[data-testid="${Qh}.close"]`)?.focus())}cancel(){let e=this.#n;e!==null&&(this.#t.querySelector(`dialog[data-testid="${Qh}"]`)?.close(),this.#n=null,this.#i(),this.#e.restoreFocus(e))}refresh(){this.#n!==null&&this.#i()}#r(e){e.key===`Escape`?(e.preventDefault(),this.cancel()):e.key===`Tab`&&kr(e)}#i(){let e=this.#e.strings(),t=this.#e.engineRev();T(this.#n===null?w:x`<dialog
            class="dialog about"
            data-testid=${Qh}
            aria-labelledby="ab-title"
            @keydown=${e=>this.#r(e)}
            @cancel=${e=>{e.preventDefault(),this.cancel()}}
          >
            <div class="dialog-head">
              <h2 id="ab-title">${e.aboutTitle}</h2>
              <button type="button" class="close" data-testid="${Qh}.close" aria-label=${e.close} @click=${()=>this.cancel()}>×</button>
            </div>
            <dl class="facts">
              <dt>${e.aboutVersion}</dt>
              <dd data-testid="${Qh}.version">${this.#e.version()}</dd>
              <dt>${e.aboutEngine}</dt>
              <dd data-testid="${Qh}.engine"><code>${t||e.aboutEngineLoading}</code></dd>
            </dl>
            <p class="notices" data-testid="${Qh}.notices">${e.aboutNotices}</p>
          </dialog>`,this.#t)}},eg=16384,tg=4096;function ng(e){let t=rg(e);return t?Math.min(2,tg/Math.max(t.width,t.height)):2}function rg(e){let t=/<svg\b[^>]*>/.exec(e)?.[0];if(!t)return null;let n=e=>{let n=RegExp(`\\s${e}=['"]([0-9.]+)(?:px)?['"]`).exec(t);return n?Number(n[1]):NaN},[r,i]=[n(`width`),n(`height`)];return Number.isFinite(r)&&Number.isFinite(i)&&r>0&&i>0?{width:r,height:i}:null}function ig(e,t=1){if(typeof t!=`number`||!Number.isFinite(t)||t<=0)throw RangeError(`png: scale must be a finite number above 0, not ${String(t)}`);let n=Math.min(t,eg/Math.max(e.width,e.height));return{width:Math.max(1,Math.ceil(e.width*n)),height:Math.max(1,Math.ceil(e.height*n)),scale:n}}function ag(e,t,n){return e.replace(/<svg\b[^>]*>/,e=>{let r=e.replace(/\swidth=(['"])[^'"]*\1/,` width='${t}px'`).replace(/\sheight=(['"])[^'"]*\1/,` height='${n}px'`);return/\swidth=/.test(r)||(r=r.replace(/^<svg/,`<svg width='${t}px'`)),/\sheight=/.test(r)||(r=r.replace(/^<svg/,`<svg height='${n}px'`)),/\sxmlns=/.test(r)||(r=r.replace(/^<svg/,`<svg xmlns='http://www.w3.org/2000/svg'`)),r})}var og=(()=>{let e=new Uint32Array(256);for(let t=0;t<256;t++){let n=t;for(let e=0;e<8;e++)n=n&1?3988292384^n>>>1:n>>>1;e[t]=n>>>0}return e})();function sg(e){let t=4294967295;for(let n of e)t=og[(t^n)&255]^t>>>8;return(t^4294967295)>>>0}function cg(e,t){if(e.length<33||[137,80,78,71,13,10,26,10].some((t,n)=>e[n]!==t))return e;let n=new DataView(e.buffer,e.byteOffset,e.byteLength),r=t=>String.fromCharCode(e[t+4],e[t+5],e[t+6],e[t+7]);if(r(8)!==`IHDR`)return e;let i=20+n.getUint32(8),a=[];for(let t=i;t+12<=e.length;){let i=n.getUint32(t),o=t+12+i;r(t)!==`pHYs`&&a.push(e.subarray(t,Math.min(o,e.length))),t=o}let o=Math.round(t/.0254),s=new Uint8Array(21),c=new DataView(s.buffer);c.setUint32(0,9),s.set([112,72,89,115],4),c.setUint32(8,o),c.setUint32(12,o),s[16]=1,c.setUint32(17,sg(s.subarray(4,17)));let l=new Uint8Array(i+s.length+a.reduce((e,t)=>e+t.length,0));l.set(e.subarray(0,i),0),l.set(s,i);let u=i+s.length;for(let e of a)l.set(e,u),u+=e.length;return l}async function lg(e,t,n=1){let r=rg(e);if(!r)return ig({width:1,height:1},n),null;let i=ig(r,n),a=t.defaultView;if(!a)throw Error(`png: the sketcher is not in a window`);let o=new a.Image(i.width,i.height);await new Promise((t,n)=>{o.onload=()=>t(),o.onerror=()=>n(Error(`png: the browser could not draw the SVG`)),o.src=`data:image/svg+xml;charset=utf-8,${encodeURIComponent(ag(e,i.width,i.height))}`});let s=t.createElement(`canvas`);s.width=i.width,s.height=i.height;let c=s.getContext(`2d`);if(!c)throw Error(`png: the browser gave no canvas to draw on`);c.drawImage(o,0,0,i.width,i.height);let l=await new Promise(e=>s.toBlob(e,`image/png`));if(s.width=s.height=0,!l)throw Error(`png: the browser could not write the PNG`);let u=cg(new Uint8Array(await l.arrayBuffer()),96*i.scale);return new a.Blob([u],{type:`image/png`})}var ug=class{#e;#t;#n=null;#r=null;#i=null;#a=null;#o=!1;#s=null;constructor(e,t){this.#e=e,this.#t=t}get busy(){return this.#i!==null||this.#a!==null}get unsettled(){return this.#o||this.#n!==null}step(e,t,n){let r=this.#u();r&&this.#c(r,je(r.zoom,e),t,n)}zoomTo(e,t,n){let r=this.#u();r&&this.#c(r,e,t,n)}#c(e,t,n,r){let i=this.#e.size();if(!i)return;this.#f();let a=r?Ue(e,r,t):u(e,i,n??y(i),t);dg(a,e)||(this.#e.show(a),this.#l())}panBy(e,t){let n=this.#u();n&&(e||t)&&(this.#f(),this.#e.show(Je(n,e,t)),this.#l())}wheel(e,t){let n=this.#u(),r=this.#e.size();n&&r&&this.#d(ni(n,r,t,e))}grab(e){this.#a=null,this.#i={...e}}drag(e){let t=this.#i,n=this.#u();t&&n&&(this.#i={...e},(e.x!==t.x||e.y!==t.y)&&this.#d(Je(n,e.x-t.x,e.y-t.y)))}release(){this.#i=null,this.flush()}pinchStart(e,t){this.#i=null,this.#a=new Map([[e[0],{...e[1]}],[t[0],{...t[1]}]])}pinchMove(e,t){let n=this.#a,r=n?.get(e),i=this.#u(),a=this.#e.size();if(!n||!r||!i||!a)return;let[o]=[...n.entries()].filter(([t])=>t!==e).map(([,e])=>e);o&&(n.set(e,{...t}),this.#d(se(i,a,[r,o],[t,o])))}pinchEnd(){this.#a=null,this.flush()}flush(){let e=this.#n,t=this.#o;this.#f(),this.#l();let n=e??(t?this.#e.view():null);n&&this.#e.show(n)}drawn(){this.#o&&this.#l()}cancel(){this.#i=null,this.#a=null,this.#f(),this.#l()}#l(){this.#o=!1,this.#s?.(),this.#s=null}#u(){return this.#n??this.#e.view()}#d(e){this.#n=e,this.#r??=this.#t.frame(()=>{this.#r=null;let e=this.#n;if(this.#n=null,e){if(!this.#e.move?.(e)){this.#e.show(e);return}this.#o=!0,this.#s?.(),this.#s=this.busy?null:this.#t.after(100,()=>{this.#s=null,this.flush()})}})}#f(){this.#n=null,this.#r?.(),this.#r=null}},dg=(e,t)=>e.zoom===t.zoom&&e.cx===t.cx&&e.cy===t.cy&&e.docLen===t.docLen,fg=Object.freeze({slots:Object.freeze([]),age:Object.freeze([])});function pg(e,t,n,r=7){if(n.has(t)||e.slots.includes(t)||r<1)return e;if(e.slots.length<r)return{slots:[...e.slots,t],age:[...e.age,t]};let[i,...a]=e.age;return{slots:e.slots.map(e=>e===i?t:e),age:[...a,t]}}var mg=class{format;errors;#e;#t=null;#n=null;constructor(e,t,n){this.format=e,this.#e=t??null,this.errors=n}get available(){return this.#e!==null}loaded(){return this.#e?.loaded()??null}get loading(){return this.#n!==null}get failed(){return this.#t}load(){let e=this.#e;return e?(this.#n??=e.load().then(e=>(this.#n=null,this.#t=null,e),e=>{throw this.#n=null,this.#t={cause:e},e}),this.#n):Promise.reject(Error(`no ${this.format} module`))}notLoaded(){return this.errors.notLoaded(this.#t?.cause??Error(this.#e?`not loaded`:`no ${this.format} module`))}forGetter(){let e=this.loaded();if(e)return e;throw this.#t&&!this.#n?this.errors.notLoaded(this.#t.cause):(this.#n||this.load().catch(()=>{}),this.errors.loading)}},hg=class{#e=new Map;on(e,t){let n=this.#e.get(e);return n||this.#e.set(e,n=new Set),n.add(t),()=>this.off(e,t)}off(e,t){this.#e.get(e)?.delete(t)}emit(e,t){let n=this.#e.get(e);if(n?.size)for(let e of[...n])try{e(t)}catch(e){queueMicrotask(()=>{throw e})}}clear(){this.#e.clear()}},gg=`Crux Sketch: this editor was destroyed`,_g=`engine`,vg=e=>Number.isInteger(e),yg=e=>e instanceof Error?`${e.name}: ${e.message}`:String(e);function bg(e){if(e&&typeof e==`object`&&!(e instanceof Error)){let t=e;if(typeof t.format==`string`&&typeof t.message==`string`){let e={format:t.format,message:t.message};return vg(t.position)&&(e.position=t.position),vg(t.line)&&(e.line=t.line),Array.isArray(t.atoms)&&(e.atoms=Object.freeze(t.atoms.filter(vg))),Object.freeze(e)}}return Object.freeze({format:`unknown`,message:e instanceof Error?e.message:String(e)})}function xg(e){return typeof WebAssembly<`u`&&e instanceof WebAssembly.RuntimeError}function Sg(e){return Object.freeze({format:_g,message:`Crux Sketch: the engine failed (${yg(e)}); reload the page to use the sketcher again`})}function Cg(e){return Object.freeze({format:_g,message:`Crux Sketch: the engine could not be loaded (${yg(e)})`})}function wg(e,t){return Object.freeze({format:e,message:`the engine could not write this drawing as ${e} (${yg(t)})`})}function Tg(e){return Object.freeze({format:`unknown`,message:`the engine could not copy this molecule to edit it (${yg(e)})`})}function Eg(){return Object.freeze({format:`unknown`,message:p.aliasLimit})}function Dg(e,t){return Object.freeze({format:e,message:`the engine could not draw this molecule (${yg(t)})`})}var Og=class{#e=[];#t=[];#n=0;#r=0;#i;budget;constructor(e,t=1/0){this.#i=e,this.budget=t}get current(){return this.#e[this.#r]??null}get index(){return this.#r}get size(){return this.#e.length}get weight(){return this.#n}get canUndo(){return this.#r>0}get canRedo(){return this.#r+1<this.#e.length}at(e){return this.#e[e]??null}restart(e,t=0){this.free(),e!==null&&(this.#e.push(e),this.#t.push(t),this.#n=t)}push(e,t=0){let n=this.#e.splice(this.#r+1);for(let e of this.#t.splice(this.#r+1))this.#n-=e;this.#e.push(e),this.#t.push(t),this.#n+=t;let r=0;for(;this.#n>this.budget&&r<this.#e.length-1;)this.#n-=this.#t[r++];let i=this.#e.splice(0,r);this.#t.splice(0,r),this.#r=this.#e.length-1;for(let e of n)this.#i(e);for(let e of i)this.#i(e)}pop(){if(this.#r===0||this.#r!==this.#e.length-1)return!1;let e=this.#e.pop();return this.#n-=this.#t.pop(),--this.#r,this.#i(e),!0}moveTo(e){if(!Number.isInteger(e)||e<0||e>=this.#e.length)throw RangeError(`the history has no state ${e} (it has ${this.#e.length})`);return this.#r=e,this.#e[e]}free(){let e=this.#e;this.#e=[],this.#t=[],this.#n=0,this.#r=0;for(let t of e)this.#i(t)}};function kg(e){try{e.free()}catch{}}function Ag(e,t){if(e.setBondTypes)e.setBondTypes(Uint32Array.from(t.map(([e])=>e)),t.map(([,e])=>e));else for(let[n,r]of t)e.setBondType(n,r)}function jg(e,t){if(e.isEdited()||!e.hasQuery())return t(e);let n=un(e);if(!n.length)return t(e);let r=e.clone();try{return Ag(r,n),t(r)}finally{kg(r)}}function Mg(e,t){let n=JSON.parse(e);return(t===`atom`?n.kind:n.type)===`custom`&&typeof n.smarts==`string`?n.smarts:null}function Ng(e,t){if(!t.isEdited()||!t.hasQuery())return t.toSmarts();let n;try{n=e.parse(t.toMolblockV3000(),`molV3000`,!0)}catch(e){if(xg(e))throw e;return t.toSmarts()}try{if(n.numAtoms()!==t.numAtoms()||n.numBonds()!==t.numBonds())return t.toSmarts();for(let e=0;e<t.numAtoms();e++){let r=Mg(t.atomQueryJson(e),`atom`);if(r!==null&&!n.setAtomQuerySmarts(e,r).ok)return t.toSmarts()}for(let e=0;e<t.numBonds();e++){let r=Mg(t.bondQueryJson(e),`bond`);if(r!==null&&!n.setBondQuerySmarts(e,r).ok)return t.toSmarts()}return n.toSmarts()}finally{kg(n)}}function Pg(e,t,n){return t===e?{format:e.includes(` |`)?`cxsmiles`:`smiles`,text:e}:n(e)===t?{format:`smiles`,text:e}:{format:`cxsmiles`,text:t}}function Fg(e,t){let n=null;try{return e.cxsmilesOf?e.cxsmilesOf(t):(n=e.parse(t,`smiles`,!0),n.toCxsmiles())}catch(e){if(xg(e))throw e;return null}finally{try{n?.free()}catch{}}}var Ig=`0.1.0`,Lg=[`smiles`,`molfile`,`molV3000`,`cxsmiles`,`smarts`],Rg={smiles:`smiles`,plainSmiles:`smiles`,molfile:`molV2000`,molV3000:`molV3000`,cxsmiles:`cxsmiles`,cxsmilesXY:`cxsmiles`,smarts:`smarts`,canonicalSmiles:`smiles`,canonicalPlain:`smiles`,canonicalCxsmiles:`cxsmiles`},zg=[`smiles`,`molV2000`,`molV3000`,`cxsmiles`,`smarts`],Bg={plainSmiles:e=>e.toSmiles(),molfile:e=>e.toMolblock(),molV3000:e=>e.toMolblockV3000(),cxsmiles:e=>e.toCxsmiles(),cxsmilesXY:e=>e.toCxsmiles(!0),smarts:e=>e.toSmarts(),canonicalPlain:e=>e.toSmilesCanonical(),canonicalCxsmiles:e=>e.toCxsmilesCanonical()},Vg=e=>{let t=e.indexOf(` |`);return t>=0?e.slice(0,t):e},Hg=new Set([`plainSmiles`,`cxsmiles`,`cxsmilesXY`,`canonicalPlain`,`canonicalCxsmiles`]),Ug=[`plainSmiles`,`molfile`,`molV3000`,`cxsmiles`,`smarts`],Wg=5e5,Gg=(e,t)=>`${e} ${t?[...new Set(t)].sort((e,t)=>e-t).join(` `):`*`}`,Kg=24,qg=e=>e.numAtoms()+e.numBonds(),Jg=!0,Yg=e=>e===`query`?`query`:`molecule`,Xg=new WeakMap;function Zg(e){return Object.freeze(e.offset===null?{format:`query`,message:e.reason}:{format:`smarts`,message:e.reason,position:e.offset})}function G(e){try{e?.free()}catch{}}function Qg(e,t){if(!t.numBonds())return[];let n=new Set(t.ringBonds().flat());if(!n.size)return[];let r=null;try{return r=e.parse(t.toMolblockV3000(),`molV3000`,!0),r.numAtoms()!==t.numAtoms()||r.numBonds()!==t.numBonds()?[]:Yn(r.geometryJson()).bonds.flatMap((e,t)=>e.type===`aromatic`&&n.has(t)?[t]:[])}catch(e){if(xg(e))throw e;return[]}finally{G(r)}}var $g=new FinalizationRegistry(e=>e.free()),e_=class{of=null;key=``;copy=null;free(){G(this.copy),this.copy=null,this.of=null,this.key=``}},t_=class{#e=`pending`;#t=null;#n=null;#r=new Og(G,Wg);#i=null;#a=null;#o=!0;#s;#c=new hg;#l;#u=``;#d=[];#f=null;#p=!0;#m=new Map;#h=null;#g=null;#_=null;#v={};#y=new Map;#b=null;#x=null;#S=``;#C=0;#w=null;#T=!1;#E;#D;#O=!1;#k=0;#A=new e_;#j=(e,t)=>{let n=this.#N(e,t)??this.#fe(e,t);return this.#W.push(...n.drops),this.#P=n.info,n.mol};#M=null;#N(e,t){let n=this.#M;return!n||n.text!==e||n.format!==t?null:(this.#M=null,n)}#P=null;#F=null;#I;#L;#R=[];#z=null;#B=null;#V=null;#H=``;#U=null;#W=[];#G=null;#K=Object.freeze([]);#q=!1;#J={kind:`write`,label:null};#Y=new WeakMap;get lastChange(){return this.#J}#X=new WeakMap;#Z=ge;constructor(e={}){this.#l=Yg(e.mode),this.#s=e.canvas??null,this.#E=e.interrupt??null,this.#D=e.hostCall??null,this.#I=new mg(`cdxml`,e.cdxml,{loading:Fa,notLoaded:Pa}),this.#L=new mg(`inchi`,e.inchi,{loading:Xa,notLoaded:Za})}get ready(){return this.#e===`ready`}get destroyed(){return this.#e===`destroyed`}start(e){if(this.#e!==`pending`)return;this.#t=e,this.#e=`ready`;let t=this.#d;this.#d=[],this.#u=``,this.#Me(()=>{this.#S=e.engineRev();let t=e.empty();this.#r.restart(t,qg(t)),this.#i=t,this.#a=t,this.#p=!0});let n=[];for(let e of t)this.#Le(e.text,e.format,!1)&&n.push(e);let r=this.#s;for(;;){let e=this.#r.current;if(!e||!r||this.#Ae())return;let t;try{t=r.draw(e,`fit`,this.#be(e,!1))}catch(e){let t=n.pop(),r=this.#je(e);if(this.#Ne(r??Dg(t?t.format??this.#Ve(t.text):`unknown`,e)),r||!t)return;this.#Be(n.at(-1));continue}this.#ze(t),r.show(t);return}}failToLoad(e){this.#e===`pending`&&(this.#e=`ready`,this.#d=[],this.#n=Cg(e),this.#Ne(this.#n))}releaseWith(e){this.#e!==`destroyed`&&($g.register(e,this.#r,this.#r),$g.register(e,this.#A,this.#A))}destroy(){if(this.#e!==`destroyed`){for(let e of[...Lg,`plainSmiles`,`canonicalSmiles`,`canonicalCxsmiles`,...this.#q?[`cxsmilesXY`]:[]])try{this.#Fe(e)}catch{}if(this.warnings,this.#re())try{this.#me()}catch{}else!this.#p&&this.#r.current&&this.#I.available&&(this.#B=Ia);if(this.#ie())try{this.#he()}catch{}else!this.#p&&this.#r.current&&this.#L.available&&(this.#V=Qa);this.#b=this.analysis,this.#x=this.hasQuery;try{this.#Q=this.#$()}catch(e){this.#Q=e}this.#C=this.renderCount,this.#o=this.isSaved,this.#e=`destroyed`,$g.unregister(this.#r),$g.unregister(this.#A),this.#r.free(),this.#A.free(),this.#d=[],this.#s?.destroy(),this.#s=null,this.#c.clear(),this.#E=null,this.#D=null}}get smiles(){return this.#Fe(`smiles`)}get plainSmiles(){let e=this.#Fe(`plainSmiles`);return this.hasQuery?Vg(e):e}get molfile(){return this.#Fe(`molfile`)}get molV3000(){return this.#Fe(`molV3000`)}get cxsmiles(){return this.#Fe(this.#q?`cxsmilesXY`:`cxsmiles`)}get canonicalSmiles(){return this.#Fe(`canonicalSmiles`)}get canonicalCxsmiles(){return this.#Fe(`canonicalCxsmiles`)}get cxsmilesCoordinates(){return this.#q}set cxsmilesCoordinates(e){this.#e!==`destroyed`&&(this.#q=Oa(e))}get readWarnings(){return this.#K}get smarts(){return this.#Fe(`smarts`)}get svg(){let e=this.#$();return e?sn(e,Mr()):``}#Q=``;#$(){if(this.#e===`destroyed`){if(typeof this.#Q==`string`)return this.#Q;throw this.#Q}let e=this.#r.current,t=this.#s;if(this.#e!==`ready`||!e||this.#p||!t?.exportSvg)return``;let n=this.#Ae();if(n)throw n;try{let n=t.exportSvg(e,this.#w?.docLen??e.meanBondLength(),{...this.#be(e,this.#O),text:!0});if(!n)throw Error(`it drew nothing`);return n}catch(e){throw this.#je(e)??wg(`svg`,e)}}get isEmpty(){return this.#p}get isSaved(){return this.#e===`destroyed`?this.#o:this.#r.current===this.#a}markSaved(){this.#e!==`destroyed`&&(this.#a=this.#r.current)}get lastError(){return this.#f}get noEzRefused(){return this.#G}noEz(e){let t=this.#r.current;return this.#e!==`ready`||!t||this.#Ae()||!e.length?!1:this.#Me(()=>{let n=t.clone();try{return e.every(e=>{try{return n.setBondStereo(e,`crossed`),!1}catch(e){if(nt(e))return!0;if(xg(e))throw e;return!1}})}finally{G(n)}})??!1}get warnings(){let e=!!this.#re(),t=!!this.#ie(),n=`${e?`cdxml`:``},${t?`inchi`:``},${this.#q?`xy`:``}`;if(this.#g&&(this.#H===n||this.#e===`destroyed`))return this.#g;let r=this.#r.current;if(this.#e!==`ready`||!r||this.#p||this.#Ae())return[];let i=this.problems.map(e=>Object.freeze({kind:e.kind,message:e.message,atoms:Object.freeze([...e.atoms])}));try{for(let e of zg){let t=e===`smiles`?this.#ne():e,n=e===`cxsmiles`&&this.#q?r.writeWarnings(t,!0):r.writeWarnings(t),a=e===`molV2000`?n.filter(e=>e.code===`v3000`):[];for(let t of a.length?a:n)i.push(Ea(e,t))}}catch(e){return this.#je(e),[]}if(e)try{i.push(...this.#me()?.warnings??[])}catch{}if(t)try{i.push(...this.#he()?.warnings??[])}catch{}return this.#g=Object.freeze(i),this.#H=n,this.#g}get problems(){if(this.#_)return this.#_;let e=this.#r.current;if(this.#e!==`ready`||!e||this.#p||this.#Ae())return[];try{let t=e.numAtoms(),n=Yn(e.geometryJson()),r=hp(n.atoms,n.bonds,this.#w?.docLen??zn(e.meanBondLength())),i=gp(n.atoms,n.bonds,t=>ip(me(e.atomJson(t)))),a=dp(this.#_e(e,!1),n.bonds,t),o=vp(e.problems(),t,t=>me(e.atomJson(t)),this.#ee(e,t),pp(e.sgroupsJson?.()),r,i,a);this.#_=Object.freeze(o.map(e=>Object.freeze({...e,atoms:Object.freeze([...e.atoms])})))}catch(e){return this.#je(e),[]}return this.#_}#ee(e,t){let n=this.#_e(e,!1).atoms;return n.some(e=>e.status===`unspecified`)?up(n,Yn(e.geometryJson()).bonds,t):[]}get analysis(){return this.#e===`destroyed`?this.#b:Lm(this.massInfo())}radicalDropped(e,t){let n=this.#r.current;if(this.#e!==`ready`||!n||this.#Ae())return null;let r=null;try{return r=n.clone(),C(r,e,this.#j)?xp(r,t,this.#j):null}catch(e){return xg(e)&&this.#je(e),null}finally{G(r)}}analysisOf(e){return e.length?Lm(this.massInfo(e)):null}massInfo(e){let t=this.#r.current;if(this.#e!==`ready`||!t||this.#Ae())return null;let n=e?[...new Set(e)].sort((e,t)=>e-t).join(` `):null;if(n===null&&this.#v.all!==void 0)return this.#v.all;if(n!==null&&this.#v.some?.key===n)return this.#v.some.info;let r;try{r=this.#p?Im(`{"query":false,"formula":"","averageMass":0,"exactMass":0,"fragments":[]}`):Im(e?t.massInfo([...e]):t.massInfo())}catch(e){this.#je(e),r=null}return n===null?this.#v.all=r:this.#v.some={key:n,info:r},r}#te(e,t,n,r){let i=this.#r.current;if(this.#e!==`ready`||!i||this.#Ae())return r;let a=Gg(e,t);if(this.#y.has(a))return this.#y.get(a);let o;try{o=n(i)}catch(e){this.#je(e),o=r}return this.#y.size>=Kg&&this.#y.clear(),this.#y.set(a,o),o}readoutKnown(e,t){return this.#p||this.#e!==`ready`||e===`inchikey`&&(!this.#L.available||!this.#ie())?!0:e===`inchikey`&&!t?this.#V!==null:this.#y.has(Gg(e,t))}descriptors(e){return this.#p?null:this.#te(`descriptors`,e,t=>typeof t.descriptors==`function`?oh(e?t.descriptors(Uint32Array.from(e)):t.descriptors()):null,null)}bondsAmong(e){return this.#p?0:this.#te(`bonds`,e,t=>{if(!e)return t.numBonds();let n=t.extract([...new Set(e)]);try{return n.numBonds()}finally{G(n)}},0)}readoutText(e,t){if(this.#p)return``;if(e===`inchikey`){if(!this.#L.available)return``;if(!this.#ie())try{this.#L.forGetter()}catch(e){return e===this.#L.errors.loading?null:``}}if(!t)try{return e===`smiles`?this.plainSmiles:e===`cxsmiles`?this.cxsmiles:e===`smarts`?this.smarts:this.inchiKey}catch{return``}return this.#te(e,t,n=>{let r=n.extract([...new Set(t)]);try{if(e===`smiles`)return r.toSmiles();if(e===`cxsmiles`)return r.toCxsmiles(this.cxsmilesCoordinates);if(e===`smarts`)return this.#Pe(`smarts`,r);let t=this.#ie();return t?oo(t,Bg.molfile(r)).key:``}finally{G(r)}},``)}#ne(){try{this.#Fe(`smiles`)}catch{}return this.#h??`smiles`}get mode(){return this.#l}set mode(e){this.#e!==`destroyed`&&(this.#l=Yg(e))}get hasQuery(){if(this.#e===`destroyed`)return this.#x??!1;if(this.#x!==null)return this.#x;let e=this.#r.current;if(this.#e!==`ready`||!e||this.#p||this.#Ae())return!1;try{this.#x=e.hasQuery()}catch(e){return this.#je(e),!1}return this.#x}hasQueryIn(e,t){let n=this.#r.current;if(this.#e!==`ready`||!n||!e.length&&!t.length||this.#Ae())return!1;try{return n.hasQueryIn(Uint32Array.from(e),Uint32Array.from(t))}catch(e){return this.#je(e),!1}}get version(){return Ig}get engineRev(){return this.#S}get renderCount(){return this.#s?.renderCount??this.#C}get value(){return this.#e===`pending`&&this.#d.length?this.#d[this.#d.length-1].text:this.#O?this.#Fe(`molfile`):this.#u}set value(e){this.setValue(e)}setValue(e,t={}){if(this.#e===`destroyed`)throw Error(gg);let n=e==null?``:String(e),r=t?.format||void 0,i=t?.view===`keep`;this.#e===`pending`?this.#d.push({text:n,format:r}):this.#R.length||this.#se(n,r)?this.#ce(n,r,i):(this.#D?.(),this.#Le(n,r,!0,i))}#re(){return this.#I.loaded()}#ie(){return this.#L.loaded()}#ae(e){return e===`cdxml`?this.#I:e===`inchi`?this.#L:null}#oe(e,t,n){return t===`cdxml`||t===`inchi`?t===`inchi`&&Ja(e)?null:t:t?null:ka(e)&&n(e)===`cdxml`?`cdxml`:Ka(e)&&n(e)===`inchi`?`inchi`:null}#se(e,t){let n=this.#oe(e,t,e=>this.#Ve(e)),r=n?this.#ae(n):null;return r&&r.available&&!r.loaded()?n:null}#ce(e,t,n){this.#R.push({text:e,format:t,keep:n,needs:this.#se(e,t),tried:!1}),!this.#le&&(this.#le=!0,this.#z=this.#ue())}#le=!1;async#ue(){try{for(;this.#R.length;){let e=this.#R[0],t=e.needs&&!e.tried?this.#ae(e.needs):null;if(t&&!t.loaded()){e.tried=!0;try{await t.load()}catch{}continue}if(this.#R.shift(),this.#e!==`ready`){this.#R=[];break}this.#D?.(),this.#Le(e.text,e.format,!0,e.keep)}}finally{this.#le=!1}}formatLoaded(e){let t=this.#ae(e);return!t||!!t.loaded()}async loadFormat(e){if(this.#e===`destroyed`)throw Error(gg);let t=this.#ae(e);if(!t||t.loaded()){this.#z&&await this.#z;return}if(!t.available)throw t.notLoaded();try{await t.load()}catch(e){throw this.#z&&await this.#z,t.errors.notLoaded(e)}this.#z&&await this.#z}async prepare(e){if(this.#e!==`pending`)return;let t=new Set;for(let n of this.#d){let r=this.#oe(n.text,n.format,t=>this.#de(e,t)),i=r?this.#ae(r):null;r&&i?.available&&!i.loaded()&&t.add(r)}await Promise.all([...t].map(e=>this.#ae(e).load().catch(()=>{})))}#de(e,t){try{return e.detectFormat(t)}catch{return`unknown`}}#fe(e,t,n=!1){let r=this.#t;if((!t||t===`inchi`)&&Ja(e))throw Ya;let i=this.#oe(e,t,e=>r.detectFormat(e));if(i===`cdxml`){let t=this.#re();if(!t)throw this.#I.notLoaded();let n=za(r,t,e,Jg);return{mol:n.mol,drops:n.drops,info:{format:`cdxml`,structures:n.structures,left:n.left,lost:Ga(n.drops,n.left)}}}if(i===`inchi`){let t=this.#ie();if(!t)throw this.#L.notLoaded();let n=ro(t,e,e=>this.#pe(e,`smiles`));return{mol:n.mol,drops:n.drops,info:{format:`inchi`,structures:0,left:[],lost:[]}}}let a=n?r.parse(e,t,Jg,void 0,!0):this.#pe(e,t),o=Da(a),s=Ga(o),c=0;if(s.length)try{let e=Yn(a.geometryJson());c=Ra(e.atoms.length,e.bonds)}catch(e){throw G(a),e}return{mol:a,drops:o,info:{format:t??`detected`,structures:c,left:[],lost:s}}}#pe(e,t){let n=this.#t,r=n.parse(e,t,Jg,wt),i;try{let{atoms:e,bonds:t}=Yn(r.geometryJson());i=Ee(e.length,t.map(e=>[e.a,e.b]))}catch(e){throw G(r),e}return i?r:(G(r),n.parse(e,t,Jg))}get lastRead(){return this.#F}get inspectError(){return this.#U}get cdxml(){return this.#me()?.text??``}#me(){let e=this.#B;if(e){if(`text`in e)return e;throw e}let t=this.#r.current;if(this.#e!==`ready`||!t||this.#p||!this.#I.available)return null;let n=this.#Ae();if(n)throw n;let r=this.#I.forGetter(),i;try{i=Ba(t)}catch(e){throw this.#B=this.#je(e)??wg(`cdxml`,e),this.#B}try{this.#B=Ha(r,i)}catch(e){throw this.#B=xg(e)?Na(e):wg(`cdxml`,e),this.#B}return this.#B}get inchi(){return this.#he()?.inchi??``}get inchiKey(){return this.#he()?.key??``}#he(){let e=this.#V;if(e){if(`inchi`in e)return e;throw e}let t=this.#r.current;if(this.#e!==`ready`||!t||this.#p||!this.#L.available)return null;let n=this.#Ae();if(n)throw n;let r=this.#L.forGetter(),i;try{i=this.#Fe(`molfile`)}catch(e){let t=e&&typeof e==`object`&&typeof e.message==`string`?e.message:String(e);throw this.#V=this.#Ae()??Object.freeze({format:`inchi`,message:`the molfile InChI is written from could not be written (${t})`}),this.#V}try{this.#V=oo(r,i)}catch(e){throw this.#V=bg(e),this.#V}return this.#V}resize(e){let t=this.#r.current,n=this.#s;if(this.#e!==`ready`||!t||!n||this.#Ae())return;let r=this.#w&&!this.#T?this.#w:`fit`,i=this.#Me(()=>n.redraw(t,e,r,this.#be(t,this.#O)));i&&(this.#ze(i),n.show(i))}get view(){return this.#w}get edited(){return this.#O}get generation(){return this.#k}geometry(){let e=this.#r.current;return this.#e!==`ready`||!e||this.#Ae()?null:this.#Me(()=>e.geometryJson())}sgroups(){let e=this.#r.current;if(this.#e!==`ready`||!e||this.#Ae()||!e.sgroupsJson)return null;try{return e.sgroupsJson()}catch(e){return this.#je(e),null}}sgroupRefusal(e){let t=this.#r.current;if(this.#e!==`ready`||!t||!e.ops.length||this.#Ae())return null;let n=null;try{return n=t.clone(),C(n,e,this.#j),null}catch(e){return e instanceof oe||e instanceof Re?e.message:(this.#je(e),null)}finally{G(n)}}contractible(e,t=[]){let n=this.#r.current;if(this.#e!==`ready`||!n||!e.length||this.#Ae()||!n.contractToAbbreviation)return!1;let r=null;try{return r=(t.length?this.#xe(n,t):n).clone(),r.contractToAbbreviation([...e])!==null}catch(e){return this.#je(e),!1}finally{G(r)}}groupPreview(e,t,n){let r=this.#r.current,i=this.#s;if(this.#e!==`ready`||!r||!e.length||!i?.exportSvg||this.#Ae())return null;let a=[...new Set(e)].sort((e,t)=>e-t),o=null;try{o=r.extract(a);let e=JSON.parse(o.sgroupsJson?.()??`[]`);for(let t of e)t.type===`SUP`&&t.contracted&&o.expandAbbreviation?.(t.index);a.forEach((e,n)=>{let r=t.get(e);r&&o.setAtomAttachmentPoints(n,r)});let s=i.exportSvg(o,n,{edited:!0});return s?sn(s,Mr()):null}catch(e){return this.#je(e),null}finally{G(o)}}stereo(e=!1){let t=this.#r.current;return this.#e!==`ready`||!t||this.#Ae()?ae:this.#_e(t,e)}aromaticity(e,t,n=[]){let r=this.#r.current;if(this.#e!==`ready`||!r||this.#p||this.#Ae()||!(e===`aromatize`?r.aromatize:r.dearomatize))return null;let i=Ot(t);if(i&&!i.length)return Se;let a=`${e}:${i===null?`*`:i.join(`,`)}`,o=this.#ge.get(r)??new Map;this.#ge.set(r,o);let s=o.get(a);if(s)return s;let c=this.#Ce(r);if(!c)return null;let l=null;try{n.length&&C(c,{label:`prepare`,ops:n});let t=e===`aromatize`?c.aromatize:c.dearomatize;l=It(t.call(c,i===null?null:Uint32Array.from(i)))}catch(e){l=this.#je(e)?null:Se}finally{try{c.free()}catch{}}return l&&o.set(a,l),l}#ge=new WeakMap;#_e(e,t){let n=this.#X.get(e);if(n&&(n.cip||!t))return n.info;if(!e.stereoInfo)return ae;let r;try{r=e.isEmpty()?ae:Te(e.stereoInfo(t))}catch(e){if(this.#je(e))throw e;r=ae}return this.#X.set(e,{cip:t,info:r}),r}stereoisomers(){let e=this.#r.current;if(this.#e!==`ready`||!e||this.#Ae()||!e.stereoisomerCount)return null;if(this.#ve.has(e))return this.#ve.get(e)??null;let t;try{t=e.isEmpty()?null:e.stereoisomerCount(),t!==null&&!(t>0)&&(t=null)}catch(e){if(this.#je(e))throw e;t=null}return this.#ve.set(e,t),t}#ve=new WeakMap;get axisGroups(){let e=this.#r.current;return!!(e?.setStereoGroup&&e.setStereoGroup.length>=4)}get stereoShown(){return this.#Z}set stereoShown(e){let t=this.#Z;(e.labels!==t.labels||e.flag!==t.flag||e.plain!==t.plain||e.registration!==t.registration)&&(this.#Z=Object.freeze({labels:e.labels,flag:e.flag,plain:e.plain,registration:e.registration}),this.refresh())}#ye=p;set strings(e){e!==this.#ye&&(this.#ye=e,this.#Z.flag&&this.refresh())}#be(e,t,n=!1){let r=n&&this.#Z.labels?{...this.#Z,labels:!1}:this.#Z,i=n?this.#r.current:e,a=i&&(r.labels||r.flag)?this.#_e(i,r.labels):null,o=fn(a,r,{absolute:this.#ye.flagAbsolute,racemic:this.#ye.flagRacemic,relative:this.#ye.flagRelative});return Object.keys(o).length?{edited:t,stereo:o}:{edited:t}}doubleBondPosition(e){let t=this.#r.current;return this.#e!==`ready`||!t||this.#Ae()?`auto`:this.#Me(()=>t.doubleBondPosition(e))??`auto`}edit(e){let t=this.#r.current,n=this.#s,r=this.#w;if(this.#G=null,this.#e!==`ready`||!t||!n||!r||!e.ops.length||this.#Ae())return!1;if(ti(e))return this.#Ne(Eg()),!1;this.#W=[],this.#P=null;let i;try{i=t.clone()}catch(e){return this.#Ne(this.#je(e)??Tg(e)),!1}let a=l(e),o=e.prepared??0,s=null,c=[];try{let t;if(a&&o?(C(i,{label:e.label,ops:e.ops.slice(0,o),prepared:o},this.#j),s=this.#Ce(i),t=C(i,{label:e.label,ops:e.ops.slice(o)},this.#j,c)):t=C(i,e,this.#j,c),!t||a&&this.#we(s,i,a))return this.#G=c[0]??null,G(i),!1}catch(e){let t=this.#je(e),n=e instanceof ue||e instanceof ve||e instanceof Re||e instanceof oe||e instanceof qe,r=e instanceof qe?Zg(e):null;if(this.#Ne(t??r??Dg(`unknown`,e instanceof ue||e instanceof ve&&e.reason!==null?e.reason:e)),t||n)return G(i),!1}finally{G(s)}let u;try{u={drawing:n.draw(i,e.fit?`fill`:r,this.#be(i,!0)),empty:i.isEmpty(),weight:qg(i)}}catch(e){let t=this.#je(e);return G(i),this.#Ne(t??Dg(`unknown`,e)),!1}return this.#r.push(i,u.weight),this.#Y.set(i,e.label),this.#O=!0,this.#Se(),this.#p=u.empty,e.ops.some(e=>e.op===`paste`)&&(this.#K=Object.freeze([...this.#W]),this.#F=this.#P),this.#W=[],this.#P=null,e.fit&&this.#ze(u.drawing),n.show(u.drawing),this.#J={kind:`edit`,label:e.label,said:e.said??null},this.#c.emit(`change`,{source:e.source??`user`}),!0}refuseRead(e){this.#e!==`destroyed`&&this.#Ne(e)}preview(e,t=!1){let n=this.#r.current,r=this.#s,i=this.#w;if(this.#e!==`ready`||!n||!r||!i||!e.ops.length||ti(e)||this.#Ae())return null;let a=null;try{if(t)return a=this.#t.empty(),C(a,e,this.#j)?r.draw(a,i,this.#be(a,!0,!0)):null;let o=e.prepared??0;return a=(o?this.#xe(n,e.ops.slice(0,o)):n).clone(),C(a,o?{label:e.label,ops:e.ops.slice(o)}:e,this.#j)?r.draw(a,i,this.#be(a,!0,!0)):null}catch(e){return this.#je(e),null}finally{G(a)}}previewText(e){let t=this.#r.current,n=this.#s,r=this.#w;if(!n?.drawText)return;if(this.#e!==`ready`||!t||!r||!e.ops.length||ti(e)||this.#Ae())return null;let i=null;try{let a=e.prepared??0;return i=(a?this.#xe(t,e.ops.slice(0,a)):t).clone(),C(i,a?{label:e.label,ops:e.ops.slice(a)}:e,this.#j)?n.drawText(i,r,this.#be(i,!0,!0)):null}catch(e){return this.#je(e),null}finally{G(i)}}partDrawing(e,t=[]){let n=this.#r.current,r=this.#s,i=this.#w;if(this.#e!==`ready`||!n||!r||!i||!e.length||this.#Ae())return null;let a=null;try{return a=(t.length?this.#xe(n,t):n).extract([...e]),r.draw(a,i,this.#be(a,!0,!0))}catch(e){return this.#je(e),null}finally{G(a)}}#xe(e,t){let n=this.#A,r=JSON.stringify(t);if(n.copy&&n.of===e&&n.key===r)return n.copy;n.free();let i=e.clone();try{C(i,{label:`prepare`,ops:t},this.#j)}catch(e){throw G(i),e}return Object.assign(n,{of:e,key:r,copy:i}),i}#Se(){this.#m.clear(),this.#B=null,this.#V=null,this.#h=null,this.#g=null,this.#_=null,this.#v={},this.#y.clear(),this.#x=null,this.#A.free()}#Ce(e){try{return e.clone()}catch(e){if(this.#je(e))throw e;return null}}#we(e,t,n){let r=e??this.#r.current;if(!r)return!1;try{return n.every(e=>_t(r.atomJson(e),t.atomJson(e)))?Ug.every(n=>(e?this.#Pe(n,e):this.#Fe(n))===this.#Pe(n,t)):!1}catch(e){if(this.#je(e))throw e;return!1}}get stateToken(){return this.#e===`ready`?this.#r.current:null}takeBack(e){let t=this.#s,n=this.#w,r=this.#r.index;if(this.#e!==`ready`||this.#r.current!==e||r===0||this.#r.canRedo||!t||!n||this.#Ae())return!1;let i=this.#r.at(r-1),a=i!==this.#i,o=this.#Me(()=>({drawing:t.draw(i,n,this.#be(i,a)),empty:i.isEmpty()}));return o?(this.#J={kind:`undo`,label:this.#Y.get(this.#r.current)??null},this.#r.pop(),this.#O=a,this.#Se(),this.#p=o.empty,t.show(o.drawing),this.#c.emit(`change`,{source:`user`}),!0):!1}abbreviationList(){if(this.#Te)return this.#Te;if(this.#e!==`ready`||!this.#t?.abbreviationsJson||this.#Ae())return[];try{let e=JSON.parse(this.#t.abbreviationsJson());this.#Te=Array.isArray(e)?e:[]}catch(e){if(this.#je(e))throw e;this.#Te=[]}return this.#Te}#Te=null;inspect(e,t,n,r=!1){if(this.#e!==`ready`||!this.#t||this.#Ae())return null;let i=null;this.#U=null;try{return i=this.#fe(e,t,r).mol,n(i)}catch(e){return this.#U=this.#je(e)??bg(e),null}finally{G(i)}}readOnce(e,t,n,r){if(this.#e!==`ready`||!this.#t||this.#Ae())return null;let i=null,a;this.#U=null;try{i=this.#fe(e,t),a=n(i.mol)}catch(e){return this.#U=this.#je(e)??bg(e),G(i?.mol??null),null}this.#Ee(),this.#M={text:e,format:t,...i};try{return r(a)}finally{this.#Ee()}}#Ee(){let e=this.#M;this.#M=null,G(e?.mol??null)}copyText(e,t=[],n,r){let i=this.#r.current;if(this.#e!==`ready`||!i||!e.length||this.#Ae())return null;let a=null;try{a=(t.length?this.#xe(i,t):i).extract([...e]);let o=(n===`smarts`?this.#Pe(`smarts`,a):n?xr(a,n):mr(a))||null;return o&&n&&r&&r(Pt(a,n)),o}catch(e){return this.#je(e),null}finally{G(a)}}copySvg(e){let t=this.#r.current,n=this.#s;if(this.#e!==`ready`||!t||!e.length||!n?.exportSvg||this.#Ae())return``;let r=null;try{r=t.extract([...e]);let i=n.exportSvg(r,this.#w?.docLen??t.meanBondLength(),{...this.#be(r,!0),text:!0});return i?sn(i,Mr()):``}catch(e){return this.#je(e),``}finally{G(r)}}copyCdxml(e,t=[],n){let r=this.#r.current;if(this.#e!==`ready`||!r||!e.length||!this.#I.available||this.#Ae())return null;let i=null,a;try{i=(t.length?this.#xe(r,t):r).extract([...e]),a=Ba(i)}catch(e){return this.#je(e),null}finally{G(i)}let o=e=>{try{let t=Ha(e,a);return t.text&&n&&n({asked:`cdxml`,as:`cdxml`,lost:[...new Set(t.warnings.flatMap(e=>e.code?[e.code]:[]))],why:[]}),t.text||null}catch{return null}},s=this.#re();return s?o(s):this.#I.load().then(o,()=>null)}copyInchi(e,t,n,r){let i=this.#r.current;if(this.#e!==`ready`||!i||!e.length||!this.#L.available||this.#Ae())return null;let a=null,o;try{a=(t.length?this.#xe(i,t):i).extract([...e]),o=Bg.molfile(a)}catch(e){return this.#je(e),null}finally{G(a)}let s=e=>{try{let t=oo(e,o),i=n===`inchikey`?t.key:t.inchi;return r?.({asked:n,as:n,lost:so(t),why:[],...i?{}:{refused:!0}}),i||null}catch{return null}},c=this.#ie();return c?s(c):this.#L.load().then(s,()=>null)}formatOf(e){return this.#e===`ready`?this.#Ve(e):`unknown`}atomRNumbers(e){return this.#De(e,t=>[...t.atomRNumbers(e)])}atomAttachmentPoints(e){return this.#De(e,t=>t.atomAttachmentPoints(e))}#De(e,t){let n=this.#r.current;if(this.#e!==`ready`||!n||this.#Ae()||!Number.isInteger(e)||e<0)return null;try{return e<n.numAtoms()?t(n):null}catch(e){return this.#je(e),null}}atomJson(e){let t=this.#r.current;if(this.#e!==`ready`||!t||this.#Ae()||!Number.isInteger(e)||e<0)return null;try{return e<t.numAtoms()?t.atomJson(e):null}catch(e){return this.#je(e),null}}bondQueryJson(e){let t=this.#r.current;if(this.#e!==`ready`||!t||this.#Ae()||!Number.isInteger(e)||e<0||typeof t.bondQueryJson!=`function`)return null;try{return e<t.numBonds()?t.bondQueryJson(e):null}catch(e){return this.#je(e),null}}aromaticRingBonds(){let e=this.#r.current;if(this.#e!==`ready`||!e||this.#Ae())return[];try{return Qg(this.#t,e)}catch(e){return this.#je(e),[]}}tryQuery(e,t=null,n=null){let r=this.#r.current;if(this.#e!==`ready`||!r||this.#Ae())return{ok:!0,query:null};let i=null;try{return i=r.clone(),e&&C(i,e,this.#j),t===null?{ok:!0,query:n!==null&&n>=0&&n<i.numBonds()?i.bondQueryJson(n):null}:{ok:!0,query:t>=0&&t<i.numAtoms()?i.atomQueryJson(t):null}}catch(e){return e instanceof qe?{ok:!1,reason:e.reason,offset:e.offset}:(this.#je(e),{ok:!0,query:null})}finally{G(i)}}showView(e){let t=this.#r.current,n=this.#s;if(this.#e!==`ready`||!t||!n||this.#Ae())return!1;let r=this.#Me(()=>n.draw(t,e,this.#be(t,this.#O)));return r?(this.#ze(r),n.show(r),!0):!1}moveView(e){return this.#e!==`ready`||!this.#r.current||!this.#s||!this.#w||this.#Ae()?!1:(this.#w=e,!0)}fit(e=!1){let t=this.#r.current,n=this.#s;if(this.#e!==`ready`||!t||!n||this.#Ae())return!1;let r=this.#Me(()=>n.draw(t,e?`fit`:`fill`,this.#be(t,this.#O)));return r?(this.#ze(r),n.show(r),!0):!1}refresh(){let e=this.#r.current,t=this.#s,n=this.#w;if(this.#e!==`ready`||!e||!t||!n||this.#Ae())return;let r=this.#Me(()=>t.draw(e,n,this.#be(e,this.#O)));r&&t.show(r)}get canUndo(){return this.#e===`ready`&&this.#r.canUndo}get canRedo(){return this.#e===`ready`&&this.#r.canRedo}undo(e=`api`){return this.#Oe(e),this.#ke(-1,e)}redo(e=`api`){return this.#Oe(e),this.#ke(1,e)}#Oe(e){e===`api`&&this.#e===`ready`&&this.#D?.()}#ke(e,t){if(this.#e===`destroyed`)throw Error(gg);let n=this.#s,r=this.#w,i=this.#r.index+e,a=(e<0?this.canUndo:this.canRedo)?this.#r.at(i):null;if(!a||!n||!r||this.#Ae())return!1;let o=a!==this.#i,s=this.#Me(()=>({drawing:n.draw(a,r,this.#be(a,o)),empty:a.isEmpty()}));return s?(this.#E?.(),this.#J={kind:e<0?`undo`:`redo`,label:this.#Y.get(e<0?this.#r.current:a)??null},this.#r.moveTo(i),this.#O=o,this.#Se(),this.#p=s.empty,n.show(s.drawing),this.#c.emit(`change`,{source:t}),!0):!1}on(e,t){return this.#e===`destroyed`?()=>{}:this.#c.on(e,t)}off(e,t){this.#c.off(e,t)}settingsApplied(e){this.#e!==`destroyed`&&this.#c.emit(`settings`,e)}rendered(){this.#e!==`destroyed`&&this.#c.emit(`render`,Object.freeze({renderCount:this.renderCount}))}selectionChanged(e){this.#e!==`destroyed`&&this.#c.emit(`selection`,e)}noticed(e){this.#e!==`destroyed`&&this.#c.emit(`notice`,e)}#Ae(){return this.#n?this.#n:this.#t?Xg.get(this.#t)??null:null}#je(e){return!xg(e)||!this.#t?null:(Xg.has(this.#t)||Xg.set(this.#t,Sg(e)),Xg.get(this.#t))}#Me(e){let t=this.#Ae();if(t)return this.#Ne(t),null;try{return e()}catch(e){return this.#Ne(this.#je(e)??Dg(`unknown`,e)),null}}#Ne(e){this.#f=e,this.#c.emit(`error`,e)}#Pe(e,t){return e===`molfile`||e===`molV3000`?jg(t,Bg[e]):e===`smarts`?Ng(this.#t,t):Bg[e](t)}#Fe(e){let t=this.#m.get(e);if(typeof t==`string`)return t;if(t)throw t;let n=this.#r.current;if(!n)return``;let r=``;if(!this.#p){let t=this.#Ae();if(t)r=t;else if(e===`smiles`)r=this.#Ie();else if(e===`canonicalSmiles`)r=this.#Ie(!0);else try{r=this.#Pe(e,n)||(Hg.has(e)&&this.hasQuery?``:wg(Rg[e],Error(`it wrote nothing`)))}catch(t){r=this.#je(t)??wg(Rg[e],t)}}if(this.#m.set(e,r),typeof r==`string`)return r;throw r}#Ie(e=!1){let t;try{t=this.#Fe(e?`canonicalPlain`:`plainSmiles`)}catch(e){return e}let n;try{n=this.#Fe(e?`canonicalCxsmiles`:`cxsmiles`)}catch{return this.#Ae()||(e||(this.#h=`smiles`),t)}try{let r=Pg(t,n,e=>Fg(this.#t,e));return e||(this.#h=r.format),r.text}catch(e){return this.#je(e)??wg(Rg.smiles,e)}}#Le(e,t,n,r=!1){let i=this.#Ae();if(i)return this.#Ne(i),!1;let a,o,s;try{({mol:a,drops:o,info:s}=this.#fe(e,t))}catch(e){return this.#Ne(this.#je(e)??bg(e)),!1}let c,l,u=null;try{c=a.isEmpty(),l=qg(a);let e=r&&this.#w&&!this.#T?this.#w:`fit`;n&&(u=this.#s?.draw(a,e,this.#be(a,!1))??null)}catch(n){let r=this.#je(n);return G(a),this.#Ne(r??Dg(t??this.#Ve(e),n)),!1}return this.#Re(a,c,l,e),n||(this.#a=a),this.#K=Object.freeze(o),this.#F=s,this.#f=null,n&&(u&&(this.#ze(u),this.#s?.show(u)),this.#J={kind:`write`,label:null},this.#c.emit(`change`,{source:`api`})),!0}#Re(e,t,n,r){this.#r.restart(e,n),this.#i=e,this.#Se(),this.#p=t,this.#u=r,this.#O=!1,this.#k+=1}#ze(e){e.view&&(this.#w=e.view,this.#T=!!e.blind)}#Be(e){try{let t=e?this.#fe(e.text,e.format):null,n=t?t.mol:this.#t.empty();this.#Re(n,n.isEmpty(),qg(n),e?.text??``),this.#a=n,this.#K=Object.freeze(t?.drops??[]),this.#F=t?.info??null}catch(t){!this.#je(t)&&e&&this.#Be(void 0)}}#Ve(e){if(this.#Ae())return`unknown`;try{return this.#t.detectFormat(e)}catch(e){return this.#je(e),`unknown`}}},K=`dialog.bond-properties`,n_=class{#e;#t;#n=null;constructor(e,t){this.#e=e,this.#t=t}get open(){return this.#n!==null}show(e,t){if(this.#n||!this.#t.isConnected)return;let n=this.#e.bondChoice(e);if(!n)return;let r=this.#e.bondForm?.(e)??null,i=n.shared===`crossed`,a=i?`double`:n.shared,o=r&&i&&r.form.type===`crossed`?{...r.form,type:`double`}:r?.form;this.#n={bondIds:[...e],opened:a,crossed:i,opener:t,query:r?.shows&&o?{opened:o,form:o,error:null}:null},this.#f();let s=this.#a(),c=this.#o();if(!s||!c)return;let l=this.#n.query?this.#n.query.form.type:a;l?c.value=l:c.selectedIndex=-1,this.#c(),s.showModal(),c.focus()}cancel(){let e=this.#n;e&&(this.#i(),this.#e.restoreFocus(e.opener))}#r(){let e=this.#n;if(!e||e.query?.error)return;let t=this.#o()?.value??``;if(e.query&&this.#e.applyForm){let n={...e.query.form,type:t||null},r=this.#e.applyForm(e.bondIds,e.query.opened,n);if(!r.ok){e.query.error={reason:r.reason,offset:r.offset},this.#f(),this.#s(`smarts`)?.focus();return}this.#i(),this.#e.restoreFocus(e.opener);return}let n=vr.includes(t)?t:null;this.#i(),n&&!(e.crossed&&n===e.opened)&&this.#e.apply(e.bondIds,n),this.#e.restoreFocus(e.opener)}#i(){this.#a()?.close(),this.#n=null,this.#f()}#a(){return this.#t.querySelector(`dialog[data-testid="${K}"]`)}#o(){return this.#t.querySelector(`[data-testid="${K}.type"]`)}#s(e){return this.#t.querySelector(`[data-testid="${K}.${e}"]`)}#c(){let e=this.#n?.query;if(!e)return;let t=this.#s(`topology`);t&&(t.value=e.form.topology??``);let n=this.#s(`custom`);n&&(n.checked=e.form.custom===!0,n.indeterminate=e.form.custom===null);let r=this.#s(`smarts`);r&&(r.value=e.form.smarts)}#l(e){let t=this.#n,n=t?.query;if(t&&n){if(e){let e={...n.form,type:this.#o()?.value||null,custom:!1},r=this.#e.bondSmarts?.(t.bondIds,n.opened,e)??``;n.form={...n.form,custom:!0,smarts:r,filled:r}}else n.form={...n.form,custom:!1};n.error=null,this.#f(),this.#c(),e&&this.#s(`smarts`)?.focus()}}#u(e){if(e.key===`Escape`){e.preventDefault(),this.cancel();return}if(e.key===`Tab`){kr(e);return}e.key===`Enter`&&(e.composedPath()[0]===this.#o()||e.composedPath()[0]===this.#s(`smarts`))&&(e.preventDefault(),this.#r())}#d(e){let t=this.#e.strings(),n=e.form.custom===!0,r=e.error?t.format(`smartsNotRead`,{position:(e.error.offset??0)+1,message:Qe(e.error.reason)}):``;return x`<label for="bp-topology">${t.topology}</label>
      <select
        id="bp-topology"
        data-testid="${K}.topology"
        ?disabled=${n}
        @change=${e=>{this.#n?.query&&(this.#n.query.form={...this.#n.query.form,topology:e.target.value||null})}}
      >
        ${e.form.topology===null?x`<option value="" selected>${t.severalValues}</option>`:w}
        ${de.map(e=>x`<option value=${e}>${jn(e,t)}</option>`)}
      </select>
      <span class="q-check bp-custom"
        ><input id="bp-custom" type="checkbox" data-testid="${K}.custom" @change=${e=>this.#l(e.target.checked)} /><label for="bp-custom"
          >${t.customQuery}</label
        ></span
      >
      <input
        id="bp-smarts"
        type="text"
        autocomplete="off"
        spellcheck="false"
        data-testid="${K}.smarts"
        aria-label=${t.customQuerySmarts}
        aria-invalid=${e.error?`true`:w}
        aria-describedby=${e.error?`bp-error`:w}
        ?disabled=${!n}
        @input=${e=>{let t=this.#n;t?.query&&(t.query.form={...t.query.form,smarts:e.target.value},t.query.error&&(t.query.error=null,this.#f()))}}
      />
      <span class="error" id="bp-error" role="alert" data-testid="${K}.error" ?hidden=${!e.error}>${r}</span>`}#f(){let e=this.#n,t=this.#e.strings(),n=e?.query??null,r=[...vr,...n?Tt:[]],i=n?.form.custom===!0;T(e?x`<dialog
            class="dialog"
            data-testid=${K}
            aria-labelledby="bp-title"
            @keydown=${e=>this.#u(e)}
            @cancel=${e=>{e.preventDefault(),this.cancel()}}
          >
            <div class="dialog-head">
              <h2 id="bp-title">${t.bondProperties}</h2>
              <button type="button" class="close" data-testid="${K}.close" aria-label=${t.close} @click=${()=>this.cancel()}>×</button>
            </div>
            <div class="fields">
              <label for="bp-type">${t.fieldType}</label>
              <select id="bp-type" class="list" size=${r.length} data-testid="${K}.type" ?disabled=${i} aria-describedby=${e.crossed?`bp-ez`:w}>
                ${r.map(e=>pr(e)?x`<option value=${e} ?disabled=${this.#e.queryOff?.(e)??!1}>
                        ${this.#e.queryOff?.(e)?t.format(`queryElementOff`,{name:ut(e,t)}):ut(e,t)}
                      </option>`:x`<option value=${e}>${D(`bond.${e}`,t)}</option>`)}
              </select>
              ${e.crossed?x`<p class="bp-ez" id="bp-ez" data-testid="${K}.ez">${t.bondEzUnknown}</p>`:w}
              ${n?this.#d(n):w}
            </div>
            <div class="dialog-foot">
              <button type="button" class="primary" data-testid="${K}.apply" ?disabled=${!!n?.error} @click=${()=>this.#r()}>${t.apply}</button>
              <button type="button" data-testid="${K}.cancel" @click=${()=>this.cancel()}>${t.cancel}</button>
            </div>
          </dialog>`:w,this.#t)}},r_=`dialog.dissolve`,i_=class{#e;#t;#n=null;constructor(e,t){this.#e=e,this.#t=t}get open(){return this.#n!==null}show(e,t){if(!this.#t.isConnected){this.#e.answer(!1);return}this.#n||(this.#n={labels:[...e],opener:t},this.#s(),this.#i()?.showModal(),this.#a(`remove`)?.focus())}cancel(){this.#r(!1)}#r(e){let t=this.#n;t&&(this.#i()?.close(),this.#n=null,this.#s(),this.#e.answer(e),this.#e.restoreFocus(t.opener))}#i(){return this.#t.querySelector(`dialog[data-testid="${r_}"]`)}#a(e){return this.#t.querySelector(`[data-testid="${r_}.${e}"]`)}#o(e){if(e.key===`Escape`){e.preventDefault(),this.cancel();return}e.key===`Tab`&&kr(e)}#s(){let e=this.#n,t=this.#e.strings();T(e?x`<dialog
            class="dialog"
            role="alertdialog"
            data-testid=${r_}
            aria-labelledby="dd-title"
            aria-describedby="dd-question dd-hint"
            @keydown=${e=>this.#o(e)}
            @cancel=${e=>{e.preventDefault(),this.cancel()}}
          >
            <div class="dialog-head">
              <h2 id="dd-title">${t.dissolveTitle}</h2>
              <button type="button" class="close" data-testid="${r_}.close" aria-label=${t.close} @click=${()=>this.cancel()}>×</button>
            </div>
            <p id="dd-question" class="dissolve-question">${t.dissolveQuestion(e.labels)}</p>
            <p id="dd-hint" class="rg-hint">${t.dissolveHint}</p>
            <div class="dialog-foot">
              <button type="button" class="primary" data-testid="${r_}.remove" @click=${()=>this.#r(!0)}>${t.dissolveRemove}</button>
              <button type="button" data-testid="${r_}.cancel" @click=${()=>this.cancel()}>${t.cancel}</button>
            </div>
          </dialog>`:w,this.#t)}},a_=class extends sr{constructor(e){if(super(e),this.it=w,e.type!==Kn.CHILD)throw Error(this.constructor.directiveName+`() can only be used in child bindings`)}render(e){if(e===w||e==null)return this._t=void 0,this.it=e;if(e===Ke)return e;if(typeof e!=`string`)throw Error(this.constructor.directiveName+`() called with a non-string value`);if(e===this.it)return this._t;this.it=e;let t=[e];return t.raw=t,this._t={_$litType$:this.constructor.resultType,strings:t,values:[]}}};a_.directiveName=`unsafeHTML`,a_.resultType=1,pn(a_);var o_=class extends a_{};o_.directiveName=`unsafeSVG`,o_.resultType=2;var s_=pn(o_),c_=new Set([`Shift`,`Control`,`Alt`,`AltGraph`,`Meta`,`OS`,`Super`,`Hyper`,`Fn`,`FnLock`,`CapsLock`,`NumLock`,`ScrollLock`,`Symbol`,`SymbolLock`]),l_=new Set([`Backspace`,`Delete`,`Enter`,`Escape`,`Tab`,`ArrowLeft`,`ArrowRight`,`ArrowUp`,`ArrowDown`,`Home`,`End`,`PageUp`,`PageDown`,`Insert`,`Dead`,`Compose`,`Process`,`Unidentified`]),u_=new Set([`a`,`c`,`v`,`x`,`z`,`y`,`Backspace`,`Delete`,`ArrowLeft`,`ArrowRight`,`ArrowUp`,`ArrowDown`,`Home`,`End`]),d_=new Set([`Enter`,` `]);function f_(e,t){if(c_.has(e.key))return!1;let n=!!(e.ctrlKey||e.metaKey),r=!!e.altKey;switch(t){case`text`:return e.isComposing?!0:e.key.length===1?!n||e.ctrlKey&&r&&!e.metaKey?!0:!r&&u_.has(e.key.toLowerCase()):n?!r&&u_.has(e.key):l_.has(e.key);case`dialog`:return!0;case`button`:return!n&&!r&&d_.has(e.key);default:return!1}}var p_=new Set([`mouseup`,`contextmenu`,`lostpointercapture`,`pointerout`,`pointerleave`,`mouseout`,`mouseleave`]),m_=new Set([`pointermove`,`pointerover`,`pointerenter`,`pointerout`,`pointerleave`,`mousemove`,`mouseover`,`mouseenter`,`mouseout`,`mouseleave`]),h_=new Set([`click`,`auxclick`,`dblclick`]),g_=e=>e.type.startsWith(`touch`)||e.pointerType===`touch`||h_.has(e.type),__=e=>e.pointerId===void 0||e.pointerType===`mouse`?`mouse`:`${e.pointerType??``} ${e.pointerId}`,v_=class{#e=new Map;owner(e){if(g_(e))return`sketcher`;let t=__(e),n=this.#e.get(t);if(n&&(n.ended?!p_.has(e.type):m_.has(e.type)&&e.buttons===0)&&(this.#e.delete(t),n=void 0),e.type===`pointerdown`||e.type===`mousedown`&&!n)n=this.#t(t,`sketcher`);else if(!n){if(e.buttons===0)return`sketcher`;n=this.#t(t,`host`)}return e.type===`gotpointercapture`?n.captured=!0:e.type===`lostpointercapture`?n.captured=!1:(e.type===`pointercancel`||(e.type===`pointerup`||e.type===`mouseup`)&&e.buttons===0||e.leaves&&!n.captured)&&(n.ended=!0),n.owner}peek(e){if(g_(e))return`sketcher`;let t=this.#e.get(__(e));return!t||(t.ended?!p_.has(e.type):e.buttons===0)?`sketcher`:t.owner}get size(){return this.#e.size}clear(){this.#e.clear()}#t(e,t){for(let[t,n]of this.#e)n.ended&&t!==e&&this.#e.delete(t);let n={owner:t,captured:!1,ended:!1};return this.#e.set(e,n),n}};function y_(e){return e.type===`pointerleave`||e.type===`mouseleave`||e.type===`pointercancel`||e.buttons===0}var b_=[`pointerdown`,`pointermove`,`pointerup`,`pointercancel`,`pointerover`,`pointerout`,`gotpointercapture`,`lostpointercapture`,`mousedown`,`mousemove`,`mouseup`,`mouseover`,`mouseout`,`click`,`dblclick`,`auxclick`,`touchstart`,`touchmove`,`touchend`,`touchcancel`],x_=[`pointerdown`,`pointermove`,`pointerup`,`pointercancel`,`pointerover`,`pointerout`,`gotpointercapture`,`lostpointercapture`,`mousedown`,`mousemove`,`mouseup`,`mouseover`,`mouseout`,`contextmenu`],S_=new Set([`pointerenter`,`pointerleave`,`mouseenter`,`mouseleave`]),C_=new WeakSet,w_=new WeakMap;function q(e){if(C_.has(e))return!0;if(!S_.has(e.type))return!1;let t=e.currentTarget??e.target,n=typeof t?.getRootNode==`function`?t.getRootNode():null;return!!(n&&w_.get(n)?.hosts(e))}var T_=[`beforeinput`,`input`,`compositionstart`,`compositionupdate`,`compositionend`],E_=[`wheel`,`copy`,`cut`,`paste`,`dragenter`,`dragover`,`dragleave`,`drop`],D_=new Set([`copy`,`cut`,`paste`]),O_=new Set([`text`,`search`,`number`,`email`,`url`,`tel`,`password`]),k_=new Set([`wheel`,`touchstart`,`touchmove`]);function A_(e){let t=e.composedPath()[0];if(!t||typeof t.closest!=`function`)return`canvas`;let n=t.localName;return n===`textarea`||n===`select`||t.isContentEditable||n===`input`&&O_.has(t.type)?`text`:t.closest(`dialog`)?`dialog`:t.closest(`button, [role="menuitem"]`)?`button`:`canvas`}function j_(e){let t=e.composedPath()[0];return t&&typeof t.closest==`function`&&t.closest(`dialog`)?`dialog`:A_(e)}var M_=e=>e.code||e.key,N_=class{#e;#t=new WeakSet;#n=new Set;#r=new v_;#i=[];#a;constructor(e,t){this.#e=e,this.#a=t??null,w_.set(e,this);let n=e=>e.stopPropagation();for(let e of x_)this.#o(e,e=>this.#s(e),!0);for(let e of b_)this.#o(e,e=>{C_.has(e)||e.stopPropagation()});for(let e of T_)this.#o(e,n);for(let e of E_)this.#o(e,e=>{(e.defaultPrevented||this.#t.has(e)||D_.has(e.type)&&A_(e)===`text`)&&e.stopPropagation()});this.#o(`keydown`,e=>this.#u(e)),this.#o(`keypress`,e=>{this.#n.has(M_(e))&&e.stopPropagation()}),this.#o(`keyup`,e=>{this.#n.delete(M_(e))&&e.stopPropagation()})}own(e){this.#t.add(e)}hosts(e){return this.#r.peek(this.#c(e))===`host`}destroy(){let e=this.#e;this.#e=null,e&&w_.delete(e);for(let[t,n,r]of this.#i)e?.removeEventListener(t,n,r);this.#i=[],this.#n.clear(),this.#r.clear(),this.#a=null}#o(e,t,n=!1){let r={passive:k_.has(e),capture:n};this.#e?.addEventListener(e,t,r),this.#i.push([e,t,r])}#s(e){this.#r.owner(this.#c(e))===`host`&&(C_.add(e),this.#a?.())}#c(e){let t=`pointerId`in e?e:null,n=e.type===`pointerout`||e.type===`mouseout`;return{type:e.type,pointerId:t?.pointerId,pointerType:t?.pointerType,buttons:e.buttons,leaves:n&&!this.#l(e.relatedTarget)}}#l(e){return e instanceof Node&&!!this.#e?.contains(e)}#u(e){let t=M_(e);e.defaultPrevented||this.#t.has(e)||f_(e,j_(e))?(e.stopPropagation(),this.#n.add(t)):this.#n.delete(t)}},P_=(e,t)=>e.x>=t.left&&e.x<=t.right&&e.y>=t.top&&e.y<=t.bottom;function F_(e,t){let{button:n,tip:r}=t;return P_(e,{left:Math.min(n.left,r.left),top:Math.min(n.top,r.top),right:Math.max(n.right,r.right),bottom:Math.max(n.bottom,r.bottom)})}var I_=class{#e;#t=null;#n=null;#r=null;#i=-1/0;constructor(e){this.#e=e}get open(){return this.#t}enter(e){if(this.#n=e,this.#s(),this.#t!==e){if(this.#t!==null||this.#e.now()-this.#i<500){this.#a(e);return}this.#r=this.#e.after(300,()=>{this.#r=null,this.#n===e&&this.#a(e)})}}leave(e,t){if(this.#n===e&&(this.#n=null),this.#s(),this.#t!==e)return;let n=this.#e.boxes(e);t&&n&&F_(t,n)||this.#o()}moved(e){if(this.#t===null)return;let t=this.#e.boxes(this.#t);(!t||!F_(e,t))&&this.#o()}pressed(){this.#n=null,this.#s(),this.#o()}gone(){this.pressed()}destroy(){this.#n=null,this.#s(),this.#t=null}#a(e){this.#t!==null&&this.#t!==e&&(this.#i=this.#e.now()),this.#t=e,this.#e.changed()}#o(){this.#t!==null&&(this.#t=null,this.#i=this.#e.now(),this.#e.changed())}#s(){this.#r?.(),this.#r=null}},J=`panel.stereo`,L_=`menu.stereo-group`,R_=6,z_=4,B_=[`remove-stereo`,`recompute-wedges`,`close`],V_=e=>`sp-tip-${e}`,H_=e=>B_.includes(e),U_=`sp-tip-row`,W_=e=>`sp-mean-${e}`,G_=e=>e===`centre`?_e`<path d="M8 8 L14 5.5 L14 10.5 Z" fill="currentColor"/><path d="M8 8 L2.5 4 M8 8 L2.5 12 M8 8 L8 2" stroke="currentColor" stroke-width="1.4" fill="none" stroke-linecap="round"/>`:e===`double`?_e`<path d="M2 6 L14 6 M2 10 L14 10" stroke="currentColor" stroke-width="1.5" stroke-linecap="round"/>`:_e`<path d="M8 1.5 L8 14.5" stroke="currentColor" stroke-width="1.5" stroke-linecap="round"/><path d="M3.2 6.2 A5 2.6 0 1 0 6 5.1" stroke="currentColor" stroke-width="1.3" fill="none"/><path d="M6.6 3.6 L6 5.2 L7.7 5.6" stroke="currentColor" stroke-width="1.3" fill="none" stroke-linecap="round" stroke-linejoin="round"/>`,K_=()=>x`<svg viewBox="0 0 16 16" width="14" height="14" aria-hidden="true" focusable="false" class="pm"><circle cx="8" cy="8" r="6" fill="none" stroke-width="1.6" stroke-dasharray="3 2" /></svg>`,q_=(e,t)=>{let n=hd(e,t);return x`${n.head?`${n.head} `:w}<span class="sp-keep">${n.tail}, <span class="sp-warn"><span class="sp-pm">${K_()}</span>${n.first}</span></span>${n.rest?x`<span class="sp-warn"> ${n.rest}</span>`:w}`},J_=e=>x`<svg class="icon" viewBox=${oc[e].viewBox} width="18" height="18" aria-hidden="true" focusable="false">${s_(oc[e].body)}</svg>`,Y_=()=>x`<svg viewBox="0 0 24 24" width="12" height="12" aria-hidden="true" focusable="false"><path d="M6 9l6 6 6-6" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round" /></svg>`,X_=e=>e.replace(/\./g,`-`),Z_=e=>{let t=e.search(/[=–≡#]/);return t>0&&t<e.length-1?x`${e.slice(0,t+1)}<wbr />${e.slice(t+1)}`:e},Q_=class{element;#e;#t=null;#n=new Set;#r=null;#i=null;#a=null;#o=null;#s=null;#c=!1;#l;#u=null;constructor(e,t){this.#e=t,this.#l=new I_({now:()=>performance.now(),boxes:e=>this.#G(e),after:(e,t)=>this.#e.after(e,t),changed:()=>this.#P()});let n=e.createElement(`aside`);n.className=`sp`,n.setAttribute(`data-testid`,J),n.setAttribute(`aria-labelledby`,`sp-title`),n.hidden=!0,this.element=n,n.addEventListener(`keydown`,this.#j),n.addEventListener(`pointerdown`,this.#b),n.addEventListener(`focusout`,this.#x)}set place(e){e!==this.#d&&(this.#d=e,this.#k(!1),this.#P())}#d=null;get open(){return this.#c}set open(e){if(e!==this.#c){if(this.#c=e,this.element.hidden=!e,!e){this.#k(!1),this.#l.gone(),this.#u=null,this.#a&&this.#e.mark(null),this.#a=null,T(w,this.element);return}this.refresh()}}refresh(){if(this.#c){if(this.#t=this.#e.read(),this.#h(),this.#o&&!this.#t?.rows.some(e=>e.key===this.#o?.key))this.#k(!1);else if(this.#o&&this.#t){let e=this.#o.at;this.#o={key:this.#o.key,menu:fd(this.#t.info,this.#t.rows,this.#n,this.#o.key,this.#e.strings()),at:e}}this.#P()}}selectionChanged(){this.#c&&this.#t&&this.#f(!0)}hot(e){let t=e&&this.#t?.rows.some(t=>t.key===`${e.kind}.${e.index}`)?`${e.kind}.${e.index}`:null;t!==this.#i&&(this.#i=t,this.#c&&this.#f(!1))}#f(e){this.#m||=e,!this.#p&&(this.#p=this.#e.later(()=>{this.#p=null;let e=this.#m;this.#m=!1,this.#c&&(e?this.refresh():this.#P())}))}#p=null;#m=!1;focusFirst(){if(!this.#c)return;let e=this.#t?.rows??[],t=e.find(e=>this.#n.has(e.key))??e[0];if(t){this.#r=t.key,this.#P(),this.#_(t.key)?.focus({preventScroll:!1});return}this.element.querySelector(`[data-testid="${J}.close"]`)?.focus()}holdsFocus(){let e=this.element.getRootNode(),t=`activeElement`in e?e.activeElement:null;return!!(t&&this.element.contains(t))}dismiss(){this.#k(!1)}pointerMoved(e,t){this.#l.moved({x:e,y:t})}pointerGone(){this.#l.gone()}escape(){return this.#l.open!==null&&(this.#l.gone(),!0)}destroy(){this.#l.destroy(),this.#p?.(),this.#p=null,this.element.removeEventListener(`keydown`,this.#j),this.element.removeEventListener(`pointerdown`,this.#b),this.element.removeEventListener(`focusout`,this.#x),T(w,this.element),this.#t=null}#h(){let e=this.#t?.rows??[],t=od(e,this.#e.selected()),n=t.size===this.#n.size&&[...t].every(e=>this.#n.has(e));return this.#n=t,this.#r&&!e.some(e=>e.key===this.#r)&&(this.#r=null),!n}#g(){return this.#t?.rows??[]}#_(e){return this.element.querySelector(`[data-testid="${J}.row.${e}"]`)}#v(e){let t=this.#g(),n=sd(t,e);this.#e.select(n),this.#h(),this.#P();let r=t.filter(t=>e.has(t.key)).flatMap(e=>e.atoms);r.length&&this.#e.reveal(r)}#y(e){let t=(e instanceof Element?e.closest(`[data-row]`):null)?.getAttribute(`data-row`);return t?this.#g().find(e=>e.key===t)??null:null}#b=e=>{this.#l.pressed(),this.#o&&!(e.target instanceof Element&&e.target.closest(`.sp-menu, .grp`))&&this.#k(!1)};#x=e=>{let t=e.relatedTarget;this.#o&&!(t&&this.element.querySelector(`.sp-menu`)?.contains(t))&&!(t instanceof Element&&t.closest(`.grp`))&&this.#k(!1)};#S(e,t){this.#r=e.key,this.#v(cd(this.#n,e.key,t.shiftKey||t.ctrlKey||t.metaKey)),this.#_(e.key)?.focus({preventScroll:!0})}#C(e,t,n){n.stopPropagation(),this.#r=e.key,!this.#t?.readOnly&&(e.segments.find(e=>e.id===t)?.pressed||this.#e.config(e,t))}#w(e,t){if(t.stopPropagation(),this.#r=e.key,e.chip&&e.chip.disabled===null&&!this.#t?.readOnly){if(this.#o?.key===e.key){this.#k(!0);return}this.#T(e)}}#T(e){let t=this.#t,n=this.element.querySelector(`[data-testid="${J}.row.${e.key}.group"]`);if(!t||!n)return;let r=fd(t.info,t.rows,this.#n,e.key,this.#e.strings()),i=this.element.getBoundingClientRect(),a=n.getBoundingClientRect();this.#o={key:e.key,menu:r,at:{left:a.right-i.left,top:a.bottom-i.top}},this.#P(),this.#s?.(),this.#s=yr(n,e=>this.#e.frame(e),()=>this.#k(this.#D()));let o=this.#O();(o.find(e=>e.getAttribute(`aria-checked`)===`true`)??o[0])?.focus({preventScroll:!0})}#E(e){let t=this.element.querySelector(`.sp-menu`);if(!t)return;En(t),t.style.maxHeight=``;let n=t.ownerDocument.defaultView,r=Wt(t),i={width:n?.innerWidth??r.width,height:n?.innerHeight??r.height},a=Yr(e,r,i,this.element.getBoundingClientRect());if(hr(t,a.x,a.y),a.maxHeight!==null){let e=getComputedStyle(t),n=e=>parseFloat(e)||0,r=n(e.paddingTop)+n(e.paddingBottom)+n(e.borderTopWidth)+n(e.borderBottomWidth);t.style.maxHeight=`${Math.floor(a.maxHeight/nr(t)-r)}px`}}#D(){let e=this.element.getRootNode(),t=`activeElement`in e?e.activeElement:null;return!!(t&&this.element.querySelector(`.sp-menu`)?.contains(t))}#O(){return[...this.element.querySelectorAll(`.sp-menu [role^="menuitem"]`)]}#k(e){let t=this.#o?.key;this.#o&&(this.#s?.(),this.#s=null,this.#o=null,this.#P(),e&&t&&this.element.querySelector(`[data-testid="${J}.row.${t}.group"]`)?.focus({preventScroll:!0}))}#A(e){let t=this.#o,n=this.#t;if(!t||!n)return;let r=t.menu.items.find(t=>t.id===e);if(!r)return;let i=ud(n.rows,this.#n,t.key);if(this.#k(!0),!r.choice){this.#v(new Set(t.menu.selectKeys));return}let a=i.filter(e=>e.kind===`centre`).map(e=>e.atom),o=i.filter(e=>e.kind===`axis`).map(e=>e.bond);(a.length||o.length)&&this.#e.group(a,o,r.choice)}#j=e=>{if(e.defaultPrevented)return;let t=e.composedPath()[0];if(t?.closest?.(`.sp-menu`)){this.#M(e,t);return}if(this.#e.history(e)){e.preventDefault();return}if(e.key===`Escape`){e.preventDefault(),this.#l.gone(),this.#o?this.#k(!0):this.#e.leave();return}let n=this.#y(t??null);if(!n||!t)return;let r=this.#g(),i=this.#_(n.key),a=t===i,o=e.ctrlKey||e.metaKey;if(o&&!e.altKey&&e.key.toLowerCase()===`a`){e.preventDefault(),this.#v(new Set(r.map(e=>e.key)));return}if(e.key===`ArrowDown`||e.key===`ArrowUp`||e.key===`Home`||e.key===`End`){e.preventDefault();let t=r[vd(r.findIndex(e=>e.key===n.key),r.length,e.key)];if(!t)return;this.#r=t.key;let i=e.shiftKey?new Set([...this.#n,n.key,t.key]):new Set([t.key]);this.#v(i),this.#_(t.key)?.focus();return}if(e.key===`ArrowRight`||e.key===`ArrowLeft`){e.preventDefault();let n=[...i?.querySelectorAll(`button`)??[]],r=a?-1:n.indexOf(t),o=e.key===`ArrowRight`?Math.min(n.length-1,r+1):r-1;o<0?i?.focus():n[o]?.focus();return}if(a&&(e.key===`Enter`||e.key===` `)){e.preventDefault(),this.#v(cd(this.#n,n.key,o));return}if(!o&&!e.altKey&&e.key.length===1){let t=_d(n,e.key);t&&(e.preventDefault(),!this.#t?.readOnly&&!n.segments.find(e=>e.id===t)?.pressed&&this.#e.config(n,t))}};#M(e,t){let n=this.#O(),r=n.indexOf(t);if(e.key===`Escape`){e.preventDefault(),this.#k(!0);return}if(e.key===`Tab`){e.preventDefault(),this.#k(!0);return}let i=e.key===`ArrowDown`?(r+1)%n.length:e.key===`ArrowUp`?(r-1+n.length)%n.length:e.key===`Home`?0:e.key===`End`?n.length-1:-1;if(i>=0){e.preventDefault(),n[i]?.focus();return}(e.key===`Enter`||e.key===` `)&&(e.preventDefault(),t.click())}#N(e){let t=e?.key??null;t!==this.#a&&(this.#a=t,this.#e.mark(e?sd(this.#g(),new Set([e.key])):null))}#P(){if(!this.#c)return;let e=this.#e.strings(),t=this.#t,n=t?.rows??[],r=this.#r&&n.some(e=>e.key===this.#r)?this.#r:n.find(e=>this.#n.has(e.key))?.key??n[0]?.key??null;this.#r=r;let i=t&&!t.query&&n.length?gd(t.info,n,t.count,e):null,a=!!t?.readOnly,o=t?.acts??{remove:!1,rewedge:!1},s=!t||t.query||!n.length?t?.query?e.noStereoForQueries:e.noStereoElements:null,c=(t,n,r,i,o,s,c)=>x`<button
        type="button"
        class="ibtn sp-cmd"
        data-testid="${J}.${t}"
        aria-label=${r}
        aria-disabled=${o&&!a?w:`true`}
        aria-describedby=${o?V_(t):`sp-why-${t}`}
        @click=${()=>{o&&!a&&c()}}
        @pointerenter=${e=>this.#I(t,e)}
        @pointerleave=${e=>this.#L(t,e)}
        @focus=${e=>this.#R(t,e)}
        @blur=${()=>this.#z(t)}
      >
        ${J_(n)}</button
      ><span class="sp-why" id="sp-why-${t}" hidden>${s}</span>${this.#B(t,e.format(`commandTip`,{name:r,tip:o?i:s}))}`,l=x`${c(`remove-stereo`,`removeStereo`,e.removeStereoButton,e.removeStereoTip,o.remove,e.noStereoToRemove,()=>this.#e.command(`remove`))}${c(`recompute-wedges`,`recomputeWedges`,e.recomputeWedgesButton,e.recomputeWedgesTip,o.rewedge,e.noWedgesToPlace,()=>this.#e.command(`rewedge`))}`;if(T(x`<div class="sp-head">
          <h2 id="sp-title">${e.stereochemistry}</h2>
          <div class="sp-cmds">${l}</div>
          <button
            type="button"
            class="ibtn sp-close"
            data-testid="${J}.close"
            aria-label=${e.closeStereoPanel}
            aria-describedby=${V_(`close`)}
            @click=${()=>this.#e.close()}
            @pointerenter=${e=>this.#I(`close`,e)}
            @pointerleave=${e=>this.#L(`close`,e)}
            @focus=${e=>this.#R(`close`,e)}
            @blur=${()=>this.#z(`close`)}
          >
            <svg viewBox="0 0 24 24" width="16" height="16" aria-hidden="true" focusable="false"><path d="M6 6l12 12M18 6L6 18" stroke="currentColor" stroke-width="2" stroke-linecap="round" /></svg></button
          >${this.#B(`close`,e.close)}
        </div>
        ${i?x`<div class="sp-sum" data-testid="${J}.summary">
              <span class="sp-count">${i.unspecified?q_(i.elements,i.unspecified):i.elements}</span
              >${i.count||i.flag?x`<span class="muted"
                    >${i.count?x`<span data-testid="${J}.count">${i.count}</span>`:w}${i.count&&i.flag?` · `:w}${i.flag?x`<span data-testid="${J}.flag">${i.flag}</span>`:w}</span
                  >`:w}
            </div>`:w}
        ${s?x`<p class="sp-empty" data-testid="${J}.empty">${s}</p>`:x`<div class="sp-list" role="grid" aria-label=${e.stereoList} aria-multiselectable="true" @pointerleave=${()=>this.#N(null)}>
              ${en(n,e=>e.key,t=>this.#X(t,t.key===r,a,e))}
            </div>`}
        ${this.#o?this.#Q(this.#o,e):w}${this.#V(this.#K(s?[]:n,e))}${this.#q(e)}`,this.element),this.#o){let e=this.element.querySelector(`.sp-menu`);if(e&&!Zn(e)){let e=this.element.querySelector(`[data-testid="${J}.row.${this.#o.key}.group"]`);e&&this.#E(e.getBoundingClientRect())}}this.#W()}#F(){return this.#l.open??this.#u}#I(e,t){q(t)||this.#l.enter(e)}#L(e,t){q(t)||this.#l.leave(e,t.pointerType===`touch`?null:{x:t.clientX,y:t.clientY})}#R(e,t){let n=t.target.matches(`:focus-visible`)?e:null;n!==this.#u&&(this.#u=n,this.#P())}#z(e){this.#u===e&&(this.#u=null,this.#P())}#B(e,t){let n=this.#F()===e;return x`<span class="tip sp-tip" popover="manual" role="tooltip" id=${V_(e)} data-testid=${n?`${J}.${e}.tooltip`:w}>${t}</span>`}#V(e){let t=this.#F(),n=t&&!H_(t)?e.get(t):void 0;return x`<span class="tip sp-tip sp-row-tip" popover="manual" role="tooltip" id=${U_} data-testid=${n?`${t}.tooltip`:w}
      >${n?(e=>e.inline?e.meaning?this.#e.strings().format(`commandTip`,{name:e.name,tip:e.meaning}):e.name:x`<span class="tip-name">${e.name}</span>${e.meaning?x`<span class="tip-meaning">${e.meaning}</span>`:w}`)(n):w}</span
    >`}#H(e){return this.element.querySelector(`[data-testid="${H_(e)?`${J}.${e}`:e}"]`)}#U(e){return this.element.querySelector(`#${H_(e)?V_(e):U_}`)}#W(){let e=this.#c?this.#F():null,t=e!==null&&!H_(e)?e:null;for(let n of[...B_,U_]){let r=this.element.querySelector(`#${n===U_?U_:V_(n)}`);if(!r)continue;let i=n===U_?t:n===e?n:null,a=i===null?null:this.#H(i);if(!a||n===U_&&!r.textContent?.trim()){wr(r);continue}En(r);let o=a.getBoundingClientRect(),s=Wt(r),c=r.ownerDocument.defaultView,[l,u]=[c?.innerWidth??1/0,c?.innerHeight??1/0],d=o.bottom+R_;d+s.height>u-z_&&(d=o.top-R_-s.height);let f=Math.max(z_,Math.min(o.left+o.width/2-s.width/2,l-z_-s.width));hr(r,f,Math.max(z_,d))}}#G(e){let t=this.#U(e),n=this.#H(e);if(!t||!n||!Zn(t))return null;let[r,i]=[n.getBoundingClientRect(),t.getBoundingClientRect()];return{button:{left:r.left,top:r.top,right:r.right,bottom:r.bottom},tip:{left:i.left,top:i.top,right:i.right,bottom:i.bottom}}}#K(e,t){let n=new Map;for(let r of e){for(let e of r.segments)n.set(`${J}.row.${r.key}.config.${e.id}`,{name:e.name,meaning:Gu(e.meaning,t),inline:!0});r.chip&&r.chip.disabled===null&&r.chip.meaning&&n.set(`${J}.row.${r.key}.group`,{name:r.chip.name,meaning:Gu(r.chip.meaning,t),inline:!1})}return n}#q(e){return x`<div class="sp-meanings" hidden>
      ${Wu.map(t=>{let n=Gu(t,e);return n?x`<span id=${W_(t)}>${n}</span>`:w})}
    </div>`}#J(e,t){let n=e.filter(e=>e!==null&&Gu(e,t)!==``).map(W_);return n.length?n.join(` `):w}#Y(e){return{enter:t=>this.#I(e,t),leave:t=>this.#L(e,t),focus:t=>this.#R(e,t),blur:()=>this.#z(e)}}#X(e,t,n,r){let i=this.#n.has(e.key),a=`${J}.row.${e.key}`,o=e.chip?.disabled??null,s=`sp-why-${X_(e.key)}`;return x`<div
      class="sp-row ${this.#i===e.key?`hover`:``}"
      role="row"
      data-row=${e.key}
      data-testid=${a}
      aria-selected=${i?`true`:`false`}
      aria-label=${e.label}
      aria-describedby=${this.#J(e.meanings,r)}
      tabindex=${t?`0`:`-1`}
      @click=${t=>this.#S(e,t)}
      @pointerenter=${t=>{t.pointerType===`mouse`&&this.#N(e)}}
      @focusin=${()=>{this.#r!==e.key&&(this.#r=e.key,this.#P())}}
    >
      <span class="who" role="gridcell"
        ><span class="kind" aria-hidden="true"><svg viewBox="0 0 16 16" width="16" height="16" focusable="false">${G_(e.kind)}</svg></span
        ><span class="nm"><span class="name" title=${e.name}>${Z_(e.name)}</span><span class="kword">${e.kindWord}</span></span>${e.problem?x`<span class="pm-wrap" title=${r.problemStereo}>${K_()}</span>`:w}</span
      ><span role="gridcell" class="cfg"
        ><span class="seg" role="group" aria-label=${r.format(`stereoConfigOf`,{name:e.name})}
          >${e.segments.map(t=>{let i=this.#Y(`${a}.config.${t.id}`);return x`<button
              type="button"
              tabindex="-1"
              class=${t.unset?`unset`:``}
              data-testid="${a}.config.${t.id}"
              aria-pressed=${t.pressed?`true`:`false`}
              aria-disabled=${n?`true`:w}
              aria-label=${t.name}
              aria-describedby=${this.#J([t.meaning],r)}
              @click=${n=>this.#C(e,t.id,n)}
              @pointerenter=${i.enter}
              @pointerleave=${i.leave}
              @focus=${i.focus}
              @blur=${i.blur}
            >
              ${t.face}
            </button>`})}</span
        ></span
      ><span role="gridcell" class="gcell"
        >${e.chip?this.#Z(e,a,o,s,n,r):x`<span class="grp-empty" aria-hidden="true"></span>`}</span
      >
    </div>`}#Z(e,t,n,r,i,a){let o=e.chip,s=this.#Y(`${t}.group`),c=n===null;return x`<button
        type="button"
        tabindex="-1"
        class="grp"
        data-testid="${t}.group"
        aria-haspopup="menu"
        aria-expanded=${this.#o?.key===e.key?`true`:`false`}
        aria-disabled=${!c||i?`true`:w}
        aria-describedby=${c?this.#J([o.meaning],a):r}
        aria-label=${o.name}
        title=${c?w:n}
        @click=${t=>this.#w(e,t)}
        @pointerenter=${c?s.enter:w}
        @pointerleave=${c?s.leave:w}
        @focus=${c?s.focus:w}
        @blur=${c?s.blur:w}
      >
        <span class=${o.group?``:`none`}>${o.face}</span>${Y_()}</button
      >${c?w:x`<span class="sp-why" id=${r} hidden>${n}</span>`}`}#Q(e,t){let n=e.menu.items,r=n.filter(e=>e.checked!==null),i=n.filter(e=>e.checked===null),a=e=>{let n=e.meaning?Gu(e.meaning,t):``;return x`<button
        type="button"
        tabindex="-1"
        role=${e.checked===null?`menuitem`:`menuitemradio`}
        class="item ${n?`has-meaning`:``}"
        data-testid="${L_}.${e.id}"
        aria-checked=${e.checked===null?w:e.checked?`true`:`false`}
        aria-describedby=${this.#J([e.meaning],t)}
        @click=${()=>this.#A(e.id)}
      >
        <span class="mark" aria-hidden="true">${e.checked?`✓`:``}</span
        ><span class="item-body"
          ><span class="item-line"><span class="item-name">${e.label}</span>${e.hint?x`<span class="hint">${e.hint}</span>`:w}</span
          >${n?x`<span class="meaning" aria-hidden="true">${n}</span>`:w}</span
        >
      </button>`};return x`<div class="sp-menu" popover="manual" role="menu" data-testid=${L_} aria-label=${e.menu.heading}>
      <div class="heading" aria-hidden="true">${e.menu.heading}</div>
      ${r.map(a)}
      <div class="sep" role="separator"></div>
      ${i.map(a)}
    </div>`}},$_=`menu.context`,ev=e=>`ctx-${(e??`main`).replace(/[^a-z0-9-]/gi,`_`)}`,tv=e=>e.flatMap(e=>e.type===`item`?[e]:e.type===`group`?[...e.items]:[]),nv=e=>`ctx-tip-${e.replace(/[^a-z0-9-]/gi,`_`)}`,rv=Zn,iv=e=>!q(e)&&(e.pointerType===`mouse`||e.pointerType===`pen`&&e.buttons===0),av=class{#e;#t;#n=null;#r=null;#i=!1;#a=null;constructor(e,t){this.#e=e,this.#t=t,t.setAttribute(`data-menus`,``),t.addEventListener(`focusout`,this.#u),t.addEventListener(`contextmenu`,this.#l)}get#o(){return this.#n?.testid??$_}get isOpen(){return this.#n!==null}open(e){if(this.#n&&this.#s(!1),!this.#t.isConnected)return;this.#n={...e,subs:[]},this.#i=e.keyboard,this.#O();let t=this.#x(null);if(!t)return;this.#T(t,null);let n=this.#C(t)[0];e.keyboard&&n?this.#y(n):t.focus({preventScroll:!0});let r=this.#t.getRootNode().host,i=this.#e.frame?.bind(this.#e);r&&i&&this.#n&&(this.#a=yr(r,i,()=>this.#s(this.#c())))}close(e=!1){this.#n&&this.#s(e)}destroy(){this.#a?.(),this.#a=null,this.#t.removeEventListener(`focusout`,this.#u),this.#t.removeEventListener(`contextmenu`,this.#l),this.#n=null,T(w,this.#t)}#s(e){let t=this.#n;if(t){this.#a?.(),this.#a=null,this.#n=null,this.#r=null;for(let e of Array.from(this.#t.querySelectorAll(`[popover], [data-popover]`)))wr(e);T(w,this.#t),t.closed?t.closed():this.#e.closed(),e&&t.back?.focus({preventScroll:!0})}}#c(){let e=this.#t.getRootNode().activeElement;return!!(e&&this.#t.contains(e))}#l=e=>{e.preventDefault(),e.stopPropagation()};#u=e=>{if(!this.#n)return;let t=e.relatedTarget;t&&this.#t.contains(t)||this.#s(!1)};#d(e,t){if(!this.#n)return;if(e.disabled){e.reason&&!t&&this.#A(e);return}if(e.submenu){this.#h(e,this.#f(e.id),t);return}let n=this.#n;if(e.keepOpen&&n.toggle){n.entries=n.toggle(e.id),this.#O();let t=this.#x(null);t&&this.#T(t,null);for(let e of n.subs){let t=this.#x(e);t&&this.#T(t,e)}this.#S(e.id)?.focus({preventScroll:!0});return}this.#s(!0),n.chose?n.chose(e.id):this.#e.chose(e.id)}#f(e){let t=this.#n;if(!t)return 0;for(let n=t.subs.length;n>0;n--)if(this.#p(n,e))return n;return 0}#p(e,t){return tv(this.#m(e)).some(e=>e.id===t)}#m(e){let t=this.#n;if(!t)return[];let n=t.entries;for(let r=0;r<e;r++)n=tv(n).find(e=>e.id===t.subs[r])?.submenu??[];return n}#h(e,t,n){let r=this.#n;if(!r)return;if(r.subs[t]!==e.id||r.subs.length>t+1){this.#g(t),r.subs=[...r.subs.slice(0,t),e.id],this.#O();let n=this.#x(e.id);n&&this.#T(n,e.id)}let i=n?this.#C(this.#x(e.id))[0]:void 0;i&&this.#y(i)}#g(e){let t=this.#n;if(t)for(let n of t.subs.slice(e)){let e=this.#x(n);e&&wr(e)}}#_(e,t){let n=this.#n;if(!n||n.subs.length<e||e<1)return;let r=n.subs[e-1],i=n.subs.slice(e-1).map(e=>this.#x(e)),a=this.#t.getRootNode().activeElement;if(t||a&&i.some(e=>e?.contains(a))){let e=this.#S(r);e?this.#y(e):this.#x(null)?.focus({preventScroll:!0})}this.#g(e-1),n.subs=n.subs.slice(0,e-1),this.#O()}#v(e,t){if(e.stopPropagation(),this.#i=!0,!this.#n)return;let n=e.currentTarget,r=this.#C(n),i=r.indexOf(e.composedPath()[0]),a=()=>e.preventDefault();switch(e.key){case`ArrowDown`:case`ArrowUp`:{if(a(),!r.length)return;let t=e.key===`ArrowDown`?1:-1,n=i<0?t>0?0:r.length-1:(i+t+r.length)%r.length;this.#y(r[n]);return}case`Home`:case`End`:a(),r.length&&this.#y(r[e.key===`Home`?0:r.length-1]);return;case`Enter`:case` `:{a();let e=i>=0?this.#w(r[i]):null;e&&this.#d(e,!0);return}case`ArrowRight`:{a();let e=i>=0?this.#w(r[i]):null;e?.submenu&&!e.disabled&&this.#h(e,t,!0);return}case`ArrowLeft`:a(),t>0&&this.#_(t,!0);return;case`Escape`:a(),t>0?this.#_(t,!0):this.#s(!0);return;case`Tab`:a(),this.#s(!0);return}}#y(e){e.focus({preventScroll:!0}),this.#b(e)}#b(e){let t=e.closest(`[role="menu"]`);if(!t||t.scrollHeight<=t.clientHeight)return;let n=t.getBoundingClientRect(),r=getComputedStyle(t),[i,a]=[parseFloat(r.paddingTop)||0,parseFloat(r.paddingBottom)||0],o={top:n.top+t.clientTop+i,height:t.clientHeight-i-a};t.scrollTop=_(t.scrollTop,o,e.getBoundingClientRect())}#x(e){return this.#t.querySelector(`#${ev(e)}`)}#S(e){return this.#t.querySelector(`[data-testid="${this.#o}.${e}"]`)}#C(e){return Array.from(e?.querySelectorAll(`[role="menuitem"], [role="menuitemradio"], [role="menuitemcheckbox"]`)??[])}#w(e){let t=e.dataset.item,n=this.#n;return!t||!n?null:xn(n.entries).find(e=>e.id===t)??null}#T(e,t){let n=this.#n;if(!n||!e.isConnected)return;En(e),e.style.maxHeight=``;let r=e.ownerDocument.defaultView,i=Wt(e),a={width:r?.innerWidth??i.width,height:r?.innerHeight??i.height},o=t?this.#S(t):null,s=o?.getBoundingClientRect()??null,c=o?.closest(`[role="menu"]`),l=this.#x(null),u=c&&l&&c!==l&&c.getBoundingClientRect().left<l.getBoundingClientRect().left?`left`:`right`,d=e.querySelector(`[aria-haspopup="menu"]`)?i.width:0,p=s?Wr(s,i,a,n.box,void 0,u,d):Yt(n.at,i,a,n.box),m=t?n.subs.indexOf(t):-1;if(o&&m>=1){let e=this.#E(n,m),t=at(o.getBoundingClientRect(),i,a,n.box,e,void 0,u);if(!t&&c&&this.#D(n,m,c,a)){let r=c.getBoundingClientRect().left<(l?.getBoundingClientRect().left??0)?`left`:`right`;t=at(o.getBoundingClientRect(),i,a,n.box,e,void 0,r),p=t??Wr(o.getBoundingClientRect(),i,a,n.box,void 0,r)}t??=f(o.getBoundingClientRect(),i,a,n.box,e),t&&(p=t)}if(hr(e,p.x,p.y),p.maxHeight!==null){let t=getComputedStyle(e),n=e=>parseFloat(e)||0,r=n(t.paddingTop)+n(t.paddingBottom)+n(t.borderTopWidth)+n(t.borderBottomWidth);e.style.maxHeight=`${Math.floor(p.maxHeight/nr(e)-r)}px`}e.tabIndex=e.scrollHeight>e.clientHeight?0:-1}#E(e,t){return t<1?[]:[null,...e.subs.slice(0,t-1)].flatMap(e=>{let t=this.#x(e);return t&&rv(t)?[t.getBoundingClientRect()]:[]})}#D(e,t,n,r){let i=this.#S(e.subs[t-1]);if(!i)return!1;let a=i.getBoundingClientRect(),o=n.getBoundingClientRect().left<a.left?`right`:`left`,s=Wt(n),c=at(a,s,r,e.box,this.#E(e,t-1),void 0,o,!0);return!c||c.x<a.left!=(o===`left`)?!1:(hr(n,c.x,c.y),!0)}#O(){let e=this.#n;if(!e){T(w,this.#t);return}let t=[],n=e.entries;e.subs.forEach((e,r)=>{let i=tv(n).find(t=>t.id===e);i?.submenu&&(t.push(this.#M(i.submenu,i.label,i.id,null,r+1)),n=i.submenu)});let r=xn(e.entries).filter(e=>e.disabled&&e.reason);T(x`${this.#M(e.entries,e.label,null,e.target,0)}${t}${r.map(e=>this.#k(e))}`,this.#t)}#k(e){let t=this.#r===e.id;return x`<span class="tip" popover="manual" role="tooltip" id=${nv(e.id)} data-testid=${t?`${this.#o}.${e.id}.tooltip`:w}
      >${e.reason}</span
    >`}#A(e){if(!this.#n||!e.reason)return;this.#r=e.id,this.#O();let t=this.#t.querySelector(`#${nv(e.id)}`),n=this.#S(e.id);if(!t||!n)return;En(t);let r=n.getBoundingClientRect(),i=t.ownerDocument.defaultView,a=Wt(t),o=i?.innerWidth??r.right+a.width,s=i?.innerHeight??r.bottom+a.height,c=r.right+4;c+a.width>o-4&&(c=r.left-4-a.width);let l=r.top+(r.height-a.height)/2;c<4&&(c=Math.max(4,Math.min(r.left,o-4-a.width)),l=r.bottom+4),hr(t,c,Math.max(4,Math.min(l,s-4-a.height)))}#j(e=null){if(this.#r===null||e!==null&&this.#r!==e)return;let t=this.#t.querySelector(`#${nv(this.#r)}`);t&&wr(t),this.#r=null,this.#n&&this.#O()}#M(e,t,n,r,i){return x`<div
      class="palette menu ctx"
      popover="manual"
      role="menu"
      tabindex="-1"
      id=${ev(n)}
      aria-label=${t}
      data-testid=${n?`${this.#o}.${n}.submenu`:this.#o}
      data-target=${r??w}
      @keydown=${e=>this.#v(e,i)}
    >
      ${e.map(e=>this.#N(e,i))}
    </div>`}#N(e,t){return e.type===`separator`?x`<div class="sep" role="separator"></div>`:e.type===`group`?x`<div class="grp" role="group" aria-label=${e.label}>
        ${e.heading?x`<div class="heading" aria-hidden="true">${e.label}</div>`:w}${e.items.map(e=>this.#F(e,t))}
      </div>`:this.#F(e,t)}#P(e,t){e.disabled&&e.reason&&this.#A(e),e.submenu&&!e.disabled?this.#h(e,t,!1):(this.#n?.subs.length??0)>t&&this.#_(t+1,!1)}#F(e,t){let n=this.#n?.subs[t]===e.id,r=e.icon?oc[e.icon]:null;return x`<button
      type="button"
      class="item"
      role=${e.radio?`menuitemradio`:e.checkbox?`menuitemcheckbox`:`menuitem`}
      tabindex="-1"
      data-testid="${this.#o}.${e.id}"
      data-item=${e.id}
      aria-label=${e.label}
      aria-checked=${e.radio||e.checkbox?String(!!e.checked):w}
      aria-disabled=${e.disabled?`true`:w}
      aria-describedby=${e.disabled&&e.reason?nv(e.id):w}
      aria-haspopup=${e.submenu?`menu`:w}
      aria-expanded=${e.submenu?String(n):w}
      aria-controls=${e.submenu&&n?ev(e.id):w}
      aria-keyshortcuts=${e.keys?.aria??w}
      @click=${t=>this.#d(e,t.detail===0)}
      @pointerenter=${n=>{!iv(n)||this.#i||this.#P(e,t)}}
      @pointermove=${n=>{this.#i&&iv(n)&&(n.movementX!==0||n.movementY!==0)&&(this.#i=!1,this.#P(e,t))}}
      @pointerleave=${t=>{q(t)||this.#j(e.id)}}
      @focus=${()=>e.disabled&&e.reason?this.#A(e):this.#j()}
      @blur=${()=>this.#j(e.id)}
    >
      <span class="mark" aria-hidden="true">${(e.radio||e.checkbox)&&e.checked?`✓`:``}</span>
      ${r?x`<svg class="icon" viewBox=${r.viewBox} aria-hidden="true" focusable="false">${s_(r.body)}</svg>`:x`<span class="icon" aria-hidden="true"></span>`}
      <span class="item-name">${e.label}</span>
      ${e.keys?x`<kbd aria-hidden="true">${e.keys.label}</kbd>`:w}
      ${e.submenu?x`<span class="chev" aria-hidden="true">▸</span>`:w}
    </button>`}},ov=`http://www.w3.org/2000/svg`,sv={formula:`panel.info.formula`,"molecular-weight":`panel.info.weight`,"exact-mass":`panel.info.exact-mass`},cv=30,lv=class{element;#e;#t;#n;#r;#i=[];#a=``;#o=`shown`;#s=ch;#c=p;#l=null;#u=null;#d=[];#f=0;#p=0;#m=[];#h=!0;#g=!1;#_;#v=null;#y=!1;#b=!1;constructor(e,t){this.#_=t;let n=e.createElement(`div`);n.className=`info`,n.setAttribute(`data-testid`,`panel.info`),n.setAttribute(`role`,`group`);let r=e.createElement(`div`);r.className=`info-body`,r.id=`crux-info-body`,r.setAttribute(`data-testid`,`panel.info.values`);let i=e.createElement(`span`);i.className=`info-caption`,i.hidden=!0,r.append(i);let a=e.createElement(`button`);a.type=`button`,a.className=`info-toggle`,a.setAttribute(`data-testid`,`panel.info.toggle`),a.setAttribute(`aria-controls`,r.id);let o=e.createElementNS(ov,`svg`);o.setAttribute(`class`,`info-chevron`),o.setAttribute(`viewBox`,oc.open.viewBox),o.setAttribute(`aria-hidden`,`true`),o.setAttribute(`focusable`,`false`),o.innerHTML=oc.open.body,a.append(o);let s=e.createElement(`button`);s.type=`button`,s.className=`info-more`,s.tabIndex=-1,s.hidden=!0,s.setAttribute(`data-testid`,`panel.info.more`),s.setAttribute(`aria-haspopup`,`menu`),s.setAttribute(`aria-expanded`,`false`),s.textContent=`…`,n.append(a,r,s),a.addEventListener(`click`,this.#x),s.addEventListener(`click`,this.#A),n.addEventListener(`click`,this.#k),n.addEventListener(`keydown`,this.#M),n.addEventListener(`contextmenu`,this.#N),n.addEventListener(`pointerdown`,this.#P),n.addEventListener(`pointermove`,this.#F),n.addEventListener(`pointerup`,this.#I),n.addEventListener(`pointercancel`,this.#I),this.element=n,this.#e=a,this.#t=r,this.#n=i,this.#r=s,this.relabel(p),this.#E()}relabel(e){this.#c=e,this.element.setAttribute(`aria-label`,e.infoPanel),this.#e.setAttribute(`aria-label`,e.infoPanel),this.#n.textContent=e.infoSelection,this.#r.setAttribute(`aria-label`,e.readoutMore),this.#r.title=e.readoutMore,this.#C()}#x=()=>{let e=this.#_;if(!e)return;let t=Km(this.#o,this.#h,this.#g);t.state!==this.#o&&e.toggled(t.state),this.#g=t.opened,this.#E()};destroy(){this.#e.removeEventListener(`click`,this.#x),this.#r.removeEventListener(`click`,this.#A);let e=this.element;e.removeEventListener(`click`,this.#k),e.removeEventListener(`keydown`,this.#M),e.removeEventListener(`contextmenu`,this.#N),e.removeEventListener(`pointerdown`,this.#P),e.removeEventListener(`pointermove`,this.#F),e.removeEventListener(`pointerup`,this.#I),e.removeEventListener(`pointercancel`,this.#I),this.#v?.cancel(),this.#v=null,this.#_=null,this.#l=null}set state(e){e!==this.#o&&(this.#g=!1),this.#o=e,this.#C()}get state(){return this.#o}fit(e){if(!(e.width>0&&e.height>0))return;this.#u={width:e.width,height:e.height},this.#w();let t=this.#T();this.#g=qm(this.#g,this.#h,t),this.#h=t,this.#E()}get text(){return Sh(this.#s,this.#c)}get line(){return this.#s}get folded(){if(this.#r.hidden)return[];let e=new Set(this.#m);return this.#i.filter((t,n)=>!e.has(n)).map(e=>e.item)}show(e){let t=JSON.stringify(e);t!==this.#a&&(this.#a=t,this.#s=e,this.#S(e),this.#C())}#S(e){let t=this.element.ownerDocument,n=e.kind===`items`?e.items:[];this.#n.hidden=!(e.kind===`items`&&e.selection);let r=this.element.getRootNode().activeElement,i=this.#i.some(e=>e.el===r);this.#i=n.map(e=>({item:e,el:dv(t,e)})),this.#t.replaceChildren(this.#n,...this.#i.flatMap(e=>[t.createTextNode(` `),e.el])),i&&this.#e.focus({preventScroll:!0}),this.element.toggleAttribute(`data-empty`,e.kind===`empty`)}#C(){this.#w(),this.#h=this.#T(),this.#E()}#w(){if(this.#o!==`shown`)return;this.element.hidden&&(this.element.hidden=!1),this.#t.hidden&&(this.#t.hidden=!1);for(let e of this.#i)e.el.hidden=!1;if(!this.#t.getClientRects().length)return;let e=4+(parseFloat(getComputedStyle(this.#t).wordSpacing)||0);this.#p=this.#n.hidden?0:this.#n.getBoundingClientRect().width+e,this.#d=this.#i.map(t=>Math.ceil(t.el.getBoundingClientRect().width)+e),this.#f=this.#t.scrollWidth}#T(){let e=this.#d.length,t=[...Array(e).keys()];if(!this.#u||!e)return this.#m=t,Wm(40+this.#p+uv(this.#d),this.#u);if(this.#u.height<200)return this.#m=[],!1;let n=Hm(this.#u.width)-1-40-this.#p,r=this.#i.map(e=>!!e.item.pointer);this.#m=this.#f<=n+this.#p?t:Ch(this.#d,r,n,cv);let i=this.#m.length<e;return Wm(40+this.#p+(this.#m.length?uv(this.#m.map(e=>this.#d[e])):this.#d[0]??0)+(i?cv:0),this.#u)}#E(){let e=Gm(this.#o,this.#h,this.#g),t=Gm(this.#o,this.#h,this.#g,this.#s.kind===`empty`)===`hidden`,n=this.element.getRootNode().activeElement;t&&!this.#b&&n&&this.element.contains(n)&&this.#_?.away(),this.#b=t,this.element.hidden!==(e===`hidden`)&&(this.element.hidden=e===`hidden`),this.#t.hidden!==(this.#o!==`shown`)&&(this.#t.hidden=this.#o!==`shown`);let r=String(e===`open`);this.#e.getAttribute(`aria-expanded`)!==r&&this.#e.setAttribute(`aria-expanded`,r),this.element.toggleAttribute(`data-collapsed`,this.#o===`collapsed`),this.element.toggleAttribute(`data-folded`,this.#o===`shown`&&e===`folded`);let i=this.#i.length,a=e===`open`&&this.#u!==null&&this.#m.length>0&&this.#m.length<i,o=new Set(a?this.#m:this.#i.keys());this.#i.forEach((e,t)=>{e.el.hidden!==!o.has(t)&&(e.el.hidden=!o.has(t))}),this.#r.hidden=!a,a||this.#r.setAttribute(`aria-expanded`,`false`);let s=e===`open`&&this.#u!==null&&i>0&&this.#m.length===0&&!Um(40+this.#f,this.#u.width)?Sh(this.#s,this.#c):``;s?this.element.title!==s&&(this.element.title=s):this.element.removeAttribute(`title`),this.#D()}yieldTo(e){this.#l=e,this.#D()}#D(){let e=this.#l;for(let t of this.#i){let n=e!==null&&!t.el.hidden&&t.el.getClientRects().length>0&&e(t.el.getBoundingClientRect());t.el.hasAttribute(`data-yields`)!==n&&t.el.toggleAttribute(`data-yields`,n)}}set moreOpen(e){this.#r.setAttribute(`aria-expanded`,String(e))}#O(e){let t=e instanceof Element?e.closest(`[data-item]`):null;return t&&this.#i.find(e=>e.el===t)||null}#k=e=>{let t=this.#O(e.target);if(t&&this.#_){if(this.#y){this.#y=!1;return}t.item.pending||this.#_.copy(t.item)}};#A=e=>{this.#_?.more(this.folded,this.#r,e.detail===0)};#j(){return[this.#e,...this.#i.filter(e=>!e.el.hidden&&!this.#t.hidden).map(e=>e.el),...this.#r.hidden?[]:[this.#r]]}#M=e=>{let t=this.#_;if(!t)return;let n=e.composedPath()[0],r=this.#j(),i=r.indexOf(n);if(i<0)return;let a=()=>{e.preventDefault(),e.stopPropagation()};if(e.key===`ArrowRight`||e.key===`ArrowLeft`||e.key===`Home`||e.key===`End`){a(),r[e.key===`Home`?0:e.key===`End`?r.length-1:Math.max(0,Math.min(r.length-1,i+(e.key===`ArrowRight`?1:-1)))].focus({preventScroll:!0});return}if(e.key===`ContextMenu`||e.key===`F10`&&e.shiftKey){a(),t.menu(null,this.#O(n)?.item.id??null,n);return}let o=this.#O(n);o&&(e.key===`Enter`||e.key===` `)&&(a(),o.item.pending||t.copy(o.item))};#N=e=>{if(e.preventDefault(),e.stopPropagation(),!this.#_||q(e))return;let t=this.#O(e.target);this.#_.menu({x:e.clientX,y:e.clientY},t?.item.id??null,t?.el??this.#e)};#P=e=>{if(this.#v?.cancel(),this.#v=null,this.#y=!1,e.pointerType!==`touch`||q(e)||!this.#_)return;let t=this.#O(e.target);if(!t)return;e.preventDefault();let n={x:e.clientX,y:e.clientY},r=this.#_.after(500,()=>{this.#v=null,this.#y=!0,this.#_?.menu(n,t.item.id,t.el)});this.#v={id:e.pointerId,...n,cancel:r}};#F=e=>{let t=this.#v;!t||t.id!==e.pointerId||Math.hypot(e.clientX-t.x,e.clientY-t.y)<4||(t.cancel(),this.#v=null)};#I=e=>{let t=this.#v;t&&t.id===e.pointerId&&(t.cancel(),this.#v=null)}},uv=e=>e.reduce((e,t)=>e+t,0);function dv(e,t){let n=e.createElement(`span`);if(n.className=`info-item`,n.setAttribute(`role`,`button`),n.tabIndex=-1,n.setAttribute(`data-testid`,`panel.info.item.${t.id}`),n.setAttribute(`data-item`,t.id),t.pending&&n.setAttribute(`aria-busy`,`true`),t.title&&(n.title=t.title),t.label){let r=e.createElement(t.labelName?`abbr`:`span`);r.className=`info-label`,r.textContent=t.label,t.labelName&&(r.title=t.labelName),n.append(r,` `)}let r=e.createElement(`span`);r.className=t.runs?`info-formula`:t.id===`smiles`||t.id===`cxsmiles`||t.id===`smarts`||t.id===`inchikey`?`info-value info-notation`:`info-value`;let i=sv[t.id];if(i&&r.setAttribute(`data-testid`,i),t.runs?r.append(...fv(e,t.runs)):r.textContent=t.value,n.append(r),t.note){let r=e.createElement(`span`);r.className=`info-note`,r.textContent=t.note,n.append(` `,r)}return n}function fv(e,t){return t.map(t=>{if(t.kind===`text`)return e.createTextNode(t.text);let n=e.createElement(t.kind);return n.textContent=t.text,n})}var pv=6e3,mv=class{#e;#t;#n=null;#r=null;constructor(e,t){this.#e=e,this.#t=t}get text(){return this.#n??this.#r?.message??``}get closable(){return this.#n===null&&this.#r!==null}refusal(e){e!==this.#n&&(this.#n=e,this.#t())}notice(e,t=pv){this.#r?.cancel(),this.#r={message:e,cancel:this.#e(t,()=>this.dismiss())},this.#t()}dismiss(){let e=this.#r;return e?(e.cancel(),this.#r=null,this.#t(),!0):!1}destroy(){this.#r?.cancel(),this.#r=null,this.#n=null}},hv=`dialog.rgroup`,gv=`dialog.attachment-points`,_v=class{#e;#t;#n=null;#r=null;constructor(e,t){this.#e=e,this.#t=t}get open(){return this.#n!==null||this.#r!==null}openRGroup(e,t){if(this.open||!this.#t.isConnected)return;let n=`atomId`in e?this.#e.rNumbers(e.atomId):[];if(n===null)return;let r=wp(n),i=[...r.pressed].sort((e,t)=>e-t)[0]??1;this.#n={target:e,choice:r,stop:i,opener:t},this.#g(),this.#m(hv)?.showModal(),this.#i(i)?.focus()}#i(e){return this.#t.querySelector(`[data-testid="${hv}.r.${e}"]`)}#a(e){let t=this.#n;t&&(t.choice=Tp(t.choice,e),t.stop=e,this.#g())}#o(e,t){let n=this.#n;n&&Op(e.key)&&(e.preventDefault(),n.stop=Dp(t,e.key),this.#g(),this.#i(n.stop)?.focus())}#s(){let e=this.#n;e&&(this.#l(),this.#e.applyRNumbers(e.target,e.choice.opened,Ep(e.choice)),this.#e.restoreFocus(e.opener))}#c(){let e=this.#n;e&&(this.#l(),this.#e.restoreFocus(e.opener))}#l(){this.#m(hv)?.close(),this.#n=null,this.#g()}openAttachmentPoints(e,t){if(this.open||!this.#t.isConnected)return;let n=this.#e.attachmentPoints(e);n!==null&&(this.#r={atomId:e,was:n,...Ip(n),opener:t},this.#g(),this.#m(gv)?.showModal(),this.#t.querySelector(`[data-testid="${gv}.primary"]`)?.focus())}#u(e,t){let n=this.#r;n&&(n[e]=t,this.#g())}#d(){let e=this.#r;e&&(this.#p(),this.#e.applyAttachmentPoints(e.atomId,e.was,Lp(e.primary,e.secondary)),this.#e.restoreFocus(e.opener))}#f(){let e=this.#r;e&&(this.#p(),this.#e.restoreFocus(e.opener))}#p(){this.#m(gv)?.close(),this.#r=null,this.#g()}cancelAll(){if(!this.open)return;let e=this.#n?.opener??this.#r?.opener??null;this.#m(hv)?.close(),this.#m(gv)?.close(),this.#n=null,this.#r=null,this.#g(),this.#e.restoreFocus(e)}#m(e){return this.#t.querySelector(`dialog[data-testid="${e}"]`)}#h(e,t){if(e.key===`Escape`){e.preventDefault(),t===`r`?this.#c():this.#f();return}if(e.key===`Tab`){kr(e);return}let n=e.composedPath()[0];t===`points`&&e.key===`Enter`&&n instanceof HTMLInputElement&&(e.preventDefault(),this.#d())}#g(){T(x`${this.#n?this.#_(this.#n):w}${this.#r?this.#v(this.#r):w}`,this.#t)}#_(e){return x`<dialog
      class="dialog"
      data-testid=${hv}
      aria-labelledby="rg-title"
      aria-describedby="rg-hint"
      @keydown=${e=>this.#h(e,`r`)}
      @cancel=${e=>{e.preventDefault(),this.#c()}}
    >
      <div class="dialog-head">
        <h2 id="rg-title">${this.#e.strings().rgroupDialog}</h2>
        <button type="button" class="close" data-testid="${hv}.close" aria-label=${this.#e.strings().close} @click=${()=>this.#c()}>×</button>
      </div>
      <p class="rg-hint" id="rg-hint">${this.#e.strings().rgroupHint}</p>
      <div class="rgrid" role="group" aria-label=${this.#e.strings().rgroupNumbers}>
        ${Sp.map(t=>x`<button
            type="button"
            class="r-toggle"
            data-testid="${hv}.r.${t}"
            aria-pressed=${String(e.choice.pressed.has(t))}
            tabindex=${e.stop===t?`0`:`-1`}
            @click=${()=>this.#a(t)}
            @keydown=${e=>this.#o(e,t)}
            @focus=${()=>{this.#n&&this.#n.stop!==t&&(this.#n.stop=t,this.#g())}}
          >
            ${this.#e.strings().rName(t)}
          </button>`)}
      </div>
      <div class="dialog-foot">
        <button type="button" class="primary" data-testid="${hv}.ok" @click=${()=>this.#s()}>${this.#e.strings().ok}</button>
        <button type="button" data-testid="${hv}.cancel" @click=${()=>this.#c()}>${this.#e.strings().cancel}</button>
      </div>
    </dialog>`}#v(e){let t=(t,n)=>x`<label class="check">
      <input
        type="checkbox"
        data-testid="${gv}.${t}"
        .checked=${e[t]}
        @change=${e=>this.#u(t,e.target.checked)}
      />
      ${n}
    </label>`;return x`<dialog
      class="dialog"
      data-testid=${gv}
      aria-labelledby="apt-title"
      @keydown=${e=>this.#h(e,`points`)}
      @cancel=${e=>{e.preventDefault(),this.#f()}}
    >
      <div class="dialog-head">
        <h2 id="apt-title">${this.#e.strings().attachmentDialog}</h2>
        <button type="button" class="close" data-testid="${gv}.close" aria-label=${this.#e.strings().close} @click=${()=>this.#f()}>×</button>
      </div>
      <div class="checks">${t(`primary`,this.#e.strings().primaryPoint)} ${t(`secondary`,this.#e.strings().secondaryPoint)}</div>
      <div class="dialog-foot">
        <button type="button" class="primary" data-testid="${gv}.ok" @click=${()=>this.#d()}>${this.#e.strings().ok}</button>
        <button type="button" data-testid="${gv}.cancel" @click=${()=>this.#f()}>${this.#e.strings().cancel}</button>
      </div>
    </dialog>`}},Y=`dialog.sgroup`,vv=class{#e;#t;#n=null;constructor(e,t){this.#e=e,this.#t=t}get open(){return this.#n!==null}openFor(e,t){if(this.open||!this.#t.isConnected)return;let n=this.#e.opening(e);n&&(this.#n={request:e,form:n.form,editing:n.editing,refused:null,told:!1,opener:t},this.#u(),this.#c()?.showModal(),this.#r()?.focus())}#r(){let e=this.#n;return e?this.#t.querySelector(`[data-testid="${Y}.${e.editing?{SUP:`name`,DAT:`field-name`,MUL:`count`,SRU:`label`}[e.form.type]:`type`}"]`):null}#i(e){let t=this.#n;t&&(t.form={...t.form,...e},t.refused=null,t.told=!0,this.#u())}#a(){let e=this.#n;if(!e)return;if(yf(e.form,this.#e.strings())){e.told||(e.told=!0,this.#u());return}let t=this.#e.apply(e.request,e.form);if(t!==null){e.refused=t,this.#u();return}this.#s(),this.#e.restoreFocus(e.opener)}#o(){let e=this.#n;e&&(this.#s(),this.#e.restoreFocus(e.opener))}#s(){this.#c()?.close(),this.#n=null,this.#u()}cancelAll(){this.open&&this.#o()}#c(){return this.#t.querySelector(`dialog[data-testid="${Y}"]`)}#l(e){if(e.key===`Escape`){e.preventDefault(),this.#o();return}if(e.key===`Tab`){kr(e);return}let t=e.composedPath()[0];e.key===`Enter`&&t instanceof HTMLInputElement&&t.type!==`radio`&&(e.preventDefault(),this.#a())}#u(){T(this.#n?this.#d(this.#n):w,this.#t)}#d(e){let t=this.#e.strings(),n=e.form,r=yf(n,t),i=e.refused??(e.told?r:null),a=e=>!!r&&(e===`count`&&n.type===`MUL`&&n.count.trim()!==``||e===`name`&&n.type===`SUP`&&n.name.trim()!==``||e===`label`&&n.type===`SRU`&&n.label.trim()!==``),o=(e,t,n,r,i={})=>x`<label for="sg-${e}">${t}</label>
      <input
        id="sg-${e}"
        type="text"
        autocomplete="off"
        spellcheck="false"
        inputmode=${i.numeric?`numeric`:w}
        maxlength=${i.max??w}
        data-testid="${Y}.${e}"
        .value=${n}
        aria-invalid=${a(e)?`true`:`false`}
        aria-describedby="sg-error"
        @input=${e=>r(e.target.value)}
      />`,s=(e,t,n,r,i,a,o=!1)=>x`<label for="sg-${e}"
        >${t}</label
      >
      <select id="sg-${e}" data-testid="${Y}.${e}" ?disabled=${o} @change=${e=>a(e.target.value)}>
        ${r.map(e=>x`<option value=${e} ?selected=${e===n}>${i(e)}</option>`)}
      </select>`,c;switch(n.type){case`SUP`:c=o(`name`,t.sgroupName,n.name,e=>this.#i({name:e}),{max:15});break;case`DAT`:c=x`${s(`context`,t.sgroupContext,n.context,cf,e=>t.sgroupContextName(e),e=>this.#i({context:e}),e.editing)}
        ${o(`field-name`,t.sgroupFieldName,n.fieldName,e=>this.#i({fieldName:e}))}
        <label for="sg-field-value">${t.sgroupFieldValue}</label>
        <textarea
          id="sg-field-value"
          rows="2"
          data-testid="${Y}.field-value"
          .value=${n.value}
          aria-describedby="sg-error"
          @input=${e=>this.#i({value:e.target.value})}
        ></textarea>
        <span class="sg-placement-label" id="sg-placement-label">${t.sgroupPlacement}</span>
        <div class="radios" role="radiogroup" aria-labelledby="sg-placement-label">
          ${lf.map(e=>x`<label class="radio"
              ><input
                type="radio"
                name="sg-placement"
                value=${e}
                data-testid="${Y}.placement.${e}"
                .checked=${n.placement===e}
                @change=${()=>this.#i({placement:e})}
              />${t.sgroupPlacementName(e)}</label
            >`)}
        </div>`;break;case`MUL`:c=o(`count`,t.sgroupCount,n.count,e=>this.#i({count:e}),{numeric:!0});break;case`SRU`:c=x`${o(`label`,t.sgroupLabel,n.label,e=>this.#i({label:e}))}
        ${s(`connect`,t.sgroupConnect,n.connect,uf,e=>t.sgroupConnectName(e),e=>this.#i({connect:e}))}`}return x`<dialog
      class="dialog sgroup"
      data-testid=${Y}
      aria-labelledby="sg-title"
      @keydown=${e=>this.#l(e)}
      @cancel=${e=>{e.preventDefault(),this.#o()}}
    >
      <div class="dialog-head">
        <h2 id="sg-title">${t.sgroupDialog}</h2>
        <button type="button" class="close" data-testid="${Y}.close" aria-label=${t.close} @click=${()=>this.#o()}>×</button>
      </div>
      <div class="fields">
        ${s(`type`,t.sgroupType,n.type,sf,e=>t.sgroupTypeName(e),e=>this.#i({type:e}),e.editing)} ${c}
        <span class="error" id="sg-error" data-testid="${Y}.error" role=${e.refused?`alert`:w} aria-live="polite" ?hidden=${!i}>${i??``}</span>
      </div>
      <div class="dialog-foot">
        <button type="button" class="primary" data-testid="${Y}.ok" aria-disabled=${r?`true`:`false`} @click=${()=>this.#a()}>${t.ok}</button>
        <button type="button" data-testid="${Y}.cancel" @click=${()=>this.#o()}>${t.cancel}</button>
      </div>
    </dialog>`}},yv=`canvas.label-editor`,bv=`label-editor-message`,xv=class{#e;#t;#n=null;#r=null;constructor(e,t,n){this.#e=e,this.#t=t,this.#i=n}#i;set label(e){this.#i=e,this.#n?.setAttribute(`aria-label`,e)}get element(){return this.#n}get value(){return this.#n?.value??``}show(e,t,n){this.hide();let r=this.#e.ownerDocument.createElement(`input`);r.type=`text`,r.className=`label-editor`,r.setAttribute(`data-testid`,yv),r.setAttribute(`aria-label`,this.#i),r.setAttribute(`autocomplete`,`off`),r.setAttribute(`spellcheck`,`false`),r.maxLength=80,r.value=t,Object.assign(r.style,{height:`${1.5*n}px`,fontSize:`${n}px`}),r.addEventListener(`keydown`,this.#u),r.addEventListener(`blur`,this.#f),r.addEventListener(`mousedown`,this.#p),r.addEventListener(`input`,this.#d),this.#e.append(r),this.#n=r,this.#a=e,this.#s(),r.focus({preventScroll:!0}),r.select()}move(e,t){let n=this.#n;n&&(this.#a=e,Object.assign(n.style,{height:`${1.5*t}px`,fontSize:`${t}px`}),this.#s(),this.#r&&this.refuse(this.#r.textContent??``))}#a=null;#o=null;#s(){let e=this.#n,t=this.#a,n=this.#e.ownerDocument.defaultView;if(!e||!t||!n)return;let r=n.getComputedStyle(e),i=parseFloat(r.fontSize)||13;this.#o??=this.#e.ownerDocument.createElement(`canvas`).getContext(`2d`);let a=i*.6*Math.max(1,e.value.length);this.#o&&(this.#o.font=`${r.fontWeight} ${i}px ${r.fontFamily}`,a=Math.max(this.#o.measureText(e.value).width,this.#o.measureText(`M`).width));let o=parseFloat(r.paddingLeft)+parseFloat(r.paddingRight)+parseFloat(r.borderLeftWidth)+parseFloat(r.borderRightWidth),s=this.#c(),c=Math.min(s.width,Math.ceil(a+o+.5*i)),l=parseFloat(e.style.height)||1.5*i,u=(e,t,n)=>Math.min(Math.max(e,t),n),d=(e,t,n)=>u(u(e-t/2,0,Math.max(0,n-t)),e-t+4,e-4);Object.assign(e.style,{left:`${d(t.x,c,s.width)}px`,top:`${d(t.y,l,s.height)}px`,width:`${c}px`})}hide(){this.#l();let e=this.#n;e&&(this.#n=null,this.#a=null,e.removeEventListener(`keydown`,this.#u),e.removeEventListener(`blur`,this.#f),e.removeEventListener(`mousedown`,this.#p),e.removeEventListener(`input`,this.#d),e.remove())}refuse(e){let t=this.#n;if(!t)return;let n=this.#r;n||(n=this.#e.ownerDocument.createElement(`div`),n.className=`label-editor-message`,n.id=bv,n.setAttribute(`role`,`alert`),this.#e.append(n),this.#r=n),n.textContent=e;let r=parseFloat(t.style.left)||0,i=parseFloat(t.style.top)||0,a=i+(parseFloat(t.style.height)||0)+4,o=this.#c(),s={width:n.offsetWidth,height:n.offsetHeight};Object.assign(n.style,{left:`${Math.max(0,Math.min(r,o.width-s.width))}px`,top:`${a+s.height<=o.height?a:Math.max(0,i-s.height-4)}px`}),t.setAttribute(`aria-invalid`,`true`),t.setAttribute(`aria-describedby`,bv)}#c(){let e=this.#e.ownerDocument.defaultView,t=e?Sr(e.getComputedStyle(this.#e)):null;if(t)return t.content;let n=this.#e.getBoundingClientRect();return{width:n.width,height:n.height}}#l(){this.#r?.remove(),this.#r=null,this.#n?.removeAttribute(`aria-invalid`),this.#n?.removeAttribute(`aria-describedby`)}#u=e=>{e.key===`Enter`||e.key===`Tab`?(e.preventDefault(),this.#t.commit()):e.key===`Escape`&&(e.preventDefault(),this.#t.cancel())};#d=()=>{this.#r&&this.#l(),this.#s()};#f=e=>{this.#n&&this.#t.blur(e.relatedTarget)};#p=e=>{e.detail>=3&&this.#t.press(e.detail)&&e.preventDefault()}},Sv=`dialog.paste-special`,Cv=class{#e;#t;#n=null;constructor(e,t){this.#e=e,this.#t=t}get open(){return this.#n!==null}show(e,t){!this.open&&this.#t.isConnected&&(this.#n={as:e,text:``,error:null,opener:t},this.#l(),this.#a()?.showModal(),this.#o(`text`)?.focus())}cancel(){let e=this.#n;e&&(this.#i(),this.#e.restoreFocus(e.opener))}#r(){let e=this.#n;if(e&&e.text.trim()){if(!r(e.as,this.#e.queries())||!this.#e.canRead(e.text,e.as)){e.error=this.#e.strings().notReadAs(e.as),this.#l(),this.#o(`text`)?.focus();return}this.#i(),this.#e.paste(e.text,e.as),this.#e.restoreFocus(e.opener)}}#i(){this.#a()?.close(),this.#n=null,this.#l()}#a(){return this.#t.querySelector(`dialog[data-testid="${Sv}"]`)}#o(e){return this.#t.querySelector(`[data-testid="${Sv}.${e}"]`)}#s(e){let t=this.#n;t&&(Object.assign(t,e,{error:null}),this.#l())}#c(e){e.key===`Escape`?(e.preventDefault(),this.cancel()):e.key===`Tab`?kr(e):e.key===`Enter`&&(e.ctrlKey||e.metaKey)&&(e.preventDefault(),this.#r())}#l(){T(this.#n?this.#u(this.#n):w,this.#t)}#u(e){let t=!e.text.trim();return x`<dialog
      class="dialog paste-special"
      data-testid=${Sv}
      aria-labelledby="ps-title"
      aria-describedby="ps-hint"
      @keydown=${e=>this.#c(e)}
      @cancel=${e=>{e.preventDefault(),this.cancel()}}
    >
      <div class="dialog-head">
        <h2 id="ps-title">${this.#e.strings().pasteSpecial}</h2>
        <button type="button" class="close" data-testid="${Sv}.close" aria-label=${this.#e.strings().close} @click=${()=>this.cancel()}>×</button>
      </div>
      <p class="rg-hint" id="ps-hint">${this.#e.strings().pasteSpecialHint(this.#e.pasteKey())}</p>
      <div class="fields">
        <label for="ps-text">${this.#e.strings().pasteText}</label>
        <textarea
          id="ps-text"
          rows="8"
          spellcheck="false"
          data-testid="${Sv}.text"
          aria-invalid=${e.error?`true`:`false`}
          aria-describedby=${e.error?`ps-error`:w}
          .value=${e.text}
          @input=${e=>this.#s({text:e.target.value})}
        ></textarea>
        <label for="ps-format">${this.#e.strings().readAs}</label>
        <select id="ps-format" data-testid="${Sv}.format" @change=${e=>this.#s({as:e.target.value})}>
          ${fe.map(t=>x`<option value=${t} ?disabled=${!r(t,this.#e.queries())} ?selected=${t===e.as}>
              ${r(t,this.#e.queries())?this.#e.strings().formatName(t):this.#e.strings().format(`queryModeOnly`,{name:this.#e.strings().formatName(t)})}
            </option>`)}
        </select>
        <p class="error" id="ps-error" data-testid="${Sv}.error" role="alert" ?hidden=${!e.error}>${e.error??``}</p>
      </div>
      <div class="dialog-foot">
        <button type="button" class="primary" data-testid="${Sv}.ok" ?disabled=${t} @click=${()=>this.#r()}>${this.#e.strings().ok}</button>
        <button type="button" data-testid="${Sv}.cancel" @click=${()=>this.cancel()}>${this.#e.strings().cancel}</button>
      </div>
    </dialog>`}};function wv(e,t){let n=e.trim();if(!/^[+-]?(\d+([.,]\d*)?|[.,]\d+)$/.test(n))return null;let r=Number(n.replace(`,`,`.`));return!Number.isFinite(r)||r<t.min||r>t.max||t.integer&&!Number.isInteger(r)?null:r}var Tv=class{#e;#t;#n=null;constructor(e,t){this.#e=e,this.#t=t}get open(){return this.#n!==null}show(e,t){if(this.open||!this.#t.isConnected)return;this.#n={request:e,text:String(e.value),direction:e.direction??null,error:null,opener:t},this.#c(),this.#a()?.showModal();let n=this.#o();n?.focus(),n?.select()}cancel(){let e=this.#n;e&&(this.#i(),this.#e.restoreFocus(e.opener))}#r(){let e=this.#n;if(!e)return;let t=wv(e.text,e.request);if(t===null){e.error=e.request.invalid,this.#c();let t=this.#o();t?.focus(),t?.select();return}this.#i(),e.request.ok(t,e.direction),this.#e.restoreFocus(e.opener)}#i(){this.#a()?.close(),this.#n=null,this.#c()}#a(){return this.#t.querySelector(`dialog`)}#o(){return this.#t.querySelector(`input[type="text"]`)}#s(e){e.key===`Escape`?(e.preventDefault(),this.cancel()):e.key===`Tab`?kr(e):e.key===`Enter`&&e.target.tagName===`INPUT`&&(e.preventDefault(),this.#r())}#c(){T(this.#n?this.#l(this.#n):w,this.#t)}#l(e){let t=e.request,n=this.#e.strings(),r=t.id,i=`${r}-title`.replace(/\./g,`-`),a=`${r}-error`.replace(/\./g,`-`),o=`${r}-value`.replace(/\./g,`-`),s=(t,n)=>x`<label class="radio"
        ><input
          type="radio"
          name="${r}-direction"
          data-testid="${r}.${t}"
          .checked=${e.direction===t}
          @change=${()=>{e.direction=t,this.#c()}}
        />${n}</label
      >`;return x`<dialog
      class="dialog number-dialog"
      data-testid=${r}
      aria-labelledby=${i}
      @keydown=${e=>this.#s(e)}
      @cancel=${e=>{e.preventDefault(),this.cancel()}}
    >
      <div class="dialog-head">
        <h2 id=${i}>${t.title}</h2>
        <button type="button" class="close" data-testid="${r}.close" aria-label=${n.close} @click=${()=>this.cancel()}>×</button>
      </div>
      <div class="fields">
        <label for=${o}>${t.label}</label>
        <input
          id=${o}
          type="text"
          inputmode=${t.integer?`numeric`:`decimal`}
          spellcheck="false"
          autocomplete="off"
          data-testid="${r}.value"
          aria-invalid=${e.error?`true`:`false`}
          aria-describedby=${e.error?a:w}
          .value=${e.text}
          @input=${t=>{e.text=t.target.value,e.error&&(e.error=null,this.#c())}}
        />
        ${e.direction?x`<span id=${`${o}-direction`}>${n.rotateDirection}</span>
              <div class="radios" role="radiogroup" aria-labelledby=${`${o}-direction`}>${s(`cw`,n.clockwise)}${s(`ccw`,n.counterClockwise)}</div>`:w}
        <p class="error" id=${a} data-testid="${r}.error" role="alert" ?hidden=${!e.error}>${e.error??``}</p>
      </div>
      <div class="dialog-foot">
        <button type="button" class="primary" data-testid="${r}.ok" @click=${()=>this.#r()}>${n.ok}</button>
        <button type="button" data-testid="${r}.cancel" @click=${()=>this.cancel()}>${n.cancel}</button>
      </div>
    </dialog>`}},Ev=`canvas.selection-actions`,Dv=[{action:`flip.horizontal`,id:`flip-horizontal`,name:`flipHorizontal`,icon:`flipH`},{action:`flip.vertical`,id:`flip-vertical`,name:`flipVertical`,icon:`flipV`},{action:`delete`,id:`delete`,name:`deleteSelection`,icon:`delete`}],Ov=`http://www.w3.org/2000/svg`,kv=6,Av=4,jv=class{element;#e;#t=new Map;#n;constructor(e,t){this.#e=t;let n=e.createElement(`div`);n.className=`actions`,n.setAttribute(`role`,`toolbar`),n.setAttribute(`data-testid`,Ev),n.hidden=!0,Object.assign(n.style,{width:`96px`,height:`32px`}),this.element=n,this.#n=new I_({now:()=>performance.now(),boxes:e=>this.#i(e),after:(e,t)=>this.#e.after(e,t),changed:()=>this.#r()});for(let t of Dv){let r=`${Ev}.${t.id}`,i=`tip-${r.replace(/\./g,`-`)}`,a=e.createElement(`button`);a.type=`button`,a.className=`act`,a.tabIndex=-1,a.setAttribute(`data-testid`,r),a.setAttribute(`aria-describedby`,i),Object.assign(a.style,{width:`28px`,height:`28px`});let o=oc[t.icon],s=e.createElementNS(Ov,`svg`);s.setAttribute(`class`,`icon`),s.setAttribute(`viewBox`,o.viewBox),s.setAttribute(`aria-hidden`,`true`),s.setAttribute(`focusable`,`false`),s.innerHTML=o.body,a.append(s),a.addEventListener(`mousedown`,e=>e.preventDefault()),a.addEventListener(`click`,()=>{this.#n.pressed(),this.#e.act(t.action)}),a.addEventListener(`pointerenter`,e=>{q(e)||this.#n.enter(r)}),a.addEventListener(`pointerleave`,e=>{q(e)||this.#n.leave(r,e.pointerType===`touch`?null:{x:e.clientX,y:e.clientY})});let c=e.createElement(`span`);c.className=`tip`,c.id=i,c.setAttribute(`role`,`tooltip`),c.setAttribute(`popover`,`manual`),n.append(a,c),this.#t.set(r,{button:a,tip:c,action:t.action,name:t.name})}this.relabel()}relabel(){let e=this.#e.strings();this.element.setAttribute(`aria-label`,e.selectionActions);for(let{button:t,tip:n,action:r,name:i}of this.#t.values()){let a=e[i];t.setAttribute(`aria-label`,a);let o=this.#e.keys(r);if(o){let t=n.ownerDocument.createElement(`kbd`);t.textContent=o;let[r,i]=e.aroundHotkey(a);n.replaceChildren(r,t,i)}else n.textContent=a}}show(e){if(!e){this.element.hidden=!0,this.#n.gone(),this.#r();return}Object.assign(this.element.style,{left:`${e.x}px`,top:`${e.y}px`}),this.element.hidden=!1,this.#r()}pointerMoved(e,t){this.#n.moved({x:e,y:t})}pressed(){this.#n.pressed()}pointerGone(){this.#n.gone()}escape(){return this.#n.open!==null&&(this.#n.gone(),!0)}destroy(){this.#n.destroy();for(let{tip:e}of this.#t.values())wr(e)}#r(){let e=this.element.hidden?null:this.#n.open;for(let[t,{button:n,tip:r}]of this.#t){if(!(t===e&&r.isConnected)){r.removeAttribute(`data-testid`),wr(r);continue}r.setAttribute(`data-testid`,`${t}.tooltip`),En(r);let i=n.getBoundingClientRect(),a=r.getBoundingClientRect(),o=r.ownerDocument.defaultView?.innerWidth??1/0,s=i.top-kv-a.height;s<Av&&(s=i.bottom+kv);let c=Math.max(Av,Math.min(i.left+i.width/2-a.width/2,o-Av-a.width));hr(r,c,s)}}#i(e){let t=this.#t.get(e);if(!t||!Zn(t.tip))return null;let[n,r]=[t.button.getBoundingClientRect(),t.tip.getBoundingClientRect()];return{button:{left:n.left,top:n.top,right:n.right,bottom:n.bottom},tip:{left:r.left,top:r.top,right:r.right,bottom:r.bottom}}}},Mv=e=>e===`top`||e===`bottom`?`x`:`y`,Nv=2;function Pv(e,t,n,r,i){return e?!t||t.key!==n?{fold:!1,measured:t}:{fold:t.width-i>1,measured:t}:{fold:r-i>1,measured:{width:r,key:n}}}function Fv(e){let t=e.closest(`.frame`),n=Array.from(e.querySelectorAll(`button[data-testid^="toolbar."]`)).filter(e=>!e.classList.contains(`scroll`)&&!e.classList.contains(`folds`)&&!e.closest(`.palette:not(.fold-palette)`)).map(e=>e.getAttribute(`data-testid`));return`${t?.getAttribute(`data-density`)??``} ${t?.getAttribute(`lang`)??``} ${n.join(` `)}`}var Iv=(e,t,n,r)=>({testid:e,name:t,split:!0,items:n,segment:e=>Rt(e),itemName:(e,t)=>t.formatName(e),acts:(e,t)=>r(e,t)}),X={copy:Iv(`menu.copy-as`,`copyAs`,Et,e=>$e(e)),paste:Iv(`menu.paste-special`,`pasteSpecial`,fe,(e,t)=>r(e,t)),clean:{testid:`menu.clean-up`,name:`cleanUpMenu`,split:!0,items:[`clean`,`layout`],segment:e=>e,itemName:(e,t)=>D(e,t),acts:()=>!0,line:(e,t)=>e===`clean`?t.cleanUpLine:t.layoutLine},zoom:{testid:`menu.zoom`,name:`zoomMenu`,split:!1,items:[`zoom.in`,`zoom.out`,`zoom.reset`,`zoom.fit`],segment:e=>e.slice(5),itemName:(e,t)=>e===`zoom.reset`?t.zoomActualSize:D(e,t),acts:()=>!0},options:{testid:`menu.options`,name:`optionsMenu`,split:!1,items:[`settings`,`mode`,`keys`,`about`],segment:e=>e,itemName:(e,t)=>e===`settings`?t.settingsItem:e===`mode`?t.queryModeItem:e===`about`?t.aboutItem:t.keysItem,acts:()=>!0,checks:e=>e===`mode`},hydrogens:{testid:`menu.hydrogens`,name:`hydrogensMenu`,split:!1,items:[`hydrogens.add`,`hydrogens.remove`],segment:e=>e.slice(10),itemName:(e,t)=>e===`hydrogens.add`?t.addExplicitHydrogens:t.removeExplicitHydrogens,acts:(e,t)=>!(t&&e===`hydrogens.add`),why:(e,t)=>t.addHydrogensQueryMode,line:(e,t)=>e===`hydrogens.add`?t.addExplicitHydrogensLine:t.removeExplicitHydrogensLine,readsDrawing:!0}},Lv=e=>Object.prototype.hasOwnProperty.call(X,e),Rv=(e,t)=>`${X[e].testid}.${X[e].segment(t)}`,zv=e=>X[e].split?`toolbar.${e}.open`:`toolbar.${e}`,Bv={undo:`undo`,redo:`redo`,clear:`clear`,"element.table":`table`,erase:`erase`,"bond.single":`single`,"bond.double":`double`,"bond.triple":`triple`,"bond.aromatic":`aromatic`,"bond.wedge":`wedge`,"bond.hash":`hash`,"bond.wavy":`wavy`,"bond.crossed":`crossed`,"bond.dative":`dative`,"bond.any":`anyBond`,"bond.single-or-double":`singleOrDouble`,"bond.single-or-aromatic":`singleOrAromatic`,"bond.double-or-aromatic":`doubleOrAromatic`,chain:`chain`,"ring.benzene":`benzene`,"ring.cyclopentadiene":`cyclopentadiene`,"ring.cyclohexane":`cyclohexane`,"ring.cyclopentane":`cyclopentane`,"ring.cyclopropane":`cyclopropane`,"ring.cyclobutane":`cyclobutane`,"ring.cycloheptane":`cycloheptane`,"ring.cyclooctane":`cyclooctane`,"template.library":`library`,hand:`hand`,"zoom.in":`zoomIn`,"zoom.out":`zoomOut`,"zoom.reset":`zoomReset`,"zoom.fit":`fit`,layout:`layout`,clean:`clean`,cut:`cut`,copy:`copy`,paste:`paste`,"structure.enhanced-stereo":`enhancedStereo`,"structure.stereo-panel":`enhancedStereo`,"structure.aromatize":`aromatize`,"structure.dearomatize":`dearomatize`,options:`settings`,hydrogens:`hydrogens`,"rgroup.label":`rgroupLabel`,"rgroup.attachment":`attachmentPoint`,sgroup:`sgroup`,"select.rect":`selectRect`,"select.lasso":`selectLasso`,"select.structure":`selectStructure`,"charge.plus":`chargePlus`,"charge.minus":`chargeMinus`,keys:`keys`,settings:`settings`},Vv={layout:`layoutTip`,clean:`cleanUpTip`,"structure.stereo-panel":`stereoPanelTip`};function Hv(e,t,n){let[r,i]=n.aroundHotkey(e);return x`${r}<kbd>${t}</kbd>${i}`}function Uv(e,t){let n=mn(e);return n?x`<span class="symbol" aria-hidden="true" ${an({"--atom-colour":Ln(n,t)})}>${n}</span>`:Wv(Bv[e]??`open`)}var Wv=e=>x`<svg class="icon ${e===`open`?`opener-icon`:``}" viewBox=${oc[e].viewBox} aria-hidden="true" focusable="false">${s_(oc[e].body)}</svg>`;function Gv(e){e?.focus(e.closest(`.palette`)?{preventScroll:!0}:void 0)}var Kv=`[role="menuitem"], [role="menuitemcheckbox"]`,Z=(e,t)=>`${e}-${t.replace(/[^a-z0-9-]/g,`_`)}`,qv=class{#e;#t;#n=null;#r=null;#i=!1;#a=new Map;#o=new Map;#s;#c=null;#l=!1;#u=!1;#d=null;#f=null;#p=null;#m=null;constructor(e,t){this.#e=e,this.#t=t,this.#s=new I_({now:()=>performance.now(),boxes:e=>this.#K(e),after:(e,t)=>this.#e.after(e,t),changed:()=>this.#e.requestUpdate()})}set config(e){this.#t=e,this.#n=null,this.#f=null,this.#o.clear(),this.#e.requestUpdate()}get narrow(){return this.#u}set narrow(e){e!==this.#u&&(this.#u=e,this.#f=null,this.#n=null,this.#o.clear(),this.#e.requestUpdate())}fit(e){if(!e||e.clientWidth<=0)return;let t=e.querySelector(`.group.folded`)!==null,n=Pv(t,this.#d,Fv(e),e.scrollWidth,e.clientWidth);this.#d=n.measured,(n.fold!==t||!t)&&(this.narrow=n.fold)}toolChanged(){this.#o.clear();let e=this.#e.activeTool();for(let t of[this.#t.top,this.#t.left,this.#t.right,this.#t.bottom])for(let n of t)for(let t of n.items)t.kind===`palette`&&t.members.includes(e)&&this.#a.set(t.id,e);this.#e.requestUpdate()}render(e){let t=this.#t[e];if(!t.length)return w;let n=Mv(e);return x`<div class="bar bar-${e}" data-testid="toolbar.${e}" @scroll=${e=>this.#h(e.currentTarget,n)}>
      ${this.#g(e,-1)} ${t.map((t,n)=>this.#_(e,t,n))} ${this.#g(e,1)}
    </div>`}#h(e,t){Jv(e,t);for(let t of Array.from(e.querySelectorAll(`[popover], [data-popover]`)))Qv(t)&&ey(t)}#g(e,t){let n=this.#e.strings(),r=t<0?n.scrollBack:n.scrollForward,i=Mv(e);return x`<button
      type="button"
      class="scroll ${t<0?`back`:`forward`}"
      data-testid="toolbar.${e}.${t<0?`scroll-back`:`scroll-forward`}"
      aria-label=${r}
      title=${r}
      tabindex="-1"
      @mousedown=${e=>e.preventDefault()}
      @click=${e=>Xv(e.currentTarget.parentElement,t,i)}
    >
      ${Wv(i===`x`?t<0?`scrollBack`:`scrollForward`:t<0?`scrollUp`:`scrollDown`)}
    </button>`}#_(e,t,n){let r=`toolbar.${e}.${n}`,i=this.#S(t),a=this.#V(i);if(e===`top`&&this.#u&&a.length>Nv)return this.#v(r,t,i);let o=this.#H(r,a);return x`<div
      class="group"
      role="toolbar"
      data-testid=${r}
      aria-label=${ta(t.label,this.#e.strings())||w}
      aria-orientation=${e===`top`||e===`bottom`?`horizontal`:`vertical`}
    >
      ${en(i,e=>e.id,e=>this.#C(e,o))}
    </div>`}#v(e,t,n){let r=this.#e.strings(),i=n.find(e=>e.kind===`button`),a=n.some(e=>e.kind===`readout`&&e.id===`zoom.level`||e.kind===`button`&&e.id===`zoom`),o=a?this.#e.readout(`zoom.level`):``,s=ta(t.label,r)||(i?D(i.id,r):``),c=a?r.format(`zoomButton`,{percent:o}):s,l=`${e}.open`,u=this.#f===e,d=Z(`fold`,e);return x`<div class="group folded" role="toolbar" data-testid=${e} aria-label=${s} aria-orientation="horizontal">
      <div class="pal fold">
        ${this.#z({testid:l,label:c,tip:c,icon:a?x`<output class="level" data-testid="${e}.level">${o}</output>${Wv(`open`)}`:x`${i?Uv(i.id,this.#e.atomColours()):Wv(`open`)}<span class="corner" aria-hidden="true">${Wv(`corner`)}</span>`,tabbable:!0,popup:{expanded:u,controls:u?d:void 0},kind:a?`folds zoom`:`folds`,noTip:u,onClick:t=>this.#y(e,t.detail===0)})}
        <div class="palette fold-palette" popover="manual" role="group" id=${d} aria-label=${s} data-testid="${e}.palette" ?data-open=${u}>
          ${en(n,e=>e.id,e=>this.#C(e,``))}
        </div>
      </div>
    </div>`}#y(e,t){if(this.#f===e){this.#b();return}this.#n=null,this.#f=e,this.#e.requestUpdate(),t&&this.#x(`${e}.palette`)}#b(){this.#f!==null&&(this.#f=null,this.#n=null,this.#e.requestUpdate())}#x(e){let t=this.#e.root().querySelector(`[data-testid="${e}"]`);t&&($v(t),Gv(t.querySelector(`button[data-testid^="toolbar."]:not(.opener)`)))}#S(e){let t=e.items.findIndex(e=>e.kind===`button`&&e.id===`element.table`);if(t<0)return e.items;let n=this.#e.recents().map(e=>({kind:`button`,id:cr(e)})),r=this.#e.queryAtom()?[{kind:`button`,id:br}]:[];return[...e.items.slice(0,t),...n,...r,...e.items.slice(t)]}#C(e,t){return e.kind===`button`&&Lv(e.id)?X[e.id].split?this.#j(e.id,t):this.#N(e.id).length?this.#M(e.id,t):w:e.kind===`button`?this.#T(e.id,`toolbar.${e.id}`,t):e.kind===`readout`?this.#w(e.id):this.#E(e,t)}#w(e){return x`<output class="readout" data-testid="toolbar.${e}" aria-label=${D(e,this.#e.strings())}>${this.#e.readout(e)}</output>`}#T(e,t,n,r){let i=gr[e],a=this.#e.strings(),o=e===`element.query`?this.#e.queryAtom():null,s=o?.name??D(e,a),c=i.tool?e:void 0,l=c===void 0?this.#e.pressed?.(e):void 0,u=c!==void 0||i.disabled?String(c!==void 0&&this.#e.activeTool()===c):l===void 0?void 0:String(l),d=this.#e.hotkey(e),f=Vv[e],p=f?a[f]:s,m=c?this.#e.toolOff(c):null,h=i.disabled?a.notAvailableYet(s):m||(d?Hv(p,d,a):p),g=mn(e);return this.#z({testid:t,label:s,tip:h,icon:o?x`<span class="symbol query-atom ${o.face.length>2?`long`:``}" aria-hidden="true">${o.face}</span>`:Uv(e,this.#e.atomColours()),pressed:u,disabled:i.disabled||m!==null||!i.tool&&!this.#e.available(e),tabbable:n===t,haspopup:i.opens??(e===`element.table`||e===`keys`||e===`settings`?`dialog`:void 0),element:g??void 0,tile:g!==null,onClick:()=>this.#Y(e,r)})}#E(e,t){let n=this.#B(e),r=`toolbar.${e.id}`,i=this.#n===e.id,a=Z(`palette`,e.id),o=this.#e.strings(),s=D(n,o),c=na(e,o),l=this.#e.hotkey(n),u=this.#e.toolOff(n);return x`<div class="pal">
      ${this.#z({testid:r,label:s,tip:u??(l?Hv(s,l,o):s),disabled:u!==null,icon:x`${Uv(n,this.#e.atomColours())}<span class="corner" aria-hidden="true">${Wv(`corner`)}</span>`,tile:mn(n)!==null,pressed:String(this.#e.activeTool()===n),tabbable:t===r,popup:{expanded:i,controls:i?a:void 0},merged:!0,noTip:i,onClick:t=>this.#D(e,n,t),onPointerDown:t=>this.#O(e,t)})}
      <div class="palette" popover="manual" role="group" id=${a} aria-label=${c} data-testid="${r}.palette" ?data-open=${i}>
        ${e.members.map(t=>this.#T(t,`toolbar.${t}`,null,e))}
      </div>
    </div>`}#D(e,t,n){if(n.detail>0&&this.#m===e.id){this.#m=null;return}if(this.#m=null,n.detail>0&&this.#e.activeTool()===t){this.#$(e,!1);return}this.#ee(),this.#e.toolOff(t)||this.#e.chooseTool(t)}#O(e,t){if(this.#A(),this.#m=null,t.button!==0)return;let n=this.#e.after(500,()=>{this.#p=null,this.#n!==e.id&&(this.#m=e.id,this.#f=this.#f&&this.#e.root().querySelector(`[data-testid="${this.#f}.palette"] [data-testid="toolbar.${e.id}"]`)?this.#f:null,this.#n=e.id,this.#e.requestUpdate())});this.#p={palette:e.id,x:t.clientX,y:t.clientY,cancel:n}}#k(e){this.#p&&Math.hypot(e.clientX-this.#p.x,e.clientY-this.#p.y)>8&&this.#A()}#A(){this.#p?.cancel(),this.#p=null}#j(e,t){let n=X[e],r=this.#e.strings()[n.name],i=`toolbar.${e}.open`,a=this.#n===e,o=Z(`menu`,n.testid);return x`<div class="pal">
      ${this.#T(e,`toolbar.${e}`,t)}
      ${this.#z({testid:i,label:r,tip:r,icon:Wv(`open`),tabbable:t===i,disabled:!this.#e.menuAvailable(e),opener:{expanded:a,controls:a?o:void 0},noTip:a,onClick:t=>this.#X(e,t.detail===0)})}
      ${this.#P(e,a)}
    </div>`}#M(e,t){let n=this.#e.strings(),r=`toolbar.${e}`,i=this.#n===e,a=e===`zoom`,o=a?this.#e.readout(`zoom.level`):``,s=a?n.format(`zoomButton`,{percent:o}):n[X[e].name],c=!a&&e!==`options`,l=Uv(e,this.#e.atomColours());return x`<div class="pal">
      ${this.#z({testid:r,label:s,tip:s,icon:a?x`<output class="level" data-testid="toolbar.zoom.level">${o}</output>${Wv(`open`)}`:c?x`${l}<span class="chevron" aria-hidden="true">${Wv(`menuMark`)}</span>`:l,tabbable:t===r,disabled:!this.#e.menuAvailable(e),popup:{expanded:i,controls:i?Z(`menu`,X[e].testid):void 0},kind:a?`zoom`:c?`menu-button chevroned`:`menu-button`,noTip:i,onClick:t=>this.#X(e,t.detail===0)})}
      ${this.#P(e,i)}
    </div>`}#N(e){let t=X[e].items;return e===`copy`||e===`paste`?[...t]:t.filter(e=>!this.#e.hidden(e))}#P(e,t){let n=X[e],r=this.#e.strings(),i=this.#N(e),a=this.#e.queryMode(),o=i.filter(e=>!n.acts(e,a)),s=x`${i.map(t=>this.#R(e,t))}`,c=e===`zoom`?x`<div class="palette zoom-pop" popover="manual" ?data-open=${t}>
            <div class="menu" role="menu" id=${Z(`menu`,n.testid)} aria-label=${r[n.name]} data-testid=${n.testid}>${s}</div>
            ${this.#F()}
          </div>`:x`<div class="palette menu" popover="manual" role="menu" id=${Z(`menu`,n.testid)} aria-label=${r[n.name]} data-testid=${n.testid} ?data-open=${t}>
            ${s}
          </div>`;return x`${c}
      ${o.map(i=>{let a=Rv(e,i),o=t&&this.#U()===a;return x`<span class="tip" popover="manual" role="tooltip" id=${Z(`tip`,a)} data-anchor=${a} data-testid=${o?`${a}.tooltip`:w} ?data-open=${o}
          >${n.why?n.why(i,r):e===`paste`?r.format(`queryModeOnly`,{name:n.itemName(i,r)}):r.notAvailableYet(n.itemName(i,r))}</span
        >`})}`}#F(){let e=this.#e.strings(),t=this.#e.readout(`zoom.level`),n=e.format(`zoomRange`,{min:xt(Hr),max:xt(h)});return x`<div class="zoom-field">
      <label class="zoom-label" for="zoom-field-input">${e.zoomField}</label
      ><input
        id="zoom-field-input"
        class="zoom-input"
        type="text"
        inputmode="decimal"
        autocomplete="off"
        spellcheck="false"
        size="6"
        data-testid="menu.zoom.field"
        aria-invalid=${this.#i?`true`:w}
        aria-describedby="zoom-field-hint"
        .value=${this.#r??t}
        @input=${e=>{this.#r=e.target.value,this.#i=!1}}
        @keydown=${e=>this.#I(e)}
        @change=${e=>this.#L(e.target,!1)}
      /><span class="zoom-hint ${this.#i?`error`:``}" id="zoom-field-hint">${n}</span>
    </div>`}#I(e){let t=e.target;if(e.key===`Enter`)e.preventDefault(),this.#L(t,!0);else if(e.key===`Escape`)e.preventDefault(),this.#r=null,t.value=this.#e.readout(`zoom.level`),this.#ee(),this.#te(zv(`zoom`));else if(e.key===`ArrowDown`||e.key===`ArrowUp`){e.preventDefault();let n=Array.from(t.closest(`.zoom-pop`)?.querySelectorAll(Kv)??[]);Gv(e.key===`ArrowDown`?n[0]:n[n.length-1])}}#L(e,t){if(this.#n!==`zoom`)return;let n=this.#r,r=n===null?null:rt(e.value);if(n!==null&&r===null){this.#i=!0,this.#e.requestUpdate();return}this.#r=null,this.#i=!1,t&&(this.#ee(),this.#te(zv(`zoom`))),r!==null&&this.#e.zoomTo(r),this.#e.requestUpdate()}#R(e,t){let n=X[e],r=Rv(e,t),i=this.#e.strings(),a=n.itemName(t,i),o=n.acts(t,this.#e.queryMode()),s=o&&(n.readsDrawing&&this.#n!==e?!0:this.#e.menuItemAvailable(e,t)),c=o?this.#e.menuKeys(e,t):void 0,l=n.line?.(t,i),u=Z(`line`,r),d=n.checks!==void 0,f=n.checks?.(t)??!1,p=f&&this.#e.menuItemChecked(e,t);return x`<button
      type="button"
      role=${f?`menuitemcheckbox`:`menuitem`}
      class="item ${l?`lined`:``} ${d?`marked`:``}"
      data-testid=${r}
      aria-label=${a}
      aria-checked=${f?String(p):w}
      aria-disabled=${s?w:`true`}
      aria-describedby=${o?l?u:w:Z(`tip`,r)}
      aria-keyshortcuts=${c?.aria??w}
      tabindex="-1"
      @click=${()=>this.#Z(e,t)}
      @pointerenter=${e=>{q(e)||(o?this.#s.gone():this.#W(r,e))}}
      @pointerleave=${e=>this.#G(r,e)}
      @focus=${e=>this.#q(r,e.target)}
      @blur=${()=>this.#J(r)}
    >
      ${d?x`<span class="mark" aria-hidden="true">${p?`✓`:``}</span>`:w}<span class="item-name"
        >${a}${l?x`<span class="item-line" id=${u}>${l}</span>`:w}</span
      >${c?x`<kbd aria-hidden="true">${c.label}</kbd>`:w}
    </button>`}#z(e){let t=Z(`tip`,e.testid),n=!e.noTip&&this.#U()===e.testid,r=e.element,i=e.opener??e.popup,a=e.onPointerDown,o=r||a?e=>{r&&this.#e.pressElement(r,e),a?.(e)}:w;return x`<span
      class="slot"
      @pointerenter=${t=>this.#W(e.testid,t)}
      @pointerleave=${t=>this.#G(e.testid,t)}
    >
      <button
        type="button"
        class="btn ${e.opener?`opener`:``} ${e.merged?`merged`:``} ${e.kind??``} ${r?`element`:``} ${e.tile?`tile`:``}"
        data-testid=${e.testid}
        aria-label=${e.label}
        aria-describedby=${t}
        aria-pressed=${e.pressed??w}
        aria-disabled=${e.disabled?`true`:w}
        aria-haspopup=${i?`true`:e.haspopup??w}
        aria-expanded=${i?String(i.expanded):w}
        aria-controls=${i?.controls??w}
        tabindex=${e.tabbable?`0`:`-1`}
        @click=${e.onClick}
        @pointerdown=${o}
        @pointermove=${a?e=>this.#k(e):w}
        @pointerup=${a?()=>this.#A():w}
        @pointercancel=${a?()=>this.#A():w}
        @pointerleave=${a?()=>this.#A():w}
        @contextmenu=${a?e=>e.preventDefault():w}
        @focus=${t=>this.#q(e.testid,t.target)}
        @blur=${()=>this.#J(e.testid)}
      >
        ${e.icon}
      </button>
      <span class="tip" popover="manual" role="tooltip" id=${t} data-testid=${n?`${e.testid}.tooltip`:w} ?data-open=${n}>${e.tip}</span>
    </span>`}#B(e){let t=this.#a.get(e.id),n=e=>Xr(e)&&!this.#e.toolOff(e);return t&&e.members.includes(t)&&n(t)?t:e.members.find(n)??e.members.find(e=>Xr(e))??e.members[0]}#V(e){return e.flatMap(e=>e.kind===`button`?Lv(e.id)&&!X[e.id].split&&!this.#N(e.id).length?[]:[{testid:`toolbar.${e.id}`,tool:Xr(e.id)?e.id:void 0},...Lv(e.id)&&X[e.id].split?[{testid:`toolbar.${e.id}.open`}]:[]]:e.kind===`readout`?[]:[{testid:`toolbar.${e.id}`,tool:this.#B(e)}])}#H(e,t){let n=this.#o.get(e);if(n&&t.some(e=>e.testid===n))return n;let r=this.#e.activeTool();return t.find(e=>e.tool===r)?.testid??t[0]?.testid??``}#U(){return this.#l?null:this.#s.open??this.#c}#W(e,t){q(t)||(this.#l=!1,this.#s.enter(e))}#G(e,t){q(t)||this.#s.leave(e,t.pointerType===`touch`?null:{x:t.clientX,y:t.clientY})}#K(e){let t=this.#e.root(),n=t.querySelector(`[data-testid="${e}"]`),r=t.querySelector(`#${Z(`tip`,e)}`);if(!n||!r||!Qv(r))return null;let[i,a]=[n.getBoundingClientRect(),r.getBoundingClientRect()];return{button:{left:i.left,top:i.top,right:i.right,bottom:i.bottom},tip:{left:a.left,top:a.top,right:a.right,bottom:a.bottom}}}pointerMoved(e,t){this.#s.moved({x:e,y:t})}pressed(){this.#s.pressed()}pointerGone(){this.#s.gone()}destroy(){this.#s.destroy(),this.#A()}#q(e,t){let n=t.closest(`[role="toolbar"]`)?.getAttribute(`data-testid`);n&&!t.closest(`.palette`)&&this.#o.set(n,e),this.#c=t.matches(`:focus-visible`)?e:null,this.#l=!1,this.#e.requestUpdate()}#J(e){this.#c===e&&(this.#c=null,this.#e.requestUpdate())}#Y(e,t){let n=gr[e];if(this.#e.wasDragged()||n.disabled||n.tool&&this.#e.toolOff(e))return;let r=this.#f;if(r&&(n.tool||n.opens||this.#e.available(e))&&(this.#b(),this.#te(`${r}.open`)),n.opens){this.#ee(),this.#e.open(e);return}if(!n.tool){if(!this.#e.available(e))return;this.#ee(),this.#e.runAction(e);return}if(t){this.#a.set(t.id,e),this.#ee(),this.#e.chooseTool(e),this.#te(`toolbar.${t.id}`);return}this.#ee(),this.#e.chooseTool(e)}#X(e,t){if(this.#n===e){this.#ee();return}if(!this.#e.menuAvailable(e))return;let n=this.#N(e)[0];n!==void 0&&(this.#n=e,this.#e.requestUpdate(),t&&this.#te(Rv(e,n),!0))}#Z(e,t){if(!X[e].acts(t,this.#e.queryMode())||this.#n!==e||!this.#e.menuItemAvailable(e,t))return;this.#ee();let n=this.#f;this.#b(),this.#te(n?`${n}.open`:zv(e)),this.#e.runMenu(e,t)}#Q(e,t){let n=t.closest(`.zoom-pop`);if(n&&/^[0-9]$/.test(e.key)){this.#r=e.key,this.#i=!1;let t=n.querySelector(`.zoom-input`);return t&&(t.value=e.key,Gv(t)),this.#e.requestUpdate(),!0}if(![`ArrowUp`,`ArrowDown`,`Home`,`End`].includes(e.key))return!1;let r=Array.from((n??t.closest(`[role="menu"]`))?.querySelectorAll(n?`${Kv}, .zoom-input`:Kv)??[]),i=r.indexOf(t),a=r.length;return a?(Gv(r[e.key===`Home`?0:e.key===`End`?a-1:(i+(e.key===`ArrowDown`?1:-1)+a)%a]),!0):!1}#$(e,t){if(this.#n===e.id){this.#ee();return}this.#n=e.id,this.#e.requestUpdate(),t&&this.#te(`toolbar.${e.members[0]}`,!0)}#ee(){this.#n!==null&&(this.#n=null,this.#r=null,this.#i=!1,this.#e.requestUpdate())}#te(e,t=!1){let n=()=>Gv(this.#e.root().querySelector(`[data-testid="${e}"]`));if(!t){n();return}let r=this.#e.root().querySelector(`[data-testid="${e}"]`)?.closest(`.palette`);r&&$v(r),n()}sync(){for(let e of Array.from(this.#e.root().querySelectorAll(`[popover], [data-popover]`)))e.closest(`[data-menus]`)||(e.hasAttribute(`data-open`)?$v(e):wr(e))}keyDown(e,t){if(e.key===`Escape`)return this.escape();let n=t.closest(`button.merged[data-testid^="toolbar."]`);if(n&&e.key===`ArrowDown`&&!e.ctrlKey&&!e.metaKey&&!e.shiftKey){let t=n.closest(`[role="toolbar"]`)?.getAttribute(`aria-orientation`)===`vertical`,r=this.#ie(n.getAttribute(`data-testid`).slice(8));if(r&&(e.altKey||!t))return this.#n===r.id?this.#te(`toolbar.${r.members[0]}`,!0):this.#$(r,!0),!0}if(e.altKey||e.ctrlKey||e.metaKey)return!1;let r=t.closest(Kv);if(r)return this.#Q(e,r);let i=t.closest(`button[data-testid^="toolbar."]`);if(!i)return!1;let a=/^(toolbar\.top\.\d+)\.open$/.exec(i.getAttribute(`data-testid`)??``);if(a&&e.key===`ArrowDown`)return this.#f===a[1]?this.#x(`${a[1]}.palette`):this.#y(a[1],!0),!0;let o=this.#ne(i.getAttribute(`data-testid`)??``);if(o&&e.key===`ArrowDown`){let e=this.#N(o)[0];return this.#n===o?e!==void 0&&this.#te(Rv(o,e),!0):this.#X(o,!0),!0}let s=i.closest(`.palette`),c=s??i.closest(`[role="toolbar"]`);if(!c)return!1;let[l,u]=s!==null||c.getAttribute(`aria-orientation`)===`vertical`?[`ArrowUp`,`ArrowDown`]:[`ArrowLeft`,`ArrowRight`];if(![l,u,`Home`,`End`].includes(e.key))return!1;let d=Array.from(c.querySelectorAll(`button[data-testid^="toolbar."]`)).filter(e=>(e.closest(`.palette`)??null)===s),f=d.indexOf(i),p=d.length;return Gv(d[e.key===`Home`?0:e.key===`End`?p-1:(f+(e.key===u?1:-1)+p)%p]),!0}#ne(e){for(let t of Object.keys(X))if(zv(t)===e)return t;return null}#re(e){return Lv(e)?zv(e):`toolbar.${e}`}#ie(e){for(let t of[this.#t.top,this.#t.left,this.#t.right,this.#t.bottom])for(let n of t)for(let t of n.items)if(t.kind===`palette`&&t.id===e)return t;return null}escape(){if(this.#n!==null){let e=this.#re(this.#n);return this.#ee(),this.#l=!0,this.#te(e),!0}if(this.#f!==null){let e=`${this.#f}.open`;return this.#b(),this.#l=!0,this.#te(e),!0}return this.#U()!==null&&(this.#l=!0,this.#e.requestUpdate(),!0)}focusOut(e){let t=e.relatedTarget;if(this.#f!==null){let e=this.#e.root().querySelector(`[data-testid="${this.#f}.open"]`)?.closest(`.pal`);if(!(t&&e?.contains(t))){this.#b();return}}if(this.#n===null)return;let n=this.#e.root().querySelector(`[data-testid="${this.#re(this.#n)}"]`)?.closest(`.pal`);t&&n?.contains(t)||this.#ee()}};function Jv(e,t=`x`){if(!e)return;let[n,r,i]=t===`x`?[e.scrollLeft,e.clientWidth,e.scrollWidth]:[e.scrollTop,e.clientHeight,e.scrollHeight],a=n>1,o=n+r<i-1;a!==e.hasAttribute(`data-back`)&&e.toggleAttribute(`data-back`,a),o!==e.hasAttribute(`data-forward`)&&e.toggleAttribute(`data-forward`,o)}var Yv=e=>Math.max(48,Math.round(e*.7));function Xv(e,t,n=`x`){if(!e)return;let r=e.ownerDocument.defaultView?.matchMedia?.(`(prefers-reduced-motion: reduce)`).matches??!1,i=t*Yv(n===`x`?e.clientWidth:e.clientHeight);e.scrollBy({...n===`x`?{left:i}:{top:i},behavior:r?`instant`:`smooth`}),Jv(e,n)}var Zv=6,Q=4,Qv=Zn;function $v(e){e.isConnected&&(En(e),ey(e))}function ey(e){let t=e.classList.contains(`palette`),n=!t&&e.dataset.anchor?e.getRootNode().querySelector?.(`[data-testid="${e.dataset.anchor}"]`):null,r=t?e.closest(`.pal`):n??e.previousElementSibling;if(!r)return;let i=r.getBoundingClientRect(),a=e.getBoundingClientRect(),o=e.ownerDocument.defaultView,[s,c]=[o?.innerWidth??1/0,o?.innerHeight??1/0],l=t?!e.closest(`.bar-top`)&&!e.closest(`.bar-bottom`):!!(n??e.closest(`.palette`)??e.closest(`.bar-left`)),u=!t&&!e.closest(`.palette`)&&!!e.closest(`.bar-right`),d=!!e.closest(`.bar-bottom`)&&(t||!e.closest(`.palette`)),f,p;if(d){f=t?i.left:i.left+i.width/2-a.width/2,p=i.top-Zv-a.height,p<Q&&(p=i.bottom+Zv);let n=e.getRootNode().querySelector?.(`.bar-left`)?.getBoundingClientRect();n&&n.width&&f<n.right&&p<n.bottom&&p+a.height>n.top&&(f=n.right+Zv),f+a.width>s-Q&&(f=s-Q-a.width)}else if(u)f=i.left-Zv-a.width,p=i.top+i.height/2-a.height/2,f<Q&&(f=i.right+Zv);else if(l)f=i.right+Zv,p=t?i.top:i.top+i.height/2-a.height/2,f+a.width>s-Q&&(f=i.left-Zv-a.width),p+a.height>c-Q&&(p=t?i.bottom-a.height:c-Q-a.height);else{f=i.left,p=i.bottom+Zv,p+a.height>c-Q&&(p=i.top-Zv-a.height);let t=e.getRootNode().querySelector?.(`.bar-left`)?.getBoundingClientRect();t&&t.width&&f<t.right&&p<t.bottom&&p+a.height>t.top&&(f=t.right+Zv),f+a.width>s-Q&&(f=s-Q-a.width)}hr(e,Math.max(Q,Math.min(f,s-Q-a.width)),Math.max(Q,Math.min(p,c-Q-a.height)))}var ty=new WeakMap;function ny(e,t){ty.set(e,t)}var ry=new WeakMap;function iy(e,t){ry.set(e,t)}function ay(e,t,n){ry.get(e)?.(t,n)}function oy(e){let t=ty.get(e);if(!t)throw Error(`Crux Sketch: not a <crux-sketch> made by this copy of the sketcher`);return t}var sy=Object.freeze({canvas:null,toolbars:Object.freeze({}),controls:Object.freeze({}),labelEditor:null}),cy=Object.freeze({atoms:Object.freeze([]),bonds:Object.freeze([])}),ly=e=>e!==null&&e.width>0&&e.height>0,uy=(e,t)=>{let[n,r]=[e.left+e.width/2,e.top+e.height/2];return n>=t.left&&n<=t.left+t.width&&r>=t.top&&r<=t.top+t.height};function dy(e){if(!e)return sy;let{origin:t}=e,n=e=>Object.freeze({x:e.left-t.left,y:e.top-t.top,width:e.width,height:e.height}),r={};for(let[t,i]of e.toolbars)ly(i)&&(r[t]=n(i));let i={};for(let t of e.controls)!ly(t.box)||!ly(t.holder)||!uy(t.box,t.holder)||t.name in i||(i[t.name]=n(t.box));return Object.freeze({canvas:ly(e.canvas)?n(e.canvas):null,toolbars:Object.freeze(r),controls:Object.freeze(i),labelEditor:ly(e.labelEditor)?n(e.labelEditor):null})}function fy(e){if(!e||!e.atoms.length&&!e.bonds.length)return cy;let t=e=>Object.freeze([...new Set(e)].sort((e,t)=>e-t));return Object.freeze({atoms:t(e.atoms),bonds:t(e.bonds)})}var py=`crux-sketch`,my=3,hy=2,gy=1,_y=0,vy=`canvas-description`,yy=[`top`,`left`,`right`,`bottom`],by=200,xy=e=>e.length>48?`${e.slice(0,47)}…`:e,Sy=`.bar button[data-testid^="toolbar."]:not(.scroll), .bar output[data-testid^="toolbar."]`,Cy=e=>e.preventDefault(),wy=[`pointermove`,`pointerup`,`pointercancel`,`pointerleave`,`mousemove`,`mouseup`,`mouseleave`];function Ty(e,t){let n=[],r=[];for(let r of Array.from(e?.querySelectorAll(`g[data-atom]`)??[]))n[Number(r.getAttribute(`data-atom`))]={id:Number(r.getAttribute(`data-id`)),label:Ps(t,r.getAttribute(`data-label`)??r.getAttribute(`data-element`)??`?`,Number(r.getAttribute(`data-charge`)??0)),x:Number(r.getAttribute(`data-x`)),y:Number(r.getAttribute(`data-y`))};for(let t of Array.from(e?.querySelectorAll(`g[data-bond]`)??[])){let[e,n]=(t.getAttribute(`data-atoms`)??``).split(` `).map(Number);r[Number(t.getAttribute(`data-bond`))]={id:Number(t.getAttribute(`data-id`)),a:e,b:n,type:Ms(t.getAttribute(`data-order`)??``,t.getAttribute(`data-stereo`)??``)}}return{atoms:n.filter(Boolean),bonds:r.filter(Boolean)}}var Ey=16777216,Dy=new Set([`clear`,`undo`,`redo`,`layout`,`clean`,`cut`,`paste`,`element.table`,`structure.aromatize`,`structure.dearomatize`]),Oy={"bond.single":`tool.bond.single`,"ring.benzene":`tool.ring.benzene`,"ring.cyclopentadiene":`tool.ring.cyclopentadiene`,chain:`tool.chain`,clear:`clear`,undo:`undo`,redo:`redo`,hand:`tool.hand`,"zoom.in":`zoom.in`,"zoom.out":`zoom.out`,"zoom.reset":`zoom.reset`,"zoom.fit":`zoom.fit`,layout:`layout`,clean:`clean`,"rgroup.label":`atom.rgroup`,cut:`cut`,copy:`copy`,paste:`paste`,"charge.plus":`atom.charge.plus`,"charge.minus":`atom.charge.minus`,keys:`keys.open`,settings:`settings.open`,"structure.aromatize":`structure.aromatize`,"structure.dearomatize":`structure.dearomatize`},ky={clean:`clean`,layout:`layout`,"zoom.in":`zoom.in`,"zoom.out":`zoom.out`,"zoom.reset":`zoom.reset`,"zoom.fit":`zoom.fit`,settings:`settings.open`,keys:`keys.open`},Ay={"zoom.in":`in`,"zoom.out":`out`,"zoom.reset":`reset`,"zoom.fit":`fit`},jy=e=>tt(e),My=(e,t)=>tt(`inset 0 0 0 ${e}px ${t}`),Ny=tt(Object.entries(Hi).map(([e,t])=>`.pt-el.${e} { background: ${t.fill}; box-shadow: inset 0 -3px 0 ${t.stripe}; }`).join(`
    `)),Py=class e extends Mn{static styles=Ae`
    :host {
      display: block;
      position: relative;
      /* Its own stacking context (spike mac-fonts, HOST-031): the parts it stacks (its bars, its readout, its panels)
         stack among themselves, never among the page's own boxes, so a box the page lays over it is above all of it (a
         host's dialog at z-index 1 beside another sketcher, api/host-drags on a Mac, 2026-10-21). What opens in the top
         layer is above every stacking context anyway; only the fixed fallback for a browser without the Popover API
         (menus-top-layer) now stays inside it. */
      isolation: isolate;
      box-sizing: border-box;
      width: 100%;
      height: 100%;
      overflow: hidden;
    }
    :host([hidden]) {
      display: none;
    }
    /* No box round the whole sketcher on focus (spike reports-oct14, report 1). */
    :host(:focus),
    :host(:focus-visible) {
      outline: none;
    }
    /* Nothing the page's styles pass down reaches the sketcher (spike host-isolation, I4; HOST-031): the boxes at the top
       of its shadow root, the frame and the dialogs' layers, start from the initial values, as if they had no parent,
       whatever font, colour, line height, spacing, direction or colour scheme the page sets. The page decides only whether
       the sketcher shows and takes the pointer. The all property leaves the direction alone, so it is set here; the rules after these
       give the frame and the dialogs their own looks. */
    .frame,
    .layer {
      all: initial;
      direction: ltr;
      visibility: inherit;
      pointer-events: inherit;
      font: 12px/1.3 system-ui, sans-serif;
      color: #1f2328;
    }
    .layer {
      display: block;
    }
    /* A host's pointer crossing the sketcher (spike polish-look, row 6): a press made outside it, a dialog's edge or a
       splitter dragged across it. Meanwhile the frame takes no pointer, the hit shapes' own pointer-events included, so that
       no button shows its hover look and the canvas not its cursor; the pointer keeps the page's. Chrome keeps the hovered
       state of what the host's pointer first landed on until the pointer moves again, so the hover looks are given only
       outside such a frame (their rules' :where(:not(.frame[data-host-pointer] *)), which weighs nothing). */
    .frame[data-host-pointer],
    .frame[data-host-pointer] * {
      pointer-events: none !important;
    }
    .frame {
      display: grid;
      grid-template-columns: auto minmax(0, 1fr) auto;
      grid-template-rows: auto minmax(0, 1fr) auto;
      width: 100%;
      height: 100%;
      font: 12px/1.3 system-ui, sans-serif;
      color: #1f2328;
      /* The bars' sizes (spike bar-fit, B1): regular here, compact below (the default since spike compact-default, in a box
         of any size). A target's short side never goes under 24 px (A11Y-015): the openers and the scroll buttons keep it
         in both. */
      --bar-pad: 4px;
      --bar-gap: 8px;
      --btn: 32px;
      --icon: 24px;
      --symbol: 14px;
      --readout: 44px;
      --target: 24px;
      --scroll-icon: 18px;
      --group-gap: 2px;
    }
    /* Compact (B1): the default density in a box of any size (spike compact-default), and the 'auto' density's in a box under
       720 by 600 px, such as Datagrok's sketcher dialog. */
    .frame[data-density='compact'] {
      --bar-pad: 3px;
      --bar-gap: 6px;
      --btn: 28px;
      --icon: 20px;
      --symbol: 13px;
      --readout: 40px;
      --scroll-icon: 16px;
      --group-gap: 1px;
    }
    /* No bar shows a scrollbar, in any browser (spike reports-oct5, R1: the product owner's ruling, over bar-fit's B2): a bar
       that overflows scrolls by its scroll buttons, the wheel, a finger and a focused control scrolling into view. With no
       scrollbar there is no gutter to keep either, so nothing moves as a bar scrolls (B4). */
    .bar {
      position: relative;
      z-index: 2;
      display: flex;
      box-sizing: border-box;
      min-width: 0;
      min-height: 0;
      gap: var(--bar-gap);
      padding: var(--bar-pad);
      /* White, as the canvas and a host's dialog are (the product owner, 2026-10-12, white-bars): the line on
         each bar's canvas side, below, keeps the drawing area's edge. */
      background: #ffffff;
      scrollbar-width: none;
    }
    /* The same in a WebKit browser that does not read scrollbar-width yet (Safari before 18.2). */
    .bar::-webkit-scrollbar {
      display: none;
    }
    /* The top bar scrolls sideways in a box narrower than its buttons, as the ring bar does (spike view:
       with the zoom group it is about 300 px wide). */
    .bar-top {
      grid-row: 1;
      grid-column: 1 / -1;
      flex-direction: row;
      border-bottom: 1px solid #c3c8d0;
      overflow-x: auto;
      overflow-y: hidden;
      /* A focused control scrolls into view clear of the scroll buttons (K8); a scroll moves nothing outside it (B4). */
      scroll-padding-inline: calc(var(--target) + var(--bar-gap));
      /* It moves smoothly, by its buttons or by focus moving along it, and at once under reduced motion (A11Y-019; the
         rule below; milestone-4-fixes, F9). */
      scroll-behavior: smooth;
    }
    /* A bar's scroll buttons (K8; the side bars' too since spike bar-fit, B3): at the ends with more, over the bar's ends,
       out of the flow (their size and the bar's gap given back, so that one coming or going moves nothing, B4), shown only
       while more lies that way. A sideways bar's are 24 px wide and a button tall, a side bar's a button wide and 24 px tall. */
    .bar .scroll {
      position: sticky;
      z-index: 1;
      flex: none;
      display: none;
      align-items: center;
      justify-content: center;
      box-sizing: border-box;
      width: var(--target);
      height: var(--btn);
      margin: 0;
      padding: 0;
      border: 0;
      border-radius: 4px;
      background: #ffffff;
      color: #1f2328;
      cursor: pointer;
      box-shadow: 0 0 6px 4px #ffffff;
    }
    .bar .scroll .icon {
      width: var(--scroll-icon);
      height: var(--scroll-icon);
    }
    .bar .scroll.back {
      left: 0;
      margin-right: calc(-1 * (var(--target) + var(--bar-gap)));
    }
    .bar .scroll.forward {
      right: 0;
      margin-left: calc(-1 * (var(--target) + var(--bar-gap)));
    }
    .bar-left .scroll,
    .bar-right .scroll {
      width: var(--btn);
      height: var(--target);
    }
    .bar-left .scroll.back,
    .bar-right .scroll.back {
      top: 0;
      left: auto;
      margin: 0 0 calc(-1 * (var(--target) + var(--bar-gap)));
    }
    .bar-left .scroll.forward,
    .bar-right .scroll.forward {
      bottom: 0;
      right: auto;
      margin: calc(-1 * (var(--target) + var(--bar-gap))) 0 0;
    }
    .bar[data-back] .scroll.back,
    .bar[data-forward] .scroll.forward {
      display: flex;
    }
    @media (hover: hover) {
      .bar .scroll:where(:not(.frame[data-host-pointer] *)):hover {
        background: #dde3ec;
      }
    }
    /* A folded group's button (K8) looks like a palette's (spike reports-oct5, R2; B10): a button's size, its first control's
       face and the same corner mark (the rules for .btn.merged below), no chevron under it. */
    .fold-palette {
      flex-direction: column;
    }
    .fold-palette .pal,
    .fold-palette .group {
      flex-direction: row;
    }
    .bar-top .slot,
    .bar-top .readout {
      flex: none;
    }
    .bar-left {
      grid-row: 2;
      grid-column: 1;
      flex-direction: column;
      border-right: 1px solid #c3c8d0;
    }
    /* The side bars scroll in their column when the canvas is shorter than their buttons (reading 30; milestone-3-fixes,
       F1), with no scrollbar and no room kept for one (R1); their scroll buttons say that more lies that way and take a
       press (B3), where a fade and a chevron took none; a focused control scrolls into view clear of them. */
    .bar-left,
    .bar-right {
      overflow-x: hidden;
      overflow-y: auto;
      scroll-padding-block: calc(var(--target) + var(--bar-gap));
      scroll-behavior: smooth;
    }
    .bar-left .slot {
      flex: none;
    }
    /* The ring bar (spike rings-chains, UI-001): a row below the canvas, spanning the box as Ketcher's
       bottom row does; it scrolls sideways in a box narrower than its buttons. */
    .bar-bottom {
      grid-row: 3;
      grid-column: 1 / -1;
      flex-direction: row;
      border-top: 1px solid #c3c8d0;
      overflow-x: auto;
      overflow-y: hidden;
      /* As the top bar's: its scroll buttons at the ends with more, no scrollbar (R1). */
      scroll-padding-inline: calc(var(--target) + var(--bar-gap));
      /* Smoothly, and at once under reduced motion, as the top bar (F9). */
      scroll-behavior: smooth;
    }
    .bar-bottom .slot {
      flex: none;
    }
    /* The element palette (spike atoms): a column on the right of the canvas, which scrolls in its
       column when the canvas is shorter than its buttons (reading 30). */
    .bar-right {
      grid-row: 2;
      grid-column: 3;
      flex-direction: column;
      border-left: 1px solid #c3c8d0;
    }
    .group,
    .pal {
      position: relative;
      display: flex;
      gap: var(--group-gap);
    }
    .bar-top .group,
    .bar-top .pal,
    .bar-bottom .group,
    .bar-bottom .pal {
      flex-direction: row;
    }
    .bar-left .group,
    .bar-left .pal,
    .bar-right .group,
    .bar-right .pal {
      flex-direction: column;
    }
    .bar-right .slot {
      flex: none;
    }
    /* A finger's vertical swipe scrolls the element palette; a drag out of the bar is a palette drag
       (fix round 4, the ui-critic's blocking finding). */
    .btn.element {
      touch-action: pan-y;
    }
    /* An element symbol is drawn in the colour the canvas draws that element in, which the toolbars hand it as
       --atom-colour (ATOM-033). */
    .symbol {
      font: 600 var(--symbol)/1 system-ui, sans-serif;
      color: var(--atom-colour);
      pointer-events: none;
    }
    /* The query atom tool's face (spike query-edit, Q8): a generic's symbol, black; a list's elements, smaller, cut at the edge. */
    .symbol.query-atom {
      max-width: 100%;
      overflow: hidden;
      color: #1f2328;
      white-space: nowrap;
    }
    .symbol.query-atom.long {
      font-size: 10px;
      letter-spacing: -0.02em;
    }
    .slot {
      position: relative;
      display: flex;
    }
    .btn {
      box-sizing: border-box;
      display: flex;
      align-items: center;
      justify-content: center;
      width: var(--btn);
      height: var(--btn);
      margin: 0;
      padding: 0;
      border: 0;
      border-radius: 4px;
      background: transparent;
      color: #1f2328;
      cursor: pointer;
    }
    .bar-left .opener {
      height: var(--target);
    }
    /* The zoom level (spike view, CANVAS-005): a readout between zoom out and zoom in, wide enough for
       "400%" so that the buttons after it never move. */
    .readout {
      box-sizing: border-box;
      display: flex;
      align-items: center;
      justify-content: center;
      min-width: var(--readout);
      height: var(--btn);
      padding: 0 2px;
      font: 12px/1 system-ui, sans-serif;
      font-variant-numeric: tabular-nums;
      color: #1f2328;
      white-space: nowrap;
    }
    .bar-left .readout,
    .bar-right .readout {
      min-width: var(--btn);
      height: var(--target);
      font-size: 11px;
    }
    .bar-top .opener,
    .bar-bottom .opener {
      width: var(--target);
    }
    .icon {
      width: var(--icon);
      height: var(--icon);
      pointer-events: none;
    }
    .opener-icon {
      width: 14px;
      height: 14px;
      transform: rotate(-90deg);
    }
    /* The top bar's openers (Copy As and Paste Special, spike copy-as) open their menus below them: they point down. */
    .bar-top .opener-icon {
      transform: none;
    }
    /* A palette's button (spike bar-fit, B10, Ketcher's multi-tool button): its member's face and a small triangle in its
       corner, which says a palette opens from it; a cue, which takes no press of its own. */
    .btn.merged,
    .btn.folds {
      position: relative;
    }
    .btn.merged .corner,
    .btn.folds .corner {
      position: absolute;
      right: 2px;
      bottom: 2px;
      display: flex;
      pointer-events: none;
    }
    .btn.merged .corner .icon,
    .btn.folds .corner .icon {
      width: 6px;
      height: 6px;
    }
    /* A menu button's chevron (spike ui-polish-oct12, row 1: the Hydrogens button): inside its square, in the strip under
       its face at the lower right, pointing down as Copy's and Paste's openers do; a cue, which takes no press of its own,
       so the button keeps bar-fit's square and the top bar its width (milestone-6-fixes' Q1). */
    .btn.chevroned {
      position: relative;
    }
    .btn.chevroned .chevron {
      position: absolute;
      right: 1px;
      bottom: 0;
      display: flex;
      pointer-events: none;
    }
    .btn.chevroned .chevron .icon {
      width: 7px;
      height: calc((var(--btn) - var(--icon)) / 2);
      max-height: 4px;
    }
    /* The zoom (B7, Ketcher's): one button showing the level, a chevron after it, as wide as "400%" needs. */
    .btn.zoom {
      width: auto;
      min-width: calc(var(--readout) + 14px);
      gap: 1px;
      padding: 0 2px 0 4px;
      font: 12px/1 system-ui, sans-serif;
      font-variant-numeric: tabular-nums;
      white-space: nowrap;
    }
    .btn.zoom .opener-icon {
      width: 12px;
      height: 12px;
      transform: none;
    }
    .bar-left .btn.zoom,
    .bar-right .btn.zoom {
      width: var(--btn);
      min-width: 0;
      flex-direction: column;
      padding: 0;
      font-size: 10px;
    }
    /* Hover only where a pointer can hover (a tap would leave it stuck), and never on a disabled button,
       which would then look like one that acts (the ui-critic's undo review). */
    @media (hover: hover) {
      .btn:not([aria-disabled='true']):where(:not(.frame[data-host-pointer] *)):hover {
        background: #dde3ec;
      }
      .btn[aria-pressed='true']:not([aria-disabled='true']):where(:not(.frame[data-host-pointer] *)):hover {
        background: #b9d2f7;
      }
    }
    .btn[aria-pressed='true'] {
      background: #cfe0fb;
      box-shadow: inset 0 0 0 1px #1c65c0;
    }
    /* The tiles of the buttons that show an element's symbol (spike palette-colours, ATOM-033, P3; the numbers
       are in config/panel-fills.ts). The symbols are drawn in the canvas's hues, and O, F, P, S and Cl keep
       4.5:1 only on a fill at least about 0.935 light (the bars are white since white-bars, but their hover and
       pressed fills are 0.76 and darker). So such a button has an opaque fill of its own in every state, white
       or a pale tint, and what the other bars show by a darker fill it shows by a ring: hovered, a thin
       grey-blue one; pressed, a 2 px blue one (never by colour alone, WCAG 1.4.1, 1.4.11). */
    .btn.tile {
      background: ${jy(Vi.rest)};
    }
    .btn.tile[aria-pressed='true'] {
      background: ${jy(Vi.lit)};
      box-shadow: ${My(Vi.pressedRingPx,Vi.pressedRing)};
    }
    @media (hover: hover) {
      .btn.tile:not([aria-disabled='true']):where(:not(.frame[data-host-pointer] *)):hover {
        background: ${jy(Vi.lit)};
        box-shadow: ${My(Vi.hoverRingPx,Vi.hoverRing)};
      }
      .btn.tile[aria-pressed='true']:not([aria-disabled='true']):where(:not(.frame[data-host-pointer] *)):hover {
        background: ${jy(Vi.lit)};
        box-shadow: ${My(Vi.pressedRingPx,Vi.pressedRing)};
      }
    }
    .btn[aria-disabled='true'] {
      color: #6b7280;
      cursor: default;
    }
    .btn:focus {
      outline: none;
    }
    .btn:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: 1px;
    }
    /* Palettes and tooltips open in the top layer (popover), placed by the toolbars' script. */
    .palette,
    .tip {
      position: fixed;
      inset: auto;
      margin: 0;
      overflow: visible;
    }
    .palette:not(:popover-open, [data-top-open]),
    .tip:not(:popover-open, [data-top-open]) {
      display: none;
    }
    /* Spike menus-top-layer: where the browser has no Popover API, what would open in the top layer is fixed over the page
       instead (top-layer.ts: data-popover, and data-top-open while it shows), where no box of the host's clips it. */
    @supports not selector(:popover-open) {
      [popover] {
        display: none;
      }
    }
    [data-popover] {
      position: fixed;
      inset: auto;
      margin: 0;
      z-index: 2147483647;
    }
    [data-popover]:not([data-top-open]) {
      display: none !important;
    }
    .palette {
      display: flex;
      flex-direction: column;
      gap: 2px;
      padding: 4px;
      color: #1f2328;
      background: #f4f5f7;
      border: 1px solid #c3c8d0;
      border-radius: 6px;
      box-shadow: 0 2px 8px rgb(0 0 0 / 20%);
    }
    /* A menu of the Copy or Paste button (spike copy-as): a column of items, each its format's name and, at its right,
       its keys, as a menu of an application shows them. */
    .menu {
      min-width: 168px;
    }
    /* The zoom's popover (spike polish-feedback, row 13): its menu, a palette's box at its top, where the menu always opened,
       and under it the field that takes a typed percentage, in a box of its own. */
    .palette.zoom-pop {
      gap: 4px;
      padding: 0;
      border: 0;
      background: transparent;
      box-shadow: none;
    }
    .zoom-pop > .menu,
    .zoom-pop > .zoom-field {
      display: flex;
      flex-direction: column;
      gap: 2px;
      padding: 4px;
      background: #f4f5f7;
      border: 1px solid #c3c8d0;
      border-radius: 6px;
      box-shadow: 0 2px 8px rgb(0 0 0 / 20%);
    }
    .zoom-field {
      display: grid;
      grid-template-columns: max-content 1fr;
      align-items: center;
      gap: 2px 8px;
      padding: 6px 8px;
      color: #1f2328;
      font: 12px/1.3 system-ui, sans-serif;
    }
    .zoom-field .zoom-input {
      box-sizing: border-box;
      min-height: 24px;
      width: 6em;
      padding: 2px 6px;
      border: 1px solid #8c959f;
      border-radius: 4px;
      font: inherit;
      color: inherit;
      background: #ffffff;
    }
    .zoom-field .zoom-input:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: 1px;
    }
    .zoom-field .zoom-input[aria-invalid='true'] {
      border-color: #b3261e;
      box-shadow: inset 0 0 0 1px #b3261e;
    }
    .zoom-field .zoom-hint {
      grid-column: 2;
      font-size: 11px;
      color: #4b5563;
    }
    .zoom-field .zoom-hint.error {
      color: #b3261e;
    }
    .menu .item {
      box-sizing: border-box;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      width: 100%;
      min-height: 28px;
      margin: 0;
      padding: 4px 8px;
      border: 0;
      border-radius: 4px;
      background: transparent;
      color: #1f2328;
      font: 12px/1.3 system-ui, sans-serif;
      text-align: left;
      white-space: nowrap;
      cursor: pointer;
    }
    .menu .item kbd {
      font: 11px/1 ui-monospace, 'Cascadia Mono', Consolas, 'DejaVu Sans Mono', 'Courier New', monospace;
      color: #4b5563;
    }
    /* A menu with a check item (spike reports-oct5, R3: the gear's Query mode): a mark's column before every item's name, a
       check in it while the item is on, as the context menu's. */
    .menu .item .mark {
      flex: none;
      width: 12px;
      margin-right: -8px;
      text-align: center;
    }
    .menu .item.marked .item-name {
      flex: 1;
    }
    /* An item with a line saying what it does (B8: Clean Up and Layout), under its name. */
    .menu .item .item-line {
      display: block;
      margin-top: 2px;
      color: #4b5563;
      font-size: 11px;
    }
    .menu .item[aria-disabled='true'] .item-line {
      color: #6b7280;
    }
    @media (hover: hover) {
      .menu .item:not([aria-disabled='true']):where(:not(.frame[data-host-pointer] *)):hover {
        background: #dde3ec;
      }
    }
    .menu .item[aria-disabled='true'] {
      color: #6b7280;
      cursor: default;
    }
    .menu .item:focus {
      outline: none;
    }
    .menu .item:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: -2px;
    }
    /* The abbreviation lookup and Choose Nickname (spike abbreviation-library, L5, L6): a box in the top layer at its point, a
       field over a list of options, each 28 px high for a finger, the active one marked as a menu's focused item is. */
    .lookup {
      position: fixed;
      inset: auto;
      margin: 0;
      box-sizing: border-box;
      width: 280px;
      padding: 6px;
      color: #1f2328;
      background: #fff;
      border: 1px solid #8a94a3;
      border-radius: 6px;
      box-shadow: 0 2px 8px rgb(0 0 0 / 20%);
      font: 13px/1.3 system-ui, sans-serif;
    }
    .lookup:not(:popover-open, [data-top-open]) {
      display: none;
    }
    .lookup-search {
      box-sizing: border-box;
      width: 100%;
      min-height: 28px;
      padding: 2px 8px;
      border: 1px solid #8a94a3;
      border-radius: 4px;
      font: inherit;
      color: #1f2328;
      background: #fff;
    }
    .lookup-search:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: 1px;
    }
    .lookup-list {
      max-height: 240px;
      margin: 6px 0 0;
      padding: 0;
      overflow-y: auto;
      list-style: none;
    }
    .lookup-option {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 8px;
      min-height: 28px;
      padding: 2px 8px;
      border-radius: 4px;
      cursor: pointer;
    }
    .lookup-option[aria-selected='true'] {
      background: #dde3ec;
      box-shadow: inset 0 0 0 1px #1a5fb4;
    }
    .lookup-label {
      overflow-wrap: anywhere;
    }
    .lookup-kind {
      flex: none;
      color: #4b5563;
      font-size: 11px;
    }
    .lookup-status {
      margin: 0;
    }
    .lookup-status span {
      display: block;
      padding: 6px 2px 2px;
    }
    /* The context menu (spike context-menus, C2, C7): rows 32 px high for a finger, a check mark beside a radio item that
       is so, the bond types' icons, a group's heading, separators drawn as borders so that forced colours keep them. */
    .menu.ctx {
      min-width: 210px;
      overflow-y: auto;
    }
    .menu.ctx:focus {
      outline: none;
    }
    .ctx .item {
      justify-content: flex-start;
      gap: 6px;
      min-height: 32px;
    }
    .ctx .item-name {
      flex: 1;
    }
    .ctx .item kbd {
      flex: none;
      margin-left: 12px;
    }
    .ctx .mark {
      flex: none;
      width: 12px;
      text-align: center;
    }
    .ctx .icon {
      flex: none;
      width: 20px;
      height: 20px;
    }
    .ctx .chev {
      flex: none;
      margin-left: 6px;
      color: #4b5563;
    }
    .ctx .sep {
      flex: none;
      margin: 3px 2px;
      border-top: 1px solid #c3c8d0;
    }
    .ctx .grp {
      display: flex;
      flex-direction: column;
      gap: 2px;
    }
    .ctx .heading {
      padding: 2px 8px 0 44px;
      font: 11px/1.3 system-ui, sans-serif;
      color: #4b5563;
    }
    /* A tooltip ignores the pointer: a drag that starts where one opens, over the canvas beside a left
       bar button, draws there and is not swallowed by it (reading 30, the ui-critic's note 4). Its
       transparent border is what forced colours draw, so that it keeps its box there (note 11); it
       takes a pixel off the padding, so that the box is the size it was. */
    .tip {
      pointer-events: none;
      white-space: nowrap;
      padding: 3px 7px;
      border: 1px solid transparent;
      border-radius: 4px;
      background: #1f2328;
      color: #fff;
      font: 12px/1.3 system-ui, sans-serif;
    }
    /* A tooltip's key, in a monospace face where l, I and 1 differ (fix round 4). */
    .tip kbd {
      font: 600 12px/1 ui-monospace, 'Cascadia Mono', Consolas, 'DejaVu Sans Mono', 'Courier New', monospace;
    }
    .stage {
      position: relative;
      grid-row: 2;
      grid-column: 2;
      min-width: 0;
      min-height: 0;
      overflow: hidden;
    }
    .surface {
      position: absolute;
      inset: 0;
      touch-action: none;
      user-select: none;
      -webkit-user-select: none;
      /* A finger held on the canvas opens the sketcher's own menu (spike context-menus, C7), never the browser's callout. */
      -webkit-touch-callout: none;
    }
    .surface > svg {
      display: block;
      width: 100%;
      height: 100%;
      /* The pointer's shape (CANVAS-019): a cross, an arrow over a bond. Set once here, never per
         hover, so moving between atoms and bonds restyles nothing (the perf review, fix round 1). */
      cursor: crosshair;
    }
    .surface .bond-selector {
      cursor: default;
    }
    /* The hand tool (spike view, CANVAS-011): an open hand over the empty canvas, and the move cursor over an atom or a bond,
       which it grabs and moves (spike polish-look, row 3; CANVAS-019); a closed one while the view or what it holds is
       dragged, by the hand or the middle button. By the hit shapes, as the bond's arrow, so that nothing restyles per hover. */
    .surface[data-tool='hand'] > svg {
      cursor: grab;
    }
    .surface[data-tool='hand'] .atom-selector,
    .surface[data-tool='hand'] .bond-selector {
      cursor: move;
    }
    /* A selection tool (spike select): the arrow, over the empty canvas and over what it selects. */
    .surface[data-tool='select'] > svg,
    .surface[data-tool='select'] .bond-selector {
      cursor: default;
    }
    .surface[data-panning] > svg,
    .surface[data-panning] .atom-selector,
    .surface[data-panning] .bond-selector,
    .surface[data-moving] > svg,
    .surface[data-moving] .atom-selector,
    .surface[data-moving] .bond-selector {
      cursor: grabbing;
    }
    /* An edit refused where the pointer is (basics-fixes F3): whatever the tool, a held atom included. */
    .surface[data-refused] > svg,
    .surface[data-refused] .atom-selector,
    .surface[data-refused] .bond-selector {
      cursor: not-allowed;
    }
    /* The canvas's foot (spike charges-checks, K12): the status line on the left and the formula and mass readout on the
       right, over the drawing, so that no box loses height. Only the readout's toggle takes the pointer. Where the two do
       not fit side by side, the foot wraps (spike ui-polish-oct12, row 6; the ui-critic, milestone 3: the status line was
       squeezed into a column beside the readout): the readout keeps the foot's right end, and the status line goes above
       it, as wide as the foot. */
    .foot {
      position: absolute;
      left: 6px;
      right: 6px;
      bottom: 6px;
      z-index: 2;
      display: flex;
      flex-wrap: wrap;
      align-items: flex-end;
      gap: 4px 8px;
      pointer-events: none;
    }
    /* Why an edit is refused (F3), and why an atom is flagged (K10): a polite live region, shown while it says something. */
    .status {
      flex: 0 1 auto;
      min-width: 0;
      box-sizing: border-box;
      padding: 2px 8px;
      border: 1px solid #a15c00;
      border-radius: 4px;
      background: #fff;
      color: #5a3300;
      font: 12px/1.4 system-ui, sans-serif;
      pointer-events: none;
    }
    .status:empty {
      display: none;
    }
    /* The status line and a notice's close button (spike keyboard-a11y, UI-020): the button alone takes the pointer, 24 px. */
    .status-line {
      display: flex;
      flex: 0 1 auto;
      align-items: center;
      gap: 2px;
      min-width: 0;
      pointer-events: none;
    }
    .dismiss {
      flex: none;
      box-sizing: border-box;
      width: 24px;
      height: 24px;
      padding: 0;
      border: 1px solid #a15c00;
      border-radius: 4px;
      background: #fff;
      color: #5a3300;
      font: 16px/1 system-ui, sans-serif;
      cursor: pointer;
      pointer-events: auto;
    }
    .dismiss[hidden] {
      display: none;
    }
    .dismiss:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: 1px;
    }
    /* What screen readers hear and the canvas's description (K3): in the page for them, shown to no one. */
    .sr-only {
      position: absolute;
      width: 1px;
      height: 1px;
      margin: -1px;
      padding: 0;
      overflow: hidden;
      clip-path: inset(50%);
      white-space: nowrap;
      border: 0;
    }
    /* The formula and mass readout (K12): a white box at the foot's right end, its toggle first and its values after it,
       which take no pointer events (the foot's), so that the canvas under them is the canvas's. At most the foot's width
       and one line tall (spike readout-fit). */
    .info {
      display: flex;
      flex: none;
      align-items: center;
      gap: 6px;
      margin-left: auto;
      max-width: 100%;
      box-sizing: border-box;
      padding: 0 8px 0 0;
      border: 1px solid #c3c8d0;
      border-radius: 4px;
      background: #fff;
      color: #1f2328;
      font: 12px/1.4 system-ui, sans-serif;
    }
    .info[hidden],
    .info [hidden] {
      display: none;
    }
    .info[data-collapsed],
    .info[data-folded] {
      padding-right: 0;
    }
    /* With nothing drawn (spike ui-polish-oct12, row 4; INFO-009; the ui-critic, milestone 6: "a lone chevron button"): no
       readout, its toggle neither, and no room taken in the foot; out of sight, so that no one reads its words and its
       toggle is no Tab stop. */
    .info[data-empty] {
      position: absolute;
      right: 0;
      bottom: 0;
      visibility: hidden;
    }
    /* Folded by its box (spike readout-fit, R1): its toggle alone, as collapsed shows it. Its values stay laid out, unseen
       and taking no room, so that their line is still measured. */
    .info[data-folded] {
      gap: 0;
    }
    .info[data-folded] .info-body {
      max-width: 0;
      visibility: hidden;
    }
    /* Its values keep to one line (spike readout-fit, R3): what does not fit is cut with an ellipsis, and the readout's
       tooltip gives the whole text. Its parts are apart by spaces, words to a screen reader, widened to the gaps they had. */
    .info-body {
      flex: 0 1 auto;
      min-width: 0;
      overflow: hidden;
      white-space: nowrap;
      text-overflow: ellipsis;
      word-spacing: 4px;
    }
    .info-body > * {
      word-spacing: normal;
    }
    .info-caption {
      font-weight: 600;
    }
    .info-formula {
      font-weight: 600;
    }
    .info-formula sub,
    .info-formula sup {
      font-size: 75%;
      line-height: 0;
    }
    .info-value {
      font-variant-numeric: tabular-nums;
    }
    .info-label {
      color: #4b5563;
      text-decoration: none;
    }
    .info-note {
      color: #4b5563;
      font-style: italic;
    }
    /* Its items (spike readout-config): each takes the pointer (a click copies it, a right-click opens the readout's menu),
       one line's height, at least 24 px wide, a notation cut with an ellipsis past 16em (its tooltip the whole). */
    .info-item {
      display: inline-block;
      box-sizing: border-box;
      min-width: 24px;
      height: 24px;
      padding: 0 2px;
      border-radius: 3px;
      line-height: 24px;
      vertical-align: middle;
      cursor: copy;
      pointer-events: auto;
      /* A finger held on an item opens the readout's menu: the browser takes none of its moves to scroll. */
      touch-action: none;
    }
    /* A value over the drawing yields the pointer to the atom or bond under it (spike readout-config; info-panel.ts
       #yield): the canvas's, as the line's gaps are. */
    .info-item[data-yields] {
      pointer-events: none;
    }
    .info-notation {
      display: inline-block;
      max-width: 16em;
      overflow: hidden;
      text-overflow: ellipsis;
      vertical-align: top;
      font-family: ui-monospace, 'Cascadia Mono', Consolas, 'DejaVu Sans Mono', monospace;
    }
    @media (hover: hover) {
      .info-item:where(:not(.frame[data-host-pointer] *)):hover {
        background: #dde3ec;
      }
    }
    .info-item:focus {
      outline: none;
    }
    .info-item:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: -2px;
    }
    /* The … button (ruling 3): the items that do not fit, listed by it. */
    .info-more {
      box-sizing: border-box;
      flex: none;
      width: 24px;
      height: 24px;
      margin: 0 -2px 0 0;
      padding: 0;
      border: 0;
      border-radius: 3px;
      background: transparent;
      color: #1f2328;
      font: 600 14px/24px system-ui, sans-serif;
      cursor: pointer;
      pointer-events: auto;
    }
    @media (hover: hover) {
      .info-more:where(:not(.frame[data-host-pointer] *)):hover {
        background: #dde3ec;
      }
    }
    .info-more:focus {
      outline: none;
    }
    .info-more:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: -2px;
    }
    .info-toggle {
      box-sizing: border-box;
      display: flex;
      flex: none;
      align-items: center;
      justify-content: center;
      width: 24px;
      height: 24px;
      margin: 0;
      padding: 0;
      border: 0;
      border-radius: 3px;
      background: transparent;
      color: #1f2328;
      cursor: pointer;
      pointer-events: auto;
    }
    .info-chevron {
      width: 14px;
      height: 14px;
    }
    .info-toggle[aria-expanded='false'] .info-chevron {
      transform: rotate(180deg);
    }
    @media (hover: hover) {
      .info-toggle:where(:not(.frame[data-host-pointer] *)):hover {
        background: #dde3ec;
      }
    }
    .info-toggle:focus {
      outline: none;
    }
    .info-toggle:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: -2px;
    }
    /* No box round the canvas, ever (spike reports-oct14, report 1; the product owner, 2026-10-14: "There never should be
       such a border"): focused by a click, a key, Tab, a menu's Escape or a dialog closing, the canvas draws no outline, and
       neither do its surface, the stage, the frame or the host. Where the keyboard is, the drawing shows: the keyboard's
       hotspot mark on its atom or bond (A11Y-006, H44). The bars' buttons, the menus' items, the fields and the dialogs keep
       their focus rings. */
    .surface > svg:focus,
    .surface > svg:focus-visible,
    .surface:focus,
    .surface:focus-visible,
    .stage:focus,
    .stage:focus-visible,
    .frame:focus,
    .frame:focus-visible {
      outline: none;
    }
    .overlay {
      position: absolute;
      inset: 0;
      width: 100%;
      height: 100%;
      overflow: visible;
      pointer-events: none;
      /* The marks tint the white and leave the drawing's black lines black (fix round 1). */
      mix-blend-mode: multiply;
    }
    /* The selection's chrome (spike rotate-flip): the rotation handle, its centre and the angle, over the marks, drawn
       as they are; the editor takes the presses on them, so the layer takes none. */
    .chrome {
      position: absolute;
      inset: 0;
      width: 100%;
      height: 100%;
      overflow: visible;
      pointer-events: none;
    }
    /* The floating buttons beside a selection (spike rotate-flip, R5): a white row, 32 px high, over the canvas. */
    .actions {
      position: absolute;
      z-index: 2;
      display: flex;
      gap: 4px;
      box-sizing: border-box;
      padding: 1px;
      border: 1px solid #c3c8d0;
      border-radius: 6px;
      background: #ffffff;
      box-shadow: 0 2px 6px rgb(0 0 0 / 18%);
    }
    .actions[hidden] {
      display: none;
    }
    .actions .act {
      box-sizing: border-box;
      display: flex;
      flex: none;
      align-items: center;
      justify-content: center;
      margin: 0;
      padding: 0;
      border: 0;
      border-radius: 4px;
      background: transparent;
      color: #1f2328;
      cursor: pointer;
      touch-action: manipulation;
    }
    .actions .act .icon {
      width: 20px;
      height: 20px;
    }
    @media (hover: hover) {
      .actions .act:where(:not(.frame[data-host-pointer] *)):hover {
        background: #dde3ec;
      }
    }
    .actions .act:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: 1px;
    }
    /* The label editor (spike atoms, ATOM-012): over the atom, in the label's font: the drawing's
       Roboto (draw's editor style) where the system has it, then faces of its metrics (fix round 4). */
    .label-editor {
      position: absolute;
      z-index: 3;
      box-sizing: border-box;
      margin: 0;
      padding: 0 4px;
      border: 1px solid #1c65c0;
      border-radius: 3px;
      background: #fff;
      color: #1f2328;
      font-family: Roboto, 'Helvetica Neue', Arial, sans-serif;
      text-align: center;
      box-shadow: 0 1px 4px rgb(0 0 0 / 25%);
    }
    .label-editor:focus-visible,
    .label-editor:focus {
      outline: 2px solid #1a5fb4;
      outline-offset: 1px;
    }
    /* Text refused (fix round 3): the input marked, and the message under it. */
    .label-editor[aria-invalid='true'] {
      border-color: #b3261e;
    }
    /* A group's tooltip (spike sgroups, ABBR-043): beside the pointer, over the canvas, taking no pointer events. */
    .sgroup-tooltip {
      position: absolute;
      z-index: 4;
      box-sizing: border-box;
      width: max-content;
      max-width: min(320px, calc(100% - 8px));
      padding: 3px 8px;
      border: 1px solid #c3c8d0;
      border-radius: 4px;
      background: #fff;
      color: #1f2328;
      font: 12px/1.4 system-ui, sans-serif;
      box-shadow: 0 1px 4px rgb(0 0 0 / 20%);
      pointer-events: none;
      overflow-wrap: anywhere;
    }
    /* The S-Group Properties dialog's placement radios (spike sgroups, G4): 24 px targets in a row. */
    .fields .radios {
      grid-column: 2;
      display: flex;
      flex-wrap: wrap;
      gap: 2px 12px;
    }
    .fields .radio {
      display: inline-flex;
      align-items: center;
      gap: 2px;
      min-height: 24px;
    }
    .fields .radio input {
      width: 18px;
      height: 18px;
      min-height: 0;
      margin: 3px;
      padding: 0;
      border: none;
    }
    .fields .sg-placement-label {
      grid-column: 1;
    }
    .dialog.sgroup .fields textarea {
      min-height: 3.5em;
    }
    .dialog.sgroup button[aria-disabled='true'] {
      opacity: 0.55;
      cursor: default;
    }
    .label-editor-message {
      position: absolute;
      z-index: 3;
      box-sizing: border-box;
      width: max-content;
      max-width: 100%;
      padding: 2px 6px;
      border: 1px solid #b3261e;
      border-radius: 3px;
      background: #fff;
      color: #b3261e;
      font: 12px/1.4 system-ui, sans-serif;
      box-shadow: 0 1px 4px rgb(0 0 0 / 25%);
    }
    /* The dialogs (spike atoms, UI-018): modal, in the top layer. */
    .dialog {
      box-sizing: border-box;
      max-width: calc(100vw - 16px);
      max-height: calc(100vh - 16px);
      overflow: auto;
      padding: 12px 16px 14px;
      border: 1px solid #c3c8d0;
      border-radius: 8px;
      background: #fff;
      color: #1f2328;
      font: 13px/1.4 system-ui, sans-serif;
      box-shadow: 0 8px 28px rgb(0 0 0 / 25%);
    }
    .dialog::backdrop {
      background: rgb(15 23 42 / 25%);
    }
    .dialog h2 {
      margin: 0;
      font-size: 15px;
      font-weight: 600;
    }
    .dialog-head {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      margin-bottom: 10px;
    }
    .dialog button {
      box-sizing: border-box;
      min-width: 24px;
      min-height: 24px;
      padding: 4px 12px;
      border: 1px solid #8a94a3;
      border-radius: 4px;
      background: #f4f5f7;
      color: #1f2328;
      font: inherit;
      cursor: pointer;
    }
    .dialog button.primary {
      border-color: #1a5fb4;
      background: #1a5fb4;
      color: #fff;
    }
    .dialog button:disabled {
      border-color: #c3c8d0;
      background: #f4f5f7;
      color: #5c6470;
      cursor: default;
    }
    .dialog button.close {
      width: 28px;
      height: 28px;
      padding: 0;
      border-color: transparent;
      background: transparent;
      font-size: 18px;
      line-height: 1;
    }
    .dialog button:focus-visible,
    .dialog input:focus-visible,
    .dialog select:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: 1px;
    }
    /* A dialog's buttons stay in view while it scrolls (spike polish-look, row 11: at 200%, Atom Properties in query mode
       outgrew a 432 px window, and Apply and Cancel scrolled away under its Query section): the foot sticks to the dialog's
       bottom edge, over the dialog's bottom padding (its 14 px, taken into the foot so that nothing moves where the dialog
       does not scroll; the sticky offset reaches into the padding, which Chrome's sticky box leaves out), on white, with
       what scrolls passing under it. The dialog keeps a field that keyboard focus moves to clear of it (scroll-padding). */
    .dialog {
      scroll-padding-bottom: 56px;
    }
    .dialog-foot {
      position: sticky;
      bottom: -14px;
      z-index: 1;
      display: flex;
      justify-content: flex-end;
      gap: 8px;
      margin: 4px 0 -14px;
      padding: 8px 0 14px;
      background: #fff;
    }
    .fields {
      display: grid;
      grid-template-columns: max-content minmax(10em, 1fr);
      align-items: center;
      gap: 6px 10px;
    }
    .fields label {
      grid-column: 1;
    }
    .fields input,
    .fields select,
    .fields output,
    .fields textarea,
    .fields .with-button {
      grid-column: 2;
    }
    /* The Paste Special dialog (spike copy-as): its field wide enough for a molfile's lines, in a face where they line up. */
    .dialog.paste-special {
      width: min(540px, calc(100vw - 16px));
    }
    .paste-special .fields label[for='ps-text'] {
      align-self: start;
      padding-top: 4px;
    }
    .fields textarea {
      box-sizing: border-box;
      width: 100%;
      min-height: 8em;
      padding: 4px 6px;
      border: 1px solid #8a94a3;
      border-radius: 4px;
      font: 12px/1.35 ui-monospace, 'Cascadia Mono', Consolas, 'DejaVu Sans Mono', 'Courier New', monospace;
      color: #1f2328;
      background: #fff;
      resize: vertical;
    }
    /* The Settings dialog (spike settings-touch, K6): a check box before its label, a list after its label. Compact since
       spike stereo-wording (the product owner, 2026-10-20: "it looks overly tall"): sections under headings, in two columns
       where the window is wide (the readout's the second), else one; rows of 26 px, the targets still 24 px. */
    .settings-fields {
      display: flex;
      flex-direction: column;
      gap: 10px;
      min-width: min(300px, calc(100vw - 48px));
      /* Taller than the window (spike readout-config): its fields scroll, its buttons stay. */
      max-height: calc(100vh - 140px);
      overflow-y: auto;
    }
    .st-column {
      display: flex;
      flex-direction: column;
      gap: 10px;
      min-width: 0;
    }
    .st-section {
      display: flex;
      flex-direction: column;
      gap: 2px;
      min-width: 0;
    }
    .dialog .st-heading {
      margin: 0 0 2px;
      padding-bottom: 2px;
      border-bottom: 1px solid #dde3ec;
      color: #4b5563;
      font-size: 12px;
      font-weight: 600;
    }
    /* Wide: the readout's section a second column, wide enough for its items in three columns. */
    @media (min-width: 720px) {
      .dialog.settings {
        width: min(840px, calc(100vw - 16px));
      }
      .settings-fields {
        display: grid;
        grid-template-columns: minmax(0, 3fr) minmax(0, 4fr);
        align-items: start;
        column-gap: 24px;
      }
    }
    .check-row,
    .select-row {
      display: flex;
      align-items: center;
      gap: 8px;
      min-height: 26px;
    }
    /* A list after its label, at the row's end; under it where both do not fit (a longer language), so that neither is
       cut. */
    .select-row {
      flex-wrap: wrap;
      justify-content: space-between;
      row-gap: 2px;
    }
    .select-row select {
      max-width: 100%;
      margin-left: auto;
    }
    .check-row input {
      flex: none;
      width: 18px;
      height: 18px;
      margin: 3px;
    }
    /* A label longer than its row (a longer language) wraps, a word longer than the row breaking (spike stereo-wording: in the
       pseudo-locale a readout item's label ran into the next column). */
    .check-row label,
    .select-row label {
      min-width: 0;
      overflow-wrap: anywhere;
    }
    /* The readout's items (spike readout-config; compact since spike stereo-wording): under their heading, each group's
       name over its check boxes, which flow in columns; no frame. */
    .readout-items {
      display: flex;
      flex-direction: column;
      gap: 4px;
      margin-top: 2px;
    }
    .readout-items-heading {
      font-weight: 600;
    }
    .readout-group {
      display: flex;
      flex-direction: column;
      min-width: 0;
    }
    .readout-group-name {
      color: #4b5563;
      font-size: 12px;
    }
    .readout-group-items {
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(132px, 1fr));
      column-gap: 8px;
      padding-left: 4px;
    }
    .readout-group .check-row {
      min-height: 24px;
      gap: 4px;
    }
    .settings select {
      min-height: 28px;
      padding: 2px 4px;
      border: 1px solid #8a94a3;
      border-radius: 4px;
      background: #fff;
      color: #1f2328;
      font: inherit;
    }
    .settings input:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: 1px;
    }
    .dialog button.reset[aria-disabled='true'] {
      border-color: #c3c8d0;
      color: #5c6470;
      cursor: default;
    }
    /* The About box (spike later-quick-2, UI-027): its facts in two columns, the notices under them. */
    .dialog.about {
      width: min(420px, calc(100vw - 16px));
    }
    .about .facts {
      display: grid;
      grid-template-columns: auto 1fr;
      gap: 4px 12px;
      margin: 8px 0;
    }
    .about dt {
      color: #4b5563;
    }
    .about dd {
      margin: 0;
      overflow-wrap: anywhere;
    }
    .about .notices {
      margin: 8px 0 0;
      font-size: 12px;
      line-height: 1.4;
    }
    /* The shortcuts sheet (K7): its tables scroll inside the dialog; each key in a monospace kbd. */
    .dialog.keys {
      width: min(560px, calc(100vw - 16px));
    }
    .keys .sheet {
      max-height: min(60vh, 520px);
      overflow: auto;
    }
    .keys .sheet:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: 1px;
    }
    .keys table {
      width: 100%;
      margin: 0 0 10px;
      border-collapse: collapse;
    }
    .keys caption {
      padding: 6px 0 4px;
      font-weight: 600;
      text-align: left;
    }
    .keys th,
    .keys td {
      padding: 3px 6px;
      border-bottom: 1px solid #e3e6ea;
      text-align: left;
      vertical-align: top;
      font-weight: normal;
    }
    .keys th[scope='col'] {
      color: #4b5563;
      font-size: 12px;
    }
    .keys kbd {
      display: inline-block;
      padding: 0 4px;
      border: 1px solid #c3c8d0;
      border-radius: 3px;
      background: #f4f5f7;
      font: 12px/1.5 ui-monospace, 'Cascadia Mono', Consolas, 'DejaVu Sans Mono', 'Courier New', monospace;
      white-space: nowrap;
    }
    .keys .or {
      color: #4b5563;
    }
    /* Bond Properties' type list (spike context-menus, C8): a list box, every type showing. */
    .fields select.list {
      min-width: 12em;
      padding: 2px;
    }
    .fields textarea[aria-invalid='true'] {
      border-color: #b3261e;
      box-shadow: inset 0 0 0 1px #b3261e;
    }
    .dialog textarea:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: 1px;
    }
    .fields input,
    .fields select {
      box-sizing: border-box;
      min-height: 26px;
      padding: 2px 6px;
      border: 1px solid #8a94a3;
      border-radius: 4px;
      font: inherit;
      color: #1f2328;
      background: #fff;
    }
    .fields input[aria-invalid='true'],
    .fields select[aria-invalid='true'] {
      border-color: #b3261e;
      box-shadow: inset 0 0 0 1px #b3261e;
    }
    .with-button {
      display: flex;
      gap: 6px;
    }
    .with-button input {
      flex: 1;
      min-width: 4em;
    }
    .fields .hint,
    .fields .error {
      grid-column: 2;
      margin-top: -4px;
      font-size: 12px;
      color: #4b5563;
    }
    .fields .error {
      color: #b3261e;
    }
    .fields .error[hidden] {
      display: none;
    }
    .ptable {
      display: grid;
      grid-template-columns: repeat(18, 30px);
      grid-template-rows: repeat(7, 30px) 10px 30px 30px;
      gap: 2px;
    }
    .dialog button.pt-el {
      width: 30px;
      height: 30px;
      min-width: 0;
      padding: 0;
      border-color: #c3c8d0;
      color: var(--atom-colour);
      font-size: 12px;
      font-weight: 600;
    }
    /* The chosen cell (ATOM-033, P3): a light fill and a ring, not today's dark blue fill with white text, which
       no symbol colour reads on. */
    .dialog button.pt-el[aria-pressed='true'] {
      border-color: ${jy(R.ring)};
      background: ${jy(R.fill)};
      box-shadow: ${My(R.ringPx,R.ring)};
    }
    /* The categories (ATOM-033, P3). The nonmetals and the halogens hold the coloured symbols, which keep 4.5:1
       only on white or a fill nearly as light, so their cells are white; the other eight keep their tints,
       where black reads. Every category is also told by a stripe along the foot of its cells, in a deeper
       shade of its hue, so that the two white categories are told apart by more than a faint tint. */
    ${Ny}
    /* Spike query-edit (Q6, Q8): the periodic table's modes and Generics row, Atom Properties' and Bond Properties' query
       fields: check boxes and radio buttons 18 px within 24 px targets, the query fields two pairs a row where there is room. */
    .q-check {
      display: inline-flex;
      align-items: center;
      gap: 4px;
      min-height: 24px;
    }
    .q-check input {
      width: 18px;
      height: 18px;
      min-height: 0;
      margin: 3px;
      padding: 0;
    }
    .pt-modes {
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 2px 14px;
      margin: 0 0 8px;
      padding: 0;
      border: none;
    }
    .pt-modes legend {
      float: left;
      margin-right: 10px;
      padding: 0;
      font-weight: 600;
    }
    .pt-generics {
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 2px;
      margin-top: 8px;
    }
    .pt-generics[data-hidden] {
      visibility: hidden;
    }
    .pt-generics-label {
      margin-right: 8px;
      font-weight: 600;
    }
    .dialog button.pt-generic {
      min-width: 36px;
      height: 30px;
      padding: 0 4px;
      border-color: #c3c8d0;
      background: #fff;
      font-weight: 600;
    }
    .dialog button.pt-generic[aria-pressed='true'] {
      border-color: ${jy(R.ring)};
      background: ${jy(R.fill)};
      box-shadow: ${My(R.ringPx,R.ring)};
    }
    .dialog button.pt-generic[aria-disabled='true'] {
      opacity: 0.55;
      cursor: default;
    }
    .dialog .hint {
      margin: 6px 0 0;
      font-size: 12px;
      color: #4b5563;
    }
    .query-section {
      margin: 10px 0 0;
      padding: 6px 8px 8px;
      border: 1px solid #c3c8d0;
      border-radius: 4px;
    }
    .query-section legend {
      padding: 0 4px;
      font-weight: 600;
    }
    .query-fields {
      display: grid;
      grid-template-columns: max-content minmax(6em, 1fr) max-content minmax(6em, 1fr);
      align-items: center;
      gap: 6px 10px;
    }
    .query-fields .q-check {
      grid-column: span 2;
    }
    @media (max-width: 520px) {
      .query-fields {
        grid-template-columns: max-content minmax(6em, 1fr);
      }
    }
    .query-fields select,
    .query-custom input[type='text'],
    .fields #bp-topology,
    .fields #bp-smarts {
      box-sizing: border-box;
      min-height: 26px;
      padding: 2px 6px;
      border: 1px solid #8a94a3;
      border-radius: 4px;
      font: inherit;
      color: #1f2328;
      background: #fff;
    }
    .query-custom {
      display: grid;
      grid-template-columns: max-content minmax(10em, 1fr);
      align-items: center;
      gap: 6px 10px;
      margin-top: 8px;
    }
    .query-custom .error,
    .fields #bp-error {
      grid-column: 2;
      font-size: 12px;
      color: #b3261e;
    }
    .query-custom .error[hidden] {
      display: none;
    }
    /* Atom Properties' Apply refused by the engine (spike polish-feedback, row 1): said above the buttons. */
    .dialog .refused {
      margin: 10px 0 0;
      font-size: 12px;
      color: #b3261e;
    }
    .dialog .refused[hidden] {
      display: none;
    }
    .query-section input[aria-invalid='true'],
    .fields #bp-smarts[aria-invalid='true'] {
      border-color: #b3261e;
      box-shadow: inset 0 0 0 1px #b3261e;
    }
    .query-section select:disabled,
    .query-section input:disabled,
    .fields select:disabled,
    .fields input:disabled {
      background: #f4f5f7;
      color: #5c6470;
    }
    .fields .bp-custom {
      grid-column: 1;
    }
    /* Bond Properties on crossed double bonds (spike ui-polish-oct12, row 5): under Type, that their E/Z is unknown. */
    .fields .bp-ez {
      grid-column: 2;
      margin: 0;
      color: #4b535e;
      font-size: 12px;
    }
    .pt-gap {
      display: flex;
      align-items: center;
      justify-content: center;
      color: #4b5563;
    }
    .details {
      min-height: 1.4em;
      margin: 10px 0 0;
    }
    /* The R-Group dialog's toggles, R1 to R32 in rows of 8, and the Attachment Points dialog's check boxes
       (spike rgroups). A pressed toggle has a light fill and a ring, as the periodic table's chosen cell. */
    .dialog .rg-hint {
      margin: 0 0 8px;
      font-size: 12px;
      color: #4b5563;
    }
    /* "Remove abbreviation?" (spike abbreviations, P5): the question, then the hint, at most a phone's width. */
    .dialog .dissolve-question {
      margin: 0 0 6px;
      max-width: min(360px, calc(100vw - 48px));
    }
    .rgrid {
      display: grid;
      grid-template-columns: repeat(8, 44px);
      gap: 4px;
    }
    .dialog button.r-toggle {
      min-width: 0;
      height: 30px;
      padding: 0;
      border-color: #c3c8d0;
      background: #fff;
      font-weight: 600;
    }
    .dialog button.r-toggle[aria-pressed='true'] {
      border-color: ${jy(R.ring)};
      background: ${jy(R.fill)};
      box-shadow: ${My(R.ringPx,R.ring)};
    }
    .checks {
      display: flex;
      flex-direction: column;
      gap: 8px;
    }
    .checks .check {
      display: flex;
      align-items: center;
      gap: 8px;
      min-height: 24px;
    }
    .checks input {
      width: 18px;
      height: 18px;
      margin: 0;
    }
    .checks input:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: 2px;
    }
    .dialog .hint {
      margin: 0 0 8px;
      font-size: 12px;
      color: #4b5563;
    }
    /* The stereochemistry panel (spike stereo-panel; the design's mockup-side-panel.html): in the stage, the canvas giving
       it room (the stage's data-panel and its --panel-w and --panel-h). On the right the full height, 216 or 288 px wide; at
       the bottom the full width, its height clamped, the list scrolling; over the whole stage where neither fits. Its width
       and its compact rules follow the sketcher's box, never the bars' density option (spike compact-default, the
       orchestrator's ruling of 2026-10-14): compact (216 px) in a box under 720 by 600 px, the frame's
       data-panel-density, regular (288 px) otherwise. */
    .stage[data-panel='right'] .surface {
      right: var(--panel-w, 0px);
    }
    .stage[data-panel='bottom'] .surface {
      bottom: var(--panel-h, 0px);
    }
    /* Spike panel-narrow (F4; the product owner: "the preview drawn pixels are off"): the editor's layers take the canvas's
       box as the surface does, so that a preview mapped through the canvas's viewBox (overlay.ts showDrawn) lies where the
       edit draws, and every mark keeps the canvas's px; their left and top edges stay (an RTL page's too). */
    /* While the drawing is moved to show a view not drawn yet (spike perf-big), the selection's floating buttons wait for the
       drawing: they are placed on it. */
    .stage[data-view-moving] .actions {
      visibility: hidden;
    }
    .stage[data-panel='right'] .overlay,
    .stage[data-panel='right'] .chrome {
      right: auto;
      width: calc(100% - var(--panel-w, 0px));
    }
    .stage[data-panel='bottom'] .overlay,
    .stage[data-panel='bottom'] .chrome {
      bottom: auto;
      height: calc(100% - var(--panel-h, 0px));
    }
    .stage[data-panel='right'] .foot {
      right: calc(var(--panel-w, 0px) + 6px);
    }
    .stage[data-panel='bottom'] .foot {
      bottom: calc(var(--panel-h, 0px) + 6px);
    }
    .sp {
      position: absolute;
      z-index: 3;
      display: flex;
      flex-direction: column;
      box-sizing: border-box;
      min-width: 0;
      min-height: 0;
      background: #ffffff;
      color: #1f2328;
      font: 12px/1.3 system-ui, sans-serif;
      /* Spike panel-narrow (C4): where the rows' contents start (their 3 px selection border and padding), which the
         title, the summary, the foot and the empty message start at too. */
      --sp-inset: 13px;
    }
    .frame[data-panel-density='compact'] .sp {
      --sp-inset: 10px;
    }
    .sp[hidden] {
      display: none;
    }
    .sp[data-place='right'] {
      top: 0;
      right: 0;
      bottom: 0;
      width: var(--panel-w, 216px);
      border-left: 1px solid #c3c8d0;
      box-shadow: -2px 0 10px rgb(15 23 42 / 8%);
    }
    .sp[data-place='bottom'] {
      left: 0;
      right: 0;
      bottom: 0;
      height: var(--panel-h, 172px);
      border-top: 1px solid #c3c8d0;
      box-shadow: 0 -2px 10px rgb(15 23 42 / 8%);
    }
    .sp[data-place='cover'],
    .sp:not([data-place]) {
      inset: 0;
    }
    .sp button {
      font: inherit;
      color: inherit;
    }
    /* Spike panel-narrow (F1, C1): Close in the head's corner, the head keeping its room; what does not fit beside the
       title (the commands, at the bottom) wraps below it, inside the panel, at any width the panel takes. */
    .sp-head {
      position: relative;
      display: flex;
      flex: none;
      flex-wrap: wrap;
      align-items: center;
      gap: 4px;
      box-sizing: border-box;
      min-height: 36px;
      padding-block: 5px;
      padding-inline: var(--sp-inset) 36px;
      border-bottom: 1px solid #e3e6ea;
    }
    .frame[data-panel-density='compact'] .sp-head {
      min-height: 32px;
      padding-block: 3px;
    }
    .sp-head .sp-close {
      position: absolute;
      top: 6px;
      inset-inline-end: 6px;
    }
    .frame[data-panel-density='compact'] .sp-head .sp-close {
      top: 4px;
    }
    .sp-cmds {
      display: flex;
      flex: 0 1 auto;
      flex-wrap: wrap;
      gap: 2px;
      min-width: 0;
    }
    /* Remove stereo and Recompute wedges as icon buttons (spike reports-oct14, report 5): 24 px, their icons 18 px, greyed
       while they cannot act. */
    .sp .sp-cmd .icon {
      display: block;
      width: 18px;
      height: 18px;
    }
    .sp .sp-cmd[aria-disabled='true'] {
      color: #8a94a3;
      cursor: default;
    }
    .sp[data-place='bottom'] .sp-sum {
      flex-direction: row;
      flex-wrap: wrap;
      column-gap: 12px;
    }
    .sp-head h2 {
      flex: 1 1 auto;
      min-width: 0;
      margin: 0;
      font-size: 13px;
      font-weight: 600;
      line-height: 24px;
      overflow-wrap: anywhere;
    }
    .sp .ibtn {
      box-sizing: border-box;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      width: 24px;
      height: 24px;
      margin: 0;
      padding: 0;
      border: 0;
      border-radius: 4px;
      background: transparent;
      cursor: pointer;
    }
    .sp-sum {
      display: flex;
      flex: none;
      flex-direction: column;
      gap: 2px;
      padding: 6px 10px 7px var(--sp-inset);
      border-bottom: 1px solid #e3e6ea;
    }
    .sp .muted {
      color: #4b5563;
    }
    /* The count's line (spike reports-oct14, report 5; spike mac-fonts, the product owner's "second line when needed"):
       its words flow as text on one baseline, on one line where they fit, else on two, the second starting where the
       first does. The elements' last word, the unspecified marker and the unspecified count's first word stay together
       (.sp-keep), so the marker is always on the first line; its box is as tall as the line, its top on the line's, so it
       is centred on the line whatever the font, at the rows' size. */
    .sp-count {
      display: block;
    }
    .sp-keep {
      white-space: nowrap;
    }
    .sp-warn {
      color: #b3261e;
    }
    .sp-pm {
      display: inline-flex;
      align-items: center;
      height: 1lh;
      margin-right: 4px;
      vertical-align: top;
    }
    .sp .pm {
      flex: none;
      stroke: #d1242f;
    }
    .sp-empty {
      flex: 1;
      margin: 0;
      padding: 10px 10px 10px var(--sp-inset);
      color: #4b5563;
    }
    /* Spike panel-narrow (F3, C3): one set of columns, each row a subgrid of the list's: the name with its marker, the
       configuration (as wide as the widest), the group chip (every chip as wide as the widest; a row without one leaves
       its cell empty). A row stays one box: its role, focus outline, selection border and hover. */
    .sp-list {
      display: grid;
      grid-template-columns: minmax(0, 1fr) auto auto;
      align-content: start;
      column-gap: 6px;
      flex: 1;
      min-height: 0;
      overflow-y: auto;
      padding: 2px 0;
      touch-action: pan-y;
      overscroll-behavior: contain;
    }
    .frame[data-panel-density='compact'] .sp-list {
      column-gap: 4px;
    }
    .frame[data-panel-density='compact'] .sp-row .who {
      gap: 4px;
    }
    .sp-row {
      display: grid;
      grid-column: 1 / -1;
      grid-template-columns: subgrid;
      align-items: center;
      box-sizing: border-box;
      min-height: 40px;
      padding: 0 8px 0 calc(var(--sp-inset) - 3px);
      border-left: 3px solid transparent;
      cursor: default;
    }
    .frame[data-panel-density='compact'] .sp-row {
      min-height: 32px;
      padding-right: 6px;
    }
    .sp-row .who {
      display: flex;
      min-width: 0;
      align-items: center;
      gap: 6px;
    }
    .sp-row .cfg {
      display: inline-flex;
      justify-self: start;
    }
    .sp-row .gcell {
      display: flex;
      justify-self: stretch;
    }
    /* Spike panel-narrow (F2, C2): the name keeps its room, whole where the row has it, else cut with an ellipsis (its
       title the whole name); the marker after it never shrinks, so nothing is drawn over the name. */
    .sp-row .nm {
      display: flex;
      flex: 0 1 auto;
      flex-direction: column;
      min-width: 0;
      overflow: hidden;
      line-height: 1.25;
    }
    .sp-row .name {
      overflow: hidden;
      font-weight: 600;
      text-overflow: ellipsis;
    }
    .sp-row .kword {
      overflow: hidden;
      font-size: 11px;
      color: #5c6470;
      white-space: nowrap;
      text-overflow: ellipsis;
    }
    .frame[data-panel-density='compact'] .sp[data-place='right'] .sp-row .kword {
      display: none;
    }
    .sp-row .kind {
      display: inline-flex;
      flex: none;
      color: #4b5563;
    }
    .sp-row .pm-wrap {
      display: inline-flex;
      flex: none;
    }
    .sp-row.hover {
      background: #eaf3ff;
    }
    .sp-row[aria-selected='true'] {
      border-left-color: #1b7a3a;
      background: #e3f5e6;
    }
    .sp-row:focus {
      outline: none;
    }
    .sp-row:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: -2px;
    }
    @media (hover: hover) {
      .sp-row:where(:not(.frame[data-host-pointer] *)):hover:not([aria-selected='true']) {
        background: #f4f5f7;
      }
      .sp .seg button:where(:not(.frame[data-host-pointer] *)):not([aria-pressed='true']):hover,
      .sp .grp:where(:not([aria-disabled='true'])):hover,
      .sp .ibtn:where(:not([aria-disabled='true'])):hover,
      .sp-menu .item:hover {
        background: #dde3ec;
      }
    }
    .sp .seg {
      display: inline-flex;
      flex: none;
      overflow: hidden;
      border: 1px solid #8a94a3;
      border-radius: 4px;
      background: #ffffff;
    }
    .sp .seg button {
      box-sizing: border-box;
      min-width: 24px;
      height: 24px;
      margin: 0;
      padding: 0 6px;
      border: 0;
      border-left: 1px solid #c3c8d0;
      background: #ffffff;
      font-weight: 600;
      cursor: pointer;
    }
    .frame[data-panel-density='compact'] .sp .seg button {
      padding: 0 5px;
    }
    .sp .seg button:first-child {
      border-left: 0;
    }
    .sp .seg button[aria-pressed='true'] {
      background: #cfe0fb;
      box-shadow: inset 0 0 0 1px #1c65c0;
    }
    .sp .seg button[aria-pressed='true'].unset {
      background: #fdecea;
      box-shadow: inset 0 0 0 1px #d1242f;
      color: #b3261e;
    }
    .sp .grp {
      box-sizing: border-box;
      display: inline-flex;
      flex: none;
      align-items: center;
      justify-content: space-between;
      gap: 2px;
      width: 100%;
      min-width: 52px;
      height: 24px;
      margin: 0;
      padding: 0 3px 0 7px;
      border: 1px solid #8a94a3;
      border-radius: 4px;
      background: #ffffff;
      cursor: pointer;
      white-space: nowrap;
    }
    .frame[data-panel-density='compact'] .sp .grp,
    .frame[data-panel-density='compact'] .sp .grp-empty {
      min-width: 46px;
    }
    .frame[data-panel-density='compact'] .sp .grp {
      padding-left: 5px;
    }
    .sp .grp[aria-disabled='true'] {
      border-color: #c3c8d0;
      background: #f4f5f7;
      color: #5c6470;
      cursor: default;
    }
    .sp .grp .none {
      color: #4b5563;
    }
    .sp .grp-empty {
      display: block;
      min-width: 52px;
    }
    .sp .seg button:focus,
    .sp .grp:focus,
    .sp .ibtn:focus {
      outline: none;
    }
    .sp .seg button:focus-visible,
    .sp .grp:focus-visible,
    .sp .ibtn:focus-visible {
      position: relative;
      z-index: 1;
      outline: 2px solid #1a5fb4;
      outline-offset: 1px;
    }
    /* The group menu opens in the page's top layer (spike menus-top-layer), placed in the window beside its chip by the
       panel's script (stereo-panel.ts), whole where the window has room: past the panel's and the sketcher's edges. */
    .sp-menu {
      position: fixed;
      inset: auto;
      margin: 0;
      display: flex;
      flex-direction: column;
      gap: 2px;
      box-sizing: border-box;
      min-width: 196px;
      overflow-x: hidden;
      overflow-y: auto;
      color: #1f2328;
      font: 12px/1.3 system-ui, sans-serif;
      padding: 4px;
      border: 1px solid #c3c8d0;
      border-radius: 6px;
      background: #f4f5f7;
      box-shadow: 0 2px 8px rgb(0 0 0 / 20%);
    }
    .sp-menu:not(:popover-open, [data-top-open]) {
      display: none;
    }
    .sp-menu .heading {
      padding: 3px 8px 1px 28px;
      font-size: 11px;
      color: #4b5563;
    }
    .sp-menu .item {
      display: flex;
      align-items: center;
      gap: 6px;
      box-sizing: border-box;
      width: 100%;
      min-height: 28px;
      margin: 0;
      padding: 4px 8px;
      border: 0;
      border-radius: 4px;
      background: transparent;
      text-align: left;
      white-space: nowrap;
      cursor: pointer;
    }
    .sp-menu .item:focus {
      outline: none;
      background: #dde3ec;
    }
    .sp-menu .item:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: -2px;
    }
    .sp-menu .mark {
      flex: none;
      width: 14px;
      text-align: center;
    }
    .sp-menu .hint {
      margin-left: auto;
      padding-left: 12px;
      color: #4b5563;
    }
    .sp-menu .sep {
      margin: 3px 2px;
      border-top: 1px solid #c3c8d0;
    }
    /* Spike stereo-wording: a group's item says what the group means on a second, quieter line, which wraps within the
       width the items' names give the menu (its own width counted as none: contain), so a long wording never widens it. */
    .sp-menu {
      min-width: 280px;
    }
    .sp-menu .item.has-meaning {
      align-items: flex-start;
      padding-top: 5px;
      padding-bottom: 5px;
    }
    .sp-menu .item-body {
      display: flex;
      flex-direction: column;
      flex: 1 1 auto;
      min-width: 0;
    }
    .sp-menu .item-line {
      display: flex;
      align-items: center;
    }
    .sp-menu .meaning {
      contain: inline-size;
      margin-top: 1px;
      color: #4b5563;
      font-size: 11px;
      line-height: 1.3;
      white-space: normal;
      /* A word longer than the line (a compound word, the pseudo-locale's tildes) breaks rather than overflows. */
      overflow-wrap: anywhere;
    }
    /* The rows' controls' tooltip (spike stereo-wording): a segment's name and meaning on one line, a chip's name over its
       group's meaning; it wraps where the bars' tooltips do not, so a long wording stays beside its control. */
    .sp-row-tip {
      max-width: 280px;
      white-space: normal;
      overflow-wrap: anywhere;
    }
    .sp-row-tip .tip-name,
    .sp-row-tip .tip-meaning {
      display: block;
    }
    .sp-row-tip .tip-meaning {
      margin-top: 1px;
      color: #d0d7de;
    }
    /* The Structure Library (spike templates): a larger dialog whose groups scroll under its search. Its
       cards sit in a grid of 120 px columns (structure-library.ts, CARD_COLUMN and CARD_GAP); a template
       too large for one at its readable size spans more columns (reading 6). */
    .dialog.library {
      width: min(760px, calc(100vw - 16px));
      height: min(640px, calc(100vh - 16px));
      overflow: hidden;
    }
    .dialog.library[open] {
      display: flex;
      flex-direction: column;
    }
    /* The library while its module loads (milestone-4-fixes, F10): its dialog at once, saying so. */
    .sl-wait {
      margin: 24px 0;
      text-align: center;
    }
    .sl-search {
      display: flex;
      align-items: center;
      gap: 8px;
    }
    .sl-search input {
      flex: 1;
      box-sizing: border-box;
      min-width: 0;
      min-height: 28px;
      padding: 2px 8px;
      border: 1px solid #8a94a3;
      border-radius: 4px;
      font: inherit;
      color: #1f2328;
      background: #fff;
    }
    .sl-search input::placeholder {
      color: #5c6470;
    }
    .sl-search input:focus-visible {
      outline: 2px solid #1a5fb4;
      outline-offset: 1px;
    }
    /* The library's tabs (spike abbreviation-library, L1): a row under the search, the chosen one underlined in the accent,
       as Ketcher's are; a flat tab's cards in the same grid as a group's. */
    .sl-tabs {
      display: flex;
      flex-wrap: wrap;
      gap: 2px;
      margin-top: 10px;
      border-bottom: 1px solid #d0d5dc;
    }
    .dialog button.sl-tab {
      min-height: 32px;
      margin-bottom: -1px;
      padding: 4px 12px;
      border: 0;
      border-bottom: 3px solid transparent;
      border-radius: 0;
      background: transparent;
    }
    .dialog button.sl-tab[aria-selected='true'] {
      border-bottom-color: #1a5fb4;
      font-weight: 600;
    }
    .sl-caption sub,
    .sl-caption sup {
      font-size: 75%;
      line-height: 0;
    }
    .sl-status p {
      margin: 12px 0 0;
    }
    /* The tab's panel and, under it while a search is typed, what the other tabs hold (spike ui-polish-oct12, row 3), which
       scroll together under the search. */
    .sl-body {
      flex: 1;
      min-height: 0;
      margin-top: 10px;
      overflow: auto;
    }
    .sl-elsewhere {
      margin-top: 12px;
      padding-top: 4px;
      border-top: 2px solid #d0d5dc;
    }
    .sl-elsewhere h3 {
      margin: 4px 0;
      font-size: 12px;
      font-weight: 600;
      color: #4b535e;
    }
    .sl-group[hidden] {
      display: none;
    }
    .sl-group h3,
    .sl-group h4 {
      margin: 0;
      font: inherit;
    }
    .dialog button.sl-header {
      display: flex;
      align-items: center;
      gap: 6px;
      width: 100%;
      min-height: 32px;
      padding: 4px 6px;
      border: 0;
      border-bottom: 1px solid #d0d5dc;
      border-radius: 0;
      background: transparent;
      font-weight: 600;
      text-align: left;
    }
    .sl-chevron {
      flex: none;
      width: 14px;
      height: 14px;
    }
    .sl-header[aria-expanded='true'] .sl-chevron {
      transform: rotate(90deg);
    }
    .sl-cards {
      display: grid;
      grid-template-columns: repeat(auto-fill, 120px);
      gap: 8px;
      padding: 8px 2px 12px;
    }
    .dialog button.sl-card {
      display: flex;
      flex-direction: column;
      align-items: center;
      gap: 4px;
      min-width: 0;
      padding: 6px;
      border: 1px solid #c3c8d0;
      border-radius: 6px;
      background: #fff;
      color: #1f2328;
      font-size: 12px;
      line-height: 1.3;
    }
    @media (hover: hover) {
      .dialog button.sl-card:hover {
        border-color: #1a5fb4;
        box-shadow: inset 0 0 0 1px #1a5fb4;
      }
    }
    .sl-thumb {
      display: flex;
      align-items: center;
      justify-content: center;
      width: 100%;
    }
    /* The drawing is cut to its ink by its viewBox, so it does not clip (an <svg>'s default): what shows
       of a card is the card itself. */
    .sl-drawing {
      display: block;
      flex: none;
      overflow: visible;
      background: #fff;
      pointer-events: none;
    }
    .sl-caption {
      display: block;
      width: 100%;
      text-align: center;
      overflow-wrap: anywhere;
    }
    /* In forced colours (Windows High Contrast) the browser drops our fills, shadows and colours, so what
       they show is drawn here in system colours, last, after the rules it must win against (fix round
       4; .btn:focus once hid a focused pressed button's outline):
       - a button disabled with aria-disabled (it keeps its Tab stop, reading 1) greys itself, since the
         browser greys only the native disabled attribute;
       - the element symbols of the palette and of the periodic table's cells take the system's button text,
         not the elements' colours (ATOM-033, P4). The browser paints them so anyway; declaring it makes
         the colour a script reads off a symbol the one it shows, since Chromium reports the author's colour
         as the symbol's text fill whatever it paints;
       - a pressed button has a border of Highlight, and the periodic table's chosen element an outline, so
         that a pressed button that keyboard focus is on shows its focus ring too: the border is the
         pressed state, the outline the focus (the ui-critic's note 10);
       - a tooltip keeps its box, drawn as a border in the text's colour (note 11);
       - an invalid field's border turns dashed, beside its message. */
    /* Reduced motion (spike settings-touch, K10, A11Y-019): nothing moves by itself, and a bar scrolls without smoothing:
       the sideways bars' smooth scrolling is turned off here (milestone-4-fixes, F9). */
    @media (prefers-reduced-motion: reduce) {
      *,
      *::before,
      *::after {
        transition: none !important;
        animation: none !important;
        scroll-behavior: auto !important;
      }
    }
    @media (forced-colors: active) {
      .bar .scroll {
        border: 1px solid ButtonText;
      }
      .btn[aria-disabled='true'],
      .menu .item[aria-disabled='true'] {
        color: GrayText;
      }
      .symbol,
      .dialog button.pt-el {
        color: ButtonText;
      }
      .btn[aria-pressed='true'] {
        border: 2px solid Highlight;
      }
      .dialog button.pt-el[aria-pressed='true'],
      .dialog button.r-toggle[aria-pressed='true'] {
        outline: 2px solid Highlight;
        outline-offset: -2px;
      }
      .tip {
        border-color: CanvasText;
      }
      .fields input[aria-invalid='true'],
      .fields select[aria-invalid='true'],
      .label-editor[aria-invalid='true'] {
        border-style: dashed;
        border-width: 2px;
      }
      /* The drawing area's marks (spike milestone-6-fixes, F1, R1): forced colours leave the drawing white (the engine's
         background is a fill, which they do not touch), so its marks are drawn for white, 3:1 or more against it, in every
         forced-colours theme, and the browser adjusts none of them. */
      .overlay [data-testid='overlay.selection'],
      .overlay [data-testid='overlay.hover'],
      .overlay [data-testid='overlay.problems'] {
        forced-color-adjust: none;
      }
      /* A light theme's system colours pair with that white, so there the marks take them: */
      @media (prefers-color-scheme: light) {
        /* The selection's marks (spike select, S5) in the system's colour for what is selected, which forced
           colours keep apart from the text and the canvas. */
        .overlay [data-testid='overlay.selection'] > * {
          fill: Highlight;
          fill-opacity: 0.3;
          stroke: Highlight;
        }
        /* A selected bond's band is its fill alone here too (spike polish-look, row 1), so that the bond's order reads. */
        .overlay [data-testid='overlay.selection'] > [data-band] {
          stroke: none;
        }
        /* The hotspot's mark on what is selected (spike v1-display, CANVAS-018): the system's Highlight, deeper than the
           selection's, so that it is told apart from a plain hotspot's; its outline as a plain hotspot's is, not the
           selection's Highlight. A band keeps its fill alone. */
        .overlay [data-testid='overlay.hover'][data-selected='true'] > * {
          fill: Highlight;
          fill-opacity: 0.6;
        }
        /* A flagged atom's ring (spike charges-checks, K10) in the system's text colour, dashed as ever, so that it is told
           apart from the selection's Highlight by its colour and its dashes. */
        .overlay [data-testid='overlay.problems'] > * {
          fill: none;
          stroke: CanvasText;
        }
      }
      /* A dark theme's are made for its dark canvas (white text, a pale Highlight: white on white, a faint cyan, the
         ui-critic's milestone-6 blocking finding), so there the marks keep the colours made for the white drawing, as
         outside forced colours: the problem ring's dashed red, the selection's green, the hotspot's blue. */
      .info {
        border-color: CanvasText;
      }
      /* The stereochemistry panel (spike stereo-panel; the design's forced colours): the pressed segment Highlight on
         HighlightText, the selected row outlined in Highlight, borders ButtonText, disabled GrayText, the problem mark
         CanvasText. Nothing is told by colour alone: ? pressed and the summary say "unspecified"; aria-selected. */
      .sp {
        border-color: CanvasText;
      }
      .sp .seg,
      .sp .grp {
        border-color: ButtonText;
      }
      .sp .seg button[aria-pressed='true'],
      .sp .seg button[aria-pressed='true'].unset {
        forced-color-adjust: none;
        background: Highlight;
        color: HighlightText;
        box-shadow: none;
      }
      .sp .grp[aria-disabled='true'],
      .sp .sp-cmd[aria-disabled='true'] {
        border-color: GrayText;
        color: GrayText;
      }
      .sp .pm {
        stroke: CanvasText;
      }
      .sp-row {
        border-left-color: Canvas;
      }
      .sp-row[aria-selected='true'] {
        border-left-color: Highlight;
        outline: 2px solid Highlight;
        outline-offset: -2px;
      }
      .sp-menu {
        border-color: CanvasText;
      }
    }
  `;static get observedAttributes(){return[...super.observedAttributes,`value`,`mode`,`atom-colours`,`keymap`,`warning-marks`,`analysis-panel`,`stereo-labels`,`stereo-label-style`,`stereo-flag`,`unspecified-double-bonds`,`locale`,`methyl-labels`,`cxsmiles-coordinates`,`density`,`persist-settings`,`readout-items`,`stereo-panel`]}ready;#e;#t=new Em(co);#n=null;#r=!1;#i;#a;#o=null;#s=!1;#c;#l=null;#u=new mv((e,t)=>this.#t.expire(e,t),()=>this.#Ot());#d=null;#f=null;#p=null;#m=null;#h=null;#g;#_;#v;#y=Zi(Ui,tn);#b=Object.freeze([]);#x=null;#S;#C=!1;#w={undo:!1,redo:!1,empty:!0};#T=!1;#E=!1;#D=null;#O=`{"atoms":[],"bonds":[]}`;#k=[];#A;#j;#M;#N;#P=null;#F;#I=!1;#L=null;#R=null;#z=`mouse`;#B=null;#V;#H=null;#U=!1;#W;#G;#K;#q;#J;#Y;#X;#Z=null;#Q=null;#$;#ee;#te;#ne;#re=null;#ie=`closed`;#ae=null;#oe=!1;#se=null;#ce=Cd;#le;#ue;#de;#fe;#pe;#me;#he;#ge;#_e=Object.freeze({});#ve=new Set;#ye=Qo.chemdraw;#be=`hidden`;#xe=`hidden`;#Se=Vh;#Ce=!0;#we=new Set;#Te=!1;#Ee=!1;#De=null;#Oe=null;#ke=null;#Ae=!1;#je=null;#Me={kind:`canvas`};#Ne=null;#Pe=!1;#Fe=new Set;#Ie=`chemdraw`;#Le=`en`;#Re=te(void 0);#ze=p;#Be=fg;#Ve=`element`;#He=null;#Ue=null;#We;#Ge=new Map;#Ke=new Set;#qe=null;#Je={x:0,y:0};#Ye=!1;#Xe=null;#Ze=1;#Qe=null;#$e=`shown`;#et=`shown`;#tt=$m;#nt=null;#rt=`null`;#it=!1;#at=null;#ot=!1;#st=`100 true true`;#ct=Gt;#lt=!0;#ut=!1;#dt=null;#ft=null;#pt(){if(this.#C)return;let e=this.getAttribute(`role`);(e===null||e===this.#dt)&&this.setAttribute(`role`,this.#dt=`group`);let t=this.getAttribute(`aria-label`);(t===null||t===this.#ft)&&this.setAttribute(`aria-label`,this.#ft=this.#ze.sketcherLabel)}#mt=()=>{this.#lt=!1,this.setAttribute(`aria-busy`,`false`)};constructor(){super(),this.#S=Uy(this,[`value`,`mode`,`disabledQueryElements`,`toolbar`,`atomColours`,`keymap`,`warningMarks`,`analysisPanel`,`readoutItems`,`stereoPanel`,`stereoLabels`,`stereoLabelStyle`,`stereoFlag`,`unspecifiedDoubleBonds`,`locale`,`strings`,`methylLabels`,`keys`,`settings`,`cxsmilesCoordinates`,`density`,`persistSettings`,`viewOnly`]);let e=this.ownerDocument??document,t=e.createElement(`div`);t.className=`stage`;let r=e.createElement(`div`);r.className=`surface`;let i=new gl(e);this.#o=i;let a=e.createElement(`div`);a.className=`status`,a.setAttribute(`role`,`status`),a.setAttribute(`aria-live`,`polite`),a.setAttribute(`data-testid`,`canvas.status`);let o=e.createElement(`button`);o.type=`button`,o.className=`dismiss`,o.tabIndex=-1,o.hidden=!0,o.setAttribute(`data-testid`,`canvas.status.close`),o.textContent=`×`,o.addEventListener(`mousedown`,Cy),o.addEventListener(`click`,this.#jt),this.#t.attach(()=>{this.#l?.removeEventListener(`mousedown`,Cy),this.#l?.removeEventListener(`click`,this.#jt)});let s=e.createElement(`div`);s.className=`status-line`,s.append(a,o);let c=e.createElement(`div`);c.className=`sr-only`,c.setAttribute(`role`,`status`),c.setAttribute(`aria-live`,`polite`),c.setAttribute(`aria-atomic`,`true`),c.setAttribute(`data-testid`,`canvas.announce`);let l=e.createElement(`div`);l.className=`sr-only`,l.id=vy,l.setAttribute(`data-testid`,`canvas.description`),this.#d=c,this.#f=l,this.#l=o;let u=new jv(e,{strings:()=>this.#ze,act:e=>{e===`delete`?this.#_?.deleteSelection():this.#_?.flip(e===`flip.horizontal`?`horizontal`:`vertical`,`mirror`)},keys:e=>ds(this.#yr,e===`delete`?`selection.delete`:e,Hy(this.ownerDocument?.defaultView),this.#ze.keysOr),after:(e,t)=>this.#t.after(e,t)});this.#j=u;let d=new lv(e,{toggled:e=>this.#Nr({...this.settings,analysisPanel:e}),away:()=>this.#g.svg?.focus({preventScroll:!0}),copy:e=>this.#St(e),menu:(e,t,n)=>this.#wt(e,n),more:(e,t,n)=>this.#Et(e,t,n),after:(e,t)=>this.#t.after(e,t)});this.#Qe=d,d.state=`shown`,this.#t.attach(()=>this.#Qe?.destroy());let f=e.createElement(`div`);f.className=`foot`,f.append(s,d.element);let p=new Q_(e,{strings:()=>this.#ze,read:()=>this.#Ur(),selected:()=>this.#_?.selected??{atoms:[],bonds:[]},select:e=>this.#_?.selectItems(e),mark:e=>this.#_?.markItems(e),config:(e,t)=>void this.#_?.panelConfig(e,t),group:(e,t,n)=>void this.#_?.panelGroup(e,t,n),command:e=>void this.#_?.panelStereoCommand(e),close:()=>this.#Lr(!1,`canvas`),leave:()=>this.#g.svg?.focus({preventScroll:!0}),history:e=>this.#Wr(e),reveal:e=>this.#Br(e),later:e=>this.#t.painted(e),after:(e,t)=>this.#t.after(e,t),frame:e=>this.#t.painted(e,!0)});this.#re=p,this.#t.attach(()=>this.#re?.destroy()),t.append(r,i.element,i.chrome,u.element,f,p.element,c,l),this.#i=t,this.#a=r,this.#c=a;let m=$(e);this.#N=m;let h=$(e);this.#F=h;let g=$(e);this.#G=g;let _=$(e);this.#q=_,this.#J=null;let v=$(e);this.#X=v;let y=$(e);this.#Q=y;let b=$(e);this.#ee=b;let ee=$(e);this.#ne=ee;let te=$(e);this.#fe=te;let ne=$(e);this.#me=ne;let ie=$(e);this.#ge=ie;let ae=$(e);this.#ue=ae;let x=$(e);x.setAttribute(`data-menus`,``),this.#V=x,this.#Fr(),this.#g=new Dr(r,this.#gt,this.#ze.canvasLabel,()=>{this.#Vt(),this.#_?.afterRender(),this.#zt(),this.#_t(),this.#Kr(),this.#Hr(),this.#e.rendered()}),this.#g.describedBy=vy,this.#e=new t_({canvas:this.#g,cdxml:{loaded:Ri,load:Li},inchi:{loaded:Bi,load:zi},interrupt:()=>this.#_?.historyMoved(),hostCall:()=>this.#Gr()}),this.#_=new Sm(this.#e,this.#g,i,{toolChanged:()=>{this.#v.toolChanged(),this.#Ut()},showLabelEditor:(e,t,n)=>this.#A?.show(e,t,n),hideLabelEditor:()=>this.#A?.hide(),labelText:()=>this.#A?.value??``,openAtomProperties:e=>this.#M?.openAtomProperties(e,this.#Qn()),focusCanvas:()=>this.#g.svg?.focus({preventScroll:!0}),refuseLabel:e=>this.#A?.refuse(e),moveLabelEditor:(e,t)=>this.#A?.move(e,t),zoom:e=>{this.#Lt(e),this.#Mt(this.#ze.format(`zoomTo`,{percent:this.#ze.percent(xt(this.#e.view?.zoom??1))}),_y)},announce:(e,t)=>this.#Mt(e,t===`hotspot`?hy:_y),notice:(e,t)=>this.#kt(e,t),relayout:e=>this.#Rt(e),openRGroup:e=>this.#W?.openRGroup(e,this.#Qn()),openSGroup:e=>this.#K?.openFor(e,this.#Qn()),sgroupTip:(e,t)=>this.#er(e,t),openAttachmentPoints:e=>this.#W?.openAttachmentPoints(e,this.#Qn()),openAddChain:e=>this.#Z?.show({id:`dialog.add-chain`,title:this.#ze.addChain,label:this.#ze.addChainLength,value:e.last??3,min:1,max:50,integer:!0,invalid:this.#ze.addChainInvalid,ok:t=>void this.#_?.addChain(e,t)},this.#Qn()),refusal:(e,t)=>this.#Dt(e,t),pan:(e,t)=>{this.#C||this.#We.panBy(e,t)},selectionActions:e=>this.#j?.show(e),selectionChanged:e=>{this.#C||e===this.#T||(this.#T=e,this.requestUpdate())},selectionUpdated:()=>{this.#hr()&&this.requestUpdate(),this.#_t(),this.#Ft(),this.#Kr(),this.#re?.selectionChanged()},hotChanged:e=>this.#re?.hot(e),noticed:(e,t)=>{this.#C||this.#e.noticed(Object.freeze({kind:t,message:e}))},copyAs:e=>this.#wn(e),pasteAs:e=>void this.#Dn(e,this.#Qn()),stereoAvailable:e=>{this.#C||e===this.#Ae||(this.#Ae=e,this.requestUpdate())},openMenu:()=>this.#Pn(`keyboard`),openBondProperties:e=>this.#$?.show(e,this.#Qn()),openKeys:()=>this.#pe?.show(this.#Qn()),openSettings:()=>this.#de?.show(this.#Qn()),openLibrary:()=>this.#nr(this.#Qn()),confirmDissolve:e=>this.#le?.show(e,this.#Qn()),openLookup:e=>this.#or(e)},void 0,this.#t),this.#_.queryOff=new Set(this.#b),this.#We=new ug({view:()=>this.#e.view,size:()=>this.#g.size,show:e=>void this.#e.showView(e),move:e=>this.#Bt(e)},this.#t),this.#A=new xv(t,{commit:()=>this.#_?.commitLabel(!0),cancel:()=>this.#_?.cancelLabel(),blur:e=>this.#_?.commitLabel(this.#$n(e),!0),press:e=>this.#_?.labelPress(e)??!1},this.#ze.labelEditor),this.#M=new Vn({strings:()=>this.#ze,atomInfo:e=>{let t=this.#_?.atomIndex(e)??-1;return t>=0?me(this.#e.atomJson(t)):null},applyAtom:(e,t)=>{let n=this.#_?.atomIndex(e)??-1;n>=0&&this.#_?.apply(S(`atom.properties`,t(n)))},applyAtoms:(e,t)=>{let n=e.map(e=>this.#_?.atomIndex(e)??-1);n.length&&n.every(e=>e>=0)&&this.#_?.apply(S(`atom.properties`,t(n)))},refuseAtoms:(e,t)=>{let n=e.map(e=>this.#_?.atomIndex(e)??-1);if(!n.length||n.some(e=>e<0))return null;let r=t(n),i=r.flatMap(e=>e.op===`setAtomRadical`&&e.radical>0&&typeof e.atom==`number`?[e.atom]:[]);return i.length?this.#_?.radicalRefusalOf(S(`atom.properties`,r),i)??null:null},activeElement:()=>this.#_?.element??null,atomColours:()=>this.#Ve,confirmElement:e=>this.#pr(e),restoreFocus:e=>this.#tr(e),queryMode:()=>this.#e.mode===`query`,queryOff:e=>this.#b.includes(e),confirmQueryAtom:e=>this.#_?.setQueryAtom(e),tryQuery:(e,t)=>{let n=e.map(e=>this.#_?.atomIndex(e)??-1);return!n.length||n.some(e=>e<0)||!this.#_?{ok:!0,smarts:null}:this.#_.tryQuery(S(`atom.properties`,t(n)),n[0])}},m),this.#W=new _v({strings:()=>this.#ze,rNumbers:e=>{let t=this.#_?.atomIndex(e)??-1;return t>=0?this.#e.atomRNumbers(t):null},attachmentPoints:e=>{let t=this.#_?.atomIndex(e)??-1;return t>=0?this.#e.atomAttachmentPoints(t):null},applyRNumbers:(e,t,n)=>{if(!(`atomId`in e)){this.#_?.apply(Mp(e.at,n));return}let r=this.#_?.atomIndex(e.atomId)??-1;r>=0&&this.#_?.apply(jp(r,t,n))},applyAttachmentPoints:(e,t,n)=>{let r=this.#_?.atomIndex(e)??-1;r>=0&&this.#_?.apply(Rp(r,t,n))},restoreFocus:e=>this.#tr(e)},g),this.#Y=new Cv({strings:()=>this.#ze,canRead:(e,t)=>this.#_?.canPasteAs(e,t)??!1,paste:(e,t)=>void this.#_?.pasteAs(e,t),restoreFocus:e=>this.#tr(e),pasteKey:()=>this.#Sn(`paste`),queries:()=>this.#e.mode===`query`},v),this.#Z=new Tv({strings:()=>this.#ze,restoreFocus:e=>this.#tr(e)},y),this.#$=new n_({strings:()=>this.#ze,bondChoice:e=>this.#_?.bondChoiceOf(e)??null,apply:(e,t)=>void this.#_?.applyBondChoice(e,t),restoreFocus:e=>this.#tr(e),bondForm:e=>this.#_?.bondFormOf(e)??null,applyForm:(e,t,n)=>this.#_?.applyBondForm(e,t,n)??{ok:!0},bondSmarts:(e,t,n)=>this.#_?.bondSmarts(e,t,n)??``,queryOff:e=>this.#b.includes(e)},b),this.#K=new vv({strings:()=>this.#ze,opening:e=>this.#_?.sgroupOpening(e)??null,apply:(e,t)=>this.#_?.applySGroup(e,t)??null,restoreFocus:e=>this.#tr(e)},_),this.#le=new i_({strings:()=>this.#ze,answer:e=>void this.#_?.dissolveAnswered(e),restoreFocus:e=>this.#tr(e)},ae),this.#de=new Yh({strings:()=>this.#ze,settings:()=>this.settings,apply:e=>this.#Nr(e),restoreFocus:e=>this.#tr(e)},te),this.#pe=new Zh({strings:()=>this.#ze,sheet:()=>Mm(this.#yr,this.#ze,Hy(this.ownerDocument?.defaultView)),restoreFocus:e=>this.#tr(e)},ne),this.#he=new $h({strings:()=>this.#ze,version:()=>this.version,engineRev:()=>this.#e.engineRev,restoreFocus:e=>this.#tr(e)},ie),this.#te=new av({closed:()=>this.#_?.menuClosed(),chose:e=>this.#Fn(e),frame:e=>this.#t.painted(e,!0)},ee),this.#t.attach(()=>this.#te?.destroy()),this.#v=new qv({strings:()=>this.#ze,requestUpdate:()=>this.requestUpdate(),root:()=>this.renderRoot,activeTool:()=>this.#_?.tool??`bond.single`,chooseTool:e=>this.#_?.setTool(e),runAction:e=>{e===`clear`?this.#_?.clear():e===`undo`?this.#e.undo(`user`):e===`redo`?this.#e.redo(`user`):e===`element.table`?this.#M?.openTable(`palette`,this.#Qn()):e===`zoom.in`?this.#Lt(`in`):e===`zoom.out`?this.#Lt(`out`):e===`zoom.reset`?this.#Lt(`reset`):e===`zoom.fit`?this.#Lt(`fit`):e===`layout`||e===`clean`?this.#Rt(e):e===`cut`||e===`copy`?this.#yn(e):e===`paste`?this.#bn():e===`structure.stereo-panel`||e===`structure.enhanced-stereo`?this.#Lr(!this.#re?.open,`panel`):e===`structure.aromatize`||e===`structure.dearomatize`?this.#_?.aromaticity(e===`structure.aromatize`?`aromatize`:`dearomatize`):e===`keys`?this.#pe?.show({kind:`control`,testid:`toolbar.${e}`}):e===`settings`&&this.#de?.show({kind:`control`,testid:`toolbar.${e}`})},open:e=>{e===`template.library`&&this.#nr({kind:`control`,testid:`toolbar.${Ar}`},this.#z===`touch`)},available:e=>this.#E&&Dy.has(e)?!1:e===`undo`?this.#e.canUndo:e===`redo`?this.#e.canRedo:e===`zoom.in`?pt(this.#Ze):e===`zoom.out`?re(this.#Ze):e===`layout`||e===`clean`?!this.#e.isEmpty:e===`cut`||e===`copy`?this.#T:e===`structure.aromatize`||e===`structure.dearomatize`?(this.#_?.aromaticityState(e===`structure.aromatize`?`aromatize`:`dearomatize`)??`none`)!==`none`:!0,pressed:e=>e===`structure.stereo-panel`||e===`structure.enhanced-stereo`?!!this.#re?.open:void 0,readout:()=>this.#ze.percent(xt(this.#Ze)),zoomTo:e=>{this.#C||this.#We.zoomTo(e,void 0,this.#_?.selectionCentre()??null)},hotkey:e=>{let t=Hy(this.ownerDocument?.defaultView),n=mn(e),r=(e,t)=>this.#ze.format(t?`cycleFirst`:`cycleLater`,{key:e});if(n)return vs(this.#yr,n,t)??(Xr(e)?bs(this.#yr,e,t,this.#ze.keysOr,r):void 0);let i=Oy[e];return[i&&!i.startsWith(`tool.`)?ds(this.#yr,i,t,this.#ze.keysOr):void 0,Xr(e)?bs(this.#yr,e,t,this.#ze.keysOr,r):void 0].filter(Boolean).join(this.#ze.keysOr)||void 0},menuAvailable:e=>this.#E&&e!==`copy`&&e!==`zoom`&&e!==`options`?!1:e===`copy`||e===`clean`||e===`hydrogens`?!this.#e.isEmpty:!0,menuItemAvailable:(e,t)=>t===`zoom.in`?pt(this.#Ze):t===`zoom.out`?re(this.#Ze):t===`clean`||t===`layout`?!this.#e.isEmpty:t===`hydrogens.add`||t===`hydrogens.remove`?!!this.#_?.hydrogensAvailable()[t===`hydrogens.add`?`add`:`remove`]:!0,runMenu:(e,t)=>this.#Cn(e,t),hidden:e=>this.#ve.has(e),queryAtom:()=>{let e=this.#e.mode===`query`?this.#_?.queryAtom:null;return e?{face:e.kind===`generic`?n(e):`${e.kind===`not-list`?`!`:``}${(e.elements??[]).join(`,`)}`,name:this.#_r(e)}:null},toolOff:e=>this.#vr(e),queryMode:()=>this.#e.mode===`query`,menuItemChecked:(e,t)=>t===`mode`&&this.#e.mode===`query`,menuKeys:(e,t)=>{let n=Hy(this.ownerDocument?.defaultView),r=n===Br.cmd;if(e===`copy`||e===`paste`){let i=e===`copy`?`copy.as`:`paste.as`,a=ps(this.#yr,i,t,n,this.#ze.keysOr);return a?{label:a,aria:ms(this.#yr,i,t,r)??``}:void 0}let i=ky[t],a=i?ds(this.#yr,i,n,this.#ze.keysOr):void 0;return i&&a?{label:a,aria:hs(this.#yr,i,r)??``}:void 0},recents:()=>this.#Be.slots,atomColours:()=>this.#Ve,pressElement:(e,t)=>this.#Wn(e,t),wasDragged:()=>{let e=this.#Ue!==null;return this.#Ue=null,e},after:(e,t)=>this.#t.after(e,t)},this.#y),ny(this,this.#e),iy(this,(e,t)=>{e===`atomNumbers`&&(this.#Tr(`atomNumbers`),this.#xr(Oh({atomNumbers:t}).atomNumbers))}),this.#e.on(`change`,e=>{this.#ht(),this.#Pt(),this.dispatchEvent(new CustomEvent(`change`,{detail:e})),e.source===`api`&&this.#e.lastChange.kind===`write`&&this.#_?.reportRead(`write`)}),this.#e.on(`error`,e=>this.dispatchEvent(new CustomEvent(`error`,{detail:e}))),this.#e.on(`render`,e=>this.dispatchEvent(new CustomEvent(`render`,{detail:e}))),this.#e.on(`selection`,e=>this.dispatchEvent(new CustomEvent(`selection`,{detail:e}))),this.#e.on(`notice`,e=>this.dispatchEvent(new CustomEvent(`notice`,{detail:e}))),this.#in(r),this.#Ut(),this.ready=this.#Yt(),this.ready.then(this.#mt,this.#mt)}#ht(){let{canUndo:e,canRedo:t,isEmpty:n}=this.#e,r=this.#w;(e!==r.undo||t!==r.redo||n!==r.empty)&&(this.#w={undo:e,redo:t,empty:n},this.requestUpdate())}#gt=()=>this.#a?er(this.#a):{width:0,height:0};#_t(){this.#C||!this.#Qe||this.#et===`hidden`||this.#nt||(this.#nt=this.#t.painted(()=>{this.#nt=null,this.#vt()},!0))}#vt(){this.#nt?.(),this.#nt=null;let e=this.#Qe;if(this.#C||!e||this.#et===`hidden`)return;let t=this.#_?.selectedAtoms??[],n=t.length?t:void 0,r=this.#e,i=!1,a=e=>this.#ot||r.readoutKnown(e,n)?!1:(i=!0,!0),o=bh(vh(this.#tt,{mass:()=>n?r.massInfo(n):r.massInfo(),text:e=>e===`inchikey`&&a(`inchikey`)?null:r.readoutText(e,n),descriptors:()=>a(`descriptors`)?`pending`:r.descriptors(n),bonds:()=>r.bondsAmong(n),pointer:()=>this.#_?.pointerFacts()??null},t.length>0,r.isEmpty,this.#ze),e.line);if(e.show(o),e.yieldTo(this.#Xn),i){this.#yt();return}o.kind===`items`&&o.items.some(e=>e.pending)&&this.#bt()}#yt(){if(this.#at||this.#C)return;let e=this.ownerDocument?.defaultView,t=()=>{this.#at=null,this.#ot=!0;try{this.#vt()}finally{this.#ot=!1}};if(typeof e?.requestIdleCallback==`function`){let n=()=>{},r=e.requestIdleCallback(()=>{n(),t()},{timeout:by});n=this.#t.attach(()=>e.cancelIdleCallback(r)),this.#at=()=>n()}else this.#at=this.#t.after(0,t)}#bt(){if(this.#it||this.#C)return;this.#it=!0;let e=()=>{this.#it=!1,this.#_t()};this.#t.alive(this.#e.loadFormat(`inchi`)).then(e,t=>{t!==wm&&e()})}#xt(){if(this.#C||!(this.#tt.includes(`under-pointer`)||this.#tt.includes(`drawing`)))return;let e=JSON.stringify(this.#_?.pointerFacts()??null);e!==this.#rt&&(this.#rt=e,this.#_t())}#St(e){this.#C||!e.copy||this.#Ln(e.copy,this.#ze.clipboardBlocked(this.#ze.copy,this.#Sn(`copy`)),()=>this.#ze.format(`readoutCopied`,{item:e.name}))}#Ct(){let e=this.#ze,t=this.#Qe?.line,n=t?.kind===`items`?t.items.filter(e=>!e.pending&&e.copy):[],r=ph(e),i=new Set(this.#tt);return[{type:`item`,id:`copy-all`,label:e.readoutCopyAll,disabled:!n.length},...n.map(t=>({type:`item`,id:`copy.${t.id}`,label:e.format(`readoutCopyItem`,{item:t.name})})),{type:`separator`},...Zm.map(t=>({type:`item`,id:t.id,label:mh(t.id,e),submenu:t.items.map(e=>({type:`item`,id:`show.${e}`,label:r[e].name,checkbox:!0,checked:i.has(e),keepOpen:!0}))}))]}#wt(e,t){let n=this.#te;if(this.#C||!n)return;this.#_?.editingLabel&&this.#_.commitLabel(!1,!0);let r=this.getBoundingClientRect(),i=t.getBoundingClientRect();n.open({entries:this.#Ct(),label:this.#ze.readoutMenu,target:`readout`,testid:`menu.readout`,at:e??{x:i.left,y:i.top},box:{left:r.left,top:r.top,right:r.right,bottom:r.bottom},keyboard:e===null,back:t,chose:e=>this.#Tt(e),closed:()=>{},toggle:e=>(this.#Tt(e),this.#Ct())})}#Tt(e){if(this.#C)return;if(e===`copy-all`){let e=this.#Qe?.text??``;e&&this.#Ln(e,this.#ze.clipboardBlocked(this.#ze.copy,this.#Sn(`copy`)),()=>this.#ze.readoutCopiedAll);return}let t=this.#Qe?.line;if(e.startsWith(`copy.`)){let n=t?.kind===`items`?t.items.find(t=>t.id===e.slice(5)):void 0;n&&this.#St(n);return}let n=e.startsWith(`show.`)?e.slice(5):``;th(n)&&(this.#Nr({...this.settings,readoutItems:ih(this.#tt,n,!this.#tt.includes(n))}),this.#vt())}#Et(e,t,n){let r=this.#te,i=this.#Qe;if(this.#C||!r||!i||!e.length)return;let a=this.getBoundingClientRect(),o=t.getBoundingClientRect();i.moreOpen=!0,r.open({entries:e.map(e=>({type:`item`,id:e.id,label:xy(xh(e))})),label:this.#ze.readoutMore,target:`readout`,testid:`menu.readout-more`,at:{x:o.left,y:o.top},box:{left:a.left,top:a.top,right:a.right,bottom:a.bottom},keyboard:n,back:t,chose:t=>{let n=e.find(e=>e.id===t);n&&this.#St(n)},closed:()=>{this.#Qe&&(this.#Qe.moreOpen=!1)}})}#Dt(e,t){this.#C||(this.#u.refusal(e),e!==null&&t?this.#a?.setAttribute(`data-refused`,``):this.#a?.removeAttribute(`data-refused`))}#Ot(){let e=this.#u.text;this.#c&&this.#c.textContent!==e&&(this.#c.textContent=e),this.#l&&(this.#l.hidden=!this.#u.closable,this.#l.setAttribute(`aria-label`,this.#ze.dismiss))}#kt(e,t=!0){this.#C||(this.#c&&this.#u.text===e&&(this.#c.textContent=``),this.#u.notice(e,t?Qp:pv))}#At(){return this.#u.dismiss()}#jt=()=>{this.#At()};#Mt(e,t){if(this.#C||!e||this.#p&&this.#p.rank>t)return;let n=this.#p===null;this.#p={message:e,rank:t},n&&queueMicrotask(()=>this.#Nt())}#Nt(){let e=this.#p;this.#p=null;let t=this.#d;if(!(!e||!t||this.#C)){if(t.textContent!==e.message){t.textContent=e.message;return}t.textContent=``,this.#t.frame(()=>{this.#d&&(this.#d.textContent=e.message)})}}#Pt(){if(this.#C)return;let e=Ty(this.#g.svg,this.#ze);this.#Mt(Hs(this.#ze,this.#e.lastChange,this.#m,e),my),this.#m=e,this.#At(),this.#h||=this.#t.painted(()=>{this.#h=null,this.#It()})}#Ft(){let e=this.#_?.selected;!e||this.#C||this.#Mt(Ls(this.#ze,e.atoms.length,e.bonds.length),gy)}#It(){let e=this.#f;if(!e||this.#C)return;let t=this.#m??Ty(this.#g.svg,this.#ze),n=``,r=null;if(t.atoms.length){try{n=this.#e.smiles}catch{n=``}r=this.#e.analysis?.formula??null}let i=Gs(this.#ze,{atoms:t.atoms.length,bonds:t.bonds.length,formula:r,smiles:n});e.textContent!==i&&(e.textContent=i)}#Lt(e){if(this.#C)return;if(e===`fit`){this.#We.flush(),this.#e.fit();return}let t=this.#_?.selectionCentre()??null;e===`reset`?this.#We.zoomTo(1,void 0,t):this.#We.step(e===`in`?1:-1,void 0,t)}#Rt(e){if(this.#C||!this.#_)return;let t=e===`layout`||e===`clean`;t&&this.#We.flush();let n=this.#e.lastError,r=this.#_.relayout(e);t&&!r&&!this.#e.isEmpty&&this.#e.lastError===n&&this.#e.fit()}#zt(){let e=this.#e.view?.zoom??1;this.#Ze=e;let t=`${xt(e)} ${pt(e)} ${re(e)}`;t!==this.#st&&(this.#st=t,this.requestUpdate())}#Bt(e){if(this.#C||!this.#e.moveView(e))return!1;let t=this.#g.move(e),n=this.#g.size;return!t||!n?!1:(this.#s||(this.#s=!0,this.#_?.viewMoving(),this.#i?.setAttribute(`data-view-moving`,``)),this.#o?.moveWith(ze(t,n)),this.#zt(),!0)}#Vt(){this.#We?.drawn(),this.#s&&(this.#s=!1,this.#o?.moveWith(null),this.#i?.removeAttribute(`data-view-moving`))}#Ht(){!this.#C&&this.#We.unsettled&&!this.#We.busy&&this.#We.flush()}#Ut(){let e=this.#_?.tool;e===`hand`?this.#a?.setAttribute(`data-tool`,`hand`):e&&Fr(e)?this.#a?.setAttribute(`data-tool`,`select`):this.#a?.removeAttribute(`data-tool`)}#Wt(e,t){this.#_?.navigating(),this.#qe=e.pointerId,this.#Je=t,this.#Ye=!1,this.#We.grab(t),this.#a?.setAttribute(`data-panning`,``),By(e)}#Gt(e){let t=performance.now(),n=this.#ke;if(n&&t-n.time<=300&&Math.hypot(e.x-n.at.x,e.y-n.at.y)<=10){this.#ke=null,this.#Lt(`fit`);return}this.#ke={at:e,time:t}}#Kt(){this.#qe=null,this.#We.release(),this.#a?.removeAttribute(`data-panning`)}#qt(e){this.#Bn();let[t,n]=[...this.#Ge.entries()];t&&n&&(this.#_?.navigating(),this.#Ke.add(t[0]).add(n[0]),this.#We.pinchStart(t,n),By(e))}#Jt(e,t,n){if(this.#Ge.delete(e.pointerId),e.pointerId===this.#Xe&&(this.#Xe=null,this.#a?.removeAttribute(`data-moving`)),e.pointerId===this.#qe){this.#qe=null,this.#We.release(),this.#a?.removeAttribute(`data-panning`),!n&&e.pointerType===`mouse`&&this.#_?.pointerMove(t,e.pointerId,Xp(e)),!n&&e.pointerType===`touch`&&!this.#Ye&&this.#Gt(t);return}if(this.#Ke.delete(e.pointerId)){this.#We.busy&&this.#We.pinchEnd();return}this.#Ht(),n?this.#_?.pointerCancel():e.isPrimary&&this.#_?.pointerUp(t,e.pointerId,Xp(e))}async#Yt(){let e;try{e=await this.#t.alive(Ni())}catch(e){if(e===wm)return;throw this.#C||this.#e.failToLoad(e),e}if(this.#rn(),this.#C)return;this.#Dr();try{await this.#t.alive(this.#e.prepare(e))}catch(e){if(e===wm)return}if(this.#C)return;if(this.isConnected){if(await this.updateComplete,this.#C)return;this.#Sr(),await this.updateComplete}if(this.#C)return;this.isConnected&&this.#Qe?.fit(this.#gt()),this.#e.start(e),this.#m=Ty(this.#g.svg,this.#ze),this.#It(),this.#ht(),this.#e.releaseWith(this);let t=new ResizeObserver(e=>{let t=e.filter(e=>e.target===this).at(-1)?.contentRect,n=!1;t&&!this.#C&&(n=this.#Cr(t),n||this.#tn());let r=this.#ae;n||this.#zr();let i=this.#ae!==r&&this.#ae?.place!==`cover`,a=e.filter(e=>e.target!==this),o=a[a.length-1]?.contentRect;o&&!n&&(this.#te?.close(),this.#re?.dismiss(),this.#C||this.#e.resize(o),this.#Qe?.fit(o),this.#oe&&!i&&(this.#oe=!1,this.#Br())),this.#nn()});this.#x=t,this.#t.attach(()=>t.disconnect()),this.#a&&t.observe(this.#a),t.observe(this),this.isConnected&&await this.updateComplete}connectedCallback(){super.connectedCallback(),this.#lt&&this.setAttribute(`aria-busy`,`true`),this.#ut=!0,this.#pt(),this.#rn(),this.#C||this.requestUpdate()}disconnectedCallback(){super.disconnectedCallback(),this.#M?.cancelAll(),this.#dr(),this.#cr(),this.#W?.cancelAll(),this.#K?.cancelAll(),this.#er(null,null),this.#Y?.cancel(),this.#Z?.cancel(),this.#$?.cancel(),this.#de?.cancel(),this.#pe?.cancel(),this.#he?.cancel(),this.#le?.cancel(),this.#te?.close(),this.#Bn(),this.#_?.cancelLabel(),this.#Yn(),this.#_?.forgetRefusal(),this.#At()}attributeChangedCallback(e,t,n){e===`value`?this.#C||(this.#e.value=n??``):e===`mode`?this.#C||(this.mode=Yg(n)):e===`atom-colours`?this.#C||(this.atomColours=n??`element`):e===`keymap`?this.#C||(this.keymap=$o(n)):e===`warning-marks`?this.#C||(this.warningMarks=n??`shown`):e===`analysis-panel`?this.#C||(this.analysisPanel=n??`shown`):e===`stereo-panel`?this.#C||(this.stereoPanel=Sd(n)):e===`readout-items`?this.#C||(this.readoutItems=nh(n??void 0)):e===`stereo-labels`?this.#C||(this.stereoLabels=n??`hidden`):e===`stereo-label-style`?this.#C||(this.stereoLabelStyle=n??`standard`):e===`stereo-flag`?this.#C||(this.stereoFlag=n??`shown`):e===`unspecified-double-bonds`?this.#C||(this.unspecifiedDoubleBonds=n??`crossed`):e===`locale`?this.#C||(this.locale=n??`en`):e===`methyl-labels`?this.#C||(this.methylLabels=n??`hidden`):e===`cxsmiles-coordinates`?this.#C||(this.#e.cxsmilesCoordinates=Oa(n??!1)):e===`density`?this.#C||(this.density=Uh(n)):e===`persist-settings`?this.#C||(this.persistSettings=Bh(n)):super.attributeChangedCallback(e,t,n)}shouldUpdate(){return!this.#C}createRenderRoot(){let e=super.createRenderRoot();if(e instanceof ShadowRoot&&(e.adoptedStyleSheets=[...e.adoptedStyleSheets,Lr(this.ownerDocument??document)]),e instanceof ShadowRoot&&!this.#n&&!this.#C){let t=new N_(e,this.#Xt);this.#n=t,this.#t.attach(()=>t.destroy());for(let e of wy)this.addEventListener(e,this.#Zt,{passive:!0});this.#t.attach(()=>{for(let e of wy)this.removeEventListener(e,this.#Zt);this.#Qt(!1)}),e.addEventListener(`keydown`,this.#$t,{capture:!0}),e.addEventListener(`keyup`,this.#en),this.#t.attach(()=>{e.removeEventListener(`keydown`,this.#$t,{capture:!0}),e.removeEventListener(`keyup`,this.#en),this.#Fe.clear()})}return e}#Xt=()=>this.#Qt(!0);#Zt=e=>{this.#r&&y_(e)&&this.#Qt(!1)};#Qt(e){e===this.#r||e&&this.#C||(this.#r=e,this.renderRoot?.querySelector(`.frame`)?.toggleAttribute(`data-host-pointer`,e))}#$t=e=>{this.#ne&&e.composedPath().includes(this.#ne)&&this.#Fe.add(e.code)};#en=e=>{let t=e.code;(this.#Fe.delete(t)||this.#ne&&e.composedPath().includes(this.#ne))&&e.stopPropagation()};updated(){this.#v.sync(),this.#nn(),this.#tn()}#tn(){this.#C||!this.renderRoot||this.#v.fit(this.renderRoot.querySelector(`.bar-top`))}#nn(){this.#C||!this.renderRoot||(Jv(this.renderRoot.querySelector(`.bar-left`),`y`),Jv(this.renderRoot.querySelector(`.bar-right`),`y`),Jv(this.renderRoot.querySelector(`.bar-top`),`x`),Jv(this.renderRoot.querySelector(`.bar-bottom`),`x`))}render(){return this.#i?x`<div
        class="frame"
        lang=${this.#ze.locale}
        data-density=${Wh(this.#Se,this.#Oe)}
        data-panel-density=${this.#wr()}
        @keydown=${this.#an}
        @keyup=${this.#on}
        @focusout=${this.#cn}
        @pointermove=${this.#Hn}
        @pointerleave=${this.#Un}
        @copy=${this.#un}
        @cut=${this.#dn}
        @paste=${this.#_n}
        @pointerdown=${this.#An}
        @dragenter=${this.#kn}
        @dragover=${this.#kn}
        @dragleave=${this.#jn}
        @drop=${this.#Mn}
      >
        ${this.#v.render(`top`)} ${this.#v.render(`left`)} ${this.#i} ${this.#v.render(`right`)} ${this.#v.render(`bottom`)}
      </div>
      ${this.#N??w}${this.#F??w}${this.#G??w}${this.#q??w}${this.#X??w}${this.#Q??w}${this.#ee??w}${this.#fe??w}${this.#me??w}${this.#ge??w}${this.#ue??w}${this.#V??w}${this.#ne??w}`:w}#rn(){let e=this.#S;e&&(this.#S=null,`toolbar`in e&&(this.toolbar=e.toolbar),`atomColours`in e&&(this.atomColours=e.atomColours),`keymap`in e&&(this.keymap=e.keymap),`warningMarks`in e&&(this.warningMarks=e.warningMarks),`analysisPanel`in e&&(this.analysisPanel=e.analysisPanel),`readoutItems`in e&&(this.readoutItems=e.readoutItems),`stereoPanel`in e&&(this.stereoPanel=e.stereoPanel),`stereoLabels`in e&&(this.stereoLabels=e.stereoLabels),`stereoLabelStyle`in e&&(this.stereoLabelStyle=e.stereoLabelStyle),`stereoFlag`in e&&(this.stereoFlag=e.stereoFlag),`unspecifiedDoubleBonds`in e&&(this.unspecifiedDoubleBonds=e.unspecifiedDoubleBonds),`locale`in e&&(this.locale=e.locale),`strings`in e&&(this.strings=e.strings),`methylLabels`in e&&(this.methylLabels=e.methylLabels),`keys`in e&&(this.keys=e.keys),`settings`in e&&(this.settings=e.settings),`cxsmilesCoordinates`in e&&(this.cxsmilesCoordinates=e.cxsmilesCoordinates),`density`in e&&(this.density=e.density),`persistSettings`in e&&(this.persistSettings=e.persistSettings),`viewOnly`in e&&(this.viewOnly=e.viewOnly),`value`in e&&(this.value=e.value),`mode`in e&&(this.mode=e.mode),`disabledQueryElements`in e&&(this.disabledQueryElements=e.disabledQueryElements))}#in(t){let n=e=>{let t=e.currentTarget;return nn({x:e.clientX,y:e.clientY},t.getBoundingClientRect(),kn(t))},r=(e,n)=>{let r=(e=>{q(e)||n(e)});t.addEventListener(e,r),this.#k.push([e,r])};r(`pointerdown`,e=>{if(this.#_&&(this.#_.pointerKind=e.pointerType),this.#v.pressed(),this.#j?.pressed(),this.#te?.isOpen){this.#te.close();return}this.#Ht();let t=n(e);if(e.pointerType===`touch`){if(this.#Ge.set(e.pointerId,t),By(e),this.#Ge.size===2&&this.#qe!==null&&this.#Ge.has(this.#qe)&&!this.#Ke.size){this.#Kt(),this.#qt(e);return}if(this.#Ke.size||this.#We.busy){this.#Ke.add(e.pointerId);return}if(this.#Ge.size===2){this.#qt(e);return}}if(Be(e,Hy(this.ownerDocument?.defaultView)===Br.cmd)||this.#We.busy)return;if(e.pointerType===`mouse`&&e.button===1){this.#Wt(e,t);return}if(!e.isPrimary||e.pointerType===`mouse`&&e.button!==0)return;let r=this.#_?.tool===`hand`;if(r&&!this.#_?.grabs(t,e.pointerType!==`mouse`,e.altKey)){this.#Wt(e,t),e.pointerType===`touch`&&this.#zn(e.pointerId,t);return}this.#_?.pointerDown(t,e.pointerId,e.pointerType===`touch`,Xp(e))&&(By(e),r&&(this.#Xe=e.pointerId,this.#a?.setAttribute(`data-moving`,``)),e.pointerType===`touch`&&this.#zn(e.pointerId,t))}),r(`mousedown`,e=>{e.button===0?this.#_?.clickCount(e.detail):e.button===1?e.preventDefault():e.button===2&&this.#_?.busy&&this.#_.cancelGesture()});let i=e=>{if(this.#C||q(e))return;e.preventDefault(),e.stopPropagation();let t=e;this.#_?.busy||this.#te?.isOpen||t.pointerType===`touch`||this.#Ge.size||(this.#Ht(),this.#Pn(`pointer`,n(t)))};t.addEventListener(`contextmenu`,i),this.#k.push([`contextmenu`,i]),r(`pointermove`,e=>{e.button===2&&this.#_?.busy&&this.#_.cancelGesture();let t=n(e);this.#Ge.has(e.pointerId)&&this.#Ge.set(e.pointerId,t);let r=this.#Ne;if(r&&r.pointerId===e.pointerId&&!In(r.at,t)&&this.#Bn(),!this.#te?.isOpen){if(this.#We.busy){e.pointerId===this.#qe?(this.#We.drag(t),Math.hypot(t.x-this.#Je.x,t.y-this.#Je.y)>=4&&(this.#Ye=!0)):this.#Ke.has(e.pointerId)&&this.#We.pinchMove(e.pointerId,t);return}this.#Ke.has(e.pointerId)||e.isPrimary&&this.#_&&(this.#Ht(),this.#_.busy||(this.#_.pointerKind=e.pointerType),this.#_.pointerMove(t,e.pointerId,Xp(e)))}}),r(`pointerup`,e=>{this.#Bn(),this.#Jt(e,n(e),!1)}),r(`pointercancel`,e=>{this.#Bn(),this.#Jt(e,n(e),!0)}),r(`lostpointercapture`,e=>{(this.#Ge.has(e.pointerId)||e.pointerId===this.#qe||e.pointerId===this.#Xe)&&this.#Jt(e,n(e),!0)}),r(`pointerleave`,()=>this.#_?.pointerLeave());for(let e of[`pointerdown`,`pointermove`,`pointerup`,`pointercancel`,`pointerleave`])r(e,()=>this.#xt());let a=e=>{let t=this.#g?.svg;this.#C||!t||e.target!==t||this.#g?.showing||t.matches(`:focus-visible`)&&this.#_?.keyboardFocus()};t.addEventListener(`focusin`,a),this.#k.push([`focusin`,a]),r(`dblclick`,e=>{this.#Ht(),!(this.#C||this.#Ye||this.#_?.tool!==`hand`||this.#_.busy||this.#_.grabs(n(e)))&&this.#Lt(`fit`)});let o=e=>{if(e.preventDefault(),this.#C||this.#_?.busy)return;this.#te?.close();let t=e.currentTarget;this.#We.wheel({deltaX:e.deltaX,deltaY:e.deltaY,deltaMode:e.deltaMode,ctrl:e.ctrlKey||e.metaKey,shift:e.shiftKey},nn({x:e.clientX,y:e.clientY},t.getBoundingClientRect(),kn(t)))};t.addEventListener(`wheel`,o,{passive:!1}),this.#k.push([`wheel`,o]),t.addEventListener(`touchmove`,Vy,{passive:!1}),this.#k.push([`touchmove`,Vy]);let s=e=>{this.#Pe&&!e.touches.length&&(this.#Pe=!1,e.cancelable&&e.preventDefault())};t.addEventListener(`touchend`,s,{passive:!1}),t.addEventListener(`touchcancel`,s,{passive:!1}),this.#k.push([`touchend`,s],[`touchcancel`,s]);let c=t=>{this.#C||!e.#On(t)||(t.preventDefault(),t.dataTransfer&&(t.dataTransfer.dropEffect=`copy`))},l=t=>{if(this.#C||!e.#On(t))return;t.preventDefault();let r=t.dataTransfer?.files[0];r&&!this.#E&&this.#Nn(r,n(t))},u=t=>{!this.#C&&e.#On(t)&&this.#n?.own(t)};t.addEventListener(`dragover`,c),t.addEventListener(`drop`,l),t.addEventListener(`dragenter`,u),t.addEventListener(`dragleave`,u),this.#k.push([`dragover`,c],[`drop`,l],[`dragenter`,u],[`dragleave`,u])}#an=e=>{if(this.#z=`key`,this.#C||e.defaultPrevented)return;if(this.#Ht(),this.#H){let t=fo(e);if(t!==null||e.key===`Escape`||e.key===`Backspace`){t===null?e.key===`Backspace`?this.#H.typed=this.#H.typed.slice(0,-1):this.#H=null:this.#H.typed+=t,e.preventDefault();return}}let t=e.composedPath()[0]??null;if(t&&this.#re?.element.contains(t)||t&&t===this.#A?.element||t instanceof HTMLInputElement&&t.classList.contains(`zoom-input`))return;let n=!!t?.closest?.(`.bar`);if(e.key===`Tab`&&this.#ln(),e.key===`Escape`&&this.#He){this.#He.cancelled=!0,this.#_?.cancelDrop(),e.preventDefault();return}if(e.key===`Escape`){if(this.#_?.keyDown(e)){e.preventDefault();return}let t=this.#v.escape()||(this.#j?.escape()??!1)||(this.#re?.escape()??!1)||this.#At(),n=t?!1:this.#_?.escape()??!1;(t||n)&&e.preventDefault();return}if(n&&t&&this.#v.keyDown(e,t)){e.preventDefault();return}(!n&&!t?.closest?.(`.info`)||e.key!==`Enter`&&e.key!==` `)&&(this.#_?.keyDown(e,!!t&&t===this.#g.svg)?e.preventDefault():this.#sn(e)&&this.#n?.own(e))};#on=e=>{this.#C||e.defaultPrevented||this.#_?.keyUp(e)&&e.preventDefault()};#sn(e){let t=rs(this.#yr,e,`generic`)?.action;return t===`paste`||(t===`cut`||t===`copy`)&&this.#T}#cn=e=>{if(this.#v.focusOut(e),this.#C||this.#g?.showing)return;let t=e.relatedTarget;t&&this.renderRoot.contains(t)||this.#_?.keyboardGone()};#ln(){this.#_?.cancelGesture(),this.#He?.dragging&&(this.#He.cancelled=!0,this.#_?.cancelDrop())}#un=e=>this.#fn(e,!1);#dn=e=>this.#fn(e,!0);#fn(e,t){if(this.#C||e.defaultPrevented||!e.clipboardData||zy(e)||this.#ve.has(t?`cut`:`copy`)||t&&this.#E)return;let n=this.#_,r=n?.copySelection()??null;n&&r!==null&&(e.clipboardData.setData(`text/plain`,r),e.preventDefault(),this.#gn(r,n.copyPicture()),this.#At(),t&&n.deleteSelection())}#pn(e,t,n){let r=this.ownerDocument?.defaultView,i=r?.ClipboardItem,a=this.#mn();if(!r||!i||!t||!a)return null;let o={"text/plain":new r.Blob([e],{type:`text/plain`})};return a.png&&(o[`image/png`]=n??lg(t,this.ownerDocument,ng(t)).then(e=>e??Promise.reject(Error(`no picture`)))),a.svg&&(o[`image/svg+xml`]=new r.Blob([t],{type:`image/svg+xml`})),new i(o)}#mn(){let e=this.ownerDocument?.defaultView?.ClipboardItem;if(!e)return null;let t=typeof e.supports==`function`,n=!t||e.supports(`image/png`),r=t&&e.supports(`image/svg+xml`);return n||r?{png:n,svg:r}:null}#hn=0;#gn(e,t){let n=this.ownerDocument?.defaultView?.navigator?.clipboard,r=this.#mn();if(!n?.write||!r||!t)return;let i=++this.#hn;(r.png?lg(t,this.ownerDocument,ng(t)):Promise.resolve(void 0)).then(r=>{if(r===null||i!==this.#hn||this.#C)return;let a=this.#pn(e,t,r??void 0);return a?n.write([a]):void 0}).catch(()=>void 0)}#_n=e=>{this.#C||e.defaultPrevented||!e.clipboardData||zy(e)||this.#ve.has(`paste`)||this.#E||(e.preventDefault(),this.#xn(e.clipboardData.getData(`text/plain`)))};async#vn(e){try{let t=this.ownerDocument?.defaultView?.navigator?.clipboard;if(!t?.writeText)throw Error(`no asynchronous Clipboard API`);await this.#t.alive(t.writeText(e))}catch{this.#C||this.#_?.notify(this.#ze.clipboardBlocked(this.#ze.copy,this.#Sn(`copy`)));return}this.#C||this.#At()}async#yn(e){let t=this.#_,n=t?.copySelection()??null;if(!t||n===null)return;let[r,i]=[this.#e.stateToken,JSON.stringify(t.selected)];try{let e=this.ownerDocument?.defaultView?.navigator?.clipboard;if(!e?.writeText)throw Error(`no asynchronous Clipboard API`);let r=null;try{r=typeof e.write==`function`?this.#pn(n,t.copyPicture()):null}catch{r=null}let i=r?await this.#t.alive(e.write([r]).then(()=>!0,()=>!1)):!1;if(this.#C)return;i||await this.#t.alive(e.writeText(n))}catch{this.#C||this.#_?.notify(this.#ze.clipboardBlocked(e===`cut`?this.#ze.cut:this.#ze.copy,this.#Sn(e)));return}this.#C||this.#At(),!(e!==`cut`||this.#C||this.#e.stateToken!==r||JSON.stringify(this.#_?.selected)!==i)&&this.#_?.deleteSelection()}async#bn(){let e;try{let t=this.ownerDocument?.defaultView?.navigator?.clipboard;if(!t?.readText)throw Error(`no asynchronous Clipboard API`);e=await this.#t.alive(t.readText())}catch{this.#C||this.#_?.notify(this.#ze.clipboardBlocked(this.#ze.paste,this.#Sn(`paste`)));return}this.#C||this.#xn(e)}#xn(e){if(e.trim()){this.#_?.paste(e);return}let t=ps(this.#yr,`copy.as`,`cdxml`,Hy(this.ownerDocument?.defaultView),this.#ze.keysOr)??ps(Qo.chemdraw,`copy.as`,`cdxml`,Hy(this.ownerDocument?.defaultView),this.#ze.keysOr)??``;this.#_?.notify(this.#ze.pasteNoText(t))}#Sn(e){return ds(this.#yr,e,Hy(this.ownerDocument?.defaultView),this.#ze.keysOr)??``}#Cn(e,t){e===`copy`?this.#wn(t):e===`paste`?this.#Dn(t,{kind:`control`,testid:`toolbar.paste.open`}):e===`clean`?this.#Rt(t===`layout`?`layout`:`clean`):e===`hydrogens`?this.#_?.hydrogens(t===`hydrogens.add`?`add`:`remove`):e===`zoom`?(this.#Lt(Ay[t]),this.#Mt(this.#ze.format(`zoomTo`,{percent:this.#ze.percent(xt(this.#e.view?.zoom??1))}),_y)):t===`settings`?this.#de?.show({kind:`control`,testid:`toolbar.options`}):t===`keys`?this.#pe?.show({kind:`control`,testid:`toolbar.options`}):t===`about`?this.#he?.show({kind:`control`,testid:`toolbar.options`}):t===`mode`&&(this.mode=this.#e.mode===`query`?`molecule`:`query`,this.#Mt(this.#e.mode===`query`?this.#ze.queryModeOn:this.#ze.queryModeOff,_y))}#wn(e){let t=null,n=this.#_?.copyAs(e,e=>t=e)??null;this.#C||this.#Ln(n,this.#ze.copyAsBlocked(e),()=>t?we(this.#ze,t):null,this.#Tn(e),this.#En(e))}#Tn(e){return e===`inchi`||e===`inchikey`?this.#ze.inchiNotLoaded:this.#ze.cdxmlNotLoaded}#En(e){return e===`inchi`||e===`inchikey`?this.#ze.inchiWriterLoading:this.#ze.cdxmlWriterLoading}async#Dn(e,t){if(this.#C||!r(e,this.#e.mode===`query`))return;let n;try{let e=this.ownerDocument?.defaultView?.navigator?.clipboard;if(!e?.readText)throw Error(`no asynchronous Clipboard API`);n=await this.#t.alive(e.readText())}catch{this.#C||this.#Y?.show(e,t);return}this.#C||this.#_?.pasteAs(n,e)}static#On(e){return Array.from(e.dataTransfer?.types??[]).includes(`Files`)}#kn=t=>{this.#C||t.defaultPrevented||!e.#On(t)||(t.preventDefault(),t.dataTransfer&&(t.dataTransfer.dropEffect=`copy`))};#An={handleEvent:e=>{this.#z=e.pointerType||`mouse`},capture:!0};#jn=t=>{!this.#C&&e.#On(t)&&this.#n?.own(t)};#Mn=t=>{this.#C||t.defaultPrevented||!e.#On(t)||(t.preventDefault(),this.#_?.notify(this.#ze.dropOnCanvas))};async#Nn(e,t){let n=``,r=!1;try{n=await this.#t.alive(e.slice(0,Ey).text())}catch{r=!0}if(this.#C)return;this.#g.svg?.focus({preventScroll:!0});let i=e.name||this.#ze.unnamedFile,a=r?{kind:`nothing`}:At(i,e.type,()=>n,e=>this.#e.formatOf(e),Ja,this.#e.mode===`query`,(e,t)=>this.#e.inspect(e,t,()=>!0)??!1);a.kind===`read`?this.#_?.pasteAt(a.text,a.format,t,this.#ze.dropNothing(i)):this.#_?.notify(a.kind===`query-only`?this.#ze.format(`dropQueryOnly`,{format:this.#ze.formatName(a.format)}):a.kind===`not-yet`?this.#ze.dropNotYet(a.format):a.kind===`reaction`?this.#ze.reactionsRefused:a.kind===`cdx`?this.#ze.cdxNotRead:a.kind===`inchikey`?this.#ze.inchiKeyRefused:this.#ze.dropNothing(i))}#Pn(e,t=null){let n=this.#_,r=this.#a;if(this.#C||!n||!r||!this.#te||this.#E)return;n.editingLabel&&n.commitLabel(!1,!0);let a=n.menuTarget(e===`pointer`?t:null);if(!a)return;let o=this.renderRoot.activeElement??this.#g.svg;this.#Me=this.#Qn()??{kind:`canvas`},this.#je=a.target,n.menuOpened(a.target);let s=this.ownerDocument?.defaultView,c=Hy(s),l=n.menuFacts(a.target),u={...l,stereo:{...l.stereo,panel:this.#ie===`hidden`?null:!!this.#re?.open}},d=zt(u,e=>_s(this.#yr,e.action,e.format,c,c===Br.cmd),this.#ze),f=r.getBoundingClientRect(),p=this.getBoundingClientRect();this.#te.open({entries:d,label:Pr(a.target,this.#ze,!!u.abbreviations?.alone),target:i(a.target),at:Wn(a.anchor,f,kn(r)),box:{left:p.left,top:p.top,right:p.right,bottom:p.bottom},keyboard:e===`keyboard`,back:o})}#Fn(e){let t=this.#_,n=this.#je;if(this.#C||!t||!n)return;let r=this.#Me;if(e.startsWith(`bond.`)){t.menuBondChoice(n,e.slice(5));return}if(e.startsWith(`position.`)){t.menuPosition(n,e.slice(9));return}if(e.startsWith(`ez.`)){let r=e.slice(3);t.menuEz(n,r===`unknown`?`unknown`:r.toUpperCase());return}if(e.startsWith(`align.`)){t.menuAlign(n,e.slice(6));return}if(e.startsWith(`distribute.`)){t.menuDistribute(n,e.slice(11));return}if(e.startsWith(`atrop.`)){t.menuAtropisomer(n,e===`atrop.m`?`M`:`P`);return}if(e.startsWith(`stereo.`)){let r=e.slice(7);t.menuStereoCentre(n,r===`r`?`R`:r===`s`?`S`:null);return}if(e.startsWith(`copy-as.`)){let r=Et.find(t=>Rt(t)===e.slice(8)),i=null;r&&this.#Ln(t.menuCopy(n,r,e=>i=e),this.#ze.copyAsBlocked(r),()=>i?we(this.#ze,i):null,this.#Tn(r),this.#En(r));return}if(e.startsWith(`paste-special.`)){let t=fe.find(t=>Rt(t)===e.slice(14));t&&this.#Dn(t,r);return}if(e.startsWith(`query.`)){let[r,i]=e.slice(6).split(`.`),a=De(r);if(r===`clear`)t.menuClearQuery(n);else if(a===`unsaturated`)t.menuQueryField(n,a,!0);else if(a&&i!==void 0){let e=Qn(a,i);e!==void 0&&t.menuQueryField(n,a,e)}return}if(e.startsWith(`bond-query.`)){let r=e.slice(11);pr(r)&&t.menuBondQuery(n,r);return}if(e.startsWith(`topology.`)){let r=de.find(t=>t===e.slice(9));r&&t.menuTopology(n,r);return}switch(e){case`edit-label`:t.menuEditLabel(n);return;case`atom-properties`:{let e=t.menuAtomIds(n);e.length&&this.#M?.openAtomProperties(e.length===1?e[0]:e,r);return}case`rgroup`:n.atomIds.length&&this.#W?.openRGroup({atomId:n.atomIds[0]},r);return;case`attachment-points`:n.atomIds.length&&this.#W?.openAttachmentPoints(n.atomIds[0],r);return;case`edit-sgroup`:t.menuSGroup(n,`edit`);return;case`remove-sgroup`:t.menuSGroup(n,`remove`);return;case`bond-properties`:this.#$?.show(n.bondIds,r);return;case`change-direction`:t.menuReverse(n);return;case`charge.plus`:case`charge.minus`:t.menuCharge(n,e===`charge.plus`?1:-1);return;case`hydrogens.add`:case`hydrogens.remove`:t.menuHydrogens(n,e===`hydrogens.add`?`add`:`remove`);return;case`structure.add-hydrogens`:case`structure.remove-hydrogens`:t.menuHydrogens(n,e===`structure.add-hydrogens`?`add`:`remove`);return;case`structure.aromatize`:case`structure.dearomatize`:t.menuAromaticity(n,e===`structure.aromatize`?`aromatize`:`dearomatize`);return;case`select-structure`:t.menuSelectStructure(n);return;case`cut`:this.#Rn(n);return;case`copy`:this.#Ln(t.menuCopy(n),this.#ze.clipboardBlocked(this.#ze.copy,this.#Sn(`copy`)));return;case`delete`:t.menuDelete(n);return;case`paste`:this.#bn();return;case`select-all`:t.selectAll();return;case`clean`:case`straighten`:this.#Rt(e);return;case`copy-readout`:{let e=this.#Qe?.text??``;e&&this.#vn(e);return}case`rotate`:this.#Z?.show({id:`dialog.rotate`,title:this.#ze.rotateDialog,label:this.#ze.rotateAngle,value:90,min:0,max:360,integer:!1,invalid:this.#ze.rotateInvalid,direction:`ccw`,ok:(e,t)=>void this.#_?.rotateSelection(e,t??`ccw`)},r);return;case`zoom.fit`:this.#Lt(`fit`);return;case`enhanced-stereo`:this.#In(n);return;case`stereo-panel`:this.#Lr(!this.#re?.open,`panel`);return;case`chiral-flag`:t.menuChiralFlag();return;case`remove-stereo`:case`recompute-wedges`:t.menuStereoCommand(n,e===`remove-stereo`?`remove`:`rewedge`);return;case`aromatize`:case`dearomatize`:t.menuAromaticity(n,e);return;case`show-stereochemistry`:this.#Nr({...this.settings,stereoLabels:this.stereoLabels===`shown`?`hidden`:`shown`});return;case`expand-abbreviation`:t.menuAbbreviation(n,`expand`);return;case`contract-abbreviation`:t.menuAbbreviation(n,`contract`);return;case`remove-abbreviation`:t.menuAbbreviation(n,`remove`);return;case`contract-label`:t.menuAbbreviation(n,`contract-label`);return}}#In(e){let t=this.#_,n=t?.stereoTargetsOf(e);if(this.#C||!t||!n||!n.atomIds.length&&!n.bondIds.length||this.#ie===`hidden`)return;let r=n.atomIds.map(e=>t.atomIndex(e)).filter(e=>e>=0),i=n.bondIds.map(e=>t.bondIndex(e)).filter(e=>e>=0);t.selectItems({atoms:r,bonds:i}),this.#Lr(!0,`panel`)}#Ln(e,t,n,r=this.#ze.cdxmlNotLoaded,i){if(this.#C)return;let a=(e,t=!0)=>{this.#C||this.#_?.notify(e,t)};if(e===null){let e=n?.();e&&a(e);return}let o=()=>{let e=n?.();e?a(e,!1):this.#C||this.#At()},s=()=>n?.()??r,c=this.ownerDocument?.defaultView,l=c?.navigator?.clipboard;if(typeof e==`string`){l?.writeText?this.#t.alive(l.writeText(e)).then(o,()=>a(t)):a(t);return}i&&a(i,!1);let u=!1,d=e.then(e=>(e===null&&(u=!0),e));if(l?.write&&c?.ClipboardItem){let e=d.then(e=>{if(e===null)throw Error(`not written`);return new c.Blob([e],{type:`text/plain`})});this.#t.alive(l.write([new c.ClipboardItem({"text/plain":e})])).then(o,()=>a(u?s():t))}else l?.writeText?this.#t.alive(d).then(e=>e===null?a(s()):l.writeText(e).then(o)).catch(()=>a(u?s():t)):a(t)}async#Rn(e){let t=this.#_,n=t?.menuCopy(e)??null;if(!t||typeof n!=`string`)return;let r=this.#e.stateToken;try{let e=this.ownerDocument?.defaultView?.navigator?.clipboard;if(!e?.writeText)throw Error(`no asynchronous Clipboard API`);await e.writeText(n)}catch{this.#C||this.#_?.notify(this.#ze.clipboardBlocked(this.#ze.cut,this.#Sn(`cut`)));return}!this.#C&&this.#e.stateToken===r&&this.#_?.menuDelete(e)}#zn(e,t){this.#Bn();let n=this.#t.after(500,()=>this.#Vn());this.#Ne={pointerId:e,at:t,cancel:n}}#Bn(){this.#Ne&&=(this.#Ne.cancel(),null)}#Vn(){let e=this.#Ne;this.#Ne=null,!(!e||this.#C||this.#Ge.size!==1||!this.#Ge.has(e.pointerId))&&(this.#qe===e.pointerId&&this.#Kt(),this.#Pn(`pointer`,e.at),this.#Pe=!!this.#te?.isOpen)}#Hn=e=>{this.#C||q(e)||(this.#v.pointerMoved(e.clientX,e.clientY),this.#j?.pointerMoved(e.clientX,e.clientY),this.#re?.pointerMoved(e.clientX,e.clientY))};#Un=e=>{this.#C||q(e)||(this.#v.pointerGone(),this.#j?.pointerGone(),this.#re?.pointerGone())};#Wn(e,t){if(this.#C||this.#He||!t.isPrimary||t.pointerType===`mouse`&&t.button!==0)return;let n=t.currentTarget;try{n.setPointerCapture(t.pointerId)}catch{return}this.#Ue=null;let r=t.pointerType===`mouse`?null:(n.closest(`.bar`)??n).getBoundingClientRect();this.#He={symbol:e,pointerId:t.pointerId,button:n,start:{x:t.clientX,y:t.clientY},bar:r,dragging:!1,cancelled:!1},n.addEventListener(`pointermove`,this.#Gn),n.addEventListener(`pointerup`,this.#qn),n.addEventListener(`pointercancel`,this.#Jn),n.addEventListener(`touchmove`,this.#Kn,{passive:!1})}#Gn=e=>{let t=this.#He;if(t&&e.pointerId===t.pointerId&&!t.cancelled){if(!t.dragging){let{bar:n}=t;t.dragging=n?e.clientX<n.left||e.clientX>n.right||e.clientY<n.top||e.clientY>n.bottom:Math.hypot(e.clientX-t.start.x,e.clientY-t.start.y)>=4}t.dragging&&this.#_?.dragElement(t.symbol,this.#Zn(e.clientX,e.clientY))}};#Kn=e=>{let t=this.#He,n=e.touches[0];if(t&&n&&e.cancelable){if(!t.intent){let e=n.clientX-t.start.x,r=n.clientY-t.start.y;if(Math.abs(e)<2&&Math.abs(r)<2)return;t.intent=Math.abs(e)>=Math.abs(r)?`drag`:`scroll`}t.intent===`drag`&&e.preventDefault()}};#qn=e=>{let t=this.#He;if(!t||e.pointerId!==t.pointerId||(this.#Yn(),!t.dragging))return;this.#Ue=t.symbol,this.#t.after(0,()=>{this.#Ue===t.symbol&&(this.#Ue=null)});let n=t.cancelled?null:this.#Zn(e.clientX,e.clientY);n?this.#_?.dropElement(t.symbol,n):this.#_?.cancelDrop()};#Jn=()=>{let e=this.#He;this.#Yn(),e?.dragging&&this.#_?.cancelDrop()};#Yn(){let e=this.#He;e&&(this.#He=null,e.button.removeEventListener(`pointermove`,this.#Gn),e.button.removeEventListener(`pointerup`,this.#qn),e.button.removeEventListener(`pointercancel`,this.#Jn),e.button.removeEventListener(`touchmove`,this.#Kn))}#Xn=e=>{let t=this.#a,n=this.#_;if(!t||!n||!this.isConnected||!(e.width>0&&e.height>0))return!1;let r=t.getBoundingClientRect(),i=kn(t),a=nn({x:e.left,y:e.top},r,i),o=nn({x:e.right,y:e.bottom},r,i);return n.drawingReaches({x:Math.min(a.x,o.x),y:Math.min(a.y,o.y),width:Math.abs(o.x-a.x),height:Math.abs(o.y-a.y)})};#Zn(e,t){let n=this.#i,r=this.#a;if(!n||!r||!this.isConnected)return null;let i=r.getBoundingClientRect();if(e<i.left||e>i.right||t<i.top||t>i.bottom)return null;let a=this.renderRoot.elementFromPoint?.(e,t);return a&&!n.contains(a)?null:nn({x:e,y:t},i,kn(r))}#Qn(){let e=this.renderRoot.activeElement;if(!e)return null;if(e===this.#g.svg)return{kind:`canvas`};let t=e.getAttribute(`data-testid`);return t?{kind:`control`,testid:t}:null}#$n(e){let t=this.#a;return!!(t&&e&&t.contains(e))}#er(e,t){let n=this.#i;if(!e||!t||!n||this.#C){this.#J?.remove(),this.#J=null;return}let r=this.#J;r||(r=n.ownerDocument.createElement(`div`),r.className=`sgroup-tooltip`,r.setAttribute(`role`,`tooltip`),r.setAttribute(`data-testid`,`canvas.sgroup-tooltip`),n.append(r),this.#J=r),r.textContent=e;let i=n.clientWidth,a=n.clientHeight,o=r.offsetWidth,s=r.offsetHeight,c=Math.max(4,Math.min(t.x+12,i-o-4)),l=t.y+18+s<=a-4?t.y+18:Math.max(4,t.y-12-s);Object.assign(r.style,{left:`${c}px`,top:`${l}px`})}#tr(e){!e||this.#C||(e.kind===`canvas`?this.#g.svg?.focus({preventScroll:!0}):this.renderRoot.querySelector(`[data-testid="${e.testid}"]`)?.focus({preventScroll:!0}))}#nr(e,t=!1){if(!this.#C){if(this.#P){this.#P.show(e,t);return}this.#R={opener:e,touch:t},this.#rr(),!this.#I&&(this.#I=!0,this.#t.alive(import(`./structure-library-j0nGc2cQ.js`)).then(({StructureLibrary:e})=>{this.#I=!1;let t=this.#F;if(this.#C||!t)return;this.#P=new e(this.#ar(),t);let n=this.#R;this.#R=null,this.#ir(),n&&this.#P.show(n.opener,n.touch)},()=>{this.#I=!1;let e=this.#R;this.#R=null,this.#ir(),!(!e||this.#C)&&(this.#tr(e.opener),this.#_?.notify(this.#ze.libraryNotLoaded))}))}}#rr(){let e=this.#F;if(!e?.isConnected||this.#L)return;let t=e.ownerDocument,n=this.#ze,r=t.createElement(`dialog`);r.className=`dialog library library-wait`,r.setAttribute(`data-testid`,`dialog.library-loading`),r.setAttribute(`aria-labelledby`,`sl-wait-title`);let i=t.createElement(`div`);i.className=`dialog-head`;let a=t.createElement(`h2`);a.id=`sl-wait-title`,a.textContent=n.structureLibrary;let o=t.createElement(`button`);o.type=`button`,o.className=`close`,o.setAttribute(`data-testid`,`dialog.library-loading.close`),o.setAttribute(`aria-label`,n.close),o.textContent=`×`,o.addEventListener(`click`,()=>this.#dr()),i.append(a,o);let s=t.createElement(`p`);s.className=`sl-wait`,s.setAttribute(`role`,`status`),s.setAttribute(`data-testid`,`dialog.library-loading.status`),s.textContent=n.libraryLoading,r.append(i,s),r.addEventListener(`keydown`,e=>{e.key===`Escape`?(e.preventDefault(),this.#dr()):e.key===`Tab`&&kr(e)}),r.addEventListener(`cancel`,e=>{e.preventDefault(),this.#dr()}),e.append(r),this.#L=r,r.showModal(),o.focus()}#ir(){let e=this.#L;this.#L=null,e&&(e.open&&e.close(),e.remove())}#ar(){return{strings:()=>this.#ze,read:(e,t)=>this.#e.inspect(e.molblock,`molV2000`,t,!0),atomColours:()=>this.#Ve,pick:e=>this.#fr(e),restoreFocus:e=>this.#tr(e),groups:()=>this.#e.abbreviationList(),inspect:(e,t,n)=>this.#e.inspect(e,t,n),pickGroup:e=>this.#_?.setGroup(e),pickSalt:e=>this.#_?.setSalt(e)}}#or(e){if(!this.#C){if(this.#H={request:e,typed:``},this.#B){this.#sr();return}this.#U||(this.#U=!0,this.#t.alive(import(`./structure-library-j0nGc2cQ.js`)).then(({AbbreviationLookup:e})=>{this.#U=!1;let t=this.#V;this.#C||!t||(this.#B=new e(this.#lr(),t),this.#sr())},()=>{this.#U=!1,this.#H=null}))}}#sr(){let e=this.#H;this.#H=null;let t=this.#g.svg;if(!e||!this.#B||!t)return;let n=t.getBoundingClientRect(),r=e.request,i=r.at?Wn(r.at,n,kn(t)):{x:n.left+n.width/2,y:n.top+n.height/2},a=this.renderRoot.querySelector(`.frame`)?.getBoundingClientRect(),o=a?{left:a.left,top:a.top,right:a.right,bottom:a.bottom}:null;this.#B.open({mode:r.mode,text:r.text+e.typed,at:i,box:o,...r.atomId===void 0?{}:{atomId:r.atomId}})}#cr(){this.#H=null,this.#B?.close()}#lr(){return{strings:()=>this.#ze,groups:()=>this.#e.abbreviationList(),pick:(e,t)=>this.#ur(e,t.mode,t.atomId),closed:e=>{e&&!this.#C&&this.#g.svg?.focus({preventScroll:!0})}}}#ur(e,t,n){if(!this.#C){if(this.#g.svg?.focus({preventScroll:!0}),t===`nickname`){n!==void 0&&this.#_?.applyNickname(n,String(e.key));return}e.kind===`element`?this.#pr(String(e.key)):e.kind===`group`?this.#_?.setGroup(String(e.key)):e.template?this.#fr(e.template):e.salt&&this.#_?.setSalt(e.salt)}}#dr(){let e=this.#R;this.#R=null,this.#L&&(this.#ir(),e&&!this.#C&&this.#tr(e.opener)),this.#P?.cancel()}#fr(e){let t=this.#e.inspect(e.molblock,`molV2000`,t=>jr(t,e),!0);t&&this.#_?.setTemplate(t)}#pr(e){this.#Be=pg(this.#Be,e,$i(this.#y)),this.#_?.setTool(cr(e)),this.requestUpdate()}get atomColours(){return this.#Ve}set atomColours(e){this.#C||(this.#Tr(`atomColours`),this.#Ve=a(e),this.#Ve!==this.#g.colours&&(this.#g.colours=this.#Ve,this.#e.refresh(),this.requestUpdate(),this.#M?.recolour(),this.#P?.recolour()))}get stereoLabels(){return this.#e.stereoShown.labels?`shown`:`hidden`}set stereoLabels(e){this.#C||(this.#Tr(`stereoLabels`),this.#e.stereoShown={...this.#e.stereoShown,labels:Fy(e)===`shown`})}get stereoLabelStyle(){return this.#e.stereoShown.registration?`registration`:`standard`}set stereoLabelStyle(e){this.#C||(this.#Tr(`stereoLabelStyle`),this.#e.stereoShown={...this.#e.stereoShown,registration:Fy(e)===`registration`})}get stereoFlag(){return this.#e.stereoShown.flag?`shown`:`hidden`}set stereoFlag(e){this.#C||(this.#Tr(`stereoFlag`),this.#e.stereoShown={...this.#e.stereoShown,flag:Fy(e)!==`hidden`})}get unspecifiedDoubleBonds(){return this.#e.stereoShown.plain?`plain`:`crossed`}set unspecifiedDoubleBonds(e){this.#C||(this.#Tr(`unspecifiedDoubleBonds`),this.#e.stereoShown={...this.#e.stereoShown,plain:Fy(e)===`plain`})}get toolbar(){return JSON.parse(JSON.stringify(this.#mr))}set toolbar(e){this.#C||(this.#mr=e&&typeof e==`object`?e:Ui,this.#gr(),this.#ve=Yi(this.#mr),this.#br())}#mr=Ui;#hr(){return[...this.#y.top,...this.#y.left,...this.#y.right,...this.#y.bottom].some(e=>e.items.some(e=>e.kind===`button`&&(e.id===`structure.aromatize`||e.id===`structure.dearomatize`)))}#gr(){let e=new Set(this.#e.mode===`query`?[]:tn);this.#ie===`hidden`&&e.add(`structure.stereo-panel`).add(`structure.enhanced-stereo`),this.#y=Zi(this.#mr,e),this.#v.config=this.#y}#_r(e){let t=e.kind===`generic`?n(e):`[${(e.elements??[]).join(`,`)}]`,r=e.kind===`list`?`queryAtomList`:e.kind===`not-list`?`queryAtomNotList`:`queryAtomGeneric`;return this.#ze.format(r,{label:t})}#vr(e){if(this.#E&&!Fr(e))return this.#ze.format(`viewOnlyOff`,{name:D(e,this.#ze)});let t=(e,t)=>e!==null&&this.#b.includes(e)?this.#ze.format(`queryElementOff`,{name:t}):null;if(e===`element.query`){let e=this.#_?.queryAtom;return e?t(Ce(e),this.#_r(e)):null}let n=e.startsWith(`bond.`)?e.slice(5):``;return pr(n)?t(n,ut(n,this.#ze)):null}get disabledQueryElements(){return this.#b}set disabledQueryElements(e){this.#C||(this.#b=Tr(e),this.#_&&(this.#_.queryOff=new Set(this.#b)),this.#te?.close(),this.requestUpdate())}get smiles(){return this.#e.smiles}get plainSmiles(){return this.#e.plainSmiles}get molfile(){return this.#e.molfile}get molV3000(){return this.#e.molV3000}get cxsmiles(){return this.#e.cxsmiles}get canonicalSmiles(){return this.#e.canonicalSmiles}get canonicalCxsmiles(){return this.#e.canonicalCxsmiles}get smarts(){return this.#e.smarts}get svg(){return this.#e.svg}png(e){try{return lg(this.svg,this.ownerDocument,e?.scale??1)}catch(e){return Promise.reject(e)}}get cdxml(){return this.#e.cdxml}get inchi(){return this.#e.inchi}get inchiKey(){return this.#e.inchiKey}loadFormat(e){return this.#e.loadFormat(e)}get keymap(){return this.#Ie}set keymap(e){if(this.#C)return;this.#Tr(`keymap`);let t=$o(e);t!==this.#Ie&&(this.#Ie=t,this.#br())}get keys(){return this.#_e}set keys(e){this.#C||(this.#_e=As(e),this.#br())}get#yr(){return this.#ye}#br(){let e=this.#e.mode===`query`?[...Qo[this.#Ie],...Zo]:Qo[this.#Ie];this.#ye=js(Ss(e,this.#ve),this.#_e),this.#_&&(this.#_.keymap=this.#ye,this.#_.eraseWhole=this.#Ie===`ketcher`,this.#_.hoverKeys=this.#Ie===`ketcher`),this.#pe?.refresh(),this.requestUpdate()}get methylLabels(){return this.#be}set methylLabels(e){if(this.#C)return;this.#Tr(`methylLabels`);let t=Fy(e),n=t===`shown`||t===`all`?t:`hidden`;n!==this.#be&&(this.#be=n,this.#g.methyl=n!==`hidden`,this.#g.carbons=n===`all`,this.#e.refresh())}#xr(e){this.#C||e===this.#xe||(this.#xe=e,this.#g.numbers=e===`shown`,this.#e.refresh())}get density(){return this.#Se}set density(e){if(this.#C)return;let t=Uh(e);t!==this.#Se&&(this.#Se=t,this.requestUpdate())}#Sr(){let e=this.ownerDocument?.defaultView;if(!e)return;let t=e.getComputedStyle(this),n=(e,t)=>(Number.parseFloat(e)||0)+(Number.parseFloat(t)||0);this.#Cr({width:this.clientWidth-n(t.paddingLeft,t.paddingRight),height:this.clientHeight-n(t.paddingTop,t.paddingBottom)})}#Cr(e){if(this.#C||!(e.width>0&&e.height>0))return!1;let t=Wh(this.#Se,this.#Oe),n=this.#wr();return this.#Oe={width:e.width,height:e.height},Wh(this.#Se,this.#Oe)===t?(this.#wr()!==n&&this.requestUpdate(),!1):(this.requestUpdate(),!0)}#wr(){return Wh(`auto`,this.#Oe)}get settings(){return Oh({atomColours:this.atomColours,warningMarks:this.warningMarks,analysisPanel:this.analysisPanel,stereoPanel:this.stereoPanel,stereoLabels:this.stereoLabels,stereoLabelStyle:this.stereoLabelStyle,stereoFlag:this.stereoFlag,unspecifiedDoubleBonds:this.unspecifiedDoubleBonds,methylLabels:this.methylLabels,atomNumbers:this.#xe,keymap:this.keymap,readoutItems:this.#tt})}set settings(e){if(this.#C)return;if(!this.#Te&&typeof e==`object`&&e)for(let t of Th)t in e&&this.#we.add(t);let t=Oh(e);this.#Er(()=>{this.atomColours=t.atomColours,this.warningMarks=t.warningMarks,this.analysisPanel=t.analysisPanel,this.stereoPanel=t.stereoPanel,this.stereoLabels=t.stereoLabels,this.stereoLabelStyle=t.stereoLabelStyle,this.stereoFlag=t.stereoFlag,this.unspecifiedDoubleBonds=t.unspecifiedDoubleBonds,this.methylLabels=t.methylLabels,this.#xr(t.atomNumbers),this.keymap=t.keymap,this.readoutItems=t.readoutItems})}#Tr(e){this.#Te||this.#we.add(e)}#Er(e){let t=this.#Te;this.#Te=!0;try{e()}finally{this.#Te=t}}get persistSettings(){return this.#Ce}set persistSettings(e){this.#C||(this.#Ce=Bh(e),this.#Ce?this.#Ee&&this.#kr():this.#jr())}#Dr(){this.#Ee||(this.#Ee=!0,this.#Ce&&(this.#Or(Rh(this.ownerDocument?.defaultView)),this.#kr()))}#Or(e){if(this.#C)return;let t=this.settings,n=Ih(t,e,this.#we);kh(t,n)||this.#Er(()=>this.settings=n)}#kr(){if(this.#De||this.#C||!this.#Ce)return;let t=this.ownerDocument?.defaultView?.BroadcastChannel;if(typeof t==`function`)try{let n=new t(Mh);n.onmessage=e.#Ar(new WeakRef(this)),this.#De=n,Ly.register(this,n,n),this.#t.attach(()=>this.#jr())}catch{}}static#Ar(e){return t=>{let n=e.deref();n&&n.#Ce&&n.#Or(Nh(t.data))}}#jr(){let e=this.#De;this.#De=null,e&&(Ly.unregister(e),Ry(e))}#Mr(e){if(this.#Ce&&Object.keys(e).length){zh(this.ownerDocument?.defaultView,e);try{this.#De?.postMessage(Ph(e))}catch{}}}#Nr(e){let t=this.settings;this.#Er(()=>this.settings=e);let n=this.settings;this.#C||kh(t,n)||(this.#Mr(Fh(t,n)),this.#e.settingsApplied(n),this.dispatchEvent(new CustomEvent(`settings`,{detail:n})))}get locale(){return this.#ze.locale}set locale(e){this.#C||(this.#Le=e,this.#Pr())}get strings(){return this.#Re}set strings(e){this.#C||(this.#Re=te(e),this.#Pr())}#Pr(){let e=dt(this.#Le,this.#Re);this.#ze=e,this.#_&&(this.#_.strings=e),this.#e.strings=e,this.#g.label=e.canvasLabel,this.#ut&&this.#pt(),this.#A&&(this.#A.label=e.labelEditor),this.#j?.relabel(),this.#Qe?.relabel(e),this.#_t(),this.#Ot(),this.#It(),this.#Fr(),this.#de?.relabel(),this.#pe?.refresh(),this.#he?.refresh(),this.#re?.refresh(),this.requestUpdate()}#Fr(){for(let e of[this.#N,this.#F,this.#G,this.#q,this.#X,this.#Q,this.#ee,this.#fe,this.#me,this.#ge,this.#ue,this.#V,this.#ne])e?.setAttribute(`lang`,this.#ze.locale)}get value(){return this.#e.value}set value(e){this.#e.value=e}setValue(e,t){this.#e.setValue(e,t)}addMolecule(e,t){if(this.#C)throw Error(gg);let n=e==null?``:String(e),r=t?.format||void 0;if(!this.#e.ready||!this.#_){this.ready.then(()=>{this.#C||this.addMolecule(n,{format:r})},()=>void 0);return}let i=this.#_.addMolecule(n,r);i&&this.#e.refuseRead(i)}get lastError(){return this.#e.lastError}get warnings(){return this.#e.warnings}get readWarnings(){return this.#e.readWarnings}get cxsmilesCoordinates(){return this.#e.cxsmilesCoordinates}set cxsmilesCoordinates(e){this.#C||(this.#e.cxsmilesCoordinates=e)}get analysis(){return this.#e.analysis}get selectionAnalysis(){let e=this.#C?[]:this.#_?.selectedAtoms??[];return e.length?this.#e.analysisOf(e):null}get hasQuery(){return this.#e.hasQuery}get selectionHasQuery(){if(this.#C||!this.#_?.hasSelection)return!1;let{atoms:e,bonds:t}=this.#_.selected;return this.#e.hasQueryIn(e,t)}get warningMarks(){return this.#$e}set warningMarks(e){this.#C||(this.#Tr(`warningMarks`),this.#$e=tp(e),this.#_&&(this.#_.problemMarks=this.#$e===`shown`))}get analysisPanel(){return this.#et}set analysisPanel(e){this.#C||(this.#Tr(`analysisPanel`),this.#et=Nm(e),this.#Qe&&(this.#Qe.state=this.#et),this.#_t())}get readoutItems(){return this.#tt}set readoutItems(e){this.#C||(this.#Tr(`readoutItems`),this.#tt=nh(e),this.#_t())}get stereoPanel(){return this.#ie}set stereoPanel(e){if(this.#C)return;this.#Tr(`stereoPanel`);let t=Sd(e);if(t===this.#ie)return;let n=this.#ie===`hidden`;this.#ie=t,(n||t===`hidden`)&&(this.#gr(),this.#br()),this.#Rr(t===`open`)}get stereo(){if(this.#C)return this.#ce;let e=this.#e.stateToken;if(e&&e===this.#Ir)return this.#ce;let t=this.hasQuery?null:this.#_?.stereoRows();return this.#ce=t?wd(t.info,t.rows,this.#e.stereoisomers()):Cd,this.#Ir=t?e:null,this.#ce}#Ir=null;#Lr(e,t){if(this.#C||this.#ie===`hidden`)return;let n=!e&&t===`canvas`&&!!this.#re?.holdsFocus();this.#Er(()=>this.stereoPanel=e?`open`:`closed`),e&&t===`panel`&&this.#re?.focusFirst(),n&&this.#g.svg?.focus({preventScroll:!0})}#Rr(e){let t=this.#re;if(!(!t||this.#C||e===t.open)){if(e)this.#ae=null,this.#zr(),this.#Vr=this.#e.stateToken,t.open=!0;else{let e=t.holdsFocus();t.open=!1,this.#ae=null,this.#oe=!1,this.#i?.removeAttribute(`data-panel`),this.#i?.style.removeProperty(`--panel-w`),this.#i?.style.removeProperty(`--panel-h`),e&&this.#g.svg?.focus({preventScroll:!0})}this.requestUpdate()}}#zr(){let e=this.#i,t=this.#re;if(!e||!t||this.#C||this.#ie!==`open`)return;let n={width:e.clientWidth,height:e.clientHeight};if(!(n.width>0&&n.height>0))return;let r=xd(n,this.#wr()===`compact`,this.#ae?.place??null),i=this.#ae;i&&i.place===r.place&&i.width===r.width&&i.height===r.height||(this.#ae=r,e.setAttribute(`data-panel`,r.place),t.element.setAttribute(`data-place`,r.place),t.place=r.place,e.style.setProperty(`--panel-w`,`${r.width}px`),e.style.setProperty(`--panel-h`,`${r.height}px`),(!i||i.place!==r.place||(r.place===`right`?i.width!==r.width:i.height!==r.height))&&(this.#oe=!0))}#Br(e){let t=this.#_;if(!t||this.#C)return;let n=e??(t.selectedAtoms.length?t.selectedAtoms:null),r=t.atomPx(n),{width:i,height:a}=this.#g.size??{width:0,height:0};if(!r.length||!(i>0)||!(a>0))return;let o=(e,t,n)=>t-e>n?Math.round(n/2-(e+t)/2):e<0?Math.round(-e):t>n?Math.round(n-t):0,s=r.map(e=>e.x),c=r.map(e=>e.y),l=o(Math.min(...s)-24,Math.max(...s)+24,i),u=o(Math.min(...c)-24,Math.max(...c)+24,a);(l||u)&&this.#We.panBy(l,u)}#Vr=null;#Hr(){!this.#re?.open||this.#se||this.#C||this.#e.stateToken!==this.#Vr&&(this.#se=this.#t.painted(()=>{this.#se=null,this.#Vr=this.#e.stateToken,this.#re?.refresh()}))}#Ur(){let e=this.#_;if(!e||this.#C)return null;let t=this.hasQuery,n=e.stereoRows();return n?{info:n.info,rows:t?[]:n.rows,count:t?null:this.#e.stereoisomers(),query:t,acts:this.#E?{remove:!1,rewedge:!1}:e.panelCommandsAct(),readOnly:this.#E}:null}#Wr(e){let t=rs(this.#yr,e,`generic`)?.action;return t!==`undo`&&t!==`redo`||this.#E||this.#C?!1:(t===`undo`?this.#e.undo(`user`):this.#e.redo(`user`),!0)}#Gr(){this.#M?.cancelAll(),this.#dr(),this.#cr(),this.#W?.cancelAll(),this.#K?.cancelAll(),this.#Y?.cancel(),this.#Z?.cancel(),this.#$?.cancel(),this.#le?.cancel(),this.#te?.close()}get viewOnly(){return this.#E}set viewOnly(e){let t=e===!0||e===`true`;if(this.#C||t===this.#E)return;this.#E=t;let n=this.#_;t&&(this.#Gr(),this.#Bn(),n&&!Fr(n.tool)&&(this.#D=n.tool,n.setTool(`select.rect`))),n&&(n.viewOnly=t),!t&&n&&this.#D&&(n.tool===`select.rect`&&n.setTool(this.#D),this.#D=null),this.#re?.refresh(),this.requestUpdate()}get isSaved(){return this.#e.isSaved}markSaved(){this.#e.markSaved()}get isEmpty(){return this.#e.isEmpty}get mode(){return this.#e.mode}set mode(e){if(this.#C)return;let t=this.#e.mode;this.#e.mode=e,this.#e.mode!==t&&(this.#te?.close(),this.#gr(),this.#br(),this.#_?.modeChanged(),this.requestUpdate())}get version(){return this.#e.version}get engineRev(){return this.#e.engineRev}get renderCount(){return this.#e.renderCount}get positions(){let e=this.#C?null:this.#g.svg;if(!e||!this.#_)return lo;let t=e.getBoundingClientRect(),n=this.getBoundingClientRect();return this.#_.positions({x:t.left-n.left,y:t.top-n.top},this.#g.moved)}get isPending(){return this.#C?!1:this.#t.busy||!!this.#_?.busy||this.#We.busy}get layout(){let e=this.renderRoot;if(this.#C||!this.#e.ready||!this.isConnected||!e)return sy;let t=e=>e?.getBoundingClientRect()??null,n=[];for(let r of Array.from(e.querySelectorAll(Sy))){let e=(r.getAttribute(`data-testid`)??``).slice(8);n.push({name:e,box:r.getBoundingClientRect(),holder:t(r.closest(`.palette, .bar`))})}return dy({origin:this.getBoundingClientRect(),canvas:t(this.#g.svg??this.#a),toolbars:yy.map(n=>[n,t(e.querySelector(`.bar-${n}`))]),controls:n,labelEditor:t(this.#A?.element)})}get selection(){return this.#C||!this.#e.ready||!this.#_?cy:fy(this.#_.selected)}#Kr(){if(this.#C||!this.#e.ready)return;let e=this.selection,t=JSON.stringify(e);t!==this.#O&&(this.#O=t,this.#e.selectionChanged(e))}get tool(){return this.#_?.tool??this.#ct}undo(){return this.#e.undo(`api`)}redo(){return this.#e.redo(`api`)}get canUndo(){return this.#e.canUndo}get canRedo(){return this.#e.canRedo}on(e,t){return this.#e.on(e,t)}off(e,t){this.#e.off(e,t)}destroy(){if(!this.#C){this.#C=!0,this.#t.end(),this.#n=null,this.#S=null,this.#x?.disconnect(),this.#x=null,this.#We.cancel(),this.#Ge.clear(),this.#Ke.clear(),this.#Yn(),this.#v.destroy(),this.#M?.cancelAll(),this.#M=null,this.#dr(),this.#P=null,this.#H=null,this.#B?.destroy(),this.#B=null,this.#W?.cancelAll(),this.#W=null,this.#K?.cancelAll(),this.#K=null,this.#er(null,null),this.#Y?.cancel(),this.#Y=null,this.#Z?.cancel(),this.#Z=null,this.#$?.cancel(),this.#$=null,this.#re=null,this.#de?.cancel(),this.#de=null,this.#pe?.cancel(),this.#pe=null,this.#he?.cancel(),this.#he=null,this.#le?.cancel(),this.#le=null,this.#te?.destroy(),this.#te=null,this.#je=null,this.#Bn(),this.#A?.hide(),this.#A=null,this.#Qe=null,this.#j?.destroy(),this.#j=null,this.#ct=this.#_?.tool??this.#ct,this.#_?.destroy(),this.#_=null;for(let[e,t]of this.#k)this.#a?.removeEventListener(e,t);this.#k=[],this.#e.destroy(),this.#a=null,this.#i=null,this.#o=null,this.#s=!1,this.#c=null,this.#l=null,this.#u.destroy(),this.#d=null,this.#f=null,this.#p=null,this.#m=null,this.#h=null,this.#N=null,this.#F=null,this.#G=null,this.#q=null,this.#X=null,this.#Q=null,this.#ee=null,this.#fe=null,this.#me=null,this.#ge=null,this.#ue=null,this.#V=null,this.#ne=null,this.renderRoot&&T(w,this.renderRoot,this.renderOptions),Iy(this)}}},Fy=e=>e==null?``:String(e).replace(/[A-Z]/g,e=>e.toLowerCase());function Iy(e){e.ownerDocument?.documentElement?.getBoundingClientRect()}var Ly=new FinalizationRegistry(Ry);function Ry(e){e.onmessage=null;try{e.close()}catch{}}function zy(e){let t=e.composedPath()[0],n=typeof t?.tagName==`string`?t.tagName.toLowerCase():``;return n===`input`||n===`textarea`||!!t?.isContentEditable}function $(e){let t=e.createElement(`div`);return t.className=`layer`,t}function By(e){try{e.currentTarget.setPointerCapture(e.pointerId)}catch{}}function Vy(e){e.cancelable&&e.preventDefault()}function Hy(e){let t=e?.navigator,n=t?.userAgentData?.platform||t?.platform||``;return/mac|iphone|ipad/i.test(n)?Br.cmd:Br.ctrl}function Uy(e,t){let n=null;for(let r of t)Object.prototype.hasOwnProperty.call(e,r)&&((n??={})[r]=e[r],delete e[r]);return n}customElements.get(`crux-sketch`)||customElements.define(py,Py);var Wy=class{#e;#t;#n;constructor(e){this.#t=e,this.#e=oy(e),this.#n=e.ownerDocument}get smiles(){return this.#e.smiles}get plainSmiles(){return this.#e.plainSmiles}get molfile(){return this.#e.molfile}get molV3000(){return this.#e.molV3000}get cxsmiles(){return this.#e.cxsmiles}get canonicalSmiles(){return this.#e.canonicalSmiles}get canonicalCxsmiles(){return this.#e.canonicalCxsmiles}get smarts(){return this.#e.smarts}get svg(){return this.#e.svg}png(e){try{let t=this.#t?.ownerDocument??this.#n;return lg(this.svg,t,e?.scale??1)}catch(e){return Promise.reject(e)}}get cdxml(){return this.#e.cdxml}get inchi(){return this.#e.inchi}get inchiKey(){return this.#e.inchiKey}loadFormat(e){return this.#e.loadFormat(e)}get value(){return this.#e.value}set value(e){this.#e.value=e}setValue(e,t){this.#e.setValue(e,t)}addMolecule(e,t){if(!this.#t)throw Error(gg);this.#t.addMolecule(e,t)}get lastError(){return this.#e.lastError}get warnings(){return this.#e.warnings}get readWarnings(){return this.#e.readWarnings}get cxsmilesCoordinates(){return this.#e.cxsmilesCoordinates}set cxsmilesCoordinates(e){this.#e.cxsmilesCoordinates=e}get analysis(){return this.#e.analysis}get selectionAnalysis(){return this.#t?.selectionAnalysis??null}get isSaved(){return this.#e.isSaved}markSaved(){this.#e.markSaved()}get isEmpty(){return this.#e.isEmpty}get hasQuery(){return this.#e.hasQuery}get selectionHasQuery(){return this.#t?.selectionHasQuery??!1}get mode(){return this.#e.mode}set mode(e){this.#t&&(this.#t.mode=e)}get disabledQueryElements(){return this.#t&&(this.#r=this.#t.disabledQueryElements),this.#r}set disabledQueryElements(e){this.#t&&(this.#t.disabledQueryElements=e,this.#r=this.#t.disabledQueryElements)}#r=Object.freeze([]);get version(){return this.#e.version}get engineRev(){return this.#e.engineRev}get renderCount(){return this.#e.renderCount}get positions(){return this.#t?.positions??lo}get isPending(){return this.#t?.isPending??!1}get layout(){return this.#t?.layout??sy}get selection(){return this.#t?.selection??cy}get tool(){return this.#t&&(this.#i=this.#t.tool),this.#i}#i=`bond.single`;undo(){return this.#e.undo(`api`)}redo(){return this.#e.redo(`api`)}get canUndo(){return this.#e.canUndo}get canRedo(){return this.#e.canRedo}get atomColours(){return this.#t&&(this.#a=this.#t.atomColours),this.#a}set atomColours(e){this.#t&&(this.#t.atomColours=e,this.#a=this.#t.atomColours)}#a=`element`;get keymap(){return this.#t&&(this.#o=this.#t.keymap),this.#o}set keymap(e){this.#t&&(this.#t.keymap=e,this.#o=this.#t.keymap)}#o=`chemdraw`;get warningMarks(){return this.#t&&(this.#s=this.#t.warningMarks),this.#s}set warningMarks(e){this.#t&&(this.#t.warningMarks=e,this.#s=this.#t.warningMarks)}#s=`shown`;get analysisPanel(){return this.#t&&(this.#c=this.#t.analysisPanel),this.#c}set analysisPanel(e){this.#t&&(this.#t.analysisPanel=e,this.#c=this.#t.analysisPanel)}#c=`shown`;get readoutItems(){return this.#t&&(this.#l=this.#t.readoutItems),this.#l}set readoutItems(e){this.#t&&(this.#t.readoutItems=e,this.#l=this.#t.readoutItems)}#l=wh.readoutItems;get stereoPanel(){return this.#t&&(this.#u=this.#t.stereoPanel),this.#u}set stereoPanel(e){this.#t&&(this.#t.stereoPanel=e,this.#u=this.#t.stereoPanel)}#u=`closed`;get stereo(){return this.#t&&(this.#d=this.#t.stereo),this.#d}#d=Cd;get stereoLabels(){return this.#t&&(this.#f=this.#t.stereoLabels),this.#f}set stereoLabels(e){this.#t&&(this.#t.stereoLabels=e,this.#f=this.#t.stereoLabels)}#f=`hidden`;get stereoLabelStyle(){return this.#t&&(this.#p=this.#t.stereoLabelStyle),this.#p}set stereoLabelStyle(e){this.#t&&(this.#t.stereoLabelStyle=e,this.#p=this.#t.stereoLabelStyle)}#p=`standard`;get stereoFlag(){return this.#t&&(this.#m=this.#t.stereoFlag),this.#m}set stereoFlag(e){this.#t&&(this.#t.stereoFlag=e,this.#m=this.#t.stereoFlag)}#m=`shown`;get unspecifiedDoubleBonds(){return this.#t&&(this.#h=this.#t.unspecifiedDoubleBonds),this.#h}set unspecifiedDoubleBonds(e){this.#t&&(this.#t.unspecifiedDoubleBonds=e,this.#h=this.#t.unspecifiedDoubleBonds)}#h=`crossed`;get locale(){return this.#t&&(this.#g=this.#t.locale),this.#g}set locale(e){this.#t&&(this.#t.locale=e,this.#g=this.#t.locale)}#g=`en`;get strings(){return this.#t&&(this.#_=this.#t.strings),this.#_}set strings(e){this.#t&&(this.#t.strings=e,this.#_=this.#t.strings)}#_=Object.freeze({});get methylLabels(){return this.#t&&(this.#v=this.#t.methylLabels),this.#v}set methylLabels(e){this.#t&&(this.#t.methylLabels=e,this.#v=this.#t.methylLabels)}#v=`hidden`;get keys(){return this.#t&&(this.#y=this.#t.keys),this.#y}set keys(e){this.#t&&(this.#t.keys=e,this.#y=this.#t.keys)}#y=Object.freeze({});get settings(){return this.#t&&(this.#b=this.#t.settings),this.#b}set settings(e){this.#t&&(this.#t.settings=e,this.#b=this.#t.settings)}#b=wh;get density(){return this.#t&&(this.#x=this.#t.density),this.#x}set density(e){this.#t&&(this.#t.density=e,this.#x=this.#t.density)}#x=Vh;get persistSettings(){return this.#t&&(this.#S=this.#t.persistSettings),this.#S}set persistSettings(e){this.#t&&(this.#t.persistSettings=e,this.#S=this.#t.persistSettings)}#S=!0;get viewOnly(){return this.#t&&(this.#C=this.#t.viewOnly),this.#C}set viewOnly(e){this.#t&&(this.#t.viewOnly=e,this.#C=this.#t.viewOnly)}#C=!1;on(e,t){return this.#e.on(e,t)}off(e,t){this.#e.off(e,t)}destroy(){let e=this.#t;e&&(this.#i=e.tool),this.#t=null,e&&(e.remove(),e.destroy())}};async function Gy(e,t={}){if(!e||typeof e.append!=`function`)throw TypeError(`Crux Sketch: createSketcher needs an element to mount into`);let n=document.createElement(py);if(!(n instanceof Py))throw Error(`Crux Sketch: <${py}> is defined by another script`);if(t.mode!==void 0&&(n.mode=t.mode),t.toolbar!==void 0&&(n.toolbar=t.toolbar),t.atomColours!==void 0&&(n.atomColours=t.atomColours),t.keymap!==void 0&&(n.keymap=t.keymap),t.warningMarks!==void 0&&(n.warningMarks=t.warningMarks),t.analysisPanel!==void 0&&(n.analysisPanel=t.analysisPanel),t.readoutItems!==void 0&&(n.readoutItems=t.readoutItems),t.stereoPanel!==void 0&&(n.stereoPanel=t.stereoPanel),t.stereoLabels!==void 0&&(n.stereoLabels=t.stereoLabels),t.stereoLabelStyle!==void 0&&(n.stereoLabelStyle=t.stereoLabelStyle),t.stereoFlag!==void 0&&(n.stereoFlag=t.stereoFlag),t.unspecifiedDoubleBonds!==void 0&&(n.unspecifiedDoubleBonds=t.unspecifiedDoubleBonds),t.locale!==void 0&&(n.locale=t.locale),t.strings!==void 0&&(n.strings=t.strings),t.methylLabels!==void 0&&(n.methylLabels=t.methylLabels),t.keys!==void 0&&(n.keys=t.keys),typeof t.settings==`object`&&t.settings!==null){let e=t.settings,r=n;for(let t of Th)t in e&&(t in n?r[t]=e[t]:ay(n,t,e[t]))}t.cxsmilesCoordinates!==void 0&&(n.cxsmilesCoordinates=t.cxsmilesCoordinates),t.density!==void 0&&(n.density=t.density),t.persistSettings!==void 0&&(n.persistSettings=t.persistSettings),t.disabledQueryElements!==void 0&&(n.disabledQueryElements=t.disabledQueryElements),t.viewOnly!==void 0&&(n.viewOnly=t.viewOnly),t.value!==void 0&&(n.value=t.value),e.append(n);try{await n.ready}catch(e){throw n.destroy(),n.remove(),e}return new Wy(n)}var Ky=new WeakMap;async function qy(){let e;try{e=await zi()}catch(e){throw Za(e)}let t=Ky.get(e);if(t)throw t;return e}function Jy(e,t){try{return t(e)}catch(t){if(typeof WebAssembly<`u`&&t instanceof WebAssembly.RuntimeError)throw Ky.has(e)||Ky.set(e,$a(t)),Ky.get(e);let n=t&&typeof t==`object`&&typeof t.message==`string`?t.message:String(t);throw Object.freeze({format:`inchi`,message:n})}}async function Yy(e){return e==null?``:Jy(await qy(),t=>t.inchiFromInput(e)).inchi}async function Xy(e){return e?Jy(await qy(),t=>t.inchiToInchiKey(e)):``}export{Py as CruxSketchElement,Gy as createSketcher,Yy as inchiFromInput,Xy as inchiKeyFromInchi,Ig as version};