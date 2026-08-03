	
	
	var UDOC = {};
	
	UDOC.G = {
		concat : function(p,r) {
			for(var i=0; i<r.cmds.length; i++) p.cmds.push(r.cmds[i]);
			for(var i=0; i<r.crds.length; i++) p.crds.push(r.crds[i]);
		},
		getBB  : function(ps) {
			var x0=1e99, y0=1e99, x1=-x0, y1=-y0;
			for(var i=0; i<ps.length; i+=2) {  var x=ps[i],y=ps[i+1];  if(x<x0)x0=x; else if(x>x1)x1=x;  if(y<y0)y0=y;  else if(y>y1)y1=y;  }
			return [x0,y0,x1,y1];
		},
		rectToPath: function(r) {  return  {cmds:["M","L","L","L","Z"],crds:[r[0],r[1],r[2],r[1], r[2],r[3],r[0],r[3]]};  },
		// a inside b
		insideBox: function(a,b) {  return b[0]<=a[0] && b[1]<=a[1] && a[2]<=b[2] && a[3]<=b[3];   },
		isBox : function(p, bb) {
			var sameCrd8 = function(pcrd, crds) {
				for(var o=0; o<8; o+=2) {  var eq = true;  for(var j=0; j<8; j++) if(Math.abs(crds[j]-pcrd[(j+o)&7])>=2) {  eq = false;  break;  }    if(eq) return true;  }
				return false;
			};
			if(p.cmds.length>10) return false;
			var cmds=p.cmds.join(""), crds=p.crds;
			var sameRect = false;
			if((cmds=="MLLLZ"  && crds.length== 8) 
			 ||(cmds=="MLLLLZ" && crds.length==10) ) {
				if(crds.length==10) crds=crds.slice(0,8);
				var x0=bb[0],y0=bb[1],x1=bb[2],y1=bb[3];
				if(!sameRect) sameRect = sameCrd8(crds, [x0,y0,x1,y0,x1,y1,x0,y1]);
				if(!sameRect) sameRect = sameCrd8(crds, [x0,y1,x1,y1,x1,y0,x0,y0]);
			}
			return sameRect;
		},
		boxArea: function(a) {  var w=a[2]-a[0], h=a[3]-a[1];  return w*h;  },
		newPath: function(gst    ) {  gst.pth = {cmds:[], crds:[]};  },
		moveTo : function(gst,x,y) {  var p=UDOC.M.multPoint(gst.ctm,[x,y]);  //if(gst.cpos[0]==p[0] && gst.cpos[1]==p[1]) return;
										gst.pth.cmds.push("M");  gst.pth.crds.push(p[0],p[1]);  gst.cpos = p;  },
		lineTo : function(gst,x,y) {  var p=UDOC.M.multPoint(gst.ctm,[x,y]);  if(gst.cpos[0]==p[0] && gst.cpos[1]==p[1]) return;
										gst.pth.cmds.push("L");  gst.pth.crds.push(p[0],p[1]);  gst.cpos = p;  },
		curveTo: function(gst,x1,y1,x2,y2,x3,y3) {   var p;  
			p=UDOC.M.multPoint(gst.ctm,[x1,y1]);  x1=p[0];  y1=p[1];
			p=UDOC.M.multPoint(gst.ctm,[x2,y2]);  x2=p[0];  y2=p[1];
			p=UDOC.M.multPoint(gst.ctm,[x3,y3]);  x3=p[0];  y3=p[1];  gst.cpos = p;
			gst.pth.cmds.push("C");  
			gst.pth.crds.push(x1,y1,x2,y2,x3,y3);  
		},
		closePath: function(gst  ) {  gst.pth.cmds.push("Z");  },
		arc : function(gst,x,y,r,a0,a1, neg) {
			
			// circle from a0 counter-clock-wise to a1
			if(neg) while(a1>a0) a1-=2*Math.PI;
			else    while(a1<a0) a1+=2*Math.PI;
			var th = (a1-a0)/4;
			
			var x0 = Math.cos(th/2), y0 = -Math.sin(th/2);
			var x1 = (4-x0)/3, y1 = y0==0 ? y0 : (1-x0)*(3-x0)/(3*y0);
			var x2 = x1, y2 = -y1;
			var x3 = x0, y3 = -y0;
			
			var p0 = [x0,y0], p1 = [x1,y1], p2 = [x2,y2], p3 = [x3,y3];
			
			var pth = {cmds:[(gst.pth.cmds.length==0)?"M":"L","C","C","C","C"], crds:[x0,y0,x1,y1,x2,y2,x3,y3]};
			
			var rot = [1,0,0,1,0,0];  UDOC.M.rotate(rot,-th);
			
			for(var i=0; i<3; i++) {
				p1 = UDOC.M.multPoint(rot,p1);  p2 = UDOC.M.multPoint(rot,p2);  p3 = UDOC.M.multPoint(rot,p3);
				pth.crds.push(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1]);
			}
			
			var sc = [r,0,0,r,x,y];  
			UDOC.M.rotate(rot, -a0+th/2);  UDOC.M.concat(rot, sc);  UDOC.M.multArray(rot, pth.crds);
			UDOC.M.multArray(gst.ctm, pth.crds);
			
			UDOC.G.concat(gst.pth, pth);
			var y=pth.crds.pop();  x=pth.crds.pop();
			gst.cpos = [x,y];
		},
		toPoly : function(p) {
			if(p.cmds[0]!="M" || p.cmds[p.cmds.length-1]!="Z") return null;
			for(var i=1; i<p.cmds.length-1; i++) if(p.cmds[i]!="L") return null;
			var out = [], cl = p.crds.length;
			if(p.crds[0]==p.crds[cl-2] && p.crds[1]==p.crds[cl-1]) cl-=2;
			for(var i=0; i<cl; i+=2) out.push([p.crds[i],p.crds[i+1]]);
			if(UDOC.G.polyArea(p.crds)<0) out.reverse();
			return out;
		},
		fromPoly : function(p) {
			var o = {cmds:[],crds:[]};
			for(var i=0; i<p.length; i++) { o.crds.push(p[i][0], p[i][1]);  o.cmds.push(i==0?"M":"L");  }
			o.cmds.push("Z");
			return o;
		},
		polyArea : function(p) {
			if(p.length <6) return 0;
			var l = p.length - 2;
			var sum = (p[0]-p[l]) * (p[l+1]+p[1]);
			for(var i=0; i<l; i+=2)
				sum += (p[i+2]-p[i]) * (p[i+1]+p[i+3]);
			return - sum * 0.5;
		},
		polyClip : function(p0, p1) {  // p0 clipped by p1
            var cp1, cp2, s, e;
            var inside = function (p) {
                return (cp2[0]-cp1[0])*(p[1]-cp1[1]) > (cp2[1]-cp1[1])*(p[0]-cp1[0]);
            };
            var isc = function () {
                var dc = [ cp1[0] - cp2[0], cp1[1] - cp2[1] ],
                    dp = [ s[0] - e[0], s[1] - e[1] ],
                    n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0],
                    n2 = s[0] * e[1] - s[1] * e[0], 
                    n3 = 1.0 / (dc[0] * dp[1] - dc[1] * dp[0]);
                return [(n1*dp[0] - n2*dc[0]) * n3, (n1*dp[1] - n2*dc[1]) * n3];
            };
            var out = p0;
            cp1 = p1[p1.length-1];
            for (j in p1) {
                var cp2 = p1[j];
                var inp = out;
                out = [];
                s = inp[inp.length - 1]; //last on the input list
                for (i in inp) {
                    var e = inp[i];
                    if (inside(e)) {
                        if (!inside(s)) {
                            out.push(isc());
                        }
                        out.push(e);
                    }
                    else if (inside(s)) {
                        out.push(isc());
                    }
                    s = e;
                }
                cp1 = cp2;
            }
            return out
        }
	}
	UDOC.M = {
		getScale : function(m) {  return Math.sqrt(Math.abs(m[0]*m[3]-m[1]*m[2]));  },
		translate: function(m,x,y) {  UDOC.M.concat(m, [1,0,0,1,x,y]);  },
		rotate   : function(m,a  ) {  UDOC.M.concat(m, [Math.cos(a), -Math.sin(a), Math.sin(a), Math.cos(a),0,0]);  },
		scale    : function(m,x,y) {  UDOC.M.concat(m, [x,0,0,y,0,0]);  },
		concat   : function(m,w  ) {  
			var a=m[0],b=m[1],c=m[2],d=m[3],tx=m[4],ty=m[5];
			m[0] = (a *w[0])+(b *w[2]);       m[1] = (a *w[1])+(b *w[3]);
			m[2] = (c *w[0])+(d *w[2]);       m[3] = (c *w[1])+(d *w[3]);
			m[4] = (tx*w[0])+(ty*w[2])+w[4];  m[5] = (tx*w[1])+(ty*w[3])+w[5]; 
		},
		invert   : function(m    ) {  
			var a=m[0],b=m[1],c=m[2],d=m[3],tx=m[4],ty=m[5], adbc=a*d-b*c;
			m[0] = d/adbc;  m[1] = -b/adbc;  m[2] =-c/adbc;  m[3] =  a/adbc;
			m[4] = (c*ty - d*tx)/adbc;  m[5] = (b*tx - a*ty)/adbc;
		},
		multPoint: function(m, p ) {  var x=p[0],y=p[1];  return [x*m[0]+y*m[2]+m[4],   x*m[1]+y*m[3]+m[5]];  },
		multArray: function(m, a ) {  for(var i=0; i<a.length; i+=2) {  var x=a[i],y=a[i+1];  a[i]=x*m[0]+y*m[2]+m[4];  a[i+1]=x*m[1]+y*m[3]+m[5];  }  }
	}
	UDOC.C = {
		srgbGamma : function(x) {  return x < 0.0031308 ? 12.92 * x : 1.055 * Math.pow(x, 1.0 / 2.4) - 0.055;  },
		cmykToRgb : function(clr) { 
			var c=clr[0], m=clr[1], y=clr[2], k=clr[3];
			// return [1-Math.min(1,c+k), 1-Math.min(1, m+k), 1-Math.min(1,y+k)];
			var r = 255
			+ c * (-4.387332384609988  * c + 54.48615194189176  * m +  18.82290502165302  * y + 212.25662451639585 * k +  -285.2331026137004) 
			+ m * ( 1.7149763477362134 * m - 5.6096736904047315 * y + -17.873870861415444 * k - 5.497006427196366) 
			+ y * (-2.5217340131683033 * y - 21.248923337353073 * k +  17.5119270841813) 
			+ k * (-21.86122147463605  * k - 189.48180835922747);
			var g = 255
			+ c * (8.841041422036149   * c + 60.118027045597366 * m +  6.871425592049007  * y + 31.159100130055922 * k +  -79.2970844816548) 
			+ m * (-15.310361306967817 * m + 17.575251261109482 * y +  131.35250912493976 * k - 190.9453302588951) 
			+ y * (4.444339102852739   * y + 9.8632861493405    * k -  24.86741582555878) 
			+ k * (-20.737325471181034 * k - 187.80453709719578);
			var b = 255
			+ c * (0.8842522430003296  * c + 8.078677503112928  * m +  30.89978309703729  * y - 0.23883238689178934 * k + -14.183576799673286) 
			+ m * (10.49593273432072   * m + 63.02378494754052  * y +  50.606957656360734 * k - 112.23884253719248) 
			+ y * (0.03296041114873217 * y + 115.60384449646641 * k + -193.58209356861505)
			+ k * (-22.33816807309886  * k - 180.12613974708367);

			return [Math.max(0, Math.min(1, r/255)), Math.max(0, Math.min(1, g/255)), Math.max(0, Math.min(1, b/255))];
			//var iK = 1-c[3];  
			//return [(1-c[0])*iK, (1-c[1])*iK, (1-c[2])*iK];  
		},
		labToRgb  : function(lab) {
			var k = 903.3, e = 0.008856, L = lab[0], a = lab[1], b = lab[2];
			var fy = (L+16)/116, fy3 = fy*fy*fy;
			var fz = fy - b/200, fz3 = fz*fz*fz;
			var fx = a/500 + fy, fx3 = fx*fx*fx;
			var zr = fz3>e ? fz3 : (116*fz-16)/k;
			var yr = fy3>e ? fy3 : (116*fy-16)/k;
			var xr = fx3>e ? fx3 : (116*fx-16)/k;
				
			var X = xr*96.72, Y = yr*100, Z = zr*81.427, xyz = [X/100,Y/100,Z/100];
			var x2s = [3.1338561, -1.6168667, -0.4906146, -0.9787684,  1.9161415,  0.0334540, 0.0719453, -0.2289914,  1.4052427];
			
			var rgb = [ x2s[0]*xyz[0] + x2s[1]*xyz[1] + x2s[2]*xyz[2],
						x2s[3]*xyz[0] + x2s[4]*xyz[1] + x2s[5]*xyz[2],
						x2s[6]*xyz[0] + x2s[7]*xyz[1] + x2s[8]*xyz[2]  ];
			for(var i=0; i<3; i++) rgb[i] = Math.max(0, Math.min(1, UDOC.C.srgbGamma(rgb[i])));
			return rgb;
		}
	}
	
	UDOC.getState = function(crds) {
		return {
			font : UDOC.getFont(),
			dd: {flat:1},  // device-dependent
			space :"/DeviceGray",
			// fill
			ca: 1,
			colr  : [0,0,0],
			sspace:"/DeviceGray",
			// stroke
			CA: 1,
			COLR : [0,0,0],
			bmode: "/Normal",
			SA:false, OPM:0, AIS:false, OP:false, op:false, SMask:"/None",
			lwidth : 1,
			lcap: 0,
			ljoin: 0,
			mlimit: 10,
			SM : 0.1,
			doff: 0,
			dash: [],
			ctm : [1,0,0,1,0,0],
			cpos: [0,0],
			pth : {cmds:[],crds:[]}, 
			cpth: crds ? UDOC.G.rectToPath(crds) : null  // clipping path
		};
	}
	
	UDOC.getFont = function() {
		return {
			Tc: 0, // character spacing
			Tw: 0, // word spacing
			Th:100, // horizontal scale
			Tl: 0, // leading
			Tf:"Helvetica-Bold", 
			Tfs:1, // font size
			Tmode:0, // rendering mode
			Trise:0, // rise
			Tk: 0,  // knockout
			Tal:0,  // align, 0: left, 1: right, 2: center
			Tun:0,  // 0: no, 1: underline
			
			Tm :[1,0,0,1,0,0],
			Tlm:[1,0,0,1,0,0],
			Trm:[1,0,0,1,0,0]
		};
	}