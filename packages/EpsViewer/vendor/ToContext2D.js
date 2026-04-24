

	function ToContext2D(needPage, scale)
	{
		this.canvas = document.createElement("canvas");
		this.ctx = this.canvas.getContext("2d");
		this.bb = null;
		this.currPage = 0;
		this.needPage = needPage;
		this.scale = scale;
	}
	ToContext2D.prototype.StartPage = function(x,y,w,h) {
		if(this.currPage!=this.needPage) return;
		this.bb = [x,y,w,h];
		var scl = this.scale, dpr = window.devicePixelRatio;
		var cnv = this.canvas, ctx = this.ctx;
		cnv.width = Math.round(w*scl);  cnv.height = Math.round(h*scl);
		ctx.translate(0,h*scl);  ctx.scale(scl,-scl);
		cnv.setAttribute("style", "border:1px solid; width:"+(cnv.width/dpr)+"px; height:"+(cnv.height/dpr)+"px");
	}
	ToContext2D.prototype.Fill = function(gst, evenOdd) {
		if(this.currPage!=this.needPage) return;
		var ctx = this.ctx;
		ctx.beginPath();
		this._setStyle(gst, ctx);
		this._draw(gst.pth, ctx);
		ctx.fill();
	}
	ToContext2D.prototype.Stroke = function(gst) {
		if(this.currPage!=this.needPage) return;
		var ctx = this.ctx;
		ctx.beginPath();
		this._setStyle(gst, ctx);
		this._draw(gst.pth, ctx);
		ctx.stroke();
	}
	ToContext2D.prototype.PutText = function(gst, str, stw) {
		if(this.currPage!=this.needPage) return;
		var scl = this._scale(gst.ctm);
		var ctx = this.ctx;
		this._setStyle(gst, ctx);
		ctx.save();
		var m = [1,0,0,-1,0,0];  this._concat(m, gst.font.Tm);  this._concat(m, gst.ctm);
		//console.log(str, m, gst);  throw "e";
		ctx.transform(m[0],m[1],m[2],m[3],m[4],m[5]);
		ctx.fillText(str,0,0);
		ctx.restore();
	}
	ToContext2D.prototype.PutImage = function(gst, buff, w, h, msk) {
		if(this.currPage!=this.needPage) return;
		var ctx = this.ctx;
		
		if(buff.length==w*h*4) {
			buff = buff.slice(0);
			if(msk && msk.length==w*h*4) for(var i=0; i<buff.length; i+=4) buff[i+3] = msk[i+1];
			
			var cnv = document.createElement("canvas"), cctx = cnv.getContext("2d");
			cnv.width = w;  cnv.height = h;
			var imgd = cctx.createImageData(w,h);
			for(var i=0; i<buff.length; i++) imgd.data[i]=buff[i];
			cctx.putImageData(imgd,0,0);
			
			ctx.save();
			var m = [1,0,0,1,0,0];  this._concat(m, [1/w,0,0,-1/h,0,1]);  this._concat(m, gst.ctm);
			ctx.transform(m[0],m[1],m[2],m[3],m[4],m[5]);
			ctx.drawImage(cnv,0,0);
			ctx.restore();
		}
	}
	ToContext2D.prototype.ShowPage = function() {  this.currPage++;  }
	ToContext2D.prototype.Done = function() {}
	
	
	ToContext2D.prototype._setStyle = function(gst, ctx) {
		var scl = this._scale(gst.ctm);
		ctx.fillStyle = this._getFill(gst.colr, gst.ca, ctx);
		ctx.strokeStyle=this._getFill(gst.COLR, gst.CA, ctx);
		
		ctx.lineCap = ["butt","round","square"][gst.lcap];
		ctx.lineJoin= ["miter","round","bevel"][gst.ljoin];
		ctx.lineWidth=gst.lwidth*scl;
		var dsh = gst.dash.slice(0);  for(var i=0; i<dsh.length; i++) dsh[i] = ToPDF._flt(dsh[i]*scl);
		ctx.setLineDash(dsh); 
		ctx.miterLimit = gst.mlimit*scl;
		
		var fn = gst.font.Tf, ln = fn.toLowerCase();
		var p0 = ln.indexOf("bold")!=-1 ? "bold " : "";
		var p1 = (ln.indexOf("italic")!=-1 || ln.indexOf("oblique")!=-1) ? "italic " : "";
		ctx.font = p0+p1 + gst.font.Tfs+"px \""+fn+"\"";
	}
	ToContext2D.prototype._getFill = function(colr, ca, ctx)
	{
		if(colr.typ==null) return this._colr(colr,ca);
		else {
			var grd = colr, crd = grd.crds, mat = grd.mat, scl=this._scale(mat), gf;
			if     (grd.typ=="lin") {
				var p0 = this._multPoint(mat,crd.slice(0,2)), p1 = this._multPoint(mat,crd.slice(2));
				gf=ctx.createLinearGradient(p0[0],p0[1],p1[0],p1[1]);
			}
			else if(grd.typ=="rad") {
				var p0 = this._multPoint(mat,crd.slice(0,2)), p1 = this._multPoint(mat,crd.slice(3));
				gf=ctx.createRadialGradient(p0[0],p0[1],crd[2]*scl,p1[0],p1[1],crd[5]*scl);
			}
			for(var i=0; i<grd.grad.length; i++)  gf.addColorStop(grd.grad[i][0],this._colr(grd.grad[i][1], ca));
			return gf;
		}
	}
	ToContext2D.prototype._colr  = function(c,a) {  return "rgba("+Math.round(c[0]*255)+","+Math.round(c[1]*255)+","+Math.round(c[2]*255)+","+a+")";  };
	ToContext2D.prototype._scale = function(m)  {  return Math.sqrt(Math.abs(m[0]*m[3]-m[1]*m[2]));  };
	ToContext2D.prototype._concat= function(m,w  ) {  
			var a=m[0],b=m[1],c=m[2],d=m[3],tx=m[4],ty=m[5];
			m[0] = (a *w[0])+(b *w[2]);       m[1] = (a *w[1])+(b *w[3]);
			m[2] = (c *w[0])+(d *w[2]);       m[3] = (c *w[1])+(d *w[3]);
			m[4] = (tx*w[0])+(ty*w[2])+w[4];  m[5] = (tx*w[1])+(ty*w[3])+w[5]; 
	}
	ToContext2D.prototype._multPoint= function(m, p) {  var x=p[0],y=p[1];  return [x*m[0]+y*m[2]+m[4],   x*m[1]+y*m[3]+m[5]];  },
	ToContext2D.prototype._draw  = function(path, ctx)
	{
		var c = 0, crds = path.crds;
		for(var j=0; j<path.cmds.length; j++) {
			var cmd = path.cmds[j];
			if     (cmd=="M") {  ctx.moveTo(crds[c], crds[c+1]);  c+=2;  }
			else if(cmd=="L") {  ctx.lineTo(crds[c], crds[c+1]);  c+=2;  }
			else if(cmd=="C") {  ctx.bezierCurveTo(crds[c], crds[c+1], crds[c+2], crds[c+3], crds[c+4], crds[c+5]);  c+=6;  }
			else if(cmd=="Q") {  ctx.quadraticCurveTo(crds[c], crds[c+1], crds[c+2], crds[c+3]);  c+=4;  }
			else if(cmd=="Z") {  ctx.closePath();  }
		}
	}