#include "com.h"
#include <vector>
// ################################################################################################
// rendering
// ################################################################################################
#define SZ 	4
#define SX 	(1.0 / IMAGESCALE)
#define SY 	(1.0 / IMAGESCALE)

#define FOREACHI(len)	for(int i = 0; i < len; i ++)
#define FOREACHJ(wid)	for(int j = 0; j < wid; j ++)
#define BLEND1(i, sz)	blend(-1, 1, i / real(sz));

vec3 campos(0.5, 0.5, -CAM_DIS);

vec2 rndvecmap[512][512];
vec2 rndvecmap1[512][512];
real rndmap[512][512];
int rndcormap[512][512];
int rndcormap1[512][512];
real rndmapx1[SZ+2];
real rndmapx2[SZ+2];
real rndmapy1[SZ+2];
real rndmapy2[SZ+2];
real alpha;
vec2 pt;
vec3 pt3;
int paint;
vec2 o, v;
real level = 0.15;
int fog = RGB(153, 190, 228);
// ------------------------------------------------
  void init()
{
	int i, j;
	for(i = 0; i < 512; i ++)
	for(j = 0; j < 512; j ++)
	{
		rndmap[i][j] = rrnd(0, 1);
		rndvecmap[i][j][0] = rrnd(0.01, .99);
		rndvecmap[i][j][1] = rrnd(0.01, .99);
		rndvecmap1[i][j][0] = rrnd(0, 1);
		rndvecmap1[i][j][1] = rrnd(0, 1);		
		
		rndcormap[i][j] = RGB(rrnd(0, 55), rrnd(0, 255), rrnd(0, 55));
		rndcormap1[i][j] = RGB(rrnd(0, 55), rrnd(0, 255), rrnd(0, 55));
	}
	rndmapx1[i] = 0;
	rndmapx2[i] = 0; 
	rndmapy1[i] = 0;
	rndmapy2[i] = 0; 
	for(int i = 0; i < SZ+2; i ++)	
	{
		rndmapx1[i] += i / float(SZ) * rrnd(0.8, 1.2);
		rndmapx2[i] += i / float(SZ) * rrnd(0.8, 1.2);
		rndmapy1[i] += i / float(SZ) * rrnd(0.8, 1.2);
		rndmapy2[i] += i / float(SZ) * rrnd(0.8, 1.2);	
	}
	fheightmap = fheightmap0;
}
// -------------------------------------------------
real dot(const vector2& v1, const vector2& v2)
{
	return v1.x * v2.x + v1.y * v2.y;
}
real fract(real v)
{
	return v - int(v);
}
float hash(vec2 p) {
	float h = dot(p, vec2(127.1, 311.7));
	return fract(sin(h)*43758.5453123);
}
// -------------------------------------------------
 void bk()
{
	fheightmap = fheightmap1;
	int i, j;
	real alpha, alpha_;
	for(j = IMAGESCALE - 1; j > 0; j --)
	//for(j = 0; j < IMAGESCALE; j ++)	
	for(i = 0; i < IMAGESCALE; i ++)	
	{
		alpha = blend2d(0, 1, i / FIMAGESCALE, j / FIMAGESCALE);
		alpha_ = blend2d(0, 1, j / FIMAGESCALE, i / FIMAGESCALE);
		alpha = blend2d(0, 1, alpha, alpha_);
		alpha = blend2d(0, 1, alpha_, alpha);
		// real time = blend(100, 1, j / FIMAGESCALE);
		// int cor = blend(0xFFff0000, 0xFF0000ff, alpha);
		// cor = blendcor(1, cor, 0.25 * alpha_);
		real hashh = hash(vec2(i, j));
		drawmap[i][j] =  blendcor(RGB(202, 199, 194), blendcor(RGB(83, 69, 56), RGB(167, 164, 159), hashh), (j / FIMAGESCALE), 0.5);
	//	drawmap[i][j] =  blendcor(RGB(216, 216, 216), blendcor(RGB(192, 192, 192), RGB(172, 172, 172), hashh), (j / FIMAGESCALE), 0.5);
		 
	//	 drawmap[i][j] = blendcor(drawmap[i][j], 0xFFFFFFFF, 0.25);
		// vec2 p(j / FIMAGESCALE, i / FIMAGESCALE);
		// real dis = level - (p.y - level);
		
		// int deta = FIMAGESCALE * 0.01 * blend(0, 1, abs(dis) / 0.1, 0.25) * sin( dis * PI * 180) * alpha;
		// if(j >= int(0.25 * IMAGESCALE))
			// drawmap[i][j] = cor;
		// else
			// drawmap[i][j] = drawmap[i][int(0.25 * IMAGESCALE + abs(j - 0.25 * IMAGESCALE)) % IMAGESCALE];
		//drawmap[i][j] = 0xFFFFFFFF;
	}
	for(j = 0; j < IMAGESCALE; j ++)	
	for(i = 0; i < IMAGESCALE; i ++)	
	{
		alpha = blend2d(0, 1, i / FIMAGESCALE, j / FIMAGESCALE);
		alpha_ = blend2d(0, 1, j / FIMAGESCALE, i / FIMAGESCALE);
		alpha = blend2d(0, 1, alpha, alpha_);
		alpha = blend2d(0, 1, alpha_, alpha);
		real pos =  0.15;
		if(j < pos * IMAGESCALE)
		{
			  drawmap[i][j] = blendcor(drawmap[i][int(pos * 2 * IMAGESCALE - j)], blendcor(1, RGB(56, 43, 34), alpha), (pos * IMAGESCALE - j) / (pos * IMAGESCALE));
			 
		}
		// pos = blend(0.25, 0.35, alpha);
		// if(j < pos * IMAGESCALE)
		// {
			 // drawmap[i][j] = blendcor(RGB(49, 34, 27), drawmap[i][j], j / (pos * IMAGESCALE));
		// }
	}
	
	fheightmap = fheightmap0;
}
// -------------------------------------------------
vec2 bezier2(vec2 cp[3], real t)
{
	return cp[0] * ((1 - t) * (1 - t)) + cp[1] * (2 * t * (1 - t)) + cp[2] * (t * t);
}
// -------------------------------------------------
vec2 bezier3(vec2 cp[4], real t)
{
	real s = 1 - t;
	return  cp[0] * (s * s * s) +
		cp[1] * (3 * t * s * s) +
		cp[2] * (3 * t * t * s) +
		cp[3] * (t * t * t);
}

real drawputaosz8 = 1;
bool drawrelfect8 = 0;
// -------------------------------------------------
void drawputao8(vec3 p, vec3 v)
{		 
	int cor0 = blendcor(0xFFFF80AC, 0xFFCC00FF, rrnd(0, 1));
	vec3 lgtdir(1, -1, 1); lgtdir.norm();
	v.norm();
	for(int i = 0; i < 100; i ++)
	{
		real ai = real(i) / 100;
		p = p + v * (0.000125 * drawputaosz8);		

		real rr = 0.008 * drawputaosz8 * CircleBlend2(0, 1, ai, 1);
		
		vec3 lstipp1[1000];
		vec3 lstiipp1;
		vec3 norm0 = v.cross(vector3::UZ);norm0.norm();
		for(int ii = 0; ii < 100; ii ++)
		{	
			real aii = real(ii) / 100;	
			real ang1 = blend(0, 2 * PI, aii);				
			quaternion q(ang1, v);			
			vec3 norm = q * norm0;
			
			vec3 pp = p + norm * rr;
			
			int cor;
			{// color
				norm = -(lstipp1[ii] - pp).cross(lstiipp1 - pp);
				norm.norm();
				real lum = -lgtdir.dot(norm);	
				
				cor = blendcor(cor0, 0xFF0000FF, blend2d(0, 1, 0.5 + 0.5 * sin(pp.x * 8*PI*ai), 0.5 + 0.5*sin(pp.y * 8*PI*aii)));
				cor = blendcor(0xFF00FF00, cor, ai, .25);
				cor = blendcor(1, cor, 0.5 + 0.5 * lum, 0.5);
				
				lstiipp1 = pp;
				lstipp1[ii] = pp;
			}
			if(i > 0 && ii > 0)
			{
				pixel(viewprj(pp), cor, 1, pp.z);
				
			}			
		}		
	}
}

// ------------------------------------------------'
real sszz = 1;
void bigleaf3(vec3 p0, vec3 v)
{
	vec3 p = p0;
	float s_t0[] = {rrnd(0, 1), rrnd(0, 1), rrnd(0, 1), rrnd(0, 1)};
	real sz = sszz*rrnd(0.4, 0.5);
	real s3d = rrnd(0.01, 0.4);
	vec3 lgtdir(1, -1, 1); lgtdir.norm();
	v.norm();
	vec3 g(0, -0.0001, 0);
	vec3 up = vec3(rrnd(-0.5, 0.5), rrnd(-0.5, 0.5), rrnd(0.25, 1));
	up.norm();
	real deta = rrnd(.5, 1);
	int corx = RGB(rrnd(128, 198), rrnd(18, 198), rrnd(0, 98));

	for(int i = 0; i < 2000; i ++)
	{
		real ai = real(i) / 2000;
		v = v - up * blend(0.04, 0, ai, 1);
		v.norm();
		p = p + v * (0.0001 * sz);
		vec3 norm0 = v.cross(up);norm0.norm();
		up = norm0.cross(v);
		
		vec3 dp = norm0;
		up.norm();
		
		real alpha = blend2d(0.5, 0.8, p.x, p.y);
		
		vec3 pp0 = p0;
		//real szi = blend(0.0, 1, ai, 1);						
		//pp0 = pp0 + up * (szi * sz * 0.02);
		
		vec3 vv = dp;		
		
		int wid0 = 1000;
		//int wid = wid0 * blend(0, 1, ai) * blend(0.5, 1, 0.5 + 0.5 * sin(ai * PI * 8));
		int wid = wid0 * blend2(0, 1, ai, 0.5);// * blend(0.5, 1, 0.5 + 0.5 * abs(sin(ai * PI * 8))); 
		for(int ii = 0; ii < wid; ii ++)
		{
			real aii = real(ii) / wid;	
			//vv.rot(-PI * blend(0, 0.05, aii, 8) * blend(0, 2, ai - 0.25, 2), v);
			//vv.rot(-PI * blend(0, 0.05, real(ii) / wid0, 8) * blend(0, 2, aii - 0.5, 2), v);
			vv.rot(-PI * blend(0, 0.005, real(ii) / wid0, 4), v);
			//vv = vv + up * 0.001;
			vv.norm();
			
			pp0 = pp0 + vv * (0.0003);
			
			real szii = blend(0, 16, real(ii) / wid0, 0.5) * blend(0.0, 2, ai, 1);
			
			vec3 pp = pp0 + up * (szii * sz * 0.005);
			
			{
				vec3 norm = (up + norm0);
				norm.norm();
				//real lum = -lgtdir.dot(norm);			
				real alpha2d = blend2d(0.0, 1.0, pp.x, pp.y);
				alpha2d = blend2d(0.0, 1.0, alpha2d, pp.y);
				alpha2d = blend2d(0.0, 1.0, pp.x, alpha2d);
				int cor00 = blendcor(0xFFFFFFFF, RGB(180, 180, 247), alpha2d);
				cor00 = blendcor(cor00, RGB(255, 255, 255), aii, 0.25);
				int cor = blendcor(cor00, RGB(180, 180, 247), 0.5 + 0.5 * sin(alpha * 1000 * ai));	
				//if(norm.z > 0)
				//	cor = blendcor(cor, RGB(1, 1, 1), aii);
				//else
					
					cor = blendcor(cor, RGB(67, 43, 93), aii, 16);
				//if(i % 8 == 0)
					pixel(pp, cor);		
			}
		}
	}
}

// ------------------------------------------------
int color0 = 0xFFFFFFF;
int color1 = RGB(188, 64, 82);
real meiguisize = 1;
void meigui(vec3 p0, vec3 v0)
{
	vec3 p = p0;
	
	real sz = sszz*rrnd(0.4, 0.5) * 1.2;
	real s3d = rrnd(0.01, 0.4);
	vec3 lgtdir(1, -1, 1); lgtdir.norm();
	
	vec3 g(0, -0.0001, 0);
	vec3 up = v0.cross(vec3::UY);
	up.norm();
	
	v0.norm();
	vec3 v = up.cross(v0);
	v.norm();
	
	real deta = rrnd(.5, 1);
	int corx = RGB(rrnd(128, 198), rrnd(18, 198), rrnd(0, 98));

	for(int i = 0; i < 1000; i ++)
	{
		real ai = real(i) / 1000;
		v = v - up * 0.01;//blend(0.04, 0, ai, 1);
		v.norm();
		p = p + v * (0.0001 * sz);
		vec3 norm0 = v.cross(up);norm0.norm();
		up = norm0.cross(v);
		
		vec3 dp = norm0;
		up.norm();
		
		//real alpha = blend2d(0.5, 0.8, p.x, p.y);
		
		vec3 pp0 = p0;
		//real szi = blend(0.0, 1, ai, 1);						
		//pp0 = pp0 + up * (szi * sz * 0.02);
		
		vec3 vv = dp;		
		
		int wid0 = 300 * meiguisize;
		//int wid = wid0 * blend(0, 1, ai) * blend(0.5, 1, 0.5 + 0.5 * sin(ai * PI * 8));
		int wid = wid0 * blend(0.25, 1, 0.4 + 0.6 * abs(sin(ai * PI * 2 * 8))); 
		for(int ii = 0; ii < wid; ii ++)
		{
			real aii = real(ii) / wid;	
			//vv.rot(-PI * blend(0, 0.05, aii, 8) * blend(0, 2, ai - 0.25, 2), v);
			//vv.rot(-PI * blend(0, 0.05, real(ii) / wid0, 8) * blend(0, 2, aii - 0.5, 2), v);
			vv.rot(-PI * blend(0, meiguisize*0.02, real(ii) / wid0, 4), v);
			//vv = vv + up * 0.001;
			vv.norm();
			
			pp0 = pp0 + vv * (0.00038);
			
			real szii = meiguisize*blend(0, 16, real(ii) / wid0, 0.35);// * blend(0.0, 2, ai, 1);
			
			vec3 pp = pp0 + up * (szii * sz * 0.01);
			
			{
				vec3 norm = (up + norm0);
				norm.norm();
				//real lum = -lgtdir.dot(norm);			
				real alpha2d = blend2d(0.0, 1.0, pp.x, pp.y);
				alpha2d = blend2d(0.0, 1.0, alpha2d, pp.y);
				alpha2d = blend2d(0.0, 1.0, pp.x, alpha2d);
				int cor00 = blendcor(0xFFFFFFFF, color0, alpha2d);
				cor00 = blendcor(cor00, color0, aii, 4);
				int cor = blendcor(cor00, color1, 0.5 + 0.5 * sin(alpha * 1000 * ai * alpha2d));	
				//if(norm.z > 0)
				//	cor = blendcor(cor, RGB(1, 1, 1), aii);
				//else
					
					cor = blendcor(cor, RGB(67, 0, 0), aii, 2);
				//if(i % 8 == 0)
					pixel(pp, cor);		
			}
		}
	}
}

// ------------------------------------------------
void leavejing(vec3 p, vec3 v, int len, real wid0)
{	
	vec3 vcp4 = v.rotcopy(rrnd(-PI, PI), vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)).normcopy());
	vcp4.z = -1;
	vcp4.norm();
	vec3 cp1[] = {
		v,
		v.rotcopy(rrnd(-PI, PI), vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)).normcopy()),
		v.rotcopy(rrnd(-PI, PI), vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)).normcopy()),
		vcp4
	};	
	
	v.norm();
	vec3 v0 = v;
	vec3 aa;
	int lst = 250;
	int step = rrnd(0.1, 0.12) * len;	
	for(int i = 0; i < len; i ++)
	{
		real ai = real(i) / len;
		
		real wid = i < len * 0.9 ? blend(wid0, wid0 * 0.25, i / (len * 0.9), 0.5) :  blend(wid0 * 0.25, wid0 * 4, (i - len * 0.9) / (0.1 * len), 8);
		
		if(i == len - 1)
		{
			sszz = blend(1, 0.5, ai);			
			meigui(p, v);
			step = rrnd(0.1, 0.12) * len;
			lst = i;
		}		
		if(rrnd(0, 15000) < 50)
		{
			int leaflen = 1000 * 0.5;
			//v = vec3(0.5, 0.5, 0.25);v.norm();
			vec3 zz = vec3::UZ;	
			
			{		
				//drawputao8(p + vec3(rrnd(-1, 1),rrnd(-1, 1), rrnd(-1, 1)).normcopy() * 0.01, vec3(rrnd(-1, 1),rrnd(-1, 1), rrnd(-1, 1)));
			}
		}
		//if(depth > 0)
		{	
			v = bezier2(cp1, ai);
			v = v + vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)) * 0.01;
		}
		v.norm();
		p = p + v * 0.00046;
		for(int ii = 0; ii < 400; ii ++)
		{
			real aii = real(ii) / 400;
			real deta1 = blend2d(0, 1, ai, aii);
			real ang1 = blend(0, 2 * PI, aii);	
			
			vec3 norm = v.cross(vec3::UZ);			
			norm.rot(ang1 < PI ? ang1 : ang1 - PI, v);
			real alpha1 = blend2d(0.95, 1.05, aii, ai);
			
			vec3 pp = p + norm * (alpha1 * wid);
			real alpha2 = blend2d(0.8, 1, aii, ai);
			//int cor = blendcor(RGB(13, 15, 40), RGB(57, 58, 80), (0.5 + 0.5 * sin(alpha2 * aii * PI * 1)) * 1);	
			
			int cor = blendcor(RGB(200, 200, 200), RGB(158, 58, 58), (0.5 + 0.5 * sin(alpha2 * aii * PI * 1)) * 1);
			cor = blendcor(cor, RGB(58, 58, 58), real(i - lst) / step, 8);
			cor = blendcor(1, cor, (0.5 + 0.5 * norm.x) * rrnd(0, 2));
			
			pixelAA(pp, cor);
		}	
	}
}

// -------------------------------------------------
void vase(vec3 p, vec3 v)
{
	vec3 lgtdir(1, -1, 1); lgtdir.norm();
	v.norm();
	real ra = 0, rv = 0.00, rr = 0.2;
	for(int i = 0; i < 1600; i ++)
	{
		p = p + v * 0.000125;
		
		real rrr = rr;
		if(i > 1550)
			rrr *= 1;// + 0.1 * abs(sin(PI * (i - 1550) / 50.0));
		else if(i < 80)
		{
			rr += rv * 0.01;
			rrr = rr;
		}
		else
		{
			real deta = blend2d(-1, 1, p.y, p.x);
			//rv += deta * 0.0000001;
			rr += rv * 0.5;
		}	
		
		vec3 norm0 = v.cross(vector3::UZ);norm0.norm();
		for(int ii = 0; ii < 2000; ii ++)
		{	
			real deta = blend2d(0.5, 1.25, i / 1600.0, ii / 2000.0);
			real rrr1 = rrr * blend2(1, 1.414 * deta, real(ii % 400) / 400);
			
			vec3 norm = norm0;
			real ang1 = blend(0, 2 * PI, real(ii) / 2000);				
			quaternion q(ang1, v);
			real ang2 = blend(PI / 2, PI / 2, real(i) / 1600, 1);
			norm = v.rotcopy(ang2, q * norm);	
			norm.norm();	
			vec3 pp = p + norm * (rrr1 * 0.32);
			
			int cor;
			{// color			rrnd(200, 255), rrnd(115, 215)	
				cor = blendcor(RGB(40, 45, 41), RGB(146, 159, 113), blend2d(0, 1, 0.5 + 0.5 * sin(2*PI*real(i) / 1600), 0.5 + 0.5*sin(2*PI*real(ii) / 2000)));
				cor = blendcor(RGB(45, 50, 45), cor, real(i) / 1600, .25);
				cor = blendcor(1, cor, abs(deta * norm.z), 2);
			}
			if(i > 0 && ii > 0)
			{
				pixel(pp, cor, 1);

				
				pp.y = - pp.y;
				pixel(pp, blendcor(1, cor, 0.5), 0.5);
			}			
		}		
	}
}

// -------------------------------------------------
void lan(vec3 p, vec3 v, real len, real wid, real sw, int cor)
{	
	vec3 cp1[] = {
		v,
		v.rotcopy(rrnd(-PI, PI), vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)).normcopy()),
		v.rotcopy(rrnd(-PI, PI), vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)).normcopy()),
		v.rotcopy(rrnd(-PI, PI), vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)).normcopy())
	};		
	
	v.norm();
	vec3 norm0 = vec3::UY;
	vec3 G = vec3(0, -0.001, 0);	
	vec3 VG;
	int cori = int(p.x * 100)%10;
	for(int i = 0; i < len; i ++)
	{
		real alpha1 = blend(0, 1, i / len);
		real a2d = blend2d(.1, 1, p.x * 0.75, p.y * 0.75);
		v = bezier2(cp1, alpha1);
		VG = VG + G;
		v = v + VG * a2d;
		v.norm();
		p = p + v * 0.00052;		
		
		real ang = blend(0, -PI, alpha1);
		vec3 N = v.cross(norm0);N.norm();
		norm0 = N.cross(v);
		
		vec3 vv = N.rotcopy(PI / 2, v);
				
		real ws = sw * blend2(0, 0.02, alpha1, 2);
		
		int cor0 = blendcor(cor, RGB(92, 123, 55), i / len);
		int black = blendcor(RGB(68, 134, 2), RGB(40, 45, 41), i / len, .2);
		for(int j = 0; j < wid; j ++)
		{
			real alpha2 = blend(1, 0, j / wid, 2);						
			
			//if(p.y > 0)
			{
				real hashh = hash(vec2(i, j));
				
				vec3 pp1 = blend(p, p - vv * (ws * blend2d(0.5, 1.5, p.x * 0.75, p.y * 0.75) + blend(0, 0.0025, hashh)), alpha2);
				vec3 pp2 = blend(p, p + vv * (ws * blend2d(0.5, 1.5, p.x * 0.75, p.y * 0.75) + blend(0, 0.0025, hashh)), alpha2);
				
				real aalpha1 = blend2d(0.75, 1.1, pp1.x, pp1.y);				
				real aalpha2 = blend2d(0.75, 1.1, pp2.x, pp2.y);				
				
				pixel(p, 1);				
				
				int cor1 = blendcor(RGB(120, 81, 82), cor0, aalpha1*aalpha1*aalpha1*aalpha1, 0.1);
				cor1 = blendcor(cor1, black, 0.5 + 0.5 * sin(aalpha1 * alpha2 * 10));
				
				int cor2 = blendcor(RGB(120, 81, 82), cor0, aalpha2*aalpha2*aalpha2*aalpha2, 0.1);
				cor2 = blendcor(cor2, black, 0.5 + 0.5 * sin(aalpha2 * alpha2 * 10));
				
				pixelAA(pp1, cor1);
				pixelAA(pp2, cor2);
			}
		}		
	}
}

// ------------------------------------------------
real leaf2size = 1;
void bigleaf2(vec3 p, vec3 v)
{
	real sz = rrnd(0.25, 0.5) * 3 * leaf2size;
	
	vec3 lgtdir(1, -1, 1); lgtdir.norm();
	v.norm();
	vec3 g(0, -0.0001, 0);
	vec3 up = vec3(rrnd(-1, 1), rrnd(0.5, 1), rrnd(-1, 1));
	up.norm();
	real deta = rrnd(.5, 1);
	
	int cor1 = RGB(rrnd(200, 255), rrnd(115, 215), rrnd(0, 5));
	int cor2 = RGB(rrnd(0, 52), rrnd(18, 28), rrnd(7, 17));
	int cor3 = RGB(rrnd(1, 82), rrnd(8, 18), rrnd(2, 4));
	int cor4 = RGB(rrnd(80, 112), rrnd(8, 80), rrnd(0, 2));
	
	real curve = rrnd(0.001, 0.005);
	for(int i = 0; i < 1500; i ++)
	{
		real ai = real(i) / 1500;
		v = v + g;
		v = v - up * blend(curve, 0, ai, 1);
		v.norm();
		p = p + v * (0.00006 * sz);
		vec3 norm0 = v.cross(up);norm0.norm();
		up = norm0.cross(v);
		real alpha = blend2d(0.4, 1.2, p.x, p.y);
		vec3 p1 = p + (up * (0.5 * alpha) + norm0) * CircleBlend2(0, 0.02 * sz * alpha, ai, 1);
		vec3 p2 = p + (up * (0.5 * alpha) - norm0) * CircleBlend2(0, 0.02 * sz * alpha, ai, 1);
		
		for(int ii = 0; ii < 250; ii += 1)
		{
			real aii = real(ii) / 250;								
			{
				vec3 pp = blend(p, p1, aii);				
				
				vec3 norm = (up * 0.5 + norm0);
				norm.norm();
				//real lum = -lgtdir.dot(norm);			
		
				int cor = blendcor(cor1, cor2, 0.5 + 0.5 * sin(alpha * 80 * (0.75 * ai - 0.25 * aii)));	
				if(norm.z > 0)
					cor = blendcor(cor, cor3, aii);
				else
					cor = blendcor(cor, cor4, aii);
				
				//cor = blendcor(1, cor, 0.5 - 0.5 * up.z, 1);					
				
				pixel(pp, cor);		
	
			}
			{
				vec3 pp = blend(p, p2, aii);		
				alpha = blend2d(0.8, 1.2, pp.x, pp.y);
				vec3 norm = (up - norm0);
				norm.norm();
				//real lum = -lgtdir.dot(norm);			

				int cor = blendcor(cor1, cor2, 0.5 + 0.5 * sin(alpha * 80 * (0.75 * ai - 0.25 * aii)));	
				if(norm.z > 0)
					cor = blendcor(cor, cor3, aii);
				else
					cor = blendcor(cor, cor4, aii);	

				// cor = blendcor(1, cor, 0.5 - 0.5 * norm.z, 1);

				pixel(pp, cor);
			 }					
		}
	}
}

void curvedleaf(vec3 p0, vec3 v0)
{
	vec3 p = p0;
	
	real sz = sszz*rrnd(0.4, 0.5) * 0.8 * leaf2size;
	real s3d = rrnd(0.01, 0.4);
	vec3 lgtdir(1, -1, 1); lgtdir.norm();
	
	vec3 g(0, -0.0001, 0);
	vec3 up = v0.cross(vec3::UZ);
	up.y = -abs(up.y);
	up.norm();

	vec3 v = v0;
	v.norm();
	
	real deta = rrnd(.5, 1);
	int corx = RGB(rrnd(128, 198), rrnd(18, 198), rrnd(0, 98));
	int cor1 = RGB(rrnd(200, 255), rrnd(115, 215), rrnd(0, 5));
	int cor2 = RGB(rrnd(0, 52), rrnd(18, 28), rrnd(7, 17));
	int cor3 = RGB(rrnd(1, 82), rrnd(8, 18), rrnd(2, 4));
	int cor4 = RGB(rrnd(80, 112), rrnd(8, 80), rrnd(0, 2));
	real cur1 = rrnd(-0.2, 0.2);
	real cur2 = rrnd(0.01, 0.05);
	for(int i = 0; i < 800; i ++)
	{
		real ai = real(i) / 800;
		v = v - up * blend(-0.0025, cur1, ai, 4) - vec3(0, 0.001, 0);
		v.norm();
		p = p + v * (0.0005 * sz);
		vec3 norm0 = v.cross(up);norm0.norm();
		up = norm0.cross(v);
		
		vec3 dp = norm0;
		up.norm();
		
		//real alpha = blend2d(0.5, 0.8, p.x, p.y);
		
		vec3 pp01 = p;
		vec3 pp02 = p;
		//real szi = blend(0.0, 1, ai, 1);						
		//pp0 = pp0 + up * (szi * sz * 0.02);
		
		vec3 vv1 = dp;		
		vec3 vv2 = -dp;	
		
		int wid0 = 300;
		//int wid = wid0 * blend(0, 1, ai) * blend(0.5, 1, 0.5 + 0.5 * sin(ai * PI * 8));
		
		int wid = wid0;
		real s1 = CircleBlend(0, 1, ai);
		real s2 = blend(1, 0, ai, 0.01);		
		wid *= blend(s1, s2, ai, 1);		
	
		
		for(int ii = 0; ii < wid; ii ++)
		{
			real aii = real(ii) / wid;	
			
			vv1.rot(-PI * blend(0, cur2, real(ii) / wid0, 4), v);
			vv1.norm();
			vv2.rot(PI * blend(0, cur2, real(ii) / wid0, 4), v);
			vv2.norm();
			
			pp01 = pp01 + vv1 * (0.0002);
			pp02 = pp02 + vv2 * (0.0002);
			
			real szii = leaf2size*blend(0, 16, real(ii) / wid0, 2);// * blend(0.0, 2, ai, 1);
			
			vec3 pp1 = pp01 - up * (szii * sz * 0.01);
			vec3 pp2 = pp02 - up * (szii * sz * 0.01);
			{
				vec3 norm = (up + norm0);
				norm.norm();
				//real lum = -lgtdir.dot(norm);			
				int cor = blendcor(cor1, cor2, 0.5 + 0.5 * sin(100 * (0.75 * ai - 0.25 * aii)));	
				//if(norm.z > 0)
					cor = blendcor(cor, cor3, aii, 4);
				//else
				//	cor = blendcor(cor, cor4, aii);	
				
					cor = blendcor(cor, 1, ai, 4);
					
				//if(i % 8 == 0)
					pixel(pp1, cor);	
					pixel(pp2, cor);		
			}
		}
	}
}

// ------------------------------------------------
void bigjing(vec3 p, vec3 v, int len, real wid0)
{	
	vec3 cp1[] = {
		v,
		v.rotcopy(rrnd(-PI, PI), vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)).normcopy()),
		v.rotcopy(rrnd(-PI, PI), vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)).normcopy()),
		v.rotcopy(rrnd(-PI, PI), vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)).normcopy())
	};	
	
	v.norm();
	vec3 v0 = v;
	vec3 aa;
	int lst = 250;
	int step = rrnd(0.1, 0.12) * len;	
	for(int i = 0; i < len; i ++)
	{
		real ai = real(i) / len;
		
		real wid = i < len * 0.9 ? blend(wid0, wid0 * 0.25, i / (len * 0.9), 0.5) :  blend(wid0 * 0.25, wid0 * 1, (i - len * 0.9) / (0.1 * len), 8);
		
		if(i == len - 1)
		{					
			bigleaf2(p, v);
			step = rrnd(0.1, 0.12) * len;
			lst = i;
		}		
		if(rrnd(0, 15000) < 50)
		{
			int leaflen = 1000 * 0.5;
			//v = vec3(0.5, 0.5, 0.25);v.norm();
			vec3 zz = vec3::UZ;	
			
			{		
				//drawputao8(p + vec3(rrnd(-1, 1),rrnd(-1, 1), rrnd(-1, 1)).normcopy() * 0.01, vec3(rrnd(-1, 1),rrnd(-1, 1), rrnd(-1, 1)));
			}
		}
		//if(depth > 0)
		{	
			v = bezier2(cp1, ai);
			v = v + vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)) * 0.01;
			//v = v + vec3(0, -0.0001, 0);
		}
		v.norm();
		p = p + v * 0.00046;
		for(int ii = 0; ii < 400; ii ++)
		{
			real aii = real(ii) / 400;
			real deta1 = blend2d(0, 1, ai, aii);
			real ang1 = blend(0, 2 * PI, aii);	
			
			vec3 norm = v.cross(vec3::UZ);			
			norm.rot(ang1 < PI ? ang1 : ang1 - PI, v);
			real alpha1 = blend2d(0.95, 1.05, aii, ai);
			
			vec3 pp = p + norm * (alpha1 * wid);
			real alpha2 = blend2d(0.8, 1, aii, ai);
			//int cor = blendcor(RGB(13, 15, 40), RGB(57, 58, 80), (0.5 + 0.5 * sin(alpha2 * aii * PI * 1)) * 1);	
			
			int cor = blendcor(RGB(200, 200, 200), RGB(58, 58, 58), (0.5 + 0.5 * sin(alpha2 * aii * PI * 1)) * 1);
			cor = blendcor(cor, RGB(58, 58, 58), real(i - lst) / step, 8);
			cor = blendcor(1, cor, (0.5 + 0.5 * norm.x) * rrnd(0, 2));
			
			pixelAA(pp, cor);
		}	
	}
}

// ------------------------------------------------
void curvejing(vec3 p, vec3 v, int len, real wid0)
{	
	vec3 cp1[] = {
		v,
		v.rotcopy(rrnd(-PI, PI), vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)).normcopy()),
		v.rotcopy(rrnd(-PI, PI), vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)).normcopy()),
		v.rotcopy(rrnd(-PI, PI), vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)).normcopy())
	};	
	
	v.norm();
	vec3 v0 = v;
	vec3 aa;
	int lst = 250;
	int step = rrnd(0.1, 0.12) * len;	
	for(int i = 0; i < len; i ++)
	{
		real ai = real(i) / len;
		
		real wid = i < len * 0.9 ? blend(wid0, wid0 * 0.25, i / (len * 0.9), 0.5) :  blend(wid0 * 0.25, wid0 * 1, (i - len * 0.9) / (0.1 * len), 8);
		
		if(i == len - 1)
		{					
			curvedleaf(p, v);
			step = rrnd(0.1, 0.12) * len;
			lst = i;
		}		
		if(rrnd(0, 15000) < 50)
		{
			int leaflen = 1000 * 0.5;
			//v = vec3(0.5, 0.5, 0.25);v.norm();
			vec3 zz = vec3::UZ;	
			
			{		
				//drawputao8(p + vec3(rrnd(-1, 1),rrnd(-1, 1), rrnd(-1, 1)).normcopy() * 0.01, vec3(rrnd(-1, 1),rrnd(-1, 1), rrnd(-1, 1)));
			}
		}
		//if(depth > 0)
		{	
			//v = bezier2(cp1, ai);
			v = v + vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)) * 0.01;
			v = v + vec3(0, -0.002, 0);
		}
		v.norm();
		p = p + v * 0.00046;
		for(int ii = 0; ii < 400; ii ++)
		{
			real aii = real(ii) / 400;
			real deta1 = blend2d(0, 1, ai, aii);
			real ang1 = blend(0, 2 * PI, aii);	
			
			vec3 norm = v.cross(vec3::UZ);			
			norm.rot(ang1 < PI ? ang1 : ang1 - PI, v);
			real alpha1 = blend2d(0.95, 1.05, aii, ai);
			
			vec3 pp = p + norm * (alpha1 * wid);
			real alpha2 = blend2d(0.8, 1, aii, ai);
			//int cor = blendcor(RGB(13, 15, 40), RGB(57, 58, 80), (0.5 + 0.5 * sin(alpha2 * aii * PI * 1)) * 1);	
			
			int cor = blendcor(RGB(200, 200, 200), RGB(58, 58, 58), (0.5 + 0.5 * sin(alpha2 * aii * PI * 1)) * 1);
			cor = blendcor(cor, RGB(58, 58, 58), real(i - lst) / step, 8);
			cor = blendcor(1, cor, (0.5 + 0.5 * norm.x) * rrnd(0, 2));
			
			pixelAA(pp, cor);
		}	
	}
}

// ------------------------------------------------
void teng(vec3 p, vec3 v, int len, real wid0, int depth = 0)
{
	vec3 cp1[] = {
		v,
		v.rotcopy(rrnd(-PI, PI) * 2.5, vec3(rrnd(-1, 1), rrnd(1, 1), rrnd(-1, 1)).normcopy()),
		v.rotcopy(rrnd(-PI, PI) * 2.5, vec3(rrnd(-1, 1), rrnd(1, 1), rrnd(-1, 1)).normcopy()),
		v.rotcopy(rrnd(-PI, PI) * 2.5, vec3(rrnd(-1, 1), rrnd(1, 1), rrnd(-1, 1)).normcopy())
	};	
	
	v.norm();
	vec3 v0 = v;
	vec3 aa;
	int last = 0;
	int step = rrnd(0.1, 0.15) * len;	
	vec3 gv0;
	real ang = 0;
	vec3 vg;
	vec3 vrot(rrnd(-1, 1), 1, rrnd(-1, 1));
	vrot.norm();
	for(int i = 0; i < len; i ++)
	{
		real ai = real(i) / len;
		
		real wid = blend(wid0, 0, ai, 0.8);
		
		if(i - last > step || i == len - 1)
		{
			step = rrnd(0.01, 0.1) * len;	
			last = i;
			
			if(rrnd(0, 1000) < (depth == 0 ? 120 : 100) && wid > 0.0001)
			{
				real deta1 = rrnd(-1, 1) * PI * (depth == 0 ? 0.1 : 0.3);	
				{
					teng(p, v.rotcopy(deta1), rrnd(0.1, 0.8) * len / 2, depth == 0 ? rrnd(0.08, 0.45) * wid : rrnd(0.5, 0.7) * wid, depth + 1);		
					//v.rot(-deta1 * rrnd(0.1, 0.5) * PI * 0.8);		
				}			
			}				
			if(rrnd(0, 1500) < 400)
			{
				//leavejing(p, vec3(rrnd(-0.1, 0.1), rrnd(-0.1, 0.1), rrnd(-0.1, 0.1)) + v, rrnd(150, 250), 0.001);
			}
			if(rrnd(0, 1500) < 150)
			{
				int leaflen = len * 0.5;
				//v = vec3(0.5, 0.5, 0.25);v.norm();
				vec3 zz = vec3::UZ;
				{		
					drawputao8(p + vec3(rrnd(-1, 1),rrnd(-1, 1), rrnd(-1, 1)).normcopy() * 0.01, vec3(rrnd(-1, 1),rrnd(-1, 1), rrnd(-1, 1)));
				}
			}
			if(rrnd(0, 200) < 16)
			{
				int cnt = 1 + rrnd(0, 3);
				for(int i = 0; i < cnt; i ++)
					curvedleaf(p, vec3(rrnd(-0.1, 0.1), rrnd(-0.1, 0.1), rrnd(-0.1, 0.1)) + v);
			}	
		}
		
		//if(p.y > 0)
		//	gv0 = gv0 + vec3(0, -0.0008, 0);
		//else
		//	gv0 = blend(gv0, vec3(0, 0, 0), 0.5);
		
		if(depth > 0)
		{	
			v = bezier2(cp1, ai);
			vec3 norm = v.cross(vec3::UX);	
			norm.norm();
			if(depth > 0)
				ang = blend(0, 4 * PI, ai, blend2d(4, 8, ai, p.y));
			else
				ang = blend(0, 4 * PI, ai, blend2d(14, 18, ai, p.y));
				
			v.rot(depth == 0 ? ang / 4 : ang, norm);
			
			vg = vg + vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)) * 0.01;
			
			v = v + vg;
		}
		else
		{
			ang = blend(.008 * PI, .0001 * PI, ai, blend2d(8, 2, ai, p.y));
			//v = v + vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)) * 0.01;
			vrot = vrot + vec3(rrnd(-1, 1), rrnd(-1, 1), rrnd(-1, 1)) * 0.01;
			vrot.norm();
			v.rot(ang, vrot);
		}
		
		v.norm();
		p = p + v * 0.0005;
		for(int ii = 0; ii < 400; ii ++)
		{
			real aii = real(ii) / 400;
			real deta1 = blend2d(0, 1, ai, aii);
			real ang1 = blend(0, 2 * PI, aii);	
			
			vec3 norm = v.cross(vec3::UZ);			
			norm.rot(ang1 < PI ? ang1 : ang1 - PI, v);
			real alpha1 = blend2d(0.95, 1.05, aii, ai);
			
			vec3 pp = p + norm * (alpha1 * wid);
			real alpha2 = blend2d(0.1, 10, aii, ai);
			//int cor = blendcor(RGB(13, 15, 40), RGB(57, 58, 80), (0.5 + 0.5 * sin(alpha2 * aii * PI * 1)) * 1);	
			if(depth > 0)
			{
				int cor = blendcor(RGB(200, 200, 200), RGB(58, 58, 58), (0.5 + 0.5 * sin(alpha2 * aii * PI * 1)) * 1);
				cor = blendcor(cor, RGB(58, 58, 58), real(i - last) / step, 8);
				cor = blendcor(1, cor, (0.5 + 0.5 * norm.x) * rrnd(0, 2));
				
				pixelAA(pp, cor);
			}
			else
			{
				int cor = blendcor(RGB(50, 50, 50), RGB(18, 18, 18), (0.5 + 0.5 * sin(alpha2 * aii * PI * 1)) * 1);
				cor = blendcor(cor, RGB(18, 18, 18), real(i - last) / step, 8);
				cor = blendcor(1, cor, (0.5 + 0.5 * norm.x) * rrnd(0, 2));
				
				pixelAA(pp, cor);				
			}
		}
	}
}

// ------------------------------------------------
void draw()
{
	fheightmap = fheightmap0;	

	leaf2size = rrnd(0.5, 0.8);	
	
	teng(vec3(0.5, 0.0, 0.5), vec3(rrnd(-1, 1), 0.16, rrnd(-1, 1)), rrnd(5000, 6000), rrnd(0.0025, 0.005));	
	
	teng(vec3(0.5, 0.0, 0.5), vec3(rrnd(-1, 1), 0.16, rrnd(-1, 1)), rrnd(3000, 4000), rrnd(0.0025, 0.005));	
	
}
// ================================================
int form()
{
	perlinmap(fheightmap0, IMAGESCALE, &fmin1, &fmax1, 12);
	perlinmap(fheightmap1, IMAGESCALE, &fmin1, &fmax1, 2);	
	fheightmap = fheightmap0;
	
	init();
	bk();
	draw();
	return 0;
}