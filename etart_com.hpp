#include <map>
#include <array>
#include <functional>
#include <string>
#include <vector>
#include <stdio.h>
#include <time.h>
#include <string>
#include <stdlib.h> 
#include <stdarg.h>
#include <math.h>
#include <tgmath.h>
#include <stdio.h>
#include <algorithm>
#include <list>
#include <sstream>
#include <fstream>
#include <iostream >
#include <set>
#include <queue>
#include <regex>
#include <numeric> 
#include <stdexcept>
#include <any>
#include <filesystem>
#include <unordered_set>
#include <random>
#include <limits>
#include <mutex>

// **********************************************************************
// 编译警告关闭
#pragma warning(disable:4244)
#pragma warning(disable:4305)
#pragma warning(disable:4267)
#pragma warning(disable:4819)
#pragma warning(disable:4018)
#pragma warning(disable:4005)
#pragma warning(disable:4190)

// **********************************************************************
#define	IMAGESCALE			1024							// 画布尺寸
#define						PMDLL
#ifdef WINDOWS
EXPORT_API    HANDLE		hConsole;						// 标准输出句柄HANDLE		
#endif

#define	DEVICE_CALLABLE
// **********************************************************************
#include "C:/Users/18858/Documents/pmsys/pmsys.hpp"

// **********************************************************************
// globles
// **********************************************************************
real				fmin1								= 0;
real				fmax1								= 0;
real				fheightmap0[IMAGESCALE][IMAGESCALE]	= {0};
real				fheightmap1[IMAGESCALE][IMAGESCALE]	= {0};
real				(*fheightmap)[IMAGESCALE]			= fheightmap0;

int					drawmap0[IMAGESCALE][IMAGESCALE]	= {0};
real				depthmap[IMAGESCALE][IMAGESCALE]	= {0};

int					(*drawmap)[IMAGESCALE]				= drawmap0;

// **********************************************************************
// projection
// **********************************************************************
#define	CAM_DIS	4.0
inline vector2 viewprj(const vector3& p)
{
	real zc = CAM_DIS;
	real vpf = zc / (zc + p.z);	
	vpf < 0.0 ? vpf = 0.0 : 0;
	return vector2(
		0.5 + (p.x - 0.5) * vpf,
		0.4 + (p.y - 0.4) * vpf
		);		
}
// ------------------------------------------------
vec3 v2tov3 (const vec2& p)
{
	real zz = blend(0, CAM_DIS, (p.y) / 0.4);
	real zc = CAM_DIS;
	real vpf = zc / (zc + zz);	
	real x3 = (p.x - 0.5) * vpf + 0.5;	
	return vec3(x3, -0.05, zz);
}
// **********************************************************************
// common draw 
// **********************************************************************
inline void pixel(int ix, int iy , int color, real alpha = 1, real depth = 0)
{	
	if(ix < 0 || iy < 0 || ix >= IMAGESCALE || iy >= IMAGESCALE)
		return;
	if(depthmap[ix][iy] == 0 || depth <= depthmap[ix][iy])
	{
		if(alpha >= 1)
		{
			drawmap[ix][iy] = color;
			depthmap[ix][iy] = depth < 0 ? 0 : depth;	
		}
		else
		{
			drawmap[ix][iy] = blendcor(drawmap[ix][iy], color, alpha, 1);
			//depthmap[ix][iy] = depth < 0 ? 0 : depth;
		}
	}	
}

inline void pixel(real x, real y , int color, real alpha = 1, real depth = 0)
{
	pixel(int(x * IMAGESCALE), int(y * IMAGESCALE), color, alpha, depth);
}
inline void pixel(const vector2& p, int color, real alpha = 1, real depth = 0)
{		
	pixel(p.x, p.y, color, alpha, depth);
}
inline void pixel(const vector3& p, int color, real alpha = 1)
{		
	pixel(viewprj(p), color, alpha, p.z);
	
}
inline void pixelAA(real x, real y, int color, real depth = 0)
{
	x = x * IMAGESCALE;	
	y = y * IMAGESCALE;
	int ix = int(x);
	int iy = int(y);
	real fx = x - ix;
	real fy = y - iy;

	if (ix < 0 || iy < 0 || ix >= IMAGESCALE || iy >= IMAGESCALE)
		return;
	if(depthmap[ix][iy] == 0 || depth <= depthmap[ix][iy])
	{
		drawmap[ix][iy] = color;
		depthmap[ix][iy] = depth < 0 ? 0 : depth;
		
		if (fx > 0 && ix + 1 < IMAGESCALE)
			drawmap[ix + 1][iy] = blendcor(drawmap[ix + 1][iy], color, fx);
		if (fy > 0 && iy + 1 < IMAGESCALE)
			drawmap[ix][iy + 1] = blendcor(drawmap[ix][iy + 1], color, fy);
	}
}
inline void pixelAA(const vector2& p, int color, real depth = 0)
{
	pixelAA(p.x, p.y, color, depth);
}
inline void pixelAA(const vector3& p, int color)
{
	pixelAA(viewprj(p), color, p.z);
}
// **********************************************************************
// write BMP
// **********************************************************************
int bmp_write( unsigned char* image, int xSize, int ySize, char *fileName)  
{  
	unsigned char header[54] = {
		0x42, 0x4d, 0, 0, 0, 0, 0, 0, 0, 0,
		54, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0, 
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		0, 0, 0, 0
	};
          
	// sizeof(long) = 4  
	long fileSize = (long)xSize * (long)ySize * 3 +54 ;     // image array size + 54  
	long temp = 0;                      // temp for byte convertion  
	long width = xSize;                 // image width  
	long height = ySize;                // image height 
	int i = 0;                          // loop variable  

	FILE *fp;  

	for( i=0 ; i != 4 ; ++i)                // write fileSize from byte2 ~ byte 6  
	{  
		header[2+i] = (unsigned char)( fileSize >> (8*i) ) & 0x000000ff;        
	}  
	for( i=0 ; i != 4 ; ++i)                // write width from byte18 ~ byte 22  
	{  
		header[18+i] = (unsigned char)( width >> (8*i) ) & 0x000000ff;          
	}  
	for( i=0 ; i != 4 ; ++i)                // write height from byte22 ~ byte 26  
	{  
		header[22+i] = (unsigned char)( height >> (8*i) ) & 0x000000ff;         
	}

	if( ! ( fp = fopen(fileName,"wb")  ) )                      
			return -1;  

	fwrite(header, sizeof(unsigned char), 54, fp);                              // write header  
	fwrite(image, sizeof(unsigned char),(size_t)(long)xSize*ySize*3, fp);       // write image array  

	fclose(fp);  
	return 0;  
} 

// **********************************************************************
// main
// **********************************************************************
extern int form_image();
unsigned char buf[IMAGESCALE * IMAGESCALE * 3];
int main(int argc, char* argv[])
{
	int i, ii, jj;
	char path[260], imgname[512];	
	//getcwd(path, 260);
	
	srand(18858146);
	for( i = 0; i < 1; i ++)
	{
		// clear
		memset(drawmap0, 0, sizeof(drawmap0));
		memset(depthmap, 0, sizeof(depthmap));

		printf("generating image ... \n");
		form_image();
		printf ("done!\n");

		printf("saving ... ");
		sprintf(imgname, "C:/Users/18858/Documents/_ETART/output/images%d.bmp", i);
		for(ii = 0; ii < IMAGESCALE; ii ++)
		for(jj = 0; jj < IMAGESCALE; jj ++)
		{
			buf[3 * (jj * IMAGESCALE + ii) + 0] = GetBValue(drawmap[ii][jj]);
			buf[3 * (jj * IMAGESCALE + ii) + 1] = GetGValue(drawmap[ii][jj]);
			buf[3 * (jj * IMAGESCALE + ii) + 2] = GetRValue(drawmap[ii][jj]);
		}
		bmp_write(buf, IMAGESCALE, IMAGESCALE, imgname);

		printf("done!\n");

		//system(imgname);
	}
	//getchar();
	return 0;
}
// **********************************************************************
#define SZ 	4
#define SX 	(1.0 / IMAGESCALE)
#define SY 	(1.0 / IMAGESCALE)
#define FIMAGESCALE 	real(IMAGESCALE)
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

real blend2d(real a,real b, real ax, real ay)
{
	return blend(a, b, ax*ay);
}
real CircleBlend2(real a,real b, real ax, real ay)
{
	return blend(a, b, ax*ay);
}
real hash(vec2 v) {
    // 一个简单的哈希算法，将 x 和 y 分量组合
    real value = v.x * 12.9898 + v.y * 78.233;
    // 使用 std::sin 函数进一步打乱值
    value = std::sin(value) * 43758.5453;
    // 取小数部分作为哈希值
    return value - std::floor(value);
}