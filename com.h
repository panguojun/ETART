#include <stdio.h>  
#include <stdlib.h>  
#include <memory.h>
#include <math.h>
#include <sstream>
// **********************************************************************
// macro
// **********************************************************************
#define PRINT(msg)		{std::stringstream ss; ss << "\n> " << msg; printf(ss.str().c_str());}
#define MYTRACE(msg)	{std::stringstream ss; ss << "\n" << msg; printf(ss.str().c_str());}// ::OutputDebugString(ss.str().c_str());}

#define real			double
#define PI				3.1415926
#define	IMAGESCALE		2048
#define	FIMAGESCALE		((real)IMAGESCALE)
#define	PSZ				(1.0f / FIMAGESCALE)

#define RGB(r, g, b)	((int)(((unsigned char)(r)) | (((short)(unsigned char)(g)) << 8) | (((int)(unsigned char)(b)) << 16)))
#define GetRValue(rgb)	((unsigned char)(0x000000ff & (rgb)))
#define GetGValue(rgb)	((unsigned char)((0x0000ff00 & (rgb))>>8))
#define GetBValue(rgb)	((unsigned char)((0x00ff0000 & (rgb))>>16))

#define rrnd(min, max)	((max) == (min) ? (min) : (min) + ( (real)(rand()) / (RAND_MAX + 1.0)) * ((max) - (min)))

#define VEC2(x, y, v)	(v[0] = x, v[1] = y)
#define vec2			vector2
#define vec3			vector3
#define vec				vector3
// **********************************************************************
// 2D
// **********************************************************************
struct vector2{
	union{
		real val[2];
		struct{
			real x;
			real y;
		};
	};
	
	static const vector2 ZERO;
	static const vector2 UX;
	static const vector2 UY;
	static const vector2 CENTER;
	
	real& operator [](int ind){
		return val[ind];
	}

	vector2(){
		x = 0;
		y = 0;
	}
	vector2(real _x, real _y){
		x = _x;
		y = _y;
	}
	vector2 operator + (const vector2& _p) const
	{
		vector2 fp;
		fp.x = x + _p.x;
		fp.y = y + _p.y;

		return fp;
	}
	vector2 operator - (const vector2& _p) const
	{
		vector2 fp;
		fp.x = x - _p.x;
		fp.y = y - _p.y;
		return fp;
	}
	vector2 operator - () const
	{
		vector2 fp;
		fp.x = - x;
		fp.y = - y;
		return fp;
	}

	vector2 operator * (real s) const
	{
		vector2 fp;
		fp.x = s * x;
		fp.y = s * y;
		return fp;		
	}
	vector2 operator * (const vector2& b) const
	{
		return vector2(x * b.x - y * b.y, x * b.y + y * b.x);	
	}
	vector2 operator / (real s) const
	{
		vector2 fp;
		fp.x = x / s;
		fp.y = y / s;
		return fp;		
	}
	real len() const
	{
		return sqrt(x * x + y * y);
	}
	real sqrlen() const
	{
		return (x * x + y * y);
	}
	real angle() const
	{
		return atan2(y, x);
	}
	void norm()
	{
		real r = len();
		if(r > 0)
		{
			x /= r;
			y /= r;
		}	
	}
	vector2 normcpy() const
	{
		real r = len();
		if(r > 0)
		{
			return vector2( x / r, y / r);
		}	
		return vector2::ZERO;
	}
	void rot(real angle)
	{
		(*this) = (*this) * vector2::fromanglelength(angle, 1);		
	}
	vector2 rotcpy(real angle) const
	{
		return (*this) * vector2::fromanglelength(angle, 1);		
	}
	void rot(real angle, const vector2& o)
	{
		vector2 v = (*this) - o;
		v = v * vector2::fromanglelength(angle, 1);
		(*this) = v + o;
	}
	vector2 rotcpy(real angle, const vector2& o) const
	{
		vector2 v = (*this) - o;
		v = v * vector2::fromanglelength(angle, 1);
		return v + o;	
	}	
	real dot(const vector2& v) const
	{
		return x * v.x + y * v.y;
	}
	real cross(const vector2& v) const
	{
		return x * v.y - y * v.x;
	}
	static vector2 fromanglelength(real _angle, real _r);
};

const vector2 vector2::ZERO = vector2(0, 0);
const vector2 vector2::UX = vector2(1, 0);
const vector2 vector2::UY = vector2(0, 1);
const vector2 vector2::CENTER = vector2(0.5f, 0.5f);

vector2 vector2::fromanglelength(real _angle, real _r)
{
	return vector2(_r * cos(_angle), _r * sin(_angle));
}
// **********************************************************************
// 3D
// **********************************************************************

struct vector3
{
	static const vector3 ZERO;
	static const vector3 UX;
	static const vector3 UY;
	static const vector3 UZ;
	static const vector3 CENTER;
	union{
		real val[3];
		struct
		{
			float x; 
			float y;
			float z;
		};
	};
	vector3()
	{
		x = 0;
		y = 0;
		z = 0;
	}
	vector3(float _x, float _y, float _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}
	real& operator [](int ind){
		return val[ind];
	}
	vector3 operator + (const vector3& _p) const
	{
		vector3 fp;
		fp.x = x + _p.x;
		fp.y = y + _p.y;
		fp.z = z + _p.z;
		return fp;
	}
	vector3 operator - (const vector3& _p) const
	{
		vector3 fp;
		fp.x = x - _p.x;
		fp.y = y - _p.y;
		fp.z = z - _p.z;
		return fp;
	}
	vector3 operator - () const
	{
		vector3 fp;
		fp.x = - x;
		fp.y = - y;
		fp.z = - z;
		return fp;
	}
	vector3 operator * (float s) const
	{
		vector3 fp;
		fp.x = s * x;
		fp.y = s * y;
		fp.z = s * z;
		return fp;		
	}
	vector3 operator / (float s) const
	{
		vector3 fp;
		fp.x = x / s;
		fp.y = y / s;
		fp.z = z / s;
		return fp;		
	}
	float len() const
	{
		return sqrt(x * x + y * y + z * z);
	}
	float sqrlen() const
	{
		return (x * x + y * y + z * z);
	}	
	void norm()
	{
		float r = len();
		if(r > 0)
		{
			x /= r;
			y /= r;
			z /= r;
		}	
	}
	vector3 normcopy()
	{
		float r = len();
		if(r > 0)
		{
			return vector3( this->x / r,
							this->y / r,
							this->z / r);
		}	
		return vector3(0, 0, 0);
	}
	void fromanglelength(float _angle, float _r)
	{
		x = _r * cos(_angle);
		y = _r * sin(_angle);
	}
	float angle() const
	{
		return atan2(y, x);
	}
	float dot(const vector3& v) const
	{
		return x * v.x + y * v.y + z * v.z;
	}

	vector3 cross(const vector3& v) const
	{
		vector3 n;
		n.x = - (y * v.z - z * v.y);
		n.y = - (z * v.x - x * v.z);
		n.z = - (x * v.y - y * v.x);
		return n;
	}
	void rot(float angle, const vector3& ax = vector3::UZ);
	vector3 rotcopy(float angle, const vector3& ax = vector3::UZ) const;

	static vector3 rnd(float min = 0, float max = 1);

	static vector3 rndrad(float r = 1);
};

const vector3 vector3::ZERO = vector3(0, 0, 0);
const vector3 vector3::UX = vector3(1, 0, 0);
const vector3 vector3::UY = vector3(0, 1, 0);
const vector3 vector3::UZ = vector3(0, 0, 1);
const vector3 vector3::CENTER = vector3(IMAGESCALE / 2, IMAGESCALE / 2, IMAGESCALE / 2);

// ----------------------------------------
vector3 vector3::rnd(float min, float max)
{
	return vector3( rrnd(min, max), rrnd(min, max), rrnd(min, max) );
}
// ----------------------------------------
vector3 vector3::rndrad(float r)
{
	return rnd(-r, r).normcopy();
}

// **********************************************************************
// 四元数 
// **********************************************************************
struct  quaternion
{
	float w, x, y, z;

    quaternion (
        float fW = 1.0,
        float fX = 0.0, float fY = 0.0, float fZ = 0.0)
	{
		w = fW;
		x = fX;
		y = fY;
		z = fZ;
	}
    quaternion (const quaternion& rkQ)
	{
		w = rkQ.w;
		x = rkQ.x;
		y = rkQ.y;
		z = rkQ.z;
	}
	quaternion(const float& rfAngle, const vector3& rkAxis)
    {
        this->FromAngleAxis(rfAngle, rkAxis);
    }
    //-----------------------------------------------------------------------
    vector3 operator* (const vector3& v) const
    {
		// nVidia SDK implementation
		vector3 uv, uuv; 
		vector3 qvec(x, y, z);
		uv = qvec.cross(v); 
		uuv = qvec.cross(uv); 
		uv = uv * (2.0f * w); 
		uuv = uuv * 2.0f; 
		
		return v + uv + uuv;
    }
	//-----------------------------------------------------------------------
    void FromAngleAxis (float rfAngle,
        const vector3& rkAxis)
    {
        // assert:  axis[] is unit length
        //
        // The quaternion representing the rotation is
        //   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

        float fHalfAngle ( 0.5*rfAngle );
        float fSin = sin(fHalfAngle);
        w = cos(fHalfAngle);
        x = fSin*rkAxis.x;
        y = fSin*rkAxis.y;
        z = fSin*rkAxis.z;
    }
};
// ----------------------------------------
void vector3::rot(float angle, const vector3& ax)
{
	quaternion q(angle, ax);
	*this = q * (*this);
}
// ----------------------------------------
vector3 vector3::rotcopy(float angle, const vector3& ax) const
{
	quaternion q(angle, ax);
	return q * (*this);
}
// **********************************************************************
// ray3d
// **********************************************************************
struct ray3
{
	vector3 o;
	vector3 dir;
	ray3(){}
	ray3(const vector3& p1, const vector3& p2)
	{
		o = p1;
		dir = (p2 - p1).normcopy();
	}
	
};
// **********************************************************************
// segment
// **********************************************************************
struct segment
{
	vector3 s;
	vector3 e;
	segment(const vector3& _s, const vector3& _e)
	{
		s = _s;
		e = _e;
	}
};
// **********************************************************************
// plane
// **********************************************************************
struct plane
{
	vector3 o;
	vector3 n;
	plane(const vector3& _o, const vector3& _n)
	{
		o = _o;
		n = _n;
		n.norm();
	}
	plane(const vector3& p1, const vector3& p2, const vector3& p3)
	{
		n = (p2 - p1).cross(p3 - p1);
		n.norm();
		o = p1;
	}
};

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
// Perlin 2d random
// **********************************************************************
// Ken Perlin's Map Gennor
// ------------------------------------------------
int myfloor(real value){return (value >= 0 ? (int)value : (int)value - 1); } // 取整
//inline real fade(real x){	return (x*x*x*(x*(6*x -15) + 10)); } // this equates to 6x^5 - 15x^4 + 10x^3  // 插值算法
real fade(real x)
{	
	if(x < 0.5)
		x = pow(x * 2.0f, 2);
	else
		x = 1 - pow(1 - (x - 0.5f) * 2, 2) + 1;
	return x / 2;
}
// ----------------------------------------
inline real dot(real x1, real y1, real x2, real y2)
{
	return x1 * x2 + y1 * y2;
}
// ------------------------------------------------
extern real blend(real h1, real h2, real alpha, real power);
void perlinmap(real map[IMAGESCALE][IMAGESCALE], int size, real* min1, real* max1, int octaves)
{
	//set up some variables
	int i,j,k,x,y,grad11,grad12,grad21,grad22;

	real pixel_value, 
		  fracX, fracY,
		  noise11, noise12, noise21, noise22,
		  interpolatedx1, interpolatedx2, interpolatedxy,
		  amplitude, frequency, 
		  gain = 0.50f, lacunarity = 2.0f;
	real gradients[8][2];
	int permutations[256];

	*min1 = 1000000;
	*max1 = -1000000;

	//梯度	
	for (i = 0; i < 8; ++i)
	{
		//gradients[i][0] = cos(0.785398163f * (real)i);
		//gradients[i][1] = sin(0.785398163f * (real)i);
		real rr = 1;//(i % 2 == 1 ? 1 : 0.01);
		gradients[i][0] = rr * cos(0.785398163f * (real)i);// PI / 4
		gradients[i][1] = rr * sin(0.785398163f * (real)i);
	}

	// 置换表	
	{
		for (i = 0; i < 256; ++ i)
			permutations[i] = i;
		for (i = 0; i < 256; ++ i)
		{
			j = rand() % 256;
			k = permutations[i];
			permutations[i] = permutations[j];
			permutations[j] = k;
		}
	}	
	
	// 生成MAP
	for (i = 0;	i < size; ++i)
	{
		for (j = 0;	j < size; ++j)
		{
			//get the value for this pixel by adding successive layers
			amplitude = 1.0f;
			frequency = 1.0f / (real)size;
			pixel_value = 0.0f;
			map[j][i] = 0;
			for (k = 0; k < octaves; ++k)
			{
				// 整数
				x = (int)((real)j * frequency);
				y = (int)((real)i * frequency);

				// 小数
				fracX = (real)j * frequency - (real)x;
				fracY = (real)i * frequency - (real)y;

				// 得到梯度索引
				grad11 = permutations[(x + permutations[y % 255]) % 255]			% 8;
				grad12 = permutations[(x + 1 + permutations[y % 255]) % 255]		% 8;
				grad21 = permutations[(x + permutations[(y + 1) % 255]) % 255]		% 8;
				grad22 = permutations[(x + 1 + permutations[(y + 1) % 255]) % 255]	% 8;

				// 四个角的梯度投影
				noise11 = dot(gradients[grad11][0], gradients[grad11][1], fracX, fracY);
				noise12 = dot(gradients[grad12][0], gradients[grad12][1], fracX - 1.0f, fracY);
				noise21 = dot(gradients[grad21][0], gradients[grad21][1], fracX, fracY - 1.0f);
				noise22 = dot(gradients[grad22][0], gradients[grad22][1], fracX - 1.0f, fracY - 1.0f);

				// 插值算法
				fracX = fade(fracX);
				fracY = fade(fracY);
				
				interpolatedx1 = blend(noise11, noise12, fracX, 0.5);
				interpolatedx2 = blend(noise21, noise22, fracX, 0.5);
				interpolatedxy = blend(interpolatedx1, interpolatedx2, fracY, 0.5);

				//	叠加
				pixel_value += interpolatedxy * amplitude;
				amplitude *= gain;

				// 缩小区域
				frequency *= lacunarity;
			}

			//put it in the map
			map[j][i] = pixel_value;

			//do some quick checks
			if (pixel_value < *min1)
				*min1 = pixel_value;
			else if (pixel_value > *max1)
				*max1 = pixel_value;
		}
	}
}
// **********************************************************************
// shape blend
// **********************************************************************
real blend(real h1, real h2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	if(power != 1.0)
		alpha = pow(alpha, power);
	
	return h1 * (1 - alpha) + h2 * alpha;
}
vector2 blend(vector2 v1, vector2 v2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	if(power != 1.0)
		alpha = pow(alpha, power);
	
	return v1 * (1 - alpha) + v2 * alpha;
}
vector3 blend(vector3 v1, vector3 v2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	if(power != 1.0)
		alpha = pow(alpha, power);
	
	return v1 * (1 - alpha) + v2 * alpha;
}

// ------------------------------------------------
real blend2(real h1, real h2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	
	if(alpha < 0.5f)
	{
		alpha = (0.5f - alpha) * 2.0f;
		if(power != 1.0)
			alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;
	}
	else
	{
		alpha = (alpha - 0.5f) * 2.0f;
		if(power != 1.0)
			alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;	
	}
}
// ------------------------------------------------
vector2 blend2(vector2 h1, vector2 h2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	
	if(alpha < 0.5f)
	{
		alpha = (0.5f - alpha) * 2.0f;
		if(power != 1.0)
			alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;
	}
	else
	{
		alpha = (alpha - 0.5f) * 2.0f;
		if(power != 1.0)
			alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;	
	}
}
// ------------------------------------------------
vector3 blend2(vector3 h1, vector3 h2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	
	if(alpha < 0.5f)
	{
		alpha = (0.5f - alpha) * 2.0f;
		if(power != 1.0)
			alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;
	}
	else
	{
		alpha = (alpha - 0.5f) * 2.0f;
		if(power != 1.0)
			alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;	
	}
}

float blend3(float h1, float h2, float alpha, float mid = 0.5, float power = 1)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	
	if(alpha < mid)
	{
		alpha = (mid - alpha) / mid;
		alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;	}
	else
	{
		alpha = (alpha - mid) / (1 - mid);
		alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;	
	}
}

float blend4(float h1, float h2, float alpha, float power = 1)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	
	float mid = 0.5;
	if(alpha < mid)
	{
		alpha = alpha / mid;
		alpha = pow(alpha, 1.0f / power);
		return h1 * (1 - alpha) + h2 * alpha;
	}
	else
	{
		alpha = (1 - alpha) / (1 - mid);
		alpha = pow(alpha, 1.0f / power);
		return h1 * (1 - alpha) + h2 * alpha;	
	}
}

float blend5(float h1, float h2, float alpha, float mid = 0.5, float power = 1)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	
	if(alpha < mid)
	{
		alpha = alpha / mid;
		alpha = pow(alpha, 1.0f / power);
		return h1 * (1 - alpha) + h2 * alpha;
	}
	else
	{
		alpha = (1 - alpha) / (1 - mid);
		alpha = pow(alpha, 1.0f / power);
		return h1 * (1 - alpha) + h2 * alpha;	
	}
}

// ------------------------------------------------------------------------------------------------
vec3 bezier2(vec3 cp[3], real t)
{
	return cp[0] * ((1 - t) * (1 - t)) + cp[1] * (2 * t * (1 - t)) + cp[2] * (t * t);
}
// ------------------------------------------------------------------------------------------------
vec3 bezier3(vec3 cp[4], real t)
{
	real s = 1 - t;
	return  cp[0] * (s * s * s) 	+ 
			cp[1] * (3 * t * s * s) +
			cp[2] * (3 * t * t * s) +
			cp[3] * (t * t * t);
}
// -----------------------------------------------------------------------
inline float CircleBlend(float h1, float h2, float alpha)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	alpha = sqrt(1 - (1-alpha)*(1-alpha));
	
	return h1 * (1 - alpha) + h2 * alpha;
}
// -----------------------------------------------------------------------
inline float CircleBlend2(float h1, float h2, float alpha, float power = 1)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	
	alpha = pow(alpha, power);

	if(alpha < 0.5)
	{
		alpha = alpha * 2;
		alpha = sqrt(1 - (1-alpha)*(1-alpha));
	
		return h1 * (1 - alpha) + h2 * alpha;
	}
	else
	{
		alpha = (alpha - 0.5) * 2;
		alpha = sqrt(1 - (alpha)*(alpha));

		return h1 * (1 - alpha) + h2 * alpha;
	}

}

// -----------------------------------------------------------------------
// 三角函数插值
inline float BlendSin(float h1, float h2, float alpha)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;

	alpha = sin(alpha * PI / 2);

	return h1 * (1 - alpha) + h2 * alpha;
}
// -----------------------------------------------------------------------
// 傅立叶级数
inline float FT(float angle, float t[] = 0, float dt = 0)
{
	if(t == 0)
	{
		static float s_t0[] = {rrnd(0, 1), rrnd(0, 1), rrnd(0, 1), rrnd(0, 1)};
		t = s_t0;
	}

	float yy = 0;
	yy += 1		* sin(1 * angle + (t[0] + dt) * PI);
	yy += 0.5	* sin(2 * angle + (t[1] + dt) * PI);
	yy += 0.25	* sin(4 * angle + (t[2] + dt) * PI);
	yy += 0.125	* sin(8 * angle + (t[3] + dt) * PI);

	return yy;
}
// -----------------------------------------------------------------------
inline float FTU(float ang, float t[] = 0, float dt = 0)
{
	float ft = FT(ang, t, dt);
	float max1 = (1 + 0.5 + 0.25 + 0.125);
	float min1 = - max1;
	return (ft - min1) / (max1 - min1);
}
// -----------------------------------------------------------------------
inline float BlendFT(float h1, float h2, float alpha, float t[] = 0, float dt = 0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;

	alpha = FTU(alpha, t, dt);

	return h1 * (1 - alpha) + h2 * alpha;
}

// ------------------------------------------------
// 2d
real blend2d(real h1, real h2, real alphaX, real alphaY)
{
	int size;
	real alpha;
	alphaX < 0 ? alphaX = 0 : 0;
	alphaX > 1 ? alphaX = 1 : 0;
	alphaY < 0 ? alphaY = 0 : 0;
	alphaY > 1 ? alphaY = 1 : 0;
	size = IMAGESCALE;
	alpha = fheightmap[(int)(alphaX * size) % size][(int)(alphaY * size) % size];
	alpha = (alpha - fmin1) / (fmax1 - fmin1);

	return h1 * (1 - alpha) + h2 * alpha;
}

// **********************************************************************
// color blend
// **********************************************************************
int blendcor(int c1, int c2, real alpha, real power = 1.0)
{
	//alpha = alpha != 1 ? abs(alpha) - (int)(alpha) : 1;	
	alpha = alpha > 1 ? 1 : alpha;
	alpha = alpha < 0 ? 0 : alpha;
	
	if(power != 1.0)
		alpha = pow(alpha, power);
	 return RGB( GetRValue(c2) * alpha + GetRValue(c1) * (1 - alpha),
				 GetGValue(c2) * alpha + GetGValue(c1) * (1 - alpha),
				 GetBValue(c2) * alpha + GetBValue(c1) * (1 - alpha)
				 );
	
}

// ------------------------------------------------
// 2d
int blendcor2d(int c1, int c2, real alphaX, real alphaY, real power = 1.0)
{
	real alpha;
	int size = IMAGESCALE;
	alphaX < 0 ? alphaX = 0 : 0;
	alphaX > 1 ? alphaX = 1 : 0;
	alphaY < 0 ? alphaY = 0 : 0;
	alphaY > 1 ? alphaY = 1 : 0;
	
	alpha = fheightmap[(int)(alphaX * size) % size][(int)(alphaY * size) % size];
	alpha = (alpha - fmin1) / (fmax1 - fmin1);

	return blendcor(c1, c2, alpha, power);
}

// **********************************************************************
// collision
// **********************************************************************
// A = y2 - y1,B = x1- x2,C = x2 * y1 - x1 * y2
//distance(P,l)=|aX+bY+c|/   
inline real distance(const vector3& p1, const vector3& p2)
{
	return (p1 - p2).len();
}
//点到线段的最短距离,x0,y0是圆心
real distance(const vector3& p, const vector3& p1,const vector3& p2)
{
    double ans = 0;
    double a, b, c;
    a = distance(p1, p2);
    b = distance(p1, p);
    c = distance(p2, p);
    if (c+b==a) {//点在线段上
      ans = 0;
      return ans;
    }
    if (a<=0.00001) {//不是线段，是一个点
      ans = b;
      return ans;
    }
    if (c*c >= a*a + b*b) { //组成直角三角形或钝角三角形，p1为直角或钝角
      ans = b;
      return ans;
    }
    if (b * b >= a * a + c * c) {// 组成直角三角形或钝角三角形，p2为直角或钝角
      ans = c;
      return ans;
    }
    // 组成锐角三角形，则求三角形的高
    double p0 = (a + b + c) / 2;// 半周长
    double s = sqrt(p0 * (p0 - a) * (p0 - b) * (p0 - c));// 海伦公式求面积
    ans = 2*s / a;// 返回点到线的距离（利用三角形面积公式求高）
    return ans;
}

real distance(const vector3& p, const segment& l)
{	
	return distance(p, l.s, l.e);
}
real distance(const vector3& p1s, const vector3& p1e,const vector3& p2s, const vector3& p2e)
{	
	segment l1(p1s, p1e);
	segment l2(p2s, p2e);
	
	return std::min(
		std::min( distance(p1s,l2),distance(p1e,l2) ),	
		std::min( distance(p2s,l1),distance(p2e,l1) )
	);
} 
// intersect3_segmentplane(): find the 3D intersection of a segment and a plane
//    Input:  S = a segment, and Pn = a plane = {Point V0;  Vector n;}
//    Output: *I0 = the intersect point (when it exists)
//    Return: 0 = disjoint (no intersection)
//            1 =  intersection in the unique point *I0
//            2 = the  segment lies in the plane
int	intersect3_segmentplane( const segment& S, const plane&  Pn, vector3& I )
{
    vector3    u = S.e - S.s;
    vector3    w = S.s - Pn.o;

    float     D = Pn.n.dot(u);
    float     N = -Pn.n.dot(w);

    if (fabs(D) < 0.0001) {           	 // segment is parallel to plane
        if (N == 0)                      // segment lies in plane
            return 2;
        else
            return 0;                    // no intersection
    }
    // they are not parallel
    // compute intersect param
    float sI = N / D;
    if (sI < 0 || sI > 1)
        return 0;                        // no intersection

    I = S.s + u * sI;                  // compute segment intersect point
    return 1;
}
int	intersect3_segmentplaneex( const segment& S, const plane&  Pn, real& dis )
{
    vector3    u = S.e - S.s;
    vector3    w = S.s - Pn.o;

    float     D = Pn.n.dot(u);
    float     N = -Pn.n.dot(w);

    if (fabs(D) < 0.0001) {           	 // segment is parallel to plane
        if (N == 0)                      // segment lies in plane
		{
			dis = 0;
            return 2;
		}
        else
            return 0;                    // no intersection
    }
    // they are not parallel
    // compute intersect param
    float sI = N / D;
    if (sI < 0 || sI > 1)
        return 0;                        // no intersection

    dis = (u * sI).len();                  // compute segment intersect point
    return 1;
}

// -----------------------------------------------------------------
// 三角形碰撞 
bool  checkPointInTriangle(
	const vector3& point, 
	const vector3& a, const vector3& b, const vector3& c,
	float tolerance = 0.005f) 
{		
	real total_angles = 0.0f;	

	// make the 3 vectors

	vector3 v1(point.x - a.x, point.y - a.y, point.z - a.z);
	v1.norm();
	vector3 v2(point.x - b.x, point.y - b.y, point.z - b.z);
	v2.norm();
	vector3 v3(point.x - c.x, point.y - c.y, point.z - c.z);	
	v3.norm();
	
	real Dot1 = v2.dot(v1);
	Dot1 < -1.0f ?	 Dot1 = -1.0f : NULL;
	Dot1 >  1.0f ?   Dot1 =  1.0f : NULL;
	total_angles += acos(Dot1); 

	real Dot2 = v3.dot(v2);
	Dot2 < -1.0f ? 	Dot2 = -1.0f : NULL;
	Dot2 >  1.0f ?  Dot2 =  1.0f : NULL;
	total_angles += acos(Dot2);

	real Dot3 = v1.dot(v3);
	Dot3 < -1.0f ?	Dot3 = -1.0f : NULL;
	Dot3 >  1.0f ?	Dot3 =  1.0f : NULL;
	total_angles += acos(Dot3); 
			
	if (fabs(total_angles - 2.0f * PI) <= tolerance)			
		return true;

	return false;
}

// -----------------------------------------------------------------
bool  intersectTrangle(
			const ray3& ray, const vector3& pt1, const vector3& pt2, const vector3& pt3, vector3& hit)
{	
	if(intersect3_segmentplane(segment(ray.o, ray.o + ray.dir * 100), plane(pt1,pt2,pt3), hit))
		return checkPointInTriangle(hit, pt1, pt2, pt3);		
	return 0;
}

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
extern int form();
unsigned char buf[IMAGESCALE * IMAGESCALE * 3];
int main(int argc, char* argv[])
{
	int i, ii, jj;
	char path[260], imgname[512];	
	//getcwd(path, 260);
	
	srand(18858146);
	for( i = 0; i < 100; i ++)
	{
		// clear
		memset(drawmap0, 0, sizeof(drawmap0));
		memset(depthmap, 0, sizeof(depthmap));

		printf("generating image ... \n");
		form();
		printf ("done!\n");

		printf("saving ... ");
		sprintf(imgname, "output/images%d.bmp", i);
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
