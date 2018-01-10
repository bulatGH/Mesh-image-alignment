#ifndef BASICIMAGES_H
#define BASICIMAGES_H

#include "AllTypes.h"

class Point3D_F;
class Point4D_F;
class Images2D_S;
class Images2D_F;
class Images3D_S;
class Images3D_F;
class sImages2D_S;
class sImages3D_S;
class GradientImage_2D;
class GradientImage_3D;

template <class T>
class Point2D
{
public:

    Point2D() { };

    ~Point2D() { };

    Point2D(T tx, T ty) :x(tx), y(ty) {};

    __declspec ( property ( put = SetX, get = GetX ) ) 
                T X ;

    __declspec ( property ( put = SetY, get = GetY ) ) 
                T Y ;

    void SetX(T tx) {x = tx;};

    void SetY(T ty) {y = ty;};

    T GetX() { return x; };

    T GetY() { return y; };

private:

    T x;
    T y;

};

template <class T>
class Point3D
{
public:

    Point3D() { };

    ~Point3D() { };

    Point3D(T tx, T ty, T tz) :x(tx), y(ty), z(tz) { };

    __declspec ( property ( put = SetX, get = GetX ) ) 
                T X ;

    __declspec ( property ( put = SetY, get = GetY ) ) 
                T Y ;

    __declspec ( property ( put = SetZ, get = GetZ ) ) 
                T Z ;

    void SetX(T tx) {x = tx;};

    void SetY(T ty) {y = ty;};

    void SetZ(T tz) {z = tz;};

    T GetX() { return x; };

    T GetY() { return y; };

    T GetZ() { return z; };

private:

    T x;
    T y;
    T z;

};

class Point3D_F
{
public:

    Point3D_F();

    ~Point3D_F() { };

    Point3D_F(float tx, float ty, float tz) :x(tx), y(ty), z(tz) { };

    __declspec ( property ( put = SetX, get = GetX ) ) 
                float X ;

    __declspec ( property ( put = SetY, get = GetY ) ) 
                float Y ;

    __declspec ( property ( put = SetZ, get = GetZ ) ) 
                float Z ;

    void SetX(float tx) {x = tx;};

    void SetY(float ty) {y = ty;};

    void SetZ(float tz) {z = tz;};

    float GetX() { return x; };

    float GetY() { return y; };

    float GetZ() { return z; };

    Point3D_F operator+(const Point3D_F& tpoint);
    Point3D_F operator-(const Point3D_F& tpoint);
    Point3D_F operator*(const float tcoef);
    Point3D_F operator/(const float tcoef);
    void operator+=(const Point3D_F& tpoint);
    void operator-=(const Point3D_F& tpoint);
    void operator*=(const float tcoef);
    void operator/=(const float tcoef);
    void operator=(const Point3D_F& tpoint);
    bool operator==(const Point3D_F& tpoint);

    static double Distance(const Point3D_F& tpoint1, const Point3D_F& tpoint2);

    static Point3D_F Transform(Point3D_F tpoint1, double** tmatr);

    static double** RotationYawPitchRoll(float dx, float dy, float dz);

    static double ComputeDeterminant3_3(double** tmatr);

    static double Dot(Point3D_F& tpoint1, Point3D_F& tpoint2);

private:

    float x;
    float y;
    float z;
};

class Point4D_F
{
public:

    Point4D_F();

    ~Point4D_F() { };

    Point4D_F(float tx, float ty, float tz, float tw) :x(tx), y(ty), z(tz), w(tw) { };

    __declspec ( property ( put = SetX, get = GetX ) ) 
                float X ;

    __declspec ( property ( put = SetY, get = GetY ) ) 
                float Y ;

    __declspec ( property ( put = SetZ, get = GetZ ) ) 
                float Z ;

    __declspec ( property ( put = SetW, get = GetW ) ) 
                float W ;

    void SetX(float tx) {x = tx;};

    void SetY(float ty) {y = ty;};

    void SetZ(float tz) {z = tz;};

    void SetW(float tw) {w = tw;};

    float GetX() { return x; };

    float GetY() { return y; };

    float GetZ() { return z; };

    float GetW() { return w; };

    Point4D_F operator+(const Point4D_F& tpoint);
    Point4D_F operator-(const Point4D_F& tpoint);
    Point4D_F operator*(const float tcoef);
    Point4D_F operator/(const float tcoef);
    void operator+=(const Point4D_F& tpoint);
    void operator-=(const Point4D_F& tpoint);
    void operator*=(const float tcoef);
    void operator/=(const float tcoef);
    void operator=(const Point4D_F& tpoint);
    bool operator==(const Point4D_F& tpoint);

    static double Distance(const Point4D_F& tpoint1, const Point4D_F& tpoint2);

private:

    float x;
    float y;
    float z;
    float w;
};

//a 2D images without scale and origin, type short
class Images2D_S
{
public :

      Images2D_S();
     ~Images2D_S();

    Images2D_S(byte* tbytes, int tstartByte, int tsize_x, int tsize_y);

    Images2D_S(int tsize_x, int tsize_y, short tinitialValue);

    bool InsideImage(int x, int y);

    void Save(ofstream& tbw);

    void Print(const char* tfileName);

    void Draw(const char* tfileName);

    short gV(int x, int y){ return image[x * size[1] + y]; }

    void sV(int x, int y, short tvalue){ image[x * size[1] + y] = tvalue; }

    __declspec ( property ( get = GetSize ) ) 
                int* Size ;

    int* GetSize(){ return size; }

    Point2D<int> GetMaxIndex();

    Point2D<int> GetMinIndex();

    Point2D<int> GetMaxValue( short& tresult);

    Point2D<int> GetMinValue( short& tresult);

    short GetInterp(double tx, double ty);

private:

    int size[2];
    short* image;
};

//a 3D images without scale and origin, type short
class Images3D_S
{
public:

     Images3D_S();
    ~Images3D_S();

    Images3D_S(std::string tfileName);

    Images3D_S(int tsizex, int tsizey, int tsizez, short tinitialValue);

    bool InsideImage(int x, int y, int z);

    bool SliceExist(int z);

    void Save(const char* tfileName);

    void PrintCrossSections(std::string tfileBase, int x, int y, int z);

    short gV(int x, int y, int z){ return image[z]->gV(x, y); }

    void sV(int x, int y, int z, short tvalue){ image[z]->sV(x, y, tvalue); }

    __declspec ( property ( get = GetSize ) ) 
                int* Size ;

    int* GetSize(){ return size; }

    void SetSlice(shared_ptr<Images2D_S> tslice, int z);

    Point3D<int> GetMaxIndex();

    Point3D<int> GetMinIndex();

    Point3D<int> GetMaxValue(short& tresult);

    Point3D<int> GetMinValue(short& tresult);

    shared_ptr<Images2D_S> SliceZ(int z);

    shared_ptr<Images2D_S> SliceX(int x);

    shared_ptr<Images2D_S> SliceY(int y);

    void CleanOrthogonal();

    shared_ptr<Images3D_S> CreateNormalizedImages(short tminValue = minValueS, short tmaxValue = maxValueS);

    short GetInterp(double tx, double ty, double tz);

private:
    
    Images3D_S(int tsizex, int tsizey, int tsizez);

    int size[3];
    std::vector<shared_ptr<Images2D_S>> image;
    std::vector<shared_ptr<Images2D_S>> imageX;
    std::vector<shared_ptr<Images2D_S>> imageY;
};

//a 2D images without scale and origin, type float
class Images2D_F
{

public :

     Images2D_F();
    ~Images2D_F();

    Images2D_F(byte* tbytes, int tstartByte, int tsize_x, int tsize_y);

    Images2D_F(int tsize_x, int tsize_y, float tinitialValue);

    bool InsideImage(int x, int y);

    void Save(ofstream& tbw);

    void Print(const char* tfileName);

    void Draw(const char* tfileName);

    float gV(int x, int y){ return image[x * size[1] + y]; }

    void sV(int x, int y, float tvalue){ image[x * size[1] + y] = tvalue; }

    __declspec ( property ( get = GetSize ) ) 
                int* Size ;

    int* GetSize(){ return size; }

    Point2D<int> GetMaxIndex();

    Point2D<int> GetMinIndex();

    Point2D<int> GetMaxValue( float& tresult);

    Point2D<int> GetMinValue( float& tresult);
    
private:

    int size[2];
    float* image;
};

//a 3D images without scale and origin, type float
class Images3D_F
{
public:

     Images3D_F();
    ~Images3D_F();

    Images3D_F(int tsizex, int tsizey, int tsizez, float tinitialValue);

    bool InsideImage(int x, int y, int z);

    float gV(int x, int y, int z){ return image[z]->gV(x, y); }

    void sV(int x, int y, int z, float tvalue){ image[z]->sV(x, y, tvalue); }

    __declspec ( property ( get = GetSize ) ) 
                int* Size ;

    int* GetSize(){ return size; }

    void SetSlice(shared_ptr<Images2D_F> tslice, int z);

    void Save(const char* tfileName);

    void PrintCrossSections(std::string tfileBase, int x, int y, int z);

    Point3D<int> GetMaxIndex();

    Point3D<int> GetMinIndex();

    Point3D<int> GetMaxValue(float& tresult);

    Point3D<int> GetMinValue(float& tresult);

    shared_ptr<Images2D_F> SliceZ(int z);

    shared_ptr<Images2D_F> SliceX(int x);

    shared_ptr<Images2D_F> SliceY(int y);

    bool SliceExist(int z);

    void CleanOrthogonal();

private:
    
    Images3D_F(int tsizex, int tsizey, int tsizez);

    int size[3];
    std::vector<shared_ptr<Images2D_F>> image;
    std::vector<shared_ptr<Images2D_F>> imageX;
    std::vector<shared_ptr<Images2D_F>> imageY;
};

class sImages2D_S
{
public :

     sImages2D_S();
    ~sImages2D_S();

    sImages2D_S(byte* tbytes, int tstartByte,
                int tsize_x, int tsize_y,
                float tscale_x,  float tscale_y,
                float torigin_x, float torigin_y);

    sImages2D_S(int tsize_x, int tsize_y, short tinitialValue,
                float tscale_x,  float tscale_y,
                float torigin_x, float torigin_y);

    void SetScaleAndOrigin(float tscale_x, float tscale_y, float torigin_x, float torigin_y);

    bool InsideImage(int x, int y);

    short gV(int x, int y){ return image[x * size[1] + y]; }

    void sV(int x, int y, short tvalue){ image[x * size[1] + y] = tvalue; }

    __declspec ( property ( get = GetSize ) ) 
                int* Size ;

    int* GetSize(){ return size; }

    __declspec ( property ( get = GetScale ) ) 
                float* Scale ;

    float* GetScale(){ return scale; }

    __declspec ( property ( get = GetOrigin ) ) 
                float* Origin ;

    float* GetOrigin(){ return origin; }

    short* Image() { return image; }

private:

    int size[2];
    float scale[2];
    float origin[2];
    short* image;
};

class sImages3D_S
{
public :

     sImages3D_S();
    ~sImages3D_S();

    sImages3D_S(int tsizex,         int tsizey,      int tsizez,
                float tscale_x,  float tscale_y,  float tscale_z,
                float torigin_x, float torigin_y, float torigin_z);

    bool SliceExist(int z);

    bool InsideImage(int x, int y, int z);

    short gV(int x, int y, int z){ return image[z]->gV(x, y); }

    void sV(int x, int y, int z, short tvalue){ image[z]->sV(x, y, tvalue); }

    __declspec ( property ( get = GetSize ) ) 
                int* Size ;

    int* GetSize(){ return size; }

    __declspec ( property ( get = GetScale ) ) 
                float* Scale ;

    float* GetScale(){ return scale; }

    __declspec ( property ( get = GetOrigin ) ) 
                float* Origin ;

    float* GetOrigin(){ return origin; }

private:
    int size[3];
    float scale[3];
    float origin[3];
    std::vector<shared_ptr<sImages2D_S>> image;
    std::vector<shared_ptr<sImages2D_S>> imageX;
    std::vector<shared_ptr<sImages2D_S>> imageY;
};

class BitmapCreator
{
public:

    static void SaveBitmapToFile( BYTE* pBitmapBits,  
                                  LONG lWidth,  
                                  LONG lHeight,  
                                  WORD wBitsPerPixel,
                                  const char* tfileName );

    static void SaveBitmapToFile(    std::vector<std::vector<Point3D<byte>>> pBitmapBits,  
                                    LONG lWidth,  
                                    LONG lHeight,  
                                    WORD wBitsPerPixel,
                                    const char* tfileName );
};

//a one slice of 3D gradient image
class GradientImage_2D
{
public :

    GradientImage_2D(int tsize_x, int tsize_y);

    void SetGradient(int x, int y, float tdx, float tdy, float tdz);

    void SetGradientDecomposed(int x, int y, float tdx, float tdy, float tdz);

    float* gV(int x, int y) { return gradients[x * size[1] + y]; }

    void sV(int x, int y, float* tvalue) { gradients[x * size[1] + y] = tvalue; }

    bool InsideImage(int x, int y);

private:

    int size[2];
    float** gradients;
};

//gradient image
class GradientImage_3D
{
public:

    GradientImage_3D(int tsizex, int tsizey, int tsizez);

    void SetGradient(int x, int y, int z, float tdx, float tdy, float tdz);

    void SetGradientDecomposed(int x, int y, int z, float tdx, float tdy, float tdz);

    float* gV(int x, int y, int z) { return gradients[z]->gV(x, y); }

    void sV(int x, int y, int z, float* tvalue) { gradients[z]->sV(x, y, tvalue); }

    float Mag(int x, int y, int z);

    __declspec ( property ( get = GetSize ) ) 
                int* Size ;

    int* GetSize(){ return size; }

    bool InsideImage(int x, int y, int z);

private:

    int size[3];
    std::vector<shared_ptr<GradientImage_2D>> gradients;
};


int roundF(float tvalue);
int roundD(double tvalue);
double randD();
float rand_thread();

#endif