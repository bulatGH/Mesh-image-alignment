#ifndef FILTERING_H
#define FILTERING_H

#include "AllTypes.h"
#include "BasicClasses.h"
#include "Object_3D.h"

class GradientFiltering;
class ImageRescaling;
class ReadStraightforwardly;

//class for genetating grandient images
class GradientFiltering
{
public:

    static shared_ptr<GradientImage_3D> GradientImage(shared_ptr<Images3D_S> timage);

	//gradient image is computed for only a specific region of the image
    static shared_ptr<GradientImage_3D> GradientImage(shared_ptr<Images3D_S> timage, shared_ptr<Images3D_S> tmask);
};


class ImageRescaling
{
public:
    
	//rescaling an sImage to Image with isotropic 1x1x1 resolution
    static shared_ptr<Images3D_S> RescaleUnit(shared_ptr<sImages3D_S> timage);

	//rescaling an sImage to Image with isotropic newScale resolution
    static shared_ptr<Images3D_S> RescaleISO(shared_ptr<sImages3D_S> timage, float tnewScale);

private:

	static short InterpolateValues(double x, double y, double z, double tscalex, double tscaley, double tscalez, shared_ptr<sImages3D_S> timage);
};

//this class is created to read images of a specific format
class ReadStraightforwardly
{
public:

    static shared_ptr<sImages3D_S> ReadOneImage(std::string tfileName, 
												int tsize_x,  int tsize_y, int tsize_z,
												float tscale_x, float tscale_y, float tscale_z,
												float torigin_x, float torigin_y, float torigin_z);
  
    static shared_ptr<sImages3D_S> ReadOneXLFormat(std::string tfileName);

private:

	static void FillOneXLFormat(std::string tline, int* tsize, float* tscale, float* torigin, int& tread);
};

class Threshold
{
public:

    static shared_ptr<Images3D_S> AboveThreshold(shared_ptr<Images3D_S> timage, short tT);

    static shared_ptr<Images3D_S> BelowThreshold(shared_ptr<Images3D_S> timage, short tT);

    static shared_ptr<Images3D_S> BetweenThreshold(shared_ptr<Images3D_S> timage, short tD, short tU);

	static shared_ptr<Images3D_S> AboveThreshold(shared_ptr<Images3D_F> timage, float tT);

    static shared_ptr<Images3D_S> BelowThreshold(shared_ptr<Images3D_F> timage, float tT);

    static shared_ptr<Images3D_S> BetweenThreshold(shared_ptr<Images3D_F> timage, float tD, float tU);
};

#endif