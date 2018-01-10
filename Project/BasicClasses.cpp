#include "AllTypes.h"
#include "BasicClasses.h" 



short bytesToShort(byte* b, int tstart)
{
    short val = 0;
    int j = 0;
    for (int i = tstart+1; i >= tstart; --i)
    {
        val += (b[i] & 0xFF) << (8*j);
        ++j;
    }
 
    return val;
}

float bytesToFloat(byte* b, int tstart)
{
    float val = 0;
    int j = 0;
    for (int i = tstart+3; i >= tstart; --i)
    {
        val += (b[i] & 0xFF) << (8*j);
        ++j;
    }
 
    return val;
}

int roundF(float tvalue)
{
	return (int)(tvalue < 0.0 ? std::ceil(tvalue - 0.5) : std::floor(tvalue + 0.5));
}

int roundD(double tvalue)
{
	return (int)(tvalue < 0.0 ? std::ceil(tvalue - 0.5) : std::floor(tvalue + 0.5));
}

double randD()
{
	return ((double) rand() / (RAND_MAX));
}

// I am not sure about this function. 
//This should be a thread safe random function, if you have a better one, please change this function.
//Bulat
float rand_thread() {
	//static __declspec( thread ) std::mt19937 generator;
	static __declspec( thread ) mt19937* generator = nullptr;
    if (!generator) generator = new mt19937(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    uniform_real_distribution<float> distribution(0, 1);
    return distribution(*generator);
}

Point3D_F::Point3D_F() { }

Point3D_F Point3D_F::operator+(const Point3D_F& tpoint)
{
	return Point3D_F(this->x + tpoint.x, this->y + tpoint.y, this->z + tpoint.z);
}

Point3D_F Point3D_F::operator-(const Point3D_F& tpoint)
{
	return Point3D_F(this->x - tpoint.x, this->y - tpoint.y, this->z - tpoint.z);
}

Point3D_F Point3D_F::operator*(const float tcoef)
{
	return Point3D_F(this->x * tcoef, this->y * tcoef, this->z * tcoef);
}

Point3D_F Point3D_F::operator/(const float tcoef)
{
	return Point3D_F(this->x / tcoef, this->y / tcoef, this->z / tcoef);
}

void Point3D_F::operator+=(const Point3D_F& tpoint)
{
	this->x += tpoint.x;
	this->y += tpoint.y;
	this->z += tpoint.z;
}

void Point3D_F::operator-=(const Point3D_F& tpoint)
{
	this->x -= tpoint.x;
	this->y -= tpoint.y;
	this->z -= tpoint.z;
}

void Point3D_F::operator*=(const float tcoef)
{
	this->x *= tcoef;
	this->y *= tcoef;
	this->z *= tcoef;
}

void Point3D_F::operator/=(const float tcoef)
{
	this->x /= tcoef;
	this->y /= tcoef;
	this->z /= tcoef;
}

void Point3D_F::operator=(const Point3D_F& tpoint)
{
	x = tpoint.x;
	y = tpoint.y;
	z = tpoint.z;
}

bool Point3D_F::operator==(const Point3D_F& tpoint)
{
	return (tpoint.x == this->x) && (tpoint.y == this->y) && (tpoint.z == this->z);
}

double Point3D_F::Distance(const Point3D_F& tpoint1, const Point3D_F& tpoint2)
{
	return std::sqrt((tpoint1.x - tpoint2.x) * (tpoint1.x - tpoint2.x) + 
			 		 (tpoint1.y - tpoint2.y) * (tpoint1.y - tpoint2.y) + 
					 (tpoint1.z - tpoint2.z) * (tpoint1.z - tpoint2.z));
}

Point3D_F Point3D_F::Transform(Point3D_F tpoint1, double** tmatr)
{
	Point3D_F tpointRes = Point3D_F(0, 0, 0);
	tpointRes.X = (float)(tmatr[0][0] * tpoint1.X + tmatr[0][1] * tpoint1.Y + tmatr[0][2] * tpoint1.Z);
	tpointRes.Y = (float)(tmatr[1][0] * tpoint1.X + tmatr[1][1] * tpoint1.Y + tmatr[1][2] * tpoint1.Z);
	tpointRes.Z = (float)(tmatr[2][0] * tpoint1.X + tmatr[2][1] * tpoint1.Y + tmatr[2][2] * tpoint1.Z);
	return tpointRes;
}

double** Point3D_F::RotationYawPitchRoll(float dx, float dz, float dy)
{
	double** tmatr = new double*[3]; 
	for(int i = 0; i < 3; ++i)
		tmatr[i] = new double[3];

	tmatr[0][0] = std::cos(dx) * std::cos(dy);
	tmatr[0][1] = std::cos(dx) * std::sin(dy) * std::sin(dz) - std::sin(dx) * std::cos(dz);
	tmatr[0][2] = std::cos(dx) * std::sin(dy) * std::cos(dz) + std::sin(dx) * std::sin(dz);

	tmatr[1][0] = std::sin(dx) * std::cos(dy);
	tmatr[1][1] = std::sin(dx) * std::sin(dy) * std::sin(dz) + std::cos(dx) * std::cos(dz);
	tmatr[1][2] = std::sin(dx) * std::sin(dy) * std::cos(dz) - std::cos(dx) * std::sin(dz);

	tmatr[2][0] = -std::sin(dy);
	tmatr[2][1] = std::cos(dy) * std::sin(dz);
	tmatr[2][2] = std::cos(dy) * std::cos(dz);

	return tmatr;
}

double Point3D_F::ComputeDeterminant3_3(double** tmatr)
{
	return tmatr[0][0] * (tmatr[1][1] * tmatr[2][2] - tmatr[1][2] * tmatr[2][1]) - 
	       tmatr[0][1] * (tmatr[1][0] * tmatr[2][2] - tmatr[1][2] * tmatr[2][0]) +
	       tmatr[0][2] * (tmatr[1][0] * tmatr[2][1] - tmatr[1][1] * tmatr[2][0]);
}

double Point3D_F::Dot(Point3D_F& tpoint1, Point3D_F& tpoint2)
{
	return tpoint1.X * tpoint2.X + tpoint1.Y * tpoint2.Y + tpoint1.Z * tpoint2.Z;
}



Point4D_F::Point4D_F() { };

Point4D_F Point4D_F::operator+(const Point4D_F& tpoint)
{
	return Point4D_F(this->x + tpoint.x, this->y + tpoint.y, this->z + tpoint.z, this->w + tpoint.w);
}

Point4D_F Point4D_F::operator-(const Point4D_F& tpoint)
{
	return Point4D_F(this->x - tpoint.x, this->y - tpoint.y, this->z - tpoint.z, this->w + tpoint.w);
}

Point4D_F Point4D_F::operator*(const float tcoef)
{
	return Point4D_F(this->x * tcoef, this->y * tcoef, this->z * tcoef, this->w * tcoef);
}

Point4D_F Point4D_F::operator/(const float tcoef)
{
	return Point4D_F(this->x / tcoef, this->y / tcoef, this->z / tcoef, this->w / tcoef);
}

void Point4D_F::operator+=(const Point4D_F& tpoint)
{
	this->x += tpoint.x;
	this->y += tpoint.y;
	this->z += tpoint.z;
	this->w += tpoint.w;
}

void Point4D_F::operator-=(const Point4D_F& tpoint)
{
	this->x -= tpoint.x;
	this->y -= tpoint.y;
	this->z -= tpoint.z;
	this->w += tpoint.w;
}

void Point4D_F::operator*=(const float tcoef)
{
	this->x *= tcoef;
	this->y *= tcoef;
	this->z *= tcoef;
	this->w *= tcoef;
}

void Point4D_F::operator/=(const float tcoef)
{
	this->x /= tcoef;
	this->y /= tcoef;
	this->z /= tcoef;
	this->w /= tcoef;
}

void Point4D_F::operator=(const Point4D_F& tpoint)
{
	x = tpoint.x;
	y = tpoint.y;
	z = tpoint.z;
	w = tpoint.w;
}

bool Point4D_F::operator==(const Point4D_F& tpoint)
{
	return (tpoint.x == this->x) && (tpoint.y == this->y) && (tpoint.z == this->z) && (tpoint.w == this->w);
}

double Point4D_F::Distance(const Point4D_F& tpoint1, const Point4D_F& tpoint2)
{
	return std::sqrt((tpoint1.x - tpoint2.x) * (tpoint1.x - tpoint2.x) + 
			 		 (tpoint1.y - tpoint2.y) * (tpoint1.y - tpoint2.y) + 
			 		 (tpoint1.z - tpoint2.z) * (tpoint1.z - tpoint2.z) + 
			 		 (tpoint1.w - tpoint2.w) * (tpoint1.w - tpoint2.w));
}



Images2D_S::Images2D_S()
{

}

Images2D_S::~Images2D_S()
{
	if(image)
	{
		delete[] image;
	}
}

Images2D_S::Images2D_S(byte* tbytes, int tstartByte, int tsize_x, int tsize_y)
{
    size[0] = tsize_x;
    size[1] = tsize_y;
    image = new short[tsize_x * tsize_y];
    for (int i = 0; i < tsize_x * tsize_y; i++)
    {
        image[i] = bytesToShort(tbytes, tstartByte);
        tstartByte += 2;
    }
}

Images2D_S::Images2D_S(int tsize_x, int tsize_y, short tinitialValue)
{
    size[0] = tsize_x;
    size[1] = tsize_y;
    image = new short[tsize_x * tsize_y];
    for (int i = 0; i < tsize_x * tsize_y; i++)
    {
        image[i] = tinitialValue;
    }
}

bool Images2D_S::InsideImage(int x, int y)
{
    if (x < 0)
    {
        return false;
    }
    if (x > size[0] - 1)
    {
        return false;
    }

    if (y < 0)
    {
        return false;
    }
    if (y > size[1] - 1)
    {
        return false;
    }
    return true;
}

void Images2D_S::Save(ofstream& tbw)
{
    for (int i = 0; i < size[0] * size[1]; i++)
    {
	tbw.write((char * )&image[i], sizeof(short));
        //tbw << image[i];
    }
}

void Images2D_S::Print(const char* tfileName)
{
	ofstream tbw;
	tbw.open(tfileName, ios::out | ios::binary);
	tbw.write((char *) &size[0], sizeof(int));
	tbw.write((char *) &size[1], sizeof(int));
	for (int i = 0; i < size[0] * size[1]; i++)
    {
        tbw.write((char * )&image[i], sizeof(short));
    }
	tbw.close();
}

void Images2D_S::Draw(const char* tfileName)
{
	int tsizeX = (int)(4 * ceil((1.0 * this->Size[0] / 4)));
	int tsizeY = (int)(4 * ceil((1.0 * this->Size[1] / 4)));
	short tmin = 0;
	short tmax = 0;
	GetMaxValue(tmax);
	GetMinValue(tmin);
	BYTE* buf = new BYTE[ tsizeX * 4 * tsizeY ];
	int c = 0;
  
	for(short j = tsizeY - 1; j >= 0; j--)
    {
        for(short i = 0; i < tsizeX; i++)
        {
			if((i < this->Size[0]) && (j < this->Size[1]))
			{
				int tcol = (int)(255.0 * (this->gV(i, j) - tmin) / (tmax - tmin + 1));
				if(tcol < 20)
				{
					tcol = 20;
				}
				buf[ c + 0 ] = (BYTE) tcol;
				buf[ c + 1 ] = (BYTE) tcol;
				buf[ c + 2 ] = (BYTE) tcol;
				buf[ c + 3 ] = (BYTE) 0;
			}
			else
			{
				buf[ c + 0 ] = (BYTE) 0;
				buf[ c + 1 ] = (BYTE) 0;
				buf[ c + 2 ] = (BYTE) 0;
				buf[ c + 3 ] = (BYTE) 0;
			}
			c += 4;
		}
	}
	BitmapCreator::SaveBitmapToFile( (BYTE*) buf,
									 tsizeX,
									 tsizeY,
									 32,
									 tfileName );
	delete [] buf;
}

Point2D<int> Images2D_S::GetMaxIndex()
{
    short tmax = minValueS;
    int tmaxInd = -1;
    for (int i = 0; i < size[0] * size[1]; i++)
    {
        if (image[i] > tmax)
        {
            tmax = image[i];
            tmaxInd = i;
        }
    }
    Point2D<int> tresPoint;
	tresPoint.X = (tmaxInd / size[1]);
    tresPoint.Y = (tmaxInd - tresPoint.X * size[1]);

    return tresPoint;
}

Point2D<int> Images2D_S::GetMinIndex()
{
    short tmin = maxValueS;
    int tminInd = -1;
    for (int i = 0; i < size[0] * size[1]; i++)
    {
        if (image[i] < tmin)
        {
            tmin = image[i];
            tminInd = i;
        }
    }
	Point2D<int> tresPoint;
    tresPoint.X = (tminInd / size[1]);
    tresPoint.Y = (tminInd - tresPoint.X * size[1]);

    return tresPoint;
}

Point2D<int> Images2D_S::GetMaxValue( short& tresult)
{
	Point2D<int> tresPoint = this->GetMaxIndex();
	int tindex = tresPoint.X * size[1] + tresPoint.Y;
	tresult = image[tindex];
	return tresPoint;
}

Point2D<int> Images2D_S::GetMinValue( short& tresult)
{
    Point2D<int> tresPoint = this->GetMinIndex();
	int tindex = tresPoint.X * size[1] + tresPoint.Y;
	tresult = image[tindex];
	return tresPoint;
}

short Images2D_S::GetInterp(double tx, double ty)
{
	double xWei = 1 - tx + (int)tx;
	double yWei = 1 - ty + (int)ty;

	int xMin = (int)std::max(0.0, std::min(size[0] - 1.0, tx));
	int yMin = (int)std::max(0.0, std::min(size[1] - 1.0, ty));
	int xMax = std::min(size[0] - 1, xMin + 1);
	int yMax = std::min(size[1] - 1, yMin + 1);

	return (short)(image[xMin * size[1] + yMin] *	   xWei	 *      yWei +
				   image[xMin * size[1] + yMax] *      xWei  * (1 - yWei) +
				   image[xMax * size[1] + yMin] * (1 - xWei) *      yWei +
				   image[xMax * size[1] + yMax] * (1 - xWei) * (1 - yWei));
}



Images2D_F::Images2D_F()
{

}

Images2D_F::~Images2D_F()
{
	if(image)
	{
		delete[] image;
	}
}

Images2D_F::Images2D_F(byte* tbytes, int tstartByte, int tsize_x, int tsize_y)
{
	size[0] = tsize_x;
    size[1] = tsize_y;
    image = new float[tsize_x * tsize_y];
    for (int i = 0; i < tsize_x * tsize_y; i++)
    {
        image[i] = bytesToFloat(tbytes, tstartByte);
        tstartByte += 4;
    }
}

Images2D_F::Images2D_F(int tsize_x, int tsize_y, float tinitialValue)
{
    size[0] = tsize_x;
    size[1] = tsize_y;
    image = new float[tsize_x * tsize_y];
    for (int i = 0; i < tsize_x * tsize_y; i++)
    {
        image[i] = tinitialValue;
    }
}

bool Images2D_F::InsideImage(int x, int y)
{
    if (x < 0)
    {
        return false;
    }
    if (x > size[0] - 1)
    {
        return false;
    }

    if (y < 0)
    {
        return false;
    }
    if (y > size[1] - 1)
    {
        return false;
    }
    return true;
}

void Images2D_F::Save(ofstream& tbw)
{
	for (int i = 0; i < size[0] * size[1]; i++)
    {
		tbw.write((char * )&image[i], sizeof(float));
    }
}

void Images2D_F::Print(const char* tfileName)
{
	ofstream tbw;
	tbw.open(tfileName, ios::out | ios::binary);
	tbw.write((char *)(&size[0]), sizeof(int));
	tbw.write((char *)(&size[1]), sizeof(int));
	for (int i = 0; i < size[0] * size[1]; i++)
    {
        tbw.write((char * )&image[i], sizeof(float));
    }
	tbw.close();
}

void Images2D_F::Draw(const char* tfileName)
{
	int tsizeX = (int)(4 * ceil((1.0 * this->Size[0] / 4)));
	int tsizeY = (int)(4 * ceil((1.0 * this->Size[1] / 4)));
	float tmin = 0;
	float tmax = 0;
	this->GetMaxValue(tmax);
	this->GetMinValue(tmin);
	BYTE* buf = new BYTE[ tsizeX * 4 * tsizeY ];
	int c = 0;
	tmin = tmin > 0 ? 0 : tmin;
	for(short j = tsizeY - 1; j >= 0; j--)
    {
        for(short i = 0; i < tsizeX; i++)
        {
			if((i < this->Size[0]) && (j < this->Size[1]))
			{
				int tcol = (int)(255.0 * (this->gV(i, j) - tmin) / (tmax - tmin));
				if(tcol < 20)
				{
					tcol = 20;
				}
				buf[ c + 0 ] = (BYTE) tcol;
				buf[ c + 1 ] = (BYTE) tcol;
				buf[ c + 2 ] = (BYTE) tcol;
				buf[ c + 3 ] = (BYTE) 0;
			}
			else
			{
				buf[ c + 0 ] = (BYTE) 0;
				buf[ c + 1 ] = (BYTE) 0;
				buf[ c + 2 ] = (BYTE) 0;
				buf[ c + 3 ] = (BYTE) 0;
			}
			c += 4;
		}
	}
	BitmapCreator::SaveBitmapToFile( (BYTE*) buf,
									 tsizeX,
									 tsizeY,
									 32,
									 tfileName );
	delete [] buf;
}

Point2D<int> Images2D_F::GetMaxIndex()
{
    float tmax = minValueF;
    int tmaxInd = -1;
    for (int i = 0; i < size[0] * size[1]; i++)
    {
        if (image[i] > tmax)
        {
            tmax = image[i];
            tmaxInd = i;
        }
    }
    Point2D<int> tresPoint;
	tresPoint.X = (tmaxInd / size[1]);
    tresPoint.Y = (tmaxInd - tresPoint.X * size[1]);

    return tresPoint;
}

Point2D<int> Images2D_F::GetMinIndex()
{
    float tmin = maxValueF;
    int tminInd = -1;
    for (int i = 0; i < size[0] * size[1]; i++)
    {
        if (image[i] < tmin)
        {
            tmin = image[i];
            tminInd = i;
        }
    }
	Point2D<int> tresPoint;
    tresPoint.X = (tminInd / size[1]);
    tresPoint.Y = (tminInd - tresPoint.X * size[1]);

    return tresPoint;
}

Point2D<int> Images2D_F::GetMaxValue( float& tresult)
{
	Point2D<int> tresPoint = this->GetMaxIndex();
	int tindex = tresPoint.X * size[1] + tresPoint.Y;
	tresult = image[tindex];
	return tresPoint;
}

Point2D<int> Images2D_F::GetMinValue( float& tresult)
{
    Point2D<int> tresPoint = this->GetMinIndex();
	int tindex = tresPoint.X * size[1] + tresPoint.Y;
	tresult = image[tindex];
	return tresPoint;
}



Images3D_S::Images3D_S()
{

}

Images3D_S::~Images3D_S()
{
	for(int k = 0; k < size[2]; k++)
	{
		image[k].reset();
	}
	for(size_t k = 0; k < imageX.size(); k++)
	{
		imageX[k].reset();
	}
	for(size_t k = 0; k < imageY.size(); k++)
	{
		imageY[k].reset();
	}
}

Images3D_S::Images3D_S(std::string tfileName)
{
	ifstream tbr;
	tbr.open(tfileName, ios::in | ios::binary);

	tbr.read((char *)&size[0], sizeof(int));
	tbr.read((char *)&size[1], sizeof(int));
	tbr.read((char *)&size[2], sizeof(int));


	for (int k = 0; k < size[2]; k++)
	{
		shared_ptr<Images2D_S> image2(new Images2D_S(size[0], size[1], 0));
		std::vector<short> tval(size[0] * size[1]);
		tbr.read((char *)&tval[0], tval.size()*sizeof(tval[0]));
		int tcur = 0;
		for (int i = 0; i < size[0]; i++)
		{
			for (int j = 0; j < size[1]; j++)
			{
				image2->sV(i, j, tval[tcur]);
				tcur++;
			}
		}

		image.push_back(image2);
	}
	for (int i = 0; i < size[0]; i++)
	{
		shared_ptr<Images2D_S> image2(new Images2D_S(1, 1, 0));
		imageX.push_back(image2);
	}
	for (int j = 0; j < size[1]; j++)
	{
		shared_ptr<Images2D_S> image2(new Images2D_S(1, 1, 0));
		imageY.push_back(image2);
	}

	tbr.close();
}

Images3D_S::Images3D_S(int tsizex, int tsizey, int tsizez)
{
    size[0] = tsizex;
    size[1] = tsizey;
    size[2] = tsizez;

	for (int k = 0; k < size[2]; k++)
    {
		shared_ptr<Images2D_S> image2(new Images2D_S(size[0], size[1], 0));
		image.push_back(image2);
	}
	for(int i = 0; i < size[0]; i++)
	{
		shared_ptr<Images2D_S> image2(new Images2D_S(1, 1, 0));
		imageX.push_back(image2);
	}
	for(int j = 0; j < size[1]; j++)
	{
		shared_ptr<Images2D_S> image2(new Images2D_S(1, 1, 0));
		imageY.push_back(image2);
	}
}

Images3D_S::Images3D_S(int tsizex, int tsizey, int tsizez, short tinitialValue)
{
    size[0] = tsizex;
    size[1] = tsizey;
    size[2] = tsizez;

    for (int i = 0; i < size[2]; i++)
    {
		shared_ptr<Images2D_S> image2(new Images2D_S(size[0], size[1], tinitialValue));
		image.push_back(image2);
    }
	for(int i = 0; i < size[0]; i++)
	{
		shared_ptr<Images2D_S> image2(new Images2D_S(1, 1, 0));
		imageX.push_back(image2);
	}
	for(int j = 0; j < size[1]; j++)
	{
		shared_ptr<Images2D_S> image2(new Images2D_S(1, 1, 0));
		imageY.push_back(image2);
	}
}

bool Images3D_S::SliceExist(int z)
{
    if (z < 0)
    {
        return false;
    }
    if (z > size[2] - 1)
    {
        return false;
    }
    return true;
}

void Images3D_S::Save(const char* tfileName)
{
	ofstream tbw;
	tbw.open(tfileName, ios::out | ios::binary);
	tbw.write((char *)&size[0], sizeof(int));
	tbw.write((char *)&size[1], sizeof(int));
	tbw.write((char *)&size[2], sizeof(int));
	for (int k = 0; k < size[2]; k++)
	{
		image[k]->Save(tbw);
	}
	tbw.close();
}

void Images3D_S::PrintCrossSections(std::string tfileBase, int x, int y, int z)
{
	if((x >= 0) && (x < size[0]))
	{
		shared_ptr<Images2D_S> tslice = this->SliceX(x);
		std::string tfileNameS = tfileBase + "X.bmp";
		tslice->Draw(tfileNameS.c_str());
	}
	if((y >= 0) && (y < size[1]))
	{
		shared_ptr<Images2D_S> tslice = this->SliceY(y);
		std::string tfileNameS = tfileBase + "Y.bmp";
		tslice->Draw(tfileNameS.c_str());
	}
	if((z >= 0) && (z < size[2]))
	{
		shared_ptr<Images2D_S> tslice = this->SliceZ(z);
		std::string tfileNameS = tfileBase + "Z.bmp";
		tslice->Draw(tfileNameS.c_str());
	}
}

bool Images3D_S::InsideImage(int x, int y, int z)
{
    if (!SliceExist(z))
    {
        return false;
    }
    return image[z]->InsideImage(x, y);
}

void Images3D_S::SetSlice(shared_ptr<Images2D_S> tslice, int z)
{
	image[z] = tslice;
}

Point3D<int> Images3D_S::GetMaxIndex()
{
    Point3D<int> tres;
	short tmax = minValueS;
    for (int k = 0; k < size[2]; k++)
    {
		short tcurrMaxValue = 0;
        Point2D<int> tcurrMax = image[k]->GetMaxValue(tcurrMaxValue);
		if (tmax < tcurrMaxValue)
        {
			tres.X = tcurrMax.X;
			tres.Y = tcurrMax.Y;
            tres.Z = k;
			tmax = tcurrMaxValue;
        }
    }
    return tres;
}

Point3D<int> Images3D_S::GetMinIndex()
{
    Point3D<int> tres;
	short tmin = maxValueS;
    for (int k = 0; k < size[2]; k++)
    {
		short tcurrMinValue = 0;
        Point2D<int> tcurrMin = image[k]->GetMinValue(tcurrMinValue);
		if (tmin > tcurrMinValue)
        {
			tres.X = tcurrMin.X;
			tres.Y = tcurrMin.Y;
            tres.Z = k;
			tmin = tcurrMinValue;
        }
    }
    return tres;
}

Point3D<int> Images3D_S::GetMaxValue(short& tresult)
{
    Point3D<int> tres;
	short tmax = minValueS;
    for (int k = 0; k < size[2]; k++)
    {
		short tcurrMaxValue = 0;
        Point2D<int> tcurrMax = image[k]->GetMaxValue(tcurrMaxValue);
		if (tmax < tcurrMaxValue)
        {
			tres.X = tcurrMax.X;
			tres.Y = tcurrMax.Y;
            tres.Z = k;
			tmax = tcurrMaxValue;
        }
    }
	tresult = tmax;
    return tres;
}

Point3D<int> Images3D_S::GetMinValue(short& tresult)
{
    Point3D<int> tres;
	short tmin = maxValueS;
    for (int k = 0; k < size[2]; k++)
    {
		short tcurrMinValue = 0;
        Point2D<int> tcurrMin = image[k]->GetMinValue(tcurrMinValue);
		if (tmin > tcurrMinValue)
        {
			tres.X = tcurrMin.X;
			tres.Y = tcurrMin.Y;
            tres.Z = k;
			tmin = tcurrMinValue;
        }
    }
	tresult = tmin;
    return tres;
}

shared_ptr<Images2D_S> Images3D_S::SliceZ(int z)
{
	return image[z];
}

shared_ptr<Images2D_S> Images3D_S::SliceX(int x)
{
	if(imageX[x]->Size[0] + imageX[x]->Size[1] == 2)
	{
		shared_ptr<Images2D_S> additionalSlice(new Images2D_S(size[1], size[2], 0));
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				additionalSlice->sV(j, k, image[k]->gV(x, j));
			}
		}
		imageX[x] = additionalSlice;
	}
	return imageX[x];
}

shared_ptr<Images2D_S> Images3D_S::SliceY(int y)
{
	if(imageY[y]->Size[0] + imageY[y]->Size[1] == 2)
	{
		shared_ptr<Images2D_S> additionalSlice(new Images2D_S(size[0], size[2], 0));
		for (int i = 0; i < size[0]; i++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				additionalSlice->sV(i, k, image[k]->gV(i, y));
			}
		}
		imageY[y] = additionalSlice;
	}
	return imageY[y];
}

void Images3D_S::CleanOrthogonal()
{

	for(int k = 0; k < size[0]; k++)
	{
		imageX[k].reset();
	}
	for(int k = 0; k < size[1]; k++)
	{
		imageY[k].reset();
	}
	imageX.clear();
	imageY.clear();
	for(int i = 0; i < size[0]; i++)
	{
		shared_ptr<Images2D_S> image2(new Images2D_S(1, 1, 0));
		imageX.push_back(image2);
	}
	for(int j = 0; j < size[1]; j++)
	{
		shared_ptr<Images2D_S> image2(new Images2D_S(1, 1, 0));
		imageY.push_back(image2);
	}
}

shared_ptr<Images3D_S> Images3D_S::CreateNormalizedImages(short tminValue, short tmaxValue)
{
    shared_ptr<Images3D_S> timageNormalized(new Images3D_S(size[0], size[1], size[2]));
    short tmin_max[2];
	GetMinValue(tmin_max[0]);
	GetMaxValue(tmin_max[1]);
	
	tminValue = std::max(tmin_max[0], tminValue);
    tmaxValue = std::min(tmin_max[1], tmaxValue);
    tmin_max[1] = tmaxValue;
    //for (int k = 0; k < size[2]; k++)
    //Parallel.For(0, size[2], k =>
	parallel_for (int(0), size[2], [&](int k)
    {
        for (int i = 0; i < size[0]; i++)
        {
            for (int j = 0; j < size[1]; j++)
            {
				if (tmaxValue < image[k]->gV(i, j))
                {
					timageNormalized->sV(i, j, k, 255);
                }
                else if (tminValue > image[k]->gV(i, j))
				{
					timageNormalized->sV(i, j, k, 0);
				}
				else
                {
                    timageNormalized->sV(i, j, k, (short)(255.0 * (image[k]->gV(i, j) - tmin_max[0]) / (tmin_max[1] - tmin_max[0])));
                }
            }
        }
    });
	//printf ("%i %i \n", (int)tmin_max[0], (int)tmin_max[1]);
    return timageNormalized;
}

short Images3D_S::GetInterp(double tx, double ty, double tz)
{
	double zWei = 1 - tz + (int)tz;

	int zMin = (int)std::max(0.0, std::min(size[2] - 1.0, tz));
	int zMax = std::min(size[2] - 1, zMin + 1);

	return (short)(image[zMin]->GetInterp(tx, ty) *      zWei +
				   image[zMax]->GetInterp(tx, ty) * (1 - zWei));
}




Images3D_F::Images3D_F()
{

}

Images3D_F::~Images3D_F()
{
	for(int k = 0; k < size[2]; k++)
	{
		image[k].reset();
	}
	for(size_t k = 0; k < imageX.size(); k++)
	{
		imageX[k].reset();
	}
	for(size_t k = 0; k < imageY.size(); k++)
	{
		imageY[k].reset();
	}
}

Images3D_F::Images3D_F(int tsizex, int tsizey, int tsizez)
{
    size[0] = tsizex;
    size[1] = tsizey;
    size[2] = tsizez;

	for (int k = 0; k < size[2]; k++)
    {
		shared_ptr<Images2D_F> image2(new Images2D_F(size[0], size[1], 0));
		image.push_back(image2);
	}
	for(int i = 0; i < size[0]; i++)
	{
		shared_ptr<Images2D_F> image2(new Images2D_F(1, 1, 0));
		imageX.push_back(image2);
	}
	for(int j = 0; j < size[1]; j++)
	{
		shared_ptr<Images2D_F> image2(new Images2D_F(1, 1, 0));
		imageY.push_back(image2);
	}
}

Images3D_F::Images3D_F(int tsizex, int tsizey, int tsizez, float tinitialValue)
{
    size[0] = tsizex;
    size[1] = tsizey;
    size[2] = tsizez;

    for (int i = 0; i < size[2]; i++)
    {
		shared_ptr<Images2D_F> image2(new Images2D_F(size[0], size[1], tinitialValue));
		image.push_back(image2);
    }
	for(int i = 0; i < size[0]; i++)
	{
		shared_ptr<Images2D_F> image2(new Images2D_F(1, 1, 0));
		imageX.push_back(image2);
	}
	for(int j = 0; j < size[1]; j++)
	{
		shared_ptr<Images2D_F> image2(new Images2D_F(1, 1, 0));
		imageY.push_back(image2);
	}
}

bool Images3D_F::SliceExist(int z)
{
    if (z < 0)
    {
        return false;
    }
    if (z > size[2] - 1)
    {
        return false;
    }
    return true;
}

void Images3D_F::PrintCrossSections(std::string tfileBase, int x, int y, int z)
{
	if((x >= 0) && (x < size[0]))
	{
		shared_ptr<Images2D_F> tslice = this->SliceX(x);
		std::string tfileNameS = tfileBase + "X.bmp";
		tslice->Draw(tfileNameS.c_str());
	}
	if((y >= 0) && (y < size[1]))
	{
		shared_ptr<Images2D_F> tslice = this->SliceY(y);
		std::string tfileNameS = tfileBase + "Y.bmp";
		tslice->Draw(tfileNameS.c_str());
	}
	if((z >= 0) && (z < size[2]))
	{
		shared_ptr<Images2D_F> tslice = this->SliceZ(z);
		std::string tfileNameS = tfileBase + "Z.bmp";
		tslice->Draw(tfileNameS.c_str());
	}
}

bool Images3D_F::InsideImage(int x, int y, int z)
{
    if (!SliceExist(z))
    {
        return false;
    }
    return image[z]->InsideImage(x, y);
}

void Images3D_F::SetSlice(shared_ptr<Images2D_F> tslice, int z)
{
	image[z] = tslice;
}

void Images3D_F::Save(const char* tfileName)
{
	ofstream tbw;
	tbw.open(tfileName, ios::out | ios::binary);
	tbw.write((char *)&size[0], sizeof(int));
	tbw.write((char *)&size[1], sizeof(int));
	tbw.write((char *)&size[2], sizeof(int));
	for (int k = 0; k < size[2]; k++)
	{
		image[k]->Save(tbw);
	}
	tbw.close();
}

Point3D<int> Images3D_F::GetMaxIndex()
{
    Point3D<int> tres;
	float tmax = minValueF;
    for (int k = 0; k < size[2]; k++)
    {
		float tcurrMaxValue = 0;
        Point2D<int> tcurrMax = image[k]->GetMaxValue(tcurrMaxValue);
		if (tmax < tcurrMaxValue)
        {
			tres.X = tcurrMax.X;
			tres.Y = tcurrMax.Y;
            tres.Z = k;
			tmax = tcurrMaxValue;
        }
    }
    return tres;
}

Point3D<int> Images3D_F::GetMinIndex()
{
    Point3D<int> tres;
	float tmin = maxValueF;
    for (int k = 0; k < size[2]; k++)
    {
		float tcurrMinValue = 0;
        Point2D<int> tcurrMin = image[k]->GetMaxValue(tcurrMinValue);
		if (tmin > tcurrMinValue)
        {
			tres.X = tcurrMin.X;
			tres.Y = tcurrMin.Y;
            tres.Z = k;
			tmin = tcurrMinValue;
        }
    }
    return tres;
}

Point3D<int> Images3D_F::GetMaxValue(float& tresult)
{
    Point3D<int> tres;
	float tmax = minValueF;
    for (int k = 0; k < size[2]; k++)
    {
		float tcurrMaxValue = 0;
        Point2D<int> tcurrMax = image[k]->GetMaxValue(tcurrMaxValue);
		if (tmax < tcurrMaxValue)
        {
			tres.X = tcurrMax.X;
			tres.Y = tcurrMax.Y;
            tres.Z = k;
			tmax = tcurrMaxValue;
        }
    }
	tresult = tmax;
    return tres;
}

Point3D<int> Images3D_F::GetMinValue(float& tresult)
{
    Point3D<int> tres;
	float tmin = maxValueF;
    for (int k = 0; k < size[2]; k++)
    {
		float tcurrMinValue = 0;
        Point2D<int> tcurrMin = image[k]->GetMaxValue(tcurrMinValue);
		if (tmin > tcurrMinValue)
        {
			tres.X = tcurrMin.X;
			tres.Y = tcurrMin.Y;
            tres.Z = k;
			tmin = tcurrMinValue;
        }
    }
	tresult = tmin;
    return tres;
}

shared_ptr<Images2D_F> Images3D_F::SliceZ(int z)
{
	return image[z];
}

shared_ptr<Images2D_F> Images3D_F::SliceX(int x)
{
	if(imageX[x]->Size[0] + imageX[x]->Size[1] == 2)
	{
		shared_ptr<Images2D_F> additionalSlice(new Images2D_F(size[1], size[2], 0));
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				additionalSlice->sV(j, k, image[k]->gV(x, j));
			}
		}
		imageX[x] = additionalSlice;
	}
	return imageX[x];
}

shared_ptr<Images2D_F> Images3D_F::SliceY(int y)
{
	if(imageY[y]->Size[0] + imageY[y]->Size[1] == 2)
	{
		shared_ptr<Images2D_F> additionalSlice(new Images2D_F(size[0], size[2], 0));
		for (int i = 0; i < size[0]; i++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				additionalSlice->sV(i, k, image[k]->gV(i, y));
			}
		}
		imageY[y] = additionalSlice;
	}
	return imageY[y];
}

void Images3D_F::CleanOrthogonal()
{
	for(int k = 0; k < size[0]; k++)
	{
		imageX[k].reset();
	}
	for(int k = 0; k < size[1]; k++)
	{
		imageY[k].reset();
	}
	imageX.clear();
	imageY.clear();
	for(int i = 0; i < size[0]; i++)
	{
		shared_ptr<Images2D_F> image2(new Images2D_F(1, 1, 0));
		imageX.push_back(image2);
	}
	for(int j = 0; j < size[1]; j++)
	{
		shared_ptr<Images2D_F> image2(new Images2D_F(1, 1, 0));
		imageY.push_back(image2);
	}
}



sImages2D_S::sImages2D_S()
 {

 }

sImages2D_S::~sImages2D_S()
{
	delete[] image;
}

sImages2D_S::sImages2D_S(byte* tbytes, int tstartByte,
						 int tsize_x, int tsize_y,
						 float tscale_x,  float tscale_y,
						 float torigin_x, float torigin_y)
{
	size[0] = tsize_x;
    size[1] = tsize_y;
    scale[0] = tscale_x;
    scale[1] = tscale_y;
    origin[0] = torigin_x;
    origin[1] = torigin_y;
    image = new short[tsize_x * tsize_y];
    for (int i = 0; i < tsize_x * tsize_y; i++)
    {
        image[i] = bytesToShort(tbytes, tstartByte);
        tstartByte += 2;
    }

	size[0] = tsize_x;
    size[1] = tsize_y;
    image = new short[tsize_x * tsize_y];
    for (int i = 0; i < tsize_x * tsize_y; i++)
    {
        image[i] = bytesToShort(tbytes, tstartByte);
        tstartByte += 2;
    }
}

sImages2D_S::sImages2D_S(int tsize_x, int tsize_y, short tinitialValue,
						 float tscale_x,  float tscale_y,
						 float torigin_x, float torigin_y)
{
	size[0] = tsize_x;
    size[1] = tsize_y;
    scale[0] = tscale_x;
    scale[1] = tscale_y;
    origin[0] = torigin_x;
    origin[1] = torigin_y;
    image = new short[tsize_x * tsize_y];
    for (int i = 0; i < tsize_x * tsize_y; i++)
    {
        image[i] = tinitialValue;
    }
}

void sImages2D_S::SetScaleAndOrigin(float tscale_x, float tscale_y, float torigin_x, float torigin_y)
{
    scale[0] = tscale_x;
    scale[1] = tscale_y;
    origin[0] = torigin_x;
    origin[1] = torigin_x;
}

bool sImages2D_S::InsideImage(int x, int y)
{
    if (x < 0)
    {
        return false;
    }
    if (x > size[0] - 1)
    {
        return false;
    }

    if (y < 0)
    {
        return false;
    }
    if (y > size[1] - 1)
    {
        return false;
    }
    return true;
}


sImages3D_S::sImages3D_S()
{

}

sImages3D_S::~sImages3D_S()
{
	for(int k = 0; k < size[2]; k++)
	{
		image[k].reset();
	}
	for(size_t k = 0; k < imageX.size(); k++)
	{
		imageX[k].reset();
	}
	for(size_t k = 0; k < imageY.size(); k++)
	{
		imageY[k].reset();
	}
}

sImages3D_S::sImages3D_S(int tsizex,	  int tsizey,	   int tsizez,
					     float tscale_x,  float tscale_y,  float tscale_z,
					     float torigin_x, float torigin_y, float torigin_z)
{
    size[0] = tsizex;
    size[1] = tsizey;
    size[2] = tsizez;

    scale[0] = tscale_x;
    scale[1] = tscale_y;
    scale[2] = tscale_z;

    origin[0] = torigin_x;
    origin[1] = torigin_y;
    origin[2] = torigin_z;

	for (int k = 0; k < size[2]; k++)
    {
		shared_ptr<sImages2D_S> image2(new sImages2D_S(size[0],   size[1], 0, 
													   scale[0],  scale[1],
													   origin[0], origin[1]));
		image.push_back(image2);
	}
	for(int i = 0; i < size[0]; i++)
	{
		shared_ptr<sImages2D_S> image2(new sImages2D_S(1, 1, 0, 
													   scale[0],  scale[1],
													   origin[0], origin[1]));
		imageX.push_back(image2);
	}
	for(int j = 0; j < size[1]; j++)
	{
		shared_ptr<sImages2D_S> image2(new sImages2D_S(1, 1, 0, 
													   scale[0],  scale[1],
													   origin[0], origin[1]));
		imageY.push_back(image2);
	}
}


bool sImages3D_S::SliceExist(int z)
{
    if (z < 0)
    {
        return false;
    }
    if (z > size[2] - 1)
    {
        return false;
    }
    return true;
}

bool sImages3D_S::InsideImage(int x, int y, int z)
{
    if (!SliceExist(z))
    {
        return false;
    }
    return image[z]->InsideImage(x, y);
}



void BitmapCreator::SaveBitmapToFile( BYTE* pBitmapBits,  
									  LONG lWidth,  
									  LONG lHeight,  
									  WORD wBitsPerPixel,
									  const char* tfileName )  
{  
	int fileSize = 54 + (4 * lWidth * lHeight);

    //The sections of the file
    unsigned char generalHeader[14] = {'B','M',0,0, 0,0,0,0, 0,0,54,0, 0,0};
    unsigned char DIBHeader[40]     = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,32,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
    unsigned char pixelArray[1000000];
    unsigned char bmpPad[3] = {0, 0, 0};

    //Set the binary portion of the generalHeader, mainly just file size
    generalHeader[2] = (unsigned char)(fileSize);
    generalHeader[3] = (unsigned char)(fileSize >> 8);
    generalHeader[4] = (unsigned char)(fileSize >> 16);
    generalHeader[5] = (unsigned char)(fileSize >> 24);

    //The binary variable portion of the DIB header
    DIBHeader[4]  = (unsigned char)(lWidth);
    DIBHeader[5]  = (unsigned char)(lWidth >> 8);
    DIBHeader[6]  = (unsigned char)(lWidth >> 16);
    DIBHeader[7]  = (unsigned char)(lWidth >> 24);
    DIBHeader[8]  = (unsigned char)(lHeight);
    DIBHeader[9]  = (unsigned char)(lHeight >> 8);
    DIBHeader[10] = (unsigned char)(lHeight >> 16);
    DIBHeader[11] = (unsigned char)(lHeight >> 24);

    //Loop through all width and height places to add all pixels

	int pad = 0; // Set pad byte count per row to zero by default.
	// Each row needs to be a multiple of 4 bytes.  
	/*if ((lWidth * 4) % 4 != 0) 
	{
		pad = 4 - ((lWidth * 4) % 4);
	}*/
    int counter = 0;
    for(short i = 0; i < lHeight; i++)
    {
        for(short j = 0; j < lWidth; j++)
        {
            //Add all 3 RGB values
            pixelArray[counter++] = pBitmapBits[counter];//pixelColour[i][j].red;
            pixelArray[counter++] = pBitmapBits[counter];//pixelColour[i][j].green;
            pixelArray[counter++] = pBitmapBits[counter];//pixelColour[i][j].blue;
			pixelArray[counter++] = pBitmapBits[counter];
        }
		for (int padVal = 0; padVal < pad; padVal++) pixelArray[counter++] = 0;
    }

    //Open it
    ofstream fileWorking(tfileName);

    //Write the sections
    fileWorking.write((const char*)generalHeader, 14);
    fileWorking.write((const char*)DIBHeader, 40);
    fileWorking.write((const char*)pixelArray, (4 * lWidth + pad) * lHeight);

    //NO MEMORY LEAKS 4 ME
    fileWorking.close();
}  

void BitmapCreator::SaveBitmapToFile( std::vector<std::vector<Point3D<byte>>> pBitmapBits,  
									  LONG lWidth,  
									  LONG lHeight,  
									  WORD wBitsPerPixel,
									  const char* tfileName )  
{  
	int fileSize = 54 + (3 * lWidth * lHeight);

    //The sections of the file
    unsigned char generalHeader[14] = {'B','M',0,0, 0,0,0,0, 0,0,54,0, 0,0};
    unsigned char DIBHeader[40]     = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,32,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
    unsigned char pixelArray[1000000];
    unsigned char bmpPad[3] = {0, 0, 0};

    //Set the binary portion of the generalHeader, mainly just file size
    generalHeader[2] = (unsigned char)(fileSize);
    generalHeader[3] = (unsigned char)(fileSize >> 8);
    generalHeader[4] = (unsigned char)(fileSize >> 16);
    generalHeader[5] = (unsigned char)(fileSize >> 24);

    //The binary variable portion of the DIB header
    DIBHeader[4]  = (unsigned char)(lWidth);
    DIBHeader[5]  = (unsigned char)(lWidth >> 8);
    DIBHeader[6]  = (unsigned char)(lWidth >> 16);
    DIBHeader[7]  = (unsigned char)(lWidth >> 24);
    DIBHeader[8]  = (unsigned char)(lHeight);
    DIBHeader[9]  = (unsigned char)(lHeight >> 8);
    DIBHeader[10] = (unsigned char)(lHeight >> 16);
    DIBHeader[11] = (unsigned char)(lHeight >> 24);

    //Loop through all width and height places to add all pixels

	int pad = 0; // Set pad byte count per row to zero by default.
	// Each row needs to be a multiple of 4 bytes.  
	if ((lWidth * 3) % 4 != 0) 
	{
		pad = 4 - ((lWidth * 3) % 4);
	}
    int counter = 0;
    for(short j = lHeight; j >= 0; j--) 
	{
		for(short i = 0; i < lWidth; i++)
		{
            //Add all 3 RGB values
			pixelArray[counter++] = pBitmapBits[i][j].X;//pixelColour[i][j].red;
            pixelArray[counter++] = pBitmapBits[i][j].Y;//pixelColour[i][j].green;
            pixelArray[counter++] = pBitmapBits[i][j].Z;//pixelColour[i][j].blue;
			pixelArray[counter++] = 0;
            //counter += 3;
        }
		for (int padVal = 0; padVal < pad; padVal++) pixelArray[counter++] = 0;
    }

    //Open it
    ofstream fileWorking(tfileName);

    //Write the sections
    fileWorking.write((const char*)generalHeader, 14);
    fileWorking.write((const char*)DIBHeader, 40);
    fileWorking.write((const char*)pixelArray, (3 * lWidth + pad) * lHeight);

    //NO MEMORY LEAKS 4 ME
    fileWorking.close();
 
}  


GradientImage_2D::GradientImage_2D(int tsize_x, int tsize_y)
{
    size[0] = tsize_x;
    size[1] = tsize_y;
    gradients = new float*[tsize_x * tsize_y];
	for(int i = 0; i < tsize_x * tsize_y; i++)
	{
		gradients[i] = new float[4];
	}
}

void GradientImage_2D::SetGradient(int x, int y, float tdx, float tdy, float tdz)
{
	float tdist = (float)std::sqrt(tdx * tdx + tdy * tdy + tdz * tdz);
    if (tdist == 0)
    {
        tdist = (float)std::sqrt(3);
    }
    int tind = x * size[1] + y;
	gradients[tind][0] = tdx / tdist;
	gradients[tind][1] = tdy / tdist;
	gradients[tind][2] = tdz / tdist;
	gradients[tind][3] = tdist;
}

void GradientImage_2D::SetGradientDecomposed(int x, int y, float tdx, float tdy, float tdz)
{
    SetGradient(x, y, tdx, tdy, tdz);
}

bool GradientImage_2D::InsideImage(int x, int y)
{
    if (x < 0)
    {
        return false;
    }
    if (y < 0)
    {
        return false;
    }
    if (x > size[0] - 1)
    {
        return false;
    }
    if (y > size[1] - 1)
    {
        return false;
    }
    return true;
}


GradientImage_3D::GradientImage_3D(int tsizex, int tsizey, int tsizez)
{
    size[0] = tsizex;
    size[1] = tsizey;
    size[2] = tsizez;
    for (int k = 0; k < size[2]; k++)
    {
		gradients.push_back(std::make_shared<GradientImage_2D>(size[0], size[1]));
	}
}

void GradientImage_3D::SetGradient(int x, int y, int z, float tdx, float tdy, float tdz)
{
    gradients[z]->SetGradient(x, y, tdx, tdy, tdz);
}

void GradientImage_3D::SetGradientDecomposed(int x, int y, int z, float tdx, float tdy, float tdz)
{
    gradients[z]->SetGradientDecomposed(x, y, tdx, tdy, tdz);
}

float GradientImage_3D::Mag(int x, int y, int z)
{
	return gradients[z]->gV(x, y)[3];
}

bool GradientImage_3D::InsideImage(int x, int y, int z)
{
    if (x < 0)
    {
        return false;
    }
    if (y < 0)
    {
        return false;
    }
    if (z < 0)
    {
        return false;
    }
    if (x > size[0] - 1)
    {
        return false;
    }
    if (y > size[1] - 1)
    {
        return false;
    }
    if (z > size[2] - 1)
    {
        return false;
    }
    return true;
}

