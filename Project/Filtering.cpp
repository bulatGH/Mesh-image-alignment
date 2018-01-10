#include "AllTypes.h"
#include "Filtering.h"


shared_ptr<GradientImage_3D> GradientFiltering::GradientImage(shared_ptr<Images3D_S> timage)
{
	shared_ptr<GradientImage_3D> tresultImage(new GradientImage_3D(timage->Size[0], timage->Size[1], timage->Size[2]));

	parallel_for (int(1), timage->Size[2] - 1, [&](int k)
	//for(int k = 1; k < timage->Size[2] - 1; k++)
	{
		for (int i = 1; i < timage->Size[0] - 1; i++)
		{
			for (int j = 1; j < timage->Size[1] - 1; j++)
			{
				tresultImage->SetGradient(i, j, k, 0.5f * (timage->gV(i + 1, j,	k) - timage->gV(i - 1, j, k)),
												   0.5f * (timage->gV(i, j + 1, k) - timage->gV(i, j - 1, k)),
												   0.5f * (timage->gV(i, j, k + 1) - timage->gV(i, j, k - 1)));
			}
		}
	});


	for (int i = 1; i < timage->Size[0] - 1; i++)
	{
		for (int j = 1; j < timage->Size[1] - 1; j++)
		{
			tresultImage->SetGradient(i, j, 0,
									  0.5f * (timage->gV(i + 1, j, 0) - timage->gV(i - 1, j, 0)),
									  0.5f * (timage->gV(i, j + 1, 0) - timage->gV(i, j - 1, 0)),
											  timage->gV(i, j, 1) - timage->gV(i, j, 0));
			tresultImage->SetGradient(i, j, timage->Size[2] - 1,
									  0.5f * (timage->gV(i + 1, j, timage->Size[2] - 1) - timage->gV(i - 1, j, timage->Size[2] - 1)),
									  0.5f * (timage->gV(i, j + 1, timage->Size[2] - 1) - timage->gV(i, j - 1, timage->Size[2] - 1)),
											  timage->gV(i, j, timage->Size[2] - 1) - timage->gV(i, j, timage->Size[2] - 2));
		}
	}
	for (int i = 1; i < timage->Size[0] - 1; i++)
	{
		for (int k = 1; k < timage->Size[2] - 1; k++)
		{
			tresultImage->SetGradient(i, 0, k,
									  0.5f * (timage->gV(i + 1, 0, k) - timage->gV(i - 1, 0, k)),
											  timage->gV(i, 1, k) - timage->gV(i, 0, k),
									  0.5f * (timage->gV(i, 0, k + 1) - timage->gV(i, 0, k - 1)));
			tresultImage->SetGradient(i, timage->Size[1] - 1, k,
									  0.5f * (timage->gV(i + 1, timage->Size[1] - 1, k) - timage->gV(i - 1, timage->Size[1] - 1, k)),
											  timage->gV(i, timage->Size[1] - 1, k) - timage->gV(i, timage->Size[1] - 2, k),
									  0.5f * (timage->gV(i, timage->Size[1] - 1, k + 1) - timage->gV(i, timage->Size[1] - 1, k - 1)));
		}
	}
	for (int j = 1; j < timage->Size[1] - 1; j++)
	{
		for (int k = 1; k < timage->Size[2] - 1; k++)
		{
			tresultImage->SetGradient(0, j, k,
											  timage->gV(1, j, k) - timage->gV(0, j, k),
									  0.5f * (timage->gV(0, j + 1, k) - timage->gV(0, j - 1, k)),
									  0.5f * (timage->gV(0, j, k + 1) - timage->gV(0, j, k - 1)));
			tresultImage->SetGradient(timage->Size[0] - 1, j, k,
											  timage->gV(timage->Size[0] - 1, j, k) - timage->gV(timage->Size[0] - 2, j, k),
									  0.5f * (timage->gV(timage->Size[0] - 1, j + 1, k) - timage->gV(timage->Size[0] - 1, j - 1, k)),
									  0.5f * (timage->gV(timage->Size[0] - 1, j, k + 1) - timage->gV(timage->Size[0] - 1, j, k - 1)));
		}
	}

	return tresultImage;
}

shared_ptr<GradientImage_3D> GradientFiltering::GradientImage(shared_ptr<Images3D_S> timage, shared_ptr<Images3D_S> tmask)
{
	shared_ptr<GradientImage_3D> tresultImage(new GradientImage_3D(timage->Size[0], timage->Size[1], timage->Size[2]));

	parallel_for (int(1), timage->Size[2] - 1, [&](int k)
	{
		for (int i = 1; i < timage->Size[0] - 1; i++)
		{
			for (int j = 1; j < timage->Size[1] - 1; j++)
			{
				if (tmask->gV(i, j, k) > 0)
				{
					tresultImage->SetGradient(i, j, k, 0.5f * (timage->gV(i + 1, j, k) - timage->gV(i - 1, j, k)),
													   0.5f * (timage->gV(i, j + 1, k) - timage->gV(i, j - 1, k)),
													   0.5f * (timage->gV(i, j, k + 1) - timage->gV(i, j, k - 1)));
				}
			}
		}
	});


	for (int i = 1; i < timage->Size[0] - 1; i++)
	{
		for (int j = 1; j < timage->Size[1] - 1; j++)
		{
			if (tmask->gV(i, j, 0) > 0)
			{
				tresultImage->SetGradient(i, j, 0,
										  0.5f * (timage->gV(i + 1, j, 0) - timage->gV(i - 1, j, 0)),
										  0.5f * (timage->gV(i, j + 1, 0) - timage->gV(i, j - 1, 0)),
												  timage->gV(i, j, 1) - timage->gV(i, j, 0));
			}
			if (tmask->gV(i, j, timage->Size[2] - 1) > 0)
			{
				tresultImage->SetGradient(i, j, timage->Size[2] - 1,
										  0.5f * (timage->gV(i + 1, j, timage->Size[2] - 1) - timage->gV(i - 1, j, timage->Size[2] - 1)),
										  0.5f * (timage->gV(i, j + 1, timage->Size[2] - 1) - timage->gV(i, j - 1, timage->Size[2] - 1)),
												  timage->gV(i, j, timage->Size[2] - 1) - timage->gV(i, j, timage->Size[2] - 2));
			}
		}
	}
	for (int i = 1; i < timage->Size[0] - 1; i++)
	{
		for (int k = 1; k < timage->Size[2] - 1; k++)
		{
			if (tmask->gV(i, 0, k) > 0)
			{
				tresultImage->SetGradient(i, 0, k,
										  0.5f * (timage->gV(i + 1, 0, k) - timage->gV(i - 1, 0, k)),
												  timage->gV(i, 1, k) - timage->gV(i, 0, k),
										  0.5f * (timage->gV(i, 0, k + 1) - timage->gV(i, 0, k - 1)));
			}
			if (tmask->gV(i, timage->Size[1] - 1, k) > 0)
			{
				tresultImage->SetGradient(i, timage->Size[1] - 1, k,
										  0.5f * (timage->gV(i + 1, timage->Size[1] - 1, k) - timage->gV(i - 1, timage->Size[1] - 1, k)),
												  timage->gV(i, timage->Size[1] - 1, k) - timage->gV(i, timage->Size[1] - 2, k),
										  0.5f * (timage->gV(i, timage->Size[1] - 1, k + 1) - timage->gV(i, timage->Size[1] - 1, k - 1)));
			}
		}
	}
	for (int j = 1; j < timage->Size[1] - 1; j++)
	{
		for (int k = 1; k < timage->Size[2] - 1; k++)
		{
			if (tmask->gV(0, j, k) > 0)
			{
				tresultImage->SetGradient(0, j, k,
												  timage->gV(1, j, k) - timage->gV(0, j, k),
										  0.5f * (timage->gV(0, j + 1, k) - timage->gV(0, j - 1, k)),
										  0.5f * (timage->gV(0, j, k + 1) - timage->gV(0, j, k - 1)));
			}
			if (tmask->gV(timage->Size[0] - 1, j, k) > 0)
			{
				tresultImage->SetGradient(timage->Size[0] - 1, j, k,
												  timage->gV(timage->Size[0] - 1, j, k) - timage->gV(timage->Size[0] - 2, j, k),
										  0.5f * (timage->gV(timage->Size[0] - 1, j + 1, k) - timage->gV(timage->Size[0] - 1, j - 1, k)),
										  0.5f * (timage->gV(timage->Size[0] - 1, j, k + 1) - timage->gV(timage->Size[0] - 1, j, k - 1)));
			}
		}
	}
	return tresultImage;
}


shared_ptr<Images3D_S> ImageRescaling::RescaleUnit(shared_ptr<sImages3D_S> timage)
{
	int tsizeX = (int)round(timage->Size[0] * timage->Scale[0]);
	int tsizeY = (int)round(timage->Size[1] * timage->Scale[1]);
	int tsizeZ = (int)round(timage->Size[2] * timage->Scale[2]);
			
	double tscalex = 1.0f * tsizeX / timage->Size[0];// / timage.Scale[0];
	double tscaley = 1.0f * tsizeY / timage->Size[1];// / timage.Scale[1];
	double tscalez = 1.0f * tsizeZ / timage->Size[2];// / timage.Scale[2];

	shared_ptr<Images3D_S> tresImage(new Images3D_S(tsizeX, tsizeY, tsizeZ, 0));

	parallel_for (int(0), tsizeX, [&](int i)
	{
		for (int j = 0; j < tsizeY; j++)
		{
			for (int k = 0; k < tsizeZ; k++)
			{
				tresImage->sV(i, j, k, InterpolateValues((i + 0.5) / tscalex, (j + 0.5) / tscaley, (k + 0.5) / tscalez,
														 tsizeX, tsizeY, tsizeZ, timage));
			}
		}
	});

	return tresImage;
}

shared_ptr<Images3D_S> ImageRescaling::RescaleISO(shared_ptr<sImages3D_S> timage, float tnewScale)
{
	double tscalex = timage->Scale[0] / tnewScale;// / timage.Scale[0];
	double tscaley = timage->Scale[1] / tnewScale;// / timage.Scale[1];
	double tscalez = timage->Scale[2] / tnewScale;// / timage.Scale[2];

	int tsizeX = (int)round(timage->Size[0] * timage->Scale[0] / tnewScale);
	int tsizeY = (int)round(timage->Size[1] * timage->Scale[1] / tnewScale);
	int tsizeZ = (int)round(timage->Size[2] * timage->Scale[2] / tnewScale);
			
	shared_ptr<Images3D_S> tresImage(new Images3D_S(tsizeX, tsizeY, tsizeZ, 0));
	parallel_for (int(0), tsizeX, [&](int i)
	{
		for (int j = 0; j < tsizeY; j++)
		{
			for (int k = 0; k < tsizeZ; k++)
			{
				float tval = InterpolateValues((i + 0.5) / tscalex, (j + 0.5) / tscaley, (k + 0.5) / tscalez,
											   tsizeX, tsizeY, tsizeZ, timage);
				tresImage->sV(i, j, k, tval);
			}
		}
	});
	return tresImage;
}

short ImageRescaling::InterpolateValues(double x, double y, double z, double tscalex, double tscaley, double tscalez, shared_ptr<sImages3D_S> timage)
{
	double tres = 0;
	double tsize = 1 / tscalex / tscaley / tscalez;
	double tx_st = x - 0.5 / tscalex;
	double tx_fi = x + 0.5 / tscalex;
	double ty_st = y - 0.5 / tscaley;
	double ty_fi = y + 0.5 / tscaley;
	double tz_st = z - 0.5 / tscalez;
	double tz_fi = z + 0.5 / tscalez;

	double tcontsumm = 0;
	double tval22 = fmod(tx_fi, 1);
	for (int i = (int)std::floor(tx_st); i < (fmod(tx_fi, 1) == 0.0 ? tx_fi + 1 : (int)std::ceil(tx_fi)); i++)
	{
		for (int j = (int)std::floor(ty_st); j < (fmod(ty_fi, 1) == 0.0 ? ty_fi + 1 : (int)std::ceil(ty_fi)); j++)
		{
			for (int k = (int)std::floor(tz_st); k < (fmod(tz_fi, 1) == 0.0 ? tz_fi + 1 : (int)std::ceil(tz_fi)); k++)
			{
				double taddx = 1;
				if (i == (int)std::floor(tx_st))
				{
					taddx = taddx - tx_st + i;
				}
				if (i == (fmod(tx_fi, 1) == 0.0 ? tx_fi + 1 : (int)std::ceil(tx_fi)) - 1)
				{
					taddx = taddx - (1 - tx_fi + i);
				}
				double taddy = 1;
				if (j == (int)std::floor(ty_st))
				{
					taddy = taddy - ty_st + j;
				}
				if (j == (fmod(ty_fi, 1) == 0.0 ? ty_fi + 1 : (int)std::ceil(ty_fi)) - 1)
				{
					taddy = taddy - (1 - ty_fi + j);
				}
				double taddz = 1;
				if (k == (int)std::floor(tz_st))
				{
					taddz = taddz - tz_st + k;
				}
				if (k == (fmod(tz_fi, 1) == 0.0 ? tz_fi + 1 : (int)std::ceil(tz_fi)) - 1)
				{
					taddz = taddz - (1 - tz_fi + k);
				}
				tres = tres + timage->gV(std::max(0, std::min(timage->Size[0] - 1, i)), std::max(0, std::min(timage->Size[1] - 1, j)), std::max(0, std::min(timage->Size[2] - 1, k))) * taddx * taddy * taddz;
				tcontsumm = tcontsumm + taddx * taddy * taddz;
			}
		}
	}
	return (short)round(tres / tsize);
}


shared_ptr<sImages3D_S> ReadStraightforwardly::ReadOneImage(std::string tfileName,
															int tsize_x,	 int tsize_y,	 int tsize_z,
															float tscale_x,  float tscale_y,  float tscale_z,
															float torigin_x, float torigin_y, float torigin_z)
{
	shared_ptr<sImages3D_S> tImage(new sImages3D_S(tsize_y,   tsize_x,   tsize_z,
												   tscale_y,  tscale_x,  tscale_z,
												   torigin_y, torigin_x, torigin_z));
	ifstream tbr; 
	tbr.open(tfileName, ios::in | ios::binary);

	std::vector<short> tval(tsize_z * tsize_x * tsize_y);
	tbr.read((char *)&tval[0], tval.size()*sizeof(tval[0]));
	int tcur = 0;
	for (int k = 0; k < tsize_z; k++)
	{
		for (int i = 0; i < tsize_x; i++)
		{
			for (int j = 0; j < tsize_y; j++)
			{
				tImage->sV(i, j, k, tval[tcur]);
				tcur++;
			}
		}
	}
	tbr.close();
	return tImage;
}
  
shared_ptr<sImages3D_S> ReadStraightforwardly::ReadOneXLFormat(std::string tfileName)
{
	std::ifstream tsr(tfileName);
	std::string tline;
	int* tsize = new int[3];
	float* tscale = new float[3];
	float* torigin = new float[3];
	int tread = 0;
	while (tread < 9)
	{
		std::getline(tsr, tline);
		FillOneXLFormat(tline, tsize, tscale, torigin, tread);
	}
	tsr.close();
	string tfileImage = tfileName.substr(0, tfileName.size() - 2) + "prinfo";
	return ReadOneImage(tfileImage, tsize[0],   tsize[1],   tsize[2],
									tscale[0],  tscale[1],  tscale[2],
									torigin[0], torigin[1], torigin[2]);
}

void ReadStraightforwardly::FillOneXLFormat(std::string tline, int* tsize, float* tscale, float* torigin, int& tread)
{
	if (tline.find("DIM_0") != std::string::npos)
	{
		size_t tsta = tline.find_first_of(">") + 1;
		size_t tfin = tline.find_last_of("<");
		int tval = std::stoi(tline.substr(tsta, tfin - tsta));
		tsize[1] = tval;
		tread++;
		return;
	}
	if (tline.find("DIM_1") != std::string::npos)
	{
		size_t tsta = tline.find_first_of(">") + 1;
		size_t tfin = tline.find_last_of("<");
		int tval = std::stoi(tline.substr(tsta, tfin - tsta));
		tsize[0] = tval;
		tread++;
		return;
	}
	if (tline.find("DIM_2") != std::string::npos)
	{
		size_t tsta = tline.find_first_of(">") + 1;
		size_t tfin = tline.find_last_of("<");
		int tval = std::stoi(tline.substr(tsta, tfin - tsta));
		tsize[2] = tval;
		tread++;
		return;
	}
	if (tline.find("MM_0") != std::string::npos)
	{
		size_t tsta = tline.find_first_of(">") + 1;
		size_t tfin = tline.find_last_of("<");
		float tval = std::stof(tline.substr(tsta, tfin - tsta));
		tscale[1] = tval;
		tread++;
		return;
	}
	if (tline.find("MM_1") != std::string::npos)
	{
		size_t tsta = tline.find_first_of(">") + 1;
		size_t tfin = tline.find_last_of("<");
		float tval = std::stof(tline.substr(tsta, tfin - tsta));
		tscale[0] = tval;
		tread++;
		return;
	}
	if (tline.find("MM_2") != std::string::npos)
	{
		size_t tsta = tline.find_first_of(">") + 1;
		size_t tfin = tline.find_last_of("<");
		float tval = std::stof(tline.substr(tsta, tfin - tsta));
		tscale[2] = tval;
		tread++;
		return;
	}
	if (tline.find("OFFSET_0") != std::string::npos)
	{
		size_t tsta = tline.find_first_of(">") + 1;
		size_t tfin = tline.find_last_of("<");
		float tval = std::stof(tline.substr(tsta, tfin - tsta));
		torigin[1] = tval;
		tread++;
		return;
	}
	if (tline.find("OFFSET_1") != std::string::npos)
	{
		size_t tsta = tline.find_first_of(">") + 1;
		size_t tfin = tline.find_last_of("<");
		float tval = std::stof(tline.substr(tsta, tfin - tsta));
		torigin[0] = tval;
		tread++;
		return;
	}
	if (tline.find("OFFSET_2") != std::string::npos)
	{
		size_t tsta = tline.find_first_of(">") + 1;
		size_t tfin = tline.find_last_of("<");
		float tval = std::stof(tline.substr(tsta, tfin - tsta));
		torigin[2] = tval;
		tread++;
		return;
	}
}



shared_ptr<Images3D_S> Threshold::AboveThreshold(shared_ptr<Images3D_S> timage, short tT)
{
	shared_ptr<Images3D_S> tmask(new Images3D_S(timage->Size[0], timage->Size[1], timage->Size[2], 0));
	for (int k = 0; k < timage->Size[2]; k++)
	{
		shared_ptr<Images2D_S> tslice = timage->SliceZ(k);
		for (int i = 0; i < timage->Size[0]; i++)
		{
			for (int j = 0; j < timage->Size[1]; j++)
			{
				if (tslice->gV(i, j) >= tT)
				{
					tmask->sV(i, j, k, 1);
				}
			}
		}
	}
	return tmask;
}

shared_ptr<Images3D_S> Threshold::BelowThreshold(shared_ptr<Images3D_S> timage, short tT)
{
	shared_ptr<Images3D_S> tmask(new Images3D_S(timage->Size[0], timage->Size[1], timage->Size[2], 0));
	for (int k = 0; k < timage->Size[2]; k++)
	{
		shared_ptr<Images2D_S> tslice = timage->SliceZ(k);
		for (int i = 0; i < timage->Size[0]; i++)
		{
			for (int j = 0; j < timage->Size[1]; j++)
			{
				if (tslice->gV(i, j) <= tT)
				{
					tmask->sV(i, j, k, 1);
				}
			}
		}
	}
	return tmask;
}

shared_ptr<Images3D_S> Threshold::BetweenThreshold(shared_ptr<Images3D_S> timage, short tD, short tU)
{
	shared_ptr<Images3D_S> tmask(new Images3D_S(timage->Size[0], timage->Size[1], timage->Size[2], 0));
	for (int k = 0; k < timage->Size[2]; k++)
	{
		shared_ptr<Images2D_S> tslice = timage->SliceZ(k);
		for (int i = 0; i < timage->Size[0]; i++)
		{
			for (int j = 0; j < timage->Size[1]; j++)
			{
				if ((tslice->gV(i, j) <= tU) && ((tslice->gV(i, j) >= tD)))
				{
					tmask->sV(i, j, k, 1);
				}
			}
		}
	}
	return tmask;
}

shared_ptr<Images3D_S> Threshold::AboveThreshold(shared_ptr<Images3D_F> timage, float tT)
{
	shared_ptr<Images3D_S> tmask(new Images3D_S(timage->Size[0], timage->Size[1], timage->Size[2], 0));
	for (int k = 0; k < timage->Size[2]; k++)
	{
		shared_ptr<Images2D_F> tslice = timage->SliceZ(k);
		for (int i = 0; i < timage->Size[0]; i++)
		{
			for (int j = 0; j < timage->Size[1]; j++)
			{
				if (tslice->gV(i, j) >= tT)
				{
					tmask->sV(i, j, k, 1);
				}
			}
		}
	}
	return tmask;
}

shared_ptr<Images3D_S> Threshold::BelowThreshold(shared_ptr<Images3D_F> timage, float tT)
{
	shared_ptr<Images3D_S> tmask(new Images3D_S(timage->Size[0], timage->Size[1], timage->Size[2], 0));
	for (int k = 0; k < timage->Size[2]; k++)
	{
		shared_ptr<Images2D_F> tslice = timage->SliceZ(k);
		for (int i = 0; i < timage->Size[0]; i++)
		{
			for (int j = 0; j < timage->Size[1]; j++)
			{
				if (tslice->gV(i, j) <= tT)
				{
					tmask->sV(i, j, k, 1);
				}
			}
		}
	}
	return tmask;
}

shared_ptr<Images3D_S> Threshold::BetweenThreshold(shared_ptr<Images3D_F> timage, float tD, float tU)
{
	shared_ptr<Images3D_S> tmask(new Images3D_S(timage->Size[0], timage->Size[1], timage->Size[2], 0));
	for (int k = 0; k < timage->Size[2]; k++)
	{
		shared_ptr<Images2D_F> tslice = timage->SliceZ(k);
		for (int i = 0; i < timage->Size[0]; i++)
		{
			for (int j = 0; j < timage->Size[1]; j++)
			{
				if ((tslice->gV(i, j) <= tU) && ((tslice->gV(i, j) >= tD)))
				{
					tmask->sV(i, j, k, 1);
				}
			}
		}
	}
	return tmask;
}