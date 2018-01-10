#include "AllTypes.h"
#include "BasicClasses.h"
#include "Object_3D.h"
#include "Procrustes.h"
#include "Segmenter.h"



double RigidRegistration::DetectOptimalPointsSelectedVertices(shared_ptr<GradientImage_3D> tgradImage, std::vector<bool>& tselectedVerts)
{
	std::vector<Point3D_F> tvertices(meshMain->Vertices.size());
	int tlen = meshMain->Vertices.size();
	parallel_for (int(0), tlen, [&](int i)
    {
        if (tselectedVerts[i])
        {
            tvertices[i] = FindOptimalPoint_LowLevel(meshMain->Vertices[i]->norm, tgradImage, meshMain->Vertices[i]->coord);
        }
    });
    return Procrustes::ProcrustesMainVertices_NoScale(meshMain, tvertices, tselectedVerts);
}

double RigidRegistration::DetectOptimalPointsSelectedVertices(shared_ptr<GradientImage_3D> tgradImage)
{
	std::vector<Point3D_F> tvertices(meshMain->Vertices.size());
	int tlen = meshMain->Vertices.size();
	parallel_for(int(0), tlen, [&](int i)
	//for (int i = 0; i < tlen; i++)
	{
		tvertices[i] = FindOptimalPoint_LowLevel(meshMain->Vertices[i]->norm, tgradImage, meshMain->Vertices[i]->coord);
	});
	return Procrustes::ProcrustesMainVertices_NoScale(meshMain, tvertices);
}

double RigidRegistration::DetectOptimalPointsSelectedVertices_NoScale(shared_ptr<GradientImage_3D> tgradImage, std::vector<std::vector<double>>& trotation, std::vector<double> & ttrans)
{
	std::vector<Point3D_F> tvertices(meshMain->Vertices.size());
	int tlen = meshMain->Vertices.size();
	parallel_for(int(0), tlen, [&](int i)
	{
		tvertices[i] = FindOptimalPoint_LowLevel(meshMain->Vertices[i]->norm, tgradImage, meshMain->Vertices[i]->coord);
	});


	double tcorrect = 0;
	tcorrect = Procrustes::ProcrustesMainVertices_NoScale(meshMain, tvertices);
	return tcorrect;
}

double RigidRegistration::DetectOptimalPointsSelectedVertices_TwoLevel(shared_ptr<GradientImage_3D> tgradImage, std::vector<bool>& tselectedVerts, double tswitchValue)
{
	std::vector<Point3D_F> tvertices(meshMain->Vertices.size());
	int tlen = meshMain->Vertices.size();
	parallel_for(int(0), tlen, [&](int i)
	{
		if (tselectedVerts[i])
		{
			if (!highLevelFinished)
				tvertices[i] = FindOptimalPoint_HighLevel(meshMain->Vertices[i]->norm, tgradImage, meshMain->Vertices[i]->coord);
			else
				tvertices[i] = FindOptimalPoint_LowLevel( meshMain->Vertices[i]->norm, tgradImage, meshMain->Vertices[i]->coord);
		}
		else
		{
			tvertices[i] = Point3D_F(0, 0, 0);
		}
	});

	double tcorrect = Procrustes::ProcrustesMainVertices_NoScale(meshMain, tvertices, tselectedVerts);
	if (tcorrect <= tswitchValue)
	{
		highLevelFinished = !highLevelFinished;
		if (highLevelFinished)
			tcorrect *= 2;
	}
	return tcorrect;
}


Point3D_F RigidRegistration::FindOptimalPoint_LowLevel( Point3D_F tnormal, shared_ptr<GradientImage_3D> tgradImage, Point3D_F tpoint)
{
    Point3D_F tres = Point3D_F(0, 0, 0);
	double tbest = minValueF;
    for (int i = -profileLen; i <= profileLen; i++)
    {
        Point3D_F tpointLoc = tpoint + tnormal * step_delta * i;
        int tx = (int)round(tpointLoc.X);
        int ty = (int)round(tpointLoc.Y);
        int tz = (int)round(tpointLoc.Z);
        double tvalue = 0;
        if (tgradImage->InsideImage(tx, ty, tz))
        {
			float* tgradDir = tgradImage->gV(tx, ty, tz);
			float tgradMag = tgradImage->Mag(tx, ty, tz);
            tvalue = ComputeGradientCost_LowLevel(tnormal, tgradDir, tgradMag, i);
        }
        else
        {
            tvalue = ComputeGradientCost_Dummy(i);
        }
		
        if (tvalue > tbest)
        {
            tbest = tvalue;
            tres = tpointLoc;
        }
    }
    return tres;
}

Point3D_F RigidRegistration::FindOptimalPoint_HighLevel(Point3D_F tnormal, shared_ptr<GradientImage_3D> tgradImage, Point3D_F tpoint)
{
	Point3D_F tres = Point3D_F(0, 0, 0);
	double tbest = minValueF;
	for (int i = -profileLen; i <= profileLen; i++)
	{
		Point3D_F tpointLoc = tpoint + tnormal * step_delta * i;
		int tx = (int)round(tpointLoc.X);
		int ty = (int)round(tpointLoc.Y);
		int tz = (int)round(tpointLoc.Z);
		double tvalue = 0;
		if (tgradImage->InsideImage(tx, ty, tz))
		{
			float* tgradDir = tgradImage->gV(tx, ty, tz);
			float tgradMag = tgradImage->Mag(tx, ty, tz);
			tvalue = ComputeGradientCost_HighLevel(tnormal, tgradDir, tgradMag, i);
		}
		else
		{
			tvalue = ComputeGradientCost_Dummy(i);
		}

		if (tvalue > tbest)
		{
			tbest = tvalue;
			tres = tpointLoc;
		}
	}
	return tres;
}

double RigidRegistration::ComputeGradientCost_LowLevel( Point3D_F tnormal, float* tgradient, float tnorm, int tstep_j)
{
    return (tgradient[0] * tnorm * tnormal.X + tgradient[1] * tnorm * tnormal.Y + tgradient[2] * tnorm * tnormal.Z) * 
		   Gmax * (Gmax + tnorm) / (Gmax2 + tnorm * tnorm) - D * step_delta2 * tstep_j * tstep_j;
}

double RigidRegistration::ComputeGradientCost_HighLevel(Point3D_F tnormal, float* tgradient, float tnorm, int tstep_j)
{
	return (tgradient[0] * tnorm * tnormal.X + tgradient[1] * tnorm * tnormal.Y + tgradient[2] * tnorm * tnormal.Z);
}

double RigidRegistration::ComputeGradientCost_Dummy(int tstep_j)
{
    return -D * step_delta2 * tstep_j * tstep_j;
}
