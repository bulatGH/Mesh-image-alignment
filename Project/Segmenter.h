#ifndef SEGMENTER_H
#define SEGMENTER_H

#include "AllTypes.h"
#include "Object_3D.h"

//this class alignes mesh with images
class RigidRegistration
{
public:

    double DetectOptimalPointsSelectedVertices(shared_ptr<GradientImage_3D> tgradImage, std::vector<bool>& tselectedVerts);

    double DetectOptimalPointsSelectedVertices(shared_ptr<GradientImage_3D> tgradImage);

    double DetectOptimalPointsSelectedVertices_NoScale(shared_ptr<GradientImage_3D> tgradImage, std::vector<std::vector<double>>& trotation, std::vector<double>& ttrans);

    double DetectOptimalPointsSelectedVertices_TwoLevel(shared_ptr<GradientImage_3D> tgradImage, std::vector<bool>& tselectedVerts, double tswitchValue);
    
    RigidRegistration(shared_ptr<Mesh> tmesh, int tprofileLen = 7, float tstep_delta = 1, float tD = 2, float tGmax = 100) 
        :meshMain(tmesh),
         profileLen(tprofileLen),
         step_delta(tstep_delta),
         D(tD),
         Gmax(tGmax),
         step_delta2(tstep_delta * tstep_delta),
         Gmax2(Gmax * Gmax)
        { };

private:

    double ComputeGradientCost_LowLevel( Point3D_F tnormal, float* tgradient, float tnorm, int tstep_j);

    static double ComputeGradientCost_HighLevel(Point3D_F tnormal, float* tgradient, float tnorm, int tstep_j);

    double ComputeGradientCost_Dummy(int tstep_j);

    Point3D_F FindOptimalPoint_LowLevel( Point3D_F tnormal, shared_ptr<GradientImage_3D> tgradImage, Point3D_F tpoint);

    Point3D_F FindOptimalPoint_HighLevel(Point3D_F tnormal, shared_ptr<GradientImage_3D> tgradImage, Point3D_F tpoint);


    shared_ptr<Mesh> meshMain;
    int profileLen;
    float step_delta;
    float D;
    float Gmax;
    float step_delta2;
    float Gmax2;
    bool highLevelFinished = false;
};

#endif