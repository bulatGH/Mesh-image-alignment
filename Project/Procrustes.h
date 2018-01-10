#ifndef PROCRUSTES_H
#define PROCRUSTES_H

#include "AllTypes.h"
#include "BasicClasses.h"
#include "Object_3D.h"

class Procrustes;

//this class is not literaly procrustes method, but rather different ways of aligning 3D meshes with 3D points.
class Procrustes
{
public:

    static double ProcrustesMainVertices_NoScale(shared_ptr<Mesh> tmesh, std::vector<Point3D_F>& tnewPoints);

    // only some mesh points are of interest (tselectedVerts).
    static double ProcrustesMainVertices_NoScale(shared_ptr<Mesh> tmesh, std::vector<Point3D_F>& tnewPoints, std::vector<bool>& tselectedVerts);

    // aligning with transformation parameters returned
    static void ProcrustesMainVertices_GetParameters(shared_ptr<Mesh> tmesh, std::vector<Point3D_F>& trefPoints, std::vector<std::vector<double>>& trotation, std::vector<double> & ttrans);

    static std::vector<Point3D_F> RotateAxes(std::vector<std::vector<double>>& trotation);

private:

    static Point3D_F FindMeanPoint(std::vector<Point3D_F>& tnewPoints);

    static Point3D_F FindMeanPoint(std::vector<Point3D_F>& tnewPoints, std::vector<bool>& tselectedVerts, int tnumSelected);

    static Point3D_F FindMeanPointMesh(std::vector<shared_ptr<Vertex>>& tnewPoints, std::vector<bool>& tselectedVerts, int tnumSelected);

    static double** RotationInternalVertices(shared_ptr<Mesh> tmesh_Norm, std::vector<Point3D_F>& tnewPoints_Norm, std::vector<bool>& tselectedVerts, int tnumSelected);

    static double** RotationInternalVerticesInverse(shared_ptr<Mesh> tmesh_Norm, std::vector<Point3D_F>& tnewPoints_Norm);

    static Point3D_F RotateOneVector(double** tmart, Point3D_F told);

    static Point3D_F RotateOneVector(std::vector<std::vector<double>> tmart, Point3D_F told);

    static double** RotationVertices(shared_ptr<Mesh> tmesh_Norm, std::vector<Point3D_F>& tnewPoints_Norm);

    static double** RotationVerticesInverse(shared_ptr<Mesh> tmesh_Norm, std::vector<Point3D_F>& tnewPoints_Norm);

    static double** RotationVertices(shared_ptr<Mesh> tmesh_Norm, std::vector<Point3D_F>& tnewPoints_Norm, std::vector<bool>& tselectedVerts, int tnumSelected);

    static bool CheckStopping(float tscale, double tdet, double tdist);

    static double CheckStopping(double tdet, double tdist);

    static Point3D_F FindMeanPoint(std::vector<std::vector<Point3D_F>>& tlineMain);
};

#endif