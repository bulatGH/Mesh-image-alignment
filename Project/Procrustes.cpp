#include "AllTypes.h"
#include "Procrustes.h"


double Procrustes::ProcrustesMainVertices_NoScale(shared_ptr<Mesh> tmesh, std::vector<Point3D_F>& tnewPoints, std::vector<bool>& tselectedVerts)
{
    int tnumSelected = 0;
    for (size_t i = 0; i < tmesh->Vertices.size(); i++)
    {
        if (tselectedVerts[i])
        {
            tnumSelected++;
        }
    }
    if (tnumSelected == 0)
    {
        return 0;
    }
    Point3D_F toldCenter = FindMeanPointMesh(tmesh->Vertices, tselectedVerts, tnumSelected);
    tmesh->MoveObject(-toldCenter.X, -toldCenter.Y, -toldCenter.Z);
    Point3D_F tcenter = FindMeanPoint(tnewPoints, tselectedVerts, tnumSelected);


    std::vector<Point3D_F> tnewPoints_Norm(tnewPoints.size());
    for (size_t i = 0; i < tnewPoints.size(); i++)
    {
        if (tselectedVerts[i])
        {
            tnewPoints_Norm[i] = tnewPoints[i] - tcenter;
        }
    }
    double** tmart = RotationVertices(tmesh, tnewPoints_Norm, tselectedVerts, tnumSelected);

    double tdet = tmart[0][0] * tmart[1][1] * tmart[2][2];//alglib.rmatrixdet(tmart);
    int tlen = tnewPoints.size();

    parallel_for (int(0), tlen, [&](int i)
    {
        tmesh->Vertices[i]->coord = RotateOneVector(tmart, tmesh->Vertices[i]->coord);
    });
    tmesh->MoveObject(tcenter.X, tcenter.Y, tcenter.Z);

    tmesh->RecomputeNormals();
    tmesh->RecomputeFaceCenters();

    return CheckStopping(tdet, std::sqrt((toldCenter.X - tcenter.X) * (toldCenter.X - tcenter.X) +
                                         (toldCenter.Y - tcenter.Y) * (toldCenter.Y - tcenter.Y) +
                                         (toldCenter.Z - tcenter.Z) * (toldCenter.Z - tcenter.Z)));
}

double Procrustes::ProcrustesMainVertices_NoScale(shared_ptr<Mesh> tmesh, std::vector<Point3D_F>& tnewPoints)
{
    int tnumSelected = tmesh->Vertices.size();
    std::vector<bool> tselectedVerts;
    for (int i = 0; i < tmesh->Vertices.size(); i++)
    {
        tselectedVerts.push_back(true);
    }
    Point3D_F toldCenter = FindMeanPointMesh(tmesh->Vertices, tselectedVerts, tnumSelected);
    tmesh->MoveObject(-toldCenter.X, -toldCenter.Y, -toldCenter.Z);
    Point3D_F tcenter = FindMeanPoint(tnewPoints, tselectedVerts, tnumSelected);


    std::vector<Point3D_F> tnewPoints_Norm(tnewPoints.size());
    for (size_t i = 0; i < tnewPoints.size(); i++)
    {
        if (tselectedVerts[i])
        {
            tnewPoints_Norm[i] = tnewPoints[i] - tcenter;
        }
    }
    double** tmart = RotationVertices(tmesh, tnewPoints_Norm, tselectedVerts, tnumSelected);

    double tdet = tmart[0][0] * tmart[1][1] * tmart[2][2];//alglib.rmatrixdet(tmart);
    int tlen = tnewPoints.size();

    parallel_for(int(0), tlen, [&](int i)
    {
        tmesh->Vertices[i]->coord = RotateOneVector(tmart, tmesh->Vertices[i]->coord);
    });
    tmesh->MoveObject(tcenter.X, tcenter.Y, tcenter.Z);

    tmesh->RecomputeNormals();
    tmesh->RecomputeFaceCenters();

    return CheckStopping(tdet, std::sqrt((toldCenter.X - tcenter.X) * (toldCenter.X - tcenter.X) +
        (toldCenter.Y - tcenter.Y) * (toldCenter.Y - tcenter.Y) +
        (toldCenter.Z - tcenter.Z) * (toldCenter.Z - tcenter.Z)));
}

void Procrustes::ProcrustesMainVertices_GetParameters(shared_ptr<Mesh> tmesh, std::vector<Point3D_F>& trefPoints, std::vector<std::vector<double>>& trotation, std::vector<double> & ttrans)
{
    Point3D_F toldCenter = tmesh->Center;

    tmesh->PositionObject(0, 0, 0);
    Point3D_F tcenter = FindMeanPoint(trefPoints);
    std::vector<Point3D_F> trefPoints_Norm;// = new Vector3[tnewPoints.GetLength(0)];
    for (int i = 0; i < trefPoints.size(); i++)
    {
        trefPoints_Norm.push_back(trefPoints[i] - tcenter);
    }

    double** tmart = RotationVerticesInverse(tmesh, trefPoints_Norm);

    tmesh->PositionObject(toldCenter.X, toldCenter.Y, toldCenter.Z);

    trotation.clear();
    for (int i = 0; i < 3; i++)
    {
        trotation.push_back(std::vector<double>());
        for (int j = 0; j < 3; j++)
        {
            trotation[i].push_back(tmart[i][j]);
        }
    }
    ttrans[0] = (toldCenter.X - tcenter.X);
    ttrans[1] = (toldCenter.Y - tcenter.Y);
    ttrans[2] = (toldCenter.Z - tcenter.Z);
}

std::vector<Point3D_F> Procrustes::RotateAxes(std::vector<std::vector<double>>& trotation)
{
    std::vector<Point3D_F> taxes;
    taxes.push_back(RotateOneVector(trotation, Point3D_F(1, 0, 0)));
    taxes.push_back(RotateOneVector(trotation, Point3D_F(0, 1, 0)));
    taxes.push_back(RotateOneVector(trotation, Point3D_F(0, 0, 1)));
    return taxes;
}

Point3D_F Procrustes::FindMeanPoint(std::vector<Point3D_F>& tnewPoints)
{
    Point3D_F tvect = Point3D_F(0, 0, 0);
    for (size_t i = 0; i < tnewPoints.size(); i++)
    {
        tvect += tnewPoints[i];
    }
    tvect /= (float)tnewPoints.size();
    return tvect;
}

Point3D_F Procrustes::FindMeanPoint(std::vector<Point3D_F>& tnewPoints, std::vector<bool>& tselectedVerts, int tnumSelected)
{
    Point3D_F tvect = Point3D_F(0, 0, 0);
    for (size_t i = 0; i < tnewPoints.size(); i++)
    {
        if (tselectedVerts[i])
        {
            tvect += tnewPoints[i];
        }
    }
    tvect /= (float)tnumSelected;
    return tvect;
}

Point3D_F Procrustes::FindMeanPointMesh(std::vector<shared_ptr<Vertex>>& tnewPoints, std::vector<bool>& tselectedVerts, int tnumSelected)
{
    Point3D_F tvect = Point3D_F(0, 0, 0);
    for (size_t i = 0; i < tnewPoints.size(); i++)
    {
        if (tselectedVerts[i])
        {
            tvect += tnewPoints[i]->coord;
        }
    }
    tvect /= (float)tnumSelected;
    return tvect;
}

double** Procrustes::RotationInternalVertices(shared_ptr<Mesh> tmesh_Norm, std::vector<Point3D_F>& tnewPoints_Norm, std::vector<bool>& tselectedVerts, int tnumSelected)
{
    double** tres = new double*[3];
    for(int i = 0; i < 3; i++)
        tres[i] = new double[3];
    tres[0][0] = 0;
    tres[0][1] = 0;
    tres[0][2] = 0;
    tres[1][0] = 0;
    tres[1][1] = 0;
    tres[1][2] = 0;
    tres[2][0] = 0;
    tres[2][1] = 0;
    tres[2][2] = 0;

    for (size_t i = 0; i < tnewPoints_Norm.size(); i++)
    {
        if (tselectedVerts[i])
        {
            tres[0][0] += tmesh_Norm->Vertices[i]->coord.X * tnewPoints_Norm[i].X;
            tres[0][1] += tmesh_Norm->Vertices[i]->coord.Y * tnewPoints_Norm[i].X;
            tres[0][2] += tmesh_Norm->Vertices[i]->coord.Z * tnewPoints_Norm[i].X;

            tres[1][0] += tmesh_Norm->Vertices[i]->coord.X * tnewPoints_Norm[i].Y;
            tres[1][1] += tmesh_Norm->Vertices[i]->coord.Y * tnewPoints_Norm[i].Y;
            tres[1][2] += tmesh_Norm->Vertices[i]->coord.Z * tnewPoints_Norm[i].Y;

            tres[2][0] += tmesh_Norm->Vertices[i]->coord.X * tnewPoints_Norm[i].Z;
            tres[2][1] += tmesh_Norm->Vertices[i]->coord.Y * tnewPoints_Norm[i].Z;
            tres[2][2] += tmesh_Norm->Vertices[i]->coord.Z * tnewPoints_Norm[i].Z;
        }
    }
    return tres;
}

double** Procrustes::RotationInternalVerticesInverse(shared_ptr<Mesh> tmesh_Norm, std::vector<Point3D_F>& tnewPoints_Norm)
{
    double** tres = new double*[3];
    for (int i = 0; i < 3; i++)
        tres[i] = new double[3];
    tres[0][0] = 0;
    tres[0][1] = 0;
    tres[0][2] = 0;
    tres[1][0] = 0;
    tres[1][1] = 0;
    tres[1][2] = 0;
    tres[2][0] = 0;
    tres[2][1] = 0;
    tres[2][2] = 0;

    for (size_t i = 0; i < tnewPoints_Norm.size(); i++)
    {
        tres[0][0] += tmesh_Norm->Vertices[i]->coord.X * tnewPoints_Norm[i].X;
        tres[1][0] += tmesh_Norm->Vertices[i]->coord.Y * tnewPoints_Norm[i].X;
        tres[2][0] += tmesh_Norm->Vertices[i]->coord.Z * tnewPoints_Norm[i].X;

        tres[0][1] += tmesh_Norm->Vertices[i]->coord.X * tnewPoints_Norm[i].Y;
        tres[1][1] += tmesh_Norm->Vertices[i]->coord.Y * tnewPoints_Norm[i].Y;
        tres[2][1] += tmesh_Norm->Vertices[i]->coord.Z * tnewPoints_Norm[i].Y;

        tres[0][2] += tmesh_Norm->Vertices[i]->coord.X * tnewPoints_Norm[i].Z;
        tres[1][2] += tmesh_Norm->Vertices[i]->coord.Y * tnewPoints_Norm[i].Z;
        tres[2][2] += tmesh_Norm->Vertices[i]->coord.Z * tnewPoints_Norm[i].Z;
    }
    return tres;
}

Point3D_F Procrustes::RotateOneVector(double** tmart, Point3D_F told)
{
    Point3D_F tres = Point3D_F(0, 0, 0);
    tres.X = (float)(tmart[0][0] * told.X + tmart[0][1] * told.Y + tmart[0][2] * told.Z);
    tres.Y = (float)(tmart[1][0] * told.X + tmart[1][1] * told.Y + tmart[1][2] * told.Z);
    tres.Z = (float)(tmart[2][0] * told.X + tmart[2][1] * told.Y + tmart[2][2] * told.Z);
    return tres;
}

Point3D_F Procrustes::RotateOneVector(std::vector<std::vector<double>> tmart, Point3D_F told)
{
    Point3D_F tres = Point3D_F(0, 0, 0);
    tres.X = (float)(tmart[0][0] * told.X + tmart[0][1] * told.Y + tmart[0][2] * told.Z);
    tres.Y = (float)(tmart[1][0] * told.X + tmart[1][1] * told.Y + tmart[1][2] * told.Z);
    tres.Z = (float)(tmart[2][0] * told.X + tmart[2][1] * told.Y + tmart[2][2] * told.Z);
    return tres;
}

double** Procrustes::RotationVertices(shared_ptr<Mesh> tmesh_Norm, std::vector<Point3D_F>& tnewPoints_Norm)
{
    int tnumSelected = tnewPoints_Norm.size();
    std::vector<bool> tselectedVerts;
    for (int i = 0; i < tnewPoints_Norm.size(); i++)
    {
        tselectedVerts.push_back(true);
    }

    double** tM = RotationInternalVertices(tmesh_Norm, tnewPoints_Norm, tselectedVerts, tnumSelected);


    alglib::real_2d_array tM_a;
    tM_a.setlength(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            tM_a(i, j) = tM[i][j];

    alglib::real_1d_array w_a;
    alglib::real_2d_array u_a;
    alglib::real_2d_array vt_a;

    alglib::rmatrixsvd(tM_a, 3, 3, 1, 1, 2, w_a, u_a, vt_a);
    double** tres = new double*[3];
    for (int i = 0; i < 3; i++)
    {
        tres[i] = new double[3];
        for (int j = 0; j < 3; j++)
        {
            tres[i][j] = 0;
            for (int t = 0; t < 3; t++)
            {
                tres[i][j] += u_a(i, t) * vt_a(t, j);
            }
        }
    }
    return tres;
}

double** Procrustes::RotationVerticesInverse(shared_ptr<Mesh> tmesh_Norm, std::vector<Point3D_F>& tnewPoints_Norm)
{
    int tnumSelected = tnewPoints_Norm.size();
    std::vector<bool> tselectedVerts;
    for (int i = 0; i < tnewPoints_Norm.size(); i++)
    {
        tselectedVerts.push_back(true);
    }

    double** tM = RotationInternalVerticesInverse(tmesh_Norm, tnewPoints_Norm);


    alglib::real_2d_array tM_a;
    tM_a.setlength(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            tM_a(i, j) = tM[i][j];

    alglib::real_1d_array w_a;
    alglib::real_2d_array u_a;
    alglib::real_2d_array vt_a;

    alglib::rmatrixsvd(tM_a, 3, 3, 1, 1, 2, w_a, u_a, vt_a);
    double** tres = new double*[3];
    for (int i = 0; i < 3; i++)
    {
        tres[i] = new double[3];
        for (int j = 0; j < 3; j++)
        {
            tres[i][j] = 0;
            for (int t = 0; t < 3; t++)
            {
                tres[i][j] += u_a(i, t) * vt_a(t, j);
            }
        }
    }
    return tres;
}


double** Procrustes::RotationVertices(shared_ptr<Mesh> tmesh_Norm, std::vector<Point3D_F>& tnewPoints_Norm, std::vector<bool>& tselectedVerts, int tnumSelected)
{
    double** tM = RotationInternalVertices(tmesh_Norm, tnewPoints_Norm, tselectedVerts, tnumSelected);


    alglib::real_2d_array tM_a;
    tM_a.setlength(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            tM_a(i, j) = tM[i][j];

    alglib::real_1d_array w_a;
    alglib::real_2d_array u_a;
    alglib::real_2d_array vt_a;

    alglib::rmatrixsvd(tM_a, 3, 3, 1, 1, 2, w_a, u_a, vt_a);
    double** tres = new double*[3];
    for (int i = 0; i < 3; i++)
    {
        tres[i] = new double[3];
        for (int j = 0; j < 3; j++)
        {
            tres[i][j] = 0;
            for (int t = 0; t < 3; t++)
            {
                tres[i][j] += u_a(i, t) * vt_a(t, j);
            }
        }
    }
    return tres;
}

bool Procrustes::CheckStopping(float tscale, double tdet, double tdist)
{
    if (10 * std::abs(1 - tscale) + 10 * std::abs(1 - tdet) + tdist > 0.02)
    {
        return false;
    }
    else
    {
        return true;
    }
}

double Procrustes::CheckStopping(double tdet, double tdist)
{
    return 33300 * std::abs(1 - tdet) + 10 * tdist;
}
