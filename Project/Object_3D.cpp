#include "AllTypes.h"
#include "BasicClasses.h"
#include "Object_3D.h"


Face::Face(Face& toriginal)
{
    A = toriginal.A;
    B = toriginal.B;
    C = toriginal.C;
	center = Point3D_F(toriginal.center.X, toriginal.center.Y, toriginal.center.Z);
    faceNormal = Point3D_F(toriginal.faceNormal.X, toriginal.faceNormal.Y, toriginal.faceNormal.Z);
}

Face::Face(Face& toriginal, const float tscale)
{
    A = toriginal.A;
    B = toriginal.B;
    C = toriginal.C;
    center = Point3D_F(toriginal.center.X, toriginal.center.Y, toriginal.center.Z);
    faceNormal = Point3D_F(toriginal.faceNormal.X, toriginal.faceNormal.Y, toriginal.faceNormal.Z);
    center *= tscale;
}



Vertex::Vertex(Vertex& toriginal)
{
    coord = Point3D_F(toriginal.coord.X, toriginal.coord.Y, toriginal.coord.Z);
    norm = Point3D_F(toriginal.norm.X, toriginal.norm.Y, toriginal.norm.Z);
    wCoord = Point3D_F(toriginal.wCoord.X, toriginal.wCoord.Y, toriginal.wCoord.Z);
}

Vertex::Vertex(Vertex& toriginal, const float tscale)
{
    coord = Point3D_F(toriginal.coord.X, toriginal.coord.Y, toriginal.coord.Z);
    norm = Point3D_F(toriginal.norm.X, toriginal.norm.Y, toriginal.norm.Z);
    wCoord = Point3D_F(toriginal.wCoord.X, toriginal.wCoord.Y, toriginal.wCoord.Z);
    coord *= tscale;
    wCoord *= tscale;
}



Mesh::Mesh(char* tname, int tverticesCount, int tfacesCount)
{
	Vertices.reserve(tverticesCount);
	Faces.reserve(tfacesCount);
    //faceNormals = new Vector3[tfacesCount];
	facesOfVerteces.reserve(tverticesCount);
    name = tname;
}

Mesh::Mesh(char* tfileName, float tscale)
{
	ifstream tbr; 
    tbr.open(tfileName, ios::in | ios::binary);

    //BinaryReader tbr = new BinaryReader(File.Open(tfileName, FileMode.Open));
    int tlen;// = tbr.ReadInt32();

	tbr.read((char *) &tlen, sizeof(int));

	float tx; float ty; float tz;
	tbr.read((char *) &tx, sizeof(float));
	tbr.read((char *) &ty, sizeof(float));
	tbr.read((char *) &tz, sizeof(float));
	Point3D_F center = Point3D_F(tx, ty, tz);


    for (int i = 0; i < tlen; i++)
    {
		std::shared_ptr<Vertex> tvert = std::shared_ptr<Vertex>();
		tbr.read((char *) &tx, sizeof(float));
		tbr.read((char *) &ty, sizeof(float));
		tbr.read((char *) &tz, sizeof(float));

		tvert->coord.X = tx;
		tvert->coord.Y = ty;
		tvert->coord.Z = tz;

		tvert->coord.X /= tscale; tvert->coord.Y /= tscale; tvert->coord.Z /= tscale;
		Vertices.push_back(tvert);
    }
	facesOfVerteces.reserve(tlen);// = new List<int>[tlen];
    for (int i = 0; i < tlen; i++)
    {
		float ttemp = 0;
		tbr.read((char *) &ttemp, sizeof(float));
		tbr.read((char *) &ttemp, sizeof(float));
		tbr.read((char *) &ttemp, sizeof(float));
    }
    tbr.read((char *) &tlen, sizeof(int));
    //faces = new Face[tlen];
    for (int i = 0; i < tlen; i++)
    {
		std::shared_ptr<Face> tface = std::shared_ptr<Face>();
		int tA;  int tB;  int tC;
		tbr.read((char *) &tA, sizeof(int));
		tbr.read((char *) &tB, sizeof(int));
		tbr.read((char *) &tC, sizeof(int));

		tface->A = tA; tface->B = tB; tface->C = tC;
        tface->center = ComputeFaceCenter(*Vertices[tface->A], *Vertices[tface->B], *Vertices[tface->C]);
		Faces.push_back(tface);
    }
    tbr.close();
    AssignVertecesToFaces();
    ComputeFaceNormals();
    ComputeVertexNormals();
    ComputeCenter();
}

Mesh::Mesh(ifstream tbr)
{
    int tlen;
	tbr.read((char *) &tlen, sizeof(int));

	float tx; float ty; float tz;
	tbr.read((char *) &tx, sizeof(float));
	tbr.read((char *) &ty, sizeof(float));
	tbr.read((char *) &tz, sizeof(float));
	Point3D_F center = Point3D_F(tx, ty, tz);

	/*Point3D_F center = Point3D_F();
	tbr.read((char *) &center.X, sizeof(float));
	tbr.read((char *) &center.Y, sizeof(float));
	tbr.read((char *) &center.Z, sizeof(float));*/

    for (int i = 0; i < tlen; i++)
    {
		std::shared_ptr<Vertex> tvert = std::shared_ptr<Vertex>();
		tbr.read((char *) &tx, sizeof(float));
		tbr.read((char *) &ty, sizeof(float));
		tbr.read((char *) &tz, sizeof(float));

		tvert->coord.X = tx;
		tvert->coord.Y = ty;
		tvert->coord.Z = tz;
		Vertices.push_back(tvert);
    }
	facesOfVerteces.reserve(tlen);// = new List<int>[tlen];
    for (int i = 0; i < tlen; i++)
    {
		float ttemp = 0;
		tbr.read((char *) &ttemp, sizeof(float));
		tbr.read((char *) &ttemp, sizeof(float));
		tbr.read((char *) &ttemp, sizeof(float));
    }
    tbr.read((char *) &tlen, sizeof(int));
    //faces = new Face[tlen];
    for (int i = 0; i < tlen; i++)
    {
		std::shared_ptr<Face> tface = std::shared_ptr<Face>();
		tbr.read((char *) &tface->A, sizeof(int));
		tbr.read((char *) &tface->B, sizeof(int));
		tbr.read((char *) &tface->C, sizeof(int));
        tface->center = ComputeFaceCenter(*Vertices[tface->A], *Vertices[tface->B], *Vertices[tface->C]);
		Faces.push_back(tface);
    }
    tbr.close();
    AssignVertecesToFaces();
    ComputeFaceNormals();
    ComputeVertexNormals();
    ComputeCenter();
}

Mesh::Mesh(std::vector<shared_ptr<Vertex>>& tvertices, std::vector<shared_ptr<Face>>& tfaces)
{
	for (size_t i = 0; i < tvertices.size(); i++)
    {
		Vertices.push_back(tvertices[i]);
    }
    for (size_t i = 0; i < tfaces.size(); i++)
    {
		Faces.push_back(tfaces[i]);
    }
	facesOfVerteces.reserve(tvertices.size());
    AssignVertecesToFaces();
    ComputeFaceNormals();
    ComputeVertexNormals();
    ComputeCenter();
}

Mesh::Mesh(std::vector<shared_ptr<Point3D_F>>& tvertices, std::vector<std::vector<int>>& tfaces, float tscale)
{
    for (size_t i = 0; i < tvertices.size(); i++)
    {
		shared_ptr<Vertex> tvert(new Vertex());
		tvert->coord.X = tvertices[i]->X * tscale;
		tvert->coord.Y = tvertices[i]->Y * tscale;
		tvert->coord.Z = tvertices[i]->Z * tscale;
		Vertices.push_back(tvert);//tvertices[i]);
    }
	//int tlen = tfaces.size();
	//Faces = std::vector<shared_ptr<Face>>(tlen);
	for (size_t i = 0; i < tfaces.size(); i++)
	//parallel_for (int(0), tlen, [&](int i)
    {
		std::shared_ptr<Face> tface(new Face());
		tface->A = tfaces[i][0];
        tface->B = tfaces[i][1];
        tface->C = tfaces[i][2];
        tface->center = ComputeFaceCenter(*Vertices[tface->A], *Vertices[tface->B], *Vertices[tface->C]);
		Faces.push_back(tface);
    }//);
    //facesOfVerteces.reserve(tvertices.size());

    AssignVertecesToFaces();
    ComputeFaceNormals();
    ComputeVertexNormals();
    ComputeCenter();
}

shared_ptr<Mesh> Mesh::ReadSTLMeshFlip(const char* tfileName, float newScale)
{
	ifstream tbr(tfileName, ios::in | ios::binary);
    //BinaryReader tbr = new BinaryReader(File.Open(tfileName, FileMode.Open));
    for (int i = 0; i < 10; i++)
    {
		double ttemp = 0;
		tbr.read((char *) &ttemp, sizeof(double));
        //tbr.ReadDouble();
    }
    int tnumFaces = 0;
	tbr.read((char *) &tnumFaces, sizeof(int));
	std::vector<Point4D_F> tvertices;// = new Vector4[3 * tnumFaces];
    std::vector<std::vector<int>> tfaces(tnumFaces, std::vector<int>(3));
    std::vector<Point3D_F> tfaceNorms;// = new Vector3[tnumFaces];
    for (int i = 0; i < tnumFaces; i++)
    {
		Point3D_F tfNorm = Point3D_F();
		float tx; float ty; float tz;
		tbr.read((char *) &tx, sizeof(float));
		tbr.read((char *) &ty, sizeof(float));
		tbr.read((char *) &tz, sizeof(float));
		tfNorm.X = tx;  tfNorm.Y = ty;  tfNorm.Z = tz;

		for(int j = 0; j < 3; j++)
		{
			Point4D_F tvert = Point4D_F();
			tbr.read((char *) &ty, sizeof(float));
			tbr.read((char *) &tx, sizeof(float));
			tbr.read((char *) &tz, sizeof(float));
			tvert.X = tx;  
			tvert.Y = ty;  
			tvert.Z = tz;

			tvert.W = 3 * i + j;
			tvertices.push_back(tvert);
		}
		
        //int tface[3];
        //tface[0] = 3 * i;
        //tface[1] = 3 * i + 1;
        //tface[2] = 3 * i + 2;
		tfaces[i][0] = 3 * i;
		tfaces[i][1] = 3 * i + 1;
		tfaces[i][2] = 3 * i + 2;

		short ttemp = 0;
		tbr.read((char *) &ttemp, sizeof(short));
        //tbr.ReadInt16();
    }
    tbr.close();

	std::sort(tvertices.begin(), tvertices.end(),
			  [](Point4D_F a, Point4D_F b) -> bool
	{
		return 10000 * b.X + 100 * b.Y + b.Z > 10000 * a.X + 100 * a.Y + a.Z;
	});

    //Array.Sort(tvertices, (o1, o2) => o1.X.CompareTo(o2.X));

	std::vector<int> tmarks;
	tmarks.reserve(tvertices.size());
	for(size_t i = 0; i < tvertices.size(); i++)
    {
        tmarks.push_back(-1);
    }

    int tselected = MarkToRemove(tvertices, tmarks);
    ReNumFaces(tfaces, tmarks, tvertices);
	std::vector<shared_ptr<Point3D_F>> tverticesCorr = ExtractCorectVertices(tmarks, tvertices, tselected);
    FlipFacesAndVertices(tfaces, tverticesCorr);
    shared_ptr<Mesh> tresult(new Mesh(tverticesCorr, tfaces, newScale));

	/*for(size_t i = 0; i < tfaces.size(); i++)
	{
		delete[] tfaces[i];
	}*/

    return tresult;
}

void Mesh::SaveAsSTL(const char* tfileBase, const char* tfileResult)
{
	ifstream tbr; 
    tbr.open(tfileBase, ios::ate | ios::binary);
	ifstream::pos_type pos = tbr.tellg();
    int length = pos;
    char *tbytes = new char[length];
    tbr.seekg(0, ios::beg);
    tbr.read(tbytes, length);
    tbr.close();
	ofstream tbw;
	tbw.open(tfileResult, ios::out | ios::binary);
    for (int i = 0; i < 80; i++)
    {
		tbw.write((char *) &tbytes[i], sizeof(char));
    }
    int tcurrent = 84;
	int tsize = Faces.size();
	tbw.write((char *) &tsize, sizeof(int));
    for (size_t i = 0; i < Faces.size(); i++)
    {
		float tx = -Faces[i]->faceNormal.X;
		float ty = -Faces[i]->faceNormal.Y;
		float tz = -Faces[i]->faceNormal.Z;
		tbw.write((char *) &ty, sizeof(float));
		tbw.write((char *) &tx, sizeof(float));
		tbw.write((char *) &tz, sizeof(float));

		tx = Vertices[Faces[i]->B]->coord.X;
		ty = Vertices[Faces[i]->B]->coord.Y;
		tz = Vertices[Faces[i]->B]->coord.Z;
		tbw.write((char *) &ty, sizeof(float));
		tbw.write((char *) &tx, sizeof(float));
		tbw.write((char *) &tz, sizeof(float));
		tx = Vertices[Faces[i]->A]->coord.X;
		ty = Vertices[Faces[i]->A]->coord.Y;
		tz = Vertices[Faces[i]->A]->coord.Z;
		tbw.write((char *) &ty, sizeof(float));
		tbw.write((char *) &tx, sizeof(float));
		tbw.write((char *) &tz, sizeof(float));
		tx = Vertices[Faces[i]->C]->coord.X;
		ty = Vertices[Faces[i]->C]->coord.Y;
		tz = Vertices[Faces[i]->C]->coord.Z;
		tbw.write((char *) &ty, sizeof(float));
		tbw.write((char *) &tx, sizeof(float));
		tbw.write((char *) &tz, sizeof(float));

        tcurrent += 48;
		tbw.write((char *) &tbytes[tcurrent],	  sizeof(char));
		tbw.write((char *) &tbytes[tcurrent + 1], sizeof(char));
        tcurrent += 2;
    }
	
    for (int i = tcurrent; i < (sizeof(tbytes)/sizeof(*tbytes)); i++)
    {
        tbw.write((char *) &tbytes[i], sizeof(char));
    }
    tbw.close();
	delete[] tbytes;
}

std::vector<int> Mesh::FOV(int i)
{
    return facesOfVerteces[i];
}

void Mesh::SetNewVerteces(std::vector<Point3D_F>& tvertices)
{
	for (size_t i = 0; i < Vertices.size(); i++)
    {
		Vertices[i]->coord = tvertices[i];
    }
    for (size_t i = 0; i < Faces.size(); i++)
    {
        Faces[i]->center = ComputeFaceCenter(*Vertices[Faces[i]->A], *Vertices[Faces[i]->B], *Vertices[Faces[i]->C]);
    }
    ComputeFaceNormals();
    ComputeVertexNormals();
    ComputeCenter();
}

void Mesh::PositionObject(float tx_g, float ty_g, float tz_g)
{
    MoveObject(tx_g - center.X, ty_g - center.Y, tz_g - center.Z);
}

void Mesh::MoveObject(float tx, float ty, float tz)
{
    Point3D_F tvect = Point3D_F(tx, ty, tz);

	int tlen = Vertices.size();
	parallel_for (int(0), tlen, [&](int t)
    {
        Vertices[t]->coord += tvect;
    });
    tlen = Faces.size();
    parallel_for (int(0), tlen, [&](int t)
    {
        Faces[t]->center += tvect;
    });
    center += tvect;
}

void Mesh::ScaleIndividualObject(float dx, float dy, float dz)
{
	int tlen = Vertices.size();
	parallel_for (int(0), tlen, [&](int i)
    {
        Vertices[i]->coord.X += ((Vertices[i]->coord.X - center.X) * (dx - 1));
		Vertices[i]->coord.Y += ((Vertices[i]->coord.Y - center.Y) * (dy - 1));
		Vertices[i]->coord.Z += ((Vertices[i]->coord.Z - center.Z) * (dz - 1));
    });
	tlen = Faces.size();
    parallel_for (int(0), tlen, [&](int i)
    {
        Faces[i]->center.X += ((Faces[i]->center.X - center.X) * (dx - 1));
		Faces[i]->center.Y += ((Faces[i]->center.Y - center.Y) * (dy - 1));
		Faces[i]->center.Z += ((Faces[i]->center.Z - center.Z) * (dz - 1));
    });
    center.X += (float)((center.X - center.X) * (dx - 1));
	center.Y += (float)((center.Y - center.Y) * (dy - 1));
	center.Z += (float)((center.Z - center.Z) * (dz - 1));
}

void Mesh::ScaleBrutally(float dx, float dy, float dz)
{
	for (int i = 0; i < (int)Vertices.size(); i++)
    {
        Vertices[i]->coord.X = (float)(Vertices[i]->coord.X * dx);
        Vertices[i]->coord.Y = (float)(Vertices[i]->coord.Y * dy);
        Vertices[i]->coord.Z = (float)(Vertices[i]->coord.Z * dz);
    }
    for (int i = 0; i < (int)Faces.size(); i++)
    {
        Faces[i]->center.X = (float)(Faces[i]->center.X * dx);
        Faces[i]->center.Y = (float)(Faces[i]->center.Y * dx);
        Faces[i]->center.Z = (float)(Faces[i]->center.Z * dx);
    }
}

void Mesh::Rotate(float dx, float dy, float dz)
{
    double** worldMatrix = Point3D_F::RotationYawPitchRoll(dy, dx, dz);// * Matrix.Translation(mesh.Position);
    //var worldMatrix = Quaternion.RotationYawPitchRoll(rotation.Y, rotation.X, rotation.Z);// * Matrix.Translation(mesh.Position);

	for (size_t i = 0; i < Vertices.size(); i++)
    {
        Vertices[i]->coord = RotateOneVertex(Vertices[i]->coord, worldMatrix);
        Vertices[i]->norm  = RotateOneVertex(Vertices[i]->norm, worldMatrix);
    }
    for (size_t i = 0; i < Faces.size(); i++)
    {
        Faces[i]->center = RotateOneVertex(Faces[i]->center, worldMatrix);
        Faces[i]->faceNormal = RotateOneVertex(Faces[i]->faceNormal, worldMatrix);
    }
	for(int i = 0; i < 3; i++)
	{
		delete[] worldMatrix[i];
	}
	delete[] worldMatrix;
}

void Mesh::RecomputeFaceCenters()
{
	int tlen = Faces.size();
	parallel_for (int(0), tlen, [&](int i)
    {
        Faces[i]->center = ComputeFaceCenter(*Vertices[Faces[i]->A], *Vertices[Faces[i]->B], *Vertices[Faces[i]->C]);
    });
}

void Mesh::RecomputeNormals()
{
    ComputeFaceNormals();
    ComputeVertexNormals();
}

void Mesh::ReNumFaces(std::vector<std::vector<int>>& tfaces, std::vector<int>& tmarks, std::vector<Point4D_F>& tvertices)
{
    for (size_t i = 0; i < tmarks.size(); i++)
    {
        int tind = (int)tvertices[i].W / 3;
        int tadd = (int)tvertices[i].W - tind * 3;
        tfaces[tind][tadd] = tmarks[i];
    }
}

std::vector<shared_ptr<Point3D_F>> Mesh::ExtractCorectVertices(std::vector<int>& tmarks, std::vector<Point4D_F>& tvertices, int tlen)
{
    std::vector<shared_ptr<Point3D_F>> tresVerts;
    int tind = 0;
    for (int i = 0; i < tlen; i++)
    {
        while(tmarks[tind] != i)
        {
            tind++;
        }
		tresVerts.push_back(make_shared<Point3D_F>(tvertices[tind].X, tvertices[tind].Y, tvertices[tind].Z));
    }
    return tresVerts;
}

int Mesh::MarkToRemove(std::vector<Point4D_F>& tvertices, std::vector<int>& tmarks)
{
    size_t i = 0;
    int tcurrent = 0;
	while (i < tvertices.size())
    {
        if (tmarks[i] > -1)
        {
            i++;
            continue;
        }

        tmarks[i] = tcurrent;
        int j = i + 1;
		if (j == tvertices.size())
        {
            tcurrent++;
            break;
        }
        while (std::abs(tvertices[i].X - tvertices[j].X) < 1e-6)
        {
            if (std::sqrt((tvertices[i].X - tvertices[j].X) * (tvertices[i].X - tvertices[j].X) +
                          (tvertices[i].Y - tvertices[j].Y) * (tvertices[i].Y - tvertices[j].Y) +
                          (tvertices[i].Z - tvertices[j].Z) * (tvertices[i].Z - tvertices[j].Z))< 1e-3)
            {
                tmarks[j] = tcurrent;
            }
            j++;
			if (j == tvertices.size())
            {
                break;
            }
        }
        i++;
        tcurrent++;
    }
    return tcurrent;
}

void Mesh::FlipFacesAndVertices(std::vector<std::vector<int>>& tfaces, std::vector<shared_ptr<Point3D_F>>& tverticesCorr)
{
    for (size_t i = 0; i < tfaces.size(); i++)
    {
        int t = tfaces[i][0];
        tfaces[i][0] = tfaces[i][1];
        tfaces[i][1] = t;
    }
    return;
    float taverZ = 0;
	for (size_t i = 0; i < tverticesCorr.size(); i++)
    {
        taverZ += tverticesCorr[i]->Z;
    }
    taverZ /= tverticesCorr.size();
    for (size_t i = 0; i < tverticesCorr.size(); i++)
    {
        tverticesCorr[i]->Z = 2 * taverZ - tverticesCorr[i]->Z;
    }
    taverZ = 0;
    for (size_t i = 0; i < tverticesCorr.size(); i++)
    {
        taverZ += tverticesCorr[i]->Z;
    }
    taverZ /= tverticesCorr.size();
}

void Mesh::ComputeCenter()
{
    float taver_x = 0;
    float taver_y = 0;
    float taver_z = 0;
    for (size_t i = 0; i < Vertices.size(); i++)
    {
        taver_x += Vertices[i]->coord.X;
        taver_y += Vertices[i]->coord.Y;
        taver_z += Vertices[i]->coord.Z;
    }
	center = Point3D_F(taver_x / Vertices.size(), taver_y / Vertices.size(), taver_z / Vertices.size());
}

Point3D_F Mesh::ComputeFaceCenter(Vertex a, Vertex b, Vertex c)
{
    Point3D_F tres = Point3D_F();
    tres.X += (a.coord.X + b.coord.X + c.coord.X) / 3.0f;
    tres.Y += (a.coord.Y + b.coord.Y + c.coord.Y) / 3.0f;
    tres.Z += (a.coord.Z + b.coord.Z + c.coord.Z) / 3.0f;
    return tres;
}

Point3D_F Mesh::VectorOffset(Point3D_F pIn, Point3D_F pOffset)
{
    Point3D_F pOut = Point3D_F();
    pOut.X = pIn.X - pOffset.X;
    pOut.Y = pIn.Y - pOffset.Y;
    pOut.Z = pIn.Z - pOffset.Z;
    return pOut;
}

Point3D_F Mesh::VectorGetNormal(Point3D_F a, Point3D_F b)
{
    Point3D_F pOut = Point3D_F();
    pOut.X = a.Y * b.Z - a.Z * b.Y;
    pOut.Y = a.Z * b.X - a.X * b.Z;
    pOut.Z = a.X * b.Y - a.Y * b.X;
    return pOut;
}

Point3D_F Mesh::VectorNormalize(Point3D_F pIn)
{
    Point3D_F pOut = Point3D_F();
    float len = (float)(std::sqrt((pIn.X * pIn.X) + (pIn.Y * pIn.Y) + (pIn.Z * pIn.Z)));
    if (len > 0)
    {
        pOut.X = pIn.X / len;
        pOut.Y = pIn.Y / len;
        pOut.Z = pIn.Z / len;
    }
    return pOut;
}

Point3D_F Mesh::ComputeFNormal(Point3D_F p1, Point3D_F p2, Point3D_F p3)
{
    // Uses p2 as a new origin for p1,p3
    Point3D_F a = VectorOffset(p3, p2);
    Point3D_F b = VectorOffset(p1, p2);
    Point3D_F pn = VectorGetNormal(a, b);
	pn.X = -pn.X; pn.Y = -pn.Y; pn.Z = -pn.Z;
    return VectorNormalize(pn);
}

void Mesh::ComputeFaceNormals()
{
	int tlen = Faces.size();
	parallel_for (int(0), tlen, [&](int i)
    {
        Faces[i]->faceNormal = ComputeFNormal(Vertices[Faces[i]->A]->coord, Vertices[Faces[i]->B]->coord, Vertices[Faces[i]->C]->coord);
    });
}

void Mesh::AssignVertecesToFaces()
{
	facesOfVerteces = std::vector<std::vector<int>>(Vertices.size(), std::vector<int>());
    /*for (size_t t = 0; t < Vertices.size(); t++)
    {
		facesOfVerteces.push_back(std::vector<int>());
    }*/
    for (size_t i = 0; i < Faces.size(); i++)
    {
		facesOfVerteces[Faces[i]->A].push_back(i);
        facesOfVerteces[Faces[i]->B].push_back(i);
        facesOfVerteces[Faces[i]->C].push_back(i);
    }
}

Point3D_F Mesh::ComputeVNormal(int tnum)
{
    Point3D_F pOut = Point3D_F(0, 0, 0);
    for (size_t i = 0; i < facesOfVerteces[tnum].size(); i++)
    {
        int j = facesOfVerteces[tnum][i];
        pOut.X += Faces[j]->faceNormal.X;// faceNormals[j].X;
        pOut.Y += Faces[j]->faceNormal.Y;//faceNormals[j].Y;
        pOut.Z += Faces[j]->faceNormal.Z;//faceNormals[j].Z;
    }
	pOut.X /= facesOfVerteces[tnum].size();
    pOut.Y /= facesOfVerteces[tnum].size();
    pOut.Z /= facesOfVerteces[tnum].size();

    return VectorNormalize(pOut);

}

void Mesh::ComputeVertexNormals()
{
	int tlen = Vertices.size();
	parallel_for (int(0), tlen, [&](int i)
    //Parallel.For(0, vertices.GetLength(0), i =>
    {
        Vertices[i]->norm = ComputeVNormal(i);
    });
}

Point3D_F Mesh::RotateOneVertex(Point3D_F tpoint, double** tMatr)
{
    Point3D_F tRes = Point3D_F::Transform(Point3D_F(tpoint.X - center.X, tpoint.Y - center.Y, tpoint.Z - center.Z), tMatr);
    tRes.X += center.X;
    tRes.Y += center.Y;
    tRes.Z += center.Z;
    return tRes;
}