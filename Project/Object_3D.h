#ifndef OBJECT_3D_H
#define OBJECT_3D_H

#include "AllTypes.h"
#include "BasicClasses.h"
#include <list>

class Face;
class Vertex;
class Mesh;

// face of 3D mesh
class Face
{
public:

    Face() { };
    ~Face() { };

    Face(Face& toriginal);
    
    Face(Face& toriginal, const float tscale);

    int A;
    int B;
    int C;
    Point3D_F center;
    Point3D_F faceNormal;

};

// vertex of 3D mesh
class Vertex
{
public:

    ~Vertex() { };
    Vertex() { };

    Vertex(Vertex& toriginal);

    Vertex(Vertex& toriginal, const float tscale);

    Point3D_F coord;
    Point3D_F norm;
    Point3D_F wCoord; // wCoord is needed for visualizing the mesh. It is not important for analysis
};

class Mesh
{
public:

    Mesh(char* tname, int tverticesCount, int tfacesCount);

    Mesh(char* tfileName, float tscale);

    Mesh(ifstream tbr);

    Mesh(std::vector<shared_ptr<Vertex>>& tvertices, std::vector<shared_ptr<Face>>& tfaces);

    Mesh(std::vector<shared_ptr<Point3D_F>>& tvertices, std::vector<std::vector<int>>& tfaces, float tscale);

    static shared_ptr<Mesh> ReadSTLMeshFlip(const char* tfileName, float newScale);

    void SaveAsSTL(const char* tfileName/*a file from which the header will be copied*/, const char* tfileResult);

    __declspec ( property ( put = SetPosition, get = GetPosition ) ) 
                Point3D_F Position ;

    __declspec ( property ( put = SetRotation, get = GetRotation ) ) 
                Point3D_F Rotation ;

    __declspec ( property ( get = GetCenter) ) 
                Point3D_F Center ;

    Point3D_F GetPosition() { return position; }

    void SetPosition(Point3D_F tvalue) { position = Point3D_F(tvalue.X, tvalue.Y, tvalue.Z); }

    Point3D_F GetRotation() { return rotation; }

    void SetRotation(Point3D_F tvalue) { rotation = Point3D_F(tvalue.X, tvalue.Y, tvalue.Z); }

    Point3D_F GetCenter() { return center; }

    std::vector<int> FOV(int i);

    void SetNewVerteces(std::vector<Point3D_F>& tvertices);

    void PositionObject(float tx_g, float ty_g, float tz_g);

    void MoveObject(float tx, float ty, float tz);

    void ScaleIndividualObject(float dx, float dy, float dz);

    // scaling, when mesh is positioned at the center of coordinate system : center = (0, 0, 0)
    void ScaleBrutally(float dx, float dy, float dz);

    void Rotate(float dx, float dy, float dz);

    void RecomputeFaceCenters();

    void RecomputeNormals();

    std::vector<shared_ptr<Vertex>> Vertices;
    std::vector<shared_ptr<Face>> Faces;

private:
    
    //the next three functions are needed during reading of stl file. The key idea is to remove duplicates
    static void ReNumFaces(std::vector<std::vector<int>>& tfaces, std::vector<int>& tmarks, std::vector<Point4D_F>& tvertices);

    static std::vector<shared_ptr<Point3D_F>> ExtractCorectVertices(std::vector<int>& tmarks, std::vector<Point4D_F>& tvertices, int tlen);

    static int MarkToRemove(std::vector<Point4D_F>& tvertices, std::vector<int>& tmarks);

    static void FlipFacesAndVertices(std::vector<std::vector<int>>& tfaces, std::vector<shared_ptr<Point3D_F>>& tverticesCorr);

    void ComputeCenter();

    static Point3D_F ComputeFaceCenter(Vertex a, Vertex b, Vertex c);

    Point3D_F VectorOffset(Point3D_F pIn, Point3D_F pOffset);

    Point3D_F VectorGetNormal(Point3D_F a, Point3D_F b);

    Point3D_F VectorNormalize(Point3D_F pIn);

    Point3D_F ComputeFNormal(Point3D_F p1, Point3D_F p2, Point3D_F p3);

    void ComputeFaceNormals();

    void AssignVertecesToFaces();

    Point3D_F ComputeVNormal(int tnum);

    void ComputeVertexNormals();

    Point3D_F RotateOneVertex(Point3D_F tpoint, double** tMatr);


    char* name;
    Point3D_F position;
    Point3D_F rotation;

    // this variable is needed to have easy access to neighboring faces and vertices
    std::vector<std::vector<int>> facesOfVerteces;
    Point3D_F center;
};

#endif