#include "general_tools.h"

struct Face{

    int ID=0;

    //bool isbound=false;  // // true if it is a boundary face

    int v0=0;
    int v1=1;

    int Lcell=0;   // Left  cell ID
    int Rcell=0;   // Right cell ID

    double nx=0.0;  // normal y component
    double ny=0.0;  // normal x component

    double Xf=0.0;  // face center x
    double Yf=0.0;  // face center y

    double Af=1;  // face area

//    Face& copy_face(const Face& RFace){

//        ID = RFace.ID;

//        Af = RFace.Af;

//        v0 = RFace.v0;
//        v1 = RFace.v1;

//        Lcell = RFace.Lcell;
//        Rcell = RFace.Rcell;

//        Xf = RFace.Xf;
//        Yf = RFace.Yf;

//        nx = RFace.nx;
//        ny = RFace.ny;
//    }

    virtual Face& operator =(const Face& RFace){

        ID = RFace.ID;

        Af = RFace.Af;

        v0 = RFace.v0;
        v1 = RFace.v1;

        Lcell = RFace.Lcell;
        Rcell = RFace.Rcell;

        Xf = RFace.Xf;
        Yf = RFace.Yf;

        nx = RFace.nx;
        ny = RFace.ny;

        return *this;

    }
};

struct BoundFace: Face{

    int bnd_type=-1; // -1: Wall , -2: FarFeild

    virtual BoundFace& operator =(const Face& RFace){

        ID = RFace.ID;

        Af = RFace.Af;

        v0 = RFace.v0;
        v1 = RFace.v1;

        Lcell = RFace.Lcell;
        Rcell = RFace.Rcell;

        Xf = RFace.Xf;
        Yf = RFace.Yf;

        nx = RFace.nx;
        ny = RFace.ny;

        return *this;
    }

    virtual BoundFace& operator =(const BoundFace& RFace){

        ID = RFace.ID;

        Af = RFace.Af;

        v0 = RFace.v0;
        v1 = RFace.v1;

        Lcell = RFace.Lcell;
        Rcell = RFace.Rcell;

        Xf = RFace.Xf;
        Yf = RFace.Yf;

        nx = RFace.nx;
        ny = RFace.ny;

        bnd_type = RFace.bnd_type;

        return *this;
    }

};

struct Elem{

    bool isbound=false;  // true if it is a boundary elem

    int ID=0;   // Elem ID

    int Nn=4;  // no. of nodes
    int Nf=4;  // no. of faces

    int *nID=nullptr;   // node ID
    int *fID=nullptr;   // face ID

    int *neigh_eID=nullptr; // neighbour elements ID

    double Vc=1;  // cell volume

    double Xc=0.; // cell center xc
    double Yc=0.; // cell center yc

};

struct BoundElem: Elem{

    int bnd_type=-1; // -1: Wall , -2: FarFeild

};

struct GhostElem{


};


struct MeshData{

    int Nnodes=0;
    int Nfaces=0;
    int Nelem =0;

    int NbndNodes=0;
    int Nbndfaces=0;
    int NbndElem=0;

    int NintElem=0;
    int Nintfaces=0;
    int NintNodes=0;

    double *Xn=nullptr; // node x coord
    double *Yn=nullptr; // node y coord

    Face *dummy_facelist=nullptr;
    Elem *dummy_elemlist=nullptr;

    Face *int_facelist=nullptr;  // array of faces
    Elem *int_elemlist=nullptr;  // array of elements

    BoundElem *bnd_elemlist=nullptr;
    BoundFace *bnd_facelist=nullptr;

    GhostElem *gh_elemlist=nullptr;

    std::map <int, int> elem_gid_to_bid;  // map from global ID to bound ID
    std::map <int, int> elem_gid_to_intid;  // map from global ID to interior ID

    void Reset_dummy(){

        emptyarray(dummy_facelist);
        emptyarray(dummy_elemlist);

    }

};

