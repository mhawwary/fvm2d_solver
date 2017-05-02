#ifndef MESHDATA_H
#define MESHDATA_H

#include "general_tools.h"

struct Face{

    int ID=0;

    bool isbound=false;  // // true if it is a boundary face

    int v0=0;
    int v1=1;

    int Lcell=0;   // Left  cell ID
    int Rcell=0;   // Right cell ID

    double nx=0.0;  // normal y component
    double ny=0.0;  // normal x component

    double Xf=0.0;  // face center x
    double Yf=0.0;  // face center y

    double Af=1.0;  // face area

    int bnd_type=0; // 0:interior , -1: Wall , -2: FarFeild

/*    Face& copy_face(const Face& RFace){

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
//    } */

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

        bnd_type = RFace.bnd_type;

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

    ~Elem(){

        emptyarray(to_face);
        emptyarray(to_node);
        emptyarray(DX);
        emptyarray(DY);
    }

    bool isbound=false;  // true if it is a boundary elem
    int bnd_type=0; // 0:interior , -1: Wall , -2: FarFeild

    int ID=0;   // Elem ID

    int *to_node=nullptr;   // elem_to_node list
    int *to_face=nullptr;   // elem_to_face list

    int n_local_faces = 4; // no. of elem faces
    int n_local_nodes = 4; // no. of elem nodes

    int Nneigh_elem=1;

    //int *to_neigh_elem=nullptr; // neighbour elements ID

    double Vc=1.0;  // cell volume

    double Xc=0.0; // cell center xc
    double Yc=0.0; // cell center yc

    double *DX=nullptr;
    double *DY=nullptr;

    virtual Elem& operator =(const Elem& RElem){

        isbound= RElem.isbound;
        bnd_type = RElem.bnd_type;

        ID = RElem.ID;

        n_local_faces = RElem.n_local_faces;
        n_local_nodes = RElem.n_local_nodes;

        to_face = new int[n_local_faces];
        to_node = new int[n_local_nodes];

        int i=0;

        for(i=0; i<n_local_faces; i++) to_face[i] = RElem.to_face[i];
        for(i=0; i<n_local_nodes; i++) to_node[i] = RElem.to_node[i];

        Nneigh_elem = RElem.Nneigh_elem;

        Vc = RElem.Vc;
        Xc = RElem.Xc;
        Yc = RElem.Yc;

        //DX = RElem.DX;
        //DY = RElem.DY;

        return *this;
    }

};

struct BoundElem: Elem{

    int bnd_type=-1; // -1: Wall , -2: FarFeild

};

struct GhostElem{

    //~GhostElem(){}

    int bnd_type=0; // 0:interior , -1: Wall , -2: FarFeild

    int ID=0;   // GhostElem ID

    double Vc=1.;  // cell volume

    double Xc=0.; // cell center xc
    double Yc=0.; // cell center yc

};

struct MeshData{

    ~MeshData(){

        emptyarray(elemlist);
        emptyarray(facelist);
        emptyarray(post_proc_elemlist);
        emptyarray(polygon_elem_origID);

        emptyarray(wall_nodelist);
        emptyarray(upper_wall_nodelist);
        emptyarray(lower_wall_nodelist);

        emptyarray(Nnode_neighElem);
        emptyarray(Nnodes,node_to_elemlist);

        emptyarray(Xn);
        emptyarray(Yn);
        emptyarray(post_Xn);
        emptyarray(post_Yn);

        emptyarray(Ixx);
        emptyarray(Iyy);
        emptyarray(Ixy);

    }

    int Nnodes=0;
    int Nfaces=0;
    int Nelem =0;
    int Nelem_extend=0;

    int NbndNodes=0;
    int Nbndfaces=0;
    int NbndElem=0;

    int Nwallnodes=0;
    int NupperWallnodes=0;
    int NlowerWallnodes=0;

    int NintElem=0;
    int Nintfaces=0;
    int NintNodes=0;

    int NtriElem=0;
    int NquadElem=0;
    int Npolygon=0;
    int NpostProc=0;
    int Nnodes_postproc=0;

    double *Xn=nullptr; // node x coord
    double *Yn=nullptr; // node y coord

    Face *facelist=nullptr;
    Elem *elemlist=nullptr;

    int *wall_nodelist=nullptr;
    int *upper_wall_nodelist=nullptr;
    int *lower_wall_nodelist=nullptr;

    int **node_to_elemlist=nullptr;
    int **node_to_facelist=nullptr;

    int *Nnode_neighElem=nullptr;

    // Least Squares Reconstruction Data:
    double *Ixx = nullptr;
    double *Iyy = nullptr;
    double *Ixy = nullptr;

    // PostProcessing ElementList all Triangles
    // and Quardilaterals with 4 sides only:

    Elem* post_proc_elemlist=nullptr;

    double* post_Xn=nullptr;
    double* post_Yn=nullptr;

    int *polygon_elem_origID=nullptr;

    std::map <int, int> uppwall_node_gid_to_bid;
    std::map <int, int> lowwall_node_gid_to_bid;
    std::map <int, int> wall_node_gid_to_bid;

//    std::map <int, int> elem_gid_to_bid;  // map from global ID to bound ID
//    std::map <int, int> elem_gid_to_intId;  // map from global ID to interior ID

//    std::map <int, int> face_gid_to_bid;  // map from global ID to bound ID
//    std::map <int, int> face_gid_to_intId;  // map from global ID to interior ID

};

#endif

