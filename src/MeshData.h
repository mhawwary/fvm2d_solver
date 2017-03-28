

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

    double Af=1;  // face area
};

struct Elem{

    bool isbound=false;  // true if it is a boundary elem

    int ID=0;   // Elem ID

    int Nn=4;  // no. of nodes
    int Nf=4;  // no. of faces

    int *nID=nullptr;   // node ID
    int *fID=nullptr;   // face ID

    int *neigh_eID=nullptr; // neighbour elements ID

    double Ve=1;  // cell volume

    double Xc=0.; // cell center xc
    double Yc=0.; // cell center yc
};

struct Boundary{

    unsigned int *gID=nullptr;  // global ID
    int *bnd_type=nullptr;    // flag for boundary types
};

struct MeshData{

    int Nnodes=0;
    int Nfaces=0;
    int Nelem =0;

    int NbndNodes=0;
    int Nbndfaces=0;
    int NbndElem=0;

    double *Xn=nullptr; // node x coord
    double *Yn=nullptr; // node y coord

    Face *facelist=nullptr;  // array of faces
    Elem *elemlist=nullptr;  // array of elements

    int* bNodelist=nullptr;     // array of all boundary nodes
    Face* bFacelist=nullptr;     // array of all boundary nodes
    Elem* bElemlist=nullptr;     // array of all boundary nodes


};
