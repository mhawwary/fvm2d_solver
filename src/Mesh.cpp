#include "Mesh.hpp"

Mesh::Mesh(void){
    Initialize();
}

Mesh::~Mesh(void){
    Reset();
}

void Mesh::Initialize(){
    Reset();
    grid_data_ = new MeshData;
}

void Mesh::Reset(){
    emptypointer(grid_data_);
    return;
}

void Mesh::Read(std::string &mesh_fname_){

    register int i;

    std::ifstream input (mesh_fname_);

    // Read in No. of nodes, faces , and elements
    //---------------------------------------------
    input >> grid_data_->Nnodes
          >> grid_data_->Nfaces
          >> grid_data_->Nelem;


    // Memory allocation:
    grid_data_->Xn = new double[grid_data_->Nnodes];
    grid_data_->Yn = new double[grid_data_->Nnodes];

    grid_data_->facelist= new Face[grid_data_->Nfaces];
    grid_data_->elemlist= new Elem[grid_data_->Nelem];

    // Fill in node coordinates list:
    //------------------------------------------------------
    for(i=0; i<grid_data_->Nnodes; i++){

        input >> grid_data_->Xn[i]
              >> grid_data_->Yn[i];
    }

    // Fill in some information for the dummy face list:
    //------------------------------------------------------
    for(i=0; i<grid_data_->Nfaces; i++){

        grid_data_->facelist[i].ID=i;

        input >> grid_data_->facelist[i].v0
              >> grid_data_->facelist[i].v1;

        grid_data_->facelist[i].v0 += -1;
        grid_data_->facelist[i].v1 += -1;
    }

    // Create a dummy element list that contains all elements:
    //--------------------------------------------------------

    for(i=0; i<grid_data_->Nelem; i++){
        grid_data_->elemlist[i].ID=i;
    }

    // Fill in more information for the dummy face list:
    //------------------------------------------------------
    grid_data_->Nbndfaces=0;

    for(i=0; i<grid_data_->Nfaces; i++){

        input >> grid_data_->facelist[i].Lcell
              >> grid_data_->facelist[i].Rcell;

        grid_data_->facelist[i].Lcell += -1;
        grid_data_->facelist[i].Rcell += -1;

        if(grid_data_->facelist[i].Rcell<0){

            grid_data_->Nbndfaces++;

            grid_data_->elemlist[grid_data_->facelist[i].Lcell].isbound=true;

            grid_data_->elemlist[grid_data_->facelist[i].Lcell].bnd_type
                                     = grid_data_->facelist[i].Rcell;

            grid_data_->facelist[i].isbound=true;
        }
    }

    // Calc. the no. of boundary elements:
    // ------------------------------------

    grid_data_->NbndElem =0;

    for(i=0; i<grid_data_->Nelem; i++)
        if(grid_data_->elemlist[i].isbound==true)
            grid_data_->NbndElem++;


    return;
}


/*void Mesh::generate_meshData(){

    register int i;

    //------------------------------------------------------------
    // Allocating arrays for boundary and interior entities
    //------------------------------------------------------------

    grid_data_->bnd_elemlist= new BoundElem[grid_data_->NbndElem];
    grid_data_->bnd_facelist= new BoundFace[grid_data_->Nbndfaces];

    grid_data_->NintElem = grid_data_->Nelem-grid_data_->NbndElem;
    grid_data_->int_elemlist= new Elem[grid_data_->NintElem];

    grid_data_->Nintfaces = grid_data_->Nfaces-grid_data_->Nbndfaces;
    grid_data_->int_facelist= new Face[grid_data_->Nintfaces];

    _print(grid_data_->Nelem,grid_data_->NbndElem);
    _(grid_data_->NintElem);
    _print(grid_data_->Nfaces,grid_data_->Nbndfaces);
    _(grid_data_->Nintfaces);

    //--------------------------------------------
    // Detecting interior and Boundary elements:
    //--------------------------------------------
    unsigned int ke=0,int_ke=0;
    for(i=0; i<grid_data_->Nelem; i++)
        if(grid_data_->elemlist[i].isbound==false){

            grid_data_->int_elemlist[int_ke].ID = grid_data_->elemlist[i].ID;

            grid_data_->elem_gid_to_intId[grid_data_->elemlist[i].ID] = int_ke;

            int_ke++; // interior Elem ID counter

        }else{

            grid_data_->bnd_elemlist[ke].ID = grid_data_->elemlist[i].ID;

            grid_data_->elem_gid_to_bid [grid_data_->elemlist[i].ID] = ke;

            ke++;  // boundary Elem ID counter
        }

    //--------------------------------------------
    // detecting boundary faces and elements:
    //--------------------------------------------
    unsigned int kf=0, int_kf=0;
    for(i=0; i<grid_data_->Nfaces; i++){

        // Boundary Face:
        if(grid_data_->facelist[i].Rcell<0){

            grid_data_->bnd_facelist[kf]=grid_data_->facelist[i];

            grid_data_->face_gid_to_bid[grid_data_->facelist[i].ID] = kf;

            if(grid_data_->facelist[i].Rcell==-1){

                //grid_data_->bnd_elemlist[ke].bnd_type=-1; // Wall boundary element

                ke=grid_data_->elem_gid_to_bid[grid_data_->facelist[i].Lcell];

                grid_data_->bnd_elemlist[ke].bnd_type=-1;
                grid_data_->bnd_facelist[kf].bnd_type=-1; // Wall boundary face

            }else {

                ke=grid_data_->elem_gid_to_bid[grid_data_->facelist[i].Lcell];
                grid_data_->bnd_elemlist[ke].bnd_type=-2; // FarFeild boundary element
                grid_data_->bnd_facelist[kf].bnd_type=-2; // FarFeild boundary face
            }

            kf++;  // boundary face ID counter

        } else {  // detecting interior faces:


            grid_data_->int_facelist[int_kf] = grid_data_->facelist[i];

            grid_data_->face_gid_to_intId[grid_data_->facelist[i].ID] = int_kf;

            int_kf++;  // interior face ID counter
        }
    }


    //grid_data_->Reset_dummy();

    compute_faceData();

    compute_elemData();


    return;

} */

void Mesh::generate_meshData(){

    //------------------------------------------------------------
    // Printing Some grid statistics
    //------------------------------------------------------------

    grid_data_->NintElem  = grid_data_->Nelem  - grid_data_->NbndElem;
    grid_data_->Nintfaces = grid_data_->Nfaces - grid_data_->Nbndfaces;

    _compare(grid_data_->NintElem,grid_data_->NbndElem);
    _(grid_data_->Nelem);
    _compare(grid_data_->Nintfaces,grid_data_->Nbndfaces);
    _(grid_data_->Nfaces);

    compute_faceData();

    compute_elemData();

    return;

}

void Mesh::compute_faceData(){

    register int i=0;

    double x0,x1,y0,y1,dx,dy,L;

    int v0,v1;

    for(i=0; i<grid_data_->Nfaces; i++){

        v0 = grid_data_->facelist[i].v0;
        v1 = grid_data_->facelist[i].v1;

        x0 = grid_data_->Xn[v0];
        y0 = grid_data_->Yn[v0];
        x1 = grid_data_->Xn[v1];
        y1 = grid_data_->Yn[v1];

        dx = x1-x0; dy = y1-y0;

        L = sqrt(pow(dx,2)+pow(dy,2));

        grid_data_->facelist[i].Af = L;

        grid_data_->facelist[i].nx =  dy/L;
        grid_data_->facelist[i].ny = -dx/L;

        grid_data_->facelist[i].Xf = 0.5*(x0+x1);
        grid_data_->facelist[i].Yf = 0.5*(y0+y1);
    }

    return;
}

void Mesh::compute_elemData(){


    std::set<int> *elem_to_face_set=nullptr;
    std::set<int> *elem_to_node_set=nullptr;

    std::set<int>::iterator elem_to_face_it;
    std::set<int>::iterator elem_to_node_it;

    elem_to_face_set = new std::set<int>[grid_data_->Nelem];
    elem_to_node_set = new std::set<int>[grid_data_->Nelem];

    register int i;

    int j,fID,iL,iR;

    for(i=0; i<grid_data_->Nfaces; i++){

        fID = grid_data_->facelist[i].ID;

        iL = grid_data_->facelist[i].Lcell;
        iR = grid_data_->facelist[i].Rcell;

        elem_to_face_set[iL].insert(fID);

        elem_to_node_set[iL].insert(grid_data_->facelist[i].v0);
        elem_to_node_set[iL].insert(grid_data_->facelist[i].v1);

        if(iR>=0){

            elem_to_face_set[iR].insert(fID);

            elem_to_node_set[iR].insert(grid_data_->facelist[i].v0);
            elem_to_node_set[iR].insert(grid_data_->facelist[i].v1);
        }
    }

    int Ntri=0,Nquad=0,Nothers=0;

    for(i=0; i<grid_data_->Nelem; i++){

        grid_data_->elemlist[i].n_local_faces = elem_to_face_set[i].size();
        grid_data_->elemlist[i].to_face = new int[grid_data_->elemlist[i].n_local_faces];
        grid_data_->elemlist[i].n_local_nodes = elem_to_face_set[i].size();
        grid_data_->elemlist[i].to_node = new int[grid_data_->elemlist[i].n_local_nodes];

        j=0;
        elem_to_node_it = elem_to_node_set[i].begin();

        for(elem_to_face_it=elem_to_face_set[i].begin();
            elem_to_face_it!=elem_to_face_set[i].end(); ++elem_to_face_it){

            grid_data_->elemlist[i].to_face[j] = *elem_to_face_it;

            grid_data_->elemlist[i].to_node[j] = *elem_to_node_it;

            j++; ++elem_to_node_it;
        }

        if(elem_to_face_set[i].size()==3){
            Ntri++;
        }else if(elem_to_face_set[i].size()==4) {
            Nquad++;
        }else if(elem_to_face_set[i].size()>=5){
            Nothers++;
        }else{
            FatalErrorST("\nProblem unidentified element type\n");
        }

        // Computing Cell Volume and Cell center position:

        compute_cell_volume_center(i); // needs centroid and volume verification
    }

    grid_data_->NtriElem= Ntri;
    grid_data_->NquadElem = Nquad;
    grid_data_->Npolygon = Nothers;

    /* Need to define a list of polygon elements (more than 4 sides) and divide them into triangles
     * for post processing.
     * */

    //printf("\nNtriangles: %d,  Nquads: %d,  Nother:  %d\n\n",Ntri,Nquad,Nothers);

    emptyarray(elem_to_face_set);
    emptyarray(elem_to_node_set);

    return;
}

void Mesh::compute_cell_volume_center(const int &ii){

    double Xf,Yf,nx,ny,A,Volume, xc,yc, test;

    int fID;

    register int i;

    Volume=0.0;
    xc =0.0;
    yc =0.0;

    test=0.0;

    for(i=0; i<grid_data_->elemlist[ii].n_local_faces; i++){

        fID = grid_data_->elemlist[ii].to_face[i] ;
        Xf = grid_data_->facelist[ fID ].Xf;
        Yf = grid_data_->facelist[ fID ].Yf;

        nx = grid_data_->facelist[ fID ].nx;
        ny = grid_data_->facelist[ fID ].ny;

        if(ii == grid_data_->facelist[fID].Rcell){
            nx = -nx;
            ny = -ny;
        }

        A = grid_data_->facelist[fID].Af;

        Volume += ( (Xf*nx) + (Yf*ny) ) * A;

        xc += ( (Xf*nx) + (Yf*ny) ) * A * Xf;
        yc += ( (Xf*nx) + (Yf*ny) ) * A * Yf;

        test += (nx+ny)*A;
    }

    if(abs(test)>0){
        FatalError("Area test Failed for some elements");
        std::cin.get();
    }

    Volume = 0.5 * Volume;

    grid_data_->elemlist[ii].Vc = Volume;

    xc = xc/(3.0 * Volume);
    yc = yc/(3.0 * Volume);

    grid_data_->elemlist[ii].Xc = xc;
    grid_data_->elemlist[ii].Yc = yc;

    /* Volume and Centroid calculations are derived from
     * "Improved Formulation for Geometric Properties of Arbitrary Polyhedra",
     * Z.J.Wang, AIAA Journal 1999
     * */

    return;
}

void Mesh::WriteMesh(std::string &write_fname_){

    std::ofstream output (write_fname_);

    output << grid_data_->Nnodes << " "
           << grid_data_->Nfaces << " "
           << grid_data_->Nelem  << "\n";

    register int i;

    for(i=0; i<grid_data_->Nnodes; i++){

        output << grid_data_->Xn[i] <<" "
               << grid_data_->Yn[i] <<"\n";
    }

    for(i=0; i<grid_data_->Nfaces; i++){

        output << grid_data_->facelist[i].ID << " "
               << grid_data_->facelist[i].v0 << " "
               << grid_data_->facelist[i].v1 << " "
               << grid_data_->facelist[i].Lcell << " "
               << grid_data_->facelist[i].bnd_type << "\n";
    }

    for(i=0; i<grid_data_->Nfaces; i++){

        output << grid_data_->facelist[i].ID << " "
               << grid_data_->facelist[i].v0 << " "
               << grid_data_->facelist[i].v1 << " "
               << grid_data_->facelist[i].Lcell << " "
               << grid_data_->facelist[i].Rcell << "\n";
    }

    for(i=0; i<grid_data_->Nelem; i++){

        output << grid_data_->elemlist[i].ID << " "
               << grid_data_->elemlist[i].bnd_type << "\n";
    }

    return;

}

MeshData* Mesh::Release_meshData(){

    MeshData* mesh_data_=grid_data_;

    grid_data_=nullptr;

    return mesh_data_;
}















