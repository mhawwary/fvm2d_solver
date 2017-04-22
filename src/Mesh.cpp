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

    if(grid_data_!=NULL)
        delete grid_data_;

    grid_data_ = NULL;

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
    for(i=0; i<grid_data_->Nfaces; i++){

        input >> grid_data_->facelist[i].Lcell
              >> grid_data_->facelist[i].Rcell;

        grid_data_->facelist[i].Lcell += -1;
        grid_data_->facelist[i].Rcell += -1;

        if(grid_data_->facelist[i].Rcell<0){

            grid_data_->Nbndfaces++;

            grid_data_->elemlist[grid_data_->facelist[i].Lcell].isbound=true;
        }
    }


    return;
}


void Mesh::generate_meshData(){

    register int i;

    // Calc. the no. of boundary elements:
    // ------------------------------------

    for(i=0; i<grid_data_->Nelem; i++)
        if(grid_data_->elemlist[i].isbound==true)
            grid_data_->NbndElem++;

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

}


void Mesh::compute_faceData(){

    register int i=0;

    double x0,x1,y0,y1,dx,dy,L;

    int v0,v1;

    for(i=0; i<grid_data_->Nbndfaces; i++){

        v0 = grid_data_->bnd_facelist[i].v0;
        v1 = grid_data_->bnd_facelist[i].v1;

        x0 = grid_data_->Xn[v0];
        y0 = grid_data_->Yn[v0];
        x1 = grid_data_->Xn[v1];
        y1 = grid_data_->Yn[v1];

        dx = x1-x0; dy = y1-y0;

        L = sqrt(pow(dx,2)+pow(dy,2));

        grid_data_->bnd_facelist[i].Af = L;

        grid_data_->bnd_facelist[i].nx =  dy/L;
        grid_data_->bnd_facelist[i].ny = -dx/L;

        grid_data_->bnd_facelist[i].Xf = 0.5(x0+x1);
        grid_data_->bnd_facelist[i].Yf = 0.5(y0+y1);
    }

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

    for(i=0; i<grid_data_->Nbndfaces; i++){

        output << grid_data_->bnd_facelist[i].ID << " "
               << grid_data_->bnd_facelist[i].v0 << " "
               << grid_data_->bnd_facelist[i].v1 << " "
               << grid_data_->bnd_facelist[i].Lcell << " "
               << grid_data_->bnd_facelist[i].bnd_type << "\n";
    }

    for(i=0; i<grid_data_->Nintfaces; i++){

        output << grid_data_->int_facelist[i].ID << " "
               << grid_data_->int_facelist[i].v0 << " "
               << grid_data_->int_facelist[i].v1 << " "
               << grid_data_->int_facelist[i].Lcell << " "
               << grid_data_->int_facelist[i].Rcell << "\n";
    }

    for(i=0; i<grid_data_->NbndElem; i++){

        output << grid_data_->bnd_elemlist[i].ID << " "
               << grid_data_->bnd_elemlist[i].bnd_type << "\n";
    }

    return;

}

MeshData* Mesh::Release_meshData(){

    MeshData* mesh_Data_=grid_data_;

    Reset();

    return mesh_Data_;
}















