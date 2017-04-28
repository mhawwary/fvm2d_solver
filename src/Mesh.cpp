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

    // Fill in more information for the dummy face list:
    //------------------------------------------------------
    grid_data_->Nbndfaces=0;

    grid_data_->Nwallnodes=0;  // No. of wall nodes = No. of wall faces

    for(i=0; i<grid_data_->Nfaces; i++){

        input >> grid_data_->facelist[i].Lcell
              >> grid_data_->facelist[i].Rcell;

        grid_data_->facelist[i].Lcell += -1;
        grid_data_->facelist[i].Rcell += -1;

        if(grid_data_->facelist[i].Rcell<0){

            grid_data_->facelist[i].bnd_type = grid_data_->facelist[i].Rcell;

            grid_data_->facelist[i].Rcell
                    = grid_data_->Nelem + grid_data_->Nbndfaces; // modifying the ID for ghost cells

            grid_data_->Nbndfaces++;

            grid_data_->facelist[i].isbound=true;

            if(grid_data_->facelist[i].bnd_type==-1)  grid_data_->Nwallnodes++;
        }
    }

    grid_data_->wall_nodelist=new int[grid_data_->Nwallnodes];

    // Constructing extended ElemList :
    //------------------------------------
    /* No. of ghost cells equal the no. of boundary faces */
    grid_data_->Nelem_extend = grid_data_->Nelem + grid_data_->Nbndfaces;
    grid_data_->elemlist= new Elem[grid_data_->Nelem_extend];

    for(i=0; i<grid_data_->Nelem_extend; i++)
        grid_data_->elemlist[i].ID=i;

    for(i=0; i<grid_data_->Nbndfaces; i++){

        grid_data_->elemlist[grid_data_->facelist[i].Lcell].isbound=true;

        grid_data_->elemlist[grid_data_->facelist[i].Lcell].bnd_type
                = grid_data_->facelist[i].bnd_type;

        grid_data_->elemlist[grid_data_->facelist[i].Rcell].isbound=true;

        grid_data_->elemlist[grid_data_->facelist[i].Rcell].bnd_type
                = grid_data_->facelist[i].bnd_type;
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

    compute_wallNodes();

    compute_LSreconstruction_metrics();

    return;

}

void Mesh::compute_wallNodes(){

    std::set<int> wall_node_set;

    register int i;

    for(i=0; i<grid_data_->Nbndfaces; i++){

        if(grid_data_->facelist[i].bnd_type==-1){
            wall_node_set.insert(grid_data_->facelist[i].v0);
            wall_node_set.insert(grid_data_->facelist[i].v1);
        }
    }

    std::vector<int> wall_node_vector (wall_node_set.begin(), wall_node_set.end());

    std::sort(wall_node_vector.begin(),wall_node_vector.end());

    int Nupper=0,Nlower=0;
    std::vector<int> upperwall_node_vector;
    std::vector<int> lowerwall_node_vector;

    double Y_;

    for(i=0; i<grid_data_->Nwallnodes; i++){
        grid_data_->wall_nodelist[i] = wall_node_vector[i];
        Y_ = grid_data_->Yn[grid_data_->wall_nodelist[i]];

        if(Y_ > 0){ // upper wall point
            Nupper++;
            upperwall_node_vector.push_back(wall_node_vector[i]);

        }else if(Y_<0) {// lower wall point
            Nlower++;
            lowerwall_node_vector.push_back(wall_node_vector[i]);

        }else{  // Ywn==0
            Nupper++; Nlower++;
            upperwall_node_vector.push_back(wall_node_vector[i]);
            lowerwall_node_vector.push_back(wall_node_vector[i]);
        }
    }

    //Nupper++;  // for adding the point 1,0;
    //upperwall_node_vector.push_back(0);

    grid_data_->NupperWallnodes = Nupper;
    grid_data_->NlowerWallnodes = Nlower;

    double *Xupp=nullptr,*Yupp=nullptr;
    Xupp = new double[Nupper];
    Yupp = new double[Nupper];

    grid_data_->upper_wall_nodelist = new int[Nupper];

    for(i=0; i<Nupper; i++){
        grid_data_->upper_wall_nodelist[i] = upperwall_node_vector[i];
        Xupp[i] = grid_data_->Xn[upperwall_node_vector[i]];
        Yupp[i] = grid_data_->Yn[upperwall_node_vector[i]];
    }

    double *Xlow=nullptr,*Ylow=nullptr;
    Xlow = new double[Nlower];
    Ylow = new double[Nlower];

    grid_data_->lower_wall_nodelist = new int[Nlower];

    for(i=0; i<Nlower; i++){
        grid_data_->lower_wall_nodelist[i] = lowerwall_node_vector[i];
        Xlow[i] = grid_data_->Xn[lowerwall_node_vector[i]];
        Ylow[i] = grid_data_->Yn[lowerwall_node_vector[i]];
    }

    QuickSort3(Xlow,Ylow,grid_data_->lower_wall_nodelist,0,Nlower-1);
    QuickSort3(Xupp,Yupp,grid_data_->upper_wall_nodelist,0,Nupper-1);

//    printf("\n Printing upper Wall node IDs: \n");
//    for(i=0; i<Nupper; i++){
//        printf("\n%d %e %e\n",grid_data_->upper_wall_nodelist[i],Xupp[i],Yupp[i]);
//    }

//    printf("\n Printing lower Wall node IDs: \n");
//    for(i=0; i<Nlower; i++){
//        printf("\n%d %e %e\n",grid_data_->lower_wall_nodelist[i],Xlow[i],Ylow[i]);
//    }

    wall_node_set.clear();

    for(i=Nupper-1; i>=0; i--){
        wall_node_set.insert(grid_data_->upper_wall_nodelist[i]);
    }

    for(i=0; i<Nlower; i++){
        wall_node_set.insert(grid_data_->lower_wall_nodelist[i]);
    }

    std::set<int>::iterator wall_node_set_it;

    i=0;
    for(wall_node_set_it=wall_node_set.begin();
        wall_node_set_it!=wall_node_set.end();
        ++wall_node_set_it) {

        grid_data_->wall_nodelist[i] = *wall_node_set_it;

        i++;
    }

    std::sort(grid_data_->wall_nodelist
              ,grid_data_->wall_nodelist+grid_data_->Nwallnodes
              ,std::greater<int>());

    emptyarray(Xupp);
    emptyarray(Yupp);

    emptyarray(Xlow);
    emptyarray(Ylow);

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
    std::set<int> *node_to_elem_set=nullptr;

    std::set<int>::iterator elem_to_face_it;
    std::set<int>::iterator elem_to_node_it;
    std::set<int>::iterator node_to_elem_it;

    elem_to_face_set = new std::set<int>[grid_data_->Nelem];
    elem_to_node_set = new std::set<int>[grid_data_->Nelem];

    node_to_elem_set = new std::set<int>[grid_data_->Nnodes];

    register int i;

    int j,fID,iL,iR;

    for(i=0; i<grid_data_->Nfaces; i++){

        fID = i;

        iL = grid_data_->facelist[i].Lcell;
        iR = grid_data_->facelist[i].Rcell;

        elem_to_face_set[iL].insert(fID);

        elem_to_node_set[iL].insert(grid_data_->facelist[i].v0);
        elem_to_node_set[iL].insert(grid_data_->facelist[i].v1);

        node_to_elem_set[grid_data_->facelist[i].v0].insert(iL);
        node_to_elem_set[grid_data_->facelist[i].v1].insert(iL);

        node_to_elem_set[grid_data_->facelist[i].v0].insert(iR);
        node_to_elem_set[grid_data_->facelist[i].v1].insert(iR);

        if(grid_data_->facelist[i].bnd_type==0){

            elem_to_face_set[iR].insert(fID);

            elem_to_node_set[iR].insert(grid_data_->facelist[i].v0);
            elem_to_node_set[iR].insert(grid_data_->facelist[i].v1);
        }
    }

    grid_data_->node_to_elemlist = new int*[grid_data_->Nnodes];

    grid_data_->Nnode_neighElem = new int [grid_data_->Nnodes];

    for(i=0; i<grid_data_->Nnodes; i++){

        grid_data_->Nnode_neighElem[i] = node_to_elem_set[i].size();
        grid_data_->node_to_elemlist[i] = new int[grid_data_->Nnode_neighElem[i] ];

        j=0;
        for(node_to_elem_it=node_to_elem_set[i].begin();
            node_to_elem_it!=node_to_elem_set[i].end();
            ++node_to_elem_it){

            grid_data_->node_to_elemlist[i][j] = *node_to_elem_it;
            j++;
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

//        for(elem_to_face_it=elem_to_face_set[i].begin();
//            elem_to_face_it!=elem_to_face_set[i].end(); ++elem_to_face_it){

//            grid_data_->elemlist[i].to_face[j] = *elem_to_face_it;

//            grid_data_->elemlist[i].to_node[j] = *elem_to_node_it;

//            j++; ++elem_to_node_it;
//        }

        elem_to_face_it=elem_to_face_set[i].begin();

        for(j=0; j<grid_data_->elemlist[i].n_local_faces; j++){

            grid_data_->elemlist[i].to_face[j] = *elem_to_face_it;

            grid_data_->elemlist[i].to_node[j] = *elem_to_node_it;

            ++elem_to_node_it; elem_to_face_it++;
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

    printf("\nNtriangles: %d,  Nquads: %d,  Nother:  %d\n",Ntri,Nquad,Nothers);

    // Calculation of Ghost cells geometric data:
    //---------------------------------------------
    double Xf,Yf,Xc,Yc,nx,ny;

    for(j=0; j<grid_data_->Nbndfaces; j++){
        iL = grid_data_->facelist[j].Lcell;
        iR = grid_data_->facelist[j].Rcell;
        Xf = grid_data_->facelist[j].Xf;
        Yf = grid_data_->facelist[j].Yf;
        nx = grid_data_->facelist[j].nx;
        ny = grid_data_->facelist[j].ny;

        Xc = grid_data_->elemlist[iL].Xc;
        Yc = grid_data_->elemlist[iL].Yc;

        grid_data_->elemlist[iR].Vc = grid_data_->elemlist[iL].Vc;

        grid_data_->elemlist[iR].Xc = Xc*((ny*ny)-(nx*nx)) -2*Yc*nx*ny - Yf*nx + Xf *(1.0+ny) ;
        grid_data_->elemlist[iR].Yc = Yc*((nx*nx)-(ny*ny)) -2*Xc*nx*ny - Xf*nx + Yf *(1.0-ny) ;
    }

    emptyarray(elem_to_face_set);
    emptyarray(elem_to_node_set);
    emptyarray(node_to_elem_set);

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

    if(abs(test)!=0){
        FatalError("Area test Failed for some elements");
        printf("Area test: %e\n",test);
        std::cin.get();
    }

    Volume = 0.5 * Volume;

    if(Volume < 0) FatalErrorST("Negative Volume");

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

void Mesh::compute_LSreconstruction_metrics(){

    register int j;

    int i=0,fID=0,iL,iR,ii;

    double Xj=0.0,Xn=0.0,Yj=0.0,Yn=0.0;

    grid_data_->Ixx = new double[grid_data_->Nelem];
    grid_data_->Iyy = new double[grid_data_->Nelem];
    grid_data_->Ixy = new double[grid_data_->Nelem];
    double IXY = 0.0;

    for(j=0; j<grid_data_->Nelem; j++){

        grid_data_->Ixx[j] = 0.0;
        grid_data_->Iyy[j] = 0.0;
        grid_data_->Ixy[j] = 0.0;
        IXY = 0.0;

        Xj = grid_data_->elemlist[j].Xc;
        Yj = grid_data_->elemlist[j].Yc;

        grid_data_->elemlist[j].DX = new double[grid_data_->elemlist[j].n_local_faces];
        grid_data_->elemlist[j].DY = new double[grid_data_->elemlist[j].n_local_faces];

        for(i=0; i<grid_data_->elemlist[j].n_local_faces; i++){

            fID = grid_data_->elemlist[j].to_face[i];
            iL = grid_data_->facelist[fID].Lcell;
            iR = grid_data_->facelist[fID].Rcell;

            if(iL!=j) ii=iL;
            else ii=iR;

            Xn = grid_data_->elemlist[ii].Xc;
            Yn = grid_data_->elemlist[ii].Yc;

            grid_data_->Ixx[j] += pow((Xn-Xj),2);
            grid_data_->Iyy[j] += pow((Yn-Yj),2);
            grid_data_->Ixy[j] += ( (Xn-Xj) * (Yn-Yj) );

            grid_data_->elemlist[j].DX[i] = (grid_data_->elemlist[ii].Xc - grid_data_->elemlist[j].Xc);
            grid_data_->elemlist[j].DY[i] = (grid_data_->elemlist[ii].Yc - grid_data_->elemlist[j].Yc);
        }

        IXY = 1.0/( (grid_data_->Ixx[j]*grid_data_->Iyy[j]) - pow(grid_data_->Ixy[j],2) ) ;
        grid_data_->Ixx[j] = IXY * grid_data_->Ixx[j];
        grid_data_->Iyy[j] = IXY * grid_data_->Iyy[j];
        grid_data_->Ixy[j] = IXY * grid_data_->Ixy[j];
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

    for(i=0; i<grid_data_->Nfaces; i++){

        output << grid_data_->facelist[i].ID << " "
               << grid_data_->facelist[i].v0 << " "
               << grid_data_->facelist[i].v1 << " "
               << grid_data_->facelist[i].Lcell << " "
               << grid_data_->facelist[i].Rcell << " "
               << grid_data_->facelist[i].bnd_type << "\n";
    }

    for(i=0; i<grid_data_->Nelem; i++){

        output << grid_data_->elemlist[i].ID << " "
               << grid_data_->elemlist[i].bnd_type << "\n";
    }

    std::ofstream output1 ("./post_process/CellVolume.dat");

    for(i=0; i<grid_data_->Nelem; i++){

        output1 << grid_data_->elemlist[i].Xc << " "
               << grid_data_->elemlist[i].Yc << " "
               << grid_data_->elemlist[i].Vc << "\n";
    }

    return;
}

void Mesh::WriteMeshTecplot(std::string &write_fname_){

    register int k,j;

    const char* fname = write_fname_.c_str();

    //fname = write_fname_.c_str();

    FILE*  outfile=fopen(fname,"wt");

//    fprintf(outfile, "VARIABLES = \"X\",\"Y\",\"M\",\"P\",\"RHO\",\"u\",\"v\"");
    fprintf(outfile, "VARIABLES = \"X\",\"Y\",\"M\"");
    fprintf(outfile, "\nZONE N=%d, E=%d, F=FEPOINT, ET=FEPOLYGON\n", grid_data_->Nnodes, grid_data_->Nelem);

    for(k=0; k<grid_data_->Nnodes; k++)
    {
        fprintf(outfile, "%e %e %e\n",
                grid_data_->Xn[k], grid_data_->Yn[k], 0.0);
    }

    int node_id;

    for(k=0; k<grid_data_->Nelem; k++)
    {
        fprintf(outfile, "\n");
        for(j=0; j<grid_data_->elemlist[k].n_local_nodes; j++){
            node_id = grid_data_->elemlist[k].to_node[j];
            fprintf(outfile, "%d ",node_id+1);
        }
    }

    fclose(outfile);

    //emptyarray(fname);

    return;
}

MeshData* Mesh::Release_meshData(){

    MeshData* mesh_data_=grid_data_;

    grid_data_=nullptr;

    return mesh_data_;
}















