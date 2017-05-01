#ifndef MESH_H
#define MESH_H

#include"MeshData.h"
//#include "general_tools.h"

class Mesh{

    public:

    Mesh(void);
    ~Mesh(void);

    void Read(std::string& mesh_fname_);
    void WriteMesh(std::string& write_fname_);
    void WriteMeshTecplot(std::string& write_fname_);
    void generate_meshData();

    MeshData* Release_meshData(void);

protected:
    void Initialize();
    void Reset();

    void compute_faceData();
    void compute_elemData();
    void compute_ghost_elemData();
    void compute_wallNodes();
    void compute_LSreconstruction_metrics();
    void compute_cell_volume_center(const int& cellID);
    void prepare_post_proc_elemlist(std::vector<int> &polygon_elemlist
                                    , const int& N_new_tri);

protected:
    MeshData *grid_data_=nullptr;

};

#endif
