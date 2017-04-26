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
    void compute_wallNodes();

    void compute_cell_volume_center(const int& cellID);

protected:
    MeshData *grid_data_=nullptr;

};

#endif
