#include"MeshData.h"
#include "general_tools.h"

class Mesh{

    public:

    Mesh(void);
    ~Mesh(void);

    void Read(std::string& mesh_fname_);
    void WriteMesh(std::string& write_fname_);
    void generate_meshData();

    MeshData* ReleaseMeshData(void);

protected:
    void Initialize();
    void Reset();


protected:
    MeshData *_grid_data_=NULL;

};
