//////////////////////////////////////////////////////////////////////
// Uses the ASSIMP library to read mesh models in one of 30+ file
// types into a structure suitable for the raytracer.
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
namespace fs = std::filesystem;

#include "geom.h"
#include "raytrace.h"

#ifndef NOASSIMP
#include <assimp/Importer.hpp>
#include <assimp/version.h>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

aiColor3D ones(1,1,1);

void recurseModelNodes(Scene* scene,
                       const  aiScene* aiscene,
                       const  aiNode* node,
                       const aiMatrix4x4& parentTr,
                       std::vector<Material*>,
                       const int level=0);

void Scene::ReadAssimpFile(const std::string& path, const mat4& M)
{
    aiMatrix4x4 modelTr(M[0][0], M[1][0], M[2][0], M[3][0],
                        M[0][1], M[1][1], M[2][1], M[3][1],
                        M[0][2], M[1][2], M[2][2], M[3][2],
                        M[0][3], M[1][3], M[2][3], M[3][3]);

    // Does the file exist?
    std::ifstream find_it(path.c_str());
    if (find_it.fail()) {
        std::cerr << "File not found: "  << path << std::endl;
        exit(-1); }

    // Invoke assimp to read the file.
    printf("Assimp %d.%d Reading %s\n", aiGetVersionMajor(), aiGetVersionMinor(), path.c_str());
    Assimp::Importer importer;
    const aiScene* aiscene = importer.ReadFile(path.c_str(),
                                               aiProcess_Triangulate|aiProcess_GenSmoothNormals);
    
    if (!aiscene) {
        printf("... Failed to read.\n");
        exit(-1); }

    if (!aiscene->mRootNode) {
        printf("Scene has no rootnode.\n");
        exit(-1); }
    
    //aiProcessPreset_TargetRealtime_Quality);
    printf("mNumCameras: %d\n", aiscene->mNumCameras);
    printf("mNumLights: %d\n", aiscene->mNumLights);
    printf("mNumMaterials: %d\n", aiscene->mNumMaterials);
    printf("mNumMeshes: %d\n", aiscene->mNumMeshes);
    printf("mNumTextures: %d\n", aiscene->mNumTextures);

    std::vector<Material*> materials;

#ifdef IGNORE_MATERIALS
    for (int i=0;  i<aiscene->mNumMaterials;  i++) {
        materials.push_back(currentMat); }
#else
    for (int i=0;  i<aiscene->mNumMaterials;  i++) {
        aiMaterial* mtl = aiscene->mMaterials[i];
        aiString name;
        mtl->Get(AI_MATKEY_NAME, name);
        aiColor3D emit(0.f,0.f,0.f); 
        aiColor3D diff(0.f,0.f,0.f), spec(0.f,0.f,0.f); 
        float alpha = 20.0;
        bool he = mtl->Get(AI_MATKEY_COLOR_EMISSIVE, emit);
        bool hd = mtl->Get(AI_MATKEY_COLOR_DIFFUSE, diff);
        bool hs = mtl->Get(AI_MATKEY_COLOR_SPECULAR, spec);
        bool ha = mtl->Get(AI_MATKEY_SHININESS, &alpha, NULL);
        Material* newmat;
        
        if (!emit.IsBlack() && diff.IsBlack()) { // An emitter
            newmat = new Light(vec3(emit.r, emit.g, emit.b)); }
        
        else if (diff==ones && spec==ones) { // The default created by assimp in none exist
            newmat = currentMat; }

        else { 
            vec3 Kd(0.5f, 0.5f, 0.5f); 
            vec3 Ks(0.03f, 0.03f, 0.03f);
            if (AI_SUCCESS == hd) Kd = vec3(diff.r, diff.g, diff.b);
            if (AI_SUCCESS == hs) Ks = vec3(spec.r, spec.g, spec.b);
            newmat = new Material(Kd, Ks, alpha, vec3(0,0,0), 1);
            
            aiString texPath;
            if (AI_SUCCESS == mtl->GetTexture(aiTextureType_DIFFUSE, 0, &texPath)) {
                fs::path fullPath = path;
                fullPath.replace_filename(texPath.C_Str());
                newmat->setTexture(fullPath.string()); } }
    
        materials.push_back(newmat); }
#endif
    recurseModelNodes(this, aiscene, aiscene->mRootNode, modelTr, materials);
}

// Recursively traverses the assimp node hierarchy, accumulating
// modeling transformations, and creating and transforming any meshes
// found.  Meshes comming from assimp can have associated surface
// properties, so each mesh *copies* the current BRDF as a starting
// point and modifies it from the assimp data structure.
void recurseModelNodes(Scene* scene,
                       const aiScene* aiscene,
                       const aiNode* node,
                       const aiMatrix4x4& parentTr,
                       std::vector<Material*> materials,
                       const int level)
{
    // Print line with indentation to show structure of the model node hierarchy.
    //for (int i=0;  i<level;  i++) printf("| ");
    //printf("%s \n", node->mName.data);

    // Accumulating transformations while traversing down the hierarchy.
    aiMatrix4x4 childTr = parentTr*node->mTransformation;
    aiMatrix3x3 normalTr = aiMatrix3x3(childTr); // Really should be inverse-transpose for full generality
     
    // Loop through this node's meshes
    for (unsigned int m=0;  m<node->mNumMeshes; ++m) {
        aiMesh* aimesh = aiscene->mMeshes[node->mMeshes[m]];
        //printf("  %d: %d:%d\n", m, aimesh->mNumVertices, aimesh->mNumFaces);

        // Arrays to hold all vertex and triangle data.
        MeshData* meshdata = new MeshData;
        
        // Loop through all vertices and record the
        // vertex/normal/texture/tangent data with the node's model
        // transformation applied.
        //printf("  mNumVertices: %d\n", aimesh->mNumVertices);
        for (unsigned int t=0;  t<aimesh->mNumVertices;  ++t) {
            aiVector3D aipnt = childTr*aimesh->mVertices[t];
            aiVector3D ainrm = aimesh->HasNormals() ? normalTr*aimesh->mNormals[t] : aiVector3D(0,0,1);
            aiVector3D aitex = aimesh->HasTextureCoords(0) ? aimesh->mTextureCoords[0][t] : aiVector3D(0,0,0);
            aiVector3D aitan = aimesh->HasTangentsAndBitangents() ? normalTr*aimesh->mTangents[t] :  aiVector3D(1,0,0);


            meshdata->vertices.push_back(VertexData(vec3(aipnt.x, aipnt.y, aipnt.z),
                                                    normalize(vec3(ainrm.x, ainrm.y, ainrm.z)),
                                                    vec2(aitex.x, aitex.y),
                                                    vec3(aitan.x, aitan.y, aitan.z))); }
        
        // Loop through all faces, recording indices
        for (unsigned int t=0;  t<aimesh->mNumFaces;  ++t) {
            aiFace* aiface = &aimesh->mFaces[t];
            for (int i=2;  i<aiface->mNumIndices;  i++) {
                meshdata->triangles.push_back(ivec3(aiface->mIndices[0],
                                                    aiface->mIndices[i-1],
                                                    aiface->mIndices[i])); } }
        meshdata->mat = materials[aimesh->mMaterialIndex];
        //printf("Mesh %d uses %d\n", m, aimesh->mMaterialIndex);
        scene->triangleMesh(meshdata); }


    // Recurse onto this node's children
    for (unsigned int i=0;  i<node->mNumChildren;  ++i)
        recurseModelNodes(scene, aiscene, node->mChildren[i], childTr, materials, level+1);
}
#endif
