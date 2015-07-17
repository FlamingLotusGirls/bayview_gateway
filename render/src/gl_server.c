/* Copyright 2013 Ka-Ping Yee

Licensed under the Apache License, Version 2.0 (the "License"); you may not
use this file except in compliance with the License.  You may obtain a copy
of the License at: http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied.  See the License for the
specific language governing permissions and limitations under the License. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <ctype.h>
#include <errno.h>
#include <jpeglib.h>
#ifdef __APPLE__
#include <OpenGL/CGLCurrent.h>
#include <OpenGL/CGLTypes.h>
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "cJSON.h"
#include "opc.h"

typedef struct {
  float ambient[4];
  float diffuse[4];
  float specular[4];
  float emission[4];
  float shininess;
  char name[32];
} Material;

// defining a cylindrical projection, for now
typedef struct {
  GLuint  id;              // opengl id
  char  name[32];        // name, for internal use/debugging
  long  textureWidth;    // width of original jpeg
  long  textureHeight;   // height of original jpeg
  double opacity;         // for overlaying with other textures
  double winding;         // How many times the texture should wrap around the axis
  int   height;          // height of texture (will scale to fit)
  float transformMat[9]; // cylinder defined by z-axis, so transform into world coordinates
  double theta;           // current rotation on axis
  double dThetadt;        // radians/second to spin the transformation
  double zOffset;          // world units/second to translate texture on local z-axis
  double dZdt;            // how fast to move the texture, in units/second  
  unsigned char* textureBytes;    // rgb triples. For now no alpha, no stride
} Texture;

typedef struct {
  long startIdx;
  Material *material;
  Texture  *texture;
  char name[16];
} PartDefinition;

typedef enum {
    MATERIALS_SECTION,
    TEXTURES_SECTION,
    PARTS_SECTION
} SectionType;


// forward references
unsigned char* readJpeg(const char *filename, long *width, long *height);
void matMultiplyf(float *mat1, float *mat2, float *out);
void matMultiply3vf(float *matrix, float *vertex, float *result);
int parseTexture(cJSON *textureSpec, Texture *texture);
int parseMaterial(cJSON *materialSpec, Material *material);
int parsePart(cJSON *partSpec, PartDefinition *part);
char* read_file(char* filename);
void init_gl_textures();


opc_source source = -1;
int verbose = 0;
int bUseJSONColor = 0; // Option set by 'jsoncolor' on command line... uses color in the JSON file rather than color from the client
int gUseNormals = 1;   // whether we use the normals in the model or not
int triangle_idx_start = 0;
int triangle_idx_end = 1000;

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef PI
#define PI 3.14159264
#endif

// Camera parameters
#define FOV_DEGREES 20
#define WINDOW_WIDTH 800
int orbiting = 0, dollying = 0;
double start_angle, start_elevation, start_distance;
int start_x, start_y;
double orbit_angle = 192.0;  // camera orbit angle, degrees
double camera_elevation = -15;  // camera elevation angle, degrees
double camera_distance = 38.0;  // distance from origin, metres
double camera_aspect = 1.0;  // will be updated to match window aspect ratio
double scale_factor = 1.0;

// Adding panning, at least in world coordinates. This is a bit of a hack
// from a UI perspective, but it works... Use keys x/X, y/Y, z/Z to move 
// the world along the axis.
double world_x = 0.0; 
double world_y = 0.0;
double world_z = 0.0;

// Extents of model, for better calculation of initial
// camera distance. These are absolute values - we're still looking 
// at the origin
/*
double model_max_x = 10.0;
double model_max_y = 10.0;
double model_max_z = 10.0;
*/

// note that textures need to be an power of two height and width


#define MAX_PARTS     1024
#define MAX_MATERIALS 256
#define MAX_TEXTURES  256
/*
#define MATERIAL_HEADER "MATERIALS"
#define TEXTURE_HEADER  "TEXTURES"
#define PART_HEADER     "PARTS"
*/

Material materials[MAX_MATERIALS]; 
Texture  textures[MAX_TEXTURES];
PartDefinition modelParts[MAX_PARTS];
int numParts = 0;
int numTextures = 0;
int numMaterials = 0;

void initIdentityMatf(float *mat)
{
    mat[0] = 1.0f;
    mat[1] = 0.0f;
    mat[2] = 0.0f;
    mat[3] = 0.0f;
    mat[4] = 1.0f;
    mat[5] = 0.0f;
    mat[6] = 0.0f;
    mat[7] = 0.0f;
    mat[8] = 1.0f;
}

void init_gl_textures()
{
    // Put the textures into OpenGL so they can be referenced later
    for (int i=0; i<numTextures; i++) {
        GLuint texId;
        glGenTextures(1, &textures[i].id); 
        glBindTexture(GL_TEXTURE_2D, textures[i].id);
        // simple texture load 
        glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, textures[i].textureWidth, textures[i].textureHeight, 0, GL_BGR, GL_UNSIGNED_BYTE, textures[i].textureBytes);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    }
}

int parseTexture(cJSON *textureSpec, Texture *texture) 
{
    if (!textureSpec || !texture) {
        return FALSE;
    }

    char *name = cJSON_GetObjectItem(textureSpec, "name")->valuestring;
    char *textureFileName = cJSON_GetObjectItem(textureSpec, "filename")->valuestring;
    if (!name || !textureFileName) {
        printf("Texture spec lacks either name or texture file, ignoring\n");
        goto ErrExit;
    }

    printf("Parsing texture - name: %s, filename: %s\n", name, textureFileName);
    if (strlen(name) > 32) {
        printf("Texture name %s is too long, truncating\n", name);
        name[31] = '\0';
    }
    strcpy(texture->name, name);
    

    texture->opacity  = cJSON_GetObjectItem(textureSpec, "opacity")->valuedouble; 
    texture->winding  = cJSON_GetObjectItem(textureSpec, "winding")->valuedouble;
    texture->height   = cJSON_GetObjectItem(textureSpec, "height")->valueint;
    texture->theta    = cJSON_GetObjectItem(textureSpec, "theta")->valuedouble;
    texture->dThetadt = cJSON_GetObjectItem(textureSpec, "dTheta")->valuedouble;
    texture->zOffset  = cJSON_GetObjectItem(textureSpec, "z")->valuedouble;
    texture->dZdt     = cJSON_GetObjectItem(textureSpec, "dZ")->valuedouble;
    
    //printf("opacity is %f, winding is %f, height is %d, theta is %f zOffset is %f\n", tmpTexture->opacity, tmpTexture->winding, tmpTexture->height, tmpTexture->theta, tmpTexture->zOffset); 

    double rotX = cJSON_GetObjectItem(textureSpec, "rotX")->valuedouble;
    double rotY = cJSON_GetObjectItem(textureSpec, "rotY")->valuedouble;
    double rotZ = cJSON_GetObjectItem(textureSpec, "rotZ")->valuedouble;

    float matRotX[9];
    float matRotY[9];
    float matRotZ[9];
    float matIdentity[9];
    
    initIdentityMatf(matRotX);
    initIdentityMatf(matRotY);
    initIdentityMatf(matRotZ);
    
    printf("rot X %f, rot Y %f, rot Z %f\n", rotX, rotY, rotZ);
    printf("matRotX is %f, %f, %f, %f,%f,%f,  %f,%f,%F\n", 
          matRotX[0],  matRotX[1],  matRotX[2],
          matRotX[3],  matRotX[4],  matRotX[5],
          matRotX[6],  matRotX[7],  matRotX[8]);
    
    matRotX[4] = cos(rotX);
    matRotX[5] = sin(rotX);
    matRotX[7] = -sin(rotX);
    matRotX[8] = cos(rotX);
    
    matRotY[0] = cos(rotY);
    matRotY[2] = sin(rotY);
    matRotY[6] = -sin(rotY);
    matRotY[8] = cos(rotY);
    
    matRotZ[0] = cos(rotZ);
    matRotZ[1] = sin(rotZ);
    matRotZ[3] = -sin(rotZ);
    matRotZ[4] = cos(rotZ);

    printf("matRotX is %f, %f, %f, %f,%f,%f,  %f,%f,%F\n", 
          matRotX[0],  matRotX[1],  matRotX[2],
          matRotX[3],  matRotX[4],  matRotX[5],
          matRotX[6],  matRotX[7],  matRotX[8]);
    
    // create rotation matrix out of rotX, rotY, and rot Z.
    float matOut[9], matOut2[9];
    matMultiplyf(matRotY, matRotX, matOut);
    matMultiplyf(matRotZ, matOut, matOut2);

    printf("matOut is %f, %f, %f, %f,%f,%f,  %f,%f,%F\n", 
          matOut[0],  matOut[1],  matOut[2],
          matOut[3],  matOut[4],  matOut[5],
          matOut[6],  matOut[7],  matOut[8]);


    printf("matOut2 is %f, %f, %f, %f,%f,%f,  %f,%f,%F\n", 
          matOut2[0],  matOut2[1],  matOut2[2],
          matOut2[3],  matOut2[4],  matOut2[5],
          matOut2[6],  matOut2[7],  matOut2[8]);
    
    memcpy(texture->transformMat, matOut2, sizeof(float) * 9);

    long width, height;
    unsigned char *textureBytes = readJpeg(textureFileName, &width, &height);
    if (!textureBytes) {
        printf("Texture file could not be read, ignoring texture %s\n", name);
        goto ErrExit;
    }
    printf("Setting up texture\n");
    texture->textureBytes = textureBytes;
    texture->textureWidth = width;
    texture->textureHeight = height;
    printf("about to call the glthingy\n");
   // glGenTextures(1, &tmpTexture->id); // XXX why is this having trouble? Probably before gl init?
    
    return TRUE;

ErrExit:
    return FALSE;
}

void setUpModelParts(char *filename)
{
   if (!filename) return;

   // clear buffers
   memset(materials,  0, MAX_MATERIALS * sizeof(Material));
   memset(textures,   0, MAX_TEXTURES  * sizeof(Texture));
   memset(modelParts, 0, MAX_PARTS     * sizeof(PartDefinition));
   numMaterials = 0;
   numTextures = 0;
   numParts = 0;

    char *bytes = read_file(filename);
    if (!bytes) {
        printf("Could not read %s\n", filename);
        return;
    }
    
    cJSON* root = cJSON_Parse(bytes);
    if (!root) {
        printf("Could not parse file %s\n", filename);
        goto ErrExit;
    }
    
    cJSON* materialsJSON = cJSON_GetObjectItem(root, "materials");
    cJSON* texturesJSON  = cJSON_GetObjectItem(root, "textures");
    cJSON* partsJSON     = cJSON_GetObjectItem(root, "parts");
    
    if (!partsJSON) {
        printf("No parts defined in model of file %s, skipping parts setup\n", filename);
        goto ErrExit;
    }
    // bring in materials
    int nMaterials = cJSON_GetArraySize(materialsJSON);
    if (nMaterials > MAX_MATERIALS) {
        printf("Too many materials, truncating\n");
        nMaterials = MAX_MATERIALS;
    }
    
    for (int i=0; i<nMaterials; i++) {
        if (parseMaterial(cJSON_GetArrayItem(materialsJSON, i), &materials[numMaterials])) {
            numMaterials++;
        }
    }
    
    // and textures
    int nTextures = cJSON_GetArraySize(texturesJSON);
    if (nTextures > MAX_TEXTURES) {
        printf("Too many textures, truncating\n");
        nTextures = MAX_TEXTURES;
    }
    
    for (int i=0; i<nTextures; i++) {
        if (parseTexture(cJSON_GetArrayItem(texturesJSON, i), &textures[numTextures])) {
            numTextures++;
        }
    }

    // and parts
    int nParts = cJSON_GetArraySize(partsJSON);
    if (nParts > MAX_PARTS) {
        printf("Too many parts, truncating\n");
        nParts = MAX_PARTS;
    }
    
    for (int i=0; i<nParts; i++) {
        if (parsePart(cJSON_GetArrayItem(partsJSON, i), &modelParts[numParts])) {
            numParts++;
        }
    }
    
ErrExit:
    if (bytes) {
        free(bytes);
    }
    if (root) {
        cJSON_Delete(root);
    }

}

// Attempt to load material. Return TRUE or FALSE depending on whether material is loaded
int parseMaterial(cJSON *materialSpec, Material *material)
{
    if (!materialSpec || !material) {
      return FALSE;
    }
    
    char *name = cJSON_GetObjectItem(materialSpec, "name")->valuestring;
    
    strncpy(material->name, name, 32); 
    material->name[31] = '\0';
    
    
    cJSON *ambient   = cJSON_GetObjectItem(materialSpec, "ambient");
    cJSON *diffuse   = cJSON_GetObjectItem(materialSpec, "diffuse");
    cJSON *emission  = cJSON_GetObjectItem(materialSpec, "emmissive");
    cJSON *specular  = cJSON_GetObjectItem(materialSpec, "specular");
    cJSON *shininess = cJSON_GetObjectItem(materialSpec, "shininess");
    
    if (ambient && cJSON_GetArraySize(ambient) >= 4) {
        material->ambient[0] = cJSON_GetArrayItem(ambient, 0)->valuedouble;
        material->ambient[1] = cJSON_GetArrayItem(ambient, 1)->valuedouble;
        material->ambient[2] = cJSON_GetArrayItem(ambient, 2)->valuedouble;
        material->ambient[3] = cJSON_GetArrayItem(ambient, 3)->valuedouble;
    } else {
        printf("Found no ambient lighting for material %s\n", name);
        material->ambient[0] = 0.0f;
        material->ambient[1] = 0.0f;
        material->ambient[2] = 0.0f;
        material->ambient[3] = 1.0f;
    }

    if (diffuse && cJSON_GetArraySize(diffuse) >= 4) {
        material->diffuse[0] = cJSON_GetArrayItem(diffuse, 0)->valuedouble;
        material->diffuse[1] = cJSON_GetArrayItem(diffuse, 1)->valuedouble;
        material->diffuse[2] = cJSON_GetArrayItem(diffuse, 2)->valuedouble;
        material->diffuse[3] = cJSON_GetArrayItem(diffuse, 3)->valuedouble;
    } else {
        printf("Found no diffuse lighting for material %s\n", name);
        material->diffuse[0] = 0.0f;
        material->diffuse[1] = 0.0f;
        material->diffuse[2] = 0.0f;
        material->diffuse[3] = 1.0f;
    }
    
    if (specular && cJSON_GetArraySize(specular) >= 4) {
        material->specular[0] = cJSON_GetArrayItem(specular, 0)->valuedouble;
        material->specular[1] = cJSON_GetArrayItem(specular, 1)->valuedouble;
        material->specular[2] = cJSON_GetArrayItem(specular, 2)->valuedouble;
        material->specular[3] = cJSON_GetArrayItem(specular, 3)->valuedouble;
    } else {
        printf("Found no specular lighting for material %s\n", name);
        material->specular[0] = 0.0f;
        material->specular[1] = 0.0f;
        material->specular[2] = 0.0f;
        material->specular[3] = 1.0f;
    }

    if (emission && cJSON_GetArraySize(emission) >= 4) {
        material->emission[0] = cJSON_GetArrayItem(emission, 0)->valuedouble;
        material->emission[1] = cJSON_GetArrayItem(emission, 1)->valuedouble;
        material->emission[2] = cJSON_GetArrayItem(emission, 2)->valuedouble;
        material->emission[3] = cJSON_GetArrayItem(emission, 3)->valuedouble;
    } else {
        printf("Found no emission lighting for material %s\n", name);
        material->emission[0] = 0.0f;
        material->emission[1] = 0.0f;
        material->emission[2] = 0.0f;
        material->emission[3] = 1.0f;
    }
    

    if (shininess) {
        material->shininess = shininess->valuedouble;
    }
    return TRUE;

ErrExit:
    return FALSE;
}

int parsePart(cJSON *partSpec, PartDefinition *part)
{
    if (!partSpec || !part) {
        return FALSE;
    }
    
    cJSON *nameJson = cJSON_GetObjectItem(partSpec, "name");
    cJSON *materialJson = cJSON_GetObjectItem(partSpec, "material");
    cJSON *textureJson  = cJSON_GetObjectItem(partSpec, "texture");
    cJSON *vertexStart = cJSON_GetObjectItem(partSpec, "startVertex");
    if (!vertexStart) {
        printf("Invalid texture\n");
        goto ErrExit;
    }
    part->startIdx = vertexStart->valueint;

    if (nameJson && nameJson->valuestring) {
        strncpy(part->name, nameJson->valuestring, 16);
        part->name[15] = '\0';
    }

    if (materialJson && materialJson->valuestring) {
        int foundMaterial = FALSE;
        char *name = materialJson->valuestring;
        for (int i=0; i<numMaterials; i++) {
           if (!strcasecmp(materials[i].name, name)) { 
             part->material = &materials[i];
             foundMaterial = TRUE;
             //printf("Found material for part %d\n", numParts);
             break;
           }
        }

        if (!foundMaterial) {
          printf("Did not find material for part %d, looking for %s\n", numParts,  name); 
          part->material = NULL;
        }
    }

    if (textureJson && textureJson->valuestring){
        int foundTexture = FALSE;
        char *name = textureJson->valuestring;
        for (int i=0; i<numTextures; i++) {
           if (!strcasecmp(textures[i].name, name)) {
             part->texture = &textures[i];
             foundTexture = TRUE;
             //printf("Found texture for part %d\n", numParts);
             break;
           }
        }

        if (!foundTexture) {
          printf("Did not find texture for part %d, looking for %s\n", numParts, name); 
          part->texture = NULL;
        }
    }
    
    /*printf("Part %s starts at vertex %ld, uses material %s and texture %s\n", 
          part->name ? part->name: "no name", 
          part->startIdx, 
          part->material ? part->material->name:"None", 
          part->texture ? part->texture->name: "None");*/
    
    
    return TRUE;
    
ErrExit:
    return FALSE;
}

// Shape parameters
#define SHAPE_THICKNESS 0.18  // thickness of points and lines, metres

// LED colours
#define MAX_PIXELS 30000
int num_pixels = 0;
pixel pixels[MAX_PIXELS];

// Floating-point colours
typedef struct {
  double r, g, b;
} colour;

colour tmp_colour;
#define set_colour(c) ((tmp_colour = c), glColor3dv(&(tmp_colour.r)))
#define set_rgb(r, g, b) (glColor3d(r, g, b))
colour xfer[256];

// Vector arithmetic
typedef struct {
  double x, y, z;
} vector;

vector tmp_vector;
#define put_vertex(v) ((tmp_vector = v), glVertex3dv(&(tmp_vector.x)))
#define put_pair(v, w) (put_vertex(v), put_vertex(w))

vector add(vector v, vector w) {
  vector result;
  result.x = v.x + w.x;
  result.y = v.y + w.y;
  result.z = v.z + w.z;
  return result;
}

vector subtract(vector v, vector w) {
  vector result;
  result.x = v.x - w.x;
  result.y = v.y - w.y;
  result.z = v.z - w.z;
  return result;
}

vector multiply(double f, vector v) {
  vector result;
  result.x = f*v.x;
  result.y = f*v.y;
  result.z = f*v.z;
  return result;
}

double length(vector v) {
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

double dot(vector v, vector w) {
  return v.x*w.x + v.y*w.y + v.z*w.z;
}

vector cross(vector v, vector w) {
  vector result;
  result.x = v.y*w.z - w.y*v.z;
  result.y = v.z*w.x - w.z*v.x;
  result.z = v.x*w.y - w.x*v.y;
  return result;
}

// row major format for the multiply....
void matMultiplyf(float *mat1, float *mat2, float *out)
{
    out[0] = mat1[0]*mat2[0] + mat1[1]*mat2[3] + mat1[2]*mat2[6];
    out[1] = mat1[0]*mat2[1] + mat1[1]*mat2[4] + mat1[2]*mat2[7];
    out[2] = mat1[0]*mat2[2] + mat1[1]*mat2[5] + mat1[2]*mat2[8];

    out[3] = mat1[3]*mat2[0] + mat1[4]*mat2[3] + mat1[5]*mat2[6];
    out[4] = mat1[3]*mat2[1] + mat1[4]*mat2[4] + mat1[5]*mat2[7];
    out[5] = mat1[3]*mat2[2] + mat1[4]*mat2[5] + mat1[5]*mat2[8];

    out[6] = mat1[6]*mat2[0] + mat1[7]*mat2[3] + mat1[8]*mat2[6];
    out[7] = mat1[6]*mat2[1] + mat1[7]*mat2[4] + mat1[8]*mat2[7];
    out[8] = mat1[6]*mat2[2] + mat1[7]*mat2[5] + mat1[8]*mat2[8];
}


// Shapes
typedef struct shape {
  void (*draw)(struct shape* this, GLUquadric* quad);
  int index;
  double red;
  double green;
  double blue;
  union {
    vector point;
    struct { vector start, end; } line;
  } g;
} shape;

#define MAX_SHAPES 30000
int num_shapes = 0;
shape shapes[MAX_SHAPES];

void draw_point(shape* this, GLUquadric* quad) {
  pixel p = pixels[this->index];
  if (bUseJSONColor) {
    glColor3d(this->red, this->green, this->blue);
  } else {
    glColor3d(xfer[p.r].r, xfer[p.g].g, xfer[p.b].b);
  }
  glPushMatrix();
  glTranslatef(this->g.point.x, this->g.point.y, this->g.point.z);
  gluSphere(quad, SHAPE_THICKNESS/2, 6, 3);
  glPopMatrix();
}

void draw_line(shape* this, GLUquadric* quad) {
  pixel p = pixels[this->index];
  vector start = this->g.line.start;
  vector delta = subtract(this->g.line.end, this->g.line.start);
  vector z = {0, 0, 1};
  vector hinge = cross(z, delta);
  double len = length(delta);
  double angle = 180./M_PI * acos(dot(z, delta) / len);
  glColor3d(xfer[p.r].r, xfer[p.g].g, xfer[p.b].b);
  glPushMatrix();
  glTranslated(start.x, start.y, start.z);
  glRotated(angle, hinge.x, hinge.y, hinge.z);
  gluSphere(quad, SHAPE_THICKNESS/2, 6, 3);
  gluCylinder(quad, SHAPE_THICKNESS/2, SHAPE_THICKNESS/2, len, 6, 1);
  glTranslated(0, 0, len);
  gluSphere(quad, SHAPE_THICKNESS/2, 6, 3);
  glPopMatrix();
}

void draw_axes() {
  vector o = {0, 0, 0};
  vector x = {1, 0, 0};
  vector y = {0, 1, 0};
  vector z = {0, 0, 1};
  vector xx = {10, 0, 0};
  vector yy = {0, 10, 0};
  vector zz = {0, 0, 10};
  glLineWidth(2);
  glBegin(GL_LINES);
  set_rgb(0.3, 0.3, 0.3);
  put_pair(o, x);
  put_pair(o, y);
  put_pair(o, z);
  set_rgb(0.3, 0, 0);
  put_pair(x, xx);
  set_rgb(0, 0.3, 0);
  put_pair(y, yy);
  set_rgb(0, 0, 0.3);
  put_pair(z, zz);
  glEnd();
}

// Mesh loaded from STL file, if any
float * mesh = NULL; // array of triangle vertices
int num_triangles = 0;


void load_stl_mesh(char* filename) {

  // check for valid file extension
  char* ext = strrchr(filename,'.');
  if (!ext || (strcmp(ext+1, "stl") && strcmp(ext+1, "STL")) || strlen(ext+1) > 3) {
    printf("Invalid filename: %s. Mesh file must be an STL file.\n", filename);
    return;
  }

  char* data = read_file(filename);
  if (!data) {
    printf("Invalid filename: %s. File does not exist or can't be read.\n", filename);
    return;
  }

  int float_size = sizeof(float);
  if (float_size == 4)
  {
    char* file_loc = data;
    int triangleSize = 9;

    // STL has 80-byte header followed by uint32 telling you the # of triangles in the file
    file_loc += 80;
    num_triangles = *(u32*)(file_loc);
    printf("Number of triangles is %d\n", num_triangles);

    file_loc += 4; //skip past triangle count
    if (!gUseNormals) {
      file_loc += 3*float_size; //skip normal of first triangle
    } else {
      triangleSize += 3;
    }

    mesh = (float*)malloc(num_triangles*triangleSize*float_size);
    if (!mesh) {
      printf("Could not allocate memory for STL mesh data.\n");
      num_triangles = 0;
    }
    else {
      float* loc = mesh;
      int i;

      // each STL triangle is 3 normal floats, then 9 vertex floats, then a uint16.
//      float *initialVtxPtr = ((float*)file_loc) + 3*gUseNormals;
      float max_X, max_Y, max_Z, min_X, min_Y, min_Z;
      float offsetX = 3450;
      float offsetY = 311;
      float offsetZ = 3468;
      float scale = 1.0f;
      for (i = 0; i < num_triangles; i++) {
        memcpy(loc, file_loc, triangleSize*float_size);
        float *vertex_ptr = loc + 3*gUseNormals;
        // and now I'm just going to do a little translation...
        *(vertex_ptr)   -= offsetX;
        *(vertex_ptr+1) -= offsetY;
        *(vertex_ptr+2) -= offsetZ;
        *(vertex_ptr+3) -= offsetX;
        *(vertex_ptr+4) -= offsetY;
        *(vertex_ptr+5) -= offsetZ;
        *(vertex_ptr+6) -= offsetX;
        *(vertex_ptr+7) -= offsetY;
        *(vertex_ptr+8) -= offsetZ;
        if (i==0) {
            max_X = min_X = *(vertex_ptr);
            max_Y = min_Y = *(vertex_ptr + 1);
            max_Z = min_Z = *(vertex_ptr + 2);
        }

        *(vertex_ptr)   = *(vertex_ptr)/scale;
        *(vertex_ptr+1) = *(vertex_ptr+1)/scale;
        *(vertex_ptr+2) = *(vertex_ptr+2)/scale;
        *(vertex_ptr+3) = *(vertex_ptr+3)/scale;
        *(vertex_ptr+4) = *(vertex_ptr+4)/scale;
        *(vertex_ptr+5) = *(vertex_ptr+5)/scale;
        *(vertex_ptr+6) = *(vertex_ptr+6)/scale;
        *(vertex_ptr+7) = *(vertex_ptr+7)/scale;
        *(vertex_ptr+8) = *(vertex_ptr+8)/scale;
        
        if (*(vertex_ptr) > max_X) max_X = *(vertex_ptr);
        if (*(vertex_ptr) < min_X) { min_X = *(vertex_ptr); /*printf("Setting min_X to %f\n", min_X);*/}
        if (*(vertex_ptr+1) > max_Y) max_Y = *(vertex_ptr+1); // Y, +1
        if (*(vertex_ptr+1) < min_Y) min_Y = *(vertex_ptr+1); // Y, +1
        if (*(vertex_ptr+2) > max_Z) max_Z = *(vertex_ptr+2);
        if (*(vertex_ptr+2) < min_Z) min_Z = *(vertex_ptr+2);
        if (*(vertex_ptr+3) > max_X) max_X = *(vertex_ptr+3);
        if (*(vertex_ptr+3) < min_X) min_X = *(vertex_ptr+3);
        if (*(vertex_ptr+4) > max_Y) max_Y = *(vertex_ptr+4); //Y, +4
        if (*(vertex_ptr+4) < min_Y) min_Y = *(vertex_ptr+4); //Y, +4
        if (*(vertex_ptr+5) > max_Z) max_Z = *(vertex_ptr+5);
        if (*(vertex_ptr+5) < min_Z) min_Z = *(vertex_ptr+5);
        if (*(vertex_ptr+6) > max_X) max_X = *(vertex_ptr+6);
        if (*(vertex_ptr+6) < min_X) min_X = *(vertex_ptr+6);
        if (*(vertex_ptr+7) > max_Y) max_Y = *(vertex_ptr+7); //Y, +7
        if (*(vertex_ptr+7) < min_Y) min_Y = *(vertex_ptr+7); //Y, +7
        if (*(vertex_ptr+8) > max_Z) max_Z = *(vertex_ptr+8);
        if (*(vertex_ptr+8) < min_Z) min_Z = *(vertex_ptr+8);
        if (i%101 == 0) {
           //printf("Uint16 on end is %f\n", *(vertex_ptr+10));
        }
        loc += triangleSize;
        file_loc += 12*float_size + 2;
      }

      // set initial camera distance so that the model fits in the view
      
   
      double model_max_x = abs(min_X) > abs(max_X) ? abs(min_X) : abs(max_X);
      double model_max_y = abs(min_Y) > abs(max_Y) ? abs(min_Y) : abs(max_Y);
      double model_max_z = abs(min_Z) > abs(max_Z) ? abs(min_Z) : abs(max_Z);
     
      double max_extents = model_max_x > model_max_y ? model_max_x : model_max_y;
      max_extents = model_max_z > max_extents ? model_max_z : max_extents;

//      max_extents = max_extents/0.8f;
      max_extents = max_extents/1.0f;
      printf("max_extents is %f\n", max_extents);

      camera_distance = 2*max_extents/tan(FOV_DEGREES/2);
      scale_factor = 10.0/max_extents;
      //printf("scale factor is %f\n", scale_factor);
      //scale_factor = 1.0f;
 
      printf("X bounds are %f, %f\n", min_X, max_X);
      printf("Y bounds are %f, %f\n", min_Y, max_Y);
      printf("Z bounds are %f, %f\n", min_Z, max_Z);
      printf("*** X range is %f ***\n", max_X-min_X);
      printf("Camera distance is %f\n", camera_distance);

    }
  }
  free(data);
}

void render_triangles(float* triangles, int nTriangles)
{
  Texture *currentTexture = NULL;
 // glEnable(GL_TEXTURE_2D);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
 
  if (triangles) {
    glBegin(GL_TRIANGLES);
 //   printf("render\n");
 
    // render everything in grey (for now)
    //glColor3d(1.0, 0.1, 0.1);
    float tempMaterial[] = {1.0f, 0.2f, 0.2f, 1.0f, 0.2f,0.2f,0.2f,1.0f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, &tempMaterial[4]);
    
    PartDefinition *nextPart;

    int partIdx = 0;
    if (numParts > partIdx) {
      nextPart = &modelParts[partIdx];
    } else {
      nextPart = NULL;
    }

    int i;
    int perTriangle  = 9;  // number of coordinates per triangle, 9 if no normals, 12 if normals
    int vertexOffset = 0;  // offset from start of triangle to vertex data. 0 if no normals, 3 if normals
    if (gUseNormals) {
      perTriangle  += 3;
      vertexOffset += 3;
    }
    
    int start=0;
    
    printf("nTriangles is %d\n", nTriangles);
    for(i = 0; i < nTriangles; i++){
      // Setup environment on part boundaries
      if (nextPart != NULL && nextPart->startIdx == i) {
        glEnd();
        
        // Set textures
        if (nextPart->texture != NULL) {
            currentTexture = nextPart->texture;
//            printf("Using texture %s on part %s\n", currentTexture->name, nextPart->name);
            glEnable(GL_TEXTURE_2D);
            //glClear(GL_COLOR_BUFFER_BIT);
            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
            glBindTexture(GL_TEXTURE_2D, currentTexture->id);
        } else {
            currentTexture = NULL;
            glDisable(GL_TEXTURE_2D);
        }
        
        // Set materials
        if (nextPart->material == NULL) {
          glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,&tempMaterial[4]);
        } else {
          glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, nextPart->material->ambient);
          glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, nextPart->material->diffuse);
          glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, nextPart->material->specular);
          glMaterialf(GL_FRONT, GL_SHININESS, nextPart->material->shininess * 128.0);
        }          

        // Prep for next part
        partIdx++;
        if (partIdx < numParts) {
          nextPart = &modelParts[partIdx];
        } else {
          nextPart = NULL;
        }
        
        // Re-init drawing with new environment
        glBegin(GL_TRIANGLES);
      } 
      
      // set the normal if we're using normals for the lighting mode
      if (gUseNormals) {
        glNormal3fv(triangles+(12*i));
      }

      if (currentTexture) {
        // coordinates of a point are...  Where are my modelworld matrices?
        // map texture projection into model space. Then v is model z / texture height + animation z factor
        float result[3];

        matMultiply3vf(currentTexture->transformMat, triangles+(perTriangle*i) + vertexOffset, result); // should do 4x4 so I can do translation XXX
        if (i==19500) {
        /*    printf("Transform mat is %f,%f,%f,   %f, %f, %f,  %f, %f,%f\n",
                   currentTexture->transformMat[0], currentTexture->transformMat[1], currentTexture->transformMat[2],
                   currentTexture->transformMat[3], currentTexture->transformMat[4], currentTexture->transformMat[5],
                   currentTexture->transformMat[6], currentTexture->transformMat[7], currentTexture->transformMat[8]); */
            printf("vertex is at %f %f %f, result is %f %f %f \n", 
                   *(triangles+(perTriangle*i) + vertexOffset),
                   *(triangles+(perTriangle*i) + vertexOffset + 1),
                   *(triangles+(perTriangle*i) + vertexOffset + 2),
                   result[0],
                   result[1],
                   result[2]);
        }
        // And u is atan2f(y, x) / 2PI * f, + animation spin factor
        float s = fmod(1,(atan2f(result[0], result[2])/(2*PI))*currentTexture->winding);
        float t = result[1]/currentTexture->height;
        if (i==19500) {
            printf("s=%f, t=%f\n", s, t);
        }
        glTexCoord2f(s, t); // no animation just yet
        //glTexCoord2f(0.01, 0.01); // no animation just yet
/*        if (start++ <10) {
            printf("tex coord is %f, %f\n", atan2f(result[1],result[0])/(2*PI*currentTexture->winding),result[2]); 
        }*/
        
      }
      glVertex3fv(triangles+(perTriangle*i) + vertexOffset);
      if (currentTexture) {
        float result[3];
        matMultiply3vf(currentTexture->transformMat, triangles+(perTriangle*i) + vertexOffset + 3, result);
//        glTexCoord2f(atan2f(result[1], result[0])/(2*PI*currentTexture->winding), result[2]); // no animation just yet
        float s = fmod(1,(atan2f(result[0], result[2])/(2*PI))*currentTexture->winding);
        float t = result[1]/currentTexture->height;
        glTexCoord2f(s, t); // no animation just yet
//       glTexCoord2f(0.01, 0.02); // no animation just yet
/*        if (start++ <10) {
            printf("tex coord is %f, %f\n", atan2f(result[1],result[0])/(2*PI*currentTexture->winding),result[2]); 
        }
*/
      }
      glVertex3fv(triangles+(perTriangle*i) + vertexOffset + 3);
      if (currentTexture) {
        float result[3];
        matMultiply3vf(currentTexture->transformMat, triangles+(perTriangle*i) + vertexOffset + 6, result); 
//        glTexCoord2f(atan2f(result[1], result[0])/(2*PI*currentTexture->winding), result[2]); // no animation just yet
        float s = fmod(1,(atan2f(result[0], result[2])/(2*PI))*currentTexture->winding);
        float t = result[1]/currentTexture->height;
//        glTexCoord2f(0.02, 0.02); // no animation just yet
/*        if (start++ <10) {
            printf("tex coord is %f, %f\n", atan2f(result[1],result[0])/(2*PI*currentTexture->winding),result[2]); 
        }
*/
      }
      glVertex3fv(triangles+(perTriangle*i) + vertexOffset + 6);
    }

    glEnd();
    glFlush();
  }
}

// multiply 3x3 matrix by a vector
void matMultiply3vf(float *matrix, float *vertex, float *result) 
{
    result[0] = matrix[0] * vertex[0] + matrix[1] * vertex[1] + matrix[2] * vertex[2];
    result[1] = matrix[3] * vertex[0] + matrix[4] * vertex[1] + matrix[5] * vertex[2];
    result[2] = matrix[6] * vertex[0] + matrix[7] * vertex[1] + matrix[5] * vertex[8];
}

// read jpg file, return allocated buffer of bytes
unsigned char* readJpeg(const char *filename, long *outwidth, long *outheight) {

  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  
  printf("Attempting to read jpeg file %s\n", filename);

  FILE * infile;  
//  JSAMPARRAY pJpegBuffer;
//  int stride;
  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return NULL;
  }
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, infile);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);
  long width  = cinfo.output_width;
  long height = cinfo.output_height;
  printf("jpeg structs inited\n");
  if (outwidth) {
    *outwidth = width;
  }
  if (outheight) {
    *outheight = height;
  }
  printf("malloc jpeg buf. width is %ld height is %ld\n", width,height);
  
  // would like to scale output to power of 2 for texture mapping purposes...

  printf("components is is %d\n", cinfo.output_components);
  // what's up with these lines... ?
  int  row_stride = width * cinfo.output_components ;
  unsigned char **pJpegBuffer = (*cinfo.mem->alloc_sarray) ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
  unsigned char *buffer = malloc(height*row_stride); // should be * components, but I don't want a stride. Deal.
  if (!buffer) {
    fprintf(stderr, "Cannot allocate memory for jpeg");
    return NULL;
  }
  
  // Note that we could invoke scaling here, but we're not going to right now.
  while (cinfo.output_scanline < cinfo.output_height) {
    int nLines = cinfo.output_scanline;
    (void) jpeg_read_scanlines(&cinfo, pJpegBuffer, 1);
    memcpy(buffer+(nLines*row_stride), pJpegBuffer[0], row_stride);
  }

  // exit now
  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(infile);
  
  return buffer;
}


void display() {
  int i;
  shape* sh;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  draw_axes();
  GLUquadric* quad = gluNewQuadric();
  for (i = 0, sh = shapes; i < num_shapes; i++, sh++) {
    sh->draw(sh, quad);
  }
  render_triangles(mesh, num_triangles);
  gluDeleteQuadric(quad);
  glutSwapBuffers();
}

void update_camera() {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(FOV_DEGREES, camera_aspect, 0.1, 1e3); // fov, aspect, zrange
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glScaled(scale_factor, scale_factor, scale_factor);
  double camera_y = -cos(camera_elevation*M_PI/180)*camera_distance;
  double camera_z = sin(camera_elevation*M_PI/180)*camera_distance;
  gluLookAt(camera_y, camera_y, camera_z, /* target */ 0, -5000, 0, /* up */ 0, -1, 0);
//  gluLookAt(0, camera_y, camera_z, /* target */ 0, 0, 0, /* up */ 0, 0, 1);
  glRotatef(orbit_angle, 0, 0, 1);
  glTranslatef(world_x, world_y, world_z);
  display();
}

// pull back - move along y/z axis
// pan left, pan right, pan up, pan down
// rotate around point at which you're looking - right, left, up, down

void reshape(int width, int height) {
  glViewport(0, 0, width, height);
  camera_aspect = ((double) width)/((double) height);
  update_camera();
}

void mouse(int button, int state, int x, int y) {
  if (state == GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
    dollying = 1;
    start_distance = camera_distance;
    start_x = x;
    start_y = y;
  } else if (state == GLUT_DOWN) {
    orbiting = 1;
    start_angle = orbit_angle;
    start_elevation = camera_elevation;
    start_x = x;
    start_y = y;
  } else {
    orbiting = 0;
    dollying = 0;
  }
}

void motion(int x, int y) {
  if (orbiting) {
    orbit_angle = start_angle + (x - start_x)*1.0;
    double elevation = start_elevation + (y - start_y)*1.0;
    camera_elevation = elevation < -89 ? -89 : elevation > 89 ? 89 : elevation;
    update_camera();
  }
  if (dollying) {
    double distance = start_distance + (y - start_y)*0.1;
    camera_distance = distance < 1.0 ? 1.0 : distance;
    update_camera();
  }
}

void keyboard(unsigned char key, int x, int y) {
  printf("have key %d\n", key);
  printf("X key is %d\n", 'X');
  if (key == '\x1b' || key == 'q') exit(0);
  if (key == 'x') {
    world_x++;
  } else if (key =='X') {
    world_x--;  
  } else if (key == 'y') {
    world_y++;
  } else if (key == 'Y') {
    world_y--;
  } else if (key == 'z') {
    world_z++;
  } else if (key == 'Z') {
    world_z--;
  } else if (key == 'i') {
    triangle_idx_start += 5;
    printf("update triangle_idx_start\n");
  } else if (key == 'I') {
    triangle_idx_start += 50;
  } else if (key == 'm') {
    triangle_idx_start -= 5;
  } else if (key == 'M') {
    triangle_idx_start -= 50;
  } else if (key == 'j') {
    triangle_idx_end += 5;
  } else if (key == 'J') {
    triangle_idx_end += 50;
  } else if (key == 'l') {
    triangle_idx_end -= 5;
  } else if (key == 'L') {
    triangle_idx_end -= 50;
  } else if (key == 'P') {
    printf("triangle start %d\n", triangle_idx_start);
    printf("triangle end %d\n", triangle_idx_end);
  }
  if (triangle_idx_start < 0) {
    triangle_idx_start = 0;
  } 
  if (triangle_idx_start > num_triangles) {
    triangle_idx_start = num_triangles;
  }
  if (triangle_idx_end < 0) {
    triangle_idx_end = 0;
  }
  if (triangle_idx_end > num_triangles) {
    triangle_idx_end = num_triangles;
  }
  if (triangle_idx_end < triangle_idx_start) {
    triangle_idx_end = triangle_idx_start;
  }
  
  update_camera();
}

void handler(u8 channel, u16 count, pixel* p) {
  int i = 0;

  if (verbose) {
    char* sep = " =";
    printf("-> channel %d: %d pixel%s", channel, count, count == 1 ? "" : "s");
    for (i = 0; i < count; i++) {
      if (i >= 4) {
        printf(", ...");
        break;
      }
      printf("%s %02x %02x %02x", sep, p[i].r, p[i].g, p[i].b);
      sep = ",";
    }
    printf("\n");
  }

  for (i = 0; i < count; i++) {
    pixels[i] = p[i];
  }
}

void idle() {
  /*
   * Receive all pending frames. We'll often draw slower than an OPC source
   * is producing pixels; to avoid runaway lag due to data buffered in the socket,
   * we want to skip frames.
   *
   * A short timeout (20 ms) on the first receive keeps us responsive to mouse events.
   * A zero timeout on subsequent receives lets us drain any queued frames without
   * waiting for them.
   */

  if (opc_receive(source, handler, 20) > 0) {

    // Drain queue
    while (opc_receive(source, handler, 0) > 0);

    // Show the last received frame
    display();
  }
}

char* read_file(char* filename) {
  FILE* fp;
  struct stat st;
  char* buffer;

  if (stat(filename, &st) != 0) {
    return strdup("");
  }
  buffer = malloc(st.st_size + 1);
  fp = fopen(filename, "r");
  fread(buffer, st.st_size, 1, fp);
  fclose(fp);
  buffer[st.st_size] = 0;
  return buffer;
}


// Read in the JSON file and ... do something to it.
void init(char* filename) {
  char* buffer;
  cJSON* json;
  cJSON* item;
  cJSON* index;
  cJSON* point;
  cJSON* color;
  cJSON* x;
  cJSON* line;
  cJSON* start;
  cJSON* x2;
  int i = 0;
  int isValidShape = 0;
  
  buffer = read_file(filename);
  json = cJSON_Parse(buffer);
  free(buffer);

  num_shapes = 0;
  for (item = json->child, i = 0; item; item = item->next, i++) {
    // Reset the index information if necessary
    index = cJSON_GetObjectItem(item, "index");
    if (index) {
      printf("Resetting index to %d\n", index->valueint);
      i = index->valueint;
    }
    // Extract point information, if it exists...
    point = cJSON_GetObjectItem(item, "point");
    x = point ? point->child : NULL;
    if (x && x->next && x->next->next) {
      shapes[num_shapes].draw = draw_point;
      shapes[num_shapes].index = i;
      shapes[num_shapes].g.point.x = x->valuedouble;
      shapes[num_shapes].g.point.y = x->next->valuedouble;
      shapes[num_shapes].g.point.z = x->next->next->valuedouble;
      isValidShape = 1;
    }

    // Extract line information, if we haven't already figured out what kind of shape
    // this is (a shape can be a point or a line, but not both)
    if (!isValidShape) {
        line = cJSON_GetObjectItem(item, "line");
        start = line ? line->child : NULL;
        x = start ? start->child : NULL;
        x2 = start && start->next ? start->next->child : NULL;
        if (x && x->next && x->next->next && x2 && x2->next && x2->next->next) {
          shapes[num_shapes].draw = draw_line;
          shapes[num_shapes].index = i;
          shapes[num_shapes].g.line.start.x = x->valuedouble;
          shapes[num_shapes].g.line.start.y = x->next->valuedouble;
          shapes[num_shapes].g.line.start.z = x->next->next->valuedouble;
          shapes[num_shapes].g.line.end.x = x2->valuedouble;
          shapes[num_shapes].g.line.end.y = x2->next->valuedouble;
          shapes[num_shapes].g.line.end.z = x2->next->next->valuedouble;
        }
    }

    // If we have a valid shape, extract color information and update index
    // into shape table.
    if (isValidShape) {
        color = cJSON_GetObjectItem(item, "color");
        x = color ? color->child : NULL;
        if (x && x->next &&  x->next->next) {
            shapes[num_shapes].red = x->valuedouble;
            shapes[num_shapes].green = x->next->valuedouble;
            shapes[num_shapes].blue = x->next->next->valuedouble;
        } else {
            shapes[num_shapes].red = 1.0;
            shapes[num_shapes].green = 1.0;
            shapes[num_shapes].blue = 1.0;
        }
        num_shapes++;
    }
  }
  
  // nb - the following clause does nothing, and we're relying
  // on the compiler to initialize the pixels array.
  // This is probably a bad idea....
  num_pixels = i;
  for (i = 0; i < num_pixels; i++) {
    pixels[i].r = pixels[i].g = pixels[i].b = 1;
  }
  
  // Initialize the channel transformation vectors. The pixel information
  // comes to us as a value between 0 and 255, and we need to transform 
  // that to something between 0 and 1.0f
  for (i = 0; i < 256; i++) {
    xfer[i].r = xfer[i].g = xfer[i].b = 0.1 + i*0.9/256;
  }
}

int main(int argc, char** argv) {
  u16 port;

  glutInitWindowSize(WINDOW_WIDTH, WINDOW_WIDTH*0.75);
  glutInit(&argc, argv);
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <gl-options> <filename.json> [<port>] [meshfile.stl] [modelfile]\n", argv[0]);
    exit(1);
  }
  init(argv[1]);
  port = argc > 2 ? strtol(argv[2], NULL, 10) : 0;
  port = port ? port : OPC_DEFAULT_PORT;
  source = opc_new_source(port);

  if (argc > 3) {
    load_stl_mesh(argv[3]);
  }
  
  // currently undocumented option to take shape color values from the initial json file.
  // Handy for some types of debugging.
 /* if (argc > 4){
    if (!strcasecmp("jsoncolor", argv[4])) {
        bUseJSONColor = 1;
    }
  }
  // XXX - this should be a different type of option!
 */ 

  if (argc > 4) {
    setUpModelParts(argv[4]);
  }

  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutCreateWindow("OPC");
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutIgnoreKeyRepeat(1);
  glutKeyboardFunc(keyboard);
  glutIdleFunc(idle);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  float sun_pos[4] = {0.0f, 1.0f, 1.0f, 0.0f};
  glLightfv(GL_LIGHT0, GL_POSITION, sun_pos);
  glEnable(GL_NORMALIZE);
  init_gl_textures();
#ifdef __APPLE__
  /* Make glutSwapBuffers wait for vertical refresh to avoid frame tearing. */
  int swap_interval = 1;
  CGLContextObj context = CGLGetCurrentContext();
  CGLSetParameter(context, kCGLCPSwapInterval, &swap_interval);
#endif

  glutMainLoop();
  return 0;
}


