#include "Lattice.h"
//#include "stdlib.h"
//#include "fcntl.h"
//#include "sys\stat.h"
//#include "io.h"





void getVoxels_d(double *Vm,  int sizeX, int sizeY, int sizeZ, double *** voxels)
{
//	printf("forming voxels\t");
    for(int x = 0; x < sizeX; x++)
    {
        for(int y = 0; y < sizeY; y++)
        {
            for (int z=0; z<sizeZ; z++)
            {
                    voxels[x][y][z] = Vm[z*sizeX*sizeY+y*sizeX+x];
            }
        }
    }
//	printf("Voxels formed\t");
}

vertex interpolate(double isolevel, vertex p1, vertex p2, int valp1, int valp2)
{
    if(fabs(isolevel - valp1) < 0.00001)
        return p1;
    if(fabs(isolevel - valp2) < 0.00001)
        return p2;
    if(fabs(valp1 - valp2) < 0.00001)
        return p1;

    vertex p;
    double diff = (double)(isolevel - valp1) / (valp2 - valp1);
    p.x = p1.x + diff * (p2.x - p1.x);
    p.y = p1.y + diff * (p2.y - p1.y);
    p.z = p1.z + diff * (p2.z - p1.z);

    p.normal_x = p1.normal_x + diff * (p2.normal_x - p1.normal_x);
    p.normal_y = p1.normal_y + diff * (p2.normal_y - p1.normal_y);
    p.normal_z = p1.normal_z + diff * (p2.normal_z - p1.normal_z);

    return p;
}

void processCube(cube cube, double isolevel, vector<vertex> *vertexList)
{
    int cubeindex = 0;
    if(cube.val[0] > isolevel) cubeindex |= 1;
    if(cube.val[1] > isolevel) cubeindex |= 2;
    if(cube.val[2] > isolevel) cubeindex |= 4;
    if(cube.val[3] > isolevel) cubeindex |= 8;
    if(cube.val[4] > isolevel) cubeindex |= 16;
    if(cube.val[5] > isolevel) cubeindex |= 32;
    if(cube.val[6] > isolevel) cubeindex |= 64;
    if(cube.val[7] > isolevel) cubeindex |= 128;

    // Cube is entirely in/out of the surface
    if(edgeTable[cubeindex] == 0 || edgeTable[cubeindex] == 255)
        return;

    vertex vertlist[12];
    // Find the vertices where the surface intersects the cube
    if(edgeTable[cubeindex] & 1)
        vertlist[0] = interpolate(isolevel,cube.p[0],cube.p[1],cube.val[0],cube.val[1]);
    if(edgeTable[cubeindex] & 2)
        vertlist[1] = interpolate(isolevel,cube.p[1],cube.p[2],cube.val[1],cube.val[2]);
    if(edgeTable[cubeindex] & 4)
        vertlist[2] = interpolate(isolevel,cube.p[2],cube.p[3],cube.val[2],cube.val[3]);
    if(edgeTable[cubeindex] & 8)
        vertlist[3] = interpolate(isolevel,cube.p[3],cube.p[0],cube.val[3],cube.val[0]);
    if(edgeTable[cubeindex] & 16)
        vertlist[4] = interpolate(isolevel,cube.p[4],cube.p[5],cube.val[4],cube.val[5]);
    if(edgeTable[cubeindex] & 32)
        vertlist[5] = interpolate(isolevel,cube.p[5],cube.p[6],cube.val[5],cube.val[6]);
    if(edgeTable[cubeindex] & 64)
        vertlist[6] = interpolate(isolevel,cube.p[6],cube.p[7],cube.val[6],cube.val[7]);
    if(edgeTable[cubeindex] & 128)
        vertlist[7] = interpolate(isolevel,cube.p[7],cube.p[4],cube.val[7],cube.val[4]);
    if(edgeTable[cubeindex] & 256)
        vertlist[8] = interpolate(isolevel,cube.p[0],cube.p[4],cube.val[0],cube.val[4]);
    if(edgeTable[cubeindex] & 512)
        vertlist[9] = interpolate(isolevel,cube.p[1],cube.p[5],cube.val[1],cube.val[5]);
    if(edgeTable[cubeindex] & 1024)
        vertlist[10] = interpolate(isolevel,cube.p[2],cube.p[6],cube.val[2],cube.val[6]);
    if(edgeTable[cubeindex] & 2048)
        vertlist[11] = interpolate(isolevel,cube.p[3],cube.p[7],cube.val[3],cube.val[7]);

    for(int i = 0; triTable[cubeindex][i] != -1; i++) {
        (*vertexList).push_back(vertlist[triTable[cubeindex][i]]);
    }
}



vector<vertex> runMarchingCubes_d(double ***voxels, double *Der, double isovalue)
{
    vector<vertex> vertexList;
    // Run the processCube function on every cube in the grid
            int sign;
    for(int x = 1; x < N-2; x++)
    {

        for(int y = 1; y < N-2; y++)
        {

            for(int z = 1; z < H-2; z++)
            {

                sign = -1;
                cube c =
                {
                    {
                    {
                        x,y,z,
                        (double)(-voxels[x+1][y][z]+voxels[x-1][y][z]) ,
                        (double)(-voxels[x][y+1][z]+voxels[x][y-1][z]) ,
                        (double)(-voxels[x][y][z+1]+voxels[x][y][z-1])
                    },
                    {
                        x+1,y,z,
                        (double)(-voxels[x+2*1][y][z]+voxels[x][y][z]) ,
                        (double)(-voxels[x+1][y+1][z]+voxels[x+1][y-1][z]),
                        (double)(-voxels[x+1][y][z+1]+voxels[x+1][y][z-1])
                    },
                    {
                        x+1,y,z+1,
                        (double)(-voxels[x+1][y][z+1]+voxels[x][y][z+1])  ,
                        (double)(-voxels[x+1][y+1][z+1]+voxels[x+1][y-1][z+1])  ,
                        (double)(-voxels[x+1][y][z+1]+voxels[x+1][y][z])
                    },
                    {
                        x,y,z+1,
                        (double)(-voxels[x+1][y][z+1]+voxels[x-1][y][z+1])  ,
                        (double)(-voxels[x][y+1][z+1]+voxels[x][y-1][z+1])  ,
                        (double)(-voxels[x][y][z+1]+voxels[x][y][z])
                    },
                    {
                        x,y+1,z,
                        (double)(-voxels[x+1][y+1][z]+voxels[x-1][y+1][z])  ,
                        (double)(-voxels[x][y+1][z]+voxels[x][y][z])  ,
                        (double)(-voxels[x][y+1][z+1]+voxels[x][y+1][z-1])
                    },
                    {
                        x+1,y+1,z,
                        (double)(-voxels[x+1][y+1][z]+voxels[x+1][y+1][z])  ,
                        (double)(-voxels[x+1][y+1][z]+voxels[x+1][y][z])  ,
                        (double)(-voxels[x+1][y+1][z+1]+voxels[x+1][y+1][z-1])
                    },
                    {x+1,y+1,z+1,
                        (double)(-voxels[x+1][y+1][z+1]+voxels[x][y+1][z+1])  ,
                        (double)(-voxels[x+1][y+1][z+1]+voxels[x+1][y][z+1])  ,
                        (double)(-voxels[x+1][y+1][z+1]+voxels[x+1][y+1][z])
                    },
                    {
                        x,y+1,z+1,
                        (double)(-voxels[x+1][y+1][z+1]+voxels[x-1][y+1][z+1])  ,
                        (double)(-voxels[x][y+1][z+1]+voxels[x][y][z+1])  ,
                        (double)(-voxels[x][y+1][z+1]+voxels[x][y+1][z])
                    }
                },
                    {
                        voxels[x][y][z],
                        voxels[x+1][y][z],
                        voxels[x+1][y][z+1],
                        voxels[x][y][z+1],
                        voxels[x][y+1][z],
                        voxels[x+1][y+1][z],
                        voxels[x+1][y+1][z+1],
                        voxels[x][y+1][z+1]
                    }
                };
                int c1,c2,c3;
                for (int i=0; i<8; i++)
                {
                    c1 = c.p[i].x;
                    c2 = c.p[i].y;
                    c3 = c.p[i].z;
                    if (Der[c3*N*N+c2*N+c1]>0) {sign = 1;break;}
                }
                if (sign == 1)
                    processCube(c, isovalue,&vertexList);
            }
        }
    }
    return vertexList;
}

