#ifndef REACTOR_WALL_H
#define REACTOR_WALL_H


const size_t reactorWallsCnt = 4;

enum ReactorWallTypes {
    NONE_WALL = -1,
    RIGHT_WALL = 0, 
    LEFT_WALL = 1, 
    BOTTOM_WALL = 2,
    TOP_WALL = 3, 
};

struct ReactorWall {    
    double measure = 0;
    double energy = 0;
    double energyTransferCoef = 0;
};


#endif // REACTOR_WALL_H