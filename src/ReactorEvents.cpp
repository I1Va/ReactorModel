#include "ReactorEvents.h"
#include "ReactorWall.h"

#include <limits>

static const double BOOM_DT = 0.001;
static const double BOOM_SPEED_COEF = 100;

MoleculeReactionEvent detectMoleculeCollision(const double endDelta, Molecule *fstMolecule, Molecule *sndMolecule) {
    assert(fstMolecule);
    assert(sndMolecule);

    if (!fstMolecule->isAlive() || !sndMolecule->isAlive()) return MoleculeReactionEvent::POISON();

    double coliisionRadius = fstMolecule->getCollideCircleRadius() + sndMolecule->getCollideCircleRadius();
    double coliisionRadius2 = coliisionRadius * coliisionRadius;

    gm_vector<double, 2> V = fstMolecule->getSpeedVector() - sndMolecule->getSpeedVector();
    gm_vector<double, 2> P = fstMolecule->getPosition() - sndMolecule->getPosition();

    double t1 = 0, t2 = 0;
    int nRoots = 0;
    double aCoef = V.get_x() * V.get_x() + V.get_y() * V.get_y();
    double bCoef = 2 * (P.get_x() * V.get_x() + P.get_y() * V.get_y());
    double cCoef = P.get_x() * P.get_x() + P.get_y() * P.get_y() - coliisionRadius2;

    solveQuadratic(aCoef, bCoef, cCoef, &t1, &t2, &nRoots);

    if (nRoots != 2) return MoleculeReactionEvent::POISON();
     
    if (t1 < 0 || t1 > endDelta) return MoleculeReactionEvent::POISON();


    return MoleculeReactionEvent(t1, endDelta, fstMolecule, sndMolecule);
}

WallCollisionEvent detectWallCollision(const double endDelta, Molecule *molecule, ReactorWall reactorWalls[reactorWallsCnt]) {
    assert(molecule);

    if (!molecule->isAlive()) return WallCollisionEvent::POISON();

    
    double collideRadius = molecule->getCollideCircleRadius();
    double speedX = molecule->getSpeedVector().get_x();
    double speedY = molecule->getSpeedVector().get_y();
    double posX = molecule->getPosition().get_x();
    double posY = molecule->getPosition().get_y();

    if (std::fabs(speedX) < std::numeric_limits<double>::epsilon() && 
        std::fabs(speedY) < std::numeric_limits<double>::epsilon()) {
            return WallCollisionEvent::POISON();
        }
    

    double tMin = INFINITY;
    ReactorWallTypes collisionWallType = NONE_WALL;

    if (std::fabs(speedX) >= std::numeric_limits<double>::epsilon()) {
        if (speedX > 0) { // right wall
            double tx = (reactorWalls[RIGHT_WALL].measure - collideRadius - posX) / speedX;
            if (tx < tMin) {
                tMin = tx;
                collisionWallType = RIGHT_WALL;
            }
        } else { // left wall
            double tx = (collideRadius - posX) / speedX; 
            if (tx < tMin) {
                tMin = tx;
                collisionWallType = LEFT_WALL;
            }
        }
    }
    if (std::fabs(speedX) >= std::numeric_limits<double>::epsilon()) {
        if (speedY > 0) { // bottom wall
            double ty = (reactorWalls[BOTTOM_WALL].measure - collideRadius - posY) / speedY;
            if (ty < tMin) {
                tMin = ty;
                collisionWallType = BOTTOM_WALL;
            }
        } else {// top wall
            double ty = (collideRadius - posY) / speedY;
            if (ty < tMin) {
                tMin = ty;
                collisionWallType = TOP_WALL;
            } 
        }                 
    }

    if (std::isinf(tMin) || tMin > endDelta || collisionWallType == NONE_WALL) return WallCollisionEvent::POISON();


    ReactorWall *collisionWall = &reactorWalls[collisionWallType];
    double moleculEnergy = molecule->getKinecticEnergy();
    double wallEnergy = collisionWall->energy;

    double targetEnergy = (1 - collisionWall->energyTransferCoef) * moleculEnergy + wallEnergy * collisionWall->energyTransferCoef;
    double deltaEnergy = targetEnergy - moleculEnergy;

    double newWallEnergy = collisionWall->energy - deltaEnergy;
    double newMoleculeEnergy = moleculEnergy + deltaEnergy;

    gm_vector<double, 2> newSpeed = {};
    switch (collisionWallType)
    {
        case RIGHT_WALL:
        case LEFT_WALL:
            newSpeed = {-speedX, speedY};
            break;

        case BOTTOM_WALL:
        case TOP_WALL:
            newSpeed = {speedX, -speedY};
            break;
        default:
            assert(0 && "None WallType");
    }        

    newSpeed = newSpeed * (1 / std::sqrt(moleculEnergy)) * std::sqrt(newMoleculeEnergy); 
    gm_vector<double, 2> newPostion = molecule->getPosition() + newSpeed * (endDelta - tMin);


    return WallCollisionEvent(tMin, endDelta, molecule, newPostion, newSpeed, collisionWall, newWallEnergy);
} 

inline gm_vector<double, 2> getCollideCenter(Molecule *fstMolecule, Molecule *sndMolecule) {
    double R1 = fstMolecule->getCollideCircleRadius();
    double R2 = sndMolecule->getCollideCircleRadius();
    gm_vector<double, 2> Pos1 = fstMolecule->getPosition();
    gm_vector<double, 2> Pos2 = sndMolecule->getPosition();

    return (Pos1 * R2 + Pos2 * R1) * (1.0 / (R1 + R2));
}

void QuadritQuadritReaction(
    std::vector<Molecule*> &molecules,
    Molecule *fstMolecule,
    Molecule *sndMolecule,
    double duration
) {

    gm_vector<double, 2> collideCenter = getCollideCenter(fstMolecule, sndMolecule);

    int boomMoleculeCnt = fstMolecule->getMass() + sndMolecule->getMass();
    double summaryEnergy = fstMolecule->getPotentialEnergy() + fstMolecule->getKinecticEnergy() +
                           sndMolecule->getPotentialEnergy() + sndMolecule->getKinecticEnergy();
    
    double speedCoef = std::sqrt(2 * summaryEnergy / (boomMoleculeCnt * /*mass=*/1));

    double boomRootationAngle = 2 * std::numbers::pi / boomMoleculeCnt;
    
    double boomRadius = std::sqrt(2 / std::sin(boomRootationAngle)) * 2;
    gm_vector<double, 2> boomCurDestVector = (gm_vector<double, 2>(0, -1));

    fstMolecule->setPhyState(DEATH);
    sndMolecule->setPhyState(DEATH);

    for (int i = 0; i < boomMoleculeCnt; i++) {
        Molecule *newMolecule = (Molecule *) new Circlit(collideCenter + boomCurDestVector * duration, boomCurDestVector * speedCoef, 1);
        molecules.push_back(newMolecule);
        boomCurDestVector = boomCurDestVector.rotate(boomRootationAngle);
    }
}

void CirclitQuadritReaction(
    std::vector<Molecule*> &molecules,
    Molecule *fstMolecule,
    Molecule *sndMolecule,
    double duration
) {
    gm_vector<double, 2> collideCenter = getCollideCenter(fstMolecule, sndMolecule);

    double summaryEnergy = fstMolecule->getPotentialEnergy() + fstMolecule->getKinecticEnergy() +
                           sndMolecule->getPotentialEnergy() + sndMolecule->getKinecticEnergy();
    
    int newMass = fstMolecule->getMass() + sndMolecule->getMass();
   
    
    gm_vector<double, 2> newspeedVector = (fstMolecule->getSpeedVector() * fstMolecule->getMass() + 
                                          sndMolecule->getSpeedVector() * sndMolecule->getMass()) * (1.0 / newMass);

    double newPotentialEnergy = std::max(0.0, summaryEnergy - newspeedVector.get_len2() * newMass / 2);
    
    fstMolecule->setPhyState(DEATH);
    sndMolecule->setPhyState(DEATH);

    Molecule *newMolecule = (Molecule *) new Quadrit(collideCenter + newspeedVector * duration, newspeedVector, newMass);
    newMolecule->setPotentialEnergy(newPotentialEnergy);
   
    
    molecules.push_back(newMolecule);
}