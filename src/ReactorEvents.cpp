#include "ReactorEvents.h"
#include <limits>

static const double BOOM_DT = 0.001;
static const double BOOM_SPEED_COEF = 100;

MoleculeReactionEvent tryMoleculeReactionEvent(Molecule *fstMolecule, Molecule *sndMolecule) {
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
     
    if (t1 < 0) return MoleculeReactionEvent::POISON();

    return MoleculeReactionEvent(t1, fstMolecule, sndMolecule);
}

WallCollisionEvent tryWallCollisionEvent(Molecule *molecule, const double width, const double height) {
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
    
        

    double tx = INFINITY;
    double ty = INFINITY;


    if (std::fabs(speedX) >= std::numeric_limits<double>::epsilon()) {
         if (speedX > 0) // right wall
            tx = (width - collideRadius - posX) / speedX; 
        else // left wall
            tx = (collideRadius - posX) / speedX; 
    }

    if (std::fabs(speedX) >= std::numeric_limits<double>::epsilon()) {
        if (speedY > 0) // bottom wall
            ty = (height - collideRadius - posY) / speedY; 
        else // top wall
            ty = (collideRadius - posY) / speedY;                 
    }

    double tMin = std::min(tx, ty);

    if (std::isinf(tMin)) return WallCollisionEvent::POISON();

    if (tMin == tx) {
        return WallCollisionEvent(tMin, molecule, {-speedX, speedY});
    } else {
        return WallCollisionEvent(tMin, molecule, {speedX, -speedY});
    }
} 

void handleReactorEvent(const WallCollisionEvent &event) {
    assert(event.molecule_);
    event.molecule_->setSpeedVector(event.newSpeedVector_);
}

void handleReactorEvent(const MoleculeReactionEvent &event, std::vector<Molecule *> &molecules) {   
    MoleculeReaction moleculeReaction = moleculeReactionsVTable[event.fstMolecule_->getType()][event.sndMolecule_->getType()];
    moleculeReaction(molecules, event.fstMolecule_, event.sndMolecule_);
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
    Molecule *sndMolecule
) {

    gm_vector<double, 2> collideCenter = getCollideCenter(fstMolecule, sndMolecule);

    int boomMoleculeCnt = fstMolecule->getMass() + sndMolecule->getMass();
    double boomRootationAngle = 2 * std::numbers::pi / boomMoleculeCnt;
    
    double boomRadius = std::sqrt(2 / std::sin(boomRootationAngle)) * 2;
    gm_vector<double, 2> boomCurOffsetVector = (gm_vector<double, 2>(0, -1)) * boomRadius;

    fstMolecule->setPhyState(DEATH);
    sndMolecule->setPhyState(DEATH);

    for (int i = 0; i < boomMoleculeCnt; i++) {
        Molecule *newMolecule = (Molecule *) new Circlit(collideCenter + boomCurOffsetVector * BOOM_DT, boomCurOffsetVector * BOOM_SPEED_COEF, 1);
        molecules.push_back(newMolecule);
        boomCurOffsetVector = boomCurOffsetVector.rotate(boomRootationAngle);
    }
}





void CirclitQuadritReaction(
    std::vector<Molecule*> &molecules,
    Molecule *fstMolecule,
    Molecule *sndMolecule
) {
    gm_vector<double, 2> collideCenter = getCollideCenter(fstMolecule, sndMolecule);


    int newMass = fstMolecule->getMass() + sndMolecule->getMass();
   
    
    gm_vector<double, 2> newspeedVector = (fstMolecule->getSpeedVector() * fstMolecule->getMass() + 
                                          sndMolecule->getSpeedVector() * sndMolecule->getMass()) * (1.0 / newMass);

    
    fstMolecule->setPhyState(DEATH);
    sndMolecule->setPhyState(DEATH);

    Molecule *newMolecule = (Molecule *) new Quadrit(collideCenter, newspeedVector, newMass);
    
    molecules.push_back(newMolecule);
}