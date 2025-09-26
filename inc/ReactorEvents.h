#ifndef REACTOR_EVENTS_H
#define REACTOR_EVENTS_H

#include "Molecule.h"
#include "gm_primitives.hpp"
#include "ReactorWall.h"

#include <vector>


typedef void (*MoleculeReaction) (
    std::vector<Molecule*> &molecules,
    Molecule *fstMolecule,
    Molecule *sndMolecule,
    double duration
);


void CirclitQuadritReaction(
    std::vector<Molecule*> &molecules,
    Molecule *fstMolecule,
    Molecule *sndMolecule,
    double duration
);

void QuadritQuadritReaction(
    std::vector<Molecule*> &molecules,
    Molecule *fstMolecule,
    Molecule *sndMolecule,
    double duration
);

const MoleculeReaction CirclitCirclitReaction = CirclitQuadritReaction;

const MoleculeReaction moleculeReactionsVTable[2][2] = 
{
    {/*[0:CIRCLIT][0:CIRCLIT]=*/CirclitCirclitReaction, /*[0:CIRCLIT][1:QUADRIT]=*/ CirclitQuadritReaction},
    {/*[1:QUADRIT][0:CIRCLIT]=*/CirclitQuadritReaction, /*[1:QUADRIT][1:QUADRIT]=*/ QuadritQuadritReaction}
};


enum ReactorEventType { NoneEvent, WallCollision, Reaction };

class ReactorEvent {
friend class WallCollisionEvent;
friend class MoleculeReactionEvent;
    double startDelta_; 
    double endDelta_;
    ReactorEventType type_;
public:
    ReactorEvent() = default;
    ReactorEvent(const double startDelta, const double endDelta, const ReactorEventType type) : 
                 startDelta_(startDelta), endDelta_(endDelta), type_(type) {}

    virtual ~ReactorEvent() = default;

    bool isPoison() const {
        return (type_ == NoneEvent);
    }

    bool operator<(const ReactorEvent& other) const {
        return startDelta_ < other.startDelta_;
    }

    double getStartDelta() const { return startDelta_; }

    virtual void handleReactorEvent(std::vector<Molecule*> &molecules) {}
};

class WallCollisionEvent : public ReactorEvent {
    Molecule* molecule_;
    gm_vector<double, 2> newPosition_;    
    gm_vector<double, 2> newSpeedVector_;
    
    ReactorWall *reactorWall_;
    double newWallEnergy_;

public:
    WallCollisionEvent() = default;
    WallCollisionEvent(
        const double startDelta, const double endDelta, Molecule * molecule, 
        const gm_vector<double, 2>& newPosition, const gm_vector<double, 2>& newSpeedVector,
        ReactorWall *reactorWall, const double newWallEnergy
    ): 
        ReactorEvent(startDelta, endDelta, ReactorEventType::WallCollision), molecule_(molecule), 
        newPosition_(newPosition), newSpeedVector_(newSpeedVector),
        reactorWall_(reactorWall), newWallEnergy_(newWallEnergy) {}
    
    static WallCollisionEvent POISON() {
        WallCollisionEvent posionEvent = {};
        posionEvent.type_ = ReactorEventType::NoneEvent;
        return posionEvent;
    }

    void handleReactorEvent(std::vector<Molecule*> &molecules) override {
        if (!molecule_->isAlive()) return;

        molecule_->setSpeedVector(newSpeedVector_);
        molecule_->setPosition(newPosition_);
        molecule_->setPhyState(MoleculePhysicalStates::UNRESPONSIVE);
        reactorWall_->energy = newWallEnergy_;
    }
};


class MoleculeReactionEvent : public ReactorEvent {
    Molecule* fstMolecule_;
    Molecule* sndMolecule_;
public:
    MoleculeReactionEvent() = default;
    MoleculeReactionEvent(const double startDelta, const double endDelta, Molecule *fstMolecule, Molecule *sndMolecule)
        : ReactorEvent(startDelta, endDelta, ReactorEventType::Reaction), fstMolecule_(fstMolecule), sndMolecule_(sndMolecule) {}

    static MoleculeReactionEvent POISON() {
        MoleculeReactionEvent posionEvent = {};
        posionEvent.type_ = ReactorEventType::NoneEvent;
        return posionEvent;
    }

    void handleReactorEvent(std::vector<Molecule *> &molecules) override {
        if (!fstMolecule_->isAlive() || !sndMolecule_->isAlive()) return;

        MoleculeReaction moleculeReaction = moleculeReactionsVTable[fstMolecule_->getType()][sndMolecule_->getType()];
        moleculeReaction(molecules, fstMolecule_, sndMolecule_, endDelta_ - startDelta_);
    }
};


MoleculeReactionEvent detectMoleculeCollision(const double endDelta, Molecule *fstMolecule, Molecule *sndMolecule);
WallCollisionEvent detectWallCollision(const double endDelta, Molecule *molecule, ReactorWall reactorWalls[reactorWallsCnt]);


#endif // REACTOR_EVENTS_H