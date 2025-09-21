#ifndef REACTOR_EVENTS_H
#define REACTOR_EVENTS_H

#include "Molecule.h"
#include "gm_primitives.hpp"

#include <vector>


typedef void (*MoleculeReaction) (
    std::vector<Molecule*> &molecules,
    Molecule *fstMolecule,
    Molecule *sndMolecule
);


void CirclitQuadritReaction(
    std::vector<Molecule*> &molecules,
    Molecule *fstMolecule,
    Molecule *sndMolecule
);

void QuadritQuadritReaction(
    std::vector<Molecule*> &molecules,
    Molecule *fstMolecule,
    Molecule *sndMolecule
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

    double deltaSecs_; 
    ReactorEventType type_;
public:
    ReactorEvent() = default;
    ReactorEvent(const double deltaSecs, const ReactorEventType type) : deltaSecs_(deltaSecs), type_(type) {}

    virtual ~ReactorEvent() = default;

    bool isPoison() const {
        return (type_ == NoneEvent);
    }

    bool operator<(const ReactorEvent& other) const {
        return deltaSecs_ < other.deltaSecs_;
    }

    double getDeltaSecs() const { return deltaSecs_; }
};

class WallCollisionEvent : public ReactorEvent {
    Molecule* molecule_;
    gm_vector<double, 2> newSpeedVector_;    
public:
    WallCollisionEvent() = default;
    WallCollisionEvent(const double deltaSecs, Molecule * molecule, const gm_vector<double, 2>& newSpeedVector)
        : ReactorEvent(deltaSecs, ReactorEventType::WallCollision), molecule_(molecule), newSpeedVector_(newSpeedVector) {}


    static WallCollisionEvent POISON() {
        WallCollisionEvent posionEvent = {};
        posionEvent.type_ = ReactorEventType::NoneEvent;
        return posionEvent;
    }

    friend void handleReactorEvent(const WallCollisionEvent &event);
};

class MoleculeReactionEvent : public ReactorEvent {
    Molecule* fstMolecule_;
    Molecule* sndMolecule_;
public:
    MoleculeReactionEvent() = default;
    MoleculeReactionEvent(const double deltaSecs, Molecule *fstMolecule, Molecule *sndMolecule)
        : ReactorEvent(deltaSecs, ReactorEventType::Reaction), fstMolecule_(fstMolecule), sndMolecule_(sndMolecule) {}

    static MoleculeReactionEvent POISON() {
        MoleculeReactionEvent posionEvent = {};
        posionEvent.type_ = ReactorEventType::NoneEvent;
        return posionEvent;
    }

    friend void handleReactorEvent(const MoleculeReactionEvent &event, std::vector<Molecule *> &molecules);
};


MoleculeReactionEvent tryMoleculeReactionEvent(Molecule *fstMolecule, Molecule *sndMolecule);
WallCollisionEvent tryWallCollisionEvent(Molecule *molecule, const double width, const double height);

void handleReactorEvent(const WallCollisionEvent &event);
void handleReactorEvent(const MoleculeReactionEvent &event, std::vector<Molecule*> &molecules);


#endif // REACTOR_EVENTS_H