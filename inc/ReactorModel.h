#ifndef REACTOR_MODEL_H
#define REACTOR_MODEL_H

#include <random>
#include <optional>
#include <algorithm>
#include <functional>

#include "Molecule.h"
#include "ReactorEvents.h"
#include "gm_primitives.hpp"

const double minSpeed = 1;
const double maxSpeed = 10;
const int initMass = 1;

class ReactorPreUpdateState {
    int newCirclitsCnt_ = 0;
    int newQuadritsCnt_ = 0;

    bool newSizeState_ = false;
    double newWidth_ = 0;
    double newHeight_ = 0;

public:
    ReactorPreUpdateState() = default;
    ~ReactorPreUpdateState() = default;

    void reset() {
        newCirclitsCnt_ = 0;
        newQuadritsCnt_ = 0;
        newSizeState_ = false;
        newWidth_ = 0;
        newHeight_ = 0;
    }

    void addCirclit() { newCirclitsCnt_++; }
    void addQuadrit() { newQuadritsCnt_++; }
    void setNewSize(const double newWidth, const double newHeight) {
        newWidth_ = newWidth;
        newHeight_ = newHeight;
        newSizeState_ = true;
    }

    int getNewCirclitsCnt() const { return newCirclitsCnt_; }
    int getNewQuadritsCnt() const { return newQuadritsCnt_; }
    bool getNewSizeState() const { return newSizeState_;}
    double getNewWidth() const { return newWidth_; }
    double getNewHeight() const { return newHeight_; }
};

class ReactorModel {
    std::mt19937 randomGenerator_ = {};

    double width_ = 0;
    double height_ = 0;

    std::vector<Molecule *> molecules_ = {};
    ReactorPreUpdateState preUpdateState_ = {};
    std::function<void()> onUpdate_ = nullptr;

public:
    ReactorModel
    (       
        const double width, const double height,
        std::function<void()> onUpdate=nullptr,
        std::optional<unsigned int> seed=std::nullopt
    ):
        width_(width), height_(height),
        onUpdate_(onUpdate),
        randomGenerator_(seed.value_or(std::random_device{}()))
    {}

    void setOnUpdate(std::function<void()> onUpdate) {
        onUpdate_ = onUpdate;
    }

    ReactorModel() = default;

    ~ReactorModel() = default;


    void update(const double deltaSecs) {
        preUpdate();

        resolveInvalides();
        
        processReactortEvent(deltaSecs);
        
        if (onUpdate_) onUpdate_();
    }
    

    // USER API
    const std::vector<Molecule*> &getMolecules() const { return molecules_; }
    
    void addCirclit() { preUpdateState_.addCirclit(); }
    void addQuadrit() { preUpdateState_.addQuadrit(); }

    void resize(const double newWidth, const double newHeight) {preUpdateState_.setNewSize(newWidth, newHeight); }

private:
    double randRange(double start, double end) {
        return start + (end - start) * randomGenerator_() / double(randomGenerator_.max());
    }
    
    void addMolecule(MoleculeTypes type) {
        Molecule *newMolecule = nullptr;

        double posX = randRange(width_ / 3, 2 * width_ / 3);
        double posY = randRange(height_ / 3, 2 * height_ / 3);

        double randomAngle = randRange(0, 2 * std::numbers::pi);
        gm_vector<double, 2> SpeedVector = (gm_vector<double, 2>(0, 1) * randRange(minSpeed, maxSpeed)).rotate(randomAngle);

        switch (type) {
            case MoleculeTypes::CIRCLIT:
                newMolecule = (Molecule *) new Circlit({posX, posY}, SpeedVector, initMass);
                break;
            case MoleculeTypes::QUADRIT:
                newMolecule = (Molecule *) new Quadrit({posX, posY}, SpeedVector, initMass);
                break;
            default:
                assert(newMolecule);
                break;
        }
        
        
        molecules_.push_back(newMolecule);
    }

    void reorderMolecules() {
        int curMoleculeIdx = molecules_.size() - 1;
        while (curMoleculeIdx - 1 >= 0) {
            Molecule *leftMolecule = molecules_[curMoleculeIdx - 1];
            Molecule *rightMolecule = molecules_[curMoleculeIdx];

            if (!leftMolecule->isAlive() && leftMolecule->isAlive()) {
                std::swap(molecules_[curMoleculeIdx - 1], molecules_[curMoleculeIdx]);
            }
            curMoleculeIdx--;
        }
    }

    void preUpdate() {
        for (int i = 0; i < preUpdateState_.getNewCirclitsCnt(); i++) addMolecule(MoleculeTypes::CIRCLIT);
        for (int i = 0; i < preUpdateState_.getNewQuadritsCnt(); i++) addMolecule(MoleculeTypes::QUADRIT);
 
        if (preUpdateState_.getNewSizeState()) {
            width_ = preUpdateState_.getNewWidth();
            height_ = preUpdateState_.getNewHeight();
        }

        preUpdateState_.reset();
    }

    void clearDeathMolecules() {
        reorderMolecules();
        int aliveMoleculesCnt = 0;
        for (auto molecule : molecules_) {
            if (molecule->isAlive()) {
                aliveMoleculesCnt += molecule->isAlive();
            } else {
                delete molecule;
            }
        }
        molecules_.resize(aliveMoleculesCnt); 
    }

    void resolveInvalides() {
        for (auto molecule : molecules_) {
            double collideRadius = molecule->getCollideCircleRadius();
            if ((0 + collideRadius > width_ - collideRadius) || (0 + collideRadius > height_ - collideRadius)) {
                molecule->setPhyState(MoleculePhysicalStates::DEATH);
            }
        }

        clearDeathMolecules();
    
        for (auto molecule : molecules_) {
            double collideRadius = molecule->getCollideCircleRadius();
            double clampedX = std::clamp(molecule->getPosition().get_x(), 0 + collideRadius, width_ - collideRadius);
            double clampedY = std::clamp(molecule->getPosition().get_y(), 0 + collideRadius, height_ - collideRadius);
            molecule->setPosition({clampedX, clampedY});
        }
    }
   
    void processReactorStableState(const double deltaSecs) {
        for (auto &molecule : molecules_)
            molecule->setPosition(molecule->getPosition() + molecule->getSpeedVector() * deltaSecs);
    }

    void processReactortEvent(const double deltaSecs) {
        WallCollisionEvent closestWallCollisionEvent = WallCollisionEvent::POISON();
        MoleculeReactionEvent closestoleculeReactionEvent = MoleculeReactionEvent::POISON();
       
        for (size_t fstMoleculeId = 0; fstMoleculeId < molecules_.size(); fstMoleculeId++) {
            WallCollisionEvent curWallCollisionEvent = tryWallCollisionEvent(molecules_[fstMoleculeId], width_, height_);
            if (closestWallCollisionEvent.isPoison() || (!curWallCollisionEvent.isPoison() && curWallCollisionEvent < closestWallCollisionEvent))
                closestWallCollisionEvent = curWallCollisionEvent;

            for (size_t sndMoleculeId = 0; sndMoleculeId < molecules_.size(); sndMoleculeId++) {
                MoleculeReactionEvent curMoleculeReactionEvent = tryMoleculeReactionEvent(molecules_[fstMoleculeId], molecules_[sndMoleculeId]);
                if (closestoleculeReactionEvent.isPoison() || (!curMoleculeReactionEvent.isPoison() && curMoleculeReactionEvent < closestoleculeReactionEvent))
                    closestoleculeReactionEvent = curMoleculeReactionEvent;
            }
        }

        if (closestWallCollisionEvent.isPoison() && closestoleculeReactionEvent.isPoison()) {
            processReactorStableState(deltaSecs);
            return;
        }
        
        bool wallCollisionEventState = (!closestWallCollisionEvent.isPoison() && closestoleculeReactionEvent.isPoison()) || (closestWallCollisionEvent < closestoleculeReactionEvent);
        double closestEventDelta = wallCollisionEventState ? closestWallCollisionEvent.getDeltaSecs() : closestoleculeReactionEvent.getDeltaSecs();

        std::cout << "closestoleculeReactionEvent : " << closestoleculeReactionEvent.getDeltaSecs() << "\n";
        
        

        if (closestEventDelta > deltaSecs) {
            processReactorStableState(deltaSecs);
            return;
        }

        processReactorStableState(closestEventDelta);
        

        if (wallCollisionEventState) {
            handleReactorEvent(closestWallCollisionEvent);
        } else {
            handleReactorEvent(closestoleculeReactionEvent, molecules_);
        }
    }
};


#endif // REACTOR_MODEL_H