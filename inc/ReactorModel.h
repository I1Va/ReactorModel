#ifndef REACTOR_MODEL_H
#define REACTOR_MODEL_H

#include <random>
#include <optional>
#include <algorithm>
#include <functional>

#include "Molecule.h"
#include "ReactorEvents.h"
#include "gm_primitives.hpp"

static const double minSpeed = 500;
static const double maxSpeed = 1000;
static const int initMass = 3;

static const double WALL_OFFSET_EPS = 0.01;

class ReactorPreUpdateState {
    std::vector<Molecule *> newMolecules_ = {};

    bool newSizeState_ = false;
    double newWidth_ = 0;
    double newHeight_ = 0;

public:
    ReactorPreUpdateState() = default;
    ~ReactorPreUpdateState() = default;

    void reset() {
        newMolecules_.clear();
        newSizeState_ = false;
        newWidth_ = 0;
        newHeight_ = 0;
    }

    void addNewMolecule(Molecule *newMolecule) { 
        assert(newMolecule);
        newMolecules_.push_back(newMolecule); 
    }
    
    std::vector<Molecule *> &getNewMolecules() {
        return newMolecules_;
    }

    void pourNewMolecules(std::vector<Molecule *> &destination) {
        destination.insert(destination.end(), newMolecules_.begin(), newMolecules_.end());
        newMolecules_.clear();        
    }

    void setNewSize(const double newWidth, const double newHeight) {
        newWidth_ = newWidth;
        newHeight_ = newHeight;
        newSizeState_ = true;
    }

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
        if (onUpdate_) onUpdate_();
    
        preUpdate();

        resolveInvalides();

        processReactortEvent(deltaSecs);

        resolveInvalides();
    }
    

    // USER API
    const std::vector<Molecule*> &getMolecules() const { return molecules_; }
    
    void addCirclit() {
        Molecule *newMolecule = genNewMolecule(MoleculeTypes::CIRCLIT); 
        preUpdateState_.addNewMolecule(newMolecule);
    }

    void addQuadrit() {
        Molecule *newMolecule = genNewMolecule(MoleculeTypes::QUADRIT); 
        preUpdateState_.addNewMolecule(newMolecule);
    }

    void resize(const double newWidth, const double newHeight) {preUpdateState_.setNewSize(newWidth, newHeight); }

private:
    double randRange(double start, double end) {
        return start + (end - start) * randomGenerator_() / double(randomGenerator_.max());
    }
    
    Molecule *genNewMolecule(MoleculeTypes type) {
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
        
        return newMolecule;
    }

    void reorderMolecules() {
        int left = 0;
        int right = (int)( molecules_.size()) - 1;
        while (right - left > 0) {
            if (!molecules_[left]->isAlive() && molecules_[right]->isAlive()) {
                std::swap(molecules_[left], molecules_[right]);
            }
            if (molecules_[left]->isAlive()) {
                left++;
            }
            if (!molecules_[right]->isAlive()) {
                right--;
            }
        }
    }

    void preUpdate() {
        preUpdateState_.pourNewMolecules(molecules_);
 
        if (preUpdateState_.getNewSizeState()) {
            width_ = preUpdateState_.getNewWidth();
            height_ = preUpdateState_.getNewHeight();
        }

        preUpdateState_.reset();
    }

    void clearDeathMolecules() {
        reorderMolecules();

        size_t aliveCnt = 0;
        for (auto molecule : molecules_) {
            aliveCnt += (molecule->isAlive());
        }


        for (auto it = molecules_.begin() + aliveCnt; it != molecules_.end();) {
            delete (*it);
            it = molecules_.erase(it);
        }
    }

    bool isMoleculeValid(const Molecule &molecule) {
        double collideRadius = molecule.getCollideCircleRadius();
        
        if (!molecule.getPosition().is_valid()) return false;
        if (!molecule.getSpeedVector().is_valid()) return false;
        if ((0 + collideRadius > width_ - collideRadius) || (0 + collideRadius > height_ - collideRadius)) return false;

        return true;
    }

    void resolveInvalides() {
        for (auto molecule : molecules_) {
            if (!isMoleculeValid(*molecule)) {
                molecule->setPhyState(MoleculePhysicalStates::DEATH);
            }
        }

        clearDeathMolecules();
    
        for (auto molecule : molecules_) {
            double x = molecule->getPosition().get_x();
            double y = molecule->getPosition().get_y();
            double r = molecule->getCollideCircleRadius();

            if (x < r) x += WALL_OFFSET_EPS;
            if (x > width_ - r) x -= WALL_OFFSET_EPS;
            if (y < r) y += WALL_OFFSET_EPS;
            if (y > height_ - r) y -= WALL_OFFSET_EPS;    

            molecule->setPosition({x, y});
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

        
        if (closestEventDelta > deltaSecs) {
            processReactorStableState(deltaSecs);
            return;
        }

        processReactorStableState(closestEventDelta);
        
        // need to process several events
        if (wallCollisionEventState) {
            handleReactorEvent(closestWallCollisionEvent);
        } else {
            handleReactorEvent(closestoleculeReactionEvent, preUpdateState_.getNewMolecules());
        }

        // processReactortEvent(deltaSecs - closestEventDelta);
    }
};


#endif // REACTOR_MODEL_H