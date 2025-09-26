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
static const int initMass = 1;
static const double WALL_OFFSET_EPS = 1;

class ReactorPreUpdateState {
    std::vector<Molecule *> newMolecules_ = {};

    bool newSizeState_ = false;
    double newWidth_ = 0;
    double newHeight_ = 0;
    bool removeMoleculeState_ = false;

public:
    ReactorPreUpdateState() = default;
    ~ReactorPreUpdateState() = default;

    void reset() {
        newMolecules_.clear();
        newSizeState_ = false;
        newWidth_ = 0;
        newHeight_ = 0;
        removeMoleculeState_ = false;
    }

    void removeMolecule() {
        removeMoleculeState_ = true;
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

    bool getRemoveMoleculeState() const { return removeMoleculeState_; }
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
    
    void removeMolecule() {
        preUpdateState_.removeMolecule();
    }

    void narrowRightWall(const double delta) {
        double newWidth = std::max(width_ - delta, 0.0);
        preUpdateState_.setNewSize(newWidth, height_);
    }

    void addCirclit() {
        Molecule *newMolecule = genNewMolecule(MoleculeTypes::CIRCLIT); 
        preUpdateState_.addNewMolecule(newMolecule);
    }

    void addQuadrit() {
        Molecule *newMolecule = genNewMolecule(MoleculeTypes::QUADRIT); 
        preUpdateState_.addNewMolecule(newMolecule);
    }


    double getWidth() const { return width_; }
    double getHeight() const { return height_; }

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
        if (preUpdateState_.getRemoveMoleculeState()) {
            if (!molecules_.size()) return;
    
            molecules_[(int) randRange(0, molecules_.size() - 1)]->setPhyState(MoleculePhysicalStates::DEATH);
        }

        preUpdateState_.reset();
    }

    void reviveUnresponsives() {
        for (auto &molecule : molecules_) {
            if (molecule->getPhyState() == MoleculePhysicalStates::UNRESPONSIVE)
                molecule->setPhyState(MoleculePhysicalStates::ALIVE);
        }
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

    void resolveInvalides() {
        for (auto molecule : molecules_) {
            double collideRadius = molecule->getCollideCircleRadius();
            if ((0 + collideRadius > width_ - collideRadius) || (0 + collideRadius > height_ - collideRadius)) {
                molecule->setPhyState(MoleculePhysicalStates::DEATH);
            }
        }
    
        reviveUnresponsives();
        clearDeathMolecules();

        for (auto molecule : molecules_) {
            double x = molecule->getPosition().get_x();
            double y = molecule->getPosition().get_y();
            double r = molecule->getCollideCircleRadius();
            
            if (x > width_)
                x = width_ - r - WALL_OFFSET_EPS;
            if (y > height_)
                y = height_ - r - WALL_OFFSET_EPS;

            molecule->setPosition({x, y});
        }
    }
   
    void processReactorStableState(const double deltaSecs) {
        for (auto &molecule : molecules_)
            molecule->setPosition(molecule->getPosition() + molecule->getSpeedVector() * deltaSecs);
    }

    void processReactortEvent(const double deltaSecs) {
        std::vector<ReactorEvent *> events = {};

        
        for (size_t fstMoleculeId = 0; fstMoleculeId < molecules_.size(); fstMoleculeId++) {
            WallCollisionEvent curWallCollisionEvent = detectWallCollision(deltaSecs, molecules_[fstMoleculeId], width_, height_);
            if (!curWallCollisionEvent.isPoison()) {
                ReactorEvent *newEvent = (ReactorEvent *) new WallCollisionEvent(curWallCollisionEvent);
                events.push_back(newEvent);
            }

            for (size_t sndMoleculeId = 0; sndMoleculeId < molecules_.size(); sndMoleculeId++) {
                MoleculeReactionEvent curMoleculeReactionEvent = detectMoleculeCollision(deltaSecs, molecules_[fstMoleculeId], molecules_[sndMoleculeId]);
                
                if (!curMoleculeReactionEvent.isPoison()) {
                    ReactorEvent *newEvent = (ReactorEvent *) new MoleculeReactionEvent(curMoleculeReactionEvent);
                    events.push_back(newEvent);
                }
            }
        }

        if (!events.size())  {
            processReactorStableState(deltaSecs);
            return;
        }

        auto cmpByTime = [](ReactorEvent* a, ReactorEvent* b) {
            return a->getStartDelta() < b->getStartDelta();
        };

        std::sort(events.begin(), events.end(), cmpByTime);
        
        for (auto event : events) {
            event->handleReactorEvent(preUpdateState_.getNewMolecules());
            delete event;
        }
    }
};


#endif // REACTOR_MODEL_H