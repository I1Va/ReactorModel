#ifndef MOLECULE_H
#define MOLECULE_H

#include "gm_primitives.hpp"

const double SQRT2 = std::sqrt(2);


enum MoleculeTypes {
    NONE_MOL_TYPE = -1,
    CIRCLIT = 0,
    QUADRIT = 1,    
};

enum MoleculePhysicalStates {
    NONE_MOL_PHYSTATE,
    DEATH,
    ALIVE,
};


class Molecule {
protected:
    MoleculeTypes type_ = NONE_MOL_TYPE;
    MoleculePhysicalStates phyState_ = NONE_MOL_PHYSTATE;

    gm_vector<double, 2> position_ = {};
    gm_vector<double, 2> speedVector_ = {};
    int mass_ = 0;
    double size_ = 0;
public:
    Molecule
    (
        const gm_vector<double, 2> &position, 
        const gm_vector<double, 2> &speedVector, 
        const int mass
    ): position_(position), speedVector_(speedVector), mass_(mass) {}

    ~Molecule() = default;

    virtual double getCollideCircleRadius() { return 0; }
    
    bool isAlive() { return (phyState_ == MoleculePhysicalStates::ALIVE); }

    MoleculeTypes getType() const { return type_; }
    MoleculePhysicalStates getPhyState() const { return phyState_; }
    gm_vector<double, 2> getPosition() const { return position_; }
    gm_vector<double, 2> getSpeedVector() const { return speedVector_; }
    
    int getMass() const { return mass_; }
    virtual double massToSize() const { return mass_; }
    double getSize() const { return size_; }
    
    void setPosition(const gm_vector<double, 2> &other) { position_ = other; }
    void setSpeedVector(const gm_vector<double, 2> &other) { speedVector_ = other; }
    void setPhyState(const MoleculePhysicalStates state) { phyState_ = state; }
};

class Circlit : public Molecule {
public:    
    Circlit
    (
        const gm_vector<double, 2> &position, 
        const gm_vector<double, 2> &speedVector, 
        const int mass
    ): 
        Molecule(position, speedVector, mass)
    { 
        size_ = massToSize();
        type_ = MoleculeTypes::CIRCLIT;
        phyState_ = MoleculePhysicalStates::ALIVE;
    }

    ~Circlit() = default;

    double massToSize() const override { return mass_; }
    double getCollideCircleRadius() override { return size_; }
};


class Quadrit : public Molecule {
public:    
    Quadrit
    (
        const gm_vector<double, 2> &position, 
        const gm_vector<double, 2> &speedVector, 
        const int mass
    ): 
        Molecule(position, speedVector, mass)
    { 
        type_ = MoleculeTypes::QUADRIT;
        phyState_ = MoleculePhysicalStates::ALIVE;
        size_ = massToSize();
    }

    ~Quadrit() = default;

    double massToSize() const override { return mass_; }
    double getCollideCircleRadius() override { return size_ / SQRT2; }
};

#endif // MOLECULE_H