
#ifndef MONOMER_H
#define MONOMER_H

#include <cstdint>
#include <array>
/*
State 0: 000
State 1: 001
State 2: 010
State 3: 011
State 4: 100
*/
struct Monomer {
    uint8_t data;
    uint16_t index;

    inline void setState(uint8_t state) {
        data = (data & 0b11111000) | (state & 0b00000111);
    }

    inline int getState() const {
        return static_cast<int>(data & 0b00000111);
    }

    inline void setIndex(uint16_t index_){
        index = index_;
    }
    inline uint16_t getIndex() const{
        return index;
    }
};

#endif