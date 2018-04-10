//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano 2018. All rights reserved.
//
// This is an abstraction on "CoherentUnit" per https://arxiv.org/abs/1710.05867
//
// "SeparatedUnit2" keeps representation of qubit states separated until explicitly
// entangled. This makes for large gains in memory and speed optimization in the
// best case scenario. "CoherentUnit" has been optimized for the worst case scenario.
//
// Licensed under the GNU General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/gpl-3.0.en.html
// for details.

#pragma once

#include "qregister.hpp"
#include <vector>

namespace Qrack {

struct QuantumBitShard {
    std::shared_ptr<CohesiveUnit *> unit;
    bitLenInt mapped;
};

class SeparatedUnit2;

class SeparatedUnit2 { // : public IQuantumComputer
public:
    SeparatedUnit(CoherentUnitEngine engine, bitLenInt qBitCount, bitCapInt initState, Complex16 phaseFac, uint32_t rand_seed);

protected:
    CoherentUnitEngine engine;
    std::vector<QuantumBitShard> shards;
    std::default_random_engine rand_generator_ptr;

    void Decompose(bitLenInt qubit);

    typedef void (CohesiveUnit::*TwoBitCall(bitLenInt, bitLenInt);
    typedef void (CohesiveUnit::*ThreeBitCall(bitLenInt, bitLenInt, bitLenInt);
    void EntangleAndCall(bitLenInt bit1, bitLenInt bit2, TwoBitCall fn);
    void EntangleAndCall(bitLenInt bit1, bitLenInt bit2, bitLenInt bit3, ThreeBitCall fn);

public:
    /* Common interfaces with CU below here. */
    double Prob(bitLenInt qubitIndex);
    double ProbAll(bitCapInt perm);
    void ProbArray(double* probArray);
    bool M(bitLenInt qubitIndex);
    bitCapInt MReg(bitLenInt start, bitLenInt length);
    void SetBit(bitLenInt qubitIndex1, bool value);
    void SetReg(bitLenInt start, bitLenInt length, bitCapInt value);
    void SetPermutation(bitCapInt value);

    void Swap(bitLenInt qubitIndex1, bitLenInt qubitIndex2);
    void Swap(bitLenInt qubitIndex1, bitLenInt qubitIndex2, bitLenInt length);

    void AND(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit);
    void OR(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit);
    void XOR(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit);
    void CLAND(bitLenInt inputQBit, bool inputClassicalBit, bitLenInt outputBit);
    void CLOR(bitLenInt inputQBit, bool inputClassicalBit, bitLenInt outputBit);
    void CLXOR(bitLenInt inputQBit, bool inputClassicalBit, bitLenInt outputBit);

    void CCNOT(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit);
    void AntiCCNOT(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit);

    void H(bitLenInt qubitIndex);
    void X(bitLenInt qubitIndex);
    void Y(bitLenInt qubitIndex);
    void Z(bitLenInt qubitIndex);

    void X(bitLenInt start, bitLenInt length);

    void CY(bitLenInt control, bitLenInt target);
    void CZ(bitLenInt control, bitLenInt target);

    void RT(double radians, bitLenInt qubitIndex);
    void RTDyad(int numerator, int denominator, bitLenInt qubitIndex);
    void RX(double radians, bitLenInt qubitIndex);
    void RXDyad(int numerator, int denominator, bitLenInt qubitIndex);
    void RY(double radians, bitLenInt qubitIndex);
    void RYDyad(int numerator, int denominator, bitLenInt qubitIndex);
    void RZ(double radians, bitLenInt qubitIndex);
    void RZDyad(int numerator, int denominator, bitLenInt qubitIndex);

    void CRT(double radians, bitLenInt control, bitLenInt target);
    void CRTDyad(int numerator, int denominator, bitLenInt control, bitLenInt target);
    void CRY(double radians, bitLenInt control, bitLenInt target);
    void CRYDyad(int numerator, int denominator, bitLenInt control, bitLenInt target);
    void CRZ(double radians, bitLenInt control, bitLenInt target);
    void CRZDyad(int numerator, int denominator, bitLenInt control, bitLenInt target);

    void ROL(bitLenInt shift, bitLenInt start, bitLenInt length);
    void ROR(bitLenInt shift, bitLenInt start, bitLenInt length);

    void INC(bitCapInt toAdd, bitLenInt start, bitLenInt length);
    void INCC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex);
    void INCS(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt overflowIndex);
    void INCSC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt overflowIndex, bitLenInt carryIndex);
    void INCSC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex);
    void INCBCD(bitCapInt toAdd, bitLenInt start, bitLenInt length);
    void INCBCDC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex);
    void DEC(bitCapInt toSub, bitLenInt start, bitLenInt length);
    void DECC(bitCapInt toSub, bitLenInt start, bitLenInt length, bitLenInt carryIndex);
    void DECS(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt overflowIndex);
    void DECSC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt overflowIndex, bitLenInt carryIndex);
    void DECSC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex);
    void DECBCD(bitCapInt toAdd, bitLenInt start, bitLenInt length);
    void DECBCDC(bitCapInt toSub, bitLenInt start, bitLenInt length, bitLenInt carryIndex);

    void QFT(bitLenInt start, bitLenInt length);
    void ZeroPhaseFlip(bitLenInt start, bitLenInt length);
    void CPhaseFlipIfLess(bitCapInt greaterPerm, bitLenInt start, bitLenInt length, bitLenInt flagIndex);
    void PhaseFlip();

    unsigned char SuperposeReg8(bitLenInt inputStart, bitLenInt outputStart, unsigned char* values);
    unsigned char AdcSuperposeReg8(
        bitLenInt inputStart, bitLenInt outputStart, bitLenInt carryIndex, unsigned char* values);
    unsigned char SbcSuperposeReg8(
        bitLenInt inputStart, bitLenInt outputStart, bitLenInt carryIndex, unsigned char* values);

} // namespace Qrack
