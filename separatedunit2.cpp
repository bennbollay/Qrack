//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano 2018. All rights reserved.
//
// This is an abstraction on "CoherentUnit" per https://arxiv.org/abs/1710.05867
//
// "SeparatedUnit" keeps representation of qubit states separated until explicitly
// entangled. This makes for large gains in memory and speed optimization in the
// best case scenario. "CoherentUnit" has been optimized for the worst case scenario.
//
// Licensed under the GNU General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/gpl-3.0.en.html
// for details.

#include "separatedunit.hpp"
#include <iostream>

namespace Qrack {

SeparatedUnit::SeparatedUnit(CoherentUnitEngine eng, bitLenInt qBitCount, bitCapInt initState, Complex16 phaseFac, uint32_t rand_seed) : engine(eng)
{
    rand_generator = std::default_random_engine();
    rand_generator->seed(rand_seed);
    shards.resize(qBitCount);

    if (phaseFac == Complex16(-999.0, -999.0)) {
        double angle = Rand() * 2.0 * M_PI;
        phaseFac = Complex16(cos(angle), sin(angle));
    }

    for (auto shard : shards) {
        shard.unit = CreateCoherentUnit(engine, 1, 0, phaseFac, rand_generator);
        shard.mapped = 0;
    }
}

void SeparatedUnit::Decompose(bitLenInt qubit)
{
    std::shared_ptr<CohesiveUnit> unit = shards[qubit].unit;
    for (auto shard : shards) {
        if (shard.unit == unit) {
            shard.unit = CreateCoherentUnit(engine, 1, 0, phaseFac, rand_generator);
            shard.unit->SetBit(0, unit->M(shard.mapped)); // Probably wrong, but YWKIM
            shard.mapped = 0;
        }
    }
}

void SeparatedUnit::EntangleAndCall(bitLenInt bit1, bitLenInt bit2, TwoBitCall fn)
{
    std::shared_ptr<CohesiveUnit> unit1 = shards[bit1].unit;
    std::shared_ptr<CohesiveUnit> unit2 = shards[bit2].unit;

    if (unit1 != unit2) {
        // Not already cohered; create a new unit and merge.
        unit1->Cohere(unit2);

        // Adjust all of the shards that referenced either of the old units.
        for (auto shard : shards) {
            if (shard.unit == unit2) {
                shard.unit = unit1;
                shard.mapped = shard.mapped + unit1->GetQubitCount();
            }
        }
    }

    (unit->*fn)(shards[bit1].mapped, shards[bit2].mapped);
}

void SeparatedUnit::EntangleAndCall(bitLenInt bit1, bitLenInt bit2, bitLenInt bit3, ThreeBitCall fn)
{
    std::shared_ptr<CohesiveUnit> unit1 = shards[bit1].unit;
    std::shared_ptr<CohesiveUnit> unit2 = shards[bit2].unit;
    std::shared_ptr<CohesiveUnit> unit3 = shards[bit3].unit;

    if (unit1 != unit2 || unit2 != unit3) {
        // Not already cohered; create a new unit and merge.
        unit1->Cohere(unit2);
        unit1->Cohere(unit3);

        // Adjust all of the shards that referenced either of the old units.
        for (auto shard : shards) {
            if (shard.unit == unit2) {
                shard.unit = unit1;
                shard.mapped = shard.mapped + unit1->GetQubitCount();
            }
            if (shard.unit == unit3) {
                shard.unit = unit1;
                shard.mapped = shard.mapped + unit1->GetQubitCount() + unit2->GetQubitCount();
        }
    }

    (unit->*fn)(shards[bit1].mapped, shards[bit2].mapped, shards[bit3].mapped);
}

double SeparatedUnit::Prob(bitLenInt qubit)
{
    QuantumBitShard &shard = shards[qubitIndex];
    return (shard.unit->Prob)(shard.mapped);
}

double SeparatedUnit::ProbAll(bitCapInt perm)
{
    double result = 1.0;

    for (auto shard : shards) {
        p = 0;
        //for (auto bit : shards[i].bits) {
        //    p |= perm & (1 << bit) ? (1 << shards[i].bits[bit]) : 0;
        //}
        // XXX: Reconstruct the perm for this particular CU's mapping.
        // result *= Call<&CoherentUnit::ProbAll>(p) << i;
    }

    return result;
}

void SeparatedUnit::ProbArray(double* probArray)
{
    for (int bit = 0; bit < shards.length(); bit++) {
        probArray[bit] = Prob(bit);
    }
}

/// Measure a bit
bool SeparatedUnit::M(bitLenInt qubitIndex)
{
    QuantumBitShard &shard = shards[qubitIndex];
    bool result = shard.unit->M(shard.mapped);

    /*
     * Decomposes all of the bits in the shard, performing M() on each one and
     * setting each new CU to the appropriate value.
     */
    Decompose(shard);

    return result;
}

/// Measure permutation state of a register
bitCapInt SeparatedUnit::MReg(bitLenInt start, bitLenInt length)
{
    bitLenInt end = start + length;
    bitCapInt result = 0;

    for (bitLenInt bit = start; bit < end; bit++) {
        result |= M(bit) << bit;
    }

    return result;
}

void SeparatedUnit::SetBit(bitLenInt qubitIndex, bool value)
{
    QuantumBitShard &shard = shards[qubitIndex];
    shard.unit->SetBit(shard.mapped, value);
    Decompose(shard);
}

/// Set register bits to given permutation
void SeparatedUnit::SetReg(bitLenInt start, bitLenInt length, bitCapInt value)
{
    for (bitLenInt bit = start; bit < length; bit++) {
        QuantumBitShard &shard = shards[bit];
        shard.unit->SetBit(shard.mapped, value & (1 << bit));
        Decompose(shard);
    }
}

void SeparatedUnit::Swap(bitLenInt qubitIndex1, bitLenInt qubitIndex2)
{
    QuantumBitShard &shard1 = shards[qubitIndex1];
    QuantumBitShard &shard2 = shards[qubitIndex2];

    QuantumBitShard tmp;

    // Swap the bit mapping.
    tmp.mapped = shard1.mapped;
    shard1.mapped = shard2.mapped;
    shard2.mapped = tmp.mapped;

    // Swap the CohesiveUnit object.
    tmp.unit = shard1.unit;
    shard1.unit = shard2.unit;
    shard2.unit = tmp.unit;
}

void SeparatedUnit::Swap(bitLenInt qubitIndex1, bitLenInt qubitIndex2, bitLenInt length)
{
    for (bitLenInt i = 0; i < length; i++) {
        Swap(qubitIndex1 + i, qubitIndex2 + i);
    }
}

void SeparatedUnit::AND(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit)
{
    EntangleAndCall(inputBit1, inputBit2, outputBit, &CohesiveUnit::AND);
}

void SeparatedUnit::OR(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit)
{
    EntangleAndCall(inputBit1, inputBit2, outputBit, &CohesiveUnit::OR);
}

void SeparatedUnit::XOR(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit)
{
    EntangleAndCall(inputBit1, inputBit2, outputBit, &CohesiveUnit::XOR);
}

void SeparatedUnit::CLAND(bitLenInt inputQBit, bool inputClassicalBit, bitLenInt outputQBit)
{
    EntangleAndCall(inputBit1, inputBit2, [&](CohesiveUnit *unit, bitLenInt b1, bitLenInt b2) {
            unit->CLAND(b1, inputClassicalBit, b2);
        });
}

void SeparatedUnit::CLOR(bitLenInt inputQBit, bool inputClassicalBit, bitLenInt outputQBit)
{
    EntangleAndCall(inputBit1, inputBit2, [&](CohesiveUnit *unit, bitLenInt b1, bitLenInt b2) {
            unit->CLOR(b1, inputClassicalBit, b2);
        });
}

void SeparatedUnit::CLXOR(bitLenInt inputQBit, bool inputClassicalBit, bitLenInt outputQBit)
{
    EntangleAndCall(inputBit1, inputBit2, [&](CohesiveUnit *unit, bitLenInt b1, bitLenInt b2) {
            unit->CLXOR(b1, inputClassicalBit, b2);
        });
}

void SeparatedUnit::CCNOT(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit)
{
    EntangleAndCall(inputBit1, inputBit2, outputBit, &CohesiveUnit::CCNOT);
}

void SeparatedUnit::AntiCCNOT(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit)
{
    EntangleAndCall(inputBit1, inputBit2, outputBit, &CohesiveUnit::AntiCCNOT);
}

void SeparatedUnit::H(bitLenInt qubit)
{
    shards[qubit]->unit->H(shards[qubit].mapped);
}

void SeparatedUnit::X(bitLenInt qubit)
{
    shards[qubit]->unit->X(shards[qubit].mapped);
}

void SeparatedUnit::Y(bitLenInt qubit)
{
    shards[qubit]->unit->Y(shards[qubit].mapped);
}

void SeparatedUnit::Z(bitLenInt qubitIndex)
{
    shards[qubit]->unit->Z(shards[qubit].mapped);
}

void SeparatedUnit::X(bitLenInt start, bitLenInt length)
{
    for (bitLenInt i = 0; i < length; i++) {
        X(start + i);
    }
}

void SeparatedUnit::CY(bitLenInt control, bitLenInt target)
{
    EntangleAndCall(control, target, &CohesiveUnit::CY);
}

void SeparatedUnit::CZ(bitLenInt control, bitLenInt target)
{
    EntangleAndCall(control, target, &CohesiveUnit::CZ);
}

void SeparatedUnit::RT(double radians, bitLenInt qubitIndex)
{
    EntangleAndCall(qubit, [&](CohesiveUnit *unit, bitLenInt q) {
            unit->RT(radians, q);
        });
}

void SeparatedUnit::RTDyad(int numerator, int denominator, bitLenInt qubit)
{
    shards[qubit]->unit->RTDyad(numerator, denominator, shards[qubit].mapped);
}

void SeparatedUnit::RX(double radians, bitLenInt qubit)
{
    shards[qubit]->unit->RX(radians, shards[qubit].mapped);
}

void SeparatedUnit::RXDyad(int numerator, int denominator, bitLenInt qubit)
{
    shards[qubit]->unit->RXDyad(numerator, denominator, shards[qubit].mapped);
}

void SeparatedUnit::RY(double radians, bitLenInt qubit)
{
    shards[qubit]->unit->RY(radians, shards[qubit].mapped);
}

void SeparatedUnit::RYDyad(int numerator, int denominator, bitLenInt qubit)
{
    shards[qubit]->unit->RYDyad(numerator, denominator, shards[qubit].mapped);
}

// XXX XXX XXX Didn't make any further changes below here...

void SeparatedUnit::RZ(double radians, bitLenInt qubitIndex)
{
    coherentUnits[qubitLookup[qubitIndex].cu]->RZ(radians, qubitLookup[qubitIndex].qb);
}

void SeparatedUnit::RZDyad(int numerator, int denominator, bitLenInt qubitIndex)
{
    coherentUnits[qubitLookup[qubitIndex].cu]->RZDyad(numerator, denominator, qubitLookup[qubitIndex].qb);
}

void SeparatedUnit::CRT(double radians, bitLenInt control, bitLenInt target)
{
    std::vector<bitLenInt> indices(2);
    indices[0] = control;
    indices[1] = target;
    EntangleIndices(indices);

    coherentUnits[qubitLookup[control].cu]->CRT(radians, qubitLookup[control].qb, qubitLookup[target].qb);
}

void SeparatedUnit::CRTDyad(int numerator, int denominator, bitLenInt control, bitLenInt target)
{
    std::vector<bitLenInt> indices(2);
    indices[0] = control;
    indices[1] = target;
    EntangleIndices(indices);

    coherentUnits[qubitLookup[control].cu]->CRTDyad(
        numerator, denominator, qubitLookup[control].qb, qubitLookup[target].qb);
}

void SeparatedUnit::CRY(double radians, bitLenInt control, bitLenInt target)
{
    std::vector<bitLenInt> indices(2);
    indices[0] = control;
    indices[1] = target;
    EntangleIndices(indices);

    coherentUnits[qubitLookup[control].cu]->CRY(radians, qubitLookup[control].qb, qubitLookup[target].qb);
}

void SeparatedUnit::CRYDyad(int numerator, int denominator, bitLenInt control, bitLenInt target)
{
    std::vector<bitLenInt> indices(2);
    indices[0] = control;
    indices[1] = target;
    EntangleIndices(indices);

    coherentUnits[qubitLookup[control].cu]->CRYDyad(
        numerator, denominator, qubitLookup[control].qb, qubitLookup[target].qb);
}

void SeparatedUnit::CRZ(double radians, bitLenInt control, bitLenInt target)
{
    std::vector<bitLenInt> indices(2);
    indices[0] = control;
    indices[1] = target;
    EntangleIndices(indices);

    coherentUnits[qubitLookup[control].cu]->CRZ(radians, qubitLookup[control].qb, qubitLookup[target].qb);
}

void SeparatedUnit::CRZDyad(int numerator, int denominator, bitLenInt control, bitLenInt target)
{
    std::vector<bitLenInt> indices(2);
    indices[0] = control;
    indices[1] = target;
    EntangleIndices(indices);

    coherentUnits[qubitLookup[control].cu]->CRZDyad(
        numerator, denominator, qubitLookup[control].qb, qubitLookup[target].qb);
}

/// "Circular shift right" - shift bits right, and carry first bits.
void SeparatedUnit::ROL(bitLenInt shift, bitLenInt start, bitLenInt length)
{
    if ((length > 0) && (shift > 0)) {
        bitLenInt end = start + length;
        if (shift >= length) {
            SetReg(start, length, 0);
        } else {
            Reverse(start, end);
            Reverse(start, start + shift);
            Reverse(start + shift, end);
        }
    }
}

/// "Circular shift right" - shift bits right, and carry first bits.
void SeparatedUnit::ROR(bitLenInt shift, bitLenInt start, bitLenInt length)
{
    if ((length > 0) && (shift > 0)) {
        bitLenInt end = start + length;
        if (shift >= length) {
            SetReg(start, length, 0);
        } else {
            Reverse(start + shift, end);
            Reverse(start, start + shift);
            Reverse(start, end);
        }
    }
}

void SeparatedUnit::INC(bitCapInt toAdd, bitLenInt start, bitLenInt length)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->INC(toAdd, qubitLookup[start].qb, length);
}

void SeparatedUnit::INCC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);
    QbListEntry carryQbe;
    carryQbe.cu = qubitLookup[carryIndex].cu;
    carryQbe.start = qubitLookup[carryIndex].qb;
    carryQbe.length = 1;
    qbList.push_back(carryQbe);
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->INCC(toAdd, qubitLookup[start].qb, length, qubitLookup[carryIndex].qb);
}

void SeparatedUnit::INCS(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt overflowIndex)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);
    QbListEntry overflowQbe;
    overflowQbe.cu = qubitLookup[overflowIndex].cu;
    overflowQbe.start = qubitLookup[overflowIndex].qb;
    overflowQbe.length = 1;
    qbList.push_back(overflowQbe);
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->INCS(toAdd, qubitLookup[start].qb, length, qubitLookup[overflowIndex].qb);
}

void SeparatedUnit::INCSC(
    bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt overflowIndex, bitLenInt carryIndex)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);
    std::vector<QbListEntry> qbExtra(2);
    QbListEntry extraQbe;
    extraQbe.cu = qubitLookup[overflowIndex].cu;
    extraQbe.start = qubitLookup[overflowIndex].qb;
    extraQbe.length = 1;
    qbExtra[0] = extraQbe;
    extraQbe.cu = qubitLookup[carryIndex].cu;
    extraQbe.start = qubitLookup[carryIndex].qb;
    extraQbe.length = 1;
    qbExtra[1] = extraQbe;
    qbList.insert(qbList.end(), qbExtra.begin(), qbExtra.end());
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->INCSC(
        toAdd, qubitLookup[start].qb, length, qubitLookup[overflowIndex].qb, qubitLookup[carryIndex].qb);
}

void SeparatedUnit::INCSC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);
    QbListEntry carryQbe;
    carryQbe.cu = qubitLookup[carryIndex].cu;
    carryQbe.start = qubitLookup[carryIndex].qb;
    carryQbe.length = 1;
    qbList.push_back(carryQbe);
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->INCSC(toAdd, qubitLookup[start].qb, length, qubitLookup[carryIndex].qb);
}

void SeparatedUnit::INCBCD(bitCapInt toAdd, bitLenInt start, bitLenInt length)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->INCBCD(toAdd, qubitLookup[start].qb, length);
}

void SeparatedUnit::INCBCDC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);
    QbListEntry carryQbe;
    carryQbe.cu = qubitLookup[carryIndex].cu;
    carryQbe.start = qubitLookup[carryIndex].qb;
    carryQbe.length = 1;
    qbList.push_back(carryQbe);
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->INCBCDC(toAdd, qubitLookup[start].qb, length, qubitLookup[carryIndex].qb);
}

void SeparatedUnit::DEC(bitCapInt toAdd, bitLenInt start, bitLenInt length)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->DEC(toAdd, qubitLookup[start].qb, length);
}

void SeparatedUnit::DECC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);
    QbListEntry carryQbe;
    carryQbe.cu = qubitLookup[carryIndex].cu;
    carryQbe.start = qubitLookup[carryIndex].qb;
    carryQbe.length = 1;
    qbList.push_back(carryQbe);
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->DECC(toAdd, qubitLookup[start].qb, length, qubitLookup[carryIndex].qb);
}

void SeparatedUnit::DECS(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt overflowIndex)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);
    QbListEntry overflowQbe;
    overflowQbe.cu = qubitLookup[overflowIndex].cu;
    overflowQbe.start = qubitLookup[overflowIndex].qb;
    overflowQbe.length = 1;
    qbList.push_back(overflowQbe);
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->DECS(toAdd, qubitLookup[start].qb, length, qubitLookup[overflowIndex].qb);
}

void SeparatedUnit::DECSC(
    bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt overflowIndex, bitLenInt carryIndex)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);
    std::vector<QbListEntry> qbExtra(2);
    QbListEntry extraQbe;
    extraQbe.cu = qubitLookup[overflowIndex].cu;
    extraQbe.start = qubitLookup[overflowIndex].qb;
    extraQbe.length = 1;
    qbExtra[0] = extraQbe;
    extraQbe.cu = qubitLookup[carryIndex].cu;
    extraQbe.start = qubitLookup[carryIndex].qb;
    extraQbe.length = 1;
    qbExtra[1] = extraQbe;
    qbList.insert(qbList.end(), qbExtra.begin(), qbExtra.end());
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->DECSC(
        toAdd, qubitLookup[start].qb, length, qubitLookup[overflowIndex].qb, qubitLookup[carryIndex].qb);
}

void SeparatedUnit::DECSC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);
    QbListEntry carryQbe;
    carryQbe.cu = qubitLookup[carryIndex].cu;
    carryQbe.start = qubitLookup[carryIndex].qb;
    carryQbe.length = 1;
    qbList.push_back(carryQbe);
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->DECSC(toAdd, qubitLookup[start].qb, length, qubitLookup[carryIndex].qb);
}

void SeparatedUnit::DECBCD(bitCapInt toAdd, bitLenInt start, bitLenInt length)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->DECBCD(toAdd, qubitLookup[start].qb, length);
}

void SeparatedUnit::DECBCDC(bitCapInt toAdd, bitLenInt start, bitLenInt length, bitLenInt carryIndex)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);
    QbListEntry carryQbe;
    carryQbe.cu = qubitLookup[carryIndex].cu;
    carryQbe.start = qubitLookup[carryIndex].qb;
    carryQbe.length = 1;
    qbList.push_back(carryQbe);
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->DECBCDC(toAdd, qubitLookup[start].qb, length, qubitLookup[carryIndex].qb);
}

void SeparatedUnit::QFT(bitLenInt start, bitLenInt length)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->QFT(qubitLookup[start].qb, length);
}

void SeparatedUnit::ZeroPhaseFlip(bitLenInt start, bitLenInt length)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->ZeroPhaseFlip(qubitLookup[start].qb, length);
}

void SeparatedUnit::CPhaseFlipIfLess(bitCapInt greaterPerm, bitLenInt start, bitLenInt length, bitLenInt flagIndex)
{
    std::vector<QbListEntry> qbList(length);
    GetParallelBitList(start, length, qbList);
    QbListEntry flagQbe;
    flagQbe.cu = qubitLookup[flagIndex].cu;
    flagQbe.start = qubitLookup[flagIndex].qb;
    flagQbe.length = 1;
    qbList.push_back(flagQbe);
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    coherentUnits[qubitLookup[start].cu]->CPhaseFlipIfLess(
        greaterPerm, qubitLookup[start].qb, length, qubitLookup[flagIndex].qb);
}

void SeparatedUnit::PhaseFlip()
{
    for (bitLenInt i = 0; i < coherentUnits.size(); i++) {
        coherentUnits[i]->PhaseFlip();
    }
}

/**
 * Set 8 bit register bits by a superposed index-offset-based read from
 * classical memory
 *
 * "inputStart" is the start index of 8 qubits that act as an index into
 * the 256 byte "values" array. The "outputStart" bits are first cleared,
 * then the separable |input, 00000000> permutation state is mapped to
 * |input, values[input]>, with "values[input]" placed in the "outputStart"
 * register.
 *
 * While a CoherentUnit represents an interacting set of qubit-based
 * registers, or a virtual quantum chip, the registers need to interact in
 * some way with (classical or quantum) RAM. SuperposeReg8 is a RAM access
 * method similar to the X addressing mode of the MOS 6502 chip, if the X
 * register can be in a state of coherent superposition when it loads from
 * RAM.
 *
 * The physical motivation for this addressing mode can be explained as
 * follows: say that we have a superconducting quantum interface device
 * (SQUID) based chip. SQUIDs have already been demonstrated passing
 * coherently superposed electrical currents. In a sufficiently
 * quantum-mechanically isolated qubit chip with a classical cache, with
 * both classical RAM and registers likely cryogenically isolated from the
 * environment, SQUIDs could (hopefully) pass coherently superposed
 * electrical currents into the classical RAM cache to load values into a
 * qubit register. The state loaded would be a superposition of the values
 * of all RAM to which coherently superposed electrical currents were
 * passed.
 *
 * In qubit system similar to the MOS 6502, say we have qubit-based
 * "accumulator" and "X index" registers, and say that we start with a
 * superposed X index register. In (classical) X addressing mode, the X
 * index register value acts an offset into RAM from a specified starting
 * address. The X addressing mode of a LoaD Accumulator (LDA) instruction,
 * by the physical mechanism described above, should load the accumulator
 * in quantum parallel with the values of every different address of RAM
 * pointed to in superposition by the X index register. The superposed
 * values in the accumulator are entangled with those in the X index
 * register, by way of whatever values the classical RAM pointed to by X
 * held at the time of the load. (If the RAM at index "36" held an unsigned
 * char value of "27," then the value "36" in the X index register becomes
 * entangled with the value "27" in the accumulator, and so on in quantum
 * parallel for all superposed values of the X index register, at once.) If
 * the X index register or accumulator are then measured, the two registers
 * will both always collapse into a random but valid key-value pair of X
 * index offset and value at that classical RAM address.
 *
 * Note that a "superposed store operation in classical RAM" is not
 * possible by analagous reasoning. Classical RAM would become entangled
 * with both the accumulator and the X register. When the state of the
 * registers was collapsed, we would find that only one "store" operation
 * to a single memory address had actually been carried out, consistent
 * with the address offset in the collapsed X register and the byte value
 * in the collapsed accumulator. It would not be possible by this model to
 * write in quantum parallel to more than one address of classical memory
 * at a time.
 */

unsigned char SeparatedUnit::SuperposeReg8(bitLenInt inputStart, bitLenInt outputStart, unsigned char* values)
{
    std::vector<QbListEntry> qbListInput(8);
    GetParallelBitList(inputStart, 8, qbListInput);
    std::vector<QbListEntry> qbListOutput(8);
    GetParallelBitList(outputStart, 8, qbListOutput);
    std::vector<QbListEntry> qbList(qbListInput.size() + qbListOutput.size());
    std::copy(qbListInput.begin(), qbListInput.end(), qbList.begin());
    std::copy(qbListOutput.begin(), qbListOutput.end(), qbList.begin() + qbListInput.size());
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    return coherentUnits[qubitLookup[inputStart].cu]->SuperposeReg8(
        qubitLookup[inputStart].qb, qubitLookup[outputStart].qb, values);
}

/**
 * Add to entangled 8 bit register state with a superposed
 * index-offset-based read from classical memory
 *
 * inputStart" is the start index of 8 qubits that act as an index into the
 * 256 byte "values" array. The "outputStart" bits would usually already be
 * entangled with the "inputStart" bits via a SuperposeReg8() operation.
 * With the "inputStart" bits being a "key" and the "outputStart" bits
 * being a value, the permutation state |key, value> is mapped to |key,
 * value + values[key]>. This is similar to classical parallel addition of
 * two arrays.  However, when either of the registers are measured, both
 * registers will collapse into one random VALID key-value pair, with any
 * addition or subtraction done to the "value." See SuperposeReg8() for
 * context.
 *
 * While a CoherentUnit represents an interacting set of qubit-based
 * registers, or a virtual quantum chip, the registers need to interact in
 * some way with (classical or quantum) RAM. SuperposeReg8 is a RAM access
 * method similar to the X addressing mode of the MOS 6502 chip, if the X
 * register can be in a state of coherent superposition when it loads from
 * RAM. "AdcSuperposReg8" and "SbcSuperposeReg8" perform add and subtract
 * (with carry) operations on a state usually initially prepared with
 * SuperposeReg8().
 */
unsigned char SeparatedUnit::AdcSuperposeReg8(
    bitLenInt inputStart, bitLenInt outputStart, bitLenInt carryIndex, unsigned char* values)
{
    QbListEntry carryQbe;
    std::vector<QbListEntry> qbListInput(8);
    GetParallelBitList(inputStart, 8, qbListInput);
    std::vector<QbListEntry> qbListOutput(8);
    GetParallelBitList(outputStart, 8, qbListOutput);
    std::vector<QbListEntry> qbList(qbListInput.size() + qbListOutput.size() + 1);
    std::copy(qbListInput.begin(), qbListInput.end(), qbList.begin());
    std::copy(qbListOutput.begin(), qbListOutput.end(), qbList.begin() + qbListInput.size());
    carryQbe.cu = qubitLookup[carryIndex].cu;
    carryQbe.start = qubitLookup[carryIndex].qb;
    carryQbe.length = 1;
    qbList[qbList.size() - 1] = carryQbe;
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    return coherentUnits[qubitLookup[inputStart].cu]->AdcSuperposeReg8(
        qubitLookup[inputStart].qb, qubitLookup[outputStart].qb, qubitLookup[carryIndex].qb, values);
}

/**
 * Subtract from an entangled 8 bit register state with a superposed
 * index-offset-based read from classical memory
 *
 * "inputStart" is the start index of 8 qubits that act as an index into
 * the 256 byte "values" array. The "outputStart" bits would usually
 * already be entangled with the "inputStart" bits via a SuperposeReg8()
 * operation.  With the "inputStart" bits being a "key" and the
 * "outputStart" bits being a value, the permutation state |key, value> is
 * mapped to |key, value - values[key]>. This is similar to classical
 * parallel addition of two arrays.  However, when either of the registers
 * are measured, both registers will collapse into one random VALID
 * key-value pair, with any addition or subtraction done to the "value."
 * See CoherentUnit::SuperposeReg8 for context.
 *
 * While a CoherentUnit represents an interacting set of qubit-based
 * registers, or a virtual quantum chip, the registers need to interact in
 * some way with (classical or quantum) RAM. SuperposeReg8 is a RAM access
 * method similar to the X addressing mode of the MOS 6502 chip, if the X
 * register can be in a state of coherent superposition when it loads from
 * RAM. "AdcSuperposReg8" and "SbcSuperposeReg8" perform add and subtract
 * (with carry) operations on a state usually initially prepared with
 * SuperposeReg8().
 */
unsigned char SeparatedUnit::SbcSuperposeReg8(
    bitLenInt inputStart, bitLenInt outputStart, bitLenInt carryIndex, unsigned char* values)
{
    QbListEntry carryQbe;
    std::vector<QbListEntry> qbListInput(8);
    GetParallelBitList(inputStart, 8, qbListInput);
    std::vector<QbListEntry> qbListOutput(8);
    GetParallelBitList(outputStart, 8, qbListOutput);
    std::vector<QbListEntry> qbList(qbListInput.size() + qbListOutput.size() + 1);
    std::copy(qbListInput.begin(), qbListInput.end(), qbList.begin());
    std::copy(qbListOutput.begin(), qbListOutput.end(), qbList.begin() + qbListInput.size());
    carryQbe.cu = qubitLookup[carryIndex].cu;
    carryQbe.start = qubitLookup[carryIndex].qb;
    carryQbe.length = 1;
    qbList[qbList.size() - 1] = carryQbe;
    OptimizeParallelBitList(qbList);

    EntangleBitList(qbList);

    return coherentUnits[qubitLookup[inputStart].cu]->SbcSuperposeReg8(
        qubitLookup[inputStart].qb, qubitLookup[outputStart].qb, qubitLookup[carryIndex].qb, values);
}

/**
 * Compile an order-preserving list of CoherentUnit bit strings for applying an register-wise operation
 *
 * This operation optimizes compiling a list out of qubit pile when bit order is important. We apply register-wise
 * operations over a pile of arbitrarily entangled and separated qubits. Entangled qubits are stored together in single
 * CoherentUnit objects, but their mapping to SeparatedUnit bit indices can be generally random. Sometimes, we must
 * preserve bit order to correctly carry out the operation, whereas sometimes our operation is bitwise parallel and does
 * not depend on the ordering of bits in the list.
 */
void SeparatedUnit::GetOrderedBitList(bitLenInt start, bitLenInt length, std::vector<QbListEntry>& qbList)
{
    // Start by getting a list (of sublists) of all the bits we need, with bit sublist length of 1.
    bitLenInt i, j;
    QbLookup qbl;
    QbListEntry qbe;
    for (i = 0; i < length; i++) {
        qbl = qubitLookup[start + i];
        qbe.cu = qbl.cu;
        qbe.start = qbl.qb;
        qbe.length = 1;
        qbList[i] = qbe;
    }

    // If contiguous sublists in the list we just made are also contiguous in the same coherent unit, we can combine
    // them to optimize with register-wise gate methods.
    j = 0;
    for (i = 0; i < (length - 1); i++) {
        if ((qbList[j].cu == qbList[j + 1].cu) && ((qbList[j].start + qbList[j].length) == qbList[j + 1].start)) {
            qbList[j].length += qbList[j + 1].length;
            qbList.erase(qbList.begin() + j + 1);
        } else {
            j++;
        }
    }
}

/// Compile a list of CoherentUnit bit strings for applying a bitwise-parallel operation
/**
 * This operation optimizes compiling a list out of qubit pile when bit order is not important. We apply register-wise
 * operations over a pile of arbitrarily entangled and separated qubits. Entangled qubits are stored together in single
 * CoherentUnit objects, but their mapping to SeparatedUnit bit indices can be generally random. Sometimes, we must
 * preserve bit order to correctly carry out the operation, whereas sometimes our operation is bitwise parallel and does
 * not depend on the ordering of bits in the list.
 */
void SeparatedUnit::GetParallelBitList(bitLenInt start, bitLenInt length, std::vector<QbListEntry>& qbList)
{
    // Start by getting a list (of sublists) of all the bits we need, with bit sublist length of 1.
    bitLenInt i, j;
    QbLookup qbl;
    QbListEntry qbe;
    for (i = 0; i < length; i++) {
        qbl = qubitLookup[start + i];
        qbe.cu = qbl.cu;
        qbe.start = qbl.qb;
        qbe.length = 1;
        qbList[i] = qbe;
    }
    // The ordering of bits returned is unimportant, so we can better optimize by sorting this list by CoherentUnit
    // index and qubit index, to maximize the reduction of the list.
    std::sort(qbList.begin(), qbList.end(), compare);
    // If contiguous sublists in the list we just sorted are also contiguous in the same coherent unit, we can combine
    // them to optimize with register-wise gate methods.
    j = 0;
    for (i = 0; i < (length - 1); i++) {
        if ((qbList[j].cu == qbList[j + 1].cu) && ((qbList[j].start + qbList[j].length) == qbList[j + 1].start)) {
            qbList[j].length += qbList[j + 1].length;
            qbList.erase(qbList.begin() + j + 1);
        } else {
            j++;
        }
    }
}

/// Combines two lists returned by GetParallelBitList() by the same logic as that algorithm
void SeparatedUnit::OptimizeParallelBitList(std::vector<QbListEntry>& qbList)
{
    if (qbList.size() < 2) {
        return;
    }

    bitLenInt i, j;
    bitLenInt length = qbList.size();
    // The ordering of bits returned is unimportant, so we can better optimize by sorting this list by CoherentUnit
    // index and qubit index, to maximize the reduction of the list.
    std::sort(qbList.begin(), qbList.end(), compare);
    // If contiguous sublists in the list we just sorted are also contiguous in the same coherent unit, we can combine
    // them to optimize with register-wise gate methods.
    j = 0;
    for (i = 0; i < (length - 1); i++) {
        if ((qbList[j].cu == qbList[j + 1].cu) && ((qbList[j].start + qbList[j].length) == qbList[j + 1].start)) {
            qbList[j].length += qbList[j + 1].length;
            qbList.erase(qbList.begin() + j + 1);
        } else {
            j++;
        }
    }
}

/// Entangle and sort the indices of a list of CoherentUnit objects
void SeparatedUnit::EntangleBitList(std::vector<QbListEntry> qbList)
{
    if (qbList.size() < 2) {
        return;
    }

    bitLenInt i, j, k;
    bitLenInt firstCu, cuLen, invLookup, cuRemoved;
    QbListEntry qbe;

    firstCu = qbList[0].cu;
    k = coherentUnits[firstCu]->GetQubitCount();
    for (i = 1; i < qbList.size(); i++) {
        qbe = qbList[i];
        cuLen = coherentUnits[qbe.cu]->GetQubitCount();
        for (j = 0; j < cuLen; j++) {
            invLookup = qubitInverseLookup[qbe.cu * qubitCount + j];
            qubitLookup[invLookup].cu = firstCu;
            qubitLookup[invLookup].qb = k + j;
            qubitInverseLookup[firstCu * qubitCount + k + j] = invLookup;
        }
        coherentUnits[firstCu]->Cohere(*(coherentUnits[qbe.cu]));
        k += cuLen;
    }

    // Swap qubits into appropriate order, then update coherentUnits list.
    cuLen = coherentUnits[firstCu]->GetQubitCount();
    QuickSortQubits(&(qubitInverseLookup[firstCu * qubitCount]), 0, cuLen - 1, coherentUnits[firstCu]);
    // Update lookup table
    for (i = 0; i < cuLen; i++) {
        invLookup = qubitInverseLookup[firstCu * qubitCount + i];
        qubitLookup[invLookup].cu = firstCu;
        qubitLookup[invLookup].qb = i;
    }

    // Update coherentUnit list and inverse lookup at end
    cuLen = qbList.size() - 1;
    if (cuLen > 0) {
        std::vector<bitLenInt> cuToDelete(cuLen);
        for (i = 0; i < cuLen; i++) {
            cuToDelete[i] = qbList[i + 1].cu;
        }
        std::sort(cuToDelete.begin(), cuToDelete.end());
        for (i = 0; i < cuLen; i++) {
            cuRemoved = cuToDelete[cuLen - i - 1];
            coherentUnits.erase(coherentUnits.begin() + cuRemoved);
            for (j = 0; j < qubitCount; j++) {
                if (qubitLookup[j].cu >= cuRemoved) {
                    qubitLookup[j].cu--;
                }
            }
            for (j = cuRemoved; j < coherentUnits.size(); j++) {
                std::copy(&(qubitInverseLookup[0]) + (j + 1) * qubitCount,
                    &(qubitInverseLookup[0]) + (j + 2) * qubitCount, &(qubitInverseLookup[0]) + j * qubitCount);
            }
        }
    }
}

void SeparatedUnit::EntangleIndices(std::vector<bitLenInt> indices)
{
    QbListEntry qbe;
    std::vector<QbListEntry> qbList(indices.size());
    for (bitLenInt i = 0; i < indices.size(); i++) {
        qbe.cu = qubitLookup[indices[i]].cu;
        qbe.start = qubitLookup[indices[i]].qb;
        qbe.length = 1;
        qbList[i] = qbe;
    }
    OptimizeParallelBitList(qbList);
    EntangleBitList(qbList);
}

void SeparatedUnit::QuickSortQubits(bitLenInt* arr, bitLenInt low, bitLenInt high, std::weak_ptr<CoherentUnit> cuWeak)
{
    std::shared_ptr<CoherentUnit> cu = cuWeak.lock();
    int i = low, j = high;
    int pivot = arr[(low + high) / 2];

    while (i <= j) {
        while (arr[i] < pivot) {
            i++;
        }
        while (arr[j] > pivot) {
            j--;
        }
        if (i <= j) {
            std::swap(arr[i], arr[j]);
            cu->Swap(i, j);
            i++;
            j--;
        }
    }
    if (low < j) {
        QuickSortQubits(arr, low, j, cuWeak);
    }
    if (i < high) {
        QuickSortQubits(arr, i, high, cuWeak);
    }
}

void SeparatedUnit::DecohereOrDispose(bool isDecohere, bitLenInt start, bitLenInt length, CoherentUnit* destination)
{
    bitLenInt i, j, k;
    std::vector<QbListEntry> qbList(length);
    GetOrderedBitList(start, length, qbList);
    EntangleBitList(qbList);

    bitLenInt cu = qubitLookup[start].cu;
    bitLenInt cuStart = qubitLookup[start].qb;
    bitLenInt cuLen = coherentUnits[cu]->GetQubitCount();
    if (cuLen == length) {
        if (isDecohere) {
            std::unique_ptr<Complex16[]> sv(new Complex16[1 << cuLen]);
            coherentUnits[cu]->CloneRawState(&(sv[0]));
            destination->SetQuantumState(&(sv[0]));
        }
        coherentUnits.erase(coherentUnits.begin() + cu);

        for (i = cu; i < (qubitCount - 1); i++) {
            std::copy(&(qubitInverseLookup[0]) + (i + 1) * qubitCount, &(qubitInverseLookup[0]) + (i + 2) * qubitCount,
                &(qubitInverseLookup[0]) + i * qubitCount);
        }
        k = 0;
        for (i = 0; i < qubitCount; i++) {
            if (qubitLookup[k].cu == cu) {
                for (j = k; j < (qubitCount - 1); j++) {
                    qubitLookup[j] = qubitLookup[j + 1];
                }
            } else {
                if (qubitLookup[k].cu > cu) {
                    qubitLookup[k].cu--;
                }
                k++;
            }
        }
    } else {
        if (isDecohere) {
            coherentUnits[cu]->Decohere(qubitLookup[start].qb, length, *destination);
        } else {
            coherentUnits[cu]->Dispose(qubitLookup[start].qb, length);
        }

        k = 0;
        for (i = 0; i < qubitCount; i++) {
            if (qubitLookup[k].cu == cu) {
                if ((qubitLookup[k].qb >= cuStart) && (qubitLookup[k].qb < (cuStart + length))) {
                    for (j = k; j < (qubitCount - 1); j++) {
                        qubitLookup[j] = qubitLookup[j + 1];
                    }
                } else {
                    if (qubitLookup[k].qb > cuStart) {
                        qubitLookup[k].qb -= length;
                    }
                    k++;
                }
            }
        }
    }

    qubitCount -= length;
    maxQPower = 1 << qubitCount;

    std::unique_ptr<QbLookup[]> ql(new QbLookup[qubitCount]);
    std::copy(&(qubitLookup[0]), &(qubitLookup[0]) + qubitCount, &(ql[0]));
    qubitLookup = std::move(ql);

    std::unique_ptr<bitLenInt[]> qil(new bitLenInt[qubitCount * qubitCount]());
    for (i = 0; i < coherentUnits.size(); i++) {
        std::copy(&(qubitInverseLookup[i * (qubitCount + length)]),
            &(qubitInverseLookup[i * (qubitCount + length)]) + qubitCount, &(qil[i * qubitCount]));
    }
    qubitInverseLookup = std::move(qil);
}

} // namespace Qrack
