# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 16:04:40 2020

@author: 41207
"""
import numpy
from rdkit import Chem
import copy
import numpy.linalg
import math
from rdkit.Chem import rdchem
from rdkit.Chem.EState import Fingerprinter as ESFP
from rdkit.Chem import Lipinski as LPK
from PyBioMed.PyMolecule.AtomProperty import GetRelativeAtomicProperty
from rdkit.Chem import MolSurf as MOE
from rdkit.Chem.EState import EState_VSA as EVSA
from rdkit.Chem import rdPartialCharges as GMCharge
from rdkit.Chem import Crippen
from rdkit.Chem import MolSurf as MS
from sklearn.externals import joblib
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import GraphDescriptors as GD
import pandas as pd
import csv
try:
    # TODO: Is this deprecated?
    # https://github.com/rdkit/rdkit/issues/2741#issuecomment-546709239
    from rdkit.Chem import pyPeriodicTable as PeriodicTable
except ImportError:
    from rdkit.Chem import PeriodicTable

def Model(current_path,model): 
    attribute= joblib.load(current_path+model)
    return attribute

def Model_selection(current_path): 
    model1 = joblib.load(current_path+'CYP_inhibitor_1A2.pkl')
    model2 = joblib.load(current_path+'CYP_inhibitor_2C9.pkl')
    model3 = joblib.load(current_path+'CYP_inhibitor_3A4.pkl')
    model4 = joblib.load(current_path+'PGP_inhibitor.pkl')
    model5 = joblib.load(current_path+'CYP_substrate_1A2.pkl')
    model6 = joblib.load(current_path+'CYP_substrate_2C9.pkl')
    model7 = joblib.load(current_path+'CYP_substrate_2C19.pkl')
    model8 = joblib.load(current_path+'CYP_substrate_3A4.pkl')
    model9 = joblib.load(current_path+'CYP_substrate_2D6.pkl')
    model10 = joblib.load(current_path+'CYP_inhibitor_2D6.pkl')
    model11 = joblib.load(current_path+'BBB.pkl')
    model12 = joblib.load(current_path+'CYP_inhibitor_2C19.pkl')
    model13 = joblib.load(current_path+'AMES_Model.pkl')
    model14 = joblib.load(current_path+'F20.pkl')
    model15 = joblib.load(current_path+'HIA.pkl')
    model16 = joblib.load(current_path+'SkinSen.pkl')
    model17 = joblib.load(current_path+'F_30.pkl')
    model18 = joblib.load(current_path+'caco2_Model.pkl')
    model19 = joblib.load(current_path+'logS_Model.pkl')
    model20 = joblib.load(current_path+'logD_Model.pkl')
    model21 = joblib.load(current_path+'PPB_Model.pkl')
    model22 = joblib.load(current_path+'VD_Model.pkl')
    model23 = joblib.load(current_path+'CL_Model.pkl')
    model24 = joblib.load(current_path+'T_Model.pkl')
    model25 = joblib.load(current_path+'hERG_Model.pkl')
    model26 = joblib.load(current_path+'HHT_Model.pkl')
    
    return model1,model2,model3,model4,model5,model6,model7,model8,model9,model10,model11,model12,model13,model14,model15,model16,model17,model18,model19,model20,model21,model22,model23,model24,model25,model26
    
#Generating function of ECFP molecular fingerprint
def ECFP_generate(filename,radius,dimension):
    df = pd.read_csv(filename,engine='python').iloc[:, 0]  
#Import the column of Smiles formula
    smiles_list = numpy.array(df)
    ECFP=[]
    for i in range(len(smiles_list)):
        Smiles=smiles_list[i]
        m= Chem.MolFromSmiles(Smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=dimension)
        fp = list(map(int,list(fp.ToBitString())))
        ECFP.append(fp)
    f = open('Molecular fingerprint.csv','w')
    csv_writer = csv.writer(f)
    b=[]
    for i in range(dimension):
        b.append(i) 
    csv_writer.writerow(b)
    for i in range(len(ECFP)):
        csv_writer.writerow(ECFP[i])
    f.close()
    fingerprint_content = pd.read_csv('Molecular fingerprint.csv')
    des_list = numpy.array(fingerprint_content)
    return des_list
#Generating function of MACCS molecular fingerprint
def MACCS_generate(filename):
    df = pd.read_csv(filename,engine='python').iloc[:, 0]  
    smiles_list = numpy.array(df)
    MACCS=[]
    for i in range(len(smiles_list)):
        Smiles=smiles_list[i]
        m= Chem.MolFromSmiles(Smiles)
        fp = AllChem.GetMACCSKeysFingerprint(m)
        fp = list(map(int,list(fp.ToBitString())))
        MACCS.append(fp)
    f = open('Molecular fingerprint.csv','w')
    csv_writer = csv.writer(f)
    b=[]
    for i in range(0,167):
        b.append(i) 
    csv_writer.writerow(b)
    for i in range(len(MACCS)):
        csv_writer.writerow(MACCS[i])
    f.close()
    fingerprint_content = pd.read_csv('Molecular fingerprint.csv')
    des_list = numpy.array(fingerprint_content)
    return des_list

#Line 90-702 is Generating function of Special 2D molecular fingerprint
periodicTable = rdchem.GetPeriodicTable()
Version = 1.0
iter_step = 12

def _CalculateElementNumber(mol, AtomicNumber=6):

    i = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == AtomicNumber:
            i = i + 1

    return i
def CalculateCarbonNumber(mol):
    return _CalculateElementNumber(mol, AtomicNumber=6)

def _CalculateEntropy(Probability):

    res = 0.0
    for i in Probability:
        if i != 0:
            res = res - i * numpy.log2(i)

    return res
def CalculateBasakIC0(mol):
    BasakIC = 0.0
    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    IC = []
    for i in range(nAtoms):
        at = Hmol.GetAtomWithIdx(i)
        IC.append(at.GetAtomicNum())
    Unique = numpy.unique(IC)
    NAtomType = len(Unique)
    NTAtomType = numpy.zeros(NAtomType, numpy.float)
    for i in range(NAtomType):
        NTAtomType[i] = IC.count(Unique[i])

    if nAtoms != 0:
        # print sum(NTAtomType/nAtoms)
        BasakIC = _CalculateEntropy(NTAtomType / nAtoms)
    else:
        BasakIC = 0.0

    return BasakIC

def _GetBurdenMatrix(mol, propertylabel="m"):
  
    mol = Chem.AddHs(mol)
    Natom = mol.GetNumAtoms()

    AdMatrix = Chem.GetAdjacencyMatrix(mol)
    bondindex = numpy.argwhere(AdMatrix)
    AdMatrix1 = numpy.array(AdMatrix, dtype=numpy.float32)

    # The diagonal elements of B, Bii, are either given by
    # the carbon normalized atomic mass,
    # van der Waals volume, Sanderson electronegativity,
    # and polarizability of atom i.

    for i in range(Natom):
        atom = mol.GetAtomWithIdx(i)
        temp = GetRelativeAtomicProperty(
            element=atom.GetSymbol(), propertyname=propertylabel
        )
        AdMatrix1[i, i] = round(temp, 3)

    # The element of B connecting atoms i and j, Bij,
    # is equal to the square root of the bond
    # order between atoms i and j.

    for i in bondindex:
        bond = mol.GetBondBetweenAtoms(int(i[0]), int(i[1]))
        if bond.GetBondType().name == "SINGLE":
            AdMatrix1[i[0], i[1]] = round(numpy.sqrt(1), 3)
        if bond.GetBondType().name == "DOUBLE":
            AdMatrix1[i[0], i[1]] = round(numpy.sqrt(2), 3)
        if bond.GetBondType().name == "TRIPLE":
            AdMatrix1[i[0], i[1]] = round(numpy.sqrt(3), 3)
        if bond.GetBondType().name == "AROMATIC":
            AdMatrix1[i[0], i[1]] = round(numpy.sqrt(1.5), 3)

    ##All other elements of B (corresponding non bonded
    # atom pairs) are set to 0.001
    bondnonindex = numpy.argwhere(AdMatrix == 0)

    for i in bondnonindex:
        if i[0] != i[1]:

            AdMatrix1[i[0], i[1]] = 0.001

    return numpy.real(numpy.linalg.eigvals(AdMatrix1))

def CalculateBurdenPolarizability(mol):
    temp = _GetBurdenMatrix(mol, propertylabel="alapha")
    temp1 = numpy.sort(temp[temp >= 0])
    temp2 = numpy.sort(numpy.abs(temp[temp < 0]))

    if len(temp1) < 8:
        temp1 = numpy.concatenate((numpy.zeros(8), temp1))
    if len(temp2) < 8:
        temp2 = numpy.concatenate((numpy.zeros(8), temp2))

    bcut = [
        "bcutp16",
        "bcutp15",
        "bcutp14",
        "bcutp13",
        "bcutp12",
        "bcutp11",
        "bcutp10",
        "bcutp9",
        "bcutp8",
        "bcutp7",
        "bcutp6",
        "bcutp5",
        "bcutp4",
        "bcutp3",
        "bcutp2",
        "bcutp1",
    ]
    bcutvalue = numpy.concatenate((temp2[-8:], temp1[-8:]))

    bcutvalue = [round(i, 3) for i in bcutvalue]
    res = dict(zip(bcut, bcutvalue))
    return res

def CalculateBurdenVDW(mol):
    temp = _GetBurdenMatrix(mol, propertylabel="V")
    temp1 = numpy.sort(temp[temp >= 0])
    temp2 = numpy.sort(numpy.abs(temp[temp < 0]))

    if len(temp1) < 8:
        temp1 = numpy.concatenate((numpy.zeros(8), temp1))
    if len(temp2) < 8:
        temp2 = numpy.concatenate((numpy.zeros(8), temp2))

    bcut = [
        "bcutv16",
        "bcutv15",
        "bcutv14",
        "bcutv13",
        "bcutv12",
        "bcutv11",
        "bcutv10",
        "bcutv9",
        "bcutv8",
        "bcutv7",
        "bcutv6",
        "bcutv5",
        "bcutv4",
        "bcutv3",
        "bcutv2",
        "bcutv1",
    ]
    bcutvalue = numpy.concatenate((temp2[-8:], temp1[-8:]))

    bcutvalue = [round(i, 3) for i in bcutvalue]
    res = dict(zip(bcut, bcutvalue))
    return res
def CalculateGutmanTopo(mol):
    nAT = mol.GetNumAtoms()
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    Distance = Chem.GetDistanceMatrix(mol)
    res = 0.0
    for i in range(nAT):
        for j in range(i + 1, nAT):
            res = res + deltas[i] * deltas[j] * Distance[i, j]

    return numpy.log10(res)

def CalculateSulfurNumber(mol):
    
    return _CalculateElementNumber(mol, AtomicNumber=16)

def _CalculateBasakICn(mol, NumPath=1):

    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    TotalPath = Chem.FindAllPathsOfLengthN(Hmol, NumPath, useBonds=0, useHs=1)
    if len(TotalPath) == 0:
        BasakIC = 0.0
    else:
        IC = {}
        for i in range(nAtoms):
            temp = []
            at = Hmol.GetAtomWithIdx(i)
            temp.append(at.GetAtomicNum())
            for index in TotalPath:
                if i == index[0]:
                    temp.append(
                        [Hmol.GetAtomWithIdx(kk).GetAtomicNum() for kk in index[1:]]
                    )
                if i == index[-1]:
                    cds = list(index)
                    cds.reverse()
                    temp.append(
                        [Hmol.GetAtomWithIdx(kk).GetAtomicNum() for kk in cds[1:]]
                    )
            # print temp

            IC[str(i)] = temp
        cds = []
        for value in IC.values():
            cds.append(value)
        kkk = list(range(len(cds)))
        aaa = copy.deepcopy(kkk)
        res = []
        for i in aaa:
            if i in kkk:
                jishu = 0
                kong = []
                temp1 = cds[i]
                for j in aaa:
                    if cds[j] == temp1:
                        jishu = jishu + 1
                        kong.append(j)
                for ks in kong:
                    kkk.remove(ks)
                res.append(jishu)

        # print res
        BasakIC = _CalculateEntropy(numpy.array(res, numpy.float) / sum(res))

    return BasakIC

def CalculateBasakIC6(mol):
    return _CalculateBasakICn(mol, NumPath=7)

def CalculateBasakCIC6(mol):
    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    IC = CalculateBasakIC6(mol)
    if nAtoms <= 1:
        BasakCIC = 0.0
    else:
        BasakCIC = numpy.log2(nAtoms) - IC

    return BasakCIC

def CalculateBurdenMass(mol):
    
    temp = _GetBurdenMatrix(mol, propertylabel="m")
    temp1 = numpy.sort(temp[temp >= 0])
    temp2 = numpy.sort(numpy.abs(temp[temp < 0]))

    if len(temp1) < 8:
        temp1 = numpy.concatenate((numpy.zeros(8), temp1))
    if len(temp2) < 8:
        temp2 = numpy.concatenate((numpy.zeros(8), temp2))

    bcut = [
        "bcutm16",
        "bcutm15",
        "bcutm14",
        "bcutm13",
        "bcutm12",
        "bcutm11",
        "bcutm10",
        "bcutm9",
        "bcutm8",
        "bcutm7",
        "bcutm6",
        "bcutm5",
        "bcutm4",
        "bcutm3",
        "bcutm2",
        "bcutm1",
    ]
    bcutvalue = numpy.concatenate((temp2[-8:], temp1[-8:]))

    bcutvalue = [round(i, 3) for i in bcutvalue]
    res = dict(zip(bcut, bcutvalue))
    return res

def _CalculateElementMaxNCharge(mol, AtomicNum=6):
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum() == AtomicNum:
            res.append(float(atom.GetProp("_GasteigerCharge")))
    if res == []:
        return 0
    else:
        return round(min(res), 3)

def CalculateNMaxNCharge(mol):

    return _CalculateElementMaxNCharge(mol, AtomicNum=7)

def CalculateMolLogP2(mol):

    res = Crippen._pyMolLogP(mol)

    return round(res ** 2, 3)

def CalculateHarmonicTopoIndex(mol):
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, "d")
    nAtoms = mol.GetNumAtoms()

    res = nAtoms / sum(1.0 / deltas)

    return res

def CalculateBalaban(mol):

    adjMat = Chem.GetAdjacencyMatrix(mol)
    Distance = Chem.GetDistanceMatrix(mol)
    Nbond = mol.GetNumBonds()
    Natom = mol.GetNumAtoms()
    S = numpy.sum(Distance, axis=1)
    mu = Nbond - Natom + 1
    sumk = 0.0
    for i in range(len(Distance)):
        si = S[i]
        for j in range(i, len(Distance)):
            if adjMat[i, j] == 1:
                sumk += 1.0 / numpy.sqrt(si * S[j])
    if mu + 1 != 0:
        J = float(Nbond) / float(mu + 1) * sumk
    else:
        J = 0
    return J

def CalculateWeiner(mol):

    return 1.0 / 2 * sum(sum(Chem.GetDistanceMatrix(mol)))

def CalculateMeanWeiner(mol):

    N = mol.GetNumAtoms()
    WeinerNumber = CalculateWeiner(mol)
    return 2.0 * WeinerNumber / (N * (N - 1))

def CalculateChi0(mol):
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, "d")
    res = sum(numpy.sqrt(1.0 / deltas))
    return res

def _HKDeltas(mol, skipHs=1):
    global periodicTable
    res = []
    for atom in mol.GetAtoms():
        n = atom.GetAtomicNum()
        if n > 1:
            nV = periodicTable.GetNOuterElecs(n)
            nHs = atom.GetTotalNumHs()
            if n < 10:
                res.append(float(nV - nHs))
            else:
                res.append(float(nV - nHs) / float(n - nV - 1))
        elif not skipHs:
            res.append(0.0)
    return res

def CalculateChiv0(mol):
    deltas = _HKDeltas(mol, skipHs=0)
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, "d")
    res = sum(numpy.sqrt(1.0 / deltas))
    return res

def CalculateDeltaChi0(mol):
   
    return abs(CalculateChiv0(mol) - CalculateChi0(mol))

def CalculateMolLogP(mol):

    return round(Crippen._pyMolLogP(mol), 3)

def CalculateTotalPCharge(mol):

    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp("_GasteigerCharge")))

    if res == []:
        return 0
    else:
        cc = numpy.array(res, "d")
        return round(sum(cc[cc > 0]), 3)

def CalculateTotalNCharge(mol):

    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp("_GasteigerCharge")))

    if res == []:
        return 0
    else:
        cc = numpy.array(res, "d")
        return round(sum(cc[cc < 0]), 3)

def CalculateTPSA(mol):

    return round(MS.TPSA(mol), 3)

def _CalculateElementSumSquareCharge(mol, AtomicNum=6):

    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum() == AtomicNum:
            res.append(float(atom.GetProp("_GasteigerCharge")))
    if res == []:
        return 0
    else:
        return round(sum(numpy.square(res)), 3)

def CalculateHSumSquareCharge(mol):


    return _CalculateElementSumSquareCharge(mol, AtomicNum=1)

def CalculateHdonorNumber(mol):

    return LPK.NumHDonors(mol)

def CalculateTPSA(mol):

    return round(MS.TPSA(mol), 3)

def CalculateBurdenElectronegativity(mol):
    """
    #################################################################
    Calculate Burden descriptors based on atomic electronegativity.

    res-->dict type with 16 descriptors
    #################################################################
    """
    temp = _GetBurdenMatrix(mol, propertylabel="En")
    temp1 = numpy.sort(temp[temp >= 0])
    temp2 = numpy.sort(numpy.abs(temp[temp < 0]))

    if len(temp1) < 8:
        temp1 = numpy.concatenate((numpy.zeros(8), temp1))
    if len(temp2) < 8:
        temp2 = numpy.concatenate((numpy.zeros(8), temp2))

    bcut = [
        "bcute16",
        "bcute15",
        "bcute14",
        "bcute13",
        "bcute12",
        "bcute11",
        "bcute10",
        "bcute9",
        "bcute8",
        "bcute7",
        "bcute6",
        "bcute5",
        "bcute4",
        "bcute3",
        "bcute2",
        "bcute1",
    ]
    bcutvalue = numpy.concatenate((temp2[-8:], temp1[-8:]))

    bcutvalue = [round(i, 3) for i in bcutvalue]
    res = dict(zip(bcut, bcutvalue))
    return res



def GetBurden(mol):

    bcut = {}
    bcut.update(CalculateBurdenMass(mol))
    bcut.update(CalculateBurdenVDW(mol))
    bcut.update(CalculateBurdenElectronegativity(mol))
    bcut.update(CalculateBurdenPolarizability(mol))
    return bcut

def CalculateEstateValue(mol):
    temp = ESFP.FingerprintMol(mol)
    res = {}
    for i, j in enumerate(temp[1]):
        res["S" + str(i + 1)] = round(j, 3)

    return res

def GetEstate(mol):
    result = {}
    result.update(CalculateEstateValue(mol))
    return result

def CalculateSLOGPVSA(mol, bins=None):
    temp = MOE.SlogP_VSA_(mol, bins, force=1)
    res = {}
    for i, j in enumerate(temp):
        res["slogPVSA" + str(i)] = round(j, 3)
    return res

def CalculateSMRVSA(mol, bins=None):
    temp = MOE.SMR_VSA_(mol, bins, force=1)
    res = {}
    for i, j in enumerate(temp):
        res["MRVSA" + str(i)] = round(j, 3)
    return res

def CalculatePEOEVSA(mol, bins=None):
    temp = MOE.PEOE_VSA_(mol, bins, force=1)
    res = {}
    for i, j in enumerate(temp):
        res["PEOEVSA" + str(i)] = round(j, 3)
    return res


def CalculateEstateVSA(mol, bins=None):
    temp = EVSA.EState_VSA_(mol, bins, force=1)
    res = {}
    for i, j in enumerate(temp):
        res["EstateVSA" + str(i)] = round(j, 3)
    return res

def CalculateVSAEstate(mol, bins=None):
    temp = EVSA.VSA_EState_(mol, bins, force=1)
    res = {}
    for i, j in enumerate(temp):
        res["VSAEstate" + str(i)] = round(j, 3)
    return res

def GetMOE(mol):
    result = {}
    result.update(CalculateSLOGPVSA(mol, bins=None))
    result.update(CalculateSMRVSA(mol, bins=None))
    result.update(CalculatePEOEVSA(mol, bins=None))
    result.update(CalculateEstateVSA(mol, bins=None))
    result.update(CalculateVSAEstate(mol, bins=None))
    return result


def GetCaco(mol):
    result = {}
    for DesLabel in _Caco.keys():
        result[DesLabel] = round(_Caco[DesLabel](mol), 3)

    return result

def _CalculateMoranAutocorrelation(mol, lag=1, propertylabel="m"):
    Natom = mol.GetNumAtoms()
    prolist = []
    for i in mol.GetAtoms():
        temp = GetRelativeAtomicProperty(i.GetSymbol(), propertyname=propertylabel)
        prolist.append(temp)

    aveweight = sum(prolist) / Natom

    tempp = [numpy.square(x - aveweight) for x in prolist]

    GetDistanceMatrix = Chem.GetDistanceMatrix(mol)
    res = 0.0
    index = 0
    for i in range(Natom):
        for j in range(Natom):
            if GetDistanceMatrix[i, j] == lag:
                atom1 = mol.GetAtomWithIdx(i)
                atom2 = mol.GetAtomWithIdx(j)
                temp1 = GetRelativeAtomicProperty(
                    element=atom1.GetSymbol(), propertyname=propertylabel
                )
                temp2 = GetRelativeAtomicProperty(
                    element=atom2.GetSymbol(), propertyname=propertylabel
                )
                res = res + (temp1 - aveweight) * (temp2 - aveweight)
                index = index + 1
            else:
                res = res + 0.0

    if sum(tempp) == 0 or index == 0:
        result = 0
    else:
        result = (res / index) / (sum(tempp) / Natom)

    return round(result, 3)


def CalculateMoranAutoMass(mol):
    res = {}
    for i in range(8):
        res["MATSm" + str(i + 1)] = _CalculateMoranAutocorrelation(
            mol, lag=i + 1, propertylabel="m"
        )

    return res


def CalculateMoranAutoVolume(mol):
    res = {}
    for i in range(8):
        res["MATSv" + str(i + 1)] = _CalculateMoranAutocorrelation(
            mol, lag=i + 1, propertylabel="V"
        )

    return res


def CalculateMoranAutoElectronegativity(mol):
    res = {}
    for i in range(8):
        res["MATSe" + str(i + 1)] = _CalculateMoranAutocorrelation(
            mol, lag=i + 1, propertylabel="En"
        )

    return res


def CalculateMoranAutoPolarizability(mol):
    res = {}
    for i in range(8):
        res["MATSp" + str(i + 1)] = _CalculateMoranAutocorrelation(
            mol, lag=i + 1, propertylabel="alapha"
        )

    return res

def GetMoranAuto(mol):
    res = {}
    res.update(CalculateMoranAutoMass(mol))
    res.update(CalculateMoranAutoVolume(mol))
    res.update(CalculateMoranAutoElectronegativity(mol))
    res.update(CalculateMoranAutoPolarizability(mol))
    return res

def CalculateBasakIC1(mol):
    return _CalculateBasakICn(mol, NumPath=2)

def CalculateAromaticBondNumber(mol):
    i = 0
    for bond in mol.GetBonds():

        if bond.GetBondType().name == "AROMATIC":
            i = i + 1

    return i

def CalculateAverageMolWeight(mol):
    MolWeight = 0
    for atom in mol.GetAtoms():
        MolWeight = MolWeight + atom.GetMass()

    return MolWeight / mol.GetNumAtoms()

def CalculateHydrophilicityFactor(mol):
    nheavy = mol.GetNumAtoms(onlyExplicit=1)
    nc = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            nc = nc + 1
    nhy = 0
    for atom in mol.GetAtoms():
        if (
            atom.GetAtomicNum() == 7
            or atom.GetAtomicNum() == 8
            or atom.GetAtomicNum() == 16
        ):
            atomn = atom.GetNeighbors()
            for i in atomn:
                if i.GetAtomicNum() == 1:
                    nhy = nhy + 1

    res = (
        (1 + nhy) * math.log((1 + nhy), 2)
        + nc * (1.0 / nheavy * math.log(1.0 / nheavy, 2))
        + math.sqrt((nhy + 0.0) / (nheavy ** 2))
    )
    return round(res, 3)

def _CalculatePathN(mol, PathLength=2):
    return len(Chem.FindAllPathsOfLengthN(mol, PathLength, useBonds=1))

def CalculatePath6(mol):
    return _CalculatePathN(mol, 6)

def _CalculateChinp(mol, NumPath=2):
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    for path in Chem.FindAllPathsOfLengthN(mol, NumPath + 1, useBonds=0):
        cAccum = 1.0
        for idx in path:
            cAccum *= deltas[idx]
        if cAccum:
            accum += 1.0 / numpy.sqrt(cAccum)
    return accum

def CalculateChi10p(mol):
    return _CalculateChinp(mol, NumPath=10)

def CalculateMolWeight(mol):
    MolWeight = 0
    for atom in mol.GetAtoms():
        MolWeight = MolWeight + atom.GetMass()

    return MolWeight

def CalculateRelativeNCharge(mol):
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp("_GasteigerCharge")))

    if res == []:
        return 0
    else:
        cc = numpy.array(res, "d")
        if sum(cc[cc < 0]) == 0:
            return 0
        else:
            return round(min(res) / sum(cc[cc < 0]), 3)

def CalculateHacceptorNumber(mol):
    return LPK.NumHAcceptors(mol)

def _CalculateChivnp(mol, NumPath=1):
    accum = 0.0
    deltas = _HKDeltas(mol, skipHs=0)
    for path in Chem.FindAllPathsOfLengthN(mol, NumPath + 1, useBonds=0):
        cAccum = 1.0
        for idx in path:
            cAccum *= deltas[idx]
        if cAccum:
            accum += 1.0 / numpy.sqrt(cAccum)
    return accum

def CalculateChiv4p(mol):
    return _CalculateChivnp(mol, NumPath=4)

def CalculateChiv1(mol):
    return _CalculateChivnp(mol, NumPath=1)

def CalculateChiv9p(mol):
    return _CalculateChivnp(mol, NumPath=9)

def _CalculateElementMaxPCharge(mol, AtomicNum=6):
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum() == AtomicNum:
            res.append(float(atom.GetProp("_GasteigerCharge")))

    if res == []:
        return 0
    else:
        return round(max(res), 3)

def CalculateOMaxPCharge(mol):
    return _CalculateElementMaxPCharge(mol, AtomicNum=8)

def CalculateBasakCIC0(mol):
    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    IC = CalculateBasakIC0(mol)
    if nAtoms <= 1:
        BasakCIC = 0.0
    else:
        BasakCIC = numpy.log2(nAtoms) - IC

    return BasakCIC

def CalculateCMaxPCharge(mol):
    return _CalculateElementMaxPCharge(mol, AtomicNum=6)

def CalculateCSumSquareCharge(mol):
    return _CalculateElementSumSquareCharge(mol, AtomicNum=6)

def CalculateGeometricTopoIndex(mol):
    nAtoms = mol.GetNumAtoms()
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, "d")

    temp = numpy.prod(deltas)
    res = numpy.power(temp, 1.0 / nAtoms)

    return res

def CalculateBasakIC2(mol):
    return _CalculateBasakICn(mol, NumPath=3)

def CalculateBasakCIC2(mol):
    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    IC = CalculateBasakIC2(mol)
    if nAtoms <= 1:
        BasakCIC = 0.0
    else:
        BasakCIC = numpy.log2(nAtoms) - IC

    return BasakCIC

def CalculateAllMaxPCharge(mol):
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp("_GasteigerCharge")))
    if res == []:
        return 0
    else:
        return round(max(res), 3)

def CalculateAllMaxNCharge(mol):
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp("_GasteigerCharge")))
    if res == []:
        return 0
    else:
        return round(min(res), 3)

def CalculateSubmolPolarityPara(mol):
    return round(CalculateAllMaxPCharge(mol) - CalculateAllMaxNCharge(mol), 3)

def CalculateOMaxNCharge(mol):
    return _CalculateElementMaxNCharge(mol, AtomicNum=8)

def CalculateCarbonNumber(mol):
    return _CalculateElementNumber(mol, AtomicNumber=6)

def CalculateArithmeticTopoIndex(mol):
    nAtoms = mol.GetNumAtoms()
    nBonds = mol.GetNumBonds()
    res = 2.0 * nBonds / nAtoms
    return res

def CalculateBertzCT(mol):
    temp = GD.BertzCT(mol)
    if temp > 0:
        return numpy.log10(temp)
    else:
        return "NaN"

def _HallKierAlpha(mol):
    alphaSum = 0.0
    rC = PeriodicTable.nameTable["C"][5]
    for atom in mol.GetAtoms():
        atNum = atom.GetAtomicNum()
        if not atNum:
            continue
        symb = atom.GetSymbol()
        alphaV = PeriodicTable.hallKierAlphas.get(symb, None)
        if alphaV is not None:
            hyb = atom.GetHybridization() - 2
            if hyb < len(alphaV):
                alpha = alphaV[hyb]
                if alpha is None:
                    alpha = alphaV[-1]
            else:
                alpha = alphaV[-1]
        else:
            rA = PeriodicTable.nameTable[symb][5]
            alpha = rA / rC - 1
        alphaSum += alpha
    return alphaSum

def CalculateKappaAlapha3(mol):
    P3 = len(Chem.FindAllPathsOfLengthN(mol, 3))
    A = mol.GetNumAtoms(onlyExplicit=1)
    alpha = _HallKierAlpha(mol)
    denom = P3 + alpha
    if denom:
        if A % 2 == 1:
            kappa = (A + alpha - 1) * (A + alpha - 3) ** 2 / denom ** 2
        else:
            kappa = (A + alpha - 3) * (A + alpha - 2) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)

def _CalculateBondNumber(mol, bondtype="SINGLE"):
    i = 0
    for bond in mol.GetBonds():

        if bond.GetBondType().name == bondtype:
            i = i + 1

    return i

def CalculateUnsaturationIndex(mol):
    nd = _CalculateBondNumber(mol, bondtype="DOUBLE")
    nt = _CalculateBondNumber(mol, bondtype="TRIPLE")
    na = _CalculateBondNumber(mol, bondtype="AROMATIC")
    res = math.log((1 + nd + nt + na), 2)

    return round(res, 3)

def CalculateLocalDipoleIndex(mol):
    GMCharge.ComputeGasteigerCharges(mol, iter_step)
    res = []
    for atom in mol.GetAtoms():
        res.append(float(atom.GetProp("_GasteigerCharge")))
    cc = [
        numpy.absolute(res[x.GetBeginAtom().GetIdx()] - res[x.GetEndAtom().GetIdx()])
        for x in mol.GetBonds()
    ]
    B = len(mol.GetBonds())

    return round(sum(cc) / B, 3)

def CalculateNSumSquareCharge(mol):
    return _CalculateElementSumSquareCharge(mol, AtomicNum=7)

def CalculateRelativePCharge(mol):
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp("_GasteigerCharge")))
    if res == []:
        return 0
    else:
        cc = numpy.array(res, "d")
        if sum(cc[cc > 0]) == 0:
            return 0
        else:
            return round(max(res) / sum(cc[cc > 0]), 3)

def CalculateMeanNCharge(mol):
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp("_GasteigerCharge")))

    if res == []:
        return 0
    else:
        cc = numpy.array(res, "d")
        return round(numpy.mean(cc[cc < 0]), 3)

def CalculateDeltaChi3(mol):
    return abs(_CalculateChivnp(mol, NumPath=3) - _CalculateChinp(mol, NumPath=3))

def CalculateBasakSIC1(mol):
    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    IC = CalculateBasakIC1(mol)
    if nAtoms <= 1:
        BasakSIC = 0.0
    else:
        BasakSIC = IC / numpy.log2(nAtoms)
    return BasakSIC

def CalculateNMaxPCharge(mol):
    return _CalculateElementMaxPCharge(mol, AtomicNum=7)

def CalculateDoubleBondNumber(mol):
    i = 0
    for bond in mol.GetBonds():

        if bond.GetBondType().name == "DOUBLE":
            i = i + 1

    return i

def CalculateSingleBondNumber(mol):
    i = 0
    for bond in mol.GetBonds():

        if bond.GetBondType().name == "SINGLE":
            i = i + 1

    return i

def CalculateHeteroNumber(mol):
    i = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 or atom.GetAtomicNum() == 1:
            i = i + 1

    return mol.GetNumAtoms() - i

def _AtomHKDeltas(atom, skipHs=0):
    global periodicTable
    res = []
    n = atom.GetAtomicNum()
    if n > 1:
        nV = periodicTable.GetNOuterElecs(n)
        nHs = atom.GetTotalNumHs()
        if n < 10:
            res.append(float(nV - nHs))
        else:
            res.append(float(nV - nHs) / float(n - nV - 1))
    elif not skipHs:
        res.append(0.0)
    return res

def CalculateChiv3c(mol):
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    patt = Chem.MolFromSmarts("*~*(~*)~*")
    HPatt = mol.GetSubstructMatches(patt)
    # print HPatt
    for cluster in HPatt:
        deltas = [_AtomHKDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas != []:
            deltas1 = numpy.array(deltas, numpy.float)
            accum = accum + 1.0 / numpy.sqrt(deltas1.prod())
    return accum

def CalculateChi3c(mol):
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    patt = Chem.MolFromSmarts("*~*(~*)~*")
    HPatt = mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas = [mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas != []:
            deltas1 = numpy.array(deltas, numpy.float)
            accum = accum + 1.0 / numpy.sqrt(deltas1.prod())
    return accum

def CalculateChi4pc(mol):
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    patt = Chem.MolFromSmarts("*~*(~*)~*~*")
    HPatt = mol.GetSubstructMatches(patt)
    # print HPatt
    for cluster in HPatt:
        deltas = [mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas != []:
            deltas1 = numpy.array(deltas, numpy.float)
            accum = accum + 1.0 / numpy.sqrt(deltas1.prod())
    return accum

def CalculateDeltaChi3c4pc(mol):
    return abs(CalculateChi3c(mol) - CalculateChi4pc(mol))

def CalculateRingNumber(mol):
    return Chem.GetSSSR(mol)

def CalculateBasakIC3(mol):
    return _CalculateBasakICn(mol, NumPath=4)

def CalculateBasakCIC3(mol):
    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    IC = CalculateBasakIC3(mol)
    if nAtoms <= 1:
        BasakCIC = 0.0
    else:
        BasakCIC = numpy.log2(nAtoms) - IC

    return BasakCIC

def CalculateHeavyAtomNumber(mol):
    return mol.GetNumAtoms(onlyExplicit=1)

def CalculateAllAtomNumber(mol):
    return Chem.AddHs(mol).GetNumAtoms()

def CalculateMZagreb2(mol):
    cc = [
        x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree()
        for x in mol.GetBonds()
    ]
    while 0 in cc:
        cc.remove(0)
    cc = numpy.array(cc, "d")
    res = sum((1.0 / cc) ** 2)
    return res

def CalculateKappaAlapha1(mol):
    P1 = mol.GetNumBonds(onlyHeavy=1)
    A = mol.GetNumAtoms(onlyExplicit=1)
    alpha = _HallKierAlpha(mol)
    denom = P1 + alpha
    if denom:
        kappa = (A + alpha) * (A + alpha - 1) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)

def CalculateKappaAlapha2(mol):
    P2 = len(Chem.FindAllPathsOfLengthN(mol, 2))
    A = mol.GetNumAtoms(onlyExplicit=1)
    alpha = _HallKierAlpha(mol)
    denom = P2 + alpha
    if denom:
        kappa = (A + alpha - 1) * (A + alpha - 2) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)

def CalculateFlexibility(mol):
    kappa1 = CalculateKappaAlapha1(mol)
    kappa2 = CalculateKappaAlapha2(mol)
    A = mol.GetNumAtoms(onlyExplicit=1)
    phi = kappa1 * kappa2 / (A + 0.0)
    return phi

def CalculateKappa2(mol):
    P2 = len(Chem.FindAllPathsOfLengthN(mol, 2))
    A = mol.GetNumAtoms(onlyExplicit=1)

    denom = P2 + 0.0
    if denom:
        kappa = (A - 1) * (A - 2) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)

def CalculateOSumSquareCharge(mol):
    return _CalculateElementSumSquareCharge(mol, AtomicNum=8)

def CalculateAllSumSquareCharge(mol):
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp("_GasteigerCharge")))

    if res == []:
        return 0
    else:
        return round(sum(numpy.square(res)), 3)

def CalculateNitrogenNumber(mol):
    return _CalculateElementNumber(mol, AtomicNumber=7)

def CalculateMZagreb1(mol):
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, "d")
    res = sum((1.0 / deltas) ** 2)
    return res

def CalculateKappa3(mol):
    P3 = len(Chem.FindAllPathsOfLengthN(mol, 3))
    A = mol.GetNumAtoms(onlyExplicit=1)

    denom = P3 + 0.0
    if denom:
        if A % 2 == 1:
            kappa = (A - 1) * (A - 3) ** 2 / denom ** 2
        else:
            kappa = (A - 3) * (A - 2) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)

def CalculateHMaxPCharge(mol):
    return _CalculateElementMaxPCharge(mol, AtomicNum=1)

def caco2(filename):
	df = pd.read_csv(filename,engine='python').iloc[:, 0]  
	#Import the column of Smiles formula
	smiles_list = numpy.array(df)
	b=[]
	for i in range(len(smiles_list)):
		m = Chem.MolFromSmiles(smiles_list[i])
		a=[]
		a.append(round(CalculateCarbonNumber(m),3))
		a.append(round(CalculateBasakIC0(m),3))
		a.append(GetBurden(m).get('bcutp1'))
		a.append(GetBurden(m).get('bcutv10'))
		a.append(round(CalculateGutmanTopo(m),3))
		a.append(round(CalculateSulfurNumber(m),3))
		a.append(round(CalculateBasakCIC6(m),3))
		a.append(GetBurden(m).get('bcutm12'))
		a.append(GetEstate(m).get('S34'))
		a.append(GetBurden(m).get('bcutp8'))
		a.append(GetMOE(m).get('slogPVSA2'))
		a.append(round(CalculateNMaxNCharge(m),3))
		a.append(round(CalculateMolLogP2(m),3))
		a.append(GetBurden(m).get('bcutm1'))
		a.append(GetMOE(m).get('EstateVSA9'))
		a.append(GetMOE(m).get('slogPVSA1'))
		a.append(round(CalculateHarmonicTopoIndex(m),3))
		a.append(round(CalculateBalaban(m),3))
		a.append(round(CalculateMeanWeiner(m),3))
		a.append(GetEstate(m).get('S7'))
		a.append(round(CalculateDeltaChi0(m),3))
		a.append(GetMOE(m).get('MRVSA1'))
		a.append(round(CalculateMolLogP(m),3))
		a.append(round(CalculateTotalPCharge(m),3))
		a.append(GetMOE(m).get('PEOEVSA0'))
		a.append(round(CalculateTotalNCharge(m),3))
		a.append(GetEstate(m).get('S13'))
		a.append(round(CalculateTPSA(m),3))
		a.append(round(CalculateHSumSquareCharge(m),3))
		a.append(round(CalculateHdonorNumber(m),3))
		b.append(a)
	des_list=numpy.array(b)
	return des_list

def logS(filename):
    df = pd.read_csv(filename, engine='python').iloc[:, 0]
    # Import the column of Smiles formula
    smiles_list = numpy.array(df)
    b = []
    for i in range(len(smiles_list)):
        m = Chem.MolFromSmiles(smiles_list[i])
        a = []
        a.append(GetMoranAuto(m).get('MATSm2'))
        a.append(100)
        a.append(round(CalculateGutmanTopo(m), 3))
        a.append(round(CalculateBasakIC1(m), 3))
        a.append(round(CalculateAromaticBondNumber(m), 3))
        a.append(GetMoranAuto(m).get('MATSm1'))
        a.append(round(CalculateSulfurNumber(m), 3))
        a.append(round(CalculateTotalPCharge(m), 3))
        a.append(GetMOE(m).get('slogPVSA7'))
        a.append(GetBurden(m).get('bcutp1'))
        a.append(round(CalculateAverageMolWeight(m), 3))
        a.append(round(CalculateTotalNCharge(m), 3))
        a.append(GetMOE(m).get('MRVSA9'))
        a.append(GetBurden(m).get('bcutp3'))
        a.append(round(CalculateBasakIC0(m), 3))
        a.append(round(CalculateMeanWeiner(m), 3))
        a.append(round(CalculateHydrophilicityFactor(m), 3))
        a.append(GetBurden(m).get('bcutv10'))
        a.append(GetMOE(m).get('MRVSA6'))
        a.append(round(CalculatePath6(m), 3))
        a.append(GetBurden(m).get('bcutm1'))
        a.append(GetBurden(m).get('bcutm8'))
        a.append(GetMOE(m).get('slogPVSA1'))
        a.append(100)
        a.append(round(CalculateChi10p(m), 3))
        a.append(round(CalculateTPSA(m), 3))
        a.append(round(CalculateMolWeight(m), 3))
        a.append(round(CalculateRelativeNCharge(m), 3))
        a.append(round(CalculateHacceptorNumber(m), 3))
        a.append(GetBurden(m).get('bcutp5'))
        a.append(round(CalculateChiv4p(m), 3))
        a.append(GetBurden(m).get('bcutm2'))
        a.append(round(CalculateChiv1(m), 3))
        a.append(GetBurden(m).get('bcutm3'))
        a.append(round(CalculateChiv9p(m), 3))
        a.append(round(CalculateCarbonNumber(m), 3))
        a.append(GetBurden(m).get('bcutm4'))
        a.append(GetMOE(m).get('PEOEVSA5'))
        a.append(round(CalculateMolLogP2(m), 3))
        a.append(round(CalculateMolLogP(m), 3))
        b.append(a)
    des_list = numpy.array(b)
    return des_list

def logD(filename):
    df = pd.read_csv(filename, engine='python').iloc[:, 0]
    # Import the column of Smiles formula
    smiles_list = numpy.array(df)
    b = []
    for i in range(len(smiles_list)):
        m = Chem.MolFromSmiles(smiles_list[i])
        a = []
        a.append(GetMoranAuto(m).get('MATSe5'))
        a.append(GetMOE(m).get('PEOEVSA9'))
        a.append(GetMOE(m).get('EstateVSA7'))
        a.append(GetEstate(m).get('S13'))
        a.append(GetMOE(m).get('EstateVSA0'))
        a.append(round(CalculateChiv4p(m), 3))
        a.append(GetEstate(m).get('S28'))
        a.append(round(CalculateMeanWeiner(m), 3))
        a.append(round(CalculateOMaxPCharge(m), 3))
        a.append(GetBurden(m).get('bcutp2'))
        a.append(GetMOE(m).get('EstateVSA4'))
        a.append(GetMoranAuto(m).get('MATSe1'))
        a.append(round(CalculatePath6(m), 3))
        a.append(round(CalculateHarmonicTopoIndex(m), 3))
        a.append(GetEstate(m).get('S24'))
        a.append(round(CalculateBasakCIC0(m), 3))
        a.append(round(CalculateCMaxPCharge(m), 3))
        a.append(round(CalculateCSumSquareCharge(m), 3))
        a.append(round(CalculateGeometricTopoIndex(m), 3))
        a.append(round(CalculateTPSA(m), 3))
        a.append(round(CalculateGeometricTopoIndex(m), 3))
        a.append(GetBurden(m).get('bcutm11'))
        a.append(round(CalculateBasakCIC2(m), 3))
        a.append(round(CalculateBalaban(m), 3))
        a.append(GetEstate(m).get('S34'))
        a.append(GetMOE(m).get('PEOEVSA5'))
        a.append(round(CalculateHydrophilicityFactor(m), 3))
        a.append(round(CalculateSubmolPolarityPara(m), 3))
        a.append(GetEstate(m).get('S36'))
        a.append(GetEstate(m).get('S9'))
        a.append(GetEstate(m).get('S16'))
        a.append(GetMOE(m).get('MRVSA4'))
        a.append(round(CalculateMolLogP2(m), 3))
        a.append(round(CalculateOMaxNCharge(m), 3))
        a.append(round(CalculateMolLogP(m), 3))
        b.append(a)
    des_list = numpy.array(b)
    return des_list

def PPB(filename):
    df = pd.read_csv(filename, engine='python').iloc[:, 0]
    # Import the column of Smiles formula
    smiles_list = numpy.array(df)
    b = []
    for i in range(len(smiles_list)):
        m = Chem.MolFromSmiles(smiles_list[i])
        a = []
        a.append(round(CalculateCarbonNumber(m), 3))
        a.append(round(CalculateAromaticBondNumber(m), 3))
        a.append(round(CalculateGutmanTopo(m), 3))
        a.append(round(CalculateMeanWeiner(m), 3))
        a.append(round(CalculateGeometricTopoIndex(m), 3))
        a.append(round(CalculateArithmeticTopoIndex(m), 3))
        a.append(round(CalculateBertzCT(m), 3))
        a.append(round(CalculateBalaban(m), 3))
        a.append(round(CalculateKappaAlapha3(m), 3))
        a.append(GetMoranAuto(m).get('MATSv1'))
        a.append(GetMoranAuto(m).get('MATSe1'))
        a.append(round(CalculateHydrophilicityFactor(m), 3))
        a.append(round(CalculateMolLogP(m), 3))
        a.append(round(CalculateMolLogP2(m), 3))
        a.append(round(CalculateUnsaturationIndex(m), 3))
        a.append(round(CalculateHSumSquareCharge(m), 3))
        a.append(round(CalculateLocalDipoleIndex(m), 3))
        a.append(round(CalculateNSumSquareCharge(m), 3))
        a.append(round(CalculateRelativePCharge(m), 3))
        a.append(round(CalculateMeanNCharge(m), 3))
        a.append(GetMOE(m).get('PEOEVSA10'))
        a.append(GetMOE(m).get('PEOEVSA0'))
        a.append(GetMOE(m).get('PEOEVSA6'))
        a.append(GetMOE(m).get('PEOEVSA5'))
        a.append(GetMOE(m).get('PEOEVSA4'))
        a.append(GetMOE(m).get('slogPVSA10'))
        a.append(GetMOE(m).get('MRVSA6'))
        a.append(GetMOE(m).get('slogPVSA0'))
        a.append(GetMOE(m).get('slogPVSA1'))
        a.append(GetMOE(m).get('slogPVSA5'))
        b.append(a)
    des_list = numpy.array(b)
    return des_list

def VD(filename):
    df = pd.read_csv(filename, engine='python').iloc[:, 0]
    # Import the column of Smiles formula
    smiles_list = numpy.array(df)
    b = []
    for i in range(len(smiles_list)):
        m = Chem.MolFromSmiles(smiles_list[i])
        a = []
        a.append(round(CalculateGutmanTopo(m), 3))
        a.append(round(CalculateUnsaturationIndex(m), 3))
        a.append(GetMoranAuto(m).get('MATSe1'))
        a.append(GetMoranAuto(m).get('MATSp1'))
        a.append(round(CalculateChiv4p(m), 3))
        a.append(GetMoranAuto(m).get('MATSm2'))
        a.append(GetEstate(m).get('S12'))
        a.append(round(CalculateDeltaChi3(m), 3))
        a.append(0)
        a.append(GetMOE(m).get('PEOEVSA7'))
        a.append(GetBurden(m).get('bcutp1'))
        a.append(GetBurden(m).get('bcutm9'))
        a.append(round(CalculateBasakSIC1(m), 3))
        a.append(GetMOE(m).get('MRVSA6'))
        a.append(round(CalculateBasakIC1(m), 3))
        a.append(round(CalculateNMaxPCharge(m), 3))
        a.append(round(CalculateBasakCIC0(m), 3))
        a.append(GetMOE(m).get('PEOEVSA6'))
        a.append(GetMoranAuto(m).get('MATSe4'))
        a.append(GetMOE(m).get('VSAEstate8'))
        a.append(round(CalculateGeometricTopoIndex(m), 3))
        a.append(GetMOE(m).get('EstateVSA3'))
        a.append(GetMOE(m).get('MRVSA5'))
        a.append(round(CalculateMolLogP2(m), 3))
        a.append(round(CalculateTotalNCharge(m), 3))
        a.append(GetEstate(m).get('S7'))
        a.append(round(CalculateSubmolPolarityPara(m), 3))
        a.append(round(CalculateOMaxNCharge(m), 3))
        a.append(GetMOE(m).get('EstateVSA7'))
        a.append(round(CalculateMolLogP(m), 3))
        a.append(round(CalculateNMaxNCharge(m), 3))
        a.append(GetMOE(m).get('MRVSA9'))
        a.append(GetEstate(m).get('S19'))
        a.append(GetMoranAuto(m).get('MATSv2'))
        a.append(round(CalculateSulfurNumber(m), 3))
        a.append(GetEstate(m).get('S17'))
        a.append(GetEstate(m).get('S9'))
        a.append(round(CalculateDoubleBondNumber(m), 3))
        a.append(round(CalculateAverageMolWeight(m), 3))
        a.append(round(CalculateCSumSquareCharge(m), 3))
        a.append(GetMOE(m).get('EstateVSA9'))
        a.append(round(CalculateHydrophilicityFactor(m), 3))
        a.append(GetEstate(m).get('S16'))
        a.append(round(CalculateBasakIC0(m), 3))
        a.append(GetEstate(m).get('S30'))
        b.append(a)
    des_list = numpy.array(b)
    return des_list

def CL(filename):
    df = pd.read_csv(filename, engine='python').iloc[:, 0]
    # Import the column of Smiles formula
    smiles_list = numpy.array(df)
    b = []
    for i in range(len(smiles_list)):
        m = Chem.MolFromSmiles(smiles_list[i])
        a = []
        a.append(round(CalculateSulfurNumber(m), 3))
        a.append(GetMOE(m).get('VSAEstate8'))
        a.append(round(CalculateNMaxNCharge(m), 3))
        a.append(0)
        a.append(round(CalculateDoubleBondNumber(m), 3))
        a.append(GetMOE(m).get('slogPVSA2'))
        a.append(GetMoranAuto(m).get('MATSv5'))
        a.append(GetEstate(m).get('S32'))
        a.append(round(CalculateCSumSquareCharge(m), 3))
        a.append(GetBurden(m).get('bcutm4'))
        a.append(GetEstate(m).get('S9'))
        a.append(GetBurden(m).get('bcutp8'))
        a.append(round(CalculateTotalNCharge(m), 3))
        a.append(round(CalculateSingleBondNumber(m), 3))
        a.append(round(CalculateGeometricTopoIndex(m), 3))
        a.append(GetBurden(m).get('bcutp11'))
        a.append(GetEstate(m).get('S7'))
        a.append(GetMoranAuto(m).get('MATSm2'))
        a.append(round(CalculateGutmanTopo(m), 3))
        a.append(round(CalculateHeteroNumber(m), 3))
        a.append(GetMoranAuto(m).get('MATSe1'))
        a.append(round(CalculateBasakCIC0(m), 3))
        a.append(GetBurden(m).get('bcutp3'))
        a.append(0)
        a.append(GetMOE(m).get('EstateVSA9'))
        a.append(GetMoranAuto(m).get('MATSe3'))
        a.append(GetMoranAuto(m).get('MATSe5'))
        a.append(round(CalculateUnsaturationIndex(m), 3))
        a.append(GetEstate(m).get('S53'))
        a.append(round(CalculateBalaban(m), 3))
        a.append(GetBurden(m).get('bcute1'))
        a.append(GetMOE(m).get('MRVSA9'))
        a.append(GetMOE(m).get('PEOEVSA0'))
        a.append(GetMoranAuto(m).get('MATSv2'))
        a.append(0)
        a.append(round(CalculateAverageMolWeight(m), 3))
        a.append(round(CalculateBasakIC0(m), 3))
        a.append(GetEstate(m).get('S16'))
        a.append(GetBurden(m).get('bcutp1'))
        a.append(GetMOE(m).get('PEOEVSA12'))
        b.append(a)
    des_list = numpy.array(b)
    return des_list

def T(filename):
    df = pd.read_csv(filename, engine='python').iloc[:, 0]
    # Import the column of Smiles formula
    smiles_list = numpy.array(df)
    b = []
    for i in range(len(smiles_list)):
        m = Chem.MolFromSmiles(smiles_list[i])
        a = []
        a.append(GetMoranAuto(m).get('MATSv5'))
        a.append(0)
        a.append(round(CalculateChiv3c(m), 3))
        a.append(GetMOE(m).get('PEOEVSA7'))
        a.append(round(CalculateDeltaChi3c4pc(m), 3))
        a.append(GetBurden(m).get('bcutp3'))
        a.append(GetBurden(m).get('bcutm9'))
        a.append(GetMOE(m).get('EstateVSA3'))
        a.append(GetMoranAuto(m).get('MATSp1'))
        a.append(GetBurden(m).get('bcutp11'))
        a.append(GetMOE(m).get('VSAEstate7'))
        a.append(round(CalculateBasakIC0(m), 3))
        a.append(round(CalculateUnsaturationIndex(m), 3))
        a.append(round(CalculateGeometricTopoIndex(m), 3))
        a.append(round(CalculateOMaxNCharge(m), 3))
        a.append(round(CalculateBasakCIC0(m), 3))
        a.append(round(CalculateDeltaChi3(m), 3))
        a.append(GetMoranAuto(m).get('MATSp4'))
        a.append(GetBurden(m).get('bcutm4'))
        a.append(round(CalculateHarmonicTopoIndex(m), 3))
        a.append(GetMoranAuto(m).get('MATSe4'))
        a.append(round(CalculateBasakCIC6(m), 3))
        a.append(round(CalculateChiv4p(m), 3))
        a.append(GetMOE(m).get('EstateVSA9'))
        a.append(GetMoranAuto(m).get('MATSv2'))
        a.append(round(CalculateRingNumber(m), 3))
        a.append(GetBurden(m).get('bcute1'))
        a.append(GetMOE(m).get('VSAEstate8'))
        a.append(GetMOE(m).get('MRVSA9'))
        a.append(GetMOE(m).get('PEOEVSA6'))
        a.append(round(CalculateBasakSIC1(m), 3))
        a.append(GetBurden(m).get('bcutp8'))
        a.append(GetMoranAuto(m).get('MATSp6'))
        a.append(round(CalculateCSumSquareCharge(m), 3))
        a.append(round(CalculateBalaban(m), 3))
        a.append(0)
        a.append(round(CalculateBasakCIC2(m), 3))
        a.append(round(CalculateHydrophilicityFactor(m), 3))
        a.append(GetMOE(m).get('MRVSA6'))
        a.append(round(CalculateAromaticBondNumber(m), 3))
        a.append(round(CalculateSubmolPolarityPara(m), 3))
        a.append(GetMOE(m).get('EstateVSA7'))
        a.append(GetBurden(m).get('bcutv10'))
        a.append(GetEstate(m).get('S12'))
        a.append(round(CalculateMolLogP2(m), 3))
        a.append(GetBurden(m).get('bcutp2'))
        a.append(round(CalculateBasakCIC3(m), 3))
        a.append(GetEstate(m).get('S17'))
        a.append(round(CalculateMolLogP(m), 3))
        a.append(GetBurden(m).get('bcutp1'))
        b.append(a)
    des_list = numpy.array(b)
    return des_list
def hERG(filename):
    df = pd.read_csv(filename, engine='python').iloc[:, 0]
    # Import the column of Smiles formula
    smiles_list = numpy.array(df)
    b = []
    for i in range(len(smiles_list)):
        m = Chem.MolFromSmiles(smiles_list[i])
        a = []
        a.append(round(CalculateDoubleBondNumber(m), 3))
        a.append(round(CalculateSingleBondNumber(m), 3))
        a.append(round(CalculateCarbonNumber(m), 3))
        a.append(round(CalculateSulfurNumber(m), 3))
        a.append(round(CalculateAromaticBondNumber(m), 3))
        a.append(round(CalculateHdonorNumber(m), 3))
        a.append(round(CalculateHeavyAtomNumber(m), 3))
        a.append(round(CalculateHacceptorNumber(m), 3))
        a.append(round(CalculateAllAtomNumber(m), 3))
        a.append(round(CalculateRingNumber(m), 3))
        a.append(round(CalculatePath6(m), 3))
        a.append(round(CalculateGutmanTopo(m), 3))
        a.append(round(CalculateMeanWeiner(m), 3))
        a.append(round(CalculateGeometricTopoIndex(m), 3))
        a.append(round(CalculateBertzCT(m), 3))
        a.append(round(CalculateBalaban(m), 3))
        a.append(round(CalculateMZagreb2(m), 3))
        a.append(round(CalculateFlexibility(m), 3))
        a.append(round(CalculateKappa2(m), 3))
        a.append(GetMoranAuto(m).get('MATSv1'))
        a.append(GetMoranAuto(m).get('MATSv5'))
        a.append(GetMoranAuto(m).get('MATSe4'))
        a.append(GetMoranAuto(m).get('MATSe5'))
        a.append(GetMoranAuto(m).get('MATSe6'))
        a.append(round(CalculateTPSA(m), 3))
        a.append(round(CalculateHydrophilicityFactor(m), 3))
        a.append(round(CalculateMolLogP(m), 3))
        a.append(round(CalculateMolLogP2(m), 3))
        a.append(round(CalculateUnsaturationIndex(m), 3))
        a.append(round(CalculateOSumSquareCharge(m), 3))
        a.append(round(CalculateSubmolPolarityPara(m), 3))
        a.append(round(CalculateLocalDipoleIndex(m), 3))
        a.append(round(CalculateAllSumSquareCharge(m), 3))
        a.append(round(CalculateOMaxNCharge(m), 3))
        a.append(round(CalculateNMaxPCharge(m), 3))
        a.append(round(CalculateAllMaxNCharge(m), 3))
        a.append(round(CalculateMeanNCharge(m), 3))
        a.append(GetMOE(m).get('EstateVSA7'))
        a.append(GetMOE(m).get('EstateVSA0'))
        a.append(GetMOE(m).get('EstateVSA3'))
        a.append(GetMOE(m).get('PEOEVSA0'))
        a.append(GetMOE(m).get('PEOEVSA6'))
        a.append(GetMOE(m).get('MRVSA5'))
        a.append(GetMOE(m).get('MRVSA4'))
        a.append(GetMOE(m).get('MRVSA3'))
        a.append(GetMOE(m).get('MRVSA6'))
        a.append(GetMOE(m).get('slogPVSA1'))
        b.append(a)
    des_list = numpy.array(b)
    return des_list


def HHT(filename):
    df = pd.read_csv(filename, engine='python').iloc[:, 0]
    # Import the column of Smiles formula
    smiles_list = numpy.array(df)
    b = []
    for i in range(len(smiles_list)):
        m = Chem.MolFromSmiles(smiles_list[i])
        a = []
        a.append(round(CalculateDoubleBondNumber(m), 3))
        a.append(round(CalculateSingleBondNumber(m), 3))
        a.append(round(CalculateNitrogenNumber(m), 3))
        a.append(round(CalculateAromaticBondNumber(m), 3))
        a.append(round(CalculateHdonorNumber(m), 3))
        a.append(round(CalculateHeteroNumber(m), 3))
        a.append(round(CalculateHeavyAtomNumber(m), 3))
        a.append(round(CalculateRingNumber(m), 3))
        a.append(round(CalculatePath6(m), 3))
        a.append(round(CalculateGutmanTopo(m), 3))
        a.append(round(CalculateGeometricTopoIndex(m), 3))
        a.append(0)
        a.append(round(CalculateArithmeticTopoIndex(m), 3))
        a.append(round(CalculateHarmonicTopoIndex(m), 3))
        a.append(round(CalculateBertzCT(m), 3))
        a.append(round(CalculateGeometricTopoIndex(m), 3))
        a.append(round(CalculateBalaban(m), 3))
        a.append(round(CalculateMZagreb2(m), 3))
        a.append(round(CalculateMZagreb1(m), 3))
        a.append(round(CalculateFlexibility(m), 3))
        a.append(round(CalculateKappa3(m), 3))
        a.append(round(CalculateKappa2(m), 3))
        a.append(round(CalculateKappaAlapha3(m), 3))
        a.append(GetMoranAuto(m).get('MATSp4'))
        a.append(GetMoranAuto(m).get('MATSp6'))
        a.append(GetMoranAuto(m).get('MATSv3'))
        a.append(GetMoranAuto(m).get('MATSv2'))
        a.append(GetMoranAuto(m).get('MATSv5'))
        a.append(GetMoranAuto(m).get('MATSv7'))
        a.append(GetMoranAuto(m).get('MATSv6'))
        a.append(GetMoranAuto(m).get('MATSm4'))
        a.append(GetMoranAuto(m).get('MATSm5'))
        a.append(GetMoranAuto(m).get('MATSm6'))
        a.append(GetMoranAuto(m).get('MATSm2'))
        a.append(GetMoranAuto(m).get('MATSm3'))
        a.append(GetMoranAuto(m).get('MATSe4'))
        a.append(GetMoranAuto(m).get('MATSe5'))
        a.append(GetMoranAuto(m).get('MATSe6'))
        a.append(GetMoranAuto(m).get('MATSe1'))
        a.append(GetMoranAuto(m).get('MATSe2'))
        a.append(GetMoranAuto(m).get('MATSe3'))
        a.append(GetMoranAuto(m).get('MATSp3'))
        a.append(GetMoranAuto(m).get('MATSp2'))
        a.append(round(CalculateTPSA(m), 3))
        a.append(round(CalculateMolLogP(m), 3))
        a.append(round(CalculateMolLogP2(m), 3))
        a.append(round(CalculateUnsaturationIndex(m), 3))
        a.append(round(CalculateNMaxNCharge(m), 3))
        a.append(round(CalculateOSumSquareCharge(m), 3))
        a.append(round(CalculateHSumSquareCharge(m), 3))
        a.append(round(CalculateSubmolPolarityPara(m), 3))
        a.append(round(CalculateLocalDipoleIndex(m), 3))
        a.append(round(CalculateAllSumSquareCharge(m), 3))
        a.append(round(CalculateCMaxPCharge(m), 3))
        a.append(round(CalculateOMaxPCharge(m), 3))
        a.append(round(CalculateTotalPCharge(m), 3))
        a.append(round(CalculateOMaxNCharge(m), 3))
        a.append(round(CalculateTotalPCharge(m), 3))
        a.append(round(CalculateHMaxPCharge(m), 3))
        a.append(round(CalculateRelativeNCharge(m), 3))
        a.append(round(CalculateRelativePCharge(m), 3))
        a.append(round(CalculateAllMaxNCharge(m), 3))
        a.append(round(CalculateMeanNCharge(m), 3))
        a.append(GetMOE(m).get('EstateVSA9'))
        a.append(GetMOE(m).get('EstateVSA4'))
        a.append(GetMOE(m).get('EstateVSA5'))
        a.append(GetMOE(m).get('EstateVSA6'))
        a.append(GetMOE(m).get('EstateVSA7'))
        a.append(GetMOE(m).get('EstateVSA0'))
        a.append(GetMOE(m).get('EstateVSA1'))
        a.append(GetMOE(m).get('EstateVSA2'))
        a.append(GetMOE(m).get('EstateVSA3'))
        a.append(GetMOE(m).get('PEOEVSA11'))
        a.append(GetMOE(m).get('PEOEVSA2'))
        a.append(GetMOE(m).get('PEOEVSA1'))
        a.append(GetMOE(m).get('PEOEVSA7'))
        a.append(GetMOE(m).get('PEOEVSA6'))
        a.append(GetMOE(m).get('PEOEVSA5'))
        a.append(GetMOE(m).get('MRVSA5'))
        a.append(GetMOE(m).get('MRVSA4'))
        a.append(GetMOE(m).get('PEOEVSA9'))
        a.append(GetMOE(m).get('PEOEVSA8'))
        a.append(GetMOE(m).get('MRVSA3'))
        a.append(GetMOE(m).get('MRVSA2'))
        a.append(GetMOE(m).get('MRVSA9'))
        a.append(GetMOE(m).get('MRVSA6'))
        a.append(GetMOE(m).get('slogPVSA2'))
        a.append(GetMOE(m).get('slogPVSA4'))
        a.append(GetMOE(m).get('slogPVSA5'))
        b.append(a)
        des_list = numpy.array(b)
    return des_list
   
def predict(filename,savename,current_path):
	model1,model2,model3,model4,model5,model6,model7,model8,model9,model10,model11,model12,model13,model14,model15,model16,model17,model18,model19,model20,model21,model22,model23,model24,model25,model26= Model_selection(current_path)
	dimension=2048
	radius=2
	
	des_list1 = ECFP_generate(filename,radius,dimension)
	y_predict_label1 = model1.predict(des_list1)

	y_predict_label2 = model2.predict(des_list1)

	y_predict_label3 = model3.predict(des_list1)

	y_predict_label4 = model4.predict(des_list1)

	dimension=1024
	radius=2
	des_list2 = ECFP_generate(filename,radius,dimension)
	y_predict_label5 = model5.predict(des_list2)

	y_predict_label6 = model6.predict(des_list2)

	y_predict_label7 = model7.predict(des_list2)

	y_predict_label8 = model8.predict(des_list2)

	y_predict_label9 = model9.predict(des_list2)

	y_predict_label10 = model10.predict(des_list2)

	dimension=2048
	radius=1
	des_list3 = ECFP_generate(filename,radius,dimension)
	y_predict_label11 = model11.predict(des_list3)

	y_predict_label12 = model12.predict(des_list3)
	
	des_list4 = MACCS_generate(filename)
	y_predict_label13 = model13.predict(des_list4)
	y_predict_label14 = model14.predict(des_list4)
	y_predict_label15 = model15.predict(des_list4)
	y_predict_label16 = model16.predict(des_list4)
	
	dimension=2048
	radius=3
	des_list5 = ECFP_generate(filename,radius,dimension)	
	y_predict_label17 = model17.predict(des_list5)
	
	des_list6 = caco2(filename)
	y_predict_label18 = model18.predict(des_list6)

	des_list7 = logS(filename)
	y_predict_label19 = model19.predict(des_list7)

	des_list8 = logD(filename)
	y_predict_label20 = model20.predict(des_list8)

	des_list9 = PPB(filename)
	y_predict_label21 = model21.predict(des_list9)

	des_list10 = VD(filename)
	y_predict_label22 = model22.predict(des_list10)

	des_list11 = CL(filename)
	y_predict_label23 = model23.predict(des_list11)

	des_list12 = T(filename)
	y_predict_label24 = model24.predict(des_list12)

	des_list13 = hERG(filename)
	y_predict_label25 = model25.predict(des_list13)

	des_list14 = HHT(filename)
	y_predict_label26 = model26.predict(des_list14)
	
	df = pd.read_csv(filename,engine='python').iloc[:, 0]  
	f = open(savename,'w')
	csv_writer = csv.writer(f)
	csv_writer.writerow(["Smiles","CYP_inhibitor_1A2","CYP_inhibitor_2C9","CYP_inhibitor_3A4","PGP_inhibitor","CYP_substrate_1A2","CYP_substrate_2C9","CYP_substrate_2C19","CYP_substrate_3A4","CYP_substrate_2D6","CYP_inhibitor_2D6","BBB","CYP_inhibitor_2C19","Ames","F_20","HIA","SkinSen","F_30","Caco2_Permeability","LogS","logD","PPB","VD","CL","T","hERG","HHT"])
	for i in range(len(df)):
		csv_writer.writerow([df[i],y_predict_label1[i],y_predict_label2[i],y_predict_label3[i],y_predict_label4[i],y_predict_label5[i],y_predict_label6[i],y_predict_label7[i],y_predict_label8[i],y_predict_label9[i],y_predict_label10[i],y_predict_label11[i],y_predict_label12[i],y_predict_label13[i],y_predict_label14[i],y_predict_label15[i],y_predict_label16[i],y_predict_label17[i],y_predict_label18[i],y_predict_label19[i],y_predict_label20[i],y_predict_label21[i],y_predict_label22[i],y_predict_label23[i],y_predict_label24[i],y_predict_label25[i],y_predict_label26[i]])
	f.close()