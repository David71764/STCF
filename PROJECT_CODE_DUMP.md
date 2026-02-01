# 1. 程序结构

```markdown
TauMuGammaAna/
├── Submission/                 # 作业提交脚本与配置文件
│   ├── Sub_Ana/                # 分析阶段 (Analysis)
│   │   ├── joboptions_ana.py   # 分析作业选项配置
│   │   ├── runANA.sh           # 运行分析的脚本
│   │   └── SubJob_Ana.sh       # 提交分析作业的脚本
│   ├── Sub_Rec/                # 重建阶段 (Reconstruction)
│   │   ├── joboptions_rec.py   # 重建作业选项配置
│   │   ├── runREC.sh           # 运行重建的脚本
│   │   └── SubJob_Rec.sh       # 提交重建作业的脚本
│   └── Sub_Sim/                # 模拟阶段 (Simulation)
│       ├── decay.dec           # 衰减卡文件 (Decay card)
│       ├── joboptions_sim.py   # 模拟作业选项配置
│       ├── runSIM.sh           # 运行模拟的脚本
│       └── SubJob_Sim.sh       # 提交模拟作业的脚本
├── TauMuGammaAna/              # 核心算法源代码
│   ├── CMakeLists.txt          # 项目构建配置文件
│   ├── include/                # 头文件目录
│   │   └── TauMuGammaAlg.h     # 算法类定义
│   ├── python/                 # Python 接口或辅助脚本
│   │   └── TauMuGammaAna/
│   │       └── __init__.py
│   └── src/                    # 源文件目录
│       └── TauMuGammaAlg.cc    # 算法逻辑实现
├── plotting/                   # 作图与复现论文图的脚本
│   ├── __init__.py
│   ├── io.py
│   ├── selection.py
│   ├── plot_paper_figs.py
│   ├── samples.json
│   └── ellipse_params.json
└── .vscode/                    # 编辑器配置
    ├── c_cpp_properties.json
    ├── launch.json
    └── settings.json
```

# 2. 具体程序内容

## Submission/Sub_Ana/joboptions_ana.py

```python
#!/usr/bin/env python3
import os
import Sniper

topdir = os.getenv("OFFLINETOP")

AlgTask = Sniper.Task("TauMuGamma_Ana")
AlgTask.setLogLevel(2)
AlgTask.setEvtMax(MAXEVT)   # MAXEVT 必须是数字占位符（不要加引号）

# ---------------- Services ----------------
import PodioDataSvc
dsvc = AlgTask.createSvc("PodioDataSvc")

import PodioSvc
Isvc = AlgTask.createSvc("PodioInputSvc/InputSvc")
Isvc.property("InputFile").set("INPUT")   # INPUT 必须在引号里（字符串）

# ✅ 分析阶段不需要 PodioOutputSvc（否则会并发写同一个 AnaTmp.root 造成冲突）
# outsvc = AlgTask.createSvc("PodioOutputSvc/OutputSvc")
# outsvc.property("OutputFile").set("AnaTmp.root")

import RootWriter
osvc = AlgTask.createSvc("RootWriter/hSvc")
osvc.property("Output").set({"Fkey": "ANAOUT"})  # ANAOUT 必须在引号里（字符串路径）

# Global PID（按需）
import GlobalPID
GlobalPID = AlgTask.createSvc("GlobalPIDSvc")
model_path = None
if topdir:
    model_path = os.path.join(topdir, "Analysis/GlobalPID/src/xgb.model")
if model_path and os.path.exists(model_path):
    GlobalPID.property("SetModelPath").set(model_path)
    GlobalPID.property("SetMethod").set("XGBoost")
else:
    print("[WARN] GlobalPID model not found, fallback to default PID")

# ---------------- Algorithm ----------------
import TauMuGammaAna
TauMuGamma = AlgTask.createAlg("TauMuGammaAlg")

TauMuGamma.property("cmsEnergy").set(ECMS)

# 下面这些 cut 你之前按 Table/正文写的逻辑保持不变
TauMuGamma.property("cosThetaAccCut").set(0.93)
TauMuGamma.property("cosThetaMuGammaCut").set(-0.35)
TauMuGamma.property("chargedPmin").set(0.2)
TauMuGamma.property("muPidMin").set(0.2)
TauMuGamma.property("muPidSlack").set(0.05)

# 运行模式与严格度（与论文对齐）
TauMuGamma.property("analysisMode").set(0)  # 0=PaperStrict, 1=FullSimRobust
TauMuGamma.property("requireExactPhotonMultiplicity").set(True)
TauMuGamma.property("useEoverPTagID").set(True)

# photon acceptance (paper)
TauMuGamma.property("photonBarrelCosMax").set(0.8)
TauMuGamma.property("photonEndcapCosMin").set(0.86)
TauMuGamma.property("photonEndcapCosMax").set(0.92)
TauMuGamma.property("photonBarrelEmin").set(0.025)
TauMuGamma.property("photonEndcapEmin").set(0.05)
TauMuGamma.property("photonMergeAngleDeg").set(2.0)

# E/p thresholds (paper)
TauMuGamma.property("eoverpElectronMin").set(0.8)
TauMuGamma.property("eoverpMuonMax").set(0.5)

# pi0 mass window (paper)
TauMuGamma.property("pi0MassWinLow").set(0.12)
TauMuGamma.property("pi0MassWinHigh").set(0.14)

# signal muon/photon kinematics (paper)
TauMuGamma.property("pMuMin").set(0.4)
TauMuGamma.property("pMuMax").set(1.7)
TauMuGamma.property("EsigGammaMin").set(0.4)
TauMuGamma.property("EsigGammaMax").set(1.7)

# e-tag
TauMuGamma.property("cut_cosSigGamma_TagCharged").set(-0.2)
TauMuGamma.property("cut_pTagMin").set(0.5)
TauMuGamma.property("cut_EmissMax_e").set(1.7)

# pi-tag
TauMuGamma.property("cut_EmissMin_pi").set(0.7)
TauMuGamma.property("cut_absCosMissMax_pi").set(0.6)
TauMuGamma.property("cut_EsigGammaMin_pi").set(0.8)
TauMuGamma.property("cut_M2missMax_pi").set(0.050)

# rho-tag
TauMuGamma.property("cut_M2missMax_rho").set(0.075)
TauMuGamma.property("cut_cosHelMax_rho").set(0.8)
TauMuGamma.property("cut_absCosMissMax_rho").set(0.9)

AlgTask.show()
AlgTask.run()
```

## Submission/Sub_Ana/runANA.sh

```bash
#!/bin/bash
set -e

cd PATH
python NAME_ID.py > ../out/NAME_ID.out 2> ../out/NAME_ID.err
```

## Submission/Sub_Ana/SubJob_Ana.sh

```bash
#!/bin/bash
set -e
cur=$(pwd)

subjob="sub_taumugamma"

rec_root="/home/qukl/workerea/offline/TauMuGammaAna/Submission/Sub_Rec/sub_taumugamma/root"

root="${cur}/${subjob}/root"
script="${cur}/${subjob}/script"
out="${cur}/${subjob}/out"
mkdir -p "$root" "$script" "$out"

maxjob=100
maxevt=2000
jobname=joboptions_ana
ECMS=4.26

# escape & to avoid sed replacement issues
escape_sed () { echo "$1" | sed 's/[&]/\\&/g'; }

for k in $(seq 1 $maxjob); do
  idx=$((k-1))

  input="${rec_root}/joboptions_rec_${idx}.root"
  anaout="${root}/ana_taumugamma_${idx}.root"

  input_esc=$(escape_sed "$input")
  anaout_esc=$(escape_sed "$anaout")

  sed "s#INPUT#${input_esc}#g; \
       s#ANAOUT#${anaout_esc}#g; \
       s#MAXEVT#${maxevt}#g; \
       s#ECMS#${ECMS}#g" \
       "${cur}/${jobname}.py" > "${script}/${jobname}_${idx}.py"

  sed "s#PATH#${script}#g; \
       s#NAME#${jobname}#g; \
       s#ID#${idx}#g" \
       "${cur}/runANA.sh" > "${script}/run_${jobname}_${idx}.sh"

  chmod +x "${script}/run_${jobname}_${idx}.sh"
done

cd "$script"
hep_sub "run_${jobname}_%{ProcId}.sh" -g STCF -n ${maxjob} \
  -o "${out}/${jobname}_%{ProcId}.out" -e "${out}/${jobname}_%{ProcId}.err" \
  -wt long -mem 20000
cd "$cur"
```

## Submission/Sub_Rec/joboptions_rec.py

```python
#!/usr/bin/env python
#******************General_begin
import os
topdir = os.getenv('OFFLINETOP')

# MultiTask for background mixing (same as tau3mu)
from MultiStreamAnalysis import MultiStreamTask
AlgTask = MultiStreamTask("task")
AlgTask.property("TansSubColCache").set(True)
AlgTask.setLogLevel(4)  # higher is less

import RandomSvc
random_svc = AlgTask.createSvc("RandomSvc")
random_svc.property("Seed").set(RANDOM)           # template
random_svc.property("GenerateSeed").set(False)

import PodioDataSvc
dsvc = AlgTask.createSvc("PodioDataSvc")

# ---------------- Input (signal hits from full simulation) ----------------
import PodioSvc
AlgTask.addInputFile("INPUT", "Sig")  # template, same style as tau3mu :contentReference[oaicite:4]{index=4}

# ---------------- Optional background streams (kept as reference style) ---
BKGDir = os.getenv("OSCARDATA") + "/BkgDATA/offline_2.60_v7c3/ITKM_BTOF_20241225"
BKGNameList = ["Tous_e", "Tous_p", "Lumi", "Beamgas_e", "Beamgas_p"]
for n in BKGNameList:
    AlgTask.addInputFile(BKGDir + "/" + n + ".root", n)

# ---------------- Output (standard podio rec file) ----------------
outsvc = AlgTask.createSvc("PodioOutputSvc/OutputSvc")
outsvc.property("OutputFile").set("SDOUTPUT")  # template

# Keep a reduced list that contains what analysis needs (tracks/showers/PID/assembler)
reduced_list = [
    "TruthEventTimeCol",
    "MCParticleCol1",
    "CanTrackCol",
    "TrackerRecTrackCol",
    "TrackerHypoTrackCol",
    "DEDXRecTrackCol",
    "RecExtTrackCol",
    "RecECALHitCol",
    "RecECALShowerCol",
    "BTOFLikelihoodCol",
    "BTOFHypoTrackCol",
    "MUDTrackCol",
    "MUDClusterCol",
    "GlobalWLPidCol",
    "ReconstructedParticleCol",
    "DTOFCNNCol",
    "DTOFCNNHypoRecCol",
    "DTOFLikelihoodCol",
    "DTOFHypoTrackCol",
]
outsvc.property("OutputCollections").set(reduced_list)  # same idea as tau3mu :contentReference[oaicite:5]{index=5}

# ---------------- RootWriter (user-defined info; can be empty for now) ----
import RootWriter
osvc = AlgTask.createSvc("RootWriter/hSvc")
osvc.property("Output").set({"Fkey": "OUTDATA"})  # template :contentReference[oaicite:6]{index=6}

# ---------------- Geometry ----------------
import DDXMLSvc
xmlsvc = AlgTask.createSvc("DDXMLSvc")
xmlsvc.property("GeoCompactFileName").set(topdir + "/Geometry/FullGeometry/compact/STCF_Default.xml")

# ---------------- Background Mixer service (kept; mixing off by default) ---
import BKGMixerSvc
BKGMix = AlgTask.createSvc("BKGMixerSvc")
BKGMix.property("BkgStreamList").set(BKGNameList)
BKGMix.property("SubDetInfo").set(["ITK", "MDC", "PIDB", "PIDE", "ECAL", "MUD"])

BKGMix.property("IsAddPhysBkg").set(False)
BKGMix.property("IsAddBeamBkg").set(False)

BKGMix.property("EnforceFixTimeWindow").set(False)
BKGMix.property("AvePhyNumInEachColl").set(1.6e-3)
BKGMix.property("BeamBkgStrength").set(BEAMBKGSTRENGTH)  # template
BKGMix.property("FixSigTime").set(-999)
BKGMix.property("SigTimeInterval").set(-1)

import EventSampling
EvtSample = AlgTask.createAlg("EventSamplingAlg/EventSamplingAlg1")
EvtSample.property("T0Resolution").set(0.040)  # 40 ps
#******************General_end


###################################### Digitization ###########################################

#******************Tracker_begin
import ItkmGeomSvc
Itkmgeom1 = AlgTask.createSvc("ItkmGeomSvc")
Itkmgeom1.property("FromXML").set(True)

import MdcGeomSvc
Mdcgeom1 = AlgTask.createSvc("MdcGeomSvc")
Mdcgeom1.property("geomDataFilePath").set(topdir + "/CommonSvc/MDCSvc/MdcGeomSvc/dat/mdcGeomData.dat")
Mdcgeom1.property("geomDataSource").set(1)

import MdcCalibSvc
mdccalib = AlgTask.createSvc("MdcCalibSvc")

import TrackInfoSvc
trInfo = AlgTask.createSvc("TrackInfoSvc")

import ITKMDigitization
ITKMDigi = AlgTask.createAlg("ITKMDigiAlg")
ITKMDigi.property("Mode").set(11)
ITKMDigi.property("SamplingMethod").set(3)
ITKMDigi.property("SimpleDigitization").set(False)
ITKMDigi.property("TimeWindow").set(1000)
ITKMDigi.property("Threshold").set(300)
ITKMDigi.property("UseClock").set(False)
ITKMDigi.property("CombineNearbyHits").set(True)
ITKMDigi.property("SaveWaveform").set(False)
ITKMDigi.property("mobility_model").set("masetti_canali")
ITKMDigi.property("recombination_model").set("srh_auger")
ITKMDigi.property("multiplication_model").set("none")
ITKMDigi.property("trapping_model").set("none")
ITKMDigi.property("detrapping_model").set("none")

import MDCDigitization
MDCDigit = AlgTask.createAlg("MDCDigiAlg/MDCDigiAlg")
MDCDigit.property("Mode").set(0)
MDCDigit.property("TimeStep").set(5)
MDCDigit.property("UseFFT").set(True)
MDCDigit.property("CalPropTime").set(True)
MDCDigit.property("SaveWaveHist").set(False)
#******************Tracker_end


#******************ECAL_begin
import ECALParaMgr
ECALPar = AlgTask.createSvc("ECALParaMgr")
ECALPar.property("ECALBarrelName").set("ECAL_Barrel")
ECALPar.property("ECALEndcapName").set("ECAL_Endcap")

import ECALHitProcSvc
ECALHitp = AlgTask.createSvc("ECALHitProcSvc")
ECALHitp.property("fMethod").set(4)

import ECALDigitization
ECALDigit = AlgTask.createAlg("ECALDigiAlg")
#******************ECAL_end


#******************DTOF_begin
import DTOFGeoSvc
DTOFgeom = AlgTask.createSvc("DTOFGeoSvc")

import DTOFDigitization
DTOFDigit = AlgTask.createAlg("DTOFDigiAlg")
DTOFDigit.property("Mode").set(1)
#******************DTOF_end


#******************BTOF_begin
import BTOFGeoSvc
BTOFgeom = AlgTask.createSvc("BTOFGeoSvc")

import BTOFDigitization
BTOFDigit = AlgTask.createAlg("BTOFDigiAlg")
BTOFDigit.property("Mode").set(1)
#******************BTOF_end


#******************MUD_begin
import MUDIDSvc
mudidsvc = AlgTask.createSvc("MUDIDSvc")

import MUDGeoSvc
mudgeosvc = AlgTask.createSvc("MUDGeoSvc")

import MUDDigitization
MUDDigi = AlgTask.createAlg("MUDDigiAlg")
MUDDigi.property("MUDDigiDataPath").set(topdir + "/Digitization/MUDDigitization/share/")
#******************MUD_end


###################################### Reconstruction ###########################################

#******************Tracker_begin
import ITKMClusterRec
ITKMRec = AlgTask.createAlg("ITKMClusterRecAlg")
ITKMRec.property("Mode").set(0)
ITKMRec.property("ToT_Threshold").set(0)
ITKMRec.property("CalibrateBias").set(False)
ITKMRec.property("UseDeltaTCut").set(True)

import EventStartTimeSvc
ESTSvc = AlgTask.createSvc("EventStartTimeSvc")

import MDCHitRec
mdcHitrec = AlgTask.createAlg("MDCHitRecAlg/MDCHitRecAlg1")
mdcHitrec.property("SelectMode").set(11)

import HoughFinder
houghfind = AlgTask.createAlg("HoughFinderAlg/HoughFinderAlg1")
houghfind.property("MapHitPreFilter").set(True)
houghfind.property("SecTrackFinder").set(True)
houghfind.property("UseStereoLayerForPeaking").set(False)

import FitTrackAlg
FitAlg1 = AlgTask.createAlg("FitTrackAlg/FitTrackAlg1")

import DeDxRec
dedxrecsvc = AlgTask.createSvc("DeDxRecSvc")
dedxrecalg = AlgTask.createAlg("DeDxRecAlg")
dedxrecalg.property("SaveITKMdEdx").set(True)

import TrackExtAlg
TrkExtAlg1 = AlgTask.createAlg("TrackExtAlg/TrackExtAlg1")
#******************Tracker_end


#******************ECAL_begin
import ECALRecAlg
ECALRec = AlgTask.createAlg("ECALRec/ECALRecAlg")
ECALRec.property("RecMode").set(2)
#******************ECAL_end


#******************DTOF_begin
import DTOFRecAlg
DTOFRec = AlgTask.createAlg("DTOFRecAlg")
DTOFRec.property("DigiMode").set(1)
DTOFRec.property("RecMode").set(0)
DTOFRec.property("SplitMom").set(1800)
DTOFRec.property("SavePDF").set(True)
#******************DTOF_end


#******************BTOF_begin
import BTOFRecAlg
BTOFRec = AlgTask.createAlg("BTOFRecAlg")
BTOFRec.property("DigiMode").set(1)
BTOFRec.property("RecMode").set(1)
BTOFRec.property("SplitMom").set(1800)
BTOFRec.property("Tmax").set(40)
BTOFRec.property("SavePDF").set(True)
#******************BTOF_end


#******************MUD_begin
import MUDRecAlg
MUDRec = AlgTask.createAlg("MUDRecIdentifyAlg")
MUDRec.property("UseTrackRec").set(True)
MUDRec.property("UseClusterRec").set(True)

HitTool = MUDRec.createTool("MUDRecHitTool")
HitTool.property("UseDigi").set(True)

TCFitTool = MUDRec.createTool("MUDRecTCFitTool")
TCFitTool.property("TrackRecVer").set(2)
TCFitTool.property("ClusterRecVer").set(1)

BDTTool = MUDRec.createTool("MUDRecBDTTool")
BDTTool.property("TrackBDTWeightsList").set(os.getenv("OSCARDATA") + "/MUD/MUDBDTData/20250103/TrackWeights/WeightsList.dat")
BDTTool.property("ClusterBDTWeightsList").set(os.getenv("OSCARDATA") + "/MUD/MUDBDTData/20250103/ClusterWeights/WeightsList.dat")
BDTTool.property("TrackBDTDataList").set(os.getenv("OSCARDATA") + "/MUD/MUDBDTData/20250103/TrackData/DataList.dat")
BDTTool.property("ClusterBDTDataList").set(os.getenv("OSCARDATA") + "/MUD/MUDBDTData/20250103/ClusterData/DataList.dat")
BDTTool.property("ReRunTrackBDTCut").set(False)
BDTTool.property("ReRunClusterBDTCut").set(False)
BDTTool.property("TrackSuppression").set(30)
BDTTool.property("ClusterSuppression").set(30)
#******************MUD_end


#******************Classic Global PID (produces GlobalWLPidCol etc.)
import GlobalWLPID
wlpid = AlgTask.createAlg("GlobalWLPID")
wlpid.property("Pprior").set([1, 1, 1, 1, 1])

#******************EventAssembler (produces ReconstructedParticleCol)
import EventAssembler
EventAs = AlgTask.createAlg("EventAssembler")

#******************SummaryWriter (optional; keep as reference)
import TrackSummaryWriter
recheck = AlgTask.createAlg("TrackSummaryWriter")
recheck.property("TrackingSummary").set(True)
recheck.property("ExtTrackRes").set(True)
recheck.property("ITKWRecCheck").set(False)
recheck.property("ITKMRecCheck").set(True)

import EMCSummaryWriter
EMCrecheck = AlgTask.createAlg("EMCSummaryWriter")

import DTOFSummaryWriter
DTOFrecheck = AlgTask.createAlg("DTOFSummaryWriter")

import BTOFSummaryWriter
alg = AlgTask.createAlg("BTOFSummaryWriter")

import MUDSummaryWriter
MUDSummaryWriterAlg = AlgTask.createAlg("MUDSummaryWriter")
#******************SummaryWriter END

AlgTask.setEvtMax(MAXEVT)  # template
AlgTask.show()
AlgTask.run()
```

## Submission/Sub_Rec/runREC.sh

```bash
#!/bin/bash

# templates in this file:
# PATH, NAME, ID
# will be replaced by SubJob_Rec.sh via sed.

scriptpath=PATH
cd $scriptpath

python NAME_ID.py
```

## Submission/Sub_Rec/SubJob_Rec.sh

```bash
#!/bin/bash
# SubJob_Rec.sh
# Full digitization + reconstruction for tau -> mu gamma full-sim samples
# NOTE: please `source hepjob` before running.

set -e
cur=$(pwd)

# If you want to test beam background mixing later, change it here.
BKGSTRENGTH=1

# Output subjob name (reco outputs will be under Sub_Rec/sub_taumugamma/)
subJobList=("sub_taumugamma")

# IMPORTANT: point to your simulation output directory (Sub_Sim)
dataList=("/home/qukl/workerea/offline/TauMuGammaAna/Submission/Sub_Sim/sub_taumugamma")

for i in "${!subJobList[@]}"; do
  subjob=${subJobList[$i]}
  datadir=${dataList[$i]}

  root=${cur}/${subjob}/root
  script=${cur}/${subjob}/script
  out=${cur}/${subjob}/out
  mkdir -p "$root" "$script" "$out"

  maxjob=100     # job number (should match your sim production)
  maxevt=2000    # events per job
  jobname=joboptions_rec

  for k in $(seq 1 $maxjob); do
    idx=$((k-1))

    input=${datadir}/root/joboptions_sim_${idx}.root
    output=${root}/${jobname}_${idx}.root
    outputRoot=${root}/TrackingInfo_${idx}.root
    random=123456$k

    # 1) generate joboptions python from template
    sed "s#RANDOM#${random}#g; \
         s#BEAMBKGSTRENGTH#${BKGSTRENGTH}#g; \
         s#INPUT#${input}#g; \
         s#SDOUTPUT#${output}#g; \
         s#OUTDATA#${outputRoot}#g; \
         s#MAXEVT#${maxevt}#g" \
         ${cur}/${jobname}.py > ${script}/${jobname}_${idx}.py

    # 2) generate run shell from template
    sed "s#PATH#${script}#g; \
         s#NAME#${jobname}#g; \
         s#ID#${idx}#g" \
         ${cur}/runREC.sh > ${script}/run_${jobname}_${idx}.sh

    chmod +x ${script}/run_${jobname}_${idx}.sh
  done

  cd "$script"
  echo "[INFO] submit ${jobname} ..."
  hep_sub "run_${jobname}_%{ProcId}.sh" \
    -g STCF \
    -n ${maxjob} \
    -o ${out}/${jobname}_%{ProcId}.out \
    -e ${out}/${jobname}_%{ProcId}.err \
    -wt long \
    -mem 20000
  echo "[INFO] submit finished."
done

cd "${cur}"
```

## Submission/Sub_Sim/decay.dec

```text
#
# taumugamma.dec
# Signal:  e+e- -> tau+ tau- (from KKMC) ; decay via EvtGen user decay table
#   - one tau forced to mu gamma (PHSP)  [paper: PHSP signal]
#   - the other tau forced into tag modes: e nu nu / pi nu / pi pi0 nu (rho)
#

Particle vpho 4.26 0.0

Decay vpho
  1.000 tau+ tau- PHSP;
Enddecay

# ---- Signal side: force tau- -> mu- gamma (PHSP)
# (If you want tau+ as signal instead, swap tau+ and tau- blocks.)
Decay tau-
  1.0000  mu-  gamma  PHSP;
Enddecay

# ---- Tag side: only keep three tag modes (e / pi / rho)
# Use relative weights (sum to 1). You may tune these to PDG values if desired.
Decay tau+
  
  # Mode 1: Electronic decay (V-A interaction)
  # 使用 TAULNUNU 模型代替 PHSP
  0.33  e+   nu_e      anti-nu_tau   TAULNUNU;

  # Mode 2: Pion decay (Scalar meson)
  # 使用 TAUSCALARNU 模型
  0.20  pi+  anti-nu_tau             TAUSCALARNU;

  # Mode 3: Rho decay (Vector meson -> pi pi0)
  # 使用 TAUVECTORNU 模型，它会自动处理 rho 的共振态结构
  # 这里的 rho+ 会进一步衰变为 pi+ pi0
  0.47  rho+ anti-nu_tau             TAUVECTORNU; 

Enddecay

# pi0 -> gamma gamma is defined in standard DECAY.DEC; no need to redefine.

End
```

## Submission/Sub_Sim/joboptions_sim.py

```python
#!/usr/bin/env python
# -*- coding:utf-8 -*-

import Sniper
import os

task = Sniper.Task("task")
task.setLogLevel(4)  # higher is less

topdir = os.getenv('OFFLINETOP')

# -----------------------
# Random
# -----------------------
import RandomSvc
random_svc = task.createSvc("RandomSvc")
random_svc.property("Seed").set(RANDOM)
random_svc.property("GenerateSeed").set(False)

# =======================
# Generator: KKMC -> tau+ tau-
# =======================
import KKMC
kkmc = task.createAlg("KKMC")
kkmc.property("CMSEnergy").set(CME)
kkmc.property("BeamEnergySpread").set(0.0013)
kkmc.property("NumberOfEventPrinted").set(0)

# Only generate tau pair
kkmc.property("GenerateDownQuark").set(False)
kkmc.property("GenerateUpQuark").set(False)
kkmc.property("GenerateStrangeQuark").set(False)
kkmc.property("GenerateCharmQuark").set(False)
kkmc.property("GenerateMuonPair").set(False)
kkmc.property("GenerateTauPair").set(True)

kkmc.property("ParticleDecayThroughEvtGen").set(True)
kkmc.property("RadiationCorrection").set(True)

# -----------------------
# EvtGen decay with user decay table (PHSP for tau->mu gamma)
# -----------------------
import StcfEvtGen
evt = task.createAlg("EvtDecay")
# Use default DECAY.DEC / pdt.table from $STCFEVTGENROOT, but override user decay table
evt.property("DecayDecDir").set("")          # "" means default: $STCFEVTGENROOT/share/DECAY.DEC
evt.property("PdtTableDir").set("")          # "" means default: $STCFEVTGENROOT/share/pdt.table
evt.property("userDecayTableName").set("DECAYFILE")

# =======================
# Data service & Output
# =======================
import PodioDataSvc
dsvc = task.createSvc("PodioDataSvc")

import PodioSvc
oSvc = task.createSvc("PodioOutputSvc/OutputSvc")
oSvc.property("OutputFile").set("OUTPUT")

# 完全照抄 tau3mu 的 reduced 输出（关键！）
reduced_coll_list = [
  "EventHeaderCol","TruthEventTimeCol","MCParticleCol",
  "ITKHitCol","MDCHitCol",
  "DTOFBarHitCol","DTOFPDHitCol",
  "BTOFBarHitCol","BTOFPDHitCol",
  "ECALHit3Col","ECALTPointCol","MUDPointCol"
]
oSvc.property("OutputCollections").set(reduced_coll_list)


# !!! IMPORTANT !!!
# 你的参考 joboptions 里 OutputCollections 被省略成了 "...".
# 为保证可用性，这里不强制裁剪集合（输出全量/默认集合），适合 full-sim -> 后续 digit/reco。
# 如果你们组有固定的 reduced_coll_list，请把那份不带省略号的列表贴过来，我可以一键对齐。
# 若你的版本必须显式设置 OutputCollections，可先注释掉下一行运行验证。
# oSvc.property("OutputCollections").set([...])

# -----------------------
# Geometry
# -----------------------
import GeometrySvc
myxmlsvc = task.createSvc("GeometrySvc")
myxmlsvc.property("GeoCompactFileName").set(topdir + "/Geometry/FullGeometry/compact/STCF_Default.xml")
myxmlsvc.property("LowLooperThreshold").set(False)
myxmlsvc.property("TightFieldParameters").set(False)

# Sub-detector geometry services (same style as your tau3mu template)
import MdcGeomSvc
mdcgeom = task.createSvc("MdcGeomSvc")
mdcgeom.property("geomDataFilePath").set(topdir + "/CommonSvc/MDCSvc/MdcGeomSvc/dat/mdcGeomData.dat")
mdcgeom.property("geomDataSource").set(0)

import ItkmGeomSvc
itkmgeom = task.createSvc("ItkmGeomSvc")
itkmgeom.property("FromXML").set(True)

import BTOFGeoSvc
btofgeom = task.createSvc("BTOFGeoSvc")

import DTOFGeoSvc
dtofgeom = task.createSvc("DTOFGeoSvc")

# -----------------------
# G4
# -----------------------
import G4Svc
g4svc = task.createSvc("G4Svc")
# g4svc.property("VisMac").set("vis.mac")
# g4svc.property("RunMac").set("run.mac")

# -----------------------
# Full detector simulation
# -----------------------
import DetSimAlg
simalg = task.createAlg("DetSimAlg/DetSimAlg")
simalg.property("DetFactory").set("FullFactory")

import FullSim
genTool = simalg.createTool("GeneratorMgr")
genTool.property("Translation").set([0, 0, 0])      # [mm]
genTool.property("ParticleSource").set("Generator") # Generator or ParticleGun
genTool.property("Boost").set(True)
genTool.property("XTSampling").set(True)

factory = task.createSvc("FullSimFactory/FullFactory")
factory.property("G4Verbose").set(0)
factory.property("OpticalSimOn").set(True)
factory.property("PAISimOn").set(False)

factory.property("ITKPAISimOn").set(True)
factory.property("ITKMStepLimiter").set(True)
factory.property("ITKMMaximumStep").set(10.0)

# AnaMgrList：你参考脚本里也被 "..." 省略；这里给一个保守可用的列表
# 若你们版本要求必须列全，请把你们 FullFactory 的默认 AnaMgrList 配置贴出来，我再对齐。
factory.property("AnaMgrList").set([
    "GeoAnaMgr", "GeneratorMgr",
    "MDCAnaMgr", "ECALAnaMgr",
    "DTOFAnaMgr", "BTOFAnaMgr",
    "MUDAnaMgr"
])

# optional: MUD PID service (same as tau3mu template)
import MUDIDSvc
mudsvc = task.createSvc("MUDIDSvc")

# -----------------------
# Run
# -----------------------
task.setEvtMax(MAXEVT)
task.show()
task.run()
```

## Submission/Sub_Sim/runSIM.sh

```bash
#!/bin/bash

scriptpath=PATH
cd $scriptpath

python NAME_ID.py
```

## Submission/Sub_Sim/SubJob_Sim.sh

```bash
#!/bin/bash
# SubJob_Sim_TauMuGamma.sh
# FullSim signal sample production for tau -> mu gamma
# NOTE: please `source hepjob` before running.

set -e

cur=$(pwd)

subJobList=("sub_taumugamma")

for subjob in "${subJobList[@]}"; do
  root=${cur}/${subjob}/root
  script=${cur}/${subjob}/script
  out=${cur}/${subjob}/out
  mkdir -p "$root" "$script" "$out"

  maxjob=100        # number of jobs
  maxevt=2000       # events per job
  jobname=joboptions_sim
  decayFile=${cur}/decay.dec
  cme=4.26

  for k in $(seq 1 $maxjob); do
    output=${root}/${jobname}_$((k-1)).root
    random=123456$k

    # 1) produce joboptions python
    sed "s#RANDOM#${random}#g; \
         s#CME#${cme}#g; \
         s#DECAYFILE#${decayFile}#g; \
         s#OUTPUT#${output}#g; \
         s#MAXEVT#${maxevt}#g;" \
         ${cur}/${jobname}.py > ${script}/${jobname}_$((k-1)).py

    # 2) produce run shell script
    sed "s#PATH#${script}#g; \
         s#NAME#${jobname}#g; \
         s#ID#$((k-1))#g;" \
         ${cur}/runSIM.sh > ${script}/run_${jobname}_$((k-1)).sh

    chmod +x ${script}/run_${jobname}_$((k-1)).sh
  done

  cd "$script"
  echo "[INFO] submit ${jobname} ..."

  # Adjust hep_sub options if your site requires different flags
  hep_sub "run_${jobname}_%{ProcId}.sh" \
    -g STCF \
    -n ${maxjob} \
    -o ${out}/${jobname}_%{ProcId}.out \
    -e ${out}/${jobname}_%{ProcId}.err \
    -wt long

  echo "[INFO] submit finished."
done

cd "${cur}"
```

## TauMuGammaAna/CMakeLists.txt

```cmake
# Analysis/TauMuGammaAna/CMakeLists.txt

# 与 Tau3mu 一致的最小做法：只声明你确定在 build/lib 中存在的库
# （否则会触发 /usr/bin/ld: cannot find -lXXX）
find_package(podio REQUIRED)
find_package(ROOT REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# 生成 libTauMuGammaAna.so
# 注意：不要写 SniperKernel / RootWriter，否则会强制 -lSniperKernel/-lRootWriter 并导致你现在的链接失败
oscar_add_pkg(TauMuGammaAna
  DEPENDS
    PodioDataSvc
    DataModel
    GlobalPID
    VertexFit
)

# 安装 python 模块（用于 joboptions: import TauMuGammaAna）
install(DIRECTORY python/ DESTINATION ${CMAKE_INSTALL_PREFIX}/python)
```

## TauMuGammaAna/include/TauMuGammaAlg.h

```cpp
#ifndef TAUMUGAMMAANA_TAUMUGAMMAALG_H
#define TAUMUGAMMAANA_TAUMUGAMMAALG_H

#include "SniperKernel/AlgBase.h"
#include <string>

class TTree;
class RootWriter;
class GlobalPIDSvc;

class TauMuGammaAlg : public AlgBase {
public:
    TauMuGammaAlg(const std::string& name);
    virtual ~TauMuGammaAlg();

    virtual bool initialize() override;
    virtual bool execute() override;
    virtual bool finalize() override;

private:
    void resetBranches();

private:
    // -------- Services --------
    RootWriter*   m_rootWriter {nullptr};
    GlobalPIDSvc* m_gPIDSvc    {nullptr};
    TTree*        m_tree       {nullptr};

    // -------- Input collection names (make them configurable) --------
    std::string m_recParColName {"ReconstructedParticleCol"};
    std::string m_showerColName {"RecECALShowerCol"};
    std::string m_treePath      {"Fkey/TauMuGamma"};

    // -------- Paper-driven properties --------
    double m_ecms {4.26};

    int  m_analysisMode {0};               // 0=PaperStrict, 1=FullSimRobust
    bool m_requireExactPhotonMultiplicity {true};
    bool m_useEoverPTagID {true};

    double m_pMu {0.0};        // |p_mu|
    double m_cosMuGamma {0.0}; // cos(angle(mu, gamma))


    // acceptance & basic
    double m_cosThetaAccCut {0.93};
    double m_cosThetaMuGammaCut {-0.35};
    double m_chargedPmin {0.2}; // GeV, clean extra soft tracks
    double m_muPidMin {0.2};    // min mu PID prob
    double m_muPidSlack {0.05}; // allow small PID overlap

    // photon acceptance & thresholds (paper)
    double m_photonBarrelCosMax {0.8};
    double m_photonEndcapCosMin {0.86};
    double m_photonEndcapCosMax {0.92};
    double m_photonBarrelEmin {0.025};
    double m_photonEndcapEmin {0.05};
    double m_photonMergeAngleDeg {2.0}; // merge split clusters within this angle

    // E/p tag ID thresholds
    double m_eoverpElectronMin {0.8};
    double m_eoverpMuonMax {0.5};

    // rho tag: pi0 mass window (paper)
    double m_pi0MassWinLow {0.12};
    double m_pi0MassWinHigh {0.14};

    // signal muon/photon kinematics (paper)
    double m_pMuMin {0.4};
    double m_pMuMax {1.7};
    double m_EsigGammaMin {0.4};
    double m_EsigGammaMax {1.7};

    // -------- Table-1 cuts (locked by joboptions_ana.py anyway) --------
    // e tag
    double m_cut_cosSigGamma_TagCharged {-0.2};
    double m_cut_pTagMin {0.5};
    double m_cut_EmissMax_e {1.7};

    // pi tag
    double m_cut_EmissMin_pi {0.7};
    double m_cut_absCosMissMax_pi {0.6};
    double m_cut_EsigGammaMin_pi {0.8};
    double m_cut_M2missMax_pi {0.050};

    // rho tag
    double m_cut_M2missMax_rho {0.075};
    double m_cut_cosHelMax_rho {0.8};
    double m_cut_absCosMissMax_rho {0.9};

    // -------- Ntuple variables --------
    int    m_nCharged {0};
    int    m_nGamma {0};
    int    m_nExtraGamma {0};
    int    m_tagMode {0};         // 1=e, 2=pi, 3=rho

    int    m_passBasic {0};
    int    m_passFinal {0};

    double m_invMass_Signal {0.0}; // M(mu,gamma)
    double m_E_mugamma {0.0};
    double m_E_sys {0.0};
    double m_P_sys {0.0};
    double m_Emiss {0.0};
    double m_E_miss {0.0};
    double m_P_miss {0.0};
    double m_M2miss {0.0};
    double m_cosTheta_miss {0.0};
    int    m_check_miss {0};

    double m_cosSigGamma_TagCharged {0.0};
    double m_pTag {0.0};
    double m_EsigGamma {0.0};
    double m_cosHel {0.0};

    // -------- counters --------
    long long m_evtAll {0};
    long long m_evtPassBasic {0};
    long long m_evtPassFinal {0};

    // -------- cutflow counters --------
    long long m_fail_noRecOrShower {0};
    long long m_fail_nCharged {0};
    long long m_fail_qsum {0};
    long long m_fail_nGamma {0};

    long long m_fail_pid_noOppCharge {0};
    long long m_fail_pid_noMuonLike {0};
    long long m_fail_pid_noMuonEop {0};
    long long m_fail_pid_noTagLike {0};
    long long m_fail_pid_noTagEop {0};
    long long m_fail_pid_noPair {0};

    long long m_fail_rho_noSigG {0};
    long long m_fail_sigG_exact {0};
    long long m_fail_sigG_none {0};

    long long m_fail_kin_cosMuGamma {0};
    long long m_fail_kin_pMu {0};
    long long m_fail_kin_EsigGamma {0};

    long long m_fail_e_cosSigGamma {0};
    long long m_fail_e_pTag {0};
    long long m_fail_e_Emiss {0};

    long long m_fail_pi_Emiss {0};
    long long m_fail_pi_cosMiss {0};
    long long m_fail_pi_EsigGamma {0};
    long long m_fail_pi_M2miss {0};

    long long m_fail_rho_M2miss {0};
    long long m_fail_rho_cosHel {0};
    long long m_fail_rho_cosMiss {0};

    long long m_fail_tagModeUnknown {0};
};

#endif
```

## TauMuGammaAna/src/TauMuGammaAlg.cc

```cpp
#include "TauMuGammaAlg.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "SniperKernel/AlgFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperPtr.h"

#include "RootWriter/RootWriter.h"
#include "GlobalPID/GlobalPIDSvc.h"

#include "DataModel/ReconstructedParticle.h"
#include "DataModel/ReconstructedParticleCollection.h"
#include "DataModel/RecECALShower.h"
#include "DataModel/RecECALShowerCollection.h"
#include "podio/ObjectID.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"

DECLARE_ALGORITHM(TauMuGammaAlg);

namespace {
    const double kMassE   = 0.00051099895;
    const double kMassMu  = 0.1056583745;
    const double kMassPi  = 0.13957039;
    const double kMassPi0 = 0.1349768;
    const double kMassTau = 1.77686;
    const double kSigmaM  = 0.06;
    const double kSigmaE  = 0.2;

    inline double cosAngle(const TVector3& a, const TVector3& b){
        double da=a.Mag(), db=b.Mag();
        if (da<=0 || db<=0) return 1.0;
        return a.Dot(b)/(da*db);
    }

    inline bool passAccCos(const TVector3& p, double cutAbsCos){
        double pp = p.Mag();
        if (pp<=0) return false;
        return std::fabs(p.Z()/pp) < cutAbsCos;
    }

    inline double calcEoverP(const ReconstructedParticle& rp){
        auto p = rp.getMomentum();
        TVector3 v(p.x, p.y, p.z);
        double pmag = v.Mag();
        if (pmag<=0) return -1.0;

        double e = 0.0;
        if (rp.getShower().isAvailable()) {
            const auto& sh = rp.getShower();
            e = sh.getEnergy();
        } else {
            e = rp.getEnergy();
        }
        if (e<=0) return -1.0;
        return e/pmag;
    }

    TLorentzVector makePhotonP4(const RecECALShower& sh){
        auto pos = sh.getPosition();
        TVector3 dir(pos.x, pos.y, pos.z);
        if (dir.Mag()<=0) dir.SetXYZ(0,0,1);
        dir = dir.Unit();
        double E = sh.getEnergy();
        return TLorentzVector(E*dir.X(), E*dir.Y(), E*dir.Z(), E);
    }

    TLorentzVector makeChargedP4(const ReconstructedParticle& rp, double mass){
        auto p = rp.getMomentum();
        TVector3 v(p.x, p.y, p.z);
        double E = std::sqrt(v.Mag2() + mass*mass);
        return TLorentzVector(v.X(), v.Y(), v.Z(), E);
    }

    inline int pickBestSignalPhoton(const RecECALShowerCollection* showers,
                                    const std::vector<int>& idxs,
                                    const TLorentzVector& p4_mu,
                                    double targetEtau){
        int best = -1;
        double bestScore = std::numeric_limits<double>::max();
        for (int id : idxs) {
            TLorentzVector p4g = makePhotonP4(showers->at(id));
            TLorentzVector p4mug = p4_mu + p4g;
            double dM = std::fabs(p4mug.M() - kMassTau);
            double dE = std::fabs(p4mug.E() - targetEtau);
            double score = (dM/kSigmaM)*(dM/kSigmaM) + (dE/kSigmaE)*(dE/kSigmaE);
            if (score < bestScore) {
                bestScore = score;
                best = id;
            }
        }
        return best;
    }

    inline double signalScore(const RecECALShowerCollection* showers,
                              int id,
                              const TLorentzVector& p4_mu,
                              double targetEtau){
        TLorentzVector p4g = makePhotonP4(showers->at(id));
        TLorentzVector p4mug = p4_mu + p4g;
        double dM = std::fabs(p4mug.M() - kMassTau);
        double dE = std::fabs(p4mug.E() - targetEtau);
        return (dM/kSigmaM)*(dM/kSigmaM) + (dE/kSigmaE)*(dE/kSigmaE);
    }
}

TauMuGammaAlg::TauMuGammaAlg(const std::string& name) : AlgBase(name) {
    // collection names & output path
    declProp("recParColName", m_recParColName);
    declProp("showerColName", m_showerColName);
    declProp("treePath",      m_treePath);

    // basic
    declProp("cmsEnergy", m_ecms);
    declProp("analysisMode", m_analysisMode);
    declProp("requireExactPhotonMultiplicity", m_requireExactPhotonMultiplicity);
    declProp("useEoverPTagID", m_useEoverPTagID);
    declProp("cosThetaAccCut", m_cosThetaAccCut);
    declProp("cosThetaMuGammaCut", m_cosThetaMuGammaCut);
    declProp("chargedPmin", m_chargedPmin);
    declProp("muPidMin", m_muPidMin);
    declProp("muPidSlack", m_muPidSlack);

    // photon acceptance & thresholds
    declProp("photonBarrelCosMax", m_photonBarrelCosMax);
    declProp("photonEndcapCosMin", m_photonEndcapCosMin);
    declProp("photonEndcapCosMax", m_photonEndcapCosMax);
    declProp("photonBarrelEmin", m_photonBarrelEmin);
    declProp("photonEndcapEmin", m_photonEndcapEmin);
    declProp("photonMergeAngleDeg", m_photonMergeAngleDeg);

    // E/p tag ID
    declProp("eoverpElectronMin", m_eoverpElectronMin);
    declProp("eoverpMuonMax", m_eoverpMuonMax);

    // rho tag: pi0 mass window
    declProp("pi0MassWinLow",  m_pi0MassWinLow);
    declProp("pi0MassWinHigh", m_pi0MassWinHigh);

    // signal muon/photon kinematics
    declProp("pMuMin", m_pMuMin);
    declProp("pMuMax", m_pMuMax);
    declProp("EsigGammaMin", m_EsigGammaMin);
    declProp("EsigGammaMax", m_EsigGammaMax);

    // Table-1 cuts
    declProp("cut_cosSigGamma_TagCharged", m_cut_cosSigGamma_TagCharged);
    declProp("cut_pTagMin", m_cut_pTagMin);
    declProp("cut_EmissMax_e", m_cut_EmissMax_e);

    declProp("cut_EmissMin_pi", m_cut_EmissMin_pi);
    declProp("cut_absCosMissMax_pi", m_cut_absCosMissMax_pi);
    declProp("cut_EsigGammaMin_pi", m_cut_EsigGammaMin_pi);
    declProp("cut_M2missMax_pi", m_cut_M2missMax_pi);

    declProp("cut_M2missMax_rho", m_cut_M2missMax_rho);
    declProp("cut_cosHelMax_rho", m_cut_cosHelMax_rho);
    declProp("cut_absCosMissMax_rho", m_cut_absCosMissMax_rho);
}

TauMuGammaAlg::~TauMuGammaAlg() {}

void TauMuGammaAlg::resetBranches(){
    m_nCharged=0; m_nGamma=0; m_nExtraGamma=0; m_tagMode=0;
    m_passBasic=0; m_passFinal=0;
    m_pMu=0; m_cosMuGamma=0;
    m_invMass_Signal=0; m_E_mugamma=0;
    m_E_sys=0; m_P_sys=0;
    m_Emiss=0; m_E_miss=0; m_P_miss=0; m_M2miss=0; m_cosTheta_miss=0; m_check_miss=0;
    m_cosSigGamma_TagCharged=0;
    m_pTag=0; m_EsigGamma=0; m_cosHel=0;
}

bool TauMuGammaAlg::initialize() {
    // GlobalPID
    SniperPtr<GlobalPIDSvc> pid(getParent(), "GlobalPIDSvc");
    if (!pid.valid()) {
        LogError << "Failed to get GlobalPIDSvc" << std::endl;
        return false;
    }
    m_gPIDSvc = pid.data();

    // RootWriter
    SniperPtr<RootWriter> rw(getParent(), "hSvc");
    if (!rw.valid()) {
        LogError << "Failed to get RootWriter (hSvc)" << std::endl;
        return false;
    }
    m_rootWriter = rw.data();

    // book tree
    m_tree = m_rootWriter->bookTree(*getParent(), m_treePath, "TauMuGamma");
    if (!m_tree) {
        LogError << "Failed to book tree: " << m_treePath << std::endl;
        return false;
    }

    m_tree->Branch("pMu", &m_pMu, "pMu/D");
    m_tree->Branch("cosMuGamma", &m_cosMuGamma, "cosMuGamma/D");

    m_tree->Branch("nCharged", &m_nCharged, "nCharged/I");
    m_tree->Branch("nGamma",   &m_nGamma,   "nGamma/I");
    m_tree->Branch("nExtraGamma", &m_nExtraGamma, "nExtraGamma/I");
    m_tree->Branch("tagMode",  &m_tagMode,  "tagMode/I");
    m_tree->Branch("passBasic", &m_passBasic, "passBasic/I");
    m_tree->Branch("passFinal", &m_passFinal, "passFinal/I");

    m_tree->Branch("M_mugamma", &m_invMass_Signal, "M_mugamma/D");
    m_tree->Branch("E_mugamma", &m_E_mugamma, "E_mugamma/D");
    m_tree->Branch("E_sys", &m_E_sys, "E_sys/D");
    m_tree->Branch("P_sys", &m_P_sys, "P_sys/D");
    m_tree->Branch("Emiss",     &m_Emiss,          "Emiss/D");
    m_tree->Branch("E_miss",    &m_E_miss,         "E_miss/D");
    m_tree->Branch("P_miss",    &m_P_miss,         "P_miss/D");
    m_tree->Branch("check_miss", &m_check_miss,    "check_miss/I");
    m_tree->Branch("M2miss",    &m_M2miss,         "M2miss/D");
    m_tree->Branch("cosTheta_miss", &m_cosTheta_miss, "cosTheta_miss/D");

    m_tree->Branch("cosSigGamma_TagCharged", &m_cosSigGamma_TagCharged, "cosSigGamma_TagCharged/D");
    m_tree->Branch("pTag", &m_pTag, "pTag/D");
    m_tree->Branch("EsigGamma", &m_EsigGamma, "EsigGamma/D");
    m_tree->Branch("cosHel", &m_cosHel, "cosHel/D");

    return true;
}

bool TauMuGammaAlg::execute() {
    ++m_evtAll;
    resetBranches();

    auto recPars = getROColl(ReconstructedParticleCollection, m_recParColName);
    auto showers = getROColl(RecECALShowerCollection, m_showerColName);
    if (!recPars || !showers) {
        ++m_fail_noRecOrShower;
        return true;
    }

    const bool requireExact = (m_analysisMode==0) ? true : m_requireExactPhotonMultiplicity;
    const bool useEoverP = (m_analysisMode==0) ? true : m_useEoverPTagID;
    const bool strictMode = (m_analysisMode==0);

    // -------- charged selection from ReconstructedParticle --------
    std::vector<ReconstructedParticle> charged;
    std::vector<podio::ObjectID> chargedShowerIds;
    int qsum=0;
    for (int i=0; i<(int)recPars->size(); ++i) {
        auto rp = recPars->at(i);
        if (rp.getCharge()==0) continue;

        if (rp.getShower().isAvailable()) {
            chargedShowerIds.push_back(rp.getShower().getObjectID());
        }

        auto p = rp.getMomentum();
        TVector3 v(p.x,p.y,p.z);
        if (v.Mag() < m_chargedPmin) continue;
        if (!passAccCos(v, m_cosThetaAccCut)) continue;

        charged.push_back(rp);
        qsum += rp.getCharge();
    }
    m_nCharged = (int)charged.size();
    if (strictMode) {
        if (m_nCharged!=2) { ++m_fail_nCharged; return true; }
        if (qsum!=0) { ++m_fail_qsum; return true; }
    } else {
        if (m_nCharged<2) { ++m_fail_nCharged; return true; }
    }

    // -------- photon candidates from ECAL shower (paper acceptance + split merge) --------
    struct PhotonCand { int idx {-1}; TVector3 dir; double E {0.0}; };
    std::vector<PhotonCand> cands;
    cands.reserve(showers->size());
    auto isChargedShower = [&](const RecECALShower& sh) {
        auto sid = sh.getObjectID();
        for (const auto& cid : chargedShowerIds) {
            if (cid == sid) return true;
        }
        return false;
    };
    for (int i=0; i<(int)showers->size(); ++i) {
        const auto& sh = showers->at(i);
        if (isChargedShower(sh)) continue;
        auto pos = sh.getPosition();
        TVector3 dir(pos.x, pos.y, pos.z);
        if (dir.Mag()<=0) continue;
        dir = dir.Unit();
        double absCos = std::fabs(dir.Z());

        bool inBarrel = (absCos < m_photonBarrelCosMax);
        bool inEndcap = (absCos > m_photonEndcapCosMin && absCos < m_photonEndcapCosMax);
        if (!inBarrel && !inEndcap) continue;

        double thr = inBarrel ? m_photonBarrelEmin : m_photonEndcapEmin;
        double E = sh.getEnergy();
        if (E < thr) continue;
        cands.push_back({i, dir, E});
    }

    std::vector<int> gamIdx;
    gamIdx.reserve(cands.size());
    if (!cands.empty()) {
        if (m_photonMergeAngleDeg > 0) {
            std::sort(cands.begin(), cands.end(), [](const PhotonCand& a, const PhotonCand& b){
                return a.E > b.E;
            });
            const double deg2rad = std::acos(-1.0) / 180.0;
            const double cosMerge = std::cos(m_photonMergeAngleDeg * deg2rad);
            std::vector<TVector3> keptDirs;
            keptDirs.reserve(cands.size());
            for (const auto& c : cands) {
                bool merged = false;
                for (const auto& d : keptDirs) {
                    if (c.dir.Dot(d) > cosMerge) { merged = true; break; }
                }
                if (!merged) {
                    gamIdx.push_back(c.idx);
                    keptDirs.push_back(c.dir);
                }
            }
        } else {
            for (const auto& c : cands) gamIdx.push_back(c.idx);
        }
    }
    m_nGamma = (int)gamIdx.size();
    if (m_nGamma<1) { ++m_fail_nGamma; return true; }

    // -------- PID (tau3mu-style + E/p refinement) --------
    auto getPID = [&](const ReconstructedParticle& rp, double& pe, double& pmu, double& ppi){
        m_gPIDSvc->calculate(rp);
        m_gPIDSvc->setmode(m_gPIDSvc->all());
        pe  = m_gPIDSvc->prob(Electron);
        pmu = m_gPIDSvc->prob(Muon);
        ppi = m_gPIDSvc->prob(Pion);
    };

    struct ChargedPID { double pe {0.0}; double pmu {0.0}; double ppi {0.0}; double eop {0.0}; int charge {0}; };
    std::vector<ChargedPID> pidInfo(charged.size());
    for (size_t i=0; i<charged.size(); ++i) {
        getPID(charged[i], pidInfo[i].pe, pidInfo[i].pmu, pidInfo[i].ppi);
        pidInfo[i].eop = calcEoverP(charged[i]);
        pidInfo[i].charge = charged[i].getCharge();
    }

    int idxMu = -1, idxTag = -1, baseTag = 0;
    double bestScore = -1.0;
    bool hasOppCharge = false, hasMuonLike = false, hasMuonEop = false, hasTagLike = false, hasTagEop = false, hasFullPair = false;

    for (size_t i=0; i<charged.size(); ++i) {
        for (size_t j=i+1; j<charged.size(); ++j) {
            if (pidInfo[i].charge * pidInfo[j].charge >= 0) continue;
            hasOppCharge = true;
            for (int assign=0; assign<2; ++assign) {
                int muIdx = (assign==0) ? (int)i : (int)j;
                int tagIdx = (assign==0) ? (int)j : (int)i;
                const auto& mu = pidInfo[muIdx];
                const auto& tag = pidInfo[tagIdx];

                double muMax = std::max(mu.pe, mu.ppi);
                bool muLike = (mu.pmu >= m_muPidMin) && (mu.pmu >= (muMax - m_muPidSlack));
                if (!muLike) continue;
                hasMuonLike = true;

                bool muEopPass = true;
                if (!muEopPass) continue;
                hasMuonEop = true;

                int candTagNoEop = 0;
                if (tag.pe>tag.pmu && tag.pe>tag.ppi) candTagNoEop = 1;
                else if (tag.ppi>tag.pe && tag.ppi>tag.pmu) candTagNoEop = 2;

                int candTagEop = 0;
                if (candTagNoEop==1) {
                    if (tag.eop <= 0) candTagEop = 1;
                    else if (tag.eop > m_eoverpElectronMin) candTagEop = 1;
                } else if (candTagNoEop==2) {
                    if (tag.eop <= 0) candTagEop = 2;
                    else if (tag.eop>0 && tag.eop < m_eoverpMuonMax) candTagEop = 2;
                }

                if (candTagNoEop!=0) hasTagLike = true;
                if (candTagEop!=0) hasTagEop = true;

                int candTag = useEoverP ? candTagEop : candTagNoEop;
                if (candTag==0) continue;
                hasFullPair = true;

                double tagScore = (candTag==1) ? tag.pe : tag.ppi;
                double score = mu.pmu + tagScore;
                if (score > bestScore) {
                    bestScore = score;
                    idxMu = muIdx;
                    idxTag = tagIdx;
                    baseTag = candTag;
                }
            }
        }
    }

    if (idxMu<0 || idxTag<0 || baseTag==0) {
        if (!hasOppCharge) ++m_fail_pid_noOppCharge;
        else if (!hasMuonLike) ++m_fail_pid_noMuonLike;
        else if (useEoverP && !hasMuonEop) ++m_fail_pid_noMuonEop;
        else if (!hasTagLike) ++m_fail_pid_noTagLike;
        else if (useEoverP && !hasTagEop) ++m_fail_pid_noTagEop;
        else ++m_fail_pid_noPair;
        return true;
    }

    const auto& rpMu  = charged[idxMu];
    const auto& rpTag = charged[idxTag];
    TLorentzVector p4_mu  = makeChargedP4(rpMu,  kMassMu);
    const double targetEtau = m_ecms * 0.5;

    // -------- pi0 search (rho-tag + veto) --------
    std::vector<char> isPi0(showers->size(), 0);
    int pi0g1 = -1, pi0g2 = -1;
    struct Pi0Pair { int a; int b; double diff; };
    std::vector<Pi0Pair> pi0Pairs;
    if (m_nGamma>=2) {
        for (size_t a=0; a<gamIdx.size(); ++a) {
            for (size_t b=a+1; b<gamIdx.size(); ++b) {
                TLorentzVector p1 = makePhotonP4(showers->at(gamIdx[a]));
                TLorentzVector p2 = makePhotonP4(showers->at(gamIdx[b]));
                double mgg = (p1+p2).M();
                if (mgg < m_pi0MassWinLow || mgg > m_pi0MassWinHigh) continue;
                double diff = std::fabs(mgg-kMassPi0);
                pi0Pairs.push_back({gamIdx[a], gamIdx[b], diff});
            }
        }
    }

    int sigG = -1;
    int tagMode = baseTag;
    bool foundPi0 = false;

    if (!pi0Pairs.empty()) {
        if (baseTag==2) {
            double bestSigScore = std::numeric_limits<double>::max();
            for (const auto& pair : pi0Pairs) {
                std::vector<int> remaining;
                for (int id : gamIdx) {
                    if (id==pair.a || id==pair.b) continue;
                    remaining.push_back(id);
                }
                if (requireExact && remaining.size()!=1) continue;
                int candSig = requireExact ? (remaining.empty() ? -1 : remaining[0])
                    : pickBestSignalPhoton(showers, remaining, p4_mu, targetEtau);
                if (candSig<0) continue;
                double score = signalScore(showers, candSig, p4_mu, targetEtau);
                if (score < bestSigScore) {
                    bestSigScore = score;
                    pi0g1 = pair.a;
                    pi0g2 = pair.b;
                }
            }
        } else {
            auto best = std::min_element(pi0Pairs.begin(), pi0Pairs.end(),
                [](const Pi0Pair& x, const Pi0Pair& y){ return x.diff < y.diff; });
            if (best != pi0Pairs.end()) { pi0g1 = best->a; pi0g2 = best->b; }
        }
        foundPi0 = (pi0g1>=0 && pi0g2>=0);
        if (foundPi0) { isPi0[pi0g1] = 1; isPi0[pi0g2] = 1; }
    }

    if (baseTag==2 && foundPi0) {
        std::vector<int> remaining;
        for (int id : gamIdx) {
            if (id==pi0g1 || id==pi0g2) continue;
            remaining.push_back(id);
        }
        if (!remaining.empty()) {
            if (requireExact) {
                if (remaining.size()==1) sigG = remaining[0];
            } else {
                sigG = pickBestSignalPhoton(showers, remaining, p4_mu, targetEtau);
                if (sigG>=0) m_nExtraGamma = (int)remaining.size() - 1;
            }
        }
        if (sigG>=0) tagMode = 3;
        else if (strictMode) { ++m_fail_rho_noSigG; return true; }
    }

    if (tagMode!=3) {
        std::vector<int> candidates;
        for (int id : gamIdx) {
            if (!isPi0[id]) candidates.push_back(id);
        }
        const std::vector<int>& useCandidates = candidates.empty() ? gamIdx : candidates;
        if (requireExact) {
            if (useCandidates.size()!=1) { ++m_fail_sigG_exact; return true; }
            sigG = useCandidates[0];
        } else {
            sigG = pickBestSignalPhoton(showers, useCandidates, p4_mu, targetEtau);
            if (sigG>=0) m_nExtraGamma = (int)useCandidates.size() - 1;
        }
    }

    if (sigG<0) { ++m_fail_sigG_none; return true; }
    m_tagMode = tagMode;

    // -------- build 4-momenta --------
    double tagMass = (m_tagMode==1) ? kMassE : kMassPi;

    TLorentzVector p4_tag = makeChargedP4(rpTag, tagMass);
    TLorentzVector p4_sigG= makePhotonP4(showers->at(sigG));

    m_pMu = p4_mu.Vect().Mag();
    m_cosMuGamma = cosAngle(p4_mu.Vect(), p4_sigG.Vect());

    m_EsigGamma = p4_sigG.E();
    bool passCosMuGamma = (m_cosMuGamma < m_cosThetaMuGammaCut);
    bool passPMu = (m_pMu >= m_pMuMin && m_pMu <= m_pMuMax);
    bool passEsig = (m_EsigGamma >= m_EsigGammaMin && m_EsigGamma <= m_EsigGammaMax);
    if (!passCosMuGamma) ++m_fail_kin_cosMuGamma;
    if (!passPMu) ++m_fail_kin_pMu;
    if (!passEsig) ++m_fail_kin_EsigGamma;
    bool passKinematics = passCosMuGamma && passPMu && passEsig;

    m_invMass_Signal = (p4_mu+p4_sigG).M();
    m_E_mugamma = (p4_mu+p4_sigG).E();
    m_pTag = p4_tag.Vect().Mag();
    m_cosSigGamma_TagCharged = cosAngle(p4_sigG.Vect(), p4_tag.Vect());

    // missing 4-momentum (simple: Pbeams=(0,0,0,Ec.m.))
    TLorentzVector P_beams(0,0,0,m_ecms);

    TLorentzVector P_sys = p4_mu + p4_tag + p4_sigG;
    if (m_tagMode==3) {
        if (pi0g1<0 || pi0g2<0) return true;
        P_sys += makePhotonP4(showers->at(pi0g1));
        P_sys += makePhotonP4(showers->at(pi0g2));
    }
    TLorentzVector P_miss = P_beams - P_sys;

    m_E_sys = P_sys.E();
    m_P_sys = P_sys.Vect().Mag();
    m_E_miss = P_miss.E();
    m_P_miss = P_miss.P();
    m_check_miss = (m_E_miss >= m_P_miss) ? 1 : 0;

    m_Emiss = m_E_miss;
    m_M2miss = P_miss.M2();
    m_cosTheta_miss = (P_miss.P()>0) ? (P_miss.Pz()/P_miss.P()) : 1.0;

    // helicity angle (approx, signal tau = mu+gamma)
    {
        TLorentzVector p4_tau = p4_mu + p4_sigG;
        TVector3 tau_dir = p4_tau.Vect();
        if (tau_dir.Mag()>0) tau_dir = tau_dir.Unit();

        TLorentzVector mu_in_tau = p4_mu;
        mu_in_tau.Boost(-p4_tau.BoostVector());
        TVector3 mu_dir = mu_in_tau.Vect();
        if (mu_dir.Mag()>0) mu_dir = mu_dir.Unit();

        m_cosHel = mu_dir.Dot(tau_dir);
    }

    if (passKinematics) ++m_evtPassBasic;
    m_passBasic = passKinematics ? 1 : 0;

    // -------- Table-1 further cuts --------
    bool passFinal = passKinematics;
    if (passKinematics) {
        if (m_tagMode==1) {
            bool passCos = (m_cosSigGamma_TagCharged < m_cut_cosSigGamma_TagCharged);
            bool passPTag = (m_pTag > m_cut_pTagMin);
            bool passEmiss = (m_Emiss < m_cut_EmissMax_e);
            if (!passCos) ++m_fail_e_cosSigGamma;
            if (!passPTag) ++m_fail_e_pTag;
            if (!passEmiss) ++m_fail_e_Emiss;
            passFinal = passFinal && passCos && passPTag && passEmiss;
        } else if (m_tagMode==2) {
            bool passEmiss = (m_Emiss > m_cut_EmissMin_pi);
            bool passCosMiss = (std::fabs(m_cosTheta_miss) < m_cut_absCosMissMax_pi);
            bool passEsigPi = (m_EsigGamma > m_cut_EsigGammaMin_pi);
            bool passM2 = (m_M2miss < m_cut_M2missMax_pi);
            if (!passEmiss) ++m_fail_pi_Emiss;
            if (!passCosMiss) ++m_fail_pi_cosMiss;
            if (!passEsigPi) ++m_fail_pi_EsigGamma;
            if (!passM2) ++m_fail_pi_M2miss;
            passFinal = passFinal && passEmiss && passCosMiss && passEsigPi && passM2;
        } else if (m_tagMode==3) {
            bool passM2 = (m_M2miss < m_cut_M2missMax_rho);
            bool passCosHel = (m_cosHel < m_cut_cosHelMax_rho);
            bool passCosMiss = (std::fabs(m_cosTheta_miss) < m_cut_absCosMissMax_rho);
            if (!passM2) ++m_fail_rho_M2miss;
            if (!passCosHel) ++m_fail_rho_cosHel;
            if (!passCosMiss) ++m_fail_rho_cosMiss;
            passFinal = passFinal && passM2 && passCosHel && passCosMiss;
        } else {
            ++m_fail_tagModeUnknown;
            passFinal = false;
        }
    }

    m_passFinal = passFinal ? 1 : 0;
    if (passFinal) ++m_evtPassFinal;
    m_tree->Fill();
    return true;
}

bool TauMuGammaAlg::finalize() {
    LogInfo << "TauMuGammaAlg summary:\n"
            << "  All events     : " << m_evtAll << "\n"
            << "  Pass basic     : " << m_evtPassBasic << "\n"
            << "  Pass final fill: " << m_evtPassFinal << "\n"
            << "  Fail no rec/shower: " << m_fail_noRecOrShower << "\n"
            << "  Fail nCharged     : " << m_fail_nCharged << "\n"
            << "  Fail qsum         : " << m_fail_qsum << "\n"
            << "  Fail nGamma       : " << m_fail_nGamma << "\n"
            << "  Fail PID no oppQ  : " << m_fail_pid_noOppCharge << "\n"
            << "  Fail PID mu-like  : " << m_fail_pid_noMuonLike << "\n"
            << "  Fail PID mu E/p   : " << m_fail_pid_noMuonEop << "\n"
            << "  Fail PID tag-like : " << m_fail_pid_noTagLike << "\n"
            << "  Fail PID tag E/p  : " << m_fail_pid_noTagEop << "\n"
            << "  Fail PID pair     : " << m_fail_pid_noPair << "\n"
            << "  Fail rho sigG     : " << m_fail_rho_noSigG << "\n"
            << "  Fail sigG exact   : " << m_fail_sigG_exact << "\n"
            << "  Fail sigG none    : " << m_fail_sigG_none << "\n"
            << "  Fail kin cosMuG   : " << m_fail_kin_cosMuGamma << "\n"
            << "  Fail kin pMu      : " << m_fail_kin_pMu << "\n"
            << "  Fail kin EsigG    : " << m_fail_kin_EsigGamma << "\n"
            << "  Fail e cosSigG    : " << m_fail_e_cosSigGamma << "\n"
            << "  Fail e pTag       : " << m_fail_e_pTag << "\n"
            << "  Fail e Emiss      : " << m_fail_e_Emiss << "\n"
            << "  Fail pi Emiss     : " << m_fail_pi_Emiss << "\n"
            << "  Fail pi cosMiss   : " << m_fail_pi_cosMiss << "\n"
            << "  Fail pi EsigG     : " << m_fail_pi_EsigGamma << "\n"
            << "  Fail pi M2miss    : " << m_fail_pi_M2miss << "\n"
            << "  Fail rho M2miss   : " << m_fail_rho_M2miss << "\n"
            << "  Fail rho cosHel   : " << m_fail_rho_cosHel << "\n"
            << "  Fail rho cosMiss  : " << m_fail_rho_cosMiss << "\n"
            << "  Fail tag unknown  : " << m_fail_tagModeUnknown << std::endl;
    return true;
}
```

## TauMuGammaAna/python/TauMuGammaAna/__init__.py

```python
import Sniper
Sniper.loadDll("libTauMuGammaAna.so")
```

## plotting/__init__.py

```python
"""plotting 包的初始化文件。

对 Python 初学者的说明：
- 这个文件存在的主要意义是告诉 Python：这个文件夹是一个可导入的包（package）。
- 当你写 `import plotting` 或 `from plotting import xxx` 时，Python 会先执行这里的代码。
- 如果这里是空的也完全没问题，只是起到"标记为包"的作用。
"""

# 当前项目不需要在导入时自动执行任何代码，所以保持为空即可。
```

## plotting/io.py

```python
"""IO 辅助工具：负责从 ROOT 文件中读取数据，并把多个文件的数组拼接起来。

对 Python 初学者的说明：
- 这里的函数只做"数据读取/合并"，不做任何绘图。
- 通过拆分成小函数，代码可复用、易测试、也更清晰。
"""

import numpy as np
import uproot


def load_arrays(file_path, tree_path, branches):
    """从 ROOT 文件读取指定分支，并返回 NumPy 数组字典。

    参数说明：
    - file_path：ROOT 文件路径。
    - tree_path：树的路径，例如 "Fkey/TauMuGamma"。
    - branches：分支名列表（要读出的变量）。

    返回：
    - dict，键是分支名，值是 NumPy 数组。
    """
    with uproot.open(file_path) as root_file:
        tree = root_file[tree_path]
        arrays = tree.arrays(branches, library="np")
    return arrays


def concat_arrays(arrays_list):
    """把多个文件读取出的数组字典按键拼接成一个大字典。"""
    if not arrays_list:
        return {}
    keys = arrays_list[0].keys()
    out = {}
    for key in keys:
        out[key] = np.concatenate([arr[key] for arr in arrays_list])
    return out
```

## plotting/selection.py

```python
"""事件选择（selection）相关的条件与掩码函数。

对 Python 初学者的说明：
- 本文件只负责"选择条件"，不负责读取数据或绘图。
- 通过布尔掩码（boolean mask）来筛选数组，这是一种常见且高效的做法。
"""

import numpy as np


DEFAULT_CUTS = {
    "cosMuGammaCut": -0.35,
    "pMuMin": 0.4,
    "pMuMax": 1.7,
    "EsigGammaMin": 0.4,
    "EsigGammaMax": 1.7,
    # e-tag
    "cosSigGamma_TagCharged": -0.2,
    "pTagMin": 0.5,
    "EmissMax_e": 1.7,
    # pi-tag
    "EmissMin_pi": 0.7,
    "absCosMissMax_pi": 0.6,
    "EsigGammaMin_pi": 0.8,
    "M2missMax_pi": 0.050,
    # rho-tag
    "M2missMax_rho": 0.075,
    "cosHelMax_rho": 0.8,
    "absCosMissMax_rho": 0.9,
}


def basic_mask(arr):
    """基础选择：只保留 `passBasic == 1` 的事件。"""
    return arr["passBasic"] == 1


def final_mask(arr):
    """最终选择：只保留 `passFinal == 1` 的事件。"""
    return arr["passFinal"] == 1


def kinematic_mask(arr, cuts=DEFAULT_CUTS, drop=None):
    """运动学选择：对 pMu、EsigGamma、cosMuGamma 施加范围限制。"""
    mask = np.ones(len(arr["pMu"]), dtype=bool)
    if drop != "pMu":
        mask &= (arr["pMu"] >= cuts["pMuMin"]) & (arr["pMu"] <= cuts["pMuMax"])
    if drop != "EsigGamma":
        mask &= (arr["EsigGamma"] >= cuts["EsigGammaMin"]) & (arr["EsigGamma"] <= cuts["EsigGammaMax"])
    if drop != "cosMuGamma":
        mask &= (arr["cosMuGamma"] < cuts["cosMuGammaCut"])
    return mask


def tag_mask(arr, tag_mode, cuts=DEFAULT_CUTS, drop=None):
    """按 tag 模式选择事件，并叠加相应的 cut。"""
    mask = (arr["tagMode"] == tag_mode)
    mask &= kinematic_mask(arr, cuts=cuts)

    if tag_mode == 1:
        if drop != "cosSigGamma_TagCharged":
            mask &= (arr["cosSigGamma_TagCharged"] < cuts["cosSigGamma_TagCharged"])
        if drop != "pTag":
            mask &= (arr["pTag"] > cuts["pTagMin"])
        if drop != "Emiss":
            mask &= (arr["Emiss"] < cuts["EmissMax_e"])
    elif tag_mode == 2:
        if drop != "Emiss":
            mask &= (arr["Emiss"] > cuts["EmissMin_pi"])
        if drop != "absCosTheta_miss":
            mask &= (np.abs(arr["cosTheta_miss"]) < cuts["absCosMissMax_pi"])
        if drop != "EsigGamma":
            mask &= (arr["EsigGamma"] > cuts["EsigGammaMin_pi"])
        if drop != "M2miss":
            mask &= (arr["M2miss"] < cuts["M2missMax_pi"])
    elif tag_mode == 3:
        if drop != "M2miss":
            mask &= (arr["M2miss"] < cuts["M2missMax_rho"])
        if drop != "cosHel":
            mask &= (arr["cosHel"] < cuts["cosHelMax_rho"])
        if drop != "absCosTheta_miss":
            mask &= (np.abs(arr["cosTheta_miss"]) < cuts["absCosMissMax_rho"])
    return mask
```

## plotting/plot_paper_figs.py

```python
#!/usr/bin/env python3
"""用于绘制论文中各类分布图的主脚本。

纵轴单位约定（与论文一致）：
- "Events / (0.05 GeV)" 等：纵轴表示事件密度 = 每 bin 计数 / bin 宽度。
- "×10³"：纵轴数值先除以 1000 再标刻度，可读性更好。
"""

import argparse
import json
import math
import os

import numpy as np
import ROOT

from plotting.io import load_arrays, concat_arrays
from plotting import selection


PLOTS_1D = [
    # Fig.2: kinematics
    {"name": "pMu", "var": "pMu", "bins": 60, "bin_width": 0.05, "range": (0.0, 2.0), "drop": "pMu",
     "xlabel": "p_{#mu} (GeV/c)", "y_unit": "0.05 GeV/c"},
    {"name": "EsigGamma", "var": "EsigGamma", "bins": 60, "bin_width": 0.05, "range": (0.0, 2.0), "drop": "EsigGamma",
     "xlabel": "E_{#gamma} (GeV)"},
    {"name": "cosMuGamma", "var": "cosMuGamma", "bins": 60, "bin_width": 0.05, "range": (-1.0, 1.0), "drop": "cosMuGamma",
     "xlabel": "cos#theta_{#gamma#mu}", "y_unit": "0.05"},
    # Fig.3: e-tag（cos#theta_{#gamma,tag} 纵轴为对数）
    {"name": "e_cosSigGamma_TagCharged", "var": "cosSigGamma_TagCharged", "bins": 60, "bin_width": 0.05, "range": (-1.0, 1.0),
     "tag_mode": 1, "drop": "cosSigGamma_TagCharged", "xlabel": "cos#theta_{#gamma,tag}", "y_unit": "0.05", "log_y": True},
    {"name": "e_pTag", "var": "pTag", "bins": 60, "bin_width": 0.05, "range": (0.0, 2.0),
     "tag_mode": 1, "drop": "pTag", "xlabel": "p_{tag} (GeV/c)", "y_unit": "0.05 GeV/c"},
    {"name": "e_Emiss", "var": "Emiss", "bins": 60, "bin_width": 0.05, "range": (0.0, 2.5),
     "tag_mode": 1, "drop": "Emiss", "xlabel": "E_{miss} (GeV)"},
    # Fig.4: pi-tag
    {"name": "pi_Emiss", "var": "Emiss", "bins": 60, "bin_width": 0.05, "range": (0.0, 2.5),
     "tag_mode": 2, "drop": "Emiss", "xlabel": "E_{miss} (GeV)"},
    {"name": "pi_absCosMiss", "var": "cosTheta_miss", "bins": 60, "bin_width": 0.05, "range": (0.0, 1.0),
     "tag_mode": 2, "drop": "absCosTheta_miss", "transform": "abs", "xlabel": "|cos#theta_{miss}|", "y_unit": "0.05", "log_y": True},
    {"name": "pi_EsigGamma", "var": "EsigGamma", "bins": 60, "bin_width": 0.05, "range": (0.0, 2.0),
     "tag_mode": 2, "drop": "EsigGamma", "xlabel": "E_{#gamma} (GeV)"},
    {"name": "pi_M2miss", "var": "M2miss", "bins": 60, "bin_width": 0.05, "range": (-0.2, 0.2),
     "tag_mode": 2, "drop": "M2miss", "xlabel": "M^{2}_{miss} (GeV^{2})", "y_unit": "0.05 GeV^{2}/c^{4}"},
    # Fig.5: rho-tag
    {"name": "rho_M2miss", "var": "M2miss", "bins": 60, "bin_width": 0.05, "range": (-0.2, 0.2),
     "tag_mode": 3, "drop": "M2miss", "xlabel": "M^{2}_{miss} (GeV^{2})", "y_unit": "0.05 GeV^{2}/c^{4}"},
    {"name": "rho_cosHel", "var": "cosHel", "bins": 60, "bin_width": 0.05, "range": (-1.0, 1.0),
     "tag_mode": 3, "drop": "cosHel", "xlabel": "cos#theta_{hel}", "y_unit": "0.05"},
    {"name": "rho_absCosMiss", "var": "cosTheta_miss", "bins": 60, "bin_width": 0.05, "range": (0.0, 1.0),
     "tag_mode": 3, "drop": "absCosTheta_miss", "transform": "abs", "xlabel": "|cos#theta_{miss}|", "y_unit": "0.05"},
]

Y_BIN_WIDTH_GEV = 0.05
Y_SCALE = 1.0e-3

PLOT_2D = {
    "name": "EvsM_mugamma",
    "xvar": "E_mugamma",
    "yvar": "M_mugamma",
    "xbins": 60,
    "ybins": 60,
    "xrange": (1.8, 2.2),
    "yrange": (1.5, 1.9),
    "xlabel": "E_{#gamma#mu} (GeV)",
    "ylabel": "M_{#gamma#mu} (GeV/c^{2})",
    "draw_style": "BOX",
}


def parse_color(color):
    if isinstance(color, (int, float)):
        return int(color)
    if not isinstance(color, str):
        return ROOT.kBlack
    if "+" in color or "-" in color:
        for op in ("+", "-"):
            if op in color:
                base, off = color.split(op, 1)
                base_val = getattr(ROOT, base, ROOT.kBlack)
                off_val = int(off)
                return base_val + off_val if op == "+" else base_val - off_val
    return getattr(ROOT, color, ROOT.kBlack)


def sample_weight(sample, cfg):
    weight = sample.get("scale")
    if weight is None:
        xsec_pb = sample.get("xsec_pb", 0.0)
        n_gen = sample.get("n_gen", 0)
        lumi_ab = cfg.get("lumi_ab", 1.0)
        if xsec_pb > 0 and n_gen > 0:
            lumi_pb = lumi_ab * 1.0e6
            weight = xsec_pb * lumi_pb / n_gen
        else:
            weight = 1.0
    if cfg.get("tag_forced", False) and sample.get("apply_tag_br", False):
        weight *= cfg.get("br_tag_sum", 0.541)
    return weight


def make_hist1d(name, title, values, weights, bins, xmin, xmax):
    hist = ROOT.TH1D(name, title, bins, xmin, xmax)
    hist.Sumw2()
    if values.size == 0:
        return hist
    counts, edges = np.histogram(values, bins=bins, range=(xmin, xmax), weights=weights)
    counts2, _ = np.histogram(values, bins=bins, range=(xmin, xmax), weights=weights * weights)
    for i in range(bins):
        hist.SetBinContent(i + 1, counts[i])
        hist.SetBinError(i + 1, math.sqrt(counts2[i]))
    return hist


def normalize_hist_y(hist, target_bin_width, y_scale):
    """按目标 bin 宽度和额外比例缩放直方图。"""
    bin_width = hist.GetXaxis().GetBinWidth(1)
    if bin_width > 0:
        hist.Scale(target_bin_width / bin_width)
    if y_scale != 1.0:
        hist.Scale(y_scale)


def make_hist2d(name, title, xvals, yvals, weights, xbins, ybins, xmin, xmax, ymin, ymax):
    hist = ROOT.TH2D(name, title, xbins, xmin, xmax, ybins, ymin, ymax)
    if xvals.size == 0:
        return hist
    counts, xedges, yedges = np.histogram2d(
        xvals, yvals, bins=(xbins, ybins), range=((xmin, xmax), (ymin, ymax)), weights=weights
    )
    for ix in range(xbins):
        for iy in range(ybins):
            hist.SetBinContent(ix + 1, iy + 1, counts[ix, iy])
    return hist


def ellipse_graph(params, npoints=200):
    if not params.get("enabled"):
        return None
    center = params.get("center")
    if not center or center[0] is None or center[1] is None:
        return None
    cx, cy = center
    angle = math.radians(params.get("angle_deg", 0.0))

    rx = params.get("rx")
    ry = params.get("ry")
    rx_neg = params.get("rx_neg", rx)
    rx_pos = params.get("rx_pos", rx)
    ry_neg = params.get("ry_neg", ry)
    ry_pos = params.get("ry_pos", ry)

    if None in (rx, ry, rx_neg, rx_pos, ry_neg, ry_pos):
        return None

    xs = []
    ys = []
    for t in np.linspace(0.0, 2.0 * math.pi, npoints, endpoint=True):
        ux = math.cos(t)
        uy = math.sin(t)
        rx_use = rx_pos if ux >= 0 else rx_neg
        ry_use = ry_pos if uy >= 0 else ry_neg
        x_local = ux * rx_use
        y_local = uy * ry_use
        x = cx + x_local * math.cos(angle) - y_local * math.sin(angle)
        y = cy + x_local * math.sin(angle) + y_local * math.cos(angle)
        xs.append(x)
        ys.append(y)
    graph = ROOT.TGraph(npoints)
    for i, (x, y) in enumerate(zip(xs, ys)):
        graph.SetPoint(i, x, y)
    graph.SetLineWidth(2)
    graph.SetLineColor(ROOT.kBlack)
    return graph


def build_mask(arr, plot_def):
    """根据图的配置决定使用哪种选择条件。"""
    tag_mode = plot_def.get("tag_mode")
    drop = plot_def.get("drop")
    if tag_mode is None:
        return selection.kinematic_mask(arr, drop=drop)
    return selection.tag_mask(arr, tag_mode, drop=drop)


def main():
    """程序入口：读取配置 → 读取数据 → 绘图输出。"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="plotting/samples.json")
    parser.add_argument("--ellipse", default="plotting/ellipse_params.json")
    parser.add_argument("--outdir", default="plotting/output")
    parser.add_argument("--root-out", default="plotting/plots.root")
    parser.add_argument("--no-pdf", action="store_true")
    parser.add_argument("--tree-path", default=None)
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)

    with open(args.config, "r") as f:
        cfg = json.load(f)
    tree_path = args.tree_path or cfg.get("tree_path", "Fkey/TauMuGamma")

    os.makedirs(args.outdir, exist_ok=True)

    branches = [
        "pMu", "EsigGamma", "cosMuGamma",
        "tagMode", "Emiss", "M2miss", "cosTheta_miss",
        "cosSigGamma_TagCharged", "pTag", "cosHel",
        "E_mugamma", "M_mugamma", "passBasic", "passFinal"
    ]

    sample_arrays = {}
    for sample in cfg["samples"]:
        files = sample["file"]
        if isinstance(files, str):
            files = [files]
        arrays_list = [load_arrays(path, tree_path, branches) for path in files]
        sample_arrays[sample["name"]] = concat_arrays(arrays_list)

    out_file = ROOT.TFile(args.root_out, "RECREATE")

    for plot_def in PLOTS_1D:
        plot_name = plot_def["name"]
        plot_dir = out_file.mkdir(plot_name)
        plot_dir.cd()

        signal_hist = None
        stack = ROOT.THStack("stack_" + plot_name, "")
        legend = ROOT.TLegend(0.65, 0.7, 0.88, 0.88)
        legend.SetBorderSize(0)

        for sample in cfg["samples"]:
            name = sample["name"]
            arr = sample_arrays[name]
            mask = build_mask(arr, plot_def)
            values = arr[plot_def["var"]][mask]
            if plot_def.get("transform") == "abs":
                values = np.abs(values)
            weight = sample_weight(sample, cfg)
            weights = np.full(values.shape, weight, dtype=float)

            bins = plot_def.get("bins")
            if plot_def.get("bin_width"):
                bw = plot_def["bin_width"]
                bins = max(1, int(round((plot_def["range"][1] - plot_def["range"][0]) / bw)))

            hist = make_hist1d(
                name,
                plot_name,
                values,
                weights,
                bins,
                plot_def["range"][0],
                plot_def["range"][1],
            )
            y_scale = 1.0 if plot_def.get("log_y") else Y_SCALE
            normalize_hist_y(hist, Y_BIN_WIDTH_GEV, y_scale)
            color = parse_color(sample.get("color", "kBlack"))
            hist.SetLineColor(color)

            if sample.get("is_signal", False):
                hist.SetLineWidth(2)
                hist.SetFillStyle(0)
                signal_hist = hist
                legend.AddEntry(hist, sample["label"], "l")
            else:
                hist.SetFillColor(color)
                hist.SetLineWidth(1)
                stack.Add(hist)
                legend.AddEntry(hist, sample["label"], "f")

            hist.Write(name)

        if not args.no_pdf:
            canvas = ROOT.TCanvas("c_" + plot_name, plot_name, 800, 600)
            if plot_def.get("log_y"):
                canvas.SetLogy(1)
            hists = stack.GetHists()
            has_bkg = bool(hists) and hists.GetSize() > 0
            y_unit = plot_def.get("y_unit", "0.05 GeV")
            y_label = "Events / (%s)" % y_unit if plot_def.get("log_y") else "Events / (%s) #times 10^{3}" % y_unit
            if has_bkg:
                stack.Draw("hist")
                stack.GetXaxis().SetTitle(plot_def["xlabel"])
                stack.GetYaxis().SetTitle(y_label)
                if signal_hist:
                    signal_hist.Draw("hist same")
            elif signal_hist:
                signal_hist.Draw("hist")
                signal_hist.GetXaxis().SetTitle(plot_def["xlabel"])
                signal_hist.GetYaxis().SetTitle(y_label)
            else:
                continue
            legend.Draw()
            canvas.SaveAs(os.path.join(args.outdir, plot_name + ".pdf"))
            canvas.SaveAs(os.path.join(args.outdir, plot_name + ".png"))

        out_file.cd()

    # 2D 图：E_mugamma vs M_mugamma（应用 final selection）
    with open(args.ellipse, "r") as f:
        ellipse_params = json.load(f)
    ellipse = ellipse_graph(ellipse_params)

    plot2d = PLOT_2D
    plot_dir = out_file.mkdir(plot2d["name"])
    plot_dir.cd()

    for sample in cfg["samples"]:
        name = sample["name"]
        arr = sample_arrays[name]
        mask = selection.final_mask(arr)
        xvals = arr[plot2d["xvar"]][mask]
        yvals = arr[plot2d["yvar"]][mask]
        weight = sample_weight(sample, cfg)
        weights = np.full(xvals.shape, weight, dtype=float)

        hist2d = make_hist2d(
            name,
            plot2d["name"],
            xvals,
            yvals,
            weights,
            plot2d["xbins"],
            plot2d["ybins"],
            plot2d["xrange"][0],
            plot2d["xrange"][1],
            plot2d["yrange"][0],
            plot2d["yrange"][1],
        )
        hist2d.SetXTitle(plot2d["xlabel"])
        hist2d.SetYTitle(plot2d["ylabel"])
        hist2d.SetStats(0)
        hist2d.Write(name)

        if not args.no_pdf:
            canvas = ROOT.TCanvas("c2d_" + plot2d["name"] + "_" + name, plot2d["name"], 800, 700)
            draw_opt = plot2d.get("draw_style", "BOX")
            hist2d.SetFillStyle(0)
            hist2d.SetLineColor(ROOT.kBlack)
            hist2d.Draw(draw_opt)
            if ellipse:
                ellipse.Draw("L same")
            canvas.SaveAs(os.path.join(args.outdir, plot2d["name"] + "_" + name + ".pdf"))
            canvas.SaveAs(os.path.join(args.outdir, plot2d["name"] + "_" + name + ".png"))

    if ellipse:
        ellipse.Write("signal_ellipse")

    out_file.Close()


if __name__ == "__main__":
    main()
```

## plotting/samples.json

```json
{
  "tree_path": "Fkey/TauMuGamma",
  "lumi_ab": 1.0,
  "tag_forced": true,
  "br_tag_sum": 0.541,
  "samples": [
    {
      "name": "signal",
      "label": "Signal",
      "file": "path/to/signal.root",
      "is_signal": true,
      "xsec_pb": 0.0,
      "n_gen": 0,
      "scale": null,
      "apply_tag_br": true,
      "color": "kRed"
    },
    {
      "name": "dimu",
      "label": "Di-#mu",
      "file": "path/to/dimu.root",
      "is_signal": false,
      "xsec_pb": 0.0,
      "n_gen": 0,
      "scale": null,
      "apply_tag_br": false,
      "color": "kBlue"
    },
    {
      "name": "ditau",
      "label": "#tau#tau",
      "file": "path/to/ditau.root",
      "is_signal": false,
      "xsec_pb": 0.0,
      "n_gen": 0,
      "scale": null,
      "apply_tag_br": true,
      "color": "kGreen+2"
    },
    {
      "name": "had",
      "label": "q#bar{q}",
      "file": "path/to/had.root",
      "is_signal": false,
      "xsec_pb": 0.0,
      "n_gen": 0,
      "scale": null,
      "apply_tag_br": false,
      "color": "kOrange+7"
    }
  ]
}
```

## plotting/ellipse_params.json

```json
{
  "enabled": false,
  "center": [null, null],
  "angle_deg": 0.0,
  "rx": null,
  "ry": null,
  "rx_neg": null,
  "rx_pos": null,
  "ry_neg": null,
  "ry_pos": null
}
```

## .vscode/c_cpp_properties.json

```json
{
  "configurations": [
    {
      "name": "macos-clang-arm64",
      "includePath": [
        "${workspaceFolder}/**"
      ],
      "compilerPath": "/usr/bin/clang",
      "cStandard": "${default}",
      "cppStandard": "${default}",
      "intelliSenseMode": "macos-clang-arm64",
      "compilerArgs": [
        ""
      ]
    }
  ],
  "version": 4
}
```

## .vscode/launch.json

```json
{
  "version": "0.2.0",
  "configurations": [
    {
      "name": "C/C++ Runner: Debug Session",
      "type": "lldb",
      "request": "launch",
      "args": [],
      "cwd": "/Users/david/Desktop/毕设/tau->gamma mu/TauMuGammaAna",
      "program": "/Users/david/Desktop/毕设/tau->gamma mu/TauMuGammaAna/build/Debug/outDebug"
    }
  ]
}
```

## .vscode/settings.json

```json
{
  "C_Cpp_Runner.cCompilerPath": "clang",
  "C_Cpp_Runner.cppCompilerPath": "clang++",
  "C_Cpp_Runner.debuggerPath": "lldb",
  "C_Cpp_Runner.cStandard": "",
  "C_Cpp_Runner.cppStandard": "",
  "C_Cpp_Runner.msvcBatchPath": "",
  "C_Cpp_Runner.useMsvc": false,
  "C_Cpp_Runner.warnings": [
    "-Wall",
    "-Wextra",
    "-Wpedantic",
    "-Wshadow",
    "-Wformat=2",
    "-Wcast-align",
    "-Wconversion",
    "-Wsign-conversion",
    "-Wnull-dereference"
  ],
  "C_Cpp_Runner.msvcWarnings": [
    "/W4",
    "/permissive-",
    "/w14242",
    "/w14287",
    "/w14296",
    "/w14311",
    "/w14826",
    "/w44062",
    "/w44242",
    "/w14905",
    "/w14906",
    "/w14263",
    "/w44265",
    "/w14928"
  ],
  "C_Cpp_Runner.enableWarnings": true,
  "C_Cpp_Runner.warningsAsError": false,
  "C_Cpp_Runner.compilerArgs": [],
  "C_Cpp_Runner.linkerArgs": [],
  "C_Cpp_Runner.includePaths": [],
  "C_Cpp_Runner.includeSearch": [
    "*",
    "**/*"
  ],
  "C_Cpp_Runner.excludeSearch": [
    "**/build",
    "**/build/**",
    "**/.*",
    "**/.*/**",
    "**/.vscode",
    "**/.vscode/**"
  ],
  "C_Cpp_Runner.useAddressSanitizer": false,
  "C_Cpp_Runner.useUndefinedSanitizer": false,
  "C_Cpp_Runner.useLeakSanitizer": false,
  "C_Cpp_Runner.showCompilationTime": false,
  "C_Cpp_Runner.useLinkTimeOptimization": false,
  "C_Cpp_Runner.msvcSecureNoWarnings": false
}
```
