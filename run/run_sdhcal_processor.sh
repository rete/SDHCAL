#!/bin/bash

singleParticleEnergy=$1
runNumberSingleParticle=""

if [ ${singleParticleEnergy} == "10" ]
then
  runNumberSingleParticle="715693"
elif [ ${singleParticleEnergy} == "20" ]
then
  runNumberSingleParticle="715675"
elif [ ${singleParticleEnergy} == "30" ]
then
  runNumberSingleParticle="715747"
elif [ ${singleParticleEnergy} == "40" ]
then
  runNumberSingleParticle="715748"
elif [ ${singleParticleEnergy} == "50" ]
then
  runNumberSingleParticle="715751"
elif [ ${singleParticleEnergy} == "60" ]
then
  runNumberSingleParticle="715753"
elif [ ${singleParticleEnergy} == "70" ]
then
  runNumberSingleParticle="715754"
elif [ ${singleParticleEnergy} == "80" ]
then
  runNumberSingleParticle="715756"
else
  echo "bad input energy for single particle !"
  exit 1
fi


collectionName="SDHCAL_HIT"
fileName="/home/rete/data/sdhcaldata/testbeam/SeptAug2012/CUT/TDHCAL_${runNumberSingleParticle}_cut_spillcut.slcio"
#fileName="/home/rete/data/sdhcaldata/testbeam/SeptAug2012/TRIVENT/TDHCAL_${runNumberSingleParticle}.slcio"
gearFile="/home/rete/soft/SDHCAL/xml/SDHCALGearFile.xml"
rootFileName="ROOTOutputFile_single_pi-_${singleParticleEnergy}GeV.root"

source /home/rete/soft/ilcsoft/v01-17-08/init_ilcsoft.sh
export MARLIN_DLL=/home/rete/soft/SDHCAL/lib/libSDHCALProcessor.so:/home/rete/soft/SDHCAL/lib/libCutProcessor.so


cat > MarlinXML.xml << EOF

<marlin>

	<execute>
	  <!-- <processor name="MyCutProcessor"/> -->
	  <processor name="MySDHCALProcessor"/>
	</execute>

	<global>
	  <parameter name="LCIOInputFiles">
		${fileName}
	  </parameter>
	  <parameter name="GearXMLFile" value="${gearFile}"/>
	  <parameter name="MaxRecordNumber" value="0"/>
	  <parameter name="SkipNEvents" value="0"/>
	  <parameter name="SupressCheck" value="false"/>
	  <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE </parameter>
	  <parameter name="RandomSeed" value="1234567890" />
	</global>

	<processor name="MySDHCALProcessor" type="SDHCALProcessor">
	
	  <!-- Collection names -->
	  <parameter name="HCalCaloHitCollection" type="String"> ${collectionName} </parameter>
	  <parameter name="RootOuputFileName" type ="string"> ${rootFileName} </parameter>
	  <parameter name="RootTreeName" type ="string"> CaloHitAnalysis </parameter>
	  <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE DEBUG </parameter>
	</processor>
	
	<processor name="MyCutProcessor" type="CutProcessor">
	
	  <parameter name="SDHCALCollectionName" type="string" lcioInType="CalorimeterHit"> ${collectionName} </parameter>
	  
	      <parameter name="NHitCut" type="int"> 10 </parameter>
          <parameter name="NHitOverNLayerCut" type="double"> 4 </parameter>
          <parameter name="LayerFractionCut" type="double"> 0.2 </parameter>
          <parameter name="RadiusOverCog2Cut" type="double"> 0.4 </parameter>
          <parameter name="ShowerStartingLayerCut" type="int"> 4 </parameter>
          <parameter name="NTouchedLayersCut" type="int"> 8 </parameter>
          <parameter name="FractalTimesCentralCellsCut" type="double"> 0 </parameter>
          <parameter name="NHolesCut" type="int"> 1 </parameter>
          <parameter name="NHitEdgePercentCut" type="double"> 0.5 </parameter>
          <parameter name="BarycenterPositionCut" type="double"> 17.0 </parameter>
          <parameter name="LargeRMSCut" type="double"> 4.0 </parameter>
          <parameter name="CosThetaCut" type="double"> 0.95 </parameter>
          <parameter name="NeutralFirstLayerCut" type="int"> 3 </parameter>
          <parameter name="SpillTimeCut" type="double"> 2 </parameter>
          
          <parameter name="DecoderString" type="string"> M:3,S-1:3,I:9,J:9,K-1:6 </parameter>
          <parameter name="IJKEncoding"   type="StringVec"> I J K-1 </parameter>
	  
		  <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE </parameter>
	  
	</processor>
	
</marlin>

EOF

Marlin MarlinXML.xml
rm MarlinXML.xml


