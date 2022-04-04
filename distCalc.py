####################
## distCalc 0.1.0 ##
####################

##Start of file##

##Importing libraries##
import pymzml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

##Options##
pd.options.mode.chained_assignment = None

##Object Defnitions##
class ExperimentalPeak:
  def __init__(self, mz, inte):
    self.mz = mz
    self.inte = inte
   
class MzmlData:
  def __init__(self, scan_id, spec_id, peaks, retention_times):
    self.__scan_id = scan_id
    self.__spec_id = spec_id
    self.__peaks = peaks
    self.__retention_times = retention_times
    

  def __len__(self):
    return len(self.__retention_times)
  
  @property
  def scan_id(self): 
    return self.__scan_id
  
  @property
  def spec_id(self): 
    return self.__spec_id

  @property
  def peaks(self): 
    return self.__peaks
  
  @property
  def retention_times(self): 
    return self.__retention_times
  
  @classmethod
  def read_mzML(cls, mzmlFile, MS_Level):
    run = pymzml.run.Reader(mzmlFile, MS1_Precision = 5e-6, MSn_Precision = 20e-6)
    retention_times = []
    spec_ids = []
    scan_ids = []
    peaks = []
    i = 0
    for s in run:
      if s.ms_level == MS_Level:
        spec_ids.append(i)
        scan_ids.append(s.id_dict['scan'])
        retention_times.append(s.scan_time[0]*60)
        peaks.append(MzmlData._get_peaks(s))
        i = i + 1  
    return cls(scan_ids, spec_ids, peaks, retention_times)
  
  @staticmethod
  def _get_peaks(spectrum):
    peaks = spectrum.peaks("raw").tolist()
    spectral_peaks = [ExperimentalPeak(tmp[0], tmp[1]) for tmp in peaks]
    return spectral_peaks

class distribution:
    def __init__(self, dataFrame, avgMass):
        self.__dataFrame=dataFrame
        self.__avgMass=avgMass
    
    @property
    def dataFrame(self):
        return self.__dataFrame
    
    @property
    def avgMass(self):
        return self.__avgMass
    
    @property
    def graph(self):
        if len(self.dataFrame) > 1:
            print("Graph of distribution:")
            plt.rcParams.update({'font.size': 7})
            plt.bar(self.dataFrame.mzr, self.dataFrame.Intensity,width=.02, linewidth=1, color='darkorchid')
            return plt.show()
        if len(self.foundDataFrame) == 1:
            print("Graph of distribution:")
            newFoundDF=pd.DataFrame({"mzr": [(self.dataFrame.mzr[0]-1), self.dataFrame.mzr[0], (self.dataFrame.mzr[0]+1)], "Intensity": [0, self.dataFrame.Intensity[0], 0]})
            plt.rcParams.update({'font.size': 7})
            plt.bar(newFoundDF.mzr, newFoundDF.Intensity,width=.02, linewidth=1, color='darkorchid')
            return plt.show()
    
    @classmethod
    def distCompile(cls, foundDataFrame):
        foundDataFrame=foundDataFrame.drop_duplicates()
        intensitySum=sum(foundDataFrame.Intensity)
        foundDataFrame['RelativeAbundance'] = (foundDataFrame.Intensity/intensitySum)
        weights = foundDataFrame.mzr * foundDataFrame.RelativeAbundance
        foundAvgMass=(sum(weights))
        return cls(foundDataFrame, foundAvgMass)
   
##Function definitions## 
   
def preprocess(mzmlFile):
  mzml_data = MzmlData.read_mzML(mzmlFile, MS_Level = 1)
  return mzml_data

def findSpecIndex(preprocessedMZMLFile, scanNumber):
    specIndex=preprocessedMZMLFile.scan_id.index(scanNumber - 1)
    return(specIndex)

def distCalc(preprocessedMZMLFile, specIndex, precursorMZR, charge):
    peakNumber=len(preprocessedMZMLFile.peaks[specIndex])
    peaks=preprocessedMZMLFile.peaks[specIndex]
    df = pd.DataFrame(columns = ['mzr', 'Intensity'])
    for x in range(peakNumber):
      i=peaks[x].inte
      mzr=float("{:.2f}".format(peaks[x].mz))
      if i>0:
        tempDF=pd.DataFrame({'mzr': [mzr], 'Intensity': [i]})
        df = pd.concat([df, tempDF], axis=0, ignore_index=True)
    df=df.apply(pd.to_numeric)
    index = abs(df.mzr - precursorMZR).idxmin()
    tempDF=pd.DataFrame(columns=['mzr', 'Intensity'])
    foundMZR=df.mzr[index]
    foundDF = pd.DataFrame({'mzr':[foundMZR], 'Intensity': [df.Intensity[index]]})
    step=1
    while step <= 40:
        y=0
        tempDF=pd.DataFrame(columns=['mzr', 'Intensity'])
        for x in np.arange((foundMZR +((1.04/charge)*step)), (foundMZR + ((1.10/charge)*step)), 0.001):
            condition = np.any(np.isclose(x, df.mzr, atol=.2))
            if condition == True:
                index = abs(df['mzr'] - x).idxmin()
                tempDF=pd.concat([tempDF, pd.DataFrame({'mzr' : [df.mzr[index]], 'Intensity' : [df.Intensity[index]]})], axis = 0, ignore_index=True)
                tempDF=tempDF.apply(pd.to_numeric)
        for x in tempDF.Intensity:
            if x>y:
                y=x
        if y >= 1:
            index = abs(tempDF['Intensity'] - y).idxmin()
            foundDF=pd.concat([foundDF, pd.DataFrame({'mzr':[tempDF.mzr[index]],'Intensity':[y]})], axis=0, ignore_index=True)
            step += 1
        else:
            break
        del(tempDF)
    dist=distribution.distCompile(foundDF)
    return dist
        
def distComp(preprocessedMZMLFile, specIndex, precursorMZR, charge, secondPreprocessedMZMLFile=None, thirdPreprocessedMZMLFile=None):
  if secondPreprocessedMZMLFile==None and thirdPreprocessedMZMLFile==None:
    firstDist=distCalc(preprocessedMZMLFile, specIndex, precursorMZR, charge)
    return firstDist
  if secondPreprocessedMZMLFile!=None and thirdPreprocessedMZMLFile==None:
    firstDist=distCalc(preprocessedMZMLFile, specIndex, precursorMZR, charge)
    retentionTime=preprocessedMZMLFile.retention_times[specIndex]
    secondIndexes=[]
    for x in range(len(secondPreprocessedMZMLFile.spec_id)):
        condition = np.any(np.isclose(retentionTime, secondPreprocessedMZMLFile.retention_times[x], atol=30))
        if condition == True:
            secondIndexes.append(x)
    dists=[]
    for x in secondIndexes:
        dist = distCalc(secondPreprocessedMZMLFile, x, precursorMZR, charge)
        dists.append(dist)
    dists.sort(key=lambda x: sum(x.dataFrame.Intensity), reverse=True)
    dists=dists[0:3]
    for x in range(3):
        cD=firstDist.avgMass-dists[x].avgMass
        print("Centroid distance:", cD)
        print("Average deuterium atoms incorporated:", cD*charge)
    distList=[firstDist, dists]
    return distList
  if secondPreprocessedMZMLFile!=None and thirdPreprocessedMZMLFile!=None:
      firstDist=distCalc(preprocessedMZMLFile, specIndex, precursorMZR, charge)
      retentionTime=preprocessedMZMLFile.retention_times[specIndex]
      secondIndexes=[]
      for x in range(len(secondPreprocessedMZMLFile.spec_id)):
          condition = np.any(np.isclose(retentionTime, secondPreprocessedMZMLFile.retention_times[x], atol=30))
          if condition == True:
              secondIndexes.append(x)
      dists1=[]
      for x in secondIndexes:
          dist = distCalc(secondPreprocessedMZMLFile, x, precursorMZR, charge)
          dists1.append(dist)
      dists1.sort(key=lambda x: sum(x.dataFrame.Intensity), reverse=True)
      dists1=dists1[0:3]
      for x in range(3):
          cD=firstDist.avgMass-dists1[x].avgMass
          print("Centroid distance to first comparison:", cD)
          print("Average deuterium atoms incorporated for first comparison:", cD*charge)
      thirdIndexes=[]
      for x in range(len(thirdPreprocessedMZMLFile.spec_id)):
          condition = np.any(np.isclose(retentionTime, thirdPreprocessedMZMLFile.retention_times[x], atol=30))
          if condition == True:
              thirdIndexes.append(x)
      dists2=[]
      for x in secondIndexes:
          dist = distCalc(thirdPreprocessedMZMLFile, x, precursorMZR, charge)
          dists2.append(dist)
      dists2.sort(key=lambda x: sum(x.dataFrame.Intensity), reverse=True)
      dists2=dists2[0:3]
      for x in range(3):
          cD=firstDist.avgMass-dists2[x].avgMass
          print("Centroid distance to second comparison:", cD)
          print("Average deuterium atoms incorporated for second comparison:", cD*charge)
      distList=[firstDist, dists1, dists2]
      return distList
        
        
def distWrite(distribution, outputPath):
    if isinstance(distribution, list)==True:
        with open(outputPath, 'w') as f:
            for x in range(len(distribution)):   
                dFStr=distribution[x].dataFrame.to_string()
                f.write("Data frame of peaks:")
                f.write('\n')
                f.write(dFStr)
                f.write('\n')
                f.write("Average mass:")
                f.write(str(distribution[x].avgMass))
                f.write('\n')
                f.write('\n')
    else:
        with open(outputPath, 'w') as f:
            dFStr=distribution.dataFrame.to_string()
            f.write("Data frame of peaks:")
            f.write('\n')
            f.write(dFStr)
            f.write('\n')
            f.write("Average mass:")
            f.write(str(distribution.avgMass))
            
##End of file##