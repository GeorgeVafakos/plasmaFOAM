import os
import platform

class bolsigExecutor:
    
    # Class constructor
    # -----------------
    def __init__(self, reactionSelectorObj, filename):
        self._bolsigInputFileName = filename
        self._crossSectionFileName = reactionSelectorObj._crossSectionFileName
        self._crossSectionFilePath = reactionSelectorObj._crossSectionFilePath
        self._selectedReactions = reactionSelectorObj._selectedReactions
        self._inputFilePath = reactionSelectorObj._inputFilePath
        self._inputDirPath = reactionSelectorObj._inputDirPath
        self._inputFileName = reactionSelectorObj._inputFileName
        self._bolsigInputDirName = reactionSelectorObj._bolsigInputDirName
        self.__bolsigInputFilePath = os.path.join(self._inputDirPath, reactionSelectorObj._bolsigInputDirName, self._bolsigInputFileName)
        self.__bolsigRunType = 0 # 0=single run, 1=SERIESRUN 1D either E/N or one species mole fraction, 2=RUN2D->E/N and one species mole fraction
        self._bolsigOuputFileName = 'bolsigOutput.dat'
        self._printSelectedReactions = reactionSelectorObj._printSelectedReactions

    # Helper method to read the entries after a keyword in input file
    # Ouputs a list of the remaining entries in the same line (after the keyword)
    # ---------------------------------------------------------------------------
    def __readInputFile(self, keyword, multipleEntries = False):
        count = 0
        entry = list()
        singleEntry = list()
        with open(self._inputFilePath, 'r') as file:
            lines = file.readlines()

            for lineNumber, line in enumerate(lines):
                # Remove white spaces at begging and end of line
                line = line.strip()

                # Ignore commented lines and empty lines
                if line.startswith('#') or not line.strip():
                    continue

                # Find the keyword
                if line.startswith(keyword):
                    count += 1
                    singleEntry = line.split()
                    singleEntry.pop(0)
                    if not multipleEntries: 
                        entry = singleEntry
                        if count>1:
                            raise Exception(f'More than one entry for keyword"{keyword}" was found in line {lineNumber+1}. File: {self._inputFileName} \nTo allow multiple entries enable the "multipleEntries" argument')
                    else:
                        entry.append(singleEntry)
        return entry # Outputs list for single entry, or list of lists for multiple entries

    @staticmethod
    def __tabs(num = 6):
        return '\t'*num
    
    def __writeElectricFieldInputEntries(self, file):
        # Write electric field
        self.__redElecField = self.__readInputFile('E/N', multipleEntries = True) 
        for i, entry in enumerate(self.__redElecField):
            if len(self.__redElecField[i]) == 1:
                self.__redElecField[i][0] = float(self.__redElecField[i][0])
            elif len(self.__redElecField[i]) == 4:
                self.__redElecField[i][0] = float(self.__redElecField[i][0])
                self.__redElecField[i][1] = float(self.__redElecField[i][1])
                self.__redElecField[i][2] = int(self.__redElecField[i][2])
                if self.__redElecField[i][3] not in ['linear', 'quadratic', 'exponential']:
                    raise Exception(f'E/N grid spacing entry must be: "linear", "quadratic" or "exponential". Found "{self.__redElecField[i][3]}". File: {self._inputFileName}')
                if self.__bolsigRunType == 0: # Increase the self.__bolsigRunType only once
                    self.__bolsigRunType += 1
            else:
                raise Exception(f'E/N keyword must be followed by 1 or 4 values. Found {len(self.__redElecField[i])}. File: {self._inputFileName}')
        file.write('{:g}{}/ Electric field / N (Td)\n'.format(self.__redElecField[0][0], self.__tabs()))

    def __writeGasCompositionInputEntries(self, file):
        # Gas composition
        gasComposition = self.__readInputFile('gasComposition') # read whole entry of input file
        gasComposition = ' '.join(gasComposition) # make all entries a single string
        gasComposition = gasComposition.split(',') # force split using commas
        gasComposition = filter(None, gasComposition) # filter out empty elements
        gasComposition = [sp.strip() for sp in gasComposition] # remove spaces at begging and end
        # Check for extra species in the gasComposition input entry 
        speciesOfGasComposition = list()
        for i, sp in enumerate(gasComposition):
            speciesOfGasComposition.append(sp.split()[0])
            if speciesOfGasComposition[i] not in self._selectedReactions.keys():
                raise Exception(f'Species "{speciesOfGasComposition[i]}" does not exist in selectedReactions. File: {self._inputFileName}')
        # Check for missing species in the gasComposition input entry 
        for sp in self._selectedReactions.keys():
            if sp not in speciesOfGasComposition:
                raise Exception(f'Species "{sp}" in not declared in selectedReactions. File: {self._inputFileName}')
        # Check for invalid mole fraction entries and write mole fractions in input file
        moleFractions = dict()
        for entry in sorted(gasComposition):
            tempList = entry.split(' ')
            sp = tempList[0]
            tempList.pop(0)
            moleFractions[sp] = tempList
            if len(moleFractions[sp]) == 1:
                moleFractions[sp][0] = float(moleFractions[sp][0])
                file.write(f'{moleFractions[sp][0]} ')
            elif len(moleFractions[sp]) == 4:
                moleFractions[sp][0] = float(moleFractions[sp][0])
                moleFractions[sp][1] = float(moleFractions[sp][1])
                moleFractions[sp][2] = int(moleFractions[sp][2])
                if moleFractions[sp][3] not in ['linear', 'quadratic', 'exponential']:
                    raise Exception(f'gasComposition grid spacing of {sp} entry must be: "linear", "quadratic" or "exponential". Found "{moleFractions[sp][3]}". File: {self._inputFileName}')
                self.__bolsigRunType += 1
                speciesRUN2D = sp
                file.write(f'VAR ')
            else:
                raise Exception(f'gasComposition keyword must be followed by 1 or 4 values. Found {len(moleFractions[sp])} entries for {sp}. File: {self._inputFileName}')
        file.write(f'{self.__tabs(3)}/ Gas composition fractions\n')

    # Write the input file for BOSIG+
    # -------------------------------
    def writeBolsigInput(self):
        with open(self.__bolsigInputFilePath, 'w') as file:
            file.write(f'! Comment\n')
            file.write(f'/ Comment\n\n')
            file.write(f'NOSCREEN\n\n')
            file.write(f'/ NOLOGFILE\n\n')
            file.write(f'/READCOLLISIONS can be called multiple times to read from different files\n\n')
            file.write(f'/CLEARCOLLISIONS\n\n')
            file.write(f'READCOLLISIONS\n')
            file.write(f'{os.path.basename(self._crossSectionFileName).split('/')[-1]}\t\t/ File\n')
            file.write(' '.join(sp for sp in sorted(list(self._selectedReactions.keys()))))
            file.write(f'{self.__tabs(4)}/ Species\n')
            file.write('{:d}{}/ Extrapolate: 0= No 1= Yes\n\n'.format(1, self.__tabs()))

            file.write(f'CONDITIONS\n')

            # # Write electric field
            # redElecField = self.__readInputFile('E/N', multipleEntries = True) 
            # for i, entry in enumerate(redElecField):
            #     if len(redElecField[i]) == 1:
            #         redElecField[i][0] = float(redElecField[i][0])
            #     elif len(redElecField[i]) == 4:
            #         redElecField[i][0] = float(redElecField[i][0])
            #         redElecField[i][1] = float(redElecField[i][1])
            #         redElecField[i][2] = int(redElecField[i][2])
            #         if redElecField[i][3] not in ['linear', 'quadratic', 'exponential']:
            #             raise Exception(f'E/N grid spacing entry must be: "linear", "quadratic" or "exponential". Found "{redElecField[i][3]}". File: {self._inputFileName}')
            #         if self.__bolsigRunType == 0: # Increase the self.__bolsigRunType only once
            #             self.__bolsigRunType += 1
            #     else:
            #         raise Exception(f'E/N keyword must be followed by 1 or 4 values. Found {len(redElecField[i])}. File: {self._inputFileName}')
            # file.write('{:g}{}/ Electric field / N (Td)\n'.format(redElecField[0][0], self.__tabs()))

            self.__writeElectricFieldInputEntries(file)

            file.write('{:g}{}/ Angular field frequency / N (m3/s)\n'.format(0, self.__tabs()))
            file.write('{:g}{}/ Cosine of E-B field angle\n'.format(0, self.__tabs()))

            # Gas temperature
            gasTemp = self.__readInputFile('gasTemp')
            if len(gasTemp) > 1:
                raise Exception(f'The gas temperature entries are more than one. File: {self._inputFileName}')
            file.write('{:g}{}/ Gas temperature (K)\n'.format(float(gasTemp[0]), self.__tabs()))

            file.write('{:g}{}/ Excitation temperature (K)\n'.format(300.0, self.__tabs()))
            file.write('{:g}{}/ Transition energy (eV)\n'.format(0, self.__tabs()))
            file.write('{:g}{}/ Ionization degree\n'.format(0, self.__tabs()))

            # Gas particle density
            gasDensity = self.__readInputFile('gasDensity')
            if len(gasDensity) > 1:
                raise Exception(f'The gas density entries are more than one. File: {self._inputFileName}')
            file.write('{:g}{}/ Gas particle density (1/m3)\n'.format(float(gasDensity[0]), self.__tabs(4)))

            file.write('{:g}{}/ Ion charge parameter\n'.format(-1, self.__tabs()))
            file.write('{:g}{}/ Ion/neutral mass ratio\n'.format(1, self.__tabs()))
            file.write('{:d}{}/ e-e momentum effects & modified Coulomb logarithm: 0=No&No; 1=Yes&No; 2=No&Yes; 3=Yes&Yes*\n'.format(0, self.__tabs()))
            file.write('{:d}{}/ Energy sharing: 1=Equal*; 2=One takes all\n'.format(1, self.__tabs()))
            file.write('{:d}{}/ Growth: 1=Temporal*; 2=Spatial; 3=Not included; 4=Grad-n expansion\n'.format(1, self.__tabs()))
            file.write('{:g}{}/ Maxwellian mean energy (eV) \n'.format(0, self.__tabs()))
            file.write('{:d}{}/ # of grid points\n'.format(100, self.__tabs()))
            file.write('{:d}{}/ Manual grid: 0=No; 1=Linear; 2=Parabolic \n'.format(0, self.__tabs()))
            file.write('{:g}{}/ Manual maximum energy (eV)\n'.format(200, self.__tabs()))
            file.write('{:g}{}/ Precision\n'.format(1e-10, self.__tabs(5)))
            file.write('{:g}{}/ Convergence\n'.format(1e-5, self.__tabs(5)))
            file.write('{:d}{}/ Maximum # of iterations\n'.format(100, self.__tabs()))

            # Gas composition
            gasComposition = self.__readInputFile('gasComposition') # read whole entry of input file
            gasComposition = ' '.join(gasComposition) # make all entries a single string
            gasComposition = gasComposition.split(',') # force split using commas
            gasComposition = filter(None, gasComposition) # filter out empty elements
            gasComposition = [sp.strip() for sp in gasComposition] # remove spaces at begging and end
            # Check for extra species in the gasComposition input entry 
            speciesOfGasComposition = list()
            for i, sp in enumerate(gasComposition):
                speciesOfGasComposition.append(sp.split()[0])
                if speciesOfGasComposition[i] not in self._selectedReactions.keys():
                    raise Exception(f'Species "{speciesOfGasComposition[i]}" does not exist in selectedReactions. File: {self._inputFileName}')
            # Check for missing species in the gasComposition input entry 
            for sp in self._selectedReactions.keys():
                if sp not in speciesOfGasComposition:
                    raise Exception(f'Species "{sp}" in not declared in selectedReactions. File: {self._inputFileName}')
            # Check for invalid mole fraction entries and write mole fractions in input file
            moleFractions = dict()
            for entry in sorted(gasComposition):
                tempList = entry.split(' ')
                sp = tempList[0]
                tempList.pop(0)
                moleFractions[sp] = tempList
                if len(moleFractions[sp]) == 1:
                    moleFractions[sp][0] = float(moleFractions[sp][0])
                    file.write(f'{moleFractions[sp][0]} ')
                elif len(moleFractions[sp]) == 4:
                    moleFractions[sp][0] = float(moleFractions[sp][0])
                    moleFractions[sp][1] = float(moleFractions[sp][1])
                    moleFractions[sp][2] = int(moleFractions[sp][2])
                    if moleFractions[sp][3] not in ['linear', 'quadratic', 'exponential']:
                        raise Exception(f'gasComposition grid spacing of {sp} entry must be: "linear", "quadratic" or "exponential". Found "{moleFractions[sp][3]}". File: {self._inputFileName}')
                    self.__bolsigRunType += 1
                    speciesRUN2D = sp
                    file.write(f'VAR ')
                else:
                    raise Exception(f'gasComposition keyword must be followed by 1 or 4 values. Found {len(moleFractions[sp])} entries for {sp}. File: {self._inputFileName}')
            file.write(f'{self.__tabs(3)}/ Gas composition fractions\n')

            file.write('{:d}{}/ Normalize composition to unity: 0=No;\n\n'.format(1, self.__tabs()))

            # Write the RUN method
            if self.__bolsigRunType == 0:
                file.write(f'RUN\n\n')
            elif self.__bolsigRunType == 1:
                for i, entry in enumerate(self.__redElecField):
                    file.write(f'RUNSERIES\n')
                    file.write('{:d}{}/ Variable: 1=E/N; 2=Mean energy; 3=Maxwellian energy\n'.format(1, self.__tabs()))
                    file.write('{:g} {:g}{}/ Min Max\n'.format(self.__redElecField[i][0], self.__redElecField[i][1], self.__tabs(5)))
                    file.write('{:d}{}/ Number\n'.format(self.__redElecField[i][2], self.__tabs()))
                    for j, gridType in enumerate(['linear','quadratic','exponential']):
                        if self.__redElecField[i][3] == gridType:
                            file.write(f'{j+1}{self.__tabs()}/ Type1 Type2: 1=Linear; 2=Quadratic; 3=Exponential\n\n')
            elif self.__bolsigRunType == 2:
                file.write(f'RUN2D\n')
                file.write('{:d}{}/ First variable: 1=E/N; 2=Mean energy; 3=Maxwellian energy \n'.format(1, self.__tabs()))
                file.write('{:g} {:g} {:g} {:g}{}/ Min1 Max1 Min2 Max2\n'.format(self.__redElecField[0][0], self.__redElecField[0][1], moleFractions[speciesRUN2D][0], moleFractions[speciesRUN2D][1], self.__tabs(4)))
                file.write('{:d} {:d}{}/ Num1 Num2\n'.format(self.__redElecField[0][2], moleFractions[speciesRUN2D][2], self.__tabs(5)))
                for i, gridType in enumerate(['linear','quadratic','exponential']):
                    if self.__redElecField[0][3] == gridType:
                        file.write(f'{i+1} ')
                for i, gridType in enumerate(['linear','quadratic','exponential']):
                    if moleFractions[speciesRUN2D][3] == gridType:
                        file.write(f'{i+1} ')
                file.write(f'{self.__tabs(5)}/ Type1 Type2: 1=Linear; 2=Quadratic; 3=Exponential\n')
                file.write(f'bolsigOutRun2D.dat{self.__tabs(2)}/ Output file of RUN2D\n\n') # A secondary output file, for additional info
            else:
                raise Exception('Error with the "__bolsigRunType" variable.')

            # Write the output file
            file.write(f'SAVERESULTS\n')
            file.write(f'{self._bolsigOuputFileName}\t\t/ File\n')
            file.write('{:d}{}/ Format: 1=Run by run; 2=Combined; 3=E/N; 4=Energy; 5=SIGLO; 6=PLASIMO\n'.format(2, self.__tabs()))
            file.write('{:d}{}/ Conditions: 0=No; 1=Yes\n'.format(1, self.__tabs()))
            file.write('{:d}{}/ Transport coefficients: 0=No; 1=Yes\n'.format(1, self.__tabs()))
            file.write('{:d}{}/ Rate coefficients: 0=No; 1=Yes\n'.format(1, self.__tabs()))
            file.write('{:d}{}/ Reverse rate coefficients: 0=No; 1=Yes\n'.format(1, self.__tabs()))
            file.write('{:d}{}/ Energy loss coefficients: 0=No; 1=Yes\n'.format(0, self.__tabs()))
            file.write('{:d}{}/ Distribution function: 0=No; 1=Yes \n'.format(0, self.__tabs()))
            file.write('{:d}{}/ Skip failed runs: 0=No; 1=Yes\n'.format(0, self.__tabs()))
            file.write('{:d}{}/ Include cross sections: 0=No; 1=Yes\n\n'.format(0, self.__tabs()))

            file.write('END')

        # Print screen available reactions according to input file
        # --------------------------------------------------------
        boolMap = { 'True': True, 'true': True, 'TRUE': True, 'False': False}
        boolValue = boolMap.get(self.__readInputFile('printScreenSelectedReactions')[0])
        if boolValue:
            self._printSelectedReactions()


    # Run the BOLSIG+ solver
    # ----------------------
    def runBolsigCode(self):
        # Find BOLSIG+ executable path
        bolsigExecutablePath = os.path.join(os.path.dirname(__file__),'bolsigExe','bolsigminus')

        # Execute BOLSIG+ code according to current OS
        try:
            if platform.system() == 'Linux':
                os.chdir(os.path.join(self._inputDirPath, self._bolsigInputDirName)) # Change direction to the case
                os.system(f'{bolsigExecutablePath} {self._bolsigInputFileName}')
            elif platform.system() == 'Windows':
                raise Exception(f'Windows OS not yet supported.')
            else:
                raise Exception(f'Current OS not supported. Available OS are: UNIX/Linux.')
        except:
            raise Exception('Problem with BOLSIG+ code execution.')


    
                
