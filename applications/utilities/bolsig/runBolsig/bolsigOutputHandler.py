import os
import csv

class bolsigOutputReader:
    
    # Class constructor
    # -----------------
    def __init__(self, bolsigExecutorObj):
        # Data from BOLSIG+ output
        self.moleFractionsTable = list()
        self.transportPropertiesTable = list()
        self.rateConstantsTable = list()
        self.bolsigCoeffsTable = list()
        self.bolsigCoeffsTableColumnNames = list()
        self._coeffsTableFileName = str()
        self._bolsigOuputFileName = bolsigExecutorObj._bolsigOuputFileName 
        self._inputDirPath = bolsigExecutorObj._inputDirPath 

    # Method to check if string is a number
    # -------------------------------------
    @staticmethod
    def __isNumber(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    # Check if the end of table has been reached
    # ------------------------------------------
    def __checkEndOfTable(self, line):
        splitLine = line.strip().split()
        if not splitLine:
            return True  # Return True to indicate the end of the table in case of empty line
        return not self.__isNumber(splitLine[0])

    # Method to find the mole fractions' names
    # ----------------------------------------
    def __findMoleFractions(self):
        moleFractionFound = False
        with open(self._bolsigOuputFileName, 'r') as file:
            for line in file:
                parts = line.split()
                if len(parts) >= 2 and parts[1] == 'Mole' and parts[2] == 'fraction':
                    moleFractionFound = True
                    self.__moleFractionsBolsigOutputEntries['bolsigOutputHeaders'].append(parts[0])
                    self.__moleFractionsBolsigOutputEntries['Species'].append(parts[3])
                    if len(parts)==4:
                        self.__moleFractionsBolsigOutputEntries['moleFractionValue'].append(float('NaN'))
                    elif len(parts)==5:
                        self.__moleFractionsBolsigOutputEntries['moleFractionValue'].append(parts[4])
                    else:
                        raise RuntimeError(f'Error in BOLSIG+ output file. Problem with entry {parts[0]} for species {part[3]}.')
            # Extend coefficients table header list with found species
            self.bolsigCoeffsTableColumnNames.extend(self.__moleFractionsBolsigOutputEntries['Species']) 
        # Check if there are species entries in BOLSIG+ output file
        if not moleFractionFound:
            raise RuntimeError('No mole fractions found.')
        

    # Method to find the reactions names
    # ----------------------------------
    def __findReactionNames(self):
        reactionNamesFound = False
        with open(self._bolsigOuputFileName, 'r') as file:
            for line in file:
                parts = line.split()
                if parts and parts[0] == 'Rate':
                    reactionNamesFound = True
                    for line in file:
                        parts = line.split()
                        if parts and parts[0] == 'R#':
                            break
                        self.bolsigCoeffsTableColumnNames.append(parts[0])
        if not reactionNamesFound:
            raise RuntimeError('No reactions found.')
        
    # Method to find the E/N value, when it is constant throughout the BOLSIG+ runs
    # -----------------------------------------------------------------------------
    def __findReducedElectricField(self):
        with open(self._bolsigOuputFileName, 'r') as file:
            for line in file:
                parts = line.split()
                if len(parts)>=2 and parts[1] == 'Electric' and parts[2] == 'field':
                    E_N = parts[-1]
                    break
        return float(E_N)

    # Helper method to process a table
    # --------------------------------
    def __readTabularData(self, file, table):
        for line in file:
            if self.__checkEndOfTable(line):
                break
            row = [float(x) for x in line.split() if self.__isNumber(x)]
            if row:
                table.append(row)

    # Public method to process the BOLSIG+ output file
    # ------------------------------------------------
    def processFile(self):
        # Define a dictionary for the mole fractions (populated in the __findMoleFractions function)
        self.__moleFractionsBolsigOutputEntries = dict()
        self.__moleFractionsBolsigOutputEntries['Species'] = list()
        self.__moleFractionsBolsigOutputEntries['bolsigOutputHeaders'] = list()
        self.__moleFractionsBolsigOutputEntries['moleFractionValue'] = list()

        # Chech if BOLSIG+ output file exists
        if not os.path.isfile(self._bolsigOuputFileName):
            raise FileNotFoundError('File not found')

        # Create bolsigCoeffsTable table columns names
        if not self.bolsigCoeffsTableColumnNames:
            self.bolsigCoeffsTableColumnNames.append('E/N_(Td)')
            self.__findMoleFractions()
            self.bolsigCoeffsTableColumnNames.extend(['MeanEnergy_(eV)', 'Mobility*N_((1/m/V/s))', 'Diffusion*N_(1/m/s)'])
            self.__findReactionNames()
        else:
            raise RuntimeError('Error when initializing bolsigCoeffsTableColumnNames vector. Vector was not initially empty.')

        # Reading table data
        with open(self._bolsigOuputFileName, 'r') as file:
            count = 0
            for line in file:
                line = line.strip()

                # Check for lines starting with 'R#'
                if line.startswith('R#'):
                    count += 1

                    # Get the E/N and mole fractions
                    if count == 1:
                        self.__readTabularData(file, self.moleFractionsTable) 
                        headerList = line.split() # Get the headers of the first table
                        headerList.pop(0) # Delete the "R#" entry for run number

                        # When BOLSIG+ performs a single RUN then the headerList will contain only the "R#" element
                        if headerList:
                            headerList.pop(0) # Delete the A1" entry for E/N in case of RUN2D or SERIESRUN

                        # Write to coeffs table the E/N value for the case of a single BOLSIG+ RUN
                        if not headerList:
                            self.bolsigCoeffsTable.append([self.__findReducedElectricField()])

                        if set(headerList) == set(self.__moleFractionsBolsigOutputEntries['bolsigOutputHeaders']):
                            for row in self.moleFractionsTable:
                                self.bolsigCoeffsTable.append(row[1:])
                        else:
                            for i, row in enumerate(self.moleFractionsTable):
                                self.bolsigCoeffsTable.append(row[1:2])
                                self.bolsigCoeffsTable[i].extend(self.__moleFractionsBolsigOutputEntries['moleFractionValue'])

                    # Get the mean energy, electron transport and diffusion coefficients
                    elif count == 2:
                        self.__readTabularData(file, self.transportPropertiesTable)
                        for i, row in enumerate(self.transportPropertiesTable):
                            self.bolsigCoeffsTable[i].extend(row[2:5])

                    # Get the reaction rates coefficients
                    elif count == 3:
                        self.__readTabularData(file, self.rateConstantsTable)
                        for i, row in enumerate(self.rateConstantsTable):
                            self.bolsigCoeffsTable[i].extend(row[3:])
                        break

    # Method to print a table
    # -----------------------
    @staticmethod
    def printTable(table):
        for row in table:
            print(' '.join(map(str, row)))

    # # Write BOLSIG+ data table
    # # ------------------------
    # def writeResults(self, filename, fileFormat='csv', headers=True):
    #     self._coeffsTableFileName = filename
    #     self._bolsigCoeffsDirName = 'bolsigGeneratedCoeffs'

    #     # Create directory to store coefficients table
    #     os.makedirs(os.path.join(self._inputDirPath, self._bolsigCoeffsDirName), exist_ok=True)

    #     # Check if the fileFormat argument is correct
    #     if fileFormat not in ["csv", "dat"]:
    #         raise ValueError('File format must be "csv" or "dat".')

    #     separator = ',' if fileFormat == 'csv' else ' '
    #     self._coeffsTableFileName = f'{self._coeffsTableFileName}.{fileFormat}'
    #     self._coeffsTableFilePath = os.path.join(self._inputDirPath, self._bolsigCoeffsDirName, self._coeffsTableFileName)

    #     with open(self._coeffsTableFilePath, mode='w', newline='') as file:
    #         writer = csv.writer(file, delimiter=separator)
    #         if headers:
    #             writer.writerow(self.bolsigCoeffsTableColumnNames)
    #         writer.writerows(self.bolsigCoeffsTable)

    #     print(f'Data written to {self._coeffsTableFileName}')



    def writeResults(self, fileFormat, filename):
        """
        Write the BOLSIG+ coefficients table to a file.
        
        Args:
            filename (str): The name of the output file.
            fileFormat (str): The format of the output file ('csv', 'dat', or 'openfoam').
            headers (bool): Whether to include headers in the output file.
        """

        self._coeffsTableFileName = filename

        # Create directory to store BOLSIG+ results
        self._bolsigOuputDirName = os.path.join(self._inputDirPath, 'bolsigOutput')
        os.makedirs(self._bolsigOuputDirName, exist_ok=True)

        # Dispatch based on format
        writers = {
            "csv": self._writeCsv,
            "dat": self._writeDat,
            "openfoam": self._writeOpenFoamDict
        }

        if fileFormat not in writers:
            raise ValueError('File format must be "csv", "dat", or "openfoam".')

        # Call the appropriate writer method
        writers[fileFormat]()


    def _writeCsv(self):
        """
        Write the BOLSIG+ coefficients table to a CSV file.
        """
        # Define the path for the CSV file
        path = os.path.join(self._bolsigOuputDirName, self._coeffsTableFileName + '.csv')

        # Write the coefficients table to the CSV file
        with open(path, mode='w', newline='') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(self.bolsigCoeffsTableColumnNames)
            writer.writerows(self.bolsigCoeffsTable)
        print(f'Data written to {path}')


    def _writeDat(self):
        """
        Write the BOLSIG+ coefficients table to a .dat file using fixed-width columns.
        """
        path = os.path.join(self._bolsigOuputDirName, self._coeffsTableFileName + '.dat')

        with open(path, mode='w') as file:
            # Write header with fixed width
            header_line = ''.join(f"{col:>26}" for col in self.bolsigCoeffsTableColumnNames)
            file.write(header_line + '\n')

            # Write rows with scientific notation (e.g., 1.23e+04) and fixed width
            for row in self.bolsigCoeffsTable:
                line = ''.join(f"{float(val):26.10e}" for val in row)
                file.write(line + '\n')

        print(f'Data written to {path}')


    def _writeOpenFoamDict(self):
        """
        Write the BOLSIG+ coefficients table to an OpenFOAM dictionary file.
        """
        path = os.path.join(self._bolsigOuputDirName, 'bolsigProperties')

        with open(path, mode='w') as file:
            file.write(r'/*--------------------------------*- C++ -*----------------------------------*\ ' + '\n')
            file.write(r'| =========                 |                                                 |' + '\n')
            file.write(r'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |' + '\n')
            file.write(r'|  \\    /   O peration     | Version:  2406                                  |' + '\n')
            file.write(r'|   \\  /    A nd           | Website:  www.openfoam.com                      |' + '\n')
            file.write(r'|    \\/     M anipulation  |                                                 |' + '\n')
            file.write(r'\*---------------------------------------------------------------------------*/' + '\n')
            file.write(r'FoamFile' + '\n')
            file.write(r'{' + '\n')
            file.write(r'    version     2.0;' + '\n')
            file.write(r'    format      ascii;' + '\n')
            file.write(r'    class       dictionary;' + '\n')
            file.write(r'    object      bolsigProperties;' + '\n')
            file.write(r'}' + '\n')
            file.write(r'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //' + '\n\n')

            # Write EN block
            self.__writeENBlock(file)

            print(type(file))

            print(f'Data written to {path}')


    def __writeENBlock(self, file):
        """
        Write the unique E/N values from the coefficients table to an OpenFOAM-style dictionary block.

        Args:
            file : File object to write the E/N values to.
        """
        uniqueEN_values = set()
        EN_values = []

        for row in self.bolsigCoeffsTable:
            val = float(row[0])
            if val not in uniqueEN_values:
                uniqueEN_values.add(val)
                EN_values.append(val)

        file.write('EN\n(\n')
        for val in EN_values:
            file.write(f'    {val:.6f}\n')
        file.write(');\n\n')