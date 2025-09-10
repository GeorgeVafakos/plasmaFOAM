import os
import re

# ---------------------------------------------------------------
# Class: reactionExtractor
# ---------------------------------------------------------------
class reactionExtractor:

    # Class constructor
    # -----------------
    def __init__(self, dataPath):
        self._dataPath = dataPath
        self._availableReactions = dict()
        self.__availableReactionsFilename = str()


    # Helper method - Extract the reaction of a file by recognizing the identifier string 'PROCESS:'
    # -------------------------------------------------------------------------------------------
    def __extractReaction(self, speciesDir, filename):
        identifier = "PROCESS:"
        with open(os.path.join(self._dataPath,speciesDir,filename), 'r') as file:
            for line in file:
                if line.startswith(identifier):
                    # Extract and return the part after 'PROCESS:'
                    return line.split(identifier)[1].strip()  # Remove leading/trailing whitespace
        # Return if 'PROCESS:' is not found
        return 'NA'
    
    # Helper method - Extract the activation energy by finding the token before "eV"
    # ------------------------------------------------------------------------------
    def __extractEnergy(self, speciesDir, filename):
        identifier = "PARAM.:"
        with open(os.path.join(self._dataPath, speciesDir, filename), 'r') as file:
            for line in file:
                if line.startswith(identifier):
                    tokens = line.split()
                    for i, token in enumerate(tokens):
                        if token.strip().lower().rstrip(',:;.') == "ev" and i > 0:
                            return tokens[i - 1]  # Return the token before "eV"
                    return '--'  # No "eV" found
        return '--'
    
    # Helper method - Extract the database source using the "DATABASE:" line
    # ----------------------------------------------------------------------
    def __extractDatabase(self, speciesDir, filename):
        identifier = "DATABASE:"
        with open(os.path.join(self._dataPath, speciesDir, filename), 'r') as file:
            for line in file:
                if line.startswith(identifier):
                    return line.split(identifier)[1].strip()
        return '--'
    
    # Helper method - Custom sorting function to extract the numeric part
    # Used in the findReactions functino to sort the .dat files in ascending order
    # ----------------------------------------------------------------------------
    @staticmethod
    def __naturalKey(string):
        return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string)]
    

    # Create the dictionary with the reactions according to main species
    # ------------------------------------------------------------------
    def findReactions(self):
        # List all directories and files
        reactionDir = sorted(os.listdir(self._dataPath))
        
        # Filter out only the directories
        reactionDir = [entry for entry in reactionDir if os.path.isdir(os.path.join(self._dataPath, entry))]

        # Loop through the species directories
        for sp in reactionDir:
            # Get reaction files for each species
            reactionsFiles = sorted(os.listdir(os.path.join(self._dataPath, sp)), key=self.__naturalKey)
            self._availableReactions[sp] = dict()
            self._availableReactions[sp]['ID'] = list()
            self._availableReactions[sp]['file'] = list()
            self._availableReactions[sp]['reaction'] = list()
            self._availableReactions[sp]['energy'] = list()
            self._availableReactions[sp]['database'] = list()

            # Append reaction name and number
            reactionNum = 0
            for r in reactionsFiles:
                reactionNum += 1
                self._availableReactions[sp]['ID'].append(reactionNum)
                self._availableReactions[sp]['file'].append(r)
                self._availableReactions[sp]['reaction'].append(self.__extractReaction(sp, r))
                self._availableReactions[sp]['energy'].append(self.__extractEnergy(sp, r))
                self._availableReactions[sp]['database'].append(self.__extractDatabase(sp, r))
        return self._availableReactions
    

    # # Print available reactions to screen
    # # -----------------------------------
    # def printAvailableReactions(self):
    #     # Define headers
    #     headers = ['ID', 'Reaction', 'Energy', 'File']

    #     # Specify the length of each header (for Reaction header find the longest string and add more 5 spaces)
    #     lengthID = 6
    #     lengthEnergy = 10
    #     lengthFile = 20
    #     longestReaction = str()
    #     for sp in self._availableReactions.values():
    #         for reaction in sp['reaction']:
    #             # Update longest_reaction if the current reaction is longer
    #             if len(reaction) > len(longestReaction):
    #                 longestReaction = reaction
    #     lengthReaction = len(longestReaction) + 5

    #     # Print info
    #     print(f'The available reactions for BOLSIG+ are:')

    #     # Print headers
    #     print(f'{headers[0]:<{lengthID}} {headers[1]:<{lengthReaction}} {headers[2]:<{lengthEnergy}} {headers[3]:<{lengthFile}}')

    #     # Print a line separator
    #     separatorChar = '='
    #     print(f'{separatorChar*len(headers[0])}{' '*(lengthID-len(headers[0]))} {separatorChar*len(headers[1])}{' '*(lengthReaction-len(headers[1]))} {separatorChar*len(headers[2])}{' '*(lengthEnergy-len(headers[2]))} {separatorChar*len(headers[3])}{' '*(lengthFile-len(headers[3]))}')

    #     # Populate rows
    #     for sp, entry in self._availableReactions.items():
    #         for i in range(len(entry['ID'])):
    #             print(f'{entry['ID'][i]:<{lengthID}} {entry['reaction'][i]:<{lengthReaction}} {entry['energy'][i]:<{lengthEnergy}} ./{os.path.join(os.path.basename(self._dataPath), sp, entry['file'][i])}')
    #         print(f'')


    # Write available reactions to file
    # ---------------------------------
    def writeAvailableReactions(self, filename):
        self.__availableReactionsFilename = filename
        # Define headers
        headers = ['ID', 'Reaction', 'Energy', 'Database', 'File']

        # Specify the length of each header (for Reaction header find the longest string and add more 5 spaces)
        lengthID = 6
        lengthEnergy = 12
        lengthFile = 10

        longestDatabase = str()
        for sp in self._availableReactions.values():
            for database in sp['database']:
                if len(database) > len(longestDatabase):
                    longestDatabase = database
        lengthDatabase = len(longestDatabase) + 5

        longestReaction = str()
        for sp in self._availableReactions.values():
            for reaction in sp['reaction']:
                if len(reaction) > len(longestReaction):
                    longestReaction = reaction
        lengthReaction = len(longestReaction) + 5

        with open(os.path.join(self._dataPath, self.__availableReactionsFilename), 'w') as file:
            # Write headers
            file.write(f'{headers[0]:<{lengthID}} {headers[1]:<{lengthReaction}} {headers[2]:<{lengthEnergy}} {headers[3]:<{lengthDatabase}} {headers[4]:<{lengthFile}}\n')
            
            # Write a line separator
            separatorChar = '='
            file.write(f'{separatorChar*len(headers[0])}{" " * (lengthID-len(headers[0]))} {separatorChar*len(headers[1])}{" " * (lengthReaction-len(headers[1]))} {separatorChar*len(headers[2])}{" " * (lengthEnergy-len(headers[2]))} {separatorChar*len(headers[3])}{" " * (lengthDatabase-len(headers[3]))} {separatorChar*len(headers[4])}{" " * (lengthFile-len(headers[4]))}\n')            
            # Populate rows
            for sp, entry in self._availableReactions.items():
                for i in range(len(entry['ID'])):
                    file.write(f'{entry["ID"][i]:<{lengthID}} {entry["reaction"][i]:<{lengthReaction}} {entry["energy"][i]:<{lengthEnergy}} {entry["database"][i]:<{lengthDatabase}} ./{os.path.join(os.path.basename(self._dataPath), sp, entry["file"][i])}\n')
                file.write(f'\n')


# ---------------------------------------------------------------
# Class: reactionSelector
# ---------------------------------------------------------------
class reactionSelector:

    # Constructor
    # -----------
    def __init__(self, reactionExtractorObj, inputFilePath):
        self._inputFilePath = inputFilePath
        self._inputFileName = os.path.basename(self._inputFilePath)
        self._dataPath = reactionExtractorObj._dataPath
        self._availableReactions = reactionExtractorObj._availableReactions
        self._inputDirPath = os.path.dirname(self._inputFilePath)
        self._selectedReactions = dict()

    # Method to extract the reactions of a literature paper
    # -----------------------------------------------------
    def __extractReactionsFromLiterature(self, filename, keyword):
        selectedReactionsLocalDict = dict()
        keywordFound = False

        # Open file
        with open(filename, 'r') as file:
            lines = file.readlines()

            for lineNumber, line in enumerate(lines):
                # Remove white spaces at begging and end of line
                line = line.strip()

                # Ignore commented lines and empty lines
                if line.startswith('#') or not line.strip():
                    continue

                # Find the desired paper reactions
                if line.startswith(keyword+':'):
                    keywordFound = True  # Boolean to check if keyword exists in file
                    tempLine = line.split(':')  # Drop the keyword string from line
                    tempLine = tempLine[1].split(',')  # Get a list with the reactions IDs as strings
                    reactionsList = [int(s.strip()) for s in tempLine]  # Get a list with the reaction ID as ints

                    for sp, values in self._availableReactions.items():
                        selectedReactionsLocalDictValues = {'ID': [], 'file': [], 'reaction': [],  'energy': [],  'database': []}
                        for i, reactionID in enumerate(values['ID']):
                            if reactionID in reactionsList:
                                selectedReactionsLocalDictValues['ID'].append(reactionID)
                                selectedReactionsLocalDictValues['file'].append(values['file'][i])
                                selectedReactionsLocalDictValues['reaction'].append(values['reaction'][i])
                                selectedReactionsLocalDictValues['energy'].append(values['energy'][i])
                                selectedReactionsLocalDictValues['database'].append(values['database'][i])
                        if selectedReactionsLocalDictValues['ID']:  # Only add species with matching IDs
                            selectedReactionsLocalDict[sp] = selectedReactionsLocalDictValues
                    break

        # Check if keyword exists in file
        if not keywordFound:
            raise Exception(f"Keyword: {keyword+':'}: was not found in file: {filename}.")

        return selectedReactionsLocalDict


    # Method to read the input file
    # -----------------------------
    def readInputFile(self):
        selectionMethod = str() # Define empty string
        with open(self._inputFilePath, 'r') as file:
            lines = file.readlines()
            
            for lineNumber, line in enumerate(lines):
                # Remove white spaces at begging and end of line
                line = line.strip()
                
                # Ignore commented lines and empty lines
                if line.startswith('#') or not line.strip():
                    continue

                # Find the selectReactionMethod keyword
                if line.startswith('selectReactionMethod'):
                    # Check for exceptions in input file
                    if len(line.split()) > 2:
                        raise Exception(f'File:"{self._inputFilePath}". Invalid input in line: {lineNumber}')
                    else:
                        selectionMethod = line.split()[1]

                        # For 'custom' keyword
                        # if selectionMethod == 'custom':
                        #     for i in range(lineNumber+1, len(lines)):
                        #         lines[i] = lines[i].strip() # Get rid of trailing newline character (\n)
                        #         tempLine = lines[i].split() # Create this list only for the following if statement
                        #         if tempLine and tempLine[0].endswith(':'):
                        #             tempLine = str() # Repurpose the previous string to assign 
                        #             tempLine = lines[i].split(':',1)
                        #             sp = tempLine[0].strip()
                        #             tempReactions = tempLine[1].strip()
                        #             self._selectedReactions[sp] = dict()
                        #             self._selectedReactions[sp]['ID'] = [int(s.strip()) for s in tempReactions.split(',')]
                        #             self._selectedReactions[sp]['file'] = [file for i, file in zip(self._availableReactions[sp]['ID'], self._availableReactions[sp]['file']) if i in self._selectedReactions[sp]['ID']]
                        #             self._selectedReactions[sp]['reaction'] = [file for i, file in zip(self._availableReactions[sp]['ID'], self._availableReactions[sp]['reaction']) if i in self._selectedReactions[sp]['ID']]
                        #             self._selectedReactions[sp]['energy'] = [file for i, file in zip(self._availableReactions[sp]['ID'], self._availableReactions[sp]['energy']) if i in self._selectedReactions[sp]['ID']]
                        #             self._selectedReactions[sp]['database'] = [file for i, file in zip(self._availableReactions[sp]['ID'], self._availableReactions[sp]['database']) if i in self._selectedReactions[sp]['ID']]
                                    


                        if selectionMethod == 'custom':
                            for i in range(lineNumber + 1, len(lines)):
                                lines[i] = lines[i].strip()  # Remove leading/trailing spaces
                                tempLine = lines[i].split()
                                if tempLine and tempLine[0].endswith(':'):
                                    spLine = lines[i].split(':', 1)
                                    sp = spLine[0].strip()
                                    idList = [int(s.strip()) for s in spLine[1].strip().split(',')]
                                    self._selectedReactions[sp] = {key: [] for key in ['ID', 'file', 'reaction', 'energy', 'database']}

                                    for rid in idList:
                                        # Find the index of the reaction ID in availableReactions
                                        try:
                                            idx = self._availableReactions[sp]['ID'].index(rid)
                                        except ValueError:
                                            raise Exception(f'Reaction ID {rid} for species "{sp}" not found in available reactions.')

                                        self._selectedReactions[sp]['ID'].append(rid)
                                        self._selectedReactions[sp]['file'].append(self._availableReactions[sp]['file'][idx])
                                        self._selectedReactions[sp]['reaction'].append(self._availableReactions[sp]['reaction'][idx])
                                        self._selectedReactions[sp]['energy'].append(self._availableReactions[sp]['energy'][idx])
                                        self._selectedReactions[sp]['database'].append(self._availableReactions[sp]['database'][idx])





                                # Uncomment the else statement if you want to stop skipping empty lines between the reactions in the input file
                                # else: 
                                #     break
                        # For 'allAvailable' keyword
                        elif selectionMethod == 'allAvailable':
                            self._selectedReactions = self._availableReactions
                        # For 'literature' keyword
                        elif selectionMethod == 'literature':
                            tempLine = lines[lineNumber+1].split()

                            # Check that the the input arguments are exactly two
                            if len(tempLine)!=2:
                                raise Exception(f'File:"{self._inputFileName}". Arguments for selectionMethod "{selectionMethod}" must be exactly 2. The filename and the paper keyword.')

                            # Select the reactions according to input argument
                            literatureDataFile = tempLine[0]
                            literaturePaperKeyword = tempLine[1]
                            self._selectedReactions = self.__extractReactionsFromLiterature(os.path.join(self._dataPath,literatureDataFile), literaturePaperKeyword)
                        else:
                            raise Exception(f'File:"{self._inputFileName}". Unknown input "{selectionMethod}". Available options are "allAvailable", "literature" or "custom"')
                    
                    # If the selectReactionMethod is found, break the loop
                    break
                else:
                    raise Exception(f'File:"{self._inputFileName}". The selectReactionMethod entry not found.')
    

    # Create cross sections file
    # --------------------------
    def createCrossSectionFile(self, filename):
        self._crossSectionFileName = filename
        self._bolsigInputDirName = 'bolsigFiles'
        os.makedirs(os.path.join(self._inputDirPath, self._bolsigInputDirName), exist_ok=True)
        self._crossSectionFilePath = os.path.join(self._inputDirPath, self._bolsigInputDirName, self._crossSectionFileName)
        fileList = list()

        with open(self._crossSectionFilePath, 'w') as fout:
            fout.write('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n')

            # Define which reaction files will be used
            for sp in self._selectedReactions.keys():
                # Compare the selectedReactions against the availableReactions and select the corresponding files
                # fileList = [file for i, file in zip(self._availableReactions[sp]['ID'], self._availableReactions[sp]['file']) if i in self._selectedReactions[sp]]
                fileList = self._selectedReactions[sp]['file']
                fout.write('********************************************************** ')
                fout.write(sp)
                fout.write(' **********************************************************\n\n')
                for file in fileList:
                    with open(os.path.join(self._dataPath, sp, file), 'r') as fin:
                        fout.write(fin.read())
                        fout.write('\n')
                        
            fout.write('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')


    # Print selected reactions to screen
    # ----------------------------------
    def _printSelectedReactions(self):
        # Define headers
        headers = ['ID', 'Reaction', 'Energy', 'Database', 'File']

        # Specify the length of each header (for Reaction header find the longest string and add more 5 spaces)
        lengthID = 6
        lengthEnergy = 12
        lengthFile = 10

        longestDatabase = str()
        for sp in self._selectedReactions.values():
            for database in sp['database']:
                if len(database) > len(longestDatabase):
                    longestDatabase = database
        lengthDatabase = len(longestDatabase) + 5

        longestReaction = str()
        for sp in self._selectedReactions.values():
            for reaction in sp['reaction']:
                if len(reaction) > len(longestReaction):
                    longestReaction = reaction
        lengthReaction = len(longestReaction) + 5

        # Print info
        print(f'The selected reactions for BOLSIG+ are:\n')

        # Print headers
        print(f'{headers[0]:<{lengthID}} {headers[1]:<{lengthReaction}} {headers[2]:<{lengthEnergy}} {headers[3]:<{lengthDatabase}} {headers[4]:<{lengthFile}}')

        # Print a line separator
        separatorChar = '='
        print(f'{separatorChar*len(headers[0])}{" "*(lengthID-len(headers[0]))} {separatorChar*len(headers[1])}{" "*(lengthReaction-len(headers[1]))} {separatorChar*len(headers[2])}{" "*(lengthEnergy-len(headers[2]))} {separatorChar*len(headers[3])}{" "*(lengthDatabase-len(headers[3]))} {separatorChar*len(headers[4])}{" "*(lengthFile-len(headers[4]))}')

        # Populate rows
        for sp, entry in self._selectedReactions.items():
            for i in range(len(entry['ID'])):
                print(f'{entry["ID"][i]:<{lengthID}} {entry["reaction"][i]:<{lengthReaction}} {entry["energy"][i]:<{lengthEnergy}} {entry["database"][i]:<{lengthDatabase}} ./{os.path.join(os.path.basename(self._dataPath), sp, entry["file"][i])}')
            print(f'')


    # Write selected reactions to file
    # -------------------------------
    def writeSelectedReactions(self, filename):
        self._selectedReactionsFileName = filename
        self._selectedReactionsFilePath = os.path.join(self._inputDirPath, self._selectedReactionsFileName)

        # Define headers
        headers = ['ID', 'Reaction', 'Energy', 'Database', 'File']

        # Specify the length of each header (for Reaction header find the longest string and add more 5 spaces)
        lengthID = 6
        lengthEnergy = 12
        lengthFile = 10

        longestDatabase = str()
        for sp in self._selectedReactions.values():
            for database in sp['database']:
                if len(database) > len(longestDatabase):
                    longestDatabase = database
        lengthDatabase = len(longestDatabase) + 5

        longestReaction = str()
        for sp in self._selectedReactions.values():
            for reaction in sp['reaction']:
                if len(reaction) > len(longestReaction):
                    longestReaction = reaction
        lengthReaction = len(longestReaction) + 5

        with open(self._selectedReactionsFilePath, 'w') as file:
            # Write info
            file.write(f'The selected reactions for BOLSIG+ are:\n\n')
            
            # Write headers
            file.write(f'{headers[0]:<{lengthID}} {headers[1]:<{lengthReaction}} {headers[2]:<{lengthEnergy}} {headers[3]:<{lengthDatabase}} {headers[4]:<{lengthFile}}\n')
            
            # Write a line separator
            separatorChar = '='
            file.write(f'{separatorChar*len(headers[0])}{" " * (lengthID-len(headers[0]))} {separatorChar*len(headers[1])}{" " * (lengthReaction-len(headers[1]))} {separatorChar*len(headers[2])}{" " * (lengthEnergy-len(headers[2]))} {separatorChar*len(headers[3])}{" " * (lengthDatabase-len(headers[3]))} {separatorChar*len(headers[4])}{" " * (lengthFile-len(headers[4]))}\n')
            
            # Populate rows
            for sp, entry in self._selectedReactions.items():
                for i in range(len(entry['ID'])):
                    file.write(f'{entry["ID"][i]:<{lengthID}} {entry["reaction"][i]:<{lengthReaction}} {entry["energy"][i]:<{lengthEnergy}} {entry["database"][i]:<{lengthDatabase}} ./{os.path.join(os.path.basename(self._dataPath), sp, entry["file"][i])}\n')
                file.write(f'\n')

                