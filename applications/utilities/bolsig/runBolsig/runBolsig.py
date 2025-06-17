import os
from LXCatDataHandler import reactionExtractor, reactionSelector
from bolsigInputHandler import bolsigExecutor 
from bolsigOutputHandler import bolsigOutputReader 

def runBolsigMain(projectPath, inputFilePath):
    # ---------------------------------------------------------------
    # --- Read reactions from database
    # ---------------------------------------------------------------
    # Constructor to initialize reactionExtractor class
    LXCatDataPath = os.path.join(projectPath, 'LXCatData')
    LXCatData = reactionExtractor(LXCatDataPath)

    # Find reactions in the directory
    LXCatData.findReactions()

    # Write all avaiable reactions to file
    LXCatData.writeAvailableReactions('availableReactions.txt')



    # ---------------------------------------------------------------
    # --- Select the desired reactions 
    # ---------------------------------------------------------------
    # Constructor to initialize reactionSelector class
    selectReactions = reactionSelector(LXCatData, inputFilePath)

    # Read input file to select reactions
    selectReactions.readInputFile()

    # Write the cross-sections file for BOLSIG+
    selectReactions.createCrossSectionFile('crossSections.txt')

    # Write to file the selected reactions
    selectReactions.writeSelectedReactions('selectedReactions.txt')


    # ---------------------------------------------------------------
    # --- Create BOLSIG+ input and run
    # ---------------------------------------------------------------
    # Constructor to initialize reactionSelector class
    bolsigInput = bolsigExecutor(selectReactions, 'bolsigInput.dat')

    # Write BOLSIG+ input file
    bolsigInput.writeBolsigInput()

    # Run BOLSIG+ code
    bolsigInput.runBolsigCode()


    # ---------------------------------------------------------------
    # --- Read BOLSIG+ output and hadle data
    # ---------------------------------------------------------------
    # Constructor to initialize readBolsigOutput class
    bolsigData = bolsigOutputReader(bolsigInput)

    # Process the output file of BOLSIG+ and extract tables
    bolsigData.processFile()

    # Write BOLSIG+ coefficients table in file
    #- 1st input: file name
    #- 2nd input (optional): choose between csv (default) and dat file formats
    #- 3rd input (optional): write the column header - true (default) or false 
    bolsigData.writeResults('openfoam', 'bolsigDict')


    print('\nEND')


