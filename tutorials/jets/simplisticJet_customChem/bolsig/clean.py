import os
import shutil
import run  # To import the inputFileName variable

def deleteFiles(directory, excludeList):
    """
    Deletes all files and directories except those in the exclude list.
    
    Args:
        directory (str): The directory to clean up.
        exclude (list): A list of item names to exclude.
    """
    try:
        # Convert exclude list to full paths for comparison
        excludePaths = [os.path.join(directory, item) for item in excludeList]

        for item in os.listdir(directory):
            itemPath = os.path.join(directory, item)

            # Skip items in the exclude list
            if itemPath in excludePaths:
                continue

            # Remove file or directory
            if os.path.isfile(itemPath):
                os.remove(itemPath)
                print(f'Deleted: {itemPath}')
            elif os.path.isdir(itemPath):
                shutil.rmtree(itemPath)
                print(f'Deleted: {itemPath}')
    except Exception as e:
        print(f'An error occurred: {e}')

# Example usage
if __name__ == "__main__":
    # Specify the directory
    currenPath = os.path.dirname(os.path.abspath(__file__))

    # Specify what to keep
    excludedItems = ['run.py', 'README.md', os.path.basename(__file__), run.inputFile]

    # Clean case
    deleteFiles(currenPath, excludedItems)
